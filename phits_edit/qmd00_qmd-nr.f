************************************************************************
*                                                                      *
      subroutine jqmdin(ityp,ktyp,eein,mmas,mchg,bmax0)
*                                                                      *
*       control routine of JQMD                                        *
*       modified by K.Niita on 2005/08/15                              *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*        icfg = 1, 2, 4                                                *
*                                                                      *
*              input       : content                     ; variables   *
*                                                                      *
*              'proj'      : projectile                  ; mstq1(1)    *
*                                                        ; mstq1(2)    *
*                                                        ; mstq1(3)    *
*              'targ'      : target                      ; mstq1(4)    *
*                                                        ; mstq1(5)    *
*                                                        ; mstq1(6)    *
*              'event'     : number of events            ; mstq1(7)    *
*              'tstep'     : total number of time step   ; mstq1(8)    *
*              'frame'     : reference frame             ; mstq1(9)    *
*                                                                      *
*                                                                      *
*              'win'       : incident energy or momentum ; parq1(1)    *
*                                                        ; parq1(2)    *
*              'bmin'      : minimum impact parameter    ; parq1(3)    *
*              'bmax'      : maximum impact parameter    ; parq1(4)    *
*              'dt'        : time step                   ; parq1(5)    *
*                                                                      *
*              'fname(i)'  : file name                   ; fname(i)    *
*              'mstq1(i)'  : integer parameters          ; mstq1(i)    *
*              'parq1(i)'  : real parameters             ; parq1(i)    *
*                                                                      *
*     output: in cldist                                                *
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
      use QMD_COOD2_MOD, only : coodalloc, lfirst

      implicit real*8(a-h,o-z)

      parameter ( nnnn = 800 )      ! S.Abe 2020/09/29

*-----------------------------------------------------------------------

      include 'param00.inc'
      include 'param01.inc'

      parameter ( pi  = 3.1415926535898d0 )

*-----------------------------------------------------------------------

      common /swich2/ icfg, imany, icpus, idatm
!$OMP THREADPRIVATE(/swich2/)
      common /input1/ mstq1(mxpa1), parq1(mxpa1)
!$OMP THREADPRIVATE(/input1/)

      common /swich3/ ielst, jelst, kelst
!$OMP THREADPRIVATE(/swich3/)

      common /clustf/ nclst, iclust(nnn)
!$OMP THREADPRIVATE(/clustf/)
      common /clustg/ jclust(0:8,nnn),  qclust(0:12,nnn)
!$OMP THREADPRIVATE(/clustg/)

      common /framtr/ betafr(0:2), gammfr(0:2)
!$OMP THREADPRIVATE(/framtr/)
      common /crshi/  bplus, icrhi, ijudg, imadj, iqmax

      common /pnint/  ipnint
      integer ipnint

      common /muflag/ imuinthit,imubrmhit,imuppdhit,imucaphit,imucapflag
!$OMP THREADPRIVATE(/muflag/)

      common /qmdscm/ qmdscm_h0, qmdscm_d, qmdscm_rcls, iqmdscm
      common /ldpart/ rlpjdg,ilpflg,idlp(nnnn)
!$OMP THREADPRIVATE(/ldpart/)
      common /lclsti/ eklcp(0:nnnn),ilclst(0:nnnn,0:nnnn)
!$OMP THREADPRIVATE(/lclsti/)

*-----------------------------------------------------------------------
*     input data and initialization
*-----------------------------------------------------------------------

               IF ( .not. lfirst ) then
                lfirst = .true.
                call coodalloc
                if( iqmdscm .eq. 1 ) then
                   do k = 1, nnnn
                      idlp(k) = 0
                   enddo
                   do i = 0, nnnn
                      eklcp(i) = 0.d0
                      do j = 0, nnnn
                         ilclst(i,j) = 0
                      enddo
                   enddo
                endif
               ENDIF

               ieo  = 6
               ierr = 0

               mstq1(1) = 0
               mstq1(2) = ibryf(ityp,ktyp)
               mstq1(3) = ichgf(ityp,ktyp)

               msp = mstq1(3)
               msn = mstq1(2) - mstq1(3)

               if( ( msp .eq. 0 .and. msn .gt. 1 ) .or.
     &             ( msp .gt. 1 .and. msn .le. 0 ) ) then

                  write(ieo,'(/
     &                '' **** Collision is skipped by JQMD ****''/
     &                '' projectile charge = '',i4/
     &                ''            mass   = '',i4)') msp, msp+msn

                  nclst = -2
                  return

               end if

               mstq1(4) = 0
               mstq1(5) = mmas
               mstq1(6) = mchg

               if( ( mchg .eq. 0 .and. mmas-mchg .gt. 1 ) .or.
     &             ( mchg .gt. 1 .and. mmas-mchg .le. 0 ) ) then

                  write(ieo,'(/
     &                '' **** Collision is skipped by JQMD ****''/
     &                ''     target charge = '',i4/
     &                ''            mass   = '',i4)') mchg, mmas

                  nclst = -2
                  return

               end if

*-----------------------------------------------------------------------
*           bmax for QMD calculation
*           ielst, moving frame adjust
*-----------------------------------------------------------------------

               ap = mstq1(2)
               at = mstq1(5)


               bmax = 1.15 * ( ap**(1./3.) + at**(1./3.) )
     &              - 0.4 + bplus



*-----------------------------------------------------------------------

               mstq1(190) = imadj    ! 0:no, 1:moving frame adjust
               mstq1(191) = ijudg    ! 0:no, 1:new ielst judge

               mstq1(17)  = 1        ! ielst for judge

*-----------------------------------------------------------------------

               if ( mstq1(2) /= 0. ) then
                  parq1(1)   = eein / dble( mstq1(2) ) / 1000.0d0
               else
                  parq1(1) = 0.0
               endif

               parq1(2)   = -1.0
               parq1(3)   = 0.0
               parq1(4)   = bmax
               mstq1(11)  = 0

               mstq1(7)   = 1
               mstq1(8)   = iqmax
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

               bmax = 1.2 * ( ap**(1./3.) + at**(1./3.) )
     &              + bplus

               parq1(4)   = bmax

            end if

*-----------------------------------------------------------------------

               parq1(120) = 1.0      ! sdmemin

               icfg = 4

               call qmdint

cABE add @2014/05/20, change @2014/08/14, projectile is photon, mod @2015/02/13
               if ( (ipnint .ge. 1 .or. imuinthit .eq. 1 .or.
     &               imucaphit .eq. 1) .and. ityp .eq. 14 ) then
                  masspr = 0
                  msprpr = 0
               endif

*-----------------------------------------------------------------------
*           QMD Events
*-----------------------------------------------------------------------

               call qmdevent(ierr,ityp)

                  if( ierr .ne. 0 )  goto 999
                  if( nclst .lt. 0 ) goto 999

*-----------------------------------------------------------------------
*           elastic judge and frame trnasform to lab system
*-----------------------------------------------------------------------

               call qmdjudge

                  if( kelst .eq. 1 ) goto 999

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

      return

*-----------------------------------------------------------------------
*     error or nothing happened
*-----------------------------------------------------------------------

  999 continue

                  nclst = -1

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine qmdevent(ierr,ityp)
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 26                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to simulate one event by QMD                            *
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8(a-h,o-z)

      include 'param00.inc'

      parameter ( nnnn = 800 )      ! S.Abe 2020/09/29

*-----------------------------------------------------------------------

      common /const1/ elab, rdist, bmin, bmax, ibch, ibin
!$OMP THREADPRIVATE(/const1/)
      common /const2/ dt, ntmax, iprun, iprun0
!$OMP THREADPRIVATE(/const2/)
      common /const3/ nfreq, nfrec, nfred
!$OMP THREADPRIVATE(/const3/)
      common /vriab1/ b, llnow, ntnow
!$OMP THREADPRIVATE(/vriab1/)
      common /swich1/ ipot, insys, irkg, icolt
!$OMP THREADPRIVATE(/swich1/)

      common /impact/ bval(1000), ibnum(1000), bweight(1000), bdef
!$OMP THREADPRIVATE(/impact/)

      common /rannum/ iseed, iseed0, iseed1
!$OMP THREADPRIVATE(/rannum/)
      common /coln01/ iccoll
!$OMP THREADPRIVATE(/coln01/)

      data ibsf /0/
      save ibsf
!$OMP THREADPRIVATE(ibsf)

      common /pnint/  ipnint
      integer ipnint

      common /pionflag/ iqdflag, ipionflag, istrflag, iphreturn
!$OMP THREADPRIVATE(/pionflag/)

      common /muflag/ imuinthit,imubrmhit,imuppdhit,imucaphit,imucapflag
!$OMP THREADPRIVATE(/muflag/)

      common /gg005/ xxx, yyy, zzz, uuu, vvv, www, tme, erg,
     &               dls, wgt, vel, dtc,
     &               icl, iii, jjj, kkk, jsu, iap, jgp, ipt,
     &               mtp, iexp, iex, idx,
     &               npa, ncp
!$OMP THREADPRIVATE(/gg005/)
      common /qmdsuc/ isucc
!$OMP THREADPRIVATE(/qmdsuc/)
      common /clustf/ nclst, iclust(nnn)
!$OMP THREADPRIVATE(/clustf/)

      common /const0/ idnta, idnpr, massta, masspr, mstapr, msprpr
!$OMP THREADPRIVATE(/const0/)

      common /ldpart/ rlpjdg,ilpflg,idlp(nnnn)
!$OMP THREADPRIVATE(/ldpart/)
      common /qmdscm/ qmdscm_h0, qmdscm_d, qmdscm_rcls, iqmdscm

      data iqmstac /0/ ! 0: stability check OFF(original),  1: stability check ON(revised)

*-----------------------------------------------------------------------
*        Initialization of one QMD event
*-----------------------------------------------------------------------

*-----------------------------------------------------------------------
*           event number, initial time, and initial randum seed
*           total collision flag
*-----------------------------------------------------------------------

                  ierr =  0

                  ntnow  = 0
                  iseed1 = iseed
                  iccoll = 0

*-----------------------------------------------------------------------
*           impact parameter and multi run control
*-----------------------------------------------------------------------

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
*           make ground state and boost
*-----------------------------------------------------------------------

               icount = 0 ! counter to record # of trials to make a stable ground state
               isucc  = 0 ! flag to record whether ground state setting is successful. 0: No, 1: Target OK, 2: Target and projectile OK
               ifinal = 0 ! flag to finish ground state setting and start final collision simulation
               do ! Loop to search stable ground + final collision

                icount = icount + 1

                if(icount .eq. 100 .or. isucc .eq. 2) then
                 ifinal = 1
                 call ground_final(ierr) ! make ground state according to previously examined stable state

cABE add @2014/05/20, change @2014/08/14, to avoid stability check
                  if ( (ipnint .ge. 1 .or. imuinthit .eq. 1) .and.
     &                ityp .eq. 14 .and. massta .ge. 2 .and.
     &                iqdflag .eq. 1 ) then
                     call qdeuteron( uuu, vvv, www, erg)
                  end if

                  if ( (ipnint .ge. 1 .or. imuinthit .eq. 1) .and.
     &                ityp .eq. 14 .and. ipionflag .eq. 1 ) then
                     call pionprod
                  end if

                  if ( imucaphit .eq. 1 .and. ityp .eq. 14 ) then
                     call pnconvert_qmd
                  end if

                  if( ierr .ne. 0 ) return

                ELSE
                 call ground_pre(ierr) ! make ground state for stability test
                 if( ierr .ne. 0 ) return ! S.Abe 2015/11/24
                ENDIF

                IF(iqmstac .eq. 0 .and. ifinal .eq. 0) goto 200  ! if iqmstac == 0, skip stability check time evolution


                  call rboost(ierr)

                     if( ierr .ne. 0 ) return

                  if( iqmdscm .eq. 1 ) call lpinit

*-----------------------------------------------------------------------
*        Time Evolution
*-----------------------------------------------------------------------

         do 100 nt = 1, ntmax

               ntnow = nt

*-----------------------------------------------------------------------
*           time integration
*-----------------------------------------------------------------------

               if( irkg .eq. 2 ) then

                  call rk12(dt)

               else if( irkg .eq. 4 ) then

                  call rkg4(dt)

               end if

*-----------------------------------------------------------------------
*           collision term
*-----------------------------------------------------------------------

               if( icolt .eq. 1 ) then

                  call pionem(dt)
                  call relcol
                  call pionab

               end if

               if( iqmdscm .eq. 1 ) call lpjudge

*-----------------------------------------------------------------------

  100    continue

*-----------------------------------------------------------------------
*        Final pion decay and Final analysis of clusters
*-----------------------------------------------------------------------

                  call fpidecay

               if( iqmdscm .ne. 1) then
                  call cldist
               else
                  call cldist3
               endif


*-----------------------------------------------------------------------

*-----------------------------------------------------------------------
*-----Judgement of ground state setting---------------------------------
  200          IF(ifinal .eq. 0) then ! this trial is ground state setting

                IF(nclst .ne. 1 .and. iqmstac .ge. 1) then ! ground state was broken!! -> try another ground state
                 nclst = -1
                 cycle
                ELSE ! this ground state is good
                 nclst = -1
cABE change @2014/05/20, change @2014/08/14, to skip projectile stability check for photonuclear reaction, mod @2015/02/13
                 if ( (ipnint .ge. 1 .or. imuinthit .eq. 1 .or.
     &                 imucaphit .eq. 1) .and. ityp .eq. 14 ) then
                    isucc = isucc + 2
                 else
                    isucc = isucc + 1
                 endif
                 cycle
                ENDIF

               ELSE ! this trial is for collision (true run)
                exit
               ENDIF

               enddo ! Loop to (search stable ground) + (collision)
*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine qmdjudge
*                                                                      *
*                                                                      *
*        Last Revised:     2008 06 19                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to determin weight of the event, [qmdfac] and           *
*              to judge the elastic or inelastic reaction type and     *
*              to sum up the reaction cross section                    *
*                                                                      *
*                                                                      *
*        Variables: in common block /swich3/                           *
*                                                                      *
*              ielst       : input flag, jelst < ielst: elastic        *
*              jelst       : elastic or inelastic flag                 *
*                     = 0  : elastic without collision                 *
*                     = 1  : elastic with collision                    *
*                     = 2  : inelastic without collision               *
*                     = 3  : inelastic with collision                  *
*              kelst       : 1 -> elastic, 0-> inelastic               *
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8(a-h,o-z)

*-----------------------------------------------------------------------

      include 'param00.inc'
      include 'param01.inc'

*-----------------------------------------------------------------------

      common /const0/ idnta, idnpr, massta, masspr, mstapr, msprpr
!$OMP THREADPRIVATE(/const0/)
      common /const1/ elab, rdist, bmin, bmax, ibch, ibin
!$OMP THREADPRIVATE(/const1/)

      common /vriab3/ qmdfac, sdmfac
!$OMP THREADPRIVATE(/vriab3/)

      common /clustf/ nclst, iclust(nnn)
!$OMP THREADPRIVATE(/clustf/)
      common /clustg/ jclust(0:8,nnn),  qclust(0:12,nnn)
!$OMP THREADPRIVATE(/clustg/)

      common /coln01/ iccoll
!$OMP THREADPRIVATE(/coln01/)
      common /swich3/ ielst, jelst, kelst
!$OMP THREADPRIVATE(/swich3/)
      common /sdmcut/ sdmemin
!$OMP THREADPRIVATE(/sdmcut/)

      common /impact/ bval(1000), ibnum(1000), bweight(1000), bdef
!$OMP THREADPRIVATE(/impact/)

      common /input1/ mstq1(mxpa1), parq1(mxpa1)
!$OMP THREADPRIVATE(/input1/)

      common /muflag/ imuinthit,imubrmhit,imuppdhit,imucaphit,imucapflag
!$OMP THREADPRIVATE(/muflag/)

      common /pionflag/ iqdflag, ipionflag, istrflag, iphreturn
!$OMP THREADPRIVATE(/pionflag/)

*-----------------------------------------------------------------------

            if( nclst .lt. 0 ) then

               kelst = 1
               return

            end if

      if( imucaphit .eq. 1 ) then
       kelst = 0
       return
      endif

      if( ipionflag .eq. 1 ) then
       kelst = 0
       return
      endif

*-----------------------------------------------------------------------

                     iselc = mstq1(191)

*-----------------------------------------------------------------------
*        judgement of elastic or inelastic reaction
*-----------------------------------------------------------------------

                     iels = 1

            if( nclst .eq. 2 ) then

               if( ( jclust(1,1) .eq. mstapr .and.
     &               jclust(2,1) .eq. massta - mstapr .and.
     &               jclust(1,2) .eq. msprpr .and.
     &               jclust(2,2) .eq. masspr - msprpr )     .or.
     &             ( jclust(1,2) .eq. mstapr .and.
     &               jclust(2,2) .eq. massta - mstapr .and.
     &               jclust(1,1) .eq. msprpr .and.
     &               jclust(2,1) .eq. masspr - msprpr ) )   then

               if( iselc .eq. 0 ) then

                     if( qclust(6,1) .lt. sdmemin .and.
     &                   qclust(6,2) .lt. sdmemin  ) then

                        iels = 0

                     end if

               else

                  if( massta .gt. 1 .and. masspr .gt. 1 .and.
     &                iccoll .eq. 0 ) then

                        remin1 = 0.3d0 * ( jclust(1,1) + jclust(2,1) )
                        remin2 = 0.3d0 * ( jclust(1,2) + jclust(2,2) )

                     if( qclust(6,1) .lt. remin1 .and.
     &                   qclust(6,2) .lt. remin2  ) then

                        iels = 0

                     end if

                  else

                     if( qclust(6,1) .lt. sdmemin .and.
     &                   qclust(6,2) .lt. sdmemin  ) then

                        iels = 0

                     end if

                  end if

               end if

               end if

            end if

*-----------------------------------------------------------------------

            if( iccoll .eq. 0 ) then

               if( iels .eq. 0 ) then

                  jelst = 0

               else

                  jelst = 2

               end if

            else

               if( iels .eq. 0 ) then

                  jelst = 1

               else

                  jelst = 3

               end if

            end if

*-----------------------------------------------------------------------
*        elastic flag
*-----------------------------------------------------------------------

            if( jelst .ge. ielst ) then

               kelst = 0

            else

               kelst = 1

            end if

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine trfram(px,py,pz,et,rm,ifrm)
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 27                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to determine energy and momentum by Lorentz transform   *
*                                                                      *
*        Variables:                                                    *
*                                                                      *
*              px, py, pz     : momentum                               *
*              et             : energy                                 *
*              rm             : rest mass                              *
*              ifrm           : 0-> to lab                             *
*                               1-> to cm                              *
*                               2-> to n-n                             *
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8(a-h,o-z)

*-----------------------------------------------------------------------

      common /framtr/ betafr(0:2), gammfr(0:2)
!$OMP THREADPRIVATE(/framtr/)

*-----------------------------------------------------------------------

         if( ifrm .lt. 0 .or. ifrm .gt. 2 ) then

            write(*,*) ' **** Error at [trfram], unrecognized frame '
            stop 999

         end if

*-----------------------------------------------------------------------

            gamm = gammfr(ifrm)
            beta = betafr(ifrm)

*-----------------------------------------------------------------------

            pz = pz * gamm - beta * gamm * et
            et = sqrt( px**2 + py**2 + pz**2 + rm**2 )

*-----------------------------------------------------------------------

      return
      end


