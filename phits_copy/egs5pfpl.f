!-----------------------------egs5_photon.f-----------------------------
! Version: 051219-1435
!          080425-1100   Add time as the time after start.
!          091105-0835   Remove redundant tvstep/ustep equivalence
! Reference: SLAC-R-730/KEK-2005-8
!-----------------------------------------------------------------------
!23456789|123456789|123456789|123456789|123456789|123456789|123456789|12

      subroutine egs5pfpl(fpl,eig,wgt,irl,totpni,totegs)
!
      USE egs5_media_mod !<- 2015.08xx allocatable
      USE egs5_edge_mod  !<- 2015.08xx allocatable
      USE egs5_photin_mod  !<- 2015.08xx allocatable
      implicit none
!      save
      include 'include/egs5_h.f'               ! Main EGS5 "header" file

      include 'include/egs5_bounds.f'    ! COMMONs required by EGS5 code
      include 'include/egs5_epcont.f'
      include 'include/egs5_misc.f'
      include 'include/egs5_stack.f'
      include 'include/egs5_uphiot.f'
      include 'include/egs5_useful.f'
      include 'include/egs5_uservr.f'

      include 'include/counters.f'       ! Additional (non-EGS5) COMMONs

      real*8 rnnow                                           ! Arguments
      integer ircode,iarg

      real*8                                           ! Local variables
     * eig,                                  ! Energy of incident photon
     * gbr1,gbr2,t,temp,
     * bexptr,dpmfp,gmfpr0,gmfp,cohfac
      integer idr,lgle,irl,iij

      real*8 EPSGMFP                                  ! Local parameters
      data EPSGMFP/1.D-6/                     ! Smallest gamma mfp value

      integer medin
      real*8 fpl     ! flight path in PHITS
      real*8 pstep   ! photon step, ne 0 after cross-surface
      real*8 wgt
      common /egs5cmn5/ircode
!$OMP THREADPRIVATE(/egs5cmn5/)
      common /egs5cmn6/pstep,dpmfp,gmfp
!$OMP THREADPRIVATE(/egs5cmn6/)

      integer ifgsq
      common /egs5cmn9/ifgsq
!$OMP THREADPRIVATE(/egs5cmn9/)


      real*8 totttl,totegs,totpni

      iphoton = iphoton + 1                ! Count entry into subroutine

      if(np.eq.0) np=1 ! not initialized yet ! T.Sato 2022/10/27

! iwase reset dE
      edep = 0d0

      ir(np) =  irl
      medium =  med(irl)
      gle    = log(eig)
      wt(np) = wgt


! T.Sato 2019/07/06, add dpmfp condition because it sometimes take over from previous history value
      if(ifgsq .eq. 1.and.dpmfp.ne.0) goto 3

      ircode = 1                              ! Set up for normal return


!                                           ! ---------------------------
!                                           ! ---------------------------

                                              ! ------------------------
 2    continue                                ! Start of NEW-ENERGY loop
                                              ! ------------------------

!     ------------------------------------------------------
!     Sample number of mfp's to transport before interacting
!     ------------------------------------------------------
      call randomset(rnnow)
      if (rnnow .eq. 0.) rnnow = 1.E-30
      dpmfp = -log(rnnow)                              ! Number of mfp's

      if (cexptr .ne. 0.) then        ! Apply exponential transformation
        if (w(np) .gt. 0.) then
          temp = cexptr*w(np)
          bexptr = 1./(1. - temp)
          dpmfp = dpmfp*bexptr               ! Number of mfp's (revised)
          wt(np) = wt(np)*bexptr*exp(-dpmfp*temp)    ! Associated weight
        end if
      end if

      irold = ir(np)                        ! Initialize previous region

                                              ! ------------------------
 3    continue                                ! Start of NEW-MEDIUM loop
                                              ! ------------------------
                                              ! Here each time we change
                                              ! medium during transport
                                              ! ------------------------

      if (medium .ne. 0) then                        ! Set PWLF interval
        lgle = ge1(medium)*gle + ge0(medium)
        iextp=0
        if (eig .lt. 0.15) then
          do iij=1,nedgb(medium)
            if (ledgb(iij,medium) .eq. lgle) then
              if (edgb(iij,medium) .le .eig) then
                iextp = 1
              else
                iextp = -1
              end if
            end if
          end do
        end if
        gmfpr0 = gmfp1(lgle+iextp,medium)*gle +
     *           gmfp0(lgle+iextp,medium)
        if(gmfpr0.le.0.d0) then
          iextp=0
          gmfpr0 = gmfp1(lgle+iextp,medium)*gle +
     *             gmfp0(lgle+iextp,medium)
        end if
      end if

 4    continue

      if (medium .le. 0) then      ! Set for large vacuum-step transport
        tstep = vacdst               ! Distance to next interaction (cm)
      else                                            ! Normal transport
        rhof = rhor(irl)/rhom(medium)                    ! Density ratio
        gmfp = gmfpr0/rhof                                  ! Scaled mfp
        if (iraylr(irl) .eq. 1) then         ! Apply Rayleigh correction
          cohfac = cohe1(lgle+iextp,medium)*gle +
     *             cohe0(lgle+iextp,medium)
          gmfp = gmfp*cohfac                             ! Corrected mfp
        end if
        if ( gmfp .gt. 0.0d0 ) then
          totegs = 1.0d0 / gmfp
        else
          totegs = 0.0d0
        endif

        totttl = totegs + totpni
        if ( totttl .gt. 0.0d0 ) then
          tstep = (1.0d0/totttl)*dpmfp           ! Distance to next interaction (cm)
        else
          tstep = 0.0d0
        end if
      endif
      fpl   = tstep
      pstep = fpl
      wgt   = wt(np)

      return
      end
!
