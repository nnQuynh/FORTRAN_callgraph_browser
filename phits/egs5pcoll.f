      subroutine egs5pcoll(eig,u0,v0,w0,wt0,irl,medin,nbeta)
      USE egs5_edge_mod   !<-
      USE egs5_photin_mod !<-
      USE egs5_thresh_mod !<-

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

      integer nnn,nomp,nclsts,iclusts,jclusts,i
      real*8 qclusts
      include 'param00.inc'
      common /clustt/ nclsts, iclusts(nnn)
!$OMP THREADPRIVATE(/clustt/)
      common /clustw/ jclusts(0:8,nnn),qclusts(0:12,nnn)
!$OMP THREADPRIVATE(/clustw/)

      real*8 x0,y0,z0,u0,v0,w0,wt0
      real*8 charge,pr,ekine,etot,rmg,pxl,pyl,pzl,rms
      integer medin,nnpp,kf,itype,id,k
      integer nbeta,mark,iegsxray
      real*8 pstep
      common /egs5cmn5/ircode
!$OMP THREADPRIVATE(/egs5cmn5/)
      common /egs5cmn6/pstep,dpmfp,gmfp
!$OMP THREADPRIVATE(/egs5cmn6/)
      common /egs5cmn8/iegsxray
!$OMP THREADPRIVATE(/egs5cmn8/)

      integer ncpls
      common /anglth/ cosph,costh,sinph,sinth
      real*8 cosph,costh,sinph,sinth,ein
!$OMP THREADPRIVATE(/anglth/)

      medium = medin

      if(medium .eq. 0) then
         write(506,*)'error in egs5pcoll, stopped'
         write(*,*)'error in egs5pcoll, stopped'
         stop
      endif
      np = 1
      iq(np) = 0
       e(np) = eig


       u(np) =  u0
       v(np) =  v0
       w(np) =  w0

      wt(np) = wt0
      ir(np) = irl
      gle    = log(eig)

      lgle = ge1(medium)*gle + ge0(medium)

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

!     ------------------------------------------------------------------
!     It is finally time to interact  ---  determine type of interaction
!     ------------------------------------------------------------------
!     First check if it is a Rayleigh scatter (provided option is ON)
!     ---------------------------------------------------------------
      if (iraylr(irl) .eq. 1) then
         cohfac = cohe1(lgle+iextp,medium)*gle +
     *             cohe0(lgle+iextp,medium)
        call randomset(rnnow)
        if (rnnow .le. (1.0 - cohfac)) then

!         ===========
          call raylei
!         ===========


          goto 10

        end if
      end if

!     ------------------------------------------------------------------
!     Otherwise determine if PAIR, COMPTON, or PHOTOELECTRIC interaction
!     ------------------------------------------------------------------
!     ---------------------
!     PAIR production check
!     ---------------------
      call randomset(rnnow)
      gbr1 = gbr11(lgle+iextp,medium)*gle + gbr10(lgle+iextp,medium)
      if (rnnow. le .gbr1 .and. e(np) .gt. RMT2) then

!       =========
        call pair               ! To determine energies and polar angles
!       =========

        ircode = 2
        goto 10

      end if

!     -------------
!     COMPTON check
!     -------------
      gbr2 = gbr21(lgle+iextp,medium)*gle + gbr20(lgle+iextp,medium)
      if (rnnow .lt. gbr2) then

!       ==========
        call compt
!       ==========


                               ! ---------------------------------------
      else                     ! Must be PHOTOELECTRIC (only thing left)
                               ! ---------------------------------------

!                                  =================
!       ==========
        call photo
!     flag x-ray event to tell PHITS
          do i=1,np
             if(iq(i).eq.0)iegsxray = 1
          enddo


!       ==========

        if (np .eq. 0) then
          ircode = 2                 ! ---------------------------------
          goto 10                    ! Stack is EMPTY - return to SHOWER
        end if                       !----------------------------------

                                                      ! ----------------
      end if                                          ! ----------------


      if (eig .lt. pcut(irl)) then
         ircode = 2
         nbeta = 3
      endif
                                                 ! ---------------------

!-----------------------last line of egs5_photon.f----------------------

 10   continue

           ncpls = 1

      call egs2phits(ncpls)

      costh = 1.0
      sinth = 0.0
      cosph = 1.0
      sinph = 0.0

         mark = 1               ! collision flags in PHITS
         nbeta = 2
         pstep = 0d0

         np = 1

      return
      end
