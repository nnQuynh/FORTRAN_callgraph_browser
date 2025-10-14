

      subroutine egs5ecoll(e0,u0,v0,w0,wt0,lelec,irl,mark,nbeta)
      use EGS5_MS_MOD !FURUTA20140825
      USE egs5_media_mod !<- 2015.08xx allocatable
      implicit none
!      save
      include 'err.inc'
      include 'include/egs5_h.f'               ! Main EGS5 "header" file

      include 'include/egs5_bounds.f'    ! COMMONs required by EGS5 code
      include 'include/egs5_epcont.f'
      include 'include/egs5_misc.f'
      include 'include/egs5_mults.f'
      include 'include/egs5_stack.f'
      include 'include/egs5_uphiin.f'
      include 'include/egs5_uphiot.f'
      include 'include/egs5_useful.f'
      include 'include/egs5_userpr.f'
      include 'include/egs5_usersc.f'

      include 'include/counters.f'       ! Additional (non-EGS5) COMMONs

      real*8 rnnow                                           ! Arguments
      integer iarg
      logical go1

      real*8                                           ! Local variables
     * eie,                        ! Energy (total) of incident electron
     * de,ekef,
     * estepe,range,range0,dedx,dedx0,sig,sig0,scpow0,
     * ams,blcold,tmxs

      integer ierust,idr,lelec,irl,lelke,ib,nk1i,dok1s0

      real*8
     * ustep0,kinit0,ktotal,detot,scpow,thard,tmscat,tinel,
     * hardstep,ecsda,k1s0

      real*8 EPSEMFP,ENEPS                            ! Local parameters

      data
     * EPSEMFP/1.E-12/,                    ! Smallest electron mfp value
     * ENEPS/0.0005/,     ! Difference between ecut and end-point energy
     * ierust/0/

      common /egs5cmn3/thard,tinel,tmscat,hardstep,sig,scpow,dedx,sig0,
     $                 ams
!$OMP THREADPRIVATE(/egs5cmn3/)

      integer nnn,nomp,nclsts,iclusts,jclusts
      include 'param00.inc'

      real*8 qclusts
      common /clustt/ nclsts, iclusts(nnn)
!$OMP THREADPRIVATE(/clustt/)
      common /clustw/ jclusts(0:8,nnn),  qclusts(0:12,nnn)
!$OMP THREADPRIVATE(/clustw/)

      integer kf,itype,id,charge
      real*8 rms,rmg,pr,etot,ekine,pxl,pyl,pzl

      integer mark,nbeta,k,i
      real*8 e0,x0,y0,z0,u0,v0,w0,wt0,k1i,k1r,k1s
      real*8 fdummy
      common /egs5cmn4/k1i,k1r,k1s
!$OMP THREADPRIVATE(/egs5cmn4/)
      integer ircode  !<- move here
      common /egs5cmn5/ircode
!$OMP THREADPRIVATE(/egs5cmn5/)

      integer ncpls

      integer numpal
      real*8 rumpal
      common /clustl/ rumpal(0:20), numpal(0:20)
!$OMP THREADPRIVATE(/clustl/)

      common /anglth/ cosph,costh,sinph,sinth
      real*8 cosph,costh,sinph,sinth,ein
!$OMP THREADPRIVATE(/anglth/)
      common /paraj/ mstz(300), parz(300)
      real*8  parz
      integer mstz
      common /egs5cmn10/denstepold,de2
      real*8 denstepold,de2
!$OMP THREADPRIVATE(/egs5cmn10/)


      np         =  1
      eke        =  e0
      eie        =  eke + RM

      e(np)      =  eie


      u(np)      =  u0
      v(np)      =  v0
      w(np)      =  w0

      wt(np)     =  wt0
      iq(np)     =  lelec
      nclsts     =  0
      do i=0,20
       numpal(i)=0
       rumpal(i)=0
      enddo

      medium     =  med(irl)

      k1init(np) = k1i
      k1rsd(np)  = k1r
      k1step(np) = k1s



       if(tmscat .eq. 0.0) then

         tmstep = (k1init(np) + k1rsd(np)) / scpow

         if(useGSD(medium).eq.0) then
            omega0 = ams*tmstep*rhof
            if (omega0 .le. 2.718282) then
               iskpms = 1
            else
               iskpms = 0
               blc = log(omega0)
               blcold = blc
               if (blc .lt. 1.306853) then
                  b = -10.27666 + blc*(17.82596 - 6.468813*blc)
               else
                  ib = b0bgb + blc*b1bgb
                  if (ib .gt. nbgb) then
                     write(6,101) ib
101                 FORMAT('electr warning: IB > NBGB =',I5,' set to 8')
                     ib = nbgb
                  end if
                  b = bgb0(ib) + blc*(bgb1(ib) + blc*bgb2(ib))
               end if
            end if
         end if


         call mscat
         call uphi(2,1)         ! Set direction cosines
         goto 99


      else if(thard .eq. 0.0) then !  hard collision

            e(np) = e(np) - deinitial + denstep

       if(e(np)-RM.lt.1.0d-3) then ! T.Sato 2023/07/16 Bugfix for electromagnetic field
        e(np)=RM+1.001d-3 ! kinetic energy should be slightly above 1 keV
       endif

            if(mstz(85).eq.99) then
               write(93,*)'===== HARDCOL ====='
               write(93,'(a10,f15.9)')'(e)',e(np)
            endif


            call egs5collis(lelec,irl,sig0,go1)

            k1i = k1init(1)
            k1r = k1rsd(1)
            k1s = k1step(1)

! energy cut treatment
            if(abs(iq(np)).eq.1 .and. e(np).le.ecut(irl))then
               thard  = 0d0
               tinel  = 0d0
               tmscat = 0d0

            goto 99
            endif

            if (iq(np) .eq. 0) then
               thard = 0.d0
               tinel = 0.d0
              tmscat = 0.d0

              goto 99

             end if


          endif

 99    continue ! move EGS5 particle informations into PHITS


           ncpls = 1

       call egs2phits(ncpls)

      costh = 1.0
      sinth = 0.0
      cosph = 1.0
      sinph = 0.0


         mark = 1               ! collision flags in PHITS
         nbeta = 2

         k1i = k1init(np)
         k1r = k1rsd(np)
         k1s = k1step(np)

         np = 1
         ircode = 2

        return

      !++++++++++++++++++++++++++++++++++++++++
      ! USER-REQUESTED ELECTRON DISCARD SECTION
      !++++++++++++++++++++++++++++++++++++++++
 16   continue

      ! adjust energy for mid-hinge discard case
      e(np) = e(np) - deinitial + denstep
      deinitial  = 0.d0
      denstep    = 0.d0
      deresid    = 0.d0
      denstepold = 0.d0

      idisc = abs(idisc)

      if (lelec .eq. -1 .or. idisc .eq. 99) then
        edep = e(np) - RM
      else
        edep = e(np) + RM     ! positron escape
      end if

      iarg = 3                   ! User requested discard
!                                =================
      if (iausfl(iarg+1) .ne. 0) call ausgab(iarg)
!                                =================

                                                      ! ----------------
      return                                          ! Return to SHOWER
                                                      ! ----------------
      end


