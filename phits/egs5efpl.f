c

      subroutine egs5efpl(fpl,ein,icl,ich)
      use EGS5_MS_MOD !FURUTA20140825
      USE egs5_elecin_mod !<- 2015.08xx allocatable
      USE egs5_media_mod !<- 2015.08xx allocatable
      USE egs5_scpw_mod !<- 2015.08xx allocatable
      implicit none
!      save
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
      integer ircode,iarg
      logical go1

      real*8                                           ! Local variables
     * eie,                        ! Energy (total) of incident electron
     * de,ekef,
     * estepe,range,range0,dedx,dedx0,sig,sig0,scpow0,
     * ams,blcold,tmxs

      integer ierust,idr,lelec,irl,lelke,ib,nk1i,dok1s0

      save eie,de,ekef,estepe,range,range0,dedx0,
     * scpow0,blcold,tmxs
      save ierust,idr,lelec,irl,lelke,ib,nk1i,dok1s0
!$OMP THREADPRIVATE(eie,de,ekef,estepe,range,range0,dedx0)
!$OMP THREADPRIVATE(scpow0,blcold,tmxs)
!$OMP THREADPRIVATE(ierust,idr,lelec,irl,lelke,ib,nk1i,dok1s0)

      real*8
     * ustep0,kinit0,ktotal,detot,scpow,thard,tmscat,tinel,
     * hardstep,ecsda,k1s0

      real*8 EPSEMFP,ENEPS                            ! Local parameters

      data
     * EPSEMFP/1.E-12/,                    ! Smallest electron mfp value
     * ENEPS/0.0005/,     ! Difference between ecut and end-point energy
     * ierust/0/

      real *8 fpl,ein,u0,v0,w0,k1i,k1r,k1s
      integer icl,ich,mk,nbeta,i
      common /egs5cmn3/thard,tinel,tmscat,hardstep,sig,scpow,dedx,sig0,
     $                 ams
!$OMP THREADPRIVATE(/egs5cmn3/)
      common /egs5cmn4/k1i,k1r,k1s
!$OMP THREADPRIVATE(/egs5cmn4/)
      common /egs5cmn5/ircode
!$OMP THREADPRIVATE(/egs5cmn5/)
      common /paraj/ mstz(300), parz(300)
      real*8  parz
      integer mstz
      common /egs5cmn10/denstepold
      real*8 denstepold
!$OMP THREADPRIVATE(/egs5cmn10/)

      integer iwarning,iwcheck ! number of estepe warning & check parameter
      save iwarning,iwcheck
      data iwarning/0/
      data iwcheck/0/

      ielectr = ielectr + 1                ! Count entry into subroutine

      eke        =  ein
      eie        =  eke+RM
cc 2014/9/19 H.Iwase bugfix
      elke       = log(eke)
      np         =  1
      e(np)      =  eie
      ir(np)     =  icl
      irl        =  ir(np)
      lelec      =  ich
      iq(np)     =  lelec
      medium     =  med(irl)

      k1init(np) =  k1i
      k1rsd(np)  =  k1r
      k1step(np) =  k1s

            edep = 0d0


! keep previous run ----------------------------------------------
!     after collision
      if( thard.eq.0 .and. tinel.ne.0 .and. tmscat.ne.0)then
         deresid    = 0.d0
         deinitial  = 0.d0
         denstep    = 0.d0
         denstepold = 0.d0
         ircode = 1
         irold  = ir(np)
         irl    = ir(np)
         goto 1
      endif



! one of (tinel,tmscat,thard) is zero after events
      if( thard.ne.0 .and.  tinel.eq.0 )goto 3
      if( thard.ne.0 .and. tmscat.eq.0 )goto 4

! after crossing surface, none of (tinel,tmscat,thard) is zero
      if( tinel.ne.0 .and. tmscat.ne.0 .and. thard.ne.0) then

!!! same treatment tinel and tmscat here?

! T.Sato 2015/08/29, use starting scattering strengh again in new medium, only imsegs = 1
            if(medium.ne.medold.and.mstz(112).eq.1) then
             if(ek1s1(1,1).ne.0.d0 .or. ek1s0(1,1).ne.0.d0) then
              dok1s0 = 1
              nk1i = 0
             endif
            endif
! End of revision 2015/08/29 by T.Sato
            if(denstep.eq.0.d0) then  ! new medium or density
               go to 2
            else
               go to 3
            end if

                                                 ! (after collis go top)
      endif

! ----------------------------------------------------------------

! new electron ---------------------------------------------------
      if( thard.eq.0 .and. tinel.eq.0 .and. tmscat.eq.0)then
         k1step(np) = 0d0
         k1rsd(np)  = 0d0
         k1init(np) = 0d0
      endif

! ----------------------------------------------------------------
c

      deresid    = 0.d0
      deinitial  = 0.d0
      denstep    = 0.d0
      hardstep   = 0.d0
      denstepold = 0.d0

      dok1s0 = 0
      if(ircode.eq.-1) then
        if(ek1s1(1,1).ne.0.d0 .or. ek1s0(1,1).ne.0.d0) then
          dok1s0 = 1
          nk1i = 0
        endif
      endif
                                            ! --------------------------
      ircode = 1                            ! Set up for normal return
       irold = ir(np)                        ! Initialize previous region
         irl = ir(np)                          ! Region number
                                            ! --------------------------
1     continue                              ! Start of NEW-ELECTRON loop
                                            ! --------------------------
      if(mstz(85).eq.99) write(93,*)'--- 1 cont.'
           eie = e(np)
        medium = med(irl)

        !-->  Vacuum
        if (medium .eq.0 ) then
          tstep = vacdst
! T.Sato 2015/08/30, goto 999 instead of direct return
           goto 999
        end if

        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
        ! Top of tracking loop - compute parameters
        !++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

        !-->  re-enter here after a hard collision or to check cutoffs
2       continue
        if(mstz(85).eq.99) write(93,*)'--- 2 cont.'
        if(mstz(85).eq.99) write(93,*)'e,ecut',e(np),ecut(irl)


! T.Sato 2015/8/30, goto 999 instead of direct return
          if(e(np) .le. ecut(irl) .and. deresid.eq.0.d0) goto 999

          !--> Sample number of mfp's to transport before interacting
          if(hardstep .eq. 0.d0) then
            call randomset(rnnow)
            hardstep = max(-log(rnnow),EPSEMFP)
          end if

          !-->  re-enter loop here after an energy hinge
3         continue
          if(mstz(85).eq.99) write(93,*)'--- 3 cont.'

          rhof=rhor(irl)/rhom(medium)

          !-->  Get energy grid parameters
          elke = log(eke)
          lelke = eke1(medium)*elke + eke0(medium)

          !-->  Get stopping, scattering power, range
          if (lelec .eq. -1) then
            dedx0 = ededx1(lelke,medium)*elke + ededx0(lelke,medium)
            scpow0 = escpw1(lelke,medium)*elke + escpw0(lelke,medium)
            range0 = erang1(lelke,medium)*elke + erang0(lelke,medium)
            range0 = range0 - ectrng(irl)
          else
            dedx0 = pdedx1(lelke,medium)*elke + pdedx0(lelke,medium)
            scpow0 = pscpw1(lelke,medium)*elke + pscpw0(lelke,medium)
            range0 = prang1(lelke,medium)*elke + prang0(lelke,medium)
            range0 = range0 - pctrng(irl)
          end if
          dedx = rhof*dedx0
          scpow = rhof*scpow0
          range = rhof*range0 + ENEPS/dedx

          !-->  use current energy to get Moliere parameters
          !-->  and set up test for step size max check
          if(useGSD(medium).eq.0) then
            ems = e(np)
            bms = (ems - RM)*(ems + RM)/ems**2
            ams = blcc(medium)/bms
            gms = rhof*(xcc(medium)/(ems*bms))**2

            tmxs = 1.d0 / (log(ams/gms) * gms)
          end if

          !-->  get hard cross section
          call hardx(lelec,eke,lelke,elke,sig0)
          sig = rhof*sig0                       ! Density-ratio scaling

          if (sig .le. 0.) then
            thard = vacdst
          else
            thard = hardstep / sig
          end if

          !-->  Get a new energy loss step, set energy hinge distance
          if(denstep.eq.0.d0) then
            denstep = deresid
            estepe = estep1(lelke,medium)*elke + estep0(lelke,medium)
            if(estepe.lt.1.0d-10) then ! T.Sato 2016/5/18, avoid infinite loop, further revision on 2017/3/7
             if(iwcheck.eq.0.and.iwarning.le.100) then ! first time in a history
              write(*,*)'warning: estepe is smaller than 1.0d-10',estepe
              if(iwarning.eq.100) then ! too much warning, better to change dmax, no more warning
               write(*,*)'**** too much warning related to estepe *****'
               write(*,*)'It is better to decrease dmax(12) & dmax(13),'
               write(*,*)'otherwise calculation may take too much time'
               write(*,*)'*********************************************'
              endif
              iwarning=iwarning+1
             endif
             estepe=1.0d-4
             iwcheck=1 ! stop to write warning
            else
             iwcheck=0
            endif
            detot = eke * estepe
            !-->  allow region dependent scaling
            if(estepr(irl) .ne. 0) detot = detot * estepr(irl)
            !-->  if this step takes us below the cutoff, adjust
            if( (e(np)-detot) .lt. ecut(irl)) then
              detot = e(np)-ecut(irl)
            end if
            call randomset(rnnow)
            deinitial = rnnow * detot
            deresid = detot - deinitial
            denstep = denstep + deinitial

            if(mstz(85).eq.99) then ! debug output
               write(93,*)'--------------------------------new denstep'
               write(93,'(a10,f15.9)')'(rnnow)'  ,rnnow
               write(93,'(a10,f15.9)')'(denstep)',denstep
               write(93,'(a10,f15.9)')'(deinit)' ,deinitial
               write(93,'(a10,f15.9)')'(deresid)',deresid
               write(93,*)'--------------------------------'
            endif

          end if

          if(dedx.le.0.) then
            tinel = vacdst
            range = vacdst
          else
            tinel = denstep / dedx
          end if

          !-->  re-enter loop here after a multiple scatter
4         continue

      if(mstz(85).eq.99) then
         write(93,*)'--- 4 cont.'
         write(93,'(a10,f15.9)')'(e)',e(np)
      endif

          !-->  Get the scattering strength, sent the hinge length
          if(k1step(np) .eq. 0.d0) then
            k1step(np) = k1rsd(np)

            !-->  Get max scattering strength
            if(lelec .eq. -1) then
              kinit0 = ekini1(lelke,medium)*elke + ekini0(lelke,medium)
            else
              kinit0 = pkini1(lelke,medium)*elke + pkini0(lelke,medium)
            end if

            if(k1Lscl(irl).ne.0.d0) then
              kinit0 = kinit0 * (k1Lscl(irl) +  k1Hscl(irl) * elke)
            end if

            !-->  Get starting scattering strength
            if(dok1s0.eq.1) then
              if(lelec .eq. -1) then
                k1s0 = ek1s1(lelke,medium)*elke + ek1s0(lelke,medium)
              else
                k1s0 = pk1s1(lelke,medium)*elke + pk1s0(lelke,medium)
              end if
              k1s0 = k1s0*(2**nk1i)
              nk1i = nk1i + 1

              if(kinit0.gt.k1s0) then
                kinit0 = k1s0
              else
                dok1s0=0
              end if
            end if

            ktotal = rhof*kinit0

            if(useGSD(medium).eq.0) then
              !->  make sure total K1 is less than that of tmxs
              if(ktotal/scpow.gt.tmxs) then
                itmxs = itmxs + 1
                if(tmxset) ktotal = tmxs*scpow
              !-> make sure that kinit gives us a valid omega0.
              !-> use 2.80 instead of e because the increase in scpow
              !-> as particle slows means tmstep will be < tmscat.
              else
                omega0 = ams * ktotal/scpow
                if(omega0 .lt. 2.80) then
                   ktotal = scpow * 2.80 / ams
                end if
              end if
            end if

            call randomset(rnnow)
            k1init(np) = rnnow * ktotal
            k1rsd(np) = ktotal - k1init(np)
            k1step(np) = k1step(np) + k1init(np)

          endif

          if (scpow .le. 0.) then
            tmscat = vacdst
          else
            tmscat = k1step(np) / scpow
          end if

          !-->  re-enter loop here after a boundary crossing or B
          !-->  field step was taken
5         continue
          if(mstz(85).eq.99) write(93,*)'--- 5 cont.'

          tstep = MIN(tmscat,tinel,thard)

          if(mstz(85).eq.99)then ! debug output
             write(93,*)'----------------------------------steps'
             write(93,'(a10,f15.9)')'(e)'     ,e(np)
             write(93,'(a10,f15.9)')'(thard)' ,thard
             write(93,'(a10,f15.9)')'(tinel)' ,tinel ! dedx
             write(93,'(a10,f15.9)')'(tmscat)',tmscat
             write(93,*)'----------------------------------'
          endif

          fpl = tstep

          k1i =    k1init(np)
          k1r =    k1rsd(np)
          k1s =    k1step(np)

! T.Sato 2015/08/30, keep old medium
  999 medold=medium

      return

      end
