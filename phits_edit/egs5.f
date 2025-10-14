      subroutine ausgab(iarg)

      implicit none

      include 'include/egs5_h.f'                ! Main EGS "header" file

      include 'include/egs5_stack.f'     ! COMMONs required by EGS5 code
      include 'include/egs5_useful.f'

      integer iarg                                          ! Arguments

      return
      end







      subroutine howfar

      implicit none

      include 'include/egs5_h.f'                ! Main EGS "header" file

      include 'include/egs5_epcont.f'    ! COMMONs required by EGS5 code
      include 'include/egs5_stack.f'

      end







      subroutine counters_out(ioflag)

      implicit none

      include 'include/counters.f'

      include 'err.inc'

      integer ioflag                                         ! Arguments


      common /paran/ icfn(100), ilfn(100), chfn(100)
      character chfn*200
      integer icfn, ilfn

      if (ioflag.eq.0) then       ! Open unit=99 and initialize counters

        iannih           = 0
        iaphi            = 0
        ibhabha          = 0
        ibrems           = 0
        icollis          = 0
        icompt           = 0
        iedgbin          = 0
        ieii             = 0
        ielectr          = 0
        ihatch           = 0
        ihardx           = 0
        ikauger          = 0
        ikshell          = 0
        ikxray           = 0
        ilauger          = 0
        ilshell          = 0
        ilxray           = 0
        imoller          = 0
        imscat           = 0
        ipair            = 0
        iphoto           = 0
        iphoton          = 0
        iraylei          = 0
        ishower          = 0
        iuphi            = 0
        itmxs            = 0
        noscat           = 0
        iblock           = 0

      else if (ioflag.eq.1) then
        iannih_sum           = 0
        iaphi_sum            = 0
        ibhabha_sum          = 0
        ibrems_sum           = 0
        icollis_sum          = 0
        icompt_sum           = 0
        iedgbin_sum          = 0
        ieii_sum             = 0
        ielectr_sum          = 0
        ihatch_sum           = 0
        ihardx_sum           = 0
        ikauger_sum          = 0
        ikshell_sum          = 0
        ikxray_sum           = 0
        ilauger_sum          = 0
        ilshell_sum          = 0
        ilxray_sum           = 0
        imoller_sum          = 0
        imscat_sum           = 0
        ipair_sum            = 0
        iphoto_sum           = 0
        iphoton_sum          = 0
        iraylei_sum          = 0
        ishower_sum          = 0
        iuphi_sum            = 0
        itmxs_sum            = 0
        noscat_sum           = 0
        iblock_sum           = 0

        iannih_sum  = iannih_sum  + iannih
        iaphi_sum   = iaphi_sum   + iaphi
        ibhabha_sum = ibhabha_sum + ibhabha
        ibrems_sum  = ibrems_sum  + ibrems
        icollis_sum = icollis_sum + icollis
        icompt_sum  = icompt_sum  + icompt
        iedgbin_sum = iedgbin_sum + iedgbin
        ieii_sum    = ieii_sum    + ieii
        ielectr_sum = ielectr_sum + ielectr
        ihatch_sum  = ihatch_sum  + ihatch
        ihardx_sum  = ihardx_sum  + ihardx
        ikauger_sum = ikauger_sum + ikauger
        ikshell_sum = ikshell_sum + ikshell
        ikxray_sum  = ikxray_sum  + ikxray
        ilauger_sum = ilauger_sum + ilauger
        ilshell_sum = ilshell_sum + ilshell
        ilxray_sum  = ilxray_sum  + ilxray
        imoller_sum = imoller_sum + imoller
        imscat_sum  = imscat_sum  + imscat
        ipair_sum   = ipair_sum   + ipair
        iphoto_sum  = iphoto_sum  + iphoto
        iphoton_sum = iphoton_sum + iphoton
        iraylei_sum = iraylei_sum + iraylei
        ishower_sum = ishower_sum + ishower
        iuphi_sum   = iuphi_sum   + iuphi
        itmxs_sum   = itmxs_sum   + itmxs
        noscat_sum  = noscat_sum  + noscat
        iblock_sum  = iblock_sum  + iblock



        OPEN(UNIT=599,FILE=chfn(23)(1:ilfn(23))//'job.out99',&
            STATUS='unknown')
        write(599,*) 'Values for subroutine-entry counters:'
        write(599,*)
        write(599,*) 'iannih =', iannih_sum !<- iannih
        write(599,*) 'iaphi  =', iaphi_sum !<- iaphi
        write(599,*) 'ibhabha=', ibhabha_sum !<- ibhabha
        write(599,*) 'ibrems =', ibrems_sum !<- ibrems
        write(599,*) 'icollis=', icollis_sum !<- icollis
        write(599,*) 'icompt =', icompt_sum !<- icompt
        write(599,*) 'iedgbin=', iedgbin_sum !<- iedgbin
        write(599,*) 'ieii   =', ieii_sum !<- ieii
        write(599,*) 'ielectr=', ielectr_sum !<- ielectr
        write(599,*) 'ihatch =', ihatch_sum !<- ihatch
        write(599,*) 'ihardx =', ihardx_sum !<- ihardx
        write(599,*) 'ikauger=', ikauger_sum !<- ikauger
        write(599,*) 'ikshell=', ikshell_sum !<- ikshell
        write(599,*) 'ikxray =', ikxray_sum !<- ikxray
        write(599,*) 'ilauger=', ilauger_sum !<- ilauger
        write(599,*) 'ilshell=', ilshell_sum !<- ilshell
        write(599,*) 'ilxray =', ilxray_sum !<- ilxray
        write(599,*) 'imoller=', imoller_sum !<- ilmoller
        write(599,*) 'imscat =', imscat_sum !<- imscat
        write(599,*) 'ipair  =', ipair_sum !<- ipair
        write(599,*) 'iphoto =', iphoto_sum !<- iphoto
        write(599,*) 'iphoton=', iphoton_sum !<- iphoton
        write(599,*) 'iraylei=', iraylei_sum !<- iraylei
        write(599,*) 'ishower=', ishower_sum !<- ishower
        write(599,*) 'iuphi  =', iuphi_sum !<- iuphi
        write(599,*) 'itmxs  =', itmxs_sum !<- itmxs
        write(599,*) 'noscat =', noscat_sum !<- noscat
        write(599,*) 'iblock =', iblock_sum !<- iblock
        CLOSE(UNIT=599)

      else
        write(506,*) ' *** Error using subroutine counters_out'
        write(506,*) '       ioflag=',ioflag
         ErrCha = ''
         ErrID = 'L:185/R:counters_out/F:egs5.f' !E86_001_001
         call ErrWrite(ErrID,ErrCha)

        write(*,*) ' *** Error using subroutine counters_out'
        write(*,*) '       ioflag=',ioflag
        CLOSE(UNIT=599)
        write(506,*) '     Program stopped'
        write(*,*) '     Program stopped'
        stop
      end if

      return
      end








      subroutine annih

      implicit none

      include 'include/egs5_h.f'               ! Main EGS5 "header" file

      include 'include/egs5_stack.f'     ! COMMONs required by EGS5 code
      include 'include/egs5_uphiot.f'
      include 'include/egs5_useful.f'

      include 'include/counters.f'       ! Additional (non-EGS5) COMMONs

      real*8 rnnow                                           ! Arguments


      common /tstara/ atmrc(10,8)      ! Ogawa 2023/7/18, (10,7) -> (10,8) for plasmon

      real*8 atmrc

      iannih = iannih + 1                  ! Count entry into subroutine

      avip = e(np) + RM
      a    = avip/RM
      g    = a - 1.0
      t    = g - 1.0
      p    = sqrt(a*t)
      pot  = p/t
      ep0  = 1.0/(a + p)

1     continue                    ! Sample 1/ep from ep = ep0 to 1 - ep0
        call randomset(rnnow)
        ep = ep0*exp(rnnow*log((1.0 - ep0)/ep0))

                                       ! Decide whether or not to accept
        call randomset(rnnow)
        rejf = 1.0 - ep + (2.0*g - 1.0/ep)/a**2
        if (rnnow .gt. rejf) go to 1

                 ! This completes sampling of a distribution which is
                 ! asymmetric about ep = 1/2, but which when symmetrized
                 ! is the symmetric annihilation distribution.

                                  ! Set up energies, place them on stack
      ep        = max(ep,1.D0 - ep)          ! Pick ep in (1/2, 1 - ep0)
      esg1      = avip*ep                 ! Energy of secondary gamma #1
      e(np)     = esg1           ! Place energy of gamma #1 on the stack
      esg2      = avip - esg1             ! Energy of secondary gamma #2
      e(np + 1) = esg2           ! Place energy of gamma #2 on the stack
      iq(np)    = 0                          ! Make particle np a photon

                             ! Set up angles for the higher energy gamma
      costhe = (esg1 - RM)*pot/esg1
      costhe = min(1.D0,costhe)                      ! Fix (860724/dwor)
      sinthe = sqrt((1.0 - costhe)*(1.0 + costhe))
      call uphi(2,1)                             ! Set direction cosines

                                             ! Set up lower energy gamma
      np     = np + 1                          ! Increase the stack size
      iq(np) = 0                                      ! Make it a photon
      costhe = (esg2 - RM)*pot/esg2
      costhe = min(1.D0,costhe)                      ! Fix (860724/dwor)
      sinthe = - sqrt((1.0 - costhe)*(1.0 + costhe))
      call uphi(3,2)                             ! Set direction cosines



      t = e(np)
      e(np) = e(np-1)
      e(np-1) = t
      t = u(np)
      u(np) = u(np-1)
      u(np-1) = t
      t = v(np)
      v(np) = v(np-1)
      v(np-1) = t
      t = w(np)
      w(np) = w(np-1)
      w(np-1) = t

      call npconv

      call npconv
      atmrc(3,5) = atmrc(3,5) + 1.d0
                                                      ! ----------------
      return                                          ! Return to ELECTR
                                                      ! ----------------
      end







      subroutine aphi(br)

      implicit none

      include 'include/egs5_h.f'               ! Main EGS5 "header" file

      include 'include/egs5_epcont.f'    ! COMMONs required by EGS5 code
      include 'include/egs5_stack.f'
      include 'include/egs5_uphiin.f'
      include 'include/egs5_uphiot.f'

      include 'include/counters.f'       ! Additional (non-EGS5) COMMONs

      real*8 br,rnnow                                        ! Arguments
      integer iarg

      integer ldpola

      iaphi = iaphi + 1                    ! Count entry into subroutine

      iarg = 21

      if (iausfl(iarg+1) .ne. 0) call ausgab(iarg)


      pnorm0 = sqrt(u(np)*u(np) + v(np)*v(np) + w(np)*w(np))
      u(np) = u(np)/pnorm0
      v(np) = v(np)/pnorm0
      w(np) = w(np)/pnorm0
      valmax = br + 1./br

 1    continue
        call randomset(rnnow)
        ph0 = rnnow*twopi
        sinph0 = sin(ph0)
        cph0 = pi5d2 - ph0
        cosph0 = sin(cph0)
        valloc = sqrt(sinph0*sinph0 + cosph0*cosph0)
        sinph0 = sinph0/valloc
        cosph0 = cosph0/valloc
        val = (valmax - 2.*sinthe*sinthe*cosph0*cosph0)/valmax
        call randomset(rnnow)
        if (rnnow .le. val) go to 2
      go to 1

 2    continue
      anorm2 =costhe*costhe*cosph0*cosph0 + sinph0*sinph0
      call randomset(rnnow)


      if (ldpola .eq. 1) then
        call randomset(rnnow)
        eta = rnnow*twopi
        sineta = sin(eta)
        ceta = pi5d2 - eta
        coseta = sin(ceta)
      else
        anormr = 1./sqrt(anorm2)
        sineta = -anormr*sinph0
        coseta = anormr*costhe*cosph0
      end if

      ufa = costhe*cosph0*coseta - sinph0*sineta
      vfa = costhe*sinph0*coseta + cosph0*sineta
      wfa = -sinthe*coseta

      asav = u(np)
      bsav = v(np)
      csav = w(np)

      sinps2 = asav*asav + bsav*bsav
      if (sinps2 .lt. 1.E-20) then
        cosomg = uf(np)
        sinomg = vf(np)
      else
        sinpsi = sqrt(sinps2)
        sindel = bsav/sinpsi
        cosdel = asav/sinpsi
        cosomg = cosdel*csav*uf(np) + sindel*csav*vf(np) - sinpsi*wf(np)
        sinomg = -sindel*uf(np) + cosdel*vf(np)
      end if

      enorm = sqrt(uf(np)*uf(np) + vf(np)*vf(np) + wf(np)*wf(np))
      if (enorm .lt. 1.E-4) then
        call randomset(rnnow)
        omg = rnnow*twopi
        sinomg = sin(omg)
        comg = pi5d2 - omg
        cosomg = sin(comg)
      end if

      cosphi = cosomg*cosph0 - sinomg*sinph0
      sinphi = sinomg*cosph0 + cosomg*sinph0

      ufb = cosomg*ufa - sinomg*vfa
      vfb = sinomg*ufa + cosomg*vfa
      wfb = wfa

      if (sinps2 .lt. 1.E-20) then
        uf(np) = ufb
        vf(np) = vfb
        wf(np) = wfb
      else
        uf(np) = cosdel*csav*ufb - sindel*vfb + asav*wfb
        vf(np) = sindel*csav*ufb + cosdel*vfb + bsav*wfb
        wf(np) = -sinpsi*ufb + csav*wfb
      end if
      enorm = sqrt(uf(np)*uf(np) + vf(np)*vf(np) + wf(np)*wf(np))
      uf(np) = uf(np)/enorm
      vf(np) = vf(np)/enorm
      wf(np) = wf(np)/enorm

      return

      end









      subroutine bhabha

      USE egs5_thresh_mod

      implicit none

      include 'include/egs5_h.f'               ! Main EGS5 "header" file


      include 'include/egs5_stack.f'     ! COMMONs required by EGS5 code
      include 'include/egs5_uphiot.f'
      include 'include/egs5_useful.f'

      include 'include/counters.f'       ! Additional (non-EGS5) COMMONs

      real*8 rnnow                                           ! Arguments



      common /tstara/ atmrc(10,8)      ! Ogawa 2023/7/18, (10,7) -> (10,8) for plasmon

      real*8 atmrc

      ibhabha = ibhabha + 1                ! Count entry into subroutine

      eip   = e(np)
      ekin  = eip - RM
      t0    = ekin/RM
      e0    = t0 + 1.0
      yy    = 1.0/(t0 + 2.0)
      e02   = e0*e0
      betai2= e02/(e02-1.0)
      ep0   = te(medium)/ekin
      ep0c  = 1.0 - ep0
      y2    = yy*yy
      yp    = 1.0 - 2.0*yy
      yp2   = yp*yp
      b4    = yp2*yp
      b3    = b4 + yp2
      b2    = yp*(3.0 + y2)
      b1    = 2.0 - y2
      br    = ep0
      re1   = ep0c*(betai2-br*(b1-br*(b2-br*(b3-br*b4))))
      br    = 1.0
      re2   = ep0c*(betai2-br*(b1-br*(b2-br*(b3-br*b4))))
      remax = max(re1,re2)

                                            ! Sample/reject to obtain br
1     continue
        call randomset(rnnow)
        br = ep0/(1.0 - ep0c*rnnow)

                                       ! Decide whether or not to accept
        call randomset(rnnow)
        rejf = ep0c*(betai2-br*(b1-br*(b2-br*(b3-br*b4))))
        rejf = rejf/remax

        if (rnnow .gt.rejf) go to 1

                            ! If e- got more energy than e+, move the e+
                            ! pointer and reflect br (this puts e+ on
                            ! top of stack if it has less energy)
      if (br .lt. 0.5) then
        iq(np+1) = -1
        k1step(np+1) = 0.
        k1init(np+1) = 0.
        k1rsd(np+1) = 0.

      else
        iq(np) = -1
        iq(np+1) = 1
        br = 1.0 - br
        k1step(np+1) = k1step(np)
        k1init(np+1) = k1init(np)
        k1rsd(np+1) = k1rsd(np)
        k1step(np) = 0.
        k1init(np) = 0.
        k1rsd(np) = 0.
      end if

                                                  ! Divide up the energy
      br = max(br,0.D0)    ! Avoids possible negative energy (round off)
      ekse2 = br*ekin          ! Kinetic energy of secondary electron #2
      ese1 = eip - ekse2         ! Total energy of secondary electron #1
      ese2 = ekse2 + RM          ! Total energy of secondary electron #2
      e(np) = ese1
      e(np+1) = ese2
                                           ! Bhabha angles are uniquely
                                           ! determined by kinematics
      h1 = (eip + RM)/ekin
      dcosth = h1*(ese1 - RM)/(ese1 + RM)
      ttt = 1.0 - dcosth
      if (ttt.le.0.0) then
        sinthe = 0.0
      else
        sinthe = sqrt(1.0 - dcosth)
      end if
      costhe = sqrt(dcosth)
      call uphi(2,1)                             ! Set direction cosines
      np = np + 1
      dcosth = h1*(ese2 - RM)/(ese2 + RM)
      ttt = 1.0 - dcosth
      if (ttt.le.0.0) then
        sinthe = 0.0
      else
        sinthe = -sqrt(1.0 - dcosth)
      end if
      costhe = sqrt(dcosth)
      call uphi(3,2)                             ! Set direction cosines


      if(iq(np).eq.1) then
        iq(np) = -1
        iq(np-1) = 1
        t=e(np)
        e(np) = e(np-1)
        e(np-1) = t
        t = u(np)
        u(np) = u(np-1)
        u(np-1) = t
        t = v(np)
        v(np) = v(np-1)
        v(np-1) = t
        t = w(np)
        w(np) = w(np-1)
        w(np-1) = t
        t=k1step(np)
        k1step(np) = k1step(np-1)
        k1step(np-1) = t
        t=k1init(np)
        k1init(np) = k1init(np-1)
        k1init(np-1) = t
        t=k1rsd(np)
        k1rsd(np) = k1rsd(np-1)
        k1rsd(np-1) = t
      end if


      call npconv
      atmrc(3,4) = atmrc(3,4) + 1.d0
                                                      ! ----------------
      return                                          ! Return to ELECTR
                                                      ! ----------------
      end






      block data egs5

      implicit none

      include 'include/egs5_h.f'               ! Main EGS5 "header" file

      include 'include/egs5_epcont.f'
      include 'include/egs5_mults.f'
      include 'include/egs5_uphiot.f'
      include 'include/egs5_useful.f'
      include 'include/randomm.f'


      end


      block data egs5_2

      include 'include/egs5_h.f'
      include 'include/egs5_stack.f'
      include 'include/egs5_uservr.f'
      end








      subroutine block_set
      use EGS5_MS_MOD !FURUTA20140825

      USE egs5_brempr_mod
      USE egs5_eiicom_mod
      USE egs5_edge_mod
      USE egs5_mscon_mod

      implicit none

      include 'include/egs5_h.f'               ! Main EGS5 "header" file






      include 'include/egs5_bounds.f'    ! COMMONs required by EGS5 code
      include 'include/egs5_misc.f'
      include 'include/egs5_usersc.f'
      include 'include/egs5_userxt.f'
      include 'include/counters.f'  !<- counters are private


      integer  mxmat,mxmat0,mxnel
      common /kmat1a/ mxmat, mxmat0, mxnel
      integer i

      iblock = iblock + 1 !<-  ada.2014.12 this is " threadprivate "


      vacdst = 1.d8
      do i = 1, MXREG
        ecut(i) = 0.d0
        pcut(i) = 0.d0
      end do


      kmpi = 512
      kmpo = 508
      dunit = 1.d0

      do i = 1, MXREG
        med(i) = 1
        rhor(i) = 0.d0
        k1Hscl(i) = 0.d0
        k1Lscl(i) = 0.d0
        ectrng(i) = 0.d0
        pctrng(i) = 0.d0
        nomsct(i) = 0
      end do
      med(1) = 0


      tmxset = .true.


      k1mine = 1.d30
      k1maxe = -1.d30
      k1minp = 1.d30
      k1maxp = -1.d30


      emaxe = 0.d0
      do i = 1, MXREG
        estepr(i) = 0.d0
        esave(i) = 0.d0
      end do


      do i = 1, MXREG
        iphter(i) = 0
      end do

      return
      end









      subroutine brems

      USE egs5_brempr_mod !<- 2015.08xx allocatable
      USE egs5_thresh_mod !<- 2015.08xx allocatable

      implicit none

      include 'include/egs5_h.f'               ! Main EGS5 "header" file


      include 'include/egs5_stack.f'
      include 'include/egs5_uphiot.f'
      include 'include/egs5_useful.f'

      include 'include/counters.f'       ! Additional (non-EGS5) COMMONs

      real*8 rnnow,rnnow1,rnnow2                             ! Arguments
      integer lvx,lvl0,lvl,idistr,i

      real*8 AILN2,AI2LN2                             ! Local parameters
     
      common /tstara/ atmrc(10,8)      ! Ogawa 2023/7/18, (10,7) -> (10,8) for plasmon

      real*8 atmrc
      integer ichkelc

      ibrems = ibrems + 1                  ! Count entry into subroutine

      if(latch(np).ne.0) then
        write(506,*) 'e(np) at the top of brems=',e(np)
      end if

      eie = e(np)
      ichkelc = iq(np)   ! S.Abe 2017/02/08
      np = np + 1
      if (eie .lt. 50.0) then        ! Choose Bethe-Heitler distribution
        lvx = 1
        lvl0 = 0
      else         ! Choose Coulomb-corrected Bethe Heitler distribution
        lvx = 2
        lvl0 = 3
      end if
      abrems = float(int(AILN2*log(eie/ap(medium))))

                                 ! Start of main sampling-rejection loop
1     continue
        call randomset(rnnow)

                                              ! Start of (1-br)/br
                                              ! subdistribution sampling
        if (0.5 .lt. (abrems*alphi(lvx,medium) + 0.5)*rnnow) then
          call randomset(rnnow)
          idistr = abrems*rnnow
          p = pwr2i(idistr+1)
          lvl = lvl0 + 1
          call randomset(rnnow)
          if (rnnow .ge. AI2LN2) then
2           continue
              call randomset(rnnow)
              call randomset(rnnow1)
              call randomset(rnnow2)
              h = max(rnnow1,rnnow2)
              br = 1.0 - 0.5*h
              if (rnnow .gt. 0.5/br) go to 2
          else
            call randomset(rnnow)
            br = rnnow*0.5
          end if
          br = br*p
                                 ! Start of 2br subdistribution sampling
        else
          call randomset(rnnow1)
          call randomset(rnnow2)
          br = max(rnnow1,rnnow2)
          lvl = lvl0 + 2
        end if

        esg = eie*br
        if (esg .lt. ap(medium)) go to 1

        ese = eie - esg
        if (ese .lt. RM) go to 1
        del = br/ese
                                                     ! Check that Adelta
                                                     ! and Bdelta > 0
        if (del .ge. delpos(lvx,medium)) go to 1
        delta = delcm(medium)*del

        call randomset(rnnow)
        if (rnnow .gt. rejf) go to 1
                                                     ! Set up new photon
      if (ibrdst .ne. 1) then             ! Polar angle is m/E (default)
        theta = RM/eie
      else                                          ! Sample polar angle
        ztarg = zbrang(medium)   ! iwase i.c.o. iprdst not set ztarg wrong
        tteie = eie/RM
        ttese = ese/RM
        esedei = ttese/tteie
        y2max = (PI*tteie)**2d0
        rjarg1 = 1.0d0 + esedei**2d0
        rjarg2 = 3.0d0*rjarg1 - 2.d0*esedei
        rjarg3 = ((1.0d0 - esedei)/(2.d0*tteie*esedei))**2d0
        y2tst1 = (1.d0 + 0.d0)**2d0
       
5       continue
          call randomset(rnnow)
          y2tst =rnnow/(1.0 - rnnow + 1.0/y2max)
          y2tst1 = (1.0 + y2tst)**2
        
          call randomset(rnnow)

          if (rnnow .gt. (rejtst/rejtop)) go to 5
        theta = sqrt(y2tst)/tteie
      end if

      call uphi(1,3)                  ! Set direction cosines for photon
                            ! Put lowest energy particle on top of stack

        iq(np) = 0
        e(np) = esg
        e(np-1) = ese
        k1step(np) = 0.
        k1init(np) = 0.
        k1rsd(np) = 0.


      call npconv
      if( ichkelc .eq. -1) atmrc(2,3) = atmrc(2,3) + 1.d0
      if( ichkelc .eq.  1) atmrc(3,3) = atmrc(3,3) + 1.d0
                                                      ! ----------------
      return                                          ! Return to ELECTR
                                                      ! ----------------
      end








      subroutine egs5collis(lelec,irl,sig0,go1)

      USE egs5_brempr_mod
      USE egs5_elecin_mod
      USE egs5_media_mod
      USE egs5_edge_mod
      USE egs5_thresh_mod

      implicit none

      include 'include/egs5_h.f'               ! Main EGS5 "header" file


      include 'include/egs5_bounds.f'    ! COMMONs required by EGS5 code
      include 'include/egs5_epcont.f'
      include 'include/egs5_misc.f'
      include 'include/egs5_stack.f'
      include 'include/egs5_uphiin.f'
      include 'include/egs5_uphiot.f'
      include 'include/egs5_useful.f'
      include 'include/egs5_userxt.f'

      include 'include/counters.f'       ! Additional (non-EGS5) COMMONs

      include 'err.inc'

      real*8 rnnow                                           ! Arguments
      integer iarg

      integer lelec
      integer irl,i
      logical go1
      real*8 sig0



      real*8 ebr1, pbr1,pbr2,frstbr,fdummy
      integer idummy,npstrt,icsplt,lelke

      common /paraj/ mstz(300), parz(300)
      real*8  parz
      integer mstz

      icollis = icollis + 1                ! Count entry into subroutine

      go1 = .false.

      eke = e(np) - RM
      elke = log(eke)
      lelke = eke1(medium)*elke + eke0(medium)





1     continue
      if (lelec .lt. 0) then                                      ! e-

        if (eke .le. ap(med(irl))) ebr1 = 0.
        call randomset(rnnow)
        if (rnnow .lt. ebr1) then                              ! Brems
         go to 9
        else                                         ! Probably Moller
          if (e(np) .le. thmoll(medium)) then             ! Not Moller
            if (ebr1 .le. 0.) then                  ! Not Brems either
              go1 = .true.
              return
            end if
            go to 9                                     ! Forced Brems
          end if

          iarg = 8                                ! BEFORE call moller

          if (iausfl(iarg+1) .ne. 0) call ausgab(iarg)


          if(mstz(85).eq.99) write(93,*)'----- MOLLER -----' ! debug output
          call moller         ! To determine energies and polar angles


          iarg = 9                                 ! AFTER call moller

          if (iausfl(iarg+1) .ne. 0) call ausgab(iarg)

          if(iq(np) .eq. 0) then
            return              ! KEK addition to prevent EII
                                ! K-Xray from being discarded
          endif                 ! in ELECTRA

        end if
        go1 = .true.
        return
      end if




      pbr1 = pbr11(lelke,medium)*elke + pbr10(lelke,medium)
      if (eke .le. ap(med(irl))) pbr1 = 0.
      call randomset(rnnow)
      if (rnnow .lt. pbr1) then                                ! Brems
        go to 9
      end if



      pbr2 = pbr21(lelke,medium)*elke + pbr20(lelke,medium)
      if (rnnow .lt. pbr2) then                               ! Bhabha

        iarg = 10                                 ! BEFORE call bhabha

        if (iausfl(iarg+1) .ne. 0) call ausgab(iarg)


        call bhabha           ! To determine energies and polar angles


        iarg = 11                                  ! AFTER call bhabha

        if (iausfl(iarg+1) .ne. 0) call ausgab(iarg)


      else                                    ! Annihilation-in-flight

        iarg = 12                                  ! BEFORE call annih

        if (iausfl(iarg+1) .ne. 0) call ausgab(iarg)


        call annih            ! To determine energies and polar angles


        uf(np) = 0.
        uf(np+1) = 0.
        vf(np) = 0.
        vf(np+1) = 0.
        wf(np) = 0.
        wf(np+1) = 0.

        iarg = 13                                   ! AFTER call annih

        if (iausfl(iarg+1) .ne. 0) call ausgab(iarg)

        go1 = .true. ! must set for PHITS-egs because photons move to 1 and 2
        go to 8
      end if
      go1 = .true.
      return

 8    continue                                        ! ----------------
      return                                          ! Return to SHOWER
                                                      ! ----------------
 9    continue




      iarg = 6                                       ! BEFORE call brems

      if (iausfl(iarg+1) .ne. 0) call ausgab(iarg)


      call brems                ! To determine energies and polar angles


      uf(np) = 0.
      uf(1) = 0.
      vf(np) = 0.
      vf(1) = 0.
      wf(np) = 0.
      wf(1) = 0.












      if (ibrspl .eq. 1) then         ! Splitting has been requested



        if(fbrspl.eq.0.d0) then
          if(nbrspl.gt.0) then
            fbrspl = 1.d0/float(nbrspl)
          else
            write(506,105)
            ErrCha = ''
            ErrID = 'L:1175/R:egs5collis/F:egs5.f' !E86_002_001
            call ErrWrite(ErrID,ErrCha)
            write(*,105)
           stop
          endif
        endif




 10       continue
            write(506,106) MXSTACK,nbrspl,(2*nbrspl + 1)/3
            nbrspl = (2*nbrspl + 1)/3
            fbrspl = 1./float(nbrspl)

            if (nbrspl .eq. 1) then
              write(506,107) MXSTACK
              ibrspl=0
            end if

            if((np+nbrspl) .lt. MXSTACK) go to 11
          go to 10
 11       continue
        end if





        if (iq(np).eq.0) then
          npstrt = np - 1
          fdummy = u(np-1)
          u(np-1) = u(np)
          u(np) = fdummy
          fdummy = v(np-1)
          v(np-1) = v(np)
          v(np) = fdummy
          fdummy = w(np-1)
          w(np-1) = w(np)
          w(np) = fdummy
          fdummy = e(np-1)
          e(np-1) = e(np)
          e(np) = fdummy
          fdummy = wt(np-1)
          wt(np-1) = wt(np)
          wt(np) = fdummy
          idummy = iq(np-1)
          iq(np-1) = iq(np)
          iq(np) = idummy
          idummy = latch(np-1)
          latch(np-1) = latch(np)
          latch(np) = idummy
          fdummy = uf(np-1)
          uf(np-1) = uf(np)
          uf(np) = fdummy
          fdummy = vf(np-1)
          vf(np-1) = vf(np)
          vf(np) = fdummy
          fdummy = wf(np-1)
          wf(np-1) = wf(np)
          wf(np) = fdummy
          k1step(np) = k1step(np-1)
          k1step(np-1) = 0.
          k1init(np) = k1init(np-1)
          k1init(np-1) = 0.
          k1rsd(np) = k1rsd(np-1)
          k1rsd(np-1) = 0.
        else
          npstrt = np
        end if

                                     !  (because interaction reduced it)

        wt(1) = wt(1)*fbrspl         ! Adjust weight of initial photon
        frstbr = e(1)                ! Store energy of initial photon

        e(np) = e(np) + e(1)         ! Restore electron's initial energy
                                     !  (because interaction reduced it)

        iarg = 7                                      ! AFTER call brems

        if (iausfl(iarg+1) .ne. 0) call ausgab(iarg)

        icsplt = 1                        ! Initialize splitting counter

 12     continue
        if (icsplt .ge. nbrspl) go to 13
          icsplt = icsplt+1

          iarg = 6                                   ! BEFORE call brems

          if (iausfl(iarg+1) .ne. 0) call ausgab(iarg)


          call brems            ! To determine energies and polar angles






          if (iq(np) .eq. 0) then
            fdummy = u(np-1)
            u(np-1) = u(np)
            u(np) = fdummy
            fdummy = v(np-1)
            v(np-1) = v(np)
            v(np) = fdummy
            fdummy = w(np-1)
            w(np-1) = w(np)
            w(np) = fdummy
            fdummy = e(np-1)
            e(np-1) = e(np)
            e(np) = fdummy
            fdummy = wt(np-1)
            wt(np-1) = wt(np)
            wt(np) = fdummy
            idummy = iq(np-1)
            iq(np-1) = iq(np)
            iq(np) = idummy
            idummy = latch(np-1)
            latch(np-1) = latch(np)
            latch(np) = idummy
            fdummy = uf(np-1)
            uf(np-1) = uf(np)
            uf(np) = fdummy
            fdummy = vf(np-1)
            vf(np-1) = vf(np)
            vf(np) = fdummy
            fdummy = wf(np-1)
            wf(np-1) = wf(np)
            wf(np) = fdummy
            k1step(np) = k1step(np-1)
            k1step(np-1) = 0.
            k1init(np) = k1init(np-1)
            k1init(np-1) = 0.
            k1rsd(np) = k1rsd(np-1)
            k1rsd(np-1) = 0.
          end if
          wt(1) = wt(1)*fbrspl           ! Adjust weight of photon
          e(np) = e(np) + e(1)    ! Restore electron's initial energy

          iarg = 7                                    ! AFTER call brems

          if (iausfl(iarg+1) .ne. 0) call ausgab(iarg)

        go to 12

 13     continue


        e(np) = e(np) - frstbr


      if(latch(np).ne.0) then
        do i=1,np
          write(506,*) 'i,iq,e,wt=',i,iq(i),e(i),wt(i)
          write(*,*) 'i,iq,e,wt=',i,iq(i),e(i),wt(i)
        end do
      end if

      iarg = 7                                        ! AFTER call brems

      if (iausfl(iarg+1) .ne. 0) call ausgab(iarg)

      if (iq(np) .eq. 0) then                         ! ----------------
        return                                        ! Return to SHOWER
      else                                            ! ----------------
        go1 = .true.
        return
      end if
      end









      subroutine compt
      use EGS5_BCOMP_MOD

      USE egs5_thresh_mod !<-

      implicit none

      include 'include/egs5_h.f'               ! Main EGS5 "header" file


      include 'include/egs5_epcont.f'   ! COMMONs required by EGS5 code
      include 'include/egs5_misc.f'
      include 'include/egs5_stack.f'
      include 'include/egs5_uphiot.f'
      include 'include/egs5_useful.f'
      include 'include/counters.f'       ! Additional (non-EGS5) COMMONs

      include 'err.inc'

      real*8 rnnow,rnnow1,rnnow2,rnnow3                      ! Arguments
      integer iarg

      integer ishell,iqtmp,lvallc,irloc,lxlv

      common /tstara/ atmrc(10,8)      ! Ogawa 2023/7/18, (10,7) -> (10,8) for plasmon

      real*8 atmrc

      icompt = icompt + 1                  ! Count entry into subroutine

      irloc = ir(np)
      if (incohr(irloc) .ne. 1 .and. iprofr(irloc) .eq. 1) then
        write(506,101)
        ErrCha = ''
        ErrID = 'L:1416/R:compt/F:egs5.f' !E86_003_001
        call ErrWrite(ErrID,ErrCha)
        write(*,101)
        stop
      end if

      eig = e(np)
      egp = eig/RM
      br0i =1. + 2.*egp
      alph1 = log(br0i)
      alph2 = egp*(br0i + 1.)/(br0i*br0i)
      sumalp = alph1 + alph2
                                    ! Start main sampling-rejection loop
1     continue
        call randomset(rnnow)
                                              ! Start of 1/br
                                              ! subdistribution sampling
        if (alph1 .ge. sumalp*rnnow) then
          call randomset(rnnow)
          br = exp(alph1*rnnow)/br0i
                                              ! Start of br
                                              ! subdistribution sampling
        else
          call randomset(rnnow1)
          call randomset(rnnow)
          if (egp .ge. (egp + 1.)*rnnow) then
            call randomset(rnnow)
            rnnow1 = max(rnnow1,rnnow)
          end if
          br = ((br0i - 1.)*rnnow1 + 1.)/br0i
        end if

        esg = br*eig                           ! Set up secondary photon

        a1mibr = 1. - br
        esedef = eig*a1mibr + RM
        temp = RM*a1mibr/esg
        sinthe = max(0.D0,temp*(2. - temp))         ! Prevent sinthe < 0
        call randomset(rnnow)
        rejf3 = 1. - br*sinthe/(1. + br*br)

        irloc = ir(np)
        if (incohr(irloc) .eq. 1) then
          medium = med(irloc)
          alamb = 0.012398520
          xval = sqrt(temp/2.)*eig/alamb
          if (xval .ge. 5.E-3 .and. xval .le. 80.) then
            xlv = log(xval)
            lxlv = sco1(medium)*xlv + sco0(medium)
            sxz = sxz1(lxlv,medium)*xlv + sxz0(lxlv,medium)
          else if (xval .gt. 80.) then
            sxz = 1.
          else
            sxz = 0.
          end if
          rejf3 = rejf3*sxz
        end if
        if (rnnow .gt. rejf3) go to 1

      sinthe = sqrt(sinthe)
      costhe = 1. - temp

      irloc = ir(np)
      if (iprofr(irloc) .eq. 1) then
        esgmax = eig - cpimev
        valloc = sqrt(eig*eig + esgmax*esgmax - 2.*eig*esgmax*costhe)
        if (qvalmx .ge. 100.) then
          esg = eig/(1. + eig/RM*(1. - costhe))
          go to 2
        end if

 3      continue
        medium = med(irloc)
        call randomset(rnnow)

        if (icprof(medium) .eq. 1) then
          lvallc = cco1(medium)*rnnow + cco0(medium)
          cpr = cpr1(lvallc,medium)*rnnow + cpr0(lvallc,medium)
        end if

        if (icprof(medium) .eq. 3) then
          call randomset(rnnow1)
          do ishell=1,mxshel(medium)
            if (rnnow1 .le. elecno(ishell,medium)) go to 4
          end do
 4        continue
          lvallc = ccos1(medium)*rnnow + ccos0(medium)
        end if

        f1 = (cpr/137.)*(cpr/137.)
        f2 = (1. - costhe)/RM
        f3 = (1. + f2*eig)*(1. + f2*eig) - f1
        f4 = (f1*costhe - 1. - f2*eig)
        f5 = f4*f4 - f3 + f1*f3
        eps1 = 0.0
        if (f5 .lt. eps1) go to 3
        esg1 = (-f4 - sqrt(f5))/f3*eig
        esg2 = (-f4 + sqrt(f5))/f3*eig
        call randomset(rnnow2)
        if (rnnow2 .lt. 0.5) then
          esg = esg1
        else
          esg = esg2
        end if
        if (icprof(medium) .eq. 3) esgmax = eig - capio(ishell,medium)
        if (esg .gt. esgmax .or. esg .lt. 0.) go to 3
        call randomset(rnnow3)
        if (esg .lt. esgmax*rnnow3) go to 3
 2      continue
      end if

      ese = eig - esg + RM

      if (iprofr(irloc) .eq. 1) then
        ese = ese - capio(ishell,medium)
        ecapbind = ecapbind + capio(ishell,medium)



      end if

      if (lpolar(irloc).eq. 0) then

        call uphi(2,1)

      else

        call aphi(br)


        call uphi(3,1)

      end if

      np = np + 1                            ! Set up secondary electron
      psq = esedef**2 - RMSQ

      if (psq .le. 0.0) then                 ! To avoid division by zero
        costhe = 0.0
        sinthe = -1.0
      else
        costhe = (ese + esg)*a1mibr/sqrt(psq)
        sinthe = -sqrt(max(0.D0,1.D0 - costhe*costhe))
      end if

      call uphi(3,2)                             ! Set direction cosines

      uf(np) = 0.
      vf(np) = 0.
      wf(np) = 0.
      k1step(np) = 0.
      k1init(np) = 0.
      k1rsd(np) = 0.
      k1step(np-1) = 0.
      k1init(np-1) = 0.
      k1rsd(np-1) = 0.


      iq(np) = -1
      e(np) = ese
      e(np-1) = esg


      call npconv
      atmrc(1,4) = atmrc(1,4) + 1.d0
                                                      ! ----------------
      return                                          ! Return to PHOTON
                                                      ! ----------------
      end








      subroutine edgbin

      USE egs5_brempr_mod
      USE egs5_media_mod
      USE egs5_edge_mod
      USE egs5_photin_mod
      USE egs5_thresh_mod

      implicit none

      include 'include/egs5_h.f'               ! Main EGS5 "header" file


      include 'include/egs5_epcont.f'
      include 'include/egs5_misc.f'
      include 'include/egs5_stack.f'
      include 'include/egs5_uphiot.f'
      include 'include/egs5_useful.f'

      include 'include/counters.f'       ! Additional (non-EGS5) COMMONs

      include 'err.inc'

      real*8 eig,eee                                   ! Local variables
      integer ii,jj,izn,ner,iz1,ikl

      iedgbin = iedgbin + 1                ! Count entry into subroutine

      do medium=1,nmed
        ner = nepm(medium)

        if (ner .gt. MXEPERMED) then
          write(506,101)

         ErrCha = ''
         ErrID = 'L:1634/R:edgbin/F:egs5.f' !E86_004_001
         call ErrWrite(ErrID,ErrCha)
          write(*,101)
 101      FORMAT(' Number of elements in medium must be <= MXEPERMED !')
          stop
        end if

        nedgb(medium) = 0
        do izn=1,ner
          iz1 = zelem(medium,izn)
          do ikl=1,4
            eee = eedge(ikl,iz1)/1000.0
            if (eee .gt. ap(medium)) then
              nedgb(medium) = nedgb(medium) + 1
              eig = log(eee)
              edgb(nedgb(medium),medium) = eee
            end if
          end do
        end do

        if (nedgb(medium) .gt. 0) then
          do ii=1,nedgb(medium)
            do jj=1,nedgb(medium)
              if (ii .ne. jj) then
                if (ledgb(ii,medium) .eq. ledgb(jj,medium)) then
                  write(506,102)medium
                end if
              end if
            end do
          end do
        end if
      end do

      return                                           ! Return to HATCH

      end










      subroutine eii

      USE egs5_eiicom_mod
      USE egs5_edge_mod

      implicit none

      include 'include/egs5_h.f'               ! Main EGS5 "header" file


      include 'include/egs5_epcont.f'
      include 'include/egs5_misc.f'
      include 'include/egs5_stack.f'
      include 'include/egs5_uphiot.f'
      include 'include/egs5_useful.f'

      include 'include/counters.f'       ! Additional (non-EGS5) COMMONs

      include 'err.inc'

      real*8 rnnow                                           ! Arguments
      integer iarg

      real*8 capbind,etmp,ese1,ese2,ratio,esesum       ! Local variables
      real*8 x_sp,y_sp,z_sp,wt_sp,time_sp,dnear_sp
      integer iqtmp,jnp,ieie
      integer ir_sp,latch_sp,i
      common /tstara/ atmrc(10,8)      ! Ogawa 2023/7/18, (10,7) -> (10,8) for plasmon

      real*8 atmrc

      ieii = ieii + 1                      ! Count entry into subroutine

      nxray = 0
      nauger = 0
      capbind = eedge(1,iz)*0.001
      ese1 = e(np)  ! higher electron at np
      ese2 = e(1)   ! lower electron at 1
      esesum = e(np) + e(1) - 2.*RM
      if (ese1 - RM .gt. capbind .or. ese2 - RM .gt. capbind) then
        call randomset(rnnow)
        if (rnnow .gt. 0.5) then
          e(np) = ese1 - capbind
          e(1) = ese2
        else
          e(np) = ese1
          e(1) = ese2 - capbind
        end if
        if (e(np) .le. RM .or. e(1) .le. RM) then
          if (rnnow .le. 0.5) then
            e(np) = ese1 - capbind
            e(1) = ese2
          else
            e(np) = ese1
            e(1) = ese2 - capbind
          end if
        end if
      else
        ratio = 1. - capbind/esesum
        e(np) = (ese1 - rm)*ratio + RM
        e(1) = (ese2 - rm)*ratio + RM
      end if

      if (ieispl .eq. 1) then



        if(feispl.eq.0.d0) then
          if(neispl.gt.0) then
            feispl = 1.d0/float(neispl)
          else
            write(506,105)

           ErrCha = ''
           ErrID = 'L:1759/R:eii/F:egs5.f' !E86_005_001
           call ErrWrite(ErrID,ErrCha)

            write(*,105)
           stop
          endif
        endif

        if (neispl .gt. 1 .and. (np + neispl) .ge. MXSTACK) then
 1        continue
            write(506,100) MXSTACK,neispl,(2*neispl+1)/3
            neispl = (2*neispl + 1)/3
            feispl = 1./float(neispl)
            if (neispl .eq. 1) then
              write(506,200) MXSTACK
              ieispl = 0
            end if
            if (np + neispl .lt. MXSTACK) go to 2
          go to 1
 2        continue
        end if
      else
        neispl = 1
        feispl = 1.
      end if


      call kshell


      if (nxray .ge. 1 .and. exray(1) .gt. eedge(2,iz)*0.001) then
        enew = exray(1)
        ieie = 1
      else
        enew = 0.0
        ieie = 0
      end if
      edep = capbind - enew
      etmp = e(np)
      e(np)= edep
      iqtmp = iq(np)
      iq(np) = -1


      ecapbind=ecapbind+edep  ! dead electron energy

      iarg=4


      if (iausfl(iarg+1) .ne. 0) call ausgab(iarg)

      e(np) = etmp
      iq(np) = iqtmp
      if (ieie .eq. 1) then
        wt_sp = wt(np)*feispl
        x_sp = x(np)
        y_sp = y(np)
        z_sp = z(np)
        ir_sp = ir(np)
        time_sp = time(np)
        dnear_sp = dnear(np)
        latch_sp = latch(np)
        do jnp=1,neispl
          np = np + 1
          iq(np) = 0
          e(np) = enew
          call randomset(rnnow)
          costhe = 2.*rnnow - 1.
          sinthe = sqrt(1. - costhe*costhe)
          u(np) = 0.
          v(np) = 0.
          w(np) = 1.

          call uphi(2,1)

          x(np) = x_sp
          y(np) = y_sp
          z(np) = z_sp
          ir(np) = ir_sp
          wt(np) = wt_sp
          time(np) = time_sp
          dnear(np) = dnear_sp
          latch(np) = latch_sp
          k1step(np) = 0.
          k1init(np) = 0.
          k1rsd(np) = 0.

          call npconv
        end do
      end if
      if( nxray .ge. 1 ) atmrc(2,1) = atmrc(2,1) + 1.d0
      if( nauger .ge. 1 ) atmrc(2,2) = atmrc(2,2) + 1.d0
                                                      ! ----------------
      return                                          ! Return to MOLLER
                                                      ! ----------------
      end









      subroutine hardx(charge,kEnergy,keIndex,keFraction,sig0)

      USE egs5_elecin_mod
      USE egs5_edge_mod
      USE egs5_thresh_mod

      implicit none

      integer charge
      integer keIndex
      double precision kEnergy
      double precision keFraction

      include 'include/egs5_h.f'               ! Main EGS5 "header" file


      include 'include/egs5_useful.f'

      include 'include/counters.f'       ! Additional (non-EGS5) COMMONs

      double precision sig0

      integer mollerIndex
      double precision mollerThresh
      double precision logMollerThresh

      ihardx = ihardx + 1                  ! Count entry into subroutine

      if (charge .lt. 0) then                                  ! e-
        iextp = 0
        !  correction for Moller threshold ?
        if (kEnergy .lt. (thmoll(medium)-RM)*1.5) then
          mollerThresh = thmoll(medium) - RM
          logMollerThresh = log(mollerThresh)
          mollerIndex = eke1(medium)*logMollerThresh + eke0(medium)
          if (mollerIndex .eq. keIndex) then
            if (thmoll(medium)-RM .le. kEnergy) then
              iextp = 1
            else
              iextp = -1
            end if
          end if
        end if
        sig0 = esig1(keIndex+iextp,medium)*keFraction +&
                                   esig0(keIndex+iextp,medium)

       else                                                    ! e+
         sig0 = psig1(keIndex,medium)*keFraction +&
                                         psig0(keIndex,medium)
       end if
      if(sig0.le.0.0)sig0=1.e-10

       return
       end


      subroutine hatch
      use EGS5_BCOMP_MOD

      USE egs5_brempr_mod
      USE egs5_eiicom_mod
      USE egs5_elecin_mod
      USE egs5_media_mod
      USE egs5_edge_mod
      USE egs5_photin_mod
      USE egs5_scpw_mod
      USE egs5_thresh_mod

      implicit none
      save
      include 'include/egs5_h.f'               ! Main EGS5 "header" file

      include 'include/egs5_bounds.f'      ! COMMONs required by EGS5 code

      include 'include/egs5_misc.f'

      include 'include/egs5_stack.f'

      include 'include/egs5_uphiin.f' ! Probably don't need this anymore
      include 'include/egs5_uphiot.f'
      include 'include/egs5_useful.f'
      include 'include/egs5_userpr.f'
      include 'include/egs5_usersc.f'
      include 'include/egs5_uservr.f'
      include 'include/egs5_userxt.f'

      include 'include/counters.f'       ! Additional (non-EGS5) COMMONs
      include 'include/randomm.f'

      include 'err.inc'

      real*8 rnnow                                           ! Arguments
      real*8, allocatable :: sxz0_W(:,:), sxz1_W(:,:)
      real*8, allocatable :: cpr0_W(:,:), cpr1_W(:,:)

      character mbuf(72),mdlabl(8)


        imxmed = nmed  ! = mxmat <= phits2 input . material .
      
      IF( flg_local .ne. 0 ) then
         ErrCha = ''
         ErrID = 'L:2130/R:hatch/F:egs5.f' !E86_006_001
         call ErrWrite(ErrID,ErrCha)

       write(*,*)'# Hatch local array can not allocate !! '
       STOP '# Hatch local array can not allocate !! '
      ENDIF


      ihatch = ihatch + 1                  ! Count entry into subroutine

      if(iblock.ne.1) then            ! Check for block_set initizations
        write(506,5000)
         ErrCha = ''
         ErrID = 'L:2143/R:hatch/F:egs5.f' !E86_007_001
         call ErrWrite(ErrID,ErrCha)
        write(*,5000)
        stop
      endif

      if (i1st .ne. 0) then
        i1st=0
        if (.not.rluxset) then
          write(506,*) 'RNG ranlux not initialized:  doing so in HATCH'
          call rluxgo
        end if
        latchi = 0.0

        nisub = mxsinc_loc - 2
        fnsss = nsinss
        wid = PI5D2/float(nisub)
        wss = wid/(fnsss - 1.0)
        zeros(1) = 0.
        zeros(2) = PI
        zeros(3) = TWOPI
        do isub=1,mxsinc_loc
          sx = 0.
          sy = 0.
          sxx = 0.
          sxy = 0.
          xs0 = wid*float(isub - 2)
          xs1 = xs0 + wid
          iz = 0
          do izz=1,3
            if (xs0 .le. zeros(izz) .and. zeros(izz) .le. xs1) then
              iz = izz
              go to 1
            end if
          end do
 1        continue
          if (iz .eq. 0) then
            xsi = xs0
          else
            xsi = zeros(iz)
          end if
          do iss=1,nsinss
            xs = wid*float(isub - 2) + wss*float(iss - 1) -xsi
            ys = sin(xs + xsi)
            sx = sx + xs
            sy = sy + ys
            sxx = sxx + xs*xs
            sxy = sxy + xs*ys
          end do
          if (iz .ne. 0) then
            sin1(isub) = sxy/sxx
            sin0(isub) = -sin1(isub)*xsi
          else
            del = fnsss*sxx-sx*sx
            sin1(isub) = (fnsss*sxy - sy*sx)/del
            sin0(isub) = (sy*sxx - sx*sxy)/del - sin1(isub)*xsi
          end if
        end do

        sinc0 = 2.
        sinc1 = 1./wid
        if (istest .ne. 0) then
          adev = 0.
          rdev = 0.
          s2c2mn = 10.
          s2c2mx = 0.
          do isub=1,nisub
            do iss=1,nsinss
              theta = wid*float(isub - 1) + wss*float(iss - 1)
              cthet = PI5D2 - theta
              sinthe = sin(theta)
              costhe = sin(cthet)
              sint = sin(theta)
              cost = cos(theta)
              asd = abs(sinthe - sint)
              acd = abs(costhe - cost)
              adev = max(adev,asd,acd)
              if (sint .ne. 0.) rdev = max(rdev,asd/abs(sint))
              if (cost .ne. 0.) rdev = max(rdev,acd/abs(cost))
              s2c2 = sinthe**2 + costhe**2
              s2c2mn = min(s2c2mn,s2c2)
              s2c2mx = max(s2c2mx,s2c2)

            end do
          end do
          write(506,1350) mxsinc_loc,nsinss
          write(506,1360) adev,rdev,s2c2mn,s2c2mx
          adev = 0.
          rdev = 0.
          s2c2mn = 10.
          s2c2mx = 0.
          do irn=1,nrna
            call randomset(rnnow)
            theta = rnnow*PI5D2
            cthet = PI5D2 - theta
            sinthe = sin(theta)
            costhe = sin(cthet)
            sint = sin(theta)
            cost = cos(theta)
            asd = abs(sinthe - sint)
            acd = abs(costhe - cost)
            adev = max(adev,asd,acd)
            if (sint .ne. 0.) rdev = max(rdev,asd/abs(sint))
            if (cost .ne. 0.) rdev = max(rdev,acd/abs(cost))
            s2c2 = sinthe**2 + costhe**2
            s2c2mn = min(s2c2mn,s2c2)
            s2c2mx = max(s2c2mx,s2c2)
          end do
          write(506,1380) nrna
          write(506,1390) adev,rdev,s2c2mn,s2c2mx
        end if
        p = 1.
        do i=1,50
          pwr2i(i) = p
          p = p/2.
        end do
      end if

      do j=1,nmed
        do i=1,nreg
          if (iraylr(i) .eq. 1 .and. med(i) .eq. j) then
            iraylm(j) = 1
            go to 2
          end if
        end do
 2      continue
      end do

      do j=1,nmed
        do i=1,nreg
          if (incohr(i) .eq. 1 .and. med(i).eq.j) then
            incohm(j) = 1
            go to 20
          end if
        end do
 20     continue
      end do

      do j=1,nmed
        do i=1,nreg
          if (iprofr(i) .eq. 1 .and. med(i) .eq. j) then
            iprofm(j) = 1
            go to 21
          end if
        end do
 21     continue
      end do

      do j=1,nmed
        do i=1,nreg
          if (impacr(i) .eq. 1 .and. med(i) .eq. j) then
            impacm(j) = 1
            go to 22
          end if
        end do
 22     continue
      end do

      rewind kmpi

      nm = 0
      do im=1,nmed
        lok(im) = 0
        if (iraylm(im) .eq. 1) write(506,1440) im
      end do

      do im=1,nmed
        if (incohm(im) .eq. 1) write(506,4090) im
      end do

      do im=1,nmed
        if (iprofm(im) .eq. 1) write(506,4110) im
      end do

      do im=1,nmed
        if (impacm(im) .eq. 1) write(506,4130) im
      end do

 5    continue

      read(kmpi,1270,end=1470) mbuf
      do ib=1,lmdl
        if (mbuf(ib) .ne. mdlabl(ib)) go to 5
      end do
      do 6 im=1,nmed
        do ib=1,lmdn
          il = lmdl + ib
          if (mbuf(il) .ne. media(ib,im)) go to 6
          if (ib .eq. lmdn) go to 7
        end do
 6    continue
      go to 5

 7    continue
      if (lok(im) .ne. 0) go to 5
      lok(im) = 1
      nm = nm + 1

      write(kmpo,1510) im,mbuf

      go to 8
 1000 backspace(kmpi)


 8    continue
      write(kmpo,1540) (mbuf(i),i=1,5),rhom(im),nne(im)


      write(kmpo,1580)
      read(kmpi,1260) rlcm(im),ae(im),ap(im),ue(im),up(im)
      write(kmpo,1260)rlcm(im),ae(im),ap(im),ue(im),up(im)

      te(im) = ae(im) - RM
      thmoll(im) = te(im)*2. + RM

      nsge = msge(im)
      nge = mge(im)
      nseke = mseke(im)
      neke = meke(im)
      nleke = mleke(im)
      ncmfp = mcmfp(im)
      nrange = mrange(im)


      write(kmpo,4340)
      read(kmpi,1260) tebinda,ge0(im),ge1(im)
      write(kmpo,1260) tebinda,ge0(im),ge1(im)

      write(kmpo,1660)
     
      if (iraylm(im) .eq. 1 .and. (irayl .lt. 1 .or. irayl .gt. 2)) then
        write(506,4360) im

        ErrCha = ''
        ErrID = 'L:2427/R:hatch/F:egs5.f' !E86_008_001
        call ErrWrite(ErrID,ErrCha)
        write(*,4360) im
        stop
      end if

      if (incohm(im) .eq. 1 .and. incoh .ne. 1) then
        write(506,4370) im
        ErrCha = ''
        ErrID = 'L:2436/R:hatch/F:egs5.f' !E86_009_001
        call ErrWrite(ErrID,ErrCha)
        write(*,4370) im
        stop
      end if

      if (incohm(im) .eq. 1 .and. ibound .ne. 1) then
        write(506,4375) im
        ErrCha = ''
        ErrID = 'L:2445/R:hatch/F:egs5.f' !E86_010_001
        call ErrWrite(ErrID,ErrCha)
        write(*,4375) im
        stop
      end if

      if (iprofm(im) .eq. 1) then
        if(icprof(im) .ne. 3) then
          write(506,4380) im
          ErrCha = ''
          ErrID = 'L:2455/R:hatch/F:egs5.f' !E86_011_001
          call ErrWrite(ErrID,ErrCha)
          write(*,4380) im
          stop
        else if(incohm(im) .eq. 0) then
          write(506,4390) im
          ErrCha = ''
          ErrID = 'L:2462/R:hatch/F:egs5.f' !E86_012_001
          call ErrWrite(ErrID,ErrCha)
          write(*,4390) im
          stop
        else if(ibound .eq. 0) then
          write(506,4400) im
          ErrCha = ''
          ErrID = 'L:2469/R:hatch/F:egs5.f' !E86_013_001
          call ErrWrite(ErrID,ErrCha)
          write(*,4400) im
          stop
 50   endif
      end if

      if (impacm(im) .eq. 1 .and. impact .eq. 0) then
        write(506,4410) im
        ErrCha = ''
        ErrID = 'L:2479/R:hatch/F:egs5.f' !E86_014_001
        call ErrWrite(ErrID,ErrCha)
        write(*,4410) im
        stop
      end if

      if (irayl .eq. 1 .or. irayl .eq. 2 .or. irayl .eq. 3) then
        write(kmpo,1680)
        read(kmpi,1250) ngr(im)
        write(kmpo,1250) ngr(im)

        ngrim = ngr(im)

        write(kmpo,1690)
        read(kmpi,1260) rco0(im),rco1(im)
        write(kmpo,1260) rco0(im),rco1(im)
        write(kmpo,1700)
        read(kmpi,1260) (rsct0(i,im),rsct1(i,im),i=1,ngrim)
        write(kmpo,1260) (rsct0(i,im),rsct1(i,im),i=1,ngrim)
        write(kmpo,1710)
        read(kmpi,1260) (cohe0(i,im),cohe1(i,im),i=1,nge)
        write(kmpo,1260) (cohe0(i,im),cohe1(i,im),i=1,nge)

        if (iraylm(im) .ne. 1) write(506,1720) im
      end if

      if (incoh .eq. 1) then
        write(kmpo,4450)
        read(kmpi,1250) ngs(im)
        write(kmpo,1250) ngs(im)
        ngsim = ngs(im)

        if( nm .eq. 1)then
         
        elseif( nm .gt. 1 )then
          if( ngsim .gt. mx_ngsim )then
           
            do j=1, im-1 !  save until previous material = im-1
            do i=1, mx_ngsim ! max ngs(j)
               sxz0_W(i, j) = sxz0(i,j)
               sxz1_W(i, j) = sxz1(i,j)
            enddo
            enddo

             deallocate(sxz0 ,sxz1  )
             
             do j=1,im - 1      !<- current material im-1
             do i=1,mx_ngsim   !<- MXNS
              sxz0(i,j) = sxz0_W(i,j)
              sxz1(i,j) = sxz1_W(i,j)
             enddo
             enddo

             deallocate(sxz0_W ,sxz1_W)
          endif
        endif

        write(kmpo,4460)
        read(kmpi,1260) sco0(im),sco1(im)
        write(kmpo,1260) sco0(im),sco1(im)
        write(kmpo,4470)
        read(kmpi,1260) (sxz0(i,im),sxz1(i,im),i=1,ngsim)
        write(kmpo,1260) (sxz0(i,im),sxz1(i,im),i=1,ngsim)
        if (incohm(im) .ne. 1) write(506,4480) im
      end if

        if (icprof(im) .eq. 1 .or. icprof(im) .eq. 2) then
          write(kmpo,4490)
          read(kmpi,1250) ngc(im)
          write(kmpo,1250) ngc(im)
          ngcim = ngc(im)

        if( nm .eq. 1 ) then
             
        elseif( nm .gt. 1)then
          if( ngcim .gt. mx_ngcim)then
             
            do j=1,im-1 !  save until previous material = im-1
            do i=1, mx_ngcim ! max ngc(j)
               cpr0_W(i,j) = cpr0(i,j)
               cpr1_W(i,j) = cpr1(i,j)
            enddo
            enddo

           deallocate( cpr0 ,cpr1)
          
           do j=1, im-1 !  restore until previous material = im-1
           do i=1,mx_ngcim   !<- MXCP
              cpr0(i,j) = cpr0_W(i,j)
              cpr1(i,j) = cpr1_W(i,j)
           enddo
           enddo

           deallocate( cpr0_W, cpr1_W)
          endif
        endif

          write(kmpo,4500)
          read(kmpi,1260) cco0(im),cco1(im),cpimev
          write(kmpo,1260) cco0(im),cco1(im),cpimev
          write(kmpo,4510)
          read(kmpi,1260) (cpr0(i,im),cpr1(i,im),i=1,ngcim)
          write(kmpo,1260) (cpr0(i,im),cpr1(i,im),i=1,ngcim)
           if (iprofm(im) .ne. 1) write(506,4520) im
        end if

        if (icprof(im) .eq. 3 .or. icprof(im) .eq. 4) then
          write(kmpo,4530)
          read(kmpi,1250) mxshel(im),ngc(im)
          write(kmpo,1250) mxshel(im),ngc(im)
          ngcim = ngc(im)
          mxsim = mxshel(im)

          if    ( nm .eq. 1)then
             
          elseif( nm .gt. 1)then
           if( (ngcim .gt. mx_ngcim3) .or. (mxsim .gt. mx_mxsim) )then
            
             do k=1, im-1 !  save until previous material = im-1
             do j=1,mx_mxsim    ! max mxshel(k) MXNS
                   elecno_W(j,k) = elecno(j,k)
                    capio_W(j,k) =  capio(j,k)
               do i=1,mx_ngcim3  ! max ngc(k)  MXCP
                    cprs0_W(i,j,k) = cprs0(i,j,k)
                    cprs1_W(i,j,k) = cprs1(i,j,k)
               enddo
             enddo
             enddo

             deallocate(cprs0,cprs1, elecno,capio)
            
             do k=1, im-1 !  restore until previous material = im-1
             do j=1,mx_mxsim    ! max mxshel(k) MXNS
                   elecno(j,k) = elecno_W(j,k)
                    capio(j,k) =  capio_W(j,k)
               do i=1,mx_ngcim3 ! max ngc(k)  MXCP
                  cprs0(i,j,k) = cprs0_W(i,j,k)
                  cprs1(i,j,k) = cprs1_W(i,j,k)
               enddo
             enddo
             enddo

             deallocate( cprs0_W, cprs1_W, elecno_W, capio_W)
           endif !<- ngcim > mx_ngcim3
          endif

          write(kmpo,4540)
          read(kmpi,1260) (elecno(i,im),i=1,mxsim)
          write(kmpo,1260) (elecno(i,im),i=1,mxsim)
          write(kmpo,4550)
          read(kmpi,1260) (capio(i,im),i=1,mxsim)
          write(kmpo,1260) (capio(i,im),i=1,mxsim)
          write(kmpo,4560)
          read(kmpi,1260) ccos0(im),ccos1(im)
          write(kmpo,1260) ccos0(im),ccos1(im)
          write(kmpo,4570)

        if (impact .ge. 1) then
          write(kmpo,4590)
          read(kmpi,1250) nepm(im)
          write(kmpo,1250) nepm(im)
          ner = nepm(im)
          write(kmpo,4610)
          read(kmpi,1250) neii(im)
          write(kmpo,1250) neii(im)
          write(kmpo,4620)
          read(kmpi,1260) eico0(im),eico1(im)
          write(kmpo,1260) eico0(im),eico1(im)
          write(kmpo,4630)
         

       if( incoh     .eq. 1 )then
          mx_ngsim = max( ngsim ,mx_ngsim) !<- max for MXSCTFF
       endif ! incoh = 1


       if( icprof(im) .eq. 1 .or. icprof(im) .eq.2)then
          mx_ngcim = max(ngcim,mx_ngcim) !<- max for MXCP
       endif ! icprof(im) = 1 or 2


       if( icprof(im) .eq. 3 .or. icprof(im) .eq.4)then
          mx_ngcim3 = max(ngcim,mx_ngcim3) !<- max for MXCP
          mx_mxsim = max(mxsim,mx_mxsim) !<- max for MXNS
       endif ! icprof(im) = 3 or 4
      if (nm .ge. nmed) go to 9
      go to 5

 9    continue








      if(emaxe .eq. 0.d0) then
        emaxe=1.d50
        do j=1,nmed
          emaxe=min(ue(j),up(j)+RM,emaxe)
        end do
        write(506,5006) emaxe
      else if(emaxe .lt. 1.d-4) then
        write(506,5007) emaxe
        ErrCha = ''
        ErrID = 'L:2719/R:hatch/F:egs5.f' !E86_015_001
        call ErrWrite(ErrID,ErrCha)
        write(*,5007) emaxe
        stop
      else
        do j=1,nmed
          if (emaxe .gt. ue(j)) then
            write(506,5008) emaxe, ue(j), j
            ErrCha = ''
            ErrID = 'L:2728/R:hatch/F:egs5.f' !E86_016_001
            call ErrWrite(ErrID,ErrCha)
            write(*,5008) emaxe, ue(j), j
            stop
          end if
          if (emaxe-RM .gt. up(j)) then
            write(506,5009) emaxe-RM, up(j), j
            ErrCha = ''
            ErrID = 'L:2736/R:hatch/F:egs5.f' !E86_017_001
            call ErrWrite(ErrID,ErrCha)
            write(*,5009) emaxe-RM, up(j), j
            stop
          end if
        end do
      end if

      dunitr = dunit
      if (dunit .lt. 0.) then

        id = max0(1,min0(nmed,int(-dunit)))
        dunit = rlcm(id)
      end if
      if (dunit .ne. 1.) write(506,1730) dunitr,dunit

      do im=1,nmed
        dfact = rlcm(im)/dunit
        dfacti = 1.0/dfact

        i = 1
        go to 11
 10     i = i + 1
 11     if (i - meke(im) .gt. 0) go to 12
        esig0(i,im) = esig0(i,im)*dfacti
        esig1(i,im) = esig1(i,im)*dfacti
        psig0(i,im) = psig0(i,im)*dfacti
        psig1(i,im) = psig1(i,im)*dfacti
        ededx0(i,im) = ededx0(i,im)*dfacti
        ededx1(i,im) = ededx1(i,im)*dfacti
        pdedx0(i,im) = pdedx0(i,im)*dfacti
        pdedx1(i,im) = pdedx1(i,im)*dfacti
        tmxs0(i,im) = tmxs0(i,im)*dfact
        tmxs1(i,im) = tmxs1(i,im)*dfact
        go to 10
 12     continue

        i = 1
        go to 14
 13     i = i + 1
 14     if (i - mleke(im) .gt. 0) go to 15
        if (dunitr.eq.1) then
          dfactr = 1.d0
        else
          dfactr = 1.d0 / dunit
        end if
        escpw0(i,im) = escpw0(i,im)*dfactr
        escpw1(i,im) = escpw1(i,im)*dfactr
        pscpw0(i,im) = pscpw0(i,im)*dfactr
        pscpw1(i,im) = pscpw1(i,im)*dfactr
        ekini0(i,im) = ekini0(i,im)*dfactr
        ekini1(i,im) = ekini1(i,im)*dfactr
        pkini0(i,im) = pkini0(i,im)*dfactr
        pkini1(i,im) = pkini1(i,im)*dfactr
        erang0(i,im) = erang0(i,im)*dfactr
        erang1(i,im) = erang1(i,im)*dfactr
        prang0(i,im) = prang0(i,im)*dfactr
        prang1(i,im) = prang1(i,im)*dfactr
        go to 13
 15     continue

        teff0(im) = teff0(im)*dfact
        blcc(im) = blcc(im)*dfacti
        xcc(im) = xcc(im)*sqrt(dfacti)
        rldu(im) = rlcm(im)/dunit

        i = 1
        go to 17
 16     i = i + 1
 17     if (i - mge(im) .gt. 0) go to 18
        gmfp0(i,im) = gmfp0(i,im)*dfact
        gmfp1(i,im) = gmfp1(i,im)*dfact
        go to 16
 18     continue
      end do

      ecutmn = 1.d20
      vacdst = vacdst*dunito/dunit
      dunito = dunit
      do jr=1,nreg
        md = med(jr)
        if (md .ge. 1 .and. md .le. nmed) then
          ecut(jr) = max(ecut(jr),ae(md))
          pcut(jr) = max(pcut(jr),ap(md))
          ecutmn = min(ecutmn,ecut(jr))
          !  get the range at cutoff for e- and e+, so we don't
          !  have to comput it continually
          eke = ecut(jr) - RM
          elke = log(eke)
          lelke = eke1(md)*elke + eke0(md)
          ectrng(jr) = erang1(lelke,md)*elke + erang0(lelke,md)
          pctrng(jr) = prang1(lelke,md)*elke + prang0(lelke,md)
          if (rhor(jr) .eq. 0.0) rhor(jr) = rhom(md)
        end if
      end do

      if (ibrdst .eq. 1 .or. iprdst .gt. 0) then
        do im=1,nmed
          zbrang(im) = 0.
          pznorm = 0.
          do ie=1,nne(im)
            zbrang(im) = zbrang(im) + pz(im,ie)*zelem(im,ie)*&
                        (zelem(im,ie) + 1.0)
            pznorm = pznorm + pz(im,ie)
          end do
          zbrang(im) = (8.116224E-05)*(zbrang(im)/pznorm)**(1./3.)
        end do
      end if

      do ii=1,nreg
        if (iedgfl(ii) .ne. 0 .or. iauger(ii) .ne. 0) then

          call edgbin

          go to 19
        end if
      end do

 19   continue

      if (nmed .eq. 1) then
        write(506,1830)
      else
        write(506,1840) nmed
      end if


      call rk1                         ! get problem scattering strength



      call rmsfit                        ! read multiple scattering data





      do i=1,nmed
        if(emaxe/ue(i).lt.0.1d0 .and. charD(i).eq.0.d0) then
          write(506,5005) (media(j,i),j=1,24)
        endif
      end do
      return



      subroutine kauger

      USE egs5_media_mod
      USE egs5_edge_mod

      implicit none

      include 'include/egs5_h.f'               ! Main EGS5 "header" file

      include 'include/egs5_bounds.f'    ! COMMONs required by EGS5 code


      include 'include/egs5_epcont.f'
      include 'include/egs5_misc.f'
      include 'include/egs5_stack.f'
      include 'include/egs5_uphiot.f'
      include 'include/egs5_useful.f'

      include 'include/counters.f'       ! Additional (non-EGS5) COMMONs

      real*8 rnnow                                           ! Arguments

      integer kaug                                     ! Local variables

      ikauger = ikauger + 1                ! Count entry into subroutine

      if (dfkaug(13,iz) .eq. 0.) return

      nauger = nauger + 1
      call randomset(rnnow)
      do kaug=1,13
        if (rnnow .le. dfkaug(kaug,iz)) then
          eauger(nauger) = ekaug(kaug,iz)*1.E-3
          go to 1
        end if
      end do
      eauger(nauger) = ekaug(14,iz)*1.E-3

 1    continue
      if (kaug .eq. 1) then

        call lshell(1)


        call lshell(1)

      else if (kaug .eq. 2) then

        call lshell(1)


        call lshell(2)

      else if (kaug .eq. 3) then

        call lshell(1)


        call lshell(3)

      else if (kaug .eq. 4) then

        call lshell(2)


        call lshell(2)

      else if (kaug.eq.5) then

        call lshell(2)


        call lshell(3)

      else if (kaug .eq. 6) then

        call lshell(3)


        call lshell(3)

      else if (kaug .eq. 7 .or. kaug .eq. 10) then

        call lshell(1)

      else if (kaug .eq. 8 .or. kaug .eq. 11) then

        call lshell(2)

      else if (kaug .eq. 9 .or. kaug .eq. 12) then

        call lshell(3)

      else
        return
      end if

      return
      end








      subroutine kshell

      USE egs5_media_mod
      USE egs5_edge_mod

      implicit none

      include 'include/egs5_h.f'               ! Main EGS5 "header" file


      include 'include/egs5_bounds.f'    ! COMMONs required by EGS5 code
      include 'include/egs5_epcont.f'
      include 'include/egs5_misc.f'
      include 'include/egs5_stack.f'
      include 'include/egs5_uphiot.f'
      include 'include/egs5_useful.f'

      include 'include/counters.f'       ! Additional (non-EGS5) COMMONs

      real*8 rnnow                                           ! Arguments

      ikshell = ikshell + 1                ! Count entry into subroutine

      call randomset(rnnow)

      if (rnnow .gt. omegak(iz)) then

        call kauger

      else

        call kxray

      end if

      return

      end








      subroutine kxray

      USE egs5_media_mod
      USE egs5_edge_mod

      implicit none

      include 'include/egs5_h.f'               ! Main EGS5 "header" file

      include 'include/egs5_bounds.f'    ! COMMONs required by EGS5 code


      include 'include/egs5_epcont.f'
      include 'include/egs5_misc.f'
      include 'include/egs5_stack.f'
      include 'include/egs5_uphiot.f'
      include 'include/egs5_useful.f'

      include 'include/counters.f'       ! Additional (non-EGS5) COMMONs

      real*8 rnnow                                           ! Arguments

      integer ik                                       ! Local variables

      ikxray = ikxray + 1                  ! Count entry into subroutine

      if (dfkx(9,iz) .eq. 0.) return

      nxray = nxray + 1
      call randomset(rnnow)
      do ik=1,9
        if (rnnow .le. dfkx(ik,iz)) then
          exray(nxray) = ekx(ik,iz)*1.E-3
          go to 1
        end if
      end do
      exray(nxray) = ekx(10,iz)*1.E-3

 1    continue

      if (ik .eq. 1) call lshell(3)


      if (ik .eq. 2) call lshell(2)


      if (ik .eq. 3) call lshell(1)

      return

      end








      subroutine lauger(ll)

      USE egs5_edge_mod
      USE egs5_media_mod

      implicit none

      include 'include/egs5_h.f'               ! Main EGS5 "header" file

      include 'include/egs5_bounds.f'    ! COMMONs required by EGS5 code


      include 'include/egs5_epcont.f'
      include 'include/egs5_misc.f'
      include 'include/egs5_stack.f'
      include 'include/egs5_uphiot.f'
      include 'include/egs5_useful.f'

      include 'include/counters.f'       ! Additional (non-EGS5) COMMONs

      real*8 rnnow                                           ! Arguments
      integer ll

      integer laug                                     ! Local variables

      ilauger = ilauger + 1                ! Count entry into subroutine

      call randomset(rnnow)
      go to (1,2,3) ll

 1    continue
      if (dfl1aug(5,iz) .eq. 0.) return

      nauger = nauger + 1
      do laug=1,5
        if (rnnow .le. dfl1aug(laug,iz)) then
          eauger(nauger) = el1aug(laug,iz)*1.E-3
          return
        end if
      end do
      eauger(nauger) = el1aug(6,iz)*1.E-3
      return

 2    continue
      if (dfl2aug(5,iz) .eq. 0.) return

      nauger = nauger + 1
      do laug=1,5
        if (rnnow .le. dfl2aug(laug,iz)) then
          eauger(nauger) = el2aug(laug,iz)*1.E-3
          return
        end if
      end do
      eauger(nauger) = el2aug(6,iz)*1.E-3
      return

 3    continue
      if (dfl3aug(5,iz) .eq. 0.) return

      nauger = nauger + 1
      do laug=1,5
        if (rnnow .le. dfl3aug(laug,iz)) then
          eauger(nauger) = el3aug(laug,iz)*1.E-3
          return
        end if
      end do
      eauger(nauger) = el3aug(6,iz)*1.E-3
      return

      end








      subroutine lshell(ll)

      USE egs5_media_mod
      USE egs5_edge_mod

      implicit none

      include 'include/egs5_h.f'               ! Main EGS5 "header" file

      include 'include/egs5_bounds.f'    ! COMMONs required by EGS5 code


      include 'include/egs5_epcont.f'
      include 'include/egs5_misc.f'
      include 'include/egs5_stack.f'
      include 'include/egs5_uphiot.f'
      include 'include/egs5_useful.f'

      include 'include/counters.f'       ! Additional (non-EGS5) COMMONs

      real*8 rnnow                                           ! Arguments
      integer ll

      integer ickflg                                   ! Local variables

      ilshell = ilshell + 1                ! Count entry into subroutine

      ickflg = 0
      if (ll .eq. 2) go to 1
      if (ll .eq. 3) go to 2
      call randomset(rnnow)
      ebind = eedge(2,iz)*1.E-3
      if (rnnow .gt. omegal1(iz)) then
        if (rnnow .le. omegal1(iz) + f12(iz)) then
          ickflg = 1
          go to 1
        else if (rnnow .le. omegal1(iz) + f12(iz) + f13(iz)) then
          ickflg = 1
          go to 2
        else

          call lauger(1)

        end if
      else

        call lxray(1)

      end if
      return

 1    continue
      call randomset(rnnow)
      if (ickflg .eq. 0) then
        ebind = eedge(3,iz)*1.E-3
      end if
      if (rnnow .gt. omegal2(iz)) then
        if (rnnow .le. omegal2(iz) + f23(iz)) then
          ickflg = 1
          go to 2
        else

          call lauger(2)

        end if
      else

        call lxray(2)

      end if
      return

 2    continue
      if (ickflg .eq. 0) ebind = eedge(4,iz)*1.E-3
      call randomset(rnnow)
      if (rnnow .gt. omegal3(iz)) then

        call lauger(3)

      else

        call lxray(3)

      end if

      return

      end








      subroutine lxray(ll)

      USE egs5_edge_mod
      USE egs5_media_mod

      implicit none

      include 'include/egs5_h.f'               ! Main EGS5 "header" file


      include 'include/egs5_bounds.f'    ! COMMONs required by EGS5 code
      include 'include/egs5_epcont.f'
      include 'include/egs5_misc.f'
      include 'include/egs5_stack.f'
      include 'include/egs5_uphiot.f'
      include 'include/egs5_useful.f'

      include 'include/counters.f'       ! Additional (non-EGS5) COMMONs

      real*8 rnnow                                           ! Arguments
      integer ll

      integer lx                                       ! Local variables

      ilxray = ilxray + 1                  ! Count entry into subroutine

      call randomset(rnnow)
      go to (1,2,3) ll

 1    continue
      if (dflx1(7,iz) .eq. 0.) return

      nxray = nxray + 1
      do lx=1,7
        if (rnnow .le. dflx1(lx,iz)) then
          exray(nxray) = elx1(lx,iz)*1.E-3
          return
        end if
      end do
      exray(nxray) = elx1(8,iz)*1.E-3
      return

 2    continue
      if (dflx2(4,iz) .eq. 0.) return

      nxray = nxray + 1
      do lx=1,4
        if (rnnow .le. dflx2(lx,iz)) then
          exray(nxray) = elx2(lx,iz)*1.E-3
          return
        end if
      end do
      exray(nxray) = elx2(5,iz)*1.E-3
      return

 3    continue
      if (dflx3(6,iz) .eq. 0.) return

      nxray = nxray + 1
      do lx=1,6
        if (rnnow .le. dflx3(lx,iz)) then
          exray(nxray) = elx3(lx,iz)*1.E-3
          return
        end if
      end do
      exray(nxray) = elx3(7,iz)*1.E-3
      return

      end










      subroutine moller

      USE egs5_brempr_mod
      USE egs5_eiicom_mod
      USE egs5_elecin_mod
      USE egs5_edge_mod
      USE egs5_thresh_mod

      implicit none
      include 'include/egs5_h.f'               ! Main EGS5 "header" file


      include 'include/egs5_epcont.f'
      include 'include/egs5_misc.f'
      include 'include/egs5_stack.f'

      include 'include/egs5_uphiot.f'
      include 'include/egs5_useful.f'

      include 'include/counters.f'       ! Additional (non-EGS5) COMMONs

      real*8 rnnow                                           ! Arguments
      integer iarg

      integer ifun,lelke2
      common /paraj/ mstz(300), parz(300)
      real*8  parz
      integer mstz,i
      common /tstara/ atmrc(10,8)      ! Ogawa 2023/7/18, (10,7) -> (10,8) for plasmon

      real*8 atmrc

      imoller = imoller + 1                ! Count entry into subroutine

      eie   = e(np)
      ekin  = eie - RM
      t0    = ekin/RM
      e0    = t0 + 1.0
      extrae= eie - thmoll(medium)
      e02   = e0*e0
      ep0   = te(medium)/ekin
      g2    = t0*t0/e02
      g3    = (2.0*t0 + 1.0)/e02
      gmax  = (1.0 + 1.25*g2)
                                            ! Sample/reject to obtain br
1     continue
        call randomset(rnnow)
        br = te(medium)/(ekin - extrae*rnnow)
         r = br/(1.0 - br)

        ! Decide whether or not to accept
        call randomset(rnnow)
        rejf = 1.0 + g2*br*br + r*(r - g3)
        rnnow = gmax*rnnow
        if (rnnow .gt. rejf) go to 1
                                                  ! Divide up the energy
      ekse2 = br*ekin
      ese1 = eie - ekse2
      ese2 = ekse2 + RM

      e(np) = ese1
      e(np+1) = ese2
                                            ! Moller angles are uniquely
                                            ! determined by kinematics
      h1 = (eie + RM)/ekin
      dcosth = h1*(ese1 - RM)/(ese1 + RM)
      sinthe = sqrt(1.0 - dcosth)
      costhe = sqrt(dcosth)
      call uphi(2,1)                             ! Set direction cosines
      np = np + 1
      iq(np)=-1
      dcosth = h1*(ese2 - RM)/(ese2 + RM)
      sinthe = - sqrt(1.0 - dcosth)
      costhe = sqrt(dcosth)
      call uphi(3,2)                             ! Set direction cosines
      k1step(np) = 0.
      k1init(np) = 0.
      k1rsd(np) = 0.


      call npconv
      atmrc(2,4) = atmrc(2,4) + 1.d0




      if (impacr(ir(np)) .eq. 1.and. iedgfl(ir(np)) .ne. 0) then
        call randomset(rnnow)
        eke = eie - RM
        elke2 = log(eke)
        lelke2 = eico1(medium)*elke2 + eico0(medium)
        do ifun=1,nepm(medium)
      
          if (rnnow .lt. eiir) then

            iarg = 25                                  ! BEFORE call eii

            if (iausfl(iarg+1) .ne. 0) call ausgab(iarg)

            iz = zelem(medium,ifun)

            call eii

            iarg = 26                                   ! AFTER call eii

            if (iausfl(iarg+1) .ne. 0) call ausgab(iarg)

            return
          end if
        end do
      end if

      if(mstz(85).eq.99)then ! debug output
         write(93,*)'---- end of moller -----'
         do i=1,np
            
         enddo
      endif
                                                      ! ----------------
      return                                          ! Return to ELECTR
                                                      ! ----------------
      end

      subroutine mscat
      use EGS5_MS_MOD !FURUTA20140825
      USE egs5_elecin_mod
      USE egs5_media_mod
      implicit none
      include 'include/egs5_h.f'               ! Main EGS5 "header" file
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
      real*8 rms1,rms2,rms3,rms4,rms5,rms6,rms7,rms8         ! Arguments
      integer i21,i22,i31,i32
      integer iegrid, iamu, iprt, ik1, im
      real*8 decade,delog,demod, xmu, xi,  b1,c1, x1,x2, fject,fmax
      real*8 ktot
      common /tstara/ atmrc(10,8)      ! Ogawa 2023/7/18, (10,7) -> (10,8) for plasmon
      real*8 atmrc
      imscat = imscat + 1                  ! Count entry into subroutine
      im = medium
      if(useGSD(im).ne.0 .and. eke.lt.msgrid(nmsgrd(im),im)) then
        if(iq(np) .eq. -1) then
          iprt = 1
        else
          iprt = 2
        endif
        !-->  find the correct energy interval
        delog = DLOG10(eke*1.d6)
        demod = MOD(delog,1.d0)
        decade = delog - demod
        iegrid = nmsdec(im) * (decade - initde(im) + demod) - jskip(im)
        if(iegrid .ge. nmsgrd(im)) iegrid = nmsgrd(im) - 1
        !-->  randommly select which energy point to use
        call randomset(xi)
        if(xi .lt. fject) iegrid = iegrid +1
        !-->  find the correct K1 interval
        ktot = k1rsd(np) + k1init(np)
        ik1 = DLOG(ktot/k1grd(iprt,1)) / dk1log(iprt) + 1
        if(ik1 .gt. NK1) then
          ik1 = NK1-1
        else if(ik1 .le. 0) then
          ik1 = 1
        endif
        call randomset(xi)
        if(xi .lt. fject) ik1 = ik1 +1
        !-->  first check for no-scatter probability
        call randomset(xi)
        if(xi.lt.pnoscat(iprt,iegrid,ik1,im)) return
        !-->  get the angular interval
        call randomset(xi)
        iamu = xi * neqp(im) + 1
        !-->  if we're in the last bin, get the sub-bin number
        if(iamu .eq. neqp(im)) then
          call randomset(xi)
          if(xi .lt. ecdf(2,iprt,iegrid,ik1,im)) then
            iamu = 1
          else
            call findi(ecdf(1,iprt,iegrid,ik1,im),xi,neqa(im)+1,iamu)
          endif
          iamu = iamu + neqp(im) - 1
        endif
        b1 = ebms(iamu,iprt,iegrid,ik1,im)
        eta = eetams(iamu,iprt,iegrid,ik1,im)
        x1 = eamu(iamu,iprt,iegrid,ik1,im)
        x2 = eamu(iamu+1,iprt,iegrid,ik1,im)
        c1 = (x2 + eta) / (x2 - x1)
        fmax = 1.d0 + .25d0 * b1 * (x2 - x1)**2
        !-->  rejection loop
 6      continue
          !-->  sample Wentzel shape part of fit
          call randomset(xi)
          xmu = ((eta * xi) + (x1 * c1)) / (c1 - xi)
          !-->  rejection test
          fject = 1.d0 + b1 * (xmu - x1) * (x2 - xmu)
          call randomset(xi)
          if(xi * fmax .gt. fject) go to 6
        costhe = 1.d0 - 2.d0 * xmu
        sinthe =  DSQRT(1.d0 - costhe * costhe)
        iskpms = 0
        if( iq(np) .eq. -1 ) atmrc(2,6) = atmrc(2,6) + 1.d0
        if( iq(np) .eq.  1 ) atmrc(3,6) = atmrc(3,6) + 1.d0
        return
      !--> user or lower limit initiated skip of MS
      else if (nomsct(ir(np)).eq.1 .or. iskpms.ne.0) then
        sinthe = 0.
        costhe = 1.
        theta = 0.
        noscat = noscat + 1
        iskpms = 0
        if( iq(np) .eq. -1 ) atmrc(2,6) = atmrc(2,6) + 1.d0
        if( iq(np) .eq.  1 ) atmrc(3,6) = atmrc(3,6) + 1.d0
        return
      end if
      xr = sqrt(gms*tmstep*b)
      if (b .gt. 2.) then
        bi = 1./b
      else
        bi = 0.5
      end if
      bmd = 1. + 1.75*bi
      bm1 = (1. - 2./b)/bmd
      bm2 = (1. + 0.025*bi)/bmd
                 ! -----------------------------------------------------
 1    continue   ! Loop for Bethe correction factor (or other) rejection
                 ! -----------------------------------------------------
        call randomset(rms1)
        if (rms1 .le. bm1) then                           ! Gaussian, F1
          call randomset(rms2)
          if (rms2 .eq. 0.) then
            rms2 = 1.E-30
          end if
          thr = sqrt(max(0.D0,-log(rms2)))
        else if (rms1 .le. bm2) then                          ! Tail, F3
          call randomset(rms3)
          call randomset(rms4)
          eta = max(rms3,rms4)
          i31 = b0g31 + eta*b1g31
          g31 = g310(i31) + eta*(g311(i31) + eta*g312(i31))
          i32 = b0g32 + eta*b1g32
          g32 = g320(i32) + eta*(g321(i32) + eta*g322(i32))
          g3 = g31 + g32*bi                         ! Rejection function
          call randomset(rms5)
          if (rms5 .gt. g3) go to 1
          thr = 1./eta
        else           ! Central correction, F2
          call randomset(rms6)
          thr = rms6
          i21 = b0g21 + thr*b1g21
          g21 = g210(i21) + thr*(g211(i21) + thr*g212(i21))
          i22 = b0g22 + thr*b1g22
          g22 = g220(i22) + thr*(g221(i22) + thr*g222(i22))
          g2  = g21 + g22*bi                        ! Rejection function
          call randomset(rms7)
          if (rms7 .gt. g2) go to 1
        end if
        theta = thr*xr           ! Real angle (thr is the reduced angle)
        if (theta .ge. PI) go to 1
        sinthe = sin(theta)
        call randomset(rms8)
        if (rms8**2*theta .le. sinthe) go to 2
      go to 1
 2    continue
      cthet = PI5D2 - theta
      costhe = sin(cthet)
      if( iq(np) .eq. -1 ) atmrc(2,6) = atmrc(2,6) + 1.d0
      if( iq(np) .eq.  1 ) atmrc(3,6) = atmrc(3,6) + 1.d0
                                                      ! ----------------
      return                                          ! Return to ELECTR
                                                      ! ----------------
      end
      subroutine npconv
      implicit none
      include 'include/egs5_h.f'               ! Main EGS5 "header" file
      include 'include/egs5_stack.f'
      integer iqdum,irdum,latchdum,i,j
      edum=e(np)
      xdum=x(np)
      ydum=y(np)
      zdum=z(np)
      udum=u(np)
      vdum=v(np)
      wdum=w(np)
      wtdum=wt(np)
      ufdum=uf(np)
      vfdum=vf(np)
      wfdum=wf(np)
      k1sdum=k1step(np)
      k1idum=k1init(np)
      k1rdum=k1rsd(np)
      dndum=dnear(np)
      tdum=time(np)
      iqdum=iq(np)
      irdum=ir(np)
      latchdum=latch(np)
      do i=1,np-1
        j=np+1-i         ! move 1-(np-1) to 2-np
        e(j)=e(j-1)
        x(j)=x(j-1)
        y(j)=y(j-1)
        z(j)=z(j-1)
        u(j)=u(j-1)
        v(j)=v(j-1)
        w(j)=w(j-1)
        wt(j)=wt(j-1)
        uf(j)=uf(j-1)
        vf(j)=vf(j-1)
        wf(j)=wf(j-1)
        k1step(j)=k1step(j-1)
        k1init(j)=k1init(j-1)
        k1rsd(j)=k1rsd(j-1)
        dnear(j)=dnear(j-1)
        time(j)=time(j-1)
        iq(j)=iq(j-1)
        ir(j)=ir(j-1)
        latch(j)=latch(j-1)
      end do
      e(1)=edum
      x(1)=xdum
      y(1)=ydum
      z(1)=zdum
      u(1)=udum
      v(1)=vdum
      w(1)=wdum
      wt(1)=wtdum
      uf(1)=ufdum
      vf(1)=vfdum
      wf(1)=wfdum
      k1step(1)=k1sdum
      k1init(1)=k1idum
      k1rsd(1)=k1rdum
      dnear(1)=dndum
      time(1)=tdum
      iq(1)=iqdum
      ir(1)=irdum
      latch(1)=latchdum
      return
      end
      subroutine npshift(nshift)
      implicit none
      include 'include/egs5_h.f'               ! Main EGS5 "header" file
      include 'include/egs5_stack.f'
      integer nshift
      integer i,j,nnp,npola,npold
      npold=np
      nnp=np-nshift
      do j=1,nnp
        e(npold+j)=e(j)
        x(npold+j)=x(j)
        y(npold+j)=y(j)
        z(npold+j)=z(j)
        u(npold+j)=u(j)
        v(npold+j)=v(j)
        w(npold+j)=w(j)
        wt(npold+j)=wt(j)
        uf(npold+j)=uf(j)
        vf(npold+j)=vf(j)
        wf(npold+j)=wf(j)
        k1step(npold+j)=k1step(j)
        k1init(npold+j)=k1init(j)
        k1rsd(npold+j)=k1rsd(j)
        dnear(npold+j)=dnear(j)
        time(npold+j)=time(j)
        iq(npold+j)=iq(j)
        ir(npold+j)=ir(j)
        latch(npold+j)=latch(j)
      end do
      do j=1,npold
        e(j)=e(j+nnp)
        x(j)=x(j+nnp)
        y(j)=y(j+nnp)
        z(j)=z(j+nnp)
        u(j)=u(j+nnp)
        v(j)=v(j+nnp)
        w(j)=w(j+nnp)
        wt(j)=wt(j+nnp)
        uf(j)=uf(j+nnp)
        vf(j)=vf(j+nnp)
        wf(j)=wf(j+nnp)
        k1step(j)=k1step(j+nnp)
        k1init(j)=k1init(j+nnp)
        k1rsd(j)=k1rsd(j+nnp)
        dnear(j)=dnear(j+nnp)
        time(j)=time(j+nnp)
        iq(j)=iq(j+nnp)
        ir(j)=ir(j+nnp)
        latch(j)=latch(j+nnp)
      end do
      return
      end
      subroutine pair
      USE egs5_brempr_mod
      implicit none
      include 'include/egs5_h.f'               ! Main EGS5 "header" file
      include 'include/egs5_stack.f'
      include 'include/egs5_uphiot.f'
      include 'include/egs5_useful.f'
      include 'include/counters.f'       ! Additional (non-EGS5) COMMONs
      real*8 rnnow,rnnow1,rnnow2,rnnow3                      ! Arguments
      integer lvx,lvl0,lvl,ichrg,it
      common /tstara/ atmrc(10,8)      ! Ogawa 2023/7/18, (10,7) -> (10,8) for plasmon
      real*8 atmrc
      ipair = ipair + 1                    ! Count entry into subroutine
      eig = e(np)
      if (eig .le. 2.1) then               ! Below 2.1 MeV (approximate)
        call randomset(rnnow)                 ! KEK method for smoothing
        ese2 = RM + rnnow*(eig/2. - RM)       !   connection at boundary
      else
        if (eig .lt. 50.) then                  ! Above 2.1 MeV - sample
          lvx = 1
          lvl0 = 0
        else
          lvx = 2
          lvl0 = 3
        end if
                                 ! Start of main sampling-rejection loop
1       continue
          call randomset(rnnow1)
          call randomset(rnnow)
                                              ! Start of 12(br-0.5)**2
                                              ! subdistribution sampling
          if (rnnow .ge. bpar(lvx,medium)) then
            lvl = lvl0 + 1
            call randomset(rnnow2)
            call randomset(rnnow3)
            br = 0.5*(1.0 - max(rnnow1,rnnow2,rnnow3))
                             ! Start of uniform subdistribution sampling
          else
            lvl = lvl0 + 3
            br = rnnow1*0.5
          end if
                                                 ! Check that br, Adelta
                                                 ! and Cdelta > 0
          if(eig*br .lt. RM) go to 1
          del = 1.0/(eig*br*(1.0 - br))
          if (del .ge. delpos(lvx,medium)) go to 1
          delta = delcm(medium)*del
          if (delta .lt. 1.0) then
          else
          end if
          call randomset(rnnow)
          if (rnnow .gt. rejf) go to 1
        ese2 = br*eig
      end if
                                        ! Set up secondary electron #1
                                        ! (electron #2 has lower energy)
      ese1 = eig - ese2
      e(np) = ese1
      e(np+1) = ese2
      k1step(np) = 0.
      k1init(np) = 0.
      k1rsd(np) = 0.
      k1step(np+1) = 0.
      k1init(np+1) = 0.
      k1rsd(np+1) = 0.
                                            ! Sample to get polar angles
                                            ! of secondary electrons
                              ! Sample lowest-order angular distribution
      if ((iprdst .eq. 1) .or.&
          ((iprdst .eq. 2) .and. (eig .lt. 4.14))) then
        do ichrg=1,2
          if (ichrg .eq. 1) then
            ese = ese2
          else
            ese = ese1
          end if
          pse = sqrt(max(0.D0,(ese - RM)*(ese + RM)))
          call randomset(rnnow)
          costhe = 1.0 - 2.0*rnnow
          sinthe = RM*sqrt((1.0 - costhe)*(1.0 + costhe))/&
                  (pse*costhe + ese)
          costhe = (ese*costhe + pse)/(pse*costhe + ese)
          if (ichrg .eq. 1) then
            call uphi(2,1)
          else
            np = np + 1
            sinthe = -sinthe
            call uphi(3,2)
          end if
        end do
        call randomset(rnnow)
        if (rnnow .le. 0.5) then
          iq(np) = 1
          iq(np-1) = -1
        else
          iq(np) = -1
          iq(np-1) = 1
        end if
        go to 100
                                           ! Sample from Motz-Olsen-Koch
                                           ! (1969) distribution
      else if ((iprdst .eq. 2) .and.&
              (eig .ge. 4.14)) then
        ztarg = zbrang(medium)
        tteig = eig/RM
        do ichrg=1,2
          if (ichrg .eq. 1) then
            ese = ese2
          else
            ese = ese1
          end if
          ttese = ese/RM
          ttpse = sqrt((ttese - 1.0)*(ttese + 1.0))
          esedei = ttese/(tteig - ttese)
          eseder = 1.0/esedei
          ximin = 1.0/(1.0 + (PI*ttese)**2)
          ya = (2.0/tteig)**2
          xitry = max(0.01D0,max(ximin,min(0.5D0,sqrt(ya/ztarg))))
          galpha = 1.0 + 0.25*log(ya + ztarg*xitry**2)
          gbeta = 0.5*ztarg*xitry/(ya + ztarg*xitry**2)
          galpha = galpha - gbeta*(xitry - 0.5)
          ximid = galpha/(3.0*gbeta)
          if (galpha .ge. 0.0) then
            ximid = 0.5 - ximid + sqrt(ximid**2 + 0.25)
          else
            ximid = 0.5 - ximid - sqrt(ximid**2+0.25)
          end if
          ximid = max(0.01D0,max(ximin,min(0.5D0,ximid)))
          rejtop = 1.02*max(rejmin,rejmid)
2         continue
            call randomset(xitst)
            call randomset(rtest)
            theta = sqrt(1.0/xitst - 1.0)/ttese
            if ((rtest .gt. (rejtst/rejtop) .or.&
               (theta .ge. PI))) go to 2
          sinthe=sin(theta)
          costhe=cos(theta)
          if (ichrg .eq. 1) then
            call uphi(2,1)
          else
            np = np+1
            sinthe = -sinthe
            call uphi(3,2)
          end if
          end do
        call randomset(rnnow)
        if (rnnow .le. 0.5) then
          iq(np) = 1
          iq(np-1) = -1
        else
          iq(np) = -1
          iq(np-1) = 1
        end if
        go to 100
      ! Polar angle is m/k (default)
      else
        theta=RM/eig
      end if
      call uphi(1,1)             ! Set direction cosines for electron #1
      np = np + 1
      sinthe = -sinthe
      call uphi(3,2)             ! Set direction cosines for electron #2
                          ! Randomly decide which particle is "positron"
      call randomset(rnnow)
      if (rnnow .le. 0.5) then
        iq(np) = 1
        iq(np-1) = -1
      else
        iq(np) = -1
        iq(np-1) = 1
      end if
100   t = e(np)
      e(np) = e(np-1)
      e(np-1) = t
      it = iq(np)
      iq(np) =iq(np-1)
      iq(np-1)=it
      t = u(np)
      u(np) = u(np-1)
      u(np-1) = t
      t = v(np)
      v(np) = v(np-1)
      v(np-1) = t
      t = w(np)
      w(np) = w(np-1)
      w(np-1) = t
      call npconv
      call npconv
      atmrc(1,5) = atmrc(1,5) + 1.d0
                                                      ! ----------------
      return                                          ! Return to PHOTON
                                                      ! ----------------
      end
      subroutine photo
      USE egs5_brempr_mod
      USE egs5_media_mod
      USE egs5_edge_mod
      implicit none
      include 'include/egs5_h.f'               ! Main EGS5 "header" file
      include 'include/egs5_bounds.f'    ! COMMONs required by EG55 code
      include 'include/egs5_epcont.f'
      include 'include/egs5_misc.f'
      include 'include/egs5_stack.f'
      include 'include/egs5_uphiot.f'
      include 'include/egs5_useful.f'
      include 'include/egs5_uservr.f'
      include 'include/egs5_userxt.f'
      include 'include/counters.f'       ! Additional (non-EGS5) COMMONs
      real*8 rnnow                                           ! Arguments
      integer iarg
      common /tstara/ atmrc(10,8)      ! Ogawa 2023/7/18, (10,7) -> (10,8) for plasmon
      real*8 atmrc
      iphoto = iphoto + 1                  ! Count entry into subroutine
      nxray = 0
      nauger = 0
      irl = ir(np)
      eig = e(np)
      phol = log(eig)
      medium = med(irl)
      eigk = eig*1000.D0
      pholk = log(eigk)
      pholk2 = pholk*pholk
      pholk3 = pholk2*pholk
      total = 0.
      do i=1,nne(medium)          ! Shell-wise photoelectric calculation
        iz = zelem(medium,i)
        if (eigk .le. eedge(1,iz)) then
          crosk(i) = 0.
        else
        end if
        if (pm0(2,iz) .eq. 0. .or. eigk .le. eedge(2,iz)) then
          crosl1(i) = 0.
        else
        end if
        if (pm0(3,iz) .eq. 0. .or. eigk .le. eedge(3,iz)) then
          crosl2(i) = 0.
        else
        end if
        if (pm0(4,iz) .eq. 0. .or. eigk .le. eedge(4,iz)) then
          crosl3(i) = 0.
        else
        end if
        if (eigk .le. embind(iz)) then
          tcros(i) = 0.
        else
          if (pm0(5,iz) .eq. 0.) then
            crosm(i) = 0.
          else
          end if
        end if
        tcros(i) = tcros(i)*pz(medium,i)
        total = total + tcros(i)
      end do
      if (total .eq. 0.) then            ! Below M-edge for all elements
        edep = eig
        iarg = 4                    ! Deposit all of the photon's energy
        if (iausfl(iarg+1) .ne. 0) call ausgab(iarg)
        e(np) = 0.
        return
      end if
      if (nne(medium) .eq. 1) then
        iz = zelem(medium,1)
        noel = 1
        go to 1
      end if
      do i=1,nne(medium)-1
        if (i .eq. 1) then
          pbran(i) = tcros(i)/total
        else
          pbran(i) = pbran(i-1) + tcros(i)/total
        end if
      end do
      call randomset(rnnow)
      do i=1,nne(medium)-1
        if (rnnow .le. pbran(i)) then
          iz = zelem(medium,i)
          noel = i
          go to 1
        end if
      end do
      iz = zelem(medium,nne(medium))
      noel = nne(medium)
 1    continue                  ! Determine K, L x-rays for each element
      if (eigk .le. eedge(4,iz)) then        ! Below L3 edge, treat as M
        ebind = embind(iz)*1.D-3
        edep = ebind
        go to 2                                          ! Below L3 edge
      end if
      call randomset(rnnow)                     ! Sample to decide shell
      if (rnnow .gt. bshl3(noel)) then   ! M,N,...absorption, treat as M
        ebind = embind(iz)*1.D-3
        edep = ebind
        go to 2
      else if(rnnow .le. bshk(noel)) then                 ! K absorption
        call kshell
        ebind = eedge(1,iz)*1.D-3
      else if(rnnow .le. bshl1(noel)) then               ! L1 absorption
        call lshell(1)
      else if(rnnow .le. bshl2(noel)) then               ! L2 absorption
        call lshell(2)
      else                                               ! L3 absorption
        call lshell(3)
      end if
      if (iedgfl(irl) .le. 0) nxray = 0
      if (iauger(irl) .le. 0) nauger = 0
      edep = ebind
      if (nxray .ge. 1) then
        do iphot=1,nxray
          edep = edep - exray(iphot)
        end do
      end if
      if (nauger .ge. 1) then
        do ielec=1,nauger
          edep = edep - eauger(ielec)
        end do
      end if
      if (edep .lt. 0.) edep = 0. ! To avoid numerical precision problem
 2    continue
      ecapbind=ecapbind+edep  ! dead electron energy
      e(np) = edep
      edep = 0.0 ! reset edep, T.Sato 2016/02/01
      iarg = 4                     ! Part of binding energy is deposited
      if (iausfl(iarg+1) .ne. 0) call ausgab(iarg)
      iq(np) = -1                        ! Photoelectron (always set up)
      e(np) = eig - ebind + RM
      if (iphter(ir(np)) .eq. 1) then   ! Select photoelectron direction
        eelec = e(np)
        if (eelec .gt. ecut(ir(np))) then
          beta = sqrt((eelec - RM)*(eelec + RM))/eelec
          gamma = eelec/RM
          alpha = 0.5D0*gamma - 0.5D0 + 1.D0/gamma
          ratio = beta/alpha
 3        continue
            call randomset(rnnow)
            rnpht = 2.D0*rnnow - 1.D0
            if (ratio .le. 0.2D0) then
              fkappa = rnpht + 0.5D0*ratio*(1.D0 - rnpht)*(1.D0 + rnpht)
              costhe = (beta + fkappa)/(1.D0 + beta*fkappa)
              xi = 1.D0/(1.D0 - beta*costhe)
            else
              costhe = (1.D0 - 1.D0/xi)/beta
            end if
            sinth2 = max(0.D0,(1.D0 - costhe)*(1.D0 + costhe))
            call randomset(rnnow)
            if(rnnow .le. 0.5D0*(1.D0 + gamma)*sinth2*xi/gamma) go to 4
          go to 3
 4        continue
          sinthe = sqrt(sinth2)
          call uphi(2,1)
        end if
      end if
      call npconv
      atmrc(1,3) = atmrc(1,3) + 1.d0
      if (nauger .ne. 0) then                   ! Set up Auger electrons
        do ielec=1,nauger
          np = np + 1
          e(np) = eauger(nauger+1-ielec) + RM
          iq(np) = -1
          call randomset(rnnow)
          costhe = 2.D0*rnnow - 1.D0
          sinthe = sqrt(1.D0 -costhe*costhe)
          u(np) = 0.
          v(np) = 0.
          w(np) = 1.D0
          call uphi(2,1)
          x(np) = x(1)
          y(np) = y(1)
          z(np) = z(1)
          uf(np) = 0.
          vf(np) = 0.
          wf(np) = 0.
          ir(np) = ir(1)
          wt(np) = wt(1)
          time(np) = time(1)
          dnear(np) = dnear(1)
          latch(np) = latch(1)
          k1step(np) = 0.
          k1init(np) = 0.
          k1rsd(np) = 0.
          call npconv
        end do
        atmrc(1,2) = atmrc(1,2) + 1.d0
      end if
      if (nxray .ne. 0) then                ! Set up fluorescent photons
        do iphot=1,nxray
          np = np + 1
          e(np) = exray(nxray+1-iphot)
          iq(np) = 0
          call randomset(rnnow)
          costhe = 2.D0*rnnow - 1.D0
          sinthe = sqrt(1.D0 - costhe*costhe)
          u(np) = 0.
          v(np) = 0.
          w(np) = 1.D0
          call uphi(2,1)
          x(np) = x(1)
          y(np) = y(1)
          z(np) = z(1)
          uf(np) = 0.
          vf(np) = 0.
          wf(np) = 0.
          ir(np) = ir(1)
          wt(np) = wt(1)
          time(np) = time(1)
          dnear(np) = dnear(1)
          latch(np) = latch(1)
          k1step(np) = 0.
          k1init(np) = 0.
          k1rsd(np) = 0.
          call npconv
        end do
        atmrc(1,1) = atmrc(1,1) + 1.d0
      end if
      if(nauger+nxray.gt.0) then
      end if
                                                      ! ----------------
      return                                          ! Return to PHOTON
                                                      ! ----------------
      end
      subroutine raylei
      USE egs5_photin_mod !<- 2015.08xx allocatable
      USE egs5_thresh_mod !<-
      implicit none
      include 'include/egs5_h.f'               ! Main EGS5 "header" file
      include 'include/egs5_misc.f'      ! COMMONs required by EGS5 code
      include 'include/egs5_stack.f'
      include 'include/egs5_uphiot.f'
      include 'include/egs5_useful.f'
      include 'include/counters.f'       ! Additional (non-EGS5) COMMONs
      real*8 rnnow                                           ! Arguments
      real*8 x2,q2,csqthe,rejf,br                      ! Local variables
      integer lxxx
      common /tstara/ atmrc(10,8)      ! Ogawa 2023/7/18, (10,7) -> (10,8) for plasmon
      real*8 atmrc
      iraylei = iraylei + 1                ! Count entry into subroutine
 1    continue
        call randomset(rnnow)
        lxxx = rco1(medium)*rnnow + rco0(medium)
        x2 = rsct1(lxxx,medium)*rnnow + rsct0(lxxx,medium)
        q2 = x2*RMSQ/(20.60744*20.60744)
        costhe = 1.-q2/(2.*e(np)*e(np))
        if (abs(costhe) .gt. 1.) go to 1
        csqthe = costhe*costhe
        rejf = (1. + csqthe)/2.
        call randomset(rnnow)
        if (rnnow .le. rejf) go to 2
      go to 1
 2    continue
      sinthe = sqrt(1. - csqthe)
      if (lpolar(ir(np)) .eq. 0) then
        call uphi(2,1)
      else
        br = 1.
        call aphi(br)
        call uphi(3,1)
      end if
      atmrc(1,6) = atmrc(1,6) + 1.d0   ! S.Abe 2017/02/08
                                                      ! ----------------
      return                                          ! Return to PHOTON
                                                      ! ----------------
      end
      subroutine rk1
      USE egs5_brempr_mod
      USE egs5_elecin_mod
      USE egs5_media_mod
      USE egs5_mscon_mod
      USE egs5_scpw_mod
      USE egs5_thresh_mod
      implicit none
      include 'include/egs5_h.f'         ! COMMONs required by EGSD5 code
      include 'include/egs5_useful.f'
      include 'include/egs5_misc.f'
      common /paran/ icfn(100), ilfn(100), chfn(100)
      character chfn*200
      integer icfn, ilfn
      integer nmatk1,nek1(MXEPERMED)
      integer nz2, idx, iek1, iedx, niek, neke, iz2
      integer i,im,ie,j,n, m
      read(517,*) nmatk1
      do i = 1, nmatk1
        read(517,*) z, z2k1w(i), watot, rhok1(i)
        z2k1w(i) = z2k1w(i) / watot
        read(517,*) nek1(i)
        do j = nek1(i), 1, -1
          ! prep for interp in D * rho
          bk1(j,i) = bk1(j,i) - slopek1(j,i) * DLOG(rhok1(i))
          dhighk1(j,i) = dhighk1(j,i) * rhok1(i)
          dlowk1(j,i) = dlowk1(j,i) * rhok1(i)
        end do
      end do
      close(517)
      do 10 im = 1, nmed
        if(charD(im).eq.0.d0) go to 10
        lcharD = dlog(charD(im)*rhom(im))
        watot = 0.d0
        z2w = 0.d0
        do ie = 1, nne(im)
          z2w = z2w + pz(im,ie) * zelem(im,ie) * (zelem(im,ie) + 1.d0)
          watot = watot + pz(im,ie) * wa(im,ie)
        end do
        z2w = z2w / watot
        if(z2w.gt.z2k1w(nmatk1)) then
          nz2 = 1
          iz2 = nmatk1
          zfrac = 0.d0
        else if(z2w.lt.z2k1w(1)) then
          nz2 = 1
          iz2 = 1
          zfrac = 0.d0
        else
          nz2 = 2
          call findi(z2k1w,z2w,nmatk1,iz2)
          zfrac = (z2w - z2k1w(iz2)) / (z2k1w(iz2+1) - z2k1w(iz2))
        end if
        !-->  get scpow at max E of the reference materials,
        !-->  in case we have to extrapolate
        do  n = 1, nz2
          idx = iz2 + n - 1
          eke = ek0k1(nek1(idx),idx)
          if(ue(im)-RM.gt.eke) then
            elke = log(eke)
            j = eke0(im) + elke * eke1(im)
            scpeMax(n) = escpw1(j,im)*elke + escpw0(j,im)
          else
            scpeMax(n) = 0.0
          end if
        end do
        neke = meke(im)
        do j = 1,neke-1
          if(j.eq.1) then
            eke = ae(im) - RM
            elke = log(eke)
          else if(j.eq.neke-1) then
            eke = ue(im) - RM
            elke = log(eke)
          else
            elke = (j + 1 - eke0(im)) / eke1(im)
            eke = DEXP(elke)
          end if
          scpe = escpw1(j,im)*elke + escpw0(j,im)
          scpp = pscpw1(j,im)*elke + pscpw0(j,im)
          scprat = scpp/scpe
          do n = 1, nz2
            extrape = 1.d0
            idx = iz2 + n - 1
            if(eke.gt.ek0k1(nek1(idx),idx)) then
              niek = 1
              iek1 = nek1(idx)
              frace = 0.d0
              extrape = scpe/scpeMax(n)
            else if(eke.lt.ek0k1(1,idx)) then
              niek = 1
              iek1 = 1
              frace = 0.d0
            else
              niek = 2
              call findi(ek0k1(1,idx),eke,nek1(idx),iek1)
            end if
            do m = 1, niek
              iedx = iek1 + m - 1
              if(charD(im)*rhom(im).gt.dhighk1(iedx,idx)) then
                k1ez(m) = k1high(iedx,idx)
              else if(charD(im)*rhom(im).lt.dlowk1(iedx,idx)) then
                k1ez(m) = k1low(iedx,idx)
              else
                k1ez(m) = DEXP(slopek1(iedx,idx)*lcharD + bk1(iedx,idx))
              end if
              k1sez(m) =  k1low(iedx,idx)
            end do
            k1z(n) = k1ez(1) + frace * (k1ez(2) - k1ez(1))
            k1z(n) = extrape * k1z(n)
            k1sz(n) = k1sez(1) + frace * (k1sez(2) - k1sez(1))
            k1sz(n) = extrape * k1sz(n)
          end do
          k1new = k1z(1) + zfrac * (k1z(2) - k1z(1))
          k1s = k1sz(1) + zfrac * (k1sz(2) - k1sz(1))
          if(k1new.gt.k1maxe) then
            k1maxe = k1new
          else if(k1new.lt.k1mine) then
            k1mine = k1new
          end if
          if(k1new*scprat.gt.k1maxp) then
            k1maxp = k1new*scprat
          else if(k1new*scprat.lt.k1minp) then
            k1minp = k1new*scprat
          end if
          if(k1s.lt.k1mine) then
            k1mine = k1s
          end if
          if(k1s*scprat.lt.k1minp) then
            k1minp = k1s*scprat
          end if
          if(j.gt.1) then
            delke = elke - elkeold
            ekini1(j,im) =  (k1new - k1old) / delke
            ekini0(j,im) = (k1old * elke - k1new * elkeold) / delke
            pkini1(j,im) = ekini1(j,im) * scprat
            pkini0(j,im) = ekini0(j,im) * scprat
            ek1s1(j,im) =  (k1s - k1sold) / delke
            ek1s0(j,im) = (k1sold * elke - k1s * elkeold) / delke
            pk1s1(j,im) = ek1s1(j,im) * scprat
            pk1s0(j,im) = ek1s0(j,im) * scprat
          end if
          elkeold = elke
          k1old = k1new
          k1sold = k1s
        end do
        ekini1(1,im) = ekini1(2,im)
        ekini0(1,im) = ekini0(2,im)
        pkini1(1,im) = pkini1(2,im)
        pkini0(1,im) = pkini0(2,im)
        ekini1(neke,im) = ekini1(neke-1,im)
        ekini0(neke,im) = ekini0(neke-1,im)
        pkini1(neke,im) = pkini1(neke-1,im)
        pkini0(neke,im) = pkini0(neke-1,im)
        ek1s1(1,im) = ek1s1(2,im)
        ek1s0(1,im) = ek1s0(2,im)
        pk1s1(1,im) = pk1s1(2,im)
        pk1s0(1,im) = pk1s0(2,im)
        ek1s1(neke,im) = ek1s1(neke-1,im)
        ek1s0(neke,im) = ek1s0(neke-1,im)
        pk1s1(neke,im) = pk1s1(neke-1,im)
        pk1s0(neke,im) = pk1s0(neke-1,im)
10    continue
      do j = 1, nreg
        if(k1Lscl(j).gt.0.d0 .and. k1Hscl(j).gt.0.d0) then
          im = med(j)
          c1 = k1Hscl(j) - dlog(ue(im) - RM) * c2
          k1Lscl(j) = c1
          k1Hscl(j) = c2
        else if(k1Lscl(j).ne.0.d0) then
          k1Lscl(j) = 0.d0
        end if
      end do
      return
      end
      subroutine rmsfit
      use EGS5_MS_MOD !FURUTA20140825
      USE egs5_media_mod
      USE egs5_mscon_mod
      USE egs5_thresh_mod
      implicit none
      include 'include/egs5_h.f'               ! Main EGS5 "header" file
      include 'include/egs5_misc.f'
      include 'include/egs5_usersc.f'
      include 'err.inc'
      integer i,j,k,ib,il,iprt,ik1,idummy, nfmeds, hasGS
      integer nm, lmdl, lmdn
      integer , allocatable :: lok(:) !<- lok(MXMED)
      integer imxmed, flg_local
      character buffer(72),mdlabl(8)
      logical openMS, useGS, doGS, GSFile
      common /paran/ icfn(100), ilfn(100), chfn(100)
      character chfn*200
      integer icfn, ilfn
      openMS = .false.
      useGS = .false.
      do k=1,nmed
        if(charD(k).eq.0.d0) then
          openMS = .true.
        endif
        if(useGSD(k).eq.0) then
          nmsgrd(k) = 1
          msgrid(1,k) = 0.d0
        else
          useGS = .true.
        endif
      end do
      if(.not.(useGS .or. openMS)) then
        return
      endif
      if(useGS) then
        do k=1,nmed
          if(useGSD(k).eq.0) then
            write(506,1001) k
            useGSD(k) = 1
          endif
        end do
      endif
      GSFile = .false.
      doGS=.false.
      GSFile = .true.
      go to 15
      !  all this is for corrupt gsdist.dat files
13    close(517)
14    open(UNIT=517,FILE=chfn(23)(1:ilfn(23))//'.msfit',STATUS='old')
15    read(517,*,end=13) nfmeds
      do 18 j = 1, nfmeds
        read(517,'(72a1)',end=1000) buffer
        read(517,*) hasGS, charD0, efrch0, efrcl0, ue0, ae0
        do 16 k=1,nmed
          do ib=1,lmdn
            il = lmdl + ib
            if (buffer(il) .ne. media(ib,k)) go to 16
            if (ib .eq. lmdn) go to 17
          end do
16      continue
        go to 18
        !  here if we are using this material
17      if(useGS .and. hasGS.eq.0) then
          write(506,1002) k
          doGS = .true.
        endif
        if(charD(k).eq.0.d0) then
          if(charD0.ne.0.d0) then
            if(useGS) then
              write(506,1003) k, charD0
            else
              write(506,1004) k
              ErrCha = ''
              ErrID = 'L:4901/R:rmsfit/F:egs5.f' !E86_018_001
              call ErrWrite(ErrID,ErrCha)
              write(*,1004) k
              stop
            endif
          else
            write(506,1005) k, efrch0, efrcl0
          endif
        endif
        !  make sure previous data was for less than or equal ue/ae
        if(useGS .and. (ue0.lt.ue(j) .or. ae0.gt.ae(j))) then
          write(506,1006) k, ae0, ue0, ae(j), ue(j)
          doGS = .true.
        endif
18    continue
      if(.not.useGS) then
        close(517)
        return
      endif
      do j = 1, nreg
        k = med(j)
        if(k.ne.0) then
          if( rhom(k).ne.rhor(j)) then
            write(506,1007) j,k
            ErrCha = ''
            ErrID = 'L:4960/R:rmsfit/F:egs5.f' !E86_019_001
            call ErrWrite(ErrID,ErrCha)
            write(*,1007) j,k
            stop
          end if
        end if
      end do
      do j = 1, nreg
        if( (k1Hscl(j) + k1Lscl(j)) .gt. 0.d0) then
          k1Hscl(j) = 0.d0
          k1Lscl(j) = 0.d0
          write(506,1008) j
        endif
      end do
      ! instead of checking charD, we can check k1max and k1min
      if(GSFile .and. .not.doGS) then
        read(517,'(72a1)',end=1000) buffer
        read(517,'(72a1)',end=1000) buffer
        do i = 1, NK1
          read(517,*) idummy, k1grd(1,i), k1grd(2,i)
        end do
        if(abs(k1mine-k1grd(1,1))/k1mine .gt. 1.d-6 .or.&
            abs(k1minp-k1grd(2,1))/k1minp .gt. 1.d-6 .or.&
              abs(k1maxe-k1grd(1,NK1))/k1maxe .gt. 1.d-6 .or.&
                abs(k1maxp-k1grd(2,NK1))/k1maxp .gt. 1.d-6) then
          write(506,1009)
          doGS = .true.
        else
          go to 30
        end if
      end if
      ! if we need to run elastino, we'll be writing the gsdist.dat
      ! file, so close whichever file we're reading on 17, open
      ! gsdist.dat, and run elastino. Then close the file and open it
      ! again for reading.  Finally, skip to the right spot
      if(doGS) then
        close(517)

        call elastino
        close(517)
     
        read(517,*) nfmeds
        do j = 1, nfmeds
          read(517,'(72a1)',end=1000) buffer
          read(517,*) hasGS, charD0, efrch0, efrcl0, ue0, ae0
        end do
      end if
      read(517,'(72a1)',end=1000) buffer
      read(517,'(72a1)',end=1000) buffer
      do i = 1, NK1
        read(517,*) idummy, k1grd(1,i), k1grd(2,i)
      end do
 30   dk1log(1) = dlog(k1grd(1,2)/k1grd(1,1))
      dk1log(2) = dlog(k1grd(2,2)/k1grd(2,1))
      imxmed = nmed !<- EGS
      allocate(lok( imxmed ) ,stat=flg_local)
      if( flg_local .ne. 0 ) then
          ErrCha = ''
          ErrID = 'L:5040/R:rmsfit/F:egs5.f' !E86_020_001
          call ErrWrite(ErrID,ErrCha)
          write(*,*)'# RMSFIT local array can not allocate !!'
          STOP '# RMSFIT local array can not allocate !!'
      endif
      nm = 0
      do k=1,nmed
        lok(k) = 0
      end do
5     continue
      read(517,'(72a1)',end=1000) buffer
      !--> look for the MEDIUM= string
      do ib=1,lmdl
        if (buffer(ib) .ne. mdlabl(ib)) go to 5
      end do
      !--> Match medium with the problem media
      do 6 k=1,nmed
        do ib=1,lmdn
          il = lmdl + ib
          if (buffer(il) .ne. media(ib,k)) go to 6
          if (ib .eq. lmdn) go to 7
        end do
 6    continue
      go to 5
 7    continue
      if (lok(k) .ne. 0) go to 5
      lok(k) = 1
      nm = nm + 1
      read(517,'(72a1)') buffer
      read(517,'(72a1)') buffer
      read(517,*) nmsgrd(k)
      read(517,'(72a1)') buffer
      read(517,*) initde(k)
      read(517,'(72a1)') buffer
      read(517,*) nmsdec(k)
      read(517,'(72a1)') buffer
      read(517,*) jskip(k)
      read(517,'(72a1)') buffer
      read(517,*) neqp(k)
      read(517,'(72a1)') buffer
      read(517,*) neqa(k)
      jskip(k) = jskip(k) - 1
      do iprt = 1,2
        read(517,'(72a1)') buffer
        do i = 1,nmsgrd(k)
          read(517,*) msgrid(i,k)
          do ik1 = 1,NK1
            read(517,'(72a1)') buffer
            read(517,*) pnoscat(iprt,i,ik1,k)
            read(517,'(72a1)') buffer
            eamu(1,iprt,i,ik1,k) = 0.0
            do j = 1, neqp(k) + neqa(k) - 1

              !--> do this now to save later
            end do
            read(517,'(72a1)') buffer
            do j = 1, neqa(k) + 1
              read(517,*) ecdf(j,iprt,i,ik1,k)
            end do
            ecdf(neqa(k)+1,iprt,i,ik1,k) = 1.d0
          end do
        end do
      end do
      if (nm .ge. nmed) go to 9
      go to 5
 9    continue
      if (nmed .eq. 1) then
        write(506,5001)
      else
        write(506,5002) nmed
      end if
      close(517)
       deallocate(lok )
      return
1000  write(506,5003)
      write(506,5004)
      do k=1,nmed
        if (lok(k) .ne. 1) write(506,'(24a1)') (media(i,k),i=1,lmdn)
      end do
      write(*,5003)
      write(*,5004)
      do k=1,nmed
        if (lok(k) .ne. 1) write(*,'(24a1)') (media(i,k),i=1,lmdn)
      end do
      close(517)
      deallocate(lok )
      stop
5001  format('Read GS Mult scat params for 1 medium')
5002  format('Read GS Mult scat params for ',i2,' media')
5003  format('End-of-file on gsdist.dat')
5004  format('The following media were not located:')
      end
      subroutine uphi(ientry,lvl)
      implicit none
      include 'include/egs5_h.f'               ! Main EGS5 "header" file
      include 'include/egs5_epcont.f'    ! COMMONs required by EGS5 code
      include 'include/egs5_stack.f'
      include 'include/egs5_uphiin.f' ! Probably don't need this anymore
      include 'include/egs5_uphiot.f'
      include 'include/counters.f'       ! Additional (non-EGS5) COMMONs
      include 'err.inc'
      integer ientry,lvl,iarg                                ! Arguments
      save usav,vsav,wsav,sinpsi,sinps2,cosdel,sindel
      save us,vs,cthet,cphi,phi,rnnow
      iuphi = iuphi + 1                    ! Count entry into subroutine
      iarg=21
      if (iausfl(iarg+1) .ne. 0) then
        call ausgab(iarg)
      end if
      go to (1,2,3),ientry
      go to 10
1     continue
      sinthe = sin(theta)
      cthet = PI5D2 - theta
      costhe = sin(cthet)
2     call randomset(rnnow)
      phi = rnnow*TWOPI
      sinphi = sin(phi)
      cphi = PI5D2 - phi
      cosphi = sin(cphi)
3     go to (4,5,6),lvl
      go to 10
4     usav = u(np)
      vsav = v(np)
      wsav = w(np)
      go to 7
6     usav = u(np-1)
      vsav = v(np-1)
      wsav = w(np-1)
5     x(np) = x(np-1)
      y(np) = y(np-1)
      z(np) = z(np-1)
      ir(np) = ir(np-1)
      wt(np) = wt(np-1)
      dnear(np) = dnear(np-1)
      latch(np) = latch(np-1)
      time(np) = time(np-1)
      k1step(np) = 0.d0
      k1init(np) = 0.d0
      k1rsd(np) = 0.d0
7     sinps2 = usav*usav + vsav*vsav
      if (sinps2 .lt. 1.e-20) then
        u(np) = sinthe*cosphi
        v(np) = sinthe*sinphi
        w(np) = wsav*costhe
      else
        sinpsi = sqrt(sinps2)
        us = sinthe*cosphi
        vs = sinthe*sinphi
        sindel = vsav/sinpsi
        cosdel = usav/sinpsi
        u(np) = wsav*cosdel*us - sindel*vs + usav*costhe
        v(np) = wsav*sindel*us + cosdel*vs + vsav*costhe
        w(np) = -sinpsi*us + wsav*costhe
      end if
      iarg=22
      if (iausfl(iarg+1) .ne. 0) then
        call ausgab(iarg)
      end if
      return
10    write(506,100) ientry,lvl
      ErrCha = ''
      ErrID = 'L:5265/R:uphi/F:egs5.f' !E86_021_001
      call ErrWrite(ErrID,ErrCha)
      write(*,100) ientry,lvl
100   format(' Stopped in uphi with ientry,lvl=',2I6)
      stop
      end
      subroutine randomset(rndum)
      implicit none
      include 'include/randomm.f'
      real*8 rndum
      integer isk
      real uni
      common /paraj/ mstz(300), parz(300)
      real*8  parz
      integer mstz

      integer iegsrand
      common /egsrand/ iegsrand
      real*8 pdummy, unirn
      if( iegsrand .eq. 1 ) then
         rndum = unirn(pdummy)
         return
      end if

      uni = seeds(j24) - seeds(i24) - carry
      if (uni .lt. 0.)  then
        uni = uni + 1.0
        carry = twom24
      else
        carry = 0.
      endif
      seeds(i24) = uni
      i24 = next(i24)
      j24 = next(j24)
      rndum = uni
      if (uni .lt. twom12)  then
        rndum = rndum + twom24*seeds(j24)
        if (rndum .eq. 0.)  rndum = twom48
      endif
      if(mstz(85).eq.99) write(93,'("(random) ",f10.5)')rndum ! debug output
      in24 = in24 + 1
      if (in24 .eq. 24)  then
        in24 = 0
        kount = kount + nskip
        do isk = 1, nskip
          uni = seeds(j24) - seeds(i24) - carry
          if (uni .lt. 0.)  then
            uni = uni + 1.0
            carry = twom24
          else
            carry = 0.
          endif
          seeds(i24) = uni
          i24 = next(i24)
          j24 = next(j24)
        end do
      endif
      ! store for restart from number of rngs used
      kount = kount + 1
      if(kount .ge. igiga) then
        mkount = mkount + 1
        kount = kount - igiga
      endif
      return
      end
      subroutine rluxinit
      implicit none
      !  inputs, from common block RLUXCOM in randomm.f
      !  integer luxlev  !  desired luxury level, default is 1,
                         !  levels 2 and higher should never be
                         !  needed
      !  integer inseed  !  seed for generating unique sequences,
                         !  in theory, any number between 1 and 2^31
                         !  "every possible value of INS gives rise
                         !  to a valid, independent sequence which will
                         !  not overlap any sequence initialized with
                         !  any other value of INS.".
      !  integer kount   !  number of randoms already used, for
      !  integer mkount  !  restart with rluxgo
      !  integer isdext  !  array of 25 ints for restart with rluxin
      include 'include/randomm.f'
      integer i, tisdext
      !  see if isdext has been set
      tisdext = 0
      do i = 1,25
        tisdext = tisdext + isdext(i)
      end do
      if(tisdext.ne.0) then
        call rluxin                    ! restart with 25 seeds
      else
        call rluxgo
      endif
      return
      end
      subroutine rluxgo
      implicit none
      include 'include/randomm.f'
      include 'err.inc'
      integer jsdflt, itwo24, icons
      parameter (jsdflt=314159265)
      parameter (itwo24=2**24, icons=2147483563)
      real uni
      integer jseed, i,k, iseeds(24), izip, izip2, isk, inner, iouter
      if (luxlev .le. 0 .or. luxlev .gt. maxlev) then
        write (506,'(a,i7)') ' illegal ranlux level in rluxgo: ',luxlev
        ErrID = 'L:5461/R:rluxgo/F:egs5.f' !E86_022_001
        call ErrWrite(ErrID,ErrCha)
        write (*,'(a,i7)') ' illegal ranlux level in rluxgo: ',luxlev
        stop
      endif
      nskip = ndskip(luxlev)
     
      in24 = 0
      if (inseed .lt. 0) then
       
        ErrID = 'L:5476/R:rluxgo/F:egs5.f' !E86_023_001
        call ErrWrite(ErrID,ErrCha)
       
        stop
      else if (inseed .gt. 0) then
        jseed = inseed
       
      else
        jseed = jsdflt
       write(506,'(a)')' ranlux initialized by rluxgo from default seed'
      endif
      inseed = jseed
      twom24 = 1.
      do i= 1, 24
        twom24 = twom24 * 0.5
        k = jseed/53668
        jseed = 40014*(jseed-k*53668) -k*12211
        if (jseed .lt. 0)  jseed = jseed+icons
        iseeds(i) = mod(jseed,itwo24)
      end do
      twom12 = twom24 * 4096.
      do i= 1,24
        seeds(i) = real(iseeds(i))*twom24
        next(i) = i-1
      end do
      next(1) = 24
      i24 = 24
      j24 = 10
      carry = 0.
      if (seeds(24) .eq. 0.) carry = twom24
      if (kount+mkount .ne. 0)  then
       
        do iouter= 1, mkount+1
          inner = igiga
          if (iouter .eq. mkount+1)  inner = kount
          do isk= 1, inner
            uni = seeds(j24) - seeds(i24) - carry
            if (uni .lt. 0.)  then
               uni = uni + 1.0
               carry = twom24
            else
               carry = 0.
            endif
            seeds(i24) = uni
            i24 = next(i24)
            j24 = next(j24)
          end do
        end do
        in24 = mod(kount, nskip+24)
        if (mkount .gt. 0)  then
           izip = mod(igiga, nskip+24)
           izip2 = mkount*izip + in24
           in24 = mod(izip2, nskip+24)
        endif
        if (in24 .gt. 23) then
          
           ErrID = 'L:5547/R:rluxgo/F:egs5.f' !E86_024_001
           call ErrWrite(ErrID,ErrCha)
          
           stop
        endif
      endif
      rluxset = .true.
      return
      end
      subroutine rluxin
      implicit none
      include 'include/randomm.f'
      include 'err.inc'
      integer i, isd
      write(506,'(a)')' full initialization of ranlux with 25 integers:'
      write(506,'(5x,5i12)') isdext
      twom24 = 1.
      do i = 1, 24
        next(i) = i-1
        twom24 = twom24 * 0.5
      end do
      next(1) = 24
      twom12 = twom24 * 4096.
      do i = 1, 24
        seeds(i) = real(isdext(i))*twom24
      end do
      carry = 0.
      if (isdext(25) .lt. 0)  carry = twom24
      isd = iabs(isdext(25))
      i24 = mod(isd,100)
      isd = isd/100
      j24 = mod(isd,100)
      isd = isd/100
      in24 = mod(isd,100)
      isd = isd/100
      luxlev = isd
      if (luxlev .le. maxlev) then
        nskip = ndskip(luxlev)
      
      else
        write (506,'(a,i5)') ' ranlux illegal luxury rluxin: ',luxlev
        ErrID = 'L:5608/R:rluxin/F:egs5.f' !E86_025_001
        call ErrWrite(ErrID,ErrCha)
        write (*,'(a,i5)') ' ranlux illegal luxury rluxin: ',luxlev
        stop
      endif
      inseed = -1
      rluxset = .true.
      return
      end
      subroutine rluxout
      implicit none
      include 'include/randomm.f'
      real twop12
      parameter (twop12=4096.)
      integer i
      do i = 1,24
        isdext(i) = int(seeds(i)*twop12*twop12)
      end do
      isdext(25) = i24 + 100*j24 + 10000*in24 + 1000000*luxlev
      if (carry .gt. 0.)  isdext(25) = -isdext(25)
      return
      end
      double precision function csdar(iq,e)
      implicit none
      include 'include/egs5_h.f'
      include 'pegscommons/rngspl.f'
      include 'pegscommons/dercon.f'
      integer iq
      double precision e
      integer i
      double precision ke
      ke = e - RM
      call findi(erng,ke,nrng,i)
      if(iq .eq. -1) then
       csdar = arnge(i)+ke*(brnge(i)+ke*(crnge(i)+ke*drnge(i)))
      else
       csdar = arngp(i)+ke*(brngp(i)+ke*(crngp(i)+ke*drngp(i)))
      endif
      return
      end
      double precision function DCSEL(RMU)
      implicit none
      include 'include/egs5_h.f'
      include 'include/egs5_cdcsep.f'
      include 'include/egs5_cdcspl.f'
      double precision rmu
      integer i
      CALL FINDI(XMU,RMU,NREDA,I)
      DCSEL=EXP(RA(I)+RMU*(RB(I)+RMU*(RC(I)+RMU*RD(I))))
      RETURN
      END
      double precision function DCSN(RMU)
      implicit none
      include 'include/egs5_csplcf.f'
      double precision  rmu, x, pl(1000)
      X=1.0D0-2.0D0*RMU
      CALL LEGENP(X,PL,IL)
      DCSN=EXP(ASPL+RMU*(BSPL+RMU*(CSPL+RMU*DSPL)))*PL(IL)
      RETURN
      END
      subroutine dcsstor(ecut,emax)
      USE PEGS_DCSSTR_MOD !<- allocatable 2015.08xx
      implicit none
      real*8 ecut, emax
      include 'include/egs5_h.f'
      include 'include/egs5_cdcsep.f'
      include 'pegscommons/molvar.f'
      include 'pegscommons/pmcons.f'
      include 'pegscommons/mxdatc.f'
      include 'pegscommons/mscom.f'
      integer i,j
      nsdcs = nsdcs + 1
      atomd(nsdcs) = an * rho / wm
      efrch(nsdcs) = efrach
      efrcl(nsdcs) = efracl
      do i = 1, 24
        mednam(nsdcs,i) = medium(i)
      end do
      do i = 1, negrid
        secs(nsdcs,i) = ecs(i)
        setcs1(nsdcs,i) = etcs1(i)
        setcs2(nsdcs,i) = etcs2(i)
        spcs(nsdcs,i) = pcs(i)
        sptcs1(nsdcs,i) = ptcs1(i)
        sptcs2(nsdcs,i) = ptcs2(i)
        do j= 1, nreda
          sedcs(nsdcs,i,j) = edcs(i,j)
          spdcs(nsdcs,i,j) = pdcs(i,j)
        end do
      end do
      egrdlo(nsdcs) = ecut
      egrdhi(nsdcs) = emax
      nlegmd(nsdcs) = nleg0
      return
      end
      subroutine dcsload(n,ecut,emax,nleg0)
      USE PEGS_DCSSTR_MOD !<- allocatable 2015.08xx
      include 'include/egs5_h.f'
      include 'include/egs5_cdcsep.f'
      include 'pegscommons/mxdatc.f'
      integer n, nleg0
      double precision ecut, emax
      integer i,j
      do i = 1, 24
        medium(i) = mednam(n,i)
      end do
      do i = 1, negrid
        ecs(i) = secs(n,i)
        etcs1(i) = setcs1(n,i)
        etcs2(i) = setcs2(n,i)
        pcs(i) = spcs(n,i)
        ptcs1(i) = sptcs1(n,i)
        ptcs2(i) = sptcs2(n,i)
        do j= 1, nreda
          edcs(i,j) = sedcs(n,i,j)
          pdcs(i,j) = spdcs(n,i,j)
        end do
      end do
      ecut = egrdlo(n)
      emax = egrdhi(n)
      nleg0 = nlegmd(n)
      return
      end
      SUBROUTINE DCSTAB(E,IELEC)
      include 'include/egs5_h.f'
      include 'include/egs5_cdcsep.f'
      include 'include/egs5_cdcspl.f'
      include 'include/egs5_coefgs.f'
      integer ielec
      double precision e
      integer ia, ie
      double precision el

      DO IE=1,NEGRID
        X(IE)=LOG(ET(IE))
      ENDDO
      EL=LOG(E)
      DO IA=1,NREDA
        DO IE=1,NEGRID
          IF(IELEC.EQ.-1) THEN
            Y(IE)=LOG(EDCS(IE,IA))
          ELSE
            Y(IE)=LOG(PDCS(IE,IA))
          ENDIF
        ENDDO
        CALL egs5SPLINE(X,Y,A,B,C,D,0.0D0,0.0D0,NEGRID)
        CALL FINDI(X,EL,NEGRID,IE)
        DCSIL(IA)=A(IE)+EL*(B(IE)+EL*(C(IE)+EL*D(IE)))
        DCSI(IA)=EXP(DCSIL(IA))
      ENDDO
      CALL egs5SPLINE(XMU,DCSIL,RA,RB,RC,RD,0.0D0,0.0D0,NREDA)
      RETURN
      END
      subroutine elastino
      USE PEGS_DCSSTR_MOD
      USE egs5_media_mod
      USE egs5_mscon_mod
      implicit none
      include 'include/egs5_h.f'
      include 'include/egs5_cdcsep.f'
      include 'include/egs5_coefgs.f'
      include 'include/egs5_useful.f'
      double precision Y(NREDA),RA(NREDA),RB(NREDA),RC(NREDA),RD(NREDA)
      double precision Y1(NREDA),Y2(NREDA)
      integer ntb
      PARAMETER (NTB=NFIT1)

      logical printDiag, longDiag, otherM, uniInt
      integer i, j, n, ik1, iener, ipart, ielec, nleg, ntermm, nterm
      integer ia, didGS
      integer nleg0
      double precision e, s, thdeg,raddeg
      double precision sumt, sum, sump, dsum, xl, xu, xx, st, dmu
      double precision k1start(2), dk1(2), k1, ecut, emax
      double precision DCSEL
      common /paran/ icfn(100), ilfn(100), chfn(100)
      character chfn*200
      integer icfn, ilfn
      raddeg = 180/3.14159
      longDiag = .false.
      printDiag = .true.
      otherM = .false.
      uniInt = .false.
      didGs = 1
      write(517,*) nsdcs
      do n = 1, nsdcs
        write(517,5001) (mednam(n,i),i=1,24)
      end do
5001  format(' MEDIUM=',24A1)
      
      dk1(1) = dexp( (1.d0/(NK1-1.d0)) * dlog(k1maxe/k1mine) )
      dk1(2) = dexp( (1.d0/(NK1-1.d0)) * dlog(k1maxp/k1minp) )
      k1start(1) = k1mine
      k1start(2) = k1minp
      write(517,*) 'Scattering strength ladder for this run'
      write(517,*) 'IK1       K1(e-)         K1(e+)'
      do ik1 = 1, NK1
        
      end do
      if(printDiag) write(526,1001)
1001  format('In Elastino to get Multiple Scattering distribution')
      do n = 1, nsdcs
      call dcsload(n,ecut,emax,nleg0)
      call inigrd(ecut,emax,0)
      do ipart = 1,2
        NLEG = MIN(nleg0,NGT)
        IF(ipart.EQ.1) THEN
          ielec = -1
          if(printDiag) WRITE(526,*) 'Projectile: electron.'
        ELSE
          ielec = +1
          if(printDiag) WRITE(526,*) 'Projectile: positron.'
        ENDIF
        do iener = 1,nmscate
          e = mscate(iener)
          if(e .gt. 1.d8) go to 1005
          CALL DCSTAB(E,IELEC)
          CALL GSCOEF(NLEG)
          if(printDiag) then
          WRITE(526,*) 'Kinetic energy (MeV) =',E*1.d-6
          WRITE(526,*) ' Total cross section (cm^-1) =',CS*atomd(n)
          WRITE(526,*) ' 1st transport cross section =',TCS1*atomd(n)
          if(longDiag) WRITE(526,*)' 2nd transport cross section =',TCS2
          endif
          do ik1 = 1, NK1
            k1 = k1start(ipart) * dk1(ipart) ** (ik1-1)
            s = k1 * cs / tcs1
            if(printDiag) then
            WRITE(526,*) ' For this step, K1 =',k1
            WRITE(526,*) '  At this K1, path length (mfps)=',S
           
            endif
            if(longDiag) then
              WRITE(526,*) '# '
              WRITE(526,*) '# Differential cross section'
              WRITE(526,*) '# theta(deg)       mu        dcs(cm^2/sr)'
              DO IA=1,NREDA
                WRITE(526,'(1P,4E14.6)') TH(IA),XMU(IA),DCSI(IA)
              ENDDO
              WRITE(526,*) '# Goudsmit-Saunderson distribution, PDF(mu)'
              IF(IELEC.EQ.1) THEN
                WRITE(526,*) '# Projectile: positron'
              ELSE
                WRITE(526,*) '# Projectile: electron'
              ENDIF
              WRITE(526,*) '# Number of terms in Legendre series =',NLEG
              WRITE(526,*) '# '

            endif
            nterm = nleg
            CALL GSDIST(S,0.0D0,Y(1),NTERM)
            NTERMM=NTERM
            DO I=1,NREDA
              CALL GSDIST(S,XMU(I),Y(I),NTERM)
              if(longDiag) then
                THDEG=ACOS(1.0D0-2.0D0*XMU(I))*RADDEG

              endif
              NTERMM=MAX(NTERMM,NTERM)
            ENDDO
    2       CONTINUE
            if(printDiag) WRITE(526,*) ' Number of terms used: ',NTERMM
            !  this is a temp fix until single scattering
            if(nterm.lt.0) nterm = ntermm
            CALL egs5SPLINE(XMU,Y,RA,RB,RC,RD,0.0D0,0.0D0,NREDA)
            CALL egs5INTEG(XMU,RA,RB,RC,RD,0.0D0,XMU(NREDA),SUMP,NREDA)
            if(s.lt.100)  then
              SUM=SUMP+EXP(-S)
            else
              sum = sump
            endif
            if(printDiag) then
            WRITE(526,*) ' Normalization =',SUM
            WRITE(526,*) ' No-scatter probability = ',sum-sump
            endif
            probns(ipart,iener,ik1) = sum-sump
            DSUM=SUMP/DBLE(NBFIT)
            XGR(1)=0.0D0
            DO I=2,NBFIT
              SUMT=(I-1)*DSUM
              XL=XGR(I-1)
              XU=1.0D0
1234          XX=0.5D0*(XL+XU)
              CALL egs5INTEG(XMU,RA,RB,RC,RD,0.0D0,XX,ST,NREDA)
              IF(ST.GT.SUMT) THEN
                XU=XX
              ELSE
                XL=XX
              ENDIF
              IF(ABS(XU-XL).GT.1.0D-6*XX) GO TO 1234
              XGR(I)=XX
            ENDDO
            dmu = (1.0d0 - xgr(NBFIT)) / NEXFIT
            do i=1,NEXFIT
              j = NBFIT + i
              xgr(j) = xgr(NBFIT) + i * dmu
            end do
            DO I=1,NREDA
              YY(I)=Y(I)/SUM
            ENDDO
            CALL egs5SPLINE(XMU,YY,TA,TB,TC,TD,0.0D0,0.0D0,NREDA)
            nterm = nleg
            DO I=1,NTB
              CALL GSDIST(S,XGR(I),YC(I),NTERM)
              if(i.le.NBFIT) then
                y0(i) = (i-1)*dsum/sum
              else
                CALL egs5INTEG(XMU,TA,TB,TC,TD,0.0D0,XGR(I),Y0(I),NREDA)
              endif
            ENDDO
            call fitms(xgr,yc,y0,ntb,ipart,iener,ik1)
           ! other moments...
            if(otherM) then
              DO I=1,NREDA
                YY(I)=YY(I)*XMU(I)
              ENDDO
              CALL egs5SPLINE(XMU,YY,TA,TB,TC,TD,0.0D0,0.0D0,NREDA)
              DO I=1,NTB
                CALL egs5INTEG(XMU,TA,TB,TC,TD,0.0D0,XGR(I),Y1(I),NREDA)
              ENDDO
              DO I=1,NREDA
                YY(I)=YY(I)*XMU(I)
              ENDDO
              CALL egs5SPLINE(XMU,YY,TA,TB,TC,TD,0.0D0,0.0D0,NREDA)
              DO I=1,NTB
                CALL egs5INTEG(XMU,TA,TB,TC,TD,0.0D0,XGR(I),Y2(I),NREDA)
              ENDDO
              if(printDiag) then


              IF(IELEC.EQ.1) THEN
                WRITE(526,*) '# Projectile: positron'
              ELSE
                WRITE(526,*) '# Projectile: electron'
              ENDIF
              WRITE(526,*) '# Kinetic energy (eV) =',E
              WRITE(526,*) '# Path length/mfp =',S
              WRITE(526,*) '# Number of terms in Legendre series =',NLEG
              if(longDiag) then
              WRITE(526,*) '# '
              WRITE(526,3001)

              DO I=1,NTB
                THDEG=ACOS(1.0D0-2.0D0*XGR(I))*RADDEG

              ENDDO
              endif
              endif
            endif
            if(uniInt) then
              DMU=1.0D0/DBLE(NTB-1)
              DO I=1,NTB
                XGR(I)=(I-1)*DMU
              ENDDO
              DO I=1,NREDA
                YY(I)=Y(I)/SUM
              ENDDO
              CALL egs5SPLINE(XMU,YY,TA,TB,TC,TD,0.0D0,0.0D0,NREDA)
              DO I=1,NTB
                CALL GSDIST(S,XGR(I),YC(I),NTERM)
                CALL egs5INTEG(XMU,TA,TB,TC,TD,0.0D0,XGR(I),Y0(I),NREDA)
              ENDDO
              call fitms(xgr,yc,y0,ntb,ipart,iener,ik1)
              if(otherM) then
                DO I=1,NREDA
                  YY(I)=YY(I)*XMU(I)
                ENDDO
                CALL egs5SPLINE(XMU,YY,TA,TB,TC,TD,0.0D0,0.0D0,NREDA)
                DO I=1,NTB
                CALL egs5INTEG(XMU,TA,TB,TC,TD,0.0D0,XGR(I),Y1(I),NREDA)
                ENDDO
                DO I=1,NREDA
                  YY(I)=YY(I)*XMU(I)
                ENDDO
                CALL egs5SPLINE(XMU,YY,TA,TB,TC,TD,0.0D0,0.0D0,NREDA)
                DO I=1,NTB
                CALL egs5INTEG(XMU,TA,TB,TC,TD,0.0D0,XGR(I),Y2(I),NREDA)
                ENDDO
                if(printDiag) then

                IF(IELEC.EQ.1) THEN
                  WRITE(526,*) '# Projectile: positron'
                ELSE
                  WRITE(526,*) '# Projectile: electron'
                ENDIF
                WRITE(526,*)'# Kinetic energy (eV) =',E
                WRITE(526,*)'# Path length/mfp =',S
                WRITE(526,*)'# Number of terms in Legendre series=',NLEG
                if(longDiag) then
                WRITE(526,*) '# '
                WRITE(526,3002)

                DO I=1,NTB
                  THDEG=ACOS(1.0D0-2.0D0*XGR(I))*RADDEG
   
                ENDDO
                endif
                endif
              endif        !  get other moments
            endif          !  use uniform ints
          end do	   !  end path length loop
        end do	           !  end energy loop
      end do               !  end particle types loop
1005  call wmsfit(k1start,dk1)
1010  end do               !  end material loop
      close(526)
      return
      END
      SUBROUTINE ELINIT(Z,STF,NELEM)
      use EGS5_MS_MOD !FURUTA20140825
      implicit none
      include 'include/egs5_h.f'
      include 'include/egs5_cdcsep.f'
      include 'include/egs5_uphiot.f'
      include 'pegscommons/scpspl.f'
      include 'pegscommons/mscom.f'
      include 'pegscommons/dercon.f'
      common /paran/ icfn(100), ilfn(100), chfn(100)
      character chfn*200
      integer icfn, ilfn
      integer nelem, IZ(MXEPERMED)
      double precision STF(NELEM), Z(NELEM), zcor(MXEPERMED)
      integer i,ie,ia, iel, ielec, izz, izr, ns, ns1,ns2,ns3
      double precision zplus1, stff, enr, csin, csin1, csin2
      CHARACTER*1 LIT10(10),LIT1,LIT2,LIT3
      DATA LIT10/'0','1','2','3','4','5','6','7','8','9'/
      CHARACTER*24 FILE1,FILE2
      double precision EGRID(16)

      integer igrid
      double precision fgrid, e
      double precision g1e
      logical printDiag
      printDiag = .false.
      if(et(1).ne.1.d2) then
      IE=0
      IGRID=0
      FGRID=100.0D0
   10 IGRID=IGRID+1
      E=EGRID(IGRID)*FGRID
      IF(IGRID.EQ.16) THEN
        IGRID=1
        FGRID=10.0D0*FGRID
      ENDIF
      IE=IE+1
      ET(IE)=E
      IF(IE.LT.NEGRID) GO TO 10
      endif
      do i=1,nelem
        iz(i) = z(i)
        zcor(i) =  (z(i) + fudgeMS) / z(i)
      end do
      DO IE=1,NEGRID
        ECS(IE)=0.0D0
        ETCS1(IE)=0.0D0
        ETCS2(IE)=0.0D0
        PCS(IE)=0.0D0
        PTCS1(IE)=0.0D0
        PTCS2(IE)=0.0D0
        DO IA=1,NREDA
          EDCS(IE,IA)=0.0D0
          PDCS(IE,IA)=0.0D0
        ENDDO
      ENDDO
      DO IEL=1,NELEM
        IZZ=IZ(IEL)
        STFF=STF(IEL)
        zplus1 = zcor(iel)
        NS=IZ(IEL)
        IF(NS.GT.999) NS=999
        NS1=NS-10*(NS/10)
        NS=(NS-NS1)/10
        NS2=NS-10*(NS/10)
        NS=(NS-NS2)/10
        NS3=NS-10*(NS/10)
        LIT1=LIT10(NS1+1)
        LIT2=LIT10(NS2+1)
        LIT3=LIT10(NS3+1)
       
        
        DO IE=1,NEGRID
          
          ECS(IE) = ECS(IE) + zplus1 * stff * csin
          ETCS1(IE) = ETCS1(IE) + zplus1 * stff * csin1
          ETCS2(IE) = ETCS2(IE) + zplus1 * stff * csin2
         
          READ(531,'(1P,10E12.5)') (DCSI(IA),IA=1,NREDA)
          DO IA=1,NREDA
            EDCS(IE,IA)=EDCS(IE,IA)+zplus1*STFF*DCSI(IA)
          ENDDO
         
          PCS(IE) = PCS(IE) + zplus1 * stff * csin
          PTCS1(IE) = PTCS1(IE) + zplus1 * stff * csin1
          PTCS2(IE) = PTCS2(IE) + zplus1 * stff * csin2
         
          READ(532,'(1P,10E12.5)') (DCSI(IA),IA=1,NREDA)
          DO IA=1,NREDA
            PDCS(IE,IA)=PDCS(IE,IA)+zplus1*STFF*DCSI(IA)
          ENDDO
        ENDDO
        CLOSE(531)
        CLOSE(532)
      ENDDO
      do ie = 1,negrid-1
        etl(ie) = dlog(et(ie))
        etcs1(ie) = dlog(etcs1(ie))
        ptcs1(ie) = dlog(ptcs1(ie))
      end do
      etcs1(negrid) = dlog(g1e(-1,110.d0+RM,-1))
      ptcs1(negrid) = dlog(g1e(+1,110.d0+RM,-1))
      etcs1(negrid+1) = dlog(g1e(-1,120.d0+RM,-1))
      ptcs1(negrid+1) = dlog(g1e(+1,120.d0+RM,-1))
      etcs1(negrid+2) = dlog(g1e(-1,130.d0+RM,-1))
      ptcs1(negrid+2) = dlog(g1e(+1,130.d0+RM,-1))
      etl(negrid) = dlog(110.d6)
      etl(negrid+1) = dlog(120.d6)
      etl(negrid+2) = dlog(130.d6)
      call egs5spline(etl,etcs1,ag1e,bg1e,cg1e,dg1e,0.d0,0.d0,negrds)
      call egs5spline(etl,ptcs1,ag1p,bg1p,cg1p,dg1p,0.d0,0.d0,negrds)
      RETURN
      END
      subroutine esteplim
      USE egs5_mscon_mod
      implicit none
      real*8 g1e, sptote, sptotp, csdar
      include 'include/egs5_h.f'
      include 'pegscommons/dercon.f'
      include 'pegscommons/molvar.f'
      include 'pegscommons/rngspl.f'
      include 'pegscommons/thres2.f'
      include 'err.inc'
      integer i, ivar
      data etol/1.d-3/
      !  load the energy array
      nrng = NESCPW
      erng(1) = log(mscate(1)*1.d-6)
      erng(nrng) = log(mscate(nmscpw)*1.d-6)
      dloge = (erng(nrng) - erng(1)) / (nrng - 1)
      do i = 2, nrng-1
        erng(i) = erng(1) + (i-1) * dloge
      end do
      do ivar = 1,4
        do i=1,nrng
          if (ivar.eq.1) then
            erng(i) = exp(erng(i))
            estepl(i) = .50
          end if
          e1t = erng(i) + RM
          if (ivar.eq.1) then
            fspe(i) = rlc / sptote(e1t,ae,ap)
          else if (ivar.eq.2) then
            fspe(i) = rlc / sptotp(e1t,ae,ap)
          else if (ivar.eq.4) then
            fspe(i) = g1e(-1,e1t,0) * (rlc / sptote(e1t,ae,ap))
          else if (ivar.eq.5) then
            fspe(i) = 1.d0
          end if
        end do
        !  spline the function for local use
        if (ivar.ne.3) then
          call egs5spline(erng,fspe,aspe,bspe,cspe,dspe,0.d0,0.d0,nrng)
        end if
        !  get the electron or positron CSDA range
        if (ivar.le.2) then
          e1 = erng(1)
          if (ivar.eq.1) then
            rnge(1) = rlc / sptote(e1t,ae,ap) * e1
          else
            rnge(1) = rlc / sptotp(e1t,ae,ap) * e1
          end if
          do i = 2, nrng
            e2 = erng(i)
            call egs5integ(erng,aspe,bspe,cspe,dspe,e1,e2,anal,nrng)
            rnge(i) = rnge(i-1) + dabs(anal)
            e1 = e2
          end do
          !  get global splines for use later is csdar
          if (ivar.eq.1) then
            call egs5spline(erng,rnge,arnge,brnge,crnge,drnge,&
                                              0.d0,0.d0,nrng)
          else
            call egs5spline(erng,rnge,arngp,brngp,crngp,drngp,&
                                              0.d0,0.d0,nrng)
          end if
        !  do the energy limit
        else
          !  start at 4 because of spline issues at ends
          do i = 4, nrng-3
            e1 = erng(i)
            ! Expand energy region to avoid unstability for second material.
            ! To avoid negative estepe.
            ! Related to the following part in the "subroutine egs5efpl":
            ! --> if(estepe.lt.1.0d-10) then
            if (ivar.eq.4 .and. e1.gt.30.d0 .and. e1.lt.300.d0) then
              if (i.gt.4) then
                estepl(i) = estepl(i-1)
                go to 200
              end if
            end if
            e2 = e1 * (1.d0 - estepl(i))
            if (e2.lt.erng(1)) then
              e2 = erng(1)
              estepl(i) = 1.d0 - e2/e1
            endif
            e1t = e1 + RM
 100        e2t = e2 + RM
            if (ivar.le.4) then
              deltaE = (e1 - e2)
              if (ivar.eq.3) then
                trap = rlc / sptote(e1t,ae,ap) + rlc / sptote(e2t,ae,ap)
                anal = csdar(-1,e1t) - csdar(-1,e2t)
              else
                trap = rlc * g1e(-1,e1t,0) / sptote(e1t,ae,ap)
                trap = trap + rlc * g1e(-1,e2t,0) / sptote(e2t,ae,ap)
                call egs5integ(erng,aspe,bspe,cspe,dspe,e2,e1,anal,nrng)
              end if
              trap = deltaE * trap / 2.d0
              !  this should never happen now except in development
              if (anal.eq.0.d0) then
                write(506,*) 'ERROR in esteplim -- integral = 0.d0'
                ErrCha = ''
                ErrID = 'L:6618/R:esteplim/F:egs5.f' !E86_026_001
                call ErrWrite(ErrID,ErrCha)
                write(*,*) 'ERROR in esteplim -- integral = 0.d0'
                stop
              ! if we're within the tolerance, get next energy
              else if (dabs((trap - anal) / anal) .le. etol) then
                ! convergence -- ignore the discreteness from using 5%
                go to 200
              else
              ! else, store previous, and try something shorter
                e2 = e1 * (1.d0 - estepl(i)*.95)
                estepl(i) = 1.d0 - e2/e1
                !  trap numerical problems -- estepe is smooth
                if(i.gt.4 .and. estepl(i).lt.estepl(i-1)/2.d0) go to 200
                go to 100
              endif
            else
          !  do the mfp limit
            end if
 200        continue
          end do
          !  skip regions where splines may be bad
          do i = 1, 3
            estepl(i) = estepl(4)
            estepl(nrng-3+i) = estepl(nrng-3)
          end do
        end if
      end do
      ! get splines for use in pegs
      call egs5spline(erng,estepl,aeste,beste,ceste,deste,&
                     0.d0,0.d0,nrng)
      return
      end
      double precision function estepmax(e)
      implicit none
      include 'include/egs5_h.f'
      include 'pegscommons/rngspl.f'
      include 'pegscommons/dercon.f'
      double precision e
      integer i
      double precision ke
      ke = e - RM
      call findi(erng,ke,nrng,i)
      estepmax = aeste(i)+ke*(beste(i)+ke*(ceste(i)+ke*deste(i)))
      return
      end
      SUBROUTINE FINDI(X,XC,N,I)
      implicit none
      integer n,i, i1, it
      double precision xc, X(N)
      IF(XC.GT.X(N)) THEN
        I=N-1
        RETURN
      ENDIF
      IF(XC.LT.X(1)) THEN
        I=1
        RETURN
      ENDIF
      I=1
      I1=N
    1 IT=(I+I1)/2
      IF(XC.GT.X(IT)) I=IT
      IF(XC.LE.X(IT)) I1=IT
      IF(I1-I.GT.1) GO TO 1
      RETURN
      END
      subroutine fitms(x,p,cdf,na, ipart,eindex,ik1)
      USE egs5_mscon_mod
      implicit none
      include 'include/egs5_h.f'
      include 'err.inc'
      integer ipart, eindex, na, ik1, jint
      double precision x(na), p(na), cdf(na)
      integer NB
      parameter (NB=NFIT+1)
      integer mesh, istart, iend, interval, i
      double precision fit(NB), r, a, b, c, n
      double precision sum, diff, biggest
      double precision p1, p2, x1, x2, n1, n2, ln21, alpha
      double precision F, G
      logical printDiag
      printDiag = .false.
      if(na.ne.NB) then
        write (506,*) 'Error:  number of MS dist points != NB.  Stop.'
        ErrID = 'L:6782/R:fitms/F:egs5.f' !E86_027_001
        call ErrWrite(ErrID,ErrCha)
        write (*,*) 'Error:  number of MS dist points != NB.  Stop.'
        stop
      endif
      mesh = NB-1
      sum = 0
      biggest = 0
      do interval = 1, mesh
        istart = interval
        iend = interval + 1
        p1 = p(istart)
        p2 = p(iend)
        x1 = x(istart)
        x2 = x(iend)
        if(p2 .ne. p1) then
          r = sqrt(p2/p1)
          n = (r*x2 - x1)/(1 - r)
          a = p1*(x1 + n)**2
        else
          r = 1
          n = 1e10
          a = p1*n**2
        endif
        n1 = x1 + n
        n2 = x2 + n
        ln21 = log(n2/n1)
        F = 1/n1 - 1/n2
        G = (n1 + n2)*ln21 - 2*(n2 - n1)
        alpha = cdf(iend) - cdf(istart) - a*F
        b = alpha/G
        c = 0
        if (p2 .ne. p1) then
          b = b/a
          c = c/a
        else
          b = 0
          c = 0
        endif
        if(printDiag) then
        do i = istart, iend
          
          diff = abs((fit(i) - p(i))/p(i))
          biggest = max(biggest,diff)
          sum = sum + diff
          write(526,*) x(i), (fit(i) - p(i))/p(i)
        enddo
        write(526,'(''  Biggest relative difference = '',g14.7)')biggest
        write(526,'(''  Goodness of fit = '',g14.7)') sum/NA
        endif
        ams(ipart,eindex,ik1,interval) = a
        bms(ipart,eindex,ik1,interval) = b
        cms(ipart,eindex,ik1,interval) = c
        etams(ipart,eindex,ik1,interval) = n
        amums(ipart,eindex,ik1,interval) = x1
        if(interval.ge.NBFIT) then
          jint = interval - NBFIT + 1
          
        endif
      enddo
      cumdist(ipart,eindex,ik1,NEXFIT) = 1.d0
      amums(ipart,eindex,ik1,mesh+1) = x2
      return
      end
      double precision function g1ededx(e)
      implicit none
      include 'pegscommons/thres2.f'
      include 'pegscommons/molvar.f'
      real*8 e, sptote, g1e
      g1ededx = g1e(-1,e,0) * (rlc / sptote(e,e,e))
      return
      end
      double precision function g1pdedx(e)
      implicit none
      include 'pegscommons/thres2.f'
      include 'pegscommons/molvar.f'
      real*8 e, sptotp, g1e
      g1pdedx = g1e(+1,e,0) * (rlc / sptotp(e,e,e))
      return
      end
      double precision function g1e(iq,e,iflag)
      implicit none
      include 'include/egs5_h.f'
      include 'pegscommons/scpspl.f'
      include 'pegscommons/molvar.f'
      include 'pegscommons/pmcons.f'
      include 'pegscommons/dercon.f'
      integer iq,iflag
      double precision e
      integer i
      double precision ke, xa2, onemb2, beta2, g1mol, logke
      ke = e - RM
      logke = dlog(ke*1.d6)
      if(logke.le.etl(negrds)) then
        call findi(etl,logke,negrds,i)
        if(iq .eq. -1) then
         
        else
      
        endif
      else
        onemb2 = (1.d0 + (ke/RM) )**(-2)
        beta2 = 1.d0 - onemb2
        xa2 = fsc*fsc * 1.13 /(.885*.885) * dexp(zx/zs) / dexp(ze/zs)
        xa2 = xa2 * onemb2 / beta2
        g1mol = dlog((pi*pi + xa2)/xa2) - pi*pi/(pi*pi + xa2)
        g1e = r0*r0 * zs * 2*pi * g1mol * onemb2 / (beta2*beta2)
      endif
      if(iflag.ne.-1) g1e = g1e * an * rho / wm
      return
      end
      SUBROUTINE GAULEG(X,W,N)
      implicit none
      integer N
      double precision  X(N),W(N)
      integer i,j,m
      double precision xm,xl, z,z1, p1,p2,p3,pp
      double precision eps
      PARAMETER (EPS=1.0D-15)
      M=(N+1)/2
      XM=0.0d0
      XL=1.0d0
      DO I=1,M
        Z=COS(3.141592654D0*(I-0.25D0)/(N+0.5D0))
    1   CONTINUE
          P1=1.0D0
          P2=0.0D0
          DO J=1,N
            P3=P2
            P2=P1
            P1=((2.0D0*J-1.0D0)*Z*P2-(J-1.0D0)*P3)/J
          ENDDO
          PP=N*(Z*P1-P2)/(Z*Z-1.0D0)
          Z1=Z
          Z=Z1-P1/PP
        IF(ABS(Z-Z1).GT.EPS*ABS(Z)+EPS) GO TO 1
        X(I)=XM-XL*Z
        X(N+1-I)=XM+XL*Z
        W(I)=2.0D0*XL/((1.0D0-Z*Z)*PP*PP)
        W(N+1-I)=W(I)
      ENDDO
      RETURN
      END
      SUBROUTINE GSCOEF(NLEG)
      implicit none
      include 'include/egs5_h.f'
      include 'include/egs5_uphiot.f'
      include 'include/egs5_cdcsep.f'
      include 'include/egs5_coefgs.f'
      integer nleg
      double precision DCSEL
      integer nxc, nlm, j,i,l
      double precision dcsl, dcsu, xi, funxi
      double precision X(NGT),W(NGT)
      double precision PL(NGT),XC(100),F0(100),F1(100)
      NLM=MAX(500,MIN(NLEG/2,NGT/2))
      CALL GAULEG(X,W,NLM)
      J=1
      XC(J)=0.0D0
      DCSL=1.0D-1*DCSI(1)
      DCSU=1.0D+1*DCSI(1)
      DO I=2,NREDA
        
      ENDDO
      IF (XC(J).NE.1.0D0) THEN
        J=J+1
        XC(J)=1.0D0
      ENDIF
      NXC=J
      DO J=2,NXC
        F0(J-1)=(XC(J)+XC(J-1))/2.0D0
        F1(J-1)=(XC(J)-XC(J-1))/2.0D0
      ENDDO
      NLEGEN=NLEG
      IF(NLEG.GT.NGT) NLEGEN=NGT
      DO L=1,NLEGEN
        GL(L)=0.0D0
      ENDDO
      DO I=1,NLM
        DO J=1,NXC-1
          XI=F0(J)+X(I)*F1(J)
          FUNXI=F1(J)*W(I)*DCSEL(XI)
          CALL LEGENP(1.0D0-2.0D0*XI,PL,NLEGEN)
          DO L=1,NLEGEN
            IF(L.EQ.1) THEN
              GL(L)=GL(L)+FUNXI
            ELSE
              GL(L)=GL(L)+(1.0D0-PL(L))*FUNXI
            ENDIF
          ENDDO
        ENDDO
      ENDDO
      CS0=GL(1)
      GL(1)=0.0D0
      DO L=2,NLEGEN
        GL(L)=GL(L)/CS0
      ENDDO
      CS0=2.0D0*CS0
      CS=CS0*TWOPI
      TCS1=CS*GL(2)
      TCS2=CS*GL(3)
      RETURN
      END
      SUBROUTINE GSDIST(S,RMU,PDF,NTERM)
      implicit none
      integer nterm
      double precision s, rmu, pdf, DCSEL
      include 'include/egs5_h.f'
      include 'include/egs5_coefgs.f'
      logical printDiag
      integer l
      double precision x, dxs, sum, sumf, f
      double precision  PL(NGT)
      printDiag = .false.
      IF(RMU.LT.0.0D0.OR.RMU.GT.1.000001D0) STOP 'GS error.'
      X=1.0D0-2.0D0*RMU
      if(nterm.lt.0) then
        SUM=S*DCSEL(RMU)/CS0
        go to 200
      end if
      CALL LEGENP(X,PL,NLEGEN)
      ! trap to prevent underflow exceptions for thick targets
      if(s.lt.100.d0) then
        DXS=EXP(-S)
      else
        dxs = 0.d0
      endif
      SUM=0.0D0
      SUMF=0.0D0
      NTERM=NLEGEN
      DO L=1,NLEGEN
        ! trap to prevent underflow exceptions for thick targets
        if(s*GL(L).lt.100.d0) then
          F=(L-0.5D0)*(EXP(-S*GL(L))-DXS)
        else
          F=-(L-0.5D0)*DXS
        endif
        SUM=SUM+F*PL(L)
        SUMF=SUMF+F
        IF(ABS(F).LT.1.0D-6*ABS(SUM)) then
          NTERM=MIN(NTERM,L)
          if(l.gt.1) go to 100
        endif
      ENDDO
  100 continue
      IF(ABS(F).GT.1.0D-2*ABS(SUM)) THEN
        if(printDiag) then
        WRITE(526,1000) RMU,SUM,F

        endif
        if(sum.gt.0d0) SUM=-SUM
      ENDIF
      IF(sum.lt.0d0 .or. ABS(SUM).LT.1.0D-4*SUMF) then
        SUM=S*DCSEL(RMU)/CS0
        !  temp fix until single scattering
        nterm = -1
      end if
 200  PDF=2.0D0*SUM
      RETURN
      END
      subroutine inigrd(ecut,emax,sflag)
      USE egs5_mscon_mod
      implicit none
      integer sflag
      double precision ecut, emax
      include 'include/egs5_h.f'
      include 'include/egs5_cdcsep.f'
      include 'include/egs5_uphiot.f'
      include 'include/egs5_useful.f'
      include 'pegscommons/mscom.f'
      integer i, j, k
      double precision e, edec, ecutev, emaxev
      ecutev = (ecut - RM) * 1.0d6
      emaxev = (emax - RM) * 1.0d6
      decade1 = dint(dlog10(ecutev))
      decade2 = dint(dlog10(emaxev))
      nmscate = 0
      nmscpw = 0
      do i = decade1,decade2
        edec = 10.d0**i
        k = 0
        do j = 1,NDEC
          k = k + 1
          e = edec * 10.d0**((j-1.d0)/NDEC)
          if(e.ge.ecutev) then
            if(nmscpw.eq.0) then
              e = ecutev
              if(k.eq.1) then
                joffset = NDEC-1
              else
                joffset = k - 2
              endif
            else
              e = edec * 10.d0**((j-2.d0)/NDEC)
            endif
            if(e.gt.emaxev) then
              if(mscate(nmscpw).ge.emaxev) go to 200
              e = emaxev
            endif
            nmscpw = nmscpw + 1
            mscate(nmscpw) = e
            if(e .le. 1.d8) nmscate = nmscate + 1
          endif
        end do
      end do
200   continue
      if(sflag.eq.1) call dcsstor(ecut,emax)
      if(th(2).eq.1.d-4) return
      I=1
      TH(I)=0.0D0
      THR(I)=TH(I)*PI/180.0D0
      XMU(I)=(1.0D0-COS(THR(I)))/2.0D0
      I=2
      TH(I)=1.0D-4
      THR(I)=TH(I)*PI/180.0D0
      XMU(I)=(1.0D0-COS(THR(I)))/2.0D0
   20 CONTINUE
      I=I+1
      IF(TH(I-1).LT.0.9999D-3) THEN
        TH(I)=TH(I-1)+2.5D-5
      ELSE IF(TH(I-1).LT.0.9999D-2) THEN
        TH(I)=TH(I-1)+2.5D-4
      ELSE IF(TH(I-1).LT.0.9999D-1) THEN
        TH(I)=TH(I-1)+2.5D-3
      ELSE IF(TH(I-1).LT.0.9999D+0) THEN
        TH(I)=TH(I-1)+2.5D-2
      ELSE IF(TH(I-1).LT.0.9999D+1) THEN
        TH(I)=TH(I-1)+1.0D-1
      ELSE IF(TH(I-1).LT.2.4999D+1) THEN
        TH(I)=TH(I-1)+2.5D-1
      ELSE
        TH(I)=TH(I-1)+5.0D-1
      ENDIF
      THR(I)=TH(I)*PI/180.0D0
      XMU(I)=(1.0D0-COS(THR(I)))/2.0D0
      IF(I.LT.NREDA) GO TO 20
      RETURN
      END
      SUBROUTINE egs5INTEG(X,A,B,C,D,XL,XU,SUM,N)
      implicit none
      integer n
      double precision  X(N),A(N),B(N),C(N),D(N), xl, xu, sum
      integer i, il, iu
      double precision x1,x2, xll, xuu, sign, sump
      IF(XU.GT.XL) THEN
        XLL=XL
        XUU=XU
        SIGN=1.0D0
      ELSE
        XLL=XU
        XUU=XL
        SIGN=-1.0D0
      ENDIF
      IF(XLL.LT.X(1).OR.XUU.GT.X(N)) THEN
        WRITE(506,10)
   10   FORMAT(5X,'Integral limits out of range. Stop.')
      ENDIF
      SUM=0.0D0
      CALL FINDI(X,XLL,N,IL)
      CALL FINDI(X,XUU,N,IU)
      IF(IL.EQ.IU) THEN
        X1=XLL
        X2=XUU
       
      ELSE
        X1=XLL
        X2=X(IL+1)
      
        IL=IL+1
        DO I=IL,IU
          X1=X(I)
          X2=X(I+1)
          IF(I.EQ.IU) X2=XUU
         
          SUM=SUM+SUMP
        ENDDO
      ENDIF
      SUM=SIGN*SUM
      RETURN
      END
      double precision function k1e(iq,e)
      implicit none
      include 'include/egs5_h.f'
      include 'pegscommons/k1spl.f'
      include 'pegscommons/dercon.f'
      integer iq
      double precision e
      integer i
      double precision ke
      ke = e - RM
      call findi(ehinge,ke,nhinge,i)
      if(iq .eq. -1) then
        k1e = ak1e(i)+ke*(bk1e(i)+ke*(ck1e(i)+ke*dk1e(i)))
      else
        k1e = ak1p(i)+ke*(bk1p(i)+ke*(ck1p(i)+ke*dk1p(i)))
      endif
      return
      end
      SUBROUTINE LEGENP(X,PL,NL)
      implicit none
      integer nl
      double precision x,  PL(NL)
      integer j
      double precision twox, f1, f2, d
      PL(1)=1.0D0
      PL(2)=X
      IF(NL.GT.2) THEN
        TWOX=2.0D0*X
        F1=X
        D=1.0D0
        DO J=3,NL
          F1=F1+TWOX
          F2=D
          D=D+1.0D0
          PL(J)=(F1*PL(J-1)-F2*PL(J-2))/D
        ENDDO
      ENDIF
      RETURN
      END
      subroutine makek1(emin,emax)
      USE PEGS_DCSSTR_MOD
      USE egs5_mscon_mod
      implicit none
      real*8 emax,emin
      real*8 g1pdedx, g1ededx, sumga
      external g1pdedx, g1ededx
      include 'include/egs5_h.f'
      include 'pegscommons/dercon.f'
      include 'pegscommons/mscom.f'
      include 'pegscommons/k1spl.f'
      integer i
      real*8 ekmax, ekmin, efrac, destp1, destp2, ctol
      real*8 e1, e2, ek1init(NESCPW), pk1init(NESCPW)
      ekmax = emax - RM
      ekmin = emin - RM
      nhinge = nmscpw
      destp2 = (efracl-efrach) / dlog(ekmin/ekmax)
      destp1 = efrach - dlog(ekmax) * destp2
      e2 = ekmax
      do i=1,nhinge-1
        efrac = destp1 + destp2 * dlog(e2)
        e1 = e2 * (1.d0 - efrac)
        if(efrac.ge.0.20) then
          ctol = 1.d-3
        else if(efrac.ge.0.05) then
          ctol = 1.d-4
        else
          ctol = 1.d-5
        endif
        pk1init(nhinge-i+1) = sumga(g1pdedx,e1+RM,e2+RM,ctol)
        ek1init(nhinge-i+1) = sumga(g1ededx,e1+RM,e2+RM,ctol)
        ehinge(nhinge-i+1) = e2
        e2 = mscate(nhinge-i)*1.d-6
      end do
      !-->  Extrapolate to get first hinge.
      ehinge(1) = ekmin
      efrac = (ehinge(2)-ehinge(1)) / (ehinge(3)-ehinge(2))
      pk1init(1) = pk1init(2) + (pk1init(2) - pk1init(3)) * efrac
      ek1init(1) = ek1init(2) + (ek1init(2) - ek1init(3)) * efrac
      do i = 1, nhinge
        if(ehinge(i) .le. 1.d8) then
          if(ek1init(i) .lt. k1mine) k1mine = ek1init(i)
          if(ek1init(i) .gt. k1maxe) k1maxe = ek1init(i)
          if(pk1init(i) .lt. k1minp) k1minp = pk1init(i)
          if(pk1init(i) .gt. k1maxp) k1maxp = pk1init(i)
        endif
      end do
      !-->  get splines for use later in PWLF1
100   call egs5spline(ehinge,ek1init,ak1e,bk1e,ck1e,dk1e,&
          0.d0,0.d0,nhinge)
      call egs5spline(ehinge,pk1init,ak1p,bk1p,ck1p,dk1p,&
          0.d0,0.d0,nhinge)
      return
      end
      double precision function addmol(x)
      implicit none
      double precision x
      include 'pegscommons/dercon.f'
      save
      addmol=1.0/(x-rm)**2
      return
      end
      double precision function adfmol(x)
      implicit none
      double precision x
      include 'pegscommons/dercon.f'
      save
      adfmol=-1.0/(x-RM)
      return
      end
      double precision function adimol(x)
      implicit none
      double precision x
      include 'pegscommons/dercon.f'
      save
      adimol=-1.0/x+RM
      return
      end


      subroutine adscpr(mxrawt,mxshet,elecnt,nshelt,capint,scprot,qcapt,&
      pz)
      implicit none
      integer iiend, mxrawt, i, j, ii, nshelt, mxshet
      double precision qcapt, elecnt, pz, capint, scprot
      include 'include/egs5_h.f'
      include 'pegscommons/cpcom.f'
      save
      iiend=MXRAW+mxrawt+1
      if (iiend.gt.MXMXRAW) then
        write(  526,100)
        write(  526,101)
        write(  526,102)
        write(  *,100)
        write(  *,101)
        write(  *,102)
100   format('Error Compton profile data  MXRAW .GT. MXMXRAW')
101   format('You have to decrease the number of elements per material')
102   format('or you have to set iprofr=0 in the [parameters] section')
        close(526)
        stop
      end if
      do i=1,31
        if (MXRAW.eq.0) then
          qcap(i)=qcapt(i)
        else
          if (qcap(i).ne.qcapt(i)) then
            write(  526,110)
            write(  *,110)
110         format(' Error Compton profile data  qcap are not agreed')
            close(526)
            stop
          end if
        end if
      end do
      do j=1,31
        if (MXRAW.eq.0) then
          scprof(j,iiend)=0.0
        else
          scprof(j,iiend)=scprof(j,MXRAW+1)
        end if
      end do
      do i=1,mxrawt
        ii=MXRAW+I
        nshell(ii)=nshelt(i)+MXSHEL
        elecni(ii)=elecnt(i)*pz
        capin(ii)=capint(i)
        do j=1,31
          scprof(j,ii)=scprot(j,i)
          scprof(j,iiend)=scprof(j,iiend)+scprot(j,i)*elecnt(i)*pz
        end do
      end do
      MXRAW=MXRAW+mxrawt
      mxshel=mxshel+mxshet
      return
      end
      double precision function affact(x)
      implicit none
      double precision aintp, x
      include 'pegscommons/cohcom.f'
      save
      affact=aintp(x,XVAL(1),100,AFAC2(1),1,.true.,.true.)
      return
      end
      double precision function aintp(x,xa,nx,ya,isk,xlog,ylog)
      implicit none
      integer nx, isk, j, i
      double precision xi, xj, xv, yi, yj
      double precision xa(nx),x
      double precision ya(isk,nx)
      logical xlog,ylog,xlogl
      save
      xlogl=xlog
      do j=2,nx
        if (x.lt.xa(j)) go to 100
      end do
      j=nx
100   i=j-1
      if (xa(i).le.0.0) then
        xlogl=.false.
      end if
      if (.not.xlogl) then
        xi=xa(i)
        xj=xa(j)
        xv=x
      else
        xi=dlog(xa(i))
        xj=dlog(xa(j))
        xv=dlog(x)
      end if
      if (ylog.and.(ya(1,i).eq.0.0.or.ya(1,j).eq.0.0)) then
        aintp=0.0
      else
        if (ylog) then
          yi=dlog(ya(1,i))
          yj=dlog(ya(1,j))
          if (xj.eq.xi) then
            aintp=yi
          else
            aintp=(yi*(xj-xv)+yj*(xv-xi))/(xj-xi)
          end if
          aintp=dexp(aintp)
        else
          yi=ya(1,i)
          yj=ya(1,j)
          if (xj.eq.xi) then
            aintp=yi
          else
            aintp=(yi*(xj-xv)+yj*(xv-xi))/(xj-xi)
          end if
        end if
      end if
      return
      end
      double precision function alin(x)
      implicit none
      double precision x
      save
      alin=x
      return
      end
      double precision function alini(x)
      implicit none
      double precision x
      save
      alini=x
      return
      end
      double precision function alke(e)
      implicit none
      double precision e
      include 'pegscommons/dercon.f'
      save
      alke=dlog(e-RM)
      return
      end
      double precision function alkei(x)
      implicit none
      double precision x, dexp
      include 'pegscommons/dercon.f'
      save
      alkei=dexp(x) + RM
      return
      end
      double precision function amoldm(en0,en)
      implicit none
      double precision en0, tm, em, betasq, amolfm, en
      include 'pegscommons/pmcons.f'
      include 'pegscommons/dercon.f'
      include 'pegscommons/molvar.f'
      include 'pegscommons/lamolm.f'
      save
      T0=en0-RM
      tm=T0/RM
      em=tm+1.
      c1=(tm/em)**2
      c2=(2.*tm+1.)/em**2
      betasq=1.-1./em**2
      cmoll=RLC*EDEN*2.*PI*R0**2/(betasq*T0*tm)
      amoldm=amolfm(en)
      return
      end
      double precision function amolfm(en)
      implicit none
      double precision t, en, eps, epsp, epsi, epspi
      include 'pegscommons/dercon.f'
      include 'pegscommons/lamolm.f'
      save
      t=en-RM
      eps=t/T0
      epsp=1.-eps
      epsi=1./eps
      epspi=1./epsp
      amolfm=CMOLL*(C1+epsi*(epsi-C2)+epspi*(epspi-C2))
      return
      end
      double precision function amolrm(en0,en1,en2)
      implicit none
      
      include 'pegscommons/pmcons.f'
      include 'pegscommons/dercon.f'
      include 'pegscommons/molvar.f'
      save
      t0=en0-RM
      t1=en1-RM
      t2=en2-RM
      tm=t0/RM
      em=tm+1.
      c1=(tm/em)**2
      c2=(2.*tm+1.)/em**2
      betasq=1.-1./em**2
      cmoll2=RLC*EDEN*2.*PI*R0**2/(betasq*tm)
      eps1=t1/t0
      epsp1=1.-eps1
      eps2=t2/t0
      epsp2=1.-eps2
      amolrm=cmoll2*(c1*(eps2-eps1)+1./eps1-1./eps2+1./epsp2-1./epsp1 -&
            c2*dlog(eps2*epsp1/(eps1*epsp2)))
      return
      end
      double precision function amoltm(e0)
      implicit none
      double precision e0, t0, amolrm
      include 'pegscommons/thres2.f'
      include 'pegscommons/dercon.f'
      save
      if (e0.le.THMOLL) then
        amoltm=0.
      else
        t0=e0-RM
        amoltm=amolrm(e0,ae,t0*0.5+RM)
      end if
      return
      end
      double precision function anihdm(e0,k)
      implicit none
      double precision gam, e0, t0p, anihfm
      include 'pegscommons/pmcons.f'
      include 'pegscommons/dercon.f'
      include 'pegscommons/molvar.f'
      include 'pegscommons/lanihm.f'
      double precision k
      save
      gam=e0/RM
      a=gam+1.
      t0p=gam-1.
      C1=RLC*EDEN*PI*R0**2/(a*t0p*RM)
      C2=A+2.0*gam/a
      anihdm=anihfm(k)
      return
      end
      double precision function anihfm(k)
      implicit none
      double precision s1
      include 'pegscommons/dercon.f'
      include 'pegscommons/lanihm.f'
      double precision k,kp,x
      save
      s1(x)=C1*(-1.+(C2-1.0/x)/x)
      kp=k/RM
      anihfm=s1(kp)+s1(A-kp)
      return
      end
      double precision function anihrm(e0,k1,k2)
      implicit none
      double precision s2, c1, c2, dlog, gam, e0, a, t0p
      include 'pegscommons/pmcons.f'
      include 'pegscommons/dercon.f'
      include 'pegscommons/molvar.f'
      double precision k1,k2,kp1,kp2
      double precision x
      save
      s2(x)=RM*c1*(-x+c2*dlog(x)+1.0/x)
      gam=e0/RM
      kp1=k1/RM
      kp2=k2/RM
      a=gam+1.
      T0P=gam-1.
      c1=RLC*EDEN*PI*R0**2/(A*T0P*RM)
      c2=A+2.*gam/A
      anihrm=s2(kp2)-s2(kp1)+s2(A-kp1)-s2(A-KP2)
      return
      end
      double precision function anihtm(e0)
      implicit none
      double precision gam, e0, p0p2, p0p, dsqrt, canih, dlog
      include 'pegscommons/pmcons.f'
      include 'pegscommons/dercon.f'
      include 'pegscommons/molvar.f'
      gam=e0/RM
      p0p2=gam*gam-1.0
      p0p=dsqrt(p0p2)
      canih=RLC*EDEN*PI*R0**2/(gam+1.)
      anihtm=canih*((gam*gam+4.*gam+1.)/p0p2*dlog(gam+p0p)-(gam+3.)/p0p)
      return
      end
      double precision function aprim(z,e)
      implicit none
      integer ie, naprz, napre, iz
      double precision e, em, aintp, z
      include 'include/egs5_h.f'
      include 'pegscommons/dercon.f'
      include 'pegscommons/epstar.f'
      double precision aprimd(115,14),eprim(115),zprim(14),aprimz(115)
      save
      if (IAPRIM.eq.0) then
        if (IAPRFL .eq. 0) then
          IAPRFL=1
          write(  526,100)
        end if
        if (e.ge.50) then
          APRIM=1.
        else
          em=e/RM
          do ie=1,18
            aprimz(ie)= aintp(z,zprim,5,aprimd(ie,1),115,.false.,.false.&
           )
          end do
          aprim=aintp(em,eprim,18,aprimz,1,.false.,.false.)
        end if
      else if(IAPRIM.eq.1) then
        if (IAPRFL.eq.0) then
          write(  526,110)
          read(522,*) naprz, napre
          read(522,*) (eprim(ie),ie=1,napre)
          do ie=1,napre
            eprim(ie)=1.+eprim(ie)/RM
          end do
          do iz=1,naprz
            read(522,*) zprim(iz),(aprimd(ie,iz),ie=1,napre)
          end do
          iaprfl=1
          rewind(522)
        end if
        em=e/RM
        do ie=1,napre
          aprimz(ie)= aintp(z,zprim,naprz,aprimd(ie,1),115,.true.,.false&
         .)
        end do
        aprim=aintp(em,eprim,napre,aprimz,1,.false.,.false.)
      else if (iaprim.eq.2) then
        if (iaprfl .eq. 0) then
          iaprfl=1
          write(  526,140)
        end if
        aprim=1.0
      else
        write(  526,150) iaprim
        write(  *,150) iaprim
150     format(//,' Illegal value for iaprim: ',I4)
        close(526)
        stop
      end if
      return
      end
      double precision function arec(x)
      implicit none
      double precision x
      save
      arec=1.0/x
      return
      end
      double precision function bhabdm(en0,en)
      implicit none
      double precision en0, tm, em, y, bhabfm, en
      include 'pegscommons/pmcons.f'
      include 'pegscommons/dercon.f'
      include 'pegscommons/molvar.f'
      include 'pegscommons/lbhabm.f'
      save
      T0=en0-RM
      tm=T0/RM
      em=tm+1.
      y=1./(tm+2.)
      betasi=1./(1.-1./em**2)
      CBHAB=RLC*EDEN*2.*PI*R0**2/(T0*tm)
      B1=2.-y**2
      B2=3.-y*(6.-y*(1.-y*2.))
      B3=2.-y*(10.-y*(16.-y*8.))
      B4=1.-y*(6.-y*(12.-y*8.))
      bhabdm=bhabfm(en)
      return
      end
      double precision function bhabfm(en)
      implicit none
      double precision t, en, eps, epsi
      include 'pegscommons/dercon.f'
      include 'pegscommons/lbhabm.f'
      save
      t=en-RM
      eps=t/T0
      epsi=1./eps
      bhabfm=CBHAB*(epsi*(epsi*BETASI-B1)+B2+EPS*(eps*B4-B3))
      return
      end
      double precision function bhabrm(en0,en1,en2)
      implicit none
      
      include 'pegscommons/pmcons.f'
      include 'pegscommons/dercon.f'
      include 'pegscommons/molvar.f'

      return
      end
      double precision function bhabtm(e0)
      implicit none
      double precision e0, bhabrm
      include 'pegscommons/thres2.f'
      include 'pegscommons/dercon.f'
      save
      if (e0.le.ae) then
        bhabtm=0.
      else
        bhabtm=bhabrm(e0,ae,e0)
      end if
      return
      end
      double precision function bremdr(ea,k)
      implicit none
      integer ls
      double precision ea, bremfr
      double precision k
      include 'pegscommons/lbremr.f'
      save
      e=ea
      if (e.ge.50.) then
        LD=2
        ls=3
      else
        LD=1
        ls=0
      end if
      LA=ls+1
      LB=ls+2
      bremdr=bremfr(k)
      return
      end
      double precision function bremdz(z,e,k)
      implicit none
      double precision brmsdz, z, e
      double precision k
      save
      bremdz=brmsdz(z,e,k)/k
      return
      end
      double precision function bremfr(k)
      implicit none
      double precision eps, del, delta, a, b
      double precision k
      include 'pegscommons/bremp2.f'
      include 'pegscommons/dbrpr.f'
      include 'pegscommons/lbremr.f'
      save
      eps=k/e
      del=eps/(e*(1-eps))
      if (del.gt.delpos(ld)) then
        bremfr=0.0
        return
      end if
      delta=DELCM*del
      if (delta.le.1.) then
        a=DL1(LA)+delta*(DL2(LA)+delta*DL3(LA))
        b =DL1(LB)+delta*(DL2(LB)+delta*DL3(LB))
      else
        a=DL4(LA)+DL5(LA)*dlog(delta+DL6(LA))
        b =DL4(LB)+DL5(LB)*dlog(delta+DL6(LB))
      end if
      bremfr=(ALPHI(LD)*(1.-eps)/eps/AL2*A+0.5*(2.*eps)*b)/e
      return
      end
      double precision function bremfz(k)
      implicit none
      double precision brmsfz
      double precision k
      save
      bremfz=brmsfz(k)/k
      return
      end
      double precision function bremrm(e,k1,k2)
      implicit none
      integer i
      double precision bremrz, e
      double precision k1,k2
      include 'include/egs5_h.f'
      include 'pegscommons/mixdat.f'
      include 'pegscommons/mxdatc.f'
      save
      bremrm=0.
      do i=1,ne
        bremrm=bremrm+pz(i)*bremrz(z(i),e,k1,k2)
      end do
      return
      end
      double precision function bremrr(e,k1,k2)
      implicit none
      double precision dummy, bremdr, e, qd
      double precision k1,k2
      external bremfr
      save
      dummy=bremdr(e,k1)
      bremrr=qd(bremfr,k1,k2,'bremfr')
      return
      end
      double precision function bremrz(z,e,k1,k2)
      implicit none
      double precision dummy, bremdz, z, e, qd
      double precision k1,k2
      external bremfz
      save
      dummy=bremdz(z,e,k1)
      bremrz=qd(bremfz,k1,k2,'bremfz')
      return
      end
      double precision function bremtm(e0)
      implicit none
      double precision e0, bremrm
      include 'pegscommons/thres2.f'
      include 'pegscommons/dercon.f'
      save
      if (e0.le.AP+RM) then
        bremtm=0.
      else
        bremtm=bremrm(e0,AP,e0-RM)
      end if
      return
      end
      double precision function bremtr(e0)
      implicit none
      double precision e0, bremrr
      include 'pegscommons/thres2.f'
      include 'pegscommons/dercon.f'
      save
      if (e0.le.AP+RM) then
        bremtr=0.
      else
        bremtr=bremrr(e0,AP,e0-RM)
      end if
      return
      end
      double precision function brmsdz(z,ea,k)
      implicit none
      double precision ea, z, aprim, xsif, dlog, fcoulc, brmsfz
      double precision k
      include 'pegscommons/pmcons.f'
      include 'pegscommons/dercon.f'
      include 'pegscommons/molvar.f'
      include 'pegscommons/lbremz.f'
      save
      e=ea
      delc=136.*z**(-1./3.)*RM/e
      const=aprim(z,e)*(AN*RHO/WM)*R0**2*FSC*z*(Z+XSIF(Z))*RLC
      XLNZ=4./3.*dlog(z)
      if (e.ge.50) XLNZ=XLNZ+4.*fcoulc(z)
      DELTAM=dexP((21.12-XLNZ)/4.184)-0.952
      brmsdz=brmsfz(k)
      return
      end
      double precision function brmsfz(k)
      implicit none
      double precision emkloc, delta, sb1, sb2, dlog, ee
      double precision k
      include 'pegscommons/lbremz.f'
      save
      emkloc=E-k
      if (emkloc.eq.0.0) then
        emkloc=1.D-25
      end if
      delta=DELC*k/emkloc
      if (delta.ge.DELTAM) then
        brmsfz=0.0
      else
        if (delta.le.1.) then
          sb1=20.867+delta*(-3.242+delta*0.625)-XLNZ
          sb2=20.209+delta*(-1.930+delta*(-0.086))-XLNZ
        else
          sb1=21.12-4.184*dlog(delta+0.952)-XLNZ
          sb2=sb1
        end if
        ee=emkloc/E
        brmsfz=CONST*((1.+ee*ee)*sb1-0.666667*ee*sb2)
      end if
      return
      end
      double precision function brmsrm(e,k1,k2)
      implicit none
      integer i
      double precision brmsrz, e
      double precision k1,k2
      include 'include/egs5_h.f'
      include 'pegscommons/mixdat.f'
      include 'pegscommons/mxdatc.f'
      save
      brmsrm=0.
      do i=1,ne
        brmsrm=brmsrm+PZ(I)*brmsrz(z(i),e,k1,k2)
      end do
      return
      end
      double precision function brmsrz(z,e,k1,k2)
      implicit none
      double precision dummy, brmsdz, z, e, qd
      double precision k1,k2
      external brmsfz
      save
      dummy=brmsdz(z,e,k1)
      brmsrz=qd(brmsfz,k1,k2,'brmsfz')
      return
      end
      double precision function brmstm(e0,eg)
      implicit none
      double precision e0, au, eg, brmsrm
      include 'pegscommons/dercon.f'
      save
      if (e0.le.RM) then
        brmstm=0.
      else
        au=dmin1(eg,e0-RM)
        brmstm=brmsrm(e0,0.D0,AU)
      end if
      return
      end
      subroutine cfuns(e,v)
      implicit none
      double precision aintp, e
      include 'include/egs5_h.f'
      include 'pegscommons/bcom.f'
      include 'pegscommons/cpcom.f'
      double precision v(1)
      save
      v(1)=aintp(e,QCAP(1),31,AVCPRF(1),1,.true.,.true.)
      return
      end
      subroutine cfuns2(e,v)
      implicit none
      double precision aintp, e
      include 'include/egs5_h.f'
      include 'pegscommons/bcom.f'
      include 'pegscommons/cpcom.f'
      double precision v(1)
      save
      v(1)=aintp(e,CPROFI(1),301,QCAP10(1),1,.true.,.false.)
      return
      end
      subroutine cfuns3(e,v)
      implicit none
      integer ishell
      double precision aintp, e
      include 'include/egs5_h.f'
      include 'pegscommons/bcom.f'
      include 'pegscommons/cpcom.f'
      double precision v(MXMXRAW)
      save
      do ishell=1,MXSHEL
      end do
      return
      end
      subroutine cfuns4(e,v)
      implicit none
      integer ishell
      double precision aintp, e
      include 'include/egs5_h.f'
      include 'pegscommons/bcom.f'
      include 'pegscommons/cpcom.f'
      double precision v(MXMXRAW)
      save
      do ishell=1,mxshel
        v(ishell)=aintp(e,QCAP(1),31,SCPSUM(1,ISHELL),1,.true.,.true.)
      end do
      return
      end
      double precision function cohetm(k)
      implicit none
      integer i
      double precision cohetz, cohetzint
      double precision k
      include 'include/egs5_h.f'
      include 'pegscommons/mixdat.f'
      include 'pegscommons/mxdatc.f'
      include 'pegscommons/cohcom.f'
      save
      cohetm=0.d0
      if (irayl.eq.1) then
        do i=1,NE
          cohetm=cohetm+PZ(I)*cohetz(z(i),k)
        end do
        return
      else if(irayl.eq.2) then
        cohetm=cohetzint(1.d0,k)
        return
      end if
      end
      double precision function cohetz(z,k)
      implicit none
      integer iz
      double precision pcon, z, aintp
      double precision k
      include 'pegscommons/molvar.f'
      include 'pegscommons/phpair.f'
      include 'pegscommons/pmcons.f'
      include 'pegscommons/cohcom.f'
      save
      pcon= 1.D-24*(AN*RHO/WM)*RLC
      iz=z
      cohetz=pcon*aintp(k,PHE(1,iz),NPHE(iz),COHE(1,iz),1,.true.,.true.)
      return
      end
      double precision function cohetzint(z,k)
      implicit none
      integer iz
      double precision pcon, z, aintp
      double precision k
      include 'pegscommons/molvar.f'
      include 'pegscommons/phpair.f'
      include 'pegscommons/pmcons.f'
      include 'pegscommons/cohcom.f'
      save
      pcon= 1.D-24*(AN*RHO/WM)*RLC
      iz=z
      return
      end
      double precision function compdm(k0a,k)
      implicit none
      include 'pegscommons/pmcons.f'
      include 'pegscommons/dercon.f'
      include 'pegscommons/molvar.f'
      include 'pegscommons/lcompm.f'
      double precision k0a,k0p,compfm, k
      save
      k0=k0a
      k0p=k0/RM
      CCOMP=RLC*EDEN*PI*R0**2/(k0*k0p)
      C1=1./k0p**2
      C2=1.-(2.+2.*k0p)/k0p**2
      C3=(1.+2.*k0p)/k0p**2
      compdm=compfm(k)
      return
      end
      double precision function compfm(k)
      implicit none
      double precision eps, epsi
      double precision k
      include 'pegscommons/lcompm.f'
      save
      eps=k/K0
      epsi=1./eps
      compfm=CCOMP*( (C1*epsi+C2)*epsi+C3+eps )
      return
      end
      double precision function comprm(k0,k1,k2)
      implicit none
      double precision ccomp2
      include 'pegscommons/pmcons.f'
      include 'pegscommons/dercon.f'
      include 'pegscommons/molvar.f'
      double precision k0,k1,k2
      real*8 c1,c2,c3,eps1,eps2,k0p
      save
      k0p=k0/RM
      ccomp2=RLC*EDEN*PI*R0**2/k0p
      c1=1./k0p**2
      c2=1.-(2.+2.*k0p)/k0p**2
      c3=(1.+2.*k0p)/k0p**2
      eps1=k1/k0
      eps2=k2/k0
      return
      end
      double precision function comptm(k0)
      implicit none
      integer i, iz
      double precision pcon, comsum, aintp, comprm
      include 'include/egs5_h.f'
      include 'pegscommons/bcom.f'
      include 'pegscommons/dercon.f'
      include 'pegscommons/mixdat.f'
      include 'pegscommons/mxdatc.f'
      include 'pegscommons/phpair.f'
      include 'pegscommons/pmcons.f'
      include 'pegscommons/molvar.f'
      double precision k0,k1
      save
      if (IBOUND.eq.1) then
        pcon=1.D-24*(AN*RHO/WM)*RLC
        comsum=0.0
        do i=1,NE
          iz=Z(i)
        end do
        comptm=pcon*comsum
      else
        K1=K0*RM/(RM+2.*K0)
        comptm=comprm(K0,K1,K0)
      end if
      return
      end
      double precision function cprfil(x)
      implicit none
      double precision aintp, x
      include 'include/egs5_h.f'
      include 'pegscommons/cpcom.f'
      save
      cprfil=aintp(x,QCAP(1),31,AVCPRF(1),1,.true.,.true.)
      return
      end
      double precision function cratio(e)
      implicit none
      double precision tot, pairtu, e, comptm, photte, cohetm
      save
      tot=pairtu(e)+comptm(e)+photte(e)
      cratio=tot/(tot+cohetm(e))
      return
      end
      double precision function dcadre(f,a,b,aerr,rerr,error,ier)
      implicit none

      external f
      dimension t(10,10),r(10),ait(10),dif(10),rn(4),ts(2049)
      dimension ibegs(30),begin(30),finis(30),est(30)
      dimension reglsv(30)
      logical h2conv,aitken,right,reglar,reglsv
      double precision t,r,ait,dif,rn,ts,begin,finis,est,aitlow
      double precision h2tol,aittol,length,jumptl,zero,p1,half,one
      double precision two,four,fourp5,ten,hun,cadre,error,a,b
      double precision aerr,rerr,stepmn,stepnm,stage,curest,fnsize
      double precision prever,beg,fbeg,end,fend,step,astep,tabs,hovn
      double precision fn,sum,sumabs,absi,vint,tabtlm,ergl,ergoal
      double precision erra,errr,fextrp,errer,diff,sing,fextm1,alg4o2
      double precision h2nxt,singnx,slope,fbeg2,alpha
      double precision erret,h2tfex,fi
      double precision rval,f
      save
      alg4o2=dlog10(TWO)
      cadre=zero
      error=zero
      curest=zero
      vint=zero
      ier=0
      length=dabs(B-A)
      if (length.eq.zero) go to 215
      if (rerr.gt.p1.or.rerr.lt.zero) go to 210
      if (aerr.eq.zero.and.(rerr+hun).le.hun) go to 210
      errr=rerr
      erra=dabs(AERR)
      stepmn=(length/float(2**mxstge))
      stepnm=dmax1(length,dabs(A),abs(B))*TEN
      stage=half
      istage=1
      fnsize=zero
      prever=zero
      reglar=.false.
      beg=A
      rval=beg
      fbeg=f(rval)*half
      ts(1)=fbeg
      ibeg=1
      end=B
      rval=end
      fend=f(rval)*half
      ts(2)=fend
      iend=2
5     right=.false.
10    step=end - beg
      astep=dabs(step)
      if (astep.lt.stepmn) go to 205
      if (stepnm+astep.eq.stepnm) go to 205
      t(1,1)=fbeg + fend
      tabs=dabs(fbeg) + dabs(fend)
      l=1
      n=1
      h2conv=.false.
      aitken=.false.
15    lm1=l
      l=l + 1
      n2=n + n
      fn=n2
      istep=(iend - ibeg)/n
      if (istep.gt.1) go to 25
      ii=iend
      iend=iend + n
      if (iend.gt.maxts) go to 200
      hovn=step/fn
      iii=iend
      fi=one
      do i=1,n2,2
        ts(iii)=ts(ii)
        rval=end-fi*hovn
        ts(iii-1)=f(rval)
        fi=fi+two
        iii=iii-2
        ii=ii-1
      end do
      istep=2
25    istep2=ibeg + istep/2
      sum=zero
      sumabs=zero
      do i=istep2,iend,istep
        sum=sum + ts(i)
        sumabs=sumabs + dabs(ts(i))
      end do
      t(l,1)=t(l-1,1)*half+sum/fn
      tabs=tabs*half+sumabs/fn
      absi=astep*tabs
      n=n2
      it=1
      vint=step*t(l,1)
      tabtlm=tabs*ten
      fnsize=dmax1(fnsize,dabs(t(l,1)))
      ergl=astep*fnsize*ten
      ergoal=stage*dmax1(erra,errr*dabs(curest+vint))
      fextrp=one
      do i=1,lm1
        fextrp=fextrp*four
        t(i,l)=t(l,i) - t(l-1,i)
        t(l,i+1)=t(l,i) + t(i,l)/(fextrp-one)
      end do
      errer=astep*dabs(t(1,l))
      if (l.gt.2) go to 40
      if (tabs+p1*dabs(t(1,2)).eq.tabs) go to 135
      go to 15
40    do i=2,lm1
      diff=zero
      if (tabtlm+dabs(t(i-1,l)).ne.tabtlm) diff=t(i-1,lm1)/t(i-1,l)
      t(i-1,lm1)=diff
      end do
      if (dabs(four-t(1,lm1)).le.h2tol) go to 60
      if (t(1,lm1).eq.zero) go to 55
      if (dabs(two-abs(t(1,lm1))).lt.jumptl) go to 130
      if (l.eq.3) go to 15
      h2conv=.false.
      if (dabs((t(1,lm1)-t(1,l-2))/t(1,lm1)).le.aittol) go to 75
50    if (reglar) go to 55
      if (l.eq.4) go to 15
55    if(errer.gt.ergoal.and.(ergl+errer).ne.ergl) go to 175
      go to 145
60    if(h2conv) go to 65
      aitken=.false.
      h2conv=.true.
65    fextrp=four
70    it=it + 1
      vint=step*t(l,it)
      errer=dabs(step/(fextrp-one)*t(it-1,l))
      if (errer.le.ergoal) go to 160
      if (ergl+errer.eq.ergl) go to 160
      if (it.eq.lm1) go to 125
      if (t(it,lm1).eq.zero) go to 70
      if (t(it,lm1).le.fextrp) go to 125
75    if(t(1,lm1).lt.aitlow) go to 175
      if (aitken) go to 80
      h2conv=.false.
      aitken=.true.
80    fextrp=t(l-2,lm1)
      if (fextrp.gt.fourp5) go to 65
      if (fextrp.lt.aitlow) go to 175
      if (dabs(fextrp-t(l-3,lm1))/t(1,lm1).gt.h2tol) go to 175
      sing=fextrp
      fextm1=one/(fextrp - one)
      ait(1)=zero
      do i=2,l
      ait(i)=t(i,1) + (t(i,1)-t(i-1,1))*fextm1
      r(i)=t(1,i-1)
      dif(i)=ait(i) - ait(i-1)
      end do
      it=2
90    vint=step*ait(l)
      errer=errer*fextm1
      if (errer.gt.ergoal.and.(ergl+errer).ne.ergl) go to 95
      alpha=dlog10(sing)/alg4o2 - one
      ier=max0(ier,65)
      go to 160
95    it=it + 1
      if (it.eq.lm1) go to 125
      if (it.gt.3) go to 100
      h2nxt=four
      singnx=sing+sing
100   if(h2nxt.lt.singnx) go to 105
      fextrp=singnx
      singnx=singnx+singnx
      go to 110
105   fextrp=h2nxt
      h2nxt=four*h2nxt
110   do i=it,lm1
      r(i+1)=zero
      if (tabtlm+dabs(dif(i+1)).ne.tabtlm) r(i+1)=dif(i)/dif(i+1)
      end do
      h2tfex=-h2tol*fextrp
      if (r(l)-fextrp.lt.h2tfex) go to 125
      if (r(l-1)-fextrp.lt.h2tfex) go to 125
      errer=astep*dabs(dif(l))
      fextm1=one/(fextrp - one)
      do i=it,l
      ait(i)=ait(i) + dif(i)*fextm1
      dif(i)=ait(i) - ait(i-1)
      end do
      go to 90
125   fextrp=dmax1(prever/errer,aitlow)
      prever=errer
      if (l.lt.5) go to 15
      if (l-it.gt.2.and.istage.lt.mxstge) go to 170
      erret=errer/(fextrp**(maxtbl-l))
      if (erret.gt.ergoal.and.(ergl+erret).ne.ergl) go to 170
      go to 15
130   if(errer.gt.ergoal.and.(ergl+errer).ne.ergl) go to 170
      diff=dabs(t(1,l))*(fn+fn)
      go to 160
135   slope=(fend-fbeg)*two
      fbeg2=fbeg+fbeg
      do i=1,4
      rval=beg+rn(i)*step
      diff=dabs(f(rval) - fbeg2-rn(i)*slope)
      if (tabtlm+diff.ne.tabtlm) go to 155
      end do
      go to 160
145   slope=(fend-fbeg)*two
      fbeg2=fbeg+fbeg
      i=1
150   rval=beg+rn(i)*step
      diff=dabs(f(rval) - fbeg2-rn(i)*slope)
155   errer=dmax1(errer,astep*diff)
      if (errer.gt.ergoal.and.(ergl+errer).ne.ergl) go to 175
      i=i+1
      if (i.le.4) go to 150
      ier=66
160   cadre=cadre + vint
      error=error + errer
      if (right) go to 165
      istage=istage - 1
      if (istage.eq.0) go to 220
      reglar=reglsv(istage)
      beg=begin(istage)
      end=finis(istage)
      curest=curest - est(istage+1) + vint
      iend=ibeg - 1
      fend=ts(iend)
      ibeg=ibegs(istage)
      go to 180
165   curest=curest + vint
      stage=stage+stage
      iend=ibeg
      ibeg=ibegs(istage)
      end=beg
      beg=begin(istage)
      fend=fbeg
      fbeg=ts(ibeg)
      go to 5
170   reglar=.true.
175   if(istage.eq.mxstge) go to 205
      if (right) go to 185
      reglsv(istage+1)=reglar
      begin(istage)=beg
      ibegs(istage)=ibeg
      stage=stage*half
180   right=.true.
      beg=(beg+end)*half
      ibeg=(ibeg+iend)/2
      ts(ibeg)=ts(ibeg)*half
      fbeg=ts(ibeg)
      go to 10
185   nnleft=ibeg - ibegs(istage)
      if (iend+nnleft.ge.maxts) go to 200
      iii=ibegs(istage)
      ii=iend
      do i=iii,ibeg
      ii=ii + 1
      ts(ii)=ts(i)
      end do
      do i=ibeg,ii
      ts(iii)=ts(i)
      iii=iii + 1
      end do
      iend=iend + 1
      ibeg=iend - nnleft
      fend=fbeg
      fbeg=ts(ibeg)
      finis(istage)=end
      end=beg
      beg=begin(istage)
      begin(istage)=end
      reglsv(istage)=reglar
      istage=istage + 1
      reglar=reglsv(istage)
      est(istage)=vint
      curest=curest + est(istage)
      go to 5
200   ier=131
      go to 215
205   ier=132
      go to 215
210   ier=133
215   cadre=curest + vint
220   dcadre=cadre
      return
      end
      subroutine differ
      implicit none
      
      include 'pegscommons/molvar.f'
      include 'pegscommons/bremp2.f'
      include 'pegscommons/dbrpr.f'
      include 'pegscommons/dercon.f'
      include 'pegscommons/radlen.f'
      save
      al2 = dlog(2.D0)
      al183= dlog(a183)
      alphi(1)= al2*(4./3. + 1./(9.*al183*(1.+ZP)))
      alphi(2)= al2*(4./3. + 1./(9.*al183*(1.+ZU)))
      alfp1(1)= 2./3. - 1./(36.*al183*(1.+ZP))
      alfp1(2)= 2./3. - 1./(36.*al183*(1.+ZU))
      alfp2(1)= (1./12.)*(4./3. + 1./(9.*al183*(1+ZP)))
      alfp2(2)= (1./12.)*(4./3. + 1./(9.*al183*(1+ZU)))
      bpar(1)= alfp1(1)/(alfp1(1)+alfp2(1))
      bpar(2)= alfp1(2)/(alfp1(2)+alfp2(2))
      delcm= 136.0*dexp(ZG)*RM
      delpos(1)= (dexp((21.12+4.*ZG)/4.184)-0.952)/DELCM
      delpos(2)= (dexp((21.12+4.*ZV)/4.184)-0.952)/DELCM
      f10=4.*al183
      f20=f10 - 2./3.
      a1den =3.0*f10- f20 + 8.0*ZG
      a2den =3.0*f10- f20 + 8.0*ZV
      b1den = f10 + 4.0*ZG
      b2den = f10 + 4.0*ZV
      c1den = 3.0*f10+ f20 + 16.0*ZG
      c2den = 3.0*f10+ f20 + 16.0*ZV
      dl1(1)= (3.0*20.867-20.209+8.0*ZG)/a1den
      dl2(1)= (3.0*(-3.242)-(-1.930))/a1den
      dl3(1)= (3.0*(0.625)-(0.086))/a1den
      dl4(1)= (2.0*21.12+8.0*ZG)/a1den
      dl5(1)= 2.0*(-4.184)/a1den
      dl6(1)= 0.952
      dl1(4)= (3.0*20.867-20.209+8.0*ZV)/a2den
      dl2(4)= (3.0*(-3.242)-(-1.930))/a2den
      dl3(4)= (3.0*(0.625)-(0.086))/a2den
      dl4(4)= (2.0*21.12+8.0*ZV)/a2den
      dl5(4)= 2.0*(-4.184)/a2den
      dl6(4)= 0.952
      dl1(2)= (20.867+4.0*ZG)/b1den
      dl2(2)= -3.242/b1den
      dl3(2)= 0.625/b1den
      dl4(2)= (21.12+4.0*ZG)/b1den
      dl5(2)= -4.184/b1den
      dl6(2)= 0.952
      dl1(5)= (20.867+4.0*ZV)/b2den
      dl2(5)= -3.242/b2den
      dl3(5)= 0.625/b2den
      dl4(5)= (21.12+4.0*ZV)/b2den
      dl5(5)= -4.184/b2den
      dl6(5)= 0.952
      dl1(3)= (3.0*20.867+20.209+16.0*ZG)/c1den
      dl2(3)= (3.0*(-3.242)+(-1.930))/c1den
      dl3(3)= (3.0*0.625+(-0.086))/c1den
      dl4(3)= (4.0*21.12+16.0*ZG)/c1den
      dl5(3)= 4.0*(-4.184)/c1den
      dl6(3)= 0.952
      dl1(6)= (3.0*20.867+20.209+16.0*ZV)/c2den
      dl2(6)= (3.0*(-3.242)+(-1.930))/c2den
      dl3(6)= (3.0*0.625+(-0.086))/c2den
      dl4(6)= (4.0*21.12+16.0*ZV)/c2den
      dl5(6)= 4.0*(-4.184)/c2den
      dl6(6)= 0.952
      write(  526,100)
      return
      end
      double precision function ebind(e)
      implicit none
      integer i, j
      double precision phottz, e, stot, photte
      include 'include/egs5_h.f'
      include 'pegscommons/mixdat.f'
      include 'pegscommons/mxdatc.f'
      include 'pegscommons/phpair.f'
      save
      ebind=0.0
      do i=1,NE
        j=z(i)
        ebind=ebind+PZ(i)*phottz(z(i),e)*EKEDGE(j)*0.001
      end do
      stot=photte(e)
      if (stot.ne.0.0) ebind=ebind/stot
      return
      end
      double precision function ebr1(e)
      implicit none
      double precision brem, bremtm, e, tebr, amoltm
      save
      brem=bremtm(e)
      tebr=brem+amoltm(e)
      if (tebr.gt.0.0) then
        ebr1=brem/tebr
      else
        ebr1=0.0
      end if
      return
      end
      double precision function ededx(e)
      implicit none
      double precision sptote, e
      include 'pegscommons/thres2.f'
      save
      ededx=sptote(e,ae,ap)
      return
      end
      subroutine efuns(e,v)
      implicit none

      double precision v(15)
      include 'pegscommons/thres2.f'
      include 'pegscommons/legacy.f'
      save
      if (iunrst.eq.0 .or. iunrst.eq.1 .or. iunrst.eq.5) then
        brem=bremtm(e)
        amoll=amoltm(e)
        bhab=bhabtm(e)
        annih=anihtm(e)
        esig=brem+amoll
        v(1)=esig
        psig=brem+bhab+annih
        v(2)=psig
        v(3)=sptote(e,AE,AP)
        v(4)=sptotp(e,AE,AP)
        if (esig.gt.0.0) then
          v(5)=brem/esig
        else
          if (thbrem.le.thmoll) then
            v(5)=1.0
          else
            v(5)=0.0
          end if
        end if
        v(6)=brem/psig
        v(7)=(brem+bhab)/psig
        v(8)=tmxs(e)
        v(9) = g1e(-1,e,0)
        v(10) = g1e(+1,e,0)
        if(oldK1run) then
          v(11) = k1e(-1,e)
          v(12) = k1e(+1,e)
        else
          v(11) = 0.d0
          v(12) = 0.d0
        endif
        v(13) = csdar(-1,e)
        v(14) = csdar(+1,e)
        v(15) = estepmax(e)
      else if (iunrst.eq.2) then
        v(1)=0.0
        v(2)=0.0
        v(5)=0.0
        v(6)=0.0
        v(7)=0.0
        v(3) = sptote(e,e,e)
        v(4) = sptotp(e,e,e)
        v(8) = tmxs(e)
        v(9) = g1e(-1,e,0)
        v(10) = g1e(+1,e,0)
        if(oldK1run) then
          v(11) = k1e(-1,e)
          v(12) = k1e(+1,e)
        else
          v(11) = 0.d0
          v(12) = 0.d0
        endif
        v(13) = csdar(-1,e)
        v(14) = csdar(+1,e)
        v(15) = estepmax(e)
      else if (iunrst.eq.3) then
        brem=bremtm(e)
        annih=anihtm(e)
        v(1)=brem
        v(2)=brem + annih
        v(3)=sptote(e,e,AP)
        v(4)=sptotp(e,e,AP)
        v(5)=1.0
        v(6)=brem/v(2)
        v(7)=v(6)
        v(8)=tmxs(e)
        v(9) = g1e(-1,e,0)
        v(10) = g1e(+1,e,0)
        if(oldK1run) then
          v(11) = k1e(-1,e)
          v(12) = k1e(+1,e)
        else
          v(11) = 0.d0
          v(12) = 0.d0
        endif
        v(13) = csdar(-1,e)
        v(14) = csdar(+1,e)
        v(15) = estepmax(e)
      else if (iunrst.eq.4) then
        v(1)=amoltm(e)
        v(2)=bhabtm(e)
        v(3)=sptote(e,AE,e)
        v(4)=sptotp(e,AE,e)
        v(5)=0.0
        v(6)=0.0
        v(7)=1.0
        v(8)=tmxs(e)
        v(9) = g1e(-1,e,0)
        v(10) = g1e(+1,e,0)
        if(oldK1run) then
          v(11) = k1e(-1,e)
          v(12) = k1e(+1,e)
        else
          v(11) = 0.d0
          v(12) = 0.d0
        endif
        v(13) = csdar(-1,e)
        v(14) = csdar(+1,e)
        v(15) = estepmax(e)
      else
        write(  526,100) iunrst
        write(  *,100) iunrst
        close(526)
        stop
      end if
      return
      end
      subroutine eiifuns(e,v)
      implicit none
      integer i
      double precision amoll, amoltm, e, eiisum, zval, eiitm
      include 'include/egs5_h.f'
      double precision v(MXEL)
      include 'pegscommons/mixdat.f'
      include 'pegscommons/mxdatc.f'
      save
      do i=1,MXEL
        v(i)=0.0
      end do
      amoll=amoltm(e)
      if (amoll.lt.1.0d-30) return
      eiisum=0.0
      do i=1,NE
        zval=Z(I)
        eiisum=eiisum+eiitm(e,zval)*PZ(i)
        v(i)=eiisum/amoll
      end do
      return
      end
      double precision function eiitm(e,zval)
      implicit none
      integer j, nismall
      
      include 'include/egs5_h.f'
      include 'pegscommons/dercon.f'
      include 'pegscommons/eimpact.f'
      include 'pegscommons/phpair.f'
      include 'pegscommons/pmcons.f'
      include 'pegscommons/mixdat.f'
      include 'pegscommons/mxdatc.f'
      include 'pegscommons/molvar.f'
      save
      J=ZVAL
      ekbmev=ekedge(j)*0.001
      if (ekbmev.eq.0.0) then
        eiitm=0.0
        return
      end if
      x=(e-RM)/ekbmev
      nismall=2
      if (x.gt.1.001) then
        capi=ekedge(j)/RM/1000.0
        cape=(e-RM)/RM
        fr1=(2.0+capi)/(2.0+cape)
        fr2=(1.0+cape)/(1.0+capi)
        fr3=(capi+cape)*(2.0+cape)*(1.0+capi)**2
        fr4=cape*(2.0+cape)*(1.0+capi)**2+capi*(2.0+capi)
        rfact=fr1*fr2**2*(fr3/fr4)**1.5
        if (impact.eq.1) then
          smalla0=5.292E3
          capi0=13.606D-3
          capu=(e-RM)/(EKEDGE(j)*0.001)
          smalld0=-0.0318
          smalld1=0.3160
          smalld2=-0.1135
          smallb0=10.57
          smallb1=-1.736
          smallb2=0.317
          sphi=(EKEDGE(j)/capi0)**(smalld0+smalld1/capu+smalld2/capu**2)
          spsi=smallb0*dexp(smallb1/capu+smallb2/capu**2)
        end if
        if (impact.eq.2) then
          cape=(e-RM)/RM
          cape1=cape+1.0
          cape2=cape+2.0
          capi=ekbmev/RM
          smalph=1.0/137.036
          eke0=0.5*(smalph*ZVAL)**2*RM*1000.0
          qcapa=cape1*cape1/capi/cape/cape2
        end if
        if (impact.eq.3) then
          cape=(e-RM)/RM
          cape1=cape+1.0
          cape2=cape+2.0
          capi=ekbmev/RM
          qcapa=cape1*cape1/capi/cape/cape2
        end if
        if (impact.eq.4) then
          beta2a=(1.0+(e-RM)/RM)**(-2)
          beta2=1.0-beta2a
          beta02=1.0-(1.0+EKEDGE(j)/(RM*1000))**(-2)
          fcap1=254.9/(EKEDGE(j)*beta2)
          fcap2=dlog(beta2/beta2a)-beta2
          fcap3=1.0-beta02/beta2
          fcap4=dlog(1.0/beta02)
          fcap5=beta02/beta2
          sma=5.14*ZVAL**(-0.48)
          smb=5.76-0.04*ZVAL
          smc=0.72+0.039*ZVAL-0.0006*ZVAL**2
          qcap=sma*fcap1*(fcap2+smb*fcap3+fcap4*fcap5**smc)
        end if
        if (impact.eq.5.or.impact.eq.6) then
          qconst=0.0656
          g1=1.0/x*((x-1.0)/(x+1.0))**1.5
          g2=1.0+0.6667*(1.0-0.5/x)*dlog(2.7+dsqrt(x-1.d0))
          qcap=qconst*nismall/ekbmev**2*g1*g2
          if (impact.eq.6) then
            qcap=qcap*rfact
          end if
        end if
        if (qcap.lt.0.0) then
          qcap=0.0
        end if
        eiitm=qcap*AN*1.0D-24/WM*RHO*RLC
      else
        eiitm=0.0
      end if
      return
      end
      double precision function esig(e)
      implicit none
      double precision bremtm, e, amoltm
      save
      esig=bremtm(e)+amoltm(e)
      return
      end
      double precision function fcoulc(z)
      implicit none
      double precision asq, z
      include 'pegscommons/dercon.f'
      save
      asq=(FSC*z)**2
      return
      end
      double precision function fi(i,x1,x2,x3,x4)
      implicit none
      integer i
     
      save
1     fi=alin(x1)
      return
2     fi=alini(x1)
      return
3     fi=adfmol(x1)
      return
4     fi=adimol(x1)
      return
5     fi=addmol(x1)
      return
6     fi=dlog(x1)
      return
7     fi=dexp(x1)
      return
8     fi=arec(x1)
      return
9     fi=alke(x1)
      return
10    fi=alkei(x1)
      return
11    fi=amoldm(x1,x2)
      return
12    fi=amolfm(x1)
      return
13    fi=amolrm(x1,x2,x3)
      return
14    fi=amoltm(x1)
      return
15    fi=anihdm(x1,x2)
      return
16    fi=anihfm(x1)
      return
17    fi=anihrm(x1,x2,x3)
      return
18    fi=anihtm(x1)
      return
19    fi=aprim(x1,x2)
      return
20    fi=bhabdm(x1,x2)
      return
21    fi=bhabfm(x1)
      return
22    fi=bhabrm(x1,x2,x3)
      return
23    fi=bhabtm(x1)
      return
24    fi=bremdr(x1,x2)
      return
25    fi=bremfr(x1)
      return
26    fi=bremdz(x1,x2,x3)
      return
27    fi=brmsdz(x1,x2,x3)
      return
28    fi=bremfz(x1)
      return
29    fi=brmsfz(x1)
      return
30    fi=bremrr(x1,x2,x3)
      return
31    fi=bremrm(x1,x2,x3)
      return
32    fi=bremrz(x1,x2,x3,x4)
      return
33    fi=bremtm(x1)
      return
34    fi=bremtr(x1)
      return
35    fi=brmsrm(x1,x2,x3)
      return
36    fi=brmsrz(x1,x2,x3,x4)
      return
37    fi=brmstm(x1,x2)
      return
38    fi=cohetm(x1)
      return
39    fi=cohetz(x1,x2)
      return
40    fi=compdm(x1,x2)
      return
41    fi=compfm(x1)
      return
42    fi=comprm(x1,x2,x3)
      return
43    fi=comptm(x1)
      return
44    fi=cratio(x1)
      return
45    fi=ebind(x1)
      return
46    fi=ebr1(x1)
      return
47    fi=ededx(x1)
      return
48    fi=eiitm(x1,x2)
      return
49    fi=esig(x1)
      return
50    fi=fcoulc(x1)
      return
51    fi=gbr1(x1)
      return
52    fi=gbr2(x1)
      return
53    fi=gmfp(x1)
      return
54    fi=pairdr(x1,x2)
      return
55    fi=pairfr(x1)
      return
56    fi=pairdz(x1,x2,x3)
      return
57    fi=pairfz(x1)
      return
58    fi=pairrm(x1,x2,x3)
      return
59    fi=pairrr(x1,x2,x3)
      return
60    fi=pairrz(x1,x2,x3,x4)
      return
61    fi=pairte(x1)
      return
62    fi=pairtm(x1)
      return
63    fi=pairtr(x1)
      return
64    fi=pairtu(x1)
      return
65    fi=pairtz(x1,x2)
      return
66    fi=pbr1(x1)
      return
67    fi=pbr2(x1)
      return
68    fi=pdedx(x1)
      return
69    fi=phottz(x1,x2)
      return
70    fi=photte(x1)
      return
71    fi=psig(x1)
      return
72    fi=spione(x1,x2)
      return
73    fi=spionp(x1,x2)
      return
74    fi=sptote(x1,x2,x3)
      return
75    fi=sptotp(x1,x2,x3)
      return
76    fi=tmxb(x1)
      return
77    fi=tmxs(x1)
      return
78    fi=tmxde2(x1)
      return
79    fi=xsif(x1)
      return
      end
      double precision function gbr1(e)
      implicit none
      double precision pair, pairtu, e, comptm, photte
      save
      pair=pairtu(e)
      gbr1=pair/(pair+comptm(e)+photte(e))
      return
      end
      double precision function gbr2(e)
      implicit none
      double precision prco, pairtu, e, comptm, photte
      save
      prco=pairtu(e)+comptm(e)
      gbr2=prco/(prco+photte(e))
      return
      end
      subroutine gfuns(e,v)
      implicit none
      
      double precision v(4)
      save
      pair=pairtu(e)
      comp=comptm(e)
      phot=photte(e)
      cohr=cohetm(e)
      tsansc=pair+comp+phot
      gmfp=1.0/tsansc
      v(1)=gmfp
      v(2)=pair*gmfp
      v(3)=(pair+comp)*gmfp
      v(4)=tsansc/(tsansc+cohr)
      return
      end
      double precision function gmfp(e)
      implicit none
      double precision pairtu, e, comptm, photte
      save
      gmfp=1.0/(pairtu(e)+comptm(e)+photte(e))
      return
      end
      subroutine hplt1(ei,el,eh,icap,ntimes,nbins,nh,idf,idsig, irsig,it,&
      sig)
      implicit none
      integer ibin, irsig, itsig, idf, nbins, ntimes, i, j, idsig, ic
      integer nh(200)
      character*4 icap(12)
      include 'pegscommons/funcs.f'
      include 'pegscommons/funcsc.f'
      character*4 l(100),cm,cr,cd,cbl
      integer ipnts
      data l/100*' '/,cm/'m'/,cr/'r'/,cd/'d'/,cbl/' '/,ipnts/10/
      save
      ibin(y)=max0(1,min0(100,idint(y/amax*100.)+1))
      rtot=fi(irsig,ei,el,eh,0.d0)
      ttot=fi(itsig,ei,0.d0,0.d0,0.d0)
      dfh=fi(idf,eh,0.d0,0.d0,0.d0)
      dfl=fi(idf,el,0.d0,0.d0,0.d0)
      deldf=(dfh-dfl)/nbins
      dnorm=rtot/(deldf*ntimes)
      amax=0.0
      eli=el
      do i=1,nbins
        ehi=fi(idf+1,dfl+deldf*i,0.d0,0.d0,0.d0)
        amax=dmax1(amax,nh(i)*dnorm,fi(irsig,ei,eli,ehi,0.d0)/deldf)
        do j=1,ipnts
        end do
        eli=ehi
      end do
      eli=el
      do i=1,nbins
        ehi=fi(idf+1,dfl+deldf*i,0.d0,0.d0,0.d0)
        v=nh(i)*dnorm
        ic=ibin(v)
        l(ic)=cm
        write(  526,130) eli,v,l
130     format(1X,1P,2E15.5,' I',100A1)
        l(ic)=cbl
        v=fi(irsig,ei,eli,ehi,0.d0)/deldf
        ic=ibin(v)
        l(ic)=cr
        write(  526,140) eli,v,l
140     format(1X,1P,2E15.5,' i',100A1)
        l(ic)=cbl
        do j=1,ipnts

          v=fi(idsig,ei,eint,0.d0,0.d0)/fi(idf+2,eint,0.d0,0.d0,0.d0)
          ic=ibin(v)
          l(ic)=cd
          write(  526,150) eint,v,l
150       format(1X,1P,2E15.5,' i',100A1)
          l(ic)=cbl
        end do
        eli=ehi
      end do
      return
      end
      integer function ifunt(name)
      implicit none
      integer if, j
      character*4 name(6)
      include 'pegscommons/funcs.f'
      include 'pegscommons/funcsc.f'
      save
      do 100 if=1,nfuns
        do j=1,6
          if (name(j).ne.fname(j,if)) go to 100
        end do
        ifunt=if
        return
100   continue
      ifunt=-1
      write(  526,110) name
110   format(' FUNC=',6A1,' not matched')
      return
      end
      subroutine lay
      implicit none
     
      include 'include/egs5_h.f'
      include 'pegscommons/bremp2.f'
      include 'pegscommons/mixdat.f'
      include 'pegscommons/mxdatc.f'
      include 'pegscommons/cohcom.f'
      include 'pegscommons/rslts.f'
      include 'pegscommons/thres2.f'
      include 'pegscommons/molvar.f'
      include 'pegscommons/eimpact.f'
      include 'pegscommons/epstar.f'
      include 'pegscommons/phpair.f'
      include 'pegscommons/bcom.f'
      include 'pegscommons/cpcom.f'
      include 'pegscommons/sfcom.f'
      save
100   FORMAT(1X,14I5)
110   FORMAT(1X,1P,5E14.5)
      ip=507
      iuecho=  526
      write(iuecho,120)
120   FORMAT(' $ECHO WRITE:MEDIUM,IDSTRN')
      write(ip,130) medium,idstrn
      write(iuecho,130) medium,idstrn
130   FORMAT(' MEDIUM=',24A1,',STERNCID=',24A1)
      if (gasp.ne.0.0) then
        write(iuecho,140)
140     FORMAT(' $ECHO WRITE:MTYP,RHO,NE,GASP, IUNRST,EPSTFL,IAPRIM')
        write(ip,150) mtyp,rho,ne,gasp, iunrst,epstfl,iaprim
        write(iuecho,150) mtyp,rho,ne,gasp, iunrst,epstfl,iaprim

      else
        write(iuecho,160)
160     FORMAT(' $ECHO WRITE:MTYP,RHO,NE,IUNRST,EPSTFL,IAPRIM')
        write(ip,170) mtyp,rho,ne,iunrst,epstfl,iaprim
        write(iuecho,170) mtyp,rho,ne,iunrst,epstfl,iaprim

      end if
      do ie=1,ne
        write(iuecho,180)
180     FORMAT(' $ECHO WRITE:ASYM(IE),Z(IE),WA(IE),PZ(IE),RHOZ(IE)')
        write(ip,190) asym(ie),z(ie),wa(ie),pz(ie),rhoz(ie)
        write(iuecho,190) asym(ie),z(ie),wa(ie),pz(ie),rhoz(ie)

      end do
      write(iuecho,200)
200   FORMAT(' $ECHO WRITE:RLC,AE,AP,UE,UP')
      write(ip,110) rlc,ae,ap,ue,up
      write(iuecho,110) rlc,ae,ap,ue,up
      nsge=0
      nseke=0
      nleke=0
      ncmfp=0
      nrange=0
      nge=ngl
      neke=nel
      write(iuecho,210)


      write(ip,110) (dl1(i),dl2(i),dl3(i),dl4(i),dl5(i),dl6(i),i=1,6)

      write(iuecho,230)
230   FORMAT(' $ECHO WRITE:DELCM,(ALPHI(I),BPAR(I),DELPOS(I),I=1,2)')
      write(ip,110) delcm,(alphi(i),bpar(i),delpos(i),i=1,2)
      write(iuecho,110) delcm,(alphi(i),bpar(i),delpos(i),i=1,2)
      write(iuecho,240)
240   FORMAT(' $ECHO WRITE:XR0,TEFF0,BLCC,XCC')
      write(ip,110) xr0,teff0,blcc,xcc
      write(iuecho,110) xr0,teff0,blcc,xcc
      write(iuecho,250)
250   FORMAT(' $Echo write:BXE,AXE')
      write(ip,110) bxe,axe
      write(iuecho,110) bxe,axe
      write(iuecho,260)

      write(ip,110) ((bfe(i,ifun),afe(i,ifun),ifun=1,15),i=1,neke)
      write(iuecho,110) ((bfe(i,ifun),afe(i,ifun),ifun=1,15),i=1,neke)
      write(iuecho,270)
270   FORMAT(' $Echo write:EBINDA,BXG,AXG')
      write(ip,110) ebinda,bxg,axg
      write(iuecho,110) ebinda,bxg,axg
      write(iuecho,280)

      write(ip,110) ((bfg(i,ifun),afg(i,ifun),ifun=1,3),i=1,nge)
      write(iuecho,110) ((bfg(i,ifun),afg(i,ifun),ifun=1,3),i=1,nge)
      if (irayl.ne.0) then
        write(iuecho,290)
290     FORMAT(' $Echo write:ngr')
        write(ip,100) ngr
        write(iuecho,100) ngr
        write(iuecho,300)
300     FORMAT(' $Echo write:bxr,axr')
        write(ip,110) bxr,axr
        write(iuecho,110) bxr,axr
        write(iuecho,310)
310     FORMAT(' $Echo write:(bfr(i),afr(i),i=1,ngr)')
        write(ip,110) (bfr(i),afr(i),i=1,ngr)
        write(iuecho,110) (bfr(i),afr(i),i=1,ngr)
        write(iuecho,320)
320     FORMAT(' $Echo write:(bfg(i,4),afg(i,4),i=1,nge)')
        write(ip,110) (bfg(i,4),afg(i,4),i=1,nge)
        write(iuecho,110) (bfg(i,4),afg(i,4),i=1,nge)
      end if
      if (incoh.eq.1) then
        write(iuecho,330)
330     FORMAT(' $Echo write:ngs')
        write(ip,100) ngs
        write(iuecho,100) ngs
        write(iuecho,340)
340     FORMAT(' $Echo write:bxs,axs')
        write(ip,110) bxs,axs
        write(iuecho,110) bxs,axs
        write(iuecho,350)
350     FORMAT(' $Echo write:(bfs(i),afs(i),i=1,ngs)')
        write(ip,110) (bfs(i),afs(i),i=1,ngs)
        write(iuecho,110) (bfs(i),afs(i),i=1,ngs)
      end if
      if (icprof.eq.1.or.icprof.eq.2) then
        write(iuecho,360)
360     FORMAT(' $Echo write:ngc')
        write(ip,100) ngc
        write(iuecho,100) ngc
        write(iuecho,370)
370     FORMAT(' $Echo write:bxc,axc,cpimev')
        write(ip,110) bxc,axc,cpimev
        write(iuecho,110) bxc,axc,cpimev
        write(iuecho,390)
390     FORMAT(' $Echo write:(bfc(i),afc(i),i=1,ngc)')
        write(ip,110) (bfc(i),afc(i),i=1,ngc)
        write(iuecho,110) (bfc(i),afc(i),i=1,ngc)
      end if
      if (icprof.eq.3.or.icprof.eq.4) then
        write(iuecho,400)
400     FORMAT(' $Echo write:mxshel,ngcs')
        write(ip,100) mxshel,ngcs
        write(iuecho,100) mxshel,ngcs
        write(iuecho,410)
410     FORMAT(' $Echo write:(elecno(ishell),ishell=1,mxshel)')
        write(ip,110) (elecno(ishell),ishell=1,mxshel)
        write(iuecho,110) (elecno(ishell),ishell=1,mxshel)
        write(iuecho,420)
420     FORMAT(' $Echo write:(capio(ishell),ishell=1,mxshel)')
        write(ip,110) (capio(ishell),ishell=1,mxshel)
        write(iuecho,110) (capio(ishell),ishell=1,mxshel)
        write(iuecho,430)
430     FORMAT(' $Echo write:bxcs,axcs')
        write(ip,110) bxcs,axcs
        write(iuecho,110) bxcs,axcs
        write(iuecho,440)

        write(ip,110) ((bfcs(i,is),afcs(i,is),is=1,mxshel),i=1,ngcs)
        write(iuecho,110) ((bfcs(i,is),afcs(i,is),is=1,mxshel),i=1,ngcs)
      end if
      if (impact.ge.1) then
        write(iuecho,450)
450     FORMAT(' $Echo write:ne')
        write(ip,100) ne
        write(iuecho,100) ne
        write(iuecho,460)
460     FORMAT(' $Echo write:neii')
        write(ip,100) neii
        write(iuecho,100) neii
        write(iuecho,470)
470     FORMAT(' $Echo write:bxeii,axeii')
        write(ip,110) bxeii,axeii
        write(iuecho,110) bxeii,axeii
        write(iuecho,480)

        write(ip,110) ((bfeii(i,ifun),afeii(i,ifun),ifun=1,ne),i=1,neii)

      end if
      return
      end
      subroutine MIX
      implicit none
      integer i, IZZ
      
      include 'include/egs5_h.f'
      include 'pegscommons/mimsd.f'
      include 'pegscommons/mixdat.f'
      include 'pegscommons/mxdatc.f'
      include 'pegscommons/molvar.f'
      include 'pegscommons/pmcons.f'
      include 'pegscommons/dercon.f'
      include 'pegscommons/radlen.f'
      include 'pegscommons/mscom.f'
      dimension XSI(MXEL),ZZX(MXEL),FZC(MXEL),FCOUL(MXEL),ZZ(MXEL)
      save
      write(  526,100)
100   format(/' In subroutine mix: '/)
      if (GASP.eq.0.0) then
        write(  526,110) NE,RHO

      else
        write(  526,120) NE,RHO,GASP
      end if
      write(  526,130) (i,Z(i),WA(i),PZ(i),RHOZ(i),i=1,NE)

      if (GASP.ne.0.0) then
        RHO=GASP*RHO
      end if
      AL183 = DLOG(A183)
      TPZ=0.0
      WM=0.0
      ZC=0.0
      ZT=0.0
      ZB=0.0
      ZF=0.0
      ZS=0.0
      ZE=0.0
      ZX=0.0
      ZAB=0.0
      do i=1,NE
        TPZ = TPZ + PZ(i)
        WM = WM + PZ(i)*WA(i)
        ZC = ZC + PZ(i)*Z(i)
        FZC(i) =(FSC*Z(i))**2
        FCOUL(i) = FCOULC(Z(i))
        XSI(i) = XSIF (Z(i))
        ZZX(i) = PZ(i)*Z(i)*(Z(i)+XSI(i))
        if (Z(i).le.4.0) then
          IZZ=Z(i)
          ZAB=ZAB+ZZX(i)*ALRAD(IZZ)
        else
          ZAB=ZAB+ZZX(i)*(AL183+DLOG(Z(i)**(-1./3.)))
        end if
        ZT = ZT + ZZX(i)
        ZB = ZB + ZZX(i)*dlog(Z(i)**(-1./3.))
        ZF = ZF + ZZX(i)*FCOUL(i)
        ZZ(i) = PZ(i)*Z(i)*(Z(i)+fudgeMS)
        ZS = ZS + ZZ(i)
        ZE = ZE + ZZ(i)*((-2./3.)*DLOG(Z(i)))
        ZX = ZX + ZZ(i)*DLOG(1.d0+3.34*FZC(i))
      end do
      EZ = ZC/TPZ
      ZA = AL183*ZT
      ZG = ZB/ZT
      ZP = ZB/ZA
      ZV = (ZB-ZF)/ZT
      ZU = (ZB-ZF)/ZA
      EDEN=AN*RHO/WM*ZC
      RLC = 1./( (AN*RHO/WM)*4.0*FSC*R0**2*(ZAB-ZF) )

      V3120=EDEN
      write(  526,150) V3120
150   FORMAT(' Eden=',1P,G15.7)
      BLCC= A6680*RHO*ZS*DEXP(ZE/ZS)*RLC / (WM*DEXP(ZX/ZS))
      TEFF0 = ( DEXP(BMIN)/BMIN )/BLCC
      XCC= (A22P9/RADDEG) * DSQRT( ZS*RHO*RLC/WM )
      XR0 = XCC*DSQRT(TEFF0*BMIN)
      write(  526,160) BLCC,XCC,TEFF0,XR0
160   format(' BLCC,XCC,TEFF0,XR0=',1P,4E14.5)
      return
      end
      subroutine molier
      implicit none
      save
     
      double precision BLCMIN, B, BLCA, BA, B1, PTOT, P, Q, PPP, PP
      include 'pegscommons/mimsd.f'

      double precision TH(29),DTH(29),F0(29),F1(29),F2(29),BOLD,BLC
      
      save
      write(  526,100) (TH(i),DTH(i),F0(i),F1(i),F2(i),i=1,29)
100   format(' Bethe table used for input'/(1X,0P,2F10.2,1P,3E18.5))
      IS=1
        go to 120
110     IS=IS+1
120     if(IS-(MSTEPS-1).gt.0) go to 160
        J=FSTEP(IS)
          go to 140
130       J=J+1
140       IF(J-(FSTEP(IS+1)-1).gt.0) go to 150
          MSMAP(J)=IS
        go to 130
150     continue
      go to 110
160   continue
      MSMAP(JRMAX)=MSTEPS
      BLCMIN = BMIN - DLOG(BMIN)
      do 260 IS=1,MSTEPS
        BLC=BLCMIN+DLOG(FSTEP(IS))
        B=BLC+DLOG(BLC)
170     continue
          BOLD=B
          B=BOLD - (BOLD-DLOG(BOLD)-BLC)/(1.0-1.0/BOLD)
          if (dabs((B-BOLD)/BOLD) .lt. 1.D-5) go to 180
        go to 170
180     continue
        BLCA(IS)=BLC
        BA(IS)=B
        FSQR(IS)=DSQRT(FSTEP(IS)*B/BMIN)
        B1=1.0/B
        PTOT=0.0
        do i=1,29
          P(i,IS)=TH(i)*DTH(i)*(F0(i)+B1*(F1(i)+B1*F2(i)))
          PTOT=PTOT+P(i,IS)
        end do
        do i=1,29
          P(i,IS)=P(i,IS)/PTOT
        end do
        do i=1,29
          Q(i,IS)=P(i,IS)
        end do
        i=29
190     continue
          l=1
200       if(Q(i,IS).ge.0.001.or.i.le.l) go to 210
            Q(i,IS)=Q(i,IS)+Q(i-l,IS)
            Q(i-l,IS)=0.0
            l=l+1
          go to 200
210       continue
          i=i-l
          if (i.le.0) go to 220
        go to 190
220     continue
        PPP=0.5
        PP=0.5
        do JLR=1,10
          ITOT=0
          do i=1,29
            IP1(i,IS)=Q(i,IS)*1000.0+PP
            ITOT=ITOT+IP1(i,IS)
          end do
          IDIF=ITOT-1000
          if (IDIF.eq.0) go to 260
          PPP=PPP*0.5
          if (IDIF.lt.0) then
            PP=PP+PPP
          else
            PP=PP-PPP
          end if
        end do
        do i=1,29
          IP2(i)=1
        end do
        n=29
230     continue
          n=n-1
          IFLG=0
          do j=1,n
            I01=IP2(j)
            I02=IP2(j+1)
            if (IP1(I01,IS).lt.IP1(I02,IS))  then
              ISWP=IP2(j)
              IP2(j)=IP2(j+1)
              IP2(j+1)=ISWP
              IFLG=1
            end if
          end do
          if (IFLG.eq.0) go to 240
        go to 230
240     continue
        write(  526,250) ITOT
250     FORMAT(' Rounding failed, itot has',I6,' entries')
        if (IDIF.lt.0) then
          IDA=-IDIF
          INC=1
        else
          IDA=IDIF
          INC=-1
        end if
        do i=1,IDA
          I01=IP2(i)
          IP1(I01,IS)=IP1(I01,IS)+INC
        end do
260   continue
      MXV1=0
      do i=1,29
        IALL(i)=IP1(i,1)
        do is=2,MSTEPS
          IALL(i)=MIN0(IALL(i),IP1(i,is))
        end do
        MXV1=MXV1+IALL(I)
      end do
      MXV2=1000-MXV1
      ii=0
      do i=1,29
        j=1
          go to 280
270       j=j+1
280       if(j-(IALL(i)).gt.0) go to 290
          ii=ii+1
          VERT1(ii)=TH(i)
        go to 270
290     continue
      end do
      do is=1,MSTEPS
        ii=0
        do i=1,29
          IXTR(i,is)=IP1(i,is)-IALL(i)
          j=1
            go to 310
300         j=j+1
310         if(j-(IXTR(i,is)).gt.0) go to 320
            ii=ii+1
            VERT2(ii,is)=TH(i)
          go to 300
320       continue
        end do
      end do
      write(  526,330) BMIN,MSTEPS,JRMAX,MXV1,MXV2
330   format(' BMIN,MSTEPS,JRMAX,MXV1,MXV2=', F11.5,4I8)
      ISU=0
340   continue
        ISL=ISU+1
        ISU=MIN0(ISL+9,MSTEPS)
        write(  526,350) ISL,ISU
350     format('  Data for steps ',I3,' to ',I3)
        write(  526,360) (IS,IS=ISL,ISU)
360     format(11X,'ISTEP',I6,9I11)
        write(  526,370) (FSTEP(IS),IS=ISL,ISU)
370     format(11X,'FSTEP',10F11.0)
        write(  526,380) (FSQR(IS),IS=ISL,ISU)
380     format(11X,'FSQR ',10F11.5)
        write(  526,390) (BLCA(IS),IS=ISL,ISU)
390     format(11X,'BLC  ',10F11.5)
        write(  526,400) (BA (IS),IS=ISL,ISU)
400     format(11X,'B    ',10F11.5)
        write(  526,410)
410     format('0I  TH IALL')
        do i=1,29
          if ((i.eq.11).or.(i.eq.23)) then
            write(  526,420)
420         format('1I  TH IALL')
          end if
          write(  526,430) i,TH(i),IALL(i),(P(i,IS),IS=ISL,ISU)
430       format(1X,I2,F5.2,I4,' PR ',10F11.8)
          write(  526,440) (Q(I,IS),IS=ISL,ISU)
440       format(11X,'  Q  ',10F11.8)
          write(  526,450) (IP1(I,IS),IS=ISL,ISU)
450       format(11X,' IP1 ',I7,9I11)
          write(  526,460) (IXTR(I,IS),IS=ISL,ISU)
460       format(11X,'EXTRA',I7,9I11)
        end do
        if (ISU.ge.MSTEPS) go to 470
      go to 340
470   continue
480   format(1X,14I5)
490   format(1X,14F5.2)
      iuecho=  526
      ipun=507
      write(iuecho,500)
500   format(' $echo write:')
      write(ipun,510)
      write(iuecho,510)
510   format(' Material independent multiple scattering data')
      write(iuecho,520)
520   format(' $echo write:JRMAX,MSTEPS,MXV1,MXV2')
      write(ipun,480) JRMAX,MSTEPS,MXV1,MXV2
      write(iuecho,480) JRMAX,MSTEPS,MXV1,MXV2
      write(iuecho,530)
530   format(' $echo write:(FSTEP(i),FSQR(i),i=1,MSTEPS)')
      write(ipun,540) (FSTEP(i),FSQR(i),i=1,MSTEPS)
      write(iuecho,540) (FSTEP(i),FSQR(i),i=1,MSTEPS)
540   format((1X,4(F5.0,F11.6)))
      write(iuecho,550)
550   format(' $echo write:(MSMAP(I),I=1,JRMAX)')
      write(ipun,480) (MSMAP(i),i=1,JRMAX)
      write(iuecho,480) (MSMAP(i),i=1,JRMAX)
      write(iuecho,560)
560   format(' $echo write:(VERT1(I),I=1,MXV1)')
      write(ipun,490) (VERT1(i),i=1,MXV1)
      write(iuecho,490) (VERT1(i),i=1,MXV1)
      do MST=1,MSTEPS
        write(iuecho,570) MST
570     format(' MST=',I5)
        write(iuecho,580)
580     format(' $echo write:(VERT2(I,MST),I=1,MXV2)')
        write(ipun,490) (VERT2(i,MST),i=1,MXV2)
        write(iuecho,490) (VERT2(i,MST),i=1,MXV2)
      end do
      return
      end
      double precision function pairdr(ka,e)
      implicit none
      integer LS
      double precision pairfr, e
      include 'pegscommons/lpairr.f'
      double precision ka
      save
      K=ka
      if (K.lt.50.) then
        LE=1
        LS=0
      else
        LE=2
        LS=3
      end if
      LA=LS+1
      LC=LS+3
      PAIRDR=PAIRFR(e)
      return
      end
      double precision function pairdz(Z,KA,E)
      implicit none
      double precision KA, Z, xsif, dlog, fcoulc, pairfz, E
      include 'pegscommons/pmcons.f'
      include 'pegscommons/dercon.f'
      include 'pegscommons/molvar.f'
      include 'pegscommons/lpairz.f'
      save
      k=ka
      DELC=136.*z**(-1./3.)*RM/K
      CONST=(AN*RHO/WM)*R0**2*FSC*z*(Z+XSIF(Z))*RLC/K**3
      XLNZ=4./3.*dlog(Z)
      if (K.ge.50) XLNZ=XLNZ+4.*FCOULC(Z)
      DELTAM=dexp((21.12-XLNZ)/4.184)-0.952
      PAIRDZ=PAIRFZ(E)
      return
      end
      double precision function pairfr(E)
      implicit none
      double precision EPS, E, DEL, DELTA, A, CC
      include 'pegscommons/bremp2.f'
      include 'pegscommons/dbrpr.f'
      include 'pegscommons/lpairr.f'
      save
      EPS=E/K
      DEL=1./(K*EPS*(1.-EPS))
      if (DEL.gt.DELPOS(LE)) then
        pairfr=0.0
      else
        DELTA=DELCM*DEL
        if (DELTA.le.1.) then
          A=DL1(LA)+DELTA*(DL2(LA)+DELTA*DL3(LA))
          CC=DL1(LC)+DELTA*(DL2(LC)+DELTA*DL3(LC))
        else
          A=DL4(LA)+DL5(LA)*DLOG(DELTA+DL6(LA))
          CC=DL4(LC)+DL5(LC)*DLOG(DELTA+DL6(LC))
        end if
        pairfr=(ALFP1(LE)*CC+ALFP2(LE)*12.*(E/K-0.5)**2*A)/K
      end if
      return
      end
      double precision function pairfz(E)
      implicit none
      double precision EPS, E, ONEEPS, DELTA, SB1, SB2, DLOG, EPLUS
      include 'pegscommons/lpairz.f'
      save
      EPS=E/K
      ONEEPS=1.-EPS
      if (ONEEPS.eq.0.0) then
        ONEEPS=1.18D-38
      end if
      DELTA=DELC/(EPS*ONEEPS)
      if (DELTA.ge.DELTAM) then
        pairfz=0.0
      else
        if (DELTA.le.1.) then
          SB1=20.867+DELTA*(-3.242+DELTA*0.625)-XLNZ
          SB2=20.209+DELTA*(-1.930+DELTA*(-0.086))-XLNZ
        else
          SB1=21.12-4.184*DLOG(DELTA+0.952)-XLNZ
          SB2=SB1
        end if
        EPLUS=K-E
        pairfz=CONST*((E**2+EPLUS**2)*SB1+0.666667*E*EPLUS*SB2 )
      end if
      return
      end
      double precision function pairrm(K,E1,E2)
      implicit none
      integer i
      double precision pairrz, E1, E2
      double precision K
      include 'include/egs5_h.f'
      include 'pegscommons/mixdat.f'
      include 'pegscommons/mxdatc.f'
      save
      pairrm=0.
      do i=1,NE
        pairrm=pairrm+PZ(i)*pairrz(Z(i),K,E1,E2)
      end do
      return
      end
      double precision function pairrr(K,E1,E2)
      implicit none
      double precision DUMMY, pairdr, E1, QD, E2
      double precision K
      external pairfr
      save
      DUMMY=pairdr(K,E1)
      PAIRRR=QD(PAIRFR,E1,E2,'PAIRFR')
      return
      end
      double precision function pairrz(Z,K,E1,E2)
      implicit none
      double precision DUMMY, pairdz, Z, E1, QD, E2
      double precision K
      external pairfz
      save
      DUMMY=pairdz(Z,K,E1)
      pairrz=QD(PAIRFZ,E1,E2,'PAIRFZ')
      return
      end
      double precision function pairte(K)
      implicit none
      integer i
      double precision pairtz
      double precision K
      include 'include/egs5_h.f'
      include 'pegscommons/dercon.f'
      include 'pegscommons/mixdat.f'
      include 'pegscommons/mxdatc.f'
      save
      pairte=0.0
      if (K.le.2.0*RM) return
      do i=1,NE
        pairte=pairte+PZ(i)*pairtz(Z(i),K)
      end do
      return
      end
      double precision function pairtm(K0)
      implicit none
      double precision pairrm
      double precision K0
      include 'pegscommons/dercon.f'
      save
      if (K0.le.2.*RM) then
        pairtm=0.0
      else
        pairtm=pairrm(K0,RM,K0-RM)
      end if
      return
      end
      double precision function pairtr(K0)
      implicit none
      double precision pairrr
      double precision K0
      include 'pegscommons/dercon.f'
      save
      if (K0.le.2.*RM) then
        pairtr=0.0
      else
        pairtr=pairrr(K0,RM,K0-RM)
      end if
      return
      end
      double precision function pairtu(K)
      implicit none
      double precision pairte, pairtm
      double precision K
      save
      if (K.lt.50) then
        pairtu=pairte(K)
      else
        pairtu=pairtm(K)
      end if
      return
      end
      double precision function pairtz(Z,K)
      implicit none
      integer IZ
      double precision PCON, Z, AINTP
      double precision K
      include 'pegscommons/dercon.f'
      include 'pegscommons/phpair.f'
      include 'pegscommons/pmcons.f'
      include 'pegscommons/molvar.f'
      save
      if (K.le.RMT2) then
        pairtz=0.0
      else
        PCON=1.D-24*(AN*RHO/WM)*RLC
        IZ=Z
        pairtz=PCON*AINTP(K,PRE,17,PRD(1,IZ),1,.TRUE.,.FALSE.)
      end if
      return
      end
      double precision function pbr1(E)
      implicit none
      double precision BREM, BREMTM, E, BHABTM, ANIHTM
      save
      BREM=BREMTM(E)
      pbr1=BREM/(BREM+BHABTM(E)+ANIHTM(E))
      return
      end
      double precision function pbr2(E)
      implicit none
      double precision BRBH, BREMTM, E, BHABTM, ANIHTM
      save
      BRBH=BREMTM(E)+BHABTM(E)
      pbr2=BRBH/(BRBH+ANIHTM(E))
      return
      end
      double precision function pdedx(E)
      implicit none
      double precision SPTOTP, E
      include 'pegscommons/thres2.f'
      save
      pdedx=SPTOTP(E,AE,AP)
      return
      end




      subroutine pegs5
      USE PEGS_DCSSTR_MOD
      USE egs5_media_mod
      implicit none
     
      include 'include/egs5_h.f'
      include 'pegscommons/bremp2.f'
      include 'pegscommons/dbrpr.f'
      include 'pegscommons/dercon.f'
      include 'pegscommons/elemtb.f'
      include 'pegscommons/elmtbc.f'
      include 'pegscommons/funcs.f'
      include 'pegscommons/funcsc.f'
      include 'pegscommons/lspion.f'
      include 'pegscommons/mimsd.f'
      include 'pegscommons/mixdat.f'
      include 'pegscommons/mxdatc.f'
      include 'pegscommons/molvar.f'
      include 'pegscommons/phpair.f'
      include 'pegscommons/pmcons.f'
      include 'pegscommons/pwlfin.f'
      include 'pegscommons/cohcom.f'
      include 'pegscommons/rslts.f'
      include 'pegscommons/thres2.f'
      include 'pegscommons/epstar.f'
      include 'pegscommons/bcom.f'
      include 'pegscommons/cpcom.f'
      include 'pegscommons/eimpact.f'
      include 'pegscommons/sfcom.f'
      include 'pegscommons/mscom.f'
      include 'pegscommons/legacy.f'
      common /paran/ icfn(100), ilfn(100), chfn(100)
      character chfn*200
      integer icfn, ilfn
      integer imategs,ilscat  ! temporary parameters for checking EGS options
      double precision XP(4),WASAV(MXEPERMED)
      logical MEDSET,ENGSET
      character*4 OPTION(4,14),OPT(4),BLKW,NAME(6)
      character*4 NAMESB(12),IDFNAM(6)
      integer NH(200)
      
      intrinsic DLOG,DEXP
     
      data NOPT/15/,BLKW/' '/
      data MEDSET/.FALSE./,ENGSET/.FALSE./
      data NPTS/50/,IDF/6/
     
      namelist/BCOMDT/BCOMP,PBC,NPBC
      namelist/ISCADT/SCATF,XSVAL
      namelist/SCPRDT/MXRAW,MXSHEL,ELECNI,NSHELL,CAPIN,SCPROF,QCAP

      integer iegsemi, iegsout
      common /egsemi/ iegsemi, iegsout
      integer npeme
      common /mpiegs/ npeme

      save


      

      open(UNIT=525,FILE=chfn(23)(1:ilfn(23))//'.inp',STATUS='old')

      if( iegsout .eq. 2 .and. npeme .eq. 0 ) then
         
      else
         open(UNIT=526,form='formatted',status='scratch')
      end if

 

      if( iegsout .eq. 2 .and. npeme .eq. 0 ) then
       open(UNIT=510,FILE=chfn(23)(1:ilfn(23))//'.err',STATUS='unknown')

      else
         open(UNIT=510,form='formatted',status='scratch')
         open(UNIT=521,form='formatted',status='scratch')
      end if
      if( iegsout .eq. 2 .and. npeme .eq. 0 ) then
       open(UNIT=515,FILE=chfn(23)(1:ilfn(23))//'.scp',STATUS='unknown')
       open(UNIT=518,FILE=chfn(23)(1:ilfn(23))//'.iff',STATUS='unknown')
       open(UNIT=519,FILE=chfn(23)(1:ilfn(23))//'.ics',STATUS='unknown')
      else
         open(UNIT=515,form='formatted',status='scratch')
         open(UNIT=518,form='formatted',status='scratch')
         open(UNIT=519,form='formatted',status='scratch')
      end if



      write(  526,100)

      call pmdcon
      write(  526,110)
110   FORMAT(/' This version reads units 8 and 9 in free format'/)
      read(509,*) XVAL, ((AFAC(in,iz),in=1,100),iz=1,100)
      do i=1,100
        XVAL(i)=XVAL(i)**2.
      end do
      ISSBS=ISSB
      AFACTS=AFACT
      CBARS=CBAR
      SKS=SK
      X0S=X0
      X1S=X1
      IEVS=IEV
      ibounds=IBOUND
      incohs=INCOH
      icprofs=ICPROF
      gasps=GASP
      irayls=IRAYL
      impacts=IMPACT
      iunrsts=IUNRST
      fudgeMSs=fudgeMS
      efracHs=efracH
      efracLs=efracL
      nleg0s=nleg0
      epstfls=EPSTFL
      iepsts=IEPST
      iaprims=IAPRIM
      iaprfls=IAPRFL
1040  continue
      do i=1,MXEPERMED
        WASAV(i)=WA(i)
        WA(i)=0.
      end do
      RHOSAV=RHO
      RHO=0.0
      read(  525,130,end=1060) OPT
130   FORMAT(4A1)
      write(  526,140) OPT

      read(  525,INP,end=8990)
      if (RHO.eq.0) then
        do ICH=1,4
          if (OPT(ICH).ne.OPTION(ICH,1)) then
            RHO=RHOSAV
            go to 150
          end if
        end do
150     continue
      end if
      do 160 IOPT=1,NOPT
        do ICH=1,4
          if (OPT(ICH).ne.OPTION(ICH,IOPT)) go to 160
        end do
        go to 1120
160   continue
      write(  526,1130)
      write(    *,1130)
1130  format(' Option not found, job aborted.')
      close(526)
      stop 16
1120  if (IOPT.gt.3) then
        RHO=RHOSAV
        do i=1,MXEPERMED
          WA(i)=WASAV(i)
        end do
      else
        do i=1,4
          MTYP(i)=OPT(i)
        end do
      end if

1160  NE=1
      PZ(1)=1
      IMIXT=0
      go to 1290
1170  if (NE.le.1) then
        write(  526,1300) NE
        write(  *,1300) NE
1300    format(//,' NE=',I6,' is improperly defined for a mixture.')
        close(526)
        stop
      end if
      IMIXT=1
      go to 1290
1180  if (NE.le.1) then
        write(  526,1310) NE
        write(  *,1310) NE
1310    format(//,' NE=',I6,' is improperly defined for a compound.')
        close(526)
        stop
      end if
      IMIXT=0
1290  read(  525,1320) MEDIUM,IDSTRN
1320  FORMAT(24A1,6X,24A1)
      read(  525,1330) (ASYM(i),i=1,NE)
      if (IDSTRN(1).eq.BLKW) then
        do i=1,LMED
          IDSTRN(i)=MEDIUM(i)
        end do
      end if
      write(  526,1350) MEDIUM,IDSTRN,(ASYM(i),i=1,NE)

      do 16 medIdx = 1, nmed
        do ib = 1, 24
          if(MEDIUM(ib) .ne. media(ib,medIdx)) go to 16
          if(ib .eq. 24) go to 17
        end do
16    continue
17    oldK1run = .true.
      if(medIDx.le.nmed) then
        if(charD(medIdx).ne.0.d0) then
          oldK1run = .false.
        endif
      endif
      if (IUNRST.eq.1) then
        write(  526,1360)

      else if (IUNRST.eq.2) then
        write(  526,1370)

      else if (IUNRST.eq.3) then
        write(  526,1380)

      else if (IUNRST.eq.4) then
        write(  526,1390)

      else if (IUNRST.eq.5) then
        write(  526,1400)

      else if (IUNRST.eq.6) then
        write(  526,1410)

      else if (IUNRST.eq.7) then
        write(  526,1420)

      end if
      do i=1,NE
        Z(i)=ZTBL(ASYM(i))
        if (Z(i).eq.0.0) then
          write(  526,1440)
          write(  *,1440)
1440      format(' Bad atomic symbol....job aborted.')
          close(526)
          stop 16
        end if
        if (WA(I).eq.0.) then
          I01=Z(i)
          WA(i)=WATBL(I01)
        end if
        if (IMIXT.ne.0) then
          PZ(i)=RHOZ(i)/WA(I)
        else
          RHOZ(i)=PZ(i)*WA(i)
        end if
      end do
      if (NE.eq.1.and.RHO.eq.0.) then
        I01=Z(1)
        RHO=RHOTBL(I01)
      end if
      call MIX
      call SPINIT
      call DIFFER
      ! Read in the DCS and get the scattering power
      call elinit(z,pz,ne)
      MEDSET=.true.
      read(511,BCOMDT)
      rewind 511
      read(  513,ISCADT)
      rewind   513
      if (ICPROF.eq.3.or.ICPROF.eq.4) then
        write(  526,1450) ICPROF
1450    format(' Reading shellwise Compton profile as ICPROF=',I3)
        read(515,SCPRDT,ERR=8991)
      end if
      if (ICPROF.eq.-3.or.ICPROF.eq.-4) then
        write(  526,1460) ICPROF
1460    format(' Reading shellwise Compton profile as ICPROF=',I3)
        MXRAW=0
        MXSHEL=0
        call RDSCPR
        open(UNIT=531,FILE=chfn(23)(1:ilfn(23))//'job.ssl',STATUS='old')
        read(531,SCPRDT)
        close(531)
        ICPROF=ABS(ICPROF)
      end if
      if (IRAYL.eq.2) then
        do i=1,100
          read(518,*,ERR=8992) XVALFAC(i),BFAC2(i)
          write(  526,1480) XVALFAC(i),BFAC2(i)
1480      format(2G10.5)
        end do
        write(  526,1490)

        do i=1,61
          read(519,*,ERR=8993) XVALINT(i),COHEINT(i,1)
          write(  526,1510) XVALINT(i),COHEINT(i,1)
1510      format(2G10.5)
        end do
        write(  526,1520)

      end if
      write(  526,1530)
1530  format(//,' End of elem, mixt, or comp option',///)
      go to 1040
1190  if (AE.lt.0) AE=-AE*RM
      if (UE.lt.0) UE=-UE*RM
      if (AP.lt.0) AP=-AP*RM
      if (UP.lt.0) UP=-UP*RM
      TE=AE-RM
      TET2=TE*2.0
      TEM=TE/RM
      THBREM=RM+AP
      THMOLL=AE+TE
      write(  526,1540) AE,UE,AP,UP,TE,TET2,TEM,THBREM,THMOLL

      call inigrd(ae,ue,1)
      ENGSET=.true.
      call esteplim
      go to 1040
1200  call MOLIER
      go to 1040
1210  if (MEDSET.and.ENGSET) go to 1550
      write(  526,1560) MEDSET,ENGSET
      write(  *,1560) MEDSET,ENGSET
1560  format(' MEDSET,ENGSET=',2L2,',PWLF req. ignored.')
      close(526)
      stop 16
1550  EBINDA=EBIND(AP)
      write(  526,PWLFNM)
      write(  526,1570) EBINDA
1570  format(' Average K-ionization energy=',F10.6,'(MeV)')
      if(efracH .gt. 0.5d0) then
        efracH = .5d0
        write(506,1541) efracH
        write(526,1541) efracH
      else if(efracH .le. 1.d-6) then
        efracH = 5.d-2
        write(506,1542) efracH
        write(526,1542) efracH
      endif
      if(efracL .gt. 0.5d0) then
        efracL = .5d0
        write(506,1543) efracL
        write(526,1543) efracL
      else if(efracL .le. 1.d-6) then
        efracL = 5.d-2
        write(506,1544) efracL
        write(526,1544) efracL
      endif
1541  format('Warning:  EFRAC_HIGH > MAX.  setting to ',F6.4)
1542  format('Warning:  EFRAC_HIGH < MIN.  setting to ',F6.4)
1543  format('Warning:  EFRAC_LOW > MAX.  setting to ',F6.4)
1544  format('Warning:  EFRAC_LOW < MIN.  setting to ',F6.4)
      if(oldK1run) call makek1(ae,ue)
      write(  526,1580)
1580  format(' Do PWLF to electron data sets.'/)
      call PWLF1(NEL,NALE,AE,UE,THMOLL,EPE,ZTHRE,ZEPE,NIPE,ALKE,ALKEI, A&
      ,XE,BXE,150,15,AFE,BFE,EFUNS)
      if (IMPACT.ge.1) then
        write(  526,1590) IMPACT
        call PWLF1(NEII,NALE,AE,UE,THMOLL,EPE,ZTHRE,ZEPE,NIPE,ALKE,ALKEI &
        , AXEII,BXEII,150,20,AFEII,BFEII,EIIFUNS)
      end if
      write(  526,1600)
1600  format(' Do pwlf to photon data sets.'/)
      if (IBOUND.eq.1) then
        write(  526,1610)

      end if
      call PWLF1(NGL,NALG,AP,UP,RMT2,EPG,ZTHRG,ZEPG,NIPG,DLOG,DEXP,AXG,B,&
       XG,1000,4,AFG,BFG,GFUNS)
      if (IRAYL.eq.1) then
        write(  526,1620)
1620    format(//,' ***** IRAYL=1: Rayleigh data included *****',//)
        write(  526,1630)
1630    format(' Do PWLF to Rayleigh distribution.'/)
        do i=1,100
          AFAC2(i)=0.0
          do in=1,NE
            IZZ=Z(in)
            AFAC2(i)=AFAC2(i)+PZ(in)*AFAC(i,IZZ)**2
          end do
        end do
        AFFI(1)=0.0
        do i=2,97
          AX=XVAL(i-1)
          BX=XVAL(i)
          AFFI(i)=QD(AFFACT,AX,BX,'AFFACT')
        end do
        do i=2,97
          AFFI(i)=AFFI(i)+AFFI(i-1)
        end do
        do i=1,97
          AFFI(i)=AFFI(i)/AFFI(97)
        end do
        call PWLF1(NGR,NALR,0.D0,1.D0,0.D0,EPR,ZTHRR,ZEPR,NIPR,ALIN,&
        ALINI,AXR,BXR,100,1,AFR,BFR,RFUNS)
      end if
      if (IRAYL.eq.2) then
        write(  526,1690)

        write(  526,1700)
1700    format(' Do PWLF to Rayleigh distribution.'/)
        do i=1,100
          AFAC2(i)=(BFAC2(i)*DSQRT(WM))**2
        end do
        do i=1,100
          XVAL(i)=XVALFAC(i)**2
        end do
        AFFI(1)=0.0
        do i=2,97
          AX=XVAL(i-1)
          BX=XVAL(i)
          AFFI(i)=QD(AFFACT,AX,BX,'AFFACT')
        end do
        do i=2,97
          AFFI(i)=AFFI(i)+AFFI(i-1)
        end do
        do i=1,97
          AFFI(i)=AFFI(i)/AFFI(97)
        end do
        call PWLF1(NGR,NALR,0.D0,1.D0,0.D0,EPR,ZTHRR,ZEPR,NIPR,ALIN,&
        ALINI,AXR,BXR,100,1,AFR,BFR,RFUNS)
      end if
      if (IRAYL.eq.3) then
        write(  526,1760)
1760    format(//,' ***** IRAYL=3: Rayleigh data included *****',//)
        write(  526,1770)

        do i=1,100
          AFAC2(i)=0.0
          do IN=1,NE
            IZZ=Z(IN)
            AFAC2(I)=AFAC2(I)+PZ(IN)*AFAC(I,IZZ)**2
          end do
        end do
        call PWLF1(NGR,NALR,1.0D-3,1.0D2,1.0D-3,EPR,ZTHRR,ZEPR,NIPR,DLOG&
        ,DEXP, AXR,BXR,100,1,AFR,BFR,RFUNS2)
      end if
      if (INCOH.eq.1) then
        write(  526,1800)
1800    format(//,' ***** INCOH=1: S(X,Z) data included *****',//)
        write(  526,1810)
1810    format(' Do PWLF to S(X,Z)'/)
        do i=1,45
          SCSUM=0.0
          ZSUM=0.0
          do in=1,NE
            IZZ=Z(in)
            SCSUM=SCSUM+PZ(in)*SCATF(i,IZZ)
            ZSUM=ZSUM+PZ(in)*Z(in)
          end do
          SCATZ(i)=SCSUM/ZSUM
        end do
        call PWLF1(NGS,NALS,5.0D-3,80.D0,80.D0,EPSF,ZTHRS,ZEPS,NIPS,&
        DLOG,DEXP,AXS,BXS,100,1,AFS,BFS,SFUNS)
      end if
      if(ICPROF.eq.1) then
        write(  526,1830) ICPROF

        ICPROF=0
      else if(ICPROF.eq.2) then
        write(  526,1830) ICPROF
        ICPROF=0
      endif
      if (ICPROF.eq.1.or.ICPROF.eq.2) then
        CPIMEV=IEV*1.0D-6
        write(  526,1840) CPIMEV
1840    format(' CPIMEV=',E12.5)
        write(  526,1850) ICPROF

        if (ICPROF.eq.1) then
          write(  526,1860)
1860      format(' Do PWLF to Q vs int J(Q)'/)
        end if
        if (ICPROF.eq.2) then
          write(  526,1870)
1870      format(' Do PWLF to J(Q)'/)
        end if
        do i=1,31
          CPSUM=0.0
          PZSUM=0.0
          do in=1,NE
            IZZ=Z(in)
            CPSUM=CPSUM+PZ(in)*CPROF(i,IZZ)
            PZSUM=PZSUM+PZ(in)
          end do
          AVCPRF(i)=CPSUM/PZSUM
          write(  526,1900) QCAP(I),AVCPRF(I)
1900      format(' QCAP= ',1P,E9.2,' AVCPRF= ',1P,E9.2)
        end do
        if (ICPROF.eq.2) then
          call PWLF1(NGC,NALC,1.D-2,1.D2,1.D2,EPCP,ZTHRC,ZEPC,NIPC,DLOG&
          ,DEXP, AXC,BXC,2000,1,AFC,BFC,CFUNS)
        end if
        if (ICPROF.eq.1) then
          CPROFI(1)=0.0
          QCAP10(1)=0.0
          do i=2,301
            ILOC=(i-2)/10+1
            STEP=(QCAP(ILOC+1)-QCAP(ILOC))*0.1
            AX=QCAP(ILOC)+STEP*FLOAT(i-2-(ILOC-1)*10)
            BX=AX+STEP
            QCAP10(i)=BX
            CPROFI(i)=QD(CPRFIL,AX,BX,'CPRFIL')
            write(  526,1920) i,CPROFI(i)
1920        format(1X,I5,'-th interval CPROFI=',1P,E12.5)
          end do
          write(  526,1930)
1930      format(' Sum of integration of CPROF')
          do i=2,301
            CPROFI(i)=CPROFI(i)+CPROFI(i-1)
            write(  526,1950)I,CPROFI(i)
1950        format(1X,I5,'-th interval. Sum CPROFI=',1P,E12.5)
          end do
          write(  526,1960)
1960      format(' Normalized sum of integration of CPROF')
          do i=1,301
            CPROFI(i)=CPROFI(i)/CPROFI(301)
            write(  526,1980) i,CPROFI(i)
1980        format(1X,I5,'-th interval. NORM.SUM CPROFI=',1P,E12.5)
          end do
          call PWLF1(NGC,NALC,0.D0,1.D0,0.D0,EPCP,ZTHRC,ZEPC,NIPC,ALIN,&
          ALINI, AXC,BXC,2000,1,AFC,BFC,CFUNS2)
        end if
        write(  526,1990) NGC
1990    format(' Compton profile was PWLF-ED, NGC= ',I5)
      end if
      if (ICPROF.eq.3.OR.ICPROF.eq.4) then
        write(  526,2000) ICPROF

        if (ICPROF.eq.3) then
          write(  526,2010)
2010      format(' Do pwlf to Q vs int J(Q)'/)
        end if
        if (ICPROF.eq.4) then
          write(  526,2020)
2020      format(' Do PWLF to J(Q)'/)
        end if
        do i=1,MXSHEL
          ELECNO(i)=0.
          ELECNJ(i)=0.
        end do
        do i=1,MXRAW
          ITEMP=NSHELL(i)
          ELECNO(ITEMP)=ELECNO(ITEMP)+ELECNI(i)
        end do
        do i=1,MXSHEL
          write(  526,2060) i,ELECNO(i)
2060      format(' Shell ,ELECNO=',I5,E12.5)
        end do
        do i=1,MXSHEL
          CAPIO(i)=0.0
          CAPILS(i)=0.0
        end do
        do i=1,MXRAW
          if (CAPIN(i).gt.0.0) then
            CAPIL=DLOG(CAPIN(i))
            J=NSHELL(i)
            CAPILS(J)=CAPILS(J)+CAPIL*ELECNI(i)
            ELECNJ(J)=ELECNJ(J)+ELECNI(i)
          end if
        end do
        do i=1,MXSHEL
          if (ELECNJ(i).gt.0.) then
            CAPIO(i)=DEXP(CAPILS(i)/ELECNJ(i))
          end if
        end do
        do i=1,MXRAW
          write(  526,2110) i,CAPIN(i)
2110      format('i,CAPIN(i)= ',I5,E12.5)
        end do
        do i=1,MXSHEL
          write(  526,2130) i,CAPIO(i)
2130      format(' i,CAPIO(i)= ',I5,E12.5)
        end do
        do i=1,31
          do ishell=1,MXSHEL
            SCPSUM(i,ishell)=0.0
          end do
          do in=1,MXRAW
            ISHELL=NSHELL(in)
            SCPSUM(i,ishell)=SCPSUM(i,ishell)+SCPROF(i,in)*ELECNI(in)
          end do
        end do
        do ishell=1,MXSHEL
          write(  526,2180) ishell
2180      format(' ',I5,'-th shell. SCPSUM')
          write(  526,2190) (SCPSUM(i,ishell),I=1,31)
2190      format(' ',1P,7E9.2)
        end do
        do i=2,MXSHEL
          ELECNO(i)=ELECNO(i)+ELECNO(i-1)
        end do
        do i=1,MXSHEL
          ELECNO(i)=ELECNO(i)/ELECNO(MXSHEL)
        end do
        do i=1,MXSHEL
          write(  526,2230) i,ELECNO(i)
2230      format(' Shell ,ELECNO=',I5,E12.5)
        end do
        if (ICPROF.eq.4) then
          call PWLF1(NGCS,NALC,1.D-2,1.D2,1.D2,EPCP,ZTHRC,ZEPC,NIPC,DLOG&
         ,DEXP, AXCS,BXCS,2000,MXSHEL,AFCS,BFCS,CFUNS4)
        end if
        if (ICPROF.eq.3) then
          do ishell=1,MXSHEL
            CPROFI(1)=0.0
            QCAP10(1)=0.0
            write(  526,2250) ishell
2250        format(' ICPROF=3,ISHELL=',I5)
            do i=1,31
              AVCPRF(i)=SCPSUM(i,ishell)
            end do
            write(  526,2270)
2270        format(' AVCPROF')
            write(  526,2280) (AVCPRF(i),i=1,31)
2280        FORMAT(' ',7E12.5)
            do i=2,301
              ILOC=(i-2)/10+1
              STEP=(QCAP(ILOC+1)-QCAP(ILOC))*0.1
              AX=QCAP(ILOC)+STEP*FLOAT(I-2-(ILOC-1)*10)
              BX=AX+STEP
              QCAP10(i)=BX
              CPROFI(i)=QD(CPRFIL,AX,BX,'CPRFIL')
            end do
            do i=2,301
              CPROFI(i)=CPROFI(i)+CPROFI(i-1)
            end do
            do i=1,301
              SCPROI(i,ishell)=CPROFI(i)/CPROFI(301)
            end do
          end do

        end if
        write(  526,2320) NGCS,NALC
2320    format(' S-Compton profile was PWLF-ED, NGCS= ',I5,' NALC=',I5)
      end if
      go to 1040
1220  call LAY
      ISSB=ISSBS
      AFACT=AFACTS
      CBAR=CBARS
      SK=SKS
      X0=X0S
      X1=X1S
      IEV=IEVS
      write(  526,2330) ISSB,AFACT,CBAR,SK,X0,X1,IEV
2330  format(' After DECK. ISSB,AFACT,CBAR,SK,X0,X1,IEV=',I5,6E11.4)
      IBOUND=ibounds
      INCOH=incohs
      ICPROF=icprofs
      GASP=gasps
      IRAYL=irayls
      IMPACT=impacts
      IUNRST=iunrsts
      fudgeMS=fudgeMSs
      efracH=efracHs
      efracL=efracLs
      nleg0=nleg0s
      EPSTFL=epstfls
      IEPST=iepsts
      IAPRIM=iaprims
      IAPRFL=iaprfls
      go to 1040
1230  continue
      call PLOT(49,XP,1,AE,UE,NPTS,9)
      call PLOT(71,XP,1,AE,UE,NPTS,9)
      call PLOT(47,XP,1,AE,UE,NPTS,9)
      call PLOT(68,XP,1,AE,UE,NPTS,9)
      call PLOT(46,XP,1,AE,UE,NPTS,9)
      call PLOT(66,XP,1,AE,UE,NPTS,9)
      call PLOT(67,XP,1,AE,UE,NPTS,9)
      call PLOT(77,XP,1,AE,UE,NPTS,9)
      call PLOT(78,XP,1,AE,UE,NPTS,9)
      call PLOT(53,XP,1,AP,UP,NPTS,6)
      call PLOT(51,XP,1,AP,UP,NPTS,6)
      call PLOT(52,XP,1,AP,UP,NPTS,6)
      call PLOT(44,XP,1,AP,UP,NPTS,6)
      write(  526,2340)
2340  format('1')
      go to 1040
1240  continue
      go to 1040
1250  read(  525,2350) NAME
2350  format(6A1)
      IFUN=IFUNT(NAME)
      if (IFUN.le.0) go to 1040
      VALUE=FI(IFUN,XP(1),XP(2),XP(3),XP(4))
      NA=NFARG(IFUN)
      write(  526,2360) VALUE,(FNAME(i,IFUN),i=1,6),(XP(i),i=1,NA)
2360  format(' Function call: ',1P,G15.6,' = ',6A1,' OF ',4G15.6)
      go to 1040
1260  read(  525,2370) NAME,IDFNAM
2370  format(12A1)
      IFUN=IFUNT(NAME)
      if (IFUN.le.0) go to 1040
      ID=IFUNT(IDFNAM)
      if (ID.lt.0) go to 1040
      if (ID.ne.0) IDF=ID
1270  call PLOT(IFUN,XP,IV,VLO,VHI,NPTS,IDF)
      go to 1040

      write(  526,2390) EI,ISUB
2390  FORMAT(' EI=',F14.3,',ISUB=',I3)
      go to (2400,2410,2420,2430,2440,2450),ISUB
2400  call HPLT1(EI,RM,EI-RM,NAMESB,NTIMES,NBINS,NH, 1,54,59,63 )
      go to 1040
2410  call HPLT1(EI,EI/(1.0+2.0*EI/RM),EI,NAMESB,NTIMES,NBINS,NH, 6,40,4&
      2,43 )
      go to 1040
2420  call HPLT1(EI,AP,EI-RM,NAMESB,NTIMES,NBINS,NH, 6,24,30,34 )
      go to 1040
2430  call HPLT1(EI,AE,RM+(EI-RM)*0.5,NAMESB,NTIMES,NBINS,NH, 3,11,13,14&
       )
      go to 1040
2440  call HPLT1(EI,AE,EI,NAMESB,NTIMES,NBINS,NH, 3,20,22,23 )
      go to 1040
2450  PINC=DSQRT(EI**2-RM**2)
      AVE=EI+RM

      go to 1040
1060  call prelastino
1070  write(  526,2470)
2470  format(///' End of file read - exit from pegs5'/'1')
      close(UNIT=525)
      close(UNIT=526)
      close(UNIT=507)
      close(UNIT=528)
      close(UNIT=509)
      close(UNIT=510)
      close(UNIT=511)
      close(UNIT=513)
      close(UNIT=514)
      close(UNIT=515)
      close(UNIT=517)
      close(UNIT=518)
      close(UNIT=519)
      close(UNIT=521)
      close(UNIT=522)
      return
8990  write(526,9990)
      write(*,9990)

      go to 9995
8991  write(526,9991)
      write(*,9991)
9991  format(//,'EOF on user supplied Compton profile data - stopping.')
      go to 9995
8992  write(526,9992)
      write(*,9992)
9992  format(//,'EOF on user supplied interference cs data - stopping.')
      go to 9995
8993  write(526,9993)
      write(*,9993)
9993  format(//,'EOF on user supplied interference ff data - stopping.')
9995  close(526)
      stop
      end
      block data PEGSDAT
      
      include 'include/egs5_h.f'
      include 'pegscommons/bremp2.f'
      include 'pegscommons/dbrpr.f'
      include 'pegscommons/dercon.f'
      include 'pegscommons/elemtb.f'
      include 'pegscommons/elmtbc.f'
      include 'pegscommons/funcs.f'
      include 'pegscommons/funcsc.f'
      include 'pegscommons/lspion.f'
      include 'pegscommons/mimsd.f'
      include 'pegscommons/mixdat.f'
      include 'pegscommons/mxdatc.f'
      include 'pegscommons/molvar.f'
      include 'pegscommons/phpair.f'
      include 'pegscommons/pmcons.f'
      include 'pegscommons/pwlfin.f'
      include 'pegscommons/radlen.f'
      include 'pegscommons/cohcom.f'
      include 'pegscommons/spcomm.f'
      include 'pegscommons/spcomc.f'
      include 'pegscommons/rslts.f'
      include 'pegscommons/thres2.f'
      include 'pegscommons/epstar.f'
      include 'pegscommons/eimpact.f'
      include 'pegscommons/bcom.f'
      include 'pegscommons/cpcom.f'
      include 'pegscommons/sfcom.f'
      include 'pegscommons/mscom.f'
     
      equivalence (MEDTBL(1,1),MEDTB1(1,1))
      equivalence (MEDTBL(1,21),MEDTB2(1,1))
      equivalence (MEDTBL(1,41),MEDTB3(1,1))
      equivalence (MEDTBL(1,61),MEDTB4(1,1))
      equivalence (MEDTBL(1,81),MEDTB5(1,1))
      equivalence (MEDTBL(1,101),MEDTB6(1,1))
      equivalence (MEDTBL(1,121),MEDTB7(1,1))
      equivalence (MEDTBL(1,131),MEDTB8(1,1))
      equivalence (MEDTBL(1,141),MEDTB9(1,1))
      equivalence (MEDTBL(1,151),MEDT10(1,1))
      equivalence (MEDTBL(1,161),MEDT11(1,1))
      equivalence (MEDTBL(1,171),MEDT12(1,1))
      equivalence (MEDTBL(1,181),MEDT13(1,1))
      equivalence (MEDTBL(1,191),MEDT14(1,1))
      equivalence (MEDTBL(1,201),MEDT15(1,1))
      equivalence (MEDTBL(1,211),MEDT16(1,1))
      equivalence (MEDTBL(1,221),MEDT17(1,1))
      equivalence (MEDTBL(1,231),MEDT18(1,1))
      equivalence (MEDTBL(1,241),MEDT19(1,1))
      equivalence (MEDTBL(1,251),MEDT20(1,1))
      equivalence (MEDTBL(1,261),MEDT21(1,1))
      equivalence (MEDTBL(1,271),MEDT22(1,1))
      
      equivalence (STDATA(1,1),STDAT1(1,1))
      equivalence (STDATA(1,21),STDAT2(1,1))
      equivalence (STDATA(1,41),STDAT3(1,1))
      equivalence (STDATA(1,61),STDAT4(1,1))
      equivalence (STDATA(1,81),STDAT5(1,1))
      equivalence (STDATA(1,101),STDAT6(1,1))
      equivalence (STDATA(1,121),STDAT7(1,1))
      equivalence (STDATA(1,141),STDAT8(1,1))
      equivalence (STDATA(1,161),STDAT9(1,1))
      equivalence (STDATA(1,181),STDA10(1,1))
      equivalence (STDATA(1,201),STDA11(1,1))
      equivalence (STDATA(1,221),STDA12(1,1))
      equivalence (STDATA(1,241),STDA13(1,1))
      equivalence (STDATA(1,261),STDA14(1,1))

      end
      double precision function photte(K)
      implicit none
      integer i
      double precision phottz
      double precision K
      include 'include/egs5_h.f'
      include 'pegscommons/mixdat.f'
      include 'pegscommons/mxdatc.f'
      save
      photte=0.0
      do i=1,NE
        photte=PHOTTE+PZ(i)*PHOTTZ(Z(i),K)
      end do
      return
      end
      double precision function phottz(Z,K)
      implicit none
      integer IZ
      double precision PCON, Z, AINTP
      double precision K
      include 'pegscommons/phpair.f'
      include 'pegscommons/pmcons.f'
      include 'pegscommons/molvar.f'
      save
      PCON=1.D-24*(AN*RHO/WM)*RLC
      IZ=Z
      phottz=PCON*AINTP(K,PHE(1,IZ),NPHE(IZ),PHD(1,IZ),1,.TRUE.,.TRUE.)
      return
      end
      subroutine plot(IFUN,XP,IV,EL,EH,NPT,IDF)
      implicit none
      
      double precision DFL, FI, EL, DFH, EH, BDF, YMAX, X, DY
      include 'include/egs5_h.f'
      include 'pegscommons/funcs.f'
      include 'pegscommons/funcsc.f'
      include 'pegscommons/molvar.f'
      include 'pegscommons/dercon.f'
      include 'pegscommons/epstar.f'
      character*4 PBUF(101),ID(5),ORDNL(3,4),ICOM,IRPAR,ICOL,IX,IBL
      double precision XTAB(200),XTABA(18),YSAV(200),XP(4),XQ(5)
      
      data IXTABF/0/
      
      data PBUF/'I',100*' '/,NMAX/200/,NUPL/  526/,IX/'X'/,IBL/' '/
      save
      if (IXTABF.eq.0) then
        itab= 0
        do ixt1=0,6
          do ixt2=1,17
            if (ixt2.ne.11) then
              itab=itab+ 1
              XTAB(itab)=XTABA(ixt2)*10.**(ixt1-3)
            else if (ixt2.eq.11 .and. ixt1.eq.4) then
              itab=itab+ 1
              XTAB(itab)=XTABA(ixt2)*10.**(ixt1-3)
            end if
          end do
        end do
        itab=itab+1
        XTAB(itab)=XTABA(18)*10.**(6-3)
        IXTABF=1
        NU=ITAB
      end if
      IDFI=IDF+1
      NA=NFARG(IFUN)
      DFL=FI(IDF,EL,0.d0,0.d0,0.d0)
      DFH=FI(IDF,EH,0.d0,0.d0,0.d0)
      BDF=(DFH-DFL)/FLOAT(NU-1)
      YMAX=0.0
      do ip=1,NU
        X=XTAB(ip)+RM
        XP(IV)=X
        YSAV(ip)=FI(IFUN,XP(1),XP(2),XP(3),XP(4))/(RLC*RHO)
        YMAX=DMAX1(YSAV(ip),YMAX)
      end do
      DY=YMAX/100.
      ja=0
      do ia=1,NA
        ja=ja+1
        if (ia.ne.IV) then
          XQ(ja)=XP(ia)
          ID(ja)=ICOM
        else
          XQ(ja)=EL
          ID(ja)=ICOL
          ja=ja+1
          XQ(ja)=EH
          ID(ja)=ICOM
        end if
      end do
      ID(ja)=IRPAR
      write(NUPL,100) (FNAME(i,IFUN),i=1,6),(XQ(i),ID(i),i=1,JA)
100   format (' Plot of function ',6A1,'(',5(1P,G15.6,1X,A1) )


      write(NUPL,120) RLC,RHO

      do ip=1,NU
        X=XTAB(ip)
        if (DY.ne.0.0) then
          IBIN=YSAV(ip)/DY+1.0
        else
          IBIN=1
        end if
        if (IBIN.ge.2) PBUF(IBIN)=IX
        write(NUPL,130) IP,X,YSAV(ip),(PBUF(j),j=1,IBIN)
130     format (1X,I3,1P,2G13.6,1X,101A1)
        if (IBIN.ge.2) PBUF(IBIN)=IBL
      end do
      write(521,*) NU
      write(521,140) (XTAB(ip),ip=1,NU)
      write(521,*) NU
      write(521,140) (YSAV(ip),ip=1,NU)
140   FORMAT(5(1P,E15.7))
      return
      end
      subroutine PLOT1(IFUN,XP,IV,EL,EH,NPT,IDF)
      implicit none
    
      double precision DFL, FI, EL, DFH, EH, BDF, YMAX, DF, X, DY
      include 'pegscommons/funcs.f'
      include 'pegscommons/funcsc.f'
      character*4 PBUF(101),ID(5),ORDNL(3,4),ICOM,IRPAR,ICOL,IX,IBL
      double precision YSAV(200),XP(4),XQ(5)
     
      data PBUF/'I',100*' '/,NMAX/200/,NUPL/506/,IX/'X'/,IBL/' '/
      save
      NU=MIN0(NPT,NMAX)
      IDFI=IDF+1
      NA=NFARG(IFUN)
      DFL=FI(IDF,EL,0.d0,0.d0,0.d0)
      DFH=FI(IDF,EH,0.d0,0.d0,0.d0)
      BDF=(DFH-DFL)/FLOAT(NU-1)
      YMAX=0.0
      do ip=1,NU
        i=ip-1
        DF=DFL+BDF*FLOAT(i)
        X=FI(IDFI,DF,0.d0,0.d0,0.d0)
        XP(IV)=X
        YSAV(ip)=FI(IFUN,XP(1),XP(2),XP(3),XP(4))
        YMAX=DMAX1(YSAV(ip),YMAX)
      end do
      DY=YMAX/100.
      ja=0
      do ia=1,NA
        ja=ja+1
        if (ia.ne.IV) then
          XQ(ja)=XP(ia)
          ID(ja)=ICOM
        else
          XQ(ja)=EL
          ID(ja)=ICOL
          ja=ja+1
          XQ(ja)=EH
          ID(ja)=ICOM
        end if
      end do
      ID(ja)=IRPAR
      write(NUPL,100) (FNAME(i,IFUN),i=1,6),(XQ(i),ID(i),i=1,ja)
100   format (' Plot of function ',6A1,'(',5(1P,G15.6,1X,A1) )

      do ip=1,NU
        i=ip-1
        X=FI(IDFI,DFL+BDF*FLOAT(i),0.d0,0.d0,0.d0)
        if (DY.ne.0.0) then
          IBIN=YSAV(ip)/DY+1.0
        else
          IBIN=1
        end if
        if (IBIN.ge.2) PBUF(IBIN)=IX
        write(NUPL,120) ip,X,YSAV(ip),(PBUF(j),j=1,IBIN)
120     format (1X,I3,1P,2G13.6,1X,101A1)
        if (IBIN.ge.2) PBUF(IBIN)=IBL
      end do
      return
      end
      subroutine PMDCON
      implicit none
      double precision FSCI
      include 'pegscommons/pmcons.f'
      include 'pegscommons/dercon.f'
      save
      PI=3.1415926535897932D+0
      C=2.99792458D+10
      RME=9.10938188D-28
      HBAR=1.054571596E-27
      ECGS=4.8032068D-10
      EMKS=1.602176462D-19
      AN=6.02214199D+23
      RADDEG=180./PI
      FSC = ECGS**2/(HBAR*C)
      FSCI=1./FSC
      ERGMEV = (1.D+6)*(EMKS*1.D+7)
      R0 = (ECGS**2)/(RME*C**2)
      RM = RME*C**2/ERGMEV
      RMT2 = RM*2.0
      RMSQ = RM*RM
      A22P9 = RADDEG*DSQRT(4.*PI*AN)*ECGS**2/ERGMEV
      A6680 = 4.0*PI*AN*(HBAR/(RME*C))**2*(0.885**2/(1.167*1.13))
      return
      end
      double precision function PSIG(E)
      implicit none
      double precision BREMTM, E, BHABTM, ANIHTM
      save
      PSIG=BREMTM(E)+BHABTM(E)+ANIHTM(E)
      return
      end
      subroutine PWLF1(NI,NIMX,XL,XU,XR,EP,ZTHR,ZEP,NIP,XFUN,XFI, AX,BX,&
      NALM,NFUN,AF,BF,VFUNS)
      implicit none
      integer NALM, NFUN, NL, NU, IPRN, NJ, NIMX, NIP, NI, NK
      double precision XL, XU, XR, EP, ZTHR, ZEP, REM, AX, BX, AF, BF
      external XFUN,XFI,VFUNS
      dimension AF(NALM,NFUN),BF(NALM,NFUN),ZTHR(NFUN),ZEP(NFUN)
      logical QFIT
      save
      NL=0
      NU=1
      IPRN=0
100   continue
        NJ=MIN0(NU,NIMX)

        if (NU.ge.NIMX) then
          write(  526,110) NIMX,EP

          NI=NJ
          return
        end if
        NL=NU
        NU=NU*2
      go to 100
120   continue
      NU=NJ
130   if (NU.le.NL+1) go to 140
        NJ=(NL+NU)/2
        NK=NJ
      go to 130
140   continue
      NI=NU
      if (NI.eq.NJ) return

150   format(' Catastrophe---does not fit when it should,NI=',I5)
      return
      end
      double precision function QD(F,A,B,MSG)
      implicit none
      integer IER
      double precision A, B
      external F
      double precision DCADRE,ADUM,BDUM,ERRDUM
      character*6 MSG
      save
      ADUM=A
      BDUM=B
      QD=DCADRE(F,ADUM,BDUM,1.D-16,1.D-5,ERRDUM,IER)
      if (IER.gt.66) then
        write(510,100) IER,MSG,A,B,QD,ERRDUM

      end if
      return
      end
      logical function QFIT(NJ,XL,XH,XR,EP,ZTHR,ZEP,REM,NJP,XFUN,XFI, AX&
      ,BX,NALM,NFUN,AF,BF,VFUNS,IPRN)
      implicit none
      
      external XFUN,XFI,VFUNS
      dimension FSXL(NALM),FSXH(NALM),FIP(NALM),FFIP(NALM),AFIP(NALM)
      dimension RE(NALM),AER(NALM)
      dimension AF(NALM,NFUN),BF(NALM,NFUN),ZTHR(NFUN),ZEP(NFUN)
      data EPS1/1.d-15/
      data NKP/3/
      save
      if(nfun.gt.500) then  ! this number is equal to NFUNS
       write(  *,1100) nfun
       write(  *,1101)
       write(  *,1102)
       stop
      endif
1100  format('Error in compton profile data NFUN (mxshel) =',i4,' > 79')
1101  format('You have to decrease the number of elements per material')
1102  format('or you have to set iprofr=0 in the [parameters] section')
      if (XH.le.XL) then
        write(  526,100) XL,XH
100     format(' QFIT error:XL should be < XH. XL,XH=',2G14.6)
        QFIT=.false.
        return
      end if
      XS=DMAX1(XL,DMIN1(XH,XR))
      NI=NJ-2
      if (((XS.eq.XL.or.XS.eq.XH).and.NI.ge.1).or.NI.ge.2) then
        XFL=XFUN(XL)
      else
        QFIT=.false.
        return
      end if
      XFH=XFUN(XH)
      XFS=XFUN(XS)
      XM=DMAX1(XFH-XFS,XFS-XFL)
      DX=XFH-XFL
      IF(DABS(XM-DX).LE.EPS1*DX) THEN
        W=XM/DMAX1(1.d0,AINT(NI*1.d0))
      else
        W=XM/DMAX1(1.d0,AINT(NI*XM/DX))
      end if
      NI=NI-AINT(NI-DX/W)
      NIP=MAX0(NKP,(NJP+NI-1)/NI)
      NIP=(NIP/2)*2+1
      if (XFH-XFS.le.XFS-XFL) then
        XLL=XFL
      else
        XLL=XFH-NI*W
      end if
      AX=1./W
      BX=2.-XLL*AX
      REM=0.0
      QFIT=.true.
      SXFL=DMAX1(XLL,XFL)
      ISUB=0
      XSXF=XFI(SXFL)
      call VFUNS(XSXF,FSXL)

110   format(' QFIT:ISUB,SXF,XSXF,FSX()=',I4,1P,9G11.4/(1X,12G11.4))
      do isub=1,NI
        JSUB=isub+1
        SXFH=DMIN1(XLL+W*isub,XH)
        XSXF=XFI(SXFH)
        call VFUNS(XSXF,FSXH)
    
        DSXF=SXFH-SXFL
        do ifun=1,NFUN
          AF(JSUB,ifun)=(FSXH(ifun)-FSXL(ifun))/DSXF
          BF(JSUB,ifun)=(FSXL(ifun)*SXFH-FSXH(ifun)*SXFL)/DSXF
        end do
        WIP=DSXF/(NIP+1)
        do ip=1,NIP
          SXFIP=SXFL+ip*WIP
          XIP=XFI(SXFIP)
          call VFUNS(XIP,FIP)
          do ifun=1,NFUN
            FFIP(ifun)=AF(JSUB,ifun)*SXFIP+BF(JSUB,ifun)
            AFIP(ifun)=dabs(FIP(ifun))
            AER(ifun)=dabs(FFIP(ifun)-FIP(ifun))
            RE(ifun)=0.0
            if (FIP(ifun).ne.0.0) then
              RE(ifun)=AER(ifun)/AFIP(ifun)
            end if
            if (AFIP(ifun).ge.ZTHR(ifun)) then
              REM=dmax1(REM,RE(ifun))
            else if (AER(ifun).gt.ZEP(ifun)) then
              QFIT=.false.
            end if
          end do
          if (IPRN.ne.0) then

          end if
        end do
        SXFL=SXFH
        do ifun=1,NFUN
          FSXL(ifun)=FSXH(ifun)
        end do
      end do
      do ifun=1,NFUN
        AF(1,ifun)=AF(2,ifun)
        BF(1,ifun)=BF(2,ifun)
        AF(NI+2,ifun)=AF(NI+1,ifun)
        BF(NI+2,ifun)=BF(NI+1,ifun)
      end do
      QFIT=QFIT.AND.REM.LE.EP
      NJ=NI+2
      return
      end
      subroutine RDSCPR
      implicit none
      integer MXRAW, MXSHEL, NSHELL, i
      double precision ELECNI, CAPIN, SCPROF, QCAP, PZT
      include 'include/egs5_h.f'
      include 'pegscommons/mixdat.f'
      include 'pegscommons/mxdatc.f'
      character FILENM*200,NOZ*3

      NAMELIST/SCPRDT/MXRAW,MXSHEL,ELECNI,NSHELL,CAPIN,SCPROF,QCAP
      common /paran/ icfn(100), ilfn(100), chfn(100)
      character chfn*200
      integer icfn, ilfn
      save
      do i=1,NE
        WRITE(NOZ,'(I3.3)') NINT(Z(i))
       
        write(  526,100) FILENM
100     format(' Reading ',A50)
        open(UNIT=530,FILE=FILENM,STATUS='old')
        PZT=PZ(i)
        read(530,SCPRDT)
        call ADSCPR(MXRAW,MXSHEL,ELECNI,NSHELL,CAPIN,SCPROF,QCAP,PZT)
        close(530)
      end do
      call WTSCPR
      return
      end
      subroutine RFUNS(E,V)
      implicit none
      double precision AINTP, E
      include 'pegscommons/cohcom.f'
      double PRECISION V(1)
      save
      V(1)=AINTP(E,AFFI(1),97,XVAL(1),1,.TRUE.,.TRUE.)
      return
      end
      subroutine RFUNS2(E,V)
      implicit none
      double precision AINTP, E
      include 'pegscommons/cohcom.f'
      double precision V(1)
      save
      V(1)=AINTP(E,XVAL(1),97,AFAC2(1),1,.TRUE.,.TRUE.)
      return
      end
      subroutine SFUNS(E,V)
      implicit none
      double precision AINTP, E
      include 'pegscommons/bcom.f'
      include 'pegscommons/sfcom.f'
      double precision V(1)
      save
      V(1)=AINTP(E,XSVAL(1),41,SCATZ(1),1,.TRUE.,.TRUE.)
      return
      end
      subroutine SPINIT
      implicit none
      integer im, j, IZ, ie, i, ICHECK, IESPEL, IPEGEL,ios
     
      include 'pegscommons/pmcons.f'
      include 'pegscommons/spcomm.f'  !<- ada.comment "nmed"in spcomm.f . set val = 278 in BLOCKDATA.PEGS .
      include 'include/egs5_h.f'
      include 'pegscommons/spcomc.f'
      include 'pegscommons/dercon.f'
      include 'pegscommons/molvar.f'
      include 'pegscommons/mixdat.f'
      include 'pegscommons/mxdatc.f'
      include 'pegscommons/elemtb.f'
      include 'pegscommons/elmtbc.f'
      include 'pegscommons/lspion.f'
      include 'pegscommons/epstar.f'
      include 'pegscommons/thres2.f'
      include 'err.inc'
      double precision IMEV
      common /paran/ icfn(100), ilfn(100), chfn(100)
      character chfn*200
      integer icfn, ilfn
      save
      TOLN10=2.0*DLOG(10.d0)
      IM=-100
      if (EPSTFL .lt. 0 .or. EPSTFL .gt. 1) then
        EPSTFL = 0
      end if
      write(526,9) EPSTFL
9     format(' EPSTFL=',I15)
      if (EPSTFL.eq.0) then
        if (ISSB.ne.0) then
          if ( AFACT.eq.0.0 .or. CBAR.le.0.0 .or. SK.eq.0.0 .or. X0.eq.&
         0.0 .or. X1.eq.0.0 .or. IEV.eq.0.0 ) then
            write(  526,100)
            write(  *,100)

            close(526)
            stop
          end if
          IMEV=IEV*1.D-6
          VPLASM=DSQRT(EDEN*R0*C**2/PI)
          IM=-1
        else
          if (ISSB.eq.0.and.(AFACT.ne.0.0.or.CBAR.ne.0.0.or.SK.ne.0.0.&
         or.X0.ne.0.0.or.X1.ne.0.0.or.IEV.ne.0.0)) then
            write(  526,110)
            write(  *,110)

            close(526)
            stop
          end if
          DO 120 im=1,NMED
            do j=1,LMED
              if (IDSTRN(j).ne.MEDTBL(j,im)) go to 120
            end do
            AFACT=STDATA(1,im)
            SK=STDATA(2,im)
            X0=STDATA(3,im)
            X1=STDATA(4,im)
            IEV=STDATA(5,im)
            CBAR=STDATA(6,im)
            DELTA0=STDATA(7,im)
            IMEV=IEV*1.0D-6
            VPLASM=DSQRT(EDEN*R0*C**2/PI)
            go to 150
120       continue
          DELTA0=0.0
          IM=0
          if (NE.eq.1) then
            IZ=Z(1)
            if (IZ.eq.1.or.IZ.eq.7.or.IZ.eq.8) then
              write(  526,130)
              write(  *,130)

              close(526)
              stop
            end if
            IEV=ITBL(IZ)
          else
            ALIADG=0.0
            do ie=1,NE
              IZ=Z(ie)
              if (IZ.eq.1) then
                IEV=19.2
              else if (IZ.eq.6) then
                if (GASP.eq.0.0) then
                  IEV=81.0
                else
                  IEV=70.0
                end if
              else if (IZ.eq.7) then
                IEV=82.0
              else if (IZ.eq.8) then
                if (GASP.eq.0.0) then
                  IEV=106.0
                else
                  IEV=97.0
                end if
              else if (IZ.eq.9) then
                IEV=112.0
              else if (IZ.eq.17) then
                IEV=180.0
              else
                IEV=1.13*ITBL(IZ)
              end if
              ALIADG=ALIADG + PZ(IE)*Z(IE)*DLOG(IEV)
            end do
            ALIADG=ALIADG/ZC
            IEV=DEXP(ALIADG)
          end if
          IMEV=IEV*1.0D-6
          if (GASP.eq.0.0) then
            EDENL=EDEN
          else
            EDENL=EDEN/GASP
          end if
          VPLASM = DSQRT(EDENL*R0*C**2/PI)
          CBAR=1. + 2.*DLOG(IMEV/(HBAR*2*PI*VPLASM/ERGMEV))
          if (NE.eq.1.and.IDINT(Z(1)).eq.2.and.GASP.ne.0.0) then
            X0=2.191
            X1=3.0
            SK=3.297
          else if&
           (NE.eq.2.and.IDINT(Z(1)).eq.1.and.IDINT(Z(2)).eq.1) then
            if (GASP.eq.0.0) then
              X0=0.425
              X1=2.0
              SK=5.949
            else
              X0=1.837
              X1=3.0
              SK=4.754
            end if
          else
            SK=3.0
            if (GASP.eq.0.0) then
              if (IEV.lt.100.0) then
                if (CBAR.lt.3.681) then
                  X0=0.2
                  X1=2.0
                else
                  X0=0.326*CBAR - 1.0
                  X1=2.0
                end if
              else
                if (CBAR.lt.5.215) then
                  X0=0.2
                  X1=3.0
                else
                  X0=0.326*CBAR - 1.5
                  X1=3.0
                end if
              end if
              if (X0.ge.X1) then
                write(  526,140) X0,X1,CBAR
                write(  *,140) X0,X1,CBAR

                close(526)
                stop
              end if
            else
              if (CBAR.lt.10.0) then
                X0=1.6
                X1=4.0
              else if (CBAR.lt.10.5) then
                X0=1.7
                X1=4.0
              else if (CBAR.lt.11.0) then
                X0=1.8
                X1=4.0
              else if (CBAR.lt.11.5) then
                X0=1.9
                X1=4.0
              else if (CBAR.lt.12.25) then
                X0=2.0
                X1=4.0
              else if (CBAR.lt.13.804) then
                X0=2.0
                X1=5.0
              else
                X0=0.326*CBAR - 2.5
                X1=5.0
              end if
            end if
          end if
        end if
150     if (GASP.ne.0.0) then
          ALGASP=DLOG(GASP)
          CBAR=CBAR - ALGASP
          X0=X0 - ALGASP/TOLN10
          X1=X1 - ALGASP/TOLN10
        end if
        if (IM.eq.0) then
          AFACT=(CBAR - TOLN10*X0)/(X1 - X0)**SK
        end if
      else
        if(IDSTRN(1).eq.'H'.and.IDSTRN(2).eq.'2'.and.IDSTRN(3).eq.'O')&
        then
       
        elseif(IDSTRN(1).eq.'C'.and.IDSTRN(2).eq.'-') then
        
        elseif(IDSTRN(1).eq.'A'.and.IDSTRN(2).eq.'I'.and.&
        IDSTRN(3).eq.'R') then
         
        else
         ErrCha = ''
         ErrID = 'L:12970/R:SPINIT/F:egs5.f' !E86_028_001
         call ErrWrite(ErrID,ErrCha)
    
         stop
        endif
        read(520,160,ERR=9991) EPSTTL
160     format(80A1)
      
        read(520,*,ERR=9991) (EPSTEN(i),EPSTD(i),i=1,NEPST)
        go to 9993
9991    write(526,9992)
        write(*,9992)
9992    format(/,/,' *****END-OF-FILE on epstar.dat ')
        close(526)
        stop
9993    close(520) ! T.Sato 2020/07/02
        if (NEPST.gt.150) then
          write(  526,170) NEPST
          write(  *,170) NEPST
170       format(//' *****NEPST=',I4,' is greater than the 150 allowed')
          close(526)
          stop
        end if
        do i=1,NEPST
          EPSTEN(i) = EPSTEN(i) + RM
        end do
        IMEV = IEV*1.D-06
        if ( AE .lt. EPSTEN(1)) then
          write(  526,180) EPSTEN(1),AE

          EPSTEN(1) = AE
        end if
        if ( UE .gt. EPSTEN(NEPST)) then
          write(  526,190) EPSTEN(NEPST),UE

          EPSTEN(NEPST) = UE
        end if
        ICHECK=0
        TLRNCE=0.01
        if (NELEPS.ne.NE) ICHECK=1

        if (EPSTRH .lt. 0.0) then
          TLRNCE=1.0
          EPSTRH=ABS(EPSTRH)
          end if
        if (NELEPS.ne.NE) ICHECK=1

        if ((ICHECK.eq.0).and.(EPSTRH.lt.(1.0-TLRNCE)*RHO)) ICHECK=2
        if ((ICHECK.eq.0).and.(EPSTRH.GT.(1.0+TLRNCE)*RHO)) ICHECK=3
        EPSTWT = 0.0
        do i=1,NE
          EPSTWT = EPSTWT + RHOZ(i)
        end do
        if (EPSTWT.eq.0.0) then
          write(  526,200)

        end if
        if (ICHECK.eq.0) then
          IESPEL=0
          ICHECK=4
210       continue
            IESPEL=IESPEL+1
            IPEGEL=0
220         continue
              IPEGEL=IPEGEL+1
              if (DINT(Z(IPEGEL)).eq.ZEPST(IESPEL)) then
                ICHECK=0
                go to 230
              end if
              if (IPEGEL.ge.NE) go to 230
            go to 220
230         continue
            if ((ICHECK.eq.0)  .and. (WEPST(IESPEL).lt.((1.0-TLRNCE)*&
           RHOZ(IPEGEL)/EPSTWT))) ICHECK=5
            if ((ICHECK.eq.0)  .and. (WEPST(IESPEL).gt.((1.0+TLRNCE)*&
           RHOZ(IPEGEL)/EPSTWT))) ICHECK=6
            if (IESPEL.GE.NELEPS) go to 240
          go to 210
240       continue
        end if
        if (ICHECK.ge.1) then
          write(  526,250)
          write(  *,250)

          write(*,'(''epstar.dat mismatches with : '',80A1)') EPSTTL
          if(icheck.eq.1) then
           write(*,*) 'NELEPS should be equal to NE',NELEPS,NE
          elseif(icheck.eq.2) then
           write(*,*) 'EPSTRH < (1.0-TLRNCE)*RHO',EPSTRH,TLRNCE,RHO
          elseif(icheck.eq.3) then
           write(*,*) 'EPSTRH > (1.0-TLRNCE)*RHO',EPSTRH,TLRNCE,RHO
          elseif(icheck.eq.5) then
      
          elseif(icheck.eq.6) then

          endif
          close(526)
          stop
        end if
      end if
      SPC1=2.*PI*R0**2*RM*EDEN*RLC
      SPC2=DLOG((IMEV/RM)**2/2.0)
      write(  526,260)
260   format(//' Parameters computed in SPINIT.'//1X,64('-'))
      if (IM.eq.0) then
        write(  526,270)

      else if (IM.gt.0) then
        write(  526,280)

      else if (IM .eq. -1) then
        write(  526,290)

      else
        write(  526,300) EPSTTL
300     format(' Density effect read in directly:'/T10,80A1)
      end if
      write(  526,310)
310   format(1X,64('-')/)
      write(  526,320) IEV
320   format(/' Adjusted mean ionization = ',F8.2,' eV'/1X,38('-')//)
      if (EPSTFL .eq. 0) then
        V4110=IEV
        write(  526,330) V4110
330     format(' IEV=',1P,G15.7)
        V4130=VPLASM
        write(  526,340) V4130
340     format(' VPLASM=',1P,G15.7)
        V4150=CBAR
        write(  526,350) V4150
350     format(' CBAR=',1P,G15.7)
        V4170=X0
        write(  526,360) V4170
360     format(' X0=',1P,G15.7)
        V4190=X1
        write(  526,370) V4190
370     format(' X1=',1P,G15.7)
        V4210=SK
        write(  526,380) V4210
380     format(' SK=',1P,G15.7)
        V4230=AFACT
        write(  526,390) V4230
390     format(' AFACT=',1P,G15.7)
        V4240=DELTA0
        WRITE(  526,400)V4240
400     FORMAT(' DELTA0=',1P,G15.7)
      end if
      V4250=SPC1
      write(  526,410) V4250
410   format(' SPC1=',1P,G15.7)
      V4270=SPC2
      write(  526,420) V4270
420   format(' SPC2=',1P,G15.7)
      return
      end
      double precision function SPIONB(E0,EE,POSITR)
      implicit none
      integer i
      
      logical POSITR
      include 'include/egs5_h.f'
      include 'pegscommons/dercon.f'
      include 'pegscommons/lspion.f'
      include 'pegscommons/epstar.f'
      save
      G=E0/RM
      EEM=EE/RM-1.
      T=G-1
      ETA2=T*(G+1.)
      BETA2=ETA2/G**2
      ALETA2=DLOG(ETA2)
      X=0.21715*ALETA2
      if (.NOT.POSITR) then
        D=DMIN1(EEM,0.5*T)
        FTERM=-1.-BETA2+DLOG((T-D)*D)+T/(T-D) +(D*D/2.+(2.*T+1.)*DLOG(1.&
       -D/T))/(G*G)
      else
        D=DMIN1(EEM,T)
        TP2=T+2.
        D2=D*D
        D3=D*D2
        D4=D*D3
        FTERM=DLOG(T*D)-(BETA2/T)*( T + 2.*D - (3.*D2/2.)/TP2 -(D-D3/3.)&
       /(TP2*TP2)-(D2/2.-T*D3/3.+D4/4.)/TP2**3)
      end if
      if (EPSTFL .eq. 0) then
        if (X.le.X0) then
          DELTA=DELTA0*10**(2.0*(X-X0))
        else if (X.lt.X1) then
          DELTA=TOLN10*X - CBAR + AFACT*(X1 - X)**SK
        else
          DELTA=TOLN10*X - CBAR
        end if
      else
        if (E0 .ge. EPSTEN(IEPST)) then
          if (E0 .eq. EPSTEN(IEPST)) then
            go to 100
          end if
          do i=IEPST,NEPST-1
            if (E0.lt.EPSTEN(i+1)) then
              IEPST = I
              go to 100
            end if
          end do
          IEPST = NEPST
          go to 100
        else
          do i=IEPST,2,-1
            if (E0 .ge. EPSTEN(i-1)) then
              IEPST = I-1
              go to 100
            end if
          end do
          IEPST = 1
        END IF
100    if (IEPST .lt. NEPST) then
          DELTA = EPSTD(IEPST) + (E0 - EPSTEN(IEPST))/ (EPSTEN(IEPST+1)&
         - EPSTEN(IEPST)) * (EPSTD(IEPST+1) - EPSTD(IEPST))
        else
          DELTA = EPSTD(NEPST)
        end if
      end if
      SPIONB=(SPC1/BETA2)*(DLOG(T + 2.) - SPC2 + FTERM - DELTA)
      return
      end
      double precision function spione(E0,EE)
      implicit none
      double precision spionb, E0, EE
      save
      spione=spionb(E0,EE,.FALSE.)
      return
      end
      double precision function spionp(E0,EE)
      implicit none
      double precision spionb, E0, EE
      save
      spionp=spionb(E0,EE,.TRUE.)
      return
      end
      double precision function sptote(E0,EE,EG)
      implicit none
      double precision spione, E0, EE, brmstm, EG
      include 'pegscommons/thres2.f'
      save
      if (IUNRST.eq.0) then
        sptote=spione(E0,EE)+brmstm(E0,EG)
      else if (IUNRST.eq.1) then
        sptote=spione(E0,E0)
      else if (IUNRST.eq.2) then
        sptote=spione(E0,E0)+brmstm(E0,E0)
      else if (IUNRST.eq.3) then
        sptote=spione(E0,E0)+brmstm(E0,EG)
      else if (IUNRST.eq.4) then
        sptote=spione(E0,EE)+brmstm(E0,E0)
      else if (IUNRST.eq.5)  then
        sptote=brmstm(E0,E0)
      else if (IUNRST.eq.6) then
        sptote=brmstm(E0,EG)
      else if (IUNRST.eq.7) then
        sptote=spione(E0,EE)
      end if
      return
      end
      double precision function sptotp(E0,EE,EG)
      implicit none
      double precision spionp, E0, EE, brmstm, EG
      include 'pegscommons/thres2.f'
      save
      if (IUNRST.eq.0) then
        sptotp=spionp(E0,EE)+brmstm(E0,EG)
      else if (IUNRST.eq.1) then
        sptotp=spionp(E0,E0)
      else if (IUNRST.eq.2) then
        sptotp=spionp(E0,E0)+brmstm(E0,E0)
      else if (IUNRST.eq.3) then
        sptotp=spionp(E0,E0)+brmstm(E0,EG)
      else if (IUNRST.eq.4) then
        sptotp=spionp(E0,EE)+brmstm(E0,E0)
      else if (IUNRST.eq.5) then
        sptotp=brmstm(E0,E0)
      else if (IUNRST.eq.6) then
        sptotp=brmstm(E0,EG)
      else if (IUNRST.eq.7) then
        sptotp=spionp(E0,EE)
      end if
      return
      end
      double precision function tmxb(E)
      implicit none
      double precision ESQ, E, BETA2, PX2
      include 'pegscommons/dercon.f'
      include 'pegscommons/molvar.f'
      save
      ESQ=E**2
      BETA2=1.0-RMSQ/ESQ
      PX2=ESQ*BETA2/XCC**2
      tmxb=PX2*BETA2/DLOG(BLCC*PX2)
      return
      end
      double precision function tmxde2(E)
      implicit none
      double precision ESQ, E, BETASQ, TMXB
      include 'pegscommons/dercon.f'
      save
      ESQ=E**2
      BETASQ=1.0-RMSQ/ESQ
      tmxde2=TMXB(E)/(ESQ*BETASQ**2)
      return
      end
      double precision function tmxs(E)
      implicit none
      double precision SAFETY, TABSMX, TMXB, E
      data SAFETY/0.8/,TABSMX/10.0/
      save
      tmxs=DMIN1(TMXB(E)*SAFETY,TABSMX)
      return
      end
      subroutine wtscpr
      implicit none
      integer i, j
      include 'include/egs5_h.f'
      include 'pegscommons/cpcom.f'
      common /paran/ icfn(100), ilfn(100), chfn(100)
      character chfn*200
      integer icfn, ilfn
      save
     
      write(531,'(1H ,A)') '&SCPRDT'
      write(531,'(1H ,A)') 'QCAP='
      write(531,'((1H ,7(F9.2,A)))') (QCAP(I),',',I=1,31)
      write(531,'(1H ,A,I5)') 'MXRAW=',MXRAW
      Write(531,'(1H ,A)') 'ELECNI='
      write(531,'((1H ,7(1PE9.3,A)))') (ELECNI(I),',',I=1,MXRAW)
      write(531,'(1H ,A,I5)') 'MXSHEL=',MXSHEL
      write(531,'(1H ,A)') 'NSHELL='
      write(531,'((1H ,14(I4,A)))') (NSHELL(I),',',I=1,MXRAW)
      write(531,'(1H ,A)') 'CAPIN='
      write(531,'((1H ,5(1PE11.5,A)))') (CAPIN(I),',',I=1,MXRAW)
      write(531,'(1H ,A)') 'SCPROF='
      do i=1,MXRAW+1
        write(531,'((1H ,7(1PE9.3,A)))') (SCPROF(j,i),',',j=1,31)
      end do
      write(531,'(1H ,A)') '/END'
      endfile 531
      close(531)
      return
      end
      double precision function xsif (Z)
      implicit none
      integer IZ
      double precision Z, FCOULC
      include 'pegscommons/radlen.f'
      save
      if (Z.le.4.0) then
        IZ=Z
        xsif=ALRADP(IZ)/(ALRAD(IZ)-FCOULC(Z))
      else
        xsif=dlog(A1440*Z**(-2./3.))/(dlog(A183*Z**(-1./3.))-FCOULC(Z))
      end if
      return
      end
      double precision function ztbl(IASYM)
      implicit none
      integer ie
      include 'pegscommons/elemtb.f'
      include 'pegscommons/elmtbc.f'
      character*4 IASYM,IA
      data IA/'A'/
      save
      if (IASYM.eq.IA) then
        ztbl=18.0
        return
      end if
      do ie=1,NET
        if (IASYM.eq.ASYMT(ie)) then
          ztbl=ie
          return
        end if
      end do
      write(  526,100) IASYM,NET
100   format(1X,A2,' Not an atomic symbol for an element with Z LE ',I3)
      ZTBL=0.0
      return
      end
      subroutine prelastino
      USE PEGS_DCSSTR_MOD
      USE egs5_media_mod
      implicit none
      include 'include/egs5_h.f'
      integer i,n,didGS
      didGS = 0
      write(517,*) nsdcs
      do n = 1, nsdcs
        write(517,5001) (mednam(n,i),i=1,24)

      end do
5001  format(' MEDIUM=',24A1)
      return
      END
      SUBROUTINE egs5SPLINE(X,Y,A,B,C,D,S1,SN,N)
      implicit none
      include 'err.inc'
      integer N
      double precision  X(N),Y(N),A(N),B(N),C(N),D(N), S1, SN
      integer i, k, n1, n2
      double precision r, h, hi, si, si1
      IF(N.LT.4) THEN
        WRITE(506,10) N
        ErrCha = ''
        ErrID = 'L:13494/R:egs5SPLINE/F:egs5.f' !E86_029_001
        call ErrWrite(ErrID,ErrCha)
        WRITE(*,10) N

        STOP
      ENDIF
      N1=N-1
      N2=N-2
      DO I=1,N1
        IF(X(I+1)-X(I).LT.1.0D-25) THEN
          WRITE(506,11)
          ErrCha = ''
          ErrID = 'L:13508/R:egs5SPLINE/F:egs5.f' !E86_030_001
          call ErrWrite(ErrID,ErrCha)
          WRITE(*,11)
   11     FORMAT(5X,'Spline x values not in increasing order. Stop.')
          STOP
        ENDIF
        A(I)=X(I+1)-X(I)
        D(I)=(Y(I+1)-Y(I))/A(I)
      ENDDO
      DO I=1,N2
        B(I)=2.0D0*(A(I)+A(I+1))
        K=N1-I+1
        D(K)=6.0D0*(D(K)-D(K-1))
      ENDDO
      D(2)=D(2)-A(1)*S1
      D(N1)=D(N1)-A(N1)*SN
      DO I=2,N2
        R=A(I)/B(I-1)
        B(I)=B(I)-R*A(I)
        D(I+1)=D(I+1)-R*D(I)
      ENDDO
      D(N1)=D(N1)/B(N2)
      DO I=2,N2
        K=N1-I+1
        D(K)=(D(K)-A(K)*D(K+1))/B(K-1)
      ENDDO
      D(N)=SN
      SI1=S1
      DO I=1,N1
        SI=SI1
        SI1=D(I+1)
        H=A(I)
        HI=1.0D0/H
        
      ENDDO
      RETURN
      END
      double precision FUNCTION SUMGA(FCT,XL,XU,TOL)
      implicit none
      integer NP, NST, NCALLS, I1, I2, I3, ICALL, LH, I, LHN
      PARAMETER(NP=10,NST=128,NCALLS=20000)
      
      external FCT
     
      CTOL=MIN(MAX(TOL,1.0D-13),1.0D-2)
      PTOL=0.01D0*CTOL
      ERR=1.0D35
      H=XU-XL
      SUMGA=0.0D0
      A=0.5D0*(XU-XL)
      B=0.5D0*(XL+XU)
      C=A*X(1)
      D=W(1)*(FCT(B+C)+FCT(B-C))
      DO I1=2,NP
        C=A*X(I1)
        D=D+W(I1)*(FCT(B+C)+FCT(B-C))
      ENDDO
      ICALL=NP+NP
      LH=1
      S(1)=D*A
      XR(1)=XL
    1 CONTINUE
      HO=H
      H=0.5D0*H
      SUMR=0.0D0
      LHN=0
      DO I=1,LH
        SI=S(I)
        XA=XR(I)
        XB=XA+H
        XC=XA+HO
        A=0.5D0*(XB-XA)
        B=0.5D0*(XB+XA)
        C=A*X(1)
        D=W(1)*(FCT(B+C)+FCT(B-C))
        DO I2=2,NP
          C=A*X(I2)
          D=D+W(I2)*(FCT(B+C)+FCT(B-C))
        ENDDO
        S1=D*A
        A=0.5D0*(XC-XB)
        B=0.5D0*(XC+XB)
        C=A*X(1)
        D=W(1)*(FCT(B+C)+FCT(B-C))
        DO I3=2,NP
          C=A*X(I3)
          D=D+W(I3)*(FCT(B+C)+FCT(B-C))
        ENDDO
        S2=D*A
        ICALL=ICALL+4*NP
        S12=S1+S2
        IF(ABS(S12-SI).LE.MAX(PTOL*ABS(S12),1.0D-25)) THEN
          SUMGA=SUMGA+S12
        ELSE
          SUMR=SUMR+S12
          LHN=LHN+2
          IF(LHN.GT.NST) GO TO 2
          SN(LHN)=S2
          XRN(LHN)=XB
          SN(LHN-1)=S1
          XRN(LHN-1)=XA
        ENDIF
        IF(ICALL.GT.NCALLS) GO TO 2
      ENDDO
      ERR=ABS(SUMR)/MAX(ABS(SUMR+SUMGA),1.0D-25)
      IF(ERR.LT.CTOL.OR.LHN.EQ.0) RETURN
      LH=LHN
      DO I=1,LH
        S(I)=SN(I)
        XR(I)=XRN(I)
      ENDDO
      GO TO 1
    2 CONTINUE
      WRITE(506,11)
   11 FORMAT(/2X,'>>> SUMGA. Gauss adaptive-bipartition quadrature.')
      WRITE(506,12) XL,XU,TOL
   12 FORMAT(2X,'XL =',1P,E19.12,',  XU =',E19.12,',  TOL =',E8.1)
      WRITE(506,13) ICALL,SUMGA,ERR,LHN
   
      RETURN
      END
      subroutine wmsfit(k1start,dk1)
      USE egs5_mscon_mod
      implicit none
      double precision k1start(2), dk1(2)
      include 'include/egs5_h.f'
      include 'include/egs5_cdcsep.f'
      include 'pegscommons/mxdatc.f'
      integer ipart, iener, iang, ik1
      write(517,5001) medium
5001  format(' MEDIUM=',24A1)
      write(517,*) 'MS fitting coefficients for this media'
      write(517,*) 'Total number of energy steps:'
      write(517,*) nmscate
      write(517,*) 'First energy decade:'
      write(517,*) decade1
      write(517,*) 'Number of energy steps inside each decade:'
      write(517,*) NDEC
      write(517,*) 'Number of energy steps to skip in first decade:'
      write(517,*) joffset
      write(517,*) 'Number of equally probably angle bins'
      write(517,*) NBFIT
      write(517,*) 'Number of equally spaced angle bins'
      write(517,*) NEXFIT
      do ipart = 1,2
        if(ipart.eq.1) then
          write(517,*) 'Electrons'
        else
          write(517,*) 'Positrons'
        endif
        do iener = 1,nmscate
          write(517,*) mscate(iener)*1.e-6, ' MeV => ladder energy'
          do ik1 = 1,NK1
          
          write(517,*) probns(ipart,iener,ik1), ' => No Scat Prob'
          write(517,*) '    amu     -       amu         b        eta'
            do iang = 1,NFIT
             
            end do
            write(517,*) 'CDF for the current region'
            write(517,*) ' 0.0000000E+00'
            do iang = 1,NEXFIT
              write(517,'(1x,1pe14.7,1x)') cumdist(ipart,iener,ik1,iang)
            end do
          end do               !-->  K1 steps
        end do                 !-->  Energy grid
      end do                   !-->  Particle types
      return
      end
