
      subroutine egs5ede(irin,ityp,ein,eout,nbeta)
c     treats dE of electron and positron
      use EGS5_MS_MOD !FURUTA20140825
      USE egs5_media_mod !<- 2015.08xx allocatable
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

      real*8
     * ustep0,kinit0,ktotal,detot,scpow,thard,tmscat,tinel,
     * hardstep,ecsda,k1s0

      real*8 EPSEMFP,ENEPS                            ! Local parameters

      data
     * EPSEMFP/1.D-12/,                    ! Smallest electron mfp value
     * ENEPS/0.0005D0/,     ! Difference between ecut and end-point energy
     * ierust/0/

      real *8 fpl,ein,eout,k1i,k1r,k1s
      integer icl,ich,mk,irin,nbeta,ityp,ncol
      common /egs5cmn3/thard,tinel,tmscat,hardstep,sig,scpow,dedx,sig0,
     $                 ams
!$OMP THREADPRIVATE(/egs5cmn3/)
      common /egs5cmn4/k1i,k1r,k1s
!$OMP THREADPRIVATE(/egs5cmn4/)
      real*8 esmax,esmin,emin
      common /eparm/  esmax, esmin, emin(20)
      real*8 eadd
      parameter(eadd=1d-15)    ! additional extra energy of 1eV
                              ! to transport electrons with e=emin
c
         np = 1
        eke = ein

        eie = eke+RM
       eout = eke      ! ein = eout for non event

        irl = irin
      lelec = iq(np)

      k1init(np) = k1i
      k1rsd(np)  = k1r
      k1step(np) = k1s

c

      if(tinel .eq. 0.0) then

         de = deinitial + deresid

         ekef = eke - de
         eold = eie
         enew = eold - de


         eie = eie - de
         e(np) = eie
         eout = eke -de

! special treatment for PHITS
         if( e(np) .le. ( ecut(irl) + eadd) ) then ! ecut

            if(deresid.eq.0.d0) then
               eout = emin(ityp) ! set energy PHITS e-cut
               thard = 0d0
               tinel = 0d0
               tmscat = 0d0
               nbeta = 3




            else

               nbeta = 2
               ! solve later (dE treatment) solved by eadd ?

               eout = emin(ityp) + eadd

               denstep = deresid


               deinitial = 0.d0
               deresid = 0.d0

               k1i = k1init(np)
               k1r = k1rsd(np)
               k1s = k1step(np)

            end if

         endif



!         eout = eout - RM       ! kinetic energy for PHITS
!         treat at each sentence


      endif

      return
      end
