!----------------------------egs5_stack.f-------------------------------
! Version: 051219-1435
!          080425-1100
!-----------------------------------------------------------------------
!23456789|123456789|123456789|123456789|123456789|123456789|123456789|12

      real*8 e,x,y,z,u,v,w,uf,vf,wf
      real*8 dnear,wt,k1step,k1rsd,k1init,time
      real*8 deinitial,deresid,denstep
      real*8 ecapbind ! T.Sato 2015/10/12 for save EII capbind
      integer  iq,ir,latch,np,latchi

      common/STACK/           ! Information kept about current particles
     * e(MXSTACK),      ! Total energy of particle (including rest mass)
     * x(MXSTACK),                                          ! X-position
     * y(MXSTACK),                                          ! Y-position
     * z(MXSTACK),                                          ! Z-position
     * u(MXSTACK),                             ! X-axis direction cosine
     * v(MXSTACK),                             ! Y-axis direction cosine
     * w(MXSTACK),                             ! Z-axis direction cosine
     * uf(MXSTACK),         ! Electric field vectors of polarized photon
     * vf(MXSTACK),
     * wf(MXSTACK),
     * dnear(MXSTACK),          ! Estimated distance to nearest boundary
     * wt(MXSTACK),                                    ! Particle weight
     * k1step(MXSTACK),      ! Scat stren to next hinge
     * k1rsd(MXSTACK),       ! Scat stren from hinge to end of step
     * k1init(MXSTACK),      ! Scat of prev hinge end of step
     * time(MXSTACK),
     * deinitial,
     * deresid,
     * denstep,
     * iq(MXSTACK),        ! Particle charge, -1(e-), 0(photons), +1(e+)
     * ir(MXSTACK),                                      ! Region number
     * latch(MXSTACK),                               ! Latching variable
     * np,                                         ! Stack pointer index
     * latchi                        ! Initialization for latch variable
      common/STACK2/ecapbind           ! T.Sato 2015/10/12, EII capbind energy
!$OMP THREADPRIVATE(/STACK/)
!$OMP THREADPRIVATE(/STACK2/)
!------------------------last line of egs5_stack.f----------------------

