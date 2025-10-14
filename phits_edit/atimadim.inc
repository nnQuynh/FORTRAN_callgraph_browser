

      parameter ( nmax=1000, kmax=10 )

      real*8 t1(nmax+kmax), bcoef1(nmax)
      real*8 t2(nmax+kmax), bcoef2(nmax)
      real*8 t3(nmax+kmax), bcoef3(nmax)
      real*8 t4(nmax+kmax), bcoef4(nmax)

      common /sys2/ nb1, kpm1, nb2, kpm2, nb3, kpm3, nb4, kpm4
!$OMP THREADPRIVATE(/sys2/)
      common /sys3/ t1, bcoef1, t2, bcoef2, t3, bcoef3, t4, bcoef4
!$OMP THREADPRIVATE(/sys3/)

