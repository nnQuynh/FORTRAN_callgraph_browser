*    projectil:

      integer           zp
      double precision  mp

*    target:

      integer    maxnuc
      parameter (maxnuc=300) !FURUTA20201203
      integer           nnuc         ! no. of nucleons
      integer           zt(maxnuc)
      double precision  mt(maxnuc)
      double precision  pot(maxnuc)
      double precision  anuc(maxnuc)
      integer           gas          ! 1=gas, 0=solid
      double precision  rho
      double precision  fntp

      integer           dedxmat

      common /sys1/ mp, mt, pot, anuc, rho, fntp, zp, nnuc, zt, gas
     &              ,dedxmat
!$OMP THREADPRIVATE(/sys1/)


