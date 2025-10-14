      module mod_ompparallel
      integer ipomp,npomp
      common /ipomp0/ipomp,npomp
!$OMP THREADPRIVATE(/ipomp0/)
      end module
