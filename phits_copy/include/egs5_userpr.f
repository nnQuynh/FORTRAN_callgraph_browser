!-----------------------------egs5_userpr.f-----------------------------
! Version: 051219-1435
!-----------------------------------------------------------------------
!23456789|123456789|123456789|123456789|123456789|123456789|123456789|12


      real*8
     z b
     z ,blc
     z ,omega0
     z ,ems
     z ,bms
     z ,gms
      integer iskpms
c
      common/USERPR/                                       ! User-MS
     * b,
     * blc,
     * omega0,
     * ems,
     * bms,
     * gms,
     * iskpms
!$OMP THREADPRIVATE ( /USERPR/ )





!-----------------------last line of egs5_userpr.f----------------------
