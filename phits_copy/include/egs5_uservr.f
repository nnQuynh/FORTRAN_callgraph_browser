!-----------------------------egs5_uservr.f-----------------------------
! Version: 051219-1435
! Reference:  Note, only cexptr is active in this release
!-----------------------------------------------------------------------
!23456789|123456789|123456789|123456789|123456789|123456789|123456789|12
      real*8 cexptr
      common/USERVR/                           ! User-Variance-Reduction
     * cexptr     ! c variable for exponential pathlength transformation
!$OMP THREADPRIVATE ( /USERVR/ )

!-----------------------last line of egs5_uservr.f----------------------
c memo
c smp loop  .  ref only "cexptr"  in egs5pfpl .
c memo

