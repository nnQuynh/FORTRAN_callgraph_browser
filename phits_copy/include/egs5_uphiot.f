!------------------------------egs5_uphiot.f----------------------------
! Version: 051219-1435
!-----------------------------------------------------------------------
!23456789|123456789|123456789|123456789|123456789|123456789|123456789|12
      real*8 theta,sinthe,costhe,sinphi,cosphi,PI,TWOPI,PI5D2
      common/UPHIOT/     ! Subroutine uphi's input/output with its users
     * theta,                                         ! Scattering angle
     * sinthe,                                   ! sin(scattering angle)
     * costhe,                                   ! cos(scattering angle)
     * sinphi,     ! sin(change in azimuthal direction after scattering)
     * cosphi,     ! cos(change in azimuthal direction after scattering)
     * PI,                                                          ! Pi
     * TWOPI,                                             ! Two times Pi
     * PI5D2                                ! Five times Pi divided by 2
!$OMP THREADPRIVATE(/UPHIOT/)
!-----------------------last line of egs5_uphiot.f---------------------
