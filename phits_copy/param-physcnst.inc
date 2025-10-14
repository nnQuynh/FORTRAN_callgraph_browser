************************************************************************
*                                                                      *
*     'param-physcnst.inc'                                             *
*        CODATA Internationally recommended values of                  *
*     fundamental physical constants obtained from                     *
*     http://physics.nist.gov/cuu/Constants/index.html                 *
*                                                                      *
************************************************************************

*----------------------------------------------------------------------*

      real*8 physc(9), rstms(0:20)

* physical constants
      data physc/
     1	   3.141592653589793d0, ! pi
     2	   2.99792458d8,        ! speed of light in vacuum in m/s
     3	   197.3269718d0,       ! hbar c in MeV fm
     4	   137.035999074d0,     ! 1/alpha (inverse fine-structure constant)
     5	   0.931494061d0,       ! atomic mass unit in GeV/(c^2)
     6	   8.854187817d-12,     ! electric constant (epsilon_0) in F/m
     7	   8.6173324d-5,        ! Boltzmann constant (k_B) in eV/K
     8	   6.02214129d23,       ! Avogadro constant (N_A) in 1/mol
     9	   1.602176565d-19      ! elementary charge in C
     &     /

* rest mass in GeV/(c^2)
* (PDG2014: Particle data group 2014 edition [Chinese Physics C38 (2014) 090001])
      data rstms/
     &     0.931494061d0,    ! a.m.u.
     1     0.938272046d0,    ! proton
     2     0.939565379d0,    ! neutron
     3     0.13957018d0,     ! pion+ (from PDG2014)
     4     0.1349766d0,      ! pion0 (from PDG2014)
     5     0.13957018d0,     ! pion- (the same as pion+)
     6     0.1056583715d0,   ! muon+
     7     0.1056583715d0,   ! muon- (the same as muon+)
     8     0.493677d0,       ! kaon+ (from PDG2014)
     9     0.497614d0,       ! kaon0 (from PDG2014)
     &     0.493677d0,       ! kaon- (the same as kaon+)
     1     0.d0,             ! other particles such as hadrons
     2     0.000510998928d0, ! electron
     3     0.000510998928d0, ! positron (the same as electron)
     4     0d0,              ! photon
     5     1.875612859d0,    ! deuteron
     6     2.808921005d0,    ! triton
     7     2.808391482d0,    ! 3He
     8     3.727379240d0,    ! alpha (4He)
     9     0d0,              ! nucleus
     &     0d0               ! all
     &     /

*-----------------------------------------------------------------------
