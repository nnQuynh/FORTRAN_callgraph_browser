************************************************************************
      module ELEDATAMOD
*                                                                      *
*     Elemental data for chemistry                                     *
*        created by T.Ogawa : 2020/06/19                               *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              setup elemental data to be shred all                    *
*              over PHITS before openmp parallelization                *
*                                                                      *
*                                                                      *
*        Data to be setup:                                             *
*                                                                      *
*           pot_elem(i,nz)          : ionization potential              *
*           dkin_elem(i,nz)         : kinetic energy                   *
*          elec_elem(i,nz)          : shell electron number            *
*                                  i = 1   K    (1s)                   *
*                                    = 2   L1   (2s)                   *
*                                    = 3   L2   (2p)                   *
*                                    = 4   L3   (2p)                   *
*                                    = 5   M1   (3s)                   *
*                                    = 6   M2   (3p)                   *
*                                    = 7   M3   (3p)                   *
*                                    = 8   M4   (3d)                   *
*                                    = 9   M5   (3d)                   *
*                                    = 10  N1   (4s)                   *
*                                    = 11  N2   (4p)                   *
*                                    = 12  N3   (4p)                   *
*                                    = 13  N4   (4d)                   *
*                                    = 14  N5   (4d)                   *
*                                    = 15  N6   (4f)                   *
*                                    = 16  N7   (4f)                   *
*                                    = 17  O1   (5s)                   *
*                                    = 18  O2   (5p)                   *
*                                    = 19  O3   (5p)                   *
*                                    = 20  O4   (5d)                   *
*                                    = 21  O5   (5d)                   *
*                                    = 22  O6   (5f)                   *
*                                    = 23  O7   (5f)                   *
*                                    = 24  P1   (6s)                   *
*                                    = 25  P2   (6p)                   *
*                                    = 26  P3   (6p)                   *
*                                    = 27  P4,5 (6d)                   *
*                                    = 28  Q1   (7s)                   *
*                                                                      *
*                                                                      *
*                                                                      *
*        Data source:                                                  *
*       1,IUPAC, Compendium of Chemical Terminology, 2nd ed.           *
*             (the "Gold Book") (1975)                                 *
*       2, UCRL--50400 Vol.30                                          *
*         Tables and Graphs of Atomic Subshell and                     *
*         Relaxation Data Derived from the LLNL                        *
*         Evaluated Atomic Data Library (EADL),                        *
*         Z = 1 -100                                                   *
*           S. T. Perkins, D. E. CuUen, M. H. Chen,                    *
*           J. H. Hubbell (National Institute of Standards and         *
*           Technology), J. Rathkopf, J. Scofield                      *
*             Part 2  EADL Atomic Subshell Parameters p2               *
*                                                                      *
*        Material ID number:                                           *
*                                                                      *
*            inorganic materials :  1    : H2O,                        *
*                                   2    : CO2,                        *
*                                   3    : NH2                         *
*                                   4    : NH3                         *
*                                   5    : SF6                         *
*                                   6    : TeF6                        *
*                                                                      *
*            organic materials   : 10001 : CH4,                        *
*                                  10002 : CH3                         *
*                                  10003 : C2H2                        *
*                                  10004 : C2H4                        *
*                                  10005 : C2H6 (ethylene)             *
*                                  10006 : C6H6 (benzen)               *
*                                  10007 : (CH3)2NH (dymethylamine)    *
*                                                                      *
*                                                                      *
*                                                                      *
*                                                                      *
*            pure materials      : 20010 : H2,                         *
*                                  20070 : N2,                         *
*                                  20080 : O2,.....                    *
*                                                                      *
*                                                                      *
*            Future rules                                              *
*             *Understand both chemical formula and name               *
*              such as C2H5OH and Ethernol                             *
*                                                                      *
*             *chem = Co means carbon mono-oxide                       *
*                                                                      *
*             *m1 Co without chem parameter is cobalt                  *
*                                                                      *
*             *currently definition is molar ratios.                   *
*              understand mass ratio in the future                     *
*                                                                      *
*                                                                      *
************************************************************************
      implicit double precision (a-h,o-z)
      include 'err.inc'

      common /paran/ icfn(100), ilfn(100), chfn(100) ! phits paths
      character chfn*200
      character FNAme*200, FNAme1*200, FNAme2*200

      allocatable ichem(:,:), frac(:,:) ! chemical id, fraction  ->  chemical composition. Arguments : Material ID number, components per material
      parameter(maxcomp = 100) ! number of components

      parameter(nelem = 108) ! number of elements
      parameter(mnum_inorg = 8)    ! Inorganic material species
      parameter(mnum_org   = 7)    ! Organic material species
      parameter(mnum_pure  = 90)   ! Pure material species covered.
      parameter(mnum_inident = 3)   ! Materials whose 1st ionisation potential is available

      character chnm_inorg (mnum_inorg)*7,
     &          chnm_org  (mnum_org)*8,
     &          chnm_pure (mnum_pure)*5,
     &          chnm_inident(mnum_inident)*5
      integer   idlen_inorg (mnum_inorg),
     &          idlen_org   (mnum_org),
     &          idlen_pure  (mnum_pure),
     &          idlen_inident(mnum_inident)

      double precision, private :: pot_inident_list(mnum_inident) ! potential of inidentified material

      public ion_potential_table, chemi_form, elec_mol, pot_mol

c ionization potential and electron number of ....
c elements
      double precision, public  :: pot_elem (31,nelem)       ! Shell-wise element ionization potential
      double precision, public  :: dkin_elem (31,nelem)       ! Shell-wise element ionization potential
      double precision,  public  :: elec_elem(31,nelem)       ! Shell-wise number of electrons

c inorganic compounds
      double precision, private :: pot_inorg  (28,mnum_inorg) ! inorganic material ionization potential
      integer,          private :: nelec_inorg(28,mnum_inorg) ! Shell-wise number of electrons

c organic compounds
      double precision, private :: pot_org  (28,mnum_org)     ! organic material ionization potential
      integer,          private :: nelec_org(28,mnum_org)     ! Shell-wise number of electrons

******** EADL Data table on Auger and X-ray ****************************
      parameter(mz =  6) ! 1st element
      parameter(nz =100) ! last element
      parameter(nd = 19) ! variety of vacant level
      parameter(nf = 29) ! variety of de-exciting level
      parameter(nk = 29) ! variety of ejected level

! Length : # of elements * variety of initial vacant levels * variety of final vacant levels * variety of kicked electron levels (0:X-ray emission)
      double precision, private :: str_ejec_s(mz:nz,nd,nf,0:nk)
      double precision, private :: eng_ejec_s(mz:nz,nd,nf,0:nk)
      integer,          private :: isp_ejec_s(mz:nz,nd,nf,0:nk)

! Take the atomic de-excitation table of corresponding part
      double precision, private :: str_ejec(nd,nf,0:nk)
      double precision, private :: eng_ejec(nd,nf,0:nk)
      integer,          private :: isp_ejec(nd,nf,0:nk)
!$OMP THREADPRIVATE(str_ejec, eng_ejec, isp_ejec)

      integer,          private :: levconv(58) ! level conversion table

! Convert : K, L1, L2, L3 ....-> 1,2,3,4.....
      data levconv/
!       1  2  3  4  5  6  7  8  9 10
     &  1, 0, 2, 0, 3, 4, 0, 5, 0, 6,
     &  7, 0, 8, 9, 0,10, 0,11,12, 0,
     & 13,14, 0,15,16, 0,17, 0,18,19,
     &  0,20,21, 0,22,23, 0, 0, 0, 0,
     & 24, 0,25,26, 0,27,28, 0, 0, 0,
     &  0, 0, 0, 0, 0, 0, 0,29/

      private call_atomic_relaxation_table

! Store the result of atomic relaxation
      integer,          public :: lng_rel ! array length
      double precision, public :: eng_rel(100)
      integer,          public :: ktp_rel(100)
!$OMP THREADPRIVATE(lng_rel, eng_rel, ktp_rel)

****************************************************************************

******** Data Taken from Table I(b) . M.E.Rudd, Review of Modern Physics, 64, (1992), 441-490 *******
*  Further data for extention is available ->  https://physics.nist.gov/PhysRefData/Ionization/molTable.html

c     H2O
      data (pot_inorg(i,1),i=   1, 28)/  ! ionization potential
     & 539.7d0, 32.2d0, 18.55d0, 14.73d0, 12.61d0, 23*0.d0/
      data (nelec_inorg(i,1),i=   1, 28)/  ! shell electron
     & 5*2, 23*0/

c     CO2
      data (pot_inorg(i,2),i=   1, 28)/  ! ionization potential
     & 2*541.d0, 297.5d0, 38.6d0, 37.0d0, 19.4d0, 18.08d0, 17.6d0,
     & 13.79d0, 19*0.d0/
      data (nelec_inorg(i,2),i=   1, 28)/  ! shell electron
     & 7*2, 2*4, 19*0/

c     NH2
      data (pot_inorg(i,3),i=   1, 28)/  ! ionization potential
     & 423.8d0, 30.09d0, 17.47d0, 13.02d0, 11.14d0, 23*0.d0/
      data (nelec_inorg(i,3),i=   1, 28)/  ! shell electron
     & 4*2, 1, 23*0/

c     NH3
      data (pot_inorg(i,4),i=   1, 28)/  ! ionization potential
     & 405.6d0, 27.77d0, 16.d0, 10.88d0, 24*0.d0/
      data (nelec_inorg(i,4),i=   1, 28)/  ! shell electron
     & 2*2, 4, 2, 24*0/

c     SF6
      data (pot_inorg(i,5),i=   1, 28)/  ! ionization potential
     & 44.2d0, 41.2d0, 39.3d0, 26.8d0, 22.5d0, 19.9d0, 18.6d0, 16.9d0,
     & 16.9d0, 15.8d0, 18*0.d0/
      data (nelec_inorg(i,5),i=   1, 28)/  ! shell electron
     & 2, 6, 4, 2, 6, 6, 4, 6, 6, 6, 18*0/

c     TeF6
      data (pot_inorg(i,6),i=   1, 28)/  ! ionization potential
     & 47.61d0, 45.85d0, 45.01d0, 27.55d0, 23.05d0, 21.17d0, 19.74d0,
     & 19.46d0, 19.31d0, 18.97d0, 18*0.d0/
      data (nelec_inorg(i,6),i=   1, 28)/  ! shell electron
     & 2, 6, 4, 2, 6, 6, 6, 6, 4, 6, 18*0/

c                *********----------************

c     CH4
      data (pot_org(i,1),i=   1, 28)/  ! ionization potential
     & 290.7d0, 23.d0, 14.35d0, 25*0.d0/
      data (nelec_org(i,1),i=   1, 28)/  ! shell electron
     & 2, 2, 6, 25*0/

c     CH3
      data (pot_org(i,2),i=   1, 28)/  ! ionization potential
     & 305.5d0, 24.57d0, 15.1d0, 9.84d0, 24*0.d0/
      data (nelec_org(i,2),i=   1, 28)/  ! shell electron
     & 2, 2, 4, 1, 24*0/

c     C2H2
      data (pot_org(i,3),i=   1, 28)/  ! ionization potential
     & 2*291.1d0, 23.5d0, 18.38d0, 16.36d0, 11.4d0, 22*0.d0/
      data (nelec_org(i,3),i=   1, 28)/  ! shell electron
     & 5*2, 4, 22*0/

c     C2H4
      data (pot_org(i,4),i=   1, 28)/  ! ionization potential
     & 2*290.9d0, 23.68d0, 19.1d0, 15.87d0, 14.66d0, 12.85d0, 10.51d0,
     & 20*0.d0/
      data (nelec_org(i,4),i=   1, 28)/  ! shell electron
     & 8*2, 20*0/

c     C2H6
      data (pot_org(i,5),i=   1, 28)/  ! ionization potential
     & 2*290.5d0, 23.6d0, 20.16d0, 15.4d0, 13.5d0, 12.36d0, 21*0.d0/
      data (nelec_org(i,5),i=   1, 28)/  ! shell electron
     & 4*2, 4, 2, 4, 21*0/

c     C6H6
      data (pot_org(i,6),i=   1, 28)/  ! ionization potential
     & 4*290.2d0, 26.d0, 22.7d0, 19.0d0, 16.9d0, 15.4d0, 14.8d0,
     & 14.d0, 12.3d0, 11.8d0, 9.24d0, 14*0.d0/
      data (nelec_org(i,6),i=   1, 28)/  ! shell electron
     & 2, 2*4, 2*2, 2*4, 3*2, 4, 2, 4, 4, 14*0/

c     (CH3)2NH
      data (pot_org(i,7),i=   1, 28)/  ! ionization potential
     & 422.7d0, 871.9d0, 305.8d0, 32.62d0, 25.79d0, 23.26d0, 16.7d0,
     & 15.49d0, 15.05d0, 13.85d0, 13.27d0, 12.64d0, 8.94d0, 15*0.d0/
      data (nelec_org(i,7),i=   1, 28)/  ! shell electron
     & 13*2, 15*0/


*****************************************************************************************************

      data (elec_elem( 1,i),i=   1, nelem)/
     & 1.d0, 107*2.d0/

      data (elec_elem( 2,i),i=   1, nelem)/
     & 2*0.d0, 1.d0, 105*2.d0/

      data (elec_elem( 3,i),i=   1, nelem)/
     &4*0.d0, 1.d0, 2.d0, 3.d0, 4.d0, 5.d0, 99*2.d0/ ! correction 2021/10/11 to match with potential table

      data (elec_elem( 4,i),i=   1, nelem)/
     &9*0.d0, 99*4.d0/

      data (elec_elem( 5,i),i=   1, nelem)/
     & 10*0.d0, 1.d0, 97*2.d0/

      data (elec_elem( 6,i),i=   1, nelem)/
     &12*0.d0, 0.33333d0, 0.66667d0, 1.d0, 1.33333d0, 1.66667d0, 91*2.d0
     &/

      data (elec_elem( 7,i),i=   1, nelem)/
     &12*0.d0, 0.66667d0, 1.33333d0, 2.d0, 2.66667d0, 3.33333d0, 91*4.d0
     &/

      data (elec_elem( 8,i),i=   1, nelem)/
     & 20*0.d0, 0.4d0, 0.8d0, 1.2d0, 2.d0, 2.d0, 2.4d0, 2.8d0, 3.2d0,
     & 80*4.d0/

      data (elec_elem( 9,i),i=   1, nelem)/
     & 20*0.d0, 0.6d0, 1.2d0, 1.8d0, 3.d0, 3.d0, 3.6d0, 4.2d0, 4.8d0,
     & 80*6.d0/

      data (elec_elem(10,i),i=   1, nelem)/
     & 18*0.d0, 1.d0, 4*2.d0, 1.d0, 4*2.d0, 1.d0, 79*2.d0/

      data (elec_elem(11,i),i=   1, nelem)/
     &30*0.d0, 0.33333d0, 0.66667d0, 1.d0, 1.33333d0, 1.66667d0, 73*2.d0
     &/

      data (elec_elem(12,i),i=   1, nelem)/
     &30*0.d0, 0.66667d0, 1.33333d0, 2.d0, 2.66667d0, 3.33333d0, 73*4.d0
     &/

      data (elec_elem(13,i),i=   1, nelem)/
     & 38*0.d0, 0.4d0, 0.8d0, 1.6d0, 2.d0, 2.d0, 2.8d0, 3.2d0, 63*4.d0/

      data (elec_elem(14,i),i=   1, nelem)/
     & 38*0.d0, 0.6d0, 1.2d0, 2.4d0, 3.d0, 3.d0, 4.2d0, 4.8d0, 63*6.d0/

      data (elec_elem(15,i),i=   1, nelem)/
     & 57*0.d0,    0.428571d0, 1.285714d0, 1.714286d0, 2.142857d0,
     & 2.571429d0,       3.d0,       3.d0, 3.857143d0, 4.285714d0,
     & 4.714286d0, 5.142857d0, 5.571429d0,    39*6.d0/

      data (elec_elem(16,i),i=   1, nelem)/
     &    57*0.d0, 0.571429d0, 1.714286d0, 2.285714d0, 2.857143d0,
     & 3.428571d0,       4.d0,       4.d0, 5.142857d0, 5.714286d0,
     & 6.285714d0, 6.857143d0, 7.428571d0,    39*8.d0/

      data (elec_elem(17,i),i=   1, nelem)/
     & 36*0.d0, 1.d0, 3*2.d0, 2*1.d0, 2.d0, 2*1.d0, 0.d0, 1.d0, 61*2.d0/

      data (elec_elem(18,i),i=   1, nelem)/
     &48*0.d0, 1.d0, 2.d0, 3.d0, 4.d0, 5.d0, 55*2.d0 ! correction 2021/10/11 to match with potential table
     &/

      data (elec_elem(19,i),i=   1, nelem)/
     &53*0.d0, 55*4.d0 ! correction 2021/10/11 to match with potential table
     &/

      data (elec_elem(20,i),i=   1, nelem)/
     & 56*0.d0, 2*0.4d0, 5*0.d0, 0.4d0, 6*0.d0,   0.4d0, 0.8d0, 1.2d0,
     &   1.6d0,    2.d0,  2.4d0, 2.8d0,  3.6d0, 30*4.d0/

      data (elec_elem(21,i),i=   1, nelem)/
     & 56*0.d0, 2*0.6d0, 5*0.d0, 0.6d0, 6*0.d0,   0.6d0, 1.2d0, 1.8d0,
     &   2.4d0,    3.d0,  3.6d0, 4.2d0,  5.4d0, 30*6.d0/

      data (elec_elem(22,i),i=   1, nelem)/
     &    90*0.d0, 0.857143d0, 1.285714d0, 1.714286d0, 2.571429d0,
     &       3.d0,       3.d0, 3.857143d0, 4.285714d0, 4.714286d0,
     & 5.142857d0, 5.571429d0,     7*6.d0/

      data (elec_elem(23,i),i=   1, nelem)/
     &    90*0.d0, 1.142857d0, 1.714286d0, 2.285714d0, 3.428571d0,
     &       4.d0,       4.d0, 5.142857d0, 5.714286d0, 6.285714d0,
     & 6.857143d0, 7.428571d0,     7*8.d0/

      data (elec_elem(24,i),i=   1, nelem)/
     & 54*0.d0, 1.d0, 22*2.d0, 2*1.d0, 29*2.d0/

      data (elec_elem(25,i),i=   1, nelem)/  ! not certain above 96th (Cm)
     &80*0.d0, 0.33333d0, 0.66667d0, 1.d0, 1.33333d0, 1.66667d0, 23*2.d0
     &/

      data (elec_elem(26,i),i=   1, nelem)/
     &80*0.d0, 0.66667d0, 1.33333d0, 2.d0, 2.66667d0, 3.33333d0, 23*4.d0
     &/

      data (elec_elem(27,i),i=   1, nelem)/  ! not certaint above 96th (Cm)
     & 88*0.d0, 0.4d0, 0.8d0, 3*0.4d0, 2*0.d0, 0.4d0, 7*0.d0, 0.8d0,
     &   1.2d0, 1.6d0,  2.d0,   2.4d0/

      data (elec_elem(28,i),i=   1, nelem)/  ! not certaint above 96th (Cm)
     & 88*0.d0, 0.6d0, 1.2d0, 3*0.6d0, 2*0.d0, 0.6d0, 7*0.d0, 1.2d0,
     &   1.8d0, 2.4d0,  3.d0,   3.6d0/

      data (elec_elem(29,i),i=   1, nelem)/  ! not certaint above 96th (Cm)
     & 86*0.d0, 1.d0, 21*2.d0/

      data (elec_elem(30,i),i=   1, nelem)/  ! not certaint above 96th (Cm)
     & 102*0.d0, 0.33333d0, 5*0.d0/

      data (elec_elem(31,i),i=   1, nelem)/  ! not certaint above 96th (Cm)
     & 102*0.d0, 0.66667d0, 5*0.d0/

!  Name data of chemicals

      data (chnm_inorg(i),i=1,mnum_inorg)/
     & 'h2o    ','co2    ','nh2    ','nh3    ','sf6    ','tef6   ',
     & 'water  ','ammonia'/

      data (idlen_inorg(i),i=1,mnum_inorg)/
     &  3,      3,      3,      3,      3,      4,      5,      7/

      data (chnm_org(i),i=1,mnum_org)/
     & 'ch4     ','ch3     ','c2h2    ','c2h4    ','c2h6    ','c6h6    '
     &,'(CH3)2NH'/

      data (idlen_org(i),i=1,mnum_org)/
     & 3,      3,      4,      4,      4,      4,      8/

      data (chnm_pure(i),i=1,mnum_pure)/ ! Table dedicated to allotropes such as ozone (O3)
     & '     ','     ','     ','     ','     ','     ','     ','     ', ! 1-9
     & '     ',
     & 'h2   ','     ','     ','     ','     ','     ','     ','     ', ! 10-19
     & '     ','     ',
     & 'he   ','     ','     ','     ','     ','     ','     ','     ', ! 20-29
     & '     ','     ',
     & 'li   ','     ','     ','     ','     ','     ','     ','     ', ! 30-39
     & '     ','     ',
     & 'be   ','     ','     ','     ','     ','     ','     ','     ', ! 40-49
     & '     ','     ',
     & 'b    ','     ','     ','     ','     ','     ','     ','     ', ! 50-59
     & '     ','     ',
     & 'c    ','     ','     ','     ','     ','     ','     ','     ', ! 60-69
     & '     ','     ',
     & 'n2   ','     ','     ','     ','     ','     ','     ','     ', ! 70-79
     & '     ','     ',
     & 'o2   ','     ','     ','     ','     ','     ','     ','     ', ! 80-90
     & '     ','     ','     '
     &/

      data (idlen_pure(i),i=1,mnum_pure)/
     &  0,      0,      0,      0,      0,      0,      0,      0,
     &  0,
     &  2,      0,      0,      0,      0,      0,      0,      0,
     &  0,      0,
     &  2,      0,      0,      0,      0,      0,      0,      0,
     &  0,      0,
     &  2,      0,      0,      0,      0,      0,      0,      0,
     &  0,      0,
     &  2,      0,      0,      0,      0,      0,      0,      0,
     &  0,      0,
     &  1,      0,      0,      0,      0,      0,      0,      0,
     &  0,      0,
     &  1,      0,      0,      0,      0,      0,      0,      0,
     &  0,      0,
     &  2,      0,      0,      0,      0,      0,      0,      0,
     &  0,      0,
     &  2,      0,      0,      0,      0,      0,      0,      0,
     &  0,      0,      0/

* Element table
      character element(118)*2,elmnt(118)*2
      integer   elmntlng(118)

      data element /
     &    'h ','he','li','be','b ','c ','n ','o ','f ','ne',
     &    'na','mg','al','si','p ','s ','cl','ar','k ','ca',
     &    'sc','ti','v ','cr','mn','fe','co','ni','cu','zn',
     &    'ga','ge','as','se','br','kr','rb','sr','y ','zr',
     &    'nb','mo','tc','ru','rh','pd','ag','cd','in','sn',
     &    'sb','te','i ','xe','cs','ba','la','ce','pr','nd',
     &    'pm','sm','eu','gd','tb','dy','ho','er','tm','yb',
     &    'lu','hf','ta','w ','re','os','ir','pt','au','hg',
     &    'tl','pb','bi','po','at','rn','fr','ra','ac','th',
     &    'pa','u ','np','pu','am','cm','bk','cf','es','fm',
     &    'md','no','lr','rf','db','sg','bh','hs','mt','ds',
     &    'rg','cn','nh','fl','mc','lv','ts','og'/

      data elmnt /
     &    'H ','He','Li','Be','B ','C ','N ','O ','F ','Ne',
     &    'Na','Mg','Al','Si','P ','S ','Cl','Ar','K ','Ca',
     &    'Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn',
     &    'Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y ','Zr',
     &    'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn',
     &    'Sb','Te','I ','Xe','Cs','Ba','La','Ce','Pr','Nd',
     &    'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb',
     &    'Lu','Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg',
     &    'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th',
     &    'Pa','U ','Np','Pu','Am','Cm','Bk','Cf','Es','Fm',
     &    'Md','No','Lr','Rf','Db','Sg','Bh','Hs','Mt','Ds',
     &    'Rg','Cn','Nh','Fl','Mc','Lv','Ts','Og'/

      data elmntlng /
     &       1,           3*2,                     5*1,
     &                    5*2,      2*1,      2*2,   1,
     &          3*2,   1,
     &                                       15*2,   1,
     &         13*2,   1,
     &              20*2,  1,
     &    17*2,   1,
     &    26*2/

***************************************************************************
*                                                                         *
*       Guidance on defining new materials                                *
*           1 : add name and length to 'chnm_inident' and 'idlen_inident' *
*           2 : add outermost potential to array 'pot_inident_list'       *
*           3 : add case to pot_inident and nelec_inident. Inner shells   *
*                are inherited from original elements. The outermost shell*
*                is replaced by pot_inident_list                          *
*                                                                         *
*          2023/9/12 potential -> band gap energy. Electrons above this   *
*                                                 become mobile.          *
*                                                                         *
*          Attention!  Don't forget to change pot_inident and             *
*                      nelec_inident when you change it                   *
*                                                                         *
*          Ionization potential can be found from here                    *
*              https://webbook.nist.gov/chemistry/ie-ser/                 *
***************************************************************************

      data (chnm_inident(i),i=1,mnum_inident)/
     & 'sio2 ','lif  ','si   '/

      data (idlen_inident(i),i=1,mnum_inident)/
     & 4,      3,      2/

      data (pot_inident_list(i),i=1,mnum_inident)/
     & 8.95d0, 12.0d0, 1.13d0/

      contains


*********  Functions/Subroutines to setup data **************************************

************************************************************************
*                                                                      *
      subroutine allocate_frac_ichem
*                                                                      *
*        Last Revised:     2020/07/04                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              Allocate frac and ichem                                 *
*                                                                      *
*                                                                      *
************************************************************************
      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      allocate(ichem(1000,maxcomp),frac(1000,maxcomp)) ! maxmat < 1000 assumed. Components in one material < 100 assumed.
      frac  = 0.d0
      ichem = 0

      return
      end subroutine


************************************************************************
*                                                                      *
      subroutine ion_potential_table(ierr)
*                                                                      *
*        Last Revised:     2020/06/22                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to read ionization energy table from external file      *
*                                                                      *
*                                                                      *
************************************************************************
      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------


*-----------------------------------------------------------------------
*           Normal Nucleus data from block data
*-----------------------------------------------------------------------
      FNAme = chfn(1)(1:ilfn(1))//'/data/'//'ion_potential.dat'

      OPEN(3000, FILE = FNAme, STATUS = 'old',ERR=201)

      i = 0

      READ(3000, *, END = 202)
      READ(3000, *, END = 202)

      do
        i = i + 1

        READ(3000, *, END = 200) pot_elem(1:28,i)

      end do

 200  CLOSE(3000)
      return
*-----------------------------------------------------------------------

 201  write(*,*) 'ion_potential.dat is missing in data folder'
      ErrID = 'L:587/R:ion_potential_table/F:elemdatamod.f'
      ErrCha = ''
      call ErrWrite(ErrID,ErrCha)
      ierr = -1
      return
 202  write(*,*) 'ion_potential.dat is empty'
      ErrID = 'L:593/R:ion_potential_table/F:elemdatamod.f'
      ErrCha = ''
      call ErrWrite(ErrID,ErrCha)
      ierr = -1
      return
      end subroutine

************************************************************************
*                                                                      *
      subroutine electron_kinetic_energy_table
*                                                                      *
*        Last Revised:     2022/10/26     by Yuho Hirata               *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to read ionization energy table from external file      *
*                                                                      *
*                                                                      *
************************************************************************
      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------


*-----------------------------------------------------------------------
*           Normal Nucleus data from block data
*-----------------------------------------------------------------------
      FNAme = chfn(1)(1:ilfn(1))//'/data/'//
     &'subshell_kinetic_energy.dat'

      OPEN(100, FILE = FNAme, STATUS = 'old',ERR=201)

      i = 0

      READ(100, *, END = 202)
      READ(100, *, END = 202)

      do
        i = i + 1

        READ(100, *, END = 200) dkin_elem(1:28,i)
      end do

 200  CLOSE(100)
      return
*-----------------------------------------------------------------------

 201  write(*,*) 'subshell_kinetic_energy.dat
     & is missing in data folder'
      return
 202  write(*,*) 'subshell_kinetic_energy.dat
     & is empty'
      return
      end subroutine
************************************************************************
*                                                                      *
      subroutine set_atomic_relaxation_table(ierr)
*                                                                      *
*        Last Revised:     2021/09/08                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to read Auger and X-ray emission table                  *
*              Data source is UCRL~50400-Vol.30 by S. T. Perkins et al,*
*                                                                      *
*              Tables and Graphs of Atomic Subshell and                *
*              Relaxation Data Derived from the LLNL                   *
*              Evaluated Atomic Data Library (EADL),                   *
*              Z = 1 -100                                              *
*              https://www-nds.iaea.org/epdl97/libsall.htm             *
************************************************************************
      implicit double precision (a-h,o-z)

      FNAme1=chfn(1)(1:ilfn(1))//'/data/'//'Relaxation-NonRadiative.dat'
      FNAme2=chfn(1)(1:ilfn(1))//'/data/'//'Relaxation-Radiative.dat'

      OPEN(3001, FILE = FNAme1, STATUS = 'old',ERR=201)
      OPEN(3002, FILE = FNAme2, STATUS = 'old',ERR=201)

      str_ejec_s = 0.d0
      eng_ejec_s = 0.d0

      do
        READ(3001, *, END = 200) lz,ld,lf,lk,str,eng ! element Z, Vacant level, Fall electron level, Kicked electron level, transition intensity, energy
        str_ejec_s(lz,levconv(ld),levconv(lf),levconv(lk)) = str
        eng_ejec_s(lz,levconv(ld),levconv(lf),levconv(lk)) = eng*1.d-6
      end do

 200  CLOSE(3001)

      do
        READ(3002, *, END = 300) lz,ld,lf,str,eng
        str_ejec_s(lz,levconv(ld),levconv(lf),0) = str
        eng_ejec_s(lz,levconv(ld),levconv(lf),0) = eng*1.d-6
      end do

 300  CLOSE(3002)

      return
*-----------------------------------------------------------------------

 201  write(*,*)'Relaxation-NonRadiative.dat or Relaxation-Radiative.dat
     & is missing in data folder'
      ErrID = 'L:696/R:set_atomic_relaxation_table/F:elemdatamod.f'
      ErrCha = ''
      call ErrWrite(ErrID,ErrCha)
      ierr = -1
      return
      end subroutine


*********  Functions/Subroutines to call data **************************************

************************************************************************
*                                                                      *
      subroutine call_atomic_relaxation_table(iz)
*                                                                      *
*        Last Revised:     2021/09/08                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to call Auger and X-ray emission table                  *
*                                                                      *
*       input                                                          *
*        iz  :  atomic number                                          *
*                                                                      *
*       output                                                         *
*                                                                      *
*                                                                      *
*                                                                      *
************************************************************************
      implicit double precision (a-h,o-z)

      integer, save :: izrem
      data izrem /0/
!$OMP THREADPRIVATE(izrem)

      if(iz .eq. izrem) return

      izrem = iz

*-----------------------------------------------------------------------

        str_ejec(1:nd,1:nf,0:nk) = str_ejec_s(iz,1:nd,1:nf,0:nk)
        eng_ejec(1:nd,1:nf,0:nk) = eng_ejec_s(iz,1:nd,1:nf,0:nk)

      return
      end subroutine

************************************************************************
*                                                                      *
      subroutine chemi_form(chemID,len,j,ierr)
*                                                                      *
*       Read chemical formula defined in [Track structure]             *
*                                                                      *
*       Data source: M.E.Rudd, Review of Modern Physics, 64 (1992)     *
*                    441-491                                           *
*                                                                      *
*     input  :                                                         *
*                                                                      *
*          chemID : Chemical form input                                *
*          len    : Length of chemID                                   *
*                                                                      *
*     output :                                                         *
*                                                                      *
*           ierr  : error flag                                         *
*             j   : Chemical form ID.                                  *
*                                                                      *
************************************************************************

      implicit double precision (a-h,o-z)

      character chemID*200

      ierr = 0
      j    = 99999999

      if(abs(j) .le. 1) return

      do i = 1, mnum_inorg
        if( chemID(1:len) .eq. chnm_inorg(i)(1:idlen_inorg(i)) ) then
          j = i
          return
        endif
      enddo

      do i = 1, mnum_org
        if( chemID(1:len) .eq. chnm_org(i)(1:idlen_org(i)) ) then
          j = i + 10000
          return
        endif
      enddo

      do i = 1, mnum_pure
        if( chemID(1:len) .eq. chnm_pure(i)(1:idlen_pure(i)) ) then
          j = i + 20000
          return
        endif
      enddo

      do i = 1, nelem
        if( chemID(1:len) .eq. element(i)(1:elmntlng(i)) ) then
          j = i + 22000
          return
        endif
      enddo

      do i = 1, mnum_inident
        if( chemID(1:len) .eq. chnm_inident(i)(1:idlen_inident(i))) then
          j = i + 30000
          return
        endif
      enddo

      if(j .eq. 99999999) ierr = 1 ! agree with nothing!

      return
      end subroutine


*-----------------------------------------------------------------------
************************************************************************
*                                                                      *
      subroutine form_chemi(chemID,len,j,ierr)
*                                                                      *
*       Read chemical ID and return chemical formula                   *
*       Inverse of chemi_form                                          *
*                                                                      *
*     input  :                                                         *
*                                                                      *
*             j   : Chemical form ID.                                  *
*                                                                      *
*     output :                                                         *
*                                                                      *
*          chemID : Chemical form input                                *
*          len    : Length of chemID                                   *
*           ierr  : error flag                                         *
*                                                                      *
************************************************************************

      implicit double precision (a-h,o-z)

      character chemID*200

      ierr = 0

      i = mod(j, 10000)

      if(    j .gt. 30000 .and. i .le. mnum_inident) then
         chemID = chnm_inident(i)(1:idlen_inident(i))
         len    = idlen_inident(i)
         return
      elseif(j .gt. 20000 .and. i .le. mnum_pure) then
         chemID = chnm_pure(i)(1:idlen_pure(i))
         len    = idlen_pure(i)
         return
      elseif(j .gt. 20000 .and. i .le. 2000+nelem) then
         chemID = element(i-2000)(1:elmntlng(i-2000))
         len    = elmntlng(i-2000)
         return
      elseif(j .gt. 10000 .and. i .le. mnum_org) then
         chemID = chnm_org(i)(1:idlen_org(i))
         len    = idlen_org(i)
         return
      elseif(i .le. mnum_inorg) then
         chemID = chnm_inorg(i)(1:idlen_inorg(i))
         len    = idlen_inorg(i)
         return
      elseif(i .le. mnum_inident) then
         chemID = chnm_inident(i)(1:idlen_inident(i))
         len    = idlen_inident(i)
         return
      endif

      ierr = 1 ! agree with nothing!

      return
      end subroutine



************************************************************************
*                                                                      *
      function pot_mol(ishell,ichemID)
*                                                                      *
*       Ionization potential of chemical compounds and elements        *
*                                                                      *
*       Data source:                                                   *
*                                                                      *
*                                                                      *
*     input  :                                                         *
*                                                                      *
*         ichemID  : Chemical compound ID                              *
*         ishell   : Shell ID  (1-28 assumed)                          *
*                                                                      *
*     output :                                                         *
*                                                                      *
*          pot_mol : Ionization potential                              *
*                                                                      *
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      SelectCase(ichemID)
         Case(    1:10000)
           pot_mol = pot_inorg(ishell, ichemID)
         Case(10001:20000)
           pot_mol = pot_org(ishell, ichemID-10000)
         Case(20001:22000)
           pot_mol = pot_elem(ishell, (ichemID-20000)/10 )
         Case(22001:30000)
           pot_mol = pot_elem(ishell, (ichemID-22000) )
         Case(30001:40000)
           pot_mol = pot_inident(ishell, ichemID-30000 )
      end Select

      return
      end function

************************************************************************
*                                                                      *
      function elec_mol(ishell,ichemID)
*                                                                      *
*       Electron number in shells                                      *
*                                                                      *
*                                                                      *
*     input  :                                                         *
*                                                                      *
*         ichemID  : Chemical compound ID                              *
*         ishell   : Shell ID  (1-28 assumed)                          *
*                                                                      *
*     output :                                                         *
*                                                                      *
*        elec_mol : Electron number in shells                         *
*                                                                      *
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      elec_mol=0  ! T.Sato 2022/11/21

      SelectCase(ichemID)
         Case(    1:10000)
           elec_mol = dble(nelec_inorg(ishell, ichemID))
         Case(10001:20000)
           elec_mol = dble(nelec_org(ishell, ichemID-10000))
         Case(20001:22000)
           elec_mol = elec_elem(ishell, (ichemID-20000)/10 )
         Case(22001:30000)
           elec_mol = elec_elem(ishell, ichemID-22000 )
         Case(30001:40000)
           elec_mol = dble(nelec_inident(ishell, ichemID-30000))
      end Select

      return
      end function


************************************************************************
*                                                                      *
      function pot_inident(ishell,ichemID)
*                                                                      *
*       Ionization potential of partially-defined chemical compounds   *
*                                                                      *
*       Note : if more chemicals are defined and this function become  *
*              too long, consider to define as an external file        *
*                                                                      *
*                                                                      *
*                                                                      *
*  Instruction on how to define new materials :                        *
*   1, Open ion_potential.dat                                          *
*   2, Look up the constituent element                                 *
*   3, First group of orbits is                                        *
*         (the number of occupied shells of first element) - 1         *
*                  (The outermost shell -1 goes to bond shell          *
*   4, Second group of orbits is                                       *
*         (the number of occupied shells of second element) - 1        *
*   5, Take bond energy for the bond shell, which is the outermost     *
*                                                                      *
*                                                                      *
*     input  :                                                         *
*                                                                      *
*         ichemID  : Chemical compound ID                              *
*         ishell   : Shell ID  (1-28 assumed)                          *
*                                                                      *
*     output :                                                         *
*                                                                      *
*          pot_mol : Ionization potential                              *
*                                                                      *
*                                                                      *
*                                                                      *
************************************************************************
      implicit real*8 (a-h,o-z)

      Select Case(ichemID)
        case(1) ! SiO2
          Select Case(ishell)
            case(1:2) ! O shell
                pot_inident = pot_elem(ishell, 8)
            case(3:8) ! Si shell
                pot_inident = pot_elem(ishell-2, 14)
            case(9)   ! Bond shell
                pot_inident = pot_inident_list(ichemID)
          end Select

        case(2) ! LiF
          Select Case(ishell)
            case(1) ! Li shell
                pot_inident = pot_elem(ishell, 3)
            case(2:3) ! F shell
                pot_inident = pot_elem(ishell-1, 9)
            case(4)   ! Bond shell
                pot_inident = pot_inident_list(ichemID)
          end Select

      end Select


      return
      end function

************************************************************************
*                                                                      *
      function nelec_inident(ishell,ichemID)
*                                                                      *
*       Electron number in shells                                      *
*                                                                      *
*                                                                      *
*     input  :                                                         *
*                                                                      *
*         ichemID  : Chemical compound ID                              *
*         ishell   : Shell ID  (1-28 assumed)                          *
*                                                                      *
*     output :                                                         *
*                                                                      *
*        elec_mol : Electron number in shells                         *
*                                                                      *
*                                                                      *
*                                                                      *
************************************************************************
      implicit real*8 (a-h,o-z)

      nelec_inident = 0

      Select Case(ichemID)
        case(1) ! SiO2
          Select Case(ishell)
            case(1:2) ! O shell
                nelec_inident = elec_elem(ishell, 8)
            case(3:8) ! Si shell
                nelec_inident = elec_elem(ishell-2, 14)
            case(9)   ! Bond shell
                nelec_inident = elec_elem(3, 8)+elec_elem(7, 14)
          end Select

        case(2) ! LiF
          Select Case(ishell)
            case(1:2) ! Li shell
                nelec_inident = elec_elem(ishell, 3)
            case(3:8) ! F shell
                nelec_inident = elec_elem(ishell-2, 9)
            case(9)   ! Bond shell
                nelec_inident = elec_elem(3, 8)+elec_elem(7, 9)
          end Select

      end Select

      return
      end function

*-----------------------------------------------------------------------

*********  Functions/Subroutines to use data **************************

************************************************************************
*                                                                      *
      subroutine do_atom_relax(iz,ld)
*                                                                      *
*        Last Revised:     2021/09/08                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              perform Auger and X-ray emission reaction               *
*                                                                      *
*       input                                                          *
*        iz  :  atomic number                                          *
*        ld  :  initial vacant level                                   *
*                                                                      *
*       output                                                         *
*        lng_rel :  length of array                                    *
*        eng_rel :  ejectile energy (MeV)                              *
*        ktp_rel :  ejectile species (11:electron, 22:photon)          *
*                                                                      *
************************************************************************
      implicit double precision (a-h,o-z)
      include 'err.inc'

      integer, save :: ierr
      data ierr /0/
!$OMP THREADPRIVATE(ierr)

      integer ivac(29) ! vacancy
*-----------------------------------------------------------------------
      lng_rel  = 0

      if(iz .lt. 6) return

      call call_atomic_relaxation_table(iz)

      eng_rel  = 0.d0
      ktp_rel  = 0
      ivac     = 0
      ivac(ld) = 1

      do k = 1, 10000

! find vacant level
       do i = ld, 29 ! no vacancy below initial vacancy
         if(ivac(i) .gt. 0) exit
       enddo

       if(i.ge.20) return ! 2023/4/7 Ogawa. No info in EADL beyond O5 shell
       lvac = i

       r = unirn(dummy)

! find deexcitation path
       do i = lvac+1, nf
        r = r - str_ejec(lvac,i,0) ! X-ray emission
        if(r .lt. 0.d0) then
            j = 0
            goto 100
        endif
        do j = i, nk               ! Auger electron emission
         r = r - str_ejec(lvac,i,j)
         if(r .lt. 0.d0) goto 100
        enddo
       enddo

       return ! no more relaxation probability

! relaxation following the path
 100   lng_rel = lng_rel + 1
       eng_rel(lng_rel) = eng_ejec(lvac,i,j)
       ivac(lvac) = ivac(lvac) - 1
       ivac(i)    = ivac(i)    + 1
       if(j .eq. 0) then
         ktp_rel(lng_rel) = 22
       else
         ivac(j) = ivac(j) + 1
         ktp_rel(lng_rel) = 11
       endif

      enddo

      if( k .ge. 10000 .and. ierr .lt. 5) then
      write(ErrCha,*) 'Strange atomic relaxation. 10000 attempts. Charg
     &e and vacant level are ', iz, ld
      ErrID = 'L:1155/R:do_atom_relax/F:elemdatamod.f'
      call ErrWrite(ErrID,ErrCha)
      ierr = ierr + 1
      endif

      end subroutine




*-----------------------------------------------------------------------

      end module
