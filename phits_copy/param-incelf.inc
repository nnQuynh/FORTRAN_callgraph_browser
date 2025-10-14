!***********************************************************************
!*                                                                     *
!*      basic parameters and global variables                          *
!*                                                                     *
!***********************************************************************

!-----------------------------------------------------------------------
!   declaration
!-----------------------------------------------------------------------

      integer               :: nnn,lll,element
      character( len = 5  ) :: version
      character( len = 14 ) :: author
      character( len = 18 ) :: fname(19)
      character( len = 18 ) :: aname(10)
      character( len = 18 ) :: rname(6)

!-----------------------------------------------------------------------
!   author and version of SPHINX
!-----------------------------------------------------------------------

      parameter ( author  = 'Hiroki IWAMOTO' ) !author of iwamoto INC
      parameter ( version = '2.06'           ) !version of iwamoto INC

!-----------------------------------------------------------------------
!   Funamental physical and mathematical constant
!-----------------------------------------------------------------------

      parameter( hbar   = 197.326960177d+00     ) !Dirac constant
      parameter( ee     = 137.035999765d+00     ) !fine-structure const
      parameter( rmass  = 938.9265d0            ) !mass of nucleon(MeV)
      parameter( ddmass = 1877.85d0             ) !mass of deuteron(MeV)
      parameter( tmass  = 2816.78d0             ) !mass of triton(MeV)
      parameter( hmass  = 2816.78d0             ) !mass of helion(MeV)
      parameter( amass  = 3727.4d0              ) !mass of alpha (MeV)
      parameter( pmass  = 0.138d0               ) !mass of pion  (GeV)
      parameter( dmass  = 1.232d0               ) !mass of Delta (GeV)
      parameter( smass  = 1.440d0               ) !mass of Nstar (GeV)
      parameter( pi     = 3.14159265358979d+00  ) !circle ratio
      parameter( mxnucl = 104                   ) !number of elements

!-----------------------------------------------------------------------
!   Simulation run constant
!-----------------------------------------------------------------------

      parameter ( mmm    = 100000         )
      parameter ( nnn    = 240            ) !maximum particle number
      parameter ( kkk    = 1000           )
      parameter ( lll    = 1              ) !event number of ground state

!-----------------------------------------------------------------------
!   file number
!-----------------------------------------------------------------------

      parameter ( ninp   = 101               ) !input file
      parameter ( nrdst  = 102, npdst  = 103 ) !spatial & momentum dist.
      parameter ( nio    = 104               ) !initial info.
      parameter ( nout   = 105               ) !INC results

      parameter ( nddxn  = 111, nddxp  = 112 ) !DDX(neutron),  DDX(proton)
      parameter ( nddxd  = 113, nddxt  = 114 ) !DDX(deuteron), DDX(triton)
      parameter ( nddxh  = 115, nddxa  = 116 ) !DDX(helium3),  DDX(alpha)

      parameter ( naddxn = 121, naddxp = 122 ) !DDX(neutron),  DDX(proton)
      parameter ( naddxd = 123, naddxt = 124 ) !DDX(deuteron), DDX(triton)
      parameter ( naddxh = 125, naddxa = 126 ) !DDX(helium3),  DDX(alpha)

      parameter ( nrddxe = 131, nrddxa = 132 ) !DDX(energy),  DDX(angle)

      parameter ( nchgms = 141               ) !charge & mass dist.
      parameter ( ntwopl = 142               ) !2-D charge & mass plot
      parameter ( ntwom  = 143               ) !2-D momentum plot

      parameter ( nchgds = 151, nmasds = 152 ) !charge & mass dist.
      parameter ( nrtwop = 153               ) !2-D charge & mass plot
      parameter ( nrtwom = 154               ) !2-D momentum plot
      parameter ( nremex = 155               ) !remnant vs excitation

!-----------------------------------------------------------------------
!   file name
!-----------------------------------------------------------------------

      data fname/
     &    'input.inp         ',   ! 1 input file name
     &    'spatial.dat       ',   ! 2 output file name of spatial dist.
     &    'momentum.dat      ',   ! 3 output file name of momentum dist.
     &    'inputdata.dat     ',   ! 4 output file name of initial echo
     &    'incresult.dat     ',   ! 5 output file name of INC results
     &    'ddxneutron.dat    ',   ! 6 output file name of DDX(neutron)
     &    'ddxproton.dat     ',   ! 7 output file name of DDX(proton)
     &    'ddxdeuteron.dat   ',   ! 8 output file name of DDX(deuteron)
     &    'ddxtriton.dat     ',   ! 9 output file name of DDX(triton)
     &    'ddxhelium3.dat    ',   !10 output file name of DDX(helium3)
     &    'ddxalpha.dat      ',   !11 output file name of DDX(alpha)
     &    '                  ',   !12 output file name of momentum dist. of emit. part.
     &    '                  ',   !13 output file name of DDX(proton)
     &    '                  ',   !14 output file name of DDX(deuteron)
     &    '                  ',   !15 output file name of DDX(triton)
     &    '                  ',   !16 output file name of DDX(helium3)
     &    '                  ',   !17 output file name of DDX(alpha)
     &    '                  ',   !18 output file name of charge&mass dist.
     &    '                  ' /  !19 output file name of 2-D plot

!-----------------------------------------------------------------------
!   angel file name
!-----------------------------------------------------------------------

      data aname/
     &    'ddxneutron.ang    ',   ! 1 output file name of DDX(neutron)
     &    'ddxproton.ang     ',   ! 2 output file name of DDX(proton)
     &    'ddxdeuteron.ang   ',   ! 3 output file name of DDX(deuteron)
     &    'ddxtriton.ang     ',   ! 4 output file name of DDX(triton)
     &    'ddxhelium3.ang    ',   ! 5 output file name of DDX(helium3)
     &    'ddxalpha.ang      ',   ! 6 output file name of DDX(alpha)
     &    'chargemass.ang    ',   ! 7 output file name of charge&mass dist.
     &    'twodimplot.ang    ',   ! 8 output file name of 2-D plot
     &    'momdist.ang       ',   ! 9 output file name of momentum dist.
     &    'remnant.ang       ' /  !10 output file name of remnant vs excitation

!-----------------------------------------------------------------------
!    element name
!-----------------------------------------------------------------------

      character( len =  3 ) :: elmnt0(mxnucl), elmnt1(mxnucl)

      data elmnt0 /
     &  'H  ','He ','Li ','Be ','B  ','C  ','N  ','O  ','F  ','Ne ',
     &  'Na ','Mg ','Al ','Si ','P  ','S  ','Cl ','Ar ','K  ','Ca ',
     &  'Sc ','Ti ','V  ','Cr ','Mn ','Fe ','Co ','Ni ','Cu ','Zn ',
     &  'Ga ','Ge ','As ','Se ','Br ','Kr ','Rb ','Sr ','Y  ','Zr ',
     &  'Nb ','Mo ','Te ','Ru ','Rh ','Pd ','Ag ','Cd ','In ','Sn ',
     &  'Sb ','Te ','J  ','Xe ','Cs ','Ba ','La ','Ce ','Pr ','Nd ',
     &  'Pm ','Sm ','Eu ','Gd ','Tb ','Dy ','Ho ','Er ','Tm ','Yb ',
     &  'Lu ','Hf ','Ta ','W  ','Re ','Os ','Ir ','Pt ','Au ','Hg ',
     &  'Tl ','Pb ','Bi ','Po ','At ','Rn ','Fr ','Ra ','Ac ','Th ',
     &  'Pa ','U  ','Np ','Pu ','Am ','Cm ','Bk ','Cf ','Es ','Fm ',
     &  'Md ','No ','Lw ','Ku ' /

      data elmnt1 /
     &  'h  ','he ','li ','be ','b  ','c  ','n  ','o  ','f  ','ne ',
     &  'na ','mg ','al ','si ','p  ','s  ','cl ','ar ','k  ','ca ',
     &  'sc ','ti ','v  ','cr ','mn ','fe ','co ','ni ','cu ','zn ',
     &  'ga ','ge ','as ','se ','br ','kr ','rb ','sr ','y  ','zr ',
     &  'nb ','mo ','te ','ru ','rh ','pd ','ag ','cd ','in ','sn ',
     &  'sb ','te ','j  ','xe ','cs ','ba ','la ','ce ','pr ','nd ',
     &  'pm ','sm ','eu ','gd ','tb ','dy ','ho ','er ','tm ','yb ',
     &  'lu ','hf ','ta ','w  ','re ','os ','ir ','pt ','au ','hg ',
     &  'tl ','pb ','bi ','po ','at ','rn ','fr ','ra ','ac ','th ',
     &  'pa ','u  ','np ','pu ','am ','cm ','bk ','cf ','es ','fm ',
     &  'md ','no ','lw ','ku ' /
