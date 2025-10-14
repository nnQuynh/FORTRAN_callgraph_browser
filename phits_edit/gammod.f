************************************************************************
*                                                                      *
      module levdat
*                                                                      *
*       Give levels, gamma energy, transition probability, etc.        *
*       By deleting "comp", and disactivating                          *
*       "! memory wasting but speedy"                                  *
************************************************************************

      implicit double precision(a-h, o-z)

      private :: levcall, levread, allocRIPL
      public  :: levset, levUNset
      
      
* Level energy
      integer                       :: nlevel
      double precision, allocatable :: elevel(:)
      double precision, allocatable :: thalf(:)
      integer, allocatable          :: parit(:)
      integer, allocatable          :: ndch(:) 
* Spin and Parity (Energy Level; Variation(1,2,3...))
      double precision, allocatable :: spinpar(:,:)
      integer, allocatable          :: nspv(:) 
* Branching ratio (Same as egamma)
      double precision, allocatable :: bratio(:,:) 
* Internal Conversion coefficient (Same as egamma)
      double precision, allocatable :: braicc(:,:)
* Final level (Same as egamma)
      integer, save, allocatable    :: nlevdn(:,:) 
* Gamma-ray energy (Energy Level; Decay channel)
      double precision, allocatable :: egamma(:,:) 

!$OMP THREADPRIVATE(elevel, thalf, parit, ndch, nlevel)
!$OMP THREADPRIVATE(spinpar, nspv, bratio, braicc, nlevdn, egamma)

* Particle energy 
      integer          ::  nph
      double precision ::  eph(1000),delay(1000)
      integer          ::  kfejec(1000)
!$OMP THREADPRIVATE(nph, eph, kfejec, delay)

      parameter (iprodsize = 3000)  ! number of produced nuclear species 

      integer, dimension(1:120,1:200) :: lev_nth_cap = 0   ! level number of n-capture state

      integer, private, dimension(0:120,0:200) :: levset_order_arr = 0   ! order where this nuclide is set. -1 : data is not available, 0 : data is not yet read, 1 >= assigned at this position.  Arguments : 1st dimension : proton number,  2nd dimension : neutron number
      integer, private :: levset_order = 1 ! number of defined nuclei

      integer, private :: ielev_adr(iprodsize,2) ! elevel address start(1) and end(2)
      integer, private :: ispin_adr(iprodsize,2) ! spin address start and end
      integer, private :: ipari_adr(iprodsize,2) ! parity address start and end
      integer, private :: ilife_adr(iprodsize,2) ! half-life address start and end
      integer, private :: inbra_adr(iprodsize,2) ! # of branches address start and end
      integer, private :: ibran_adr(iprodsize,2) ! branching ratio address start and end
      integer, private :: ibicc_adr(iprodsize,2) ! internal conversion ratio address start and end
      integer, private :: ilvdn_adr(iprodsize,2) ! branch destination level start and end
      
      double precision, allocatable, private :: elev_arr(:)  ! elevel array
      double precision, allocatable, private :: spin_arr(:)  ! spin array
      integer,          allocatable, private :: pari_arr(:)  ! parity array
      double precision, allocatable, private :: life_arr(:)  ! half-life array
      integer,          allocatable, private :: nbra_arr(:)  ! # of branches array
      double precision, allocatable, private :: bran_arr(:)  ! branching ratio array
      double precision, allocatable, private :: bicc_arr(:)  ! branching ratio array
      integer,          allocatable, private :: lvdn_arr(:)  ! branch destination array
!$OMP THREADPRIVATE(lev_nth_cap, levset_order_arr, levset_order)         ! memory wasting but speedy 
!$OMP THREADPRIVATE(ielev_adr,ispin_adr,ipari_adr,ilife_adr,inbra_adr)   ! memory wasting but speedy 
!$OMP THREADPRIVATE(ibran_adr,ibicc_adr,ilvdn_adr)                       ! memory wasting but speedy 
!$OMP THREADPRIVATE(elev_arr,spin_arr,pari_arr,life_arr)                 ! memory wasting but speedy 
!$OMP THREADPRIVATE(nbra_arr,bran_arr,bicc_arr,lvdn_arr)                 ! memory wasting but speedy 

      parameter (maxzab=556)
      integer, private, dimension(maxzab)    :: ismzab ! Isomer definition table
      integer, public                        :: isolev1(120,200)
      integer, public                        :: isolev2(120,200)


c Isomer definition table taken from DCHAIN
c          ismzab .... stored za having isomeric states
c             positive : only 1st meta-stable
c             negative : 1st and 2nd meta-stables
c             (note) below za truncated because of unknown nuclide as
c                    decay data:
c                    28068, 65145, 67156, 72184, 79184, 80183
      data ismzab
     & /  11024,  13024,  13026,  17034,  17038,  19038,          21042,
     &    21044,  21045,  21046,  21050,  23046,          25050,  25052,
     &    25058,  25060,  26052,  26053,  27053,  27054,  27058,  27060,
     &    27062,          29068,  29070,  29076,  30069,  30071,  30073,
     &    30077,  31072,  31074,  32071,  32073,  32075,  32077,  32079,
     &    33082,  34073,  34077,  34079,  34081,  34083,  35070,  35072,
     &    35074,  35076,  35077,  35079,  35080,  35082,  35084,  36079,
     &    36081,  36083,  36085,  37078,  37081,  37082,  37084,  37086,
     &            37090,  37098,  38083,  38085,  38087,  39083,  39084,
     &    39085,  39086,  39087,  39089,  39090,  39091,  39093,  39096,
     &   -39097,  39098,  39100,  39102,  40083,  40085,  40087,  40089,
     &    40090,  41086,  41087,  41088,  41089,  41090,  41091,  41092,
     &    41093,  41094,  41095,  41097,  41098,  41099,  41100,  41102,
     &    41104,  42087,  42089,  42091,  42093,  43088,  43089,  43090,
     &    43091,  43093,  43094,  43095,  43096,  43097,  43099,  43102,
     &    44091,  44093,  44103,  45094,  45095,  45096,  45097,  45098,
     &    45099,  45100,  45101,  45102,  45103,  45104,  45105,  45106,
     &    45108,  45110,  45112,  45116,  46107,  46109,  46111,  46113,
     &    46115,  46117,  47094,  47099,  47100,  47101,  47102,  47103,
     &    47104,  47105,  47106,  47107,  47108,  47109,  47110,  47111,
     &    47113,  47114,  47115,  47116,  47117,  47118,  47119,  47120,
     &    48111,  48113,  48115,  48117,  48119,  48121,  48123,  48125,
     &    49104,  49105,  49106,  49107,  49108, -49109,  49110,  49111,
     &    49112,  49113, -49114,  49115, -49116,  49117, -49118,  49119,
     &   -49120,  49121, -49122,  49123,  49124,  49125,  49126,  49127,
     &    49128,  49129,  49130, -49131,  50113,  50117,  50119,  50121,
     &    50123,  50125,  50127,  50128,  50129,  50130,  50131,  51116,
     &    51118,  51120,  51122,  51124,  51126,  51128,  51129,  51130,
     &    51132,  51134,  52115,  52117,  52119,  52121,  52123,  52125,
     &    52127,          52129,  52131,  52133,  53114,  53118,  53120,
     &    53130,  53132,  53133,  53134,  53136,  54125,  54127,  54129,
     &    54131,  54132,  54133,  54134,  54135,  55116,  55117,  55118,
     &    55119,  55120,  55121, -55122,  55123,  55124,  55130,  55134,
     &    55135,  55136,  55138,  56127,  56129,  56130,  56131,  56133,
     &    56135,  56136,  56137,  57127,  57129,  57132,  57136,
     &    57146,  58131,  58132,  58133,  58135,  58137,  58138,  58139,
     &    59131,  59134,  59138,  59142,  59144,  59148,  60133,  60135,
     &    60137,  60139,          60141,          61134,  61136,  61138,
     &    61139,  61140,  61148, -61152,  61154,  62139,  62141, -62143,
     &                    63136,  63140,  63141,  63142,  63150, -63152,
     &    63154,  64141,  64143,  64145,          64155,  65141,  65142,
     &    65143,  65144,          65146,  65147,  65148,  65149,  65150,
     &    65151,  65152,  65154,  65156,  65158,  66145,  66146,  66147,
     &    66149,  66165,  67148,  67149,  67150,  67151,  67152,  67153,
     &    67154,         -67158,  67159, -67160,  67161,  67162,  67163,
     &    67164,  67166,  67168,  67170,  68149,  68151,  68167,  69146,
     &    69147, -69151,  69152,  69153,  69154,  69155,  69156,  69160,
     &    69162,  69164,  70151,  70169,  70175,  70176,  70177,  71151,
     &   -71155,  71156,  71157,  71160,  71161, -71162,  71165, -71166,
     &    71168,  71169,  71170,  71171,  71172,  71174,  71176,  71177,
     &    71178, -72156,         -72177, -72178, -72179,  72180,  72182,
     &            73178,  73180, -73182,  74158,  74179,  74180,  74183,
     &    74185,  75167,  75168,  75172,  75182,  75184,  75186,
     &    75188,  75190,  76181,  76183,          76189,  76190,  76191,
     &    76192,  77172,  77173,  77174,  77186,  77187, -77190, -77191,
     &   -77192,  77193,  77194,  77195,  77196,  77197,  78183,  78185,
     &            78193,  78195,  78197,  78199,          79185,  79187,
     &    79189,  79191,  79192,  79193, -79194,  79195, -79196,  79197,
     &    79198,  79200,          80185,  80187,  80189,  80191,  80193,
     &    80195,  80197,  80199,  81179,  81183,  81185,  81186,  81187,
     &    81189,  81190,  81192,  81193,  81194,  81195,  81196,  81197,
     &    81198,  81200,  81206,  81207,  82187,  82191,  82195,  82197,
     &    82199,  82201,  82202, -82203,  82204,  82207,  83187,  83188,
     &    83189,  83190,  83191,  83192,  83193, -83194,  83195, -83196,
     &    83197, -83198,  83199, -83200,  83201, -83204,  83210, -83212,
     &    84193,  84195,  84197,  84199,  84201,  84203,  84207,  84211,
     &    84212,  85197,  85198, -85200, -85202,  85204,  85212, -85214,
     &    86197,  86199,  86201,  86203,  87200, -87204,  87206,  87213,
     &    87214,  87218,  88203,  88205,  88207,  88213,  88216,  89208,
     &    89217,  89222,  90216,          91217,  91234,  92234,  92235,
     &    92238,  93236,  93240,  93242,  94237,  95242,  95244,  95246,
     &    97248,  99254,  99256, 100247, 101254, 101258, 102254, 105258,
     &   107262  /
     
      common /qparm/  ielas,icasc,iqstep,lvlopt,igamma ! use igamma
      
      contains

************************************************************************
*                                                                      *
      subroutine allocRIPL(id_alloc)
*                                                                      *
*        allocation of arrays for EBITEM                               *
*        created by T.Ogawa on 2020/02/24                              *
*                                                                      *
*        parameters (transferred as arguments) :                       *
*                                                                      *
*        Input :                                                       *
*        id_alloc : which data to be allocated                         *
*                 0-> initialization                                   *
*                 1-> elev_arr(:)                                      *
*                 2-> spin_arr(:)                                      *
*                 3-> pari_arr(:)                                      *
*                 4-> life_arr(:)                                      *
*                 5-> nbra_arr(:)                                      *
*                 6-> bran_arr(:)                                      *
*                 7-> bicc_arr(:)                                      *
*                 8-> lvdn_arr(:)                                      *
*                                                                      *
*        Output : none                                                 *
*                                                                      *
************************************************************************

      implicit double precision(a-h, o-z)

      ! temporary arrays for reallocation
      double precision, allocatable :: elev_arr_tmp(:)
      double precision, allocatable :: spin_arr_tmp(:)
      integer         , allocatable :: pari_arr_tmp(:)
      double precision, allocatable :: life_arr_tmp(:)
      integer         , allocatable :: nbra_arr_tmp(:)
      double precision, allocatable :: bran_arr_tmp(:)
      double precision, allocatable :: bicc_arr_tmp(:)
      integer         , allocatable :: lvdn_arr_tmp(:)
      

      Selectcase(id_alloc) ! Extend the array if not enough
      Case(0)
         allocate(elev_arr(10000))
         allocate(spin_arr(10000))
         allocate(pari_arr(10000))
         allocate(life_arr(10000))
         allocate(nbra_arr(100000))
         allocate(bran_arr(100000))
         allocate(bicc_arr(100000))
         allocate(lvdn_arr(100000))
         
         ielev_adr(1,1) = 1
         ispin_adr(1,1) = 1
         ipari_adr(1,1) = 1
         ilife_adr(1,1) = 1
         inbra_adr(1,1) = 1
         ibran_adr(1,1) = 1
         ibicc_adr(1,1) = 1
         ilvdn_adr(1,1) = 1
         return
      Case(1)
         length = size(elev_arr)
         allocate(elev_arr_tmp(length)) 
         elev_arr_tmp = elev_arr
comp !$OMP CRITICAL (RIPL_elev_crit) ! prohibit access to elev_arr
         deallocate(elev_arr)
         allocate(elev_arr(length*2))
         elev_arr(1:length) = elev_arr_tmp(1:length)
         deallocate(elev_arr_tmp)
comp !$OMP END CRITICAL (RIPL_elev_crit)
      Case(2)
         length = size(spin_arr)
         allocate(spin_arr_tmp(length)) 
         spin_arr_tmp = spin_arr
comp !$OMP CRITICAL (RIPL_spin_crit) ! prohibit access 
         deallocate(spin_arr)
         allocate(spin_arr(length*2))
         spin_arr(1:length) = spin_arr_tmp(1:length)
         deallocate(spin_arr_tmp)
comp !$OMP END CRITICAL (RIPL_spin_crit)
      Case(3)
         length = size(pari_arr)
         allocate(pari_arr_tmp(length)) 
         pari_arr_tmp = pari_arr
comp !$OMP CRITICAL (RIPL_pari_crit) ! prohibit access 
         deallocate(pari_arr)
         allocate(pari_arr(length*2))
         pari_arr(1:length) = pari_arr_tmp(1:length)
         deallocate(pari_arr_tmp)
comp !$OMP END CRITICAL (RIPL_pari_crit)
      Case(4)
         length = size(life_arr)
         allocate(life_arr_tmp(length)) 
         life_arr_tmp = life_arr
comp !$OMP CRITICAL (RIPL_life_crit) ! prohibit access
         deallocate(life_arr)
         allocate(life_arr(length*2))
         life_arr(1:length) = life_arr_tmp(1:length)
         deallocate(life_arr_tmp)
comp !$OMP END CRITICAL (RIPL_life_crit)
      Case(5)
         length = size(nbra_arr)
         allocate(nbra_arr_tmp(length)) 
         nbra_arr_tmp = nbra_arr
comp !$OMP CRITICAL (RIPL_nbra_crit) ! prohibit access 
         deallocate(nbra_arr)
         allocate(nbra_arr(length*2))
         nbra_arr(1:length) = nbra_arr_tmp(1:length)
         deallocate(nbra_arr_tmp)
comp !$OMP END CRITICAL (RIPL_nbra_crit)
      Case(6)
         length = size(bran_arr)
         allocate(bran_arr_tmp(length)) 
         bran_arr_tmp = bran_arr
comp !$OMP CRITICAL (RIPL_bran_crit) ! prohibit access 
         deallocate(bran_arr)
         allocate(bran_arr(length*2))
         bran_arr(1:length) = bran_arr_tmp(1:length)
         deallocate(bran_arr_tmp)
comp !$OMP END CRITICAL (RIPL_bran_crit)
      Case(7)
         length = size(bicc_arr)
         allocate(bicc_arr_tmp(length)) 
         bicc_arr_tmp = bicc_arr
comp !$OMP CRITICAL (RIPL_bicc_crit) ! prohibit access 
         deallocate(bicc_arr)
         allocate(bicc_arr(length*2))
         bicc_arr(1:length) = bicc_arr_tmp(1:length)
         deallocate(bicc_arr_tmp)
comp !$OMP END CRITICAL (RIPL_bicc_crit)
      Case(8)
         length = size(lvdn_arr)
         allocate(lvdn_arr_tmp(length)) 
         lvdn_arr_tmp = lvdn_arr
comp !$OMP CRITICAL (RIPL_lvdn_crit) ! prohibit access 
         deallocate(lvdn_arr)
         allocate(lvdn_arr(length*2))
         lvdn_arr(1:length) = lvdn_arr_tmp(1:length)
         deallocate(lvdn_arr_tmp)
comp !$OMP END CRITICAL (RIPL_lvdn_crit)
      ENDSELECT
      
      return
      end subroutine
      
      
************************************************************************
*                                                                      *
      subroutine levset(iares, izres, lflg)
*                                                                      *
*       nuclear structure setup based on RIPL                          *
*       created by T.Ogawa on 2020/02/24                               *
*                                                                      *
*       ENSDF Based Isomeric Transition/ isomEr production Model       *
*               -->  EBITEM                                            *
*                                                                      *
*        parameters (transferred as arguments) :                       *
*                                                                      *
*        Input :                                                       *
*        iares  : mass number of residual nucleus                      *
*        izres  : charge number of residual nucleus                    *
*                                                                      *
*        Output :                                                      *
*        lflg  : Level structure setup                                 *
*                      true : setup success   false : setup failed     *
*                                                                      *
************************************************************************

      implicit double precision(a-h, o-z)
      logical lflg

      if( izres .lt. 0 .or. izres .gt. 120 .or. 
     &    iares .lt. 0 .or. iares .gt. 294                  ) then
           lflg = .false.
           return ! structure unavailable. 
      elseif( levset_order_arr(izres, iares - izres) .eq. -1) then
           lflg = .false.
           return ! structure unavailable. 
      endif
      
comp !$OMP CRITICAL (RIPL_crit)
      if( levset_order_arr(izres, iares - izres) .eq. 0) then ! level not defined. read RIPL data.
           call levread(izres, iares, levset_order, lflg)
           if(lflg ) then
               levset_order_arr(izres, iares - izres) = levset_order
               levset_order = levset_order + 1 
           else
               levset_order_arr(izres, iares - izres) = -1 ! Data missing, error, or only ground state is available
               return
           endif
      endif
comp !$OMP END CRITICAL (RIPL_crit)

      lflg = .true.
      call levcall(izres, iares)  ! call structure data based on stored array

      return
      end subroutine
      
      
************************************************************************
*                                                                      *
      subroutine levread(izres, iares, levset_order, lflg)
*                                                                      *
*       read RIPL nuclear structure.                                   *
*       modified version of read-levels.f90 provided by IAEA           *
*       created by T.Ogawa on 2020/02/24                               *
*                                                                      *
*        parameters (transferred as arguments) :                       *
*                                                                      *
*        Input :                                                       *
*        iares  : mass number of residual nucleus                      *
*        izres  : charge number of residual nucleus                    *
*                                                                      *
*        Output :                                                      *
*        lflg  : Level structure setup                                 *
*                      true : setup success   false : setup failed     *
*                                                                      *
************************************************************************
!$    use omp_lib
      use NGSDATAMOD, only : bindeg

      implicit double precision(a-h, o-z)
      
      include 'param-physcnst.inc'
      
      logical lflg     ! this is the nucleus I want
      logical lisosrch ! looking for isomer?
      common /ipomp0/ipomp,npomp
!$OMP THREADPRIVATE(/ipomp0/)
      
      common /paran/ icfn(100), ilfn(100), chfn(100)
      character chfn*200
      character hdpth*300
! PARAMETER definitions
!
      INTEGER, PARAMETER :: LEVMAX = 1000, GAMMAX = 10000, BRMAX = 100
!
! COMMON variables
!
      INTEGER :: A, I, J, JJ, K, M, N, NC, NGAm, NLEv, NMAx, W, Z
      LOGICAL :: BOTh, EVEn, HASq
      double precision :: DIFf, EG, IT_percent, L1, LL1, PF, PI, RISum, 
     & RITot, SF, SI, SN, SP, SUM, TOTsum
      double precision :: DMIx, EGD, ICC
      CHARACTER(300) :: FNAme
      CHARACTER(1) :: IMPos
      CHARACTER(10) :: STRmul, TEXt
      CHARACTER(5) :: SYMb

      double precision, allocatable :: CC_gam(:), E_Gam(:), RI1_gam(:), 
     & RI_gam(:), MIXr_gam(:), PROb_gam(:)
      INTEGER, allocatable :: FINal_gam(:), INItial_gam(:)
      CHARACTER(10), allocatable :: MULt_gam(:)
      CHARACTER(7), allocatable :: DMOdes_lev(:,:)
      double precision, allocatable :: DPErcent_lev(:,:)
      double precision, allocatable :: E_Lev(:), JPMy_lev(:), T_Lev(:), 
     & JSMy_lev(:)
      CHARACTER(1), allocatable :: JPEstimate_lev(:)
      CHARACTER(18), allocatable :: JPText_lev(:)
      INTEGER, allocatable :: NOBr_lev(:), NOG_lev(:)
      CHARACTER(4), allocatable :: UNCertaint_lev(:)
      CHARACTER(2), allocatable :: PFIx_lev(:,:)
      
! Local variables
      INTEGER :: aa, zz
!
      lflg = .false.
!Check if data request is in range
      IF(izres .LT. 0 .OR. izres .GE. 118)THEN
         WRITE(*, *)'Data with this charge number Z does not exists Z=',
     & izres
         RETURN
      ENDIF
      IF(iares .LT. 1 .OR. iares .GT. 294)THEN
         WRITE(*, *)'Data with this mass number A does not exists A=', 
     & iares
         RETURN
      ENDIF
!open file with the give Z
      hdpth=chfn(1)(1:ilfn(1))//'/data/RIPL/' ! default value of datapath in xsdir
      IF(izres .LT. 10)WRITE(FNAme, '(a3,i1,a4)')'z00', izres, '.dat'
      IF(izres .LT. 100 .AND. izres .GE.10) WRITE(FNAme, '(a2,i2,a4)')'z
     &0', izres,'.dat'
      IF(izres .GT. 99) WRITE(FNAme, '(a1,i3,a4)')'z', izres, '.dat'
      
      FNAme = hdpth(1:leng(hdpth))//FNAme
      
      allocate(CC_gam(GAMMAX), E_Gam(GAMMAX), RI1_gam(GAMMAX), 
     & RI_gam(GAMMAX), MIXr_gam(GAMMAX), PROb_gam(GAMMAX),
     & FINal_gam(GAMMAX), INItial_gam(GAMMAX), MULt_gam(GAMMAX))
      allocate(DMOdes_lev(LEVMAX, BRMAX), DPErcent_lev(LEVMAX, BRMAX), 
     & PFIx_lev(LEVMAX, BRMAX))
      allocate(E_Lev(LEVMAX), JPMy_lev(LEVMAX), JSMy_lev(LEVMAX), 
     & T_Lev(LEVMAX), JPEstimate_lev(LEVMAX), JPText_lev(LEVMAX),
     & NOBr_lev(LEVMAX), NOG_lev(LEVMAX), UNCertaint_lev(LEVMAX))
      
      OPEN(301+ipomp, FILE = FNAme, STATUS = 'old',ERR=200)
      
      do while (.not. lflg) ! scroll down the file 
      READ(301+ipomp, '(a5,6i5,2f12.6)', END = 200)SYMb, aa, zz, NLEv, 
     & NGAm, NMAx, NC, SN, SP

      J = 0
      JJ = 0
      IF(iares.EQ.aa .AND. izres.EQ.zz) lflg = .true. ! data exists
          
          if(lflg .and. NLEv .eq. 1) then
             lflg = .false.
             CLOSE(301+ipomp)
             return
          endif
          
          lisosrch = .true.
          
          DO JJ = 1, NLEv
          READ(301+ipomp, 
     &    '(i3,1x,f10.6,1x,f5.1,i3,1x,(1pe10.2),i3,1x,a1,1x,a4,1x,a18,i3
     &    ,10(1x,a2,1x,0pf10.4,1x,a7))')I, E_Lev(I), JSMy_lev(I), W, 
     &    T_Lev(I),  NOG_lev(I), JPEstimate_lev(I), UNCertaint_lev(I), 
     &    JPText_lev(I),NOBr_lev(I), (PFIx_lev(I, M), DPErcent_lev(I, M)
     &    , DMOdes_lev(I, M), M = 1, NOBr_lev(I))

              if( lflg ) then ! correction of data if nucleus of interest
              
              if( abs( bindeg(izres, iares - izres)
     &               - bindeg(izres, iares - izres - 1) - E_Lev(I) )
     &               .lt. 1.d-3 .and. NOG_lev(I) .gt. 1)
     &        lev_nth_cap(izres, iares - izres) = I-1 ! neutron absorption state
              
              if( W .eq. 0 ) then ! parity is unknown
                    if( JPText_lev(I)(17:17) .eq. '-' .or. 
     &                  JPText_lev(I)(18:18) .eq. '-' ) then
                        W = -1
                    elseif( JPText_lev(I)(17:17) .eq. '+' .or.  
     &                      JPText_lev(I)(18:18) .eq. '+' ) then
                        W = 1
                    else
                        W = 1 ! assumed to be +. Keep coherence between threads in OpenMP mode
	              endif
              endif

              if( JSMy_lev(I) .eq. -1.d0 ) then ! spin is unknown
                 if( JPText_lev(I) .eq. '                 +' .or. 
     &               JPText_lev(I) .eq. '                 -' .or. 
     &               JPText_lev(I) .eq. '               (+)' .or. 
     &               JPText_lev(I) .eq. '               (-)' .or. 
     &               JPText_lev(I) .eq. '           NATURAL' .or. 
     &               JPText_lev(I) .eq. '         UNNATURAL' .or. 
     &               JPText_lev(I) .eq. '                  ' ) then
                    JSMy_lev(I) = -1.d0
                 elseif( JPText_lev(I)(1:1) .eq. 'J' .or. 
     &               JPText_lev(I)(2:2) .eq. 'J' .or. 
     &               JPText_lev(I)(3:3) .eq. 'J' .or. 
     &               JPText_lev(I)(4:4) .eq. 'J' .or. 
     &               JPText_lev(I)(5:5) .eq. 'J' .or. 
     &               JPText_lev(I)(6:6) .eq. 'J' .or. 
     &               JPText_lev(I)(7:7) .eq. 'J' .or. 
     &               JPText_lev(I)(8:8) .eq. 'J' .or. 
     &               JPText_lev(I)(9:9) .eq. 'J' .or. 
     &               JPText_lev(I)(10:10) .eq. 'J' .or. 
     &               JPText_lev(I)(11:11) .eq. 'J' .or. 
     &               JPText_lev(I)(12:12) .eq. 'J' .or. 
     &               JPText_lev(I)(13:13) .eq. 'J' .or. 
     &               JPText_lev(I)(14:14) .eq. 'J' .or. 
     &               JPText_lev(I)(15:15) .eq. 'J' .or. 
     &               JPText_lev(I)(16:16) .eq. 'J' .or. 
     &               JPText_lev(I)(17:17) .eq. 'J' .or. 
     &               JPText_lev(I)(18:18) .eq. 'J' ) then
                    JSMy_lev(I) = -1.d0
                 elseif( JPText_lev(I) .ne. '                  ' ) then
                    if( mod(iares, 2) .eq. 0 ) then
                         kk = 1
                         do while ( JPText_lev(I)(kk:kk) .ne. '0'
     &                   .and. JPText_lev(I)(kk:kk) .ne. '1' 
     &                   .and. JPText_lev(I)(kk:kk) .ne. '2' 
     &                   .and. JPText_lev(I)(kk:kk) .ne. '3' 
     &                   .and. JPText_lev(I)(kk:kk) .ne. '4' 
     &                   .and. JPText_lev(I)(kk:kk) .ne. '5' 
     &                   .and. JPText_lev(I)(kk:kk) .ne. '6' 
     &                   .and. JPText_lev(I)(kk:kk) .ne. '7' 
     &                   .and. JPText_lev(I)(kk:kk) .ne. '8' 
     &                   .and. JPText_lev(I)(kk:kk) .ne. '9')
                            kk = kk + 1
                         enddo
                         if(kk .eq. 18 .or. 
     &                    JPText_lev(I)(kk+1:kk+1) .eq. '(' .or.
     &                    JPText_lev(I)(kk+1:kk+1) .eq. '[' .or.
     &                    JPText_lev(I)(kk+1:kk+1) .eq. ',' .or.
     &                    JPText_lev(I)(kk+1:kk+1) .eq. '+' .or.
     &                    JPText_lev(I)(kk+1:kk+1) .eq. '-') then ! pickup the first spin specification 
                            read( JPText_lev(I)(kk:kk), *) iii
                            JSMy_lev(I) = iii
                         elseif(
     &                    JPText_lev(I)(kk+1:kk+1) .eq. '1' .or.
     &                    JPText_lev(I)(kk+1:kk+1) .eq. '2' .or.
     &                    JPText_lev(I)(kk+1:kk+1) .eq. '3' .or.
     &                    JPText_lev(I)(kk+1:kk+1) .eq. '4' .or.
     &                    JPText_lev(I)(kk+1:kk+1) .eq. '5' .or.
     &                    JPText_lev(I)(kk+1:kk+1) .eq. '6' .or.
     &                    JPText_lev(I)(kk+1:kk+1) .eq. '7' .or.
     &                    JPText_lev(I)(kk+1:kk+1) .eq. '8' .or.
     &                    JPText_lev(I)(kk+1:kk+1) .eq. '9' .or.
     &                    JPText_lev(I)(kk+1:kk+1) .eq. '0') then ! pickup the first spin specification 
                            read( JPText_lev(I)(kk:kk+1), *) iii
                            JSMy_lev(I) = iii
                         elseif(JPText_lev(I)(kk-3:kk-3) .eq. '('
     &                     .or. JPText_lev(I)(kk-3:kk-3) .eq. '['
     &                     .or. JPText_lev(I)(kk-3:kk-3) .eq. ',')
     &                       then ! in case spin is 11/2, 13/2 etc
                            read( JPText_lev(I)(kk:kk+1), *) iii
                            JSMy_lev(I) = iii
                         endif
                    else ! odd A nucleus
                           kk = 1
                         do while( JPText_lev(I)(kk:kk) .ne. '/' .and. 
     &                    kk .lt. 18)
                            kk = kk + 1
                         enddo
                         if(kk .eq. 2) then                        ! string is full 
                            read( JPText_lev(I)(kk-1:kk-1), *) iii
                            JSMy_lev(I) = dble(iii/2.d0)
                         elseif(JPText_lev(I)(kk-2:kk-2) .eq. '('
     &                     .or. JPText_lev(I)(kk-2:kk-2) .eq. '['
     &                     .or. JPText_lev(I)(kk-2:kk-2) .eq. ' ') then! pickup the first spin specification 
                            read( JPText_lev(I)(kk-1:kk-1), *) iii
                            JSMy_lev(I) = dble(iii/2.d0)
                         elseif(JPText_lev(I)(kk-3:kk-3) .eq. '(' 
     &                     .or. JPText_lev(I)(kk-3:kk-3) .eq. '['
     &                     .or. JPText_lev(I)(kk-3:kk-3) .eq. ' ')
     &                     then ! in case spin is 11/2, 13/2 etc
                            read( JPText_lev(I)(kk-2:kk-1), *) iii
                            JSMy_lev(I) = dble(iii/2.d0)
                         endif
                    endif
                 endif ! if completely unknown, do not use this level afterwards
              endif

              if( NOBr_lev(I) .gt. 0 .and. lisosrch .and. 
     &         abs(igamma) .ge. 2 .and. I .ne. 1 ) then ! Isomer level definition

                  nisomer = 0
                  do ii = 1, maxzab
                        if(ismzab(ii) .eq. izres*1000+iares) then
                              nisomer = 1
                        elseif(ismzab(ii) .eq. -(izres*1000+iares)) then
                              nisomer = 2
                        endif
                  enddo
     
                  if(nisomer .eq. 1) then                                ! 1 isomer
                      isolev1(izres, iares - izres) = I - 1
                      lisosrch = .false.
                  elseif(nisomer .eq. 2 .and. isolev1(izres, iares 
     &             - izres) .ne. 0) then                                 ! 2nd isomer of 2 isomers
                      isolev2(izres, iares - izres) = I - 1
                      lisosrch = .false.
                  elseif(nisomer .eq. 2 .and. isolev1(izres, iares 
     &             - izres) .eq. 0) then                                 ! 1st isomer of 2 isomers
                      isolev1(izres, iares - izres) = I - 1
                  else                                                   ! nucleus without isomer 
                      lisosrch = .false.
                  endif
              endif
              endif
              
              JPMy_lev(I) = W
              DO K = 1, NOG_lev(I)
                  J = J + 1
                  READ(301+ipomp, '(39x,i4,1x,f10.3,3(1x,e10.3))')
     &           FINal_gam(J),E_Gam(J), RI_gam(J), RI1_gam(J), CC_gam(J)
              ENDDO
          ENDDO
      enddo
      
      ! first allocation
      if(.not. allocated(elev_arr)) then 
          call allocRIPL(0)
      endif
      
      ! data address specification
      
      if(levset_order .gt. iprodsize) write(*,*) 'Error : in gammod.f 
     &increase iprodsize'
      
      ielev_adr(levset_order,2)   = ielev_adr(levset_order,1)+NLEv - 1
      ielev_adr(levset_order+1,1) = ielev_adr(levset_order,2) + 1
      
      ispin_adr(levset_order,2)   = ispin_adr(levset_order,1)+NLEv - 1 ! spin バラエティ考慮しない。
      ispin_adr(levset_order+1,1) = ispin_adr(levset_order,2) + 1

      ipari_adr(levset_order,2)   = ipari_adr(levset_order,1)+NLEv - 1
      ipari_adr(levset_order+1,1) = ipari_adr(levset_order,2) + 1

      ilife_adr(levset_order,2)   = ilife_adr(levset_order,1)+NLEv - 1
      ilife_adr(levset_order+1,1) = ilife_adr(levset_order,2) + 1

      inbra_adr(levset_order,2)   = inbra_adr(levset_order,1)+NLEv - 1
      inbra_adr(levset_order+1,1) = inbra_adr(levset_order,2) + 1
      
      ibran_adr(levset_order,2)   = ibran_adr(levset_order,1)+ J - 1
      ibran_adr(levset_order+1,1) = ibran_adr(levset_order,2) + 1

      ibicc_adr(levset_order,2)   = ibicc_adr(levset_order,1)+ J - 1
      ibicc_adr(levset_order+1,1) = ibicc_adr(levset_order,2) + 1

      ilvdn_adr(levset_order,2)   = ilvdn_adr(levset_order,1)+ J - 1
      ilvdn_adr(levset_order+1,1) = ilvdn_adr(levset_order,2) + 1
      
      ! Data input
      
comp !$OMP CRITICAL (RIPL_elev_crit) ! prohibit access to elev_arr
      if(ielev_adr(levset_order,2) .gt. size(elev_arr) ) 
     & call allocRIPL(1)
      elev_arr(ielev_adr(levset_order,1):ielev_adr(levset_order,2))
     & = E_Lev(1:NLEv)
comp !$OMP END CRITICAL (RIPL_elev_crit)
comp !$OMP CRITICAL (RIPL_spin_crit) ! prohibit access 
      if(ispin_adr(levset_order,2) .gt. size(spin_arr) ) 
     & call allocRIPL(2)
      spin_arr(ispin_adr(levset_order,1):ispin_adr(levset_order,2))
     & = JSMy_lev(1:NLEv)
comp !$OMP END CRITICAL (RIPL_spin_crit)
comp !$OMP CRITICAL (RIPL_pari_crit) ! prohibit access 
      if(ipari_adr(levset_order,2) .gt. size(pari_arr) ) 
     & call allocRIPL(3)
      pari_arr(ipari_adr(levset_order,1):ipari_adr(levset_order,2))
     & = JPMy_lev(1:NLEv)
comp !$OMP END CRITICAL (RIPL_pari_crit)
comp !$OMP CRITICAL (RIPL_life_crit) ! prohibit access 
      if(ilife_adr(levset_order,2) .gt. size(life_arr) ) 
     & call allocRIPL(4)
      life_arr(ilife_adr(levset_order,1):ilife_adr(levset_order,2))
     & = T_Lev(1:NLEv)
comp !$OMP END CRITICAL (RIPL_life_crit)
comp !$OMP CRITICAL (RIPL_nbra_crit) ! prohibit access 
      if(inbra_adr(levset_order,2) .gt. size(nbra_arr) ) 
     & call allocRIPL(5)
      nbra_arr(inbra_adr(levset_order,1):inbra_adr(levset_order,2))
     & = NOG_lev(1:NLEv)
comp !$OMP END CRITICAL (RIPL_nbra_crit)
comp !$OMP CRITICAL (RIPL_bran_crit) ! prohibit access 
      if(ibran_adr(levset_order,2) .gt. size(bran_arr) ) 
     & call allocRIPL(6)
      bran_arr(ibran_adr(levset_order,1):ibran_adr(levset_order,2)) ! total branching ratio
     & = RI1_gam(1:J)
comp !$OMP END CRITICAL (RIPL_bran_crit)
comp !$OMP CRITICAL (RIPL_bicc_crit) ! prohibit access 
      if(ibicc_adr(levset_order,2) .gt. size(bicc_arr) ) 
     & call allocRIPL(7)
      do ii = 1, J
          if(RI1_gam(ii) .eq. 0.d0) then
              bicc_arr(ibicc_adr(levset_order,1) + ii - 1) = 0.d0 ! workaroud for unknown branching ratio
          else
              bicc_arr(ibicc_adr(levset_order,1) + ii - 1) ! ratio of ICC over total
     &        = RI_gam(ii) * CC_gam(ii) / RI1_gam(ii)
          endif
      enddo
comp !$OMP END CRITICAL (RIPL_bicc_crit)
comp !$OMP CRITICAL (RIPL_lvdn_crit) ! prohibit access 
      if(ilvdn_adr(levset_order,2) .gt. size(lvdn_arr) ) 
     & call allocRIPL(8)
      lvdn_arr(ilvdn_adr(levset_order,1):ilvdn_adr(levset_order,2))
     & = FINal_gam(1:J)
comp !$OMP END CRITICAL (RIPL_lvdn_crit)

200   CLOSE(301+ipomp) ! Data is not existing
      deallocate(CC_gam, E_Gam, RI1_gam, RI_gam, MIXr_gam, PROb_gam,
     & FINal_gam, INItial_gam, MULt_gam, DMOdes_lev, DPErcent_lev, 
     & E_Lev, JPMy_lev, JSMy_lev, T_Lev, JPEstimate_lev, JPText_lev,
     & NOBr_lev, NOG_lev, UNCertaint_lev, PFIx_lev)
      return     ! return with lflg .false.
      
      end subroutine
      

************************************************************************
*                                                                      *
      subroutine levcall(izres, iares)
*                                                                      *
*       call stored RIPL nuclear structure.                            *
*       created by T.Ogawa on 2020/02/24                               *
*                                                                      *
*        parameters (transferred as arguments) :                       *
*                                                                      *
*        Input :                                                       *
*        iares  : mass number of residual nucleus                      *
*        izres  : charge number of residual nucleus                    *
*                                                                      *
************************************************************************

      implicit double precision(a-h, o-z)
      
      id = levset_order_arr(izres, iares-izres) ! nuclear storage id
      
      nlevel = ielev_adr(id,2) - ielev_adr(id,1) + 1
      allocate(elevel(0:nlevel-1))
comp !$OMP CRITICAL (RIPL_elev_crit) ! prohibit access to elev_arr
      elevel(0:nlevel-1) = elev_arr(ielev_adr(id,1):ielev_adr(id,2)) ! shift by 1 because EBITEM ground state is 0th, not 1st
comp !$OMP END CRITICAL (RIPL_elev_crit)
      
      allocate( thalf(0:nlevel-1))
comp !$OMP CRITICAL (RIPL_life_crit) ! prohibit access
      thalf(0:nlevel-1) = life_arr(ilife_adr(id,1):ilife_adr(id,2))
comp !$OMP END CRITICAL (RIPL_life_crit)

      allocate( nspv(0:nlevel-1)) 
      nspv(0:nlevel-1) = 1 ! temporal workaround
      
      length = ispin_adr(id,2) - ispin_adr(id,1) + 1 
      allocate(spinpar(0:length-1, 1)) 
!      spinpar = spin_arr(ispin_adr(id,1):ispin_adr(id,2)) ! スピンの広がりを考慮しない

comp !$OMP CRITICAL (RIPL_spin_crit) ! prohibit access 
      k = 0
      do i = 0, ielev_adr(id,2) - ielev_adr(id,1) ! assign spin variety for each level
         do j = 1, nspv(i) ! **** spin variety of the level
           spinpar(i,j) = spin_arr(ispin_adr(id,1)+k)
           k = k + 1
         enddo
      enddo
comp !$OMP END CRITICAL (RIPL_spin_crit)

      
      allocate( parit(0:nlevel-1))
comp !$OMP CRITICAL (RIPL_pari_crit) ! prohibit access 
      parit(0:nlevel-1) = pari_arr(ipari_adr(id,1):ipari_adr(id,2))
comp !$OMP END CRITICAL (RIPL_pari_crit)

      nlevel = inbra_adr(id,2) - inbra_adr(id,1) + 1 
      allocate(  ndch(0:nlevel-1))
comp !$OMP CRITICAL (RIPL_nbra_crit) ! prohibit access 
      ndch(0:nlevel-1) = nbra_arr(inbra_adr(id,1):inbra_adr(id,2))
comp !$OMP END CRITICAL (RIPL_nbra_crit)

      maxbr = maxval(nbra_arr(inbra_adr(id,1):inbra_adr(id,2)))
      
      allocate(bratio(0:nlevel-1, maxbr)) 
      bratio = 0.d0
      allocate(braicc(0:nlevel-1, maxbr)) 
      braicc = 0.d0
      allocate(nlevdn(0:nlevel-1, maxbr)) 
      nlevdn = 0

comp !$OMP CRITICAL (RIPL_bran_crit) ! prohibit access 
comp !$OMP CRITICAL (RIPL_bicc_crit) ! prohibit access 
comp !$OMP CRITICAL (RIPL_lvdn_crit) ! prohibit access 
      k = 0
      do i = 1, nlevel-1 ! assign branch scheme
         do j = 1, ndch(i)
           bratio(i,j) = bran_arr(ibran_adr(id,1)+k)
           braicc(i,j) = bicc_arr(ibicc_adr(id,1)+k)
           nlevdn(i,j) = lvdn_arr(ilvdn_adr(id,1)+k) - 1 ! definition of level is shifted (ground = 0 in EBITEM, 1 in RIPL)
           k = k + 1
         enddo
      enddo
comp !$OMP END CRITICAL (RIPL_lvdn_crit)
comp !$OMP END CRITICAL (RIPL_bicc_crit)
comp !$OMP END CRITICAL (RIPL_bran_crit)
      

      return
      end subroutine
      
************************************************************************
*                                                                      *
      subroutine levUNset
*                                                                      *
*       deallocate arrays for nuclear structure                        *
*       created by T.Ogawa on 2020/04/03                               *
*                                                                      *
************************************************************************

      deallocate( elevel, thalf, spinpar, parit, ndch, bratio, braicc, 
     & nlevdn, nspv)  

      return
      end subroutine

************************************************************************
*                                                                      *
      subroutine neufdc(iares, izres, eexe, nsta) 
*                                                                      *
*      first decay for neutron capture                                 *
*      Data : https://www-nds.iaea.org/pgaa/PGAAdatabase/LANL/lanl.htm *
************************************************************************

      implicit double precision(a-h, o-z)

      br = unirn(dummy)

      SelectCase(izres*1000+iares)
         Case(3007) 
           IF(br .le. 0.618401207) then
            nsta = 0
           ELSE
            nsta = 1
           ENDIF

         Case( 3008) 
           IF(br .le. 0.894004283) then
            nsta = 0
           ELSE
            nsta = 1
           ENDIF

         Case( 4010) 
           IF(br .le. 0.655) then
            nsta = 0
           ELSEIF(br .le. 0.761) then
            nsta = 1
           ELSEIF(br .le. 0.996) then
            nsta = 2
           ELSEIF(br .le. 0.9984) then
            nsta = 4
           ELSE
            nsta = 5
           ENDIF

         Case( 5011 ) 
           IF(br .le. 0.047) then
            nsta = 0
           ELSEIF(br .le. 0.601) then
            nsta = 2
           ELSEIF(br .le. 0.857) then
            nsta = 4
           ELSE
            nsta = 9
           ENDIF

         Case( 6013 ) 
           IF(br .le. 0.6747) then
            nsta = 0
           ELSEIF(br .le. 0.6763) then
            nsta = 1
           ELSE
            nsta = 2
           ENDIF

         Case( 6014 ) 
           IF(br .le. 0.84) then
            nsta = 0
           ELSEIF(br .le. 0.865) then
            nsta = 1
           ELSEIF(br .le. 0.95) then
            nsta = 2
           ELSE
            nsta = 4
           ENDIF

         Case( 7015 ) 
           IF(br .le. 0.142076503) then
            nsta = 0
           ELSEIF(br .le. 0.248385494) then
            nsta = 1
           ELSEIF(br .le. 0.443119722) then
            nsta = 2
           ELSEIF(br .le. 0.609041232) then
            nsta = 3
           ELSEIF(br .le. 0.753104819) then
            nsta = 4
           ELSEIF(br .le. 0.84199) then
            nsta = 5
           ELSEIF(br .le. 0.89738) then
            nsta = 7
           ELSEIF(br .le. 0.89814) then
            nsta = 8
           ELSEIF(br .le. 0.90061) then
            nsta = 9
           ELSEIF(br .le. 0.91693) then
            nsta = 10
           ELSEIF(br .le. 0.99593) then
            nsta = 11
           ELSEIF(br .le. 0.99665) then
            nsta = 12
           ELSEIF(br .le. 0.99753) then
            nsta = 13
           ELSEIF(br .le. 0.99913) then
            nsta = 15
           ELSEIF(br .le. 0.99974) then
            nsta = 16
           ELSEIF(br .le. 0.99981) then
            nsta = 17
           ELSE
            nsta = 20
           ENDIF


         Case( 7016 ) 
           IF(br .le. 0.179538616) then
            nsta = 1
           ELSE
            nsta = 2
           ENDIF

         Case( 8017 ) 
           IF(br .le. 0.180) then
            nsta = 1
           ELSE
            nsta = 2
           ENDIF

         Case( 8018 ) 
           IF(br .le. 0.378d0) then
            nsta = 1
           ELSEIF(br .le. 0.489d0) then
            nsta = 1
           ELSE ! not yet known
            nsta = 18 ! 20
           ENDIF

         Case( 9020 ) 
           IF(br .le. 0.0985d0) then
            nsta = 0      
           ELSEIF(br .le. 0.113d0) then
            nsta = 3      
           ELSEIF(br .le. 0.1559d0) then
            nsta = 4      
           ELSEIF(br .le. 0.1807d0) then
            nsta = 5      
           ELSEIF(br .le. 0.2005d0) then
            nsta = 7      
           ELSEIF(br .le. 0.2011d0) then
            nsta = 8      
           ELSEIF(br .le. 0.2558d0) then
            nsta = 9      
           ELSEIF(br .le. 0.281d0) then
            nsta = 15     
           ELSEIF(br .le. 0.3008d0) then
            nsta = 16     
           ELSEIF(br .le. 0.3432d0) then
            nsta = 17     
           ELSEIF(br .le. 0.3531d0) then
            nsta = 20     
           ELSEIF(br .le. 0.3632d0) then
            nsta = 22     
           ELSEIF(br .le. 0.3706d0) then
            nsta = 23     
           ELSEIF(br .le. 0.3828d0) then
            nsta = 26     
           ELSEIF(br .le. 0.3883d0) then
            nsta = 28     
           ELSEIF(br .le. 0.3932d0) then
            nsta = 32     
           ELSEIF(br .le. 0.3959d0) then
            nsta = 37     
           ELSEIF(br .le. 0.3965d0) then
            nsta = 42     
           ELSEIF(br .le. 0.3989d0) then
            nsta = 44     
           ELSEIF(br .le. 0.4079d0) then
            nsta = 45     
           ELSEIF(br .le. 0.4088d0) then
            nsta = 51     
           ELSEIF(br .le. 0.4274d0) then
            nsta = 52     
           ELSEIF(br .le. 0.4338d0) then
            nsta = 55     
           ELSEIF(br .le. 0.4342d0) then
            nsta = 62     
           ELSEIF(br .le. 0.5903d0) then
            nsta = 63     
           ELSEIF(br .le. 0.601d0) then
            nsta = 64     
           ELSEIF(br .le. 0.9783d0) then
            nsta = 67     
           ELSEIF(br .le. 0.9995d0) then
            nsta = 68     
           ELSE ! 79 is out of level data
            nsta = 68 !79
           ENDIF

         Case(10021 ) 
           IF(br .le. 0.0589d0) then
            nsta = 0
           ELSEIF(br .le. 0.0640d0) then
            nsta = 1
           ELSEIF(br .le. 0.0749d0) then
            nsta = 3
           ELSEIF(br .le. 0.0788d0) then
            nsta = 4
           ELSEIF(br .le. 0.0811d0) then
            nsta = 6
           ELSEIF(br .le. 0.0834d0) then
            nsta = 7
           ELSEIF(br .le. 0.0915d0) then
            nsta = 11
           ELSEIF(br .le. 0.8362d0) then
            nsta = 12
           ELSEIF(br .le. 0.09723d0) then
            nsta = 17
           ELSEIF(br .le. 0.09919d0) then
            nsta = 21
           ELSE
            nsta = 23
           ENDIF

         Case(10022 ) 
           IF(br .le. 15.64d0 /85.94d0) then
            nsta = 1
           ELSEIF(br .le. 36.44d0 /85.94d0) then
            nsta = 7
           ELSE
            nsta = 13
           ENDIF

         Case(10023 ) 
           IF(br .le. 0.0247d0) then
            nsta = 1
           ELSEIF(br .le. 0.0432d0) then
            nsta = 3
           ELSEIF(br .le. 0.7100d0) then
            nsta = 6
           ELSEIF(br .le. 0.8642d0) then
            nsta = 10
           ELSE
            nsta = 16
           ENDIF


         Case(11024 ) 
           IF(br .le. 0.00003d0) then
            nsta = 0      
           ELSEIF(br .le. 0.00428d0) then
            nsta = 1      
           ELSEIF(br .le. 0.20195d0) then
            nsta = 2      
           ELSEIF(br .le. 0.24715d0) then
            nsta = 3      
           ELSEIF(br .le. 0.25202d0) then
            nsta = 4      
           ELSEIF(br .le. 0.25693d0) then
            nsta = 5      
           ELSEIF(br .le. 0.26183d0) then
            nsta = 7      
           ELSEIF(br .le. 0.26593d0) then
            nsta = 8      
           ELSEIF(br .le. 0.27028d0) then
            nsta = 11     
           ELSEIF(br .le. 0.27257d0) then
            nsta = 13     
           ELSEIF(br .le. 0.40120d0) then
            nsta = 14     
           ELSEIF(br .le. 0.51643d0) then
            nsta = 16     
           ELSEIF(br .le. 0.52527d0) then
            nsta = 17     
           ELSEIF(br .le. 0.55224d0) then
            nsta = 18     
           ELSEIF(br .le. 0.55452d0) then
            nsta = 19     
           ELSEIF(br .le. 0.55673d0) then
            nsta = 20     
           ELSEIF(br .le. 0.56408d0) then
            nsta = 21     
           ELSEIF(br .le. 0.57441d0) then
            nsta = 23     
           ELSEIF(br .le. 0.60133d0) then
            nsta = 28     
           ELSEIF(br .le. 0.60200d0) then
            nsta = 29     
           ELSEIF(br .le. 0.61224d0) then
            nsta = 32     
           ELSEIF(br .le. 0.71653d0) then
            nsta = 33     
           ELSEIF(br .le. 0.86255d0) then
            nsta = 35     
           ELSEIF(br .le. 0.87547d0) then
            nsta = 38     
           ELSEIF(br .le. 0.87703d0) then
            nsta = 39     
           ELSEIF(br .le. 0.87762d0) then
            nsta = 40     
           ELSEIF(br .le. 0.93186d0) then
            nsta = 41     
           ELSEIF(br .le. 0.94378d0) then
            nsta = 48     
           ELSEIF(br .le. 0.95614d0) then
            nsta = 49     
           ELSEIF(br .le. 0.96061d0) then
            nsta = 50     
           ELSEIF(br .le. 0.96114d0) then
            nsta = 53     
           ELSEIF(br .le. 0.96710d0) then
            nsta = 56     
           ELSEIF(br .le. 0.97058d0) then
            nsta = 61     
           ELSEIF(br .le. 0.97972d0) then
            nsta = 69     
           ELSEIF(br .le. 0.98009d0) then
            nsta = 71     
           ELSEIF(br .le. 0.98297d0) then
            nsta = 73     
           ELSEIF(br .le. 0.99087d0) then
            nsta = 76     
           ELSEIF(br .le. 0.99106d0) then
            nsta = 80
           ELSE
            nsta = 81
           ENDIF


         Case(12025 ) 
           IF(br .le. 0.0004d0) then
            nsta = 0
           ELSEIF(br .le. 0.0037d0) then
            nsta = 1
           ELSEIF(br .le. 0.0278d0) then
            nsta = 2
           ELSEIF(br .le. 0.0354d0) then
            nsta = 5
           ELSEIF(br .le. 0.0438d0) then
            nsta = 7
           ELSEIF(br .le. 0.7993d0) then
            nsta = 9
           ELSEIF(br .le. 0.9910d0) then
            nsta = 13
           ELSEIF(br .le. 0.9926d0) then
            nsta = 14
           ELSE
            nsta = 18
           ENDIF


         Case(12026 ) 
           IF(br .le. 0.001407672d0) then
            nsta = 0      
           ELSEIF(br .le. 0.024282339d0) then
            nsta = 1      
           ELSEIF(br .le. 0.174098839d0) then
            nsta = 2      
           ELSEIF(br .le. 0.185611583d0) then
            nsta = 4      
           ELSEIF(br .le. 0.189080489d0) then
            nsta = 5      
           ELSEIF(br .le. 0.200593233d0) then
            nsta = 6      
           ELSEIF(br .le. 0.221758584d0) then
            nsta = 7      
           ELSEIF(br .le. 0.225177216d0) then
            nsta = 8      
           ELSEIF(br .le. 0.227841737d0) then
            nsta = 9      
           ELSEIF(br .le. 0.230355437d0) then
            nsta = 10     
           ELSEIF(br .le. 0.246191745d0) then
            nsta = 11     
           ELSEIF(br .le. 0.248906541d0) then
            nsta = 12     
           ELSEIF(br .le. 0.25116887d0) then
            nsta = 13     
           ELSEIF(br .le. 0.251822432d0) then
            nsta = 14     
           ELSEIF(br .le. 0.337288221d0) then
            nsta = 15     
           ELSEIF(br .le. 0.342013976d0) then
            nsta = 18     
           ELSEIF(br .le. 0.345382334d0) then
            nsta = 19     
           ELSEIF(br .le. 0.420793324d0) then
            nsta = 20     
           ELSEIF(br .le. 0.42994319d0) then
            nsta = 22     
           ELSEIF(br .le. 0.438791413d0) then
            nsta = 23     
           ELSEIF(br .le. 0.657986024d0) then
            nsta = 26     
           ELSEIF(br .le. 0.708762757d0) then
            nsta = 27     
           ELSEIF(br .le. 0.780654567d0) then
            nsta = 28     
           ELSEIF(br .le. 0.784525665d0) then
            nsta = 29     
           ELSEIF(br .le. 0.840832537d0) then
            nsta = 32     
           ELSEIF(br .le. 0.842039113d0) then
            nsta = 34     
           ELSEIF(br .le. 0.846412951d0) then
            nsta = 35     
           ELSEIF(br .le. 0.847770348d0) then
            nsta = 43     
           ELSEIF(br .le. 0.861847067d0) then
            nsta = 44     
           ELSEIF(br .le. 0.870393645d0) then
            nsta = 46     
           ELSEIF(br .le. 0.882208034d0) then
            nsta = 47     
           ELSEIF(br .le. 0.886330501d0) then
            nsta = 49     
           ELSEIF(br .le. 0.891307627d0) then
            nsta = 52     
           ELSEIF(br .le. 0.899351465d0) then
            nsta = 53     
           ELSEIF(br .le. 0.905082701d0) then
            nsta = 57     
           ELSEIF(br .le. 0.936051481d0) then
            nsta = 59     
           ELSEIF(br .le. 0.937660248d0) then
            nsta = 61     
           ELSEIF(br .le. 0.938716002d0) then
            nsta = 63     
           ELSEIF(br .le. 0.939922578d0) then
            nsta = 68     
           ELSEIF(br .le. 0.955708612d0) then
            nsta = 74     
           ELSEIF(br .le. 0.961842039d0) then
            nsta = 77     
           ELSEIF(br .le. 0.975013825d0) then
            nsta = 81     
           ELSEIF(br .le. 0.980845609d0) then
            nsta = 91     
           ELSEIF(br .le. 0.982353828d0) then
            nsta = 101    
           ELSEIF(br .le. 0.984314514d0) then
            nsta = 102    
           ELSEIF(br .le. 0.984716706d0) then
            nsta = 107    
           ELSEIF(br .le. 0.985822734d0) then
            nsta = 113    
           ELSEIF(br .le. 0.993816299d0) then
            nsta = 114    
           ELSEIF(br .le. 0.997536574d0) then
            nsta = 124    
           ELSEIF(br .le. 0.997888492d0) then
            nsta = 127    
           ELSEIF(br .le. 0.999044794d0) then
            nsta = 130    
           ELSEIF(br .le. 0.999648082d0) then
            nsta = 132    
           ELSE
            nsta = 134
           ENDIF


         Case(12027 ) 
           IF(br .le. 0.090745715d0) then
            nsta = 0      
           ELSEIF(br .le. 0.115263844d0) then
            nsta = 1      
           ELSEIF(br .le. 0.136750533d0) then
            nsta = 6      
           ELSEIF(br .le. 0.139278379d0) then
            nsta = 7      
           ELSEIF(br .le. 0.786406904d0) then
            nsta = 8      
           ELSEIF(br .le. 0.822803934d0) then
            nsta = 10     
           ELSEIF(br .le. 0.989641757d0) then
            nsta = 16     
           ELSEIF(br .le. 0.99393712 d0) then
            nsta = 18
           ELSE
            nsta = 30
           ENDIF


         Case(13028 ) 
           IF(br .le. 0.262946253d0) then
            nsta = 0      
           ELSEIF(br .le. 0.295017654d0) then
            nsta = 1      
           ELSEIF(br .le. 0.29560612 d0) then
            nsta = 2      
           ELSEIF(br .le. 0.304433111d0) then
            nsta = 3      
           ELSEIF(br .le. 0.305511965d0) then
            nsta = 4      
           ELSEIF(br .le. 0.330816006d0) then
            nsta = 6      
           ELSEIF(br .le. 0.341604551d0) then
            nsta = 7      
           ELSEIF(br .le. 0.342193017d0) then
            nsta = 8      
           ELSEIF(br .le. 0.343860337d0) then
            nsta = 9      
           ELSEIF(br .le. 0.346410357d0) then
            nsta = 10     
           ELSEIF(br .le. 0.346606512d0) then
            nsta = 11     
           ELSEIF(br .le. 0.348273833d0) then
            nsta = 12     
           ELSEIF(br .le. 0.352687328d0) then
            nsta = 13     
           ELSEIF(br .le. 0.36063162 d0) then
            nsta = 16     
           ELSEIF(br .le. 0.36484896 d0) then
            nsta = 17     
           ELSEIF(br .le. 0.431541781d0) then
            nsta = 18     
           ELSEIF(br .le. 0.499215379d0) then
            nsta = 20     
           ELSEIF(br .le. 0.500588466d0) then
            nsta = 22     
           ELSEIF(br .le. 0.507748137d0) then
            nsta = 23     
           ELSEIF(br .le. 0.538348372d0) then
            nsta = 25     
           ELSEIF(br .le. 0.543840722d0) then
            nsta = 26     
           ELSEIF(br .le. 0.55237348 d0) then
            nsta = 27     
           ELSEIF(br .le. 0.553550412d0) then
            nsta = 30     
           ELSEIF(br .le. 0.554531189d0) then
            nsta = 33     
           ELSEIF(br .le. 0.556885053d0) then
            nsta = 35     
           ELSEIF(br .le. 0.643193409d0) then
            nsta = 36     
           ELSEIF(br .le. 0.73734798 d0) then
            nsta = 38     
           ELSEIF(br .le. 0.775304041d0) then
            nsta = 40     
           ELSEIF(br .le. 0.782071401d0) then
            nsta = 44     
           ELSEIF(br .le. 0.82346018 d0) then
            nsta = 45     
           ELSEIF(br .le. 0.824931346d0) then
            nsta = 47     
           ELSEIF(br .le. 0.826500588d0) then
            nsta = 52     
           ELSEIF(br .le. 0.873087485d0) then
            nsta = 54     
           ELSEIF(br .le. 0.883581797d0) then
            nsta = 56     
           ELSEIF(br .le. 0.895351118d0) then
            nsta = 58     
           ELSEIF(br .le. 0.899862691d0) then
            nsta = 60     
           ELSEIF(br .le. 0.903687721d0) then
            nsta = 67     
           ELSEIF(br .le. 0.921341703d0) then
            nsta = 71     
           ELSEIF(br .le. 0.951745783d0) then
            nsta = 73     
           ELSEIF(br .le. 0.953609259d0) then
            nsta = 75     
           ELSEIF(br .le. 0.964790114d0) then
            nsta = 76     
           ELSEIF(br .le. 0.971949784d0) then
            nsta = 86     
           ELSEIF(br .le. 0.977540212d0) then
            nsta = 87     
           ELSEIF(br .le. 0.981659474d0) then
            nsta = 90     
           ELSEIF(br .le. 0.994899961d0) then
            nsta = 97     
           ELSEIF(br .le. 0.997155747d0) then
            nsta = 105
           ELSE
            nsta = 109
           ENDIF


         Case(14029 ) 
           IF(br .le. 0.021609709d0) then
            nsta = 0      
           ELSEIF(br .le. 0.0918637  d0) then
            nsta = 1      
           ELSEIF(br .le. 0.093040417d0) then
            nsta = 2      
           ELSEIF(br .le. 0.096291347d0) then
            nsta = 3      
           ELSEIF(br .le. 0.096650346d0) then
            nsta = 4      
           ELSEIF(br .le. 0.79621855 d0) then
            nsta = 10     
           ELSEIF(br .le. 0.991035012d0) then
            nsta = 18     
           ELSEIF(br .le. 0.991443872d0) then
            nsta = 25     
           ELSEIF(br .le. 0.996579543d0) then
            nsta = 27     
           ELSEIF(br .le. 0.998703617d0) then
            nsta = 30     
           ELSEIF(br .le. 0.999411642d0) then
            nsta = 35
           ELSE
            nsta = 42
           ENDIF


         Case(14030 ) 
           IF(br .le. 0.0685467448d0) then
            nsta = 0      
           ELSEIF(br .le. 0.0753117616d0) then
            nsta = 1      
           ELSEIF(br .le. 0.1264874894d0) then
            nsta = 2      
           ELSEIF(br .le. 0.1639701694d0) then
            nsta = 3      
           ELSEIF(br .le. 0.1749225693d0) then
            nsta = 4      
           ELSEIF(br .le. 0.1776632162d0) then
            nsta = 11     
           ELSEIF(br .le. 0.2178865432d0) then
            nsta = 15     
           ELSEIF(br .le. 0.5217010351d0) then
            nsta = 17     
           ELSEIF(br .le. 0.5329081423d0) then
            nsta = 19     
           ELSEIF(br .le. 0.774248105 d0) then
            nsta = 26     
           ELSEIF(br .le. 0.7761329379d0) then
            nsta = 30     
           ELSEIF(br .le. 0.7787003838d0) then
            nsta = 34     
           ELSEIF(br .le. 0.858719133 d0) then
            nsta = 35     
           ELSEIF(br .le. 0.9015099031d0) then
            nsta = 50     
           ELSEIF(br .le. 0.904250551 d0) then
            nsta = 51     
           ELSEIF(br .le. 0.911351781 d0) then
            nsta = 52     
           ELSEIF(br .le. 0.949944983 d0) then
            nsta = 57     
           ELSEIF(br .le. 0.965431168 d0) then
            nsta = 62     
           ELSEIF(br .le. 0.973816122 d0) then
            nsta = 71     
           ELSEIF(br .le. 0.984595322 d0) then
            nsta = 73     
           ELSEIF(br .le. 0.995720923 d0) then
            nsta = 80     
           ELSEIF(br .le. 0.998889478 d0) then
            nsta = 92
           ELSE
            nsta = 94
           ENDIF


         Case(14031) 
            nsta = 37

         Case(15032) 
            nsta = 113

         Case(16033) 
            nsta = 116

         Case(16034) 
            nsta = 172

         Case(16035) 
            nsta = 56

         Case(16037 )
            nsta = 23
      
         Case(17036) 
            nsta = 145

         Case(17038 ) 
           IF(br .le. 0.036144578d0) then
            nsta = 0
           ELSEIF(br .le. 0.119386637d0) then
            nsta = 2
           ELSEIF(br .le. 0.284775465d0) then
            nsta = 4
           ELSEIF(br .le. 0.424972618d0) then
            nsta = 5
           ELSEIF(br .le. 0.47645126 d0) then
            nsta = 6
           ELSEIF(br .le. 0.487404162d0) then
            nsta = 9
           ELSEIF(br .le. 0.671412924d0) then
            nsta = 10
           ELSEIF(br .le. 0.764512596d0) then
            nsta = 12
           ELSEIF(br .le. 0.775465498d0) then
            nsta = 14
           ELSEIF(br .le. 0.779846659d0) then
            nsta = 19
           ELSEIF(br .le. 0.800657174d0) then
            nsta = 21
           ELSEIF(br .le. 0.812705367d0) then
            nsta = 22
           ELSEIF(br .le. 0.837897043d0) then
            nsta = 23
           ELSEIF(br .le. 0.851040526d0) then
            nsta = 24
           ELSEIF(br .le. 0.91456736 d0) then
            nsta = 25
           ELSEIF(br .le. 0.946330778d0) then
            nsta = 27
           ELSEIF(br .le. 0.966046002d0) then
            nsta = 28
           ELSEIF(br .le. 0.976998905d0) then
            nsta = 30
           ELSEIF(br .le. 0.986856517d0) then
            nsta = 31
           ELSEIF(br .le. 0.992332968d0) then
            nsta = 33
           ELSE                                           
            nsta = 37
           ENDIF

         Case(18037 ) 
            nsta = 104
      
         Case(18041 ) 
           IF(br .le. 0.128383758d0) then
            nsta = 3      
           ELSEIF(br .le. 0.133658439d0) then
            nsta = 4      
           ELSEIF(br .le. 0.685509554d0) then
            nsta = 5      
           ELSEIF(br .le. 0.778065287d0) then
            nsta = 10     
           ELSEIF(br .le. 0.816878981d0) then
            nsta = 13     
           ELSEIF(br .le. 0.859175955d0) then
            nsta = 15     
           ELSEIF(br .le. 0.868630573d0) then
            nsta = 16     
           ELSEIF(br .le. 0.951731688d0) then
            nsta = 19     
           ELSEIF(br .le. 0.990047771d0) then
            nsta = 31
           ELSE
            nsta = 36
           ENDIF


         Case(19040 ) 
           IF(br .le. 0.062613907d0) then
            nsta = 1      
           ELSEIF(br .le. 0.086653175d0) then
            nsta = 2      
           ELSEIF(br .le. 0.088889386d0) then
            nsta = 5      
           ELSEIF(br .le. 0.150385187d0) then
            nsta = 6      
           ELSEIF(br .le. 0.175877992d0) then
            nsta = 7      
           ELSEIF(br .le. 0.238491899d0) then
            nsta = 8      
           ELSEIF(br .le. 0.273935843d0) then
            nsta = 10     
           ELSEIF(br .le. 0.362266176d0) then
            nsta = 13     
           ELSEIF(br .le. 0.366481434d0) then
            nsta = 16     
           ELSEIF(br .le. 0.39219786 d0) then
            nsta = 17     
           ELSEIF(br .le. 0.406174178d0) then
            nsta = 18     
           ELSEIF(br .le. 0.426076456d0) then
            nsta = 20     
           ELSEIF(br .le. 0.43915829 d0) then
            nsta = 21     
           ELSEIF(br .le. 0.46353299 d0) then
            nsta = 23     
           ELSEIF(br .le. 0.464114405d0) then
            nsta = 29     
           ELSEIF(br .le. 0.471493901d0) then
            nsta = 30     
           ELSEIF(br .le. 0.477308049d0) then
            nsta = 31     
           ELSEIF(br .le. 0.483904872d0) then
            nsta = 35     
           ELSEIF(br .le. 0.488600915d0) then
            nsta = 36     
           ELSEIF(br .le. 0.505037065d0) then
            nsta = 37     
           ELSEIF(br .le. 0.553451033d0) then
            nsta = 38     
           ELSEIF(br .le. 0.556581728d0) then
            nsta = 40     
           ELSEIF(br .le. 0.561613203d0) then
            nsta = 42     
           ELSEIF(br .le. 0.586546955d0) then
            nsta = 43     
           ELSEIF(br .le. 0.594485504d0) then
            nsta = 44     
           ELSEIF(br .le. 0.632612901d0) then
            nsta = 45     
           ELSEIF(br .le. 0.649719915d0) then
            nsta = 46     
           ELSEIF(br .le. 0.652190928d0) then
            nsta = 48     
           ELSEIF(br .le. 0.670192426d0) then
            nsta = 49     
           ELSEIF(br .le. 0.684615987d0) then
            nsta = 50     
           ELSEIF(br .le. 0.701163948d0) then
            nsta = 51     
           ELSEIF(br .le. 0.718606393d0) then
            nsta = 52     
           ELSEIF(br .le. 0.729340206d0) then
            nsta = 53     
           ELSEIF(br .le. 0.732470901d0) then
            nsta = 55     
           ELSEIF(br .le. 0.742869282d0) then
            nsta = 57     
           ELSEIF(br .le. 0.758858191d0) then
            nsta = 59     
           ELSEIF(br .le. 0.775517962d0) then
            nsta = 60     
           ELSEIF(br .le. 0.800339904d0) then
            nsta = 61     
           ELSEIF(br .le. 0.808949316d0) then
            nsta = 62     
           ELSEIF(br .le. 0.811375605d0) then
            nsta = 63     
           ELSEIF(br .le. 0.863926563d0) then
            nsta = 64     
           ELSEIF(br .le. 0.87566667 d0) then
            nsta = 65     
           ELSEIF(br .le. 0.886847725d0) then
            nsta = 70     
           ELSEIF(br .le. 0.889307557d0) then
            nsta = 71     
           ELSEIF(br .le. 0.898140591d0) then
            nsta = 72     
           ELSEIF(br .le. 0.925310554d0) then
            nsta = 74     
           ELSEIF(br .le. 0.957288371d0) then
            nsta = 78     
           ELSEIF(br .le. 0.962878898d0) then
            nsta = 80     
           ELSEIF(br .le. 0.971041068d0) then
            nsta = 83     
           ELSEIF(br .le. 0.990719724d0) then
            nsta = 87
           ELSE
            nsta = 89
           ENDIF

         Case(19041 ) 
           IF(br .le. 0.119  d0/93.383d0) then
            nsta = 0
           ELSEIF(br .le. 0.42   d0/93.383d0) then
            nsta = 1
           ELSEIF(br .le. 0.462  d0/93.383d0) then
            nsta = 3
           ELSEIF(br .le. 0.608  d0/93.383d0) then
            nsta = 4
           ELSEIF(br .le. 0.763  d0/93.383d0) then
            nsta = 7
           ELSEIF(br .le. 1.103  d0/93.383d0) then
            nsta = 8
           ELSEIF(br .le. 1.241  d0/93.383d0) then
            nsta = 9
           ELSEIF(br .le. 1.292  d0/93.383d0) then
            nsta = 11
           ELSEIF(br .le. 1.355  d0/93.383d0) then
            nsta = 12
           ELSEIF(br .le. 1.538  d0/93.383d0) then
            nsta = 13
           ELSEIF(br .le. 2.002  d0/93.383d0) then
            nsta = 15
           ELSEIF(br .le. 2.134  d0/93.383d0) then
            nsta = 16
           ELSEIF(br .le. 2.464  d0/93.383d0) then
            nsta = 17
           ELSEIF(br .le. 2.691  d0/93.383d0) then
            nsta = 19
           ELSEIF(br .le. 2.992  d0/93.383d0) then
            nsta = 24
           ELSEIF(br .le. 3.325  d0/93.383d0) then
            nsta = 25
           ELSEIF(br .le. 3.445  d0/93.383d0) then
            nsta = 26
           ELSEIF(br .le. 3.642  d0/93.383d0) then
            nsta = 27
           ELSEIF(br .le. 3.771  d0/93.383d0) then
            nsta = 28
           ELSEIF(br .le. 3.896  d0/93.383d0) then
            nsta = 31
           ELSEIF(br .le. 4.024  d0/93.383d0) then
            nsta = 34
           ELSEIF(br .le. 4.643  d0/93.383d0) then
            nsta = 35
           ELSEIF(br .le. 5.753  d0/93.383d0) then
            nsta = 37
           ELSEIF(br .le. 5.907  d0/93.383d0) then
            nsta = 38
           ELSEIF(br .le. 6.021  d0/93.383d0) then
            nsta = 39
           ELSEIF(br .le. 6.119  d0/93.383d0) then
            nsta = 40
           ELSEIF(br .le. 6.352  d0/93.383d0) then
            nsta = 41
           ELSEIF(br .le. 6.598  d0/93.383d0) then
            nsta = 44
           ELSEIF(br .le. 7.029  d0/93.383d0) then
            nsta = 45
           ELSEIF(br .le. 7.39   d0/93.383d0) then
            nsta = 46
           ELSEIF(br .le. 7.925  d0/93.383d0) then
            nsta = 47
           ELSEIF(br .le. 8.72   d0/93.383d0) then
            nsta = 48
           ELSEIF(br .le. 8.933  d0/93.383d0) then
            nsta = 51
           ELSEIF(br .le. 9.793  d0/93.383d0) then
            nsta = 52
           ELSEIF(br .le. 11.052 d0/93.383d0) then
            nsta = 53
           ELSEIF(br .le. 11.334 d0/93.383d0) then
            nsta = 54
           ELSEIF(br .le. 11.454 d0/93.383d0) then
            nsta = 55
           ELSEIF(br .le. 11.62  d0/93.383d0) then
            nsta = 56
           ELSEIF(br .le. 11.835 d0/93.383d0) then
            nsta = 57
           ELSEIF(br .le. 12.487 d0/93.383d0) then
            nsta = 60
           ELSEIF(br .le. 14.557 d0/93.383d0) then
            nsta = 61
           ELSEIF(br .le. 14.789 d0/93.383d0) then
            nsta = 62
           ELSEIF(br .le. 15.039 d0/93.383d0) then
            nsta = 63
           ELSEIF(br .le. 15.199 d0/93.383d0) then
            nsta = 64
           ELSEIF(br .le. 15.671 d0/93.383d0) then
            nsta = 66
           ELSEIF(br .le. 17.116 d0/93.383d0) then
            nsta = 67
           ELSEIF(br .le. 17.455 d0/93.383d0) then
            nsta = 68
           ELSEIF(br .le. 17.791 d0/93.383d0) then
            nsta = 69
           ELSEIF(br .le. 17.942 d0/93.383d0) then
            nsta = 73
           ELSEIF(br .le. 18.671 d0/93.383d0) then
            nsta = 74
           ELSEIF(br .le. 19.537 d0/93.383d0) then
            nsta = 75
           ELSEIF(br .le. 19.654 d0/93.383d0) then
            nsta = 76
           ELSEIF(br .le. 19.855 d0/93.383d0) then
            nsta = 77
           ELSEIF(br .le. 20.622 d0/93.383d0) then
            nsta = 78
           ELSEIF(br .le. 20.724 d0/93.383d0) then
            nsta = 79
           ELSEIF(br .le. 21.687 d0/93.383d0) then
            nsta = 80
           ELSEIF(br .le. 22.208 d0/93.383d0) then
            nsta = 81
           ELSEIF(br .le. 23.009 d0/93.383d0) then
            nsta = 82
           ELSEIF(br .le. 23.359 d0/93.383d0) then
            nsta = 83
           ELSEIF(br .le. 23.463 d0/93.383d0) then
            nsta = 85
           ELSEIF(br .le. 23.721 d0/93.383d0) then
            nsta = 86
           ELSEIF(br .le. 24.309 d0/93.383d0) then
            nsta = 89
           ELSEIF(br .le. 26.027 d0/93.383d0) then
            nsta = 90
           ELSEIF(br .le. 26.348 d0/93.383d0) then
            nsta = 91
           ELSEIF(br .le. 26.574 d0/93.383d0) then
            nsta = 92
           ELSEIF(br .le. 26.784 d0/93.383d0) then
            nsta = 94
           ELSEIF(br .le. 27.003 d0/93.383d0) then
            nsta = 95
           ELSEIF(br .le. 27.123 d0/93.383d0) then
            nsta = 96
           ELSEIF(br .le. 28.037 d0/93.383d0) then
            nsta = 97
           ELSEIF(br .le. 28.527 d0/93.383d0) then
            nsta = 98
           ELSEIF(br .le. 28.987 d0/93.383d0) then
            nsta = 99
           ELSEIF(br .le. 29.814 d0/93.383d0) then
            nsta = 100
           ELSEIF(br .le. 30.249 d0/93.383d0) then
            nsta = 101
           ELSEIF(br .le. 31.506 d0/93.383d0) then
            nsta = 102
           ELSEIF(br .le. 32.144 d0/93.383d0) then
            nsta = 103
           ELSEIF(br .le. 33.036 d0/93.383d0) then
            nsta = 104
           ELSEIF(br .le. 33.553 d0/93.383d0) then
            nsta = 105
           ELSEIF(br .le. 33.993 d0/93.383d0) then
            nsta = 108
           ELSEIF(br .le. 34.267 d0/93.383d0) then
            nsta = 109
           ELSEIF(br .le. 34.845 d0/93.383d0) then
            nsta = 110
           ELSEIF(br .le. 35.31  d0/93.383d0) then
            nsta = 111
           ELSEIF(br .le. 39.42  d0/93.383d0) then
            nsta = 112
           ELSEIF(br .le. 40.881 d0/93.383d0) then
            nsta = 113
           ELSEIF(br .le. 40.943 d0/93.383d0) then
            nsta = 114
           ELSEIF(br .le. 41.886 d0/93.383d0) then
            nsta = 115
           ELSEIF(br .le. 42.473 d0/93.383d0) then
            nsta = 116
           ELSEIF(br .le. 43.065 d0/93.383d0) then
            nsta = 117
           ELSEIF(br .le. 43.433 d0/93.383d0) then
            nsta = 118
           ELSEIF(br .le. 45.443 d0/93.383d0) then
            nsta = 119
           ELSEIF(br .le. 47.453 d0/93.383d0) then
            nsta = 120
           ELSEIF(br .le. 51.143 d0/93.383d0) then
            nsta = 121
           ELSEIF(br .le. 52.173 d0/93.383d0) then
            nsta = 122
           ELSEIF(br .le. 58.733 d0/93.383d0) then
            nsta = 124
           ELSEIF(br .le. 59.273 d0/93.383d0) then
            nsta = 125
           ELSEIF(br .le. 61.923 d0/93.383d0) then
            nsta = 126
           ELSEIF(br .le. 62.048 d0/93.383d0) then
            nsta = 128
           ELSEIF(br .le. 62.106 d0/93.383d0) then
            nsta = 129
           ELSEIF(br .le. 62.489 d0/93.383d0) then
            nsta = 130
           ELSEIF(br .le. 65.009 d0/93.383d0) then
            nsta = 131
           ELSEIF(br .le. 66.827 d0/93.383d0) then
            nsta = 132
           ELSEIF(br .le. 68.192 d0/93.383d0) then
            nsta = 133
           ELSEIF(br .le. 71.902 d0/93.383d0) then
            nsta = 134
           ELSEIF(br .le. 74.162 d0/93.383d0) then
            nsta = 135
           ELSEIF(br .le. 74.349 d0/93.383d0) then
            nsta = 136
           ELSEIF(br .le. 74.504 d0/93.383d0) then
            nsta = 137
           ELSEIF(br .le. 75.137 d0/93.383d0) then
            nsta = 138
           ELSEIF(br .le. 75.447 d0/93.383d0) then
            nsta = 139
           ELSEIF(br .le. 75.817 d0/93.383d0) then
            nsta = 140
           ELSEIF(br .le. 76.217 d0/93.383d0) then
            nsta = 142
           ELSEIF(br .le. 77.118 d0/93.383d0) then
            nsta = 143
           ELSEIF(br .le. 77.487 d0/93.383d0) then
            nsta = 144
           ELSEIF(br .le. 78.272 d0/93.383d0) then
            nsta = 147
           ELSEIF(br .le. 79.808 d0/93.383d0) then
            nsta = 151
           ELSEIF(br .le. 80.014 d0/93.383d0) then
            nsta = 154
           ELSEIF(br .le. 80.649 d0/93.383d0) then
            nsta = 155
           ELSEIF(br .le. 81.434 d0/93.383d0) then
            nsta = 159
           ELSEIF(br .le. 87.064 d0/93.383d0) then
            nsta = 161
           ELSEIF(br .le. 88.316 d0/93.383d0) then
            nsta = 164
           ELSEIF(br .le. 91.526 d0/93.383d0) then
            nsta = 167
           ELSE
            nsta = 169
           ENDIF
      
         Case(19042 ) 
           IF(br .le. 0.019885846d0) then
            nsta = 0
           ELSEIF(br .le. 0.054260411d0) then
            nsta = 1
           ELSEIF(br .le. 0.067667563d0) then
            nsta = 4
           ELSEIF(br .le. 0.127956909d0) then
            nsta = 5
           ELSEIF(br .le. 0.13636315 d0) then
            nsta = 8
           ELSEIF(br .le. 0.143163102d0) then
            nsta = 9
           ELSEIF(br .le. 0.145369071d0) then
            nsta = 10
           ELSEIF(br .le. 0.161817461d0) then
            nsta = 12
           ELSEIF(br .le. 0.170223702d0) then
            nsta = 13
           ELSEIF(br .le. 0.174175171d0) then
            nsta = 14
           ELSEIF(br .le. 0.183823608d0) then
            nsta = 16
           ELSEIF(br .le. 0.188107043d0) then
            nsta = 17
           ELSEIF(br .le. 0.19206922 d0) then
            nsta = 18
           ELSEIF(br .le. 0.195388882d0) then
            nsta = 20
           ELSEIF(br .le. 0.202167418d0) then
            nsta = 25
           ELSEIF(br .le. 0.203827249d0) then
            nsta = 26
           ELSEIF(br .le. 0.205829755d0) then
            nsta = 27
           ELSEIF(br .le. 0.230031162d0) then
            nsta = 28
           ELSEIF(br .le. 0.269760021d0) then
            nsta = 29
           ELSEIF(br .le. 0.284377242d0) then
            nsta = 32
           ELSEIF(br .le. 0.303845454d0) then
            nsta = 35
           ELSEIF(br .le. 0.330831093d0) then
            nsta = 36
           ELSEIF(br .le. 0.334932482d0) then
            nsta = 38
           ELSEIF(br .le. 0.33555358 d0) then
            nsta = 39
           ELSEIF(br .le. 0.336892154d0) then
            nsta = 40
           ELSEIF(br .le. 0.36526991 d0) then
            nsta = 41
           ELSEIF(br .le. 0.371052547d0) then
            nsta = 42
           ELSEIF(br .le. 0.441300879d0) then
            nsta = 45
           ELSEIF(br .le. 0.452330724d0) then
            nsta = 46
           ELSEIF(br .le. 0.491845411d0) then
            nsta = 47
           ELSEIF(br .le. 0.513369671d0) then
            nsta = 48
           ELSEIF(br .le. 0.534893931d0) then
            nsta = 49
           ELSEIF(br .le. 0.538834692d0) then
            nsta = 52
           ELSEIF(br .le. 0.545174175d0) then
            nsta = 53
           ELSEIF(br .le. 0.551460116d0) then
            nsta = 54
           ELSEIF(br .le. 0.561558314d0) then
            nsta = 55
           ELSEIF(br .le. 0.562222246d0) then
            nsta = 56
           ELSEIF(br .le. 0.577867492d0) then
            nsta = 57
           ELSEIF(br .le. 0.621879785d0) then
            nsta = 58
           ELSEIF(br .le. 0.626859278d0) then
            nsta = 61
           ELSEIF(br .le. 0.633048842d0) then
            nsta = 62
           ELSEIF(br .le. 0.635982995d0) then
            nsta = 63
           ELSEIF(br .le. 0.640694773d0) then
            nsta = 64
           ELSEIF(br .le. 0.646231113d0) then
            nsta = 65
           ELSEIF(br .le. 0.655783173d0) then
            nsta = 67
           ELSEIF(br .le. 0.662615251d0) then
            nsta = 68
           ELSEIF(br .le. 0.676075945d0) then
            nsta = 69
           ELSEIF(br .le. 0.680734181d0) then
            nsta = 70
           ELSEIF(br .le. 0.689590182d0) then
            nsta = 74
           ELSEIF(br .le. 0.694516132d0) then
            nsta = 75
           ELSEIF(br .le. 0.69976334 d0) then
            nsta = 76
           ELSEIF(br .le. 0.709550989d0) then
            nsta = 77
           ELSEIF(br .le. 0.710836019d0) then
            nsta = 78
           ELSEIF(br .le. 0.7131812  d0) then
            nsta = 79
           ELSEIF(br .le. 0.715430003d0) then
            nsta = 80
           ELSEIF(br .le. 0.717850144d0) then
            nsta = 81
           ELSEIF(br .le. 0.721287601d0) then
            nsta = 82
           ELSEIF(br .le. 0.739684953d0) then
            nsta = 83
           ELSEIF(br .le. 0.745981603d0) then
            nsta = 84
           ELSEIF(br .le. 0.748744418d0) then
            nsta = 86
           ELSEIF(br .le. 0.749858111d0) then
            nsta = 87
           ELSEIF(br .le. 0.753606117d0) then
            nsta = 90
           ELSEIF(br .le. 0.762183695d0) then
            nsta = 91
           ELSEIF(br .le. 0.767762869d0) then
            nsta = 92
           ELSEIF(br .le. 0.778075239d0) then
            nsta = 93
           ELSEIF(br .le. 0.779167514d0) then
            nsta = 95
           ELSEIF(br .le. 0.787381001d0) then
            nsta = 96
           ELSEIF(br .le. 0.789533427d0) then
            nsta = 97
           ELSEIF(br .le. 0.790786332d0) then
            nsta = 98
           ELSEIF(br .le. 0.800059968d0) then
            nsta = 99
           ELSEIF(br .le. 0.807866528d0) then
            nsta = 100
           ELSEIF(br .le. 0.809483525d0) then
            nsta = 101
           ELSEIF(br .le. 0.81308161 d0) then
            nsta = 102
           ELSEIF(br .le. 0.816711821d0) then
            nsta = 103
           ELSEIF(br .le. 0.83218573 d0) then
            nsta = 104
           ELSEIF(br .le. 0.837240183d0) then
            nsta = 105
           ELSEIF(br .le. 0.838953557d0) then
            nsta = 106
           ELSEIF(br .le. 0.841630704d0) then
            nsta = 109
           ELSEIF(br .le. 0.844115096d0) then
            nsta = 110
           ELSEIF(br .le. 0.866281871d0) then
            nsta = 111
           ELSEIF(br .le. 0.87326387 d0) then
            nsta = 112
           ELSEIF(br .le. 0.875566217d0) then
            nsta = 113
           ELSEIF(br .le. 0.877343842d0) then
            nsta = 114
           ELSEIF(br .le. 0.878628872d0) then
            nsta = 115
           ELSEIF(br .le. 0.881648694d0) then
            nsta = 116
           ELSEIF(br .le. 0.895130805d0) then
            nsta = 117
           ELSEIF(br .le. 0.90434019 d0) then
            nsta = 118
           ELSEIF(br .le. 0.906621119d0) then
            nsta = 119
           ELSEIF(br .le. 0.915134446d0) then
            nsta = 120
           ELSEIF(br .le. 0.92086354 d0) then
            nsta = 121
           ELSEIF(br .le. 0.92472934 d0) then
            nsta = 122
           ELSEIF(br .le. 0.929344742d0) then
            nsta = 123
           ELSEIF(br .le. 0.931979054d0) then
            nsta = 125
           ELSEIF(br .le. 0.934474155d0) then
            nsta = 126
           ELSEIF(br .le. 0.935523596d0) then
            nsta = 127
           ELSEIF(br .le. 0.936744375d0) then
            nsta = 129
           ELSEIF(br .le. 0.938393498d0) then
            nsta = 130
           ELSEIF(br .le. 0.95028003 d0) then
            nsta = 131
           ELSEIF(br .le. 0.956908645d0) then
            nsta = 132
           ELSEIF(br .le. 0.958279344d0) then
            nsta = 133
           ELSEIF(br .le. 0.959617918d0) then
            nsta = 134
           ELSEIF(br .le. 0.960999325d0) then
            nsta = 135
           ELSEIF(br .le. 0.963108917d0) then
            nsta = 136
           ELSEIF(br .le. 0.964393948d0) then
            nsta = 137
           ELSEIF(br .le. 0.967959907d0) then
            nsta = 138
           ELSEIF(br .le. 0.971183192d0) then
            nsta = 139
           ELSEIF(br .le. 0.973614041d0) then
            nsta = 140
           ELSEIF(br .le. 0.977147875d0) then
            nsta = 141
           ELSEIF(br .le. 0.978561408d0) then
            nsta = 143
           ELSEIF(br .le. 0.983530193d0) then
            nsta = 144
           ELSEIF(br .le. 0.985489864d0) then
            nsta = 145
           ELSEIF(br .le. 0.986164505d0) then
            nsta = 148
           ELSEIF(br .le. 0.986710643d0) then
            nsta = 152
           ELSEIF(br .le. 0.988188428d0) then
            nsta = 155
           ELSEIF(br .le. 0.991829348d0) then
            nsta = 156
           ELSEIF(br .le. 0.993489179d0) then
            nsta = 160
           ELSEIF(br .le. 0.995052633d0) then
            nsta = 162
           ELSEIF(br .le. 0.995502393d0) then
            nsta = 165
           ELSEIF(br .le. 0.998725678d0) then
            nsta = 168
           ELSE                                           
            nsta = 170
           ENDIF


         Case(20041 ) 
           IF(br .le. 0.443396691d0) then
            nsta = 1
           ELSEIF(br .le. 0.451269621d0) then
            nsta = 2
           ELSEIF(br .le. 0.530087597d0) then
            nsta = 3
           ELSEIF(br .le. 0.547814006d0) then
            nsta = 6
           ELSEIF(br .le. 0.551755397d0) then
            nsta = 9
           ELSEIF(br .le. 0.573433051d0) then
            nsta = 13
           ELSEIF(br .le. 0.583286529d0) then
            nsta = 15
           ELSEIF(br .le. 0.60594953 d0) then
            nsta = 16
           ELSEIF(br .le. 0.609890922d0) then
            nsta = 19
           ELSEIF(br .le. 0.622690591d0) then
            nsta = 22
           ELSEIF(br .le. 0.839417856d0) then
            nsta = 25
           ELSEIF(br .le. 0.863076059d0) then
            nsta = 39
           ELSEIF(br .le. 0.941884182d0) then
            nsta = 43
           ELSEIF(br .le. 0.974380955d0) then
            nsta = 44
           ELSEIF(br .le. 0.980293042d0) then
            nsta = 52
           ELSEIF(br .le. 0.983249086d0) then
            nsta = 55
           ELSEIF(br .le. 0.992117217d0) then
            nsta = 68
           ELSEIF(br .le. 0.996058608d0) then
            nsta = 70
           ELSE
            nsta = 71
           ENDIF

         Case(20043 ) 
           IF(br .le. 7.2  d0/ 98.5 d0) then
            nsta = 2
           ELSEIF(br .le. 7.9  d0/ 98.5 d0) then
            nsta = 8
           ELSEIF(br .le. 65.9 d0/ 98.5 d0) then
            nsta = 9
           ELSEIF(br .le. 67   d0/ 98.5 d0) then
            nsta = 12
           ELSEIF(br .le. 67.2 d0/ 98.5 d0) then
            nsta = 15
           ELSEIF(br .le. 72   d0/ 98.5 d0) then
            nsta = 18
           ELSEIF(br .le. 74.6 d0/ 98.5 d0) then
            nsta = 25
           ELSEIF(br .le. 78.2 d0/ 98.5 d0) then
            nsta = 26
           ELSEIF(br .le. 80.8 d0/ 98.5 d0) then
            nsta = 37
           ELSEIF(br .le. 84.2 d0/ 98.5 d0) then
            nsta = 44
           ELSEIF(br .le. 93.4 d0/ 98.5 d0) then
            nsta = 67
           ELSEIF(br .le. 93.6 d0/ 98.5 d0) then
            nsta = 81
           ELSEIF(br .le. 94.1 d0/ 98.5 d0) then
            nsta = 83
           ELSEIF(br .le. 96.4 d0/ 98.5 d0) then
            nsta = 93
           ELSE
            nsta = 98
           ENDIF
      
         Case(20044 ) 
            nsta  = 151

         Case(20045 ) 
           IF(br .le. 0.002071259d0) then
            nsta = 1
           ELSEIF(br .le. 0.108651491d0) then
            nsta = 2
           ELSEIF(br .le. 0.64936966 d0) then
            nsta = 8
           ELSEIF(br .le. 0.745981956d0) then
            nsta = 11
           ELSEIF(br .le. 0.786769831d0) then
            nsta = 19
           ELSEIF(br .le. 0.810688893d0) then
            nsta = 26
           ELSEIF(br .le. 0.913285933d0) then
            nsta = 31
           ELSEIF(br .le. 0.922258071d0) then
            nsta = 41
           ELSEIF(br .le. 0.939206548d0) then
            nsta = 42
           ELSEIF(br .le. 0.942193942d0) then
            nsta = 53
           ELSEIF(br .le. 0.976080938d0) then
            nsta = 56
           ELSEIF(br .le. 0.998008405d0) then
            nsta = 65
           ELSE
            nsta = 71
           ENDIF


         Case(21046 ) 
           IF(br .le. 0.009380476d0) then
            nsta = 0
           ELSEIF(br .le. 0.009206799d0) then
            nsta = 0
           ELSEIF(br .le. 0.011634966d0) then
            nsta = 2
           ELSEIF(br .le. 0.06384055 d0) then
            nsta = 3
           ELSEIF(br .le. 0.071529745d0) then
            nsta = 5
           ELSEIF(br .le. 0.091865641d0) then
            nsta = 6
           ELSEIF(br .le. 0.193241603d0) then
            nsta = 7
           ELSEIF(br .le. 0.220153784d0) then
            nsta = 8
           ELSEIF(br .le. 0.225819506d0) then
            nsta = 10
           ELSEIF(br .le. 0.228955888d0) then
            nsta = 14
           ELSEIF(br .le. 0.254148118d0) then
            nsta = 16
           ELSEIF(br .le. 0.258599757d0) then
            nsta = 18
           ELSEIF(br .le. 0.260623229d0) then
            nsta = 21
           ELSEIF(br .le. 0.26315257 d0) then
            nsta = 23
           ELSEIF(br .le. 0.269020639d0) then
            nsta = 26
           ELSEIF(br .le. 0.290165925d0) then
            nsta = 27
           ELSEIF(br .le. 0.293707001d0) then
            nsta = 30
           ELSEIF(br .le. 0.295022258d0) then
            nsta = 32
           ELSEIF(br .le. 0.297450425d0) then
            nsta = 33
           ELSEIF(br .le. 0.305240793d0) then
            nsta = 36
           ELSEIF(br .le. 0.33903278 d0) then
            nsta = 37
           ELSEIF(br .le. 0.362707406d0) then
            nsta = 38
           ELSEIF(br .le. 0.380109268d0) then
            nsta = 39
           ELSEIF(br .le. 0.381525698d0) then
            nsta = 41
           ELSEIF(br .le. 0.382436261d0) then
            nsta = 42
           ELSEIF(br .le. 0.387494941d0) then
            nsta = 43
           ELSEIF(br .le. 0.394071226d0) then
            nsta = 44
           ELSEIF(br .le. 0.419061109d0) then
            nsta = 46
           ELSEIF(br .le. 0.42189397 d0) then
            nsta = 47
           ELSEIF(br .le. 0.42664913 d0) then
            nsta = 48
           ELSEIF(br .le. 0.432112505d0) then
            nsta = 49
           ELSEIF(br .le. 0.438486443d0) then
            nsta = 50
           ELSEIF(br .le. 0.439194658d0) then
            nsta = 51
           ELSEIF(br .le. 0.445872117d0) then
            nsta = 54
           ELSEIF(br .le. 0.473796034d0) then
            nsta = 55
           ELSEIF(br .le. 0.475010117d0) then
            nsta = 56
           ELSEIF(br .le. 0.507588021d0) then
            nsta = 57
           ELSEIF(br .le. 0.51092675 d0) then
            nsta = 58
           ELSEIF(br .le. 0.520335896d0) then
            nsta = 59
           ELSEIF(br .le. 0.521347633d0) then
            nsta = 62
           ELSEIF(br .le. 0.524585188d0) then
            nsta = 64
           ELSEIF(br .le. 0.529745042d0) then
            nsta = 65
           ELSEIF(br .le. 0.554127883d0) then
            nsta = 66
           ELSEIF(br .le. 0.556252529d0) then
            nsta = 67
           ELSEIF(br .le. 0.561007689d0) then
            nsta = 68
           ELSEIF(br .le. 0.563941724d0) then
            nsta = 69
           ELSEIF(br .le. 0.579117766d0) then
            nsta = 70
           ELSEIF(br .le. 0.585289357d0) then
            nsta = 71
           ELSEIF(br .le. 0.589235127d0) then
            nsta = 74
           ELSEIF(br .le. 0.591056253d0) then
            nsta = 76
           ELSEIF(br .le. 0.596620801d0) then
            nsta = 77
           ELSEIF(br .le. 0.619283691d0) then
            nsta = 78
           ELSEIF(br .le. 0.621610684d0) then
            nsta = 79
           ELSEIF(br .le. 0.623937677d0) then
            nsta = 80
           ELSEIF(br .le. 0.626871712d0) then
            nsta = 81
           ELSEIF(br .le. 0.632537434d0) then
            nsta = 82
           ELSEIF(br .le. 0.63354917 d0) then
            nsta = 83
           ELSEIF(br .le. 0.644475921d0) then
            nsta = 84
           ELSEIF(br .le. 0.647207608d0) then
            nsta = 85
           ELSEIF(br .le. 0.648624039d0) then
            nsta = 86
           ELSEIF(br .le. 0.651861594d0) then
            nsta = 87
           ELSEIF(br .le. 0.659550789d0) then
            nsta = 88
           ELSEIF(br .le. 0.669465803d0) then
            nsta = 90
           ELSEIF(br .le. 0.676143262d0) then
            nsta = 91
           ELSEIF(br .le. 0.679785512d0) then
            nsta = 92
           ELSEIF(br .le. 0.683528936d0) then
            nsta = 93
           ELSEIF(br .le. 0.68757588 d0) then
            nsta = 94
           ELSEIF(br .le. 0.691420478d0) then
            nsta = 96
           ELSEIF(br .le. 0.699008499d0) then
            nsta = 97
           ELSEIF(br .le. 0.70690004 d0) then
            nsta = 98
           ELSEIF(br .le. 0.711149332d0) then
            nsta = 100
           ELSEIF(br .le. 0.715095103d0) then
            nsta = 101
           ELSEIF(br .le. 0.720356131d0) then
            nsta = 102
           ELSEIF(br .le. 0.730574666d0) then
            nsta = 103
           ELSEIF(br .le. 0.733306354d0) then
            nsta = 104
           ELSEIF(br .le. 0.740084986d0) then
            nsta = 105
           ELSEIF(br .le. 0.759105625d0) then
            nsta = 106
           ELSEIF(br .le. 0.765479563d0) then
            nsta = 108
           ELSEIF(br .le. 0.773775799d0) then
            nsta = 109
           ELSEIF(br .le. 0.776709834d0) then
            nsta = 110
           ELSEIF(br .le. 0.778328612d0) then
            nsta = 111
           ELSEIF(br .le. 0.78358964 d0) then
            nsta = 112
           ELSEIF(br .le. 0.785916633d0) then
            nsta = 113
           ELSEIF(br .le. 0.791582355d0) then
            nsta = 114
           ELSEIF(br .le. 0.794212869d0) then
            nsta = 115
           ELSEIF(br .le. 0.795730473d0) then
            nsta = 116
           ELSEIF(br .le. 0.804836099d0) then
            nsta = 117
           ELSEIF(br .le. 0.831647106d0) then
            nsta = 118
           ELSEIF(br .le. 0.833569405d0) then
            nsta = 119
           ELSEIF(br .le. 0.838526912d0) then
            nsta = 120
           ELSEIF(br .le. 0.844597329d0) then
            nsta = 121
           ELSEIF(br .le. 0.850870093d0) then
            nsta = 122
           ELSEIF(br .le. 0.855422906d0) then
            nsta = 123
           ELSEIF(br .le. 0.85694051 d0) then
            nsta = 124
           ELSEIF(br .le. 0.866248482d0) then
            nsta = 128
           ELSEIF(br .le. 0.868676649d0) then
            nsta = 129
           ELSEIF(br .le. 0.874140024d0) then
            nsta = 130
           ELSEIF(br .le. 0.877175233d0) then
            nsta = 131
           ELSEIF(br .le. 0.880716309d0) then
            nsta = 132
           ELSEIF(br .le. 0.88577499 d0) then
            nsta = 133
           ELSEIF(br .le. 0.892755969d0) then
            nsta = 138
           ELSEIF(br .le. 0.896802914d0) then
            nsta = 139
           ELSEIF(br .le. 0.899736949d0) then
            nsta = 140
           ELSEIF(br .le. 0.904997977d0) then
            nsta = 144
           ELSEIF(br .le. 0.908134359d0) then
            nsta = 147
           ELSEIF(br .le. 0.910562525d0) then
            nsta = 148
           ELSEIF(br .le. 0.914204775d0) then
            nsta = 149
           ELSEIF(br .le. 0.920072845d0) then
            nsta = 151
           ELSEIF(br .le. 0.923310401d0) then
            nsta = 152
           ELSEIF(br .le. 0.930696074d0) then
            nsta = 155
           ELSEIF(br .le. 0.935248887d0) then
            nsta = 156
           ELSEIF(br .le. 0.940408741d0) then
            nsta = 159
           ELSEIF(br .le. 0.949919061d0) then
            nsta = 160
           ELSEIF(br .le. 0.952751922d0) then
            nsta = 161
           ELSEIF(br .le. 0.957203561d0) then
            nsta = 163
           ELSEIF(br .le. 0.960542291d0) then
            nsta = 164
           ELSEIF(br .le. 0.96388102 d0) then
            nsta = 166
           ELSEIF(br .le. 0.965499798d0) then
            nsta = 170
           ELSEIF(br .le. 0.971772562d0) then
            nsta = 171
           ELSEIF(br .le. 0.977741805d0) then
            nsta = 174
           ELSEIF(br .le. 0.98067584 d0) then
            nsta = 178
           ELSEIF(br .le. 0.984014569d0) then
            nsta = 180
           ELSEIF(br .le. 0.993019021d0) then
            nsta = 190
           ELSEIF(br .le. 0.997268312d0) then
            nsta = 192
           ELSE
            nsta = 202
           ENDIF


         Case(22047 ) 
           IF(br .le. 0.223043876d0) then
            nsta = 5
           ELSEIF(br .le. 0.504753856d0) then
            nsta = 7     
           ELSEIF(br .le. 0.619761955d0) then
            nsta = 20     
           ELSEIF(br .le. 0.720684555d0) then
            nsta = 30    
           ELSEIF(br .le. 0.741812804d0) then
            nsta = 40
           ELSEIF(br .le. 0.807521656d0) then
            nsta = 52
           ELSEIF(br .le. 0.89203465 d0) then
            nsta = 57
           ELSE
            nsta = 65
           ENDIF


         Case(22048 ) 
           IF(br .le. 0.023183297d0) then
            nsta = 1
           ELSEIF(br .le. 0.02720987 d0) then
            nsta = 2
           ELSEIF(br .le. 0.115062364d0) then
            nsta = 3
           ELSEIF(br .le. 0.119997289d0) then
            nsta = 7
           ELSEIF(br .le. 0.136483189d0) then
            nsta = 8
           ELSEIF(br .le. 0.140686009d0) then
            nsta = 10
           ELSEIF(br .le. 0.146691974d0) then
            nsta = 11
           ELSEIF(br .le. 0.225596529d0) then
            nsta = 13
           ELSEIF(br .le. 0.22720987 d0) then
            nsta = 15
           ELSEIF(br .le. 0.228294469d0) then
            nsta = 17
           ELSEIF(br .le. 0.261917028d0) then
            nsta = 18
           ELSEIF(br .le. 0.278592733d0) then
            nsta = 21
           ELSEIF(br .le. 0.323603579d0) then
            nsta = 22
           ELSEIF(br .le. 0.335100325d0) then
            nsta = 27
           ELSEIF(br .le. 0.336456074d0) then
            nsta = 28
           ELSEIF(br .le. 0.38539859 d0) then
            nsta = 34
           ELSEIF(br .le. 0.616960412d0) then
            nsta = 38
           ELSEIF(br .le. 0.623603579d0) then
            nsta = 44
           ELSEIF(br .le. 0.642041757d0) then
            nsta = 46
           ELSEIF(br .le. 0.646827549d0) then
            nsta = 48
           ELSEIF(br .le. 0.670824295d0) then
            nsta = 49
           ELSEIF(br .le. 0.681941432d0) then
            nsta = 50
           ELSEIF(br .le. 0.692245119d0) then
            nsta = 54
           ELSEIF(br .le. 0.696922451d0) then
            nsta = 56
           ELSEIF(br .le. 0.699742408d0) then
            nsta = 57
           ELSEIF(br .le. 0.755056941d0) then
            nsta = 63
           ELSEIF(br .le. 0.767800976d0) then
            nsta = 65
           ELSEIF(br .le. 0.777020065d0) then
            nsta = 75
           ELSEIF(br .le. 0.781954989d0) then
            nsta = 79
           ELSEIF(br .le. 0.794292299d0) then
            nsta = 88
           ELSEIF(br .le. 0.828321584d0) then
            nsta = 90
           ELSEIF(br .le. 0.850013557d0) then
            nsta = 101
           ELSEIF(br .le. 0.899905098d0) then
            nsta = 113
           ELSEIF(br .le. 0.902169197d0) then
            nsta = 115
           ELSEIF(br .le. 0.903660521d0) then
            nsta = 133
           ELSEIF(br .le. 0.913557484d0) then
            nsta = 137
           ELSEIF(br .le. 0.917218004d0) then
            nsta = 143
           ELSEIF(br .le. 0.923671367d0) then
            nsta = 146
           ELSEIF(br .le. 0.933256508d0) then
            nsta = 152
           ELSEIF(br .le. 0.941906182d0) then
            nsta = 161
           ELSEIF(br .le. 0.947085141d0) then
            nsta = 167
           ELSEIF(br .le. 0.949376356d0) then
            nsta = 173
           ELSEIF(br .le. 0.95210141 d0) then
            nsta = 179
           ELSEIF(br .le. 0.954948482d0) then
            nsta = 183
           ELSEIF(br .le. 0.958500542d0) then
            nsta = 190
           ELSEIF(br .le. 0.959802061d0) then
            nsta = 195
           ELSEIF(br .le. 0.961239154d0) then
            nsta = 199
           ELSEIF(br .le. 0.965170824d0) then
            nsta = 207
           ELSEIF(br .le. 0.983337852d0) then
            nsta = 231
           ELSEIF(br .le. 0.989845445d0) then
            nsta = 244
           ELSEIF(br .le. 0.992136659d0) then
            nsta = 247
           ELSE                                           
            nsta = 250
           ENDIF


         Case(22049 ) 
           IF(br .le. 0.467605918d0) then
            nsta = 1
           ELSEIF(br .le. 0.519012271d0) then
            nsta = 3
           ELSEIF(br .le. 0.827046407d0) then
            nsta = 6
           ELSEIF(br .le. 0.831368984d0) then
            nsta = 10
           ELSEIF(br .le. 0.833449477d0) then
            nsta = 14
           ELSEIF(br .le. 0.862636974d0) then
            nsta = 20
           ELSEIF(br .le. 0.908185628d0) then
            nsta = 21
           ELSEIF(br .le. 0.917840731d0) then
            nsta = 23
           ELSEIF(br .le. 0.92349644 d0) then
            nsta = 25
           ELSEIF(br .le. 0.927081755d0) then
            nsta = 34
           ELSEIF(br .le. 0.928455285d0) then
            nsta = 36
           ELSEIF(br .le. 0.941473514d0) then
            nsta = 44
           ELSEIF(br .le. 0.94738171 d0) then
            nsta = 49
           ELSEIF(br .le. 0.954289754d0) then
            nsta = 55
           ELSEIF(br .le. 0.970145937d0) then
            nsta = 57
           ELSEIF(br .le. 0.970691309d0) then
            nsta = 60
           ELSEIF(br .le. 0.9747614  d0) then
            nsta = 63
           ELSEIF(br .le. 0.995667323d0) then
            nsta = 65
           ELSEIF(br .le. 0.997535727d0) then
            nsta = 67
           ELSEIF(br .le. 0.99902035 d0) then
            nsta = 75
           ELSE
            nsta = 83
           ENDIF


         Case(23051 ) 
           IF(br .le. 4.8   d0/74.875d0) then
            nsta = 6
           ELSEIF(br .le. 7.48  d0/74.875d0) then
            nsta = 7
           ELSEIF(br .le. 9.3   d0/74.875d0) then
            nsta = 13
           ELSEIF(br .le. 11.29 d0/74.875d0) then
            nsta = 24
           ELSEIF(br .le. 21.1  d0/74.875d0) then
            nsta = 27
           ELSEIF(br .le. 21.163d0/74.875d0) then
            nsta = 28
           ELSEIF(br .le. 21.4  d0/74.875d0) then
            nsta = 32
           ELSEIF(br .le. 28.7  d0/74.875d0) then
            nsta = 37
           ELSEIF(br .le. 29.055d0/74.875d0) then
            nsta = 50
           ELSEIF(br .le. 34.355d0/74.875d0) then
            nsta = 53
           ELSEIF(br .le. 37.725d0/74.875d0) then
            nsta = 55
           ELSEIF(br .le. 40.825d0/74.875d0) then
            nsta = 56
           ELSEIF(br .le. 41.566d0/74.875d0) then
            nsta = 58
           ELSEIF(br .le. 41.659d0/74.875d0) then
            nsta = 59
           ELSEIF(br .le. 44.749d0/74.875d0) then
            nsta = 63
           ELSEIF(br .le. 46.669d0/74.875d0) then
            nsta = 66
           ELSEIF(br .le. 49.489d0/74.875d0) then
            nsta = 72
           ELSEIF(br .le. 50.344d0/74.875d0) then
            nsta = 76
           ELSEIF(br .le. 51.29 d0/74.875d0) then
            nsta = 81
           ELSEIF(br .le. 51.622d0/74.875d0) then
            nsta = 82
           ELSEIF(br .le. 52.338d0/74.875d0) then
            nsta = 101
           ELSEIF(br .le. 52.481d0/74.875d0) then
            nsta = 103
           ELSEIF(br .le. 52.796d0/74.875d0) then
            nsta = 110
           ELSEIF(br .le. 53.086d0/74.875d0) then
            nsta = 111
           ELSEIF(br .le. 54.086d0/74.875d0) then
            nsta = 117
           ELSEIF(br .le. 56.356d0/74.875d0) then
            nsta = 123
           ELSEIF(br .le. 57.286d0/74.875d0) then
            nsta = 124
           ELSEIF(br .le. 58.826d0/74.875d0) then
            nsta = 126
           ELSEIF(br .le. 59.145d0/74.875d0) then
            nsta = 129
           ELSEIF(br .le. 59.59 d0/74.875d0) then
            nsta = 131
           ELSEIF(br .le. 61.32 d0/74.875d0) then
            nsta = 132
           ELSEIF(br .le. 62.076d0/74.875d0) then
            nsta = 138
           ELSEIF(br .le. 63.766d0/74.875d0) then
            nsta = 153
           ELSEIF(br .le. 64.35 d0/74.875d0) then
            nsta = 158
           ELSEIF(br .le. 64.911d0/74.875d0) then
            nsta = 159
           ELSEIF(br .le. 65.889d0/74.875d0) then
            nsta = 161
           ELSEIF(br .le. 66.583d0/74.875d0) then
            nsta = 163
           ELSEIF(br .le. 66.919d0/74.875d0) then
            nsta = 164
           ELSEIF(br .le. 67.398d0/74.875d0) then
            nsta = 165
           ELSEIF(br .le. 67.585d0/74.875d0) then
            nsta = 168
           ELSEIF(br .le. 67.998d0/74.875d0) then
            nsta = 170
           ELSEIF(br .le. 68.629d0/74.875d0) then
            nsta = 173
           ELSEIF(br .le. 68.876d0/74.875d0) then
            nsta = 175
           ELSEIF(br .le. 69.806d0/74.875d0) then
            nsta = 176
           ELSEIF(br .le. 70.189d0/74.875d0) then
            nsta = 185
           ELSEIF(br .le. 70.982d0/74.875d0) then
            nsta = 186
           ELSEIF(br .le. 71.196d0/74.875d0) then
            nsta = 187
           ELSEIF(br .le. 71.712d0/74.875d0) then
            nsta = 188
           ELSEIF(br .le. 72.088d0/74.875d0) then
            nsta = 192
           ELSEIF(br .le. 72.617d0/74.875d0) then
            nsta = 194
           ELSEIF(br .le. 72.792d0/74.875d0) then
            nsta = 196
           ELSEIF(br .le. 73.2  d0/74.875d0) then
            nsta = 197
           ELSEIF(br .le. 73.501d0/74.875d0) then
            nsta = 198
           ELSEIF(br .le. 73.777d0/74.875d0) then
            nsta = 201
           ELSEIF(br .le. 73.984d0/74.875d0) then
            nsta = 202
           ELSEIF(br .le. 74.146d0/74.875d0) then
            nsta = 206
           ELSEIF(br .le. 74.458d0/74.875d0) then
            nsta = 209
           ELSEIF(br .le. 74.628d0/74.875d0) then
            nsta = 211
           ELSE
            nsta = 212
           ENDIF

         Case(23052 ) 
           IF(br .le. 5.02 d0/96.38d0) then
            nsta = 0
           ELSEIF(br .le. 7.09 d0/96.38d0) then
            nsta = 1
           ELSEIF(br .le. 8.35 d0/96.38d0) then
            nsta = 2
           ELSEIF(br .le. 21.05d0/96.38d0) then
            nsta = 4
           ELSEIF(br .le. 31.45d0/96.38d0) then
            nsta = 5
           ELSEIF(br .le. 47.95d0/96.38d0) then
            nsta = 6
           ELSEIF(br .le. 56.91d0/96.38d0) then
            nsta = 7
           ELSEIF(br .le. 59.36d0/96.38d0) then
            nsta = 10
           ELSEIF(br .le. 66.87d0/96.38d0) then
            nsta = 12
           ELSEIF(br .le. 66.99d0/96.38d0) then
            nsta = 13
           ELSEIF(br .le. 67.42d0/96.38d0) then
            nsta = 15
           ELSEIF(br .le. 67.97d0/96.38d0) then
            nsta = 16
           ELSEIF(br .le. 68.02d0/96.38d0) then
            nsta = 17
           ELSEIF(br .le. 75.93d0/96.38d0) then
            nsta = 18
           ELSEIF(br .le. 80.74d0/96.38d0) then
            nsta = 20
           ELSEIF(br .le. 84.6 d0/96.38d0) then
            nsta = 22
           ELSEIF(br .le. 85.31d0/96.38d0) then
            nsta = 23
           ELSEIF(br .le. 86.72d0/96.38d0) then
            nsta = 26
           ELSEIF(br .le. 87.05d0/96.38d0) then
            nsta = 28
           ELSEIF(br .le. 87.13d0/96.38d0) then
            nsta = 30
           ELSEIF(br .le. 87.35d0/96.38d0) then
            nsta = 33
           ELSEIF(br .le. 87.51d0/96.38d0) then
            nsta = 34
           ELSEIF(br .le. 87.86d0/96.38d0) then
            nsta = 35
           ELSEIF(br .le. 88.99d0/96.38d0) then
            nsta = 36
           ELSEIF(br .le. 89.17d0/96.38d0) then
            nsta = 39
           ELSEIF(br .le. 89.31d0/96.38d0) then
            nsta = 40
           ELSEIF(br .le. 89.49d0/96.38d0) then
            nsta = 41
           ELSEIF(br .le. 89.62d0/96.38d0) then
            nsta = 43
           ELSEIF(br .le. 91.5 d0/96.38d0) then
            nsta = 44
           ELSEIF(br .le. 91.57d0/96.38d0) then
            nsta = 45
           ELSEIF(br .le. 91.77d0/96.38d0) then
            nsta = 48
           ELSEIF(br .le. 91.97d0/96.38d0) then
            nsta = 49
           ELSEIF(br .le. 92.16d0/96.38d0) then
            nsta = 50
           ELSEIF(br .le. 92.46d0/96.38d0) then
            nsta = 52
           ELSEIF(br .le. 92.77d0/96.38d0) then
            nsta = 53
           ELSEIF(br .le. 92.96d0/96.38d0) then
            nsta = 56
           ELSEIF(br .le. 93.72d0/96.38d0) then
            nsta = 57
           ELSEIF(br .le. 94.26d0/96.38d0) then
            nsta = 58
           ELSEIF(br .le. 94.95d0/96.38d0) then
            nsta = 59
           ELSEIF(br .le. 95.11d0/96.38d0) then
            nsta = 64
           ELSEIF(br .le. 95.38d0/96.38d0) then
            nsta = 66
           ELSEIF(br .le. 95.52d0/96.38d0) then
            nsta = 67
           ELSEIF(br .le. 95.72d0/96.38d0) then
            nsta = 69
           ELSEIF(br .le. 95.89d0/96.38d0) then
            nsta = 71
           ELSEIF(br .le. 96.00d0/96.38d0) then
            nsta = 72
           ELSEIF(br .le. 96.07d0/96.38d0) then
            nsta = 75
           ELSEIF(br .le. 96.12d0/96.38d0) then
            nsta = 77
           ELSE
            nsta = 82
           ENDIF


         Case(24051 ) 
           nsta = 272

         Case(24053 ) 
           nsta = 171

         Case(24054 ) 
           IF(br .le. 14.65 d0/96.04d0) then
            nsta = 0
           ELSEIF(br .le. 58.65 d0/96.04d0) then
            nsta = 1
           ELSEIF(br .le. 66.85 d0/96.04d0) then
            nsta = 3
           ELSEIF(br .le. 69.2  d0/96.04d0) then
            nsta = 4
           ELSEIF(br .le. 79.49 d0/96.04d0) then
            nsta = 5
           ELSEIF(br .le. 80.68 d0/96.04d0) then
            nsta = 8
           ELSEIF(br .le. 82.71 d0/96.04d0) then
            nsta = 9
           ELSEIF(br .le. 87.5  d0/96.04d0) then
            nsta = 13
           ELSEIF(br .le. 88.71 d0/96.04d0) then
            nsta = 16
           ELSEIF(br .le. 88.88 d0/96.04d0) then
            nsta = 18
           ELSEIF(br .le. 89.34 d0/96.04d0) then
            nsta = 19
           ELSEIF(br .le. 90.69 d0/96.04d0) then
            nsta = 21
           ELSEIF(br .le. 90.82 d0/96.04d0) then
            nsta = 23
           ELSEIF(br .le. 90.95 d0/96.04d0) then
            nsta = 27
           ELSEIF(br .le. 91.24 d0/96.04d0) then
            nsta = 30
           ELSEIF(br .le. 91.46 d0/96.04d0) then
            nsta = 36
           ELSEIF(br .le. 93.42 d0/96.04d0) then
            nsta = 42
           ELSEIF(br .le. 93.61 d0/96.04d0) then
            nsta = 52
           ELSEIF(br .le. 93.74 d0/96.04d0) then
            nsta = 55
           ELSEIF(br .le. 94.19 d0/96.04d0) then
            nsta = 56
           ELSEIF(br .le. 94.69 d0/96.04d0) then
            nsta = 59
           ELSEIF(br .le. 95.17 d0/96.04d0) then
            nsta = 67
           ELSEIF(br .le. 95.28 d0/96.04d0) then
            nsta = 73
           ELSEIF(br .le. 95.67 d0/96.04d0) then
            nsta = 74
           ELSEIF(br .le. 95.87 d0/96.04d0) then
            nsta = 80
           ELSE
            nsta = 85
           ENDIF
      
         Case(24055 ) 
           IF(br .le. 55.5  d0/89.79d0) then
            nsta = 0
           ELSEIF(br .le. 68.9  d0/89.79d0) then
            nsta = 1
           ELSEIF(br .le. 71.04 d0/89.79d0) then
            nsta = 3
           ELSEIF(br .le. 84.14 d0/89.79d0) then
            nsta = 8
           ELSEIF(br .le. 84.6  d0/89.79d0) then
            nsta = 11
           ELSEIF(br .le. 84.95 d0/89.79d0) then
            nsta = 13
           ELSEIF(br .le. 86.38 d0/89.79d0) then
            nsta = 20
           ELSEIF(br .le. 89.1  d0/89.79d0) then
            nsta = 25
           ELSE
            nsta = 48
           ENDIF
       
         Case(25056 ) 
           IF(br .le. 3.5   d0/98.195d0) then
            nsta = 0
           ELSEIF(br .le. 15.8  d0/98.195d0) then
            nsta = 1
           ELSEIF(br .le. 21.76 d0/98.195d0) then
            nsta = 2
           ELSEIF(br .le. 32.76 d0/98.195d0) then
            nsta = 3
           ELSEIF(br .le. 35.1  d0/98.195d0) then
            nsta = 6
           ELSEIF(br .le. 38.52 d0/98.195d0) then
            nsta = 8
           ELSEIF(br .le. 38.71 d0/98.195d0) then
            nsta = 11
           ELSEIF(br .le. 39.44 d0/98.195d0) then
            nsta = 13
           ELSEIF(br .le. 41.44 d0/98.195d0) then
            nsta = 18
           ELSEIF(br .le. 41.96 d0/98.195d0) then
            nsta = 22
           ELSEIF(br .le. 42.14 d0/98.195d0) then
            nsta = 23
           ELSEIF(br .le. 42.96 d0/98.195d0) then
            nsta = 24
           ELSEIF(br .le. 44.7  d0/98.195d0) then
            nsta = 30
           ELSEIF(br .le. 44.87 d0/98.195d0) then
            nsta = 36
           ELSEIF(br .le. 51.59 d0/98.195d0) then
            nsta = 37
           ELSEIF(br .le. 52.86 d0/98.195d0) then
            nsta = 39
           ELSEIF(br .le. 53.25 d0/98.195d0) then
            nsta = 40
           ELSEIF(br .le. 53.305d0/98.195d0) then
            nsta = 45
           ELSEIF(br .le. 54.545d0/98.195d0) then
            nsta = 46
           ELSEIF(br .le. 55.395d0/98.195d0) then
            nsta = 48
           ELSEIF(br .le. 59.185d0/98.195d0) then
            nsta = 49
           ELSEIF(br .le. 59.485d0/98.195d0) then
            nsta = 52
           ELSEIF(br .le. 62.215d0/98.195d0) then
            nsta = 53
           ELSEIF(br .le. 63.235d0/98.195d0) then
            nsta = 54
           ELSEIF(br .le. 69.355d0/98.195d0) then
            nsta = 55
           ELSEIF(br .le. 69.715d0/98.195d0) then
            nsta = 57
           ELSEIF(br .le. 71.725d0/98.195d0) then
            nsta = 58
           ELSEIF(br .le. 72.085d0/98.195d0) then
            nsta = 59
           ELSEIF(br .le. 72.765d0/98.195d0) then
            nsta = 60
           ELSEIF(br .le. 73.485d0/98.195d0) then
            nsta = 61
           ELSEIF(br .le. 74.735d0/98.195d0) then
            nsta = 64
           ELSEIF(br .le. 74.76 d0/98.195d0) then
            nsta = 66
           ELSEIF(br .le. 77.52 d0/98.195d0) then
            nsta = 67
           ELSEIF(br .le. 78.52 d0/98.195d0) then
            nsta = 69
           ELSEIF(br .le. 78.62 d0/98.195d0) then
            nsta = 73
           ELSEIF(br .le. 79.05 d0/98.195d0) then
            nsta = 74
           ELSEIF(br .le. 80.99 d0/98.195d0) then
            nsta = 75
           ELSEIF(br .le. 81.67 d0/98.195d0) then
            nsta = 76
           ELSEIF(br .le. 81.77 d0/98.195d0) then
            nsta = 78
           ELSEIF(br .le. 82.54 d0/98.195d0) then
            nsta = 79
           ELSEIF(br .le. 82.74 d0/98.195d0) then
            nsta = 80
           ELSEIF(br .le. 82.88 d0/98.195d0) then
            nsta = 81
           ELSEIF(br .le. 83.42 d0/98.195d0) then
            nsta = 82
           ELSEIF(br .le. 83.73 d0/98.195d0) then
            nsta = 83
           ELSEIF(br .le. 83.84 d0/98.195d0) then
            nsta = 85
           ELSEIF(br .le. 84.53 d0/98.195d0) then
            nsta = 86
           ELSEIF(br .le. 84.74 d0/98.195d0) then
            nsta = 87
           ELSEIF(br .le. 85.47 d0/98.195d0) then
            nsta = 88
           ELSEIF(br .le. 85.64 d0/98.195d0) then
            nsta = 90
           ELSEIF(br .le. 85.99 d0/98.195d0) then
            nsta = 93
           ELSEIF(br .le. 86.38 d0/98.195d0) then
            nsta = 94
           ELSEIF(br .le. 86.51 d0/98.195d0) then
            nsta = 95
           ELSEIF(br .le. 86.68 d0/98.195d0) then
            nsta = 97
           ELSEIF(br .le. 86.73 d0/98.195d0) then
            nsta = 99
           ELSEIF(br .le. 87.11 d0/98.195d0) then
            nsta = 100
           ELSEIF(br .le. 87.17 d0/98.195d0) then
            nsta = 102
           ELSEIF(br .le. 87.265d0/98.195d0) then
            nsta = 104
           ELSEIF(br .le. 87.395d0/98.195d0) then
            nsta = 106
           ELSEIF(br .le. 87.995d0/98.195d0) then
            nsta = 107
           ELSEIF(br .le. 88.575d0/98.195d0) then
            nsta = 109
           ELSEIF(br .le. 89.645d0/98.195d0) then
            nsta = 110
           ELSEIF(br .le. 89.965d0/98.195d0) then
            nsta = 112
           ELSEIF(br .le. 90.085d0/98.195d0) then
            nsta = 113
           ELSEIF(br .le. 90.515d0/98.195d0) then
            nsta = 114
           ELSEIF(br .le. 91.125d0/98.195d0) then
            nsta = 119
           ELSEIF(br .le. 91.565d0/98.195d0) then
            nsta = 122
           ELSEIF(br .le. 91.805d0/98.195d0) then
            nsta = 123
           ELSEIF(br .le. 92.075d0/98.195d0) then
            nsta = 124
           ELSEIF(br .le. 92.615d0/98.195d0) then
            nsta = 126
           ELSEIF(br .le. 92.775d0/98.195d0) then
            nsta = 130
           ELSEIF(br .le. 95.945d0/98.195d0) then
            nsta = 131
           ELSEIF(br .le. 96.065d0/98.195d0) then
            nsta = 137
           ELSEIF(br .le. 96.385d0/98.195d0) then
            nsta = 139
           ELSEIF(br .le. 96.555d0/98.195d0) then
            nsta = 143
           ELSEIF(br .le. 96.885d0/98.195d0) then
            nsta = 178
           ELSEIF(br .le. 97.105d0/98.195d0) then
            nsta = 181
           ELSEIF(br .le. 97.505d0/98.195d0) then
            nsta = 184
           ELSEIF(br .le. 97.795d0/98.195d0) then
            nsta = 187
           ELSE
            nsta = 189
           ENDIF


         Case(26055 ) 
           IF(br .le. 65.997 d0/ 97.586d0) then
            nsta = 0
           ELSEIF(br .le. 78.297 d0/ 97.586d0) then
            nsta = 1
           ELSEIF(br .le. 80.297 d0/ 97.586d0) then
            nsta = 7
           ELSEIF(br .le. 82.194 d0/ 97.586d0) then
            nsta = 12
           ELSEIF(br .le. 85.39  d0/ 97.586d0) then
            nsta = 22
           ELSEIF(br .le. 87.689 d0/ 97.586d0) then
            nsta = 31
           ELSEIF(br .le. 90.091 d0/ 97.586d0) then
            nsta = 36
           ELSEIF(br .le. 91.287 d0/ 97.586d0) then
            nsta = 42
           ELSEIF(br .le. 94.184 d0/ 97.586d0) then
            nsta = 53
           ELSEIF(br .le. 96.782 d0/ 97.586d0) then
            nsta = 58
           ELSE
            nsta = 71
           ENDIF


         Case(26057 ) 
           nsta = 177

         Case(26058 ) 
           nsta = 189

         Case(26059 ) 
           IF(br .le. 3.4    d0/92.189d0) then
            nsta = 0
           ELSEIF(br .le. 50.4   d0/92.189d0) then
            nsta = 1
           ELSEIF(br .le. 50.94  d0/92.189d0) then
            nsta = 3
           ELSEIF(br .le. 65.84  d0/92.189d0) then
            nsta = 6
           ELSEIF(br .le. 72.939 d0/92.189d0) then
            nsta = 9
           ELSEIF(br .le. 73.269 d0/92.189d0) then
            nsta = 10
           ELSEIF(br .le. 73.639 d0/92.189d0) then
            nsta = 12
           ELSEIF(br .le. 81.839 d0/92.189d0) then
            nsta = 15
           ELSEIF(br .le. 85.239 d0/92.189d0) then
            nsta = 16
           ELSEIF(br .le. 85.729 d0/92.189d0) then
            nsta = 17
           ELSEIF(br .le. 86.119 d0/92.189d0) then
            nsta = 20
           ELSEIF(br .le. 89.699 d0/92.189d0) then
            nsta = 24
           ELSEIF(br .le. 90.409 d0/92.189d0) then
            nsta = 26
           ELSEIF(br .le. 90.539 d0/92.189d0) then
            nsta = 28
           ELSEIF(br .le. 90.869 d0/92.189d0) then
            nsta = 29
           ELSEIF(br .le. 90.909 d0/92.189d0) then
            nsta = 34
           ELSEIF(br .le. 90.969 d0/92.189d0) then
            nsta = 35
           ELSEIF(br .le. 91.029 d0/92.189d0) then
            nsta = 36
           ELSEIF(br .le. 92.089 d0/92.189d0) then
            nsta = 40
           ELSE
            nsta = 43
           ENDIF

         Case(27060 ) 
           IF(br .le. 2.88  d0/90.26d0) then
            nsta = 0
           ELSEIF(br .le. 3.1   d0/90.26d0) then
            nsta = 1
           ELSEIF(br .le. 6.8   d0/90.26d0) then
            nsta = 2
           ELSEIF(br .le. 7.65  d0/90.26d0) then
            nsta = 3
           ELSEIF(br .le. 9.41  d0/90.26d0) then
            nsta = 4
           ELSEIF(br .le. 12.17 d0/90.26d0) then
            nsta = 5
           ELSEIF(br .le. 12.86 d0/90.26d0) then
            nsta = 6
           ELSEIF(br .le. 21.28 d0/90.26d0) then
            nsta = 7
           ELSEIF(br .le. 28.82 d0/90.26d0) then
            nsta = 9
           ELSEIF(br .le. 35.07 d0/90.26d0) then
            nsta = 12
           ELSEIF(br .le. 35.1  d0/90.26d0) then
            nsta = 13
           ELSEIF(br .le. 35.18 d0/90.26d0) then
            nsta = 14
           ELSEIF(br .le. 35.79 d0/90.26d0) then
            nsta = 15
           ELSEIF(br .le. 36.43 d0/90.26d0) then
            nsta = 17
           ELSEIF(br .le. 36.9  d0/90.26d0) then
            nsta = 18
           ELSEIF(br .le. 37.44 d0/90.26d0) then
            nsta = 20
           ELSEIF(br .le. 37.96 d0/90.26d0) then
            nsta = 21
           ELSEIF(br .le. 44.86 d0/90.26d0) then
            nsta = 24
           ELSEIF(br .le. 46.67 d0/90.26d0) then
            nsta = 25
           ELSEIF(br .le. 46.96 d0/90.26d0) then
            nsta = 26
           ELSEIF(br .le. 47.01 d0/90.26d0) then
            nsta = 29
           ELSEIF(br .le. 49.18 d0/90.26d0) then
            nsta = 30
           ELSEIF(br .le. 49.65 d0/90.26d0) then
            nsta = 31
           ELSEIF(br .le. 56.79 d0/90.26d0) then
            nsta = 34
           ELSEIF(br .le. 57.75 d0/90.26d0) then
            nsta = 36
           ELSEIF(br .le. 58.79 d0/90.26d0) then
            nsta = 37
           ELSEIF(br .le. 59.97 d0/90.26d0) then
            nsta = 38
           ELSEIF(br .le. 60.07 d0/90.26d0) then
            nsta = 39
           ELSEIF(br .le. 60.57 d0/90.26d0) then
            nsta = 40
           ELSEIF(br .le. 60.66 d0/90.26d0) then
            nsta = 43
           ELSEIF(br .le. 60.69 d0/90.26d0) then
            nsta = 44
           ELSEIF(br .le. 61.19 d0/90.26d0) then
            nsta = 45
           ELSEIF(br .le. 61.63 d0/90.26d0) then
            nsta = 46
           ELSEIF(br .le. 61.76 d0/90.26d0) then
            nsta = 49
           ELSEIF(br .le. 62.81 d0/90.26d0) then
            nsta = 51
           ELSEIF(br .le. 62.97 d0/90.26d0) then
            nsta = 53
           ELSEIF(br .le. 63.18 d0/90.26d0) then
            nsta = 54
           ELSEIF(br .le. 65.85 d0/90.26d0) then
            nsta = 55
           ELSEIF(br .le. 66.13 d0/90.26d0) then
            nsta = 57
           ELSEIF(br .le. 66.27 d0/90.26d0) then
            nsta = 58
           ELSEIF(br .le. 66.36 d0/90.26d0) then
            nsta = 59
           ELSEIF(br .le. 66.84 d0/90.26d0) then
            nsta = 60
           ELSEIF(br .le. 67.15 d0/90.26d0) then
            nsta = 61
           ELSEIF(br .le. 67.19 d0/90.26d0) then
            nsta = 62
           ELSEIF(br .le. 67.43 d0/90.26d0) then
            nsta = 63
           ELSEIF(br .le. 68.16 d0/90.26d0) then
            nsta = 65
           ELSEIF(br .le. 68.35 d0/90.26d0) then
            nsta = 66
           ELSEIF(br .le. 69.27 d0/90.26d0) then
            nsta = 69
           ELSEIF(br .le. 70.97 d0/90.26d0) then
            nsta = 70
           ELSEIF(br .le. 71.58 d0/90.26d0) then
            nsta = 71
           ELSEIF(br .le. 72.29 d0/90.26d0) then
            nsta = 72
           ELSEIF(br .le. 72.4  d0/90.26d0) then
            nsta = 73
           ELSEIF(br .le. 72.51 d0/90.26d0) then
            nsta = 74
           ELSEIF(br .le. 72.63 d0/90.26d0) then
            nsta = 76
           ELSEIF(br .le. 72.67 d0/90.26d0) then
            nsta = 78
           ELSEIF(br .le. 72.92 d0/90.26d0) then
            nsta = 79
           ELSEIF(br .le. 73.06 d0/90.26d0) then
            nsta = 81
           ELSEIF(br .le. 73.42 d0/90.26d0) then
            nsta = 82
           ELSEIF(br .le. 73.67 d0/90.26d0) then
            nsta = 86
           ELSEIF(br .le. 73.91 d0/90.26d0) then
            nsta = 87
           ELSEIF(br .le. 74.21 d0/90.26d0) then
            nsta = 88
           ELSEIF(br .le. 75.05 d0/90.26d0) then
            nsta = 90
           ELSEIF(br .le. 75.18 d0/90.26d0) then
            nsta = 92
           ELSEIF(br .le. 75.59 d0/90.26d0) then
            nsta = 97
           ELSEIF(br .le. 75.77 d0/90.26d0) then
            nsta = 98
           ELSEIF(br .le. 75.89 d0/90.26d0) then
            nsta = 100
           ELSEIF(br .le. 76.08 d0/90.26d0) then
            nsta = 101
           ELSEIF(br .le. 76.17 d0/90.26d0) then
            nsta = 102
           ELSEIF(br .le. 76.2  d0/90.26d0) then
            nsta = 103
           ELSEIF(br .le. 76.27 d0/90.26d0) then
            nsta = 104
           ELSEIF(br .le. 76.56 d0/90.26d0) then
            nsta = 106
           ELSEIF(br .le. 76.91 d0/90.26d0) then
            nsta = 108
           ELSEIF(br .le. 77.14 d0/90.26d0) then
            nsta = 109
           ELSEIF(br .le. 77.46 d0/90.26d0) then
            nsta = 111
           ELSEIF(br .le. 77.56 d0/90.26d0) then
            nsta = 112
           ELSEIF(br .le. 77.87 d0/90.26d0) then
            nsta = 113
           ELSEIF(br .le. 77.89 d0/90.26d0) then
            nsta = 114
           ELSEIF(br .le. 78.06 d0/90.26d0) then
            nsta = 115
           ELSEIF(br .le. 78.16 d0/90.26d0) then
            nsta = 117
           ELSEIF(br .le. 78.34 d0/90.26d0) then
            nsta = 118
           ELSEIF(br .le. 78.46 d0/90.26d0) then
            nsta = 119
           ELSEIF(br .le. 78.73 d0/90.26d0) then
            nsta = 121
           ELSEIF(br .le. 79.48 d0/90.26d0) then
            nsta = 122
           ELSEIF(br .le. 79.92 d0/90.26d0) then
            nsta = 124
           ELSEIF(br .le. 80.12 d0/90.26d0) then
            nsta = 125
           ELSEIF(br .le. 80.3  d0/90.26d0) then
            nsta = 128
           ELSEIF(br .le. 80.97 d0/90.26d0) then
            nsta = 130
           ELSEIF(br .le. 81.8  d0/90.26d0) then
            nsta = 131
           ELSEIF(br .le. 82    d0/90.26d0) then
            nsta = 133
           ELSEIF(br .le. 82.21 d0/90.26d0) then
            nsta = 134
           ELSEIF(br .le. 83.11 d0/90.26d0) then
            nsta = 137
           ELSEIF(br .le. 83.3  d0/90.26d0) then
            nsta = 138
           ELSEIF(br .le. 83.55 d0/90.26d0) then
            nsta = 139
           ELSEIF(br .le. 83.93 d0/90.26d0) then
            nsta = 142
           ELSEIF(br .le. 84.05 d0/90.26d0) then
            nsta = 145
           ELSEIF(br .le. 85.44 d0/90.26d0) then
            nsta = 147
           ELSEIF(br .le. 85.64 d0/90.26d0) then
            nsta = 149
           ELSEIF(br .le. 85.86 d0/90.26d0) then
            nsta = 150
           ELSEIF(br .le. 86.05 d0/90.26d0) then
            nsta = 152
           ELSEIF(br .le. 86.21 d0/90.26d0) then
            nsta = 153
           ELSEIF(br .le. 86.37 d0/90.26d0) then
            nsta = 154
           ELSEIF(br .le. 86.53 d0/90.26d0) then
            nsta = 158
           ELSEIF(br .le. 86.66 d0/90.26d0) then
            nsta = 164
           ELSEIF(br .le. 87.36 d0/90.26d0) then
            nsta = 165
           ELSEIF(br .le. 87.74 d0/90.26d0) then
            nsta = 167
           ELSEIF(br .le. 87.83 d0/90.26d0) then
            nsta = 169
           ELSEIF(br .le. 87.98 d0/90.26d0) then
            nsta = 171
           ELSEIF(br .le. 88.17 d0/90.26d0) then
            nsta = 172
           ELSEIF(br .le. 88.24 d0/90.26d0) then
            nsta = 173
           ELSEIF(br .le. 88.37 d0/90.26d0) then
            nsta = 176
           ELSEIF(br .le. 88.7  d0/90.26d0) then
            nsta = 177
           ELSEIF(br .le. 88.97 d0/90.26d0) then
            nsta = 188
           ELSEIF(br .le. 89.17 d0/90.26d0) then
            nsta = 190
           ELSEIF(br .le. 89.3  d0/90.26d0) then
            nsta = 194
           ELSEIF(br .le. 89.61 d0/90.26d0) then
            nsta = 201
           ELSEIF(br .le. 89.82 d0/90.26d0) then
            nsta = 205
           ELSEIF(br .le. 90.06 d0/90.26d0) then
            nsta = 206
           ELSE
            nsta = 210
           ENDIF                                           


         Case(28059 ) 
           IF(br .le. 50.1   d0/97.65d0) then
            nsta = 0
           ELSEIF(br .le. 73.9   d0/97.65d0) then
            nsta = 2
           ELSEIF(br .le. 78.2   d0/97.65d0) then
            nsta = 3
           ELSEIF(br .le. 79.45  d0/97.65d0) then
            nsta = 5
           ELSEIF(br .le. 79.686 d0/97.65d0) then
            nsta = 9
           ELSEIF(br .le. 82.296 d0/97.65d0) then
            nsta = 15
           ELSEIF(br .le. 84.646 d0/97.65d0) then
            nsta = 27
           ELSEIF(br .le. 85.527 d0/97.65d0) then
            nsta = 28
           ELSEIF(br .le. 89.127 d0/97.65d0) then
            nsta = 35
           ELSEIF(br .le. 89.253 d0/97.65d0) then
            nsta = 42
           ELSEIF(br .le. 89.319 d0/97.65d0) then
            nsta = 48
           ELSEIF(br .le. 89.398 d0/97.65d0) then
            nsta = 51
           ELSEIF(br .le. 90.029 d0/97.65d0) then
            nsta = 55
           ELSEIF(br .le. 91.769 d0/97.65d0) then
            nsta = 58
           ELSEIF(br .le. 92.106 d0/97.65d0) then
            nsta = 60
           ELSEIF(br .le. 92.237 d0/97.65d0) then
            nsta = 65
           ELSEIF(br .le. 92.263 d0/97.65d0) then
            nsta = 66
           ELSEIF(br .le. 92.371 d0/97.65d0) then
            nsta = 67
           ELSEIF(br .le. 92.624 d0/97.65d0) then
            nsta = 71
           ELSEIF(br .le. 94.084 d0/97.65d0) then
            nsta = 77
           ELSEIF(br .le. 94.201 d0/97.65d0) then
            nsta = 84
           ELSEIF(br .le. 94.321 d0/97.65d0) then
            nsta = 87
           ELSEIF(br .le. 94.747 d0/97.65d0) then
            nsta = 103
           ELSEIF(br .le. 95.244 d0/97.65d0) then
            nsta = 115
           ELSEIF(br .le. 95.636 d0/97.65d0) then
            nsta = 120
           ELSEIF(br .le. 95.813 d0/97.65d0) then
            nsta = 134
           ELSEIF(br .le. 96.041 d0/97.65d0) then
            nsta = 142
           ELSEIF(br .le. 96.301 d0/97.65d0) then
            nsta = 143
           ELSEIF(br .le. 96.474 d0/97.65d0) then
            nsta = 146
           ELSEIF(br .le. 96.548 d0/97.65d0) then
            nsta = 148
           ELSEIF(br .le. 96.703 d0/97.65d0) then
            nsta = 159
           ELSEIF(br .le. 96.815 d0/97.65d0) then
            nsta = 162
           ELSEIF(br .le. 97.135 d0/97.65d0) then
            nsta = 164
           ELSEIF(br .le. 97.254 d0/97.65d0) then
            nsta = 168
           ELSEIF(br .le. 97.381 d0/97.65d0) then
            nsta = 169
           ELSEIF(br .le. 97.602 d0/97.65d0) then
            nsta = 189
           ELSE
             nsta = 221
           ENDIF

         Case(28061 ) 
           nsta = 258

         Case(28062 ) 
           nsta = 145

         Case(28063 ) 
           nsta = 99

         Case(28065 ) 
           IF(br .le. 0.19   d0/100.245d0) then
            nsta = 1
           ELSEIF(br .le. 0.45   d0/100.245d0) then
            nsta = 2
           ELSEIF(br .le. 1.28   d0/100.245d0) then
            nsta = 3
           ELSEIF(br .le. 1.58   d0/100.245d0) then
            nsta = 7
           ELSEIF(br .le. 2.71   d0/100.245d0) then
            nsta = 11
           ELSEIF(br .le. 2.83   d0/100.245d0) then
            nsta = 13
           ELSEIF(br .le. 6.641  d0/100.245d0) then
            nsta = 17
           ELSEIF(br .le. 15.541 d0/100.245d0) then
            nsta = 24
           ELSEIF(br .le. 33.242 d0/100.245d0) then
            nsta = 29
           ELSE
            nsta = 34
           ENDIF
      
         Case(29064 ) 
           IF(br .le. 33.1   d0/98.869d0) then
            nsta =0
           ELSEIF(br .le. 34.68  d0/98.869d0) then
            nsta =1
           ELSEIF(br .le. 50.88  d0/98.869d0) then
            nsta =2
           ELSEIF(br .le. 52.62  d0/98.869d0) then
            nsta =3
           ELSEIF(br .le. 52.704 d0/98.869d0) then
            nsta =4
           ELSEIF(br .le. 61.664 d0/98.869d0) then
            nsta =6
           ELSEIF(br .le. 65.814 d0/98.869d0) then
            nsta =7
           ELSEIF(br .le. 68.374 d0/98.869d0) then
            nsta =8
           ELSEIF(br .le. 68.677 d0/98.869d0) then
            nsta =9
           ELSEIF(br .le. 69.065 d0/98.869d0) then
            nsta =10
           ELSEIF(br .le. 72.575 d0/98.869d0) then
            nsta =12
           ELSEIF(br .le. 74.565 d0/98.869d0) then
            nsta =13
           ELSEIF(br .le. 74.664 d0/98.869d0) then
            nsta =15
           ELSEIF(br .le. 75.794 d0/98.869d0) then
            nsta =18
           ELSEIF(br .le. 76.425 d0/98.869d0) then
            nsta =19
           ELSEIF(br .le. 76.597 d0/98.869d0) then
            nsta =22
           ELSEIF(br .le. 76.692 d0/98.869d0) then
            nsta =24
           ELSEIF(br .le. 78.084 d0/98.869d0) then
            nsta =25
           ELSEIF(br .le. 78.19  d0/98.869d0) then
            nsta =26
           ELSEIF(br .le. 78.551 d0/98.869d0) then
            nsta =28
           ELSEIF(br .le. 78.718 d0/98.869d0) then
            nsta =29
           ELSEIF(br .le. 78.827 d0/98.869d0) then
            nsta =33
           ELSEIF(br .le. 78.993 d0/98.869d0) then
            nsta =38
           ELSEIF(br .le. 79.19  d0/98.869d0) then
            nsta =40
           ELSEIF(br .le. 79.797 d0/98.869d0) then
            nsta =41
           ELSEIF(br .le. 80.026 d0/98.869d0) then
            nsta =43
           ELSEIF(br .le. 81.616 d0/98.869d0) then
            nsta =44
           ELSEIF(br .le. 81.76  d0/98.869d0) then
            nsta =52
           ELSEIF(br .le. 82.269 d0/98.869d0) then
            nsta =62
           ELSEIF(br .le. 82.679 d0/98.869d0) then
            nsta =71
           ELSEIF(br .le. 83.173 d0/98.869d0) then
            nsta =72
           ELSEIF(br .le. 83.27  d0/98.869d0) then
            nsta =74
           ELSEIF(br .le. 83.544 d0/98.869d0) then
            nsta =78
           ELSEIF(br .le. 83.668 d0/98.869d0) then
            nsta =82
           ELSEIF(br .le. 83.826 d0/98.869d0) then
            nsta =86
           ELSEIF(br .le. 85.676 d0/98.869d0) then
            nsta =89
           ELSEIF(br .le. 86.075 d0/98.869d0) then
            nsta =90
           ELSEIF(br .le. 86.086 d0/98.869d0) then
            nsta =92
           ELSEIF(br .le. 86.141 d0/98.869d0) then
            nsta =95
           ELSEIF(br .le. 86.391 d0/98.869d0) then
            nsta =98
           ELSEIF(br .le. 86.471 d0/98.869d0) then
            nsta =99
           ELSEIF(br .le. 87.502 d0/98.869d0) then
            nsta =100
           ELSEIF(br .le. 88.17  d0/98.869d0) then
            nsta =106
           ELSEIF(br .le. 88.536 d0/98.869d0) then
            nsta =107
           ELSEIF(br .le. 88.805 d0/98.869d0) then
            nsta =108
           ELSEIF(br .le. 88.984 d0/98.869d0) then
            nsta =109
           ELSEIF(br .le. 89.312 d0/98.869d0) then
            nsta =111
           ELSEIF(br .le. 89.495 d0/98.869d0) then
            nsta =114
           ELSEIF(br .le. 89.685 d0/98.869d0) then
            nsta =115
           ELSEIF(br .le. 89.887 d0/98.869d0) then
            nsta =116
           ELSEIF(br .le. 90.038 d0/98.869d0) then
            nsta =119
           ELSEIF(br .le. 90.077 d0/98.869d0) then
            nsta =120
           ELSEIF(br .le. 90.107 d0/98.869d0) then
            nsta =124
           ELSEIF(br .le. 90.226 d0/98.869d0) then
            nsta =125
           ELSEIF(br .le. 90.333 d0/98.869d0) then
            nsta =126
           ELSEIF(br .le. 90.418 d0/98.869d0) then
            nsta =131
           ELSEIF(br .le. 91.188 d0/98.869d0) then
            nsta =133
           ELSEIF(br .le. 91.731 d0/98.869d0) then
            nsta =135
           ELSEIF(br .le. 91.831 d0/98.869d0) then
            nsta =136
           ELSEIF(br .le. 92.144 d0/98.869d0) then
            nsta =137
           ELSEIF(br .le. 92.631 d0/98.869d0) then
            nsta =138
           ELSEIF(br .le. 93.109 d0/98.869d0) then
            nsta =139
           ELSEIF(br .le. 93.281 d0/98.869d0) then
            nsta =140
           ELSEIF(br .le. 93.476 d0/98.869d0) then
            nsta =141
           ELSEIF(br .le. 93.616 d0/98.869d0) then
            nsta =142
           ELSEIF(br .le. 93.92  d0/98.869d0) then
            nsta =143
           ELSEIF(br .le. 94.053 d0/98.869d0) then
            nsta =144
           ELSEIF(br .le. 95.403 d0/98.869d0) then
            nsta =145
           ELSEIF(br .le. 95.683 d0/98.869d0) then
            nsta =146
           ELSEIF(br .le. 96.02  d0/98.869d0) then
            nsta =147
           ELSEIF(br .le. 96.274 d0/98.869d0) then
            nsta =149
           ELSEIF(br .le. 96.655 d0/98.869d0) then
            nsta =151
           ELSEIF(br .le. 96.909 d0/98.869d0) then
            nsta =154
           ELSEIF(br .le. 97.055 d0/98.869d0) then
            nsta =159
           ELSEIF(br .le. 97.545 d0/98.869d0) then
            nsta =160
           ELSEIF(br .le. 97.747 d0/98.869d0) then
            nsta =161
           ELSEIF(br .le. 97.809 d0/98.869d0) then
            nsta =163
           ELSEIF(br .le. 98.149 d0/98.869d0) then
            nsta =165
           ELSEIF(br .le. 98.264 d0/98.869d0) then
            nsta =167
           ELSEIF(br .le. 98.447 d0/98.869d0) then
            nsta =168
           ELSEIF(br .le. 98.591 d0/98.869d0) then
            nsta =169
           ELSE
             nsta = 172
           ENDIF       

         Case(29066 ) 
           IF(br .le. 2.23    d0/92.617d0) then
              nsta = 0
           ELSEIF(br .le. 2.581   d0/92.617d0) then
              nsta = 1
           ELSEIF(br .le. 2.603   d0/92.617d0) then
              nsta = 2
           ELSEIF(br .le. 5.233   d0/92.617d0) then
              nsta = 3
           ELSEIF(br .le. 19.033  d0/92.617d0) then
              nsta = 4
           ELSEIF(br .le. 33.433  d0/92.617d0) then
              nsta = 5
           ELSEIF(br .le. 33.577  d0/92.617d0) then
              nsta = 8
           ELSEIF(br .le. 36.017  d0/92.617d0) then
              nsta = 9
           ELSEIF(br .le. 37.727  d0/92.617d0) then
              nsta = 13
           ELSEIF(br .le. 37.847  d0/92.617d0) then
              nsta = 14
           ELSEIF(br .le. 37.882  d0/92.617d0) then
              nsta = 16
           ELSEIF(br .le. 38.994  d0/92.617d0) then
              nsta = 17
           ELSEIF(br .le. 39.377  d0/92.617d0) then
              nsta = 19
           ELSEIF(br .le. 39.422  d0/92.617d0) then
              nsta = 21
           ELSEIF(br .le. 40.081  d0/92.617d0) then
              nsta = 23
           ELSEIF(br .le. 40.176  d0/92.617d0) then
              nsta = 24
           ELSEIF(br .le. 40.269  d0/92.617d0) then
              nsta = 25
           ELSEIF(br .le. 40.801  d0/92.617d0) then
              nsta = 27
           ELSEIF(br .le. 46.951  d0/92.617d0) then
              nsta = 31
           ELSEIF(br .le. 54.251  d0/92.617d0) then
              nsta = 32
           ELSEIF(br .le. 54.514  d0/92.617d0) then
              nsta = 36
           ELSEIF(br .le. 58.014  d0/92.617d0) then
              nsta = 38
           ELSEIF(br .le. 63.874  d0/92.617d0) then
              nsta = 39
           ELSEIF(br .le. 64.246  d0/92.617d0) then
              nsta = 40
           ELSEIF(br .le. 65.016  d0/92.617d0) then
              nsta = 41
           ELSEIF(br .le. 65.336  d0/92.617d0) then
              nsta = 42
           ELSEIF(br .le. 65.568  d0/92.617d0) then
              nsta = 44
           ELSEIF(br .le. 65.619  d0/92.617d0) then
              nsta = 47
           ELSEIF(br .le. 65.815  d0/92.617d0) then
              nsta = 48
           ELSEIF(br .le. 65.987  d0/92.617d0) then
              nsta = 49
           ELSEIF(br .le. 66.369  d0/92.617d0) then
              nsta = 50
           ELSEIF(br .le. 66.704  d0/92.617d0) then
              nsta = 51
           ELSEIF(br .le. 66.748  d0/92.617d0) then
              nsta = 52
           ELSEIF(br .le. 66.803  d0/92.617d0) then
              nsta = 54
           ELSEIF(br .le. 67.873  d0/92.617d0) then
              nsta = 56
           ELSEIF(br .le. 68.14   d0/92.617d0) then
              nsta = 57
           ELSEIF(br .le. 68.363  d0/92.617d0) then
              nsta = 58
           ELSEIF(br .le. 68.418  d0/92.617d0) then
              nsta = 59
           ELSEIF(br .le. 69.448  d0/92.617d0) then
              nsta = 61
           ELSEIF(br .le. 72.948  d0/92.617d0) then
              nsta = 62
           ELSEIF(br .le. 73.718  d0/92.617d0) then
              nsta = 63
           ELSEIF(br .le. 74.564  d0/92.617d0) then
              nsta = 65
           ELSEIF(br .le. 74.791  d0/92.617d0) then
              nsta = 66
           ELSEIF(br .le. 74.906  d0/92.617d0) then
              nsta = 67
           ELSEIF(br .le. 75.034  d0/92.617d0) then
              nsta = 68
           ELSEIF(br .le. 75.384  d0/92.617d0) then
              nsta = 69
           ELSEIF(br .le. 76.474  d0/92.617d0) then
              nsta = 70
           ELSEIF(br .le. 76.804  d0/92.617d0) then
              nsta = 72
           ELSEIF(br .le. 77.684  d0/92.617d0) then
              nsta = 73
           ELSEIF(br .le. 78.304  d0/92.617d0) then
              nsta = 74
           ELSEIF(br .le. 78.438  d0/92.617d0) then
              nsta = 75
           ELSEIF(br .le. 78.686  d0/92.617d0) then
              nsta = 76
           ELSEIF(br .le. 79.436  d0/92.617d0) then
              nsta = 77
           ELSEIF(br .le. 79.535  d0/92.617d0) then
              nsta = 78
           ELSEIF(br .le. 79.865  d0/92.617d0) then
              nsta = 79
           ELSEIF(br .le. 79.994  d0/92.617d0) then
              nsta = 80
           ELSEIF(br .le. 80.385  d0/92.617d0) then
              nsta = 81
           ELSEIF(br .le. 80.813  d0/92.617d0) then
              nsta = 82
           ELSEIF(br .le. 80.91   d0/92.617d0) then
              nsta = 83
           ELSEIF(br .le. 80.947  d0/92.617d0) then
              nsta = 84
           ELSEIF(br .le. 81.013  d0/92.617d0) then
              nsta = 85
           ELSEIF(br .le. 82.223  d0/92.617d0) then
              nsta = 86
           ELSEIF(br .le. 82.641  d0/92.617d0) then
              nsta = 87
           ELSEIF(br .le. 83.146  d0/92.617d0) then
              nsta = 88
           ELSEIF(br .le. 83.628  d0/92.617d0) then
              nsta = 89
           ELSEIF(br .le. 83.873  d0/92.617d0) then
              nsta = 90
           ELSEIF(br .le. 84.343  d0/92.617d0) then
              nsta = 91
           ELSEIF(br .le. 84.489  d0/92.617d0) then
              nsta = 92
           ELSEIF(br .le. 84.582  d0/92.617d0) then
              nsta = 93
           ELSEIF(br .le. 84.724  d0/92.617d0) then
              nsta = 94
           ELSEIF(br .le. 84.799  d0/92.617d0) then
              nsta = 95
           ELSEIF(br .le. 84.915  d0/92.617d0) then
              nsta = 96
           ELSEIF(br .le. 85.335  d0/92.617d0) then
              nsta = 97
           ELSEIF(br .le. 86.045  d0/92.617d0) then
              nsta = 98
           ELSEIF(br .le. 86.375  d0/92.617d0) then
              nsta = 100
           ELSEIF(br .le. 87.015  d0/92.617d0) then
              nsta = 101
           ELSEIF(br .le. 88.065  d0/92.617d0) then
              nsta = 102
           ELSEIF(br .le. 88.217  d0/92.617d0) then
              nsta = 104
           ELSEIF(br .le. 88.547  d0/92.617d0) then
              nsta = 106
           ELSEIF(br .le. 88.877  d0/92.617d0) then
              nsta = 107
           ELSEIF(br .le. 89.117  d0/92.617d0) then
              nsta = 108
           ELSEIF(br .le. 89.317  d0/92.617d0) then
              nsta = 111
           ELSEIF(br .le. 89.49   d0/92.617d0) then
              nsta = 112
           ELSEIF(br .le. 89.72   d0/92.617d0) then
              nsta = 113
           ELSEIF(br .le. 90.017  d0/92.617d0) then
              nsta = 114
           ELSEIF(br .le. 90.207  d0/92.617d0) then
              nsta = 116
           ELSEIF(br .le. 90.327  d0/92.617d0) then
              nsta = 119
           ELSEIF(br .le. 90.557  d0/92.617d0) then
              nsta = 120
           ELSEIF(br .le. 90.627  d0/92.617d0) then
              nsta = 121
           ELSEIF(br .le. 92.257  d0/92.617d0) then
              nsta = 122
           ELSE
              nsta = 123
           ENDIF

         Case(30065 ) 
           IF(br .le. 0.47  d0/96.16d0) then
             nsta = 1
           ELSEIF(br .le. 48.27 d0/96.16d0) then
             nsta = 2
           ELSEIF(br .le. 54.99 d0/96.16d0) then
             nsta = 6
           ELSEIF(br .le. 61.91 d0/96.16d0) then
             nsta = 7
           ELSEIF(br .le. 65.95 d0/96.16d0) then
             nsta = 15
           ELSEIF(br .le. 66.38 d0/96.16d0) then
             nsta = 21
           ELSEIF(br .le. 71.44 d0/96.16d0) then
             nsta = 23
           ELSEIF(br .le. 71.97 d0/96.16d0) then
             nsta = 29
           ELSEIF(br .le. 76.57 d0/96.16d0) then
             nsta = 33
           ELSEIF(br .le. 77.53 d0/96.16d0) then
             nsta = 34
           ELSEIF(br .le. 81.3  d0/96.16d0) then
             nsta = 40
           ELSEIF(br .le. 82.7  d0/96.16d0) then
             nsta = 42
           ELSEIF(br .le. 83.21 d0/96.16d0) then
             nsta = 44
           ELSEIF(br .le. 84.1  d0/96.16d0) then
             nsta = 51
           ELSEIF(br .le. 84.67 d0/96.16d0) then
             nsta = 63
           ELSEIF(br .le. 85.27 d0/96.16d0) then
             nsta = 64
           ELSEIF(br .le. 86.56 d0/96.16d0) then
             nsta = 68
           ELSEIF(br .le. 88.85 d0/96.16d0) then
             nsta = 69
           ELSEIF(br .le. 89.72 d0/96.16d0) then
             nsta = 71
           ELSEIF(br .le. 90.15 d0/96.16d0) then
             nsta = 76
           ELSEIF(br .le. 90.59 d0/96.16d0) then
             nsta = 81
           ELSEIF(br .le. 91.4  d0/96.16d0) then
             nsta = 85
           ELSEIF(br .le. 92.21 d0/96.16d0) then
             nsta = 89
           ELSEIF(br .le. 92.73 d0/96.16d0) then
             nsta = 97
           ELSEIF(br .le. 93.17 d0/96.16d0) then
             nsta = 98
           ELSEIF(br .le. 93.83 d0/96.16d0) then
             nsta = 105
           ELSEIF(br .le. 94.31 d0/96.16d0) then
             nsta = 110
           ELSEIF(br .le. 94.9  d0/96.16d0) then
             nsta = 112
           ELSE
             nsta = 131
           ENDIF

         Case(30067 ) 
           IF(br .le. 36.3  d0/86.25d0) then
             nsta = 1
           ELSEIF(br .le. 60.05 d0/86.25d0) then
             nsta = 2
           ELSEIF(br .le. 73.7  d0/86.25d0) then
             nsta = 3
           ELSEIF(br .le. 83.65 d0/86.25d0) then
             nsta = 11
           ELSE
             nsta = 40
           ENDIF


         Case(30068 ) 
           IF(br .le. 2.893  d0/41.205d0) then
              nsta = 1
           ELSEIF(br .le. 4.936  d0/41.205d0) then
              nsta = 3
           ELSEIF(br .le. 6.623  d0/41.205d0) then
              nsta = 4
           ELSEIF(br .le. 6.645  d0/41.205d0) then
              nsta = 6
           ELSEIF(br .le. 6.856  d0/41.205d0) then
              nsta = 8
           ELSEIF(br .le. 7.063  d0/41.205d0) then
              nsta = 11
           ELSEIF(br .le. 9.682  d0/41.205d0) then
              nsta = 12
           ELSEIF(br .le. 10.064 d0/41.205d0) then
              nsta = 17
           ELSEIF(br .le. 14.564 d0/41.205d0) then
              nsta = 20
           ELSEIF(br .le. 15.136 d0/41.205d0) then
              nsta = 25
           ELSEIF(br .le. 17.462 d0/41.205d0) then
              nsta = 26
           ELSEIF(br .le. 18.785 d0/41.205d0) then
              nsta = 30
           ELSEIF(br .le. 19.028 d0/41.205d0) then
              nsta = 31
           ELSEIF(br .le. 19.338 d0/41.205d0) then
              nsta = 35
           ELSEIF(br .le. 19.536 d0/41.205d0) then
              nsta = 36
           ELSEIF(br .le. 20.031 d0/41.205d0) then
              nsta = 40
           ELSEIF(br .le. 20.647 d0/41.205d0) then
              nsta = 42
           ELSEIF(br .le. 20.984 d0/41.205d0) then
              nsta = 45
           ELSEIF(br .le. 21.348 d0/41.205d0) then
              nsta = 46
           ELSEIF(br .le. 21.825 d0/41.205d0) then
              nsta = 47
           ELSEIF(br .le. 23.58  d0/41.205d0) then
              nsta = 49
           ELSEIF(br .le. 23.872 d0/41.205d0) then
              nsta = 53
           ELSEIF(br .le. 24.029 d0/41.205d0) then
              nsta = 59
           ELSEIF(br .le. 24.627 d0/41.205d0) then
              nsta = 61
           ELSEIF(br .le. 25.216 d0/41.205d0) then
              nsta = 65
           ELSEIF(br .le. 26.125 d0/41.205d0) then
              nsta = 72
           ELSEIF(br .le. 27.169 d0/41.205d0) then
              nsta = 79
           ELSEIF(br .le. 28.222 d0/41.205d0) then
              nsta = 80
           ELSEIF(br .le. 29.545 d0/41.205d0) then
              nsta = 81
           ELSEIF(br .le. 30.08  d0/41.205d0) then
              nsta = 89
           ELSEIF(br .le. 30.381 d0/41.205d0) then
              nsta = 90
           ELSEIF(br .le. 31.15  d0/41.205d0) then
              nsta = 93
           ELSEIF(br .le. 31.487 d0/41.205d0) then
              nsta = 94
           ELSEIF(br .le. 32.45  d0/41.205d0) then
              nsta = 97
           ELSEIF(br .le. 33.566 d0/41.205d0) then
              nsta = 98
           ELSEIF(br .le. 34.016 d0/41.205d0) then
              nsta = 99
           ELSEIF(br .le. 34.448 d0/41.205d0) then
              nsta = 101
           ELSEIF(br .le. 35.001 d0/41.205d0) then
              nsta = 106
           ELSEIF(br .le. 35.892 d0/41.205d0) then
              nsta = 108
           ELSEIF(br .le. 36.963 d0/41.205d0) then
              nsta = 109
           ELSEIF(br .le. 37.57  d0/41.205d0) then
              nsta = 110
           ELSEIF(br .le. 38.038 d0/41.205d0) then
              nsta = 111
           ELSEIF(br .le. 38.78  d0/41.205d0) then
              nsta = 112
           ELSEIF(br .le. 39.68  d0/41.205d0) then
              nsta = 113
           ELSEIF(br .le. 40.373 d0/41.205d0) then
              nsta = 115
           ELSE                            
              nsta = 118
           ENDIF

         Case(30069 ) 
           IF(br .le. 4.9  d0/55.8d0) then
             nsta = 0
           ELSEIF(br .le. 8.7  d0/55.8d0) then
             nsta = 1
           ELSEIF(br .le. 29.5 d0/55.8d0) then
             nsta = 6
           ELSEIF(br .le. 31.4 d0/55.8d0) then
             nsta = 10
           ELSEIF(br .le. 34.7 d0/55.8d0) then
             nsta = 15
           ELSEIF(br .le. 36.9 d0/55.8d0) then
             nsta = 22
           ELSEIF(br .le. 38.2 d0/55.8d0) then
             nsta = 27
           ELSEIF(br .le. 40.4 d0/55.8d0) then
             nsta = 29
           ELSEIF(br .le. 49.7 d0/55.8d0) then
             nsta = 35
           ELSEIF(br .le. 53.3 d0/55.8d0) then
             nsta = 36
           ELSEIF(br .le. 54.5 d0/55.8d0) then
             nsta = 37
           ELSEIF(br .le. 55.0 d0/55.8d0) then
             nsta = 40
           ELSE
             nsta = 41
           ENDIF

         Case(32071 ) 
              nsta = 193

         Case(32073 ) 
           IF(br .le. 1.43  d0/54.76d0) then
             nsta = 2
           ELSEIF(br .le. 6.23  d0/54.76d0) then
             nsta = 5
           ELSEIF(br .le. 16.75 d0/54.76d0) then
             nsta = 6
           ELSEIF(br .le. 18.29 d0/54.76d0) then
             nsta = 10
           ELSEIF(br .le. 19.04 d0/54.76d0) then
             nsta = 22
           ELSEIF(br .le. 21.25 d0/54.76d0) then
             nsta = 23
           ELSEIF(br .le. 26.07 d0/54.76d0) then
             nsta = 28
           ELSEIF(br .le. 30.21 d0/54.76d0) then
             nsta = 31
           ELSEIF(br .le. 39.56 d0/54.76d0) then
             nsta = 35
           ELSEIF(br .le. 39.98 d0/54.76d0) then
             nsta = 39
           ELSEIF(br .le. 40.62 d0/54.76d0) then
             nsta = 41
           ELSEIF(br .le. 41.24 d0/54.76d0) then
             nsta = 49
           ELSEIF(br .le. 42.02 d0/54.76d0) then
             nsta = 58
           ELSEIF(br .le. 42.82 d0/54.76d0) then
             nsta = 60
           ELSEIF(br .le. 43.34 d0/54.76d0) then
             nsta = 61
           ELSEIF(br .le. 44.04 d0/54.76d0) then
             nsta = 64
           ELSEIF(br .le. 44.7  d0/54.76d0) then
             nsta = 69
           ELSEIF(br .le. 45.6  d0/54.76d0) then
             nsta = 70
           ELSEIF(br .le. 46.18 d0/54.76d0) then
             nsta = 73
           ELSEIF(br .le. 47.84 d0/54.76d0) then
             nsta = 75
           ELSEIF(br .le. 49.64 d0/54.76d0) then
             nsta = 80
           ELSEIF(br .le. 51.36 d0/54.76d0) then
             nsta = 81
           ELSEIF(br .le. 52.4  d0/54.76d0) then
             nsta = 83
           ELSEIF(br .le. 54.1  d0/54.76d0) then
             nsta = 86
           ELSE
             nsta = 88
           ENDIF

         Case(32074 ) 
           IF(br .le. 0.06    d0/ 26.24d0) then
             nsta = 1
           ELSEIF(br .le. 0.07    d0/ 26.24d0) then
             nsta = 2
           ELSEIF(br .le. 1.09    d0/ 26.24d0) then
             nsta = 3
           ELSEIF(br .le. 1.89    d0/ 26.24d0) then
             nsta = 5
           ELSEIF(br .le. 2.96    d0/ 26.24d0) then
             nsta = 9
           ELSEIF(br .le. 3.23    d0/ 26.24d0) then
             nsta = 16
           ELSEIF(br .le. 3.9     d0/ 26.24d0) then
             nsta = 21
           ELSEIF(br .le. 4.33    d0/ 26.24d0) then
             nsta = 25
           ELSEIF(br .le. 6.87    d0/ 26.24d0) then
             nsta = 32
           ELSEIF(br .le. 7.65    d0/ 26.24d0) then
             nsta = 37
           ELSEIF(br .le. 7.92    d0/ 26.24d0) then
             nsta = 41
           ELSEIF(br .le. 8.24    d0/ 26.24d0) then
             nsta = 42
           ELSEIF(br .le. 8.42    d0/ 26.24d0) then
             nsta = 44
           ELSEIF(br .le. 9.88    d0/ 26.24d0) then
             nsta = 46
           ELSEIF(br .le. 10.34   d0/ 26.24d0) then
             nsta = 50
           ELSEIF(br .le. 10.83   d0/ 26.24d0) then
             nsta = 55
           ELSEIF(br .le. 10.96   d0/ 26.24d0) then
             nsta = 58
           ELSEIF(br .le. 11.68   d0/ 26.24d0) then
             nsta = 64
           ELSEIF(br .le. 12.39   d0/ 26.24d0) then
             nsta = 66
           ELSEIF(br .le. 13.28   d0/ 26.24d0) then
             nsta = 69
           ELSEIF(br .le. 14.     d0/ 26.24d0) then
             nsta = 71
           ELSEIF(br .le. 14.11   d0/ 26.24d0) then
             nsta = 86
           ELSEIF(br .le. 14.47   d0/ 26.24d0) then
             nsta = 90
           ELSEIF(br .le. 14.85   d0/ 26.24d0) then
             nsta = 99
           ELSEIF(br .le. 15.105  d0/ 26.24d0) then
             nsta = 100
           ELSEIF(br .le. 15.555  d0/ 26.24d0) then
             nsta = 103
           ELSEIF(br .le. 15.915  d0/ 26.24d0) then
             nsta = 104
           ELSEIF(br .le. 16.055  d0/ 26.24d0) then
             nsta = 117
           ELSEIF(br .le. 17.095  d0/ 26.24d0) then
             nsta = 120
           ELSEIF(br .le. 17.855  d0/ 26.24d0) then
             nsta = 125
           ELSEIF(br .le. 18.035  d0/ 26.24d0) then
             nsta = 127
           ELSEIF(br .le. 18.065  d0/ 26.24d0) then
             nsta = 138
           ELSEIF(br .le. 18.38   d0/ 26.24d0) then
             nsta = 139
           ELSEIF(br .le. 18.59   d0/ 26.24d0) then
             nsta = 146
           ELSEIF(br .le. 18.71   d0/ 26.24d0) then
             nsta = 147
           ELSEIF(br .le. 19.04   d0/ 26.24d0) then
             nsta = 150
           ELSEIF(br .le. 19.46   d0/ 26.24d0) then
             nsta = 162
           ELSEIF(br .le. 19.71   d0/ 26.24d0) then
             nsta = 165
           ELSEIF(br .le. 19.93   d0/ 26.24d0) then
             nsta = 166
           ELSEIF(br .le. 20.23   d0/ 26.24d0) then
             nsta = 167
           ELSEIF(br .le. 20.66   d0/ 26.24d0) then
             nsta = 171
           ELSEIF(br .le. 21.26   d0/ 26.24d0) then
             nsta = 179
           ELSEIF(br .le. 21.59   d0/ 26.24d0) then
             nsta = 186
           ELSEIF(br .le. 21.94   d0/ 26.24d0) then
             nsta = 196
           ELSEIF(br .le. 22.135  d0/ 26.24d0) then
             nsta = 211
           ELSEIF(br .le. 22.285  d0/ 26.24d0) then
             nsta = 215
           ELSEIF(br .le. 22.465  d0/ 26.24d0) then
             nsta = 238
           ELSEIF(br .le. 22.795  d0/ 26.24d0) then
             nsta = 241
           ELSEIF(br .le. 23.515  d0/ 26.24d0) then
             nsta = 246
           ELSEIF(br .le. 23.765  d0/ 26.24d0) then
             nsta = 247
           ELSEIF(br .le. 24.875  d0/ 26.24d0) then
             nsta = 253
           ELSEIF(br .le. 25.565  d0/ 26.24d0) then
             nsta = 254
           ELSE
             nsta = 256
           ENDIF

         Case(32075 ) 
             nsta = 122

         Case(32077 ) 
             nsta = 64

         Case(33076)
               IF(br .le. 0.0043/3.0705) then
             nsta = 0
           ELSEIF(br .le. 0.0403/3.0705) then
             nsta = 1
           ELSEIF(br .le. 0.057/3.0705) then
             nsta = 2
           ELSEIF(br .le. 0.0697/3.0705) then
             nsta = 3
           ELSEIF(br .le. 0.0878/3.0705) then
             nsta = 5
           ELSEIF(br .le. 0.0915/3.0705) then
             nsta = 6
           ELSEIF(br .le. 0.1365/3.0705) then
             nsta = 8
           ELSEIF(br .le. 0.1468/3.0705) then
             nsta = 9
           ELSEIF(br .le. 0.2002/3.0705) then
             nsta = 12
           ELSEIF(br .le. 0.3042/3.0705) then
             nsta = 13
           ELSEIF(br .le. 0.3064/3.0705) then
             nsta = 14
           ELSEIF(br .le. 0.3194/3.0705) then
             nsta = 15
           ELSEIF(br .le. 0.3242/3.0705) then
             nsta = 16
           ELSEIF(br .le. 0.3291/3.0705) then
             nsta = 18
           ELSEIF(br .le. 0.3901/3.0705) then
             nsta = 19
           ELSEIF(br .le. 0.3925/3.0705) then
             nsta = 20
           ELSEIF(br .le. 0.4087/3.0705) then
             nsta = 21
           ELSEIF(br .le. 0.4255/3.0705) then
             nsta = 24
           ELSEIF(br .le. 0.4416/3.0705) then
             nsta = 26
           ELSEIF(br .le. 0.4549/3.0705) then
             nsta = 27
           ELSEIF(br .le. 1.0149/3.0705) then
             nsta = 29
           ELSEIF(br .le. 1.1749/3.0705) then
             nsta = 30
           ELSEIF(br .le. 1.1882/3.0705) then
             nsta = 31
           ELSEIF(br .le. 1.2025/3.0705) then
             nsta = 32
           ELSEIF(br .le. 1.2109/3.0705) then
             nsta = 34
           ELSEIF(br .le. 1.221/3.0705) then
             nsta = 35
           ELSEIF(br .le. 1.2319/3.0705) then
             nsta = 37
           ELSEIF(br .le. 1.2565/3.0705) then
             nsta = 38
           ELSEIF(br .le. 1.2792/3.0705) then
             nsta = 40
           ELSEIF(br .le. 1.2842/3.0705) then
             nsta = 42
           ELSEIF(br .le. 1.3146/3.0705) then
             nsta = 44
           ELSEIF(br .le. 1.317/3.0705) then
             nsta = 45
           ELSEIF(br .le. 1.3542/3.0705) then
             nsta = 46
           ELSEIF(br .le. 1.3992/3.0705) then
             nsta = 48
           ELSEIF(br .le. 1.4262/3.0705) then
             nsta = 50
           ELSEIF(br .le. 1.4299/3.0705) then
             nsta = 52
           ELSEIF(br .le. 1.4354/3.0705) then
             nsta = 53
           ELSEIF(br .le. 1.4762/3.0705) then
             nsta = 54
           ELSEIF(br .le. 1.5078/3.0705) then
             nsta = 55
           ELSEIF(br .le. 1.5201/3.0705) then
             nsta = 56
           ELSEIF(br .le. 1.5312/3.0705) then
             nsta = 57
           ELSEIF(br .le. 1.5337/3.0705) then
             nsta = 58
           ELSEIF(br .le. 1.5647/3.0705) then
             nsta = 59
           ELSEIF(br .le. 1.5867/3.0705) then
             nsta = 60
           ELSEIF(br .le. 1.6187/3.0705) then
             nsta = 61
           ELSEIF(br .le. 1.6516/3.0705) then
             nsta = 62
           ELSEIF(br .le. 1.6573/3.0705) then
             nsta = 63
           ELSEIF(br .le. 1.6847/3.0705) then
             nsta = 64
           ELSEIF(br .le. 1.6906/3.0705) then
             nsta = 65
           ELSEIF(br .le. 1.711/3.0705) then
             nsta = 66
           ELSEIF(br .le. 1.721/3.0705) then
             nsta = 67
           ELSEIF(br .le. 1.806/3.0705) then
             nsta = 69
           ELSEIF(br .le. 1.83/3.0705) then
             nsta = 70
           ELSEIF(br .le. 1.894/3.0705) then
             nsta = 71
           ELSEIF(br .le. 1.9039/3.0705) then
             nsta = 72
           ELSEIF(br .le. 1.9071/3.0705) then
             nsta = 73
           ELSEIF(br .le. 1.9484/3.0705) then
             nsta = 74
           ELSEIF(br .le. 1.9604/3.0705) then
             nsta = 75
           ELSEIF(br .le. 1.9764/3.0705) then
             nsta = 76
           ELSEIF(br .le. 2.0028/3.0705) then
             nsta = 77
           ELSEIF(br .le. 2.0133/3.0705) then
             nsta = 78
           ELSEIF(br .le. 2.0273/3.0705) then
             nsta = 79
           ELSEIF(br .le. 2.0357/3.0705) then
             nsta = 80
           ELSEIF(br .le. 2.0444/3.0705) then
             nsta = 82
           ELSEIF(br .le. 2.0531/3.0705) then
             nsta = 83
           ELSEIF(br .le. 2.0791/3.0705) then
             nsta = 84
           ELSEIF(br .le. 2.0991/3.0705) then
             nsta = 85
           ELSEIF(br .le. 2.1152/3.0705) then
             nsta = 86
           ELSEIF(br .le. 2.1376/3.0705) then
             nsta = 87
           ELSEIF(br .le. 2.1673/3.0705) then
             nsta = 88
           ELSEIF(br .le. 2.1763/3.0705) then
             nsta = 89
           ELSEIF(br .le. 2.1962/3.0705) then
             nsta = 90
           ELSEIF(br .le. 2.2172/3.0705) then
             nsta = 91
           ELSEIF(br .le. 2.2202/3.0705) then
             nsta = 92
           ELSEIF(br .le. 2.2321/3.0705) then
             nsta = 93
           ELSEIF(br .le. 2.2464/3.0705) then
             nsta = 94
           ELSEIF(br .le. 2.2539/3.0705) then
             nsta = 95
           ELSEIF(br .le. 2.2667/3.0705) then
             nsta = 96
           ELSEIF(br .le. 2.3171/3.0705) then
             nsta = 97
           ELSEIF(br .le. 2.3447/3.0705) then
             nsta = 98
           ELSEIF(br .le. 2.3597/3.0705) then
             nsta = 99
           ELSEIF(br .le. 2.3815/3.0705) then
             nsta = 100
           ELSEIF(br .le. 2.391/3.0705) then
             nsta = 101
           ELSEIF(br .le. 2.412/3.0705) then
             nsta = 102
           ELSEIF(br .le. 2.4139/3.0705) then
             nsta = 103
           ELSEIF(br .le. 2.4386/3.0705) then
             nsta = 104
           ELSEIF(br .le. 2.4461/3.0705) then
             nsta = 105
           ELSEIF(br .le. 2.4548/3.0705) then
             nsta = 106
           ELSEIF(br .le. 2.4808/3.0705) then
             nsta = 107
           ELSEIF(br .le. 2.529/3.0705) then
             nsta = 108
           ELSEIF(br .le. 2.544/3.0705) then
             nsta = 109
           ELSEIF(br .le. 2.5483/3.0705) then
             nsta = 110
           ELSEIF(br .le. 2.5507/3.0705) then
             nsta = 111
           ELSEIF(br .le. 2.5667/3.0705) then
             nsta = 112
           ELSEIF(br .le. 2.6146/3.0705) then
             nsta = 113
           ELSEIF(br .le. 2.6376/3.0705) then
             nsta = 114
           ELSEIF(br .le. 2.679/3.0705) then
             nsta = 115
           ELSEIF(br .le. 2.705/3.0705) then
             nsta = 116
           ELSEIF(br .le. 2.724/3.0705) then
             nsta = 117
           ELSEIF(br .le. 2.7412/3.0705) then
             nsta = 118
           ELSEIF(br .le. 2.7531/3.0705) then
             nsta = 119
           ELSEIF(br .le. 2.7585/3.0705) then
             nsta = 120
           ELSEIF(br .le. 2.7766/3.0705) then
             nsta = 121
           ELSEIF(br .le. 2.7836/3.0705) then
             nsta = 122
           ELSEIF(br .le. 2.7986/3.0705) then
             nsta = 123
           ELSEIF(br .le. 2.8089/3.0705) then
             nsta = 124
           ELSEIF(br .le. 2.8227/3.0705) then
             nsta = 125
           ELSEIF(br .le. 2.8417/3.0705) then
             nsta = 126
           ELSEIF(br .le. 2.8771/3.0705) then
             nsta = 127
           ELSEIF(br .le. 2.8952/3.0705) then
             nsta = 128
           ELSEIF(br .le. 2.9083/3.0705) then
             nsta = 129
           ELSEIF(br .le. 3.0593/3.0705) then
             nsta = 130
           ELSE
             nsta = 131
           ENDIF
         Case(34077 ) 
             nsta = 163
         Case(34078)
               IF(br .le. 0.0221/1.4014) then
             nsta = 0
           ELSEIF(br .le. 0.2421/1.4014) then
             nsta = 1
           ELSEIF(br .le. 0.3921/1.4014) then
             nsta = 2
           ELSEIF(br .le. 0.3991/1.4014) then
             nsta = 3
           ELSEIF(br .le. 0.401/1.4014) then
             nsta = 5
           ELSEIF(br .le. 0.449/1.4014) then
             nsta = 7
           ELSEIF(br .le. 0.503/1.4014) then
             nsta = 11
           ELSEIF(br .le. 0.561/1.4014) then
             nsta = 12
           ELSEIF(br .le. 0.5769/1.4014) then
             nsta = 15
           ELSEIF(br .le. 0.5801/1.4014) then
             nsta = 21
           ELSEIF(br .le. 0.5962/1.4014) then
             nsta = 26
           ELSEIF(br .le. 0.615/1.4014) then
             nsta = 27
           ELSEIF(br .le. 0.6292/1.4014) then
             nsta = 30
           ELSEIF(br .le. 0.6587/1.4014) then
             nsta = 34
           ELSEIF(br .le. 0.66/1.4014) then
             nsta = 40
           ELSEIF(br .le. 0.6698/1.4014) then
             nsta = 47
           ELSEIF(br .le. 0.6818/1.4014) then
             nsta = 49
           ELSEIF(br .le. 0.6919/1.4014) then
             nsta = 50
           ELSEIF(br .le. 0.7479/1.4014) then
             nsta = 51
           ELSEIF(br .le. 0.7849/1.4014) then
             nsta = 57
           ELSEIF(br .le. 0.7957/1.4014) then
             nsta = 61
           ELSEIF(br .le. 0.8058/1.4014) then
             nsta = 62
           ELSEIF(br .le. 0.8153/1.4014) then
             nsta = 68
           ELSEIF(br .le. 0.8387/1.4014) then
             nsta = 73
           ELSEIF(br .le. 0.8545/1.4014) then
             nsta = 75
           ELSEIF(br .le. 0.8802/1.4014) then
             nsta = 78
           ELSEIF(br .le. 0.895/1.4014) then
             nsta = 81
           ELSEIF(br .le. 0.9115/1.4014) then
             nsta = 87
           ELSEIF(br .le. 0.9131/1.4014) then
             nsta = 89
           ELSEIF(br .le. 0.9601/1.4014) then
             nsta = 91
           ELSEIF(br .le. 0.975/1.4014) then
             nsta = 92
           ELSEIF(br .le. 0.9782/1.4014) then
             nsta = 96
           ELSEIF(br .le. 0.9909/1.4014) then
             nsta = 101
           ELSEIF(br .le. 1.0349/1.4014) then
             nsta = 103
           ELSEIF(br .le. 1.0779/1.4014) then
             nsta = 109
           ELSEIF(br .le. 1.0842/1.4014) then
             nsta = 111
           ELSEIF(br .le. 1.0947/1.4014) then
             nsta = 112
           ELSEIF(br .le. 1.1039/1.4014) then
             nsta = 114
           ELSEIF(br .le. 1.1143/1.4014) then
             nsta = 116
           ELSEIF(br .le. 1.1434/1.4014) then
             nsta = 120
           ELSEIF(br .le. 1.1602/1.4014) then
             nsta = 122
           ELSEIF(br .le. 1.1697/1.4014) then
             nsta = 126
           ELSEIF(br .le. 1.1741/1.4014) then
             nsta = 134
           ELSEIF(br .le. 1.201/1.4014) then
             nsta = 135
           ELSEIF(br .le. 1.2127/1.4014) then
             nsta = 136
           ELSEIF(br .le. 1.22/1.4014) then
             nsta = 137
           ELSEIF(br .le. 1.2401/1.4014) then
             nsta = 138
           ELSEIF(br .le. 1.2566/1.4014) then
             nsta = 141
           ELSEIF(br .le. 1.2683/1.4014) then
             nsta = 142
           ELSEIF(br .le. 1.2759/1.4014) then
             nsta = 143
           ELSEIF(br .le. 1.2781/1.4014) then
             nsta = 151
           ELSEIF(br .le. 1.2794/1.4014) then
             nsta = 152
           ELSEIF(br .le. 1.2921/1.4014) then
             nsta = 153
           ELSEIF(br .le. 1.3016/1.4014) then
             nsta = 154
           ELSEIF(br .le. 1.3048/1.4014) then
             nsta = 155
           ELSEIF(br .le. 1.313/1.4014) then
             nsta = 156
           ELSEIF(br .le. 1.3158/1.4014) then
             nsta = 158
           ELSEIF(br .le. 1.3177/1.4014) then
             nsta = 160
           ELSEIF(br .le. 1.3357/1.4014) then
             nsta = 162
           ELSEIF(br .le. 1.3454/1.4014) then
             nsta = 164
           ELSEIF(br .le. 1.3521/1.4014) then
             nsta = 165
           ELSEIF(br .le. 1.3559/1.4014) then
             nsta = 169
           ELSEIF(br .le. 1.3629/1.4014) then
             nsta = 170
           ELSEIF(br .le. 1.3711/1.4014) then
             nsta = 171
           ELSEIF(br .le. 1.3758/1.4014) then
             nsta = 172
           ELSEIF(br .le. 1.3831/1.4014) then
             nsta = 173
           ELSEIF(br .le. 1.3913/1.4014) then
             nsta = 175
           ELSEIF(br .le. 1.3951/1.4014) then
             nsta = 176
           ELSE
             nsta = 178
           ENDIF
         Case(34079)
               IF(br .le. 0.00085/0.002185) then
             nsta = 1
           ELSEIF(br .le. 0.000875/0.002185) then
             nsta = 5
           ELSEIF(br .le. 0.001045/0.002185) then
             nsta = 6
           ELSEIF(br .le. 0.001059/0.002185) then
             nsta = 22
           ELSEIF(br .le. 0.001359/0.002185) then
             nsta = 23
           ELSEIF(br .le. 0.001426/0.002185) then
             nsta = 39
           ELSEIF(br .le. 0.001441/0.002185) then
             nsta = 41
           ELSEIF(br .le. 0.001581/0.002185) then
             nsta = 47
           ELSEIF(br .le. 0.001741/0.002185) then
             nsta = 51
           ELSEIF(br .le. 0.00181/0.002185) then
             nsta = 83
           ELSEIF(br .le. 0.001843/0.002185) then
             nsta = 84
           ELSEIF(br .le. 0.001865/0.002185) then
             nsta = 89
           ELSEIF(br .le. 0.001933/0.002185) then
             nsta = 90
           ELSEIF(br .le. 0.001979/0.002185) then
             nsta = 92
           ELSEIF(br .le. 0.002006/0.002185) then
             nsta = 94
           ELSEIF(br .le. 0.002032/0.002185) then
             nsta = 95
           ELSEIF(br .le. 0.00206/0.002185) then
             nsta = 97
           ELSEIF(br .le. 0.002097/0.002185) then
             nsta = 98
           ELSEIF(br .le. 0.00212/0.002185) then
             nsta = 103
           ELSEIF(br .le. 0.002154/0.002185) then
             nsta = 104
           ELSE
             nsta = 106
           ENDIF
         Case(34081)
               IF(br .le. 0.1/0.1402) then
             nsta = 3
           ELSEIF(br .le. 0.1013/0.1402) then
             nsta = 11
           ELSEIF(br .le. 0.1019/0.1402) then
             nsta = 13
           ELSEIF(br .le. 0.1032/0.1402) then
             nsta = 15
           ELSEIF(br .le. 0.1182/0.1402) then
             nsta = 17
           ELSEIF(br .le. 0.1236/0.1402) then
             nsta = 24
           ELSEIF(br .le. 0.1294/0.1402) then
             nsta = 26
           ELSEIF(br .le. 0.1318/0.1402) then
             nsta = 30
           ELSEIF(br .le. 0.1368/0.1402) then
             nsta = 40
           ELSE
             nsta = 56
           ENDIF
         Case(35080)
               IF(br .le. 0.006/0.5248) then
             nsta = 0
           ELSEIF(br .le. 0.0121/0.5248) then
             nsta = 3
           ELSEIF(br .le. 0.0189/0.5248) then
             nsta = 4
           ELSEIF(br .le. 0.0282/0.5248) then
             nsta = 5
           ELSEIF(br .le. 0.1362/0.5248) then
             nsta = 8
           ELSEIF(br .le. 0.1405/0.5248) then
             nsta = 10
           ELSEIF(br .le. 0.1513/0.5248) then
             nsta = 14
           ELSEIF(br .le. 0.2008/0.5248) then
             nsta = 20
           ELSEIF(br .le. 0.2033/0.5248) then
             nsta = 25
           ELSEIF(br .le. 0.2136/0.5248) then
             nsta = 37
           ELSEIF(br .le. 0.2176/0.5248) then
             nsta = 39
           ELSEIF(br .le. 0.233/0.5248) then
             nsta = 42
           ELSEIF(br .le. 0.2896/0.5248) then
             nsta = 46
           ELSEIF(br .le. 0.3343/0.5248) then
             nsta = 51
           ELSEIF(br .le. 0.3391/0.5248) then
             nsta = 53
           ELSEIF(br .le. 0.3501/0.5248) then
             nsta = 54
           ELSEIF(br .le. 0.3602/0.5248) then
             nsta = 59
           ELSEIF(br .le. 0.3654/0.5248) then
             nsta = 65
           ELSEIF(br .le. 0.3822/0.5248) then
             nsta = 73
           ELSEIF(br .le. 0.4143/0.5248) then
             nsta = 75
           ELSEIF(br .le. 0.4449/0.5248) then
             nsta = 77
           ELSEIF(br .le. 0.4767/0.5248) then
             nsta = 78
           ELSEIF(br .le. 0.5052/0.5248) then
             nsta = 84
           ELSE
             nsta = 87
           ENDIF
         Case(35082)
               IF(br .le. 0.0101/0.2471) then
             nsta = 3
           ELSEIF(br .le. 0.0351/0.2471) then
             nsta = 4
           ELSEIF(br .le. 0.0589/0.2471) then
             nsta = 6
           ELSEIF(br .le. 0.0597/0.2471) then
             nsta = 7
           ELSEIF(br .le. 0.0656/0.2471) then
             nsta = 8
           ELSEIF(br .le. 0.0682/0.2471) then
             nsta = 9
           ELSEIF(br .le. 0.0745/0.2471) then
             nsta = 11
           ELSEIF(br .le. 0.0759/0.2471) then
             nsta = 13
           ELSEIF(br .le. 0.0816/0.2471) then
             nsta = 14
           ELSEIF(br .le. 0.0878/0.2471) then
             nsta = 17
           ELSEIF(br .le. 0.1264/0.2471) then
             nsta = 18
           ELSEIF(br .le. 0.1368/0.2471) then
             nsta = 23
           ELSEIF(br .le. 0.1428/0.2471) then
             nsta = 24
           ELSEIF(br .le. 0.1508/0.2471) then
             nsta = 26
           ELSEIF(br .le. 0.1547/0.2471) then
             nsta = 27
           ELSEIF(br .le. 0.1588/0.2471) then
             nsta = 31
           ELSEIF(br .le. 0.1625/0.2471) then
             nsta = 32
           ELSEIF(br .le. 0.1953/0.2471) then
             nsta = 33
           ELSEIF(br .le. 0.2089/0.2471) then
             nsta = 34
           ELSEIF(br .le. 0.2153/0.2471) then
             nsta = 37
           ELSEIF(br .le. 0.2303/0.2471) then
             nsta = 38
           ELSE
             nsta = 39
           ENDIF
         Case(36084)
               IF(br .le. 0.068/5.256) then
             nsta = 1
           ELSEIF(br .le. 0.083/5.256) then
             nsta = 3
           ELSEIF(br .le. 0.273/5.256) then
             nsta = 4
           ELSEIF(br .le. 0.285/5.256) then
             nsta = 5
           ELSEIF(br .le. 0.341/5.256) then
             nsta = 8
           ELSEIF(br .le. 0.451/5.256) then
             nsta = 10
           ELSEIF(br .le. 0.527/5.256) then
             nsta = 13
           ELSEIF(br .le. 0.579/5.256) then
             nsta = 17
           ELSEIF(br .le. 0.829/5.256) then
             nsta = 21
           ELSEIF(br .le. 0.883/5.256) then
             nsta = 26
           ELSEIF(br .le. 1.423/5.256) then
             nsta = 29
           ELSEIF(br .le. 2.723/5.256) then
             nsta = 30
           ELSEIF(br .le. 3.363/5.256) then
             nsta = 33
           ELSEIF(br .le. 4.243/5.256) then
             nsta = 40
           ELSEIF(br .le. 4.723/5.256) then
             nsta = 44
           ELSEIF(br .le. 4.993/5.256) then
             nsta = 45
           ELSEIF(br .le. 5.173/5.256) then
             nsta = 50
           ELSE
             nsta = 51
           ENDIF
         Case(37086)
               IF(br .le. 0.0022/0.06804) then
             nsta = 0
           ELSEIF(br .le. 0.00255/0.06804) then
             nsta = 1
           ELSEIF(br .le. 0.00466/0.06804) then
             nsta = 3
           ELSEIF(br .le. 0.00483/0.06804) then
             nsta = 5
           ELSEIF(br .le. 0.00558/0.06804) then
             nsta = 6
           ELSEIF(br .le. 0.01698/0.06804) then
             nsta = 7
           ELSEIF(br .le. 0.01797/0.06804) then
             nsta = 10
           ELSEIF(br .le. 0.01854/0.06804) then
             nsta = 13
           ELSEIF(br .le. 0.02444/0.06804) then
             nsta = 15
           ELSEIF(br .le. 0.02528/0.06804) then
             nsta = 17
           ELSEIF(br .le. 0.02657/0.06804) then
             nsta = 19
           ELSEIF(br .le. 0.02712/0.06804) then
             nsta = 20
           ELSEIF(br .le. 0.02873/0.06804) then
             nsta = 26
           ELSEIF(br .le. 0.03513/0.06804) then
             nsta = 29
           ELSEIF(br .le. 0.03574/0.06804) then
             nsta = 30
           ELSEIF(br .le. 0.03622/0.06804) then
             nsta = 32
           ELSEIF(br .le. 0.03698/0.06804) then
             nsta = 33
           ELSEIF(br .le. 0.03781/0.06804) then
             nsta = 34
           ELSEIF(br .le. 0.04421/0.06804) then
             nsta = 37
           ELSEIF(br .le. 0.04586/0.06804) then
             nsta = 38
           ELSEIF(br .le. 0.05076/0.06804) then
             nsta = 39
           ELSEIF(br .le. 0.05224/0.06804) then
             nsta = 41
           ELSEIF(br .le. 0.05397/0.06804) then
             nsta = 43
           ELSEIF(br .le. 0.05504/0.06804) then
             nsta = 44
           ELSEIF(br .le. 0.05564/0.06804) then
             nsta = 45
           ELSEIF(br .le. 0.05617/0.06804) then
             nsta = 47
           ELSEIF(br .le. 0.05977/0.06804) then
             nsta = 50
           ELSEIF(br .le. 0.06029/0.06804) then
             nsta = 51
           ELSEIF(br .le. 0.06161/0.06804) then
             nsta = 52
           ELSEIF(br .le. 0.06258/0.06804) then
             nsta = 54
           ELSEIF(br .le. 0.06728/0.06804) then
             nsta = 55
           ELSE
             nsta = 56
           ENDIF
         Case(37088)
               IF(br .le. 0.00097/0.01233) then
             nsta = 0
           ELSEIF(br .le. 0.00173/0.01233) then
             nsta = 1
           ELSEIF(br .le. 0.0039/0.01233) then
             nsta = 2
           ELSEIF(br .le. 0.00459/0.01233) then
             nsta = 5
           ELSEIF(br .le. 0.00505/0.01233) then
             nsta = 6
           ELSEIF(br .le. 0.00681/0.01233) then
             nsta = 9
           ELSEIF(br .le. 0.00973/0.01233) then
             nsta = 20
           ELSEIF(br .le. 0.01049/0.01233) then
             nsta = 22
           ELSE
             nsta = 32
           ENDIF
         Case(38087)
               IF(br .le. 0.026/0.0897209) then
             nsta = 1
           ELSEIF(br .le. 0.0260205/0.0897209) then
             nsta = 3
           ELSEIF(br .le. 0.0260259/0.0897209) then
             nsta = 4
           ELSEIF(br .le. 0.0282959/0.0897209) then
             nsta = 9
           ELSEIF(br .le. 0.0297659/0.0897209) then
             nsta = 11
           ELSEIF(br .le. 0.0335759/0.0897209) then
             nsta = 14
           ELSEIF(br .le. 0.0357659/0.0897209) then
             nsta = 25
           ELSEIF(br .le. 0.0372759/0.0897209) then
             nsta = 33
           ELSEIF(br .le. 0.0373399/0.0897209) then
             nsta = 37
           ELSEIF(br .le. 0.0375499/0.0897209) then
             nsta = 38
           ELSEIF(br .le. 0.0377109/0.0897209) then
             nsta = 40
           ELSEIF(br .le. 0.0437109/0.0897209) then
             nsta = 41
           ELSEIF(br .le. 0.0541109/0.0897209) then
             nsta = 44
           ELSEIF(br .le. 0.0555309/0.0897209) then
             nsta = 48
           ELSEIF(br .le. 0.0568009/0.0897209) then
             nsta = 50
           ELSEIF(br .le. 0.0570379/0.0897209) then
             nsta = 53
           ELSEIF(br .le. 0.0574279/0.0897209) then
             nsta = 57
           ELSEIF(br .le. 0.0576449/0.0897209) then
             nsta = 58
           ELSEIF(br .le. 0.0576639/0.0897209) then
             nsta = 61
           ELSEIF(br .le. 0.0578549/0.0897209) then
             nsta = 62
           ELSEIF(br .le. 0.0578979/0.0897209) then
             nsta = 63
           ELSEIF(br .le. 0.0580029/0.0897209) then
             nsta = 66
           ELSEIF(br .le. 0.0580869/0.0897209) then
             nsta = 71
           ELSEIF(br .le. 0.0581189/0.0897209) then
             nsta = 72
           ELSEIF(br .le. 0.0581749/0.0897209) then
             nsta = 77
           ELSEIF(br .le. 0.0583519/0.0897209) then
             nsta = 81
           ELSEIF(br .le. 0.0587219/0.0897209) then
             nsta = 83
           ELSEIF(br .le. 0.0587409/0.0897209) then
             nsta = 86
           ELSEIF(br .le. 0.0606409/0.0897209) then
             nsta = 87
           ELSEIF(br .le. 0.0608859/0.0897209) then
             nsta = 90
           ELSEIF(br .le. 0.0610669/0.0897209) then
             nsta = 96
           ELSEIF(br .le. 0.0613219/0.0897209) then
             nsta = 97
           ELSEIF(br .le. 0.0614419/0.0897209) then
             nsta = 100
           ELSEIF(br .le. 0.0619619/0.0897209) then
             nsta = 103
           ELSEIF(br .le. 0.0626819/0.0897209) then
             nsta = 104
           ELSEIF(br .le. 0.0654919/0.0897209) then
             nsta = 109
           ELSEIF(br .le. 0.0665919/0.0897209) then
             nsta = 111
           ELSEIF(br .le. 0.0673919/0.0897209) then
             nsta = 113
           ELSEIF(br .le. 0.0685219/0.0897209) then
             nsta = 114
           ELSEIF(br .le. 0.0687079/0.0897209) then
             nsta = 115
           ELSEIF(br .le. 0.0700279/0.0897209) then
             nsta = 120
           ELSEIF(br .le. 0.0712279/0.0897209) then
             nsta = 125
           ELSEIF(br .le. 0.0712539/0.0897209) then
             nsta = 133
           ELSEIF(br .le. 0.0726939/0.0897209) then
             nsta = 135
           ELSEIF(br .le. 0.0732339/0.0897209) then
             nsta = 139
           ELSEIF(br .le. 0.0744039/0.0897209) then
             nsta = 142
           ELSEIF(br .le. 0.0746789/0.0897209) then
             nsta = 143
           ELSEIF(br .le. 0.0749409/0.0897209) then
             nsta = 146
           ELSEIF(br .le. 0.0760909/0.0897209) then
             nsta = 152
           ELSEIF(br .le. 0.0773209/0.0897209) then
             nsta = 153
           ELSEIF(br .le. 0.0777009/0.0897209) then
             nsta = 156
           ELSEIF(br .le. 0.0780829/0.0897209) then
             nsta = 166
           ELSEIF(br .le. 0.0791729/0.0897209) then
             nsta = 167
           ELSEIF(br .le. 0.0799529/0.0897209) then
             nsta = 172
           ELSEIF(br .le. 0.0806829/0.0897209) then
             nsta = 173
           ELSEIF(br .le. 0.0809729/0.0897209) then
             nsta = 174
           ELSEIF(br .le. 0.0813029/0.0897209) then
             nsta = 175
           ELSEIF(br .le. 0.0826229/0.0897209) then
             nsta = 176
           ELSEIF(br .le. 0.0828809/0.0897209) then
             nsta = 180
           ELSEIF(br .le. 0.0834809/0.0897209) then
             nsta = 182
           ELSEIF(br .le. 0.0837809/0.0897209) then
             nsta = 183
           ELSEIF(br .le. 0.0841009/0.0897209) then
             nsta = 184
           ELSEIF(br .le. 0.0844209/0.0897209) then
             nsta = 185
           ELSEIF(br .le. 0.0848109/0.0897209) then
             nsta = 186
           ELSEIF(br .le. 0.0854009/0.0897209) then
             nsta = 189
           ELSEIF(br .le. 0.0857909/0.0897209) then
             nsta = 190
           ELSEIF(br .le. 0.0868009/0.0897209) then
             nsta = 192
           ELSEIF(br .le. 0.0872609/0.0897209) then
             nsta = 193
           ELSEIF(br .le. 0.0880309/0.0897209) then
             nsta = 195
           ELSEIF(br .le. 0.0885009/0.0897209) then
             nsta = 196
           ELSEIF(br .le. 0.0889309/0.0897209) then
             nsta = 197
           ELSEIF(br .le. 0.0894009/0.0897209) then
             nsta = 198
           ELSE
             nsta = 199
           ENDIF
         Case(38088 )
             nsta = 233
         !Case(38088)
         !      IF(br .le. 0.00023/0.58136) then
         !    nsta = 1
         !  ELSEIF(br .le. 0.01993/0.58136) then
         !    nsta = 2
         !  ELSEIF(br .le. 0.02038/0.58136) then
         !    nsta = 4
         !  ELSEIF(br .le. 0.08908/0.58136) then
         !    nsta = 8
         !  ELSEIF(br .le. 0.09678/0.58136) then
         !    nsta = 9
         !  ELSEIF(br .le. 0.09975/0.58136) then
         !    nsta = 10
         !  ELSEIF(br .le. 0.10021/0.58136) then
         !    nsta = 14
         !  ELSEIF(br .le. 0.15041/0.58136) then
         !    nsta = 16
         !  ELSEIF(br .le. 0.19821/0.58136) then
         !    nsta = 19
         !  ELSEIF(br .le. 0.20781/0.58136) then
         !    nsta = 20
         !  ELSEIF(br .le. 0.21023/0.58136) then
         !    nsta = 21
         !  ELSEIF(br .le. 0.22293/0.58136) then
         !    nsta = 24
         !  ELSEIF(br .le. 0.23613/0.58136) then
         !    nsta = 25
         !  ELSEIF(br .le. 0.30053/0.58136) then
         !    nsta = 26
         !  ELSEIF(br .le. 0.37753/0.58136) then
         !    nsta = 40
         !  ELSEIF(br .le. 0.42523/0.58136) then
         !    nsta = 46
         !  ELSEIF(br .le. 0.42543/0.58136) then
         !    nsta = 48
         !  ELSEIF(br .le. 0.43633/0.58136) then
         !    nsta = 51
         !  ELSEIF(br .le. 0.45593/0.58136) then
         !    nsta = 58
         !  ELSEIF(br .le. 0.45943/0.58136) then
         !    nsta = 62
         !  ELSEIF(br .le. 0.47253/0.58136) then
         !    nsta = 64
         !  ELSEIF(br .le. 0.47586/0.58136) then
         !    nsta = 68
         !  ELSEIF(br .le. 0.49046/0.58136) then
         !    nsta = 77
         !  ELSEIF(br .le. 0.49576/0.58136) then
         !    nsta = 84
         !  ELSEIF(br .le. 0.50276/0.58136) then
         !    nsta = 85
         !  ELSEIF(br .le. 0.51656/0.58136) then
         !    nsta = 89
         !  ELSEIF(br .le. 0.51941/0.58136) then
         !    nsta = 91
         !  ELSEIF(br .le. 0.52262/0.58136) then
         !    nsta = 93
         !  ELSEIF(br .le. 0.53072/0.58136) then
         !    nsta = 100
         !  ELSEIF(br .le. 0.53407/0.58136) then
         !    nsta = 111
         !  ELSEIF(br .le. 0.55097/0.58136) then
         !    nsta = 126
         !  ELSEIF(br .le. 0.55667/0.58136) then
         !    nsta = 132
         !  ELSEIF(br .le. 0.56127/0.58136) then
         !    nsta = 134
         !  ELSEIF(br .le. 0.56419/0.58136) then
         !    nsta = 141
         !  ELSEIF(br .le. 0.56764/0.58136) then
         !    nsta = 147
         !  ELSEIF(br .le. 0.56992/0.58136) then
         !    nsta = 155
         !  ELSEIF(br .le. 0.57262/0.58136) then
         !    nsta = 166
         !  ELSEIF(br .le. 0.57656/0.58136) then
         !    nsta = 169
         !  ELSE
         !    nsta = 182
         !  ENDIF
         Case(38089)
               IF(br .le. 0.0002/0.007314) then
             nsta = 0
           ELSEIF(br .le. 0.00045/0.007314) then
             nsta = 1
           ELSEIF(br .le. 0.000456/0.007314) then
             nsta = 3
           ELSEIF(br .le. 0.000486/0.007314) then
             nsta = 4
           ELSEIF(br .le. 0.005986/0.007314) then
             nsta = 8
           ELSEIF(br .le. 0.00604/0.007314) then
             nsta = 9
           ELSEIF(br .le. 0.00648/0.007314) then
             nsta = 10
           ELSEIF(br .le. 0.006529/0.007314) then
             nsta = 20
           ELSEIF(br .le. 0.006584/0.007314) then
             nsta = 27
           ELSEIF(br .le. 0.006734/0.007314) then
             nsta = 29
           ELSEIF(br .le. 0.006765/0.007314) then
             nsta = 37
           ELSEIF(br .le. 0.006872/0.007314) then
             nsta = 64
           ELSEIF(br .le. 0.006932/0.007314) then
             nsta = 69
           ELSEIF(br .le. 0.007152/0.007314) then
             nsta = 70
           ELSEIF(br .le. 0.007217/0.007314) then
             nsta = 76
           ELSEIF(br .le. 0.007263/0.007314) then
             nsta = 84
           ELSE
             nsta = 93
           ENDIF
         Case(39090)
               IF(br .le. 0.00279/1.21915) then
             nsta = 0
           ELSEIF(br .le. 0.00638/1.21915) then
             nsta = 1
           ELSEIF(br .le. 0.76638/1.21915) then
             nsta = 3
           ELSEIF(br .le. 0.79538/1.21915) then
             nsta = 7
           ELSEIF(br .le. 0.79958/1.21915) then
             nsta = 9
           ELSEIF(br .le. 0.80248/1.21915) then
             nsta = 12
           ELSEIF(br .le. 0.80279/1.21915) then
             nsta = 13
           ELSEIF(br .le. 0.80859/1.21915) then
             nsta = 14
           ELSEIF(br .le. 0.81419/1.21915) then
             nsta = 16
           ELSEIF(br .le. 0.81513/1.21915) then
             nsta = 20
           ELSEIF(br .le. 0.82393/1.21915) then
             nsta = 22
           ELSEIF(br .le. 0.8277/1.21915) then
             nsta = 25
           ELSEIF(br .le. 0.837/1.21915) then
             nsta = 28
           ELSEIF(br .le. 0.8455/1.21915) then
             nsta = 30
           ELSEIF(br .le. 0.84606/1.21915) then
             nsta = 31
           ELSEIF(br .le. 0.86676/1.21915) then
             nsta = 32
           ELSEIF(br .le. 0.86716/1.21915) then
             nsta = 36
           ELSEIF(br .le. 0.86758/1.21915) then
             nsta = 37
           ELSEIF(br .le. 0.86998/1.21915) then
             nsta = 38
           ELSEIF(br .le. 0.87208/1.21915) then
             nsta = 41
           ELSEIF(br .le. 0.93908/1.21915) then
             nsta = 43
           ELSEIF(br .le. 0.94988/1.21915) then
             nsta = 45
           ELSEIF(br .le. 0.95467/1.21915) then
             nsta = 48
           ELSEIF(br .le. 0.96357/1.21915) then
             nsta = 49
           ELSEIF(br .le. 0.96469/1.21915) then
             nsta = 50
           ELSEIF(br .le. 0.97359/1.21915) then
             nsta = 55
           ELSEIF(br .le. 0.97473/1.21915) then
             nsta = 56
           ELSEIF(br .le. 0.97546/1.21915) then
             nsta = 58
           ELSEIF(br .le. 0.97684/1.21915) then
             nsta = 59
           ELSEIF(br .le. 0.98464/1.21915) then
             nsta = 64
           ELSEIF(br .le. 0.99844/1.21915) then
             nsta = 65
           ELSEIF(br .le. 1.01474/1.21915) then
             nsta = 69
           ELSEIF(br .le. 1.01695/1.21915) then
             nsta = 72
           ELSEIF(br .le. 1.01865/1.21915) then
             nsta = 75
           ELSEIF(br .le. 1.03455/1.21915) then
             nsta = 76
           ELSEIF(br .le. 1.06215/1.21915) then
             nsta = 81
           ELSEIF(br .le. 1.08135/1.21915) then
             nsta = 82
           ELSEIF(br .le. 1.09325/1.21915) then
             nsta = 85
           ELSEIF(br .le. 1.10485/1.21915) then
             nsta = 87
           ELSEIF(br .le. 1.11235/1.21915) then
             nsta = 90
           ELSEIF(br .le. 1.12435/1.21915) then
             nsta = 91
           ELSEIF(br .le. 1.13335/1.21915) then
             nsta = 100
           ELSEIF(br .le. 1.13915/1.21915) then
             nsta = 102
           ELSEIF(br .le. 1.14875/1.21915) then
             nsta = 108
           ELSEIF(br .le. 1.15905/1.21915) then
             nsta = 112
           ELSEIF(br .le. 1.16605/1.21915) then
             nsta = 115
           ELSEIF(br .le. 1.17325/1.21915) then
             nsta = 118
           ELSEIF(br .le. 1.18695/1.21915) then
             nsta = 121
           ELSEIF(br .le. 1.18875/1.21915) then
             nsta = 124
           ELSEIF(br .le. 1.19825/1.21915) then
             nsta = 130
           ELSEIF(br .le. 1.20435/1.21915) then
             nsta = 133
           ELSEIF(br .le. 1.21515/1.21915) then
             nsta = 135
           ELSE
             nsta = 136
           ENDIF
         Case(40091 ) 
             nsta = 165
         Case(40092 ) 
             nsta = 135
         !Case(40092)
         !      IF(br .le. 0.0279/0.0454) then
         !    nsta = 7
         !  ELSEIF(br .le. 0.0295/0.0454) then
         !    nsta = 25
         !  ELSEIF(br .le. 0.0359/0.0454) then
         !    nsta = 32
         !  ELSEIF(br .le. 0.0378/0.0454) then
         !    nsta = 37
         !  ELSEIF(br .le. 0.0427/0.0454) then
         !    nsta = 45
         !  ELSE
         !    nsta = 47
         !  ENDIF
         Case(40093)
               IF(br .le. 0.00091/0.00785) then
             nsta = 0
           ELSEIF(br .le. 0.00182/0.00785) then
             nsta = 1
           ELSEIF(br .le. 0.00422/0.00785) then
             nsta = 7
           ELSEIF(br .le. 0.00569/0.00785) then
             nsta = 29
           ELSEIF(br .le. 0.00601/0.00785) then
             nsta = 30
           ELSEIF(br .le. 0.00648/0.00785) then
             nsta = 32
           ELSE
             nsta = 50
           ENDIF
         Case(40095)
               IF(br .le. 0.0026/0.0055) then
             nsta = 0
           ELSE
             nsta = 19
           ENDIF
         Case(41094)
               IF(br .le. 0.00137/0.82845) then
             nsta = 0
           ELSEIF(br .le. 0.01027/0.82845) then
             nsta = 1
           ELSEIF(br .le. 0.01207/0.82845) then
             nsta = 2
           ELSEIF(br .le. 0.01239/0.82845) then
             nsta = 3
           ELSEIF(br .le. 0.01399/0.82845) then
             nsta = 4
           ELSEIF(br .le. 0.01639/0.82845) then
             nsta = 7
           ELSEIF(br .le. 0.01712/0.82845) then
             nsta = 8
           ELSEIF(br .le. 0.03462/0.82845) then
             nsta = 9
           ELSEIF(br .le. 0.03662/0.82845) then
             nsta = 11
           ELSEIF(br .le. 0.03712/0.82845) then
             nsta = 12
           ELSEIF(br .le. 0.04182/0.82845) then
             nsta = 15
           ELSEIF(br .le. 0.04276/0.82845) then
             nsta = 16
           ELSEIF(br .le. 0.04566/0.82845) then
             nsta = 17
           ELSEIF(br .le. 0.04896/0.82845) then
             nsta = 20
           ELSEIF(br .le. 0.05033/0.82845) then
             nsta = 22
           ELSEIF(br .le. 0.05134/0.82845) then
             nsta = 23
           ELSEIF(br .le. 0.05154/0.82845) then
             nsta = 24
           ELSEIF(br .le. 0.052/0.82845) then
             nsta = 27
           ELSEIF(br .le. 0.05308/0.82845) then
             nsta = 31
           ELSEIF(br .le. 0.05376/0.82845) then
             nsta = 32
           ELSEIF(br .le. 0.05636/0.82845) then
             nsta = 33
           ELSEIF(br .le. 0.05796/0.82845) then
             nsta = 35
           ELSEIF(br .le. 0.06126/0.82845) then
             nsta = 40
           ELSEIF(br .le. 0.06416/0.82845) then
             nsta = 41
           ELSEIF(br .le. 0.06966/0.82845) then
             nsta = 43
           ELSEIF(br .le. 0.07216/0.82845) then
             nsta = 44
           ELSEIF(br .le. 0.07666/0.82845) then
             nsta = 45
           ELSEIF(br .le. 0.07802/0.82845) then
             nsta = 46
           ELSEIF(br .le. 0.09632/0.82845) then
             nsta = 47
           ELSEIF(br .le. 0.09982/0.82845) then
             nsta = 48
           ELSEIF(br .le. 0.10004/0.82845) then
             nsta = 49
           ELSEIF(br .le. 0.10194/0.82845) then
             nsta = 50
           ELSEIF(br .le. 0.10734/0.82845) then
             nsta = 54
           ELSEIF(br .le. 0.10814/0.82845) then
             nsta = 57
           ELSEIF(br .le. 0.10864/0.82845) then
             nsta = 58
           ELSEIF(br .le. 0.11014/0.82845) then
             nsta = 59
           ELSEIF(br .le. 0.11204/0.82845) then
             nsta = 60
           ELSEIF(br .le. 0.11464/0.82845) then
             nsta = 63
           ELSEIF(br .le. 0.11834/0.82845) then
             nsta = 65
           ELSEIF(br .le. 0.12244/0.82845) then
             nsta = 66
           ELSEIF(br .le. 0.13044/0.82845) then
             nsta = 67
           ELSEIF(br .le. 0.13182/0.82845) then
             nsta = 68
           ELSEIF(br .le. 0.13552/0.82845) then
             nsta = 69
           ELSEIF(br .le. 0.1361/0.82845) then
             nsta = 71
           ELSEIF(br .le. 0.1388/0.82845) then
             nsta = 72
           ELSEIF(br .le. 0.1441/0.82845) then
             nsta = 73
           ELSEIF(br .le. 0.1482/0.82845) then
             nsta = 74
           ELSEIF(br .le. 0.1687/0.82845) then
             nsta = 75
           ELSEIF(br .le. 0.17027/0.82845) then
             nsta = 76
           ELSEIF(br .le. 0.17081/0.82845) then
             nsta = 77
           ELSEIF(br .le. 0.17611/0.82845) then
             nsta = 78
           ELSEIF(br .le. 0.17871/0.82845) then
             nsta = 79
           ELSEIF(br .le. 0.17929/0.82845) then
             nsta = 80
           ELSEIF(br .le. 0.18101/0.82845) then
             nsta = 81
           ELSEIF(br .le. 0.18223/0.82845) then
             nsta = 82
           ELSEIF(br .le. 0.18723/0.82845) then
             nsta = 83
           ELSEIF(br .le. 0.19113/0.82845) then
             nsta = 84
           ELSEIF(br .le. 0.19843/0.82845) then
             nsta = 85
           ELSEIF(br .le. 0.20663/0.82845) then
             nsta = 86
           ELSEIF(br .le. 0.20837/0.82845) then
             nsta = 87
           ELSEIF(br .le. 0.21467/0.82845) then
             nsta = 90
           ELSEIF(br .le. 0.21777/0.82845) then
             nsta = 91
           ELSEIF(br .le. 0.21997/0.82845) then
             nsta = 93
           ELSEIF(br .le. 0.22497/0.82845) then
             nsta = 94
           ELSEIF(br .le. 0.22622/0.82845) then
             nsta = 95
           ELSEIF(br .le. 0.22734/0.82845) then
             nsta = 96
           ELSEIF(br .le. 0.22948/0.82845) then
             nsta = 97
           ELSEIF(br .le. 0.23748/0.82845) then
             nsta = 98
           ELSEIF(br .le. 0.23868/0.82845) then
             nsta = 100
           ELSEIF(br .le. 0.2389/0.82845) then
             nsta = 101
           ELSEIF(br .le. 0.23974/0.82845) then
             nsta = 102
           ELSEIF(br .le. 0.2417/0.82845) then
             nsta = 103
           ELSEIF(br .le. 0.2489/0.82845) then
             nsta = 104
           ELSEIF(br .le. 0.2603/0.82845) then
             nsta = 105
           ELSEIF(br .le. 0.2675/0.82845) then
             nsta = 106
           ELSEIF(br .le. 0.26788/0.82845) then
             nsta = 107
           ELSEIF(br .le. 0.26894/0.82845) then
             nsta = 108
           ELSEIF(br .le. 0.2696/0.82845) then
             nsta = 109
           ELSEIF(br .le. 0.273/0.82845) then
             nsta = 110
           ELSEIF(br .le. 0.27358/0.82845) then
             nsta = 111
           ELSEIF(br .le. 0.29678/0.82845) then
             nsta = 112
           ELSEIF(br .le. 0.29732/0.82845) then
             nsta = 113
           ELSEIF(br .le. 0.30032/0.82845) then
             nsta = 114
           ELSEIF(br .le. 0.30074/0.82845) then
             nsta = 115
           ELSEIF(br .le. 0.31094/0.82845) then
             nsta = 116
           ELSEIF(br .le. 0.31434/0.82845) then
             nsta = 117
           ELSEIF(br .le. 0.3152/0.82845) then
             nsta = 118
           ELSEIF(br .le. 0.3174/0.82845) then
             nsta = 119
           ELSEIF(br .le. 0.3178/0.82845) then
             nsta = 121
           ELSEIF(br .le. 0.3236/0.82845) then
             nsta = 122
           ELSEIF(br .le. 0.32436/0.82845) then
             nsta = 123
           ELSEIF(br .le. 0.32516/0.82845) then
             nsta = 124
           ELSEIF(br .le. 0.32628/0.82845) then
             nsta = 125
           ELSEIF(br .le. 0.32958/0.82845) then
             nsta = 126
           ELSEIF(br .le. 0.33094/0.82845) then
             nsta = 127
           ELSEIF(br .le. 0.33874/0.82845) then
             nsta = 128
           ELSEIF(br .le. 0.33924/0.82845) then
             nsta = 129
           ELSEIF(br .le. 0.34434/0.82845) then
             nsta = 131
           ELSEIF(br .le. 0.34724/0.82845) then
             nsta = 132
           ELSEIF(br .le. 0.34834/0.82845) then
             nsta = 133
           ELSEIF(br .le. 0.34854/0.82845) then
             nsta = 133
           ELSEIF(br .le. 0.35124/0.82845) then
             nsta = 134
           ELSEIF(br .le. 0.35199/0.82845) then
             nsta = 135
           ELSEIF(br .le. 0.35979/0.82845) then
             nsta = 136
           ELSEIF(br .le. 0.36172/0.82845) then
             nsta = 137
           ELSEIF(br .le. 0.36245/0.82845) then
             nsta = 138
           ELSEIF(br .le. 0.36269/0.82845) then
             nsta = 139
           ELSEIF(br .le. 0.36317/0.82845) then
             nsta = 140
           ELSEIF(br .le. 0.36472/0.82845) then
             nsta = 141
           ELSEIF(br .le. 0.36625/0.82845) then
             nsta = 142
           ELSEIF(br .le. 0.36719/0.82845) then
             nsta = 143
           ELSEIF(br .le. 0.36761/0.82845) then
             nsta = 144
           ELSEIF(br .le. 0.36815/0.82845) then
             nsta = 145
           ELSEIF(br .le. 0.37385/0.82845) then
             nsta = 146
           ELSEIF(br .le. 0.37425/0.82845) then
             nsta = 147
           ELSEIF(br .le. 0.37531/0.82845) then
             nsta = 148
           ELSEIF(br .le. 0.3763/0.82845) then
             nsta = 149
           ELSEIF(br .le. 0.3834/0.82845) then
             nsta = 151
           ELSEIF(br .le. 0.38396/0.82845) then
             nsta = 152
           ELSEIF(br .le. 0.38436/0.82845) then
             nsta = 153
           ELSEIF(br .le. 0.38886/0.82845) then
             nsta = 155
           ELSEIF(br .le. 0.39276/0.82845) then
             nsta = 156
           ELSEIF(br .le. 0.39656/0.82845) then
             nsta = 157
           ELSEIF(br .le. 0.41186/0.82845) then
             nsta = 158
           ELSEIF(br .le. 0.41317/0.82845) then
             nsta = 159
           ELSEIF(br .le. 0.4142/0.82845) then
             nsta = 160
           ELSEIF(br .le. 0.4194/0.82845) then
             nsta = 161
           ELSEIF(br .le. 0.41996/0.82845) then
             nsta = 162
           ELSEIF(br .le. 0.4206/0.82845) then
             nsta = 163
           ELSEIF(br .le. 0.4265/0.82845) then
             nsta = 164
           ELSEIF(br .le. 0.433/0.82845) then
             nsta = 165
           ELSEIF(br .le. 0.4358/0.82845) then
             nsta = 166
           ELSEIF(br .le. 0.43653/0.82845) then
             nsta = 167
           ELSEIF(br .le. 0.43712/0.82845) then
             nsta = 168
           ELSEIF(br .le. 0.44182/0.82845) then
             nsta = 169
           ELSEIF(br .le. 0.44672/0.82845) then
             nsta = 170
           ELSEIF(br .le. 0.44865/0.82845) then
             nsta = 171
           ELSEIF(br .le. 0.45325/0.82845) then
             nsta = 172
           ELSEIF(br .le. 0.45795/0.82845) then
             nsta = 173
           ELSEIF(br .le. 0.45837/0.82845) then
             nsta = 174
           ELSEIF(br .le. 0.45861/0.82845) then
             nsta = 175
           ELSEIF(br .le. 0.45899/0.82845) then
             nsta = 176
           ELSEIF(br .le. 0.46389/0.82845) then
             nsta = 178
           ELSEIF(br .le. 0.46719/0.82845) then
             nsta = 179
           ELSEIF(br .le. 0.46818/0.82845) then
             nsta = 180
           ELSEIF(br .le. 0.47398/0.82845) then
             nsta = 181
           ELSEIF(br .le. 0.47488/0.82845) then
             nsta = 182
           ELSEIF(br .le. 0.47868/0.82845) then
             nsta = 183
           ELSEIF(br .le. 0.48158/0.82845) then
             nsta = 185
           ELSEIF(br .le. 0.48718/0.82845) then
             nsta = 186
           ELSEIF(br .le. 0.48828/0.82845) then
             nsta = 187
           ELSEIF(br .le. 0.48892/0.82845) then
             nsta = 188
           ELSEIF(br .le. 0.48908/0.82845) then
             nsta = 189
           ELSEIF(br .le. 0.49238/0.82845) then
             nsta = 190
           ELSEIF(br .le. 0.49518/0.82845) then
             nsta = 191
           ELSEIF(br .le. 0.49818/0.82845) then
             nsta = 192
           ELSEIF(br .le. 0.50088/0.82845) then
             nsta = 193
           ELSEIF(br .le. 0.50147/0.82845) then
             nsta = 194
           ELSEIF(br .le. 0.50577/0.82845) then
             nsta = 197
           ELSEIF(br .le. 0.50617/0.82845) then
             nsta = 198
           ELSEIF(br .le. 0.50787/0.82845) then
             nsta = 201
           ELSEIF(br .le. 0.50916/0.82845) then
             nsta = 202
           ELSEIF(br .le. 0.51356/0.82845) then
             nsta = 203
           ELSEIF(br .le. 0.51552/0.82845) then
             nsta = 204
           ELSEIF(br .le. 0.51842/0.82845) then
             nsta = 205
           ELSEIF(br .le. 0.51938/0.82845) then
             nsta = 206
           ELSEIF(br .le. 0.52018/0.82845) then
             nsta = 207
           ELSEIF(br .le. 0.52171/0.82845) then
             nsta = 208
           ELSEIF(br .le. 0.52441/0.82845) then
             nsta = 209
           ELSEIF(br .le. 0.52535/0.82845) then
             nsta = 210
           ELSEIF(br .le. 0.52965/0.82845) then
             nsta = 211
           ELSEIF(br .le. 0.53235/0.82845) then
             nsta = 212
           ELSEIF(br .le. 0.53455/0.82845) then
             nsta = 213
           ELSEIF(br .le. 0.53945/0.82845) then
             nsta = 214
           ELSEIF(br .le. 0.54082/0.82845) then
             nsta = 216
           ELSEIF(br .le. 0.54245/0.82845) then
             nsta = 217
           ELSEIF(br .le. 0.54425/0.82845) then
             nsta = 218
           ELSEIF(br .le. 0.54785/0.82845) then
             nsta = 219
           ELSEIF(br .le. 0.54831/0.82845) then
             nsta = 220
           ELSEIF(br .le. 0.54997/0.82845) then
             nsta = 221
           ELSEIF(br .le. 0.55197/0.82845) then
             nsta = 222
           ELSEIF(br .le. 0.55253/0.82845) then
             nsta = 223
           ELSEIF(br .le. 0.55363/0.82845) then
             nsta = 225
           ELSEIF(br .le. 0.55653/0.82845) then
             nsta = 226
           ELSEIF(br .le. 0.55923/0.82845) then
             nsta = 227
           ELSEIF(br .le. 0.56119/0.82845) then
             nsta = 228
           ELSEIF(br .le. 0.56227/0.82845) then
             nsta = 229
           ELSEIF(br .le. 0.56367/0.82845) then
             nsta = 230
           ELSEIF(br .le. 0.56537/0.82845) then
             nsta = 231
           ELSEIF(br .le. 0.56688/0.82845) then
             nsta = 233
           ELSEIF(br .le. 0.56968/0.82845) then
             nsta = 234
           ELSEIF(br .le. 0.57178/0.82845) then
             nsta = 235
           ELSEIF(br .le. 0.57344/0.82845) then
             nsta = 236
           ELSEIF(br .le. 0.57974/0.82845) then
             nsta = 237
           ELSEIF(br .le. 0.58167/0.82845) then
             nsta = 238
           ELSEIF(br .le. 0.58427/0.82845) then
             nsta = 239
           ELSEIF(br .le. 0.58697/0.82845) then
             nsta = 240
           ELSEIF(br .le. 0.58847/0.82845) then
             nsta = 241
           ELSEIF(br .le. 0.59057/0.82845) then
             nsta = 242
           ELSEIF(br .le. 0.59103/0.82845) then
             nsta = 243
           ELSEIF(br .le. 0.59179/0.82845) then
             nsta = 244
           ELSEIF(br .le. 0.5925/0.82845) then
             nsta = 245
           ELSEIF(br .le. 0.59314/0.82845) then
             nsta = 246
           ELSEIF(br .le. 0.59356/0.82845) then
             nsta = 247
           ELSEIF(br .le. 0.59539/0.82845) then
             nsta = 249
           ELSEIF(br .le. 0.59661/0.82845) then
             nsta = 250
           ELSEIF(br .le. 0.59737/0.82845) then
             nsta = 251
           ELSEIF(br .le. 0.5989/0.82845) then
             nsta = 252
           ELSEIF(br .le. 0.6044/0.82845) then
             nsta = 253
           ELSEIF(br .le. 0.6077/0.82845) then
             nsta = 254
           ELSEIF(br .le. 0.60836/0.82845) then
             nsta = 255
           ELSEIF(br .le. 0.61166/0.82845) then
             nsta = 256
           ELSEIF(br .le. 0.61296/0.82845) then
             nsta = 257
           ELSEIF(br .le. 0.61536/0.82845) then
             nsta = 258
           ELSEIF(br .le. 0.61836/0.82845) then
             nsta = 259
           ELSEIF(br .le. 0.61954/0.82845) then
             nsta = 260
           ELSEIF(br .le. 0.62147/0.82845) then
             nsta = 261
           ELSEIF(br .le. 0.62317/0.82845) then
             nsta = 262
           ELSEIF(br .le. 0.62647/0.82845) then
             nsta = 263
           ELSEIF(br .le. 0.62887/0.82845) then
             nsta = 264
           ELSEIF(br .le. 0.63147/0.82845) then
             nsta = 265
           ELSEIF(br .le. 0.63313/0.82845) then
             nsta = 266
           ELSEIF(br .le. 0.63693/0.82845) then
             nsta = 268
           ELSEIF(br .le. 0.63913/0.82845) then
             nsta = 269
           ELSEIF(br .le. 0.6412/0.82845) then
             nsta = 270
           ELSEIF(br .le. 0.6424/0.82845) then
             nsta = 271
           ELSEIF(br .le. 0.6463/0.82845) then
             nsta = 272
           ELSEIF(br .le. 0.6514/0.82845) then
             nsta = 273
           ELSEIF(br .le. 0.65317/0.82845) then
             nsta = 274
           ELSEIF(br .le. 0.65797/0.82845) then
             nsta = 275
           ELSEIF(br .le. 0.65863/0.82845) then
             nsta = 276
           ELSEIF(br .le. 0.66123/0.82845) then
             nsta = 277
           ELSEIF(br .le. 0.66205/0.82845) then
             nsta = 278
           ELSEIF(br .le. 0.66385/0.82845) then
             nsta = 279
           ELSEIF(br .le. 0.66444/0.82845) then
             nsta = 280
           ELSEIF(br .le. 0.66644/0.82845) then
             nsta = 281
           ELSEIF(br .le. 0.66726/0.82845) then
             nsta = 282
           ELSEIF(br .le. 0.66863/0.82845) then
             nsta = 283
           ELSEIF(br .le. 0.66949/0.82845) then
             nsta = 284
           ELSEIF(br .le. 0.67085/0.82845) then
             nsta = 285
           ELSEIF(br .le. 0.67275/0.82845) then
             nsta = 286
           ELSEIF(br .le. 0.67313/0.82845) then
             nsta = 287
           ELSEIF(br .le. 0.67423/0.82845) then
             nsta = 288
           ELSEIF(br .le. 0.6755/0.82845) then
             nsta = 289
           ELSEIF(br .le. 0.67686/0.82845) then
             nsta = 290
           ELSEIF(br .le. 0.67852/0.82845) then
             nsta = 291
           ELSEIF(br .le. 0.67884/0.82845) then
             nsta = 292
           ELSEIF(br .le. 0.68334/0.82845) then
             nsta = 293
           ELSEIF(br .le. 0.6845/0.82845) then
             nsta = 294
           ELSEIF(br .le. 0.6865/0.82845) then
             nsta = 295
           ELSEIF(br .le. 0.68801/0.82845) then
             nsta = 296
           ELSEIF(br .le. 0.69131/0.82845) then
             nsta = 298
           ELSEIF(br .le. 0.69341/0.82845) then
             nsta = 299
           ELSEIF(br .le. 0.69425/0.82845) then
             nsta = 300
           ELSEIF(br .le. 0.69587/0.82845) then
             nsta = 301
           ELSEIF(br .le. 0.69917/0.82845) then
             nsta = 302
           ELSEIF(br .le. 0.70025/0.82845) then
             nsta = 303
           ELSEIF(br .le. 0.70084/0.82845) then
             nsta = 304
           ELSEIF(br .le. 0.70214/0.82845) then
             nsta = 305
           ELSEIF(br .le. 0.70282/0.82845) then
             nsta = 306
           ELSEIF(br .le. 0.7047/0.82845) then
             nsta = 307
           ELSEIF(br .le. 0.7075/0.82845) then
             nsta = 308
           ELSEIF(br .le. 0.7103/0.82845) then
             nsta = 309
           ELSEIF(br .le. 0.7126/0.82845) then
             nsta = 310
           ELSEIF(br .le. 0.7149/0.82845) then
             nsta = 311
           ELSEIF(br .le. 0.7171/0.82845) then
             nsta = 312
           ELSEIF(br .le. 0.71792/0.82845) then
             nsta = 313
           ELSEIF(br .le. 0.72062/0.82845) then
             nsta = 314
           ELSEIF(br .le. 0.72133/0.82845) then
             nsta = 315
           ELSEIF(br .le. 0.72191/0.82845) then
             nsta = 316
           ELSEIF(br .le. 0.72351/0.82845) then
             nsta = 317
           ELSEIF(br .le. 0.72482/0.82845) then
             nsta = 318
           ELSEIF(br .le. 0.72604/0.82845) then
             nsta = 319
           ELSEIF(br .le. 0.72636/0.82845) then
             nsta = 320
           ELSEIF(br .le. 0.72906/0.82845) then
             nsta = 321
           ELSEIF(br .le. 0.73069/0.82845) then
             nsta = 322
           ELSEIF(br .le. 0.73246/0.82845) then
             nsta = 323
           ELSEIF(br .le. 0.73456/0.82845) then
             nsta = 324
           ELSEIF(br .le. 0.73654/0.82845) then
             nsta = 325
           ELSEIF(br .le. 0.73819/0.82845) then
             nsta = 326
           ELSEIF(br .le. 0.73899/0.82845) then
             nsta = 327
           ELSEIF(br .le. 0.73991/0.82845) then
             nsta = 328
           ELSEIF(br .le. 0.74184/0.82845) then
             nsta = 329
           ELSEIF(br .le. 0.74594/0.82845) then
             nsta = 330
           ELSEIF(br .le. 0.74814/0.82845) then
             nsta = 331
           ELSEIF(br .le. 0.74952/0.82845) then
             nsta = 332
           ELSEIF(br .le. 0.75068/0.82845) then
             nsta = 333
           ELSEIF(br .le. 0.75206/0.82845) then
             nsta = 334
           ELSEIF(br .le. 0.75389/0.82845) then
             nsta = 335
           ELSEIF(br .le. 0.75479/0.82845) then
             nsta = 336
           ELSEIF(br .le. 0.75729/0.82845) then
             nsta = 337
           ELSEIF(br .le. 0.76029/0.82845) then
             nsta = 338
           ELSEIF(br .le. 0.76088/0.82845) then
             nsta = 339
           ELSEIF(br .le. 0.76172/0.82845) then
             nsta = 340
           ELSEIF(br .le. 0.76248/0.82845) then
             nsta = 341
           ELSEIF(br .le. 0.76342/0.82845) then
             nsta = 342
           ELSEIF(br .le. 0.76642/0.82845) then
             nsta = 343
           ELSEIF(br .le. 0.76842/0.82845) then
             nsta = 344
           ELSEIF(br .le. 0.77025/0.82845) then
             nsta = 345
           ELSEIF(br .le. 0.77171/0.82845) then
             nsta = 346
           ELSEIF(br .le. 0.77291/0.82845) then
             nsta = 347
           ELSEIF(br .le. 0.77403/0.82845) then
             nsta = 348
           ELSEIF(br .le. 0.77549/0.82845) then
             nsta = 349
           ELSEIF(br .le. 0.77889/0.82845) then
             nsta = 350
           ELSEIF(br .le. 0.78109/0.82845) then
             nsta = 351
           ELSEIF(br .le. 0.78159/0.82845) then
             nsta = 352
           ELSEIF(br .le. 0.78359/0.82845) then
             nsta = 353
           ELSEIF(br .le. 0.78629/0.82845) then
             nsta = 354
           ELSEIF(br .le. 0.78909/0.82845) then
             nsta = 355
           ELSEIF(br .le. 0.79051/0.82845) then
             nsta = 356
           ELSEIF(br .le. 0.79281/0.82845) then
             nsta = 357
           ELSEIF(br .le. 0.79418/0.82845) then
             nsta = 358
           ELSEIF(br .le. 0.79608/0.82845) then
             nsta = 359
           ELSEIF(br .le. 0.79798/0.82845) then
             nsta = 360
           ELSEIF(br .le. 0.80078/0.82845) then
             nsta = 361
           ELSEIF(br .le. 0.80205/0.82845) then
             nsta = 362
           ELSEIF(br .le. 0.80271/0.82845) then
             nsta = 363
           ELSEIF(br .le. 0.80327/0.82845) then
             nsta = 364
           ELSEIF(br .le. 0.80402/0.82845) then
             nsta = 365
           ELSEIF(br .le. 0.80456/0.82845) then
             nsta = 366
           ELSEIF(br .le. 0.8051/0.82845) then
             nsta = 367
           ELSEIF(br .le. 0.8069/0.82845) then
             nsta = 368
           ELSEIF(br .le. 0.80838/0.82845) then
             nsta = 369
           ELSEIF(br .le. 0.80918/0.82845) then
             nsta = 370
           ELSEIF(br .le. 0.81019/0.82845) then
             nsta = 371
           ELSEIF(br .le. 0.81439/0.82845) then
             nsta = 372
           ELSEIF(br .le. 0.81849/0.82845) then
             nsta = 373
           ELSEIF(br .le. 0.8198/0.82845) then
             nsta = 374
           ELSEIF(br .le. 0.82102/0.82845) then
             nsta = 375
           ELSEIF(br .le. 0.82239/0.82845) then
             nsta = 376
           ELSEIF(br .le. 0.82499/0.82845) then
             nsta = 377
           ELSEIF(br .le. 0.82652/0.82845) then
             nsta = 378
           ELSE
             nsta = 379
           ENDIF
         Case(42101)
               IF(br .le. 0.00008/0.004833) then
             nsta = 0
           ELSEIF(br .le. 0.00118/0.004833) then
             nsta = 9
           ELSEIF(br .le. 0.001279/0.004833) then
             nsta = 19
           ELSEIF(br .le. 0.001353/0.004833) then
             nsta = 29
           ELSEIF(br .le. 0.001409/0.004833) then
             nsta = 31
           ELSEIF(br .le. 0.001829/0.004833) then
             nsta = 36
           ELSEIF(br .le. 0.001888/0.004833) then
             nsta = 39
           ELSEIF(br .le. 0.002558/0.004833) then
             nsta = 45
           ELSEIF(br .le. 0.002798/0.004833) then
             nsta = 48
           ELSEIF(br .le. 0.003068/0.004833) then
             nsta = 49
           ELSEIF(br .le. 0.003178/0.004833) then
             nsta = 50
           ELSEIF(br .le. 0.003298/0.004833) then
             nsta = 51
           ELSEIF(br .le. 0.003448/0.004833) then
             nsta = 52
           ELSEIF(br .le. 0.003648/0.004833) then
             nsta = 54
           ELSEIF(br .le. 0.003838/0.004833) then
             nsta = 55
           ELSEIF(br .le. 0.004148/0.004833) then
             nsta = 56
           ELSEIF(br .le. 0.004348/0.004833) then
             nsta = 57
           ELSEIF(br .le. 0.004407/0.004833) then
             nsta = 58
           ELSEIF(br .le. 0.004657/0.004833) then
             nsta = 59
           ELSEIF(br .le. 0.004713/0.004833) then
             nsta = 60
           ELSE
             nsta = 63
           ENDIF
         Case(42093)
               IF(br .le. 0.00044/0.00464) then
             nsta = 0
           ELSEIF(br .le. 0.00364/0.00464) then
             nsta = 1
           ELSE
             nsta = 4
           ENDIF
         Case(42095)
               IF(br .le. 0.0027/0.00386) then
             nsta = 2
           ELSEIF(br .le. 0.00305/0.00386) then
             nsta = 8
           ELSEIF(br .le. 0.00326/0.00386) then
             nsta = 13
           ELSEIF(br .le. 0.00333/0.00386) then
             nsta = 36
           ELSEIF(br .le. 0.00364/0.00386) then
             nsta = 46
           ELSE
             nsta = 54
           ENDIF
         Case(42096)
               IF(br .le. 0.0189/0.36575) then
             nsta = 1
           ELSEIF(br .le. 0.0231/0.36575) then
             nsta = 4
           ELSEIF(br .le. 0.0495/0.36575) then
             nsta = 5
           ELSEIF(br .le. 0.05035/0.36575) then
             nsta = 8
           ELSEIF(br .le. 0.05103/0.36575) then
             nsta = 9
           ELSEIF(br .le. 0.05253/0.36575) then
             nsta = 10
           ELSEIF(br .le. 0.15853/0.36575) then
             nsta = 11
           ELSEIF(br .le. 0.16103/0.36575) then
             nsta = 13
           ELSEIF(br .le. 0.1625/0.36575) then
             nsta = 14
           ELSEIF(br .le. 0.1755/0.36575) then
             nsta = 17
           ELSEIF(br .le. 0.1791/0.36575) then
             nsta = 19
           ELSEIF(br .le. 0.1818/0.36575) then
             nsta = 20
           ELSEIF(br .le. 0.1847/0.36575) then
             nsta = 23
           ELSEIF(br .le. 0.1901/0.36575) then
             nsta = 24
           ELSEIF(br .le. 0.19088/0.36575) then
             nsta = 27
           ELSEIF(br .le. 0.19328/0.36575) then
             nsta = 32
           ELSEIF(br .le. 0.21678/0.36575) then
             nsta = 33
           ELSEIF(br .le. 0.21802/0.36575) then
             nsta = 36
           ELSEIF(br .le. 0.22042/0.36575) then
             nsta = 38
           ELSEIF(br .le. 0.22176/0.36575) then
             nsta = 40
           ELSEIF(br .le. 0.23796/0.36575) then
             nsta = 43
           ELSEIF(br .le. 0.23906/0.36575) then
             nsta = 44
           ELSEIF(br .le. 0.24526/0.36575) then
             nsta = 46
           ELSEIF(br .le. 0.24926/0.36575) then
             nsta = 48
           ELSEIF(br .le. 0.25125/0.36575) then
             nsta = 50
           ELSEIF(br .le. 0.25235/0.36575) then
             nsta = 51
           ELSEIF(br .le. 0.26725/0.36575) then
             nsta = 52
           ELSEIF(br .le. 0.27085/0.36575) then
             nsta = 57
           ELSEIF(br .le. 0.28045/0.36575) then
             nsta = 60
           ELSEIF(br .le. 0.29355/0.36575) then
             nsta = 65
           ELSEIF(br .le. 0.34155/0.36575) then
             nsta = 71
           ELSE
             nsta = 78
           ENDIF
         Case(42097)
               IF(br .le. 0.0007/0.0371) then
             nsta = 8
           ELSEIF(br .le. 0.0015/0.0371) then
             nsta = 9
           ELSEIF(br .le. 0.0027/0.0371) then
             nsta = 36
           ELSEIF(br .le. 0.0048/0.0371) then
             nsta = 49
           ELSEIF(br .le. 0.0055/0.0371) then
             nsta = 75
           ELSEIF(br .le. 0.0155/0.0371) then
             nsta = 83
           ELSEIF(br .le. 0.0158/0.0371) then
             nsta = 86
           ELSEIF(br .le. 0.0264/0.0371) then
             nsta = 97
           ELSEIF(br .le. 0.0272/0.0371) then
             nsta = 99
           ELSEIF(br .le. 0.0278/0.0371) then
             nsta = 103
           ELSEIF(br .le. 0.0286/0.0371) then
             nsta = 104
           ELSEIF(br .le. 0.0333/0.0371) then
             nsta = 107
           ELSEIF(br .le. 0.0353/0.0371) then
             nsta = 117
           ELSE
             nsta = 130
           ENDIF
         Case(42098)
               IF(br .le. 0.00007/0.04744) then
             nsta = 1
           ELSEIF(br .le. 0.00037/0.04744) then
             nsta = 2
           ELSEIF(br .le. 0.00064/0.04744) then
             nsta = 3
           ELSEIF(br .le. 0.00107/0.04744) then
             nsta = 4
           ELSEIF(br .le. 0.00153/0.04744) then
             nsta = 5
           ELSEIF(br .le. 0.00168/0.04744) then
             nsta = 7
           ELSEIF(br .le. 0.02868/0.04744) then
             nsta = 9
           ELSEIF(br .le. 0.02918/0.04744) then
             nsta = 11
           ELSEIF(br .le. 0.03008/0.04744) then
             nsta = 12
           ELSEIF(br .le. 0.0302/0.04744) then
             nsta = 14
           ELSEIF(br .le. 0.03035/0.04744) then
             nsta = 16
           ELSEIF(br .le. 0.03195/0.04744) then
             nsta = 22
           ELSEIF(br .le. 0.03206/0.04744) then
             nsta = 23
           ELSEIF(br .le. 0.03256/0.04744) then
             nsta = 28
           ELSEIF(br .le. 0.03362/0.04744) then
             nsta = 30
           ELSEIF(br .le. 0.03384/0.04744) then
             nsta = 34
           ELSEIF(br .le. 0.03494/0.04744) then
             nsta = 37
           ELSEIF(br .le. 0.03834/0.04744) then
             nsta = 40
           ELSEIF(br .le. 0.04104/0.04744) then
             nsta = 50
           ELSEIF(br .le. 0.04131/0.04744) then
             nsta = 51
           ELSEIF(br .le. 0.04172/0.04744) then
             nsta = 55
           ELSEIF(br .le. 0.04187/0.04744) then
             nsta = 56
           ELSEIF(br .le. 0.04241/0.04744) then
             nsta = 60
           ELSEIF(br .le. 0.04331/0.04744) then
             nsta = 61
           ELSEIF(br .le. 0.04611/0.04744) then
             nsta = 64
           ELSEIF(br .le. 0.04676/0.04744) then
             nsta = 66
           ELSE
             nsta = 68
           ENDIF
         Case(42099)
               IF(br .le. 0.0012/0.00684) then
             nsta = 0
           ELSEIF(br .le. 0.0023/0.00684) then
             nsta = 3
           ELSEIF(br .le. 0.00312/0.00684) then
             nsta = 15
           ELSEIF(br .le. 0.00399/0.00684) then
             nsta = 19
           ELSEIF(br .le. 0.00467/0.00684) then
             nsta = 44
           ELSEIF(br .le. 0.00607/0.00684) then
             nsta = 74
           ELSE
             nsta = 101
           ENDIF
         Case(44100)
               IF(br .le. 0.00009/0.38709) then
             nsta = 0
           ELSEIF(br .le. 0.00389/0.38709) then
             nsta = 1
           ELSEIF(br .le. 0.00849/0.38709) then
             nsta = 3
           ELSEIF(br .le. 0.01529/0.38709) then
             nsta = 4
           ELSEIF(br .le. 0.02029/0.38709) then
             nsta = 6
           ELSEIF(br .le. 0.15229/0.38709) then
             nsta = 7
           ELSEIF(br .le. 0.16339/0.38709) then
             nsta = 9
           ELSEIF(br .le. 0.16819/0.38709) then
             nsta = 12
           ELSEIF(br .le. 0.18039/0.38709) then
             nsta = 14
           ELSEIF(br .le. 0.18589/0.38709) then
             nsta = 21
           ELSEIF(br .le. 0.19229/0.38709) then
             nsta = 25
           ELSEIF(br .le. 0.19759/0.38709) then
             nsta = 26
           ELSEIF(br .le. 0.19999/0.38709) then
             nsta = 27
           ELSEIF(br .le. 0.20069/0.38709) then
             nsta = 30
           ELSEIF(br .le. 0.20219/0.38709) then
             nsta = 31
           ELSEIF(br .le. 0.22019/0.38709) then
             nsta = 32
           ELSEIF(br .le. 0.22329/0.38709) then
             nsta = 34
           ELSEIF(br .le. 0.23049/0.38709) then
             nsta = 40
           ELSEIF(br .le. 0.23249/0.38709) then
             nsta = 44
           ELSEIF(br .le. 0.23559/0.38709) then
             nsta = 45
           ELSEIF(br .le. 0.23759/0.38709) then
             nsta = 50
           ELSEIF(br .le. 0.23939/0.38709) then
             nsta = 53
           ELSEIF(br .le. 0.24639/0.38709) then
             nsta = 55
           ELSEIF(br .le. 0.25339/0.38709) then
             nsta = 60
           ELSEIF(br .le. 0.26329/0.38709) then
             nsta = 66
           ELSEIF(br .le. 0.26789/0.38709) then
             nsta = 68
           ELSEIF(br .le. 0.27069/0.38709) then
             nsta = 75
           ELSEIF(br .le. 0.27439/0.38709) then
             nsta = 77
           ELSEIF(br .le. 0.28539/0.38709) then
             nsta = 86
           ELSEIF(br .le. 0.29739/0.38709) then
             nsta = 89
           ELSEIF(br .le. 0.32139/0.38709) then
             nsta = 90
           ELSEIF(br .le. 0.33249/0.38709) then
             nsta = 91
           ELSEIF(br .le. 0.33899/0.38709) then
             nsta = 94
           ELSEIF(br .le. 0.35119/0.38709) then
             nsta = 99
           ELSEIF(br .le. 0.35989/0.38709) then
             nsta = 111
           ELSEIF(br .le. 0.37289/0.38709) then
             nsta = 112
           ELSEIF(br .le. 0.37599/0.38709) then
             nsta = 114
           ELSEIF(br .le. 0.37749/0.38709) then
             nsta = 115
           ELSEIF(br .le. 0.38029/0.38709) then
             nsta = 118
           ELSEIF(br .le. 0.38139/0.38709) then
             nsta = 119
           ELSEIF(br .le. 0.38289/0.38709) then
             nsta = 121
           ELSEIF(br .le. 0.38469/0.38709) then
             nsta = 129
           ELSE
             nsta = 134
           ENDIF
         Case(44101)
               IF(br .le. 0.0033/0.08662) then
             nsta = 1
           ELSEIF(br .le. 0.0213/0.08662) then
             nsta = 7
           ELSEIF(br .le. 0.02182/0.08662) then
             nsta = 15
           ELSEIF(br .le. 0.02412/0.08662) then
             nsta = 24
           ELSEIF(br .le. 0.03482/0.08662) then
             nsta = 41
           ELSEIF(br .le. 0.03582/0.08662) then
             nsta = 57
           ELSEIF(br .le. 0.04212/0.08662) then
             nsta = 60
           ELSEIF(br .le. 0.04342/0.08662) then
             nsta = 63
           ELSEIF(br .le. 0.04442/0.08662) then
             nsta = 64
           ELSEIF(br .le. 0.04832/0.08662) then
             nsta = 65
           ELSEIF(br .le. 0.04912/0.08662) then
             nsta = 66
           ELSEIF(br .le. 0.05692/0.08662) then
             nsta = 70
           ELSEIF(br .le. 0.06262/0.08662) then
             nsta = 79
           ELSEIF(br .le. 0.06632/0.08662) then
             nsta = 80
           ELSEIF(br .le. 0.07052/0.08662) then
             nsta = 90
           ELSEIF(br .le. 0.07152/0.08662) then
             nsta = 92
           ELSEIF(br .le. 0.07312/0.08662) then
             nsta = 93
           ELSEIF(br .le. 0.07732/0.08662) then
             nsta = 94
           ELSEIF(br .le. 0.08082/0.08662) then
             nsta = 98
           ELSEIF(br .le. 0.08292/0.08662) then
             nsta = 106
           ELSEIF(br .le. 0.08502/0.08662) then
             nsta = 122
           ELSE
             nsta = 124
           ENDIF
         Case(44102)
               IF(br .le. 0.00033/0.22876) then
             nsta = 0
           ELSEIF(br .le. 0.00144/0.22876) then
             nsta = 1
           ELSEIF(br .le. 0.00834/0.22876) then
             nsta = 4
           ELSEIF(br .le. 0.01274/0.22876) then
             nsta = 5
           ELSEIF(br .le. 0.01397/0.22876) then
             nsta = 6
           ELSEIF(br .le. 0.0148/0.22876) then
             nsta = 8
           ELSEIF(br .le. 0.027/0.22876) then
             nsta = 13
           ELSEIF(br .le. 0.02953/0.22876) then
             nsta = 14
           ELSEIF(br .le. 0.07053/0.22876) then
             nsta = 17
           ELSEIF(br .le. 0.07893/0.22876) then
             nsta = 18
           ELSEIF(br .le. 0.08034/0.22876) then
             nsta = 19
           ELSEIF(br .le. 0.08059/0.22876) then
             nsta = 20
           ELSEIF(br .le. 0.08859/0.22876) then
             nsta = 26
           ELSEIF(br .le. 0.18159/0.22876) then
             nsta = 29
           ELSEIF(br .le. 0.18659/0.22876) then
             nsta = 30
           ELSEIF(br .le. 0.18962/0.22876) then
             nsta = 31
           ELSEIF(br .le. 0.19154/0.22876) then
             nsta = 33
           ELSEIF(br .le. 0.19535/0.22876) then
             nsta = 36
           ELSEIF(br .le. 0.20025/0.22876) then
             nsta = 38
           ELSEIF(br .le. 0.20156/0.22876) then
             nsta = 39
           ELSEIF(br .le. 0.21856/0.22876) then
             nsta = 48
           ELSEIF(br .le. 0.22556/0.22876) then
             nsta = 53
           ELSE
             nsta = 54
           ENDIF
         Case(44103)
               IF(br .le. 0.0008/0.15486) then
             nsta = 2
           ELSEIF(br .le. 0.0018/0.15486) then
             nsta = 3
           ELSEIF(br .le. 0.0037/0.15486) then
             nsta = 7
           ELSEIF(br .le. 0.00443/0.15486) then
             nsta = 9
           ELSEIF(br .le. 0.00486/0.15486) then
             nsta = 10
           ELSEIF(br .le. 0.00606/0.15486) then
             nsta = 14
           ELSEIF(br .le. 0.00726/0.15486) then
             nsta = 15
           ELSEIF(br .le. 0.00946/0.15486) then
             nsta = 17
           ELSEIF(br .le. 0.01046/0.15486) then
             nsta = 19
           ELSEIF(br .le. 0.01116/0.15486) then
             nsta = 22
           ELSEIF(br .le. 0.01236/0.15486) then
             nsta = 25
           ELSEIF(br .le. 0.01496/0.15486) then
             nsta = 28
           ELSEIF(br .le. 0.01586/0.15486) then
             nsta = 32
           ELSEIF(br .le. 0.01666/0.15486) then
             nsta = 33
           ELSEIF(br .le. 0.01766/0.15486) then
             nsta = 48
           ELSEIF(br .le. 0.01806/0.15486) then
             nsta = 53
           ELSEIF(br .le. 0.01906/0.15486) then
             nsta = 63
           ELSEIF(br .le. 0.03326/0.15486) then
             nsta = 65
           ELSEIF(br .le. 0.04326/0.15486) then
             nsta = 68
           ELSEIF(br .le. 0.04596/0.15486) then
             nsta = 69
           ELSEIF(br .le. 0.06466/0.15486) then
             nsta = 75
           ELSEIF(br .le. 0.06896/0.15486) then
             nsta = 79
           ELSEIF(br .le. 0.07376/0.15486) then
             nsta = 80
           ELSEIF(br .le. 0.08276/0.15486) then
             nsta = 85
           ELSEIF(br .le. 0.09776/0.15486) then
             nsta = 86
           ELSEIF(br .le. 0.10576/0.15486) then
             nsta = 88
           ELSEIF(br .le. 0.12026/0.15486) then
             nsta = 89
           ELSEIF(br .le. 0.13336/0.15486) then
             nsta = 99
           ELSEIF(br .le. 0.14536/0.15486) then
             nsta = 101
           ELSE
             nsta = 114
           ENDIF
         Case(44105)
               IF(br .le. 0.0077/0.1061) then
             nsta = 0
           ELSEIF(br .le. 0.019/0.1061) then
             nsta = 3
           ELSEIF(br .le. 0.0199/0.1061) then
             nsta = 9
           ELSEIF(br .le. 0.0242/0.1061) then
             nsta = 11
           ELSEIF(br .le. 0.026/0.1061) then
             nsta = 13
           ELSEIF(br .le. 0.0293/0.1061) then
             nsta = 14
           ELSEIF(br .le. 0.0332/0.1061) then
             nsta = 25
           ELSEIF(br .le. 0.0345/0.1061) then
             nsta = 26
           ELSEIF(br .le. 0.0349/0.1061) then
             nsta = 27
           ELSEIF(br .le. 0.0489/0.1061) then
             nsta = 30
           ELSEIF(br .le. 0.0689/0.1061) then
             nsta = 34
           ELSEIF(br .le. 0.081/0.1061) then
             nsta = 42
           ELSEIF(br .le. 0.0815/0.1061) then
             nsta = 43
           ELSEIF(br .le. 0.0912/0.1061) then
             nsta = 45
           ELSEIF(br .le. 0.0967/0.1061) then
             nsta = 46
           ELSEIF(br .le. 0.0971/0.1061) then
             nsta = 47
           ELSEIF(br .le. 0.1025/0.1061) then
             nsta = 49
           ELSEIF(br .le. 0.1038/0.1061) then
             nsta = 52
           ELSEIF(br .le. 0.1041/0.1061) then
             nsta = 55
           ELSEIF(br .le. 0.1046/0.1061) then
             nsta = 56
           ELSE
             nsta = 58
           ENDIF
         Case(45104)
               IF(br .le. 0.17/7.6857) then
             nsta = 0
           ELSEIF(br .le. 0.1974/7.6857) then
             nsta = 1
           ELSEIF(br .le. 0.2025/7.6857) then
             nsta = 2
           ELSEIF(br .le. 0.2195/7.6857) then
             nsta = 5
           ELSEIF(br .le. 0.2234/7.6857) then
             nsta = 6
           ELSEIF(br .le. 0.6934/7.6857) then
             nsta = 8
           ELSEIF(br .le. 0.7024/7.6857) then
             nsta = 10
           ELSEIF(br .le. 0.7584/7.6857) then
             nsta = 18
           ELSEIF(br .le. 0.7661/7.6857) then
             nsta = 30
           ELSEIF(br .le. 0.7833/7.6857) then
             nsta = 34
           ELSEIF(br .le. 0.7936/7.6857) then
             nsta = 36
           ELSEIF(br .le. 0.7975/7.6857) then
             nsta = 37
           ELSEIF(br .le. 0.8835/7.6857) then
             nsta = 41
           ELSEIF(br .le. 0.8874/7.6857) then
             nsta = 43
           ELSEIF(br .le. 1.3474/7.6857) then
             nsta = 49
           ELSEIF(br .le. 1.4194/7.6857) then
             nsta = 51
           ELSEIF(br .le. 1.6184/7.6857) then
             nsta = 53
           ELSEIF(br .le. 1.6248/7.6857) then
             nsta = 56
           ELSEIF(br .le. 1.6466/7.6857) then
             nsta = 56
           ELSEIF(br .le. 1.6633/7.6857) then
             nsta = 57
           ELSEIF(br .le. 1.6761/7.6857) then
             nsta = 59
           ELSEIF(br .le. 2.5661/7.6857) then
             nsta = 59
           ELSEIF(br .le. 2.6971/7.6857) then
             nsta = 60
           ELSEIF(br .le. 2.701/7.6857) then
             nsta = 61
           ELSEIF(br .le. 2.7049/7.6857) then
             nsta = 62
           ELSEIF(br .le. 2.8369/7.6857) then
             nsta = 63
           ELSEIF(br .le. 2.8889/7.6857) then
             nsta = 64
           ELSEIF(br .le. 2.9129/7.6857) then
             nsta = 67
           ELSEIF(br .le. 3.6629/7.6857) then
             nsta = 68
           ELSEIF(br .le. 3.6783/7.6857) then
             nsta = 69
           ELSEIF(br .le. 3.7053/7.6857) then
             nsta = 69
           ELSEIF(br .le. 3.7237/7.6857) then
             nsta = 72
           ELSEIF(br .le. 3.7477/7.6857) then
             nsta = 73
           ELSEIF(br .le. 4.0257/7.6857) then
             nsta = 75
           ELSEIF(br .le. 4.0477/7.6857) then
             nsta = 76
           ELSEIF(br .le. 4.0527/7.6857) then
             nsta = 79
           ELSEIF(br .le. 4.6327/7.6857) then
             nsta = 80
           ELSEIF(br .le. 4.6867/7.6857) then
             nsta = 81
           ELSEIF(br .le. 4.7077/7.6857) then
             nsta = 82
           ELSEIF(br .le. 5.5877/7.6857) then
             nsta = 84
           ELSEIF(br .le. 5.7217/7.6857) then
             nsta = 91
           ELSEIF(br .le. 5.7857/7.6857) then
             nsta = 93
           ELSEIF(br .le. 7.0957/7.6857) then
             nsta = 96
           ELSE
             nsta = 105
           ENDIF
         Case(46106)
               IF(br .le. 0.034/0.1238) then
             nsta = 13
           ELSEIF(br .le. 0.0423/0.1238) then
             nsta = 27
           ELSEIF(br .le. 0.0569/0.1238) then
             nsta = 28
           ELSEIF(br .le. 0.0664/0.1238) then
             nsta = 30
           ELSEIF(br .le. 0.0727/0.1238) then
             nsta = 32
           ELSEIF(br .le. 0.0897/0.1238) then
             nsta = 61
           ELSEIF(br .le. 0.1077/0.1238) then
             nsta = 64
           ELSE
             nsta = 76
           ENDIF
         Case(46109)
               IF(br .le. 0.0018/0.7297) then
             nsta = 0
           ELSEIF(br .le. 0.0128/0.7297) then
             nsta = 5
           ELSEIF(br .le. 0.0448/0.7297) then
             nsta = 9
           ELSEIF(br .le. 0.0768/0.7297) then
             nsta = 16
           ELSEIF(br .le. 0.0798/0.7297) then
             nsta = 17
           ELSEIF(br .le. 0.0813/0.7297) then
             nsta = 18
           ELSEIF(br .le. 0.0959/0.7297) then
             nsta = 24
           ELSEIF(br .le. 0.1129/0.7297) then
             nsta = 27
           ELSEIF(br .le. 0.1739/0.7297) then
             nsta = 35
           ELSEIF(br .le. 0.1789/0.7297) then
             nsta = 39
           ELSEIF(br .le. 0.1859/0.7297) then
             nsta = 41
           ELSEIF(br .le. 0.19/0.7297) then
             nsta = 44
           ELSEIF(br .le. 0.204/0.7297) then
             nsta = 45
           ELSEIF(br .le. 0.216/0.7297) then
             nsta = 46
           ELSEIF(br .le. 0.248/0.7297) then
             nsta = 48
           ELSEIF(br .le. 0.255/0.7297) then
             nsta = 50
           ELSEIF(br .le. 0.2562/0.7297) then
             nsta = 51
           ELSEIF(br .le. 0.3682/0.7297) then
             nsta = 58
           ELSEIF(br .le. 0.3802/0.7297) then
             nsta = 60
           ELSEIF(br .le. 0.3946/0.7297) then
             nsta = 63
           ELSEIF(br .le. 0.4053/0.7297) then
             nsta = 64
           ELSEIF(br .le. 0.4133/0.7297) then
             nsta = 66
           ELSEIF(br .le. 0.4233/0.7297) then
             nsta = 68
           ELSEIF(br .le. 0.4343/0.7297) then
             nsta = 69
           ELSEIF(br .le. 0.4543/0.7297) then
             nsta = 71
           ELSEIF(br .le. 0.4613/0.7297) then
             nsta = 75
           ELSEIF(br .le. 0.4713/0.7297) then
             nsta = 77
           ELSEIF(br .le. 0.4853/0.7297) then
             nsta = 78
           ELSEIF(br .le. 0.4923/0.7297) then
             nsta = 81
           ELSEIF(br .le. 0.5073/0.7297) then
             nsta = 85
           ELSEIF(br .le. 0.5148/0.7297) then
             nsta = 90
           ELSEIF(br .le. 0.5292/0.7297) then
             nsta = 91
           ELSEIF(br .le. 0.5353/0.7297) then
             nsta = 96
           ELSEIF(br .le. 0.5503/0.7297) then
             nsta = 98
           ELSEIF(br .le. 0.5683/0.7297) then
             nsta = 100
           ELSEIF(br .le. 0.5843/0.7297) then
             nsta = 104
           ELSEIF(br .le. 0.5963/0.7297) then
             nsta = 107
           ELSEIF(br .le. 0.6053/0.7297) then
             nsta = 109
           ELSEIF(br .le. 0.6183/0.7297) then
             nsta = 111
           ELSEIF(br .le. 0.6403/0.7297) then
             nsta = 111
           ELSEIF(br .le. 0.6486/0.7297) then
             nsta = 120
           ELSEIF(br .le. 0.6766/0.7297) then
             nsta = 125
           ELSEIF(br .le. 0.6816/0.7297) then
             nsta = 130
           ELSEIF(br .le. 0.6906/0.7297) then
             nsta = 130
           ELSEIF(br .le. 0.6996/0.7297) then
             nsta = 130
           ELSEIF(br .le. 0.7056/0.7297) then
             nsta = 131
           ELSEIF(br .le. 0.7186/0.7297) then
             nsta = 131
           ELSE
             nsta = 131
           ENDIF
         Case(46111)
               IF(br .le. 0.012/0.0205) then
             nsta = 4
           ELSE
             nsta = 6
           ENDIF
         Case(47108)
               IF(br .le. 0.284/1.596) then
             nsta = 0
           ELSEIF(br .le. 0.575/1.596) then
             nsta = 4
           ELSEIF(br .le. 0.678/1.596) then
             nsta = 5
           ELSEIF(br .le. 0.799/1.596) then
             nsta = 9
           ELSEIF(br .le. 0.808/1.596) then
             nsta = 12
           ELSEIF(br .le. 0.887/1.596) then
             nsta = 16
           ELSEIF(br .le. 0.908/1.596) then
             nsta = 21
           ELSEIF(br .le. 0.991/1.596) then
             nsta = 31
           ELSEIF(br .le. 1.001/1.596) then
             nsta = 32
           ELSEIF(br .le. 1.0059/1.596) then
             nsta = 40
           ELSEIF(br .le. 1.0539/1.596) then
             nsta = 46
           ELSEIF(br .le. 1.0689/1.596) then
             nsta = 48
           ELSEIF(br .le. 1.0745/1.596) then
             nsta = 49
           ELSEIF(br .le. 1.085/1.596) then
             nsta = 50
           ELSEIF(br .le. 1.098/1.596) then
             nsta = 56
           ELSEIF(br .le. 1.111/1.596) then
             nsta = 58
           ELSEIF(br .le. 1.124/1.596) then
             nsta = 62
           ELSEIF(br .le. 1.234/1.596) then
             nsta = 63
           ELSEIF(br .le. 1.263/1.596) then
             nsta = 66
           ELSEIF(br .le. 1.284/1.596) then
             nsta = 68
           ELSEIF(br .le. 1.43/1.596) then
             nsta = 71
           ELSEIF(br .le. 1.489/1.596) then
             nsta = 73
           ELSEIF(br .le. 1.552/1.596) then
             nsta = 77
           ELSE
             nsta = 83
           ENDIF
         Case(47110)
               IF(br .le. 0.083/4.918) then
             nsta = 1
           ELSEIF(br .le. 0.116/4.918) then
             nsta = 5
           ELSEIF(br .le. 0.161/4.918) then
             nsta = 7
           ELSEIF(br .le. 0.42/4.918) then
             nsta = 10
           ELSEIF(br .le. 0.479/4.918) then
             nsta = 16
           ELSEIF(br .le. 0.522/4.918) then
             nsta = 18
           ELSEIF(br .le. 0.551/4.918) then
             nsta = 27
           ELSEIF(br .le. 0.607/4.918) then
             nsta = 29
           ELSEIF(br .le. 0.615/4.918) then
             nsta = 38
           ELSEIF(br .le. 0.669/4.918) then
             nsta = 39
           ELSEIF(br .le. 0.687/4.918) then
             nsta = 40
           ELSEIF(br .le. 0.707/4.918) then
             nsta = 41
           ELSEIF(br .le. 0.787/4.918) then
             nsta = 47
           ELSEIF(br .le. 0.844/4.918) then
             nsta = 50
           ELSEIF(br .le. 1.507/4.918) then
             nsta = 53
           ELSEIF(br .le. 1.587/4.918) then
             nsta = 56
           ELSEIF(br .le. 1.837/4.918) then
             nsta = 57
           ELSEIF(br .le. 1.991/4.918) then
             nsta = 61
           ELSEIF(br .le. 2.075/4.918) then
             nsta = 67
           ELSEIF(br .le. 2.092/4.918) then
             nsta = 75
           ELSEIF(br .le. 2.114/4.918) then
             nsta = 76
           ELSEIF(br .le. 2.627/4.918) then
             nsta = 77
           ELSEIF(br .le. 2.852/4.918) then
             nsta = 78
           ELSEIF(br .le. 3.081/4.918) then
             nsta = 79
           ELSEIF(br .le. 3.797/4.918) then
             nsta = 80
           ELSEIF(br .le. 3.996/4.918) then
             nsta = 81
           ELSEIF(br .le. 4.045/4.918) then
             nsta = 81
           ELSEIF(br .le. 4.085/4.918) then
             nsta = 82
           ELSEIF(br .le. 4.293/4.918) then
             nsta = 82
           ELSEIF(br .le. 4.595/4.918) then
             nsta = 83
           ELSEIF(br .le. 4.706/4.918) then
             nsta = 84
           ELSEIF(br .le. 4.812/4.918) then
             nsta = 84
           ELSE
             nsta = 84
           ENDIF
         Case(48111)
               IF(br .le. 2.2/6.91) then
             nsta = 100
           ELSEIF(br .le. 2.81/6.91) then
             nsta = 103
           ELSEIF(br .le. 5.21/6.91) then
             nsta = 111
           ELSE
             nsta = 131
           ENDIF
         Case(48114)
               IF(br .le. 7.7/240.8) then
             nsta = 0
           ELSEIF(br .le. 19.6/240.8) then
             nsta = 1
           ELSEIF(br .le. 26.8/240.8) then
             nsta = 3
           ELSEIF(br .le. 35.6/240.8) then
             nsta = 5
           ELSEIF(br .le. 52.6/240.8) then
             nsta = 6
           ELSEIF(br .le. 53.4/240.8) then
             nsta = 9
           ELSEIF(br .le. 56.7/240.8) then
             nsta = 15
           ELSEIF(br .le. 72.7/240.8) then
             nsta = 18
           ELSEIF(br .le. 83.3/240.8) then
             nsta = 27
           ELSEIF(br .le. 90.3/240.8) then
             nsta = 28
           ELSEIF(br .le. 94.1/240.8) then
             nsta = 33
           ELSEIF(br .le. 99.5/240.8) then
             nsta = 36
           ELSEIF(br .le. 100.3/240.8) then
             nsta = 39
           ELSEIF(br .le. 102.9/240.8) then
             nsta = 40
           ELSEIF(br .le. 104.7/240.8) then
             nsta = 48
           ELSEIF(br .le. 109./240.8) then
             nsta = 51
           ELSEIF(br .le. 111.9/240.8) then
             nsta = 69
           ELSEIF(br .le. 115.5/240.8) then
             nsta = 74
           ELSEIF(br .le. 125.5/240.8) then
             nsta = 78
           ELSEIF(br .le. 144.8/240.8) then
             nsta = 78
           ELSEIF(br .le. 148.2/240.8) then
             nsta = 87
           ELSEIF(br .le. 153.2/240.8) then
             nsta = 93
           ELSEIF(br .le. 222.3/240.8) then
             nsta = 96
           ELSE
             nsta = 104
           ENDIF
         Case(49114)
               IF(br .le. 0.00029/0.0393) then
             nsta = 1
           ELSEIF(br .le. 0.00074/0.0393) then
             nsta = 2
           ELSEIF(br .le. 0.00164/0.0393) then
             nsta = 4
           ELSEIF(br .le. 0.00177/0.0393) then
             nsta = 7
           ELSEIF(br .le. 0.00397/0.0393) then
             nsta = 14
           ELSEIF(br .le. 0.00446/0.0393) then
             nsta = 16
           ELSEIF(br .le. 0.00716/0.0393) then
             nsta = 19
           ELSEIF(br .le. 0.00745/0.0393) then
             nsta = 22
           ELSEIF(br .le. 0.00778/0.0393) then
             nsta = 23
           ELSEIF(br .le. 0.01428/0.0393) then
             nsta = 26
           ELSEIF(br .le. 0.01513/0.0393) then
             nsta = 28
           ELSEIF(br .le. 0.01583/0.0393) then
             nsta = 29
           ELSEIF(br .le. 0.01693/0.0393) then
             nsta = 32
           ELSEIF(br .le. 0.01718/0.0393) then
             nsta = 37
           ELSEIF(br .le. 0.01808/0.0393) then
             nsta = 38
           ELSEIF(br .le. 0.0188/0.0393) then
             nsta = 39
           ELSEIF(br .le. 0.01934/0.0393) then
             nsta = 40
           ELSEIF(br .le. 0.02134/0.0393) then
             nsta = 43
           ELSEIF(br .le. 0.02274/0.0393) then
             nsta = 44
           ELSEIF(br .le. 0.02564/0.0393) then
             nsta = 48
           ELSEIF(br .le. 0.02584/0.0393) then
             nsta = 51
           ELSEIF(br .le. 0.02604/0.0393) then
             nsta = 52
           ELSEIF(br .le. 0.02644/0.0393) then
             nsta = 54
           ELSEIF(br .le. 0.02734/0.0393) then
             nsta = 56
           ELSEIF(br .le. 0.02796/0.0393) then
             nsta = 58
           ELSEIF(br .le. 0.02854/0.0393) then
             nsta = 60
           ELSEIF(br .le. 0.02899/0.0393) then
             nsta = 61
           ELSEIF(br .le. 0.02932/0.0393) then
             nsta = 63
           ELSEIF(br .le. 0.02968/0.0393) then
             nsta = 65
           ELSEIF(br .le. 0.03148/0.0393) then
             nsta = 67
           ELSEIF(br .le. 0.0317/0.0393) then
             nsta = 69
           ELSEIF(br .le. 0.0323/0.0393) then
             nsta = 70
           ELSEIF(br .le. 0.03244/0.0393) then
             nsta = 70
           ELSEIF(br .le. 0.03344/0.0393) then
             nsta = 71
           ELSEIF(br .le. 0.03435/0.0393) then
             nsta = 71
           ELSEIF(br .le. 0.035/0.0393) then
             nsta = 72
           ELSEIF(br .le. 0.03516/0.0393) then
             nsta = 73
           ELSEIF(br .le. 0.03552/0.0393) then
             nsta = 74
           ELSEIF(br .le. 0.03574/0.0393) then
             nsta = 75
           ELSEIF(br .le. 0.03574/0.0393) then
             nsta = 76
           ELSEIF(br .le. 0.03612/0.0393) then
             nsta = 78
           ELSEIF(br .le. 0.03634/0.0393) then
             nsta = 81
           ELSEIF(br .le. 0.03703/0.0393) then
             nsta = 83
           ELSEIF(br .le. 0.03719/0.0393) then
             nsta = 84
           ELSEIF(br .le. 0.03739/0.0393) then
             nsta = 85
           ELSEIF(br .le. 0.03819/0.0393) then
             nsta = 86
           ELSEIF(br .le. 0.03879/0.0393) then
             nsta = 87
           ELSE
             nsta = 88
           ENDIF
         Case(49116)
               IF(br .le. 0.025/6.957) then
             nsta = 2
           ELSEIF(br .le. 0.036/6.957) then
             nsta = 5
           ELSEIF(br .le. 0.186/6.957) then
             nsta = 8
           ELSEIF(br .le. 0.306/6.957) then
             nsta = 11
           ELSEIF(br .le. 0.313/6.957) then
             nsta = 13
           ELSEIF(br .le. 0.398/6.957) then
             nsta = 14
           ELSEIF(br .le. 0.487/6.957) then
             nsta = 16
           ELSEIF(br .le. 0.514/6.957) then
             nsta = 17
           ELSEIF(br .le. 0.586/6.957) then
             nsta = 21
           ELSEIF(br .le. 0.611/6.957) then
             nsta = 27
           ELSEIF(br .le. 2.711/6.957) then
             nsta = 32
           ELSEIF(br .le. 2.782/6.957) then
             nsta = 36
           ELSEIF(br .le. 2.945/6.957) then
             nsta = 39
           ELSEIF(br .le. 3.045/6.957) then
             nsta = 40
           ELSEIF(br .le. 3.347/6.957) then
             nsta = 41
           ELSEIF(br .le. 3.439/6.957) then
             nsta = 43
           ELSEIF(br .le. 3.531/6.957) then
             nsta = 44
           ELSEIF(br .le. 3.58/6.957) then
             nsta = 46
           ELSEIF(br .le. 3.605/6.957) then
             nsta = 48
           ELSEIF(br .le. 3.621/6.957) then
             nsta = 52
           ELSEIF(br .le. 3.645/6.957) then
             nsta = 54
           ELSEIF(br .le. 3.669/6.957) then
             nsta = 56
           ELSEIF(br .le. 3.685/6.957) then
             nsta = 57
           ELSEIF(br .le. 3.905/6.957) then
             nsta = 58
           ELSEIF(br .le. 4.045/6.957) then
             nsta = 59
           ELSEIF(br .le. 4.063/6.957) then
             nsta = 60
           ELSEIF(br .le. 4.318/6.957) then
             nsta = 62
           ELSEIF(br .le. 4.553/6.957) then
             nsta = 64
           ELSEIF(br .le. 4.627/6.957) then
             nsta = 66
           ELSEIF(br .le. 5.157/6.957) then
             nsta = 67
           ELSEIF(br .le. 5.195/6.957) then
             nsta = 68
           ELSEIF(br .le. 5.705/6.957) then
             nsta = 69
           ELSEIF(br .le. 6.067/6.957) then
             nsta = 70
           ELSE
             nsta = 71
           ENDIF
         Case(50116)
               IF(br .le. 0.00025/0.04263) then
             nsta = 0
           ELSEIF(br .le. 0.00045/0.04263) then
             nsta = 1
           ELSEIF(br .le. 0.00055/0.04263) then
             nsta = 2
           ELSEIF(br .le. 0.00192/0.04263) then
             nsta = 4
           ELSEIF(br .le. 0.00258/0.04263) then
             nsta = 5
           ELSEIF(br .le. 0.00361/0.04263) then
             nsta = 10
           ELSEIF(br .le. 0.00379/0.04263) then
             nsta = 11
           ELSEIF(br .le. 0.00424/0.04263) then
             nsta = 12
           ELSEIF(br .le. 0.00434/0.04263) then
             nsta = 14
           ELSEIF(br .le. 0.00448/0.04263) then
             nsta = 16
           ELSEIF(br .le. 0.00616/0.04263) then
             nsta = 18
           ELSEIF(br .le. 0.00639/0.04263) then
             nsta = 23
           ELSEIF(br .le. 0.00779/0.04263) then
             nsta = 31
           ELSEIF(br .le. 0.01009/0.04263) then
             nsta = 33
           ELSEIF(br .le. 0.01023/0.04263) then
             nsta = 39
           ELSEIF(br .le. 0.01182/0.04263) then
             nsta = 40
           ELSEIF(br .le. 0.01218/0.04263) then
             nsta = 45
           ELSEIF(br .le. 0.01288/0.04263) then
             nsta = 48
           ELSEIF(br .le. 0.0135/0.04263) then
             nsta = 51
           ELSEIF(br .le. 0.01373/0.04263) then
             nsta = 59
           ELSEIF(br .le. 0.01397/0.04263) then
             nsta = 60
           ELSEIF(br .le. 0.0162/0.04263) then
             nsta = 65
           ELSEIF(br .le. 0.0165/0.04263) then
             nsta = 67
           ELSEIF(br .le. 0.01668/0.04263) then
             nsta = 73
           ELSEIF(br .le. 0.01682/0.04263) then
             nsta = 77
           ELSEIF(br .le. 0.01705/0.04263) then
             nsta = 84
           ELSEIF(br .le. 0.01746/0.04263) then
             nsta = 85
           ELSEIF(br .le. 0.01796/0.04263) then
             nsta = 86
           ELSEIF(br .le. 0.01812/0.04263) then
             nsta = 89
           ELSEIF(br .le. 0.02022/0.04263) then
             nsta = 92
           ELSEIF(br .le. 0.02124/0.04263) then
             nsta = 93
           ELSEIF(br .le. 0.02151/0.04263) then
             nsta = 96
           ELSEIF(br .le. 0.02165/0.04263) then
             nsta = 98
           ELSEIF(br .le. 0.02174/0.04263) then
             nsta = 99
           ELSEIF(br .le. 0.02365/0.04263) then
             nsta = 101
           ELSEIF(br .le. 0.02397/0.04263) then
             nsta = 104
           ELSEIF(br .le. 0.02827/0.04263) then
             nsta = 108
           ELSEIF(br .le. 0.02856/0.04263) then
             nsta = 109
           ELSEIF(br .le. 0.02929/0.04263) then
             nsta = 110
           ELSEIF(br .le. 0.0297/0.04263) then
             nsta = 112
           ELSEIF(br .le. 0.03002/0.04263) then
             nsta = 113
           ELSEIF(br .le. 0.03013/0.04263) then
             nsta = 117
           ELSEIF(br .le. 0.03073/0.04263) then
             nsta = 120
           ELSEIF(br .le. 0.03115/0.04263) then
             nsta = 121
           ELSEIF(br .le. 0.03223/0.04263) then
             nsta = 123
           ELSEIF(br .le. 0.03238/0.04263) then
             nsta = 126
           ELSEIF(br .le. 0.03296/0.04263) then
             nsta = 127
           ELSEIF(br .le. 0.03347/0.04263) then
             nsta = 128
           ELSEIF(br .le. 0.03465/0.04263) then
             nsta = 129
           ELSEIF(br .le. 0.03495/0.04263) then
             nsta = 134
           ELSEIF(br .le. 0.03518/0.04263) then
             nsta = 137
           ELSEIF(br .le. 0.03532/0.04263) then
             nsta = 140
           ELSEIF(br .le. 0.03554/0.04263) then
             nsta = 143
           ELSEIF(br .le. 0.03568/0.04263) then
             nsta = 147
           ELSEIF(br .le. 0.03595/0.04263) then
             nsta = 156
           ELSEIF(br .le. 0.03622/0.04263) then
             nsta = 162
           ELSEIF(br .le. 0.03651/0.04263) then
             nsta = 165
           ELSEIF(br .le. 0.03665/0.04263) then
             nsta = 167
           ELSEIF(br .le. 0.03742/0.04263) then
             nsta = 180
           ELSEIF(br .le. 0.03774/0.04263) then
             nsta = 181
           ELSEIF(br .le. 0.03795/0.04263) then
             nsta = 189
           ELSEIF(br .le. 0.03902/0.04263) then
             nsta = 190
           ELSEIF(br .le. 0.03965/0.04263) then
             nsta = 207
           ELSEIF(br .le. 0.04013/0.04263) then
             nsta = 210
           ELSEIF(br .le. 0.04113/0.04263) then
             nsta = 214
           ELSEIF(br .le. 0.04213/0.04263) then
             nsta = 236
           ELSEIF(br .le. 0.04247/0.04263) then
             nsta = 257
           ELSE
             nsta = 258
           ENDIF
         Case(51122)
               IF(br .le. 0.0102/0.42427) then
             nsta = 0
           ELSEIF(br .le. 0.0192/0.42427) then
             nsta = 1
           ELSEIF(br .le. 0.0632/0.42427) then
             nsta = 2
           ELSEIF(br .le. 0.06357/0.42427) then
             nsta = 3
           ELSEIF(br .le. 0.06437/0.42427) then
             nsta = 6
           ELSEIF(br .le. 0.06767/0.42427) then
             nsta = 7
           ELSEIF(br .le. 0.06917/0.42427) then
             nsta = 8
           ELSEIF(br .le. 0.07137/0.42427) then
             nsta = 9
           ELSEIF(br .le. 0.07387/0.42427) then
             nsta = 10
           ELSEIF(br .le. 0.14887/0.42427) then
             nsta = 14
           ELSEIF(br .le. 0.15647/0.42427) then
             nsta = 15
           ELSEIF(br .le. 0.15907/0.42427) then
             nsta = 16
           ELSEIF(br .le. 0.16287/0.42427) then
             nsta = 17
           ELSEIF(br .le. 0.16687/0.42427) then
             nsta = 19
           ELSEIF(br .le. 0.16857/0.42427) then
             nsta = 28
           ELSEIF(br .le. 0.17517/0.42427) then
             nsta = 32
           ELSEIF(br .le. 0.18727/0.42427) then
             nsta = 33
           ELSEIF(br .le. 0.18897/0.42427) then
             nsta = 34
           ELSEIF(br .le. 0.20897/0.42427) then
             nsta = 40
           ELSEIF(br .le. 0.21297/0.42427) then
             nsta = 42
           ELSEIF(br .le. 0.21527/0.42427) then
             nsta = 43
           ELSEIF(br .le. 0.21747/0.42427) then
             nsta = 44
           ELSEIF(br .le. 0.27147/0.42427) then
             nsta = 46
           ELSEIF(br .le. 0.27587/0.42427) then
             nsta = 47
           ELSEIF(br .le. 0.28017/0.42427) then
             nsta = 48
           ELSEIF(br .le. 0.28127/0.42427) then
             nsta = 49
           ELSEIF(br .le. 0.28217/0.42427) then
             nsta = 50
           ELSEIF(br .le. 0.29287/0.42427) then
             nsta = 52
           ELSEIF(br .le. 0.30217/0.42427) then
             nsta = 53
           ELSEIF(br .le. 0.31317/0.42427) then
             nsta = 54
           ELSEIF(br .le. 0.31467/0.42427) then
             nsta = 55
           ELSEIF(br .le. 0.31677/0.42427) then
             nsta = 56
           ELSEIF(br .le. 0.32067/0.42427) then
             nsta = 57
           ELSEIF(br .le. 0.33477/0.42427) then
             nsta = 59
           ELSEIF(br .le. 0.34217/0.42427) then
             nsta = 60
           ELSEIF(br .le. 0.34837/0.42427) then
             nsta = 61
           ELSEIF(br .le. 0.35317/0.42427) then
             nsta = 62
           ELSEIF(br .le. 0.36817/0.42427) then
             nsta = 63
           ELSEIF(br .le. 0.38417/0.42427) then
             nsta = 66
           ELSEIF(br .le. 0.38837/0.42427) then
             nsta = 67
           ELSEIF(br .le. 0.40937/0.42427) then
             nsta = 68
           ELSE
             nsta = 69
           ENDIF
         Case(51124)
               IF(br .le. 0.021/0.28921) then
             nsta = 0
           ELSEIF(br .le. 0.0287/0.28921) then
             nsta = 1
           ELSEIF(br .le. 0.02908/0.28921) then
             nsta = 3
           ELSEIF(br .le. 0.07308/0.28921) then
             nsta = 5
           ELSEIF(br .le. 0.09808/0.28921) then
             nsta = 6
           ELSEIF(br .le. 0.11508/0.28921) then
             nsta = 8
           ELSEIF(br .le. 0.11938/0.28921) then
             nsta = 12
           ELSEIF(br .le. 0.11974/0.28921) then
             nsta = 15
           ELSEIF(br .le. 0.12064/0.28921) then
             nsta = 16
           ELSEIF(br .le. 0.13864/0.28921) then
             nsta = 19
           ELSEIF(br .le. 0.14034/0.28921) then
             nsta = 21
           ELSEIF(br .le. 0.15834/0.28921) then
             nsta = 22
           ELSEIF(br .le. 0.16534/0.28921) then
             nsta = 23
           ELSEIF(br .le. 0.16714/0.28921) then
             nsta = 25
           ELSEIF(br .le. 0.16814/0.28921) then
             nsta = 26
           ELSEIF(br .le. 0.20214/0.28921) then
             nsta = 31
           ELSEIF(br .le. 0.20274/0.28921) then
             nsta = 37
           ELSEIF(br .le. 0.20317/0.28921) then
             nsta = 39
           ELSEIF(br .le. 0.20537/0.28921) then
             nsta = 41
           ELSEIF(br .le. 0.21737/0.28921) then
             nsta = 42
           ELSEIF(br .le. 0.22257/0.28921) then
             nsta = 43
           ELSEIF(br .le. 0.22293/0.28921) then
             nsta = 44
           ELSEIF(br .le. 0.22383/0.28921) then
             nsta = 45
           ELSEIF(br .le. 0.22415/0.28921) then
             nsta = 47
           ELSEIF(br .le. 0.22515/0.28921) then
             nsta = 49
           ELSEIF(br .le. 0.23315/0.28921) then
             nsta = 50
           ELSEIF(br .le. 0.23525/0.28921) then
             nsta = 51
           ELSEIF(br .le. 0.24245/0.28921) then
             nsta = 52
           ELSEIF(br .le. 0.25645/0.28921) then
             nsta = 53
           ELSEIF(br .le. 0.25855/0.28921) then
             nsta = 54
           ELSEIF(br .le. 0.25941/0.28921) then
             nsta = 56
           ELSEIF(br .le. 0.26721/0.28921) then
             nsta = 57
           ELSEIF(br .le. 0.26941/0.28921) then
             nsta = 58
           ELSEIF(br .le. 0.27211/0.28921) then
             nsta = 59
           ELSEIF(br .le. 0.27311/0.28921) then
             nsta = 60
           ELSEIF(br .le. 0.27741/0.28921) then
             nsta = 61
           ELSEIF(br .le. 0.28571/0.28921) then
             nsta = 62
           ELSEIF(br .le. 0.28591/0.28921) then
             nsta = 63
           ELSEIF(br .le. 0.28781/0.28921) then
             nsta = 64
           ELSE
             nsta = 65
           ENDIF
         Case(52123)
               IF(br .le. 0.000038/0.022502) then
             nsta = 0
           ELSEIF(br .le. 0.000113/0.022502) then
             nsta = 1
           ELSEIF(br .le. 0.000141/0.022502) then
             nsta = 4
           ELSEIF(br .le. 0.000221/0.022502) then
             nsta = 8
           ELSEIF(br .le. 0.000282/0.022502) then
             nsta = 20
           ELSEIF(br .le. 0.000343/0.022502) then
             nsta = 31
           ELSEIF(br .le. 0.000943/0.022502) then
             nsta = 34
           ELSEIF(br .le. 0.003743/0.022502) then
             nsta = 58
           ELSEIF(br .le. 0.009243/0.022502) then
             nsta = 61
           ELSEIF(br .le. 0.012243/0.022502) then
             nsta = 62
           ELSEIF(br .le. 0.012455/0.022502) then
             nsta = 65
           ELSEIF(br .le. 0.012818/0.022502) then
             nsta = 70
           ELSEIF(br .le. 0.013252/0.022502) then
             nsta = 72
           ELSEIF(br .le. 0.013652/0.022502) then
             nsta = 78
           ELSEIF(br .le. 0.013855/0.022502) then
             nsta = 80
           ELSEIF(br .le. 0.017755/0.022502) then
             nsta = 83
           ELSEIF(br .le. 0.018208/0.022502) then
             nsta = 85
           ELSEIF(br .le. 0.021308/0.022502) then
             nsta = 113
           ELSEIF(br .le. 0.021393/0.022502) then
             nsta = 116
           ELSEIF(br .le. 0.021663/0.022502) then
             nsta = 121
           ELSEIF(br .le. 0.022093/0.022502) then
             nsta = 124
           ELSEIF(br .le. 0.022233/0.022502) then
             nsta = 133
           ELSEIF(br .le. 0.022412/0.022502) then
             nsta = 146
           ELSE
             nsta = 147
           ENDIF
         Case(52124)
               IF(br .le. 0.001/0.43519) then
             nsta = 0
           ELSEIF(br .le. 0.0125/0.43519) then
             nsta = 1
           ELSEIF(br .le. 0.0295/0.43519) then
             nsta = 3
           ELSEIF(br .le. 0.02955/0.43519) then
             nsta = 4
           ELSEIF(br .le. 0.02987/0.43519) then
             nsta = 5
           ELSEIF(br .le. 0.03297/0.43519) then
             nsta = 9
           ELSEIF(br .le. 0.05997/0.43519) then
             nsta = 10
           ELSEIF(br .le. 0.06077/0.43519) then
             nsta = 11
           ELSEIF(br .le. 0.06119/0.43519) then
             nsta = 12
           ELSEIF(br .le. 0.06989/0.43519) then
             nsta = 19
           ELSEIF(br .le. 0.07059/0.43519) then
             nsta = 23
           ELSEIF(br .le. 0.07109/0.43519) then
             nsta = 28
           ELSEIF(br .le. 0.07199/0.43519) then
             nsta = 29
           ELSEIF(br .le. 0.07339/0.43519) then
             nsta = 35
           ELSEIF(br .le. 0.07589/0.43519) then
             nsta = 38
           ELSEIF(br .le. 0.07859/0.43519) then
             nsta = 44
           ELSEIF(br .le. 0.07979/0.43519) then
             nsta = 53
           ELSEIF(br .le. 0.08739/0.43519) then
             nsta = 59
           ELSEIF(br .le. 0.09069/0.43519) then
             nsta = 61
           ELSEIF(br .le. 0.09159/0.43519) then
             nsta = 67
           ELSEIF(br .le. 0.09239/0.43519) then
             nsta = 88
           ELSEIF(br .le. 0.09449/0.43519) then
             nsta = 91
           ELSEIF(br .le. 0.10019/0.43519) then
             nsta = 100
           ELSEIF(br .le. 0.10419/0.43519) then
             nsta = 108
           ELSEIF(br .le. 0.20319/0.43519) then
             nsta = 109
           ELSEIF(br .le. 0.22939/0.43519) then
             nsta = 125
           ELSEIF(br .le. 0.23409/0.43519) then
             nsta = 127
           ELSEIF(br .le. 0.24309/0.43519) then
             nsta = 130
           ELSEIF(br .le. 0.25509/0.43519) then
             nsta = 142
           ELSEIF(br .le. 0.26509/0.43519) then
             nsta = 153
           ELSEIF(br .le. 0.27429/0.43519) then
             nsta = 154
           ELSEIF(br .le. 0.28339/0.43519) then
             nsta = 161
           ELSEIF(br .le. 0.28689/0.43519) then
             nsta = 171
           ELSEIF(br .le. 0.30069/0.43519) then
             nsta = 174
           ELSEIF(br .le. 0.33469/0.43519) then
             nsta = 176
           ELSEIF(br .le. 0.33819/0.43519) then
             nsta = 183
           ELSEIF(br .le. 0.34189/0.43519) then
             nsta = 184
           ELSEIF(br .le. 0.35869/0.43519) then
             nsta = 192
           ELSEIF(br .le. 0.39569/0.43519) then
             nsta = 195
           ELSEIF(br .le. 0.40079/0.43519) then
             nsta = 196
           ELSEIF(br .le. 0.40459/0.43519) then
             nsta = 209
           ELSEIF(br .le. 0.42359/0.43519) then
             nsta = 212
           ELSEIF(br .le. 0.42559/0.43519) then
             nsta = 213
           ELSEIF(br .le. 0.43299/0.43519) then
             nsta = 224
           ELSE
             nsta = 253
           ENDIF
         Case(52125)
               IF(br .le. 0.0053/0.1979) then
             nsta = 0
           ELSEIF(br .le. 0.0161/0.1979) then
             nsta = 1
           ELSEIF(br .le. 0.0193/0.1979) then
             nsta = 5
           ELSEIF(br .le. 0.0201/0.1979) then
             nsta = 8
           ELSEIF(br .le. 0.026/0.1979) then
             nsta = 13
           ELSEIF(br .le. 0.0292/0.1979) then
             nsta = 28
           ELSEIF(br .le. 0.0403/0.1979) then
             nsta = 32
           ELSEIF(br .le. 0.0518/0.1979) then
             nsta = 41
           ELSEIF(br .le. 0.0548/0.1979) then
             nsta = 42
           ELSEIF(br .le. 0.0584/0.1979) then
             nsta = 44
           ELSEIF(br .le. 0.0703/0.1979) then
             nsta = 45
           ELSEIF(br .le. 0.0716/0.1979) then
             nsta = 48
           ELSEIF(br .le. 0.0754/0.1979) then
             nsta = 56
           ELSEIF(br .le. 0.0792/0.1979) then
             nsta = 58
           ELSEIF(br .le. 0.0914/0.1979) then
             nsta = 63
           ELSEIF(br .le. 0.0937/0.1979) then
             nsta = 65
           ELSEIF(br .le. 0.0998/0.1979) then
             nsta = 69
           ELSEIF(br .le. 0.1058/0.1979) then
             nsta = 75
           ELSEIF(br .le. 0.1085/0.1979) then
             nsta = 78
           ELSEIF(br .le. 0.1111/0.1979) then
             nsta = 81
           ELSEIF(br .le. 0.1187/0.1979) then
             nsta = 88
           ELSEIF(br .le. 0.1278/0.1979) then
             nsta = 102
           ELSEIF(br .le. 0.1325/0.1979) then
             nsta = 104
           ELSEIF(br .le. 0.1377/0.1979) then
             nsta = 112
           ELSEIF(br .le. 0.1466/0.1979) then
             nsta = 117
           ELSEIF(br .le. 0.1503/0.1979) then
             nsta = 120
           ELSEIF(br .le. 0.1522/0.1979) then
             nsta = 127
           ELSEIF(br .le. 0.1543/0.1979) then
             nsta = 129
           ELSEIF(br .le. 0.1585/0.1979) then
             nsta = 131
           ELSEIF(br .le. 0.1646/0.1979) then
             nsta = 135
           ELSEIF(br .le. 0.1697/0.1979) then
             nsta = 146
           ELSEIF(br .le. 0.1748/0.1979) then
             nsta = 147
           ELSEIF(br .le. 0.1789/0.1979) then
             nsta = 167
           ELSEIF(br .le. 0.1847/0.1979) then
             nsta = 175
           ELSEIF(br .le. 0.1908/0.1979) then
             nsta = 179
           ELSEIF(br .le. 0.1958/0.1979) then
             nsta = 185
           ELSE
             nsta = 208
           ENDIF
         Case(52129)
               IF(br .le. 0.0067/0.06195) then
             nsta = 47
           ELSEIF(br .le. 0.0092/0.06195) then
             nsta = 62
           ELSEIF(br .le. 0.0301/0.06195) then
             nsta = 67
           ELSEIF(br .le. 0.0409/0.06195) then
             nsta = 69
           ELSEIF(br .le. 0.0514/0.06195) then
             nsta = 76
           ELSEIF(br .le. 0.05275/0.06195) then
             nsta = 99
           ELSEIF(br .le. 0.05835/0.06195) then
             nsta = 102
           ELSE
             nsta = 112
           ENDIF
         Case(53128)
               IF(br .le. 0.00081/0.9766) then
             nsta = 0
           ELSEIF(br .le. 0.00227/0.9766) then
             nsta = 1
           ELSEIF(br .le. 0.00507/0.9766) then
             nsta = 2
           ELSEIF(br .le. 0.01007/0.9766) then
             nsta = 3
           ELSEIF(br .le. 0.04707/0.9766) then
             nsta = 4
           ELSEIF(br .le. 0.05117/0.9766) then
             nsta = 5
           ELSEIF(br .le. 0.06017/0.9766) then
             nsta = 6
           ELSEIF(br .le. 0.06113/0.9766) then
             nsta = 7
           ELSEIF(br .le. 0.06316/0.9766) then
             nsta = 8
           ELSEIF(br .le. 0.07016/0.9766) then
             nsta = 10
           ELSEIF(br .le. 0.07051/0.9766) then
             nsta = 11
           ELSEIF(br .le. 0.07078/0.9766) then
             nsta = 13
           ELSEIF(br .le. 0.07151/0.9766) then
             nsta = 18
           ELSEIF(br .le. 0.07351/0.9766) then
             nsta = 19
           ELSEIF(br .le. 0.0737/0.9766) then
             nsta = 19
           ELSEIF(br .le. 0.07597/0.9766) then
             nsta = 20
           ELSEIF(br .le. 0.08117/0.9766) then
             nsta = 21
           ELSEIF(br .le. 0.08367/0.9766) then
             nsta = 22
           ELSEIF(br .le. 0.08547/0.9766) then
             nsta = 24
           ELSEIF(br .le. 0.08608/0.9766) then
             nsta = 25
           ELSEIF(br .le. 0.0865/0.9766) then
             nsta = 26
           ELSEIF(br .le. 0.0965/0.9766) then
             nsta = 28
           ELSEIF(br .le. 0.09715/0.9766) then
             nsta = 31
           ELSEIF(br .le. 0.12115/0.9766) then
             nsta = 32
           ELSEIF(br .le. 0.12284/0.9766) then
             nsta = 35
           ELSEIF(br .le. 0.12357/0.9766) then
             nsta = 36
           ELSEIF(br .le. 0.13167/0.9766) then
             nsta = 38
           ELSEIF(br .le. 0.13271/0.9766) then
             nsta = 41
           ELSEIF(br .le. 0.13481/0.9766) then
             nsta = 42
           ELSEIF(br .le. 0.13504/0.9766) then
             nsta = 44
           ELSEIF(br .le. 0.13734/0.9766) then
             nsta = 45
           ELSEIF(br .le. 0.13876/0.9766) then
             nsta = 46
           ELSEIF(br .le. 0.14039/0.9766) then
             nsta = 48
           ELSEIF(br .le. 0.14097/0.9766) then
             nsta = 49
           ELSEIF(br .le. 0.14181/0.9766) then
             nsta = 54
           ELSEIF(br .le. 0.14212/0.9766) then
             nsta = 56
           ELSEIF(br .le. 0.14262/0.9766) then
             nsta = 58
           ELSEIF(br .le. 0.14277/0.9766) then
             nsta = 59
           ELSEIF(br .le. 0.14327/0.9766) then
             nsta = 60
           ELSEIF(br .le. 0.14557/0.9766) then
             nsta = 61
           ELSEIF(br .le. 0.14967/0.9766) then
             nsta = 63
           ELSEIF(br .le. 0.15457/0.9766) then
             nsta = 65
           ELSEIF(br .le. 0.15684/0.9766) then
             nsta = 67
           ELSEIF(br .le. 0.15964/0.9766) then
             nsta = 68
           ELSEIF(br .le. 0.16083/0.9766) then
             nsta = 70
           ELSEIF(br .le. 0.16613/0.9766) then
             nsta = 71
           ELSEIF(br .le. 0.16903/0.9766) then
             nsta = 73
           ELSEIF(br .le. 0.16934/0.9766) then
             nsta = 73
           ELSEIF(br .le. 0.17155/0.9766) then
             nsta = 75
           ELSEIF(br .le. 0.18455/0.9766) then
             nsta = 76
           ELSEIF(br .le. 0.18524/0.9766) then
             nsta = 77
           ELSEIF(br .le. 0.19094/0.9766) then
             nsta = 78
           ELSEIF(br .le. 0.19148/0.9766) then
             nsta = 80
           ELSEIF(br .le. 0.19209/0.9766) then
             nsta = 80
           ELSEIF(br .le. 0.19419/0.9766) then
             nsta = 80
           ELSEIF(br .le. 0.19503/0.9766) then
             nsta = 82
           ELSEIF(br .le. 0.19545/0.9766) then
             nsta = 84
           ELSEIF(br .le. 0.20195/0.9766) then
             nsta = 86
           ELSEIF(br .le. 0.20455/0.9766) then
             nsta = 88
           ELSEIF(br .le. 0.2052/0.9766) then
             nsta = 89
           ELSEIF(br .le. 0.20566/0.9766) then
             nsta = 90
           ELSEIF(br .le. 0.21466/0.9766) then
             nsta = 91
           ELSEIF(br .le. 0.21636/0.9766) then
             nsta = 92
           ELSEIF(br .le. 0.23636/0.9766) then
             nsta = 93
           ELSEIF(br .le. 0.23682/0.9766) then
             nsta = 96
           ELSEIF(br .le. 0.23782/0.9766) then
             nsta = 97
           ELSEIF(br .le. 0.23886/0.9766) then
             nsta = 98
           ELSEIF(br .le. 0.2397/0.9766) then
             nsta = 99
           ELSEIF(br .le. 0.2429/0.9766) then
             nsta = 100
           ELSEIF(br .le. 0.2456/0.9766) then
             nsta = 102
           ELSEIF(br .le. 0.2496/0.9766) then
             nsta = 103
           ELSEIF(br .le. 0.2501/0.9766) then
             nsta = 105
           ELSEIF(br .le. 0.2526/0.9766) then
             nsta = 106
           ELSEIF(br .le. 0.25387/0.9766) then
             nsta = 106
           ELSEIF(br .le. 0.26587/0.9766) then
             nsta = 109
           ELSEIF(br .le. 0.28687/0.9766) then
             nsta = 110
           ELSEIF(br .le. 0.28725/0.9766) then
             nsta = 111
           ELSEIF(br .le. 0.33125/0.9766) then
             nsta = 112
           ELSEIF(br .le. 0.33365/0.9766) then
             nsta = 113
           ELSEIF(br .le. 0.34865/0.9766) then
             nsta = 115
           ELSEIF(br .le. 0.35165/0.9766) then
             nsta = 116
           ELSEIF(br .le. 0.35605/0.9766) then
             nsta = 118
           ELSEIF(br .le. 0.36305/0.9766) then
             nsta = 119
           ELSEIF(br .le. 0.38105/0.9766) then
             nsta = 120
           ELSEIF(br .le. 0.38163/0.9766) then
             nsta = 121
           ELSEIF(br .le. 0.39963/0.9766) then
             nsta = 123
           ELSEIF(br .le. 0.40093/0.9766) then
             nsta = 124
           ELSEIF(br .le. 0.4022/0.9766) then
             nsta = 124
           ELSEIF(br .le. 0.40426/0.9766) then
             nsta = 125
           ELSEIF(br .le. 0.40676/0.9766) then
             nsta = 125
           ELSEIF(br .le. 0.40852/0.9766) then
             nsta = 126
           ELSEIF(br .le. 0.40994/0.9766) then
             nsta = 127
           ELSEIF(br .le. 0.41182/0.9766) then
             nsta = 127
           ELSEIF(br .le. 0.41402/0.9766) then
             nsta = 128
           ELSEIF(br .le. 0.42102/0.9766) then
             nsta = 128
           ELSEIF(br .le. 0.42217/0.9766) then
             nsta = 129
           ELSEIF(br .le. 0.42432/0.9766) then
             nsta = 130
           ELSEIF(br .le. 0.42932/0.9766) then
             nsta = 131
           ELSEIF(br .le. 0.42982/0.9766) then
             nsta = 131
           ELSEIF(br .le. 0.44082/0.9766) then
             nsta = 132
           ELSEIF(br .le. 0.44252/0.9766) then
             nsta = 132
           ELSEIF(br .le. 0.4429/0.9766) then
             nsta = 133
           ELSEIF(br .le. 0.4456/0.9766) then
             nsta = 133
           ELSEIF(br .le. 0.44648/0.9766) then
             nsta = 134
           ELSEIF(br .le. 0.45348/0.9766) then
             nsta = 135
           ELSEIF(br .le. 0.45436/0.9766) then
             nsta = 135
           ELSEIF(br .le. 0.45584/0.9766) then
             nsta = 137
           ELSEIF(br .le. 0.45695/0.9766) then
             nsta = 137
           ELSEIF(br .le. 0.45726/0.9766) then
             nsta = 137
           ELSEIF(br .le. 0.48826/0.9766) then
             nsta = 138
           ELSEIF(br .le. 0.48906/0.9766) then
             nsta = 139
           ELSEIF(br .le. 0.49336/0.9766) then
             nsta = 140
           ELSEIF(br .le. 0.49496/0.9766) then
             nsta = 141
           ELSEIF(br .le. 0.50396/0.9766) then
             nsta = 142
           ELSEIF(br .le. 0.50646/0.9766) then
             nsta = 143
           ELSEIF(br .le. 0.50711/0.9766) then
             nsta = 144
           ELSEIF(br .le. 0.50757/0.9766) then
             nsta = 145
           ELSEIF(br .le. 0.5096/0.9766) then
             nsta = 146
           ELSEIF(br .le. 0.51005/0.9766) then
             nsta = 147
           ELSEIF(br .le. 0.5119/0.9766) then
             nsta = 147
           ELSEIF(br .le. 0.5129/0.9766) then
             nsta = 147
           ELSEIF(br .le. 0.5219/0.9766) then
             nsta = 148
           ELSEIF(br .le. 0.5242/0.9766) then
             nsta = 149
           ELSEIF(br .le. 0.5562/0.9766) then
             nsta = 150
           ELSEIF(br .le. 0.55838/0.9766) then
             nsta = 151
           ELSEIF(br .le. 0.5602/0.9766) then
             nsta = 152
           ELSEIF(br .le. 0.56074/0.9766) then
             nsta = 152
           ELSEIF(br .le. 0.56247/0.9766) then
             nsta = 153
           ELSEIF(br .le. 0.56468/0.9766) then
             nsta = 154
           ELSEIF(br .le. 0.56506/0.9766) then
             nsta = 154
           ELSEIF(br .le. 0.56826/0.9766) then
             nsta = 155
           ELSEIF(br .le. 0.56972/0.9766) then
             nsta = 155
           ELSEIF(br .le. 0.57042/0.9766) then
             nsta = 156
           ELSEIF(br .le. 0.57115/0.9766) then
             nsta = 156
           ELSEIF(br .le. 0.57333/0.9766) then
             nsta = 156
           ELSEIF(br .le. 0.57673/0.9766) then
             nsta = 157
           ELSEIF(br .le. 0.60073/0.9766) then
             nsta = 160
           ELSEIF(br .le. 0.61573/0.9766) then
             nsta = 161
           ELSEIF(br .le. 0.61718/0.9766) then
             nsta = 162
           ELSEIF(br .le. 0.61768/0.9766) then
             nsta = 163
           ELSEIF(br .le. 0.61837/0.9766) then
             nsta = 165
           ELSEIF(br .le. 0.62077/0.9766) then
             nsta = 165
           ELSEIF(br .le. 0.62123/0.9766) then
             nsta = 166
           ELSEIF(br .le. 0.62286/0.9766) then
             nsta = 167
           ELSEIF(br .le. 0.62886/0.9766) then
             nsta = 167
           ELSEIF(br .le. 0.63046/0.9766) then
             nsta = 168
           ELSEIF(br .le. 0.63436/0.9766) then
             nsta = 168
           ELSEIF(br .le. 0.65836/0.9766) then
             nsta = 170
           ELSEIF(br .le. 0.66446/0.9766) then
             nsta = 170
           ELSEIF(br .le. 0.66481/0.9766) then
             nsta = 171
           ELSEIF(br .le. 0.67481/0.9766) then
             nsta = 172
           ELSEIF(br .le. 0.67596/0.9766) then
             nsta = 172
           ELSEIF(br .le. 0.67707/0.9766) then
             nsta = 174
           ELSEIF(br .le. 0.67887/0.9766) then
             nsta = 174
           ELSEIF(br .le. 0.67937/0.9766) then
             nsta = 175
           ELSEIF(br .le. 0.68347/0.9766) then
             nsta = 177
           ELSEIF(br .le. 0.72047/0.9766) then
             nsta = 178
           ELSEIF(br .le. 0.73047/0.9766) then
             nsta = 179
           ELSEIF(br .le. 0.73547/0.9766) then
             nsta = 180
           ELSEIF(br .le. 0.73686/0.9766) then
             nsta = 181
           ELSEIF(br .le. 0.73886/0.9766) then
             nsta = 182
           ELSEIF(br .le. 0.74086/0.9766) then
             nsta = 182
           ELSEIF(br .le. 0.74171/0.9766) then
             nsta = 183
           ELSEIF(br .le. 0.74259/0.9766) then
             nsta = 184
           ELSEIF(br .le. 0.74989/0.9766) then
             nsta = 184
           ELSEIF(br .le. 0.75069/0.9766) then
             nsta = 185
           ELSEIF(br .le. 0.75669/0.9766) then
             nsta = 186
           ELSEIF(br .le. 0.76029/0.9766) then
             nsta = 186
           ELSEIF(br .le. 0.76106/0.9766) then
             nsta = 186
           ELSEIF(br .le. 0.76263/0.9766) then
             nsta = 187
           ELSEIF(br .le. 0.76301/0.9766) then
             nsta = 187
           ELSEIF(br .le. 0.76408/0.9766) then
             nsta = 187
           ELSEIF(br .le. 0.76818/0.9766) then
             nsta = 187
           ELSEIF(br .le. 0.76958/0.9766) then
             nsta = 188
           ELSEIF(br .le. 0.77108/0.9766) then
             nsta = 188
           ELSEIF(br .le. 0.77368/0.9766) then
             nsta = 189
           ELSEIF(br .le. 0.77588/0.9766) then
             nsta = 190
           ELSEIF(br .le. 0.77798/0.9766) then
             nsta = 192
           ELSEIF(br .le. 0.77817/0.9766) then
             nsta = 192
           ELSEIF(br .le. 0.78417/0.9766) then
             nsta = 193
           ELSEIF(br .le. 0.78647/0.9766) then
             nsta = 193
           ELSEIF(br .le. 0.78927/0.9766) then
             nsta = 194
           ELSEIF(br .le. 0.79337/0.9766) then
             nsta = 194
           ELSEIF(br .le. 0.79477/0.9766) then
             nsta = 195
           ELSEIF(br .le. 0.79777/0.9766) then
             nsta = 195
           ELSEIF(br .le. 0.79865/0.9766) then
             nsta = 196
           ELSEIF(br .le. 0.80325/0.9766) then
             nsta = 196
           ELSEIF(br .le. 0.81725/0.9766) then
             nsta = 197
           ELSEIF(br .le. 0.81825/0.9766) then
             nsta = 197
           ELSEIF(br .le. 0.82145/0.9766) then
             nsta = 198
           ELSEIF(br .le. 0.82336/0.9766) then
             nsta = 198
           ELSEIF(br .le. 0.83336/0.9766) then
             nsta = 198
           ELSEIF(br .le. 0.83451/0.9766) then
             nsta = 198
           ELSEIF(br .le. 0.83539/0.9766) then
             nsta = 198
           ELSEIF(br .le. 0.83769/0.9766) then
             nsta = 198
           ELSEIF(br .le. 0.83811/0.9766) then
             nsta = 199
           ELSEIF(br .le. 0.83893/0.9766) then
             nsta = 199
           ELSEIF(br .le. 0.84473/0.9766) then
             nsta = 199
           ELSEIF(br .le. 0.84553/0.9766) then
             nsta = 199
           ELSEIF(br .le. 0.84637/0.9766) then
             nsta = 200
           ELSEIF(br .le. 0.84668/0.9766) then
             nsta = 200
           ELSEIF(br .le. 0.84748/0.9766) then
             nsta = 200
           ELSEIF(br .le. 0.85058/0.9766) then
             nsta = 200
           ELSEIF(br .le. 0.85138/0.9766) then
             nsta = 201
           ELSEIF(br .le. 0.85253/0.9766) then
             nsta = 201
           ELSEIF(br .le. 0.85403/0.9766) then
             nsta = 202
           ELSEIF(br .le. 0.85813/0.9766) then
             nsta = 202
           ELSEIF(br .le. 0.86713/0.9766) then
             nsta = 203
           ELSEIF(br .le. 0.86893/0.9766) then
             nsta = 203
           ELSEIF(br .le. 0.87253/0.9766) then
             nsta = 204
           ELSEIF(br .le. 0.87295/0.9766) then
             nsta = 204
           ELSEIF(br .le. 0.87665/0.9766) then
             nsta = 204
           ELSEIF(br .le. 0.87703/0.9766) then
             nsta = 205
           ELSEIF(br .le. 0.88003/0.9766) then
             nsta = 205
           ELSEIF(br .le. 0.88273/0.9766) then
             nsta = 205
           ELSEIF(br .le. 0.88473/0.9766) then
             nsta = 206
           ELSEIF(br .le. 0.88853/0.9766) then
             nsta = 207
           ELSEIF(br .le. 0.89083/0.9766) then
             nsta = 207
           ELSEIF(br .le. 0.89363/0.9766) then
             nsta = 208
           ELSEIF(br .le. 0.89478/0.9766) then
             nsta = 208
           ELSEIF(br .le. 0.8966/0.9766) then
             nsta = 208
           ELSEIF(br .le. 0.9096/0.9766) then
             nsta = 208
           ELSEIF(br .le. 0.9133/0.9766) then
             nsta = 208
           ELSEIF(br .le. 0.9148/0.9766) then
             nsta = 208
           ELSEIF(br .le. 0.9166/0.9766) then
             nsta = 208
           ELSEIF(br .le. 0.9185/0.9766) then
             nsta = 208
           ELSEIF(br .le. 0.9245/0.9766) then
             nsta = 208
           ELSEIF(br .le. 0.925/0.9766) then
             nsta = 208
           ELSEIF(br .le. 0.927/0.9766) then
             nsta = 208
           ELSEIF(br .le. 0.9286/0.9766) then
             nsta = 208
           ELSEIF(br .le. 0.92945/0.9766) then
             nsta = 208
           ELSEIF(br .le. 0.93042/0.9766) then
             nsta = 209
           ELSEIF(br .le. 0.93139/0.9766) then
             nsta = 209
           ELSEIF(br .le. 0.93299/0.9766) then
             nsta = 209
           ELSEIF(br .le. 0.93699/0.9766) then
             nsta = 209
           ELSEIF(br .le. 0.94099/0.9766) then
             nsta = 209
           ELSEIF(br .le. 0.94281/0.9766) then
             nsta = 209
           ELSEIF(br .le. 0.94346/0.9766) then
             nsta = 209
           ELSEIF(br .le. 0.94676/0.9766) then
             nsta = 209
           ELSEIF(br .le. 0.94767/0.9766) then
             nsta = 209
           ELSEIF(br .le. 0.95077/0.9766) then
             nsta = 209
           ELSEIF(br .le. 0.95577/0.9766) then
             nsta = 209
           ELSEIF(br .le. 0.95817/0.9766) then
             nsta = 209
           ELSEIF(br .le. 0.96097/0.9766) then
             nsta = 209
           ELSEIF(br .le. 0.97497/0.9766) then
             nsta = 209
           ELSE
             nsta = 209
           ENDIF
         Case(54130)
               IF(br .le. 0.0005/0.5141) then
             nsta = 0
           ELSEIF(br .le. 0.0017/0.5141) then
             nsta = 1
           ELSEIF(br .le. 0.0091/0.5141) then
             nsta = 2
           ELSEIF(br .le. 0.014/0.5141) then
             nsta = 6
           ELSEIF(br .le. 0.0148/0.5141) then
             nsta = 7
           ELSEIF(br .le. 0.0207/0.5141) then
             nsta = 16
           ELSEIF(br .le. 0.0293/0.5141) then
             nsta = 23
           ELSEIF(br .le. 0.0321/0.5141) then
             nsta = 29
           ELSEIF(br .le. 0.034/0.5141) then
             nsta = 35
           ELSEIF(br .le. 0.0407/0.5141) then
             nsta = 37
           ELSEIF(br .le. 0.0533/0.5141) then
             nsta = 40
           ELSEIF(br .le. 0.0603/0.5141) then
             nsta = 45
           ELSEIF(br .le. 0.0659/0.5141) then
             nsta = 47
           ELSEIF(br .le. 0.0929/0.5141) then
             nsta = 49
           ELSEIF(br .le. 0.1239/0.5141) then
             nsta = 52
           ELSEIF(br .le. 0.1599/0.5141) then
             nsta = 53
           ELSEIF(br .le. 0.1769/0.5141) then
             nsta = 54
           ELSEIF(br .le. 0.3369/0.5141) then
             nsta = 58
           ELSEIF(br .le. 0.3469/0.5141) then
             nsta = 61
           ELSEIF(br .le. 0.3522/0.5141) then
             nsta = 62
           ELSEIF(br .le. 0.3762/0.5141) then
             nsta = 63
           ELSEIF(br .le. 0.3802/0.5141) then
             nsta = 65
           ELSEIF(br .le. 0.3901/0.5141) then
             nsta = 69
           ELSEIF(br .le. 0.4045/0.5141) then
             nsta = 70
           ELSEIF(br .le. 0.4203/0.5141) then
             nsta = 70
           ELSEIF(br .le. 0.4248/0.5141) then
             nsta = 75
           ELSEIF(br .le. 0.4468/0.5141) then
             nsta = 75
           ELSEIF(br .le. 0.4908/0.5141) then
             nsta = 75
           ELSEIF(br .le. 0.4948/0.5141) then
             nsta = 77
           ELSEIF(br .le. 0.4988/0.5141) then
             nsta = 77
           ELSEIF(br .le. 0.5071/0.5141) then
             nsta = 77
           ELSE
             nsta = 78
           ENDIF
         Case(54132)
               IF(br .le. 1.33/2.2584) then
             nsta = 20
           ELSEIF(br .le. 1.54/2.2584) then
             nsta = 23
           ELSEIF(br .le. 1.57/2.2584) then
             nsta = 29
           ELSEIF(br .le. 1.635/2.2584) then
             nsta = 59
           ELSEIF(br .le. 1.645/2.2584) then
             nsta = 60
           ELSEIF(br .le. 1.674/2.2584) then
             nsta = 62
           ELSEIF(br .le. 1.78/2.2584) then
             nsta = 63
           ELSEIF(br .le. 1.7885/2.2584) then
             nsta = 64
           ELSEIF(br .le. 1.8405/2.2584) then
             nsta = 66
           ELSEIF(br .le. 1.8645/2.2584) then
             nsta = 67
           ELSEIF(br .le. 1.9215/2.2584) then
             nsta = 69
           ELSEIF(br .le. 2.0285/2.2584) then
             nsta = 70
           ELSEIF(br .le. 2.0755/2.2584) then
             nsta = 71
           ELSEIF(br .le. 2.0864/2.2584) then
             nsta = 72
           ELSEIF(br .le. 2.1444/2.2584) then
             nsta = 73
           ELSEIF(br .le. 2.2154/2.2584) then
             nsta = 75
           ELSE
             nsta = 76
           ENDIF
         Case(54137)
               IF(br .le. 0.0048/0.009089) then
             nsta = 1
           ELSEIF(br .le. 0.0063/0.009089) then
             nsta = 2
           ELSEIF(br .le. 0.006345/0.009089) then
             nsta = 6
           ELSEIF(br .le. 0.007545/0.009089) then
             nsta = 17
           ELSEIF(br .le. 0.008095/0.009089) then
             nsta = 24
           ELSEIF(br .le. 0.008255/0.009089) then
             nsta = 38
           ELSEIF(br .le. 0.008299/0.009089) then
             nsta = 50
           ELSEIF(br .le. 0.008989/0.009089) then
             nsta = 52
           ELSE
             nsta = 55
           ENDIF
         Case(55134)
               IF(br .le. 0.008/3.584) then
             nsta = 0
           ELSEIF(br .le. 0.0135/3.584) then
             nsta = 1
           ELSEIF(br .le. 0.0485/3.584) then
             nsta = 2
           ELSEIF(br .le. 0.0555/3.584) then
             nsta = 4
           ELSEIF(br .le. 0.1455/3.584) then
             nsta = 6
           ELSEIF(br .le. 0.1535/3.584) then
             nsta = 7
           ELSEIF(br .le. 0.3775/3.584) then
             nsta = 8
           ELSEIF(br .le. 0.3835/3.584) then
             nsta = 9
           ELSEIF(br .le. 0.3848/3.584) then
             nsta = 10
           ELSEIF(br .le. 0.3928/3.584) then
             nsta = 11
           ELSEIF(br .le. 0.3968/3.584) then
             nsta = 14
           ELSEIF(br .le. 0.3972/3.584) then
             nsta = 15
           ELSEIF(br .le. 0.4412/3.584) then
             nsta = 17
           ELSEIF(br .le. 0.5232/3.584) then
             nsta = 21
           ELSEIF(br .le. 0.5262/3.584) then
             nsta = 22
           ELSEIF(br .le. 0.5267/3.584) then
             nsta = 23
           ELSEIF(br .le. 0.5377/3.584) then
             nsta = 24
           ELSEIF(br .le. 0.5412/3.584) then
             nsta = 25
           ELSEIF(br .le. 0.5426/3.584) then
             nsta = 26
           ELSEIF(br .le. 0.5926/3.584) then
             nsta = 27
           ELSEIF(br .le. 0.5966/3.584) then
             nsta = 28
           ELSEIF(br .le. 0.6406/3.584) then
             nsta = 29
           ELSEIF(br .le. 0.6414/3.584) then
             nsta = 30
           ELSEIF(br .le. 0.6442/3.584) then
             nsta = 31
           ELSEIF(br .le. 0.6448/3.584) then
             nsta = 32
           ELSEIF(br .le. 0.6828/3.584) then
             nsta = 34
           ELSEIF(br .le. 0.684/3.584) then
             nsta = 39
           ELSEIF(br .le. 0.691/3.584) then
             nsta = 40
           ELSEIF(br .le. 0.726/3.584) then
             nsta = 41
           ELSEIF(br .le. 0.917/3.584) then
             nsta = 42
           ELSEIF(br .le. 1.169/3.584) then
             nsta = 43
           ELSEIF(br .le. 1.207/3.584) then
             nsta = 45
           ELSEIF(br .le. 1.268/3.584) then
             nsta = 47
           ELSEIF(br .le. 1.2695/3.584) then
             nsta = 49
           ELSEIF(br .le. 1.2727/3.584) then
             nsta = 53
           ELSEIF(br .le. 1.2787/3.584) then
             nsta = 54
           ELSEIF(br .le. 1.2887/3.584) then
             nsta = 55
           ELSEIF(br .le. 1.5287/3.584) then
             nsta = 56
           ELSEIF(br .le. 1.5291/3.584) then
             nsta = 61
           ELSEIF(br .le. 1.6281/3.584) then
             nsta = 62
           ELSEIF(br .le. 1.6551/3.584) then
             nsta = 63
           ELSEIF(br .le. 1.6591/3.584) then
             nsta = 65
           ELSEIF(br .le. 1.7041/3.584) then
             nsta = 66
           ELSEIF(br .le. 1.7061/3.584) then
             nsta = 67
           ELSEIF(br .le. 1.7531/3.584) then
             nsta = 68
           ELSEIF(br .le. 1.8691/3.584) then
             nsta = 70
           ELSEIF(br .le. 1.873/3.584) then
             nsta = 75
           ELSEIF(br .le. 1.993/3.584) then
             nsta = 77
           ELSEIF(br .le. 2.006/3.584) then
             nsta = 78
           ELSEIF(br .le. 2.143/3.584) then
             nsta = 79
           ELSEIF(br .le. 2.289/3.584) then
             nsta = 83
           ELSEIF(br .le. 2.376/3.584) then
             nsta = 85
           ELSEIF(br .le. 2.384/3.584) then
             nsta = 93
           ELSEIF(br .le. 2.661/3.584) then
             nsta = 95
           ELSEIF(br .le. 2.772/3.584) then
             nsta = 97
           ELSEIF(br .le. 3.021/3.584) then
             nsta = 104
           ELSEIF(br .le. 3.354/3.584) then
             nsta = 111
           ELSE
             nsta = 113
           ENDIF
         Case(56135)
               IF(br .le. 0.000173/0.011032) then
             nsta = 0
           ELSEIF(br .le. 0.000251/0.011032) then
             nsta = 4
           ELSEIF(br .le. 0.000294/0.011032) then
             nsta = 16
           ELSEIF(br .le. 0.000363/0.011032) then
             nsta = 18
           ELSEIF(br .le. 0.002023/0.011032) then
             nsta = 21
           ELSEIF(br .le. 0.002084/0.011032) then
             nsta = 22
           ELSEIF(br .le. 0.002364/0.011032) then
             nsta = 23
           ELSEIF(br .le. 0.002496/0.011032) then
             nsta = 24
           ELSEIF(br .le. 0.002586/0.011032) then
             nsta = 32
           ELSEIF(br .le. 0.002756/0.011032) then
             nsta = 35
           ELSEIF(br .le. 0.003196/0.011032) then
             nsta = 38
           ELSEIF(br .le. 0.003237/0.011032) then
             nsta = 46
           ELSEIF(br .le. 0.003887/0.011032) then
             nsta = 49
           ELSEIF(br .le. 0.003987/0.011032) then
             nsta = 53
           ELSEIF(br .le. 0.004187/0.011032) then
             nsta = 58
           ELSEIF(br .le. 0.004337/0.011032) then
             nsta = 60
           ELSEIF(br .le. 0.004487/0.011032) then
             nsta = 62
           ELSEIF(br .le. 0.005027/0.011032) then
             nsta = 63
           ELSEIF(br .le. 0.005217/0.011032) then
             nsta = 68
           ELSEIF(br .le. 0.005287/0.011032) then
             nsta = 69
           ELSEIF(br .le. 0.005507/0.011032) then
             nsta = 70
           ELSEIF(br .le. 0.005667/0.011032) then
             nsta = 73
           ELSEIF(br .le. 0.005787/0.011032) then
             nsta = 74
           ELSEIF(br .le. 0.006067/0.011032) then
             nsta = 79
           ELSEIF(br .le. 0.006132/0.011032) then
             nsta = 81
           ELSEIF(br .le. 0.010632/0.011032) then
             nsta = 88
           ELSEIF(br .le. 0.010767/0.011032) then
             nsta = 95
           ELSEIF(br .le. 0.010793/0.011032) then
             nsta = 97
           ELSEIF(br .le. 0.010826/0.011032) then
             nsta = 107
           ELSEIF(br .le. 0.010956/0.011032) then
             nsta = 111
           ELSE
             nsta = 112
           ENDIF
         Case(56136)
               IF(br .le. 0.00635/0.04564) then
             nsta = 0
           ELSEIF(br .le. 0.00984/0.04564) then
             nsta = 1
           ELSEIF(br .le. 0.01104/0.04564) then
             nsta = 2
           ELSEIF(br .le. 0.01404/0.04564) then
             nsta = 3
           ELSEIF(br .le. 0.0143/0.04564) then
             nsta = 6
           ELSEIF(br .le. 0.0158/0.04564) then
             nsta = 7
           ELSEIF(br .le. 0.01635/0.04564) then
             nsta = 8
           ELSEIF(br .le. 0.01665/0.04564) then
             nsta = 10
           ELSEIF(br .le. 0.01925/0.04564) then
             nsta = 13
           ELSEIF(br .le. 0.01945/0.04564) then
             nsta = 15
           ELSEIF(br .le. 0.01968/0.04564) then
             nsta = 19
           ELSEIF(br .le. 0.02077/0.04564) then
             nsta = 21
           ELSEIF(br .le. 0.02122/0.04564) then
             nsta = 22
           ELSEIF(br .le. 0.02273/0.04564) then
             nsta = 24
           ELSEIF(br .le. 0.02347/0.04564) then
             nsta = 28
           ELSEIF(br .le. 0.02436/0.04564) then
             nsta = 30
           ELSEIF(br .le. 0.02462/0.04564) then
             nsta = 31
           ELSEIF(br .le. 0.02486/0.04564) then
             nsta = 33
           ELSEIF(br .le. 0.02508/0.04564) then
             nsta = 36
           ELSEIF(br .le. 0.0253/0.04564) then
             nsta = 44
           ELSEIF(br .le. 0.02626/0.04564) then
             nsta = 46
           ELSEIF(br .le. 0.03142/0.04564) then
             nsta = 47
           ELSEIF(br .le. 0.03242/0.04564) then
             nsta = 50
           ELSEIF(br .le. 0.03413/0.04564) then
             nsta = 59
           ELSEIF(br .le. 0.0353/0.04564) then
             nsta = 72
           ELSEIF(br .le. 0.0435/0.04564) then
             nsta = 73
           ELSEIF(br .le. 0.04475/0.04564) then
             nsta = 81
           ELSE
             nsta = 82
           ENDIF
         Case(56137)
               IF(br .le. 0.0001/0.048184) then
             nsta = 0
           ELSEIF(br .le. 0.0035/0.048184) then
             nsta = 1
           ELSEIF(br .le. 0.00376/0.048184) then
             nsta = 13
           ELSEIF(br .le. 0.003935/0.048184) then
             nsta = 14
           ELSEIF(br .le. 0.030335/0.048184) then
             nsta = 18
           ELSEIF(br .le. 0.030375/0.048184) then
             nsta = 26
           ELSEIF(br .le. 0.039075/0.048184) then
             nsta = 36
           ELSEIF(br .le. 0.039295/0.048184) then
             nsta = 44
           ELSEIF(br .le. 0.039605/0.048184) then
             nsta = 49
           ELSEIF(br .le. 0.040455/0.048184) then
             nsta = 51
           ELSEIF(br .le. 0.040695/0.048184) then
             nsta = 55
           ELSEIF(br .le. 0.042215/0.048184) then
             nsta = 60
           ELSEIF(br .le. 0.043465/0.048184) then
             nsta = 63
           ELSEIF(br .le. 0.043555/0.048184) then
             nsta = 72
           ELSEIF(br .le. 0.044485/0.048184) then
             nsta = 74
           ELSEIF(br .le. 0.044695/0.048184) then
             nsta = 76
           ELSEIF(br .le. 0.045235/0.048184) then
             nsta = 79
           ELSEIF(br .le. 0.045284/0.048184) then
             nsta = 83
           ELSEIF(br .le. 0.045354/0.048184) then
             nsta = 83
           ELSEIF(br .le. 0.045484/0.048184) then
             nsta = 87
           ELSEIF(br .le. 0.046024/0.048184) then
             nsta = 88
           ELSEIF(br .le. 0.046774/0.048184) then
             nsta = 90
           ELSEIF(br .le. 0.047194/0.048184) then
             nsta = 91
           ELSEIF(br .le. 0.047414/0.048184) then
             nsta = 92
           ELSEIF(br .le. 0.047564/0.048184) then
             nsta = 93
           ELSEIF(br .le. 0.047974/0.048184) then
             nsta = 94
           ELSE
             nsta = 95
           ENDIF
         Case(56138)
               IF(br .le. 0.00032/0.22984) then
             nsta = 0
           ELSEIF(br .le. 0.00134/0.22984) then
             nsta = 1
           ELSEIF(br .le. 0.00471/0.22984) then
             nsta = 4
           ELSEIF(br .le. 0.00523/0.22984) then
             nsta = 6
           ELSEIF(br .le. 0.00545/0.22984) then
             nsta = 7
           ELSEIF(br .le. 0.00623/0.22984) then
             nsta = 10
           ELSEIF(br .le. 0.01553/0.22984) then
             nsta = 12
           ELSEIF(br .le. 0.01993/0.22984) then
             nsta = 13
           ELSEIF(br .le. 0.0203/0.22984) then
             nsta = 14
           ELSEIF(br .le. 0.02054/0.22984) then
             nsta = 15
           ELSEIF(br .le. 0.08224/0.22984) then
             nsta = 17
           ELSEIF(br .le. 0.08294/0.22984) then
             nsta = 19
           ELSEIF(br .le. 0.08353/0.22984) then
             nsta = 20
           ELSEIF(br .le. 0.08516/0.22984) then
             nsta = 21
           ELSEIF(br .le. 0.08812/0.22984) then
             nsta = 22
           ELSEIF(br .le. 0.09342/0.22984) then
             nsta = 24
           ELSEIF(br .le. 0.09384/0.22984) then
             nsta = 26
           ELSEIF(br .le. 0.09426/0.22984) then
             nsta = 27
           ELSEIF(br .le. 0.10306/0.22984) then
             nsta = 29
           ELSEIF(br .le. 0.10425/0.22984) then
             nsta = 32
           ELSEIF(br .le. 0.10517/0.22984) then
             nsta = 33
           ELSEIF(br .le. 0.10631/0.22984) then
             nsta = 35
           ELSEIF(br .le. 0.11231/0.22984) then
             nsta = 37
           ELSEIF(br .le. 0.11446/0.22984) then
             nsta = 40
           ELSEIF(br .le. 0.12426/0.22984) then
             nsta = 45
           ELSEIF(br .le. 0.12503/0.22984) then
             nsta = 49
           ELSEIF(br .le. 0.13133/0.22984) then
             nsta = 50
           ELSEIF(br .le. 0.14533/0.22984) then
             nsta = 53
           ELSEIF(br .le. 0.14699/0.22984) then
             nsta = 54
           ELSEIF(br .le. 0.14939/0.22984) then
             nsta = 55
           ELSEIF(br .le. 0.15379/0.22984) then
             nsta = 57
           ELSEIF(br .le. 0.1546/0.22984) then
             nsta = 58
           ELSEIF(br .le. 0.1574/0.22984) then
             nsta = 60
           ELSEIF(br .le. 0.16/0.22984) then
             nsta = 62
           ELSEIF(br .le. 0.1604/0.22984) then
             nsta = 64
           ELSEIF(br .le. 0.16091/0.22984) then
             nsta = 65
           ELSEIF(br .le. 0.16221/0.22984) then
             nsta = 68
           ELSEIF(br .le. 0.16911/0.22984) then
             nsta = 69
           ELSEIF(br .le. 0.17811/0.22984) then
             nsta = 70
           ELSEIF(br .le. 0.18401/0.22984) then
             nsta = 71
           ELSEIF(br .le. 0.18781/0.22984) then
             nsta = 72
           ELSEIF(br .le. 0.19351/0.22984) then
             nsta = 73
           ELSEIF(br .le. 0.19871/0.22984) then
             nsta = 74
           ELSEIF(br .le. 0.20191/0.22984) then
             nsta = 75
           ELSEIF(br .le. 0.20441/0.22984) then
             nsta = 76
           ELSEIF(br .le. 0.20589/0.22984) then
             nsta = 77
           ELSEIF(br .le. 0.20819/0.22984) then
             nsta = 78
           ELSEIF(br .le. 0.21199/0.22984) then
             nsta = 80
           ELSEIF(br .le. 0.21541/0.22984) then
             nsta = 82
           ELSEIF(br .le. 0.21831/0.22984) then
             nsta = 83
           ELSEIF(br .le. 0.21964/0.22984) then
             nsta = 86
           ELSEIF(br .le. 0.22224/0.22984) then
             nsta = 87
           ELSEIF(br .le. 0.22314/0.22984) then
             nsta = 88
           ELSEIF(br .le. 0.22734/0.22984) then
             nsta = 91
           ELSE
             nsta = 92
           ENDIF
         Case(56139)
               IF(br .le. 0.155/0.2622) then
             nsta = 1
           ELSEIF(br .le. 0.2112/0.2622) then
             nsta = 2
           ELSEIF(br .le. 0.2127/0.2622) then
             nsta = 3
           ELSEIF(br .le. 0.2314/0.2622) then
             nsta = 31
           ELSEIF(br .le. 0.2404/0.2622) then
             nsta = 32
           ELSEIF(br .le. 0.2506/0.2622) then
             nsta = 36
           ELSE
             nsta = 46
           ENDIF
         Case(57140)
               IF(br .le. 0.089/5.75656) then
             nsta = 0
           ELSEIF(br .le. 0.1049/5.75656) then
             nsta = 1
           ELSEIF(br .le. 0.2189/5.75656) then
             nsta = 2
           ELSEIF(br .le. 0.8989/5.75656) then
             nsta = 5
           ELSEIF(br .le. 0.9134/5.75656) then
             nsta = 9
           ELSEIF(br .le. 1.0634/5.75656) then
             nsta = 10
           ELSEIF(br .le. 1.7244/5.75656) then
             nsta = 12
           ELSEIF(br .le. 1.7732/5.75656) then
             nsta = 19
           ELSEIF(br .le. 1.9372/5.75656) then
             nsta = 20
           ELSEIF(br .le. 1.9447/5.75656) then
             nsta = 23
           ELSEIF(br .le. 2.1917/5.75656) then
             nsta = 24
           ELSEIF(br .le. 2.4467/5.75656) then
             nsta = 26
           ELSEIF(br .le. 2.45103/5.75656) then
             nsta = 28
           ELSEIF(br .le. 2.45209/5.75656) then
             nsta = 31
           ELSEIF(br .le. 2.45237/5.75656) then
             nsta = 36
           ELSEIF(br .le. 2.47067/5.75656) then
             nsta = 39
           ELSEIF(br .le. 2.49447/5.75656) then
             nsta = 41
           ELSEIF(br .le. 2.49557/5.75656) then
             nsta = 44
           ELSEIF(br .le. 2.52527/5.75656) then
             nsta = 45
           ELSEIF(br .le. 2.55497/5.75656) then
             nsta = 46
           ELSEIF(br .le. 2.55596/5.75656) then
             nsta = 48
           ELSEIF(br .le. 2.55816/5.75656) then
             nsta = 49
           ELSEIF(br .le. 2.57016/5.75656) then
             nsta = 51
           ELSEIF(br .le. 2.57186/5.75656) then
             nsta = 52
           ELSEIF(br .le. 2.57441/5.75656) then
             nsta = 54
           ELSEIF(br .le. 2.59421/5.75656) then
             nsta = 55
           ELSEIF(br .le. 2.59481/5.75656) then
             nsta = 58
           ELSEIF(br .le. 2.64791/5.75656) then
             nsta = 59
           ELSEIF(br .le. 2.65351/5.75656) then
             nsta = 62
           ELSEIF(br .le. 2.6564/5.75656) then
             nsta = 64
           ELSEIF(br .le. 2.6695/5.75656) then
             nsta = 67
           ELSEIF(br .le. 2.6742/5.75656) then
             nsta = 68
           ELSEIF(br .le. 2.6976/5.75656) then
             nsta = 73
           ELSEIF(br .le. 2.7328/5.75656) then
             nsta = 75
           ELSEIF(br .le. 2.7498/5.75656) then
             nsta = 77
           ELSEIF(br .le. 2.8228/5.75656) then
             nsta = 80
           ELSEIF(br .le. 2.8612/5.75656) then
             nsta = 82
           ELSEIF(br .le. 2.8962/5.75656) then
             nsta = 86
           ELSEIF(br .le. 2.9284/5.75656) then
             nsta = 87
           ELSEIF(br .le. 3.0674/5.75656) then
             nsta = 89
           ELSEIF(br .le. 3.2024/5.75656) then
             nsta = 94
           ELSEIF(br .le. 3.20396/5.75656) then
             nsta = 98
           ELSEIF(br .le. 3.20766/5.75656) then
             nsta = 99
           ELSEIF(br .le. 3.21316/5.75656) then
             nsta = 100
           ELSEIF(br .le. 3.26796/5.75656) then
             nsta = 101
           ELSEIF(br .le. 3.32356/5.75656) then
             nsta = 102
           ELSEIF(br .le. 3.33926/5.75656) then
             nsta = 103
           ELSEIF(br .le. 3.33968/5.75656) then
             nsta = 105
           ELSEIF(br .le. 3.35258/5.75656) then
             nsta = 106
           ELSEIF(br .le. 3.36558/5.75656) then
             nsta = 108
           ELSEIF(br .le. 3.36978/5.75656) then
             nsta = 110
           ELSEIF(br .le. 3.37311/5.75656) then
             nsta = 111
           ELSEIF(br .le. 3.37951/5.75656) then
             nsta = 113
           ELSEIF(br .le. 3.38671/5.75656) then
             nsta = 115
           ELSEIF(br .le. 3.38801/5.75656) then
             nsta = 116
           ELSEIF(br .le. 3.39031/5.75656) then
             nsta = 117
           ELSEIF(br .le. 3.39591/5.75656) then
             nsta = 118
           ELSEIF(br .le. 3.41291/5.75656) then
             nsta = 119
           ELSEIF(br .le. 3.41561/5.75656) then
             nsta = 120
           ELSEIF(br .le. 3.46001/5.75656) then
             nsta = 121
           ELSEIF(br .le. 3.46271/5.75656) then
             nsta = 123
           ELSEIF(br .le. 3.48261/5.75656) then
             nsta = 124
           ELSEIF(br .le. 3.48401/5.75656) then
             nsta = 125
           ELSEIF(br .le. 3.52501/5.75656) then
             nsta = 126
           ELSEIF(br .le. 3.54971/5.75656) then
             nsta = 127
           ELSEIF(br .le. 3.60771/5.75656) then
             nsta = 129
           ELSEIF(br .le. 3.63091/5.75656) then
             nsta = 130
           ELSEIF(br .le. 3.64901/5.75656) then
             nsta = 131
           ELSEIF(br .le. 3.65581/5.75656) then
             nsta = 132
           ELSEIF(br .le. 3.67291/5.75656) then
             nsta = 133
           ELSEIF(br .le. 3.68901/5.75656) then
             nsta = 134
           ELSEIF(br .le. 3.69211/5.75656) then
             nsta = 135
           ELSEIF(br .le. 3.71631/5.75656) then
             nsta = 136
           ELSEIF(br .le. 3.71951/5.75656) then
             nsta = 137
           ELSEIF(br .le. 3.72171/5.75656) then
             nsta = 138
           ELSEIF(br .le. 3.73371/5.75656) then
             nsta = 139
           ELSEIF(br .le. 3.73933/5.75656) then
             nsta = 140
           ELSEIF(br .le. 3.74216/5.75656) then
             nsta = 142
           ELSEIF(br .le. 3.75116/5.75656) then
             nsta = 143
           ELSEIF(br .le. 3.75755/5.75656) then
             nsta = 144
           ELSEIF(br .le. 3.75918/5.75656) then
             nsta = 145
           ELSEIF(br .le. 3.79108/5.75656) then
             nsta = 146
           ELSEIF(br .le. 3.79488/5.75656) then
             nsta = 147
           ELSEIF(br .le. 3.79646/5.75656) then
             nsta = 148
           ELSEIF(br .le. 3.80026/5.75656) then
             nsta = 149
           ELSEIF(br .le. 3.80586/5.75656) then
             nsta = 150
           ELSEIF(br .le. 3.80883/5.75656) then
             nsta = 151
           ELSEIF(br .le. 3.81303/5.75656) then
             nsta = 152
           ELSEIF(br .le. 3.86363/5.75656) then
             nsta = 153
           ELSEIF(br .le. 3.91683/5.75656) then
             nsta = 154
           ELSEIF(br .le. 3.92363/5.75656) then
             nsta = 155
           ELSEIF(br .le. 3.95363/5.75656) then
             nsta = 156
           ELSEIF(br .le. 3.96803/5.75656) then
             nsta = 157
           ELSEIF(br .le. 3.9724/5.75656) then
             nsta = 158
           ELSEIF(br .le. 3.9937/5.75656) then
             nsta = 159
           ELSEIF(br .le. 4.0475/5.75656) then
             nsta = 160
           ELSEIF(br .le. 4.061/5.75656) then
             nsta = 161
           ELSEIF(br .le. 4.0683/5.75656) then
             nsta = 162
           ELSEIF(br .le. 4.1007/5.75656) then
             nsta = 163
           ELSEIF(br .le. 4.1097/5.75656) then
             nsta = 164
           ELSEIF(br .le. 4.1417/5.75656) then
             nsta = 165
           ELSEIF(br .le. 4.1656/5.75656) then
             nsta = 166
           ELSEIF(br .le. 4.1689/5.75656) then
             nsta = 167
           ELSEIF(br .le. 4.1807/5.75656) then
             nsta = 168
           ELSEIF(br .le. 4.1983/5.75656) then
             nsta = 169
           ELSEIF(br .le. 4.2303/5.75656) then
             nsta = 170
           ELSEIF(br .le. 4.2494/5.75656) then
             nsta = 171
           ELSEIF(br .le. 4.2608/5.75656) then
             nsta = 172
           ELSEIF(br .le. 4.4008/5.75656) then
             nsta = 173
           ELSEIF(br .le. 4.4138/5.75656) then
             nsta = 174
           ELSEIF(br .le. 4.418/5.75656) then
             nsta = 175
           ELSEIF(br .le. 4.4237/5.75656) then
             nsta = 176
           ELSEIF(br .le. 4.4431/5.75656) then
             nsta = 177
           ELSEIF(br .le. 4.4614/5.75656) then
             nsta = 178
           ELSEIF(br .le. 4.4652/5.75656) then
             nsta = 179
           ELSEIF(br .le. 4.4946/5.75656) then
             nsta = 180
           ELSEIF(br .le. 4.5464/5.75656) then
             nsta = 181
           ELSEIF(br .le. 4.5794/5.75656) then
             nsta = 182
           ELSEIF(br .le. 4.6465/5.75656) then
             nsta = 183
           ELSEIF(br .le. 4.6498/5.75656) then
             nsta = 185
           ELSEIF(br .le. 4.652/5.75656) then
             nsta = 186
           ELSEIF(br .le. 4.6656/5.75656) then
             nsta = 187
           ELSEIF(br .le. 4.7114/5.75656) then
             nsta = 188
           ELSEIF(br .le. 4.7214/5.75656) then
             nsta = 189
           ELSEIF(br .le. 4.7378/5.75656) then
             nsta = 190
           ELSEIF(br .le. 4.7787/5.75656) then
             nsta = 191
           ELSEIF(br .le. 4.7864/5.75656) then
             nsta = 192
           ELSEIF(br .le. 4.8126/5.75656) then
             nsta = 193
           ELSEIF(br .le. 4.81628/5.75656) then
             nsta = 195
           ELSEIF(br .le. 4.81904/5.75656) then
             nsta = 196
           ELSEIF(br .le. 4.82502/5.75656) then
             nsta = 197
           ELSEIF(br .le. 4.86852/5.75656) then
             nsta = 198
           ELSEIF(br .le. 4.87712/5.75656) then
             nsta = 199
           ELSEIF(br .le. 4.89012/5.75656) then
             nsta = 200
           ELSEIF(br .le. 4.90252/5.75656) then
             nsta = 201
           ELSEIF(br .le. 4.91372/5.75656) then
             nsta = 202
           ELSEIF(br .le. 4.92182/5.75656) then
             nsta = 203
           ELSEIF(br .le. 4.92633/5.75656) then
             nsta = 204
           ELSEIF(br .le. 4.93133/5.75656) then
             nsta = 205
           ELSEIF(br .le. 4.94143/5.75656) then
             nsta = 206
           ELSEIF(br .le. 4.94323/5.75656) then
             nsta = 207
           ELSEIF(br .le. 5.01623/5.75656) then
             nsta = 208
           ELSEIF(br .le. 5.03013/5.75656) then
             nsta = 209
           ELSEIF(br .le. 5.03413/5.75656) then
             nsta = 210
           ELSEIF(br .le. 5.04053/5.75656) then
             nsta = 211
           ELSEIF(br .le. 5.06003/5.75656) then
             nsta = 212
           ELSEIF(br .le. 5.06513/5.75656) then
             nsta = 213
           ELSEIF(br .le. 5.07143/5.75656) then
             nsta = 214
           ELSEIF(br .le. 5.07463/5.75656) then
             nsta = 215
           ELSEIF(br .le. 5.08125/5.75656) then
             nsta = 216
           ELSEIF(br .le. 5.10155/5.75656) then
             nsta = 217
           ELSEIF(br .le. 5.11245/5.75656) then
             nsta = 218
           ELSEIF(br .le. 5.1172/5.75656) then
             nsta = 219
           ELSEIF(br .le. 5.1459/5.75656) then
             nsta = 220
           ELSEIF(br .le. 5.1748/5.75656) then
             nsta = 221
           ELSEIF(br .le. 5.2263/5.75656) then
             nsta = 222
           ELSEIF(br .le. 5.2461/5.75656) then
             nsta = 223
           ELSEIF(br .le. 5.2661/5.75656) then
             nsta = 224
           ELSEIF(br .le. 5.2849/5.75656) then
             nsta = 225
           ELSEIF(br .le. 5.28786/5.75656) then
             nsta = 226
           ELSEIF(br .le. 5.30296/5.75656) then
             nsta = 227
           ELSEIF(br .le. 5.31706/5.75656) then
             nsta = 228
           ELSEIF(br .le. 5.32876/5.75656) then
             nsta = 229
           ELSEIF(br .le. 5.33966/5.75656) then
             nsta = 230
           ELSEIF(br .le. 5.35816/5.75656) then
             nsta = 231
           ELSEIF(br .le. 5.36966/5.75656) then
             nsta = 232
           ELSEIF(br .le. 5.39506/5.75656) then
             nsta = 233
           ELSEIF(br .le. 5.40506/5.75656) then
             nsta = 234
           ELSEIF(br .le. 5.42976/5.75656) then
             nsta = 236
           ELSEIF(br .le. 5.45606/5.75656) then
             nsta = 237
           ELSEIF(br .le. 5.47206/5.75656) then
             nsta = 238
           ELSEIF(br .le. 5.48046/5.75656) then
             nsta = 239
           ELSEIF(br .le. 5.48196/5.75656) then
             nsta = 240
           ELSEIF(br .le. 5.49686/5.75656) then
             nsta = 241
           ELSEIF(br .le. 5.50546/5.75656) then
             nsta = 242
           ELSEIF(br .le. 5.53986/5.75656) then
             nsta = 243
           ELSEIF(br .le. 5.56296/5.75656) then
             nsta = 244
           ELSEIF(br .le. 5.56846/5.75656) then
             nsta = 245
           ELSEIF(br .le. 5.60576/5.75656) then
             nsta = 246
           ELSEIF(br .le. 5.63166/5.75656) then
             nsta = 247
           ELSEIF(br .le. 5.65476/5.75656) then
             nsta = 248
           ELSEIF(br .le. 5.66666/5.75656) then
             nsta = 249
           ELSEIF(br .le. 5.68546/5.75656) then
             nsta = 250
           ELSEIF(br .le. 5.72076/5.75656) then
             nsta = 251
           ELSEIF(br .le. 5.74016/5.75656) then
             nsta = 252
           ELSE
             nsta = 253
           ENDIF
         Case(58137)
               IF(br .le. 0.000016/0.0001294) then
             nsta = 0
           ELSEIF(br .le. 0.000047/0.0001294) then
             nsta = 1
           ELSEIF(br .le. 0.0000524/0.0001294) then
             nsta = 6
           ELSEIF(br .le. 0.0000549/0.0001294) then
             nsta = 10
           ELSEIF(br .le. 0.0000629/0.0001294) then
             nsta = 15
           ELSEIF(br .le. 0.0000739/0.0001294) then
             nsta = 21
           ELSEIF(br .le. 0.0000758/0.0001294) then
             nsta = 22
           ELSEIF(br .le. 0.000078/0.0001294) then
             nsta = 26
           ELSEIF(br .le. 0.0000802/0.0001294) then
             nsta = 39
           ELSEIF(br .le. 0.0000882/0.0001294) then
             nsta = 40
           ELSEIF(br .le. 0.0001042/0.0001294) then
             nsta = 45
           ELSEIF(br .le. 0.0001202/0.0001294) then
             nsta = 46
           ELSEIF(br .le. 0.0001259/0.0001294) then
             nsta = 52
           ELSE
             nsta = 58
           ENDIF
         Case(58141)
               IF(br .le. 0.113/0.2208) then
             nsta = 1
           ELSEIF(br .le. 0.166/0.2208) then
             nsta = 2
           ELSEIF(br .le. 0.1755/0.2208) then
             nsta = 11
           ELSEIF(br .le. 0.1794/0.2208) then
             nsta = 17
           ELSEIF(br .le. 0.186/0.2208) then
             nsta = 27
           ELSEIF(br .le. 0.1932/0.2208) then
             nsta = 36
           ELSEIF(br .le. 0.2046/0.2208) then
             nsta = 39
           ELSEIF(br .le. 0.215/0.2208) then
             nsta = 40
           ELSE
             nsta = 45
           ENDIF
         Case(58143)
               IF(br .le. 0.00061/0.04692) then
             nsta = 0
           ELSEIF(br .le. 0.02571/0.04692) then
             nsta = 6
           ELSEIF(br .le. 0.02941/0.04692) then
             nsta = 8
           ELSEIF(br .le. 0.03321/0.04692) then
             nsta = 11
           ELSEIF(br .le. 0.03438/0.04692) then
             nsta = 14
           ELSEIF(br .le. 0.03571/0.04692) then
             nsta = 21
           ELSEIF(br .le. 0.037/0.04692) then
             nsta = 23
           ELSEIF(br .le. 0.0378/0.04692) then
             nsta = 31
           ELSEIF(br .le. 0.03849/0.04692) then
             nsta = 35
           ELSEIF(br .le. 0.04089/0.04692) then
             nsta = 42
           ELSEIF(br .le. 0.04379/0.04692) then
             nsta = 48
           ELSEIF(br .le. 0.04425/0.04692) then
             nsta = 52
           ELSEIF(br .le. 0.04505/0.04692) then
             nsta = 54
           ELSEIF(br .le. 0.04597/0.04692) then
             nsta = 58
           ELSE
             nsta = 61
           ENDIF
         Case(59142)
               IF(br .le. 0.147/3.65423) then
             nsta = 0
           ELSEIF(br .le. 0.187/3.65423) then
             nsta = 2
           ELSEIF(br .le. 0.2241/3.65423) then
             nsta = 4
           ELSEIF(br .le. 0.2259/3.65423) then
             nsta = 5
           ELSEIF(br .le. 0.2376/3.65423) then
             nsta = 8
           ELSEIF(br .le. 0.6166/3.65423) then
             nsta = 10
           ELSEIF(br .le. 0.61988/3.65423) then
             nsta = 11
           ELSEIF(br .le. 0.65288/3.65423) then
             nsta = 14
           ELSEIF(br .le. 0.92188/3.65423) then
             nsta = 17
           ELSEIF(br .le. 1.01988/3.65423) then
             nsta = 18
           ELSEIF(br .le. 1.22788/3.65423) then
             nsta = 19
           ELSEIF(br .le. 1.26078/3.65423) then
             nsta = 21
           ELSEIF(br .le. 1.27428/3.65423) then
             nsta = 22
           ELSEIF(br .le. 1.28548/3.65423) then
             nsta = 25
           ELSEIF(br .le. 1.42548/3.65423) then
             nsta = 26
           ELSEIF(br .le. 1.43048/3.65423) then
             nsta = 27
           ELSEIF(br .le. 1.44538/3.65423) then
             nsta = 28
           ELSEIF(br .le. 1.44647/3.65423) then
             nsta = 29
           ELSEIF(br .le. 1.52947/3.65423) then
             nsta = 30
           ELSEIF(br .le. 1.82047/3.65423) then
             nsta = 31
           ELSEIF(br .le. 1.82244/3.65423) then
             nsta = 32
           ELSEIF(br .le. 1.82386/3.65423) then
             nsta = 33
           ELSEIF(br .le. 1.82856/3.65423) then
             nsta = 34
           ELSEIF(br .le. 1.83776/3.65423) then
             nsta = 35
           ELSEIF(br .le. 1.85426/3.65423) then
             nsta = 39
           ELSEIF(br .le. 1.86686/3.65423) then
             nsta = 40
           ELSEIF(br .le. 1.87316/3.65423) then
             nsta = 41
           ELSEIF(br .le. 1.97116/3.65423) then
             nsta = 42
           ELSEIF(br .le. 1.97476/3.65423) then
             nsta = 44
           ELSEIF(br .le. 1.99756/3.65423) then
             nsta = 45
           ELSEIF(br .le. 2.02276/3.65423) then
             nsta = 46
           ELSEIF(br .le. 2.02516/3.65423) then
             nsta = 49
           ELSEIF(br .le. 2.03006/3.65423) then
             nsta = 50
           ELSEIF(br .le. 2.05696/3.65423) then
             nsta = 51
           ELSEIF(br .le. 2.07356/3.65423) then
             nsta = 52
           ELSEIF(br .le. 2.07816/3.65423) then
             nsta = 53
           ELSEIF(br .le. 2.08106/3.65423) then
             nsta = 54
           ELSEIF(br .le. 2.09346/3.65423) then
             nsta = 55
           ELSEIF(br .le. 2.10226/3.65423) then
             nsta = 56
           ELSEIF(br .le. 2.10746/3.65423) then
             nsta = 57
           ELSEIF(br .le. 2.10855/3.65423) then
             nsta = 58
           ELSEIF(br .le. 2.11405/3.65423) then
             nsta = 60
           ELSEIF(br .le. 2.15805/3.65423) then
             nsta = 61
           ELSEIF(br .le. 2.19005/3.65423) then
             nsta = 62
           ELSEIF(br .le. 2.19213/3.65423) then
             nsta = 63
           ELSEIF(br .le. 2.19552/3.65423) then
             nsta = 64
           ELSEIF(br .le. 2.20272/3.65423) then
             nsta = 65
           ELSEIF(br .le. 2.20792/3.65423) then
             nsta = 66
           ELSEIF(br .le. 2.21632/3.65423) then
             nsta = 67
           ELSEIF(br .le. 2.25502/3.65423) then
             nsta = 68
           ELSEIF(br .le. 2.29002/3.65423) then
             nsta = 69
           ELSEIF(br .le. 2.29362/3.65423) then
             nsta = 70
           ELSEIF(br .le. 2.29802/3.65423) then
             nsta = 71
           ELSEIF(br .le. 2.33882/3.65423) then
             nsta = 72
           ELSEIF(br .le. 2.35182/3.65423) then
             nsta = 73
           ELSEIF(br .le. 2.36022/3.65423) then
             nsta = 74
           ELSEIF(br .le. 2.36472/3.65423) then
             nsta = 75
           ELSEIF(br .le. 2.37842/3.65423) then
             nsta = 76
           ELSEIF(br .le. 2.37984/3.65423) then
             nsta = 77
           ELSEIF(br .le. 2.38148/3.65423) then
             nsta = 78
           ELSEIF(br .le. 2.39008/3.65423) then
             nsta = 79
           ELSEIF(br .le. 2.40338/3.65423) then
             nsta = 80
           ELSEIF(br .le. 2.41148/3.65423) then
             nsta = 81
           ELSEIF(br .le. 2.41323/3.65423) then
             nsta = 82
           ELSEIF(br .le. 2.41563/3.65423) then
             nsta = 83
           ELSEIF(br .le. 2.42553/3.65423) then
             nsta = 84
           ELSEIF(br .le. 2.42826/3.65423) then
             nsta = 85
           ELSEIF(br .le. 2.43396/3.65423) then
             nsta = 86
           ELSEIF(br .le. 2.46096/3.65423) then
             nsta = 87
           ELSEIF(br .le. 2.47966/3.65423) then
             nsta = 88
           ELSEIF(br .le. 2.4813/3.65423) then
             nsta = 90
           ELSEIF(br .le. 2.4882/3.65423) then
             nsta = 91
           ELSEIF(br .le. 2.48951/3.65423) then
             nsta = 92
           ELSEIF(br .le. 2.50641/3.65423) then
             nsta = 93
           ELSEIF(br .le. 2.52171/3.65423) then
             nsta = 94
           ELSEIF(br .le. 2.53121/3.65423) then
             nsta = 95
           ELSEIF(br .le. 2.55421/3.65423) then
             nsta = 96
           ELSEIF(br .le. 2.55961/3.65423) then
             nsta = 97
           ELSEIF(br .le. 2.56234/3.65423) then
             nsta = 98
           ELSEIF(br .le. 2.60434/3.65423) then
             nsta = 99
           ELSEIF(br .le. 2.61604/3.65423) then
             nsta = 100
           ELSEIF(br .le. 2.65504/3.65423) then
             nsta = 101
           ELSEIF(br .le. 2.66084/3.65423) then
             nsta = 102
           ELSEIF(br .le. 2.66694/3.65423) then
             nsta = 103
           ELSEIF(br .le. 2.66814/3.65423) then
             nsta = 104
           ELSEIF(br .le. 2.68454/3.65423) then
             nsta = 105
           ELSEIF(br .le. 2.68585/3.65423) then
             nsta = 106
           ELSEIF(br .le. 2.70575/3.65423) then
             nsta = 107
           ELSEIF(br .le. 2.71315/3.65423) then
             nsta = 108
           ELSEIF(br .le. 2.71825/3.65423) then
             nsta = 109
           ELSEIF(br .le. 2.72335/3.65423) then
             nsta = 110
           ELSEIF(br .le. 2.73305/3.65423) then
             nsta = 111
           ELSEIF(br .le. 2.73795/3.65423) then
             nsta = 112
           ELSEIF(br .le. 2.73948/3.65423) then
             nsta = 113
           ELSEIF(br .le. 2.76258/3.65423) then
             nsta = 114
           ELSEIF(br .le. 2.76367/3.65423) then
             nsta = 115
           ELSEIF(br .le. 2.77807/3.65423) then
             nsta = 116
           ELSEIF(br .le. 2.78327/3.65423) then
             nsta = 117
           ELSEIF(br .le. 2.92327/3.65423) then
             nsta = 118
           ELSEIF(br .le. 2.92747/3.65423) then
             nsta = 119
           ELSEIF(br .le. 2.93537/3.65423) then
             nsta = 120
           ELSEIF(br .le. 2.94707/3.65423) then
             nsta = 121
           ELSEIF(br .le. 2.97007/3.65423) then
             nsta = 122
           ELSEIF(br .le. 2.97247/3.65423) then
             nsta = 123
           ELSEIF(br .le. 2.99017/3.65423) then
             nsta = 124
           ELSEIF(br .le. 2.99387/3.65423) then
             nsta = 125
           ELSEIF(br .le. 2.99507/3.65423) then
             nsta = 126
           ELSEIF(br .le. 3.01417/3.65423) then
             nsta = 127
           ELSEIF(br .le. 3.02037/3.65423) then
             nsta = 128
           ELSEIF(br .le. 3.02365/3.65423) then
             nsta = 129
           ELSEIF(br .le. 3.02474/3.65423) then
             nsta = 130
           ELSEIF(br .le. 3.07174/3.65423) then
             nsta = 131
           ELSEIF(br .le. 3.07684/3.65423) then
             nsta = 132
           ELSEIF(br .le. 3.08274/3.65423) then
             nsta = 133
           ELSEIF(br .le. 3.08405/3.65423) then
             nsta = 134
           ELSEIF(br .le. 3.08547/3.65423) then
             nsta = 135
           ELSEIF(br .le. 3.09617/3.65423) then
             nsta = 136
           ELSEIF(br .le. 3.13017/3.65423) then
             nsta = 137
           ELSEIF(br .le. 3.13279/3.65423) then
             nsta = 138
           ELSEIF(br .le. 3.15209/3.65423) then
             nsta = 139
           ELSEIF(br .le. 3.15849/3.65423) then
             nsta = 140
           ELSEIF(br .le. 3.16013/3.65423) then
             nsta = 141
           ELSEIF(br .le. 3.22013/3.65423) then
             nsta = 142
           ELSEIF(br .le. 3.23283/3.65423) then
             nsta = 143
           ELSEIF(br .le. 3.29383/3.65423) then
             nsta = 144
           ELSEIF(br .le. 3.30883/3.65423) then
             nsta = 145
           ELSEIF(br .le. 3.31673/3.65423) then
             nsta = 146
           ELSEIF(br .le. 3.32173/3.65423) then
             nsta = 147
           ELSEIF(br .le. 3.34173/3.65423) then
             nsta = 148
           ELSEIF(br .le. 3.36573/3.65423) then
             nsta = 149
           ELSEIF(br .le. 3.37173/3.65423) then
             nsta = 150
           ELSEIF(br .le. 3.38173/3.65423) then
             nsta = 151
           ELSEIF(br .le. 3.43573/3.65423) then
             nsta = 152
           ELSEIF(br .le. 3.45913/3.65423) then
             nsta = 153
           ELSEIF(br .le. 3.47303/3.65423) then
             nsta = 154
           ELSEIF(br .le. 3.48583/3.65423) then
             nsta = 155
           ELSEIF(br .le. 3.50943/3.65423) then
             nsta = 156
           ELSEIF(br .le. 3.5114/3.65423) then
             nsta = 157
           ELSEIF(br .le. 3.5159/3.65423) then
             nsta = 158
           ELSEIF(br .le. 3.5221/3.65423) then
             nsta = 159
           ELSEIF(br .le. 3.5261/3.65423) then
             nsta = 160
           ELSEIF(br .le. 3.5371/3.65423) then
             nsta = 161
           ELSEIF(br .le. 3.5498/3.65423) then
             nsta = 162
           ELSEIF(br .le. 3.5786/3.65423) then
             nsta = 163
           ELSEIF(br .le. 3.6046/3.65423) then
             nsta = 164
           ELSEIF(br .le. 3.6306/3.65423) then
             nsta = 165
           ELSEIF(br .le. 3.6375/3.65423) then
             nsta = 166
           ELSEIF(br .le. 3.6455/3.65423) then
             nsta = 167
           ELSEIF(br .le. 3.64703/3.65423) then
             nsta = 168
           ELSE
             nsta = 169
           ENDIF
         Case(60143)
               IF(br .le. 0.49/2.1218) then
             nsta = 1
           ELSEIF(br .le. 0.625/2.1218) then
             nsta = 3
           ELSEIF(br .le. 0.65/2.1218) then
             nsta = 9
           ELSEIF(br .le. 0.691/2.1218) then
             nsta = 10
           ELSEIF(br .le. 0.881/2.1218) then
             nsta = 13
           ELSEIF(br .le. 1.026/2.1218) then
             nsta = 14
           ELSEIF(br .le. 1.086/2.1218) then
             nsta = 16
           ELSEIF(br .le. 1.093/2.1218) then
             nsta = 17
           ELSEIF(br .le. 1.0966/2.1218) then
             nsta = 19
           ELSEIF(br .le. 1.1146/2.1218) then
             nsta = 23
           ELSEIF(br .le. 1.1356/2.1218) then
             nsta = 35
           ELSEIF(br .le. 1.1516/2.1218) then
             nsta = 38
           ELSEIF(br .le. 1.1595/2.1218) then
             nsta = 57
           ELSEIF(br .le. 1.1775/2.1218) then
             nsta = 61
           ELSEIF(br .le. 1.1829/2.1218) then
             nsta = 67
           ELSEIF(br .le. 1.2019/2.1218) then
             nsta = 71
           ELSEIF(br .le. 1.2129/2.1218) then
             nsta = 76
           ELSEIF(br .le. 1.2249/2.1218) then
             nsta = 77
           ELSEIF(br .le. 1.2989/2.1218) then
             nsta = 79
           ELSEIF(br .le. 1.3031/2.1218) then
             nsta = 85
           ELSEIF(br .le. 1.3067/2.1218) then
             nsta = 88
           ELSEIF(br .le. 1.3146/2.1218) then
             nsta = 89
           ELSEIF(br .le. 1.4286/2.1218) then
             nsta = 90
           ELSEIF(br .le. 1.4526/2.1218) then
             nsta = 92
           ELSEIF(br .le. 1.5096/2.1218) then
             nsta = 95
           ELSEIF(br .le. 1.5286/2.1218) then
             nsta = 99
           ELSEIF(br .le. 1.5371/2.1218) then
             nsta = 106
           ELSEIF(br .le. 1.5432/2.1218) then
             nsta = 108
           ELSEIF(br .le. 1.5682/2.1218) then
             nsta = 109
           ELSEIF(br .le. 1.6312/2.1218) then
             nsta = 116
           ELSEIF(br .le. 1.6442/2.1218) then
             nsta = 118
           ELSEIF(br .le. 1.6503/2.1218) then
             nsta = 120
           ELSEIF(br .le. 1.6576/2.1218) then
             nsta = 121
           ELSEIF(br .le. 1.6618/2.1218) then
             nsta = 123
           ELSEIF(br .le. 1.6654/2.1218) then
             nsta = 124
           ELSEIF(br .le. 1.6739/2.1218) then
             nsta = 126
           ELSEIF(br .le. 1.6939/2.1218) then
             nsta = 129
           ELSEIF(br .le. 1.7069/2.1218) then
             nsta = 131
           ELSEIF(br .le. 1.7189/2.1218) then
             nsta = 134
           ELSEIF(br .le. 1.7399/2.1218) then
             nsta = 138
           ELSEIF(br .le. 1.7539/2.1218) then
             nsta = 144
           ELSEIF(br .le. 1.8729/2.1218) then
             nsta = 146
           ELSEIF(br .le. 1.9169/2.1218) then
             nsta = 150
           ELSEIF(br .le. 1.9299/2.1218) then
             nsta = 155
           ELSEIF(br .le. 1.9599/2.1218) then
             nsta = 156
           ELSEIF(br .le. 1.9678/2.1218) then
             nsta = 157
           ELSEIF(br .le. 1.9918/2.1218) then
             nsta = 165
           ELSEIF(br .le. 2.0058/2.1218) then
             nsta = 171
           ELSEIF(br .le. 2.0608/2.1218) then
             nsta = 172
           ELSE
             nsta = 194
           ENDIF
         Case(60144)
               IF(br .le. 3.18/4.68) then
             nsta = 2
           ELSE
             nsta = 4
           ENDIF
         Case(60145)
               IF(br .le. 0.32/0.48) then
             nsta = 6
           ELSE
             nsta = 41
           ENDIF
         Case(60146)
               IF(br .le. 0.368/1.5561) then
             nsta = 1
           ELSEIF(br .le. 0.563/1.5561) then
             nsta = 3
           ELSEIF(br .le. 0.5791/1.5561) then
             nsta = 4
           ELSEIF(br .le. 0.8751/1.5561) then
             nsta = 7
           ELSEIF(br .le. 0.8799/1.5561) then
             nsta = 13
           ELSEIF(br .le. 0.9749/1.5561) then
             nsta = 15
           ELSEIF(br .le. 0.9858/1.5561) then
             nsta = 17
           ELSEIF(br .le. 0.9898/1.5561) then
             nsta = 18
           ELSEIF(br .le. 1.0328/1.5561) then
             nsta = 22
           ELSEIF(br .le. 1.0458/1.5561) then
             nsta = 23
           ELSEIF(br .le. 1.0688/1.5561) then
             nsta = 24
           ELSEIF(br .le. 1.0722/1.5561) then
             nsta = 29
           ELSEIF(br .le. 1.0818/1.5561) then
             nsta = 32
           ELSEIF(br .le. 1.1418/1.5561) then
             nsta = 37
           ELSEIF(br .le. 1.1527/1.5561) then
             nsta = 39
           ELSEIF(br .le. 1.1578/1.5561) then
             nsta = 41
           ELSEIF(br .le. 1.2028/1.5561) then
             nsta = 53
           ELSEIF(br .le. 1.2057/1.5561) then
             nsta = 56
           ELSEIF(br .le. 1.2128/1.5561) then
             nsta = 57
           ELSEIF(br .le. 1.2225/1.5561) then
             nsta = 60
           ELSEIF(br .le. 1.2334/1.5561) then
             nsta = 64
           ELSEIF(br .le. 1.2567/1.5561) then
             nsta = 65
           ELSEIF(br .le. 1.2947/1.5561) then
             nsta = 66
           ELSEIF(br .le. 1.301/1.5561) then
             nsta = 67
           ELSEIF(br .le. 1.373/1.5561) then
             nsta = 70
           ELSEIF(br .le. 1.3863/1.5561) then
             nsta = 72
           ELSEIF(br .le. 1.3914/1.5561) then
             nsta = 74
           ELSEIF(br .le. 1.4193/1.5561) then
             nsta = 76
           ELSEIF(br .le. 1.4242/1.5561) then
             nsta = 84
           ELSEIF(br .le. 1.4283/1.5561) then
             nsta = 91
           ELSEIF(br .le. 1.4534/1.5561) then
             nsta = 93
           ELSEIF(br .le. 1.4631/1.5561) then
             nsta = 96
           ELSEIF(br .le. 1.4671/1.5561) then
             nsta = 101
           ELSEIF(br .le. 1.4692/1.5561) then
             nsta = 104
           ELSEIF(br .le. 1.474/1.5561) then
             nsta = 105
           ELSEIF(br .le. 1.4955/1.5561) then
             nsta = 110
           ELSEIF(br .le. 1.5016/1.5561) then
             nsta = 113
           ELSEIF(br .le. 1.5456/1.5561) then
             nsta = 114
           ELSEIF(br .le. 1.5507/1.5561) then
             nsta = 117
           ELSE
             nsta = 122
           ENDIF
         Case(60147)
               IF(br .le. 0.00108/0.030524) then
             nsta = 4
           ELSEIF(br .le. 0.00416/0.030524) then
             nsta = 5
           ELSEIF(br .le. 0.00557/0.030524) then
             nsta = 6
           ELSEIF(br .le. 0.00599/0.030524) then
             nsta = 7
           ELSEIF(br .le. 0.01349/0.030524) then
             nsta = 10
           ELSEIF(br .le. 0.02139/0.030524) then
             nsta = 11
           ELSEIF(br .le. 0.021555/0.030524) then
             nsta = 14
           ELSEIF(br .le. 0.021884/0.030524) then
             nsta = 15
           ELSEIF(br .le. 0.022544/0.030524) then
             nsta = 42
           ELSEIF(br .le. 0.023394/0.030524) then
             nsta = 43
           ELSEIF(br .le. 0.024974/0.030524) then
             nsta = 48
           ELSEIF(br .le. 0.025704/0.030524) then
             nsta = 49
           ELSEIF(br .le. 0.026154/0.030524) then
             nsta = 50
           ELSEIF(br .le. 0.028974/0.030524) then
             nsta = 66
           ELSE
             nsta = 73
           ENDIF
         Case(60149)
               IF(br .le. 0.00085/0.02104) then
             nsta = 3
           ELSEIF(br .le. 0.00111/0.02104) then
             nsta = 6
           ELSEIF(br .le. 0.00371/0.02104) then
             nsta = 8
           ELSEIF(br .le. 0.01671/0.02104) then
             nsta = 14
           ELSEIF(br .le. 0.01722/0.02104) then
             nsta = 15
           ELSEIF(br .le. 0.01756/0.02104) then
             nsta = 19
           ELSEIF(br .le. 0.01794/0.02104) then
             nsta = 23
           ELSEIF(br .le. 0.01944/0.02104) then
             nsta = 31
           ELSEIF(br .le. 0.01971/0.02104) then
             nsta = 32
           ELSEIF(br .le. 0.02043/0.02104) then
             nsta = 35
           ELSE
             nsta = 39
           ENDIF
         Case(60151)
               IF(br .le. 0.00054/0.0334) then
             nsta = 0
           ELSEIF(br .le. 0.00254/0.0334) then
             nsta = 2
           ELSEIF(br .le. 0.0029/0.0334) then
             nsta = 7
           ELSEIF(br .le. 0.0043/0.0334) then
             nsta = 13
           ELSEIF(br .le. 0.0053/0.0334) then
             nsta = 14
           ELSEIF(br .le. 0.00577/0.0334) then
             nsta = 18
           ELSEIF(br .le. 0.00599/0.0334) then
             nsta = 22
           ELSEIF(br .le. 0.01099/0.0334) then
             nsta = 29
           ELSEIF(br .le. 0.01189/0.0334) then
             nsta = 31
           ELSEIF(br .le. 0.01279/0.0334) then
             nsta = 32
           ELSEIF(br .le. 0.01311/0.0334) then
             nsta = 34
           ELSEIF(br .le. 0.01451/0.0334) then
             nsta = 35
           ELSEIF(br .le. 0.01601/0.0334) then
             nsta = 36
           ELSEIF(br .le. 0.01635/0.0334) then
             nsta = 43
           ELSEIF(br .le. 0.01735/0.0334) then
             nsta = 44
           ELSEIF(br .le. 0.01795/0.0334) then
             nsta = 45
           ELSEIF(br .le. 0.01885/0.0334) then
             nsta = 46
           ELSEIF(br .le. 0.01955/0.0334) then
             nsta = 47
           ELSEIF(br .le. 0.02015/0.0334) then
             nsta = 51
           ELSEIF(br .le. 0.02075/0.0334) then
             nsta = 52
           ELSEIF(br .le. 0.02115/0.0334) then
             nsta = 59
           ELSEIF(br .le. 0.02145/0.0334) then
             nsta = 61
           ELSEIF(br .le. 0.02305/0.0334) then
             nsta = 64
           ELSEIF(br .le. 0.0235/0.0334) then
             nsta = 66
           ELSEIF(br .le. 0.0249/0.0334) then
             nsta = 68
           ELSEIF(br .le. 0.0319/0.0334) then
             nsta = 73
           ELSE
             nsta = 78
           ENDIF
         Case(62148)
               IF(br .le. 5.2/50.7) then
             nsta = 10
           ELSEIF(br .le. 11./50.7) then
             nsta = 54
           ELSEIF(br .le. 20./50.7) then
             nsta = 96
           ELSEIF(br .le. 25./50.7) then
             nsta = 116
           ELSEIF(br .le. 30.2/50.7) then
             nsta = 117
           ELSEIF(br .le. 32.7/50.7) then
             nsta = 125
           ELSEIF(br .le. 41.4/50.7) then
             nsta = 158
           ELSEIF(br .le. 45.6/50.7) then
             nsta = 162
           ELSE
             nsta = 163
           ENDIF
         Case(62150)
               IF(br .le. 40.2/292.) then
             nsta = 3
           ELSEIF(br .le. 43.4/292.) then
             nsta = 13
           ELSEIF(br .le. 46./292.) then
             nsta = 15
           ELSEIF(br .le. 46.8/292.) then
             nsta = 18
           ELSEIF(br .le. 47.6/292.) then
             nsta = 27
           ELSEIF(br .le. 48.8/292.) then
             nsta = 32
           ELSEIF(br .le. 53.6/292.) then
             nsta = 34
           ELSEIF(br .le. 64.4/292.) then
             nsta = 38
           ELSEIF(br .le. 71.7/292.) then
             nsta = 42
           ELSEIF(br .le. 75.2/292.) then
             nsta = 44
           ELSEIF(br .le. 79./292.) then
             nsta = 50
           ELSEIF(br .le. 81.4/292.) then
             nsta = 54
           ELSEIF(br .le. 83./292.) then
             nsta = 59
           ELSEIF(br .le. 91.2/292.) then
             nsta = 61
           ELSEIF(br .le. 92.9/292.) then
             nsta = 66
           ELSEIF(br .le. 97.7/292.) then
             nsta = 71
           ELSEIF(br .le. 128.1/292.) then
             nsta = 75
           ELSEIF(br .le. 136.1/292.) then
             nsta = 81
           ELSEIF(br .le. 144.7/292.) then
             nsta = 88
           ELSEIF(br .le. 148./292.) then
             nsta = 89
           ELSEIF(br .le. 153.2/292.) then
             nsta = 94
           ELSEIF(br .le. 156.4/292.) then
             nsta = 97
           ELSEIF(br .le. 165.4/292.) then
             nsta = 99
           ELSEIF(br .le. 167.8/292.) then
             nsta = 102
           ELSEIF(br .le. 174.2/292.) then
             nsta = 112
           ELSEIF(br .le. 179.6/292.) then
             nsta = 115
           ELSEIF(br .le. 188.6/292.) then
             nsta = 121
           ELSEIF(br .le. 204.5/292.) then
             nsta = 122
           ELSEIF(br .le. 210.3/292.) then
             nsta = 130
           ELSEIF(br .le. 216.2/292.) then
             nsta = 133
           ELSEIF(br .le. 221.7/292.) then
             nsta = 137
           ELSEIF(br .le. 231.7/292.) then
             nsta = 138
           ELSEIF(br .le. 238.3/292.) then
             nsta = 140
           ELSEIF(br .le. 242.8/292.) then
             nsta = 141
           ELSEIF(br .le. 249.5/292.) then
             nsta = 142
           ELSEIF(br .le. 253.5/292.) then
             nsta = 145
           ELSEIF(br .le. 256.8/292.) then
             nsta = 147
           ELSEIF(br .le. 261.2/292.) then
             nsta = 148
           ELSEIF(br .le. 267.3/292.) then
             nsta = 149
           ELSEIF(br .le. 269.7/292.) then
             nsta = 150
           ELSEIF(br .le. 274./292.) then
             nsta = 152
           ELSEIF(br .le. 277.2/292.) then
             nsta = 153
           ELSEIF(br .le. 279.6/292.) then
             nsta = 154
           ELSEIF(br .le. 283.8/292.) then
             nsta = 157
           ELSEIF(br .le. 287./292.) then
             nsta = 159
           ELSE
             nsta = 160
           ENDIF
         Case(62151)
               IF(br .le. 0.26/3.397) then
             nsta = 1
           ELSEIF(br .le. 0.484/3.397) then
             nsta = 5
           ELSEIF(br .le. 0.905/3.397) then
             nsta = 13
           ELSEIF(br .le. 0.944/3.397) then
             nsta = 16
           ELSEIF(br .le. 0.997/3.397) then
             nsta = 17
           ELSEIF(br .le. 1.487/3.397) then
             nsta = 18
           ELSEIF(br .le. 1.521/3.397) then
             nsta = 20
           ELSEIF(br .le. 1.551/3.397) then
             nsta = 21
           ELSEIF(br .le. 1.561/3.397) then
             nsta = 25
           ELSEIF(br .le. 1.758/3.397) then
             nsta = 32
           ELSEIF(br .le. 1.784/3.397) then
             nsta = 37
           ELSEIF(br .le. 1.811/3.397) then
             nsta = 39
           ELSEIF(br .le. 1.832/3.397) then
             nsta = 45
           ELSEIF(br .le. 1.939/3.397) then
             nsta = 48
           ELSEIF(br .le. 2.015/3.397) then
             nsta = 50
           ELSEIF(br .le. 2.071/3.397) then
             nsta = 53
           ELSEIF(br .le. 2.127/3.397) then
             nsta = 55
           ELSEIF(br .le. 2.205/3.397) then
             nsta = 60
           ELSEIF(br .le. 2.264/3.397) then
             nsta = 68
           ELSEIF(br .le. 2.397/3.397) then
             nsta = 70
           ELSEIF(br .le. 2.422/3.397) then
             nsta = 76
           ELSEIF(br .le. 2.592/3.397) then
             nsta = 83
           ELSEIF(br .le. 2.692/3.397) then
             nsta = 93
           ELSEIF(br .le. 2.726/3.397) then
             nsta = 99
           ELSEIF(br .le. 2.761/3.397) then
             nsta = 110
           ELSEIF(br .le. 2.802/3.397) then
             nsta = 114
           ELSEIF(br .le. 2.842/3.397) then
             nsta = 120
           ELSEIF(br .le. 2.982/3.397) then
             nsta = 121
           ELSEIF(br .le. 3.041/3.397) then
             nsta = 124
           ELSEIF(br .le. 3.08/3.397) then
             nsta = 129
           ELSEIF(br .le. 3.34/3.397) then
             nsta = 131
           ELSE
             nsta = 132
           ENDIF
         Case(62153)
               IF(br .le. 0.092/9.6945) then
             nsta = 0
           ELSEIF(br .le. 0.098/9.6945) then
             nsta = 1
           ELSEIF(br .le. 0.708/9.6945) then
             nsta = 2
           ELSEIF(br .le. 3.408/9.6945) then
             nsta = 9
           ELSEIF(br .le. 3.46/9.6945) then
             nsta = 20
           ELSEIF(br .le. 3.525/9.6945) then
             nsta = 21
           ELSEIF(br .le. 3.5293/9.6945) then
             nsta = 22
           ELSEIF(br .le. 3.5355/9.6945) then
             nsta = 23
           ELSEIF(br .le. 3.7255/9.6945) then
             nsta = 27
           ELSEIF(br .le. 3.8655/9.6945) then
             nsta = 29
           ELSEIF(br .le. 4.1155/9.6945) then
             nsta = 36
           ELSEIF(br .le. 4.1365/9.6945) then
             nsta = 41
           ELSEIF(br .le. 4.7665/9.6945) then
             nsta = 44
           ELSEIF(br .le. 4.7765/9.6945) then
             nsta = 45
           ELSEIF(br .le. 5.1365/9.6945) then
             nsta = 48
           ELSEIF(br .le. 5.1615/9.6945) then
             nsta = 53
           ELSEIF(br .le. 5.4015/9.6945) then
             nsta = 55
           ELSEIF(br .le. 5.4765/9.6945) then
             nsta = 67
           ELSEIF(br .le. 5.6965/9.6945) then
             nsta = 71
           ELSEIF(br .le. 5.9265/9.6945) then
             nsta = 73
           ELSEIF(br .le. 5.9825/9.6945) then
             nsta = 76
           ELSEIF(br .le. 6.0385/9.6945) then
             nsta = 83
           ELSEIF(br .le. 6.6185/9.6945) then
             nsta = 92
           ELSEIF(br .le. 6.7055/9.6945) then
             nsta = 98
           ELSEIF(br .le. 7.2655/9.6945) then
             nsta = 111
           ELSEIF(br .le. 7.4355/9.6945) then
             nsta = 112
           ELSEIF(br .le. 7.8655/9.6945) then
             nsta = 116
           ELSEIF(br .le. 7.9335/9.6945) then
             nsta = 119
           ELSEIF(br .le. 8.0145/9.6945) then
             nsta = 120
           ELSEIF(br .le. 8.5545/9.6945) then
             nsta = 136
           ELSEIF(br .le. 8.7245/9.6945) then
             nsta = 143
           ELSEIF(br .le. 8.9345/9.6945) then
             nsta = 192
           ELSEIF(br .le. 9.2245/9.6945) then
             nsta = 193
           ELSEIF(br .le. 9.5745/9.6945) then
             nsta = 238
           ELSE
             nsta = 246
           ENDIF
         Case(63152)
               IF(br .le. 1.81/116.08) then
             nsta = 0
           ELSEIF(br .le. 5.91/116.08) then
             nsta = 3
           ELSEIF(br .le. 6.71/116.08) then
             nsta = 6
           ELSEIF(br .le. 7.43/116.08) then
             nsta = 14
           ELSEIF(br .le. 10.11/116.08) then
             nsta = 19
           ELSEIF(br .le. 10.25/116.08) then
             nsta = 21
           ELSEIF(br .le. 11.35/116.08) then
             nsta = 22
           ELSEIF(br .le. 11.95/116.08) then
             nsta = 31
           ELSEIF(br .le. 13.75/116.08) then
             nsta = 38
           ELSEIF(br .le. 21.95/116.08) then
             nsta = 44
           ELSEIF(br .le. 23.05/116.08) then
             nsta = 45
           ELSEIF(br .le. 25.35/116.08) then
             nsta = 56
           ELSEIF(br .le. 27.74/116.08) then
             nsta = 58
           ELSEIF(br .le. 28.84/116.08) then
             nsta = 61
           ELSEIF(br .le. 29.34/116.08) then
             nsta = 69
           ELSEIF(br .le. 30.24/116.08) then
             nsta = 76
           ELSEIF(br .le. 30.84/116.08) then
             nsta = 82
           ELSEIF(br .le. 33.14/116.08) then
             nsta = 94
           ELSEIF(br .le. 36.64/116.08) then
             nsta = 96
           ELSEIF(br .le. 38.14/116.08) then
             nsta = 102
           ELSEIF(br .le. 39.04/116.08) then
             nsta = 110
           ELSEIF(br .le. 40.44/116.08) then
             nsta = 118
           ELSEIF(br .le. 41.34/116.08) then
             nsta = 123
           ELSEIF(br .le. 45.04/116.08) then
             nsta = 124
           ELSEIF(br .le. 46.24/116.08) then
             nsta = 127
           ELSEIF(br .le. 48.14/116.08) then
             nsta = 129
           ELSEIF(br .le. 49.54/116.08) then
             nsta = 133
           ELSEIF(br .le. 50.24/116.08) then
             nsta = 134
           ELSEIF(br .le. 51.77/116.08) then
             nsta = 137
           ELSEIF(br .le. 52.67/116.08) then
             nsta = 138
           ELSEIF(br .le. 54.14/116.08) then
             nsta = 139
           ELSEIF(br .le. 56.83/116.08) then
             nsta = 140
           ELSEIF(br .le. 57.49/116.08) then
             nsta = 143
           ELSEIF(br .le. 58.71/116.08) then
             nsta = 144
           ELSEIF(br .le. 59.73/116.08) then
             nsta = 145
           ELSEIF(br .le. 61.13/116.08) then
             nsta = 149
           ELSEIF(br .le. 62.57/116.08) then
             nsta = 150
           ELSEIF(br .le. 65.27/116.08) then
             nsta = 152
           ELSEIF(br .le. 67.67/116.08) then
             nsta = 154
           ELSEIF(br .le. 72.97/116.08) then
             nsta = 155
           ELSEIF(br .le. 74.27/116.08) then
             nsta = 159
           ELSEIF(br .le. 75.37/116.08) then
             nsta = 160
           ELSEIF(br .le. 78.27/116.08) then
             nsta = 162
           ELSEIF(br .le. 79.67/116.08) then
             nsta = 164
           ELSEIF(br .le. 81.57/116.08) then
             nsta = 165
           ELSEIF(br .le. 83.07/116.08) then
             nsta = 167
           ELSEIF(br .le. 90.07/116.08) then
             nsta = 168
           ELSEIF(br .le. 90.66/116.08) then
             nsta = 170
           ELSEIF(br .le. 94.26/116.08) then
             nsta = 171
           ELSEIF(br .le. 95.56/116.08) then
             nsta = 172
           ELSEIF(br .le. 98.46/116.08) then
             nsta = 173
           ELSEIF(br .le. 101.36/116.08) then
             nsta = 174
           ELSEIF(br .le. 102.89/116.08) then
             nsta = 175
           ELSEIF(br .le. 104.38/116.08) then
             nsta = 176
           ELSEIF(br .le. 106.08/116.08) then
             nsta = 177
           ELSEIF(br .le. 115.28/116.08) then
             nsta = 180
           ELSE
             nsta = 181
           ENDIF
         Case(63154)
               IF(br .le. 1./9.23) then
             nsta = 110
           ELSEIF(br .le. 2.2/9.23) then
             nsta = 118
           ELSEIF(br .le. 2.93/9.23) then
             nsta = 141
           ELSEIF(br .le. 4.23/9.23) then
             nsta = 147
           ELSEIF(br .le. 5.33/9.23) then
             nsta = 149
           ELSEIF(br .le. 7.13/9.23) then
             nsta = 150
           ELSEIF(br .le. 8.03/9.23) then
             nsta = 150
           ELSE
             nsta = 151
           ENDIF
         Case(64153)
               IF(br .le. 0.01423/0.09104) then
             nsta = 0
           ELSEIF(br .le. 0.01452/0.09104) then
             nsta = 4
           ELSEIF(br .le. 0.01524/0.09104) then
             nsta = 10
           ELSEIF(br .le. 0.01583/0.09104) then
             nsta = 15
           ELSEIF(br .le. 0.02164/0.09104) then
             nsta = 16
           ELSEIF(br .le. 0.02458/0.09104) then
             nsta = 19
           ELSEIF(br .le. 0.02837/0.09104) then
             nsta = 22
           ELSEIF(br .le. 0.02896/0.09104) then
             nsta = 28
           ELSEIF(br .le. 0.03132/0.09104) then
             nsta = 31
           ELSEIF(br .le. 0.03299/0.09104) then
             nsta = 34
           ELSEIF(br .le. 0.03722/0.09104) then
             nsta = 37
           ELSEIF(br .le. 0.03781/0.09104) then
             nsta = 39
           ELSEIF(br .le. 0.03807/0.09104) then
             nsta = 41
           ELSEIF(br .le. 0.03846/0.09104) then
             nsta = 44
           ELSEIF(br .le. 0.03877/0.09104) then
             nsta = 45
           ELSEIF(br .le. 0.03947/0.09104) then
             nsta = 52
           ELSEIF(br .le. 0.03968/0.09104) then
             nsta = 55
           ELSEIF(br .le. 0.0403/0.09104) then
             nsta = 56
           ELSEIF(br .le. 0.04136/0.09104) then
             nsta = 57
           ELSEIF(br .le. 0.04172/0.09104) then
             nsta = 62
           ELSEIF(br .le. 0.04219/0.09104) then
             nsta = 67
           ELSEIF(br .le. 0.04247/0.09104) then
             nsta = 68
           ELSEIF(br .le. 0.04293/0.09104) then
             nsta = 70
           ELSEIF(br .le. 0.04384/0.09104) then
             nsta = 71
           ELSEIF(br .le. 0.04599/0.09104) then
             nsta = 74
           ELSEIF(br .le. 0.04682/0.09104) then
             nsta = 78
           ELSEIF(br .le. 0.04723/0.09104) then
             nsta = 80
           ELSEIF(br .le. 0.049/0.09104) then
             nsta = 85
           ELSEIF(br .le. 0.04937/0.09104) then
             nsta = 87
           ELSEIF(br .le. 0.0497/0.09104) then
             nsta = 89
           ELSEIF(br .le. 0.05025/0.09104) then
             nsta = 92
           ELSEIF(br .le. 0.05152/0.09104) then
             nsta = 97
           ELSEIF(br .le. 0.06305/0.09104) then
             nsta = 100
           ELSEIF(br .le. 0.06362/0.09104) then
             nsta = 104
           ELSEIF(br .le. 0.06409/0.09104) then
             nsta = 106
           ELSEIF(br .le. 0.06542/0.09104) then
             nsta = 107
           ELSEIF(br .le. 0.06623/0.09104) then
             nsta = 112
           ELSEIF(br .le. 0.06691/0.09104) then
             nsta = 113
           ELSEIF(br .le. 0.06798/0.09104) then
             nsta = 114
           ELSEIF(br .le. 0.06834/0.09104) then
             nsta = 121
           ELSEIF(br .le. 0.06875/0.09104) then
             nsta = 123
           ELSEIF(br .le. 0.06953/0.09104) then
             nsta = 124
           ELSEIF(br .le. 0.07003/0.09104) then
             nsta = 125
           ELSEIF(br .le. 0.07024/0.09104) then
             nsta = 126
           ELSEIF(br .le. 0.07626/0.09104) then
             nsta = 128
           ELSEIF(br .le. 0.07725/0.09104) then
             nsta = 130
           ELSEIF(br .le. 0.08125/0.09104) then
             nsta = 134
           ELSEIF(br .le. 0.08208/0.09104) then
             nsta = 136
           ELSEIF(br .le. 0.08315/0.09104) then
             nsta = 138
           ELSEIF(br .le. 0.08426/0.09104) then
             nsta = 139
           ELSEIF(br .le. 0.08462/0.09104) then
             nsta = 145
           ELSEIF(br .le. 0.08613/0.09104) then
             nsta = 150
           ELSEIF(br .le. 0.08683/0.09104) then
             nsta = 152
           ELSE
             nsta = 171
           ENDIF
         Case(64155)
               IF(br .le. 0.0041/0.37363) then
             nsta = 0
           ELSEIF(br .le. 0.0174/0.37363) then
             nsta = 3
           ELSEIF(br .le. 0.01842/0.37363) then
             nsta = 11
           ELSEIF(br .le. 0.02452/0.37363) then
             nsta = 14
           ELSEIF(br .le. 0.02762/0.37363) then
             nsta = 18
           ELSEIF(br .le. 0.06262/0.37363) then
             nsta = 24
           ELSEIF(br .le. 0.06323/0.37363) then
             nsta = 31
           ELSEIF(br .le. 0.07763/0.37363) then
             nsta = 34
           ELSEIF(br .le. 0.11463/0.37363) then
             nsta = 36
           ELSEIF(br .le. 0.21063/0.37363) then
             nsta = 38
           ELSEIF(br .le. 0.21913/0.37363) then
             nsta = 63
           ELSEIF(br .le. 0.25813/0.37363) then
             nsta = 64
           ELSEIF(br .le. 0.26043/0.37363) then
             nsta = 68
           ELSEIF(br .le. 0.27133/0.37363) then
             nsta = 78
           ELSEIF(br .le. 0.27483/0.37363) then
             nsta = 82
           ELSEIF(br .le. 0.28403/0.37363) then
             nsta = 85
           ELSEIF(br .le. 0.28593/0.37363) then
             nsta = 89
           ELSEIF(br .le. 0.28883/0.37363) then
             nsta = 91
           ELSEIF(br .le. 0.29173/0.37363) then
             nsta = 95
           ELSEIF(br .le. 0.29673/0.37363) then
             nsta = 96
           ELSEIF(br .le. 0.30043/0.37363) then
             nsta = 103
           ELSEIF(br .le. 0.30433/0.37363) then
             nsta = 106
           ELSEIF(br .le. 0.30623/0.37363) then
             nsta = 110
           ELSEIF(br .le. 0.30913/0.37363) then
             nsta = 112
           ELSEIF(br .le. 0.31143/0.37363) then
             nsta = 115
           ELSEIF(br .le. 0.31953/0.37363) then
             nsta = 120
           ELSEIF(br .le. 0.32393/0.37363) then
             nsta = 122
           ELSEIF(br .le. 0.32683/0.37363) then
             nsta = 124
           ELSEIF(br .le. 0.36583/0.37363) then
             nsta = 127
           ELSE
             nsta = 129
           ENDIF
         Case(64156)
               IF(br .le. 24/1283.4) then
             nsta = 34
           ELSEIF(br .le. 31.4/1283.4) then
             nsta = 47
           ELSEIF(br .le. 35.7/1283.4) then
             nsta = 56
           ELSEIF(br .le. 37.5/1283.4) then
             nsta = 61
           ELSEIF(br .le. 41.6/1283.4) then
             nsta = 62
           ELSEIF(br .le. 43.3/1283.4) then
             nsta = 70
           ELSEIF(br .le. 92.3/1283.4) then
             nsta = 75
           ELSEIF(br .le. 110.3/1283.4) then
             nsta = 78
           ELSEIF(br .le. 168.3/1283.4) then
             nsta = 82
           ELSEIF(br .le. 220.3/1283.4) then
             nsta = 97
           ELSEIF(br .le. 239.3/1283.4) then
             nsta = 101
           ELSEIF(br .le. 244.2/1283.4) then
             nsta = 104
           ELSEIF(br .le. 246.3/1283.4) then
             nsta = 106
           ELSEIF(br .le. 253./1283.4) then
             nsta = 109
           ELSEIF(br .le. 255.3/1283.4) then
             nsta = 112
           ELSEIF(br .le. 269.3/1283.4) then
             nsta = 116
           ELSEIF(br .le. 270.3/1283.4) then
             nsta = 118
           ELSEIF(br .le. 272./1283.4) then
             nsta = 121
           ELSEIF(br .le. 274.1/1283.4) then
             nsta = 123
           ELSEIF(br .le. 275.6/1283.4) then
             nsta = 126
           ELSEIF(br .le. 279.8/1283.4) then
             nsta = 129
           ELSEIF(br .le. 285.3/1283.4) then
             nsta = 132
           ELSEIF(br .le. 299.3/1283.4) then
             nsta = 133
           ELSEIF(br .le. 308.3/1283.4) then
             nsta = 136
           ELSEIF(br .le. 321.3/1283.4) then
             nsta = 139
           ELSEIF(br .le. 334.3/1283.4) then
             nsta = 141
           ELSEIF(br .le. 347.3/1283.4) then
             nsta = 143
           ELSEIF(br .le. 350.8/1283.4) then
             nsta = 144
           ELSEIF(br .le. 353.7/1283.4) then
             nsta = 147
           ELSEIF(br .le. 394.7/1283.4) then
             nsta = 151
           ELSEIF(br .le. 397.3/1283.4) then
             nsta = 154
           ELSEIF(br .le. 399.5/1283.4) then
             nsta = 166
           ELSEIF(br .le. 403.7/1283.4) then
             nsta = 169
           ELSEIF(br .le. 431.7/1283.4) then
             nsta = 176
           ELSEIF(br .le. 438.3/1283.4) then
             nsta = 183
           ELSEIF(br .le. 443.3/1283.4) then
             nsta = 184
           ELSEIF(br .le. 466.3/1283.4) then
             nsta = 191
           ELSEIF(br .le. 478.3/1283.4) then
             nsta = 193
           ELSEIF(br .le. 482.8/1283.4) then
             nsta = 196
           ELSEIF(br .le. 486.4/1283.4) then
             nsta = 200
           ELSEIF(br .le. 502.4/1283.4) then
             nsta = 204
           ELSEIF(br .le. 529.4/1283.4) then
             nsta = 205
           ELSEIF(br .le. 536.6/1283.4) then
             nsta = 207
           ELSEIF(br .le. 565.6/1283.4) then
             nsta = 208
           ELSEIF(br .le. 574.6/1283.4) then
             nsta = 210
           ELSEIF(br .le. 582.6/1283.4) then
             nsta = 213
           ELSEIF(br .le. 588.6/1283.4) then
             nsta = 216
           ELSEIF(br .le. 612.6/1283.4) then
             nsta = 219
           ELSEIF(br .le. 622.6/1283.4) then
             nsta = 220
           ELSEIF(br .le. 634.6/1283.4) then
             nsta = 221
           ELSEIF(br .le. 642.1/1283.4) then
             nsta = 222
           ELSEIF(br .le. 648.6/1283.4) then
             nsta = 223
           ELSEIF(br .le. 656.3/1283.4) then
             nsta = 223
           ELSEIF(br .le. 670.3/1283.4) then
             nsta = 224
           ELSEIF(br .le. 683.3/1283.4) then
             nsta = 230
           ELSEIF(br .le. 688./1283.4) then
             nsta = 231
           ELSEIF(br .le. 691.6/1283.4) then
             nsta = 233
           ELSEIF(br .le. 697.6/1283.4) then
             nsta = 236
           ELSEIF(br .le. 712.6/1283.4) then
             nsta = 237
           ELSEIF(br .le. 734.6/1283.4) then
             nsta = 237
           ELSEIF(br .le. 743.6/1283.4) then
             nsta = 238
           ELSEIF(br .le. 754.6/1283.4) then
             nsta = 239
           ELSEIF(br .le. 765.6/1283.4) then
             nsta = 239
           ELSEIF(br .le. 773.6/1283.4) then
             nsta = 240
           ELSEIF(br .le. 780.1/1283.4) then
             nsta = 240
           ELSEIF(br .le. 785.3/1283.4) then
             nsta = 240
           ELSEIF(br .le. 791.7/1283.4) then
             nsta = 241
           ELSEIF(br .le. 806.7/1283.4) then
             nsta = 242
           ELSEIF(br .le. 827.7/1283.4) then
             nsta = 242
           ELSEIF(br .le. 842.7/1283.4) then
             nsta = 243
           ELSEIF(br .le. 856.7/1283.4) then
             nsta = 244
           ELSEIF(br .le. 874.7/1283.4) then
             nsta = 245
           ELSEIF(br .le. 880.6/1283.4) then
             nsta = 246
           ELSEIF(br .le. 893.6/1283.4) then
             nsta = 246
           ELSEIF(br .le. 899.4/1283.4) then
             nsta = 247
           ELSEIF(br .le. 906.4/1283.4) then
             nsta = 248
           ELSEIF(br .le. 919.4/1283.4) then
             nsta = 248
           ELSEIF(br .le. 941.4/1283.4) then
             nsta = 249
           ELSEIF(br .le. 959.4/1283.4) then
             nsta = 250
           ELSEIF(br .le. 976.4/1283.4) then
             nsta = 251
           ELSEIF(br .le. 985.4/1283.4) then
             nsta = 251
           ELSEIF(br .le. 998.4/1283.4) then
             nsta = 251
           ELSEIF(br .le. 1012.4/1283.4) then
             nsta = 251
           ELSEIF(br .le. 1022.4/1283.4) then
             nsta = 253
           ELSEIF(br .le. 1030.4/1283.4) then
             nsta = 253
           ELSEIF(br .le. 1037.4/1283.4) then
             nsta = 253
           ELSEIF(br .le. 1048.4/1283.4) then
             nsta = 253
           ELSEIF(br .le. 1068.4/1283.4) then
             nsta = 253
           ELSEIF(br .le. 1079.4/1283.4) then
             nsta = 254
           ELSEIF(br .le. 1092.4/1283.4) then
             nsta = 254
           ELSEIF(br .le. 1102.4/1283.4) then
             nsta = 254
           ELSEIF(br .le. 1117.4/1283.4) then
             nsta = 255
           ELSEIF(br .le. 1130.4/1283.4) then
             nsta = 255
           ELSEIF(br .le. 1157.4/1283.4) then
             nsta = 256
           ELSEIF(br .le. 1165.4/1283.4) then
             nsta = 256
           ELSEIF(br .le. 1182.4/1283.4) then
             nsta = 256
           ELSEIF(br .le. 1195.4/1283.4) then
             nsta = 256
           ELSEIF(br .le. 1215.4/1283.4) then
             nsta = 257
           ELSEIF(br .le. 1233.4/1283.4) then
             nsta = 259
           ELSEIF(br .le. 1248.4/1283.4) then
             nsta = 260
           ELSEIF(br .le. 1268.4/1283.4) then
             nsta = 260
           ELSE
             nsta = 260
           ENDIF
         Case(64158)
               IF(br .le. 6.2/6416.1) then
             nsta = 1
           ELSEIF(br .le. 14.5/6416.1) then
             nsta = 5
           ELSEIF(br .le. 68.5/6416.1) then
             nsta = 6
           ELSEIF(br .le. 1033.5/6416.1) then
             nsta = 10
           ELSEIF(br .le. 1116.5/6416.1) then
             nsta = 14
           ELSEIF(br .le. 1247.5/6416.1) then
             nsta = 27
           ELSEIF(br .le. 1307.5/6416.1) then
             nsta = 39
           ELSEIF(br .le. 1322.5/6416.1) then
             nsta = 44
           ELSEIF(br .le. 1334.5/6416.1) then
             nsta = 45
           ELSEIF(br .le. 1367.5/6416.1) then
             nsta = 48
           ELSEIF(br .le. 1392.4/6416.1) then
             nsta = 52
           ELSEIF(br .le. 1449.4/6416.1) then
             nsta = 54
           ELSEIF(br .le. 1478.4/6416.1) then
             nsta = 64
           ELSEIF(br .le. 1935.4/6416.1) then
             nsta = 65
           ELSEIF(br .le. 1957.4/6416.1) then
             nsta = 68
           ELSEIF(br .le. 1983.4/6416.1) then
             nsta = 70
           ELSEIF(br .le. 1999.4/6416.1) then
             nsta = 71
           ELSEIF(br .le. 2017.4/6416.1) then
             nsta = 73
           ELSEIF(br .le. 2122.4/6416.1) then
             nsta = 75
           ELSEIF(br .le. 2149.4/6416.1) then
             nsta = 79
           ELSEIF(br .le. 2179.4/6416.1) then
             nsta = 80
           ELSEIF(br .le. 2196.4/6416.1) then
             nsta = 82
           ELSEIF(br .le. 2334.4/6416.1) then
             nsta = 83
           ELSEIF(br .le. 2372.4/6416.1) then
             nsta = 85
           ELSEIF(br .le. 2496.4/6416.1) then
             nsta = 86
           ELSEIF(br .le. 2542.4/6416.1) then
             nsta = 87
           ELSEIF(br .le. 2556.4/6416.1) then
             nsta = 91
           ELSEIF(br .le. 2631.4/6416.1) then
             nsta = 93
           ELSEIF(br .le. 2722.4/6416.1) then
             nsta = 96
           ELSEIF(br .le. 2877.4/6416.1) then
             nsta = 97
           ELSEIF(br .le. 2888.9/6416.1) then
             nsta = 98
           ELSEIF(br .le. 3000.9/6416.1) then
             nsta = 100
           ELSEIF(br .le. 3017.9/6416.1) then
             nsta = 102
           ELSEIF(br .le. 3039.9/6416.1) then
             nsta = 104
           ELSEIF(br .le. 3103.9/6416.1) then
             nsta = 109
           ELSEIF(br .le. 3223.9/6416.1) then
             nsta = 110
           ELSEIF(br .le. 3247.9/6416.1) then
             nsta = 112
           ELSEIF(br .le. 3276.9/6416.1) then
             nsta = 113
           ELSEIF(br .le. 3309.9/6416.1) then
             nsta = 115
           ELSEIF(br .le. 3360.9/6416.1) then
             nsta = 117
           ELSEIF(br .le. 3376.1/6416.1) then
             nsta = 120
           ELSEIF(br .le. 3413.1/6416.1) then
             nsta = 121
           ELSEIF(br .le. 3435.1/6416.1) then
             nsta = 122
           ELSEIF(br .le. 3538.1/6416.1) then
             nsta = 125
           ELSEIF(br .le. 3621.1/6416.1) then
             nsta = 126
           ELSEIF(br .le. 3657.1/6416.1) then
             nsta = 129
           ELSEIF(br .le. 3767.1/6416.1) then
             nsta = 130
           ELSEIF(br .le. 3829.1/6416.1) then
             nsta = 133
           ELSEIF(br .le. 3847.1/6416.1) then
             nsta = 134
           ELSEIF(br .le. 3875.1/6416.1) then
             nsta = 137
           ELSEIF(br .le. 3909.1/6416.1) then
             nsta = 138
           ELSEIF(br .le. 3938.1/6416.1) then
             nsta = 141
           ELSEIF(br .le. 4043.1/6416.1) then
             nsta = 144
           ELSEIF(br .le. 4073.1/6416.1) then
             nsta = 146
           ELSEIF(br .le. 4126.1/6416.1) then
             nsta = 147
           ELSEIF(br .le. 4154.1/6416.1) then
             nsta = 149
           ELSEIF(br .le. 4180.1/6416.1) then
             nsta = 152
           ELSEIF(br .le. 4223.1/6416.1) then
             nsta = 153
           ELSEIF(br .le. 4295.1/6416.1) then
             nsta = 154
           ELSEIF(br .le. 4530.1/6416.1) then
             nsta = 156
           ELSEIF(br .le. 4561.1/6416.1) then
             nsta = 157
           ELSEIF(br .le. 4600.1/6416.1) then
             nsta = 159
           ELSEIF(br .le. 4671.1/6416.1) then
             nsta = 160
           ELSEIF(br .le. 4729.1/6416.1) then
             nsta = 162
           ELSEIF(br .le. 4772.1/6416.1) then
             nsta = 168
           ELSEIF(br .le. 4806.1/6416.1) then
             nsta = 169
           ELSEIF(br .le. 4831.1/6416.1) then
             nsta = 171
           ELSEIF(br .le. 4878.1/6416.1) then
             nsta = 173
           ELSEIF(br .le. 4922.1/6416.1) then
             nsta = 174
           ELSEIF(br .le. 4957.1/6416.1) then
             nsta = 175
           ELSEIF(br .le. 4993.1/6416.1) then
             nsta = 177
           ELSEIF(br .le. 5033.1/6416.1) then
             nsta = 179
           ELSEIF(br .le. 5064.1/6416.1) then
             nsta = 180
           ELSEIF(br .le. 5104.1/6416.1) then
             nsta = 181
           ELSEIF(br .le. 5132.1/6416.1) then
             nsta = 184
           ELSEIF(br .le. 5181.1/6416.1) then
             nsta = 185
           ELSEIF(br .le. 5241.1/6416.1) then
             nsta = 187
           ELSEIF(br .le. 5308.1/6416.1) then
             nsta = 188
           ELSEIF(br .le. 5346.1/6416.1) then
             nsta = 190
           ELSEIF(br .le. 5411.1/6416.1) then
             nsta = 193
           ELSEIF(br .le. 5487.1/6416.1) then
             nsta = 195
           ELSEIF(br .le. 5553.1/6416.1) then
             nsta = 196
           ELSEIF(br .le. 5619.1/6416.1) then
             nsta = 197
           ELSEIF(br .le. 5679.1/6416.1) then
             nsta = 198
           ELSEIF(br .le. 5702.1/6416.1) then
             nsta = 199
           ELSEIF(br .le. 5724.1/6416.1) then
             nsta = 200
           ELSEIF(br .le. 5788.1/6416.1) then
             nsta = 201
           ELSEIF(br .le. 5855.1/6416.1) then
             nsta = 203
           ELSEIF(br .le. 5929.1/6416.1) then
             nsta = 204
           ELSEIF(br .le. 5974.1/6416.1) then
             nsta = 206
           ELSEIF(br .le. 6077.1/6416.1) then
             nsta = 207
           ELSEIF(br .le. 6114.1/6416.1) then
             nsta = 208
           ELSEIF(br .le. 6146.1/6416.1) then
             nsta = 209
           ELSEIF(br .le. 6188.1/6416.1) then
             nsta = 210
           ELSEIF(br .le. 6254.1/6416.1) then
             nsta = 211
           ELSEIF(br .le. 6317.1/6416.1) then
             nsta = 212
           ELSE
             nsta = 213
           ENDIF
         Case(65160)
               IF(br .le. 0.0063/2.0876) then
             nsta = 0
           ELSEIF(br .le. 0.0343/2.0876) then
             nsta = 1
           ELSEIF(br .le. 0.0412/2.0876) then
             nsta = 4
           ELSEIF(br .le. 0.0702/2.0876) then
             nsta = 5
           ELSEIF(br .le. 0.1422/2.0876) then
             nsta = 7
           ELSEIF(br .le. 0.1622/2.0876) then
             nsta = 9
           ELSEIF(br .le. 0.3522/2.0876) then
             nsta = 10
           ELSEIF(br .le. 0.3559/2.0876) then
             nsta = 12
           ELSEIF(br .le. 0.4659/2.0876) then
             nsta = 17
           ELSEIF(br .le. 0.4769/2.0876) then
             nsta = 20
           ELSEIF(br .le. 0.4859/2.0876) then
             nsta = 21
           ELSEIF(br .le. 0.5999/2.0876) then
             nsta = 28
           ELSEIF(br .le. 0.7029/2.0876) then
             nsta = 29
           ELSEIF(br .le. 0.7259/2.0876) then
             nsta = 31
           ELSEIF(br .le. 0.8629/2.0876) then
             nsta = 33
           ELSEIF(br .le. 0.8989/2.0876) then
             nsta = 36
           ELSEIF(br .le. 0.9529/2.0876) then
             nsta = 38
           ELSEIF(br .le. 0.9689/2.0876) then
             nsta = 39
           ELSEIF(br .le. 0.9859/2.0876) then
             nsta = 41
           ELSEIF(br .le. 1.0269/2.0876) then
             nsta = 44
           ELSEIF(br .le. 1.1469/2.0876) then
             nsta = 45
           ELSEIF(br .le. 1.1779/2.0876) then
             nsta = 46
           ELSEIF(br .le. 1.2069/2.0876) then
             nsta = 49
           ELSEIF(br .le. 1.2409/2.0876) then
             nsta = 50
           ELSEIF(br .le. 1.2679/2.0876) then
             nsta = 52
           ELSEIF(br .le. 1.3049/2.0876) then
             nsta = 53
           ELSEIF(br .le. 1.3078/2.0876) then
             nsta = 54
           ELSEIF(br .le. 1.3328/2.0876) then
             nsta = 56
           ELSEIF(br .le. 1.3748/2.0876) then
             nsta = 58
           ELSEIF(br .le. 1.3777/2.0876) then
             nsta = 60
           ELSEIF(br .le. 1.3867/2.0876) then
             nsta = 62
           ELSEIF(br .le. 1.4157/2.0876) then
             nsta = 63
           ELSEIF(br .le. 1.4667/2.0876) then
             nsta = 66
           ELSEIF(br .le. 1.4857/2.0876) then
             nsta = 67
           ELSEIF(br .le. 1.4997/2.0876) then
             nsta = 68
           ELSEIF(br .le. 1.5077/2.0876) then
             nsta = 69
           ELSEIF(br .le. 1.5367/2.0876) then
             nsta = 71
           ELSEIF(br .le. 1.5507/2.0876) then
             nsta = 73
           ELSEIF(br .le. 1.5647/2.0876) then
             nsta = 75
           ELSEIF(br .le. 1.5737/2.0876) then
             nsta = 76
           ELSEIF(br .le. 1.5977/2.0876) then
             nsta = 77
           ELSEIF(br .le. 1.6137/2.0876) then
             nsta = 78
           ELSEIF(br .le. 1.6166/2.0876) then
             nsta = 79
           ELSEIF(br .le. 1.6196/2.0876) then
             nsta = 80
           ELSEIF(br .le. 1.6326/2.0876) then
             nsta = 81
           ELSEIF(br .le. 1.6486/2.0876) then
             nsta = 82
           ELSEIF(br .le. 1.6696/2.0876) then
             nsta = 84
           ELSEIF(br .le. 1.6966/2.0876) then
             nsta = 86
           ELSEIF(br .le. 1.7096/2.0876) then
             nsta = 87
           ELSEIF(br .le. 1.7316/2.0876) then
             nsta = 88
           ELSEIF(br .le. 1.7956/2.0876) then
             nsta = 89
           ELSEIF(br .le. 1.8566/2.0876) then
             nsta = 90
           ELSEIF(br .le. 1.8826/2.0876) then
             nsta = 91
           ELSEIF(br .le. 1.9346/2.0876) then
             nsta = 92
           ELSEIF(br .le. 1.9746/2.0876) then
             nsta = 93
           ELSEIF(br .le. 1.9916/2.0876) then
             nsta = 94
           ELSEIF(br .le. 2.0316/2.0876) then
             nsta = 96
           ELSEIF(br .le. 2.0646/2.0876) then
             nsta = 97
           ELSE
             nsta = 99
           ENDIF
         Case(66161)
               IF(br .le. 0.0083/1.9041) then
             nsta = 3
           ELSEIF(br .le. 0.8583/1.9041) then
             nsta = 13
           ELSEIF(br .le. 0.9833/1.9041) then
             nsta = 15
           ELSEIF(br .le. 1.0513/1.9041) then
             nsta = 28
           ELSEIF(br .le. 1.2013/1.9041) then
             nsta = 44
           ELSEIF(br .le. 1.2153/1.9041) then
             nsta = 48
           ELSEIF(br .le. 1.2423/1.9041) then
             nsta = 55
           ELSEIF(br .le. 1.2674/1.9041) then
             nsta = 74
           ELSEIF(br .le. 1.3013/1.9041) then
             nsta = 75
           ELSEIF(br .le. 1.3323/1.9041) then
             nsta = 78
           ELSEIF(br .le. 1.3379/1.9041) then
             nsta = 79
           ELSEIF(br .le. 1.3492/1.9041) then
             nsta = 80
           ELSEIF(br .le. 1.3683/1.9041) then
             nsta = 81
           ELSEIF(br .le. 1.3771/1.9041) then
             nsta = 83
           ELSEIF(br .le. 1.3837/1.9041) then
             nsta = 85
           ELSEIF(br .le. 1.4135/1.9041) then
             nsta = 86
           ELSEIF(br .le. 1.427/1.9041) then
             nsta = 88
           ELSEIF(br .le. 1.494/1.9041) then
             nsta = 98
           ELSEIF(br .le. 1.577/1.9041) then
             nsta = 99
           ELSEIF(br .le. 1.599/1.9041) then
             nsta = 100
           ELSEIF(br .le. 1.611/1.9041) then
             nsta = 100
           ELSEIF(br .le. 1.814/1.9041) then
             nsta = 100
           ELSEIF(br .le. 1.88/1.9041) then
             nsta = 103
           ELSE
             nsta = 104
           ENDIF
         Case(66163)
               IF(br .le. 0.0007/13.9221) then
             nsta = 3
           ELSEIF(br .le. 0.0247/13.9221) then
             nsta = 7
           ELSEIF(br .le. 0.5047/13.9221) then
             nsta = 8
           ELSEIF(br .le. 0.9847/13.9221) then
             nsta = 11
           ELSEIF(br .le. 1.1197/13.9221) then
             nsta = 31
           ELSEIF(br .le. 1.2127/13.9221) then
             nsta = 33
           ELSEIF(br .le. 1.2178/13.9221) then
             nsta = 34
           ELSEIF(br .le. 1.2356/13.9221) then
             nsta = 35
           ELSEIF(br .le. 2.3356/13.9221) then
             nsta = 37
           ELSEIF(br .le. 2.4136/13.9221) then
             nsta = 41
           ELSEIF(br .le. 2.4966/13.9221) then
             nsta = 43
           ELSEIF(br .le. 2.6126/13.9221) then
             nsta = 49
           ELSEIF(br .le. 2.6131/13.9221) then
             nsta = 51
           ELSEIF(br .le. 3.0731/13.9221) then
             nsta = 61
           ELSEIF(br .le. 3.4531/13.9221) then
             nsta = 62
           ELSEIF(br .le. 3.5191/13.9221) then
             nsta = 66
           ELSEIF(br .le. 3.6461/13.9221) then
             nsta = 75
           ELSEIF(br .le. 3.7911/13.9221) then
             nsta = 77
           ELSEIF(br .le. 4.0611/13.9221) then
             nsta = 79
           ELSEIF(br .le. 4.1381/13.9221) then
             nsta = 96
           ELSEIF(br .le. 4.4481/13.9221) then
             nsta = 100
           ELSEIF(br .le. 4.5191/13.9221) then
             nsta = 103
           ELSEIF(br .le. 4.5501/13.9221) then
             nsta = 109
           ELSEIF(br .le. 4.6401/13.9221) then
             nsta = 112
           ELSEIF(br .le. 4.9001/13.9221) then
             nsta = 117
           ELSEIF(br .le. 5.0001/13.9221) then
             nsta = 127
           ELSEIF(br .le. 5.0821/13.9221) then
             nsta = 131
           ELSEIF(br .le. 5.8521/13.9221) then
             nsta = 136
           ELSEIF(br .le. 5.9921/13.9221) then
             nsta = 151
           ELSEIF(br .le. 6.0421/13.9221) then
             nsta = 154
           ELSEIF(br .le. 6.4421/13.9221) then
             nsta = 161
           ELSEIF(br .le. 6.5621/13.9221) then
             nsta = 162
           ELSEIF(br .le. 6.8321/13.9221) then
             nsta = 166
           ELSEIF(br .le. 7.3221/13.9221) then
             nsta = 169
           ELSE
             nsta = 313
           ENDIF
         Case(66164)
               IF(br .le. 0.066/1.552) then
             nsta = 1
           ELSEIF(br .le. 0.164/1.552) then
             nsta = 2
           ELSEIF(br .le. 0.216/1.552) then
             nsta = 4
           ELSEIF(br .le. 0.356/1.552) then
             nsta = 5
           ELSEIF(br .le. 0.388/1.552) then
             nsta = 8
           ELSEIF(br .le. 0.408/1.552) then
             nsta = 10
           ELSEIF(br .le. 0.477/1.552) then
             nsta = 27
           ELSEIF(br .le. 0.727/1.552) then
             nsta = 35
           ELSEIF(br .le. 0.753/1.552) then
             nsta = 38
           ELSEIF(br .le. 0.913/1.552) then
             nsta = 43
           ELSEIF(br .le. 1.133/1.552) then
             nsta = 45
           ELSEIF(br .le. 1.323/1.552) then
             nsta = 49
           ELSEIF(br .le. 1.379/1.552) then
             nsta = 70
           ELSEIF(br .le. 1.422/1.552) then
             nsta = 74
           ELSE
             nsta = 86
           ENDIF
         Case(66165)
               IF(br .le. 35.9/242.235) then
             nsta = 2
           ELSEIF(br .le. 64.6/242.235) then
             nsta = 3
           ELSEIF(br .le. 71.2/242.235) then
             nsta = 17
           ELSEIF(br .le. 79.6/242.235) then
             nsta = 18
           ELSEIF(br .le. 95.3/242.235) then
             nsta = 19
           ELSEIF(br .le. 95.316/242.235) then
             nsta = 20
           ELSEIF(br .le. 101.416/242.235) then
             nsta = 21
           ELSEIF(br .le. 101.545/242.235) then
             nsta = 36
           ELSEIF(br .le. 101.805/242.235) then
             nsta = 41
           ELSEIF(br .le. 104.405/242.235) then
             nsta = 45
           ELSEIF(br .le. 110.105/242.235) then
             nsta = 47
           ELSEIF(br .le. 112.005/242.235) then
             nsta = 48
           ELSEIF(br .le. 112.235/242.235) then
             nsta = 51
           ELSEIF(br .le. 113.835/242.235) then
             nsta = 57
           ELSEIF(br .le. 115.235/242.235) then
             nsta = 66
           ELSEIF(br .le. 115.545/242.235) then
             nsta = 67
           ELSEIF(br .le. 116.755/242.235) then
             nsta = 69
           ELSEIF(br .le. 117.705/242.235) then
             nsta = 71
           ELSEIF(br .le. 117.862/242.235) then
             nsta = 74
           ELSEIF(br .le. 119.262/242.235) then
             nsta = 76
           ELSEIF(br .le. 119.409/242.235) then
             nsta = 80
           ELSEIF(br .le. 119.575/242.235) then
             nsta = 83
           ELSEIF(br .le. 121.675/242.235) then
             nsta = 84
           ELSEIF(br .le. 134.775/242.235) then
             nsta = 86
           ELSEIF(br .le. 135.295/242.235) then
             nsta = 88
           ELSEIF(br .le. 139.595/242.235) then
             nsta = 89
           ELSEIF(br .le. 140.235/242.235) then
             nsta = 90
           ELSEIF(br .le. 142.735/242.235) then
             nsta = 92
           ELSEIF(br .le. 147.435/242.235) then
             nsta = 97
           ELSEIF(br .le. 149.635/242.235) then
             nsta = 98
           ELSEIF(br .le. 154.835/242.235) then
             nsta = 103
           ELSEIF(br .le. 159.735/242.235) then
             nsta = 106
           ELSEIF(br .le. 162.835/242.235) then
             nsta = 111
           ELSEIF(br .le. 164.735/242.235) then
             nsta = 120
           ELSEIF(br .le. 167.835/242.235) then
             nsta = 121
           ELSEIF(br .le. 172.535/242.235) then
             nsta = 123
           ELSEIF(br .le. 175.735/242.235) then
             nsta = 124
           ELSEIF(br .le. 186.335/242.235) then
             nsta = 130
           ELSEIF(br .le. 191.035/242.235) then
             nsta = 137
           ELSEIF(br .le. 194.335/242.235) then
             nsta = 139
           ELSEIF(br .le. 200.135/242.235) then
             nsta = 142
           ELSEIF(br .le. 204.635/242.235) then
             nsta = 146
           ELSEIF(br .le. 207.335/242.235) then
             nsta = 147
           ELSEIF(br .le. 212.435/242.235) then
             nsta = 151
           ELSEIF(br .le. 216.235/242.235) then
             nsta = 152
           ELSEIF(br .le. 223.135/242.235) then
             nsta = 158
           ELSEIF(br .le. 226.135/242.235) then
             nsta = 160
           ELSEIF(br .le. 229.735/242.235) then
             nsta = 165
           ELSEIF(br .le. 232.535/242.235) then
             nsta = 166
           ELSEIF(br .le. 235.835/242.235) then
             nsta = 167
           ELSEIF(br .le. 238.435/242.235) then
             nsta = 171
           ELSE
             nsta = 172
           ENDIF
         Case(67166)
               IF(br .le. 0.0032/9.8189) then
             nsta = 2
           ELSEIF(br .le. 0.0322/9.8189) then
             nsta = 5
           ELSEIF(br .le. 0.0396/9.8189) then
             nsta = 6
           ELSEIF(br .le. 0.2276/9.8189) then
             nsta = 7
           ELSEIF(br .le. 0.3206/9.8189) then
             nsta = 8
           ELSEIF(br .le. 0.3238/9.8189) then
             nsta = 13
           ELSEIF(br .le. 0.328/9.8189) then
             nsta = 14
           ELSEIF(br .le. 0.524/9.8189) then
             nsta = 15
           ELSEIF(br .le. 0.748/9.8189) then
             nsta = 16
           ELSEIF(br .le. 0.7612/9.8189) then
             nsta = 21
           ELSEIF(br .le. 0.7644/9.8189) then
             nsta = 21
           ELSEIF(br .le. 1.3044/9.8189) then
             nsta = 24
           ELSEIF(br .le. 1.3197/9.8189) then
             nsta = 27
           ELSEIF(br .le. 1.4097/9.8189) then
             nsta = 28
           ELSEIF(br .le. 1.4467/9.8189) then
             nsta = 29
           ELSEIF(br .le. 1.5537/9.8189) then
             nsta = 30
           ELSEIF(br .le. 1.5827/9.8189) then
             nsta = 32
           ELSEIF(br .le. 1.6001/9.8189) then
             nsta = 34
           ELSEIF(br .le. 1.6381/9.8189) then
             nsta = 35
           ELSEIF(br .le. 1.7351/9.8189) then
             nsta = 37
           ELSEIF(br .le. 1.7562/9.8189) then
             nsta = 38
           ELSEIF(br .le. 1.8092/9.8189) then
             nsta = 42
           ELSEIF(br .le. 1.8472/9.8189) then
             nsta = 45
           ELSEIF(br .le. 1.8498/9.8189) then
             nsta = 46
           ELSEIF(br .le. 1.8524/9.8189) then
             nsta = 47
           ELSEIF(br .le. 1.8566/9.8189) then
             nsta = 48
           ELSEIF(br .le. 1.883/9.8189) then
             nsta = 49
           ELSEIF(br .le. 1.8909/9.8189) then
             nsta = 52
           ELSEIF(br .le. 1.9094/9.8189) then
             nsta = 53
           ELSEIF(br .le. 1.93/9.8189) then
             nsta = 55
           ELSEIF(br .le. 1.9543/9.8189) then
             nsta = 56
           ELSEIF(br .le. 1.9622/9.8189) then
             nsta = 58
           ELSEIF(br .le. 2.0132/9.8189) then
             nsta = 59
           ELSEIF(br .le. 2.0158/9.8189) then
             nsta = 62
           ELSEIF(br .le. 2.2078/9.8189) then
             nsta = 65
           ELSEIF(br .le. 2.2099/9.8189) then
             nsta = 68
           ELSEIF(br .le. 2.2157/9.8189) then
             nsta = 70
           ELSEIF(br .le. 2.232/9.8189) then
             nsta = 71
           ELSEIF(br .le. 2.2436/9.8189) then
             nsta = 72
           ELSEIF(br .le. 2.2816/9.8189) then
             nsta = 74
           ELSEIF(br .le. 2.2832/9.8189) then
             nsta = 78
           ELSEIF(br .le. 2.2937/9.8189) then
             nsta = 80
           ELSEIF(br .le. 2.5167/9.8189) then
             nsta = 82
           ELSEIF(br .le. 2.5487/9.8189) then
             nsta = 84
           ELSEIF(br .le. 2.5787/9.8189) then
             nsta = 85
           ELSEIF(br .le. 2.5919/9.8189) then
             nsta = 92
           ELSEIF(br .le. 2.5966/9.8189) then
             nsta = 93
           ELSEIF(br .le. 2.6696/9.8189) then
             nsta = 95
           ELSEIF(br .le. 2.678/9.8189) then
             nsta = 97
           ELSEIF(br .le. 2.72/9.8189) then
             nsta = 98
           ELSEIF(br .le. 2.83/9.8189) then
             nsta = 101
           ELSEIF(br .le. 2.8316/9.8189) then
             nsta = 104
           ELSEIF(br .le. 2.8479/9.8189) then
             nsta = 108
           ELSEIF(br .le. 2.8526/9.8189) then
             nsta = 109
           ELSEIF(br .le. 2.8769/9.8189) then
             nsta = 111
           ELSEIF(br .le. 2.8811/9.8189) then
             nsta = 112
           ELSEIF(br .le. 2.8843/9.8189) then
             nsta = 114
           ELSEIF(br .le. 2.8933/9.8189) then
             nsta = 115
           ELSEIF(br .le. 2.9513/9.8189) then
             nsta = 117
           ELSEIF(br .le. 2.9587/9.8189) then
             nsta = 119
           ELSEIF(br .le. 2.9761/9.8189) then
             nsta = 120
           ELSEIF(br .le. 2.9814/9.8189) then
             nsta = 121
           ELSEIF(br .le. 2.9825/9.8189) then
             nsta = 122
           ELSEIF(br .le. 2.9878/9.8189) then
             nsta = 123
           ELSEIF(br .le. 3.2478/9.8189) then
             nsta = 126
           ELSEIF(br .le. 3.2648/9.8189) then
             nsta = 130
           ELSEIF(br .le. 3.5178/9.8189) then
             nsta = 132
           ELSEIF(br .le. 3.5738/9.8189) then
             nsta = 135
           ELSEIF(br .le. 3.6318/9.8189) then
             nsta = 138
           ELSEIF(br .le. 3.8028/9.8189) then
             nsta = 139
           ELSEIF(br .le. 3.8368/9.8189) then
             nsta = 141
           ELSEIF(br .le. 3.8426/9.8189) then
             nsta = 142
           ELSEIF(br .le. 4.1726/9.8189) then
             nsta = 143
           ELSEIF(br .le. 4.2046/9.8189) then
             nsta = 144
           ELSEIF(br .le. 4.212/9.8189) then
             nsta = 145
           ELSEIF(br .le. 4.2157/9.8189) then
             nsta = 147
           ELSEIF(br .le. 4.2547/9.8189) then
             nsta = 149
           ELSEIF(br .le. 4.3927/9.8189) then
             nsta = 151
           ELSEIF(br .le. 4.3948/9.8189) then
             nsta = 153
           ELSEIF(br .le. 4.4768/9.8189) then
             nsta = 154
           ELSEIF(br .le. 4.4789/9.8189) then
             nsta = 155
           ELSEIF(br .le. 4.5089/9.8189) then
             nsta = 156
           ELSEIF(br .le. 4.5258/9.8189) then
             nsta = 157
           ELSEIF(br .le. 4.5432/9.8189) then
             nsta = 158
           ELSEIF(br .le. 4.5569/9.8189) then
             nsta = 159
           ELSEIF(br .le. 4.569/9.8189) then
             nsta = 160
           ELSEIF(br .le. 4.652/9.8189) then
             nsta = 162
           ELSEIF(br .le. 4.6647/9.8189) then
             nsta = 163
           ELSEIF(br .le. 4.7117/9.8189) then
             nsta = 164
           ELSEIF(br .le. 4.7787/9.8189) then
             nsta = 165
           ELSEIF(br .le. 4.7957/9.8189) then
             nsta = 166
           ELSEIF(br .le. 4.8094/9.8189) then
             nsta = 167
           ELSEIF(br .le. 4.8242/9.8189) then
             nsta = 168
           ELSEIF(br .le. 4.8752/9.8189) then
             nsta = 169
           ELSEIF(br .le. 4.8826/9.8189) then
             nsta = 170
           ELSEIF(br .le. 4.8984/9.8189) then
             nsta = 173
           ELSEIF(br .le. 4.9344/9.8189) then
             nsta = 174
           ELSEIF(br .le. 5.1484/9.8189) then
             nsta = 175
           ELSEIF(br .le. 5.2134/9.8189) then
             nsta = 176
           ELSEIF(br .le. 5.2297/9.8189) then
             nsta = 177
           ELSEIF(br .le. 5.2397/9.8189) then
             nsta = 178
           ELSEIF(br .le. 5.2455/9.8189) then
             nsta = 179
           ELSEIF(br .le. 5.2529/9.8189) then
             nsta = 180
           ELSEIF(br .le. 5.2899/9.8189) then
             nsta = 181
           ELSEIF(br .le. 5.2925/9.8189) then
             nsta = 182
           ELSEIF(br .le. 5.3855/9.8189) then
             nsta = 184
           ELSEIF(br .le. 5.4103/9.8189) then
             nsta = 185
           ELSEIF(br .le. 5.4703/9.8189) then
             nsta = 186
           ELSEIF(br .le. 5.5333/9.8189) then
             nsta = 187
           ELSEIF(br .le. 5.5486/9.8189) then
             nsta = 189
           ELSEIF(br .le. 5.5613/9.8189) then
             nsta = 190
           ELSEIF(br .le. 5.5629/9.8189) then
             nsta = 191
           ELSEIF(br .le. 5.6459/9.8189) then
             nsta = 192
           ELSEIF(br .le. 5.6669/9.8189) then
             nsta = 193
           ELSEIF(br .le. 5.8129/9.8189) then
             nsta = 194
           ELSEIF(br .le. 5.8319/9.8189) then
             nsta = 195
           ELSEIF(br .le. 5.8649/9.8189) then
             nsta = 196
           ELSEIF(br .le. 5.8834/9.8189) then
             nsta = 197
           ELSEIF(br .le. 5.8887/9.8189) then
             nsta = 198
           ELSEIF(br .le. 6.0187/9.8189) then
             nsta = 199
           ELSEIF(br .le. 6.0308/9.8189) then
             nsta = 200
           ELSEIF(br .le. 6.0898/9.8189) then
             nsta = 201
           ELSEIF(br .le. 6.1083/9.8189) then
             nsta = 202
           ELSEIF(br .le. 6.1723/9.8189) then
             nsta = 204
           ELSEIF(br .le. 6.1781/9.8189) then
             nsta = 205
           ELSEIF(br .le. 6.1901/9.8189) then
             nsta = 206
           ELSEIF(br .le. 6.2711/9.8189) then
             nsta = 207
           ELSEIF(br .le. 6.2753/9.8189) then
             nsta = 208
           ELSEIF(br .le. 6.2837/9.8189) then
             nsta = 209
           ELSEIF(br .le. 6.2895/9.8189) then
             nsta = 210
           ELSEIF(br .le. 6.3053/9.8189) then
             nsta = 211
           ELSEIF(br .le. 6.318/9.8189) then
             nsta = 212
           ELSEIF(br .le. 6.3301/9.8189) then
             nsta = 213
           ELSEIF(br .le. 6.3354/9.8189) then
             nsta = 214
           ELSEIF(br .le. 6.3412/9.8189) then
             nsta = 215
           ELSEIF(br .le. 6.3852/9.8189) then
             nsta = 216
           ELSEIF(br .le. 6.3889/9.8189) then
             nsta = 217
           ELSEIF(br .le. 6.3984/9.8189) then
             nsta = 218
           ELSEIF(br .le. 6.4594/9.8189) then
             nsta = 219
           ELSEIF(br .le. 6.4868/9.8189) then
             nsta = 220
           ELSEIF(br .le. 6.5238/9.8189) then
             nsta = 221
           ELSEIF(br .le. 6.5264/9.8189) then
             nsta = 222
           ELSEIF(br .le. 6.5784/9.8189) then
             nsta = 223
           ELSEIF(br .le. 6.5916/9.8189) then
             nsta = 224
           ELSEIF(br .le. 6.6426/9.8189) then
             nsta = 225
           ELSEIF(br .le. 6.6556/9.8189) then
             nsta = 227
           ELSEIF(br .le. 6.6593/9.8189) then
             nsta = 228
           ELSEIF(br .le. 6.7363/9.8189) then
             nsta = 229
           ELSEIF(br .le. 6.7693/9.8189) then
             nsta = 230
           ELSEIF(br .le. 6.7856/9.8189) then
             nsta = 231
           ELSEIF(br .le. 6.7972/9.8189) then
             nsta = 232
           ELSEIF(br .le. 6.8572/9.8189) then
             nsta = 233
           ELSEIF(br .le. 6.9022/9.8189) then
             nsta = 234
           ELSEIF(br .le. 6.9522/9.8189) then
             nsta = 235
           ELSEIF(br .le. 6.9632/9.8189) then
             nsta = 236
           ELSEIF(br .le. 7.0192/9.8189) then
             nsta = 237
           ELSEIF(br .le. 7.025/9.8189) then
             nsta = 238
           ELSEIF(br .le. 7.082/9.8189) then
             nsta = 239
           ELSEIF(br .le. 7.107/9.8189) then
             nsta = 240
           ELSEIF(br .le. 7.175/9.8189) then
             nsta = 241
           ELSEIF(br .le. 7.1898/9.8189) then
             nsta = 242
           ELSEIF(br .le. 7.2003/9.8189) then
             nsta = 243
           ELSEIF(br .le. 7.2061/9.8189) then
             nsta = 244
           ELSEIF(br .le. 7.2188/9.8189) then
             nsta = 245
           ELSEIF(br .le. 7.2283/9.8189) then
             nsta = 246
           ELSEIF(br .le. 7.2573/9.8189) then
             nsta = 247
           ELSEIF(br .le. 7.2913/9.8189) then
             nsta = 248
           ELSEIF(br .le. 7.3082/9.8189) then
             nsta = 249
           ELSEIF(br .le. 7.3156/9.8189) then
             nsta = 250
           ELSEIF(br .le. 7.3256/9.8189) then
             nsta = 251
           ELSEIF(br .le. 7.3293/9.8189) then
             nsta = 252
           ELSEIF(br .le. 7.3703/9.8189) then
             nsta = 254
           ELSEIF(br .le. 7.3983/9.8189) then
             nsta = 255
           ELSEIF(br .le. 7.4083/9.8189) then
             nsta = 256
           ELSEIF(br .le. 7.4268/9.8189) then
             nsta = 257
           ELSEIF(br .le. 7.4379/9.8189) then
             nsta = 258
           ELSEIF(br .le. 7.44/9.8189) then
             nsta = 259
           ELSEIF(br .le. 7.4585/9.8189) then
             nsta = 260
           ELSEIF(br .le. 7.4733/9.8189) then
             nsta = 261
           ELSEIF(br .le. 7.4791/9.8189) then
             nsta = 262
           ELSEIF(br .le. 7.4838/9.8189) then
             nsta = 263
           ELSEIF(br .le. 7.4943/9.8189) then
             nsta = 264
           ELSEIF(br .le. 7.5523/9.8189) then
             nsta = 265
           ELSEIF(br .le. 7.5618/9.8189) then
             nsta = 266
           ELSEIF(br .le. 7.6778/9.8189) then
             nsta = 267
           ELSEIF(br .le. 7.6836/9.8189) then
             nsta = 268
           ELSEIF(br .le. 7.6978/9.8189) then
             nsta = 270
           ELSEIF(br .le. 7.7025/9.8189) then
             nsta = 271
           ELSEIF(br .le. 7.7083/9.8189) then
             nsta = 272
           ELSEIF(br .le. 7.7463/9.8189) then
             nsta = 273
           ELSEIF(br .le. 7.7793/9.8189) then
             nsta = 274
           ELSEIF(br .le. 7.7883/9.8189) then
             nsta = 275
           ELSEIF(br .le. 7.8183/9.8189) then
             nsta = 277
           ELSEIF(br .le. 7.8209/9.8189) then
             nsta = 278
           ELSEIF(br .le. 7.8709/9.8189) then
             nsta = 279
           ELSEIF(br .le. 7.8783/9.8189) then
             nsta = 280
           ELSEIF(br .le. 7.9036/9.8189) then
             nsta = 281
           ELSEIF(br .le. 7.9426/9.8189) then
             nsta = 282
           ELSEIF(br .le. 7.9452/9.8189) then
             nsta = 283
           ELSEIF(br .le. 7.9499/9.8189) then
             nsta = 284
           ELSEIF(br .le. 7.9839/9.8189) then
             nsta = 285
           ELSEIF(br .le. 7.9939/9.8189) then
             nsta = 286
           ELSEIF(br .le. 8.0187/9.8189) then
             nsta = 287
           ELSEIF(br .le. 8.1027/9.8189) then
             nsta = 288
           ELSEIF(br .le. 8.1387/9.8189) then
             nsta = 289
           ELSEIF(br .le. 8.1767/9.8189) then
             nsta = 290
           ELSEIF(br .le. 8.1878/9.8189) then
             nsta = 291
           ELSEIF(br .le. 8.2218/9.8189) then
             nsta = 292
           ELSEIF(br .le. 8.2678/9.8189) then
             nsta = 293
           ELSEIF(br .le. 8.3008/9.8189) then
             nsta = 294
           ELSEIF(br .le. 8.3135/9.8189) then
             nsta = 295
           ELSEIF(br .le. 8.3445/9.8189) then
             nsta = 296
           ELSEIF(br .le. 8.3582/9.8189) then
             nsta = 297
           ELSEIF(br .le. 8.4202/9.8189) then
             nsta = 298
           ELSEIF(br .le. 8.4239/9.8189) then
             nsta = 299
           ELSEIF(br .le. 8.4439/9.8189) then
             nsta = 300
           ELSEIF(br .le. 8.4839/9.8189) then
             nsta = 301
           ELSEIF(br .le. 8.4944/9.8189) then
             nsta = 302
           ELSEIF(br .le. 8.4991/9.8189) then
             nsta = 303
           ELSEIF(br .le. 8.5096/9.8189) then
             nsta = 304
           ELSEIF(br .le. 8.536/9.8189) then
             nsta = 305
           ELSEIF(br .le. 8.556/9.8189) then
             nsta = 306
           ELSEIF(br .le. 8.5755/9.8189) then
             nsta = 307
           ELSEIF(br .le. 8.5892/9.8189) then
             nsta = 309
           ELSEIF(br .le. 8.6762/9.8189) then
             nsta = 310
           ELSEIF(br .le. 8.701/9.8189) then
             nsta = 311
           ELSEIF(br .le. 8.735/9.8189) then
             nsta = 312
           ELSEIF(br .le. 8.771/9.8189) then
             nsta = 313
           ELSEIF(br .le. 8.7837/9.8189) then
             nsta = 314
           ELSEIF(br .le. 8.8287/9.8189) then
             nsta = 315
           ELSEIF(br .le. 8.8497/9.8189) then
             nsta = 316
           ELSEIF(br .le. 8.8907/9.8189) then
             nsta = 317
           ELSEIF(br .le. 8.9247/9.8189) then
             nsta = 318
           ELSEIF(br .le. 8.9374/9.8189) then
             nsta = 319
           ELSEIF(br .le. 8.9448/9.8189) then
             nsta = 320
           ELSEIF(br .le. 8.9688/9.8189) then
             nsta = 321
           ELSEIF(br .le. 8.9888/9.8189) then
             nsta = 322
           ELSEIF(br .le. 9.0268/9.8189) then
             nsta = 323
           ELSEIF(br .le. 9.0331/9.8189) then
             nsta = 324
           ELSEIF(br .le. 9.0841/9.8189) then
             nsta = 325
           ELSEIF(br .le. 9.1181/9.8189) then
             nsta = 326
           ELSEIF(br .le. 9.1286/9.8189) then
             nsta = 327
           ELSEIF(br .le. 9.1686/9.8189) then
             nsta = 328
           ELSEIF(br .le. 9.1923/9.8189) then
             nsta = 329
           ELSEIF(br .le. 9.2192/9.8189) then
             nsta = 330
           ELSEIF(br .le. 9.2271/9.8189) then
             nsta = 331
           ELSEIF(br .le. 9.2524/9.8189) then
             nsta = 332
           ELSEIF(br .le. 9.2664/9.8189) then
             nsta = 333
           ELSEIF(br .le. 9.2774/9.8189) then
             nsta = 334
           ELSEIF(br .le. 9.2837/9.8189) then
             nsta = 335
           ELSEIF(br .le. 9.2927/9.8189) then
             nsta = 336
           ELSEIF(br .le. 9.3143/9.8189) then
             nsta = 337
           ELSEIF(br .le. 9.3206/9.8189) then
             nsta = 338
           ELSEIF(br .le. 9.3285/9.8189) then
             nsta = 339
           ELSEIF(br .le. 9.347/9.8189) then
             nsta = 340
           ELSEIF(br .le. 9.3713/9.8189) then
             nsta = 341
           ELSEIF(br .le. 9.3863/9.8189) then
             nsta = 342
           ELSEIF(br .le. 9.3993/9.8189) then
             nsta = 343
           ELSEIF(br .le. 9.4373/9.8189) then
             nsta = 344
           ELSEIF(br .le. 9.4531/9.8189) then
             nsta = 345
           ELSEIF(br .le. 9.4991/9.8189) then
             nsta = 346
           ELSEIF(br .le. 9.5461/9.8189) then
             nsta = 347
           ELSEIF(br .le. 9.5619/9.8189) then
             nsta = 348
           ELSEIF(br .le. 9.5851/9.8189) then
             nsta = 349
           ELSEIF(br .le. 9.6001/9.8189) then
             nsta = 350
           ELSEIF(br .le. 9.6541/9.8189) then
             nsta = 351
           ELSEIF(br .le. 9.6641/9.8189) then
             nsta = 352
           ELSEIF(br .le. 9.6794/9.8189) then
             nsta = 353
           ELSEIF(br .le. 9.6989/9.8189) then
             nsta = 354
           ELSE
             nsta = 355
           ENDIF
         Case(68167)
               IF(br .le. 1.41/3.2816) then
             nsta = 3
           ELSEIF(br .le. 1.76/3.2816) then
             nsta = 4
           ELSEIF(br .le. 1.7735/3.2816) then
             nsta = 12
           ELSEIF(br .le. 1.8095/3.2816) then
             nsta = 27
           ELSEIF(br .le. 1.8185/3.2816) then
             nsta = 28
           ELSEIF(br .le. 1.8365/3.2816) then
             nsta = 31
           ELSEIF(br .le. 1.8965/3.2816) then
             nsta = 56
           ELSEIF(br .le. 1.9645/3.2816) then
             nsta = 63
           ELSEIF(br .le. 2.2645/3.2816) then
             nsta = 71
           ELSEIF(br .le. 2.3215/3.2816) then
             nsta = 85
           ELSEIF(br .le. 2.3362/3.2816) then
             nsta = 97
           ELSEIF(br .le. 2.3612/3.2816) then
             nsta = 101
           ELSEIF(br .le. 2.3832/3.2816) then
             nsta = 107
           ELSEIF(br .le. 2.4432/3.2816) then
             nsta = 109
           ELSEIF(br .le. 2.5332/3.2816) then
             nsta = 111
           ELSEIF(br .le. 2.7432/3.2816) then
             nsta = 115
           ELSEIF(br .le. 2.8062/3.2816) then
             nsta = 118
           ELSEIF(br .le. 2.8152/3.2816) then
             nsta = 122
           ELSEIF(br .le. 2.8342/3.2816) then
             nsta = 124
           ELSEIF(br .le. 2.9442/3.2816) then
             nsta = 131
           ELSEIF(br .le. 2.9556/3.2816) then
             nsta = 135
           ELSEIF(br .le. 3.1956/3.2816) then
             nsta = 138
           ELSEIF(br .le. 3.2316/3.2816) then
             nsta = 146
           ELSEIF(br .le. 3.2456/3.2816) then
             nsta = 147
           ELSE
             nsta = 148
           ENDIF
         Case(68168)
               IF(br .le. 0.115/62.218) then
             nsta = 1
           ELSEIF(br .le. 0.355/62.218) then
             nsta = 2
           ELSEIF(br .le. 0.474/62.218) then
             nsta = 3
           ELSEIF(br .le. 0.784/62.218) then
             nsta = 4
           ELSEIF(br .le. 0.871/62.218) then
             nsta = 5
           ELSEIF(br .le. 0.982/62.218) then
             nsta = 7
           ELSEIF(br .le. 2.002/62.218) then
             nsta = 8
           ELSEIF(br .le. 2.075/62.218) then
             nsta = 9
           ELSEIF(br .le. 2.373/62.218) then
             nsta = 10
           ELSEIF(br .le. 2.526/62.218) then
             nsta = 14
           ELSEIF(br .le. 2.806/62.218) then
             nsta = 18
           ELSEIF(br .le. 4.106/62.218) then
             nsta = 19
           ELSEIF(br .le. 4.376/62.218) then
             nsta = 21
           ELSEIF(br .le. 4.428/62.218) then
             nsta = 24
           ELSEIF(br .le. 5.968/62.218) then
             nsta = 26
           ELSEIF(br .le. 6.698/62.218) then
             nsta = 27
           ELSEIF(br .le. 6.988/62.218) then
             nsta = 28
           ELSEIF(br .le. 8.488/62.218) then
             nsta = 30
           ELSEIF(br .le. 9.058/62.218) then
             nsta = 34
           ELSEIF(br .le. 9.348/62.218) then
             nsta = 35
           ELSEIF(br .le. 9.498/62.218) then
             nsta = 36
           ELSEIF(br .le. 9.678/62.218) then
             nsta = 37
           ELSEIF(br .le. 10.028/62.218) then
             nsta = 38
           ELSEIF(br .le. 10.438/62.218) then
             nsta = 39
           ELSEIF(br .le. 11.308/62.218) then
             nsta = 49
           ELSEIF(br .le. 12.258/62.218) then
             nsta = 50
           ELSEIF(br .le. 12.408/62.218) then
             nsta = 52
           ELSEIF(br .le. 12.461/62.218) then
             nsta = 53
           ELSEIF(br .le. 13.241/62.218) then
             nsta = 56
           ELSEIF(br .le. 14.011/62.218) then
             nsta = 59
           ELSEIF(br .le. 14.351/62.218) then
             nsta = 60
           ELSEIF(br .le. 14.439/62.218) then
             nsta = 62
           ELSEIF(br .le. 14.619/62.218) then
             nsta = 69
           ELSEIF(br .le. 15.039/62.218) then
             nsta = 71
           ELSEIF(br .le. 15.249/62.218) then
             nsta = 72
           ELSEIF(br .le. 15.599/62.218) then
             nsta = 73
           ELSEIF(br .le. 16.039/62.218) then
             nsta = 75
           ELSEIF(br .le. 16.189/62.218) then
             nsta = 76
           ELSEIF(br .le. 16.339/62.218) then
             nsta = 77
           ELSEIF(br .le. 16.539/62.218) then
             nsta = 79
           ELSEIF(br .le. 16.769/62.218) then
             nsta = 80
           ELSEIF(br .le. 16.969/62.218) then
             nsta = 82
           ELSEIF(br .le. 17.499/62.218) then
             nsta = 83
           ELSEIF(br .le. 17.829/62.218) then
             nsta = 85
           ELSEIF(br .le. 17.999/62.218) then
             nsta = 92
           ELSEIF(br .le. 18.179/62.218) then
             nsta = 97
           ELSEIF(br .le. 18.253/62.218) then
             nsta = 98
           ELSEIF(br .le. 18.613/62.218) then
             nsta = 102
           ELSEIF(br .le. 18.697/62.218) then
             nsta = 106
           ELSEIF(br .le. 18.867/62.218) then
             nsta = 107
           ELSEIF(br .le. 19.227/62.218) then
             nsta = 114
           ELSEIF(br .le. 19.331/62.218) then
             nsta = 115
           ELSEIF(br .le. 19.461/62.218) then
             nsta = 119
           ELSEIF(br .le. 20.121/62.218) then
             nsta = 121
           ELSEIF(br .le. 20.271/62.218) then
             nsta = 123
           ELSEIF(br .le. 20.337/62.218) then
             nsta = 126
           ELSEIF(br .le. 20.477/62.218) then
             nsta = 130
           ELSEIF(br .le. 21.207/62.218) then
             nsta = 131
           ELSEIF(br .le. 21.417/62.218) then
             nsta = 132
           ELSEIF(br .le. 21.827/62.218) then
             nsta = 134
           ELSEIF(br .le. 21.856/62.218) then
             nsta = 135
           ELSEIF(br .le. 22.056/62.218) then
             nsta = 139
           ELSEIF(br .le. 22.196/62.218) then
             nsta = 141
           ELSEIF(br .le. 22.356/62.218) then
             nsta = 145
           ELSEIF(br .le. 23.156/62.218) then
             nsta = 146
           ELSEIF(br .le. 23.356/62.218) then
             nsta = 149
           ELSEIF(br .le. 23.476/62.218) then
             nsta = 151
           ELSEIF(br .le. 24.276/62.218) then
             nsta = 154
           ELSEIF(br .le. 25.176/62.218) then
             nsta = 156
           ELSEIF(br .le. 25.466/62.218) then
             nsta = 158
           ELSEIF(br .le. 26.086/62.218) then
             nsta = 159
           ELSEIF(br .le. 26.566/62.218) then
             nsta = 163
           ELSEIF(br .le. 26.816/62.218) then
             nsta = 164
           ELSEIF(br .le. 26.926/62.218) then
             nsta = 165
           ELSEIF(br .le. 27.016/62.218) then
             nsta = 166
           ELSEIF(br .le. 27.396/62.218) then
             nsta = 169
           ELSEIF(br .le. 27.476/62.218) then
             nsta = 170
           ELSEIF(br .le. 27.536/62.218) then
             nsta = 173
           ELSEIF(br .le. 28.136/62.218) then
             nsta = 174
           ELSEIF(br .le. 28.766/62.218) then
             nsta = 176
           ELSEIF(br .le. 29.156/62.218) then
             nsta = 177
           ELSEIF(br .le. 29.396/62.218) then
             nsta = 180
           ELSEIF(br .le. 29.646/62.218) then
             nsta = 181
           ELSEIF(br .le. 29.986/62.218) then
             nsta = 183
           ELSEIF(br .le. 30.336/62.218) then
             nsta = 184
           ELSEIF(br .le. 30.586/62.218) then
             nsta = 185
           ELSEIF(br .le. 30.806/62.218) then
             nsta = 188
           ELSEIF(br .le. 31.286/62.218) then
             nsta = 192
           ELSEIF(br .le. 32.686/62.218) then
             nsta = 193
           ELSEIF(br .le. 32.976/62.218) then
             nsta = 195
           ELSEIF(br .le. 33.646/62.218) then
             nsta = 196
           ELSEIF(br .le. 33.856/62.218) then
             nsta = 200
           ELSEIF(br .le. 34.006/62.218) then
             nsta = 201
           ELSEIF(br .le. 34.566/62.218) then
             nsta = 203
           ELSEIF(br .le. 34.636/62.218) then
             nsta = 207
           ELSEIF(br .le. 34.724/62.218) then
             nsta = 210
           ELSEIF(br .le. 35.274/62.218) then
             nsta = 211
           ELSEIF(br .le. 36.294/62.218) then
             nsta = 213
           ELSEIF(br .le. 36.784/62.218) then
             nsta = 215
           ELSEIF(br .le. 36.924/62.218) then
             nsta = 216
           ELSEIF(br .le. 37.304/62.218) then
             nsta = 217
           ELSEIF(br .le. 37.644/62.218) then
             nsta = 219
           ELSEIF(br .le. 37.824/62.218) then
             nsta = 220
           ELSEIF(br .le. 38.024/62.218) then
             nsta = 222
           ELSEIF(br .le. 38.134/62.218) then
             nsta = 224
           ELSEIF(br .le. 38.514/62.218) then
             nsta = 227
           ELSEIF(br .le. 38.864/62.218) then
             nsta = 228
           ELSEIF(br .le. 39.704/62.218) then
             nsta = 230
           ELSEIF(br .le. 39.884/62.218) then
             nsta = 233
           ELSEIF(br .le. 40.014/62.218) then
             nsta = 234
           ELSEIF(br .le. 40.066/62.218) then
             nsta = 236
           ELSEIF(br .le. 40.946/62.218) then
             nsta = 238
           ELSEIF(br .le. 41.366/62.218) then
             nsta = 239
           ELSEIF(br .le. 41.616/62.218) then
             nsta = 241
           ELSEIF(br .le. 41.906/62.218) then
             nsta = 244
           ELSEIF(br .le. 42.146/62.218) then
             nsta = 246
           ELSEIF(br .le. 42.366/62.218) then
             nsta = 247
           ELSEIF(br .le. 42.516/62.218) then
             nsta = 249
           ELSEIF(br .le. 42.567/62.218) then
             nsta = 253
           ELSEIF(br .le. 43.177/62.218) then
             nsta = 254
           ELSEIF(br .le. 43.427/62.218) then
             nsta = 257
           ELSEIF(br .le. 43.504/62.218) then
             nsta = 260
           ELSEIF(br .le. 43.714/62.218) then
             nsta = 261
           ELSEIF(br .le. 43.994/62.218) then
             nsta = 263
           ELSEIF(br .le. 44.304/62.218) then
             nsta = 264
           ELSEIF(br .le. 44.584/62.218) then
             nsta = 265
           ELSEIF(br .le. 44.944/62.218) then
             nsta = 266
           ELSEIF(br .le. 45.214/62.218) then
             nsta = 268
           ELSEIF(br .le. 45.564/62.218) then
             nsta = 269
           ELSEIF(br .le. 45.694/62.218) then
             nsta = 270
           ELSEIF(br .le. 45.924/62.218) then
             nsta = 271
           ELSEIF(br .le. 46.274/62.218) then
             nsta = 273
           ELSEIF(br .le. 46.324/62.218) then
             nsta = 276
           ELSEIF(br .le. 46.37/62.218) then
             nsta = 278
           ELSEIF(br .le. 47.77/62.218) then
             nsta = 280
           ELSEIF(br .le. 47.98/62.218) then
             nsta = 281
           ELSEIF(br .le. 48.26/62.218) then
             nsta = 283
           ELSEIF(br .le. 48.5/62.218) then
             nsta = 285
           ELSEIF(br .le. 48.73/62.218) then
             nsta = 286
           ELSEIF(br .le. 48.94/62.218) then
             nsta = 287
           ELSEIF(br .le. 49.23/62.218) then
             nsta = 288
           ELSEIF(br .le. 49.97/62.218) then
             nsta = 290
           ELSEIF(br .le. 50.55/62.218) then
             nsta = 291
           ELSEIF(br .le. 51.85/62.218) then
             nsta = 292
           ELSEIF(br .le. 52.23/62.218) then
             nsta = 294
           ELSEIF(br .le. 52.41/62.218) then
             nsta = 295
           ELSEIF(br .le. 52.66/62.218) then
             nsta = 296
           ELSEIF(br .le. 53.14/62.218) then
             nsta = 299
           ELSEIF(br .le. 54.12/62.218) then
             nsta = 301
           ELSEIF(br .le. 54.61/62.218) then
             nsta = 302
           ELSEIF(br .le. 54.99/62.218) then
             nsta = 304
           ELSEIF(br .le. 55.21/62.218) then
             nsta = 305
           ELSEIF(br .le. 55.37/62.218) then
             nsta = 308
           ELSEIF(br .le. 55.468/62.218) then
             nsta = 309
           ELSEIF(br .le. 56.418/62.218) then
             nsta = 312
           ELSEIF(br .le. 56.628/62.218) then
             nsta = 313
           ELSEIF(br .le. 56.778/62.218) then
             nsta = 314
           ELSEIF(br .le. 57.958/62.218) then
             nsta = 317
           ELSEIF(br .le. 58.828/62.218) then
             nsta = 318
           ELSEIF(br .le. 60.528/62.218) then
             nsta = 320
           ELSEIF(br .le. 61.018/62.218) then
             nsta = 321
           ELSEIF(br .le. 61.198/62.218) then
             nsta = 322
           ELSE
             nsta = 324
           ENDIF
         Case(68169)
               IF(br .le. 0.0048/2.2196) then
             nsta = 0
           ELSEIF(br .le. 0.2648/2.2196) then
             nsta = 1
           ELSEIF(br .le. 0.3848/2.2196) then
             nsta = 15
           ELSEIF(br .le. 0.3896/2.2196) then
             nsta = 20
           ELSEIF(br .le. 0.4196/2.2196) then
             nsta = 30
           ELSEIF(br .le. 0.4306/2.2196) then
             nsta = 43
           ELSEIF(br .le. 0.8406/2.2196) then
             nsta = 45
           ELSEIF(br .le. 0.8526/2.2196) then
             nsta = 47
           ELSEIF(br .le. 0.8596/2.2196) then
             nsta = 50
           ELSEIF(br .le. 1.0996/2.2196) then
             nsta = 65
           ELSEIF(br .le. 1.1176/2.2196) then
             nsta = 71
           ELSEIF(br .le. 1.1236/2.2196) then
             nsta = 72
           ELSEIF(br .le. 1.1356/2.2196) then
             nsta = 73
           ELSEIF(br .le. 1.1496/2.2196) then
             nsta = 75
           ELSEIF(br .le. 1.1636/2.2196) then
             nsta = 78
           ELSEIF(br .le. 1.1806/2.2196) then
             nsta = 85
           ELSEIF(br .le. 1.1826/2.2196) then
             nsta = 87
           ELSEIF(br .le. 1.187/2.2196) then
             nsta = 89
           ELSEIF(br .le. 1.337/2.2196) then
             nsta = 91
           ELSEIF(br .le. 1.353/2.2196) then
             nsta = 98
           ELSEIF(br .le. 1.36/2.2196) then
             nsta = 101
           ELSEIF(br .le. 1.362/2.2196) then
             nsta = 102
           ELSEIF(br .le. 1.3656/2.2196) then
             nsta = 103
           ELSEIF(br .le. 1.3736/2.2196) then
             nsta = 104
           ELSEIF(br .le. 1.3856/2.2196) then
             nsta = 105
           ELSEIF(br .le. 1.3976/2.2196) then
             nsta = 106
           ELSEIF(br .le. 1.4086/2.2196) then
             nsta = 108
           ELSEIF(br .le. 1.5786/2.2196) then
             nsta = 110
           ELSEIF(br .le. 1.6686/2.2196) then
             nsta = 114
           ELSEIF(br .le. 1.6742/2.2196) then
             nsta = 115
           ELSEIF(br .le. 1.6778/2.2196) then
             nsta = 116
           ELSEIF(br .le. 1.6868/2.2196) then
             nsta = 119
           ELSEIF(br .le. 1.7198/2.2196) then
             nsta = 120
           ELSEIF(br .le. 1.9298/2.2196) then
             nsta = 123
           ELSEIF(br .le. 1.9338/2.2196) then
             nsta = 124
           ELSEIF(br .le. 1.9468/2.2196) then
             nsta = 126
           ELSEIF(br .le. 1.9538/2.2196) then
             nsta = 129
           ELSEIF(br .le. 2.0138/2.2196) then
             nsta = 130
           ELSEIF(br .le. 2.0166/2.2196) then
             nsta = 131
           ELSEIF(br .le. 2.0186/2.2196) then
             nsta = 133
           ELSEIF(br .le. 2.0866/2.2196) then
             nsta = 134
           ELSEIF(br .le. 2.1766/2.2196) then
             nsta = 135
           ELSEIF(br .le. 2.1986/2.2196) then
             nsta = 137
           ELSEIF(br .le. 2.2056/2.2196) then
             nsta = 138
           ELSE
             nsta = 139
           ENDIF
         Case(68171)
               IF(br .le. 0.0015/0.5942) then
             nsta = 3
           ELSEIF(br .le. 0.0105/0.5942) then
             nsta = 4
           ELSEIF(br .le. 0.0124/0.5942) then
             nsta = 17
           ELSEIF(br .le. 0.0284/0.5942) then
             nsta = 24
           ELSEIF(br .le. 0.0984/0.5942) then
             nsta = 31
           ELSEIF(br .le. 0.1044/0.5942) then
             nsta = 32
           ELSEIF(br .le. 0.1072/0.5942) then
             nsta = 34
           ELSEIF(br .le. 0.11/0.5942) then
             nsta = 35
           ELSEIF(br .le. 0.14/0.5942) then
             nsta = 39
           ELSEIF(br .le. 0.37/0.5942) then
             nsta = 41
           ELSEIF(br .le. 0.392/0.5942) then
             nsta = 46
           ELSEIF(br .le. 0.3942/0.5942) then
             nsta = 49
           ELSEIF(br .le. 0.3992/0.5942) then
             nsta = 51
           ELSEIF(br .le. 0.4052/0.5942) then
             nsta = 53
           ELSEIF(br .le. 0.4152/0.5942) then
             nsta = 57
           ELSEIF(br .le. 0.4282/0.5942) then
             nsta = 58
           ELSEIF(br .le. 0.4372/0.5942) then
             nsta = 59
           ELSEIF(br .le. 0.4552/0.5942) then
             nsta = 60
           ELSEIF(br .le. 0.4642/0.5942) then
             nsta = 61
           ELSEIF(br .le. 0.5442/0.5942) then
             nsta = 62
           ELSE
             nsta = 63
           ENDIF
         Case(69170)
               IF(br .le. 0.03/16.152) then
             nsta = 0
           ELSEIF(br .le. 0.68/16.152) then
             nsta = 1
           ELSEIF(br .le. 1.15/16.152) then
             nsta = 3
           ELSEIF(br .le. 2.63/16.152) then
             nsta = 6
           ELSEIF(br .le. 2.715/16.152) then
             nsta = 7
           ELSEIF(br .le. 3.135/16.152) then
             nsta = 8
           ELSEIF(br .le. 4.125/16.152) then
             nsta = 29
           ELSEIF(br .le. 5.155/16.152) then
             nsta = 38
           ELSEIF(br .le. 6.665/16.152) then
             nsta = 39
           ELSEIF(br .le. 7.155/16.152) then
             nsta = 43
           ELSEIF(br .le. 7.505/16.152) then
             nsta = 45
           ELSEIF(br .le. 7.915/16.152) then
             nsta = 51
           ELSEIF(br .le. 8.062/16.152) then
             nsta = 60
           ELSEIF(br .le. 9.482/16.152) then
             nsta = 72
           ELSEIF(br .le. 10.652/16.152) then
             nsta = 73
           ELSEIF(br .le. 10.912/16.152) then
             nsta = 75
           ELSEIF(br .le. 11.016/16.152) then
             nsta = 79
           ELSEIF(br .le. 11.176/16.152) then
             nsta = 99
           ELSEIF(br .le. 11.324/16.152) then
             nsta = 105
           ELSEIF(br .le. 11.474/16.152) then
             nsta = 107
           ELSEIF(br .le. 11.704/16.152) then
             nsta = 108
           ELSEIF(br .le. 11.944/16.152) then
             nsta = 109
           ELSEIF(br .le. 12.334/16.152) then
             nsta = 110
           ELSEIF(br .le. 12.477/16.152) then
             nsta = 112
           ELSEIF(br .le. 12.657/16.152) then
             nsta = 114
           ELSEIF(br .le. 12.847/16.152) then
             nsta = 120
           ELSEIF(br .le. 13.027/16.152) then
             nsta = 123
           ELSEIF(br .le. 13.119/16.152) then
             nsta = 135
           ELSEIF(br .le. 13.589/16.152) then
             nsta = 139
           ELSEIF(br .le. 13.899/16.152) then
             nsta = 142
           ELSEIF(br .le. 14.179/16.152) then
             nsta = 146
           ELSEIF(br .le. 14.569/16.152) then
             nsta = 152
           ELSEIF(br .le. 14.672/16.152) then
             nsta = 155
           ELSEIF(br .le. 14.832/16.152) then
             nsta = 162
           ELSEIF(br .le. 15.092/16.152) then
             nsta = 167
           ELSEIF(br .le. 15.252/16.152) then
             nsta = 177
           ELSEIF(br .le. 15.832/16.152) then
             nsta = 181
           ELSE
             nsta = 187
           ENDIF
         Case(70169)
               IF(br .le. 0.0133/0.19332) then
             nsta = 1
           ELSEIF(br .le. 0.0713/0.19332) then
             nsta = 3
           ELSEIF(br .le. 0.08/0.19332) then
             nsta = 21
           ELSEIF(br .le. 0.08116/0.19332) then
             nsta = 24
           ELSEIF(br .le. 0.09146/0.19332) then
             nsta = 32
           ELSEIF(br .le. 0.09196/0.19332) then
             nsta = 36
           ELSEIF(br .le. 0.09282/0.19332) then
             nsta = 47
           ELSEIF(br .le. 0.09334/0.19332) then
             nsta = 57
           ELSEIF(br .le. 0.09397/0.19332) then
             nsta = 62
           ELSEIF(br .le. 0.09425/0.19332) then
             nsta = 65
           ELSEIF(br .le. 0.09453/0.19332) then
             nsta = 66
           ELSEIF(br .le. 0.09903/0.19332) then
             nsta = 69
           ELSEIF(br .le. 0.10203/0.19332) then
             nsta = 77
           ELSEIF(br .le. 0.10553/0.19332) then
             nsta = 78
           ELSEIF(br .le. 0.11183/0.19332) then
             nsta = 80
           ELSEIF(br .le. 0.11324/0.19332) then
             nsta = 89
           ELSEIF(br .le. 0.11376/0.19332) then
             nsta = 90
           ELSEIF(br .le. 0.11489/0.19332) then
             nsta = 91
           ELSEIF(br .le. 0.11672/0.19332) then
             nsta = 92
           ELSEIF(br .le. 0.11775/0.19332) then
             nsta = 94
           ELSEIF(br .le. 0.11844/0.19332) then
             nsta = 97
           ELSEIF(br .le. 0.13254/0.19332) then
             nsta = 100
           ELSEIF(br .le. 0.13594/0.19332) then
             nsta = 110
           ELSEIF(br .le. 0.13704/0.19332) then
             nsta = 114
           ELSEIF(br .le. 0.13817/0.19332) then
             nsta = 117
           ELSEIF(br .le. 0.1388/0.19332) then
             nsta = 118
           ELSEIF(br .le. 0.1419/0.19332) then
             nsta = 122
           ELSEIF(br .le. 0.145/0.19332) then
             nsta = 123
           ELSEIF(br .le. 0.14683/0.19332) then
             nsta = 124
           ELSEIF(br .le. 0.14838/0.19332) then
             nsta = 125
           ELSEIF(br .le. 0.15008/0.19332) then
             nsta = 127
           ELSEIF(br .le. 0.1506/0.19332) then
             nsta = 128
           ELSEIF(br .le. 0.15201/0.19332) then
             nsta = 130
           ELSEIF(br .le. 0.15242/0.19332) then
             nsta = 132
           ELSEIF(br .le. 0.15332/0.19332) then
             nsta = 133
           ELSEIF(br .le. 0.15782/0.19332) then
             nsta = 134
           ELSEIF(br .le. 0.16062/0.19332) then
             nsta = 138
           ELSEIF(br .le. 0.16372/0.19332) then
             nsta = 140
           ELSEIF(br .le. 0.16424/0.19332) then
             nsta = 141
           ELSEIF(br .le. 0.16579/0.19332) then
             nsta = 144
           ELSEIF(br .le. 0.16637/0.19332) then
             nsta = 145
           ELSEIF(br .le. 0.16712/0.19332) then
             nsta = 148
           ELSEIF(br .le. 0.16792/0.19332) then
             nsta = 149
           ELSEIF(br .le. 0.16962/0.19332) then
             nsta = 150
           ELSEIF(br .le. 0.17372/0.19332) then
             nsta = 151
           ELSEIF(br .le. 0.17419/0.19332) then
             nsta = 153
           ELSEIF(br .le. 0.1754/0.19332) then
             nsta = 156
           ELSEIF(br .le. 0.17737/0.19332) then
             nsta = 157
           ELSEIF(br .le. 0.17877/0.19332) then
             nsta = 158
           ELSEIF(br .le. 0.18187/0.19332) then
             nsta = 159
           ELSEIF(br .le. 0.18487/0.19332) then
             nsta = 160
           ELSEIF(br .le. 0.18545/0.19332) then
             nsta = 161
           ELSEIF(br .le. 0.18614/0.19332) then
             nsta = 162
           ELSEIF(br .le. 0.18811/0.19332) then
             nsta = 163
           ELSEIF(br .le. 0.18981/0.19332) then
             nsta = 164
           ELSEIF(br .le. 0.19102/0.19332) then
             nsta = 165
           ELSE
             nsta = 166
           ENDIF
         Case(70172)
               IF(br .le. 0.0088/0.51833) then
             nsta = 0
           ELSEIF(br .le. 0.01029/0.51833) then
             nsta = 1
           ELSEIF(br .le. 0.01499/0.51833) then
             nsta = 5
           ELSEIF(br .le. 0.01779/0.51833) then
             nsta = 6
           ELSEIF(br .le. 0.02029/0.51833) then
             nsta = 7
           ELSEIF(br .le. 0.02269/0.51833) then
             nsta = 9
           ELSEIF(br .le. 0.02989/0.51833) then
             nsta = 17
           ELSEIF(br .le. 0.03179/0.51833) then
             nsta = 18
           ELSEIF(br .le. 0.03729/0.51833) then
             nsta = 19
           ELSEIF(br .le. 0.03776/0.51833) then
             nsta = 28
           ELSEIF(br .le. 0.03804/0.51833) then
             nsta = 29
           ELSEIF(br .le. 0.04294/0.51833) then
             nsta = 45
           ELSEIF(br .le. 0.05014/0.51833) then
             nsta = 53
           ELSEIF(br .le. 0.05504/0.51833) then
             nsta = 58
           ELSEIF(br .le. 0.05539/0.51833) then
             nsta = 61
           ELSEIF(br .le. 0.05645/0.51833) then
             nsta = 64
           ELSEIF(br .le. 0.07125/0.51833) then
             nsta = 68
           ELSEIF(br .le. 0.07182/0.51833) then
             nsta = 71
           ELSEIF(br .le. 0.07285/0.51833) then
             nsta = 75
           ELSEIF(br .le. 0.07457/0.51833) then
             nsta = 78
           ELSEIF(br .le. 0.09177/0.51833) then
             nsta = 93
           ELSEIF(br .le. 0.09457/0.51833) then
             nsta = 101
           ELSEIF(br .le. 0.09677/0.51833) then
             nsta = 111
           ELSEIF(br .le. 0.09828/0.51833) then
             nsta = 112
           ELSEIF(br .le. 0.11828/0.51833) then
             nsta = 113
           ELSEIF(br .le. 0.12038/0.51833) then
             nsta = 116
           ELSEIF(br .le. 0.12478/0.51833) then
             nsta = 123
           ELSEIF(br .le. 0.12512/0.51833) then
             nsta = 124
           ELSEIF(br .le. 0.12682/0.51833) then
             nsta = 131
           ELSEIF(br .le. 0.20982/0.51833) then
             nsta = 133
           ELSEIF(br .le. 0.21292/0.51833) then
             nsta = 136
           ELSEIF(br .le. 0.21389/0.51833) then
             nsta = 139
           ELSEIF(br .le. 0.21442/0.51833) then
             nsta = 140
           ELSEIF(br .le. 0.21485/0.51833) then
             nsta = 141
           ELSEIF(br .le. 0.21527/0.51833) then
             nsta = 143
           ELSEIF(br .le. 0.21605/0.51833) then
             nsta = 145
           ELSEIF(br .le. 0.21751/0.51833) then
             nsta = 148
           ELSEIF(br .le. 0.21794/0.51833) then
             nsta = 149
           ELSEIF(br .le. 0.21828/0.51833) then
             nsta = 150
           ELSEIF(br .le. 0.21862/0.51833) then
             nsta = 151
           ELSEIF(br .le. 0.22062/0.51833) then
             nsta = 154
           ELSEIF(br .le. 0.22492/0.51833) then
             nsta = 157
           ELSEIF(br .le. 0.22612/0.51833) then
             nsta = 162
           ELSEIF(br .le. 0.22695/0.51833) then
             nsta = 163
           ELSEIF(br .le. 0.22856/0.51833) then
             nsta = 166
           ELSEIF(br .le. 0.23046/0.51833) then
             nsta = 169
           ELSEIF(br .le. 0.23466/0.51833) then
             nsta = 173
           ELSEIF(br .le. 0.235/0.51833) then
             nsta = 174
           ELSEIF(br .le. 0.2387/0.51833) then
             nsta = 175
           ELSEIF(br .le. 0.2402/0.51833) then
             nsta = 176
           ELSEIF(br .le. 0.2424/0.51833) then
             nsta = 179
           ELSEIF(br .le. 0.2436/0.51833) then
             nsta = 181
           ELSEIF(br .le. 0.245/0.51833) then
             nsta = 182
           ELSEIF(br .le. 0.2477/0.51833) then
             nsta = 184
           ELSEIF(br .le. 0.2505/0.51833) then
             nsta = 186
           ELSEIF(br .le. 0.2519/0.51833) then
             nsta = 188
           ELSEIF(br .le. 0.2638/0.51833) then
             nsta = 190
           ELSEIF(br .le. 0.26458/0.51833) then
             nsta = 192
           ELSEIF(br .le. 0.26568/0.51833) then
             nsta = 194
           ELSEIF(br .le. 0.26768/0.51833) then
             nsta = 195
           ELSEIF(br .le. 0.27028/0.51833) then
             nsta = 196
           ELSEIF(br .le. 0.27188/0.51833) then
             nsta = 198
           ELSEIF(br .le. 0.27418/0.51833) then
             nsta = 200
           ELSEIF(br .le. 0.27678/0.51833) then
             nsta = 201
           ELSEIF(br .le. 0.27848/0.51833) then
             nsta = 205
           ELSEIF(br .le. 0.28188/0.51833) then
             nsta = 207
           ELSEIF(br .le. 0.28398/0.51833) then
             nsta = 212
           ELSEIF(br .le. 0.28708/0.51833) then
             nsta = 216
           ELSEIF(br .le. 0.29168/0.51833) then
             nsta = 219
           ELSEIF(br .le. 0.29908/0.51833) then
             nsta = 220
           ELSEIF(br .le. 0.30268/0.51833) then
             nsta = 222
           ELSEIF(br .le. 0.30356/0.51833) then
             nsta = 228
           ELSEIF(br .le. 0.30455/0.51833) then
             nsta = 230
           ELSEIF(br .le. 0.30735/0.51833) then
             nsta = 235
           ELSEIF(br .le. 0.31675/0.51833) then
             nsta = 237
           ELSEIF(br .le. 0.31885/0.51833) then
             nsta = 238
           ELSEIF(br .le. 0.32075/0.51833) then
             nsta = 239
           ELSEIF(br .le. 0.32295/0.51833) then
             nsta = 240
           ELSEIF(br .le. 0.32383/0.51833) then
             nsta = 241
           ELSEIF(br .le. 0.32603/0.51833) then
             nsta = 244
           ELSEIF(br .le. 0.33103/0.51833) then
             nsta = 245
           ELSEIF(br .le. 0.33383/0.51833) then
             nsta = 246
           ELSEIF(br .le. 0.33533/0.51833) then
             nsta = 247
           ELSEIF(br .le. 0.34133/0.51833) then
             nsta = 248
           ELSEIF(br .le. 0.34653/0.51833) then
             nsta = 249
           ELSEIF(br .le. 0.34793/0.51833) then
             nsta = 252
           ELSEIF(br .le. 0.34983/0.51833) then
             nsta = 253
           ELSEIF(br .le. 0.35263/0.51833) then
             nsta = 255
           ELSEIF(br .le. 0.35553/0.51833) then
             nsta = 257
           ELSEIF(br .le. 0.36383/0.51833) then
             nsta = 258
           ELSEIF(br .le. 0.36553/0.51833) then
             nsta = 259
           ELSEIF(br .le. 0.36823/0.51833) then
             nsta = 260
           ELSEIF(br .le. 0.37103/0.51833) then
             nsta = 262
           ELSEIF(br .le. 0.37553/0.51833) then
             nsta = 263
           ELSEIF(br .le. 0.37953/0.51833) then
             nsta = 264
           ELSEIF(br .le. 0.38143/0.51833) then
             nsta = 268
           ELSEIF(br .le. 0.38483/0.51833) then
             nsta = 269
           ELSEIF(br .le. 0.38783/0.51833) then
             nsta = 271
           ELSEIF(br .le. 0.39093/0.51833) then
             nsta = 272
           ELSEIF(br .le. 0.39213/0.51833) then
             nsta = 273
           ELSEIF(br .le. 0.39623/0.51833) then
             nsta = 274
           ELSEIF(br .le. 0.40043/0.51833) then
             nsta = 275
           ELSEIF(br .le. 0.40473/0.51833) then
             nsta = 276
           ELSEIF(br .le. 0.41413/0.51833) then
             nsta = 277
           ELSEIF(br .le. 0.42243/0.51833) then
             nsta = 278
           ELSEIF(br .le. 0.42453/0.51833) then
             nsta = 279
           ELSEIF(br .le. 0.42983/0.51833) then
             nsta = 280
           ELSEIF(br .le. 0.43143/0.51833) then
             nsta = 282
           ELSEIF(br .le. 0.43823/0.51833) then
             nsta = 283
           ELSEIF(br .le. 0.43943/0.51833) then
             nsta = 284
           ELSEIF(br .le. 0.44213/0.51833) then
             nsta = 286
           ELSEIF(br .le. 0.44573/0.51833) then
             nsta = 288
           ELSEIF(br .le. 0.44743/0.51833) then
             nsta = 291
           ELSEIF(br .le. 0.45103/0.51833) then
             nsta = 292
           ELSEIF(br .le. 0.45283/0.51833) then
             nsta = 293
           ELSEIF(br .le. 0.45653/0.51833) then
             nsta = 294
           ELSEIF(br .le. 0.46023/0.51833) then
             nsta = 295
           ELSEIF(br .le. 0.46853/0.51833) then
             nsta = 296
           ELSEIF(br .le. 0.47533/0.51833) then
             nsta = 297
           ELSEIF(br .le. 0.47823/0.51833) then
             nsta = 298
           ELSEIF(br .le. 0.48203/0.51833) then
             nsta = 299
           ELSEIF(br .le. 0.48483/0.51833) then
             nsta = 300
           ELSEIF(br .le. 0.48623/0.51833) then
             nsta = 301
           ELSEIF(br .le. 0.48863/0.51833) then
             nsta = 302
           ELSEIF(br .le. 0.49343/0.51833) then
             nsta = 303
           ELSEIF(br .le. 0.50003/0.51833) then
             nsta = 304
           ELSEIF(br .le. 0.51153/0.51833) then
             nsta = 305
           ELSE
             nsta = 306
           ENDIF
         Case(70173)
               IF(br .le. 0.012/0.03951) then
             nsta = 5
           ELSEIF(br .le. 0.0136/0.03951) then
             nsta = 8
           ELSEIF(br .le. 0.0245/0.03951) then
             nsta = 23
           ELSEIF(br .le. 0.0338/0.03951) then
             nsta = 25
           ELSEIF(br .le. 0.0362/0.03951) then
             nsta = 33
           ELSEIF(br .le. 0.03721/0.03951) then
             nsta = 37
           ELSE
             nsta = 44
           ENDIF
         Case(70174)
               IF(br .le. 0.014/0.157) then
             nsta = 91
           ELSEIF(br .le. 0.027/0.157) then
             nsta = 96
           ELSEIF(br .le. 0.054/0.157) then
             nsta = 114
           ELSE
             nsta = 187
           ENDIF
         Case(70175)
               IF(br .le. 0.02/4.2879) then
             nsta = 7
           ELSEIF(br .le. 1.42/4.2879) then
             nsta = 10
           ELSEIF(br .le. 1.6/4.2879) then
             nsta = 18
           ELSEIF(br .le. 1.6046/4.2879) then
             nsta = 22
           ELSEIF(br .le. 1.8546/4.2879) then
             nsta = 26
           ELSEIF(br .le. 1.8623/4.2879) then
             nsta = 36
           ELSEIF(br .le. 1.9023/4.2879) then
             nsta = 50
           ELSEIF(br .le. 1.9283/4.2879) then
             nsta = 51
           ELSEIF(br .le. 1.9345/4.2879) then
             nsta = 55
           ELSEIF(br .le. 1.9437/4.2879) then
             nsta = 56
           ELSEIF(br .le. 1.9499/4.2879) then
             nsta = 60
           ELSEIF(br .le. 2.0079/4.2879) then
             nsta = 65
           ELSEIF(br .le. 2.0959/4.2879) then
             nsta = 68
           ELSEIF(br .le. 2.1099/4.2879) then
             nsta = 70
           ELSEIF(br .le. 2.1329/4.2879) then
             nsta = 72
           ELSEIF(br .le. 2.1589/4.2879) then
             nsta = 75
           ELSEIF(br .le. 2.1789/4.2879) then
             nsta = 89
           ELSEIF(br .le. 2.5589/4.2879) then
             nsta = 97
           ELSEIF(br .le. 3.2789/4.2879) then
             nsta = 102
           ELSEIF(br .le. 3.3639/4.2879) then
             nsta = 105
           ELSEIF(br .le. 3.4379/4.2879) then
             nsta = 106
           ELSEIF(br .le. 3.4639/4.2879) then
             nsta = 108
           ELSEIF(br .le. 3.5209/4.2879) then
             nsta = 111
           ELSEIF(br .le. 3.5609/4.2879) then
             nsta = 112
           ELSEIF(br .le. 3.5719/4.2879) then
             nsta = 114
           ELSEIF(br .le. 3.6149/4.2879) then
             nsta = 115
           ELSEIF(br .le. 3.8449/4.2879) then
             nsta = 117
           ELSEIF(br .le. 3.8879/4.2879) then
             nsta = 124
           ELSE
             nsta = 126
           ENDIF
         Case(70177)
               IF(br .le. 0.0005/0.03027) then
             nsta = 6
           ELSEIF(br .le. 0.0007/0.03027) then
             nsta = 8
           ELSEIF(br .le. 0.00087/0.03027) then
             nsta = 17
           ELSEIF(br .le. 0.00095/0.03027) then
             nsta = 21
           ELSEIF(br .le. 0.00133/0.03027) then
             nsta = 23
           ELSEIF(br .le. 0.00353/0.03027) then
             nsta = 26
           ELSEIF(br .le. 0.00389/0.03027) then
             nsta = 29
           ELSEIF(br .le. 0.00423/0.03027) then
             nsta = 37
           ELSEIF(br .le. 0.00623/0.03027) then
             nsta = 39
           ELSEIF(br .le. 0.01123/0.03027) then
             nsta = 43
           ELSEIF(br .le. 0.0114/0.03027) then
             nsta = 45
           ELSEIF(br .le. 0.0122/0.03027) then
             nsta = 47
           ELSEIF(br .le. 0.01239/0.03027) then
             nsta = 57
           ELSEIF(br .le. 0.01699/0.03027) then
             nsta = 61
           ELSEIF(br .le. 0.01829/0.03027) then
             nsta = 65
           ELSEIF(br .le. 0.01863/0.03027) then
             nsta = 69
           ELSEIF(br .le. 0.0189/0.03027) then
             nsta = 75
           ELSEIF(br .le. 0.01902/0.03027) then
             nsta = 82
           ELSEIF(br .le. 0.01929/0.03027) then
             nsta = 83
           ELSEIF(br .le. 0.01968/0.03027) then
             nsta = 87
           ELSEIF(br .le. 0.02068/0.03027) then
             nsta = 92
           ELSEIF(br .le. 0.02128/0.03027) then
             nsta = 94
           ELSEIF(br .le. 0.02488/0.03027) then
             nsta = 96
           ELSEIF(br .le. 0.02568/0.03027) then
             nsta = 100
           ELSEIF(br .le. 0.02614/0.03027) then
             nsta = 109
           ELSEIF(br .le. 0.02659/0.03027) then
             nsta = 111
           ELSEIF(br .le. 0.02809/0.03027) then
             nsta = 118
           ELSEIF(br .le. 0.02879/0.03027) then
             nsta = 138
           ELSEIF(br .le. 0.02939/0.03027) then
             nsta = 145
           ELSEIF(br .le. 0.02984/0.03027) then
             nsta = 149
           ELSE
             nsta = 156
           ENDIF
         Case(71176)
               IF(br .le. 0.0107/1.0259) then
             nsta = 11
           ELSEIF(br .le. 0.0227/1.0259) then
             nsta = 15
           ELSEIF(br .le. 0.0907/1.0259) then
             nsta = 18
           ELSEIF(br .le. 0.1027/1.0259) then
             nsta = 19
           ELSEIF(br .le. 0.1103/1.0259) then
             nsta = 21
           ELSEIF(br .le. 0.1194/1.0259) then
             nsta = 22
           ELSEIF(br .le. 0.1344/1.0259) then
             nsta = 26
           ELSEIF(br .le. 0.1451/1.0259) then
             nsta = 29
           ELSEIF(br .le. 0.1801/1.0259) then
             nsta = 32
           ELSEIF(br .le. 0.2621/1.0259) then
             nsta = 44
           ELSEIF(br .le. 0.2771/1.0259) then
             nsta = 47
           ELSEIF(br .le. 0.3121/1.0259) then
             nsta = 49
           ELSEIF(br .le. 0.3411/1.0259) then
             nsta = 56
           ELSEIF(br .le. 0.4321/1.0259) then
             nsta = 66
           ELSEIF(br .le. 0.4561/1.0259) then
             nsta = 71
           ELSEIF(br .le. 0.4622/1.0259) then
             nsta = 72
           ELSEIF(br .le. 0.5052/1.0259) then
             nsta = 76
           ELSEIF(br .le. 0.5352/1.0259) then
             nsta = 80
           ELSEIF(br .le. 0.5562/1.0259) then
             nsta = 82
           ELSEIF(br .le. 0.5802/1.0259) then
             nsta = 87
           ELSEIF(br .le. 0.7702/1.0259) then
             nsta = 88
           ELSEIF(br .le. 0.9302/1.0259) then
             nsta = 88
           ELSEIF(br .le. 0.9592/1.0259) then
             nsta = 95
           ELSEIF(br .le. 0.9699/1.0259) then
             nsta = 102
           ELSEIF(br .le. 0.9959/1.0259) then
             nsta = 105
           ELSE
             nsta = 109
           ENDIF
         Case(71177)
               IF(br .le. 0.38/5.755) then
             nsta = 3
           ELSEIF(br .le. 0.388/5.755) then
             nsta = 4
           ELSEIF(br .le. 0.409/5.755) then
             nsta = 6
           ELSEIF(br .le. 0.529/5.755) then
             nsta = 11
           ELSEIF(br .le. 0.607/5.755) then
             nsta = 20
           ELSEIF(br .le. 0.619/5.755) then
             nsta = 28
           ELSEIF(br .le. 0.664/5.755) then
             nsta = 30
           ELSEIF(br .le. 0.685/5.755) then
             nsta = 31
           ELSEIF(br .le. 0.693/5.755) then
             nsta = 45
           ELSEIF(br .le. 0.78/5.755) then
             nsta = 46
           ELSEIF(br .le. 0.792/5.755) then
             nsta = 47
           ELSEIF(br .le. 0.809/5.755) then
             nsta = 48
           ELSEIF(br .le. 0.85/5.755) then
             nsta = 51
           ELSEIF(br .le. 0.879/5.755) then
             nsta = 55
           ELSEIF(br .le. 1.063/5.755) then
             nsta = 57
           ELSEIF(br .le. 1.104/5.755) then
             nsta = 58
           ELSEIF(br .le. 1.137/5.755) then
             nsta = 61
           ELSEIF(br .le. 1.216/5.755) then
             nsta = 63
           ELSEIF(br .le. 1.446/5.755) then
             nsta = 65
           ELSEIF(br .le. 1.52/5.755) then
             nsta = 68
           ELSEIF(br .le. 1.578/5.755) then
             nsta = 72
           ELSEIF(br .le. 1.607/5.755) then
             nsta = 75
           ELSEIF(br .le. 1.665/5.755) then
             nsta = 76
           ELSEIF(br .le. 1.992/5.755) then
             nsta = 77
           ELSEIF(br .le. 2.058/5.755) then
             nsta = 79
           ELSEIF(br .le. 2.083/5.755) then
             nsta = 80
           ELSEIF(br .le. 2.468/5.755) then
             nsta = 82
           ELSEIF(br .le. 2.493/5.755) then
             nsta = 88
           ELSEIF(br .le. 2.501/5.755) then
             nsta = 92
           ELSEIF(br .le. 2.518/5.755) then
             nsta = 94
           ELSEIF(br .le. 2.736/5.755) then
             nsta = 96
           ELSEIF(br .le. 2.748/5.755) then
             nsta = 96
           ELSEIF(br .le. 2.773/5.755) then
             nsta = 97
           ELSEIF(br .le. 2.835/5.755) then
             nsta = 101
           ELSEIF(br .le. 2.847/5.755) then
             nsta = 103
           ELSEIF(br .le. 2.868/5.755) then
             nsta = 105
           ELSEIF(br .le. 2.984/5.755) then
             nsta = 109
           ELSEIF(br .le. 3.025/5.755) then
             nsta = 111
           ELSEIF(br .le. 3.107/5.755) then
             nsta = 113
           ELSEIF(br .le. 3.367/5.755) then
             nsta = 115
           ELSEIF(br .le. 3.379/5.755) then
             nsta = 117
           ELSEIF(br .le. 3.524/5.755) then
             nsta = 120
           ELSEIF(br .le. 3.691/5.755) then
             nsta = 121
           ELSEIF(br .le. 3.703/5.755) then
             nsta = 127
           ELSEIF(br .le. 3.794/5.755) then
             nsta = 130
           ELSEIF(br .le. 3.815/5.755) then
             nsta = 131
           ELSEIF(br .le. 3.84/5.755) then
             nsta = 134
           ELSEIF(br .le. 3.852/5.755) then
             nsta = 137
           ELSEIF(br .le. 3.873/5.755) then
             nsta = 139
           ELSEIF(br .le. 3.947/5.755) then
             nsta = 140
           ELSEIF(br .le. 4.025/5.755) then
             nsta = 141
           ELSEIF(br .le. 4.162/5.755) then
             nsta = 144
           ELSEIF(br .le. 4.174/5.755) then
             nsta = 145
           ELSEIF(br .le. 4.195/5.755) then
             nsta = 146
           ELSEIF(br .le. 4.294/5.755) then
             nsta = 147
           ELSEIF(br .le. 4.41/5.755) then
             nsta = 148
           ELSEIF(br .le. 4.447/5.755) then
             nsta = 150
           ELSEIF(br .le. 4.468/5.755) then
             nsta = 154
           ELSEIF(br .le. 4.48/5.755) then
             nsta = 156
           ELSEIF(br .le. 4.505/5.755) then
             nsta = 157
           ELSEIF(br .le. 4.63/5.755) then
             nsta = 160
           ELSEIF(br .le. 4.647/5.755) then
             nsta = 161
           ELSEIF(br .le. 4.684/5.755) then
             nsta = 164
           ELSEIF(br .le. 4.717/5.755) then
             nsta = 165
           ELSEIF(br .le. 4.738/5.755) then
             nsta = 166
           ELSEIF(br .le. 4.759/5.755) then
             nsta = 167
           ELSEIF(br .le. 4.78/5.755) then
             nsta = 168
           ELSEIF(br .le. 4.813/5.755) then
             nsta = 169
           ELSEIF(br .le. 4.838/5.755) then
             nsta = 170
           ELSEIF(br .le. 5.014/5.755) then
             nsta = 171
           ELSEIF(br .le. 5.229/5.755) then
             nsta = 173
           ELSEIF(br .le. 5.344/5.755) then
             nsta = 177
           ELSEIF(br .le. 5.389/5.755) then
             nsta = 178
           ELSEIF(br .le. 5.414/5.755) then
             nsta = 180
           ELSEIF(br .le. 5.505/5.755) then
             nsta = 181
           ELSE
             nsta = 189
           ENDIF
         Case(72175)
               IF(br .le. 0.0099/0.1151) then
             nsta = 2
           ELSEIF(br .le. 0.0259/0.1151) then
             nsta = 4
           ELSEIF(br .le. 0.0849/0.1151) then
             nsta = 27
           ELSEIF(br .le. 0.0901/0.1151) then
             nsta = 29
           ELSEIF(br .le. 0.0934/0.1151) then
             nsta = 32
           ELSEIF(br .le. 0.0976/0.1151) then
             nsta = 50
           ELSEIF(br .le. 0.1053/0.1151) then
             nsta = 52
           ELSE
             nsta = 63
           ENDIF
         Case(72178)
               IF(br .le. 0.122/4.696) then
             nsta = 6
           ELSEIF(br .le. 0.442/4.696) then
             nsta = 9
           ELSEIF(br .le. 0.592/4.696) then
             nsta = 12
           ELSEIF(br .le. 0.659/4.696) then
             nsta = 15
           ELSEIF(br .le. 0.739/4.696) then
             nsta = 16
           ELSEIF(br .le. 1.659/4.696) then
             nsta = 26
           ELSEIF(br .le. 1.771/4.696) then
             nsta = 27
           ELSEIF(br .le. 1.946/4.696) then
             nsta = 30
           ELSEIF(br .le. 2.116/4.696) then
             nsta = 35
           ELSEIF(br .le. 2.466/4.696) then
             nsta = 52
           ELSEIF(br .le. 2.606/4.696) then
             nsta = 69
           ELSEIF(br .le. 2.686/4.696) then
             nsta = 70
           ELSEIF(br .le. 2.836/4.696) then
             nsta = 72
           ELSEIF(br .le. 3.246/4.696) then
             nsta = 77
           ELSEIF(br .le. 3.506/4.696) then
             nsta = 87
           ELSEIF(br .le. 3.866/4.696) then
             nsta = 101
           ELSEIF(br .le. 4.186/4.696) then
             nsta = 116
           ELSEIF(br .le. 4.386/4.696) then
             nsta = 120
           ELSEIF(br .le. 4.556/4.696) then
             nsta = 124
           ELSE
             nsta = 135
           ENDIF
         Case(72179)
               IF(br .le. 1.97/7.5041) then
             nsta = 5
           ELSEIF(br .le. 2.26/7.5041) then
             nsta = 6
           ELSEIF(br .le. 2.2613/7.5041) then
             nsta = 12
           ELSEIF(br .le. 2.3953/7.5041) then
             nsta = 16
           ELSEIF(br .le. 2.6053/7.5041) then
             nsta = 17
           ELSEIF(br .le. 2.6213/7.5041) then
             nsta = 19
           ELSEIF(br .le. 2.6264/7.5041) then
             nsta = 35
           ELSEIF(br .le. 2.655/7.5041) then
             nsta = 54
           ELSEIF(br .le. 2.6604/7.5041) then
             nsta = 58
           ELSEIF(br .le. 2.6664/7.5041) then
             nsta = 59
           ELSEIF(br .le. 2.6868/7.5041) then
             nsta = 61
           ELSEIF(br .le. 2.6898/7.5041) then
             nsta = 80
           ELSEIF(br .le. 2.6984/7.5041) then
             nsta = 85
           ELSEIF(br .le. 2.7124/7.5041) then
             nsta = 92
           ELSEIF(br .le. 2.739/7.5041) then
             nsta = 96
           ELSEIF(br .le. 2.7539/7.5041) then
             nsta = 102
           ELSEIF(br .le. 2.7794/7.5041) then
             nsta = 108
           ELSEIF(br .le. 2.7911/7.5041) then
             nsta = 111
           ELSEIF(br .le. 2.9411/7.5041) then
             nsta = 115
           ELSEIF(br .le. 2.975/7.5041) then
             nsta = 117
           ELSEIF(br .le. 3.226/7.5041) then
             nsta = 118
           ELSEIF(br .le. 3.33/7.5041) then
             nsta = 119
           ELSEIF(br .le. 3.77/7.5041) then
             nsta = 122
           ELSEIF(br .le. 4.12/7.5041) then
             nsta = 125
           ELSEIF(br .le. 4.1301/7.5041) then
             nsta = 127
           ELSEIF(br .le. 4.1356/7.5041) then
             nsta = 128
           ELSEIF(br .le. 4.1491/7.5041) then
             nsta = 129
           ELSEIF(br .le. 4.2531/7.5041) then
             nsta = 130
           ELSEIF(br .le. 4.2921/7.5041) then
             nsta = 132
           ELSEIF(br .le. 4.2967/7.5041) then
             nsta = 135
           ELSEIF(br .le. 4.3058/7.5041) then
             nsta = 136
           ELSEIF(br .le. 4.3379/7.5041) then
             nsta = 139
           ELSEIF(br .le. 4.3402/7.5041) then
             nsta = 142
           ELSEIF(br .le. 4.3454/7.5041) then
             nsta = 143
           ELSEIF(br .le. 4.4054/7.5041) then
             nsta = 145
           ELSEIF(br .le. 4.4236/7.5041) then
             nsta = 147
           ELSEIF(br .le. 4.4556/7.5041) then
             nsta = 149
           ELSEIF(br .le. 4.4711/7.5041) then
             nsta = 150
           ELSEIF(br .le. 4.5171/7.5041) then
             nsta = 161
           ELSEIF(br .le. 4.5561/7.5041) then
             nsta = 163
           ELSEIF(br .le. 4.6281/7.5041) then
             nsta = 164
           ELSEIF(br .le. 4.7311/7.5041) then
             nsta = 165
           ELSEIF(br .le. 4.7701/7.5041) then
             nsta = 168
           ELSEIF(br .le. 4.8091/7.5041) then
             nsta = 169
           ELSEIF(br .le. 4.8431/7.5041) then
             nsta = 170
           ELSEIF(br .le. 4.9061/7.5041) then
             nsta = 171
           ELSEIF(br .le. 5.0661/7.5041) then
             nsta = 173
           ELSEIF(br .le. 5.1771/7.5041) then
             nsta = 175
           ELSEIF(br .le. 5.2131/7.5041) then
             nsta = 180
           ELSEIF(br .le. 5.3731/7.5041) then
             nsta = 182
           ELSEIF(br .le. 5.4161/7.5041) then
             nsta = 184
           ELSEIF(br .le. 5.5361/7.5041) then
             nsta = 187
           ELSEIF(br .le. 5.6341/7.5041) then
             nsta = 189
           ELSEIF(br .le. 5.8141/7.5041) then
             nsta = 192
           ELSEIF(br .le. 5.9041/7.5041) then
             nsta = 194
           ELSEIF(br .le. 6.0141/7.5041) then
             nsta = 195
           ELSEIF(br .le. 6.3241/7.5041) then
             nsta = 199
           ELSEIF(br .le. 6.5541/7.5041) then
             nsta = 200
           ELSEIF(br .le. 6.6141/7.5041) then
             nsta = 204
           ELSEIF(br .le. 6.7041/7.5041) then
             nsta = 206
           ELSEIF(br .le. 6.7641/7.5041) then
             nsta = 215
           ELSEIF(br .le. 6.8941/7.5041) then
             nsta = 219
           ELSEIF(br .le. 6.9841/7.5041) then
             nsta = 221
           ELSEIF(br .le. 7.1341/7.5041) then
             nsta = 223
           ELSEIF(br .le. 7.2741/7.5041) then
             nsta = 227
           ELSE
             nsta = 231
           ENDIF
         Case(72180)
               IF(br .le. 0.0035/2.00847) then
             nsta = 2
           ELSEIF(br .le. 0.00393/2.00847) then
             nsta = 3
           ELSEIF(br .le. 0.00417/2.00847) then
             nsta = 10
           ELSEIF(br .le. 0.00531/2.00847) then
             nsta = 13
           ELSEIF(br .le. 0.00569/2.00847) then
             nsta = 14
           ELSEIF(br .le. 0.00819/2.00847) then
             nsta = 17
           ELSEIF(br .le. 0.01019/2.00847) then
             nsta = 18
           ELSEIF(br .le. 0.01189/2.00847) then
             nsta = 22
           ELSEIF(br .le. 0.01709/2.00847) then
             nsta = 24
           ELSEIF(br .le. 0.01859/2.00847) then
             nsta = 27
           ELSEIF(br .le. 0.03659/2.00847) then
             nsta = 31
           ELSEIF(br .le. 0.03713/2.00847) then
             nsta = 33
           ELSEIF(br .le. 0.05613/2.00847) then
             nsta = 36
           ELSEIF(br .le. 0.06313/2.00847) then
             nsta = 38
           ELSEIF(br .le. 0.06352/2.00847) then
             nsta = 40
           ELSEIF(br .le. 0.06912/2.00847) then
             nsta = 48
           ELSEIF(br .le. 0.06957/2.00847) then
             nsta = 50
           ELSEIF(br .le. 0.44957/2.00847) then
             nsta = 51
           ELSEIF(br .le. 0.47657/2.00847) then
             nsta = 57
           ELSEIF(br .le. 0.48447/2.00847) then
             nsta = 66
           ELSEIF(br .le. 0.59447/2.00847) then
             nsta = 76
           ELSEIF(br .le. 0.59697/2.00847) then
             nsta = 78
           ELSEIF(br .le. 0.85697/2.00847) then
             nsta = 80
           ELSEIF(br .le. 0.86677/2.00847) then
             nsta = 87
           ELSEIF(br .le. 1.15677/2.00847) then
             nsta = 88
           ELSEIF(br .le. 1.30677/2.00847) then
             nsta = 93
           ELSEIF(br .le. 1.39677/2.00847) then
             nsta = 95
           ELSEIF(br .le. 1.39947/2.00847) then
             nsta = 102
           ELSEIF(br .le. 1.40217/2.00847) then
             nsta = 103
           ELSEIF(br .le. 1.75217/2.00847) then
             nsta = 116
           ELSEIF(br .le. 1.75847/2.00847) then
             nsta = 124
           ELSEIF(br .le. 1.80847/2.00847) then
             nsta = 170
           ELSE
             nsta = 219
           ENDIF
         Case(72181)
               IF(br .le. 1.09/3.1072) then
             nsta = 0
           ELSEIF(br .le. 1.42/3.1072) then
             nsta = 1
           ELSEIF(br .le. 1.48/3.1072) then
             nsta = 4
           ELSEIF(br .le. 1.57/3.1072) then
             nsta = 26
           ELSEIF(br .le. 1.591/3.1072) then
             nsta = 28
           ELSEIF(br .le. 1.6/3.1072) then
             nsta = 29
           ELSEIF(br .le. 1.635/3.1072) then
             nsta = 33
           ELSEIF(br .le. 1.654/3.1072) then
             nsta = 37
           ELSEIF(br .le. 1.734/3.1072) then
             nsta = 43
           ELSEIF(br .le. 1.7382/3.1072) then
             nsta = 52
           ELSEIF(br .le. 1.7452/3.1072) then
             nsta = 55
           ELSEIF(br .le. 1.7702/3.1072) then
             nsta = 58
           ELSEIF(br .le. 1.7782/3.1072) then
             nsta = 65
           ELSEIF(br .le. 1.8382/3.1072) then
             nsta = 66
           ELSEIF(br .le. 1.8572/3.1072) then
             nsta = 71
           ELSEIF(br .le. 2.0472/3.1072) then
             nsta = 75
           ELSEIF(br .le. 2.0882/3.1072) then
             nsta = 77
           ELSEIF(br .le. 2.1092/3.1072) then
             nsta = 81
           ELSEIF(br .le. 2.1422/3.1072) then
             nsta = 84
           ELSEIF(br .le. 2.2722/3.1072) then
             nsta = 88
           ELSEIF(br .le. 2.2942/3.1072) then
             nsta = 90
           ELSEIF(br .le. 2.3442/3.1072) then
             nsta = 93
           ELSEIF(br .le. 2.3712/3.1072) then
             nsta = 96
           ELSEIF(br .le. 2.4812/3.1072) then
             nsta = 97
           ELSEIF(br .le. 2.5102/3.1072) then
             nsta = 99
           ELSEIF(br .le. 2.5242/3.1072) then
             nsta = 102
           ELSEIF(br .le. 2.5392/3.1072) then
             nsta = 104
           ELSEIF(br .le. 2.6192/3.1072) then
             nsta = 111
           ELSEIF(br .le. 2.6462/3.1072) then
             nsta = 115
           ELSEIF(br .le. 2.6732/3.1072) then
             nsta = 117
           ELSEIF(br .le. 2.6842/3.1072) then
             nsta = 118
           ELSEIF(br .le. 2.7192/3.1072) then
             nsta = 123
           ELSEIF(br .le. 2.7532/3.1072) then
             nsta = 124
           ELSEIF(br .le. 2.7702/3.1072) then
             nsta = 126
           ELSEIF(br .le. 2.7792/3.1072) then
             nsta = 129
           ELSEIF(br .le. 2.7872/3.1072) then
             nsta = 130
           ELSEIF(br .le. 2.7932/3.1072) then
             nsta = 130
           ELSEIF(br .le. 2.8012/3.1072) then
             nsta = 131
           ELSEIF(br .le. 2.8512/3.1072) then
             nsta = 132
           ELSEIF(br .le. 2.8842/3.1072) then
             nsta = 136
           ELSEIF(br .le. 2.9002/3.1072) then
             nsta = 137
           ELSEIF(br .le. 2.9162/3.1072) then
             nsta = 140
           ELSEIF(br .le. 2.9272/3.1072) then
             nsta = 141
           ELSEIF(br .le. 2.9872/3.1072) then
             nsta = 142
           ELSEIF(br .le. 2.9962/3.1072) then
             nsta = 149
           ELSEIF(br .le. 3.0072/3.1072) then
             nsta = 150
           ELSEIF(br .le. 3.0222/3.1072) then
             nsta = 154
           ELSEIF(br .le. 3.0452/3.1072) then
             nsta = 156
           ELSEIF(br .le. 3.0662/3.1072) then
             nsta = 158
           ELSE
             nsta = 159
           ENDIF
         Case(73182) 
             nsta = 179
         !Case(73182)
         !      IF(br .le. 0.087/1.5374) then
         !    nsta = 0
         !  ELSEIF(br .le. 0.225/1.5374) then
         !    nsta = 2
         !  ELSEIF(br .le. 0.245/1.5374) then
         !    nsta = 3
         !  ELSEIF(br .le. 0.2516/1.5374) then
         !    nsta = 6
         !  ELSEIF(br .le. 0.2596/1.5374) then
         !    nsta = 9
         !  ELSEIF(br .le. 0.2936/1.5374) then
         !    nsta = 11
         !  ELSEIF(br .le. 0.3063/1.5374) then
         !    nsta = 12
         !  ELSEIF(br .le. 0.3233/1.5374) then
         !    nsta = 28
         !  ELSEIF(br .le. 0.3433/1.5374) then
         !    nsta = 30
         !  ELSEIF(br .le. 0.3603/1.5374) then
         !    nsta = 32
         !  ELSEIF(br .le. 0.3713/1.5374) then
         !    nsta = 38
         !  ELSEIF(br .le. 0.38/1.5374) then
         !    nsta = 41
         !  ELSEIF(br .le. 0.3887/1.5374) then
         !    nsta = 45
         !  ELSEIF(br .le. 0.4367/1.5374) then
         !    nsta = 49
         !  ELSEIF(br .le. 0.4481/1.5374) then
         !    nsta = 55
         !  ELSEIF(br .le. 0.4991/1.5374) then
         !    nsta = 58
         !  ELSEIF(br .le. 0.5082/1.5374) then
         !    nsta = 69
         !  ELSEIF(br .le. 0.5282/1.5374) then
         !    nsta = 71
         !  ELSEIF(br .le. 0.5452/1.5374) then
         !    nsta = 76
         !  ELSEIF(br .le. 0.5872/1.5374) then
         !    nsta = 78
         !  ELSEIF(br .le. 0.6202/1.5374) then
         !    nsta = 79
         !  ELSEIF(br .le. 0.63/1.5374) then
         !    nsta = 81
         !  ELSEIF(br .le. 0.6424/1.5374) then
         !    nsta = 83
         !  ELSEIF(br .le. 0.6654/1.5374) then
         !    nsta = 85
         !  ELSEIF(br .le. 0.6954/1.5374) then
         !    nsta = 90
         !  ELSEIF(br .le. 0.7324/1.5374) then
         !    nsta = 92
         !  ELSEIF(br .le. 0.7804/1.5374) then
         !    nsta = 93
         !  ELSEIF(br .le. 0.8854/1.5374) then
         !    nsta = 94
         !  ELSEIF(br .le. 0.9254/1.5374) then
         !    nsta = 105
         !  ELSEIF(br .le. 0.9384/1.5374) then
         !    nsta = 107
         !  ELSEIF(br .le. 0.9824/1.5374) then
         !    nsta = 112
         !  ELSEIF(br .le. 1.0174/1.5374) then
         !    nsta = 116
         !  ELSEIF(br .le. 1.0494/1.5374) then
         !    nsta = 118
         !  ELSEIF(br .le. 1.0814/1.5374) then
         !    nsta = 119
         !  ELSEIF(br .le. 1.1234/1.5374) then
         !    nsta = 127
         !  ELSEIF(br .le. 1.1344/1.5374) then
         !    nsta = 129
         !  ELSEIF(br .le. 1.1654/1.5374) then
         !    nsta = 131
         !  ELSEIF(br .le. 1.1884/1.5374) then
         !    nsta = 141
         !  ELSEIF(br .le. 1.2104/1.5374) then
         !    nsta = 142
         !  ELSEIF(br .le. 1.2944/1.5374) then
         !    nsta = 148
         !  ELSEIF(br .le. 1.3314/1.5374) then
         !    nsta = 156
         !  ELSEIF(br .le. 1.3454/1.5374) then
         !    nsta = 157
         !  ELSEIF(br .le. 1.3584/1.5374) then
         !    nsta = 158
         !  ELSEIF(br .le. 1.3764/1.5374) then
         !    nsta = 159
         !  ELSEIF(br .le. 1.3874/1.5374) then
         !    nsta = 161
         !  ELSEIF(br .le. 1.4214/1.5374) then
         !    nsta = 163
         !  ELSEIF(br .le. 1.4514/1.5374) then
         !    nsta = 164
         !  ELSEIF(br .le. 1.4614/1.5374) then
         !    nsta = 165
         !  ELSEIF(br .le. 1.4934/1.5374) then
         !    nsta = 168
         !  ELSEIF(br .le. 1.5164/1.5374) then
         !    nsta = 170
         !  ELSE
         !    nsta = 171
         !  ENDIF
         Case(74183)
               IF(br .le. 0.45/1.23498) then
             nsta = 0
           ELSEIF(br .le. 0.624/1.23498) then
             nsta = 1
           ELSEIF(br .le. 0.624191/1.23498) then
             nsta = 3
           ELSEIF(br .le. 0.636391/1.23498) then
             nsta = 25
           ELSEIF(br .le. 0.826391/1.23498) then
             nsta = 30
           ELSEIF(br .le. 0.832091/1.23498) then
             nsta = 42
           ELSEIF(br .le. 0.832501/1.23498) then
             nsta = 52
           ELSEIF(br .le. 0.834301/1.23498) then
             nsta = 55
           ELSEIF(br .le. 0.853201/1.23498) then
             nsta = 55
           ELSEIF(br .le. 0.868201/1.23498) then
             nsta = 60
           ELSEIF(br .le. 0.871201/1.23498) then
             nsta = 62
           ELSEIF(br .le. 0.874101/1.23498) then
             nsta = 64
           ELSEIF(br .le. 0.900101/1.23498) then
             nsta = 65
           ELSEIF(br .le. 0.905301/1.23498) then
             nsta = 66
           ELSEIF(br .le. 0.908401/1.23498) then
             nsta = 69
           ELSEIF(br .le. 0.947401/1.23498) then
             nsta = 71
           ELSEIF(br .le. 0.949201/1.23498) then
             nsta = 76
           ELSEIF(br .le. 0.951301/1.23498) then
             nsta = 77
           ELSEIF(br .le. 0.963701/1.23498) then
             nsta = 78
           ELSEIF(br .le. 0.965601/1.23498) then
             nsta = 79
           ELSEIF(br .le. 0.966501/1.23498) then
             nsta = 80
           ELSEIF(br .le. 0.968401/1.23498) then
             nsta = 82
           ELSEIF(br .le. 0.985401/1.23498) then
             nsta = 84
           ELSEIF(br .le. 1.0114/1.23498) then
             nsta = 87
           ELSEIF(br .le. 1.014/1.23498) then
             nsta = 88
           ELSEIF(br .le. 1.034/1.23498) then
             nsta = 91
           ELSEIF(br .le. 1.0383/1.23498) then
             nsta = 91
           ELSEIF(br .le. 1.0418/1.23498) then
             nsta = 91
           ELSEIF(br .le. 1.0436/1.23498) then
             nsta = 92
           ELSEIF(br .le. 1.0866/1.23498) then
             nsta = 93
           ELSEIF(br .le. 1.08668/1.23498) then
             nsta = 97
           ELSEIF(br .le. 1.09068/1.23498) then
             nsta = 98
           ELSEIF(br .le. 1.10288/1.23498) then
             nsta = 101
           ELSEIF(br .le. 1.10898/1.23498) then
             nsta = 103
           ELSEIF(br .le. 1.11638/1.23498) then
             nsta = 106
           ELSEIF(br .le. 1.13438/1.23498) then
             nsta = 107
           ELSEIF(br .le. 1.13818/1.23498) then
             nsta = 110
           ELSEIF(br .le. 1.15718/1.23498) then
             nsta = 110
           ELSEIF(br .le. 1.16118/1.23498) then
             nsta = 110
           ELSEIF(br .le. 1.21118/1.23498) then
             nsta = 111
           ELSEIF(br .le. 1.21818/1.23498) then
             nsta = 112
           ELSEIF(br .le. 1.22598/1.23498) then
             nsta = 113
           ELSE
             nsta = 114
           ENDIF
         Case(74184)
               IF(br .le. 0.071/0.58979) then
             nsta = 0
           ELSEIF(br .le. 0.0869/0.58979) then
             nsta = 1
           ELSEIF(br .le. 0.0967/0.58979) then
             nsta = 4
           ELSEIF(br .le. 0.1397/0.58979) then
             nsta = 5
           ELSEIF(br .le. 0.1632/0.58979) then
             nsta = 7
           ELSEIF(br .le. 0.1648/0.58979) then
             nsta = 8
           ELSEIF(br .le. 0.1912/0.58979) then
             nsta = 10
           ELSEIF(br .le. 0.1919/0.58979) then
             nsta = 15
           ELSEIF(br .le. 0.2279/0.58979) then
             nsta = 18
           ELSEIF(br .le. 0.2316/0.58979) then
             nsta = 20
           ELSEIF(br .le. 0.2477/0.58979) then
             nsta = 30
           ELSEIF(br .le. 0.2707/0.58979) then
             nsta = 31
           ELSEIF(br .le. 0.2733/0.58979) then
             nsta = 32
           ELSEIF(br .le. 0.2744/0.58979) then
             nsta = 38
           ELSEIF(br .le. 0.2754/0.58979) then
             nsta = 43
           ELSEIF(br .le. 0.2771/0.58979) then
             nsta = 45
           ELSEIF(br .le. 0.2881/0.58979) then
             nsta = 48
           ELSEIF(br .le. 0.28824/0.58979) then
             nsta = 52
           ELSEIF(br .le. 0.28944/0.58979) then
             nsta = 53
           ELSEIF(br .le. 0.29104/0.58979) then
             nsta = 55
           ELSEIF(br .le. 0.29364/0.58979) then
             nsta = 56
           ELSEIF(br .le. 0.29409/0.58979) then
             nsta = 58
           ELSEIF(br .le. 0.29679/0.58979) then
             nsta = 59
           ELSEIF(br .le. 0.29969/0.58979) then
             nsta = 63
           ELSEIF(br .le. 0.30189/0.58979) then
             nsta = 64
           ELSEIF(br .le. 0.30289/0.58979) then
             nsta = 65
           ELSEIF(br .le. 0.30939/0.58979) then
             nsta = 67
           ELSEIF(br .le. 0.32089/0.58979) then
             nsta = 69
           ELSEIF(br .le. 0.32689/0.58979) then
             nsta = 70
           ELSEIF(br .le. 0.32809/0.58979) then
             nsta = 74
           ELSEIF(br .le. 0.32884/0.58979) then
             nsta = 75
           ELSEIF(br .le. 0.33124/0.58979) then
             nsta = 76
           ELSEIF(br .le. 0.34264/0.58979) then
             nsta = 77
           ELSEIF(br .le. 0.34334/0.58979) then
             nsta = 79
           ELSEIF(br .le. 0.34396/0.58979) then
             nsta = 82
           ELSEIF(br .le. 0.34636/0.58979) then
             nsta = 83
           ELSEIF(br .le. 0.35086/0.58979) then
             nsta = 84
           ELSEIF(br .le. 0.36706/0.58979) then
             nsta = 87
           ELSEIF(br .le. 0.37026/0.58979) then
             nsta = 89
           ELSEIF(br .le. 0.37063/0.58979) then
             nsta = 91
           ELSEIF(br .le. 0.37303/0.58979) then
             nsta = 92
           ELSEIF(br .le. 0.38203/0.58979) then
             nsta = 93
           ELSEIF(br .le. 0.38313/0.58979) then
             nsta = 97
           ELSEIF(br .le. 0.38473/0.58979) then
             nsta = 99
           ELSEIF(br .le. 0.38753/0.58979) then
             nsta = 101
           ELSEIF(br .le. 0.38853/0.58979) then
             nsta = 105
           ELSEIF(br .le. 0.38915/0.58979) then
             nsta = 108
           ELSEIF(br .le. 0.38959/0.58979) then
             nsta = 110
           ELSEIF(br .le. 0.39359/0.58979) then
             nsta = 111
           ELSEIF(br .le. 0.39469/0.58979) then
             nsta = 112
           ELSEIF(br .le. 0.39659/0.58979) then
             nsta = 113
           ELSEIF(br .le. 0.39819/0.58979) then
             nsta = 114
           ELSEIF(br .le. 0.40139/0.58979) then
             nsta = 120
           ELSEIF(br .le. 0.40349/0.58979) then
             nsta = 122
           ELSEIF(br .le. 0.40689/0.58979) then
             nsta = 123
           ELSEIF(br .le. 0.40959/0.58979) then
             nsta = 125
           ELSEIF(br .le. 0.41839/0.58979) then
             nsta = 126
           ELSEIF(br .le. 0.41969/0.58979) then
             nsta = 128
           ELSEIF(br .le. 0.42609/0.58979) then
             nsta = 129
           ELSEIF(br .le. 0.43189/0.58979) then
             nsta = 131
           ELSEIF(br .le. 0.43359/0.58979) then
             nsta = 135
           ELSEIF(br .le. 0.43459/0.58979) then
             nsta = 137
           ELSEIF(br .le. 0.43579/0.58979) then
             nsta = 138
           ELSEIF(br .le. 0.43749/0.58979) then
             nsta = 140
           ELSEIF(br .le. 0.43979/0.58979) then
             nsta = 142
           ELSEIF(br .le. 0.44269/0.58979) then
             nsta = 147
           ELSEIF(br .le. 0.44689/0.58979) then
             nsta = 148
           ELSEIF(br .le. 0.44769/0.58979) then
             nsta = 150
           ELSEIF(br .le. 0.44879/0.58979) then
             nsta = 151
           ELSEIF(br .le. 0.45259/0.58979) then
             nsta = 152
           ELSEIF(br .le. 0.45859/0.58979) then
             nsta = 154
           ELSEIF(br .le. 0.46229/0.58979) then
             nsta = 156
           ELSEIF(br .le. 0.46499/0.58979) then
             nsta = 159
           ELSEIF(br .le. 0.47339/0.58979) then
             nsta = 163
           ELSEIF(br .le. 0.47599/0.58979) then
             nsta = 165
           ELSEIF(br .le. 0.47969/0.58979) then
             nsta = 169
           ELSEIF(br .le. 0.48339/0.58979) then
             nsta = 172
           ELSEIF(br .le. 0.48999/0.58979) then
             nsta = 174
           ELSEIF(br .le. 0.49149/0.58979) then
             nsta = 176
           ELSEIF(br .le. 0.49349/0.58979) then
             nsta = 177
           ELSEIF(br .le. 0.49539/0.58979) then
             nsta = 179
           ELSEIF(br .le. 0.49769/0.58979) then
             nsta = 181
           ELSEIF(br .le. 0.49969/0.58979) then
             nsta = 183
           ELSEIF(br .le. 0.50459/0.58979) then
             nsta = 184
           ELSEIF(br .le. 0.50819/0.58979) then
             nsta = 186
           ELSEIF(br .le. 0.51199/0.58979) then
             nsta = 189
           ELSEIF(br .le. 0.51339/0.58979) then
             nsta = 192
           ELSEIF(br .le. 0.51409/0.58979) then
             nsta = 194
           ELSEIF(br .le. 0.51569/0.58979) then
             nsta = 197
           ELSEIF(br .le. 0.51689/0.58979) then
             nsta = 201
           ELSEIF(br .le. 0.51759/0.58979) then
             nsta = 202
           ELSEIF(br .le. 0.51879/0.58979) then
             nsta = 204
           ELSEIF(br .le. 0.52039/0.58979) then
             nsta = 205
           ELSEIF(br .le. 0.52309/0.58979) then
             nsta = 206
           ELSEIF(br .le. 0.52449/0.58979) then
             nsta = 212
           ELSEIF(br .le. 0.52689/0.58979) then
             nsta = 215
           ELSEIF(br .le. 0.52919/0.58979) then
             nsta = 216
           ELSEIF(br .le. 0.53029/0.58979) then
             nsta = 219
           ELSEIF(br .le. 0.53229/0.58979) then
             nsta = 220
           ELSEIF(br .le. 0.54329/0.58979) then
             nsta = 224
           ELSEIF(br .le. 0.54479/0.58979) then
             nsta = 226
           ELSEIF(br .le. 0.54749/0.58979) then
             nsta = 228
           ELSEIF(br .le. 0.54969/0.58979) then
             nsta = 230
           ELSEIF(br .le. 0.55179/0.58979) then
             nsta = 234
           ELSEIF(br .le. 0.55879/0.58979) then
             nsta = 235
           ELSEIF(br .le. 0.56979/0.58979) then
             nsta = 237
           ELSEIF(br .le. 0.57249/0.58979) then
             nsta = 239
           ELSEIF(br .le. 0.57459/0.58979) then
             nsta = 240
           ELSEIF(br .le. 0.57589/0.58979) then
             nsta = 241
           ELSEIF(br .le. 0.57829/0.58979) then
             nsta = 242
           ELSEIF(br .le. 0.58179/0.58979) then
             nsta = 244
           ELSEIF(br .le. 0.58289/0.58979) then
             nsta = 245
           ELSEIF(br .le. 0.58549/0.58979) then
             nsta = 246
           ELSEIF(br .le. 0.58709/0.58979) then
             nsta = 247
           ELSE
             nsta = 248
           ENDIF
         Case(74185)
               IF(br .le. 0.0112/0.4856) then
             nsta = 0
           ELSEIF(br .le. 0.0147/0.4856) then
             nsta = 1
           ELSEIF(br .le. 0.0165/0.4856) then
             nsta = 3
           ELSEIF(br .le. 0.0865/0.4856) then
             nsta = 15
           ELSEIF(br .le. 0.0935/0.4856) then
             nsta = 18
           ELSEIF(br .le. 0.1125/0.4856) then
             nsta = 19
           ELSEIF(br .le. 0.1295/0.4856) then
             nsta = 20
           ELSEIF(br .le. 0.1414/0.4856) then
             nsta = 23
           ELSEIF(br .le. 0.2014/0.4856) then
             nsta = 34
           ELSEIF(br .le. 0.2094/0.4856) then
             nsta = 40
           ELSEIF(br .le. 0.3134/0.4856) then
             nsta = 54
           ELSEIF(br .le. 0.3934/0.4856) then
             nsta = 55
           ELSEIF(br .le. 0.3975/0.4856) then
             nsta = 57
           ELSEIF(br .le. 0.4195/0.4856) then
             nsta = 59
           ELSEIF(br .le. 0.4275/0.4856) then
             nsta = 68
           ELSEIF(br .le. 0.4345/0.4856) then
             nsta = 78
           ELSEIF(br .le. 0.4685/0.4856) then
             nsta = 84
           ELSEIF(br .le. 0.4745/0.4856) then
             nsta = 88
           ELSEIF(br .le. 0.4786/0.4856) then
             nsta = 97
           ELSE
             nsta = 129
           ENDIF
         Case(74187)
               IF(br .le. 0.023/4.0371) then
             nsta = 0
           ELSEIF(br .le. 0.628/4.0371) then
             nsta = 2
           ELSEIF(br .le. 1.488/4.0371) then
             nsta = 4
           ELSEIF(br .le. 1.638/4.0371) then
             nsta = 23
           ELSEIF(br .le. 1.69/4.0371) then
             nsta = 30
           ELSEIF(br .le. 1.814/4.0371) then
             nsta = 31
           ELSEIF(br .le. 1.966/4.0371) then
             nsta = 38
           ELSEIF(br .le. 1.991/4.0371) then
             nsta = 40
           ELSEIF(br .le. 2.027/4.0371) then
             nsta = 46
           ELSEIF(br .le. 2.075/4.0371) then
             nsta = 52
           ELSEIF(br .le. 2.132/4.0371) then
             nsta = 59
           ELSEIF(br .le. 2.172/4.0371) then
             nsta = 64
           ELSEIF(br .le. 2.287/4.0371) then
             nsta = 71
           ELSEIF(br .le. 2.33/4.0371) then
             nsta = 80
           ELSEIF(br .le. 2.364/4.0371) then
             nsta = 85
           ELSEIF(br .le. 2.423/4.0371) then
             nsta = 86
           ELSEIF(br .le. 2.474/4.0371) then
             nsta = 91
           ELSEIF(br .le. 2.503/4.0371) then
             nsta = 97
           ELSEIF(br .le. 2.537/4.0371) then
             nsta = 102
           ELSEIF(br .le. 2.554/4.0371) then
             nsta = 105
           ELSEIF(br .le. 2.578/4.0371) then
             nsta = 108
           ELSEIF(br .le. 2.592/4.0371) then
             nsta = 109
           ELSEIF(br .le. 2.643/4.0371) then
             nsta = 115
           ELSEIF(br .le. 2.663/4.0371) then
             nsta = 121
           ELSEIF(br .le. 2.689/4.0371) then
             nsta = 125
           ELSEIF(br .le. 2.715/4.0371) then
             nsta = 126
           ELSEIF(br .le. 2.784/4.0371) then
             nsta = 131
           ELSEIF(br .le. 2.818/4.0371) then
             nsta = 136
           ELSEIF(br .le. 2.834/4.0371) then
             nsta = 151
           ELSEIF(br .le. 2.894/4.0371) then
             nsta = 154
           ELSEIF(br .le. 2.957/4.0371) then
             nsta = 156
           ELSEIF(br .le. 2.997/4.0371) then
             nsta = 157
           ELSEIF(br .le. 3.03/4.0371) then
             nsta = 162
           ELSEIF(br .le. 3.081/4.0371) then
             nsta = 165
           ELSEIF(br .le. 3.184/4.0371) then
             nsta = 168
           ELSEIF(br .le. 3.239/4.0371) then
             nsta = 170
           ELSEIF(br .le. 3.278/4.0371) then
             nsta = 172
           ELSEIF(br .le. 3.308/4.0371) then
             nsta = 176
           ELSEIF(br .le. 3.349/4.0371) then
             nsta = 182
           ELSEIF(br .le. 3.364/4.0371) then
             nsta = 192
           ELSEIF(br .le. 3.3741/4.0371) then
             nsta = 197
           ELSEIF(br .le. 3.4161/4.0371) then
             nsta = 203
           ELSEIF(br .le. 3.4461/4.0371) then
             nsta = 206
           ELSEIF(br .le. 3.4831/4.0371) then
             nsta = 210
           ELSEIF(br .le. 3.5441/4.0371) then
             nsta = 219
           ELSEIF(br .le. 3.6301/4.0371) then
             nsta = 219
           ELSEIF(br .le. 3.6551/4.0371) then
             nsta = 223
           ELSEIF(br .le. 3.6701/4.0371) then
             nsta = 226
           ELSEIF(br .le. 3.6991/4.0371) then
             nsta = 232
           ELSEIF(br .le. 3.7071/4.0371) then
             nsta = 242
           ELSEIF(br .le. 3.7211/4.0371) then
             nsta = 251
           ELSEIF(br .le. 3.7541/4.0371) then
             nsta = 267
           ELSEIF(br .le. 3.7771/4.0371) then
             nsta = 270
           ELSEIF(br .le. 3.8101/4.0371) then
             nsta = 282
           ELSEIF(br .le. 3.8421/4.0371) then
             nsta = 287
           ELSEIF(br .le. 3.8631/4.0371) then
             nsta = 290
           ELSEIF(br .le. 3.8891/4.0371) then
             nsta = 292
           ELSEIF(br .le. 3.9131/4.0371) then
             nsta = 294
           ELSEIF(br .le. 3.9441/4.0371) then
             nsta = 302
           ELSEIF(br .le. 3.9651/4.0371) then
             nsta = 304
           ELSEIF(br .le. 3.9961/4.0371) then
             nsta = 306
           ELSEIF(br .le. 4.0261/4.0371) then
             nsta = 307
           ELSE
             nsta = 308
           ENDIF
         Case(75186)
               IF(br .le. 0.01/6.6) then
             nsta = 0
           ELSEIF(br .le. 0.192/6.6) then
             nsta = 1
           ELSEIF(br .le. 0.347/6.6) then
             nsta = 2
           ELSEIF(br .le. 0.437/6.6) then
             nsta = 3
           ELSEIF(br .le. 0.518/6.6) then
             nsta = 5
           ELSEIF(br .le. 0.531/6.6) then
             nsta = 7
           ELSEIF(br .le. 1.131/6.6) then
             nsta = 8
           ELSEIF(br .le. 1.271/6.6) then
             nsta = 13
           ELSEIF(br .le. 1.285/6.6) then
             nsta = 17
           ELSEIF(br .le. 1.296/6.6) then
             nsta = 19
           ELSEIF(br .le. 1.394/6.6) then
             nsta = 22
           ELSEIF(br .le. 1.509/6.6) then
             nsta = 23
           ELSEIF(br .le. 1.597/6.6) then
             nsta = 27
           ELSEIF(br .le. 1.706/6.6) then
             nsta = 31
           ELSEIF(br .le. 1.82/6.6) then
             nsta = 39
           ELSEIF(br .le. 1.839/6.6) then
             nsta = 48
           ELSEIF(br .le. 1.865/6.6) then
             nsta = 51
           ELSEIF(br .le. 1.995/6.6) then
             nsta = 53
           ELSEIF(br .le. 2.029/6.6) then
             nsta = 56
           ELSEIF(br .le. 2.055/6.6) then
             nsta = 60
           ELSEIF(br .le. 2.195/6.6) then
             nsta = 62
           ELSEIF(br .le. 2.255/6.6) then
             nsta = 64
           ELSEIF(br .le. 2.272/6.6) then
             nsta = 66
           ELSEIF(br .le. 2.29/6.6) then
             nsta = 69
           ELSEIF(br .le. 2.36/6.6) then
             nsta = 70
           ELSEIF(br .le. 2.54/6.6) then
             nsta = 74
           ELSEIF(br .le. 2.93/6.6) then
             nsta = 78
           ELSEIF(br .le. 2.953/6.6) then
             nsta = 82
           ELSEIF(br .le. 3.047/6.6) then
             nsta = 84
           ELSEIF(br .le. 3.143/6.6) then
             nsta = 87
           ELSEIF(br .le. 3.433/6.6) then
             nsta = 89
           ELSEIF(br .le. 3.703/6.6) then
             nsta = 92
           ELSEIF(br .le. 3.733/6.6) then
             nsta = 93
           ELSEIF(br .le. 3.776/6.6) then
             nsta = 95
           ELSEIF(br .le. 3.819/6.6) then
             nsta = 101
           ELSEIF(br .le. 3.837/6.6) then
             nsta = 104
           ELSEIF(br .le. 3.947/6.6) then
             nsta = 106
           ELSEIF(br .le. 4.317/6.6) then
             nsta = 107
           ELSEIF(br .le. 4.338/6.6) then
             nsta = 109
           ELSEIF(br .le. 4.364/6.6) then
             nsta = 112
           ELSEIF(br .le. 4.544/6.6) then
             nsta = 117
           ELSEIF(br .le. 4.657/6.6) then
             nsta = 121
           ELSEIF(br .le. 4.692/6.6) then
             nsta = 125
           ELSEIF(br .le. 4.725/6.6) then
             nsta = 126
           ELSEIF(br .le. 4.965/6.6) then
             nsta = 127
           ELSEIF(br .le. 5.05/6.6) then
             nsta = 130
           ELSEIF(br .le. 5.131/6.6) then
             nsta = 132
           ELSEIF(br .le. 5.191/6.6) then
             nsta = 135
           ELSEIF(br .le. 5.219/6.6) then
             nsta = 136
           ELSEIF(br .le. 5.231/6.6) then
             nsta = 137
           ELSEIF(br .le. 5.262/6.6) then
             nsta = 138
           ELSEIF(br .le. 5.313/6.6) then
             nsta = 139
           ELSEIF(br .le. 5.343/6.6) then
             nsta = 141
           ELSEIF(br .le. 5.394/6.6) then
             nsta = 142
           ELSEIF(br .le. 5.504/6.6) then
             nsta = 146
           ELSEIF(br .le. 5.584/6.6) then
             nsta = 147
           ELSEIF(br .le. 5.644/6.6) then
             nsta = 149
           ELSEIF(br .le. 5.683/6.6) then
             nsta = 154
           ELSEIF(br .le. 5.734/6.6) then
             nsta = 155
           ELSEIF(br .le. 5.806/6.6) then
             nsta = 156
           ELSEIF(br .le. 5.84/6.6) then
             nsta = 157
           ELSEIF(br .le. 5.928/6.6) then
             nsta = 158
           ELSEIF(br .le. 6.028/6.6) then
             nsta = 159
           ELSEIF(br .le. 6.046/6.6) then
             nsta = 160
           ELSEIF(br .le. 6.176/6.6) then
             nsta = 161
           ELSEIF(br .le. 6.316/6.6) then
             nsta = 162
           ELSEIF(br .le. 6.376/6.6) then
             nsta = 163
           ELSEIF(br .le. 6.399/6.6) then
             nsta = 164
           ELSEIF(br .le. 6.43/6.6) then
             nsta = 165
           ELSEIF(br .le. 6.52/6.6) then
             nsta = 166
           ELSE
             nsta = 167
           ENDIF
         Case(75188)
               IF(br .le. 0.299/4.4823) then
             nsta = 0
           ELSEIF(br .le. 0.329/4.4823) then
             nsta = 2
           ELSEIF(br .le. 0.415/4.4823) then
             nsta = 5
           ELSEIF(br .le. 0.515/4.4823) then
             nsta = 6
           ELSEIF(br .le. 0.635/4.4823) then
             nsta = 8
           ELSEIF(br .le. 0.675/4.4823) then
             nsta = 9
           ELSEIF(br .le. 0.706/4.4823) then
             nsta = 11
           ELSEIF(br .le. 0.798/4.4823) then
             nsta = 12
           ELSEIF(br .le. 0.806/4.4823) then
             nsta = 14
           ELSEIF(br .le. 0.853/4.4823) then
             nsta = 15
           ELSEIF(br .le. 0.9/4.4823) then
             nsta = 16
           ELSEIF(br .le. 0.966/4.4823) then
             nsta = 17
           ELSEIF(br .le. 1.038/4.4823) then
             nsta = 18
           ELSEIF(br .le. 1.046/4.4823) then
             nsta = 20
           ELSEIF(br .le. 1.06/4.4823) then
             nsta = 22
           ELSEIF(br .le. 1.091/4.4823) then
             nsta = 24
           ELSEIF(br .le. 1.166/4.4823) then
             nsta = 25
           ELSEIF(br .le. 1.177/4.4823) then
             nsta = 28
           ELSEIF(br .le. 1.19/4.4823) then
             nsta = 30
           ELSEIF(br .le. 1.1953/4.4823) then
             nsta = 31
           ELSEIF(br .le. 1.2063/4.4823) then
             nsta = 32
           ELSEIF(br .le. 1.2963/4.4823) then
             nsta = 33
           ELSEIF(br .le. 1.3223/4.4823) then
             nsta = 34
           ELSEIF(br .le. 1.5223/4.4823) then
             nsta = 35
           ELSEIF(br .le. 1.6053/4.4823) then
             nsta = 38
           ELSEIF(br .le. 1.6753/4.4823) then
             nsta = 39
           ELSEIF(br .le. 1.6893/4.4823) then
             nsta = 40
           ELSEIF(br .le. 1.7003/4.4823) then
             nsta = 41
           ELSEIF(br .le. 1.7103/4.4823) then
             nsta = 42
           ELSEIF(br .le. 1.7733/4.4823) then
             nsta = 43
           ELSEIF(br .le. 1.8543/4.4823) then
             nsta = 44
           ELSEIF(br .le. 1.9003/4.4823) then
             nsta = 45
           ELSEIF(br .le. 2.0403/4.4823) then
             nsta = 46
           ELSEIF(br .le. 2.0683/4.4823) then
             nsta = 48
           ELSEIF(br .le. 2.3183/4.4823) then
             nsta = 49
           ELSEIF(br .le. 2.3433/4.4823) then
             nsta = 51
           ELSEIF(br .le. 2.4193/4.4823) then
             nsta = 52
           ELSEIF(br .le. 2.4423/4.4823) then
             nsta = 54
           ELSEIF(br .le. 2.5403/4.4823) then
             nsta = 55
           ELSEIF(br .le. 2.9703/4.4823) then
             nsta = 56
           ELSEIF(br .le. 3.0053/4.4823) then
             nsta = 57
           ELSEIF(br .le. 3.1653/4.4823) then
             nsta = 59
           ELSEIF(br .le. 3.2123/4.4823) then
             nsta = 60
           ELSEIF(br .le. 3.3103/4.4823) then
             nsta = 61
           ELSEIF(br .le. 3.4903/4.4823) then
             nsta = 62
           ELSEIF(br .le. 3.6603/4.4823) then
             nsta = 63
           ELSEIF(br .le. 3.8303/4.4823) then
             nsta = 64
           ELSEIF(br .le. 3.9803/4.4823) then
             nsta = 65
           ELSEIF(br .le. 4.1203/4.4823) then
             nsta = 66
           ELSEIF(br .le. 4.1583/4.4823) then
             nsta = 69
           ELSEIF(br .le. 4.2603/4.4823) then
             nsta = 70
           ELSEIF(br .le. 4.3413/4.4823) then
             nsta = 72
           ELSE
             nsta = 73
           ENDIF
         Case(76185)
               IF(br .le. 0.093/0.2349) then
             nsta = 1
           ELSEIF(br .le. 0.1014/0.2349) then
             nsta = 4
           ELSEIF(br .le. 0.1051/0.2349) then
             nsta = 12
           ELSEIF(br .le. 0.1491/0.2349) then
             nsta = 15
           ELSEIF(br .le. 0.1556/0.2349) then
             nsta = 23
           ELSEIF(br .le. 0.1584/0.2349) then
             nsta = 33
           ELSEIF(br .le. 0.1631/0.2349) then
             nsta = 35
           ELSEIF(br .le. 0.1721/0.2349) then
             nsta = 37
           ELSEIF(br .le. 0.1768/0.2349) then
             nsta = 40
           ELSEIF(br .le. 0.1858/0.2349) then
             nsta = 43
           ELSEIF(br .le. 0.1932/0.2349) then
             nsta = 45
           ELSEIF(br .le. 0.1997/0.2349) then
             nsta = 50
           ELSEIF(br .le. 0.2044/0.2349) then
             nsta = 52
           ELSEIF(br .le. 0.2109/0.2349) then
             nsta = 56
           ELSEIF(br .le. 0.2165/0.2349) then
             nsta = 59
           ELSEIF(br .le. 0.2239/0.2349) then
             nsta = 62
           ELSE
             nsta = 64
           ENDIF
         Case(76187)
               IF(br .le. 0.00032/0.05279) then
             nsta = 0
           ELSEIF(br .le. 0.00064/0.05279) then
             nsta = 1
           ELSEIF(br .le. 0.00204/0.05279) then
             nsta = 3
           ELSEIF(br .le. 0.05204/0.05279) then
             nsta = 23
           ELSE
             nsta = 24
           ENDIF
         Case(76188)
               IF(br .le. 0.0208/0.9007) then
             nsta = 0
           ELSEIF(br .le. 0.0455/0.9007) then
             nsta = 1
           ELSEIF(br .le. 0.059/0.9007) then
             nsta = 3
           ELSEIF(br .le. 0.0648/0.9007) then
             nsta = 7
           ELSEIF(br .le. 0.0713/0.9007) then
             nsta = 9
           ELSEIF(br .le. 0.0747/0.9007) then
             nsta = 10
           ELSEIF(br .le. 0.0887/0.9007) then
             nsta = 14
           ELSEIF(br .le. 0.1027/0.9007) then
             nsta = 16
           ELSEIF(br .le. 0.1123/0.9007) then
             nsta = 21
           ELSEIF(br .le. 0.1138/0.9007) then
             nsta = 25
           ELSEIF(br .le. 0.1189/0.9007) then
             nsta = 26
           ELSEIF(br .le. 0.1261/0.9007) then
             nsta = 28
           ELSEIF(br .le. 0.1435/0.9007) then
             nsta = 30
           ELSEIF(br .le. 0.145/0.9007) then
             nsta = 31
           ELSEIF(br .le. 0.1607/0.9007) then
             nsta = 32
           ELSEIF(br .le. 0.1733/0.9007) then
             nsta = 39
           ELSEIF(br .le. 0.1812/0.9007) then
             nsta = 41
           ELSEIF(br .le. 0.1932/0.9007) then
             nsta = 47
           ELSEIF(br .le. 0.2372/0.9007) then
             nsta = 50
           ELSEIF(br .le. 0.2461/0.9007) then
             nsta = 52
           ELSEIF(br .le. 0.2567/0.9007) then
             nsta = 54
           ELSEIF(br .le. 0.2647/0.9007) then
             nsta = 55
           ELSEIF(br .le. 0.2821/0.9007) then
             nsta = 59
           ELSEIF(br .le. 0.3321/0.9007) then
             nsta = 63
           ELSEIF(br .le. 0.3551/0.9007) then
             nsta = 68
           ELSEIF(br .le. 0.3732/0.9007) then
             nsta = 69
           ELSEIF(br .le. 0.4252/0.9007) then
             nsta = 72
           ELSEIF(br .le. 0.4632/0.9007) then
             nsta = 75
           ELSEIF(br .le. 0.4704/0.9007) then
             nsta = 76
           ELSEIF(br .le. 0.5194/0.9007) then
             nsta = 77
           ELSEIF(br .le. 0.529/0.9007) then
             nsta = 78
           ELSEIF(br .le. 0.557/0.9007) then
             nsta = 83
           ELSEIF(br .le. 0.5727/0.9007) then
             nsta = 85
           ELSEIF(br .le. 0.5888/0.9007) then
             nsta = 87
           ELSEIF(br .le. 0.6178/0.9007) then
             nsta = 88
           ELSEIF(br .le. 0.6393/0.9007) then
             nsta = 89
           ELSEIF(br .le. 0.6703/0.9007) then
             nsta = 91
           ELSEIF(br .le. 0.6846/0.9007) then
             nsta = 93
           ELSEIF(br .le. 0.7216/0.9007) then
             nsta = 96
           ELSEIF(br .le. 0.7387/0.9007) then
             nsta = 101
           ELSEIF(br .le. 0.7597/0.9007) then
             nsta = 103
           ELSEIF(br .le. 0.7967/0.9007) then
             nsta = 104
           ELSEIF(br .le. 0.8097/0.9007) then
             nsta = 105
           ELSEIF(br .le. 0.8587/0.9007) then
             nsta = 108
           ELSE
             nsta = 110
           ENDIF
         Case(76189)
               IF(br .le. 0.041/0.315) then
             nsta = 2
           ELSEIF(br .le. 0.09/0.315) then
             nsta = 14
           ELSEIF(br .le. 0.143/0.315) then
             nsta = 17
           ELSEIF(br .le. 0.166/0.315) then
             nsta = 19
           ELSEIF(br .le. 0.194/0.315) then
             nsta = 20
           ELSEIF(br .le. 0.214/0.315) then
             nsta = 46
           ELSEIF(br .le. 0.234/0.315) then
             nsta = 76
           ELSEIF(br .le. 0.244/0.315) then
             nsta = 77
           ELSEIF(br .le. 0.25/0.315) then
             nsta = 84
           ELSEIF(br .le. 0.302/0.315) then
             nsta = 88
           ELSEIF(br .le. 0.308/0.315) then
             nsta = 90
           ELSE
             nsta = 90
           ENDIF
         Case(76190)
               IF(br .le. 0.034/0.7364) then
             nsta = 0
           ELSEIF(br .le. 0.0503/0.7364) then
             nsta = 1
           ELSEIF(br .le. 0.0943/0.7364) then
             nsta = 3
           ELSEIF(br .le. 0.1073/0.7364) then
             nsta = 4
           ELSEIF(br .le. 0.1142/0.7364) then
             nsta = 8
           ELSEIF(br .le. 0.1402/0.7364) then
             nsta = 13
           ELSEIF(br .le. 0.1542/0.7364) then
             nsta = 15
           ELSEIF(br .le. 0.1802/0.7364) then
             nsta = 20
           ELSEIF(br .le. 0.1852/0.7364) then
             nsta = 22
           ELSEIF(br .le. 0.2042/0.7364) then
             nsta = 27
           ELSEIF(br .le. 0.2232/0.7364) then
             nsta = 29
           ELSEIF(br .le. 0.2432/0.7364) then
             nsta = 35
           ELSEIF(br .le. 0.2612/0.7364) then
             nsta = 39
           ELSEIF(br .le. 0.3572/0.7364) then
             nsta = 41
           ELSEIF(br .le. 0.3922/0.7364) then
             nsta = 46
           ELSEIF(br .le. 0.4222/0.7364) then
             nsta = 47
           ELSEIF(br .le. 0.4322/0.7364) then
             nsta = 49
           ELSEIF(br .le. 0.4345/0.7364) then
             nsta = 51
           ELSEIF(br .le. 0.4525/0.7364) then
             nsta = 53
           ELSEIF(br .le. 0.4765/0.7364) then
             nsta = 54
           ELSEIF(br .le. 0.5025/0.7364) then
             nsta = 56
           ELSEIF(br .le. 0.5075/0.7364) then
             nsta = 57
           ELSEIF(br .le. 0.5235/0.7364) then
             nsta = 60
           ELSEIF(br .le. 0.5685/0.7364) then
             nsta = 63
           ELSEIF(br .le. 0.5708/0.7364) then
             nsta = 68
           ELSEIF(br .le. 0.5731/0.7364) then
             nsta = 69
           ELSEIF(br .le. 0.5971/0.7364) then
             nsta = 70
           ELSEIF(br .le. 0.6021/0.7364) then
             nsta = 71
           ELSEIF(br .le. 0.6071/0.7364) then
             nsta = 72
           ELSEIF(br .le. 0.6101/0.7364) then
             nsta = 73
           ELSEIF(br .le. 0.6551/0.7364) then
             nsta = 74
           ELSEIF(br .le. 0.6761/0.7364) then
             nsta = 75
           ELSEIF(br .le. 0.6791/0.7364) then
             nsta = 77
           ELSEIF(br .le. 0.6821/0.7364) then
             nsta = 78
           ELSEIF(br .le. 0.6961/0.7364) then
             nsta = 80
           ELSEIF(br .le. 0.7101/0.7364) then
             nsta = 84
           ELSEIF(br .le. 0.7124/0.7364) then
             nsta = 89
           ELSE
             nsta = 91
           ENDIF
         Case(76191)
               IF(br .le. 0.167/0.817) then
             nsta = 1
           ELSEIF(br .le. 0.205/0.817) then
             nsta = 2
           ELSEIF(br .le. 0.279/0.817) then
             nsta = 14
           ELSEIF(br .le. 0.688/0.817) then
             nsta = 28
           ELSEIF(br .le. 0.729/0.817) then
             nsta = 37
           ELSEIF(br .le. 0.758/0.817) then
             nsta = 39
           ELSEIF(br .le. 0.782/0.817) then
             nsta = 64
           ELSE
             nsta = 74
           ENDIF
         Case(76193)
               IF(br .le. 0.076/0.342) then
             nsta = 0
           ELSEIF(br .le. 0.192/0.342) then
             nsta = 6
           ELSEIF(br .le. 0.217/0.342) then
             nsta = 13
           ELSEIF(br .le. 0.307/0.342) then
             nsta = 16
           ELSE
             nsta = 36
           ENDIF
         Case(77192)
               IF(br .le. 0.072/57.364) then
             nsta = 3
           ELSEIF(br .le. 0.102/57.364) then
             nsta = 4
           ELSEIF(br .le. 0.662/57.364) then
             nsta = 7
           ELSEIF(br .le. 3.282/57.364) then
             nsta = 8
           ELSEIF(br .le. 3.572/57.364) then
             nsta = 9
           ELSEIF(br .le. 3.62/57.364) then
             nsta = 11
           ELSEIF(br .le. 3.674/57.364) then
             nsta = 13
           ELSEIF(br .le. 3.731/57.364) then
             nsta = 16
           ELSEIF(br .le. 3.988/57.364) then
             nsta = 19
           ELSEIF(br .le. 4.038/57.364) then
             nsta = 20
           ELSEIF(br .le. 4.127/57.364) then
             nsta = 21
           ELSEIF(br .le. 4.381/57.364) then
             nsta = 24
           ELSEIF(br .le. 5.131/57.364) then
             nsta = 26
           ELSEIF(br .le. 6.921/57.364) then
             nsta = 27
           ELSEIF(br .le. 7.151/57.364) then
             nsta = 35
           ELSEIF(br .le. 7.601/57.364) then
             nsta = 36
           ELSEIF(br .le. 7.658/57.364) then
             nsta = 38
           ELSEIF(br .le. 8.448/57.364) then
             nsta = 39
           ELSEIF(br .le. 9.178/57.364) then
             nsta = 40
           ELSEIF(br .le. 9.223/57.364) then
             nsta = 42
           ELSEIF(br .le. 9.383/57.364) then
             nsta = 45
           ELSEIF(br .le. 9.496/57.364) then
             nsta = 46
           ELSEIF(br .le. 9.976/57.364) then
             nsta = 47
           ELSEIF(br .le. 11.316/57.364) then
             nsta = 50
           ELSEIF(br .le. 11.806/57.364) then
             nsta = 54
           ELSEIF(br .le. 11.996/57.364) then
             nsta = 56
           ELSEIF(br .le. 12.266/57.364) then
             nsta = 57
           ELSEIF(br .le. 12.388/57.364) then
             nsta = 59
           ELSEIF(br .le. 14.118/57.364) then
             nsta = 60
           ELSEIF(br .le. 14.283/57.364) then
             nsta = 62
           ELSEIF(br .le. 16.963/57.364) then
             nsta = 64
           ELSEIF(br .le. 17.343/57.364) then
             nsta = 65
           ELSEIF(br .le. 17.733/57.364) then
             nsta = 66
           ELSEIF(br .le. 17.813/57.364) then
             nsta = 67
           ELSEIF(br .le. 18.873/57.364) then
             nsta = 70
           ELSEIF(br .le. 19.593/57.364) then
             nsta = 71
           ELSEIF(br .le. 19.693/57.364) then
             nsta = 72
           ELSEIF(br .le. 20.363/57.364) then
             nsta = 73
           ELSEIF(br .le. 22.073/57.364) then
             nsta = 74
           ELSEIF(br .le. 22.236/57.364) then
             nsta = 75
           ELSEIF(br .le. 22.286/57.364) then
             nsta = 76
           ELSEIF(br .le. 23.676/57.364) then
             nsta = 77
           ELSEIF(br .le. 24.436/57.364) then
             nsta = 79
           ELSEIF(br .le. 24.656/57.364) then
             nsta = 81
           ELSEIF(br .le. 24.846/57.364) then
             nsta = 82
           ELSEIF(br .le. 25.016/57.364) then
             nsta = 83
           ELSEIF(br .le. 25.326/57.364) then
             nsta = 85
           ELSEIF(br .le. 25.926/57.364) then
             nsta = 86
           ELSEIF(br .le. 26.436/57.364) then
             nsta = 87
           ELSEIF(br .le. 27.216/57.364) then
             nsta = 89
           ELSEIF(br .le. 27.417/57.364) then
             nsta = 90
           ELSEIF(br .le. 27.817/57.364) then
             nsta = 93
           ELSEIF(br .le. 28.041/57.364) then
             nsta = 94
           ELSEIF(br .le. 28.329/57.364) then
             nsta = 95
           ELSEIF(br .le. 29.359/57.364) then
             nsta = 97
           ELSEIF(br .le. 29.539/57.364) then
             nsta = 98
           ELSEIF(br .le. 30.079/57.364) then
             nsta = 99
           ELSEIF(br .le. 30.789/57.364) then
             nsta = 100
           ELSEIF(br .le. 30.939/57.364) then
             nsta = 102
           ELSEIF(br .le. 31.669/57.364) then
             nsta = 103
           ELSEIF(br .le. 32.519/57.364) then
             nsta = 105
           ELSEIF(br .le. 33.029/57.364) then
             nsta = 107
           ELSEIF(br .le. 33.12/57.364) then
             nsta = 108
           ELSEIF(br .le. 33.32/57.364) then
             nsta = 109
           ELSEIF(br .le. 33.35/57.364) then
             nsta = 110
           ELSEIF(br .le. 34.07/57.364) then
             nsta = 112
           ELSEIF(br .le. 34.45/57.364) then
             nsta = 113
           ELSEIF(br .le. 34.79/57.364) then
             nsta = 114
           ELSEIF(br .le. 34.99/57.364) then
             nsta = 115
           ELSEIF(br .le. 35.95/57.364) then
             nsta = 117
           ELSEIF(br .le. 36.21/57.364) then
             nsta = 118
           ELSEIF(br .le. 37.5/57.364) then
             nsta = 119
           ELSEIF(br .le. 37.89/57.364) then
             nsta = 120
           ELSEIF(br .le. 38.79/57.364) then
             nsta = 121
           ELSEIF(br .le. 38.99/57.364) then
             nsta = 122
           ELSEIF(br .le. 39.1/57.364) then
             nsta = 123
           ELSEIF(br .le. 39.24/57.364) then
             nsta = 124
           ELSEIF(br .le. 39.61/57.364) then
             nsta = 125
           ELSEIF(br .le. 39.876/57.364) then
             nsta = 126
           ELSEIF(br .le. 40.026/57.364) then
             nsta = 127
           ELSEIF(br .le. 40.286/57.364) then
             nsta = 128
           ELSEIF(br .le. 40.406/57.364) then
             nsta = 129
           ELSEIF(br .le. 40.976/57.364) then
             nsta = 130
           ELSEIF(br .le. 41.196/57.364) then
             nsta = 131
           ELSEIF(br .le. 41.866/57.364) then
             nsta = 132
           ELSEIF(br .le. 42.526/57.364) then
             nsta = 133
           ELSEIF(br .le. 42.876/57.364) then
             nsta = 134
           ELSEIF(br .le. 43.276/57.364) then
             nsta = 135
           ELSEIF(br .le. 43.856/57.364) then
             nsta = 136
           ELSEIF(br .le. 44.676/57.364) then
             nsta = 137
           ELSEIF(br .le. 45.026/57.364) then
             nsta = 138
           ELSEIF(br .le. 45.226/57.364) then
             nsta = 139
           ELSEIF(br .le. 45.376/57.364) then
             nsta = 140
           ELSEIF(br .le. 45.686/57.364) then
             nsta = 141
           ELSEIF(br .le. 46.206/57.364) then
             nsta = 142
           ELSEIF(br .le. 46.456/57.364) then
             nsta = 143
           ELSEIF(br .le. 46.566/57.364) then
             nsta = 144
           ELSEIF(br .le. 46.626/57.364) then
             nsta = 145
           ELSEIF(br .le. 46.916/57.364) then
             nsta = 146
           ELSEIF(br .le. 46.986/57.364) then
             nsta = 147
           ELSEIF(br .le. 47.396/57.364) then
             nsta = 148
           ELSEIF(br .le. 47.746/57.364) then
             nsta = 149
           ELSEIF(br .le. 48.076/57.364) then
             nsta = 150
           ELSEIF(br .le. 48.756/57.364) then
             nsta = 151
           ELSEIF(br .le. 49.206/57.364) then
             nsta = 152
           ELSEIF(br .le. 49.486/57.364) then
             nsta = 153
           ELSEIF(br .le. 49.636/57.364) then
             nsta = 154
           ELSEIF(br .le. 49.786/57.364) then
             nsta = 155
           ELSEIF(br .le. 49.886/57.364) then
             nsta = 156
           ELSEIF(br .le. 49.996/57.364) then
             nsta = 157
           ELSEIF(br .le. 50.276/57.364) then
             nsta = 158
           ELSEIF(br .le. 50.716/57.364) then
             nsta = 159
           ELSEIF(br .le. 50.836/57.364) then
             nsta = 160
           ELSEIF(br .le. 50.936/57.364) then
             nsta = 161
           ELSEIF(br .le. 51.256/57.364) then
             nsta = 162
           ELSEIF(br .le. 51.501/57.364) then
             nsta = 163
           ELSEIF(br .le. 51.607/57.364) then
             nsta = 164
           ELSEIF(br .le. 51.997/57.364) then
             nsta = 165
           ELSEIF(br .le. 52.377/57.364) then
             nsta = 166
           ELSEIF(br .le. 52.827/57.364) then
             nsta = 167
           ELSEIF(br .le. 52.994/57.364) then
             nsta = 168
           ELSEIF(br .le. 53.164/57.364) then
             nsta = 169
           ELSEIF(br .le. 53.224/57.364) then
             nsta = 170
           ELSEIF(br .le. 53.284/57.364) then
             nsta = 171
           ELSEIF(br .le. 53.374/57.364) then
             nsta = 172
           ELSEIF(br .le. 53.734/57.364) then
             nsta = 173
           ELSEIF(br .le. 53.914/57.364) then
             nsta = 174
           ELSEIF(br .le. 54.174/57.364) then
             nsta = 175
           ELSEIF(br .le. 54.324/57.364) then
             nsta = 176
           ELSEIF(br .le. 54.424/57.364) then
             nsta = 177
           ELSEIF(br .le. 54.534/57.364) then
             nsta = 178
           ELSEIF(br .le. 54.754/57.364) then
             nsta = 179
           ELSEIF(br .le. 55.324/57.364) then
             nsta = 180
           ELSEIF(br .le. 55.554/57.364) then
             nsta = 181
           ELSEIF(br .le. 55.694/57.364) then
             nsta = 182
           ELSEIF(br .le. 55.874/57.364) then
             nsta = 183
           ELSEIF(br .le. 55.934/57.364) then
             nsta = 184
           ELSEIF(br .le. 55.994/57.364) then
             nsta = 185
           ELSEIF(br .le. 56.604/57.364) then
             nsta = 186
           ELSEIF(br .le. 56.724/57.364) then
             nsta = 187
           ELSEIF(br .le. 56.924/57.364) then
             nsta = 188
           ELSE
             nsta = 189
           ENDIF
         Case(77194)
               IF(br .le. 0.171/14.402) then
             nsta = 1
           ELSEIF(br .le. 0.383/14.402) then
             nsta = 2
           ELSEIF(br .le. 1.123/14.402) then
             nsta = 4
           ELSEIF(br .le. 1.453/14.402) then
             nsta = 5
           ELSEIF(br .le. 1.793/14.402) then
             nsta = 8
           ELSEIF(br .le. 2.273/14.402) then
             nsta = 14
           ELSEIF(br .le. 2.703/14.402) then
             nsta = 18
           ELSEIF(br .le. 3.123/14.402) then
             nsta = 20
           ELSEIF(br .le. 4.273/14.402) then
             nsta = 23
           ELSEIF(br .le. 4.566/14.402) then
             nsta = 32
           ELSEIF(br .le. 4.881/14.402) then
             nsta = 33
           ELSEIF(br .le. 5.002/14.402) then
             nsta = 35
           ELSEIF(br .le. 5.582/14.402) then
             nsta = 44
           ELSEIF(br .le. 5.648/14.402) then
             nsta = 45
           ELSEIF(br .le. 6.238/14.402) then
             nsta = 46
           ELSEIF(br .le. 6.298/14.402) then
             nsta = 57
           ELSEIF(br .le. 6.498/14.402) then
             nsta = 62
           ELSEIF(br .le. 6.598/14.402) then
             nsta = 66
           ELSEIF(br .le. 6.698/14.402) then
             nsta = 73
           ELSEIF(br .le. 7.038/14.402) then
             nsta = 78
           ELSEIF(br .le. 7.378/14.402) then
             nsta = 79
           ELSEIF(br .le. 7.738/14.402) then
             nsta = 81
           ELSEIF(br .le. 7.858/14.402) then
             nsta = 86
           ELSEIF(br .le. 8.048/14.402) then
             nsta = 87
           ELSEIF(br .le. 8.568/14.402) then
             nsta = 90
           ELSEIF(br .le. 8.848/14.402) then
             nsta = 91
           ELSEIF(br .le. 9.048/14.402) then
             nsta = 93
           ELSEIF(br .le. 9.258/14.402) then
             nsta = 97
           ELSEIF(br .le. 9.438/14.402) then
             nsta = 108
           ELSEIF(br .le. 9.918/14.402) then
             nsta = 114
           ELSEIF(br .le. 10.328/14.402) then
             nsta = 116
           ELSEIF(br .le. 10.528/14.402) then
             nsta = 118
           ELSEIF(br .le. 10.988/14.402) then
             nsta = 126
           ELSEIF(br .le. 11.268/14.402) then
             nsta = 127
           ELSEIF(br .le. 11.312/14.402) then
             nsta = 128
           ELSEIF(br .le. 11.492/14.402) then
             nsta = 130
           ELSEIF(br .le. 11.822/14.402) then
             nsta = 132
           ELSEIF(br .le. 11.952/14.402) then
             nsta = 134
           ELSEIF(br .le. 12.142/14.402) then
             nsta = 135
           ELSEIF(br .le. 12.352/14.402) then
             nsta = 137
           ELSEIF(br .le. 12.512/14.402) then
             nsta = 138
           ELSEIF(br .le. 12.692/14.402) then
             nsta = 146
           ELSEIF(br .le. 12.822/14.402) then
             nsta = 147
           ELSEIF(br .le. 12.962/14.402) then
             nsta = 148
           ELSEIF(br .le. 13.192/14.402) then
             nsta = 150
           ELSEIF(br .le. 13.542/14.402) then
             nsta = 152
           ELSEIF(br .le. 13.932/14.402) then
             nsta = 153
           ELSEIF(br .le. 14.042/14.402) then
             nsta = 154
           ELSEIF(br .le. 14.182/14.402) then
             nsta = 155
           ELSE
             nsta = 156
           ENDIF
         Case(78195)
               IF(br .le. 0.01/0.0763) then
             nsta = 24
           ELSEIF(br .le. 0.02/0.0763) then
             nsta = 46
           ELSEIF(br .le. 0.0219/0.0763) then
             nsta = 47
           ELSEIF(br .le. 0.0288/0.0763) then
             nsta = 49
           ELSEIF(br .le. 0.0508/0.0763) then
             nsta = 50
           ELSEIF(br .le. 0.0609/0.0763) then
             nsta = 53
           ELSEIF(br .le. 0.0655/0.0763) then
             nsta = 60
           ELSE
             nsta = 66
           ENDIF
         Case(78196)
               IF(br .le. 0.053/2.7903) then
             nsta = 0
           ELSEIF(br .le. 0.0735/2.7903) then
             nsta = 1
           ELSEIF(br .le. 0.1275/2.7903) then
             nsta = 2
           ELSEIF(br .le. 0.1385/2.7903) then
             nsta = 5
           ELSEIF(br .le. 0.1558/2.7903) then
             nsta = 8
           ELSEIF(br .le. 0.1629/2.7903) then
             nsta = 10
           ELSEIF(br .le. 0.172/2.7903) then
             nsta = 15
           ELSEIF(br .le. 0.1782/2.7903) then
             nsta = 17
           ELSEIF(br .le. 0.1903/2.7903) then
             nsta = 21
           ELSEIF(br .le. 0.1983/2.7903) then
             nsta = 24
           ELSEIF(br .le. 0.2114/2.7903) then
             nsta = 25
           ELSEIF(br .le. 0.2404/2.7903) then
             nsta = 27
           ELSEIF(br .le. 0.2584/2.7903) then
             nsta = 28
           ELSEIF(br .le. 0.3674/2.7903) then
             nsta = 30
           ELSEIF(br .le. 0.4404/2.7903) then
             nsta = 33
           ELSEIF(br .le. 0.4614/2.7903) then
             nsta = 34
           ELSEIF(br .le. 0.5474/2.7903) then
             nsta = 36
           ELSEIF(br .le. 0.559/2.7903) then
             nsta = 37
           ELSEIF(br .le. 0.61/2.7903) then
             nsta = 44
           ELSEIF(br .le. 0.657/2.7903) then
             nsta = 46
           ELSEIF(br .le. 0.666/2.7903) then
             nsta = 49
           ELSEIF(br .le. 0.6898/2.7903) then
             nsta = 53
           ELSEIF(br .le. 0.7268/2.7903) then
             nsta = 56
           ELSEIF(br .le. 0.8108/2.7903) then
             nsta = 58
           ELSEIF(br .le. 0.8338/2.7903) then
             nsta = 60
           ELSEIF(br .le. 0.8838/2.7903) then
             nsta = 61
           ELSEIF(br .le. 0.9548/2.7903) then
             nsta = 62
           ELSEIF(br .le. 0.9748/2.7903) then
             nsta = 63
           ELSEIF(br .le. 0.979/2.7903) then
             nsta = 64
           ELSEIF(br .le. 1.001/2.7903) then
             nsta = 67
           ELSEIF(br .le. 1.014/2.7903) then
             nsta = 69
           ELSEIF(br .le. 1.154/2.7903) then
             nsta = 74
           ELSEIF(br .le. 1.186/2.7903) then
             nsta = 76
           ELSEIF(br .le. 1.198/2.7903) then
             nsta = 77
           ELSEIF(br .le. 1.23/2.7903) then
             nsta = 78
           ELSEIF(br .le. 1.247/2.7903) then
             nsta = 79
           ELSEIF(br .le. 1.277/2.7903) then
             nsta = 83
           ELSEIF(br .le. 1.289/2.7903) then
             nsta = 91
           ELSEIF(br .le. 1.367/2.7903) then
             nsta = 93
           ELSEIF(br .le. 1.391/2.7903) then
             nsta = 94
           ELSEIF(br .le. 1.407/2.7903) then
             nsta = 96
           ELSEIF(br .le. 1.52/2.7903) then
             nsta = 98
           ELSEIF(br .le. 1.5312/2.7903) then
             nsta = 100
           ELSEIF(br .le. 1.6492/2.7903) then
             nsta = 108
           ELSEIF(br .le. 1.7462/2.7903) then
             nsta = 113
           ELSEIF(br .le. 2.1562/2.7903) then
             nsta = 114
           ELSEIF(br .le. 2.2412/2.7903) then
             nsta = 120
           ELSEIF(br .le. 2.3772/2.7903) then
             nsta = 121
           ELSEIF(br .le. 2.4772/2.7903) then
             nsta = 128
           ELSEIF(br .le. 2.4962/2.7903) then
             nsta = 130
           ELSEIF(br .le. 2.5122/2.7903) then
             nsta = 130
           ELSEIF(br .le. 2.5252/2.7903) then
             nsta = 131
           ELSEIF(br .le. 2.5391/2.7903) then
             nsta = 132
           ELSEIF(br .le. 2.6081/2.7903) then
             nsta = 132
           ELSEIF(br .le. 2.6671/2.7903) then
             nsta = 133
           ELSEIF(br .le. 2.7231/2.7903) then
             nsta = 133
           ELSEIF(br .le. 2.7581/2.7903) then
             nsta = 134
           ELSEIF(br .le. 2.7841/2.7903) then
             nsta = 135
           ELSE
             nsta = 138
           ENDIF
         Case(78197)
               IF(br .le. 0.0269/0.1579) then
             nsta = 0
           ELSEIF(br .le. 0.0519/0.1579) then
             nsta = 5
           ELSEIF(br .le. 0.0573/0.1579) then
             nsta = 13
           ELSEIF(br .le. 0.0649/0.1579) then
             nsta = 18
           ELSE
             nsta = 20
           ENDIF
         Case(79198)
               IF(br .le. 0.3/12.796) then
             nsta = 32
           ELSEIF(br .le. 0.46/12.796) then
             nsta = 36
           ELSEIF(br .le. 0.7/12.796) then
             nsta = 40
           ELSEIF(br .le. 0.77/12.796) then
             nsta = 41
           ELSEIF(br .le. 0.86/12.796) then
             nsta = 42
           ELSEIF(br .le. 1.41/12.796) then
             nsta = 46
           ELSEIF(br .le. 2.68/12.796) then
             nsta = 48
           ELSEIF(br .le. 2.73/12.796) then
             nsta = 52
           ELSEIF(br .le. 2.79/12.796) then
             nsta = 53
           ELSEIF(br .le. 3.13/12.796) then
             nsta = 54
           ELSEIF(br .le. 3.3/12.796) then
             nsta = 64
           ELSEIF(br .le. 4.1/12.796) then
             nsta = 66
           ELSEIF(br .le. 4.52/12.796) then
             nsta = 68
           ELSEIF(br .le. 4.73/12.796) then
             nsta = 70
           ELSEIF(br .le. 4.8/12.796) then
             nsta = 72
           ELSEIF(br .le. 5.201/12.796) then
             nsta = 84
           ELSEIF(br .le. 5.391/12.796) then
             nsta = 90
           ELSEIF(br .le. 5.915/12.796) then
             nsta = 91
           ELSEIF(br .le. 6.295/12.796) then
             nsta = 92
           ELSEIF(br .le. 6.805/12.796) then
             nsta = 94
           ELSEIF(br .le. 7.225/12.796) then
             nsta = 96
           ELSEIF(br .le. 7.497/12.796) then
             nsta = 97
           ELSEIF(br .le. 7.657/12.796) then
             nsta = 101
           ELSEIF(br .le. 7.991/12.796) then
             nsta = 105
           ELSEIF(br .le. 8.11/12.796) then
             nsta = 106
           ELSEIF(br .le. 8.57/12.796) then
             nsta = 107
           ELSEIF(br .le. 8.965/12.796) then
             nsta = 108
           ELSEIF(br .le. 9.126/12.796) then
             nsta = 112
           ELSEIF(br .le. 9.282/12.796) then
             nsta = 114
           ELSEIF(br .le. 10.152/12.796) then
             nsta = 116
           ELSEIF(br .le. 10.652/12.796) then
             nsta = 118
           ELSEIF(br .le. 10.804/12.796) then
             nsta = 119
           ELSEIF(br .le. 10.863/12.796) then
             nsta = 123
           ELSEIF(br .le. 11.043/12.796) then
             nsta = 124
           ELSEIF(br .le. 11.223/12.796) then
             nsta = 125
           ELSEIF(br .le. 11.336/12.796) then
             nsta = 126
           ELSEIF(br .le. 11.449/12.796) then
             nsta = 128
           ELSEIF(br .le. 11.759/12.796) then
             nsta = 129
           ELSEIF(br .le. 11.849/12.796) then
             nsta = 131
           ELSEIF(br .le. 12.01/12.796) then
             nsta = 132
           ELSEIF(br .le. 12.64/12.796) then
             nsta = 134
           ELSE
             nsta = 135
           ENDIF
         Case(80200)
               IF(br .le. 0.1/264.98) then
             nsta = 0
           ELSEIF(br .le. 0.56/264.98) then
             nsta = 1
           ELSEIF(br .le. 23.66/264.98) then
             nsta = 6
           ELSEIF(br .le. 27.36/264.98) then
             nsta = 10
           ELSEIF(br .le. 31.36/264.98) then
             nsta = 14
           ELSEIF(br .le. 93.86/264.98) then
             nsta = 26
           ELSEIF(br .le. 94.86/264.98) then
             nsta = 36
           ELSEIF(br .le. 96.64/264.98) then
             nsta = 42
           ELSEIF(br .le. 124.14/264.98) then
             nsta = 46
           ELSEIF(br .le. 125.15/264.98) then
             nsta = 51
           ELSEIF(br .le. 142.65/264.98) then
             nsta = 59
           ELSEIF(br .le. 144.39/264.98) then
             nsta = 72
           ELSEIF(br .le. 144.87/264.98) then
             nsta = 73
           ELSEIF(br .le. 164.87/264.98) then
             nsta = 76
           ELSEIF(br .le. 170.09/264.98) then
             nsta = 80
           ELSEIF(br .le. 174.1/264.98) then
             nsta = 82
           ELSEIF(br .le. 194.1/264.98) then
             nsta = 89
           ELSEIF(br .le. 197.8/264.98) then
             nsta = 92
           ELSEIF(br .le. 199.53/264.98) then
             nsta = 93
           ELSEIF(br .le. 211.93/264.98) then
             nsta = 95
           ELSEIF(br .le. 242.03/264.98) then
             nsta = 96
           ELSEIF(br .le. 255.03/264.98) then
             nsta = 99
           ELSEIF(br .le. 259.26/264.98) then
             nsta = 106
           ELSEIF(br .le. 260.24/264.98) then
             nsta = 109
           ELSEIF(br .le. 261.28/264.98) then
             nsta = 113
           ELSE
             nsta = 118
           ENDIF
         Case(80202)
               IF(br .le. 0.067/0.2855) then
             nsta = 0
           ELSEIF(br .le. 0.089/0.2855) then
             nsta = 1
           ELSEIF(br .le. 0.11/0.2855) then
             nsta = 2
           ELSEIF(br .le. 0.123/0.2855) then
             nsta = 4
           ELSEIF(br .le. 0.1307/0.2855) then
             nsta = 6
           ELSEIF(br .le. 0.1321/0.2855) then
             nsta = 7
           ELSEIF(br .le. 0.1451/0.2855) then
             nsta = 11
           ELSEIF(br .le. 0.1504/0.2855) then
             nsta = 15
           ELSEIF(br .le. 0.1577/0.2855) then
             nsta = 18
           ELSEIF(br .le. 0.1917/0.2855) then
             nsta = 20
           ELSEIF(br .le. 0.2257/0.2855) then
             nsta = 22
           ELSEIF(br .le. 0.2355/0.2855) then
             nsta = 27
           ELSEIF(br .le. 0.2485/0.2855) then
             nsta = 33
           ELSE
             nsta = 39
           ENDIF
         Case(81204 ) 
             nsta = 124
         !Case(81204)
         !      IF(br .le. 0.0104/3.1229) then
         !    nsta = 0
         !  ELSEIF(br .le. 0.1394/3.1229) then
         !    nsta = 1
         !  ELSEIF(br .le. 0.1639/3.1229) then
         !    nsta = 4
         !  ELSEIF(br .le. 0.2289/3.1229) then
         !    nsta = 10
         !  ELSEIF(br .le. 0.3099/3.1229) then
         !    nsta = 11
         !  ELSEIF(br .le. 0.4759/3.1229) then
         !    nsta = 12
         !  ELSEIF(br .le. 0.4991/3.1229) then
         !    nsta = 13
         !  ELSEIF(br .le. 0.5213/3.1229) then
         !    nsta = 15
         !  ELSEIF(br .le. 0.6053/3.1229) then
         !    nsta = 18
         !  ELSEIF(br .le. 0.612/3.1229) then
         !    nsta = 19
         !  ELSEIF(br .le. 0.928/3.1229) then
         !    nsta = 24
         !  ELSEIF(br .le. 1.21/3.1229) then
         !    nsta = 26
         !  ELSEIF(br .le. 1.341/3.1229) then
         !    nsta = 30
         !  ELSEIF(br .le. 1.3593/3.1229) then
         !    nsta = 31
         !  ELSEIF(br .le. 1.4383/3.1229) then
         !    nsta = 33
         !  ELSEIF(br .le. 1.5853/3.1229) then
         !    nsta = 34
         !  ELSEIF(br .le. 1.7923/3.1229) then
         !    nsta = 36
         !  ELSEIF(br .le. 1.8763/3.1229) then
         !    nsta = 37
         !  ELSEIF(br .le. 1.8919/3.1229) then
         !    nsta = 40
         !  ELSEIF(br .le. 2.0329/3.1229) then
         !    nsta = 44
         !  ELSEIF(br .le. 2.0909/3.1229) then
         !    nsta = 48
         !  ELSEIF(br .le. 2.1489/3.1229) then
         !    nsta = 51
         !  ELSEIF(br .le. 2.1849/3.1229) then
         !    nsta = 54
         !  ELSEIF(br .le. 2.3489/3.1229) then
         !    nsta = 59
         !  ELSEIF(br .le. 2.3563/3.1229) then
         !    nsta = 60
         !  ELSEIF(br .le. 2.4463/3.1229) then
         !    nsta = 61
         !  ELSEIF(br .le. 2.4601/3.1229) then
         !    nsta = 64
         !  ELSEIF(br .le. 2.6081/3.1229) then
         !    nsta = 65
         !  ELSEIF(br .le. 2.6212/3.1229) then
         !    nsta = 66
         !  ELSEIF(br .le. 2.6792/3.1229) then
         !    nsta = 67
         !  ELSEIF(br .le. 2.7772/3.1229) then
         !    nsta = 68
         !  ELSEIF(br .le. 2.8064/3.1229) then
         !    nsta = 72
         !  ELSEIF(br .le. 2.8244/3.1229) then
         !    nsta = 73
         !  ELSEIF(br .le. 2.8657/3.1229) then
         !    nsta = 74
         !  ELSEIF(br .le. 2.9087/3.1229) then
         !    nsta = 78
         !  ELSEIF(br .le. 2.9181/3.1229) then
         !    nsta = 81
         !  ELSEIF(br .le. 2.9389/3.1229) then
         !    nsta = 83
         !  ELSEIF(br .le. 2.9729/3.1229) then
         !    nsta = 86
         !  ELSEIF(br .le. 2.9939/3.1229) then
         !    nsta = 90
         !  ELSEIF(br .le. 2.9996/3.1229) then
         !    nsta = 91
         !  ELSEIF(br .le. 3.0026/3.1229) then
         !    nsta = 92
         !  ELSEIF(br .le. 3.0476/3.1229) then
         !    nsta = 93
         !  ELSEIF(br .le. 3.0849/3.1229) then
         !    nsta = 94
         !  ELSEIF(br .le. 3.1071/3.1229) then
         !    nsta = 98
         !  ELSEIF(br .le. 3.1157/3.1229) then
         !    nsta = 99
         !  ELSE
         !    nsta = 100
         !  ENDIF
         Case(81206)
               IF(br .le. 0.004/0.0312) then
             nsta = 0
           ELSEIF(br .le. 0.0149/0.0312) then
             nsta = 2
           ELSEIF(br .le. 0.024/0.0312) then
             nsta = 3
           ELSE
             nsta = 4
           ENDIF
         Case(82205)
               IF(br .le. 0.0032/0.005908) then
             nsta = 1
           ELSEIF(br .le. 0.00342/0.005908) then
             nsta = 6
           ELSEIF(br .le. 0.003828/0.005908) then
             nsta = 12
           ELSEIF(br .le. 0.003899/0.005908) then
             nsta = 18
           ELSEIF(br .le. 0.004509/0.005908) then
             nsta = 26
           ELSEIF(br .le. 0.004652/0.005908) then
             nsta = 31
           ELSEIF(br .le. 0.004872/0.005908) then
             nsta = 36
           ELSEIF(br .le. 0.005552/0.005908) then
             nsta = 37
           ELSEIF(br .le. 0.005705/0.005908) then
             nsta = 41
           ELSEIF(br .le. 0.00583/0.005908) then
             nsta = 44
           ELSE
             nsta = 53
           ENDIF
         Case(82207)
               IF(br .le. 0.00691/0.006957) then
             nsta = 0
           ELSE
             nsta = 2
           ENDIF
         Case(82208)
               IF(br .le. 0.137/0.137455) then
             nsta = 0
           ELSEIF(br .le. 0.137233/0.137455) then
             nsta = 1
           ELSEIF(br .le. 0.137318/0.137455) then
             nsta = 11
           ELSEIF(br .le. 0.137354/0.137455) then
             nsta = 17
           ELSEIF(br .le. 0.137384/0.137455) then
             nsta = 43
           ELSEIF(br .le. 0.137422/0.137455) then
             nsta = 94
           ELSEIF(br .le. 0.137433/0.137455) then
             nsta = 151
           ELSE
             nsta = 204
           ENDIF
         Case(83210)
               IF(br .le. 0.00042/0.061334) then
             nsta = 3
           ELSEIF(br .le. 0.00282/0.061334) then
             nsta = 4
           ELSEIF(br .le. 0.01992/0.061334) then
             nsta = 5
           ELSEIF(br .le. 0.02165/0.061334) then
             nsta = 6
           ELSEIF(br .le. 0.03055/0.061334) then
             nsta = 7
           ELSEIF(br .le. 0.04425/0.061334) then
             nsta = 8
           ELSEIF(br .le. 0.04561/0.061334) then
             nsta = 13
           ELSEIF(br .le. 0.04582/0.061334) then
             nsta = 14
           ELSEIF(br .le. 0.04607/0.061334) then
             nsta = 18
           ELSEIF(br .le. 0.04777/0.061334) then
             nsta = 19
           ELSEIF(br .le. 0.04944/0.061334) then
             nsta = 20
           ELSEIF(br .le. 0.04966/0.061334) then
             nsta = 23
           ELSEIF(br .le. 0.049707/0.061334) then
             nsta = 25
           ELSEIF(br .le. 0.049917/0.061334) then
             nsta = 26
           ELSEIF(br .le. 0.050527/0.061334) then
             nsta = 28
           ELSEIF(br .le. 0.050937/0.061334) then
             nsta = 29
           ELSEIF(br .le. 0.052387/0.061334) then
             nsta = 33
           ELSEIF(br .le. 0.053187/0.061334) then
             nsta = 36
           ELSEIF(br .le. 0.054977/0.061334) then
             nsta = 38
           ELSEIF(br .le. 0.055067/0.061334) then
             nsta = 44
           ELSEIF(br .le. 0.056607/0.061334) then
             nsta = 47
           ELSEIF(br .le. 0.056877/0.061334) then
             nsta = 50
           ELSEIF(br .le. 0.058537/0.061334) then
             nsta = 52
           ELSEIF(br .le. 0.058584/0.061334) then
             nsta = 53
           ELSEIF(br .le. 0.058894/0.061334) then
             nsta = 55
           ELSEIF(br .le. 0.060994/0.061334) then
             nsta = 58
           ELSE
             nsta = 60
           ENDIF
         Case(90233)
               IF(br .le. 0.0037/0.9092) then
             nsta = 0
           ELSEIF(br .le. 0.0084/0.9092) then
             nsta = 3
           ELSEIF(br .le. 0.0101/0.9092) then
             nsta = 13
           ELSEIF(br .le. 0.0144/0.9092) then
             nsta = 19
           ELSEIF(br .le. 0.0237/0.9092) then
             nsta = 34
           ELSEIF(br .le. 0.027/0.9092) then
             nsta = 35
           ELSEIF(br .le. 0.038/0.9092) then
             nsta = 36
           ELSEIF(br .le. 0.0397/0.9092) then
             nsta = 39
           ELSEIF(br .le. 0.0415/0.9092) then
             nsta = 43
           ELSEIF(br .le. 0.0475/0.9092) then
             nsta = 48
           ELSEIF(br .le. 0.0593/0.9092) then
             nsta = 51
           ELSEIF(br .le. 0.063/0.9092) then
             nsta = 56
           ELSEIF(br .le. 0.0671/0.9092) then
             nsta = 61
           ELSEIF(br .le. 0.0939/0.9092) then
             nsta = 63
           ELSEIF(br .le. 0.0954/0.9092) then
             nsta = 66
           ELSEIF(br .le. 0.0977/0.9092) then
             nsta = 72
           ELSEIF(br .le. 0.1034/0.9092) then
             nsta = 74
           ELSEIF(br .le. 0.1059/0.9092) then
             nsta = 75
           ELSEIF(br .le. 0.113/0.9092) then
             nsta = 80
           ELSEIF(br .le. 0.1228/0.9092) then
             nsta = 84
           ELSEIF(br .le. 0.13/0.9092) then
             nsta = 85
           ELSEIF(br .le. 0.1415/0.9092) then
             nsta = 88
           ELSEIF(br .le. 0.1501/0.9092) then
             nsta = 89
           ELSEIF(br .le. 0.155/0.9092) then
             nsta = 90
           ELSEIF(br .le. 0.1584/0.9092) then
             nsta = 95
           ELSEIF(br .le. 0.1657/0.9092) then
             nsta = 96
           ELSEIF(br .le. 0.1714/0.9092) then
             nsta = 98
           ELSEIF(br .le. 0.1833/0.9092) then
             nsta = 99
           ELSEIF(br .le. 0.1871/0.9092) then
             nsta = 104
           ELSEIF(br .le. 0.2268/0.9092) then
             nsta = 105
           ELSEIF(br .le. 0.2388/0.9092) then
             nsta = 107
           ELSEIF(br .le. 0.2558/0.9092) then
             nsta = 108
           ELSEIF(br .le. 0.2607/0.9092) then
             nsta = 109
           ELSEIF(br .le. 0.3177/0.9092) then
             nsta = 112
           ELSEIF(br .le. 0.3246/0.9092) then
             nsta = 113
           ELSEIF(br .le. 0.3479/0.9092) then
             nsta = 115
           ELSEIF(br .le. 0.369/0.9092) then
             nsta = 117
           ELSEIF(br .le. 0.3881/0.9092) then
             nsta = 120
           ELSEIF(br .le. 0.3925/0.9092) then
             nsta = 121
           ELSEIF(br .le. 0.406/0.9092) then
             nsta = 123
           ELSEIF(br .le. 0.4111/0.9092) then
             nsta = 124
           ELSEIF(br .le. 0.4279/0.9092) then
             nsta = 127
           ELSEIF(br .le. 0.4381/0.9092) then
             nsta = 128
           ELSEIF(br .le. 0.4432/0.9092) then
             nsta = 130
           ELSEIF(br .le. 0.4597/0.9092) then
             nsta = 131
           ELSEIF(br .le. 0.466/0.9092) then
             nsta = 132
           ELSEIF(br .le. 0.4716/0.9092) then
             nsta = 134
           ELSEIF(br .le. 0.4746/0.9092) then
             nsta = 136
           ELSEIF(br .le. 0.4869/0.9092) then
             nsta = 137
           ELSEIF(br .le. 0.504/0.9092) then
             nsta = 140
           ELSEIF(br .le. 0.5119/0.9092) then
             nsta = 141
           ELSEIF(br .le. 0.5208/0.9092) then
             nsta = 142
           ELSEIF(br .le. 0.5416/0.9092) then
             nsta = 143
           ELSEIF(br .le. 0.5503/0.9092) then
             nsta = 144
           ELSEIF(br .le. 0.5561/0.9092) then
             nsta = 145
           ELSEIF(br .le. 0.5601/0.9092) then
             nsta = 146
           ELSEIF(br .le. 0.5687/0.9092) then
             nsta = 148
           ELSEIF(br .le. 0.5726/0.9092) then
             nsta = 149
           ELSEIF(br .le. 0.581/0.9092) then
             nsta = 152
           ELSEIF(br .le. 0.5841/0.9092) then
             nsta = 154
           ELSEIF(br .le. 0.5892/0.9092) then
             nsta = 155
           ELSEIF(br .le. 0.5958/0.9092) then
             nsta = 157
           ELSEIF(br .le. 0.6042/0.9092) then
             nsta = 158
           ELSEIF(br .le. 0.6106/0.9092) then
             nsta = 159
           ELSEIF(br .le. 0.6139/0.9092) then
             nsta = 160
           ELSEIF(br .le. 0.6221/0.9092) then
             nsta = 163
           ELSEIF(br .le. 0.6246/0.9092) then
             nsta = 164
           ELSEIF(br .le. 0.6339/0.9092) then
             nsta = 165
           ELSEIF(br .le. 0.6416/0.9092) then
             nsta = 166
           ELSEIF(br .le. 0.6475/0.9092) then
             nsta = 168
           ELSEIF(br .le. 0.6619/0.9092) then
             nsta = 169
           ELSEIF(br .le. 0.6729/0.9092) then
             nsta = 170
           ELSEIF(br .le. 0.6839/0.9092) then
             nsta = 172
           ELSEIF(br .le. 0.6914/0.9092) then
             nsta = 175
           ELSEIF(br .le. 0.6944/0.9092) then
             nsta = 176
           ELSEIF(br .le. 0.7007/0.9092) then
             nsta = 177
           ELSEIF(br .le. 0.7088/0.9092) then
             nsta = 178
           ELSEIF(br .le. 0.716/0.9092) then
             nsta = 179
           ELSEIF(br .le. 0.724/0.9092) then
             nsta = 180
           ELSEIF(br .le. 0.74/0.9092) then
             nsta = 181
           ELSEIF(br .le. 0.753/0.9092) then
             nsta = 183
           ELSEIF(br .le. 0.767/0.9092) then
             nsta = 184
           ELSEIF(br .le. 0.775/0.9092) then
             nsta = 185
           ELSEIF(br .le. 0.7835/0.9092) then
             nsta = 187
           ELSEIF(br .le. 0.7965/0.9092) then
             nsta = 188
           ELSEIF(br .le. 0.8065/0.9092) then
             nsta = 189
           ELSEIF(br .le. 0.8175/0.9092) then
             nsta = 190
           ELSEIF(br .le. 0.8246/0.9092) then
             nsta = 192
           ELSEIF(br .le. 0.8364/0.9092) then
             nsta = 196
           ELSEIF(br .le. 0.8433/0.9092) then
             nsta = 197
           ELSEIF(br .le. 0.8502/0.9092) then
             nsta = 198
           ELSEIF(br .le. 0.8602/0.9092) then
             nsta = 199
           ELSEIF(br .le. 0.8678/0.9092) then
             nsta = 200
           ELSEIF(br .le. 0.8808/0.9092) then
             nsta = 201
           ELSEIF(br .le. 0.8895/0.9092) then
             nsta = 202
           ELSEIF(br .le. 0.9002/0.9092) then
             nsta = 203
           ELSE
             nsta = 204
           ENDIF
         Case(92239)
               IF(br .le. 0.00007/0.4877) then
             nsta = 0
           ELSEIF(br .le. 0.00123/0.4877) then
             nsta = 4
           ELSEIF(br .le. 0.00463/0.4877) then
             nsta = 5
           ELSEIF(br .le. 0.00773/0.4877) then
             nsta = 7
           ELSEIF(br .le. 0.00921/0.4877) then
             nsta = 18
           ELSEIF(br .le. 0.01067/0.4877) then
             nsta = 21
           ELSEIF(br .le. 0.01189/0.4877) then
             nsta = 22
           ELSEIF(br .le. 0.01203/0.4877) then
             nsta = 23
           ELSEIF(br .le. 0.01933/0.4877) then
             nsta = 24
           ELSEIF(br .le. 0.20533/0.4877) then
             nsta = 25
           ELSEIF(br .le. 0.20563/0.4877) then
             nsta = 26
           ELSEIF(br .le. 0.22973/0.4877) then
             nsta = 29
           ELSEIF(br .le. 0.25563/0.4877) then
             nsta = 30
           ELSEIF(br .le. 0.25614/0.4877) then
             nsta = 32
           ELSEIF(br .le. 0.25634/0.4877) then
             nsta = 35
           ELSEIF(br .le. 0.25642/0.4877) then
             nsta = 36
           ELSEIF(br .le. 0.25906/0.4877) then
             nsta = 38
           ELSEIF(br .le. 0.26586/0.4877) then
             nsta = 40
           ELSEIF(br .le. 0.26785/0.4877) then
             nsta = 41
           ELSEIF(br .le. 0.26911/0.4877) then
             nsta = 43
           ELSEIF(br .le. 0.27077/0.4877) then
             nsta = 44
           ELSEIF(br .le. 0.27095/0.4877) then
             nsta = 46
           ELSEIF(br .le. 0.27109/0.4877) then
             nsta = 48
           ELSEIF(br .le. 0.27489/0.4877) then
             nsta = 49
           ELSEIF(br .le. 0.27633/0.4877) then
             nsta = 52
           ELSEIF(br .le. 0.27833/0.4877) then
             nsta = 53
           ELSEIF(br .le. 0.28523/0.4877) then
             nsta = 54
           ELSEIF(br .le. 0.29743/0.4877) then
             nsta = 55
           ELSEIF(br .le. 0.31203/0.4877) then
             nsta = 56
           ELSEIF(br .le. 0.35403/0.4877) then
             nsta = 59
           ELSEIF(br .le. 0.35472/0.4877) then
             nsta = 61
           ELSEIF(br .le. 0.35586/0.4877) then
             nsta = 62
           ELSEIF(br .le. 0.35627/0.4877) then
             nsta = 63
           ELSEIF(br .le. 0.36047/0.4877) then
             nsta = 64
           ELSEIF(br .le. 0.3614/0.4877) then
             nsta = 65
           ELSEIF(br .le. 0.36193/0.4877) then
             nsta = 66
           ELSEIF(br .le. 0.36258/0.4877) then
             nsta = 67
           ELSEIF(br .le. 0.36422/0.4877) then
             nsta = 68
           ELSEIF(br .le. 0.36507/0.4877) then
             nsta = 70
           ELSEIF(br .le. 0.36781/0.4877) then
             nsta = 71
           ELSEIF(br .le. 0.3686/0.4877) then
             nsta = 72
           ELSEIF(br .le. 0.36943/0.4877) then
             nsta = 73
           ELSEIF(br .le. 0.37393/0.4877) then
             nsta = 75
           ELSEIF(br .le. 0.37624/0.4877) then
             nsta = 76
           ELSEIF(br .le. 0.37646/0.4877) then
             nsta = 77
           ELSEIF(br .le. 0.3765/0.4877) then
             nsta = 78
           ELSEIF(br .le. 0.37898/0.4877) then
             nsta = 79
           ELSEIF(br .le. 0.37912/0.4877) then
             nsta = 80
           ELSEIF(br .le. 0.38003/0.4877) then
             nsta = 81
           ELSEIF(br .le. 0.38157/0.4877) then
             nsta = 82
           ELSEIF(br .le. 0.38198/0.4877) then
             nsta = 83
           ELSEIF(br .le. 0.38222/0.4877) then
             nsta = 84
           ELSEIF(br .le. 0.38372/0.4877) then
             nsta = 85
           ELSEIF(br .le. 0.38772/0.4877) then
             nsta = 86
           ELSEIF(br .le. 0.39142/0.4877) then
             nsta = 87
           ELSEIF(br .le. 0.39292/0.4877) then
             nsta = 88
           ELSEIF(br .le. 0.39992/0.4877) then
             nsta = 89
           ELSEIF(br .le. 0.4002/0.4877) then
             nsta = 90
           ELSEIF(br .le. 0.4042/0.4877) then
             nsta = 91
           ELSEIF(br .le. 0.4142/0.4877) then
             nsta = 92
           ELSEIF(br .le. 0.4262/0.4877) then
             nsta = 93
           ELSEIF(br .le. 0.4422/0.4877) then
             nsta = 94
           ELSEIF(br .le. 0.4469/0.4877) then
             nsta = 95
           ELSEIF(br .le. 0.4536/0.4877) then
             nsta = 96
           ELSEIF(br .le. 0.4616/0.4877) then
             nsta = 99
           ELSEIF(br .le. 0.4686/0.4877) then
             nsta = 100
           ELSEIF(br .le. 0.4757/0.4877) then
             nsta = 101
           ELSE
             nsta = 104
           ENDIF
      EndSelect
      
      if(nsta .eq. -1) return ! 2015/7/29 neutron capture scheme unavailable
       
      nph      = nph + 1
      eph(nph) = eexe - elevel(nsta)
      kfejec(nph) = 22 ! photon
      eexe     = elevel(nsta)

      return
      end subroutine

      
      end module levdat
************************************************************************

