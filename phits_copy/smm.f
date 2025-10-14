************************************************************************
*                                                                      *
      module smmmod
*                                                                      *
*                                                                      *
*        Last Revised:     2014 03 29  by T.Ogawa                      *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              Run the statistical multifragmentation model developed  *
*              by Bondorf dedicated to heavy nuclei OR run the Fermi   *
*              breakup model by Fermi for light nuclei                 *
*                                                                      *
*                                                                      *
************************************************************************

      use NGSDATAMOD, only : spin, bindeg
      implicit double precision(a-h, o-z)
      private :: eps
      private :: vf
      private :: xi
      private :: tlam
      private :: sint
      private :: fact

      logical,private:: smmsuc
      integer,private:: multf
      integer,private:: ibar(500),icha(500),nchart(0:260,0:100) ! T.Sato 2017/1/11, 250 -> 260
      integer,private:: nchf(0:18,0:18), nchfd(0:18,0:18)
      integer,private:: ibarec(500,10000),ichrec(500,10000)
      integer,private:: nlcrec(10000),mulrec(10000)
      double precision,private:: etran(500),trec(10000),wrec(10000)
      double precision,private:: ekinfbm(10000) ! used only in fbmexec
!$OMP THREADPRIVATE(ibar,icha,nchart,ibarec,ichrec,nlcrec,mulrec)
!$OMP THREADPRIVATE(nchf, nchfd)
!$OMP THREADPRIVATE(multf)
!$OMP THREADPRIVATE(etran,trec,wrec,smmsuc)

      double precision, allocatable :: exmfrg(:)
      double precision, allocatable :: pxmfrg(:)
      double precision, allocatable :: pymfrg(:)
      double precision, allocatable :: pzmfrg(:)
      double precision, allocatable :: frgrm (:)
      integer, allocatable :: ibarf(:)
      integer, allocatable :: ichaf(:)
!$OMP THREADPRIVATE(exmfrg,pxmfrg,pymfrg,pzmfrg,frgrm,ibarf,ichaf)


* to change maxp, increase size of **rec vectors and matrice
      integer,private,parameter:: maxp=10000
      double precision,private,parameter::  tc=16.0d0
      double precision,private,parameter::  vp=8.85418782d-12
      double precision,private,parameter::  ec=1.60217646d-19
      double precision,private,parameter:: temmin=1.d-1
      double precision,private,parameter:: r0=1.2d-15
      double precision,private,parameter:: r0f=1.4d-15
      double precision,private,parameter:: syec=28.1d0
      double precision,private,parameter:: benm=-15.68d0
      double precision,private,parameter:: srec=18.56d0
      double precision,private,parameter:: sol=2.99792458d8
      double precision,private,parameter:: hbar =197.3269718
      double precision,private,parameter:: rpmas=9.38272046d2
      double precision,private,parameter:: rnmas=9.39565379d2

      logical,public,save:: lfirst=.false.

      integer,private:: maxa(0:100)
      integer,private:: mina(0:100)
      data mina /1, 1, 3, 3, 5, 6, 8, 10, 12, 14, 16, 18, 19, 21, 22,
     &         24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 45, 47, 48,
     &         52, 54, 56, 58, 60, 65, 67, 69, 71, 73, 76, 78, 81, 83,
     &         85, 87, 89, 91, 93, 95, 97, 99, 103, 105, 108, 110, 112,
     &         114, 117, 119, 121, 124, 126, 128, 130, 134, 136, 138,
     &         140, 143, 145, 148, 150, 153, 155, 158, 160, 162, 164,
     &         166, 169, 171, 176, 178, 184, 188, 193, 195, 199, 202,
     &         206, 209, 212, 217, 225, 228, 231, 233, 235, 237, 240,
     &         242/
      data maxa /1, 7, 10, 12, 16, 19, 22, 25, 28, 31, 34, 37, 40, 42,
     &         44, 46, 49, 51, 53, 55, 57, 60, 63, 65, 67, 69, 72, 75,
     &         78, 80, 83, 86, 89, 92, 94, 97, 100, 102, 105, 108, 110,
     &         113, 115, 118, 120, 122, 124, 130, 132, 135, 137, 139,
     &         142, 144, 147, 151, 153, 155, 157, 159, 161, 163, 165,
     &         167, 169, 171, 173, 175, 177, 179, 181, 184, 188, 190,
     &         192, 194, 196, 199, 202, 205, 210, 212, 215, 218, 220,
     &         223, 228, 232, 234, 236, 238, 240, 242, 244, 247, 249,
     &         250, 250, 250, 250, 250/

      contains


************************************************************************
*                                                                      *
         subroutine smmexec(ia0, iz0, px, py, pz, ex, multf )
*                                                                      *
*        Last Revised:     2014 03 29                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to execute statistical multi-fragmentation              *
*                                                                      *
*        Variables:                                                    *
*                                                                      *
*        ia0, iz0   : mass and charge number of mother (input)         *
*        ex         : excitation energy of mother (MeV)   (input)      *
*        multf      : multiplicity (output)                            *
*        ibarf(i)   : barion number of i-th fragment (i = 1- multf)    *
*        ichaf(i)   : charge of i-th fragment        (i = 1- multf)    *
*        exmfrg(i)  : excitation energy of fragment  (i = 1- multf)    *
*        p*mfrg(i)  : momentum of fragment           (i = 1- multf)    *
*                                                                      *
************************************************************************

      implicit doubleprecision(a-h,o-z)

*-----------------------------------------------------------------------

      dimension enecn0(500), epsi(500)

      smmsuc = .false.


*-----------------------------------------------------------------------

*  Loop for different partition
      do 100 mp = 1, maxp

*---------------------------------------------------------------
      wrec(mp) = 0.d0
      trec(mp) = 0.d0
      mulrec(mp) = 0
      nlcrec(mp) = 0

      If(mp .eq. 1) then
       ibar(1) = ia0
       icha(1) = iz0
       mult    = 1
       goto 110
      ENDIF

*  Partition of nucleons (i.e. determine mass of fragments)

      jcont = 1
      ipart = 0

* Sample a large fragment once per 5 trials
       IF( mod(mp, 5) .eq. 0) then
        ibar(1)  = idnint( dble(ia0) * (0.3d0 + 0.5d0 * unirn(dummy) ) )
        ipart    = ibar(1)
        jcont    = 2
       ENDIF

* Expected multiplicity
       mulmea = idint( dint( dmax1(
     & 2.0d0**(2.d0*unirn(dummy)-1.d0) * ( (7.d0 / 7.5d4 * dble(ia0)
     & **2 + 1.d-2 * dble(ia0) + 4.d0 / 15.d0 ) * ex / dble(ia0)
     & + (-1.d0/3.d3 * dble(ia0) **2) - 3.d0 / 1.d2 * dble(ia0)
     & + 1.d0 / 3.d0 ) ,
     & + 1.d0 + 4.d0 * unirn(dummy)  ))+1.d0)

* Nucleon sampling routine #1
       IF( mod(mp, 2) .eq. 0 ) then
        do 200 icont = 1 + ipart, ia0 -1
         IF( unirn(dummy) .lt. (dble(mulmea) -1.d0)/ dble(ia0) ) then
          ibar(jcont) = icont - ipart
          ipart = icont
          jcont = jcont + 1
         ENDIF
 200    continue
       mult   = jcont
       ibar(mult)= ia0 - ipart

* Nucleon sampling routine #2
       ELSE
        mult = mulmea
        do 205 icont = jcont, mult - 1
         ibar(icont) = int( dint( dble( ia0 - ipart) * unirn(dummy)) )
         ipart = ipart + ibar(icont)
 205    continue
       ibar(mult)= ia0 - ipart
       ENDIF

       IF(jcont .eq. 1 .or. ibar(mult) .le. 0)  goto 99
*---------------------------------------------------------------

*---------------------------------------------------------------
* Partition of protons. Giving protons to partitions at random
        ic0    = iz0
        ib0    = ia0
       do 210 jcont = 1, mult-1
        kcont = 0
        do 220 icont = 1, ibar(jcont)
         IF( dble(ic0 - kcont) / dble(max0(1,ib0-icont))
     &                                .gt. unirn(dummy) ) then
          kcont = kcont + 1
         ENDIF
 220    continue
         icha(jcont) = kcont
         ib0 = ib0 - ibar(jcont)
         ic0 = ic0 - icha(jcont)
 210   continue
       icha(mult)= ic0
*---------------------------------------------------------------

*---------------------------------------------------------------
* Destroy neutronium
       do 230 jcont = 1, mult
        IF(icha(jcont) .eq. 0 .and. ibar(jcont) .gt. 1) then
         do 240 icont = mult + 1, mult + ibar(jcont) - 1
          icha(icont) = 0
          ibar(icont) = 1
 240     continue
         mult        = mult + ibar(jcont) - 1
         icha(jcont) = 0
         ibar(jcont) = 1
        ENDIF
 230   continue

* Destroy proton matter
       do 250 jcont = 1, mult
        IF(icha(jcont) .eq. ibar(jcont) .and. ibar(jcont) .gt. 1) then
         do 260 icont = mult + 1, mult + ibar(jcont) - 1
          icha(icont) = 1
          ibar(icont) = 1
 260     continue
         mult        = mult + ibar(jcont) - 1
         icha(jcont) = 1
         ibar(jcont) = 1
        ENDIF
 250   continue
*---------------------------------------------------------------

*---------------------------------------------------------------
* Label #110 is for "no partitioning" check

* chekck if nuclei are in "normal" zone in nuclear chart
 110   do 300 jcont = 1, mult
        IF(ibar(jcont) .gt. maxa(icha(jcont)) .or.
     &     ibar(jcont) .lt. mina(icha(jcont))) then
         goto 99
        ENDIF
 300    continue
        IF(icha(mult) .lt. 0)  goto 99
*---------------------------------------------------------------

*---------------------------------------------------------------
*------------- isotopic multiplicity----------------------------
*  initialize dimension
        do 309 jcont = 0 , 100
        do 308 icont = 0 , 250
          nchart(icont, jcont) = 0
 308    continue
 309    continue
*  count fragments
        do 310 jcont = 1, mult
         nchart(ibar(jcont),icha(jcont))
     &     = min(30, nchart(ibar(jcont),icha(jcont)) + 1) ! 2017/11/13 ogawa. Avoid Overflow

 310    continue



*---------------------------------------------------------------
      do jcont = 1, mult
        epsi(jcont) = 1.d0 / eps(ibar(jcont))
      end do

*---------------------------------------------------------------
*--------temperature determination routine----------------------
      do 320 levt = 1, idint(dint(tc/5.d-2))

      tem = dble(levt) * 5.d-2

* Energy of compound nucleus at temperature tem
      enecn = 6.d-1 * dble(iz0) **2 * ec / (r0 * dble(ia0)
     &   ** 0.33333d0 * (1.d0 + xi(mult,ia0) )**0.33333d0 )
     &  / (4.d0 * 3.141d0 * vp ) * 1.d-6 + 1.5d0 * tem * dble(mult)

*---------------------------------------------------------------
* Sum-up internal energy of fragmented system

      do jcont = 1, mult
        enecn0(jcont) = (benm + tem**2 * epsi(jcont))
     &   * dble(ibar(jcont)) + syec
     &   * dble( (ibar(jcont) - 2 * icha(jcont))**2) / dble(ibar(jcont))
     &   + (srec * ((tc**2 - tem**2) / (tc**2 + tem**2))
     &   ** 1.25d0 + srec * 5.d0 * tem**2 * tc**2
     &   * ((tc**2 - tem**2) / (tc**2 + tem**2))
     &   ** 0.25d0 / (tc**2 + tem**2) **2 )
     &   * dble(ibar(jcont)) ** 0.666666d0
     &   + 0.6d0 * dble(icha(jcont) **2) * ec
     &   / (r0 * dble(ibar(jcont)) ** 0.3333d0)
     &   / (4.d0 * 3.141d0 * vp)
     &   * (1.d0 - 1.d0 / (1.d0 + xi(mult,ia0))** 0.333333d0 ) * 1.d-6
      enddo
      do jcont = 1, mult
       IF(ibar(jcont) .eq. 1 ) then
        enecn0(jcont) = 0.0d0
       ELSEIF (ibar(jcont) .le. 3) then
        enecn0(jcont) = -bindeg(icha(jcont),ibar(jcont)-icha(jcont))
     & * dble(ibar(jcont))
       ELSEIF (ibar(jcont) .eq. 4 .and. icha(jcont) .eq. 2 ) then
        enecn0(jcont) = -bindeg(icha(jcont),ibar(jcont)-icha(jcont))
     & * dble(ibar(jcont)) + tem**2 / epsi(jcont) * dble(ibar(jcont))
       ENDIF
      enddo

      enecn = enecn + sum(enecn0(1:mult))

*---------------------------------------------------------------

* Enecn mononeneously increases with temperature
       IF ( enecn .gt. (ex - bindeg(iz0,ia0-iz0) * dble(ia0)) ) goto 340
 320  continue

 340  tem = ( dble(levt) - 5.d-1 ) * 5.d-2

*---------------------------------------------------------------
*--------temperature limitation---------------------------------
      IF(levt .eq. 1 .or. levt .ge. idint(dint(tc/5.d-2))) goto 99
*---------------------------------------------------------------

*---------------------------------------------------------------
*--------remember this partition--------------------------

       smmsuc   = .true.

       IF(mult .gt. 1) then
        w =1.d0
       do 400 jcont=1, mult
        w = w * ( ( (2.d0 * spin(19,icha(jcont)*1000000+ibar(jcont))
     &   + 1.d0)
     &   * vf(ia0, mult) / tlam(tem)**3 * dble(ibar(jcont)) ** 1.5d0
     &   * exp(1.5d0 + sint(tem,ibar(jcont),icha(jcont),ia0,iz0)) )
     &   ** nchart(ibar(jcont),icha(jcont))
     &   /  fact( nchart(ibar(jcont),icha(jcont)))  )
     &   ** (1.d0/ dble( nchart(ibar(jcont),icha(jcont)) ) )
 400   continue
       ELSE
        w = exp( sint (tem, ia0, iz0, ia0, iz0 ) )
       ENDIF

       do 450 icont=1, mult
        ibarec(icont, mp) = ibar(icont)
        ichrec(icont, mp) = icha(icont)

 450   continue

       nlcrec(mp) = nchart(1,0) + nchart(1,1) + nchart(2,1)
     &  + nchart(3,1) + nchart(3,2)
       mulrec(mp) = mult
       trec  (mp) = tem
       wrec  (mp) = w

*---------------------------------------------------------------

 99   do 600 icont = 1, mult
* Initialize barion, charge assignment not to affect next partition
       ibar(icont) = 0
       icha(icont) = 0
 600  continue


*---------------------------------------------------------------
 100  continue

      IF( .not. smmsuc) then
           tfin = 0.d0
          multf = 1

         allocate(ibarf(multf), ichaf(multf), exmfrg(multf),
     &    frgrm(multf), pxmfrg(multf), pymfrg(multf), pzmfrg(multf) )

       ibarf(1) = ia0
       ichaf(1) = iz0
       goto 685
      ENDIF

*----------------------------------------------------------------
*--Determine partition-------------------------------------------

       ransel = unirn(dummy)

        wsum = 0.d0
       do 650 icont = 1, maxp
        wsum = wrec(icont) + wsum
 650   continue

       wrecpt = 0.d0
       do icont = 1, maxp
         wrecpt = wrec(icont) / wsum + wrecpt
         IF(wrecpt .ge. ransel) then
           ifinal =  icont
         exit
         ENDIF
       end do

*    Score selected partition
       multf = mulrec(ifinal)
       tfin  = trec(ifinal)

      allocate( ibarf(multf), ichaf(multf), exmfrg(multf),
     &  frgrm(multf), pxmfrg(multf), pymfrg(multf), pzmfrg(multf) )

       do 680 icont = 1, mulrec(ifinal)
        ibarf(icont) = ibarec(icont, ifinal)
        ichaf(icont) = ichrec(icont, ifinal)
 680   continue


*----------------------------------------------------------------

* ----SMM finalization. Determine energy, momentum, rest mass...-


 685   totmas = 0.d0
       do 690 icont = 1, multf
*  Rest mass
        frgrm(icont) = rnmas * dble(ibarf(icont) - ichaf(icont))
     &              + rpmas * dble(ichaf(icont))
     &              - bindeg(ichaf(icont),ibarf(icont)-ichaf(icont))
 690   continue
*
       do 700 icont = 1, multf
*  Excitation energy
       IF(multf .eq. 1) then
        exmfrg(icont) = ex
       ELSEIF(ibarf(icont) .le. 3) then
        exmfrg(icont) = 0.d0
       ELSE
        exmfrg(icont) = tfin **2 / eps(ibarf(icont)) * ibarf(icont)
     &   + (srec * ((tc**2 - tfin **2) / (tc**2 + tfin**2))
     &   ** 1.25d0 + srec * 5.d0 * tfin**2 * tc**2
     &   * ((tc**2 - tfin**2) / (tc**2 + tfin**2)) **0.25d0
     &   / (tc**2 + tfin **2) **2 - srec) * dble(ibar(jcont))
     &   ** 0.666666d0
       ENDIF

 700   continue

*  Translational motion kinetic energy in CMS frame
       IF(multf .eq. 1) then
        etran (1)     = 0.d0
*  Analytical solution of relativistic kinematics for multiplicity = 2
       ELSEIF(multf .eq. 2) then
        etran (1)     = 1.5d0 * tfin * (frgrm(2) + 7.5d-1 * tfin)
     & / ( frgrm(1) + frgrm(2) + 1.5d0 * tfin )
        etran (2)     = 1.5d0 * tfin * (frgrm(1) + 7.5d-1 * tfin)
     & / ( frgrm(1) + frgrm(2) + 1.5d0 * tfin )

       ELSE
* Approximative solution for multiplicity > 2
       do 705 icont = 1, multf
        etran (icont) = - 1.5d0 * tfin * dlog(1.d0 - unirn(dummy))
 705   continue
       ENDIF


*---------------------------------------------------------------
*-----------   Correction for energy conservation --------------
       IF(multf .gt. 1) then
         etot = 0
        do 710 icont = 1, multf
         etot = etot + etran (icont)
 710    continue
*   Correction factor to conserve energy
         etocof = 1.5d0 * dble(multf) * tfin / etot
        do 720 icont = 1, multf
         etran(icont) = etran(icont) * etocof
 720    continue
       ENDIF
*---------------------------------------------------------------


*---------------------------------------------------------------
* ----------   Correction for momentum conservation ------------
       IF(multf .eq. 1) then
        pxmfrg(1) = px * 1.d3
        pymfrg(1) = py * 1.d3
        pzmfrg(1) = pz * 1.d3
       ELSEIF(multf .eq. 2) then
*    Sample angular cosines
        randx  = 2.d0 * unirn(dummy) -1.d0
        randy  = 2.d0 * unirn(dummy) -1.d0
        randz  = 2.d0 * unirn(dummy) -1.d0
*    Momenta of fragments
        ptmfrg = dsqrt( (etran(1) + frgrm(1))**2 - frgrm(1)**2 )
        pxmfrg(1) = ptmfrg * randx / dsqrt(randx**2 + randy**2
     & + randz**2 )
        pymfrg(1) = ptmfrg * randy / dsqrt(randx**2 + randy**2
     & + randz**2 )
        pzmfrg(1) = ptmfrg * randz / dsqrt(randx**2 + randy**2
     & + randz**2 )
        pxmfrg(2) = - pxmfrg(1)
        pymfrg(2) = - pymfrg(1)
        pzmfrg(2) = - pzmfrg(1)
* shift to laboratory frame   vx, vy, vz are velocity in lab frame
        do 730 icont = 1, multf
         vx = (px * 1.d3 /dsqrt( totmas**2 + (px * 1.d3) **2) +
     & pxmfrg(icont) /dsqrt( frgrm(icont)**2 + pxmfrg(icont)**2))
     & / (1.d0 + px * 1.d3 /dsqrt( totmas**2 + (px * 1.d3) **2) *
     & pxmfrg(icont) /dsqrt( frgrm(icont)**2 + pxmfrg(icont)**2))
         vy = (py * 1.d3 /dsqrt( totmas**2 + (py * 1.d3) **2) +
     & pymfrg(icont) /dsqrt( frgrm(icont)**2 + pymfrg(icont)**2))
     & / (1.d0 + py * 1.d3 /dsqrt( totmas**2 + (py * 1.d3) **2) *
     & pymfrg(icont) /dsqrt( frgrm(icont)**2 + pymfrg(icont)**2))
         vz = (pz * 1.d3 /dsqrt( totmas**2 + (pz * 1.d3) **2) +
     & pzmfrg(icont) /dsqrt( frgrm(icont)**2 + pzmfrg(icont)**2))
     & / (1.d0 + pz * 1.d3 /dsqrt( totmas**2 + (pz * 1.d3) **2) *
     & pzmfrg(icont) /dsqrt( frgrm(icont)**2 + pzmfrg(icont)**2))
*
         pxmfrg(icont) = frgrm(icont) * vx / dsqrt(1.d0 - vx**2
     & - vy**2 - vz**2)
         pymfrg(icont) = frgrm(icont) * vy / dsqrt(1.d0 - vx**2
     & - vy**2 - vz**2)
         pzmfrg(icont) = frgrm(icont) * vz / dsqrt(1.d0 - vx**2
     & - vy**2 - vz**2)
 730    continue

       ELSE
       do 740 icont = 1, multf
*    Sample angular cosines
        randx  = 2.d0 * unirn(dummy) -1.d0
        randy  = 2.d0 * unirn(dummy) -1.d0
        randz  = 2.d0 * unirn(dummy) -1.d0
*    Momenta of fragments
        ptmfrg = dsqrt( (etran(icont) + frgrm(icont))**2
     & - frgrm(icont)**2 )
        pxmfrg(icont) = ptmfrg * randx / dsqrt(randx**2
     & + randy**2 + randz**2 )
        pymfrg(icont) = ptmfrg * randy / dsqrt(randx**2
     & + randy**2 + randz**2 )
        pzmfrg(icont) = ptmfrg * randz / dsqrt(randx**2
     & + randy**2 + randz**2 )
 740   continue
*       Momentum of center-of-mass
        pcmsx = 0.d0
        pcmsy = 0.d0
        pcmsz = 0.d0
       do 750 icont = 1, multf
        pcmsx = pcmsx + pxmfrg(icont)
        pcmsy = pcmsy + pymfrg(icont)
        pcmsz = pcmsz + pzmfrg(icont)
 750   continue

* Correct for CMS translational motion
       do 760 icont = 1, multf
        pxmfrg(icont) = pcmsx * ibarf(icont) / ia0 + pxmfrg(icont)
        pymfrg(icont) = pcmsy * ibarf(icont) / ia0 + pymfrg(icont)
        pzmfrg(icont) = pcmsz * ibarf(icont) / ia0 + pzmfrg(icont)
 760   continue

* Shift to laboratory frame   vx, vy, vz are velocity in lab frame

        do 790 icont = 1, multf
         vx = (px * 1.d3 /dsqrt( totmas**2 + (px * 1.d3) **2) +
     & pxmfrg(icont) /dsqrt( frgrm(icont)**2 + pxmfrg(icont)**2))
     & / (1.d0 + px * 1.d3 /dsqrt( totmas**2 + (px * 1.d3) **2) *
     & pxmfrg(icont) /dsqrt( frgrm(icont)**2 + pxmfrg(icont)**2))
         vy = (py * 1.d3 /dsqrt( totmas**2 + (py * 1.d3) **2) +
     & pymfrg(icont) /dsqrt( frgrm(icont)**2 + pymfrg(icont)**2))
     & / (1.d0 + py * 1.d3 /dsqrt( totmas**2 + (py * 1.d3) **2) *
     & pymfrg(icont) /dsqrt( frgrm(icont)**2 + pymfrg(icont)**2))
         vz = (pz * 1.d3 /dsqrt( totmas**2 + (pz * 1.d3) **2) +
     & pzmfrg(icont) /dsqrt( frgrm(icont)**2 + pzmfrg(icont)**2))
     & / (1.d0 + pz * 1.d3 /dsqrt( totmas**2 + (pz * 1.d3) **2) *
     & pzmfrg(icont) /dsqrt( frgrm(icont)**2 + pzmfrg(icont)**2))
*
         pxmfrg(icont) = frgrm(icont) * vx / dsqrt(1.d0 - vx**2
     & - vy**2 - vz**2)
         pymfrg(icont) = frgrm(icont) * vy / dsqrt(1.d0 - vx**2
     & - vy**2 - vz**2)
         pzmfrg(icont) = frgrm(icont) * vz / dsqrt(1.d0 - vx**2
     & - vy**2 - vz**2)
 790    continue
      ENDIF
*---------------------------------------------------------------

      return
      end subroutine smmexec

************************************************************************
*                                                                      *
         subroutine fbmexec(ia0, iz0, px, py, pz, ex, multf )
*                                                                      *
*        Last Revised:     2014 03 29                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to execute Fermi breakup model                          *
*                                                                      *
*        Variables:                                                    *
*                                                                      *
*        ia0, iz0   : mass and charge number of mother (input)         *
*        ex         : excitation energy of mother (MeV)   (input)      *
*        multf      : multiplicity (output)                            *
*        ibarf(i)   : barion number of i-th fragment (i = 1- multf)    *
*        ichaf(i)   : charge of i-th fragment        (i = 1- multf)    *
*        exmfrg(i)  : excitation energy of fragment  (i = 1- multf)    *
*        p*mfrg(i)  : momentum of fragment           (i = 1- multf)    *
*                                                                      *
************************************************************************

      implicit doubleprecision(a-h,o-z)

*-----------------------------------------------------------------------

      dimension enecn0(500), epsi(500)

      smmsuc = .false.


*-----------------------------------------------------------------------

      maxpfb = int(10.d0**(dble(ia0)/5.d0))

*  Loop for different partition
      do 100 mp = 1, maxpfb

*---------------------------------------------------------------
      wrec(mp) = 0.d0
      trec(mp) = 0.d0
      mulrec(mp) = 0
      nlcrec(mp) = 0

      If(mp .eq. 1) then
       ibar(1) = ia0
       icha(1) = iz0
       mult    = 1
       goto 110
      ENDIF

*  Partition of nucleons (i.e. determine mass of fragments)

      jcont = 1
      ipart = 0

* Sample a large fragment once per 5 trials
       IF( mod(mp, 5) .eq. 0) then
        ibar(1)  = idnint( dble(ia0) * (0.3d0 + 0.5d0 * unirn(dummy) ) )
        ipart    = ibar(1)
        jcont    = 2
       ENDIF

* Expected multiplicity
       mulmea = idint( 5.d0 * unirn(dummy) +1.5d0)

* Nucleon sampling routine #1
       IF( mod(mp, 2) .eq. 0 ) then
        do 200 icont = 1 + ipart, ia0 -1
         IF( unirn(dummy) .lt. (dble(mulmea) -1.d0)/ dble(ia0) ) then
          ibar(jcont) = icont - ipart
          ipart = icont
          jcont = jcont + 1
         ENDIF
 200    continue
       mult   = jcont
       ibar(mult)= ia0 - ipart

* Nucleon sampling routine #2
       ELSE
        mult = mulmea
        do 205 icont = jcont, mult - 1
         ibar(icont) = max(1,int(dint(dble(ia0 - ipart)* unirn(dummy))))
         ipart = ipart + ibar(icont)
 205    continue
       ibar(mult) = ia0 - ipart
       ENDIF

       IF(jcont .eq. 1 .or. ibar(mult) .le. 0)  goto 99
*---------------------------------------------------------------

*---------------------------------------------------------------
* Partition of protons. Giving protons to partitions at random
        ic0    = iz0
        ib0    = ia0
       do 210 jcont = 1, mult-1
        kcont = 0
        do 220 icont = 1, ibar(jcont)
         IF( dble(ic0 - kcont) / dble(max0(1,ib0-icont))
     &                                .gt. unirn(dummy) ) then
          kcont = kcont + 1
         ENDIF
 220    continue
         icha(jcont) = kcont
         ib0 = ib0 - ibar(jcont)
         ic0 = ic0 - icha(jcont)
 210   continue
       icha(mult)= ic0
*---------------------------------------------------------------

*---------------------------------------------------------------
* Destroy neutronium
       do 230 jcont = 1, mult
        IF(icha(jcont) .eq. 0 .and. ibar(jcont) .gt. 1) then
         do 240 icont = mult + 1, mult + ibar(jcont) - 1
          icha(icont) = 0
          ibar(icont) = 1
 240     continue
         mult        = mult + ibar(jcont) - 1
         icha(jcont) = 0
         ibar(jcont) = 1
        ENDIF
 230   continue

* Destroy proton matter
       do 250 jcont = 1, mult
        IF(icha(jcont) .eq. ibar(jcont) .and. ibar(jcont) .gt. 1) then
         do 260 icont = mult + 1, mult + ibar(jcont) - 1
          icha(icont) = 1
          ibar(icont) = 1
 260     continue
         mult        = mult + ibar(jcont) - 1
         icha(jcont) = 1
         ibar(jcont) = 1
        ENDIF
 250   continue
*---------------------------------------------------------------

*---------------------------------------------------------------
* Label #110 is for "no partitioning" check

* chekck if nuclei are between proton/neutron drip lines
 110   do 300 jcont = 1, mult
        IF(ibar(jcont) .gt. maxa(icha(jcont)) .or.
     &     ibar(jcont) .lt. mina(icha(jcont))) then
         goto 99
        ENDIF
 300    continue
        IF(icha(mult) .lt. 0)  goto 99
*---------------------------------------------------------------

*---------------------------------------------------------------
*------------- isotopic multiplicity----------------------------
*  initialize dimension
          nchf  = 0
          nchfd = 0

*  count fragments
        do 310 jcont = 1, mult
        nchf(ibar(jcont),icha(jcont))
     &          = nchf(ibar(jcont),icha(jcont)) + 1

 310    continue

*---------------------------------------------------------------
*------------- Partition duplication check----------------------------

        do icont = 1, mp - 1
         nchfd = nchf
         if(mult .ne. mulrec(icont)) cycle
         do jcont = 1, mulrec(icont)
          nchfd(ibarec(jcont,icont),ichrec(jcont,icont))
     &              = nchfd(ibarec(jcont,icont),ichrec(jcont,icont)) - 1
          if(nchfd(ibarec(jcont,icont),ichrec(jcont,icont)) .lt. 0) exit
         enddo
         if(jcont .eq. mulrec(icont) + 1) goto 99 ! This partition is duplicated
        enddo

*---------------------------------------------------------------
*--------Kinetic energy determination routine-------------------

* Kinetic energy of compound nucleus
      enecn = ex - bindeg( iz0, ia0-iz0) * dble(ia0)
     &  - 6.d-1 * dble(iz0 **2) * ec / (r0f * dble(ia0) ** 0.33333d0
     &   * 2.d0** 0.33333d0 ) / (4.d0 * 3.141d0 * vp ) * 1.d-6

*---------------------------------------------------------------
* Sum-up internal energy of fragmented system

      do jcont = 1, mult
       IF(ibar(jcont) .eq. 1 ) then
        enecn0(jcont) = 0.0d0
       ELSEIF (ibar(jcont) .le. 4) then
        enecn0(jcont) = -bindeg(icha(jcont),ibar(jcont)-icha(jcont))
     &  * dble(ibar(jcont))
       ELSE
        enecn0(jcont) = - bindeg(icha(jcont),ibar(jcont)-icha(jcont))
     &  * dble(ibar(jcont)) - 6.d-1 * dble(icha(jcont)) **2 * ec / (r0f
     &  * dble(ibar(jcont)) ** 0.33333d0 * 2.d0 ** 0.33333d0)
     &  / (4.d0 * 3.141d0 * vp ) * 1.d-6

       ENDIF
      enddo

      ekinfbm(mp) = enecn - sum(enecn0(1:mult))

      if(ekinfbm(mp) .lt. 0.d0) goto 99

*---------------------------------------------------------------


*---------------------------------------------------------------
*--------remember this partition--------------------------

       smmsuc   = .true.

       IF(mult .gt. 1) then
        w = 1.d0
        do jcont = 1, mult
        w = w * (4.d0 * 3.14d0 / 3.d0 * (r0f * 1.d15) **3 * dble(
     &   ibar(jcont) * 2)) / ( (2.d0 * 3.14d0) ** 1.5d0 * hbar **3 )
     &   ** (mult - 1) * Sqrt(ekinfbm(mp)) ** (3 * mult - 5)
     &   / gammaf(1.5d0 * dble(mult - 1) )
     &   / fact( nchf(ibar(jcont),icha(jcont)))
     &   ** (1.d0/ dble( nchf(ibar(jcont),icha(jcont)) ) )

         w = w * (2.d0 * spin(19,icha(jcont)*1000000+ibar(jcont))
     &   + 1.d0) * ibar(jcont) **1.5d0
        enddo
       ELSE
        w = 1.d-50  ! FBM do not allow mult = 1. mult can be 1 only when excitation energy was insufficient
       ENDIF

       do 450 icont=1, mult
        ibarec(icont, mp) = ibar(icont)
        ichrec(icont, mp) = icha(icont)

 450   continue

       mulrec(mp) = mult
       wrec  (mp) = w

*---------------------------------------------------------------

 99   do 600 icont = 1, mult
* Initialize barion, charge assignment not to affect next partition
       ibar = 0
       icha = 0
 600  continue


*---------------------------------------------------------------
 100  continue

      IF( .not. smmsuc) then
          multf = 1

         allocate(ibarf(multf), ichaf(multf), exmfrg(multf),
     &    frgrm(multf), pxmfrg(multf), pymfrg(multf), pzmfrg(multf) )

       ibarf(1) = ia0
       ichaf(1) = iz0
       goto 685
      ENDIF

*----------------------------------------------------------------
*--Determine partition-------------------------------------------

       ransel = unirn(dummy)

       wsum = sum(wrec(1:maxpfb))

       wrecpt = 0.d0
       do icont = 1, maxpfb
         wrecpt = wrec(icont) / wsum + wrecpt
         IF(wrecpt .ge. ransel) then
           ifinal =  icont
         exit
         ENDIF
       end do

*    Score selected partition
       multf = mulrec(ifinal)
       tfin  = 0.666666d0 * ekinfbm(ifinal) / dble(multf)

      allocate( ibarf(multf), ichaf(multf), exmfrg(multf),
     &  frgrm(multf), pxmfrg(multf), pymfrg(multf), pzmfrg(multf) )

       do 680 icont = 1, mulrec(ifinal)
        ibarf(icont) = ibarec(icont, ifinal)
        ichaf(icont) = ichrec(icont, ifinal)
 680   continue


*----------------------------------------------------------------

* ----SMM finalization. Determine energy, momentum, rest mass...-


 685   totmas = 0.d0
       do 690 icont = 1, multf
*  Rest mass
        frgrm(icont) = rnmas * (ibarf(icont) - ichaf(icont))
     &              + rpmas * ichaf(icont)
     &              - bindeg(ichaf(icont),ibarf(icont)-ichaf(icont))
     &              * ibarf(icont)
        totmas       = totmas + frgrm(icont)
 690   continue
*
       do 700 icont = 1, multf
*  Excitation energy
       IF(multf .eq. 1) then
        exmfrg(icont) = ex
       ELSEIF(ibarf(icont) .le. 3) then
        exmfrg(icont) = 0.d0
       ELSE
        exmfrg(icont) = 0.d0  ! problem in FBM. FBM assume 0 but non-zero in reality
       ENDIF

 700   continue

*  Translational motion kinetic energy in CM frame
       IF(multf .eq. 1) then
        etran (1)     = 0.d0
*  Analytical solution of relativistic kinematics for multiplicity = 2
       ELSEIF(multf .eq. 2) then
        etran (1)     = 1.5d0 * tfin * (frgrm(2) + 7.5d-1 * tfin)
     & / ( frgrm(1) + frgrm(2) + 1.5d0 * tfin )
        etran (2)     = 1.5d0 * tfin * (frgrm(1) + 7.5d-1 * tfin)
     & / ( frgrm(1) + frgrm(2) + 1.5d0 * tfin )

       ELSE
* Approximative solution for multiplicity > 2
       do 705 icont = 1, multf
        etran (icont) = - 1.5d0 * tfin * dlog(1.d0 - unirn(dummy))
 705   continue
       ENDIF


*---------------------------------------------------------------
*-----------   Correction for energy conservation --------------
       IF(multf .gt. 1) then
         etot = 0
        do 710 icont = 1, multf
         etot = etot + etran (icont)
 710    continue
*   Correction factor to conserve energy
         etocof = 1.5d0 * dble(multf) * tfin / etot
        do 720 icont = 1, multf
         etran(icont) = etran(icont) * etocof
 720    continue
       ENDIF
*---------------------------------------------------------------


*---------------------------------------------------------------
* ----------   Correction for momentum conservation ------------
       IF(multf .eq. 1) then
        pxmfrg(1) = px * 1.d3
        pymfrg(1) = py * 1.d3
        pzmfrg(1) = pz * 1.d3
       ELSEIF(multf .eq. 2) then
*    Sample angular cosines
        randx  = 2.d0 * unirn(dummy) -1.d0
        randy  = 2.d0 * unirn(dummy) -1.d0
        randz  = 2.d0 * unirn(dummy) -1.d0
*    Momenta of fragments
        ptmfrg = dsqrt( (etran(1) + frgrm(1))**2 - frgrm(1)**2 )
        pxmfrg(1) = ptmfrg * randx / dsqrt(randx**2 + randy**2
     & + randz**2 )
        pymfrg(1) = ptmfrg * randy / dsqrt(randx**2 + randy**2
     & + randz**2 )
        pzmfrg(1) = ptmfrg * randz / dsqrt(randx**2 + randy**2
     & + randz**2 )
        pxmfrg(2) = - pxmfrg(1)
        pymfrg(2) = - pymfrg(1)
        pzmfrg(2) = - pzmfrg(1)
* shift to laboratory frame   vx, vy, vz are velocity in lab frame
        do 730 icont = 1, multf
         vx = (px * 1.d3 /dsqrt( totmas**2 + (px * 1.d3) **2) +
     & pxmfrg(icont) /dsqrt( frgrm(icont)**2 + pxmfrg(icont)**2))
     & / (1.d0 + px * 1.d3 /dsqrt( totmas**2 + (px * 1.d3) **2) *
     & pxmfrg(icont) /dsqrt( frgrm(icont)**2 + pxmfrg(icont)**2))
         vy = (py * 1.d3 /dsqrt( totmas**2 + (py * 1.d3) **2) +
     & pymfrg(icont) /dsqrt( frgrm(icont)**2 + pymfrg(icont)**2))
     & / (1.d0 + py * 1.d3 /dsqrt( totmas**2 + (py * 1.d3) **2) *
     & pymfrg(icont) /dsqrt( frgrm(icont)**2 + pymfrg(icont)**2))
         vz = (pz * 1.d3 /dsqrt( totmas**2 + (pz * 1.d3) **2) +
     & pzmfrg(icont) /dsqrt( frgrm(icont)**2 + pzmfrg(icont)**2))
     & / (1.d0 + pz * 1.d3 /dsqrt( totmas**2 + (pz * 1.d3) **2) *
     & pzmfrg(icont) /dsqrt( frgrm(icont)**2 + pzmfrg(icont)**2))
*
         pxmfrg(icont) = frgrm(icont) * vx / dsqrt(1.d0 - vx**2
     & - vy**2 - vz**2)
         pymfrg(icont) = frgrm(icont) * vy / dsqrt(1.d0 - vx**2
     & - vy**2 - vz**2)
         pzmfrg(icont) = frgrm(icont) * vz / dsqrt(1.d0 - vx**2
     & - vy**2 - vz**2)
 730    continue

       ELSE
       do 740 icont = 1, multf
*    Sample angular cosines
        randx  = 2.d0 * unirn(dummy) -1.d0
        randy  = 2.d0 * unirn(dummy) -1.d0
        randz  = 2.d0 * unirn(dummy) -1.d0
*    Momenta of fragments
        ptmfrg = dsqrt( (etran(icont) + frgrm(icont))**2
     & - frgrm(icont)**2 )
        pxmfrg(icont) = ptmfrg * randx / dsqrt(randx**2
     & + randy**2 + randz**2 )
        pymfrg(icont) = ptmfrg * randy / dsqrt(randx**2
     & + randy**2 + randz**2 )
        pzmfrg(icont) = ptmfrg * randz / dsqrt(randx**2
     & + randy**2 + randz**2 )
 740   continue
*       Momentum of center-of-mass
        pcmsx = 0.d0
        pcmsy = 0.d0
        pcmsz = 0.d0
       do 750 icont = 1, multf
        pcmsx = pcmsx + pxmfrg(icont)
        pcmsy = pcmsy + pymfrg(icont)
        pcmsz = pcmsz + pzmfrg(icont)
 750   continue

* Correct for CMS translational motion
       do 760 icont = 1, multf
        pxmfrg(icont) = pcmsx * ibarf(icont) / ia0 + pxmfrg(icont)
        pymfrg(icont) = pcmsy * ibarf(icont) / ia0 + pymfrg(icont)
        pzmfrg(icont) = pcmsz * ibarf(icont) / ia0 + pzmfrg(icont)
 760   continue

* Shift to laboratory frame   vx, vy, vz are velocity in lab frame

        do 790 icont = 1, multf
         vx = (px * 1.d3 /dsqrt( totmas**2 + (px * 1.d3) **2) +
     & pxmfrg(icont) /dsqrt( frgrm(icont)**2 + pxmfrg(icont)**2))
     & / (1.d0 + px * 1.d3 /dsqrt( totmas**2 + (px * 1.d3) **2) *
     & pxmfrg(icont) /dsqrt( frgrm(icont)**2 + pxmfrg(icont)**2))
         vy = (py * 1.d3 /dsqrt( totmas**2 + (py * 1.d3) **2) +
     & pymfrg(icont) /dsqrt( frgrm(icont)**2 + pymfrg(icont)**2))
     & / (1.d0 + py * 1.d3 /dsqrt( totmas**2 + (py * 1.d3) **2) *
     & pymfrg(icont) /dsqrt( frgrm(icont)**2 + pymfrg(icont)**2))
         vz = (pz * 1.d3 /dsqrt( totmas**2 + (pz * 1.d3) **2) +
     & pzmfrg(icont) /dsqrt( frgrm(icont)**2 + pzmfrg(icont)**2))
     & / (1.d0 + pz * 1.d3 /dsqrt( totmas**2 + (pz * 1.d3) **2) *
     & pzmfrg(icont) /dsqrt( frgrm(icont)**2 + pzmfrg(icont)**2))
*
         pxmfrg(icont) = frgrm(icont) * vx / dsqrt(1.d0 - vx**2
     & - vy**2 - vz**2)
         pymfrg(icont) = frgrm(icont) * vy / dsqrt(1.d0 - vx**2
     & - vy**2 - vz**2)
         pzmfrg(icont) = frgrm(icont) * vz / dsqrt(1.d0 - vx**2
     & - vy**2 - vz**2)
 790    continue
      ENDIF
*---------------------------------------------------------------

      return
      end subroutine fbmexec

***********************************************************************
*                                                                      *
      function eps(iae)
*    Level density parameter of nucleus with mass "iae"

      implicit doubleprecision(a-h,o-z)


*     s is 0.5 when iae .gt. 4, -0.5 when iae .le. 4
      if(iae .eq. 1) then
       eps = 16.d0
      return
      endif

      s = sign(0.5d0, dble(iae - 5))
      e = (0.5d0 - s)
     &  + (0.5d0 + s) * (1.d0 + 3.d0 / (dble(iae) - 1.d0) )
      eps = e * 16.d0
      return
      end function eps
************************************************************************
*                                                                      *


************************************************************************
*                                                                      *
      function vf(iav,m)
*    volume of expanded nucleus iav: nuclear mass, m: multiplicity

      implicit doubleprecision(a-h,o-z)

      vf = 4.d0/3.d0 * dble(iav) * 3.141d0 * r0 **3.d0 * xi(m,iav)

      return
      end function vf
************************************************************************
*                                                                      *

************************************************************************
*                                                                      *
      function xi(m, ia)
*    nuclear volume expansion factor, m: multiplicity, ia: mother nuclei mass

      implicit doubleprecision(a-h,o-z)

      xi = (1.d0 + 1.4d-15 * ( dble(m) **0.33333d0 -1.d0) /
     & (r0 * dble(ia) ** 0.3333333d0) ) **3.d0 - 1.d0

      return
      end function xi
************************************************************************
*                                                                      *

************************************************************************
*                                                                      *
      function tlam(t)
*    de Broglie wave length (fm) of nucleons at temperature t (MeV)

      implicit doubleprecision(a-h,o-z)

      tlam = dsqrt(6.626068d-34 **2 /(2.d0 * 3.141d0) /
     & (1.67262158d-27 * t * 1.0d6 * ec))

      return
      end function tlam
************************************************************************
*                                                                      *

************************************************************************
*                                                                      *
      function sint(t,ia,iz,ia0,iz0)
*    Partial entropy of fragments

*=====================================================================
* <variables>
*    ia   :   mass of fragment                      (IN)
*    iz   :   charge  of fragment                   (IN)
*     t   :   temperature of compound nucleus       (IN)
*   ia0   :   mass  of original nucleus             (IN)
*   iz0   :   charge of original nucleus            (IN)
*   sint  :   Entropy           {MeV]              (OUT)
*/////////////////////////////////////////////////////////////////////

      implicit doubleprecision(a-h,o-z)

      IF(ia .le. 3) then
       sint = 0.d0
      ELSEIF(ia .eq. 4 .and. iz .eq. 2) then
       sint = 2.d0 * t / eps(ia) * dble(ia)
      ELSE
       sint = 2.d0 * t / eps(ia) * dble(ia)
     &  + 18.d0 * 5.d0 * t * tc**2
     &  * ((tc**2 - t**2)/(tc**2 + t**2)) ** 0.25d0
     &  / (t**2 + tc**2)**2 * dble(ia) ** 0.6666666d0
      ENDIF

      return
      end function sint
************************************************************************
*                                                                      *

************************************************************************
*                                                                      *
      function fact(num)
*     factorial funtion

      implicit doubleprecision(a-h,o-z)

      fact = 1.d0
      do 10 icont = 0, num
       fact = fact * dble( max0(1, icont ) )
 10   continue

      return
      end function fact
************************************************************************
*                                                                      *

************************************************************************
*                                                                      *
      function gammaf(arg)
*     gamma funtion

      implicit doubleprecision(a-h,o-z)

      gammaf = sqrt( 3.141592653589793d0 )
      do icont = 1, idnint(arg) - 2 ! round up
       gammaf = gammaf * (arg - dble(icont) )
      enddo

      return
      end function gammaf
************************************************************************
*                                                                      *

      end module smmmod

