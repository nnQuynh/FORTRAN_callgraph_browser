************************************************************************
*                                                                      *
        subroutine aamain(itz,ita)
*                                                                      *
*       control routine of AAMA program                                *
*       original file was made by V.R.Akylas and P.Vogel               *
*       ref.) Muonic atom cascade program,                             *
*             Computer Physics Communications, vol.15, p.291, 1978     *
*                                                                      *
*       modified by S.Abe on 2015/03/05                                *
*       last modified by S.Abe on 2018/10/30                           *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      common /paran/ icfn(100), ilfn(100), chfn(100)
      character chfn*200
      integer icfn,ilfn

      logical   exex

*-----------------------------------------------------------------------

      common /muint/ imuint, imubrm, imuppd, imucap

      common /aama01/ energ, econs, econst, d2p1s, d2p1sm
!$OMP THREADPRIVATE(/aama01/)
      common /aama02/ zsa(3), be(3), bem(3)
!$OMP THREADPRIVATE(/aama02/)
      common /aama03/ k0, k1, k2, k3
!$OMP THREADPRIVATE(/aama03/)
      common /aama04/ nn0(3), nn1(7), nn2(7), nn3(7)
!$OMP THREADPRIVATE(/aama04/)
      common /aama05/ ip1(7), ip2(7), ip3(7), iq1(7), iq2(7), iq3(7)
!$OMP THREADPRIVATE(/aama05/)
      common /aama06/ z, zsk, zsl, zsm, zskz, zslz, zsmz
!$OMP THREADPRIVATE(/aama06/)
      common /aama11/ rr(18), rau, rad, ra(4), rd(4), rsa(4)
!$OMP THREADPRIVATE(/aama11/)
      common /aama12/ pop(6)
!$OMP THREADPRIVATE(/aama12/)
      common /aama13/ e(1000), ai(1000), energya(20,40), m, ia(1000)
!$OMP THREADPRIVATE(/aama13/)
      common /aama14/ dza, dza2, dredm
!$OMP THREADPRIVATE(/aama14/)
      common /aama15/ pl(20), pln(210)
!$OMP THREADPRIVATE(/aama15/)
      common /aama16/ a, cfm
!$OMP THREADPRIVATE(/aama16/)
      common /aama17/ yc(4)
!$OMP THREADPRIVATE(/aama17/)

      common /aama99/ popul(20),ptran(210,210),cctr(210,210),
     &                esp0,energy0(20,40)
!$OMP THREADPRIVATE(/aama99/)
      common /aamai99/ nmax0
!$OMP THREADPRIVATE(/aamai99/)

      common /aamac01/ iread,iw,ipunch,iprint
!$OMP THREADPRIVATE(/aamac01/)
      common /aamac02/ ffact(60),ffactd
!$OMP THREADPRIVATE(/aamac02/)
      common /aamac03/ pi,coeff,picoef,amassm,amasse,alfa,ne(3),jk(3)
!$OMP THREADPRIVATE(/aamac03/)
      common /aamac04/ amassa,amassn
!$OMP THREADPRIVATE(/aamac04/)
      common /aamac05/ m1(7),m2(7),m3(7),idb
!$OMP THREADPRIVATE(/aamac05/)
      common /aamac06/ zmk,zml,zmm,zmkm,zmlm,zmmm,cl1,cl2,ivers
!$OMP THREADPRIVATE(/aamac06/)
      common /aamac07/ jtm(6),jtd(6),jtq(6),jto(6)
!$OMP THREADPRIVATE(/aamac07/)
      common /aamac08/ coemon(30),expmon(30),ifm(6),jm(10)
!$OMP THREADPRIVATE(/aamac08/)
      common /aamac09/ coedip(42),expdip(42),coed(4),coedp(9),
     &                 ifd(9),jd(14)
!$OMP THREADPRIVATE(/aamac09/)
      common /aamac10/ coequa(45),expqua(45),coeq(11),ifq(10),jq(15)
!$OMP THREADPRIVATE(/aamac10/)
      common /aamac11/ coeoct(45),expoct(45),coeo(11),ifo(10),jo(15)
!$OMP THREADPRIVATE(/aamac11/)
      common /aamac12/ hbar,widthk,npol(20),ipol,ide,ip8,ipc(3),nmax
!$OMP THREADPRIVATE(/aamac12/)
      common /aamac13/ ehigh,elow,climit,eres,esp,espm,cd(5),ea,eb,
     &                 icc,mpu,icpu(200),ipn
!$OMP THREADPRIVATE(/aamac13/)
      common /aamac14/ tfm,step,rmatch,alexp,nopt,idir,iyc
!$OMP THREADPRIVATE(/aamac14/)
      common /aamac15/ ij(4),yj(4),jj1(4)
!$OMP THREADPRIVATE(/aamac15/)
      common /aamac16/ ic,ll(20)
!$OMP THREADPRIVATE(/aamac16/)
      common /aamac17/ elecbe(104,3)
      common /aamadbg/ idebug

*-----------------------------------------------------------------------
c  ***note -- non-standard way of overprinting -- change as needed
c  ***standard fortran requires data for i9 /1h ,1h+,1h+,1h+/

      dimension i9(4)
      character i9*1
      data i9 / '+', '+', '+', ' '/

      dimension m9(4,7), k9(4,4)
      character m9*2, l9*2, k9*2
      data l9 / '--'/
      data k9 / 'rd', '1s', '2t', '3t', 'rd',
     &          '1s', '2s', '3s', 'rd', '**',
     &          '2p', '3p', 'rd', '**', '**',
     &          '3d'/

      dimension j4(9,5,4), jx(4)
      character j4*1, jx*1, jb*1
      data jx / 'x', 'i', 'h', 'o'/
      data jb / ' '/

      dimension j9(9,10), id9(5), ixx(3,20), plx(20)
      data j9 / 124, 254, 455, 387, 387, 387, 455, 254, 124, 024,
     &          056, 120, 024, 024, 024, 024, 126, 126, 124, 254,
     &          455, 391, 030, 056, 112, 255, 511, 511, 510, 028,
     &          112, 056, 014, 391, 254, 124, 014, 030, 054, 102,
     &          255 ,511, 006, 015, 015, 511, 511, 384, 508, 254,
     &          007, 391, 254, 124, 124, 254, 385 ,508, 510, 387,
     &          455, 254, 124, 511, 511, 007, 014, 028, 056, 112,
     &          224, 448, 124 ,254, 387, 455, 254, 455, 387, 254,
     &          124, 124, 254, 455, 387, 255, 127, 259, 254 ,124/

*-----------------------------------------------------------------------

      character chln*200

      real*8,allocatable,save :: 
     &           d2p1s_udf(:),esp_udf(:),
     &           yc_udf(:,:),pln_udf(:,:),pl_udf(:,:),
     &           cfm_udf(:),tfm_udf(:),
     &           zsk_udf(:),zsl_udf(:),zsm_udf(:),
     &           be_udf(:,:),pop_udf(:,:),alexp_udf(:),
     &           cl1_udf(:),cl2_udf(:),ehigh_udf(:),
     &           elow_udf(:),climit_udf(:),eres_udf(:),
     &           amassa_udf(:),zmk_udf(:),zml_udf(:),
     &           zmm_udf(:),zmkm_udf(:),zmlm_udf(:),
     &           zmmm_udf(:),yj_udf(:,:)

      integer,allocatable,save ::
     &           k0_udf(:),k1_udf(:),
     &           k2_udf(:),k3_udf(:),
     &           nn0_udf(:,:),nn1_udf(:,:),
     &           nn2_udf(:,:),nn3_udf(:,:),
     &           m1_udf(:,:),m2_udf(:,:), m3_udf(:,:),
     &           ip1_udf(:,:),ip2_udf(:,:),ip3_udf(:,:),
     &           iq1_udf(:,:),iq2_udf(:,:),iq3_udf(:,:),
     &           jtm_udf(:,:),jtd_udf(:,:),jtq_udf(:,:),
     &           jto_udf(:,:),widthk_udf(:),ipc_udf(:,:),
     &           nopt_udf(:),nmax_udf(:),ip8_udf(:),
     &           npol_udf(:,:),ipol_udf(:),
     &           jm_udf(:,:),jd_udf(:,:),jq_udf(:,:),
     &           jo_udf(:,:),iyc_udf(:),
     &           ifm_udf(:,:),ifd_udf(:,:),ifq_udf(:,:),
     &           ifo_udf(:,:),ij_udf(:,:),jj1_udf(:,:),
     &           iz_udf(:),ia_udf(:)

      integer,save:: ireadflg, nelmnt
      data ireadflg /0/
      data nelmnt /0/

*-----------------------------------------------------------------------
*     open file(28) (D=aama.dat) and allocate udf arrays
*-----------------------------------------------------------------------

!$OMP CRITICAL (aamadat_read)
      if( ireadflg.eq. 0 ) then

       ireadflg = 1
       nelmnt = 0

*-----------------------------------------------------------------------

       if( imucap .eq. 2 ) then

        inquire( file = chfn(28)(1:ilfn(28)), exist = exex )
        if( exex .eqv. .false. ) then
         write(*,'(/'' Error: aama parameter file does not exist.''/
     &              '' file name = '',200a1)')
     &    ( chfn(28)(i:i), i=1,ilfn(28))
         stop
        endif

        iotaama = 999
        open(iotaama, file=chfn(28)(1:ilfn(28)), form = 'formatted')

 1000   read(iotaama,'(a200)', end=1998) chln
        if( chln(1:3) .eq. 'tot' ) then
         read(chln(11:200),'(i190)') nelmnt
         goto 1999
        else
         goto 1000
        endif

 1998   continue
        write(*,'('' Warning: parameter "tot" is not found in'',
     &            '' aama parameter file.''/
     &            '' default parameter set is applied.'')')

 1999   continue

       endif

*-----------------------------------------------------------------------

       allocate( d2p1s_udf(0:nelmnt),esp_udf(0:nelmnt),
     &           k0_udf(0:nelmnt),k1_udf(0:nelmnt),
     &           k2_udf(0:nelmnt),k3_udf(0:nelmnt),
     &           nn0_udf(0:nelmnt,3),nn1_udf(0:nelmnt,7),
     &           nn2_udf(0:nelmnt,7),nn3_udf(0:nelmnt,7),
     &           m1_udf(0:nelmnt,7),m2_udf(0:nelmnt,7),
     &           m3_udf(0:nelmnt,7),
     &           ip1_udf(0:nelmnt,7),ip2_udf(0:nelmnt,7),
     &           ip3_udf(0:nelmnt,7),
     &           iq1_udf(0:nelmnt,7),iq2_udf(0:nelmnt,7),
     &           iq3_udf(0:nelmnt,7),
     &           yc_udf(0:nelmnt,4),
     &           pln_udf(0:nelmnt,210),pl_udf(0:nelmnt,20),
     &           cfm_udf(0:nelmnt),tfm_udf(0:nelmnt),
     &           zsk_udf(0:nelmnt),zsl_udf(0:nelmnt),zsm_udf(0:nelmnt),
     &           be_udf(0:nelmnt,3),pop_udf(0:nelmnt,6),
     &           jtm_udf(0:nelmnt,6),jtd_udf(0:nelmnt,6),
     &           jtq_udf(0:nelmnt,6),jto_udf(0:nelmnt,6),
     &           widthk_udf(0:nelmnt),ipc_udf(0:nelmnt,3),
     &           nopt_udf(0:nelmnt),nmax_udf(0:nelmnt),
     &           alexp_udf(0:nelmnt),
     &           cl1_udf(0:nelmnt),cl2_udf(0:nelmnt),
     &           ip8_udf(0:nelmnt),
     &           npol_udf(0:nelmnt,20),ipol_udf(0:nelmnt),
     &           ehigh_udf(0:nelmnt),elow_udf(0:nelmnt),
     &           climit_udf(0:nelmnt),eres_udf(0:nelmnt),
     &           jm_udf(0:nelmnt,10),jd_udf(0:nelmnt,14),
     &           jq_udf(0:nelmnt,15),jo_udf(0:nelmnt,15),
     &           iyc_udf(0:nelmnt),amassa_udf(0:nelmnt),
     &           ifm_udf(0:nelmnt,6),ifd_udf(0:nelmnt,9),
     &           ifq_udf(0:nelmnt,10),ifo_udf(0:nelmnt,10),
     &           zmk_udf(0:nelmnt),zml_udf(0:nelmnt),
     &           zmm_udf(0:nelmnt),zmkm_udf(0:nelmnt),
     &           zmlm_udf(0:nelmnt),zmmm_udf(0:nelmnt),
     &           ij_udf(0:nelmnt,4),yj_udf(0:nelmnt,4),
     &           jj1_udf(0:nelmnt,4),
     &           iz_udf(0:nelmnt),ia_udf(0:nelmnt) )

*-----------------------------------------------------------------------
*     set default value for each parameter
*-----------------------------------------------------------------------

       z = 0.d0
       a = 0.d0

       d2p1s = 0.d0

       k0 = 3
       k1 = 4
       k2 = 4
       k3 = 4

       nn0(1) = 1
       nn0(2) = 2
       nn0(3) = 3

       nn1(1) = 0
       nn1(2) = 1
       nn1(3) = 2
       nn1(4) = 3
       nn1(5) = 3
       nn1(6) = 3
       nn1(7) = 3

       ip1(1) = 0
       do i = 2, 7
        ip1(i) = 1
       enddo

       do i = 1, 7
        iq1(i) = 0
       enddo

       do i = 1, 7
        nn2(i) = nn1(i)
        nn3(i) = nn1(i)
        ip2(i) = ip1(i)
        ip3(i) = ip1(i)
        iq2(i) = iq1(i)
        iq3(i) = iq1(i)
       enddo

       do i = 1, 4
        yc(i) = 1.d0
       enddo

       do i = 1, 210
        pln(i) = 0.d0
       enddo

       do i = 1, 20
        pl(i) = 0.d0
       enddo

       cfm = 0.d0

       zsk = -1.d0
       zsl = -1.d0
       zsm = -1.d0

       do i = 1, 3
        be(i) = -1.d0
       enddo

       do i = 1, 6
        pop(i) = -1.d0
       enddo

*-----------------------------------------------------------------------

       do iel = 0, nelmnt

*-----------------------------------------------------------------------
*     assign input parameters to udf array, iel=0: default parameters
*-----------------------------------------------------------------------

        iz_udf(iel) = idint(z)
        ia_udf(iel) = idint(a)

        d2p1s_udf(iel) = d2p1s
        esp_udf(iel) = esp

        k0_udf(iel) = k0
        k1_udf(iel) = k1
        k2_udf(iel) = k2
        k3_udf(iel) = k3

        do i = 1, 3
         nn0_udf(iel,i) = nn0(i)
        enddo

        do i = 1, 7
         nn1_udf(iel,i) = nn1(i)
         nn2_udf(iel,i) = nn2(i)
         nn3_udf(iel,i) = nn3(i)
         ip1_udf(iel,i) = ip1(i)
         ip2_udf(iel,i) = ip2(i)
         ip3_udf(iel,i) = ip3(i)
         iq1_udf(iel,i) = iq1(i)
         iq2_udf(iel,i) = iq2(i)
         iq3_udf(iel,i) = iq3(i)
        enddo

        do i = 1, 4
         yc_udf(iel,i) = yc(i)
        enddo

        do i = 1, 7
         m1_udf(iel,i) = m1(i)
         m2_udf(iel,i) = m2(i)
         m3_udf(iel,i) = m3(i)
        enddo

        do i = 1, 210
         pln_udf(iel,i) = pln(i)
        enddo

        do i = 1, 20
         pl_udf(iel,i) = pl(i)
        enddo

        cfm_udf(iel) = cfm
        tfm_udf(iel) = tfm

        zsk_udf(iel) = zsk
        zsl_udf(iel) = zsl
        zsm_udf(iel) = zsm

        do i = 1, 3
         be_udf(iel,i) = be(i)
        enddo

        do i = 1, 6
         pop_udf(iel,i) = pop(i)
        enddo

        do i = 1, 6
         jtm_udf(iel,i) = jtm(i)
         jtd_udf(iel,i) = jtd(i)
         jtq_udf(iel,i) = jtq(i)
         jto_udf(iel,i) = jto(i)
        enddo

        widthk_udf(iel) = widthk
        do i = 1, 3
         ipc_udf(iel,i) = ipc(i)
        enddo
        nopt_udf(iel) = nopt
        nmax_udf(iel) = nmax
        alexp_udf(iel) = alexp
        cl1_udf(iel) = cl1
        cl2_udf(iel) = cl2
        ip8_udf(iel) = ip8
        do i = 1, 20
         npol_udf(iel,i) = npol(i)
        enddo
        ipol_udf(iel) = ipol
        ehigh_udf(iel) = ehigh
        elow_udf(iel) = elow
        climit_udf(iel) = climit
        eres_udf(iel) = eres

        do i = 1, 10
         jm_udf(iel,i) = jm(i)
        enddo
        do i = 1, 14
         jd_udf(iel,i) = jd(i)
        enddo
        do i = 1, 15
         jq_udf(iel,i) = jq(i)
        enddo
        do i = 1, 15
         jo_udf(iel,i) = jo(i)
        enddo

        iyc_udf(iel) = iyc
        amassa_udf(iel) = amassa

        do i = 1, 6
         ifm_udf(iel,i) = ifm(i)
        enddo
        do i = 1, 9
         ifd_udf(iel,i) = ifd(i)
        enddo
        do i = 1, 10
         ifq_udf(iel,i) = ifq(i)
        enddo
        do i = 1, 10
         ifo_udf(iel,i) = ifo(i)
        enddo

        zmk_udf(iel) = zmk
        zml_udf(iel) = zml
        zmm_udf(iel) = zmm

        zmkm_udf(iel) = zmkm
        zmlm_udf(iel) = zmlm
        zmmm_udf(iel) = zmmm

        do i = 1, 4
         ij_udf(iel,i) = ij(i)
         yj_udf(iel,i) = yj(i)
         jj1_udf(iel,i) = jj1(i)
        enddo

*-----------------------------------------------------------------------
*     exit when all udf parameter sets are stored
*-----------------------------------------------------------------------

        if( iel .eq. nelmnt ) exit

*-----------------------------------------------------------------------
*     initialize input parameters before loading file(28)
*-----------------------------------------------------------------------

        d2p1s = d2p1s_udf(0)
        esp = esp_udf(0)

        k0 = k0_udf(0)
        k1 = k1_udf(0)
        k2 = k2_udf(0)
        k3 = k3_udf(0)

        do i = 1, 3
         nn0(i) = nn0_udf(0,i)
        enddo

        do i = 1, 7
         nn1(i) = nn1_udf(0,i)
         nn2(i) = nn2_udf(0,i)
         nn3(i) = nn3_udf(0,i)
         ip1(i) = ip1_udf(0,i)
         ip2(i) = ip2_udf(0,i)
         ip3(i) = ip3_udf(0,i)
         iq1(i) = iq1_udf(0,i)
         iq2(i) = iq2_udf(0,i)
         iq3(i) = iq3_udf(0,i)
        enddo

        do i = 1, 4
         yc(i) = yc_udf(0,i)
        enddo

        do i = 1, 7
         m1(i) = m1_udf(0,i)
         m2(i) = m2_udf(0,i)
         m3(i) = m3_udf(0,i)
        enddo

        do i = 1, 210
         pln(i) = pln_udf(0,i)
        enddo

        do i = 1, 20
         pl(i) = pl_udf(0,i)
        enddo

        cfm = cfm_udf(0)
        tfm = tfm_udf(0)

        zsk = zsk_udf(0)
        zsl = zsl_udf(0)
        zsm = zsm_udf(0)

        do i = 1, 3
         be(i) = be_udf(0,i)
        enddo

        do i = 1, 6
         pop(i) = pop_udf(0,i)
        enddo

        do i = 1, 6
         jtm(i) = jtm_udf(0,i)
         jtd(i) = jtd_udf(0,i)
         jtq(i) = jtq_udf(0,i)
         jto(i) = jto_udf(0,i)
        enddo

        widthk = widthk_udf(0)
        do i = 1, 3
         ipc(i) = ipc_udf(0,i)
        enddo
        nopt = nopt_udf(0)
        nmax = nmax_udf(0)
        alexp = alexp_udf(0)
        cl1 = cl1_udf(0)
        cl2 = cl2_udf(0)
        ip8 = ip8_udf(0)
        do i = 1, 20
         npol(i) = npol_udf(0,i)
        enddo
        ipol = ipol_udf(0)

        ehigh = ehigh_udf(0)
        elow = elow_udf(0)
        climit = climit_udf(0)
        eres = eres_udf(0)

        do i = 1, 10
         jm(i) = jm_udf(0,i)
        enddo
        do i = 1, 14
         jd(i) = jd_udf(0,i)
        enddo
        do i = 1, 15
         jq(i) = jq_udf(0,i)
        enddo
        do i = 1, 15
         jo(i) = jo_udf(0,i)
        enddo

        iyc = iyc_udf(0)
        amassa = amassa_udf(0)

        do i = 1, 6
         ifm(i) = ifm_udf(0,i)
        enddo
        do i = 1, 9
         ifd(i) = ifd_udf(0,i)
        enddo
        do i = 1, 10
         ifq(i) = ifq_udf(0,i)
        enddo
        do i = 1, 10
         ifo(i) = ifo_udf(0,i)
        enddo

        zmk = zmk_udf(0)
        zml = zml_udf(0)
        zmm = zmm_udf(0)

        zmkm = zmkm_udf(0)
        zmlm = zmlm_udf(0)
        zmmm = zmmm_udf(0)

        do i = 1, 4
         ij(i) = ij_udf(0,i)
         yj(i) = yj_udf(0,i)
         jj1(i) = jj1_udf(0,i)
        enddo

*-----------------------------------------------------------------------
*     read one parameter set from file(28)
*-----------------------------------------------------------------------

 1100   read(iotaama, '(a200)', end=1199) chln
        if( chln(1:3) .ne. 'iel' ) goto 1100

        call aamaread(iotaama)

       enddo

*-----------------------------------------------------------------------
*     check the number of reading parameter sets
*-----------------------------------------------------------------------

 1199  continue

       if( imucap .eq. 2 ) close(iotaama)
       if( iel .lt. nelmnt ) then
        write(*,'('' Warning: lack of parameter set in'',
     &            '' aama parameter file.''/
     &            '' default parameters fill deficient sets.''/
     &            '' iel, nelmnt = '',2i8)') iel, nelmnt
        nelmnt = iel
       endif

      endif
!$OMP END CRITICAL (aamadat_read)

*-----------------------------------------------------------------------
*     search iel and load parameters from udf array
*-----------------------------------------------------------------------

      if( nelmnt .gt. 0 ) then
       do iel = 1, nelmnt
         if( itz .eq. iz_udf(iel) .and. ita .eq. ia_udf(iel) ) goto 100
       enddo
      endif
      iel = 0
  100 continue

*-----------------------------------------------------------------------

      econs = 5.505355d-3
      d2p1s = d2p1s_udf(iel)
      esp = esp_udf(iel)

      k0 = k0_udf(iel)
      k1 = k1_udf(iel)
      k2 = k2_udf(iel)
      k3 = k3_udf(iel)

      do i = 1, 3
       nn0(i) = nn0_udf(iel,i)
      enddo

      do i = 1, 7
       nn1(i) = nn1_udf(iel,i)
       nn2(i) = nn2_udf(iel,i)
       nn3(i) = nn3_udf(iel,i)
       ip1(i) = ip1_udf(iel,i)
       ip2(i) = ip2_udf(iel,i)
       ip3(i) = ip3_udf(iel,i)
       iq1(i) = iq1_udf(iel,i)
       iq2(i) = iq2_udf(iel,i)
       iq3(i) = iq3_udf(iel,i)
      enddo

      do i = 1, 4
       yc(i) = yc_udf(iel,i)
      enddo

      do i = 1, 7
       m1(i) = m1_udf(iel,i)
       m2(i) = m2_udf(iel,i)
       m3(i) = m3_udf(iel,i)
      enddo

      do i = 1, 210
       pln(i) = pln_udf(iel,i)
      enddo

      do i = 1, 20
       pl(i) = pl_udf(iel,i)
      enddo

      m = 1

      do i = 1, 20
       do j = 1, 40
        energya(i,j) = 0.0d0
       enddo
      enddo

      cfm = cfm_udf(iel)
      tfm = tfm_udf(iel)

*-----------------------------------------------------------------------
*     set parameters from PHITS
*-----------------------------------------------------------------------

      z = dble(itz)
      a = dble(ita)

      zsk = z - 1
      zsl = z - 3.0d0
      zsm = z - 11.0d0

      if( zsk_udf(iel) .ge. 0.d0 ) zsk = zsk_udf(iel)
      if( zsl_udf(iel) .ge. 0.d0 ) zsl = zsl_udf(iel)
      if( zsm_udf(iel) .ge. 0.d0 ) zsm = zsm_udf(iel)

      if( zsl .lt. 0.d0 ) zsl = 0.0d0
      if( zsm .lt. 0.d0 ) zsm = 0.0d0

      do i = 1, 3
       if( be_udf(iel,i) .ge. 0.d0 ) then
        be(i) = be_udf(iel,i)
       else
cABE 2022/11/04
        if( itz .eq. 1 ) then
         be(i) = 0.d0
        else
         be(i) = elecbe(itz-1,i)
        endif
       endif
      enddo

      pop(1) = z / 2.0d0
      pop(2) = (z-2.0d0)  /  2.0d0
      pop(3) = (z-4.0d0)  /  6.0d0
      pop(4) = (z-10.0d0) /  2.0d0
      pop(5) = (z-12.0d0) /  6.0d0
      pop(6) = (z-18.0d0) / 10.0d0

      do i = 1, 6
       if( pop_udf(iel,i) .ge. 0.d0 ) pop(i) = pop_udf(iel,i)
       if( pop(i) .ge. 1.0d0 ) pop(i) = 1.0d0
       if( pop(i) .le. 0.0d0 ) pop(i) = 0.0d0
      enddo

      if( itz .eq. 19 .or. itz .eq. 20 ) then
       pop(6) = 0.0d0
      elseif( itz .eq. 21 ) then
       pop(6) = 0.1d0
      elseif( itz .eq. 22 ) then
       pop(6) = 0.2d0
      elseif( itz .eq. 23 ) then
       pop(6) = 0.3d0
      elseif( itz .eq. 24 .or. itz .eq. 25 ) then
       pop(6) = 0.5d0
      elseif( itz .eq. 26 ) then
       pop(6) = 0.6d0
      elseif( itz .eq. 27 ) then
       pop(6) = 0.7d0
      elseif( itz .eq. 28 ) then
       pop(6) = 0.8d0
      elseif( itz .eq. 29 ) then
       pop(6) = 1.0d0
      endif

      do i = 1, 6
       jtm(i) = jtm_udf(iel,i)
       jtd(i) = jtd_udf(iel,i)
       jtq(i) = jtq_udf(iel,i)
       jto(i) = jto_udf(iel,i)
      enddo

*-----------------------------------------------------------------------
*     set parameters defined in data block
*-----------------------------------------------------------------------

      widthk = widthk_udf(iel)
      do i = 1, 3
       ipc(i) = ipc_udf(iel,i)
      enddo
      nopt = nopt_udf(iel)
      nmax = nmax_udf(iel)
      alexp = alexp_udf(iel)
      cl1 = cl1_udf(iel)
      cl2 = cl2_udf(iel)
      ip8 = ip8_udf(iel)
      do i = 1, 20
       npol(i) = npol_udf(iel,i)
      enddo
      ipol = ipol_udf(iel)

      ehigh = ehigh_udf(iel)
      elow = elow_udf(iel)
      climit = climit_udf(iel)
      eres = eres_udf(iel)

*-----------------------------------------------------------------------

      if(idebug.eq.1)then
20000  if(z.le.0.0d0) write(iw,20100)
       if(z.gt.99.0d0) write(iw,20200)
       if(z.gt.137.04d0) write(iw,20300)
       if(dmod(z+0.001d0,1.0d0).gt.10.002) write(iw,20400)
       if(z.le.0.0d0) stop
20100  format(53h *** error *** atomic number z zero or negative *** e,
     &60hxecution terminated *** ....................................)
20200  format(53h *** warning *** atomic number z too big *** last two,
     &60hdigits will be printed in block letters in the upper left of,
     &8htable **)
20300  format(53h *** warning *** atomic number z too big *** point li,
     &60hke dirac formulae have problems *** no attempt to rectify pr,
     &8hoblem **)
20400  format(53h *** warning *** atomic number z not close to an inte,
     &60hger value *** integer part printed in table,but actual value,
     &8h used **)

       if(z.lt.30.0d0 .and. (z-20.0d0)*0.1d0.gt.pop(6)) write(iw,27500)
       if(z.lt.18.0d0 .and. (z-12.0d0)/6.0d0.gt.pop(5)) write(iw,27500)
       if(z.lt.12.0d0 .and. (z-10.0d0)*0.5d0.gt.pop(4)) write(iw,27500)
       if(z.lt.10.0d0 .and. (z-4.0d0)/6.0d0.gt.pop(3)) write(iw,27500)
       if(z.lt.4.0d0 .and. (z-2.0d0)*0.5d0.gt.pop(2)) write(iw,27500)
       if((z.le.20.0d0 .and. pop(6).ne.0.0d0) .or.
     &    (z.le.12.0d0 .and. pop(5).ne.0.0d0) .or.
     &    (z.le.10.0d0 .and. pop(4).ne.0.0d0) .or.
     &    (z.le.4.0d0 .and. pop(3).ne.0.0d0) .or.
     &    (z.le.2.0d0 .and. pop(2).ne.0.0d0)) write(iw,27500)
27500  format(53h *** error *** z is too low to support specified popu,
     &60hlation of electronic shells *** possible erroneous results *,
     &8h********)

       if( nopt .ne. 0 ) then
        write(6,*)"error: l,pl(l) or cl1,cl2 are necessary for nop =/=0"
        stop
       endif
      endif

*-----------------------------------------------------------------------
*     set unknown input parameters
*-----------------------------------------------------------------------

      do i = 1, 10
       jm(i) = jm_udf(iel,i)
      enddo
      do i = 1, 14
       jd(i) = jd_udf(iel,i)
      enddo
      do i = 1, 15
       jq(i) = jq_udf(iel,i)
      enddo
      do i = 1, 15
       jo(i) = jo_udf(iel,i)
      enddo

      iyc = iyc_udf(iel)
      amassa = amassa_udf(iel)

      do i = 1, 6
       ifm(i) = ifm_udf(iel,i)
      enddo
      do i = 1, 9
       ifd(i) = ifd_udf(iel,i)
      enddo
      do i = 1, 10
       ifq(i) = ifq_udf(iel,i)
      enddo
      do i = 1, 10
       ifo(i) = ifo_udf(iel,i)
      enddo

      zmk = zmk_udf(iel)
      zml = zml_udf(iel)
      zmm = zmm_udf(iel)

      zmkm = zmkm_udf(iel)
      zmlm = zmlm_udf(iel)
      zmmm = zmmm_udf(iel)

      do i = 1, 4
       ij(i) = ij_udf(iel,i)
       yj(i) = yj_udf(iel,i)
       jj1(i) = jj1_udf(iel,i)
      enddo

*-----------------------------------------------------------------------

      call aamaffix

*-----------------------------------------------------------------------
*     replace energya(n,l)
*-----------------------------------------------------------------------

      do n = 1, nmax

       k = 2*n - 1

       do l = 1, k

        l1 = (l-1)/2 + 1
        nj = n
        pt = point(nj,l1)

        if(idebug.eq.1)then
        if ( energya(n,l) .le. 0.0d0 .and. mod(iprint,2) .eq. 0 )
     &   write(iw,1200) n, l, pt
 1200   format(27h no input data for state n=,i2,7h, (lj)=,i2,
     &         29h,  point-like dirac energy = ,f10.6,5h(MeV))
        endif

        if ( energya(n,l) .le. 0.0d0 )
     &   energya(n,l) = pt - energya(n,l)

       enddo

      enddo

*-----------------------------------------------------------------------
*     nmax: starting principal quantum number for the cascade
*-----------------------------------------------------------------------

      if ( nmax .lt. 20 ) then

       mx1 = nmax + 1

       do i = mx1, 20
        pl(i) = 0.0d0
       enddo

      endif

*-----------------------------------------------------------------------
*     modified angular momentum distribution, pl(i)
*     nop = -1: use set parameter, no modification
*         =  0: modified statistical
*         =  2: modified quadratic
*-----------------------------------------------------------------------

      if ( nopt .eq. 0 ) then
       do i = 1, nmax
        pl(i) = dble(2*i-1) * dexp(alexp*dble(i-1))
       enddo
      elseif ( nopt .eq. 2 ) then
       do i = 1, nmax
        pl(i) = 1.0d0 + cl1*dble(i-1) + cl2*dble((i-1)**2)
       enddo
      endif

      ss = 0.0d0

      do i = 1, nmax
       ss = ss + pl(i)
      enddo

      if(idebug.eq.1)then
       if(ss.le.0.0d0) write(iw,2100) (pl(i),i=1,nmax)
 2100  format(47h *** error *** initial l distribution wrong ***/
     &        (1x,5f10.6))
      endif

      if(ss.le.0.0d0) ss = 1.0d0

      do i = 1, nmax
       pl(i) = pl(i)/ss
      enddo

*-----------------------------------------------------------------------
*     set character in m9(i,j)
*-----------------------------------------------------------------------

      if(idebug.eq.1)then
       if(mpu.gt.1 .and. ipn.eq.0) write(ipunch,2300) mpu, ide
 2300  format(14(1h*),19h new fit -- at most,i4,15h lines punched ,
     &        20(1h*),i5)
      endif

      do i = 1, 4
       do j = 1, 7
        m9(i,j) = l9
       enddo
      enddo


      if ( k0 .ne. 0 ) then
       do i = 1, k0
        j = nn0(i)
        m9(1,i) = k9(j+1,1)
       enddo
      endif

      if ( k1 .ne. 0 ) then
       do i = 1, k1
        j = nn1(i)
        k = m1(i)
        m9(2,i) = k9(j+1,k+1)
       enddo
      endif

      if ( k2 .ne. 0 ) then
       do i = 1, k2
        j = nn2(i)
        k = m2(i)
        m9(3,i) = k9(j+1,k+1)
       enddo
      endif

      if ( k3 .ne. 0 ) then
       do i = 1, k3
        j = nn3(i)
        k = m3(i)
        m9(4,i) = k9(j+1,k+1)
       enddo
      endif

*-----------------------------------------------------------------------
*     set integer number of atomic number, z
*     and punch identification number, ide
*-----------------------------------------------------------------------

      iz9 = idint(z+0.5d0)
      iz1 = mod(iz9,10) + 1
      iz2 = iz9/10 + 1
      id1 = mod(ide,10) + 1
      id2 = mod(ide/10,10) + 1
      id3 = mod(ide/100,10) + 1
      id4 = mod(ide/1000,10) + 1
      id5 = ide/10000 + 1

*-----------------------------------------------------------------------
*     set l-distribution 
*-----------------------------------------------------------------------

      if( ip8 .ne. 0 ) then

       ss = 0.0d0
       nu = nmax*(nmax+1)/2

       do i = 1, nu
        ss = ss + pln(i)
       enddo

       if(idebug.eq.1)then
        if(ss.le.0.0d0) write(iw,2100)
       endif
       if(ss.le.0.0d0) ss = 1.0d0

       do i = 1, nu
        pln(i) = pln(i)/ss 
       enddo

      endif

*-----------------------------------------------------------------------
*     for printing input echo in output file
*     it does not use in muonab routine, thus skip to 9999
*-----------------------------------------------------------------------

      if(idebug.ne.1) goto 9999

      id9(1) = iz2
      id9(2) = iz1

      write(iw,4000) i9(1)

      do i = 2, 4
       write(iw,4010) i9(i)
      enddo

      do i = 1, 4
       write(iw,4100) i9(i)
      enddo

      do i = 1, 9
       do j = 1, 2
        do k = 1, 4
         j4(i,j,k) = jb
         i99 = id9(j)
         if( mod(j9(1,i99)/2**(9-i),2) .eq. 1 ) j4(i,j,k) = jx(k)
        enddo
       enddo
      enddo

      write(iw,4200) i9(1), ((j4(i,j,1),i=1,9),j=1,2)

      do k = 2, 4
       write(iw,4210) i9(k), ((j4(i,j,k),i=1,9),j=1,2)
      enddo

      do i = 1, 9
       do j = 1, 2
        do k = 1, 4
         j4(i,j,k) = jb
         i99 = id9(j)
         if ( mod(j9(2,i99)/2**(9-i),2) .eq. 1 ) j4(i,j,k) = jx(k)
        enddo
       enddo
      enddo

      write(iw,4300) i9(1), ((j4(i,j,1),i=1,9),j=1,2)

      do k = 2, 4
       write(iw,4310) i9(k), ((j4(i,j,k),i=1,9),j=1,2)
      enddo

      do i = 1, 9
       do j = 1, 2
        do k = 1, 4
         j4(i,j,k) = jb
         i99 = id9(j)
         if ( mod(j9(3,i99)/2**(9-i),2) .eq. 1 ) j4(i,j,k) = jx(k)
        enddo
       enddo
      enddo

      write(iw,4400) i9(1), ((j4(i,j,1),i=1,9),j=1,2)

      do k = 2, 4
       write(iw,4410)i9(k),((j4(i,j,k),i=1,9),j=1,2)
      enddo

      do i = 1, 9
       do j = 1, 2
        do k = 1, 4
         j4(i,j,k) = jb
         i99 = id9(j)
         if ( mod(j9(4,i99)/2**(9-i),2) .eq. 1 ) j4(i,j,k) = jx(k)
        enddo
       enddo
      enddo

      write(iw,4500) i9(1), ((j4(i,j,1),i=1,9),j=1,2), nmax, nopt

      do k = 2, 4
       write(iw,4510) i9(k), ((j4(i,j,k),i=1,9),j=1,2), nmax, nopt
      enddo

      do i = 1, 9
       do j = 1, 2
        do k = 1, 4
         j4(i,j,k) = jb
         i99 = id9(j)
         if ( mod(j9(5,i99)/2**(9-i),2) .eq. 1 ) j4(i,j,k) = jx(k)
        enddo
       enddo
      enddo

      write(iw,4600) i9(1), ((j4(i,j,1),i=1,9),j=1,2), alexp, cl1, cl2

      do k = 2, 4
       write(iw,4610) i9(k), ((j4(i,j,k),i=1,9),j=1,2), alexp, cl1, cl2
      enddo

      do i = 1, 9
       do j = 1, 2
        do k = 1, 4
         j4(i,j,k) = jb
         i99 = id9(j)
         if ( mod(j9(6,i99)/2**(9-i),2) .eq. 1 ) j4(i,j,k) = jx(k)
        enddo
       enddo
      enddo

      if ( ip8 .eq. 0 ) then 
       write(iw,4700) i9(1), ((j4(i,j,1),i=1,9),j=1,2), (pl(i),i=1,10)
      else 
       write(iw,4720) i9(1), ((j4(i,j,1),i=1,9),j=1,2)
      endif

      do k = 2, 4
       if ( ip8 .eq. 0 ) then
        write(iw,4710) i9(k), ((j4(i,j,k),i=1,9),j=1,2), (pl(i),i=1,10)
       else
        write(iw,4730) i9(k), ((j4(i,j,k),i=1,9),j=1,2)
       endif
      enddo

      do i = 1, 9
       do j = 1, 2
        do k = 1, 4
         j4(i,j,k) = jb
         i99 = id9(j)
         if ( mod(j9(7,i99)/2**(9-i),2) .eq. 1 ) j4(i,j,k) = jx(k)
        enddo
       enddo
      enddo

      if ( ip8 .eq. 0 ) then
       write(iw,4800) i9(1), ((j4(i,j,1),i=1,9),j=1,2), (pl(i),i=11,20)
      else
       write(iw,4820) i9(1), ((j4(i,j,1),i=1,9),j=1,2)
      endif

      do k = 2, 4
       if ( ip8 .eq. 0 ) then
        write(iw,4810) i9(k), ((j4(i,j,k),i=1,9),j=1,2), (pl(i),i=11,20)
       else
        write(iw,4830) i9(k), ((j4(i,j,k),i=1,9),j=1,2)
       endif
      enddo

      do i = 1, 9
       do j = 1, 2
        do k = 1, 4
         j4(i,j,k) = jb
         i99 = id9(j)
         if ( mod(j9(8,i99)/2**(9-i),2) .eq. 1 ) j4(i,j,k) = jx(k)
        enddo
       enddo
      enddo

      write(iw,4900) i9(1), ((j4(i,j,1),i=1,9),j=1,2)

      do k = 2, 4
       write(iw,4910) i9(k), ((j4(i,j,k),i=1,9),j=1,2)
      enddo

      do i = 1, 9
       do j = 1, 2
        do k = 1, 4
         j4(i,j,k) = jb
         i99 = id9(j)
         if ( mod(j9(9,i99)/2**(9-i),2) .eq. 1 ) j4(i,j,k) = jx(k)
        enddo
       enddo
      enddo

      write(iw,5000) i9(1), ((j4(i,j,1),i=1,9),j=1,2)

      do k = 2, 4
       write(iw,5010) i9(k), ((j4(i,j,k),i=1,9),j=1,2)
      enddo

      write(iw,5100) i9(1)

      do k = 2, 4
       write(iw,5110) i9(k)
      enddo

      write(iw,5200) i9(1), zsk, zsl, zsm, (pop(i),i=1,3)

      do k = 2, 4
       write(iw,5210) i9(k), zsk, zsl, zsm, (pop(i),i=1,3)
      enddo

      write(iw,5300) i9(1), (pop(i),i=4,6)

      do k = 2, 4
       write(iw,5310) i9(k), (pop(i),i=4,6)
      enddo

      write(iw,5400) i9(1)

      do k = 2, 4
       write(iw,5410) i9(k)
      enddo

      write(iw,5500) i9(1), ipc, widthk

      do k = 2, 4
       write(iw,5510) i9(k), ipc, widthk
      enddo

      write(iw,5600) i9(1), be

      do k = 2, 4
       write(iw,5610) i9(k), be
      enddo

      write(iw,5700) i9(1)

      do k = 2, 4
       write(iw,5710) i9(k)
      enddo

*-----------------------------------------------------------------------
*     tfm
*-----------------------------------------------------------------------

      a9 = a
      if ( amassa .ne. 0.0d0 ) a9 = amassa
      tqm = dabs(tfm-2.3001d0)

      if ( tqm .ge. 1.0d-10 ) then
       write(iw,5800) i9(1), a9, cfm, tfm
      else
       write(iw,5820) i9(1), a9
      endif

      do k = 2, 4
       if( tqm .ge. 1.0d-10 ) then
        write(iw,5810) i9(k), a9, cfm, tfm
       else
        write(iw,5830) i9(k), a9
       endif
      enddo

      if ( step+rmatch .gt. 1.0d-20 ) then
       write(iw,5900) i9(1),step,rmatch
      else
       write(iw,5920) i9(1)
      endif

      do k = 2, 4
       if(step+rmatch.gt.1.0d-20) then
        write(iw,5910) i9(k), step, rmatch
       else
        write(iw,5930) i9(k)
       endif
      enddo

*-----------------------------------------------------------------------

      write(iw,6000) i9(1), amassm, amasse, amassn
      do k = 2, 4
       write(iw,6010) i9(k), amassm, amasse, amassn
      enddo

      write(iw,6100) i9(1), d2p1s, esp

      do k = 2, 4
       write(iw,6110) i9(k), d2p1s, esp
      enddo

      write(iw,6200) i9(1)

      do k = 2, 4
       write(iw,6210) i9(k)
      enddo

      write(iw,6300) i9(1), ehigh, climit

      do k = 2, 4
       write(iw,6310) i9(k), ehigh, climit
      enddo

      write(iw,6400) i9(1), elow, icc

      do k = 2, 4
       write(iw,6410) i9(k), elow, icc
      enddo

*-----------------------------------------------------------------------

      eab = (ea-99.0d0)**2 + (eb-99.0d0)**2

      if ( eab .gt. 1.0d-20 ) then
       write(iw,6500) i9(1), eres, ea, eb
      else
       write(iw,6520) i9(1), eres
      endif

      do k = 2, 4
       if ( eab .gt. 1.0d-20 ) then
        write(iw,6510) i9(k), eres, ea, eb
       else
        write(iw,6530) i9(k), eres
       endif
      enddo

      write(iw,6600) i9(1), cd

      do k = 2, 4
       write(iw,6610) i9(k), cd
      enddo

      write(iw,6700) i9(1), npol, ipol

      do k = 2, 4
       write(iw,6710) i9(k), npol, ipol
      enddo

      do k = 1, 4
       write(iw,6800) i9(k)
      enddo

      do k = 1, 4
       write(iw,6900) i9(k)
      enddo

      do k = 1, 4
       write(iw,7000) i9(k)
      enddo

      write(iw,7100) i9(1)

      do k = 2, 4
       write(iw,7110) i9(k)
      enddo

      write(iw,7200) i9(1)

      do k = 2, 4
       write(iw,7210) i9(k)
      enddo

      write(iw,7300) i9(1), k0, k1, k2, k3, iread, iw, ipunch

      do k = 2, 4
       write(iw,7310) i9(k), k0, k1, k2, k3, iread, iw, ipunch
      enddo

      write(iw,7400) i9(1)

      do k = 2, 4
       write(iw,7410) i9(k)
      enddo

      do k = 1, 4
       write(iw,7500) i9(k)
      enddo

      write(iw,7600) i9(1), iprint

      do k = 2, 4
       write(iw,7610) i9(k), iprint
      enddo

      write(iw,7700) i9(1), (m9(1,j),j=1,3), idb

      do k = 2, 4
       write(iw,7710) i9(k), (m9(1,j),j=1,3), idb
      enddo

      write(iw,7800) i9(1), (m9(2,j),j=1,7), ic

      do k = 2, 4
       write(iw,7810) i9(k), (m9(2,j),j=1,7), ic
      enddo

      write(iw,7900) i9(1), (m9(3,j),j=1,7)

      do k = 2, 4
       write(iw,7910) i9(k), (m9(3,j),j=1,7)
      enddo

      write(iw,8000) i9(1), (m9(4,j),j=1,7)
      do k = 2, 4
       write(iw,8010) i9(k), (m9(4,j),j=1,7)
      enddo

      write(iw,8100) i9(1), ffactd
      do k = 2, 4
       write(iw,8110) i9(k), ffactd
      enddo

*-----------------------------------------------------------------------

      mdir = nmax**2

      write(iw,8200) i9(1), idir, mdir

      do k = 2, 4
       write(iw,8210) i9(k), idir, mdir
      enddo

      write(iw,8300) i9(1), mpu

      do k = 2, 4
       write(iw,8310) i9(k), mpu
      enddo

      write(iw,8400) i9(1), ipn

      do k = 2, 4
       write(iw,8410) i9(k), ipn
      enddo

      write(iw,8500) i9(1)

      do k = 2, 4
       write(iw,8510) i9(k)
      enddo

      do k = 1, 4
       write(iw,8600) i9(k)
      enddo

      write(iw,8700) i9(1)

      do k = 2, 4
       write(iw,8710) i9(k)
      enddo

      write(iw,8800) i9(1)

      do k = 2, 4
       write(iw,8810) i9(k)
      enddo

*-----------------------------------------------------------------------

      id9(1) = id5
      id9(2) = id4
      id9(3) = id3
      id9(4) = id2
      id9(5) = id1

      do i = 1, 9
       do j = 1, 5
        do k = 1, 4
         j4(i,j,k) = jb
         i99 = id9(j)
         if ( mod(j9(1,i99)/2**(9-i),2) .eq. 1 ) j4(i,j,k) = jx(k)
        enddo
       enddo
      enddo

      write(iw,8900) i9(1), ((j4(i,j,1),i=1,9),j=1,5)

      do k = 2, 4
       write(iw,8910) i9(k), ((j4(i,j,k),i=1,9),j=1,5)
      enddo

      do i = 1, 9
       do j = 1, 5
        do k = 1, 4
         j4(i,j,k) = jb
         i99 = id9(j)
         if ( mod(j9(2,i99)/2**(9-i),2) .eq. 1 ) j4(i,j,k) = jx(k)
        enddo
       enddo
      enddo

      write(iw,9000) i9(1), ((j4(i,j,1),i=1,9),j=1,5)

      do k = 2, 4
       write(iw,9010) i9(k), ((j4(i,j,k),i=1,9),j=1,5)
      enddo

      do i = 1, 9
       do j = 1, 5
        do k = 1, 4
         j4(i,j,k) = jb
         i99 = id9(j)
         if ( mod(j9(3,i99)/2**(9-i),2) .eq. 1) j4(i,j,k) = jx(k)
        enddo
       enddo
      enddo

      write(iw,9100) i9(1), jtm, ((j4(i,j,1),i=1,9),j=1,5)

      do k = 2, 4
       write(iw,9110) i9(k), jtm, ((j4(i,j,k),i=1,9),j=1,5)
      enddo

      do i = 1, 9
       do j = 1, 5
        do k = 1, 4
         j4(i,j,k) = jb
         i99 = id9(j)
         if ( mod(j9(4,i99)/2**(9-i),2) .eq. 1 ) j4(i,j,k) = jx(k)
        enddo
       enddo
      enddo

      write(iw,9200) i9(1), jtd, ip1, ((j4(i,j,1),i=1,9),j=1,5)

      do k = 2, 4
       write(iw,9210) i9(k), jtd, ip1, ((j4(i,j,k),i=1,9),j=1,5)
      enddo

      do i = 1, 9
       do j = 1, 5
        do k = 1, 4
         j4(i,j,k) = jb
         i99 = id9(j)
         if ( mod(j9(5,i99)/2**(9-i),2) .eq. 1 ) j4(i,j,k) = jx(k)
        enddo
       enddo
      enddo

      write(iw,9300) i9(1), jtq, ip2, ((j4(i,j,1),i=1,9),j=1,5)

      do k = 2, 4
       write(iw,9310) i9(k), jtq, ip2, ((j4(i,j,k),i=1,9),j=1,5)
      enddo

      do i = 1, 9
       do j = 1, 5
        do k = 1, 4
         j4(i,j,k) = jb
         i99 = id9(j)
         if ( mod(j9(6,i99)/2**(9-i),2) .eq. 1 ) j4(i,j,k) = jx(k)
        enddo
       enddo
      enddo

      write(iw,9400) i9(1), jto, ip3, ((j4(i,j,1),i=1,9),j=1,5)

      do k = 2, 4
       write(iw,9410) i9(k), jto, ip3, ((j4(i,j,k),i=1,9),j=1,5)
      enddo

      do i = 1, 9
       do j = 1, 5
        do k = 1, 4
         j4(i,j,k) = jb
         i99 = id9(j)
         if ( mod(j9(7,i99)/2**(9-i),2) .eq. 1 ) j4(i,j,k) = jx(k)
        enddo
       enddo
      enddo

      write(iw,9500) i9(1), ((j4(i,j,1),i=1,9),j=1,5)

      do k = 2, 4
       write(iw,9510) i9(k), ((j4(i,j,k),i=1,9),j=1,5)
      enddo

      do i = 1, 9
       do j = 1, 5
        do k = 1, 4
         j4(i,j,k) = jb
         i99 = id9(j)
         if ( mod(j9(8,i99)/2**(9-i),2) .eq. 1 ) j4(i,j,k) = jx(k)
        enddo
       enddo
      enddo

      do k = 1, 4
       write(iw,9600) i9(k), ((j4(i,j,k),i=1,9),j=1,5)
      enddo

      do i = 1, 9
       do j = 1, 5
        do k = 1, 4
         j4(i,j,k) = jb
         i99 = id9(j)
         if ( mod(j9(9,i99)/2**(9-i),2) .eq. 1) j4(i,j,k) = jx(k)
        enddo
       enddo
      enddo

      write(iw,9700) i9(1), yc, ((j4(i,j,1),i=1,9),j=1,5)

      do k = 2, 4
       write(iw,9710)i9(k),yc,((j4(i,j,k),i=1,9),j=1,5)
      enddo

      do k = 1, 4
       write(iw,9800) i9(k)
      enddo

      do k = 1, 4
       write(iw,9900)i9(k)
      enddo

      write(iw,9910)

*-----------------------------------------------------------------------

 4000 format(1h1/a1,120(1h*))
 4010 format(a1,120(1h*))
 4100 format(a1,1h*,22x,1hi,95x,1h*)
 4200 format(a1,2h* ,9a1,2x,9a1,34h i table of all input parameters -,
     &60h-- defaults (if appropriate) follow the values, in parenthes,
     &4hes *)                                                         
 4210 format(a1,1h*,1x,9a1,2x,9a1,1x,1hi,95x,1h*)
 4300 format(a1,2h* ,9a1,2x,9a1,2h i,95(1h.),1h*)
 4310 format(a1,1h*,1x,9a1,2x,9a1,1x,1hi,95x,1h*)
 4400 format(a1,2h* ,9a1,2x,9a1,2h i,95x,1h*)
 4410 format(a1,1h*,1x,9a1,2x,9a1,1x,1hi,95x,1h*)
 4500 format(a1,2h* ,9a1,2x,9a1,15h i e12 init. n=,i2,12h(max=20,def=,
     &23h15) e11 l-dist. option=,i2,30h(def=0 /-1=inputed,0=statist.,,
     &14h2=quadratic) *)
 4510 format(a1,1h*,1x,9a1,2x,9a1,1x,1hi,13x,i2,35x,i2,43x,1h*)
 4600 format(a1)
 4610 format(a1,2h* ,9a1,2x,9a1,2h i,22x,f7.5,25x,f7.5,8x,f7.5,19x,
     &1h*)                                                            
 4700 format(a1,2h* ,9a1,2x,9a1,17h i e13 norm.init/,f7.6,9(1x,f7.6),
     &2h *)
 4710 format(a1,2h* ,9a1,2x,9a1,2h i,15x,10(f7.6,1x),1h*)
 4720 format(a1,2h* ,9a1,2x,9a1,7h i e13 ,24h  l-distribution extends,
     &60h beyond starting n.  see next page for complete distribution,
     &1h.,5x,1h*)
 4730 format(a1,1h*,1x,9a1,2x,9a1,2h i,95x,1h*)
 4800 format(a1,2h* ,9a1,2x,9a1,17h i l-dist.(0-19)/,10(f7.6,1x),1h*)
 4810 format(a1,2h* ,9a1,2x,9a1,2h i,15x,10(f7.6,1x),1h*)
 4820 format(a1,2h* ,9a1,2x,9a1,3h i ,4x,8(9x,1h*),10x,1h*)
 4830 format(a1,2h* ,9a1,2x,9a1,2h i,95x,1h*)
 4900 format(a1,2h* ,9a1,2x,9a1,2h i,95(1h.),1h*)
 4910 format(a1,2h* ,9a1,2x,9a1,2h i,95x,1h*)
 5000 format(a1,2h* ,9a1,2x,9a1,2h i,45x,1hi,49x,1h*)
 5010 format(a1,2h* ,9a1,2x,9a1,2h i,45x,1hi,49x,1h*)
 5100 format(a1,1h*,22x,42hi e07 effective charge for electronic shel,
     &55hls  i e08  popul. of el. subshells (fraction of full) *)
 5110 format(a1,1h*,22x,1hi,45x,1hi,49x,1h*)
 5200 format(a1,27h*  e01  z /  (needed)  i k/,f6.3,4h, l/,f6.3,4h, m/
     &,f6.3,21h  (all needed)  i 1s=,f5.3,11h(1.000) 2s=,f5.3,
     &11h(1.000) 2p=,f5.3,9h(1.000) *)
 5210 format(a1,1h*,22x,1hi,3x,f6.3,4x,f6.3,4x,f6.3,16x,1hi,4x,f5.3,
     &11x,f5.3,11x,f5.3,8x,1h*)
 5300 format(a1,1h*,22(1h-),1h+,17(1h-),1h+,27(1h-),5h+ 3s=,f5.3,
     &11h(1.000) 3p=,f5.3,11h(1.000) 3d=,f5.3,9h(1.000) *)
 5310 format(a1,1h*,22x,1h+,17x,1h+,27x,1h+,4x,f5.3,11x,f5.3,11x,f5.3,
     &8x,1h*)
 5400 format(a1,50h*   e10 depletion of electronic shells   i e09 ele,
     &20hct. 1s width in eV +,49(1h-),1h*)
 5410 format(a1,1h*,40x,1hi,27x,1h+,49x,1h*)
 5500 format(a1,1h*,5x,2hk/,i1,3h(0),6x,2hl/,i1,3h(0),6x,2hm/,i1,3h(0)
     &,5x,1hi,6x,f7.3,9h(000.000),5x,29hi e02 av. el. binding ener. f,
     &22hor atom z-1 (needed) *)
 5510 format(a1,1h*,7x,i1,11x,i1,11x,i1,8x,1hi,6x,f7.3,14x,1hi,49x,1h*
     &)

 5600 format(a1,50h* (0=yes,1=no - if k/0,1s width is used) i (exper.,
     &23h or inputed value) i k/,f9.2,8h(eV)  l/,f8.2,8h(ev)  m/,f8.2,
     &6h(eV) *)
 5610 format(a1,1h*,40x,1hi,27x,1hi,3x,f9.2,8x,f8.2,8x,f8.2,5x,1h*)
 5700 format(a1,1h*,40(1h-),1h+,27(1h-),1h+,49(1h-),1h*)
 5710 format(a1,1h*,40x,1h+,27x,1h+,49x,1h*)
 5800 format(a1,20h* e19 atomic weight=,f6.2,21h(140.00) e23 fermi pa,
     &34hrameters (not used in program)  c=,f7.5,17h(fm),skin thick. ,
     &2ht=,f7.5,6h(fm) *)
 5810 format(a1,1h*,19x,f6.2,55x,f7.5,19x,f7.5,5x,1h*)
 5820 format(a1,20h* e19 atomic weight=,f6.2,21h(140.00) e23 fermi pa,
     &60hrameters (not used in program)    * *  n o t   s p e c i f i,
     &13h e d  * *   *)
 5830 format(a1,1h*,19x,f6.2,93x,1h*)
 5900 format(a1,50h* e24 /dirac/ program parameters(not used)/ step i,
     &14hn integration=,1pe9.3,17h matching radius=,e9.3,
     &21h(fm)(for reference) *)
 5910 format(a1,1h*,63x,1pe9.3,17x,e9.3,20x,1h*)
 5920 format(a1,50h* e24 /dirac/ program parameters(not used)/   *  *,
     &46h    n  o  t     s  p  e  c  i  f  i  e  d    *,8(2x,1h*))
 5930 format(a1,1h*,118x,1h*)
 6000 format(a1,23h* e30 masses/ particle=,f9.4,18h(206.7686)(elec. m,
     &17hasses)  electron=,f8.1,24h(511003.4)(eV)  nucleon=,f6.2,
     &15h(931.48)(MeV) *)
 6010 format(a1,1h*,22x,f9.4,35x,f8.1,24x,f6.2,14x,1h*)
 6100 format(a1)
 6110 format(a1)
 6200 format(a1,1h*,12(1h-),1h+,31(1h-),1h+,55(1h-),1h+,17(1h-),1h*)
 6210 format(a1,1h*,12x,1h+,31x,1h+,55x,1h+,17x,1h*)
 6300 format(a1,1h*,12x,14h/ e20 hi. cut=,f6.3,17h(20.000)MeV i e21,
     &16h intens. cutoff=,1pe9.3,33h(1.000e-06)(per particle) i e33 s,
     &12htar option *)
 6310 format(a1,1h*,12x,1h/,13x,f6.3,12x,1hi,20x,1pe9.3,26x,1hi,17x,
     &1h*)
 6400 format(a1)
 6410 format(a1,1h*,12x,1h/,13x,f6.3,12x,1hi,55x,1hi,12x,i1,4x,1h*)
 6500 format(a1,27h* parameters / e22 resol. =,f6.5,13h(.00030)MeV i,
     &5x,2ha=,f7.2,11h (keV) , b=,f7.3,27h  channel no =(e-a)/b) i 0=,
     &15hdef, 1=readin *)
 6510 format(a1,1h*,12x,1h/,13x,f6.5,12x,1hi,7x,f7.2,11x,f7.3,23x,
     &1hi,17x,1h*)
 6520 format(a1,27h* parameters / e22 resol. =,f6.5,13h(.00030)MeV i,
     &5x,57h  *  *   n o t   s p e c i f i e d   *  *  *  *   i 0=def,
     &12h, 1=readin *)
 6530 format(a1,1h*,12x,1h/,13x,f6.5,12x,1hi,55x,1hi,17x,1h*)
 6600 format(a1,1h*,12x,24h/ e34 intensities /  5*=,1pe7.1,8h(.1) 4*=,
     &e7.1,9h(.01) 3*=,e7.1,10h(.001) 2*=,e7.1,11h(.0001) 1*=,e7.1,
     &10h(.00001) *)
 6610 format(a1,1h*,12x,1h/,23x,1pe7.1,8x,e7.1,9x,e7.1,10x,e7.1,11x,
     &e7.1,9x,1h*)
 6700 format(a1,18h* e17,e18 quan.dep,4(1h/,i2,1h,,i2,1h,,i2,1h,,i2,
     11h,,i2),28h(all/-1=start ran) do depol=,i1,13h(0/0=y,1=n) *)
 6710 format(a1,1h*,17x,20(1x,i2),28x,i1,12x,1h*)
 6800 format(a1,1h*,118x,1h*)
 6900 format(a1,120(1h*))
 7000 format(a1,1h*,62x,1hi,55x,1h*)
 7100 format(a1,50h* s h e l l  a n d  s u b s h e l l   c o m b i n ,
     &60ha t i o n s  i      b o o k k e e p i n g    p a r a m e t e,
     &10h r s     *)
 7110 format(a1,1h*,62x,1hi,55x,1h*)
 7200 format(a1,1h*,62x,1hi,55x,1h*)
 7210 format(a1,1h*,62x,1hi,55x,1h*)
 7300 format(a1,16h* e38 cases/  m/,i1,16h(max=3,def=3),d/,i1,3h,q/,
     &i1,3h,o/,i1,48h (d,q,o/max=7,def=4) i e31 logical unit no.s rea,
     &2hd=,i1,10h(5),write=,i1,10h(6),punch=,i1,5h(7) *)
 7310 format(a1,1h*,15x,i1,16x,i1,3x,i1,3x,i1,21x,1hi,28x,i1,10x,i1,
     &10x,i1,4x,1h*)
 7400 format(a1,1h*,62(1h.),1hi,55(1h.),1h*)
 7410 format(a1,1h*,62x,1hi,55x,1h*)
 7500 format(a1,1h*,62x,1hi,55x,1h*)
 7600 format(a1,8h* e39,40,1x,2hc1,6x,2hc2,6x,2hc3,6x,2hc4,6x,2hc5,6x,
     &2hc6,6x,2hc7,4x,20hi e35 print option =,i2,17h(0/print all) sel,
     &18hect codes 0 - 63 *)
 7610 format(a1,1h*,62x,1hi,19x,i2,34x,1h*)
 7700 format(a1,7h* mon/ ,a2,6h(1s)  ,a2,6h(2t)  ,a2,13h(3t)  ------ ,
     &46h ------  ------  ------  i e36 debug option = ,i1,
     &35h(0)(0=n,1=y) deb prints all rates *)
 7710 format(a1,1h*,6x,a2,6x,a2,6x,a2,38x,1hi,20x,i1,34x,1h*)
 7800 format(a1)
 7810 format(a1,1h*,6x,7(a2,6x),1hi,20x,i1,34x,1h*)
 7900 format(a1)
 7910 format(a1,1h*,6x,7(a2,6x),1hi,55x,1h*)
 8000 format(a1)
 8010 format(a1,1h*,6x,7(a2,6x),1hi,55x,1h*)
 8100 format(a1)
 8110 format(a1,1h*,62x,1hi,16x,f5.2,34x,1h*)
 8200 format(a1)
 8210 format(a1,1h*,62x,1hi,35x,i3,8x,i3,6x,1h*)
 8300 format(a1)
 8310 format(a1,63(1h*),1hi,36x,i3,16x,1h*)
 8400 format(a1,1h*,62x,39hi e29 punch specified transitions /    ,i1,
     &17h(0 =yes, 1 =no) *)
 8410 format(a1,1h*,62x,1hi,38x,i1,16x,1h*)
 8500 format(a1,1h*,6x,43hs e l e c t i o n    o f    p e n e t r a t,
     &6h i o n,7x,1hi,55(1h.),1h*)
 8510 format(a1,1h*,62x,1hi,55x,1h*)
 8600 format(a1,1h*,62x,1hi,55x,1h*)
 8700 format(a1,50h* e42 max number of terms in i e41 penetration sel,
     &60hection codes i e28 punched card identity no. in col.s 73-78 ,
     &10h/(10000) *)
 8710 format(a1,1h*,28x,1hi,33x,1hi,55x,1h*)
 8800 format(a1,50h* penetration (max=3,def.=1) i     (1=y,0=n *=must,
     &14h be 0,def=1) i,55x,1h*)
 8810 format(a1,1h*,28x,1hi,33x,1hi,55x,1h*)
 8900 format(a1,1h*,28x,17hi cases as in e26,17x,1hi,5(1x,9a1,1x),1h*)
 8910 format(a1,1h*,28x,1hi,33x,1hi,5(1x,9a1,1x),1h*)
 9000 format(a1,50h*    1s  2s  2p  3s  3p  3d  i      c1  c2  c3  c4,
     &14h  c5  c6  c7 i,5(1x,9a1,1x),1h*)
 9010 format(a1,1h*,28x,1hi,33x,1hi,5(1x,9a1,1x),1h*)
 9100 format(a1,4h* m/,1x,6(1x,i1,2x),5hi  m/,5h  -  ,6(2x,1h-,1x),1hi
     &,5(1x,9a1,1x),1h*)
 9110 format(a1,1h*,4x,6(1x,i1,2x),1hi,33x,1hi,5(1x,9a1,1x),1h*)
 9200 format(a1,5h* d/ ,6(1x,i1,2x),7hi  d/  ,i1,1h*,6(3x,i1),2h i,
     &5(1x,9a1,1x),1h*)
 9210 format(a1,1h*,4x,6(1x,i1,2x),1hi,6x,i1,1x,6(3x,i1),2h i,
     &5(1x,9a1,1x),1h*)
 9300 format(a1,5h* q/ ,6(1x,i1,2x),7hi  q/  ,i1,1h*,6(3x,i1),2h i,
     &5(1x,9a1,1x),1h*)
 9310 format(a1,1h*,4x,6(1x,i1,2x),1hi,6x,i1,1x,6(3x,i1),2h i,
     &5(1x,9a1,1x),1h*)
 9400 format(a1,5h* o/ ,6(1x,i1,2x),7hi  o/  ,i1,1h*,6(3x,i1),2h i,
     &5(1x,9a1,1x),1h*)
 9410 format(a1,1h*,4x,6(1x,i1,2x),1hi,6x,i1,1x,6(3x,i1),2h i,
     &5(1x,9a1,1x),1h*)
 9500 format(a1,1h*,28(1h.),1hi,33(1h.),1hi,5(1x,9a1,1x),1h*)
 9510 format(a1,1h*,28x,1hi,33x,1hi,5(1x,9a1,1x),1h*)
 9600 format(a1,1h*,62x,1hi,5(1x,9a1,1x),1h*)
 9700 format(a1,17h* e43 y-cutoff m/,f5.2,3h d/,f5.2,3h q/,f5.2,3h o/,
     &f5.2,18h (all def/ 1.00) i,5(1x,9a1,1x),1h*)
 9710 format(a1,1h*,16x,f5.2,3(3x,f5.2),17x,1hi,5(1x,9a1,1x),1h*)
 9800 format(a1,1h*,62x,1hi,55x,1h*)
 9900 format(a1,120(1h*))
 9910 format(1h1)

*-----------------------------------------------------------------------

 9999 continue

      if( ip8 .eq. 0 ) goto 11100

*-----------------------------------------------------------------------

      if(idebug.eq.1)then
       write(iw,10000)
       write(iw,10100)
       write(iw,10200)
      endif

      do i = 1, nmax

       n = nmax + 1 -i
       ss = 0.0d0

       do j = 1, n
        k = n*(n-1)/2 + j
        plx(j) = pln(k)
        ss = ss + plx(j)
        k = idint(plx(j)*1.0d6+0.5d0)
        ixx(1,j) = mod(k/100,10)
        ixx(2,j) = mod(k/10,10)
        ixx(3,j) = mod(k,10)
        if(dmod(1000.0d0*plx(j),1.0d0).gt..499999d0)
     &   plx(j)=plx(j)-.000499999d0
        if(plx(j).lt.0.0d0) plx(j) = 0.0d0
       enddo

       if(idebug.eq.1)then

        if ( n .gt. 2 ) then
         write(iw,10300) n, ss, (plx(j),j=1,n)
         write(iw,10400) ((ixx(ii,j),ii=1,3),j=1,n)
        elseif ( n .eq. 2 ) then
         write(iw,10300) n, ss, (plx(j),j=1,n)
         write(iw,10700) ((ixx(ii,j),ii=1,3),j=1,2)
        elseif ( n .eq. 1 ) then
         write(iw,10900) n, ss, plx(1)
         write(iw,11000) (ixx(ii,1),ii=1,3)
        endif

        if ( n .gt. 15 ) then
         write(iw,10450)
        elseif ( n .gt. 10 .and. n .le. 15 ) then
         write(iw,10460)
        elseif ( n .gt. 10 .and. n .le. 15 ) then
         write(iw,10470)
        elseif ( n .gt. 10 .and. n .le. 15 ) then
         write(iw,10480)
        elseif ( n .eq. 3 ) then
         write(iw,10500)
        elseif ( n .eq. 2 ) then
         write(iw,10800)
        endif

       endif

      enddo

*-----------------------------------------------------------------------

11100 call aamacascad
      call aamasort

cabe add
      esp0 = esp
      nmax0 = nmax
      do i = 1, 20
       do j = 1, 40
        energy0(i,j) = energya(i,j)
       enddo
      enddo
cabe end

*-----------------------------------------------------------------------

10000 format(20x,49h*  *   n o r m a l i z e d   i n i t i a l   l - ,
     &30hd i s t r i b u t i o n   *  */)
10100 format(7x,50htotal  i  l=0  l=1  l=2  l=3  l=4 i  l=5  l=6  l=7,
     &60h  l=8  l=9 i l=10 l=11 l=12 l=13 l=14 i l=15 l=16 l=17 l=18 ,
     &4hl=19)
10200 format(1x,13(1h-),1h+,26(1h-),1h+,26(1h-),1h+,26(1h-),1h+,26(1h-))
10300 format(3h n=,i2,1x,f7.5,4(2h i,1x,f4.3,1x,f4.3,1x,f4.3,1x,f4.3,
     &1x,f4.3))
10400 format(13x,4(2h i,2x,3i1,2x,3i1,2x,3i1,2x,3i1,2x,3i1))
10450 format(14x,1hi,26x,1hi,26x,1hi,26x,1hi)
10460 format(14x,1hi,26x,1hi,26x,1hi)
10470 format(14x,1hi,26x,1hi)
10480 format(14x,1hi)
10500 format(14x,1hi,56x,1h+,32(1h-),1h+)
10600 format(3h n=,i2,1x,f7.5,3h i ,f4.3,1x,f4.3,46x,1hi,32x,1hi)
10700 format(14x,1hi,2x,3i1,2x,3i1,46x,27hi entries folded  .abc = .a,
     &7hbcdef i)
10800 format(14x,1hi,56x,34hi to save space    def           i)
10900 format(3h n=,i2,f8.5,3h i ,f4.3,51x,1hi,32x,1hi)
11000 format(14x,1hi,2x,3i1,51x,1h+,32(1h-),1h+/1h1)

*-----------------------------------------------------------------------
      return
      end


************************************************************************
*                                                                      *
      subroutine aamaread(iotaama)
*                                                                      *
*     reads input cards and sets parameters according to codes         *
*                                                                      *
*       last modified by S.Abe on 2018/10/30                           *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      common /aama01/ energ, econs, econst, d2p1s, d2p1sm
!$OMP THREADPRIVATE(/aama01/)
      common /aama02/ zsa(3), be(3), bem(3)
!$OMP THREADPRIVATE(/aama02/)
      common /aama03/ k0, k1, k2, k3
!$OMP THREADPRIVATE(/aama03/)
      common /aama04/ nn0(3), nn1(7), nn2(7), nn3(7)
!$OMP THREADPRIVATE(/aama04/)
      common /aama05/ ip1(7), ip2(7), ip3(7), iq1(7), iq2(7), iq3(7)
!$OMP THREADPRIVATE(/aama05/)
      common /aama06/ z, zsk, zsl, zsm, zskz, zslz, zsmz
!$OMP THREADPRIVATE(/aama06/)
      common /aama07/ amzz(3)
!$OMP THREADPRIVATE(/aama07/)
      common /aama08/ angd, aid, aidsq
!$OMP THREADPRIVATE(/aama08/)
      common /aama09/ angq, aiq, aiqsq
!$OMP THREADPRIVATE(/aama09/)
      common /aama10/ ango, aio, aiosq
!$OMP THREADPRIVATE(/aama10/)
      common /aama11/ rr(18), rau, rad, ra(4), rd(4), rsa(4)
!$OMP THREADPRIVATE(/aama11/)
      common /aama12/ pop(6)
!$OMP THREADPRIVATE(/aama12/)
      common /aama13/ e(1000), ai(1000), energya(20,40), m, ia(1000)
!$OMP THREADPRIVATE(/aama13/)
      common /aama14/ dza, dza2, dredm
!$OMP THREADPRIVATE(/aama14/)
      common /aama15/ pl(20), pln(210)
!$OMP THREADPRIVATE(/aama15/)
      common /aama16/ a, cfm
!$OMP THREADPRIVATE(/aama16/)
      common /aama17/ yc(4)
!$OMP THREADPRIVATE(/aama17/)

      common /aamac01/ iread,iw,ipunch,iprint
!$OMP THREADPRIVATE(/aamac01/)
      common /aamac02/ ffact(60),ffactd
!$OMP THREADPRIVATE(/aamac02/)
      common /aamac03/ pi,coeff,picoef,amassm,amasse,alfa,ne(3),jk(3)
!$OMP THREADPRIVATE(/aamac03/)
      common /aamac04/ amassa,amassn
!$OMP THREADPRIVATE(/aamac04/)
      common /aamac05/ m1(7),m2(7),m3(7),idb
!$OMP THREADPRIVATE(/aamac05/)
      common /aamac06/ zmk,zml,zmm,zmkm,zmlm,zmmm,cl1,cl2,ivers
!$OMP THREADPRIVATE(/aamac06/)
      common /aamac07/ jtm(6),jtd(6),jtq(6),jto(6)
!$OMP THREADPRIVATE(/aamac07/)
      common /aamac08/ coemon(30),expmon(30),ifm(6),jm(10)
!$OMP THREADPRIVATE(/aamac08/)
      common /aamac09/ coedip(42),expdip(42),coed(4),coedp(9),
     &                 ifd(9),jd(14)
!$OMP THREADPRIVATE(/aamac09/)
      common /aamac10/ coequa(45),expqua(45),coeq(11),ifq(10),jq(15)
!$OMP THREADPRIVATE(/aamac10/)
      common /aamac11/ coeoct(45),expoct(45),coeo(11),ifo(10),jo(15)
!$OMP THREADPRIVATE(/aamac11/)
      common /aamac12/ hbar,widthk,npol(20),ipol,ide,ip8,ipc(3),nmax
!$OMP THREADPRIVATE(/aamac12/)
      common /aamac13/ ehigh,elow,climit,eres,esp,espm,cd(5),ea,eb,
     &                 icc,mpu,icpu(200),ipn
!$OMP THREADPRIVATE(/aamac13/)
      common /aamac14/ tfm,step,rmatch,alexp,nopt,idir,iyc
!$OMP THREADPRIVATE(/aamac14/)
      common /aamac15/ ij(4),yj(4),jj1(4)
!$OMP THREADPRIVATE(/aamac15/)
      common /aamac16/ ic,ll(20)
!$OMP THREADPRIVATE(/aamac16/)
      common /aamac17/ elecbe(104,3)
      common /aamadbg/ idebug

*-----------------------------------------------------------------------

      character chln*200

      dimension chead(79)
      character chead*3
      data chead /'ijk', 'ene', 'ama', 'amn', 'd21',
     &            'be ', 'k  ', 'nn0', 'nn1', 'nn2',
     &            'nn3', 'ip1', 'ip2', 'ip3', 'ic ',
     &            'ire', 'iwr', 'ipu', 'fd ', 'amm',
     &            'ame', 'z  ', 'zs ', 'ifm', 'ifd',
     &            'ifq', 'ifo', 'm1 ', 'm2 ', 'm3 ',
     &            'pop', 'jm ', 'jd ', 'jq ', 'jo ',
     &            'ehi', 'elo', 'clm', 'ers', 'icc',
     &            'cd ', 'eab', 'sto', 'c  ', 'xeq',
     &            'a  ', 'ct ', 'stp', 'kwd', 'nop',
     &            'nmx', 'pl ', 'dir', 'ipr', 'npl',
     &            'ipl', 'zm ', 'zmm', 'iq1', 'iq2',
     &            'iq3', 'idb', 'yc ', 'iyc', 'ij ',
     &            'yj ', 'jj1', 'ipc', 'cl ', 'esp',
     &            'pun', 'ide', 'ipn', 'jtm', 'jtd',
     &            'jtq', 'jto', 'pln', 'ip '/
      data t8 / 0.0d0 /
      data irc / 0 /

      dimension s(70)
      character s*1, b1*3, b2*3

      character c*1

*-----------------------------------------------------------------------

  100 irc = irc+1
      read(iotaama,'(a200)') chln

      do i = 1, 200
         c = chln(i:i)
         if( c .ge. 'A' .and. c .le. 'Z' ) 
     &       c = char(ichar(c)+ichar('a')-ichar('A'))
         chln(i:i) = c
      end do

      if( chln(1:4) .eq. 'stop' ) return
      read(chln,'(a3,a3,a4,70a1)') b, b1, b2, s

      j = 0

      if(idebug.eq.1)then
       if(mod(iprint/2,2).eq.0) write(iw,300) irc, chln
  300  format(1x,14hinput card no.,i3,5x,3h---,5x,a)
      endif

      do 400 i = 1, 89
       if( chln(1:3) .eq. chead(i) ) then
         j = i
         goto 600
       endif
  400 continue

      if(idebug.eq.1)then
       write(iw,500) chln(1:3)
  500  format(/46h *** error *** input action code illegal *** =,a3,
     &         21h *** card ignored ***)
      endif
      goto 100

  600 jt = j/10 + 1
      ju = mod(j,10) + 1
      ni = 1

      goto(999,1999,2999,3999,4999,5999,6999,7999,8999), jt

  900 format(53h *** warning *** unimplemented user code for this ver,
     &29hsion (no effect produced) ***)

*-----------------------------------------------------------------------
*     read parameter 00-09 ...
*-----------------------------------------------------------------------

  999 goto (100,1100,1200,1300,1400,1500,1600,1700,1800,1900), ju

 1100 continue
      if(idebug.eq.1)then
       if(ivers.ne.1) write(iw,900)
      endif
      if(ivers.ne.1) goto 100
      call aamadcyphr(s,ni,0,dum,ijk)
      goto 100

 1200 continue
      if(idebug.eq.1)then
       if(ivers.ne.1) write(iw,900)
      endif
      if(ivers.ne.1) goto 100
      call aamadcyphr(s,ni,1,energ,idum)
      goto 100

 1300 continue
      if(idebug.eq.1)then
       if(ivers.ne.1) write(iw,900)
      endif
      if(ivers.ne.1) goto 100
      call aamadcyphr(s,ni,1,amassa,idum)
      if(idebug.eq.1)then
       write(iw,1350)
 1350  format(53h *** warning *** amassa has been provided -- if new a,
     &        51h new amassa must be given (or zero for default) ***)
      endif
      goto 100

 1400 call aamadcyphr(s,ni,1,amassn,idum)
      goto 100

 1500 call aamadcyphr(s,ni,1,d2p1s,idum)
      goto 100

 1600 do 1650 k = 1, 3
       call aamadcyphr(s,ni,1,be(k),idum)
 1650 continue
      goto 100

 1700 call aamadcyphr(s,ni,0,dum,k0)
      call aamadcyphr(s,ni,0,dum,k1)
      call aamadcyphr(s,ni,0,dum,k2)
      call aamadcyphr(s,ni,0,dum,k3)
      goto 100

 1800 do 1850 k = 1, k0
       call aamadcyphr(s,ni,0,dum,nn0(k))
 1850 continue
      goto 100

 1900 do 1950 k = 1, k1
       call aamadcyphr(s,ni,0,dum,nn1(k))
 1950 continue
      if(idebug.eq.1)then
       if(nn1(1).ne.0) write(iw,1975)
 1975  format(53h *** error *** radiarion is not computed as the first,
     &        50h item in the multipolarity *** no action taken but,
     &        18h results are bad *)
      endif
      goto 100

*-----------------------------------------------------------------------
*     read parameter 10-19 ...
*-----------------------------------------------------------------------

 1999 goto (2000,2100,2200,2300,2400,2500,2600,2700,2800,2900), ju

 2000 do 2050 k = 1, k2
       call aamadcyphr(s,ni,0,dum,nn2(k))
 2050 continue
      if(idebug.eq.1)then
       if(nn2(1).ne.0) write(iw,1975)
      endif
      goto 100

 2100 do 2150 k = 1, k3
       call aamadcyphr(s,ni,0,dum,nn3(k))
 2150 continue
      if(idebug.eq.1)then
       if(nn3(1).ne.0) write(iw,1975)
      endif
      goto 100

 2200 do 2250 k = 1, k1
       call aamadcyphr(s,ni,0,dum,ip1(k))
 2250 continue
      if(idebug.eq.1)then
       if(ip1(1).ne.0) write(iw,2275)
 2275  format(46h *** error *** illegal code for penetration in,
     &        46h radiation (1st entry must be 0) *** no action,
     &        30h taken, but results a,re bad *)
      endif
      goto 100

 2300 do 2350 k = 1, k2
       call aamadcyphr(s,ni,0,dum,ip2(k))
 2350 continue
      if(idebug.eq.1)then
       if(ip2(1).ne.0) write(iw,2275)
      endif
      goto 100

 2400 do 2450 k = 1, k3
       call aamadcyphr(s,ni,0,dum,ip3(k))
 2450 continue
      if(idebug.eq.1)then
       if(ip3(1).ne.0) write(iw,2275)
      endif
      goto 100

 2500 call aamadcyphr(s,ni,0,dum,ic)
      goto 100

 2600 call aamadcyphr(s,ni,0,dum,iread)
      goto 100

 2700 call aamadcyphr(s,ni,0,dum,iw)
      goto 100

 2800 call aamadcyphr(s,ni,0,dum,ipunch)
      goto 100

 2900 call aamadcyphr(s,ni,1,ffactd,idum)
      goto 100

*-----------------------------------------------------------------------
*     read parameter 20-29 ...
*-----------------------------------------------------------------------

 2999 goto (3000,3100,3200,3300,3400,3500,3600,3700,3800,3900), ju

 3000 call aamadcyphr(s,ni,1,amassm,idum)
      goto 100

 3100 call aamadcyphr(s,ni,1,amasse,idum)
      goto 100

 3200 call aamadcyphr(s,ni,1,z,idum)
      goto 100

 3300 call aamadcyphr(s,ni,1,zsk,idum)
      call aamadcyphr(s,ni,1,zsl,idum)
      call aamadcyphr(s,ni,1,zsm,idum)
      goto 100

 3400 continue
      if(idebug.eq.1)then
       if(ivers.ne.1) write(iw,900)
      endif
      if(ivers.ne.1) goto 100
      do 3450 k = 1, 6
       call aamadcyphr(s,ni,0,dum,ifm(k))
 3450 continue
      goto 100

 3500 continue
      if(idebug.eq.1)then
       if(ivers.ne.1) write(iw,900)
      endif
      if(ivers.ne.1) goto 100
      do 3550 k = 1, 9
       call aamadcyphr(s,ni,0,dum,ifd(k))
 3550 continue
      goto 100

 3600 continue
      if(idebug.eq.1)then
       if(ivers.ne.1) write(iw,900)
      endif
      if(ivers.ne.1) goto 100
      do 3650 k = 1, 10
       call aamadcyphr(s,ni,0,dum,ifq(k))
 3650 continue
      goto 100

 3700 continue
      if(idebug.eq.1)then
       if(ivers.ne.1) write(iw,900)
      endif
      if(ivers.ne.1) goto 100
      do 3750 k = 1, 10
       call aamadcyphr(s,ni,0,dum,ifo(k))
 3750 continue
      goto 100

 3800 do 3850 k = 1, k1
       call aamadcyphr(s,ni,0,dum,m1(k))
 3850 continue
      goto 100

 3900 do 3950 k = 1, k2
       call aamadcyphr(s,ni,0,dum,m2(k))
 3950 continue
      goto 100

*-----------------------------------------------------------------------
*     read parameter 30-39 ...
*-----------------------------------------------------------------------

 3999 goto (4000,4100,4200,4300,4400,4500,4600,4700,4800,4900), ju

 4000 do 4050 k = 1, k3
       call aamadcyphr(s,ni,0,dum,m3(k))
 4050 continue
      goto 100

 4100 do 4150 k = 1, 6
       call aamadcyphr(s,ni,1,pop(k),idum)
       if(idebug.eq.1)then
        if(pop(k).lt.0.0d0) write(iw,4110) k, pop(k)
 4110   format(21h *** warning *** pop(,i1,3h) =,f6.3,14h *** set to 1.,
     &          7h000 ***)
        if(pop(k).gt.1.0d0) write(iw,4120) k, pop(k)
 4120   format(21h *** warning *** pop(,i1,3h) =,f6.3,14h *** set to 0.,
     &          7h000 ***)
       endif
       if(pop(k).lt.0.0d0) pop(k) = 0.0d0
       if(pop(k).gt.1.0d0) pop(k) = 1.0d0
 4150 continue
      goto 100

 4200 continue
      if(idebug.eq.1)then
       if(ivers.ne.1) write(iw,900)
      endif
      if(ivers.ne.1) goto 100
      do 4250 k = 1, 10
       call aamadcyphr(s,ni,0,dum,jm(k))
 4250 continue
      goto 100

 4300 continue
      if(idebug.eq.1)then
       if(ivers.ne.1) write(iw,900)
      endif
      if(ivres.ne.1) goto 100
      do 4350 k=1,14
       call aamadcyphr(s,ni,0,dum,jd(k))
 4350 continue
      goto 100

 4400 continue
      if(idebug.eq.1)then
       if(ivers.ne.1) write(iw,900)
      endif
      if(ivers.ne.1) goto 100
      do 4450 k = 1, 15
       call aamadcyphr(s,ni,0,dum,jq(k))
 4450 continue
      goto 100

 4500 continue
      if(idebug.eq.1)then
       if(ivers.ne.1) write(iw,900)
      endif
      if(ivers.ne.1) goto 100
      do 4550 k = 1, 15
       call aamadcyphr(s,ni,0,dum,jo(k))
 4550 continue
      goto 100

 4600 call aamadcyphr(s,ni,1,ehigh,idum)
      goto 100

 4700 call aamadcyphr(s,ni,1,elow,idum)
      goto 100

 4800 call aamadcyphr(s,ni,1,climit,idum)
      goto 100

 4900 call aamadcyphr(s,ni,1,eres,idum)
      goto 100

*-----------------------------------------------------------------------
*     read parameter 40-49 ...
*-----------------------------------------------------------------------

 4999 goto (5000,5100,5200,5300,100,5500,5600,5700,5800,5900), ju

 5000 call aamadcyphr(s,ni,0,dum,icc)
      goto 100

 5100 do 5150 k = 1, 5
      call aamadcyphr(s,ni,1,cd(k),idum)
 5150 continue
      goto 100

 5200 call aamadcyphr(s,ni,1,ea,idum)
      call aamadcyphr(s,ni,1,eb,idum)
      goto 100

 5300 write(iw,5350)
      stop
 5350 format(53h *** stop card encountered in the input file *** exec,
     &60htion terminated normally *** no subsequent input cards (if a,
     &8hny) read)

 5500 goto 10000

 5600 call aamadcyphr(s,ni,1,a,idum)
      goto 100

 5700 call aamadcyphr(s,ni,1,cfm,idum)
      call aamadcyphr(s,ni,1,tfm,idum)
      goto 100

 5800 call aamadcyphr(s,ni,1,step,idum)
      call aamadcyphr(s,ni,1,rmatch,idum)
      goto 100

 5900 call aamadcyphr(s,ni,1,widthk,idum)
      goto 100

*-----------------------------------------------------------------------
*     read parameter 50-59 ...
*-----------------------------------------------------------------------

 5999 goto (6000,6100,6200,6300,6400,6500,6600,6700,6800,6900), ju

 6000 call aamadcyphr(s,ni,0,dum,nopt)
      goto 100

 6100 call aamadcyphr(s,ni,0,dum,nmax)
      if(nopt.lt.1) call aamadcyphr(s,ni,1,alexp,idum)
      goto 100

 6200 call aamadcyphr(s,ni,0,dum,l)
      call aamadcyphr(s,ni,1,pl(l+1),idum)
      if(idebug.eq.1)then
       if(l.lt.0 .or. l.gt.19 .or. pl(l+1).lt.0.0d0) write(iw,6250)
 6250  format(53h *** warning *** specifications for inputed initial l,
     & 60h-distribution wrong *** l out of range or population negativ,
     & 8he ******)
      endif
      goto 100

 6300 call aamadcyphr(s,ni,0,dum,nstate)
      call aamadcyphr(s,ni,0,dum,kappa)
      call aamadcyphr(s,ni,1,evacp,idum)
      call aamadcyphr(s,ni,1,ebind,idum)
      if(idebug.eq.1)then
       if(nstate.gt.nmax .or. kappa.gt.2*nstate-1 .or. nstate.le.0 .or.
     &    kappa.le.0) write(iw,6325)
 6325  format(53h *** error *** dirac state indecies negative, zero or,
     & 60h out of limits *** check the q.numbers of the states inputed)
      endif
      if(evacp.lt.0.0d0) ebind = ebind + evacp
      idr = idr + 1
      if(idebug.eq.1)then
       if(evacp.lt.0.0d0) write(iw,6350) nstate, kappa
 6350  format(46h *** vacuum polarization less than zero for n=,i2,
     &        12h, and kappa=,i2,4h ***)
      endif
      energya(nstate,kappa) = ebind + evacp
      if(abs(ebind).le.1.0d-20) energya(nstate,kappa) = ebind + evacp
      goto 100

 6400 call aamadcyphr(s,ni,0,dum,iprint)
      goto 100

 6500 do 6550 k=1,20
       call aamadcyphr(s,ni,0,dum,npol(k))
 6550 continue
      goto 100

 6600 call aamadcyphr(s,ni,0,dum,ipol)
      goto 100

 6700 if(ivers.ne.1) write(iw,900)
      if(ivers.ne.1) goto 100
      call aamadcyphr(s,ni,1,zmk,idum)
      call aamadcyphr(s,ni,1,zml,idum)
      call aamadcyphr(s,ni,1,zmm,idum)
      goto 100

 6800 if(ivres.ne.1) write(iw,900)
      if(ivers.ne.1) goto 100
      call aamadcyphr(s,ni,1,zmkm,idum)
      call aamadcyphr(s,ni,1,zmlm,idum)
      call aamadcyphr(s,ni,1,zmmm,idum)
      goto 100

 6900 if(ivers.ne.1) write(iw,900)
      if(ivers.ne.1) goto 100
      do 6950 k = 1, k1
       call aamadcyphr(s,ni,0,dum,iq1(k))
 6950 continue
      goto 100

*-----------------------------------------------------------------------
*     read parameter 60-69 ...
*-----------------------------------------------------------------------

 6999 goto(7000,7100,7200,7300,7400,7500,7600,7700,7800,7900), ju

 7000 continue
      if(idebug.eq.1)then
       if(ivers.ne.1) write(iw,900)
      endif
      if(ivers.ne.1) goto 100
      do 7050 k = 1, k2
       call aamadcyphr(s,ni,0,dum,iq2(k))
 7050 continue
      goto 100

 7100 continue
      if(idebug.eq.1)then
       if(ivers.ne.1) write(iw,900)
      endif
      if(ivers.ne.1) goto 100
      do 7150 k = 1, k3
       call aamadcyphr(s,ni,0,dum,iq3(k))
 7150 continue
      goto 100

 7200 call aamadcyphr(s,ni,0,dum,idb)
      goto 100

 7300 do 7350 k = 1, 4
      call aamadcyphr(s,ni,1,yc(k),idum)
 7350 continue
      goto 100

 7400 continue
      if(idebug.eq.1)then
       if(ivers.ne.1) write(iw,900)
      endif
      if(ivers.ne.1) goto 100
      call aamadcyphr(s,ni,0,dum,iyc)
      goto 100

 7500 continue
      if(idebug.eq.1)then
       if(ivers.ne.1) write(iw,900)
      endif
      if(ivers.ne.1) goto 100
      do 7550 k = 1, 4
      call aamadcyphr(s,ni,0,dum,ij(k))
 7550 continue
      goto 100

 7600 continue
      if(idebug.eq.1)then
       if(ivers.ne.1) write(iw,900)
      endif
      if(ivers.ne.1) goto 100
      do 7650 k = 1, 4
       call aamadcyphr(s,ni,1,yj(k),idum)   ! S.Abe 2018/12/05, dum -> idum
 7650 continue
      goto 100

 7700 continue
      if(idebug.eq.1)then
       if(ivers.ne.1) write(iw,900)
      endif
      if(ivers.ne.1) goto 100
      do 7750 k = 1, 4
       call aamadcyphr(s,ni,0,dum,jj1(k))
 7750 continue
      goto 100

 7800 do 7850 k = 1, 3
       call aamadcyphr(s,ni,0,dum,ipc(k))
 7850 continue
      goto 100

 7900 call aamadcyphr(s,ni,1,cl1,idum)
      call aamadcyphr(s,ni,1,cl2,idum)
      goto 100

*-----------------------------------------------------------------------
*     read parameter 70-79 ...
*-----------------------------------------------------------------------

 7999 goto(8000,8100,8200,8300,8400,8500,8600,8700,8800,8900), ju

 8000 call aamadcyphr(s,ni,1,esp,idum)
      goto 100

 8100 mpu = mpu + 1
      call aamadcyphr(s,ni,0,dum,n1j)
      call aamadcyphr(s,ni,0,dum,l1j)
      call aamadcyphr(s,ni,0,dum,j1j)
      call aamadcyphr(s,ni,0,dum,n2j)
      call aamadcyphr(s,ni,0,dum,l2j)
      call aamadcyphr(s,ni,0,dum,j2j)
      call aamadcyphr(s,ni,0,dum,irs)
      if(idebug.eq.1)then
       if(n1j.gt.20 .or. n1j.le.1 .or. n2j.gt.n1j .or. n2j.le.0 .or.
     &    l1j.ge.n1j .or. l1j.lt.0 .or. l2j.ge.n2j .or. l2j.lt.0 .or.
     &    j1j*j1j.ne.j1j .or. j2j*j2j.ne.j2j .or. iabs(l1j-l2j).gt.3.or.
     &    irs.gt.2 .or. irs.lt.0 .or. (n2j.eq.n1j.and.n2j.ne.2))
     &  write(iw,8150)
 8150   format(53h *** warning *** specifications for transition to be ,
     &  60hpunched are wrong *** no such line or group exists *** no pu,
     &  8hnch ****)
      endif
      icpu(mpu) = n1j + 32*l1j + 1024*j1j + 2048*n2j + 65536*l2j +
     &            2097152*j2j + 4194304*irs
      goto 100

 8200 call aamadcyphr(s,ni,0,dum,ide)
      goto 100

 8300 call aamadcyphr(s,ni,0,dum,ipn)
      goto 100

 8400 do 8450 k = 1, 6
       call aamadcyphr(s,ni,0,dum,jtm(k))
 8450 continue
      goto 100

 8500 do 8550 k = 1, 6
       call aamadcyphr(s,ni,0,dum,jtd(k))
 8550 continue
      goto 100

 8600 do 8650 k = 1, 6
       call aamadcyphr(s,ni,0,dum,jtq(k))
 8650 continue
      go to 100

 8700 do 8750 k = 1, 6
       call aamadcyphr(s,ni,0,dum,jto(k))
 8750 continue
      goto 100

 8800 call aamadcyphr(s,ni,0,dum,n8)
      call aamadcyphr(s,ni,0,dum,ld8)
      call aamadcyphr(s,ni,0,dum,lu8)
      call aamadcyphr(s,ni,1,a8,idum)
      l8 = lu8 - ld8 + 1
      s8 = 0.0d0
      do 8850 k = 1, l8
       call aamadcyphr(s,ni,1,b8,idum)
       if(idebug.eq.1)then

        if(n8.lt.1 .or. ld8.lt.0 .or. lu8.lt.ld8 .or. n8.gt.20 .or.
     &     lu8.ge.n8 .or. a8.lt.0.0d0 .or. b8.lt.0.0d0)
     &  write(iw,8880)
 8880   format(53h *** warning *** spec.s for the initial l-distributio,
     &  60hn given are wrong *** integers out of limits or reals negati,
     &  8hve *****)
       endif
       k8 = n8*(n8-1)/2 + ld8 + k
       pln(k8) = b8
       s8 = s8 + b8
 8850 continue
      do 8875 k = 1, l8
      k8 = n8*(n8-1)/2 + ld8 + k
      if(a8.le.0.0d0) a8 = s8
      pln(k8) = pln(k8)*a8/dmax1(s8,1.d-20)
 8875 continue
      goto 100

 8900 call aamadcyphr(s,ni,0,dum,ip8)
      goto 100

*-----------------------------------------------------------------------
*     read parameter 80-89 ...
*-----------------------------------------------------------------------

 8999 goto( 100, 100, 100, 100, 100, 100, 100, 100, 100, 100), ju

*-----------------------------------------------------------------------
*     this section is devoted to the discovery and recovery of errors.
*     all possible input data are screened for errors (within reason)
*     this portion may be removed if you are confident that you make
*     no mistakes in the input specifications
*-----------------------------------------------------------------------

10000 continue

      if(idebug.ne.1) goto 9999

      if(z.le.0.0d0) write(iw,10100)
      if(z.gt.99.0d0) write(iw,10200)
      if(z.gt.137.04d0) write(iw,10300)
      if(dmod(z+0.001d0,1.0d0).gt.10.002) write(iw,10400)
10100 format(53h *** error *** atomic number z zero or negative *** e,
     &60hxecution terminated *** ....................................)
10200 format(53h *** warning *** atomic number z too big *** last two,
     &60hdigits will be printed in block letters in the upper left of,
     &8htable **)
10300 format(53h *** warning *** atomic number z too big *** point li,
     &60hke dirac formulae have problems *** no attempt to rectify pr,
     &8hoblem **)
10400 format(53h *** warning *** atomic number z not close to an inte,
     &60hger value *** integer part printed in table,but actual value,
     &8h used **)
      if(z.le.0.0d0) stop

      do 10500 k = 1, 3
       if(be(k).lt.0.0d0) write(iw,10600)
10500 continue
      if(be(1).lt.4.0*be(2) .or. be(2).lt.2.0*be(3)) write(iw,10700)
      if(be(1).gt.15.0d0*z*z .or. be(1).lt.5.0d0*z*z) write(iw,10800)
      if(be(3).gt.1000.0d0) write(iw,10900)
10600 format(53h *** error *** binding energy(ies) negative or zero *,
     &60h** program will halt later *** no attempt to correct problem,
     &8h *******)
10700 format(53h *** warning *** binding energies not in any reasonab,
     &60hle proportion *** possibly enetred out of sequence (must be ,
     &8hk,l,m **)
10800 format(53h *** warning *** binding energies not in any reasonab,
     &60hle range *** possibly in the wrong units (must be in ev) or ,
     &8hz ******)
10900 format(53h *** warning *** binding energy of m shell unreasonab,
     &60hly high *** check your source or disregard if intentionally ,
     &8hset ****)

      if(k0.gt.3 .or. k1.gt.7 .or. k2.gt.7 .or. k3.gt.7)
     & write(iw,11000)   
11000 format(53h *** error *** too many cases specified for the multi,
     &60hpolarities *** max are m/3 d,q,o/7 *** results will be incor,
     &8hrect ***)

      if(k0.lt.0 .or. k1.lt.0 .or. k2.lt.0 .or. k3.lt.0)
     & write(iw,11100)
11100 format(53h *** error *** negative number of cases specified for,
     &60h the multipolarities *** must be positive or zero to skip an,
     &8hyone ***)

      if((k0.eq.2.and.nn0(1).eq.nn0(2)) .or.
     &   (k0.eq.3.and.(nn0(1).eq.nn0(2) .or.
     &                 nn0(1).eq.nn0(3) .or.
     &                 nn0(2).eq.nn0(3)))) write(iw,11200)
11200 format(53h *** error *** duplicate shell specification in monop,
     &60hole rate calculation *** will not abort but results will be ,
     &8hwrong **)

      if((k0.eq.1.and.nn0(1).eq.0) .or. 
     &   (k0.eq.1.and.nn0(1)*nn0(2).eq.0) .or.
     &   (k0.eq.3.and.nn0(1)*nn0(2)*nn0(3).eq.0)) write(iw,11300)
11300 format(53h *** error *** radiation specified for monopole cases,
     &60h *** no such rate exists, so the program might do strange th,
     &8hings ***)                                                     

      if(k1.le.1) goto 12000
      do 11500 k = 2, k1
       if(nn1(k).le.0) write(iw,11600)
       if(m1(k).lt.0) write(iw,11700)
       if(nn1(k).gt.3) write(iw,11900)
       if(m1(k).ge.nn1(k) .and. nn1(k).ne.1) write(iw,11900)
       if(k.eq.2) goto 11500
       k4 = k-1
       do 11400 i = 2, k4
        if(nn1(i).eq.nn1(k) .and. m1(k).eq.m1(i)) write(iw,11800)
        if(nn1(i).eq.nn1(k) .and. m1(i)*m1(k).eq.0) write(iw,11800)
11400  continue
11500 continue
11600 format(53h *** error *** negative shell q. n. specified or dupl,
     &60hicate radiation calculation *** no attempt to fix *** check ,
     &8hnnj ****)
11700 format(53h *** error *** negative subshell code specification *,
     &60h** program will take an unpredictable branch or abort *** ch,
     &8heck mj *)
11800 format(53h *** error *** duplicate or overlapping subshell code,
     &60h specified *** will not abort, but rates will be done twice ,
     &8h********)
11900 format(53h *** error *** non existent shell (.gt.3) or subshell,
     &60h combination (e.g. 2d) *** outcome unpredictable *** check n,
     &8hnj,mj **)

12000 if(k2.le.1) goto 12300
      do 12200 k=2,k2
       if(nn2(k).le.0) write(iw,11600)
       if(m2(k).lt.0) write(iw,11700)
       if(nn2(k).gt.3) write(iw,11900)
       if(m2(k).ge.nn2(k) .and. nn2(k).ne.1) write(iw,11900)
       if(k.eq.2) goto 12200
       k4 = k-1
       do 12100 i = 2, k4
        if(nn2(i).eq.nn2(k) .and. m2(k).eq.m2(i)) write(iw,11800)
        if(nn2(i).eq.nn2(k) .and. m2(i)*m2(k).eq.0) write(iw,11800)
12100  continue
12200 continue

12300 if(k3.le.1) goto 12600
      do 12500 k = 2, k3
       if(nn3(k).le.0) write(iw,11600)
       if(m3(k).lt.0) write(iw,11700)
       if(nn3(k).gt.3) write(iw,11900)
       if(m3(k).ge.nn3(k) .and. nn3(k).ne.1) write(iw,11900)
       if(k.eq.2) goto 12500
       k4 = k-1
       do 12400 i = 2, k4
        if(nn3(i).eq.nn3(k).and.m3(k).eq.m3(i))  write(iw,11800)
        if(nn3(i).eq.nn3(k).and.m3(i)*m3(k).eq.0)  write(iw,11800)
12400  continue
12500 continue

12600 if(k1.le.1) goto 12900
      do 12700 k=2,k1
       if(ip1(k)*ip1(k).ne.ip1(k)) write(iw,12800)
12700 continue
12800 format(53h *** warning *** penetration codes for subshells not ,
     &60hzero or one *** possible erroneous result *** check arrays i,
     &8hpj(k) **)

12900 if(k2.le.1) goto 13100
      do 13000 k = 2, k2
       if(ip2(k)*ip2(k).ne.ip2(k)) write(iw,12800)
13000 continue

13100 if(k3.le.1) goto 13300
      do 13200 k = 2, k3
       if(ip3(k)*ip3(k).ne.ip3(k)) write(iw,12800)
13200 continue

13300 do 13400 i = 1, 6
       if(jtm(i).gt.3 .or. jtm(i).lt.1 .or.
     &    jtd(i).lt.1 .or. jtd(i).gt.3 .or.
     &    jtq(i).lt.1 .or. jtq(i).gt.3 .or.
     &    jto(i).gt.3 .or. jto(i).lt.1) write(iw,13500)
13400 continue
13500 format(53h *** warning *** number of terms in penetration calcu,
     &60hlation outside range (1-3) *** disregard if penetration not ,
     &8hused ***)

      if(ic.lt.0 .or. ic.gt.3) write(iw,13600)
13600 format(53h *** error *** accuracy control option code outside p,
     &60hermissible range (0-3) *** unpredictable results can happen ,
     &8h********)

      if(iread.le.0 .or. iw.le.0 .or. ipunch.le.0) write(iw,13700)
13700 format(53h *** warning *** i/o unit number negative or zero ***,
     &60h unlikely to be read or write unit *** if punch...who knows ,
     &8hresult *)

      if(iread.gt.99 .or. iw.gt.99 .or. ipunch.gt.99) write(iw,13800)
13800 format(53h *** warning *** i/o unit number exceeding 99 *** non,
     &60h standard fortran assignment *** disregard if intentionally ,
     &8hset ****)

      if(iprint.lt.0) write(iw,13900)
      if(iprint.gt.63) write(iw,14000)
13900 format(53h *** error *** print selection option code negative *,
     &60h** modulo routine will figure options erroneously *** check ,
     &8hipr ****)
14000 format(53h *** warning *** print selection option code .gt. 63 ,
     &60h*** last 6 bits of number will be used in print (mod(ipr,64),
     &8h) ******)

      if(ffactd.gt.99.99d0 .or. ffactd.lt.0.01d0)  write(iw,14100)
14100 format(53h *** warning *** factorial divider not in any reasona,
     160hble range (0.01-99.99) *** could cause severe arithmetic pro,
     28hblems **)

      if(idb*idb.ne.idb) write(iw,14200)
14200 format(53h *** warning *** debug option selection switch not ze,
     &60hro or one *** could cause errors in the printing of detailed,
     &8h rates *)

      if(ipn*ipn.ne.ipn) write(iw,14300)
14300 format(53h *** warning *** punch selection switch not zero or o,
     &60hne *** could result in unintentional inclusion or omission o,
     &8hf punch*)

      if(ide.lt.0 .or. ide.gt.99999) write(iw,14400)
14400 format(53h *** warning *** punch card identification number neg,
     &60hative ot .gt. 99999 *** will punch 5 stars instead, if .gt. ,
     &8h5 digits)

      if(widthk.lt.0.0d0.or.widthk.gt.999.999d0) write(iw,14500)
14500 format(53h *** warning *** refilling width of k-electron shell ,
     &60hnegative or unreasonably large *** check for proper units (e,
     &8hv) *****)

      if(d2p1s.lt.z*z*amassm .or. d2p1s.gt.10.2d0*z*z*amassm)
     & write(iw,14600)
14600 format(53h *** warning *** energy of the 2p-1s muonic transitio,
     &60hn negative, too low or unreasonably high *** check for units,
     &8h (ev) **)

      if(esp.lt.0.0d0.or.esp.gt.1.0d7) write(iw,14700)
14700 format(53h *** warning *** energy of the 2s-2p muonic transitio,
     &60hn negative or too large *** set to zero if transition to be ,
     &8hskipped*)

      if(nopt.lt.-1 .or. nopt.gt.2) write(iw,14800)
14800 format(53h *** warning *** initial l-distribution option code u,
     &60hnrecognizable *** could cause unexpected complications if of,
     &8h limits*)

      if(nmax.lt.2 .or. nmax.gt.20) write(iw,14900)
14900 format(53h *** error *** starting n quantum number of the casca,
     &60hde not in the range 2-20 *** program will abort if .gt. 20 o,
     &8hr .lt. 1)

      if(nopt.eq.0 .and. abs(alexp).gt.1.000) write(iw,15000)
15000 format(53h *** warning *** modified statistical l-distribution ,
     &60hexponent too high or too low (neg) *** could cause arith. ov,
     &8herflow *)

      if(ip8*ip8.ne.ip8) write(iw,15100)
15100 format(53h *** warning *** l-distribution table selection code ,
     &60h(top n only or full n-l) not zero or one *** could result in,
     &8h errors*)

      if(abs(cl1).gt.10.0d0 .or. abs(cl2).gt.10.0d0) write(iw,15200)
15200 format(53h *** warning *** initial quadratic l-distribution par,
     &60hameters unreasonably high or low (abs .gt. 10.0) *** possibl,
     &8he errors)

      if(ipc(1)*ipc(1).ne.ipc(1) .or.
     &   ipc(2)*ipc(2).ne.ipc(2) .or.
     &   ipc(3)*ipc(3).ne.ipc(3)) write(iw,15300)
15300 format(53h *** warning *** electron refilling control codes not,
     &60h equal to zero or one *** refilling might not be done proper,
     &8hly *****)

      do 15400 i = 1, nmax
       if(npol(i).lt.-1.or.npol(i).gt.nmax-1) write(iw,15500)
15400 continue
      if(ipol*ipol.ne.ipol)  write(iw,15600)
15500 format(53h *** warning *** polarization code n.s for each l out,
     &60h of range *** polarization might be wrong or program will ab,
     &8hort ****)
15600 format(53h *** warning *** polarization calculation selection s,
     &60hwitch not equal to zero or one *** possible undesired result,
     &8hs ******)

      if(yc(1).lt.0.0d0 .or. yc(1).gt.20.0d0 .or.
     &   yc(2).lt.0.0d0 .or. yc(2).gt.20.0d0 .or.
     &   yc(3).lt.0.0d0 .or. yc(3).gt.20.0d0 .or.
     &   yc(4).lt.0.0d0 .or. yc(4).gt.20.0d0) write(iw,15700)
15700 format(53h *** warning *** cutoff y.s for the multipolarities n,
     &60hegative or unreasonably high *** if high check if so desired,
     &8h *******)

      if(abs(amassm-206.7686d0).gt.1.0d-10) write(iw,15800)
15800 format(53h *** warning *** nonstandard mass of particle (not mu,
     &60hon) *** check if this is intentional and the mass is in elec,
     &8ht. m.s *)

      if(abs(amasse-511003.4d0).gt.1.0d-03) write(iw,15900)
15900 format(53h *** warning *** nonstandard mass for the electron **,
     &60h* check if this is intentional and that the mass in in elect,
     &8h. volts*)

      if(abs(amassn-931.48d0).gt.1.0d-10) write(iw,16000)
16000 format(53h *** warning nonstandard mass for the average nucleon,
     &60h bound mass *** check if this is intentional and the mass in,
     &8h mev ***)

      if(abs(a-140.0d0).gt.1.0d-10 .and.
     &   (a.lt.1.3d0*z.or.a.gt.2.7d0*z)) write(iw,16100)
16100 format(53h *** warning *** atomic weight a not changed from def,
     &60hault or unreasonably high or low *** reduced mass calculatio,
     &8hn off **)

      if(cfm.lt.0.0d0 .or. cfm.gt.7.5d0 .or.
     &   tfm.lt.0.5d0 .or. tfm.gt.3.0d0) write(iw,16200)
16200 format(53h *** warning *** fermi distribution parameters unreas,
     &60honably high or low *** not used in any calculation in this p,
     &8hrogram *)

      if(step.lt.0.0d0 .or. rmatch.lt.0.0d0 .or. rmatch.gt.1.0d4)
     & write(iw,16300)
16300 format(53h *** warning *** step in integration negative or matc,
     &60hhing radius negative or unreasonably large *** not used in p,
     &8hrogram *)

      if(ehigh.gt.30.0d0 .or. ehigh.lt.0.001d0) write(iw,16400)
16400 format(53h *** warning *** high cut of energy in x-ray catalogu,
     &60he too high or too low *** check units (mev) *** ok if intent,
     &8hional **)

      if(elow.gt.ehigh .or. elow.lt.0.001d0) write(iw,16500)
16500 format(53h *** warning *** low cut of energy in x-ray catalogue,
     &60h more than the high cut or too low *** check units (mev) and,
     &8h ehi ***)

      if(climit.gt.0.5d0 .or. climit.lt.1.0d-7) write(iw,16600)
16600 format(53h *** warning *** intensity limit cutoff for the x-ray,
     &60h catalogue too high or too low *** too fow or too many lines,
     &8hwritten*)

      if(eres.lt.1.0d-6 .or. eres.gt.0.05d0) write(iw,16700)
16700 format(53h *** warning *** energy resolution in the x-ray catal,
     &60hogue too high or too low *** possibly wrong units (must be i,
     &8hn mev **)

      if(icc*icc.ne.icc) write(iw,16800)
16800 format(53h *** warning *** inputed dividing points switch is no,
     &60ht zero or one *** possible unwanted dividing points in catal,
     &8hogue ***)

      if(cd(1).gt.0.5d0 .or. cd(2).le.cd(3) .or. cd(3).le.cd(4) .or.
     &   cd(4).lt.cd(5) .or. cd(5).lt.1.2d0*climit) write(iw,16900)
16900 format(53h *** warning *** new star dividing points for the x-r,
     &60hay catalogue too high, too low or unreasonably close spaced ,
     &8h********)

      if(eb.le.0.0d0) write(iw,17000)
17000 format(53h *** warning *** calibration parameter b for the conv,
     &60hersion of energy to channel number negative or 0 *** if 0 pr,
     &8hog. halt)

      if(mpu.gt.200) write(iw,17100)
17100 format(53h *** error *** too many lines specified to be punched,
     &60h *** if you want more than 200 lines increase the dim of icp,
     &8hu in l41)

      if(ipc(1)*ipc(2)*ipc(3).eq.0 .and. ip8.eq.1) write(iw,17200)
      if(ipc(1)*ipc(2)*ipc(3).eq.0 .and. ip8.eq.1) stop
17200 format(53h *** error *** full (n,l) l-distribution requires ref,
     &60hilling of all shells (ipc 1 1 1) *** execution terminated **,
     &8h********)

      if(ipol.eq.0 .and. ip8.eq.1) write(iw,17300)
17300 format(53h *** warning *** full n-l distribution specified and ,
     &60hdepolarization calculation *** depolarization may be wrong *,
     &8h********)

      if(ea.gt.2.0d4 .or. ea.lt.-2.0d4) write(iw,17400)
17400 format(53h *** warning *** calibration energy point ea unreason,
     &60hably high or low *** check for units (must be in kev) ******,
     &8h********)

      if(z.lt.30.0d0 .and. (z-20.0d0)*0.1d0.gt.pop(6)) write(iw,17500)
      if(z.lt.18.0d0 .and. (z-12.0d0)/6.0d0.gt.pop(5)) write(iw,17500)
      if(z.lt.12.0d0 .and. (z-10.0d0)*0.5d0.gt.pop(4)) write(iw,17500)
      if(z.lt.10.0d0 .and. (z-4.0d0)/6.0d0.gt.pop(3)) write(iw,17500)
      if(z.lt.4.0d0 .and. (z-2.0d0)*0.5d0.gt.pop(2)) write(iw,17500)
      if((z.le.20.0d0 .and. pop(6).ne.0.0d0) .or.
     &   (z.le.12.0d0 .and. pop(5).ne.0.0d0) .or.
     &   (z.le.10.0d0 .and. pop(4).ne.0.0d0) .or.
     &   (z.le.4.0d0 .and. pop(3).ne.0.0d0) .or.
     &   (z.le.2.0d0 .and. pop(2).ne.0.0d0)) write(iw,17500)
17500 format(53h *** error *** z is too low to support specified popu,
     &60hlation of electronic shells *** possible erroneous results *,
     &8h********)

*-----------------------------------------------------------------------
 9999 continue

      return
      end


************************************************************************
*                                                                      *
      subroutine aamadcyphr(a,nj,itype,r,i)
*                                                                      *
*     decyphers the numbers for routine rread. finds input errors      *
*                                                                      *
*       last modified by S.Abe on 2018/10/30                           *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

c      integer a, bb, aa, bl, se, co, pl, da, po, ee

*-----------------------------------------------------------------------

      common /aamac01/ iread,iw,ipunch,iprint
!$OMP THREADPRIVATE(/aamac01/)
      common /aamadbg/ idebug

      dimension a(70), aa(20)
      character a*1, aa*1

      dimension bb(10)
      character bb*1
      data bb / '0', '1', '2', '3', '4', '5', '6', '7', '8', '9'/

      character bl*1, se*1, co*1, pl*1, da*1, po*1, ee*1
      data bl,se,co,pl,da,po,ee / ' ', '/', ',', '+', '-', '.', 'e'/

*-----------------------------------------------------------------------

      ig = 0

      do 100 j = nj, 70
       ni1 = j+1
       if((a(j).eq.bl .or. a(j).eq.se .or. a(j).eq.co) .and. ig.eq.0)
     &  goto 100
       if((a(j).eq.bl .or .a(j).eq.se .or. a(j).eq.co) .and. ig.ne.0)
     &  goto 200
       ig = ig + 1
       aa(ig) = a(j)
  100 continue

  200 nj = ni1

      if(itype.eq.1) goto 1000

*-----------------------------------------------------------------------
*     get integer parameter
*-----------------------------------------------------------------------

      do 600 j = 1, ig

       ifl = 0
       if(j.ne.1) goto 300
       if(aa(j).eq.pl .or. aa(j).eq.da) ifl = 1

  300  do 400 k = 1, 10
        if(aa(j).eq.bb(k)) ifl = 1
  400  continue

       if(idebug.eq.1)then
        if(ifl.eq.0) write(iw,500) aa(j)
  500   format(/38h *** error *** illegal data in input /,a1,5h/ ***/
     &  60h *** standard fixup taken (taken as 0), execution continuing)
       endif
       if(ifl.eq.0) aa(j) = bb(1)

  600 continue

      is = 1
      if(aa(1).eq.pl .or. aa(1).eq.da) is = 2
      n = 0

      do 700 j = is, ig
       k = ig - j
       read(aa(j),'(i1)') idum
       n = n + idum*10**k
  700 continue

      i = n
      if(aa(1).eq.da) i = -i
      return

*-----------------------------------------------------------------------
*     get real parameter that not containe "e"
*-----------------------------------------------------------------------

 1000 ie = 0

      do 1100 j = 1, ig
       if(aa(j).eq.ee) ie = 1
 1100 continue
      if(ie.eq.1) goto 2000

      ih = ig + 1
      ifl = 0
      do 1200 j = 1, ig
       if(aa(j).eq.po) ih = j
       if(aa(j).eq.po) ifl = ifl + 1
 1200 continue

      if(ifl.eq.1) goto 1500

      if(idebug.eq.1)then
       if(ifl.eq.0) write(iw,1300)
 1300  format(/48h *** warning *** no decimal point in real number,
     & 39h *** assumed to be at the right end ***)

       if(ifl.gt.1) write(iw,1400)
 1400  format(/53h *** error *** too many decimal points in real number
     & ,57h *** last encountered assumed, others changed to zero ***)
      endif

*-----------------------------------------------------------------------

 1500 do 1800 j = 1, ig
       if(j.eq.ih) goto 1800
       ifl = 0
       if(j.ne.1) goto 1600
       if(aa(j).eq.pl .or. aa(j).eq.da) ifl = 1

 1600  do 1700 k = 1, 10
        if(aa(j).eq.bb(k)) ifl = 1
 1700  continue

       if(idebug.eq.1)then
        if(ifl.eq.0) write(iw,500) aa(j)
       endif
       if(ifl.eq.0) aa(j) = bb(1)

 1800 continue

      is = 1
      if(aa(1).eq.pl .or. aa(1).eq.da) is=2
      re = 0.0d0

      do 1900 j=is,ig
       if(j.eq.ih) goto 1900
       k = ih-j-1
       if(k.lt.0)  k = k + 1
       read(aa(j),'(i1)') idum
       re = re + dble(idum)*10.0d0**k
 1900 continue

      r = re
      if(aa(1).eq.da)  r = -r
      return

*-----------------------------------------------------------------------
*     get real parameter that containe "e"
*-----------------------------------------------------------------------

 2000 ih = 0
      ifl = 0

      do 2100 j = 1, ig
       if(aa(j).eq.po) ih = j
       if(aa(j).eq.po) ifl = ifl + 1
       if(aa(j).eq.ee) ie = j
 2100 continue

      ig1 = ie - 1
      if(ifl.eq.0) ih = ig1 + 1
      if(idebug.eq.1)then
       if(ifl.eq.0) write(iw,1300)
       if(ifl.gt.1) write(iw,1400)
      endif
      if(ie.eq.ig) aa(ig+1) = bb(1)
      if(ie.eq.ig) ig = ig + 1

      do 2400 j = 1, ig1
       if(j.eq.ih) goto 2400
       ifl = 0
       if(j.ne.1) goto 2200
       if(aa(j).eq.pl .or. aa(j).eq.da) ifl = 1
 2200  do 2300 k = 1, 10
        if(aa(j).eq.bb(k)) ifl = 1
 2300  continue
       if(idebug.eq.1)then
        if(ifl.eq.0) write(iw,500) aa(j)
       endif
       if(ifl.eq.0) aa(j) = bb(1)
 2400 continue

      is = 1
      if(aa(1).eq.pl .or. aa(1).eq.da) is=2
      re = 0.0d0
      do 2500 j = is, ig1
       if(j.eq.ih) goto 2500
       k = ih - j - 1
       if(k.lt.0) k = k + 1
       read(aa(j),'(i1)') idum
       re = re + dble(idum)*10.0d0**k
 2500 continue

      if(aa(1).eq.da) re = -re
      is = ie + 1
      if(aa(ie+1).eq.pl .or. aa(ie+1).eq.da) is = ie + 2
      n = 0
      do 2600 j=is,ig
       k = ig-j
       read(aa(j),'(i1)') idum
       n = n + dble(idum)*10**k
 2600 continue

      if(aa(ie+1).eq.da) n = -n
      r = re*10.0d0**n

*-----------------------------------------------------------------------
      return
      end                                                             


************************************************************************
*                                                                      *
        subroutine aamain_old(itz,ita)
*                                                                      *
*       control routine of AAMA program                                *
*       original file was made by V.R.Akylas and P.Vogel               *
*       ref.) Muonic atom cascade program,                             *
*             Computer Physics Communications, vol.15, p.291, 1978     *
*                                                                      *
*       modified by S.Abe on 2015/03/05                                *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      common /aama01/ energ, econs, econst, d2p1s, d2p1sm
!$OMP THREADPRIVATE(/aama01/)
      common /aama02/ zsa(3), be(3), bem(3)
!$OMP THREADPRIVATE(/aama02/)
      common /aama03/ k0, k1, k2, k3
!$OMP THREADPRIVATE(/aama03/)
      common /aama04/ nn0(3), nn1(7), nn2(7), nn3(7)
!$OMP THREADPRIVATE(/aama04/)
      common /aama05/ ip1(7), ip2(7), ip3(7), iq1(7), iq2(7), iq3(7)
!$OMP THREADPRIVATE(/aama05/)
      common /aama06/ z, zsk, zsl, zsm, zskz, zslz, zsmz
!$OMP THREADPRIVATE(/aama06/)
      common /aama11/ rr(18), rau, rad, ra(4), rd(4), rsa(4)
!$OMP THREADPRIVATE(/aama11/)
      common /aama12/ pop(6)
!$OMP THREADPRIVATE(/aama12/)
      common /aama13/ e(1000), ai(1000), energya(20,40), m, ia(1000)
!$OMP THREADPRIVATE(/aama13/)
      common /aama14/ dza, dza2, dredm
!$OMP THREADPRIVATE(/aama14/)
      common /aama15/ pl(20), pln(210)
!$OMP THREADPRIVATE(/aama15/)
      common /aama16/ a, cfm
!$OMP THREADPRIVATE(/aama16/)
      common /aama17/ yc(4)
!$OMP THREADPRIVATE(/aama17/)

      common /aama99/ popul(20),ptran(210,210),cctr(210,210),
     &                esp0,energy0(20,40)
!$OMP THREADPRIVATE(/aama99/)
      common /aamai99/ nmax0
!$OMP THREADPRIVATE(/aamai99/)

cABE 2018/10/30, add threadprivate aamac
      common /aamac01/ iread,iw,ipunch,iprint
!$OMP THREADPRIVATE(/aamac01/)
      common /aamac02/ ffact(60),ffactd
!$OMP THREADPRIVATE(/aamac02/)
      common /aamac03/ pi,coeff,picoef,amassm,amasse,alfa,ne(3),jk(3)
!$OMP THREADPRIVATE(/aamac03/)
      common /aamac04/ amassa,amassn
!$OMP THREADPRIVATE(/aamac04/)
      common /aamac05/ m1(7),m2(7),m3(7),idb
!$OMP THREADPRIVATE(/aamac05/)
      common /aamac06/ zmk,zml,zmm,zmkm,zmlm,zmmm,cl1,cl2,ivers
!$OMP THREADPRIVATE(/aamac06/)
      common /aamac07/ jtm(6),jtd(6),jtq(6),jto(6)
!$OMP THREADPRIVATE(/aamac07/)
      common /aamac12/ hbar,widthk,npol(20),ipol,ide,ip8,ipc(3),nmax
!$OMP THREADPRIVATE(/aamac12/)
      common /aamac13/ ehigh,elow,climit,eres,esp,espm,cd(5),ea,eb,
     &                 icc,mpu,icpu(200),ipn
!$OMP THREADPRIVATE(/aamac13/)
      common /aamac14/ tfm,step,rmatch,alexp,nopt,idir,iyc
!$OMP THREADPRIVATE(/aamac14/)
      common /aamac16/ ic,ll(20)
!$OMP THREADPRIVATE(/aamac16/)
      common /aamac17/ elecbe(104,3)
      common /aamadbg/ idebug

*-----------------------------------------------------------------------

      dimension i9(4), k9(4,4), m9(4,7), j9(9,10), j4(9,5,4), jx(4),
     &          id9(5), ixx(3,20)
      dimension plx(20)

c  ***note -- non-standard way of overprinting -- change as needed
c  ***standard fortran requires data for i9 /1h ,1h+,1h+,1h+/
      data i9, l9 / a, b, c, d , e /

      data j9 / 124, 254, 455, 387, 387, 387, 455, 254, 124, 024,
     &          056, 120, 024, 024, 024, 024, 126, 126, 124, 254,
     &          455, 391, 030, 056, 112, 255, 511, 511, 510, 028,
     &          112, 056, 014, 391, 254, 124, 014, 030, 054, 102,
     &          255 ,511, 006, 015, 015, 511, 511, 384, 508, 254,
     &          007, 391, 254, 124, 124, 254, 385 ,508, 510, 387,
     &          455, 254, 124, 511, 511, 007, 014, 028, 056, 112,
     &          224, 448, 124 ,254, 387, 455, 254, 455, 387, 254,
     &          124, 124, 254, 455, 387, 255, 127, 259, 254 ,124/
      data jx, jb / 1a, 1a, 1a, 1a, 1a /

*-----------------------------------------------------------------------
*     initialize common parameter for private
*-----------------------------------------------------------------------
      econs = 5.505355d-3
      d2p1s = 0.d0

      k0 = 3
      k1 = 4
      k2 = 4
      k3 = 4

      nn0(1) = 1
      nn0(2) = 2
      nn0(3) = 3

      nn1(1) = 0
      nn1(2) = 1
      nn1(3) = 2
      nn1(4) = 3
      nn1(5) = 3
      nn1(6) = 3
      nn1(7) = 3

      ip1(1) = 0
      do i = 2, 7
       ip1(i) = 1
      enddo

      do i = 1, 7
       iq1(i) = 0
      enddo

      do i = 1, 7
       nn2(i) = nn1(i)
       nn3(i) = nn1(i)
       ip2(i) = ip1(i)
       ip3(i) = ip1(i)
       iq2(i) = iq1(i)
       iq3(i) = iq1(i)
      enddo

      yc(1) = 1.d0
      yc(2) = 1.d0
      yc(3) = 1.d0
      yc(4) = 1.d0

      do i = 1, 210
       pln(i) = 0.d0
      enddo

*-----------------------------------------------------------------------
*     initialize
*-----------------------------------------------------------------------

      do i = 1, 20
       pl(i) = 0.0d0
      enddo

      m = 1

      do i = 1, 20
       do j = 1, 40
        energya(i,j) = 0.0d0
       enddo
      enddo

      cfm = 0.d0

*-----------------------------------------------------------------------
*     set parameter from PHITS
*-----------------------------------------------------------------------

      z = dble(itz)
      a = dble(ita)

      zsk = z - 1
      zsl = z - 3.0d0
      zsm = z - 11.0d0
      if( zsl .lt. 0.d0 ) zsl = 0.0d0
      if( zsm .lt. 0.d0 ) zsm = 0.0d0

      be(1) = elecbe(itz-1,1)
      be(2) = elecbe(itz-1,2)
      be(3) = elecbe(itz-1,3)

      pop(1) = z / 2.0d0
      pop(2) = (z-2.0d0)  /  2.0d0
      pop(3) = (z-4.0d0)  /  6.0d0
      pop(4) = (z-10.0d0) /  2.0d0
      pop(5) = (z-12.0d0) /  6.0d0
      pop(6) = (z-18.0d0) / 10.0d0

      do i = 1, 6
       if( pop(i) .ge. 1.0d0 ) pop(i) = 1.0d0
c       if( pop(i) .le. 0.0d0 ) pop(i) = 0.0d0
      enddo

      if( itz .eq. 19 .or. itz .eq. 20 ) then
       pop(6) = 0.0d0
      elseif( itz .eq. 21 ) then
       pop(6) = 0.1d0
      elseif( itz .eq. 22 ) then
       pop(6) = 0.2d0
      elseif( itz .eq. 23 ) then
       pop(6) = 0.3d0
      elseif( itz .eq. 24 .or. itz .eq. 25 ) then
       pop(6) = 0.5d0
      elseif( itz .eq. 26 ) then
       pop(6) = 0.6d0
      elseif( itz .eq. 27 ) then
       pop(6) = 0.7d0
      elseif( itz .eq. 28 ) then
       pop(6) = 0.8d0
      elseif( itz .eq. 29 ) then
       pop(6) = 1.0d0
      endif

*-----------------------------------------------------------------------

      if(idebug.eq.1)then
20000 if(z.le.0.0d0) write(iw,20100)
      if(z.gt.99.0d0) write(iw,20200)
      if(z.gt.137.04d0) write(iw,20300)
      if(dmod(z+0.001d0,1.0d0).gt.10.002) write(iw,20400)
      if(z.le.0.0d0) stop
20100 format(53h *** error *** atomic number z zero or negative *** e,
     &60hxecution terminated *** ....................................)
20200 format(53h *** warning *** atomic number z too big *** last two,
     &60hdigits will be printed in block letters in the upper left of,
     &8htable **)
20300 format(53h *** warning *** atomic number z too big *** point li,
     &60hke dirac formulae have problems *** no attempt to rectify pr,
     &8hoblem **)
20400 format(53h *** warning *** atomic number z not close to an inte,
     &60hger value *** integer part printed in table,but actual value,
     &8h used **)

      if(z.lt.30.0d0 .and. (z-20.0d0)*0.1d0.gt.pop(6)) write(iw,27500)
      if(z.lt.18.0d0 .and. (z-12.0d0)/6.0d0.gt.pop(5)) write(iw,27500)
      if(z.lt.12.0d0 .and. (z-10.0d0)*0.5d0.gt.pop(4)) write(iw,27500)
      if(z.lt.10.0d0 .and. (z-4.0d0)/6.0d0.gt.pop(3)) write(iw,27500)
      if(z.lt.4.0d0 .and. (z-2.0d0)*0.5d0.gt.pop(2)) write(iw,27500)
      if((z.le.20.0d0 .and. pop(6).ne.0.0d0) .or.
     &   (z.le.12.0d0 .and. pop(5).ne.0.0d0) .or.
     &   (z.le.10.0d0 .and. pop(4).ne.0.0d0) .or.
     &   (z.le.4.0d0 .and. pop(3).ne.0.0d0) .or.
     &   (z.le.2.0d0 .and. pop(2).ne.0.0d0)) write(iw,27500)
27500 format(53h *** error *** z is too low to support specified popu,
     &60hlation of electronic shells *** possible erroneous results *,
     &8h********)

      if( nopt .ne. 0 ) then
       write(6,*) "error: l,pl(l) or cl1,cl2 are necessary for nop =/=0"
       stop
      endif
      endif

*-----------------------------------------------------------------------

      call aamaffix

*-----------------------------------------------------------------------
*     replace energya(n,l)
*-----------------------------------------------------------------------

      do n = 1, nmax

       k = 2*n - 1

       do l = 1, k

        l1 = (l-1)/2 + 1
        nj = n
        pt = point(nj,l1)

        if(idebug.eq.1)then
        if ( energya(n,l) .le. 0.0d0 .and. mod(iprint,2) .eq. 0 )
     &   write(iw,1200) n, l, pt
 1200   format(27h no input data for state n=,i2,7h, (lj)=,i2,
     &         29h,  point-like dirac energy = ,f10.6,5h(MeV))
        endif

        if ( energya(n,l) .le. 0.0d0 )
     &   energya(n,l) = pt - energya(n,l)

       enddo

      enddo

*-----------------------------------------------------------------------
*     nmax: starting principal quantum number for the cascade
*-----------------------------------------------------------------------

      if ( nmax .lt. 20 ) then

       mx1 = nmax + 1

       do i = mx1, 20
        pl(i) = 0.0d0
       enddo

      endif

*-----------------------------------------------------------------------
*     modified angular momentum distribution, pl(i)
*     nop = -1: use set parameter, no modification
*         =  0: modified statistical
*         =  2: modified quadratic
*-----------------------------------------------------------------------

      if ( nopt .eq. 0 ) then
       do i = 1, nmax
        pl(i) = dble(2*i-1) * dexp(alexp*dble(i-1))
       enddo
      elseif ( nopt .eq. 2 ) then
       do i = 1, nmax
        pl(i) = 1.0d0 + cl1*dble(i-1) + cl2*dble((i-1)**2)
       enddo
      endif

      ss = 0.0d0

      do i = 1, nmax
       ss = ss + pl(i)
      enddo

      if(idebug.eq.1)then
      if(ss.le.0.0d0) write(iw,2100) (pl(i),i=1,nmax)
 2100 format(47h *** error *** initial l distribution wrong ***/
     &       (1x,5f10.6))
      endif

      if(ss.le.0.0d0) ss = 1.0d0

      do i = 1, nmax
       pl(i) = pl(i)/ss
      enddo

*-----------------------------------------------------------------------
*     set character in m9(i,j)
*-----------------------------------------------------------------------

      if(idebug.eq.1)then
      if(mpu.gt.1 .and. ipn.eq.0) write(ipunch,2300) mpu, ide
 2300 format(14(1h*),19h new fit -- at most,i4,15h lines punched ,
     &       20(1h*),i5)
      endif

      do i = 1, 4
       do j = 1, 7
        m9(i,j) = l9
       enddo
      enddo


      if ( k0 .ne. 0 ) then
       do i = 1, k0
        j = nn0(i)
        m9(1,i) = k9(j+1,1)
       enddo
      endif

      if ( k1 .ne. 0 ) then
       do i = 1, k1
        j = nn1(i)
        k = m1(i)
        m9(2,i) = k9(j+1,k+1)
       enddo
      endif

      if ( k2 .ne. 0 ) then
       do i = 1, k2
        j = nn2(i)
        k = m2(i)
        m9(3,i) = k9(j+1,k+1)
       enddo
      endif

      if ( k3 .ne. 0 ) then
       do i = 1, k3
        j = nn3(i)
        k = m3(i)
        m9(4,i) = k9(j+1,k+1)
       enddo
      endif

*-----------------------------------------------------------------------
*     set integer number of atomic number, z
*     and punch identification number, ide
*-----------------------------------------------------------------------

      iz9 = idint(z+0.5d0)
      iz1 = mod(iz9,10) + 1
      iz2 = iz9/10 + 1
      id1 = mod(ide,10) + 1
      id2 = mod(ide/10,10) + 1
      id3 = mod(ide/100,10) + 1
      id4 = mod(ide/1000,10) + 1
      id5 = ide/10000 + 1

*-----------------------------------------------------------------------
*     set l-distribution 
*-----------------------------------------------------------------------

      if ( ip8 .ne. 0) then

       ss = 0.0d0
       nu = nmax*(nmax+1)/2

       do i = 1, nu
        ss = ss + pln(i)
       enddo

       if(idebug.eq.1)then
       if(ss.le.0.0d0) write(iw,2100)
       endif
       if(ss.le.0.0d0) ss = 1.0d0

       do i = 1, nu
        pln(i) = pln(i)/ss 
       enddo

      endif

*-----------------------------------------------------------------------
*     for printing input echo in output file
*     it does not use in muonab routine, thus skip to 9999
*-----------------------------------------------------------------------

      if(idebug.ne.1) goto 9999

      id9(1) = iz2
      id9(2) = iz1

      write(iw,4000) i9(1)

      do i = 2, 4
       write(iw,4010) i9(i)
      enddo

      do i = 1, 4
       write(iw,4100) i9(i)
      enddo

      do i = 1, 9
       do j = 1, 2
        do k = 1, 4
         j4(i,j,k) = jb
         i99 = id9(j)
         if( mod(j9(1,i99)/2**(9-i),2) .eq. 1 ) j4(i,j,k) = jx(k)
        enddo
       enddo
      enddo

      write(iw,4200) i9(1), ((j4(i,j,1),i=1,9),j=1,2)

      do k = 2, 4
       write(iw,4210) i9(k), ((j4(i,j,k),i=1,9),j=1,2)
      enddo

      do i = 1, 9
       do j = 1, 2
        do k = 1, 4
         j4(i,j,k) = jb
         i99 = id9(j)
         if ( mod(j9(2,i99)/2**(9-i),2) .eq. 1 ) j4(i,j,k) = jx(k)
        enddo
       enddo
      enddo

      write(iw,4300) i9(1), ((j4(i,j,1),i=1,9),j=1,2)

      do k = 2, 4
       write(iw,4310) i9(k), ((j4(i,j,k),i=1,9),j=1,2)
      enddo

      do i = 1, 9
       do j = 1, 2
        do k = 1, 4
         j4(i,j,k) = jb
         i99 = id9(j)
         if ( mod(j9(3,i99)/2**(9-i),2) .eq. 1 ) j4(i,j,k) = jx(k)
        enddo
       enddo
      enddo

      write(iw,4400) i9(1), ((j4(i,j,1),i=1,9),j=1,2)

      do k = 2, 4
       write(iw,4410)i9(k),((j4(i,j,k),i=1,9),j=1,2)
      enddo

      do i = 1, 9
       do j = 1, 2
        do k = 1, 4
         j4(i,j,k) = jb
         i99 = id9(j)
         if ( mod(j9(4,i99)/2**(9-i),2) .eq. 1 ) j4(i,j,k) = jx(k)
        enddo
       enddo
      enddo

      write(iw,4500) i9(1), ((j4(i,j,1),i=1,9),j=1,2), nmax, nopt

      do k = 2, 4
       write(iw,4510) i9(k), ((j4(i,j,k),i=1,9),j=1,2), nmax, nopt
      enddo

      do i = 1, 9
       do j = 1, 2
        do k = 1, 4
         j4(i,j,k) = jb
         i99 = id9(j)
         if ( mod(j9(5,i99)/2**(9-i),2) .eq. 1 ) j4(i,j,k) = jx(k)
        enddo
       enddo
      enddo

      write(iw,4600) i9(1), ((j4(i,j,1),i=1,9),j=1,2), alexp, cl1, cl2

      do k = 2, 4
       write(iw,4610) i9(k), ((j4(i,j,k),i=1,9),j=1,2), alexp, cl1, cl2
      enddo

      do i = 1, 9
       do j = 1, 2
        do k = 1, 4
         j4(i,j,k) = jb
         i99 = id9(j)
         if ( mod(j9(6,i99)/2**(9-i),2) .eq. 1 ) j4(i,j,k) = jx(k)
        enddo
       enddo
      enddo

      if ( ip8 .eq. 0 ) then 
       write(iw,4700) i9(1), ((j4(i,j,1),i=1,9),j=1,2), (pl(i),i=1,10)
      else 
       write(iw,4720) i9(1), ((j4(i,j,1),i=1,9),j=1,2)
      endif

      do k = 2, 4
       if ( ip8 .eq. 0 ) then
        write(iw,4710) i9(k), ((j4(i,j,k),i=1,9),j=1,2), (pl(i),i=1,10)
       else
        write(iw,4730) i9(k), ((j4(i,j,k),i=1,9),j=1,2)
       endif
      enddo

      do i = 1, 9
       do j = 1, 2
        do k = 1, 4
         j4(i,j,k) = jb
         i99 = id9(j)
         if ( mod(j9(7,i99)/2**(9-i),2) .eq. 1 ) j4(i,j,k) = jx(k)
        enddo
       enddo
      enddo

      if ( ip8 .eq. 0 ) then
       write(iw,4800) i9(1), ((j4(i,j,1),i=1,9),j=1,2), (pl(i),i=11,20)
      else
       write(iw,4820) i9(1), ((j4(i,j,1),i=1,9),j=1,2)
      endif

      do k = 2, 4
       if ( ip8 .eq. 0 ) then
        write(iw,4810) i9(k), ((j4(i,j,k),i=1,9),j=1,2), (pl(i),i=11,20)
       else
        write(iw,4830) i9(k), ((j4(i,j,k),i=1,9),j=1,2)
       endif
      enddo

      do i = 1, 9
       do j = 1, 2
        do k = 1, 4
         j4(i,j,k) = jb
         i99 = id9(j)
         if ( mod(j9(8,i99)/2**(9-i),2) .eq. 1 ) j4(i,j,k) = jx(k)
        enddo
       enddo
      enddo

      write(iw,4900) i9(1), ((j4(i,j,1),i=1,9),j=1,2)

      do k = 2, 4
       write(iw,4910) i9(k), ((j4(i,j,k),i=1,9),j=1,2)
      enddo

      do i = 1, 9
       do j = 1, 2
        do k = 1, 4
         j4(i,j,k) = jb
         i99 = id9(j)
         if ( mod(j9(9,i99)/2**(9-i),2) .eq. 1 ) j4(i,j,k) = jx(k)
        enddo
       enddo
      enddo

      write(iw,5000) i9(1), ((j4(i,j,1),i=1,9),j=1,2)

      do k = 2, 4
       write(iw,5010) i9(k), ((j4(i,j,k),i=1,9),j=1,2)
      enddo

      write(iw,5100) i9(1)

      do k = 2, 4
       write(iw,5110) i9(k)
      enddo

      write(iw,5200) i9(1), zsk, zsl, zsm, (pop(i),i=1,3)

      do k = 2, 4
       write(iw,5210) i9(k), zsk, zsl, zsm, (pop(i),i=1,3)
      enddo

      write(iw,5300) i9(1), (pop(i),i=4,6)

      do k = 2, 4
       write(iw,5310) i9(k), (pop(i),i=4,6)
      enddo

      write(iw,5400) i9(1)

      do k = 2, 4
       write(iw,5410) i9(k)
      enddo

      write(iw,5500) i9(1), ipc, widthk

      do k = 2, 4
       write(iw,5510) i9(k), ipc, widthk
      enddo

      write(iw,5600) i9(1), be

      do k = 2, 4
       write(iw,5610) i9(k), be
      enddo

      write(iw,5700) i9(1)

      do k = 2, 4
       write(iw,5710) i9(k)
      enddo

*-----------------------------------------------------------------------
*     tfm
*-----------------------------------------------------------------------

      a9 = a
      if ( amassa .ne. 0.0d0 ) a9 = amassa
      tqm = dabs(tfm-2.3001d0)

      if ( tqm .ge. 1.0d-10 ) then
       write(iw,5800) i9(1), a9, cfm, tfm
      else
       write(iw,5820) i9(1), a9
      endif

      do k = 2, 4
       if( tqm .ge. 1.0d-10 ) then
        write(iw,5810) i9(k), a9, cfm, tfm
       else
        write(iw,5830) i9(k), a9
       endif
      enddo

      if ( step+rmatch .gt. 1.0d-20 ) then
       write(iw,5900) i9(1),step,rmatch
      else
       write(iw,5920) i9(1)
      endif

      do k = 2, 4
       if(step+rmatch.gt.1.0d-20) then
        write(iw,5910) i9(k), step, rmatch
       else
        write(iw,5930) i9(k)
       endif
      enddo

*-----------------------------------------------------------------------

      write(iw,6000) i9(1), amassm, amasse, amassn
      do k = 2, 4
       write(iw,6010) i9(k), amassm, amasse, amassn
      enddo

      write(iw,6100) i9(1), d2p1s, esp

      do k = 2, 4
       write(iw,6110) i9(k), d2p1s, esp
      enddo

      write(iw,6200) i9(1)

      do k = 2, 4
       write(iw,6210) i9(k)
      enddo

      write(iw,6300) i9(1), ehigh, climit

      do k = 2, 4
       write(iw,6310) i9(k), ehigh, climit
      enddo

      write(iw,6400) i9(1), elow, icc

      do k = 2, 4
       write(iw,6410) i9(k), elow, icc
      enddo

*-----------------------------------------------------------------------

      eab = (ea-99.0d0)**2 + (eb-99.0d0)**2

      if ( eab .gt. 1.0d-20 ) then
       write(iw,6500) i9(1), eres, ea, eb
      else
       write(iw,6520) i9(1), eres
      endif

      do k = 2, 4
       if ( eab .gt. 1.0d-20 ) then
        write(iw,6510) i9(k), eres, ea, eb
       else
        write(iw,6530) i9(k), eres
       endif
      enddo

      write(iw,6600) i9(1), cd

      do k = 2, 4
       write(iw,6610) i9(k), cd
      enddo

      write(iw,6700) i9(1), npol, ipol

      do k = 2, 4
       write(iw,6710) i9(k), npol, ipol
      enddo

      do k = 1, 4
       write(iw,6800) i9(k)
      enddo

      do k = 1, 4
       write(iw,6900) i9(k)
      enddo

      do k = 1, 4
       write(iw,7000) i9(k)
      enddo

      write(iw,7100) i9(1)

      do k = 2, 4
       write(iw,7110) i9(k)
      enddo

      write(iw,7200) i9(1)

      do k = 2, 4
       write(iw,7210) i9(k)
      enddo

      write(iw,7300) i9(1), k0, k1, k2, k3, iread, iw, ipunch

      do k = 2, 4
       write(iw,7310) i9(k), k0, k1, k2, k3, iread, iw, ipunch
      enddo

      write(iw,7400) i9(1)

      do k = 2, 4
       write(iw,7410) i9(k)
      enddo

      do k = 1, 4
       write(iw,7500) i9(k)
      enddo

      write(iw,7600) i9(1), iprint

      do k = 2, 4
       write(iw,7610) i9(k), iprint
      enddo

      write(iw,7700) i9(1), (m9(1,j),j=1,3), idb

      do k = 2, 4
       write(iw,7710) i9(k), (m9(1,j),j=1,3), idb
      enddo

      write(iw,7800) i9(1), (m9(2,j),j=1,7), ic

      do k = 2, 4
       write(iw,7810) i9(k), (m9(2,j),j=1,7), ic
      enddo

      write(iw,7900) i9(1), (m9(3,j),j=1,7)

      do k = 2, 4
       write(iw,7910) i9(k), (m9(3,j),j=1,7)
      enddo

      write(iw,8000) i9(1), (m9(4,j),j=1,7)
      do k = 2, 4
       write(iw,8010) i9(k), (m9(4,j),j=1,7)
      enddo

      write(iw,8100) i9(1), ffactd
      do k = 2, 4
       write(iw,8110) i9(k), ffactd
      enddo

*-----------------------------------------------------------------------

      mdir = nmax**2

      write(iw,8200) i9(1), idir, mdir

      do k = 2, 4
       write(iw,8210) i9(k), idir, mdir
      enddo

      write(iw,8300) i9(1), mpu

      do k = 2, 4
       write(iw,8310) i9(k), mpu
      enddo

      write(iw,8400) i9(1), ipn

      do k = 2, 4
       write(iw,8410) i9(k), ipn
      enddo

      write(iw,8500) i9(1)

      do k = 2, 4
       write(iw,8510) i9(k)
      enddo

      do k = 1, 4
       write(iw,8600) i9(k)
      enddo

      write(iw,8700) i9(1)

      do k = 2, 4
       write(iw,8710) i9(k)
      enddo

      write(iw,8800) i9(1)

      do k = 2, 4
       write(iw,8810) i9(k)
      enddo

*-----------------------------------------------------------------------

      id9(1) = id5
      id9(2) = id4
      id9(3) = id3
      id9(4) = id2
      id9(5) = id1

      do i = 1, 9
       do j = 1, 5
        do k = 1, 4
         j4(i,j,k) = jb
         i99 = id9(j)
         if ( mod(j9(1,i99)/2**(9-i),2) .eq. 1 ) j4(i,j,k) = jx(k)
        enddo
       enddo
      enddo

      write(iw,8900) i9(1), ((j4(i,j,1),i=1,9),j=1,5)

      do k = 2, 4
       write(iw,8910) i9(k), ((j4(i,j,k),i=1,9),j=1,5)
      enddo

      do i = 1, 9
       do j = 1, 5
        do k = 1, 4
         j4(i,j,k) = jb
         i99 = id9(j)
         if ( mod(j9(2,i99)/2**(9-i),2) .eq. 1 ) j4(i,j,k) = jx(k)
        enddo
       enddo
      enddo

      write(iw,9000) i9(1), ((j4(i,j,1),i=1,9),j=1,5)

      do k = 2, 4
       write(iw,9010) i9(k), ((j4(i,j,k),i=1,9),j=1,5)
      enddo

      do i = 1, 9
       do j = 1, 5
        do k = 1, 4
         j4(i,j,k) = jb
         i99 = id9(j)
         if ( mod(j9(3,i99)/2**(9-i),2) .eq. 1) j4(i,j,k) = jx(k)
        enddo
       enddo
      enddo

      write(iw,9100) i9(1), jtm, ((j4(i,j,1),i=1,9),j=1,5)

      do k = 2, 4
       write(iw,9110) i9(k), jtm, ((j4(i,j,k),i=1,9),j=1,5)
      enddo

      do i = 1, 9
       do j = 1, 5
        do k = 1, 4
         j4(i,j,k) = jb
         i99 = id9(j)
         if ( mod(j9(4,i99)/2**(9-i),2) .eq. 1 ) j4(i,j,k) = jx(k)
        enddo
       enddo
      enddo

      write(iw,9200) i9(1), jtd, ip1, ((j4(i,j,1),i=1,9),j=1,5)

      do k = 2, 4
       write(iw,9210) i9(k), jtd, ip1, ((j4(i,j,k),i=1,9),j=1,5)
      enddo

      do i = 1, 9
       do j = 1, 5
        do k = 1, 4
         j4(i,j,k) = jb
         i99 = id9(j)
         if ( mod(j9(5,i99)/2**(9-i),2) .eq. 1 ) j4(i,j,k) = jx(k)
        enddo
       enddo
      enddo

      write(iw,9300) i9(1), jtq, ip2, ((j4(i,j,1),i=1,9),j=1,5)

      do k = 2, 4
       write(iw,9310) i9(k), jtq, ip2, ((j4(i,j,k),i=1,9),j=1,5)
      enddo

      do i = 1, 9
       do j = 1, 5
        do k = 1, 4
         j4(i,j,k) = jb
         i99 = id9(j)
         if ( mod(j9(6,i99)/2**(9-i),2) .eq. 1 ) j4(i,j,k) = jx(k)
        enddo
       enddo
      enddo

      write(iw,9400) i9(1), jto, ip3, ((j4(i,j,1),i=1,9),j=1,5)

      do k = 2, 4
       write(iw,9410) i9(k), jto, ip3, ((j4(i,j,k),i=1,9),j=1,5)
      enddo

      do i = 1, 9
       do j = 1, 5
        do k = 1, 4
         j4(i,j,k) = jb
         i99 = id9(j)
         if ( mod(j9(7,i99)/2**(9-i),2) .eq. 1 ) j4(i,j,k) = jx(k)
        enddo
       enddo
      enddo

      write(iw,9500) i9(1), ((j4(i,j,1),i=1,9),j=1,5)

      do k = 2, 4
       write(iw,9510) i9(k), ((j4(i,j,k),i=1,9),j=1,5)
      enddo

      do i = 1, 9
       do j = 1, 5
        do k = 1, 4
         j4(i,j,k) = jb
         i99 = id9(j)
         if ( mod(j9(8,i99)/2**(9-i),2) .eq. 1 ) j4(i,j,k) = jx(k)
        enddo
       enddo
      enddo

      do k = 1, 4
       write(iw,9600) i9(k), ((j4(i,j,k),i=1,9),j=1,5)
      enddo

      do i = 1, 9
       do j = 1, 5
        do k = 1, 4
         j4(i,j,k) = jb
         i99 = id9(j)
         if ( mod(j9(9,i99)/2**(9-i),2) .eq. 1) j4(i,j,k) = jx(k)
        enddo
       enddo
      enddo

      write(iw,9700) i9(1), yc, ((j4(i,j,1),i=1,9),j=1,5)

      do k = 2, 4
       write(iw,9710)i9(k),yc,((j4(i,j,k),i=1,9),j=1,5)
      enddo

      do k = 1, 4
       write(iw,9800) i9(k)
      enddo

      do k = 1, 4
       write(iw,9900)i9(k)
      enddo

      write(iw,9910)

*-----------------------------------------------------------------------

 4000 format(1h1/a1,120(1h*))
 4010 format(a1,120(1h*))
 4100 format(a1,1h*,22x,1hi,95x,1h*)
 4200 format(a1,2h* ,9a1,2x,9a1,34h i table of all input parameters -,
     &60h-- defaults (if appropriate) follow the values, in parenthes,
     &4hes *)                                                         
 4210 format(a1,1h*,1x,9a1,2x,9a1,1x,1hi,95x,1h*)
 4300 format(a1,2h* ,9a1,2x,9a1,2h i,95(1h.),1h*)
 4310 format(a1,1h*,1x,9a1,2x,9a1,1x,1hi,95x,1h*)
 4400 format(a1,2h* ,9a1,2x,9a1,2h i,95x,1h*)
 4410 format(a1,1h*,1x,9a1,2x,9a1,1x,1hi,95x,1h*)
 4500 format(a1,2h* ,9a1,2x,9a1,15h i e12 init. n=,i2,12h(max=20,def=,
     &23h15) e11 l-dist. option=,i2,30h(def=0 /-1=inputed,0=statist.,,
     &14h2=quadratic) *)
 4510 format(a1,1h*,1x,9a1,2x,9a1,1x,1hi,13x,i2,35x,i2,43x,1h*)
 4600 format(a1,2h* 9a1,2x,9a1,24h i e12 stat. dist. exp.=,f7.5,
     &25h(0.0) e14 quad. param./a=,f7.5,8h(0.0),b=,f7.5,
     &20h(0.0) 1+a*l+b*l**2 *)
 4610 format(a1,2h* ,9a1,2x,9a1,2h i,22x,f7.5,25x,f7.5,8x,f7.5,19x,
     &1h*)                                                            
 4700 format(a1,2h* ,9a1,2x,9a1,17h i e13 norm.init/,f7.6,9(1x,f7.6),
     &2h *)
 4710 format(a1,2h* ,9a1,2x,9a1,2h i,15x,10(f7.6,1x),1h*)
 4720 format(a1,2h* ,9a1,2x,9a1,7h i e13 ,24h  l-distribution extends,
     &60h beyond starting n.  see next page for complete distribution,
     &1h.,5x,1h*)
 4730 format(a1,1h*,1x,9a1,2x,9a1,2h i,95x,1h*)
 4800 format(a1,2h* ,9a1,2x,9a1,17h i l-dist.(0-19)/,10(f7.6,1x),1h*)
 4810 format(a1,2h* ,9a1,2x,9a1,2h i,15x,10(f7.6,1x),1h*)
 4820 format(a1,2h* ,9a1,2x,9a1,3h i ,4x,8(9x,1h*),10x,1h*)
 4830 format(a1,2h* ,9a1,2x,9a1,2h i,95x,1h*)
 4900 format(a1,2h* ,9a1,2x,9a1,2h i,95(1h.),1h*)
 4910 format(a1,2h* ,9a1,2x,9a1,2h i,95x,1h*)
 5000 format(a1,2h* ,9a1,2x,9a1,2h i,45x,1hi,49x,1h*)
 5010 format(a1,2h* ,9a1,2x,9a1,2h i,45x,1hi,49x,1h*)
 5100 format(a1,1h*,22x,42hi e07 effective charge for electronic shel,
     &55hls  i e08  popul. of el. subshells (fraction of full) *)
 5110 format(a1,1h*,22x,1hi,45x,1hi,49x,1h*)
 5200 format(a1,27h*  e01  z /  (needed)  i k/,f6.3,4h, l/,f6.3,4h, m/
     &,f6.3,21h  (all needed)  i 1s=,f5.3,11h(1.000) 2s=,f5.3,
     &11h(1.000) 2p=,f5.3,9h(1.000) *)
 5210 format(a1,1h*,22x,1hi,3x,f6.3,4x,f6.3,4x,f6.3,16x,1hi,4x,f5.3,
     &11x,f5.3,11x,f5.3,8x,1h*)
 5300 format(a1,1h*,22(1h-),1h+,17(1h-),1h+,27(1h-),5h+ 3s=,f5.3,
     &11h(1.000) 3p=,f5.3,11h(1.000) 3d=,f5.3,9h(1.000) *)
 5310 format(a1,1h*,22x,1h+,17x,1h+,27x,1h+,4x,f5.3,11x,f5.3,11x,f5.3,
     &8x,1h*)
 5400 format(a1,50h*   e10 depletion of electronic shells   i e09 ele,
     &20hct. 1s width in eV +,49(1h-),1h*)
 5410 format(a1,1h*,40x,1hi,27x,1h+,49x,1h*)
 5500 format(a1,1h*,5x,2hk/,i1,3h(0),6x,2hl/,i1,3h(0),6x,2hm/,i1,3h(0)
     &,5x,1hi,6x,f7.3,9h(000.000),5x,29hi e02 av. el. binding ener. f,
     &22hor atom z-1 (needed) *)
 5510 format(a1,1h*,7x,i1,11x,i1,11x,i1,8x,1hi,6x,f7.3,14x,1hi,49x,1h*
     &)

 5600 format(a1,50h* (0=yes,1=no - if k/0,1s width is used) i (exper.,
     &23h or inputed value) i k/,f9.2,8h(eV)  l/,f8.2,8h(ev)  m/,f8.2,
     &6h(eV) *)
 5610 format(a1,1h*,40x,1hi,27x,1hi,3x,f9.2,8x,f8.2,8x,f8.2,5x,1h*)
 5700 format(a1,1h*,40(1h-),1h+,27(1h-),1h+,49(1h-),1h*)
 5710 format(a1,1h*,40x,1h+,27x,1h+,49x,1h*)
 5800 format(a1,20h* e19 atomic weight=,f6.2,21h(140.00) e23 fermi pa,
     &34hrameters (not used in program)  c=,f7.5,17h(fm),skin thick. ,
     &2ht=,f7.5,6h(fm) *)
 5810 format(a1,1h*,19x,f6.2,55x,f7.5,19x,f7.5,5x,1h*)
 5820 format(a1,20h* e19 atomic weight=,f6.2,21h(140.00) e23 fermi pa,
     &60hrameters (not used in program)    * *  n o t   s p e c i f i,
     &13h e d  * *   *)
 5830 format(a1,1h*,19x,f6.2,93x,1h*)
 5900 format(a1,50h* e24 /dirac/ program parameters(not used)/ step i,
     &14hn integration=,1pe9.3,17h matching radius=,e9.3,
     &21h(fm)(for reference) *)
 5910 format(a1,1h*,63x,1pe9.3,17x,e9.3,20x,1h*)
 5920 format(a1,50h* e24 /dirac/ program parameters(not used)/   *  *,
     &46h    n  o  t     s  p  e  c  i  f  i  e  d    *,8(2x,1h*))
 5930 format(a1,1h*,118x,1h*)
 6000 format(a1,23h* e30 masses/ particle=,f9.4,18h(206.7686)(elec. m,
     &17hasses)  electron=,f8.1,24h(511003.4)(eV)  nucleon=,f6.2,
     &15h(931.48)(MeV) *)
 6010 format(a1,1h*,22x,f9.4,35x,f8.1,24x,f6.2,14x,1h*)
 6100 format(a1,50h* e06,e07 special experim. transition energies/  2,
     &5hp-1s=,-6pf8.6,25h(empir. fit)(MeV)  2s-2p=,f8.6,10h(0.0/no tr,
     &14hansit.)(MeV) *)
 6110 format(a1,1h*,54x,-6pf8.6,25x,f8.6,23x,1h*)
 6200 format(a1,1h*,12(1h-),1h+,31(1h-),1h+,55(1h-),1h+,17(1h-),1h*)
 6210 format(a1,1h*,12x,1h+,31x,1h+,55x,1h+,17x,1h*)
 6300 format(a1,1h*,12x,14h/ e20 hi. cut=,f6.3,17h(20.000)MeV i e21,
     &16h intens. cutoff=,1pe9.3,33h(1.000e-06)(per particle) i e33 s,
     &12htar option *)
 6310 format(a1,1h*,12x,1h/,13x,f6.3,12x,1hi,20x,1pe9.3,26x,1hi,17x,
     &1h*)
 6400 format(a1,27h* catalogue  / e20 low cut=,f6.3,14h( 0.040)MeV i ,
     &60he25 calibration parameters (conversion to channel no) i in c,
     &7hatalog/,i1,5h(0) *)
 6410 format(a1,1h*,12x,1h/,13x,f6.3,12x,1hi,55x,1hi,12x,i1,4x,1h*)
 6500 format(a1,27h* parameters / e22 resol. =,f6.5,13h(.00030)MeV i,
     &5x,2ha=,f7.2,11h (keV) , b=,f7.3,27h  channel no =(e-a)/b) i 0=,
     &15hdef, 1=readin *)
 6510 format(a1,1h*,12x,1h/,13x,f6.5,12x,1hi,7x,f7.2,11x,f7.3,23x,
     &1hi,17x,1h*)
 6520 format(a1,27h* parameters / e22 resol. =,f6.5,13h(.00030)MeV i,
     &5x,57h  *  *   n o t   s p e c i f i e d   *  *  *  *   i 0=def,
     &12h, 1=readin *)
 6530 format(a1,1h*,12x,1h/,13x,f6.5,12x,1hi,55x,1hi,17x,1h*)
 6600 format(a1,1h*,12x,24h/ e34 intensities /  5*=,1pe7.1,8h(.1) 4*=,
     &e7.1,9h(.01) 3*=,e7.1,10h(.001) 2*=,e7.1,11h(.0001) 1*=,e7.1,
     &10h(.00001) *)
 6610 format(a1,1h*,12x,1h/,23x,1pe7.1,8x,e7.1,9x,e7.1,10x,e7.1,11x,
     &e7.1,9x,1h*)
 6700 format(a1,18h* e17,e18 quan.dep,4(1h/,i2,1h,,i2,1h,,i2,1h,,i2,
     11h,,i2),28h(all/-1=start ran) do depol=,i1,13h(0/0=y,1=n) *)
 6710 format(a1,1h*,17x,20(1x,i2),28x,i1,12x,1h*)
 6800 format(a1,1h*,118x,1h*)
 6900 format(a1,120(1h*))
 7000 format(a1,1h*,62x,1hi,55x,1h*)
 7100 format(a1,50h* s h e l l  a n d  s u b s h e l l   c o m b i n ,
     &60ha t i o n s  i      b o o k k e e p i n g    p a r a m e t e,
     &10h r s     *)
 7110 format(a1,1h*,62x,1hi,55x,1h*)
 7200 format(a1,1h*,62x,1hi,55x,1h*)
 7210 format(a1,1h*,62x,1hi,55x,1h*)
 7300 format(a1,16h* e38 cases/  m/,i1,16h(max=3,def=3),d/,i1,3h,q/,
     &i1,3h,o/,i1,48h (d,q,o/max=7,def=4) i e31 logical unit no.s rea,
     &2hd=,i1,10h(5),write=,i1,10h(6),punch=,i1,5h(7) *)
 7310 format(a1,1h*,15x,i1,16x,i1,3x,i1,3x,i1,21x,1hi,28x,i1,10x,i1,
     &10x,i1,4x,1h*)
 7400 format(a1,1h*,62(1h.),1hi,55(1h.),1h*)
 7410 format(a1,1h*,62x,1hi,55x,1h*)
 7500 format(a1,1h*,62x,1hi,55x,1h*)
 7600 format(a1,8h* e39,40,1x,2hc1,6x,2hc2,6x,2hc3,6x,2hc4,6x,2hc5,6x,
     &2hc6,6x,2hc7,4x,20hi e35 print option =,i2,17h(0/print all) sel,
     &18hect codes 0 - 63 *)
 7610 format(a1,1h*,62x,1hi,19x,i2,34x,1h*)
 7700 format(a1,7h* mon/ ,a2,6h(1s)  ,a2,6h(2t)  ,a2,13h(3t)  ------ ,
     &46h ------  ------  ------  i e36 debug option = ,i1,
     &35h(0)(0=n,1=y) deb prints all rates *)
 7710 format(a1,1h*,6x,a2,6x,a2,6x,a2,38x,1hi,20x,i1,34x,1h*)
 7800 format(a1,7h* dip/ ,a2,6h(rd)* ,a2,6h(1s)  ,a2,6h(2t)  ,a2,
     &6h(3t)  ,3(a2,6h(--)  ),21hi e32 accuracy check=,i1,
     &35h(3)(0=no, 1 or 2=partial, 3=full) *)
 7810 format(a1,1h*,6x,7(a2,6x),1hi,20x,i1,34x,1h*)
 7900 format(a1,7h* qua/ ,a2,6h(rd)* ,a2,6h(1s)  ,a2,6h(2t)  ,a2,
     &6h(3t)  ,3(a2,6h(--)  ),1hi,55(1h.),1h*)
 7910 format(a1,1h*,6x,7(a2,6x),1hi,55x,1h*)
 8000 format(a1,7h* oct/ ,a2,6h(rd)* ,a2,6h(1s)  ,a2,6h(2t)  ,a2,
     &6h(3t)  ,3(a2,6h(--)  ),1hi,55x,1h*)
 8010 format(a1,1h*,6x,7(a2,6x),1hi,55x,1h*)
 8100 format(a1,50h* key/ rd=radiat., t=total shell, --=not applic., ,
     &30h*=must be rd i e37 fact. div.=,f5.2,21h(15.00) fact. stored ,
     &14hfac(n)/fd**n *)
 8110 format(a1,1h*,62x,1hi,16x,f5.2,34x,1h*)
 8200 format(a1,1h*,62x,36hi e36 no. of dirac energies inputed=,i3,
     &8h out of ,i3,7h max. *)
 8210 format(a1,1h*,62x,1hi,35x,i3,8x,i3,6x,1h*)
 8300 format(a1,63(1h*),37hi e27 max no. of transitions punched=,i3,
     &17h(max=200,def=0) *)
 8310 format(a1,63(1h*),1hi,36x,i3,16x,1h*)
 8400 format(a1,1h*,62x,39hi e29 punch specified transitions /    ,i1,
     &17h(0 =yes, 1 =no) *)
 8410 format(a1,1h*,62x,1hi,38x,i1,16x,1h*)
 8500 format(a1,1h*,6x,43hs e l e c t i o n    o f    p e n e t r a t,
     &6h i o n,7x,1hi,55(1h.),1h*)
 8510 format(a1,1h*,62x,1hi,55x,1h*)
 8600 format(a1,1h*,62x,1hi,55x,1h*)
 8700 format(a1,50h* e42 max number of terms in i e41 penetration sel,
     &60hection codes i e28 punched card identity no. in col.s 73-78 ,
     &10h/(10000) *)
 8710 format(a1,1h*,28x,1hi,33x,1hi,55x,1h*)
 8800 format(a1,50h* penetration (max=3,def.=1) i     (1=y,0=n *=must,
     &14h be 0,def=1) i,55x,1h*)
 8810 format(a1,1h*,28x,1hi,33x,1hi,55x,1h*)
 8900 format(a1,1h*,28x,17hi cases as in e26,17x,1hi,5(1x,9a1,1x),1h*)
 8910 format(a1,1h*,28x,1hi,33x,1hi,5(1x,9a1,1x),1h*)
 9000 format(a1,50h*    1s  2s  2p  3s  3p  3d  i      c1  c2  c3  c4,
     &14h  c5  c6  c7 i,5(1x,9a1,1x),1h*)
 9010 format(a1,1h*,28x,1hi,33x,1hi,5(1x,9a1,1x),1h*)
 9100 format(a1,4h* m/,1x,6(1x,i1,2x),5hi  m/,5h  -  ,6(2x,1h-,1x),1hi
     &,5(1x,9a1,1x),1h*)
 9110 format(a1,1h*,4x,6(1x,i1,2x),1hi,33x,1hi,5(1x,9a1,1x),1h*)
 9200 format(a1,5h* d/ ,6(1x,i1,2x),7hi  d/  ,i1,1h*,6(3x,i1),2h i,
     &5(1x,9a1,1x),1h*)
 9210 format(a1,1h*,4x,6(1x,i1,2x),1hi,6x,i1,1x,6(3x,i1),2h i,
     &5(1x,9a1,1x),1h*)
 9300 format(a1,5h* q/ ,6(1x,i1,2x),7hi  q/  ,i1,1h*,6(3x,i1),2h i,
     &5(1x,9a1,1x),1h*)
 9310 format(a1,1h*,4x,6(1x,i1,2x),1hi,6x,i1,1x,6(3x,i1),2h i,
     &5(1x,9a1,1x),1h*)
 9400 format(a1,5h* o/ ,6(1x,i1,2x),7hi  o/  ,i1,1h*,6(3x,i1),2h i,
     &5(1x,9a1,1x),1h*)
 9410 format(a1,1h*,4x,6(1x,i1,2x),1hi,6x,i1,1x,6(3x,i1),2h i,
     &5(1x,9a1,1x),1h*)
 9500 format(a1,1h*,28(1h.),1hi,33(1h.),1hi,5(1x,9a1,1x),1h*)
 9510 format(a1,1h*,28x,1hi,33x,1hi,5(1x,9a1,1x),1h*)
 9600 format(a1,1h*,62x,1hi,5(1x,9a1,1x),1h*)
 9700 format(a1,17h* e43 y-cutoff m/,f5.2,3h d/,f5.2,3h q/,f5.2,3h o/,
     &f5.2,18h (all def/ 1.00) i,5(1x,9a1,1x),1h*)
 9710 format(a1,1h*,16x,f5.2,3(3x,f5.2),17x,1hi,5(1x,9a1,1x),1h*)
 9800 format(a1,1h*,62x,1hi,55x,1h*)
 9900 format(a1,120(1h*))
 9910 format(1h1)

*-----------------------------------------------------------------------

 9999 continue

      if(ip8.eq.0) goto 11100

*-----------------------------------------------------------------------

      if(idebug.eq.1)then
      write(iw,10000)
      write(iw,10100)
      write(iw,10200)
      endif

      do i = 1, nmax

       n = nmax + 1 -i
       ss = 0.0d0

       do j = 1, n
        k = n*(n-1)/2 + j
        plx(j) = pln(k)
        ss = ss + plx(j)
        k = idint(plx(j)*1.0d6+0.5d0)
        ixx(1,j) = mod(k/100,10)
        ixx(2,j) = mod(k/10,10)
        ixx(3,j) = mod(k,10)
        if(dmod(1000.0d0*plx(j),1.0d0).gt..499999d0)
     &   plx(j)=plx(j)-.000499999d0
        if(plx(j).lt.0.0d0) plx(j) = 0.0d0
       enddo

       if(idebug.eq.1)then
       if ( n .gt. 2 ) then
        write(iw,10300) n, ss, (plx(j),j=1,n)
        write(iw,10400) ((ixx(ii,j),ii=1,3),j=1,n)
       elseif ( n .eq. 2 ) then
        write(iw,10300) n, ss, (plx(j),j=1,n)
        write(iw,10700) ((ixx(ii,j),ii=1,3),j=1,2)
       elseif ( n .eq. 1 ) then
        write(iw,10900) n, ss, plx(1)
        write(iw,11000) (ixx(ii,1),ii=1,3)
       endif

       if ( n .gt. 15 ) then
        write(iw,10450)
       elseif ( n .gt. 10 .and. n .le. 15 ) then
        write(iw,10460)
       elseif ( n .gt. 10 .and. n .le. 15 ) then
        write(iw,10470)
       elseif ( n .gt. 10 .and. n .le. 15 ) then
        write(iw,10480)
       elseif ( n .eq. 3 ) then
        write(iw,10500)
       elseif ( n .eq. 2 ) then
        write(iw,10800)
       endif
       endif

      enddo

*-----------------------------------------------------------------------

11100 call aamacascad
      call aamasort

cabe add
      esp0 = esp
      nmax0 = nmax
      do i = 1, 20
       do j = 1, 40
        energy0(i,j) = energya(i,j)
       enddo
      enddo
cabe end

*-----------------------------------------------------------------------

10000 format(20x,49h*  *   n o r m a l i z e d   i n i t i a l   l - ,
     &30hd i s t r i b u t i o n   *  */)
10100 format(7x,50htotal  i  l=0  l=1  l=2  l=3  l=4 i  l=5  l=6  l=7,
     &60h  l=8  l=9 i l=10 l=11 l=12 l=13 l=14 i l=15 l=16 l=17 l=18 ,
     &4hl=19)
10200 format(1x,13(1h-),1h+,26(1h-),1h+,26(1h-),1h+,26(1h-),1h+,26(1h-))
10300 format(3h n=,i2,1x,f7.5,4(2h i,1x,f4.3,1x,f4.3,1x,f4.3,1x,f4.3,
     &1x,f4.3))
10400 format(13x,4(2h i,2x,3i1,2x,3i1,2x,3i1,2x,3i1,2x,3i1))
10450 format(14x,1hi,26x,1hi,26x,1hi,26x,1hi)
10460 format(14x,1hi,26x,1hi,26x,1hi)
10470 format(14x,1hi,26x,1hi)
10480 format(14x,1hi)
10500 format(14x,1hi,56x,1h+,32(1h-),1h+)
10600 format(3h n=,i2,1x,f7.5,3h i ,f4.3,1x,f4.3,46x,1hi,32x,1hi)
10700 format(14x,1hi,2x,3i1,2x,3i1,46x,27hi entries folded  .abc = .a,
     &7hbcdef i)
10800 format(14x,1hi,56x,34hi to save space    def           i)
10900 format(3h n=,i2,f8.5,3h i ,f4.3,51x,1hi,32x,1hi)
11000 format(14x,1hi,2x,3i1,51x,1h+,32(1h-),1h+/1h1)

*-----------------------------------------------------------------------
      return
      end


************************************************************************
*                                                                      *
      subroutine aamaffix
*                                                                      *
*     initializes variables that cannot be simply assigned in data     *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      common /aama01/ energ, econs, econst, d2p1s, d2p1sm
!$OMP THREADPRIVATE(/aama01/)
      common /aama02/ zsa(3), be(3), bem(3)
!$OMP THREADPRIVATE(/aama02/)
      common /aama06/ z, zsk, zsl, zsm, zskz, zslz, zsmz
!$OMP THREADPRIVATE(/aama06/)
      common /aama07/ amzz(3)
!$OMP THREADPRIVATE(/aama07/)
      common /aama12/ pop(6)
!$OMP THREADPRIVATE(/aama12/)
      common /aama14/ dza, dza2, dredm
!$OMP THREADPRIVATE(/aama14/)
      common /aama16/ a, cfm
!$OMP THREADPRIVATE(/aama16/)
      common /aama17/ yc(4)
!$OMP THREADPRIVATE(/aama17/)

cABE 2018/10/30, add threadprivate aamac
      common /aamac01/ iread,iw,ipunch,iprint
!$OMP THREADPRIVATE(/aamac01/)
      common /aamac02/ ffact(60),ffactd
!$OMP THREADPRIVATE(/aamac02/)
      common /aamac03/ pi,coeff,picoef,amassm,amasse,alfa,ne(3),jk(3)
!$OMP THREADPRIVATE(/aamac03/)
      common /aamac04/ amassa,amassn
!$OMP THREADPRIVATE(/aamac04/)
      common /aamac05/ m1(7),m2(7),m3(7),idb
!$OMP THREADPRIVATE(/aamac05/)
      common /aamac06/ zmk,zml,zmm,zmkm,zmlm,zmmm,cl1,cl2,ivers
!$OMP THREADPRIVATE(/aamac06/)
      common /aamac07/ jtm(6),jtd(6),jtq(6),jto(6)
!$OMP THREADPRIVATE(/aamac07/)
      common /aamac08/ coemon(30),expmon(30),ifm(6),jm(10)
!$OMP THREADPRIVATE(/aamac08/)
      common /aamac09/ coedip(42),expdip(42),coed(4),coedp(9),
     &                 ifd(9),jd(14)
!$OMP THREADPRIVATE(/aamac09/)
      common /aamac10/ coequa(45),expqua(45),coeq(11),ifq(10),jq(15)
!$OMP THREADPRIVATE(/aamac10/)
      common /aamac11/ coeoct(45),expoct(45),coeo(11),ifo(10),jo(15)
!$OMP THREADPRIVATE(/aamac11/)
      common /aamac12/ hbar,widthk,npol(20),ipol,ide,ip8,ipc(3),nmax
!$OMP THREADPRIVATE(/aamac12/)
      common /aamac13/ ehigh,elow,climit,eres,esp,espm,cd(5),ea,eb,
     &                 icc,mpu,icpu(200),ipn
!$OMP THREADPRIVATE(/aamac13/)
      common /aamac14/ tfm,step,rmatch,alexp,nopt,idir,iyc
!$OMP THREADPRIVATE(/aamac14/)
      common /aamadbg/ idebug

*-----------------------------------------------------------------------

      data amasss/206.7686d0/

*-----------------------------------------------------------------------
*     ffact
*-----------------------------------------------------------------------

      ffact(1) = 1.0d0
      do 100 i = 2, 60
       ffact(i) = ffact(i-1)*dble(i-1)/ffactd
  100 continue

*-----------------------------------------------------------------------

      if(idebug.eq.1)then
      if(z.lt.1.0d-20) write(iw,150)
  150 format(52h0*** error *** z not given, execution terminated ***)
      if(z.lt.1.0d-20) stop
      endif

*-----------------------------------------------------------------------
*     zsk, zsl, zsm, zskz, zslz, zsmz, amzz
*-----------------------------------------------------------------------

      if(ivers.ne.1) goto 160

      if(z-zsk.lt.0.0d0 .or. z-zsk.gt.zmkm) zsk = z - zmk
      if(z-zsl.lt.0.0d0 .or. z-zsl.gt.zmlm) zsl = z - zml
      if(z-zsm.lt.0.0d0 .or. z-zsm.gt.zmmm) zsm = z - zmm

  160 continue

      if(idebug.eq.1)then
      if(z-zsk.lt.0.0d0 .or. z-zsk.gt.zmkm) write(iw,170)
      if(z-zsl.lt.0.0d0 .or. z-zsl.gt.zmlm) write(iw,180)
      if(z-zsm.lt.0.0d0 .or. z-zsm.gt.zmmm) write(iw,190)
  170 format(54h *** warning *** effective charge for k shell too high,
     &35h or too low *** no action taken ***)
  180 format(54h *** warning *** effective charge for l shell too high,
     &35h or too low *** no action taken ***)
  190 format(54h *** warning *** effective charge for m shell too high,
     &35h or too low *** no action taken ***)
      endif

      if(zsk.lt.1.0d0) zsk = 1.0d0
      if(zsl.lt.1.0d0) zsl = 1.0d0
      if(zsm.lt.1.0d0) zsm = 1.0d0

      zskz = zsk/z
      zslz = zsl/z
      zsmz = zsm/z

      amzz(1) = zskz/amassm
      amzz(2) = zslz/amassm
      amzz(3) = zsmz/amassm

*-----------------------------------------------------------------------
*     econst
*-----------------------------------------------------------------------

      econst = econs * z**2

*-----------------------------------------------------------------------
*     bem, zsa
*-----------------------------------------------------------------------

      do 200 i = 1, 3
       bem(i) = be(i)/amasse
  200 continue

      if(idebug.eq.1)then
      bm = bem(1)*bem(2)*bem(3)
      if(bm.lt.1.0d-20) write(iw,250)be
  250 format(53h0*** error *** undefined or zero binding energies ***,
     &3f12.3,30h  *** execution terminated ***)
      if(bm.lt.1.0d-20) stop
      endif

      zsa(1) = zsk*alfa
      zsa(2) = zsl*alfa
      zsa(3) = zsm*alfa

*-----------------------------------------------------------------------
*     d2p1sm, espm
*-----------------------------------------------------------------------

      d2p1sm = d2p1s/amasse
      espm = esp/amasse

*-----------------------------------------------------------------------
*     expmon, expdip, expqua, expoct
*-----------------------------------------------------------------------

      if(abs(amasss-amassm).lt.1.0d-20) goto 600

      do 300 i = 1, 30
       expmon(i) = expmon(i)/amassm*amasss
  300 continue
      do 400 i = 1, 42
       expdip(i) = expdip(i)/amassm*amasss
  400 continue
      do 500 i = 1, 45
       expqua(i) = expqua(i)/amassm*amasss
       expoct(i) = expoct(i)/amassm*amasss
  500 continue

      amasss = amassm

*-----------------------------------------------------------------------
*     dza, dza2, dredm
*-----------------------------------------------------------------------

  600 continue

      dza = z*alfa
      dza2 = dza**2
      ame = amasse*1.0d-06
      dredm = a*amassn*amassm*ame/(a*amassn+amassm*ame)
      amasst = amassa*amassn

      if(amassa.gt.1.0d-20)
     & dredm=amasst*amassm*ame/(amasst+amassm*ame)

*-----------------------------------------------------------------------
*     yc(k)
*-----------------------------------------------------------------------

      if(iyc.eq.0) goto 700

      ya = 0.0297d0*z**0.666667d0
      yb = 0.0667d0*dsqrt(z)
      yk = 0.0758d0*dsqrt(z)
      yd = 0.0850d0*dsqrt(z)
      yc(1) = dmin1(yc(1),ya)
      yc(2) = dmin1(yc(2),yb)
      yc(3) = dmin1(yc(3),yk)
      yc(4) = dmin1(yc(4),yd)

*-----------------------------------------------------------------------
*     cfm
*-----------------------------------------------------------------------

  700 continue

      if(cfm.lt.1.0d-20) cfm = 1.100*a**0.333333

      if(amassa.gt.1.0d-20 .and. dabs(cfm-1.1d0*a**0.3333d0).lt.1.0d-20)
     & cfm = 1.1d0*amassa**0.333333d0

*-----------------------------------------------------------------------
*     d2p1s
*-----------------------------------------------------------------------

      if(d2p1s.gt.1.0d-20) goto 800

      r1 = 1.2d0*a**0.3333d0
      x = 2.0d-05*z*r1*amassm/0.529d0
      d2p1sm = econst*(0.75d0 + 3.0d0/x**3*(x*x-4.0d0-x*x*x/3.0d0
     &                                  + dexp(-x)*(x*x+4.0d0+4.0d0*x)))

*-----------------------------------------------------------------------

  800 continue

*-----------------------------------------------------------------------
*     jm, jd, jq, jo
*-----------------------------------------------------------------------

      if(ivers.ne.1) goto 900

      jm(1)=jtm(1)
      jd(1)=jtd(1)
      jq(1)=jtq(1)
      jo(1)=jto(1)
      jm(2)=jtm(2)
      jd(2)=jtd(2)
      jq(2)=jtq(2)
      jo(2)=jto(2)
      jm(3)=jtm(2)
      jd(3)=jtd(2)
      jq(3)=jtq(2)
      jo(3)=jto(2)
      jm(4)=jtm(3)
      jd(4)=jtd(3)
      jq(4)=jtq(3)
      jo(4)=jto(3)
      jd(5)=jtd(3)
      jq(5)=jtq(3)
      jo(5)=jto(3)
      jm(5)=jtm(4)
      jd(6)=jtd(4)
      jq(6)=jtq(4)
      jo(6)=jto(4)
      jm(6)=jtm(4)
      jd(7)=jtd(4)
      jq(7)=jtq(4)
      jo(7)=jto(4)
      jm(7)=jtm(4)
      jd(8)=jtd(4)
      jq(8)=jtq(4)
      jo(8)=jto(4)
      jm(8)=jtm(5)
      jd(9)=jtd(5)
      jq(9)=jtq(5)
      jo(9)=jtq(5)
      jm(9)=jtm(5)
      jd(10)=jtd(5)
      jq(10)=jtq(5)
      jo(10)=jto(5)
      jd(11)=jtd(5)
      jq(11)=jtq(5)
      jo(11)=jto(5)
      jd(12)=jtd(5)
      jq(12)=jtq(5)
      jo(12)=jto(5)
      jm(10)=jtm(6)
      jd(13)=jtd(6)
      jq(13)=jtq(6)
      jo(13)=jto(6)
      jd(14)=jtd(6)
      jq(14)=jtq(6)
      jo(14)=jto(6)
      jq(15)=jtq(6)
      jo(15)=jto(6)

  900 continue

*-----------------------------------------------------------------------
      return
      end


************************************************************************
*                                                                      *
      subroutine aamacascad
*                                                                      *
*     main cascade routine -- does all bookkeeping...                  *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      common /aama00/ ijk
!$OMP THREADPRIVATE(/aama00/)
      common /aama01/ energ, econs, econst, d2p1s, d2p1sm
!$OMP THREADPRIVATE(/aama01/)
      common /aama06/ z, zsk, zsl, zsm, zskz, zslz, zsmz
!$OMP THREADPRIVATE(/aama06/)
      common /aama11/ rr(18), rau, rad, ra(4), rd(4), rsa(4)
!$OMP THREADPRIVATE(/aama11/)
      common /aama12/ pop(6)
!$OMP THREADPRIVATE(/aama12/)
      common /aama15/ pl(20), pln(210)
!$OMP THREADPRIVATE(/aama15/)
      common /aama16/ a, cfm
!$OMP THREADPRIVATE(/aama16/)

      common /aama99/ popul(20),ptran(210,210),cctr(210,210),
     &                esp0,energy0(20,40)
!$OMP THREADPRIVATE(/aama99/)
      common /aamai99/ nmax0
!$OMP THREADPRIVATE(/aamai99/)

cABE 2018/10/30, add threadprivate aamac
      common /aamac01/ iread,iw,ipunch,iprint
!$OMP THREADPRIVATE(/aamac01/)
      common /aamac03/ pi,coeff,picoef,amassm,amasse,alfa,ne(3),jk(3)
!$OMP THREADPRIVATE(/aamac03/)
      common /aamac12/ hbar,widthk,npol(20),ipol,ide,ip8,ipc(3),nmax
!$OMP THREADPRIVATE(/aamac12/)
      common /aamac13/ ehigh,elow,climit,eres,esp,espm,cd(5),ea,eb,
     &                 icc,mpu,icpu(200),ipn
!$OMP THREADPRIVATE(/aamac13/)
      common /aamadbg/ idebug

*-----------------------------------------------------------------------

cABE 2018/12/05 p(3) -> p(4)
      dimension popt(3), p(4), pc(3,210), pnl(210), polpos(210),
     &          polneg(210), width(210), convc(210), sporb(210),
     &          radnt(20), zt(130), zk(3,130), zr(130), ene1(19),
     &          u(4), za0(130) ,za1(130), za2(130), za(130),
     &          pc0(210), pc1(210), pc2(210), pop1(6), pop2(6)

      data u/3hm+q,3hd+o,3h q ,3h o /

*-----------------------------------------------------------------------
*     initialize
*-----------------------------------------------------------------------

      ierror = 0
      maxerr = 99

      do i = 1, 6
       pop1(i) = pop(i)
       pop2(i) = 1.0d0
      enddo

      popt(1) = 2.0d0*pop(1)
      popt(2) = 2.0d0*pop(2) + 6.0d0*pop(3)
      popt(3) = 2.0d0*pop(4) + 6.0d0*pop(5) + 10.0d0*pop(6)

      pop2(1) = pop1(1)

      if ( popt(2) .gt. 1.0d-20 ) then
       pop2(2) = pop2(2)/popt(2)
       pop2(3) = pop2(3)/popt(2)
      endif

      if ( popt(3) .gt. 1.0d-20 ) then
       pop2(4) = pop2(4)/popt(3)
       pop2(5) = pop2(5)/popt(3)
       pop2(6) = pop2(6)/popt(3)
      endif

*-----------------------------------------------------------------------

      rlyman = 0.0d0
      rk = widthk/hbar

      nu = nmax*(nmax+1)/2
      ms = 3

      do j = 1, nu

       if( ip8 .eq. 0 ) then
        pnl(j) = 0.0d0
       else
        pnl(j) = pln(j)
       endif

       polpos(j) = 0.0d0
       polneg(j) = 0.0d0
       width(j) = 0.0d0
       pc0(j) = 0.0d0
       pc1(j) = 0.0d0
       pc2(j) = 0.0d0
       convc(j) = 0.0d0
       sporb(j) = 0.0d0

       do is = 1, ms
        pc(is,j) = 0.0d0
       enddo

      enddo

*-----------------------------------------------------------------------

      mu = nmax*(nmax-1)/2

      do j = 1, nmax

       jj = mu + j

       if ( ip8 .eq. 0 ) then
        pnl(jj) = pl(j)
       endif

       pc1(jj) = 2.0d0 - 2.0d0*pop(1)
       pc2(jj) = 2.0d0*pop(1) - 1.0d0

       if ( pop(1) .le. 0.5d0 ) then
        pc0(jj)=1.0d0-2.0d0*pop(1)
        pc1(jj)=2.0d0*pop(1)
        pc2(jj)=0.0d0
       endif

       if ( pnl(jj) .gt. 1.0d-20 ) then
        polpos(jj) = (1.0d0 + 2.0d0/dble(2*j-1))/3.0d0
        if ( j .ne. 1 ) then
         polneg(jj) = (1.0d0 - 2.0d0/dble(2*j-1))/3.0d0
        endif
       endif

       do is = 1, ms
        pc(is,jj) = popt(is)
       enddo

      enddo

*-----------------------------------------------------------------------

      ipr1 = mod(iprint/4,2)
      ipr2 = mod(iprint/8,2)
      ipr3 = mod(iprint/16,2)
      ipr4 = mod(iprint/32,2)

*-----------------------------------------------------------------------
*     nmax loop
*-----------------------------------------------------------------------

      do i1 = 1, nmax
       n1 = nmax + 1 - i1
       if ( n1 .eq. 1 ) cycle

       if(idebug.eq.1)then
       if ( ipr1 .eq. 0 .and. n1 .gt. 3) write(iw,300)
  300  format(1h1)
       if ( ipr1 .eq. 1 ) write(iw,400)
  400  format(/1x,120(1h*)/)
       if ( ipr1 .eq. 0 .and. n1 .le. 3 ) write(iw,450)
  450  format(/////)
       endif

       call aamacheck(n1)

       do i = 1, 20
        radnt(i) = 0.0d0
       enddo

*-----------------------------------------------------------------------
*     n1 loop
*-----------------------------------------------------------------------

       do i2 = 1, n1

        l1 = i2 - 1
        rategt = 0.0d0
        rate0 = 2.0d0*rk
        rate1 = rk
        raterd = 0.0d0
        k1 = n1*(n1-1)/2 + i2

        if ( i1 .gt. 1 ) then

         if ( pnl(k1) .lt. 1.0d-20 .and. i1 .gt. 1) then

          pc(2,k1) = popt(2)
          pc(3,k1) = popt(3)

          if ( pop2(1) .ge. 0.5d0 ) then
           pc0(k1) = 0.0d0
           pc1(k1) = 2.0d0 - 2.0d0*pop2(1)
           pc2(k1) = 2.0d0*pop2(1) - 1.0d0
          else
           pc0(k1) = 1.0d0 - 2.0d0*pop2(1)
           pc1(k1) = 2.0d0*pop2(1)
           pc2(k1) = 0.0d0
          endif

          pc(1,k1) = pc1(k1) + 2.0d0*pc2(k1)

         else

          xnorm = 1.0d0/pnl(k1)
          pc0(k1) = pc0(k1)*xnorm
          pc1(k1) = pc1(k1)*xnorm
          pc2(k1) = pc2(k1)*xnorm

          do is = 2, ms
           pc(is,k1) = pc(is,k1)*xnorm
           if( pc(is,k1) .le. 0.0d0 ) pc(is,k1) = 0.0d0
          enddo

          pc(1,k1) = pc1(k1) + 2.0d0*pc2(k1)
          if ( pc(1,k1) .lt. 0.0d0 ) then
           pc(1,k1) = 0.0d0
          elseif ( pc(1,k1) .gt. 2.0d0 ) then
           pc(1,k1) = 2.0d0
          endif

          do is = 1, ms
           if ( ipc(is) .ne. 0 ) pc(is,k1) = popt(is)
          enddo

          if ( ipc(1) .ne. 0 ) then
           if ( pop1(1) .ge. 0.5d0 ) then
            pc2(k1) = 2.0d0*pop1(1) - 1.0d0
            pc1(k1) = 2.0d0 - 2.0d0*pop1(1)
            pc0(k1)=0.0d0
           else
            pc2(k1)=0.0d0
            pc1(k1)=2.0d0*pop1(1)
            pc0(k1)=1.0d0 - 2.0d0*pop1(1)
           endif
          endif

          polpos(k1) = polpos(k1)*xnorm*dble(2*l1+1)/dble(l1+1)
          if ( npol(n1)-l1 .le. 1 .and. npol(n1)-l1 .ge. 0 ) then
           polpos(k1) = (1.0d0 + 2.0d0/dble(2*l1+1))/3.0d0
          endif

          if ( l1 .gt. 0 ) then
           polneg(k1) = polneg(k1)*xnorm*dble(2*l1+1)/dble(l1)
           if ( npol(n1)-l1 .le. 1 .and. npol(n1)-l1 .ge. 0 ) then
            polneg(k1) = (1.0d0 - 2.0d0/dble(2*l1+1))/3.0d0
           endif
          endif

         endif

        endif

*-----------------------------------------------------------------------
      
        if ( n1 .eq. 1 ) cycle

        k = 0

        if(idebug.eq.1)then
        if(ipr2.eq.0) write(iw,800)
  800   format(//)
        endif

        pc012 = pc0(k1) + pc1(k1) + pc2(k1)

        if(idebug.eq.1)then
        if ( abs(pc012-1.0d0) .gt. 1.0d-3*dble(i1) )
     &   write(iw,850)n1, l1, pc0(k1), pc1(k1), pc2(k1)
  850  format(53h0*** warning *** probabilities of the population of t,
     & 60hhe k-shell are significantly off from being normalized prope,
     & 8hrly. ***/20h *** happened at n1=,i2,4h l1=,i2,5h pc0=,1pe12.5,
     & 5h pc1=,e12.5,5h pc2=,e12.5,24h normalization forced.../)
        endif

        if ( dabs(pc012) .lt. 1.0d-20 ) pc012 = 1.0d0

        pc0(k1) = pc0(k1)/pc012
        pc1(k1) = pc1(k1)/pc012
        pc2(k1) = pc2(k1)/pc012

        pop(1) = 0.5d0*pc(1,k1)
        pop(2) = pop2(2)*pc(2,k1)
        pop(3) = pop2(3)*pc(2,k1)
        pop(4) = pop2(4)*pc(3,k1)
        pop(5) = pop2(5)*pc(3,k1)
        pop(6) = pop2(6)*pc(3,k1)

*-----------------------------------------------------------------------
*     l2 and n2 loop
*-----------------------------------------------------------------------

        do i3 = 1, 7

         l2 = l1 - 4 + i3

         if ( l2 .lt. 0 ) cycle

         do i4 = 1, n1

          n2 = n1 - i4 + 1
          if ( ( n2 .le. l2 ) .or.
     &         ( n2 .eq. n1 .and. (n2.ne.2.or.l1.ne.0.or.l2.ne.1) ) .or.
     &         ( n2 .eq. n1 .and. espm .le. 1.0d-20 ) ) cycle

          k = k + 1

          if ( n1 .ne. n2 ) then
           ijk = 0
          else
           ijk = 1
           energ = espm
          endif

          popq = pop(1)
          pop(1) = 1.0d0
          zt(k) = rate(n1,l1,n2,l2) + 1.0d-10

          rategt = rategt + zt(k)

cabe add for save ptran
          kk1 = n1*(n1-1)/2 + l1 + 1
          kk2 = n2*(n2-1)/2 + l2 + 1
          ptran(kk1,kk2) = (zt(k) - 1.0d-10) * 1.d-12
cabe end

          if( pop(1) .lt. 1.0d-20 ) pop(1)=1.0d0

          rsa(1) = rsa(1)/pop(1)
          pop(1) = popq
          zr(k) = rad
          raterd = raterd + rad
          za(k) = rsa(1)
          xr = zt(k) - rsa(1)
          t0 = xr + 0.5d0*rsa(1) + rk
          g0 = xr + 2.0d0*rk + 1.0d-10

          rate1 = rate1 + t0 - rk + 1.0d-10
          rate0 = rate0 + xr + 1.0d-10

          xm = zt(k) - za(k)*(pc0(k1) + 0.5d0*pc1(k1))

          ierr = 0

          if ( dabs(xm) .lt. 1.0d-10 .or. dabs(t0) .lt. 1.0d-10 .or.
     &         dabs(g0) .lt. 1.0d-10. or. dabs(zt(k)) .lt. 1.0d-10)
     &     ierr = 1

          if ( ierr .ne. 0 ) then
           ierror = ierror + 1

           if(idebug.eq.1)then
           write(iw,960)
           write(iw,990) n1, l1, n2, l2, xm, t0, g0, k, zt(k),
     &                   pc0(k1), pc1(k1), pc2(k1), zr(k), xr,
     &                   za(k), za0(k), za1(k), za2(k), ym, xl,
     &                   zk(1,k), zk(2,k), zk(3,k)
           write(iw,970)
           endif

           if ( dabs(xm) .lt. 1.0d-10 ) xm=1.0d0
           if ( dabs(t0) .lt. 1.0d-10 ) t0=1.0d0
           if ( dabs(g0) .lt. 1.0d-10 ) g0=1.0d0
           if ( dabs(zt(k)) .lt. 1.0d-10 ) zt(k)=1.0d0

           if(idebug.eq.1)then
           if ( ierror .gt. maxerr ) then
            write(iw,980)
            stop
           endif
           endif

  960   format(//51h *** internal program error in routine cascad *** b,
     &  60had choice of parameters has resulted in a division by zero *,
     &  2h**/55h *** one of the variables /xm,t0,g0,zt(k),rate0,rate1,r,
     &  58hategt/ is zero *** next lines give more information... ***//)
  970   format(//51h *** standard fixup taken (zero variable put equal ,
     &  35hto one), execution continuing ***//)                         
  980   format(//51h *** too many divide checks *** execution terminate,
     &  17hd --- no dump ***//)
  990   format(8h *** n1=,i2,4h l1=,i2,5h, n2=,i2,4h l2=,i2/5x,3hxm=,1pe
     &  12.5,4h t0=,e12.5,4h g0=,e12.5,5h zt(,i3,2h)=,e12.5/5x,6hpc0(k1,
     &  2h)=,e12.5,9h pc1(k1)=,e12.5,9h pc2(k1)=,e12.5/5x,6hzr(k)=,e12.5
     &  ,4h xr=,e12.5,7h za(k)=,e12.5/5x,4hza0=,e12.5,5h za1=,e12.5,    
     &  5h za2=,e12.5,4h ym=,e12.5,4h xl=,e12.5/5x,8hzk(1,k)=,e12.5,    
     &  9h zk(2,k)=,e12.5,9h zk(3,k)=,e12.5/)                           
          endif

          zk(3,k) = pc(3,k1) - rsa(3)/xm
          if( zk(3,k) .lt. 0.0d0 ) zk(3,k) = 0.0d0

          za2(k) = xr/zt(k)*(pc2(k1) + rk*pc1(k1)/t0 +
     &                       2.0d0*rk*rk*pc0(k1)/t0/g0)

          za1(k) = za(k)*pc2(k1)/zt(k) + 
     &             (xr + rk*za(k)/zt(k))/t0*(pc1(k1) +
     &              2.0d0*rk*pc0(k1)/g0)

          za0(k) = pc1(k1)*za(k)/(2.0d0*t0) +
     &             pc0(k1)/g0*(xr + rk*za(k)/t0)

          ym = za(k)*(pc2(k1)/zt(k) + pc1(k1)*(0.5d0 + rk/zt(k))/t0 +
     &                pc0(k1)*rk/t0/g0*(1.0d0 + 2.0d0*rk/zt(k)))

          xl = rk*(pc1(k1)/t0 + 2.0d0*pc0(k1)/t0/g0*(t0 + rk))

          zk(1,k) = pc(1,k1) - ym + xl

          zk(2,k) = pc(2,k1) - rsa(2)/xm - xl

          if( zk(1,k) .lt. 0.0d0 ) zk(1,k) = 0.0d0
          if( zk(2,k) .lt. 0.0d0 ) zk(2,k) = 0.0d0

         enddo ! n2 loop

        enddo  ! l2 loop

*-----------------------------------------------------------------------

        width(k1) = (rategt*pc2(k1) + (rate1-rk)*pc1(k1) +
     &              (rate0-2.0d0*rk)*pc0(k1))*hbar
        if ( width(k1) .lt. 0.0d0 ) width(k1) = 0.0d0

        ierr = 0

        if ( dabs(rate0) .lt. 1.0d-10 .or.
     &       dabs(rate1) .lt. 1.0d-10 .or.
     &       dabs(rategt).lt. 1.0d-10 ) ierr = 1

        if ( ierr .ne. 0 ) then

         ierror = ierror + 1

         if(idebug.eq.1)then
         write(iw,960)
         write(iw,1150) n1, l1, rate0, rate1, rategt
         write(iw,970)          
 1150  format(8h *** n1=,i2,4h l1=,i2,7h rate0=,1pe12.5,7h rate1=,e12.5
     &        ,8h rategt=,e12.5,4h ***/)
         endif

         if ( dabs(rate0) .lt. 1.0d-10 ) rate0 = 1.0d0
         if ( dabs(rate1) .lt. 1.0d-10 ) rate1 = 1.0d0
         if ( dabs(rategt).lt. 1.0d-10 ) rategt = 1.0d0

         if(idebug.eq.1)then
         if ( ierror .gt. maxerr ) then
          write(iw,980)
          stop
         endif
         endif

        endif

        convc(k1) = raterd/(width(k1)/hbar - raterd + 1.0d-10)
        if ( convc(k1) .lt. 0.0d0 ) convc(k1) = 9.999d+99

        if ( l1 .ne. 0 ) sporb(k1) = 0.15d0*z**4/dble(n1**3*l1*(l1+1))

*-----------------------------------------------------------------------
*     l2 and n2 loop
*-----------------------------------------------------------------------

        k = 0

        do i3 = 1, 7

         l2 = l1 - 4 + i3
         if ( l2 .lt. 0 ) cycle

         do i4 = 1, n1

          n2 = n1 - i4 + 1
          if ( ( n2 .le. l2) .or.
     &         ( n2 .eq. n1 .and. (n2.ne.2.or.l1.ne.0.or.l2.ne.1) ) .or.
     &         ( n2 .eq. n1 .and. espm .le. 1.0d-20 ) ) cycle

          k = k + 1
          k2 = n2*(n2-1)/2 + l2 + 1

          bnorm = pnl(k1)
     &            *(pc2(k1)*zt(k)/rategt +
     &              pc1(k1)*(zt(k)-0.5d0*za(k)+rk*zt(k)/rategt)/rate1 +
     &              pc0(k1)*(zt(k)-za(k)+2.0d0*rk/rate1
     &                                   *(zt(k)-0.5d0*za(k) +
     &                                     rk*zt(k)/rategt))/rate0)

          do is = 1, ms
           pc(is,k2) = pc(is,k2) + zk(is,k)*bnorm
          enddo

          pnl(k2) = pnl(k2) + bnorm
          pc2(k2) = pc2(k2) + za2(k)*bnorm
          pc1(k2) = pc1(k2) + za1(k)*bnorm
          pc0(k2) = pc0(k2) + za0(k)*bnorm

          radint = (pc2(k1)/rategt +
     &              pc1(k1)*(1.0d0+rk/rategt)/rate1 +
     &              pc0(k1)*(1.0d0+2.0d0*rk/rate1*(1.0d0+rk/rategt))/
     &              rate0)*pnl(k1)*zr(k)

cabe add for save ccbtr
          kk1 = n1*(n1-1)/2 + l1 + 1
          kk2 = n2*(n2-1)/2 + l2 + 1
          cctr(kk1,kk2) = radint
cabe end

          radnt(n2) = radnt(n2) + radint

          ene1(n2) = econst*(1.0d0/dble(n2*n2)-1.0d0/dble(n1*n1))

          if ( n2 .eq. n1 ) ene1(n2) = espm
          if ( n2 .eq. 1 .and. d2p1sm .gt. 1.0d-20 )
     &     ene1(n2) = ene1(n2) + d2p1sm - 0.75d0*econst

          ene1(n2) = ene1(n2)*amasse

          ll = iabs(l1-l2)

          call popj(l1,l2,ll,p)

*-----------------------------------------------------------------------

          if ( ipol .eq. 0 ) then

           li = ll + 1

           j1u = 2*l1 + 1
           j2u = 2*l2 + 1

           j1d = j1u - 2
           j2d = j2u - 2


           if ( li .gt. 1 ) then

            if ( l2 .gt. l1 ) then

             polpos(k2) = polpos(k2) +
     &                    bnorm*p(1)*polpos(k1)
     &                    *aamabeta(l1,j1u,l2,j2u,ll)
             polneg(k2) = polneg(k2) +
     &                    bnorm*(p(2)*polpos(k1)
     &                    *aamabeta(l1,j1u,l2,j2d,ll)+
     &                    p(3)*polneg(k1)*aamabeta(l1,j1d,l2,j2d,ll))

            else

             polpos(k2) = polpos(k2) +
     &                    bnorm*(p(1)*polpos(k1)
     &                    *aamabeta(l1,j1u,l2,j2u,ll)+
     &                    p(2)*polneg(k1)*aamabeta(l1,j1d,l2,j2u,ll))
             if ( j1d .ne. -1 .and. j2d .ne. -1 ) then
              polneg(k2) = polneg(k2) +
     &                     bnorm*p(3)*polneg(k1)
     &                     *aamabeta(l1,j1d,l2,j2d,ll)
             endif

            endif

           else

            polpos(k2) = polpos(k2) + bnorm*p(1)*polpos(k1)
            polneg(k2) = polneg(k2) + bnorm*p(2)*polneg(k1)

           endif

          endif

*-----------------------------------------------------------------------

          li = ll
          if ( ll .eq. 0 ) li = 2

          lk = ll + 1

          call aamacode(n1,l1,n2,l2,li,radint)

*-----------------------------------------------------------------------

          if(idebug.eq.1)then
          if ( ipr2 .eq. 0 ) write(iw,1600) n1, l1, n2, l2,
     &                                     u(lk), ene1(n2), radint
 1600   format(4h n1=,i2,5h, l1=,i2,5h, n2=,i2,5h, l2=,i2,6h, mul=,a3,
     &   4h, e=,-6pf11.8,5h(MeV),6h, rad=,1pe12.4,10h(per muon))
          endif

          if ( l1 .eq. 1 .and. l2 .eq. 0 .and. n2 .eq. 1)
     &     rlyman = rlyman + radint

          if ( mpu .le. 0 .or. ipn .ne. 0 ) cycle

          do it = 1, mpu

           if ( icpu(it)/4194304 .ne. 1 ) cycle

           n1j = mod(icpu(it),32)
           n2j = mod(icpu(it)/2048,32)
           l1j = mod(icpu(it)/32,32)
           l2j = mod(icpu(it)/65536,32)

           if ( n1j .ne. n1 .or. n2j .ne. n2 .or.
     &          l1j .ne. l1 .or. l2j .ne. l2) cycle

           if(idebug.eq.1)then
           write(ipunch,1700) n1, l1, n2, l2, ene1(n2), radint, ide
 1700      format(2h 1,2i3,4x,2i3,7x,-6pf9.6,1pe12.4,26x,i5)
           endif

          enddo

*-----------------------------------------------------------------------

         enddo ! n2 loop
        enddo  ! l2 loop

       enddo   ! l1 loop

*-----------------------------------------------------------------------

       if(idebug.eq.1)then
       if ( ipr3 .eq. 0 ) write(iw,3100)
 3100 format(//)
       endif

       no = n1 - 1

       if(idebug.eq.1)then
       if( ipr3 .eq. 0 )
     &  write(iw,3200) (n1,n2,ene1(n2),radnt(n2),n2=1,no)
 3200 format(4h n1=,i2,5h, n2=,i2,4h, e=,-6pf11.8,5h(MeV),6h, rad=,
     &1pe12.4,12h(normalized)/)
       endif

       if ( mpu .eq. 0 .or. ipn .ne. 0) cycle

       do it = 1, mpu

        if ( icpu(it)/4194304 .ne. 2) cycle

        n1j = mod(icpu(it),32)
        n2 = mod(icpu(it)/2048,32)

        if(idebug.eq.1)then
        if ( n1j .eq. n1)
     &   write(ipunch,3300) n1,n2,ene1(n2),radnt(n2),ide
 3300  format(2h 2,i3,7x,i3,10x,-6pf9.6,1pe12.4,26x,i5)
        endif

       enddo

*-----------------------------------------------------------------------

      enddo

*-----------------------------------------------------------------------

      if(idebug.eq.1)then
      write(iw,4025) rlyman
 4025 format(//1h ,60(1h*)/40h lyman series (np-1s) sum of intensities,
     &3h = ,f7.5,10h(per muon)/
     &60h deviation from unity is the sum of transition intensities  /
     &60h ending in the 1s state not through an np-1s radiative tran-/
     &60h sition.  experimentally set to unity for normalization     /
     &1h ,60(1h*)//)
      endif

      if ( ipr4 .ne. 0 ) return

      if(idebug.eq.1)then
      write(iw,4050)
 4050 format(1h1)

      write(iw,4100)
 4100 format(//51h n1 l1  population  polar.up   polar.dn   wid (eV) ,
     &60h  rad/aug    s-o(eV)    k-elect    l-elect    m-elect    ***,
     &6h****  /1h ,120(1h-))
      endif

      pc(1,1) = pc1(1) + 2.0d0*pc2(1)

*-----------------------------------------------------------------------

      do m1 = 1, nmax

       n1 = nmax + 1 - m1

       do ll1 = 1, n1

        l1 = ll1 - 1
        k1 = n1*(n1-1)/2 + ll1
        polneg(k1) = -polneg(k1)*(dble(l1)+0.5d0)/(dble(l1)-0.5d0)

cabe add for get initial population
          if(n1 .eq. nmax) then
           popul(ll1) = pnl(k1)
          endif
cabe end

        if(idebug.eq.1)then
        if ( ipol .eq. 1 ) then

         if ( n1 .eq. 1 ) then
          write(iw,4400) n1, l1, pnl(k1), (pc(is,1),is=1,3)
         elseif( l1 .eq. 0 .and. pnl(k1) .gt. 1.0d-20 ) then
          write(iw,4500) n1, l1, pnl(k1), width(k1), convc(k1),
     &                   (pc(is,k1),is=1,3)
         elseif ( pnl(k1) .le. 1.0d-20 ) then
          write(iw,4130) n1, l1
         else
          write(iw,4600) n1, l1, pnl(k1), width(k1), convc(k1),
     &                   sporb(k1), (pc(is,k1),is=1,ms)
         endif

        else

         if ( n1 .eq. 1 ) then
          write(iw,4110) n1, l1, pnl(k1), polpos(1), (pc(is,1),is=1,3)
         elseif ( l1 .eq. 0 .and. pnl(k1) .gt. 1.0d-20 ) then
          write(iw,4120) n1, l1, pnl(k1), polpos(k1), width(k1), 
     &                   convc(k1), (pc(is,k1),is=1,3)
         elseif( pnl(k1) .le. 1.0d-20 ) then
          write(iw,4130) n1, l1
         else
          write(iw,4200) n1, l1, pnl(k1), polpos(k1), polneg(k1),
     &                   width(k1), convc(k1), sporb(k1),
     &                   (pc(is,k1),is=1,ms)
         endif

        endif

 4110   format(1x,i2,1h,,i2,1x,1p2e11.3,4(5x,3h***,3x),3e11.3)
 4120   format(1x,i2,1h,,i2,1x,1p2e11.3,5x,3h***,
     &         3x,2e11.3,5x,3h***,3x,3e11.3)
 4130   format(1x,i2,1h,,i2,1x,9(5x,3h***,3x),2x,13hnot populated)
 4200   format(1x,i2,1h,,i2,1x,1p10e11.3)
 4400   format(1x,i2,1h,,i2,1x,1pe11.3,5(5x,3h***,3x),3e11.3)
 4500   format(1x,i2,1h,,i2,1x,1pe11.3,2(5x,3h***,3x),2e11.3,
     &         5x,3h***,3x,3e11.3)
 4600   format(1x,i2,1h,,i2,1x,1pe11.3,2(5x,3h***,3x),7e11.3)
       endif

       enddo

       if(idebug.eq.1)then
       write(iw,4950)
 4950  format(1h ,120(1h-))
       endif

      enddo
*-----------------------------------------------------------------------

      do i = 1, 6
       pop(i) = pop1(i)
      enddo

*-----------------------------------------------------------------------
      return
      end


************************************************************************
*                                                                      *
      subroutine aamacheck(n)
*                                                                      *
c  ***  numerical accuracy control -- computes diagonal m.e.
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      real*8 matelu

cABE 2018/12/05, add threadprivate aamac
      common /aamac01/ iread,iw,ipunch,iprint
!$OMP THREADPRIVATE(/aamac01/)
      common /aamac16/ ic,ll(20)
!$OMP THREADPRIVATE(/aamac16/)
      common /aamadbg/ idebug

*-----------------------------------------------------------------------

      dimension aa(3,20)

*-----------------------------------------------------------------------

      if(ic.le.0)  return
      go to (1000,2000,3000),ic

*-----------------------------------------------------------------------

 1000 a = (matelu(n,0,n,0,1)/rav(n,0,1))**2.0d0

      if(idebug.eq.1)then
      write(iw,1100) n, a
 1100 format(/29h *** accuracy control at n = ,i2,8h, dip = ,f13.9,
     &       22h (should be unity) ***)
      endif

      return

*-----------------------------------------------------------------------

 2000 a = (matelu(n,0,n,0,1)/rav(n,0,1))**2.0d0
      b = (matelu(n,0,n,0,2)/rav(n,0,2))**2.0d0
      c = (matelu(n,0,n,0,3)/rav(n,0,3))**2.0d0

      if(idebug.eq.1)then
      write(iw,2100)n,a,b,c   
 2100 format(/29h *** accuracy control at n = ,i2,8h, dip = ,f13.9,
     &       8h, qua = ,f13.9,8h, oct = ,f13.9,
     &       22h (should be unity) ***)
      endif

      return

*-----------------------------------------------------------------------

 3000 do 3100 i = 1, 3
       do 3100 j = 1, n
        i1 = i
        j1 = j - 1
        aa(i,j) = (matelu(n,j1,n,j1,i1)/rav(n,j1,i1))**2.0d0
 3100 continue

      if(idebug.eq.1)then
      write(iw,3200) n, (ll(j),(aa(i,j),i=1,3),j=1,n)
 3200 format(/29h *** accuracy control at n = ,i2,16h (should be unit,
     &       6hy) ***/41h  l    dipoles    quadrupoles   octupoles/
     &       (1x,i2,3f13.9))
      endif

*-----------------------------------------------------------------------

      return
      end                                                    


************************************************************************
*                                                                      *
      subroutine popj(l1,l2,ll,p)
*                                                                      *
c  ***  finds relative population of all possible j-states
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)
      dimension p(4)

*-----------------------------------------------------------------------

      li = ll + 1
      l = max0(l1,l2)
      p(3) = 0.d0
      p(4) = 0.d0

      goto (1000,2000,3000,4000), li

*-----------------------------------------------------------------------

 1000 p(1) = dble(l+1)/dble(2*l+1)
      p(2) = 1.d0 - p(1)

      return

*-----------------------------------------------------------------------

 2000 p(1) = dble((l+1)*(2*l-1))/dble(4*l*l-1)
      p(2) = 1.d0/dble(4*l*l-1)
      p(3) = 1.d0 - (p(1)+p(2))

      return

*-----------------------------------------------------------------------

 3000 if(l1.eq.l2) goto 3500
      p(1) = dble(l+1)/dble(2*l+1)
      p(2) = 2.d0/dble((2*l-3)*(2*l+1))
      p(3) = 1.d0 - (p(1)+p(2))

      return

*-----------------------------------------------------------------------

 3500 d = dble((2*l+1)**2)
      p(1) = dble((l+2)*(2*l-1))/d
      p(2) = dble((l-1)*(2*l+3))/d
      p(3) = 3.d0/d
      p(4) = p(3)

      return

*-----------------------------------------------------------------------

 4000 if(iabs(l1-l2).eq.1) goto 4500
      p(1) = dble(l+1)/dble(2*l+1)
      p(2) = 3.d0/dble((2*l-5)*(2*l+1))
      p(3) = 1.d0 - (p(1)+p(2))

      return

*-----------------------------------------------------------------------

 4500 d = dble(4*l*l-1)
      p(1) = dble((2*l-3)*(l+2))/d
      p(2) = 5.d0/d
      p(3) = 6.d0/d
      p(4) = 1.d0 - (p(1)+p(2)+p(3)) 

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine aamacode(n1,l1,n2,l2,l,cc)
*                                                                      *
c  ***  puts information for line intensities in compact form to save
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      common /aama13/ e(1000), ai(1000), energya(20,40), m, ia(1000)
!$OMP THREADPRIVATE(/aama13/)

cABE 2018/10/30, add threadprivate aamac
      common /aamac01/ iread,iw,ipunch,iprint
!$OMP THREADPRIVATE(/aamac01/)
      common /aamac13/ ehigh,elow,climit,eres,esp,espm,cd(5),ea,eb,
     &                 icc,mpu,icpu(200),ipn
!$OMP THREADPRIVATE(/aamac13/)
      common /aamadbg/ idebug

*-----------------------------------------------------------------------

      if(idebug.eq.1)then
      if((l+l1-l2)/2*2.ne.l+l1-l2) write(iw,100) n1, l1, n2, l2, l
  100 format(/50h *** error *** angular momenta do not match at n1=,i2
     &,5h, l1=,i2,5h, n2=,i2,5h, l2=,i2,16h, multipolarity=,i1,4h ***)
      endif

      if(cc.le.climit) return

      if(idebug.eq.1)then
      if(m.gt.0995) write(iw,200)
  200 format(53h *** attention *** overflowing capability of sorting ,
     &50h*** please restrict criteria for lines to pass ***)
      endif

      e12 = 0.0d0
      e22 = 0.0d0

      e11 = energya(n1,2*l1+1)
      if(l1.ne.0) e12 = energya(n1,2*l1)

      e21 = energya(n2,2*l2+1)
      if(l2.ne.0)  e22 = energya(n2,2*l2)

      if(e21-e11.le.elow .and. n1.ne.n2) return
      if(e21-e11.ge.ehigh .and. n1.ne.n2) return

      if(n1.ne.n2) goto 250

      e11 = 0.0d0
      ez = 0.5d0*(energya(2,2)-energya(2,3))

      e21 = esp*1.0d-6 - ez
      e22 = e21 + 2.0d0*ez

  250 ia0 = 4194304*l + 65536*l2 + 2048*n2 + 32*l1 + n1

      if(idebug.eq.1)then
      if(l.le.0 .or. l.gt.3) write(iw,300) n1, l1, n2, l2, l
  300 format(48h *** error *** unexpected quantum numbers at n1=,i2,
     &5h, l1=,i2,5h, n2=,i2,5h, l2=,i2,4h, l=,i1,4h ***)
      endif

*-----------------------------------------------------------------------

      goto(1000,2000,3000), l

*-----------------------------------------------------------------------
*     abs(l2-l1)  = 1 
*-----------------------------------------------------------------------

 1000 ll = max0(l1,l2)

      ia(m) = ia0 + 2098176
      if(l1.eq.ll) ia(m+1) = ia0 + 2097152
      if(l2.eq.ll) ia(m+1) = ia0 + 1024
      ia(m+2) = ia0

      e(m) = e21-e11
      if(l1.eq.ll) e(m+1) = e21 - e12
      if(l2.eq.ll) e(m+1) = e22 - e11
      e(m+2) = e22 - e12

      an = cc/dble(4*ll*ll-1)
      ai(m) = an*dble((ll+1)*(2*ll-1))
      ai(m+1) = an
      ai(m+2) = an*dble((ll-1)*(2*ll+1))

      m = m + 3
      if(ll.eq.1) m = m - 1
      return

*-----------------------------------------------------------------------
*     abs(l2-l1)  = 0 or 2 
*-----------------------------------------------------------------------

 2000 if(l1.eq.l2) goto 2500

      ll = max0(l1,l2)

      ia(m) = ia0 + 2098176
      if(l1.eq.ll) ia(m+1) = ia0 + 2097152
      if(l2.eq.ll) ia(m+1) = ia0 + 1024
      ia(m+2) = ia0

      e(m) = e21 - e11
      if(l1.eq.ll) e(m+1) = e21 - e12
      if(l2.eq.ll) e(m+1) = e22 - e11
      e(m+2) = e22 - e12

      an = cc/dble((2*ll-3)*(2*ll+1))
      ai(m) = an*dble((ll+1)*(2*ll-3))
      ai(m+1) = 2.0d0*an
      ai(m+2) = an*dble((ll-2)*(2*ll+1))

      m = m + 3
      if(ll.le.2) m = m - 1

      return

*-----------------------------------------------------------------------

 2500 if(l1.eq.0) return

      ia(m) = ia0 + 20978176
      ia(m+1) = ia0 + 1024
      ia(m+2) = ia0 + 2097152
      ia(m+3) = ia0

      e(m) = e21 - e11
      e(m+1) = e22 - e11
      e(m+2) = e21 - e12
      e(m+3) = e22 - e12

      an = cc/dble((2*l1+1)**2)
      ai(m) = an*dble((l1+2)*(2*l1-1))
      ai(m+1) = 3.0d0*an
      ai(m+2) = 3.0d0*an
      ai(m+3) = an*dble((l1-1)*(2*l1+3))

      m = m + 4
      if(l1.eq.1)  m = m - 1

      return

*-----------------------------------------------------------------------
*     abs(l2-l1)  = 3
*-----------------------------------------------------------------------

 3000 if(iabs(l1-l2).eq.1) goto 3500

      ll = max0(l1,l2)

      ia(m) = ia0 + 2098176
      if(l1.eq.ll) ia(m+1) = ia0 + 2097152
      if(l2.eq.ll) ia(m+1) = ia0 + 1024
      ia(m+2) = ia0

      e(m) = e21 - e11
      if(l1.eq.ll) e(m+1) = e21 - e12
      if(l2.eq.ll) e(m+1) = e22 - e11
      e(m+2) = e22 - e12

      an = cc/dble((2*ll-5)*(2*ll+1))
      ai(m) = an*dble((2*ll-5)*(ll+1))
      ai(m+1) = 3.000*an          
      ai(m+2) = an*dble((2*ll+1)*(ll-3))

      m = m + 3
      if(ll.le.3)  m = m - 1

      return

*-----------------------------------------------------------------------

 3500 ll = max0(l1,l2)


      ia(m) = ia0 + 20978176
      ia(m+1) = ia0 + 1024
      ia(m+2) = ia0 + 2097152
      ia(m+3) = ia0

      e(m) = e21 - e11
      e(m+1) = e22 - e11
      e(m+2) = e21 - e12
      e(m+3) = e22 - e12

      an = cc/dble(4*ll*ll-1)
      ai(m) = an*dble((2*ll-3)*(ll+2))
      ai(m+1) = 5.d0*an
      ai(m+2) = 6.d0*an
      ai(m+3) = an*dble((2*ll+3)*(ll-2))

      m = m + 4
      if(ll.le.2)  m = m - 1

*-----------------------------------------------------------------------
      return
      end


************************************************************************
*                                                                      *
      subroutine aamasort
*                                                                      *
c  ***  arranges line intensities in energy for the x-ray table         
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      common /aama06/ z, zsk, zsl, zsm, zskz, zslz, zsmz
!$OMP THREADPRIVATE(/aama06/)
      common /aama13/ e(1000), ai(1000), energya(20,40), m, ia(1000)
!$OMP THREADPRIVATE(/aama13/)
      common /aama15/ pl(20), pln(210)
!$OMP THREADPRIVATE(/aama15/)
      common /aama16/ a, cfm
!$OMP THREADPRIVATE(/aama16/)

cABE 2018/10/30, add threadprivate aamac
      common /aamac01/ iread,iw,ipunch,iprint
!$OMP THREADPRIVATE(/aamac01/)
      common /aamac12/ hbar,widthk,npol(20),ipol,ide,ip8,ipc(3),nmax
!$OMP THREADPRIVATE(/aamac12/)
      common /aamac13/ ehigh,elow,climit,eres,esp,espm,cd(5),ea,eb,
     &                 icc,mpu,icpu(200),ipn
!$OMP THREADPRIVATE(/aamac13/)
      common /aamadbg/ idebug

*-----------------------------------------------------------------------

      dimension st(5),c(5),ll(3)

      data star,blank/1h*,1h /
      data ll/3hdip,3hqua,3hoct/

*-----------------------------------------------------------------------

      ich = 0
      if((ea-99.)**2 + (eb-99.)**2.lt.1.0d-20) ich = 1

      cc = climit**(1.0d0/6.0d0)

      do 100 i = 1, 5
       c(i) = cc**i
       if(icc.ne.0) c(i) = cd(i)
  100 continue

      m = m - 1
      if(m.le.0) return

      m1 = m

      do 170 i = 1, m

  120  if(ai(i).gt.climit) goto 170
       if(i.gt.m1) goto 170

       m1 = m1 - 1

       do 140 j = i, m1
        e(j) = e(j+1)
        ai(j) = ai(j+1)
        ia(j) = ia(j+1)
  140  continue

       goto 120

  170 continue

      m = m1

      if(idebug.eq.1)then
      write(iw,200) z, a, m
  200 format(1h1/1h ,120(1h*)/17h0atomic number = ,f5.1,4x,6hatomic
     &,10h weight = ,f7.3,4x,18hnumber of lines = ,i5//1x,120(1h*)/)
      endif

      if(m.eq.1) goto 500

      mm = m - 1

      do 400 i = 1, mm

       i1 = i + 1

       do 300 j = i1, m
        if(e(j).ge.e(i))  go to 300
        e0 = e(j)
        ai0 = ai(j)
        ia0 = ia(j)
        e(j) = e(i)
        ai(j) = ai(i)
        ia(j) = ia(i)
        e(i) = e0
        ai(i) = ai0
        ia(i) = ia0
  300  continue

  400 continue

*-----------------------------------------------------------------------

  500 n = 0

      ee = 0.0d0
      aa = 0.0d0
      lc = 7


      do 1500 i = 1, m

  600  ii = i + n

       if(ia(ii)/16777216.eq.1) goto 1500

       ia0 = ia(ii)

       n1 = mod(ia0,32)
       l1 = mod(ia0/32,32)
       j1 = (2*l1-1) + 2*mod(ia0/1024,2)

       n2 = mod(ia0/2048,32)
       l2 = mod(ia0/65536,32)
       j2 = (2*l2-1) + 2*mod(ia0/2097152,2)

       l = mod(ia0/4194304,4)
       l = ll(l)

       en = 1000.0d0*e(ii)
       a1 = ai(ii)
       ch = (en-ea)/dmin1(eb,1.0d20)

       do 700 k = 1, 5
        st(k) = blank
        if(a1.gt.c(k)) st(k) = star
  700  continue

       if(idebug.eq.1)then
       if(ich.eq.0) write(iw,800) n1, l1, j1, n2, l2, j2, l, en, a1,
     &                            st, ch
       if(ich.eq.1) write(iw,805) n1, l1, j1, n2, l2, j2, l, en, a1, st
  800  format(4h n1=,i2,4h,l1=,i2,4h,j1=,i2,9h/2    n2=,i2,4h,l2=,i2,
     &        4h,j2=,i2,8h/2    l=,a3,7h    en=,f12.6,13h(keV)    int=,
     &        1pe11.4,4x,5a1,4x,3hch=,0pf9.3)
  805  format(4h n1=,i2,4h,l1=,i2,4h,j1=,i2,9h/2    n2=,i2,4h,l2=,i2,  
     &        4h,j2=,i2,8h/2    l=,a3,7h    en=,f12.6,13h(keV)    int=,
     &        1pe11.4,4x,5a1)
       endif

       lc = lc + 1

       if(mpu.le.0) goto 830
       if(ipn.ne.0) goto 830 

       do 820 it = 1, mpu
        if(icpu(it)/4194304.ne.0) goto 820
        n1j = mod(icpu(it),32)
        n2j = mod(icpu(it)/2048,32)
        l1j = mod(icpu(it)/32,32)
        l2j = mod(icpu(it)/65536,32)
        j1j = mod(icpu(it)/1024,2)
        j2j = mod(icpu(it)/2097152,2)
        j1k = 2*l1j - 1 + 2*j1j
        j2k = 2*l2j - 1 + 2*j2j

        if(n1j.ne.n1 .or. n2j.ne.n2 .or.
     &     l1j.ne.l1 .or. l2j.ne.l2 .or.
     &     j1k.ne.j1 .or. j2k.ne.j2) goto 820

        if(idebug.eq.1)then
        write(ipunch,810) n1, l1, j1j, n2, l2, j2j, en, a1, ide
  810   format(2h 0,2i3,i2,2x,2i3,i2,5x,-3pf9.6,1pe12.4,26x,i5)
        endif

  820  continue

  830  if(idebug.eq.1)then
       if(lc.ge.60) write(iw,850)
  850  format(1h1)
       endif

       if(lc.ge.60) lc = 0
       ia(ii) = ia(ii) + 16777216
       ee = ee + e(ii)*a1
       aa = aa + a1

  900  n = n + 1
       ipn1 = i + n

       if(e(ipn1)-e(i).ge.eres) goto 1000
       if(ia(ipn1)/16777216.eq.1) goto 900

       n11 = mod(ia(ipn1),32)
       n22 = mod(ia(ipn1)/2048,32)

       if(n1.eq.n11 .and. n2.eq.n22) goto 600
 1000  if(dabs(a1-aa).lt.1.0d-20) goto 1300

       en = 1000.0d0*ee/aa
       ch = (en-ea)/dmin1(eb,1.0d20)

       do 1100 k = 1, 5
        st(k) = blank
        if(aa.gt.c(k)) st(k) = star
 1100  continue

       if(idebug.eq.1)then
       if(ich.eq.0) write(iw,1200) n1, n2, en, aa, st, ch
       if(ich.eq.1) write(iw,1250) n1, n2, en, aa, st
 1200  format(4h n1=,i2,1x,13(1h-),4x,3hn2=,i2,1x,13(1h-),10x,6hav.en=,
     &        f12.6,13h(keV) tot.in=,1pe11.4,4x,5a1,7h av.ch=,0pf9.3)
 1250  format(4h n1=,i2,1x,13(1h-),4x,3hn2=,i2,1x,13(1h-),10x,6hav.en=,
     &        f12.6,13h(keV) tot.in=,1pe11.4,4x,5a1)
       endif

       lc = lc + 1

       if(idebug.eq.1)then
       if(lc.ge.60) write(iw,850)
       endif
       if(lc.ge.60) lc = 0

 1300  if(idebug.eq.1)then
       write(iw,1400)
 1400  format(1h )
       endif

       lc = lc + 1

       if(idebug.eq.1)then
       if(lc.ge.60) write(iw,850)
       endif
       if(lc.ge.60) lc = 0  

       n = 0
       aa = 0.0d0
       ee = 0.0d0

 1500 continue

*-----------------------------------------------------------------------
      return
      end


************************************************************************
*                                                                      *
      double precision function rate(n1,l1,n2,l2)
*                                                                      *
c  ***  master rate routine -- interprets options and calls other rates
c  ***  point dirac or inputed energies used exclusively in rates ***
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      common /aama00/ ijk
!$OMP THREADPRIVATE(/aama00/)
      common /aama01/ energ, econs, econst, d2p1s, d2p1sm
!$OMP THREADPRIVATE(/aama01/)
      common /aama02/ zsa(3), be(3), bem(3)
!$OMP THREADPRIVATE(/aama02/)
      common /aama03/ k0, k1, k2, k3
!$OMP THREADPRIVATE(/aama03/)
      common /aama04/ nn0(3), nn1(7), nn2(7), nn3(7)
!$OMP THREADPRIVATE(/aama04/)
      common /aama05/ ip1(7), ip2(7), ip3(7), iq1(7), iq2(7), iq3(7)
!$OMP THREADPRIVATE(/aama05/)
      common /aama11/ rr(18), rau, rad, ra(4), rd(4), rsa(4)
!$OMP THREADPRIVATE(/aama11/)
      common /aama13/ e(1000), ai(1000), energya(20,40), m, ia(1000)
!$OMP THREADPRIVATE(/aama13/)
      common /aama17/ yc(4)
!$OMP THREADPRIVATE(/aama17/)

cABE 2018/10/30, add threadprivate aamac
      common /aamac01/ iread,iw,ipunch,iprint
!$OMP THREADPRIVATE(/aamac01/)
      common /aamac03/ pi,coeff,picoef,amassm,amasse,alfa,ne(3),jk(3)
!$OMP THREADPRIVATE(/aamac03/)
      common /aamac05/ m1(7),m2(7),m3(7),idb
!$OMP THREADPRIVATE(/aamac05/)
      common /aamadbg/ idebug

*-----------------------------------------------------------------------

      dimension y(3), iy(3), jp1(7), jp2(7), jp3(7), idm(3), idd(5,5),
     &          idq(4,4), ido(4,4), idr(18)
      dimension r0(3), r1(7), r2(7), r3(7)

      data idm/4hm-1t,4hm-2t,4hm-3t/
      data idd/4hd-rt,4hd-1t,4hd-2t,4hd-3t,4h****,4h****,4hd-1s,4hd-2s,
     &4hd-3s,4h****,2*4h****,4hd-2p,4hd-3p,4h****,3*4h****,4hd-3d,
     &4h****,4*4h****,4h****/
      data idq/4hq-rt,4hq-1t,4hq-2t,4hq-3t,4h****,4hq-1s,4hq-2s,4hq-3s,
     &2*4h****,4hq-2p,4hq-3p,3*4h****,4hq-3d/
      data ido/4h8-rt,4h8-1t,4h8-2t,4h8-3t,4h****,4h8-1s,4h8-2s,4h8-3s,
     &2*4h****,4h8-2p,4h8-3p,3*4h****,4h8-3d/

*-----------------------------------------------------------------------

      irr = 0
      n12 = (n1+n2+1)/2
      rate = 0.0d0
      rau = 0.0d0
      rad = 0.0d0

      do 50 i = 1, 4
       ra(i) = 0.0d0
       rd(i) = 0.0d0
       rsa(i) = 0.0d0
   50 continue

*-----------------------------------------------------------------------

      l = iabs(l1-l2) + 1
      if(l.gt.4) return

      enem = energ
      if(ijk.ne.0) goto 100

      ll1 = 1
      ll2 = 1

      if(l1.eq.0) ll1 = 0
      if(l2.eq.0) ll2 = 0

      lj11 = 2*l1 + 1
      lj12 = lj11 - ll1
      lj21 = 2*l2 + 1
      lj22 = lj21 - ll2

      ene1 = 0.5d0*(energya(n1,lj11) + energya(n1,lj12))
      ene2 = 0.5d0*(energya(n2,lj21) + energya(n2,lj22))

      if(n2.eq.1 .and. d2p1sm.gt.1.0d-20)
     & ene2 = d2p1sm*amasse*1.0d-6 + 0.5d0*(energya(2,2) + energya(2,3))

      enem = (ene2-ene1)*1.0d6/amasse

  100 do 200 i = 1, 3
       iy(i) = 1
       t = enem - bem(i)
cABE change @2015/06/18
       if(t.le.0.0d0) goto 200
c       iy(i) = 0
       y(i) = zsa(i)/dsqrt(t*t+2.000*t)
       if(2.d0*pi*y(i) .le. 700.d0)then
        iy(i) = 0
       endif
  200 continue

      goto (1000,2000,3000,4000), l

*-----------------------------------------------------------------------
*     monopole rate
*-----------------------------------------------------------------------

 1000 if(k0.eq.0) goto 3000

      do 1100 i = 1, k0
       r0(i) = 0.0d0
       nn = nn0(i)
       if(iy(nn).ne.0) goto 1100
       if(y(nn).lt.yc(nn)) goto 1100
       r0(i) = rmon(n1,l1,n2,l2,nn,y(nn))
       rate = rate + r0(i)
       rau = rau + r0(i)
       irr = irr + 1
       rr(irr) = r0(i)
       if(nn.eq.0) goto 1100
       rsa(nn) = rsa(nn) + r0(i)
       idr(irr) = idm(nn)
 1100 continue

      ra(1) = rau

      goto 3000

*-----------------------------------------------------------------------
*     dipole unpenetrated rates
*-----------------------------------------------------------------------

 2000 if(k1.eq.0)  go to 4000

      do 2100 i = 1, k1
       r1(i) = 0.0d0
       nn = nn1(i)                                                      
       mm = nn
       if(mm.eq.0) mm = 1
       jp1(i) = 0
       if(nn.eq.0) goto 2050
       if(iy(nn).ne.0) goto 2100
       if(y(nn).gt.yc(nn)) jp1(i) = ip1(i)
 2050  if(jp1(i).eq.0) r1(i) = rdipu(n1,l1,n2,l2,nn,enem,y(mm),m1(i))
       if(nn.eq.0) goto 2075
       if(jp1(i).ne.0 .and. n12.ge.iq1(i))
     &  r1(i) = rdip(n1,l1,n2,l2,nn,m1(i),y(mm))
       if(jp1(i).ne.0 .and. n12.lt.iq1(i))
     &  r1(i) = rdipu(n1,l1,n2,l2,nn,enem,y(mm),m1(i))
 2075  rate = rate + r1(i)
       if(nn.eq.0) rad = rad + r1(i)
       if(nn.ne.0) rau = rau + r1(i)
       irr = irr + 1
       rr(irr) = r1(i)
       mn = m1(i) + 1
       if(nn.eq.0) mn = 1
       idr(irr) = idd(nn+1,mn)
       if(nn.eq.0) goto 2100
       rsa(nn) = rsa(nn) + r1(i)
 2100 continue

      rd(2) = rad
      ra(2) = rau

      goto 4000

*-----------------------------------------------------------------------
*     quadrupole unpenetrates rates
*-----------------------------------------------------------------------

 3000 if(k2.eq.0 .or. l1+l2.eq.0) goto 5000

       do 3100 i = 1, k2
       r2(i) = 0.0d0
       nn = nn2(i)
       mm = nn
       if(mm.eq.0) mm = 1
       jp2(i) = 0
       if(nn.eq.0) goto 3050
       if(iy(nn).ne.0) goto 3100
       if(y(nn).gt.yc(nn)) jp2(i) = ip2(i)
 3050  if(jp2(i).eq.0) r2(i) = rquau(n1,l1,n2,l2,nn,enem,y(mm),m2(i))
       if(nn.eq.0) goto 3075
       if(jp2(i).ne.0 .and. n12.ge.iq2(i))
     &  r2(i) = rqua(n1,l1,n2,l2,nn,m2(i),y(mm))
       if(jp2(i).ne.0 .and. n12.lt.iq2(i))
     &  r2(i) = rquau(n1,l1,n2,l2,nn,enem,y(mm),m2(i))
 3075  rate = rate+r2(i)
       if(nn.eq.0) rd(3) = rd(3) + r2(i)
       if(nn.eq.0) rad = rad + r2(i)
       if(nn.ne.0) ra(3) = ra(3) + r2(i)
       if(nn.ne.0) rau = rau + r2(i)
       irr = irr + 1
       rr(irr) = r2(i)
       mn = m2(i) + 1
       if(nn.eq.0) mn = 1
       idr(irr) = idq(nn+1,mn)
       if(nn.eq.0) goto 3100
       rsa(nn) = rsa(nn) + r2(i)
 3100 continue

      goto 5000

*-----------------------------------------------------------------------
*     octupole unpenetrated rates
*-----------------------------------------------------------------------

 4000 if(k3.eq.0 .or. l1+l2.eq.1) goto 5000

       do 4100 i = 1, k3
       r3(i) = 0.0d0
       nn = nn3(i)
       mm = nn
       if(mm.eq.0) mm = 1
       jp3(i) = 0
       if(nn.eq.0) goto 4050
       if(iy(nn).ne.0) goto 4100
       if(y(mm).gt.yc(nn)) jp3(i) = ip3(i)
 4050  if(jp3(i).eq.0) r3(i) = roctu(n1,l1,n2,l2,nn,enem,y(mm),m3(i))
       if(nn.eq.0) goto 4075
       if(jp3(i).ne.0 .and. n12.ge.iq3(i))
     &  r3(i)=roct(n1,l1,n2,l2,nn,m3(i),y(mm))
       if(jp3(i).ne.0 .and. n12.lt.iq3(i))
     &  r3(i)=roctu(n1,l1,n2,l2,nn,enem,y(mm),m3(i))
 4075  if(nn.eq.0) rd(4) = rd(4) + r3(i)
       if(nn.eq.0) rad = rad + r3(i)
       if(nn.ne.0) ra(4) = ra(4) + r3(i)
       if(nn.ne.0) rau = rau + r3(i)
       rate = rate + r3(i)
       irr = irr + 1
       rr(irr) = r3(i)
       mn = m3(i) + 1
       if(nn.eq.0) mn = 1
       idr(irr) = ido(nn+1,mn)
       if(nn.eq.0) goto 4100
       rsa(nn) = rsa(nn) + r3(i)
 4100 continue

*-----------------------------------------------------------------------

 5000 if(idb.eq.0) goto 5400
      if(irr.eq.0) goto 5200

      if(idebug.eq.1)then
      write(iw,5100) n1, l1, n2, l2, y, rate, (idr(i),rr(i),i=1,irr)
 5100 format(1x,i2,1h,,i2,3h - ,i2,1h,,i2,3f7.3,1pe10.3 /7(1x,a4,
     &1pe12.4))
      endif

 5200 continue
      if(idebug.eq.1)then
      if(irr.eq.0) write(iw,5300) n1, l1, n2, l2, y
 5300 format(1x,i2,1h,,i2,3h - ,i2,1h,,i2,3f7.3,14h ***no rate***)
      endif

 5400 if(irr.eq.0) return

      if(idebug.eq.1)then
      do 5500 i = 1, irr
       if(rr(i).lt.0.0d0) write(iw,5600) n1, l1, n2, l2, i, rr(i)
 5500 continue
 5600 format(53h *** error *** in internal calculation of transition ,
     &12hrates at n1=,i2,4h l1=,i2,5h, n2=,i2,4h l2=,i2,9h rate no=,i2,
     &7h rate =,1pe13.5,4h ***)
      endif

*-----------------------------------------------------------------------
      return
      end


************************************************************************
*                                                                      *
      double precision function rmon(n1,l1,n2,l2,n,y)
*                                                                      *
c  ***  monopole rate routine (penetration only)                        
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      real*8 matel

      common /aama06/ z, zsk, zsl, zsm, zskz, zslz, zsmz
!$OMP THREADPRIVATE(/aama06/)
      common /aama12/ pop(6)
!$OMP THREADPRIVATE(/aama12/)

cABE 2018/10/30, add threadprivate aamac
      common /aamac03/ pi,coeff,picoef,amassm,amasse,alfa,ne(3),jk(3)
!$OMP THREADPRIVATE(/aamac03/)
      common /aamac08/ coemon(30),expmon(30),ifm(6),jm(10)
!$OMP THREADPRIVATE(/aamac08/)
      common /aamac15/ ij(4),yj(4),jj1(4)
!$OMP THREADPRIVATE(/aamac15/)

*-----------------------------------------------------------------------

      dimension km(10)

*-----------------------------------------------------------------------

      epiy = dexp(pi*y)
      coeffloc = picoef*2.0d0*epiy/(epiy-1.0d0/epiy)

      do 500 jj = 1, 10
       km(jj) = jm(jj)
  500 continue

      if(ij(1).eq.0) goto 900

      do 600 jj=1,10
       if(y.gt.yj(1)) km(jj) = min0(jm(jj),jj1(1))
  600 continue

  900 goto (1000,2000,3000), n

*-----------------------------------------------------------------------

 1000 a1 = 0.0d0

      if(ifm(1).eq.1) return

      do 1100 jj = 1, 3
       if(jj.gt.km(1)) goto 1100
       j = jk(jj)
       b = y**dble(ne(jj))
       a1 = a1 + coemon(jj)*matel(n1,l1,n2,l2,j+2,expmon(jj),1)/b
 1100 continue

      rmon = a1*a1*coeffloc*pop(1)

      return

*-----------------------------------------------------------------------

 2000 a1 = 0.0d0
      a2 = 0.0d0
      a3 = 0.0d0
      yy = (1.0d0+y*y)/(y*y)
      if(ifm(2)+ifm(3).eq.2) return
      do 2200 jj = 1, 3
       j = jk(jj)
       b = y**dble(ne(jj))
       if(ifm(2).eq.1) goto 2100
       if(jj.gt.km(2)) goto 2100
       a1 = a1 + coemon(jj+3)*matel(n1,l1,n2,l2,j+2,expmon(jj+3),2)/b
 2100  if(ifm(3).eq.1) goto 2200
       if(jj.gt.km(3)) goto 2150
       a2 = a2 + coemon(jj+ 6)*matel(n1,l1,n2,l2,j+3,expmon(jj+ 6),2)/b
 2150  if(jj.gt.km(4)) goto 2200
       a3 = a3 + coemon(jj+ 9)*matel(n1,l1,n2,l2,j+4,expmon(jj+ 9),2)/b
 2200 continue

      rmon=((a1+a2)**2*pop(2)+a3*a3*yy*pop(3))*coeffloc

      return

*-----------------------------------------------------------------------

 3000 a1 = 0.0d0
      a2 = 0.0d0
      a3 = 0.0d0
      a4 = 0.0d0
      a5 = 0.0d0
      a6 = 0.0d0

      yy1 = (1.0d0+y*y)/(y*y)
      yy2 = (1.0d0+y*y)*(4.0d0+y*y)/y**4

      if(ifm(4)+ifm(5)+ifm(6).eq.3) return
      do 3300 jj = 1, 3
       j = jk(jj)
       b = y**dble(ne(jj))
       if(ifm(4).eq.1) goto 3100
       if(jj.gt.km(5)) goto 3100
       a1 = a1 + coemon(jj+12)*matel(n1,l1,n2,l2,j+2,expmon(jj+12),3)/b
 3100  if(ifm(5).eq.1) goto 3200
       if(jj.gt.km(6)) goto 3150
       a2 = a2 + coemon(jj+15)*matel(n1,l1,n2,l2,j+3,expmon(jj+15),3)/b
 3150  if(jj.gt.km(7)) goto 3200
       a3 = a3 + coemon(jj+18)*matel(n1,l1,n2,l2,j+4,expmon(jj+18),3)/b
 3200  if(ifm(6).eq.1) goto 3300
       if(jj.gt.km(8)) goto 3230
       a4 = a4 + coemon(jj+21)*matel(n1,l1,n2,l2,j+4,expmon(jj+21),3)/b
 3230  if(jj.gt.km(9)) goto 3260
       a5 = a5 + coemon(jj+24)*matel(n1,l1,n2,l2,j+5,expmon(jj+24),3)/b
 3260  if(jj.gt.km(10)) goto 3300
       a6 = a6 + coemon(jj+27)*matel(n1,l1,n2,l2,j+6,expmon(jj+27),3)/b
 3300 continue

      rmon = ((a1+a2+a3)**2*pop(4) +
     &        yy1*(a4+a5)**2*pop(5) + 
     &        yy2*a6*a6*pop(6) )*coeffloc

*-----------------------------------------------------------------------
      return
      end


************************************************************************
*                                                                      *
      double precision function rdip(n1,l1,n2,l2,n,m,y)
*                                                                      *
c  ***  dipole penetration routine                                      
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      real*8 matel

      common /aama06/ z, zsk, zsl, zsm, zskz, zslz, zsmz
!$OMP THREADPRIVATE(/aama06/)
      common /aama07/ amzz(3)
!$OMP THREADPRIVATE(/aama07/)
      common /aama08/ angd, aid, aidsq
!$OMP THREADPRIVATE(/aama08/)
      common /aama12/ pop(6)
!$OMP THREADPRIVATE(/aama12/)

cABE 2018/10/30, add threadprivate aamac
      common /aamac03/ pi,coeff,picoef,amassm,amasse,alfa,ne(3),jk(3)
!$OMP THREADPRIVATE(/aamac03/)
      common /aamac09/ coedip(42),expdip(42),coed(4),coedp(9),
     &                 ifd(9),jd(14)
!$OMP THREADPRIVATE(/aamac09/)
      common /aamac15/ ij(4),yj(4),jj1(4)
!$OMP THREADPRIVATE(/aamac15/)

*-----------------------------------------------------------------------

      dimension a(5), kd(14)

*-----------------------------------------------------------------------

      etpy = dexp(2.0d0*pi*y)
      extpy = etpy/(etpy-1.0d0)
      yy = y*y
      p = dexp(y*(2.0d0*datan(y/dble(n))-pi))

      do 500 i = 1, 5
       a(i) = 0.0d0
  500 continue

      do 550 jj = 1, 14
       kd(jj) = jd(jj)
  550 continue

      if(ij(2).eq.0) goto 900

      do 600 jj = 1, 14
       if(y.gt.yj(2)) kd(jj) = min0(jd(jj),jj1(2))
  600 continue

  900 goto (1000,2000,3000), n

*-----------------------------------------------------------------------

 1000 a1 = 0.0d0

      if(ifd(1).eq.1) goto 1200
      do 1100 jj = 1, 3
       if(jj.gt.kd(1)) goto 1100
       j = jk(jj)
       b = y**dble(ne(jj))
       a1 = a1 + coedip(jj)*matel(n1,l1,n2,l2,j+3,expdip(jj),1)/b
 1100 continue                                                        

 1200 yf = yy/(1.0d0+yy)
      yg = (1.0d0+yy)/yy

      rdip = coedp(1)*extpy*picoef*angd*yg*(yf*p*aid*amzz(1)-a1)**2
      rdip = rdip*pop(1)

      return

*-----------------------------------------------------------------------

 2000 a1 = 0.0d0
      a2 = 0.0d0
      a3 = 0.0d0

      if(m.eq.2) goto 2300
      if(ifd(2).eq.1) goto 2200

      do 2100 jj = 1, 3
       j = jk(jj)
       b = y**dble(ne(jj))
       if(jj.gt.kd(2)) goto 2050
       a1 = a1 + coedip(jj+3)*matel(n1,l1,n2,l2,j+3,expdip(jj+3),2)/b
 2050  if(jj.gt.kd(3)) goto 2100
       a1 = a1 + coedip(jj+6)*matel(n1,l1,n2,l2,j+4,expdip(jj+6),2)/b
 2100 continue

 2200 yf = yy/(4.0d0+yy)
      yg = (1.0d0+yy)/yy

      a(1) = coedp(2)*pop(2)*extpy*yg*(yf*p*aid*amzz(2)-a1)**2

 2300 if(m.eq.1) goto 2600

      if(ifd(3)+ifd(4).eq.2) goto 2500

      do 2400 jj = 1, 3
       j = jk(jj)
       b = y**dble(ne(jj))
       if(jj.gt.kd(4)) goto 2350
       a2 = a2 + coedip(jj+ 9)*matel(n1,l1,n2,l2,j+3,expdip(jj+ 9),2)/b
 2350  if(jj.gt.kd(5)) goto 2400
       a3 = a3 + coedip(jj+12)*matel(n1,l1,n2,l2,j+5,expdip(jj+12),2)/b
 2400 continue

 2500 yf = yy/(4.0d0+yy)

      a(2) = coedp(3)*pop(3)*extpy*(yf*p*aid*amzz(2)-a2)**2
      yf = yy*yy/(4.0d0+yy)**2
      yg = (1.0d0+yy)*(4.0d0+yy)/(yy*yy)
      a(3) = coedp(4)*pop(3)*extpy*yg*(yf*p*aid*amzz(2)-a3)**2

 2600 at = a(1) + a(2) + a(3)
      rdip = picoef*angd*at

      return

*-----------------------------------------------------------------------

 3000 a1 = 0.0d0
      a2 = 0.0d0
      a3 = 0.0d0
      a4 = 0.0d0
      a5 = 0.0d0

      if(m.gt.1) goto 3300
      if(ifd(5).eq.1) goto 3200

      do 3100 jj = 1, 3
       j = jk(jj)
       b = y**dble(ne(jj))
       if(jj.gt.kd(6)) goto 3030
       a1 = a1 + coedip(jj+15)*matel(n1,l1,n2,l2,j+3,expdip(jj+15),3)/b
 3030  if(jj.gt.kd(7)) goto 3060
       a1 = a1 + coedip(jj+18)*matel(n1,l1,n2,l2,j+4,expdip(jj+18),3)/b
 3060  if(jj.gt.kd(8)) goto 3100
       a1 = a1 + coedip(jj+21)*matel(n1,l1,n2,l2,j+5,expdip(jj+21),3)/b
 3100 continue

 3200 yf = yy*(27.0d0+7.0d0*yy)/(9.0d0+yy)**2
      yg = (1.0d0+yy)/yy
      a(1) = coedp(5)*pop(4)*extpy*yg*(yf*p*aid*amzz(3)-a1)**2

 3300 if(m.ne.0 .and. m.ne.2) goto 3600
      if(ifd(6)+ifd(7).eq.2) goto 3500

      do 3400 jj = 1, 3
       j = jk(jj)
       b = y**dble(ne(jj))
       if(jj.gt.kd(9)) goto 3325
       a2 = a2 + coedip(jj+24)*matel(n1,l1,n2,l2,j+3,expdip(jj+24),3)/b
 3325  if(jj.gt.kd(10)) goto 3350
       a2 = a2 + coedip(jj+27)*matel(n1,l1,n2,l2,j+4,expdip(jj+27),3)/b
 3350  if(jj.gt.kd(11)) goto 3375
       a3 = a3 + coedip(jj+30)*matel(n1,l1,n2,l2,j+5,expdip(jj+30),3)/b
 3375  if(jj.gt.kd(12)) goto 3400
       a3 = a3 + coedip(jj+33)*matel(n1,l1,n2,l2,j+6,expdip(jj+33),3)/b
 3400 continue

 3500 yf = yy*(3.0d0+yy)/(9.0d0+yy)**2

      a(2) = coedp(6)*pop(5)*extpy*(yf*p*aid*amzz(3)-a2)**2
      yf = yy*yy/(9.0d0+yy)**2
      yg = (1.0d0+yy)*(4.0d0+yy)/(yy*yy)
      a(3) = coedp(7)*pop(5)*extpy*yg*(yf*p*aid*amzz(3)-a3)**2

 3600 if(m.ne.0.and.m.ne.3) goto 3900
      if(ifd(8)+ifd(9).eq.2) goto 3800

      do 3700 jj = 1, 3
       j = jk(jj)
       b = y**dble(ne(jj))
       if(jj.gt.kd(13)) goto 3650
       a4 = a4 + coedip(jj+36)*matel(n1,l1,n2,l2,j+5,expdip(jj+36),3)/b
 3650  if(jj.gt.kd(14)) goto 3700
       a5 = a5 + coedip(jj+39)*matel(n1,l1,n2,l2,j+7,expdip(jj+39),3)/b
 3700 continue

 3800 yf = yy*yy/(9.0d0+yy)**2
      yg = (1.0d0+yy)/yy
      a(4) = coedp(8)*pop(6)*extpy*yg*(yf*p*aid*amzz(3)-a4)**2
      yf = (yy/(9.0d0+yy))**3
      yg = (1.0d0+yy)*(4.0d0+yy)*(9.0d0+yy)/yy**3
      a(5) = coedp(9)*pop(5)*extpy*yg*(yf*p*aid*amzz(3)-a5)**2

 3900 at = a(1) + a(2) + a(3) + a(4) + a(5)
      rdip = picoef*angd*at

*-----------------------------------------------------------------------
      return
      end


************************************************************************
*                                                                      *
      double precision function rqua(n1,l1,n2,l2,n,m,y)
*                                                                      *
c  ***quadrupole penetration routine
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      real*8 matel

      common /aama06/ z, zsk, zsl, zsm, zskz, zslz, zsmz
!$OMP THREADPRIVATE(/aama06/)
      common /aama07/ amzz(3)
!$OMP THREADPRIVATE(/aama07/)
      common /aama09/ angq, aiq, aiqsq
!$OMP THREADPRIVATE(/aama09/)
      common /aama12/ pop(6)
!$OMP THREADPRIVATE(/aama12/)

cABE 2018/10/30, add threadprivate aamac
      common /aamac03/ pi,coeff,picoef,amassm,amasse,alfa,ne(3),jk(3)
!$OMP THREADPRIVATE(/aamac03/)
      common /aamac10/ coequa(45),expqua(45),coeq(11),ifq(10),jq(15)
!$OMP THREADPRIVATE(/aamac10/)
      common /aamac15/ ij(4),yj(4),jj1(4)
!$OMP THREADPRIVATE(/aamac15/)

*-----------------------------------------------------------------------

      dimension a(6), kq(15)

*-----------------------------------------------------------------------

      etpy = dexp(2.0d0*pi*y)
      extpy = etpy/(etpy-1.0d0)
      yy = y*y
      p = dexp(y*(2.0d0*datan(y/dble(n))-pi))

      do 500 i = 1, 6
       a(i)=0.0d0
  500 continue

      do 550 jj = 1, 15
       kq(jj) = jq(jj)
  550 continue

      if(ij(3).eq.0) goto 900

      do 600 jj = 1, 15
       if(y.gt.yj(3)) kq(jj) = min0(jq(jj),jj1(3))
  600 continue

  900 go to (1000,2000,3000), n

*-----------------------------------------------------------------------

 1000 a1 = 0.0d0

      if(ifq(1).eq.1) goto 1200

      do 1100 jj = 1, 3
       if(jj.gt.kq(1)) goto 1100
       j = jk(jj)
       b = y**dble(ne(jj))
       a1 = a1 + coequa(jj)*matel(n1,l1,n2,l2,j+4,expqua(jj),1)/b
 1100 continue

 1200 yf = yy/(4.0d0+yy)
      yg = (1.0d0+yy)*(4.0d0+yy)/(yy*yy)
      pf = 9.0d0*p-1.0d0

      rqua = coeq(2)*extpy*picoef*angq*yg*(yf*pf*aiq*amzz(1)**2-a1)**2
      rqua = rqua*pop(1)

      return

*-----------------------------------------------------------------------

 2000 a1 =0.0d0
      a2 = 0.0d0
      a3 = 0.0d0

      if(m.eq.2) goto 2300
      if(ifq(2).eq.1) goto 2200

      do 2100 jj = 1, 3
       j = jk(jj)
       b = y**dble(ne(jj))
       if(jj.gt.kq(2)) goto 2050
       a1 = a1 + coequa(jj+3)*matel(n1,l1,n2,l2,j+4,expqua(jj+3),2)/b
 2050  if(jj.gt.kq(3)) goto 2100
       a1 = a1 + coequa(jj+6)*matel(n1,l1,n2,l2,j+5,expqua(jj+6),2)/b
 2100 continue

 2200 yf = yy/(1.0d0+yy)
      yg = (1.0d0+yy)*(4.0d0+yy)/(yy*yy)
      pf = 9.0d0*(4.0d0+5.0d0*yy)/(4.0d0+yy)*p - 1.0d0
      a(1) = coeq(3)*pop(2)*extpy*yg*(yf*pf*aiq*amzz(2)**2-a1)**2

 2300 if(m.eq.1) goto 2600
      if(ifq(3)+ifq(4).eq.2) goto 2500

      do 2400 jj = 1, 3
       j = jk(jj)
       b = y**dble(ne(jj))
       if(jj.gt.kq(4)) goto 2350
       a2 = a2 + coequa(jj+9)*matel(n1,l1,n2,l2,j+4,expqua(jj+9),2)/b
 2350  if(jj.gt.kq(5)) goto 2400
       a3 = a3 + coequa(jj+12)*matel(n1,l1,n2,l2,j+6,expqua(jj+12),2)/b
 2400 continue

 2500 yf = yy/(1.0d0+yy)
      yg = (1.0d0+yy)/yy
      pf = 3.0d0*p+1.0d0
      a(2) = coeq(4)*pop(3)*extpy*yg*(yf*pf*aiq*amzz(2)**2-a2)**2
      yf = yy*yy/((1.0d0+yy)*(9.0d0+yy))
      yg = (1.0d0+yy)*(4.0d0+yy)*(9.0d0+yy)/yy**3
      pf = (68.0d0+77.0d0*yy)/(4.0d0+yy)*p-1.0d0
      a(3) = coeq(5)*pop(3)*extpy*yg*(yf*pf*aiq*amzz(2)**2-a3)**2

 2600 at = a(1) + a(2) + a(3)
      rqua = picoef*angq*at
      return

*-----------------------------------------------------------------------

 3000 a1 = 0.0d0
      a2 = 0.0d0
      a3 = 0.0d0
      a4 = 0.0d0
      a5 = 0.0d0
      a6 = 0.0d0

      if(m.gt.1) goto 3300
      if(ifq(5).eq.1) goto 3200

      do 3100 jj = 1, 3
       j = jk(jj)
       b = y**dble(ne(jj))
       if(jj.gt.kq(6)) goto 3030
       a1 = a1 + coequa(jj+15)*matel(n1,l1,n2,l2,j+4,expqua(jj+15),3)/b
 3030  if(jj.gt.kq(7)) goto 3060
       a1 = a1 + coequa(jj+18)*matel(n1,l1,n2,l2,j+5,expqua(jj+18),3)/b
 3060  if(jj.gt.kq(8)) goto 3100
       a1 = a1 + coequa(jj+21)*matel(n1,l1,n2,l2,j+6,expqua(jj+21),3)/b
 3100 continue

 3200 yf = yy*(9.0d0+yy)/((1.0d0+yy)*(4.0d0+yy))
      yg = (1.0d0+yy)*(4.0d0+yy)/(yy*yy)
      pf = (729.0d0+1134.0d0*yy+277.0d0*yy*yy)/(9.0d0+yy)**2*p-1.0d0
      a(1) = coeq(6)*pop(4)*extpy*yg*(yf*pf*aiq*amzz(3)**2-a1)**2

 3300 if(m.ne.0 .and. m.ne.2) goto 3600
      if(ifq(6)+ifq(7).eq.2) goto 3500

      do 3400 jj = 1, 3
       j = jk(jj)
       b = y**dble(ne(jj))
       if(jj.gt.kq(9)) goto 3325
       a2 = a2 + coequa(jj+24)*matel(n1,l1,n2,l2,j+4,expqua(jj+24),3)/b
 3325  if(jj.gt.kq(10)) goto 3350
       a2 = a2 + coequa(jj+27)*matel(n1,l1,n2,l2,j+5,expqua(jj+27),3)/b
 3350  if(jj.gt.kq(11)) goto 3375
       a3 = a3 + coequa(jj+30)*matel(n1,l1,n2,l2,j+6,expqua(jj+30),3)/b
 3375  if(jj.gt.kq(12)) goto 3400
       a3 = a3 + coequa(jj+33)*matel(n1,l1,n2,l2,j+7,expqua(jj+33),3)/b
 3400 continue

 3500 yf = yy/(1.0d0+yy)
      yg = (1.0d0+yy)/yy
      pf = (27.0d0+11.0d0*yy)/(9.0d0+yy)*p+1.0d0
      a(2) = coeq(7)*pop(5)*extpy*yg*(yf*pf*aiq*amzz(3)**2-a2)**2

      yf = yy*yy/((1.0d0+yy)*(4.0d0+yy))
      yg = (1.0d0+yy)*(4.0d0+yy)*(9.0d0+yy)/yy**3
      pf = (1377.0d0+1944.0d0*yy+439.0d0*yy*yy)/(9.0d0+yy)**2*p-1.0d0
      a(3) = coeq(8)*pop(5)*extpy*yg*(yf*pf*aiq*amzz(3)**2-a3)**2

 3600 if(m.ne.0.and.m.ne.3) goto 3900
      if(ifq(8)+ifq(9)+ifq(10).eq.3) goto 3800

      do 3700 jj = 1, 3
       j = jk(jj)
       b = y**dble(ne(jj))
       if(jj.gt.kq(13)) goto 3630
       a4 = a4 + coequa(jj+36)*matel(n1,l1,n2,l2,j+4,expqua(jj+36),3)/b
 3630  if(jj.gt.kq(14)) goto 3660
       a5 = a5 + coequa(jj+39)*matel(n1,l1,n2,l2,j+6,expqua(jj+39),3)/b
 3660  if(jj.gt.kq(15)) goto 3700
       a6 = a6 + coequa(jj+42)*matel(n1,l1,n2,l2,j+8,expqua(jj+42),3)/b
 3700 continue

 3800 yf = yy/(9.0d0+yy)
      a(4) = coeq(9)*pop(6)*extpy*(yf*p*aiq*amzz(3)**2-a4)**2
      yf = yy*yy/((1.0d0+yy)*(4.0d0+yy))
      yg = (1.0d0+yy)*(4.0d0+yy)/(yy*yy)
      pf = (63.0d0+47.0d0*yy)/(9.0d0+yy)*p+1.0d0
      a(5) = coeq(10)*pop(6)*extpy*yg*(yf*pf*aiq*amzz(3)**2-a5)**2
      yf = yy**3/((1.0d0+yy)*(4.0d0+yy)*(16.0d0+yy))
      yg = (1.0d0+yy)*(4.0d0+yy)*(9.0d0+yy)*(16.0d0+yy)/yy**4
      pf = (10773.0d0+14580.0d0*yy+3167.0d0*yy*yy)/(9.0d0+yy)**2*p-5.0d0
      a(6) = coeq(11)*pop(6)*extpy*yg*(yf*pf*aiq*amzz(3)**2-a6)**2

 3900 at = a(1) + a(2) + a(3) + a(4) + a(5) + a(6)
      rqua = picoef*angq*at

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      double precision function roct(n1,l1,n2,l2,n,m,y)
*                                                                      *
c  ***  octupole penetration routine                                    
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      real*8 matel

      common /aama06/ z, zsk, zsl, zsm, zskz, zslz, zsmz
!$OMP THREADPRIVATE(/aama06/)
      common /aama07/ amzz(3)
!$OMP THREADPRIVATE(/aama07/)
      common /aama10/ ango, aio, aiosq
!$OMP THREADPRIVATE(/aama10/)
      common /aama12/ pop(6)
!$OMP THREADPRIVATE(/aama12/)

cABE 2018/10/30, add threadprivate aamac
      common /aamac03/ pi,coeff,picoef,amassm,amasse,alfa,ne(3),jk(3)
!$OMP THREADPRIVATE(/aamac03/)
      common /aamac11/ coeoct(45),expoct(45),coeo(11),ifo(10),jo(15)
!$OMP THREADPRIVATE(/aamac11/)
      common /aamac15/ ij(4),yj(4),jj1(4)
!$OMP THREADPRIVATE(/aamac15/)

*-----------------------------------------------------------------------

      dimension a(6), ko(15)

*-----------------------------------------------------------------------

      etpy = dexp(2.0d0*pi*y)
      extpy = etpy/(etpy-1.0d0)
      yy = y*y
      p = dexp(y*(2.0d0*datan(y/dble(n))-pi))

      do 500 i = 1, 6
       a(i) = 0.0d0
  500 continue

      do 550 jj = 1, 15
       ko(jj) = jo(jj)
  550 continue

      if(ij(4).eq.0) goto 900

      do 600 jj=1,15
       if(y.gt.yj(4)) ko(jj) = min0(jo(jj),jj1(4))
  600 continue

  900 goto (1000,2000,3000), n

*-----------------------------------------------------------------------

 1000 a1 = 0.0d0

      if(ifo(1).eq.1) goto 1200

      do 1100 jj = 1, 3
       if(jj.gt.ko(1)) goto 1100
       j = jk(jj)
       b = y**dble(ne(jj))
       a1 = a1 + coeoct(jj)*matel(n1,l1,n2,l2,j+5,expoct(jj),1)/b
 1100 continue

 1200 yf = yy*(3.0d0+2.0d0*yy)/((4.0d0+yy)*(9.0d0+yy))
      yg = (1.0d0+yy)*(4.0d0+yy)*(9.0d0+yy)/yy**3
      pf = 15.0d0*(1.0d0+yy)/(3.0d0+2.0d0*yy)*p-1.0d0
      roct = coeo(2)*picoef*ango*extpy*yg*(yf*pf*aio*amzz(1)**3-a1)**2
      roct = roct*pop(1)

      return

*-----------------------------------------------------------------------

 2000 a1 =0.0d0
      a2 = 0.0d0
      a3 = 0.0d0

      if(m.eq.2) goto 2300
      if(ifo(2).eq.1) goto 2200

      do 2100 jj = 1, 3
       j = jk(jj)
       b = y**dble(ne(jj))
       if(jj.gt.ko(2)) goto 2050
       a1 = a1 + coeoct(jj+3)*matel(n1,l1,n2,l2,j+5,expoct(jj+3),2)/b
 2050  if(jj.gt.ko(3)) goto 2100
       a1 = a1 + coeoct(jj+6)*matel(n1,l1,n2,l2,j+6,expoct(jj+6),2)/b
 2100 continue

 2200 yf = yy*(6.0d0+yy)/((1.0d0+yy)*(9.0d0+yy))
      yg = (1.0d0+yy)*(4.0d0+yy)*(9.0d0+yy)/yy**3
      pf = 15.0d0*(2.0d0+3.0d0*yy)/(6.0d0+yy)*p-1.0d0
      a(1) = coeo(3)*pop(2)*extpy*yg*(yf*pf*aio*amzz(2)**3-a1)**2

 2300 if(m.eq.1) goto 2600
      if(ifo(3)+ifo(4).eq.2) goto 2500

      do 2400 jj = 1,3
       j = jk(jj)
       b = y**dble(ne(jj))
       if(jj.gt.ko(4)) goto 2350
       a2 = a2 + coeoct(jj+ 9)*matel(n1,l1,n2,l2,j+5,expoct(jj+ 9),2)/b
 2350  if(jj.gt.ko(5)) goto 2400
       a3 = a3 + coeoct(jj+12)*matel(n1,l1,n2,l2,j+7,expoct(jj+12),2)/b
 2400 continue

 2500 yf = yy/(1.0d0+yy)
      yg = (1.0d0+yy)*(4.0d0+yy)/(yy*yy)
      pf = 3.0d0*p+1.0d0
      a(2) = coeo(4)*pop(3)*extpy*yg*(yf*pf*aio*amzz(2)**3-a2)**2
      yf = yy*yy*(68.0d0+13.0d0*yy)/((1.0d0+yy)*(9.0d0+yy)*(16.0d0+yy))
      yg = (1.0d0+yy)*(4.0d0+yy)*(9.0d0+yy)*(16.0d0+yy)/yy**4
      pf = 5.0d0*(116.0d0+149.0d0*yy)/(68.0d0+13.0d0*yy)*p-1.0d0
      a(3) = coeo(5)*pop(3)*extpy*yg*(yf*pf*aio*amzz(2)**3-a3)**2

 2600 at = a(1)+a(2)+a(3)

      roct = picoef*ango*at

      return

*-----------------------------------------------------------------------

 3000 a1 = 0.0d0
      a2 = 0.0d0
      a3 = 0.0d0
      a4 = 0.0d0
      a5 = 0.0d0
      a6 = 0.0d0

      if(m.gt.1) goto 3300
      if(ifo(5).eq.1) goto 3200

      do 3100 jj = 1, 3
       j = jk(jj)
       b = y**dble(ne(jj))
       if(jj.gt.ko(6)) goto 3030
       a1 = a1 + coeoct(jj+15)*matel(n1,l1,n2,l2,j+5,expoct(jj+15),3)/b
 3030  if(jj.gt.ko(7)) goto 3060
       a1 = a1 + coeoct(jj+18)*matel(n1,l1,n2,l2,j+6,expoct(jj+18),3)/b
 3060  if(jj.gt.ko(8)) goto 3100
       a1 = a1 + coeoct(jj+21)*matel(n1,l1,n2,l2,j+7,expoct(jj+21),3)/b
 3100 continue

 3200 yf = yy*(27.0d0+2.0d0*yy)/((1.0d0+yy)*(4.0d0+yy))
      yg = (1.0d0+yy)*(4.0d0+yy)*(9.0d0+yy)/yy**3
      pf = 2.5d0*(405.0d0+900.0d0*yy+254.0d0*yy*yy)/
     &     ((9.0d0+yy)*(27.0d0+2.0d0*yy))*p-1.0d0
      a(1) = coeo(6)*pop(4)*extpy*yg*(yf*pf*aio*amzz(3)**3-a1)**2

 3300 if(m.ne.0 .and. m.ne.2) goto 3600
      if(ifo(6)+ifo(7).eq.2) goto 3500

      do 3400 jj = 1, 3
       j = jk(jj)
       b = y**dble(ne(jj))
       if(jj.gt.ko(9)) goto 3325
       a2 = a2 + coeoct(jj+24)*matel(n1,l1,n2,l2,j+5,expoct(jj+24),3)/b
 3325  if(jj.gt.ko(10)) goto 3350
       a2 = a2 + coeoct(jj+27)*matel(n1,l1,n2,l2,j+6,expoct(jj+27),3)/b
 3350  if(jj.gt.ko(11)) goto 3375
       a3 = a3 + coeoct(jj+30)*matel(n1,l1,n2,l2,j+7,expoct(jj+30),3)/b
 3375  if(jj.gt.ko(12)) goto 3400
       a3 = a3 + coeoct(jj+33)*matel(n1,l1,n2,l2,j+8,expoct(jj+33),3)/b
 3400 continue

 3500 yf = yy*(9.0d0+2.0d0*yy)/((1.0d0+yy)*(4.0d0+yy))
      yg = (1.0d0+yy)*(4.0d0+yy)/(yy*yy)
      pf = (27.0d0+13.0d0*yy)/(9.0d0+2.0d0*yy)*p+1.0d0
      a(2) = coeo(7)*pop(5)*extpy*yg*(yf*pf*aio*amzz(3)**3-a2)**2
      yf = yy*yy*(153.0d0+13.0d0*yy)/((1.0d0+yy)*(4.0d0+yy)*(16.0d0+yy))
      yg = (1.0d0+yy)*(4.0d0+yy)*(9.0d0+yy)*(16.0d0+yy)/yy**4
      pf = 5.0d0*(2349.0d0+3744.0d0*yy+947.0d0*yy*yy)/
     &     ((9.0d0+yy)*(153.0d0+13.0d0*yy))*p-1.0d0
      a(3) = coeo(8)*pop(5)*extpy*yg*(yf*pf*aio*amzz(3)**3-a3)**2

 3600 if(m.ne.0 .and. m.ne.3) goto 3900
      if(ifo(8)+ifo(9)+ifo(10).eq.3) goto 3800

      do 3700 jj = 1, 3
       j = jk(jj)
       b = y**dble(ne(jj))
       if(jj.gt.ko(13)) goto 3630
       a4 = a4 + coeoct(jj+36)*matel(n1,l1,n2,l2,j+5,expoct(jj+36),3)/b
 3630  if(jj.gt.ko(14)) goto 3660
       a5 = a5 + coeoct(jj+39)*matel(n1,l1,n2,l2,j+7,expoct(jj+39),3)/b
 3660  if(jj.gt.ko(15)) goto 3700
       a6 = a6 + coeoct(jj+42)*matel(n1,l1,n2,l2,j+9,expoct(jj+42),3)/b
 3700 continue

 3800 yf = yy/(1.0d0+yy)
      yg = (1.0d0+yy)/yy
      pf = 2.0d0*p+1.0d0
      a(4) = coeo(9)*pop(6)*extpy*yg*(yf*pf*aio*amzz(3)**3-a4)**2
      yf = yy*yy/((1.0d0+yy)*(4.0d0+yy))
      yg = (1.0d0+yy)*(4.0d0+yy)*(9.0d0+yy)/yy**3
      pf = (63.0d0+47.0d0*yy)/(9.0d0+yy)*p+1.0d0
      a(5) = coeo(10)*pop(6)*extpy*yg*(yf*pf*aio*amzz(3)**3-a5)**2
      yf = yy**3*(11.0d0+yy)/
     &     ((1.0d0+yy)*(4.0d0+yy)*(16.0d0+yy)*(25.0d0+yy))
      yg = (1.0d0+yy)*(4.0d0+yy)*(9.0d0+yy)*(16.0d0+yy)*(25.0d0+yy)/
     &     yy**5
      pf = (1251.0d0+1850.0d0*yy+439.0d0*yy*yy)/
     &     ((9.0d0+yy)*(11.0d0+yy))*p - 1.0d0
      a(6) = coeo(11)*pop(6)*extpy*yg*(yf*pf*aio*amzz(3)**3-a6)**2

 3900 at = a(1)+a(2)+a(3)+a(4)+a(5)+a(6)

      roct = picoef*ango*at

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      double precision function rdipu(n1,l1,n2,l2,n,enem,y,mm)
*                                                                      *
c  ***  dipole unpenetrated rates routine                               
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      real*8 matelu

      common /aama06/ z, zsk, zsl, zsm, zskz, zslz, zsmz
!$OMP THREADPRIVATE(/aama06/)
      common /aama07/ amzz(3)
!$OMP THREADPRIVATE(/aama07/)
      common /aama08/ angd, aid, aidsq
!$OMP THREADPRIVATE(/aama08/)
      common /aama12/ pop(6)
!$OMP THREADPRIVATE(/aama12/)

cABE 2018/10/30, add threadprivate aamac
      common /aamac03/ pi,coeff,picoef,amassm,amasse,alfa,ne(3),jk(3)
!$OMP THREADPRIVATE(/aamac03/)
      common /aamac09/ coedip(42),expdip(42),coed(4),coedp(9),
     &                 ifd(9),jd(14)
!$OMP THREADPRIVATE(/aamac09/)

*-----------------------------------------------------------------------

      if(n.gt.0) goto 1000

      aid = matelu(n1,l1,n2,l2,1)
      aidsq = aid*aid
      l = (l1+l2+1)/2
      angd = dble((2*l2+1)*l)/dble((2*l-1)*(2*l+1))
      rdipu = coed(1)*coeff/(z*z*amassm*amassm)*(enem/alfa)**3
     &        *angd*aidsq
      return

*-----------------------------------------------------------------------

 1000 epiy = dexp(pi*y)
      expiy = 1.0d0/(epiy-1.0d0/epiy)
      yy = y*y
      p2 = dexp(y*(4.0d0*datan(y/dble(n))-pi))
      m = mm + 1

      goto (2000,3000,4000), n

*-----------------------------------------------------------------------

 2000 yf = yy/(1.0d0+yy)
      rdipu = coed(2)*pop(1)*picoef*angd*p2*expiy*yf*aidsq*amzz(1)**2

      return

*-----------------------------------------------------------------------

 3000 goto (3100,3200,3300), m

 3100 yf = yy*(4.0d0+3.0d0*yy)*(4.0d0+5.0d0*yy)/(4.0d0+yy)**3*
     &     (pop(2)+3.0d0*pop(3))/4.0d0
      goto 3400

 3200 yf = 4.0d0*yy*(1.0d0+yy)/(4.0d0+yy)**2*pop(2)
      goto 3400

 3300 yf = yy*yy*(12.0d0+11.0d0*yy)/(4.0d0+yy)**3*pop(3)

 3400 rdipu = coed(3)*picoef*angd*p2*expiy*yf*aidsq*amzz(2)**2

      return

*-----------------------------------------------------------------------

 4000 go to (4100,4200,4300,4400),m

 4100 yf = yy*(81.0d0+78.0d0*yy+13.0d0*yy*yy)*
     &     (81.0d0+126.0d0*yy+29.0d0*yy*yy)/
     &     (9.0d0+yy)**5*(pop(4)+3.0d0*pop(5)+5.0d0*pop(6))/9.0d0
      goto 4500                                                      

 4200 yf = yy*(1.0d0+yy)*(27.0d0+7.0d0*yy)**2/(9.0d0+yy)**4*pop(4)
      goto 4500

 4300 yf = 8.0d0*yy*yy*(81.0d0+96.0d0*yy+19.0d0*yy*yy)/
     &     (9.0d0+yy)**4*pop(5)
      goto 4500

 4400 yf = 16.0d0*yy**3*(45.0d0+11.0d0*yy)*(1.0d0+yy)/(9.0d0+yy)**5*
     &     pop(6)

 4500 rdipu = coed(4)*picoef*angd*p2*expiy*yf*aidsq*amzz(3)**2

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      double precision function rquau(n1,l1,n2,l2,n,enem,y,m)
*                                                                      *
c  ***  quadrupole unpenetrates rates routine
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      real*8 matelu

      common /aama06/ z, zsk, zsl, zsm, zskz, zslz, zsmz
!$OMP THREADPRIVATE(/aama06/)
      common /aama07/ amzz(3)
!$OMP THREADPRIVATE(/aama07/)
      common /aama09/ angq, aiq, aiqsq
!$OMP THREADPRIVATE(/aama09/)
      common /aama12/ pop(6)
!$OMP THREADPRIVATE(/aama12/)

cABE 2018/10/30, add threadprivate aamac
      common /aamac03/ pi,coeff,picoef,amassm,amasse,alfa,ne(3),jk(3)
!$OMP THREADPRIVATE(/aamac03/)
      common /aamac10/ coequa(45),expqua(45),coeq(11),ifq(10),jq(15)
!$OMP THREADPRIVATE(/aamac10/)

*-----------------------------------------------------------------------

      dimension a(6)

*-----------------------------------------------------------------------

      if(n.gt.0) goto 1000

      aiq = matelu(n1,l1,n2,l2,2)
      aiqsq = aiq*aiq
      l = (l1+l2)/2

      angq = 1.5d0
      if(l2.eq.l1) angq=1.0d0

      angq = angq*dble((2*l2+1)*l*(l+1))/dble((2*l-1)*(2*l+1)*(2*l+3))
      rquau = coeq(1)*coeff/(z*amassm)**4*(enem/alfa)**5*angq*aiqsq

      return

*-----------------------------------------------------------------------

 1000 etpy = dexp(2.0d0*pi*y)
      extpy = etpy/(etpy-1.0d0)
      yy = y*y
      p = dexp(y*(2.0d0*datan(y/dble(n))-pi))

      do 1100 i = 1, 6
       a(i) = 0.0d0
 1100 continue

      goto (2000,3000,4000), n

*-----------------------------------------------------------------------

 2000 yf = (1.0d0+yy)/(4.0d0+yy)
      pf = (9.0d0*p-1.0d0)**2
      rquau = coeq(2)*pop(1)*picoef*angq*pf*extpy*yf*aiqsq*amzz(1)**4

      return

*-----------------------------------------------------------------------

 3000 if(m.eq.2) goto 3100
      yf = (4.0d0+yy)/(1.0d0+yy)
      pf = (9.0d0*(4.0d0+5.0d0*yy)/(4.0d0+yy)*p-1.0d0)**2
      a(1) = coeq(3)*pop(2)*yf*pf

 3100 if(m.eq.1) goto 3200

      yf = yy/(1.0d0+yy)
      pf = (3.0d0*p+1.0d0)**2
      a(2) = coeq(4)*pop(3)*yf*pf
      yf = yy*(4.0d0+yy)/((1.0d0+yy)*(9.0d0+yy))
      pf = ((68.0d0+77.0d0*yy)/(4.0d0+yy)*p-1.0d0)**2
      a(3) = coeq(5)*pop(3)*yf*pf

 3200 at = a(1) + a(2) + a(3)
      rquau = picoef*angq*at*extpy*aiqsq*amzz(2)**4

      return

*-----------------------------------------------------------------------

 4000 if(m.gt.1) goto 4100

      yf = (9.0d0+yy)**2/((1.0d0+yy)*(4.0d0+yy))
      pf = ((729.0d0+1134.0d0*yy+277.0d0*yy*yy)/
     &     (9.0d0+yy)**2*p-1.0d0)**2
      a(1) = coeq(6)*pop(4)*yf*pf

 4100 if(m.ne.0.and.m.ne.2)  go to 4200

      yf = yy/(1.0d0+yy)            
      pf = ((27.0d0+11.0d0*yy)/(9.0d0+yy)*p+1.0d0)**2
      a(2) = coeq(7)*pop(5)*yf*pf
      yf = yy*(9.0d0+yy)/((1.0d0+yy)*(4.0d0+yy))
      pf = ((1377.0d0+1944.0d0*yy+439.0d0*yy*yy)/
     &     (9.0d0+yy)**2*p-1.0d0)**2
      a(3) = coeq(8)*pop(5)*yf*pf

 4200 if(m.ne.0 .and. m.ne.3) goto 4300

      yf = yy*yy/(9.0d0+yy)**2
      pf = p**2
      a(4) = coeq(9)*pop(6)*yf*pf
      yf = yy*yy/((1.0d0+yy)*(4.0d0+yy))
      pf = ((63.0d0+47.0d0*yy)/(9.0d0+yy)*p+1.0d0)**2
      a(5) = coeq(10)*pop(6)*yf*pf
      yf = yy*yy*(9.0d0+yy)/((1.0d0+yy)*(4.0d0+yy)*(16.0d0+yy))
      pf = ((10773.0d0+14580.0d0*yy+3167.0d0*yy*yy)/
     &     (9.0d0+yy)**2*p-5.0d0)**2
      a(6) = coeq(11)*pop(6)*yf*pf

 4300 at = a(1) + a(2) + a(3) +a (4) + a(5) + a(6)

      rquau = picoef*angq*at*extpy*aiqsq*amzz(3)**4.0d0

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      double precision function roctu(n1,l1,n2,l2,n,enem,y,m)
*                                                                      *
c  ***  octupole unpenetrated rates routine                             
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      real*8 matelu

      common /aama06/ z, zsk, zsl, zsm, zskz, zslz, zsmz
!$OMP THREADPRIVATE(/aama06/)
      common /aama07/ amzz(3)
!$OMP THREADPRIVATE(/aama07/)
      common /aama10/ ango, aio, aiosq
!$OMP THREADPRIVATE(/aama10/)
      common /aama12/ pop(6)
!$OMP THREADPRIVATE(/aama12/)

cABE 2018/10/30, add threadprivate aamac
      common /aamac03/ pi,coeff,picoef,amassm,amasse,alfa,ne(3),jk(3)
!$OMP THREADPRIVATE(/aamac03/)
      common /aamac11/ coeoct(45),expoct(45),coeo(11),ifo(10),jo(15)
!$OMP THREADPRIVATE(/aamac11/)

*-----------------------------------------------------------------------

      dimension a(6)

*-----------------------------------------------------------------------

      if(n.gt.0) goto 1000

      aio = matelu(n1,l1,n2,l2,3)
      aiosq = aio*aio
      l = (l1+l2+1)/2 
      ango = 2.5d0
      if(iabs(l1-l2).eq.1) ango = 1.5d0
      ango = ango*dble((2*l2+1)*(l-1)*l*(l+1))/
     &       dble((2*l-3)*(2*l-1)*(2*l+1)*(2*l+3))
      roctu = coeo(1)*coeff/(z*amassm)**6.0d0*
     &       (enem/alfa)**7.0d0*ango*aiosq

      return

*-----------------------------------------------------------------------

 1000 etpy = dexp(2.0d0*pi*y)

      extpy = etpy/(etpy-1.0d0)
      yy = y*y
      p = dexp(y*(2.0d0*datan(y/dble(n))-pi))

      do 1100 i = 1, 6
       a(i) = 0.0d0
 1100 continue

      goto (2000,3000,4000), n

*-----------------------------------------------------------------------

 2000 yf = (1.0d0+yy)*(3.0d0+2.0d0*yy)**2.0d0/(yy*(4.0d0+yy)*(9.0d0+yy))
      pf = (15.0d0*(1.0d0+yy)/(3.0d0+2.0d0*yy)*p-1.0d0)**2.0d0
      roctu = coeo(2)*pop(1)*picoef*ango*pf*extpy*yf*
     &        aiosq*amzz(1)**6.0d0

      return

*-----------------------------------------------------------------------

 3000 if(m.eq.2) goto 3100

      yf = (4.0d0+yy)*(6.0d0+yy)**2.0d0/(yy*(1.0d0+yy)*(9.0d0+yy))
      pf = (15.0d0*(2.0d0+3.0d0*yy)/(6.0d0+yy)*p-1.0d0)**2.0d0
      a(1) = coeo(3)*pop(2)*yf*pf

 3100 if(m.eq.1) goto 3200

      yf = (4.0d0+yy)/(1.0d0+yy)
      pf = (3.0d0*p+1.0d0)**2.0d0
      a(2) = coeo(4)*pop(3)*yf*pf
      yf = (4.0d0+yy)*(68.0d0+13.0d0*yy)**2.0d0/((1.0d0+yy)*
     &     (9.0d0+yy)*(16.0d0+yy))
      pf = (5.0d0*(116.0d0+149.0d0*yy)/
     &     (68.0d0+13.0d0*yy)*p-1.0d0)**2.0d0
      a(3) = coeo(5)*pop(3)*yf*pf 

 3200 at = a(1) + a(2) + a(3)
      roctu = picoef*ango*at*extpy*aiosq*amzz(2)**6.0d0

      return

*-----------------------------------------------------------------------

 4000 if(m.gt.1) goto 4100

      yf = (9.0d0+yy)*(27.0d0+2.0d0*yy)**2.0d0/
     &     (yy*(1.0d0+yy)*(4.0d0+yy))
      pf = (2.5d0*(405.0d0+900.0d0*yy+254.0d0*yy*yy)/
     &     ((9.0d0+yy)*(27.0d0+2.0d0*yy))*p-1.0d0)**2.0d0
      a(1) = coeo(6)*pop(4)*yf*pf

 4100 if(m.ne.0 .and. m.ne.2) goto 4200

      yf = (9.0d0+2.0d0*yy)**2.0d0/((1.0d0+yy)*(4.0d0+yy))
      pf = ((27.0d0+13.0d0*yy)/(9.0d0+2.0d0*yy)*p+1.0d0)**2.0d0
      a(2) = coeo(7)*pop(5)*yf*pf
      yf = (9.0d0+yy)*(153.0d0+13.0d0*yy)**2.0d0/
     &     ((1.0d0+yy)*(4.0d0+yy)*(16.0d0+yy))
      pf = (5.0d0*(2349.0d0+3744.0d0*yy+947.0d0*yy*yy)/
     &     ((9.d00+yy)*(153.0d0+13.0d0*yy))*p-1.0d0)**2.0d0
      a(3) = coeo(8)*pop(5)*yf*pf

 4200 if(m.ne.0 .and. m.ne.3) goto 4300

      yf = yy/(1.0d0+yy)
      pf = (2.0*p+1.0d0)**2.0d0
      a(4) = coeo(9)*pop(6)*yf*pf
      yf = yy*(9.0d0+yy)/((1.0d0+yy)*(4.0d0+yy))
      pf = ((63.0d0+47.0d0*yy)/(9.0d0+yy)*p+1.0d0)**2.0d0
      a(5) =coeo(10)*pop(6)*yf*pf
      yf = yy*(9.0d0+yy)*(11.0d0+yy)**2.0d0/
     &     ((1.0d0+yy)*(4.0d0+yy)*(16.0d0+yy)*(25.0d0+yy))
      pf = ((1251.0d0+1850.0d0*yy+439.0d0*yy*yy)/
     &     ((9.0d0+yy)*(11.0d0+yy))*p-1.0d0)**2.0d0
      a(6) = coeo(11)*pop(6)*yf*pf

 4300 at = a(1) + a(2) + a(3) + a(4) + a(5) + a(6)

      roctu = picoef*ango*at*extpy*aiosq*amzz(3)**6.0d0

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      double precision function matel(n1,l1,n2,l2,l,a,n)
*                                                                      *
c  ***  general dimensionless muonic matrix element for penetration
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

cABE 2018/12/05, add threadprivate aamac
      common /aama07/ amzz(3)
!$OMP THREADPRIVATE(/aama07/)
      common /aamac02/ ffact(60),ffactd
!$OMP THREADPRIVATE(/aamac02/)

*-----------------------------------------------------------------------

      an = dble(n1+n2)
      p = 1.0d0 + 1.0d-3*a*dble(n1*n2)/an
      a1 = -2.0d0*dble(n2)/an/p
      a2 = -2.0d0*dble(n1)/an/p

      m1 = n1 - l1
      m2 = n2 - l2
      m3 = n1 + l1 + 1
      m4 = n2 + l2 + 1
      m5 = 2*l1 + 2
      m6 = 2*l2 + 2
      mm = l1 + l2 + l + 3

      s1 = 0.0d0

      do 200 i1 = 1, m1
       k1 = i1 - 1
       s2 = 0.0d0
       do 100 i2 = 1, m2
        k2 = i2 - 1
        la = mm + k1 + k2
        lb = m2 - k2
        lc = m6 + k2
        s2 = s2 + a2**dble(k2)*ffact(la)/(ffact(lb)*ffact(lc)*ffact(i2))
  100  continue
       ld = m1 - k1
       le = m5 + k1
       s1 = s1 + s2*a1**dble(k1)/(ffact(ld)*ffact(le)*ffact(i1))
  200 continue

      aq = 2.0d0/(an*p)
      t1 = s1*dsqrt(ffact(m1)*ffact(m2)/
     &     aq**dble(l-1)*ffact(m3)*ffact(m4))

      matel = t1*(aq*dble(n1))**dble(l+l2+1)*
     &        (aq*dble(n2))**dble(l+l1+1)*
     &        0.5d0**dble(l+1)*ffactd**dble(l+1)/dsqrt(aq**dble(l-1))*
     &        (amzz(n))**dble(l)

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      double precision function matelu(n1,l1,n2,l2,l)
*                                                                      *
c  ***  general dimensionless muonic matrix element for nonpenetration
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

cABE 2018/12/05, add threadprivate aamac
      common /aamac02/ ffact(60),ffactd
!$OMP THREADPRIVATE(/aamac02/)

*-----------------------------------------------------------------------

      an = dble(n1+n2)
      a1 = -2.0d0*dble(n2)/an
      a2 = -2.0d0*dble(n1)/an
      m1 = n1 - l1
      m2 = n2 - l2
      m3 = n1 + l1 + 1
      m4 = n2 + l2 + 1
      m5 = 2*l1 + 2
      m6 = 2*l2 + 2
      mm = l1 + l2 + l + 3
      s1 = 0.0d0

      do 200 i1 = 1, m1
       k1 = i1 - 1
       s2 = 0.0d0
       do 100 i2 = 1, m2
        k2 = i2 - 1
        la = mm + k1 + k2
        lb = m2 - k2
        lc = m6 + k2
        s2 = s2 + a2**dble(k2)*ffact(la)/(ffact(lb)*ffact(lc)*ffact(i2))
  100  continue
       ld = m1 - k1
       le = m5 + k1
       s1 = s1 + s2*a1**dble(k1)/(ffact(ld)*ffact(le)*ffact(i1))
  200 continue

      aq = 2.0d0/an
      t1 = s1*dsqrt(ffact(m1)*ffact(m2)/
     &     aq**dble(l-1)*ffact(m3)*ffact(m4))

      matelu = t1*(aq*dble(n1))**dble(l+l2+1)*
     &         (aq*dble(n2))**dble(l+l1+1)*
     &         0.5d0**dble(l+1)*ffactd**dble(l+1)/dsqrt(aq**dble(l-1))

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      double precision function aamabeta(l1,j1,l2,j2,l)
*                                                                      *
c  ***  depolarization factor
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      a1 = 0.25d0*dble(j1*(j1+2))
      a2 = 0.25d0*dble(j2*(j2+2))

      aamabeta = (a2-dble(l2*(l2+1))+0.75d0)/(a1-dble(l1*(l1+1))+0.75d0)
     &           * (a1+a2-dble(l*(l+1)))/(2.0d0*a2)

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      double precision function point(n,j)
*                                                                      *
c  ***  point-like dirac energy function
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      common /aama14/ dza, dza2, dredm
!$OMP THREADPRIVATE(/aama14/)

*-----------------------------------------------------------------------

      dn = dble(n)
      dj = dble(j)
      d = dza/(dn-dj+dsqrt(dj*dj-dza2))
      d = dredm/dsqrt(1.d0+d*d)
      point = dredm - d

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      function id(a)
*                                                                      *
c  ***  used in decyphering to convert characters into numbers          
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      integer a, b
      dimension b(10)
      data b /1h0,1h1,1h2,1h3,1h4,1h5,1h6,1h7,1h8,1h9/

*-----------------------------------------------------------------------

      id = 0
      do 100 j = 1, 10
       if(a.eq.b(j)) id = j - 1
  100 continue

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      double precision function rav(n,l,m)
*                                                                      *
c  ***used by check to find the exact expectation values
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      goto (100,200,300),m

  100 rav = 0.5d0*dble(3*n*n-l*(l+1))
      return

  200 rav = 0.5d0*dble(n*n*(5*n*n+1-3*l*(l+1)))
      return

  300 rav = 0.125d0*dble(n*n*(35*n*n*(n*n-1)
     &                   - 30*n*n*(l+2)*(l-1)
     &                   + 3*(l+2)*(l+1)*l*(l-1)))

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      block data aama
*                                                                      *
c  ***  block data with all intrinsic parameters and default values     
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

cABE 2018/10/30, add threadprivate aamac
      common /aamac01/ iread,iw,ipunch,iprint
!$OMP THREADPRIVATE(/aamac01/)
      common /aamac02/ ffact(60),ffactd
!$OMP THREADPRIVATE(/aamac02/)
      common /aamac03/ pi,coeff,picoef,amassm,amasse,alfa,ne(3),jk(3)
!$OMP THREADPRIVATE(/aamac03/)
      common /aamac04/ amassa,amassn
!$OMP THREADPRIVATE(/aamac04/)
      common /aamac05/ m1(7),m2(7),m3(7),idb
!$OMP THREADPRIVATE(/aamac05/)
      common /aamac06/ zmk,zml,zmm,zmkm,zmlm,zmmm,cl1,cl2,ivers
!$OMP THREADPRIVATE(/aamac06/)
      common /aamac07/ jtm(6),jtd(6),jtq(6),jto(6)
!$OMP THREADPRIVATE(/aamac07/)
      common /aamac08/ coemon(30),expmon(30),ifm(6),jm(10)
!$OMP THREADPRIVATE(/aamac08/)
      common /aamac09/ coedip(42),expdip(42),coed(4),coedp(9),
     &                 ifd(9),jd(14)
!$OMP THREADPRIVATE(/aamac09/)
      common /aamac10/ coequa(45),expqua(45),coeq(11),ifq(10),jq(15)
!$OMP THREADPRIVATE(/aamac10/)
      common /aamac11/ coeoct(45),expoct(45),coeo(11),ifo(10),jo(15)
!$OMP THREADPRIVATE(/aamac11/)
      common /aamac12/ hbar,widthk,npol(20),ipol,ide,ip8,ipc(3),nmax
!$OMP THREADPRIVATE(/aamac12/)
      common /aamac13/ ehigh,elow,climit,eres,esp,espm,cd(5),ea,eb,
     &                 icc,mpu,icpu(200),ipn
!$OMP THREADPRIVATE(/aamac13/)
      common /aamac14/ tfm,step,rmatch,alexp,nopt,idir,iyc
!$OMP THREADPRIVATE(/aamac14/)
      common /aamac15/ ij(4),yj(4),jj1(4)
!$OMP THREADPRIVATE(/aamac15/)
      common /aamac16/ ic,ll(20)
!$OMP THREADPRIVATE(/aamac16/)
      common /aamac17/ elecbe(104,3)
      common /aamadbg/ idebug

*-----------------------------------------------------------------------
*     constant value
*-----------------------------------------------------------------------

        data pi /3.1415926535d0/
        data picoef /1.298778d17/
        data coeff /4.134139d16/
        data alfa /7.297353d-3/
        data amasse,amassm /5.110034d5,206.7686d0/
        data amassn,hbar /931.48d0,6.582173d-16/

*-----------------------------------------------------------------------
c     default parameter in input card
*-----------------------------------------------------------------------

        data esp /0.d0/
        data widthk /3.d0/
        data ipc /0,1,1/
        data nopt,nmax,alexp /0,16,0.d0/
        data ip8 /0/
        data npol /20*-1/
        data ipol /0/
        data ehigh,elow,climit,eres /1.d1,1.d-3,1.d-2,1.d-3/
        data tfm /2.3001d0/
        data step,rmatch /0.d0,0.d0/
        data ea,eb /99.d0,99.d0/
        data ide /10000/
        data ipn /1/
        data ic,iread,iw,ipunch,iprint /0,5,6,6,6/
        data cd /1.d-1,1.d-2,1.d-3,1.d-4,1.d-5/
        data idb /0/
        data ffactd /90.0d0/

        data m1 /0,0,0,0,1,2,3/
        data m2 /0,0,0,0,1,2,3/
        data m3 /0,0,0,0,1,2,3/
        data jtm,jtd,jtq,jto /6*1,6*1,6*1,6*1/

*-----------------------------------------------------------------------
*     default parameter in input card but no description in manual
*-----------------------------------------------------------------------

        data amassa /0.d0/
        data ifm /6*0/
        data ifd /9*0/
        data ifq /10*0/
        data ifo /10*0/
        data jm,jd,jq,jo /10*1,14*1,15*1,15*1/
        data zmk,zml,zmm /2.d0,4.d0,9.d0/
        data zmkm,zmlm,zmmm /4.d0,8.d0,18.d0/
        data iyc /0/
        data ij /1,1,1,1/
        data yj /0.d0,0.d0,0.d0,0.d0/
        data jj1 /1,1,1,1/
        data mpu /0/

*-----------------------------------------------------------------------
*     default parameter
*-----------------------------------------------------------------------

        data ne,jk /0,2,2,0,2,3/
        data cl1,cl2 /0.d0,0.d0/

        data coedp /
     &    2.133333d+01,  4.266667d+01,  3.555556d+00,  1.137778d+02,
     &    7.111111d+00,  5.688889d+01,  1.024000d+03,  2.275556d+01,
     &    1.228800d+03/
        data coeq /
     &    6.666667d-02,  8.888889d-02,  6.944444d-04,  1.000000d-02,
     &    9.375000d-04,  4.064421d-05,  3.511660d-03,  6.503074d-05,
     &    2.107000d-02,  9.290105d-05,  2.064468d-06/
        data coeo /
     &    1.693122d-03,  1.015873d-02,  1.984127d-05,  2.125850d-04,
     &    3.985969d-07,  5.734633d-08,  1.474620d-05,  5.461556d-09,
     &    2.654316d-05,  7.646177d-07,  1.365389d-06/
        data coed /
     &    1.333333d+00,  2.133333d+01,  1.066667d+01,  7.111111d+00/

        data idir /0/
        data ivers /0/
        data ll /0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19/

        data expmon /
     &    4.531799,  4.706453,  4.736786,  3.471234,  3.150354,
     &    3.052991,  4.703517,  3.459110,  3.249048,  3.150354,
     &    2.977822,  2.917974,  3.312419,  2.695684,  2.533755,
     &    4.914818,  2.970119,  2.700216,  6.548446,  3.199174,
     &    2.830713,  2.695684,  2.415030,  2.323143,  2.974074,
     &    2.526711,  2.402950,  2.415030,  2.249824,  2.188693/
        data expdip /
     &    3.766185,  4.397381,  4.511374,  2.429543,  2.723240,
     &    2.747787,  2.976262,  2.927947,  2.886848,  4.076930,
     &    3.354237,  3.192531,  2.723240,  2.747093,  2.736171,
     &    2.064114,  2.209579,  2.187568,  2.546045,  2.374398,
     &    2.298389,  2.876576,  2.498902,  2.385350,  3.993301,
     &    2.876576,  2.652421,  5.778643,  3.134388,  2.797892,
     &    2.209579,  2.154546,  2.118914,  2.374398,  2.234588,
     &    2.178575,  2.876576,  2.498902,  2.385350,  2.154546,
     &    2.084353,  2.051376/
        data expqua /
     &    3.531162,  4.240205,  4.388131,  2.543817,  2.683365,
     &    2.694792,  2.131239,  2.520683,  2.579500,  2.566813,
     &    2.816281,  2.819550,  2.520683,  2.607425,  2.618563,
     &    1.726067,  1.977932,  1.998559,  2.061514,  2.104925,
     &    2.087356,  2.283507,  2.199915,  2.155173,  2.183688,
     &    2.283507,  2.245050,  2.663082,  2.438142,  2.347529,
     &    1.977932,  1.998217,  1.987446,  2.104925,  2.062287,
     &    2.036313,  4.340990,  2.986053,  2.726274,  2.283507,
     &    2.199915,  2.155173,  1.998217,  1.972396,  1.956081/
        data expoct /
     &    3.429305,  4.153371,  4.322206,  1.982441,  2.395708,
     &    2.470619,  2.341001,  2.538774,  2.572215,  2.195239,
     &    2.571516,  2.623072,  2.395708,  2.512490,  2.537093,
     &    1.558161,  1.837921,  1.876801,  1.841217,  1.946305,
     &    1.953014,  2.017891,  2.026300,  2.011369,  1.778407,
     &    2.017891,  2.032372,  2.117493,  2.140989,  2.117166,
     &    1.837921,  1.892340,  1.895838,  1.946305,  1.948185,
     &    1.938261,  2.252893,  2.331797,  2.283271,  2.017891,
     &    2.026300,  2.011369,  1.892340,  1.921936,  1.898511/
        data coemon /
     &    9.278462d-01, -4.686876d-02,  2.606690d-03,  3.258145d-01,
     &   -1.650726d-02,  9.192955d-04, -8.728775d-02,  5.557533d-03,
     &   -3.299306d-04,  1.650726d-02, -7.889860d-04,  1.644805d-05,
     &    1.837462d-01, -9.081039d-03,  5.039261d-04, -7.674657d-02,
     &    4.090333d-03, -2.414850d-04,  6.976183d-03, -3.308179d-04,
     &    2.022719d-05,  9.886183d-03, -4.700158d-04,  9.786634d-06,
     &   -1.115213d-03,  5.895603d-05, -1.271168d-06,  4.522730d-05,
     &   -1.882362d-06,  2.089122d-08/
        data coedip /
     &    9.864278d-02, -3.552196d-03,  6.908681d-05,  2.415974d-02,
     &   -8.812569d-04,  1.721286d-05, -6.819236d-03,  3.101739d-04,
     &   -6.390954d-06,  7.532858d-02, -4.445121d-03,  2.592835d-04,
     &    2.203142d-04, -8.210685d-06,  8.801341d-08,  3.235460d-02,
     &   -1.181545d-03,  2.303800d-05, -1.225112d-02,  5.542601d-04,
     &   -1.140757d-05,  8.892663d-04, -4.580119d-05,  9.800079d-07,
     &    1.211143d-02, -6.669497d-04,  3.872454d-05, -1.512719d-03,
     &    7.933719d-05, -4.808937d-06,  4.376093d-05, -1.627056d-06,
     &    1.743224d-08, -5.132038d-06,  2.095740d-07, -2.313305d-09,
     &    2.223166d-04, -1.145030d-05,  2.450020d-07,  9.039198d-08,
     &   -3.079553d-09,  2.088144d-11/
        data coequa /
     &    1.411675d-01, -3.944496d-03,  3.945040d-05,  3.270497d-01,
     &   -1.133197d-02,  1.194093d-04, -2.776085d-01,  7.846003d-03,
     &   -7.866167d-05,  4.601989d-01, -1.831626d-02,  3.675330d-04,
     &    3.487113d-03, -1.060614d-04,  6.944398d-07,  1.253120d-00,
     &   -3.541941d-02,  3.549579d-04, -4.926037d-01,  1.705249d-02,
     &   -1.796112d-04,  3.681114d-02, -1.438755d-03,  1.569441d-05,
     &    4.626170d-01, -1.840557d-02,  3.690030d-04, -4.608811d-02,
     &    2.220260d-03, -4.671028d-05,  7.870980d-03, -2.392062d-04,
     &    1.565542d-06, -9.473607d-04,  3.134314d-05, -2.107953d-07,
     &    4.492433d-02, -2.608238d-03,  1.552426d-04,  9.202786d-03,
     &   -3.596887d-04,  3.923602d-06,  5.980154d-05, -1.722622d-06,
     &    7.889632d-09/
        data coeoct /
     &    1.835095d-02, -4.182848d-04,  2.558331d-06,  1.444508d-01,
     &   -3.328672d-03,  2.037957d-05, -4.367688d-02,  1.223204d-03,
     &   -7.847840d-06,  3.021241d-01, -8.987655d-03,  9.179069d-05,
     &    4.438300d-03, -1.140520d-04,  5.042711d-07,  1.467168d-00,
     &   -3.379181d-02,  2.068045d-04, -5.924030d-01,  1.655453d-02,
     &   -1.061699d-04,  4.508567d-02, -1.416945d-03,  9.394653d-06,
     &    6.820167d-01, -2.028855d-02,  2.071180d-04, -6.910950d-02,
     &    2.487891d-03, -2.660315d-05,  2.252870d-02, -5.785850d-04,
     &    2.556806d-06, -2.759089d-03,  7.677244d-05, -3.479514d-07,
     &    3.778402d-01, -1.582449d-02,  3.229544d-04,  4.508567d-03,
     &   -1.416945d-04,  9.394653d-07,  1.285744d-06, -3.257990d-08,
     &    1.067759d-10/

*-----------------------------------------------------------------------
      end

************************************************************************
*                                                                      *
      block data elecbind
*                                                                      *
*     average value of atomic electron binding energy [eV]             *
*     for k, l, and m shell                                            *
*     taken from below reference                                       *
*                                                                      *
*     [1] F.B.Larkins, At. Data and Nucl.Data Tables 20, 313 (1977)    *
*     [2] K.D.Sevier, Low energy enectron spectrometry (1972)          *
*     [3] F.T.Porter+, J. Phys. Chem. Ref. Data 7, 1267 (1978)         *
*     [4] D.A.Shirley+,Phys. Rev. B15, 544 (1977)                      *
*     [5] J.A.Bearden+, Rev. Mod. Phys., 39, 125 (1967)                *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     elecbe(iz,j)): atomic electron binding energy [eV]               *
*                    iz: atomic number                                 *
*                    j = 1: k-shell, = 2: l-shell, = 3: m-shell        *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      common /aamac17/ elecbe(104,3)
      common /aamadbg/ idebug

*-----------------------------------------------------------------------

      data idebug /0/

      data ( ( elecbe(iz,j), j = 1, 3 ), iz = 1, 104 ) /
*             k-shell       l-shell       m-shell
     &    1.360000d+1,  0.000000d+0,  0.000000d+0,
     &    2.460000d+1,  0.000000d+0,  0.000000d+0,
     &    5.480000d+1,  5.300000d+0,  0.000000d+0,
     &    1.121000d+2,  8.000000d+0,  0.000000d+0,
     &    1.880000d+2,  9.966667d+0,  0.000000d+0,
     &    2.838000d+2,  1.220000d+1,  0.000000d+0,
     &    4.016000d+2,  1.528000d+1,  0.000000d+0,
     &    5.320000d+2,  1.423333d+1,  0.000000d+0,
     &    6.854000d+2,  1.585714d+1,  0.000000d+0,
     1    8.701000d+2,  2.835000d+1,  0.000000d+0,
     &    1.072100d+3,  3.915000d+1,  7.000000d-1,
     &    1.305000d+3,  6.090000d+1,  2.100000d+0,
     &    1.559600d+3,  8.407500d+1,  2.300000d+0,
     &    1.838900d+3,  1.115000d+2,  5.300000d+0,
     &    2.145500d+3,  1.490250d+2,  1.242000d+1,
     &    2.472000d+3,  1.807500d+2,  1.060000d+1,
     &    2.822400d+3,  2.179500d+2,  9.857143d+0,
     &    3.206000d+3,  2.685500d+2,  1.917500d+1,
     &    3.607400d+3,  3.151500d+2,  1.940000d+1,
     2    4.038100d+3,  3.701500d+2,  2.398000d+1,
     &    4.492800d+3,  4.278750d+2,  2.920000d+1,
     &    4.966400d+3,  4.840500d+2,  2.666923d+1,
     &    5.465100d+3,  5.436250d+2,  2.852308d+1,
     &    5.989200d+3,  6.068250d+2,  2.978571d+1,
     &    6.539000d+3,  6.752500d+2,  3.216667d+1,
     &    7.112000d+3,  7.458500d+2,  3.366250d+1,
     &    7.708900d+3,  8.191000d+2,  3.263333d+1,
     &    8.332800d+3,  8.973500d+2,  3.712222d+1,
     &    8.978900d+3,  9.773250d+2,  3.873333d+1,
     3    9.658600d+3,  1.068950d+3,  4.846667d+1,
     &    1.036710d+4,  1.167700d+3,  6.196667d+1,
     &    1.110310d+4,  1.273875d+3,  7.700000d+1,
     &    1.186670d+4,  1.382825d+3,  9.298889d+1,
     &    1.265780d+4,  1.500425d+3,  1.118889d+2,
     &    1.347370d+4,  1.619450d+3,  1.284444d+2,
     &    1.432560d+4,  1.749500d+3,  1.571444d+2,
     &    1.519970d+4,  1.884450d+3,  1.778889d+2,
     &    1.610460d+4,  2.025575d+3,  2.049778d+2,
     &    1.703840d+4,  2.172000d+3,  2.331111d+2,
     4    1.799760d+4,  2.320725d+3,  2.600333d+2,
     &    1.898560d+4,  2.475850d+3,  2.890444d+2,
     &    1.999950d+4,  2.632750d+3,  3.156111d+2,
     &    2.104400d+4,  2.797375d+3,  3.456000d+2,
     &    2.211720d+4,  2.966675d+3,  3.771556d+2,
     &    2.321990d+4,  3.141400d+3,  4.094333d+2,
     &    2.435030d+4,  3.320300d+3,  4.417889d+2,
     &    2.551400d+4,  3.507925d+3,  4.787111d+2,
     &    2.671120d+4,  3.705000d+3,  5.206667d+2,
     &    2.793990d+4,  3.908925d+3,  5.652556d+2,
     5    2.920010d+4,  4.119600d+3,  6.122222d+2,
     &    3.049120d+4,  4.335775d+3,  6.603444d+2,
     &    3.181380d+4,  4.558500d+3,  7.104889d+2,
     &    3.316940d+4,  4.788600d+3,  7.636222d+2,
     &    3.456440d+4,  5.030225d+3,  8.267667d+2,
     &    3.598460d+4,  5.274375d+3,  8.814222d+2,
     &    3.744060d+4,  5.526600d+3,  9.431333d+2,
     &    3.892460d+4,  5.780575d+3,  1.000511d+3,
     &    4.044300d+4,  6.039950d+3,  1.058967d+3,
     &    4.199060d+4,  6.300950d+3,  1.114222d+3,
     6    4.356890d+4,  6.565825d+3,  1.167311d+3,
     &    4.518400d+4,  6.839825d+3,  1.224322d+3,
     &    4.683420d+4,  7.120250d+3,  1.283967d+3,
     &    4.851900d+4,  7.405725d+3,  1.343222d+3,
     &    5.023910d+4,  7.697875d+3,  1.405233d+3,
     &    5.199570d+4,  7.996900d+3,  1.470156d+3,
     &    5.378850d+4,  8.301650d+3,  1.532167d+3,
     &    5.561770d+4,  8.613550d+3,  1.596744d+3,
     &    5.748550d+4,  8.932850d+3,  1.663378d+3,
     &    5.938960d+4,  9.257150d+3,  1.733100d+3,
     7    6.133230d+4,  9.587950d+3,  1.800744d+3,
     &    6.331380d+4,  9.926800d+3,  1.871800d+3,
     &    6.535080d+4,  1.028288d+4,  1.955489d+3,
     &    6.741640d+4,  1.064495d+4,  2.039600d+3,
     &    6.952500d+4,  1.101435d+4,  2.125256d+3,
     &    7.167640d+4,  1.138900d+4,  2.210489d+3,
     &    7.387080d+4,  1.177370d+4,  2.299667d+3,
     &    7.611100d+4,  1.216825d+4,  2.393022d+3,
     &    7.839480d+4,  1.257018d+4,  2.486833d+3,
     &    8.072490d+4,  1.298095d+4,  2.584222d+3,
     8    8.310230d+4,  1.340395d+4,  2.687644d+3,
     &    8.553040d+4,  1.383990d+4,  2.796789d+3,
     &    8.800450d+4,  1.428280d+4,  2.906767d+3,
     &    9.052590d+4,  1.473395d+4,  3.018133d+3,
     &    9.310000d+4,  1.519625d+4,  3.134444d+3,
     &    9.572400d+4,  1.566800d+4,  3.251889d+3,
     &    9.839700d+4,  1.614900d+4,  3.373000d+3,
     &    1.011300d+5,  1.664575d+4,  3.501222d+3,
     &    1.039150d+5,  1.715100d+4,  3.631778d+3,
     &    1.067560d+5,  1.766675d+4,  3.764556d+3,
     9    1.096500d+5,  1.819125d+4,  3.898111d+3,
     &    1.125960d+5,  1.872125d+4,  4.027556d+3,
     &    1.156020d+5,  1.926050d+4,  4.159889d+3,
     &    1.186690d+5,  1.981175d+4,  4.296111d+3,
     &    1.217910d+5,  2.037100d+4,  4.430111d+3,
     &    1.249820d+5,  2.094500d+4,  4.570000d+3,
     &    1.282410d+5,  2.152925d+4,  4.713889d+3,
     &    1.315560d+5,  2.212425d+4,  4.858000d+3,
     &    1.349390d+5,  2.273300d+4,  5.005222d+3,
     &    1.383960d+5,  2.335375d+4,  5.155000d+3,
     1    1.419260d+5,  2.398775d+4,  5.307889d+3,
     &    1.465260d+5,  2.463425d+4,  5.463111d+3,
     &    1.492080d+5,  2.529450d+4,  5.621333d+3,
     &    1.529700d+5,  2.597600d+4,  5.789333d+3,
     &    1.562880d+5,  2.667025d+4,  5.959000d+3/

*-----------------------------------------------------------------------
      end
