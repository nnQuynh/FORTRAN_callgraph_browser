************************************************************************
*                                                                      *
      subroutine ovly14
*                                                                      *
*        modified by K.Niita on 2007/03/05                             *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              calculate charged particle spectrum by nuclear data     *
*                                                                      *
************************************************************************
C for USE_MOD_COUNTER
      use mod_counter, only: rncnt,rnint,rnintr,rnpnt,rnpntr
C for  REDUCTION_COUNTER
!$   &                      ,rncnt2,rnint2,rnintr2,rnpnt2,rnpntr2
      use moddas_material

*-----------------------------------------------------------------------

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param00.inc'
      include 'param02.inc'
      include 'param.inc'

*-----------------------------------------------------------------------

      parameter ( pir = 3.1415926535897932d0 )

*-----------------------------------------------------------------------

      common /inout/  ins, ios

*-----------------------------------------------------------------------

      common /mpi00/ npe, me

*-----------------------------------------------------------------------


*-----------------------------------------------------------------------

      common /clustf/ nclst, iclust(nnn)
!$OMP THREADPRIVATE(/clustf/)
      common /clustg/ jclust(0:8,nnn), qclust(0:12,nnn)
!$OMP THREADPRIVATE(/clustg/)
      common /clustp/ rumpat(0:20), numpat(0:20)
!$OMP THREADPRIVATE(/clustp/)

      common /clustt/ nclsts, iclusts(nnn)
!$OMP THREADPRIVATE(/clustt/)
      common /clustw/ jclusts(0:8,nnn),  qclusts(0:12,nnn)
!$OMP THREADPRIVATE(/clustw/)

      common /clustv/ kdecay(4)
!$OMP THREADPRIVATE(/clustv/)

      common /swich3/ ielst, jelst, kelst
!$OMP THREADPRIVATE(/swich3/)
      common /engch/  ejamnu, ejampi, eisobar, eqmdnu, eqmdmn, ejamqmd

      common /paraj/  mstz(300), parz(300)

      common /paran/  icfn(100), ilfn(100), chfn(100)
      character chfn*200

*-----------------------------------------------------------------------

      common /kmat1g/ kmat(kvlmax)
      common /xgeosm/ ksig(kvlmax)

      common /ndemax/ dnmax(20)
      common /regdc/  idrg(kvlmax), idgr(kvmmax)
      common /regcm/  icmg(kvlmax)
      common /celdn/  denr(kvlmax), denm(kvlmax), denc(kvlmax)


*-----------------------------------------------------------------------


      integer nrandgen
      common /randn/ nrandgen ! S.H. xorshift (2020.5.29)
      integer*8 :: iranji64 ! S.H. xorshift (2020.5.29)
      common /randm4/ rnfb,rnfs,rngb,rngs,rnmult,ranj,rani,
     &                rnrtc,nstrid,inif, iranji64
      integer*8 :: iransb64 ! S.H. xorshift (2020.2.6)
      common /randtp/ rijk,rans,ranb, iransb64
!$OMP THREADPRIVATE(/randtp/)

*-----------------------------------------------------------------------

      common /isomul/ smlwt(isrc), totfact, imsrc
      common /isorst/ jstyp(isrc), istyp(isrc), inkf0(isrc), lstyp(isrc)
!$OMP THREADPRIVATE(/isorst/)
      common /isorsp/ sx0(isrc), sy0(isrc), sz0(isrc), sx1(isrc),
     &                sy1(isrc), sz1(isrc), sr0(isrc), se0(isrc),
     &                sdir(isrc), srx(isrc), sry(isrc), swem(isrc),
     &                sphi(isrc), sdom(isrc), swt0(isrc)

*-----------------------------------------------------------------------

      common /eparm/  esmax, esmin, emin(20)
      common /cparm/  maxbch,maxcas
      common /qparm/  ielas,icasc,iqstep,lvlopt,igamma
      common /bparm/  andt,jevap,npidk
      common /tcntl/  icntl, inucr
      common /spred/ nspred, nwsprd, nedisp, itstep, ndedx

*-----------------------------------------------------------------------

      common /comps1/ mstapr, massta, msprpr, masspr
      common /comps2/ sigela, signon, fissx

*-----------------------------------------------------------------------

      dimension cngd(50,20,104)
      equivalence ( das, cngd )

      dimension izms(50)
      dimension izmr(50)
      dimension heat(50)
      dimension xtot(50)
      dimension xmtt(50)
      dimension angd(19)
      dimension angm(19,2)
      dimension engm(104)
      dimension elet(43)
      dimension emet(43)

*-----------------------------------------------------------------------

      data engm / 2.1197E+1,1.9180E+1,1.7355E+1,1.5703E+1,1.4209E+1,
     &            1.2857E+1,1.1633E+1,1.0526E+1,9.5242E+0,8.6179E+0,
     &            7.7978E+0,7.0557E+0,6.3843E+0,5.7767E+0,5.2270E+0,
     &            4.7296E+0,4.2795E+0,3.8723E+0,3.5038E+0,3.1703E+0,
     &            2.8686E+0,2.5957E+0,2.3487E+0,2.1252E+0,1.9229E+0,
     &            1.7399E+0,1.5744E+0,1.4246E+0,1.2890E+0,1.1663E+0,
     &            1.0553E+0,9.5489E-1,8.6402E-1,7.8180E-1,7.0740E-1,
     &            6.4008E-1,5.7917E-1,5.2405E-1,4.7418E-1,4.2906E-1,
     &            3.8823E-1,3.5128E-1,3.1785E-1,2.8761E-1,2.6024E-1,
     &            2.3548E-1,2.1307E-1,1.9279E-1,1.7395E-1,1.5735E-1,
     &            1.4283E-1,1.2923E-1,1.1693E-1,9.8804E-2,7.6949E-2,
     &            5.9928E-2,4.6672E-2,3.6348E-2,2.8308E-2,2.2047E-2,
     &            1.7170E-2,1.3372E-2,1.0414E-2,8.1103E-3,6.3164E-3,
     &            4.9192E-3,3.8310E-3,2.9836E-3,2.3237E-3,1.8097E-3,
     &            1.4094E-3,1.0976E-3,8.5482E-4,6.6574E-4,5.1848E-4,
     &            4.0379E-4,3.1448E-4,2.4491E-4,1.9074E-4,1.4855E-4,
     &            1.1569E-4,9.0097E-5,7.0168E-5,5.4647E-5,4.2559E-5,
     &            3.3145E-5,2.5813E-5,2.0104E-5,1.5657E-5,1.2194E-5,
     &            9.4962E-6,7.3957E-6,5.7598E-6,4.4857E-6,3.4935E-6,
     &            2.7207E-6,2.1184E-6,1.6497E-6,1.2852E-6,1.0009E-6,
     &            7.7951E-7,6.0708E-7,4.7280E-7,2.0705E-7/

      data angd /   0.0,  10.0,  20.0,  30.0,  40.0,  50.0,  60.0,
     &             70.0,  80.0,  90.0, 100.0, 110.0, 120.0, 130.0,
     &            140.0, 150.0, 160.0, 170.0, 180.0 /

      data elet / 1.00E-01,1.50E-01,2.00E-01,3.00E-01,4.00E-01,
     &            6.00E-01,8.00E-01,1.00E+00,1.50E+00,2.00E+00,
     &            3.00E+00,4.00E+00,6.00E+00,8.00E+00,1.00E+01,
     &            1.50E+01,2.00E+01,3.00E+01,4.00E+01,6.00E+01,
     &            8.00E+01,1.00E+02,1.50E+02,2.00E+02,3.00E+02,
     &            4.00E+02,6.00E+02,8.00E+02,1.00E+03,1.50E+03,
     &            2.00E+03,3.00E+03,4.00E+03,6.00E+03,8.00E+03,
     &            1.00E+04,1.50E+04,2.00E+04,3.00E+04,4.00E+04,
     &            6.00E+04,8.00E+04,1.00E+05/

*-----------------------------------------------------------------------

      character cha*1
      data cha /"'"/

      character yen*1
      yen  = char(92)

*-----------------------------------------------------------------------
*        only for executable PE
*-----------------------------------------------------------------------

            if( me .eq. 0 .and. npe .gt. 1 ) return

*-----------------------------------------------------------------------
*        output unit ( 61 )
*-----------------------------------------------------------------------

            io = 61

            open(io, file = chfn(11), status = 'unknown' )

*-----------------------------------------------------------------------
*        random number generator
*-----------------------------------------------------------------------

               if( inif .eq. 0 ) call advijk

               inif = 0
            if ( nrandgen .eq. 0 ) then ! S.H. xorshift (2020.5.29)
               ranb = rani
               rans = ranj
            else
               iransb64 = iranji64
            end if

*-----------------------------------------------------------------------
*        mevent ; number of events
*-----------------------------------------------------------------------

               mevent = maxcas

*-----------------------------------------------------------------------
*        pick up incident particle and energy from source
*-----------------------------------------------------------------------

               ityp = istyp(imsrc)
               ktyp = inkf0(imsrc)
               jtyp = ichgf(ityp,ktyp)
               mtyp = ibryf(ityp,ktyp)
               rtyp = rmtyp(ityp,ktyp)

               ein  = se0(imsrc)

            if( ityp .ne. 2 ) then

               write(6,*) ' **** Particle should be neutron. *****'
               ierr = 1
               goto 999

            end if

            if( ein .gt. dnmax(2) ) then

               write(6,*) ' **** Energy should be less than '//
     &                    'dnmax(2). *****', ein, dnmax(2)
               ierr = 1
               goto 999

            end if

*-----------------------------------------------------------------------
*        min and max energy for range
*-----------------------------------------------------------------------

               esmax = parz(170)
               esmin = parz(162)

*-----------------------------------------------------------------------
*        pick up the projectile nucleus from the first material
*-----------------------------------------------------------------------

               icl   = 1
               mat   = 1
               ireg  = 1

               lemm  = nint( dnel_das(kmat0+mat) )
               hydro = denh_das(kmat0+mat)
               jimat = lemm

               if( hydro .gt. 0.0 ) jimat = jimat + 1

            if( jimat .ne. 1 ) then

               write(6,*) ' **** Target should be only one. *****'
               ierr = 1
               goto 999

            end if

            if( lemm .eq. 1 ) then

               zpr = zz_das(kmat(mat)+lemm)
               apr = a_das(kmat(mat)+lemm)

            else

               apr = 1.0
               zpr = 1.0

            end if

               nta = nint( apr )
               ntz = nint( zpr )

               ipty = 1
               if( nta .eq. 1 .and. ntz .eq. 1 ) ipty = 2

            if( nta .eq. 1 .and. ntz .eq. 1 ) then

               ipart = 1
               kpart = 2212

            else

               ipart = 19
               kpart = 1000000 * ntz + nta

            end if

               jpart = ntz
               rpart = rmtyp(ipart,kpart)

*-----------------------------------------------------------------------
*        check whether second material is water or not
*-----------------------------------------------------------------------

               mat2 = 2
               lem2 = nint( dnel_das(kmat0+mat2) )
               hyd2 = denh_das(kmat0+mat2)
               iz2 = zz_das(kmat(mat2)+lem2)
               ia2 = a_das(kmat(mat2)+lem2)

            if( lem2 .ne. 1 .or. hyd2 .le. 0.0d0 .or.
     &          iz2 .ne. 8 .or. ia2 .ne. 16 ) then

               write(6,*) ' **** second material should be water *****'
               ierr = 1
               goto 999

            end if

*-----------------------------------------------------------------------
*        set array and data point
*-----------------------------------------------------------------------

               icn   = mmmax
               mmmax = mmmax + 20*18*104

            if( mmmax .gt. mdas ) then

               write(6,*) ' **** memory is over than mdas *****'
               ierr = 1
               goto 999

            end if

*-----------------------------------------------------------------------
*        mesh points
*-----------------------------------------------------------------------

               do i = 1, 52

                  engr = engm(105-i)
                  engm(105-i) = engm(i)
                  engm(i) = engr

               end do

               if( nta .gt. 1 ) then

                  anga = 5.d0

               else if( nta .eq. 1 ) then

                  anga = 2.5d0

                  do i = 1, 19

                     angd(i) = angd(i) / 2.0d0

                  end do

               end if

               do i = 1, 19

                  angm(i,1) = cos( ( angd(i) - anga ) / 180.d0 * pir )
                  angm(i,2) = cos( ( angd(i) + anga ) / 180.d0 * pir )

               end do

                  angm(1,1)  =  1.0d0
                  angm(19,2) = -1.0d0
                  if( nta .eq. 1 ) angm(19,2) = 0.0d0

               do i =1, 50

                  izms(i) = 0
                  izmr(i) = i
                  heat(i) = 0.0d0
                  xtot(i) = 0.0d0
                  xmtt(i) = 0.0d0

               end do

                  iznm = 0

               do i = 1, 43

                  emet(i) = 0.0d0

               end do

               do i = 1, 50
               do j = 1, 20
               do k = 1, 104

                  cngd(icn+i,j,k) = 0.0d0

               end do
               end do
               end do

*-----------------------------------------------------------------------
*        neutron for nuclear data
*-----------------------------------------------------------------------

                  mk = mat
                  rh = denm(mat)
                  tme = 0.0

            call xstneu(0,sigt,sigaa,icl,ein,tme,mk)

                  sigtot = sigt * 1000.d0
                  signon = sigt * 1000.d0
                  sigela = 0.0

                  bmax0  = sqrt( signon / 10.d0 / pir )

            call heatn(icl,ein,heatr,heatf,mk,rh,0,0,mtdum)

                  heat0 = ( heatr + heatf ) / rh / sigt

*-----------------------------------------------------------------------
*     Nuclear Reactions
*-----------------------------------------------------------------------

                  heat1 = 0.0d0

      do 1000 irunp = 1, mevent

*-----------------------------------------------------------------------
*        inelastic nuclear reaction
*-----------------------------------------------------------------------

                  wgti = 1.0d0

                  call sctneut(ein,wgti,ireg,mat,ffac1,ffac2,0,elrt)

*-----------------------------------------------------------------------
*        final summary
*-----------------------------------------------------------------------

         if( nclsts .gt. 0 ) then

            do i = 1, nclsts

                     iz = jclusts(1,i)
                     in = jclusts(2,i)

               if( iz .gt. 0 ) then

                     px = qclusts(1,i)
                     py = qclusts(2,i)
                     pz = qclusts(3,i)
                     ek = qclusts(7,i)

                     ip = jclusts(3,i)
                     ik = jclusts(7,i)
                     rm = rmtyp(ip,ik)

                     heat1 = heat1 + ek

                     pl2 = px**2 + py**2 + pz**2
                     pla = sqrt(pl2)

                     if( pl2 .eq. 0 ) then
                        cosa = 1.0
                     else
                        cosa = pz / sqrt(pl2)
                     end if

                     izin = iz * 1000 + in + iz

*-----------------------------------------------------------------------

                  if( iznm .eq. 0 ) then

                        iznm = 1
                        izhh = 1
                        izms(1) = izin

                  else

                     do j = 1, iznm

                        if( izin .eq. izms(j) ) goto 200

                     end do

                        iznm = iznm + 1

                        if( iznm .gt. 50 ) then
                           write(6,*) ' **** iznm is over than 50 *****'
                           write(6,'(2i8)') (k,izms(k),k=1,50)

                           ierr = 1
                           goto 999

                        end if

                        j = iznm
                        izms(j) = izin

  200                   izhh = j

                  end if

*-----------------------------------------------------------------------
*                 DDX of Charge Particle
*-----------------------------------------------------------------------

                  do j = 1, 19

                     if( cosa .le. angm(j,1) .and.
     &                   cosa .gt. angm(j,2) ) goto 300

                  end do

                        j = 19

  300                   iang = j

                        dang = angm(j,1) - angm(j,2)

*-----------------------------------------------------------------------

                  do j = 1, 104

                     if( ek .le. engm(j) ) goto 400

                  end do

                        j = 104

  400                   ieng = j

                     if( j .eq. 1 ) then

                        deng = engm(j)

                     else

                        deng = engm(j) - engm(j-1)

                     end if

*-----------------------------------------------------------------------

                     heat(izhh) = heat(izhh) + ek
                     xtot(izhh) = xtot(izhh) + 1.d0

                     cngd(icn+izhh,iang,ieng) =
     &               cngd(icn+izhh,iang,ieng) + 1.d0 / dang / deng

                     cngd(icn+izhh,20,ieng) =
     &               cngd(icn+izhh,20,ieng) + 1.d0 / deng

*-----------------------------------------------------------------------
*                 LET Distribution
*-----------------------------------------------------------------------

                     apo = dble( iz + in )
                     zpo = dble( iz )

                     ipty = 1
                     if( iz .eq. 1 .and. in .eq. 0 ) ipty = 2

                     inmr = 200

                     edma = ek
                     edmi = 1.d-8
                     edif = log( edma / edmi ) / dble( inmr )

                     ek1  = ek

                     if( ndedx .eq. 0 .or. ndedx .eq. 2 ) then

                        call spar(ipty,apo,zpo,ek1,rm,mat2,rng1,ecc,1)

                     else

                        call atima(apo,zpo,ek1,rm,mat2,rng1,
     &                             delt,ecc,1)

                     end if

                  do j = inmr - 1, 1, -1

                        ek2 = edmi * exp( edif * dble(j) )

                     if( ndedx .eq. 0 .or. ndedx .eq. 2 ) then

                        call spar(ipty,apo,zpo,ek2,rm,mat2,rng2,ecc,1)

                     else

                        call atima(apo,zpo,ek2,rm,mat2,rng2,
     &                             delt,ecc,1)

                     end if

                     rlet = ( ek1 - ek2 ) / ( rng1 - rng2 ) / 10.0

                     do k = 1, 43

                        if( rlet .le. elet(k) ) goto 500

                     end do

                        k = 43

  500                continue

                     plet = elet(k)
                     if( k .gt. 1 ) plet = plet - elet(k-1)

                     emet(k) = emet(k) + ( rng1 - rng2 ) / plet

                     rng1 = rng2
                     ek1  = ek2

                     if( j .eq. 1 ) then

                        rlet = ek1 / rng1 / 10.0

                        do k = 1, 43

                           if( rlet .le. elet(k) ) goto 501

                        end do

                           k = 43

  501                   continue

                        plet = elet(k)
                        if( k .gt. 1 ) plet = plet - elet(k-1)

                        emet(k) = emet(k) + rng1 / plet

                     end if


                  end do

*-----------------------------------------------------------------------

               end if

            end do

         end if

*-----------------------------------------------------------------------

 1000 continue

               rncnt(5) = mevent

               heat1 = heat1 / dble(mevent)

               heatt = 0.0
               xtott = 0.0
               xmttt = 0.0

            do i = 1, iznm

               heat(i) = heat(i) / dble(mevent)
               xmtt(i) = xtot(i) / dble(mevent)
               xtot(i) = xtot(i) / dble(mevent) * signon

               heatt = heatt + heat(i)
               xtott = xtott + xtot(i)
               xmttt = xmttt + xmtt(i)

            do j = 1, 20
            do k = 1, 104

               if( j .le. 19 ) then

                  cngd(icn+i,j,k) = cngd(icn+i,j,k) / dble(mevent)
     &                            * signon / 2.d0 / pir

               else

                  cngd(icn+i,j,k) = cngd(icn+i,j,k) / dble(mevent)
     &                            * signon
               end if

            end do
            end do
            end do

            do i = 1, iznm - 1
                     izmin = izms(izmr(i))
               do j = i + 1, iznm
                  if( izms(izmr(j)) .lt. izmin ) then
                     izmin = izms(izmr(j))
                     ii = izmr(i)
                     izmr(i) = izmr(j)
                     izmr(j) = ii
                  end if
               end do
            end do

            do i = 1, 43

                  emet(i) = emet(i) / dble(mevent)

            end do

*-----------------------------------------------------------------------

            write(io,'(''*** DDX'',
     &                 '' of Charged Particle and LET Distribution'')')
            write(io,'(''    calculated by PHITS with Nuclear Data'')')
            write(io,*)
            write(io,'(''    Target (Z,A) = '',2i6)') ntz, nta
            write(io,'(''         Energy  = '',1pg15.6,'' MeV'')') ein
            write(io,'('' total x-section = '',1pg15.6,'' mb'')') signon
            write(io,'(''          events = '',i9)') mevent

            write(io,*)
            write(io,'('' heat from nuclear data  = '',
     &                 1pg15.6,'' MeV'')') heat0
            write(io,'('' heat from simulation    = '',
     &                 1pg15.6,'' MeV'')') heat1

            write(io,'('' number of charged part. = '',i6)') iznm
            write(io,*)
            write(io,'(6x,''num    z    a     heat(MeV)      %'',
     &                 ''       x-sect(mb)   multiplicity'')')
            write(io,'(i8,3x,i3,i5,1pe15.4,3x,0pf5.3,1pe15.4,
     &                 1pe15.4)')
     &      (k,izms(izmr(k))/1000,
     &         izms(izmr(k))-izms(izmr(k))/1000*1000,
     &         heat(izmr(k)),heat(izmr(k))/heat1,
     &         xtot(izmr(k)),xmtt(izmr(k)),k=1,iznm)

            write(io,'(72(''-''))')
            write(io,'(11x,''total''3x,1pe15.4,8x,1pe15.4,1pe15.4)')
     &         heatt, xtott, xmttt

            write(io,*)

*-----------------------------------------------------------------------
*        target nucleus or proton
*-----------------------------------------------------------------------

               izz = ntz
               iaa = nta

            if( izz .eq. 1 .and. iaa .eq. 1 ) then
               ip = 1
            else if( izz .eq. 1 .and. iaa .eq. 2 ) then
               ip = 15
            else if( izz .eq. 1 .and. iaa .eq. 3 ) then
               ip = 16
            else if( izz .eq. 2 .and. iaa .eq. 3 ) then
               ip = 17
            else if( izz .eq. 2 .and. iaa .eq. 4 ) then
               ip = 18
            else
               ip = 19
            end if

               ik = kfft(ip)
               if( ip .eq. 19 ) ik = izz * 1000000 + iaa
               iz = izz
               rm = rmtyp(ip,ik)

*-----------------------------------------------------------------------
*        LET vs Energy
*-----------------------------------------------------------------------

               write(io,'(/''#'',78(''-''))')
               write(io,'( ''# LET vs Energy'')')
               write(io,'( ''#'',78(''-''))')

               write(io,'(/''x: Energy (MeV)'')')
               write(io,'( ''y: LET (keV/'',a1,''mu m)'')') yen
               write(io,'( ''p: xlog ylog afac(0.8) form(0.9)'')')

               write(io,'(/a1,''LET vs Energy for (Z,A) = ('',i3,'','',
     &                    i3,'' )'',a1)') cha, izz, iaa, cha

               write(io,'(/''h: x  y,l0  y,d0r  n'')')
               write(io,'( ''#  Energy       LET          DEDX'',
     &                     ''         Range'')')

                  nmes = 100

                  edmx = 1.e+4
                  edmi = 1.e-6
                  edlg = log( edmx / edmi ) / dble( nmes )

               do i = 1, nmes

                  ekin = edmi * exp( edlg * dble( i ) )

                  if( ndedx .eq. 0 ) then

                    call rainge(ekin*1.05,rng1,mat2,
     &                          ipart,kpart,jpart,rpart)
                    call rainge(ekin*0.95,rng2,mat2,
     &                          ipart,kpart,jpart,rpart)

                  else if( ndedx .eq. 1 .or. ndedx .eq. 3 ) then   ! S.Abe 2016/08/08

                     call atima(apr,zpr,ekin*1.05,rm,mat2,rng1,
     &                          delt,ecc,1)
                     call atima(apr,zpr,ekin*0.95,rm,mat2,rng2,
     &                          delt,ecc,1)

                  else if( ndedx .eq. 2 ) then

                     call spar(ipty,apr,zpr,ekin*1.05,rm,mat2,
     &                         rng1,ecc,1)
                     call spar(ipty,apr,zpr,ekin*0.95,rm,mat2,
     &                         rng2,ecc,1)

                  end if

                     if( rng1 .gt. rng2 ) then
                        rlet =  ekin * 0.1 / ( rng1 - rng2 ) / 10.0
                     else
                        rlet =  0.d0
                     end if

                     call dedxas(ekin,stpi,mat2,
     &                           ipart,kpart,jpart,rpart)

                     rletd = stpi / 10.0

                  write(io,'(6(1pe13.4))') ekin, rlet, rletd, rng1

               end do

*-----------------------------------------------------------------------
*        LET Distribution
*-----------------------------------------------------------------------

               write(io,'(/''#'',78(''-''))')
               write(io,'( ''# LET Distribution'')')
               write(io,'( ''#'',78(''-''))')

               write(io,'( '' newpage:'')')
               write(io,'(/''x: LET (keV/'',a1,''mu m)'')') yen
               write(io,'( ''y: Absorbed Dose (Gy '',
     &                                      a1,''mu m/ keV)'')') yen
               write(io,'( ''p: xlog ylog afac(0.8) form(0.9)'')')

               write(io,'(/a1,''LET Distribution'',a1)') cha, cha

               write(io,'(/''h: x  y(LET),hh0'')')

               sek = emet(1) * elet(1)

            do i = 2, 43

               sek = sek + emet(i) * ( elet(i) - elet(i-1) )

            end do

            do i = 1, 43

               emet(i) = emet(i) / sek

            end do

            do i = 1, 43

               write(io,'(6(1pe13.4))') elet(i), emet(i)

            end do

*-----------------------------------------------------------------------
*     DDX of Charged Particles
*-----------------------------------------------------------------------

         do k = 1, iznm

               izz = izms(izmr(k))/1000
               iaa = izms(izmr(k))-izms(izmr(k))/1000*1000

*-----------------------------------------------------------------------

               write(io,'(/''#'',78(''-''))')

            if( k .eq. 1 ) then
               write(io,'( ''# DDX of Charged Particle'')')
               write(io,'( ''#'',78(''-''))')
            end if

               write(io,'( '' newpage:'')')
               write(io,'(''#   no. ='',i5,3x,''z ='',i3,3x,
     &                     ''a ='',i3)') k, izz, iaa

*-----------------------------------------------------------------------

               write(io,'(/''x: Energy (MeV)'')')
               write(io,'( ''y: d'',a1,''sigma  / dE(mb/MeV)'')') yen
               write(io,'( ''p: xlog ylog afac(0.8) form(0.9)'')')

               ien = 104
               iem = 1

               write(io,'(/a1,''no. ='',i3,'',  (Z,A) = ('',i3,'','',
     &                    i3,'' )'',a1)') cha, k, izz, iaa, cha

               write(io,'(/''h: x  y(all angle),hh0'')')

            do ie = 1, 104

               write(io,'(6(1pe13.4))') engm(ie),
     &                    cngd(icn+izmr(k),20,ie)

               if( cngd(icn+izmr(k),20,ie) .gt. 0.0d0 .and.
     &             ie .lt. ien ) ien = ie
               if( cngd(icn+izmr(k),20,ie) .gt. 0.0d0 .and.
     &             ie .gt. iem ) iem = ie

            end do

               if( ien .gt. 1 ) then

                  enmin = engm(ien-1)

               else

                  enmin = engm(ien)

               end if

                  enmax = engm(iem)

            if( enmax .gt. enmin ) then

               write(io,'(/''p: xmin('',1pe13.4,'') xmax('',
     &               1pe13.4,'')'')') enmin, enmax

            end if

*-----------------------------------------------------------------------

               write(io,'(/''#'',78(''-''))')
               write(io,'( '' newpage:'')')
               write(io,'(/''x: Energy (MeV)'')')
               write(io,'( ''y: d^2'',a1,
     &         ''sigma  / dEd'',a1,''Omega(mb/MeV sr)'')') yen,yen
               write(io,'( ''p: xlog ylog afac(0.8) form(0.9)'')')

               ien = 104
               iem = 1

               write(io,'(/a1,''no. ='',i3,'',  (Z,A) = ('',i3,'','',
     &                    i3,'' )'',a1)') cha, k, izz, iaa, cha

*-----------------------------------------------------------------------

               write(io,'(/''h: x n n n y('',f5.1''),hh0r n'')')
     &                    angd(4)
               write(io,'( ''#  energy   ang ='',5(f5.1,8x))')
     &                   (angd(j),j=1,5)

            do ie = 1, 104

               write(io,'(6(1pe13.4))') engm(ie),
     &                   (cngd(icn+izmr(k),j,ie),j=1,5)

               if( cngd(icn+izmr(k),4,ie) .gt. 0.0d0 .and.
     &             ie .lt. ien ) ien = ie
               if( cngd(icn+izmr(k),4,ie) .gt. 0.0d0 .and.
     &             ie .gt. iem ) iem = ie

            end do

*-----------------------------------------------------------------------

               write(io,'(/''h: x n n n n y('',f5.1''),hh0b'')')
     &                    angd(10)
               write(io,'( ''#  energy   ang ='',5(f5.1,8x))')
     &                   (angd(j),j=6,10)

            do ie = 1, 104

               write(io,'(6(1pe13.4))') engm(ie),
     &                   (cngd(icn+izmr(k),j,ie),j=6,10)

               if( cngd(icn+izmr(k),10,ie) .gt. 0.0d0 .and.
     &             ie .lt. ien ) ien = ie
               if( cngd(icn+izmr(k),10,ie) .gt. 0.0d0 .and.
     &             ie .gt. iem ) iem = ie

            end do

*-----------------------------------------------------------------------

               write(io,'(/''c: x n n n n n'')')
               write(io,'( ''#  energy   ang ='',5(f5.1,8x))')
     &                   (angd(j),j=11,15)

            do ie = 1, 104

               write(io,'(6(1pe13.4))') engm(ie),
     &                   (cngd(icn+izmr(k),j,ie),j=11,15)

            end do

*-----------------------------------------------------------------------

               write(io,'(/''h: x y('',f5.1''),hh0g n n n'')')
     &                    angd(16)
               write(io,'( ''#  energy   ang ='',5(f5.1,8x))')
     &                   (angd(j),j=16,19)

            do ie = 1, 104

               write(io,'(6(1pe13.4))') engm(ie),
     &                   (cngd(icn+izmr(k),j,ie),j=16,19)

               if( cngd(icn+izmr(k),16,ie) .gt. 0.0d0 .and.
     &             ie .lt. ien ) ien = ie
               if( cngd(icn+izmr(k),16,ie) .gt. 0.0d0 .and.
     &             ie .gt. iem ) iem = ie

            end do

*-----------------------------------------------------------------------

               if( ien .gt. 1 ) then

                  enmin = engm(ien-1)

               else

                  enmin = engm(ien)

               end if

                  enmax = engm(iem)

            if( enmax .gt. enmin ) then

               write(io,'(/''p: xmin('',1pe13.4,'') xmax('',
     &               1pe13.4,'')'')') enmin, enmax

            end if

*-----------------------------------------------------------------------

         end do

*-----------------------------------------------------------------------

      return

  999 continue

*-----------------------------------------------------------------------

      return
      end


