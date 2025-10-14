************************************************************************
*                                                                      *
      subroutine ovly13
*                                                                      *
*        modified by K.Niita on 2005/11/01                             *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              calculate thin target reactions                         *
*                                                                      *
*                                                                      *
************************************************************************
!$    use omp_lib
C for USE_MOD_COUNTER
      use mod_counter, only: ALLOCATE_EVTS,DEALLOCATE_EVTS,INIT_EVTS
     &     ,rncnt,rnint,rnintr,rnpnt,rnpntr
C for REDUCTION_COUNTER
!$   &                      ,rncnt2,rnint2,rnintr2,rnpnt2,rnpntr2
      use MMBANKMOD  !FURUTA
      use MEMBANKMOD !FURUTA

      use GGBANKMOD  !FURUTA
      use GGMBANKMOD !FURUTA
      use moddas_material
      use moddas_source ! T.Sato 2022/12/07
      use liboutmod, only: acewrite ! frtati 2022/12/28
*-----------------------------------------------------------------------

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param00.inc'
      include 'param02.inc'
      include 'param.inc'
      include 'param-physcnst.inc'

*-----------------------------------------------------------------------

      common /inout/  ins, ios

*-----------------------------------------------------------------------

      common /mpi00/ npe, me

*-----------------------------------------------------------------------


*-----------------------------------------------------------------------
      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
!$OMP THREADPRIVATE(/icomon/)

      common /ccggg/  icgg
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

      common /isorse/ ngrp(isrc), ngei(isrc), ngea(isrc), ngfe(isrc),
     &                ngft(isrc), ngll(isrc), ngpi(isrc), ngpw(isrc)

*-----------------------------------------------------------------------

      common /eparm/  esmax, esmin, emin(20)
      common /cparm/  maxbch,maxcas
      common /qparm/  ielas,icasc,iqstep,lvlopt,igamma
      common /bparm/  andt,jevap,npidk
      common /tcntl/  icntl, inucr
      common /emode/  emodem, ge1, ge2, iemode
      common /ccxsm/  icxsni, icxspi
      common /pnint/  ipnint    ! S.H. 2020.11.27
      dimension sigphoton(5)

*-----------------------------------------------------------------------

      common /comps1/ mstapr, massta, msprpr, masspr
      common /comps2/ sigela, signon, fissx

      common /ptname/ pname(20), ipln(20)
      character       pname*8

      common /analevt/ irunp

*-----------------------------------------------------------------------

      common /mathzn/ mathz, mathn, jcoll, kcoll
!$OMP THREADPRIVATE(/mathzn/)

*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
      common /tall80/ itdpa
      common /tall81/ iteth
*-----------------------------------------------------------------------

      common /cMeVperu/ iMeVperu

      common /dedxfac/ dedxfd
!$OMP THREADPRIVATE(/dedxfac/)

      dimension evle(11)
      data evle /150.0, 250.0, 350.0, 500.0, 600.0, 700.0, 800.0,
     &           1000.0, 1500.0, 2000.0, 3000.0/

      dimension angds(1000)
      dimension dsiga(0:2,1000), dsigp(1000)
      dimension sek(0:2)

      character yen*1
      yen  = char(92)

*-----------------------------------------------------------------------
*        only for executable PE
*-----------------------------------------------------------------------

         if( me .eq. 0 .and. npe .gt. 1 ) return

*-----------------------------------------------------------------------
*     choice of nuclear reaction ( inucr )
*-----------------------------------------------------------------------
*        inucr = 11 ; for COMPAS
*-----------------------------------------------------------------------
*                 anal-005.f for COMPAS out put
*                            only for proton or neutron incident
*
*                    [parameters]
*                       maxcas : number of events
*                    [source]    s-type = 1
*                       proj   : projectile
*                       e0     : incident energy
*                    [material]
*                     mat[1]
*                       first nucleus : target
*
*-----------------------------------------------------------------------
*        inucr = 1 ; double differrential cross section
*-----------------------------------------------------------------------
*                 anal-002.f for double differential cross section
*                 anal-004.f for COMPAS out put
*                            only for proton or neutron incident
*
*                    [parameters]
*                       maxcas : number of events
*                    [source]    s-type = 1
*                       proj   : projectile
*                       e0     : incident energy
*                    [material]
*                     mat[1]
*                       first nucleus : target
*
*-----------------------------------------------------------------------
*        inucr = 2 ; total reaction, elastic and inelastic
*-----------------------------------------------------------------------
*
*                    [parameters]
*                       inucl  : =0 linear energy mesh
*                                =1 logarithm energy mesh
*                       maxbch : number of energy mesh
*                       emin(itype) : min energy mesh
*                    [source]    s-type = 1
*                       proj   : projectile
*                       e0     : max energy mesh
*                    [material]
*                     mat[1]
*                       first nucleus : target
*
*-----------------------------------------------------------------------
*        inucr = 3 ; inelastic cross section in simulation
*-----------------------------------------------------------------------
*
*                    [parameters]
*                       maxcas : number of events
*                       inucl  : =0 linear energy mesh
*                                =1 logarithm energy mesh
*                       maxbch : number of energy mesh
*                       emin(itype) : min energy mesh
*                    [source]    s-type = 1
*                       proj   : projectile
*                       e0     : max energy mesh
*                    [material]
*                     mat[1]
*                       first nucleus : target
*
*-----------------------------------------------------------------------
*        inucr = 4 ; angular distribution of elastic collisions
*-----------------------------------------------------------------------
*
*                    [parameters]
*                       maxcas : number of events for simulation
*                       ielms  : mesh point of angular distribution
*                       inucl  : =0 x-axis cm angle in degree
*                                =1 x-axis con(theta) in cm
*                    [source]    s-type = 1
*                       proj   : projectile ( proton or neutron )
*                       e0     : incident energy
*                    [material]
*                     mat[1]
*                       first nucleus : target
*
*-----------------------------------------------------------------------
*        inucr = 5 ; pp, np, pi-p cross section
*-----------------------------------------------------------------------
*
*                    [parameters]
*                       inucl  : =0 linear energy mesh
*                                =1 logarithm energy mesh
*                       maxbch : number of energy mesh
*                       emin(itype) : min energy mesh
*                    [source]    s-type = 1
*                       proj   : projectile
*                       e0     : max energy mesh
*
*-----------------------------------------------------------------------
*        inucr = 6 ; pp, np, pi-p cross section in simulation
*-----------------------------------------------------------------------
*
*                    [parameters]
*                       maxcas : number of events for simulation
*                       inucl  : =0 linear energy mesh
*                                =1 logarithm energy mesh
*                       maxbch : number of energy mesh
*                       emin(itype) : min energy mesh
*                    [source]    s-type = 1
*                       proj   : projectile
*                       e0     : max energy mesh
*
*-----------------------------------------------------------------------
*        inucr = 7-9 ; elastic + inelastic nuclear reactions
*        inucr = 7 ; transverse energy
*        inucr = 8 ; recoil energy
*        inucr = 9 ; DPA cross section
*-----------------------------------------------------------------------
*        inucr = 12; Kerma
*        inucr = 13; Photon energy
*        inucr = 14; Photon Kerma
*        inucr = 16; dE/dx
*        inucr = 100; nuclear data library
*-----------------------------------------------------------------------

*        inucr = 15; inelastic
*-----------------------------------------------------------------------
*
*                    [parameters]
*                       maxcas : number of events for simulation
*                       inucl  : =0 linear energy mesh
*                                =1 logarithm energy mesh
*                       maxbch : number of energy mesh
*                       emin(itype) : min energy mesh
*                    [source]    s-type = 1
*                       proj   : projectile
*                       e0     : max energy mesh
*
*-----------------------------------------------------------------------

            if( inucr .lt. 1 .or. (inucr.gt.16.and.inucr.ne.100) ) then
               write(6,*) ' ******* inucr is wrong ***'
               ierr = 1
               goto 999
            end if

*-----------------------------------------------------------------------
*     output unit ( 61 - 65 )
*-----------------------------------------------------------------------

         io = 61

         if( inucr.ne.100 ) then ! frtati 2022/12/28
           open(io, file = chfn(11), status = 'unknown' )
         end if

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
*        dynamical allocattion and initialization for GGBANK & GGMBANK
*-----------------------------------------------------------------------

               call ALLOCATE_MMBANK   !FURUTA
               call ALLOCATE_MEMBANK  !FURUTA

               if(icgg.ne.0)then
                 call ALLOCATE_GGBANK !FURUTA
                 call INIT_GGBANK     !FURUTA
               endif
               call ALLOCATE_GGMBANK  !FURUTA
               call INIT_GGMBANK      !FURUTA
               call ALLOCATE_EVTS     !FURUTA20210506
               call INIT_EVTS         !FURUTA20210506

*-----------------------------------------------------------------------
               if( inucr.eq.100 ) then
                 call acewrite
                 return
               end if

*-----------------------------------------------------------------------
*        pick up the target nucleus from the first material
*-----------------------------------------------------------------------

               icl   = 1
               mat   = 1
               ireg  = 1

               lemm  = nint( dnel_das(kmat0+mat) )
               hydro = denh_das(kmat0+mat)
               jimat = lemm

               if( hydro .gt. 0.0 ) jimat = jimat + 1

            if( jimat .eq. 0 ) then

               write(6,*) ' **** There is no target'
               ierr = 1
               goto 999

            end if

            if( inucr .ne. 11 ) jimat = 1

*-----------------------------------------------------------------------
*        for varius nucleus  for inucr = 11
*-----------------------------------------------------------------------

            inumc = 0

      do 5200 imat = 1, jimat

               lem = imat

            if( lem .le. lemm ) then

               zpr = zz_das(kmat(mat)+lem)
               apr = a_das(kmat(mat)+lem)

            else

               apr = 1.0
               zpr = 1.0

            end if

               nta = nint( apr )
               ntz = nint( zpr )

*-----------------------------------------------------------------------
*        for proton and neutron incident for inucr = 11
*-----------------------------------------------------------------------

               jipnc = 1

               if( inucr .eq. 11 ) jipnc = 2

      do 5100 ipnc = 1, jipnc

*-----------------------------------------------------------------------
*        incident particle
*-----------------------------------------------------------------------

            if( inucr .ne. 11 ) then

             ityp = istyp(imsrc)
             ktyp = inkf0(imsrc)
             if(jstyp(imsrc).eq.4) then ! energy directly specified by user, T.Sato 2022/12/07
              maxbch = ngrp(imsrc)+1
              emin(ityp) = egmin(1)
              emax = egmax(ngrp(imsrc))
             else
              emax = se0(imsrc)
             endif
             jpz = 0
             jpa = 0

             if( ityp .ge. 15 .and. ityp .le. 19 ) then

                jpz = ichgf(ityp,ktyp)
                jpa = ibryf(ityp,ktyp)

                if(iMeVperu.eq.1) then ! T.Sato 2021/10/03
                 emin(ityp) = emin(ityp) * jpa
                 emax = emax * jpa
                endif

             end if

            else

               if( ipnc .eq. 1 ) then

                  ityp = 1
                  ktyp = 2212

               else

                  ityp = 2
                  ktyp = 2112

               end if

            end if

*-----------------------------------------------------------------------

         if( jpa .eq. 0 ) then

            if( nta .gt. 1 ) then

               siggin  = sigg(ksig(mat)+lem)
               siggeo  = siggin/den_das(kmat(mat)+lem)*1000.0
               bmax10  = sqrt( siggeo / 10.0 / pi )

            else

               siggeo  = 200.0
               bmax10  = sqrt( siggeo / 10.0 / pi )

            end if

         else

               ap = jpa
               zp = jpz
               at = nta
               zt = ntz

               call sighi(ap,zp,emax,at,zt,ssigne,ssigel,bmax10)

               siggeo  = ssigne * 1000.0

         end if

*-----------------------------------------------------------------------
*        nejj   ; energy mesh
*        eini   ; minimum energy
*        ielg   ; 1-log, 0-liniear
*        mevent ; number of events
*-----------------------------------------------------------------------

            mevent = maxcas

         if( inucr .eq. 11 ) then

            nejj = 11

         end if

         if( inucr .eq. 1 ) then

            ielg = mstz(25)
            nejj = 1
            eini = emax
            emax = emax

         else if( inucr .eq. 2 .or. inucr .eq. 3 ) then

            ielg = mstz(25)
            nejj = max( 2, maxbch )
            eini = emin(ityp)
            emax = emax

         else if( inucr .ge. 5 .and. inucr .le. 9 .or.
     &            inucr .eq. 12 .or. inucr .eq. 13 .or.
     &            inucr .eq. 14 .or. inucr .eq. 15 .or.
     &            inucr .eq. 16) then  ! T.Sato 2022/12/06, dE/dx mode

            ielg = mstz(25)
            nejj = max( 2, maxbch )
            eini = emin(ityp)
            emax = emax

         else if( inucr .eq. 4 ) then

            if( ityp .gt. 2 ) then

               write(6,*) ' **** Elastic is only for nucleon'
               ierr = 1
               goto 999

            end if

            ielg = 0
            nejj = 1
            eini = emax
            emax = emax

         end if

*-----------------------------------------------------------------------

         if( inucr .eq. 1 ) then

            write(6,*)
            write(6,'(''*** Double Differential Cross Sections'')')
            write(6,*)
            write(6,'(''          ityp   = '',2i6)') ityp
            write(6,'(''          ktyp   = '',2i9)') ktyp
            if(ityp.eq.19)
     &      write(6,'('' Projectile(A,Z) = '',2i6)') jpa,jpz
            if(ityp.ge.15.and.ityp.le.19.and.iMeVperu.eq.1) then
             write(6,'(''    Max Energy   = '',1pg15.6,'' MeV/n'')')
     &       emax/jpa
            else
             write(6,'(''    Max Energy   = '',1pg15.6,'' MeV'')') emax
            endif
            write(6,'(''    Target (A,Z) = '',2i6)') nta, ntz

         end if

*-----------------------------------------------------------------------

         if( inucr .eq. 2 ) then

            write(io,*)
            write(io,'(''*** Total, Elastic and Non-Elastic Cross '',
     &                 ''Sections ***'')')
            write(io,*)
            write(io,'(''  Particle type  =    '',a8)') pname(ityp)
            if(ityp.eq.19)
     &      write(io,'('' Projectile(A,Z) = '',2i6)') jpa,jpz
            write(io,'(''    Target (A,Z) = '',2i6)') nta, ntz
            if(ityp.ge.15.and.ityp.le.19.and.iMeVperu.eq.1) then
             write(io,'(''    Min Energy   = '',1pg15.6,'' MeV/n'')')
     &       eini/jpa
             write(io,'(''    Max Energy   = '',1pg15.6,'' MeV/n'')')
     &       emax/jpa
            else
             write(io,'(''    Min Energy   = '',1pg15.6,'' MeV'')') eini
             write(io,'(''    Max Energy   = '',1pg15.6,'' MeV'')') emax
            endif
            write(io,'(''    energy mesh  = '',i6)') nejj
            if( ielg .eq. 0 ) then
            write(io,'(''    energy mesh  =    linear'')')
            else
            write(io,'(''    energy mesh  =    log'')')
            end if

         end if

*-----------------------------------------------------------------------

         if( inucr .eq. 15 ) then

            write(io,*)
            write(io,'(''*** Total, Elastic, Non-Elastic '',
     &                 ''and Inelastic Cross Sections ***'')')
            write(io,*)
            write(io,'(''  Particle type  =    '',a8)') pname(ityp)
            if(ityp.eq.19)
     &      write(io,'('' Projectile(A,Z) = '',2i6)') jpa,jpz
            write(io,'(''    Target (A,Z) = '',2i6)') nta, ntz
            if(ityp.ge.15.and.ityp.le.19.and.iMeVperu.eq.1) then
             write(io,'(''    Min Energy   = '',1pg15.6,'' MeV/n'')')
     &       eini/jpa
             write(io,'(''    Max Energy   = '',1pg15.6,'' MeV/n'')')
     &       emax/jpa
            else
             write(io,'(''    Min Energy   = '',1pg15.6,'' MeV'')') eini
             write(io,'(''    Max Energy   = '',1pg15.6,'' MeV'')') emax
            endif
            write(io,'(''    energy mesh  = '',i6)') nejj
            if( ielg .eq. 0 ) then
            write(io,'(''    energy mesh  =    linear'')')
            else
            write(io,'(''    energy mesh  =    log'')')
            end if

         end if

*-----------------------------------------------------------------------

         if( inucr .eq. 3 ) then

            write(io,*)
            write(io,'(''*** Total, Elastic and Non-Elastic Cross '',
     &                 ''Sections by Simulations ***'')')
            write(io,*)
            write(io,'(''  Particle type  =    '',a8)') pname(ityp)
            if(ityp.eq.19)
     &      write(io,'('' Projectile(A,Z) = '',2i6)') jpa,jpz
            write(io,'(''    Target (A,Z) = '',2i6)') nta, ntz
            if(ityp.ge.15.and.ityp.le.19.and.iMeVperu.eq.1) then
             write(io,'(''    Min Energy   = '',1pg15.6,'' MeV/n'')')
     &       eini/jpa
             write(io,'(''    Max Energy   = '',1pg15.6,'' MeV/n'')')
     &       emax/jpa
            else
             write(io,'(''    Min Energy   = '',1pg15.6,'' MeV'')') eini
             write(io,'(''    Max Energy   = '',1pg15.6,'' MeV'')') emax
            endif
            write(io,'(''    energy mesh  = '',i6)') nejj
            if( ielg .eq. 0 ) then
            write(io,'(''    energy mesh  =    linear'')')
            else
            write(io,'(''    energy mesh  =    log'')')
            end if

         end if

*-----------------------------------------------------------------------

         if( inucr .eq. 4 ) then

            write(io,*)
            write(io,'(''*** Angular Distribution of  '',
     &                 ''Elastic Collisions ***'')')
            write(io,*)
            write(io,'(''  Particle type  =    '',a8)') pname(ityp)
            if(ityp.eq.19)
     &      write(io,'('' Projectile(A,Z) = '',2i6)') jpa,jpz
            write(io,'(''    Target (A,Z) = '',2i6)') nta, ntz
            if(ityp.ge.15.and.ityp.le.19.and.iMeVperu.eq.1) then
             write(io,'(''    Max Energy   = '',1pg15.6,'' MeV/n'')')
     &       emax/jpa
            else
             write(io,'(''    Max Energy   = '',1pg15.6,'' MeV'')') emax
            endif
            write(io,*)

            write(io,'(''p: form(0.8) xfac(1.2)'')')
            write(io,'(''p: scal(0.8) xorg(0.2) yorg(.1)'')')
            write(io,'(''p: nofr noms nocn'')')

            if( mstz(25) .eq. 0 ) then

               write(io,'(''p: ylog xmax(180) ymin(1.e-3)'')')
               write(io,'(''p: ymax(1.e5) nosp'')')
               write(io,'(''x: CM-angle (degree)'')')
               write(io,'(''y: Cross Section (mb/sr)'')')

            else

               write(io,'(''p: ylog xmin(-1) xmax(1) nosx'')')
               write(io,'(''x: cos('',a1,''theta) in CM'')') yen
               write(io,'(''y: Probability'')')

            end if

         end if

*-----------------------------------------------------------------------

         if( inucr .eq. 5 ) then

            write(io,*)
            write(io,'(''*** Total, Elastic and Non-Elastic Cross '',
     &                 ''Sections on Hydrogen by JAM ***'')')
            write(io,*)
            write(io,'(''  Particle type  =    '',a8)') pname(ityp)
            if(ityp.eq.19)
     &      write(io,'('' Projectile(A,Z) = '',2i6)') jpa,jpz
            write(io,'(''    target       =    '',a8)') pname(1)
            if(ityp.ge.15.and.ityp.le.19.and.iMeVperu.eq.1) then
             write(io,'(''    Min Energy   = '',1pg15.6,'' MeV/n'')')
     &       eini/jpa
             write(io,'(''    Max Energy   = '',1pg15.6,'' MeV/n'')')
     &       emax/jpa
            else
             write(io,'(''    Min Energy   = '',1pg15.6,'' MeV'')') eini
             write(io,'(''    Max Energy   = '',1pg15.6,'' MeV'')') emax
            endif
            write(io,'(''    energy mesh  = '',i6)') nejj
            if( ielg .eq. 0 ) then
            write(io,'(''    energy mesh  =    linear'')')
            else
            write(io,'(''    energy mesh  =    log'')')
            end if

            write(io,*)
            write(io,'(''x: E (MeV)'')')
            write(io,'(''y: Cross Section (mb)'')')
            if( ielg .eq. 1 ) write(io,'(''p: xlog'')')
            write(io,'(''h: x y(total),l0 y(elast),d0 y(nonela),m0'')')
            write(io,'(''#   energy     JAM total    JAM elast  '',
     &                 ''  JAM inela'')')


         end if

*-----------------------------------------------------------------------

         if( inucr .eq. 6 ) then

            write(io,*)
            write(io,'(''*** Total, Elastic and Non-Elastic Cross '',
     &                 ''Sections on Hydrogen by JAM '',
     &                 ''and Simulation ***'')')
            write(io,*)
            write(io,'(''    particle     =    '',a8)') pname(ityp)
            write(io,'(''    target       =    '',a8)') pname(1)
            if(ityp.ge.15.and.ityp.le.19.and.iMeVperu.eq.1) then
             write(io,'(''    Min Energy   = '',1pg15.6,'' MeV/n'')')
     &       eini/jpa
             write(io,'(''    Max Energy   = '',1pg15.6,'' MeV/n'')')
     &       emax/jpa
            else
             write(io,'(''    Min Energy   = '',1pg15.6,'' MeV'')') eini
             write(io,'(''    Max Energy   = '',1pg15.6,'' MeV'')') emax
            endif
            write(io,'(''    energy mesh  = '',i6)') nejj
            if( ielg .eq. 0 ) then
            write(io,'(''    energy mesh  =    linear'')')
            else
            write(io,'(''    energy mesh  =    log'')')
            end if

            write(io,*)
            write(io,'(''x: E (MeV)'')')
            write(io,'(''y: Cross Section (mb)'')')
            if( ielg .eq. 1 ) write(io,'(''p: xlog'')')
            write(io,'(''h: x y(total),l0 y(inela),d0'',
     &                 '' y(Simu.inela),ur0 n n n'')')
            write(io,'(''#energy/srt    JAM total    JAM inela  '',
     &                 ''  Simu. inela  elas(%)      inela1(%)=3'',
     &                 ''  inela2(%)>3'')')

         end if

*-----------------------------------------------------------------------
*        for anal-006.f
*-----------------------------------------------------------------------

         if( inucr .eq. 7 ) then

            write(io,*)
            write(io,'(''*** Transverse Energy '')')
            write(io,*)
            write(io,'(''  Particle type  =    '',a8)') pname(ityp)
            if(ityp.eq.19)
     &      write(io,'('' Projectile(A,Z) = '',2i6)') jpa,jpz
            write(io,'(''    Target (A,Z) = '',2i6)') nta, ntz
            if(ityp.ge.15.and.ityp.le.19.and.iMeVperu.eq.1) then
             write(io,'(''    Min Energy   = '',1pg15.6,'' MeV/n'')')
     &       eini/jpa
             write(io,'(''    Max Energy   = '',1pg15.6,'' MeV/n'')')
     &       emax/jpa
            else
             write(io,'(''    Min Energy   = '',1pg15.6,'' MeV'')') eini
             write(io,'(''    Max Energy   = '',1pg15.6,'' MeV'')') emax
            endif
            write(io,'(''    energy mesh  = '',i6)') nejj
            if( ielg .eq. 0 ) then
            write(io,'(''    energy mesh  =    linear'')')
            else
            write(io,'(''    energy mesh  =    log'')')
            end if

            write(io,*)
            write(io,'(''x: E (MeV)'')')
            write(io,'(''y: Transverse Energy / E'')')
            if( ielg .eq. 1 ) write(io,'(''p: xlog'')')
            write(io,'(''h: x y(els),l0 y(nucleus),d0 y(nucleon),m0 '',
     &                      ''y(tot),l0'')')

         end if

*-----------------------------------------------------------------------
*           for anal-007.f
*-----------------------------------------------------------------------

         if( inucr .eq. 8 ) then

            write(io,*)
            write(io,'(''*** Recoil Energy '')')
            write(io,*)
            write(io,'(''  Particle type  =    '',a8)') pname(ityp)
            if(ityp.eq.19)
     &      write(io,'('' Projectile(A,Z) = '',2i6)') jpa,jpz
            write(io,'(''    Target (A,Z) = '',2i6)') nta, ntz
            if(ityp.ge.15.and.ityp.le.19.and.iMeVperu.eq.1) then
             write(io,'(''    Min Energy   = '',1pg15.6,'' MeV/n'')')
     &       eini/jpa
             write(io,'(''    Max Energy   = '',1pg15.6,'' MeV/n'')')
     &       emax/jpa
            else
             write(io,'(''    Min Energy   = '',1pg15.6,'' MeV'')') eini
             write(io,'(''    Max Energy   = '',1pg15.6,'' MeV'')') emax
            endif
            write(io,'(''    energy mesh  = '',i6)') nejj
            if( ielg .eq. 0 ) then
            write(io,'(''    energy mesh  =    linear'')')
            else
            write(io,'(''    energy mesh  =    log'')')
            end if

            write(io,*)
            write(io,'(''x: E (MeV)'')')
            write(io,'(''y: Recoil Energy (MeV)'')')
            if( ielg .eq. 1 ) write(io,'(''p: xlog'')')
            write(io,'(''h: x y(d),d0 y(t),p0 y(h3),m0 '',
     &      ''y(alpha),l0 y(tot),l0 y(recoil),l0r y(ex*0.1),l0b'')')

         end if

*-----------------------------------------------------------------------
*           for anal-008
*-----------------------------------------------------------------------

         if( inucr .eq. 9 ) then

            write(io,*)
            write(io,'(''*** DPA cross section '')')
            write(io,*)
            write(io,'(''  Particle type  =    '',a8)') pname(ityp)
            if(ityp.eq.19)
     &      write(io,'('' Projectile(A,Z) = '',2i6)') jpa,jpz
            write(io,'(''    Target (A,Z) = '',2i6)') nta, ntz
            if(ityp.ge.15.and.ityp.le.19.and.iMeVperu.eq.1) then
             write(io,'(''    Min Energy   = '',1pg15.6,'' MeV/n'')')
     &       eini/jpa
             write(io,'(''    Max Energy   = '',1pg15.6,'' MeV/n'')')
     &       emax/jpa
            else
             write(io,'(''    Min Energy   = '',1pg15.6,'' MeV'')') eini
             write(io,'(''    Max Energy   = '',1pg15.6,'' MeV'')') emax
            endif
            write(io,'(''    energy mesh  = '',i6)') nejj
            if( ielg .eq. 0 ) then
            write(io,'(''    energy mesh  =    linear'')')
            else
            write(io,'(''    energy mesh  =    log'')')
            end if
            if( itdpa .eq. 0 ) then
            write(io,'(''    damage approximation  = NRT'')')
            else
            write(io,'(''    damage approximation  = arc'')')
            end if
            if( iteth .eq. 0 ) then
            write(io,'(''    threshold energy  = original'')')
            else
            write(io,'(''    threshold energy  = new'')')
            end if

            write(io,*)
            if(ityp.ge.15.and.ityp.le.19.and.iMeVperu.eq.1) then
             write(io,'(''x: E (MeV/n)'')')
            else
             write(io,'(''x: E (MeV)'')')
            endif
            write(io,'(''y: DPA Cross Section (b)'')')
            if( ielg .eq. 1 ) write(io,'(''p: xlog'')')
            write(io,'(''h: x y(qoul),l0b y(nucl),l0r y(total),l0'')')
         end if

*-----------------------------------------------------------------------

         if( inucr .eq. 12 .or. inucr. eq. 13 .or. inucr. eq. 14 ) then

               write(io,'(''x: E (MeV)'')')

            if( ielg .eq. 1 ) then
               write(io,'(''p: xlog ylog'')')
            else
               write(io,'(''p: xlin ylin'')')
            end if

            if( inucr .eq. 12 ) then

               write(io,'(''y: heat (MeV or MeV*b)'')')
               write(io,'(''h: x y(calc-MeV),l0btt y(data-MeV),l0b'',
     &                    '' y(calc-MeV*b),l0rtt y(data-MeV*b),l0r'')')

            else if( inucr .eq. 13 ) then

               write(io,'(''y: Photon (MeV)'')')
               write(io,'(''h: x y(calc),l0btt y(data),l0r'')')

            else if( inucr .eq. 14 ) then

               write(io,'(''y: heat (MeV or MeV*b)'')')
               write(io,'(''h: x y(data-MeV),l0b y(data-MeV*b),l0r'',
     &                    '' ny(sigma)'')')

            end if

         end if

         if( inucr .eq. 16 ) then ! T.Sato 2022/12/06, dE/dx mode
            if(ityp.eq.12.or.ityp.eq.13) then
             write(*,'("Error: dE/dx cannot be output for electron",
     &       " & positron because it depends on cut-off energy")')
             stop
            endif
            write(io,*)
            write(io,'(''*** dE/dx in keV/um *** '')')
            write(io,*)
            write(io,'(''  Particle type  =    '',a8)') pname(ityp)
            if(ityp.eq.19)
     &      write(io,'('' Projectile(A,Z) = '',2i6)') jpa,jpz
            if(ityp.ge.15.and.ityp.le.19.and.iMeVperu.eq.1) then
             write(io,'(''    Min Energy   = '',1pg15.6,'' MeV/n'')')
     &       eini/jpa
             write(io,'(''    Max Energy   = '',1pg15.6,'' MeV/n'')')
     &       emax/jpa
            else
             write(io,'(''    Min Energy   = '',1pg15.6,'' MeV'')') eini
             write(io,'(''    Max Energy   = '',1pg15.6,'' MeV'')') emax
            endif
            write(io,'(''    energy mesh  = '',i6)') nejj
            if( ielg .eq. 0 ) then
            write(io,'(''    energy mesh  =    linear'')')
            else
            write(io,'(''    energy mesh  =    log'')')
            end if
            write(io,*)
            if(ityp.ge.15.and.ityp.le.19.and.iMeVperu.eq.1) then
             write(io,'(''x: E (MeV/n)'')')
            else
             write(io,'(''x: E (MeV)'')')
            endif
            write(io,'(''y: dE/dx (keV/um)'')')
            if( ielg .eq. 1 ) write(io,'(''p: xlog ylog'')')
            write(io,'(''h: x y(dE/dx),l0'')')
         end if


*-----------------------------------------------------------------------

         if( inucr .eq. 1 ) then

            write(6,*)
            write(6,'(''   Geom. inelastic x-section  = '',
     &                      1pg15.6)') siggeo
            write(6,'(''   Geom. max impact parameter = '',
     &                      1pg15.6)') bmax10
            write(6,*)

         else if( inucr .eq. 2 ) then

          if( ityp .eq. 14 ) then

            write(io,*)
            write(io,'(''x: E (MeV)'')')
            write(io,'(''y: Cross Section (mb)'')')
            if( ielg .eq. 1 ) write(io,'(''p: xlog ylog'')')

           if ( ipnint.gt.0 ) then

            write(io,'(''h: x y(total-nuclear),l0 y(GDR),dr0'',
     &'' y(quasi-deu),ub0 y(pion-pro),mg0 y(NRF),qrr0'')')
            write(io,'(''# energy       total-nucl   GDR          '',
     &''quasi-deu    pion-pro     NRF'')')

           else

            write(io,'(''h: x y(total-atomic),l0 y(Incoherent),dr0'',
     &'' y(Coherent),ub0 y(ph-el),mg0 y(pair-pro),qrr0'')')
            write(io,'(''# energy       total-atomic Incoherent   '',
     &''Coherent     photo-elect  pair-production'')')

           end if


          else

            write(io,*)
            write(io,'(''   Geom. inelastic x-section  = '',
     &                      1pg15.6)') siggeo
            write(io,'(''   Geom. max impact parameter = '',
     &                      1pg15.6)') bmax10
            write(io,*)
            if(ityp.ge.15.and.ityp.le.19.and.iMeVperu.eq.1) then
             write(io,'(''x: E (MeV/n)'')')
            else
             write(io,'(''x: E (MeV)'')')
            endif
            write(io,'(''y: Cross Section (mb)'')')
            if( ielg .eq. 1 ) write(io,'(''p: xlog'')')
            write(io,'(''h: x y(total),l0 y(nonela),d0 y(elast),m0'')')
            write(io,'(''#   energy           total    inelastic'',
     &                 ''      elastic'')')

          end if

         else if( inucr .eq. 15 ) then
            write(io,*)
            write(io,'(''   Geom. inelastic x-section  = '',
     &                      1pg15.6)') siggeo
            write(io,'(''   Geom. max impact parameter = '',
     &                      1pg15.6)') bmax10
            write(io,*)
            if(ityp.ge.15.and.ityp.le.19.and.iMeVperu.eq.1) then
             write(io,'(''x: E (MeV/n)'')')
            else
             write(io,'(''x: E (MeV)'')')
            endif
            write(io,'(''y: Cross Section (mb)'')')
            if( ielg .eq. 1 ) write(io,'(''p: xlog'')')
            write(io,'(''h: x y(total),l0 y(nonela),d0 y(elast),m0'',
     &                 '' y(inela),d0r'')')
            write(io,'(''#   energy     total        nonela'',
     &                 ''       elast        inela'')')

         else if( inucr .eq. 3 ) then

            write(io,*)
            write(io,'(''   Geom. inelastic x-section  = '',
     &                      1pg15.6)') siggeo
            write(io,'(''   Geom. max impact parameter = '',
     &                      1pg15.6)') bmax10
            write(io,*)
            if(ityp.ge.15.and.ityp.le.19.and.iMeVperu.eq.1) then
             write(io,'(''x: E (MeV/n)'')')
            else
             write(io,'(''x: E (MeV)'')')
            endif
            write(io,'(''y: Cross Section (mb)'')')
            if( ielg .eq. 1 ) write(io,'(''p: xlog'')')
            write(io,'(''h: x y(total),l0 y(nonela),d0 y(elast),m0'',
     &                 '' y(Simulation),ur0'')')
            write(io,'(''#   energy     Niita total  Niita elast'',
     &                 ''  Niita inela  Simu.inela'')')

         end if

*-----------------------------------------------------------------------
*     energy mesh and do loop for energy
*-----------------------------------------------------------------------

            if( inucr .ne. 1 .and. inucr .ne. 4 .and.
     &          inucr .ne. 11 ) then

               if( ielg .eq. 0 ) then

                  edif = ( emax - eini ) / dble( nejj - 1 )

               else

                  edif = log( emax / eini ) / dble( nejj - 1 )

               end if

            else

                  edif = 0.0

            end if

*-----------------------------------------------------------------------

      do 5000 iejj = 1, nejj

*-----------------------------------------------------------------------
        if( inucr .eq. 15 )
     &           write(6,*) 'e-mesh=', iejj

*-----------------------------------------------------------------------

         if( inucr .ne. 11 ) then
          if(jstyp(imsrc).eq.4) then ! energy directly specified by user, T.Sato 2022/12/07
           if(iejj.eq.nejj) then ! last energy bin
            ein = emax ! emax has already been multiplied by mass
           else
            ein = egmin(iejj)
            if( ityp.ge.15.and.ityp.le.19.and.iMeVperu.eq.1) then
             ein = ein * jpa ! convert MeV/n to MeV
            end if
           endif
          else
            if( ielg .eq. 0 ) then

               ein = eini + dble( iejj - 1 ) * edif

            else

               ein = eini * exp( dble( iejj - 1 ) * edif )

            end if
          endif
         else

               ein = evle( iejj )

         end if

*-----------------------------------------------------------------------

            if( inucr .eq. 1 .or. inucr .eq. 3 ) then

               write(6,'(''   energy = '',1pg15.6,'' MeV'')') ein

            end if

*-----------------------------------------------------------------------

         if( inucr .eq. 11 ) then

            inumc = inumc + 1

            call date_a_time(iyer1,imon1,iday1,
     &                       ihor1,imin1,isec1)

            if( ityp .eq. 1 )
     &      write(ios,'(i3,'': e ='',f7.1,'' MeV, (z,a)=('',i3,'','',
     &                  i3,'') : proton    '',
     &            i4,''/'',i2.2,''/'',i2.2,3x,
     &            i2.2,'':'',i2.2,'':'',i2.2)')
     &            inumc, ein, ntz, nta,
     &            iyer1, imon1, iday1, ihor1, imin1, isec1

            if( ityp .eq. 2 )
     &      write(ios,'(i3,'': e ='',f7.1,'' MeV, (z,a)=('',i3,'','',
     &                  i3,'') : neutron   '',
     &            i4,''/'',i2.2,''/'',i2.2,3x,
     &            i2.2,'':'',i2.2,'':'',i2.2)')
     &            inumc, ein, ntz, nta,
     &            iyer1, imon1, iday1, ihor1, imin1, isec1

         end if

*-----------------------------------------------------------------------
*        neutron for nuclear data
*-----------------------------------------------------------------------

         if( ityp .eq. 2 .and. ein .le. dnmax(2) ) then
                     mk = mat
                     rh = denm(mat)
                     tme = 0.0

               call xstneu(0,sigt,sigaa,icl,ein,tme,mk)

                  wgti = 1.0d0
                  ffac = 0.0d0

               call sctneu(ein,wgti,ireg,mat,ffac,0,elrt)

               elrt=min(1.0d0,elrt) ! T.Sato 2020/5/31 avoid negative inelastic cross section

                     sigtot = sigt * 1000.0
                     sigela = sigtot * elrt
                     signon = sigtot - sigela


               if( inucr .eq. 12 ) then

                  call heatn(icl,ein,heatr,heatf,mk,rh,0,0,mtdum)
                        heat = ( heatr + heatf ) / rh / sigt

               end if

*-----------------------------------------------------------------------
*        photon for nuclear data
*-----------------------------------------------------------------------

         else if( ityp .eq. 14 .and. ein .le. dnmax(14) ) then

                     mk = mat
                     rh = denm(mat)

               call xstgam(0,sigt,sigaa,ein,mk)

               if( inucr .eq. 2 ) then !S.H. added for photo-nuclear cross section (2020.11.27)
                  if ( ipnint.gt.0 ) then
                     xsc = phxs(ntz,nta-ntz,ein)
                     sigtot = xsc * 1000.0
                     signon = xsc * 1000.0
                     sigela = 0d0
                     sigphoton(1) = getGDRxsec(ntz, nta-ntz, ein) *1d3
                     sigphoton(2) = sigma_qd(ntz, nta-ntz, ein) *1d3
                     sigphoton(3)
     &                    = sigma_delta(dble(nta),dble(ein),1,5) *1d3
                     sigphoton(4) = sigmaNRF(ntz, nta-ntz, ein) *1d3
                     sigphoton(5) = sigphoton(1) + sigphoton(2)
     &                    + sigphoton(3) + sigphoton(4)
                  else
                     sigtot = sigt * 1000.0
                     signon = sigaa * 1000.0
                     sigela = sigtot - signon
                     sigphoton(1) = rtc(1,2) *1d3
                     sigphoton(2) = (rtc(2,2) -  rtc(1,2)) *1d3
                     sigphoton(3) = (rtc(3,2) -  rtc(2,2)) *1d3
                     sigphoton(4) = (rtc(4,2) -  rtc(3,2)) *1d3
                     sigphoton(5) = rtc(4,2) *1d3
                  end if
               end if

               if( inucr .eq. 14 ) then

                  call heatp(icl,ein,heatr,mk,rh,0,0,mtdum)
                        heat =  heatr / rh / sigt

               end if

*-----------------------------------------------------------------------
*        charged particle for nuclear data ! T.Sato 2021/08/14
*-----------------------------------------------------------------------

         else if( (ityp .eq. 1 .and. ein .le. dnmax(1) ) .or.
     &       (ityp .eq. 15 .and. ein/2.0 .le. dnmax(15)) .or.
     &       (ityp .eq. 18 .and. ein/4.0 .le. dnmax(18))) then
                     mk = mat
                     rh = denm(mat)

               call sig_tot(ityp,sigt,sigaa,ein,mk)

                  wgti = 1.0d0
                  ffac = 0.0d0
                  call cp_coll(ein,wgti,ireg,mat,ffac,rh,0,elrt,ityp)

                     sigtot = sigt * 1000.0
                     sigela = sigtot * elrt
                     signon = sigtot - sigela


*-----------------------------------------------------------------------
*        pion for cross section model (S.H. added on 2016.6.22)
*-----------------------------------------------------------------------

         else if( ityp .ge. 3 .and. ityp .le. 5 ) then

          if ( icxspi .eq. 1 ) then ! PHITS original model

            if ( ityp .eq. 3 ) then
               pichg = 1d0
            else if ( ityp .eq. 5 ) then
               pichg = -1d0
            else
               pichg = 0d0
            end if

            at = dble(nta)
            zt = dble(ntz)

            call pionXS(pichg,ein,at,zt,sigtot,signon,sigela,bmax)

            sigtot = sigtot * 1000d0
            signon = signon * 1000d0
            sigela = sigela * 1000d0

          else

            sigtot = siggeo
            signon = 0d0
            sigela = 0d0

          end if

*-----------------------------------------------------------------------
*        hydrogen cross section
*-----------------------------------------------------------------------

         else

            if( inucr .eq. 5 .or. inucr .eq. 6 ) then

                     kf1 = kfft(ityp)
                     kf2 = 2212

                  call sigjam(kf1,kf2,ein,sig,sigel,sigin)

                     sigthyd = sig   * 1000.0
                     sigehyd = sigel * 1000.0
                     sigihyd = sigin * 1000.0

            end if

            if( inucr .eq. 5 ) then

               write(io,'(4(1pe13.5))') ein, sigthyd, sigehyd, sigihyd

            end if

*-----------------------------------------------------------------------
*        Parametrization of Cross sections
*-----------------------------------------------------------------------

            if( ityp .le. 2 ) then

               if( nta .gt. 1 ) then

                     call sigrc(ityp,ein,nta,ntz,sigt,signe,sigel)

               else

                        kf1 = kfft(ityp)
                        kf2 = 2212

                     call sigjam(kf1,kf2,ein,sigt,sigel,signe)

               end if

                        sigtot = sigt  * 1000.0
                        signon = signe * 1000.0
                        sigela = sigel * 1000.0

            else if( ityp .ge. 15 ) then

               call sighi(ap,zp,ein,at,zt,ssigne,ssigel,bmax10)

                        sigtot = ( ssigne + ssigel ) * 1000.0
                        signon = ssigne * 1000.0
                        sigela = ssigel * 1000.0

            else

                        sigtot = siggeo
                        signon = ssigne
                        sigela = ssigel

            end if

         end if

*-----------------------------------------------------------------------
*        for anal-004.f / anal-005.f
*        for anal-006.f / anal-007.f / anal-008.f
*-----------------------------------------------------------------------

         if( inucr .eq. 1 .or. inucr .eq. 11 .or. inucr .eq. 7 .or.
     &       inucr .eq. 8 .or. inucr .eq. 9 .or.
     &       inucr .eq. 12 .or. inucr. eq. 13 .or. inucr .eq. 16) then ! T.Sato 2022/12/06

               mstapr = ntz
               massta = nta

            if( ityp .eq. 1 ) then
               msprpr = 1
            else
               msprpr = 0
            end if
               masspr = 1

         end if

*-----------------------------------------------------------------------
         if(inucr.eq.16) then
          dedxfd=1.0d0
          etmp=ein
          if(etmp/dble(max(1,jpa)).eq.1.0d-3) etmp=etmp*1.00000001d0 ! above just 1 MeV/n
          call dedxas(etmp,dedx,mat,ityp,ktyp,ichgf(ityp,ktyp),
     &    rmtyp(ityp,ktyp))
          if(ityp.ge.15.and.ityp.le.19.and.iMeVperu.eq.1) then
           einp=ein/dble(max(1,jpa))
          else
           einp=ein
          endif
          write(io,'(2(1pe13.5))') einp, dedx*0.1d0 ! convert (MeV/cm) to (keV/um)
         endif

         if( inucr .eq. 2 ) then

          if( ityp .eq. 14 .and. ein .le. dnmax(14) ) then
           if( inucr .eq. 2 ) then

             if ( ipnint.gt.0 ) then
                write(io,'(6(1pe13.5))')
     &               ein, sigphoton(5), (sigphoton(isigph),isigph=1,4)
             else
                write(io,'(6(1pe13.5))')
     &               ein, sigphoton(5), (sigphoton(isigph),isigph=1,4)
             end if

           end if

          else

            if(ityp.ge.15.and.ityp.le.19.and.iMeVperu.eq.1) then
             einp=ein/jpa
            else
             einp=ein
            endif

            write(io,'(4(1pe13.5))') einp, sigtot, signon, sigela

          end if
         end if

*-----------------------------------------------------------------------
*     elastic angular distribution
*-----------------------------------------------------------------------

      if( inucr .eq. 4 ) then

         if( mstz(25) .eq. 0 ) then

               ianm  = mstz(24)
               andif = 180.0 / dble(ianm)

            do ian = 1, ianm

               angds(ian) = 0.0

            end do

            do ian = 1, mevent

               call seldsd(ityp,ein,apr,zpr,csthcm,ian)

               andeg = acos( csthcm ) / pi *180.0
               iangv = int( andeg / andif ) + 1

               if( iangv .ge. 1 .and. iangv .le. ianm ) then

                  angds(iangv) = angds(iangv) + 1.0

               end if

            end do

               write(io,*)
               write(io,'(''c:            energy = '',1pg15.6)') ein
               write(io,'(''c: elastic x-section = '',1pg15.6)') sigela
               write(io,*)
               write(io,'(''h: x y(Pearlstein),d0b '',
     &                    ''y(Pearl + Hut),m0c '',
     &                    ''y(Niita),l0rtt '',
     &                    ''y(Niita Cal.),hhh0r'')')

               sek(0) = 0.0
               sek(1) = 0.0
               sek(2) = 0.0

            do ian = 1, ianm

               aphig = dble(ian-1) * andif + andif
               aplow = dble(ian-1) * andif
               apoin = ( aphig + aplow ) / 2.0
               aprad = apoin * pi / 180.0
               apsin = sin( aprad )
               adrad = ( aphig - aplow ) * pi / 180.0

               angds(ian) = angds(ian) / dble(mevent)
     &                      * sigela / 2.0 / pi / apsin / adrad

               dsigp(ian) = apoin

               icm   = 0

               do ich = 0, 2

                  icct  = ich

                  call dsdarc(icct,icm,ityp,
     &            ein,apr,apoin,signon/1000.,dsigs,escat,etarg,angle)

                  dsiga(ich,ian) = dsigs * 1000.0

                  sek(ich) = sek(ich)
     &                      + dsiga(ich,ian) * 2. * pi * apsin * adrad

               end do

            end do

                  fnor0 = sigela / sek(0)
                  fnor1 = sigela / sek(1)
                  fnor2 = sigela / sek(2)

            do ian = 1, ianm

               write(io,'(1p6e13.5)')
     &           dsigp(ian), dsiga(0,ian)*fnor0, dsiga(1,ian)*fnor1,
     &           dsiga(2,ian)*fnor2, angds(ian)

            end do

*-----------------------------------------------------------------------

         else

               ianm  = mstz(24)
               andif = 180.0 / dble( ianm - 1 )

               write(io,*)
               write(io,'(''c:            energy = '',1pg15.6)') ein
               write(io,'(''c: elastic x-section = '',1pg15.6)') sigela
               write(io,*)
               write(io,'(''h:  n          x        y,l0rtt'')')
               write(io,'(''#  angle      cos       probability'')')

               sek(2) = 0.0

            do ian = 1, ianm

               aphig = min( 180.d0, dble(ian-1) * andif + andif / 2.0 )
               aplow = max(  0.0d0, dble(ian-1) * andif - andif / 2.0 )
               apoin = ( aphig + aplow ) / 2.0
               aprad = apoin * pi / 180.0
               apsin = sin( aprad )
               adrad = ( aphig - aplow ) * pi / 180.0

               dsigp(ian) = apoin

               icm  = 0
               ich  = 2
               icct = 2

               call dsdarc(icct,icm,ityp,
     &         ein,apr,apoin,signon/1000.,dsigs,escat,etarg,angle)

                  dsiga(ich,ian) = dsigs * 1000.0

                  sek(ich) = sek(ich)
     &                      + dsiga(ich,ian) * apsin * adrad

            end do

            sek(0) = 0.0
            sek(1) = 0.0

            do ian = ianm, 1, -1

               aphig = min( 180.d0, dble(ian-1) * andif + andif / 2.0 )
               aplow = max(  0.0d0, dble(ian-1) * andif - andif / 2.0 )
               apoin = ( aphig + aplow ) / 2.0
               aprad = apoin * pi / 180.0
               apsin = sin( aprad )
               adrad = ( aphig - aplow ) * pi / 180.0

               ahcos = cos( aphig * pi / 180.0 )
               alcos = cos( aplow * pi / 180.0 )

                  sek(0) = sek(0)
     &                      + dsiga(ich,ian) / sek(2) * apsin * adrad

                  sek(1) = sek(1)
     &                   + dsiga(ich,ian) / sek(2) * ( alcos - ahcos )

               apoin = dsigp(ian)
               if( ian .eq. ianm ) apoin = 180.0
               if( ian .eq.    1 ) apoin = 0.0

               capon = cos( apoin * pi / 180.0 )
               dsigs = dsiga(ich,ian) / sek(2)

               write(io,'(f7.2,2x,1p2e13.5)') apoin, capon, dsigs

            end do

            write(io,'(/''#normalization :''1p2e13.5)') sek(0), sek(1)

         end if

      end if

*-----------------------------------------------------------------------

         if( inucr .eq. 2 .or.
     &       inucr .eq. 4 .or.
     &       inucr .eq. 5 .or. inucr.eq.16) goto 5000  ! T.Sato 2022/12/16

*-----------------------------------------------------------------------
*     Nuclear Reactions
*-----------------------------------------------------------------------

         if( inucr .eq. 1 .or. inucr .eq. 11 .or. inucr .eq. 7 .or.
     &       inucr .eq. 8 .or. inucr .eq. 9 .or.
     &       inucr .eq. 12 .or. inucr .eq. 13 .or.
     &       inucr .eq. 15) then

            engin  = ein / 1000.0

            if( jpa .ne. 0 ) engin = engin / dble( jpa )

            if( inucr .lt. 12) then
               call anal_int(engin,mevent)
            end if

         end if

*-----------------------------------------------------------------------

               infiss = 0
               inelsn = 0
               inelsm = 0
               inelsl = 0

               sekh = 0.0
               seki = 0.0


               inelst = 0

*-----------------------------------------------------------------------

      do 1000 irunp = 1, mevent

               kelst  = 0

        if(mod(irunp,max(1,int(mevent/100))).eq.0 .and. inucr .ne. 11 )
     &           write(6,*) 'event=', irunp

*-----------------------------------------------------------------------
*        inelastic nuclear reaction
*-----------------------------------------------------------------------

         if( inucr .eq. 1 .or. inucr .eq. 3 .or. inucr .eq. 11 .or.
     &   ( ( inucr .eq. 9 .or. inucr .eq. 12 .or. inucr .eq. 13 .or.
     &       inucr .eq. 15 ) .and.
     &   ( ( ityp .eq. 1 .and. ein .lt. dnmax(1) ) .or.
     &     ( ityp .eq. 2 .and. ein .lt. dnmax(2) ) ) ) ) then

               call cputime(14)

*-----------------------------------------------------------------------
*           Kerma and Photon energy by nucr = 12, 13
*-----------------------------------------------------------------------

            if( ityp .eq. 2 .and. ein .lt. dnmax(2) ) then

                  wgti = 1.0d0
                  ffac = 1.d0

               if( iemode .eq. 0 .or. ein .gt. emodem ) then

                  call sctneu(ein,wgti,ireg,mat,ffac,0,elrt)

                  jcoll = 6

               else

                  call sctneut(ein,wgti,ireg,mat,ffac1,ffac2,0,elrt)
                  sekh = sekh + ffac1
                  seki = seki + ffac2

               end if


*-----------------------------------------------------------------------
*       Charged particle by nuclear data ! T.Sato 2021/08/14
*-----------------------------------------------------------------------

         else if( (ityp .eq. 1 .and. ein .le. dnmax(1) ) .or.
     &       (ityp .eq. 15 .and. ein/2.0 .le. dnmax(15)) .or.
     &       (ityp .eq. 18 .and. ein/4.0 .le. dnmax(18))) then

                  wgti = 1.0d0
                  ffac = 0.0d0
                  rh = denm(mat)

                  call cp_coll(ein,wgti,ireg,mat,ffac,rh,0,elrt,ityp)

                  jcoll = 9

*-----------------------------------------------------------------------

            else

*                  inverse kinematics for proton and light ion targets
                   if( ityp.ge.15 .and. (
     &              ( nta.eq.1 .and. ntz.eq.1 )
     &              .or. ( nta.eq.2 .and. ntz.eq.1 )
     &              .or. ( nta.eq.3 .and. ntz.eq.1 )
     &              .or. ( nta.eq.3 .and. ntz.eq.2 )
     &              .or. ( nta.eq.4 .and. ntz.eq.2 ) ) ) then

                      invpl  = 1
                      iprj = ityp
                      kprj = ktyp

                    if ( nta.eq.1 .and. ntz.eq.1 ) then

                       ityp = 1
                       ktyp = 2212

                    else if ( nta.eq.2 .and. ntz.eq.1 ) then

                       ityp = 15
                       ktyp = 1000002

                    else if ( nta.eq.3 .and. ntz.eq.1 ) then

                       ityp = 16
                       ktyp = 1000003

                    else if ( nta.eq.3 .and. ntz.eq.2 ) then

                       ityp = 17
                       ktyp = 2000003

                    else if ( nta.eq.4 .and. ntz.eq.2 ) then

                       ityp = 18
                       ktyp = 2000004

                    end if

                      nta = ibryf(iprj,kprj)
                      ntz = ichgf(iprj,kprj)
                      rprj = rmtyp(iprj,kprj)

                      rlmass = rstms(ityp)*1000d0
                      ein = rlmass / rprj * ein
                      e2p  = ein + rlmass
                      p2p  = dsqrt( ein**2 + 2.0 * rlmass * ein )
                      beta = p2p / e2p
                      gamm = e2p / rlmass

                   else

                      invpl  = 0

                   end if

                  ipim = 0

   22          continue

                  call ncasc(0,ityp,ktyp,ein,nta,ntz,bmax10)

                  if( nclst .lt. 0 .and.
     &              ( ityp .le. 2 .or. ityp .ge. 15 ) ) then

                     ipim = ipim + 1

                     if( ipim .le. 20 ) goto 22

                  end if


*                  inverse kinematics for proton and light ion targets
                   if( invpl .eq. 1 ) then

                      invpl = 0
                      nta = ibryf(ityp,ktyp)
                      ntz = ichgf(ityp,ktyp)
                      ityp = iprj
                      ktyp = kprj
                      ein = rprj / rlmass * ein

                    if( nclst .gt. 0 ) then

                     do i = 1, nclst

                        px = qclust(1,i)
                        py = qclust(2,i)
                        pz = qclust(3,i)
                        et = qclust(4,i)
                        rm = qclust(5,i)

                        pz = - ( gamm * pz - beta * gamm * et )
                        et = dsqrt( px**2 + py**2 + pz**2 + rm**2 )

                        qclust(3,i) = pz
                        qclust(4,i) = et
                        qclust(7,i) = ( et - rm ) * 1000.d0

                     end do

                    end if

                   end if

            end if

*-----------------------------------------------------------------------

!$OMP CRITICAL (rncnt14_crit)
               rncnt(14) = rncnt(14) + 1.0
!$OMP END CRITICAL (rncnt14_crit)
               call cputime(14)

               if( nclst .lt. 0 ) then

                  inelsn = inelsn + 1

               end if

               if( inucr .eq. 3 ) goto 1000

         end if

*-----------------------------------------------------------------------
*        p-p n-p pi-p reactions
*-----------------------------------------------------------------------

         if( inucr .eq. 6 ) then

               call cputime(8)

                  mmas  = 1
                  mchg  = 1
                  bmax  = 0.0d0

                  call jamin(0,ityp,ktyp,ein,mmas,mchg,bmax)

               rncnt(8) = rncnt(8) + 1.0
               call cputime(8)

            if( nclst .eq. 2 .and.
     &        ( ( ityp .eq. 1 .and. numpat(1) .eq. 2 ) .or.
     &          ( ityp .ne. 1 .and.
     &            numpat(1) .eq. 1 .and. numpat(ityp) .eq. 1 ) ) ) then

                  inelsn = inelsn + 1

            else if( nclst .eq. 3 ) then

                  inelsm = inelsm + 1

            else if( nclst .gt. 3 ) then

                  inelsl = inelsl + 1

            end if

                  goto 1000

         end if

*-----------------------------------------------------------------------
*        elastic + inelastic nuclear reaction
*-----------------------------------------------------------------------

         if( inucr .eq. 7 .or. inucr .eq. 8 .or.
     &   ( ( inucr .eq. 9 .or. inucr .eq. 15 ) .and.
     &   ( ( ityp  .eq. 1 .and. ein .ge. dnmax(1) ) .or.
     &     ( ityp  .eq. 2 .and. ein .ge. dnmax(2) ) ) ) ) then

            if( nta .eq. 1 .and. ntz .eq. 1 ) then

                  call cputime(8)

                  bmax  = 0.0d0

                  call jamin(0,ityp,ktyp,ein,nta,ntz,bmax)

               rncnt(8) = rncnt(8) + 1.0
                  call cputime(8)

               if( nclst .eq. 2 .and.
     &           ( ( ityp .eq. 1 .and. numpat(1) .eq. 2 ) .or.
     &             ( ityp .ne. 1 .and.
     &               numpat(1) .eq. 1 .and.
     &               numpat(ityp) .eq. 1 ) ) ) then

                  kelst = 1

               else

                  kelst = 0

               end if

            else if( unirn(dummy) .le. signon / sigtot ) then

                  call cputime(14)

                  ipim = 0

   23          continue

                  call ncasc(0,ityp,ktyp,ein,nta,ntz,bmax10)

                  if( nclst .lt. 0 .and. ityp .le. 2 ) then

                     ipim = ipim + 1

                     if( ipim .le. 20 ) goto 23

                  end if

                  rncnt(14) = rncnt(14) + 1.0
                  call cputime(14)

                  kelst = 0

            else

                  call cputime(13)

                  call nelst(ityp,ktyp,apr,zpr,ein)

                  rncnt(13) = rncnt(13) + 1.0
                  call cputime(13)

                  kelst = 1

            end if

         end if

*-----------------------------------------------------------------------
*     after ncasc
*-----------------------------------------------------------------------

      if( nclst .gt. 0 .and. inucr .lt. 12 ) then

         do 100 i = 1, nclst

                  ik = iclust(i)

                  jj = jclust(0,i)
                  iz = jclust(1,i)
                  in = jclust(2,i)
                  id = jclust(3,i)
                  is = jclust(4,i)
                  ic = jclust(5,i)
                  iq = jclust(6,i)
                  im = jclust(7,i)

                  bi = qclust(0,i)
                  px = qclust(1,i)
                  py = qclust(2,i)
                  pz = qclust(3,i)
                  et = qclust(4,i)
                  rm = qclust(5,i)
                  ex = qclust(6,i)
                  ek = qclust(7,i)
                  we = qclust(8,i)

               call anal_qmd(ik,jj,iz,in,id,is,ic,iq,im,
     &                       bi,px,py,pz,et,rm,ex,ek,we)

  100    continue

      end if

*-----------------------------------------------------------------------
*     evaporation and fission
*-----------------------------------------------------------------------

            if( ityp .eq. 2 .and. ein .le. dnmax(2) ) then

            else if( ityp .eq. 1 .and. ein .le. dnmax(1) ) then

            else

               call nevap(0)

            end if

*-----------------------------------------------------------------------
*     after nevap
*-----------------------------------------------------------------------

               if( kdecay(4) .eq. 1 ) then

                  infiss = infiss + 1

               end if

      if( nclsts .gt. 0 .and.
     &  ( inucr .lt. 12 .or. inucr .eq. 15 )  ) then

            itrstar = 0
            nmpro = 0
            nmnut = 0

         do i = 1, nclsts

                  ik = iclusts(i)

                  jj = jclusts(0,i)
                  iz = jclusts(1,i)
                  in = jclusts(2,i)
                  id = jclusts(3,i)
                  is = jclusts(4,i)
                  ic = jclusts(5,i)
                  iq = jclusts(6,i)
                  im = jclusts(7,i)

                  bi = qclusts(0,i)
                  px = qclusts(1,i)
                  py = qclusts(2,i)
                  pz = qclusts(3,i)
                  et = qclusts(4,i)
                  rm = qclusts(5,i)
                  ex = qclusts(6,i)
                  ek = qclusts(7,i)
                  we = qclusts(8,i)

            if( inucr .lt. 12)
     &         call anal_sdm(ik,jj,iz,in,id,is,ic,iq,im,
     &                       bi,px,py,pz,et,rm,ex,ek,we)

*-----------------------------------------------------------------------

            if( kelst .eq. 0 .and.
     &          ntz .eq. jclusts(1,i) .and.
     &          nta-ntz .eq. jclusts(2,i) )  itrstar = 1

               if( iclusts(i) .eq. 1 ) nmpro = nmpro + 1
               if( iclusts(i) .eq. 2 ) nmnut = nmnut + 1

         end do

*-----------------------------------------------------------------------

         if( itrstar .eq. 0 ) then

             if( jcoll .eq. 6 .and.
     &           kcoll .ne. 3 .and.
     &         ( nmpro .eq. 0 .and. nmnut .eq. 1 ) ) then
                 itrstar = 1
             end if

             if( jcoll .eq. 9 .and.
     &           kcoll .ne. 3 .and.
     &         ( nmpro .eq. 1 .and. nmnut .eq. 0 ) ) then
                 itrstar = 1
             end if

         end if

         if( itrstar .eq. 1 ) inelst = inelst + 1

*-----------------------------------------------------------------------

      end if

*-----------------------------------------------------------------------

 1000 continue

*-----------------------------------------------------------------------

         if( inucr .eq. 15 ) then

            sigie = sigtot * dble(inelst) / dble(mevent)

            if(ityp.ge.15.and.ityp.le.19.and.iMeVperu.eq.1) then
             einp=ein/jpa
            else
             einp=ein
            endif

            write(io,'(5(1pe13.5))') einp, sigtot, signon, sigela, sigie

         end if

*-----------------------------------------------------------------------

         if( inucr .eq. 1 .or. inucr .eq. 3 .or. inucr .eq. 11 ) then

            if( ityp .le. 2 ) then

               siginp = signon

            else

               siginp = siggeo

            end if

            reaxrs = siginp * dble(mevent-inelsn) / dble(mevent)
            bmax0  = sqrt( siginp / 10.0 / pi )
            bmax1  = sqrt( reaxrs / 10.0 / pi )

            fissx  = siginp * dble(infiss) / dble(mevent)

         end if

*-----------------------------------------------------------------------

         if( inucr .eq. 1 ) then

            write(6,*)
            write(6,'('' initial inelastic x-section  = '',1pg15.6)')
     &                siginp
            write(6,'('' initial max impact parameter = '',1pg15.6)')
     &                bmax0
            write(6,*)

            write(6,'(''   final inelastic x-section  = '',1pg15.6)')
     &                reaxrs
            write(6,'(''   final max impact parameter = '',1pg15.6)')
     &                bmax1
            write(6,*)

            write(6,'(''   Niita inelastic x-section  = '',1pg15.6)')
     &                signon

            write(6,*)
            write(6,'(''           fission x-section  = '',1pg15.6)')
     &                fissx

         end if

*-----------------------------------------------------------------------

         if( inucr .eq. 7 .or. inucr .eq. 8 .or. inucr .eq. 9 ) then

            bmax0  = sqrt( sigtot / 10.0 / pi )

         end if

*-----------------------------------------------------------------------

         if( inucr .eq. 1 .or. inucr .eq. 11 .or.
     &       inucr .eq. 7 .or. inucr .eq. 8 ) then

            call anal_fin(io,bmax0)

         end if

*-----------------------------------------------------------------------
         if( inucr .eq. 9 ) then
              app = 0.0
              zpp = 0.0
              app = jpa
              zpp = jpz

            if(ktyp .eq. 2212) then
              app = 1.0
              zpp = 1.0
            else if(ktyp .eq. 2112) then
              app = 1.0
              zpp = 0.0
            endif

            call anal_fin2(io,bmax0,itdpa,iteth,zpp,app)

         end if

*-----------------------------------------------------------------------

         if( inucr .eq. 3 ) then

            if(ityp.ge.15.and.ityp.le.19.and.iMeVperu.eq.1) then
             einp=ein/jpa
            else
             einp=ein
            endif

            write(io,'(5(1pe13.5))') einp, sigtot, sigela, signon,
     &                               reaxrs

         end if

*-----------------------------------------------------------------------

         if( inucr. eq. 6 ) then

            siginp = sigthyd * dble(mevent-inelsn) / dble(mevent)
            siginn = dble(inelsn) / dble(mevent) * 100.0
            siginm = dble(inelsm) / dble(mevent) * 100.0
            siginl = dble(inelsl) / dble(mevent) * 100.0

            if( ityp .ge. 3 ) then

               em1  = rmtyp(ityp,kf)
               em2  = rmtyp(1,kf)
               elab = ein

               pproj = sqrt(elab*(2.0*em1+elab))
               srt   = sqrt((elab+em1+em2)**2-pproj**2)

               write(io,'(4(1pe13.5),3(1pg13.5))')
     &                                  srt, sigthyd, sigihyd, siginp,
     &                                  siginn, siginm, siginl

            else

            if(ityp.ge.15.and.ityp.le.19.and.iMeVperu.eq.1) then
             einp=ein/jpa
            else
             einp=ein
            endif

               write(io,'(4(1pe13.5),3(1pg13.5))')
     &                                  einp, sigthyd, sigihyd, siginp,
     &                                  siginn, siginm, siginl

            end if

         end if

*-----------------------------------------------------------------------

      if( inucr. eq. 12 .or. inucr .eq. 13 ) then

*-----------------------------------------------------------------------
*        neutron for nuclear data
*-----------------------------------------------------------------------
         if( ityp .eq. 2 .and. ein .le. dnmax(2) ) then

                     mk = mat
                     rh = denm(mat)
                     tme = 0.0

               call xstneu(0,sigt,sigaa,icl,ein,tme,mk)

                  heat1 = sekh / dble( mevent )
                  heat2 = seki / dble( mevent )

            if( inucr. eq. 12 ) then
               write(io,'(5e13.4)') ein, heat1, heat,
     &                              heat1*sigt, heat*sigt
            else
               write(io,'(3e13.4)') ein, heat1, heat2
            end if

         end if

      end if

*-----------------------------------------------------------------------

      if( inucr. eq. 14 ) then

*-----------------------------------------------------------------------
*        photon for nuclear data
*-----------------------------------------------------------------------

         if( ityp .eq. 14 .and. ein .le. dnmax(14) ) then

                     mk = mat
                     rh = denm(mat)

               call xstgam(0,sigt,sigaa,ein,mk)

               write(io,'(5e13.4)') ein, heat, heat*sigt, sigt

         end if

      end if

*-----------------------------------------------------------------------

 5000 continue

*-----------------------------------------------------------------------

 5100 continue
 5200 continue

      call DEALLOCATE_MMBANK    !FURUTA
      call DEALLOCATE_MEMBANK   !FURUTA

      if(icgg.ne.0) call DEALLOCATE_GGBANK !FURUTA
      call DEALLOCATE_GGMBANK   !FURUTA
      call DEALLOCATE_EVTS      !FURUTA20210506
      return

  999 continue

*-----------------------------------------------------------------------

      return
      end


