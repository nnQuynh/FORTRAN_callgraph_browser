************************************************************************
*                                                                      *
      subroutine nreac
*                                                                      *
*        control of intranuclear reaction part                         *
*        last modified by K.Niita on 2011/05/17                        *
*                                                                      *
*---- in common -------------------------------------------------------*
*                                                                      *
*        nclst   : total number of out going particles and nuclei      *
*        iclust(nclst)                                                 *
*        jclust(i,nclst)                                               *
*        qclust(i,nclst)                                               *
*        numpat(i) : total number of out going particles or nuclei     *
*        mathz, mathn : z and n of mather nucleus                      *
*                                                                      *
*---- in common for nuclear data process ------------------------------*
*                                                                      *
*        nclsts  : total number of out going particles and nuclei      *
*        iclusts(nclsts)                                               *
*        jclusts(i,nclsts)                                             *
*        qclusts(i,nclsts)                                             *
*                                                                      *
*----------------------------------------------------------------------*
*        jcoll reaction type identifier                                *
*----------------------------------------------------------------------*
*                                                                      *
*        jcoll : =  0, nothing happen                                  *
*                =  1, Hydrogen collisions                             *
*                =  2, Particle Decays                                 *
*                =  3, Elastic collisions                              *
*                =  4, High Energy Nuclear collisions                  *
*                =  5, Heavy Ion reactions                             *
*                =  6, Neutron reactions by data                       *
*                =  7, Photon reactions by data                        *
*                =  8, Electron reactions by data                      *
*                =  9, P,d,a, and photo-nuclear reactions by data      *
*                = 10, Neutron event mode                              *
*                = 11, Delta Ray production                            *
*                = 12, Muon atomic interaction                         *
*                = 13, Photon by EGS5                                  *
*                = 14, Electron by EGS5                                *
*                = 15, Photon photonuclear interaction                 *
*                = 16, Negative muon captured by nucleon               *
*                = 17, Muon photonuclear interaction                   *
*                = 18, Electron recoil by track strcuture mode         *
*                = 19, Muon pair production (photon -> mu+ mu-)        *
*                = 20, User defined interaction                        *
*                                                                      *
*----------------------------------------------------------------------*
*        kcoll reaction type identifier                                *
*----------------------------------------------------------------------*
*                                                                      *
*        kcoll : =  0, normal                                          *
*                =  1, high energy fission                             *
*                =  2, high energy absorption                          *
*                =  3, low energy n elastic                            *
*                =  4, low energy n non-elastic                        *
*                =  5, low energy n fission                            *
*                =  6, low energy n absorption                         *
*                                                                      *
*----------------------------------------------------------------------*
*        itcoll Track structrue identifier                              *
*----------------------------------------------------------------------*
*                                                                      *
*       itcoll : =  1, electron track-structure by Kai code            *
*                =  2, Proton track structure by Kurbuc                *
*                =  3, Carbon ion track structure by Kurbuc            *
*                =  4, Generalized target-ion track structure mode     *
*                                                                      *
*----------------------------------------------------------------------*
*        regionwise and medium count of collisions and decays          *
*----------------------------------------------------------------------*
*                                                                      *
*            n : =  0, nothing happen                                  *
*                =  1, Hydrogen collisions                             *
*                =  2, Particle Decays                                 *
*                =  3, Elastic collisions                              *
*                =  4, Nuclear collisions                              *
*                                                                      *
*         ihvi : =  5, deuteron induced                                *
*                =  6, tritium induced                                 *
*                =  7, 3He induced                                     *
*                =  8, Alpha induced                                   *
*                =  9, Nucleus induced                                 *
*                                                                      *
*            n : = 10, neutron by nuclear data                         *
*                = 11, photon by nuclear data                          *
*                = 12, electron by nuclear data                        *
*                = 57, proton by nuclear data                          *
*                = 63, neutron event mode                              *
*                                                                      *
*----------------------------------------------------------------------*
*        in neutron collision by nuclear data                          *
*----------------------------------------------------------------------*
*                                                                      *
*            n : = 13, n-photon product weight                         *
*                = 14, n-photon product energy                         *
*                                                                      *
*                = 15, n-capture loss weight                           *
*                = 16, n-capture loss energy                           *
*                                                                      *
*                = 17, elastic event weight                            *
*                                                                      *
*                = 18, (N,xN) event weight                             *
*                = 19, fission event weight                            *
*                                                                      *
*                = 20, (N,xN) neutron creation weight                  *
*                = 21, fission neutron creation weight                 *
*                                                                      *
*                = 22, upscattering neutron gain energy                *
*                = 23, downscattering neutron loss energy              *
*                = 24, (N,N') event weight                             *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*                = 25, High Energy Fission                             *
*                                                                      *
*----------------------------------------------------------------------*
*        in photon collision by nuclear data                           *
*----------------------------------------------------------------------*
*                                                                      *
*            n : = 26, b-photon product weight                         *
*                = 27, b-photon product energy                         *
*                                                                      *
*                = 28, p-capture loss weight                           *
*                = 29, p-capture loss energy                           *
*                                                                      *
*                = 30, p-annihilation product weight                   *
*                = 31, p-annihilation product energy                   *
*                                                                      *
*                = 32, electron x-ray product weight                   *
*                = 33, electron x-ray product energy                   *
*                                                                      *
*                = 34, fluorescence product weight                     *
*                = 35, fluorescence product energy                     *
*                                                                      *
*                = 36, pair product loss weight                        *
*                = 37, pair product loss energy                        *
*                                                                      *
*                = 62, compton loss energy                             *
*                                                                      *
*----------------------------------------------------------------------*
*        in electron collision by nuclear data                         *
*----------------------------------------------------------------------*
*                                                                      *
*            n : = 38, pair product weight                             *
*                = 39, pair product energy                             *
*                                                                      *
*                = 40, compton recoil weight                           *
*                = 41, compton recoil energy                           *
*                                                                      *
*                = 42, photo-electric weight                           *
*                = 43, photo-electric energy                           *
*                                                                      *
*                = 44, photo and electron auger weight                 *
*                = 45, photo and electron auger energy                 *
*                                                                      *
*                = 46, knock-on weight                                 *
*                = 47  knock-on energy                                 *
*                                                                      *
*                = 48  scattering loss energy                          *
*                = 49, brems loss energy                               *
*                                                                      *
*----------------------------------------------------------------------*
*        high energy library                                           *
*----------------------------------------------------------------------*
*                                                                      *
*                = 50, (n,xp) event weight                             *
*                = 51, (n,pion) event weight                           *
*                = 52, (n,other) event weight                          *
*                                                                      *
*                = 59, p-capture event weight                          *
*                = 60, p-elastic event weight                          *
*                = 61, p-nonelastic event weight                       *
*                                                                      *
*                = 53, (p,xp) event weight                             *
*                = 54, (p,xn) event weight                             *
*                = 58, (p,photon) event weight                         *
*                = 55, (p,pion) event weight                           *
*                = 56, (p,other) event weight                          *
*                                                                      *
************************************************************************
C for USE_MOD_COUNTER
      use mod_counter, only: rncnt,rnint,rnintr,rnpnt,rnpntr
     &                      ,aevts,aevtr,bevts,bevtr
     &                      ,iaevt, ibevt, jaevt, jbevt
!$   &                      ,rncnt2,rnint2,rnintr2,rnpnt2,rnpntr2
!$   &                      ,aevts2,aevtr2,bevts2,bevtr2

      use MMBANKMOD !FURUTA
      use MEMBANKMOD !FURUTA
      use neutrino_mod, only : neutrino_kinem1, neutrino_kinem2
      use ion_track_structure, only : itsreac
      use ets_art, only : etsreac_art
      use ELEDATAMOD, only : ichem

!<-20211126murofushi add
      use moddas_material
      use udm_Parameter
      use udm_Utility
      use udm_Manager

      implicit real*8(a-h,o-z)

*-----------------------------------------------------------------------

      include 'param00.inc'
      include 'param.inc'
      include 'param-physcnst.inc'

      parameter ( rpmass = 938.27, rnmass = 939.58 )

*-----------------------------------------------------------------------

      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
!$OMP THREADPRIVATE(/icomon/)
      common /tlgeom/ iblz1,iblz2
!$OMP THREADPRIVATE(/tlgeom/)
      common /tlglat/ ilev1,ilev2,ilat1(5,10),ilat2(5,10)
!$OMP THREADPRIVATE(/tlglat/)
      common /mathzn/ mathz, mathn, jcoll, kcoll
!$OMP THREADPRIVATE(/mathzn/)
      common /preeq/  npcle, nhole, efermi, atar, ztar
!$OMP THREADPRIVATE(/preeq/)
      common /poabs/  iabsms
!$OMP THREADPRIVATE(/poabs/)
      common /ndemax/ dnmax(20)
      common /wtsave/ oldwt
!$OMP THREADPRIVATE(/wtsave/)
      common /eparm/  esmax, esmin, emin(20)

      common /clustv/ kdecay(4)
!$OMP THREADPRIVATE(/clustv/)

      common /delesg/ rdels(kvlmax), mndel, ndels(kvlmax)
      common /delreg/ delm(kvlmax), kdelt


      common /anglth/ cosphi,costh,sinphi,sinth
!$OMP THREADPRIVATE(/anglth/)

      common /kmat1d/ idmn(0:kvlmax), idnm(kvmmax)
      common /xtotal/ totttl,tothyd,totdcy,totdel,tottsm
!$OMP THREADPRIVATE(/xtotal/)
      common /geosig/ geosig(250)
      common /regdc/  idrg(kvlmax), idgr(kvmmax)

      common /fgdata/ ifrgd, ifgm, ifgdf(kvlmax,5), frgfl(kvlmax)
      character frgfl*200


      common /kmat1g/ kmat(kvlmax)
      common /xinels/ ksige, ksign

      common /kmat1k/ kmatd(kvlmax), kmate(kvlmax)


*-----------------------------------------------------------------------

      common /clustf/ nclst, iclust(nnn)
!$OMP THREADPRIVATE(/clustf/)
      common /clustg/ jclust(0:8,nnn), qclust(0:12,nnn)
!$OMP THREADPRIVATE(/clustg/)
      common /clustp/ rumpat(0:20), numpat(0:20)
!$OMP THREADPRIVATE(/clustp/)
      common /tcntl/  icntl, inucr
      common /kmatad/ matadd
      common /celdn/  denr(kvlmax), denm(kvlmax), denc(kvlmax)

      common /clustt/ nclsts, iclusts(nnn)
!$OMP THREADPRIVATE(/clustt/)
      common /clustw/ jclusts(0:8,nnn),  qclusts(0:12,nnn)
!$OMP THREADPRIVATE(/clustw/)

      common /elpos/ ielpos
!$OMP THREADPRIVATE(/elpos/)
      common /emode/  emodem, ge1, ge2, iemode

      common /cdwba/  bdwba, idwbahit
!$OMP THREADPRIVATE(/cdwba/)

      dimension xv(3)

      dimension data(11)

*-----------------------------------------------------------------------
      common /muint/ imuint, imubrm, imuppd, imucap
      common /muflag/ imuinthit,imubrmhit,imuppdhit,imucaphit,imucapflag
!$OMP THREADPRIVATE(/muflag/)
      common /sigmuh/ sigmuv(200),sigmub(200),sigmup(200),
     &                nmuh,itzmuh(200),itamuh(200)
!$OMP THREADPRIVATE(/sigmuh/)

*-----------------------------------------------------------------------

      integer iegsemi, iegsout
      common /egsemi/ iegsemi, iegsout

      common /pnint/  ipnint
      common /gmuppd/ igmuppd ! S.Abe 2019/11/07
      common /xphotn/ totgam,totegs,totpni,totgmm ! y.sakaki 2019/09/26
!$OMP THREADPRIVATE(/xphotn/)

      common /impactcommon/ oldimpact
!$OMP THREADPRIVATE(/impactcommon/)

      common /pionflag/ iqdflag, ipionflag, istrflag, iphreturn
!$OMP THREADPRIVATE(/pionflag/)

      common /pnmul/ pnimul
      common /gmmul/ gmumul ! S.Abe 2019/11/07

      common /tstara/ atmrc(10,8)      ! Ogawa 2023/7/18, (10,7) -> (10,8) for plasmon
!$OMP THREADPRIVATE(/tstara/)
      common /tpdcta/ nflumu
!$OMP THREADPRIVATE(/tpdcta/)

      common /  tscmsg   / ktsc(kvlmax), mntsc, ntsc(kvlmax)
      common /  tscreg   / ntscell(kvlmax)
      common / etsminmax / etsmin, etsmax
      common / ptsminmax / ptsmin, ptsmax
      common / ctsminmax / ctsmin, ctsmax

      common / tsminmax  / tsmax ! T.Ogawa 2020/06/15

      common /errorcom/ncascerr ! T.Sato 2018/06/26
!$OMP THREADPRIVATE(/errorcom/)

      common /cascid/ jcasc
!$OMP THREADPRIVATE(/cascid/)

      common /scinswt/iswitch
!$OMP THREADPRIVATE(/scinswt/)

      common /dpnmaxcom/ dpnmax ! photon nuclear library maximum energy

      common /kmat1ka/ kmathd(kvlmax), kmathe(kvlmax)

*-----------------------------------------------------------------------

            call cputime(6)

*-----------------------------------------------------------------------
*        zero set
*-----------------------------------------------------------------------

                  invp  = 0
                  invpl = 0
                  ihvi  = 0

                  npcle = 0
                  nhole = 0
                  atar  = 0.0
                  ztar  = 0.0

                  mathz = 0
                  mathn = 0

                  jcoll =  0
                  kcoll =  0
                  itcoll=  0
                  nclst = -1
                 nclsts =  0
                 nflumu =  0    ! S.Abe 2017/02/08

               do i = 0, 20

                  numpat(i) = 0
                  rumpat(i) = 0

               end do

               do i = 1, 4
                  kdecay(i) = 0
               end do

                  icl = idgr(iblz1)

               iqdflag = 0
               ipionflag = 0
               istrflag = 0
               iphreturn = 0
               imuinthit = 0
               imubrmhit = 0
               imuppdhit = 0
               imucaphit = 0
               imucapflag = 0

               oldimpact = 0.0d0

! T.Sato 2018/06/27, ncascerr is ocasionally not initialized by previous reaction
               ncascerr = 0

               atmrc  = 0.d0
               nflumu = 0

*-----------------------------------------------------------------------

               if( iabsms .eq. 1 ) then

                  totdcy = 0.0d0

                  if( tothyd + totttl + totdel + tottsm .le. 0.0d0 )
     &                totdcy = 1.0d0


                  if( (ityp .eq. 6 .or. ityp .eq. 7) .and.
     &                imuint .eq. 1 ) then
                      totdcy = 1.0d0
                      totttl = 0.0d0
                      tothyd = 0.0d0
                      totdel = 0.0d0
                      tottsm = 0.0d0
                  endif

               else if( iabsms .eq. 2 ) then

                  totdcy = 1.0d0
                  tothyd = 0.0d0
                  totttl = 0.0d0
                  totdel = 0.0d0
                  tottsm = 0.0d0

               end if

*-----------------------------------------------------------------------
*        region and material
*-----------------------------------------------------------------------

                  ireg = idgr( iblz1 )

               if( mat .le. kvlmax .and. mat .gt. 0 ) then

                  if( matadd .ne. 0 ) then

                     imat = mat

                  else

                     imat = idnm( idmn(mat) )

                  end if

               else

                     imat = mat

               end if

*-----------------------------------------------------------------------
*        for test or nuclear reactions
*-----------------------------------------------------------------------

               if( icntl .eq. 1 ) then

                  lem = 1

                  goto 13

               end if

*-----------------------------------------------------------------------
* (Takeshi Kai 2017/05/31) for electron track structure mode
* (Takeshi Kai 2018/12/10) for ion track structure mode
*-----------------------------------------------------------------------
               mat1 = 0
               do i=1,mntsc
                if(ntsc(i).eq.iblz(ibkblz+no,ipomp+1))then
                 mat1 = ntscell(idgr(ntsc(i)))
                 exit
                endif
               enddo

*-----------------------------------------------------------------------
*     total flight length
*     random number for the reaction channels
*-----------------------------------------------------------------------

               delsig = totttl + tothyd + totdcy + totdel + tottsm

                  if( delsig .le. 0.0d0 ) goto 1000

 2000 continue   ! S.Abe 2021/06/03, avoid to go to ncasc

               r = unirn(dummy)

*-----------------------------------------------------------------------
*     incident direction of particle
*-----------------------------------------------------------------------

                  costh = w(ibkw+no,ipomp+1)
                  rt2 = u(ibku+no,ipomp+1)**2 + v(ibkv+no,ipomp+1)**2

               if( rt2 .eq. 0.0d0 ) then

                  sinth = 0.0d0
                  cosphi= 1.0d0
                  sinphi= 0.0d0

               else

                  rt = sqrt(rt2)
                  sinth  = rt
                  cosphi = u(ibku+no,ipomp+1) / rt
                  sinphi = v(ibkv+no,ipomp+1) / rt

               end if

*-----------------------------------------------------------------------
*     Particle Decay
*-----------------------------------------------------------------------

               r = r - totdcy / delsig

         if( r .le. 0.0d0 .or. totdcy .gt. 1.0d+09 ) then

            if( imat .ne. 0 .and. ityp .eq. 7 .and. imucap .gt. 0 ) then

               call capmujudge(ireg,imat,imucapflag,itz,ita)

               wgti  = wt(ibkwt+no,ipomp+1)
               call capmuxray(imucapflag,wgti,itz,ita,emubind)

            endif

            if( imucapflag .eq. 1 ) then

               call cputime(26)

               jcoll = 16   ! S.Abe 2017/01/17
               call capmu(wgti,itz,ita)

               rncnt(26) = rncnt(26) + 1.0
               call cputime(26)

            else

               call cputime(7)

               jcoll = 2

               iprj = ityp
               kprj = ktyp
               eein = ec(ibkec+no,ipomp+1)

               call dklos(iprj,kprj,eein)

               rncnt(7) = rncnt(7) + 1.0
               call cputime(7)

            endif

            goto 1000

         end if

*-----------------------------------------------------------------------
*     for void
*-----------------------------------------------------------------------

         if( imat .le. 0 ) goto 1000

*-----------------------------------------------------------------------
*     Hydrogen Collision calculated by JAM
*-----------------------------------------------------------------------

               r = r - tothyd / delsig

         if( r .le. 0.0d0 ) then

                  jcoll = 1
                  ipim  = 0

   21          continue

*-----------------------------------------------------------------------
*           particle - hydrogen
*-----------------------------------------------------------------------

            if( ityp .le. 14 ) then

               call cputime(8)

                  mmas  = 1
                  mchg  = 1
                  mathz = 1
                  mathn = 0
                  iprj  = ityp
                  kprj  = ktyp
                  bmax  = 0.0d0

                  eein  = ec(ibkec+no,ipomp+1)
                  if( iabsms .eq. 1 .and. eein .lt. emin(ityp) )
     &            eein = emin(ityp)


! T.Sato 2016/03/08 determine elastic-inelastic ratio
           if(kprj.eq.2212.or.kprj.eq.2112) then ! neutron or proton
            kf2 = 2212
            call sigjam(kprj,kf2,eein,sig,sigel,signo)
            tmp = unirn(dummy)
            if(tmp.lt.sigel/sig) then ! should be elastic
             ielastic=1
            else  ! should be non-elastic
             ielastic=0
            endif
            do itatus=1,20 ! try 20 time, if 20 time is not enough, just use the last reaction data
             jcasc = 3      ! S.Abe 2019/03/13, JAM flag for therdec
             call jamin(0,iprj,kprj,eein,mmas,mchg,bmax)
             icheck=0 ! basically inelastic
             if(nclst.eq.2) then ! only two particle emitted, possibility of elastic
              if(kprj.eq.2212) then ! proton incidence
               if(jclust(7,1).eq.2212.and.jclust(7,2).eq.2212) icheck=1 ! p-p elastic
              else ! neutron incidence
               if(jclust(7,1).eq.2212.and.jclust(7,2).eq.2112) icheck=1 ! p-n elastic
               if(jclust(7,1).eq.2112.and.jclust(7,2).eq.2212) icheck=1 ! p-n elastic
              endif
             endif
             if(nclst.ne.-1.and.ielastic.eq.icheck) exit ! reaction check is OK
            enddo
            if(icheck.eq.1) jcoll=3

c S.H. determine elastic-inelastic ratio for pion+,0,- (2016.7.8)
           else if( iprj .ge. 3 .and. iprj .le. 5 ) then
            if ( iprj .eq. 3 ) then
             pichg = 1d0
            else if ( iprj .eq. 5 ) then
             pichg = -1d0
            else
             pichg = 0d0
            end if
            call pionXS(pichg,eein,1d0,1d0,sigt,signe,sigel,bmax)
            tmp = unirn(dummy)
            if(sigt.le.0d0) then ! to avoid division by zero (2017.2.1)
               sigt = 1d-6 ! 0.001 mb
               sigel = sigt
            end if
            if(tmp.lt.sigel/sigt) then ! should be elastic
             ielastic=1
            else  ! should be non-elastic
             ielastic=0
            endif
            do itatus=1,20 ! try 20 time, if 20 time is not enough, just use the last reaction data
             jcasc = 3      ! S.Abe 2019/03/13, JAM flag for therdec
             call jamin(0,iprj,kprj,eein,mmas,mchg,bmax)
             icheck=0 ! basically inelastic
             if(nclst.eq.2) then ! only two particle emitted, possibility of elastic
              if ( ( jclust(3,1).eq.iprj .and. jclust(3,2).eq.1 ) .or.
     &             ( jclust(3,1).eq.1 .and. jclust(3,2).eq.iprj ) ) then
               icheck=1 ! elastic
              end if
             end if
             if(nclst.ne.-1.and.ielastic.eq.icheck) exit ! reaction check is OK
            enddo
            if(icheck.eq.1) jcoll=3

           elseif (abs(ktyp) .eq. 12 .or. abs(ktyp) .eq. 14 .or.
     &      abs(ktyp) .eq. 16 ) then ! neutrino-electron scattering or neutrino-proton scattering
              call neutrino_kinem1(ktyp,1,1,ec(ibkec+no,ipomp+1))
           else ! other particle, do not have to consider elasitc/inelasitc ratio
                  jcasc = 3      ! S.Abe 2019/03/13, JAM flag for therdec
                  call jamin(0,iprj,kprj,eein,mmas,mchg,bmax)

           endif
           if(jcoll.eq.3) then
            rncnt(13) = rncnt(13) + 1.0
           else
            rncnt(8) = rncnt(8) + 1.0
           endif

               call cputime(8)

*-----------------------------------------------------------------------
*              absorption of negative charge mesons
*-----------------------------------------------------------------------

               if( nclst .lt. 0 .and. iabsms .eq. 1 ) then

                  ipim = ipim + 1

                  if( ipim .le. 10 ) goto 21

               else if( nclst .gt. 0 .and. iabsms .eq. 1 ) then

                  ec(ibkec+no,ipomp+1) = emin(ityp)

               end if

                  goto 1000

*-----------------------------------------------------------------------
*           nucleus - hydrogen
*-----------------------------------------------------------------------

            else if( ityp .ge. 15 ) then

                  invp  = 1
                  ihvi  = ityp - 10

                  jcoll = 4

                  goto 13

            end if

*-----------------------------------------------------------------------

         end if

*-----------------------------------------------------------------------
*     neutron for nuclear data
*-----------------------------------------------------------------------

                  eein  = ec(ibkec+no,ipomp+1)

*-----------------------------------------------------------------------
*     user defined interaction
*-----------------------------------------------------------------------

      if( udm_int_num .gt. 0) then

        if( ityp .ne. 14 ) then
C         --------------------------------------------------------------
          if(abs(ktyp)==11) then
C           udm_p1  !        Probability for the UDM Process(es)
C           udm_q1  ! Biased Probability for the UDM Process(es)
C           udm_p2  !        Probability for the normal Process
C           udm_q2  ! Biased Probability for the normal Process
C           udm_w1i ! Weight Factor for the UDM Process i
C           udm_w2  ! Weight Factor for the Normal Process
C           -----
            udm_p1 = totudm*udm_fpl ! totudm[1/cm] * udm_fpl[cm]
            udm_p2 = 1d0-udm_p1
            udm_tmp = totudm_mul*udm_fpl + udm_p2
            udm_q1  = totudm_mul*udm_fpl / udm_tmp
            udm_q2  = udm_p2             / udm_tmp
            udm_w2  = udm_p2 / udm_q2
            if(udm_p1>=1d0) udm_w2=0d0
          else
            totmul = (totttl-totudm) + totudm_mul
            udm_q1 = totudm_mul / totmul
            udm_w2 = totmul / delsig ! Weight Factor for the Normal Process
          endif
          ! ------------------------------------------------------------
          if(rn(0) .lt. udm_q1) then
          ! ------------------------------------------------------------
          ! Choose a process from user-defined processes ( totudm_mul_(i) )
          ! Calculate weight for the choosed process
          tmp=0d0
          r_udp=rn(0)
          call fill_mat_info(imat)
          do i = 1, udm_int_num
            tmp=tmp + (totudm_(i) * udm_bias(i)) / totudm_mul
            udm_Kin=eein
            udm_kf_incident=ktyp
            udm_sigt=0d0
            call user_defined_interaction(11,i) ! Calculate: 'udm_sigt'
            if(udm_sigt.eq.0d0) cycle
            ! If the charged particle is an incident particle, the energy may have
            ! changed compared to that used to calculate totudm_(i) in getflt.f.
            ! Processes with a cross section of 0 for the changed energy are eliminated.
            ! For example, this can be a problem for production processes in narrow
            ! resonance processes. This could be addressed by
            ! (1) Shortening the following step lengths
            !     chard (Default value = 0.1. For electrons and positrons)
            !     deltm (Default value = 20.12345. For others)
            ! (2) Comment out 'if(udm_sigt==0d0) cycle' and address the problem in the udm_int* file.
            if(r_udp .lt. tmp) then ! index=i is accepted.
              jcoll = 20
C             ----------
              if(abs(ktyp)==11) then
                udm_w1i = udm_tmp / udm_bias(i)
              else
                udm_w1i = totmul / delsig / udm_bias(i)
              endif
C             ----------
              wt(ibkwt+no,ipomp+1) = wt(ibkwt+no,ipomp+1) * udm_w1i
              oldwt = wt(ibkwt+no,ipomp+1)
              wgti  = wt(ibkwt+no,ipomp+1)
              ! generate final states ------------
              udm_Kin=eein
              udm_kf_incident=ktyp
              call user_defined_interaction(21,i)
              ! ----------------------------------
              goto 1000
            endif
          enddo
          print*,"Skipped: User Defined Interaction"
          ! ------------------------------------------------------------
          endif
          ! Change the weight of the normal process
          wt(ibkwt+no,ipomp+1)=wt(ibkwt+no,ipomp+1)*udm_w2

        endif

      endif

*-----------------------------------------------------------------------
*     proton reactions
*-----------------------------------------------------------------------

      if( ityp .eq. 1 ) then

                  dmaxn = das_kmatd(kmatd(mat)+1)
                  dminn = das_kmatd(kmatd(mat)+2)

*-----------------------------------------------------------------------
*        mixed material
*-----------------------------------------------------------------------

         if( eein .gt. dminn .and. eein .le. dmaxn ) then

            do lema = 1, isigd(0)

                  r = r - siggd(lema) / delsig

               if( r .le. 0.0d0 ) then

                     icm = isigd(lema)

*-----------------------------------------------------------------------
*              nuclear data
*-----------------------------------------------------------------------

               if( icm .le. 0 ) then

                   if( iemode.eq.0 ) then ! frtati 2021/12/17

                     call cputime(9)

                     jcoll = 9 ! T.Sato 2020/09/12
                     wgti  = wt(ibkwt+no,ipomp+1)

                     ciim = aimp(ityp,iblz1,ilev1,ilat1,ii1)
                     csim = wtin(ibkwin+no,ipomp+1)
                     ffac = 0.0
                     if( ciim .gt. 0.0 ) ffac = csim / ciim


                     rh = denm(imat)
               call cp_coll(eein,wgti,ireg,imat,ffac,rh,-icm,elrt,ityp)

                  rncnt(12) = rncnt(12) + 1.0
                  call cputime(12)

                     goto 1000

                   else ! high-energy emode
                     lemc = -icm
                     rr = r + siggd(lema) / delsig
                     rr = rr - siggcn(lemc) / delsig
                     if( rr.le.0.d0 ) then
                       jcoll = 4
                       nmm = nint( dnel_das(kmat0+mat) )
                       do lem = 1, nmm
                         ztar = zz_das(kmat(mat)+lem)
                         atar = a_das(kmat(mat)+lem)
                         it_za = isigza(lemc)
                         if( it_za.eq.6000 ) it_za = 6012
                         if( 1000*ztar+atar .eq. it_za ) then
                           goto 13
                         end if
                       end do
                     else
                       call cputime(13)
                       jcoll = 3
                       iprj  = ityp
                       kprj  = ktyp
                       nmm = nint( dnel_das(kmat0+mat) )
                       do lem = 1, nmm
                         ztar = zz_das(kmat(mat)+lem)
                         atar = a_das(kmat(mat)+lem)
                         it_za = isigza(lemc)
                         if( it_za.eq.6000 ) it_za = 6012
                         if( 1000*ztar+atar .eq. it_za ) then
                           ztar = zz_das(kmat(mat)+lem)
                           atar = a_das(kmat(mat)+lem)
                           mathz = nint( ztar )
                           mathn = nint( atar - ztar )
                           call nelst(iprj,kprj,atar,ztar
     &                               ,ec(ibkec+no,ipomp+1))
                           rncnt(13) = rncnt(13) + 1.0
                           call cputime(13)
                           goto 1000
                         end if
                       end do
                     end if
                   end if

*-----------------------------------------------------------------------
*              high energy : elastic and non-elastic
*-----------------------------------------------------------------------

               else if( icm .gt. 0 ) then

                     rr = r + siggd(lema) / delsig
                     rr = rr - sigge(ksige+icm) / delsig

*-----------------------------------------------------------------------
*                 high energy elastic
*-----------------------------------------------------------------------

                  if( rr .le. 0.0d+0 ) then

                     call cputime(13)

                     jcoll = 3

                     iprj  = ityp
                     kprj  = ktyp

                     ztar = zz_das(kmat(mat)+icm)
                     atar = a_das(kmat(mat)+icm)
                     mathz = nint( ztar )
                     mathn = nint( atar - ztar )

                    call nelst(iprj,kprj,atar,ztar,ec(ibkec+no,ipomp+1))

                     rncnt(13) = rncnt(13) + 1.0
                     call cputime(13)

                     goto 1000

*-----------------------------------------------------------------------
*                 high energy non-elastic
*-----------------------------------------------------------------------

                   else

                     lem = icm
                     jcoll = 4

                     goto 13

                  end if

               end if

*-----------------------------------------------------------------------

               end if

            end do

*-----------------------------------------------------------------------
*        nuclear data for proton
*-----------------------------------------------------------------------

         else if( eein .le. dminn ) then

                if( iemode.eq.0 ) then ! frtati 2021/12/17

                  call cputime(9)

                  jcoll = 9
                  wgti  = wt(ibkwt+no,ipomp+1)

                  ciim = aimp(ityp,iblz1,ilev1,ilat1,ii1)
                  csim = wtin(ibkwin+no,ipomp+1)
                  ffac = 0.0
                  if( ciim .gt. 0.0 ) ffac = csim / ciim

                  rh = denm(imat)

                  call cp_coll(eein,wgti,ireg,imat,ffac,rh,0,elrt,ityp)

                  rncnt(12) = rncnt(12) + 1.0
                  call cputime(12)

                  goto 1000

                else ! high-energy emode
                  do lemc = 1, isigc(0)
                    r = r - siggcn(lemc) / delsig
                    if( r.le.0.d0 ) then
                      jcoll = 4
                      nmm = nint( dnel_das(kmat0+mat) )
                      do lem = 1, nmm
                        ztar = zz_das(kmat(mat)+lem)
                        atar = a_das(kmat(mat)+lem)
                        it_za = isigza(lemc)
                        if( it_za.eq.6000 ) it_za = 6012
                        if( 1000*ztar+atar .eq. it_za ) then
                          goto 13
                        end if
                      end do
                    end if
                    r = r - siggce(lemc) / delsig
                    if( r.le.0.d0 ) then
                      call cputime(13)
                      jcoll = 3
                      iprj  = ityp
                      kprj  = ktyp
                      nmm = nint( dnel_das(kmat(mat)+1) )
                      do lem = 1, nmm
                        ztar = zz_das(kmat(mat)+lem)
                        atar = a_das(kmat(mat)+lem)
                        it_za = isigza(lemc)
                        if( it_za.eq.6000 ) it_za = 6012
                        if( 1000*ztar+atar .eq. it_za ) then
                          mathz = nint( ztar )
                          mathn = nint( atar - ztar )
                          call nelst(iprj,kprj,atar,ztar
     &                               ,ec(ibkec+no,ipomp+1))
                          rncnt(13) = rncnt(13) + 1.0
                          call cputime(13)
                          goto 1000
                        end if
                      end do
                    end if
                  end do
                end if

*-----------------------------------------------------------------------
*         high energy elastic collisions only for proton
*-----------------------------------------------------------------------

         else if( eein .gt. dmaxn ) then

                  nmm = nint( dnel_das(kmat0+mat) )

            do lem = 1, nmm

                  r = r - sigge(ksige+lem) / delsig

               if( r .le. 0.0d+0 ) then

                  call cputime(13)

                  jcoll = 3

                  iprj  = ityp
                  kprj  = ktyp

                  ztar = zz_das(kmat(mat)+lem)
                  atar = a_das(kmat(mat)+lem)
                  mathz = nint( ztar )
                  mathn = nint( atar - ztar )

                  call nelst(iprj,kprj,atar,ztar,ec(ibkec+no,ipomp+1))

                  rncnt(13) = rncnt(13) + 1.0
                  call cputime(13)

                  goto 1000

               end if

            end do

*-----------------------------------------------------------------------
*           high energy proton non-elastic
*-----------------------------------------------------------------------

                  nmm = nint( dnel_das(kmat0+mat) )

            do lem = 1, nmm

                  r = r - siggn(ksign+lem) / delsig

               if( r .le. 0.0d+0 ) then

                  jcoll = 4

                  goto 13

               end if

            end do

         end if

      end if

*-----------------------------------------------------------------------
*     Charged particle by nuclear data
*-----------------------------------------------------------------------
      if(ityp .eq. 15 .or. ityp .eq. 18) then ! only deuteron and alpha
       if(dnmax(ityp) .gt. emin(ityp)) then ! dmax is specified
        if(ityp.eq.15) then
         eeinpern = eein/2.0
         dmaxn = das_kmathd(kmathd(mat)+3)
         dminn = das_kmathd(kmathd(mat)+4)
        else
         eeinpern = eein/4.0
         dmaxn = das_kmathd(kmathd(mat)+5)
         dminn = das_kmathd(kmathd(mat)+6)
        endif
        if( eeinpern.gt.dminn .and. eeinpern.le.dmaxn ) then
         do lema = 1, isigd(0)
          r = r - siggd(lema)/delsig
          if( r.le.0.d0 ) then
           icm = isigd(lema)
           if( icm.le.0 ) then ! library
            if( iemode.eq.0 ) then
             call cputime(9)
             jcoll = 9
             wgti  = wt(ibkwt+no,ipomp+1)
             ciim = aimp(ityp,iblz1,ilev1,ilat1,ii1)
             csim = wtin(ibkwin+no,ipomp+1)
             ffac = 0.0
             if( ciim .gt. 0.0 ) ffac = csim / ciim
             rh = denm(imat)
             call cp_coll(eein,wgti,ireg,imat,ffac,rh,0,elrt,ityp)
             rncnt(12) = rncnt(12) + 1.0
             call cputime(12)
             goto 1000
            else
             lemc = -icm
             jcoll = 4
             nmm = nint( dnel_das(kmat0+mat) )
             do lem = 1, nmm
              ztar = zz_das(kmat(mat)+lem)
              atar = a_das(kmat(mat)+lem)
              it_za = isigza(lemc)
              if( it_za.eq.6000 ) it_za = 6012
              if( 1000*ztar+atar .eq. it_za ) then
               goto 13
              end if
             end do
            end if
           else if( icm.gt.0 ) then ! model non-elastic
            rr = r + siggd(lema)/delsig
            rr = rr - siggn(ksign+icm)/delsig
            if( rr .le. 0.0d+0 ) then
             lem = icm
             jcoll = 4
             goto 13
            end if
           end if
          end if
         end do
        else if( eeinpern.le.dminn ) then
         if( iemode.eq.0 ) then
          call cputime(9)
          jcoll = 9
          wgti  = wt(ibkwt+no,ipomp+1)
          ciim = aimp(ityp,iblz1,ilev1,ilat1,ii1)
          csim = wtin(ibkwin+no,ipomp+1)
          ffac = 0.0
          if( ciim .gt. 0.0 ) ffac = csim / ciim
          rh = denm(imat)
          call cp_coll(eein,wgti,ireg,imat,ffac,rh,0,elrt,ityp)
          rncnt(12) = rncnt(12) + 1.0
          call cputime(12)
          goto 1000
         else
          do lemc = 1, isigc(0)
           r = r - siggcn(lemc) / delsig
           if( r .le. 0.0d+0 ) then
            jcoll = 4
            nmm = nint( dnel_das(kmat0+mat) )
            do lem = 1, nmm
             ztar = zz_das(kmat(mat)+lem)
             atar = a_das(kmat(mat)+lem)
             it_za = isigza(lemc)
             if( it_za.eq.6000 ) it_za = 6012
             if( 1000*ztar+atar .eq. it_za ) then
              goto 13
             end if
            end do
           end if
          end do
         end if
        else if( eeinpern.gt.dmaxn ) then
         nmm = nint( dnel_das(kmat0+mat) )
         do lem = 1, nmm
          r = r - siggn(ksign+lem) / delsig
          if( r .le. 0.0d+0 ) then
           jcoll = 4
           goto 13
          end if
         end do
        end if
       end if
      end if
*-----------------------------------------------------------------------
*     neutron reactions
*-----------------------------------------------------------------------

      if( ityp .eq. 2 ) then

                  dmaxn = das_kmatd(kmatd(mat)+3)
                  dminn = das_kmatd(kmatd(mat)+4)

*-----------------------------------------------------------------------
*        mixed material
*-----------------------------------------------------------------------

         if( eein .gt. dminn .and. eein .le. dmaxn ) then

            do lema = 1, isigd(0)

                  r = r - siggd(lema) / delsig

               if( r .le. 0.0d0 ) then

                     icm = isigd(lema)

*-----------------------------------------------------------------------
*              nuclear data
*-----------------------------------------------------------------------

               if( icm .le. 0 ) then

                   if( iemode.eq.0 .or. eein.le.emodem ) then ! frtati 2021/12/17

                     call cputime(9)

                     wgti  = wt(ibkwt+no,ipomp+1)

                     ciim = aimp(ityp,iblz1,ilev1,ilat1,ii1)
                     csim = wtin(ibkwin+no,ipomp+1)
                     ffac = 0.0
                     if( ciim .gt. 0.0 ) ffac = csim / ciim

                     if( iemode .eq. 0 .or. eein .gt. emodem ) then

                        jcoll = 6

                        call sctneu(eein,wgti,ireg,imat,ffac,-icm,elrt)

                     else

                        jcoll = 10

                        call sctneut(eein,wgti,ireg,imat,ffac1,
     &                                              ffac2,-icm,elrt)

                     end if

                     rncnt(9) = rncnt(9) + 1.0
                     call cputime(9)

                     goto 1000

                   else ! high-energy emode
                     lemc = -icm
                     rr = r + siggd(lema) / delsig
                     rr = rr - siggcn(lemc) / delsig
                     if( rr.le.0.d0 ) then
                       jcoll = 4
                       nmm = nint( dnel_das(kmat0+mat) )
                       do lem = 1, nmm
                         ztar = zz_das(kmat(mat)+lem)
                         atar = a_das(kmat(mat)+lem)
                         it_za = isigza(lemc)
                         if( it_za.eq.6000 ) it_za = 6012
                         if( 1e3*ztar+atar .eq. it_za ) then
                           goto 13
                         end if
                       end do
                     else
                       call cputime(13)
                       jcoll = 3
                       iprj  = ityp
                       kprj  = ktyp
                       nmm = nint( dnel_das(kmat0+mat) )
                       do lem = 1, nmm
                         ztar = zz_das(kmat(mat)+lem)
                         atar = a_das(kmat(mat)+lem)
                         it_za = isigza(lemc)
                         if( it_za.eq.6000 ) it_za = 6012
                         if( 1000*ztar+atar .eq. it_za ) then
                           ztar = zz_das(kmat(mat)+lem)
                           atar = a_das(kmat(mat)+lem)
                           mathz = nint( ztar )
                           mathn = nint( atar - ztar )
                           call nelst(iprj,kprj,atar,ztar
     &                               ,ec(ibkec+no,ipomp+1))
                           rncnt(13) = rncnt(13) + 1.0
                           call cputime(13)
                           goto 1000
                         end if
                       end do
                     end if
                   end if

*-----------------------------------------------------------------------
*              high energy : elastic and non-elastic
*-----------------------------------------------------------------------

               else if( icm .gt. 0 ) then

                     rr = r + siggd(lema) / delsig
                     rr = rr - sigge(ksige+icm) / delsig

*-----------------------------------------------------------------------
*                 high energy elastic
*-----------------------------------------------------------------------

                  if( rr .le. 0.0d+0 ) then

                     call cputime(13)

                     jcoll = 3

                     iprj  = ityp
                     kprj  = ktyp

                     ztar = zz_das(kmat(mat)+icm)
                     atar = a_das(kmat(mat)+icm)
                     mathz = nint( ztar )
                     mathn = nint( atar - ztar )

                    call nelst(iprj,kprj,atar,ztar,ec(ibkec+no,ipomp+1))

                     rncnt(13) = rncnt(13) + 1.0
                     call cputime(13)

                     goto 1000

*-----------------------------------------------------------------------
*                 high energy non-elastic
*-----------------------------------------------------------------------

                   else

                     lem = icm
                     jcoll = 4

                     goto 13

                  end if

               end if

*-----------------------------------------------------------------------

               end if

            end do

*-----------------------------------------------------------------------
*        nuclear data for neutron
*-----------------------------------------------------------------------

         else if( eein .le. dminn ) then

                if( iemode.eq.0 .or. eein.le.emodem ) then ! frtati 2021/12/17

                  call cputime(9)

                  wgti  = wt(ibkwt+no,ipomp+1)

                  ciim = aimp(ityp,iblz1,ilev1,ilat1,ii1)
                  csim = wtin(ibkwin+no,ipomp+1)
                  ffac = 0.0
                  if( ciim .gt. 0.0 ) ffac = csim / ciim

               if( iemode .eq. 0 .or. eein .gt. emodem ) then

                  jcoll = 6

                  call sctneu(eein,wgti,ireg,imat,ffac,0,elrt)

               else

                  jcoll = 10

                  call sctneut(eein,wgti,ireg,imat,ffac1,ffac2,0,elrt)

               end if

                  rncnt(9) = rncnt(9) + 1.0
                  call cputime(9)

                  goto 1000

                else ! high-energy emode
                  do lemc = 1, isigc(0)
                    r = r - siggcn(lemc) / delsig
                    if( r.le.0.d0 ) then
                      jcoll = 4
                      nmm = nint( dnel_das(kmat0+mat) )
                      do lem = 1, nmm
                        ztar = zz_das(kmat(mat)+lem)
                        atar = a_das(kmat(mat)+lem)
                        it_za = isigza(lemc)
                        if( it_za.eq.6000 ) it_za = 6012
                        if( 1000*ztar+atar .eq. it_za ) then
                          goto 13
                        end if
                      end do
                    end if
                    r = r - siggce(lemc) / delsig
                    if( r.le.0.d0 ) then
                      call cputime(13)
                      jcoll = 3
                      iprj  = ityp
                      kprj  = ktyp
                      nmm = nint( dnel_das(kmat0+mat) )
                      do lem = 1, nmm
                        ztar = zz_das(kmat(mat)+lem)
                        atar = a_das(kmat(mat)+lem)
                        it_za = isigza(lemc)
                        if( it_za.eq.6000 ) it_za = 6012
                        if( 1000*ztar+atar .eq. it_za ) then
                          mathz = nint( ztar )
                          mathn = nint( atar - ztar )
                          call nelst(iprj,kprj,atar,ztar
     &                               ,ec(ibkec+no,ipomp+1))
                          rncnt(13) = rncnt(13) + 1.0
                          call cputime(13)
                          goto 1000
                        end if
                      end do
                    end if
                  end do
                end if

*-----------------------------------------------------------------------
*         high energy elastic collisions only for nuetron
*-----------------------------------------------------------------------

         else if( eein .gt. dmaxn ) then

                  nmm = nint( dnel_das(kmat0+mat) )

            do lem = 1, nmm

                  r = r - sigge(ksige+lem) / delsig

               if( r .le. 0.0d+0 ) then

                  call cputime(13)

                  jcoll = 3

                  iprj  = ityp
                  kprj  = ktyp

                  ztar = zz_das(kmat(mat)+lem)
                  atar = a_das(kmat(mat)+lem)
                  mathz = nint( ztar )
                  mathn = nint( atar - ztar )

                  call nelst(iprj,kprj,atar,ztar,ec(ibkec+no,ipomp+1))

                  rncnt(13) = rncnt(13) + 1.0
                  call cputime(13)

                  goto 1000

               end if

            end do

*-----------------------------------------------------------------------
*           high energy neutron non-elastic
*-----------------------------------------------------------------------

                  nmm = nint( dnel_das(kmat0+mat) )

            do lem = 1, nmm

                  r = r - siggn(ksign+lem) / delsig

               if( r .le. 0.0d+0 ) then

                  jcoll = 4

                  goto 13

               end if

            end do

         end if

      end if

*-----------------------------------------------------------------------
*     photon
*-----------------------------------------------------------------------

      if( ityp .eq. 14 .and. eein .le. dnmax(14) ) then

*-----------------------------------------------------------------------
*     manage reaction channel of photon interaction
*        pnimul: multiplying factor for photonuclear interaction CS
*        gmumul: multiplying factor for gamma-induced muon pair production CS
*-----------------------------------------------------------------------

         if( ( ipnint  .ge. 1 .and. pnimul .gt. 1.d0 ) .or.
     &       ( igmuppd .ge. 1 .and. gmumul .gt. 1.0d0 ) .or.
     &         udm_int_num .gt. 0  ) then

            totmul = totpni * pnimul + totgmm * gmumul + totgam + totegs
     &             + totudm_mul

            wtmg1 = totmul / delsig / pnimul
            wtmg2 = totmul / delsig
            wtmg3 = totmul / delsig / gmumul

            prbg1 = totpni * pnimul / totmul
            prbg2 = ( totgam + totegs ) / totmul
            prbg3 = totgmm * gmumul / totmul
            prbg_udp = totudm_mul / totmul


*-----------------------------------------------------------------------
*           1 = prbg1 * wtmg1 + prbg2 * wtmg2 + prbg3 * wtmg3
*-----------------------------------------------------------------------

            if( r. lt. prbg1 ) then


cfrtati 2021/12/17 added photo-nuclear library part (cg_coll needed)
              wt(ibkwt+no,ipomp+1) = wt(ibkwt+no,ipomp+1) * wtmg1
              oldwt = wt(ibkwt+no,ipomp+1)
              wgti  = wt(ibkwt+no,ipomp+1)

              dmaxn = das_kmathd(kmathd(mat)+1)
              dminn = das_kmathd(kmathd(mat)+2)
              if( eein.gt.dminn .and. eein.le.dmaxn ) then
                rr = r
                do lema = 1, isigd(0)
                  rr = rr - siggd(lema)*pnimul/totmul
                  if( rr.le.0.d0 ) then
                    icm = isigd(lema)
                    if( icm.le.0 ) then
                      if( iemode.eq.0 ) then
                        jcoll = 9
                        call cg_coll(eein,wgti,ireg,imat,-icm)
                        goto 1000
                      else
                        jcoll = 15
                        call sctpni(eein,wgti,ireg,imat,-icm)
                        goto 1000
                      end if
                    else if( icm.gt.0 ) then
                      jcoll = 15
                      call sctpni(eein,wgti,ireg,imat,icm)
                      goto 1000
                    end if
                  endif
                end do
              else if( eein.le.dminn ) then
                if( iemode.eq.0 ) then
                  jcoll = 9
                  rr = r
                  do lem = 1, isigc(0)
                    rr = rr - siggc(lem)*pnimul/totmul
                    if( rr.le.0.d0 ) then
                      call cg_coll(eein,wgti,ireg,imat,lem)
                      goto 1000
                    end if
                  end do
                else
                  jcoll = 15
                  call sctpni(eein,wgti,ireg,imat,0)
                  goto 1000
                end if
              else if( eein.gt.dmaxn ) then

                  jcoll = 15   ! S.Abe 2017/01/17
                  call sctpni(eein,wgti,ireg,imat,0) ! frtati 2021/12/17 added icm
                  goto 1000

              end if

*-----------------------------------------------------------------------

            elseif( r. lt. prbg3 + prbg1 ) then      ! y.sakaki 2021/05/06

! y.sakaki 2023/07/19 -->
                  jcoll = 19
                  wt(ibkwt+no,ipomp+1) = wt(ibkwt+no,ipomp+1) * wtmg3
                  oldwt = wt(ibkwt+no,ipomp+1)                 ! cKN 2019/05/27
                  wgti  = wt(ibkwt+no,ipomp+1)
                  call sctgmm(eein,wgti,ireg,imat)
                  goto 1000

*-----------------------------------------------------------------------
            elseif( r .lt. prbg3 + prbg1 + prbg_udp ) then
            ! Choose a process from user-defined processes ( totudm_mul_(i) )
            ! Calculate weight for the choosed process
              tmp=0d0
              r_udp=rn(0)
              call fill_mat_info(imat)
              do i = 1, udm_int_num
                tmp=tmp + (totudm_(i) * udm_bias(i)) / totudm_mul
C       ! If the charged particle is an incident particle, the energy may have
C       ! changed compared to that used to calculate totudm_(i) in getflt.f.
C       ! Processes with a cross section of 0 for the changed energy are eliminated.
C       ! For example, this can be a problem for production processes in narrow
C       ! resonance processes. This could be addressed by
C       ! (1) Shortening the following step lengths
C       !     chard (Default value = 0.1. For electrons and positrons)
C       !     deltm (Default value = 20.12345. For others)
C       ! (2) Comment out 'if(udm_sigt==0d0) cycle' and address the problem in the udm_int* file.
                if(r_udp .lt. tmp) then ! index=i is accepted.
                  jcoll = 20
                  wtmg_udp = totmul / delsig / udm_bias(i)
                  wt(ibkwt+no,ipomp+1) = wt(ibkwt+no,ipomp+1) * wtmg_udp
                  oldwt = wt(ibkwt+no,ipomp+1)
                  wgti  = wt(ibkwt+no,ipomp+1)
                  ! generate final states ------------
                  udm_Kin=eein
                  udm_kf_incident=ktyp
                  call user_defined_interaction(21,i)
                  ! ----------------------------------
                  goto 1000
                endif
              enddo
              print*,"Skipped: User Defined Interaction"

*-----------------------------------------------------------------------

            else

               if( iegsemi .eq. 0 ) then

                  jcoll = 7
                  wt(ibkwt+no,ipomp+1) = wt(ibkwt+no,ipomp+1) * wtmg2
                  oldwt = wt(ibkwt+no,ipomp+1)                 ! cKN 2019/05/27
                  wgti  = wt(ibkwt+no,ipomp+1)

                  ielpos = 0                         ! S.Abe 2017/02/08

                  call sctgam(eein,wgti,ireg,imat)

                  goto 1000

               else

                  call cputime(24)
                  jcoll = 13 ! HI defined
                   u0 =  u(ibku + no,ipomp+1)
                   v0 =  v(ibkv + no,ipomp+1)
                   w0 =  w(ibkw + no,ipomp+1)
                  wt(ibkwt+no,ipomp+1) = wt(ibkwt+no,ipomp+1) * wtmg2
                  oldwt = wt(ibkwt+no,ipomp+1)                 ! cKN 2019/05/27
                  wt0 = wt(ibkwt+no,ipomp+1)
                  call egs5pcoll(eein,u0,v0,w0,wt0,ireg,imat,nbeta)
                  rncnt(24) = rncnt(24) + 1.0
                  call cputime(24)
                  goto 1000

               endif

            endif

*-----------------------------------------------------------------------
*     default process
*-----------------------------------------------------------------------

         else

*-----------------------------------------------------------------------
*     photon for nuclear data
*-----------------------------------------------------------------------

            if ( iegsemi .eq. 0 ) then

               r = r - totgam / delsig

               if( r .le. 0.0d0 ) then

                  jcoll = 7
                  wgti  = wt(ibkwt+no,ipomp+1)

                  ielpos = 0  ! S.Abe 2017/02/08

                  call sctgam(eein,wgti,ireg,imat)

                  goto 1000

               endif

            end if

*---------------------------------------------------------------------
*     photon for EGS5
*---------------------------------------------------------------------

            if( iegsemi .ne. 0 ) then

               r = r - totegs / delsig

               if( r .le. 0.0d0 ) then

                  call cputime(24)

                  jcoll = 13 ! HI defined

                   u0 =  u(ibku + no,ipomp+1)
                   v0 =  v(ibkv + no,ipomp+1)
                   w0 =  w(ibkw + no,ipomp+1)
                  wt0 = wt(ibkwt+ no,ipomp+1)

                  call egs5pcoll(eein,u0,v0,w0,wt0,ireg,imat,nbeta)


                  rncnt(24) = rncnt(24) + 1.0
                  call cputime(24)

                  goto 1000

               endif

            endif

*-----------------------------------------------------------------------
*        photon for photonuclear reaction
*-----------------------------------------------------------------------

            if( ipnint .ge. 1 ) then

               wgti  = wt(ibkwt+no,ipomp+1)

               r = r - totpni / delsig

               if( r .le. 0.0d0 ) then

cfrtati 2021/12/17 added photo-nuclear library part
                dmaxn = das_kmathd(kmathd(mat)+1)
                dminn = das_kmathd(kmathd(mat)+2)
                if( eein.gt.dminn .and. eein.le.dmaxn ) then
                  rr = r + totpni / delsig
                  do lema = 1, isigd(0)
                    rr = rr - siggd(lema)/delsig
                    if( rr.le.0.d0 ) then
                      icm = isigd(lema)
                      if( icm.le.0 ) then ! library
                        if( iemode.eq.0 ) then
                          jcoll = 9
                          call cg_coll(eein,wgti,ireg,imat,-icm)
                          goto 1000
                        else
                          jcoll = 15
                          call sctpni(eein,wgti,ireg,imat,-icm)
                          goto 1000
                        end if
                      else if( icm.gt.0 ) then ! model
                        jcoll = 15
                        call sctpni(eein,wgti,ireg,imat,icm)
                        goto 1000
                      end if
                    end if
                  end do

                else if( eein.le.dminn ) then
                  if( iemode.eq.0 ) then
                    jcoll = 9
                    rr = r + totpni / delsig
                    do lem = 1, isigc(0)
                      rr = rr - siggc(lem)/delsig
                      if( rr.le.0.d0 ) then
                        call cg_coll(eein,wgti,ireg,imat,lem)
                        goto 1000
                      end if
                    end do
                  else
                    jcoll = 15
                    call sctpni(eein,wgti,ireg,imat,0)
                    goto 1000
                  end if
                else if( eein.gt.dmaxn ) then

                  jcoll = 15   ! S.Abe 2017/01/17
                  wgti  = wt(ibkwt+no,ipomp+1)

                  call sctpni(eein,wgti,ireg,imat,0) ! frtati 2021/12/17 added icm

                  goto 1000

                end if

               endif

            endif

*-----------------------------------------------------------------------
*        photon for muon pair production
*-----------------------------------------------------------------------

            if( igmuppd .ge. 1 ) then

               r = r - totgmm / delsig

               if( r .le. 0.0d0 ) then

                  jcoll = 19   ! y.sakaki 2019/09/26
                  wgti  = wt(ibkwt+no,ipomp+1) ! basically 1. not understood yet.

                  call sctgmm(eein,wgti,ireg,imat)

                  goto 1000

               endif

            endif

*-----------------------------------------------------------------------

         endif

         goto 2000   ! S.Abe 2021/06/03, avoid to go to ncasc

      end if

*-----------------------------------------------------------------------
*     electron for track structure
*-----------------------------------------------------------------------

         if(  ( ityp .eq. 12 .or. ityp .eq. 13 ).and.
     &       mat1 .ne. 0                        .and.
     &       e(ibke+no,ipomp+1) .le. etsmax  ) then  ! use e instead of ec because cross section is not ready
           jcoll = 18 ; kcoll = 0 ; itcoll = 1

chirata etsmode02 20210915
          if(mat1.eq.-1 .and. 
     &       ichem(mat,1) .eq. 22014 .and.ichem(mat,2) .eq. 0) then 
           call etsreac02    ! ETS for Si
          else if(mat1.eq.-1 .and. 
     &       ichem(mat,1) .eq. 1 .and.ichem(mat,2) .eq. 0) then 
           call etsreac    ! water-ETS with mID=-1, chem=water
          else if(mat1.eq.-1) then ! hirata 20220829 arbitary ETS mode
           call etsreac_art
          else
           call etsreac
          endif

           rncnt(29) = rncnt(29) + 1.0
           call cputime(29)
           goto 1000
         end if

*-----------------------------------------------------------------------
*     electron for nuclear data
*-----------------------------------------------------------------------

         if( ( ityp .eq. 12 .or. ityp .eq. 13 ) .and.
     &         eein .le. dnmax(12) .and. iegsemi .eq. 0 ) then

            if( eein .le. dnmax(12) ) then

                  call cputime(11)

                  jcoll = 8
                  wgti  = wt(ibkwt+no,ipomp+1)

                     xv(1) = x(ibkx+no,ipomp+1) - xc(ibkxc+no,ipomp+1)
                     xv(2) = y(ibky+no,ipomp+1) - yc(ibkyc+no,ipomp+1)
                     xv(3) = z(ibkz+no,ipomp+1) - zc(ibkzc+no,ipomp+1)

                     xvs = sqrt( xv(1)**2 + xv(2)**2 + xv(3)**2 )

                  if( xvs .gt. 0.0d0 ) then

                     xv(1) = xv(1) / xvs
                     xv(2) = xv(2) / xvs
                     xv(3) = xv(3) / xvs

                  end if

                  if( ityp .eq. 12 ) then

                     ielpos = -1

                  else if( ityp .eq. 13 ) then

                     ielpos = 1

                  end if

                     rh = denm(imat)

                  eezo  = e(ibke+no,ipomp+1)

                  call sctelc(eein,wgti,ireg,imat,xv,rh,eezo)

                     rncnt(11) = rncnt(11) + 1.0
                     call cputime(11)

            end if

                  goto 1000

         end if


*---------------------------------------------------------------------
*     electron for EGS5
*---------------------------------------------------------------------

         if( ( ityp .eq. 12 .or. ityp .eq. 13 ) .and.
     &         eein .le. dnmax(12) .and. iegsemi .ne. 0 ) then


            if( eein .le. dnmax(14) ) then

                  call cputime(38)

                  if(ityp.eq.12)lelec = -1
                  if(ityp.eq.13)lelec =  1

                   e0 = eein
                   u0 =  u(ibku + no,ipomp+1)
                   v0 =  v(ibkv + no,ipomp+1)
                   w0 =  w(ibkw + no,ipomp+1)
                  wt0 = wt(ibkwt+ no,ipomp+1)

                  call egs5ecoll(e0,u0,v0,w0,wt0,lelec,ireg,mark,nbeta)
                  jcoll = 14 ! iwase defined

                     rncnt(25) = rncnt(25) + 1.0
                     call cputime(25)

            end if

                  goto 1000

         end if


*-----------------------------------------------------------------------
*     proton for nuclear data before introducing [data max]. Used only for proton library with EG mode
*-----------------------------------------------------------------------
         if( ityp .eq. 1 .and. eein .lt. dnmax(1) ) then
                  call cputime(12)

                  jcoll = 9
                  wgti  = wt(ibkwt+no,ipomp+1)
                  ciim = aimp(ityp,iblz1,ilev1,ilat1,ii1)
                  csim = wtin(ibkwin+no,ipomp+1)
                  ffac = 0.0
                  if( ciim .gt. 0.0 ) ffac = csim / ciim
                  rh = denm(imat)
                  call cp_coll(eein,wgti,ireg,imat,ffac,rh,0,elrt,ityp) ! T.Sato 2021/08/14
                  rncnt(12) = rncnt(12) + 1.0
                  call cputime(12)

                  goto 1000

         end if


*-----------------------------------------------------------------------
*     muon interaction
*-----------------------------------------------------------------------

         if( (ityp .eq. 6 .or. ityp .eq. 7) .and.
     &       totttl .gt. 0.d0 ) then

          do i = 1, nmuh

           r = r - sigmuv(i) / delsig
           if( r .le. 0.d0 ) then
            call cputime(23)
            jcoll = 17   ! S.Abe 2017/01/17
            eein = ec(ibkec+no,ipomp+1)
            wgti = wt(ibkwt+no,ipomp+1)
            ztar = itzmuh(i)
            atar = itamuh(i)
            call vpimu(eein,wgti,ztar,atar,ireg,imat)
                  rncnt(23) = rncnt(23) + 1.0
                  call cputime(23)
            goto 1000
           endif

           r = r - sigmub(i) / delsig
           if( r .le. 0.d0 ) then
            call cputime(27)
            jcoll = 12
            eein = ec(ibkec+no,ipomp+1)
            wgti = wt(ibkwt+no,ipomp+1)
            ztar = itzmuh(i)
            atar = itamuh(i)
            call brmmu(eein,wgti,ztar,atar)
                  rncnt(27) = rncnt(27) + 1.0
                  call cputime(27)
            goto 1000
           endif

           r = r - sigmup(i) / delsig
           if( r .le. 0.d0 ) then
            call cputime(28)
            jcoll = 12
            eein = ec(ibkec+no,ipomp+1)
            wgti = wt(ibkwt+no,ipomp+1)
            ztar = itzmuh(i)
            atar = itamuh(i)
            call ppdmu(eein,wgti,ztar,atar,ireg,imat)
                  rncnt(28) = rncnt(28) + 1.0
                  call cputime(28)
            goto 1000
           endif

          enddo

          if( jtyp .ne. 0 .and. rtyp .gt. 1.d0 .and.
     &        mndel .gt. 0 .and. delm(icl) .lt. 1.d+10 ) then
          else
           goto 2000
          endif

         endif

*-----------------------------------------------------------------------
*     Reaction of Neutrino + Deuteron
*-----------------------------------------------------------------------

         if( abs(ktyp) .eq. 12 .or. abs(ktyp) .eq. 14 .or.
     &    abs(ktyp) .eq. 16 ) then ! 2017/11/20 Ogawa Neutrino reaction
                  nmm = nint( dnel_das(kmat0+mat) )

            do lem = 1, nmm

                  r = r - siggn(ksige+lem) / delsig

               if( r .le. 0.0d+0 ) then

                  call cputime(14)

                  jcoll = 4

                  ztar = zz_das(kmat(mat)+lem)
                  atar = a_das(kmat(mat)+lem)
                  mathz = nint( ztar )
                  mathn = nint( atar - ztar )

                   call neutrino_kinem2(ktyp,int(atar),int(ztar),
     &              ec(ibkec+no,ipomp+1))

                  rncnt(14) = rncnt(14) + 1.0  ! define new category sometime
                  call cputime(14)  ! define new category sometime

                  goto 1000

               end if

            end do

            goto 2000   ! S.Abe 2021/06/03, avoid to go to ncasc

         endif

*-----------------------------------------------------------------------
*     proton and ion for track structure
*-----------------------------------------------------------------------

         if(ibryf(ityp,ktyp) .ge. 1) then
          r = r - tottsm / delsig
          if(r .lt. 0.d0) then

           if    ( ityp .eq. 1        .and.
     &       mat1 .eq. 1              .and.
     &       e(ibke+no,ipomp+1) .le. ptsmax  )then
             jcoll = 18 ; kcoll = 0 ; itcoll = 2
             call pts_reac

           elseif( ktyp .eq. 6000012  .and.
     &       mat1 .eq. 1              .and.
     &       e(ibke+no,ipomp+1) .le. ctsmax  )then
             jcoll = 18 ; kcoll = 0; itcoll = 3
             call cts_reac

           elseif( ibryf(ityp,ktyp) .ge. 1 .and. ! Rudd based track-structure mode
     &       mat1 .ne. 0                   .and.
     &       e(ibke+no,ipomp+1) .le. tsmax*dble(ibryf(ityp,ktyp)))then
             jcoll = 18 ; kcoll = 0 ; itcoll = 4 ; ierr = 0
             call itsreac(mat1,ierr)
             if(ierr .eq. 1) then
                 jcoll = 0 ; kcoll = 0 ; itcoll = 0
             endif
           endif

           if( jcoll .eq. 18) then
             rncnt(30) = rncnt(30) + 1.0
             call cputime(30)
             goto 1000
           endif

          endif

         endif
*-----------------------------------------------------------------------
*     Delta electron from charged particles
*-----------------------------------------------------------------------

         if( jtyp .ne. 0 .and. rtyp .gt. 1.d0 .and.
     &       mndel .gt. 0 .and. delm(icl) .lt. 1.d+10 ) then

                  r = r - totdel / delsig

               if( r .le. 0.0d+0 ) then

                  jcoll = 11

                  iprj  = ityp
                  kprj  = ktyp
                  rmass = rtyp
                  chag  = dabs(dble( jtyp ))
                  eezo  = e(ibke+no,ipomp+1)
                  eein  = ec(ibkec+no,ipomp+1)
                  emind = delm(icl)

                  call delprd(iprj,kprj,rmass,chag,eein,emind,eezo)


                  goto 1000

               end if

         end if

*-----------------------------------------------------------------------
*     Non-Nucleons, Non-Elastic Collisions
*-----------------------------------------------------------------------

         if( ityp .gt. 2 ) then

                  nmm = nint( dnel_das(kmat0+mat) )

            do lem = 1, nmm

                  r = r - siggn(ksign+lem) / delsig

               if( r .le. 0.0d+0 ) then

                  jcoll = 4

                  goto 13

               end if

            end do

         end if

*-----------------------------------------------------------------------

               goto 1000

*-----------------------------------------------------------------------
*     Nuclear Non-Elastic
*     =====  ncasc  =====
*-----------------------------------------------------------------------

   13 continue

                  ipim = 0
                  idwbahit = 1

   22          continue

                  call cputime(14)

*-----------------------------------------------------------------------

               if( invp .eq. 0 ) then

                     ztar = zz_das(kmat(mat)+lem)
                     atar = a_das(kmat(mat)+lem)

                     iprj  = ityp
                     kprj  = ktyp
                     mmas  = nint( atar )
                     mchg  = nint( ztar )

                     eein  = ec(ibkec+no,ipomp+1)
                     if( iabsms .eq. 1 .and. eein .lt. emin(ityp) )
     &               eein  = emin(ityp)

                  if( ityp .lt. 15 ) then

                     bmax = sqrt( geosig(mmas) * 100.0 / 3.1415926 )

                  else if( ityp .ge. 15 ) then

                     ihvi = ityp - 10

                     ap = dble( ibryf(ityp,ktyp) )
                     zp = dble( ichgf(ityp,ktyp) )

                     call sighi(ap,zp,eein,atar,ztar,signe,sigel,bmax)

*                  inverse kinematics for light ion targets
                   if( (( atar.eq.2 .and. ztar.eq.1 )
     &                    .or. ( atar.eq.3 .and. ztar.eq.1 )
     &                    .or. ( atar.eq.3 .and. ztar.eq.2 )
     &                    .or. ( atar.eq.4 .and. ztar.eq.2 ))
     &                    .and. atar.lt.ap ) then ! T.Sato 2016/03/14, not inverse for lighter target

                      invpl  = 1

                    if ( atar.eq.2 .and. ztar.eq.1 ) then

                       iprj = 15

                    else if ( atar.eq.3 .and. ztar.eq.1 ) then

                       iprj = 16

                    else if ( atar.eq.3 .and. ztar.eq.2 ) then

                       iprj = 17

                    else if ( atar.eq.4 .and. ztar.eq.2 ) then

                       iprj = 18

                    end if

                      kprj = idnint(ztar)*1000000+idnint(atar)
                      mmas = ibryf(ityp,ktyp)
                      mchg = ichgf(ityp,ktyp)

                      rlmass = rstms(iprj)*1000d0
                      eein = rlmass / rtyp * ec(ibkec+no,ipomp+1)
                      e2p  = eein + rlmass
                      p2p  = dsqrt( eein**2 + 2.0 * rlmass * eein )
                      beta = p2p / e2p
                      gamm = e2p / rlmass

                   end if

                  end if

               else if( invp .ne. 0 ) then

                     atar = 1.0
                     ztar = 1.0

                     iprj = 1
                     kprj = 2212
                     mmas = ibryf(ityp,ktyp)
                     mchg = ichgf(ityp,ktyp)

                     bmax = sqrt( geosig(mmas) * 100.0 / 3.1415926 )

                     eein = rpmass / rtyp * ec(ibkec+no,ipomp+1)
                     e2p  = eein + rpmass
                     p2p  = sqrt( eein**2 + 2.0 * rpmass * eein )
                     beta = p2p / e2p
                     gamm = e2p / rpmass

               end if

                     mathz = nint( ztar )
                     mathn = nint( atar - ztar )

                     if( e(ibke+no,ipomp+1) .lt. 0.0 ) bmax = 0.0d0

*-----------------------------------------------------------------------
*           Reactions by Fragment Data : opt = 1
*-----------------------------------------------------------------------

            if( ifrgd .gt. 0 ) then

             do i = 1, ifrgd

              if( ifgdf(i,5) .ne. 1 ) then

CS.Hashimoto revised for extension of kind of particles. (2014.9.24)
              if( ifgdf(i,1) .eq. 1 .or. ifgdf(i,1) .eq. 3
     &                .or. ifgdf(i,1) .eq. 4
     &                .or. ifgdf(i,1) .eq. 5 ) then

               if( ityp .eq. kftp(ifgdf(i,2)) ) then

                if ( ifgdf(i,3) .eq. mathz*1000000+mathn+mathz .or.
     &               ifgdf(i,3) .eq. 2212 ) then

                   ifg = i
                   efg = ec(ibkec+no,ipomp+1) / mtyp

                   call nfrgdat(ifg,efg)
                  rncnt(31) = rncnt(31) + 1.0

                   goto 500

                end if

               end if

              end if

              end if

             end do

            end if

*-----------------------------------------------------------------------
! SCINFUL mode, T.Sato 2020/02/16
! modified by d.satoh (2021.08.26)
*-----------------------------------------------------------------------
         if(iswitch.ne.0) then
           if(iprj.eq.2.and.mchg.eq.6.and.mmas.eq.12) then ! neutron on
             call scinmode(eein) ! SCINFUL mode
             rncnt(32) = rncnt(32) + 1.0
             goto 500
           endif
         endif
*-----------------------------------------------------------------------
*           Reactions by Cascade or QMD
*-----------------------------------------------------------------------


               call ncasc(0,iprj,kprj,eein,mmas,mchg,bmax)

*-----------------------------------------------------------------------
*              final transform for inverse proton-nucleus reaction
*-----------------------------------------------------------------------

               if( ( invp .ne. 0 .and. nclst .gt. 0 )
     &              .or. ( invpl .ne. 0 .and. nclst .gt. 0 ) ) then

                  do i = 1, nclst

                     px = qclust(1,i)
                     py = qclust(2,i)
                     pz = qclust(3,i)
                     et = qclust(4,i)
                     rm = qclust(5,i)

                     pz = - ( gamm * pz - beta * gamm * et )
                     et = sqrt( px**2 + py**2 + pz**2 + rm**2 )

                     qclust(3,i) = pz
                     qclust(4,i) = et
                     qclust(7,i) = ( et - rm ) * 1000.

                  end do

               end if

*-----------------------------------------------------------------------
*           Add Fragment Data : opt = 2
*-----------------------------------------------------------------------

            if( ifrgd .gt. 0 ) then

               do i = 1, ifrgd

              if( ifgdf(i,5) .ne. 1 ) then

CS.Hashimoto revised for extension of kind of particles. (2014.9.24)
              if( ifgdf(i,1) .eq. 2 ) then

               if( ityp .eq. kftp(ifgdf(i,2)) ) then

                if ( ifgdf(i,3) .eq. mathz*1000000+mathn+mathz .or.
     &               ifgdf(i,3) .eq. 2212 ) then

                   ifg = i
                   efg = ec(ibkec+no,ipomp+1) / mtyp

                   call nfrgdat(ifg,efg)
                  rncnt(31) = rncnt(31) + 1.0

                   goto 500

                end if

               end if

              end if

              end if

             end do

            end if

  500       continue

*-----------------------------------------------------------------------

                  rncnt(14) = rncnt(14) + 1.0
                  call cputime(14)

*-----------------------------------------------------------------------

               if( nclst .lt. 0 ) then

                  if( iabsms .ne. 1 ) then ! for all particles (S.H.20140814)

                     ipim = ipim + 1

                     if( ipim .le. 100 .and. nclst .gt. -2 ) goto 22 ! increase from 20 to 100 T.Sato 2022/2/20

                  else if( iabsms .eq. 1 ) then

                     ipim = ipim + 1

                   if ( ityp .eq. 5 ) then ! for low-energy absorption of pion, S.H. 2022.1.14
                     if( ipim .le. 1000 ) goto 22
                   else
                     if( ipim .le. 10 ) goto 22
                   end if

                  end if

               end if

*-----------------------------------------------------------------------
*              absorption of negative charge mesons
*-----------------------------------------------------------------------

               if( nclst .gt. 0 .and. iabsms .eq. 1 ) then

                  ec(ibkec+no,ipomp+1) = emin(ityp)

               end if

*-----------------------------------------------------------------------
*     booking of collision events
*-----------------------------------------------------------------------

 1000 continue

*-----------------------------------------------------------------------
*        jcoll = 0 : path through
*-----------------------------------------------------------------------

            if( jcoll .eq. 0 ) then

                     nclsts = -1

            else if( jcoll .ge. 1 .and. jcoll .le. 4 ) then

               if( nclst .gt. 0 ) then
                     aevts(iaevt+jcoll,ireg) =
     &               aevts(iaevt+jcoll,ireg) + oldwt
                  if( imat .gt. 0 ) then
                     bevts(ibevt+jcoll,imat) =
     &               bevts(ibevt+jcoll,imat) + oldwt
                  end if

                  if( ihvi .ne. 0 ) then

                     jcoll = 5
                     aevts(iaevt+ihvi,ireg) =
     &               aevts(iaevt+ihvi,ireg) + oldwt
                  if( imat .gt. 0 ) then
                     bevts(ibevt+ihvi,imat) =
     &               bevts(ibevt+ihvi,imat) + oldwt
                  end if

                  end if

               end if

            else if( jcoll .ge. 6 .and. jcoll .le. 8 ) then
                     aevts(iaevt+jcoll+4,ireg) =
     &               aevts(iaevt+jcoll+4,ireg) + oldwt
                  if( imat .gt. 0 ) then
                     bevts(ibevt+jcoll+4,imat) =
     &               bevts(ibevt+jcoll+4,imat) + oldwt
                  end if

            else if( jcoll .eq. 9 ) then
                     aevts(iaevt+57,ireg) =
     &               aevts(iaevt+57,ireg) + oldwt
                  if( imat .gt. 0 ) then
                     bevts(ibevt+57,imat) =
     &               bevts(ibevt+57,imat) + oldwt
                  end if

            else if( jcoll .eq. 10 ) then
                     aevts(iaevt+63,ireg) =
     &               aevts(iaevt+63,ireg) + oldwt
                  if( imat .gt. 0 ) then
                     bevts(ibevt+63,imat) =
     &               bevts(ibevt+63,imat) + oldwt
                  end if

            else if( jcoll .eq. 12 .or.
     %               jcoll .eq. 16 .or. jcoll .eq. 17 ) then   ! S.Abe 2017/01/17
               if( imubrmhit .eq. 1 ) then
                     aevts(iaevt+66,ireg) =
     &               aevts(iaevt+66,ireg) + oldwt
                  if( imat .gt. 0 ) then
                     bevts(ibevt+66,imat) =
     &               bevts(ibevt+66,imat) + oldwt
                  end if
               elseif( imuppdhit .eq. 1 ) then
                     aevts(iaevt+67,ireg) =
     &               aevts(iaevt+67,ireg) + oldwt
                  if( imat .gt. 0 ) then
                     bevts(ibevt+67,imat) =
     &               bevts(ibevt+67,imat) + oldwt
                  end if
               elseif( imucapflag .ne. 1 ) then
                     aevts(iaevt+64,ireg) =
     &               aevts(iaevt+64,ireg) + oldwt
                  if( imat .gt. 0 ) then
                     bevts(ibevt+64,imat) =
     &               bevts(ibevt+64,imat) + oldwt
                  end if
               else
                     aevts(iaevt+65,ireg) =
     &               aevts(iaevt+65,ireg) + oldwt
                  if( imat .gt. 0 ) then
                     bevts(ibevt+65,imat) =
     &               bevts(ibevt+65,imat) + oldwt
                  end if
               endif

            else if( jcoll .eq. 18 ) then
                     aevts(iaevt+67+itcoll,ireg) =
     &               aevts(iaevt+67+itcoll,ireg) + oldwt
                  if( imat .gt. 0 ) then
                     bevts(ibevt+67+itcoll,imat) =
     &               bevts(ibevt+67+itcoll,imat) + oldwt
                  end if

            end if

*-----------------------------------------------------------------------

            rncnt(6) = rncnt(6) + 1.0
            call cputime(6)

*-----------------------------------------------------------------------

      return
      end
