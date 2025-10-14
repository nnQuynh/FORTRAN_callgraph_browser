************************************************************************
*                                                                      *
      subroutine dataup(ncol,ncol14)
*                                                                      *
*       booking the particles into datahi and datalo                   *
*       modified by K.Niita on 2011/05/17                              *
*                                                                      *
*           ncol = 13 : projectile disappeared reactions               *
*           ncol = 14 : projectile appeared reactions                  *
*                       except for delta ray production                *
*           ncol14 = 0 : always ncol=13 from repeated collisions       *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*        qclusts(i,nclst) : after nevap or low energy reactions        *
*                                                                      *
*                i = 0, impact parameter                               *
*                  = 1, px (GeV/c)                                     *
*                  = 2, py (GeV/c)                                     *
*                  = 3, pz (GeV/c)                                     *
*                  = 4, etot = sqrt( p**2 + rm**2 ) (GeV)              *
*                  = 5, rest mass (GeV)                                *
*                  = 6, excitation energy (MeV)                        *
*                  = 7, kinetic energy (MeV)                           *
*                  = 8, weight change                                  *
*                  = 9, delay time                                     *
*                  = 10, x-displace                                    *
*                  = 11, y-displace                                    *
*                  = 12, z-displace                                    *
*                                                                      *
*        qclusts(i,nclst) : here changed                               *
*                                                                      *
*                i = 0, impact parameter                               *
*                  = 1, u in lab                                       *
*                  = 2, v in lab                                       *
*                  = 3, w in lab                                       *
*                  = 4, etot = sqrt( p**2 + rm**2 ) (GeV)              *
*                  = 5, rest mass (GeV)                                *
*                  = 6, excitation energy (MeV)                        *
*                  = 7, kinetic energy (MeV)                           *
*                  = 8, real weight                                    *
*                  = 9, real time                                      *
*                  = 10, x coordinate in lab                           *
*                  = 11, y coordinate in lab                           *
*                  = 12, z coordinate in lab                           *
*                                                                      *
************************************************************************
      use MMBANKMOD !FURUTA
      use moddas_region
*-----------------------------------------------------------------------

      implicit real*8(a-h,o-z)

*-----------------------------------------------------------------------

      include 'param00.inc'
      include 'param01.inc'
      include 'param.inc'

*-----------------------------------------------------------------------

      common /bnkmem/ maxbnk, maxbn2, rtrckflp
      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
!$OMP THREADPRIVATE(/icomon/)
      common /eparm/  esmax, esmin, emin(20)
      common /ndemax/ dnmax(20)
      common /tparm/  tmax(20)
      common /jcomon/ nabov,nobch,nocas,nomax
!$OMP THREADPRIVATE(/jcomon/)
      common /anglth/ cosphi,costh,sinphi,sinth
!$OMP THREADPRIVATE(/anglth/)
      common /mathzn/ mathz, mathn, jcoll, kcoll
!$OMP THREADPRIVATE(/mathzn/)

      common /trstar/ itrstar
!$OMP THREADPRIVATE(/trstar/)

      common /pnsave/ egs, uus, vvs, wws, wts, tms, nms, nct(3)
!$OMP THREADPRIVATE(/pnsave/)
      common /tlgeom/ iblz1,iblz2
!$OMP THREADPRIVATE(/tlgeom/)
      common /tlglat/ ilev1,ilev2,ilat1(5,10),ilat2(5,10)
!$OMP THREADPRIVATE(/tlglat/)

*-----------------------------------------------------------------------

      common /inout/  in,io
      common /neulo/  reutn

*-----------------------------------------------------------------------

      common /clustt/ nclsts, iclusts(nnn)
!$OMP THREADPRIVATE(/clustt/)
      common /clustw/ jclusts(0:8,nnn),  qclusts(0:12,nnn)
!$OMP THREADPRIVATE(/clustw/)
      common /otheid/ idpat(20), idoth(200), idono
      common /cntcls/ jcount(3,nnn)
!$OMP THREADPRIVATE(/cntcls/)

      common /cntreg/ ncntc(3), ncreg(3), incrc(3), incrt(3), kcont(3)
      common /cntpat/ icpan(3), icpat(3,20,2)
      common /tmtreg/ ntmrg, intmc, intmt, ktime
      common /fsmul/  ridfs, iidfs

*-----------------------------------------------------------------------

      dimension     idas(1)
      equivalence ( das, idas )

*-----------------------------------------------------------------------

      dimension data(11)
      dimension jcnt(3)
      dimension lcnt(22)   ! S.Abe 2018/02/26, 20->21, S.H. 21->22 (2022.2.28)

*-----------------------------------------------------------------------

      integer iegsemi, iegsout
      common /egsemi/ iegsemi, iegsout
      common /egs5cmn8/iegsxray
!$OMP THREADPRIVATE(/egs5cmn8/)

*-----------------------------------------------------------------------
      common /pointt/ itpont

*-----------------------------------------------------------------------
      common /trskip/ ntrskip
!$OMP THREADPRIVATE(/trskip/)
      common /tstara/ atmrc(10,8)      ! Ogawa 2023/7/18, (10,7) -> (10,8) for plasmon
!$OMP THREADPRIVATE(/tstara/)
      common /tpdcta/ nflumu
!$OMP THREADPRIVATE(/tpdcta/)

      common /regdc/  idrg(kvlmax), idgr(kvmmax)
      common /  tscmsg   / ktsc(kvlmax), mntsc, ntsc(kvlmax)
      common /  tscreg   / ntscell(kvlmax)
      common / etsminmax / etsmin, etsmax
      common / ptsminmax / ptsmin, ptsmax
      common / ctsminmax / ctsmin, ctsmax
      common / tsminmax  / tsmax

      common /etsion/ icoll
!$OMP THREADPRIVATE(/etsion/)
      common / ele2ndc / e2ndc(2),u2ndc(2),v2ndc(2),w2ndc(2),icollc
!$OMP THREADPRIVATE(/ele2ndc/)


      common /cntmx/ ncntmx(3), nctmxsr(3)
!$OMP THREADPRIVATE(/cntmx/)

*-----------------------------------------------------------------------

         idcol = 0
         ncol  = 13

*-----------------------------------------------------------------------
         itrstar = 0

         ntrskip = 0    ! S.Abe 2017/02/08

         if( nclsts .le. 0 ) return

*-----------------------------------------------------------------------


*-----------------------------------------------------------------------
*        timer after collision
*-----------------------------------------------------------------------

            if( ntmrg .gt. 0 ) then

                        idsm = intmc
                        jdsm = 0

                        kdsm = ktime
                        ldsm = 0

                  do m = 1, ntmrg

                        jdsm = jdsm + 1
                        ntrn = idas_intmc(idsm+jdsm)
                        jdsm = jdsm + 1
                        mtrn = idas_intmc(idsm+jdsm)

                        ldsm  = ldsm + 1
                        inin  = idas_ktime(kdsm+ldsm)
                        ldsm  = ldsm + 1
                        inout = idas_ktime(kdsm+ldsm)
                        ldsm  = ldsm + 1
                        incol = idas_ktime(kdsm+ldsm)
                        ldsm  = ldsm + 1
                        inref = idas_ktime(kdsm+ldsm)

                     if( incol .ne. 0 ) then

                           jj  = 0

                        do ii = 1, ntrn

                           call tregck(iblz1,ilev1,ilat1,mtrn
     &                                 ,idas_intmc(idsm+jdsm+1),jj,icc)

                           if( icc .ne. 0 ) then

                              if( incol .eq. 1 ) then

                                 tc(ibktc+no,ipomp+1) =
     &                                  -abs(tc(ibktc+no,ipomp+1))

                              else if( incol .eq. -1 ) then

                                 tc(ibktc+no,ipomp+1) = 0.d0

                              end if

                           end if

                        end do

                     end if

                        jdsm = jdsm + mtrn

                  end do

            end if

*-----------------------------------------------------------------------
*           counter after collision
*-----------------------------------------------------------------------

                  jcnt(1) = ncnt(ibknct+1,no,ipomp+1)
                  jcnt(2) = ncnt(ibknct+2,no,ipomp+1)
                  jcnt(3) = ncnt(ibknct+3,no,ipomp+1)

      if( ncntc(1) .eq. 1 .or. ncntc(2) .eq. 1 .or.
     &    ncntc(3) .eq. 1 ) then

         do k = 1, 3
            if( ncntc(k) .eq. 1 ) then

               knn  = ityp
               kcg  = jtyp
               kkf  = ktyp

               call pcchck(k,knn,kkf,kcg,icpan,icpat,icc)

               if( icc .eq. 1 ) then

                  idsm = incrt(k)
                  jdsm = 0

                  kdsm = kcont(k)
                  ldsm = 0

                  do m = 1, ncreg(k)

                     jj  = 0

                     jdsm = jdsm + 1
                     ntrn = idas_incrt(idsm+jdsm)
                     jdsm = jdsm + 1
                     mtrn = idas_incrt(idsm+jdsm)

                     do ll = 1, maxcntr-1    ! S.H. set maxcntr-1(reg case) (2022.3.24)
                        lcnt(ll) = 0
                     enddo

                     sumatmrc = 0.d0
                     do ia = 1, 5
                     do ja = 1, 5    ! except mscat
                        sumatmrc = sumatmrc + atmrc(ia,ja)
                     enddo
                     enddo
                     sumatmrc = sumatmrc + atmrc(1,6)   ! S.Abe 2018/10/03, add rayleigh scattering

                     ldsm  = ldsm + 1
                     lcnt(1) = idas_kcont(kdsm+ldsm)                 ! in
                     ldsm  = ldsm + 1
                     lcnt(2) = idas_kcont(kdsm+ldsm)                 ! out
                     ldsm  = ldsm + 1
                     if( ( jcoll .ne. 7 .and. jcoll .ne. 13 .and.
     &                     jcoll .ne. 8 .and. jcoll .ne. 14 ) .or.
     &                   ( ( jcoll .eq. 7 .or. jcoll .eq. 13 .or.
     &                       jcoll .eq. 8 .or. jcoll .eq. 14 ) .and.
     &                     sumatmrc .gt. 0.d0 ) )
     &               lcnt(3) = idas_kcont(kdsm+ldsm)                 ! collision
                     ldsm  = ldsm + 1
                     lcnt(4) = idas_kcont(kdsm+ldsm)                 ! refrect
                     ldsm  = ldsm + 1
                     if( kcoll .eq. 1  .or. kcoll .eq. 5 )
     &                lcnt(5) = idas_kcont(kdsm+ldsm)                ! fission
                     ldsm  = ldsm + 1
                     if( jcoll .eq. 3  .or. kcoll .eq. 3 )
     &                lcnt(6) = idas_kcont(kdsm+ldsm)                ! elastic
                     ldsm  = ldsm + 1
                     if( ( jcoll .eq. 1  .or.
     &                     jcoll .eq. 4  .or. jcoll .eq. 5  .or.
     &                     jcoll .eq. 10 .or. jcoll .eq. 15 .or.
     &                     jcoll .eq. 16 .or. jcoll .eq. 17 ) .and.
     &                   kcoll .ne. 3 )
     &                lcnt(7) = idas_kcont(kdsm+ldsm)                ! nonelastic
                     ldsm  = ldsm + 1
                     if( jcoll .eq. 1  .or. jcoll .eq. 3  .or.
     &                   jcoll .eq. 4  .or. jcoll .eq. 5  .or.
     &                   jcoll .eq. 10 .or. jcoll .eq. 15 .or.
     &                   jcoll .eq. 16 .or. jcoll .eq. 17 )
     &                lcnt(8) = idas_kcont(kdsm+ldsm)                ! nuclear reaction
                     ldsm  = ldsm + 1
                     if( jcoll .eq. 2 .and. ityp .ne. 13 )
     &                lcnt(9) = idas_kcont(kdsm+ldsm)                ! decay
                     ldsm  = ldsm + 1
                     if( jcoll .eq. 11 .or. sumatmrc .gt. 0.d0 )
     &                lcnt(10) = idas_kcont(kdsm+ldsm)                ! atomic reaction
                     ldsm  = ldsm + 1
                     if( jcoll .eq. 11 )
     &                lcnt(11) = idas_kcont(kdsm+ldsm)               ! delta ray production
                     ldsm  = ldsm + 1
                     if( atmrc(1,1).gt.0.d0 .or. atmrc(2,1).gt.0.d0 .or.
     &                   atmrc(3,1).gt.0.d0 .or. atmrc(4,1).gt.0.d0 )
     &                lcnt(12) = idas_kcont(kdsm+ldsm)               ! atomic fluorescence
                     ldsm  = ldsm + 1
                     if( atmrc(1,2).gt.0.d0 .or. atmrc(2,2).gt.0.d0 .or.
     &                   atmrc(3,2).gt.0.d0 )
     &                lcnt(13) = idas_kcont(kdsm+ldsm)               ! auger electron production
                     ldsm  = ldsm + 1
                     if( atmrc(2,3).gt.0.d0 .or. atmrc(3,3).gt.0.d0 .or.
     &                   atmrc(4,3).gt.0.d0 .or. atmrc(5,3).gt.0.d0 .or.
     &                   atmrc(2,5).gt.0.d0 )
     &                lcnt(14) = idas_kcont(kdsm+ldsm)               ! bremsstrahlung
                     ldsm  = ldsm + 1
                     if( atmrc(1,3).gt.0.d0 )
     &                lcnt(15) = idas_kcont(kdsm+ldsm)               ! photoelectric effect
                     ldsm  = ldsm + 1
                     if( atmrc(1,4).gt.0.d0 )
     &                lcnt(16) = idas_kcont(kdsm+ldsm)               ! compton scattering
                     ldsm  = ldsm + 1
                     if( atmrc(1,5).gt.0.d0 .or. atmrc(4,5).gt.0.d0 .or.
     &                   atmrc(5,5).gt.0.d0 )
     &                lcnt(17) = idas_kcont(kdsm+ldsm)               ! pair production
                     ldsm  = ldsm + 1
                     if( atmrc(3,5).gt.0.d0 )
     &                lcnt(18) = idas_kcont(kdsm+ldsm)               ! positron annihilation
                     ldsm  = ldsm + 1
                     if( atmrc(2,6).gt.0.d0 .or. atmrc(3,6).gt.0.d0 )
     &                lcnt(19) = idas_kcont(kdsm+ldsm)               ! mscat
                     ldsm  = ldsm + 1
                     if( atmrc(1,6).gt.0.d0 )
     &                lcnt(20) = idas_kcont(kdsm+ldsm)               ! rayleigh scattering
                     ldsm  = ldsm + 1
                     if( atmrc(2,4).gt.0.d0 .or. atmrc(3,4).gt.0.d0 )
     &                lcnt(21) = idas_kcont(kdsm+ldsm)               ! knock-on electron prodcution
                     ldsm  = ldsm + 1
                     lcnt(22) = idas_kcont(kdsm+ldsm)                ! ndata

                     do ii = 1, ntrn

                        call tregck(iblz1,ilev1,ilat1,
     &                              mtrn,idas_incrt(idsm+jdsm+1),jj,icc)

                        if( icc .ne. 0 ) then
                           do ll = 1, maxcntr-1   ! S.H. set maxcntr-1(reg case) (2022.3.24)
                              if( ll .eq. 1 .or. ll .eq. 2 .or.
     &                            ll .eq. 4 .or. ll .eq. 22 ) cycle
                              if( lcnt(ll) .eq. 10000 ) jcnt(k) = 0
                           enddo

                           jcntbk = jcnt(k)
                           do ll = 1, maxcntr-1   ! S.H. set maxcntr-1(reg case) (2022.3.24)
                              if( ll .eq. 1 .or. ll .eq. 2 .or.
     &                            ll .eq. 4 .or. ll .eq. 22 ) cycle
                              if( abs(lcnt(ll)) .le. 9999 )
     &                         jcnt(k) = jcnt(k) + lcnt(ll)
                           enddo
                           if( abs(jcnt(k)) .gt. 9999 )
     &                      jcnt(k) = jcntbk

                        endif

                        jdsm = jdsm + mtrn

                     end do

                  end do

               end if

            end if

         end do

      end if

*-----------------------------------------------------------------------
                   nmpro = 0
                   nmnut = 0

*-----------------------------------------------------------------------

      do i = 1, nclsts

*-----------------------------------------------------------------------

                if( iclusts(i) .eq. 1 ) nmpro = nmpro + 1
                if( iclusts(i) .eq. 2 ) nmnut = nmnut + 1

*-----------------------------------------------------------------------
*        energy, angle and angle transform to lab system
*-----------------------------------------------------------------------

                  xmd = qclusts(10,i)
                  ymd = qclusts(11,i)
                  zmd = qclusts(12,i)

                  qclusts(7,i)  = max( 0.0d0, qclusts(7,i) )
                  qclusts(8,i)  = wt(ibkwt+no,ipomp+1) * qclusts(8,i)

               if( tc(ibktc+no,ipomp+1) .gt. 0 ) then
                  qclusts(9,i)  = tc(ibktc+no,ipomp+1) + qclusts(9,i)
               else
                  qclusts(9,i)  = tc(ibktc+no,ipomp+1)
               end if

                  qclusts(10,i) = xc(ibkxc+no,ipomp+1) + qclusts(10,i)
                  qclusts(11,i) = yc(ibkyc+no,ipomp+1) + qclusts(11,i)
                  qclusts(12,i) = zc(ibkzc+no,ipomp+1) + qclusts(12,i)

                  pab = sqrt( qclusts(1,i)**2
     &                      + qclusts(2,i)**2
     &                      + qclusts(3,i)**2 )

               if( pab .gt. 0.0d0 ) then

                  adc = qclusts(1,i) / pab
                  bdc = qclusts(2,i) / pab
                  gdc = qclusts(3,i) / pab

               else

                  adc = 0.0
                  bdc = 0.0
                  gdc = 1.0

               end if

                  t1   = costh  * adc + sinth  * gdc
                  adcc = cosphi * t1  - sinphi * bdc
                  bdcc = sinphi * t1  + cosphi * bdc
                  gdcc = costh  * gdc - sinth  * adc

                  qclusts(1,i) = adcc
                  qclusts(2,i) = bdcc
                  qclusts(3,i) = gdcc

*-----------------------------------------------------------------------
*           counter after collision (all ejectiles are the same counter)
*-----------------------------------------------------------------------

                  jcount(1,i) = jcnt(1)
                  jcount(2,i) = jcnt(2)
                  jcount(3,i) = jcnt(3)

*-----------------------------------------------------------------------
*           particle id and decay particle or not
*-----------------------------------------------------------------------

               knn  = jclusts(3,i)
               kf   = jclusts(7,i)
               idcy = idcyc(kf)

cABE 2022/02/24, keep decay flag for positron with the energies below cutoff
               if( kf .eq. -11 .and. jclusts(4,i) .ne. 2 ) idcy = 0
               if( jclusts(4,i) .eq. 2 ) jclusts(4,i) = 0

*-----------------------------------------------------------------------
*           minimum energy per particle
*-----------------------------------------------------------------------

            if( knn .lt. 15 ) then

             do i1=1,mntsc
             if(ntsc(i1).eq.iblz(ibkblz+no,ipomp+1))then
                mat1 = ntscell(idgr(ntsc(i1)))
                goto 1112
             endif
             enddo
                mat1 = 0
 1112 continue

             if( ( knn .eq. 12 .or. knn .eq. 13 ) .and.
     .             mat1 .ne. 0                    .and.
     .           ( qclusts(7,i) .le. etsmax )     .and.
     .           ( qclusts(7,i) .ge. etsmin )     )then
               emint = etsmin
             else
               emint = emin(knn)
             endif

            else

               emint = emin(knn) * dble( ibryf(knn,kf) )

            end if

*-----------------------------------------------------------------------
*           knn = 11 : kf-code and booking
*-----------------------------------------------------------------------

                     kid = 0

               if( knn .eq. 11 ) then

                  do k = 1, idono

                     if( kf .eq. idoth(k) ) then
                        kid = k
                        exit
                     endif

                  end do

               end if

*-----------------------------------------------------------------------
*           neutron below dnmax(2) or emin(2)
*-----------------------------------------------------------------------

               if( knn .eq. 2 .and. jcoll .le. 5 .and.
     &             qclusts(7,i) .le. max(emin(2),dnmax(2)) ) then
!$OMP ATOMIC
                    reutn = reutn + 1.0d0

               end if

*-----------------------------------------------------------------------
*           iidfs : induced fission option : use spx variable !!
*-----------------------------------------------------------------------

               if( iidfs .eq. 1 .and. ridfs .gt. 0.0d0 ) then

                  spx(ibkspx+no,ipomp+1) = spx(ibkspx+no,ipomp+1)+1.0d0

                  ridfs = -1.0d0

               end if

*-----------------------------------------------------------------------
*        jelas option for t-yield, itratar
*-----------------------------------------------------------------------

         if( mathz .eq. jclusts(1,i) .and.
     &       mathn .eq. jclusts(2,i) )  itrstar = 1

*-----------------------------------------------------------------------
*        exclude non-produced particle for t-product,
*        i.e., EGS capbind electron and annihilation positron
*-----------------------------------------------------------------------

         if( jclusts(4,i) .eq. -11 .or.
     &       ( atmrc(3,5) .gt. 0.d0 .and. jclusts(4,i) .eq. -1 ) ) then
            ntrskip = i
         endif

*-----------------------------------------------------------------------
*        datahi booking for transprot particles
*        and energy, time and weight check
*        if decay particle, put in datahi,
*        photon : check photon from electron : name > 1
*-----------------------------------------------------------------------
         if( ( ( knn .ne. 11 .and. idpat(knn) .eq. 1 ) .or.
     &         ( knn .eq. 11 .and. idpat(knn) .eq. 1 .and.
     &           ( ( 900000.le.abs(kf) .and. abs(kf).le.999999 ) .or.
     &             kid .ne. 0 ) ) )
     &          .and.   qclusts(8,i) .gt. 0.0
     &          .and.   abs(qclusts(9,i)) .lt. tmax(knn)
     &          .and. ( qclusts(7,i) .gt. emint .or. idcy .eq. 1 )
     &          .and. ( jclusts(4,i) .ne. -11 )
     &      ) then

*-----------------------------------------------------------------------
*           ncol = 14 : particle for sequential run
*-----------------------------------------------------------------------
cc H.Iwase 2014/1/9 (excluding photoelectric X-rays)

cKN 2015/03/24 for point estimator !!!

            if( idcol .eq. 0     .and.
cKN 2019/10/20 for repeated collisions
     &          ncol14 .ne. 0    .and.
     &          xmd   .eq. 0.0d0 .and.
     &          ymd   .eq. 0.0d0 .and.
     &          zmd   .eq. 0.0d0 .and.
cKN 2016/05/11 exclude decay event
     &          jcoll .ne. 2     .and.
     &          jcoll .ne. 11    .and.
     &          kf    .eq. ktyp  .and.
     &         itpont .eq. 0     .and.
     &       ( ( iegsemi .ne. 0 .and. iegsxray .ne. 1 ) .or.
     &           iegsemi .eq. 0 ) ) then

                  ncol = 14

                  uus = qclusts(1,i)
                  vvs = qclusts(2,i)
                  wws = qclusts(3,i)
                  egs = qclusts(7,i)
                  wts = qclusts(8,i)
                  tms = qclusts(9,i)

                  nms = name(ibknam+no,ipomp+1)
                  if( name(ibknam+no,ipomp+1) .lt. 999 .and.
     &                 jcoll .ne. 11 .and.
     &              ( ( ityp .ne. 12 .and. ityp .ne. 13 ) .or.
     &                ( knn .ne. 12 .and. knn .ne. 13 ) .or.
     &                  name(ibknam+no,ipomp+1) .ne. 1 ) )
     &            nms = name(ibknam+no,ipomp+1) + 1

                  nct(1) = jcount(1,i)
                  nct(2) = jcount(2,i)
                  nct(3) = jcount(3,i)

                  idcol = idcol + 1

                  do itt = 1, 3
                    if ( nct(itt).gt.ncntmx(itt) ) then
                      ncntmx(itt) = nct(itt)
                    end if
                  end do

*-----------------------------------------------------------------------
*              photon heating flag : jclusts(4,nclsts) = 1
*-----------------------------------------------------------------------
               if( knn .eq. 14 .and.
     &             ( ( jcoll .ne. 7 .and. jcoll .ne. 15 .and.      ! S.Abe 2017/01/17
     &                 jcoll .ne. 8 ) .or.
     &               ( ( jcoll .eq. 7 .or. jcoll .eq. 15 ) .and.   ! S.Abe 2017/01/17
     &                 jclusts(4,i) .eq. 1 .and.
     &                 name(ibknam+no,ipomp+1) .eq. 1 ) .or.
     &               ( jcoll .eq. 8 .and.
     &                 name(ibknam+no,ipomp+1) .eq. 1 ) ) ) then

                  nms = 1

               end if

*-----------------------------------------------------------------------
*           data up in bank
*-----------------------------------------------------------------------

            else

                  nabov = nabov + 1

                  if( nabov .gt. maxbn2 ) goto 9010

                  ntya(iaknty+nabov,ipomp+1) = jclusts(3,i)
                  nkfa(iaknkf+nabov,ipomp+1) = jclusts(7,i)

                  ea(iake+nabov,ipomp+1)   = qclusts(7,i)


                  ta(iakt+nabov,ipomp+1)   = qclusts(9,i)
                  wta(iakwt+nabov,ipomp+1) = qclusts(8,i)

                  ua(iaku+nabov,ipomp+1) = qclusts(1,i)
                  va(iakv+nabov,ipomp+1) = qclusts(2,i)
                  wa(iakw+nabov,ipomp+1) = qclusts(3,i)

                  xa(iakx+nabov,ipomp+1) = qclusts(10,i)
                  ya(iaky+nabov,ipomp+1) = qclusts(11,i)
                  za(iakz+nabov,ipomp+1) = qclusts(12,i)

                  iblza(iakblz+nabov,ipomp+1) = iblz(ibkblz+no,ipomp+1)
                  nmeda(iaknmd+nabov,ipomp+1) = nmed(ibknmd+no,ipomp+1)
                  wtina(iakwin+nabov,ipomp+1) = wtin(ibkwin+no,ipomp+1)
                  wtnza(iakwnz+nabov,ipomp+1) = wtnz(ibkwnz+no,ipomp+1)

                  spxa(iakspx+nabov,ipomp+1) = spx(ibkspx+no,ipomp+1)
                  spya(iakspy+nabov,ipomp+1) = spy(ibkspy+no,ipomp+1)
                  spza(iakspz+nabov,ipomp+1) = spz(ibkspz+no,ipomp+1)

                  nzsta(iakzst+nabov,ipomp+1) = jclusts(5,i)
                  nsosa(iaksos+nabov,ipomp+1) = nsos(ibksos+no,ipomp+1)

                  nfcsa(iaknfc+nabov,ipomp+1) = 0

                  namea(iaknam+nabov,ipomp+1) = name(ibknam+no,ipomp+1)
                  itetposa(ibtetposa+nabov,ipomp+1) =
     &                                      itetpos(ibtetpos+no,ipomp+1)
                  if( name(ibknam+no,ipomp+1) .lt. 999 .and.
     &                  jcoll .ne. 11 .and.
     &              ( ( ityp .ne. 12 .and. ityp .ne. 13 ) .or.
     &                ( knn .ne. 12 .and. knn .ne. 13 ) .or.
     &                  name(ibknam+no,ipomp+1) .ne. 1 ) )
     &            namea(iaknam+nabov,ipomp+1) =
     &                              name(ibknam+no,ipomp+1) + 1

                  ncnta(iaknct+1,nabov,ipomp+1) = jcount(1,i)
                  ncnta(iaknct+2,nabov,ipomp+1) = jcount(2,i)
                  ncnta(iaknct+3,nabov,ipomp+1) = jcount(3,i)

*-----------------------------------------------------------------------
*              photon heating flag :
*-----------------------------------------------------------------------
               if( knn .eq. 14 .and.
     &             ( ( jcoll .ne. 7 .and. jcoll .ne. 15 .and.      ! S.Abe 2017/01/17
     &                 jcoll .ne. 8 ) .or.
     &               ( ( jcoll .eq. 7 .or. jcoll .eq. 15 ) .and.   ! S.Abe 2017/01/17
     &                 jclusts(4,i) .eq. 1 .and.
     &                 name(ibknam+no,ipomp+1) .eq. 1 ) .or.
     &               ( jcoll .eq. 8 .and.
     &                 name(ibknam+no,ipomp+1) .eq. 1 ) ) ) then

                  namea(iaknam+nabov,ipomp+1) = 1

               end if

*-----------------------------------------------------------------------

            end if

*-----------------------------------------------------------------------
*        stop or dead particles jclusts(4,i) = -1
*-----------------------------------------------------------------------

         else

               jclusts(4,i) = -1

         end if

*-----------------------------------------------------------------------

      end do

*-----------------------------------------------------------------------

         if( itrstar .eq. 0 ) then

             if( jcoll .eq. 6 .and.
     &         ( kcoll .eq. 3 .or.
     &         ( nmpro .eq. 0 .and. nmnut .eq. 1 ) ) ) then
                 itrstar = 1
             end if

             if( jcoll .eq. 9 .and.
     &         ( kcoll .eq. 3 .or.
     &         ( nmpro .eq. 1 .and. nmnut .eq. 0 ) ) ) then
                 itrstar = 1
             end if

          end if

*-----------------------------------------------------------------------

      return

*-----------------------------------------------------------------------
*     error message
*-----------------------------------------------------------------------

 9010 continue

         write(io,12) nabov,maxbn2,nocas
   12    format(/' *** error messge from s.dataup ***'
     &   /' storage size (nabov) of particles is greater ',
     &    'than maxbn2.'
     &   /' nabov  =',i10
     &   /' maxbn2 =',i10
     &   /' nocas  =',i10)
         call parastop( 837 )

*-----------------------------------------------------------------------

      end


