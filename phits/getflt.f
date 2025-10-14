************************************************************************
*                                                                      *
      subroutine getflt(fpl)
*                                                                      *
*       flight path length sampling                                    *
*       last modified by K.Niita on 2011/05/12                         *
*       for proton, total cross section is the max value.              *
*                                                                      *
*     output:                                                          *
*                                                                      *
*       fpl     : sampling flight path length (cm)                     *
*                                                                      *
*       totttl  : total flight length (1/cm)                           *
*       tothyd  : Hydrogen flight length (1/cm)                        *
*       totdcy  : decay fligt length (1/cm)                            *
*       totdel  : delta ray collision length (1/cm)                    *
*       tottsm  : track-structure collision length (1/cm)              *
*                                                                      *
*       siggn(ksign+lem) : nonelastic flight length for (lem,mat)      *
*       sigge(ksige+lem) : elastic flight length for (lem,mat)         *
*                                                                      *
************************************************************************
      use MMBANKMOD !FURUTA
      use MEMBANKMOD !FURUTA
      use ELEDATAMOD, only : ichem
      use neutrino_mod, only : neutrino_Xsec
      use ion_track_structure, only : pts_flt2
      use ets_art, only : etsflt_art
      use udm_Parameter
      use moddas_material
      implicit real*8(a-h,o-z)

*-----------------------------------------------------------------------

      include 'param00.inc'
      include 'param.inc'

      parameter ( rpmass = 938.27, rnmass = 939.58 )
      parameter ( nemd = 20 )
      parameter ( rlit = 2.997925d+1 )    ! S.Abe 2015/08/10, clight [cm/ns]

*-----------------------------------------------------------------------

      common /paraj/ mstz(300), parz(300) ! T.Sato 2020/09/12

      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
!$OMP THREADPRIVATE(/icomon/)

      common /regdc/  idrg(kvlmax), idgr(kvmmax)
      common /kmat1d/ idmn(0:kvlmax), idnm(kvmmax)
      common /tlgeom/ iblz1,iblz2
!$OMP THREADPRIVATE(/tlgeom/)
      common /mathzn/ mathz, mathn, jcoll, kcoll
!$OMP THREADPRIVATE(/mathzn/)
      common /qparm/  ielas,icasc,iqstep,lvlopt,igamma
      common /xtotal/ totttl,tothyd,totdcy,totdel,tottsm
!$OMP THREADPRIVATE(/xtotal/)
      common /eparm/  esmax, esmin, emin(20)
      common /ndemax/ dnmax(20)
      common /rcmini/ rcmin(20)
      common /emode/  emodem, ge1, ge2, iemode

      common /celdn/  denr(kvlmax), denm(kvlmax), denc(kvlmax)
      common /regcm/  icmg(kvlmax)

      common /delesg/ rdels(kvlmax), mndel, ndels(kvlmax)
      common /delreg/ delm(kvlmax), kdelt

      common /kmat1g/ kmat(kvlmax)
      common /xgeosm/ ksig(kvlmax)
      common /xinels/ ksige, ksign


*-----------------------------------------------------------------------

      common /clustt/ nclsts, iclusts(nnn)
!$OMP THREADPRIVATE(/clustt/)
      common /clustw/ jclusts(0:8,nnn),  qclusts(0:12,nnn)
!$OMP THREADPRIVATE(/clustw/)

      dimension crsen(kvlmax), crstl(kvlmax)

*-----------------------------------------------------------------------

      common /muint/ imuint, imubrm, imuppd, imucap
      common /sigmuh/ sigmuv(200),sigmub(200),sigmup(200),
     &                nmuh,itzmuh(200),itamuh(200)
!$OMP THREADPRIVATE(/sigmuh/)

*-----------------------------------------------------------------------

      data rcc / 2.997925d+10 /
      data pi / 3.1415926535897932d0 /

*-----------------------------------------------------------------------

      integer iegsemi, iegsout
      common /egsemi/ iegsemi, iegsout

      common /pnint/ ipnint
      common /gmuppd/ igmuppd ! S.Abe 2019/11/07
      common /xphotn/ totgam,totegs,totpni,totgmm ! y.sakaki 2019/09/26
!$OMP THREADPRIVATE(/xphotn/)

      common /tdcyl/ tdcylife, itdcyw
!$OMP THREADPRIVATE(/tdcyl/)

      common /kmat1k/ kmatd(kvlmax), kmate(kvlmax)
      common /kmat1ka/ kmathd(kvlmax), kmathe(kvlmax)

*-----------------------------------------------------------------------

      common /ccxsm/  icxsni, icxspi

*-----------------------------------------------------------------------
      common /fgdata/ ifrgd, ifgm, ifgdf(kvlmax,5), frgfl(kvlmax)
      character frgfl*200

      common /  tscmsg   / ktsc(kvlmax), mntsc, ntsc(kvlmax)
      common /  tscreg   / ntscell(kvlmax)
      common / etsminmax / etsmin, etsmax
      common / ptsminmax / ptsmin, ptsmax
      common / ctsminmax / ctsmin, ctsmax
      common / tsminmax  / tsmax
      common /  tsxcl    / itsxcl ! Ogawa 2023/06/02

      common /scinswt/iswitch
!$OMP THREADPRIVATE(/scinswt/)
      iswitch = 0

*-----------------------------------------------------------------------
*        initial values and constants
*-----------------------------------------------------------------------
               fpl    = 1.0d+20

               totttl = 0.d0
               tothyd = 0.d0
               totdcy = 0.d0
               totdel = 0.d0
               tottsm = 0.d0

               totgam = 0.d0
               totegs = 0.d0
               totpni = 0.d0
               totgmm = 0.d0 ! y.sakaki 2019/09/26

               ein    = e(ibke+no,ipomp+1)
               kf     = nkf(ibknkf+no,ipomp+1)

               tme    = abs(t(ibkt+no,ipomp+1)) / 10.0d0
               icl    = idgr(iblz1)

               ieletr = 0

               tdcylife = 0.d0   ! S.Abe 2015/08/10

               tottnd = 0d0 ! initialization (S.Hashimoto, 2016.6.29)

*-----------------------------------------------------------------------
*        forced collision or reaction cutoff option
*-----------------------------------------------------------------------

            if( nfcs(ibknfc+no,ipomp+1) .eq. 1 ) return

            emint = rcmin(ityp) * dble(max(1,mtyp))

*-----------------------------------------------------------------------
*        totdcy : decay flight legth
*        decay particle is considered even if emin is less than cmin
*-----------------------------------------------------------------------

            totdcy = dcyfl(kf,ein,rtyp)

            if( totdcy .gt. 1.0d+10 ) then

               fpl = 1.0d-10
               return

            end if

*-----------------------------------------------------------------------
*        Negative muon capture in material
*-----------------------------------------------------------------------

         if( ityp .eq. 7 .and. imucap .gt. 0 .and. mat .gt. 0 ) then

            lem   = nint( dnel_das(kmat0+mat) )
            hydro = denh_das(kmat0+mat)

            totdcy = 0.d0
            dentot = 0.d0

            do i = 1, lem

               itz = idnint( zz_das(kmat(mat)+i) )
               ita = idnint( a_das(kmat(mat)+i) )
               capr = dmax1(0.d0, caprmu(itz,ita)*1.d3)
     &                *den_das(kmat(mat)+i)
               dcyr = dmax1(0.d0, dcyrmu(itz)*1.d3)
     &                *den_das(kmat(mat)+i)
               dentot = dentot + den_das(kmat(mat)+i)
               totdcy = totdcy + capr + dcyr

            enddo

            if( hydro .gt. 0.d0 ) then

               itz = 1
               ita = 1

               capr = dmax1(0.d0, caprmu(itz,ita)*1.d3) * hydro
               dcyr = dmax1(0.d0, dcyrmu(itz)*1.d3) * hydro

               dentot = dentot + hydro
               totdcy = totdcy + capr + dcyr

            endif

            totdcy = totdcy / dentot
     &              / rcc / dsqrt(((ein/rtyp)+2.d0)*(ein/rtyp))

         endif

*-----------------------------------------------------------------------
*        tottsm : track structure flight length     moved to here 2021/11/09 T.Ogawa
*-----------------------------------------------------------------------

cc KURBUC
         if( mntsc .gt. 0) then
          sig_macro = 0.d0
          if    (ntscell( idgr(iblz(ibkblz+no,ipomp+1)) ) .ne. 0)then
            if    (ntscell( idgr(iblz(ibkblz+no,ipomp+1)) ) .eq. 1 .and.
     &             ityp .eq. 1                    .and.
     &             e(ibke+no,ipomp+1) .le. ptsmax .and.
     &             e(ibke+no,ipomp+1) .ge. ptsmin ) then

              call pts_flt(sig_macro)
              if(sig_macro .gt. 0.d0) then ! .and. itsxcl .ge. 1) then ! Ogawa 2023/06/02. KURBUC exclusive mode without nuclear reactions
                  tottsm = tottsm + sig_macro
                  goto 100
              endif

            elseif(ntscell( idgr(iblz(ibkblz+no,ipomp+1)) ) .eq. 1 .and.
     &             ktyp .eq. 6000012               .and.
     &             e(ibke+no,ipomp+1) .le. ctsmax  .and.
     &             e(ibke+no,ipomp+1) .ge. ctsmin ) then

              call cts_flt(sig_macro)
              if(sig_macro .gt. 0.d0) then ! .and. itsxcl .ge. 1) then ! Ogawa 2023/06/02. KURBUC exclusive mode without nuclear reactions.
                  tottsm = tottsm + sig_macro
                  goto 100
              endif

c hirata 20220602 for etsmode02
        elseif(ntscell( idgr(iblz(ibkblz+no,ipomp+1)) ) .eq. -1 .and.
     &              ichem(mat,1) .eq. 22014 .and.
     &              ichem(mat,2) .eq. 0     .and.
     &              (ityp .eq. 12 .or. ityp .eq. 13 ) .and.
     &              e(ibke+no,ipomp+1) .le. etsmax   .and.
     &              e(ibke+no,ipomp+1) .ge. etsmin ) then

               call etsflt02(sig_macro)  ! for silicon track structure

        elseif(ntscell( idgr(iblz(ibkblz+no,ipomp+1)) ) .eq. -1 .and.
     &              ichem(mat,1) .eq. 1 .and.   ! water-ETS can be called with mID=-1 and chem=H2O
     &              ichem(mat,2) .eq. 0     .and.
     &              (ityp .eq. 12 .or. ityp .eq. 13 ) .and.
     &              e(ibke+no,ipomp+1) .le. etsmax   .and.
     &              e(ibke+no,ipomp+1) .ge. etsmin ) then

               call etsflt(sig_macro)  ! for water track structure
c hirata 20220829 for ETSART
        elseif(ntscell( idgr(iblz(ibkblz+no,ipomp+1)) ) .eq. -1 .and.
     &              (ityp .eq. 12 .or. ityp .eq. 13 ) .and.
     &              e(ibke+no,ipomp+1) .le. etsmax   .and.
     &              e(ibke+no,ipomp+1) .ge. etsmin ) then

               call etsflt_art(sig_macro)  ! for Electron trac structure for arbitary mode

            else
             if   ((ityp .eq. 12 .or. ityp .eq. 13 ) .and.
     &              e(ibke+no,ipomp+1) .le. etsmax   .and.
     &              e(ibke+no,ipomp+1) .ge. etsmin ) then

cc ETS
              call etsflt(sig_macro)

             elseif(ibryf(ityp,ktyp) .ge. 1 .and.
     &              ichgf(ityp,ktyp) .ne. 0 .and. ityp .ne. 11 .and.
     &  e(ibke+no,ipomp+1) .le.      tsmax*dble(ibryf(ityp,ktyp)) .and.
     &  e(ibke+no,ipomp+1) .ge. emin(ityp)*dble(ibryf(ityp,ktyp))) then

cc ITSART
              call pts_flt2(e(ibke+no,ipomp+1)* 1.d6, sig_macro)
              if(sig_macro .gt. 0.d0 .and. itsxcl .eq. 0) then ! Ogawa 2023/08/18. ITSART exclusive mode without nuclear reactions
                  tottsm = tottsm + sig_macro
                  goto 100
              endif

             end if
            end if
          end if
          tottsm = tottsm + sig_macro
         end if

*-----------------------------------------------------------------------
*
*-----------------------------------------------------------------------
            if( ein .lt. emint .and. totdcy+tottsm .gt. 0.0d0 ) goto 100
            if( ein .lt. emint ) return

*-----------------------------------------------------------------------
*     for no void
*-----------------------------------------------------------------------

      if( mat .gt. 0 ) then

               lem   = nint( dnel_das(kmat0+mat) )
               hydro = denh_das(kmat0+mat)

*-----------------------------------------------------------------------
*        Delta Ray production
*-----------------------------------------------------------------------

         if( jtyp .ne. 0 .and. rtyp .gt. 1.d0 .and.
     &       mndel .gt. 0 .and. delm(icl) .lt. 1.d+10 ) then

               rmass = rtyp
               chag  = dabs(dble( jtyp ))
               edens = edns(kdelt+mat)
               ein   = ein
               emind = delm(icl)

            totdel = delmfp(rmass,chag,edens,ein,emind)

         end if

*-----------------------------------------------------------------------
*        user defined interaction
*-----------------------------------------------------------------------
         if(udm_int_num > 0) then
           call calc_udm_macro_xs(ein)
           if(900000 .le. abs(ktyp) .and. abs(ktyp) .le. 999999)goto 100 ! To avoid user defined particles going to JAM.
         endif

*-----------------------------------------------------------------------
*        non nucleus
*-----------------------------------------------------------------------

         if( ityp .lt. 15 ) then


*-----------------------------------------------------------------------
*           for proton
*-----------------------------------------------------------------------

            if( ityp .eq. 1 ) then

                  dmaxn = das_kmatd(kmatd(mat)+1)
                  dminn = das_kmatd(kmatd(mat)+2)

*-----------------------------------------------------------------------

                  if( hydro .gt. 0.0d0 ) then

                        dsmax = das_kmatd(kmatd(mat)+6)

                     if( ein .gt. dsmax ) then

                        kf1 = kf
                        kf2 = 2212

                        call sigjam(kf1,kf2,ein,sig,sigel,signo)

                        tothyd = tothyd + sig * hydro

                     end if

                  end if

*-----------------------------------------------------------------------

               if( nfcs(ibknfc+no,ipomp+1) .ne. 2 ) then
*                    --- proton pseudo ---

                     do i = 1, lem

                        siggn(ksign+i) = 0.0
                        sigge(ksige+i) = 0.0

                        totttl = totttl + siggm(ksig(mat)+i)

                     end do
*-----------------------------------------------------------------------
cKN 2013/08/20 proton for Coulomb scattering,

                  if( ein .lt. dminn ) then

                        mk = mat
                        rh = denm(mat)

                        eintem = min(ein/5.0,parz(197))  ! use epseudo, T.Sato 2020/09/11
                        eintem = max(eintem,1.0d0)

                        call sig_tot(ityp,sigt,sigaa,eintem,mk)

                      if( iemode.eq.0 ) then ! frtati 2021/12/17

                        tottnd = sigt * rh

                      else
                        tottnd = 0.d0
                        do i = 1, isigc(0)
                          if( ielas.eq.0 ) then
                            tottnd = tottnd + siggcmx(1,i) * rh
                          else
                            tottnd = tottnd + siggcmx(2,i) * rh
                          end if
                        end do
                      end if

                        if( tottnd .gt. totttl ) totttl = tottnd

                  else if( ein .gt. dminn .and. ein .le. dmaxn ) then

                        mk = mat
                        rh = denm(mat)

                        eintem = min(ein/5.0,parz(197))  ! use epseudo, T.Sato 2020/09/11
                        eintem = max(eintem,1.0d0)

                        call sig_tot(ityp,sigt,sigaa,eintem,mk)

                           tottnd = 0.d0

                   if( iemode.eq.0 ) then ! frtati 2021/12/17

                     do i = 1, isigc(0)

                        dsmax = das_kmate(kmate(mk)+(i-1)*5+11)

                        if( ein .le. dsmax ) then

                           tottnd = tottnd + siggc(i) * rh

                        endif

                     enddo

                     do i = 1, lem ! T.Sato 2022/03/14

                        dsmax = das_kmate(kmate(mk)+(i-1)*5+11)

                        if( ein .gt. dsmax ) then

                           tottnd = tottnd + siggm(ksig(mk)+i)

                        endif

                     enddo

                   else
                     do i = 1, isigc(0)
                       if( ielas.eq.0 ) then
                         tottnd = tottnd + siggcmx(1,i) * rh
                       else
                         tottnd = tottnd + siggcmx(2,i) * rh
                       end if
                     end do
                   end if

                        if( tottnd .gt. totttl ) totttl = tottnd

                  end if

*-----------------------------------------------------------------------

               else if( nfcs(ibknfc+no,ipomp+1) .eq. 2 ) then
*                    --- proton real ---

*-----------------------------------------------------------------------
*              mixed material
*-----------------------------------------------------------------------

               if( ein .gt. dminn .and. ein .le. dmaxn ) then

                     lema = 0

                     mk = mat
                     rh = denm(mat)

                     call sig_tot(ityp,sigt,sigaa,ein,mk)

                  do i = 1, isigc(0)

                        dsmax = das_kmate(kmate(mk)+(i-1)*5+11)

                     if( ein .le. dsmax ) then

                        jcoll = 9

                        siggc(i) = siggc(i) * rh

                        lema = lema + 1
                        isigd(lema) = -i
                        siggd(lema) = siggc(i)

                        totttl = totttl + siggc(i)

                     end if

                  end do

*-----------------------------------------------------------------------

                     do i = 1, lem

                        dsmax = das_kmatd(kmatd(mat)+(i-1)*5+11)

                     if( ein .gt. dsmax ) then

                        itz = nint( zz_das(kmat(mat)+i) )
                        ita = nint( a_das(kmat(mat)+i) )

                        call sigrc(ityp,ein,ita,itz,sigt,signe,sigel)

                        if( ielas .eq. 0 ) sigel = 0.0

                        signe = signe*den_das(kmat(mat)+i)
                        sigel = sigel*den_das(kmat(mat)+i)
                        sigt  = signe + sigel

                        siggn(ksign+i) = signe
                        sigge(ksige+i) = sigel

                        lema = lema + 1
                        isigd(lema) = i
                        siggd(lema) = sigt

                        totttl = totttl + sigt

                     end if
                     end do

                        isigd(0) = lema

*-----------------------------------------------------------------------
*              nuclear data
*-----------------------------------------------------------------------

               else if( ein .le. dminn ) then

                     jcoll = 9

                     mk = mat
                     rh = denm(mat)

                  call sig_tot(ityp,sigt,sigaa,ein,mk)

                     do i = 1, isigc(0)

                        if( iemode.eq.0 ) then ! S.H. 2024.12.16
                           siggc(i) = siggc(i) * rh
                           totttl = totttl + siggc(i)
                        else
                           siggcn(i) = siggcn(i) * rh
                           siggce(i) = siggce(i) * rh
                           totttl = totttl + siggcn(i) + siggce(i)
                        end if

                     end do

*-----------------------------------------------------------------------
*              high energy by sigrc
*-----------------------------------------------------------------------

               else if( ein .gt. dmaxn ) then

                     do i = 1, lem

                        itz = nint( zz_das(kmat(mat)+i) )
                        ita = nint( a_das(kmat(mat)+i) )

                        call sigrc(ityp,ein,ita,itz,sigt,signe,sigel)

                        if( ielas .eq. 0 ) sigel = 0.0

                        signe = signe*den_das(kmat(mat)+i)
                        sigel = sigel*den_das(kmat(mat)+i)
                        sigt  = signe + sigel

                        siggn(ksign+i) = signe
                        sigge(ksige+i) = sigel

                        totttl = totttl + sigt

                     end do

               end if

               end if

*-----------------------------------------------------------------------
*           for neutron
*-----------------------------------------------------------------------

            else if( ityp .eq.2 ) then

                  dmaxn = das_kmatd(kmatd(mat)+3)
                  dminn = das_kmatd(kmatd(mat)+4)

*-----------------------------------------------------------------------

                  if( hydro .gt. 0.0d0 ) then

                        dsmax = das_kmatd(kmatd(mat)+7)

                     if( ein .gt. dsmax ) then

                        kf1 = kf
                        kf2 = 2212

                        call sigjam(kf1,kf2,ein,sig,sigel,signo)

                        tothyd = tothyd + sig * hydro

                     end if

                  end if

*-----------------------------------------------------------------------
*              mixed material
*-----------------------------------------------------------------------

               if( ein .gt. dminn .and. ein .le. dmaxn ) then

                     lema = 0

                     mk = mat
                     rh = denm(mat)

                  call xstneu(0,sigt,sigaa,icl,ein,tme,mk)

                  do i = 1, isigc(0)

                        dsmax = das_kmate(kmate(mk)+(i-1)*5+12)

                     if( ein .le. dsmax ) then

                        siggc(i) = siggc(i) * rh

                        lema = lema + 1
                        isigd(lema) = -i

                      if( iemode.eq.0 .or. ein.le.emodem ) then ! frtati 2021/12/17
                        if( iemode .eq. 0 ) then ! frtati 2021/12/17
                           jcoll = 6
                        else if( ein .le. emodem ) then
                           jcoll = 10
                        end if

                        siggd(lema) = siggc(i)
                        totttl = totttl + siggc(i)

                      else
                        siggcn(i) = siggcn(i) * rh
                        siggce(i) = siggce(i) * rh
                        if( ielas.eq.0 ) siggce(i) = 0.d0
                        siggd(lema) = siggcn(i) + siggce(i)
                        totttl = totttl + siggcn(i) + siggce(i)
                      end if

                     end if

                  end do

*-----------------------------------------------------------------------

                     do i = 1, lem

                        dsmax = das_kmatd(kmatd(mat)+(i-1)*5+12)

                     if( ein .gt. dsmax ) then

                        itz = nint( zz_das(kmat(mat)+i) )
                        ita = nint( a_das(kmat(mat)+i) )

                        call sigrc(ityp,ein,ita,itz,sigt,signe,sigel)

                        if( ielas .eq. 0 ) sigel = 0.0

                        signe = signe*den_das(kmat(mat)+i)
                        sigel = sigel*den_das(kmat(mat)+i)
                        sigt  = signe + sigel

                        siggn(ksign+i) = signe
                        sigge(ksige+i) = sigel

                        lema = lema + 1
                        isigd(lema) = i
                        siggd(lema) = sigt

                        totttl = totttl + sigt

                     end if
                     end do

                        isigd(0) = lema

*-----------------------------------------------------------------------
*              nuclear data
*-----------------------------------------------------------------------

               else if( ein .le. dminn ) then

                  if( iemode .eq. 0 .or. ein .gt. emodem ) then
                     jcoll = 6
                  else
                     jcoll = 10
                  end if

                     mk = mat
                     rh = denm(mat)

                  call xstneu(0,sigt,sigaa,icl,ein,tme,mk)

*-----------------------------------------------------------------------

                  do i = 1, isigc(0)

                     siggc(i) = siggc(i) * rh

                   if( iemode.eq.0 .or. ein.le.emodem ) then ! frtati 2021/12/17

                     totttl = totttl + siggc(i)

                   else
                     siggcn(i) = siggcn(i) * rh
                     siggce(i) = siggce(i) * rh
                     if( ielas.eq.0 ) siggce(i) = 0.d0
                     totttl = totttl + siggcn(i) + siggce(i)
                   end if

                  end do

*-----------------------------------------------------------------------
*              high energy by sigrc
*-----------------------------------------------------------------------

               else if( ein .gt. dmaxn ) then

                     do i = 1, lem

                        itz = nint( zz_das(kmat(mat)+i) )
                        ita = nint( a_das(kmat(mat)+i) )

                        call sigrc(ityp,ein,ita,itz,sigt,signe,sigel)

                        if( ielas .eq. 0 ) sigel = 0.0

                        signe = signe*den_das(kmat(mat)+i)
                        sigel = sigel*den_das(kmat(mat)+i)
                        sigt  = signe + sigel

                        siggn(ksign+i) = signe
                        sigge(ksige+i) = sigel

                        totttl = totttl + sigt

                     end do

               end if

*-----------------------------------------------------------------------
*        photon
*-----------------------------------------------------------------------

            else if( ityp .eq. 14 .and. ein .le. dnmax(14) ) then

               mk = mat
               rh = denm(mat)

*-----------------------------------------------------------------------
*        photon for photonuclear reaction
*-----------------------------------------------------------------------
cfrtati 2021/12/17 added library part with [data max]
               if( ipnint .ge. 1 ) then

                 dmaxn = das_kmathd(kmathd(mat)+1)
                 dminn = das_kmathd(kmathd(mat)+2)

                 if( ein.gt.dminn .and. ein.le.dmaxn ) then
                   lema = 0
                   call pnctot(mk,totmpni,ein)
                   do i = 1, isigc(0)
                     dsmax = das_kmate(kmate(mk)+(i-1)*5+13)
                     if( ein.le.dsmax ) then
                       siggc(i) = siggc(i) * rh
                       lema = lema + 1
                       isigd(lema) = -i
                       siggd(lema) = siggc(i)
                       totpni = totpni + siggd(lema)
                     end if
                   end do
                   call xstpni(0,sigt,ein,mk)
                   do i = 1, isigpn0
                     dsmax = das_kmate(kmate(mat)+(i-1)*5+13)
                     if( ein.gt.dsmax ) then
                       siggpn(i) = siggpn(i) * rh
                       lema = lema + 1
                       isigd(lema) = i
                       siggd(lema) = siggpn(i)
                       totpni = totpni + siggd(lema)
                     end if
                   end do
                   isigd(0) = lema
                   totttl = totttl + totpni
                 else if( ein.le.dminn ) then
                   call pnctot(mk,totmpni,ein)
                   do i = 1, isigc(0)
                     siggc(i) = siggc(i) * rh
                     totpni = totpni + siggc(i)
                   end do
                   totttl = totttl + totpni
                 else if( ein.gt.dmaxn ) then
                   call xstpni(0,sigt,ein,mk)
                   totpni = sigt * rh
                   totttl = totttl + totpni
                 end if

               endif

*-----------------------------------------------------------------------
*        photon for muon pair production
*-----------------------------------------------------------------------

               if( igmuppd .ge. 1 ) then

                  call xstgmm(0,sigt,ein,mk)

                  totgmm = sigt * rh
                  totttl = totttl + totgmm

               endif

*-----------------------------------------------------------------------
*        photon for nuclear data
*-----------------------------------------------------------------------

               if( iegsemi .eq. 0 ) then

                  jcoll = 7

                  call xstgam(0,sigt,sigaa,ein,mk)
                  totgam = sigt * rh
                  totttl = totttl + sigt * rh

*---------------------------------------------------------------------
*        photon by EGS5
*---------------------------------------------------------------------

               else if( iegsemi .ne. 0 ) then

                  jcoll = 13 ! Iwase defined

                  call egs5pfpl(tstep,ein,wt(ibkwt+no,ipomp+1),icl,
     &                          totpni,totegs)

                  totttl = totttl + totegs
                  fpl = tstep
                  if( nfcs(ibknfc+no,ipomp+1).ne.0 ) fpl =
     &                                       xfcs(ibkxfc+no,ipomp+1)   ! S.Abe 2015/08/28
                  return
               endif

*-----------------------------------------------------------------------
*        electron for nuclear data
*-----------------------------------------------------------------------

            else if( ( ityp .eq. 12 .or. ityp .eq. 13 ) .and.
     &                 ein .le. dnmax(12) .and. iegsemi .eq. 0 .and.
     &                 tottsm .eq. 0.d0 ) then

                     jcoll = 8

                     mk = mat
                     rh = denm(mat)

               call eletot(pmff,icl,ein,mk,rh)

                     ieletr = 1

                     totttl = totttl + 1.0d0 / pmff

*---------------------------------------------------------------------
*        electron by EGS5
*---------------------------------------------------------------------

            else if( ( ityp .eq. 12 .or. ityp .eq. 13 ) .and.
     &                 ein .le. dnmax(12) .and. iegsemi .ne. 0 .and.
     &                 tottsm .eq. 0.d0 ) then

               jcoll = 14 ! Iwase defined

               mk = mat
               rh = denm(mat)
               tstep = 0d0

               if( ityp.eq.12 ) ich = -1
               if( ityp.eq.13 ) ich =  1

               call egs5efpl(tstep,ein,icl,ich)

               if( tstep.ne.0 ) totttl = totttl + 1.0 / tstep     ! S.Hashimoto, 2016.5.11
               fpl = tstep
               if( nfcs(ibknfc+no,ipomp+1).ne.0 ) fpl =
     &                                          xfcs(ibkxfc+no,ipomp+1)   ! S.Abe 2015/08/28
C              --------------------
               if(udm_int_num > 0) then
C                Calcurate distance to the boundary (dpr)
                 call gomdis(dpr,x(ibkx+no,ipomp+1),y(ibky+no,ipomp+1),
     &                       z(ibkz+no,ipomp+1),
     &                       u(ibku+no,ipomp+1),v(ibkv+no,ipomp+1),
     &                       w(ibkw+no,ipomp+1),
     &                       mark,markp,nmed(ibknmd+no,ipomp+1),
     &                       iblz(ibkblz+no,ipomp+1))
                 if(fpl > dpr) then
                    if(.not. udm_has_step_backed) then
C                   The user-defined interaction occurs with a probability of
C                   udm_sigt (macro XS) * fpl (step length). In steps that cross a boundary,
C                   the process does not go to nreac, causing the reaction probability to be ignored.
C                   To compensate for this, the step is stopped once just before the boundary.
                       fpl=dpr*0.99999d0
                       udm_has_step_backed = .true.
                    endif
                 else
                    udm_has_step_backed = .false.
                 endif

                 if(udm_int_num > 0) udm_fpl = fpl
                 udm_tmp_save=udm_tmp_save+fpl
               endif
C              --------------------
               return

*-----------------------------------------------------------------------
*        for high energy part
*-----------------------------------------------------------------------

            else

*-----------------------------------------------------------------------
*              Hydrogen cross section by JAM for non-nucleon
*-----------------------------------------------------------------------
                     sig   = 0.0
                     sigel = 0.0

                  if( hydro .gt. 0.0d0 ) then

                     kf1 = kf
                     kf2 = 2212

                     call sigjam(kf1,kf2,ein,sig,sigel,signo)

                  end if

                     tothyd = sig * hydro

*-----------------------------------------------------------------------
*                       --- muon ---
*-----------------------------------------------------------------------
               if( ityp .eq. 6 .or. ityp .eq. 7 ) then

                     nmuh = 0

                     do i = 1, lem

                        siggn(ksign+i) = 0.0
                        sigge(ksige+i) = 0.0

                        itz = idnint( zz_das(kmat(mat)+i) )
                        ita = idnint( a_das(kmat(mat)+i) )

                        call sigmu(ein,ita,itz,sigvpi,sigbrm,sigppd)

                        sigvpi = sigvpi*den_das(kmat(mat)+i)
                        sigbrm = sigbrm*den_das(kmat(mat)+i)
                        sigppd = sigppd*den_das(kmat(mat)+i)

                        nmuh = nmuh + 1
                        itzmuh(nmuh) = itz
                        itamuh(nmuh) = ita
                        sigmuv(nmuh) = sigvpi
                        sigmub(nmuh) = sigbrm
                        sigmup(nmuh) = sigppd

                        totttl = totttl + sigvpi + sigbrm + sigppd

                     enddo

                     if( hydro .gt. 0.d0 ) then

                        itz = 1
                        ita = 1

                        call sigmu(ein,ita,itz,sigvpi,sigbrm,sigppd)

                        sigvpi = sigvpi * hydro
                        sigbrm = sigbrm * hydro
                        sigppd = sigppd * hydro

                        nmuh = nmuh + 1
                        itzmuh(nmuh) = itz
                        itamuh(nmuh) = ita
                        sigmuv(nmuh) = sigvpi
                        sigmub(nmuh) = sigbrm
                        sigmup(nmuh) = sigppd

                        totttl = totttl + sigvpi + sigbrm + sigppd

                     endif

*-----------------------------------------------------------------------

               else if( ityp .eq. 12 .or. ityp .eq. 13 ) then
*                       --- electron or positron ---

                     do i = 1, lem

                        siggn(ksign+i) = 0.0
                        sigge(ksige+i) = 0.0

                     end do

               else if( ityp .eq. 14 ) then
*                       --- photon ---

                     do i = 1, lem

                        siggn(ksign+i) = 0.0
                        sigge(ksige+i) = 0.0

                     end do

               else if( ityp .eq. 11 .and. abs(kf) .le. 16 ) then

                if( abs(kf) .eq. 15 ) then
*                       --- tau lepton ---
                else
*                       --- neutrino ---
                call neutrino_Xsec(kf, ein, hydro, mat, lem, sigt, sigh)
                totttl = totttl + sigt
                tothyd = tothyd + sigh
                endif

*-----------------------------------------------------------------------
*                       --- pion ---
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

*           --- pion pseudo ---
                 if( nfcs(ibknfc+no,ipomp+1) .ne. 2 ) then


*              --- for hydrogen ---
                  sigt = 0d0

                  if( hydro .gt. 0.0d0 ) then

                   zt = 1d0
                   at = 1d0
                   call getpseudocs(ein,-1,nint(pichg),
     &                              nint(at),nint(zt),pseudone,pseudoel)
                   sigt = pseudone + pseudoel

                  end if

                  tothyd = sigt * hydro ! consider both elastic and inelastic

*              --- for nucleus ---
                  do i = 1, lem

                   zt = zz_das(kmat(mat)+i)
                   at = a_das(kmat(mat)+i)

                   signe = 0d0
                   sigel = 0d0

                   call getpseudocs(ein,-1,nint(pichg),
     &                              nint(at),nint(zt),pseudone,pseudoel)
                   pseudone = pseudone * den_das(kmat(mat)+i)
                   pseudoel = pseudoel * den_das(kmat(mat)+i)
                   siggn(ksign+i) = pseudone
                   sigge(ksign+i) = pseudoel

                   totttl = totttl + siggn(ksign+i)

                  end do


*           --- pion real ---
                 else

*              --- for hydrogen ---
                  sigt = 0d0

                  if( hydro .gt. 0.0d0 ) then

                   zt = 1d0
                   at = 1d0

                   call pionXS(pichg,ein,at,zt,sigt,signe,sigel,bmax)

                  end if

                  tothyd = sigt * hydro ! consider both elastic and inelastic

*              --- for nucleus ---
                  do i = 1, lem

                   zt = zz_das(kmat(mat)+i)
                   at = a_das(kmat(mat)+i)

                   call pionXS(pichg,ein,at,zt,sigt,signe,sigel,bmax)

                   signe = signe*den_das(kmat(mat)+i)
                   sigel = sigel*den_das(kmat(mat)+i)

                   siggn(ksign+i) = signe
                   sigge(ksige+i) = sigel

                   totttl = totttl + siggn(ksign+i)

                  end do

                 end if

                else ! geometrical model

                 do i = 1, lem

                  siggn(ksign+i) = sigg(ksig(mat)+i)
                  sigge(ksige+i) = 0.0

                  totttl = totttl + siggn(ksign+i)

                 end do

                end if

*-----------------------------------------------------------------------

               else
*                       --- other hadron ---

                  do i = 1, lem

                        siggn(ksign+i) = sigg(ksig(mat)+i)
                        sigge(ksige+i) = 0.0

                        totttl = totttl + siggn(ksign+i)

                  end do

               end if

            end if

*-----------------------------------------------------------------------
*        transport particle is nucleus
*-----------------------------------------------------------------------

         else if( ityp .ge. 15 .and. ityp .le. 19 ) then

C S.H. added two lines to avoid error of [forced collisions] (2023.2.13)
            jta = mtyp
            jtz = jtyp

*-----------------------------------------------------------------------
*           --- nucleus pseudo ---
*-----------------------------------------------------------------------

            if( nfcs(ibknfc+no,ipomp+1) .ne. 2 ) then

           if( ityp.eq.15 .or. ityp.eq.18 ) then
            if(ityp.eq.15) then
             einmevpern = ein/2.0
             dmaxn = das_kmathd(kmathd(mat)+3)
             dminn = das_kmathd(kmathd(mat)+4)
            else if(ityp.eq.18) then
             einmevpern = ein/4.0
             dmaxn = das_kmathd(kmathd(mat)+5)
             dminn = das_kmathd(kmathd(mat)+6)
            endif

*                     --- charged particle pseudo ---
            do i = 1, lem
             siggn(ksign+i) = 0.0
             sigge(ksige+i) = 0.0
            end do
            if( einmevpern .lt. dminn ) then
             mk = mat
             rh = denm(mat)
             eintem = min(ein/5.0,parz(197))  ! use epseudo, T.Sato 2020/09/11
             eintem = max(eintem,1.0d0)
             call sig_tot(ityp,sigt,sigaa,eintem,mk)
             if( iemode.eq.0 ) then
              tottnd = sigt * rh
             else
              tottnd = 0.d0
              do i = 1, isigc(0)
               tottnd = tottnd + siggcmx(1,i) * rh
              end do
             end if
             if( tottnd .gt. totttl ) totttl = tottnd
            else if( einmevpern .ge. dminn .and.  ! T.Sato 2023/04/29 .gt. is changed to .ge.
     &               einmevpern .le. dmaxn ) then
             mk = mat
             rh = denm(mat)
             eintem = min(ein/5.0,parz(197))  ! use epseudo, T.Sato 2020/09/11
             eintem = max(eintem,1.0d0)
             call sig_tot(ityp,sigt,sigaa,eintem,mk)
             tottnd = 0.d0
             do i = 1, isigc(0) ! library part
              if( iemode.eq.0 ) then
               if(ityp.eq.15) then
                dsmax = das_kmate(kmate(mk)+(i-1)*5+14)
               else if (ityp.eq.18) then
                dsmax = das_kmate(kmate(mk)+(i-1)*5+15)
               end if
               if( einmevpern .le. dsmax ) then
                tottnd = tottnd + siggc(i) * rh
               else
                tottnd = tottnd + siggcmx(1,i) * rh
               end if
              else
               tottnd = tottnd + siggcmx(1,i) * rh
              end if
             end do
             if( tottnd .gt. totttl ) totttl = tottnd
            end if
           else
            dminn = 0.d0
            dmaxn = 0.d0
           end if


*-----------------------------------------------------------------------
*              nucleus - proton
*-----------------------------------------------------------------------

                     signe = 0.0
                     sigel = 0.0

               if( hydro .gt. 0.0d0 ) then

                     itypn = 1
                     einn  = rpmass / rtyp * ein
                if( ityp.eq.15 ) then
                  dsmax = das_kmatd(kmatd(mat)+9)
                else if( ityp.eq.18 ) then
                  dsmax = das_kmatd(kmatd(mat)+10)
                else
                  dsmax = 0.d0
                end if
                if( einn.gt.dsmax ) then

                  if( einn .gt. esmin ) then
                    call getpseudocs(einn,jta,jtz,1,1,pseudone,pseudoel)
                    signe = pseudone

                  end if

                end if

               end if

                     tothyd = signe * hydro

*-----------------------------------------------------------------------
*              nucleus - nucleus
*-----------------------------------------------------------------------

                     ap = dble( jta )
                     zp = dble( jtz )

               if( ein/ap.gt.esmin .and. ein/ap.gt.dmaxn ) then


                  do i = 1, lem

                     zt = zz_das(kmat(mat)+i)
                     at = a_das(kmat(mat)+i)

                     signe = 0.d0
                     sigel = 0.d0

                     call getpseudocs(ein/ap,jta,jtz,nint(at),nint(zt),
     &                                pseudone,pseudoel)
                     pseudone = pseudone * den_das(kmat(mat)+i)
                     pseudoel = pseudoel * den_das(kmat(mat)+i)
                     sigt  = pseudone + pseudoel
                     siggn(ksign+i) = pseudone
                     sigge(ksige+i) = pseudoel
                     totttl = totttl + sigt

                  end do

               end if

*-----------------------------------------------------------------------
*           --- nucleus real ---
*-----------------------------------------------------------------------

            else

cfrtati 2023/11/28
           iuselib = 0
           if( ityp.eq.15 .or. ityp.eq.18 ) then
            if(ityp.eq.15) then
             einmevpern = ein/2.0
             dmaxn = das_kmathd(kmathd(mat)+3)
             dminn = das_kmathd(kmathd(mat)+4)
            else if(ityp.eq.18) then
             einmevpern = ein/4.0
             dmaxn = das_kmathd(kmathd(mat)+5)
             dminn = das_kmathd(kmathd(mat)+6)
            endif
            do i = 1, lem
             siggn(ksign+i) = 0.0
             sigge(ksige+i) = 0.0
            end do
            if( einmevpern .lt. dminn ) then
             mk = mat
             rh = denm(mat)
             call sig_tot(ityp,sigt,sigaa,ein,mk)
             do i = 1, isigc(0)
              siggc(i) = siggc(i) * rh
              jcoll = 9
              totttl = totttl + siggc(i)
             end do
             iuselib = 1
            else if( einmevpern .gt. dminn .and.
     &               einmevpern .le. dmaxn ) then
             lema = 0
             mk = mat
             rh = denm(mat)
             call sig_tot(ityp,sigt,sigaa,ein,mk)
             do i = 1, isigc(0)
              if( ityp.eq.15 ) then
               dsmax = das_kmate(kmate(mk)+(i-1)*5+14)
              else if ( ityp.eq.18 ) then
               dsmax = das_kmate(kmate(mk)+(i-1)*5+15)
              end if
              if( einpern.le.dsmax ) then
               lema = lema + 1
               isigd(lema) = -i
               siggc(i) = siggc(i) * rh
               jcoll = 9
               siggd(lema) = siggc(i)
               totttl = totttl + siggc(i)
              end if
             end do
             do i = 1, lem
              if( ityp.eq.15 ) then
               dsmax = das_kmatd(kmatd(mk)+(i-1)*5+14)
              else if ( ityp.eq.18 ) then
               dsmax = das_kmatd(kmatd(mk)+(i-1)*5+15)
              end if
              if( einpern.gt.dsmax ) then
               zt = zz_das(kmat(mat)+i)
               at = a_das(kmat(mat)+i)
               call sighi(ap,zp,ein,at,zt,signe,sigel,bmax)
               if( ielas .le. 2 ) sigel = 0.0
               signe = signe*den_das(kmat(mat)+i)
               sigel = sigel*den_das(kmat(mat)+i)
               sigt  = signe + sigel
               siggn(ksign+i) = signe
               sigge(ksige+i) = sigel
               lema = lema + 1
               isigd(lema) = i
               siggd(lema) = sigt
               totttl = totttl + sigt
              end if
             end do
             isigd(0) = lema
             iuselib = 1
            end if
           end if

*-----------------------------------------------------------------------
*              nucleus - proton
*-----------------------------------------------------------------------

                     signe = 0.0
                     sigel = 0.0

               if( hydro .gt. 0.0d0 ) then

                     itypn = 1
                     einn  = rpmass / rtyp * ein

                     call sigrc(itypn,einn,jta,jtz,sigt,signe,sigel)

               end if

                     tothyd = signe * hydro

*-----------------------------------------------------------------------
*              nucleus - nucleus
*-----------------------------------------------------------------------

                     ap = dble( jta )
                     zp = dble( jtz )

                if( iuselib.eq.0 ) then ! frtati 2023/11/28
                  do i = 1, lem

                     zt = zz_das(kmat(mat)+i)
                     at = a_das(kmat(mat)+i)

                     call sighi(ap,zp,ein,at,zt,signe,sigel,bmax)

                     if( ielas .le. 2 ) sigel = 0.0

                     signe = signe*den_das(kmat(mat)+i)
                     sigel = sigel*den_das(kmat(mat)+i)
                     sigt  = signe + sigel

                     siggn(ksign+i) = signe
                     sigge(ksige+i) = sigel

                     totttl = totttl + signe

                  end do
                end if

            end if

*-----------------------------------------------------------------------

         end if

      end if

*-----------------------------------------------------------------------
*     flight path length by exp. random number
*     distance mesh for electrons
*-----------------------------------------------------------------------

  100 continue

      if( ieletr .eq. 0 ) then

         if(totttl + tothyd + totdcy + totdel + tottsm .le. 0.0d0)return

         rnexp = exprnf(dummy)

         if( nfcs(ibknfc+no,ipomp+1) .eq. 0 ) then

            fpl = rnexp / (totttl + tothyd + totdcy + totdel + tottsm)

         else

            fpl = xfcs(ibkxfc+no,ipomp+1)

         end if

         if( totdcy .gt. 0.d0 .and. ityp .ne. 13 ) then
          tdcylife = rnexp / totdcy / rlit
     &              / dsqrt( ((ein/rtyp)+2.d0) * (ein/rtyp) )
     &              + t(ibkt+no,ipomp+1)
         endif

      else

            fpl = pmff

      end if

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      function delmfp(rmass,chag,edens,ein,emind)
*                                                                      *
*       delta ray flight length                                        *
*       last modified by K.Niita on 2010/07/30                         *
*                                                                      *
*     input:                                                           *
*       rmass : mass of projectile                                     *
*       chag  : charge of projectile                                   *
*       edens : electron density of the material (1/cm^3)              *
*       ein   : energy of projectile                                   *
*       emind : min energy of the delta ray                            *
*                                                                      *
*     output:                                                          *
*       delmfp : delta ray flight length (1/cm)                        *
*                                                                      *
************************************************************************

      implicit real*8(a-h,o-z)

*-----------------------------------------------------------------------

      parameter(e4=2.063d-29) ! e^4 in MeV^2*cm^2
      parameter(rmase=0.511)  ! electron mass
      parameter(pi=3.1415926535897932d0)

*-----------------------------------------------------------------------

               delmfp = 0.0d0

*-----------------------------------------------------------------------
*        Delta Ray production
*-----------------------------------------------------------------------

            emaxd = 4. * rmase / rmass * ein * ( 1. + ein / 2. / rmass )
     &            / ( ( 1. + rmase /rmass )**2
     &                + 2. * rmase * ein / rmass**2 )

            if( emind .lt. emaxd ) then

               beta = sqrt( ein * ( ein + 2. * rmass ) )
     &              / ( ein + rmass )

               smx = 100. * beta / chag**(2./3.)
               zstar = chag * ( 1. - exp( -1.316 * smx
     &               + 0.112 * smx**2 - 0.0650 * smx**3 ) )
               cpara = 2. * pi * edens * e4
     &               * 1.d+24 / rmase / 1.d4

               rate = cpara * ( zstar / beta )**2
     &              * ( 1. / emind - 1. / emaxd )

               delmfp = rate * 1.d+7

            end if

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      function dcyfl(kf,ein,rtyp)
*                                                                      *
*       decay flight length                                            *
*       last modified by K.Niita on 2002/10/15                         *
*                                                                      *
*     output:                                                          *
*                                                                      *
*       dcyfl : decay flight length (cm)                               *
*                                                                      *
************************************************************************

      use udm_Parameter
      use udm_Manager

      implicit real*8(a-h,o-z)

*-----------------------------------------------------------------------

      data rcc / 2.997925d+10 /

      common /dcyd/ rtdcy(50), kfdcy(50), ndcy

*-----------------------------------------------------------------------

               dcyfl = 0.0d0
               tlife = -1.0

*-----------------------------------------------------------------------
*     userdefined particle lifetime
*     "user defined particle" function can modify any particle lifetime.
*     So, this function is written at the beginning of "dcyfl" function.
*-----------------------------------------------------------------------

         if(udm_part_num .gt. 0) then
           udm_kf_for_11=kf
           udm_Kin=ein
           udm_lifetime=0d0
           call user_defined_particle(11) ! get 'udm_lifetime' for 'kf'
           if(udm_lifetime .gt. 0d0) then
             tlife = udm_lifetime
             goto 100
           elseif(900000 .le. abs(kf) .and. abs(kf) .le. 999999) then
             tlife = 1d+30
             goto 100
           endif ! otherwise go back to the normal mode
         endif

*-----------------------------------------------------------------------


         do k = 1, ndcy

            if( kf .eq. kfdcy(k) ) then

               tlife = rtdcy(k)
               goto 100

            end if

         end do

  100    continue

*-----------------------------------------------------------------------
*        decay width
*-----------------------------------------------------------------------

         if( tlife .gt. 0.0 ) then

               dcyfl = 1.0 / tlife / rcc
     &               / dsqrt((( ein / rtyp ) + 2.0 )
     &               * ( ein / rtyp ) )

         end if


*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      function idcyc(kf)
*                                                                      *
*       check decay particle                                           *
*       idcyc = 1 : decay particle                                     *
*       idcyc = 0 : no decay particle                                  *
*       last modified by K.Niita on 2002/10/15                         *
*                                                                      *
************************************************************************

      implicit real*8(a-h,o-z)

*-----------------------------------------------------------------------

      common /dcyd/ rtdcy(50), kfdcy(50), ndcy

*-----------------------------------------------------------------------

         idcyc = 0

*-----------------------------------------------------------------------

         do k = 1, ndcy

            if( kf .eq. kfdcy(k) ) then

               idcyc = 1
               return

            end if

         end do

      return
      end


************************************************************************
*                                                                      *
      function pseudpr(ityp)
*                                                                      *
*       for proton, check pseudo collision                             *
*       last modified by K.Niita on 01/06/2001                         *
*                                                                      *
*     output:                                                          *
*                                                                      *
*       pseudpr  = .true.  : pseudo collision                          *
*                  .false. : real collision                            *
*                                                                      *
************************************************************************

      implicit real*8(a-h,o-z)

*-----------------------------------------------------------------------
      logical pseudpr

      common /xtotal/ totttl,tothyd,totdcy,totdel,tottsm
!$OMP THREADPRIVATE(/xtotal/)

      common /ccxsm/  icxsni, icxspi

*-----------------------------------------------------------------------

       pseudpr = .false.

       if( totttl + tothyd + totdcy + totdel .le. 0.d0 ) return ! do not count tottsm

       tottts = totttl
       tothys = tothyd
       totdes = totdel
       tottss = tottsm

*-----------------------------------------------------------------------
*        proton
*-----------------------------------------------------------------------

         if( ityp .eq. 1 ) then

            call getfltp

*-----------------------------------------------------------------------
*        nucleus
*-----------------------------------------------------------------------

         else if( ityp .ge. 15 ) then

            call getfltn

*-----------------------------------------------------------------------
*        pion
*-----------------------------------------------------------------------

         else if( ityp .ge. 3 .and. ityp .le. 5 .and. icxspi .eq. 1)then  ! PHITS original model

            call getfltpi

*-----------------------------------------------------------------------
*        muon
*-----------------------------------------------------------------------

         else if( ityp .eq. 6 .or. ityp .eq. 7 ) then

            call getfltmu

*-----------------------------------------------------------------------

         else

            return

         end if

*-----------------------------------------------------------------------

         pseud = ( totttl + tothyd + totdcy + totdel )!+ tottsm )
     &         / ( tottts + tothys + totdcy + totdes )!+ tottss )

         if( unirn(dummy) .gt. pseud ) pseudpr = .true.

      return
      end

************************************************************************
      subroutine getrealflt(ityp) ! frtati 2023/12/22
************************************************************************
      if( ityp.eq.1 ) then
        call getfltp
      else if( ityp.ge.15 .and. ityp.le.19 ) then
        call getfltn
      end if

      return
      end

************************************************************************
*                                                                      *
      subroutine getfltp
*                                                                      *
*       for proton, get total flight path length at ec(no)             *
*       last modified by K.Niita on 2010/07/30                         *
*                                                                      *
*     output:                                                          *
*                                                                      *
*       totttl  : total flight length (1/cm)                           *
*                                                                      *
*       siggn(ksign+lem) : nonelastic flight length for (lem,mat)      *
*       sigge(ksige+lem) : elastic flight length for (lem,mat)         *
*                                                                      *
************************************************************************
      use MMBANKMOD !FURUTA
      use MEMBANKMOD !FURUTA
      use ion_track_structure, only : pts_flt2
      use udm_Parameter, only: udm_int_num
      use moddas_material
      implicit real*8(a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'

*-----------------------------------------------------------------------

      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
!$OMP THREADPRIVATE(/icomon/)
      common /tlgeom/ iblz1,iblz2
!$OMP THREADPRIVATE(/tlgeom/)
      common /mathzn/ mathz, mathn, jcoll, kcoll
!$OMP THREADPRIVATE(/mathzn/)
      common /qparm/  ielas,icasc,iqstep,lvlopt,igamma
      common /xtotal/ totttl,tothyd,totdcy,totdel,tottsm
!$OMP THREADPRIVATE(/xtotal/)
      common /kmat1g/ kmat(kvlmax)
      common /regdc/  idrg(kvlmax), idgr(kvmmax)
      common /xinels/ ksige, ksign
      common /ndemax/ dnmax(20)

      common /delesg/ rdels(kvlmax), mndel, ndels(kvlmax)
      common /delreg/ delm(kvlmax), kdelt

      common /celdn/  denr(kvlmax), denm(kvlmax), denc(kvlmax)
      common /regcm/  icmg(kvlmax)


      common /kmat1k/ kmatd(kvlmax), kmate(kvlmax)

      common /  tscmsg   / ktsc(kvlmax), mntsc, ntsc(kvlmax)
      common /  tscreg   / ntscell(kvlmax)
      common / ptsminmax / ptsmin, ptsmax
      common / tsminmax  / tsmax
      common /emode/  emodem, ge1, ge2, iemode

*-----------------------------------------------------------------------
*        only for proton
*-----------------------------------------------------------------------

         if( ityp .ne. 1 ) return

*-----------------------------------------------------------------------
*        initial values and constants
*-----------------------------------------------------------------------

               totttl = 0.0d0
               tothyd = 0.0d0
               totdel = 0.0d0
               tottsm = 0.0d0

               if( mat .le. 0 ) return

*-----------------------------------------------------------------------

               ein = ec(ibkec+no,ipomp+1)
               icl = idgr(iblz1)

               lem   = nint( dnel_das(kmat0+mat) )
               hydro = denh_das(kmat0+mat)

*-----------------------------------------------------------------------
*           user defined interaction
*-----------------------------------------------------------------------
            if(udm_int_num > 0) then
              call calc_udm_macro_xs(ein)
            endif

*-----------------------------------------------------------------------
* Track structure reaction
*-----------------------------------------------------------------------

cc KURBUC
            if( mntsc .gt. 0) then
             sig_macro = 0.d0
             if   (ntscell( idgr(iblz(ibkblz+no,ipomp+1)) ) .eq. 1 .and.
     &          (ityp .eq. 1 )                 .and.
     &          e(ibke+no,ipomp+1) .le. ptsmax .and.
     &          e(ibke+no,ipomp+1) .ge. ptsmin ) then

                call pts_flt(sig_macro)

* ITSART
             elseif(ntscell( idgr(iblz(ibkblz+no,ipomp+1))) .ne. 0 .and.
     &         e(ibke+no,ipomp+1) .le. tsmax*dble(ibryf(ityp,ktyp)).and.
     &         ibryf(ityp,ktyp) .ge. 1    ) then

               call pts_flt2(e(ibke+no,ipomp+1)* 1.d6, sig_macro)

             end if
             tottsm = tottsm + sig_macro
            end if

*-----------------------------------------------------------------------
*           Delta Ray production
*-----------------------------------------------------------------------

            if( jtyp .ne. 0 .and. rtyp .gt. 1.d0 .and.
     &          mndel .gt. 0 .and. delm(icl) .lt. 1.d+10 ) then

                  rmass = rtyp
                  chag  = dabs(dble( jtyp ))
                  edens = edns(kdelt+mat)
                  ein   = ein
                  emind = delm(icl)

               totdel = delmfp(rmass,chag,edens,ein,emind)

            end if

*-----------------------------------------------------------------------
*           for proton
*-----------------------------------------------------------------------

                  dmaxn = das_kmatd(kmatd(mat)+1)
                  dminn = das_kmatd(kmatd(mat)+2)

*-----------------------------------------------------------------------

                  if( hydro .gt. 0.0d0 ) then

                        dsmax = das_kmatd(kmatd(mat)+6)

                     if( ein .gt. dsmax ) then

                        kf1 = ktyp
                        kf2 = 2212

                        call sigjam(kf1,kf2,ein,sig,sigel,signo)

                        tothyd = tothyd + sig * hydro

                     end if

                  end if

*-----------------------------------------------------------------------
*              mixed material
*-----------------------------------------------------------------------

               if( ein .gt. dminn .and. ein .le. dmaxn ) then

                     lema = 0

                     mk = mat
                     rh = denm(mat)

                  call sig_tot(ityp,sigt,sigaa,ein,mk)

                  do i = 1, isigc(0)

                        dsmax = das_kmate(kmate(mk)+(i-1)*5+11)

                     if( ein .le. dsmax ) then

                        siggc(i) = siggc(i) * rh

                        lema = lema + 1
                        isigd(lema) = -i

                      if( iemode.eq.0 ) then ! frtati 2021/12/17

                        jcoll = 9

                        siggd(lema) = siggc(i)
                        totttl = totttl + siggc(i)

                      else
                        siggcn(i) = siggcn(i) * rh
                        siggce(i) = siggce(i) * rh
                        if( ielas.eq.0 ) siggce(i) = 0.d0
                        siggd(lema) = siggcn(i) + siggce(i)
                        totttl = totttl + siggcn(i) + siggce(i)
                      end if

                     end if

                  end do

*-----------------------------------------------------------------------

                     do i = 1, lem

                        dsmax = das_kmatd(kmatd(mat)+(i-1)*5+11)

                     if( ein .gt. dsmax ) then

                        itz = nint( zz_das(kmat(mat)+i) )
                        ita = nint( a_das(kmat(mat)+i) )

                        call sigrc(ityp,ein,ita,itz,sigt,signe,sigel)

                        if( ielas .eq. 0 ) sigel = 0.0

                        signe = signe*den_das(kmat(mat)+i)
                        sigel = sigel*den_das(kmat(mat)+i)
                        sigt  = signe + sigel

                        siggn(ksign+i) = signe
                        sigge(ksige+i) = sigel

                        lema = lema + 1
                        isigd(lema) = i
                        siggd(lema) = sigt

                        totttl = totttl + sigt

                     end if
                     end do

                        isigd(0) = lema

*-----------------------------------------------------------------------
*              nuclear data
*-----------------------------------------------------------------------

               else if( ein .le. dminn ) then

                     mk = mat
                     rh = denm(mat)

                  call sig_tot(ityp,sigt,sigaa,ein,mk)

                  do i = 1, isigc(0)

                     siggc(i) = siggc(i) * rh

                   if( iemode.eq.0 ) then ! frtati 2021/12/17

                     jcoll = 9

                     totttl = totttl + siggc(i)

                   else
                     siggcn(i) = siggcn(i) * rh
                     siggce(i) = siggce(i) * rh
                     totttl = totttl + siggcn(i) + siggce(i)
                   end if

                  end do

*-----------------------------------------------------------------------
*              high energy by sigrc
*-----------------------------------------------------------------------

               else if( ein .gt. dmaxn ) then

                     do i = 1, lem

                        itz = nint( zz_das(kmat(mat)+i) )
                        ita = nint( a_das(kmat(mat)+i) )

                        call sigrc(ityp,ein,ita,itz,sigt,signe,sigel)

                        if( ielas .eq. 0 ) sigel = 0.0

                        signe = signe*den_das(kmat(mat)+i)
                        sigel = sigel*den_das(kmat(mat)+i)
                        sigt  = signe + sigel

                        siggn(ksign+i) = signe
                        sigge(ksige+i) = sigel

                        totttl = totttl + sigt

                     end do

                  jcoll = 0

               end if

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine getfltn
*                                                                      *
*       for nucleus, get total flight path length at ec(no)            *
*       last modified by K.Niita on 2010/07/30                         *
*                                                                      *
*     output:                                                          *
*                                                                      *
*       totttl  : total flight length (1/cm)                           *
*                                                                      *
*       siggn(ksign+lem) : nonelastic flight length for (lem,mat)      *
*       sigge(ksige+lem) : elastic flight length for (lem,mat)         *
*                                                                      *
************************************************************************
      use MMBANKMOD !FURUTA
      use MEMBANKMOD !FURUTA
      use ion_track_structure, only : pts_flt2
      use udm_Parameter, only: udm_int_num
      use moddas_material
      implicit real*8(a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'

      parameter ( rpmass = 938.27, rnmass = 939.58 )

*-----------------------------------------------------------------------

      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
!$OMP THREADPRIVATE(/icomon/)

      common /qparm/  ielas,icasc,iqstep,lvlopt,igamma
      common /xtotal/ totttl,tothyd,totdcy,totdel,tottsm
!$OMP THREADPRIVATE(/xtotal/)

      common /kmat1g/ kmat(kvlmax)
      common /xinels/ ksige, ksign

      common /tlgeom/ iblz1,iblz2
!$OMP THREADPRIVATE(/tlgeom/)
      common /regdc/  idrg(kvlmax), idgr(kvmmax)
      common /delesg/ rdels(kvlmax), mndel, ndels(kvlmax)
      common /delreg/ delm(kvlmax), kdelt

      common /celdn/  denr(kvlmax), denm(kvlmax), denc(kvlmax) ! T.Sato 2021/08/13
      common /ndemax/ dnmax(20)  ! T.Sato 2021/08/13


      common /  tscmsg   / ktsc(kvlmax), mntsc, ntsc(kvlmax)
      common /  tscreg   / ntscell(kvlmax)
      common / ctsminmax / ctsmin, ctsmax
      common / tsminmax  / tsmax
      common /emode/  emodem, ge1, ge2, iemode
      common /kmat1k/ kmatd(kvlmax), kmate(kvlmax)
      common /kmat1ka/ kmathd(kvlmax), kmathe(kvlmax)

*-----------------------------------------------------------------------
*        only for nucleus
*-----------------------------------------------------------------------

         if( ityp .lt. 15 ) return

*-----------------------------------------------------------------------
*        initial values and constants
*-----------------------------------------------------------------------

               totttl = 0.0d0
               tothyd = 0.0d0
               totdel = 0.0d0
               tottsm = 0.0d0

               if( mat .le. 0 ) return

               ein    = ec(ibkec+no,ipomp+1)
               lem    = nint( dnel_das(kmat0+mat) )
               hydro  = denh_das(kmat0+mat)
               icl    = idgr(iblz1)

*-----------------------------------------------------------------------

               jta = mtyp
               jtz = jtyp

*-----------------------------------------------------------------------
*           user defined interaction
*-----------------------------------------------------------------------
            if(udm_int_num > 0) then
              call calc_udm_macro_xs(ein)
            endif

*-----------------------------------------------------------------------
* Track structure reaction
*-----------------------------------------------------------------------

cc KURBUC
            if( mntsc .gt. 0) then
             sig_macro = 0.d0
             if(ntscell( idgr(iblz(ibkblz+no,ipomp+1)) ) .eq. 1 .and.
     &          ktyp .eq. 6000012              .and.
     &          e(ibke+no,ipomp+1) .le. ctsmax .and.
     &          e(ibke+no,ipomp+1) .ge. ctsmin ) then

                call cts_flt(sig_macro)

             elseif(ntscell( idgr(iblz(ibkblz+no,ipomp+1))) .ne. 0 .and.
     &         e(ibke+no,ipomp+1) .le. tsmax*dble(ibryf(ityp,ktyp)).and.
     &         ibryf(ityp,ktyp) .ge. 1    ) then

                call pts_flt2(e(ibke+no,ipomp+1)* 1.d6, sig_macro)

             end if
             tottsm = tottsm + sig_macro
            end if

*-----------------------------------------------------------------------
*           Delta Ray production
*-----------------------------------------------------------------------

            if( jtyp .ne. 0 .and. rtyp .gt. 1.d0 .and.
     &          mndel .gt. 0 .and. delm(icl) .lt. 1.d+10 ) then

                  rmass = rtyp
                  chag  = dabs(dble( jtyp ))
                  edens = edns(kdelt+mat)
                  ein   = ein
                  emind = delm(icl)

               totdel = delmfp(rmass,chag,edens,ein,emind)

            end if


*-----------------------------------------------------------------------
*              nucleus - proton
*-----------------------------------------------------------------------

               dmaxn = dnmax(ityp) ! T.Sato 2021/06/13
               dminn = dnmax(ityp) ! T.Sato 2021/06/13
               if( ityp.eq.15 ) then
                 dmaxn = das_kmathd(kmathd(mat)+3)
                 dminn = das_kmathd(kmathd(mat)+4)
               else if( ityp.eq.18 ) then
                 dmaxn = das_kmathd(kmathd(mat)+5)
                 dminn = das_kmathd(kmathd(mat)+6)
               end if

                     signe = 0.0
                     sigel = 0.0

               if( hydro .gt. 0.0d0 ) then

                     itypn = 1
                     einn  = rpmass / rtyp * ein
                if( ityp.eq.15 ) then
                  dsmax = das_kmatd(kmatd(mat)+9)
                else if( ityp.eq.18 ) then
                  dsmax = das_kmatd(kmatd(mat)+10)
                else
                  dsmax = 0.d0
                end if
                if( einn.gt.dsmax ) then

                     call sigrc(itypn,einn,jta,jtz,sigt,signe,sigel)

                end if

               end if

                     tothyd = signe * hydro

*-----------------------------------------------------------------------
*              nucleus - nucleus
*-----------------------------------------------------------------------

*-----------------------------------------------------------------------
*              mixed material setting for composite use of lib. and model ! frtati 2021/12/17
*-----------------------------------------------------------------------

           einpern = ein
           if(ityp.eq.15) then
            einpern = ein/2.0
           elseif(ityp.eq.18) then
            einpern = ein/4.0
           endif
           ap = dble( jta )
           zp = dble( jtz )

               if( einpern.gt.dminn .and. einpern.le.dmaxn ) then
                 lema = 0
                 mk = mat
                 rh = denm(mat)
                 call sig_tot(ityp,sigt,sigaa,ein,mk)
                 do i = 1, isigc(0)
                   if( ityp.eq.15 ) then
                     dsmax = das_kmate(kmate(mk)+(i-1)*5+14)
                   else if ( ityp.eq.18 ) then
                     dsmax = das_kmate(kmate(mk)+(i-1)*5+15)
                   end if
                   if( einpern.le.dsmax ) then
                     lema = lema + 1
                     isigd(lema) = -i
                     siggc(i) = siggc(i) * rh
                     if( iemode.eq.0 ) then
                       jcoll = 9
                       siggd(lema) = siggc(i)
                       totttl = totttl + siggc(i)
                     else
                       siggcn(i) = siggcn(i) * rh
                       siggd(lema) = siggcn(i)
                       totttl = totttl + siggcn(i)
                     end if
                   end if
                 end do
                 do i = 1, lem
                   if( ityp.eq.15 ) then
                     dsmax = das_kmatd(kmatd(mk)+(i-1)*5+14)
                   else if ( ityp.eq.18 ) then
                     dsmax = das_kmatd(kmatd(mk)+(i-1)*5+15)
                   end if
                   if( einpern.gt.dsmax ) then
                     zt = zz_das(kmat(mat)+i)
                     at = a_das(kmat(mat)+i)
                     call sighi(ap,zp,ein,at,zt,signe,sigel,bmax)
                     if( ielas .le. 2 ) sigel = 0.0
                     signe = signe*den_das(kmat(mat)+i)
                     sigel = sigel*den_das(kmat(mat)+i)
                     sigt  = signe + sigel
                     siggn(ksign+i) = signe
                     sigge(ksige+i) = sigel
                     lema = lema + 1
                     isigd(lema) = i
                     siggd(lema) = sigt
                     totttl = totttl + sigt
                   end if
                 end do

                 isigd(0) = lema

*-----------------------------------------------------------------------
*              nuclear data ! T.Sato 2021/08/14
*-----------------------------------------------------------------------

               else if( einpern .le. dminn ) then

                     mk = mat
                     rh = denm(mat)

                  call sig_tot(ityp,sigt,sigaa,ein,mk)

                  do i = 1, isigc(0)

                     siggc(i) = siggc(i) * rh

                     if( iemode.eq.0 ) then ! frtati 2021/12/17

                       jcoll = 9

                       totttl = totttl + siggc(i)

                     else
                       siggcn(i) = siggcn(i) * rh
                       totttl = totttl + siggcn(i)
                     end if

                  end do

               else if( einpern.gt.dmaxn ) then ! frtati 2021/12/17

                     ap = dble( jta )
                     zp = dble( jtz )

                  do i = 1, lem

                     zt = zz_das(kmat(mat)+i)
                     at = a_das(kmat(mat)+i)

                     call sighi(ap,zp,ein,at,zt,signe,sigel,bmax)

                     if( ielas .le. 2 ) sigel = 0.0

                     signe = signe*den_das(kmat(mat)+i)
                     sigel = sigel*den_das(kmat(mat)+i)
                     sigt  = signe + sigel

                     siggn(ksign+i) = signe
                     sigge(ksige+i) = sigel

                     totttl = totttl + sigt

                  end do

               endif

*-----------------------------------------------------------------------

      return
      end



************************************************************************
*                                                                      *
      subroutine getfltpi
*                                                                      *
*       for pion, get total flight path length at ec(no)               *
*       last modified by S.Hashimoto on 2016/07/21                     *
*                                                                      *
*     output:                                                          *
*                                                                      *
*       totttl  : total flight length (1/cm)                           *
*                                                                      *
*       siggn(ksign+lem) : nonelastic flight length for (lem,mat)      *
*       sigge(ksige+lem) : elastic flight length for (lem,mat)         *
*                                                                      *
************************************************************************
      use MMBANKMOD !FURUTA
      use MEMBANKMOD !FURUTA
      use udm_Parameter, only: udm_int_num
      use moddas_material
      implicit real*8(a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'

*-----------------------------------------------------------------------

      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
!$OMP THREADPRIVATE(/icomon/)
      common /tlgeom/ iblz1,iblz2
!$OMP THREADPRIVATE(/tlgeom/)
      common /mathzn/ mathz, mathn, jcoll, kcoll
!$OMP THREADPRIVATE(/mathzn/)
      common /qparm/  ielas,icasc,iqstep,lvlopt,igamma
      common /xtotal/ totttl,tothyd,totdcy,totdel,tottsm
!$OMP THREADPRIVATE(/xtotal/)
      common /kmat1g/ kmat(kvlmax)
      common /regdc/  idrg(kvlmax), idgr(kvmmax)
      common /xinels/ ksige, ksign
      common /ndemax/ dnmax(20)

      common /delesg/ rdels(kvlmax), mndel, ndels(kvlmax)
      common /delreg/ delm(kvlmax), kdelt

      common /celdn/  denr(kvlmax), denm(kvlmax), denc(kvlmax)
      common /regcm/  icmg(kvlmax)


*-----------------------------------------------------------------------
*        only for pion
*-----------------------------------------------------------------------

      if( ityp .lt. 3 .or. ityp .gt. 5 ) return

*-----------------------------------------------------------------------
*        initial values and constants
*-----------------------------------------------------------------------

      totttl = 0.0d0
      tothyd = 0.0d0
      totdel = 0.0d0
      tottsm = 0.0d0

      if( mat .le. 0 ) return

*-----------------------------------------------------------------------

      ein = ec(ibkec+no,ipomp+1)
      icl = idgr(iblz1)

*-----------------------------------------------------------------------
*           user defined interaction
*-----------------------------------------------------------------------
            if(udm_int_num > 0) then
              call calc_udm_macro_xs(ein)
            endif

*-----------------------------------------------------------------------
*           Delta Ray production
*-----------------------------------------------------------------------

            if( jtyp .ne. 0 .and. rtyp .gt. 1.d0 .and.
     &          mndel .gt. 0 .and. delm(icl) .lt. 1.d+10 ) then

                  rmass = rtyp
                  chag  = dabs(dble( jtyp ))
                  edens = edns(kdelt+mat)
                  ein   = ein
                  emind = delm(icl)

               totdel = delmfp(rmass,chag,edens,ein,emind)

            end if

*-----------------------------------------------------------------------
*     cross section by pionXS
*-----------------------------------------------------------------------

      if ( ityp .eq. 3 ) then
       pichg = 1d0
      else if ( ityp .eq. 5 ) then
       pichg = -1d0
      else
       pichg = 0d0
      end if

*-----------------------------------------------------------------------
*     Hydrogen cross section
*-----------------------------------------------------------------------

      sig   = 0d0
      hydro = denh_das(kmat0+mat)

      if( hydro .gt. 0.0d0 ) then

       zt = 1d0
       at = 1d0

       call pionXS(pichg,ein,at,zt,sig,signe,sigel,bmax)

      end if

      tothyd = sig * hydro     ! consider both elastic and inelastic

*-----------------------------------------------------------------------
*     total, elastic and nonelastic cross sections
*-----------------------------------------------------------------------

      lem  = nint( dnel_das(kmat0+mat) )

      do i = 1, lem

       zt = zz_das(kmat(mat)+i)
       at = a_das(kmat(mat)+i)

       call pionXS(pichg,ein,at,zt,sigt,signe,sigel,bmax)

       signe = signe*den_das(kmat(mat)+i)
       sigel = sigel*den_das(kmat(mat)+i)

       siggn(ksign+i) = signe
       sigge(ksige+i) = sigel

       totttl = totttl + siggn(ksign+i)

      end do

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine getfltmu
*                                                                      *
*       for muon, get total flight path length at ec(no)               *
*       last modified by S.Abe on 2016/07/21                           *
*                                                                      *
*     output:                                                          *
*                                                                      *
*       totttl  : total flight length (1/cm)                           *
*                                                                      *
*       siggn(ksign+lem) : nonelastic flight length for (lem,mat)      *
*       sigge(ksige+lem) : elastic flight length for (lem,mat)         *
*                                                                      *
************************************************************************

      use MMBANKMOD !FURUTA
      use MEMBANKMOD !FURUTA
      use udm_Parameter, only: udm_int_num
      use moddas_material
      implicit real*8(a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'

*-----------------------------------------------------------------------

      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
!$OMP THREADPRIVATE(/icomon/)

      common /xtotal/ totttl,tothyd,totdcy,totdel,tottsm
!$OMP THREADPRIVATE(/xtotal/)

      common /kmat1g/ kmat(kvlmax)
      common /xinels/ ksige, ksign

      common /tlgeom/ iblz1,iblz2
!$OMP THREADPRIVATE(/tlgeom/)
      common /regdc/  idrg(kvlmax), idgr(kvmmax)
      common /delesg/ rdels(kvlmax), mndel, ndels(kvlmax)
      common /delreg/ delm(kvlmax), kdelt


*-----------------------------------------------------------------------
*        only for muon
*-----------------------------------------------------------------------

         if( ityp .ne. 6 .and. ityp .ne. 7 ) return

*-----------------------------------------------------------------------
*        initial values and constants
*-----------------------------------------------------------------------

               totttl = 0.0d0
               tothyd = 0.0d0
               totdel = 0.0d0
               tottsm = 0.0d0

               if( mat .le. 0 ) return

               ein    = ec(ibkec+no,ipomp+1)
               lem   = nint( dnel_das(kmat0+mat) )
               hydro = denh_das(kmat0+mat)
               icl    = idgr(iblz1)

*-----------------------------------------------------------------------
*           user defined interaction
*-----------------------------------------------------------------------
            if(udm_int_num > 0) then
              call calc_udm_macro_xs(ein)
            endif

*-----------------------------------------------------------------------
*           Delta Ray production
*-----------------------------------------------------------------------

            if( jtyp .ne. 0 .and. rtyp .gt. 1.d0 .and.
     &          mndel .gt. 0 .and. delm(icl) .lt. 1.d+10 ) then

                  rmass = rtyp
                  chag  = dabs(dble( jtyp ))
                  edens = edns(kdelt+mat)
                  ein   = ein
                  emind = delm(icl)

               totdel = delmfp(rmass,chag,edens,ein,emind)

            end if

*-----------------------------------------------------------------------
*           Muon reaction
*-----------------------------------------------------------------------

      do i = 1, lem

         siggn(ksign+i) = 0.0
         sigge(ksige+i) = 0.0

         itz = idnint( zz_das(kmat(mat)+i) )
         ita = idnint( a_das(kmat(mat)+i) )

         call sigmu(ein,ita,itz,sigvpi,sigbrm,sigppd)

         sigvpi = sigvpi*den_das(kmat(mat)+i)
         sigbrm = sigbrm*den_das(kmat(mat)+i)
         sigppd = sigppd*den_das(kmat(mat)+i)

         totttl = totttl + sigvpi + sigbrm + sigppd

      enddo

      if( hydro .gt. 0.d0 ) then

         itz = 1
         ita = 1

         call sigmu(ein,ita,itz,sigvpi,sigbrm,sigppd)

         sigvpi = sigvpi * hydro
         sigbrm = sigbrm * hydro
         sigppd = sigppd * hydro

         totttl = totttl + sigvpi + sigbrm + sigppd

      endif

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine delprd(iprj,kprj,rmass,chag,eein,emind,eezo)
*                                                                      *
*                                                                      *
*       delta electron production                                      *
*       modified by K.Niita on 2011/05/17                              *
*                                                                      *
*     input  :                                                         *
*                                                                      *
*       iprj  : typr of projectile                                     *
*       kprj  : kf code of projectile                                  *
*       rmass : mass of projectile                                     *
*       chag  : charge of projectile                                   *
*       eein  : initial energy of projectile (MeV)                     *
*       emind : min energy of the delta ray                            *
*       eezo  : start energy of proj at this step (MeV)                *
*                                                                      *
*     output :                                                         *
*                                                                      *
*       dedxfd : dedx factor scaled by delta energy loss               *
*                                                                      *
*---- in common -------------------------------------------------------*
*                                                                      *
*        nclst  = 2  : total number of out going particles and nuclei  *
*                                                                      *
*        iclust(nclst)                                                 *
*                                                                      *
*                i = 0, nucleus                                        *
*                  = 1, proton                                         *
*                  = 2, neutron                                        *
*                  = 3, pion                                           *
*                  = 4, photon                                         *
*                  = 5, kaon                                           *
*                  = 6, muon                                           *
*                  = 7, others                                         *
*                                                                      *
*        jclust(i,nclst)                                               *
*                                                                      *
*                i = 0, angular momentum                               *
*                  = 1, proton number                                  *
*                  = 2, neutron number                                 *
*                  = 3, ip, see below                                  *
*                  = 4,                                                *
*                  = 5, charge                                         *
*                  = 6, baryon number                                  *
*                  = 7, kf code                                        *
*                  = 8, isomer level (0: Ground, 1,2: 1st, 2nd isomer) *
*                                                                      *
*        qclust(i,nclst)                                               *
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
*        numpat(i) : total number of out going particles or nuclei     *
*                                                                      *
*                i =  0, nuclei                                        *
*                  =  1, proton                                        *
*                  =  2, neutron                                       *
*                  =  3, pi+                                           *
*                  =  4, pi0                                           *
*                  =  5, pi-                                           *
*                  =  6, mu+                                           *
*                  =  7, mu-                                           *
*                  =  8, K+                                            *
*                  =  9, K0                                            *
*                  = 10, K-                                            *
*                  = 11, other particles                               *
*                  = 14, gammma                                        *
*                  = 15, deuteron                                      *
*                  = 16, triton                                        *
*                  = 17, 3He                                           *
*                  = 18, Alpha                                         *
*                                                                      *
************************************************************************
      use MMBANKMOD !FURUTA
      implicit real*8(a-h,o-z)

*-----------------------------------------------------------------------

      parameter ( pi = 3.1415926535898 )
      parameter ( rmase = 0.511)  ! electron mass

*-----------------------------------------------------------------------

      include 'param00.inc'
      include 'param.inc'

*-----------------------------------------------------------------------

      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
!$OMP THREADPRIVATE(/icomon/)

*-----------------------------------------------------------------------

      common /clustf/ nclst, iclust(nnn)
!$OMP THREADPRIVATE(/clustf/)
      common /clustg/ jclust(0:8,nnn), qclust(0:12,nnn)
!$OMP THREADPRIVATE(/clustg/)
      common /clustp/ rumpat(0:20), numpat(0:20)
!$OMP THREADPRIVATE(/clustp/)

*-----------------------------------------------------------------------
      common /dedxfac/ dedxfd
!$OMP THREADPRIVATE(/dedxfac/)
      common /reslet/ irlet

*-----------------------------------------------------------------------
*        maxmum energy of delta ray
*-----------------------------------------------------------------------

               emaxd = 4. * rmase / rmass * eein
     &               * ( 1. + eein / 2. / rmass )
     &               / ( ( 1. + rmase /rmass )**2
     &                 + 2. * rmase * eein / rmass**2 )

              if( emaxd .le. emind ) return

*-----------------------------------------------------------------------
*        determine the delta energy
*-----------------------------------------------------------------------

               r = unirn(dummy)

               edel = 1.d0 / (
     &                1.d0 / emind * ( 1.d0 - r )
     &              + 1.d0 / emaxd * r )

*-----------------------------------------------------------------------
*        edel < eezo - eein ; modify dedxfd
*        edel > eezo - eein ; modify dedxfd and efin
*-----------------------------------------------------------------------

               ddex = eezo - eein

            if( irlet .eq. 1 ) then

               dedxfd = 1.d0
               efin = max( 0.d0, eein - edel )

            else

               if( edel .lt. ddex ) then

                  dedxfd = ( ddex - edel ) / ddex
                  efin = eein

               else

                  dedxfd = 1.0d-40
                  efin = max( 0.d0, eein - ( edel - ddex ) )

               end if

            endif

*-----------------------------------------------------------------------
*        cosine
*-----------------------------------------------------------------------

               pprj = sqrt( efin**2 + 2.0 * rmass * efin )
               pdel = sqrt( edel**2 + 2.0 * rmase * edel )
               etot = efin + rmass - rmase

               cos1 = edel * etot / pprj / pdel
               sin1 = sqrt( 1.0 - cos1**2 )
               phi1 = 2.0 * pi * unirn(dummy)

               pxe1 = pdel * sin1 * cos(phi1)
               pye1 = pdel * sin1 * sin(phi1)
               pze1 = pdel * cos1

*-----------------------------------------------------------------------
*        booking of outgoing two particles ( 1:elctron, 2:projectile )
*-----------------------------------------------------------------------

               nclst = 2

               numpat(12)   = numpat(12) + 1
               rumpat(12)   = rumpat(12) + 1

               iclust(1)    = ipatf(12,11)

               jclust(0,1)  = 0
               jclust(1,1)  = 0
               jclust(2,1)  = 0
               jclust(3,1)  = 12
               jclust(4,1)  = 0
               jclust(5,1)  = ichgf(12,11)
               jclust(6,1)  = ibryf(12,11)
               jclust(7,1)  = 11
               jclust(8,1)  = 0

               qclust(0,1)  = 0.0
               qclust(1,1)  = pxe1 / 1000.
               qclust(2,1)  = pye1 / 1000.
               qclust(3,1)  = pze1 / 1000.
               qclust(4,1)  = ( edel + rmase ) / 1000.
               qclust(5,1)  = rmase / 1000.
               qclust(6,1)  = 0.0
               qclust(7,1)  = edel
               qclust(8,1)  = 1.0
               qclust(9,1)  = 0.0
               qclust(10,1) = 0.0d0
               qclust(11,1) = 0.0d0
               qclust(12,1) = 0.0d0

               numpat(iprj) = numpat(iprj) + 1
               rumpat(iprj) = rumpat(iprj) + 1

               iclust(2)    = ipatf(iprj,kprj)

               jclust(0,2)  = 0
               jclust(1,2)  = ichgf(iprj,kprj)
               jclust(2,2)  = ibryf(iprj,kprj) - ichgf(iprj,kprj)
               jclust(3,2)  = iprj
               jclust(4,2)  = 0
               jclust(5,2)  = ichgf(iprj,kprj)
               jclust(6,2)  = ibryf(iprj,kprj)
               jclust(7,2)  = kprj
               jclust(8,2)  = 0

               qclust(0,2)  = 0.0
               qclust(1,2)  = 0.0
               qclust(2,2)  = 0.0
               qclust(3,2)  = pprj / 1000.
               qclust(4,2)  = ( efin + rmass ) / 1000.
               qclust(5,2)  = rmass / 1000.
               qclust(6,2)  = 0.0
               qclust(7,2)  = efin
               qclust(8,2)  = 1.0
               qclust(9,2)  = 0.0
               qclust(10,2) = 0.0d0
               qclust(11,2) = 0.0d0
               qclust(12,2) = 0.0d0

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine sigrc(incp,emev,ia,iz,sigt,sigr,sigs)
*                                                                      *
*        calculates total, nonelastic and elastic cross-sections       *
*                                                                      *
*     input:                                                           *
*        incp   : =1, proton, =2, neutron                              *
*        emev   : incident nucleon energy (MeV)                        *
*        ia     : mass number                                          *
*        iz     : charge number                                        *
*                                                                      *
*     output:                                                          *
*       sigt    : total cross-section (b)                              *
*       sigr    : nonelastic cross-section (b)                         *
*       sigs    : sigt-sigr=elastic scattering cross-section (b)       *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      include 'param.inc'
      include 'err.inc'

      common /ccxsm/  icxsni, icxspi
      common /ckurotama/  dsck

      common /scincom/iscinful  ! use SCINFUL or not
      common /scinswt/iswitch
!$OMP THREADPRIVATE(/scinswt/)

      dimension       par(4),c(4,3)

      data (c(1,k),k=1,3)/1.353,0.4993,1.076/
      data (c(2,k),k=1,3)/0.6653,0.9900,1.040/
      data (c(3,k),k=1,3)/2.456,2.060,0.9314/
      data (c(4,k),k=1,3)/0.1167,1.623,0.9736/

*-----------------------------------------------------------------------
C for [frag data] (S.Hashimoto, 2014.12.08)
      common /fgdata/ ifrgd, ifgm, ifgdf(kvlmax,5), frgfl(kvlmax)
      character frgfl*200
C S.Hashimoto added a new flag for [Frag Data]. (2015.8.25)
      common /initfgdata/ init_frag ! S.H. corrected. (2016.9.25)
      integer   init_frag
      data init_frag/1/

*-----------------------------------------------------------------------

      if ( init_frag .eq. 0 .and. ifrgd .gt. 0 ) then

            ap = 1d0
         if (incp .eq. 1) then
            zp = 1d0
            ktypfrg = 2212
         else
            zp = 0d0
            ktypfrg = 2112
         end if
            at = dble(ia)
            zt = dble(iz)
            call fragdataXS(ktypfrg,ap,zp,emev,at,zt,signe,sigel,bmax)

            sigr = signe
            sigs = sigel
            sigt = sigr + sigs
            if ( sigr .gt. 0d0 ) return

      end if

! SCINFUL mode, T.Sato 2020/02/16
! revised  by d.satoh (2020.08.25)
         if(iscinful.ne.0) then
           if(incp.eq.2.and.iz.eq.6.and.ia.eq.12) then ! neutron on 12C

             iswitch = 0

             elow = 8.0d+1             ! (MeV)
             eupp = 1.1d+2             ! (MeV)

             a = 0.d0
             b = 0.d0

             if( emev .le. elow ) then
               prb = 0.0d0
             else if( emev .gt. elow .and. emev .le. eupp ) then
!            Linear function-------------------------------------------
             a = 1.0d+0/(eupp-elow)
             b = 1.0d+0 - a*eupp
             prb = a*emev + b
!            Linear function-------------------------------------------

             else
               prb = 1.0d0
             end if

             rn = unirn(dummy)

*            -----------------------------------------------------------
             if( rn .ge. prb ) then     ! SCINFUL mode
               iswitch = 1
               sigt=totalx(emev)
               sigs=0.0d0
               sigr=sigt
               if(sigr.eq.0.0) then
                 write(ErrCha,'(''for SCINFUL mode, dmax for carbon'',
     &           '' should be greater than 0.10 MeV'')')
                 ErrID = 'L:3026/R:sigrc/F:getflt.f' !E00_000_000
                 call ErrWrite(ErrID,ErrCha)
                 call parastop(122)
               endif
               return
*            -----------------------------------------------------------
             else if( rn .lt. prb .and. emev .le. 3.0d+3 ) then
               iswitch = 0
               sigt=carbonsigma(1,emev)
               sigs=carbonsigma(2,emev)
               sigr=carbonsigma(3,emev)
               if(sigr.eq.0.0) then
                 write(ErrCha,'(''for SCINFUL mode, dmax for carbon'',
     &           '' should be greater than 0.10 MeV'')')
                 ErrID = 'L:3040/R:sigrc/F:getflt.f' !E00_000_000
                 call ErrWrite(ErrID,ErrCha)
                 call parastop(122)
               endif
               return
*            -----------------------------------------------------------
             else
               iswitch = 0
*            -----------------------------------------------------------
             endif

           endif
         endif
*-----------------------------------------------------------------------
C S.H. to avoid an error occured when 1H target (2022.7.1)

         if ( iz.eq.1 .and. ia.eq.1 ) then
            sigt = 0d0
            sigr = 0d0
            sigs = 0d0
            return
         end if

*-----------------------------------------------------------------------

         rz = 0.14

         a = dble(ia)
         z = dble(iz)

         a13 = a**0.3333333
         a23 = a**0.6666667
         alg = log(a)

         ecst = 82.0
         acst = 238.0

*-----------------------------------------------------------------------

         c1 =  0.0825
         c2 = -0.0057
         c3 =  0.14
         c4 = -0.2
         c5 =  2.72
         c6 =  1.62
         c7 = -5.3

*-----------------------------------------------------------------------

         ap = 1.0
         pi = 3.1415927

*-----------------------------------------------------------------------
*     modified by Niita
*-----------------------------------------------------------------------

            cmev = 0.0575 * a + 12.31

            facp = min( 1.0d0, 0.684 + 1.327e-3 * a )
            facc = 1.0 - facp

*-----------------------------------------------------------------------

            g1 = 1.0 - 0.62 * exp( -cmev / 200. )
     &                      * sin( 10.9 / cmev**0.28 )
            g1 = g1 * 1.1
            g2 = 1.0 + 0.016 * sin( 5.3 - 2.63 * log(a) )

            sigpa = 0.045 * a**0.7 * g1 * g2

*-----------------------------------------------------------------------

            eroot = sqrt(cmev)
            rad   = rz * a**0.3333
            wvl   = 0.1 * 1.22 * ( a + ap ) / a * sqrt(14.1) / eroot

            sigcr = pi * ( rad + wvl )**2

*-----------------------------------------------------------------------

            sigpc = facp * sigcr + facc * sigpa

            ftpa = sigpc / sigpa
            ftcr = sigpc / sigcr

            enpa = ftpa * 1.1 - 1.0

*-----------------------------------------------------------------------

         do k=1,4

            par(k) = log(c(k,1))+log(c(k,2))*alg+log(c(k,3))*alg**2
            par(k) = exp(par(k))

         end do

            epk   = par(3) * a13
            fexp1 = par(2) * log(epk/cmev)

            sigtp  = sigpc * ( 1.0 + par(4) )
     &             + par(1) * a13 * exp( -fexp1**2 )

            esub = ecst * a13 / acst**0.3333
            epk2 = epk - esub

         if( epk2 .gt. 0.0 ) then

            fexp2 = par(2) * log( epk2 / cmev )

            sigtp = sigtp + par(1) * a13 * exp( -fexp2**2 )

         end if

*-----------------------------------------------------------------------

            if( eroot .lt. 3.0 ) as = 1.0 - c1 * ( 3.0 - eroot )**2
            if( eroot .ge. 3.0 ) as = 1.0

            p = c2 * eroot + c3
            q = c4 * eroot + c5

            if( eroot .ge. 4.0 ) r = 1.2
            if( eroot .lt. 4.0 ) r = c6 * eroot + c7

            sigtc = 2.0 * sigpc * ( as - p * cos( q * a**0.3333 - r ) )

*-----------------------------------------------------------------------

            facp = min( 1.0d0, 0.578 + 1.77e-3 * a )
            facc = 1.0 - facp

            sigtpc = facp * sigtc + facc * sigtp

            fttp = sigtpc / sigtp
            fttc = sigtpc / sigtc

*-----------------------------------------------------------------------
*     calculation of sigr as
*     J.Letaw etal., Astrophys. J. Supp. series 51,271-276(1983).
*     sigrb = asymptotic high energy reaction cross-section
*-----------------------------------------------------------------------

         sigrb = 0.045 * a**0.7

*-----------------------------------------------------------------------
*     energy dependent factor
*     modified by T.Sato, 2010/8/3
*-----------------------------------------------------------------------

         if (icxsni .eq. 2) then

            f1 = 1.0 - 0.75 * exp( -emev / 200. )
     &                   * sin( 10.9 / emev**0.28 )
     &      +0.09*exp(-((log10(emev)-3.1)**2)/0.3**2)
     &      +0.12*exp(-((log10(emev)-3.6)**2)/0.5**2)

         else

            f1 = 1.0 - 0.62 * exp( -emev / 200. )
     &                   * sin( 10.9 / emev**0.28 )

         endif
*-----------------------------------------------------------------------
*     Pearlstein added low energy enhancement factor of 10%
*-----------------------------------------------------------------------

         f1 = f1 * ( 1.0 + enpa
     &      * exp( -min( 50.d0, ( emev - cmev ) / 10. )) )



*-----------------------------------------------------------------------
*     mass dependent factor
*-----------------------------------------------------------------------

         f2 = 1.0 + 0.016 * sin( 5.3 - 2.63 * log(a) )

*-----------------------------------------------------------------------

         sigr = sigrb * f1 * f2

*-----------------------------------------------------------------------
*     sigt according to fit by S. Pearlstein, June 86.
*-----------------------------------------------------------------------

      if( emev .ge. cmev ) then

         do 210 k=1,4

            par(k) = log(c(k,1))+log(c(k,2))*alg+log(c(k,3))*alg**2
            par(k) = exp(par(k))

  210    continue

            epk   = par(3) * a13
            fexp1 = par(2) * log(epk/emev)

            sigt  = sigr * ( 1.0 + par(4) )
     &            + par(1) * a13 * exp( -fexp1**2 )

            ecst = 82.0
            acst = 238.
            esub = ecst * a13 / acst**0.3333
            epk2 = epk - esub

         if( epk2 .gt. 0.0 ) then

            fexp2 = par(2) * log( epk2 / emev )

            sigt = sigt + par(1) * a13 * exp( -fexp2**2 )

         end if

            factp = 1.0 + ( fttp - 1.0 )
     &            * exp( - min( 50.d0,( emev - cmev ) / 10.0 ))

            sigt = sigt * factp

*-----------------------------------------------------------------------
*     calculates xsects according to Ramsauer effect
*     ref: Angeli and Csikai, NP A170,577-583(1971)
*     sigt/signe=as-p*cos(q*a**0.3333-r)
*-----------------------------------------------------------------------

      else

            eroot = sqrt(emev)

            if( eroot .lt. 3.0 ) as = 1.0 - c1 * ( 3.0 - eroot )**2
            if( eroot .ge. 3.0 ) as = 1.0

            p = c2 * eroot + c3
            q = c4 * eroot + c5

            if( eroot .ge. 4.0 ) r = 1.2
            if( eroot .lt. 4.0 ) r = c6 * eroot + c7

            rad = rz * a**0.3333
            wvl = 0.1 * 1.22 * ( a + ap ) / a * sqrt(14.1) / eroot

         facr = 1.0 + ( ftcr - 1.0 )
     &        * exp( - min( 50.d0,( cmev - emev ) / 10.0 ))

            sigr = pi * ( rad + wvl )**2 * facr
            sigt = 2.0 * sigr * ( as - p * cos( q * a**0.3333 - r ) )

            factc = 1.0 + ( fttc - 1.0 )
     &            * exp( - min( 50.d0, ( cmev - emev ) / 10.0 ))

            sigt = sigt * factc

      end if

*-----------------------------------------------------------------------
*     coulomb factor
*     ec, coulomb barrier half-height, similar to s. pearlstein,
*     j. nuc. energy 23,87(1975) using optical model inverse p cs's.
*-----------------------------------------------------------------------

      if( incp .eq. 1 .and. emev .lt. 200.0 ) then

            ec = 1.44 * z / ( 8.2 + 0.68 * a13 )

            w1    = 3.816 + 0.1974 * z
            fexc1 = exp( max( -50.d0, ( ec - emev ) / w1 ))
            qfac1 = 1.0 / ( 1.0 + fexc1 )

            w2    = 0.07246 * z + 6.058
            if( z .lt. 10.0 ) w2 = 12.0
            qfac2 =  1.0 - exp( - min( 50.d0, ( emev / w2 )**2 ))

            w3    = 2.0
            fexc3 = exp( max( -50.d0, ( ec - emev ) / w3 ))
            qfac3 = 1.0 / ( 1.0 + fexc3 )

            fcoul = qfac1 * qfac2 * qfac3

            sigt = sigt * fcoul
            sigr = sigr * fcoul

      end if

*-----------------------------------------------------------------------
*     Use KUROTAMA model for reaction cross section
*-----------------------------------------------------------------------

      if ( icxsni .eq. 1 .and. (( incp .eq. 1 )
     &     .or. ( incp .eq. 2 .and. emev .ge. 20.d0 ))) then

           if (incp .eq. 1) then
               zp = 1.0
            else
               zp = 0.0
            end if
            at = dble(ia)
            zt = dble(iz)
            call kurotama0(ap,zp,emev,at,zt,dsck,sigr,bmax)
            sigr = sigr / 100.0d0 !(fm^2 -> b)

            if ( sigt .lt. sigr ) sigt = sigr

      end if

*-----------------------------------------------------------------------
*   Increase the ratio of elastic cross section
      if (icxsni .eq. 2) then
       sigr=sigr*1.03               ! T.Sato 2015/2/11
       if (sigr.gt.sigt) sigr=sigt  ! T.Sato 2015/2/11
      endif

*-----------------------------------------------------------------------

            sigs = sigt - sigr

*-----------------------------------------------------------------------

       return

            sigr = sigr * 1.5
            sigt = sigr + sigs

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine seldsd(incp,ex,ata,atz,csthcm,ifirst)
*                                                                      *
*        determine csthcm according to the elastic angular             *
*        distribution by K. Niita's systematics                        *
*                                                                      *
*     input:                                                           *
*        incp     : 1->proton, 2->neutron                              *
*        ex       : kinetic energy of incident nucleon in lab (MeV)    *
*        ata      : target mass number                                 *
*        atz      : target charge                                      *
*        ifirst   : frag of multiple calculations for the same system  *
*                                                                      *
*     output:                                                          *
*        csthcm   : cosine of scattering angle in cm                   *
*                                                                      *
************************************************************************

      implicit real*8(a-h,o-z)

*-----------------------------------------------------------------------

      parameter ( pi = 3.1415926535897932 )

*-----------------------------------------------------------------------

      common /paraj/ mstz(300), parz(300)

      dimension       dsiga(1000), dsigp(1000)

      save            dsiga, dsigp
!$OMP THREADPRIVATE(dsiga,dsigp)
*-----------------------------------------------------------------------
*        mesh for angular distribution
*-----------------------------------------------------------------------

              ianm  = mstz(24)
              andif = 180.0 / ianm

*-----------------------------------------------------------------------
*        ifirst = 1 ; normal,  ne= 1 ; repeatition of the same system
*-----------------------------------------------------------------------

         if( ifirst .ne. 1 ) goto 200

*-----------------------------------------------------------------------

               ia = nint(ata)
               iz = nint(atz)

            call sigrc(incp,ex,ia,iz,sigt,sign,sige)

               sek = 0.0

            do ian = 1, ianm

               aplow = dble(ian-1) * andif
               aphig = dble(ian-1) * andif + andif
               apoin = ( aplow + aphig ) /2.0
               aprad = apoin * pi / 180.0
               apsin = sin( aprad )
               adrad = ( aphig - aplow ) * pi / 180.0

               sigs  = sign
               icct  = 2
               icm   = 0

               call dsdarc(icct,icm,incp,
     &                     ex,ata,apoin,sigs,dsigs,escat,etarg,angle)

               dsigp(ian) = aplow
               dsiga(ian) = dsigs * 2. * pi * apsin * adrad

               sek = sek + dsiga(ian)

            end do

               fnorm = sek

               sek = 0.0

           if ( fnorm .ne. 0d0 ) then ! S.H. to avoid zero division error (2021.1.28)
            do ian = 1, ianm

               sek = sek + dsiga(ian) / fnorm

               dsiga(ian) = sek

            end do
           end if

*-----------------------------------------------------------------------
*        random number 0 < ramx < 1
*-----------------------------------------------------------------------

  200 continue

               ram1 = unirn(dummy)
               ram2 = unirn(dummy)

*-----------------------------------------------------------------------

            do ian = 1, ianm

               if( dsiga(ian) .gt. ram1 ) goto 100

            end do

  100       continue

               angcm = ( dsigp(ian) + ram2 * andif ) * pi / 180.0

               csthcm = cos( angcm )

*-----------------------------------------------------------------------
*           high energy > 1000 GeV  => forward only
*-----------------------------------------------------------------------

                  ee1 =    1000.0
                  ee2 = 1000000.0

                  ec1 = log(ee1)
                  ec2 = log(ee2)

               if( ex .gt. ee1 .and. ex .le. ee2 ) then

                  aa  = 1.0 / ( ec1 - ec2 )
                  bb  = - ec2 * aa

                  csthcm = 1.0 - ( 1.0 - csthcm )
     &                   * ( aa * log(ex) + bb )**8

               else if( ex .gt. ee2 ) then

                  csthcm = 1.0

               end if

*-----------------------------------------------------------------------

      return
                  csthcm = 2. * unirn(dummy) - 1.0

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine dsdarc(icc,icm,incp,
     &                  pnrg,atarg,angle,sigs,dsigs,escat,etarg,angcm)
*                                                                      *
*        compute elastic scattering cross-section using                *
*        Fraunhoffer diffraction (b) for scattering formulae,          *
*        see S. Pearlstein, Nuc. Sci. Eng., 49, p.162-171(1972)        *
*                                                                      *
*        Hutcheon PRL 47,315(81) data 200-400 mev p+pb data            *
*                                                                      *
*     input:                                                           *
*        icc      : 0 -> Pearlstein                                    *
*                 : 1 -> Pearlstein + Hutcheon                         *
*                 : 2 -> Niita                                         *
*        icm      : 0 -> cm, 1 -> lab  input and output                *
*        incp     : 1->proton, 2->neutron                              *
*        pnrg     : kinetic energy of incident nucleon in lab (MeV)    *
*        atarg    : target mass number                                 *
*        angle    : scattering angle (deg) cm or lab                   *
*        sigs     : non-elastic scattering cross-section (b)           *
*                                                                      *
*     output:                                                          *
*        dsigs    : differential elastic cross-section (b/sr)          *
*        escat    : lab kinetic energy of scattered nucleon (MeV)      *
*        etarg    : =pnrg-escat, lab k.e. of scattered target (MeV)    *
*        angcm    : cm scattering angle (deg)                          *
*                                                                      *
************************************************************************

      implicit real*8(a-h,o-z)

*-----------------------------------------------------------------------
*     constants
*-----------------------------------------------------------------------

            rz = 0.14
            pi = 3.1415927

*-----------------------------------------------------------------------
*     transformation lab to cm
*-----------------------------------------------------------------------

            call syscrc(icm,incp,pnrg,atarg,angle,vk,angzlm,escat,conv)

            angcm = angzlm
            etarg = pnrg - escat

            a  = atarg
            ia = nint( a )

*-----------------------------------------------------------------------
*        original elastic cross section
*-----------------------------------------------------------------------

         if( icc .ne. 2 ) then

            r = ( rz * a**0.3333
     &          + 0.122 * ( a + 1.0 ) / a ) * 1.d-12

            del = 0.05 * r

            area = pi * r * r * 1.d+24
            vkr2 = 2.0 * vk * r

            call bessel(0,vkr2,besz)
            call bessel(1,vkr2,bes1)

            sig = area * ( 1.0 - besz**2 - bes1**2 )

*-----------------------------------------------------------------------
*        new parametrization by K. Niita
*-----------------------------------------------------------------------

         else

            r = ( sqrt( sigs * 100. ) / pi + 0.0 ) * 1.d-13

            r = r * ( -43.78 / ( a + 136.6 ) + 2.23 )

            r = r * ( 1.0
     &            + ( ( 0.8 + 0.2
     &            *  exp( - min( 50.d0,( a / 100.0 )**4 )) ) - 1.0 )
     &            *  exp( - min( 50.d0,( pnrg / 40.0 )**2 )) )

            del = 0.45 * 1.d-13

            rfac = 0.66 * exp( - min( 50.d0,( 541.4 / pnrg )**2) )
            dfac = 1.80 * exp( - min( 50.d0,( 61.43 / a )**2) )

            r = r * ( 1.0 + rfac )

            del = del * ( 1.0 + dfac * rfac * 1.35 )

         end if

*-----------------------------------------------------------------------
*     differential cross section
*-----------------------------------------------------------------------

            y0 = ( vk * r * r )**2 * 1.d+24

            u = cos( angzlm / 180.0 * pi )
            s = sqrt( ( 1.0 - u ) / 2.0 )

            if( u .eq. 1.0 ) s = 0.0

            xav = 2.0 * vk * r * s

         if( xav .le. 1.d-2 ) then

            xj1x = 0.5
            y1   = xj1x * xj1x
            y2   = 1.0

         else

            xpl = 2.0 * vk * ( r + del ) * s
            xmi = 2.0 * vk * ( r - del ) * s

            call bessel(1,xpl,xj1pl)
            call bessel(1,xmi,xj1mi)

            if( icc .ne. 2 ) then

               fp1 = 1.0
               fp2 = 1.0

            else

               fp1 = 1.4
               fp2 = 2.0 - fp1

            end if

            y1 = 0.5 * (
     &           fp1 * ( xj1pl / xpl )**2 +
     &           fp2 * ( xj1mi / xmi )**2 )

*-----------------------------------------------------------------------
*           following fits Hutcheon PRL 47,315(81)
*           data 200-400 mev p+pb data
*           icc = 0 -> Pearlstein
*                 1 -> Pearlstein + Hutcheon
*                 2 -> Niita
*-----------------------------------------------------------------------

            if( icc .eq. 0 ) then

               fac0 = 1.0

               y2 = fac0

            else if( icc .eq. 1 ) then

               cf1  = 0.2 * ( 208. / a )**0.3333
               fac1 = exp( - min( 50.d0,
     &                            cf1 * xav + cf1 / 40. * xav**2 ) )

               y2 = fac1

            else if( icc .eq. 2 ) then

               cf2 = 0.1 * ( 100. / a )**0.3333
               facc = exp( - min( 50.d0,
     &                            cf2 * xav + cf2 / 40. * xav**2 ) )

               emas = ( -14512.0 / ( a + 103.5 ) + 146.7 )
               fac2 = 1.0 + ( facc - 1.0 )
     &              * ( 1.0 - exp( - min( 50.d0,
     &                                  ( pnrg / emas )**4 ) ) )

               cf3 = 2.5 * ( 80. / a )**0.3333
               angc = 160.0 * 50.0 / ( pnrg + 50.0 ) + 20.0

               angd = - 15.0 / 190.0 * a + 20.79
               ange =   15.0 / 190.0 * a + 49.2

              angc = ( 180.0 - angd ) * ange / ( pnrg + ange ) + angd

               faca = exp( - min( 50.d0, cf3 * ( angzlm / angc )**2 ) )

               fac3 = 1.0 + ( faca - 1.0 )
     &              * ( 1.0 - exp( - min( 50.d0,
     &                                  ( pnrg / 60.0 )**2 ) ) )

               y2 = fac2 * fac3

            end if


         end if

*-----------------------------------------------------------------------

            dsigs = conv * y0 * y1 * y2

            if( abs(dsigs) .le. 1.d-12 ) dsigs = 0.0


*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine syscrc(icm,incp,e,a,angle,vk,angzlm,escat,conv)
*                                                                      *
*        to convert from lab parameters                                *
*        to zero linear momentum (zlm) frame                           *
*        for relativistic formulae,                                    *
*        see J.L. Fowler, J.E. Brolley, Jr.,                           *
*        Rev. Mod. Phys. 28,103-134(1956)                              *
*                                                                      *
*     input:                                                           *
*        icm      : 0 -> cm, 1 -> lab  input and output                *
*        incp     : 1->proton, 2->neutron                              *
*        e        : incident energy in lab (MeV)                       *
*        angle    : scattering angle (deg) cm or lab                   *
*                                                                      *
*     output:                                                          *
*        vk       : wave number (1/cm) in lab                          *
*        angzlm   : scattering angle in cm system (deg)                *
*        escat    : elastic scattering energy in lab (MeV)             *
*        conv     : factor to convert cross-section from c.m.          *
*                   to lab system                                      *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      pi = 3.1415927

*-----------------------------------------------------------------------

            if( incp .eq. 1 ) then

               rmass = 938.3

            else if( incp .eq. 2 ) then

               rmass = 939.58

            end if

            ang = angle * pi / 180.
            rat = e / rmass
            gz  = 1. + rat
            s   = 1. + 2. * a * gz + a**2

*-----------------------------------------------------------------------
*     if nonrelativistic, wave number in lab, vk
*-----------------------------------------------------------------------

            vk = 2.197e+12 * dsqrt(e)

*-----------------------------------------------------------------------
*     corrected for relativistic and zlm motion
*-----------------------------------------------------------------------

            add = 1.
            vk  = vk * a / ( a + add )
     &          * dsqrt( ( 1. + rat / 2. ) / ( 1. + rat )**2 )

*-----------------------------------------------------------------------

         if( icm .eq. 0 ) then

            angzlm = angle

            cosa = cos( angzlm * pi / 180. )
            conv = 1.0

*-----------------------------------------------------------------------

         else

               gcg=(a+gz)/dsqrt(s)
               g1=(1.+a*gz)/dsqrt(s)
               a1=(1.+a*gz)/(a+gz)/a
               c3=gcg**2
               c1=-a1*c3
               c2=(a1**2-1.)*c3
               if(dabs(ang-pi/2.).ge.1.e-2) go to 10
               cosa=c1/c3

               go to 20

   10          continue
               tanx=dsin(ang)/dcos(ang)
               sign=1.
               if(ang.gt.(pi/2.)) sign=-1.
               x2=tanx**2
               cosa=(c1*x2+sign*dsqrt(1.-c2*x2))/(1.+c3*x2)
               if(angle.ge.179.) cosa=-1.
               if(angle.le.0.1) cosa=1.

   20          continue

               angzlm=acos(cosa)

               angzlm=angzlm*180./pi

               if(dabs(ang-pi/2.).ge.0.02) go to 30

               ang1=ang-0.1
               ang2=ang+0.1
               cosx1=dcos(ang1)
               cosx2=dcos(ang2)
               x31=cosx1**3
               x32=cosx2**3
               x21=(dsin(ang1)/dcos(ang1))**2
               x22=(dsin(ang2)/dcos(ang2))**2
               conv=1./2.*((-2.*c1/x31+sign*c2/x31
     &              /dsqrt(1.-c2*x21))/(1.+c3*x21)
     &              +2.*c3/x31*(c1*x21
     &              +sign*dsqrt(1.-c2*x21))/(1.+c3*x21)**2
     &              +(-2.*c1/x32-sign*c2/x32/dsqrt(1.-c2*x22))
     &              /(1.+c3*x22)+2.*c3/x32
     &              *(c1*x22+sign*dsqrt(1.-c2*x22))/(1.+c3*x22)**2)

               go to 40

   30          continue

               cosx=dcos(ang)
               x3=cosx**3

               conv=(-2.*c1/x3+sign*c2/x3/dsqrt(1.-c2*x2))/(1.+c3*x2)+
     &           2.*c3/x3*(c1*x2+sign*dsqrt(1.-c2*x2))/(1.+c3*x2)**2

   40          continue

               if(dabs(conv).le.1.e-6) conv=0.
               if(conv.le.0.0) conv=1.
               if(dabs(a-1.).le.0.01.and.ang.ge.(pi/2.-1.e-4)) conv=0.

         end if

*-----------------------------------------------------------------------

         escat = rmass * ( gz - 1. )
     &         / s * ( 1. + a * ( gz - 1. + ( gz + 1. ) * cosa )
     &                    + a**2 )

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine bessel(n,x,xjn)
c ----------------------------------------------------------------------
c     bessel integer functions, see ams-55, handbook of
c     mathematical functions, chapter 9.
c ----------------------------------------------------------------------
      implicit real*8(a-h,o-z)
c
      dimension a(7),b(7),c(7),d(7),e(7),g(7)
      data a/1.0,-2.2499997, 1.2656208,-.3163866, .0444479,
     *  -.0039444, .0002100/
      data b/ .79788456,-.00000077,-.00552740,-.00009512, .00137237,
     *  -.00072805, .00014476/
      data c/-.78539816,-.04166397,-.00003954, .00262573,-.00054125,
     *  -.00029333, .00013558/
      data d/ .5,-.56249985, .21093573,-.03954289, .00443319,
     * -.00031761, .00001109/
      data e/ .79788456, .0000156, .01659667, .00017105,-.00249511,
     *  .00113653,-.00020033/
      data g/-2.35619449, .12499612, .00005650,-.00637879, .00074348,
     *  .00079824,-.00029166/
c ----------------------------------------------------------------------
c
      if(n-1) 25,210,25
   25 if(x-3.0) 50,50,100
   50 if(x-0.) 55,55,60
   55 xjzr=1.
      go to 200
   60 xjzr=a(1)
      do 70 m=2,7
      k=2*(m-1)
   70 xjzr=xjzr+a(m)*(x/3.)**k
      go to 200
  100 fzer=b(1)
      do 110 m=2,7
  110 fzer=fzer+b(m)*(3./x)**(m-1)
      thetz=x
      do 120 m=1,7
  120 thetz=thetz+c(m)*(3./x)**(m-1)
      xjzr=fzer*dcos(thetz)/dsqrt(x)
  200 if(n-1) 201,210,210
  201 xjn=xjzr
      return
  210 if(x-3.) 250,250,300
  250 if(x-1.e-6) 255,255,260
  255 xjne=0.
      go to 400
  260 xjne=d(1)
      do 270 m=2,7
      k=2*(m-1)
  270 xjne=xjne+d(m)*(x/3.)**k
      xjne=x*xjne
      go to 400
  300 fzne=e(1)
      do 310 m=2,7
  310 fzne=fzne+e(m)*(3./x)**(m-1)
      thetn=x
      do 320 m=1,7
  320 thetn=thetn+g(m)*(3./x)**(m-1)
      xjne=fzne*dcos(thetn)/dsqrt(x)
  400 if(n-1) 401,401,408
  408 if(x-1.e-6) 409,409,410
  401 xjn=xjne
      return
  409 xjn=0.
      return
  410 do 420 i=2,n
      xn=i-1
      xjn=2.0*xn*xjne/x-xjzr
      xjzr=xjne
      xjne=xjn
  420 continue
      return
      end


************************************************************************
*                                                                      *
      subroutine sighi(ap,zp,ep,at,zt,signe,sigel,bmax)
*                                                                      *
*        calculates reaction cross-sections of nucleus-nucleus         *
*        choose the models                                             *
*        last modified by K.Niita on 2016/08/12                        *
*                                                                      *
*     common:                                                          *
*       icrhi   : choice of cross section formula                      *
*               : 0; Shen's formula                                    *
*               : 1; NASA's formula                                    *
*               : 2; KUROTAMA model                                    *
*                                                                      *
*       icrdm = 1  : MWO formula for deuteron                          *
*                                                                      *
*     input:                                                           *
*        ap     : mass number of projectile                            *
*        zp     : charge number of projectile                          *
*        ep     : incident nucleus total energy (MeV)                  *
*        at     : mass number of target                                *
*        zt     : charge number of target                              *
*                                                                      *
*     output:                                                          *
*       signe   : reaction cross-section (b)                           *
*       sigel   : elastic cross-section (b)                            *
*       bmax    : correspond impact parameter (fm)                     *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      include 'param.inc'

      common /crshi/  bplus, icrhi, ijudg, imadj, iqmax
      common /ckurotama/  dsck

      dimension att(10), ztt(10)

      common /paraj/ mstz(300), parz(300)

*-----------------------------------------------------------------------
C for [frag data] (S.Hashimoto, 2014.12.08)
      common /fgdata/ ifrgd, ifgm, ifgdf(kvlmax,5), frgfl(kvlmax)
      character frgfl*200

*-----------------------------------------------------------------------

      if ( ifrgd .gt. 0 ) then

         ktypfrg = idnint(zt)*1000000+idnint(at)
         eppap = ep / ap !(MeV -> MeV/nucleon)
         call fragdataXS(ktypfrg,ap,zp,eppap,at,zt,signe,sigel,bmax)

         if ( signe .gt. 0d0 ) return

      end if

*-----------------------------------------------------------------------

         if( mstz(119) .eq. 1 .and. nint(ap) .eq. 2 .and.
     &       nint(zp) .eq. 1 ) then

            call deumino(ep,at,zt,signe,sigel,bmax)

            return

         end if

*-----------------------------------------------------------------------


         if( icrhi .eq. 0 ) then

            call shen(ap,zp,ep,at,zt,signe,sigel,bmax)

         else if( icrhi .eq. 2 ) then

            eppap = ep / ap !(MeV -> MeV/nucleon)
            call kurotama0(ap,zp,eppap,at,zt,dsck,signe,bmax)
            signe = signe / 100.0d0 !(fm^2 -> b)
            sigel = 0.0d0

         else if( icrhi .eq. 3 ) then ! special setting for deuteron

            if ( ( ap.eq.2 .and. zp.eq.1 .and. zt.ge.3 )
     1           .or. ( at.eq.2 .and. zt.eq.1 .and. zp.ge.3 ) ) then

               call deuXS(ap,zp,ep,at,zt,signe,sigel,bmax)

            else

               call nasa(ap,zp,ep,at,zt,signe,sigel,bmax)

            end if

         else

            call nasa(ap,zp,ep,at,zt,signe,sigel,bmax)

         end if

*-----------------------------------------------------------------------
*        for debug to comment out the next return
*-----------------------------------------------------------------------

         return

              signe = signe * 1.1
              sigel = sigel * 1.1
               bmax = bmax * sqrt(1.1)

         return

              bmax3 = 1.15 * ( ap**(1./3.) + at**(1./3.) )
     &              - 0.4 + bplus

             write(6,'(5e13.5)') ep/ap, bmax, bmax3


*-----------------------------------------------------------------------

         ap = 56
         zp = 26

         ap = 20
         zp = 10

         ap = 20
         zp = 10

         at = 208
         zt = 82

         att(1) =  1
         ztt(1) =  1
         att(2) = 12
         ztt(2) =  6
         att(3) = 27
         ztt(3) = 13
         att(4) = 64
         ztt(4) = 29
         att(5) = 119
         ztt(5) =  50
         att(6) = 181
         ztt(6) =  73
         att(7) = 208
         ztt(7) =  82


         write(6,'(''wt:'')')
         write(6,'('' ap ='',i4,'',  zp ='',i3)') nint(ap), nint(zp)
         write(6,'(''e:'')')
         write(6,'(''x: mass number'')')
         write(6,'(''y: bmax (fm)'')')
         write(6,'(''p: ylin xlin afac(0.8)'')')
         write(6,'(''h: x y(SHEN),l0r y(NASA),l0b y(Sato),l0g'',
     &             '' y(bmax),l0'')')

         emin = 100
         emax = 400.0
         nemd = 16
         edel = log( emax / emin ) / ( nemd - 1 )
         edel = ( emax - emin ) / ( nemd - 1 )

         nemd = 1
         nemd = 16

         do i = 1, nemd

           ep = ( ( i - 1 ) * edel + emin ) * ap


           at = 12
           zt = 6


            call shen(ap,zp,ep,at,zt,signe,sigel,bmax)

               sigsh = signe
               bmax1 = bmax

            call nasa(ap,zp,ep,at,zt,signe,sigel,bmax)

               signa = signe
               bmax2 = bmax

            bmax3 = 1.2 * ( ap**(1./3.) + at**(1./3.) ) - 0.5

             write(6,'(3e13.5)') ep/ap, sigsh*1000, signa*1000

         end do

        stop 666

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine shen(ap0,zp0,ep0,at0,zt0,signe,sigel,bmax)
*                                                                      *
*        calculates reaction cross-sections of nucleus-nucleus         *
*        this parametrization is taken from                            *
*        Nucl. Phys. A491 (1989) 130 by SHEN Wen-qing, et.al.          *
*                                                                      *
*        coded by H. Iwase, on 2002/01/10                              *
*                                                                      *
*     input:                                                           *
*        ap     : mass number of projectile                            *
*        zp     : charge number of projectile                          *
*        ep     : incident nucleus total energy (MeV)                  *
*        at     : mass number of target                                *
*        zt     : charge number of target                              *
*                                                                      *
*     output:                                                          *
*       signe   : reaction cross-section (b)                           *
*       sigel   : elastic cross-section (b)                            *
*       bmax    : correspond impact parameter (fm)                     *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      parameter ( pi = 3.1415926535897932d0 )

c----------------------------------------------------------------------
c     this routine requires that at > ap
c----------------------------------------------------------------------

      if( ap0 .gt. at0 ) then

         ep = ep0 / ap0 * at0
         ap = at0
         zp = zt0
         at = ap0
         zt = zp0

      else

         ep = ep0
         ap = ap0
         zp = zp0
         at = at0
         zt = zt0

      end if

*-----------------------------------------------------------------------

      signe = 0.0
      sigel = 0.0
      bmax  = 0.0

      ecm = ep * at / ( ap + at )

      call keibeta(ap,zp,at,zt,bbb)

      if( ecm .le. bbb ) return

      call keirrr(ap,zp,at,zt,ep,rrr)

      signe = 10.0 * pi * rrr**2 * ( 1.0 - bbb / ecm )
      signe = signe / 1000.

      bmax = rrr * sqrt( 1.0 - bbb / ecm )

      return
      end


************************************************************************
*                                                                      *
        subroutine keibeta(ap,zp,at,zt,bbb)
*                                                                      *
*        ap,zp : projectile mass and charge                            *
*        at,zt : target mass and charge                                *
*        bbb   : barrier                                               *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      b   = 1.0
      rrp = (1.12*(ap**(1.0/3.0)))-(0.94*(ap**(-1.0/3.0)))
      rrt = (1.12*(at**(1.0/3.0)))-(0.94*(at**(-1.0/3.0)))
      r   = rrp+rrt+3.2

      bbb = ((1.44*zt*zp)/r)-b*((rrt*rrp)/(rrt+rrp))

      return
      end


************************************************************************
*                                                                      *
      subroutine keirrr(ap,zp,at,zt,ep0,rrr)
*                                                                      *
*        ap,zp : projectile mass and charge                            *
*        at,zt : target mass and charge                                *
*        ep0   : projectile total energy                               *
*        rrr   : reaction radius                                       *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      common /cedata/ ee(55),ce(55)

*-----------------------------------------------------------------------

      r0=1.1
      be=0.176

      ep = ep0 / ap

      do i=1,55
         if(ee(i).ge.ep) then
            ced=ce(i)
            go to 10
         end if
      end do
            ced=ce(55)
 10   continue

      za=1.85*(at**(1.0/3.0))*(ap**(1.0/3.0))
      zb=(at**(1.0/3.0)+ap**(1.0/3.0))

      dd=za/zb
      dda=at**(1.0/3.0)+ap**(1.0/3.0)+dd-ced
      ddb=((at-2.0*zt)*zp)/(at*ap)
      ecm=ep0*at/(ap+at)

      zc=(ecm**(-1.0/3.0))*(ap**(1.0/3.0))*(at**(1.0/3.0))
      zd=(at**(1.0/3.0)+ap**(1.0/3.0))
      ddc=zc/zd
      rrr=dda*r0+ddb+be*ddc

      return
      end


************************************************************************
*                                                                      *
      block data readce
*                                                                      *
*       transparency coefficients for Shen formula                     *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      common /cedata/ ee(55),ce(55)
      data ee/1.71343E+00, 2.33903E+00, 2.92134E+00, 3.77688E+00,
     +        4.35853E+00, 5.28443E+00, 6.31251E+00, 7.42979E+00,
     +        8.78789E+00, 1.01410E+01, 1.11927E+01, 1.27256E+01,
     +        1.46115E+01, 1.59679E+01, 1.72794E+01, 1.88830E+01,
     +        2.02318E+01, 2.13589E+01, 2.28843E+01, 2.36854E+01,
     +        2.54983E+01, 2.70476E+01, 2.92652E+01, 3.14451E+01,
     +        3.49550E+01, 3.77599E+01, 4.07804E+01, 4.40631E+01,
     +        5.02993E+01, 5.53659E+01, 6.06408E+01, 6.51767E+01,
     +        7.17530E+01, 7.90173E+01, 8.95879E+01, 1.02118E+02,
     +        1.16990E+02, 1.34726E+02, 1.52204E+02, 1.78908E+02,
     +        2.07299E+02, 2.31066E+02, 2.66530E+02, 2.99893E+02,
     +        3.29279E+02, 3.56168E+02, 3.92895E+02, 4.29124E+02,
     +        4.27030E+02, 4.75548E+02, 5.37448E+02, 6.22204E+02,
     +        7.13148E+02, 8.13585E+02, 9.23209E+02/
      data ce/2.77909E-02, 4.44954E-02, 5.85097E-02, 8.40984E-02,
     +        1.01228E-01, 1.24064E-01, 1.58599E-01, 1.87341E-01,
     +        2.24817E-01, 2.50693E-01, 2.94175E-01, 3.28836E-01,
     +        3.78048E-01, 4.30302E-01, 4.68003E-01, 5.26087E-01,
     +        5.81306E-01, 6.24900E-01, 6.83035E-01, 7.38341E-01,
     +        8.48942E-01, 9.42087E-01, 1.02060E+00, 1.15135E+00,
     +        1.26087E+00, 1.33189E+00, 1.41181E+00, 1.47392E+00,
     +        1.52465E+00, 1.59636E+00, 1.67101E+00, 1.71886E+00,
     +        1.78464E+00, 1.83855E+00, 1.90151E+00, 1.94670E+00,
     +        1.98597E+00, 2.01341E+00, 2.02590E+00, 2.01784E+00,
     +        1.99782E+00, 1.96277E+00, 1.92790E+00, 1.90774E+00,
     +        1.88745E+00, 1.87895E+00, 1.86759E+00, 1.86211E+00,
     +        1.86208E+00, 1.86264E+00, 1.86030E+00, 1.87293E+00,
     +        1.89440E+00, 1.90695E+00, 1.93727E+00/
      end


************************************************************************
*                                                                      *
      subroutine deumino(ep,at,zt,signe,sigel,bmax)
*                                                                      *
*        calculates reaction cross-sections of deuteron-nucleus        *
*                                                                      *
*        this parametrization is taken from                            *
*        K. Minomo, K. Washiyama, K. Ogata                             *
*        J. Nucl. Sci. Tech. DOI:10.1080/00223131.2016.1213672         *
*                                                                      *
*        coded by K. Niita, on 2016/08/12                              *
*                                                                      *
*     input:                                                           *
*        ep    : incident nucleus total energy (MeV)                  *
*        at     : mass number of target                                *
*        zt     : charge number of target                              *
*                                                                      *
*     output:                                                          *
*       signe   : reaction cross-section (b)                           *
*       sigel   : elastic cross-section (b)                            *
*       bmax    : correspond impact parameter (fm)                     *
*                                                                      *
************************************************************************

      implicit none

      real*8 ep, at, zt, signe, sigel, bmax
      real*8 pi, a1, a2, a3, b1, b2, b3, c1, c2, c3, d1

      parameter (pi = 3.14159265358979d0)
      parameter (a1 =  0.306)   ! (fm)
      parameter (a2 = -0.923)
      parameter (a3 =  590.d0)  ! (MeV)
      parameter (b1 =  1.33)    ! (fm)
      parameter (b2 = -0.112)
      parameter (b3 =  248.d0)  ! (MeV)
      parameter (c1 =  0.00204) ! (fm)
      parameter (c2 = -0.788)
      parameter (c3 =  453.d0)  ! (MeV)
      parameter (d1 =  0.272d0) ! (MeV)

*-----------------------------------------------------------------------

         sigel = 0.0d0

         signe = ( a1 / ( 1.d0 + a2 * exp( - ep / a3 ) )
     &         + b1 / ( 1.d0 + b2 * exp( - ep / b3 ) ) * at**(1.d0/3.d0)
     &         + c1 / ( 1.d0 + c2 * exp( - ep / c3 ) ) * at )**2
     &         * exp( - d1 * zt / ep )

         bmax = sqrt( signe )

         signe = signe * pi / 100.d0

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine nasa(ap0,zp0,ep0,at0,zt0,signe,sigel,bmax)
*                                                                      *
*        calculates reaction cross-sections of nucleus-nucleus         *
*                                                                      *
*        this parametrization is taken from                            *
*        NIM B 117 (1996) 347 by R.K. Tripathi, et.al.,                *
*        NIM B 129 (1997)  11 by R.K. Tripathi, et.al.,,               *
*        and                                                           *
*        NIM B 155 (1999) 349 by R.K. Tripathi, et.al.                 *
*                                                                      *
*        coded by H. Iwase, on 2004/02/06                              *
*                                                                      *
*     input:                                                           *
*        ap     : mass number of projectile                            *
*        zp     : charge number of projectile                          *
*        ep    : incident nucleus total energy (MeV)                  *
*        at     : mass number of target                                *
*        zt     : charge number of target                              *
*                                                                      *
*     output:                                                          *
*       signe   : reaction cross-section (b)                           *
*       sigel   : elastic cross-section (b)                            *
*       bmax    : correspond impact parameter (fm)                     *
*                                                                      *
************************************************************************

      implicit none

      real*8 signe, pi, r0, deltaE, B, Ecm, R, rp, rt, rmsp, rmst,
     $     radius, S, CE, rmsc, rca, E, rap, rat, rac, D, Rc, fact,
     $     G, T1, SL, X1, Xm, sigel, bmax, ep

      real*8 Ap, At, Zp, Zt, Ac
      real*8 Ap0, At0, Zp0, Zt0, ep0
      integer i

      parameter (pi    = 3.14159265358979d0)
      parameter (r0    = 1.1)   ! (fm)
      parameter (Ac    = 12.)

c----------------------------------------------------------------------
      real*8 eep, scal
      integer mstz
      real*8 parz
      common /paraj/ mstz(300), parz(300)

c----------------------------------------------------------------------
c     this routine requires that at > ap
c----------------------------------------------------------------------

      if( ap0 .gt. at0 ) then

         ep = ep0 / ap0 * at0
         ap = at0
         zp = zt0
         at = ap0
         zt = zp0

      else

         ep = ep0
         ap = ap0
         zp = zp0
         at = at0
         zt = zt0

      end if

c----------------------------------------------------------------------

      signe = 0.0
      sigel = 0.0
      bmax  = 0.0

      E    = Ep / Ap
      Ecm  = Ep * At / ( Ap + At )

      rmsp  = radius(Ap,Zp)
      rmst  = radius(At,Zt)

      fact = sqrt(5./3.) ! = 1.29

      rp  = fact * rmsp
      rt  = fact * rmst
      rca = fact * 2.471

      rap = 3. * Ap  / ( 4. * pi * rp**3. )
      rat = 3. * At  / ( 4. * pi * rt**3. )
      rac = 3. * 12. / ( 4. * pi * rca**3. )

c---------------------------------------------------- D selection start

      if(ap.eq.1.and.zp.eq.0)then
c     neutron + X systems
         T1 = 18.
         D  = 1.85 + ( 0.16 ) / ( 1. + exp( (500.-E) / 200. ))

      elseif(ap.eq.1.and.zp.eq.1.and.at.eq.4.and.zt.eq.2)then
c     proton + alpha
         T1 = 40.
         D  = 2.05

      elseif(ap.eq.1.and.zp.eq.1.and.at.eq.3.and.zt.eq.2)then
c     proton + 3He
         T1 = 58.
         D  = 1.70

      elseif(ap.eq.1.and.zp.eq.1.and.at.eq.6.and.zt.eq.3)then
c     proton + 6Li
         T1 = 40.
         D  = 2.05

      elseif(ap.eq.1.and.zp.eq.1.and.at.eq.7.and.zt.eq.3)then
c     proton + 7Li
         T1 = 37.
         D  = 2.15

      elseif((ap.eq.1.and.zp.eq.1).and.At.le.7)then
c     proton + light nucleus
         T1 = 23.
         D  = 1.85 + ( 0.16 ) / ( 1. + exp( (500.-E) / 200. ))

      elseif(ap.eq.1.and.zp.eq.1)then
c     proton + others
         T1 = 40
         D  = 2.05

      elseif(ap.eq.2.and.zp.eq.1.and.at.eq.4.and.zt.eq.2)then
c     deuteron + alpha
         T1 = 23
         D  = 1.65 + ( 0.22 ) / ( 1. + exp( (500.-E) / 200.  ))

      elseif(ap.eq.2.and.zp.eq.1)then ! deuteron + X systems
         T1 = 23.
         D  = 1.65 + ( 0.1 ) / ( 1. + exp( (500.-E) / 200.  ))

      elseif(ap.eq.3.and.zp.eq.2 .or. at.eq.3.and.zt.eq.2)then
         T1 = 40.
         D  = 1.55


      elseif((ap.eq.4.and.zp.eq.2).and.(At.eq.4.and.Zt.eq.2))then
c     alhpa + alpha
         T1 = 40.
         G  = 300.
         D = 2.77 - 8.0E-3 * At + 1.8E-5 * At * At
     $        -0.8 / ( 1. + exp( (250.-E) /G ) )

      elseif((ap.eq.4.and.zp.eq.2).and.(zt.eq.4))then
c     alpha + Be
         T1 = 25.
         G  = 300.
         D = 2.77 - 8.0E-3 * At + 1.8E-5 * At * At
     $        -0.8 / ( 1. + exp( (250.-E) /G ) )

      elseif((ap.eq.4.and.zp.eq.2).and.(zt.eq.7))then
c     alpha + N
         T1 = 40.
         G  = 500.
         D = 2.77 - 8.0E-3 * At + 1.8E-5 * At * At
     $        -0.8 / ( 1. + exp( (250.-E) /G ) )

      elseif((ap.eq.4.and.zp.eq.2).and.(zt.eq.13))then
c     alpha + Al
         T1 = 25.
         G  = 300.
         D = 2.77 - 8.0E-3 * At + 1.8E-5 * At * At
     $        -0.8 / ( 1. + exp( (250.-E) /G ) )

      elseif((ap.eq.4.and.zp.eq.2).and.(zt.eq.26))then
c     alpha + Fe
         T1 = 40.
         G  = 300.
         D = 2.77 - 8.0E-3 * At + 1.8E-5 * At * At
     $        -0.8 / ( 1. + exp( (250.-E) /G ) )

      elseif(ap.eq.4.and.zp.eq.2 .or. at.eq.4.and.zt.eq.2)then
c     alpha + other systems
         T1 = 40.
         G  = 75.
         D = 2.77 - 8.0E-3 * At + 1.8E-5 * At * At
     $        -0.8 / ( 1. + exp( (250.-E) /G ) )

      elseif(zp.eq.3 .or. zt.eq.3)then
c     Li + X
         T1 = 40.
         D = 1.75 * ( rap + rat ) / ( rac + rac ) /3

      else
c     others
         T1 = 40.
         D = 1.75 * ( rap + rat ) / ( rac + rac )

      endif

c---------------------------------------------------- D selection end


c----------------------------------------------------------------------

      CE = D * ( 1 - exp( -E / T1) )  - 0.292 * exp( -E / 792 ) *
     $     cos( 0.229 * E**(0.453) )

      S = ( Ap**(1./3.) * At**(1./3.) ) / ( Ap**(1./3.) + At**(1./3.) )

      deltaE = 1.85 * S + ( 0.16 * S / Ecm **(1./3.) ) - CE + 0.91 * (
     $     At - 2 * Zt ) * Zp / ( At * Ap )


      R = rp + rt + 1.2 * ( Ap**(1./3.) + At**(1./3.) ) / ( Ecm**(1./3.)
     $     )


      B = 1.44 * Zp * Zt / R

c----------------------------------------------------------------------


c--------------------------------------------------- Rc selection start

      if((ap.eq.1.and.zp.eq.1).and.(at.eq.2.and.zt.eq.1))then
         Rc = 13.5              ! p + d

      elseif((ap.eq.1.and.zp.eq.1).and.(at.eq.3.and.zt.eq.2))then
         Rc = 21                ! p + 3He

      elseif((ap.eq.1.and.zp.eq.1).and.(at.eq.4.and.zt.eq.2))then
         Rc = 27                ! p + 4He

      elseif((ap.eq.1.and.zp.eq.1).and.(zt.eq.3))then
         Rc = 2.2               ! p + Li

      elseif((ap.eq.1.and.zp.eq.1).and.(zt.eq.6))then
         Rc = 3.5               ! p + C

      elseif((ap.eq.2.and.zp.eq.1).and.(at.eq.2.and.zt.eq.1))then
         Rc = 13.5              ! d + d

      elseif((ap.eq.2.and.zp.eq.1).and.(at.eq.4.and.zt.eq.2))then
         Rc = 13.5              ! d + 4He

      elseif((ap.eq.2.and.zp.eq.1).and.(zt.eq.6))then
         Rc = 6.0               ! d + C

      elseif((ap.eq.4.and.zp.eq.2).and.(zt.eq.73))then
         Rc = 0.6               ! alpha + Ta

      elseif((ap.eq.4.and.zp.eq.2).and.(zt.eq.79))then
         Rc = 0.6               ! alpha + Au

      else
         Rc = 1.0
      endif

c------------------------------------------------- Rc selection end


c---------------------------------------------- Xm calculation start

      if (ap.eq.1.and.zp.eq.0)then ! neutron + X systems

         SL = 1.2 + 1.6 * ( 1.0 - exp( -E/15. ) )


         if((at.eq.4.and.zt.eq.2) )then
            X1 = 5.2

         else
c        others
            X1 = 2.83 - 3.1E-2 * At + 1.7E-4 * At * At

         endif

         Xm = 1. - X1 * exp( -E / ( X1 * SL ) )

      else
         Xm = 1.0
      endif

c---------------------------------------------- Xm calculation end


c------------------------------------------- signe calculation start

      signe = pi * r0 * r0 * ( Ap**(1./3.) + At**(1./3.) + deltaE )**(2.
     $     ) * ( 1 - Rc * B / Ecm ) * Xm

      if(signe.le.0) then
         signe=0
         return
      endif

      signe = signe * 10        ! convert the unit from [fm**2] to [mb]
      signe = signe / 1000      ! convert the unit from [mb] to [b]

      bmax =  r0 * ( Ap**(1./3.) + At**(1./3.) + deltaE )
     $        * sqrt(1. -  Rc * B / Ecm)

c------------------------------------------- signe calculation end

c----------------------------------------------------------------------
cKN 2015/12/08 for test
c----------------------------------------------------------------------

      if( mstz(199) .ne. 0 .and. nint(ap0) .eq. 2 ) then

         eep = ep0 / ap0

        if( eep .lt. 500. ) then
           scal = ( 1.d0 + 0.45 * (( 500. - eep ) / 500.)**2 )
        else
           scal = 1.d0
        end if

           signe = signe * scal

      end if

c----------------------------------------------------------------------

      end


c----------------------------------------------------------------------
      function radius(a,z)
c----------------------------------------------------------------------
c Purpose   : to obtain the "r_rms,i" values
c References: Atomic Data adn Nuclear Data Tables 36, 495-536 (1987)
c             and
c             NIM-B 152(1999)425-431
c             by H.Iwase Fri Feb  6 2004
c----------------------------------------------------------------------

      implicit none

      real*8 rms(300,2)
      real*8 a, z, radius
      integer i
      integer irnm

c----------------------------------------------------------------------
c     these data are refered from
c     Atomic Data adn Nuclear Data Tables 36, 495-536 (1987)
c     data from page 503 to 510 are used
c----------------------------------------------------------------------

       data irnm / 135 /

       data rms(  1,1), rms(  1,2) /     1, 0.3407/
       data rms(  2,1), rms(  2,2) /  1001, 0.8507/
       data rms(  3,1), rms(  3,2) /  1002, 2.1055/
       data rms(  4,1), rms(  4,2) /  1003, 1.7200/
       data rms(  5,1), rms(  5,2) /  2003, 1.8990/
       data rms(  6,1), rms(  6,2) /  2004, 1.6810/
       data rms(  7,1), rms(  7,2) /  3006, 2.5567/
       data rms(  8,1), rms(  8,2) /  3007, 2.4000/
       data rms(  9,1), rms(  9,2) /  4009, 2.5095/
       data rms( 10,1), rms( 10,2) /  5010, 2.4500/
       data rms( 11,1), rms( 11,2) /  5011, 2.3950/
       data rms( 12,1), rms( 12,2) /  6012, 2.4690/
       data rms( 13,1), rms( 13,2) /  6013, 2.4400/
       data rms( 14,1), rms( 14,2) /  6014, 2.5600/
       data rms( 15,1), rms( 15,2) /  7014, 2.5480/
       data rms( 16,1), rms( 16,2) /  7015, 2.6537/
       data rms( 17,1), rms( 17,2) /  8016, 2.7283/
       data rms( 18,1), rms( 18,2) /  8017, 2.6620/
       data rms( 19,1), rms( 19,2) /  8018, 2.7270/
       data rms( 20,1), rms( 20,2) /  9019, 2.9000/
       data rms( 21,1), rms( 21,2) / 10020, 3.0120/
       data rms( 22,1), rms( 22,2) / 10022, 2.9690/
       data rms( 23,1), rms( 23,2) / 11023, 2.9400/
       data rms( 24,1), rms( 24,2) / 12024, 3.0467/
       data rms( 25,1), rms( 25,2) / 12025, 3.0565/
       data rms( 26,1), rms( 26,2) / 12026, 3.0600/
       data rms( 27,1), rms( 27,2) / 13027, 3.0483/
       data rms( 28,1), rms( 28,2) / 14028, 3.1140/
       data rms( 29,1), rms( 29,2) / 14029, 3.1045/
       data rms( 30,1), rms( 30,2) / 14030, 3.1760/
       data rms( 31,1), rms( 31,2) / 15031, 3.1880/
       data rms( 32,1), rms( 32,2) / 16032, 3.2423/
       data rms( 33,1), rms( 33,2) / 16034, 3.2810/
       data rms( 34,1), rms( 34,2) / 17035, 3.3880/
       data rms( 35,1), rms( 35,2) / 16036, 3.2780/
       data rms( 36,1), rms( 36,2) / 18036, 3.3270/
       data rms( 37,1), rms( 37,2) / 17037, 3.3840/
       data rms( 38,1), rms( 38,2) / 19039, 3.4040/
       data rms( 39,1), rms( 39,2) / 18040, 3.4320/
       data rms( 40,1), rms( 40,2) / 20040, 3.4703/
       data rms( 41,1), rms( 41,2) / 20048, 3.4605/
       data rms( 42,1), rms( 42,2) / 22048, 3.6550/
       data rms( 43,1), rms( 43,2) / 22050, 3.5730/
       data rms( 44,1), rms( 44,2) / 24050, 3.6690/
       data rms( 45,1), rms( 45,2) / 23051, 3.5975/
       data rms( 46,1), rms( 46,2) / 24052, 3.6467/
       data rms( 47,1), rms( 47,2) / 24053, 3.7260/
       data rms( 48,1), rms( 48,2) / 24054, 3.7127/
       data rms( 49,1), rms( 49,2) / 26054, 3.6957/
       data rms( 50,1), rms( 50,2) / 25055, 3.6800/
       data rms( 51,1), rms( 51,2) / 26056, 3.7503/
       data rms( 52,1), rms( 52,2) / 26058, 3.7750/
       data rms( 53,1), rms( 53,2) / 28058, 3.7683/
       data rms( 54,1), rms( 54,2) / 27059, 3.8130/
       data rms( 55,1), rms( 55,2) / 28060, 3.7953/
       data rms( 56,1), rms( 56,2) / 28061, 3.8060/
       data rms( 57,1), rms( 57,2) / 28062, 3.8263/
       data rms( 58,1), rms( 58,2) / 29063, 3.9187/
       data rms( 59,1), rms( 59,2) / 28064, 3.8673/
       data rms( 60,1), rms( 60,2) / 30064, 3.9370/
       data rms( 61,1), rms( 61,2) / 29065, 3.9440/
       data rms( 62,1), rms( 62,2) / 30066, 3.9580/
       data rms( 63,1), rms( 63,2) / 30068, 3.9667/
       data rms( 64,1), rms( 64,2) / 30070, 4.0077/
       data rms( 65,1), rms( 65,2) / 32070, 4.0565/
       data rms( 66,1), rms( 66,2) / 32072, 4.0550/
       data rms( 67,1), rms( 67,2) / 32074, 4.0750/
       data rms( 68,1), rms( 68,2) / 32076, 4.0810/
       data rms( 69,1), rms( 69,2) / 38088, 4.2060/
       data rms( 70,1), rms( 70,2) / 39089, 4.2500/
       data rms( 71,1), rms( 71,2) / 40090, 4.2707/
       data rms( 72,1), rms( 72,2) / 40091, 4.3090/
       data rms( 73,1), rms( 73,2) / 40092, 4.2970/
       data rms( 74,1), rms( 74,2) / 42092, 4.3047/
       data rms( 75,1), rms( 75,2) / 41093, 4.3205/
       data rms( 76,1), rms( 76,2) / 40094, 4.3235/
       data rms( 77,1), rms( 77,2) / 42094, 4.3340/
       data rms( 78,1), rms( 78,2) / 40096, 4.3960/
       data rms( 79,1), rms( 79,2) / 42096, 4.3640/
       data rms( 80,1), rms( 80,2) / 42098, 4.3880/
       data rms( 81,1), rms( 81,2) / 42100, 4.4300/
       data rms( 82,1), rms( 82,2) / 46104, 4.4370/
       data rms( 83,1), rms( 83,2) / 46106, 4.4670/
       data rms( 84,1), rms( 84,2) / 46108, 4.5240/
       data rms( 85,1), rms( 85,2) / 46110, 4.5900/
       data rms( 86,1), rms( 86,2) / 48110, 4.5780/
       data rms( 87,1), rms( 87,2) / 48112, 4.6080/
       data rms( 88,1), rms( 88,2) / 50112, 4.6205/
       data rms( 89,1), rms( 89,2) / 48114, 4.6305/
       data rms( 90,1), rms( 90,2) / 50114, 4.6020/
       data rms( 91,1), rms( 91,2) / 49115, 4.6460/
       data rms( 92,1), rms( 92,2) / 48116, 4.6390/
       data rms( 93,1), rms( 93,2) / 50116, 4.6240/
       data rms( 94,1), rms( 94,2) / 50117, 4.6250/
       data rms( 95,1), rms( 95,2) / 50118, 4.6630/
       data rms( 96,1), rms( 96,2) / 50119, 4.6390/
       data rms( 97,1), rms( 97,2) / 50120, 4.6430/
       data rms( 98,1), rms( 98,2) / 50122, 4.6580/
       data rms( 99,1), rms( 99,2) / 51122, 4.6300/
       data rms(100,1), rms(100,2) / 50124, 4.6807/
       data rms(101,1), rms(101,2) / 56138, 4.8360/
       data rms(102,1), rms(102,2) / 57139, 4.8500/
       data rms(103,1), rms(103,2) / 60142, 4.9253/
       data rms(104,1), rms(104,2) / 60144, 4.9260/
       data rms(105,1), rms(105,2) / 62144, 4.9470/
       data rms(106,1), rms(106,2) / 60146, 4.9815/
       data rms(107,1), rms(107,2) / 60148, 5.0020/
       data rms(108,1), rms(108,2) / 62148, 4.9890/
       data rms(109,1), rms(109,2) / 60150, 5.0037/
       data rms(110,1), rms(110,2) / 62150, 5.0450/
       data rms(111,1), rms(111,2) / 62152, 5.0947/
       data rms(112,1), rms(112,2) / 62154, 5.1259/
       data rms(113,1), rms(113,2) / 64154, 5.1240/
       data rms(114,1), rms(114,2) / 64156, 5.0680/
       data rms(115,1), rms(115,2) / 64158, 5.1720/
       data rms(116,1), rms(116,2) / 67165, 5.2100/
       data rms(117,1), rms(117,2) / 68166, 5.2593/
       data rms(118,1), rms(118,2) / 70174, 5.4100/
       data rms(119,1), rms(119,2) / 71175, 5.3700/
       data rms(120,1), rms(120,2) / 70176, 5.3790/
       data rms(121,1), rms(121,2) / 73181, 5.4800/
       data rms(122,1), rms(122,2) / 74184, 5.4200/
       data rms(123,1), rms(123,2) / 74186, 5.4000/
       data rms(124,1), rms(124,2) / 76192, 5.4130/
       data rms(125,1), rms(125,2) / 78196, 5.3800/
       data rms(126,1), rms(126,2) / 79197, 5.3000/
       data rms(127,1), rms(127,2) / 81203, 5.4630/
       data rms(128,1), rms(128,2) / 82204, 5.4790/
       data rms(129,1), rms(129,2) / 81205, 5.4745/
       data rms(130,1), rms(130,2) / 82206, 5.4963/
       data rms(131,1), rms(131,2) / 82207, 5.5050/
       data rms(132,1), rms(132,2) / 82208, 5.5016/
       data rms(133,1), rms(133,2) / 83209, 5.5167/
       data rms(134,1), rms(134,2) / 90232, 5.7087/
       data rms(135,1), rms(135,2) / 92238, 5.8470/

c----------------------------------------------------------------------

      do i = 1, irnm

         if( nint( z*1000 + a ) .eq. nint( rms(i,1) ) ) then

            radius = rms(i,2)
            return

         end if

      end do

c----------------------------------------------------------------------
c     if there is no data, a fitted values is used.
c     the formula is cited from NIM-B 152(1999)425-431
c----------------------------------------------------------------------

      radius = ( 0.84 * a**(1./3.) + 0.55 )

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      block data decayd
*                                                                      *
*       life time of decay particles                                   *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      common /dcyd/ rtdcy(50), kfdcy(50), ndcy

      data ndcy / 50 /
      data ( kfdcy(i), rtdcy(i), i = 1, 50 ) /
     &       211,   2.6029d-8,      ! pion+
     &       111,   1.0d-30,        ! pion0
     &      -211,   2.6029d-8,      ! pion-
     &       -13,   2.19703d-6,     ! muon+
     &        13,   2.19703d-6,     ! muon-
     &       321,   1.2371d-8,      ! kaon+
     &       311,   0.8922d-10,     ! kaon0
     &      -321,   1.2371d-8,      ! kaon+
     &       221,   1.0d-30,        ! eta
     &       331,   1.0d-30,        ! eta'
     &      -311,   0.8922d-10,     ! k0bar
     &      3122,   2.631d-10,      ! Lambda0
     &      3222,   0.799d-10,      ! Sigma+
     &      3212,   1.0d-30,        ! Sigma0
     &      3112,   1.479d-10,      ! Sigma-
     &      3322,   2.90d-10,       ! Xi0
     &      3312,   1.639d-10,      ! Xi-
     &      3334,   0.822d-10,      ! Omega-
     &   6000009,   1.265d-1,       ! 9C
     &       -11,   1.0d+30,        ! e+
     &      2112,   886.7d0,        ! neutron, add by T.Sato 2016/03/02
     &       130,   5.116d-8,       ! K_L0, K-long , y.sakaki 2021/08
     &       310,   0.8954d-10,     ! K_S0, K-short, y.sakaki 2021/08
     &        15,   290.3d-15,      !      tau- (1777.0)
     &       -15,   290.3d-15,      !      tau+ (1777.0)
     &       113,   1d-30,          !      rho0 (768.50)
     &       213,   1d-30,          !      rho+ (766.90)
     &      -213,   1d-30,          !      rho- (766.90)
     &       223,   1d-30,          !      omega(781.94)
     &     20213,   1d-30,          !      a_1+ (1230.0)
     &    -20213,   1d-30,          !      a_1- (1230.0)
     &      -313,   1d-30,          ! anti-K*0  (896.10)
     &       313,   1d-30,          !      K*0  (896.10)
     &       323,   1d-30,          !      K*+  (891.60)
     &      -323,   1d-30,          !      K*-  (891.60)
     &       333,   1d-30,          !      phi  (1019.4)
     &       411,   1033d-15,       !      D+   (1869.3)
     &      -411,   1033d-15,       !      D-   (1869.3)
     &       421,   410.3d-15,      !      D0   (1864.5)
     &      -421,   410.3d-15,      ! anti-D0   (1864.5)
     &       431,   504d-15,        !      D_s+ (1968.5)
     &      -431,   504d-15,        !      D_s- (1968.5)
     &       511,   1.519d-12,      !      B0   (5279.2)
     &      -511,   1.519d-12,      ! anti-B0   (5279.2)
     &       521,   1.638d-12,      !      B+   (5278.9)
     &      -521,   1.638d-12,      !      B-   (5278.9)
     &       531,   1.521d-12,      !      B_s0 (5369.3)
     &      -531,   1.521d-12,      ! anti-B_s0 (5369.3)
     &       541,   0.510d-12,      !      B_c+ (6594.0)
     &      -541,   0.510d-12/      !      B_c- (6594.0)
      end


************************************************************************
*                                                                      *
      subroutine deuXS(ap0,zp0,ep0,at0,zt0,signe,sigel,bmax)
*                                                                      *
*     calculation of reaction cross-sections of deuteron-nucleus       *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------

      implicit none

      real*8 ap0, zp0, ep0, at0, zt0, signe, sigel, bmax
      real*8 ep, ap, at, zp, zt, Ecm
      real*8 Rcoul, Bcoul, Eth, Etmp, fact

      real*8 RM, fk, eta
      real*8 func, aaa, bbb
      real*8 rmsp, rmst, factr, rp, rt, radius

      real*8 r0
      parameter (r0    = 1.1)   ! (fm)
      real*8 pi, HC, EC, AC
      DATA pi/3.141592653589793D0/
      DATA HC/197.32705D0/, EC/137.03599D0/, AC/938.91897D0/

*-----------------------------------------------------------------------
*     this routine requires that at > ap
*-----------------------------------------------------------------------

      if( ap0 .gt. at0 ) then

         ep = ep0 / ap0 * at0
         ap = at0
         zp = zt0
         at = ap0
         zt = zp0

      else

         ep = ep0
         ap = ap0
         zp = zp0
         at = at0
         zt = zt0

      end if

*-----------------------------------------------------------------------

      Ecm = ep * at / ( ap + at )

      rmsp  = radius(Ap,Zp)
      rmst  = radius(At,Zt)
      factr = sqrt(5./3.) ! = 1.29
      rp  = factr * rmsp
      rt  = factr * rmst

      Rcoul = rp + rt

      Bcoul = 1.44d0 * zp * zt / Rcoul
      Eth = Bcoul * 0.8d0

      if ( Ecm .ge. Eth ) then

         if ( zt.eq.6 ) then

            call shen(ap,zp,ep,at,zt,signe,sigel,bmax)
            return

         else

            call nasa(ap,zp,ep,at,zt,signe,sigel,bmax)
            return

         end if

      else

         Etmp = Eth * ( ap + at ) / at

         if ( zt.eq.6 ) then

            call shen(ap,zp,Etmp,at,zt,signe,sigel,bmax)

         else

            call nasa(ap,zp,Etmp,at,zt,signe,sigel,bmax)

         end if

         fact = signe

      end if

*-----------------------------------------------------------------------

      bbb = ( 5.3d0 - 0.8d0 ) / ( 82 - 13 ) * ( zt - 13 ) + 0.8d0
      RM = Ap * At / ( Ap + At ) * AC

*-----------------------------------------------------------------------

      fk = dsqrt( 2d0 * RM * Eth ) / HC
      eta = Zp * Zt * RM / ( EC * HC * fk )

      signe = 1d0 / Eth * dexp( -2d0 * pi * eta )

      func = 1d0 / Eth**bbb
      signe = signe * func

      fact = fact / signe

*-----------------------------------------------------------------------

      fk = dsqrt( 2d0 * RM * Ecm ) / HC
      eta = Zp * Zt * RM / ( EC * HC * fk )

      signe = 1d0 / Ecm * dexp( -2d0 * pi * eta )

      func = 1d0 / Ecm**bbb
      signe = signe * func * fact

      if(signe.le.0) then
         signe=0
         return
      endif

      bmax = dsqrt( signe / pi )
      sigel = 0d0

      return

*-----------------------------------------------------------------------

      end

************************************************************************
*                                                                      *
      subroutine fragdataXS(ktypfrg,ap,zp,ep,at,zt,signe,sigel,bmax)
*                                                                      *
*     reaction cross-sections given by [frag data]                     *
*                                                                      *
************************************************************************
      use moddas_fragdata
*-----------------------------------------------------------------------

      implicit real*8 (a-h,o-z)

      include 'param.inc'

      real*8 pi
      parameter ( pi  = 3.1415926535898d0 )

*-----------------------------------------------------------------------

      integer ktypfrg
      real*8 ep, at, zt, signe, sigel, bmax
      integer ifg, iene, nei, kne, kxs
      integer iproj, itarg, izpfrg, iapfrg, iztfrg, iatfrg
      real*8 ee1, ee2, eed
      character frgfl*200

*-----------------------------------------------------------------------

      integer ifrgd, ifgdf, ifgm

      common /fgdata/ ifrgd, ifgm, ifgdf(kvlmax,5), frgfl(kvlmax)

*-----------------------------------------------------------------------

      signe = 0d0

*-----------------------------------------------------------------------

      do ifg = 1, ifrgd

         iconsistent = 0

         iproj = kftp(ifgdf(ifg,2))
         itarg = kftp(ifgdf(ifg,3))

*-----------------------------------------------------------------------
C check consistency

C case that projectile written in [frag data] is not nucleus
       if ( iproj .ge. 1 .and. iproj .le. 14 ) then

        if( ifgdf(ifg,2) .eq. ktypfrg ) then

         if ( ifgdf(ifg,3) .eq. idnint(zt)*1000000+idnint(at) ) then
            iconsistent = 1

         else if ( ifgdf(ifg,3) .eq. 2212 ) then

          if ( idnint(zt)*1000000+idnint(at) .eq. 1000001 ) then
            iconsistent = 1
          end if

         end if

        end if

C case that projectile written in [frag data] is nucleus
       else if ( iproj .ge. 15 .and. iproj .le. 19 ) then

C           and target is proton
        if ( itarg .eq. 1 ) then

         if( ktypfrg .eq. 2212 .and.
     &       ifgdf(ifg,2) .eq. idnint(zt)*1000000+idnint(at) ) then
            iconsistent = 1
         end if

C           and target is also nucleus
        else if ( itarg .ge. 15 .and. itarg .le. 19 ) then

         if ( ( ifgdf(ifg,2) .eq. idnint(zp)*1000000+idnint(ap) .and.
     &          ifgdf(ifg,3) .eq. idnint(zt)*1000000+idnint(at) ) .or.
     &        ( ifgdf(ifg,2) .eq. idnint(zt)*1000000+idnint(at) .and.
     &        ifgdf(ifg,3) .eq. idnint(zp)*1000000+idnint(ap) ) ) then
            iconsistent = 1
         end if

        end if

       end if


       if( iconsistent .eq. 1 ) then

          nei = ifgs01_nei(ifg)
          kne = ifgs02_kne(ifg)
          kxs = ifgs03_kxs(ifg)

        do iene = 1, nei

           ee1 = frgne(kne+iene)
           ee2 = frgne(kne+iene+1)
           eed = ee2 - ee1

         if( ep .ge. ee1 .and. ep .le. ee2 ) then

            signe = ( ( ee2 - ep ) * frgxs(kxs+iene)
     &           - ( ee1 - ep ) * frgxs(kxs+iene+1) )
     &           / eed

            go to 100

         else if( ep .lt. ee1 ) then

          if ( ifgdf(ifg,1) .eq. 4 ) then

             signe = frgxs(kxs+1)
             go to 100

          end if

         else if( ep .gt. ee2 ) then

          if ( ifgdf(ifg,1) .eq. 4 ) then

             signe = frgxs(kxs+nei+1)
             go to 100

          end if


         end if

        end do

       end if

      end do

 100  continue

*-----------------------------------------------------------------------

      if( signe .lt. 0d0 ) signe = 0d0

      signe = signe / 1000d0 ! mb -> b
      bmax = dsqrt( signe * 100d0 / pi ) ! fm
      sigel = 0d0

      return

*-----------------------------------------------------------------------

      end


************************************************************************
*                                                                      *
      subroutine txsmodel(iproj,itarg,ein,signe)
*                                                                      *
*     calculation of total reaction cross-sections using models        *
*     modified by S.Hashimoto on 2016/09/24                            *
*                                                                      *
*        calculates total reaction cross sections                      *
*                                                                      *
*     input:                                                           *
*        iproj  : incident particle                                    *
*        itarg  : target nucleus                                       *
*        ein    : incident energy (MeV/u)                              *
*                                                                      *
*     output:                                                          *
*       signe   : nonelastic cross-section (b)                         *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------

      implicit none

      integer iproj, itarg
      real*8  ein, signe
      integer iaproj, izproj, iatarg, iztarg
      real*8  aproj, zproj, atarg, ztarg, epin
      integer incp
      real*8  sigt, sigel, bmax, dsck

*-----------------------------------------------------------------------

      common /paraj/ mstz(300), parz(300)
      integer mstz
      real*8  parz

      common /crshi/  bplus, icrhi, ijudg, imadj, iqmax
      real*8  bplus
      integer icrhi, ijudg, imadj, iqmax

*-----------------------------------------------------------------------

      iatarg = iabs( itarg - itarg / 1000000 * 1000000 )
      iztarg = iabs( itarg / 1000000 )
      atarg = dble(iatarg)
      ztarg = dble(iztarg)

*-----------------------------------------------------------------------
*     for proton or neutron
*-----------------------------------------------------------------------

      if ( iproj.eq.2212 .or. iproj.eq.2112 ) then

       if ( iproj.eq.2212 ) then ! proton
        incp = 1
       else if ( iproj.eq.2112 ) then ! neutron
        incp = 2
       end if

       call sigrc(incp,ein,iatarg,iztarg,sigt,signe,sigel)

       if ( signe .gt. 0d0 ) return

      end if

*-----------------------------------------------------------------------
*     for deuteron using MWO formula
*-----------------------------------------------------------------------

      if ( iproj.eq.1000002 .and. mstz(119).eq.1 ) then

       epin = ein * 2d0

       call deumino(epin,atarg,ztarg,signe,sigel,bmax)

       if ( signe .gt. 0d0 ) return

      end if

*-----------------------------------------------------------------------
*     for nucleus
*-----------------------------------------------------------------------

      if ( iproj .ge. 1000002 ) then

       iaproj = iabs( iproj - iproj / 1000000 * 1000000 )
       izproj = iabs( iproj / 1000000 )
       aproj = dble(iaproj)
       zproj = dble(izproj)

       epin = ein * aproj

       if( icrhi .eq. 0 ) then

        call shen(aproj,zproj,epin,atarg,ztarg,signe,sigel,bmax)

       else if( icrhi .eq. 2 ) then

        call kurotama0(aproj,zproj,ein,atarg,ztarg,dsck,signe,bmax)
        signe = signe / 100.0d0 !(fm^2 -> b)
        sigel = 0.0d0

       else if( icrhi .eq. 3 ) then ! special setting for deuteron

        if ( ( aproj.eq.2 .and. zproj.eq.1 .and. ztarg.ge.3 ) ) then

         call deuXS(aproj,zproj,epin,atarg,ztarg,signe,sigel,bmax)

        else

         call nasa(aproj,zproj,epin,atarg,ztarg,signe,sigel,bmax)

        end if

       else

        call nasa(aproj,zproj,epin,atarg,ztarg,signe,sigel,bmax)

       end if

       if ( signe .gt. 0d0 ) return

      end if

*-----------------------------------------------------------------------
*     for the others
*-----------------------------------------------------------------------

      if ( signe .le. 0d0 ) signe = 1d-99

      return

*-----------------------------------------------------------------------

      end


************************************************************************
*                                                                      *
      subroutine pionXS(pichg,ep,at,zt,sigtot,signe,sigel,bmax)
*                                                                      *
*     calculation of reaction cross-sections for pion-nucleus          *
*     modified by S.Hashimoto on 2016/06/22                            *
*                                                                      *
*        calculates total, nonelastic and elastic cross-sections       *
*                                                                      *
*     input:                                                           *
*        pichg  : =-1: pi-, 0: pi0, 1: pi+                             *
*        ep     : incident pion energy (MeV)                           *
*        at     : target mass number                                   *
*        zt     : target charge number                                 *
*                                                                      *
*     output:                                                          *
*       sigtot  : total cross-section (b)                              *
*       signe   : nonelastic cross-section (b)                         *
*       signel  : sigt-sigr=elastic scattering cross-section (b)       *
*       bmax    : impact parameter calculated from signe (fm)          *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------

      use NGSDATAMOD, only : bindeg
      implicit none

      include 'param-physcnst.inc'

      real*8 pichg,ep,at,zt,sigtot,signe,sigel,bmax
      real*8 Aproj, Atarg, Tplab, Eplab, Pplab, beta, gamma
      real*8 PpCM, EpCM, EtCM, totEcm, fk, Rat, blambda, fLzero
      real*8 bene
      integer iatarg, iptarg, intarg
      real*8 srt, pr
      real*8 gamma_1, gamma_0, Gamma_el, Gamma_t
      real*8 Ezero, delE
      real*8 sigtot1, signe1, sigel1, sigtot2, signe2, sigel2, ratio_tr
      real*8 sigtotp, signep, sigelp, sigtotm, signem, sigelm
      real*8 cut1, cut2, ecut, width, xene, scale_el

*-----------------------------------------------------------------------

*-----------------------------------------------------------------------
*     for pi0 induced reaction
*     cross sections are given as 0 because pi0 is a short-lived particle
*-----------------------------------------------------------------------

      if ( dabs(pichg) .lt. 0.5d0 ) then

         sigtot = 0d0
         sigel = 0d0
         signe = 0d0
         bmax = 0d0

         return

      end if


*-----------------------------------------------------------------------
*     kinematic parameters
*-----------------------------------------------------------------------

      if ( pichg .gt. 0d0 ) then
         Aproj = rstms(3)*1d3
      else if ( pichg .lt. 0d0 ) then
         Aproj = rstms(5)*1d3
      else
         Aproj = rstms(4)*1d3
      end if

      iatarg = idnint(at)
      iptarg = idnint(zt)
      intarg = idnint(at-zt)
      bene = bindeg(iptarg,intarg)
      Atarg = rstms(1)*1d3 * iptarg + rstms(2)*1d3 * intarg - bene

      Tplab = ep
      Eplab = Tplab+Aproj
      Pplab = dsqrt(Tplab*(Tplab+2d0*Aproj))
      beta = Pplab/(Eplab+Atarg)
      gamma = 1d0/dsqrt(1d0-beta**2)
      PpCM = gamma*Atarg/(Eplab+Atarg)*Pplab
      EpCM = gamma*(Aproj**2+Atarg*Eplab)/(Eplab+Atarg)
      EtCM = gamma*Atarg

      totEcm = EpCM + EtCM
      fk = PpCM / physc(3)

*     no reaction
      if( totEcm .lt. Atarg+Aproj+0.1d0) then
        sigtot = 0d0
        sigel = 0d0
        signe = 0d0
        bmax = 0d0
        return
      end if


*-----------------------------------------------------------------------
*     for proton target
*-----------------------------------------------------------------------
c...srt     : invariant mass in GeV                     (input)
c...pr      : relative momentum(GeV/c) in two-body c.m. (input)

      if ( iptarg.eq.1 .and. iatarg.eq.1 ) then !for proton

         srt = totEcm*1d-3      ! MeV -> GeV
         pr = PpCM*1d-3         ! MeV/c -> GeV/c
         call jamxpin(idnint(pichg),iptarg,srt,pr,sigtot,sigel)

         sigtot = sigtot *1d-3 ! mb -> b
         sigel = sigel *1d-3 ! mb -> b

         signe = sigtot - sigel
         bmax = dsqrt( signe *1d2 / physc(1) ) ! fm

         return

      end if


*-----------------------------------------------------------------------
*     total, elastic, and reaction cross sections for resonance region
*     referred to D.B.Ion et al.[Rom. Journ. Phys. 54, 601 (2009)]
*-----------------------------------------------------------------------

      if ( (iptarg.eq.1 .and. iatarg.eq.2)
     &     .or. (iptarg.eq.1 .and. iatarg.eq.3) ) then !for d and t
         Rat = 2.8d0
      else if ( (iptarg.eq.2 .and. iatarg.eq.3)
     &     .or. (iptarg.eq.2 .and. iatarg.eq.4) ) then !for 3He and 4He
         Rat = 2.21d0
      else if ( iptarg.eq.3 .and. iatarg.eq.6 ) then !for 6Li
         Rat = 3.23d0
      else if ( iptarg.eq.3 .and. iatarg.eq.7 ) then !for 7Li
         Rat = 3.15d0
      else if ( iptarg.eq.4 .and. iatarg.eq.9 ) then !for 9Be
         Rat = 3.12d0
      else
         Rat = 0.10854E+01*at**(1d0/3d0) + 0.57253E+00
      end if

      blambda = 1d0 / fk
      fLzero = fk * Rat

C for pion+
      if ( pichg .gt. 0d0 ) then

         if ( (iptarg.eq.1 .and. iatarg.eq.2)
     &        .or. (iptarg.eq.1 .and. iatarg.eq.3) ) then !for d and t
            gamma_0 = 0.58107E+00
            gamma_1 = 0.81284E-01
            Gamma_t = 0.13637E+03
         else if ( (iptarg.eq.2 .and. iatarg.eq.3)
     &           .or. (iptarg.eq.2 .and. iatarg.eq.4) ) then !for 3He and 4He
            gamma_0 = 0.78000E+00
            gamma_1 = 0.24254E+00
            Gamma_t = 0.20678E+03
         else if ( iptarg.eq.3 .and. iatarg.eq.6 ) then !for 6Li
            gamma_0 = 0.80427E+00
            gamma_1 = 0.20367E+00
            Gamma_t = 0.23052E+03
         else if ( iptarg.eq.3 .and. iatarg.eq.7 ) then !for 7Li
            gamma_0 = 0.81653E+00
            gamma_1 = 0.24080E+00
            Gamma_t = 0.25446E+03
         else if ( iptarg.eq.4 .and. iatarg.eq.9 ) then !for 9Be
            gamma_0 = 0.11403E+01
            gamma_1 = 0.32713E+00
            Gamma_t = 0.29625E+03
         else
            gamma_0 = 0.43090E-01*at**(2d0/3d0)
     &           + 0.60398E+00*at**(1d0/3d0) -0.47695E+00
            gamma_1 = 0.57599E-01*at**(2d0/3d0)
     &           -0.67610E-02*at**(1d0/3d0) + 0.15836E+00
            Gamma_t = 0.31570E+02*at**(2d0/3d0)
     &           -0.49310E+02*at**(1d0/3d0) + 0.26864E+03
         end if

         scale_el = 0.40823E-01*at**(2d0/3d0)
     &        -0.40049E+00*at**(1d0/3d0) + 0.16752E+01

C for pion-
      else if ( pichg .lt. 0d0 ) then

         if ( (iptarg.eq.1 .and. iatarg.eq.2)
     &        .or. (iptarg.eq.1 .and. iatarg.eq.3) ) then !for d and t
            gamma_0 = 0.58107E+00
            gamma_1 = 0.81284E-01
            Gamma_t = 0.13637E+03
         else if ( (iptarg.eq.2 .and. iatarg.eq.3)
     &           .or. (iptarg.eq.2 .and. iatarg.eq.4) ) then !for 3He and 4He
            gamma_0 = 0.78000E+00
            gamma_1 = 0.24254E+00
            Gamma_t = 0.20678E+03
         else if ( iptarg.eq.3 .and. iatarg.eq.6 ) then !for 6Li
            gamma_0 = 0.10061E+01
            gamma_1 = 0.21359E+00
            Gamma_t = 0.23639E+03
         else if ( iptarg.eq.3 .and. iatarg.eq.7 ) then !for 7Li
            gamma_0 = 0.13667E+01
            gamma_1 = 0.28338E+00
            Gamma_t = 0.25615E+03
         else if ( iptarg.eq.4 .and. iatarg.eq.9 ) then !for 9Be
            gamma_0 = 0.12134E+01
            gamma_1 = 0.34981E+00
            Gamma_t = 0.27329E+03
         else
            gamma_0 = 0.59714E-01*at**(2d0/3d0)
     &           + 0.22185E+00*at**(1d0/3d0) + 0.93526E+00
            gamma_1 = -0.49579E-02*at**(2d0/3d0)
     &           + 0.45493E+00*at**(1d0/3d0) -0.48445E+00
            Gamma_t = -0.13142E+02*at**(2d0/3d0)
     &           + 0.24062E+03*at**(1d0/3d0) -0.10202E+03
         end if

         scale_el = -0.49187E-02*at**(2d0/3d0)
     &        -0.68791E-01*at**(1d0/3d0) + 0.12043E+01

C for pion0
      else

         gamma_1 = 0d0
         gamma_0 = 0d0
         Gamma_t = 0d0
         scale_el = 0d0

      end if

      Gamma_el = fk*physc(3) * gamma_1
      Ezero = Atarg + 1236d0 - physc(5)*1d3
      delE = Ezero - totEcm

      sigtot1 = physc(1) * blambda**2 * (fLzero+1d0)**2
     &     * Gamma_el *(Gamma_t - gamma_0 *delE )
     &     /( delE**2 +(Gamma_t - gamma_0 *delE)**2/4d0 )
      sigtot1 = sigtot1 * 10d0 ! fm^2 -> mb
      if ( sigtot1 .lt. 0d0 ) sigtot1 = 0d0

      sigel1 = physc(1) * blambda**2 * (fLzero+1d0)**2
     &     * Gamma_el**2
     &     /( delE**2 +(Gamma_t - gamma_0 *delE)**2/4d0 )
      sigel1 = sigel1 * 10d0 ! fm^2 -> mb
      sigel1 = sigel1 * scale_el
      if ( sigel1 .lt. 0d0 ) sigel1 = 0d0
      if ( sigel1 .gt. sigtot1 ) sigel1 = sigtot1

      signe1 = sigtot1 - sigel1

      ratio_tr = -0.26296E-01*at**(2d0/3d0)
     &     + 0.34726E+00*at**(1d0/3d0) + 0.64584E+00


*-----------------------------------------------------------------------
*     total, elastic, and reaction cross sections for high energy region
*-----------------------------------------------------------------------
*     for deuteron target using jamxpin
      if ( iptarg.eq.1 .and. iatarg.eq.2 ) then !for deuteron

         Atarg = rstms(1)*1d3
         beta = Pplab/(Eplab+Atarg)
         gamma = 1d0/dsqrt(1d0-beta**2)
         PpCM = gamma*Atarg/(Eplab+Atarg)*Pplab
         EpCM = gamma*(Aproj**2+Atarg*Eplab)/(Eplab+Atarg)
         EtCM = gamma*Atarg
         totEcm = EpCM + EtCM
         srt = totEcm*1d-3      ! MeV -> GeV
         pr = PpCM*1d-3         ! MeV/c -> GeV/c

         call jamxpin(1,1,srt,pr,sigtotp,sigelp)
         signep = sigtotp - sigelp

         call jamxpin(-1,1,srt,pr,sigtotm,sigelm)
         signem = sigtotm - sigelm

         sigtot2 = sigtotp + sigtotm
         sigel2 = sigelp + sigelm
         signe2 = signep + signem

         signe1 = sigtot1 * signe2 / sigtot2


*     for the other targets (except proton and deuteron)
*     referred to B.W.Allardyce et al.[Nulc. Phys. A209, 1 (1973)]
      else

         signe2 = 45.0d0*at**0.689d0
         sigtot2 = signe2 * ratio_tr

      end if


*-----------------------------------------------------------------------
*     combination of two regions using Fermi funcions
*-----------------------------------------------------------------------

      ecut = 500d0
      width = 100d0
      xene = (ep-ecut)/width

      if ( xene .lt. 100d0*dlog(10d0) ) then
         cut1 = 1.0d0/(1.0d0 + dexp((ep - ecut)/width))
      else
         cut1 = 0d0
      end if
      if ( xene .gt. -100d0*dlog(10d0) ) then
         cut2 = 1.0d0/(1.0d0 + dexp((-ep + ecut)/width))
      else
         cut2 = 0d0
      end if

      signe = cut1 * signe1 + cut2 * signe2
      sigtot = cut1 * sigtot1 + cut2 * sigtot2
      sigel = sigtot - signe
      if ( sigel .lt. 0d0) sigel = 0d0


*-----------------------------------------------------------------------

      sigtot = sigtot *1d-3     ! mb -> b
      signe  = signe *1d-3      ! mb -> b
      sigel  = sigel *1d-3      ! mb -> b

      bmax = dsqrt( signe *1d2 / physc(1) ) ! fm

      return

*-----------------------------------------------------------------------

      end

************************************************************************
      subroutine getpseudocs(einpu,iap,izp,iat,izt,pseudone,pseudoel) ! frtati 2022/12/20
************************************************************************
      use membankmod

      implicit real*8(a-h,o-z)
      parameter ( epseudomx = 3000.d0 )
      parameter ( nemd = 20, nemd_db = 100 )
      common /qparm/  ielas,icasc,iqstep,lvlopt,igamma
      common /eparm/  esmax, esmin, emin(20)

 100  pseudone = 0.d0
      pseudoel = 0.d0
      epeak = 0.d0

      if( ldodbps ) then
        do idb = 1, ndbcurrent
          if( (izp.eq.dbpscs_iza(1,idb) .and.
     &         iap.eq.dbpscs_iza(2,idb) .and.
     &         izt.eq.dbpscs_iza(3,idb) .and.
     &         iat.eq.dbpscs_iza(4,idb)) .or.
     &        (izp.eq.dbpscs_iza(3,idb) .and.
     &         iap.eq.dbpscs_iza(4,idb) .and.
     &         izt.eq.dbpscs_iza(1,idb) .and.
     &         iat.eq.dbpscs_iza(2,idb)) ) then
            if( einpu.lt.dbpscs_epeak(idb) ) then
              if( iat.eq.1 ) then
                call sigrc(1,einpu,iap,izp,sigt,signe,sigel)
                pseudone = signe
                pseudoel = sigel
              else
                call sighi(dble(iap),dble(izp),einpu*dble(iap),
     &                     dble(iat),dble(izt),signe,sigel,bmax)
                if( ielas.le.2 ) sigel = 0.d0
                pseudone = signe
                pseudoel = sigel
              end if
            else
              pseudone = dbpscs_signe(idb)
              pseudoel = dbpscs_sigel(idb)
            end if
            return
          end if
        end do
      end if

      if( iap.lt.0 ) then ! pion-nucleus
        if( laddbps ) then
          edel = log( epseudomx / esmin ) / ( nemd_db - 1 )
          do k = 1, nemd_db
            ei = exp( ( k - 1 ) * edel ) * esmin
            call pionXS(dble(izp),ei,dble(iat),dble(izt),sigt,signe,
     &                  sigel,bmax)
            if( iat.eq.1 )then
              if( sigt.gt.(pseudone+pseudoel) ) then
                pseudone = signe
                pseudoel = sigel
              end if
            else
              if( signe.gt.pseudone ) then
                pseudone = signe
                pseudoel = sigel
              end if
            end if
          end do
!$omp critical (dbpseudo)
          ndbcurrent = ndbcurrent + 1
          dbpscs_iza(1,ndbcurrent) = izp
          dbpscs_iza(2,ndbcurrent) = iap
          dbpscs_iza(3,ndbcurrent) = izt
          dbpscs_iza(4,ndbcurrent) = iat
          dbpscs_signe(ndbcurrent) = pseudone
          dbpscs_sigel(ndbcurrent) = pseudoel
          if( ndbcurrent.eq.mdbpseud ) then
            laddbps = .false.
          end if
!$omp end critical (dbpseudo)
        else
          edel = log( einpu / esmin ) / ( nemd - 1 )
          do k = 1, nemd
            ei = exp( ( k - 1 ) * edel ) * esmin
            call pionXS(dble(izp),ei,dble(iat),dble(izt),sigt,signe,
     &                  sigel,bmax)
            if( iat.eq.1 )then
              if( sigt.gt.(pseudone+pseudoel) ) then
                pseudone = signe
                pseudoel = sigel
              end if
            else
              if( signe.gt.pseudone ) then
                pseudone = signe
                pseudoel = sigel
              end if
            end if
          end do
        end if
      else if( iat.eq.1 ) then ! nucleus-proton
        if( laddbps ) then
          edel = log( epseudomx / esmin ) / ( nemd_db - 1 )
          do k = 1, nemd_db
            ei = exp( ( k - 1 ) * edel ) * esmin
            call sigrc(1,ei,iap,izp,sigt,signe,sigel)
            if( signe.gt.pseudone ) then
              pseudone = signe
              pseudoel = sigel
              epeak = ei
            end if
          end do
!$omp critical (dbpseudo)
          ndbcurrent = ndbcurrent + 1
          dbpscs_iza(1,ndbcurrent) = izp
          dbpscs_iza(2,ndbcurrent) = iap
          dbpscs_iza(3,ndbcurrent) = 1
          dbpscs_iza(4,ndbcurrent) = 1
          dbpscs_signe(ndbcurrent) = pseudone
          dbpscs_sigel(ndbcurrent) = pseudoel
          if( epeak.lt.200.d0 ) then
            dbpscs_epeak(ndbcurrent) = epeak
          end if
          if( ndbcurrent.eq.mdbpseud ) then
            laddbps = .false.
          end if
          if( einpu.lt.epeak ) then ! frtati 2023/4/6
            call sigrc(1,einpu,iap,izp,sigt,signe,sigel)
            pseudone = signe
            pseudoel = sigel
          end if
!$omp end critical (dbpseudo)
        else
          edel = log( einpu / esmin ) / ( nemd - 1 )
          do k = 1, nemd
            ei = exp( ( k - 1 ) * edel ) * esmin
            call sigrc(1,ei,iap,izp,sigt,signe,sigel)
            if( signe.gt.pseudone ) then
              pseudone = signe
              pseudoel = sigel
            end if
          end do
        end if
      else ! nucleus-nucleus
        if( laddbps ) then
          edel = log( epseudomx / esmin ) / ( nemd_db - 1 )
          do k = 1, nemd_db
            ei = exp( ( k - 1 ) * edel ) * esmin
            call sighi(dble(iap),dble(izp),ei*dble(iap),
     &                 dble(iat),dble(izt),signe,sigel,bmax)
            if( ielas.le.2 ) sigel = 0.d0
            sigt = signe + sigel
            if( sigt.gt.(pseudone+pseudoel) ) then
              pseudone = signe
              pseudoel = sigel
              epeak = ei
            end if
          end do
!$omp critical (dbpseudo)
          ndbcurrent = ndbcurrent + 1
          dbpscs_iza(1,ndbcurrent) = izp
          dbpscs_iza(2,ndbcurrent) = iap
          dbpscs_iza(3,ndbcurrent) = izt
          dbpscs_iza(4,ndbcurrent) = iat
          dbpscs_signe(ndbcurrent) = pseudone
          dbpscs_sigel(ndbcurrent) = pseudoel
          if( epeak.lt.200.d0 ) then
            dbpscs_epeak(ndbcurrent) = epeak
          end if
          if( ndbcurrent.eq.mdbpseud ) then
            laddbps = .false.
          end if
!$omp end critical (dbpseudo)
          if( einpu.lt.epeak ) then
            call sighi(dble(iap),dble(izp),einpu*dble(iap),
     &                 dble(iat),dble(izt),signe,sigel,bmax)
            if( ielas.le.2 ) sigel = 0.d0
            pseudone = signe
            pseudoel = sigel
          end if
        else
          edel = log( einpu / esmin ) / ( nemd - 1 )
          do k = 1, nemd
            ei = exp( ( k - 1 ) * edel ) * esmin
            call sighi(dble(iap),dble(izp),ei*dble(iap),
     &                 dble(iat),dble(izt),signe,sigel,bmax)
            if( ielas.le.2 ) sigel = 0.d0
            sigt = signe + sigel
            if( sigt.gt.(pseudone+pseudoel) ) then
              pseudone = signe
              pseudoel = sigel
            end if
          end do
        end if
      end if

      return
      end subroutine getpseudocs

************************************************************************
      subroutine calc_udm_macro_xs(ein)
************************************************************************
      use udm_Parameter
      use udm_Utility, only : fill_mat_info
      use udm_Manager, only : user_defined_interaction
      implicit real*8(a-h,o-z)
      include 'param.inc'
      include 'ggsparam.inc'
      include 'ggmparam.inc'
      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
!$OMP THREADPRIVATE(/icomon/)
      common /celdn/  denr(kvlmax), denm(kvlmax), denc(kvlmax)
      common /regcm/  icmg(kvlmax)
      common /xtotal/ totttl,tothyd,totdcy,totdel,tottsm
!$OMP THREADPRIVATE(/xtotal/)
      ! -----------------------------
      totudm    =0d0
      totudm_mul=0d0
      ! -----------------------------
      call fill_mat_info(mat)
      do i = 1, udm_int_num
        udm_Kin=ein
        udm_kf_incident=ktyp
        udm_sigt=0d0
        call user_defined_interaction(11,i) ! Calculate: 'udm_sigt'
        totudm_this = udm_sigt * denm(mat) ! Sigma=sigma*N(Atom). udm_sigt=XS[barn], denm=(nAtom*rho/A)[1/cm3]*1e+24
        totudm_(i) = totudm_this           ! Sigma(i)
        totudm     = totudm_this               + totudm
        totudm_mul = totudm_this * udm_bias(i) + totudm_mul
      enddo
      if(abs(ktyp)/=11) totttl = totttl + totudm ! totttl = totttl + 1.0d0 / fpl
      end subroutine calc_udm_macro_xs
