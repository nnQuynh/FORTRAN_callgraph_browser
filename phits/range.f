************************************************************************
*                                                                      *
      subroutine dedxas(eng,dedx,lmat,ityp,ktyp,jtyp,rtyp)
*                                                                      *
*       get dedx of charged particles                                  *
*       modified by K.Niita on 2005/08/25                              *
*       last modified by T.Sato on 2011/09/5                           *
*                                                                      *
*       input                                                          *
*         eng : initial energy (MeV)                                   *
*         lmat : material number                                        *
*         ityp, ktyp, jtyp : particle id                               *
*         rtyp : particle mass (MeV)                                   *
*                                                                      *
*       output                                                         *
*         dedx : dedx of particle (MeV/cm)                             *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'

*-----------------------------------------------------------------------

      common /rngmem/ nrange, krnge, krngn !FURUTA
      common /eparm/  esmax, esmin, emin(20)
      common /ndemax/ dnmax(20)
      common /celdn/  denr(kvlmax), denm(kvlmax), denc(kvlmax)
      common /elect/  uint(3), qs, qo, eint, delc, am, qsex,
     &                nq, ns, n1, noz, mtel
!$OMP THREADPRIVATE(/elect/)
      common /spred/ nspred, nwsprd, nedisp, itstep, ndedx

      common /tcntl/  icntl, inucr


*-----------------------------------------------------------------------
      common /dedxfac/ dedxfd
!$OMP THREADPRIVATE(/dedxfac/)

*-----------------------------------------------------------------------

      integer iegsemi, iegsout
      common /egsemi/ iegsemi, iegsout

*-----------------------------------------------------------------------

            dedx = 0.d0

*-----------------------------------------------------------------------
*        for void nor neutral particles
*-----------------------------------------------------------------------

         if( lmat .eq. 0 .or. icntl .eq. 5 .or.
     &       icntl .eq. 14 .or.
     &       jtyp .eq. 0 ) return

*-----------------------------------------------------------------------

         mat = abs( lmat )

*-----------------------------------------------------------------------
*        for electron
*  Modified by T.Sato on 2011/9/5 to use expected dEdx for electron
*-----------------------------------------------------------------------

         if( ( ityp .eq. 12 .or. ityp .eq. 13 ) .and.
     &          eng .le. dnmax(12) .and. iegsemi .eq. 0 ) then

            if( lmat .gt. 0 ) then

               if(  qs .gt. 0. )
     &         dedx = qsex / denc(mtel) * denc(mat)

             else

               dedx = getdEdxH2O(eng)

             end if

               return

*-----------------------------------------------------------------------

         else if( ( ityp .eq. 12 .or. ityp .eq. 13 ) .and.
     &               eng .le. dnmax(12) .and. iegsemi .ne. 0 ) then

            if( lmat .gt. 0 ) then

               call egs5ededx(eng,mat,ityp,dedx)

             else

               dedx = getdEdxH2O(eng)

             end if

               return

*-----------------------------------------------------------------------
*        for nucleus by SPAR or ATIMA
*-----------------------------------------------------------------------

         else if( ityp .ge. 15 ) then

                  ene = eng
                  ap  = dble( ktyp - ktyp / 1000000 * 1000000 )
                  zp  = dble( ktyp / 1000000 )

            if( ndedx .eq. 0 .or. ndedx .eq. 2 ) then

                  ipty = 1
                  iway = 3
                  call spar(ipty,ap,zp,ene,rtyp,mat,rng,dedx,iway)

            else if( ndedx .eq. 1 .or. ndedx. eq. 3 ) then   ! S.Abe 2016/08/08

                  iway = 5
                  call atima(ap,zp,ene,rtyp,mat,rng,delt,dedx,iway)

            end if

*-----------------------------------------------------------------------
*        for proton by SPAR or ATIMA
*-----------------------------------------------------------------------

         else if( ityp .eq. 1 ) then

            if( ndedx .eq. 0 .or. ndedx .eq. 2 ) then

                  ene = eng
                  ap  = 1.d0
                  zp  = 1.d0

                  ipty = 2
                  iway = 3
                  call spar(ipty,ap,zp,ene,rtyp,mat,rng,dedx,iway)

            else if( ndedx .eq. 1 .or. ndedx .eq. 3 ) then   ! S.Abe 2016/08/08

                  ene = eng
                  ap  = 1.d0
                  zp  = 1.d0

                  iway = 5
                  call atima(ap,zp,ene,rtyp,mat,rng,delt,dedx,iway)

            end if

*-----------------------------------------------------------------------
*        for pion and muon by SPAR
*-----------------------------------------------------------------------

         else if( ityp .eq. 3 .or. ityp .eq. 5 .or.
     &            ityp .eq. 6 .or. ityp .eq. 7 ) then

            if(  ndedx .eq. 0 .or. ndedx .eq. 1 .or. ndedx .eq. 2 ) then

                  ene = eng
                  ap  = 1.d0
                  zp  = 1.d0

                  ipty = 3
                  if( ityp .eq. 6 .or. ityp .eq. 7 ) ipty = 4

                  iway = 3
                  call spar(ipty,ap,zp,ene,rtyp,mat,rng,dedx,iway)

            elseif( ndedx .eq. 3 ) then

                  ene = eng
                  ap = -1.d0 * dble(iabs(ktyp))
                  zp  = 1.d0

                  iway = 5
                  call atima(ap,zp,ene,rtyp,mat,rng,delt,dedx,iway)

            endif

*-----------------------------------------------------------------------

         elseif(  ndedx .eq. 3 .and.
     &            ( ityp .eq. 8 .or. ityp .eq. 10 ) ) then

                  ene = eng
                  ap = -1.d0 * dble(iabs(ktyp))
                  zp  = 1.d0

                  iway = 5
                  call atima(ap,zp,ene,rtyp,mat,rng,delt,dedx,iway)

*-----------------------------------------------------------------------
*        Other particles
*-----------------------------------------------------------------------

         else

                 call rainge(eng*1.05,rng1,mat,ityp,ktyp,jtyp,rtyp)
                 call rainge(eng*0.95,rng2,mat,ityp,ktyp,jtyp,rtyp)

                 if( rng1 - rng2 .gt. 0.d0 ) then

                    dedx = eng * 0.1 / ( rng1 - rng2 )

                 else

                    dedx = 0.0d0

                 end if

         end if

*-----------------------------------------------------------------------

                 dedx = dedx * dedxfd

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine rainge(eng,rng,mat,ityp,ktyp,jtyp,rtyp)
*                                                                      *
*       get range of charged particles                                 *
*       modified by K.Niita on 2005/08/25                              *
*                                                                      *
*       input                                                          *
*         eng : initial energy (MeV)                                   *
*         mat : material number                                        *
*         ityp, ktyp, jtyp : particle id                               *
*         rtyp : particle mass (MeV)                                   *
*                                                                      *
*       output                                                         *
*         rng : range of particle (cm)                                 *
*                                                                      *
************************************************************************
      use MEMBANKMOD !FURUTA
*-----------------------------------------------------------------------

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      common /rngmem/ nrange, krnge, krngn !FURUTA
      common /rnginit/initr                !FURUTA
!$OMP THREADPRIVATE(/rnginit/)
      common /eparm/  esmax, esmin, emin(20)
      common /ndemax/ dnmax(20)
      common /elect/  uint(3), qs, qo, eint, delc, am, qsex,
     &                nq, ns, n1, noz, mtel
!$OMP THREADPRIVATE(/elect/)
      common /spred/ nspred, nwsprd, nedisp, itstep, ndedx

      common /tcntl/  icntl, inucr


*-----------------------------------------------------------------------
      common /dedxfac/ dedxfd
!$OMP THREADPRIVATE(/dedxfac/)

*-----------------------------------------------------------------------
      data rprt / 938.27 /
      data initr /0/

*-----------------------------------------------------------------------

      integer iegsemi, iegsout
      common /egsemi/ iegsemi, iegsout

      data ierrms / 0 /    !S.Abe 2016/05/16
      save ierrms          !S.Abe 2016/05/16

*-----------------------------------------------------------------------
*     initialization
*-----------------------------------------------------------------------

         if( initr .eq. 0 ) then

               initr = initr + 1

               call rangei

         end if

*-----------------------------------------------------------------------
*        for void
*-----------------------------------------------------------------------

         if( mat .le. 0 .or. icntl .eq. 5 .or. icntl .eq. 14 ) then

            rng = 1.d+40
            return

         end if

*-----------------------------------------------------------------------
*        for neutral particles
*-----------------------------------------------------------------------

         if( iegsemi .eq. 0 ) then

            if( jtyp .eq. 0 ) then

               rng = 1.d+40
               return

            end if

         else

            if( jtyp .eq. 0 .and. ityp .ne. 14 ) then

               rng = 1.d+40
               return

            end if

         end if

*-----------------------------------------------------------------------
*        for electron and photon of EGS
*-----------------------------------------------------------------------

         if( ( ityp .eq. 12 .or. ityp .eq. 13 ) .and.
     &         eng .le. dnmax(12) .and. iegsemi .eq. 0 ) then

            if( qo .gt. 0. ) then

                  rng = eng / qo / dedxfd

            else

                  rng = 1.d+40

            end if

                  return

cKN H.Iwase for EGS5
         else if( ( ityp .eq. 12 .or. ityp .eq. 13 .or.
     &              ityp .eq. 14 ) .and.
     &              eng .le. dnmax(12) .and. iegsemi .ne. 0 ) then

                  rng = 1.d+29
                  return

         else if( ( ityp .eq. 12 .or. ityp .eq. 13 .or.
     &              ityp .eq. 14 ) .and.
     &              eng .gt. dnmax(12) ) then

                  rng = ( eng - emin(12) ) * 1.d+40 + 1.d+40
                  return

         end if

*-----------------------------------------------------------------------
*        for nucleus by SPAR or ATIMA
*-----------------------------------------------------------------------

         if( ityp .ge. 15 ) then

                  ene = eng
                  ap  = dble( ktyp - ktyp / 1000000 * 1000000 )
                  zp  = dble( ktyp / 1000000 )

            if( ndedx .eq. 0 .or. ndedx .eq. 2 ) then

                  ipty = 1
                  iway = 1
                  call spar(ipty,ap,zp,ene,rtyp,mat,rng,ecc,iway)

            else if( ndedx .eq. 1 .or. ndedx .eq. 3 ) then   ! S.Abe 2016/08/08

                  iway = 1
                  call atima(ap,zp,ene,rtyp,mat,rng,delt,ecc,iway)

            end if

                  rng = rng / dedxfd
                  return

         end if

*-----------------------------------------------------------------------
*        for proton by ATIMA or SPAR
*-----------------------------------------------------------------------

         if( ityp .eq. 1 .and. ndedx .ne. 0 ) then

            if( ndedx .eq. 1 .or. ndedx .eq. 3 ) then   ! S.Abe 2016/08/08

                  ene = eng
                  ap  = 1.d0
                  zp  = 1.d0

                  iway = 1
                  call atima(ap,zp,ene,rtyp,mat,rng,delt,ecc,iway)

            else if( ndedx .eq. 2 ) then

                  ene = eng
                  ap  = 1.d0
                  zp  = 1.d0

                  ipty = 2
                  iway = 1
                  call spar(ipty,ap,zp,ene,rtyp,mat,rng,ecc,iway)

            end if

                  rng = rng / dedxfd
                  return

         end if

*-----------------------------------------------------------------------
*        for pion and muon by SPAR
*-----------------------------------------------------------------------

         if( ndedx .eq. 2 .and.
     &     ( ityp .eq. 3 .or. ityp .eq. 5 .or.
     &       ityp .eq. 6 .or. ityp .eq. 7 ) ) then

                  ene = eng
                  ap  = 1.d0
                  zp  = 1.d0

                  ipty = 3
                  if( ityp .eq. 6 .or. ityp .eq. 7 ) ipty = 4

                  iway = 1
                  call spar(ipty,ap,zp,ene,rtyp,mat,rng,ecc,iway)

                  rng = rng / dedxfd
                  return

         elseif( ndedx .eq. 3 .and.
     &     ( ityp .eq. 3 .or. ityp .eq. 5 .or.
     &       ityp .eq. 6 .or. ityp .eq. 7 .or.
     &       ityp .eq. 8 .or. ityp .eq. 10 ) ) then

                  ene = eng
                  ap = -1.d0 * dble(iabs(ktyp))
                  zp  = 1.d0

                  iway = 1
                  call atima(ap,zp,ene,rtyp,mat,rng,delt,ecc,iway)

                  rng = rng / dedxfd
                  return

         end if

*-----------------------------------------------------------------------
*        for other particles
*-----------------------------------------------------------------------

         if( ndedx .eq. 3 .and. ityp .eq. 11 ) then

                  ene = eng
                  ap = -1.d0 * dble(iabs(ktyp))
                  zp  = dble(iabs(jtyp))

                  iway = 1
                  call atima(ap,zp,ene,rtyp,mat,rng,delt,ecc,iway)

                  rng = rng / dedxfd
                  return
         endif

*-----------------------------------------------------------------------
*        for zero density matter
*-----------------------------------------------------------------------

         if( rnge(krngn+1,mat) .gt. 1.d+40 ) then

            rng = 1.d+40
            return

         end if

*-----------------------------------------------------------------------
*     Original dE/dx
*-----------------------------------------------------------------------
*        mass and charge corrections
*-----------------------------------------------------------------------

            rmsi = rprt / rtyp * dble( jtyp**2 )
            rmsr = rtyp / rprt / dble( jtyp**2 )

*-----------------------------------------------------------------------
*        get range
*-----------------------------------------------------------------------

            en2 = eng * rmsi

            if( en2 .gt. esmax .and. ierrms .eq. 0 ) then
             ErrCha = ''
             ErrID = 'L:512/R:rainge/F:range.f' !W00_005_001
             call ErrWrite(ErrID,ErrCha)
             write(*,'(a)')
     &        '*** Warning: NMTC-org, energy*rmsi is larger than esmax'
             write(*,'(a,1p3d14.7)')
     &        '    Energy, rmsi, esmax: ',eng,rmsi,esmax
             write(*,'(a)')'    !!! Data is extrapolated !!!'
             ierrms = 1
             i = nrange
             goto 200
            endif

*-----------------------------------------------------------------------

            j = nrange / 2
            i = nrange / 2

  100    continue
               j = j / 2
         if( erng(krnge+i) .lt. en2 ) then
            if( j .eq. 0 ) then
               i = i + 1
               goto 200
            end if
               i = i + j
         else
            if( j .eq. 0 ) goto 200
               i = i - j
         end if
               go to 100
  200    continue

*-----------------------------------------------------------------------

         if( i .gt. 1 ) then

            e   = erng(krnge+i-1)
            rng = rnge(krngn+i-1,mat)

         else

            e   = 0.d+0
            rng = 0.d+0

         end if

*-----------------------------------------------------------------------
*        interpolation of range
*-----------------------------------------------------------------------

            rng = rng
     &          + ( en2 - e ) * ( rnge(krngn+i,mat) - rng )
     &                        / ( erng(krnge+i) - e )

            rng = rng * rmsr

            rng = rng / dedxfd

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine ecol(ee,delt0,e1,rng0,
     &                mat,ityp,ktyp,jtyp,rtyp)
*                                                                      *
*       get energy of charged particle at delt from rng0, e1           *
*       modified by K.Niita on 2013/05/21                              *
*                                                                      *
*       input                                                          *
*         delt0: moving length (cm)                                    *
*         e1   : initial energy (MeV)                                  *
*         rng0 : range of initial energy e1                            *
*                                                                      *
*         mat  : material number                                       *
*         ityp, ktyp, jtyp : particle id                               *
*         rtyp : particle mass (MeV)                                   *
*                                                                      *
*       output                                                         *
*         ee   : energy at delt (MeV)                                  *
*         delt0 : = rng0 if delt0 > rng0                               *
*                                                                      *
************************************************************************
      use MEMBANKMOD !FURUTA
      use MMBANKMOD
      use ion_track_structure, only: e_mean, pts_flt2, no, iflnel, eIlow
     & , lflgTS
      use moddas_material
*-----------------------------------------------------------------------

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'

      include 'err.inc'
      include 'param00.inc'

*-----------------------------------------------------------------------

      common /rngmem/ nrange, krnge, krngn !FURUTA
      common /rnginit/initr                !FURUTA
!$OMP THREADPRIVATE(/rnginit/)
      common /estrag/ kesta, kestz, kestd
      common /eparm/  esmax, esmin, emin(20)
      common /ndemax/ dnmax(20)
      common /elect/  uint(3), qs, qo, eint, delc, am, qsex,
     &                nq, ns, n1, noz, mtel
!$OMP THREADPRIVATE(/elect/)
      common /spred/  nspred, nwsprd, nedisp, itstep, ndedx

      common /tcntl/  icntl, inucr
      common /kmat1g/ kmat(kvlmax)




*-----------------------------------------------------------------------
      common /dedxfac/ dedxfd
!$OMP THREADPRIVATE(/dedxfac/)
      common / etsexe  / dexc_ene
!$OMP THREADPRIVATE(/etsexe/)

      data trec / 0.d0 /
      data erat / 0.d0 /
      data r1   / 0.d0 /
      save trec ! variable to check if this is first attempt
      save erat ! variable to conserve probability
      save r1   ! random number to keep in-event coherence
!$OMP THREADPRIVATE(trec, erat, r1)
*-----------------------------------------------------------------------

      data rprt / 938.27 /

*-----------------------------------------------------------------------

      integer iegsemi, iegsout
      common /egsemi/ iegsemi, iegsout

      data ierrms / 0 /    !S.Abe 2016/05/16
      save ierrms          !S.Abe 2016/05/16

      data iwarnstop /0/  ! T.Sato 2022/11/29
      save iwarnstop      ! T.Sato 2022/11/29
!$OMP THREADPRIVATE(iwarnstop)

*-----------------------------------------------------------------------
      common /reslet/ irlet
      common /delesg/ rdels(kvlmax), mndel, ndels(kvlmax)
      common /delreg/ delm(kvlmax), kdelt
      common /regdc/  idrg(kvlmax), idgr(kvmmax)
      common /tlgeom/ iblz1,iblz2
!$OMP THREADPRIVATE(/tlgeom/)
      parameter(rmase=0.511d0)  ! electron mass
      common /  tscreg   / ntscell(kvlmax)

      common /adjoint/ iadjnt ! T.Sato 2024/01/07

*-----------------------------------------------------------------------
               rng1 = rng0  * dedxfd
               delt = delt0 * dedxfd

*-----------------------------------------------------------------------

               ee = e1

*-----------------------------------------------------------------------
*        for void  and  for neutral particles
*-----------------------------------------------------------------------

         if( mat .le. 0 .or. icntl .eq. 5 .or. icntl .eq. 14 ) return
         if( jtyp .eq. 0 ) return

*-----------------------------------------------------------------------
*     initialization
*-----------------------------------------------------------------------

         if( initr .eq. 0 ) then

               initr = initr + 1

               call rangei

         end if

*-----------------------------------------------------------------------
*        range of the particle and stopped particle
*-----------------------------------------------------------------------

               rr = rng1 - delt

         if( rr .le. 0.0d0 ) then

               delt0 = rng0
               delt  = rng0
               dexc_ene = ee - emin(ityp) * ibryf(ityp,ktyp) ! in case of ITSART, energy deposition at end of track

               ee    = 0.0d0

               goto 900

         end if

*-----------------------------------------------------------------------
*        for electron
*-----------------------------------------------------------------------

         if( ( ityp .eq. 12 .or. ityp .eq. 13 ) .and.
     &         e1 .le. dnmax(12) .and. iegsemi .eq. 0 ) then

               ee = rr * qo
               return

         else if( ( ityp .eq. 12 .or. ityp .eq. 13 ) .and.
     &              e1 .le. dnmax(12) .and. iegsemi .ne. 0 ) then

               call egs5edxde(ee,delt0,e1,mat,ityp)

               return

         else if( ityp .eq. 14 .and.
     &              e1 .le. dnmax(14) .and. iegsemi .ne. 0 ) then

               call egs5edxde(ee,delt0,e1,mat,ityp)

               return

         else if( ( ityp .eq. 12 .or. ityp .eq. 13 ) .and.
     &              e1 .gt. dnmax(12) ) then

               ee = e1
               return

         end if

*-----------------------------------------------------------------------
*        for nucleus by SPAR or ATIMA
*-----------------------------------------------------------------------

         if( ityp .ge. 15 ) then

               ene = e1
               rng = rr
               ap  = dble( ktyp - ktyp / 1000000 * 1000000 )
               zp  = dble( ktyp / 1000000 )

            if( ndedx .eq. 0 .or. ndedx .eq. 2 ) then

                  ipty = 1
                  iway = 2
                  call spar(ipty,ap,zp,ene,rtyp,mat,rng,ee,iway)

            else if( ndedx .eq. 1 .or. ndedx .eq. 3 ) then   ! S.Abe 2016/08/08

               if( nedisp .ne. 10 ) then

                  iway = 2
                  call atima(ap,zp,ene,rtyp,mat,rng1,delt,ee,iway)

               else

                  iway = 3
                  call atima(ap,zp,ene,rtyp,mat,rng1,delt,ee,iway)

                  goto 900

               end if

            end if

                  goto 900

         end if

*-----------------------------------------------------------------------
*        for proton by ATIMA or SPAR
*-----------------------------------------------------------------------

         if( ityp .eq. 1 .and. ndedx .ne. 0 ) then

            if( ndedx .eq. 1 .or. ndedx .eq. 3 ) then   ! S.Abe 2016/08/08

                  ene = e1
                  ap  = 1.d0
                  zp  = 1.d0

               if( nedisp .ne. 10 ) then

                  iway = 2
                  call atima(ap,zp,ene,rtyp,mat,rng1,delt,ee,iway)

                  goto 900

               else

                  iway = 3
                  call atima(ap,zp,ene,rtyp,mat,rng1,delt,ee,iway)

                  goto 900

               end if

            else if( ndedx .eq. 2 ) then

                  ene = e1
                  rng = rr
                  ap  = 1.d0
                  zp  = 1.d0

                  ipty = 2
                  iway = 2
                  call spar(ipty,ap,zp,ene,rtyp,mat,rng,ee,iway)

                  goto 900

            end if

         end if

*-----------------------------------------------------------------------
*        for pion and muon by SPAR
*-----------------------------------------------------------------------

         if( ndedx .eq. 2 .and.
     &     ( ityp .eq. 3 .or. ityp .eq. 5 .or.
     &       ityp .eq. 6 .or. ityp .eq. 7 ) ) then

                  ene = e1
                  rng = rr
                  ap  = 1.d0
                  zp  = 1.d0

                  ipty = 3
                  if( ityp .eq. 6 .or. ityp .eq. 7 ) ipty = 4

                  iway = 2
                  call spar(ipty,ap,zp,ene,rtyp,mat,rng,ee,iway)

                  goto 900

         elseif( ndedx .eq. 3 .and.
     &     ( ityp .eq. 3 .or. ityp .eq. 5 .or.
     &       ityp .eq. 6 .or. ityp .eq. 7 .or.
     &       ityp .eq. 8 .or. ityp .eq. 10 ) ) then

                  ene = e1
                  ap = -1.d0 * dble(iabs(ktyp))
                  zp  = 1.d0

                  if( nedisp .ne. 10 ) then

                     iway = 2
                     call atima(ap,zp,ene,rtyp,mat,rng1,delt,ee,iway)

                     goto 900

                  else

                     iway = 3
                     call atima(ap,zp,ene,rtyp,mat,rng1,delt,ee,iway)

                     goto 900

                  end if

         end if

*-----------------------------------------------------------------------
*        for other particles
*-----------------------------------------------------------------------

         if( ndedx .eq. 3 .and. ityp .eq. 11 ) then

                  ene = e1
                  ap = -1.d0 * dble(iabs(ktyp))
                  zp  = dble(iabs(jtyp))

                  if( nedisp .ne. 10 ) then

                     iway = 2
                     call atima(ap,zp,ene,rtyp,mat,rng1,delt,ee,iway)

                     goto 900

                  else

                     iway = 3
                     call atima(ap,zp,ene,rtyp,mat,rng1,delt,ee,iway)

                     goto 900

                  end if

         endif

*-----------------------------------------------------------------------
*        for zero density matter
*-----------------------------------------------------------------------

         if( rnge(krngn+1,mat) .gt. 1.d+40 ) then

            ee = e1
            return

         end if

*-----------------------------------------------------------------------
*     Original dE/dx
*-----------------------------------------------------------------------
*        mass and charge corrections
*-----------------------------------------------------------------------

            rmsi = rprt / rtyp * dble( jtyp**2 )
            rmsr = rtyp / rprt / dble( jtyp**2 )

*-----------------------------------------------------------------------
*        get energy
*-----------------------------------------------------------------------

            rngp = rr * rmsi

            if( rngp.gt.rnge(krngn+nrange,mat) .and. ierrms.eq.0 ) then
             ErrCha = ''
             ErrID = 'L:939/R:ecol/F:range.f' !W00_005_002
             call ErrWrite(ErrID,ErrCha)
             write(*,'(a)')
     &        '*** Warning: NMTC-org, energy*rmsi is larger than esmax'
             write(*,'(a,1p3d14.7)')
     &        '    Energy, rmsi, esmax: ',e1,rmsi,esmax
             write(*,'(a)')'    !!! Data is extrapolated !!!'
             ierrms = 1
             i = nrange
             goto 200
            endif

*-----------------------------------------------------------------------

            j = nrange / 2
            i = nrange / 2

  100    continue
               j = j / 2
         if( rnge(krngn+i,mat) .le. rngp ) then
            if( j .eq. 0 ) then
               i = i + 1
               goto 200
            end if
               i = i + j
         else
            if( j .eq. 0 ) goto 200
               i = i - j
         end if
               go to 100
  200    continue

*-----------------------------------------------------------------------

         if( i .gt. 1 ) then

            ee = erng(krnge+i-1)
            r  = rnge(krngn+i-1,mat)

         else

            ee = 0.d+0
            r  = 0.d+0

         end if

*-----------------------------------------------------------------------
*        interpolation of energy
*-----------------------------------------------------------------------

            ee = ee
     &         + ( rngp - r ) * ( erng(krnge+i) - ee )
     &                        / ( rnge(krngn+i,mat) - r )

            ee = ee * rmsr

*-----------------------------------------------------------------------
*     without energy straggling
*-----------------------------------------------------------------------

  900 continue

         if( (irlet .eq. 1 .or. lflgTS .eq. 2) .and. lflgTS .ne. 1) then ! if lflgTS .eq. 1, this is called from tally. Energy should be recovered later in transport phase.

          icl = idgr(iblz1)

            lflg1 = 0
            if(rtyp .gt. 1.d0 .and. mndel .gt. 0 .and.
     &          delm(icl) .gt. 1.0d-3) lflg1 = 1   ! Delta-ray on, changed by T.Sato 2022/11/24

            if( jtyp .ne. 0 .and. lflgTS .eq. 2 .and. ee .eq. 0.d0) then   ! Stopped in ITSART

               ee = ibryf(ityp,ktyp) * emin(ityp)

            elseif( jtyp  .ne. 0 .and.
     &         (lflg1  .eq. 1 .or.
     &          lflgTS .eq. 2)
     &        ) then   ! Kurbuc not supported.

               rmass = rtyp
               chag  = dabs(dble( jtyp ))
               ein   = e1

               if(lflg1 .eq. 1) then
                  edens = edns(kdelt+mat) ! T.Sato 2022/11/21, move to here to avoid access violation

                  emind = delm(icl)
                  emaxd = 4.d0 * rmase / rmass * ein
     &                   * ( 1.d0 + ein / 2.d0 / rmass )
     &                   / ( ( 1.d0 + rmase /rmass )**2
     &                       + 2.d0 * rmase * ein / rmass**2 )
                  edmean = emind * emaxd / (emaxd - emind)
     &                    * dlog(emaxd/emind)                   ! derived from integration of edel in subroutine delprd

                  totdel = delmfp(rmass,chag,edens,ein,emind)   ! delta ray flight length (1/cm)

               elseif(lflgTS .eq. 2) then
                  edmean = e_mean(e1*1.d6) * 1.d-6
                  if(edmean .gt. 0.d0) then
                    call pts_flt2(e1*1.d6, sig_macro)
                    totdel = sig_macro
                  else
                    totdel = 0.d0
                  endif
               endif

               if( totdel .gt. 0.d0 ) then
                  edgain = edmean * totdel * delt0
               else
                  edgain = 0.d0
               endif

!               if(irlet.eq.1) then ! used to use irlet as a flag for delta-ray production, but this condition also meets the ITSART mode
               if(lflg1.eq.1) then ! delta-ray mode, T.Sato 2023/12/30
                edelt = e1 - ee
               else
!  T.Sato 2022/11/22, ee is sometimes stragne because of unstability of e_out in atima
                iway=5
                call atima(ap,zp,ene,rtyp,mat,rngtmp,delt0,dedxtmp,iway)
                edelt = dedxtmp*delt0
                ee = e1-edelt
               endif

               if( edgain .gt. edelt .and. lflg1 .eq. 1) edgain = edelt


       if(edgain.gt.edelt.and.iwarnstop.eq.0) then
        write(*,'("Warning: dE/dx(ITSART)",es12.4," > dE/dx(ATIMA)",
     &  es12.4," at E(MeV/n) = ",es12.4)') edgain/delt0,edelt/delt0,
     &  ene/ap
        iwarnstop=1
        edgain=edelt
       endif


               ee = ee + edgain

            endif
            lflgTS = 0

         endif

         if( e1 .le. 0.d0 ) goto 999
         if( nedisp .eq. 0 ) goto 999
         if( nedisp .eq. 10 ) goto 999   ! S.Abe 2017/06/05

*-----------------------------------------------------------------------
*     energy straggling
*-----------------------------------------------------------------------

               step  = delt
               zp    = dble( jtyp )
               xmass = rtyp
               tin   = e1

               seka = avam(kesta+mat)
               sekz = avaz(kestz+mat)
               sekm = avad(kestd+mat)

               call glando(step,zp,xmass,tin,sekz,seka,sekm,de,iflag)

               ee = max( 0.0d0, ee - de )

*-----------------------------------------------------------------------
  999 if(iadjnt.eq.2) ee=2.0d0*e1-ee ! T.Sato 2024/01/07 adjoint mode

      return
      end


************************************************************************
*                                                                      *
      subroutine rangei
*                                                                      *
*       initialization of range                                        *
*       modified by K.Niita on 23/12/1999                              *
*                                                                      *
*                                                                      *
************************************************************************
      use MEMBANKMOD !FURUTA
      use NGSDATAMOD, only : weitn
      use moddas_material
*-----------------------------------------------------------------------

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'

*-----------------------------------------------------------------------

      common /eparm/  esmax, esmin, emin(20)
      common /cparm/  maxbch,maxcas
      common /rngmem/ nrange, krnge, krngn !FURUTA
      common /spred/  nspred, nwsprd, nedisp, itstep, ndedx
      common /estrag/ kesta, kestz, kestd

      common /kmat1a/ mxmat, mxmat0, mxnel
      common /kmat1g/ kmat(kvlmax)




*-----------------------------------------------------------------------
*        constants
*-----------------------------------------------------------------------

            scale  = 1.0d0

*-----------------------------------------------------------------------

            emx  = esmax * scale
            umax = log( emx / esmin )
            fgrp = nrange

            du   = umax / fgrp
            u1   = umax - du

*-----------------------------------------------------------------------

            erng(krnge+1) = emx * exp(-u1)

         do i = 2, nrange

            u1 = u1 - du

            erng(krnge+i) = emx * exp(-u1)

         end do

*-----------------------------------------------------------------------

      do n = 1, mxmat

               u1   = umax - du
               sum1 = denh_das(kmat0+n)
               sum2 = denh_das(kmat0+n) * (-10.956)
               mmn  = nint( dnel_das(kmat0+n) )

            do i = 1, mmn

               t = den_das(kmat(n)+i)*zz_das(kmat(n)+i)

               sum1 = sum1 + t
               sum2 = sum2 + t*log( zfoi(zz_das(kmat(n)+i))*1.d-6 )

            end do

         if( sum1 .gt. 0.0 ) then

               call dxde(erng(krnge+1),sum1,sum2,f1)
               call dxde(esmin,sum1,sum2,fmin)

               rnge(krngn+1,n) = ( f1 * erng(krnge+1) + fmin * esmin )
     &                 * du / 2.0

               fi = f1

            do i = 2, nrange

               u1  = u1 - du
               fii = fi

               call dxde(erng(krnge+i),sum1,sum2,fi)

               rnge(krngn+i,n) = rnge(krngn+i-1,n)
     &                 + ( fi  * erng(krnge+i)
     &                   + fii * erng(krnge+i-1) ) * du / 2.0

            end do

         else

            do i = 1, nrange

               rnge(krngn+i,n) = 1.0d+61 * i

            end do

         end if

      end do

*-----------------------------------------------------------------------
*     energy straggling,  average z, a and density
*-----------------------------------------------------------------------

      if( nedisp .ne. 0 ) then

         do m = 1, mxmat

               sekm = 0.0d0
               seka = 0.0d0
               sekz = 0.0d0
               sekn = 0.0d0

               lem   = nint( dnel_das(kmat0+m) )
               hydro = denh_das(kmat0+m)

            if( hydro .gt. 0.0d0 ) then

               zin = 1.0
               ain = 1.0
               itz = nint( zin )
               itn = nint( ain - zin )

               rho = hydro * weitn(itz,itn)

               sekm = sekm + rho
               seka = seka + ain * hydro
               sekz = sekz + zin * hydro
               sekn = sekn + hydro

            end if

         do i = 1, lem

               zin = zz_das(kmat(m)+i)
               ain = a_das(kmat(m)+i)
               itz = nint( zin )
               itn = nint( ain - zin )

               rho = den_das(kmat(m)+i)*weitn(itz,itn)

               sekm = sekm + rho
               seka = seka + ain*den_das(kmat(m)+i)
               sekz = sekz + zin*den_das(kmat(m)+i)
               sekn = sekn + den_das(kmat(m)+i)

         end do

               seka = seka / sekn
               sekz = sekz / sekn

               avam(kesta+m) = seka
               avaz(kestz+m) = sekz
               avad(kestd+m) = sekm

         end do

      end if

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      subroutine dxde(ee,sum1,sum2,x1)
*                                                                      *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      a2  = ( 1.d0 + ee / 938.232 ) * ( 1.d0 + ee / 938.232 )
      a3  =  a2 - 1.d0
      a4  =  a3 / a2
      a6  = log(a3) - a4 + 0.0217615
      x1  = a4 / 5.0985d-1 / ( a6 * sum1 - sum2 )

      return
      end


************************************************************************
*                                                                      *
      function zfoi(z)
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      dimension fli(4)

      data fli/17.5,44.,34.5,59./

*-----------------------------------------------------------------------

      zfoi=9.*z
      if(z.gt.30.)return
      if(z.lt.5.)go to 40
   30 zfoi = z*(12.8-(3.8*(z-5.)/25.))
      return
   40 j=z
      zfoi = fli(j)
      return
      end


************************************************************************
*                                                                      *
      subroutine atima(ap,cp,ene,rms,mat,rng,delt,ecc,iway)
*                                                                      *
*       control routine of ATIMA                                       *
*       last modified by K.Niita on 2006/03/29                         *
*                                                                      *
*       ap   : mass number of projectile                               *
*       cp   : charge of projectile                                    *
*       ene  : initial energy of projectile (MeV)                      *
*       rms  : mass of projectile  (MeV)                               *
*       mat  : material number                                         *
*       rng  : range of initial energy ene                             *
*       delt : distance                                                *
*       ecc  : final energy at rng or distance, or angle straggling    *
*                                                                      *
*       iway : 1=> range, 2=>energy,                                   *
*              3=> energy straggling,                                  *
*              4=> angle straggling                                    *
*              5=> dE/dx at the initial ene, ecc is the output         *
*                                                                      *
************************************************************************
      use MEMBANKMOD,only:rhoi !FURUTA
      use MEMBANKMOD,only:kspc,ispc,mdbatima,ndbatima,dbcutoff !Wada20160118
      use MEMBANKMOD,only:dbfspace  ,dbispace  ,dbfispace      !Wada20160118
      use MEMBANKMOD,only:dbfspacene,dbispacene,dbfispacene    !Wada20160118
      use MEMBANKMOD,only:dbfspacens,dbispacens,dbfispacens    !Wada20160118
      use NGSDATAMOD, only : weitn
      use mod_ompparallel
      use moddas_material

      implicit real*8 (a-h,o-z)

      include 'param.inc'
      include 'param-physcnst.inc' !FURUTA20200219
      include 'atimasys.inc'
      include 'atimadim.inc'
      include 'err.inc'

      common /icomon/ no  ! T.Sato 2019/02/18
!$OMP THREADPRIVATE(/icomon/)
      common /fixcharge/ifixchg     !T.Sato 2019/02/17

*-----------------------------------------------------------------------

      dimension fspace(kspc), ispace(ispc)
      save fspace, ispace           !Wada20160118
!$OMP THREADPRIVATE(fspace, ispace) !Wada20160118
      external ffrange, ffrstr, ffastr, fftof

*-----------------------------------------------------------------------

      common /kmat1a/ mxmat, mxmat0, mxnel
      common /kmat1g/ kmat(kvlmax)

      common /kmat1h/ kmatg(kvlmax)

      common /atima01/ katima
      common /spred/  nspred, nwsprd, nedisp, itstep, ndedx
      common /eparm/  esmax, esmin, emin(20)
      common /paraj/  mstz(300), parz(300)
      common /ionpot/ rh2o, ih2o
!$OMP THREADPRIVATE(/ionpot/)


*-----------------------------------------------------------------------

      data ifirst / 0 /
      data matsav / 0 /
      data jzp / 0 /
      data jap / 0 /
      save ifirst,matsav,jzp,jap !FURUTA
!$OMP THREADPRIVATE(ifirst,matsav,jzp,jap)
      data mdbexceed / 0 / !FURUTA20160128
      save mdbexceed       !FURUTA20160128

      data ierrms,ierrms3,ierrms4 / 0, 0, 0 /  !S.Abe 2016/05/16
      save ierrms,ierrms3,ierrms4              !S.Abe 2016/05/16

      common /paraspar/ mat_spar !S.Abe 2017/09/13

      common /iaoru/  iflagAorU !FURUTA20200219 E/A <-> E/u
      integer ipompl
      ipompl = ipomp+1
*-----------------------------------------------------------------------

      mat_spar = mat

      if( mxnel + 1 .gt. maxnuc ) then

         ErrCha = ''
         ErrID = 'L:1430/R:atima/F:range.f' !E00_013_001
         call ErrWrite(ErrID,ErrCha)

         write(*,*) '*** Error: maxnuc is smaller than mxnel'
         write(*,'(''mxnel ='',i3,''maxnuc ='',i3)') mxnel, maxnuc
         call parastop( 540 )

      end if

*-----------------------------------------------------------------------
*     esmin and esmax
*-----------------------------------------------------------------------

      if( ap .ge. 0.d0 ) then
       esminap = esmin*ap
       if(iflagAorU.eq.0)then
        esminu = esmin*ap/rms*physc(5)*1000.0d0 !FURUTA20200219 E/A > E/u
       else
        esminu = esmin
       endif
      else
       esminap = esmin
       esminu  = esmin !FURUTA20200226
      endif
      if( ene .le. esminap ) then
         ecc = 1.d-11 !FURUTA20160422
         rng = 1.d-10 ! make it same as SPAR
         return
      end if

      if( iway .eq. 2 .and. rng .lt. 1.0d-10 ) then
         ecc = 1.d-11
         return
      end if


cFURUTA20160422---------------------------------------------------------
      if( ap .ge. 0.d0 ) then
       esmaxap = esmax*ap
       if(iflagAorU.eq.0)then
        esmaxu = esmax*ap/rms*physc(5)*1000.0d0 !FURUTA20200219 E/A > E/u
       else
        esmaxu = esmax
       endif
       eneap = ene/dble(ap)
      else
       esmaxap = esmax
       eneap = ene
       esmaxu = esmax !FURUTA20200226
      endif
      if( ene .gt. esmaxap*1.001 .and. ierrms.eq.0 ) then !FURUTA20160422, T.Sato 2024/07/03 to avoid significant digit
       ErrCha = ''
       ErrID = 'L:1482/R:atima/F:range.f' !W00_006_001
       call ErrWrite(ErrID,ErrCha)

       write(*,'(a)')'*** Warning: ATIMA, energy is larger than esmax'
       write(*,'(a,1p2d14.7)')
     &      '    Energy/A, esmax: ',eneap,esmax

       write(*,'(a)')'    !!! Data is extrapolated !!!'
       ierrms = 1
      endif
c-----------------------------------------------------------------------

*-----------------------------------------------------------------------
*     some constants
*-----------------------------------------------------------------------

        iprint = 1
        ndimf  = kspc
        ndimi  = ispc

*-----------------------------------------------------------------------
*     store density of material (g/cm^3)
*-----------------------------------------------------------------------

      if( ifirst .eq. 0 ) then

         do m = 1, mxmat

               sekm = 0.0d0

               lem   = nint(dnel_das(kmat0+m) )
               hydro = denh_das(kmat0+m)

            if( hydro .gt. 0.0d0 ) then

               zin = 1.0
               ain = 1.0
               itz = nint( zin )
               itn = nint( ain - zin )

               rhom = hydro * weitn(itz,itn)

               sekm = sekm + rhom

            end if

            do i = 1, lem

               zin = zz_das(kmat(m)+i)
               ain = a_das(kmat(m)+i)
               itz = nint( zin )
               itn = nint( ain - zin )

               rhom = den_das(kmat(m)+i)*weitn(itz,itn)

               sekm = sekm + rhom

            end do

               rhoi(katima+m) = sekm

         end do

            ifirst = 1

      end if

*-----------------------------------------------------------------------
*     store new material information and coefficients
*-----------------------------------------------------------------------

         if( mat .ne. matsav ) then

               k = 0

               lem = nint( dnel_das(kmat0+mat) )
               hydro = denh_das(kmat0+mat)

            if( hydro .gt. 0.0d0 ) then

               k = k + 1

               zt(k) = 1
               mt(k) = 1.d0

               anuc(k) = hydro

            end if

            do i = 1, lem

               k = k + 1

               zt(k) = nint( zz_das(kmat(mat)+i) )
               mt(k) = a_das(kmat(mat)+i)
               anuc(k) = den_das(kmat(mat)+i)

            end do


              gas = nint( das_kmatg(kmatg(mat)+2) )
              rho = rhoi(katima+mat)

              fntp = 0
              nnuc = k

*-----------------------------------------------------------------------
*           modified ionizaton potential for water
*           iH2O (parz(174)) <0: normal(69eV), >0 potential (eV)
*-----------------------------------------------------------------------

                     ih2o = 0

cFURUTA20190717 ! The function is activated also for natural hydrogen
            anuc_H=0.0d0
            anuc_O=0.0d0
            iflag=0
            do i=1,nnuc
             if(zt(i).eq.1)then
              anuc_H=anuc_H+anuc(i)
             elseif(zt(i).eq.8)then
              anuc_O=anuc_O+anuc(i)
             else
              iflag=1 ! Containing other elements
              exit
             endif
            enddo
            if(iflag.eq.0.and.anuc_H.gt.0.0d0.and.anuc_O.gt.0.0d0)then
             if(anuc_H/anuc_O.gt.1.9d0.and.anuc_H/anuc_O.lt.2.1d0)then
              ! H/O ratio is about 2.0 criteria changed to 5% on 20200324
              ih2o = 1
              rh2o = 75.d0 / 69.d0
              if( parz(174) .gt. 0.d0 ) then
               rh2o = parz(174) / 69.d0
              end if
             end if
            end if

*-----------------------------------------------------------------------
*        pot(j) : ionization potential
*-----------------------------------------------------------------------

            do j = 1, nnuc

               call scoef(zt(j),dummy1,dummy2,dummy3,pot(j))

            end do

*-----------------------------------------------------------------------

         end if

*-----------------------------------------------------------------------
*        calculate tables
*-----------------------------------------------------------------------

            iap = nint( ap )
            izp = nint( cp )

        if(ifixchg.ge.1.and.ctyp.ne.0) izp = abs(ctyp) ! T.Sato 2019/02/18, use actual charge.  2022/7/21 Ogawa replace nstz with ctyp

            zp = izp

         if( ap .gt. 0.d0 ) then
            ep = ene/ap
cFURUTA20200219 MeV/n -> MeV/u
          if(iflagAorU.eq.0)then
            mp = rms/physc(5)*0.001d0
          else
            mp = ap
          endif
         elseif( ap .le. 0.d0 ) then   ! pion, muon, etc...
            mp = rms / physc(5) *0.001d0
            ep = ene
         endif

C for dedxfile read mode
         dedxmat = mat



         if( mat .ne. matsav .or. iap .ne. jap .or. izp .ne. jzp ) then
!$omp critical (dbspace)
               do i = 1, mdbatima
                 if( dbfispace(1,i) .eq. 0 ) exit
                 if( dbfispace(1,i) .eq. iap .and.
     &               dbfispace(2,i) .eq. izp .and.
     &               dbfispace(3,i) .eq. mat       ) then
                   ispace(1:8        ) = dbispace(1:8        ,i)
                   fspace(1:ispace(7)) = dbfspace(1:ispace(7),i)
                   exit
                 end if
               end do
               if(i.eq.mdbatima+1.and.mdbexceed.eq.0)then
                write(6,*)'**** WARNING: MDBATIMA overflowed ',
     &               '********************************'
                write(6,*)' calculation may be speeded up',
     &               ' by increasing MDBATIMA or DBCUTOFF'
                write(6,*)'**********************************',
     &               '********************************'
                mdbexceed=1
               endif
               if( i .eq. mdbatima+1 )then
                iflag=1
               else
                if(dbfispace(1,i) .eq. 0)then
                 iflag=1
                else
                 iflag=0
                endif
               endif
               if( iflag.eq.1 )then
c----------------------------------------------------------------------
                 a = esminu      !FURUTA20200219 E/A -> E/u
                 b = esmaxu      !FURUTA20200219 E/A -> E/u
                 c = a
                 tol = 1.d-6



                 call integr(ffrange,a,b,c,tol,iprint,fspace,ispace,
     &                ndimf,ndimi,iflag)
                 if( i .ne. mdbatima+1 .and. ep .ge. dbcutoff ) then
                   dbfispace(1,i) = iap
                   dbfispace(2,i) = izp
                   dbfispace(3,i) = mat
                   dbispace(1:8        ,i) = ispace(1:8        )
                   dbfspace(1:ispace(7),i) = fspace(1:ispace(7))
                   ndbatima(1)=i
                 end if
               endif
!$omp end critical (dbspace)
c----------------------------------------------------------------------
               call bscoef(1,fspace,ispace,t1,bcoef1,nb1,kpm1) !Wada20160118

            if( nedisp .eq. 10 ) then
!$omp critical (dbspacene)
               do i = 1, mdbatima
                 if( dbfispacene(1,i) .eq. 0 ) exit
                 if( dbfispacene(1,i) .eq. iap .and.
     &               dbfispacene(2,i) .eq. izp .and.
     &               dbfispacene(3,i) .eq. mat       ) then
                   ispace(1:8        ) = dbispacene(1:8        ,i)
                   fspace(1:ispace(7)) = dbfspacene(1:ispace(7),i)
                   exit
                 end if
               end do
               if(i.eq.mdbatima+1.and.mdbexceed.eq.0)then
                write(6,*)'**** WARNING: MDBATIMA overflowed ',
     &               '********************************'
                write(6,*)' calculation may be speeded up',
     &               ' by increasing MDBATIMA or DBCUTOFF'
                write(6,*)'**********************************',
     &               '********************************'
                mdbexceed=1
               endif
               if( i .eq. mdbatima+1 )then
                iflag=1
               else
                if(dbfispacene(1,i) .eq. 0)then
                 iflag=1
                else
                 iflag=0
                endif
               endif
               if( iflag.eq.1 )then
c----------------------------------------------------------------------
                 a = esminu      !FURUTA20200219 E/A -> E/u
                 b = esmaxu      !FURUTA20200219 E/A -> E/u
                 c = a
                 tol = 1.d-6


                 call integr(ffrstr,a,b,c,tol,iprint,fspace,ispace,
     &                ndimf,ndimi,iflag)
                 if( i .ne. mdbatima+1 .and. ep .ge. dbcutoff ) then
                   dbfispacene(1,i) = iap
                   dbfispacene(2,i) = izp
                   dbfispacene(3,i) = mat
                   dbispacene(1:8        ,i) = ispace(1:8        )
                   dbfspacene(1:ispace(7),i) = fspace(1:ispace(7))
                   ndbatima(2)=i
                 end if
               endif
!$omp end critical (dbspacene)
c----------------------------------------------------------------------
 20            call bscoef(1,fspace,ispace,t2,bcoef2,nb2,kpm2) !Wada20160118

            end if

            if( nspred .eq. 10 ) then
!$omp critical (dbspacens)
               do i = 1, mdbatima
                 if( dbfispacens(1,i) .eq. 0 ) exit
                 if( dbfispacens(1,i) .eq. iap .and.
     &               dbfispacens(2,i) .eq. izp .and.
     &               dbfispacens(3,i) .eq. mat       ) then
                   ispace(1:8        ) = dbispacens(1:8        ,i)
                   fspace(1:ispace(7)) = dbfspacens(1:ispace(7),i)
                   exit
                 end if
               end do
               if(i.eq.mdbatima+1.and.mdbexceed.eq.0)then
                write(6,*)'**** WARNING: MDBATIMA overflowed ',
     &               '********************************'
                write(6,*)' calculation may be speeded up',
     &               ' by increasing MDBATIMA or DBCUTOFF'
                write(6,*)'**********************************',
     &               '********************************'
                mdbexceed=1
               endif
               if( i .eq. mdbatima+1 )then
                iflag=1
               else
                if(dbfispacens(1,i) .eq. 0)then
                 iflag=1
                else
                 iflag=0
                endif
               endif
               if( iflag.eq.1 )then
c----------------------------------------------------------------------
                 a = esminu      !FURUTA20200219 E/A -> E/u
                 b = esmaxu      !FURUTA20200219 E/A -> E/u
                 c = b
                 tol = 1.d-8


                 call integr(ffastr,a,b,c,tol,iprint,fspace,ispace,
     &                ndimf,ndimi,iflag)
                 if( i .ne. mdbatima+1 .and. ep .ge. dbcutoff ) then
                   dbfispacens(1,i) = iap
                   dbfispacens(2,i) = izp
                   dbfispacens(3,i) = mat
                   dbispacens(1:8        ,i) = ispace(1:8        )
                   dbfspacens(1:ispace(7),i) = fspace(1:ispace(7))
                   ndbatima(3)=i
                 end if
               endif
!$omp end critical (dbspacens)
c----------------------------------------------------------------------
 30            call bscoef(1,fspace,ispace,t3,bcoef3,nb3,kpm3) !Wada20160118

            end if

         end if

            matsav = mat
            jap = iap
            jzp = izp

*-----------------------------------------------------------------------
*        range
*-----------------------------------------------------------------------

         if( iway .eq. 1 ) then

            ein = ene / mp   ! S.Abe 2016/08/08

            if( ein .le. esmaxu ) then !FURUTA20200219 E/A -> E/u
             range = bvalue(t1,bcoef1,nb1,kpm1,ein,0)  ! range (mg/cm2)
            else
             ene1 = esmaxu * 0.99d0 !FURUTA20200219 E/A -> E/u
             ene2 = esmaxu          !FURUTA20200219 E/A -> E/u
             range1 = bvalue(t1,bcoef1,nb1,kpm1,ene1,0)  ! range (mg/cm2)
             range2 = bvalue(t1,bcoef1,nb1,kpm1,ene2,0)  ! range (mg/cm2)

             range = range1 + (ein-ene1) / (ene2-ene1) * (range2-range1)
            endif

            rng   = range / rho / 1000.d0             ! (cm)

*-----------------------------------------------------------------------

         else if( iway .eq. 2 ) then

            thick = delt * rho * 1000.d0

            ein = ene / mp   ! S.Abe 2016/08/08

            if( ein .le. esmaxu ) then !FURUTA20200219 E/A -> E/u
             range = bvalue(t1,bcoef1,nb1,kpm1,ein,0)  ! range (mg/cm2)
             ecc = e_out(t1,bcoef1,nb1,kpm1,range,ein,thick) ! (MeV/u)
     &           * mp   ! S.Abe 2016/08/08
            else
             ene1 = esmaxu * 0.99d0 !FURUTA20200219 E/A -> E/u
             ene2 = esmaxu          !FURUTA20200219 E/A -> E/u
             range1 = bvalue(t1,bcoef1,nb1,kpm1,ene1,0)  ! range (mg/cm2)
             range2 = bvalue(t1,bcoef1,nb1,kpm1,ene2,0)  ! range (mg/cm2)
             ecc1 = e_out(t1,bcoef1,nb1,kpm1,range1,ene1,thick) ! (MeV/u)
     &            * mp   ! S.Abe 2016/08/08
             ecc2 = e_out(t1,bcoef1,nb1,kpm1,range2,ene2,thick) ! (MeV/u)
     &            * mp   ! S.Abe 2016/08/08

             range = range1 + (ein-ene1) / (ene2-ene1) * (range2-range1)
             ecc = ecc1 + (range-range1) / (range2-range1) * (ecc2-ecc1)
            endif

*-----------------------------------------------------------------------

         else if( iway .eq. 3 ) then

            thick = delt * rho * 1000.d0

            ein = ene / mp   ! S.Abe 2016/08/08
            if( ein .le. esmaxu ) then !FURUTA20200219 E/A -> E/u
             range = bvalue(t1,bcoef1,nb1,kpm1,ein,0)  ! range (mg/cm2)
             ecc = e_out(t1,bcoef1,nb1,kpm1,range,ein,thick) ! (MeV/u)
            else
             ene1 = esmaxu * 0.99d0 !FURUTA20200219 E/A -> E/u
             ene2 = esmaxu          !FURUTA20200219 E/A -> E/u
             range1 = bvalue(t1,bcoef1,nb1,kpm1,ene1,0)  ! range (mg/cm2)
             range2 = bvalue(t1,bcoef1,nb1,kpm1,ene2,0)  ! range (mg/cm2)
             ecc1 = e_out(t1,bcoef1,nb1,kpm1,range1,ene1,thick) ! (MeV/u)
             ecc2 = e_out(t1,bcoef1,nb1,kpm1,range2,ene2,thick) ! (MeV/u)

             range = range1 + (ein-ene1) / (ene2-ene1) * (range2-range1)
             ecc = ecc1 + (range-range1) / (range2-range1) * (ecc2-ecc1)
            endif

            if(ecc.gt.esminu)then !FURUTA20200219
             if( ein .le. esmaxu ) then !FURUTA20200219 E/A -> E/u
              dedxout = mp / bvalue(t1,bcoef1,nb1,kpm1,ecc,1)
              r1 = bvalue(t2,bcoef2,nb2,kpm2,ein,0)
              r2 = bvalue(t2,bcoef2,nb2,kpm2,ecc,0)
             else
              if( ierrms3.eq.0 ) then
               ErrCha =""
               ErrID = 'L:1910/R:atima/F:range.f' !W00_006_002
               call ErrWrite(ErrID,ErrCha)

               write(*,'(a)')
     &          '*** Warning: ATIMA, energy is larger than esmax'
               write(*,'(a,1p2d14.7)')
     &          '    Energy/A, esmax: ',ene/dble(mp),esmax   ! S.Abe 2016/08/08
               write(*,'(a)')
     &          '    !!! Data is same as esmax for nedisp !!!'
               ierrms3 = 1
              endif

              dedxout = mp / bvalue(t1,bcoef1,nb1,kpm1,ecc2,1)
              r1 = bvalue(t2,bcoef2,nb2,kpm2,ene2,0)
              r2 = bvalue(t2,bcoef2,nb2,kpm2,ecc2,0)
             endif

             if( r1-r2 .lt. 0.d0 ) write(*,*) ' *** ein, ecc ',ein, ecc,
     &                                                       r1, r2
             estragg = dedxout * sqrt ( abs( r1 - r2 ) ) ! (MeV)

             edif = estragg * gaurn(dummy)

            else        !FURUTA20160422
             edif=0.0d0 ! Energy straggling to be zero if ecc<=esmin
            endif       !FURUTA20160422

            ecc = max( 0.0d0, ecc * mp - edif )   ! S.Abe 2016/08/08

*-----------------------------------------------------------------------

         else if( iway .eq. 4 ) then

            thick = delt * rho * 1000.d0

            ein = ene / mp   ! S.Abe 2016/08/08
            if( ein .le. esmaxu ) then !FURUTA20200219 E/A -> E/u
             range = bvalue(t1,bcoef1,nb1,kpm1,ein,0)  ! range (mg/cm2)
             ecc = e_out(t1,bcoef1,nb1,kpm1,range,ein,thick) ! (MeV/u)

             s3 = bvalue(t3,bcoef3,nb3,kpm3,ein,0)
             s4 = bvalue(t3,bcoef3,nb3,kpm3,ecc,0)
            else

             if( ierrms4.eq.0 ) then
              ErrCha =""
              ErrID = 'L:1956/R:atima/F:range.f' !W00_006_003
              call ErrWrite(ErrID,ErrCha)

              write(*,'(a)')
     &         '*** Warning: ATIMA, energy is larger than esmax'
              write(*,'(a,1p2d14.7)')
     &         '    Energy/A, esmax: ',ene/dble(mp),esmax   ! S.Abe 2016/08/08
              write(*,'(a)')
     &         '    !!! Data is same as esmax for nspred !!!'
              ierrms4 = 1
             endif

             ene1 = esmaxu * 0.99d0 !FURUTA20200219 E/A -> E/u
             ene2 = esmaxu          !FURUTA20200219 E/A -> E/u
             range1 = bvalue(t1,bcoef1,nb1,kpm1,ene1,0)  ! range (mg/cm2)
             range2 = bvalue(t1,bcoef1,nb1,kpm1,ene2,0)  ! range (mg/cm2)
             ecc1 = e_out(t1,bcoef1,nb1,kpm1,range1,ene1,thick) ! (MeV/u)
             ecc2 = e_out(t1,bcoef1,nb1,kpm1,range2,ene2,thick) ! (MeV/u)

             range = range1 + (ein-ene1) / (ene2-ene1) * (range2-range1)
             ecc = ecc1 + (range-range1) / (range2-range1) * (ecc2-ecc1)

             s3 = bvalue(t3,bcoef3,nb3,kpm3,ene2,0)
             s4 = bvalue(t3,bcoef3,nb3,kpm3,ecc2,0)
            endif

            astragg = sqrt( abs( s3 - s4 ) )            ! (rad)

            sigx = delt * sin( astragg )
  560       r1  = gaurn(dummy)
            r2  = gaurn(dummy)
            rm = sigx * sqrt( r1**2 + r2**2 )
            am = delt**2 - rm**2
            if( am .lt. 0.0d0 ) goto 560

            ecc = rm

*-----------------------------------------------------------------------

         else if( iway .eq. 5 ) then

            ein   = ene / mp   ! S.Abe 2016/08/08

            if( ein .le. esmaxu ) then !FURUTA20200219 E/A -> E/u
             dedx  = mp / bvalue(t1,bcoef1,nb1,kpm1,ein,1)
             range = bvalue(t1,bcoef1,nb1,kpm1,ein,0)  ! range (mg/cm2)
            else
             ene1 = esmaxu * 0.99d0 !FURUTA20200219 E/A -> E/u
             ene2 = esmaxu          !FURUTA20200219 E/A -> E/u
             dedx1 = mp / bvalue(t1,bcoef1,nb1,kpm1,ene1,1)
             dedx2 = mp / bvalue(t1,bcoef1,nb1,kpm1,ene2,1)
             range1 = bvalue(t1,bcoef1,nb1,kpm1,ene1,0)  ! range (mg/cm2)
             range2 = bvalue(t1,bcoef1,nb1,kpm1,ene2,0)  ! range (mg/cm2)

             dedx = dedx1 + (ein-ene1) / (ene2-ene1) * (dedx2-dedx1)
             range = range1 + (ecc-ene1) / (ene2-ene1) * (range2-range1)
            endif

            rng = range / rho / 1000.d0             ! (cm)
            ecc = dedx * 1000.d0 * rho              ! (MeV/cm)

         end if

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine spar(ipty,ap,zp,ene,rms,mat,rng,ecc,iway)
*                                                                      *
*       control routine of SPAR                                        *
*       last modified by K.Niita on 2005/08/17                         *
*                                                                      *
*       ipty : projectile id, 1:Nucleus, 2:Proton, 3:Pion, 4:Muon      *
*       ap   : mass number of projectile                               *
*       zp   : charge of projectile                                    *
*       ene  : initial energy of projectile (MeV)                      *
*       rms  : mass of projectile  (MeV)                               *
*       mat  : material number                                         *
*       rng  : range of initial energy ene                             *
*       ecc  : final energy at rng                                     *
*                                                                      *
*       iway : 1=> range, 2=>energy                                    *
*              3=> dE/dx at the initial ene, ecc is output             *
*                                                                      *
*                                                                      *
************************************************************************
      use MEMBANKMOD !FURUTA
      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'err.inc'

      parameter (npts=32768)
      parameter ( amu = 931.141 )

*-----------------------------------------------------------------------

      real*8 mother

*-----------------------------------------------------------------------

      common /kmat1a/ mxmat, mxmat0, mxnel
      common /eparm/  esmax, esmin, emin(20)
      common /spred/  nspred, nwsprd, nedisp, itstep, ndedx
      common /sparcmn/ khb, kh1, kh2, kh3, khe(kvlmax), khp,
     &                 kdf, kdg, kro

      common /kmat1g/ kmat(kvlmax)

      DATA NFIRST/0/
      save NFIRST !FURUTA
!$OMP THREADPRIVATE(NFIRST)
*-----------------------------------------------------------------------

      dimension eeee(0:npts), rrrr(0:npts)


      data ierrms / 0 /    !S.Abe 2016/05/16
      save ierrms          !S.Abe 2016/05/16

*-----------------------------------------------------------------------

         elog(i) = exp( ( i - 1 ) * xlinc ) * esmin

*-----------------------------------------------------------------------

      if( ene .le. esmin*ap ) then !FURUTA20160422
         ecc = 1.d-11
         rng = 1.d-10
         return
      end if

      if( iway .eq. 2 .and. rng .lt. 1.0d-10 ) then
         ecc = 1.d-11
         return
      end if

cFURUTA20160422---------------------------------------------------------
      if( ene .ge. esmax*ap .and. ierrms.eq.0 ) then ! make it same as ATIMA
       ErrCha = ''
       ErrID = 'L:2102/R:spar/F:range.f' !W00_006_004
       call ErrWrite(ErrID,ErrCha)

       write(*,'(a)')'*** Warning: SPAR, energy is larger than esmax'
       write(*,'(a,1p2d14.7)')'    Energy/A, esmax: ',ene/dble(ap),esmax
       write(*,'(a)')'    !!! Data is extrapolated !!!'
       ierrms = 1
      endif
c-----------------------------------------------------------------------

*-----------------------------------------------------------------------

      if( nfirst .eq. 0 ) then

         call prepspar

         nfirst = nfirst + 1

      end if

*-----------------------------------------------------------------------

            med  = mat
            eng  = ene
            path = rng

            z0   = zp
            xm   = rms / amu

*-----------------------------------------------------------------------

         if( iway .eq. 3 ) then

*-----------------------------------------------------------------------

                  ei = eng

                  call dedxspar(ipty,xm,z0,ei,med,stpi)
                  ecc = stpi

*-----------------------------------------------------------------------

         else if( iway .eq. 1 ) then

*-----------------------------------------------------------------------

                  ei = eng

                  call rangespar(ipty,xm,z0,ei,med,rngi)
                  rng = rngi

*-----------------------------------------------------------------------

         else if( iway .eq. 2 ) then

*-----------------------------------------------------------------------

            xlinc = log( eng / esmin ) / ( npts  - 1 )

                  j = npts / 2
                  i = npts / 2

  100       continue

                  j = j / 2

                  ei = elog(i)

                  call rangespar(ipty,xm,z0,ei,med,rngi)

                  eeee(i) = ei
                  rrrr(i) = rngi

            if( rrrr(i) .le. path ) then

               if( j .eq. 0 ) then

                  i = i + 1
                  ei = elog(i)
                  call rangespar(ipty,xm,z0,ei,med,rngi)
                  eeee(i) = ei
                  rrrr(i) = rngi

                  goto 200

               end if

                  i = i + j

            else

               if( j .eq. 0 ) then

                  ei = elog(i-1)
                  call rangespar(ipty,xm,z0,ei,med,rngi)
                  eeee(i-1) = ei
                  rrrr(i-1) = rngi

                  goto 200

               end if

                  i = i - j

            end if

                  go to 100

  200       continue

*-----------------------------------------------------------------------

            if( i .gt. 1 ) then

               ee = eeee(i-1)
               r  = rrrr(i-1)

            else

               ee = 0.d0
               r  = 0.d0

            end if

*-----------------------------------------------------------------------

               ecc = ee
     &             + ( path - r ) * ( eeee(i) - ee )
     &                            / ( rrrr(i) - r )

         end if

*-----------------------------------------------------------------------

      return
      END


************************************************************************
*                                                                      *
      SUBROUTINE RANGEspar(ITYP,XM,Z,E,MED,RNG)
*                                                                      *
*        Last Revised:     2012 10 01                                  *
*                                                                      *
************************************************************************
      use MEMBANKMOD !FURUTA
*-----------------------------------------------------------------------

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'

*-----------------------------------------------------------------------

      common /kmat1a/ mxmat, mxmat0, mxnel
      common /sparcmn/ khb, kh1, kh2, kh3, khe(kvlmax), khp,
     &                 kdf, kdg, kro

      common /kmat1g/ kmat(kvlmax)

*-----------------------------------------------------------------------

      DATA Q2D3/.6666667/
      DATA NFIRST/1/
      integer,parameter:: imax=500 !FURUTA20160422

      DIMENSION ZA(IMAX),XMA(IMAX),MEDA(IMAX),ECUT(IMAX),RNGI(125,IMAX),
     X  RECUT(IMAX),PRECUT(IMAX),EINC2(IMAX),RNGILOG(125,IMAX) !FURUTA

      DIMENSION EN(625),FEN(625),AEN(625),X(1701),FX(1701),AX(1701)

*-----------------------------------------------------------------------

      save einc, einc1, ipart, elolog, el1log, el2log
!$OMP THREADPRIVATE(einc, einc1, ipart, elolog, el1log, el2log)
      save NFIRST,ZA,XMA,MEDA,ECUT,RNGI,RNGILOG,RECUT,PRECUT,EINC2,EN !FURUTA
!$OMP THREADPRIVATE(NFIRST,ZA,XMA,MEDA,ECUT,RNGI,RNGILOG,RECUT,PRECUT)
!$OMP THREADPRIVATE(EINC2,EN)
*-----------------------------------------------------------------------

      parameter (ELO   = 1.0d-10)
      parameter (NPTS  = 1701   )
      parameter (NPTS1 = NPTS-1 )
      parameter (EL1   = 1.0d-10)
      parameter (EL2   = 1.0d-5 )
      parameter (NP1   = 26     )
      parameter (NP2   = 100    )
      parameter (N1    = (NP1-1) * 5    )
      parameter (N2    = (NP2-1) * 5 + 1)
      parameter (NMAX  = NP1 + NP2 - 1  )
      parameter (JJ    = N1+N2  )

      common /eparm/  esmax, esmin, emin(20) !FURUTA20160422
*-----------------------------------------------------------------------

      IF(NFIRST.eq.1)then
        EHI = esmax !FURUTA20160422
        NFIRST = 0
        IPART = 0
        ELOLOG = LOG(ELO)
        EINC = LOG(EHI/ELO) / NPTS1
        EL1LOG = LOG(EL1)
        EL2LOG = LOG(EL2)
        EINC1 = LOG(EL2/EL1) / N1
        do N=1,N1
          EN(N) = EL1 * EXP(EINC1 * (N-1))
        enddo
        do I=1,NPTS
          X(I) = ELO * EXP((I-1) * EINC)
        enddo
        do M = 1,MXMAT
          CALL DEDXspar(2,1.d0,1.d0,ELO,M,STP)
          RNGELO = .5 * ELO / STP
          do I = 1,NPTS
            CALL DEDXspar(2,1.d0,1.d0,X(I),M,STP)
            FX(I) = 1./ STP
          enddo
          CALL SIMP2spar(X,FX,AX,NPTS)
          J = 0
          do I = 1, NPTS, 5
            J = J + 1
            PRTRNG(khp+J,M) = AX(I) + RNGELO
            PRTRNGLOG(khp+J,M) = LOG(AX(I) + RNGELO)
          enddo
        enddo
        EINC = EINC * 5.0
        EINC1 = EINC1 * 5
      ENDIF
*-----------------------------------------------------------------------

      IF(E.LE.0.0d0)then
        RNG = 0.0d0
        return
      ENDIF
      ELOG = LOG(E)

      if(ITYP.eq.2)then
        if( E .gt. esmax ) then
         I = NPTS / 5
        else
         I = 1 + (ELOG-ELOLOG) / EINC
        endif
        if(I.GE.1)then
          E1LOG = ELOLOG + (I-1)*EINC
          RNG = PRTRNGLOG(khp+I,MED) + (ELOG-E1LOG) / EINC *
     X         (PRTRNGLOG(khp+I+1,MED) - PRTRNGLOG(khp+I,MED))
          RNG = EXP(RNG)
        else
          RNG = E/ELO * PRTRNG(khp+1,MED)
        endif
        return
      elseif(ITYP.eq.1)then
        ZI = Z
        XMI = XM
      else
        ZI = 1.0
        XMI = .1488
        IF(ITYP.EQ.4) XMI = .1129
      endif

      do I = 1,IPART
        IF(XMI.NE.XMA(I))  cycle
        IF(ZI.NE.ZA(I))    cycle
        IF(MED.EQ.MEDA(I)) exit
      enddo

      if(I.gt.IPART)then
        IPART = IPART + 1
        IF(IPART.GT.IMAX) IPART = IMAX
        I = IPART
        ZA(I) = ZI
        XMA(I) = XMI
        MEDA(I) = MED

        if( zi .le. 31.0 ) then
          B = 0.07 * ZI ** Q2D3
          ESTAR = SQRT(1./(1. - B*B)) - 1.0
        else
          B = 0.09 * ZI ** Q2D3
          ESTAR = b**4 / (1.+b**4)
        end if

        FACT = 931.141
        IF(ITYP.NE.1) FACT = 938.232
        ECUT(I) = ESTAR * FACT * XMI
        EH = ECUT(I)
        EINC2(I) = LOG(EH/EL2) / (N2 - 2)
        do N = 1,N2
          EN(N1+N) = EL2 * EXP(EINC2(I)*(N-1))
        enddo
        EN(JJ-1) = ECUT(I)
        EP = ECUT(I) / XMI
        EPLOG = LOG(EP)
        K = 1 + (EPLOG-ELOLOG) / EINC
        E1LOG = ELOLOG + (K-1)*EINC
        M = MED
        IP = I
        CALL DEDXspar(ITYP,XM,Z,EL1,M ,STP)
        RNGEL1 = .5 * EL1 / STP
        do J = 1,JJ
          CALL DEDXspar(ITYP,XM,Z,EN(J),M,STP)
          FEN(J) = 1./ STP
        enddo
        CALL SIMP2spar(EN,FEN,AEN,JJ)
        L = 0

        do J = 1,JJ,5
          L = L + 1
          RNGI(L,IP) = AEN(J) + RNGEL1
          RNGILOG(L,IP) = LOG(AEN(J) + RNGEL1)
        enddo

        RECUT(IP) = AEN(JJ-1) + RNGEL1
        PRE = PRTRNGLOG(khp+K,M) + (EPLOG-E1LOG) / EINC *
     X       (PRTRNGLOG(khp+K+1,M) - PRTRNGLOG(khp+K,M))
        PRECUT(IP) = EXP(PRE)
        EINC2(I) = EINC2(I) * 5.0
      endif

      IP = I
      IF(E.GT.ECUT(I))THEN
        EP = E / XMI
        EPLOG = LOG(EP)
        if( EP .gt. esmax ) then
         K = NPTS / 5
        else
         K = 1 + (EPLOG-ELOLOG) / EINC
        endif
        E1LOG = ELOLOG + (K-1)*EINC
        PRE = PRTRNGLOG(khp+K,MED) + (EPLOG-E1LOG) / EINC *
     X       (PRTRNGLOG(khp+K+1,MED) - PRTRNGLOG(khp+K,MED))
        PRE = EXP(PRE)
        RNG = RECUT(IP) + XMI/ (ZI*ZI) * (PRE - PRECUT(IP))
      ELSE
        IF(E.LE.EL2)THEN
          K = 1 + (ELOG-EL1LOG) / EINC1
          IF(K.LT.1)THEN
            RNG = E/EL1 * RNGI(1,IP)
            return
          ENDIF
          E1LOG = EL1LOG + (K-1)*EINC1
          EINC3 = EINC1
        ELSE
          K = 1 + (ELOG-EL2LOG) / EINC2(I)
          E1LOG = EL2LOG + (K-1)*EINC2(I)
          EINC3 = EINC2(I)
          K = NP1 + K -1
        ENDIF
        RNG = RNGILOG(K,IP) + (ELOG-E1LOG) / EINC3 *
     X       (RNGILOG(K+1,IP) - RNGILOG(K,IP))
        RNG = EXP(RNG)
      ENDIF

      RETURN
      END

************************************************************************
*                                                                      *
      SUBROUTINE DEDXspar(ITYP,XM,Z,E,MED,STP)
C
C       DEDX FOR HEAVY ION, PROTONS,PIONS,MUONS
C
*                                                                      *
************************************************************************
      use MEMBANKMOD !FURUTA
      use moddas_material

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'

*-----------------------------------------------------------------------

      common /kmat1a/ mxmat, mxmat0, mxnel
      common /sparcmn/ khb, kh1, kh2, kh3, khe(kvlmax), khp,
     &                 kdf, kdg, kro


      common /kmat1g/ kmat(kvlmax)



*-----------------------------------------------------------------------


      real(8),parameter:: Q2D3=2.0d0/3.0d0
      real(8),parameter:: Q1D3=1.0d0/3.0d0

      DATA MOLD/0/,ZOLD/0./,XMOLD/0./,IOLD/0/
      DIMENSION XMF(4)
      DATA XMF/0.,1.,0.1488,0.1129/
      save MOLD,ZOLD,XMOLD,IOLD !FURUTA
!$OMP THREADPRIVATE(MOLD,ZOLD,XMOLD,IOLD)
       save XKI             !FURUTA
!$OMP THREADPRIVATE(XKI)
       DATA XKI/0./         !FURUTA

      common /paraj/  mstz(300), parz(300)
      ismode=mstz(116)
      ROsum=0.0d0

*-----------------------------------------------------------------------

       MOLD = 0
       ZOLD = 0.
       XMOLD = 0.
       IOLD = 0

*-----------------------------------------------------------------------

      IF(E.EQ.0.) GO TO 1001
      IF(ITYP.EQ.1) GO TO 120
      ZI = 1.0
      XMI = XMF(ITYP)
      XMAS = 938.232 * XMI
      BCUT1 = .07
      BCUT2 = .0046
      GO TO 130
  120 ZI = Z
      XMI = XM
      XMAS = 931.141 * XMI
      BCUT1 = .07 * Z**Q2D3
      BCUT2 = .0046 * Z**Q1D3
  130 ESTAR = E / XMAS
      TEMP1 = (ESTAR + 1.D0)**2
      ESQ = TEMP1 - 1.D0
      BSQ = ESQ / TEMP1
      B = DSQRT(BSQ)
      IF(B.LT.BCUT1) GO TO 135
      ZSQ = ZI * ZI
      SNUC = 0.
      IFROM = 1
      if(ismode.ne.1) GO TO 150 ! T.Sato 2016/2/17, should not goto 150 for shielding calculation mode
  135 CONTINUE
      ISAME = 1

      IF(ITYP.EQ.IOLD .AND.
     &   XM.EQ.XMOLD.AND.Z.EQ.ZOLD.AND.MED.EQ.MOLD) GO TO 144

  137 XMOLD = XM
      ZOLD = Z
  138 IOLD = ITYP
      MOLD = MED
      ISAME = 0

      NELM = nint( dnel_das(kmat0+med) )

      ZIQ2D3 = ZI**Q2D3
      DO 140 J = 1,NELM
      Q4 = a_das(kmat(med)+j)
      Q1 = SQRT(ZIQ2D3 + zz_das(kmat(med)+j)**Q2D3)
      Q2 = Q4 + XMI
      Q3 = ZI*zz_das(kmat(med)+j)

      F(kdf+J) = 3.255D4 * Q4 / (Q2 * Q3 * Q1)
      G(kdg+J) = 1.96D-4 * Q4 * Q2 * Q1 / (Q3 * XMI)
      RO(kro+J) = Q4*den_das(kmat(med)+j)/0.602
      ROsum=ROsum+RO(kro+J) ! T.Sato 2016/2/17
  140 CONTINUE

      IF( denh_das(kmat0+med) .EQ.0.) GO TO 142

      NELM = NELM + 1
      J = NELM

      Q2 = 1.0 + XMI
      Q1 = SQRT(ZIQ2D3 + 1.0)

      F(kdf+J) = 3.255D4 / (Q2 * ZI * Q1)
      G(kdg+J) = 1.96D-4 * Q2 * Q1 / (ZI * XMI)
      RO(kro+J) = denh_das(kmat0+med)/0.602
      ROsum=ROsum+RO(kro+J) ! T.Sato 2016/2/17
  142 CONTINUE

      if(ismode.eq.1) then ! density calculation mode, T.Sato 2016/2/17
       STP=ROsum
       return
      endif

  144 CONTINUE
      IF(B.LT.BCUT2) GO TO 200
      SNUC = 0.

      DO 148 J = 1,NELM
      EJ = F(kdf+J) * E
      IF(EJ.GT.1000.) GO TO 148
      IF(EJ.LE.4.) GO TO 146
      DEDP = 0.5455 * LOG(EJ) / (EJ * (1.0 - 0.9988 * EJ**(-1.5391)))
      GO TO 147
  146 DEDP = 4.46426 * SQRT(EJ) * EXP(-2.542 * EJ**0.277)
  147 SNUC = SNUC + RO(kro+J) * DEDP / G(kdg+J)
  148 CONTINUE

      IFROM = 1
  145 ZSTAR = ZI * (1.0 - EXP(-125.0*B/ZIQ2D3))
      ZSQ = ZSTAR * ZSTAR
  150 A6 = DLOG(ESQ) - BSQ + 0.0217615D0
      DEL = DLOG(1.378D-9 * S1(kh1+MED) * ESQ)
     &    - 2. * S2(kh2+MED)/S1(kh1+MED) - 1.
      IF(DEL.LT.0.) DEL = 0.
      CALL SHELLspar(BSQ,S1(kh1+MED),S2(kh2+MED),S3(kh3+MED),COZ)
      A7 = 5.0985D-1 * ((A6-DEL/2.-COZ)*S1(kh1+MED)- S2(kh2+MED))/ BSQ
      IF(IFROM.EQ.2) GO TO 215
      STP = ZSQ * A7 + SNUC
      GO TO 999
  200 CONTINUE
      IF(ISAME.EQ.1) GO TO 220
      IF(ITYP.EQ.1) GO TO 205
      B = 0.0046
      GO TO 210

  205 B = 0.0046 * Z ** Q1D3
  210 BSQ = B * B

      ESQ = BSQ / (1. - BSQ)
      EI = (DSQRT(1./(1.-BSQ)) -1.) * XMAS

      IFROM = 2
      GO TO 145
  215 SI = A7 * ZSQ
      XKI = SI / SQRT(EI)
  220 CONTINUE
      SE = XKI * SQRT(E)
      SNUC = 0.

      DO 248 J = 1,NELM
      EJ = F(kdf+J) * E
      IF(EJ.GT.1000.) GO TO 248
      IF(EJ.LE.4.) GO TO 246
      DEDP = 0.5455 * LOG(EJ) / (EJ * (1.0 - 0.9988 * EJ**(-1.5391)))
      GO TO 247
  246 DEDP = 4.46426 * SQRT(EJ) * EXP(-2.542 * EJ**0.277)
  247 SNUC = SNUC + RO(kro+J) * DEDP / G(kdg+J)
  248 CONTINUE

      STP = SE + SNUC

  999 RETURN
 1001 STP = 0.
      RETURN
      END


************************************************************************
*                                                                      *
      SUBROUTINE PREPspar
C
C       SETS UP VARIABLES NEEDED BY RANGE-ENERGY ROUTINES
*                                                                      *
************************************************************************
      use MEMBANKMOD !FURUTA
      use moddas_material
*-----------------------------------------------------------------------

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'

*-----------------------------------------------------------------------

      common /kmat1a/ mxmat, mxmat0, mxnel
      common /sparcmn/ khb, kh1, kh2, kh3, khe(kvlmax), khp,
     &                 kdf, kdg, kro

      common /kmat1g/ kmat(kvlmax)

*-----------------------------------------------------------------------

      DO 40 M=1,MXMAT

      NELM = nint( dnel_das(kmat0+m) )

      DO 40 L=1,NELM

      eion(khe(m)+l) = ZFOIspar(zz_das(kmat(m)+l))
      eion(khe(m)+l) = eion(khe(m)+l) * 1.d-6

   40 CONTINUE

      DO 50 M=1,MXMAT

      S1(kh1+M) = denh_das(kmat0+m)




cKN from NMTC/JAM see above !!!
      S2(kh2+M) = denh_das(kmat0+m)*(-10.956)

      S3(kh3+M) = denh_das(kmat0+m)
      NELM  = nint( dnel_das(kmat0+m) )

      DO 50 L=1,NELM
      T = den_das(kmat(m)+l)*zz_das(kmat(m)+l)

      S1(kh1+M)= S1(kh1+M) + T
      S2(kh2+M) = S2(kh2+M) + T*LOG(eion(khe(m)+l))
      S3(kh3+M) = S3(kh3+M) + den_das(kmat(m)+l)

   50 CONTINUE

      DO 60 M=1,MXMAT

         BARI(khb+M) = EXP(S2(kh2+M)/S1(kh1+M))

   60 CONTINUE

      RETURN
      END


************************************************************************
*                                                                      *
      SUBROUTINE SHELLspar(E2,SUM1,SUM2,SUM3,COZ)

C **    COMPUTES C/Z, SHELL CORRECTION TO DE/DX FORMULA
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      DATA A1/ .42237700E-6/ ,  A2/3.85801900E-9/ ,
     X     B1/ .03040430E-6/ ,  B2/-.16679890E-9/ ,
     X     C1/-.00038106E-6/ ,  C2/ .00157955E-9/
      DATA D1/0.0217615/, RFSC/137.0/
      DATA ZAL/12.95/, ZH2O/3.34/
      DATA P1,P2,P3,P4,P5/4.774248E-4, 1.143478E-3, -5.633920E-2,
     X        4.763953E-1, 4.844536E-1/
      DATA W1,W2,W3,W4,W5,W6,W7,W8 /-1.819954E-6, -2.232760E-5,
     X  1.219912E-4, 1.837873E-3, -4.457574E-3, -6.837103E-2,
     X  5.266586E-1, 3.743715E-1/

*-----------------------------------------------------------------------

      X2 = E2 / (1.D0 - E2)
      E = (DSQRT(X2 + 1.D0) - 1.D0) * 938.232D0
      BARI = EXP(SUM2/SUM1) * 1.d6
      ZBAR = SUM1/SUM3
      IF(E.LT.8.0) GO TO 30
      F1 = A1/X2 + B1/(X2*X2) + C1/X2**3
      F2 = A2/X2 + B2/(X2*X2) + C2/X2**3
      COZ =(F1*BARI*BARI + F2*BARI**3) / ZBAR
      X = E2 * RFSC * RFSC / ZBAR
      XL = DLOG(E2) + D1 - COZ - SUM2/SUM1
      GO TO 99
   30 CONTINUE
      X = E2 * RFSC * RFSC / ZBAR
      XLOG = LOG(X)
      XLOG2 = XLOG * XLOG
      XLOG3 = XLOG2 * XLOG
      XLOG4 = XLOG3 * XLOG
      IF(ZBAR.GT.ZAL) GO TO 150
      XLOG5 = XLOG4 * XLOG
      XLOG6 = XLOG5 * XLOG
      XLOG7 = XLOG6 * XLOG
      XL1 = W1*XLOG7 + W2*XLOG6 + W3*XLOG5 + W4*XLOG4 + W5*XLOG3 +
     X  W6*XLOG2 + W7*XLOG + W8
      XL1 = EXP(XL1)
      IF(ZBAR.GT.ZH2O) GO TO 150
      XL = XL1
      GO TO 80
  150 CONTINUE
      XL = P1*XLOG4 + P2*XLOG3 + P3*XLOG2 + P4*XLOG + P5
      XL = EXP(XL)
      IF(ZBAR.GT.ZAL) GO TO 80
      XL2 = XL
      XL = XL1 + (ZBAR-ZH2O) / (ZAL-ZH2O) * (XL2 - XL1)
   80 COZ = DLOG(E2) + D1 -SUM2/SUM1 - XL
   99 CONTINUE
      RETURN
      END


************************************************************************
*                                                                      *
      SUBROUTINE SIMP2spar(XX,FX,AX,MM)

C     PARABOLIC INTEGRATION EVEN OR UNEVEN SPACING
C     MM NEGATIVE RETURNS INTEGRAL FROM XX(I-1) TO XX(I) IN AX(I)
C     MM POSITIVE RETURNS INTEGRAL FROM XX(1) TO XX(I) IN AX(I)
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      DIMENSION XX(MM),FX(MM),AX(MM) ! S.H. corrected (2022.10.31)
      NX=MM-1
      IF(MM)10,10,11
   10 NX=-MM-1
   11 AX(1)=0.0
      I=2
      D1=XX(I)-XX(I-1)
      D2=XX(I+1)-XX(I)
      D3=XX(I+1)-XX(I-1)
      A3=-0.16666667*D1/D3*D1/D2*D1
      A2=-A3/D1*D3+0.5*D1
      A1=-A2-A3+D1
      AX(I)=AX(I-1)+FX(I-1)*A1+FX(I)*A2+FX(I+1)*A3
      DO 20 I=2,NX
      D1=XX(I)-XX(I-1)
      D2=XX(I+1)-XX(I)
      D3=XX(I+1)-XX(I-1)
      A1=-0.16666667*D2/D3*D2/D1*D2
      A2=-A1/D2*D3+0.5*D2
      A3=-A1-A2+D2
   18 AX(I+1)=FX(I-1)*A1+FX(I)*A2+FX(I+1)*A3
      IF(MM)20,20,19
   19 AX(I+1)=AX(I+1)+AX(I)
   20 CONTINUE
      RETURN
      END


************************************************************************
*                                                                      *
      FUNCTION ZFOIspar(Z)
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------
      DIMENSION FLI(13)
      DATA FLI/18.7,42.0,39.0,60.0,68.0,78.0,99.5,98.5,117.0,140.0,
     1150.0,157.0,163.0/,A1/9.76/,A2/58.8/,XP/-0.19/
      IF(Z.GT.13.) GO TO 10
      IZ=NINT(Z)
      ZFOIspar=FLI(IZ)
      RETURN
   10 ZFOIspar=A1*Z + A2 *(Z**XP)
      RETURN
      END


************************************************************************
*                                                                      *
      subroutine glando(step,zp,xmass,t,z,a,rho,de,iflag)
*                                                                      *
*       modified by K.Niita on 2003/10/07                              *
*                                                                      *
************************************************************************
*                                                                      *
*        glando                                                        *
*        ------                                                        *
*        energy straggling using gaussian, landau & vavilov theories.  *
*                                                                      *
*        input                                                         *
*        -----                                                         *
*                                                                      *
*        step  =  current step-length (cm)                             *
*                                                                      *
*        zp    = atomic number ( +/- charge ) of projectile particle   *
*                                                                      *
*        xmass = projectile mass (MeV)                                 *
*                                                                      *
*        t     = kinetic energy (MeV)                                  *
*                                                                      *
*        z     = atomic number of medium                               *
*                                                                      *
*        a     = atomic weight of medium                               *
*                                                                      *
*        rho   = mass density of medium (gm/cm**3)                     *
*                                                                      *
*        output                                                        *
*        ------                                                        *
*        de     =  ( de/dx - <de/dx> ) * step   (MeV)                  *
*                                                                      *
*        iflag  =  0   step too small, CSDA assumed                    *
*               =  4   landau sampling used (cappa < 0.01)             *
*               =  5   vavilov sampling used (0.01 <= cappa <= 10.0)   *
*               =  6   gaussian sampling used (cappa > 10.0)           *
*                                                                      *
*        warning (as stated, but unconfirmed locally)                  *
*        --------------------------------------------                  *
*        only landau sampling should be used since this has been well  *
*        tested whereas both vavilov and gaussian sampling are being   *
*        developed.                                                    *
*                                                                      *
*        author      : g.n. patrick                                    *
*        date        : 03.05.1985                                      *
*        last update : 09.09.1985                                      *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      parameter (euler=.5772156649015328606065120900824d0)
      parameter (pi=3.14159265358979323846264338328d0)
      parameter (p1=.60715d+0,p2=.67794d+0,p3=.52382d-1,p4=.94753d+0,
     &           p5=.74442d+0,p6=1.1934d+0)

      parameter (emass=0.5109990615d+0)

*-----------------------------------------------------------------------

      fland(x) = p1 + p6 * x + ( p2 + p3 * x ) * exp( p4 * x + p5 )

*-----------------------------------------------------------------------

      if(step*rho.lt.0.0000001)then
          de = 0.
          iflag = 0
          return
      endif

*-----------------------------------------------------------------------
*     calculate xi factor (kev).
*-----------------------------------------------------------------------

      s     = t / xmass
      eta   = sqrt( s * ( s + 2. ) )
      gamma = s + 1.
      beta  = eta / gamma
      xi    = ( 153.5 * zp * zp * z * step * rho )
     &      / ( a * beta * beta )

*-----------------------------------------------------------------------
*     maximum energy transfer to atomic electron (kev).
*-----------------------------------------------------------------------

      etasq = eta * eta
      ratio = emass / xmass
      f1    = 2. * emass * etasq
      f2    = 1. + 2. * ratio * gamma + ratio * ratio
      emax  = 1000. * f1 / f2

*-----------------------------------------------------------------------
*     calculate kappa significance ratio.
*-----------------------------------------------------------------------

      cappa  = xi / emax

*-----------------------------------------------------------------------
*     choose correct straggling function.
*-----------------------------------------------------------------------

      if (cappa.lt.0.01) then

*-----------------------------------------------------------------------
*     sample lambda variable from james & hancock landau distribution
*-----------------------------------------------------------------------

          iflag = 4
          xmean = -beta**2-log(xi/emax)+euler-1.
          xlamx = fland(xmean)
  25      call glandg(xlamb)
          if(xlamb.gt.xlamx) go to 25

*-----------------------------------------------------------------------
*     sample lambda variable (landau not vavilov) from
*     rotondi & montagna & kolbig vavilov distribution
*-----------------------------------------------------------------------

      else if (cappa.le.10.) then

          iflag = 5
          rka = cappa
          be2 = beta*beta
          ra = unirn(dummy)
          xlamb = gvaviv(rka,be2,ra)

*-----------------------------------------------------------------------
*     sample from gaussian distribution
*-----------------------------------------------------------------------

      else

          iflag = 6
          sigma  = xi*emax*(1.-(beta*beta/2.))
          sigma  = sqrt(sigma)
   30     r1 = unirn(dummy)
          r2 = unirn(dummy)
          if(r1.le.0.)go to 30
          f1     = -2. * log(r1)
          dekev  = sigma * sqrt(f1) * cos(2.*pi*r2)
          goto 40
      endif

*-----------------------------------------------------------------------
*     calculate ( de/dx - <de/dx> ) * step
*-----------------------------------------------------------------------

      dekev = xi * ( xlamb + beta *beta + log( xi / emax ) - euler + 1.)

   40 de  =0.001 * dekev

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine glandg(yran)
*                                                                      *
*       modified by K.Niita on 2003/10/01                              *
*                                                                      *
************************************************************************
*                                                                      *
*       copy of the cern library routine genlan                        *
*       generation of landau-distributed random numbers by             *
*       4-point interpolation in the previously-tabulated              *
*       inverse cumulative distribution.                               *
*                                                                      *
*       ==>called by : glando                                          *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      parameter (range1=0.807069d0, range2=0.994869d0)
      parameter (x1bot =0.0d0, x2bot =0.79124d0)
      parameter (gp1inv=118.836d0, gp2inv=486.178d0)

*-----------------------------------------------------------------------

      dimension xcum1(100), xcum2(100)

*-----------------------------------------------------------------------
*     table 1 of inverse cumulative landau points (0 < p < .833)
*-----------------------------------------------------------------------

      data (xcum1(i),i=1,100) /
     & -2.5700d0,-2.15412d0,-1.94167d0,-1.79583d0,-1.67975d0,-1.58090d0,
     &-1.49341d0,-1.41397d0,-1.34057d0,-1.27185d0,-1.20686d0,-1.14490d0,
     &-1.08545d0,-1.02809d0,-0.97249d0,-0.91839d0,-0.86557d0,-0.81386d0,
     &-0.76308d0,-0.71311d0,-0.66384d0,-0.61515d0,-0.56696d0,-0.51919d0,
     &-0.47176d0,-0.42462d0,-0.37769d0,-0.33092d0,-0.28426d0,-0.23767d0,
     &-0.19109d0,-0.14447d0,-0.09779d0,-0.05100d0,-0.00405d0, 0.04308d0,
     & 0.09044d0, 0.13806d0, 0.18598d0, 0.23423d0, 0.28285d0, 0.33187d0,
     & 0.38134d0, 0.43128d0, 0.48174d0, 0.53275d0, 0.58435d0, 0.63658d0,
     & 0.68948d0, 0.74309d0, 0.79746d0, 0.85263d0, 0.90865d0, 0.96556d0,
     & 1.02342d0, 1.08228d0, 1.14219d0, 1.20322d0, 1.26542d0, 1.32887d0,
     & 1.39362d0, 1.45976d0, 1.52736d0, 1.59650d0, 1.66727d0, 1.73976d0,
     & 1.81407d0, 1.89032d0, 1.96860d0, 2.04905d0, 2.13180d0, 2.21699d0,
     & 2.30477d0, 2.39531d0, 2.48878d0, 2.58540d0, 2.68535d0, 2.78889d0,
     & 2.89626d0, 3.00775d0, 3.12364d0, 3.24429d0, 3.37005d0, 3.50136d0,
     & 3.63866d0, 3.78246d0, 3.93334d0, 4.09194d0, 4.25900d0, 4.43533d0,
     & 4.62186d0, 4.81960d0, 5.02974d0, 5.25368d0, 5.49312d0, 5.74987d0,
     & 6.02605d0, 6.32428d0, 6.64773d0, 7.00000d0/

*-----------------------------------------------------------------------
*     table 2 of inverse cumulative landau points (.791 < p < .995)
*-----------------------------------------------------------------------

      data (xcum2(i),i=1,100) /
     & 5.50000d0, 5.56120d0, 5.62347d0, 5.68684d0, 5.75133d0, 5.81699d0,
     & 5.88383d0, 5.95191d0, 6.02125d0, 6.09190d0, 6.16391d0, 6.23732d0,
     & 6.31219d0, 6.38855d0, 6.46646d0, 6.54597d0, 6.62712d0, 6.70998d0,
     & 6.79460d0, 6.88103d0, 6.96935d0, 7.05962d0, 7.15192d0, 7.24633d0,
     & 7.34294d0, 7.44182d0, 7.54306d0, 7.64676d0, 7.75300d0, 7.86188d0,
     & 7.97351d0, 8.08800d0, 8.20548d0, 8.32610d0, 8.44997d0, 8.57725d0,
     & 8.70808d0, 8.84262d0, 8.98103d0, 9.12349d0, 9.27021d0, 9.42141d0,
     & 9.57730d0, 9.73812d0, 9.90410d0,10.07552d0,10.25265d0,10.43584d0,
     &10.62540d0,10.82169d0,11.02508d0,11.23598d0,11.45487d0,11.68222d0,
     &11.91855d0,12.16441d0,12.42045d0,12.68734d0,12.96580d0,13.25663d0,
     &13.56075d0,13.87912d0,14.21280d0,14.56296d0,14.93088d0,15.31799d0,
     &15.72593d0,16.15650d0,16.61170d0,17.09378d0,17.60523d0,18.14884d0,
     &18.72792d0,19.34622d0,20.00796d0,20.71792d0,21.48176d0,22.30618d0,
     &23.19862d0,24.16809d0,25.22547d0,26.38320d0,27.65658d0,29.06490d0,
     &30.63048d0,32.38211d0,34.35555d0,36.59607d0,39.16212d0,42.13195d0,
     &45.61204d0,49.74758d0,54.74189d0,60.89935d0,68.67370d0,78.81480d0,
     &92.61047d0,112.50807d0,143.78539d0,200.0d0/

*-----------------------------------------------------------------------

      x = unirn(dummy)
      if(x .lt. 0.004) then

*-----------------------------------------------------------------------
*         extreme left-hand tail
*-----------------------------------------------------------------------

          yran = -sqrt(abs(log(x)))

*-----------------------------------------------------------------------
*         4-point interpolation in the first cumulative table
*-----------------------------------------------------------------------

      elseif(x.le.range1) then

          tabpo1 = (x-x1bot)*gp1inv
          j = tabpo1 + 1
          j = max(j,2)
          j = min(j,98)
          p = tabpo1 - (j-1)
          a = (p+1.0) * xcum1(j+2) - (p-2.0) * xcum1(j-1)
          b = (p-1.0) * xcum1(j)   - p       * xcum1(j+1)
          yran = a*p*(p-1.0)*0.16666667+b*(p+1.0)*(p-2.0)*0.5

*-----------------------------------------------------------------------
*         4-point interpolation in the first cumulative table
*-----------------------------------------------------------------------

      elseif(x.le.range2) then

          tabpo2 = (x-x2bot)*gp2inv
          j = tabpo2 + 1
          j = max(j,2)
          j = min(j,98)
          p = tabpo2 - (j-1)
          a = (p+1.0) * xcum2(j+2) - (p-2.0) * xcum2(j-1)
          b = (p-1.0) * xcum2(j)   - p       * xcum2(j+1)
          yran = a*p*(p-1.0)*0.16666667+b*(p+1.0)*(p-2.0)*0.5

*-----------------------------------------------------------------------
*         1/x**2 sampling for extreme landau tail
*-----------------------------------------------------------------------

      else

          yran = 200. / unirn(dummy)

      endif

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      function gvaviv(rkappa,beta2,ran)
*                                                                      *
*       modified by K.Niita on 2003/10/01                              *
*                                                                      *
************************************************************************
*                                                                      *
*       gvaviv                                                         *
*       ------                                                         *
*                                                                      *
*       initializes the parameters for a vavilov distribution          *
*                                                                      *
*       beta2 is the particle velocity squared                         *
*       rkappa the usual straggling parameter                          *
*                                                                      *
*       this routine has been extracted from a proposed cernlib        *
*       set of routines (vavcoe,vavfsm).                               *
*       the authors of these routines have submitted an article        *
*       to nim:                                                        *
*                                                                      *
*         alberto rotondi, paolo montagna                              *
*         fast calculation of vavilov distribution nimb 1990           *
*                                                                      *
*            author    k.s.koelbig     ********                        *
*                                                                      *
*        ==> called by : glando                                        *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      dimension ac(0:12),hc(0:8),h(9)

      parameter
     &(bkmnx1=0.02d0, bkmny1=0.05d0, bkmnx2=0.12d0, bkmny2=0.05d0,
     & bkmnx3=0.22d0, bkmny3=0.05d0, bkmxx1=0.1d0,  bkmxy1=1.0d0,
     & bkmxx2=0.2d0,  bkmxy2=1d0,    bkmxx3=0.3d0,  bkmxy3=1.0d0,
     & zero=0.0d0,two=2.0d0)

      parameter
     &(fbkx1=two/(bkmxx1-bkmnx1), fbkx2=two/(bkmxx2-bkmnx2),
     & fbkx3=two/(bkmxx3-bkmnx3), fbky1=two/(bkmxy1-bkmny1),
     & fbky2=two/(bkmxy2-bkmny2), fbky3=two/(bkmxy3-bkmny3))

      dimension edgec(2:7),fninv(5),drk(5),dsigm(5),alfa(2:5)
      dimension u1(13),u2(13),u3(13),u4(12),u5(13),u6(13),u7( 8),u8(13)
      dimension v1(12),v2(12),v3(13),v4(12),v5(13),v6(13),v7(11),v8(11)
      dimension w1(13),w2(11),w3(13),w4(13),w5(13),w6(13),       w8( 8)

      data fninv /1.0d0,0.5d0,0.33333333d0,0.25d0,0.2d0/

      data (edgec(j),j=2,7)
     &/ 0.16666667d+0, 0.41666667d-1, 0.83333333d-2,
     &  0.13888889d-1, 0.69444444d-2, 0.77160493d-3/
      data (u1(k),k=1,13)
     &/ 0.25850868d+0,  0.32477982d-1, -0.59020496d-2,
     &  0.0d0         , 0.24880692d-1,  0.47404356d-2,
     & -0.74445130d-3,  0.73225731d-2,  0.0d0        ,
     &  0.11668284d-2,  0.0d0        , -0.15727318d-2,-0.11210142d-2/
      data (u2(k),k=1,13)
     &/ 0.43142611d+0,  0.40797543d-1, -0.91490215d-2,
     &  0.0d0        ,  0.42127077d-1,  0.73167928d-2,
     & -0.14026047d-2,  0.16195241d-1,  0.24714789d-2,
     &  0.20751278d-2,  0.0d0        , -0.25141668d-2,-0.14064022d-2/
      data (u3(k),k=1,13)
     &/ 0.25225955d+0,  0.64820468d-1, -0.23615759d-1,
     &  0.0d0        ,  0.23834176d-1,  0.21624675d-2,
     & -0.26865597d-2, -0.54891384d-2,  0.39800522d-2,
     &  0.48447456d-2, -0.89439554d-2, -0.62756944d-2,-0.24655436d-2/
      data (u4(k),k=1,12)
     &/ 0.12593231d+1, -0.20374501d+0,  0.95055662d-1,
     & -0.20771531d-1, -0.46865180d-1, -0.77222986d-2,
     &  0.32241039d-2,  0.89882920d-2, -0.67167236d-2,
     & -0.13049241d-1,  0.18786468d-1,  0.14484097d-1/
      data (u5(k),k=1,13)
     &/-0.24864376d-1, -0.10368495d-2,  0.14330117d-2,
     &  0.20052730d-3,  0.18751903d-2,  0.12668869d-2,
     &  0.48736023d-3,  0.34850854d-2,  0.0d0        ,
     & -0.36597173d-3,  0.19372124d-2,  0.70761825d-3, 0.46898375d-3/
      data (u6(k),k=1,13)
     &/ 0.35855696d-1, -0.27542114d-1,  0.12631023d-1,
     & -0.30188807d-2, -0.84479939d-3,  0.0d0        ,
     &  0.45675843d-3, -0.69836141d-2,  0.39876546d-2,
     & -0.36055679d-2,  0.0d0        ,  0.15298434d-2, 0.19247256d-2/
      data (u7(k),k=1,8)
     &/ 0.10234691d+2, -0.35619655d+1,  0.69387764d+0,
     & -0.14047599d+0, -0.19952390d+1, -0.45679694d+0,
     &  0.0d0        ,  0.50505298d+0/
      data (u8(k),k=1,13)
     &/ 0.21487518d+2, -0.11825253d+2,  0.43133087d+1,
     & -0.14500543d+1, -0.34343169d+1, -0.11063164d+1,
     & -0.21000819d+0,  0.17891643d+1, -0.89601916d+0,
     &  0.39120793d+0,  0.73410606d+0,  0.0d0        ,-0.32454506d+0/
      data (v1(k),k=1,12)
     &/ 0.27827257d+0, -0.14227603d-2,  0.24848327d-2,
     &  0.0d0        ,  0.45091424d-1,  0.80559636d-2,
     & -0.38974523d-2,  0.0d0        , -0.30634124d-2,
     &  0.75633702d-3,  0.54730726d-2,  0.19792507d-2/
      data (v2(k),k=1,12)
     &/ 0.41421789d+0, -0.30061649d-1,  0.52249697d-2,
     &  0.0d0        ,  0.12693873d+0,  0.22999801d-1,
     & -0.86792801d-2,  0.31875584d-1, -0.61757928d-2,
     &  0.0d0        ,  0.19716857d-1,  0.32596742d-2/
      data (v3(k),k=1,13)
     &/ 0.20191056d+0, -0.46831422d-1,  0.96777473d-2,
     & -0.17995317d-2,  0.53921588d-1,  0.35068740d-2,
     & -0.12621494d-1, -0.54996531d-2, -0.90029985d-2,
     &  0.34958743d-2,  0.18513506d-1,  0.68332334d-2,-0.12940502d-2/
      data (v4(k),k=1,12)
     &/ 0.13206081d+1,  0.10036618d+0, -0.22015201d-1,
     &  0.61667091d-2, -0.14986093d+0, -0.12720568d-1,
     &  0.24972042d-1, -0.97751962d-2,  0.26087455d-1,
     & -0.11399062d-1, -0.48282515d-1, -0.98552378d-2/
      data (v5(k),k=1,13)
     &/ 0.16435243d-1,  0.36051400d-1,  0.23036520d-2,
     & -0.61666343d-3, -0.10775802d-1,  0.51476061d-2,
     &  0.56856517d-2, -0.13438433d-1,  0.0d0        ,
     &  0.0d0        , -0.25421507d-2,  0.20169108d-2,-0.15144931d-2/
      data (v6(k),k=1,13)
     &/ 0.33432405d-1,  0.60583916d-2, -0.23381379d-2,
     &  0.83846081d-3, -0.13346861d-1, -0.17402116d-2,
     &  0.21052496d-2,  0.15528195d-2,  0.21900670d-2,
     & -0.13202847d-2, -0.45124157d-2, -0.15629454d-2, 0.22499176d-3/
      data (v7(k),k=1,11)
     &/ 0.54529572d+1, -0.90906096d+0,  0.86122438d-1,
     &  0.0d0        , -0.12218009d+1, -0.32324120d+0,
     & -0.27373591d-1,  0.12173464d+0,  0.0d0        ,
     &  0.0d0        ,  0.40917471d-1/
      data (v8(k),k=1,11)
     &/ 0.93841352d+1, -0.16276904d+1,  0.16571423d+0,
     &  0.0d0        , -0.18160479d+1, -0.50919193d+0,
     & -0.51384654d-1,  0.21413992d+0,  0.0d0        ,
     &  0.0d0        ,  0.66596366d-1/
      data (w1(k),k=1,13)
     &/ 0.29712951d+0,  0.97572934d-2,  0.0d0        ,
     & -0.15291686d-2,  0.35707399d-1,  0.96221631d-2,
     & -0.18402821d-2, -0.49821585d-2,  0.18831112d-2,
     &  0.43541673d-2,  0.20301312d-2, -0.18723311d-2,-0.73403108d-3/
      data (w2(k),k=1,11)
     &/ 0.40882635d+0,  0.14474912d-1,  0.25023704d-2,
     & -0.37707379d-2,  0.18719727d+0,  0.56954987d-1,
     &  0.0d0        ,  0.23020158d-1,  0.50574313d-2,
     &  0.94550140d-2,  0.19300232d-1/
      data (w3(k),k=1,13)
     &/ 0.16861629d+0,  0.0d0        ,  0.36317285d-2,
     & -0.43657818d-2,  0.30144338d-1,  0.13891826d-1,
     & -0.58030495d-2, -0.38717547d-2,  0.85359607d-2,
     &  0.14507659d-1,  0.82387775d-2, -0.10116105d-1,-0.55135670d-2/
      data (w4(k),k=1,13)
     &/ 0.13493891d+1, -0.26863185d-2, -0.35216040d-2,
     &  0.24434909d-1, -0.83447911d-1, -0.48061360d-1,
     &  0.76473951d-2,  0.24494430d-1, -0.16209200d-1,
     & -0.37768479d-1, -0.47890063d-1,  0.17778596d-1, 0.13179324d-1/
      data (w5(k),k=1,13)
     &/ 0.10264945d+0,  0.32738857d-1,  0.0d0        ,
     &  0.43608779d-2, -0.43097757d-1, -0.22647176d-2,
     &  0.94531290d-2, -0.12442571d-1, -0.32283517d-2,
     & -0.75640352d-2, -0.88293329d-2,  0.52537299d-2, 0.13340546d-2/
      data (w6(k),k=1,13)
     &/ 0.29568177d-1, -0.16300060d-2, -0.21119745d-3,
     &  0.23599053d-2, -0.48515387d-2, -0.40797531d-2,
     &  0.40403265d-3,  0.18200105d-2, -0.14346306d-2,
     & -0.39165276d-2, -0.37432073d-2,  0.19950380d-2, 0.12222675d-2/
      data (w8(k),k=1,8)
     &/ 0.66184645d+1, -0.73866379d+0,  0.44693973d-1,
     &  0.0d0        , -0.14540925d+1, -0.39529833d+0,
     & -0.44293243d-1,  0.88741049d-1/

*-----------------------------------------------------------------------

      gvaviv=0

      if(rkappa .lt. 0.01 .or. rkappa .gt. 12) return

      if(rkappa .ge. 0.23) then
         itype=1
         npt=100
         wk=1/sqrt(rkappa)
         ac(0)=(-0.032227*beta2-0.074275)*rkappa+
     &    (0.24533*beta2+0.070152)*wk+(-0.55610*beta2-3.1579)
         ac(8)=(-0.013483*beta2-0.048801)*rkappa+
     &    (-1.6921*beta2+8.3656)*wk+(-0.73275*beta2-3.5226)
         drk(1)=wk**2
         dsigm(1)=sqrt(rkappa/(1-0.5*beta2))
         do 1 j = 1,4
            drk(j+1)=drk(1)*drk(j)
            dsigm(j+1)=dsigm(1)*dsigm(j)
    1       alfa(j+1)=(fninv(j)-beta2*fninv(j+1))*drk(j)
         hc(0)=log(rkappa)+beta2+0.42278434
         hc(1)=dsigm(1)
         hc(2)=alfa(3)*dsigm(3)
         hc(3)=(3*alfa(2)**2+alfa(4))*dsigm(4)-3
         hc(4)=(10*alfa(2)*alfa(3)+alfa(5))*dsigm(5)-10*hc(2)
         hc(5)=hc(2)**2
         hc(6)=hc(2)*hc(3)
         hc(7)=hc(2)*hc(5)
         do 2 j = 2,7
    2       hc(j)=edgec(j)*hc(j)
         hc(8)=0.39894228*hc(1)

      elseif(rkappa .ge. 0.22) then

         itype=2
         npt=150
         x=1+(rkappa-bkmxx3)*fbkx3
         y=1+(sqrt(beta2)-bkmxy3)*fbky3
         xx=2*x
         yy=2*y
         x2=xx*x-1
         x3=xx*x2-x
         y2=yy*y-1
         y3=yy*y2-y
         xy=x*y
         p2=x2*y
         p3=x3*y
         q2=y2*x
         q3=y3*x
         pq=x2*y2
         ac(1)=w1(1)+w1(2)*x+w1(4)*x3+w1(5)*y+w1(6)*y2+w1(7)*y3+
     &    w1(8)*xy+w1(9)*p2+w1(10)*p3+w1(11)*q2+w1(12)*q3+w1(13)*pq
         ac(2)=w2(1)+w2(2)*x+w2(3)*x2+w2(4)*x3+w2(5)*y+w2(6)*y2+
     &    w2(8)*xy+w2(9)*p2+w2(10)*p3+w2(11)*q2
         ac(3)=w3(1)+w3(3)*x2+w3(4)*x3+w3(5)*y+w3(6)*y2+w3(7)*y3+
     &    w3(8)*xy+w3(9)*p2+w3(10)*p3+w3(11)*q2+w3(12)*q3+w3(13)*pq
         ac(4)=w4(1)+w4(2)*x+w4(3)*x2+w4(4)*x3+w4(5)*y+w4(6)*y2+w4(7)*
     &    y3+w4(8)*xy+w4(9)*p2+w4(10)*p3+w4(11)*q2+w4(12)*q3+w4(13)*pq
         ac(5)=w5(1)+w5(2)*x+w5(4)*x3+w5(5)*y+w5(6)*y2+w5(7)*y3+
     &    w5(8)*xy+w5(9)*p2+w5(10)*p3+w5(11)*q2+w5(12)*q3+w5(13)*pq
         ac(6)=w6(1)+w6(2)*x+w6(3)*x2+w6(4)*x3+w6(5)*y+w6(6)*y2+w6(7)*
     &    y3+w6(8)*xy+w6(9)*p2+w6(10)*p3+w6(11)*q2+w6(12)*q3+w6(13)*pq
         ac(8)=w8(1)+w8(2)*x+w8(3)*x2+w8(5)*y+w8(6)*y2+w8(7)*y3+w8(8)*xy
         ac(0)=-3.05

       elseif(rkappa .ge. 0.10) then

         itype=3
         npt=200
         x=1+(rkappa-bkmxx2)*fbkx2
         y=1+(sqrt(beta2)-bkmxy2)*fbky2
         xx=2*x
         yy=2*y
         x2=xx*x-1
         x3=xx*x2-x
         y2=yy*y-1
         y3=yy*y2-y
         xy=x*y
         p2=x2*y
         p3=x3*y
         q2=y2*x
         q3=y3*x
         pq=x2*y2
         ac(1)=v1(1)+v1(2)*x+v1(3)*x2+v1(5)*y+v1(6)*y2+v1(7)*y3+
     &    v1(9)*p2+v1(10)*p3+v1(11)*q2+v1(12)*q3
         ac(2)=v2(1)+v2(2)*x+v2(3)*x2+v2(5)*y+v2(6)*y2+v2(7)*y3+
     &    v2(8)*xy+v2(9)*p2+v2(11)*q2+v2(12)*q3
         ac(3)=v3(1)+v3(2)*x+v3(3)*x2+v3(4)*x3+v3(5)*y+v3(6)*y2+v3(7)*
     &    y3+v3(8)*xy+v3(9)*p2+v3(10)*p3+v3(11)*q2+v3(12)*q3+v3(13)*pq
         ac(4)=v4(1)+v4(2)*x+v4(3)*x2+v4(4)*x3+v4(5)*y+v4(6)*y2+v4(7)*
     &    y3+v4(8)*xy+v4(9)*p2+v4(10)*p3+v4(11)*q2+v4(12)*q3
         ac(5)=v5(1)+v5(2)*x+v5(3)*x2+v5(4)*x3+v5(5)*y+v5(6)*y2+v5(7)*
     &    y3+v5(8)*xy+v5(11)*q2+v5(12)*q3+v5(13)*pq
         ac(6)=v6(1)+v6(2)*x+v6(3)*x2+v6(4)*x3+v6(5)*y+v6(6)*y2+v6(7)*
     &    y3+v6(8)*xy+v6(9)*p2+v6(10)*p3+v6(11)*q2+v6(12)*q3+v6(13)*pq
         ac(7)=v7(1)+v7(2)*x+v7(3)*x2+v7(5)*y+v7(6)*y2+v7(7)*y3+
     &    v7(8)*xy+v7(11)*q2
         ac(8)=v8(1)+v8(2)*x+v8(3)*x2+v8(5)*y+v8(6)*y2+v8(7)*y3+
     &    v8(8)*xy+v8(11)*q2
         ac(0)=-3.04

      else

         itype=4
         if(rkappa .ge. 0.02) itype=3
         npt=200
         x=1+(rkappa-bkmxx1)*fbkx1
         y=1+(sqrt(beta2)-bkmxy1)*fbky1
         xx=2*x
         yy=2*y
         x2=xx*x-1
         x3=xx*x2-x
         y2=yy*y-1
         y3=yy*y2-y
         xy=x*y
         p2=x2*y
         p3=x3*y
         q2=y2*x
         q3=y3*x
         pq=x2*y2

         if(itype .eq. 4) go to 4

         ac(1)=u1(1)+u1(2)*x+u1(3)*x2+u1(5)*y+u1(6)*y2+u1(7)*y3+
     &    u1(8)*xy+u1(10)*p3+u1(12)*q3+u1(13)*pq
         ac(2)=u2(1)+u2(2)*x+u2(3)*x2+u2(5)*y+u2(6)*y2+u2(7)*y3+
     &    u2(8)*xy+u2(9)*p2+u2(10)*p3+u2(12)*q3+u2(13)*pq
         ac(3)=u3(1)+u3(2)*x+u3(3)*x2+u3(5)*y+u3(6)*y2+u3(7)*y3+
     &    u3(8)*xy+u3(9)*p2+u3(10)*p3+u3(11)*q2+u3(12)*q3+u3(13)*pq
         ac(4)=u4(1)+u4(2)*x+u4(3)*x2+u4(4)*x3+u4(5)*y+u4(6)*y2+u4(7)*
     &    y3+u4(8)*xy+u4(9)*p2+u4(10)*p3+u4(11)*q2+u4(12)*q3
         ac(5)=u5(1)+u5(2)*x+u5(3)*x2+u5(4)*x3+u5(5)*y+u5(6)*y2+u5(7)*
     &    y3+u5(8)*xy+u5(10)*p3+u5(11)*q2+u5(12)*q3+u5(13)*pq
         ac(6)=u6(1)+u6(2)*x+u6(3)*x2+u6(4)*x3+u6(5)*y+u6(7)*y3+
     &    u6(8)*xy+u6(9)*p2+u6(10)*p3+u6(12)*q3+u6(13)*pq
    4    ac(7)=u7(1)+u7(2)*x+u7(3)*x2+u7(4)*x3+u7(5)*y+u7(6)*y2+u7(8)*xy
         ac(8)=u8(1)+u8(2)*x+u8(3)*x2+u8(4)*x3+u8(5)*y+u8(6)*y2+u8(7)*
     &    y3+u8(8)*xy+u8(9)*p2+u8(10)*p3+u8(11)*q2+u8(13)*pq

         ac(0)=-3.03

      endif

      ac(9)=(ac(8)-ac(0))/npt

      if(itype .eq. 3) then
         x=(ac(7)-ac(8))/(ac(7)*ac(8))
         y=1/log(ac(8)/ac(7))
         p2=ac(7)**2
         ac(11)=p2*(ac(1)*exp(-ac(2)*(ac(7)+ac(5)*p2)-
     &    ac(3)*exp(-ac(4)*(ac(7)+ac(6)*p2)))-0.045*y/ac(7))/
     &    (1+x*y*ac(7))
         ac(12)=(0.045+x*ac(11))*y
      endif

      if(itype .eq. 4) ac(10)=0.995/glands(ac(8))
      t=2*ran/ac(9)
      rlam=ac(0)
      fl=0
      s=0

      do 21 n = 1,npt
         rlam=rlam+ac(9)
         if(itype .eq. 1) then
            fn=1
            x=(rlam+hc(0))*hc(1)
            h(1)=x
            h(2)=x**2-1
            do 31 k = 2,8
               fn=fn+1
   31          h(k+1)=x*h(k)-fn*h(k-1)
            y=1+hc(7)*h(9)
            do 32 k = 2,6
   32          y=y+hc(k)*h(k+1)
            fu=hc(8)*exp(-0.5*x**2)*max(y,zero)
         elseif(itype .eq. 2) then
            x=rlam**2
            fu=ac(1)*exp(-ac(2)*(rlam+ac(5)*x)-
     &       ac(3)*exp(-ac(4)*(rlam+ac(6)*x)))
         elseif(itype .eq. 3) then
            if(rlam .lt. ac(7)) then
               x=rlam**2
               fu=ac(1)*exp(-ac(2)*(rlam+ac(5)*x)-
     &          ac(3)*exp(-ac(4)*(rlam+ac(6)*x)))
            else
               x=1/rlam
               fu=(ac(11)*x+ac(12))*x
            endif
         else
            fu=ac(10)*glande(rlam)
         endif
         s=s+fl+fu
         if(s .gt. t) go to 22
   21    fl=fu

   22 s0=s-fl-fu

      gvaviv=rlam-ac(9)

      if(s .gt. s0) gvaviv=gvaviv+ac(9)*(t-s0)/(s-s0)

      return
      end


************************************************************************
*                                                                      *
      function glands(x)
*                                                                      *
*       modified by K.Niita on 2003/10/01                              *
*                                                                      *
************************************************************************
*                                                                      *
*       copy of the cern library routine dstlan (g110)                 *
*                                                                      *
*         ==>called by : gvaviv                                        *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      dimension p1(0:4),p2(0:3),p3(0:3),p4(0:3),p5(0:3),p6(0:3)
      dimension q1(0:4),q2(0:3),q3(0:3),q4(0:3),q5(0:3),q6(0:3)
      dimension a1(1:3),a2(1:3)

      data (p1(i),i=0,4),(q1(j),j=0,4)
     & /0.2514091491d+0,-0.6250580444d-1, 0.1458381230d-1,
     & -0.2108817737d-2, 0.7411247290d-3,
     &  1.0d0          ,-0.5571175625d-2, 0.6225310236d-1,
     & -0.3137378427d-2, 0.1931496439d-2/

      data (p2(i),i=0,3),(q2(j),j=0,3)
     & /0.2868328584d+0, 0.3564363231d+0, 0.1523518695d+0,
     &  0.2251304883d-1,
     &  1.0d0          , 0.6191136137d+0, 0.1720721448d+0,
     &  0.2278594771d-1/

      data (p3(i),i=0,3),(q3(j),j=0,3)
     & /0.2868329066d+0, 0.3003828436d+0, 0.9950951941d-1,
     &  0.8733827185d-2,
     &  1.0d0          , 0.4237190502d+0, 0.1095631512d+0,
     &  0.8693851567d-2/

      data (p4(i),i=0,3),(q4(j),j=0,3)
     & /0.1000351630d+1, 0.4503592498d+1, 0.1085883880d+2,
     &  0.7536052269d+1,
     &  1.0d0          , 0.5539969678d+1, 0.1933581111d+2,
     &  0.2721321508d+2/

      data (p5(i),i=0,3),(q5(j),j=0,3)
     & /0.1000006517d+1, 0.4909414111d+2, 0.8505544753d+2,
     &  0.1532153455d+3,
     &  1.0d0          , 0.5009928881d+2, 0.1399819104d+3,
     &  0.4200002909d+3/

      data (p6(i),i=0,3),(q6(j),j=0,3)
     & /0.1000000983d+1, 0.1329868456d+3, 0.9162149244d+3,
     & -0.9605054274d+3,
     &  1.0d0          , 0.1339887843d+3, 0.1055990413d+4,
     &  0.5532224619d+3/

      data (a1(i),i=1,3)
     & /-0.4583333333d0, 0.6675347222d+0,-0.1641741416d+1/

      data (a2(i),i=1,3)
     & /1.0d0          ,-0.4227843351d+0,-0.2043403138d+1/

*-----------------------------------------------------------------------

      v=x

      if(v .lt. -5.5) then
          u=exp(v+1.0)
          glands=0.3989422803*exp(-1.0/u)*sqrt(u)*
     &     (1.0+(a1(1)+(a1(2)+a1(3)*u)*u)*u)
      else if(v .lt. -1.0) then
          u=exp(-v-1.0)
          glands=(exp(-u)/sqrt(u))*
     &     (p1(0)+(p1(1)+(p1(2)+(p1(3)+p1(4)*v)*v)*v)*v)/
     &     (q1(0)+(q1(1)+(q1(2)+(q1(3)+q1(4)*v)*v)*v)*v)
      else if(v .lt. 1.0) then
          glands=(p2(0)+(p2(1)+(p2(2)+p2(3)*v)*v)*v)/
     &     (q2(0)+(q2(1)+(q2(2)+q2(3)*v)*v)*v)
      else if(v .lt. 4.0) then
          glands=(p3(0)+(p3(1)+(p3(2)+p3(3)*v)*v)*v)/
     &     (q3(0)+(q3(1)+(q3(2)+q3(3)*v)*v)*v)
      else if(v .lt. 12.0) then
          u=1.0/v
          glands=(p4(0)+(p4(1)+(p4(2)+p4(3)*u)*u)*u)/
     &     (q4(0)+(q4(1)+(q4(2)+q4(3)*u)*u)*u)
      else if(v .lt. 50.0) then
          u=1.0/v
          glands=(p5(0)+(p5(1)+(p5(2)+p5(3)*u)*u)*u)/
     &     (q5(0)+(q5(1)+(q5(2)+q5(3)*u)*u)*u)
      else if(v .lt. 300.0) then
          u=1.0/v
          glands=(p6(0)+(p6(1)+(p6(2)+p6(3)*u)*u)*u)/
     &     (q6(0)+(q6(1)+(q6(2)+q6(3)*u)*u)*u)
      else
          u=1.0/(v-v*log(v)/(v+1.0))
          glands=1.0-(a2(1)+(a2(2)+a2(3)*u)*u)*u
      end if

      return
      end


************************************************************************
*                                                                      *
      function glande(x)
*                                                                      *
*       modified by K.Niita on 2003/10/01                              *
*                                                                      *
************************************************************************
*                                                                      *
*       copy of the cern library routine denlan (g110)                 *
*                                                                      *
*         ==>called by : gvaviv                                        *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      dimension p1(0:4),p2(0:4),p3(0:4),p4(0:4),p5(0:4),p6(0:4)
      dimension q1(0:4),q2(0:4),q3(0:4),q4(0:4),q5(0:4),q6(0:4)
      dimension a1(1:3),a2(1:2)

      data (p1(i),i=0,4),(q1(j),j=0,4)
     & /0.4259894875d+0,-0.1249762550d+0, 0.3984243700d-1,
     & -0.6298287635d-2, 0.1511162253d-2,
     &  1.0d0          ,-0.3388260629d+0, 0.9594393323d-1,
     & -0.1608042283d-1, 0.3778942063d-2/

      data (p2(i),i=0,4),(q2(j),j=0,4)
     & /0.1788541609d+0, 0.1173957403d+0, 0.1488850518d-1,
     & -0.1394989411d-2, 0.1283617211d-3,
     &  1.0d0          , 0.7428795082d+0, 0.3153932961d+0,
     &  0.6694219548d-1, 0.8790609714d-2/

      data (p3(i),i=0,4),(q3(j),j=0,4)
     & /0.1788544503d+0, 0.9359161662d-1, 0.6325387654d-2,
     &  0.6611667319d-4,-0.2031049101d-5,
     &  1.0d0          , 0.6097809921d+0, 0.2560616665d+0,
     &  0.4746722384d-1, 0.6957301675d-2/

      data (p4(i),i=0,4),(q4(j),j=0,4)
     & /0.9874054407d+0, 0.1186723273d+3, 0.8492794360d+3,
     & -0.7437792444d+3, 0.4270262186d+3,
     &  1.0d0          , 0.1068615961d+3, 0.3376496214d+3,
     &  0.2016712389d+4, 0.1597063511d+4/

      data (p5(i),i=0,4),(q5(j),j=0,4)
     & /0.1003675074d+1, 0.1675702434d+3, 0.4789711289d+4,
     &  0.2121786767d+5,-0.2232494910d+5,
     &  1.0d0          , 0.1569424537d+3, 0.3745310488d+4,
     &  0.9834698876d+4, 0.6692428357d+5/

      data (p6(i),i=0,4),(q6(j),j=0,4)
     & /0.1000827619d+1, 0.6649143136d+3, 0.6297292665d+5,
     &  0.4755546998d+6,-0.5743609109d+7,
     &  1.0d0          , 0.6514101098d+3, 0.5697473333d+5,
     &  0.1659174725d+6,-0.2815759939d+7/

      data (a1(i),i=1,3)
     & /0.4166666667d-1,-0.1996527778d-1, 0.2709538966d-1/

      data (a2(i),i=1,2)
     & /-0.1845568670d+1,-0.4284640743d+1/

*-----------------------------------------------------------------------

      v=x

      if(v .lt. -5.5) then
          u=exp(v+1.0)
          glande=0.3989422803*(exp(-1.0/u)/sqrt(u))*
     &     (1.0+(a1(1)+(a1(2)+a1(3)*u)*u)*u)
      else if(v .lt. -1.0) then
          u=exp(-v-1.0)
          glande=exp(-u)*sqrt(u)*
     &     (p1(0)+(p1(1)+(p1(2)+(p1(3)+p1(4)*v)*v)*v)*v)/
     &     (q1(0)+(q1(1)+(q1(2)+(q1(3)+q1(4)*v)*v)*v)*v)
      else if(v .lt. 1.0) then
          glande=(p2(0)+(p2(1)+(p2(2)+(p2(3)+p2(4)*v)*v)*v)*v)/
     &     (q2(0)+(q2(1)+(q2(2)+(q2(3)+q2(4)*v)*v)*v)*v)
      else if(v .lt. 5.0) then
          glande=(p3(0)+(p3(1)+(p3(2)+(p3(3)+p3(4)*v)*v)*v)*v)/
     &     (q3(0)+(q3(1)+(q3(2)+(q3(3)+q3(4)*v)*v)*v)*v)
      else if(v .lt. 12.0) then
          u=1.0/v
          glande=u**2*(p4(0)+(p4(1)+(p4(2)+(p4(3)+p4(4)*u)*u)*u)*u)/
     &     (q4(0)+(q4(1)+(q4(2)+(q4(3)+q4(4)*u)*u)*u)*u)
      else if(v .lt. 50.0) then
          u=1.0/v
          glande=u**2*(p5(0)+(p5(1)+(p5(2)+(p5(3)+p5(4)*u)*u)*u)*u)/
     &     (q5(0)+(q5(1)+(q5(2)+(q5(3)+q5(4)*u)*u)*u)*u)
      else if(v .lt. 300.0) then
          u=1.0/v
          glande=u**2*(p6(0)+(p6(1)+(p6(2)+(p6(3)+p6(4)*u)*u)*u)*u)/
     &     (q6(0)+(q6(1)+(q6(2)+(q6(3)+q6(4)*u)*u)*u)*u)
      else
          u=1.0/(v-v*log(v)/(v+1.0))
          glande=u**2*(1.0+(a2(1)+a2(2)*u)*u)
      end if

      return
      end

************************************************************************
*                                                                      *
      function getdEdxH2O(E)
*                                                                      *
*       modified by T.Sato on 2014/08/21                               *
*                                                                      *
*       Return collision stopping power (MeV/cm)                       *
*       of electron with energy E(MeV)                                 *
*       Data table are taken from ESTAR                                *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      parameter(nebin=92)
      dimension Ene(nebin),dEdx(nebin)
      save Ene,dEdx

*-----------------------------------------------------------------------
      data (Ene(i),dEdx(i),i=   1,  nebin)/
     &   1.000E-03,   1.198E+02,
     &   1.500E-03,   9.177E+01,
     &   2.000E-03,   7.521E+01,
     &   2.500E-03,   6.416E+01,
     &   3.000E-03,   5.621E+01,
     &   3.500E-03,   5.018E+01,
     &   4.000E-03,   4.543E+01,
     &   4.500E-03,   4.159E+01,
     &   5.000E-03,   3.841E+01,
     &   6.000E-03,   3.344E+01,
     &   7.000E-03,   2.971E+01,
     &   8.000E-03,   2.681E+01,
     &   9.000E-03,   2.448E+01,
     &   1.000E-02,   2.256E+01,
     &   1.500E-02,   1.647E+01,
     &   2.000E-02,   1.317E+01,
     &   2.500E-02,   1.109E+01,
     &   3.000E-02,   9.653E+00,
     &   3.500E-02,   8.592E+00,
     &   4.000E-02,   7.777E+00,
     &   4.500E-02,   7.130E+00,
     &   5.000E-02,   6.603E+00,
     &   6.000E-02,   5.797E+00,
     &   7.000E-02,   5.207E+00,
     &   8.000E-02,   4.757E+00,
     &   9.000E-02,   4.402E+00,
     &   1.000E-01,   4.115E+00,
     &   1.500E-01,   3.238E+00,
     &   2.000E-01,   2.793E+00,
     &   2.500E-01,   2.528E+00,
     &   3.000E-01,   2.355E+00,
     &   3.500E-01,   2.235E+00,
     &   4.000E-01,   2.148E+00,
     &   4.500E-01,   2.083E+00,
     &   5.000E-01,   2.034E+00,
     &   6.000E-01,   1.963E+00,
     &   7.000E-01,   1.917E+00,
     &   8.000E-01,   1.886E+00,
     &   9.000E-01,   1.864E+00,
     &   1.000E+00,   1.849E+00,
     &   1.500E+00,   1.822E+00,
     &   2.000E+00,   1.824E+00,
     &   2.500E+00,   1.834E+00,
     &   3.000E+00,   1.846E+00,
     &   3.500E+00,   1.858E+00,
     &   4.000E+00,   1.870E+00,
     &   4.500E+00,   1.882E+00,
     &   5.000E+00,   1.892E+00,
     &   6.000E+00,   1.911E+00,
     &   7.000E+00,   1.928E+00,
     &   8.000E+00,   1.943E+00,
     &   9.000E+00,   1.956E+00,
     &   1.000E+01,   1.968E+00,
     &   1.500E+01,   2.014E+00,
     &   2.000E+01,   2.046E+00,
     &   2.500E+01,   2.070E+00,
     &   3.000E+01,   2.089E+00,
     &   3.500E+01,   2.105E+00,
     &   4.000E+01,   2.118E+00,
     &   4.500E+01,   2.129E+00,
     &   5.000E+01,   2.139E+00,
     &   6.000E+01,   2.156E+00,
     &   7.000E+01,   2.170E+00,
     &   8.000E+01,   2.182E+00,
     &   9.000E+01,   2.193E+00,
     &   1.000E+02,   2.202E+00,
     &   1.500E+02,   2.238E+00,
     &   2.000E+02,   2.263E+00,
     &   2.500E+02,   2.282E+00,
     &   3.000E+02,   2.297E+00,
     &   3.500E+02,   2.311E+00,
     &   4.000E+02,   2.322E+00,
     &   4.500E+02,   2.332E+00,
     &   5.000E+02,   2.341E+00,
     &   6.000E+02,   2.357E+00,
     &   7.000E+02,   2.370E+00,
     &   8.000E+02,   2.381E+00,
     &   9.000E+02,   2.391E+00,
     &   1.000E+03,   2.400E+00,
     &   1.500E+03,   2.435E+00,
     &   2.000E+03,   2.459E+00,
     &   2.500E+03,   2.479E+00,
     &   3.000E+03,   2.494E+00,
     &   3.500E+03,   2.507E+00,
     &   4.000E+03,   2.519E+00,
     &   4.500E+03,   2.529E+00,
     &   5.000E+03,   2.538E+00,
     &   6.000E+03,   2.553E+00,
     &   7.000E+03,   2.566E+00,
     &   8.000E+03,   2.578E+00,
     &   9.000E+03,   2.588E+00,
     &   1.000E+04,   2.597E+00/
*-----------------------------------------------------------------------

      do ie=1,nebin
       if(E.le.Ene(ie)) exit
      enddo

      if(ie.eq.1) then
       getdEdxH2O=dEdx(1)
      elseif(ie.le.nebin) then
       ratio=(E-Ene(ie-1))/(Ene(ie)-Ene(ie-1))
       getdEdxH2O=dEdx(ie-1)+(dEdx(ie)-dEdx(ie-1))*ratio
      else
       getdEdxH2O=dEdx(nebin)
      endif

      return

      end



*-----------------------------------------------------------------------

