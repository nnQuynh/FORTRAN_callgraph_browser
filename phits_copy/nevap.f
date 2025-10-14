************************************************************************
*                                                                      *
      subroutine nevap(ipos)
*                                                                      *
*       main control routine of Nuclear Evaporation                    *
*       modified by K.Niita on 2005/12/22                              *
*                                                                      *
*     input:                                                           *
*                                                                      *
*        ipos  : = 0, call from ovly12 and 13 ,=1 from tally           *
*                = 2, call from sctneut                                *
*                                                                      *
*     output:                                                          *
*                                                                      *
*---- in common -------------------------------------------------------*
*                                                                      *
*        nclsts   : total number of out going particles and nuclei     *
*                                                                      *
*        iclusts(nclsts)                                               *
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
*        jclusts(i,nclsts)                                             *
*                                                                      *
*                i = 0, angular momentum                               *
*                  = 1, proton number                                  *
*                  = 2, neutron number                                 *
*                  = 3, ip, see below                                  *
*                  = 4, status of the particle 0: real, <0 : dead      *
*                  = 5, charge                                         *
*                  = 6, baryon number                                  *
*                  = 7, kf code                                        *
*                  = 8, isomer level (0: Ground, 1,2: 1st, 2nd isomer) *
*                                                                      *
*        qclusts(i,nclsts)                                             *
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
*        numpam(i) : after evaporation                                 *
*                                                                      *
*        numpal(i) : event of final = cascade + evaoparation           *
*        rumpal(i) : weight                                            *
*                  : total number of out going particles or nuclei     *
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
*                                                                      *
*                  = 11, other particles                               *
*                                                                      *
*                  = 12, electron                                      *
*                  = 13, positron                                      *
*                  = 14, photon                                        *
*                                                                      *
*                  = 15, deuteron                                      *
*                  = 16, triton                                        *
*                  = 17, 3He                                           *
*                  = 18, Alpha                                         *
*                  = 19, residual nucleus                              *
*                                                                      *
*        kdecay(4) = 0 : no fission                                    *
*                  = 1 : with fission                                  *
*                                                                      *
************************************************************************
C for USE_MOD_COUNTER
      use mod_counter, only: rncnt,rnint,rnintr,rnpnt,rnpntr
     &                      ,aevts,aevtr,bevts,bevtr
     &                      ,iaevt, ibevt, jaevt, jbevt
!$   &                      ,rncnt2,rnint2,rnintr2,rnpnt2,rnpntr2
!$   &                      ,aevts2,aevtr2,bevts2,bevtr2

      use levdat, only : nph, eph, kfejec
      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param00.inc'
      include 'param02.inc'
      include 'err.inc'
      include 'param.inc'

*-----------------------------------------------------------------------

      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
!$OMP THREADPRIVATE(/icomon/)
      common /qparm/  ielas,icasc,iqstep,lvlopt,igamma
      common /bparm/  andt,jevap,npidk

*-----------------------------------------------------------------------

      common /clustf/ nclst, iclust(nnn)
!$OMP THREADPRIVATE(/clustf/)
      common /clustg/ jclust(0:8,nnn),  qclust(0:12,nnn)
!$OMP THREADPRIVATE(/clustg/)
      common /clustp/ rumpat(0:20), numpat(0:20)
!$OMP THREADPRIVATE(/clustp/)

*-----------------------------------------------------------------------

      common /clusts/ nclust, kclust(3,nnn)
!$OMP THREADPRIVATE(/clusts/)
      common /clustu/ lclust(0:8,nnn), sclust(0:12,nnn)
!$OMP THREADPRIVATE(/clustu/)
      common /clustm/ rumpam(0:20), numpam(0:20)
!$OMP THREADPRIVATE(/clustm/)
      common /clustv/ kdecay(4)
!$OMP THREADPRIVATE(/clustv/)
      common /mathzn/ mathz, mathn, jcoll, kcoll
!$OMP THREADPRIVATE(/mathzn/)

*-----------------------------------------------------------------------

      common /clustt/ nclsts, iclusts(nnn)
!$OMP THREADPRIVATE(/clustt/)
      common /clustw/ jclusts(0:8,nnn),  qclusts(0:12,nnn)
!$OMP THREADPRIVATE(/clustw/)
      common /clustl/ rumpal(0:20), numpal(0:20)
!$OMP THREADPRIVATE(/clustl/)
      common /emode1/ nevhin, nlowlv, nefiss, ntwidt, mtprec,ipcnt(70)
!$OMP THREADPRIVATE(/emode1/)
      common /nrfmem/ spis, spgr, levabs, lflgnrf, mpole
!$OMP THREADPRIVATE(/nrfmem/)

*-----------------------------------------------------------------------


      common /wtsave/ oldwt
!$OMP THREADPRIVATE(/wtsave/)

*-----------------------------------------------------------------------

      common /muflag/ imuinthit,imubrmhit,imuppdhit,imucaphit,imucapflag
!$OMP THREADPRIVATE(/muflag/)

*-----------------------------------------------------------------------

      dimension numpal0(20), rumpal0(20)

*-----------------------------------------------------------------------
      common /jcomon/ nabov,nobch,nocas,nomax
!$OMP THREADPRIVATE(/jcomon/)
      common /tlgeom/ iblz1,iblz2
!$OMP THREADPRIVATE(/tlgeom/)
      common /regdc/  idrg(kvlmax), idgr(kvmmax)
      common /kmat1d/ idmn(0:kvlmax), idnm(kvmmax)
      common /kmatad/ matadd


*-----------------------------------------------------------------------
*        do loop for the clusters from CASCAD
*-----------------------------------------------------------------------

               if( ipos .eq. 2 ) then

               else if( imuinthit.eq.1 .or. imucapflag.gt.0 ) then

                  do i = 1, 20
                     numpal0(i) = numpal(i)
                     rumpal0(i) = rumpal(i)
                  end do

               else

                  nclsts = 0

               end if
*-----------------------------------------------------------------------

                  numnut = 0

                  kcoll = 0

               do i = 1, 4

                  kdecay(i) = 0

               end do

               do i = 0, 20

                     numpam(i) = 0
                     rumpam(i) = 0.0d0

                  if( ipos .ne. 2 ) then

                     numpal(i) = 0
                     rumpal(i) = 0.0d0

                  end if

               end do

*-----------------------------------------------------------------------

      if( nclst .gt. 0 ) then

      do 100 i = 1, nclst

                  ex = qclust(6,i)

*-----------------------------------------------------------------------
*        statistical decay of cluster
*-----------------------------------------------------------------------

         if( iclust(i) .eq. 0 .and. jevap .gt. 0 ) then

*-----------------------------------------------------------------------

                  iz = jclust(1,i)
                  in = jclust(2,i)

               if( iz .lt. 0 .or. in .lt. 0 .or.
     &             iz + in .eq. 0 ) goto 100

                  if( ipos .eq. 0 .or. ipos .eq. 2 ) call cputime(4)

                  jj = jclust(0,i)
                  bi = qclust(0,i)

                  px = qclust(1,i)
                  py = qclust(2,i)
                  pz = qclust(3,i)
                  et = qclust(4,i)
                  rm = qclust(5,i)
                  ex = qclust(6,i)
                  ek = qclust(7,i)
                  wt = qclust(8,i)

                  pt = sqrt( px**2 + py**2 + pz**2 )

                  ex0 = ex
                  er0 = ek

*-----------------------------------------------------------------------
*           statistical decay, call sdmexec
*-----------------------------------------------------------------------

            if( iqstep .eq. 4 ) then

               call sdmexec(iz,in,jj,ex,px,py,pz,wt,ierr)


            else if( iqstep .eq. 2 .or. iqstep .eq. 3 .or.
     &               iqstep .eq. 5 ) then

               call erupin(iz,in,ex,px,py,pz,pt,et,rm,wt)


            else if( iqstep .eq. 6 ) then

               call gemexec(iz,in,ex,jj,px,py,pz,pt,et,rm,wt,ierr,ipos)
               mtprec = 0

            end if

*-----------------------------------------------------------------------
*           new booking
*-----------------------------------------------------------------------

            do j = 1, nclust

             if(nclust .gt. 1) lflgnrf = 0 ! 2018/10/9 Ogawa. Forget NRF if evaporated. Nuclear structure (level scheme) changes

               if( kclust(1,j) .ge. 100 ) then

                     nclsts = nclsts + 1

                  if( lclust(1,j) .eq. 0 .and.
     &                lclust(2,j) .eq. 0 ) then

                     kf    = 22
                     ibary = 0
                     ipid  = 4
                     ippad = 14

                  else if( lclust(1,j) .eq. 1 .and.
     &                     lclust(2,j) .eq. 0 ) then

                     kf    = 2212
                     ibary = 1
                     ipid  = 1
                     ippad = 1

                  else if( lclust(1,j) .eq. 0 .and.
     &                     lclust(2,j) .eq. 1 ) then

                     kf    = 2112
                     ibary = 1
                     ipid  = 2
                     ippad = 2

                  else

                     ipid  = 0
                     ibary = lclust(1,j) + lclust(2,j)

                     kf    = lclust(1,j) * 1000000
     &                     + lclust(1,j) + lclust(2,j)

                     if( lclust(1,j) .eq. 1 .and.
     &                   lclust(2,j) .eq. 1 ) then

                        ippad  = 15

                     else if( lclust(1,j) .eq. 1 .and.
     &                        lclust(2,j) .eq. 2 ) then

                        ippad  = 16

                     else if( lclust(1,j) .eq. 2 .and.
     &                        lclust(2,j) .eq. 1 ) then

                        ippad  = 17

                     else if( lclust(1,j) .eq. 2 .and.
     &                        lclust(2,j) .eq. 2 ) then

                        ippad  = 18

                     else

                        ippad  = 19

                     end if

                  end if

                     iclusts(nclsts)   = ipid

                     jclusts(0,nclsts) = lclust(0,j)
                     jclusts(1,nclsts) = lclust(1,j)
                     jclusts(2,nclsts) = lclust(2,j)
                     jclusts(3,nclsts) = ippad
                     jclusts(4,nclsts) = 0
                     jclusts(5,nclsts) = lclust(1,j)
                     jclusts(6,nclsts) = ibary
                     jclusts(7,nclsts) = kf
                     jclusts(8,nclsts) = 0

                     qclusts(0,nclsts) = bi


                  do l = 1, 12

                     qclusts(l,nclsts) = sclust(l,j)

                  end do

                     if( qclusts(6,nclsts) .lt. 0.0 )
     &                   qclusts(6,nclsts) = 0.0

                     numpam(ippad) = numpam(ippad) + 1
                     rumpam(ippad) = rumpam(ippad) + sclust(8,j)

               end if

            end do

*-----------------------------------------------------------------------

               if( ipos .eq. 0 .or. ipos .eq. 2 ) then
                  rncnt(4) = rncnt(4) + 1.0
                  call cputime(4)
               end if

*-----------------------------------------------------------------------
*        particle or without evaporation : pass through
*-----------------------------------------------------------------------

         else

*-----------------------------------------------------------------------

               if( iclust(i) .eq. 0 ) then

                     if( jclust(1,i) .eq. 1 .and.
     &                   jclust(2,i) .eq. 1 ) then

                        ippad  = 15

                     else if( jclust(1,i) .eq. 1 .and.
     &                        jclust(2,i) .eq. 2 ) then

                        ippad  = 16

                     else if( jclust(1,i) .eq. 2 .and.
     &                        jclust(2,i) .eq. 1 ) then

                        ippad  = 17

                     else if( jclust(1,i) .eq. 2 .and.
     &                        jclust(2,i) .eq. 2 ) then

                        ippad  = 18

                     else

                        ippad  = 19

                     end if

                        jclust(3,i) = ippad

               end if

*-----------------------------------------------------------------------

                     nclsts = nclsts + 1

                     iclusts(nclsts)   = iclust(i)

                  do j = 0, 8

                     jclusts(j,nclsts) = jclust(j,i)

                  end do

                  do j = 0, 12

                     qclusts(j,nclsts) = qclust(j,i)

                  end do

                     ippad = jclust(3,i)

                     numpam(ippad) = numpam(ippad) + 1
                     rumpam(ippad) = rumpam(ippad) + qclust(8,i)

         end if

*-----------------------------------------------------------------------

  100 continue

      end if

      jclusts(8,1:nclsts) = 0   ! nullify isomer level just in case

*-----------------------------------------------------------------------
*        gamma production from deexcited residual nuclei
*-----------------------------------------------------------------------

      if( igamma .ne. 0 .and. nclsts .gt. 0 ) then

         do i = 1, nclsts

            if( iclusts(i) .eq. 0 .and.
     &          qclusts(6,i) .gt. 0.0 .and.
     &          jclusts(1,i) .gt. 0 .and.
     &          jclusts(2,i) .gt. 0 ) then

                     iz  = jclusts(1,i)
                     in  = jclusts(2,i)

                     pxn = qclusts(1,i)
                     pyn = qclusts(2,i)
                     pzn = qclusts(3,i)
                     emn = qclusts(5,i) + qclusts(6,i) * 1.d-3   !  2014/8/28 ogawa. Rest mass at excited state
                     ern = sqrt(pxn**2 + pyn**2 + pzn**2 + emn **2)
                     exn = qclusts(6,i)
                     ekn = qclusts(7,i)
                     wt  = qclusts(8,i)

                     bex = pxn / sqrt(ern**2 + (exn * 1.d-3)**2) / 2.d0  ! devided by true total energy
                     bey = pyn / sqrt(ern**2 + (exn * 1.d-3)**2) / 2.d0
                     bez = pzn / sqrt(ern**2 + (exn * 1.d-3)**2) / 2.d0
                     gam = ern / emn

                  if( ipos .eq. 0 .or. ipos .eq. 2 ) call cputime(5)
                     apr = dble( iz + in )
                     zpr = dble( iz )
                     exi = exn
                     spin = dble(jclusts(0,i)) ! 2020/3/26 spin is half integer. correct upstream later

                  call dexgam(apr,zpr,exi,spin,isomle)   ! 2013/4/17 Ogawa
                  jclusts(8,i) = isomle  ! record isomeric identity

                  if( ipos .eq. 0 .or. ipos .eq. 2 ) then
                     rncnt(5) = rncnt(5) + 1.0
                     call cputime(5)
                  end if

               if( nph .gt. 0 ) then

                     pxg = 0.d0
                     pyg = 0.d0
                     pzg = 0.d0

                  do j = 1, nph

                     pr   = eph(j) / 1000.0

                     if(lflgnrf .eq. 1) then ! 2014/9/4 ogawa. NRF photons are anisotropic
                      beta = sqrt(bex**2 + bey**2 + bez**2)
                      call NRF_ang(axf, ayf, azf, iz, in)
                      pxr = pr * axf
                      pyr = pr * ayf
                      pzr = pr * azf
                     else
                      cos1 = 1.0 - 2.0 * rn(0)
                      sin1 = sqrt( 1.0 - cos1**2 )
                      phi1 = 2.0 * pi * rn(0)
                      pxr = pr * sin1 * cos(phi1)
                      pyr = pr * sin1 * sin(phi1)
                      pzr = pr * cos1
                     endif

                     if(igamma .ge. 1) then
                      pxn = pxn - pxr
                      pyn = pyn - pyr
                      pzn = pzn - pzr

                      pcs = pxn * pxr + pyn * pyr + pzn * pzr

                      trans = ( pcs / ( ern + emn ) + pr ) / emn

                      pxl = pxr + trans * pxn
                      pyl = pyr + trans * pyn
                      pzl = pzr + trans * pzn

                      exn = (qclusts(6,i)-sum(eph(1:j))) * 1.d-3
                      emn = qclusts(5,i) + exn
                      ern = sqrt(pxn**2 + pyn**2 + pzn**2 + emn **2)
                      ekns= ekn
                      ekn = (ern - emn)*1.d3

                      egl = sqrt(pxr**2 + pyr**2 + pzr**2) +(ekns - ekn)
     &                * 1.d-3

                     else ! No doppler mode
                      pxl = pxr
                      pyl = pyr
                      pzl = pzr
                      egl = pr
                     endif

                     nclsts = nclsts + 1

                     if( kfejec(j) .eq. 22 ) then     ! gamma emmision
                         iclusts(nclsts)    = 4
                         jclusts(3,nclsts)  = 14
                     elseif( kfejec(j) .eq. 11 ) then ! internal conversion electron
                         iclusts(nclsts)    = 7
                         jclusts(3,nclsts)  = 12
                     else
                         write(ErrCha,*) 'warning: particle species is
     & strange in electro-magnetic deexcitation in nevap. kf=',kfejec(j)
                         ErrID = 'L:575/R:nevap/F:nevap.f'
                         call ErrWrite(ErrID,ErrCha)
c                         stop  ! 2022/12/26 Ogawa. Comment out and following emergency recovery
                         kfejec(j)          = 22
                         iclusts(nclsts)    = 4
                         jclusts(3,nclsts)  = 14
                     endif

                     jclusts(0,nclsts)  = 0
                     jclusts(1,nclsts)  = 0
                     jclusts(2,nclsts)  = 0
                     jclusts(4,nclsts)  = 0
                     jclusts(5,nclsts)  = 0
                     jclusts(6,nclsts)  = 0
                     jclusts(7,nclsts)  = kfejec(j) ! 22(photon) or 11(electron)
                     jclusts(8,nclsts)  = 0

                     qclusts(1,nclsts)  = pxl
                     qclusts(2,nclsts)  = pyl
                     qclusts(3,nclsts)  = pzl
                     qclusts(4,nclsts)  = egl
                     qclusts(5,nclsts)  = 0.0d0
                     qclusts(6,nclsts)  = 0.0d0
                     qclusts(7,nclsts)  = egl * 1000.
                     qclusts(8,nclsts)  = wt
                     qclusts(9,nclsts)  = 0.0
                     qclusts(10,nclsts) = 0.0d0
                     qclusts(11,nclsts) = 0.0d0
                     qclusts(12,nclsts) = 0.0d0

                     numpam(14) = numpam(14) + 1
                     rumpam(14) = rumpam(14) + wt

                  end do

                     lflgnrf = 0 ! 2014/9/4 ogawa.  Clear NRF flag
                     mpole = 0

                     qclusts(1,i) = pxn
                     qclusts(2,i) = pyn
                     qclusts(3,i) = pzn
                     qclusts(4,i) = ern
                     qclusts(5,i) = emn
                     qclusts(6,i) = max( 0.0d0, exn )*1.d3
                     qclusts(7,i) = (ern - emn) * 1000.d0

               end if

            end if

         end do

      end if

*-----------------------------------------------------------------------
*     total number of particles and nuclei
*-----------------------------------------------------------------------

            if( ipos .ne. 2 ) then

               do i = 0, 19

                  numpal(i) = numpam(i)
                  rumpal(i) = rumpam(i) * oldwt

               end do

               if( numpal(ityp) .gt. 0 .and.
     &           ( jcoll .eq. 8 .or. jcoll .eq. 11 ) ) then

                  numpal(ityp) = numpal(ityp) - 1
                  rumpal(ityp) = rumpal(ityp) - oldwt

               end if

               if( kdecay(4) .gt. 0 ) then

                  numpal(20) = 1
                  rumpal(20) = oldwt

               end if

               if(imuinthit.eq.1 .or.
     &            (imucaphit.eq.1 .and. imucapflag.eq.1)) then
                do i = 1, 20
                 numpal(i) = numpal(i) + numpal0(i)
                 rumpal(i) = rumpal(i) + rumpal0(i)
                enddo
               endif

            end if

*-----------------------------------------------------------------------
*        collision type check, apsorption or fission
*-----------------------------------------------------------------------

            if( jcoll .ne. 6 .and. jcoll .ne. 10 ) then

               if( nclst .eq. 1 .and.
     &             kdecay(4) .eq. 0 .and.
     &             nabov .eq. 0 ) then

                     kcoll = 2

               elseif( kdecay(4) .gt. 0 ) then

                     kcoll = 1

               endif

            endif

*-----------------------------------------------------------------------

            if( kcoll .eq. 1 ) then

               ireg = idgr( iblz1 )

               if( mat .le. kvlmax .and. mat .gt. 0 ) then
                  if( matadd .ne. 0 ) then
                     imat = mat
                  else
                     imat = idnm( idmn(mat) )
                  endif
               end if

                  aevts(iaevt+25,ireg) =
     &            aevts(iaevt+25,ireg) + oldwt
               if( imat .gt. 0 ) then
                  bevts(ibevt+25,imat) =
     &            bevts(ibevt+25,imat) + oldwt
               endif

            endif

*-----------------------------------------------------------------------

      return
      end

