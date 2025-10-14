************************************************************************
*                                                                      *
*        PART 6: Statistical Decay Model                               *
*                                                                      *
*                                                                      *
*   List of subprograms in rough order of relevance with main purpose  *
*      ( s = subroutine, f = function, b = block data, e = entry )     *
*                                                                      *
*                                                                      *
*  s  sdmtest   a short main for test of SDM                           *
*  s  sdment    to entry of SDM from QMD results                       *
*  s  sdmint    to initialize SDM                                      *
*  s  sdmexec   to execute statistical particle decay / fission.       *
*  s  sdmwid0   to calculate decay width and decay products            *
*               without angular momentum.                              *
*  s  sdmwid1   to calculate decay width and decay products            *
*               with angular momentum.                                 *
*  f  sdmlev    to calculate level density                             *
*  s  sdjsum    to sum up angular momentum.                            *
*  s  sdmfisw   to calculate fission width                             *
*  s  sdmfiss   to determine the mass and charge of fission fragment   *
*  f  sdmfisb   to determine the fission barrier.                      *
*  f  bndeng    to give binding energy per baryon (MeV)                *
*  f  eliq      to calculate liquid drop binding energy (MeV)          *
*                                                                      *
************************************************************************
*                                                                      *
************************************************************************
*                                                                      *
      subroutine sdmset
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              initialize SDM for PHITS                                *
*                                                                      *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------

      implicit double precision(a-h, o-z)

*-----------------------------------------------------------------------

      include 'param00.inc'
      include 'param01.inc'

*-----------------------------------------------------------------------

      common /input1/ mstq1(mxpa1), parq1(mxpa1)
!$OMP THREADPRIVATE(/input1/)
      common /vriab3/ qmdfac, sdmfac
!$OMP THREADPRIVATE(/vriab3/)

*-----------------------------------------------------------------------

               mstq1(121) = 0

               mstq1(120) = 1

               mstq1(122) = 1
               mstq1(123) = 1
               mstq1(124) = 1

               parq1(120) = 0.1

               mstq1(125) = 30
               mstq1(126) = 1
               mstq1(127) = 1
               mstq1(128) = 1
               mstq1(129) = 1

               qmdfac = 1.0

               call sdmint

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      subroutine sdmint
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 27                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to initialize SDM                                       *
*                                                                      *
*----------------------------------------------------------------------*
*     switch for SDM                                                   *
*----------------------------------------------------------------------*
*                                                                      *
*        issdm      : mstq1(120)                                       *
*                                                                      *
*----------------------------------------------------------------------*
*     minimum energy cut off energy for evapolation and fission        *
*----------------------------------------------------------------------*
*                                                                      *
*        sdmemin    : parq1(161)                                       *
*                     ( D = 1.0 MeV )                                  *
*                                                                      *
*----------------------------------------------------------------------*
*     switches for decay mode                                          *
*----------------------------------------------------------------------*
*                                                                      *
*        iswids     : mstq1(121)                                       *
*               = 0 : simple decay width without gamma nor angular mom.*
*                 1 : decay width with gammma and angular momentum.    *
*                                                                      *
*        isevap     : mstq1(122)                                       *
*               = 0 : without particle evapolation                     *
*                 1 : with particle evapolation                        *
*                                                                      *
*        isfiss     : mstq1(123)                                       *
*               = 0 : without fission decay                            *
*                 1 : with fission decay                               *
*                                                                      *
*        isgrnd     : mstq1(124)                                       *
*               = 0 : without ground state decay                       *
*                 1 : with ground state decay                          *
*                                                                      *
*----------------------------------------------------------------------*
*     parameters for sdmwid1                                           *
*----------------------------------------------------------------------*
*                                                                      *
*        imengb     : mstq1(125)                                       *
*                   : (D=30) number of energy bin                      *
*                                                                      *
*        imbarr     : mstq1(126)                                       *
*               = 0 : simple barrier                                   *
*                 1 : modified barrier of sepc                         *
*                                                                      *
*        imangm     : mstq1(127)                                       *
*               = 0 : without angular momentum                         *
*                 1 : with angular momentum                            *
*                                                                      *
*        imlevd     : mstq1(128)                                       *
*               = 0 : level density without angular momentum           *
*                 1 : level density with angular momentum              *
*                                                                      *
*        imgamm     : mstq1(129)                                       *
*               = 0 : without gamma decay                              *
*                 1 : with gamma decay                                 *
*                                                                      *
************************************************************************

      implicit double precision(a-h, o-z)
      include 'param00.inc'
      include 'param01.inc'

*-----------------------------------------------------------------------

      common /input1/ mstq1(mxpa1), parq1(mxpa1)
!$OMP THREADPRIVATE(/input1/)

      common /vriab3/ qmdfac, sdmfac
!$OMP THREADPRIVATE(/vriab3/)

      common /sdmcut/ sdmemin
!$OMP THREADPRIVATE(/sdmcut/)

      common /sdmsw0/ issdm
!$OMP THREADPRIVATE(/sdmsw0/)
      common /sdmsw1/ iswids, isevap, isfiss, isgrnd
!$OMP THREADPRIVATE(/sdmsw1/)
      common /sdmsw2/ imengb, imbarr, imangm, imlevd, imgamm
!$OMP THREADPRIVATE(/sdmsw2/)

      common /clustx/ tdecay(4)
!$OMP THREADPRIVATE(/clustx/)

*-----------------------------------------------------------------------
*        switches for decay modes
*-----------------------------------------------------------------------

            issdm  = mstq1(120)

               sdmfac = 1.0

               if( issdm .gt. 0 ) sdmfac = 1.0 / dble(issdm)

            iswids = mstq1(121)
            isevap = mstq1(122)
            isfiss = mstq1(123)
            isgrnd = mstq1(124)

*-----------------------------------------------------------------------
*        minimum energy cut off energy for evapolation and fission
*-----------------------------------------------------------------------

            sdmemin = parq1(120)

*-----------------------------------------------------------------------
*        parameters for sdmwid1
*-----------------------------------------------------------------------

            imengb = mstq1(125)
            imbarr = mstq1(126)
            imangm = mstq1(127)
            imlevd = mstq1(128)
            imgamm = mstq1(129)

            if( imlevd .eq. 0 ) imangm = 0

*-----------------------------------------------------------------------
*        total decay mode counter
*-----------------------------------------------------------------------

            do i = 1, 4

               tdecay(i) = 0

            end do

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine sdmexec(iz,in,jang,ex,px,py,pz,wt,ierr)
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 27                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to execute statistical particle decay / fission.        *
*                                                                      *
*                                                                      *
*        Variables:                                                    *
*                                                                      *
*              iz          : proton number                             *
*              in          : neutron number                            *
*              jang        : angular momentum                          *
*              px,py,pz    : momentum in GeV                           *
*              wt          : weight change                             *
*              ierr        : error flag                                *
*                                                                      *
*                                                                      *
************************************************************************
      use NGSDATAMOD, only : bindeg, ngsdcy
      implicit double precision(a-h, o-z)
      include 'param00.inc'
      include 'param02.inc'

      include 'err.inc'

*-----------------------------------------------------------------------

      common /vriab3/ qmdfac, sdmfac
!$OMP THREADPRIVATE(/vriab3/)
      common /sdmsw0/ issdm
!$OMP THREADPRIVATE(/sdmsw0/)

      common /clusts/ nclust, kclust(3,nnn)
!$OMP THREADPRIVATE(/clusts/)
      common /clustu/ lclust(0:8,nnn), sclust(0:12,nnn)
!$OMP THREADPRIVATE(/clustu/)
      common /clustv/ kdecay(4)
!$OMP THREADPRIVATE(/clustv/)

      common /evapol/ ievmax, izip(0:100), inip(0:100), jip(0:100)

*-----------------------------------------------------------------------

      common /masmx0/ imax0

      common /sdmcut/ sdmemin
!$OMP THREADPRIVATE(/sdmcut/)

      common /sdmsw1/ iswids, isevap, isfiss, isgrnd
!$OMP THREADPRIVATE(/sdmsw1/)
      common /sdmsw2/ imengb, imbarr, imangm, imlevd, imgamm
!$OMP THREADPRIVATE(/sdmsw2/)

*-----------------------------------------------------------------------

      dimension pwid(0:100)
      dimension pcm(5)

      data initsd /0/
      save initsd !FURUTA
!$OMP THREADPRIVATE(initsd)
*-----------------------------------------------------------------------
*     initialization
*-----------------------------------------------------------------------

         if( initsd .eq. 0 ) then

            initsd = initsd + 1

            call sdmset

         end if

*-----------------------------------------------------------------------
*        Error flag
*-----------------------------------------------------------------------

            ierr = 0

*-----------------------------------------------------------------------
*        Decay mode counter
*-----------------------------------------------------------------------

            do i = 1, 4

               kdecay(i) = 0

            end do

*-----------------------------------------------------------------------
*        check of initial nucleus
*-----------------------------------------------------------------------

            if( iz .gt. maxpt .or. in .gt. maxnt ) then

               ErrCha = ''
               ErrID = 'L:326/R:sdmexec/F:sdml.f' !E81_001_001
               call ErrWrite(ErrID,ErrCha)

               write(*,*) ' **** Error at sdmexec, too large iz, in'
               write(*,*) ' ======================================='
               write(*,*) ' iz    = ', iz,    ' in    = ', in
               write(*,*) ' maxpt = ', maxpt, ' maxnt = ', maxnt

               ierr = 1
               call parastop( 999 )

            end if

*-----------------------------------------------------------------------
*        Booking.  Store initial nucleus in K, L and S vectors.
*-----------------------------------------------------------------------

            nclust = 1

            kclust(1,1) = 0
            kclust(2,1) = 0
            kclust(3,1) = 0

            lclust(0,1) = jang

            lclust(1,1) = iz
            lclust(2,1) = in

            sclust(1,1) = px
            sclust(2,1) = py
            sclust(3,1) = pz

            be = bindeg(iz,in)

            pm = ( iz * rpmass + in * rnmass - be + ex ) / 1000.0

            sclust(4,1)  = sqrt( pm**2 + px**2 + py**2 + pz**2 )
            sclust(5,1)  = pm

            sclust(6,1)  = ex
            sclust(7,1)  = ( sclust(4,1) - sclust(5,1) ) * 1000.
            sclust(8,1)  = wt
            sclust(9,1)  = 0.0
            sclust(10,1) = 0.0d0
            sclust(11,1) = 0.0d0
            sclust(12,1) = 0.0d0

*-----------------------------------------------------------------------
*        minimum energy fot particle decay / fission
*-----------------------------------------------------------------------

               emin = sdmemin

*-----------------------------------------------------------------------
*     Start decay.
*-----------------------------------------------------------------------

               num    = 1
               inum   = 0

*-----------------------------------------------------------------------
 2000    continue
*-----------------------------------------------------------------------

               inum = inum + 1
               i1   = inum

*-----------------------------------------------------------------------
*           check dimension
*-----------------------------------------------------------------------

            if( num + 2 .gt. nnn ) then

               ErrCha = ''
               ErrID = 'L:400/R:sdmexec/F:sdml.f' !E81_001_002
               call ErrWrite(ErrID,ErrCha)

               write(*,*) ' **** Error at sdmexec, too many products'
               write(*,*) ' ========================================='
               write(*,*) ' num + 2 = ', num + 2, ' nnn  = ', nnn

               ierr = 1
               call parastop( 999 )

            end if

*-----------------------------------------------------------------------
*        mother id
*-----------------------------------------------------------------------

               iz1 = lclust(1,i1)
               in1 = lclust(2,i1)
                m1 = iz1 + in1

               e1     = sclust(6,i1)

               em1 = iz1 * rpmass
     &             + in1 * rnmass
     &             - bindeg(iz1,in1)
     &             + e1

               j1     = lclust(0,i1)

               pcm(1) = sclust(1,i1)
               pcm(2) = sclust(2,i1)
               pcm(3) = sclust(3,i1)
               pcm(4) = sclust(4,i1)
               pcm(5) = sclust(5,i1)

*-----------------------------------------------------------------------
*           Gamma case
*-----------------------------------------------------------------------

            if( m1 .eq. 0 ) then

                  kclust(1,i1) = kclust(1,i1) + 100

                  goto 3000

            end if

*-----------------------------------------------------------------------
*        Reset fission and evaporation width.
*-----------------------------------------------------------------------

               tcs  = 0.0
               tcsf = 0.0

*-----------------------------------------------------------------------
*     Decay only if minimam excitation energy has.
*-----------------------------------------------------------------------

      if( e1 .gt. emin ) then

*-----------------------------------------------------------------------
*        particle evapolation width
*-----------------------------------------------------------------------

          if( isevap .eq. 1 ) then

               if( iswids .eq. 0 ) then

                  call sdmwid0(iz1,in1,e1,j1,
     &                         iz2,in2,e2,j2,
     &                         iz3,in3,e3,j3,
     &                         tcs,pwid,0,ipdec)

               else if( iswids .eq. 1 ) then

                  call sdmwid1(iz1,in1,e1,j1,
     &                         iz2,in2,e2,j2,
     &                         iz3,in3,e3,j3,
     &                         tcs,pwid,0,ipdec)

               end if

                  tcn = pwid(1)

          else

               do i = 0, ievmax

                  pwid(i) = 0.0

               end do

                  tcn = 1.0

          end if

*-----------------------------------------------------------------------
*        fission width,
*-----------------------------------------------------------------------

         if( isfiss .eq. 1 ) then

                  call sdmfisw(iz1,in1,e1,j1,tcn,tcsf,barf,saf)

          end if

*-----------------------------------------------------------------------

      end if

*-----------------------------------------------------------------------
*        Total decay width and branches
*-----------------------------------------------------------------------

               totwid = tcs + tcsf

               if( totwid .gt. 0.0 ) then

                  jdec = 1

               else if( ngsdcy(iz1,in1+iz1) .ne. 0 .and.
     &                  isgrnd .eq. 1 ) then

                  jdec = 2

               else

                  jdec = 0

               end if

*-----------------------------------------------------------------------
*     Decay channel
*-----------------------------------------------------------------------


      if( jdec .eq. 1 .or. jdec .eq. 2) then

               ifiss = 0

*-----------------------------------------------------------------------
*        evaporation or fission
*-----------------------------------------------------------------------

         if( jdec .eq. 1 ) then

*-----------------------------------------------------------------------
*           decide evaporation or fission
*-----------------------------------------------------------------------

               xran  = rn(0) * totwid

*-----------------------------------------------------------------------
*           fission channnel
*-----------------------------------------------------------------------

            if( tcsf .ge. xran ) then

               call sdmfiss(iz1,in1,e1,j1,
     &                      iz2,in2,e2,j2,
     &                      iz3,in3,e3,j3,
     &                      barf,saf,ierr)

               if( ierr .ne. 0 ) return

               ifiss = 1

               fisfac = qmdfac * sdmfac


            else

*-----------------------------------------------------------------------
*           Evaporation channels, find decay branch.
*-----------------------------------------------------------------------

                  tcip0 = tcsf

               do ipdec = 0, ievmax

                  tcip0 = tcip0 + pwid(ipdec)

                  if( tcip0 .ge. xran ) goto 101

               end do

*-----------------------------------------------------------------------
*                 Decay branch could not be found.
*                 Choice maxmum width mode
*-----------------------------------------------------------------------

                        ipdmx = 0
                        tcmx  = 0

                  do ipdec = 0, ievmax

                     if( tcmx .lt. pwid(ipdec) ) then

                        ipdmx = ipdec
                        tcmx  = pwid(ipdec)

                     end if

                  end do

                        ipdec = ipdmx

*-----------------------------------------------------------------------

  101          continue

*-----------------------------------------------------------------------
*              Determine final nuclei.
*-----------------------------------------------------------------------

               if( iswids .eq. 0 ) then

                  call sdmwid0(iz1,in1,e1,j1,
     &                         iz2,in2,e2,j2,
     &                         iz3,in3,e3,j3,
     &                         tcs,pwid,1,ipdec)

               else if( iswids .eq. 1 ) then

                  call sdmwid1(iz1,in1,e1,j1,
     &                         iz2,in2,e2,j2,
     &                         iz3,in3,e3,j3,
     &                         tcs,pwid,1,ipdec)


               end if

*-----------------------------------------------------------------------

            end if

*-----------------------------------------------------------------------
*           Mass of produced particle in MeV.
*-----------------------------------------------------------------------

            if( iz3 + in3 .gt. 0 ) then

               m2  = iz2 + in2
               m3  = iz3 + in3

               em2 = iz2 * rpmass
     &             + in2 * rnmass
     &             - bindeg(iz2,in2)
     &             + e2

               em3 = iz3 * rpmass
     &             + in3 * rnmass
     &             - bindeg(iz3,in3)
     &             + e3

            else

               m2  = m1
               m3  = 0

               em2 = iz2 * rpmass
     &             + in2 * rnmass
     &             - bindeg(iz2,in2)
     &             + e2

               em3 = 0.0

            end if

*-----------------------------------------------------------------------
*        ground state decay
*-----------------------------------------------------------------------

         else if( jdec .eq. 2 ) then

*-----------------------------------------------------------------------

               iz2 = ngsdcy(iz1,in1+iz1)/1000
               in2 = ngsdcy(iz1,in1+iz1) - iz2*1000
                m2 = iz2 + in2
                j2 = 0
                e2 = 0

               iz3 = iz1 - iz2
               in3 = in1 - in2
                m3 = iz3 + in3
                j3 = 0
                e3 = 0

               em2 = iz2 * rpmass
     &             + in2 * rnmass
     &             - bindeg(iz2,in2)
     &             + e2

               em3 = iz3 * rpmass
     &             + in3 * rnmass
     &             - bindeg(iz3,in3)
     &             + e3

*-----------------------------------------------------------------------

         end if

*-----------------------------------------------------------------------
*        line number for two daughters
*-----------------------------------------------------------------------

               i2 = num + 1
               i3 = num + 2

               num    = num + 2
               nclust = nclust + 2

*-----------------------------------------------------------------------
*        Determine momenta of decay products.
*        Lorentz-transformation into reference frame.
*-----------------------------------------------------------------------

               em1  = em1 / 1000.
               em2  = em2 / 1000.
               em3  = em3 / 1000.

               pr   = pcmsr( pcm(5), em2, em3 )

               cos1 = 1.0 - 2.0 * rn(0)
               sin1 = sqrt( 1.0 - cos1**2 )
               phi1 = 2.0 * pi * rn(0)

               pxr  = pr * sin1 * cos(phi1)
               pyr  = pr * sin1 * sin(phi1)
               pzr  = pr * cos1

               pcs    = pcm(1) * pxr + pcm(2) * pyr + pcm(3) * pzr

               ecm1   = sqrt( em2**2 + pxr**2 + pyr**2+ pzr**2 )
               trans1 = ( pcs / ( pcm(4) + pcm(5) ) + ecm1 ) / pcm(5)

               ecm2   = sqrt( em3**2 + pxr**2 + pyr**2 + pzr**2 )
               trans2 = ( -pcs / ( pcm(4) + pcm(5) ) + ecm2 ) / pcm(5)

*-----------------------------------------------------------------------
*           counting decay mode and booking
*-----------------------------------------------------------------------

            if( ifiss .eq. 0 ) then

               if( jdec .eq. 1 ) then

                  if( m3 .ne. 0 ) then

                     kdecay(1) = kdecay(1) + 1

                     kclust(1,i3) = kclust(1,i1) / 10 * 10 + 1
                     kclust(1,i2) = kclust(1,i1) / 10 * 10 + 2

                  else

                     kdecay(2) = kdecay(2) + 1

                     kclust(1,i3) = kclust(1,i1) / 10 * 10 + 3
                     kclust(1,i2) = kclust(1,i1) / 10 * 10 + 4

                  end if

               else

                     kdecay(3) = kdecay(3) + 1

                     kclust(1,i3) = kclust(1,i1) / 10 * 10 + 5
                     kclust(1,i2) = kclust(1,i1) / 10 * 10 + 6

               end if


            else

                     kdecay(4) = kdecay(4) + 1

                     kclust(1,i2) = kclust(1,i1) / 10 * 10 + 10
                     kclust(1,i3) = kclust(1,i1) / 10 * 10 + 10

                     if( kclust(1,i2) .ge. 100 ) kclust(1,i2) = 90
                     if( kclust(1,i3) .ge. 100 ) kclust(1,i3) = 90

            end if

*-----------------------------------------------------------------------
*           determine the quantities of two daughters
*-----------------------------------------------------------------------

               lclust(0,i2)  = j2
               lclust(1,i2)  = iz2
               lclust(2,i2)  = in2

               kclust(2,i2)  = i1
               kclust(3,i2)  = kclust(3,i1) + 1

               sclust(1,i2)  = pxr + pcm(1) * trans1
               sclust(2,i2)  = pyr + pcm(2) * trans1
               sclust(3,i2)  = pzr + pcm(3) * trans1
               sclust(4,i2)  = sqrt( em2**2
     &                       + sclust(1,i2)**2
     &                       + sclust(2,i2)**2
     &                       + sclust(3,i2)**2 )
               sclust(5,i2)  = em2
               sclust(6,i2)  = e2
               sclust(7,i2)  = ( sclust(4,i2) - sclust(5,i2) ) * 1000.
               sclust(8,i2)  = wt
               sclust(9,i2)  = 0.0
               sclust(10,i2) = 0.0d0
               sclust(11,i2) = 0.0d0
               sclust(12,i2) = 0.0d0

*-----------------------------------------------------------------------

               lclust(0,i3)  = j3
               lclust(1,i3)  = iz3
               lclust(2,i3)  = in3

               kclust(2,i3)  = i1
               kclust(3,i3)  = kclust(3,i1) + 1

               sclust(1,i3)  = -pxr + pcm(1) * trans2
               sclust(2,i3)  = -pyr + pcm(2) * trans2
               sclust(3,i3)  = -pzr + pcm(3) * trans2
               sclust(4,i3)  = sqrt( em3**2
     &                       + sclust(1,i3)**2
     &                       + sclust(2,i3)**2
     &                       + sclust(3,i3)**2 )
               sclust(5,i3)  = em3
               sclust(6,i3)  = e3
               sclust(7,i3)  = ( sclust(4,i3) - sclust(5,i3) ) * 1000.
               sclust(8,i3)  = wt
               sclust(9,i3)  = 0.0
               sclust(10,i3) = 0.0d0
               sclust(11,i3) = 0.0d0
               sclust(12,i3) = 0.0d0

*-----------------------------------------------------------------------
*     final state
*-----------------------------------------------------------------------

      else if( jdec .eq. 0 ) then


                  kclust(1,i1) = kclust(1,i1) + 100


      end if

*-----------------------------------------------------------------------

 3000    if( num .gt. inum ) goto 2000

*-----------------------------------------------------------------------


      return
      end


************************************************************************
*                                                                      *
      subroutine sdmwid0(iz1,in1,ex1,j1,
     &                   iz2,in2,ex2,j2,
     &                   iz3,in3,ex3,j3,
     &                   totwid,pwid,ifdec,ipdec)
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 27                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to calculate decay width and decay products             *
*              without angular momentum.                               *
*                                                                      *
*              statistical decay of light particles                    *
*              according to the temperature.                           *
*              neutron proton deuteron triton 3he 4he                  *
*              are considered                                          *
*                                                                      *
*                                                                      *
*        Variables:                                                    *
*                                                                      *
*              izi         : charge of fragment                        *
*              ini         : nuetron number of fragment                *
*              exi         : excitation energy                         *
*              ji          : angular momentum                          *
*                                                                      *
*              totwid      : total decay width                         *
*              pwid        : partial decay width                       *
*                                                                      *
*              ifdec       : switch for this subroutine                *
*                      = 0 : calcuate decay width                      *
*                      = 1 : calculate decay product                   *
*                                                                      *
*              ipdec       : decay channel for the decay products      *
*                                                                      *
*                                                                      *
************************************************************************

      use NGSDATAMOD, only : bindeg
      implicit double precision(a-h, o-z)
      include 'param02.inc'

*-----------------------------------------------------------------------

      parameter      ( ddct = - 30.0 )
      parameter      ( mcrit = 1000 )

      parameter      ( denpa = 8.0 )

*-----------------------------------------------------------------------

      common /evapol/ ievmax, izip(0:100), inip(0:100), jip(0:100)
      common /evaval/ sepc(2,0:100)

      dimension   pwid(0:100)
      dimension   qvalp(0:100), eclbp(0:100)
      dimension   smi(0:100), ey1(0:100), s(0:100)

      save        qvalp, eclbp, smi, ey1, s
!$OMP THREADPRIVATE(qvalp, eclbp, smi, ey1, s)
*-----------------------------------------------------------------------

               e1 = ex1
               m1 = iz1 + in1

*-----------------------------------------------------------------------

         if( ifdec .eq. 1 ) goto 3000

*-----------------------------------------------------------------------

                pwid(0) = 0.0

            do ipid = 1, ievmax

               eclbp(ipid) = 0.0
               qvalp(ipid) = 0.0
                pwid(ipid) = 0.0

            end do

*-----------------------------------------------------------------------
*     Calculate total and partial decay width.
*-----------------------------------------------------------------------

               s0 = 2.0 * sqrt( dble(m1) / denpa * e1 )

               totwid = 0.0

*-----------------------------------------------------------------------

      do 100 ipid = 1, ievmax

*-----------------------------------------------------------------------

               s(ipid) = 0.0

*-----------------------------------------------------------------------
*         decay of 1 => 2 + 3
*-----------------------------------------------------------------------

*-----------------------------------------------------------------------
*           Set emission particle
*-----------------------------------------------------------------------

               iz3 = izip(ipid)
               in3 = inip(ipid)
                m3 = iz3 + in3

*-----------------------------------------------------------------------
*           Set residual particle
*-----------------------------------------------------------------------

               iz2 = iz1 - iz3
               in2 = in1 - in3
                m2 = iz2 + in2

*-----------------------------------------------------------------------
*           check
*-----------------------------------------------------------------------

               if( iz3 .gt. iz2 .or.
     &              m3 .gt.  m2 .or.
     &             iz2 .ge.  m2 .or.
     &             iz2 .lt.   0 .or.
     &             in2 .lt.   0       ) goto 100

*-----------------------------------------------------------------------
*           Q-value and Coulomb barrier
*-----------------------------------------------------------------------

               qvalp(ipid) = bindeg(iz1,in1)
     &                     - bindeg(iz2,in2)
     &                     - bindeg(iz3,in3)

               a23 = dble( m2 )**(1./3.)

               eclbp(ipid) =  ccoul * 1000.0
     &                     * dble( iz2 ) * dble( iz3 )
     &                     / ( sepc(1,ipid) + sepc(2,ipid) * a23 )

*-----------------------------------------------------------------------
*           Level density
*-----------------------------------------------------------------------

               ss = e1 - eclbp(ipid) - qvalp(ipid)

                  if( ss .le. 0.0 ) goto 100

               sa = dble(m2) / denpa

               s(ipid)  = 2.0 * sqrt( sa * ss )

            if( s(ipid) - s0 .lt. ddct ) then

                spl = 0.0

            else

                spl = exp( s(ipid) - s0 )

            end if

            if( s(ipid) .gt. -ddct ) then

                smi(ipid) = 0.0

            else

                smi(ipid) = exp( - s(ipid) )

            end if


               ey1(ipid)  = ( ( 2.0 * s(ipid)**2 - 6.0 * s(ipid) + 6.0 )
     &                      + ( s(ipid)**2 - 6.0 ) * smi(ipid) )

*-----------------------------------------------------------------------
*           Partial decay width and total decay width
*-----------------------------------------------------------------------

               pwid(ipid) = dble(m3) * a23**2
     &                    / sa**2 * spl * ey1(ipid)

               if( pwid(ipid) .lt. 0.0 ) pwid(ipid) = 0.0


               totwid = totwid + pwid(ipid)

*-----------------------------------------------------------------------

  100 continue

*-----------------------------------------------------------------------

         return



*-----------------------------------------------------------------------
*     Set decay particles:  ifdec = 1
*-----------------------------------------------------------------------

 3000 continue

*-----------------------------------------------------------------------

               iz3 = izip(ipdec)
               in3 = inip(ipdec)
               ex3 = 0.0

               iz2 = iz1 - iz3
               in2 = in1 - in3

               m2 = iz2 + in2
               m3 = iz3 + in3

               j3 = jip(ipdec)
               j2 = max( 0, j1 - j3 )

*-----------------------------------------------------------------------
*        Determine the energy and momentum.
*-----------------------------------------------------------------------

            if( e1 - qvalp(ipdec) - eclbp(ipdec) .lt. 0.1 .or.
     &          m2 .eq. 1 ) then

               erel = max( 0.0d0, e1 - qvalp(ipdec) )

            else

               sa = dble(m2) / denpa

               t  = ( ( s(ipdec)**3 - 6.0 * s(ipdec)**2
     &                + 15.0 * s(ipdec) -15.0 )
     &              + smi(ipdec) / 8.0
     &              * ( s(ipdec)**4 - 12.0 * s(ipdec)**2
     &                + 15.0 * s(ipdec) ) )
     &              / sa / ey1(ipdec)

               if( t .le. 0.0 ) then

                  erel = max( 0.0d0, e1 - qvalp(ipdec) )
                  goto 330

               end if

               ermax = e1 - qvalp(ipdec)
               smax  = t / 2.718282


               itry = 0
   20          itry = itry + 1

               if( itry .ge. mcrit ) goto 330

               erel = eclbp(ipdec) + ( ermax - eclbp(ipdec) ) * rn(0)
               erc  = erel - eclbp(ipdec)

               if( -erc / t .lt. ddct ) goto 20

               if( erc * exp( -erc / t ) .lt. smax * rn(0) ) goto 20

            end if

  330          continue

               ex2  = max( 0.0d0, e1 - erel - qvalp(ipdec) )


*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine sdmwid1(iz1,in1,e14,j1,
     &                   iz2,in2,e24,j2,
     &                   iz3,in3,e34,j3,
     &                   tcs04,tcip04,ifdec,ipdec)
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 27                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to calculate decay width and decay products             *
*              with angular momentum.                                  *
*                                                                      *
*              statistical decay of light particles                    *
*              according to the temperature.                           *
*              neutron proton deuteron triton 3he 4he                  *
*              are considered                                          *
*                                                                      *
*                                                                      *
*        Variables:                                                    *
*                                                                      *
*              izi         : charge of fragment                        *
*              ini         : nuetron number of fragment                *
*              exi         : excitation energy                         *
*              ji          : angular momentum                          *
*                                                                      *
*              tcs04       : total decay width                         *
*              tcip04      : partial decay width                       *
*                                                                      *
*              ifdec       : switch for this subroutine                *
*                      = 0 : calcuate decay width                      *
*                      = 1 : calculate decay product                   *
*                                                                      *
*              ipdec       : decay channel for the decay products      *
*                                                                      *
*                                                                      *
************************************************************************
      use NGSDATAMOD, only : bindeg

      implicit double precision(a-h, o-z)

      include 'err.inc'

      dimension  tcip04(0:100)

*-----------------------------------------------------------------------

      common /evapol/ ievmax, izip(0:100), inip(0:100), jip(0:100)
      common /evaval/ sepc(2,0:100)

      common /sdmcut/ sdmemin
!$OMP THREADPRIVATE(/sdmcut/)

      common /sdmsw2/ imengb, imbarr, imangm, imlevd, imgamm
!$OMP THREADPRIVATE(/sdmsw2/)

*-----------------------------------------------------------------------
*     Some constants
*-----------------------------------------------------------------------

      parameter ( rmass = 938.3d0, hc  = 197.3d0 )
      parameter ( ccoul = 1.439767d0 )
      parameter ( pi = 3.141592653589793d0 )
      parameter ( xi1 = 0.5d-8, xi2 = 0.3d-9 )

*-----------------------------------------------------------------------

      parameter( maxtc = 5200 )

      dimension tcsave(maxtc), tcip(0:100), ntcoffset(0:100)
      dimension j2idx(maxtc),ieidx(maxtc),lmaxsv(maxtc)

      dimension de(1000), ee(1000), we(1000)
      dimension itmp1(1000),itmp2(1000),ltmp(1000),iesv(1000)

*-----------------------------------------------------------------------

      save     tcsave, ntcoffset
      save     j2idx,ieidx,lmaxsv,ng
!$OMP THREADPRIVATE(tcsave, ntcoffset,j2idx,ieidx,lmaxsv,ng)
*-----------------------------------------------------------------------
*           iemax : number of energy bin
*           ifj   : switch of angular momentum
*           emin  : minimum energy of decay particles
*-----------------------------------------------------------------------

            iemax = imengb
            ifj   = imangm
            emin  = sdmemin

*-----------------------------------------------------------------------

            e1    = e14
            j1    = j1 * ifj

*-----------------------------------------------------------------------

            if( ifdec .eq. 1 ) then

                  tcs    = tcip(ipdec)
                  ifsave = 0
                  ifcalc = 0

            else

                  ifsave = 1
                  ifcalc = 1

               do ip = 0, ievmax

                  tcip(ip) = 0

               end do

            end if

*-----------------------------------------------------------------------

               ntc    = 0
               random = rn(0)

*-----------------------------------------------------------------------

               m1     = iz1 + in1
               r1     = 1.2 * dble( m1 )**(1./3.)
               de1    = emin
               dl1    = sdmlev(iz1,in1,e1,de1,j1)
               denomi = dl1 * 2 * pi * de1

*-----------------------------------------------------------------------
*        dl1 = 0.0 => parent is in the wrong state. but continue.
*-----------------------------------------------------------------------

            if( denomi .lt. 0.0001 ) then

               goto 888

            end if

*-----------------------------------------------------------------------
*        gamma deacy
*-----------------------------------------------------------------------

               ip    = 0
               tcip0 = 0

               ntcoffset(ip) = ntc

         if( imgamm .eq. 1 .and.
     &     ( ifdec .eq. 0 .or. ip .eq. ipdec ) ) then

               iz2 = iz1
               in2 = in1
                m2 = m1

               iz3 = 0
               in3 = 0
                m3 = 0
                e3 = 0
                j3 = 0

*-----------------------------------------------------------------------
*           make energy bin
*-----------------------------------------------------------------------

                  emax = e1
                  pow  = 1.0 / ( 0.4 * sqrt( emax * m1 / 8.0 ) )

               do ie2 = 1, iemax

                  yy = ( ie2 - 0.5 ) / ( iemax - 1 )
                  if( ie2 .eq. iemax ) yy = 1.0

                  we(ie2) = emax * yy**pow

               end do

                  ee(1) = we(1) / 2.0
                  de(1) = we(1)

               do ie2 = 2, iemax

                  ee(ie2) = ( we(ie2) + we(ie2-1) ) / 2.0
                  de(ie2) =   we(ie2) - we(ie2-1)

               end do

*-----------------------------------------------------------------------
*           gamma width
*-----------------------------------------------------------------------

                  itmp  = ifj * ( j1 + 2 ) - ifj * max( 0, j1 - 2 ) + 1
                  max10 = iemax * itmp

*-----------------------------------------------------------------------
            if( ifdec .eq. 0 ) then
*-----------------------------------------------------------------------

               do 10 i = 1, max10

                  ie2 = ( i - 1 ) / itmp + 1
                  j2  = ifj * max( 0, j1 - 2 ) + mod( i - 1 , itmp)

                  ieidx(i) = ie2
                  j2idx(i) = j2

                  de2  = de(ie2)
                  e2   = ee(ie2)
                  erel = e1 - e2

                  dl = sdmlev(iz1,in1,e2,de2,j2)

                  if( ( abs( j2 - j1 ) .le. 2 ) .and.
     &                ( j2 + j1 + mod(m1,2) .ge. 2 ) ) then

                        tcoeff = xi2 * erel**5 * dl * de2 / denomi

                     if( ( abs( j2 - j1 ) .le. 1 ) .and.
     &                   ( j2 + j1 + mod(m1,2) .ge. 1 ) ) then

                        tcoeff = xi2 * erel**5 * dl * de2 / denomi
     &                         + xi1 * erel**3 * dl * de2 / denomi

                     end if

                  else

                        tcoeff = 0.0

                  end if

                        tcsave(i) = tcoeff
                        tcip0     = tcip0 + tcoeff

   10          continue

                        ntc = max10

*-----------------------------------------------------------------------
            else
*-----------------------------------------------------------------------

              do 101 i = 1, max10

                    ie2   = ieidx(i)
                    j2    = j2idx(i)
                    de2   = de(ie2)
                    e2    = ee(ie2)
                    erel  = e1 - e2
                    tcip0 = tcip0 + tcsave(i)

                    if( tcip0 .ge. tcs * random ) goto 999

  101         continue

            end if

*-----------------------------------------------------------------------
         end if
*-----------------------------------------------------------------------

                     tcip(ip) = tcip0


*-----------------------------------------------------------------------
*     particle decay
*-----------------------------------------------------------------------

      do 30 ip = 1, ievmax

               if( ifdec .eq. 0 ) then

                  ntcoffset(ip) = ntc

               else if( ifdec .eq. 1 ) then

                  ntc = ntcoffset(ip)

               end if

*-----------------------------------------------------------------------

         if( ifdec .eq. 0 .or. ip .eq. ipdec ) then

*-----------------------------------------------------------------------

                  iz3 = izip(ip)
                  in3 = inip(ip)
                   j3 = jip(ip) * ifj

                   m3 = iz3 + in3
                   e3 = 0

                  iz2 = iz1 - iz3
                  in2 = in1 - in3
                   m2 = iz2 + in2

                  tcip0 = 0

*-----------------------------------------------------------------------

            if(  m2 .gt. 0 .and.
     &          iz2 .ge. 0 .and.
     &          in2 .ge. 0 ) then

                  r2 = 1.2 * dble( m2 )**(1./3.)
                  r3 = 1.2 * dble( m3 )**(1./3.)

                  qval12 =     ( bindeg(iz1,in1)
     &                         - bindeg(iz2,in2)
     &                         - bindeg(iz3,in3) )

               if( imbarr .eq. 1 ) then

                  a23 = dble( m2 )**(1./3.)
                  b23 = iz2 * iz3 * ccoul
     &                / ( sepc(1,ip) + sepc(2,ip) * a23 )

               else

                  b23 = iz2 * iz3 * ccoul / ( r2 + r3 )

               end if

                  e2max = e1 - qval12 - b23

*-----------------------------------------------------------------------

            if( e2max .gt. 0 ) then

*-----------------------------------------------------------------------
*              make energy bin
*-----------------------------------------------------------------------

                  emax = e2max
                  pow  = 1.0 / ( 0.4 * sqrt( emax * m1 / 8.0 ) )

               do ie2 = 1, iemax

                  yy = ( ie2 - 0.5 ) / ( iemax - 1 )
                  if( ie2 .eq. iemax ) yy = 1.0

                  we(ie2) = emax * yy**pow

               end do

                  ee(1) = we(1) / 2.0
                  de(1) = we(1)

               do ie2 = 2, iemax

                  ee(ie2) = ( we(ie2) + we(ie2-1) ) / 2.0
                  de(ie2) =   we(ie2) - we(ie2-1)

               end do


*-----------------------------------------------------------------------
            if( ifdec .eq. 0 ) then
*-----------------------------------------------------------------------

                           ntcold = ntc
                           ntctmp = 0

                  do 40 ie2 = 1, iemax

                           de2  = de(ie2)
                           e2   = ee(ie2)
                           erel = e1 - e2 - qval12

                           angmx2 = 2.0 * rmass * m3 * m2 / m1
     &                            * ( r2 + r3 )**2
     &                            * ( erel - b23 )

                        if( angmx2 .ge. 0 ) then

                           iesv(ie2)  = 1
                           ltmp(ie2)  = int( sqrt(angmx2) / hc )
                           itmp1(ie2) = ifj
     &                                * max( 0, j1 - ltmp(ie2) - 1 )
                           itmp2(ie2) = ifj * ( J1 + ltmp(ie2) + 1 )

                           ntctmp = ntctmp + itmp2(ie2) - itmp1(ie2) + 1

                        else

                           iesv(ie2) = 0

                        end if

  40              continue

                     if( ntctmp + ntc .gt. maxtc ) then

                           ng = ip


                           ErrCha = ''
                           ErrID = 'L:1539/R:sdmwid1/F:sdml.f' !W00_007_001
                           call ErrWrite(ErrID,ErrCha)

                           write(*,*)
     &                     ' **** Warning at sdmwid1, over maxtc'
                           write(*,*)
     &                     ' ==================================='
                           write(*,*) ' over maxtc : ntc = ', ntc

                     else

                           ng = -1

                     end if


               if( ng .lt. 0 ) then

                  do 41 ie2 = 1, iemax

                     if( iesv(ie2) .eq. 1 ) then

                        do 50 j2 = itmp1(ie2), itmp2(ie2)

                           ntc = ntc + 1

                           j2idx(ntc)  = j2
                           ieidx(ntc)  = ie2
                           lmaxsv(ntc) = ltmp(ie2)

  50                    continue

                     end if

  41              continue

                  do 444 i = ntcold + 1, ntc

                           ie2  = ieidx(i)
                           j2   = j2idx(i)
                           de2  = de(ie2)
                           e2   = ee(ie2)
                           lmax = lmaxsv(i)

                        if( ifj .eq. 1 ) then

                           call sdjsum( 2 * J1 + mod(m1,2),
     &                                  2 * J2 + mod(m2,2),
     &                                  2 * J3 + mod(m3,2),
     &                                  lmax, w )

                        else

                           w = 1

                        end if

                           dl = sdmlev(iz2,in2,e2,de2,j2)

                           tcoeff    = w * dl * de2 / denomi
                           tcip0     = tcip0 + tcoeff
                           tcsave(i) = tcoeff

 444              continue


               else


                  do 411 ie2 = 1, iemax

                     if( iesv(ie2) .eq. 1 ) then

                        do 501 j2 = itmp1(ie2), itmp2(ie2)

                           ntc  = ntc + 1
                           lmax = ltmp(ie2)
                           de2  = de(ie2)
                           e2   = ee(ie2)

                           if( ifj .eq. 1 ) then

                              call sdjsum( 2 * j1 + mod(m1,2),
     &                                     2 * J2 + mod(m2,2),
     &                                     2 * J3 + mod(m3,2),
     &                                     lmax, w )

                           else

                              W = 1

                           end if

                           dl = sdmlev(iz2,in2,e2,de2,j2)

                           tcoeff = w * dl * de2 / denomi
                           tcip0  = tcip0 + tcoeff

 501                    continue

                     end if

 411              continue

               end if

                           ntcoffset(ip+1) = ntc

*-----------------------------------------------------------------------
            else
*-----------------------------------------------------------------------

               if( ng .lt. ip ) then

*vocl loop,scalar

                  do 441 i = ntcoffset(ip) + 1, ntcoffset(ip+1)

                           tcip0 = tcip0 + tcsave(i)

                     if( tcip0 .ge. tcs * random ) then

                           ie2  = ieidx(i)
                           de2  = de(ie2)
                           e2   = ee(ie2)
                           erel = e1 - e2 - qval12
                           lmax = lmaxsv(i)
                           j2   = j2idx(i)

                           goto 999

                     end if

  441             continue

               else

*vocl loop,scalar

                  do 44 ie2 = 1, iemax

                           de2  = de(ie2)
                           e2   = ee(ie2)
                           erel = e1 - e2 - qval12

                           angmx2 = 2.0 * rmass * m3 * m2 / m1
     &                            * ( r2 + r3 )**2
     &                            * ( erel - b23 )

                     if( angmx2 .ge. 0 ) then

                           lmax = int( sqrt(angmx2) / hc )
*vocl loop,scalar

                        do 55 J2 = ifj * max( 0, j1 - lmax - 1 ),
     &                             ifj * ( j1 + lmax + 1 )

                           ntc = ntc + 1

                           if( ifj .eq. 1 ) then

                              call sdjsum( 2 * j1 + mod(m1,2),
     &                                     2 * j2 + mod(m2,2),
     &                                     2 * j3 + mod(m3,2),
     &                                     lmax, w )

                           else

                              w = 1

                           end if

                           dl = sdmlev(iz2,in2,e2,de2,j2)

                           tcoeff = w * dl * de2 / denomi
                           tcip0  = tcip0 + tcoeff

                           if( tcip0 .ge. tcs * random ) goto 999

  55                    continue

                     end if

  44              continue

               end if

            end if


*-----------------------------------------------------------------------

            end if

            end if

                  tcip(ip) = tcip0

*-----------------------------------------------------------------------

         end if

*-----------------------------------------------------------------------
   30 continue
*-----------------------------------------------------------------------


            if( ifdec .eq. 1 ) then

             ErrCha = ''
             ErrID = 'L:1749/R:sdmwid1/F:sdml.f' !W00_008_001
             call ErrWrite(ErrID,ErrCha)

               write(*,*) ' **** Warning at sdmwid1, no decay channel'
               write(*,*) ' ========================================='
               write(*,*) '  iz1,in1,e1,j1,tcs = ', iz1,in1,e1,j1,tcs
               write(*,*) '  ipdec,tcs,tcip(ipdec) = ',
     &                       ipdec,tcs,tcip(ipdec)

            end if

*-----------------------------------------------------------------------

  888    continue

                  tcs   = 0.0
                  tcs04 = 0.0

               do 100 ip = 0, ievmax

                  tcs = tcs + tcip(ip)
                  tcip04(ip) = (tcip(ip))
                  tcs04 = tcs04 + tcip04(ip)

  100          continue

               return

*-----------------------------------------------------------------------

  999 continue

               e24   = (e2)
               e34   = (e3)


*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      function sdmlev(iz,in,e,de,j)
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 27                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to calculate level density according to                 *
*              Ref. F.Puhlhofer, Nucl. Phys. A280(1977)267.            *
*                                                                      *
*                                                                      *
*        Variables:                                                    *
*                                                                      *
*              iz          : proton number                             *
*              in          : neutron number                            *
*              e           : excitation energy (MeV)                   *
*              de          : energy bin (MeV)                          *
*              j           : angluar momentum                          *
*                                                                      *
*                                                                      *
************************************************************************

      implicit double precision(a-h, o-z)

*-----------------------------------------------------------------------

      common /sdmsw2/ imengb, imbarr, imangm, imlevd, imgamm
!$OMP THREADPRIVATE(/sdmsw2/)

      common /sdmcut/ sdmemin
!$OMP THREADPRIVATE(/sdmcut/)

*-----------------------------------------------------------------------

      parameter( rmass = 938.3d0, hc  = 197.3d0 )
      parameter( denpa = 8.0d0 )

*-----------------------------------------------------------------------
*              initial value
*-----------------------------------------------------------------------

               denl = 0

               if( e .lt. sdmemin .and.
     &             j .eq. 0 ) denl = 1.d0 / de

               m    = iz + in
               amas = m * 1.0d0
               a    = amas / denpa

*-----------------------------------------------------------------------
*     a =< 4
*-----------------------------------------------------------------------

            if( m .le. 4 ) then

                  sdmlev = denl

                  return

            end if

*-----------------------------------------------------------------------
*     a > 4
*-----------------------------------------------------------------------
*        simple level density
*-----------------------------------------------------------------------

         if( imlevd .eq. 0 ) then

               denl = exp( 2 * sqrt( a * e ) ) / 1000

*-----------------------------------------------------------------------
*        Level density from F. Puhlhofer.
*-----------------------------------------------------------------------

         else

            if( ( m / 2 ) * 2 .ne. m ) then

               delta = -0.7d0

            else if( ( iz / 2 ) * 2. ne. iz ) then

               delta = -2.0d0

            else

               delta = 0.7d0

            end if

               r = 4 * ( amas**(5.d0/3.d0) )
     &           * rmass * ( 1.2d0**2 ) / 5.0d0 / hc / hc / a

               u1 = e - ( j + mod(m,2) / 2.0d0 )**2 / a / r - delta
               u2 = e - ( j + mod(m,2) / 2.0d0 + 1 )**2 / a / r - delta

            if( u2 .gt. 0 ) then

               t1 = ( 1.5d0 + sqrt( 1.5d0 * 1.5d0 + 4 * a * u1 ) )
     &            / 2.0d0 / a

               t2 = ( 1.5d0 + sqrt( 1.5d0 * 1.5d0 + 4 * a * u2 ) )
     &            / 2.0d0 / a

               denl = 1.0d0 / 12.0d0 / sqrt(r) / a / a
     &              * ( exp( 2.0d0 * sqrt( a * u1 ) )
     &              / t1 / t1 / t1
     &              - exp( 2.0d0 * sqrt( a * u2 ) )
     &              / t2 / t2 / t2 )

            end if

         end if


               sdmlev = denl

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine sdjsum(j1,j2,j3,lmax,w)
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 27                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to sum up angular momentum.                             *
*                                                                      *
*        Variables:                                                    *
*                                                                      *
*              j1, j2, j3  : angular momentum                          *
*              lmax        : maximum relative angular momentum         *
*              w           : number of states                          *
*                                                                      *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------

      implicit real*8(a-h,o-z)

*-----------------------------------------------------------------------

            w = 0.0

            if( mod( j1 + j2 + j3, 2 ) .eq. 1 ) return

            jj1 = abs( j2 - j3 )
            jj2 = j2 + j3

*-----------------------------------------------------------------------

            if( jj2 - jj1 .eq. 0 ) then

               l1 = abs( j1 - jj1 ) / 2
               l2 = min( ( j1 + jj1 ) / 2, lmax )

               w  = dim( l2 + 1 - l1, 0 )

            else if( jj2 - jj1 .eq. 2 ) then

               l1 = abs( j1 - jj1 ) / 2
               l2 = min( ( j1 + jj1 ) / 2, lmax )
               l3 = abs( j1 - ( jj1 + 2 ) ) / 2
               l4 = min( ( j1 + jj1 + 2 ) / 2, lmax )

               w  = dim( l2 + 1 - l1, 0 ) + dim( l4 + 1 - l3, 0 )

            else

               l1 = abs( j1 - jj1 ) / 2
               l2 = min( ( j1 + jj1 ) / 2, lmax )
               l3 = abs( j1 - ( jj1 + 2 ) ) / 2
               l4 = min( ( j1 + jj1 + 2 ) / 2, lmax )
               l5 = abs( j1 - ( jj1 + 4 ) ) / 2
               l6 = min(( j1 + jj1 + 4 ) / 2, lmax )

               w  = dim( l2 + 1 - l1, 0 )
     &            + dim( l4 + 1 - l3, 0 )
     &            + dim( l6 + 1 - l5, 0 )

            end if

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      subroutine sdmfisw(iz1,in1,e1,j1,tcn,tcsf,barf,saf)
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 27                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to calculate fission width                              *
*                                                                      *
*                                                                      *
************************************************************************
      use NGSDATAMOD, only : bindeg

      implicit double precision(a-h, o-z)

*-----------------------------------------------------------------------

      parameter      ( ddct = - 30.0 )

      parameter      ( sanf  =  10.0,
     &                 saa   =  7.9169e-08,
     &                 sbb   = -1.6143e-04,
     &                 scc   =  1.0922,
     &                 rk0   =  13.0 )

*-----------------------------------------------------------------------

               tcsf  = 0.0

               m1 = iz1 + in1

               if( m1 .lt. 50 ) return

*-----------------------------------------------------------------------

                     san   =  dble(m1) / sanf
                     saf   = ( saa * e1**2 + sbb * e1 + scc ) * san

                     qvaln = bindeg(iz1,in1) - bindeg(iz1,in1-1)
                     ss    = e1 - qvaln

                     barf  = sdmfisb(m1,iz1)

               if( e1 - qvaln .gt. 0.0 .and.
     &             e1 - barf  .gt. 0.0 ) then

                        ssn  = 2.0 * sqrt( san * ( e1 - qvaln ) )
                        ssf  = 2.0 * sqrt( saf * ( e1 - barf  ) )

                        ssfn = ssf - ssn

                     if( ssfn .gt. -ddct ) ssfn = -ddct

                        ratfn = 0.0

                     if( ssfn .ge. ddct ) then

                        ratfn = rk0 * san * ( ssf - 1.0 )
     &                        / ( 4.0 * dble(m1)**(2./3.)
     &                          * saf * ( e1 - qvaln )  )
     &                        * exp( ssfn )
                     end if

                        if( ratfn .gt. 0.0 ) tcsf = tcn * ratfn

               end if

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine sdmfiss(iz1,in1,e1,j1,
     &                   iz2,in2,e2,j2,
     &                   iz3,in3,e3,j3,
     &                   barf,saf,ierr)
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 27                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to determine the mass and charge of fission fragment    *
*                                                                      *
*        Variables:                                                    *
*                                                                      *
*              iz1         : charge of fragment                        *
*              in1         : nuetron number of fragment                *
*              e1          : excitation energy                         *
*              j1          : angular momentum                          *
*                                                                      *
*              (iz2,in2,e2,i2)                                         *
*              (iz3,in3,e3,j3)  : fragments' quantities                *
*                                                                      *
*              barf        : fission barrier                           *
*              saf         : level density parameter of fission        *
*              ierr        : error flag                                *
*                                                                      *
*                                                                      *
************************************************************************
      use NGSDATAMOD, only : bindeg

      implicit double precision(a-h, o-z)

      include 'param02.inc'

      include 'err.inc'

*-----------------------------------------------------------------------

      parameter      ( ddct = - 30.0 )
      parameter      ( bnorm = 1.2011 )
      parameter      ( mcrit = 1000 )

*-----------------------------------------------------------------------
*           Error Frag
*-----------------------------------------------------------------------

               ierr = 0

*-----------------------------------------------------------------------
*           temperature and width of the mass distribution
*-----------------------------------------------------------------------

            m1   = iz1 + in1
            rm1  = dble(m1)
            rz1  = dble(iz1)

            temp = dsqrt( dabs( ( e1 - barf ) / saf ) )

            width = ( e1 - barf + 7.0 ) * bnorm

*     -------------------------------------------------------
*     this width is questionable for high excitation energy !
*     -------------------------------------------------------

*-----------------------------------------------------------------------
*     mass of fission fragment for subactinides
*-----------------------------------------------------------------------

      if( iz1 .lt. 89 ) then

            imsa = 0
   50       imsa = imsa + 1

            if( imsa .lt. mcrit ) then

               rm2  = rn(0) * ( rm1 - 1.0 )
               rarg = - ( rm2 - rm1/2.0 )**2 / width**2
               m2   = nint( rm2 )

               if( rarg .lt. ddct )            goto 50
               if( rn(0) .gt. exp( rarg ) )    goto 50
               if( m2 .eq. 0 .or. m2 .eq. m1 ) goto 50

            else

               m2 = m1 / 2

            end if


*-----------------------------------------------------------------------
*     mass of fission fragment for actinides
*-----------------------------------------------------------------------

      else

            bm1 = 0.4 * rm1
            bm2 = 0.5 * rm1
            bm3 = 0.6 * rm1

            exm = e1 + 6.0

            alpa = 19.98160
            beta = 78.61184

         if( exm .le. 25.0 ) then

            alpa = exp( 0.5991 * exm - 13.1869 )
            beta = exp( 0.7013 * exm - 17.5325 )

         else if( exm .le. 40.0 ) then

            alpa = exp( 0.2008 * exm**0.8 - 0.8451 )
            beta = exp( 2.2672 * sqrt(exm)- 11.3431 )

         else if( exm .le. 48.0 ) then

            beta = exp( 2.2672 * sqrt(exm)- 11.3431 )

         end if

            rdm1 = alpa + beta * exp( - ( bm1 - bm2 )**2 / width**2 )
     &           + alpa * exp( - ( bm1 - bm3 )**2 / width**2 )

            rdm2 = beta
     &           + 2.0 * alpa * exp( - ( bm1 - bm2 )**2 / width**2 )

            rdmax = max( rdm1, rdm2 ) * 1.1


            imsa = 0
   60       imsa = imsa + 1

            if( imsa .lt. mcrit ) then

               rm2  = rn(0) * ( rm1 - 1.0 )

                  ear1 = 0.0
                  ear2 = 0.0
                  ear3 = 0.0

                  rar1 = - ( rm2 - bm1 )**2 / width**2
                  if( rar1 .gt. ddct ) ear1 = exp( rar1 )

                  rar2 = - ( rm2 - bm2 )**2 / width**2
                  if( rar2 .gt. ddct ) ear2 = exp( rar2 )

                  rar3 = - ( rm2 - bm3 )**2 / width**2
                  if( rar3 .gt. ddct ) ear3 = exp( rar3 )

               m2 = nint( rm2 )

               if( rn(0) * rdmax .gt. ear1 + ear2 + ear3 ) goto 60
               if( m2 .eq. 0 .or. m2 .eq. m1 )             goto 60

          else

               m2 = m1 / 2

          end if

      end if

*-----------------------------------------------------------------------
*     deetermine the charge of fission fragment
*-----------------------------------------------------------------------

         rho  = 1.1
         r0   = 1.2
         bets = -31.4506
         phi  = 44.2355

         s    = 0.1 * ccoul * 1000.0 / r0 / bets
     &          * ( rm1 / 2. )**(2./3.) * ( 1.0 - 5. / 8. / rho )

         zbar = - s * rz1 / 2. + ( 1. + s ) * dble( m2 ) * rz1 / rm1


         widi = - 16.0 * bets / rm1 / temp
     &        * ( 1. + phi / bets * ( 2. / rm1 )**(1./3.)
     &            - 0.055 * ccoul * 1000.0 / r0 / bets * rm1**(2./3.) )

            imsa = 0
   70       imsa = imsa + 1

            if( imsa .lt. mcrit ) then

               rz2  = rn(0) * ( rz1 - 1.0 )
               rarg = - ( rz2 - zbar )**2 * widi

               iz2  = nint( rz2 )

               if( rarg .lt. ddct )         goto 70
               if( rn(0) .gt. exp( rarg ) ) goto 70
               if( iz2 .eq.  0 .or.
     &             iz2 .ge. m2 .or.
     &             iz2 .eq. iz1 )           goto 70

          else

               iz2 = min( m2 / 2 , iz1 / 2 )

          end if


*-----------------------------------------------------------------------
*     now fission fragments are determined
*-----------------------------------------------------------------------

            m3  = m1 - m2

            iz2 = iz2 - 1
   80       iz2 = iz2 + 1

            if( iz2 .lt. m2 .and. iz2 .le. iz1 ) then

               iz3 = iz1 - iz2
               if( iz3 .ge. m3 ) goto 80

            else

               ErrCha = ''
               ErrID = 'L:2288/R:sdmfiss/F:sdml.f' !E81_002_001
               call ErrWrite(ErrID,ErrCha)

               write(*,*) ' **** Error at sdmfiss,',
     &                    ' iz3 can not be determined'
               write(*,*) ' ======================',
     &                    '=========================='
               write(*,*) ' iz1 = ', iz1, ' in1 = ', in1

               ierr = 1
               return

            end if

            in2 = m2 - iz2
            in3 = m3 - iz3

            rm = dble( m2 * m3 ) / dble( m1 )


*-----------------------------------------------------------------------
*     determine the kinetic energy and excitation energy
*-----------------------------------------------------------------------

            rad2 = 1.2 * dble(m2)**(1./3.)
            rad3 = 1.2 * dble(m3)**(1./3.)

            rad  = rad2 + rad3

            amoi2  = 2./5. * m2 * rad2**2
            amoi3  = 2./5. * m3 * rad3**2
            amoi23 = rm * rad**2 + amoi2 + amoi3

            xj2 = j1 * amoi2 / amoi23
            xj3 = j1 * amoi3 / amoi23

            j2  = xj2
            j3  = xj3

*-----------------------------------------------------------------------
*           rotation energy and relative motion
*-----------------------------------------------------------------------

            er = ( j1 - xj2 - xj3)**2 / 2.0
     &         / ( amoi23 - amoi2 - amoi3 )

            er = er + 0.1071 * rz1**2 / rm1**(1./3.) + 22.2

            qvalf = bindeg(iz1,in1)
     &            - bindeg(iz2,in2)
     &            - bindeg(iz3,in3)

            etf =  max( 0.0d0, e1 - qvalf )
            efx = etf - er

            if( efx .le. 0.0 ) then

                er  = etf
                efx = 0.0

            end if

            e2 = efx / rm1 * m2
            e3 = efx / rm1 * m3

*-----------------------------------------------------------------------

       return
       end


************************************************************************
*                                                                      *
      function sdmfisb(m1,iz1)
*                                                                      *
*                                                                      *
*        Last Revised:     1998 11 27                                  *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              to determine the fission barrier.                       *
*                                                                      *
*        Variables:                                                    *
*                                                                      *
*              m1          : mass of parent nucleus                    *
*              iz1         : charge of parent nucleus                  *
*              sdmfisb     : fission barrier                           *
*                                                                      *
*                                                                      *
*        Comments:                                                     *
*                                                                      *
*              These double-humped fission barrier data are            *
*              taken from  Kupriyanof et al.                           *
*              Sov. J. Nucl. Phys. 32 (1980) 184                       *
*                                                                      *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------
C...Double precision.
      implicit double precision(a-h, o-z)

      dimension  fba(64), fbb(64)

*-----------------------------------------------------------------------

      data fba / 5.69, 5.50, 5.08, 5.07, 5.69, 5.42, 5.26, 5.45, 5.59,
     &           5.42, 5.48, 5.45, 5.41, 5.29, 5.24, 5.57, 5.45, 5.65,
     &           5.48, 5.67, 5.79, 5.97, 5.92, 5.91, 5.83, 5.89, 5.65,
     &           5.63, 5.42, 5.69, 6.10, 6.02, 6.17, 5.96, 6.17, 5.82,
     &           6.25, 6.22, 6.40, 6.16, 6.17, 5.94, 5.92, 5.71, 5.67,
     &           6.40, 6.59, 6.34, 6.44, 6.09, 6.26, 5.82, 5.92, 5.37,
     &           6.56, 6.45, 6.53, 6.41, 6.54, 6.32, 6.32, 6.10, 5.89,
     &           5.48 /

      data fbb / 7.89, 7.68, 7.30, 7.34, 7.35, 7.14, 7.04, 6.58, 6.79,
     &           6.68, 6.80, 6.84, 6.86, 6.79, 6.80, 6.26, 6.21, 6.48,
     &           6.38, 5.75, 5.95, 6.20, 6.23, 6.29, 6.28, 6.40, 6.23,
     &           6.26, 6.12, 5.21, 5.69, 5.68, 5.93, 5.79, 6.08, 5.79,
     &           5.34, 5.39, 5.65, 5.48, 5.46, 5.41, 5.52, 5.32, 5.34,
     &           4.87, 5.15, 4.98, 5.16, 4.89, 5.13, 4.77, 4.94, 4.45,
     &           4.50, 4.38, 4.54, 4.50, 4.72, 4.57, 4.65, 4.50, 4.36,
     &           4.02 /

*-----------------------------------------------------------------------

      if( m1 .le. 90 ) then

            sdmfisb = 52.0 * exp( - ( ( m1 - 90.0 ) / 84.7  )**2 )

      else if( m1 .le. 200 ) then

            sdmfisb = 52.0 * exp( - ( ( m1 - 90.0 ) / 110.0 )**2 )

      else if( m1 .le. 224 ) then

            sdmfisb = 23.0 * exp( - ( ( m1 - 210.0 ) / 13.2 )**2 )

      else

            sdmfisb = 6.0

            ifis = 0

         if( iz1 .eq. 88 )      then

            if( m1 .le. 228 )                   ifis = 1  + m1 - 225

         else if( iz1 .eq. 89 ) then

            if( m1 .ge. 226 .and. m1 .le. 228 ) ifis = 5  + m1 - 226

         else if( iz1 .eq. 90 ) then

            if( m1 .ge. 227 .and. m1 .le. 234 ) ifis = 8  + m1 - 227

         else if( iz1 .eq. 91 ) then

            if( m1 .ge. 230 .and. m1 .le. 233 ) ifis = 16 + m1 - 230

         else if( iz1 .eq. 92 ) then

            if( m1 .ge. 231 .and. m1 .le. 240 ) ifis = 20 + m1 - 231

         else if( iz1 .eq. 93 ) then

            if( m1 .ge. 233 .and. m1 .le. 239 ) ifis = 30 + m1 - 233

         else if( iz1 .eq. 94 ) then

            if( m1 .ge. 237 .and. m1 .le. 245 ) ifis = 37 + m1 - 237

         else if( iz1 .eq. 95 ) then

            if( m1 .ge. 239 .and. m1 .le. 247 ) ifis = 46 + m1 - 239

         else if( iz1 .eq. 96 ) then

            if( m1 .ge. 241 .and. m1 .le. 250 ) ifis = 55 + m1 - 241

         end if


         if( ifis .ne. 0 ) then

            sdmfisb = max( fba(ifis), fbb(ifis) )

         end if


      end if

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      block data evtable
*                                                                      *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              data table of evapolation particle                      *
*                                                                      *
*                                                                      *
************************************************************************

      implicit double precision(a-h, o-z)

*-----------------------------------------------------------------------

      common /evapol/ ievmax, izip(0:100), inip(0:100), jip(0:100)
      common /evaval/ sepc(2,0:100)

*-----------------------------------------------------------------------
*     data for evapolation particle
*-----------------------------------------------------------------------

      data ievmax / 6 /

      data ( izip(i), i = 1, 6 ) / 0,1,1,1,2,2 /
      data ( inip(i), i = 1, 6 ) / 1,0,1,2,1,2 /
      data (  jip(i), i = 1, 6 ) / 0,0,1,0,0,0 /

      data ( sepc(1,i), sepc(2,i), i = 1, 6 ) /
     &            1.0000, 1.00000,  !neutron DUMMY
     &            8.0606, 0.50836,  !proton
     &            7.9869, 0.62142,  !deuteron
     &            7.9869, 0.62142,  !trition
     &            6.4355, 0.78601,  !3-helium
     &            6.1333, 0.88761/  !alpha particle

c  original

c  modified

*-----------------------------------------------------------------------

      end


