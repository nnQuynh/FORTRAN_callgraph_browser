************************************************************************
*                                                                      *
      subroutine erupin(iz,in,ex,px,py,pz,pt,et,rm,wt)
*                                                                      *
*                                                                      *
*       main control routine of Evaporation Code ERUP and Fission      *
*       and booking produced particles and nucleus                     *
*       last modified by K.Niita on 01/06/2001                         *
*                                                                      *
*     input:                                                           *
*                                                                      *
*        iz, in     : proton and neutron number of mather              *
*        ex         : excitation energy of mather (MeV)                *
*        px,py,pz   : momentum vector of mather (GeV)                  *
*        pt         : absolute value of momentum of mather (GeV)       *
*        et         : energy of mather (sqrt(p**2+m**2) GeV)           *
*        rm         : rest mass of mather (GeV)                        *
*        wt         : weight change                                    *
*                                                                      *
*     output:                                                          *
*                                                                      *
*---- in common -------------------------------------------------------*
*                                                                      *
*        nclust   : total number of out going particles and nuclei     *
*                                                                      *
*        kclust(3,nclust)                                              *
*                                                                      *
*                   kclust(1,i) = 101 : final output                   *
*                   kclust(2,i) = 0                                    *
*                   kclust(3,i) = 0                                    *
*                                                                      *
*        lclust(i,nclust)                                              *
*                                                                      *
*                i = 0, angular momentum                               *
*                  = 1, proton number                                  *
*                  = 2, neutron number                                 *
*                  = 3,                                                *
*                  = 4,                                                *
*                  = 5, charge                                         *
*                  = 6,                                                *
*                  = 7,                                                *
*                  = 8, isomer level (0: Ground, 1,2: 1st, 2nd isomer) *
*                                                                      *
*        sclust(i,nclust)                                              *
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
*        kdecay(4) = 0 : no fission                                    *
*                  = 1 : with fission                                  *
*                                                                      *
*                                                                      *
************************************************************************
      use NGSDATAMOD, only : bindeg

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------
*     common for output of evaporation and fission
*-----------------------------------------------------------------------

      include 'param00.inc'

      common /clusts/ nclust, kclust(3,nnn)
!$OMP THREADPRIVATE(/clusts/)
      common /clustu/ lclust(0:8,nnn), sclust(0:12,nnn)
!$OMP THREADPRIVATE(/clustu/)
      common /clustv/ kdecay(4)
!$OMP THREADPRIVATE(/clustv/)

*-----------------------------------------------------------------------

      common /qparm/  ielas,icasc,iqstep,lvlopt,igamma
      common /preeq/  npcle, nhole, efermi, atar, ztar
!$OMP THREADPRIVATE(/preeq/)

*-----------------------------------------------------------------------

      dimension erffg(2),affg(2),zffg(2),exffg(2)
      dimension erfrg(2),afrg(2),zfrg(2),exfrg(2)
      dimension alp0(2),bet0(2),gam0(2),ptt0(2)

      dimension izp(6), inp(6), imp(6)

      dimension npart(6),epart(100,2),hepart(100,4)
      dimension cosq2(3,100,6)

*-----------------------------------------------------------------------

      data izp/0,1,1,1,2,2/
      data inp/1,0,1,2,1,2/
      data imp/1,1,2,3,3,4/
      data rms/931.504d0/

*-----------------------------------------------------------------------
*        initial values
*-----------------------------------------------------------------------

               apr  = dble( iz + in )
               zpr  = dble( iz )
               erec = ( et - rm ) * 1000.0

*-----------------------------------------------------------------------
*        first mather
*-----------------------------------------------------------------------

               ifssev   = 1

               afrg(1)  = apr
               zfrg(1)  = zpr
               exfrg(1) = ex
               erfrg(1) = erec

            if( pt .gt. 0.0d0 ) then

               alp0(1)  = px / pt
               bet0(1)  = py / pt
               gam0(1)  = pz / pt

               ptt0(1)  = pt

            else

               alp0(1)  = 0.0
               bet0(1)  = 0.0
               gam0(1)  = 1.0

               ptt0(1)  = 0.0

            end if

*-----------------------------------------------------------------------
*        fission ( ifssev = 2 )
*-----------------------------------------------------------------------

               call fissin(ifssev,apr,zpr,ex,px,py,pz,
     &                     afrg,zfrg,exfrg,erfrg,
     &                     alp0,bet0,gam0,ptt0)

*-----------------------------------------------------------------------
*        zero set for booking
*-----------------------------------------------------------------------

               nclust = 0

*-----------------------------------------------------------------------
*        evaporation
*-----------------------------------------------------------------------

         do k = 1, ifssev

               erec = erfrg(k)
               apr  = afrg(k)
               zpr  = zfrg(k)
               ex   = exfrg(k)

               alp00 = alp0(k)
               bet00 = bet0(k)
               gam00 = gam0(k)

*-----------------------------------------------------------------------
*        subroutine ERUP
*-----------------------------------------------------------------------

               call erup(iqstep,lvlopt,
     &                   apr,zpr,ex,erec,alp00,bet00,gam00,
     &                   atar,ztar,npcle,nhole,efermi,
     &                   npart,epart,hepart,cosq2)

*-----------------------------------------------------------------------
*        booking of evapolration particles
*-----------------------------------------------------------------------

                        pxp = 0.0
                        pyp = 0.0
                        pzp = 0.0

            do i = 1, 6

               if( npart(i) .gt. 0 ) then

                        rmsp = ( rms * imp(i) - bindeg(izp(i),inp(i)) )
     &                       / 1000.0

                  do j = 1, npart(i)

                        nclust = nclust + 1

                        kclust(1,nclust) = 101
                        kclust(2,nclust) = 0
                        kclust(3,nclust) = 0

                     if( i .le. 2 ) then

                        ekin = epart(j,i) / 1000.0

                     else if( i .gt. 2 ) then

                        ekin = hepart(j,i-2) / 1000.0

                     end if


                        pabs = sqrt( ekin**2 + 2.0 * rmsp * ekin )

                        lclust(0,nclust)  = 0
                        lclust(1,nclust)  = izp(i)
                        lclust(2,nclust)  = inp(i)
                        lclust(3,nclust)  = 0
                        lclust(4,nclust)  = 0
                        lclust(5,nclust)  = izp(i)
                        lclust(6,nclust)  = 0
                        lclust(7,nclust)  = 0
                        lclust(8,nclust)  = 0

                        sclust(0,nclust)  = 0.0
                        sclust(1,nclust)  = pabs * cosq2(1,j,i)
                        sclust(2,nclust)  = pabs * cosq2(2,j,i)
                        sclust(3,nclust)  = pabs * cosq2(3,j,i)
                        sclust(4,nclust)  = ekin + rmsp
                        sclust(5,nclust)  = rmsp
                        sclust(6,nclust)  = 0.0
                        sclust(7,nclust)  = ekin * 1000.
                        sclust(8,nclust)  = wt
                        sclust(9,nclust)  = 0.0
                        sclust(10,nclust) = 0.0d0
                        sclust(11,nclust) = 0.0d0
                        sclust(12,nclust) = 0.0d0

                        pxp = pxp + sclust(1,nclust)
                        pyp = pyp + sclust(2,nclust)
                        pzp = pzp + sclust(3,nclust)

                  end do

               end if

            end do

*-----------------------------------------------------------------------
*        booking of residual nucleus
*           ( momentum is not conserved ?? , we have to use erec )
*-----------------------------------------------------------------------

            if( apr .gt. 0.0 .and. apr .ge. zpr ) then

                        nclust = nclust + 1

                        kclust(1,nclust) = 101
                        kclust(2,nclust) = 0
                        kclust(3,nclust) = 0

                        rmsp = ( rms * apr
     &                         - bindeg(nint(zpr),nint(apr-zpr)) )
     &                       / 1000.0

                        ekin = max( 0.0d0, erec / 1000.0 )
                        pabs = sqrt( ekin**2 + 2.0 * rmsp * ekin )
                        etot = sqrt( pabs**2 + rmsp**2 )

                        pab0 = sqrt( pxp**2 + pyp**2 + pzp**2 )

                     if( pab0 .gt. 0.0 ) then

                        pxp = - pabs * pxp / pab0
                        pyp = - pabs * pyp / pab0
                        pzp = - pabs * pzp / pab0

                     else

                        pxp = pabs * alp00
                        pyp = pabs * bet00
                        pzp = pabs * gam00

                     end if

                        lclust(0,nclust)  = 0
                        lclust(1,nclust)  = nint( zpr )
                        lclust(2,nclust)  = nint( apr - zpr )
                        lclust(3,nclust)  = 0
                        lclust(4,nclust)  = 0
                        lclust(5,nclust)  = nint( zpr )
                        lclust(6,nclust)  = 0
                        lclust(7,nclust)  = 0
                        lclust(8,nclust)  = 0

                        sclust(0,nclust)  = 0.0
                        sclust(1,nclust)  = pxp
                        sclust(2,nclust)  = pyp
                        sclust(3,nclust)  = pzp
                        sclust(4,nclust)  = etot
                        sclust(5,nclust)  = rmsp
                        sclust(6,nclust)  = ex
                        sclust(7,nclust)  = ekin * 1000.
                        sclust(8,nclust)  = wt
                        sclust(9,nclust)  = 0.0
                        sclust(10,nclust) = 0.0d0
                        sclust(11,nclust) = 0.0d0
                        sclust(12,nclust) = 0.0d0

            end if

*-----------------------------------------------------------------------

         end do

*-----------------------------------------------------------------------
*        fission frag
*-----------------------------------------------------------------------

               kdecay(4) = ifssev - 1

*-----------------------------------------------------------------------

      return
      end

