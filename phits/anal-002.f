************************************************************************
*                                                                      *
      subroutine anal_int(elabin,iprunin)
*                                                                      *
*        anal-002 : Ishibashi KEK experiments and Meiyer data          *
*                   and mass distribution                              *
*                   and excitation energy distribution of residual     *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              user subroutine for analysis                            *
*                                                                      *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------

c     implicit double precision(a-h, o-z)

c     above statement is already included in jam1.inc
*-----------------------------------------------------------------------

      include 'jam1.inc'
      include 'jam2.inc'
      include 'param02.inc'

*-----------------------------------------------------------------------

      common /vriab3/ qmdfac, sdmfac
!$OMP THREADPRIVATE(/vriab3/)
      common /swich3/ ielst, jelst, kelst
!$OMP THREADPRIVATE(/swich3/)

      common /clustv/ kdecay(4)
!$OMP THREADPRIVATE(/clustv/)

*-----------------------------------------------------------------------

      dimension wi(20)
      dimension wt(20)

      save elab, iprun

*-----------------------------------------------------------------------

         elab  = elabin
         iprun = iprunin

         qmdfac = 1.0d0
         sdmfac = 1.0d0

*-----------------------------------------------------------------------
*     initialize mass distribution
*-----------------------------------------------------------------------


*-----------------------------------------------------------------------
*     initialize jbook
*-----------------------------------------------------------------------

            nw   = 7
            ipx  = 1


            wi(1) = cos(   5.0 * pi / 180.0 )
            wt(1) = cos(  10.0 * pi / 180.0 )

            wi(2) = cos(  12.5 * pi / 180.0 )
            wt(2) = cos(  17.5 * pi / 180.0 )

            wi(3) = cos(  27.5 * pi / 180.0 )
            wt(3) = cos(  32.5 * pi / 180.0 )

            wi(4) = cos(  57.5 * pi / 180.0 )
            wt(4) = cos(  62.5 * pi / 180.0 )

            wi(5) = cos(  87.5 * pi / 180.0 )
            wt(5) = cos(  92.5 * pi / 180.0 )

            wi(6) = cos( 117.5 * pi / 180.0 )
            wt(6) = cos( 122.5 * pi / 180.0 )

            wi(7) = cos( 147.5 * pi / 180.0 )
            wt(7) = cos( 152.5 * pi / 180.0 )

*-----------------------------------------------------------------------

            inum = 0

            ifac = 1
            tfac = 1.0

            ilog = 0
            ix   = 30

            xmin  = 0.0
            xmax0 = elab * 1000.0
            xmax1 = xmax0 * 1.3

            dx    = ( xmax0 - xmin ) / dble( ix )
            nx    = nint( ( xmax1 - xmin ) / dx )
            xmax  = xmin + dx * dble( nx )

*-----------------------------------------------------------------------

            call jbook1(1,'neutron EVAP linear',
     &                  tfac,ilog,inum,ifac,nx,xmin,xmax,ipx,nw,wi,wt)

*-----------------------------------------------------------------------

            inum = 0

            ifac = 1
            tfac = 1.0

            ilog = 1
            ix   = 30

            xmin = 1.0
            xmax0 = elab * 1000.0
            xmax1 = xmax0 * 3.0

            dx    = log( xmax0 / xmin ) / dble( ix )
            nx    = nint( log( xmax1 / xmin ) / dx )
            xmax  = xmin + exp( dx * dble( nx ) )

*-----------------------------------------------------------------------

            call jbook1(2,'neutron EVAP log',
     &                  tfac,ilog,inum,ifac,nx,xmin,xmax,ipx,nw,wi,wt)

*-----------------------------------------------------------------------

      return


************************************************************************
*                                                                      *
      entry anal_qmd(ik,jj,iz,in,id,is,ic,iq,im,
     &               bi,px,py,pz,et,rm,ex,ek,we)
*                                                                      *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              user entry for analysis of QMD                          *
*                                                                      *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------
*        inelastic frag
*-----------------------------------------------------------------------



*-----------------------------------------------------------------------
*        nucleus
*-----------------------------------------------------------------------

            if( ik .eq. 0 ) then

                  excit = ex / dble( iz + in )

               if( excit .gt. 0.0 ) then

                  gfac = 1.0


               end if

                  call sm_mass(1,iz,in,qmdfac)


            end if

*-----------------------------------------------------------------------
*        proton
*-----------------------------------------------------------------------

            if( ik .eq. 1 ) then

                  ekin = ( et - rm ) * 1000.0

                  pl2 = px**2 + py**2 + pz**2
                  pla = sqrt(pl2)

                  if( pl2 .eq. 0 ) then
                     cosa = 1.0
                  else
                     cosa = pz / sqrt(pl2)
                  end if

           end if

*-----------------------------------------------------------------------
*        neutron
*-----------------------------------------------------------------------

            if( ik .eq. 2 ) then

                  ekin = ( et - rm ) * 1000.0

                  pl2 = px**2 + py**2 + pz**2
                  pla = sqrt(pl2)

                  if( pl2 .eq. 0 ) then
                     cosa = 1.0
                  else
                     cosa = pz / sqrt(pl2)
                  end if


            end if

*-----------------------------------------------------------------------
*        pions
*-----------------------------------------------------------------------

            if( ik .eq. 3 ) then

                  ekin = ( et - rm ) * 1000.0

                  pl2 = px**2 + py**2 + pz**2
                  pla = sqrt(pl2)

                  if( pl2 .eq. 0 ) then
                     cosa = 1.0
                  else
                     cosa = pz / sqrt(pl2)
                  end if

               if( ic .eq. 1 ) then

               else if( ic .eq. 0 ) then

               else if( ic .eq. -1 ) then

               end if

            end if


*-----------------------------------------------------------------------
*        Kaon
*-----------------------------------------------------------------------

            if( ik .eq. 5 ) then

                  rap = 0.5d0 * log( max(et+pz,1.d-8)
     &                             / max(et-pz,1.d-8) )

                  ekin = ( et - rm ) * 1000.0

                  ptsq = px**2 + py**2
                  pt   = sqrt( max(ptsq,1.d-8) )
                  emt  = sqrt( pt**2 + rm )
                  emt0 = emt - rm

               if( ic .eq. 1 ) then

               else if( ic .eq. 0 ) then

               else if( ic .eq. -1 ) then

               end if

            end if


*-----------------------------------------------------------------------

      return


************************************************************************
*                                                                      *
      entry anal_sdm(ik,jj,iz,in,id,is,ic,iq,im,
     &               bi,px,py,pz,et,rm,ex,ek,we)
*                                                                      *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              user entry for analysis of SDM                          *
*                                                                      *
*                                                                      *
************************************************************************


*-----------------------------------------------------------------------
*        inelastic frag
*-----------------------------------------------------------------------



*-----------------------------------------------------------------------
*        nucleus
*-----------------------------------------------------------------------

            if( ik .eq. 0 ) then

                  call sm_mass(2,iz,in,sdmfac)

               if( kdecay(4) .eq. 1 ) then

                  call sm_mass(3,iz,in,sdmfac)

               end if

               if( iz .eq. 3 .and. in .eq. 4 ) then

                  ekin = ( et - rm ) * 1000.0

                  pl2 = px**2 + py**2 + pz**2
                  pla = sqrt(pl2)

                  if( pl2 .eq. 0 ) then
                     cosa = 1.0
                  else
                     cosa = pz / sqrt(pl2)
                  end if


               end if

            end if

*-----------------------------------------------------------------------
*        proton
*-----------------------------------------------------------------------

            if( ik .eq. 1 ) then

                  ekin = ( et - rm ) * 1000.0

                  pl2 = px**2 + py**2 + pz**2
                  pla = sqrt(pl2)

                  if( pl2 .eq. 0 ) then
                     cosa = 1.0
                  else
                     cosa = pz / sqrt(pl2)
                  end if


            end if

*-----------------------------------------------------------------------
*        neutron
*-----------------------------------------------------------------------

            if( ik .eq. 2 ) then

                  ekin = ( et - rm ) * 1000.0

                  pl2 = px**2 + py**2 + pz**2
                  pla = sqrt(pl2)

                  if( pl2 .eq. 0 ) then
                     cosa = 1.0
                  else
                     cosa = pz / sqrt(pl2)
                  end if

                  call jfill1(1,ekin,cosa,sdmfac,1.0/2.0/pi*qmdfac)
                  call jfill1(2,ekin,cosa,sdmfac,1.0/2.0/pi*qmdfac)


            end if

*-----------------------------------------------------------------------

      return


************************************************************************
*                                                                      *
      entry anal_fin( io, bmax0 )
*                                                                      *
*                                                                      *
*        Purpose:                                                      *
*                                                                      *
*              user entry for analysis of weight                       *
*                                                                      *
*                                                                      *
************************************************************************


*-----------------------------------------------------------------------

         fac = 1.0 / dble(iprun) * pi * bmax0**2 * 10.0

         call jscale1(1,fac)
         call jscale1(2,fac)

         call jprint1(1,90)
         call jprint1(2,91)

         call sm_norm(fac)
         call sm_mdis(94,31)
         call pr_mdis(98,31)


*-----------------------------------------------------------------------
      return

************************************************************************
*                                                                      *
       entry anal_fin2( io, bmax0,itdpa,iteth,ap,zp )
*       entry anal_fin2( io, bmax0 )
*                                                                      *
*                                                                      *
*        Purpose:         dummy for anal-008.f                         *
*                                                                      *
*              DPA cross section                                       *
*                                                                      *
*       2021/06/23 Y. Iwamoto                                          *
************************************************************************

*-----------------------------------------------------------------------

*-----------------------------------------------------------------------

      return
      end
