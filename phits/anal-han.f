************************************************************************
*                                                                      *
      subroutine anal_int(elabin,iprunin)
*                                                                      *
*        anal-han : p+Al at 180MeV reaction                            *
*                   mass distribution                                  *
*                   and ddx for A=22 and angular distribution          *
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

            nw   = 9
            ipx  = 1

            wi(1) = cos(   0.0 * pi / 180.0 )
            wt(1) = cos(  15.0 * pi / 180.0 )

            wi(2) = cos(  15.0 * pi / 180.0 )
            wt(2) = cos(  25.0 * pi / 180.0 )

            wi(3) = cos(  25.0 * pi / 180.0 )
            wt(3) = cos(  35.0 * pi / 180.0 )

            wi(4) = cos(  35.5 * pi / 180.0 )
            wt(4) = cos(  45.5 * pi / 180.0 )

            wi(5) = cos(  45.0 * pi / 180.0 )
            wt(5) = cos(  55.0 * pi / 180.0 )

            wi(6) = cos(  55.0 * pi / 180.0 )
            wt(6) = cos(  65.0 * pi / 180.0 )

            wi(7) = cos(  65.0 * pi / 180.0 )
            wt(7) = cos(  75.0 * pi / 180.0 )

            wi(8) = cos(  75.0 * pi / 180.0 )
            wt(8) = cos(  85.0 * pi / 180.0 )

            wi(9) = cos(  85.0 * pi / 180.0 )
            wt(9) = cos(  95.0 * pi / 180.0 )


*-----------------------------------------------------------------------

            inum = 0

            ifac = 1
            tfac = 1.0

            ilog = 0

            nx   = 25
            xmin = 0.0
            xmax = 50.0

*-----------------------------------------------------------------------

            call jbook1(1,'A=24 DDX',
     &                  tfac,ilog,inum,ifac,nx,xmin,xmax,ipx,nw,wi,wt)

            call jbook1(2,'A=22 DDX',
     &                  tfac,ilog,inum,ifac,nx,xmin,xmax,ipx,nw,wi,wt)

            call jbook1(3,'A=16 DDX',
     &                  tfac,ilog,inum,ifac,nx,xmin,xmax,ipx,nw,wi,wt)

            call jbook1(4,'A=12 DDX',
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

                  ekin = ( et - rm ) * 1000.0

                  pl2 = px**2 + py**2 + pz**2
                  pla = sqrt(pl2)

                  if( pl2 .eq. 0 ) then
                     cosa = 1.0
                  else
                     cosa = pz / sqrt(pl2)
                  end if

               if( iz+in .eq. 24 ) then

                  call jfill1(1,ekin,cosa,sdmfac,1.0/2.0/pi*qmdfac)

               else if( iz+in .eq. 22 ) then

                  call jfill1(2,ekin,cosa,sdmfac,1.0/2.0/pi*qmdfac)

               else if( iz+in .eq. 16 ) then

                  call jfill1(3,ekin,cosa,sdmfac,1.0/2.0/pi*qmdfac)

               else if( iz+in .eq. 12 ) then

                  call jfill1(4,ekin,cosa,sdmfac,1.0/2.0/pi*qmdfac)

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
         call jscale1(3,fac)
         call jscale1(4,fac)

         call jprint1(1,90)
         call jprint1(2,91)
         call jprint1(3,92)
         call jprint1(4,93)

         call sm_norm(fac)
         call sm_mdis(95,31)
         call pr_mdis(96,31)


*-----------------------------------------------------------------------

      return
      end

