************************************************************************
*                                                                      *
      subroutine anal_int(elabin,iprunin)
*                                                                      *
*        anal-004 : Nuclear Data for Konno                             *
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
      include 'param01.inc'
      include 'param02.inc'

*-----------------------------------------------------------------------

      common /vriab3/ qmdfac, sdmfac
!$OMP THREADPRIVATE(/vriab3/)
      common /swich3/ ielst, jelst, kelst
!$OMP THREADPRIVATE(/swich3/)

      common /clustv/ kdecay(4)
!$OMP THREADPRIVATE(/clustv/)

*-----------------------------------------------------------------------

      common /comps1/ mstapr, massta, msprpr, masspr
      common /comps2/ sigela, signon, fissx

      common /summas/ sumas(3,0:maxpt,0:maxnt)

*-----------------------------------------------------------------------

      save elab, iprun

*-----------------------------------------------------------------------

      dimension wi(20), wt(20)
      dimension pxstt(10)

      dimension angvl(20)
      dimension engvl(40)

      dimension pxsdd(40)

      save angvl, engvl
      save nang, neng

*-----------------------------------------------------------------------

         elab  = elabin
         iprun = iprunin

         qmdfac = 1.0d0
         sdmfac = 1.0d0

*-----------------------------------------------------------------------
*     initialize mass distribution
*-----------------------------------------------------------------------

            call sm_mini
            call jbkreset

*-----------------------------------------------------------------------
*     initialize jbook
*-----------------------------------------------------------------------

            ipx  = 1

            inum = 0
            ifac = 1
            tfac = 1.0

            ilog = 1
            xmin = 1.0
            xmax = elab * 1000.0

            nx   = 30
            neng = nx + 2

                     engvl(1) = 0.0

               if( ilog .eq. 0 ) then

                     bmi = xmin
                     bma = xmax

               else

                     bmi = log(xmin)
                     bma = log(xmax)

               end if

                     bin = ( bma - bmi ) / nx

               do i = 1, nx + 1

                  if( ilog .eq. 0 ) then

                     engvl(i+1) = xmin + bin * ( i - 1 )

                  else

                     engvl(i+1) = xmin * exp( bin * ( i - 1 ) )

                  end if

               end do


            nw   = 19
            nang = nw

            wi(1) = cos(  0.0 * pi / 180.0 )
            wt(1) = cos(  5.0 * pi / 180.0 )

            angvl(1) = 1.0

         do ia = 2, 18

            angli = dble(ia-1) * 10.0 - 5.0
            anglt = dble(ia-1) * 10.0 + 5.0
            anglc = dble(ia-1) * 10.0

            wi(ia) = cos( angli * pi / 180.0 )
            wt(ia) = cos( anglt * pi / 180.0 )

            angvl(ia) = cos( anglc * pi / 180.0 )

         end do

            wi(19) = cos( 175.0 * pi / 180.0 )
            wt(19) = cos( 180.0 * pi / 180.0 )

            angvl(19) = -1.0

*-----------------------------------------------------------------------

            call jbook1(1,'neutron',
     &                  tfac,ilog,inum,ifac,nx,xmin,xmax,ipx,nw,wi,wt)

            call jbook1(3,'proton',
     &                  tfac,ilog,inum,ifac,nx,xmin,xmax,ipx,nw,wi,wt)

            call jbook1(4,'deuteron',
     &                  tfac,ilog,inum,ifac,nx,xmin,xmax,ipx,nw,wi,wt)

            call jbook1(5,'triton',
     &                  tfac,ilog,inum,ifac,nx,xmin,xmax,ipx,nw,wi,wt)

            call jbook1(6,'helium3',
     &                  tfac,ilog,inum,ifac,nx,xmin,xmax,ipx,nw,wi,wt)

            call jbook1(7,'alpha',
     &                  tfac,ilog,inum,ifac,nx,xmin,xmax,ipx,nw,wi,wt)

            call jbook1(8,'pi+',
     &                  tfac,ilog,inum,ifac,nx,xmin,xmax,ipx,nw,wi,wt)

            call jbook1(9,'pi0',
     &                  tfac,ilog,inum,ifac,nx,xmin,xmax,ipx,nw,wi,wt)

            call jbook1(10,'pi-',
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

         return

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

            if( ik .eq. 7 ) then

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

                  ekin = ( et - rm ) * 1000.0
                  pl2 = px**2 + py**2 + pz**2
                  if( pl2 .eq. 0 ) then
                     cosa = 1.0
                  else
                     cosa = pz / sqrt(pl2)
                  end if

            if( iz .eq. 1 .and. in .eq. 1 ) then

                  call jfill1(4,ekin,cosa,sdmfac,1.0/2.0/pi*qmdfac)

            else if( iz .eq. 1 .and. in .eq. 2 ) then

                  call jfill1(5,ekin,cosa,sdmfac,1.0/2.0/pi*qmdfac)

            else if( iz .eq. 2 .and. in .eq. 1 ) then

                  call jfill1(6,ekin,cosa,sdmfac,1.0/2.0/pi*qmdfac)

            else if( iz .eq. 2 .and. in .eq. 2 ) then

                  call jfill1(7,ekin,cosa,sdmfac,1.0/2.0/pi*qmdfac)

            end if

         end if

*-----------------------------------------------------------------------
*        nucleons
*-----------------------------------------------------------------------

         if( ik .eq. 1 .or. ik .eq. 2 ) then

                  ekin = ( et - rm ) * 1000.0
                  pl2 = px**2 + py**2 + pz**2
                  if( pl2 .eq. 0 ) then
                     cosa = 1.0
                  else
                     cosa = pz / sqrt(pl2)
                  end if

            if( ik .eq. 1 ) then

                  call jfill1(3,ekin,cosa,sdmfac,1.0/2.0/pi*qmdfac)

            else if( ik .eq. 2 ) then

                  call jfill1(1,ekin,cosa,sdmfac,1.0/2.0/pi*qmdfac)

            end if

         end if

*-----------------------------------------------------------------------
*        pions
*-----------------------------------------------------------------------

         if( ik .eq. 3 ) then

                  ekin = ( et - rm ) * 1000.0
                  pl2 = px**2 + py**2 + pz**2
                  if( pl2 .eq. 0 ) then
                     cosa = 1.0
                  else
                     cosa = pz / sqrt(pl2)
                  end if

            if( ic .eq. 1 ) then

                  call jfill1(8,ekin,cosa,sdmfac,1.0/2.0/pi*qmdfac)

            else if( ic .eq. 0 ) then

                  call jfill1(9,ekin,cosa,sdmfac,1.0/2.0/pi*qmdfac)

            else if( ic .eq. -1 ) then

                  call jfill1(10,ekin,cosa,sdmfac,1.0/2.0/pi*qmdfac)

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
         call jscale1(3,fac)
         call jscale1(4,fac)
         call jscale1(5,fac)
         call jscale1(6,fac)
         call jscale1(7,fac)
         call jscale1(8,fac)
         call jscale1(9,fac)
         call jscale1(10,fac)

         call sm_norm(fac)

*-----------------------------------------------------------------------

         write(io,'(i5)') 1

         write(io,'(4i5,1pe11.4)') mstapr, massta, msprpr, masspr,
     &                             elab*1.e+9

*-----------------------------------------------------------------------

         write(io,'(i5)') 2

         write(io,'(3(1pe11.4))')
     &              sigela*1.e-3, signon*1.e-3, fissx*1.e-3

*-----------------------------------------------------------------------

         write(io,'(i5)') 3

               call jftot(1,pxst,pxsu,pxso,2.0*pi)
            pxstt(1) = pxst + pxsu

            pxstt(2) = 0.0

               call jftot(3,pxst,pxsu,pxso,2.0*pi)
            pxstt(3) = pxst + pxsu
               call jftot(4,pxst,pxsu,pxso,2.0*pi)
            pxstt(4) = pxst + pxsu
               call jftot(5,pxst,pxsu,pxso,2.0*pi)
            pxstt(5) = pxst + pxsu
               call jftot(6,pxst,pxsu,pxso,2.0*pi)
            pxstt(6) = pxst + pxsu
               call jftot(7,pxst,pxsu,pxso,2.0*pi)
            pxstt(7) = pxst + pxsu
               call jftot(8,pxst,pxsu,pxso,2.0*pi)
            pxstt(8) = pxst + pxsu
               call jftot(9,pxst,pxsu,pxso,2.0*pi)
            pxstt(9) = pxst + pxsu
               call jftot(10,pxst,pxsu,pxso,2.0*pi)
            pxstt(10) = pxst + pxsu

         write(io,'(10(1pe11.4))') ( pxstt(i)*1.e-3, i = 1, 10 )

*-----------------------------------------------------------------------

         write(io,'(i5)') 5

            iiso = 0

            do izz = 2, maxpt
            do inn = 1, maxnt

               if( ( izz .eq. 2 .and. inn .gt. 2 ) .or.
     &             ( izz .gt. 2 .and.
     &               ( masspr .eq. 1 .and.
     &               ( izz .ne. mstapr .or.
     &                 inn .ne. massta - mstapr ) ) .or.
     &               ( masspr .gt. 1 ) ) ) then

               if( sumas(2,izz,inn) .gt. 0.0 ) then

                     iiso = iiso + 1

               end if
               end if

            end do
            end do

         write(io,'(2i5)') iiso, 0

            do izz = 2, maxpt
            do inn = 1, maxnt

               if( ( izz .eq. 2 .and. inn .gt. 2 ) .or.
     &             ( izz .gt. 2 .and.
     &               ( masspr .eq. 1 .and.
     &               ( izz .ne. mstapr .or.
     &                 inn .ne. massta - mstapr ) ) .or.
     &               ( masspr .gt. 1 ) ) ) then

               if( sumas(2,izz,inn) .gt. 0.0 ) then

                     write(io,'(2(1pe11.4),i5)')
     &               1000.0 * dble( izz ) + dble( izz + inn ),
     &               sumas(2,izz,inn)*1.e-3, 0

               end if
               end if

            end do
            end do

*-----------------------------------------------------------------------

         write(io,'(i5)') 6

         do ipk = 1, 10

            write(io,'(2i5)') nang, neng

            write(io,'(10(1pe11.4))') ( engvl(i)*1.e+6, i = 1, neng )

            do jpk = nang, 1, -1

               write(io,'(1pe11.4)') angvl(jpk)

               call jfddx(ipk,jpk,pxsdd)

               write(io,'(10(1pe11.4))') ( pxsdd(i)*1.e-9, i = 1, neng )

            end do

         end do

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
