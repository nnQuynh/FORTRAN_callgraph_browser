************************************************************************
*                                                                      *
      subroutine anal_int(elabin,iprunin)
*                                                                      *
*        anal-011 : dump information of nuclei after QMD               *
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

      common /analevt/ irunp

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

         call sm_mini

*-----------------------------------------------------------------------

         open(89, file = 'dumpqmd.dat',
     &                     form='formatted',status = 'unknown' )

         write(89,'(''  IZ   N     Px(MeV/c/A) Py'',
     &              ''    Pz     Ex(MeV)   L x y z  BIMP(fm)'',
     &              ''    ISIM'')')

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
*        nucleus
*-----------------------------------------------------------------------

            if( ik .eq. 0 ) then

                  call sm_mass(1,iz,in,qmdfac)

                  pxm = px * 1000.0 / dble( iz + in )
                  pym = py * 1000.0 / dble( iz + in )
                  pzm = pz * 1000.0 / dble( iz + in )

              write(89,'(2i4,3f9.2,f10.2,2x,4i2,f9.3,i8)')
     &           iz,in,pxm,pym,pzm,ex,0,0,0,0,bi,irunp

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
*        nucleus
*-----------------------------------------------------------------------

            if( ik .eq. 0 ) then

                  call sm_mass(2,iz,in,sdmfac)

               if( kdecay(4) .eq. 1 ) then

                  call sm_mass(3,iz,in,sdmfac)

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

         call sm_norm(fac)


         open(90, file = 'masdisqmd.dat',
     &                     form='formatted',status = 'unknown' )

         open(91, file = 'dchainqmd.dat',
     &                     form='formatted',status = 'unknown' )

         call sm_mdis(90,30)
         call pr_mdis(91,30)

         open(92, file = 'masdisgem.dat',
     &                     form='formatted',status = 'unknown' )

         open(93, file = 'dchaingem.dat',
     &                     form='formatted',status = 'unknown' )

         call sm_mdis(92,31)
         call pr_mdis(93,31)


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
