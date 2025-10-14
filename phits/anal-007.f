************************************************************************
*                                                                      *
      subroutine anal_int(elabin,iprunin)
*                                                                      *
*        anal-007 : for recoil energy                                  *
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

      dimension pxstt(10)
      save pxstt
      save eccm

*-----------------------------------------------------------------------

         elab  = elabin
         iprun = iprunin

         qmdfac = 1.0d0
         sdmfac = 1.0d0

*-----------------------------------------------------------------------
*        CM energy
*-----------------------------------------------------------------------

         n1 = masspr
         n2 = massta

         prmas  = rmass
         tamas  = rmass

         plab = sqrt( elab * ( 2.0 * prmas + elab ) )

         ptot = plab * n1
         etot = elab * n1 + prmas * n1 + tamas * n2
         stot = sqrt( etot**2 - ptot**2 )
         eccm = stot - ( prmas * n1 + tamas * n2 )

*-----------------------------------------------------------------------

         pxstt(1) = 0.0
         pxstt(2) = 0.0
         pxstt(3) = 0.0
         pxstt(4) = 0.0
         pxstt(5) = 0.0
         pxstt(6) = 0.0

*-----------------------------------------------------------------------

      return


************************************************************************
*                                                                      *
      entry anal_qmd(ik,jj,iz,in,id,is,ic,iq,im,
     &               bi,px,py,pz,et,rm,ex,ek,we)
*                                                                      *
*              user entry for analysis of QMD                          *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------

         if( ik .eq. 0 ) then

               pxstt(6) = pxstt(6) + et - rm
               pxstt(7) = pxstt(7) + ex

         end if

*-----------------------------------------------------------------------

         return

************************************************************************
*                                                                      *
      entry anal_sdm(ik,jj,iz,in,id,is,ic,iq,im,
     &               bi,px,py,pz,et,rm,ex,ek,we)
*                                                                      *
*        Purpose:                                                      *
*              user entry for analysis of SDM                          *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------

         if( ik .eq. 0 ) then

            if( iz .eq. 1 .and. in .eq. 1 ) then

               kk = 1

            else if( iz .eq. 1 .and. in .eq. 2 ) then

               kk = 2

            else if( iz .eq. 2 .and. in .eq. 1 ) then

               kk = 3

            else if( iz .eq. 2 .and. in .eq. 2 ) then

               kk = 4

            else if( iz .gt. 0 .and. in .gt. 0 ) then

               kk = 5

            end if

               pxstt(kk) = pxstt(kk) + et - rm

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

         fac = 1.0 / dble(iprun) * 1000.0

         pxstt(1) = pxstt(1) * fac
         pxstt(2) = pxstt(2) * fac
         pxstt(3) = pxstt(3) * fac
         pxstt(4) = pxstt(4) * fac
         pxstt(5) = pxstt(5) * fac
         pxstt(6) = pxstt(6) * fac
         pxstt(7) = pxstt(7) * fac * 0.0001

         pxstt(2) = pxstt(1) + pxstt(2)
         pxstt(3) = pxstt(2) + pxstt(3)
         pxstt(4) = pxstt(3) + pxstt(4)
         pxstt(5) = pxstt(4) + pxstt(5)

         write(io,'(10(1pe11.4))') elab*1000.0, ( pxstt(i), i = 1, 7 )

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
