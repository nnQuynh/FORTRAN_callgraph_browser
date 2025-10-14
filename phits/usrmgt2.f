************************************************************************
*                                                                      *
      subroutine usrmgt2(b_mag,ctime,ptime,o_mag)
*                                                                      *
*        sample subroutine for user defined time dependence of         *
*        magnetic field.                                               *
*                                                                      *
*        input :                                                       *
*           b_mag : initial magnetic field                             *
*           ctime : critical time (msec) given by input                *
*           ptime : present time (msec)                                *
*                                                                      *
*        output :                                                      *
*           o_mag : final magnetic field                               *
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'

*-----------------------------------------------------------------------

      common /verjam/ versn, lastr, iyeav, imonv, idayv

*-----------------------------------------------------------------------

      character filename*100
      data iout / 0 /
      character rpa*1
      data rpa /'}'/
      character yen*1
      save iout !FURUTA
!$OMP THREADPRIVATE(iout)

*-----------------------------------------------------------------------
*     User Define
*     Time dependence of Magnetic Field
*     a : critical time, b : present time ( msec )
*-----------------------------------------------------------------------


*-----------------------------------------------------------------------

         tfunc(a,b) = 0.33 * ( b - a )
     &              * exp ( - 0.2 * ( b - a ) ) + 0.02

*-----------------------------------------------------------------------

            yen  = char(92)
            iout = iout + 1

*-----------------------------------------------------------------------
*        example for user defined time dependence of magnetic field
*-----------------------------------------------------------------------

         if( ptime .lt. ctime ) then

            o_mag = 0.0d0

         else

            o_mag = b_mag * tfunc( ctime, ptime )

         end if

*-----------------------------------------------------------------------
*     output of time dependence of magnetic field
*-----------------------------------------------------------------------

      if( iout .eq. 1 ) then

*-----------------------------------------------------------------------

         tini = 0.0
         tfin = 50.0
         ntms = 100

         filename = 'tdep01.ang'

*-----------------------------------------------------------------------

         iot = 31
         open(iot, file = filename, status = 'unknown' )

*-----------------------------------------------------------------------

         write(iot,'(''# Time Dependent of Magnetic Field'')')

               write(iot,'(/''msuc: {'',a1,
     &             ''huge Time Dep. of Mag. Field'',a1)') yen, rpa

               write(iot,'(''msdl: {'',a1,''it plotted by '',
     &                  a1,''ANGEL '',a1,''version}'')') yen, yen, yen
               write(iot,'(''msdr: {'',a1,''it calculated by '',
     &                  a1,''PHITS '',f5.2''}'')') yen, yen, versn

               write(iot,'(/''p: afac(0.8) form(0.9)'')')

               write(iot,'(/''x: Time [msec]'')')
               write(iot,'( ''y: G / G_0'')')

               write(iot,'(/''h: x  y,l0'')')

               tdef = ( tfin - tini ) / dble( ntms )

            do i = 1, ntms + 1

               tvar = tini + dble( i - 1 ) * tdef

               if( tvar .lt. ctime ) then

                  smag = 0.0d0

               else

                  smag = tfunc( ctime, tvar )

               end if

               write(iot,'(1p2e13.4)') tvar, smag

            end do

*-----------------------------------------------------------------------

         idasa = mmmax + 1

         open(iot, file = filename, status = 'unknown' )
         call a_angel(idasa,filename)

      end if

*-----------------------------------------------------------------------

      return
      end

