************************************************************************
*                                                                      *
      subroutine usrelst2(icnt,icels,ioel,erg,costh,sigtt)
*                                                                      *
*        sample subroutine for user defined elastic option             *
*                                                                      *
*        input :                                                       *
*           icnt  : =1 total x section =2 elastic costh                *
*           icels : =0 this els, =1 back to normal elastic             *
*           ioel  : region id                                          *
*           erg   : neutron energy in MeV                              *
*                                                                      *
*        output :                                                      *
*           costh : cos(the) of outgoing angle in CM system            *
*           sigtt : total cross section (barn)                         *
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'

      parameter ( pi  = 3.1415926535898d0 )

*-----------------------------------------------------------------------

      common /verjam/ versn, lastr, iyeav, imonv, idayv
      common /elreg/  elarg(kvlmax,4), mlrgn, melrg(kvlmax,2), ielusr
      common /regdc/  idrg(kvlmax), idgr(kvmmax)
      common /paraj/ mstz(300), parz(300)

      dimension       dsiga(1000), dsigp(1000)
      save            dsiga, dsigp
!$OMP THREADPRIVATE(dsiga, dsigp)
*-----------------------------------------------------------------------

      character filename*100
      data iout / 0 /
      character rpa*1
      data rpa /'}'/
      character yen*1

*-----------------------------------------------------------------------

               icels = 0
               ssr = 0.d0

*-----------------------------------------------------------------------
*     User Define angular distribution for low energy neutron
*     elastic scattering in CM system.
*-----------------------------------------------------------------------

               rlam = 2.86e-2 / sqrt( erg * 1.e+6 )   ! nm

               c1 = elarg(ioel,1)
               c2 = elarg(ioel,2)
               c3 = elarg(ioel,3)
               c4 = elarg(ioel,4)

*-----------------------------------------------------------------------
*     total cross section by data
*-----------------------------------------------------------------------

         if( icnt .ne. 1 ) then

*-----------------------------------------------------------------------

            return

         end if

*-----------------------------------------------------------------------
*     angular distribution
*-----------------------------------------------------------------------

               ianm  = mstz(24) * 2
               andif = 180.0 / ianm

*-----------------------------------------------------------------------

               sek = 0.0

            do ian = 1, ianm

               aplow = dble(ian-1) * andif
               aphig = dble(ian-1) * andif + andif
               apoin = ( aplow + aphig ) /2.0
               aprad = apoin * pi / 180.0
               apsin = sin( aprad )
               adrad = ( aphig - aplow ) * pi / 180.0

               call elasang(rlam,c1,c2,c3,c4,aprad,dsigs)

               dsigp(ian) = aplow

               dsiga(ian) = dsigs * 2. * pi * apsin * adrad

               sek = sek + dsiga(ian)

            end do

               fnorm = sek

               sek = 0.0

            do ian = 1, ianm

               sek = sek + dsiga(ian) / fnorm

               dsiga(ian) = sek

            end do

*-----------------------------------------------------------------------
*        random number 0 < ramx < 1
*-----------------------------------------------------------------------

               ram1 = unirn(dummy)
               ram2 = unirn(dummy)

*-----------------------------------------------------------------------

            do ian = 1, ianm

               if( dsiga(ian) .gt. ram1 ) goto 100

            end do

  100       continue

               angcm = ( dsigp(ian) + ram2 * andif ) * pi / 180.0

               costh = cos( angcm )

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      subroutine elasang(rlam,c1,c2,c3,c4,aprad,dsigs)
*                                                                      *
*        sample subroutine for user defined elastic angular dist.      *
*                                                                      *
*        input :                                                       *
*           rlam : lambda (nm)                                         *
*           c1,c2,c3,c4  : user defined constants                      *
*           aprad : angle (rad)                                        *
*                                                                      *
*        output :                                                      *
*           disgs : probability                                        *
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      parameter ( pi  = 3.1415926535898d0 )

*-----------------------------------------------------------------------

         dq = 4.d0 * pi * sin( aprad / 2.d0 ) / rlam

         dsigs = 1.d0 / ( dq**2 + c1**2 )**2

*-----------------------------------------------------------------------

      return
      end

