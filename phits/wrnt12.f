************************************************************************
*                                                                      *
      subroutine wrnt12(data)
*                                                                      *
*       output the neutron information to the neutron cut-off file.    *
*       modified by K.Niita on 26/03/2000                              *
*                                                                      *
*     input : structure of data(i) array                               *
*                                                                      *
*       i=1 : particle energy [MeV]                                    *
*         2 : x position [cm]                                          *
*         3 : y position [cm]                                          *
*         4 : z position [cm]                                          *
*         5 : directional cosine of x                                  *
*         6 : directional cosine of y                                  *
*         7 : directional cosine of z                                  *
*         8 : particle weight                                          *
*         9 : particle type (1=neutron, 2=proton, 3=photon)            *
*        10 : region number                                            *
*        11 : time [ns]                                                *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'

*-----------------------------------------------------------------------

      common /ngcut/  incut, igcut, ipcut

      common /wtcntl/ iwt, icimp(20), ifcls(20), iwwin(20), ircls(20)

      dimension bufr(3000),data(11)

      data nn, nd / 0, 0 /
      data ncall / 0 /
      data rflag / 1.0d0 /

      save bufr,ni,nd,nn
      save ncall, rflag

*-----------------------------------------------------------------------

         if( incut .eq. 0 ) return

*-----------------------------------------------------------------------

         if( incut .eq. 2 .and. ncall .eq. 0 ) then

            write(12) real(rflag), real(rflag), real(rflag),
     &                real(rflag), real(rflag)

            ncall = 1

         end if

*-----------------------------------------------------------------------
*        end of file 12 ( data(1) < 0.0 )
*-----------------------------------------------------------------------

         if( data(1) .lt. 0.0 ) then

            if( nd .gt. 0 ) then

               rd = dble(nd)
               if( iwt .eq. 0 ) rd = -rd

               write(12) real(rd),real(nn),(real(bufr(i)),i=1,nd)

            end if

               endfile 12

               return

         end if

*-----------------------------------------------------------------------
*        only neutron
*-----------------------------------------------------------------------

               n = nint(data(9))

               if( n .ne. 1 ) return

*-----------------------------------------------------------------------
*     check values
*-----------------------------------------------------------------------

         do i = 1, 10

            if( abs( data(i) ) .gt. 1.d+61 ) then

               write(6,'(/''*** Error in wrnt12, data('',i2,
     &         '') is greater than 1d+61; skip this particles'')') i

               return

            end if

         end do

*-----------------------------------------------------------------------
*     check positions
*-----------------------------------------------------------------------

         if( nn .eq. 0 ) then

                  ipost = 1

         else if( nn .gt. 0 ) then

                  ipost = 0

            do i = 1, 3

               if( data(i+1) .eq. 0.0d0 ) then

                  if( bufr(ni+i) .ne. 0.0d0 ) ipost = ipost + 1

               else

                  r = abs( bufr(ni+i) / data(i+1) - 1.d+0 )
                  if( r .gt. 1.0d-4 ) ipost = ipost + 1

               end if

                  if( ipost .gt. 0 ) goto 100

            end do

  100       continue

         end if

*-----------------------------------------------------------------------
*     the same position as previous one and nd =< 2000
*-----------------------------------------------------------------------

         if( ipost .eq. 0 .and. nd .le. 2000 ) then

               bufr(ni)   = bufr(ni) + 1.0d0
               bufr(nd+1) = data(1)
               bufr(nd+2) = data(5)
               bufr(nd+3) = data(6)
               bufr(nd+4) = data(7)

               nd = nd + 4

            if( iwt .ne. 0 ) then

               bufr(nd+1) = data(8)
               nd = nd + 1

            end if

            if( incut .eq. 2 ) then

               bufr(nd+1) = data(11)
               nd = nd + 1

            end if

*-----------------------------------------------------------------------
*     new position :
*         write information on buffer
*         write information on file 12 if nd > 2000
*-----------------------------------------------------------------------

         else

            if( nd .gt. 2000 ) then

               rd = dble(nd)
               if( iwt .eq. 0 ) rd = -rd

               write(12) real(rd),real(nn),(real(bufr(i)),i=1,nd)

               nn = 0
               nd = 0

            end if

               ni = nd + 1
               nd = nd + 1

               bufr(ni)   = 1.0d0
               bufr(nd+1) = data(2)
               bufr(nd+2) = data(3)
               bufr(nd+3) = data(4)
               bufr(nd+4) = data(1)
               bufr(nd+5) = data(5)
               bufr(nd+6) = data(6)
               bufr(nd+7) = data(7)

               nd = nd + 7

            if( iwt .ne. 0 ) then

               bufr(nd+1) = data(8)
               nd = nd + 1

            end if

            if( incut .eq. 2 ) then

               bufr(nd+1) = data(11)
               nd = nd + 1

            end if

         end if

*-----------------------------------------------------------------------

            nn = nn + 1

*-----------------------------------------------------------------------

      return
      end


