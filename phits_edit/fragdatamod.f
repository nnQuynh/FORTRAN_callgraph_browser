************************************************************************
      module fragdatamod
*                                                                      *
*                                                                      *
      implicit real*8 (a-h,o-z)

      real*8, allocatable :: abcd(:,:,:,:,:)


*-----------------------------------------------------------------------

      contains

************************************************************************
*                                                                      *
      subroutine FragData_record(io,jo,ierr)
*                                                                      *
*       record differential cross section data of [Frag Data]          *
*       modified by S.Hashimoto on 2019/12/10                          *
*                                                                      *
************************************************************************
      use moddas_fragdata

      implicit real*8 (a-h,o-z)

      include 'param.inc'
      include 'ggsparam.inc'

*-----------------------------------------------------------------------

      common /fgdata/ ifrgd, ifgm, ifgdf(kvlmax,5), frgfl(kvlmax)
      character frgfl*200
      integer nfdopt5
      common /nfdopt5/ nfdopt5
      data nfdopt5/0/

*-----------------------------------------------------------------------

      dimension     idas(1)
      equivalence ( das, idas )


*-----------------------------------------------------------------------
C allocate

      maxneo = 1
      maxnag = 1
      maxnfrg = 1
      maxnei = 1
      do i = 1, ifrgd
       if( ifgdf(i,5) .eq. 1 ) exit
       if( ifgdf(i,1) .eq. 5 ) then
          nfdopt5 = nfdopt5 + 1

        nei = ifgs01_nei(i)
        maxnei = max(maxnei,nei)
        neo = iabs(ifgs04_neo(i))
        maxneo = max(maxneo,neo)
        nag = iabs(ifgs06_nag(i))
        maxnag = max(maxnag,nag)
        nfrg = ifgs08_nfrg(i)
        maxnfrg = max(maxnfrg,nfrg)

       end if
      end do

      if ( nfdopt5 .gt. 0 ) then
         allocate( abcd(maxneo,maxnag,maxnfrg,maxnei+1,ifrgd) )
      end if

*-----------------------------------------------------------------------
C do-loop for frag data files
      do i = 1, ifrgd
       if( ifgdf(i,5) .eq. 1 ) exit
       if( ifgdf(i,1) .eq. 5 ) then

        nei = ifgs01_nei(i)
        kne = ifgs02_kne(i)
        kxs = ifgs03_kxs(i)
        neo = ifgs04_neo(i)
        kef = ifgs05_kef(i)
        nag = ifgs06_nag(i)
        kaf = ifgs07_kaf(i)
        nfrg = ifgs08_nfrg(i)
        kim = ifgs09_kim(i)
        ks0 = ifgs10_ks0(i)

C do-loop for incident energy mesh points
        do iei = 1, nei + 1

         ksf = ifrge_ksf(ks0+iei)
         kk1 = ifrge_ks1(ks0+iei)
         kk2 = ifrge_ks2(ks0+iei)
         kk3 = ifrge_ks3(ks0+iei)

*-----------------------------------------------------------------------
C case 1
         if( neo .gt. 0 .and. nag .ne. 0 ) then

          do ifrg = 1, nfrg
           do iag = 1, iabs(nag)
            do ieo = 1, neo
               abcd(ieo,iag,ifrg,iei,i)
     &              = frgdd(kk2+(ifrg-1)*neo*iabs(nag)
     &              +(ieo-1)*iabs(nag)+iag)
            end do
           end do
          end do


C case 2
         else if ( neo .eq. 0 .and. nag .ne. 0 ) then

          ieo = 1
          do ifrg = 1, nfrg
           do iag = 1, iabs(nag)
              abcd(ieo,iag,ifrg,iei,i)
     &             = frgdx(kk3+(ifrg-1)*iabs(nag)+iag)
           end do
          end do


C case 4
         else if ( neo .gt. 0 .and. nag .eq. 0 ) then

          iag = 1
          do ifrg = 1, nfrg
           do ieo = 1, neo
              abcd(ieo,iag,ifrg,iei,i) = frgdx(kk3+(ifrg-1)*neo+ieo)
           end do
          end do


         end if

        end do ! end of do-loop for incident energy mesh points

*-----------------------------------------------------------------------

       end if
      end do ! end of do-loop for frag data files

*-----------------------------------------------------------------------

      return
      end subroutine FragData_record
*-----------------------------------------------------------------------


************************************************************************
*                                                                      *
      subroutine FragData_cumulative(io,jo,ierr)
*                                                                      *
*       calculation of cumulative distribution function of [Frag Data] *
*       modified by S.Hashimoto on 2019/12/10                          *
*                                                                      *
************************************************************************
      use moddas_fragdata

      implicit real*8 (a-h,o-z)

      include 'param.inc'
      include 'ggsparam.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      common /fgdata/ ifrgd, ifgm, ifgdf(kvlmax,5), frgfl(kvlmax)
      character frgfl*200

*-----------------------------------------------------------------------

      dimension     idas(1)
      equivalence ( das, idas )


*-----------------------------------------------------------------------
C do-loop for frag data files
      do i = 1, ifrgd
       if( ifgdf(i,5) .eq. 1 ) exit
       if( ifgdf(i,1) .gt. 0 ) then

        nei = ifgs01_nei(i)
        kne = ifgs02_kne(i)
        kxs = ifgs03_kxs(i)
        neo = ifgs04_neo(i)
        kef = ifgs05_kef(i)
        nag = ifgs06_nag(i)
        kaf = ifgs07_kaf(i)
        nfrg = ifgs08_nfrg(i)
        kim = ifgs09_kim(i)
        ks0 = ifgs10_ks0(i)

        if ( neo .gt. 0 ) then
           if ( ifgdf(i,1) .eq. 5 ) then
              nef = neo
           else
              nef = neo + 1
           end if
        else
           nef = neo
        end if
        if ( iabs(nag) .gt. 0 ) then
           if ( ifgdf(i,1) .eq. 5 ) then
              naf = iabs(nag)
           else
              naf = iabs(nag) + 1
           end if
        else
           naf = iabs(nag)
        end if

C do-loop for incident energy mesh points
        do j = 1, nei + 1

         ksf = ifrge_ksf(ks0+j)
         kk1 = ifrge_ks1(ks0+j)
         kk2 = ifrge_ks2(ks0+j)
         kk3 = ifrge_ks3(ks0+j)

*-----------------------------------------------------------------------
C case 1
         if( neo .gt. 0 .and. nag .ne. 0 ) then

          do k = 1, nfrg

           sek = 0.0d0
           do l = 1, nef-1
            do m = 1, naf-1

            if ( ifgdf(i,1) .ne. 5 ) then

               sek = sek + frgdd(kk2+(k-1)*neo*iabs(nag)
     &              +(l-1)*iabs(nag)+m)
     &              * 2.d0 * pi
     &              * ( dcos( frgaf(kaf+m) )
     &              - dcos( frgaf(kaf+m+1) ) )
     &              * ( frgef(kef+l+1) - frgef(kef+l) )

            else

              xx1 = frgef(kef+l)
              xx2 = frgef(kef+l+1)
              yy1 = dcos( frgaf(kaf+m) )
              yy2 = dcos( frgaf(kaf+m+1) )
              zz1 = frgdd(kk2+(k-1)*neo*iabs(nag)+(l-1)*iabs(nag)+m)
              zz2 = frgdd(kk2+(k-1)*neo*iabs(nag)+l*iabs(nag)+m)
              zz3 = frgdd(kk2+(k-1)*neo*iabs(nag)+(l-1)*iabs(nag)+m+1)
              zz4 = frgdd(kk2+(k-1)*neo*iabs(nag)+l*iabs(nag)+m+1)

              aaa = (zz2-zz1)/(xx2-xx1)
              bbb = (xx2*zz1-xx1*zz2)/(xx2-xx1)
              cc1 = aaa*(xx2**2-xx1**2)/2 + bbb*(xx2-xx1)

              aaa = (zz4-zz3)/(xx2-xx1)
              bbb = (xx2*zz3-xx1*zz4)/(xx2-xx1)
              cc2 = aaa*(xx2**2-xx1**2)/2 + bbb*(xx2-xx1)

              aaa = (cc2-cc1)/(yy2-yy1)
              bbb = (yy2*cc1-yy1*cc2)/(yy2-yy1)
              ccc = aaa*(yy2**2-yy1**2)/2 + bbb*(yy2-yy1)
              ccc = -ccc

              sek = sek + ccc *2.d0*pi

            end if

               frgdd(kk2+(k-1)*neo*iabs(nag)+(l-1)*iabs(nag)+m)
     &              = sek

            end do
           end do

           if( sek .gt. 0.0d0 ) then
            do l = 1, nef-1
             do m = 1, naf-1

                frgdd(kk2+(k-1)*neo*iabs(nag)+(l-1)*iabs(nag)+m)
     &               = frgdd(kk2+(k-1)*neo*iabs(nag)
     &               +(l-1)*iabs(nag)+m) / sek

             end do
            end do
           end if

           if ( dabs(frgsf(ksf+k)) .lt. 1d-9 ) frgsf(ksf+k) = sek

           if ( frgxs(kxs+j) .gt. 0d0 ) then
              frgsf(ksf+k) = frgsf(ksf+k) / frgxs(kxs+j)
              if ( frgsf(ksf+k) .gt. 1d3 ) frgsf(ksf+k) = 1d0
           else
              frgsf(ksf+k) = 0d0 ! when frgxs(kxs+j) = 0
           end if

          end do


C case 2
         else if ( neo .eq. 0 .and. nag .ne. 0 ) then

          do k = 1, nfrg

           sek = 0.0d0
           do m = 1, naf-1

            if ( ifgdf(i,1) .ne. 5 ) then

              sek = sek + frgdx(kk3+(k-1)*iabs(nag)+m)
     &             * 2.d0 * pi
     &             * ( dcos( frgaf(kaf+m) )
     &             - dcos( frgaf(kaf+m+1) ) )

            else

              yy1 = dcos( frgaf(kaf+m) )
              yy2 = dcos( frgaf(kaf+m+1) )
              zz1 = frgdx(kk3+(k-1)*iabs(nag)+m)
              zz2 = frgdx(kk3+(k-1)*iabs(nag)+m+1)
              aaa = (zz2-zz1)/(yy2-yy1)
              bbb = (yy2*zz1-yy1*zz2)/(yy2-yy1)
              ccc = aaa*(yy2**2-yy1**2)/2 + bbb*(yy2-yy1)
              ccc = -ccc

              sek = sek + ccc

            end if

              frgdx(kk3+(k-1)*iabs(nag)+m) = sek

           end do

           if( sek .gt. 0.0d0 ) then
            do m = 1, iabs(nag)

               frgdx(kk3+(k-1)*iabs(nag)+m) =
     &              frgdx(kk3+(k-1)*iabs(nag)+m) / sek

            end do
           end if

           if ( dabs(frgsf(ksf+k)) .lt. 1d-9 ) frgsf(ksf+k) = sek

           if ( frgxs(kxs+j) .gt. 0d0 ) then
              frgsf(ksf+k) = frgsf(ksf+k) / frgxs(kxs+j)
              if ( frgsf(ksf+k) .gt. 1d3 ) frgsf(ksf+k) = 1d0
           else
              frgsf(ksf+k) = 0d0 ! when frgxs(kxs+j) = 0
           end if

          end do


C case 3
         else if ( neo .eq. 0 .and. nag .eq. 0 ) then

          do k = 1, nfrg

           if ( dabs(frgsf(ksf+k)) .lt. 1d-9 ) goto 340

           if ( frgxs(kxs+j) .gt. 0d0 ) then
              frgsf(ksf+k) = frgsf(ksf+k) / frgxs(kxs+j)
              if ( frgsf(ksf+k) .gt. 1d3 ) frgsf(ksf+k) = 1d0
           else
              frgsf(ksf+k) = 0d0 ! when frgxs(kxs+j) = 0
           end if

          end do


C case 4
         else if ( neo .gt. 0 .and. nag .eq. 0 ) then

          iag = 1
          do k = 1, nfrg

           sek = 0.0d0
           do l = 1, nef-1

            if ( ifgdf(i,1) .ne. 5 ) then

              sek = sek + frgdx(kk3+(k-1)*nef+l)
     &             * ( frgef(kef+l+1) - frgef(kef+l) )

            else

              xx1 = frgef(kef+l)
              xx2 = frgef(kef+l+1)
              zz1 = frgdx(kk3+(k-1)*nef+l)
              zz2 = frgdx(kk3+(k-1)*nef+l+1)
              aaa = (zz2-zz1)/(xx2-xx1)
              bbb = (xx2*zz1-xx1*zz2)/(xx2-xx1)
              ccc = aaa*(xx2**2-xx1**2)/2 + bbb*(xx2-xx1)

              sek = sek + ccc

            end if

              frgdx(kk3+(k-1)*nef+l) = sek

           end do

           if( sek .gt. 0.0d0 ) then
            do l = 1, nef-1

               frgdx(kk3+(k-1)*nef+l) =
     &              frgdx(kk3+(k-1)*nef+l) / sek

            end do
           end if

           if ( dabs(frgsf(ksf+k)) .lt. 1d-9 ) frgsf(ksf+k) = sek

           if ( frgxs(kxs+j) .gt. 0d0 ) then
              frgsf(ksf+k) = frgsf(ksf+k) / frgxs(kxs+j)
              if ( frgsf(ksf+k) .gt. 1d3 ) frgsf(ksf+k) = 1d0
           else
              frgsf(ksf+k) = 0d0 ! when frgxs(kxs+j) = 0
           end if

          end do


C case 5
         else if ( neo .lt. 0 .and. nag .ne. 0 ) then

          do k = 1, nfrg

           sek = 0.0d0
           do l = 1, iabs(neo)
            do m = 1, iabs(nag)

               sek = sek + frgdd(kk2+(k-1)*iabs(neo)*iabs(nag)
     &              +(l-1)*iabs(nag)+m)
     &              * 2.d0 * pi
     &              * ( dcos( frgaf(kaf+m) )
     &              - dcos( frgaf(kaf+m+1) ) )

               frgdd(kk2+(k-1)*iabs(neo)*iabs(nag)+(l-1)*iabs(nag)+m)
     &              = sek

            end do
           end do

           if( sek .gt. 0.0d0 ) then
            do l = 1, iabs(neo)
             do m = 1, iabs(nag)

                frgdd(kk2+(k-1)*iabs(neo)*iabs(nag)+(l-1)*iabs(nag)+m)
     &               = frgdd(kk2+(k-1)*iabs(neo)*iabs(nag)
     &               +(l-1)*iabs(nag)+m) / sek

             end do
            end do
           end if

           if ( dabs(frgsf(ksf+k)) .lt. 1d-9 ) frgsf(ksf+k) = sek

           if ( frgxs(kxs+j) .gt. 0d0 ) then
              frgsf(ksf+k) = frgsf(ksf+k) / frgxs(kxs+j)
              if ( frgsf(ksf+k) .gt. 1d3 ) frgsf(ksf+k) = 1d0
           else
              frgsf(ksf+k) = 0d0 ! when frgxs(kxs+j) = 0
           end if

          end do


C case 6
         else if ( neo .lt. 0 .and. nag .eq. 0 ) then

          do k = 1, nfrg

           sek = 0.0d0
           do l = 1, iabs(neo)

              sek = sek + frgdx(kk3+(k-1)*iabs(neo)+l)

              frgdx(kk3+(k-1)*iabs(neo)+l) = sek

           end do

           if( sek .gt. 0.0d0 ) then
            do l = 1, iabs(neo)

               frgdx(kk3+(k-1)*iabs(neo)+l) =
     &              frgdx(kk3+(k-1)*iabs(neo)+l) / sek

            end do
           end if

           if ( dabs(frgsf(ksf+k)) .lt. 1d-9 ) frgsf(ksf+k) = sek

           if ( frgxs(kxs+j) .gt. 0d0 ) then
              frgsf(ksf+k) = frgsf(ksf+k) / frgxs(kxs+j)
              if ( frgsf(ksf+k) .gt. 1d3 ) frgsf(ksf+k) = 1d0
           else
              frgsf(ksf+k) = 0d0 ! when frgxs(kxs+j) = 0
           end if

          end do

         end if

        end do ! end of do-loop for incident energy mesh points

*-----------------------------------------------------------------------

        goto 399

*-----------------------------------------------------------------------

 340    continue
 1340   format('** Error : when neo=0 and nag=0, production XS '
     &       /'must not be 0 in ',a)
        write(io,1340) frgfl(i)(1:ifgdf(i,4))
        ErrCha = ''
        MsgID = 'L:517/R:FragData_cumulative/F:fragdatamod.f'
        call ErrWrite(MsgID, ErrCha)
        write(jo,1340) frgfl(i)(1:ifgdf(i,4))
        ierr = ierr + 1

*-----------------------------------------------------------------------

 399    continue
       end if
      end do ! end of do-loop for frag data files

*-----------------------------------------------------------------------

      return
      end subroutine FragData_cumulative
*-----------------------------------------------------------------------



************************************************************************
*                                                                      *
      subroutine FragData_engang(icase,ieo,iag,ifrg,iei,i
     &                  ,eout,csth,wxs)
*                                                                      *
*       calculation of cross section weight at eout and csth           *
*       modified by S.Hashimoto on 2019/12/10                          *
*                                                                      *
************************************************************************
      use moddas_fragdata

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'ggsparam.inc'

*-----------------------------------------------------------------------

      common /fgdata/ ifrgd, ifgm, ifgdf(kvlmax,5), frgfl(kvlmax)
      character frgfl*200


*-----------------------------------------------------------------------

      if ( icase .eq. 1 ) then

         kef = ifgs05_kef(i)
         kaf = ifgs07_kaf(i)
         xx1 = frgef(kef+ieo)
         xx2 = frgef(kef+ieo+1)
         yy1 = dcos( frgaf(kaf+iag) )
         yy2 = dcos( frgaf(kaf+iag+1) )
         zz1 = abcd(ieo,iag,ifrg,iei,i)
         zz2 = abcd(ieo+1,iag,ifrg,iei,i)
         zz3 = abcd(ieo,iag+1,ifrg,iei,i)
         zz4 = abcd(ieo+1,iag+1,ifrg,iei,i)
C S.H. revised the above calculations to avoid division by zero when e.g. eout=xx1 (2021.3.4)
         zbar = (zz1+zz2+zz3+zz4)/4d0
         if ( zbar .ne. 0d0 ) then
            wxs = ( (xx2-eout)*zz1 + (eout-xx1)*zz2 ) * (yy2-csth)
     1           + ( (xx2-eout)*zz3 + (eout-xx1)*zz4 ) * (csth-yy1)
            wxs = wxs /(xx2-xx1) /(yy2-yy1) /zbar
         else
            wxs = 0d0
         end if

*-----------------------------------------------------------------------

      else if ( icase .eq. 2 ) then

         kaf = ifgs07_kaf(i)
         yy1 = dcos( frgaf(kaf+iag) )
         yy2 = dcos( frgaf(kaf+iag+1) )
         zz1 = abcd(ieo,iag,ifrg,iei,i)
         zz2 = abcd(ieo,iag+1,ifrg,iei,i)
         zbar = (zz1+zz2)/2d0
         if ( zbar .ne. 0d0 ) then
            wxs = (yy2-csth)*zz1 + (csth-yy1)*zz2
            wxs = wxs /(yy2-yy1) /zbar
         else
            wxs = 0d0
         end if

*-----------------------------------------------------------------------

      else if ( icase .eq. 4 ) then

         kef = ifgs05_kef(i)
         xx1 = frgef(kef+ieo)
         xx2 = frgef(kef+ieo+1)
         zz1 = abcd(ieo,iag,ifrg,iei,i)
         zz2 = abcd(ieo+1,iag,ifrg,iei,i)
         zbar = (zz1+zz2)/2d0
         if ( zbar .ne. 0d0 ) then
            wxs = (xx2-eout)*zz1 + (eout-xx1)*zz2
            wxs = wxs /(xx2-xx1) /zbar
         else
            wxs = 0d0
         end if

*-----------------------------------------------------------------------

      end if

      return
      end subroutine FragData_engang
*-----------------------------------------------------------------------


************************************************************************
*                                                                      *
      subroutine deallocate_FragData
*                                                                      *
*       deallocate memory of frag data cross sections                  *
*       modified by S.Hashimoto on 2019/12/10                          *
*                                                                      *
************************************************************************

      deallocate( abcd )

      end subroutine deallocate_FragData
*-----------------------------------------------------------------------


*                                                                      *
*                                                                      *
      end module fragdatamod
************************************************************************
