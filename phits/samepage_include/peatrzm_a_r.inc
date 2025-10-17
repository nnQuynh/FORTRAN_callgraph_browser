      if( itaxs(m,iax) .eq. 8 ) then
        if ( iteps(m) .ne. 2 ) then
          write(iot,'( "h: n",12x,"x",12x,
     &             1000(a," n "))')
     &           ( trim(chpchl_n(i,1)), i = 1, np_mxang )

        else

          write(iot,'( "h: n",12x,"x",12x,
     &         1000(a
     &         ," ny",i0," dy",i0,"=[y",i0,"*y",i0,"] "))')
     &         ( trim(chpchl_n(i,1))
     &         ,i+np_mxang,i,i,i+np_mxang, i = 1, np_mxang )
        end if

      else

        if ( iteps(m) .ne. 2 ) then
           write(iot,'( "h: x",12x,"n",12x,
     &           1000(a," n "))')
     &         ( trim(chpchl_n(i,2)) , i = 1, np_mxang )
        else
           write(iot,'( "h: x",12x,"n",12x,
     &         1000(a
     &         ," ny",i0," dy",i0,"=[y",i0,"*y",i0,"] "))')
     &         ( trim(chpchl_n(i,2))
     &         ,i+np_mxang,i,i,i+np_mxang, i = 1, np_mxang ) ! frtati 2021/10/05
        end if
      end if

      write(iot,'( "#  a-lower      a-upper",3x,
     &           1000(2x,a11,"  r.err "))')
     &          ( chp_n(i), i = 1, nsame )


      if( ( itaxs(m,iax) .eq. 8 .and.
     &      itaty(m) .gt. 0 ) .or.
     &    ( itaxs(m,iax) .eq. 10 .and.
     &      itaty(m) .lt. 0 ) ) then

         ias = 1
         iaf = na
         iad = 1

      else if( ( itaxs(m,iax) .eq. 8 .and.
     &           itaty(m) .lt. 0 ) .or.
     &         ( itaxs(m,iax) .eq. 10 .and.
     &           itaty(m) .gt. 0 ) ) then

         ias = na
         iaf = 1
         iad = -1

      end if

      tott_n(:,:) = 0.0d+0

      samewtt(:) = 0.0d0

      nastepi = 1

      do iai = ias, iaf, iad

        if( ( itaxs(m,iax) .eq. 8 .and.
     &        itaty(m) .gt. 0 ) .or.
     &      ( itaxs(m,iax) .eq. 10 .and.
     &        itaty(m) .lt. 0 ) ) then

          abf = ab(iai)
          abt = ab(iai+1)

        else if( itaxs(m,iax) .eq. 8 .and.
     &           itaty(m) .lt. 0 ) then

          abf = cos( ab(iai+1) / 180.d0 * pi )
          abt = cos( ab(iai) / 180.d0 * pi )

        else if( itaxs(m,iax) .eq. 10 .and.
     &           itaty(m) .gt. 0 ) then

          abf = acos( ab(iai+1) ) * 180.d0 / pi
          abt = acos( ab(iai) ) * 180.d0 / pi

        end if

        write(iot,'(1p2e13.4,1000(1pe13.4,0pf8.4))')
     &      abf, abt,
     &      (((((((tr(ipi+iploop-1,iei+ieloop-1,
     &             iai+ialoop-1,iti+itloop-1,
     &             irf(iri+irloop-1,izi+izloop-1),
     &             imi,k),k=1,2),
     &      iploop=1,npstepi),ieloop=1,nestepi),
     &      ialoop=1,nastepi),itloop=1,ntstepi),
     &      irloop=1,nrstepi),izloop=1,nzstepi)


            i = 0
              do izloop=1,nzstepi
                do irloop=1,nrstepi
                  do itloop=1,ntstepi
                    do ialoop=1,nastepi
                        do ieloop=1,nestepi
                          do iploop=1,npstepi
                            i = i + 1
                            vn = tr(ipi+iploop-1,iei+ieloop-1,
     &                              iai+ialoop-1,iti+itloop-1,
     &                              irf(iri+irloop-1,izi+izloop-1),
     &                              imi,1)
     &                         * aw(iai+ialoop-1)

                            tott_n(i,1) = tott_n(i,1) + vn

                            samewtt(i) = samewtt(i) + aw(iai+ialoop-1)

                          end do
                        end do
                    end do
                  end do
                end do
              end do

      end do

*-----------------------------------------------------------------------

      do i = 1, nsame
        if( itunt(m) .eq.  4 .or. itunt(m) .eq.  5 .or.
     &      itunt(m) .eq.  6 .or.
     &      itunt(m) .eq. 14 .or. itunt(m) .eq. 15 .or.
     &      itunt(m) .eq. 16 ) then
                     tott_n(i,1) = tott_n(i,1) / samewtt(i)
        end if
      end do

      call  psufrz_sumover_getput(0,m,iax,
     &    npstepi,nestepi,nastepi,ntstepi,nrstepi,nzstepi,nmstepi,
     &    ipi,iei,1,iti,iri,izi,imi,maxtott_n,tott_n)    ! iai = 1

      write(iot,'(/"#   sum over",14x,1000(1pe13.4,0pf8.4))')
     &           (tott_n(i,1),tott_n(i,2),i=1,nsame)

