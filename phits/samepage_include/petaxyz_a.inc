c
c nais samepage
      if ( iteps(m) .ne. 2 ) then

        write(iot,'( "h: n",12x,"x",12x,
     &             1000(a," n "))')
     &            ( trim(chpchl_n(i,1)) , i = 1, np_mxang )

      else
       if(iloopmode .eq. 0) then
        write(iot,'( "h: n",12x,"x",12x,
     &         1000(a,
     &         " ny",i0," dy",i0,"=[y",i0,"*y",i0,"] "))')
     &         ( trim(chpchl_n(i,1))
     &         ,i+np_mxang,i,i,i+np_mxang, i = 1, np_mxang )
       else
        write(iot,'( "h: n",12x,"x",12x,
     &         1000(a,
     &         " ny",i4.4," dy",i4.4,"=[y",i4.4,"*y",i4.4,"] "))')
     &         ( trim(chpchl_n(i,1))
     &         ,i+np_mxang,i,i,i+np_mxang, i = 1, np_mxang )
       end if
      end if

      write(iot,'( "#  a-lower      a-upper",3x,
     &             1000(2x,a11,"  r.err "))')
     &             ( chp_n(i), i = 1, nsame )

      tott_n(:,:) = 0.0d+0

      samewtt(:) = 0.0d0

      if( ( itaxs(m,iax) .eq. 12 .and.
     &      itaty(m) .gt. 0 ) .or.
     &    ( itaxs(m,iax) .eq. 13 .and.
     &      itaty(m) .lt. 0 ) ) then

       ias = 1
       iaf = na
       iad = 1

      else if( ( itaxs(m,iax) .eq. 12 .and.
     &          itaty(m) .lt. 0 ) .or.
     &        ( itaxs(m,iax) .eq. 13 .and.
     &          itaty(m) .gt. 0 ) ) then

        ias = na
        iaf = 1
        iad = -1

      end if

      nastepi = 1
      do iai = ias, iaf, iad

        if( ( itaxs(m,iax) .eq. 12 .and.
     &        itaty(m) .gt. 0 ) .or.
     &      ( itaxs(m,iax) .eq. 13 .and.
     &        itaty(m) .lt. 0 ) ) then

          abf = ab(iai)
          abt = ab(iai+1)

        else if( itaxs(m,iax) .eq. 12 .and.
     &           itaty(m) .lt. 0 ) then

          abf = cos( ab(iai+1) / 180.d0 * pi )
          abt = cos( ab(iai) / 180.d0 * pi )

        else if( itaxs(m,iax) .eq. 13 .and.
     &           itaty(m) .gt. 0 ) then

          abf = acos( ab(iai+1) ) * 180.d0 / pi
          abt = acos( ab(iai) ) * 180.d0 / pi

        end if

        write(iot,'(1p2e13.4,1000(1pe13.4,0pf8.4))')
     &      abf,abt,
     &      ((((((((tr(ipi+iploop-1,iei+ieloop-1,
     &             iti+itloop-1,iai+ialoop-1,
     &       icf(ixi+ixloop-1,iyi+iyloop-1,izi+izloop-1),
     &      k),k=1,2),
     &      iploop=1,npstepi),ieloop=1,nestepi),
     &      itloop=1,ntstepi),ialoop=1,nastepi),
     &      ixloop=1,nxstepi),iyloop=1,nystepi),izloop=1,nzstepi)


            i = 0
            do izloop=1,nzstepi
              do iyloop=1,nystepi
                do ixloop=1,nxstepi
                  do ialoop=1,nastepi
                    do itloop=1,ntstepi
                        do ieloop=1,nestepi
                          do iploop=1,npstepi
                            i = i + 1
                            vn = tr(ipi+iploop-1,iei+ieloop-1,
     &                              iti+itloop-1,iai+ialoop-1,
     &                      icf(ixi+ixloop-1,iyi+iyloop-1,izi+izloop-1),
     &                              1)
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

      end do

      do i = 1, nsame
        if( itunt(m) .gt.  20 ) then
          tott_n(i,1) = tott_n(i,1) / samewtt(i)
        end if
      end do

      call  ppdctxyz_sumover_getput(0,m,iax,
     &    npstepi,nestepi,ntstepi,nastepi,nxstepi,nystepi,nzstepi,
     &    ipi,iei,iti,1,ixi,iyi,izi,maxtott_n,tott_n)  ! iai

      write(iot,'(/"#   sum over ",13x,
     &              1000(1pe13.4,0pf8.4))')
     &              (tott_n(i,1),tott_n(i,2),i=1,nsame)


      if( itunt(m) .eq.  2 .or. itunt(m) .eq.  4 .or.
     &    itunt(m) .eq.  6 .or. itunt(m) .eq. 12 .or.
     &    itunt(m) .eq. 14 .or. itunt(m) .eq. 16 .or.
     &    itunt(m) .eq. 22 .or. itunt(m) .eq. 24 .or.
     &    itunt(m) .eq. 26 .or. itunt(m) .eq. 32 .or.
     &    itunt(m) .eq. 34 .or. itunt(m) .eq. 36 ) then

         vm = vl(ixi,iyi,izi)
      else
        vm = 1.d+0
      end if

