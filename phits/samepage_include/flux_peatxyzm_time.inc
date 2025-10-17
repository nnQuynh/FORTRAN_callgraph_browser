c
c nais samepage
      write(iot,'( "h: n",12x,"x",12x,
     &             1000(a," n  "))')
     &            ( trim(chpchl_n(i,1)) , i = 1, nsame )

      write(iot,'( "#  t-lower      t-upper",3x,
     &             1000("f",1x,a11,"  r.err "))')
     &             ( chp_n(i), i = 1, nsame )

      tott_n(:,:) = 0.0d+0

      samewtt(:) = 0.0d0

      do iti = 1, nt, ntstepi
          write(iot,'(1p2e13.4,1000(1pe13.4,0pf8.4))')
     &      tb(iti),tb(iti+ntstepi),
     &      (((((((((tr(ipi+iploop-1,iei+ieloop-1,iai+ialoop-1,
     &      iti+itloop-1,
     &       icf(ixi+ixloop-1,iyi+iyloop-1,izi+izloop-1),
     &      imi+imloop-1,k),k=1,2),
     &      iploop=1,npstepi),ieloop=1,nestepi),ialoop=1,nastepi),
     &      itloop=1,ntstepi),
     &      ixloop=1,nxstepi),iyloop=1,nystepi),izloop=1,nzstepi),
     &      imloop=1,nmstepi)


            i = 0
            do imloop=1,nmstepi
              do izloop=1,nzstepi
                do iyloop=1,nystepi
                  do ixloop=1,nxstepi
                    do itloop=1,ntstepi
                      do ialoop=1,nastepi
                        do ieloop=1,nestepi
                          do iploop=1,npstepi
                            i = i + 1
                            vn = tr(ipi+iploop-1,iei+ieloop-1,
     &                              iai+ialoop-1,iti+itloop-1,
     &                      icf(ixi+ixloop-1,iyi+iyloop-1,izi+izloop-1),
     &                              imi+imloop-1,1)
     &                         * tw(iti+itloop-1)

                            tott_n(i,1) = tott_n(i,1) + vn

                            samewtt(i) = samewtt(i) + tw(iti+itloop-1)

                          end do
                        end do
                      end do
                    end do
                  end do
                end do
              end do
            end do

        end do

        do i = 1, nsame
            if( itunt(m) .gt.  10 ) then
               tott_n(i,1) = tott_n(i,1) / samewtt(i)
            end if

        end do

        write(iot,'(/"#   sum over",14x,1000(1pe13.4,8x))')
     &       (tott_n(i,1),i=1,nsame)

