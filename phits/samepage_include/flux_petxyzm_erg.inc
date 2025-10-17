c
c nais samepage
      write(iot,'( "h: n",12x,"x",12x,
     &             1000(a," n  "))')
     &            ( trim(chpchl_n(i,1)) , i = 1, nsame )

      write(iot,'( "#  e-lower      e-upper",3x,
     &             1000("f",1x,a11,"  r.err "))')
     &             ( chp_n(i), i = 1, nsame )

      tott_n(:,:) = 0.0d+0

      samewtt(:) = 0.0d0

      do iei = 1, ne, nestepi
          write(iot,'(1p2e13.4,1000(1pe13.4,0pf8.4))')
     &      eb(iei),eb(iei+nestepi),
     &      ((((((((tr(ipi+iploop-1,iei+ieloop-1,iti+itloop-1,
     &       icf(ixi+ixloop-1,iyi+iyloop-1,izi+izloop-1),
     &      imi+imloop-1,k),k=1,2),
     &      iploop=1,npstepi),ieloop=1,nestepi),itloop=1,ntstepi),
     &      ixloop=1,nxstepi),iyloop=1,nystepi),izloop=1,nzstepi),
     &      imloop=1,nmstepi)


            i = 0
            do imloop=1,nmstepi
              do izloop=1,nzstepi
                do iyloop=1,nystepi
                  do ixloop=1,nxstepi
                    do itloop=1,ntstepi
                        do ieloop=1,nestepi
                          do iploop=1,npstepi
                            i = i + 1
                            vn = tr(ipi+iploop-1,iei+ieloop-1,
     &                              iti+itloop-1,
     &                      icf(ixi+ixloop-1,iyi+iyloop-1,izi+izloop-1),
     &                              imi+imloop-1,1)
     &                         * ew(iei+ieloop-1)

                            tott_n(i,1) = tott_n(i,1) + vn

                            samewtt(i) = samewtt(i) + ew(iei+ieloop-1)

                          end do
                        end do
                    end do
                  end do
                end do
              end do
            end do

      end do

      do i = 1, nsame
        if( itunt(m) .eq.  2 .or. itunt(m) .eq.  3 .or.
     &      itunt(m) .eq. 12 .or. itunt(m) .eq. 13 ) then
          tott_n(i,1) = tott_n(i,1) / samewtt(i)
        end if
      end do

      call ptracxyz_sumover_getput(0,m,iax,
     &    npstepi,nestepi,ntstepi,nxstepi,nystepi,nzstepi,nmstepi,      ! nestepi = 1
     &    ipi,1,iti,ixi,iyi,izi,imi,maxtott_n,tott_n)                   !  iei = 1

      write(iot,'(/"#   sum over ",13x,1000(1pe13.4,0pf8.4))')
     &     (tott_n(i,1),tott_n(i,2),i=1,nsame)

