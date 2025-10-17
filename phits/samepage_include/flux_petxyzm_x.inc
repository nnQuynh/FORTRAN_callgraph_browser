c
c nais samepage
      write(iot,'( "h: n",12x,"x",12x,
     &             1000(a," n  "))')
     &            ( trim(chpchl_n(i,1)) , i = 1, nsame )

      write(iot,'( "#  x-lower      x-upper",3x,
     &             1000("f",1x,a11,"  r.err "))')
     &             ( chp_n(i), i = 1, nsame )

      tott_n(:,:) = 0.0d+0

      samewtt(:) = 0.0

      do ixi = 1, nx, nxstepi
          write(iot,'(1p2e13.4,1000(1pe13.4,0pf8.4))')
     &      xm(ixi),xm(ixi+nxstepi),
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
     &                      * vl(ixi+ixloop-1,iyi+iyloop-1,izi+izloop-1)

                            tott_n(i,1) = tott_n(i,1) + vn

                            samewtt(i) = samewtt(i)
     &                      + vl(ixi+ixloop-1,iyi+iyloop-1,izi+izloop-1)
                          end do
                        end do
                    end do
                  end do
                end do
              end do
            end do

      end do

      do i = 1, nsame
        if( itunt(m) .ne. 4 .and. itunt(m) .ne. 14 ) then
          tott_n(i,1) = tott_n(i,1) / samewtt(i)
        end if
      end do

      call ptracxyz_sumover_getput(0,m,iax,
     &    npstepi,nestepi,ntstepi,nxstepi,nystepi,nzstepi,nmstepi,      ! nxstepi = 1
     &    ipi,iei,iti,1,iyi,izi,imi,maxtott_n,tott_n)                   !  ixi = 1


      if( itunt(m) .ne. 4 .and. itunt(m) .ne. 14 ) then
        if(iloopmode .eq. 0) then
          write(iot,'(/"#   sum over ",1pe13.4,
     &              1000(1pe13.4,0pf8.4))')
     &              samewtt(1), (tott_n(i,1),tott_n(i,2),i=1,nsame)
        else
          write(iot,'(/"#   sum over ",13x,
     &              1000(1pe13.4,0pf8.4))')
     &              (tott_n(i,1),tott_n(i,2),i=1,nsame)
          write(iot,'(/"#   sum over ",5x,"wtt",5x,
     &              1000(1pe13.4,8x))')
     &              (samewtt(i),i=1,nsame)
        endif
      else
        write(iot,'(/"#   sum over ",13x,
     &              1000(1pe13.4,0pf8.4))')
     &              (tott_n(i,1),tott_n(i,2),i=1,nsame)
      end if

