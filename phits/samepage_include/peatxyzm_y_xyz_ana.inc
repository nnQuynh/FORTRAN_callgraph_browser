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
     &         " ny2",i1," dy",i1,"=[y",i1,"*y2",i1,"] "))')
     &         ( trim(chpchl_n(i,1))
     &         ,i,i,i,i, i = 1, np_mxang )
       else
        write(iot,'( "h: n",12x,"x",12x,
     &         1000(a,
     &         " ny2",i4.4," dy",i4.4,"=[y",i4.4,"*y2",i4.4,"] "))')
     &         ( trim(chpchl_n(i,1))
     &         ,i,i,i,i, i = 1, np_mxang )
       end if

      end if

      write(iot,'( "#  y-lower      y-upper",3x,
     &             1000(2x,a11,"  r.err "))')
     &             ( chp_n(i), i = 1, nsame )

      tott_n(:,:) = 0.0d+0

      samewtt(:) = 0.0d0

      nystepi = 1
      do iyi = 1, ny, nystepi
          write(iot,'(1p2e13.4,1000(1pe13.4,0pf8.4))')
     &      ym(iyi),ym(iyi+nystepi),
     &      ((((((((tr(ipi+iploop-1,iei+ieloop-1,
     &             iai+ialoop-1,iti+itloop-1,
     &             icf(ixi+ixloop-1,iyi+iyloop-1,izi+izloop-1),
     &             imi,k),k=1,2),
     &      iploop=1,npstepi),ieloop=1,nestepi),
     &      ialoop=1,nastepi),itloop=1,ntstepi),
     &      ixloop=1,nxstepi),iyloop=1,nystepi),izloop=1,nzstepi)


            i = 0
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
     &                          icf(ixi+ixloop-1,iyi+iyloop-1,
     &                              izi+izloop-1),imi,1)
     &                   * ax(ixi+ixloop-1,iyi+iyloop-1,izi+izloop-1)

                            tott_n(i,1) = tott_n(i,1) + vn

                            samewtt(i) = samewtt(i)
     &                   + ax(ixi+ixloop-1,iyi+iyloop-1,izi+izloop-1)

                          end do
                        end do
                    end do
                  end do
                end do
              end do
            end do
      end do

      do i = 1, nsame
          tott_n(i,1) = tott_n(i,1) / samewtt(i)
      end do

      call psufxyz_sumover_getput(0,m,iax,
     &    npstepi,nestepi,nastepi,ntstepi,nxstepi,nystepi,nzstepi,
     &    nmstepi,
     &    ipi,iei,iai,iti,ixi,1,izi,imi,maxtott_n,tott_n)   ! iy = 1

      if(iloopmode .eq. 0) then
        write(iot,'(/"#   sum over ",1pe13.4,
     &                1000(1pe13.4,0pf8.4))')
     &               samewtt(1), (tott_n(i,1),tott_n(i,2),i=1,nsame)
      else
        write(iot,'(/"#   sum over ",13x,
     &                1000(1pe13.4,0pf8.4))')
     &               (tott_n(i,1),tott_n(i,2),i=1,nsame)
        write(iot,'(/"#   sum over ",5x,"wtt",5x,
     &                1000(1pe13.4,8x))')
     &               (samewtt(i),i=1,nsame)
      end if
