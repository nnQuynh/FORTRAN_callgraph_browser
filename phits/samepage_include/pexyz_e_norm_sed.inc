c
c nais samepage
      if ( iteps(m) .ne. 2 ) then
        write(iot,'( "h: n",12x,"x",12x,
     &               1000(a," n  "))')
     &              ( trim(chpchl_n(i,1)) , i = 1, np_mxang )
      else
       if(iloopmode .eq. 0) then
        write(iot,'( "h: n",12x,"x",12x,
     &               1000(a," ny",i0," dy",i0,"=[y",i0,"*y",i0,"] "))')
     &             ( trim(chpchl_n(i,1))
     &              ,i+np_mxang,i,i,i+np_mxang, i = 1, np_mxang )
       else
        write(iot,'( "h: n",12x,"x",12x,
     &               1000(a,
     &              " ny",i4.4," dy",i4.4,"=[y",i4.4,"*y",i4.4,"] "))')
     &             ( trim(chpchl_n(i,1))
     &              ,i+np_mxang,i,i,i+np_mxang, i = 1, np_mxang )
       end if
      end if

               write(iot,'( "#  sed-lwr      sed-upr",3x,
     &             1000(2x,a11,"  r.err "))')
     &             ( chp_n(i), i = 1, nsame )

      tott_n(:,:) = 0.0d+0

      samewtt(:) = 0.0d0

      if(itunt(m).ge.7) then

        totnorm_n(:)=0.0

        i = 0
        do izloop=1,nzstepi
        do iyloop=1,nystepi
        do ixloop=1,nxstepi
        do iploop=1,npstepi
          i = i + 1
          do ie=1,ne
               totnorm_n(i)=totnorm_n(i)
     &         + tr(ipi+iploop-1,ie,
     &         icf(ixi+ixloop-1,iyi+iyloop-1,izi+izloop-1),1)
     &         *(eb(ie+1)-eb(ie))/((eb(ie+1)+eb(ie))*0.5d0)

               tott_n(i,1)=tott_n(i,1)+tr(ipi+iploop-1,ie,
     &         icf(ixi+ixloop-1,iyi+iyloop-1,izi+izloop-1),1)
     &         *(eb(ie+1)-eb(ie))

               tott_n(i,2)=tott_n(i,2)
     &                    + (tr(ipi+iploop-1,ie,
     &         icf(ixi+ixloop-1,iyi+iyloop-1,izi+izloop-1),1)
     &         * tr(ipi+iploop-1,ie,
     &         icf(ixi+ixloop-1,iyi+iyloop-1,izi+izloop-1),2))
     &         **2
     &         *(eb(ie+1)-eb(ie))**2

          enddo
          if(totnorm_n(i).eq.0.0) totnorm_n(i)=1.0d0

          if(tott_n(i,1).gt. 0.0d0) then
            tott_n(i,2) = sqrt(tott_n(i,2))/tott_n(i,1)
          else
            tott_n(i,2) = 0.0d0
          endif

          tott_n(i,1) = tott_n(i,1)/totnorm_n(i)

        enddo
        enddo
        enddo
        enddo

      else ! conventional mode
       totnorm_n(:)=1.0d0
      endif

      do iei = 1, ne, nestepi

        i = 0
        do izloop=1,nzstepi
        do iyloop=1,nystepi
        do ixloop=1,nxstepi
          do ieloop=1,nestepi
            do iploop=1,npstepi
              i = i + 1
              if(itunt(m).ge.7) then
                tr_n(i,1) = tr(ipi+iploop-1,iei+ieloop-1,
     &              icf(ixi+ixloop-1,iyi+iyloop-1,izi+izloop-1),1)
     &                    / totnorm_n(i)
!     &                    * (eb(iei+1)-eb(iei))
              else
                tr_n(i,1) = tr(ipi+iploop-1,iei+ieloop-1,
     &              icf(ixi+ixloop-1,iyi+iyloop-1,izi+izloop-1),1)
     &                  /totnorm_n(i)
              endif
              tr_n(i,2) = tr(ipi+iploop-1,iei+ieloop-1,
     &              icf(ixi+ixloop-1,iyi+iyloop-1,izi+izloop-1),2)

            end do
          end do
        end do
        end do
        end do


          write(iot,'(1p2e13.4,1000(1pe13.4,0pf8.4))')
     &      eb(iei),eb(iei+nestepi),
     &       (tr_n(i,1),tr_n(i,2),i=1,nsame)


          if( itunt(m) .le. 6 ) then
            i = 0
              do izloop=1,nzstepi
                do iyloop=1,nystepi
                  do ixloop=1,nxstepi
                        do ieloop=1,nestepi
                          do iploop=1,npstepi
                            i = i + 1
                            vn = tr(ipi+iploop-1,iei+ieloop-1,
     &                      icf(ixi+ixloop-1,iyi+iyloop-1,izi+izloop-1),
     &                      1)
     &                      * ew(iei+ieloop-1)

                            tott_n(i,1) = tott_n(i,1) + vn

                            samewtt(i) = samewtt(i) + ew(iei+ieloop-1)
                          end do
                        end do
                  end do
                end do
              end do
             end if

      end do

      if( itunt(m) .le. 4 ) then
        do i = 1, nsame
          tott_n(i,1) = tott_n(i,1) / samewtt(i)
        end do
      end if

      if(itunt(m).le.6) then
        call pletxyz_sumover_getput(0,m,iax,
     &    npstepi,nestepi,nxstepi,nystepi,nzstepi,
     &    ipi,1,ixi,iyi,izi,maxtott_n,tott_n)  ! iei = 1
      endif

      zero=0.0  ! T.Sato 2025/03/13 no uncertainty unit = 13 & 14
        if(itunt(m).eq.7) then
               write(iot,'(/"#          Frequency mean ",
     &         1000(1pe13.4,0pf8.4))')
     &         (tott_n(i,1),zero,i=1,nsame)
        elseif(itunt(m).eq.8) then
               write(iot,'(/"#               Dose mean ",
     &         1000(1pe13.4,0pf8.4))')
     &         (tott_n(i,1),zero,i=1,nsame)
        else
               write(iot,'(/"#   sum over in cm or MeV ",
     &         1000(1pe13.4,0pf8.4))')
     &         (tott_n(i,1),tott_n(i,2),i=1,nsame)
        endif

! T.Sato 2025/03/12 for output normalization factor for restart calculation
               write(iot,'("#    Normalization factor ",
     &         1000(1pe13.4,8x))') (totnorm_n(i),i=1,nsame)
