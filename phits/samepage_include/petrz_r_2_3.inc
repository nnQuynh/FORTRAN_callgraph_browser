c
c nais samepage
      if ( iteps(m) .ne. 2 ) then

        write(iot,'( "h: n",12x,"x",12x,
     &             1000(a," n "))')
     &            ( trim(chpchl_n(i,1)) , i = 1, np_mxang )

      else

        write(iot,'( "h: n",12x,"x",12x,
     &         1000(a,
     &         " ny",i0," dy",i0,"=[y",i0,"*y",i0,"] "))')
     &         ( trim(chpchl_n(i,1))
     &         ,i+np_mxang,i,i,i+np_mxang, i = 1, np_mxang )

      end if

      write(iot,'( "#  r-lower      r-upper",3x,
     &             1000(2x,a11,"  r.err "))')
     &             ( chp_n(i), i = 1, nsame )

      tott_n(:,:) = 0.0d+0

      samewtt(:) = 0.0

      do iri = 1, nr, nrstepi
          write(iot,'(1p2e13.4,1000(1pe13.4,0pf8.4))')
     &      rm(iri),rm(iri+nrstepi),
     &      ((((((tr(ipi+iploop-1,iei+ieloop-1,iti+itloop-1,
     &                iri+irloop-1,izi+izloop-1,
     &      k),k=1,2),
     &      iploop=1,npstepi),ieloop=1,nestepi),itloop=1,ntstepi),
     &      irloop=1,nrstepi),izloop=1,nzstepi)

            i = 0
              do izloop=1,nzstepi
                do irloop=1,nrstepi
                    do itloop=1,ntstepi
                        do ieloop=1,nestepi
                          do iploop=1,npstepi
                            i = i + 1
                            if( itunt(m) .le. 2 ) then
                              vm = 1.d0
                              samewtt(i) = 1.d0
                            else if( itunt(m) .ge. 3 ) then
                              vm =
     &                        vl(iri+irloop-1,izi+izloop-1)
                              samewtt(i) = samewtt(i) + vm
                            end if

                            vn = tr(ipi+iploop-1,iei+ieloop-1,
     &                              iti+itloop-1,iri+irloop-1,
     &                              izi+izloop-1,1) * vm

                            tott_n(i,1) = tott_n(i,1) + vn

                          end do

                        end do
                    end do
                end do
              end do

      end do

      do i = 1, nsame
        if( itunt(m) .ge. 3 ) then
          tott_n(i,1) = tott_n(i,1) / samewtt(i)
        endif
      end do

      call ptime_star_rz_sumover_getput(0,m,iax,
     &    npstepi,nestepi,ntstepi,nrstepi,nzstepi,
     &    ipi,iei,iti,1,izi,maxtott_n,tott_n)   ! iri = 1

      if( itunt(m) .le. 2 ) then
        write(iot,'(/"#   sum over ",13x,
     &              1000(1pe13.4,0pf8.4))')
     &             (tott_n(i,1),tott_n(i,2),i=1,nsame)
      else
        if(iloopmode .eq. 0) then
          write(iot,'(/"#   sum over ",1pe13.4,
     &                1000(1pe13.4,0pf8.4))')
     &                samewtt(1),(tott_n(i,1),tott_n(i,2),i=1,nsame)
        else
          write(iot,'(/"#   sum over ",13x,
     &                1000(1pe13.4,0pf8.4))')
     &                (tott_n(i,1),tott_n(i,2),i=1,nsame)

          write(iot,'(/"#   sum over ",5x,"wtt",5x,
     &                1000(1pe13.4,8x))')
     &                (samewtt(i),i=1,nsame)

        end if
      end if
