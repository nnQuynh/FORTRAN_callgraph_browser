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

      write(iot,'( "#  r-lower      r-upper",3x,
     &             1000(2x,a11,"  r.err "))')
     &             ( chp_n(i), i = 1, nsame )

      tott_n(:,:) = 0.0d+0

      samewtt(:) = 0.0d0

      nrstepi = 1
      do iri = 1, nr, nrstepi
          write(iot,'(1p2e13.4,1000(1pe13.4,0pf8.4))')
     &      rm(iri),rm(iri+nrstepi),
     &      ((((((tr(ipi+iploop-1,iei+ieloop-1,iti+itloop-1,
     &                iri+irloop-1,izi+izloop-1,k),k=1,2),
     &      iploop=1,npstepi),ieloop=1,nestepi),itloop=1,ntstepi),
     &      irloop=1,nrstepi),izloop=1,nzstepi)


            i = 0
              do izloop=1,nzstepi
                do irloop=1,nrstepi
                    do itloop=1,ntstepi
                        do ieloop=1,nestepi
                          do iploop=1,npstepi
                            i = i + 1
                            vn = tr(ipi+iploop-1,iei+ieloop-1,
     &                              iti+itloop-1,iri+irloop-1,
     &                              izi+izloop-1,1)
     &                         * vl(iri+irloop-1,izi+izloop-1)

                            tott_n(i,1) = tott_n(i,1) + vn

            samewtt(i) = samewtt(i) + vl(iri+irloop-1,izi+izloop-1)

                          end do
                        end do
                    end do
                end do
              end do

      end do

      do i = 1, nsame
        if( itunt(m) .le. 1 .or. itunt(m) .eq. 5 ) then
          tott_n(i,1) = tott_n(i,1) / samewtt(i)
        end if
      end do

      call  pdepstrz_sumover_getput(0,m,iax,
     &    npstepi,nestepi,ntstepi,nrstepi,nzstepi,
     &    ipi,iei,iti,1,izi,maxtott_n,tott_n)  ! iri

      if( itunt(m) .le. 1 .or. itunt(m) .eq. 5 ) then
        if(iloopmode .eq. 0) then
          write(iot,'(/"#   sum over ",1pe13.4,
     &        1000(1pe13.4,0pf8.4))')
     &        samewtt(1), (tott_n(i,1),tott_n(i,2),i=1,nsame) ! frtati
        else
          write(iot,'(/"#   sum over ",13x,
     &        1000(1pe13.4,0pf8.4))')
     &        (tott_n(i,1),tott_n(i,2),i=1,nsame) ! frtati
          write(iot,'(/"#   sum over ",5x,"wtt",5x,
     &        1000(1pe13.4,8x))')
     &        (samewtt(i),i=1,nsame) ! frtati
        end if
      else
        write(iot,'(/"#   sum over ",13x,
     &        1000(1pe13.4,0pf8.4))')
     &        (tott_n(i,1),tott_n(i,2),i=1,nsame) ! frtati 2021/10/05
      end if

