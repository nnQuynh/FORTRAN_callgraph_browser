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
     &             " ny",i4.4," dy",i4.4,"=[y",i4.4,"*y",i4.4,"] "))')
     &             ( trim(chpchl_n(i,1))
     &              ,i+np_mxang,i,i,i+np_mxang, i = 1, np_mxang )
       end if
      end if

      write(iot,'( "#  r-lower      r-upper",3x,
     &             1000(2x,a11,"  r.err "))')
     &             ( chp_n(i), i = 1, nsame )

      tott_n(:,:) = 0.0d+0

      samewtt(:) = 0.0d0

      do iri = 1, nr, nrstepi
          write(iot,'(1p2e13.4,1000(1pe13.4,0pf8.4))')
     &      rm(iri),rm(iri+nrstepi),
     &      (((((tr(ipi+iploop-1,iei+ieloop-1,
     &               iri+irloop-1,izi+izloop-1,k),k=1,2),
     &      iploop=1,npstepi),ieloop=1,nestepi),
     &      irloop=1,nrstepi),izloop=1,nzstepi)


            i = 0
              do izloop=1,nzstepi
                do irloop=1,nrstepi
                  do ieloop=1,nestepi
                    do iploop=1,npstepi
                      i = i + 1
                      vn = tr(ipi+iploop-1,iei+ieloop-1,
     &                        iri+irloop-1,izi+izloop-1,1)
     &                   * vl(iri+irloop-1,izi+izloop-1)

                      tott_n(i,1) = tott_n(i,1) + vn

            samewtt(i) = samewtt(i) + vl(iri+irloop-1,izi+izloop-1)

                    end do
                  end do
                end do
              end do

      end do

      call pletrz_sumover_getput(0,m,iax,
     &    npstepi,nestepi,nrstepi,nzstepi,
     &    ipi,iei,1,izi,maxtott_n,tott_n)  ! iri = 1

      write(iot,'(/"#   sum over ",13x,
     &              1000(1pe13.4,0pf8.4))')
     &              (tott_n(i,1),tott_n(i,2),i=1,nsame)

