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
     &      ((((((tr(ipi+iploop-1,iei+ieloop-1,iti+itloop-1,
     &      iri+irloop-1,imi+imloop-1,k),k=1,2),
     &      iploop=1,npstepi),ieloop=1,nestepi),itloop=1,ntstepi),
     &      irloop=1,nrstepi),imloop=1,nmstepi)


            i = 0
            do imloop=1,nmstepi
                    do irloop=1,nrstepi
                      do itloop=1,ntstepi
                        do ieloop=1,nestepi
                          do iploop=1,npstepi
                            i = i + 1
                            vn = tr(ipi+iploop-1,iei+ieloop-1,
     &                              iti+itloop-1,iri+irloop-1,
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

        do i = 1, nsame
            if( itunt(m) .gt.  10 ) then
               tott_n(i,1) = tott_n(i,1) / samewtt(i)
            end if

        end do

      call ptracreg_sumover_getput(0,m,iax,
     &    npstepi,nestepi,ntstepi,nrstepi,nmstepi, ! ntstepi = 1
     &    ipi,iei,1,iri,imi,maxtott_n,tott_n)      ! iti = 1

        write(iot,'(/"#   sum over",14x,1000(1pe13.4,0pf8.4))')
     &       (tott_n(i,1),tott_n(i,2),i=1,nsame)

