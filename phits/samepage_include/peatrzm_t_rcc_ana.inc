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

      write(iot,'( "#  t-lower      t-upper",3x,
     &             1000(2x,a11,"  r.err "))')
     &             ( chp_n(i), i = 1, nsame )

      tott_n(:,:) = 0.0d+0

      samewtt(:) = 0.0d0

      ntstepi = 1
      do iti = 1, nt, ntstepi
          write(iot,'(1p2e13.4,1000(1pe13.4,0pf8.4))')
     &      tb(iti),tb(iti+ntstepi),
     &      (((((((tr(ipi+iploop-1,iei+ieloop-1,
     &             iai+ialoop-1,iti+itloop-1,
     &             irf(iri+irloop-1,izi+izloop-1),
     &             imi,k),k=1,2),
     &      iploop=1,npstepi),ieloop=1,nestepi),
     &      ialoop=1,nastepi),itloop=1,ntstepi),
     &      irloop=1,nrstepi),izloop=1,nzstepi)


            i = 0
              do izloop=1,nzstepi
                do irloop=1,nrstepi
                  do itloop=1,ntstepi
                    do ialoop=1,nastepi
                        do ieloop=1,nestepi
                          do iploop=1,npstepi
                            i = i + 1
                            vn = tr(ipi+iploop-1,iei+ieloop-1,
     &                              iai+ialoop-1,iti+itloop-1,
     &                              irf(iri+irloop-1,izi+izloop-1),
     &                              imi,1)
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

      do i = 1, nsame
        if( itunt(m) .gt.  10 ) then
          tott_n(i,1) = tott_n(i,1) / samewtt(i)
        endif
      end do

      call psufrz_sumover_getput(0,m,iax,
     &    npstepi,nestepi,nastepi,ntstepi,nrstepi,nzstepi,nmstepi,
     &    ipi,iei,iai,1,iri,izi,imi,maxtott_n,tott_n)   !iti = 1

        write(iot,'(/"#   sum over ",13x,
     &                1000(1pe13.4,0pf8.4))')
     &                (tott_n(i,1),tott_n(i,2),i=1,nsame)

