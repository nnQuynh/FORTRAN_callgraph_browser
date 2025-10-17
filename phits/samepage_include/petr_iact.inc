      if ( iteps(m) .ne. 2 ) then

        write(iot,'( "h: x",12x,"n",12x,
     &             1000(a," n "))')
     &             ( trim(chpchl_n(i,2)), i = 1, np_mxang ) ! frtati 2021/10/05 np -> np_mxang

      else
       if(iloopmode .eq. 0) then
        write(iot,'( "h: x",12x,"n",12x,
     &          1000(a,
     &          " ny",i0," dy",i0,"=[y",i0,"*y",i0,"] "))')
     &          ( trim(chpchl_n(i,2))
     &          ,i+np_mxang,i,i,i+np_mxang, i = 1, np_mxang ) ! frtati 2021/10/05
       else
        write(iot,'( "h: x",12x,"n",12x,
     &          1000(a,
     &          " ny",i4.4," dy",i4.4,"=[y",i4.4,"*y",i4.4,"] "))')
     &          ( trim(chpchl_n(i,2))
     &          ,i+np_mxang,i,i,i+np_mxang, i = 1, np_mxang ) ! frtati 2021/10/05
       end if
      end if

      write(iot,'( "#  act-lower    act-upper",1x,
     &             1000(2x,a11,"  r.err "))')
     &            ( chp_n(i), i = 1, nsame ) ! frtati 2021/10/05

      tott_n(:,:) = 0.d+0

      nactstepi = 1
      do iacti = 1, nact, nactstepi
        iact = iacti

        if( iacti .ne. nact ) then

           write(iot,'(1p2e13.4,1000(1pe13.4,0pf8.4))')
     &           actb(iacti),actb(iacti+nactstepi),
     &            ((((((tr(ipi+iploop-1,iei+ieloop-1,iti+itloop-1,
     &                     iri+irloop-1,iacti+iactloop-1,k),k=1,2),
     &                     iploop=1,npstepi),ieloop=1,nestepi),
     &                     itloop=1,ntstepi),irloop=1,nrstepi),
     &                     iactloop=1,nactstepi)

         elseif( iacti .eq. nact ) then

           trchk = 0.d0

           i = 0
           do iactloop=1,nactstepi
            do irloop=1,nrstepi
              do itloop=1,ntstepi
                do ieloop=1,nestepi
                  do iploop=1,npstepi
                    i = i + 1
                    vn = tr(ipi+iploop-1,iei+ieloop-1,
     &                  iti+itloop-1,iri+irloop-1,iacti+iactloop-1,1)

                    trchk = trchk + vn

                  end do
                end do
              end do
            end do
          end do

          if( trchk .gt. 0.d0 ) then
             write(iot,'("#   overflow",14x,
     &                     1000(1pe13.4,0pf8.4))')
     &            ((((((tr(ipi+iploop-1,iei+ieloop-1,iti+itloop-1,
     &                     iri+irloop-1,iacti+iactloop-1,k),k=1,2),
     &                     iploop=1,npstepi),ieloop=1,nestepi),
     &                     itloop=1,ntstepi),irloop=1,nrstepi),
     &                     iactloop=1,nactstepi)

           endif

         endif

         i = 0
         do iactloop=1,nactstepi
           do irloop=1,nrstepi
              do itloop=1,ntstepi
                do ieloop=1,nestepi
                  do iploop=1,npstepi
                    i = i + 1
                    vn = tr(ipi+iploop-1,iei+ieloop-1,
     &                   iti+itloop-1,iri+irloop-1,iacti+iactloop-1,1)

                    tott_n(i,1) = tott_n(i,1) + vn

                  end do
                end do
              end do
            end do
          end do

      end do

      do i = 1, nsame
         tott_n(i,1) = tott_n(i,1)
      end do

      call pstarreg_sumover_getput(0,m,iax,
     &    npstepi,nestepi,ntstepi,nrstepi,nactstepi,
     &    ipi,iei,iti,iri,1,maxtott_n,tott_n)  ! iacti = 1

      write(iot,'(/"#   sum over",14x,1000(1pe13.4,0pf8.4))')
     &                     (tott_n(i,1),tott_n(i,2),i=1,nsame)
