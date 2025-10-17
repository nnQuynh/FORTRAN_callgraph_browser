      if ( iteps(m) .ne. 2 ) then
        write(iot,'( "h: n",12x,"x",12x,
     &               1000(a," n  "))')
     &             ( trim(chpchl_n(i,1)) , i = 1, np_mxang )

      else

       if(iloopmode .eq. 0) then
        write(iot,'( "h: n",12x,"x",12x,
     &         1000(a
     &         ," ny",i0," dy",i0,"=[y",i0,"*y",i0,"] "))')
     &         ( trim(chpchl_n(i,1))
     &         ,i+np_mxang,i,i,i+np_mxang, i = 1, np_mxang )
       else
        write(iot,'( "h: n",12x,"x",12x,
     &         1000(a
     &         ," ny",i4.4," dy",i4.4,"=[y",i4.4,"*y",i4.4,"] "))')
     &         ( trim(chpchl_n(i,1))
     &         ,i+np_mxang,i,i,i+np_mxang, i = 1, np_mxang )
       end if
      end if

      write(iot,'( "#  let-lwr      let-upr",3x,
     &             1000(2x,a11,"  r.err "))')
     &     ( chp_n(i), i = 1, nsame )

      tott_n(:,:) = 0.d+0

      samewtt(:) = 0.0d0
      tr_n(:,:) =0.0d0


      if(itunt(m).ge.13) then ! yf(y) or yd(y) mode

        totnorm_n(:)=0.0

        i = 0
        do irloop=1,nrstepi
        do iploop=1,npstepi
          i = i + 1
          do ie=1,ne
             if(eb(ie).eq.0.0) then ! lower boundary is 0
               totnorm_n(i)=totnorm_n(i)
     &         + tr(ipi+iploop-1,ie,iri+irloop-1,1)
     &         *(eb(ie+1)-eb(ie))/(eb(ie+1)/2.0d0)
             else ! normal case
               totnorm_n(i)=totnorm_n(i)
     &         + tr(ipi+iploop-1,ie,iri+irloop-1,1)
     &         *(eb(ie+1)-eb(ie))/sqrt(eb(ie+1)*eb(ie))
            endif
            tott_n(i,1)=tott_n(i,1)+tr(ipi+iploop-1,ie,
     &                       iri+irloop-1,1)
     &      *(eb(ie+1)-eb(ie))

            tott_n(i,2)=tott_n(i,2)
     &      +(tr(ipi+iploop-1,ie,
     &                       iri+irloop-1,1)
     &      *tr(ipi+iploop-1,ie,
     &                       iri+irloop-1,2))**2
     &      *(eb(ie+1)-eb(ie))**2

          enddo
          if(totnorm_n(i).eq.0.0) totnorm_n(i)=1.0d0
          if(tott_n(i,1) .ne. 0.0) then
            tott_n(i,2) = sqrt(tott_n(i,2))/tott_n(i,1)
          else
            tott_n(i,2) = 0.0d0
          endif

          tott_n(i,1) = tott_n(i,1)/totnorm_n(i)
        enddo
        enddo

      else ! conventional mode
       totnorm_n(:)=1.0d0

      endif

      do iei = 1, ne, nestepi

        i = 0
        do irloop=1,nrstepi
          do ieloop=1,nestepi
            do iploop=1,npstepi
              i = i + 1

              if(itunt(m).ge.13) then
                tr_n(i,1) = tr(ipi+iploop-1,iei+ieloop-1,iri+irloop-1,1)
     &                    / totnorm_n(i)
!     &                    * (eb(iei+1)-eb(iei))
              else
                tr_n(i,1) = tr(ipi+iploop-1,iei+ieloop-1,iri+irloop-1,1)
     &                  /totnorm_n(i)
              endif
              tr_n(i,2) = tr(ipi+iploop-1,iei+ieloop-1,iri+irloop-1,2)

             if( itunt(m) .le. 12 ) then
              vn = tr(ipi+iploop-1,iei+ieloop-1,iri+irloop-1,1)
     &           /totnorm_n(i) * ew(iei+ieloop-1)

              tott_n(i,1) = tott_n(i,1) + vn

              samewtt(i) = samewtt(i) + ew(iei+ieloop-1)
             endif

            end do
          end do
        end do

        write(iot,'(1p2e13.4,1000(1pe13.4,0pf8.4))')
     &        eb(iei),eb(iei+1),
     &       (tr_n(i,1),tr_n(i,2),i=1,nsame)

      end do

      do i = 1, nsame
        if( itunt(m) .ne.  5 .and. itunt(m) .ne. 6 .and.
     &      itunt(m) .ne. 11 .and. itunt(m) .ne. 12 .and.
     &      itunt(m) .ne. 13 .and. itunt(m) .ne. 14 ) then

          tott_n(i,1) = tott_n(i,1) / samewtt(i)
        end if
      end do

      if(itunt(m).le.12) then
        call pletreg_sumover_getput(0,m,iax,
     &    npstepi,nestepi,nrstepi,
     &    ipi,1,iri,maxtott_n,tott_n)  ! iei = 1
      endif

      zero=0.0  ! T.Sato 2025/03/13 no uncertainty unit = 13 & 14
        if(itunt(m).eq.13) then
               write(iot,'(/"#          Frequency mean ",
     &         1000(1pe13.4,0pf8.4))') 
     &         (tott_n(i,1),zero,i=1,nsame)
        elseif(itunt(m).eq.14) then
               write(iot,'(/"#               Dose mean ",
     &         1000(1pe13.4,0pf8.4))')
     &         (tott_n(i,1),zero,i=1,nsame)
        else
          write(iot,'(/"#   sum over in cm or MeV ",
     &         1000(1pe13.4,0pf8.4))')
     &         (tott_n(i,1),tott_n(i,2),i=1,nsame)
        end if

! T.Sato 2025/03/12 for output normalization factor for restart calculation
               write(iot,'("#    Normalization factor ",
     &         1000(1pe13.4,8x))') (totnorm_n(i),i=1,nsame)
