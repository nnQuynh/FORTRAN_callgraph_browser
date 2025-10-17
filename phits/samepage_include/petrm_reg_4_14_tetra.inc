c
c nais samepage
      if ( iteps(m) .ne. 2 ) then

               write(iot,'( "h:   x      n",5x,"n",10x,
     &         1000(a," n  "))')
     &         ( trim(chpchm_n(i,1)) , i = 1, np_mxang )

       else
             if(iloopmode .eq. 0) then  
               write(iot,'( "h:   x      n",5x,"n",10x,
     &         1000(a,
     &         " ny",i0," dy",i0,"=[y",i0,"*y",i0,"] "))')
     &         ( trim(chpchm_n(i,1))
     &         ,i+np_mxang,i,i,i+np_mxang, i = 1, np_mxang )
             else
               write(iot,'( "h:   x      n",5x,"n",10x,
     &         1000(a,
     &         " ny",i4.4," dy",i4.4,"=[y",i4.4,"*y",i4.4,"] "))')
     &         ( trim(chpchm_n(i,1))
     &         ,i+np_mxang,i,i,i+np_mxang, i = 1, np_mxang )
             end if
      end if

      write(iot,'( "#    num    tetra   volume    ",
     &             1000(2x,a11,"  r.err "))')
     &             ( chp_n(i), i = 1, nsame )

      tott_n(:,:) = 0.0d+0

      samewtt(:) = 0.0d0

      do iri = 1, nr, nrstepi

          ireg = lr(iri)

          write(iot,'(i8,1x,i8,1pe13.4,1000(1pe13.4,0pf8.4))') 
     &      iri, ireg, vl(iri),
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
                     if( itunt(m) .ne. 4 .and. itunt(m) .ne. 14 ) then
                           vn = vn * vl(iri+irloop-1)
                     end if
                            tott_n(i,1) = tott_n(i,1) + vn

                            samewtt(i) = samewtt(i) + vl(iri+irloop-1)

                          end do
                        end do
                      end do
                    end do
            end do

       end do

       if( itunt(m) .ne. 4 .and. itunt(m) .ne. 14 ) then
         do i = 1, nsame
               tott_n(i,1) = tott_n(i,1) / samewtt(i)
         end do
       end if

      call pwwgreg_sumover_getput(0,m,iax,
     &    npstepi,nestepi,ntstepi,nrstepi,nmstepi,
     &    ipi,iei,iti,1,imi,maxtott_n,tott_n)  ! iri

       if( itunt(m) .ne. 4 .and. itunt(m) .ne. 14 ) then
         if(iloopmode .eq. 0) then
           write(iot,'(/"#   sum over ",1pe13.4,
     &               1000(1pe13.4,0pf8.4))')
     &               samewtt(1), (tott_n(i,1),tott_n(i,2),i=1,nsame)
         else
           write(iot,'(/"#   sum over ",13x,
     &               1000(1pe13.4,0pf8.4))')
     &               (tott_n(i,1),tott_n(i,2),i=1,nsame)
           write(iot,'(/"#   sum over ",5x, "wtt",5x,
     &               1000(1pe13.4,8x))')
     &               (samewtt(i),i=1,nsame)
         end if
      else
           write(iot,'(/"#   sum over ",13x,
     &               1000(1pe13.4,0pf8.4))')
     &               (tott_n(i,1),tott_n(i,2),i=1,nsame)
      end if
