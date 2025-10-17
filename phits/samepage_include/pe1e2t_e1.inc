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

      write(iot,'( "#  e-lower      e-upper",3x,
     &             1000(2x,a11,"  r.err "))')
     &             ( chp_n(i), i = 1, nsame )

      tott_n(:,:) = 0.0d+0

      samewtt(:) = 0.0

      ne1stepi = 1
      do ie1i = 1, ne1, ne1stepi
          write(iot,'(1p2e13.4,1000(1pe13.4,0pf8.4))')
     &      eb1(ie1i),eb1(ie1i+ne1stepi),
     &      (((((tr(ipi+iploop-1,ie1i+ie1loop-1,ie2i+ie2loop-1,
     &              iti+itloop-1,k),k=1,2),
     &      iploop=1,npstepi),ie1loop=1,ne1stepi),ie2loop=1,ne2stepi),
     &      itloop=1,ntstepi)

            i = 0
                    do itloop=1,ntstepi
                      do ie2loop=1,ne2stepi
                        do ie1loop=1,ne1stepi
                          do iploop=1,npstepi
                            i = i + 1
                            vn = tr(ipi+iploop-1,ie1i+ie1loop-1,
     &                              ie2i+ie2loop-1,iti+itloop-1,
     &                              1)

                            tott_n(i,1) = tott_n(i,1) + vn

                          end do
                        end do
                    end do
                end do

      end do

      call pdpst2reg_sumover_getput(0,m,iax,
     &    npstepi,ne1stepi,ne2stepi,ntstepi,
     &    ipi,1,ie2i,iti,maxtott_n,tott_n)  ! ie1i

      if( itunt(m) .ne.  0 ) then
        write(iot,'(/"#   sum over ",13x,
     &              1000(1pe13.4,0pf8.4))')
     &             (tott_n(i,1),tott_n(i,2),i=1,nsame)
      end if
