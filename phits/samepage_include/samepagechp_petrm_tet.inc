c
      call getsamepage_iloopmode(m,iax,iloopmode)
c

      npstep = 1
      nestep = 1
      nrstep = 1
      nxstep = 1
      nystep = 1
      nzstep = 1
      nastep = 1
      ntstep = 1
      ne1step = 1
      ne2step = 1
      nactstep = 1


      nsame = np + np_0
      np_mxang = nsame

      if(iloopmode .eq. 0) then
        nsame = np + np_0
        npstep = nsame
      elseif(iloopmode .eq. 1) then
        nsame = ne + ne_0
        nestep = nsame
      elseif(iloopmode .eq. 2) then
        nsame = nr + nr_0
        nrstep = nsame
      elseif(iloopmode .eq. 3) then
        nsame = nx + nx_0
        nxstep = nsame
      elseif(iloopmode .eq. 4) then
        nsame = ny + ny_0
        nystep = nsame
      elseif(iloopmode .eq. 5) then
        nsame = nz + nz_0
        nzstep = nsame
      elseif(iloopmode .eq. 6) then
        nsame = nr + nr_0
        nrstep = nsame
c tet
      elseif(iloopmode .eq. 7) then
        nsame = nr + nr_0
        nrstep = nsame
      elseif(iloopmode .eq. 8) then
        nsame = na + na_0
        nastep = nsame
      elseif(iloopmode .eq. 9) then
        nsame = na + na_0
        nastep = nsame
      elseif(iloopmode .eq. 10) then
        nsame = nt + nt_0
        ntstep = nsame
      elseif(iloopmode .eq. 12) then
        nsame = ne1 + ne1_0
        ne1step = nsame
      elseif(iloopmode .eq. 13) then
        nsame = ne2 + ne2_0
        ne2step = nsame
      endif

      np_mxang = nsame
      if ( itmxang(m).ne.0 .and. iabs(itmxang(m)).lt.
     &     np_mxang ) then
        np_mxang = iabs(itmxang(m))
      end if

c
c particle name
      if(iloopmode .eq. 0) then
        do i=1,np
          ichl = mod(i,ndatamax)
          if(ichl .eq. 0) then
            ichl = ndatamax
          endif
          if ( iteps(m) .ne. 2 ) then
            chpchl_n(i,1)= chp(i)(1:10)//"),hh0"//chld(ichl)//" "
            chpchl_n(i,2)= chp(i)(1:10)//"),h0"//chld(ichl)//" "
            chpchm_n(i,1)= chp(i)(1:10)//"),"//chmd(ichl)(1:4)//" "
            chpchm_n(i,2)= chp(i)(1:10)//"),"//chmd(ichl)(1:3)//" "
          else
            write(chpchl_n(i,1),'(a1,i0,a9,"),hh0",a3)')
     &            chp(i)(1:1),i,chp(i)(2:10), chld(ichl)
            write(chpchl_n(i,2),'(a1,i0,a9,"),h0",a3)')
     &            chp(i)(1:1),i,chp(i)(2:10), chld(ichl)
            write(chpchm_n(i,1),'(a1,i0,a9,"),",a4)')
     &            chp(i)(1:1),i,chp(i)(2:10), chmd(ichl)
            write(chpchm_n(i,2),'(a1,i0,a9,"),",a3)')
     &            chp(i)(1:1),i,chp(i)(2:10), chmd(ichl)
          endif
        enddo
        do i=1,np
          chp_n(i)= chb(i)//'  '//chp(i)(3:10)
        enddo
c energy
      elseif(iloopmode .eq. 1) then
        do i = 1, ne
          ichl = mod(i,ndatamax)
          if(ichl .eq. 0) then
            ichl = ndatamax
          endif
          if ( iteps(m) .ne. 2 ) then
            write(chpchl_n(i,1),'("y(emin",1pe11.3,"),hh0",a3)')
     &            eb(i), chld(ichl)
            write(chpchl_n(i,2),'("y(emin",1pe11.3,"),h0",a3)')
     &            eb(i), chld(ichl)
            write(chpchm_n(i,1),'("y(emin",1pe11.3,"),",a4)')
     &            eb(i), chmd(ichl)
            write(chpchm_n(i,2),'("y(emin",1pe11.3,"),",a3)')
     &            eb(i), chmd(ichl)
          else
            write(chpchl_n(i,1),'("y",i0,"(emin",1pe11.3,"),hh0",a3)')
     &            i,eb(i), chld(ichl)
            write(chpchl_n(i,2),'("y",i0,"(emin",1pe11.3,"),h0",a3)')
     &            i,eb(i), chld(ichl)
            write(chpchm_n(i,1),'("y",i0,"(emin",1pe11.3,"),",a4)')
     &            i,eb(i), chmd(ichl)
            write(chpchm_n(i,2),'("y",i0,"(emin",1pe11.3,"),",a3)')
     &            i,eb(i), chmd(ichl)
          endif
          write(chp_n(i),'(1pe11.3)') eb(i)
        enddo
c reg
      elseif(iloopmode .eq. 2) then
c********
c
c

c x
      elseif(iloopmode .eq. 3) then
c y
      elseif(iloopmode .eq. 4) then
cc z
      elseif(iloopmode .eq. 5) then
cc r
      elseif(iloopmode .eq. 6) then
c tet
      elseif(iloopmode .eq. 7) then
         do i = 1, nr
          ichl = mod(i,ndatamax)
          if(ichl .eq. 0) then
            ichl = ndatamax
          endif
          if ( iteps(m) .ne. 2 ) then
            write(chpchl_n(i,1),'("y(tet",i7,"),hh0",a3)')
     &            i, chld(ichl)
            write(chpchl_n(i,2),'("y(tet",i7,"),h0",a3)')
     &            i, chld(ichl)
            write(chpchm_n(i,1),'("y(tet",i7,"),",a4)')
     &            i, chmd(ichl)
            write(chpchm_n(i,2),'("y(tet",i7,"),",a3)')
     &            i, chmd(ichl)
          else
            write(chpchl_n(i,1),'("y",i0,"(tet",i7,"),hh0",a3)')
     &            i,i, chld(ichl)
            write(chpchl_n(i,2),'("y",i0,"(tet",i7,"),h0",a3)')
     &            i,i, chld(ichl)
            write(chpchm_n(i,1),'("y",i0,"(tet",i7,"),",a4)')
     &            i,i, chmd(ichl)
            write(chpchm_n(i,2),'("y",i0,"(tet",i7,"),",a3)')
     &            i,i, chmd(ichl)
          endif
          write(chp_n(i),'("tet ",i7)') i
        enddo

c rad deg
      elseif(iloopmode .eq. 8 .or. iloopmode .eq. 9) then

c time
      elseif(iloopmode .eq. 10) then
        if( ittty(m) .ne. 0 ) then
          do i = 1, nt
            ichl = mod(i,ndatamax)
            if(ichl .eq. 0) then
              ichl = ndatamax
            endif
            if ( iteps(m) .ne. 2 ) then
              write(chpchl_n(i,1),'("y(tmin",1pe11.3,"),hh0",a3)')
     &             tb(i), chld(ichl)
              write(chpchl_n(i,2),'("y(tmin",1pe11.3,"),h0",a3)')
     &              tb(i), chld(ichl)
              write(chpchm_n(i,1),'("y(tmin",1pe11.3,"),hh0",a3)')
     &             tb(i), chmd(ichl)
              write(chpchm_n(i,2),'("y(tmin",1pe11.3,"),h0",a3)')
     &              tb(i), chmd(ichl)
            else
              write(chpchl_n(i,1),'("y",i0,"(tmin",1pe11.3,"),hh0",a3)')
     &              i,tb(i), chld(ichl)
              write(chpchl_n(i,2),'("y",i0,"(tmin",1pe11.3,"),h0",a3)')
     &              i,tb(i), chld(ichl)
              write(chpchm_n(i,1),'("y",i0,"(tmin",1pe11.3,"),hh0",a3)')
     &              i,tb(i), chmd(ichl)
              write(chpchm_n(i,2),'("y",i0,"(tmin",1pe11.3,"),h0",a3)')
     &              i,tb(i), chmd(ichl)
            endif
            write(chp_n(i),'(1pe11.3)') tb(i)
          enddo
        endif
c let
      elseif(iloopmode .eq. 14) then
c act
      elseif(iloopmode .eq. 15) then
c mset
      elseif(iloopmode .eq. 16) then
        if( itmlp(m) .gt. 0 ) then
          do i = 1, nm
            ichl = mod(i,ndatamax)
            if(ichl .eq. 0) then
              ichl = ndatamax
            endif
            if ( iteps(m) .ne. 2 ) then
              write(chpchl_n(i,1),'("y(mset",i3,"),hh0",a3)')
     &              itmnt(m,i), chld(ichl)
              write(chpchl_n(i,2),'("y(mset",i3,"),h0",a3)')
     &              itmnt(m,i), chld(ichl)
              write(chpchm_n(i,1),'("y(mset",i3,"),",a4)')
     &              itmnt(m,i), chmd(ichl)
              write(chpchm_n(i,2),'("y(mset",i3,"),",a3)')
     &              itmnt(m,i), chmd(ichl)
            else
              write(chpchl_n(i,1),'("y",i0,"(mset",i3,"),hh0",a3)')
     &              i,itmnt(m,i), chld(ichl)
              write(chpchl_n(i,2),'("y",i0,"(mset",i3,"),h0",a3)')
     &              i,itmnt(m,i), chld(ichl)
              write(chpchm_n(i,1),'("y",i0,"(mset",i3,"),",a4)')
     &              i,itmnt(m,i), chmd(ichl)
              write(chpchm_n(i,2),'("y",i0,"(mset",i3,"),",a3)')
     &              i,itmnt(m,i), chmd(ichl)
            endif
            write(chp_n(i),'("mset",i3)') itmnt(m,i)
          enddo
        endif
      endif

      npstepi = npstep
      nestepi = nestep
      nrstepi = nrstep
      nxstepi = nxstep
      nystepi = nystep
      nzstepi = nzstep
      nastepi = nastep
      ntstepi = ntstep
      ne1stepi = ne1step
      ne2stepi = ne2step
      nactstepi = nactstep


