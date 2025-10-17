c
      call getsamepage_iloopmode(m,iax,iloopmode)
c
      npstepi = 1
      nestepi = 1
      nrstepi = 1
      nzstepi = 1
      nastepi = 1
      ntstepi = 1

      if(iloopmode .eq. 0) then
        npstepi = np + np_0
        nsame = npstepi
      elseif(iloopmode .eq. 1) then
        nestepi = ne + ne_0
        nsame = nestepi
c reg
      elseif(iloopmode .eq. 5) then
        nzstepi = nz + nz_0
        nsame = nzstepi
      elseif(iloopmode .eq. 6) then
        nrstepi = nr + nr_0
        nsame = nrstepi
c tet
      elseif(iloopmode .eq. 8) then
        nastepi = na + na_0
        nsame = nastepi
      elseif(iloopmode .eq. 9) then
        nastepi = na + na_0
        nsame = nastepi
      elseif(iloopmode .eq. 10) then
        ntstepi = nt + nt_0
        nsame = ntstepi
c let
      else
        npstepi = np + np_0
        nsame = npstepi
      endif
