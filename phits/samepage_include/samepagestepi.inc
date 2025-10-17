c
      call getsamepage_iloopmode(m,iax,iloopmode)
c
      nmstepi = 1
      npstepi = 1
      nestepi = 1
      nrstepi = 1
      nxstepi = 1
      nystepi = 1
      nzstepi = 1
      nastepi = 1
      ntstepi = 1
      ne1stepi = 1
      ne2stepi = 1
      nactstepi = 1

      if(iloopmode .eq. 0) then
        npstepi = np + np_0
        nsame = npstepi
      elseif(iloopmode .eq. 1) then
        nestepi = ne + ne_0
        nsame = nestepi
      elseif(iloopmode .eq. 2) then
c reg
        nrstepi = nr + nr_0
        nsame = nrstepi
      elseif(iloopmode .eq. 3) then
        nxstepi = nx + nx_0
        nsame = nxstepi
      elseif(iloopmode .eq. 4) then
        nystepi = ny + ny_0
        nsame = nystepi
      elseif(iloopmode .eq. 5) then
        nzstepi = nz + nz_0
        nsame = nzstepi
      elseif(iloopmode .eq. 6) then
        nrstepi = nr + nr_0
        nsame = nrstepi
c tet
      elseif(iloopmode .eq. 7) then
        nrstepi = nr + nr_0
        nsame = nrstepi
      elseif(iloopmode .eq. 8) then
        nastepi = na + na_0
        nsame = nastepi
      elseif(iloopmode .eq. 9) then
        nastepi = na + na_0
        nsame = nastepi
      elseif(iloopmode .eq. 10) then
        ntstepi = nt + nt_0
        nsame = ntstepi
      elseif(iloopmode .eq. 12) then
        ne1stepi = ne1 + ne1_0
        nsame = ne1stepi
      elseif(iloopmode .eq. 13) then
        ne2stepi = ne2 + ne2_0
        nsame = ne2stepi
c let
      elseif(iloopmode .eq. 14) then
        nestepi = ne + ne_0
        nsame = nestepi
      elseif(iloopmode .eq. 15) then
        nactstepi = nact + nact_0
        nsame = nactstepi
      else
        npstepi = np + np_0
        nsame = npstepi
      endif
