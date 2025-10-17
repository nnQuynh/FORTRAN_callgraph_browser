c
      call getsamepage_iloopmode(m,iax,iloopmode)
c
      nmstepi = 1
      npstepi = 1
      nestepi = 1
      nxstepi = 1
      nystepi = 1
      nzstepi = 1

      if(iloopmode .eq. 0) then
        npstepi = np + np_0
        nsame = npstepi
      elseif(iloopmode .eq. 1) then
        nestepi = ne + ne_0
        nsame = nestepi
c reg
      elseif(iloopmode .eq. 3) then
        nxstepi = nx + nx_0
        nsame = nxstepi
      elseif(iloopmode .eq. 4) then
        nystepi = ny + ny_0
        nsame = nystepi
      elseif(iloopmode .eq. 5) then
        nzstepi = nz + nz_0
        nsame = nzstepi
c tet
c let
      elseif(iloopmode .eq. 14) then
        nestepi = ne + ne_0
        nsame = nestepi
      else
        npstepi = np + np_0
        nsame = npstepi
      endif
