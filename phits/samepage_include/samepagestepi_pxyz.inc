c
      call getsamepage_iloopmode(m,iax,iloopmode)
c
      npstepi = 1
      nxstepi = 1
      nystepi = 1
      nzstepi = 1

      if(iloopmode .eq. 0) then
        npstepi = np + np_0
        nsame = npstepi
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
      else
        npstepi = np + np_0
        nsame = npstepi
      endif
