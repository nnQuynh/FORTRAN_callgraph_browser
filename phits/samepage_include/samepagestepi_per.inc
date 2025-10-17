c
      call getsamepage_iloopmode(m,iax,iloopmode)
c
      npstepi = 1
      nestepi = 1
      nrstepi = 1

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
c tet
c let
      elseif(iloopmode .eq. 14) then
        nestepi = ne + ne_0
        nsame = nestepi
      else
        npstepi = np + np_0
        nsame = npstepi
      endif
