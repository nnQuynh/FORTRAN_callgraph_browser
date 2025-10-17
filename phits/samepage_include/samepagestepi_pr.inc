c
      call getsamepage_iloopmode(m,iax,iloopmode)
c
      npstepi = 1
      nrstepi = 1

      if(iloopmode .eq. 0) then
        npstepi = np + np_0
        nsame = npstepi
      elseif(iloopmode .eq. 2) then
c reg
        nrstepi = nr + nr_0
        nsame = nrstepi
c tet
c let
      else
        npstepi = np + np_0
        nsame = npstepi
      endif
