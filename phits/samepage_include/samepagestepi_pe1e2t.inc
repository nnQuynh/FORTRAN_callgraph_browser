c
      call getsamepage_iloopmode(m,iax,iloopmode)
c
      npstepi = 1
      ntstepi = 1
      ne1stepi = 1
      ne2stepi = 1

      if(iloopmode .eq. 0) then
        npstepi = np + np_0
        nsame = npstepi
c reg
c tet
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
      else
        npstepi = np + np_0
        nsame = npstepi
      endif
