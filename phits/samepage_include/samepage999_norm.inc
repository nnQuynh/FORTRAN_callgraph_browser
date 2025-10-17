c end_pro
      if(samepage_flag) then
        deallocate(tott_n)
        deallocate(samewtt)
        deallocate(totnorm_n,tr_n)
        deallocate(chp_n)
        deallocate(chpchl_n,chpchm_n)
      end if
      if(nmnr_flag) then
        deallocate(dum1_n,dum2_n)
        deallocate(lng1_n,lng2_n)
      end if
