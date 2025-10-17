c step1
      samepage_flag = .false.
      nmnr_flag = .false.

      np_n = 0
      nm_n = 0
      na_n = 0
      nt_n = 0
      nr_n = 0
      nx_n = 0
      ny_n = 0
      nz_n = 0
      ne_n = 0
      ne1_n = 0
      ne2_n = 0
      nact_n = 0

      chld(1:6)  = (/'l  ','dr ','ub ','mg ','qrr','vbb'/)
      chld(7:12) = (/'lr ','db ','ug ','mrr','qbb','v  '/)
      chld(13:18)= (/'lb ','dg ','urr','mbb','q  ','vr '/)
      chld(19:24)= (/'lg ','drr','ubb','m  ','qr ','vb '/)
      chld(25:30)= (/'lrr','dbb','u  ','mr ','qb ','vg '/)
      chld(31:36)= (/'lbb','d  ','ur ','mb ','qg ','vrr'/)
      chmd(1:6)  = (/'l3   ','d4r  ','u5b  ','m6g  ','q7rr ','v8bb '/)
      chmd(7:12) = (/'l9r  ','d10b ','u11g ','m12rr','q13bb','v14  '/)
      chmd(13:18)= (/'l3b  ','d4g  ','u5rr ','m6bb ','q7   ','v8r  '/)
      chmd(19:24)= (/'l9g  ','d10rr','u11bb','m12  ','q13r ','v14b '/)
      chmd(25:30)= (/'l3rr ','d4bb ','u5   ','m6r  ','q7b  ','v8g  '/)
      chmd(31:36)= (/'l9bb ','d10  ','u11r ','m12b ','q13g ','v14rr'/)

