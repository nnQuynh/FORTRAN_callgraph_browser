c step1
      np_0 = 0
      ne_0 = 0
      nr_0 = 0
      nx_0 = 0
      ny_0 = 0
      nz_0 = 0
      ntet_0 = 0
      na_0 = 0
      nt_0 = 0
      ne1_0 = 0
      ne2_0 = 0
      nact_0 = 0
      np_n = np + np_0
      nm_n = nm
      ne_n = ne + ne_0
      nr_n = nr + nr_0
      nx_n = nx + nx_0
      ny_n = ny + ny_0
      nz_n = nz + nz_0
      na_n = na + na_0
      nt_n = nt + nt_0
      ne1_n = ne1 + ne1_0
      ne2_n = ne2 + ne2_0
      nact_n = nact + nact_0
      maxtott_n = max(np_n, nm_n, ne_n+1, nr_n+1,
     &                na_n+1, nt_n+1)
c samepage kind
      iloopmode = 0
c part
      npstep = np_n
c e
      nestep = 1
c reg
      nrstep = 1
c x
      nxstep = 1
c y
      nystep = 1
c z
      nzstep = 1
c angle
      nastep = 1
c time
      ntstep = 1
c e1
      ne1step = 1
c e2
      ne2step = 1
c act
      nactstep = 1
c mset
      nmstep = 1
c samepage counts
      nsame =  np_n

      nmnr_flag = .false.
      samepage_flag = .true.
      allocate(tott_n(maxtott_n,2))
      allocate(samewtt(maxtott_n))
      allocate(chp_n(maxtott_n))
      allocate(chpchl_n(maxtott_n,2))
      allocate(chpchm_n(maxtott_n,2))
      if(nr_n*nm_n >= 1) then
         nmnr_flag = .true.
      endif
      if(nmnr_flag) then
        allocate(dum1_n(nr_n*nm_n),dum2_n(nr_n*nm_n))
        allocate(lng1_n(nr_n*nm_n),lng2_n(nr_n*nm_n))
      endif

      tott_n(:,:) = 0.0
      samewtt(:) = 0.0
      chp_n(:) = ' '
      chpchl_n(:,:) = ' '
      chpchm_n(:,:) = ' '
      if(nmnr_flag) then
        dum1_n(:) = ' '
        dum2_n(:) = ' '
        lng1_n(:) = 0
        lng2_n(:) = 0
      end if
