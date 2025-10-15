      subroutine jxsdir(jc,iom,ierr)
      use GGMARRAYMOD !2020ASTOM
      use moddas_material
      implicit real*8 (a-h,o-z)
      include 'param.inc'
      include 'ggsparam.inc'
      include 'ggmparam.inc'
      include 'err.inc'
! T.Sato 2018/08/28, default value of hdpth
      common /paran/ icfn(100), ilfn(100), chfn(100)
      character chfn*200
      common /kmat1o/ iom1, iom2, iom3
      common /dircha/ idirch
      character yen*1
      common /kmat1k/ kmatd(kvlmax), kmate(kvlmax) ! frtati 2022/03/28
      integer, allocatable :: inlibflg(:) ! frtati 2023/3/6
      common /nlibcom/ nlibdef ! T.Sato 2023/12/26, read from natural_abundance.dat
      character*2 nlibdef
!      character*2 :: nlibdef = '50' ! frtati 2021/12/17 20MeV neutron default
      character ha*1,hb*30,hc*181,hd*10,hl*281,hn*13,hs*10,ht*10,hz*10 ! T.Sato 2023/08/22
      yen  = char(92)
      hdpth=chfn(1)(1:ilfn(1))//'/XS' ! default value of datapath in xsdir
      allocate(inlibflg(mxe)) ! frtati 2023/3/6
      iuo = iom1
      jc = 0
      nx = 0
      do 90 ie = 1, nxsc
         if( nty(ie) .gt. 0 ) goto 90
   10    if( nx .eq. nxsc ) goto 92
         hc = ' '
         call zaid(2,hc,ixc(1,ie))
         do 20 j = 11, 181
   20    hc(j:j) = char(mod(ixc((j+1)/3,ie)/256**mod(j+1,3),256))
         hz = hc(1:10)
         ha = hc(10:10)
         nt = index(htn,ha)
         do 30 je = 1, mxe
            if( nty(je) .ne. 0 ) goto 30
            call zaid(2,ht,ixl(1,je))
            if( ht(1:7) .ne. hz(1:7) ) goto 30
            if( ht(8:9) .ne. hz(8:9) .and.
     &          ht(8:9) .ne. ' ' ) goto 30
            if( ht(10:10) .eq. ha ) goto 40
            if( nt .eq. 1 .and.
     &          ht(10:10) .eq. 'd' ) nty(je) = -1
            if( nty(je) .lt. 0 ) goto 50
   30    continue
         if(ixc(61,ie).gt.0)
     &      write(iom1,'(
     &      '' ## warning. table on xs'',i3,
     &      '' card is not used in this problem.'')')
     &      ixc(61,ie)
         ixc(61,ie) = -abs(ixc(61,ie))
         goto 90
   40    call zaid(1,hz,ixl(1,je))
         nty(je) = nt
         nx = nx + 1
         if( kdata(hc(1:6)) .ne. 2 ) goto 50
         read(hc,'(bn,i6)') iz
         call nxtsym(hc,' ',11,it,iu,1)
         if( iu .eq. 0 ) goto 50
         hb = hc(it:iu)
         read(hb,'(bn,e30.0)') a
         do 48 im = 1, mix
   48    if( iza(im) .eq. iz ) awc(im) = a
   50    if( je .eq. ie ) goto 90
         do 60 i = 1, 61
            j = ixc(i,je)
            ixc(i,je) = ixc(i,ie)
   60       ixc(i,ie) = j
         if( je .gt. ie .and. ixc(61,ie) .gt. 0 ) goto 10
   90 continue
   92 if( mxe .gt. 0 .and. nx .eq. mxe ) goto 490
      jc = 1
      hl = xsdir
      if( inqire(hl) .eq. 1 ) goto 95
      if( idirch .eq. 0 ) then
         if( leng(hdpath) .gt. 0 )
     &         hl = hdpath(1:leng(hdpath))//'/'//xsdir
      else
         if( leng(hdpath) .gt. 0 )
     &         hl = hdpath(1:leng(hdpath))//yen//xsdir
      end if
      if( inqire(hl) .eq. 1 ) goto 95
         hdpath = hdpth
      if( idirch .eq. 0 ) then
         if( leng(hdpth) .gt. 0 )
     &         hl = hdpth(1:leng(hdpth))//'/'//xsdir
      else
         if( leng(hdpth) .gt. 0 )
     &         hl = hdpth(1:leng(hdpth))//yen//xsdir
      end if
      if( inqire(hl) .eq. 1 ) goto 95
       write(ErrCha,'(''Error: cannot find file(7): '',a)')
     & chfn(7)(1:ilfn(7))
       ErrID = 'L:1035/R:jxsdir/F:ggm02.f' !E02_001_001
       call ErrWriteIO(ErrID,ErrCha,iom)
       call ErrWrite(ErrID,ErrCha)
         ierr = 1
         return
   95 jc = 2
      open(iud,file=hl,status='old')
      rewind iud
   96 read(iud,'(a132)',err=98,iostat=ios) hc
      if( ios .eq. -1 ) goto 98
      if( hc(131:132) .eq. ' ' ) goto 96
      write(iom,'(a70)')hc(1:70)
         write(ErrCha,'(''Error: xsdir file data beyond column 131.'')')
         ErrID = 'L:1059/R:jxsdir/F:ggm02.f'
         call ErrWrite(ErrID,ErrCha)
         ierr = 1
         return
   98 rewind iud
  100    read(iud,'(a130)',iostat=ios) klin
         if( ios .eq. -1 ) goto 160
         hn = klin(1:13)
         call nxtsym(hn,' ',1,i,j,2)
         if( hn .eq. 'atomic weight' ) goto 105
         if( hn(1:8) .ne. 'datapath' ) goto 100
         call nxtsym(klin,' =',9,i,j,0)
         if( j .ne. 0 .and. klin(i:i) .eq. '=' )
     &   call nxtsym(klin,' =',j+1,i,j,0)
         if( j .ne. 0 ) hdpth = klin(i:j)
         goto 100
  105    iu = 0
  110    read(iud,'(a130)',iostat=ios) klin
         if( ios .eq. -1 ) goto 170
  120    call nxtsym(klin,' ',iu+1,it,iu,1)
         if( iu .eq. 0 ) goto 110
         if( kdata(klin(it:iu)) .eq. 0 ) goto 170
         hb = klin(it:iu)
         read(hb,'(bn,i30)') i
  130    call nxtsym(klin,' ',iu+1,it,iu,1)
         if( iu .ne. 0 ) goto 140
         read(iud,'(a130)',iostat=ios) klin
         if( ios .eq. -1 ) goto 360
         goto 130
  140    hb = klin(it:iu)
         read(hb,'(bn,e30.0)') a
      do 150 j = 1, mix
  150    if( iza(j) .eq. i ) awc(j) = a
         goto 120
  160    continue
         write(iom,'(''Error: cross-section directory file '',
     &               '' has no atomic weights table.'')')
         ierr = 1
         return
  170    if( mxe .eq. 0 ) goto 470
         if( klin .eq. 'directory' ) goto 190
         rewind iud
  180    read(iud,'(a130)',iostat=ios) klin
         if( ios .eq. -1 ) goto 360
         call nxtsym(klin,' ',1,i,j,2)
         if( klin .ne. 'directory' ) goto 180
  190    na = 0
         nb = 0
         hl = ' '
         hd = ' '
  200    iu = 0
         read(iud,'(a130)',iostat=ios) klin
         if( ios .eq. -1 ) goto 370
  210    if( na .eq. 0 ) call nxtsym(klin,' ',iu+1,it,iu,1)
         if( na .ne. 0 ) call nxtsym(klin,' ',iu+1,it,iu,0)
         if( iu .eq. 0 ) goto 240
         if( klin(it:iu) .eq. 'obsolete' ) nb = 1
         if( klin(it:iu) .eq. 'obsolete' .or.
     &       klin(it:iu) .eq. '+' ) goto 200
         na = na + 1
         if( na .gt. 1 ) goto 220
         i = index(klin(it:iu),'.')
         if( i .eq. 0 ) goto 190
         hl(8-i:10) = klin(it:iu)
         ic = 10
         goto 210
  220    if( na .gt. 10 .and. klin(it:iu) .ne. 'ptable' ) goto 230
         hl(ic+2:ic+iu-it+2) = klin(it:iu)
         ic = ic + iu - it + 2
         goto 210
  230    i = index(klin(it:iu),'/')
         j = index(klin(min(iu,it+i):iu),'/')
         if( ( i .eq. 2 .or. i .eq. 3 ) .and.
     &       ( j .eq. 2 .or. j .eq. 3 ) ) hd = klin(it:iu)
         goto 210
  240    hz = hl(1:10)
         ha = hz(10:10)
         inlibflg = 0 ! frtati 2021/12/17
  250    do 260 je = 1, mxe
            if( nty(je) .gt. 0 .or. inlibflg(je).eq.1 ) goto 260 ! frtati 2021/12/17
            call zaid(2,hs,ixl(1,je))
            if( hs(1:7) .ne. hz(1:7) .or. hs(10:10) .ne. ha ) goto 260
            if( hs(8:9) .eq. hz(8:9) .or. hs(8:9) .eq. ' ' )  goto 290
            if( ha.eq.'c' .and. hz(8:9).eq.nlibdef ) then
              nty(je) = -10
              inlibflg(je) = 1
              goto 300
            end if
  260    continue
         if( kdr(1) .eq. 0 .or. ha .ne. 'c' ) goto 190
         do 270 je = 1, mxe
            if( nty(je) .ne. 0 ) goto 270
            call zaid(2,ht,ixl(1,je))
            if( ht(1:7) .ne. hz(1:7) .or. ht(10:10) .ne. 'd' ) goto 270
            if( ht(8:9) .eq. hz(8:9) .or. ht(8:9) .eq. ' ' )   goto 280
  270    continue
         goto 190
  280    nty(je) = -1
         goto 300
  290    nty(je) = index(htn,ha)
         call zaid(1,hz,ixl(1,je))
         nx = nx + 1
  300    do 310 i = 1, 10
            n = index('0123456789/',hd(i:i))
  310       if( n .ne. 0 ) kxd(je) = 11 * kxd(je) + n
         kxs(je) = nb
         do 320 i = 1, 61
  320       ixc(i,je) = 0
         call zaid(1,hl,ixc(1,je))
         do 330 i = 11, 181
  330       ixc((i+1)/3,je) = ixc((i+1)/3,je)
     &                           + ichar(hl(i:i)) * 256**mod(i+1,3)
         if( nx .eq. mxe ) goto 385
         goto 250
  360    write(iom,'(
     & ''Error: bad data in cross-section directory file'')')
         ierr = 1
         return
  370    do 380 ie = 1, mxe
            if( nty(ie) .ge. 0 .or. nty(ie).eq.-10 ) goto 380 ! frtati 2021/12/17
            call zaid(2,ht,ixl(1,ie))
            call zaid(2,hs,ixc(1,ie))
            call zaid(1,hs,ixl(1,ie))
           nty(ie) = 1
           nx = nx + 1
           write(iom1,*) '## warning. '//
     &     'continuous-energy cross-section table used for '//ht
  380    continue
  385    ie = 0
  390    if( ie .eq. mxe ) goto 470
         ie = ie + 1
         call zaid(2,ht,ixl(1,ie))
  400    do 410 je = ie + 1, mxe
            call zaid(2,hz,ixl(1,je))
  410       if( hz .eq. ht ) goto 420
         goto 390
  420    do 430 m = 1, mix
         do 430 i = 1, mipt
           if( lme(i,m) .eq. je ) lme(i,m) = ie
  430      if( lme(i,m) .gt. je ) lme(i,m) = lme(i,m) - 1
         do 435 m = 1, mix
            if( lmn(m) .eq. je ) lmn(m) = ie
  435       if( lmn(m) .gt. je ) lmn(m) = lmn(m) - 1
         do 460 ke = je, mxe - 1
         do 440 i = 1, 3
  440       ixl(i,ke) = ixl(i,ke+1)
         do 450 i = 1, 61
  450       ixc(i,ke) = ixc(i,ke+1)
            kxd(ke) = kxd(ke+1)
            kxs(ke) = kxs(ke+1)
  460       nty(ke) = nty(ke+1)
            mxe = mxe - 1
            nx = nx - 1
         goto 400
  470    if( nx .eq. mxe ) goto 490
         i1_erflg = 0
         do ie = 1, mxe
          if( nty(ie) .eq. 0 ) then
           call zaid(2,ht,ixl(1,ie))
           if( ht(10:10).eq.'c' ) then
            write(ErrCha,'(1x,a10,1x,
     &      ": no library for low-energy neutron")') ht
            ErrID = 'L:1326/R:jxsdir/F:ggm02.f'
            call ErrWriteIO(ErrID,ErrCha,iom)
            call ErrWrite(ErrID,ErrCha)
            ierr = 1
            return
           else
            write(*,'(1x,a10)') ht
            i1_erflg = 1
           end if
          endif
         end do
         if(i1_erflg.eq.1) then
           write(*,'("These nuclear data libraries are missing. ",
     &     "Physical models are used for corresponding ",
     &     "nuclear reactions.")')
         end if
         do ie = 1, mxe
          if( nty(ie) .eq. -10 ) then
            mn = 1
            do km = 1, mix
              if( km .ge. jmd(1+mn+1) ) mn = mn + 1
              lem_tem = km-jmd(1+mn)+1
              if( lme(1,km).eq.ie ) then
                t0_dmax = das_kmate(kmate(mn)+(lem_tem-1)*5+12)
                exit
              end if
            end do
            if( t0_dmax.gt.20.d0 ) then
              i1_erflg = 2
              call zaid(2,ht,ixl(1,ie))
              write(*,'(1x,a10)') ht
            end if
          endif
         enddo
         if(i1_erflg.eq.2) then
           write(*,'("These nuclear data libraries are missing. ",
     &     " JENDL-4.0 is used for these nuclei up to 20MeV.")')
         end if
         np = 0
         do 475 i = 1, mix
  475       if( lmn(i) .eq. ie ) izn(i) = 0
  480    continue
  490    do 500 ie = 1, mxe
            if( kxs(ie) .eq. 0 ) goto 500
            call zaid(2,ht,ixl(1,ie))
            write(iom1,*)
     &      '## warning. '//ht//
     &      ' is an obsolete table to be eliminated soon.'
  500    continue
         if( jc .ne. 0 ) close(iud)
      deallocate(inlibflg) ! frtati 2023/3/6
      return
      end
      