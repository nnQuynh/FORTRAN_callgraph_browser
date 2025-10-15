      subroutine setgg(iom,jo,ierr)
      use GGBANKMOD !FURUTA
      use LATDATAMOD
      use LAFDATAMOD !FURUTA20201127
      use moddas
      use moddas_character
      use moddas_ggs
      implicit real*8 (a-h,o-z)
      include 'param.inc'
      include 'ggsparam.inc'
      include 'err.inc'
      parameter ( ibmt = 40 )
      common /ccggg/  icgg
      common /ggmes/  iggcm
      common /cggmm/  ngstar, ngfini, ngfin0
      common /tcntl/  icntl, inucr
      common /inpec/  ititl, ipara, ibody, iregn, llarr, itby, itar
      common /inggs/  iog, igcel, ioa, igsuf, iob, igtrs
      common /celdb/  idsn(kvlmax), idtn(kvlmax)
      common /regda/  ichl(kvlmax), chsm(kvlmax), ichmx, iod
      character       chsm*10
      common /regdc/  idrg(kvlmax), idgr(kvmmax)
      common /regdf/  ioc, ifilt, ifil(kvlmax)
      common /regdu/  iuni(kvlmax)
      common /regde/  ichp(kvlmax), ilat(kvlmax), idct(kvlmax)
      common /regdm/  idmg(kvlmax)
      common /regdl/  ilike(kvlmax)
      common /celda/  deng(kvlmax)
      common /volreg/ dvol(kvlmax)
      common /tmpreg/ dtmp(kvlmax)
      common /regdd/  ivolm, iimpo
      integer nlat3,itetcl(10),itetvol,itetauto,itgchk
      common /tetf1/ nlat3,itetcl,itetvol,itetauto,itgchk
      real(8) tetsfac(10)
      common /tetf2/ tetsfac,ltfile,itfform,tfilename
      common /taliin/ rsouin, nzztin, nrgnin
      common /mpi00/ npe, me
      dimension vtrs(13)
      logical deqn5
      dimension bval(50)
      dimension nr(39)
      dimension ns(39)
      character ksf(50)*3
      dimension kfil(6)
      real*8 xangle, yangle, zangle
      integer inonzero
      real*8 rmatrix(3,3,3)
      common /geomemcom/ igeomem  ! T.Sato 2024/12/25
      integer,allocatable :: lattmp(:),lattmp2(:)
      common /ggcell/ icells, iobo
      character chcfg*100
      character chtrs*100
      common /paran/  icfn(100), ilfn(100), chfn(100)
      character       chfn*200
      logical   exex
               ierr  = 0
         if( icells .eq. 1 .or. icells .eq. 2 ) then
               chcfg = chfn(19)(1:ilfn(19))//'.cfg'
               ichcfg = ilfn(19) + 4
            if( icells .eq. 1 ) then
                  ierror = 0
                  inquire( file = chfn(19), exist = exex )
               if( exex .eqv. .false. ) then
                  ierror = ierror + 1
               end if
                  inquire( file = chcfg, exist = exex )
                  if( ierror .gt. 0 ) goto 999
               iod = 26
            end if
               ioj = 203
         end if
             if( icells .eq. 1 ) then
                read(ioj) igtrs, iobo
             else if( icells .eq. 2 ) then
                write(ioj) igtrs, iobo
             end if
         if( igtrs .gt. 0 ) then
               call moddas_allocate_dbl2( 17, 0, igtrs, trf )
            do j = 0, igtrs
               do i = 1, 14
                  trf(i,j) = 0.0d0
               end do
            end do
                  trf( 5,0) = 1.0
                  trf( 9,0) = 1.0
                  trf(13,0) = 1.0
            do j = 1, igtrs
               do i = 5, 13
                  trf(i,j) = huge
               end do
            end do
            if( igtrs .gt. 0 ) then
                  if( icells .ne. 1 ) then
                        rewind iob
                  else if( icells .eq. 1 ) then
                     if( iobo .eq. 1 ) then
                        iob = 71
                        open(iob,status='scratch',form='unformatted')
                     end if
                  end if
               do i = 1, igtrs
                  if( icells .eq. 1 ) then
                     read(ioj) itrs, ( vtrs(j), j = 1, 13 ), idtn(i)
                     if( iobo .eq. 1 ) then
                        write(iob) itrs, ( vtrs(j), j = 1, 13 )
                     end if
                  else
                     read(iob) itrs, ( vtrs(j), j = 1, 13 )
                  end if
                  if( icells .eq. 2 ) then
                     write(ioj) itrs, ( vtrs(j), j = 1, 13 ), idtn(i)
                  end if
                        trf( 1, i ) = - idtn(i)
                        if ( iabs(idnint(vtrs(13))) .eq. 2 ) then
                         if( itrs .ne. 0 ) then
                          zangle = vtrs(4) * pi / 180.
                          yangle = vtrs(5) * pi / 180.
                          xangle = vtrs(6) * pi / 180.
                         else
                          zangle = vtrs(4)
                          yangle = vtrs(5)
                          xangle = vtrs(6)
                         end if
                         if ( iabs(idnint(vtrs(13))) .ne. 1 ) then
                          inonzero = 0
                          do j = 7, 12
                           if ( vtrs(j) .ne. 0d0 ) inonzero = 1
                          end do
                          if ( inonzero .eq. 1 ) then
                           ErrCha = ''
                           MsgID = 'L:278/R:setgg/F:ggs00.f'
                           call ErrWrite(MsgID, ErrCha)
                          end if
                         end if
                        elseif ( iabs(idnint(vtrs(13))) .eq. 3 ) then
                         rmatrix(1,1,3)=1.0d0
                         rmatrix(1,2,3)=0.0d0
                         rmatrix(1,3,3)=0.0d0
                         rmatrix(2,1,3)=0.0d0
                         rmatrix(2,2,3)=1.0d0
                         rmatrix(2,3,3)=0.0d0
                         rmatrix(3,1,3)=0.0d0
                         rmatrix(3,2,3)=0.0d0
                         rmatrix(3,3,3)=1.0d0
                         do im=1,3 ! up to 3 matrix
                          iraxis=nint(vtrs(5+im*2)) ! data ID of 7, 9, 11
                          if(iraxis.ge.4) then
       write(*,*) "Error in [transform] definition"
                           MsgID = 'L:320/R:setgg/F:ggs00.f'
                           call ErrWrite(MsgID, ErrCha)
                           stop
                          endif
                          if(iraxis.ge.1) then ! iraxis = 1,2, or 3
                           if( itrs .ne. 0 ) then
                            arot = vtrs(6+im*2) * pi / 180.
                           else
                            arot = vtrs(6+im*2)
                           end if
                           if(iraxis.eq.1) then     ! rotate around X-axis
                            rmatrix(1,1,2)=1.0d0
                            rmatrix(1,2,2)=0.0
                            rmatrix(1,3,2)=0.0
                            rmatrix(2,1,2)=0.0
                            rmatrix(2,2,2)=cos(arot)
                            rmatrix(2,3,2)=-sin(arot)
                            rmatrix(3,1,2)=0.0
                            rmatrix(3,2,2)=sin(arot)
                            rmatrix(3,3,2)=cos(arot)
                           elseif(iraxis.eq.2) then ! rotate around Y-axis
                            rmatrix(1,1,2)=cos(arot)
                            rmatrix(1,2,2)=0.0
                            rmatrix(1,3,2)=sin(arot)
                            rmatrix(2,1,2)=0.0
                            rmatrix(2,2,2)=1.0d0
                            rmatrix(2,3,2)=0.0
                            rmatrix(3,1,2)=-sin(arot)
                            rmatrix(3,2,2)=0.0
                            rmatrix(3,3,2)=cos(arot)
                           else                    ! rotate around Z-axis
                            rmatrix(1,1,2)=cos(arot)
                            rmatrix(1,2,2)=-sin(arot)
                            rmatrix(1,3,2)=0.0
                            rmatrix(2,1,2)=sin(arot)
                            rmatrix(2,2,2)=cos(arot)
                            rmatrix(2,3,2)=0.0
                            rmatrix(3,1,2)=0.0
                            rmatrix(3,2,2)=0.0
                            rmatrix(3,3,2)=1.0d0
                           endif
                           rmatrix(:,:,1)=rmatrix(:,:,3) ! copy current matrix (=3) to original maxtrix (=1)
                           rmatrix(:,:,3)=0.0d0
                           do j=1,3
                            do k=1,3
                             do l=1,3
                             enddo
                            enddo
                           enddo
                          endif
                         enddo
                         if(vtrs(4).ne.0) then
                          if(idnint(vtrs(13)).eq.3) then ! X0=X0+XR-(R3R2R1)XR
                          endif
                         endif
                         vtrs(4)= rmatrix(1,1,3)
                         vtrs(5)= rmatrix(2,1,3)
                         vtrs(6)= rmatrix(3,1,3)
                         vtrs(7)= rmatrix(1,2,3)
                         vtrs(8)= rmatrix(2,2,3)
                         vtrs(9)= rmatrix(3,2,3)
                         vtrs(10)=rmatrix(1,3,3)
                         vtrs(11)=rmatrix(2,3,3)
                         vtrs(12)=rmatrix(3,3,3)
                        end if
                  do j = 1, 13
                        trf( 1 + j, i ) = vtrs(j)
                     if( j .ge. 4 .and. j .le. 12 ) then
                      if ( iabs(idnint(vtrs(13))) .eq. 1 ) then
                        if( itrs .ne. 0 )&
                       trf( 1 + j, i ) = cos( vtrs(j) * pi / 180. )
                      end if
                        if( abs(  trf(1+j,i)&
                         - anint(trf(1+j,i)) ) .gt. 1e-10 )&
                       trf( ltrf + 1, i ) = idtn(i)
                     end if
                     if ( j .eq. 13 )&
                         trf( 1 + j, i ) = dsign(1d0,vtrs(j))
                  end do
                     imxt = i
                     call trfdef(iom,jo,ierr,imxt)
               end do
                     if( ierr .ne. 0 ) then
                        iog = 20
                        io  = iog
                        open(io,form='formatted',status='scratch')
                        goto 999
                     end if
            end if
         end if
            if( icgg  .eq. 0 ) return
               iog = 20
               io  = iog
               open(io,form='formatted',status='scratch')
               mcmx = ( mdas / 2  - 1 ) * 8 + 1
               call moddas_allocate_cha(MAX_NUM_CHRG, chrg)
               mci = 0
               ngstar = mmmax
            if( icells .eq. 1 ) then
                  read(ioj)  igcel, igsuf, ifilt, nlat3
               do i = 1, igcel
               end do
               do i = 1, nlat3
                  read(ioj)  itetcl(i), tetsfac(i)
               end do
               if( ifilt .gt. 0 ) then
                  inquire( file = chfn(18), exist = exex )
                  if( exex .eqv. .false. ) then
                     goto 999
                  end if
                  ioc = 66
               end if
            else if( icells .eq. 2 ) then
                  write(ioj) igcel, igsuf, ifilt, nlat3
               do i = 1, igcel
               end do
               do i = 1, nlat3
                  write(ioj) itetcl(i), tetsfac(i)
               end do
            end if
               mxa    = igcel
               mxafs  = mxa
               mxtr   = igtrs
               nrgnin = mxa
            if( icells .eq. 1 ) goto 10000
               call moddas_allocate_int(igsuf, kst)
               nsc = 0
               mxj = 0
               rewind  ioa
         do i = 1, igsuf
                  kst(i) = 0
            read(ioa) idrf, idtr, idsf, igkst, ( bval(j), j = 1, igkst )
               mxj = mxj + 1
               n = 10
               if( idsf .eq. 1 .and. igkst .gt. 8 ) n = 10
               if( idsf .ge. 30 ) n = nr(idsf)
               nsc = nsc + n
            if( idsf .ge. 30 .and. idsf .le. 39 ) then
                  mxj = mxj + jn
                  kst(i) = jn
            end if
         end do
            ncomp = 0
            mlja  = 0
            mxit  = 0
            rewind iod
         k_all_total=0  ! T.Sato 2024/12/26
         do i = 1, mxa
! Nais_2024
            k_all_check = 0
            k_4_check = 0
            k_back = 0
! Nais_2024
            m1c = 0
            ic  = 0
            img = 0
            read(iod) (chrg(k:k),k=1,ichp(i))
  451       ic = ic + 1
            if( ic .gt. ichl(i) ) goto 461
               k = index('():#',chrg(ic:ic))
! Nais_2024 >>>
               k_all = 0
               if(k == 4) then
                  if(ic + 3 < ichl(i)) then
                    k_all = index(chrg(ic:),'#all ')
                  else
                    k_all = index(chrg(ic:),'#all')
                  end if
                  if(k_all == 1) then
                     k_all_check = 1
                     k_all_total=k_all_total+1
                     if(k_all_total.ge.2) then ! more than two cells include #all
                      MsgID = 'L:670/R:setgg/F:ggs00.f' !E06_001_001
                      call ErrWrite(MsgID, ErrCha)
                      stop
                     endif
                  else
                     k_4_check = 1
                  end if
               end if
! Nais_2024 <<<
               if( deqn5( chrg(ic:ic) ) ) k = 5
               if( k .eq. 5 .and. img .eq. 0 ) then
                     img = 1
                     ini = ic
               end if
               if( img .eq. 1 .and.&
                ( ( k .eq. 5 .and. ic .eq. ichl(i) ) .or.&
                  ( k .ne. 5 ) ) ) then
                     img = 0
                     m1c = m1c + 1
                  if( k .eq. 5 ) then
                     ifi = ic
                  else
                     ifi = ic - 1
                  end if
                     call onum(chrg(1:ifi),ini,ifi,cvvv,ierr)
                     if( ierr .ne. 0 ) goto 999
                     idsuf = abs( nint( cvvv ) )
! Nais_2024 >>>
                  if(k_back /= 4) then
! Nais_2024 <<<
                    do j = 1, igsuf
                     if( idsuf .eq. idsn(j) ) then
                           m1c = m1c + 2 * kst(j)
                        if( idct(i) .ne. 0 ) then
                           mxj = mxj + kst(j) + 1
                           nsc = nsc + max( 10, 5 * kst(j) ) ! 4*kst is occasionally insufficient
                        end if
                     end if
                    end do
                 end if
                 k_back = k
! Nais_2024 <<<
               end if
! Nais_2024 >>>
               if(k == 4 .and. k_all == 1) then
                  if(ic + 3 <  ichl(i)) then
                     ic = ic + 4
                  else
                     ic = ic + 3
                  endif
                  k_back = 0
                  goto 451
               end if
! Nais_2024 <<<
               if( k .ge. 1 .and. k .le. 4 ) then
                     m1c = m1c + 1
                     if( k .eq. 4 ) ncomp = ncomp + 1
                     k_back = k
               end if
               goto 451
  461       continue
! Nais_2024 >>>
            if(k_all_check == 1 .and. k_4_check == 1) then
             MsgID = 'L:781/R:setgg/F:ggs00.f' !E06_001_001
             call ErrWrite(MsgID, ErrCha)
            endif
! Nais_2024 <<<
! #all: get cells
            if(k_all_check == 1) then
                 do i_all = 1,igcel
                    if(i /= i_all) then
                      if(idmg(i_all) == -1) then
                      else
                        if(iuni(i_all) == 0) then
                           m1c = m1c + 2
                           ncomp = ncomp + 1
                        end if
                      end if
                   end if
                 end do
            end if
            mlja = mlja + m1c
            mxit = max( mxit, m1c )
         end do
               nlat = 0
            do i = 1, mxa
               if( iuni(i) .ne. 0 .or. ilat(i) .ne. 0 .or.&
                  idct(i) .ne. 0 ) junf = 1
               if( ilat(i) .ne. 0 ) nlat = nlat + 1
            end do
               if( ifilt .ne. 0 ) junf = 1
         if( junf .ne. 0 ) then
               call moddas_allocate_int(mxa, jun)
               call moddas_allocate_int2(1, 3, mxa, mfl)
               call moddas_allocate_int2(1, 3, mxa, mazp)
               llaf = 1 !FURUTA20201127 0 is avoided due to special fucntion
            do i = 1, mxa
               jun(i) = iuni(i)
            do j = 1, 3
               mazp(j,i) = 0
                mfl(j,i) = 0
            end do
            end do
         end if
               mlaf = 0
               nmzu = 0
               nmaz = 0
            if( ifilt .ne. 0 ) then
                  rewind  ioc
                  if(nlatind.gt.0)call REWIND_ldata
               do i = 1, mxa
                  if( ifil(i) .eq. 1 ) then
                        read(ioc) klat, jlat
                        mfl(1,i) = klat
                  else if( ifil(i) .eq. 6 ) then
                        mfl(1,i) = -llaf - mlaf
                        read(ioc) ( kfil(j), j = 1, 6 ), latot
                        call EXTEND_laf(abs(latot)+2) !FURUTA20201127
                    if(latot.gt.0)then !FURUTA20181009
                     do j = 1, latot
                        read(ioc) klat, jlat
                        laf(1,llaf+mlaf+2+j) = klat !FURUTA20201127
                     end do
                    else               !FURUTA20181009
                     latot=-latot
                     allocate(lattmp(latot))
                     if(nlatind.gt.0)then
                      ilatind=ilatind+1
                      latot2=ldata(ilatind)
                      allocate(lattmp2(latot2))
                      lattmp2(1:latot2)=ldata(ilatind+1:ilatind+latot2)
                      ilatind=ilatind+latot2
                     else
                      read(ioc)latot2
                      allocate(lattmp2(latot2))
                      read(ioc)lattmp2(1:latot2)
                     endif
                     call latdecode(latot2,lattmp2,latot,lattmp,&
                         icount,ierr)
                     do j=1,latot
                      laf(1,llaf+mlaf+2+j)=lattmp(j) !FURUTA20201127
                     enddo
                     deallocate(lattmp,lattmp2)
                     if(ierr.ne.0)then
                      goto 999
                     endif
                    endif              !FURUTA20181009
                     do k = 1, 6
                        if( mod(k,2) .ne. 0 ) then
                           laf(1+k/2,llaf+mlaf+1) = kfil(k) !FURUTA20201127
                        else if( mod(k,2) .eq. 0 ) then
                        end if
                     end do
                     mlaf = mlaf + ( latot + 2 ) !FURUTA20201127
                  end if
               end do
                  if( ierr .ne. 0 ) goto 999
            end if
            if( icells .eq. 2 ) then
                  rewind  ioa
               do i = 1, igsuf
               end do
            end if
10000       continue
            if( icells .eq. 1 ) then
                  rewind  ioa
               do i = 1, igsuf
               end do
            end if
            mlja0 = mlja + 2 * mxit * ncomp
            mxj0  = mxj
               igcnt = 0
            if( icells .eq. 1 ) then
               igcnt = 1
               read(ioj)  mxj, nljc
            end if
 2000 continue
            igcnt = igcnt + 1
         if( igcnt .eq. 1 ) then
            mlja = mlja + 2 * mxit * ncomp
         else
            mlja  = nljc + 1
         end if
            mlaj = ( 12 + 50 * junf ) * mxa + 50
            mtasks = 1
            nlse = 0
            call moddas_allocate_dbl(nsc+igeomem, scf) ! T.Sato 2024/12/25 change +2000 to +igeomem
            call moddas_allocate_dbl3(3, 7, nlat, vcl)
               call moddas_allocate_int(mxj, ksc)
               call moddas_allocate_int(mxa, nlv)
               call moddas_allocate_int(mxa*mtasks, lse)
               call moddas_allocate_int(mxa*junf, jun)
               call moddas_allocate_int2(1, 2, mxa*junf, lat)
               call moddas_allocate_int(mxa*junf, ktr)
               call moddas_allocate_int2(1, 3, mxa*junf, mfl)
               call moddas_allocate_int2(1, 3, mxa*junf, mazp)
               if( nmaz .ne. 0 ) call moddas_allocate_int(nmzu, mazu)
               call moddas_allocate_int(mxj, kst)
               call moddas_allocate_int(mxj, ksu)
               call moddas_allocate_int(mlja+igeomem, idna) ! T.Sato 2024/12/25 change +1000 to +igeomem
               call moddas_allocate_int(mxj+igeomem, idne)  ! T.Sato 2024/12/25 change +1000 to +igeomem
               call moddas_allocate_int(mxj, idns)
               call moddas_allocate_int(mxj, idnt)
               call moddas_allocate_int(mxafs+1, lca)
               call moddas_allocate_int(mlja+igeomem, lja)  ! T.Sato 2024/12/25 change +1000 to +geomem
               call moddas_allocate_int(mxj+1, lsc)
               call moddas_allocate_int(mxa+2, ncl)
               call moddas_allocate_int(mxj, nsfm)
               call moddas_allocate_int(mxj, ksm)
               call moddas_allocate_int(mxj, jtr)
            nlaj_bank  = ( mlaj + mxa ) * mtasks !FURUTA
            nlcaj_bank = ( mlja + 1 ) * mtasks   !FURUTA
            kdb  = 0
            mnax = mmmax
            mmmax  = mnax + 1
            ngfini = mnax + 1
            if( igcnt .eq. 1 ) ngfin0 = mnax
            call INIT_laf !FURUTA20201127
            lsc(1) = 0
            nlja = 0
            nljc = 0
            rewind iod
         do i = 1, mxa
               lca( i ) = nlja + 1
               ncl( i ) = idrg(i)
            read(iod) (chrg(k:k),k=1,ichp(i))
            ic  = 0
            img = 0
  551       ic = ic + 1
            if( ic .gt. ichl(i) ) goto 561
               k = index('():#',chrg(ic:ic))
! Nais_2024 >>>
               k_all = 0
               if(k == 4) then
                  if(ic + 3 < ichl(i)) then
                    k_all = index(chrg(ic:),'#all ')
                  else
                    k_all = index(chrg(ic:),'#all')
                  end if
               endif
! Nais_2024 <<<
               if( deqn5( chrg(ic:ic) ) ) k = 5
               if( k .eq. 5 .and. img .eq. 0 ) then
                     img = 1
                     ini = ic
               end if
               if( img .eq. 1 .and.&
                ( ( k .eq. 5 .and. ic .eq. ichl(i) ) .or.&
                  ( k .ne. 5 ) ) ) then
                     img = 0
                  if( k .eq. 5 ) then
                     ifi = ic
                  else
                     ifi = ic - 1
                  end if
                     call onum(chrg(1:ifi),ini,ifi,cvvv,ierr)
                     if( ierr .ne. 0 ) goto 999
                     nlja = nlja + 1
                     nljc = nljc + 1
                     lja( nlja ) = nint( cvvv )
                  if( abs( cvvv - anint(cvvv)).gt.&
                     5e-14 * abs(cvvv) ) then
                  end if
               end if
! Nais_2024 >>>
               if( k == 4 .and. k_all == 1 ) then
                 if(ic + 3 < ichl(i)) then
                    ic = ic + 4
                 else
                   ic = ic + 3
                 end if
                 lca( i )  = -abs( lca( i ) )
! #all: get cells
                 do i_all = 1,igcel
                    if(i /= i_all) then
                      if(idmg(i_all) == -1) then
                      else
                        if(iuni(i_all) == 0) then
                           nlja = nlja + 1
                           nljc = nljc + 1
                           lja( nlja ) = 1000000 + k
                           nlja = nlja + 1
                           nljc = nljc + 1
                           lja( nlja ) = idrg(i_all)
!
!                           write(*,'(1x,i8)') idrg(i_all)
!
                        end if
                      end if
                   end if
                 end do
                 go to 551
               end if
! Nais_2024 <<<
               if( k .ge. 1 .and. k .le. 4 ) then
                     nlja = nlja + 1
                     nljc = nljc + 1
                     lca( i )  = -abs( lca( i ) )
                     lja( nlja ) = 1000000 + k
               end if
               goto 551
  561       continue
         end do
            lca( mxa + 1 ) = nlja + 1
         if( junf .ne. 0 ) then
            do i = 1, mxa
               jun(i) = iuni(i)
            end do
               mlaf = 0
               nlat = 0
               if( ifilt .ne. 0 )then
                rewind  ioc
                if(nlatind.gt.0)call REWIND_ldata
               endif
               do i = 1, mxa
                  if( ifil(i) .eq. 1 ) then
                        read(ioc) klat, jlat
                        mfl(1,i) = klat
                        mfl(3,i) = jlat
                  else if( ifil(i) .eq. 6 ) then
                        mfl(1,i) = -llaf - mlaf
                        read(ioc) ( kfil(j), j = 1, 6 ), latot
                    if(latot.gt.0)then !FURUTA20181009
                     do j = 1, latot
                        read(ioc) klat, jlat
                        laf(1,llaf+mlaf+2+j) = klat !FURUTA20201127
                        laf(3,llaf+mlaf+2+j) = jlat !FURUTA20201127
                     end do
                    else               !FURUTA20181009
                     latot=-latot
                     allocate(lattmp(latot))
                     if(nlatind.gt.0)then
                      ilatind=ilatind+1
                      latot2=ldata(ilatind)
                      allocate(lattmp2(latot2))
                      lattmp2(1:latot2)=ldata(ilatind+1:ilatind+latot2)
                      ilatind=ilatind+latot2
                     else
                      read(ioc)latot2
                      allocate(lattmp2(latot2))
                      read(ioc)lattmp2(1:latot2)
                     endif
                     call latdecode(latot2,lattmp2,latot,lattmp,&
                         icount,ierr)
                     do j=1,latot
                      laf(1,llaf+mlaf+2+j)=lattmp(j) !FURUTA20201127
                     enddo
                     do j=1,latot
                      laf(3,llaf+mlaf+2+j)=0 !FURUTA20201127
                     enddo
                     deallocate(lattmp,lattmp2)
                     if(ierr.ne.0)then
                      goto 999
                     endif
                    endif              !FURUTA20181009
                     do k = 1, 6
                        if( mod(k,2) .ne. 0 ) then
                           laf(1+k/2,llaf+mlaf+1) = kfil(k) !FURUTA20201127
                        else if( mod(k,2) .eq. 0 ) then
                        end if
                     end do
                        mlaf = mlaf + ( latot + 2 ) !FURUTA20201127
                  end if
                  if( ilat(i) .ne. 0 ) then
                     lat(1,i) = ilat(i)
                     nlat = nlat + 1
                     lat(2,i) = nlat
                  end if
                  if( idct(i) .ne. 0 ) then
                     ktr( i ) = idct(i)
                  end if
               end do
         end if
               mxj = 0
            rewind  ioa
      do 200 i = 1, igsuf
            read(ioa) idrf, idtr, idsf, igkst, ( bval(j), j = 1, igkst )
! T.Sato 2025/03/11, avoid 0 for TRC
      if(idsf.eq.36.and.bval(8).eq.0.0) then
       ErrID = 'L:1397/R:setgg/F:ggs00.f'
       call ErrWrite(ErrID,ErrCha)
       bval(8)=1.0e-5
      endif
               mxj = mxj + 1
               nsfm( mxj ) = idsn(i)
            if( idrf .ne. 0 ) ksu( mxj ) = idrf
            if( idtr .ne. 0 ) then
               if( idtr .gt. 0 ) jtr( mxj ) =  idtr
               if( idtr .lt. 0 ) ksu( mxj ) = -idtr
            end if
               m0c = idsf
               m1c = idsf
            do j = 1, igkst
               scf( lsc( mxj ) + j ) = bval(j)
            end do
               ix = lsc( mxj )
            if( ( m1c .ge. 27 .and. m1c .le. 29 ) .or.&
               ( m1c .eq.  1 .and. igkst .gt. 4 ) ) then
                  call pdefs(io,ierr,ix,m1c,igkst,idsn(i))
                  goto 180
            end if
               nx = ix    + ns(m1c)
               n  = igkst - ns(m1c)
            if( ( m1c .ge. 30 .and. m1c. le. 39 .and. n .eq. 0 ) .or.&
               ( m1c .eq. 30 .and. igkst .eq. 9 ) .or.&
               ( m1c .eq. 39 .and.&
                 igkst .ge. 7 .and. igkst .le. 15 ) .or.&
               ( m1c .eq. 32 .and. igkst .eq. 1 ) ) then
               call mbody1(io,ierr,ix,igkst,m1c,idsn(i))
               goto 200
            end if
            if( n .eq. 0 .or.&
             ( m1c .ge. 16 .and.&
    &          m1c .le. 21 .and. n .eq. -1 ) ) then
               call sufchk(io,ierr,ix,nx,m0c,m1c,idsn(i))
               goto 180
            end if
               ierr = ierr + 1
               goto 180
  180    continue
            lsc(mxj+1) = ix + ns(m1c)
            kst(mxj)   = m1c
  200 continue
         if( ierr .ne. 0 ) goto 999
            call mbody2(io,ierr)
            if( ierr .ne. 0 ) goto 999
         do js = 1, mxj
            if( jtr( js ) .ne. 0 ) call trfsuf(io,ierr,js)
         end do
            if( ierr .ne. 0 ) goto 999
      if( junf .ne. 0 ) then
         do ic = 1, mxa
            if( ktr(ic) .ne. 0 ) call addsuf(io,ierr,ic)
         end do
            if( ierr .ne. 0 ) goto 999
      end if
      if( junf .ne. 0 ) then
         do 70 ic = 1, mxa
               if(mfl(1,ic).eq.0) goto 70
               if(mfl(1,ic).lt.0) goto 40
               if(mfl(3,ic).eq.0) mfl(3,ic) = ktr(ic)
               if(mfl(3,ic).eq.0) goto 70
               l = 0
            do m = 1, mxtr
               if( abs(trf(1,m)) .eq. mfl(3,ic) ) l = m
            end do
            if( l .eq. 0 ) then
               ierr = ierr + 1
               goto 999
            end if
               mfl(3,ic)=l
               goto 70
   40          lp=-mfl(1,ic)
            do 60 j=3,2+laf(1,lp+2)*laf(2,lp+2)*laf(3,lp+2)
               if(laf(1,lp+j).eq.0.or.&
                 laf(1,lp+j).eq.abs(jun(ic))) goto 60
               if(laf(3,lp+j).eq.0) laf(3,lp+j)=ktr(ic)
               if(laf(3,lp+j).eq.0) goto 60
               l=0
               do m=1,mxtr
                  if(abs(trf(1,m)).eq.laf(3,lp+j))l=m
               end do
               if(l.eq.0) then
                               ierr = ierr + 1
                  goto 999
               end if
               laf(3,lp+j)=l
   60       continue
   70      continue
      end if
            call chkcdf(io,ierr)
            if( ierr .ne. 0 ) goto 999
            if( nlat .ne. 0 ) call callat(io,ierr)
            if( ierr .ne. 0 ) goto 999
            if( igcnt.gt.1)then !FURUTA20190110
             if( nlat3 .gt. 0 ) call settetra(io,ierr)
            endif
            if( ierr .ne. 0 ) goto 999
            if( ierr .ne. 0 ) goto 999
            call celsuf(io,ierr)
            if( ierr .ne. 0 ) goto 999
            if( igcnt .eq. 1 ) then
                  close( io )
               if( icells .eq. 2 ) then
                  write(ioj) mxj, nljc
                  close(ioj)
                  stop
               else
                  open(io,form='formatted',status='scratch')
                  goto 2000
               end if
            end if
                  goto 998
  999 continue
            ierr = 1
  998 continue
            if( me .ne. 0 ) then
                  close( io )
            else
               if( iggcm .eq. 0 .and. ierr .eq. 0 ) then
                  close( io )
                  open(io,form='formatted',status='scratch')
               else
               end if
                  endfile( io )
            end if
            if( ifilt .ne. 0 )then
             close( ioc )
             if(nlatind.gt.0)call DEALLOCATE_ldata
            endif
      call moddas_deallocate_cha(chrg)
      return
      end

      subroutine pdefs(io,ierr,ix,m1c,nc,idsn)
      use moddas_ggs !frtati20220905
      implicit real*8 (a-h,o-z)
      include 'param.inc'
      include 'ggsparam.inc'
      dimension c(2,3),a(4,4)
      character hs(9)*24
      return
      end


      subroutine trfdef(io,jo,ierr,jt)
      use moddas_ggs !frtati20220905
      implicit real*8 (a-h,o-z)
      include 'param.inc'
      include 'ggsparam.inc'
      include 'err.inc'
      dimension tr(17),iy(2,3)
         call defbas(tr(3*j+2),tr(3*(mod(j,3)+1)+2),&
         tr(3*(mod(j+1,3)+1)+2),i)
         call crspro(tr(3*(mod(j,3)+1)+2),&
                    tr(3*(mod(j+1,3)+1)+2),tr(3*j+2))
         if(r.gt.2.e-6) then
               call ErrWrite(MsgID, ErrCha)
         end if
         call crspro(tr(5),tr(8),tpp)
         if(tr(14).ne.-1.) call matmpy(3,3,1,tr(5),3,tpp,3,tr(2),3,0)
         if(tr(14).eq.-1.) call matmpy(3,3,1,tr(5),3,tpp,3,tr(15),3,1)
      return
      end
      subroutine trfsuf(io,ierr,js)
      use moddas_ggs !frtati20220905
      implicit real*8 (a-h,o-z)
      include 'param.inc'
      include 'ggsparam.inc'
      dimension tm(4,4),aa(4,4),nc(23)
      data tm/1.,15*0./,nc/4,4*1,4,3*2,3*3,3*1,3*5,3*3,2*10/
      do 10 jt = 1, mxtr
   10    if(abs(trf(1,jt)).eq.jtr(ljtr+js)) goto 20
         ierr = ierr + 1
         return
   20    if(trf(14,jt).eq.-2.) return
         ks = kst(js)
         if(abs(ks-25).le.1) goto 180
      do 30 i = 1, 4
      do 30 j = 2, 4
   30    tm(j,i) = trf(3*i+j-3,jt)
         call sufdef(js,aa)
         tpp(18) = aa(2,2)
         tpp(19) = aa(3,3)
         tpp(20) = aa(4,4)
         call matmpy(4,4,4,aa,4,tm,4,tpp,4,0)
         call matmpy(4,4,4,tm,4,tpp,4,aa,4,4)
         ii = 3
         ic = 0
         im = 1
      do 40 i = 1, 3
         if(tpp(i+17).eq.0.) ic = ic+1
         if(abs(aa(i+1,i+1)).gt.abs(aa(ii+1,ii+1))) ii = i
   40    if(abs(aa(i+1,i+1)).lt.abs(aa(im+1,im+1))) im = i
         if(trf(1,jt).gt.0..and.ic.ne.3) ic = 0
         if(ic.gt.0) aa(im+1,im+1) = 0.
         if(ic.gt.1) aa(7-im-ii,7-im-ii) = 0.
         if(ic.gt.2) aa(ii+1,ii+1) = 0.
         tpp(17) = aa(1,1)
         tpp(10) = aa(1,1)
      do 50 i = 1, 3
         tpp(i+13) = 0.
         tpp(i+17) = 0.
         tpp(i) = aa(i+1,i+1)
         tpp(i+10) = aa(i+1,i+1)
         tpp(i+6) = 2.*aa(1,i+1)
         if(aa(i+1,i+1).ne.0.) tpp(i+17) = -aa(1,i+1)/aa(i+1,i+1)
         if(aa(i+1,i+1).eq.0.) tpp(i+13) = aa(1,i+1)
   50    tpp(17) = tpp(17)-aa(i+1,i+1)*tpp(i+17)**2
         if(ks.lt.10.or.ks.lt.23.and.trf(1,jt).lt.0.) goto 60
         if(ks.ge.16.and.ks.le.21.and.&
           scf(lsc(js)+nc(ks)).ne.0.) then
           
            ierr = ierr + 1
         end if
         ks = 23
         tpp(4) = 2.*aa(2,3)
         tpp(5) = 2.*aa(3,4)
         tpp(6) = 2.*aa(2,4)
         if(aa(2,3).ne.0..or.aa(3,4).ne.0..or.aa(2,4).ne.0.) goto 130
         ks = 22
         goto 110
   60    if(ks.ge.10) goto 80
         kk = 1
      do 70 i= 1 , 3
         if(ks.le.4) tpp(i) = 2.*tpp(i+13)
         if(ks.ge.5) tpp(i) = tpp(i+17)
   70    if(abs(tpp(i)).gt.abs(tpp(kk))) kk = i
         if(abs(tpp(1))+abs(tpp(2))+abs(tpp(3))-abs(tpp(kk)).gt.&
           epss*(abs(tpp(kk))+abs(tpp(17))).or.&
           ks.le.4.and.tpp(kk).lt.0.) kk = 0
         if(ks.ge.5) kk = kk+5
         ks = kk+1
         if(ks.gt.6)then
           if(abs(tpp(kk-5)).lt.-epss*tpp(17)) ks = 5
         endif
         if(abs(ks-3).le.1) tpp(17) = tpp(17)/tpp(kk)
         if(ks.gt.6) tpp(1) = tpp(kk-5)
         tpp(nc(ks)) = -tpp(17)
         goto 130
   80    if(ks.ge.16) goto 90
         tpp(1) = tpp(18)
         tpp(2) = tpp(20)
         if(im.eq.1) tpp(1) = tpp(19)
         if(im.eq.3) tpp(2) = tpp(19)
         ks = im+9
         if(abs(tpp(1))+abs(tpp(2)).lt.-epss*tpp(17)) ks = im + 12
         tpp(nc(ks)) = -tpp(17)
         goto 130
   90    if(ks.eq.22) goto 110
         k = lsc(js)+nc(ks)
         ic = 1
         if(tpp(1)*tpp(2).gt.0.) ic = 3
         if(tpp(1)*tpp(3).gt.0.) ic = 2
         j = mod(ks-16,3)+3*ic+2
         ks = ic+15
      do 100 i = 1, 3
  100    tpp(i) = tpp(17+i)
         if(abs(tpp(1))+abs(tpp(2))+abs(tpp(3))-abs(tpp(ic)).lt.&
           epss*(abs(tpp(ic))+1.)) ks = ic+18
         if(ks.gt.18) tpp(1) = tpp(ic)
         tpp(nc(ks)-1) = scf(k-1)
         tpp(nc(ks)) = scf(k)*sign(one,trf(j,jt))
         goto 130
  110 do 120 i = 1, 10
  120    tpp(i) = tpp(10+i)
  130    kst(js) = ks
         n = lsc(js)+nc(ks)-lsc(js+1)
         if(n.eq.0) goto 160
      do 140 iz = lsc(js+1)+1,lsc(mxj+1),1
         i = lsc(mxj+1)+lsc(js+1)+1-iz
         j = i
         if(n.lt.0) j = lsc(js+1)+1+lsc(mxj+1)-i
  140    scf(j+n) = scf(j)
      do 150 j = js, mxj
  150    lsc(j+1) = lsc(j+1)+n
  160 do 170 i = 1, nc(ks)
  170    scf(lsc(js)+i) = tpp(i)
         return
  180    continue
      return
      end
      subroutine sufdef(js,aa)
      use moddas_ggs !frtati20220905
      implicit real*8 (a-h,o-z)
      include 'param.inc'
      include 'ggsparam.inc'
      dimension aa(4,4)
      do 10 j = 1, 4
      do 10 i = 1, 4
   10    aa(i,j) = 0.
         l = lsc(js)+1
         k = kst(js)
      do 20 i = 2, 4
   20    if(k.ge.5.and.k.le.21) aa(i,i) = 1.

   30    aa(1,1) = -scf(l+3)
         aa(1,2) = .5*scf(l)
         aa(1,3) = .5*scf(l+1)
         aa(1,4) = .5*scf(l+2)
         goto 140
   40    aa(1,1) = -scf(l)
         aa(1,k) = .5
         goto 140
   50    aa(1,1) = -scf(l)
         goto 140
   60    aa(1,1) = scf(l)**2+scf(l+1)**2+scf(l+2)**2-scf(l+3)
         aa(1,2) = -scf(l)
         aa(1,3) = -scf(l+1)
         aa(1,4) = -scf(l+2)
         goto 140
   70    aa(1,1) = scf(l)**2-scf(l+1)
         aa(1,k-5) = -scf(l)
         goto 140
   80    aa(1,1) = scf(l)**2+scf(l+1)**2-scf(l+2)
         aa(1,2) = -scf(l)
         aa(1,3) = -scf(l+k/12)
         aa(1,4) = -scf(l+1)
         aa(1,k-8) = 0.
         aa(k-8,k-8) = 0.
         goto 140
   90    aa(1,1) = -scf(l)
         aa(k-11,k-11) = 0.
         goto 140
         aa(k-14,k-14) = -scf(l+3)
         goto 140
  110    aa(1,1) = -scf(l+1)*scf(l)**2
         aa(1,k-17) = scf(l)*scf(l+1)
         aa(k-17,k-17) = -scf(l+1)

  130    aa(1,1) = scf(l+9)
         aa(1,2) = .5*scf(l+6)
         aa(1,3) = .5*scf(l+7)
         aa(1,4) = .5*scf(l+8)
         aa(2,2) = scf(l)
         aa(2,3) = .5*scf(l+3)
         aa(2,4) = .5*scf(l+5)
         aa(3,3) = scf(l+1)
         aa(3,4) = .5*scf(l+4)
         aa(4,4) = scf(l+2)
  140 do 150 i = 2, 4
      do 150 j = 2, i
  150    aa(i,j-1) = aa(j-1,i)
      return
      end
      subroutine chkcdf(io,ierr)
      use GGBANKMOD  
      use LAFDATAMOD  
      use moddas_ggs  
          MsgID = 'L:2863/R:chkcdf/F:ggs00.f' !E06_001_001
          call ErrWrite(MsgID, ErrCha)
         if(kst(ja).le.21.or.kst(ja).ge.24) goto 830
      do 825 i = 1, na
         if(scf(lsc(ja)+i).eq.0.) goto 825
         if(scf(lsc(ja)+i)*scf(lsc(j)+i).lt.0.) is = -1
         goto 830
  825 continue
  830 do 835 i = 1, na
  835    if(abs(scf(lsc(ja)+i)*is-scf(lsc(j)+i)).gt.&
         1.d-12*abs(scf(lsc(ja)+i)*is&
         +scf(lsc(j)+i))) goto 870
         if( kst(ja).ge.24 .and. kst(ja).le.26&
            .and. jtr(ja).ne.jtr(j) ) goto 870 ! frtati 2022/04/08
  840    hb = ' '
         ka = nsf(ja)
         mq = kfq
         if(kfq.ne.0) write(hb,'(1h.,i1)') kfq
         hc = ' '
         k = nsf(j)
         if(kfq.ne.0) write(hc,'(1h.,i1)') kfq
         write(io,850) ka,hb,k,hc,k,hc

         if(idnt(ja).ne.0) goto 860
         in = in+1
         idnt(ja) = ja
         if(idne(1).eq.0) idne(2) = 2
         idne(1) = idne(1)+1
         idne(idne(2)+1) = 1
         idns(ja) = idne(2)+1
         idne(2) = idne(2)+2
         idne(idne(2)) = ja
         if(ksu(ja).ne.0) then
           
         end if
  860    idnt(j) = idnt(ja)*is
         jq = idns(ja)
         idns(j) = jq
         idne(jq) = idne(jq)+1
         idne(2) = idne(2)+1
         idne(idne(2)) = j
         if(ksu(j).ne.0) then
           
         end if
         mc = mc + 1
      ErrID = 'L:3194/R:chkcdf/F:ggs00.f' !E06_001_001
      call ErrWriteIO(ErrID,ErrCha,io)
  490    if(mfl(1,ic).eq.0) nlv(ic) = 1
  500 do 530 ic = 1, mxa
         if(nlv(ic).ne.0) goto 530
         l = 0
         lp = -mfl(1,ic)
         na = 1
         if(lp.gt.0) na = laf(1,lp+2)*laf(2,lp+2)*laf(3,lp+2) !FURUTA20201127
      do 520 j = 3, na+2
         lu = -lp
         if(lp.gt.0) lu = laf(1,lp+j) !FURUTA20201127
         if(lu.eq.0.or.lu.eq.abs(jun(ic))) goto 520
      do 510 i = 1, mxa
         if(abs(jun(i)).ne.lu) goto 510
         if(nlv(i).eq.0) goto 530
         l = max(l,nlv(i)+1)
  510 continue
         if(l.eq.0) then
           
            ierr = ierr + 1
            return
         end if
  520 continue
         nlv(ic) = max(l,1)
         nlev = max(nlev,l)
         goto 500
  530    continue
         if(nlev.gt.mxlv) then
           
            ierr = ierr + 1
            return
         end if
         goto 670
         nt = 0
         lv = 0
         iu(0) = 0
  540    il(lv) = 1
  550    if(jun(il(lv)).ne.iu(lv)) goto 640
      do 580 l = 0, lv-1
      do 570 k = abs(lca(il(lv))),abs(lca(il(lv)+1))-1
         jk = abs(lja(k))
         if(jk.gt.1000000) goto 570
         ks = 1
         if(kst(jk).le.4) ks = 2
      do 565 j = abs(lca(il(l))),abs(lca(il(l)+1))-1
         jj = abs(lja(j))
         if(jk.ne.jj) goto 565
      do 555 jt = 1, nt
  555    if(nint(tpp(jt)).eq.jj) goto 565
         ka = nsf(jj)
         hb = ' '
         if(kfq.ne.0) write(hb,'(1h.,i1)') kfq
      if( ks .eq. 1 ) then
        
         ierr = ierr + 1
         return
      else
        
      end if
        
         nt = nt + 1
         if(nt.gt.20) goto 670
         tpp(nt) = jj
  565 continue
  570 continue
  580 continue
         mu(lv) = mfl(1,il(lv))
         if(mu(lv).eq.0) goto 640
         iu(lv+1) = mu(lv)
         if(mu(lv).gt.0) goto 610
         jm(lv) = 3
         nu(lv) = laf(1,-mu(lv)+2)*laf(2,-mu(lv)+2)*laf(3,-mu(lv)+2) !FURUTA20201127
  590    iu(lv+1) = laf(1,-mu(lv)+jm(lv)) !FURUTA20201127
         if(iu(lv+1).eq.0.or.iu(lv+1).eq.iu(lv)) goto 630
      do 600 j = 3, jm(lv)-1
  600    if(iu(lv+1).eq.laf(1,-mu(lv)+j)) goto 630 !FURUTA20201127
  610    lv = lv+1
         goto 540
  620    lv = lv-1
         if(mu(lv).ge.0) goto 640
  630    jm(lv) = jm(lv)+1
         if(jm(lv).le.nu(lv)+2) goto 590
  640    il(lv) = il(lv)+1
         if(il(lv).le.mxa) goto 550
         if(lv.ne.0) goto 620
  670    nt= 4
         coincd = .0001
         coincd = 1.d-08
      do 690 js= 1, mxj
         if(kst(js).gt.4) goto 690
         if(kst(js).gt.1) ksc(js) = kst(js)
         if(ksc(js).gt.1) goto 690
         nt = nt + 1
         ksc(js) = nt
         is = lsc(js)+1
         ts = scf(is)**2+scf(is+1)**2+scf(is+2)**2
         if(scf(is+1)**2+scf(is+2)**2.lt.epss*ts) ksc(js) = 2
         if(scf(is)**2+scf(is+2)**2.lt.epss*ts) ksc(js) = 3
         if(scf(is)**2+scf(is+1)**2.lt.epss*ts) ksc(js) = 4
      do 680 jk = js+1,mxj
         if(kst(jk).gt.1) goto 680
         if(ksc(jk).ne.0) goto 680
         ik = lsc(jk)+1
         tk = scf(ik)**2+scf(ik+1)**2+scf(ik+2)**2
         if((scf(is)*scf(ik)+scf(is+1)*scf(ik+1)&
           +scf(is+2)*scf(ik+2))**2.gt..9998*ts*tk)&
           ksc(jk) = ksc(js)
  680 continue
  690 continue
      return
      end
      function namchg(mm,ji)
      use moddas_ggs !frtati20220905
      implicit real*8 (a-h,o-z)
      include 'param.inc'
      include 'ggsparam.inc'
      goto (10,30,50) mm
   10 do 20 namchg = 1, mxa
   20    if(ncl(namchg).eq.ji) return
         goto 70
   30 do 40 namchg = 1,mxj
   40    if(nsf(namchg).eq.abs(ji)) return
         goto 70
   50 do 60 kk = 1, nlja
         namchg = abs(lja(kk))
         if(namchg.gt.1000000) goto 60
         if(nsf(namchg).eq.abs(ji)) return
   60 continue
   70    namchg = 0
      return
      end
      subroutine defbas(a,b,c,ie)
      implicit double precision (a-h,o-z)
      dimension a(3),b(3),c(3)
         ie = 1
         j = 1
         if(abs(a(2)).gt.abs(a(1))) j = 2
         if(abs(a(3)).gt.abs(a(j))) j = 3
         k = mod(j,3)+1
         l = 6-j-k
         t = a(j)**2+a(k)**2
         if(t.eq.0.) return
         b(l) = .0637922
         b(k) = (a(j)*sqrt(t-(t+a(l)**2)*b(l)**2)-a(k)*a(l)*b(l))/t
         b(j) = -(b(k)*a(k)+b(l)*a(l))/a(j)
         call crspro(a,b,c)
         ie = 0
      return
      end
      subroutine crspro(v1,v2,v3)
      implicit double precision (a-h,o-z)
      dimension v1(3),v2(3),v3(3)
         v3(1) = v1(2)*v2(3)-v1(3)*v2(2)
         v3(2) = v1(3)*v2(1)-v1(1)*v2(3)
         v3(3) = v1(1)*v2(2)-v1(2)*v2(1)
      return
      end
      subroutine matmpy(m,n,l,a,ia,b,ib,c,ic,jt)
      implicit double precision (a-h,o-z)
      dimension a(ia,*),b(ib,*),c(ic,*)
         if(jt.ne.0) goto 20
      do 10 j = 1,l
      do 10 i = 1,m
         c(i,j) = 0.
      do 10 k = 1, n
   10    c(i,j) = c(i,j)+a(i,k)*b(k,j)
         return
   20 do 30 j = 1,l
      do 30 i = 1,m
         c(i,j) = 0.
      do 30 k = 1,n
   30    c(i,j) = c(i,j)+a(k,i)*b(k,j)
      return
      end
      function lgeval(lg,n)
      dimension lg(*)
         lg(n+1) = 1000002
         i = 0
         l = 0
   10    lgeval = 1
   20    i = i+1
         if(lg(i).gt.1) goto 40
         lgeval = min(lgeval,lg(i))
         if(lgeval.ne.0) goto 20
   30    i = i + 1
         if(lg(i).lt.2) goto 30
   40    if(lg(i).ne.1000003) goto 50
         if(lgeval.eq.0) goto 10
         if(l.eq.0) return
         goto 60
   50    if(lg(i).ne.1000001) goto 80
         l = l+1
         if(lgeval.ne.0) goto 20
   60    m = 1
   70    i = i+1
         if(lg(i).eq.1000001) m = m+1
         if(lg(i).eq.1000002) m = m-1
         if(m.ne.0) goto 70
   80    l = l-1
         if(i.ne.n+1) goto 20
      return
      end
      subroutine quart(n,c,r4,jj,mm)
      implicit double precision (a-h,o-z)
      dimension c(5),r4(4),tp(16)
      parameter (pf=(2.*3.1415926535898d0)/3.,th=1d0/3.)
      parameter ( epss = 1.0d-10 )
         jj = 0
         if(n.eq.3) goto 90
         tp(1) = .5*c(2)
         tp(2) = .25*c(2)
         tp(3) = c(3)-1.5*tp(1)**2
         tp(4) = c(4)+tp(1)*(tp(1)**2-c(3))
         tp(5) = c(5)-.0625*tp(1)*(5.*c(4)-c(3)*tp(1)+3.*tp(4))
         if(abs(tp(4)).le.epss*max(abs(c(4)),tp(1)**2,abs(tp(1)*c(3))))&
         goto 40
         tp(6) = 2.*th*tp(3)
         tp(7) = tp(3)**2-4.*tp(5)
         tp(8) = -tp(4)**2
         tp(9) = th*tp(7)-tp(6)**2
         tp(10) = .5*(tp(8)-tp(6)*(tp(7)-2.*tp(6)**2))
         tp(11) = tp(9)**3+tp(10)**2
         if(abs(tp(11)).le.epss*max(abs(tp(9)**3),tp(10)**2)) goto 10
         if(tp(11).gt.0.) goto 30
         tp(12) = sqrt(-tp(9))
         if(tp(6).ge.0..or.tp(7).le.0.) return
         tp(13) = -tp(10)/tp(12)**3
         if(abs(tp(13)).ge.1.) goto 10
         tp(13) = acos(tp(13))
         tp(14) = 2.*tp(12)*cos(th*tp(13))-tp(6)
         if(tp(14).lt.epss*abs(tp(6))) goto 40
         tp(15) = .5*sqrt(tp(14))
         tp(14) = -tp(14)-3.*tp(6)
         tp(16) = tp(4)/tp(15)
         if(abs(tp(16)).ge.tp(14)) goto 10
         jj = 4
         tp(13) = .5*sqrt(tp(14)+tp(16))
         r4(1) = -tp(13)-tp(15)-tp(2)
         r4(2) = tp(13)-tp(15)-tp(2)
         tp(13) = .5*sqrt(tp(14)-tp(16))
         r4(3) = tp(15)-tp(13)-tp(2)
         r4(4) = tp(15)+tp(13)-tp(2)
         return
   10    tp(12) = sign(sqrt(abs(tp(9))),tp(10))
         tp(13) = -2.*tp(12)-tp(6)
         tp(14) = tp(12)-tp(6)
         tp(15) = sign(.5*sqrt(abs(tp(13))),tp(4))
         if(tp(14).le.0.) goto 20
         tp(16) = sqrt(tp(14))
         jj = 2
         r4(1) = -tp(15)-tp(2)-tp(16)
         r4(2) = tp(16)-tp(15)-tp(2)
         if(mm.eq.0) return
         if(tp(10).eq.0.) return
         jj = 3
         r4(3) = tp(15)-tp(2)
         return
   20    if(mm.eq.0) return
         jj = 1
         r4(1) = -tp(15)-tp(2)
         return
   30    tp(12) = -tp(10)-sign(sqrt(tp(11)),tp(10))
         tp(12) = sign(abs(tp(12))**th,tp(12))
         tp(14) = tp(12)-tp(9)/tp(12)-tp(6)
         if(tp(14).lt.epss*abs(tp(6))) goto 40
         tp(15) = .5*sign(sqrt(tp(14)),tp(4))
         tp(16) = abs(tp(4)/tp(15))-tp(14)-3.*tp(6)
         if(tp(16).le.0.) return
         tp(16) = .5*sqrt(tp(16))
         jj = 2
         r4(1) = -tp(16)-tp(15)-tp(2)
         r4(2) = tp(16)-tp(15)-tp(2)
         return
   40    tp(6) = .5*tp(3)
         tp(7) = tp(6)**2-tp(5)
         tp(9) = 1.e-8*max(tp(6)**2,abs(tp(5)))
         if(tp(7).le.tp(9)) goto 60
         tp(8) = sqrt(tp(7))
         tp(16) = tp(8)-tp(6)
         if(tp(16).le.0.) goto 70
         tp(15) = -tp(6)-tp(8)
         tp(16) = sqrt(tp(16))
         if(tp(15).gt.0.) goto 50
         jj = 2
         r4(1) = -tp(16)-tp(2)
         r4(2) = tp(16)-tp(2)
         if(mm.eq.0) return
         if(tp(15).ne.0.) return
         jj = 3
         r4(3) = -tp(2)
         return
   50    tp(15) = sqrt(tp(15))
         jj = 4
         r4(1) = -tp(16)-tp(2)
         r4(2) = -tp(15)-tp(2)
         r4(3) = tp(15)-tp(2)
         r4(4) = tp(16)-tp(2)
         return
   60    if(mm.eq.0) return
         if(tp(7).lt.-tp(9).or.tp(6).gt.0.) return
         if(tp(6).eq.0.) goto 80
         jj = 2
         tp(16) = sqrt(-tp(6))
         r4(1) = tp(16)-tp(2)
         r4(2) = -tp(16)-tp(2)
         return
   70    if(mm.eq.0) return
         if(tp(16).lt.0.) return
   80    jj= 1
         r4(1) = -tp(2)
         return
   90    tp(1) = th*c(2)
         tp(2) = th*c(3)-tp(1)**2
         tp(3) = .5*(c(4)-tp(1)*(c(3)-2.*tp(1)**2))
         tp(4) = tp(2)**3+tp(3)**2
         if(tp(4).gt.0.) goto 110
         if(tp(4).eq.0.) goto 100
         tp(5) = 2.*sqrt(-tp(2))
         tp(6) = 8.*tp(3)/tp(5)**3
         if(abs(tp(6)).ge.1.) goto 100
         tp(7) = th*acos(-tp(6))
         r4(1) = tp(5)*cos(pf+tp(7))-tp(1)
         r4(2) = tp(5)*cos(pf-tp(7))-tp(1)
         r4(3) = tp(5)*cos(tp(7))-tp(1)
         jj = 3
         return
  100    r4(1) = -2.*sign(sqrt(-tp(2)),tp(3))-tp(1)
         jj = 1
         if(mm.eq.0) return
         if(r4(1).eq.0.) return
         jj = 2
         r4(2) = -.5*r4(1)
         return
  110    tp(5) = -tp(3)-sign(sqrt(tp(4)),tp(3))
         tp(5) = sign(abs(tp(5))**th,tp(5))
         r4(1) = tp(5)-tp(2)/tp(5)-tp(1)
         jj= 1
      return
      end
