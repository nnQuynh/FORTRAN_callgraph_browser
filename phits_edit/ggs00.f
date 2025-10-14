************************************************************************
*                                                                      *
      subroutine setgg(iom,jo,ierr)
*                                                                      *
*       setup GG geometry                                              *
*       Last modified by K.Niita on 2009/09/30                              *
*                                                                      *
************************************************************************

      use GGBANKMOD !FURUTA
      use LATDATAMOD
      use LAFDATAMOD !FURUTA20201127
      use moddas
      use moddas_character
      use moddas_ggs
      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'ggsparam.inc'
      include 'err.inc'

      parameter ( ibmt = 40 )

*-----------------------------------------------------------------------

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

*-----------------------------------------------------------------------
cFURUTA20150714 TETRA
      integer nlat3,itetcl(10),itetvol,itetauto,itgchk
      common /tetf1/ nlat3,itetcl,itetvol,itetauto,itgchk
      real(8) tetsfac(10)
      common /tetf2/ tetsfac,ltfile,itfform,tfilename
*-----------------------------------------------------------------------

      common /taliin/ rsouin, nzztin, nrgnin

      common /mpi00/ npe, me

*-----------------------------------------------------------------------

      dimension vtrs(13)

      logical deqn5
      dimension bval(50)

      dimension nr(39)
      dimension ns(39)

      data nr/4,4*1,4,3*2,3*3,3*1,3*5,3*3,2*10,3*7,3*0,24,15,4,18,18,
     &        10,18,20,24,32/

      data ns/4,4*1,4,3*2,3*3,3*1,3*5,3*3,2*10,3*6,3*0,12,6,4,7,
     &        12,7,8,12,24,15/

      character ksf(50)*3

      data ( ksf(i), i = 1, ibmt ) /
     &         'p  ','px ','py ','pz ','so ',
     &         's  ','sx ','sy ','sz ','c/x',
     &         'c/y','c/z','cx ','cy ','cz ',
     &         'k/x','k/y','k/z','kx ','ky ',
     &         'kz ','sq ','gq ','tx ','ty ',
     &         'tz ','x  ','y  ','z  ','box',
     &         'rpp','sph','rcc','rec','ell',
     &         'trc','wed','arb','rhp','hex'/

      dimension kfil(6)

*-----------------------------------------------------------------------


*-----------------------------------------------------------------------
      real*8 xangle, yangle, zangle
      integer inonzero
      real*8 rmatrix(3,3,3)

      common /geomemcom/ igeomem  ! T.Sato 2024/12/25
*-----------------------------------------------------------------------
      integer,allocatable :: lattmp(:),lattmp2(:)
*-----------------------------------------------------------------------
      common /ggcell/ icells, iobo
      character chcfg*100
      character chtrs*100

      common /paran/  icfn(100), ilfn(100), chfn(100)
      character       chfn*200

      logical   exex

*-----------------------------------------------------------------------

               ierr  = 0

*-----------------------------------------------------------------------

         if( icells .eq. 1 .or. icells .eq. 2 ) then

               chcfg = chfn(19)(1:ilfn(19))//'.cfg'
               ichcfg = ilfn(19) + 4

            if( icells .eq. 1 ) then

                  ierror = 0

                  inquire( file = chfn(19), exist = exex )

               if( exex .eqv. .false. ) then
                  write(*,'(/" Error : icells data file",
     &            " does not exist."/
     &            " file name = ",200a1)')
     &            ( chfn(19)(i:i),i=1, ilfn(19) )
                  ierror = ierror + 1
               end if

                  inquire( file = chcfg, exist = exex )

               if( exex .eqv. .false. ) then
                  write(*,'(/" Error : icells cfg data file",
     &            " does not exist."/
     &            " file name = ",200a1)')
     &            ( chcfg(i:i),i=1, ichcfg )
                  ierror = ierror + 1
               end if

                  if( ierror .gt. 0 ) goto 999

               iod = 26
               open(iod,status='unknown',file= chfn(19),
     &              form='unformatted')

            end if

               ioj = 203
               open(ioj,status='unknown',file= chcfg,
     &              form='unformatted')

         end if

*-----------------------------------------------------------------------
*     initialization of transform
*-----------------------------------------------------------------------

             if( icells .eq. 1 ) then

                read(ioj) igtrs, iobo

             else if( icells .eq. 2 ) then

                write(ioj) igtrs, iobo

             end if

*-----------------------------------------------------------------------

         if( igtrs .gt. 0 ) then

c allocate ggs array trf
               call moddas_allocate_dbl2( 17, 0, igtrs, trf )

*-----------------------------------------------------------------------

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

*-----------------------------------------------------------------------
*           re-read transform information and set valiables
*-----------------------------------------------------------------------

            if( igtrs .gt. 0 ) then

                  if( icells .ne. 1 ) then

                        rewind iob

                  else if( icells .eq. 1 ) then

                     if( iobo .eq. 1 ) then
                        iob = 71
                        open(iob,status='scratch',form='unformatted')
                     end if

                  end if

*-----------------------------------------------------------------------

               do i = 1, igtrs

*-----------------------------------------------------------------------

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

*-----------------------------------------------------------------------

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
                           write(iom,'(a,i8,a,/a,a)')
     &       '** warning : in definition of transform id', iabs(idtn(i))
     &      ,' in [transform] section,'
     &      ,'non-zero values of the 7th - 12th parameters are ignored,'
     &      ,' when M= -2 or 2.'
                           ErrCha = ''
                           MsgID = 'L:278/R:setgg/F:ggs00.f'
                           call ErrWrite(MsgID, ErrCha)
                           write(jo,'(a,i8,a,/a,a)')
     &       '** warning : in definition of transform id', iabs(idtn(i))
     &      ,' in [transform] section,'
     &      ,'non-zero values of the 7th - 12th parameters are ignored,'
     &      ,' when M= -2 or 2.'
                          end if

                         end if

                         vtrs(4) = dcos(zangle)*dcos(yangle)
                         vtrs(5) = dsin(zangle)*dcos(xangle)
     &                        +dcos(zangle)*dsin(yangle)*dsin(xangle)
                         vtrs(6) = dsin(zangle)*dsin(xangle)
     &                        -dcos(zangle)*dsin(yangle)*dcos(xangle)
                         vtrs(7) = -dsin(zangle)*dcos(yangle)
                         vtrs(8) = dcos(zangle)*dcos(xangle)
     &                        -dsin(zangle)*dsin(yangle)*dsin(xangle)
                         vtrs(9) = dcos(zangle)*dsin(xangle)
     &                        +dsin(zangle)*dsin(yangle)*dcos(xangle)
                         vtrs(10) = dsin(yangle)
                         vtrs(11) = -dcos(yangle)*dsin(xangle)
                         vtrs(12) = dcos(yangle)*dcos(xangle)

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
       write(ErrCha,'("For M=+-3,",
     & " transform axis should be either 1(=x), 2(=y), or 3(=z)",
     & " in ",i2,"th data of transform id",i8)') 5+im*2,iabs(idtn(i))
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
                              rmatrix(j,k,3)=rmatrix(j,k,3)
     &                             +rmatrix(j,l,2)*rmatrix(l,k,1)
                             enddo
                            enddo
                           enddo
                          endif
                         enddo
                         if(vtrs(4).ne.0
     &                        .or.vtrs(5).ne.0
     &                        .or.vtrs(6).ne.0) then ! rotation origin is defined
                          if(idnint(vtrs(13)).eq.3) then ! X0=X0+XR-(R3R2R1)XR
                           vtrs(1)=vtrs(1)+vtrs(4)
     &                          -rmatrix(1,1,3)*vtrs(4)
     &                          -rmatrix(1,2,3)*vtrs(5)
     &                          -rmatrix(1,3,3)*vtrs(6)
                           vtrs(2)=vtrs(2)+vtrs(5)
     &                          -rmatrix(2,1,3)*vtrs(4)
     &                          -rmatrix(2,2,3)*vtrs(5)
     &                          -rmatrix(2,3,3)*vtrs(6)
                           vtrs(3)=vtrs(3)+vtrs(6)
     &                          -rmatrix(3,1,3)*vtrs(4)
     &                          -rmatrix(3,2,3)*vtrs(5)
     &                          -rmatrix(3,3,3)*vtrs(6)
                          else ! X0=X0+XR-(R3R2R1)^T XR
                           vtrs(1)=vtrs(1)+vtrs(4)
     &                          -rmatrix(1,1,3)*vtrs(4)
     &                          -rmatrix(2,1,3)*vtrs(5)
     &                          -rmatrix(3,1,3)*vtrs(6)
                           vtrs(2)=vtrs(2)+vtrs(5)
     &                          -rmatrix(1,2,3)*vtrs(4)
     &                          -rmatrix(2,2,3)*vtrs(5)
     &                          -rmatrix(3,2,3)*vtrs(6)
                           vtrs(3)=vtrs(3)+vtrs(6)
     &                          -rmatrix(1,3,3)*vtrs(4)
     &                          -rmatrix(2,3,3)*vtrs(5)
     &                          -rmatrix(3,3,3)*vtrs(6)
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

                        if( itrs .ne. 0 )
     &                  trf( 1 + j, i ) = cos( vtrs(j) * pi / 180. )

                      end if

                        if( abs(  trf(1+j,i)
     &                    - anint(trf(1+j,i)) ) .gt. 1e-10 )
     &                  trf( ltrf + 1, i ) = idtn(i)

                     end if

                     if ( j .eq. 13 )
     &                    trf( 1 + j, i ) = dsign(1d0,vtrs(j))

                  end do

                     imxt = i

                     call trfdef(iom,jo,ierr,imxt)

               end do

*-----------------------------------------------------------------------

                     if( ierr .ne. 0 ) then
                        iog = 20
                        io  = iog
                        open(io,form='formatted',status='scratch')
                        goto 999
                     end if

*-----------------------------------------------------------------------

            end if

*-----------------------------------------------------------------------

         end if

*-----------------------------------------------------------------------
*     start of GG initialization
*-----------------------------------------------------------------------

            if( icgg  .eq. 0 ) return

*-----------------------------------------------------------------------
*        open temporary file for errors
*-----------------------------------------------------------------------

               iog = 20
               io  = iog

               open(io,form='formatted',status='scratch')

*-----------------------------------------------------------------------
*        set initial pointer
*-----------------------------------------------------------------------
               mcmx = ( mdas / 2  - 1 ) * 8 + 1
               call moddas_allocate_cha(MAX_NUM_CHRG, chrg)
               mci = 0

               ngstar = mmmax

*-----------------------------------------------------------------------
*        read data for icells = 1
*-----------------------------------------------------------------------

            if( icells .eq. 1 ) then

                  read(ioj)  igcel, igsuf, ifilt, nlat3

               do i = 1, igcel

                  read(ioj)  ichl(i), chsm(i), ichp(i), ilat(i),
     &                       idct(i), ifil(i), ilike(i), idmg(i),
     &                       iuni(i), deng(i), dvol(i), dtmp(i),
     &                       idrg(i), idgr(idrg(i))

               end do

               do i = 1, nlat3

                  read(ioj)  itetcl(i), tetsfac(i)

               end do

               if( ifilt .gt. 0 ) then

                  inquire( file = chfn(18), exist = exex )

                  if( exex .eqv. .false. ) then
                     write(*,'(/" Error : voxel data file",
     &               " does not exist."/
     &               " file name = ",200a1)')
     &               ( chfn(18)(i:i),i=1, ilfn(18) )
                     goto 999
                  end if

                  ioc = 66
                  open(ioc,status='unknown',file= chfn(18),
     &                 form='unformatted')

               end if

*-----------------------------------------------------------------------
*        write data for icells = 2
*-----------------------------------------------------------------------

            else if( icells .eq. 2 ) then

                  write(ioj) igcel, igsuf, ifilt, nlat3

               do i = 1, igcel

                  write(ioj) ichl(i), chsm(i), ichp(i), ilat(i),
     &                       idct(i), ifil(i), ilike(i), idmg(i),
     &                       iuni(i), deng(i), dvol(i), dtmp(i),
     &                       idrg(i), idgr(idrg(i))

               end do

               do i = 1, nlat3

                  write(ioj) itetcl(i), tetsfac(i)

               end do

            end if

*-----------------------------------------------------------------------
*        set number of cell and transforms
*-----------------------------------------------------------------------

               mxa    = igcel
               mxafs  = mxa
               mxtr   = igtrs

               nrgnin = mxa

*-----------------------------------------------------------------------
*        skip for icells = 1
*-----------------------------------------------------------------------

            if( icells .eq. 1 ) goto 10000

*-----------------------------------------------------------------------
*        set temporary pointer for normal case
*-----------------------------------------------------------------------

c allocate kst temporary
               call moddas_allocate_int(igsuf, kst)

*-----------------------------------------------------------------------
*        re-read surface information and get dimension
*-----------------------------------------------------------------------

               nsc = 0
               mxj = 0

               rewind  ioa

         do i = 1, igsuf

                  kst(i) = 0

            read(ioa) idrf, idtr, idsf, igkst, ( bval(j), j = 1, igkst )

               mxj = mxj + 1

               n = 10

               if( idtr .eq. 0  .and.
     &             idsf .le. 26 .and. idsf .ge. 1 ) n = nr(idsf)

               if( idtr .ne. 0  .and.
     &             idsf .le. 4 .and. idsf .ge. 2 ) n = 11 ! 2016/8/16 Avoid error in transform of PX, PY, or PZ

               if( idsf .eq. 1 .and. igkst .gt. 8 ) n = 10
               if( idsf .ge. 30 ) n = nr(idsf)

               nsc = nsc + n

            if( idsf .ge. 30 .and. idsf .le. 39 ) then

               jn = 6
               if( idsf .eq. 33 .or.
     &             idsf .eq. 34 .or.
     &             idsf .eq. 36 ) jn = 3
               if( idsf .eq. 39 ) jn = 8
               if( idsf .eq. 32 .or.
     &             idsf .eq. 35 ) jn = 0
               if( idsf .eq. 37 ) jn = 5

                  mxj = mxj + jn

                  kst(i) = jn

            end if

         end do

*-----------------------------------------------------------------------
*        re-read cell information and get dimension
*-----------------------------------------------------------------------

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
                      write(ErrCha,'("More than two cells include #all."
     &                ," 2nd cell with #all is ",i8)') idrg(i)
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

               if( img .eq. 1 .and.
     &           ( ( k .eq. 5 .and. ic .eq. ichl(i) ) .or.
     &             ( k .ne. 5 ) ) ) then

                     img = 0
                     m1c = m1c + 1

*-----------------------------------------------------------------------

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

*-----------------------------------------------------------------------
                           m1c = m1c + 2 * kst(j)
*-----------------------------------------------------------------------

                        if( idct(i) .ne. 0 ) then

*-----------------------------------------------------------------------
                           mxj = mxj + kst(j) + 1
                           nsc = nsc + max( 10, 5 * kst(j) ) ! 4*kst is occasionally insufficient
*-----------------------------------------------------------------------
                        end if

                     end if

                    end do

                 end if

                 k_back = k

! Nais_2024 <<<
*-----------------------------------------------------------------------

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
             write(ErrCha,'("Warning: #all and #cell coexist in cell",
     &       " no =",i8)') idrg(i)
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

*-----------------------------------------------------------------------
*        check repeated structure
*-----------------------------------------------------------------------

               nlat = 0

            do i = 1, mxa

               if( iuni(i) .ne. 0 .or. ilat(i) .ne. 0 .or.
     &             idct(i) .ne. 0 ) junf = 1

               if( ilat(i) .ne. 0 ) nlat = nlat + 1

            end do

               if( ifilt .ne. 0 ) junf = 1

*-----------------------------------------------------------------------
*        set up array for repeated structures / lattice universe map.
*-----------------------------------------------------------------------

         if( junf .ne. 0 ) then

c allocate jun, mfl, mazp
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

*-----------------------------------------------------------------------
*        re-read filling cell information
*-----------------------------------------------------------------------

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
cFURUTA20181009---------------------------------------------------------
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
                     call latdecode(latot2,lattmp2,latot,lattmp,
     &                    icount,ierr)
                     do j=1,latot
                      laf(1,llaf+mlaf+2+j)=lattmp(j) !FURUTA20201127
                     enddo
                     deallocate(lattmp,lattmp2)
                     if(ierr.ne.0)then
                      write(io,'(/"** ERROR : in [cell] section, ",
     &                   "compressed lattice element is inconsistent."/
     &                     "   latot =",i9,/
     &                     "   recoded number=",i9,/)')  latot,icount
                      goto 999
                     endif
*-----------------------------------------------------------------------
                    endif              !FURUTA20181009

                     do k = 1, 6

                        if( mod(k,2) .ne. 0 ) then

                           laf(1+k/2,llaf+mlaf+1) = kfil(k) !FURUTA20201127

                        else if( mod(k,2) .eq. 0 ) then

                           laf(k/2,llaf+mlaf+2) = kfil(k)       !FURUTA20201127
     &                                          - kfil(k-1) + 1 !FURUTA20201127

                        end if

                     end do

                     mlaf = mlaf + ( latot + 2 ) !FURUTA20201127

                  end if

               end do


                  if( ierr .ne. 0 ) goto 999

            end if

*-----------------------------------------------------------------------
*        write for icells = 2
*-----------------------------------------------------------------------

            if( icells .eq. 2 ) then

                  rewind  ioa

                  write(ioj) mxj, nsc, ncomp, mlja, mxit,
     &                       nlat, junf, mlaf, nmzu, nmaz

               do i = 1, igsuf

                  read(ioa)  idrf, idtr, idsf, igkst,
     &                       ( bval(j), j = 1, igkst )

                  write(ioj) idrf, idtr, idsf, igkst,
     &                       ( bval(j), j = 1, igkst )

               end do

            end if

*-----------------------------------------------------------------------
*        skip for icells = 1
*-----------------------------------------------------------------------

10000       continue

            if( icells .eq. 1 ) then

                  read(ioj)  mxj, nsc, ncomp, mlja, mxit,
     &                       nlat, junf, mlaf, nmzu, nmaz

                  rewind  ioa

               do i = 1, igsuf

                  read(ioj)  idrf, idtr, idsf, igkst,
     &                       ( bval(j), j = 1, igkst )

                  read(ioa)  idrf, idtr, idsf, igkst,
     &                       ( bval(j), j = 1, igkst )

               end do

            end if

*-----------------------------------------------------------------------
*     set dimension pointers
*-----------------------------------------------------------------------

            mlja0 = mlja + 2 * mxit * ncomp
            mxj0  = mxj

*-----------------------------------------------------------------------

               igcnt = 0

*-----------------------------------------------------------------------

            if( icells .eq. 1 ) then

               igcnt = 1

               read(ioj)  mxj, nljc

            end if

*-----------------------------------------------------------------------

 2000 continue

            igcnt = igcnt + 1

         if( igcnt .eq. 1 ) then

            mlja = mlja + 2 * mxit * ncomp

         else

            mlja  = nljc + 1

         end if

*-----------------------------------------------------------------------

            mlaj = ( 12 + 50 * junf ) * mxa + 50

            mtasks = 1
            nlse = 0

*-----------------------------------------------------------------------
*        double variables
*-----------------------------------------------------------------------

            call moddas_allocate_dbl(nsc+igeomem, scf) ! T.Sato 2024/12/25 change +2000 to +igeomem
            call moddas_allocate_dbl3(3, 7, nlat, vcl)

*-----------------------------------------------------------------------
*        integer variables
*-----------------------------------------------------------------------

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

*-----------------------------------------------------------------------


            nlaj_bank  = ( mlaj + mxa ) * mtasks !FURUTA
            nlcaj_bank = ( mlja + 1 ) * mtasks   !FURUTA

            kdb  = 0

*-----------------------------------------------------------------------

            mnax = mmmax

            mmmax  = mnax + 1
            ngfini = mnax + 1

            if( igcnt .eq. 1 ) ngfin0 = mnax


*-----------------------------------------------------------------------
*        initialize of dimensions
*-----------------------------------------------------------------------


            call INIT_laf !FURUTA20201127

            lsc(1) = 0

*-----------------------------------------------------------------------
*        re-read cell information and set valiables
*-----------------------------------------------------------------------

            nlja = 0
            nljc = 0

            rewind iod

         do i = 1, mxa

               lca( i ) = nlja + 1
               ncl( i ) = idrg(i)

            read(iod) (chrg(k:k),k=1,ichp(i))

*-----------------------------------------------------------------------

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

               if( img .eq. 1 .and.
     &           ( ( k .eq. 5 .and. ic .eq. ichl(i) ) .or.
     &             ( k .ne. 5 ) ) ) then

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

*-----------------------------------------------------------------------
*                 macrobody ( cvvv is real )
*-----------------------------------------------------------------------

                  if( abs( cvvv - anint(cvvv) ) .gt.
     &                5e-14 * abs(cvvv) ) then

                     lja( nlja ) =
     &                            nint( sign( aint(abs(cvvv)), cvvv) )
                     idna( nlja ) =
     &                    - nint( 10. * ( abs(cvvv)-aint(abs(cvvv)) ) )

                  end if

*-----------------------------------------------------------------------

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

*-----------------------------------------------------------------------
*        set up array for repeated structures / lattice universe map.
*-----------------------------------------------------------------------

         if( junf .ne. 0 ) then

            do i = 1, mxa

               jun(i) = iuni(i)

            end do

*-----------------------------------------------------------------------
*        re-read filling cell information
*-----------------------------------------------------------------------

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
cFURUTA20181009---------------------------------------------------------
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
                     call latdecode(latot2,lattmp2,latot,lattmp,
     &                    icount,ierr)
                     do j=1,latot
                      laf(1,llaf+mlaf+2+j)=lattmp(j) !FURUTA20201127
                     enddo
                     do j=1,latot
                      laf(3,llaf+mlaf+2+j)=0 !FURUTA20201127
                     enddo
                     deallocate(lattmp,lattmp2)
                     if(ierr.ne.0)then
                      write(io,'(/"** ERROR : in [cell] section, ",
     &                   "compressed lattice element is inconsistent."/
     &                     "   latot =",i9,/
     &                     "   recoded number=",i9,/)')  latot,icount
                      goto 999
                     endif
*-----------------------------------------------------------------------
                    endif              !FURUTA20181009

                     do k = 1, 6

                        if( mod(k,2) .ne. 0 ) then

                           laf(1+k/2,llaf+mlaf+1) = kfil(k) !FURUTA20201127

                        else if( mod(k,2) .eq. 0 ) then

                           laf(k/2,llaf+mlaf+2) = kfil(k)       !FURUTA20201127
     &                                          - kfil(k-1) + 1 !FURUTA20201127

                        end if

                     end do

                        mlaf = mlaf + ( latot + 2 ) !FURUTA20201127

                  end if

*-----------------------------------------------------------------------
*              set up lattice type
*-----------------------------------------------------------------------

                  if( ilat(i) .ne. 0 ) then

                     lat(1,i) = ilat(i)

                     nlat = nlat + 1

                     lat(2,i) = nlat

                  end if

*-----------------------------------------------------------------------
*              set up cell transformation
*-----------------------------------------------------------------------

                  if( idct(i) .ne. 0 ) then

                     ktr( i ) = idct(i)

                  end if

*-----------------------------------------------------------------------

               end do

         end if

*-----------------------------------------------------------------------
*        re-read surface information and set valiables
*-----------------------------------------------------------------------

               mxj = 0

            rewind  ioa

      do 200 i = 1, igsuf

            read(ioa) idrf, idtr, idsf, igkst, ( bval(j), j = 1, igkst )

! T.Sato 2025/03/11, avoid 0 for TRC
      if(idsf.eq.36.and.bval(8).eq.0.0) then
       write(ErrCha,'("Top radius of TRC surface (R2) should not be 0.",
     & " It is automatically changed to 1.0e-5")')
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

*-----------------------------------------------------------------------

               ix = lsc( mxj )

*-----------------------------------------------------------------------
*        translate point-defined surfaces into standard surfaces.
*-----------------------------------------------------------------------

            if( ( m1c .ge. 27 .and. m1c .le. 29 ) .or.
     &          ( m1c .eq.  1 .and. igkst .gt. 4 ) ) then

                  call pdefs(io,ierr,ix,m1c,igkst,idsn(i))

                  goto 180

            end if

*-----------------------------------------------------------------------
*        check of surface items
*-----------------------------------------------------------------------

               nx = ix    + ns(m1c)
               n  = igkst - ns(m1c)

*-----------------------------------------------------------------------
*        set up surface for macrobody
*-----------------------------------------------------------------------

            if( ( m1c .ge. 30 .and. m1c. le. 39 .and. n .eq. 0 ) .or.
     &          ( m1c .eq. 30 .and. igkst .eq. 9 ) .or.
     &          ( m1c .eq. 39 .and.
     &            igkst .ge. 7 .and. igkst .le. 15 ) .or.
     &          ( m1c .eq. 32 .and. igkst .eq. 1 ) ) then

               call mbody1(io,ierr,ix,igkst,m1c,idsn(i))

               goto 200

            end if

*-----------------------------------------------------------------------
*        set up surfaces and check them
*-----------------------------------------------------------------------

            if( n .eq. 0 .or.
     &        ( m1c .ge. 16 .and.
     &          m1c .le. 21 .and. n .eq. -1 ) ) then

               call sufchk(io,ierr,ix,nx,m0c,m1c,idsn(i))

               goto 180

            end if

*-----------------------------------------------------------------------
*           input error
*-----------------------------------------------------------------------

               write(io,'(/"** ERROR : in [surface] section, ",
     &                     "number of entries is incorrect."/
     &                     "   surface id =",i8)')  idsn(i)

               ierr = ierr + 1

               goto 180

*-----------------------------------------------------------------------
*     increment coefficients count and set up the final surface type.
*-----------------------------------------------------------------------

  180    continue

            lsc(mxj+1) = ix + ns(m1c)
            kst(mxj)   = m1c

*-----------------------------------------------------------------------

  200 continue

         if( ierr .ne. 0 ) goto 999

*-----------------------------------------------------------------------
*     expand macrobody
*-----------------------------------------------------------------------

            call mbody2(io,ierr)

            if( ierr .ne. 0 ) goto 999

*-----------------------------------------------------------------------
*     transform the coefficients of any surfaces that need it.
*-----------------------------------------------------------------------

         do js = 1, mxj

            if( jtr( js ) .ne. 0 ) call trfsuf(io,ierr,js)

         end do

            if( ierr .ne. 0 ) goto 999

*-----------------------------------------------------------------------
*     create surfaces for any cells that have trcl entries
*-----------------------------------------------------------------------

      if( junf .ne. 0 ) then

         do ic = 1, mxa

            if( ktr(ic) .ne. 0 ) call addsuf(io,ierr,ic)

         end do

            if( ierr .ne. 0 ) goto 999

      end if

*-----------------------------------------------------------------------
*     endow each deprived universe with its cell's transformation.
*     convert mfl(3, ) entries from names to indexes.
*-----------------------------------------------------------------------

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

               write(io,'(/"** ERROR : in setup GG"/
     &                     "   transform = ",i6/
     &                     " for fill of cell = ",i7,
     &                     " is absent."/)')
     &                     mfl(3,ic), ncl(ic)

               ierr = ierr + 1
               goto 999

            end if

               mfl(3,ic)=l
               goto 70

   40          lp=-mfl(1,ic)
            do 60 j=3,2+laf(1,lp+2)*laf(2,lp+2)*laf(3,lp+2)
               if(laf(1,lp+j).eq.0.or.
     &            laf(1,lp+j).eq.abs(jun(ic))) goto 60
               if(laf(3,lp+j).eq.0) laf(3,lp+j)=ktr(ic)
               if(laf(3,lp+j).eq.0) goto 60
               l=0

               do m=1,mxtr
                  if(abs(trf(1,m)).eq.laf(3,lp+j))l=m

               end do

               if(l.eq.0) then

                  write(io,'(/"** ERROR : in setup GG"/
     &                        "   transform = ",i6/
     &                        " for fill of cell = ",i7,
     &                        " is absent."/)')
     &                        laf(3,lp+j), ncl(ic)

                  ierr = ierr + 1
                  goto 999

               end if

               laf(3,lp+j)=l

   60       continue
c-----------------------------------------------------------
   70      continue

      end if

*-----------------------------------------------------------------------
*     expand # operators and check the cells and surfaces.
*-----------------------------------------------------------------------

            call chkcdf(io,ierr)

            if( ierr .ne. 0 ) goto 999

*-----------------------------------------------------------------------
*     calculate the constants of any lattices in the geometry.
*-----------------------------------------------------------------------

            if( nlat .ne. 0 ) call callat(io,ierr)

            if( ierr .ne. 0 ) goto 999

*-----------------------------------------------------------------------
*     set up geometry of tetrahedrons. (FURUTA20150714 TETRA)
*-----------------------------------------------------------------------

            if( igcnt.gt.1)then !FURUTA20190110
             if( nlat3 .gt. 0 ) call settetra(io,ierr)
            endif

            if( ierr .ne. 0 ) goto 999

*-----------------------------------------------------------------------
*     set up array for repeated structures / lattice universe map.
*-----------------------------------------------------------------------


            if( ierr .ne. 0 ) goto 999

*-----------------------------------------------------------------------
*     finish up the cells and surfaces.
*-----------------------------------------------------------------------

            call celsuf(io,ierr)

            if( ierr .ne. 0 ) goto 999

*-----------------------------------------------------------------------
*     set again for saving the memory
*-----------------------------------------------------------------------

            if( igcnt .eq. 1 ) then

                  close( io )

               if( icells .eq. 2 ) then

                  write(ioj) mxj, nljc

                  write(*,*) 'Binary file was successfully generated'
     &            ,' by icells = 2'

                  close(ioj)

                  stop

               else

                  open(io,form='formatted',status='scratch')
                  goto 2000

               end if

            end if

                  goto 998

*-----------------------------------------------------------------------

  999 continue

            ierr = 1

*-----------------------------------------------------------------------

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

*-----------------------------------------------------------------------
*     close temporary files
*-----------------------------------------------------------------------

            if( ifilt .ne. 0 )then
             close( ioc )
             if(nlatind.gt.0)call DEALLOCATE_ldata
            endif

*-----------------------------------------------------------------------

      call moddas_deallocate_cha(chrg)
      return
      end

************************************************************************
*                                                                      *
      subroutine pdefs(io,ierr,ix,m1c,nc,idsn)
*                                                                      *
*       translate point-defined surface mxj into a standard surface.   *
*       Last modified by K.Niita on 2009/09/30                         *
*                                                                      *
************************************************************************
      use moddas_ggs !frtati20220905

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'ggsparam.inc'

*-----------------------------------------------------------------------

      dimension c(2,3),a(4,4)

      character hs(9)*24

      data hs/'plane','cylinder','one-sheet cone','sphere','paraboloid',
     1 'hyperboloid of one sheet','hyperboloid of 2 sheets','ellipsoid',
     2 'two-sheet cone'/

*-----------------------------------------------------------------------

      if( m1c .eq. 1 ) goto 170

*-----------------------------------------------------------------------
*           check for an acceptable number of input data items.
*-----------------------------------------------------------------------

               if( nc .ne. 2 .and. nc .ne. 4 .and. nc .ne. 6 ) then

                  write(io,'(/"** ERROR : in [surface] section,"/
     &            "   number of entries is incorrect."/
     &            "   surface id =",i7)')  idsn

                  ierr = ierr + 1

                  return

               end if

*-----------------------------------------------------------------------
*        make it three points if only one or two are given.
*-----------------------------------------------------------------------

            do j = 1, 3
            do i = 1, 2

               c(i,j) = scf(ix+i+2*(j-1))

            end do
            end do

            do i = ( nc + 2 ) / 2, 3
            do j = 1, 2

               c(j,i) = c(j,i-1)

            end do
            end do

            do i = 1, nc

               scf(ix+i) = 0.

            end do

*-----------------------------------------------------------------------
*        set up to detect special cases.
*-----------------------------------------------------------------------

         t = sqrt(max((c(1,1)-c(1,2))**2+(c(2,1)-c(2,2))**2,
     &           (c(1,2)-c(1,3))**2+(c(2,2)-c(2,3))**2,
     &           (c(1,3)-c(1,1))**2+(c(2,3)-c(2,1))**2))
         t1 = min(c(1,1),c(1,2),c(1,3))
         t2 = max(c(1,1),c(1,2),c(1,3))

*-----------------------------------------------------------------------
*        plane if it is narrow or cylinder if it is horizontal.
*-----------------------------------------------------------------------

         ks = 3
         if(max(abs(c(2,1)-c(2,2)),abs(c(2,2)-c(2,3)),
     &          abs(c(2,3)-c(2,1))).le.1.e-12*t) ks = 2
         if(t2-t1.le.1.e-12*t) ks = 1
         if(ks.eq.3) goto 50
         m1c = m1c-36+11*ks
         scf(ix+1) = ((c(ks,1)+c(ks,2)+c(ks,3))/3.)**ks
         goto 160

*-----------------------------------------------------------------------
*        one-sheet cone if it is a straight line.
*-----------------------------------------------------------------------

   50    if(abs((c(1,2)-c(1,1))*(c(2,3)-c(2,2))-(c(2,2)-c(2,1))*
     &    (c(1,3)-c(1,2)))*1.e12.gt.t**2) goto 60
         m1c = m1c-8
         t1 = c(1,1)+c(1,2)+c(1,3)
         t2 = c(2,1)+c(2,2)+c(2,3)
         scf(ix+3) = (3.*(c(1,1)*c(2,1)+c(1,2)*c(2,2)+c(1,3)*c(2,3))-
     &    t1*t2)/(3.*(c(1,1)**2+c(1,2)**2+c(1,3)**2)-t1**2)
         scf(ix+2) = scf(ix+3)**2
         scf(ix+1) = (t1-t2/scf(ix+3))/3.
         goto 160

*-----------------------------------------------------------------------
*        check for a plane of two sheets.
*-----------------------------------------------------------------------

   60    if(abs(c(2,1)-c(2,2)).le.1.e12*abs(c(1,1)-c(1,2)).and.
     1    abs(c(2,2)-c(2,3)).le.1.e12*abs(c(1,2)-c(1,3)).and.
     2    abs(c(2,3)-c(2,1)).le.1.e12*abs(c(1,3)-c(1,1))) goto 70

                  write(io,'(/"** ERROR : in [surface] section,"/
     &            "   surface would be two parallel planes."/
     &            "   surface id =",i7)')  idsn

                  ierr = ierr + 1

         m1c = m1c-25
         ks = 0
         goto 160

*-----------------------------------------------------------------------
*        find coefficients of symmetric gq equation by elimination.
*        r**2 + a(4,1)*x**2 + a(4,2)*x + a(4,3) = 0.
*-----------------------------------------------------------------------

   70    i1 = 1

      do 80 i = 1, 3
         a(i,1) = c(1,i)**2
         if(a(i,1).gt.a(i1,1)) i1 = i
         a(i,2) = c(1,i)
         a(i,3) = 1.
   80    a(i,4) = -c(2,i)**2

         i2 = max(1,3-i1)

      do 100 i = 1, 3
         if(i.eq.i1) goto 100
      do 90 j = 2, 4
   90    a(i,j) = a(i,j)-a(i,1)*a(i1,j)/a(i1,1)
         if(abs(a(i,2)).gt.abs(a(i2,2))) i2 = i
  100 continue

         i3 = 6-i1-i2
         a(i3,3) = a(i3,3)-a(i3,2)*a(i2,3)/a(i2,2)
         a(4,3) = (a(i3,4)-a(i3,2)*a(i2,4)/a(i2,2))/a(i3,3)
         a(4,2) = (a(i2,4)-a(i2,3)*a(4,3))/a(i2,2)
         a(4,1) = (a(i1,4)-a(i1,3)*a(4,3)-a(i1,2)*a(4,2))/a(i1,1)

*-----------------------------------------------------------------------
*        set up coefficients for simplest applicable standard surface.
*        sphere if x**2 coefficient is 1.
*-----------------------------------------------------------------------

         if(abs(a(4,1)-1.)*1.e12.gt.1.) goto 120
         ks = 4

*-----------------------------------------------------------------------
*        sphere at origin if linear term is zero.
*-----------------------------------------------------------------------

         if(abs(a(4,2))*1.e12.gt.t) goto 110
         m1c = 5
         scf(ix+1) = -a(4,3)
         goto 160

*-----------------------------------------------------------------------
*        sphere, not at origin.
*-----------------------------------------------------------------------

  110    m1c = m1c-20
         scf(ix+1) = -.5*a(4,2)
         scf(ix+2) = scf(ix+1)**2-a(4,3)
         goto 160

*-----------------------------------------------------------------------
*        sq surface of some sort.
*-----------------------------------------------------------------------

  120 do 130 i = 1, 10
  130    scf(ix+11-i) = i/8

*-----------------------------------------------------------------------
*        paraboloid if it has no x**2 term.  let constant term be zero.
*-----------------------------------------------------------------------

         if(abs(a(4,1))*1.e12.gt.1.) goto 140
         scf(ix+m1c-26) = 0.
         scf(ix+m1c-23) = .5*a(4,2)
         scf(ix+m1c-19) = -a(4,3)/a(4,2)
         m1c = 22
         ks = 5
         goto 160

*-----------------------------------------------------------------------
*        hyperboloid, ellipsoid, or cone.  let linear term be zero.
*-----------------------------------------------------------------------

  140    scf(ix+m1c-26) = a(4,1)
         scf(ix+7) = a(4,3)-.25*a(4,2)**2/a(4,1)
         scf(ix+m1c-19) = -.5*a(4,2)/a(4,1)

*-----------------------------------------------------------------------
*        two-sheet cone if constant term is zero.
*-----------------------------------------------------------------------

         if(abs(scf(ix+7))*1.e12.gt.t**2.or.a(4,1).gt.0.) goto 150
         scf(ix+1) = scf(ix+m1c-19)
         scf(ix+2) = -a(4,1)
         scf(ix+3) = 0.
         m1c = m1c-8
         ks = 9
         goto 160

*-----------------------------------------------------------------------
*        hyperboloid or ellipsoid.
*-----------------------------------------------------------------------

  150    ks = 6
         if(scf(ix+7).gt.0.) ks = 7
         if(a(4,1).gt.0.) ks = 8

*-----------------------------------------------------------------------
*        check that all points are on the same sheet of a hyperboloid.
*-----------------------------------------------------------------------

      if(scf(ix+m1c-19).lt.t2.and.
     &   scf(ix+m1c-19).gt.t1.and.ks.eq.7) then

               write(io,'(/"** ERROR : in [surface] section,"/
     &         "   the points are on different sheets of a ",
     &             "hyperboloid."/
     &         "   surface id =",i7)')  idsn
               ierr = ierr + 1

      end if

         m1c = 22

*-----------------------------------------------------------------------
*        print the shape of the translated surface.
*-----------------------------------------------------------------------

  160    if(ks.ne.0) then

               write(io,'(/"* surface id =",i7,": -> ",a24)')
     &               idsn, hs(ks)

         end if

      return

*-----------------------------------------------------------------------
*        convert a point-defined plane to a p-plane.
*-----------------------------------------------------------------------

  170    tpp(4) = 0.
      do 180 i = 1, 3
         j = mod(i,3)+1
         k = 6-i-j
         tpp(i) = scf(ix+j)*(scf(ix+k+3)-scf(ix+k+6))
     &          + scf(ix+j+3)*(scf(ix+k+6)-scf(ix+k))
     &          + scf(ix+j+6)*(scf(ix+k)-scf(ix+k+3))
  180    tpp(4) = tpp(4)+scf(ix+i)*(scf(ix+j+3)*scf(ix+k+6)-
     &    scf(ix+j+6)*scf(ix+k+3))

*-----------------------------------------------------------------------
*        check for the points being in line.  put the results in scf.
*-----------------------------------------------------------------------

         xm = 0.

      do 190 i = 1, 4
         if(xm.eq.0..and.tpp(5-i).ne.0.) xm = 1./tpp(5-i)
  190    scf(ix+5-i) = tpp(5-i)*xm

         if(xm.eq.0.) then

               write(io,'(/"** ERROR : in [surface] section,"/
     &         "   cannot convert surface to a plane."/
     &         "   surface id =",i7)')  idsn
               ierr = ierr + 1

         end if

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine trfdef(io,jo,ierr,jt)
*                                                                      *
*       complete the transformation matrix trf(jt).                    *
*       Last modified by K.Niita on 2009/09/30                         *
*                                                                      *
************************************************************************
      use moddas_ggs !frtati20220905

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'ggsparam.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      dimension tr(17),iy(2,3)

      data iy/6,8,7,11,10,12/

      dp(i,j)=tr(i)*tr(j)+tr(i+1)*tr(j+1)+tr(i+2)*tr(j+2)

*-----------------------------------------------------------------------

         nb = 0
         k  = 0

      do i = 1, 14

         tr(i) = trf(i,jt)

         if( tr(i) .ne. huge .and.
     &       i .ge. 5 .and. i .le. 13 ) nb = nb + 1

      end do

         if( tr(14) .eq. 0. ) tr(14) = 1.

      goto (20,320,320,40,320,160,90,320,320,230) nb + 1

*-----------------------------------------------------------------------
*     nb=0 -- no vectors.  simple translation.
*-----------------------------------------------------------------------

   20 do 30 i = 5, 13
   30    tr(i) = trf(i,0)
         goto 260

*-----------------------------------------------------------------------
*     nb=3 -- one vector.
*-----------------------------------------------------------------------

   40 do 60 k = 1, 2
      do 50 j = 1, 3
   50    if(tr(3*j+2).ne.huge.and.tr(3*j+3).ne.huge.and.
     &    tr(3*j+4).ne.huge) goto 70

      do 60 i = 1, 3
         t = tr(iy(1,i))
         tr(iy(1,i)) = tr(iy(2,i))
   60    tr(iy(2,i)) = t
         goto 320

   70    r = sqrt(tr(3*j+2)**2+tr(3*j+3)**2+tr(3*j+4)**2)
         if(r.eq.0.) goto 320

      do 80 i = 2, 4
   80    tr(3*j+i) = tr(3*j+i)/r

         call defbas(tr(3*j+2),tr(3*(mod(j,3)+1)+2),
     1    tr(3*(mod(j+1,3)+1)+2),i)

         go to 140

*-----------------------------------------------------------------------
*     nb=6 -- two vectors in the same system.
*-----------------------------------------------------------------------

   90 do 110 k = 1, 2
      do 100 j = 1, 3
  100    if(tr(3*j+2).eq.huge.and.tr(3*j+3).eq.huge.and.
     &      tr(3*j+4).eq.huge) goto 120
      do 110 i = 1, 3
         t = tr(iy(1,i))
         tr(iy(1,i)) = tr(iy(2,i))
  110    tr(iy(2,i)) = t
         goto 320

  120 do 130 l = 1, 2
         m = mod(j+l-1,3)+1
         r = sqrt(tr(3*m+2)**2+tr(3*m+3)**2+tr(3*m+4)**2)
         if(r.eq.0.) goto 320

      do 130 i=2,4
  130    tr(3*m+i) = tr(3*m+i)/r

         call crspro(tr(3*(mod(j,3)+1)+2),
     &               tr(3*(mod(j+1,3)+1)+2),tr(3*j+2))

  140    if(k.ne.2) goto 260

      do 150 i = 1, 3
         t = tr(iy(1,i))
         tr(iy(1,i)) = tr(iy(2,i))
  150    tr(iy(2,i)) = t
         goto 260

*-----------------------------------------------------------------------
*     nb=5 -- one vector in each system.
*-----------------------------------------------------------------------

  160 do 170 i1 = 1, 3
  170    if(tr(3*i1+2).ne.huge.and.tr(3*i1+3).ne.huge.and.
     &    tr(3*i1+4).ne.huge) goto 180
         goto 320

  180 do 190 j1 = 1, 3
  190    if(tr(j1+4).ne.huge.and.tr(j1+7).ne.huge.and.
     &    tr(j1+10).ne.huge) goto 200
         goto 320

  200    r = tr(3*i1+j1+1)**2
      do 210 i = 5, 13
  210    if(tr(i).ne.huge) r = r+tr(i)**2
         r = sqrt(.5*r)
         if(r.eq.0.) goto 320

      do 220 i = 5, 13
  220    if(tr(i).ne.huge) tr(i) = tr(i)/r
         if(tr(3*i1+j1+1)**2.eq.1.) goto 320
         i2 = mod(i1,3)+1
         i3 = mod(i2,3)+1
         j2 = mod(j1,3)+1
         j3 = mod(j2,3)+1
         r = 1./(1.-tr(3*i1+j1+1)**2)
         a = tr(3*i1+j3+1)*tr(3*i3+j1+1)*r
         b = tr(3*i1+j2+1)*tr(3*i2+j1+1)*r
         c = tr(3*i2+j1+1)*tr(3*i1+j3+1)*r
         d = tr(3*i3+j1+1)*tr(3*i1+j2+1)*r
         tr(3*i2+j2+1) = -a-b*tr(3*i1+j1+1)
         tr(3*i2+j3+1) = d-c*tr(3*i1+j1+1)
         tr(3*i3+j2+1) = c-d*tr(3*i1+j1+1)
         tr(3*i3+j3+1) = -b-a*tr(3*i1+j1+1)
         goto 260

*-----------------------------------------------------------------------
*     nb=9 -- entire b matrix.
*-----------------------------------------------------------------------

  230 do 250 j = 1, 3
         r = 0.
      do 235 i = 1, 3
  235    r = r+tr(3*j+i+1)**2
         if(r.eq.0.) goto 320
         s = 1./sqrt(r)
      do 240 i = 1, 3
  240    tr(3*j+i+1) = tr(3*j+i+1)*s
  250    continue

*-----------------------------------------------------------------------
*     set flag (trf(ltrf+1,jt) > 0 to indicate skew rotation.
*-----------------------------------------------------------------------

  260 do 265 i = 5, 13
  265    if(tr(i).gt.epss.and.abs(1.-tr(i)).gt.epss) tr(1) = abs(tr(1))

*-----------------------------------------------------------------------
*     transpose matrix and check orthogonality.
*-----------------------------------------------------------------------

         r = 0.
      do 270 i = 1, 3
  270    r = max(r,abs(dp(3*i+2,3*(mod(i,3)+1)+2)),
     &             abs(dp(3*i+2,3*i+2)-1.))
         if(k.eq.2) goto 290

      do 280 i = 1, 3
         t = tr(iy(1,i))
         tr(iy(1,i)) = tr(iy(2,i))
  280    tr(iy(2,i)) = t

  290    continue

         if(r.gt.2.e-6) then

               write(io,'("** warning : in [transform] section,"/
     &         "   non-orthogonality of transformation."/
     &         "   transform id =",i8,":  r =",e10.2,
     &         " > 2.e-6")')
     &         int(abs(tr(1))), r

               ErrCha = ''
               MsgID = 'L:2222/R:trfdef/F:ggs00.f'
               call ErrWrite(MsgID, ErrCha)
               write(jo,'("** warning : in [transform] section,"/
     &         "   non-orthogonality of transformation."/
     &         "   transform id =",i8,":  r =",e10.2,
     &         " > 2.e-6")')
     &         int(abs(tr(1))), r

         end if

         if(r.gt.2.e-6) tr(1) = abs(tr(1))

*-----------------------------------------------------------------------
*     check for degenerate rotation matrix.
*-----------------------------------------------------------------------

         rr = dp(5,5)
         rs = dp(5,8)
         ra = sqrt(rr)

         rb = sqrt(max(zero,dp(8,8)*rr**2-rr*rs**2))

         if(ra.lt.epss.or.rb.lt.epss) tr(14) = -2.

         if(tr(14).eq.-2.) then

               write(io,'(/"** ERROR : in [transform] section,"/
     &         "   transformation is degenerate."/
     &         "   transform id =",i8)')
     &         int(tr(1))

               ErrCha = ''
               MsgID = 'L:2254/R:trfdef/F:ggs00.f'
               call ErrWrite(MsgID, ErrCha)
               write(jo,'(/"** ERROR : in [transform] section,"/
     &         "   transformation is degenerate."/
     &         "   transform id =",i8)')
     &         int(tr(1))

               ierr = ierr + 1

         end if

         if(tr(14).eq.-2.) goto 330

*-----------------------------------------------------------------------
*     ensure rotation matrix is orthonormal.
*-----------------------------------------------------------------------

      do 300 i = 1, 3
         tr(i+7) = (tr(i+7)*rr-tr(i+4)*rs)/rb
  300    tr(i+4) = tr(i+4)/ra

         call crspro(tr(5),tr(8),tpp)

         a = sign(one,tpp(1)*tr(11)+tpp(2)*tr(12)+tpp(3)*tr(13))
         tr(11) = a*tpp(1)
         tr(12) = a*tpp(2)
         tr(13) = a*tpp(3)
      do 305 i = 5, 13
  305    if(abs(tr(i)).lt.epss) tr(i) = 0.

*-----------------------------------------------------------------------
*     generate the other inter-origin vector.
*-----------------------------------------------------------------------

      do 310 i = 1, 3
         tr(14+i) = tr(1+i)
  310    tpp(i) = -tr(1+i)

         if(tr(14).ne.-1.) call matmpy(3,3,1,tr(5),3,tpp,3,tr(2),3,0)
         if(tr(14).eq.-1.) call matmpy(3,3,1,tr(5),3,tpp,3,tr(15),3,1)

         goto 330

*-----------------------------------------------------------------------
*     bad b matrix.
*-----------------------------------------------------------------------

  320 continue

               write(io,'(/"** ERROR : in [transform] section,"/
     &         "   transformation is incorrectly defined."/
     &         "   transform id =",i8)')
     &         int(abs(tr(1)))

               ErrCha = ''
               MsgID = 'L:2309/R:trfdef/F:ggs00.f'
               call ErrWrite(MsgID, ErrCha)
               write(jo,'(/"** ERROR : in [transform] section,"/
     &         "   transformation is incorrectly defined."/
     &         "   transform id =",i8)')
     &         int(abs(tr(1)))

               ierr = ierr + 1

         tr(14) = -2.

*-----------------------------------------------------------------------

  330 do 340 i=1,17
  340    trf(i,jt) = tr(i)

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine trfsuf(io,ierr,js)
*                                                                      *
*       transform surface js into the xyz coordinate system.           *
*       Last modified by K.Niita on 2009/09/30                              *
*                                                                      *
************************************************************************
      use moddas_ggs !frtati20220905

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'ggsparam.inc'

*-----------------------------------------------------------------------

      dimension tm(4,4),aa(4,4),nc(23)
      data tm/1.,15*0./,nc/4,4*1,4,3*2,3*3,3*1,3*5,3*3,2*10/

*-----------------------------------------------------------------------

      do 10 jt = 1, mxtr
   10    if(abs(trf(1,jt)).eq.jtr(ljtr+js)) goto 20

         write(io,'(/"** ERROR : in GG setup,"/
     &   "transformation",i6," for surface ",i7," is absent.")')
     &   jtr(js),nsf(js)

         ierr = ierr + 1
         return

   20    if(trf(14,jt).eq.-2.) return
         ks = kst(js)
         if(abs(ks-25).le.1) goto 180

*-----------------------------------------------------------------------
*     generate and transform the aa matrix.
*-----------------------------------------------------------------------

      do 30 i = 1, 4
      do 30 j = 2, 4
   30    tm(j,i) = trf(3*i+j-3,jt)

         call sufdef(js,aa)

         tpp(18) = aa(2,2)
         tpp(19) = aa(3,3)
         tpp(20) = aa(4,4)
         call matmpy(4,4,4,aa,4,tm,4,tpp,4,0)
         call matmpy(4,4,4,tm,4,tpp,4,aa,4,4)

*-----------------------------------------------------------------------
*        set up to convert to the simplest possible surface.
*        zero out diagonal terms where possible.
*-----------------------------------------------------------------------

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

*-----------------------------------------------------------------------
*        put gq coefficients into tpp(1-10) and attempt to put
*        sq coefficients into tpp(11-20)
*-----------------------------------------------------------------------

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

*-----------------------------------------------------------------------
*        use gq surface if necessary.
*-----------------------------------------------------------------------

         if(ks.lt.10.or.ks.lt.23.and.trf(1,jt).lt.0.) goto 60

         if(ks.ge.16.and.ks.le.21.and.
     &      scf(lsc(js)+nc(ks)).ne.0.) then

            write(io,'(/"** ERROR : in GG setup,"/
     &      "one-sheet cone surface ",i7,
     &      " has a skew transformation.")')
     &      nsf(js)

            ierr = ierr + 1

         end if

         ks = 23
         tpp(4) = 2.*aa(2,3)
         tpp(5) = 2.*aa(3,4)
         tpp(6) = 2.*aa(2,4)
         if(aa(2,3).ne.0..or.aa(3,4).ne.0..or.aa(2,4).ne.0.) goto 130
         ks = 22
         goto 110

*-----------------------------------------------------------------------
*        produce a planar or spherical surface.
*-----------------------------------------------------------------------

   60    if(ks.ge.10) goto 80
         kk = 1

      do 70 i= 1 , 3
         if(ks.le.4) tpp(i) = 2.*tpp(i+13)
         if(ks.ge.5) tpp(i) = tpp(i+17)
   70    if(abs(tpp(i)).gt.abs(tpp(kk))) kk = i
         if(abs(tpp(1))+abs(tpp(2))+abs(tpp(3))-abs(tpp(kk)).gt.
     &      epss*(abs(tpp(kk))+abs(tpp(17))).or.
     &      ks.le.4.and.tpp(kk).lt.0.) kk = 0
         if(ks.ge.5) kk = kk+5
         ks = kk+1
         if(ks.gt.6)then
           if(abs(tpp(kk-5)).lt.-epss*tpp(17)) ks = 5
         endif
         if(abs(ks-3).le.1) tpp(17) = tpp(17)/tpp(kk)
         if(ks.gt.6) tpp(1) = tpp(kk-5)
         tpp(nc(ks)) = -tpp(17)
         goto 130

*-----------------------------------------------------------------------
*        produce a cylindrical surface.
*-----------------------------------------------------------------------

   80    if(ks.ge.16) goto 90
         tpp(1) = tpp(18)
         tpp(2) = tpp(20)
         if(im.eq.1) tpp(1) = tpp(19)
         if(im.eq.3) tpp(2) = tpp(19)
         ks = im+9
         if(abs(tpp(1))+abs(tpp(2)).lt.-epss*tpp(17)) ks = im + 12
         tpp(nc(ks)) = -tpp(17)
         goto 130

*-----------------------------------------------------------------------
*        produce a conical surface.
*-----------------------------------------------------------------------

   90    if(ks.eq.22) goto 110
         k = lsc(js)+nc(ks)
         ic = 1
         if(tpp(1)*tpp(2).gt.0.) ic = 3
         if(tpp(1)*tpp(3).gt.0.) ic = 2
         j = mod(ks-16,3)+3*ic+2
         ks = ic+15

      do 100 i = 1, 3
  100    tpp(i) = tpp(17+i)

         if(abs(tpp(1))+abs(tpp(2))+abs(tpp(3))-abs(tpp(ic)).lt.
     &      epss*(abs(tpp(ic))+1.)) ks = ic+18
         if(ks.gt.18) tpp(1) = tpp(ic)
         tpp(nc(ks)-1) = scf(k-1)
         tpp(nc(ks)) = scf(k)*sign(one,trf(j,jt))
         goto 130

*-----------------------------------------------------------------------
*        produce a special quadratic.
*-----------------------------------------------------------------------

  110 do 120 i = 1, 10
  120    tpp(i) = tpp(10+i)

*-----------------------------------------------------------------------
*        load surface coefficients into the scf array.
*-----------------------------------------------------------------------

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

*-----------------------------------------------------------------------
*        transform torus if transformation is not skewed.
*-----------------------------------------------------------------------

  180    continue

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      subroutine sufdef(js,aa)
*                                                                      *
*       put surface js into the form  (1,x,y,z)*aa(4x4)*(1,x,y,z)t=0   *
*       Last modified by K.Niita on 2001/02/28                         *
*                                                                      *
************************************************************************
      use moddas_ggs !frtati20220905

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'ggsparam.inc'

*-----------------------------------------------------------------------

      dimension aa(4,4)

*-----------------------------------------------------------------------

      do 10 j = 1, 4
      do 10 i = 1, 4
   10    aa(i,j) = 0.

         l = lsc(js)+1
         k = kst(js)

      do 20 i = 2, 4
   20    if(k.ge.5.and.k.le.21) aa(i,i) = 1.

         goto (30,40,40,40,50,60,70,70,70,80,80,80,90,90,90,
     &    100,100,100,110,110,110,120,130) k

*-----------------------------------------------------------------------
* >>>>>  p
*-----------------------------------------------------------------------

   30    aa(1,1) = -scf(l+3)
         aa(1,2) = .5*scf(l)
         aa(1,3) = .5*scf(l+1)
         aa(1,4) = .5*scf(l+2)
         goto 140

*-----------------------------------------------------------------------
* >>>>>  px py pz
*-----------------------------------------------------------------------

   40    aa(1,1) = -scf(l)
         aa(1,k) = .5
         goto 140

*-----------------------------------------------------------------------
* >>>>>  so
*-----------------------------------------------------------------------

   50    aa(1,1) = -scf(l)
         goto 140

*-----------------------------------------------------------------------
* >>>>>  s
*-----------------------------------------------------------------------

   60    aa(1,1) = scf(l)**2+scf(l+1)**2+scf(l+2)**2-scf(l+3)
         aa(1,2) = -scf(l)
         aa(1,3) = -scf(l+1)
         aa(1,4) = -scf(l+2)
         goto 140

*-----------------------------------------------------------------------
* >>>>>  sx sy sz
*-----------------------------------------------------------------------

   70    aa(1,1) = scf(l)**2-scf(l+1)
         aa(1,k-5) = -scf(l)
         goto 140

*-----------------------------------------------------------------------
* >>>>>  c/x c/y c/z
*-----------------------------------------------------------------------

   80    aa(1,1) = scf(l)**2+scf(l+1)**2-scf(l+2)
         aa(1,2) = -scf(l)
         aa(1,3) = -scf(l+k/12)
         aa(1,4) = -scf(l+1)
         aa(1,k-8) = 0.
         aa(k-8,k-8) = 0.
         goto 140

*-----------------------------------------------------------------------
* >>>>>  cx cy cz
*-----------------------------------------------------------------------

   90    aa(1,1) = -scf(l)
         aa(k-11,k-11) = 0.
         goto 140

*-----------------------------------------------------------------------
* >>>>>  k/x k/y k/z
*-----------------------------------------------------------------------

  100    aa(1,1) = -(1.+scf(l+3))*scf(l+k-16)**2
     &           + scf(l)**2+scf(l+1)**2+scf(l+2)**2
         aa(1,2) = -scf(l)
         aa(1,3) = -scf(l+1)
         aa(1,4) = -scf(l+2)
         aa(1,k-14) = scf(l+k-16)*scf(l+3)
         aa(k-14,k-14) = -scf(l+3)
         goto 140

*-----------------------------------------------------------------------
* >>>>>  kx ky kz
*-----------------------------------------------------------------------

  110    aa(1,1) = -scf(l+1)*scf(l)**2
         aa(1,k-17) = scf(l)*scf(l+1)
         aa(k-17,k-17) = -scf(l+1)
         goto 140

*-----------------------------------------------------------------------
* >>>>>  sq
*-----------------------------------------------------------------------

  120    aa(1,1) = scf(l)*scf(l+7)**2+scf(l+1)*scf(l+8)**2+scf(l+2)*
     &    scf(l+9)**2-2.*(scf(l+3)*scf(l+7)+scf(l+4)*scf(l+8)+
     &    scf(l+5)*scf(l+9))+scf(l+6)

         aa(1,2) = -scf(l)*scf(l+7)+scf(l+3)
         aa(1,3) = -scf(l+1)*scf(l+8)+scf(l+4)
         aa(1,4) = -scf(l+2)*scf(l+9)+scf(l+5)
         aa(2,2) = scf(l)
         aa(3,3) = scf(l+1)
         aa(4,4) = scf(l+2)
         goto 140

*-----------------------------------------------------------------------
* >>>>>  gq
*-----------------------------------------------------------------------

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

*-----------------------------------------------------------------------
*        complete the symmetric aa matrix.
*-----------------------------------------------------------------------

  140 do 150 i = 2, 4
      do 150 j = 2, i
  150    aa(i,j-1) = aa(j-1,i)

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      subroutine chkcdf(io,ierr)
*                                                                      *
*       simplify and check the cell descriptions.                      *
*       Last modified by K.Niita on 2009/09/30                              *
*                                                                      *
************************************************************************
      use GGBANKMOD !FURUTA
      use LAFDATAMOD !FURUTA20201127
      use moddas_ggs !frtati20220905

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'ggsparam.inc'

      include 'err.inc'

*-----------------------------------------------------------------------

      common /inggs/ iog, igcel, ioa, igsuf, iob, igtrs
      common /geomemcom/ igeomem  ! T.Sato 2024/12/25
*-----------------------------------------------------------------------

      dimension iu(0:mxlv),il(0:mxlv),mu(0:mxlv),jm(0:mxlv),nu(0:mxlv)
      character ha*1,hb*2,hc*2
      data ha/'<'/

      nlev=0 !FURUTA

*-----------------------------------------------------------------------
*        simplify and expand the cell descriptions.
*-----------------------------------------------------------------------

      do 260 ic = 1, mxa
         m1 = abs(lca(ic))
   10    m2 = abs(lca(ic+1))-1

*-----------------------------------------------------------------------
*        find and check the next #n descriptor.
*-----------------------------------------------------------------------

      do 30 j = m1, m2
         j1 = lja(j+1)
         if(lja(j).ne.1000004.or.j1.eq.1000001) goto 30

         if(j1.le.1000000.and.j.ne.m2.and.j1.ne.ncl(ic)) goto 20

            write(io,'(/"** ERROR : in GG setup,"/
     &      "# used incorrectly in description of cell",i7)')
     &      ncl(ic)

            ierr = ierr + 1

         goto 250

   20    n = namchg(1,j1)

         if(n.ne.0) goto 40

            write(io,'(/"** ERROR : in GG setup,"/
     &      "cell",i7," in description of cell",i7,
     &      " is not defined.")')
     &      j1,ncl(ic)

            ierr = ierr + 1

         goto 250

   30    continue
         goto 80

*-----------------------------------------------------------------------
*        expand the #n descriptor.
*-----------------------------------------------------------------------

   40    lh = abs(lca(n+1))-abs(lca(n))

      do 50 iz = 1, nlja - j
         i = nlja-j+1-iz
         if(j+lh+2+i.gt.mlja+igeomem) then
          write(ErrCha,'("Too complicated cell definition. ",
     &    " Reduce # in a cell or increase igeomem up to",i8,
     &    " or more")') nlja
          MsgID = 'L:2797/R:chkcdf/F:ggs00.f' !E06_001_001
          call ErrWrite(MsgID, ErrCha)
          stop
         endif
         idna(j+lh+2+i) = idna(j+1+i)
   50    lja(j+lh+2+i) = lja(j+1+i)
         lja(j+1) = 1000001
         lja(j+lh+2) = 1000002
         idna(j+1) = 0
         idna(j+lh+2) = 0
         nlja = nlja+lh+1
         nljc = nljc+lh+1

      do 60 i = ic, mxa
   60    lca(i+1) = lca(i+1)+sign(lh+1,lca(i+1))

      do 70 i = 1, lh
         idna(j+1+i) = idna(abs(lca(n))+i-1)
   70    lja(j+1+i) = lja(abs(lca(n))+i-1)
         goto 10

*-----------------------------------------------------------------------
*        find the next # complementation operator.
*-----------------------------------------------------------------------

   80 do 100 j = m1, m2
         if(lja(j).ne.1000004) goto 100
         np = 0

      do 90 k = j, m2
         if(lja(k+1).eq.1000001) np = np+1
         if(lja(k+1).eq.1000002) np = np-1
         if(np.eq.0) goto 110
   90    if(lja(k+1).eq.1000004) goto 100
         goto 250
  100 continue
         goto 180

*-----------------------------------------------------------------------
*        multiply the complementation operator out.
*-----------------------------------------------------------------------

  110 do 120 i = 1, nlja - j + 1
         idna(j+i-1) = idna(j+i)
  120    lja(j+i-1) = lja(j+i)
         nlja = nlja-1

      do 130 i = ic, mxa
  130    lca(i+1) = lca(i+1)-sign(1,lca(i+1))

         np = 0
  140    if(lja(j).eq.1000001) np = np+1
         if(lja(j).eq.1000002) np = np-1
         if(np.eq.0) goto 10

         if(lja(j).le.1000000) lja(j) = -lja(j)

         if(lja(j).eq.1000001.or.lja(j).ne.1000003.and.
     &    lja(j+1).gt.1000001) goto 170

      do 150 iz = 1, nlja - j + 1
         i = nlja-j+2-iz
         if(j+i+1.gt.mlja+igeomem) then
          write(ErrCha,'("Too complicated cell definition. ",
     &    " Reduce # in a cell or increase igeomem up to",i8,
     &    " or more")') nlja
          MsgID = 'L:2863/R:chkcdf/F:ggs00.f' !E06_001_001
          call ErrWrite(MsgID, ErrCha)
          stop
         endif
         idna(j+i+1) = idna(j+i)
  150    lja(j+1+i) = lja(j+i)

      do 160 i = ic, mxa
  160    lca(i+1) = lca(i+1)+sign(1,lca(i+1))
         idna(j+1) = 0
         if(lja(j).ne.1000003) lja(j+1) = 1000003
         if(lja(j).eq.1000003) lja(j+1) = 1000001
         if(lja(j).eq.1000003) lja(j) = 1000002
         nlja = nlja+1
         nljc = nljc+1
         m2 = m2+1
         j = j+1
  170    j = j+1
         if(j.le.m2) goto 140
         goto 250

*-----------------------------------------------------------------------
*        remove any redundant parentheses from the cell description.
*-----------------------------------------------------------------------

  180 do 200 ip = m1, m2
         if(lja(ip).ne.1000001) goto 200
         ng = 0
         np = 0

         do 190 l = ip, m2
         if(lja(l).eq.1000001) np = np + 1
         if(lja(l).eq.1000003.and.np.eq.1) ng = 1
         if(lja(l).eq.1000002) np = np - 1
         if(np.eq.0.and.(ng.eq.0.or.ip.eq.m1.and.l.eq.m2)) goto 210

  190    if(np.eq.0) goto 200
         goto 250

  200    continue
         goto 250

  210 do 220 i = 1, nlja-ip+1
         if(ip+i.gt.mlja+igeomem) then
          write(ErrCha,'("Too complicated cell definition. ",
     &    " Reduce # in a cell or increase igeomem up to",i8,
     &    " or more")') nlja
          MsgID = 'L:2910/R:chkcdf/F:ggs00.f' !E06_001_001
          call ErrWrite(MsgID, ErrCha)
          stop
         endif
         idna(ip+i-1) = idna(ip+i)
  220    lja(ip+i-1) = lja(ip+i)

      do 230 i = 1, nlja-l+1
         idna(l+i-2) = idna(l+i-1)
  230    lja(l+i-2) = lja(l+i-1)
         nlja = nlja-2

      do 240 i = ic, mxa
  240    lca(i+1) = lca(i+1)-sign(2,lca(i+1))
         goto 10

  250    if(m2-m1.ge.mlgc) then

             write(io,'(/"** ERROR : in GG setup,"/
     &      "description of cell",i7," uses",i4," words.  ",i4,
     &      " is max.")')
     &      ncl(ic),m2-m1+1,mlgc

            ierr = ierr + 1

         end if

  260    continue

*-----------------------------------------------------------------------
*        check the definitions of the surfaces of all the cells.
*-----------------------------------------------------------------------

      do 290 ic = 1, mxa
      do 280 la = abs(lca(ic)), abs(lca(ic+1))-1

         if(lja(la).gt.1000000) goto 280

      do 270 jb = 1, mxj
  270    if(abs(lja(la)).eq.nsf(jb)) goto 280

            write(io,'(/"** ERROR : in GG setup,"/
     &      "surface",i7," of cell ",i7," is not defined.")')
     &      lja(la),ncl(ic)

            ierr = ierr + 1

  280 continue
  290 continue

*-----------------------------------------------------------------------
*        change the contents of lja from names to indexes.
*-----------------------------------------------------------------------

      do 300 j = 1, nlja
         j1 = lja(j)

         if(j1.gt.1000000) goto 300

         j2 =namchg(2,abs(j1))
         lja(j) = sign(j2+abs(idna(j)),j1)
  300 continue

*-----------------------------------------------------------------------
*        check for unused tr cards.
*-----------------------------------------------------------------------

      do 305 m = 1, mxtr
         i = nint(abs(trf(1,m)))
         if(i.gt.1000) goto 305

      do 301 j = 1, mxj
  301    if(jtr(j).eq.i) goto 305
         if(junf.eq.0) goto 304

      do 303 j = 1, mxa
         k = -mfl(1,j)
         if(k.le.0) goto 303

      do 302 n = 3, 2+laf(1,k+2)*laf(2,k+2)*laf(3,k+2) !FURUTA20201127
  302    if(laf(3,k+n).eq.m) goto 305                  !FURUTA20201127

  303    if(mfl(3,j).eq.m.or.ktr(j).eq.i) goto 305
  304       write(io,'("** warning : in GG setup,"/
     &      "tr",i6," card unused.")')
     &      i

  305 continue

*-----------------------------------------------------------------------
*        prepare to delete all identical surfaces.  the surface is
*        not physically deleted, but is removed from cell cards and
*        tallies so it is not used anywhere.  also deleted if the
*        same except for sign, then sign is swiched on cell cards.
*-----------------------------------------------------------------------
*        set up arrays involving identical surfaces.
*-----------------------------------------------------------------------

         mc = 0
         in = 0

      do 895 ja = 1, mxj - 1
         if(ksm(ja).lt.0) goto 895
         if(idnt(ja).ne.0) goto 895
         na = lsc(ja+1)-lsc(ja)

      do 870 j = ja+1,mxj
         if(ksm(j).lt.0) goto 870
         if(idnt(j).ne.0) goto 870
         is = 1
         if(kst(ja).gt.4.or.kst(j).gt.4) goto 820

*-----------------------------------------------------------------------
*        check for identical planes.
*-----------------------------------------------------------------------

         if(kst(ja).eq.1.and.kst(j).eq.1) goto 810
         if(kst(ja).eq.1.or.kst(j).eq.1) goto 800
         if(kst(ja).ne.kst(j)) goto 870
         if(abs(scf(lsc(ja)+1)-scf(lsc(j)+1)).gt.
     &    1e-12*abs(scf(lsc(ja)+1)+scf(lsc(j)+1))) goto 870
         goto 840

  800    j1 = ja
         j2 = j
         if(kst(ja).ne.1) j1 = j
         if(kst(ja).ne.1) j2 = ja
         j3 = kst(j2)-1
         if(scf(lsc(j1)+j3).eq.0.) goto 870

      do 805 k = 1, 3
  805    if(k.ne.j3.and.
     &      abs(scf(lsc(j1)+k)/scf(lsc(j1)+j3)).gt.
     &      1.d-12) goto 870
         d = scf(lsc(j1)+4)/scf(lsc(j1)+j3)
         if(abs(scf(lsc(j2)+1)-d).gt.
     &      1.d-12*abs(scf(lsc(j2)+1)+d)) goto 870
         if(scf(lsc(j1)+j3).lt.0.) is = -1
         goto 840

  810    j3 = 1
         g = abs(scf(lsc(ja)+1))
         if(abs(scf(lsc(ja)+2)).gt.g) j3 = 2
         if(abs(scf(lsc(ja)+2)).gt.g) g = abs(scf(lsc(ja)+2))
         if(abs(scf(lsc(ja)+3)).gt.g) j3 = 3
         if(abs(scf(lsc(ja)+3)).gt.g) g = abs(scf(lsc(ja)+3))
         if(g.eq.0.) goto 870
         g = scf(lsc(j)+j3)/scf(lsc(ja)+j3)

      do 815 k = 1, 4
         if(k.eq.j3) goto 815
         d = g*scf(lsc(ja)+k)
         if(abs(scf(lsc(j)+k)-d).gt.
     &    1.d-12*abs(scf(lsc(j)+k)+d)) goto 870
  815 continue

         if(g.lt.0.) is = -1
         goto 840

  820    if(kst(ja).ne.kst(j)) goto 870

*-----------------------------------------------------------------------
*        check surfaces other than planes.
*-----------------------------------------------------------------------

         if(kst(ja).le.21.or.kst(ja).ge.24) goto 830

      do 825 i = 1, na
         if(scf(lsc(ja)+i).eq.0.) goto 825
         if(scf(lsc(ja)+i)*scf(lsc(j)+i).lt.0.) is = -1
         goto 830
  825 continue

  830 do 835 i = 1, na
  835    if(abs(scf(lsc(ja)+i)*is-scf(lsc(j)+i)).gt.
     &    1.d-12*abs(scf(lsc(ja)+i)*is
     &    +scf(lsc(j)+i))) goto 870
         if( kst(ja).ge.24 .and. kst(ja).le.26
     &       .and. jtr(ja).ne.jtr(j) ) goto 870 ! frtati 2022/04/08
  840    hb = ' '
         ka = nsf(ja)
         mq = kfq
         if(kfq.ne.0) write(hb,'(1h.,i1)') kfq
         hc = ' '
         k = nsf(j)
         if(kfq.ne.0) write(hc,'(1h.,i1)') kfq

         write(io,850) ka,hb,k,hc,k,hc
  850    format(8h surface,i8,a2,12h and surface,i8,a2,
     &    14h are the same.,i7,a2,17h will be deleted.)

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

            write(io,'(/"C** Warning : in GG setup,"/
     &      "   boundary condition on identical surface (1).")')

         end if

  860    idnt(j) = idnt(ja)*is
         jq = idns(ja)
         idns(j) = jq
         idne(jq) = idne(jq)+1
         idne(2) = idne(2)+1
         idne(idne(2)) = j

         if(ksu(j).ne.0) then

            write(io,'(/"C** Warning : in GG setup,"/
     &      "   boundary condition on identical surface (2).")')

         end if

         mc = mc + 1

*-----------------------------------------------------------------------
*        change deleted surface on cell cards.
*-----------------------------------------------------------------------

      do 865 k = 1, mlja
         if(abs(lja(k)).ne.j) goto 865
         lja(k) = sign(ja,is*lja(k))
         idna(k) = mq
  865 continue
  870 continue

*-----------------------------------------------------------------------
*     make master identical surface a non-facet if possible.
*-----------------------------------------------------------------------

         if(idnt(ja).eq.0) goto 895
         if(ksm(ja).eq.0) goto 895

      do 875 jq = 2, idne(idns(ja))
         j1 = idne(idns(ja)+jq)
  875    if(ksm(j1).eq.0) goto 880
         goto 895

  880    is = 1
         if(idnt(j1).lt.0) is = -1
         iq = j1*is
         idnt(j1) = j1
         idnt(ja) = j1*is
         idne(idns(ja)+jq) = ja
         idne(idns(ja)+1) = j1

      do 885 k = 1, mlja
         if(abs(lja(k)).ne.ja) goto 885
         lja(k) = sign(j1,is*lja(k))
         idna(k) = 0
  885 continue

         if(idne(idns(ja)).le.2) goto 895

      do 890 j2 = 2,idne(idns(ja))
         j3 = idne(idns(ja)+j2)
         if(ja.eq.j3) goto 890
         j4 = idnt(j3)
         if(j4.gt.0) idnt(j3) = iq
         if(j4.lt.0) idnt(j3) = -iq
  890 continue

*-----------------------------------------------------------------------

  895 continue

      if(mc.ne.0) then

         write(io,'("** warning : in GG setup,"/
     &   i4," surfaces were deleted for being the same as others.")')
     &   mc

      write (io,*)

      ErrCha = ''
      ErrID = 'L:3194/R:chkcdf/F:ggs00.f' !E06_001_001
      call ErrWriteIO(ErrID,ErrCha,io)

      write(io,*) 'Please check the section [cell]'
      write(io,*) 'There is a possibility of cross referencing cells'
      write(io,*) 'Please check the cells using #'
c ----------------------------------------------------------------------

      end if


*-----------------------------------------------------------------------
*     check whether there are too many levels of universes.
*-----------------------------------------------------------------------

         if(junf.eq.0) return

      do 490 ic = 1, mxa

         if(mfl(1,ic).ne.0.and.mfl(1,ic).eq.abs(jun(ic)))
     &    then

            write(io,'(/"** ERROR : in GG setup,"/
     &      "   cell = ",i7," is filled with its own universe.")')
     &           ncl(ic)

            ierr = ierr + 1
            return

          end if

         if(lat(1,ic).eq.0.and.mfl(1,ic).lt.0) then

            write(io,'(/"** ERROR : in GG setup,"/
     &      "   non lattice cell = ",i7," has fill array.")')
     &           ncl(lncl+ic)

            ierr = ierr + 1
            return

         end if

*-----------------------------------------------------------------------

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

            write(io,'(/"** ERROR : in GG setup,"/
     &      "   universe = ",i7," which fills cell = ",i7,
     &      " has no cells.")')
     &           lu, ncl(ic)

            ierr = ierr + 1
            return

         end if

  520 continue

         nlv(ic) = max(l,1)
         nlev = max(nlev,l)
         goto 500

  530    continue

         if(nlev.gt.mxlv) then

            write(io,'(/"** ERROR : in GG setup,"/
     &      "   level of univers = ",i7," but max. is = ",i7)')
     &           nlev,mxlv

            ierr = ierr + 1
            return

         end if

*-----------------------------------------------------------------------
*     check that no surface appears more than once in any chain.
*-----------------------------------------------------------------------
*     Special for reducing cpu time for huge lattice
*     by KN
*-----------------------------------------------------------------------

         goto 670

*-----------------------------------------------------------------------

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

         write(io,'(/"** ERROR : in GG setup,"/
     &   "   surface = ",i8,a2,
     &   " appears more than once in a chain.")')
     &        ka, hb

         ierr = ierr + 1
         return

      else

         write(io,'(/"** Warning : in GG setup,"/
     &   "   surface = ",i8,a2,
     &   " appears more than once in a chain.")')
     &        ka, hb

      end if

         write(io,560) ka,hb,ncl(il(lv)),ncl(il(l)),
     &    ncl(il(lv)),('<',ncl(il(lv-i)),i=1,lv)
  560    format(/8h surface,i7,a2,12h is in cells,i6,4h and,i6,
     &    9h in chain/ i6,15(2x,a1,i5))

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

*-----------------------------------------------------------------------
*     set flags (ksc) for parallel and possibly coincident surfaces.
*     ksc(lksc+j)=0 nonplanar; =2/3/4 for px/py/pz; >4 for p.
*-----------------------------------------------------------------------

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
         if((scf(is)*scf(ik)+scf(is+1)*scf(ik+1)
     &      +scf(is+2)*scf(ik+2))**2.gt..9998*ts*tk)
     &      ksc(jk) = ksc(js)
  680 continue
  690 continue

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      function namchg(mm,ji)
*                                                                      *
*       return the array index corresponding to the name ji.           *
*       Last modified by K.Niita on 2009/09/30                         *
*                                                                      *
************************************************************************
      use moddas_ggs !frtati20220905

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'ggsparam.inc'

*-----------------------------------------------------------------------

      goto (10,30,50) mm

*-----------------------------------------------------------------------
* >>>>>  mm=1 -- cell
*-----------------------------------------------------------------------

   10 do 20 namchg = 1, mxa
   20    if(ncl(namchg).eq.ji) return
         goto 70

*-----------------------------------------------------------------------
* >>>>>  mm=2 -- surface
*-----------------------------------------------------------------------

   30 do 40 namchg = 1,mxj
   40    if(nsf(namchg).eq.abs(ji)) return
         goto 70

*-----------------------------------------------------------------------
* >>>>>  mm=3 -- bounding surface
*-----------------------------------------------------------------------

   50 do 60 kk = 1, nlja
         namchg = abs(lja(kk))

         if(namchg.gt.1000000) goto 60

         if(nsf(namchg).eq.abs(ji)) return
   60 continue

*-----------------------------------------------------------------------
*     namchg=0 means that no valid name was found.
*-----------------------------------------------------------------------

   70    namchg = 0

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      subroutine defbas(a,b,c,ie)
*                                                                      *
*       provide 2 arbitrary base vectors b and c perpendicular to a.   *
*       the length of a is assumed to be 1.                            *
*       make abc a right-handed system.  return with ie=1 if a=0.      *
*       b is made catawampous to avoid singularities.                  *
*                                                                      *
*       Last modified by K.Niita on 2009/09/30                              *
*                                                                      *
************************************************************************

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      dimension a(3),b(3),c(3)

*-----------------------------------------------------------------------

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

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      subroutine crspro(v1,v2,v3)
*                                                                      *
*       calculate the cross product   v3 = v1 x v2                     *
*                                                                      *
*       Last modified by K.Niita on 2009/09/30                         *
*                                                                      *
************************************************************************

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      dimension v1(3),v2(3),v3(3)

*-----------------------------------------------------------------------

         v3(1) = v1(2)*v2(3)-v1(3)*v2(2)
         v3(2) = v1(3)*v2(1)-v1(1)*v2(3)
         v3(3) = v1(1)*v2(2)-v1(2)*v2(1)

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      subroutine matmpy(m,n,l,a,ia,b,ib,c,ic,jt)
*                                                                      *
*       multiply a (or a(transpose)) and b to give c.                  *
*       m = number of rows of a (or a(t)) and number of rows of c.     *
*       n = number of columns of a (or a(t)) and number of rows of b.  *
*       l = number of columns of b and number of columns of c.         *
*       jt = 0 for c=a*b    jt.ne.0 for c=a(t)*b                       *
*                                                                      *
*       Last modified by K.Niita on 2009/09/30                         *
*                                                                      *
************************************************************************

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      dimension a(ia,*),b(ib,*),c(ic,*)

*-----------------------------------------------------------------------

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

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      function lgeval(lg,n)
*                                                                      *
*       evaluate logical expression lg of length n.                    *
*       representation of logical elements:                            *
*          (      1000001       true            1                      *
*          )      1000002       false           0                      *
*          union  1000003       intersection   implicit                *
*                                                                      *
*       Lasr modified by K.Niita on 2009/09/30                              *
*                                                                      *
************************************************************************

      dimension lg(*)

*-----------------------------------------------------------------------

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

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      subroutine quart(n,c,r4,jj,mm)
*                                                                      *
*       use the procedure of la-4299 to solve the quartic equation     *
*          x**4+b*x**3+c*x**2+d*x+e=0  (c(i),i=1,5)=1,b,c,d,e          *
*       the jj real roots are returned in r4.  the roots with odd      *
*       multiplicity are listed first, in ascending order.  if mm=0,   *
*       the roots with even multiplicity are omitted.                  *
*                                                                      *
*       Last modified by K.Niita on 2009/09/30                         *
*                                                                      *
************************************************************************

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      dimension c(5),r4(4),tp(16)
      parameter (pf=(2.*3.1415926535898d0)/3.,th=1d0/3.)
      parameter ( epss = 1.0d-10 )

*-----------------------------------------------------------------------

         jj = 0
         if(n.eq.3) goto 90

*-----------------------------------------------------------------------
*     let y=x+b/4 to reduce the quartic equation to
*       y**4+q*y**2+r*y+s=0  (tp(i),i=1,5)=b/2,b/4,q,r,s
*-----------------------------------------------------------------------

         tp(1) = .5*c(2)
         tp(2) = .25*c(2)
         tp(3) = c(3)-1.5*tp(1)**2
         tp(4) = c(4)+tp(1)*(tp(1)**2-c(3))
         tp(5) = c(5)-.0625*tp(1)*(5.*c(4)-c(3)*tp(1)+3.*tp(4))

*-----------------------------------------------------------------------
*     steps 2-7 (la-4299) (tp(i),i=6,11)=b/3,c,d,p/3,q/2,w
*-----------------------------------------------------------------------

         if(abs(tp(4)).le.epss*max(abs(c(4)),tp(1)**2,abs(tp(1)*c(3))))
     &    goto 40
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

*-----------------------------------------------------------------------
*     steps 8-11 (la-4299) double or triple roots.
*-----------------------------------------------------------------------

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

*-----------------------------------------------------------------------
*     step 12 (la-4299)
*-----------------------------------------------------------------------

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

*-----------------------------------------------------------------------
*     steps 13-24 (la-4299) trivial case; r=tp(4)=0
*-----------------------------------------------------------------------

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

*-----------------------------------------------------------------------
*     cubic equation.  x**3+b*x**2+c*x+d=0  (c(i),i=1,4)=1,b,c,d
*     let y=x+b/3.  then  y**3+p*y+q=0   (tp(i),i=1,4)=b/3,p/3,q/2,w
*-----------------------------------------------------------------------

   90    tp(1) = th*c(2)
         tp(2) = th*c(3)-tp(1)**2
         tp(3) = .5*(c(4)-tp(1)*(c(3)-2.*tp(1)**2))
         tp(4) = tp(2)**3+tp(3)**2
         if(tp(4).gt.0.) goto 110
         if(tp(4).eq.0.) goto 100

*-----------------------------------------------------------------------
*     case i (la-4299)  3 distinct roots
*-----------------------------------------------------------------------

         tp(5) = 2.*sqrt(-tp(2))
         tp(6) = 8.*tp(3)/tp(5)**3
         if(abs(tp(6)).ge.1.) goto 100
         tp(7) = th*acos(-tp(6))
         r4(1) = tp(5)*cos(pf+tp(7))-tp(1)
         r4(2) = tp(5)*cos(pf-tp(7))-tp(1)
         r4(3) = tp(5)*cos(tp(7))-tp(1)
         jj = 3
         return

*-----------------------------------------------------------------------
*     case ii (la-4299) 1 distinct and 1 double root
*-----------------------------------------------------------------------

  100    r4(1) = -2.*sign(sqrt(-tp(2)),tp(3))-tp(1)
         jj = 1
         if(mm.eq.0) return
         if(r4(1).eq.0.) return
         jj = 2
         r4(2) = -.5*r4(1)
         return

*-----------------------------------------------------------------------
*     case iii (la-4299)  1 distinct root
*-----------------------------------------------------------------------

  110    tp(5) = -tp(3)-sign(sqrt(tp(4)),tp(3))
         tp(5) = sign(abs(tp(5))**th,tp(5))
         r4(1) = tp(5)-tp(2)/tp(5)-tp(1)
         jj= 1

*-----------------------------------------------------------------------

      return
      end

