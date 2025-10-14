************************************************************************
*                                                                      *
      subroutine gcell(jsn,jsi,dsin,idsi,ill,ilf,
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr) ! 20220906frtati iods deleted
*                                                                      *
*       read [cell] section of GG input files                          *
*       modified by K.Niita on 2009/03/02                              *
*                                                                      *
************************************************************************
      use LATDATAMOD
      use moddas
      use moddas_character

      implicit real*8 (a-h,o-z)

      include 'param.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      character m_err*200
      common /error/ m_err, l_err, k_err

*-----------------------------------------------------------------------

      common /regdm/  idmg(kvlmax)
      common /regda/  ichl(kvlmax), chsm(kvlmax), ichmx, iod
      character       chsm*10
      common /regdc/  idrg(kvlmax), idgr(kvmmax)
      common /regdu/  iuni(kvlmax)
      common /regde/  ichp(kvlmax), ilat(kvlmax), idct(kvlmax)
      common /regdf/  ioc, ifilt, ifil(kvlmax)
      common /regdl/  ilike(kvlmax)
      common /volreg/ dvol(kvlmax)
      common /tmpmsg/ rtmps(kvlmax), mntmp, ntmps(kvlmax)
      common /tmpreg/ dtmp(kvlmax)
      common /regdd/  ivolm, iimpo

      common /paraj/  mstz(300), parz(300)
      common /voxel/  ivoxel

      common /ggcell/ icells, iobo
      character chcfg*100

*-----------------------------------------------------------------------

      common /inpec/  ititl, ipara, ibody, iregn, llarr, itby, itar
      common /inggs/  iog, igcel, ioa, igsuf, iob, igtrs
      common /celda/  deng(kvlmax)
      common /celdb/  idsn(kvlmax), idtn(kvlmax)

      common /paran/  icfn(100), ilfn(100), chfn(100)
      character       chfn*200

*-----------------------------------------------------------------------



      integer,allocatable :: klat(:),jlat(:) !FURUTA20201127


      character chdm*1

*-----------------------------------------------------------------------

      character chin*200, chlw*200, chcm*200

      character dsin(0:9)*200
      dimension idsi(0:9)

      dimension ill(0:9), ilf(0:9)

      character dkam*6

      logical deqn5

      dimension vtrs(13)
      dimension kfil(6)

      dimension lfil(6)

      logical   exex

*-----------------------------------------------------------------------

      character chss(4)*2
      data chss / '( ', ') ',': ','# ' /

*-----------------------------------------------------------------------

      dimension lschn(20), ischn(20)
      character schan(20)*9
      data icsu / 14 /

      data ( schan(i), i = 1, 14 ) /
     &    'u       ','trcl    ','*trcl   ','lat     ','fill    ',
     &    '*fill   ','mat     ','rho     ','vol     ','tmp     ',
     &    'tfile   ','tsfac   ','nfile   ','hfile   '/

      data ( lschn(i), i = 1, 14 ) /
     &     1,         4,         5,         3,         4,
     &     5,         3,         3,         3,         3,
     &     5,         5,         5,         5/
*-----------------------------------------------------------------------
cFURUTA20150714 TETRA !FURUTA20190208
      integer ntfile,nlat0
      logical itsfac,itrcl
      real(8) tsfactor
      integer nlat3,itetcl(10),itetvol,itetauto,itgchk
      integer ltfile(10),itfform(10)
      real(8) tetsfac(10)
      character(200) tfilename(10)
      common /tetf1/ nlat3,itetcl,itetvol,itetauto,itgchk
      common /tetf2/ tetsfac,ltfile,itfform,tfilename
*-----------------------------------------------------------------------
c     FURUTA20181009 compressed voxel
      integer,allocatable :: lattmp(:),lattmp2(:)
*-----------------------------------------------------------------------

            ierr  = 0

            katot = 0 ! T.Sato 2021/09/20 initialization

            if( igcel .ne. 0 ) goto 991

            if( iregn .ne. 0 .or. ibody .ne. 0 ) goto 990

            icelp = 0

            latmax=1                  !FURUTA20201127
            allocate(jlat(1),klat(1)) !FURUTA20201127

cFURUTA20150714 TETRA---------------------------------------------------

            nlat3 = 0
            ntfile= 0
            itetcl(1:10)=0
            tetsfac(1:10)=1.0d0
            nlat0 = 0
            itsfac=.false.
            itrcl=.false.

*-----------------------------------------------------------------------

*-----------------------------------------------------------------------
*     set flag for cell information
*     icells = 0 : Read from input file without echo of [cell]
*            = 1 : Read from binary file 'gcell.bin'    -> ivoxel=1
*            = 2 : Write a binary file from phits input -> ivoxel=2
*            = 3 : default: normal with echo
*-----------------------------------------------------------------------

                  icells = mstz( 142 )

               if( icells .lt. 0 .or. icells .ge. 4 ) then

                     write(ErrCha,'('' Error : icells should be'',
     &               '' 0,  1, 2 or 3.'')')
                     ErrID = 'L:159/R:gcell/F:read03.f' !E00_015_002
                     call ErrWrite(ErrID,ErrCha)

                     goto 999

               end if

*-----------------------------------------------------------------------

            if( icells .eq. 1 ) then

                  jpn = 4

  150          continue
               call readl(jsn,jsi,dsin,idsi,ill,ilf,'$',
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

                  if( ierr .ne. 0 ) return
                  if( jpn  .eq. 3 ) return
                  if( iskip .ne. 0 ) goto 150

               if( i1 .le. 5 .and. chlw(i1:i1) .eq. '[' ) then

*-----------------------------------------------------------------------

                     chcfg = chfn(19)(1:ilfn(19))//'.cfg'
                     ichcfg = ilfn(19) + 4

                     inquire( file = chcfg, exist = exex )

                  if( exex .eqv. .false. ) then
                     write(*,'(/'' Error : icells cfg data file'',
     &               '' does not exist.''/
     &               '' file name = '',200a1)')
     &               ( chcfg(i:i),i=1, ichcfg )
                     goto 999
                  end if

                     ioj = 203
                     open(ioj,status='unknown',file= chcfg,
     &                    form='unformatted')

                        read(ioj) igtrs, iobo
                     do i = 1, igtrs
                        read(ioj) itrs, ( vtrs(j), j = 1, 13 ), idtn(i)
                     end do

                        read(ioj)  igcel, igsuf0, ifilt, nlat3

                     do i = 1, igcel
                        read(ioj)  ichl(i), chsm(i), ichp(i), ilat(i),
     &                             idct(i), ifil(i), ilike(i), idmg(i),
     &                             iuni(i), deng(i), dvol(i), dtmp(i),
     &                             idrg(i), idgr(idrg(i))
                     end do
                     do i = 1, nlat3
                        read(ioj)  itetcl(i), tetsfac(i)
                     end do

                     close(ioj)

*-----------------------------------------------------------------------

                  jpn = 1
                  return

               end if

                  goto 150

*-----------------------------------------------------------------------

            else if( icells .eq. 2 ) then

               iod = 26
               open(iod,status='unknown',file= chfn(19),
     &              form='unformatted')

            end if

*-----------------------------------------------------------------------
*     set flag for filling cell information
*     ivoxel = 0 : Data read from input file of phits (normal)
*            = 1 : Data read from external binary file 'voxel.bin'
*            = 2 : Generate a binary file (iov) from phits input
*            = 3 : 0 + echo of voxel data in phits.out
*     added by daiki on Feb.2008
*     modified by T. Furuta on 20200515
*-----------------------------------------------------------------------

               if( icells .eq. 2 ) then

                     ivoxel = icells

               else

                     ivoxel = mstz( 68 )

               end if

               if( ivoxel .lt. 0 .or. ivoxel .ge. 4 ) then !FURUTA20200515

                  write(ErrCha,'('' Error : ivoxel should be'',
     &            '' 0,  1, 2 or 3.'')')
                  ErrID = 'L:263/R:gcell/F:read03.f' !E00_015_001
                  call ErrWrite(ErrID,ErrCha)
                  goto 999

               end if

*-----------------------------------------------------------------------
*     assign the pointers for klat, jlat and chrg
*-----------------------------------------------------------------------


               klt=0   !FURUTA20201127
               jlt=0   !FURUTA20201127
               mgmax=( ((mmmax-1)*2+1)/2+1 ) !FURUTA20201207

                  if( mgmax .gt. mdas ) goto 978


                  mci  = ( mgmax - 1.0 ) * 8 + 1
                  mcmx = ( mdas  - 1 ) * 8 + 1 - mci
                  call moddas_allocate_cha(MAX_NUM_CHRG, chrg)
                  mci = 0

*-----------------------------------------------------------------------
*     file (iod)
*-----------------------------------------------------------------------

               if( icells .ne. 2 ) then

                  iod = 26
                  open(iod,status='scratch',form='unformatted')

               end if

*-----------------------------------------------------------------------

            ichmx = 0

*-----------------------------------------------------------------------
*     read one line from jsi
*-----------------------------------------------------------------------

  140 continue

            call readl(jsn,jsi,dsin,idsi,ill,ilf,'$',
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

               if( ierr .ne. 0 ) goto 1000
               if( jpn  .eq. 3 ) goto 1000

               if( iskip .ne. 0 ) goto 140

               if( jpn .eq. 4) jpn = 0  !FURUTA20200515 skip infl by ivoxel=1

*-----------------------------------------------------------------------
*     end of cell section
*-----------------------------------------------------------------------

            if( i1 .le. 5 .and. chlw(i1:i1) .eq. '[' ) then

               jpn = 1
               goto 1000

            end if

*-----------------------------------------------------------------------
*     start new cell line
*-----------------------------------------------------------------------

               ic = i1

      if( i1 .le. 5 ) then

*-----------------------------------------------------------------------
*     modified by daiki on Feb.2008
*-----------------------------------------------------------------------

         if( igcel .gt. 0 ) then

               ichl(igcel) = igc
               ichp(igcel) = max(igc,igp)

               ichmx = max(ichmx,ichp(igcel))

               write(iod) (chrg(mci+k:mci+k),k=1,ichp(igcel))

               if( latot .gt. 0 .and. iulat .ne. latot .and.
     &             ivoxel .ne. 1 ) goto 985

*-----------------------------------------------------------------------
               if(nlat3.ne.ntfile)then
                write(ErrCha,'('' Error : LAT=3 and TFILE should'',
     &            '' defined togheter.'')')
                 ErrID = 'L:356/R:gcell/F:read03.f' !E07_002_001
                 call ErrWrite(ErrID,ErrCha)
                goto 999
               endif

cFURUTA20150714 TETRA---------------------------------------------------
               if(nlat3.gt.nlat0)then
                if(itsfac)then
                 tetsfac(nlat3)=tsfactor
                 itsfac=.false.
                endif
                if(itrcl)then

                 ErrCha = ''
                 ErrID = 'L:370/R:gcell/F:read03.f' !E07_003_001
                 call ErrWrite(ErrID,ErrCha)

                 write(*,'(/'' Error : trcl is not compatible'',
     &                '' with LAT=3.'')')
                 write(*,'('' Apply trcl on higher rank universe'',
     &                '' instead.'')')
                  goto 999
                endif
                nlat0=nlat3
               else
                itsfac=.false.
                itrcl=.false.
               endif

*-----------------------------------------------------------------------
*           write fill information on temporary file
*           modified by daiki on Feb.2008
*-----------------------------------------------------------------------

            if( ifill .ne. 0 .or. jfill .ne. 0 ) then

                  if( ifill .eq. 0 ) then

                     ifill = jfill
                     latot = katot

                  end if

                     ifil(igcel) = ifill

                     icelp = 0

*-----------------------------------------------------------------------
               if( ivoxel .ne. 1 ) then
*-----------------------------------------------------------------------

                  if( ifill .eq. 1 ) then

                        write(ioc) klat(klt+1), jlat(jlt+1)

                  else if( ifill .eq. 6 ) then

                   jcount=maxval(jlat(jlt+1:jlt+latot)) !FURUTA20181009
                   if(jcount.gt.0)then                  !FURUTA20181009

                        write(ioc) ( kfil(i), i = 1, 6 ), latot

                     do i = 1, latot

                        write(ioc) klat(klt+i), jlat(jlt+i)

                     end do

                   else                                 !FURUTA20181009
cFURUTA20181009---------------------------------------------------------
                    write(ioc) ( kfil(i), i = 1, 6 ), -latot
                    allocate(lattmp(latot),lattmp2(latot))
                    do i=1,latot
                     lattmp(i)=klat(klt+i)
                    enddo
                    call latencode(latot,lattmp,latot2,lattmp2)
cFURUTA20200515---------------------------------------------------------

                    if( ivoxel .eq. 0 .or. ivoxel .eq. 3 ) then

                       if(nlatind.eq.0)then
                          call ALLOCATE_ldata(latot2+1)
                       else
                          call RESIZE_ldata(latot2+1)
                        end if

                     ilatind=ilatind+1
                     ldata(ilatind)=latot2
                     ldata(ilatind+1:ilatind+latot2)=lattmp2(1:latot2)
                     ilatind=ilatind+latot2

                    else if( ivoxel .eq. 2 ) then

                       write(ioc)latot2
                       write(ioc)lattmp2(1:latot2)

                    endif

*-----------------------------------------------------------------------

                    deallocate(lattmp,lattmp2)

                   endif                                !FURUTA20181009
                  end if

*-----------------------------------------------------------------------
               end if
*-----------------------------------------------------------------------

            end if

         end if

*-----------------------------------------------------------------------

               if( icelp .ne. 0 ) goto 983

               ifill = 0
               latot = 0
               iulat = 0
               ifat  = 0
               itprs = 0

               igcel = igcel + 1
               if(igcel.gt.kvlmax) goto 975
               igc   = 0
               igp   = 0
               igpn  = 0

               jlike = 0
               klike = 0
               jfill = 0
               ilike(igcel) = 0

*-----------------------------------------------------------------------
*           read cell number
*-----------------------------------------------------------------------

               call snum(chlw,ic,i3,ic2,cvvv,ierr)
               if( ierr .ne. 0 ) goto 999

                  idrg(igcel) = nint( cvvv )

                  if( idrg(igcel) .le. 0 .or.
     &                idrg(igcel) .ge. kvmmax ) goto 994

                  if( idgr(idrg(igcel)) .ne. 0 ) goto 993

                  idgr(idrg(igcel)) = igcel

                  ic = jnumc(chlw,ic2,i3)

*-----------------------------------------------------------------------
*           like but
*-----------------------------------------------------------------------

            if( chlw(ic:ic+3) .eq. 'like' ) then

                  ic = ic + 4
                  ic = jnumc(chlw,ic,i3)

               call snum(chlw,ic,i3,ic2,cvvv,ierr)
               if( ierr .ne. 0 ) goto 999

                  jlike = nint( cvvv )

               if( jlike .le. 0 .or. jlike .ge. kvmmax ) goto 982
               if( idgr(jlike) .eq. 0 ) goto 982

                  ic = jnumc(chlw,ic2,i3)

               if( chlw(ic:ic+2) .ne. 'but' ) goto 981

                  ic = ic + 3

                  ilike(igcel) = jlike
                  klike = idgr(jlike)

            end if

*-----------------------------------------------------------------------
*           read material number
*-----------------------------------------------------------------------

            if( klike .eq. 0 ) then

               call snum(chlw,ic,i3,ic2,cvvv,ierr)
               if( ierr .ne. 0 ) goto 999

                  idmg(igcel) = nint( cvvv )

                  if( idmg(igcel) .ge. kvmmax ) goto 992
                  if( idmg(igcel) .lt. -1 ) goto 992

                  ic = jnumc(chlw,ic2,i3)

*-----------------------------------------------------------------------
*           read particle density
*-----------------------------------------------------------------------

                     deng(igcel) = 0.0d0

               if( idmg(igcel) .gt. 0 ) then

                  call snum(chlw,ic2,i3,ic2,cvvv,ierr)
                  if( ierr .ne. 0 ) goto 999
                  if( cvvv .eq. 0.0d0) goto 976 !FURUTA 2014/02/04

                     deng(igcel) = cvvv

               end if

                  ic = jnumc(chlw,ic2,i3)

            end if

      end if

*-----------------------------------------------------------------------
*        rest part or sequential line
*-----------------------------------------------------------------------

            if( igcel .eq. 0 ) goto 999

               img = 0

            ic = ic - 1

  100       ic = ic + 1

            if( ic .gt. i3 ) goto 140

               k = index('():#',chlw(ic:ic))

! Nais_2024 >>>
               k_all = 0
               if(k == 4) then
                  if(ic + 3 < i3) then
                    k_all = index(chlw(ic:i3),'#all ')
                  else
                    k_all = index(chlw(ic:i3),'#all')
                  end if
!                  if(k_all == 1) then
!                    write(*,*) 'read03.f,chlw(ic:ic+3)=',chlw(ic:ic+3)
!                  endif
               endif
! Nais_2024 <<<

               if( deqn5( chlw(ic:ic) ) ) k = 5

            if( klike .ne. 0 .and. igpn .eq. 0 .and. icelp .eq. 0 .and.
     &          k .ge. 1 .and. k .le. 5 ) goto 980

            if( igpn .ne. 0 .and. icelp .eq. 0 .and.
     &          k .ge. 1 .and. k .le. 5 ) goto 999

            if( icelp .ne. 0 .and.
     &          icelp .ne. 5 .and. icelp .ne. 6 ) goto 300

            if( icelp .eq. 5 .or. icelp .eq. 6 ) then

               if( k .eq. 1 .or. k .eq. 2 .or.
     &             k .eq. 3 .or. k .eq. 5 ) goto 300

               if( k .eq. 4 ) goto 999
               if( k .eq. 0 .and. chlw(ic:ic) .eq. ' ' ) goto 100

               icelp = 0

               goto 250

            end if

*-----------------------------------------------------------------------
*           cell definitions
*-----------------------------------------------------------------------

               if( k .eq. 5 .and. img .eq. 0 ) then

                     img = 1
                     ini = ic

               end if

               if( img .eq. 1 .and.
     &           ( ( k .eq. 5 .and. ic .eq. i3 ) .or.
     &             ( k .ne. 5 ) ) ) then

                     img = 0

                  if( k .eq. 5 ) then

                     ifi = ic

                  else

                     ifi = ic - 1

                  end if

                     call onum(chlw,ini,ifi,cvvv,ierr)
                     if( ierr .ne. 0 ) goto 999


cFURUTA20131226 optimization bug? in gfortran 4.8 can be fixed
cFURUTA20140213 but it does not work in BX900
                  do i = ini, ifi

                     j = igc + i - ini + 1

                     chrg(mci+j:mci+j) = chin(i:i)

                  end do

                     igc = igc + ifi - ini + 2

                     chrg(mci+igc:mci+igc) = ' '

               end if

               if( k .ge. 1 .and. k .le. 3 ) then

                     chrg(mci+igc+1:mci+igc+2) = chss(k)

                     igc = igc + 2

               else if( k .eq. 4 ) then

                     chrg(mci+igc+1:mci+igc+1) = chss(k)(1:1)

                     igc = igc + 1

! Nais_2024 >>>
                  if(k_all == 1) then
!            write(*,*) 'read03.f,ic,i3,k,k_all,igc,chss(k)(1:1)=',
!     &                           ic,i3,k,k_all,igc,chss(k)(1:1)
                       if(ic + 3 <  i3) then
                         chrg(mci+igc+1:mci+igc+4) = 'all '
                         igc = igc + 4
                         ic = ic + 4
                       else
                         chrg(mci+igc+1:mci+igc+3) = 'all'
                         igc = igc + 3
                         ic = ic + 3
                       endif

                     end if
! Nais_2024 <<<

               end if

               if( k .ne. 0 .or.
     &           ( k .eq. 0 .and. chlw(ic:ic) .eq. ' ' ) ) goto 100

*-----------------------------------------------------------------------
*        re-read cell definition and parameters from like-but cell
*-----------------------------------------------------------------------

         if( klike .ne. 0 .and. igpn .eq. 0 ) then

                  idmg(igcel) = idmg(klike)
                  deng(igcel) = deng(klike)
                  iuni(igcel) = iuni(klike)
                  ilat(igcel) = ilat(klike)
                  idct(igcel) = idct(klike)

*-----------------------------------------------------------------------
*              read cell information for like-but
*-----------------------------------------------------------------------

                     rewind iod

                     igc = ichl(klike)

                  do k = 1, igcel - 1

                     if( k .eq. klike ) then

                        read(iod) (chrg(mci+i:mci+i),i=1,ichl(k))

                     else

                        read(iod) (chdm,i=1,ichl(k))

                     end if

                  end do

*-----------------------------------------------------------------------
*           read fill information for like-but
*-----------------------------------------------------------------------

            if( ifil(klike) .ne. 0 ) then

                     rewind ioc
                  if( ivoxel .eq. 0 .or.
     &                ivoxel .eq. 3 )  call REWIND_ldata

*-----------------------------------------------------------------------

               do k = 1, igcel - 1

*-----------------------------------------------------------------------

               if( ifil(k) .eq. 1 ) then

                           if( k .eq. klike ) then

                              read(ioc)  klat(klt+1), jlat(jlt+1)

                           else

                              read(ioc) idmm, idmm

                           end if

               else if( ifil(k) .eq. 6 ) then

                           if( k .eq. klike ) then

                              read(ioc) ( kfil(i), i = 1, 6 ), katot

                              iatot = katot

                           else

                              read(ioc) ( lfil(i), i = 1, 6 ), iatot

                           end if

*-----------------------------------------------------------------------
                  if( iatot .gt. 0 ) then       !FURUTA20181009
*-----------------------------------------------------------------------

                        do i = 1, iatot

                           if( k .eq. klike ) then

                              read(ioc)  klat(klt+i), jlat(jlt+i)

                           else

                              read(ioc) idmm, idmm

                           end if

                        end do

*-----------------------------------------------------------------------
                  else                          !FURUTA20181009
*-----------------------------------------------------------------------


cFURUTA20181009---------------------------------------------------------
                     iatot = -iatot
                     allocate(lattmp(iatot))
cFURUTA20200515---------------------------------------------------------

                  if( ivoxel .eq. 0 .or. ivoxel .eq. 3 ) then

                      ilatind=ilatind+1
                      iatot2=ldata(ilatind)
                      allocate(lattmp2(iatot2))
                      lattmp2(1:iatot2)=ldata(ilatind+1:ilatind+iatot2)
                      ilatind=ilatind+iatot2

*-----------------------------------------------------------------------

                  else

                      read(ioc)iatot2
                      allocate(lattmp2(iatot2))
                      read(ioc)lattmp2(1:iatot2)

                  end if

*-----------------------------------------------------------------------

                     call latdecode(iatot2,lattmp2,iatot,lattmp,
     &                    icount,ierr)
                     do i=1,iatot
                      klat(klt+i)=lattmp(i)
                     enddo
                     do i=1,iatot
                      jlat(jlt+i)=0
                     enddo

                     deallocate(lattmp,lattmp2)

                     if(ierr.ne.0)then
                      write(io,'(/''** ERROR : in [cell] section, '',
     &                   ''compressed lattice element is inconsistent.''/
     &                     ''   latot ='',i9,/
     &                     ''   recoded number='',i9,/)')  iatot,icount
                      goto 999
                     end if

*-----------------------------------------------------------------------
                  end if              !FURUTA20181009
*-----------------------------------------------------------------------

               end if
               end do

                           jfill = ifil(klike)

            end if

         end if

*-----------------------------------------------------------------------
*        start cell parameter part
*-----------------------------------------------------------------------

  250       continue

            do i = 1, icsu

               il = ic + lschn(i) - 1

               if( chlw(ic:il) .eq. schan(i)(1:lschn(i)) ) goto 200

            end do

               goto 999

  200       continue

*-----------------------------------------------------------------------

               if( itprs .ne. 0 ) goto 999

                  icelp = i
                  igpn  = igpn + 1
                  ieqc  = 0
                  itprs = 0
                  igkst = 0

               if( igpn .eq. 1 ) then

                  igp = igc

                  do i = 1, icsu

                     ischn(i) = 0

                  end do

               end if

               ischn( icelp ) = ischn( icelp ) + 1
               if( ischn(icelp) .gt. 1 ) goto 987

*-----------------------------------------------------------------------

  300    continue

            if( ieqc .eq. 0 ) then

               ic = inumc(chlw,ic,i3,'=')
               if( ic .gt. i3 ) goto 983

               ieqc = 1
               ic = jnumc(chlw,ic+1,i3)
               if( ic .gt. i3 )  goto 140

            end if

*-----------------------------------------------------------------------
*        u : universe
*-----------------------------------------------------------------------

         if( icelp .eq. 1 ) then

               call snum(chlw,ic,i3,icl,cvvv,ierr)
               if( ierr .ne. 0 ) goto 999

               iuni(igcel) = nint( cvvv )

               chrg(mci+igp+1:mci+igp+3) = ' u='
               igp = igp + 3
               chrg(mci+igp+1:mci+igp+icl-ic) = chlw(ic:icl-1)
               igp = igp + icl - ic

               ic = jnumc(chlw,icl,i3) - 1

               icelp = 0

*-----------------------------------------------------------------------
*        lat : lattice
*-----------------------------------------------------------------------

         else if( icelp .eq. 4 ) then

               call snum(chlw,ic,i3,icl,cvvv,ierr)
               if( ierr .ne. 0 ) goto 999

               ilat(igcel) = nint( cvvv )

               if( ilat(igcel) .ne. 1 .and.
     &             ilat(igcel) .ne. 2 .and.
     &             ilat(igcel) .ne. 3) goto 986 !FURUTA20150714 TETRA

               chrg(mci+igp+1:mci+igp+5) = ' lat='
               igp = igp + 5
               chrg(mci+igp+1:mci+igp+icl-ic) = chlw(ic:icl-1)
               igp = igp + icl - ic

               ic = jnumc(chlw,icl,i3) - 1

               icelp = 0

               if( ilat(igcel) .eq. 3)then !FURUTA20150714 TETRA
                  nlat3 = nlat3 + 1
                  itetcl(nlat3)=igcel
               endif

*-----------------------------------------------------------------------
*        trcl, *trcl : cell transform
*-----------------------------------------------------------------------

         else if( icelp .eq. 2 .or. icelp .eq. 3 ) then

               if( icelp .eq. 2 ) then

                  itrs = 0

               else if( icelp .eq. 3 ) then

                  itrs = 1

               end if

               itrcl=.true. !FURUTA20150714 TETRA

*-----------------------------------------------------------------------

            if( chlw(ic:ic) .eq. '(' ) then

                  itprs = itprs + 1
                  if( itprs .gt. 1 ) goto 999

                  vtrs( 1) = 0.0
                  vtrs( 2) = 0.0
                  vtrs( 3) = 0.0
                  vtrs( 4) = 1.0
                  vtrs( 5) = 0.0
                  vtrs( 6) = 0.0
                  vtrs( 7) = 0.0
                  vtrs( 8) = 1.0
                  vtrs( 9) = 0.0
                  vtrs(10) = 0.0
                  vtrs(11) = 0.0
                  vtrs(12) = 1.0
                  vtrs(13) = 1.0

               if( icelp .eq. 3 ) then

                  vtrs( 4) =  0.0
                  vtrs( 5) = 90.0
                  vtrs( 6) = 90.0
                  vtrs( 7) = 90.0
                  vtrs( 8) =  0.0
                  vtrs( 9) = 90.0
                  vtrs(10) = 90.0
                  vtrs(11) = 90.0
                  vtrs(12) =  0.0

               end if

               if( icelp .eq. 2 ) then

                  chrg(mci+igp+1:mci+igp+7) = ' trcl=('
                  igp = igp + 7

               else if( icelp .eq. 3 ) then

                  chrg(mci+igp+1:mci+igp+8) = ' *trcl=('
                  igp = igp + 8

               end if


                  ic = jnumc(chlw,ic+1,i3)
                  if( ic .gt. i3 )  goto 140

            end if

*-----------------------------------------------------------------------

            if( itprs .gt. 0 ) then

  400          continue

               if( chlw(ic:ic) .eq. ')' .and. igkst .eq. 0 ) then

                  goto 999

               else if( chlw(ic:ic) .eq. ')' .and. igkst .eq. 1 ) then

                  idct(igcel) = nint( vtrs(1) )

                  chrg(mci+igp+1:mci+igp+2) = ' )'
                  igp = igp + 2

                  ic = jnumc(chlw,ic+1,i3) - 1

                  itprs = 0
                  igkst = 0
                  icelp = 0

               else if( chlw(ic:ic) .eq. ')' .and. igkst .gt. 1 ) then

                  if( igtrs .eq. 0 ) then
                     iobo = 1
                     iob = 71
                     open(iob,status='scratch',form='unformatted')

                  end if

                  igtrs = igtrs + 1

                  idtn(igtrs) = 1000000 + igtrs
                  idct(igcel) = 1000000 + igtrs

                  write(iob) itrs, ( vtrs(i), i = 1, 13 )

                  chrg(mci+igp+1:mci+igp+2) = ' )'
                  igp = igp + 2

                  ic = jnumc(chlw,ic+1,i3) - 1

                  itprs = 0
                  igkst = 0
                  icelp = 0

               else

                  if( ic .gt. i3 ) goto 140

                  call unum(chlw,ic,i3,icl,cvvv,ierr)
                  if( ierr .ne. 0 ) goto 999

                  igkst = igkst + 1

                  if( igkst .gt. 13 ) goto 984

                  vtrs(igkst) = cvvv

                  chrg(mci+igp+1:mci+igp+icl-ic+1) = ' '//chlw(ic:icl-1)
                  igp = igp + icl - ic + 1

                  ic = jnumc(chlw,icl,i3)

                  goto 400

               end if

*-----------------------------------------------------------------------
*              trcl = 3 type
*-----------------------------------------------------------------------

            else if( itprs .eq. 0 ) then

                  call snum(chlw,ic,i3,icl,cvvv,ierr)
                  if( ierr .ne. 0 ) goto 999

                  idct(igcel) = nint( cvvv )

               if( icelp .eq. 2 ) then

                  chrg(mci+igp+1:mci+igp+6) = ' trcl='
                  igp = igp + 6

               else if( icelp .eq. 3 ) then

                  chrg(mci+igp+1:mci+igp+7) = ' *trcl='
                  igp = igp + 7

               end if

                  chrg(mci+igp+1:mci+igp+icl-ic) = chlw(ic:icl-1)
                  igp = igp + icl - ic

                  ic = jnumc(chlw,icl,i3) - 1

                  icelp = 0

            end if

*-----------------------------------------------------------------------
*        fill, *fill : cell-filling universes, with transformations
*-----------------------------------------------------------------------

         else if( icelp .eq. 5 .or. icelp .eq. 6 ) then

               if( ifill .eq. 0 .and. itprs .ne. 0 ) goto 999

               if( icelp .eq. 5 ) then

                  itrs = 0

               else if( icelp .eq. 6 ) then

                  itrs = 1

               end if

*-----------------------------------------------------------------------
*           modified by daiki on Feb.2008
*-----------------------------------------------------------------------

            if( ifill .eq. 0 ) then

                  call unum(chlw,ic,i3,icl,cvvv,ierr)
                  if( ierr .ne. 0 ) goto 999

                  kfil(1) = nint( cvvv )
                  klat(klt+1) = nint( cvvv )
                  jlat(jlt+1) = 0

                  ifill = ifill + 1

                  if( ifilt .eq. 0 ) then

                     ioc = 66 ! T.Sato 2018/09/10, avoid overlap with multiplier

*-----------------------------------------------------------------------
                     if( ivoxel .eq. 0 .or. ivoxel .eq. 3 ) then

                        open(ioc,status='scratch',form='unformatted')

*-----------------------------------------------------------------------
cKN 2020/06/23  for ivoxel = 1, 2

                     else if( ivoxel .eq. 1 .or. ivoxel .eq. 2 ) then

*-----------------------------------------------------------------------

                        if( ivoxel .eq. 1 ) then
                           inquire( file = chfn(18), exist = exex )

                           if( exex .eqv. .false. ) then

                              write(*,'(/'' Error : voxel data file'',
     &                        '' does not exist.''/
     &                        '' file name = '',200a1)')
     &                        ( chfn(18)(i:i),i=1, ilfn(18) )

                              goto 999

                           end if
                        end if
*-----------------------------------------------------------------------

                        open(ioc,status='unknown',file= chfn(18),
     &                       form='unformatted')

                     end if

*-----------------------------------------------------------------------

                     call INIT_ldata !FURUTA20200515

                  end if

                  ifilt = ifilt + 1

                  if( icelp .eq. 5 ) then

                     chrg(mci+igp+1:mci+igp+6) = ' fill='
                     igp = igp + 6

                  else if( icelp .eq. 6 ) then

                     chrg(mci+igp+1:mci+igp+7) = ' *fill='
                     igp = igp + 7

                  end if

                  chrg(mci+igp+1:mci+igp+icl-ic) = chlw(ic:icl-1)
                  igp = igp + icl - ic
                  ic = jnumc(chlw,icl,i3) - 1

                  goto 100

            end if

            if( chlw(ic:ic) .eq. ':' ) then

                  ifat = ifat + 1
                  if( ifat .gt. 3 ) goto 999
                  if( ifill .eq. 0 .or. ifill .ge. 6  ) goto 999
                  chrg(mci+igp+1:mci+igp+1) = ':'
                  igp = igp + 1

                  goto 100

            end if

            if( k .eq. 5 .and. ifat .gt. 0 .and.
     &          ifill .ge. 1 .and. ifill .le. 5 ) then

                  ifill = ifill + 1

                  call unum(chlw,ic,i3,icl,cvvv,ierr)
                  if( ierr .ne. 0 ) goto 999
                  kfil(ifill) = nint( cvvv )

               if( ifill .eq. 6 ) then

                  if( kfil(1) .gt. kfil(2) .or.
     &                kfil(3) .gt. kfil(4) .or.
     &                kfil(5) .gt. kfil(6) ) goto 999

                  latot = ( kfil(2) - kfil(1) + 1 ) *
     &                    ( kfil(4) - kfil(3) + 1 ) *
     &                    ( kfil(6) - kfil(5) + 1 )

cFURUTA20201127---------------------------------------------------------
                  if( latot .gt. latmax )then
                   deallocate(jlat,klat)
                   allocate(jlat(latot),klat(latot))
                   latmax=latot
                  endif
c-----------------------------------------------------------------------
!                  if( latot .gt. mdas ) goto 977 !FURUTA20201127

*-----------------------------------------------------------------------
cFURUTA20200515 skip infl by ivoxel=1
*-----------------------------------------------------------------------

                  if( ivoxel .eq. 1 ) jpn = 4

*-----------------------------------------------------------------------

               end if

               if( ifill .eq. 3 .or. ifill .eq. 5 ) then

                  chrg(mci+igp+1:mci+igp+icl-ic+1) = ' '//chlw(ic:icl-1)
                  igp = igp + icl - ic + 1

               else

                  chrg(mci+igp+1:mci+igp+icl-ic) = chlw(ic:icl-1)
                  igp = igp + icl - ic

               end if

                  ic = jnumc(chlw,icl,i3) - 1

                  goto 100

            end if

            if( k .eq. 5 .and. ifat .gt. 0 .and.
     &          ifill .eq. 6 .and. itprs .eq. 0 .and.
     &          iulat .lt. latot ) then

                  call unum(chlw,ic,i3,icl,cvvv,ierr)
                  if( ierr .ne. 0 ) goto 999

cFURUTA20181009---------------------------------------------------------
                  itmp=nint(cvvv)
                  if(itmp.gt.0)then
                   iulat=iulat+1
                   jlat(jlt+iulat)=0
                   klat(klt+iulat)=itmp
                   klatold=itmp
                  else
                   do i=1,-itmp
                    klat(klt+iulat+i)=klatold
                   enddo
                   do i=1,-itmp
                    jlat(jlt+iulat+i)=0
                   enddo
                   iulat=iulat-itmp
                  endif

*-----------------------------------------------------------------------
cFURUTA20200515 not to show voxels
*-----------------------------------------------------------------------
               if( ivoxel .eq. 3 ) then
                  chrg(mci+igp+1:mci+igp+icl-ic+1) = ' '//chlw(ic:icl-1)
                  igp = igp + icl - ic + 1
               end if

*-----------------------------------------------------------------------

                  ic = jnumc(chlw,icl,i3) - 1

                  goto 100

            end if

*-----------------------------------------------------------------------

            if( chlw(ic:ic) .eq. '(' ) then

                  if( ifill .eq. 0 ) goto 999
                  if( ifat .eq. 3 .and. ifill .ne. 6 ) goto 999
                  if( ifat .eq. 0 .and. ifill .ne. 1 ) goto 999
                  if( ifat .eq. 1 .or. ifat .eq. 2 .or. ifat .gt. 3 )
     &                                                 goto 999

                  itprs = itprs + 1
                  if( itprs .gt. 1 ) goto 999

                  vtrs( 1) = 0.0
                  vtrs( 2) = 0.0
                  vtrs( 3) = 0.0
                  vtrs( 4) = 1.0
                  vtrs( 5) = 0.0
                  vtrs( 6) = 0.0
                  vtrs( 7) = 0.0
                  vtrs( 8) = 1.0
                  vtrs( 9) = 0.0
                  vtrs(10) = 0.0
                  vtrs(11) = 0.0
                  vtrs(12) = 1.0
                  vtrs(13) = 1.0

               if( icelp .eq. 6 ) then

                  vtrs( 4) =  0.0
                  vtrs( 5) = 90.0
                  vtrs( 6) = 90.0
                  vtrs( 7) = 90.0
                  vtrs( 8) =  0.0
                  vtrs( 9) = 90.0
                  vtrs(10) = 90.0
                  vtrs(11) = 90.0
                  vtrs(12) =  0.0

               end if

               if( ivoxel .eq. 3 ) then
                  chrg(mci+igp+1:mci+igp+2) = ' ('
                  igp = igp + 2
               end if

                  ic = jnumc(chlw,ic+1,i3)
                  if( ic .gt. i3 )  goto 140

            end if

*-----------------------------------------------------------------------

            if( itprs .gt. 0 ) then

  500          continue

               if( chlw(ic:ic) .eq. ')' .and. igkst .eq. 0 ) then

                  goto 999

               else if( chlw(ic:ic) .eq. ')' .and. igkst .eq. 1 ) then

                     iult = iulat
                     if( ifat .eq. 0 ) iult = 1

                     jlat(jlt+iult) = nint( vtrs(1) )

                  if( ivoxel .eq. 3 ) then
                     chrg(mci+igp+1:mci+igp+2) = ' )'
                     igp = igp + 2
                  end if

                     ic = jnumc(chlw,ic+1,i3) - 1

                     itprs = 0
                     igkst = 0
                     if( ifat .eq. 0 ) icelp = 0

               else if( chlw(ic:ic) .eq. ')' .and. igkst .gt. 1 ) then

                  if( igtrs .eq. 0 ) then
                     iobo = 1
                     iob = 71
                     open(iob,status='scratch',form='unformatted')

                  end if

                     igtrs = igtrs + 1

                     idtn(igtrs) = 1000000 + igtrs

                     iult = iulat
                     if( ifat .eq. 0 ) iult = 1

                     jlat(jlt+iult) = 1000000 + igtrs

                     write(iob) itrs, ( vtrs(i), i = 1, 13 )

                  if( ivoxel .eq. 3 ) then
                     chrg(mci+igp+1:mci+igp+2) = ' )'
                     igp = igp + 2
                  end if

                     ic = jnumc(chlw,ic+1,i3) - 1

                     itprs = 0
                     igkst = 0
                     if( ifat .eq. 0 ) icelp = 0

               else

                     if( ic .gt. i3 ) goto 140

                     call unum(chlw,ic,i3,icl,cvvv,ierr)

                     if( ierr .ne. 0 ) goto 999

                     igkst = igkst + 1

                     if( igkst .gt. 13 ) goto 984

                     vtrs(igkst) = cvvv

                  if( ivoxel .eq. 3 ) then
                     chrg(mci+igp+1:mci+igp+icl-ic+1) =
     &                                          ' '//chlw(ic:icl-1)
                     igp = igp + icl - ic + 1
                  end if

                     ic = jnumc(chlw,icl,i3)

                     goto 500

               end if

                     goto 100

            end if

*-----------------------------------------------------------------------
*        mat : material
*-----------------------------------------------------------------------

         else if( icelp .eq. 7 ) then

               if( klike .eq. 0 ) goto 979

               call snum(chlw,ic,i3,icl,cvvv,ierr)
               if( ierr .ne. 0 ) goto 999

                  idmg(igcel) = nint( cvvv )

                  if( idmg(igcel) .ge. kvmmax ) goto 992
                  if( idmg(igcel) .lt. -1 ) goto 992

               chrg(mci+igp+1:mci+igp+5) = ' mat='
               igp = igp + 5
               chrg(mci+igp+1:mci+igp+icl-ic) = chlw(ic:icl-1)
               igp = igp + icl - ic

               ic = jnumc(chlw,icl,i3) - 1

               icelp = 0

*-----------------------------------------------------------------------
*        rho : density
*-----------------------------------------------------------------------

         else if( icelp .eq. 8 ) then

               if( klike .eq. 0 ) goto 979

               call snum(chlw,ic,i3,icl,cvvv,ierr)
               if( ierr .ne. 0 ) goto 999

                  deng(igcel) = cvvv

               chrg(mci+igp+1:mci+igp+5) = ' rho='
               igp = igp + 5
               chrg(mci+igp+1:mci+igp+icl-ic) = chlw(ic:icl-1)
               igp = igp + icl - ic

               ic = jnumc(chlw,icl,i3) - 1

               icelp = 0

*-----------------------------------------------------------------------
*        vol : volume
*-----------------------------------------------------------------------

         else if( icelp .eq. 9 ) then

               call snum(chlw,ic,i3,icl,cvvv,ierr)
               if( ierr .ne. 0 ) goto 999

                  dvol(igcel) = cvvv

                  ivolm = 1

               chrg(mci+igp+1:mci+igp+5) = ' vol='
               igp = igp + 5
               chrg(mci+igp+1:mci+igp+icl-ic) = chlw(ic:icl-1)
               igp = igp + icl - ic

               ic = jnumc(chlw,icl,i3) - 1

               icelp = 0

*-----------------------------------------------------------------------
*        tmp : temperature
*-----------------------------------------------------------------------

         else if( icelp .eq. 10 ) then

               call snum(chlw,ic,i3,icl,cvvv,ierr)
               if( ierr .ne. 0 ) goto 999

                  dtmp(igcel) = cvvv

               chrg(mci+igp+1:mci+igp+5) = ' tmp='
               igp = igp + 5
               chrg(mci+igp+1:mci+igp+icl-ic) = chlw(ic:icl-1)
               igp = igp + icl - ic

               ic = jnumc(chlw,icl,i3) - 1

               icelp = 0

*-----------------------------------------------------------------------
*        tfile : tetrahedron file (FURUTA20150714 TETRA)
*-----------------------------------------------------------------------

         else if( icelp .eq. 11 ) then

               icl=inumc(chlw,ic,i3,' ')
               icl=min(icl,inumc(chlw,ic,i3,char(9)))

               chrg(mci+igp+1:mci+igp+7) = ' tfile='
               igp = igp + 7
               chrg(mci+igp+1:mci+igp+icl-ic) = chin(ic:icl-1)
               ! Capital letters in filename are allowed
               igp = igp +icl - ic

               ntfile = ntfile + 1
               ltfile(ntfile) = icl-ic
               tfilename(ntfile)(1:icl-ic)=chin(ic:icl-1)
               itfform(ntfile)=0 !FURUTA20190208

               ic = jnumc(chlw,icl,i3) - 1

               icelp=0

*-----------------------------------------------------------------------
*        tsfac : tetrahedron scaling factor (FURUTA20150714 TETRA)
*-----------------------------------------------------------------------

         else if( icelp .eq. 12 ) then

               call snum(chlw,ic,i3,icl,cvvv,ierr)
               if( ierr .ne. 0 ) goto 999

               tsfactor=cvvv

               chrg(mci+igp+1:mci+igp+7) = ' tsfac='
               igp = igp + 7
               chrg(mci+igp+1:mci+igp+icl-ic) = chlw(ic:icl-1)
               igp = igp + icl - ic

               itsfac=.true.

               ic = jnumc(chlw,icl,i3) - 1

               icelp = 0

*-----------------------------------------------------------------------
*        nfile : NASTRAN file (FURUTA20190208 TETRA)
*-----------------------------------------------------------------------

         else if( icelp .eq. 13 ) then

               icl=inumc(chlw,ic,i3,' ')
               icl=min(icl,inumc(chlw,ic,i3,char(9)))

               chrg(mci+igp+1:mci+igp+7) = ' nfile='
               igp = igp + 7
               chrg(mci+igp+1:mci+igp+icl-ic) = chin(ic:icl-1)
               ! Capital letters in filename are allowed
               igp = igp +icl - ic

               ntfile = ntfile + 1
               ltfile(ntfile) = icl-ic
               tfilename(ntfile)(1:icl-ic)=chin(ic:icl-1)
               itfform(ntfile) = 1

               ic = jnumc(chlw,icl,i3) - 1

               icelp=0

*-----------------------------------------------------------------------
*        hfile : HDF5 file (FURUTA20240709 TETRA)
*-----------------------------------------------------------------------

         else if( icelp .eq. 14 ) then

               icl=inumc(chlw,ic,i3,' ')
               icl=min(icl,inumc(chlw,ic,i3,char(9)))

               chrg(mci+igp+1:mci+igp+7) = ' hfile='
               igp = igp + 7
               chrg(mci+igp+1:mci+igp+icl-ic) = chin(ic:icl-1)
               ! Capital letters in filename are allowed
               igp = igp +icl - ic

               ntfile = ntfile + 1
               ltfile(ntfile) = icl-ic
               tfilename(ntfile)(1:icl-ic)=chin(ic:icl-1)
               itfform(ntfile) = 2

               ic = jnumc(chlw,icl,i3) - 1

               icelp=0

*-----------------------------------------------------------------------

         end if

*-----------------------------------------------------------------------

      goto 100

*-----------------------------------------------------------------------
*     end of this section
*     modified by daiki on Feb.2008
*-----------------------------------------------------------------------

 1000 continue

               ichl(igcel) = igc
               ichp(igcel) = max(igc,igp)

               ichmx = max(ichmx,ichp(igcel))

               write(iod) (chrg(mci+k:mci+k),k=1,ichp(igcel))

               if( latot .gt. 0 .and. iulat .ne. latot .and.
     &             ivoxel .ne. 1 ) goto 985

*-----------------------------------------------------------------------
*           write fill information on temporary file
*           modified by daiki on Feb.2008
*-----------------------------------------------------------------------

            if( ifill .ne. 0 .or. jfill .ne. 0 ) then

                  if( ifill .eq. 0 ) then

                     ifill = jfill
                     latot = katot

                  end if

                     ifil(igcel) = ifill

               if( ivoxel .ne. 1 ) then

                  if( ifill .eq. 1 ) then

                        write(ioc) klat(klt+1), jlat(jlt+1)

                  else if( ifill .eq. 6 ) then

                   jcount=maxval(jlat(jlt+1:jlt+latot)) !FURUTA20181009
                   if(jcount.gt.0)then                  !FURUTA20181009

                        write(ioc) ( kfil(i), i = 1, 6 ), latot

                     do i = 1, latot

                        write(ioc) klat(klt+i), jlat(jlt+i)

                     end do

                   else                                 !FURUTA20181009
cFURUTA20181009---------------------------------------------------------
                    write(ioc) ( kfil(i), i = 1, 6 ), -latot
                    allocate(lattmp(latot),lattmp2(latot))
                    do i=1,latot
                     lattmp(i)=klat(klt+i)
                    enddo
                    call latencode(latot,lattmp,latot2,lattmp2)
cFURUTA20200515---------------------------------------------------------
                    if(ivoxel.eq.0.or.ivoxel.eq.3)then
                     if(nlatind.eq.0)then
                      call ALLOCATE_ldata(latot2+1)
                     else
                      call RESIZE_ldata(latot2+1)
                     endif
                     ilatind=ilatind+1
                     ldata(ilatind)=latot2
                     ldata(ilatind+1:ilatind+latot2)=lattmp2(1:latot2)
                     ilatind=ilatind+latot2
*-----------------------------------------------------------------------
                    else
                     write(ioc)latot2
                     write(ioc)lattmp2(1:latot2)
                    endif
                    deallocate(lattmp,lattmp2)
*-----------------------------------------------------------------------
                   endif                                !FURUTA20181009

                  end if

               end if

            end if

*-----------------------------------------------------------------------

            if( ivoxel .eq. 2 .and.
     &          icells .ne. 1 .and. icells .ne. 2 ) then

                write(*,*) 'Binary file was successfully generated',
     &          ' by ivoxel = 2'
                stop

            end if

*-----------------------------------------------------------------------
            deallocate(jlat,klat) !FURUTA20201127

         if( ichmx > MAX_NUM_CHRG ) then
            write(ErrCha,'(a,a,i5,a,a,i5,a)')
     &           'sub.gcell@read03.f'
     &              //' ?dimension over chrg?'
     &              //' ichmx > MAX_NUM_CHRG'
     &           ,' (ichmx=',ichmx,')'
     &           ,' (MAX_NUM_CHRG@moddas.f=',MAX_NUM_CHRG,')'
            ErrID = 'L:1790/R:gcell/F:read03.f'
            call ErrWrite(ErrID,ErrCha)
         endif
         call moddas_deallocate_cha(chrg)
      return

*-----------------------------------------------------------------------
*     errors
*-----------------------------------------------------------------------

 975  continue

         write(dkam,'(i6)') kvlmax
         m_err = '# of cells exceeds kvlmax = '//dkam//
     &        '. Increase kvlmax.'
         ErrCha = ''
         ErrID = 'L:1806/R:gcell/F:read03.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------
 976  continue

         m_err = 'density should not be 0. Use inner void instead.'
         ErrCha = ''
         ErrID = 'L:1817/R:gcell/F:read03.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  978 continue

         m_err = 'memory not enough, in gcell at latmax.'
         ErrCha = ''
         ErrID = 'L:1829/R:gcell/F:read03.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  979 continue

         m_err = 'mat or rho is only available in like-but cell'
         ErrCha = ''
         ErrID = 'L:1841/R:gcell/F:read03.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  980 continue

         m_err = 'like but description is wrong'
         ErrCha = ''
         ErrID = 'L:1853/R:gcell/F:read03.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  981 continue

         m_err = 'but is missing in like-but cell'
         ErrCha = ''
         ErrID = 'L:1865/R:gcell/F:read03.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  982 continue

         write(dkam,'(i6)') jlike
         m_err = 'The like-but cell '// dkam //' does not exist'
         ErrCha = ''
         ErrID = 'L:1878/R:gcell/F:read03.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  983 continue

         m_err = 'value of cell parameter (=number) is missing'
         ErrCha = ''
         ErrID = 'L:1890/R:gcell/F:read03.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  984 continue

         m_err = 'maximum number of transform parameters is 13'
         ErrCha = ''
         ErrID = 'L:1902/R:gcell/F:read03.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  985 continue

         m_err = 'number of filling cell does not match'
         ErrCha = ''
         ErrID = 'L:1914/R:gcell/F:read03.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  986 continue

         m_err = 'LAT should be 1 or 2 or 3.'
         ErrCha = ''
         ErrID = 'L:1926/R:gcell/F:read03.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  987 continue

         m_err = 'This cell parameter is specified twice.'
         ErrCha = ''
         ErrID = 'L:1938/R:gcell/F:read03.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  988 continue

         m_err = 'cell definition is too long in characters'
         ErrCha = ''
         ErrID = 'L:1950/R:gcell/F:read03.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  990 continue

         m_err = 'CG section is specified twice.'
         ErrCha = ''
         ErrID = 'L:1962/R:gcell/F:read03.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  991 continue

         m_err = '[cell] section is specified twice.'
         ErrCha = ''
         ErrID = 'L:1974/R:gcell/F:read03.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  992 continue

         write(dkam,'(i6)') idmg(igcel)
         m_err = 'ID number of material should be'//
     &           ' -1, 0, 1 - kvmmax-1, = '// dkam
         ErrCha = ''
         ErrID = 'L:1988/R:gcell/F:read03.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  993 continue

         write(dkam,'(i6)') idrg(igcel)
         m_err = 'ID number of cell is duplicated. = '// dkam
         ErrCha = ''
         ErrID = 'L:2001/R:gcell/F:read03.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  994 continue

         write(dkam,'(i6)') idrg(igcel)
         m_err = 'ID number of cell should be'//
     &           ' 1 - kvmmax-1. = '// dkam
         ErrCha = ''
         ErrID = 'L:2015/R:gcell/F:read03.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  999 continue

         m_err = 'Description of [cell] is wrong.'
         ErrCha = ''
         ErrID = 'L:2027/R:gcell/F:read03.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

      end


************************************************************************
*                                                                      *
      subroutine gsurf(jsn,jsi,dsin,idsi,ill,ilf,
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)
*                                                                      *
*       read [surface] section of GG input files                       *
*       modified by K.Niita on 2001/02/22                              *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      include 'param.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      parameter ( ibmt = 40 )

*-----------------------------------------------------------------------

      character m_err*200
      common /error/ m_err, l_err, k_err

*-----------------------------------------------------------------------


*-----------------------------------------------------------------------

      common /inggs/  iog, igcel, ioa, igsuf, iob, igtrs
      common /inpec/  ititl, ipara, ibody, iregn, llarr, itby, itar

      common /celdb/  idsn(kvlmax), idtn(kvlmax)

*-----------------------------------------------------------------------

      character chin*200, chlw*200, chcm*200

      character dsin(0:9)*200
      dimension idsi(0:9)

      dimension ill(0:9), ilf(0:9)

      character dkam*6

      logical deqn5

*-----------------------------------------------------------------------

      dimension vdsn(50)

*-----------------------------------------------------------------------

      dimension iblg(50)
      character chsf(50)*3

      data ( chsf(i), i = 1, ibmt ) /
     &         'p  ','px ','py ','pz ','so ',
     &         's  ','sx ','sy ','sz ','c/x',
     &         'c/y','c/z','cx ','cy ','cz ',
     &         'k/x','k/y','k/z','kx ','ky ',
     &         'kz ','sq ','gq ','tx ','ty ',
     &         'tz ','x  ','y  ','z  ','box',
     &         'rpp','sph','rcc','rec','ell',
     &         'trc','wed','arb','rhp','hex'/

      data ( iblg(i), i = 1, ibmt ) /
     &           2,    2,    2,    2,    2,
     &           2,    2,    2,    2,    3,
     &           3,    3,    2,    2,    3,
     &           3,    3,    3,    2,    2,
     &           2,    2,    2,    2,    2,
     &           2,    2,    2,    2,    3,
     &           3,    3,    3,    3,    3,
     &           3,    3,    3,    3,    3/

*-----------------------------------------------------------------------

            ierr  = 0
            idrf = 0 ! zero set, S.H. 2022.3.24

            if( igsuf .ne. 0 ) goto 991

            if( iregn .ne. 0 .or. ibody .ne. 0 ) goto 990

*-----------------------------------------------------------------------
*     open temporary file : unit 23
*-----------------------------------------------------------------------

            ioa = 16
            open(ioa,status='scratch',form='unformatted')

*-----------------------------------------------------------------------
*     read one line from jsi
*-----------------------------------------------------------------------

  140 continue

            call readl(jsn,jsi,dsin,idsi,ill,ilf,'$',
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

               if( ierr .ne. 0 ) goto 1000
               if( jpn  .eq. 3 ) goto 1000

               if( iskip .ne. 0 ) goto 140

*-----------------------------------------------------------------------
*        end of region section
*-----------------------------------------------------------------------

            if( i1 .le. 5 .and. chlw(i1:i1) .eq. '[' ) then

               jpn = 1
               goto 1000

            end if

*-----------------------------------------------------------------------
*        start new surface line
*-----------------------------------------------------------------------

               ic = i1

         if( i1 .le. 5 ) then

*-----------------------------------------------------------------------
*        write information on temporary file unit = 23
*-----------------------------------------------------------------------

            if( igsuf .gt. 0 ) then

               write(ioa) idrf, idtr, idsf, igkst,
     &                    ( vdsn(i), i = 1, igkst )

            end if

*-----------------------------------------------------------------------

               igsuf = igsuf + 1

                  if( igsuf .ge. kvlmax ) goto 998

               igkst = 0

*-----------------------------------------------------------------------
*           reflection (*) or white surface (+)
*-----------------------------------------------------------------------

               idrf = 0

            if( chlw(ic:ic) .eq. '*' ) then

               idrf = -1

               ic = ic + 1

            end if

            if( chlw(ic:ic) .eq. '+' ) then

               idrf = -2

               ic = ic + 1

            end if

*-----------------------------------------------------------------------
*           surface number
*-----------------------------------------------------------------------

               call snum(chlw,ic,i3,ic2,cvvv,ierr)
               if( ierr .ne. 0 ) goto 999

                     nsufn = nint( cvvv )

                  do k = 1, igsuf - 1

                     if( nsufn .eq. idsn(k) ) goto 989

                  end do

                  idsn(igsuf) = nint( cvvv )

                  ic = jnumc(chlw,ic2,i3)

*-----------------------------------------------------------------------
*           surface transform if exists
*-----------------------------------------------------------------------

            if( deqn5( chlw(ic:ic) ) ) then

               call snum(chlw,ic,i3,ic2,cvvv,ierr)
               if( ierr .ne. 0 ) goto 999

                  idtr = nint( cvvv )

            else

                  ic2  = ic
                  idtr = 0

            end if

*-----------------------------------------------------------------------
*           kind of surface
*-----------------------------------------------------------------------

               ic = jnumc(chlw,ic2,i3)

            do i = 1, ibmt

               il = ic + iblg(i) - 1

               if( chlw(ic:il) .eq. chsf(i)(1:iblg(i)) ) goto 300

            end do

               goto 999

  300       continue

               idsf = i

               if( idsf .eq. 40 ) idsf = 39

               ic2 = il + 1

               ic = jnumc(chlw,ic2,i3)

         end if

*-----------------------------------------------------------------------
*        rest part or sequential line
*-----------------------------------------------------------------------

            if( igsuf .eq. 0 ) goto 999

  100       continue

            if( ic .gt. i3 ) goto 140

                  call snum(chlw,ic,i3,ic2,cvvv,ierr)

                  if( ierr .ne. 0 ) goto 999

                  igkst = igkst + 1

                  vdsn(igkst) = cvvv

               ic = jnumc(chlw,ic2,i3)

         goto 100

*-----------------------------------------------------------------------
*     errors
*-----------------------------------------------------------------------

  989 continue

         m_err = 'Surface number is duplicated.'
         ErrCha = ''
         ErrID = 'L:2299/R:gsurf/F:read03.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  990 continue

         m_err = 'CG section is specified twice.'
         ErrCha = ''
         ErrID = 'L:2311/R:gsurf/F:read03.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  991 continue

         m_err = '[surface] section is specified twice.'
         ErrCha = ''
         ErrID = 'L:2323/R:gsurf/F:read03.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  998 continue

         m_err = 'Number of surfaces exceeds kvlmax in param.inc.'
         ErrCha = ''
         ErrID = 'L:2335/R:gsurf/F:read03.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  999 continue

         m_err = 'Description of [surface] is wrong.'
         ErrCha = ''
         ErrID = 'L:2347/R:gsurf/F:read03.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

 1000 continue

*-----------------------------------------------------------------------
*        write information on temporary file unit = 23
*-----------------------------------------------------------------------

               write(ioa) idrf, idtr, idsf, igkst,
     &                    ( vdsn(i), i = 1, igkst )

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine gtrans(jsn,jsi,dsin,idsi,ill,ilf,
     &                  jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr)
*                                                                      *
*       read [transform] section of GG input files                     *
*       modified by K.Niita on 2010/01/14                              *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      include 'param.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      character m_err*200
      common /error/ m_err, l_err, k_err

*-----------------------------------------------------------------------

      common /inggs/ iog, igcel, ioa, igsuf, iob, igtrs

      common /celdb/ idsn(kvlmax), idtn(kvlmax)

*-----------------------------------------------------------------------

      character chin*200, chlw*200, chcm*200

      character dsin(0:9)*200
      dimension idsi(0:9)

      dimension ill(0:9), ilf(0:9)

*-----------------------------------------------------------------------

      dimension vtrs(13)

*-----------------------------------------------------------------------

            ierr  = 0

            igtms = 0

*-----------------------------------------------------------------------
*     open temporary file : unit 71
*-----------------------------------------------------------------------

         if( igtrs .eq. 0 ) then

            iob = 71

            open(iob,status='scratch',form='unformatted')

         end if

*-----------------------------------------------------------------------
*     read one line from jsi
*-----------------------------------------------------------------------

  140 continue

            call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

               if( ierr .ne. 0 ) goto 1000
               if( jpn  .eq. 3 ) goto 1000

               if( iskip .ne. 0 ) goto 140

*-----------------------------------------------------------------------
*        end of region section
*-----------------------------------------------------------------------

            if( i1 .le. 5 .and. chlw(i1:i1) .eq. '[' ) then

               jpn = 1
               goto 1000

            end if

               ic = i1

*-----------------------------------------------------------------------
*        write information on temporary file unit = 71
*-----------------------------------------------------------------------

         if( i1 .le. 5 ) then

            if( igtms .gt. 0 ) then

               write(iob) itrs, ( vtrs(i), i = 1, 13 )

            end if

*-----------------------------------------------------------------------

               vtrs( 1) = 0.0
               vtrs( 2) = 0.0
               vtrs( 3) = 0.0
               vtrs( 4) = 1.0
               vtrs( 5) = 0.0
               vtrs( 6) = 0.0
               vtrs( 7) = 0.0
               vtrs( 8) = 1.0
               vtrs( 9) = 0.0
               vtrs(10) = 0.0
               vtrs(11) = 0.0
               vtrs(12) = 1.0
               vtrs(13) = 1.0

               igtrs = igtrs + 1
               igtms = igtms + 1
               igkst = 0

            if( chlw(ic:ic+2) .ne. '*tr' .and.
     &          chlw(ic:ic+1) .ne. 'tr' ) goto 999

            if( chlw(ic:ic+2) .eq. '*tr' ) then

               ic = ic + 3
               itrs = 1

               vtrs( 4) =  0.0
               vtrs( 5) = 90.0
               vtrs( 6) = 90.0
               vtrs( 7) = 90.0
               vtrs( 8) =  0.0
               vtrs( 9) = 90.0
               vtrs(10) = 90.0
               vtrs(11) = 90.0
               vtrs(12) =  0.0

            else if( chlw(ic:ic+1) .eq. 'tr' ) then

               ic = ic + 2
               itrs = 0

            end if

               icf = inumc(chlw,ic,i3,' ') - 1

               call snum(chlw,ic,icf,ic2,cvvv,ierr)

               if( ierr .ne. 0 ) goto 999
               if( nint(cvvv) .gt. 999999 ) goto 998

               idtn(igtrs) = nint( cvvv )

               ic = jnumc(chlw,ic2,i3)

         end if

*-----------------------------------------------------------------------
*        rest part or sequential line
*-----------------------------------------------------------------------

            if( igtms .eq. 0 ) goto 999

  100       continue

            if( ic .gt. i3 ) goto 140

                  call snum(chlw,ic,i3,ic2,cvvv,ierr)

                  if( ierr .ne. 0 ) goto 999

                  igkst = igkst + 1

                  if( igkst .gt. 13 ) goto 999

                  vtrs(igkst) = cvvv

               ic = jnumc(chlw,ic2,i3)

*-----------------------------------------------------------------------

      goto 100

 1000 continue

*-----------------------------------------------------------------------
*        write information on temporary file unit = 71
*-----------------------------------------------------------------------

            write(iob) itrs, ( vtrs(i), i = 1, 13 )

            return

*-----------------------------------------------------------------------
*     errors
*-----------------------------------------------------------------------

  999 continue

         m_err = 'Description of [transform] is wrong.'
         ErrCha = ''
         ErrID = 'L:2568/R:gtrans/F:read03.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

  998 continue

         m_err = 'ID of transform should be 1 - 999999.'
         ErrCha = ''
         ErrID = 'L:2578/R:gtrans/F:read03.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

      end

************************************************************************
*                                                                      *
      subroutine latencode(latot,lattmp,latot2,lattmp2)
*                                                                      *
*       encode voxel elements; lattmp, to compressed format; lattmp2   *
*       created by T.Furuta on 2018/10/09                              *
*                                                                      *
************************************************************************
      implicit none
      integer,intent(in) :: latot,lattmp(latot)
      integer,intent(out) :: latot2,lattmp2(latot)
      integer i,icount,jcount
      integer latold
*-----------------------------------------------------------------------
      latold=lattmp(1)
      lattmp2(1)=lattmp(1)
      icount=0
      jcount=1
      do i=2,latot
       if(lattmp(i).eq.latold)then
        icount=icount+1
       else
        if(icount.gt.0)then
         jcount=jcount+1
         lattmp2(jcount)=-icount
        endif
        jcount=jcount+1
        lattmp2(jcount)=lattmp(i)
        latold=lattmp(i)
        icount=0
       endif
      enddo
      if(icount.gt.0)then
       jcount=jcount+1
       lattmp2(jcount)=-icount
      endif
      latot2=jcount
*-----------------------------------------------------------------------
      return
      end

************************************************************************
*                                                                      *
      subroutine latdecode(latot2,lattmp2,latot,lattmp,icount,ierr)
*                                                                      *
*       decode compressed format; lattmp2, into voxel elements; lattmp *
*       created by T.Furuta on 2018/10/09                              *
*                                                                      *
************************************************************************
      implicit none
      integer,intent(in) :: latot2,lattmp2(latot2),latot
      integer,intent(out) :: lattmp(latot),icount,ierr
      integer i,j
      integer latold
*-----------------------------------------------------------------------
      icount=0
      do i=1,latot2
       if(lattmp2(i).gt.0)then
        icount=icount+1
        latold=lattmp2(i)
        lattmp(icount)=lattmp2(i)
       else
        do j=1,-lattmp2(i)
         lattmp(icount+j)=latold
        enddo
        icount=icount-lattmp2(i)
       endif
      enddo
      if(icount.eq.latot)then
       ierr=0
      else
       ierr=1
      endif
*-----------------------------------------------------------------------
      return
      end

