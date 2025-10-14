!***********************************************************************
!                                                                      *
      subroutine read_sangel(jsn,jsi,dsin,idsi,ill,ilf,
     &                       jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr
     &                      ,ic ,icl ,tname)
!                                                                      *
!     read sangel parameter from input file                            *
!     and write for scrach file                                        *
!     create by T.Miura on 2017/09/30                                  *
!                                                                      *
!***********************************************************************
      use sangelmod

      implicit none

      include 'param.inc'
!-----------------------------------------------------------------------
      integer   jsn                ! include level
      include 'err.inc'
      integer   jsi                ! I/O machine number
      character dsin(0:9)*200      ! input file name
      integer   idsi(0:9)          ! length of input file name
      integer   ill (0:9)          ! input file lines
      integer   ilf (0:9)          ! final line to read
      integer   jpn
      character chin*200           ! read data (raw letters)
      character chlw*200           ! read data (lower letters)
      character chcm*200           ! read data (removing the blank)
      integer   i1       ! starting column (non-blank)
      integer   i2       ! ending column   (non-blank)
      integer   i3       ! ending column (w/o comment & last part blank)
      integer   i4       ! ending column (compress letters)
      integer   iskip    ! 0:normal line|1:blank line|2:comment line
      integer   ierr               ! return code (=0:normal retrun)

      integer   ic
      integer   icl
      character*(*),intent(in) :: tname    ! tally name
!-----------------------------------------------------------------------
      character m_err*200
      integer   l_err, k_err
      common /error/ m_err, l_err, k_err
!-----------------------------------------------------------------------
      integer   itnm
      integer   ital
      integer   itals
      integer   italm
      common /tall00/ itnm, ital(itlmax), itals(itlmax), italm(itlmax)
!-----------------------------------------------------------------------
      integer   :: isc = 202               ! input unit number for scrach file
      integer   :: ios

      integer   :: i
      integer   :: isang
      double precision :: cvvv
!-----------------------------------------------------------------------
       integer    inumc, jnumc
       external   inumc, jnumc
!-----------------------------------------------------------------------
      ierr = 0
      jpn  = 0
!-----------------------------------------------------------------------
!     Open scratch file for sangel parameter
!-----------------------------------------------------------------------
      if (itsanm == 0) then
         open(isc,form='unformatted',status='scratch')
      end if

  100 continue
!-----------------------------------------------------------------------
!     Get number of line for sangel parameters
!-----------------------------------------------------------------------
      call onum(chlw,ic,icl,cvvv,ierr)
      if( ierr .ne. 0 ) goto 997

      isang = nint( cvvv )
      if( isang <= 0 ) goto 951
      itsanm = itsanm + isang                 ! total sangels line
!-----------------------------------------------------------------------
!     Loop for sangel parameter lines
!-----------------------------------------------------------------------
      do i=1, isang
!.......................................................................
!     read one line from jsi
!.......................................................................
  140    ill(jsn) = ill(jsn)+1
         read(jsi,'(a200)',iostat=ios) chin   ! read from input file
         if ( ios == -1 ) goto 251            ! EOF
  240    continue
         if ( ill(jsn) <= ilf(jsn) ) goto 250

!-----------------------------------------------------------------------
!        close files
!-----------------------------------------------------------------------
  251    continue
         if ( jsn <= 0 ) then
            jpn = 3
            goto 954
         else if ( jsn > 0 ) then
            call closef(jsi,jsn)
            goto 140
         end if

!-----------------------------------------------------------------------
!        end of the section
!-----------------------------------------------------------------------
  250    continue
         chlw = adjustl(chin(1:5))
         if ( chlw(1:1) .eq. '[' ) then
            jpn = 1
            goto 954
         end if

!-----------------------------------------------------------------------

         write(isc) itnm+1, i, chin     ! write for scrach file

      end do

      return
!-----------------------------------------------------------------------
!     errors
!-----------------------------------------------------------------------
  951 m_err = 'Number of sangel is little equal 0 in tally '//tname
      ErrCha = ''
      ErrID = 'L:126/R:read_sangel/F:sangel.f'
      goto 999


  953 m_err = 'Input data ended when reading the sangel parameter '
     &      //'in tally '//tname
      ErrCha = ''
      ErrID = 'L:133/R:read_sangel/F:sangel.f'
      goto 999

  954 m_err = 'The number of parameters of sangel is wrong in tally '
     &       //tname
      ErrCha = ''
      ErrID = 'L:139/R:read_sangel/F:sangel.f'
      goto 999

  997 m_err = 'Description of parameter is wrong in tally'//tname
      ErrCha = ''
      ErrID = 'L:144/R:read_sangel/F:sangel.f'
      goto 999
!-----------------------------------------------------------------------
  999 continue
      l_err = ill(jsn)
      k_err = jsn
      ierr  = 1
!-----------------------------------------------------------------------
 9000 continue
      close(isc)
      return

      end subroutine read_sangel



!***********************************************************************
!                                                                      *
      subroutine set_sangel
!                                                                      *
!     read sangel parameter from scrach file                           *
!     create by T.Miura on 2017/09/30                                  *
!                                                                      *
!***********************************************************************
      use sangelmod

      implicit none
!-----------------------------------------------------------------------
      include 'param.inc'
!-----------------------------------------------------------------------
      integer   :: ital_now
      integer   :: iseq
      character(len=200) :: chin           ! read data (raw letters)

      integer   :: ital_old

      integer   :: isc = 202               ! input unit number for scrach file
      integer   :: i
!-----------------------------------------------------------------------
      if ( itsanm <= 0 ) return
      endfile(isc)
      rewind (isc)
!-----------------------------------------------------------------------
      ital_old  = 0
      do i=1, itsanm
         read(isc) ital_now, iseq, chin     ! read from scrach file

         itsang(i) = trim(adjustl(chin))    ! sangel parameter
         ltsang(i) = len_trim(itsang(i))    ! sangel parameter length

         if ( ital_old /= ital_now ) then
            ltsans(ital_now) = i
            ital_old = ital_now
         end if
         itsans(ital_now) = i - ltsans(ital_now) + 1
      end do
!-----------------------------------------------------------------------
 9000 continue
      close(isc)
      return

      end subroutine set_sangel



!***********************************************************************
!                                                                      *
      subroutine write_sangel(iot,m,iecho)
!                                                                      *
!     write sangel parameter                                           *
!     create by T.Miura on 2017/09/30                                  *
!                                                                      *
!***********************************************************************
      use sangelmod

      implicit none
!-----------------------------------------------------------------------
      integer ,intent(in) :: iot      ! write i/o
      integer ,intent(in) :: m        ! tally
      integer ,intent(in) :: iecho    ! echo flag (=1:echo,=0:not echo)

      integer   :: i
!-----------------------------------------------------------------------
      if ( itsans(m) <= 0 ) return

      if ( iecho == 1 .and. iot /= 28) then   ! Echo output
         do i = ltsans(m), ltsans(m)+itsans(m)-1
            write(iot,'("# ",a)') itsang(i)(1:ltsang(i))
         end do

      else   ! tally data or phits.out output
         do i = ltsans(m), ltsans(m)+itsans(m)-1
            write(iot,'(a)') itsang(i)(1:ltsang(i))
         end do
      end if
!-----------------------------------------------------------------------
      return

      end subroutine write_sangel
