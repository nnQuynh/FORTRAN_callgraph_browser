!***********************************************************************
!                                                                      *
      subroutine init_sumtal
!                                                                      *
!     create by T.Miura on 2014/11/30                                  *
!     initialize for sumtally variables                                *
!     Modified for responding to multiple tallies. (T.Miura 2015/07/31)*
!                                                                      *
!***********************************************************************

       use sumtallymod

       implicit none

!-----------------------------------------------------------------------

       include 'param.inc'

!-----------------------------------------------------------------------

       integer   i

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!    allocated variable arrays
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

       if( itlmax > 0 ) then
        allocate (  nsumtalRead(itlmax) )

        allocate( isumtally  (itlmax))
        allocate( nfile      (itlmax))
        allocate( sumWR      (itlmax))
        allocate(  sfile     (itlmax))
        allocate( lsfile     (itlmax))
        allocate( sumfactor  (itlmax))
        allocate(  sumang    (itlmax))
        allocate( lsumang    (itlmax))

        allocate( resc2SUMTAL(itlmax))
        allocate( resc3SUMTAL(itlmax))
       end if

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!    initilized variable arrays
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

 1000  continue
       do i = 1, itlmax
        nsumtalRead(i) = 0
        isumtally  (i) = 0
         nfile     (i) = 0
         sumWR     (i) = 0.0d0

         sfile     (i) = ' '
        lsfile     (i) = 0

         sumfactor (i) = 0.0d0

         sumang    (i) = ' '
        lsumang    (i) = 0


        resc2SUMTAL(i) = 0.0d0
        resc3SUMTAL(i) = 0.0d0
       end do

!-----------------------------------------------------------------------

 9999  continue
       return

      end subroutine init_sumtal
!***********************************************************************


!***********************************************************************
!                                                                      *
      subroutine read_sumtal(jsn,jsi,dsin,idsi,ill,ilf,
     &                       jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)
!                                                                      *
!     read sumtally sub section from input file                        *
!     and write for scrach file                                        *
!     create by T.Miura on 2014/11/30                                  *
!                                                                      *
!***********************************************************************

       use sumtallymod

       implicit none

!-----------------------------------------------------------------------

       include 'param.inc'
       include 'err.inc'

       integer   mh
       integer   mxcval
       integer   icolm
       integer   ichrl
       integer   inig
       integer   mc
       integer   numtic
       integer   ipsm
       integer   nr
       integer   nc
       double precision r0max
       double precision r1max
       double precision r2max
       double precision r3max ! S.H. xorshift (2020.2.6)
       double precision r9max
       double precision r9min
       double precision pi
       integer   jol
       integer   jil
       integer   jhc
       integer   jhl
       integer   jer
       integer   jps
       integer   jab
       integer   jht
       integer   jhs
       integer   jwt
       integer   jhm
       integer   jha
       integer   jhy
       integer   jhb
       integer   jhd
       integer   jhp
       integer   jhq
       integer   jhr

       include 'angel01.inc'

!-----------------------------------------------------------------------

       common /inout/  in,io
       integer   in,io

       integer   jsn
       integer   jsi
       character dsin(0:9)*200
       integer   idsi(0:9)
       integer   ill (0:9)
       integer   ilf (0:9)
       integer   jpn
       character chin*200
       character chlw*200
       character chcm*200
       integer   i1, i2, i3, i4
       integer   iskip
       integer   ierr           ! return code (=0:normal retrun)

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

       double precision cval
       double precision aval

       common /rval1/ cval(mxcval), aval(mxcval)

!-----------------------------------------------------------------------

       character tname*12
       data      tname /'sumtally end'/

       integer   isc            ! output unit number for scrach file

       integer   nsumt
       integer   i

!-----------------------------------------------------------------------

       ierr = 0
       jpn  = 0

       isc = 151
       nsumt = 0
       do i = 1, itnm+1
        nsumt = nsumt + nsumtalRead(i)
       end do
       if( nsumt == 0 ) then
        open(isc,form='unformatted',status='scratch')
       end if
       write(isc) chlw       ! sumtally start
       write(isc) jsn        ! include level
       write(isc) dsin(jsn)  ! input file name
       write(isc) idsi(jsn)  ! length of input file name
       write(isc) ill(jsn)   ! input file lines
       write(isc) ilf(jsn)   ! final line to read
       write(isc) itnm       ! tally naumber - 1
       write(isc) cval       ! saved user defined variables. c1-c100
       nsumtalRead(itnm+1) = nsumtalRead(itnm+1) + 1 ! number of sumtally

!-----------------------------------------------------------------------
!     read one line from jsi
!-----------------------------------------------------------------------

  100  continue

       call sumtal_readlt(isc,jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &            jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

       if( ierr .ne. 0 ) go to 9000
       if( iskip /= 0 ) then
        write(isc) chin
        goto 100
       end if
       if( jpn == 3 ) then        ! q:
        write(isc) 'q:'
        goto 9000
       else if( jpn == 1 ) then        ! qp:
        write(isc) 'qp:'
        goto 9000
       end if

!-----------------------------------------------------------------------
!     end of the section (error)
!-----------------------------------------------------------------------

       if( chlw(i1:i1) .eq. '[' ) then
        m_err = 'Not found sumtally end'
        ErrCha = ''
        ErrID = 'L:237/R:read_sumtal/F:sumtally.f'
        l_err = ill(jsn)
        k_err = jsn
        ierr = 1
        goto 9000
       end if

       write(isc) chin
       if ( chlw(i1:i3) .ne. tname ) go to 100    ! end of sumtally sub section?

!-----------------------------------------------------------------------

 9000  continue
       return

      end subroutine read_sumtal
!***********************************************************************


!***********************************************************************
!                                                                      *
      subroutine read_sumtal_dchain(jsn,jsi,dsin,idsi,ill,ilf,
     &                       jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)
!                                                                      *
!     read sumtally sub section of the t-dchain tally from input file  *
!     and write to scratch file                                        *
!     create by T.Miura on 2016/08/31                                  *
!                                                                      *
!***********************************************************************

       use sumtallymod

       implicit none

!-----------------------------------------------------------------------

       include 'param.inc'
       include 'err.inc'

       integer   mh
       integer   mxcval
       integer   icolm
       integer   ichrl
       integer   inig
       integer   mc
       integer   numtic
       integer   ipsm
       integer   nr
       integer   nc
       double precision r0max
       double precision r1max
       double precision r2max
       double precision r3max ! S.H. xorshift (2020.2.6)
       double precision r9max
       double precision r9min
       double precision pi
       integer   jol
       integer   jil
       integer   jhc
       integer   jhl
       integer   jer
       integer   jps
       integer   jab
       integer   jht
       integer   jhs
       integer   jwt
       integer   jhm
       integer   jha
       integer   jhy
       integer   jhb
       integer   jhd
       integer   jhp
       integer   jhq
       integer   jhr

       include 'angel01.inc'

!-----------------------------------------------------------------------

       common /inout/  in,io
       integer   in,io

       integer   jsn
       integer   jsi
       character dsin(0:9)*200
       integer   idsi(0:9)
       integer   ill (0:9)
       integer   ilf (0:9)
       integer   jpn
       character chin*200
       character chlw*200
       character chcm*200
       integer   i1, i2, i3, i4
       integer   iskip
       integer   ierr           ! return code (=0:normal retrun)

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

       double precision cval
       double precision aval

       common /rval1/ cval(mxcval), aval(mxcval)

!-----------------------------------------------------------------------

       character tname*12
       data      tname /'sumtally end'/

       integer   isc             ! output unit number for scrach file
       integer   isct            ! in/out unit number for temp. scrach file

       integer   nsumt
       integer   i
       integer   itnmw

!-----------------------------------------------------------------------

       ierr = 0
       jpn  = 0

       isc = 151
       isct= 151+1
       nsumt = 0
       do i = 1, itnm+1
        nsumt = nsumt + nsumtalRead(i)
       end do
       if( nsumt == 0 ) then
        open(isc,form='unformatted',status='scratch')
       end if
       open(isct,form='unformatted',status='scratch',action='readwrite')
       write(isct) chlw       ! sumtally start
       write(isct) jsn        ! include level
       write(isct) dsin(jsn)  ! input file name
       write(isct) idsi(jsn)  ! length of input file name
       write(isct) ill(jsn)   ! input file lines
       write(isct) ilf(jsn)   ! final line to read
       write(isct) itnm       ! tally naumber - 1
       write(isct) cval       ! saved user defined variables. c1-c100
       nsumtalRead(itnm+1) = nsumtalRead(itnm+1) + 1 ! t-yield

!-----------------------------------------------------------------------
!     read one line from jsi
!-----------------------------------------------------------------------

  100  continue

       call sumtal_readlt(isct,jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &            jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

       if( ierr .ne. 0 ) go to 9000
       if( iskip /= 0 ) then
        write(isct) chin
        goto 100
       end if
       if( jpn == 3 ) then        ! q:
        write(isct) 'q:'
        goto 9000
       else if( jpn == 1 ) then        ! qp:
        write(isct) 'qp:'
        goto 9000
       end if

!-----------------------------------------------------------------------
!     end of the section (error)
!-----------------------------------------------------------------------

       if( chlw(i1:i1) .eq. '[' ) then
        m_err = 'Not found sumtally end'
        ErrCha = ''
        ErrID = 'L:421/R:read_sumtal_dchain/F:sumtally.f'
        l_err = ill(jsn)
        k_err = jsn
        ierr = 1
        goto 9000
       end if


       write(isct) chin
       if ( chlw(i1:i3) .ne. tname ) go to 100    ! end of sumtally sub section?

!-----------------------------------------------------------------------
!     temp. scrach file ===> scrach file
!-----------------------------------------------------------------------

       endfile(isct)

       do i=1, 3               ! t-yield / t-track / t-dchain
        rewind (isct)

        do
         read(isct,end=999) chlw
         call chlngt(chlw,200,i1,i2)
         call chcaps(chlw,i1,i2,i3,'#!$')
         if ( chlw(i1:i3) == 'sumtally start' ) then
          read(isct) jsn        ! include level
          read(isct) dsin(jsn)  ! input file name
          read(isct) idsi(jsn)  ! length of input file name
          read(isct) ill(jsn)   ! input file lines
          read(isct) ilf(jsn)   ! final line to read
          read(isct) itnmw      ! tally number - 1
          read(isct) cval       ! set user defined variables. c1-c100

          write(isc) chlw
          write(isc) jsn        ! include level
          write(isc) dsin(jsn)  ! input file name
          write(isc) idsi(jsn)  ! length of input file name
          write(isc) ill(jsn)   ! input file lines
          write(isc) ilf(jsn)   ! final line to read
          write(isc) itnmw+i-1  ! tally number - 1
          write(isc) cval       ! set user defined variables. c1-c100

         else if ( chlw == 'incinc' .or.
     &             chlw == 'incdec' ) then
          read(isct) jsn        ! include level
          read(isct) dsin(jsn)  ! input file name
          read(isct) idsi(jsn)  ! length of input file name
          read(isct) ill(jsn)   ! input file lines
          read(isct) ilf(jsn)   ! final line to read

          write(isc) chlw
          write(isc) jsn        ! include level
          write(isc) dsin(jsn)  ! input file name
          write(isc) idsi(jsn)  ! length of input file name
          write(isc) ill(jsn)   ! input file lines
          write(isc) ilf(jsn)   ! final line to read

         else
          write(isc) chlw
         end if
        end do
  999   continue

       end do

       nsumtalRead(itnm+2) = nsumtalRead(itnm+2) + 1 ! t-track
       nsumtalRead(itnm+3) = nsumtalRead(itnm+3) + 1 ! t-dchain

!-----------------------------------------------------------------------

 9000  continue
       close(isct)
       return

      end subroutine read_sumtal_dchain
!***********************************************************************


!***********************************************************************
!                                                                      *
      subroutine set_sumtal(io,jo,ierr)
!                                                                      *
!     read sumtally sub section data from scrach file                  *
!     create by T.Miura on 2014/11/30                                  *
!     modified by T.Miura on 2015/07/31 (multiple tally inputs)        *
!                                                                      *
!***********************************************************************

       use sumtallymod

       implicit none

!-----------------------------------------------------------------------

       include 'param.inc'
       include 'err.inc'

       integer   mh
       integer   mxcval
       integer   icolm
       integer   ichrl
       integer   inig
       integer   mc
       integer   numtic
       integer   ipsm
       integer   nr
       integer   nc
       double precision r0max
       double precision r1max
       double precision r2max
       double precision r3max ! S.H. xorshift (2020.2.6)
       double precision r9max
       double precision r9min
       double precision pi
       integer   jol
       integer   jil
       integer   jhc
       integer   jhl
       integer   jer
       integer   jps
       integer   jab
       integer   jht
       integer   jhs
       integer   jwt
       integer   jhm
       integer   jha
       integer   jhy
       integer   jhb
       integer   jhd
       integer   jhp
       integer   jhq
       integer   jhr

       include 'angel01.inc'

!-----------------------------------------------------------------------

       integer   io
       integer   jo

!-----------------------------------------------------------------------

       integer   jsn
       character dsin(0:9)*200
       integer   idsi(0:9)
       integer   ill (0:9)
       integer   ilf (0:9)
       integer   jpn
       character chin*200
       character chlw*200
       character chcm*200
       integer   i1, i2, i3, i4, i5
       integer   iskip
       integer   ierr           ! return code (=0:normal retrun)

       integer   itnmw          ! tally number - 1
       integer   iditnmw        ! for do-loop

!-----------------------------------------------------------------------

       integer   itnm
       integer   ital
       integer   itals
       integer   italm

       common /tall00/ itnm, ital(itlmax), itals(itlmax), italm(itlmax)

!-----------------------------------------------------------------------

       integer   itfln
       integer   itfll
       character ctfln*100

       common /tall11/ itfln(itlmax), itfll(itlmax,6), ctfln(itlmax,6)

!-----------------------------------------------------------------------

       integer   icntl
       integer   inucr

       common /tcntl/  icntl, inucr

!-----------------------------------------------------------------------

       character m_err*200
       integer   l_err, k_err
       common /error/ m_err, l_err, k_err

!-----------------------------------------------------------------------

       double precision cval
       double precision aval

       common /rval1/ cval(mxcval), aval(mxcval)

!-----------------------------------------------------------------------

       integer, parameter :: icsu = 5
       integer lschn(0:icsu+1), ischn(0:icsu+1)
       character schan(0:icsu+1)*13

       integer i
       integer irtal
       integer j

       data ( schan(i), i = 0, icsu+1 )              /'sumtallystart',
     &'isumtally    ','nfile        ','sfile        ','sumfactor    ',
     &'angel        ',
     &'sumtallyend  '/

       data ( lschn(i), i = 0, icsu+1 ) /                  13,
     &            9,            5,            5,            9,
     &            5,
     &           11 /

       character tname*8
       data      tname /'Sumtally'/

       integer   isc            ! input unit number for scrach file
       integer   ipm
       integer   ic
       integer   ic2
       integer   icl
       integer   il
       integer   icf
       integer   ict
       character chlc*200
       double precision cvvv

       integer infile          ! file counter for work
       integer nfilew          ! for work (nfile)

!-----------------------------------------------------------------------

       integer    inumc, jnumc
       external   inumc, jnumc

!-----------------------------------------------------------------------

       do irtal=1, itnm
        if( nsumtalRead(irtal) > 0 ) goto 8   ! check sumtally sub section
       end do
       return   ! no sumtally sub section

    8  continue
       ierr = 0

       jsn  = 0
       jpn  = 0
       ilf(jsn) = 10000000

       m_err = ' Error !!'
       ErrCha = ''
       ErrID = 'L:674/R:set_sumtal/F:sumtally.f'
       l_err = 1
       k_err = 0

       dsin(k_err) = 'Error Line'
       idsi(k_err) = 12
!-----------------------------------------------------------------------
       isc = 151
       endfile(isc)                 !Writes the file end, corresponding to the 'q:' & 'qp:' command
       rewind(isc)

       ischn(0)      = 0            ! count of sumtally start
       ischn(icsu+1) = 0            ! count of sumtally end

       nfilew        = 0            ! max. number of tally files

    1  continue                     ! next sumtally sub-section
       read(isc,end=99) chin        ! sumtally start
       call chlngt(chin,200,i1,i2)
       chlw = chin
       call chcaps(chlw,i1,i2,i3,'#!$')
       if( chlw(i1:i3) .eq. 'sumtally start' ) then
        read(isc) jsn        ! include level
        read(isc) dsin(jsn)  ! input file name
        read(isc) idsi(jsn)  ! length of input file name
        read(isc) ill(jsn)   ! input file lines
        read(isc) ilf(jsn)   ! final line to read
        read(isc) itnmw      ! tally number - 1
        read(isc) cval       ! set user defined variables. c1-c100
        ischn(0) = ischn(0) + 1
       end if

       do i = 1, icsu
        ischn(i) = 0
       end do

!=======================================================================
!<<<<< pre read (get allocate array size) >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!=======================================================================

   10  continue

       call sumtal_readl(jsn,isc,dsin,idsi,ill,ilf,'#!$',
     &                   jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

       if( ierr .ne. 0 ) goto 9000
       if( chlw(i1:i1+1) .eq. 'q:' ) goto 99
       if( chlw(i1:i1+2) .eq. 'qp:' ) goto 99
       if( jpn  .eq. 3 ) goto 99      ! scratch file end

       if( iskip .ne. 0 ) goto 10

       if( chlw(i1:i3) .eq. 'sumtally end' ) then
        ischn(icsu+1) = ischn(icsu+1) + 1
        if( ischn(0) == ischn(icsu+1) ) then
         goto 1     ! next sumtally sub sction data
        else
         goto 901   ! sumtally start/end is not paired
        end if

       else if( chlw .eq. 'incinc' ) then
        read(isc) jsn        ! include level
        read(isc) dsin(jsn)  ! input file name
        read(isc) idsi(jsn)  ! length of input file name
        read(isc) ill(jsn)   ! input file lines
        read(isc) ilf(jsn)   ! final line to read
        goto 10              ! next Line

       else if( chlw .eq. 'incdec' ) then
        read(isc) jsn        ! include level
        read(isc) dsin(jsn)  ! input file name
        read(isc) idsi(jsn)  ! length of input file name
        read(isc) ill(jsn)   ! input file lines
        read(isc) ilf(jsn)   ! final line to read
        goto 10              ! next Line

       else if( chlw(i1:i1+3) .eq. 'set:' ) then
        write(io,'(/" *** Warning : set: definition in"
     &             " sumtally subsection is ignored."/)')
        write(io,*) dsin(jsn)(1:idsi(jsn)),ill(jsn),':'

        ErrCha = ''
        MsgID = 'L:756/R:set_sumtal/F:sumtally.f'
        call ErrWrite(MsgID, ErrCha)
        write(jo,'(/" *** Warning : set: definition in"
     &             " sumtally subsection is ignored."/)')
        write(jo,*) dsin(jsn)(1:idsi(jsn)),ill(jsn),':'
        goto 10              ! next Line

       end if

   15  continue
       if( ierr .ne. 0 ) goto 9000
       if( jpn  .eq. 3 ) goto 99      ! scratch file end

!-----------------------------------------------------------------------
!        identify the parameters
!-----------------------------------------------------------------------

       icl = i1
   20  continue
       chlc = chlw
       call chcomp(chlc,icl,i3,i5)
       do i = 1, icsu
        il = icl + lschn(i) - 1
        if( chlc(icl:il) .eq. schan(i)(1:lschn(i)) ) goto 30
       end do

       goto 987

!-----------------------------------------------------------------------
!        read value of parameters
!-----------------------------------------------------------------------

   30  continue

       ipm = i
       ischn( ipm ) = ischn( ipm ) + 1

       ic = inumc(chlw,il+1,i3,'=') + 1
       ic = jnumc(chlw,ic,i3)

       if( ic .gt. i3 ) goto 997

       icl = inumc(chlw,ic,i3,';') - 1

!-----------------------------------------------------------------------
!      isumtally
!-----------------------------------------------------------------------

       if( ipm == 1 ) then
        call onum(chlw,ic,icl,cvvv,ierr)
        if( ierr .ne. 0 ) goto 997

        icl = jnumc(chlw,icl+2,i3)
        if( icl .le. i3 ) goto 15

!-----------------------------------------------------------------------
!      nfile, file name and weight rate
!-----------------------------------------------------------------------

       else if( ipm == 2 ) then
        call onum(chlw,ic,icl,cvvv,ierr)
        if( ierr .ne. 0 ) goto 997

        if( chlw(icl+1:icl+1) .eq. ';' ) goto 993

        nfile(itnmw+1) = nint( cvvv )
        if( nfile(itnmw+1) <= 0 ) then             ! nfile <= 0
         goto 951
        else if( nfile(itnmw+1) > nfilew ) then
         nfilew = nfile(itnmw+1)                   ! Max nfile
        end if

!-----------------------------------------------------------------------
   50   continue
        infile = 0

        ic = i3+1                       ! set the value to read next line

        do i=1, nfile(itnmw+1)

         if( ic .gt. i3 ) then
   52     call sumtal_readl(jsn,isc,dsin,idsi,ill,ilf,'#!$',
     &                      jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)
          if( ierr .ne. 0 ) goto 9000
          if( jpn  .eq. 3 ) goto 979
          if( iskip .ne. 0 ) goto 52
          if( chlw .eq. 'incinc' ) then
           read(isc) jsn      ! include level
           read(isc) dsin(jsn)  ! input file name
           read(isc) idsi(jsn)  ! length of input file name
           read(isc) ill(jsn)   ! input file lines
           read(isc) ilf(jsn)   ! final line to read
           goto 52              ! next Line

          else if( chlw .eq. 'incdec' ) then
           read(isc) jsn      ! include level
           read(isc) dsin(jsn)  ! input file name
           read(isc) idsi(jsn)  ! length of input file name
           read(isc) ill(jsn)   ! input file lines
           read(isc) ilf(jsn)   ! final line to read
           goto 52              ! next Line

          else if( chlw(i1:i1+3) .eq. 'set:' ) then
           write(io,'(/" *** Warning : set: definition in"
     &                 " sumtally subsection is ignored."/)')
           write(io,*) dsin(jsn)(1:idsi(jsn)),ill(jsn),':'

           ErrCha = ''
           MsgID = 'L:864/R:set_sumtal/F:sumtally.f'
           call ErrWrite(MsgID, ErrCha)
           write(jo,'(/" *** Warning : set: definition in"
     &                " sumtally subsection is ignored."/)')
           write(jo,*) dsin(jsn)(1:idsi(jsn)),ill(jsn),':'
           goto 52              ! next Line

          end if
          ic = i1
         end if


         ic = jnumc(chlw,ic,i3)                           ! file name start col
         icf = min( inumc(chlw,ic,i3,' ') - 1, i3 )
!-----------------------------------------------------------------------

         ic = jnumc(chlw,icf+2,i3)                       ! weight val start col
         call snum(chlw,ic,i3,ic2,cvvv,ierr)
         if( ierr .ne. 0 ) goto 979

         ic = ic2                                        ! next col
        end do

!-----------------------------------------------------------------------
!      sfile
!-----------------------------------------------------------------------

       else if( ipm == 3 ) then
        icf = min( inumc(chlw,ic,icl,' ') - 1, icl )

        icl = jnumc(chlw,icl+2,i3)
        if( icl .le. i3 ) goto 15

!-----------------------------------------------------------------------
!      sumfactor
!-----------------------------------------------------------------------

       else if( ipm == 4 ) then
        call onum(chlw,ic,icl,cvvv,ierr)
        if( ierr .ne. 0 ) goto 997

        icl = jnumc(chlw,icl+2,i3)
        if( icl .le. i3 ) goto 15

!-----------------------------------------------------------------------
!      angel parameter
!-----------------------------------------------------------------------

       else if( ipm == 5 ) then
        ict = min( ic + 199, i2 )

       else
       end if

       go to 10   ! read next line

!=======================================================================

!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!    allocate and initilized variable arrays
!>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

   99  continue
       if ( icntl /= 13 ) return   ! Not sumally calculation
       if( ischn(0) /= ischn(icsu+1) ) then
        goto 901   ! sumtally start/end is not paired
       end if

       rewind(isc)

       if( nfilew > 0 ) then
        allocate (  tallyfname(nfilew,itlmax) )
        allocate ( ltallyfname(nfilew,itlmax) )
        allocate (  weightRate(nfilew,itlmax) )
       else
        goto 951
       end if

         tallyfname(1:nfilew,1:itlmax) = ' '
        ltallyfname(1:nfilew,1:itlmax) = 0
         weightRate(1:nfilew,1:itlmax) = 0.0d0

       ischn(0)      = 0
       ischn(icsu+1) = 0

   98  continue                      ! next sumtally sub-section
       read(isc,end=800) chin        ! sumtally start
       call chlngt(chin,200,i1,i2)
       chlw = chin
       call chcaps(chlw,i1,i2,i3,'#!$')
       if( chlw(i1:i3) .eq. 'sumtally start' ) then
        read(isc) jsn        ! include level
        read(isc) dsin(jsn)  ! input file name
        read(isc) idsi(jsn)  ! length of input file name
        read(isc) ill(jsn)   ! input file lines
        read(isc) ilf(jsn)   ! final line to read
        read(isc) itnmw      ! tally number - 1
        read(isc) cval       ! set user defined variables. c1-c100
        ischn(0) = ischn(0) + 1

        if( nsumtalRead(itnmw+1) /= 0 ) then   ! replacing flag to tally number
         nsumtalRead(itnmw+1) = ital(itnmw+1)  ! number of sumtally
        else
         nsumtalRead(itnmw+1) = 0
        end if
       end if

       do i = 1, icsu
        ischn(i) = 0
       end do

!=======================================================================
!<<<<< Production read >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
!=======================================================================
!      set default or initial value
!-----------------------------------------------------------------------

       isumtally  (itnmw+1) = 1       ! default value
        nfile     (itnmw+1) = 0
        sumWR     (itnmw+1) = 0.0d0
        sfile     (itnmw+1) = ' '
       lsfile     (itnmw+1) = 0
       sumfactor  (itnmw+1) = 1.0d0   ! default value
        sumang    (itnmw+1) = ' '
       lsumang    (itnmw+1) = 0

       sumWR      (itnmw+1) = 0.0d0   ! Sum of the weight rates

!-----------------------------------------------------------------------
  100  continue

       call sumtal_readl(jsn,isc,dsin,idsi,ill,ilf,'#!$',
     &                   jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

       if( ierr .ne. 0 ) goto 9000
       if( chlw(i1:i1+1) .eq. 'q:' ) goto 800
       if( chlw(i1:i1+2) .eq. 'qp:' ) goto 800
       if( jpn  .eq. 3 ) goto 800              ! go to parameters check

       if( iskip .ne. 0 ) goto 100

       if( chlw(i1:i3) .eq. 'sumtally end' ) then
        ischn(icsu+1) = ischn(icsu+1) + 1
        if( ischn(0) == ischn(icsu+1) ) then
         goto 98    ! next sumtally sub sction data
        else
         goto 901   ! sumtally start/end is not paired
        end if

       else if( chlw .eq. 'incinc' ) then
        read(isc) jsn        ! include level
        read(isc) dsin(jsn)  ! input file name
        read(isc) idsi(jsn)  ! length of input file name
        read(isc) ill(jsn)   ! input file lines
        read(isc) ilf(jsn)   ! final line to read
        goto 100             ! next Line

       else if( chlw .eq. 'incdec' ) then
        read(isc) jsn        ! include level
        read(isc) dsin(jsn)  ! input file name
        read(isc) idsi(jsn)  ! length of input file name
        read(isc) ill(jsn)   ! input file lines
        read(isc) ilf(jsn)   ! final line to read
        goto 100             ! next Line
       else if( chlw(i1:i1+3) .eq. 'set:' ) then
        goto 100             ! next Line

       end if

  150  continue
       if( ierr .ne. 0 ) goto 9000
       if( jpn  .eq. 3 ) goto 800              ! go to parameters check

!-----------------------------------------------------------------------
!        identify the parameters
!-----------------------------------------------------------------------

       icl = i1
  200  continue
       chlc = chlw
       call chcomp(chlc,icl,i3,i5)
       do i = 1, icsu
        il = icl + lschn(i) - 1
        if( chlc(icl:il) .eq. schan(i)(1:lschn(i)) ) goto 300
       end do

       goto 987

!-----------------------------------------------------------------------
!        read value of parameters
!-----------------------------------------------------------------------

  300  continue

       ipm = i
       ischn( ipm ) = ischn( ipm ) + 1

       ic = inumc(chlw,il+1,i3,'=') + 1
       ic = jnumc(chlw,ic,i3)

       if( ic .gt. i3 ) goto 997

       icl = inumc(chlw,ic,i3,';') - 1

!-----------------------------------------------------------------------
!      isumtally
!-----------------------------------------------------------------------

       if( ipm == 1 ) then
        call onum(chlw,ic,icl,cvvv,ierr)
        if( ierr .ne. 0 ) goto 997
        isumtally(itnmw+1) = nint( cvvv )
        if(isumtally(itnmw+1) < 0 .or. isumtally(itnmw+1) > 3) goto 931

        icl = jnumc(chlw,icl+2,i3)
        if( icl .le. i3 ) goto 150

!-----------------------------------------------------------------------
!      nfile, file name and weight rate
!-----------------------------------------------------------------------

       else if( ipm == 2 ) then
        call onum(chlw,ic,icl,cvvv,ierr)
        if( ierr .ne. 0 ) goto 997

        if( chlw(icl+1:icl+1) .eq. ';' ) goto 993

        nfile(itnmw+1) = nint( cvvv )
        if( nfile(itnmw+1) <= 0 .or. nfile(itnmw+1) > nfilew ) goto 952

!-----------------------------------------------------------------------
  500   continue
        infile = 0

        ic = i3+1                       ! set the value to read next line

        do i=1, nfile(itnmw+1)

         if( ic .gt. i3 ) then
  520     call sumtal_readl(jsn,isc,dsin,idsi,ill,ilf,'#!$',
     &                      jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)
          if( ierr .ne. 0 ) goto 9000
          if( jpn  .eq. 3 ) goto 979
          if( iskip .ne. 0 ) goto 520
          if( chlw .eq. 'incinc' ) then
           read(isc) jsn      ! include level
           read(isc) dsin(jsn) ! input file name
           read(isc) idsi(jsn)  ! length of input file name
           read(isc) ill(jsn)   ! input file lines
           read(isc) ilf(jsn)   ! final line to read
           goto 520             ! next Line

          else if( chlw .eq. 'incdec' ) then
           read(isc) jsn      ! include level
           read(isc) dsin(jsn)  ! input file name
           read(isc) idsi(jsn)  ! length of input file name
           read(isc) ill(jsn)   ! input file lines
           read(isc) ilf(jsn)   ! final line to read
           goto 520             ! next Line
          else if( chlw(i1:i1+3) .eq. 'set:' ) then
           goto 520             ! next Line
          end if
          ic = i1
         end if


         ic = jnumc(chlw,ic,i3)                           ! file name start col
         icf = min( inumc(chlw,ic,i3,' ') - 1, i3 )
         infile = infile + 1
         ltallyfname(infile,itnmw+1) = icf - ic + 1             ! character length
         tallyfname (infile,itnmw+1)(1:icf-ic+1) = chin(ic:icf) ! tally file name(s)
         if ( ltallyfname(infile,itnmw+1) < 1 ) goto 957

!.......................................................................
!     Rename (input) restart file name
!.......................................................................

         ! for [T-Yield] tally
         if ( ital(itnmw+1) ==  3 .and.
     &        ital(itnmw+3) == 16 ) then        ! t-yield, t-dchain
          j = index(tallyfname(infile,itnmw+1), ".", back=.true.)
          if ( j /= 0 ) then
           j = j - 1
          else
           j = ltallyfname(infile,itnmw+1)
          end if
           tallyfname(infile,itnmw+1) =
     &     tallyfname(infile,itnmw+1)(1:j) // '.dyld'
          ltallyfname(infile,itnmw+1) = j + 5

         ! for [T-Track] tally
         else if ( ital(itnmw+1) ==  1 .and.
     &             ital(itnmw+2) == 16 ) then   ! t-track, t-dchain
          j = index(tallyfname(infile,itnmw+1), ".", back=.true.)
          if ( j /= 0 ) then
           j = j - 1
          else
           j = ltallyfname(infile,itnmw+1)
          end if
           tallyfname(infile,itnmw+1) =
     &     tallyfname(infile,itnmw+1)(1:j) // '.dtrk'
          ltallyfname(infile,itnmw+1) = j + 5

         end if

!-----------------------------------------------------------------------

         ic = jnumc(chlw,icf+2,i3)                       ! weight val start col
         call snum(chlw,ic,i3,ic2,cvvv,ierr)
         if( ierr .ne. 0 ) goto 979
         weightRate(infile,itnmw+1) = cvvv               ! weighting rate

         ic = ic2                                        ! next col
        end do

!-----------------------------------------------------------------------
!      sfile
!-----------------------------------------------------------------------

       else if( ipm == 3 ) then
        icf = min( inumc(chlw,ic,icl,' ') - 1, icl )
        lsfile(itnmw+1)             = icf - ic + 1
         sfile(itnmw+1)(1:icf-ic+1) = chin(ic:icf)
        if ( lsfile(itnmw+1) < 1 ) goto 958

        if ( itnmw+1 .gt. 1 ) then
         do iditnmw = 0, itnmw-1
          if ( sfile(itnmw+1)(1:lsfile(itnmw+1))
     &           .eq. sfile(iditnmw+1)(1:lsfile(iditnmw+1)) ) goto 960
         end do
        end if

!.......................................................................
!     Rename (output) restart file name
!.......................................................................

        ! for [T-Yield] tally
        if ( ital(itnmw+1) ==  3 .and.
     &       ital(itnmw+3) == 16 ) then        ! t-yield, t-dchain
         j = index(sfile(itnmw+1), ".", back=.true.)
         if ( j /= 0 ) then
          j = j - 1
         else
          j = lsfile(itnmw+1)
         end if
          sfile(itnmw+1) = sfile(itnmw+1)(1:j) // '.dyld'
         lsfile(itnmw+1) = j + 5

        ! for [T-Track] tally
        else if ( ital(itnmw+1) ==  1 .and.
     &            ital(itnmw+2) == 16 ) then   ! t-track, t-dchain
         j = index(sfile(itnmw+1), ".", back=.true.)
         if ( j /= 0 ) then
          j = j - 1
         else
          j = lsfile(itnmw+1)
         end if
          sfile(itnmw+1) = sfile(itnmw+1)(1:j) // '.dtrk'
         lsfile(itnmw+1) = j + 5

        end if

        icl = jnumc(chlw,icl+2,i3)
        if( icl .le. i3 ) goto 150

!-----------------------------------------------------------------------
!      sumfactor
!-----------------------------------------------------------------------

       else if( ipm == 4 ) then
        call onum(chlw,ic,icl,cvvv,ierr)
        if( ierr .ne. 0 ) goto 997
        sumfactor(itnmw+1) = cvvv
        if ( sumfactor(itnmw+1) <= 0.0d0 ) goto 956

        icl = jnumc(chlw,icl+2,i3)
        if( icl .le. i3 ) goto 150

!-----------------------------------------------------------------------
!      angel parameter
!-----------------------------------------------------------------------

       else if( ipm == 5 ) then
        ict = min( ic + 199, i2 )
         sumang(itnmw+1)(1:ict-ic+1) = chin(ic:ict)
        lsumang(itnmw+1) = ict - ic + 1
        if ( lsumang(itnmw+1) < 1 ) goto 959

       else
       end if

       go to 100   ! read next line

!-----------------------------------------------------------------------
!     check
!-----------------------------------------------------------------------

  800  continue
       if( ischn(0) /= ischn(icsu+1) ) then
        goto 901   ! sumtally start/end is not paired
       end if

!-----------------------------------------------------------------------
! single sumtally (Turn off the following if block, if multiple tally version)

c to use sumtally option without input of isumtally (D=1)
!-----------------------------------------------------------------------

       do irtal=1, itnm

        if ( itfln(irtal).gt.0 ) then
c to use sumtally option in the case a usual tally name is defined.

        else
         if ( icntl.ne.13 ) then
          m_err = 'File is not defined in this tally '
          ErrCha = ''
          ErrID = 'L:1281/R:set_sumtal/F:sumtally.f'
          goto 999
         end if

        end if

        if( nsumtalRead(irtal) > 0 ) then   ! used sumtally function

         ! input tally files
         if( nfile(irtal) <= 0 ) then
          write(m_err,'(a,i3,a)')
     &     'Number of nfiles is little equal 0 in ',irtal,
     &     ' th tally section.'
          ErrCha = ''
          ErrID = 'L:1295/R:set_sumtal/F:sumtally.f'
          goto 999
         end if

         ! When using the weighted average
         if ( isumtally(irtal) == 2 ) then
          ! Sum of the weight rates
          sumWR(irtal) = 0.0d0
          do i = 1, nfile(irtal)
            sumWR(irtal) = sumWR(irtal) + weightRate(i,irtal)
          end do
          if ( sumWR(irtal) == 0.0d0 ) goto 955
         else if ( isumtally(irtal) == 3 ) then ! S.H. added (2020.12.15)
          sumWR(irtal) = 1.0d0
         end if

         ! output tally file name
         if( lsfile(irtal) < 1 ) then
          write(m_err,'(a,i3,a)')
     &     'sumtally output file name is undefined. ',irtal,
     &     ' th tally section'
          ErrCha = ''
          ErrID = 'L:1317/R:set_sumtal/F:sumtally.f'
          goto 999
         end if

        end if

       end do

!-----------------------------------------------------------------------

 9000  continue
       close(isc)
       return

!-----------------------------------------------------------------------
!     errors
!-----------------------------------------------------------------------

  901  m_err = 'sumtally start/end is not paired in '//tname
       ErrCha = ''
       ErrID = 'L:1337/R:set_sumtal/F:sumtally.f'
       goto 999

  931  m_err = 'Isumtally can only be set to 0, 1, 2, or 3 in '//tname
       ErrCha = ''
       ErrID = 'L:1342/R:set_sumtal/F:sumtally.f'
       goto 999

  951  m_err = 'Number of nfiles is little equal 0 in '//tname
       ErrCha = ''
       ErrID = 'L:1347/R:set_sumtal/F:sumtally.f'
       goto 999

  952  m_err = 'nfiles is inconsistent in '//tname
       ErrCha = ''
       ErrID = 'L:1352/R:set_sumtal/F:sumtally.f'
       goto 999

  955  m_err = 'Total weight is zero in '//tname
       ErrCha = ''
       ErrID = 'L:1357/R:set_sumtal/F:sumtally.f'
       goto 999

  956  m_err = 'Sumfactor is inconsistent in  '//tname
       ErrCha = ''
       ErrID = 'L:1362/R:set_sumtal/F:sumtally.f'
       goto 999

  957  m_err = 'Illegal tally file name in '//tname
       ErrCha = ''
       ErrID = 'L:1367/R:set_sumtal/F:sumtally.f'
       goto 999

  958  m_err = 'Illegal sfile name in '//tname
       ErrCha = ''
       ErrID = 'L:1372/R:set_sumtal/F:sumtally.f'
       goto 999

  960  m_err = 'the same sfile name is used in '//tname
       ErrCha = ''
       ErrID = 'L:1377/R:set_sumtal/F:sumtally.f'
       goto 999

  959  m_err = 'Illegal ANGEL parameter in '//tname
       ErrCha = ''
       ErrID = 'L:1382/R:set_sumtal/F:sumtally.f'
       goto 999

  979  m_err = 'Description of parameter is wrong in '//tname
       ErrCha = ''
       ErrID = 'L:1387/R:set_sumtal/F:sumtally.f'
       goto 999

  987  m_err = 'Unknown parameter is found in '//tname
       ErrCha = ''
       ErrID = 'L:1392/R:set_sumtal/F:sumtally.f'
       goto 999

  989  m_err = 'Double definition of the parameter is found in '//tname
       ErrCha = ''
       ErrID = 'L:1397/R:set_sumtal/F:sumtally.f'
       goto 999

  993  m_err = 'In this line, [ ; ] cannot be used.'
       ErrCha = ''
       ErrID = 'L:1402/R:set_sumtal/F:sumtally.f'
       goto 999

  994  m_err = 'In this line, continuation line format cannot used.'
       ErrCha = ''
       ErrID = 'L:1407/R:set_sumtal/F:sumtally.f'
       goto 999

  997  m_err = 'Description of parameter is wrong in '//tname
       ErrCha = ''
       ErrID = 'L:1412/R:set_sumtal/F:sumtally.f'
       goto 999

!-----------------------------------------------------------------------

  999  continue
       l_err = ill(jsn)
       k_err = jsn
       ierr  = 1

*-----------------------------------------------------------------------

 9999  continue
       close(isc)

       ErrCha = ''
       MsgID = 'L:1428/R:set_sumtal/F:sumtally.f'
       call ErrWrite(MsgID, ErrCha)
       write(jo,'(/" ***** Error in Sumtally Subsection *****"/)')
       write(jo,*) dsin(k_err)(1:idsi(k_err)),l_err,':'
       write(io,'(/" ***** Error in Sumtally Subsection *****"/)')
       write(io,*) dsin(k_err)(1:idsi(k_err)),l_err,':'

       icf = 200

       do 910 i = 200, 1, -1

        if( m_err(i:i) .ne. ' ' ) goto 911

 910   continue

 911   icf = i

       call ErrWrite(ErrID, ErrCha)
       write(jo,'(" error = ",200A1/)') ( m_err(i:i), i=1,icf )
       write(io,'(" error = ",200A1/)') ( m_err(i:i), i=1,icf )

       return

      end subroutine set_sumtal


!***********************************************************************
!                                                                      *
      subroutine sumtal_readlt(isc,jsn,jsi,dsin,idsi,ill,ilf,cmmt,
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)
!                                                                      *
!     read lines of input files (text file)                            *
!     original routine is readl in utl3.f                              *
!     Change the file closing process and infl processing              *
!         of the original code.                                        *
!                                                                      *
!***********************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      character m_err*200
      common /error/ m_err, l_err, k_err

*-----------------------------------------------------------------------

      character chin*200, chlw*200, chcm*200
      character chins*200, chlws*200, chcms*200
      parameter ( icolms = 200 )

      character dsin(0:9)*200
      dimension idsi(0:9)

      dimension ill(0:9), ilf(0:9)

      character yen*1

      character*(*) cmmt

      integer   isc            ! output unit number for scrach file

*-----------------------------------------------------------------------

            yen  = char(92)
            isql = 0
            i3s  = 1
            i4s  = 1

*-----------------------------------------------------------------------

            ierr = 0

*-----------------------------------------------------------------------
*     read one line from jsi
*-----------------------------------------------------------------------

  140       ill(jsn) = ill(jsn) + 1

*-----------------------------------------------------------------------

            read(jsi,'(a200)', iostat=ios ) chin
            if( ios .eq. -1 ) goto 251

*-----------------------------------------------------------------------

  240       continue

               if( ill(jsn) .le. ilf(jsn) ) goto 250

*-----------------------------------------------------------------------
*           close files
*-----------------------------------------------------------------------

  251       continue

                  call closef(jsi,jsn)

               write(isc) 'incdec'   ! "incdec" key word
               write(isc) jsn        ! include level
               write(isc) dsin(jsn)  ! input file name
               write(isc) idsi(jsn)  ! length of input file name
               write(isc) ill(jsn)   ! input file lines
               write(isc) ilf(jsn)   ! final line to read

               if( jsn .le. 0 ) then

                  jpn = 3
                  return

               else if( jsn .gt. 0 ) then

                  goto 140

               end if

*-----------------------------------------------------------------------
*        one line edit
*-----------------------------------------------------------------------

  250    continue

                  call chlngt(chin,200,i1,i2)

                  chlw = chin
                  call chcaps(chlw,i1,i2,i3,cmmt)

                  chcm = chlw
                  call chcomp(chcm,i1,i3,i4)

                  k = i1

*-----------------------------------------------------------------------
*        comment line
*        normal line  iskip = 0
*        blank line   iskip = 1
*        comment line iskip = 2
*             ( comment character or 'c ' within 5 column )
*-----------------------------------------------------------------------

            if( i1 .eq. 0 .and. i2 .eq. 0 ) then

               iskip = 1
               return

            else if( index(cmmt,chin(i1:i1)) .ne. 0 ) then

               iskip = 2
               return

            else if( chlw(i1:i1+1) .eq. 'c ' .and. i1 .le. 5  ) then

               iskip = 2
               return

            else

               iskip = 0

            end if

*-----------------------------------------------------------------------
*        sequential line
*-----------------------------------------------------------------------

            if( chlw(i3:i3) .eq. yen ) then

               isql = isql + 1

  150          continue

               if( i3s + i3 - i1 + 1 .gt. 200 ) goto 998

               chin(i3:i3) = ' '
               chlw(i3:i3) = ' '
               chcm(i4:i4) = ' '

               chins(i3s:i3s+i3-i1) = chin(i1:i3)
               chlws(i3s:i3s+i3-i1) = chlw(i1:i3)
               chcms(i4s:i4s+i4-i1) = chlw(i1:i4)

               i3s = i3s + i3 - i1 + 1
               i4s = i4s + i4 - i1 + 1

                  ill(jsn) = ill(jsn) + 1

                  read(jsi,'(a200)', iostat = ios ) chin
                  if( ios .eq. -1 ) goto 999

                  if( ill(jsn) .gt. ilf(jsn) ) goto 999

                  call chlngt(chin,200,i1,i2)

                  chlw = chin
                  call chcaps(chlw,i1,i2,i3,cmmt)

                  chcm = chlw
                  call chcomp(chcm,i1,i3,i4)

                  if( chlw(i3:i3) .eq. yen ) goto 150

               if( i3s + i3 - i1 + 1 .gt. 200 ) goto 998

               chins(i3s:i3s+i3-i1) = chin(i1:i3)
               chlws(i3s:i3s+i3-i1) = chlw(i1:i3)
               chcms(i4s:i4s+i4-i1) = chlw(i1:i4)

               i3s = i3s + i3 - i1 + 1
               i4s = i4s + i4 - i1 + 1

            end if

            if( isql .gt. 0  ) then

               i1 = 1
               i2 = i3s - 1
               i3 = i3s - 1
               i4 = i4s - 1

               chin(i1:i3) = chins(i1:i3)
               chlw(i1:i3) = chlws(i1:i3)
               chcm(i1:i4) = chcms(i1:i4)

               k = i1

            end if

*-----------------------------------------------------------------------
*        infl: {filename} [1-15]  ; include file
*-----------------------------------------------------------------------

            if( jpn .ne. 2. and. chlw(i1:i1+4) .eq. 'infl:' ) then

                  k = k + 4

c  *** T.Sato 2014/2/20, avoid error when -, [, or ] are used in the comment
               do iii=i3+1,i2
                chin(iii:iii)=' '
               enddo
c  *****************************

               call inclf(jsn,jsi,dsin,chin,k,200,ill,ilf,idsi,ierr)

               write(isc) 'incinc'   ! "incinc" key word
               write(isc) jsn        ! include level
               write(isc) dsin(jsn)  ! input file name
               write(isc) idsi(jsn)  ! length of input file name
               write(isc) ill(jsn)   ! input file lines
               write(isc) ilf(jsn)   ! final line to read

                  if( ierr .ne. 0 ) return

                  goto 140

            end if

*-----------------------------------------------------------------------
*       stop reading or stop reading the sections
*-----------------------------------------------------------------------

            if( jpn .le. 1 ) then

               if( chlw(i1:i1+1) .eq. 'q:' ) then

                  jsn0 = jsn

                  do j = jsn0, 1, -1

                     call closef(jsi,jsn)

                  end do

                  jpn = 3
                  return

               else if( chlw(i1:i1+2) .eq. 'qp:' ) then

                  jpn = 2
                  goto 140

               end if

            end if

*-----------------------------------------------------------------------
*       skip lines under qp: up to [sections]
*-----------------------------------------------------------------------

         if( jpn .eq. 2 ) then

            if( i1 .gt. 5 .or. chlw(i1:i1) .ne. '[' ) goto 140

            jpn = 1

         end if

*-----------------------------------------------------------------------

      return

*-----------------------------------------------------------------------

  998 continue

         m_err = 'Sequential line include too more characters: max(200)'
         ErrCha = ''
         ErrID = 'L:1739/R:sumtal_readlt/F:sumtally.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  999 continue

         m_err = 'Sequential line meets the end of file.'
         ErrCha = ''
         ErrID = 'L:1751/R:sumtal_readlt/F:sumtally.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

      end subroutine sumtal_readlt
!***********************************************************************


!***********************************************************************
!                                                                      *
      subroutine sumtal_readl(jsn,jsi,dsin,idsi,ill,ilf,cmmt,
     &                        jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)
!                                                                      *
!     read lines of input files(sumtally sub-section data, binary file)*
!     original routine is readl in utl3.f                              *
!                                                                      *
!***********************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      character m_err*200
      common /error/ m_err, l_err, k_err

*-----------------------------------------------------------------------

      character chin*200, chlw*200, chcm*200
      character chins*200, chlws*200, chcms*200
      parameter ( icolms = 200 )

      character dsin(0:9)*200
      dimension idsi(0:9)

      dimension ill(0:9), ilf(0:9)

      character yen*1

      character*(*) cmmt

*-----------------------------------------------------------------------

            yen  = char(92)
            isql = 0
            i3s  = 1
            i4s  = 1

*-----------------------------------------------------------------------

            ierr = 0

*-----------------------------------------------------------------------
*     read one line from jsi
*-----------------------------------------------------------------------

  140       ill(jsn) = ill(jsn) + 1
*-----------------------------------------------------------------------

            do iii=1,200
               chin(iii:iii)=' '
            enddo
            read(jsi, iostat=ios ) chin
            if( ios .eq. -1 ) then

               jpn = 3
               return

            end if

*-----------------------------------------------------------------------

  240       continue

               if( ill(jsn) .le. ilf(jsn) ) goto 250

*-----------------------------------------------------------------------
*        one line edit
*-----------------------------------------------------------------------

  250    continue

                  call chlngt(chin,200,i1,i2)

                  chlw = chin
                  call chcaps(chlw,i1,i2,i3,cmmt)

                  chcm = chlw
                  call chcomp(chcm,i1,i3,i4)

                  k = i1

*-----------------------------------------------------------------------
*        comment line
*        normal line  iskip = 0
*        blank line   iskip = 1
*        comment line iskip = 2
*             ( comment character or 'c ' within 5 column )
*-----------------------------------------------------------------------

            if( i1 .eq. 0 .and. i2 .eq. 0 ) then

               iskip = 1
               return

            else if( index(cmmt,chin(i1:i1)) .ne. 0 ) then

               iskip = 2
               return

            else if( chlw(i1:i1+1) .eq. 'c ' .and. i1 .le. 5  ) then

               iskip = 2
               return

            else

               iskip = 0

            end if

*-----------------------------------------------------------------------
*        sequential line
*-----------------------------------------------------------------------

            if( chlw(i3:i3) .eq. yen ) then

               isql = isql + 1

  150          continue

               if( i3s + i3 - i1 + 1 .gt. 200 ) goto 998

               chin(i3:i3) = ' '
               chlw(i3:i3) = ' '
               chcm(i4:i4) = ' '

               chins(i3s:i3s+i3-i1) = chin(i1:i3)
               chlws(i3s:i3s+i3-i1) = chlw(i1:i3)
               chcms(i4s:i4s+i4-i1) = chlw(i1:i4)

               i3s = i3s + i3 - i1 + 1
               i4s = i4s + i4 - i1 + 1

                  ill(jsn) = ill(jsn) + 1

                  read(jsi, iostat = ios ) chin
                  if( ios .eq. -1 ) goto 999

                  if( ill(jsn) .gt. ilf(jsn) ) goto 999

                  call chlngt(chin,200,i1,i2)

                  chlw = chin
                  call chcaps(chlw,i1,i2,i3,cmmt)

                  chcm = chlw
                  call chcomp(chcm,i1,i3,i4)

                  if( chlw(i3:i3) .eq. yen ) goto 150

               if( i3s + i3 - i1 + 1 .gt. 200 ) goto 998

               chins(i3s:i3s+i3-i1) = chin(i1:i3)
               chlws(i3s:i3s+i3-i1) = chlw(i1:i3)
               chcms(i4s:i4s+i4-i1) = chlw(i1:i4)

               i3s = i3s + i3 - i1 + 1
               i4s = i4s + i4 - i1 + 1

            end if

            if( isql .gt. 0  ) then

               i1 = 1
               i2 = i3s - 1
               i3 = i3s - 1
               i4 = i4s - 1

               chin(i1:i3) = chins(i1:i3)
               chlw(i1:i3) = chlws(i1:i3)
               chcm(i1:i4) = chcms(i1:i4)

               k = i1

            end if

*-----------------------------------------------------------------------
*       stop reading or stop reading the sections
*-----------------------------------------------------------------------

            if( jpn .le. 1 ) then

               if( chlw(i1:i1+1) .eq. 'q:' ) then

                     jpn = 3
                     return

               else if( chlw(i1:i1+2) .eq. 'qp:' ) then

                     jpn = 2
                     return

               end if

            end if

*-----------------------------------------------------------------------

      return

*-----------------------------------------------------------------------

  998 continue

         m_err = 'Sequential line include too more characters: max(200)'
         ErrCha = ''
         ErrID = 'L:1976/R:sumtal_readl/F:sumtally.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  999 continue

         m_err = 'Sequential line meets the end of file.'
         ErrCha = ''
         ErrID = 'L:1988/R:sumtal_readl/F:sumtally.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

      end subroutine sumtal_readl
!***********************************************************************


!***********************************************************************
!                                                                      *
      subroutine sumtal_echo(iot,m)
!                                                                      *
!     echo sumtally sub section data                                   *
!     create by S.Hashimoto on 2014/11/30                              *
!                                                                      *
!***********************************************************************

      use sumtallymod

      implicit none

*-----------------------------------------------------------------------

      integer iot, m, maxltf, infile, k, icm, ic
      character ch200*400
      character ch100*400

*-----------------------------------------------------------------------

      write(iot,'("sumtally start")')

      write(iot,'("isumtally =",i5,11x,
     &" # (D=1) sumtally option, 1:integration, 2,3:weighted sum")')
     &     isumtally(m)

      write(iot,'("    nfile =",i5,11x,
     &" # number of tally files")')
     &     nfile(m)

      maxltf = 0
      do infile = 1, nfile(m)
         maxltf = max ( maxltf, ltallyfname(infile,m) )
      end do

      if ( maxltf .le. 20 ) then
       write(iot,'("#   file name",15x,"weight")')
       do infile = 1, nfile(m)
          write(iot,'(4x,a20,3x,1p1e12.5)')
     &         tallyfname(infile,m),weightRate(infile,m)
       end do

      else if ( maxltf .le. 30 ) then
       write(iot,'("#   file name",25x,"weight")')
       do infile = 1, nfile(m)
          write(iot,'(4x,a30,3x,1p1e12.5)')
     &         tallyfname(infile,m),weightRate(infile,m)
       end do

      else if ( maxltf .gt. 30 ) then
       write(iot,'("#   file name / weight")')
       do infile = 1, nfile(m)
          write(iot,'(4x,200a1)')
     &         (tallyfname(infile,m)(k:k),k=1,ltallyfname(infile,m))
          write(iot,'(4x,1p1e12.5)') weightRate(infile,m)
       end do

      end if

      icm = 41
      ch100(1:icm) = ' # file name of output by sumtally option'
      do k = 1, 27
         ch200(k:k) = ' '
      end do
      ic = 12
      ch200(1:ic) = '    sfile = '
      ch200(ic+1:ic+lsfile(m)) = sfile(m)(1:lsfile(m))
      if ( ic+lsfile(m) .le. 27 ) then
         ch200(27+1:27+icm) = ch100(1:icm)
         icm = 27+icm
      else
         ch200(ic+lsfile(m)+1:ic+lsfile(m)+icm) = ch100(1:icm)
         icm = ic+lsfile(m)+icm
      end if
      write(iot,'(200a1)') (ch200(k:k),k=1,icm)

      write(iot,'("sumfactor =",1p1e12.5,4x,
     &" # (D=1.0) normalization factor for sumtally")') sumfactor(m)

      if ( lsumang(m) > 0 ) then
       write(iot,'("angel = ",a)') sumang(m)(1:lsumang(m))
      end if

      write(iot,'("sumtally end")')

*-----------------------------------------------------------------------

      return

      end subroutine sumtal_echo
!***********************************************************************


!***********************************************************************
!                                                                      *
      subroutine sumtally(io,jo,ivers,ierr)
!                                                                      *
!     create by T.Miura on 2014/11/30                                  *
!                                                                      *
!***********************************************************************
       use GGBANKMOD
       use TALMOD,   only : tr0
       use sumtallymod
       use TDCHAINMOD, only : pdchreg2
       use partmod, only: itpan, itpat, jtpat ! frtati 2021/10/05

       implicit none

!-----------------------------------------------------------------------

       include 'param.inc'
       include 'err.inc'

!-----------------------------------------------------------------------

       integer   itnm
       integer   ital
       integer   itals
       integer   italm
       common /tall00/ itnm, ital(itlmax), itals(itlmax), italm(itlmax)

       integer itxty
       integer itxnm
       integer itxrg
       double precision rtxmi
       double precision rtxma
       double precision rtxdl
       common /tall04/ itxty(itlmax), itxnm(itlmax), itxrg(itlmax),
     &                 rtxmi(itlmax), rtxma(itlmax), rtxdl(itlmax)

       integer ityty
       integer itynm
       integer ityrg
       double precision rtymi
       double precision rtyma
       double precision rtydl
       common /tall05/ ityty(itlmax), itynm(itlmax), ityrg(itlmax),
     &                 rtymi(itlmax), rtyma(itlmax), rtydl(itlmax)

       integer itzty
       integer itznm
       integer itzrg
       double precision rtzmi
       double precision rtzma
       double precision rtzdl
       common /tall06/ itzty(itlmax), itznm(itlmax), itzrg(itlmax),
     &                 rtzmi(itlmax), rtzma(itlmax), rtzdl(itlmax)

       integer itety
       integer itenm
       integer iterg
       double precision rtemi
       double precision rtema
       double precision rtedl
       common /tall07/ itety(itlmax), itenm(itlmax), iterg(itlmax),
     &                 rtemi(itlmax), rtema(itlmax), rtedl(itlmax)


       integer   itfln
       integer   itfll
       character ctfln*100
       common /tall11/ itfln(itlmax), itfll(itlmax,6), ctfln(itlmax,6)

       integer ittty
       integer ittnm
       integer ittrg
       double precision rttmi
       double precision rttma
       double precision rttdl
       common /tall20/ ittty(itlmax), ittnm(itlmax), ittrg(itlmax),
     &                 rttmi(itlmax), rttma(itlmax), rttdl(itlmax)

       integer itaty
       integer itanm
       integer itarg
       double precision rtami
       double precision rtama
       double precision rtadl
       common /tall23/ itaty(itlmax), itanm(itlmax), itarg(itlmax),
     &                 rtami(itlmax), rtama(itlmax), rtadl(itlmax)

       integer itmcn
       integer itmtn
       integer itmtt
       common /tall22/ itmcn(itlmax), itmtn(itlmax), itmtt(itlmax)

       integer itmlp
       integer itmln
       integer itmst
       integer itmli
       double precision rtmme
       integer itmnt
       integer itmpn
       integer itmpt
       common /tall45/ itmlp(itlmax), itmln(itlmax,6), itmst(itlmax),
     &               itmli(itlmax,6), rtmme(itlmax,6), itmnt(itlmax,6),
     &               itmpn(itlmax,6), itmpt(itlmax,6,6,2)

!-----------------------------------------------------------------------

       integer   itanl
       character itang*200
       common /tall15/ itanl(itlmax), itang(itlmax)

!-----------------------------------------------------------------------

       integer, intent(in)   :: io             !
       integer, intent(in)   :: jo             !
       integer, intent(inout):: ivers          !
       integer, intent(out)  :: ierr           ! return code (=0:normal retrun)

       integer :: i,m
       integer :: iax
       integer :: ncol

!-----------------------------------------------------------------------

       integer           ntf
       integer           np
       integer           ne
       integer           nt
       integer           nx
       integer           ny
       integer           nz
       integer           nm

       integer           na

!-----------------------------------------------------------------------
       ierr = 0

       do m=1, itnm
        if( nsumtalRead(m) > 0 ) goto 8   ! check sumtally sub section
       end do
       return   ! no sumtally sub section

    8  continue

       call ALLOCATE_SUMTAL

       do m=1, itnm

        if( nsumtalRead(m) .lt. 1 ) cycle   ! check sumtally sub section
        if( ital(m) == 16 )         cycle   ! t-dchain

        ntf = 1

!-----------------------------------------------------------------------
!        read tally file.
!-----------------------------------------------------------------------

        call sumtal_read_resfiles(m,io,jo,ivers,ntf,ierr)
        if ( ierr /= 0 ) goto 9999
        call sumtal_read_talls(m,io,ntf,ierr)
        if ( ierr /= 0 ) goto 9999

!-----------------------------------------------------------------------

        if( isumtally(m) == 1 ) then
         call sumtal_calc_stdev(m,ntf,ierr) ! integration

        else if( isumtally(m) == 2 .or. isumtally(m) == 3 ) then
         call sumtal_calc_average(m,ntf,ierr) ! weighted average

        end if
        if ( ierr /= 0 ) goto 9999

        if( nfile(m) > 1 ) then
         do ntf=2, nfile(m)
          call sumtal_read_resfiles(m,io,jo,ivers,ntf,ierr)
          if ( ierr /= 0 ) goto 9999
          call sumtal_read_talls(m,io,ntf,ierr)
          if ( ierr /= 0 ) goto 9999

         if( isumtally(m) == 1 ) then
           call sumtal_calc_stdev(m,ntf,ierr) ! integration
          else if( isumtally(m) == 2 .or. isumtally(m) == 3 ) then
           call sumtal_calc_average(m,ntf,ierr) ! weighted average

          end if
          if ( ierr /= 0 ) goto 9999

         end do   ! ntf(2-nfile) loop end


        end if

       end do   ! m(tally) loop end

       call COPY_SUMTAL

!-----------------------------------------------------------------------

       call ALLOCATE_GGBANK
       call INIT_GGBANK

!-----------------------------------------------------------------------
!        write tally file.
!-----------------------------------------------------------------------

       do m=1, itnm
        if( nsumtalRead(m) .lt. 1 ) cycle   ! check sumtally sub section
        if( ital(m) > 0 ) then
         itfln(m)     = 1                    ! Number of files
         do iax=1, itfln(m)
          itfll(m,iax) = lsfile(m)                  ! file name length
          ctfln(m,iax) =  sfile(m)(1:lsfile(m))     ! file name
         end do
        end if

        ! add to the ANGEL parameter
        if ( itanl(m)+lsumang(m)+1 <= len(itang(m)) ) then
         if ( itanl(m) >= 1 ) then
          itang(m) = itang(m)(1:itanl(m)) // ' ' //
     &              sumang(m)(1:lsumang(m))
          itanl(m) = itanl(m)+lsumang(m)+1
         else
          itang(m) = sumang(m)(1:lsumang(m))
          itanl(m) =  itanl(m)+lsumang(m)
         end if
        else

         ErrCha = ''
         ErrID = 'L:2324/R:sumtally/F:sumtally.f' !E85_001_001
         call ErrWrite(ErrID,ErrCha)

         write(*,'(/" ***** Error in Sumtally Subsection *****"/)')
         write(*,'(" error = ",i3, "-th tally",
     &           " of angel parameter, it could not be added.")') m
         ierr = 1
         go to 9999
        end if
       end do

        ncol = 3
        call talls01(ncol)

!.......................................................................
!      t-dchain tally
!.......................................................................

       do m=1, itnm
        if( nsumtalRead(m) .lt. 1 ) cycle   ! check sumtally sub section
        if( ital (m) == 16 .and.
     &       ( itals(m) == 38 .or. itals(m) == 39
     &       .or. itals(m) == 40 .or. itals(m) == 54 ) ) then ! t-dchain
         call pdchreg2(m)
        end if
       end do

       call DEALLOCATE_GGBANK
       call DEALLOCATE_SUMTAL

!-----------------------------------------------------------------------

 9999  continue
       return

      end subroutine sumtally
!***********************************************************************


!***********************************************************************
!                                                                      *
      subroutine sumtal_read_resfiles(m,io,jo,ivers,ntf,ierr)
!                                                                      *
!     original routine is read_resfiles in restart.f                   *
!     modified by  T.Miura on 2014/11/30                               *
!                                                                      *
!***********************************************************************

        use sumtallymod

!-----------------------------------------------------------------------

        implicit real*8 (a-h,o-z)

        include 'param.inc'
        include 'err.inc'

*-----------------------------------------------------------------------

        common /cparm/  maxbch,maxcas
        common /restart/ resc2(itlmax), resc3(itlmax)
        common /stat / istdev, irestart, ireschk
        common /paraj/  mstz(300), parz(300)

*-----------------------------------------------------------------------

        common /tall00/ itnm, ital(itlmax), itals(itlmax), italm(itlmax)

        common /tall11/ itfln(itlmax), itfll(itlmax,6), ctfln(itlmax,6)
        character ctfln*100

        common /tall59/ irfll(itlmax), crfln(itlmax), itrff(itlmax)
        character crfln*100

*-----------------------------------------------------------------------

        common /res01/ istdevres,maxcasres,rijklstres,irdrf
        common /res02/ crdrfln, irdrfll
        character crdrfln*100

        common /res03/ lextrf(itlmax)
        logical lextrf

        common /res04/ lrijkeqrf
        logical lrijkeqrf

        dimension iristdev(itlmax)
        dimension irmaxcas(itlmax)
        dimension rrijklst(itlmax)

        common /sumtal01/ isistdev, ismaxcas, rsijklst
        integer   :: isistdev(itlmax), ismaxcas(itlmax)
        dimension rsijklst(itlmax)
        character csumfile*200
        integer isjsn, isidsi

*-----------------------------------------------------------------------

        character dsin(0:9)*200
        dimension idsi(0:9)
        dimension ill(0:9), ilf(0:9)

*-----------------------------------------------------------------------

        integer   ntf   ! Number of tallyfname (tally file name)

*-----------------------------------------------------------------------

        iristdev(:) = -1
        irmaxcas(:) = maxcas
        rrijklst(:)  = 0.d0

        lextrf(:)    = .false.

        inoprm = 0

*-----------------------------------------------------------------------
* read each restart files
*-----------------------------------------------------------------------

          jsn     = 0
          jsi     = 32
          idsi(:) = 0
          ill(:)  = 1
          ilf(:)  = 10000000

          ierr = 0

*-----------------------------------------------------------------------
* check the tally is restartable
*-----------------------------------------------------------------------

          if ( .not. any( ital(m) .eq.
     &      (/ 1, 2, 3, 4, 5, 6, 7, 8, 12, 13, 14, 15, 17,
     &         18, 21 /) ) ) then

            goto 100

          end if

*-----------------------------------------------------------------------
* read restart prameters in read resfile
*-----------------------------------------------------------------------

          if ( irfll(m) .eq. 0 ) then
            irfll(m) = itfll(m,1)
            crfln(m)(1:) = ctfln(m,1)(1:itfll(m,1))
          end if


          idsi(jsn) = ltallyfname(ntf,m)

          dsin(jsn)(1:idsi(jsn)) = tallyfname(ntf,m)(1:idsi(jsn))

          isidsi = idsi(jsn)
          csumfile(1:isidsi) = dsin(jsn)(1:idsi(jsn))

          open(jsi,
     &         file=dsin(jsn)(1:idsi(jsn)), status='old', action='read',
     &         form='formatted', iostat=ierr)

          if( ierr .eq. 0 ) then

            call rrestart(jsn,jsi,dsin,idsi,ill,ilf,
     &                    iristdev(m), resc2(m), resc3(m),
     &                    irmaxcas(m), rrijklst(m),
     &                    ierr)

            close(jsi)

            if( ierr .eq. 0 ) then
              lextrf(m) = .true.
            else
              inoprm = inoprm + 1
            end if

          else

            resc2(m) = 0.d0
            resc3(m) = 0.d0

          end if

  100   continue

        ierr = 0

*-----------------------------------------------------------------------
* check exit restart prameters in read resfile
*-----------------------------------------------------------------------

        if( inoprm .gt. 0 ) then
          !! no restart parameter in read reifles
          ErrCha = ''
          MsgID = 'L:2518/R:sumtal_read_resfiles/F:sumtally.f'
          call ErrWrite(MsgID, ErrCha)
          write(jo,
     &      '("Error: Some restart file have no restart parameter.")')
          ierr = 1
          return
        end if

*-----------------------------------------------------------------------
* check exit restart files
*-----------------------------------------------------------------------

        if ( .not. any( lextrf(1:itnm) ) ) then
          !! all restart files are not found
          ErrCha = ''
          MsgID = 'L:2533/R:sumtal_read_resfiles/F:sumtally.f'
          call ErrWrite(MsgID, ErrCha)
          write(jo,'("Error: Tally file(s) written in sumtally ",
     &             "subsection do not exist.")')
          ierr = 1
          return
        end if

*-----------------------------------------------------------------------
* tall number of primery restart file
*-----------------------------------------------------------------------

        irdrf = 0
          if ( lextrf(m) ) then
            irdrf = m
            go to 1000
          end if
 1000     continue


*-----------------------------------------------------------------------
* check all read istdev values are same in one tally
*-----------------------------------------------------------------------

          if ( isumtally(m) .eq. 1 ) then
           if ( ntf .gt. 1 ) then
            if ( isistdev(m) .ne. iristdev(m) ) then
             ErrCha = ''
             MsgID = 'L:2561/R:sumtal_read_resfiles/F:sumtally.f'
             call ErrWrite(MsgID, ErrCha)
             write(jo,'(/"Error: an inconsistent istdev is found "
     &                 ,"in tally file for sumtally.")')
             write(jo,'(9x,a8,i2,a4,a)')
     &            'istdev= ',iristdev(m),' in ',csumfile(1:isidsi)
             write(io,'(/"Error: an inconsistent istdev is found "
     &                 ,"in tally file for sumtally.")')
             write(io,'(9x,a8,i2,a4,a)')
     &            'istdev= ',iristdev(m),' in ',csumfile(1:isidsi)
             ierr = 1
             return
            end if
           end if
          end if
          isistdev(m) = iristdev(m)

*-----------------------------------------------------------------------
* check all read maxcas values are same
*-----------------------------------------------------------------------

          if ( isumtally(m) .eq. 1 ) then
           if ( ntf .gt. 1 ) then
            if ( ismaxcas(m) .ne. irmaxcas(m) ) then
             ErrCha = ''
             MsgID = 'L:2586/R:sumtal_read_resfiles/F:sumtally.f'
             call ErrWrite(MsgID, ErrCha)
             write(jo,'(/"Error: an inconsistent maxcas is found "
     &                 ,"in tally file for sumtally.")')
             write(jo,'(9x,a8,i11,a4,a)')
     &            'maxcas= ',irmaxcas(m),' in ',csumfile(1:isidsi)
             write(io,'(/"Error: an inconsistent maxcas is found "
     &                 ,"in tally file for sumtally.")')
             write(io,'(9x,a8,i11,a4,a)')
     &            'maxcas= ',irmaxcas(m),' in ',csumfile(1:isidsi)
             ierr = 1
             return
            end if
           end if
          end if
          ismaxcas(m) = irmaxcas(m)



          rsijklst(m) = rrijklst(m)

      end subroutine sumtal_read_resfiles
!***********************************************************************


!***********************************************************************
!                                                                      *
      subroutine sumtal_read_talls(m,io,ntf,ierr)
!                                                                      *
!     original routine is read_talls in restart.f                      *
!     modified by  T.Miura on 2014/11/30                               *
!                                                                      *
!***********************************************************************

        use GGBANKMOD
        use sumtallymod

        implicit real*8 (a-h,o-z)

        include 'param.inc'

*-----------------------------------------------------------------------

        common /cparm/  maxbch,maxcas
        common /restart/ resc2(itlmax), resc3(itlmax)
        common /stat / istdev, irestart, ireschk

*-----------------------------------------------------------------------

        common /tall00/ itnm, ital(itlmax), itals(itlmax), italm(itlmax)

        common /tall11/ itfln(itlmax), itfll(itlmax,6), ctfln(itlmax,6)
        character ctfln*100

        common /tall59/ irfll(itlmax), crfln(itlmax), itrff(itlmax)
        character crfln*100

*-----------------------------------------------------------------------

        character dsin(0:9)*200
        dimension idsi(0:9)
        dimension ill(0:9), ilf(0:9)

*-----------------------------------------------------------------------

        integer   ntf   ! Number of tallyfname (tally file name)

*-----------------------------------------------------------------------

        ierr    = 0

        call ALLOCATE_GGBANK
        call INIT_GGBANK

*-----------------------------------------------------------------------
* read each tally files
*-----------------------------------------------------------------------

         if( nsumtalRead(m) == 0 ) return     ! skipped, next do loop

        do 800 iax = 1, 1

*-----------------------------------------------------------------------
*         [t-track]
*-----------------------------------------------------------------------

          if ( ital(m) .eq. 1 ) then

            call sumtal_read_ttrack(m,iax,ntf,ierr)

            if( ierr .ne. 0 ) then
              goto 900
            end if

*-----------------------------------------------------------------------
*         [t-adjoint]
*-----------------------------------------------------------------------

          else if ( ital(m) .eq. 19 ) then

            call sumtal_read_tadjnt(m,iax,ntf,ierr)
            if( ierr .ne. 0 ) then
              goto 900
            end if

*-----------------------------------------------------------------------
*         [t-cross]
*-----------------------------------------------------------------------

          else if ( ital(m) .eq. 2 ) then

            call sumtal_read_tcross(m,iax,ntf,ierr)
            if( ierr .ne. 0 ) then
              goto 900
            end if

*-----------------------------------------------------------------------
*         [t-yield]
*-----------------------------------------------------------------------

          else if ( ital(m) .eq. 3 ) then

            call sumtal_read_tyield(m,iax,ntf,ierr)
            if( ierr .ne. 0 ) then
              goto 900
            end if

*-----------------------------------------------------------------------
*         [t-heat]
*-----------------------------------------------------------------------

          else if ( ital(m) .eq. 4 ) then

            call sumtal_read_theat(m,iax,ntf,ierr)
            if( ierr .ne. 0 ) then
              goto 900
            end if

*-----------------------------------------------------------------------
*         [t-star]
*-----------------------------------------------------------------------

          else if ( ital(m) .eq. 5 ) then

            call sumtal_read_tstar(m,iax,ntf,ierr)
            if( ierr .ne. 0 ) then
              goto 900
            end if

*-----------------------------------------------------------------------
*         [t-time]
*-----------------------------------------------------------------------

          else if ( ital(m) .eq. 6 ) then

            call sumtal_read_ttime(m,iax,ntf,ierr)
            if( ierr .ne. 0 ) then
              goto 900
            end if

*-----------------------------------------------------------------------
*         [t-dpa]
*-----------------------------------------------------------------------

          else if ( ital(m) .eq. 7 ) then

            call sumtal_read_tdpa(m,iax,ntf,ierr)
            if( ierr .ne. 0 ) then
              goto 900
            end if

*-----------------------------------------------------------------------
*         [t-product]
*-----------------------------------------------------------------------

          else if ( ital(m) .eq. 8 ) then

            call sumtal_read_tproduct(m,iax,ntf,ierr)
            if( ierr .ne. 0 ) then
              goto 900
            end if

*-----------------------------------------------------------------------
*         [t-let]
*-----------------------------------------------------------------------

          else if ( ital(m) .eq. 12 ) then

            call sumtal_read_tlet(m,iax,ntf,ierr)
            if( ierr .ne. 0 ) then
              goto 900
            end if

*-----------------------------------------------------------------------
*         [t-deposit]
*-----------------------------------------------------------------------

          else if ( ital(m) .eq. 13 ) then

            call sumtal_read_tdeposit(m,iax,ntf,ierr)
            if( ierr .ne. 0 ) then
              goto 900
            end if

*-----------------------------------------------------------------------
*         [t-deposit2]
*-----------------------------------------------------------------------

          else if ( ital(m) .eq. 14 ) then

            call sumtal_read_tdeposit2(m,iax,ntf,ierr)
            if( ierr .ne. 0 ) then
              goto 900
            end if

*-----------------------------------------------------------------------
*         [t-sed]
*-----------------------------------------------------------------------

          else if ( ital(m) .eq. 15 ) then

            call sumtal_read_tsed(m,iax,ntf,ierr)
            if( ierr .ne. 0 ) then
              goto 900
            end if

*-----------------------------------------------------------------------
*         [t-point]
*-----------------------------------------------------------------------
          else if ( ital(m) .eq. 17 ) then

            call sumtal_read_tpoint(m,iax,ntf,ierr)

            if( ierr .ne. 0 ) then
              goto 900
            end if

*-----------------------------------------------------------------------
*         [t-wwg]
*-----------------------------------------------------------------------
          else if ( ital(m) .eq. 18 ) then

            call sumtal_read_twwg(m,iax,ntf,ierr)

            if( ierr .ne. 0 ) then
              goto 900
            end if

*-----------------------------------------------------------------------
*         [t-volume]
*-----------------------------------------------------------------------
          else if ( ital(m) .eq. 21 ) then

            call sumtal_read_tvlm(m,iax,ntf,ierr)

            if( ierr .ne. 0 ) then
              goto 900
            end if

*-----------------------------------------------------------------------
          end if

 800    continue

        call trans_trRES2tr0()

 900    continue
        call DEALLOCATE_GGBANK

      end subroutine sumtal_read_talls
!***********************************************************************


!***********************************************************************
!                                                                      *
      subroutine sumtal_open_resfile(m,noe,ntf,
     &                               jsn,jsi,dsin,idsi,ill,ilf,
     &                               newtall,ierr)
!                                                                      *
!  subroutine open resfiles *.out, *_err.out.                          *
!                                                                      *
!     original routine is open_resfile in resutl.f                     *
!                                                                      *
!***********************************************************************

        use sumtallymod

        implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

        include 'param.inc'
        include 'err.inc'

*-----------------------------------------------------------------------

        common /tall59/ irfll(itlmax), crfln(itlmax), itrff(itlmax)
        character crfln*100

*-----------------------------------------------------------------------

        integer   ntf   ! Number of tallyfname (tally file name)

*-----------------------------------------------------------------------

        character dsin(0:9,2)*200
        dimension idsi(0:9,2)
        dimension ill(0:9,2), ilf(0:9,2)
        dimension jsn(2), jsi(2)
        dimension ierrs(2)

*-----------------------------------------------------------------------

        ierr = 0
        newtall = 0

*-----------------------------------------------------------------------
*   open restart file
*-----------------------------------------------------------------------

        do ioe = 1, noe

          jsn(ioe)    = 0
          jsi(ioe)    = 31 + ioe
          idsi(:,ioe) = 0
          ill(:,ioe)  = 1
          ilf(:,ioe)  = 10000000
          ierrs(ioe)  = 0

          jsni = jsn(ioe)
          jsii = jsi(ioe)

          if ( ioe .eq. 1 ) then
            idsi(jsni,ioe) = ltallyfname(ntf,m)
            dsin(jsni,ioe)(1:idsi(jsni,ioe))
     &                     = tallyfname(ntf,m)(1:idsi(jsni,ioe))
          else
            idsi(jsni,ioe) = ltallyfname(ntf,m) + 4
            call mk_2derrfn(tallyfname(ntf,m),dsin(jsni,ioe),
     &                     ltallyfname(ntf,m))
          end if

          open(jsii,
     &         file=dsin(jsni,ioe)(1:idsi(jsni,ioe)),
     &         status='old', action='read', form='formatted',
     &         iostat=ierrs(ioe))

        end do

*-----------------------------------------------------------------------
*   check dose restart file exit
*-----------------------------------------------------------------------

        if( all(ierrs(1:noe) .ne. 0) .and. itrff(m) .eq. 0 ) then
          !! it's newly tally
          ierr = 0
          newtall = 1
          goto 900
        end if

        if( noe .eq. 1 .and. ierrs(1) .ne. 0 .and. itrff(m) .ne. 0 )then
          ierr = 1
       write(ErrCha,'("Error: resfile = ", a, " does not exist.")')
     &          dsin(jsn(1),1)(1:idsi(jsn(1),1))
       ErrID = 'L:2950/R:sumtal_open_resfile/F:sumtally.f' !E85_002_001
       call ErrWrite(ErrID,ErrCha)

       write(ErrCha,'("If you want to create new tally",
     &              " in the restart calculation,",
     &              " you should not specify resfile parameter.")')
       ErrID = 'L:2956/R:sumtal_open_resfile/F:sumtally.f'
       call ErrWrite(ErrID,ErrCha)

          goto 900
        end if

        if( noe .eq. 2 .and. any(ierrs(1:noe) .ne. 0 ) ) then
          ierr = 1
      write(ErrCha,'("Error: Both ", a," and ", a," should exist",
     &              " when you want to restart the calculation.")')
     &          dsin(jsn(1),1)(1:idsi(jsn(1),1)),
     &          dsin(jsn(2),2)(1:idsi(jsn(2),2))
       ErrID = 'L:2968/R:sumtal_open_resfile/F:sumtally.f' !E85_003_001
       call ErrWrite(ErrID,ErrCha)

          goto 900
        end if

*-----------------------------------------------------------------------
  900   continue

        return

      end subroutine sumtal_open_resfile
!***********************************************************************


!***********************************************************************
!                                                                      *
      subroutine sumtal_calc_stdev(m,ntf,ierr)
!                                                                      *
!***********************************************************************

      use TALMOD,   only : tr0
      use RESTALMOD,   only : irestalm
      use sumtallymod
      use partmod, only: itpan, itpat, jtpat ! frtati 2021/10/05
      implicit none

!-----------------------------------------------------------------------

      include 'param.inc'

!-----------------------------------------------------------------------

      integer itnm
      integer ital
      integer itals
      integer italm
      common /tall00/ itnm, ital(itlmax), itals(itlmax), italm(itlmax)

      integer itmsh
      integer itunt
      integer itspc
      integer itout
      integer ittwo
      common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &     itout(itlmax), ittwo(itlmax)

      integer itrgn
      integer itrgm
      integer itreg
      common /tall02/ itrgn(itlmax), itrgm(itlmax), itreg(itlmax)

      integer itrty
      integer itrnm
      integer itrrg
      double precision rtrmi
      double precision rtrma
      double precision rtrdl
      common /tall03/ itrty(itlmax), itrnm(itlmax), itrrg(itlmax),
     &     rtrmi(itlmax), rtrma(itlmax), rtrdl(itlmax)

      integer itxty
      integer itxnm
      integer itxrg
      double precision rtxmi
      double precision rtxma
      double precision rtxdl
      common /tall04/ itxty(itlmax), itxnm(itlmax), itxrg(itlmax),
     &     rtxmi(itlmax), rtxma(itlmax), rtxdl(itlmax)

      integer ityty
      integer itynm
      integer ityrg
      double precision rtymi
      double precision rtyma
      double precision rtydl
      common /tall05/ ityty(itlmax), itynm(itlmax), ityrg(itlmax),
     &     rtymi(itlmax), rtyma(itlmax), rtydl(itlmax)

      integer itzty
      integer itznm
      integer itzrg
      double precision rtzmi
      double precision rtzma
      double precision rtzdl
      common /tall06/ itzty(itlmax), itznm(itlmax), itzrg(itlmax),
     &     rtzmi(itlmax), rtzma(itlmax), rtzdl(itlmax)

      integer itety
      integer itenm
      integer iterg
      double precision rtemi
      double precision rtema
      double precision rtedl
      common /tall07/ itety(itlmax), itenm(itlmax), iterg(itlmax),
     &     rtemi(itlmax), rtema(itlmax), rtedl(itlmax)


      integer itman
      integer itmat
      integer itmct
      integer itnun
      integer itnuc
      integer itndz
      integer itndn
      integer itnkz
      integer itnkn
      common /tall12/ itman(itlmax), itmat(itlmax), itmct(itlmax),
     &                itnun(itlmax), itnuc(itlmax),
     &                itndz(itlmax), itndn(itlmax),
     &                itnkz(itlmax), itnkn(itlmax)

      integer itrcn
      integer itrcr
      integer itrca
      common /tall13/ itrcn(itlmax), itrcr(itlmax), itrca(itlmax)

      integer itndy
      common /tall17/ itndy(itlmax)

      integer ittty
      integer ittnm
      integer ittrg
      double precision rttmi
      double precision rttma
      double precision rttdl
      common /tall20/ ittty(itlmax), ittnm(itlmax), ittrg(itlmax),
     &     rttmi(itlmax), rttma(itlmax), rttdl(itlmax)

      integer itaty
      integer itanm
      integer itarg
      double precision rtami
      double precision rtama
      double precision rtadl
      common /tall23/ itaty(itlmax), itanm(itlmax), itarg(itlmax),
     &                rtami(itlmax), rtama(itlmax), rtadl(itlmax)

      integer itmlp
      integer itmln
      integer itmst
      integer itmli
      double precision rtmme
      integer itmnt
      integer itmpn
      integer itmpt
      common /tall45/ itmlp(itlmax), itmln(itlmax,6), itmst(itlmax),
     &     itmli(itlmax,6), rtmme(itlmax,6), itmnt(itlmax,6),
     &     itmpn(itlmax,6), itmpt(itlmax,6,6,2)

      integer itety2
      integer itenm2
      integer iterg2
      double precision rtemi2
      double precision rtema2
      double precision rtedl2
      common /tall56/ itety2(itlmax), itenm2(itlmax), iterg2(itlmax),
     &                rtemi2(itlmax), rtema2(itlmax), rtedl2(itlmax)

      integer itactnm
      integer itactrg
      integer itactmax
      common /tall72/ itactnm(itlmax), itactrg(itlmax), itactmax(itlmax)

      integer itnzn
      integer itndm
      common /tall83/ itnzn(itlmax), itndm(itlmax)

!-----------------------------------------------------------------------

       integer, intent(out)  :: ierr           ! return code (=0:normal retrun)

       integer :: m

!-----------------------------------------------------------------------

       integer           ntf
       integer           np
       integer           ne, nei
       integer           nt
       integer           nx
       integer           ny
       integer           nz
       integer           nr
       integer           nact      ! S.Abe 2018/02/15
       integer           nm
       integer           na
       integer           nd
       integer           mn
       integer           mz
       integer           mm ! frtati 2022/02/18
       integer           ne1, ne2
       integer           idasz

!-----------------------------------------------------------------------

       ierr = 0

!-----------------------------------------------------------------------

*-----------------------------------------------------------------------
*         [t-track]
*-----------------------------------------------------------------------

      if ( ital(m) .eq. 1 ) then

* region mesh
       if( itmsh(m) .eq. 1 ) then

        np = itpan(m)
        nr = itrgn(m)
        ne = itenm(m)
        nm = itmst(m)
        nt = ittnm(m)

        call sumtal_calc_stdev_tracreg(m,ntf,
     &       np,ne,nt,nr,nm,
     &       sumfactor(m),weightRate(1,m),
     &       tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &       resc2SUMTAL(m),resc3SUMTAL(m),
     &       ierr)              ! integration

* r-z mesh
       else if ( itmsh(m) .eq. 2 ) then

        np = itpan(m)
        nr = itrnm(m)
        nz = itznm(m)
        ne = itenm(m)
        nm = itmst(m)
        nt = ittnm(m)
        na = itanm(m)

        call sumtal_calc_stdev_tracrz(m,ntf,
     &       np,ne,nt,na,nr,nz,nm,
     &       sumfactor(m),weightRate(1,m),
     &       tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &       resc2SUMTAL(m),resc3SUMTAL(m),
     &       ierr)              ! integration

* xyz mesh
       else if ( itmsh(m) .eq. 3 ) then

        np = itpan(m)
        nx = itxnm(m)
        ny = itynm(m)
        nz = itznm(m)
        ne = itenm(m)
        nm = itmst(m)
        nt = ittnm(m)

        call sumtal_calc_stdev_tracxyz(m,ntf,
     &       np,ne,nt,nx,ny,nz,nm,
     &       sumfactor(m),weightRate(1,m),
     &       tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &       resc2SUMTAL(m),resc3SUMTAL(m),
     &       ierr)              ! integration

* tetra mesh !FURUTA20190110
       else if( itmsh(m) .eq. 4 ) then
! mesh=tet can be treated same as mesh=reg
        np = itpan(m)
        nr = itrgn(m)
        ne = itenm(m)
        nm = itmst(m)
        nt = ittnm(m)

        call sumtal_calc_stdev_tracreg(m,ntf,
     &       np,ne,nt,nr,nm,
     &       sumfactor(m),weightRate(1,m),
     &       tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &       resc2SUMTAL(m),resc3SUMTAL(m),
     &       ierr)              ! integration

       end if

*-----------------------------------------------------------------------
*         [t-adjoint]
*-----------------------------------------------------------------------

      else if ( ital(m) .eq. 19 ) then

* region mesh
       if( itmsh(m) .eq. 1 ) then

        np = itpan(m)
        nr = itrgn(m)
        ne = itenm(m)
        nm = itmst(m)
        nt = ittnm(m)

        call sumtal_calc_stdev_tadjreg(m,ntf,
     &       np,ne,nt,nr,nm,
     &       sumfactor(m),weightRate(1,m),
     &       tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &       resc2SUMTAL(m),resc3SUMTAL(m),
     &       ierr)              ! integration

* r-z mesh
       else if ( itmsh(m) .eq. 2 ) then

        np = itpan(m)
        nr = itrnm(m)
        nz = itznm(m)
        ne = itenm(m)
        nm = itmst(m)
        nt = ittnm(m)
        na = itanm(m)

        call sumtal_calc_stdev_tadjrz(m,ntf,
     &       np,ne,nt,na,nr,nz,nm,
     &       sumfactor(m),weightRate(1,m),
     &       tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &       resc2SUMTAL(m),resc3SUMTAL(m),
     &       ierr)              ! integration

* xyz mesh
       else if ( itmsh(m) .eq. 3 ) then

        np = itpan(m)
        nx = itxnm(m)
        ny = itynm(m)
        nz = itznm(m)
        ne = itenm(m)
        nm = itmst(m)
        nt = ittnm(m)

        call sumtal_calc_stdev_tadjxyz(m,ntf,
     &       np,ne,nt,nx,ny,nz,nm,
     &       sumfactor(m),weightRate(1,m),
     &       tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &       resc2SUMTAL(m),resc3SUMTAL(m),
     &       ierr)              ! integration

       end if

*-----------------------------------------------------------------------
*         [t-cross]
*-----------------------------------------------------------------------
      else if ( ital(m) .eq. 2 ) then

* region mesh
       if( itmsh(m) .eq. 1 ) then

        np = itpan(m)
        ne = itenm(m)
        nt = ittnm(m)
        na = itanm(m)
        nr = itrcn(m)
        nm = itmst(m)

        call sumtal_calc_stdev_crsreg(m,ntf,
     &                              np,ne,nt,na,nr,nm,
     &                              sumfactor(m),weightRate(1,m),
     &       tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &                              resc2SUMTAL(m),resc3SUMTAL(m),
     &                              ierr)

* r-z mesh
       else if ( itmsh(m) .eq. 2 ) then

        np = itpan(m)
        ne = itenm(m)
        na = itanm(m)
        nt = ittnm(m)
        nr = itrnm(m)
        nz = itznm(m)
        nm = itmst(m)

        idasz = itenm(m) * itpan(m) * 2
     &        * ( itrnm(m) + 1 ) * itznm(m)
     &        * itanm(m)
     &        * ittnm(m)

        call sumtal_calc_stdev_crsrz(m,ntf,
     &                             np,ne,na,nt,nr,nz,nm,
     &                             sumfactor(m),weightRate(1,m),
     &       tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &       tr0(irestalm(m)+idasz),trSUMTAL(irestalm(m)+idasz),
     &                             resc2SUMTAL(m),resc3SUMTAL(m),
     &                             ierr)

* xyz mesh
       else if ( itmsh(m) .eq. 3 ) then

        nt = ittnm(m)
        np = itpan(m)
        na = itanm(m)
        ne = itenm(m)
        ny = itynm(m)
        nx = itxnm(m)
        nz = itznm(m)
        nm = itmst(m)

        call sumtal_calc_stdev_crsxyz(m,ntf,
     &                             nt,np,na,ne,ny,nx,nz,nm,
     &                             sumfactor(m),weightRate(1,m),
     &       tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &                             resc2SUMTAL(m),resc3SUMTAL(m),
     &                             ierr)

       end if

*-----------------------------------------------------------------------
*         [t-yield]
*-----------------------------------------------------------------------

      else if ( ital(m) .eq. 3 ) then

* region mesh
       if( itmsh(m) .eq. 1 ) then

        mn = itndn(m)
        mz = itndz(m)
        mm = itndm(m)
        nr = itrgn(m)

        call sumtal_calc_stdev_yieldreg(m,ntf,
     &                             mn,mz,mm,nr,
     &                             sumfactor(m),weightRate(1,m),
     &       tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &                             resc2SUMTAL(m),resc3SUMTAL(m),
     &                             ierr)

* r-z mesh
       else if ( itmsh(m) .eq. 2 ) then

        mn = itndn(m)
        mz = itndz(m)
        mm = itndm(m)
        nz = itznm(m)
        nr = itrnm(m)

        call sumtal_calc_stdev_yieldrz(m,ntf,
     &                             mn,mz,mm,nz,nr,
     &                             sumfactor(m),weightRate(1,m),
     &       tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &                             resc2SUMTAL(m),resc3SUMTAL(m),
     &                             ierr)

* xyz mesh
       else if ( itmsh(m) .eq. 3 ) then

        mn = itndn(m)
        mz = itndz(m)
        mm = itndm(m)
        nz = itznm(m)
        ny = itynm(m)
        nx = itxnm(m)

        call sumtal_calc_stdev_yieldxyz(m,ntf,
     &                             mn,mz,mm,nz,ny,nx,
     &                             sumfactor(m),weightRate(1,m),
     &       tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &                             resc2SUMTAL(m),resc3SUMTAL(m),
     &                             ierr)

* tetra mesh !FURUTA20190121
       else if( itmsh(m) .eq. 4 ) then
! mesh=tet can be treated same as mesh=reg
        mn = itndn(m)
        mz = itndz(m)
        mm = itndm(m)
        nr = itrgn(m)

        call sumtal_calc_stdev_yieldreg(m,ntf,
     &                             mn,mz,mm,nr,
     &                             sumfactor(m),weightRate(1,m),
     &       tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &                             resc2SUMTAL(m),resc3SUMTAL(m),
     &                             ierr)

       end if

*-----------------------------------------------------------------------
*         [t-heat]
*-----------------------------------------------------------------------

      else if ( ital(m) .eq. 4 ) then

* region mesh
       if( itmsh(m) .eq. 1 ) then

        ne = itenm(m)
        nd = itndy(m)
        nr = itrgn(m)

        call sumtal_calc_stdev_hetreg(m,ntf,
     &                             ne,nd,nr,
     &                             sumfactor(m),weightRate(1,m),
     &       tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &                             resc2SUMTAL(m),resc3SUMTAL(m),
     &                             ierr)

* r-z mesh
       else if ( itmsh(m) .eq. 2 ) then

        ne = itenm(m)
        nd = itndy(m)
        nz = itznm(m)
        nr = itrnm(m)

        call sumtal_calc_stdev_hetrz(m,ntf,
     &                             ne,nd,nz,nr,
     &                             sumfactor(m),weightRate(1,m),
     &       tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &                             resc2SUMTAL(m),resc3SUMTAL(m),
     &                             ierr)

* xyz mesh
       else if ( itmsh(m) .eq. 3 ) then

        ne = itenm(m)
        nd = itndy(m)
        nz = itznm(m)
        ny = itynm(m)
        nx = itxnm(m)

        call sumtal_calc_stdev_hetxyz(m,ntf,
     &                             ne,nd,nz,ny,nx,
     &                             sumfactor(m),weightRate(1,m),
     &       tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &                             resc2SUMTAL(m),resc3SUMTAL(m),
     &                             ierr)

       end if

*-----------------------------------------------------------------------
*         [t-star]
*-----------------------------------------------------------------------

      else if ( ital(m) .eq. 5 ) then

* region mesh
       if( itmsh(m) .eq. 1 ) then

        np = itpan(m)
        ne = itenm(m)
        nt = ittnm(m)
        nr = itrgn(m)
        nact = itactnm(m)      ! S.Abe 2018/02/15

        call sumtal_calc_stdev_starreg(m,ntf,
     &                             np,ne,nt,nr,nact,
     &                             sumfactor(m),weightRate(1,m),
     &       tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &                             resc2SUMTAL(m),resc3SUMTAL(m),
     &                             ierr)

* r-z mesh
       else if ( itmsh(m) .eq. 2 ) then

        np = itpan(m)
        ne = itenm(m)
        nt = ittnm(m)
        nr = itrnm(m)
        nz = itznm(m)

        call sumtal_calc_stdev_starrz(m,ntf,
     &                             np,ne,nt,nr,nz,
     &                             sumfactor(m),weightRate(1,m),
     &       tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &                             resc2SUMTAL(m),resc3SUMTAL(m),
     &                             ierr)

* xyz mesh
       else if ( itmsh(m) .eq. 3 ) then

        np = itpan(m)
        ne = itenm(m)
        nt = ittnm(m)
        nx = itxnm(m)
        ny = itynm(m)
        nz = itznm(m)

        call sumtal_calc_stdev_starxyz(m,ntf,
     &                             np,ne,nt,nx,ny,nz,
     &                             sumfactor(m),weightRate(1,m),
     &       tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &                             resc2SUMTAL(m),resc3SUMTAL(m),
     &                             ierr)

       end if

*-----------------------------------------------------------------------
*         [t-time]
*-----------------------------------------------------------------------

      else if ( ital(m) .eq. 6 ) then

* region mesh
       if( itmsh(m) .eq. 1 ) then

        np = itpan(m)
        nt = ittnm(m)
        ne = itenm(m)
        nr = itrgn(m)

        call sumtal_calc_stdev_timereg(m,ntf,
     &                             np,nt,ne,nr,
     &                             sumfactor(m),weightRate(1,m),
     &       tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &                             resc2SUMTAL(m),resc3SUMTAL(m),
     &                             ierr)

* r-z mesh
       else if ( itmsh(m) .eq. 2 ) then

        np = itpan(m)
        nt = ittnm(m)
        ne = itenm(m)
        nr = itrnm(m)
        nz = itznm(m)

        call sumtal_calc_stdev_timerz(m,ntf,
     &                             np,nt,ne,nr,nz,
     &                             sumfactor(m),weightRate(1,m),
     &       tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &                             resc2SUMTAL(m),resc3SUMTAL(m),
     &                             ierr)

* xyz mesh
       else if ( itmsh(m) .eq. 3 ) then

        np = itpan(m)
        nt = ittnm(m)
        ne = itenm(m)
        nx = itxnm(m)
        ny = itynm(m)
        nz = itznm(m)

        call sumtal_calc_stdev_timexyz(m,ntf,
     &                             np,nt,ne,nx,ny,nz,
     &                             sumfactor(m),weightRate(1,m),
     &       tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &                             resc2SUMTAL(m),resc3SUMTAL(m),
     &                             ierr)

       end if

*-----------------------------------------------------------------------
*         [t-dpa]
*-----------------------------------------------------------------------

      else if ( ital(m) .eq. 7 ) then

* region mesh
       if( itmsh(m) .eq. 1 ) then

        nd = itndy(m)
        np = itpan(m)
        nr = itrgn(m)

        call sumtal_calc_stdev_dpareg(m,ntf,
     &                             nd,np,nr,
     &                             sumfactor(m),weightRate(1,m),
     &       tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &                             resc2SUMTAL(m),resc3SUMTAL(m),
     &                             ierr)

* r-z mesh
       else if ( itmsh(m) .eq. 2 ) then

        nd = itndy(m)
        np = itpan(m)
        nr = itrnm(m)
        nz = itznm(m)

        call sumtal_calc_stdev_dparz(m,ntf,
     &                             nd,np,nr,nz,
     &                             sumfactor(m),weightRate(1,m),
     &       tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &                             resc2SUMTAL(m),resc3SUMTAL(m),
     &                             ierr)

* xyz mesh
       else if ( itmsh(m) .eq. 3 ) then

        nd = itndy(m)
        np = itpan(m)
        nx = itxnm(m)
        ny = itynm(m)
        nz = itznm(m)

        call sumtal_calc_stdev_dpaxyz(m,ntf,
     &                             nd,np,nx,ny,nz,
     &                             sumfactor(m),weightRate(1,m),
     &       tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &                             resc2SUMTAL(m),resc3SUMTAL(m),
     &                             ierr)

* tetra mesh !FURUTA20190121
       else if( itmsh(m) .eq. 4 ) then
! mesh=tet can be treated same as mesh=reg
        nd = itndy(m)
        np = itpan(m)
        nr = itrgn(m)

        call sumtal_calc_stdev_dpareg(m,ntf,
     &                             nd,np,nr,
     &                             sumfactor(m),weightRate(1,m),
     &       tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &                             resc2SUMTAL(m),resc3SUMTAL(m),
     &                             ierr)

       end if

*-----------------------------------------------------------------------
*         [t-product]
*-----------------------------------------------------------------------

      else if ( ital(m) .eq. 8 ) then

* region mesh
       if( itmsh(m) .eq. 1 ) then

        np = itpan(m)
        na = itanm(m)
        nt = ittnm(m)
        ne = itenm(m)
        nr = itrgn(m)

        call sumtal_calc_stdev_productreg(m,ntf,
     &                             np,na,nt,ne,nr,
     &                             sumfactor(m),weightRate(1,m),
     &       tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &                             resc2SUMTAL(m),resc3SUMTAL(m),
     &                             ierr)

* r-z mesh
       else if ( itmsh(m) .eq. 2 ) then

        np = itpan(m)
        na = itanm(m)
        nt = ittnm(m)
        ne = itenm(m)
        nr = itrnm(m)
        nz = itznm(m)

        call sumtal_calc_stdev_productrz(m,ntf,
     &                             np,na,nt,ne,nr,nz,
     &                             sumfactor(m),weightRate(1,m),
     &       tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &                             resc2SUMTAL(m),resc3SUMTAL(m),
     &                             ierr)

* xyz mesh
       else if ( itmsh(m) .eq. 3 ) then

        np = itpan(m)
        na = itanm(m)
        nt = ittnm(m)
        ne = itenm(m)
        nx = itxnm(m)
        ny = itynm(m)
        nz = itznm(m)

        call sumtal_calc_stdev_productxyz(m,ntf,
     &                             np,na,nt,ne,nx,ny,nz,
     &                             sumfactor(m),weightRate(1,m),
     &       tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &                             resc2SUMTAL(m),resc3SUMTAL(m),
     &                             ierr)

* tetra mesh !FURUTA20190121
       else if( itmsh(m) .eq. 4 ) then
! mesh=tet can be treated same as mesh=reg
        np = itpan(m)
        na = itanm(m)
        nt = ittnm(m)
        ne = itenm(m)
        nr = itrgn(m)

        call sumtal_calc_stdev_productreg(m,ntf,
     &                             np,na,nt,ne,nr,
     &                             sumfactor(m),weightRate(1,m),
     &       tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &                             resc2SUMTAL(m),resc3SUMTAL(m),
     &                             ierr)

       end if

*-----------------------------------------------------------------------
*         [t-let]
*-----------------------------------------------------------------------

      else if ( ital(m) .eq. 12 ) then

* region mesh
       if( itmsh(m) .eq. 1 ) then

        np = itpan(m)
        ne = itenm(m)
        nr = itrgn(m)

        call sumtal_calc_stdev_letreg(m,ntf,
     &                             np,ne,nr,
     &                             sumfactor(m),weightRate(1,m),
     &       tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &                             resc2SUMTAL(m),resc3SUMTAL(m),
     &                             ierr)

* r-z mesh
       else if ( itmsh(m) .eq. 2 ) then

        np = itpan(m)
        ne = itenm(m)
        nr = itrnm(m)
        nz = itznm(m)

        call sumtal_calc_stdev_letrz(m,ntf,
     &                             np,ne,nr,nz,
     &                             sumfactor(m),weightRate(1,m),
     &       tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &                             resc2SUMTAL(m),resc3SUMTAL(m),
     &                             ierr)

* xyz mesh
       else if ( itmsh(m) .eq. 3 ) then

        np = itpan(m)
        ne = itenm(m)
        nx = itxnm(m)
        ny = itynm(m)
        nz = itznm(m)

        call sumtal_calc_stdev_letxyz(m,ntf,
     &                             np,ne,nx,ny,nz,
     &                             sumfactor(m),weightRate(1,m),
     &       tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &                             resc2SUMTAL(m),resc3SUMTAL(m),
     &                             ierr)

       end if

*-----------------------------------------------------------------------
*         [t-deposit]
*-----------------------------------------------------------------------

      else if ( ital(m) .eq. 13 ) then

       if ( itout(m).eq.1 ) then   ! output=Dose

* region mesh
        if( itmsh(m) .eq. 1 ) then

         np = itpan(m)
         nr = itrgn(m)
         ne = itenm(m)
         nt = ittnm(m)

         if( itout(m) .le. 1 ) then
           nei = 0
         else
           nei = 1
         end if

         call sumtal_calc_stdev_depstreg(m,ntf,
     &        np,nei,ne,nt,nr,
     &        sumfactor(m),weightRate(1,m),
     &        tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &        resc2SUMTAL(m),resc3SUMTAL(m),
     &        ierr)              ! integration

* r-z mesh
        else if ( itmsh(m) .eq. 2 ) then

         np = itpan(m)
         ne = itenm(m)
         nt = ittnm(m)
         nr = itrnm(m)
         nz = itznm(m)

         if( itout(m) .le. 1 ) then
           nei = 0
         else
           nei = 1
         end if

         call sumtal_calc_stdev_depstrz(m,ntf,
     &        np,nei,ne,nt,nr,nz,
     &        sumfactor(m),weightRate(1,m),
     &        tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &        resc2SUMTAL(m),resc3SUMTAL(m),
     &        ierr)              ! integration

* xyz mesh
        else if ( itmsh(m) .eq. 3 ) then

         np = itpan(m)
         ne = itenm(m)
         nt = ittnm(m)
         nx = itxnm(m)
         ny = itynm(m)
         nz = itznm(m)

         if( itout(m) .le. 1 ) then
           nei = 0
         else
           nei = 1
         end if

         call sumtal_calc_stdev_depstxyz(m,ntf,
     &        np,nei,ne,nt,nx,ny,nz,
     &        sumfactor(m),weightRate(1,m),
     &        tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &        resc2SUMTAL(m),resc3SUMTAL(m),
     &        ierr)              ! integration

* tetra mesh !FURUTA20190121
        else if( itmsh(m) .eq. 4 ) then
! mesh=tet can be treated same as mesh=reg
         np = itpan(m)
         nr = itrgn(m)
         ne = itenm(m)
         nt = ittnm(m)

         if( itout(m) .le. 1 ) then
           nei = 0
         else
           nei = 1
         end if

         call sumtal_calc_stdev_depstreg(m,ntf,
     &        np,nei,ne,nt,nr,
     &        sumfactor(m),weightRate(1,m),
     &        tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &        resc2SUMTAL(m),resc3SUMTAL(m),
     &        ierr)              ! integration

        end if

*.......................................................................

       else if ( itout(m).eq.2 ) then   ! output=Deposit
* region mesh
        if( itmsh(m) .eq. 1 ) then

         np = itpan(m)
         nr = itrgn(m)
         ne = itenm(m)
         nt = ittnm(m)

         if( itout(m) .le. 1 ) then
           nei = 0
         else
           nei = 1
         end if

         call sumtal_calc_stdev2_depstreg(m,ntf,nfile(m),
     &        np,nei,ne,nt,nr,
     &        sumfactor(m),weightRate(1,m),
     &        tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &        resc2SUMTAL(m),resc3SUMTAL(m),
     &        ierr)              ! integration

* r-z mesh
        else if ( itmsh(m) .eq. 2 ) then

         np = itpan(m)
         ne = itenm(m)
         nt = ittnm(m)
         nr = itrnm(m)
         nz = itznm(m)

         if( itout(m) .le. 1 ) then
           nei = 0
         else
           nei = 1
         end if

         call sumtal_calc_stdev2_depstrz(m,ntf,nfile(m),
     &        np,nei,ne,nt,nr,nz,
     &        sumfactor(m),weightRate(1,m),
     &        tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &        resc2SUMTAL(m),resc3SUMTAL(m),
     &        ierr)              ! integration

* xyz mesh
        else if ( itmsh(m) .eq. 3 ) then

         np = itpan(m)
         ne = itenm(m)
         nt = ittnm(m)
         nx = itxnm(m)
         ny = itynm(m)
         nz = itznm(m)

         if( itout(m) .le. 1 ) then
           nei = 0
         else
           nei = 1
         end if

         call sumtal_calc_stdev2_depstxyz(m,ntf,nfile(m),
     &        np,nei,ne,nt,nx,ny,nz,
     &        sumfactor(m),weightRate(1,m),
     &        tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &        resc2SUMTAL(m),resc3SUMTAL(m),
     &        ierr)              ! integration

* tetra mesh !FURUTA20190121
!nais   else if( itmsh(m) .eq. 1 ) then
        else if( itmsh(m) .eq. 4 ) then
! mesh=tet can be treated same as mesh=reg
         np = itpan(m)
         nr = itrgn(m)
         ne = itenm(m)
         nt = ittnm(m)

         if( itout(m) .le. 1 ) then
           nei = 0
         else
           nei = 1
         end if

         call sumtal_calc_stdev2_depstreg(m,ntf,nfile(m),
     &        np,nei,ne,nt,nr,
     &        sumfactor(m),weightRate(1,m),
     &        tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &        resc2SUMTAL(m),resc3SUMTAL(m),
     &        ierr)              ! integration

        end if

       else                     ! illegal output parameter
        write(*,'(/" ***** Error: Sumtally, illegal output parameter"
     &            ," in [t-deposit] *****"/)')
        ierr = 1
        return                  ! error return
       end if

*-----------------------------------------------------------------------
*         [t-deposit2]
*-----------------------------------------------------------------------

      else if ( ital(m) .eq. 14 ) then

       np  = itpan(m)
       nt  = ittnm(m)
       ne1 = itenm(m)
       ne2 = itenm2(m)

       call sumtal_calc_stdev2_deposit2reg(m,ntf,nfile(m),
     &                             np,nt,ne1,ne2,
     &                             sumfactor(m),weightRate(1,m),
     &       tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &                             resc2SUMTAL(m),resc3SUMTAL(m),
     &                             ierr)

*-----------------------------------------------------------------------
*         [t-sed]
*-----------------------------------------------------------------------

      else if ( ital(m) .eq. 15 ) then

* region mesh
       if( itmsh(m) .eq. 1 ) then

        np = itpan(m)
        ne = itenm(m)
        nr = itrgn(m)

        call sumtal_calc_stdev_sedreg(m,ntf,
     &                             np,ne,nr,
     &                             sumfactor(m),weightRate(1,m),
     &       tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &                             resc2SUMTAL(m),resc3SUMTAL(m),
     &                             ierr)

* r-z mesh
       else if ( itmsh(m) .eq. 2 ) then

        np = itpan(m)
        ne = itenm(m)
        nr = itrnm(m)
        nz = itznm(m)

        call sumtal_calc_stdev_sedrz(m,ntf,
     &                             np,ne,nr,nz,
     &                             sumfactor(m),weightRate(1,m),
     &       tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &                             resc2SUMTAL(m),resc3SUMTAL(m),
     &                             ierr)

* xyz mesh
       else if ( itmsh(m) .eq. 3 ) then

        np = itpan(m)
        ne = itenm(m)
        nx = itxnm(m)
        ny = itynm(m)
        nz = itznm(m)

        call sumtal_calc_stdev_sedxyz(m,ntf,
     &                             np,ne,nx,ny,nz,
     &                             sumfactor(m),weightRate(1,m),
     &       tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &                             resc2SUMTAL(m),resc3SUMTAL(m),
     &                             ierr)

       end if

*-----------------------------------------------------------------------
*         [t-point]
*-----------------------------------------------------------------------

      else if ( ital(m) .eq. 17 ) then

        np = itpan(m)
        nr = itmsh(m)
        ne = itenm(m)
        nm = itmst(m)
        nt = ittnm(m)

        call sumtal_calc_stdev_point(m,ntf,
     &                             np,nr,ne,nm,nt,
     &                             sumfactor(m),weightRate(1,m),
     &       tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &                             resc2SUMTAL(m),resc3SUMTAL(m),
     &                             ierr)

*-----------------------------------------------------------------------
*         [t-wwg]
*-----------------------------------------------------------------------

      else if ( ital(m) .eq. 18 ) then

* region mesh
       if( itmsh(m) .eq. 1 ) then

        np = itpan(m)
        nr = itrgn(m)
        ne = itenm(m)
        nm = itmst(m)
        nt = ittnm(m)

        call sumtal_calc_stdev_wwgreg(m,ntf,
     &       np,ne,nt,nr,nm,
     &       sumfactor(m),weightRate(1,m),
     &       tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &       resc2SUMTAL(m),resc3SUMTAL(m),
     &       ierr)              ! integration

* xyz mesh
       else if ( itmsh(m) .eq. 3 ) then

        np = itpan(m)
        nx = itxnm(m)
        ny = itynm(m)
        nz = itznm(m)
        ne = itenm(m)
        nm = itmst(m)
        nt = ittnm(m)

        call sumtal_calc_stdev_wwgxyz(m,ntf,
     &       np,ne,nt,nx,ny,nz,nm,
     &       sumfactor(m),weightRate(1,m),
     &       tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &       resc2SUMTAL(m),resc3SUMTAL(m),
     &       ierr)              ! integration

* tetra mesh !FURUTA20240110
       else if( itmsh(m) .eq. 4 ) then
! mesh=tet can be treated same as mesh=reg
        np = itpan(m)
        nr = itrgn(m)
        ne = itenm(m)
        nm = itmst(m)
        nt = ittnm(m)

        call sumtal_calc_stdev_wwgreg(m,ntf,
     &       np,ne,nt,nr,nm,
     &       sumfactor(m),weightRate(1,m),
     &       tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &       resc2SUMTAL(m),resc3SUMTAL(m),
     &       ierr)              ! integration

       end if

*-----------------------------------------------------------------------
*         [t-volume]
*-----------------------------------------------------------------------

      else if ( ital(m) .eq. 21 ) then

* region mesh
       if( itmsh(m) .eq. 1 ) then

        nr = itrgn(m)

        call sumtal_calc_stdev_vlmreg(m,ntf,nr,
     &       sumfactor(m),weightRate(1,m),
     &       tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &       resc2SUMTAL(m),resc3SUMTAL(m),
     &       ierr)              ! integration


       end if

*-----------------------------------------------------------------------

      end if

!-----------------------------------------------------------------------
      return

      end subroutine sumtal_calc_stdev
!***********************************************************************


!***********************************************************************
!                                                                      *
      subroutine sumtal_calc_average(m,ntf,ierr)
!                                                                      *
!***********************************************************************
      use TALMOD,   only : tr0
      use RESTALMOD,   only : irestalm
      use sumtallymod
      use partmod, only: itpan, itpat, jtpat ! frtati 2021/10/05
      implicit none
!-----------------------------------------------------------------------
      include 'param.inc'
!-----------------------------------------------------------------------

      integer itnm
      integer ital
      integer itals
      integer italm
      common /tall00/ itnm, ital(itlmax), itals(itlmax), italm(itlmax)

      integer itmsh
      integer itunt
      integer itspc
      integer itout
      integer ittwo
      common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &     itout(itlmax), ittwo(itlmax)

      integer itrgn
      integer itrgm
      integer itreg
      common /tall02/ itrgn(itlmax), itrgm(itlmax), itreg(itlmax)

      integer itrty
      integer itrnm
      integer itrrg
      double precision rtrmi
      double precision rtrma
      double precision rtrdl
      common /tall03/ itrty(itlmax), itrnm(itlmax), itrrg(itlmax),
     &     rtrmi(itlmax), rtrma(itlmax), rtrdl(itlmax)

      integer itxty
      integer itxnm
      integer itxrg
      double precision rtxmi
      double precision rtxma
      double precision rtxdl
      common /tall04/ itxty(itlmax), itxnm(itlmax), itxrg(itlmax),
     &     rtxmi(itlmax), rtxma(itlmax), rtxdl(itlmax)

      integer ityty
      integer itynm
      integer ityrg
      double precision rtymi
      double precision rtyma
      double precision rtydl
      common /tall05/ ityty(itlmax), itynm(itlmax), ityrg(itlmax),
     &     rtymi(itlmax), rtyma(itlmax), rtydl(itlmax)

      integer itzty
      integer itznm
      integer itzrg
      double precision rtzmi
      double precision rtzma
      double precision rtzdl
      common /tall06/ itzty(itlmax), itznm(itlmax), itzrg(itlmax),
     &     rtzmi(itlmax), rtzma(itlmax), rtzdl(itlmax)

      integer itety
      integer itenm
      integer iterg
      double precision rtemi
      double precision rtema
      double precision rtedl
      common /tall07/ itety(itlmax), itenm(itlmax), iterg(itlmax),
     &     rtemi(itlmax), rtema(itlmax), rtedl(itlmax)


      integer itman
      integer itmat
      integer itmct
      integer itnun
      integer itnuc
      integer itndz
      integer itndn
      integer itnkz
      integer itnkn
      common /tall12/ itman(itlmax), itmat(itlmax), itmct(itlmax),
     &                itnun(itlmax), itnuc(itlmax),
     &                itndz(itlmax), itndn(itlmax),
     &                itnkz(itlmax), itnkn(itlmax)

      integer itrcn
      integer itrcr
      integer itrca
      common /tall13/ itrcn(itlmax), itrcr(itlmax), itrca(itlmax)

      integer itndy
      common /tall17/ itndy(itlmax)

      integer ittty
      integer ittnm
      integer ittrg
      double precision rttmi
      double precision rttma
      double precision rttdl
      common /tall20/ ittty(itlmax), ittnm(itlmax), ittrg(itlmax),
     &     rttmi(itlmax), rttma(itlmax), rttdl(itlmax)

      integer itaty
      integer itanm
      integer itarg
      double precision rtami
      double precision rtama
      double precision rtadl
      common /tall23/ itaty(itlmax), itanm(itlmax), itarg(itlmax),
     &                rtami(itlmax), rtama(itlmax), rtadl(itlmax)

      integer itmlp
      integer itmln
      integer itmst
      integer itmli
      double precision rtmme
      integer itmnt
      integer itmpn
      integer itmpt
      common /tall45/ itmlp(itlmax), itmln(itlmax,6), itmst(itlmax),
     &     itmli(itlmax,6), rtmme(itlmax,6), itmnt(itlmax,6),
     &     itmpn(itlmax,6), itmpt(itlmax,6,6,2)

      integer itety2
      integer itenm2
      integer iterg2
      double precision rtemi2
      double precision rtema2
      double precision rtedl2
      common /tall56/ itety2(itlmax), itenm2(itlmax), iterg2(itlmax),
     &                rtemi2(itlmax), rtema2(itlmax), rtedl2(itlmax)

      integer itactnm
      integer itactrg
      integer itactmax
      common /tall72/ itactnm(itlmax), itactrg(itlmax), itactmax(itlmax)

      integer itnzn
      integer itndm
      common /tall83/ itnzn(itlmax), itndm(itlmax)

!-----------------------------------------------------------------------

       integer, intent(out)  :: ierr           ! return code (=0:normal retrun)

       integer :: m

!-----------------------------------------------------------------------

       integer           ntf
       integer           np
       integer           ne, nei
       integer           nt
       integer           nx
       integer           ny
       integer           nz
       integer           nr
       integer           nact      ! S.Abe 2018/02/15
       integer           nm
       integer           na
       integer           nd
       integer           mn
       integer           mz
       integer           mm ! frtati 2022/02/18
       integer           ne1, ne2
       integer           idasz

!-----------------------------------------------------------------------

       ierr = 0

!-----------------------------------------------------------------------

*-----------------------------------------------------------------------
*         [t-track]
*-----------------------------------------------------------------------

      if ( ital(m) .eq. 1 ) then

* region mesh
       if( itmsh(m) .eq. 1 ) then

        np = itpan(m)
        nr = itrgn(m)
        ne = itenm(m)
        nm = itmst(m)
        nt = ittnm(m)

        call sumtal_calc_average_tracreg(m,ntf,nfile(m),
     &       np,ne,nt,nr,nm,
     &       sumfactor(m),weightRate(1,m),sumWR(m),
     &       tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &       resc2SUMTAL(m),resc3SUMTAL(m),
     &       ierr)              ! weighted average

* r-z mesh
       else if ( itmsh(m) .eq. 2 ) then

        np = itpan(m)
        nr = itrnm(m)
        nz = itznm(m)
        ne = itenm(m)
        nm = itmst(m)
        nt = ittnm(m)
        na = itanm(m)


        call sumtal_calc_average_tracrz(m,ntf,nfile(m),
     &       np,ne,nt,na,nr,nz,nm,
     &       sumfactor(m),weightRate(1,m),sumWR(m),
     &       tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &       resc2SUMTAL(m),resc3SUMTAL(m),
     &       ierr)              ! weighted average

* xyz mesh
       else if ( itmsh(m) .eq. 3 ) then

        np = itpan(m)
        nx = itxnm(m)
        ny = itynm(m)
        nz = itznm(m)
        ne = itenm(m)
        nm = itmst(m)
        nt = ittnm(m)

        call sumtal_calc_average_tracxyz(m,ntf,nfile(m),
     &       np,ne,nt,nx,ny,nz,nm,
     &       sumfactor(m),weightRate(1,m),sumWR(m),
     &       tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &       resc2SUMTAL(m),resc3SUMTAL(m),
     &       ierr)              ! weighted average

* tetra mesh !FURUTA20190110
       else if( itmsh(m) .eq. 4 ) then
! mesh=tet can be treated same as mesh=reg
        np = itpan(m)
        nr = itrgn(m)
        ne = itenm(m)
        nm = itmst(m)
        nt = ittnm(m)

        call sumtal_calc_average_tracreg(m,ntf,nfile(m),
     &       np,ne,nt,nr,nm,
     &       sumfactor(m),weightRate(1,m),sumWR(m),
     &       tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &       resc2SUMTAL(m),resc3SUMTAL(m),
     &       ierr)              ! weighted average

       end if

*-----------------------------------------------------------------------
*         [t-adjoint]
*-----------------------------------------------------------------------

      else if ( ital(m) .eq. 19 ) then

* region mesh
       if( itmsh(m) .eq. 1 ) then

        np = itpan(m)
        nr = itrgn(m)
        ne = itenm(m)
        nm = itmst(m)
        nt = ittnm(m)

        call sumtal_calc_average_tadjreg(m,ntf,nfile(m),
     &       np,ne,nt,nr,nm,
     &       sumfactor(m),weightRate(1,m),sumWR(m),
     &       tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &       resc2SUMTAL(m),resc3SUMTAL(m),
     &       ierr)              ! weighted average

* r-z mesh
       else if ( itmsh(m) .eq. 2 ) then

        np = itpan(m)
        nr = itrnm(m)
        nz = itznm(m)
        ne = itenm(m)
        nm = itmst(m)
        nt = ittnm(m)
        na = itanm(m)


        call sumtal_calc_average_tadjrz(m,ntf,nfile(m),
     &       np,ne,nt,na,nr,nz,nm,
     &       sumfactor(m),weightRate(1,m),sumWR(m),
     &       tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &       resc2SUMTAL(m),resc3SUMTAL(m),
     &       ierr)              ! weighted average

* xyz mesh
       else if ( itmsh(m) .eq. 3 ) then

        np = itpan(m)
        nx = itxnm(m)
        ny = itynm(m)
        nz = itznm(m)
        ne = itenm(m)
        nm = itmst(m)
        nt = ittnm(m)

        call sumtal_calc_average_tadjxyz(m,ntf,nfile(m),
     &       np,ne,nt,nx,ny,nz,nm,
     &       sumfactor(m),weightRate(1,m),sumWR(m),
     &       tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &       resc2SUMTAL(m),resc3SUMTAL(m),
     &       ierr)              ! weighted average

       end if

*-----------------------------------------------------------------------
*         [t-cross]
*-----------------------------------------------------------------------

      else if ( ital(m) .eq. 2 ) then

* region mesh
       if( itmsh(m) .eq. 1 ) then

        np = itpan(m)
        ne = itenm(m)
        nt = ittnm(m)
        na = itanm(m)
        nr = itrcn(m)
        nm = itmst(m)

        call sumtal_calc_average_crsreg(m,ntf,nfile(m),
     &                               np,ne,nt,na,nr,nm,
     &       sumfactor(m),weightRate(1,m),sumWR(m),
     &       tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &                               resc2SUMTAL(m),resc3SUMTAL(m),
     &                               ierr)

* r-z mesh
       else if ( itmsh(m) .eq. 2 ) then

        np = itpan(m)
        ne = itenm(m)
        na = itanm(m)
        nt = ittnm(m)
        nr = itrnm(m)
        nz = itznm(m)
        nm = itmst(m)

        idasz = itenm(m) * itpan(m) * 2
     &        * ( itrnm(m) + 1 ) * itznm(m)
     &        * itanm(m)
     &        * ittnm(m)

        call sumtal_calc_average_crsrz(m,ntf,nfile(m),
     &                               np,ne,na,nt,nr,nz,nm,
     &       sumfactor(m),weightRate(1,m),sumWR(m),
     &       tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &       tr0(irestalm(m)+idasz),trSUMTAL(irestalm(m)+idasz),
     &                               resc2SUMTAL(m),resc3SUMTAL(m),
     &                               ierr)

* xyz mesh
       else if ( itmsh(m) .eq. 3 ) then

        nt = ittnm(m)
        np = itpan(m)
        na = itanm(m)
        ne = itenm(m)
        ny = itynm(m)
        nx = itxnm(m)
        nz = itznm(m)
        nm = itmst(m)

        call sumtal_calc_average_crsxyz(m,ntf,nfile(m),
     &                               nt,np,na,ne,ny,nx,nz,nm,
     &       sumfactor(m),weightRate(1,m),sumWR(m),
     &       tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &                               resc2SUMTAL(m),resc3SUMTAL(m),
     &                               ierr)

       end if

*-----------------------------------------------------------------------
*         [t-yield]
*-----------------------------------------------------------------------

      else if ( ital(m) .eq. 3 ) then

* region mesh
       if( itmsh(m) .eq. 1 ) then

        mn = itndn(m)
        mz = itndz(m)
        mm = itndm(m)
        nr = itrgn(m)

        call sumtal_calc_average_yieldreg(m,ntf,nfile(m),
     &                               mn,mz,mm,nr,
     &       sumfactor(m),weightRate(1,m),sumWR(m),
     &       tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &                               resc2SUMTAL(m),resc3SUMTAL(m),
     &                               ierr)

* r-z mesh
       else if ( itmsh(m) .eq. 2 ) then

        mn = itndn(m)
        mz = itndz(m)
        mm = itndm(m)
        nz = itznm(m)
        nr = itrnm(m)

        call sumtal_calc_average_yieldrz(m,ntf,nfile(m),
     &                               mn,mz,mm,nz,nr,
     &       sumfactor(m),weightRate(1,m),sumWR(m),
     &       tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &                               resc2SUMTAL(m),resc3SUMTAL(m),
     &                               ierr)

* xyz mesh
       else if ( itmsh(m) .eq. 3 ) then

        mn = itndn(m)
        mz = itndz(m)
        mm = itndm(m)
        nz = itznm(m)
        ny = itynm(m)
        nx = itxnm(m)

        call sumtal_calc_average_yieldxyz(m,ntf,nfile(m),
     &                               mn,mz,mm,nz,ny,nx,
     &       sumfactor(m),weightRate(1,m),sumWR(m),
     &       tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &                               resc2SUMTAL(m),resc3SUMTAL(m),
     &                               ierr)

* tetra mesh !FURUTA20190121
       else if( itmsh(m) .eq. 4 ) then
! mesh=tet can be treated same as mesh=reg
        mn = itndn(m)
        mz = itndz(m)
        mm = itndm(m)
        nr = itrgn(m)

        call sumtal_calc_average_yieldreg(m,ntf,nfile(m),
     &                               mn,mz,mm,nr,
     &       sumfactor(m),weightRate(1,m),sumWR(m),
     &       tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &                               resc2SUMTAL(m),resc3SUMTAL(m),
     &                               ierr)

       end if

*-----------------------------------------------------------------------
*         [t-heat]
*-----------------------------------------------------------------------

      else if ( ital(m) .eq. 4 ) then

* region mesh
       if( itmsh(m) .eq. 1 ) then

        ne = itenm(m)
        nd = itndy(m)
        nr = itrgn(m)

        call sumtal_calc_average_hetreg(m,ntf,nfile(m),
     &                               ne,nd,nr,
     &       sumfactor(m),weightRate(1,m),sumWR(m),
     &       tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &                               resc2SUMTAL(m),resc3SUMTAL(m),
     &                               ierr)

* r-z mesh
       else if ( itmsh(m) .eq. 2 ) then

        ne = itenm(m)
        nd = itndy(m)
        nz = itznm(m)
        nr = itrnm(m)

        call sumtal_calc_average_hetrz(m,ntf,nfile(m),
     &                               ne,nd,nz,nr,
     &       sumfactor(m),weightRate(1,m),sumWR(m),
     &       tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &                               resc2SUMTAL(m),resc3SUMTAL(m),
     &                               ierr)

* xyz mesh
       else if ( itmsh(m) .eq. 3 ) then

        ne = itenm(m)
        nd = itndy(m)
        nz = itznm(m)
        ny = itynm(m)
        nx = itxnm(m)

        call sumtal_calc_average_hetxyz(m,ntf,nfile(m),
     &                               ne,nd,nz,ny,nx,
     &       sumfactor(m),weightRate(1,m),sumWR(m),
     &       tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &                               resc2SUMTAL(m),resc3SUMTAL(m),
     &                               ierr)

       end if

*-----------------------------------------------------------------------
*         [t-star]
*-----------------------------------------------------------------------

      else if ( ital(m) .eq. 5 ) then

* region mesh
       if( itmsh(m) .eq. 1 ) then

        np = itpan(m)
        ne = itenm(m)
        nt = ittnm(m)
        nr = itrgn(m)
        nact = itactnm(m)      ! S.Abe 2018/02/15

        call sumtal_calc_average_starreg(m,ntf,nfile(m),
     &                               np,ne,nt,nr,nact,      ! S.Abe 2018/02/15
     &       sumfactor(m),weightRate(1,m),sumWR(m),
     &       tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &                               resc2SUMTAL(m),resc3SUMTAL(m),
     &                               ierr)

* r-z mesh
       else if ( itmsh(m) .eq. 2 ) then

        np = itpan(m)
        ne = itenm(m)
        nt = ittnm(m)
        nr = itrnm(m)
        nz = itznm(m)

        call sumtal_calc_average_starrz(m,ntf,nfile(m),
     &                               np,ne,nt,nr,nz,
     &       sumfactor(m),weightRate(1,m),sumWR(m),
     &       tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &                               resc2SUMTAL(m),resc3SUMTAL(m),
     &                               ierr)

* xyz mesh
       else if ( itmsh(m) .eq. 3 ) then

        np = itpan(m)
        ne = itenm(m)
        nt = ittnm(m)
        nx = itxnm(m)
        ny = itynm(m)
        nz = itznm(m)

        call sumtal_calc_average_starxyz(m,ntf,nfile(m),
     &                               np,ne,nt,nx,ny,nz,
     &       sumfactor(m),weightRate(1,m),sumWR(m),
     &       tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &                               resc2SUMTAL(m),resc3SUMTAL(m),
     &                               ierr)

       end if

*-----------------------------------------------------------------------
*         [t-time]
*-----------------------------------------------------------------------

      else if ( ital(m) .eq. 6 ) then

* region mesh
       if( itmsh(m) .eq. 1 ) then

        np = itpan(m)
        nt = ittnm(m)
        ne = itenm(m)
        nr = itrgn(m)

        call sumtal_calc_average_timereg(m,ntf,nfile(m),
     &                               np,nt,ne,nr,
     &       sumfactor(m),weightRate(1,m),sumWR(m),
     &       tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &                               resc2SUMTAL(m),resc3SUMTAL(m),
     &                               ierr)

* r-z mesh
       else if ( itmsh(m) .eq. 2 ) then

        np = itpan(m)
        nt = ittnm(m)
        ne = itenm(m)
        nr = itrnm(m)
        nz = itznm(m)

        call sumtal_calc_average_timerz(m,ntf,nfile(m),
     &                               np,nt,ne,nr,nz,
     &       sumfactor(m),weightRate(1,m),sumWR(m),
     &       tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &                               resc2SUMTAL(m),resc3SUMTAL(m),
     &                               ierr)

* xyz mesh
       else if ( itmsh(m) .eq. 3 ) then

        np = itpan(m)
        nt = ittnm(m)
        ne = itenm(m)
        nx = itxnm(m)
        ny = itynm(m)
        nz = itznm(m)

        call sumtal_calc_average_timexyz(m,ntf,nfile(m),
     &                               np,nt,ne,nx,ny,nz,
     &       sumfactor(m),weightRate(1,m),sumWR(m),
     &       tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &                               resc2SUMTAL(m),resc3SUMTAL(m),
     &                               ierr)

       end if

*-----------------------------------------------------------------------
*         [t-dpa]
*-----------------------------------------------------------------------

      else if ( ital(m) .eq. 7 ) then

* region mesh
       if( itmsh(m) .eq. 1 ) then

        nd = itndy(m)
        np = itpan(m)
        nr = itrgn(m)

        call sumtal_calc_average_dpareg(m,ntf,nfile(m),
     &                               nd,np,nr,
     &       sumfactor(m),weightRate(1,m),sumWR(m),
     &       tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &                               resc2SUMTAL(m),resc3SUMTAL(m),
     &                               ierr)

* r-z mesh
       else if ( itmsh(m) .eq. 2 ) then

        nd = itndy(m)
        np = itpan(m)
        nr = itrnm(m)
        nz = itznm(m)

        call sumtal_calc_average_dparz(m,ntf,nfile(m),
     &                               nd,np,nr,nz,
     &       sumfactor(m),weightRate(1,m),sumWR(m),
     &       tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &                               resc2SUMTAL(m),resc3SUMTAL(m),
     &                               ierr)

* xyz mesh
       else if ( itmsh(m) .eq. 3 ) then

        nd = itndy(m)
        np = itpan(m)
        nx = itxnm(m)
        ny = itynm(m)
        nz = itznm(m)

        call sumtal_calc_average_dpaxyz(m,ntf,nfile(m),
     &                               nd,np,nx,ny,nz,
     &       sumfactor(m),weightRate(1,m),sumWR(m),
     &       tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &                               resc2SUMTAL(m),resc3SUMTAL(m),
     &                               ierr)

       end if

*-----------------------------------------------------------------------
*         [t-product]
*-----------------------------------------------------------------------

      else if ( ital(m) .eq. 8 ) then

* region mesh
       if( itmsh(m) .eq. 1 ) then

        np = itpan(m)
        na = itanm(m)
        nt = ittnm(m)
        ne = itenm(m)
        nr = itrgn(m)

        call sumtal_calc_average_productreg(m,ntf,nfile(m),
     &                               np,na,nt,ne,nr,
     &       sumfactor(m),weightRate(1,m),sumWR(m),
     &       tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &                               resc2SUMTAL(m),resc3SUMTAL(m),
     &                               ierr)

* r-z mesh
       else if ( itmsh(m) .eq. 2 ) then

        np = itpan(m)
        na = itanm(m)
        nt = ittnm(m)
        ne = itenm(m)
        nr = itrnm(m)
        nz = itznm(m)

        call sumtal_calc_average_productrz(m,ntf,nfile(m),
     &                               np,na,nt,ne,nr,nz,
     &       sumfactor(m),weightRate(1,m),sumWR(m),
     &       tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &                               resc2SUMTAL(m),resc3SUMTAL(m),
     &                               ierr)

* xyz mesh
       else if ( itmsh(m) .eq. 3 ) then

        np = itpan(m)
        na = itanm(m)
        nt = ittnm(m)
        ne = itenm(m)
        nx = itxnm(m)
        ny = itynm(m)
        nz = itznm(m)

        call sumtal_calc_average_productxyz(m,ntf,nfile(m),
     &                               np,na,nt,ne,nx,ny,nz,
     &       sumfactor(m),weightRate(1,m),sumWR(m),
     &       tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &                               resc2SUMTAL(m),resc3SUMTAL(m),
     &                               ierr)

       end if

*-----------------------------------------------------------------------
*         [t-let]
*-----------------------------------------------------------------------

      else if ( ital(m) .eq. 12 ) then

* region mesh
       if( itmsh(m) .eq. 1 ) then

        np = itpan(m)
        ne = itenm(m)
        nr = itrgn(m)

        call sumtal_calc_average_letreg(m,ntf,nfile(m),
     &                               np,ne,nr,
     &       sumfactor(m),weightRate(1,m),sumWR(m),
     &       tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &                               resc2SUMTAL(m),resc3SUMTAL(m),
     &                               ierr)

* r-z mesh
       else if ( itmsh(m) .eq. 2 ) then

        np = itpan(m)
        ne = itenm(m)
        nr = itrnm(m)
        nz = itznm(m)

        call sumtal_calc_average_letrz(m,ntf,nfile(m),
     &                               np,ne,nr,nz,
     &       sumfactor(m),weightRate(1,m),sumWR(m),
     &       tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &                               resc2SUMTAL(m),resc3SUMTAL(m),
     &                               ierr)

* xyz mesh
       else if ( itmsh(m) .eq. 3 ) then

        np = itpan(m)
        ne = itenm(m)
        nx = itxnm(m)
        ny = itynm(m)
        nz = itznm(m)

        call sumtal_calc_average_letxyz(m,ntf,nfile(m),
     &                               np,ne,nx,ny,nz,
     &       sumfactor(m),weightRate(1,m),sumWR(m),
     &       tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &                               resc2SUMTAL(m),resc3SUMTAL(m),
     &                               ierr)

       end if

*-----------------------------------------------------------------------
*         [t-deposit]
*-----------------------------------------------------------------------

      else if ( ital(m) .eq. 13 ) then

       if ( itout(m).eq.1 ) then   ! output=Dose

* region mesh
        if( itmsh(m) .eq. 1 ) then

         np = itpan(m)
         nr = itrgn(m)
         ne = itenm(m)
         nt = ittnm(m)

         if( itout(m) .le. 1 ) then
           nei = 0
         else
           nei = 1
         end if

         call sumtal_calc_average_depstreg(m,ntf,nfile(m),
     &        np,nei,ne,nt,nr,
     &        sumfactor(m),weightRate(1,m),sumWR(m),
     &        tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &        resc2SUMTAL(m),resc3SUMTAL(m),
     &        ierr)              ! integration

* r-z mesh
        else if ( itmsh(m) .eq. 2 ) then

         np = itpan(m)
         nr = itrnm(m)
         nz = itznm(m)
         ne = itenm(m)
         nt = ittnm(m)

         if( itout(m) .le. 1 ) then
           nei = 0
         else
           nei = 1
         end if

         call sumtal_calc_average_depstrz(m,ntf,nfile(m),
     &        np,nei,ne,nt,nr,nz,
     &        sumfactor(m),weightRate(1,m),sumWR(m),
     &        tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &        resc2SUMTAL(m),resc3SUMTAL(m),
     &        ierr)              ! integration

* xyz mesh
        else if ( itmsh(m) .eq. 3 ) then

         np = itpan(m)
         nx = itxnm(m)
         ny = itynm(m)
         nz = itznm(m)
         ne = itenm(m)
         nt = ittnm(m)

         if( itout(m) .le. 1 ) then
           nei = 0
         else
           nei = 1
         end if

         call sumtal_calc_average_depstxyz(m,ntf,nfile(m),
     &        np,nei,ne,nt,nx,ny,nz,
     &        sumfactor(m),weightRate(1,m),sumWR(m),
     &        tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &        resc2SUMTAL(m),resc3SUMTAL(m),
     &        ierr)              ! weighted average

* tetra mesh !FURUTA20190121
        else if( itmsh(m) .eq. 4 ) then
! mesh=tet can be treated same as mesh=reg
         np = itpan(m)
         nr = itrgn(m)
         ne = itenm(m)
         nt = ittnm(m)

         if( itout(m) .le. 1 ) then
           nei = 0
         else
           nei = 1
         end if

         call sumtal_calc_average_depstreg(m,ntf,nfile(m),
     &        np,nei,ne,nt,nr,
     &        sumfactor(m),weightRate(1,m),sumWR(m),
     &        tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &        resc2SUMTAL(m),resc3SUMTAL(m),
     &        ierr)              ! integration

        end if

*.......................................................................

       else if ( itout(m).eq.2 ) then   ! output=Deposit
* region mesh
        if( itmsh(m) .eq. 1 ) then

         np = itpan(m)
         nr = itrgn(m)
         ne = itenm(m)
         nt = ittnm(m)

         if( itout(m) .le. 1 ) then
           nei = 0
         else
           nei = 1
         end if

         call sumtal_calc_average2_depstreg(m,ntf,nfile(m),
     &        np,nei,ne,nt,nr,
     &        sumfactor(m),weightRate(1,m),sumWR(m),
     &        tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &        resc2SUMTAL(m),resc3SUMTAL(m),
     &        ierr)              ! integration

* r-z mesh
        else if ( itmsh(m) .eq. 2 ) then

         np = itpan(m)
         nr = itrnm(m)
         nz = itznm(m)
         ne = itenm(m)
         nt = ittnm(m)

         if( itout(m) .le. 1 ) then
           nei = 0
         else
           nei = 1
         end if

         call sumtal_calc_average2_depstrz(m,ntf,nfile(m),
     &        np,nei,ne,nt,nr,nz,
     &        sumfactor(m),weightRate(1,m),sumWR(m),
     &        tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &        resc2SUMTAL(m),resc3SUMTAL(m),
     &        ierr)              ! integration

* xyz mesh
        else if ( itmsh(m) .eq. 3 ) then

         np = itpan(m)
         nx = itxnm(m)
         ny = itynm(m)
         nz = itznm(m)
         ne = itenm(m)
         nt = ittnm(m)

         if( itout(m) .le. 1 ) then
           nei = 0
         else
           nei = 1
         end if

         call sumtal_calc_average2_depstxyz(m,ntf,nfile(m),
     &        np,nei,ne,nt,nx,ny,nz,
     &        sumfactor(m),weightRate(1,m),sumWR(m),
     &        tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &        resc2SUMTAL(m),resc3SUMTAL(m),
     &        ierr)              ! weighted average

* tetra mesh !FURUTA20190121
        else if( itmsh(m) .eq. 4 ) then
! mesh=tet can be treated same as mesh=reg
         np = itpan(m)
         nr = itrgn(m)
         ne = itenm(m)
         nt = ittnm(m)

         if( itout(m) .le. 1 ) then
           nei = 0
         else
           nei = 1
         end if

         call sumtal_calc_average2_depstreg(m,ntf,nfile(m),
     &        np,nei,ne,nt,nr,
     &        sumfactor(m),weightRate(1,m),sumWR(m),
     &        tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &        resc2SUMTAL(m),resc3SUMTAL(m),
     &        ierr)              ! integration

        end if

       else                     ! illegal output parameter
        write(*,'(/" ***** Error: Sumtally, illegal output parameter"
     &            ," in [t-deposit] *****"/)')
        ierr = 1
        return                  ! error return
       end if

*-----------------------------------------------------------------------
*         [t-deposit2]
*-----------------------------------------------------------------------

      else if ( ital(m) .eq. 14 ) then

        np  = itpan(m)
        nt  = ittnm(m)
        ne1 = itenm(m)
        ne2 = itenm2(m)

       call sumtal_calc_averag2_deposit2reg(m,ntf,nfile(m),
     &                               np,nt,ne1,ne2,
     &       sumfactor(m),weightRate(1,m),sumWR(m),
     &       tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &                               resc2SUMTAL(m),resc3SUMTAL(m),
     &                               ierr)

*-----------------------------------------------------------------------
*         [t-sed]
*-----------------------------------------------------------------------

      else if ( ital(m) .eq. 15 ) then

* region mesh
       if( itmsh(m) .eq. 1 ) then

        np = itpan(m)
        ne = itenm(m)
        nr = itrgn(m)

        call sumtal_calc_average_sedreg(m,ntf,nfile(m),
     &                               np,ne,nr,
     &       sumfactor(m),weightRate(1,m),sumWR(m),
     &       tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &                               resc2SUMTAL(m),resc3SUMTAL(m),
     &                               ierr)

* r-z mesh
       else if ( itmsh(m) .eq. 2 ) then

        np = itpan(m)
        ne = itenm(m)
        nr = itrnm(m)
        nz = itznm(m)

        call sumtal_calc_average_sedrz(m,ntf,nfile(m),
     &                               np,ne,nr,nz,
     &       sumfactor(m),weightRate(1,m),sumWR(m),
     &       tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &                               resc2SUMTAL(m),resc3SUMTAL(m),
     &                               ierr)

* xyz mesh
       else if ( itmsh(m) .eq. 3 ) then

        np = itpan(m)
        ne = itenm(m)
        nx = itxnm(m)
        ny = itynm(m)
        nz = itznm(m)

        call sumtal_calc_average_sedxyz(m,ntf,nfile(m),
     &                               np,ne,nx,ny,nz,
     &       sumfactor(m),weightRate(1,m),sumWR(m),
     &       tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &                               resc2SUMTAL(m),resc3SUMTAL(m),
     &                               ierr)

       end if

*-----------------------------------------------------------------------
*         [t-point]
*-----------------------------------------------------------------------

      else if ( ital(m) .eq. 17 ) then

        np = itpan(m)
        nr = itmsh(m)
        ne = itenm(m)
        nm = itmst(m)
        nt = ittnm(m)

        call sumtal_calc_average_point(m,ntf,nfile(m),
     &                               np,nr,ne,nm,nt,
     &       sumfactor(m),weightRate(1,m),sumWR(m),
     &       tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &                               resc2SUMTAL(m),resc3SUMTAL(m),
     &                               ierr)

*-----------------------------------------------------------------------
*         [t-wwg]
*-----------------------------------------------------------------------

      else if ( ital(m) .eq. 18 ) then

* region mesh
       if( itmsh(m) .eq. 1 ) then

        np = itpan(m)
        nr = itrgn(m)
        ne = itenm(m)
        nm = itmst(m)
        nt = ittnm(m)

        call sumtal_calc_average_wwgreg(m,ntf,nfile(m),
     &       np,ne,nt,nr,nm,
     &       sumfactor(m),weightRate(1,m),sumWR(m),
     &       tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &       resc2SUMTAL(m),resc3SUMTAL(m),
     &       ierr)              ! weighted average

* xyz mesh
       else if ( itmsh(m) .eq. 3 ) then

        np = itpan(m)
        nx = itxnm(m)
        ny = itynm(m)
        nz = itznm(m)
        ne = itenm(m)
        nm = itmst(m)
        nt = ittnm(m)

        call sumtal_calc_average_wwgxyz(m,ntf,nfile(m),
     &       np,ne,nt,nx,ny,nz,nm,
     &       sumfactor(m),weightRate(1,m),sumWR(m),
     &       tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &       resc2SUMTAL(m),resc3SUMTAL(m),
     &       ierr)              ! weighted average

* tetra mesh !FURUTA20240110
       else if( itmsh(m) .eq. 4 ) then
! mesh=tet can be treated same as mesh=reg
        np = itpan(m)
        nr = itrgn(m)
        ne = itenm(m)
        nm = itmst(m)
        nt = ittnm(m)

        call sumtal_calc_average_wwgreg(m,ntf,nfile(m),
     &       np,ne,nt,nr,nm,
     &       sumfactor(m),weightRate(1,m),sumWR(m),
     &       tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &       resc2SUMTAL(m),resc3SUMTAL(m),
     &       ierr)              ! weighted average

       end if

*-----------------------------------------------------------------------
*         [t-volume]
*-----------------------------------------------------------------------

      else if ( ital(m) .eq. 21 ) then

* region mesh
       if( itmsh(m) .eq. 1 ) then

        nr = itrgn(m)

        call sumtal_calc_average_vlmreg(m,ntf,nfile(m),nr,
     &       sumfactor(m),weightRate(1,m),sumWR(m),
     &       tr0(irestalm(m)),trSUMTAL(irestalm(m)),
     &       resc2SUMTAL(m),resc3SUMTAL(m),
     &       ierr)              ! weighted average

       end if

*-----------------------------------------------------------------------

      end if

!-----------------------------------------------------------------------
      return

      end subroutine sumtal_calc_average
!***********************************************************************


!***********************************************************************
!                                                                      *
! sumover subroutine                                                   *
!                                                                      *
!***********************************************************************

!***********************************************************************
!                                                                      *
      subroutine calc_stdev_sum(kind, m, calcfact)
!                                                                      *
!***********************************************************************

      use TALMOD
!$      use TALMOD0

      implicit none

C for nonshared_tally option
      integer italsh
      common /talsh/ italsh

      integer :: m, iax, kind

      double precision :: calcfact

      real(8),pointer :: p_sum(:)

      integer :: ln_sum

      do iax =1,6

C for nonshared_tally option
!$       if(italsh .eq. 0) then
!$          call GET_TR_HEAD_POINTER0_SUM(p_sum,m,iax)
!$          ln_sum = italsize0_sum(m,iax)/2
!$       else
            call GET_TR_HEAD_POINTER_SUM(p_sum,m,iax)
            ln_sum = italsize_sum(m,iax)/2
C for nonshared_tally option
!$       end if

        if(ln_sum > 0 ) then

          call calc_stdev_sum_sub(m,p_sum,ln_sum,calcfact)
        endif
      enddo

      return
      end

!***********************************************************************
!                                                                      *
      subroutine calc_stdev_tr_sum(m, calcfact)
!                                                                      *
!***********************************************************************

      use TALMOD
!$      use TALMOD0

      implicit none

C for nonshared_tally option
      integer italsh
      common /talsh/ italsh

      integer :: kind, m, iax

      double precision :: calcfact

      real(8),pointer :: p_sum(:)

      integer :: ln_sum

      do iax =1,6

C for nonshared_tally option
!$       if(italsh .eq. 0) then
!$          call GET_TR_HEAD_POINTER0_SUM(p_sum,m,iax)
!$          ln_sum = italsize0_sum(m,iax)/2
!$       else
            call GET_TR_HEAD_POINTER_SUM(p_sum,m,iax)
            ln_sum = italsize_sum(m,iax)/2
C for nonshared_tally option
!$       end if

        if(ln_sum > 0 ) then

          call calc_stdev_sum_sub(m,p_sum,ln_sum,calcfact)
        endif
      enddo

      return
      end

!***********************************************************************
!                                                                      *
      subroutine calc_stdev_trtz_sum(kind, m, calcfact)
!                                                                      *
!***********************************************************************

      use TALMOD
!$      use TALMOD0

      implicit none

C for nonshared_tally option
      integer italsh
      common /talsh/ italsh

      integer :: kind, m, iax

      double precision :: calcfact

      real(8),pointer :: p_sum(:)

      integer :: ln_sum

      do iax =1,6

        if(kind == 0) then
C for nonshared_tally option
!$       if(italsh .eq. 0) then
!$          call GET_TR_HEAD_POINTER0_SUM(p_sum,m,iax)
!$          ln_sum = italsize0_2_sum(m,iax)/2
!$       else
            call GET_TR_HEAD_POINTER_SUM(p_sum,m,iax)
            ln_sum = italsize_2_sum(m,iax)/2
C for nonshared_tally option
!$       end if
        else
C for nonshared_tally option
!$       if(italsh .eq. 0) then
!$          call GET_TZ_HEAD_POINTER0_SUM(p_sum,m,iax)
!$          ln_sum=(italsize0_sum(m,iax) - italsize0_2_sum(m,iax))/2
!$       else
            call GET_TZ_HEAD_POINTER_SUM(p_sum,m,iax)
            ln_sum=(italsize_sum(m,iax)  - italsize_2_sum(m,iax))/2
C for nonshared_tally option
!$       end if
        endif

        if(ln_sum > 0 ) then

          call calc_stdev_sum_sub(m,p_sum,ln_sum,calcfact)
        endif
      enddo

      return
      end

!***********************************************************************
!                                                                      *
      subroutine calc_stdev_sum_sub(m,tr0_sum,ln_sum,calcfact)
!                                                                      *
!***********************************************************************

      implicit none

      integer :: m,ln_sum,i
      real(8) :: tr0_sum(ln_sum,2),calcfact,Xa,sigx

      do i=1,ln_sum
         if( tr0_sum(i,1) .gt. 0.d0 ) then
            call calc_stdev(m,Xa,sigx,
     &             tr0_sum(i,1), tr0_sum(i,2),
     &             calcfact)
            tr0_sum(i,1) = Xa
            tr0_sum(i,2) = sigx
         else
            tr0_sum(i,2) = 0.0d0
         end if
      enddo

      return
      end

!***********************************************************************
!                                                                      *
      subroutine calc_deposit_stdev_tr_sum(m, calcfact, np)
!                                                                      *
!***********************************************************************

      use TALMOD
!$      use TALMOD0

      implicit none

C for nonshared_tally option
      integer italsh
      common /talsh/ italsh

      integer :: kind, m, iax, np

      double precision :: calcfact

      real(8),pointer :: p_sum(:)

      integer :: ln_sum

      do iax =1,6

C for nonshared_tally option
!$       if(italsh .eq. 0) then
!$          call GET_TR_HEAD_POINTER0_SUM(p_sum,m,iax)
!$          ln_sum = italsize0_sum(m,iax)/2
!$       else
            call GET_TR_HEAD_POINTER_SUM(p_sum,m,iax)
            ln_sum = italsize_sum(m,iax)/2
C for nonshared_tally option
!$       end if

        if(ln_sum > 0 ) then

          call calc_deposit_stdev_sum_sub(m,p_sum,ln_sum,calcfact,
     &                                   np)
        endif
      enddo

      return
      end

!***********************************************************************
!                                                                      *
      subroutine calc_deposit_stdev_sum_sub(m,tr0_sum,ln_sum,calcfact,
     &                                      np)
!                                                                      *
!***********************************************************************

      implicit none

      integer :: m,ln_sum,i,np,ip
      real(8) :: tr0_sum(ln_sum,2),calcfact,Xa,sigx

      ip = 0
      do i=1,ln_sum
         ip = ip + 1
         if(ip > np) then
           ip = 1
         endif
         if( tr0_sum(i,1) .gt. 0.d0 ) then
            call calc_deposit_stdev(m,Xa,sigx,
     &             tr0_sum(i,1), tr0_sum(i,2),
     &             calcfact,ip)
            tr0_sum(i,1) = Xa
            tr0_sum(i,2) = sigx
         else
            tr0_sum(i,2) = 0.0d0
         end if
      enddo

      return
      end

!***********************************************************************
!                                                                      *
      subroutine sumtal_calc_average_sub(m, ntf,
     &                                      calcfact)
!                                                                      *
!***********************************************************************

      use RESTALMOD, only: irestalm_sum,lrestalm_sum
      use TALMOD
!$      use TALMOD0
      use sumtallymod

      implicit none

C for nonshared_tally option
      integer italsh
      common /talsh/ italsh

      integer :: m, iax, ntf

      double precision :: calcfact

      real(8),pointer :: p_sum(:)
!      real(8),pointer :: resta_sum(:)
      real(8),pointer :: sumtal_sum(:)

      integer :: ln_sum

      do iax =1,6
!       resta_sum => trRES_sum(irestalm_sum(m,iax):)

C for nonshared_tally option
!$       if(italsh .eq. 0) then
!$          call GET_TR_HEAD_POINTER0_SUM(p_sum,m,iax)
!$       else
            call GET_TR_HEAD_POINTER_SUM(p_sum,m,iax)
C for nonshared_tally option
!$       end if

        sumtal_sum => trSUMTAL_SUM(irestalm_sum(m,iax):)
        ln_sum = lrestalm_sum(m,iax)/2

        if(ln_sum > 0 ) then

          call calc_average_restal2sumtally(p_sum,sumtal_sum,
     &       ln_sum, ntf, calcfact)
        endif
      enddo

      return
      end

!***********************************************************************
!                                                                      *
      subroutine calc_average_restal2sumtally(trRES_SUM,trSUMTAL_SUM,
     &       ln_sum, ntf, calcfact)
!                                                                      *
!***********************************************************************

      implicit none

      integer :: ln_sum,ntf
      real(8) :: trRES_sum(ln_sum,2),trSUMTAL_SUM(ln_sum,2),calcfact

      if( ntf > 1 ) then

            ! X_bar = F Sigma rj/r Xj_bar
          trSUMTAL_SUM(1:ln_sum,1) =  trSUMTAL_SUM(1:ln_sum,1)
     &                     +  calcfact * trRES_SUM(1:ln_sum,1)


            ! sig_x = F sqrt{ Sigma (rj/r)**2 (sig_xj)**2 }
          trSUMTAL_SUM(1:ln_sum,2) =  trSUMTAL_SUM(1:ln_sum,2)
     &            + ( trRES_SUM(1:ln_sum,2)
     &            *   trRES_SUM(1:ln_sum,1) )**2
     &                *  calcfact**2

           ! ntf = 1
       else
            ! X_bar = F Sigma rj/r Xj_bar
            trSUMTAL_SUM(1:ln_sum,1) =
     &               calcfact * trRES_SUM(1:ln_sum,1)

            ! sig_x = F sqrt{ Sigma (rj/r)**2 (sig_xj)**2 }
            trSUMTAL_SUM(1:ln_sum,2) =
     &           ( trRES_SUM(1:ln_sum,2)
     &             * trRES_SUM(1:ln_sum,1) )**2
     &            *  calcfact**2
       end if

      return
      end

!***********************************************************************
!                                                                      *
      subroutine sumtal_calc_average_sub2(m, ntf,
     &           sumfactor1,sumWR1,resc2SUMTAL1,resc3SUMTAL1)
!                                                                      *
!***********************************************************************

      use RESTALMOD, only: irestalm_sum,lrestalm_sum
      use sumtallymod

      implicit none

C for nonshared_tally option
      integer italsh
      common /talsh/ italsh

      integer :: m, iax, ntf

      real(8) :: sumfactor1,sumWR1,resc2SUMTAL1,resc3SUMTAL1

      real(8),pointer :: sumtal_sum(:)

      integer :: ln_sum

      do iax =1,6

        sumtal_sum => trSUMTAL_SUM(irestalm_sum(m,iax):)
        ln_sum = lrestalm_sum(m,iax)/2

        if(ln_sum > 0 ) then

          call calc_average_sumtally(sumtal_sum,
     &       ln_sum, ntf,
     &       sumfactor1,sumWR1,resc2SUMTAL1,resc3SUMTAL1)
        endif
      enddo

      return
      end

!***********************************************************************
!                                                                      *
      subroutine calc_average_sumtally(trSUMTAL_SUM,ln_sum, ntf,
     &                  sumfactor,sumWR,resc2SUMTAL,resc3SUMTAL)
!                                                                      *
!***********************************************************************

      implicit none

      integer :: ln_sum,ntf,i
      real(8) :: trSUMTAL_SUM(ln_sum,2)
      real(8) :: sumfactor,sumWR,resc2SUMTAL,resc3SUMTAL

      do i=1,ln_sum
            ! X_bar
            trSUMTAL_sum(i,1) =
     &             trSUMTAL_sum(i,1) * sumfactor/sumWR

            ! sig_x
            trSUMTAL_sum(i,2) = sqrt(
     &           trSUMTAL_sum(i,2)) * sumfactor/sumWR


            ! Sigma(xi wi)**2 = (sig_X**2 N(N-1)+N X_bar**2)(W/N)**2
            trSUMTAL_sum(i,2) =
     &           (trSUMTAL_sum(i,2)**2 *
     &           resc3SUMTAL * (resc3SUMTAL-1.0d0) +
     &           resc3SUMTAL *
     &           trSUMTAL_sum(i,1)**2) *
     &           (resc2SUMTAL/resc3SUMTAL)**2

            ! Sigma xi wi = X_bar W
            trSUMTAL_sum(i,1) =
     &           trSUMTAL_sum(i,1) * resc2SUMTAL

            if( trSUMTAL_sum(i,1)  < 0d0 ) then
                trSUMTAL_sum(i,1)= 0d0
                trSUMTAL_sum(i,2)= 0d0
            end if
      enddo
      return
      end


!***********************************************************************
!                                                                      *
      subroutine sumtal_calc_average2_sub2(m, ntf, sumfactor1)
!                                                                      *
!***********************************************************************

      use RESTALMOD, only: irestalm_sum,lrestalm_sum
      use sumtallymod

      implicit none

C for nonshared_tally option
      integer italsh
      common /talsh/ italsh

      integer :: m, iax, ntf

      real(8) :: sumfactor1

      real(8),pointer :: sumtal_sum(:)

      integer :: ln_sum

      do iax =1,6

        sumtal_sum => trSUMTAL_SUM(irestalm_sum(m,iax):)
        ln_sum = lrestalm_sum(m,iax)/2

        if(ln_sum > 0 ) then

          call calc_average2_sumtally(sumtal_sum,
     &       ln_sum, ntf, sumfactor1)
        endif
      enddo

      return
      end

!***********************************************************************
!                                                                      *
      subroutine calc_average2_sumtally(trSUMTAL_SUM,ln_sum, ntf,
     &                  sumfactor)
!                                                                      *
!***********************************************************************

      implicit none

      integer :: ln_sum,ntf,i
      real(8) :: trSUMTAL_SUM(ln_sum,2)
      real(8) :: sumfactor,sumWR,resc2SUMTAL,resc3SUMTAL

      do i=1,ln_sum
             ! X_bar_k = F Sigma (rj/r) Xk_bar_j
            trSUMTAL_sum(i,1) =  sumfactor * trSUMTAL_sum(i,1)

C S.H. added the following processes for negative values when isumtally=3.
C This is tentative. Correct calculation can be done if processes of
C subroutine calc_stdev is performed here, but preparation of the last
C parameter of calc_stdev is needed. (2020.12.16)

            if( trSUMTAL_sum(i,1)  < 0d0 ) then
                trSUMTAL_sum(i,1)= 0d0
                trSUMTAL_sum(i,2)= 0d0
            end if


             if( trSUMTAL_sum(i,1) > 0.0d0 ) then
              ! F sqrt{ Sigma (rj/r)**2 (Rk_bar_j Xk_bar_j)**2} / X_bar_k
              trSUMTAL_sum(i,2) = sumfactor * sqrt(trSUMTAL_sum(i,2))
     &                          / trSUMTAL_sum(i,1)
             end if

      enddo
      return
      end
