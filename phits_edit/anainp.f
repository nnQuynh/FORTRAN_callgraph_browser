************************************************************************
*                                                                      *
      subroutine anainp(jsi,jsn,ivers,ierr)
*                                                                      *
*       make temporary input files when icntl=16                       *
*       modified by S.Hashimoto on 2020/02/20                          *
*                                                                      *
*       output :                                                       *
*           jsi      unit of reading file                              *
*         ivers      version of input file                             *
*          ierr      error flag                                        *
*                                                                      *
************************************************************************
      use CHARVARMOD, only: ErrLine_Adjust
      implicit real*8 (a-h,o-z)
      include 'err.inc'
*-----------------------------------------------------------------------

      common /error/  m_err, l_err, k_err
      character       m_err*200

      common /paraj/  mstz(300), parz(300)
      common /tcntl/  icntl, inucr

*-----------------------------------------------------------------------

      character chin*200, chlw*200, chcm*200
      character chlc*200

      character dsin(0:9)*200
      dimension idsi(0:9)

      dimension ill(0:9), ilf(0:9)

      logical   exex

      dimension ild(0:9)

*-----------------------------------------------------------------------

      character cunder*5,cnumber*20
      character ctallyfile*200
      integer ltallyfile

*-----------------------------------------------------------------------
      integer itnmw  ! tally number (with anatally start)
*-----------------------------------------------------------------------

      ioptmp = 799
      ioanat = 798
      ioanattmp = 797
      open(ioptmp, file='phits_tmp.inp', status='unknown')
      open(ioanat, file='anatally.inp', status='unknown')
      open(ioanattmp, file='anatally_tmp.inp', status='unknown')
      icntl = mstz(1) ! set icntl value

      itnmw = 0
*-----------------------------------------------------------------------

      if( ivers .eq. 0 ) then
         ierr = 1
         m_err = 'icntl=16 is not available in input of old version'
         ErrCha = ''
         ErrID = 'L:64/R:anainp/F:anainp.f'
         return
      end if

      dsin(0) = 'Error Line'
      idsi(0) = 12
      m_err = ' Error !!'
      ErrCha = ''
      ErrID = 'L:72/R:anainp/F:anainp.f'
      l_err = 1
      k_err = 0

*-----------------------------------------------------------------------

      write(*,*) 'icntl = 16: Making temporary files for anatally mode'

*-----------------------------------------------------------------------
*     read first non comment line of unit 5
*-----------------------------------------------------------------------

 514  continue

      read(jsi,'(a200)', iostat = ios ) chin
      if( ios .eq. -1 ) goto 998

      call chlngt(chin,200,i1,i2)
      chlw = chin
      call chcaps(chlw,i1,i2,i3,'#!$')
      chcm = chlw
      call chcomp(chcm,i1,i3,i4)

*-----------------------------------------------------------------------

      if( i1 .eq. 0 .and. i2 .eq. 0 ) then
         iskip = 1
      else if( index('#!$',chin(i1:i1)) .ne. 0 ) then
         iskip = 2
      else
         iskip = 0
      end if

      if( iskip .ne. 0 ) goto 514

*-----------------------------------------------------------------------
*     check unit 5 file
*     if first line is started as 'file =' -> this is input file name
*                                  other   -> read unit jsi = 5
*-----------------------------------------------------------------------

      if( chcm(i1:i1+4) .ne. 'file=' ) then
         rewind jsi
         goto 100
      end if

*-----------------------------------------------------------------------
*     get input file name from file =
*-----------------------------------------------------------------------

      ic1 = inumc(chlw,i1,i3,'=') + 1
      ic1 = jnumc(chlw,ic1,i3)
      ic2 = min( i3, inumc(chlw,ic1+1,i3,' ') )

      iname = ic2 - ic1 + 1

      do i = 1, iname
         dsin(jsn+1)(i:i) = chin(ic1+i-1:ic1+i-1)
      end do

      do i = iname + 1, 200
         dsin(jsn+1)(i:i) = ' '
      end do

      idsi(jsn+1) = iname

*-----------------------------------------------------------------------
*     input file dsin exist ?
*     open first input file
*-----------------------------------------------------------------------

      inquire( file = dsin(jsn+1), exist = exex )

      if( exex .eqv. .false. ) then

         m_err = 'Input File Name Error. File not exist'//
     &        ' ->> '//dsin(jsn+1)(1:iname)
         ErrCha = ''
         ErrID = 'L:150/R:anainp/F:anainp.f'
         l_err = 1
         k_err = jsn+1

         dsin(k_err) = 'Error Line'
         idsi(k_err) = 12

         goto 999

      end if

      call openf(jsi,jsn,dsin)

*-----------------------------------------------------------------------
*     check input file : whether [  ] is exist or not
*-----------------------------------------------------------------------

 100  continue

      ivers = 0

      read(jsi, '(a10)', iostat = ios ) chin
      if( ios .eq. -1 ) goto 200

      call chlngt(chin,10,i1,i2)

      if( chin(i1:i1) .ne. '[' ) goto 100

      ivers = 1

 200  continue

      rewind jsi

*-----------------------------------------------------------------------
*     initial values for include file
*-----------------------------------------------------------------------

      do i = 0, 9
         ill(i) = 0
         ilf(i) = 10000000
         ild(i) = 0
      end do
      jpn  = 0

*-----------------------------------------------------------------------
*     get information on parallel calculation
*-----------------------------------------------------------------------

      numdmpi = -1
      numdomp = -1
      do while ( jpn .ne. 3 ) ! search $MPI or $OMP
       call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!',
     &        jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

       if (  jpn.ne.3 .and. iskip.eq.0 ) then
        if( i1 .le. 5 .and. chlw(i1:i1) .eq. '[' ) then
         exit
        else
         if( chcm(i1:i1+3) .eq. '$mpi' ) then
          il = i1 + 3
          ivl = i4
          do k = i1+4, i4
           if ( chcm(k:k) .eq. '=' ) then
              il = k + 1
              exit
           end if
          end do
          call onum(chcm,il,i4,prn,ierr)
          numdmpi = int(prn)
         else if( chcm(i1:i1+3) .eq. '$omp' ) then
          il = i1 + 3
          ivl = i4
          do k = i1+4, i4
           if ( chcm(k:k) .eq. '=' ) then
              il = k + 1
              exit
           end if
          end do
          call onum(chcm,il,i4,prn,ierr)
          numdomp = int(prn)
         end if
        end if
       end if
      end do

      if ( numdmpi .ge. 0 ) write(ioptmp,'(a7,i7)') ' $MPI = ',numdmpi
      if ( numdomp .ge. 0 ) write(ioptmp,'(a7,i7)') ' $OMP = ',numdomp

      write(ioptmp,'(a19)') 'infl: {varfile.inp}'

*-----------------------------------------------------------------------
*     make varfile
*-----------------------------------------------------------------------

      call make_varfile(nfile,iorder,ierr,ioanat)
      if ( ierr .ne. 0 ) return

*-----------------------------------------------------------------------
*     reset initial values for include file
*-----------------------------------------------------------------------

      rewind jsi
      do i = 0, 9
         ill(i) = 0
         ilf(i) = 10000000
         ild(i) = 0
      end do
      jpn  = 0

*-----------------------------------------------------------------------
*     read and write to make phits_tmp.inp
*-----------------------------------------------------------------------

      do while ( jpn .ne. 3 ) ! search [parameters]

       if ( jpn .ne. 1 ) then
        call fnd_setc(jsi,ioptmp,ioanat)
        call readl(jsn,jsi,dsin,idsi,ill,ilf,'%!$',
     &       jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)
       end if

       if (  jpn.ne.3 .and. iskip.eq.0 ) then
        if( chcm(i1:i1+11) .eq. '[parameters]' ) then
         if( chcm(i1+12:i1+14) .eq. 'off' ) then
          jpn = 2

         else ! when available [parameters] is found
           write(ioptmp,'(/200a1)') (chin(i:i),i=1,i2)
           write(ioanat,'(/200a1)') (chin(i:i),i=1,i2)

          do while ( jpn .ne. 3 ) ! search icntl

           call fnd_setc(jsi,ioptmp,ioanat)
           call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &            jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

           if (  jpn.ne.3 .and. iskip.eq.0 ) then
*     end of [parameters]
            if( i1 .le. 5 .and. chlw(i1:i1) .eq. '[' ) then
             jpn = 1
             exit

*     identify the parameters
            else
             icl = i1
             chlc = chlw
             call chcomp(chlc,icl,i3,i5)

             if( chlc(icl:icl+4) .eq. 'icntl' ) then
              write(ioptmp,'(a9)') 'icntl = 0'
              write(ioanat,'(a10)') 'icntl = 17'
              icl = inumc(chlc,icl+4,i5,';') - 1
              icl = jnumc(chlc,icl+2,i3)
              if( icl .le. i5 ) then ! when ; is included in the line of icntl
               i1 = inumc(chlw,i1,i3,';')+1
               write(ioptmp,'(200a1)') (chlw(i:i),i=i1,i3)
               write(ioanat,'(200a1)') (chlw(i:i),i=i1,i3)
              end if

             else ! except icntl
              write(ioptmp,'(200a1)') (chin(i:i),i=1,i2)
              write(ioanat,'(200a1)') (chin(i:i),i=1,i2)
             end if
            end if
           end if

          end do
         end if

        else if( chcm(i1:i1+2) .eq. '[t-' ) then ! tally is found
         if( chcm(i4-2:i4) .eq. 'off' ) then
          jpn = 2

         else ! when available tally section is found
          jpn = 0
          ianataldchaint = 0
          if ( chcm(i1:i1+8) .eq. '[t-dchain' ) ianataldchaint = 1
          write(ioptmp,'(/200a1)') (chin(i:i),i=1,i2)
          write(ioanat,'(/200a1)') (chin(i:i),i=1,i2)
          itnmw = itnmw + 1
          do while ( jpn .ne. 1 ) ! until the next section is found

           call fnd_setc(jsi,ioptmp,ioanat)
           call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &            jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

           if (  jpn.ne.3 .and. iskip.eq.0 ) then
            if( chcm(i1:i1+12) .eq. 'anatallystart' ) then
             write(ioanat,'(200a1)') (chlw(i:i),i=i1,i3)
             write(chlw(i3+1:i3+5),'(i5)') itnmw
             i3 = i3 + 5
             write(ioanattmp,'(200a1)') (chlw(i:i),i=i1,i3)
             if ( ianataldchaint .eq. 1 )
     &            write(ioanattmp,'(a17)') 'ianataldchain = 1'
             ianatalend = 0
             do while ( ianatalend .eq. 0 )
              call fnd_setc(jsi,ioptmp,ioanat)
              call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &               jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)
              write(ioanattmp,'(200a1)') (chlw(i:i),i=i1,i3)
              if( chcm(i1:i1+10) .eq. 'anatallyend' ) then
               ianatalend = 1
               backspace(ioanattmp)
               write(ioanattmp,'(a9,i7)') ' nfile = ',nfile
               do ifile = 1, nfile
                cunder(1:5) = '_'//'0000'
                iread = ifile
                write(cnumber,*) iread
                call chcomp(cnumber,1,20,i4)
                cunder(iorder+2-i4:iorder+2) = cnumber(1:i4)
                write(ioanattmp,'(a10,a,a1,a)')
     &               ' outfiles/',cunder(2:iorder+1)
     &               ,'/',ctallyfile(1:ltallyfile)
               end do
               write(ioanat,'(200a1)') (chlw(i:i),i=i1,i3)
               write(chlw(i3+1:i3+5),'(i5)') itnmw
               i3 = i3 + 5
               write(ioanattmp,'(200a1)') (chlw(i:i),i=i1,i3)
              end if
             end do

*     end of anatally subsection
            else if( i1 .le. 5 .and. chlw(i1:i1) .eq. '[' ) then
             jpn = 1

            else
             if( chcm(i1:i1+4) .eq. 'file=' ) then
                ltallyfile = i3-i1
                ctallyfile(1:ltallyfile) = chcm(i1+5:i1+5+ltallyfile)
             end if
             write(ioptmp,'(200a1)') (chlw(i:i),i=i1,i3)
             write(ioanat,'(200a1)') (chlw(i:i),i=i1,i3)

            end if

           end if
          end do

         end if


        else if( chcm(i1:i1+4) .eq. '[end]' ) then ! [end] is found
         if( chcm(i4-2:i4) .eq. 'off' ) then
          jpn = 2

         else ! write [anatally] section before [end]
          write(ioptmp,'(/200a1)') (chin(i:i),i=1,i2)
          rewind(ioanattmp)
          write(ioanat,'(/a10)') '[anatally]'
          ios = 1
          do while ( ios .ge. 0 )
           read(ioanattmp,'(a200)',iostat=ios) chin
           if ( ios .lt. 0 ) exit
           call chlngt(chin,200,i1,i2)
           chlw = chin
           call chcaps(chlw,i1,i2,i3,'#!$')
           write(ioanat,'(200a1)') (chlw(i:i),i=i1,i3)
          end do
           write(ioanat,'(/a5)') '[end]'
           jpn = 3

         end if


        else ! the other section is found

         jpn = 0
         if ( chcm(i1:i1) .eq. '[' ) then
          write(ioptmp,'(/200a1)') (chin(i:i),i=1,i2)
          write(ioanat,'(/200a1)') (chin(i:i),i=1,i2)
         else
          write(ioptmp,'(200a1)') (chin(i:i),i=1,i2)
          write(ioanat,'(200a1)') (chin(i:i),i=1,i2)
         end if

        end if
       end if

      end do


*-----------------------------------------------------------------------
      return

*-----------------------------------------------------------------------

 997  continue

      m_err = 'Unknown section name'
      ErrCha = ''
      ErrID = 'L:441/R:anainp/F:anainp.f'
      l_err = ill(jsn)
      k_err = jsn

      goto 999

*-----------------------------------------------------------------------

 998  continue

      dsin(0) = 'Error Line'
      idsi(0) = 12
      m_err = 'There is nothing in the normal input'
      ErrCha = ''
      ErrID = 'L:455/R:anainp/F:anainp.f'
      l_err = 1
      k_err = 0
      goto 999

*-----------------------------------------------------------------------

 999  continue


      write(*,'(/'' ***** Error Message from Input File *****''/)')
      call ErrLine_Adjust(dsin(k_err)(1:idsi(k_err)),l_err)
      write(*,*) dsin(k_err)(1:idsi(k_err)),l_err,':'

      icf = 200

      do 910 i = 200, 1, -1
         if( m_err(i:i) .ne. ' ' ) goto 911
 910  continue

 911  icf = i
      call ErrWrite(ErrID, ErrCha)
      write(*,'('' error = '',200A1/)') ( m_err(i:i), i=1,icf )

      ierr = 1

*-----------------------------------------------------------------------

      return
      end subroutine anainp


!***********************************************************************
!                                                                      *
      subroutine make_varfile(nfile,iorder,ierr,ioanat)
!                                                                      *
*-----------------------------------------------------------------------

      implicit none

*-----------------------------------------------------------------------

      integer nfile, ierr, ionc, iovar, ioanat
      character cncinfofile*11
      data cncinfofile / 'nc_info.inp' /
      logical exex
      integer nvar, ivar

*-----------------------------------------------------------------------

      character chin*200, chlw*200, chcm*200
      character chlc*200
      character dsin(0:9)*200
      integer idsi(0:9)
      integer ill(0:9), ilf(0:9)
      integer ios,i1,i2,i3,i4,ic,ic2,inumc
      real*8 cvvv

      integer iread,icv

      integer maxnvar, maxcvar
      parameter (maxnvar=9, maxcvar=10000)
      integer iduc(maxnvar),ictp(maxnvar),inc(maxnvar)
      real*8 cmin(maxnvar),cmax(maxnvar),cdel(maxnvar)
      real*8 cvalue(maxcvar,maxnvar)

      common /cnvar/ nvar,iduc
      common /ccvar/ cvalue ! S.H. 2021.10.5

*-----------------------------------------------------------------------

      character cvarfile*7
      data cvarfile / 'varfile' /
      character cunder*5,cnumber*20
      integer iorder

*-----------------------------------------------------------------------

      integer icntl, inucr
      common /tcntl/  icntl, inucr ! S.H. 2021.10.5

*-----------------------------------------------------------------------

      ionc = 796
      iovar = 795

*-----------------------------------------------------------------------
*     read nc_info_inp
*-----------------------------------------------------------------------

      inquire( file = cncinfofile, exist = exex )

      if( exex .eqv. .false. ) then

         nvar = 1
         iduc(1) = 1
         ictp(1) = 1
         inc(1) = 1
         cmin(1) = 0
         cmax(1) = 0
         cdel(1) = 0
         cvalue(1,1) = 0

      else

       open(ionc, file=cncinfofile, status='old')
       read(ionc,*) nvar

       do ivar = 1, nvar

        read(ionc,'(a200)', iostat=ios ) chin
        call chlngt(chin,200,i1,i2)
        chlw = chin
        call chcaps(chlw,i1,i2,i3,'#!$')
        chcm = chlw
        call chcomp(chcm,i1,i3,i4)
        if(chcm(i1:i1+3).eq.'set:')read(chcm(i1+5:i1+7),*)iduc(ivar)

        read(ionc,'(a200)', iostat=ios ) chin
        call chlngt(chin,200,i1,i2)
        chlw = chin
        call chcaps(chlw,i1,i2,i3,'#!$')
        chcm = chlw
        call chcomp(chcm,i1,i3,i4)
        if( chcm(i1:i1+6) .eq. 'c-type=' ) then
         ic = inumc(chlw,i1,i3,'=') + 1
         call snum(chlw,ic,i3,ic2,cvvv,ierr)
         if( ierr .ne. 0 ) goto 999
         ictp(ivar) = nint( cvvv )
         if ( ictp(ivar) .eq. 1 ) then ! c-type=1
          read(ionc,'(a200)', iostat=ios ) chin
          call chlngt(chin,200,i1,i2)
          chlw = chin
          call chcaps(chlw,i1,i2,i3,'#!$')
          chcm = chlw
          call chcomp(chcm,i1,i3,i4)
          if(chcm(i1:i1+2).eq.'nc=') read(chcm(i1+3:i4),*)inc(ivar)
          do iread=1,inc(ivar)
             read(ionc,*) cvalue(iread,ivar)
          end do

         else ! c-type=2,3,4,5
          do iread=1,3
           read(ionc,'(a200)', iostat=ios ) chin
           call chlngt(chin,200,i1,i2)
           chlw = chin
           call chcaps(chlw,i1,i2,i3,'#!$')
           chcm = chlw
           call chcomp(chcm,i1,i3,i4)
           if(chcm(i1:i1+2).eq.'nc=') then
            read(chcm(i1+3:i4),*)inc(ivar)
           else if(chcm(i1:i1+4).eq.'cmin=') then
            read(chcm(i1+5:i4),*) cmin(ivar)
           else if(chcm(i1:i1+4).eq.'cmax=') then
            read(chcm(i1+5:i4),*) cmax(ivar)
           else if(chcm(i1:i1+4).eq.'cdel=') then
            read(chcm(i1+5:i4),*) cdel(ivar)
           end if
          end do

          if ( ictp(ivar) .eq. 2 ) then ! c-type=2
           if ( inc(ivar) .gt. 1 ) then
            cdel(ivar) = (cmax(ivar)-cmin(ivar))/(inc(ivar)-1)
           else
            cdel(ivar) = 0
           end if
           do icv=1,inc(ivar)
              cvalue(icv,ivar) = cmin(ivar) + (icv-1)*cdel(ivar)
           end do

          else if ( ictp(ivar) .eq. 3 ) then ! c-type=3
           if ( inc(ivar) .gt. 1 ) then
            cdel(ivar) = dlog(cmax(ivar)/cmin(ivar))/(inc(ivar)-1)
           else
            cdel(ivar) = 0
           end if
           do icv=1,inc(ivar)
            cvalue(icv,ivar)
     &             = cmin(ivar) * dexp((icv-1)*cdel(ivar))
           end do

          else if ( ictp(ivar) .eq. 4 ) then ! c-type=4
           if ( cdel(ivar) .gt. 0d0 ) then
            inc(ivar) = (cmax(ivar)-cmin(ivar))/cdel(ivar) +1
           else
            inc(ivar) = 1
           end if
           do icv=1,inc(ivar)
              cvalue(icv,ivar) = cmin(ivar) + (icv-1)*cdel(ivar)
           end do

          else if ( ictp(ivar) .eq. 5 ) then ! c-type=5
           if ( cdel(ivar) .gt. 0d0 ) then
            inc(ivar)=idnint(dlog((cmax(ivar)/cmin(ivar)))/cdel(ivar))+1
           else
            inc(ivar) = 1
           end if
           do icv=1,inc(ivar)
            cvalue(icv,ivar)
     &             = cmin(ivar) * dexp((icv-1)*cdel(ivar))
           end do

          end if

         end if
        end if

       end do

      end if

      if ( icntl .eq. 17 ) return ! S.H. 2021.10.5

*-----------------------------------------------------------------------
*     make varfile
*-----------------------------------------------------------------------

      nfile = 1
      do ivar = 1, 1!nvar

         nfile = nfile*inc(ivar)

       if ( inc(ivar).ge.1 .and. inc(ivar).lt.10 ) then
          iorder = 1
       else if ( inc(ivar).ge.10 .and. inc(ivar).lt.100 ) then
          iorder = 2
       else if ( inc(ivar).ge.100 .and. inc(ivar).lt.1000 ) then
          iorder = 3
       else if ( inc(ivar).ge.1000 .and. inc(ivar).lt.10000 ) then
          iorder = 4
       else
          ierr = 1
          goto 999
       end if
       do icv=1,inc(ivar)
          cunder(1:5) = '_'//'0000'
          iread = icv
          write(cnumber,*) iread
          call chcomp(cnumber,1,20,i4)
          cunder(iorder+2-i4:iorder+2) = cnumber(1:i4)
          open(iovar,
     &      file=cvarfile//cunder(1:iorder+1)//'.inp', status='unknown')

          if ( iduc(ivar) .lt. 10 ) then
             write(iovar,'(a5,i1,a1,1p1e25.16,a1)')
     &            'set:c',iduc(ivar),'[',cvalue(icv,ivar),']'
          else if ( iduc(ivar).ge.10 .and. iduc(ivar).lt.100 ) then
             write(iovar,'(a5,i2,a1,1p1e25.16,a1)')
     &            'set:c',iduc(ivar),'[',cvalue(icv,ivar),']'
          else
             write(iovar,'(a5,i3,a1,1p1e25.16,a1)')
     &            'set:c',iduc(ivar),'[',cvalue(icv,ivar),']'
          end if
       end do
          if ( iduc(ivar) .lt. 10 ) then
             write(ioanat,'(a5,i1,a1,1p1e25.16,a1)')
     &            'set:c',iduc(ivar),'[',cvalue(1,ivar),']'
          else if ( iduc(ivar).ge.10 .and. iduc(ivar).lt.100 ) then
             write(ioanat,'(a5,i2,a1,1p1e25.16,a1)')
     &            'set:c',iduc(ivar),'[',cvalue(1,ivar),']'
          else
             write(ioanat,'(a5,i3,a1,1p1e25.16,a1)')
     &            'set:c',iduc(ivar),'[',cvalue(1,ivar),']'
          end if
      end do

*-----------------------------------------------------------------------

 999  continue
      return

      end subroutine make_varfile
!***********************************************************************


!***********************************************************************
!                                                                      *
      subroutine fnd_setc(jsi,ioptmp,ioanat)
!                                                                      *
*-----------------------------------------------------------------------

      implicit none

      integer jsi, ios, ioptmp,ioanat
      character chin*200, chlw*200, chcm*200
      integer i1,i2,i3,i4, isetcommand, i, iskip

      integer nvar,iduc(9), ivar,iducread, ihitiduc
      common /cnvar/ nvar,iduc


*-----------------------------------------------------------------------

      isetcommand = 1
      do while ( isetcommand .eq. 1 )

         read(jsi,'(a200)', iostat=ios ) chin
         call ZspcReplace(chin)

         call chlngt(chin,200,i1,i2)
         chlw = chin
         call chcaps(chlw,i1,i2,i3,'#!$')
         chcm = chlw
         call chcomp(chcm,i1,i3,i4)
         if( i1 .eq. 0 .and. i2 .eq. 0 ) then
            iskip = 1
         else if( index('#!$',chin(i1:i1)) .ne. 0 ) then
            iskip = 2
         else
            iskip = 0
         end if
         if( iskip .eq. 0 ) then
         if( chcm(i1:i1+3) .eq. 'set:' ) then
            if ( chcm(i1+6:i1+6) .eq. '[' ) then
               read(chcm(i1+5:i1+5),*) iducread
            else if ( chcm(i1+7:i1+7) .eq. '[' ) then
               read(chcm(i1+5:i1+6),*) iducread
            else if ( chcm(i1+8:i1+8) .eq. '[' ) then
               read(chcm(i1+5:i1+7),*) iducread
            end if
            ihitiduc = 0
            do ivar = 1, nvar
               if ( iducread .eq. iduc(ivar) ) ihitiduc=1
            end do
            if ( ihitiduc .eq. 1 ) then
               chlw(i1:i3+1) = '$'//chlw(i1:i3)
               write(ioptmp,'(200a1)') (chlw(i:i),i=i1,i3+1)
               write(ioanat,'(200a1)') (chlw(i:i),i=i1,i3+1)
            else
               write(ioptmp,'(200a1)') (chlw(i:i),i=i1,i3)
               write(ioanat,'(200a1)') (chlw(i:i),i=i1,i3)
            end if

         else
            isetcommand = 0
         end if

         end if ! S.H. 2020.7.14

      end do
      backspace(jsi)

*-----------------------------------------------------------------------

      return

      end subroutine fnd_setc
!***********************************************************************
