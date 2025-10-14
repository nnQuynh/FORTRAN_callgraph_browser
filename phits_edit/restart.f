************************************************************************
*                                                                      *
      subroutine prestart(m,iot)
*                                                                      *
*   print restart data, independ from tally.                           *
*                                                                      *
*   m:   [in] tally number                                             *
*   iot: [in] tally file number                                        *
************************************************************************

        use sumtallymod, only : isumtally

        implicit double precision (a-h,o-z)

        include 'param.inc'

*-----------------------------------------------------------------------

        common /taliin/ rsouin, nzztin, nrgnin
        common /cparm/  maxbch,maxcas
        common /jcomon/ nabov,nobch,nocas,nomax
!$OMP THREADPRIVATE(/jcomon/)
        common /restart/ resc2(itlmax), resc3(itlmax)
        common /stat / istdev, irestart, ireschk
        common /randm5/ rijklst, rijkinit
        integer nrandgen ! S.H. xorshift (2020.2.6)
        common /randn/ nrandgen

        common /tcntl/ icntl, inucr
        common /res01/ istdevres,maxcasres,rijklstres,irdrf

        common /sumtal01/ isistdev, ismaxcas, rsijklst
        integer   :: isistdev(itlmax), ismaxcas(itlmax)
        dimension rsijklst(itlmax)

        common /talout/ itall
        common /mpi00/ npe, me

*-----------------------------------------------------------------------

        common /tall59/ irfll(itlmax), crfln(itlmax), itrff(itlmax)
        character crfln*100

        common /res02/ crdrfln, irdrfll
        character crdrfln*100

        common /res03/ lextrf(itlmax)
        logical lextrf

*-----------------------------------------------------------------------

        if ( icntl .eq. 13 ) then

           istdev = isistdev(m)
           maxcas = ismaxcas(m)
           rijklst = rsijklst(m)

        end if

*-----------------------------------------------------------------------

        if ( istdev .eq. 2 ) then
          c2 = rsouin + resc2(m)
          c3 = dble(nobch) * dble(maxcas) + resc3(m)
        else
          c2 = rsouin / maxcas + resc2(m)
          c3 = dble(nobch) + resc3(m)
        end if

        if ( itall.eq.4 .and. istdev.eq.2 ) then
          if ( npe .gt. 1 ) then
            c2 = rsouin/dble(nobch/(npe-1))+resc2(m)
            c3 = dble(npe-1)*dble(maxcas)+resc3(m)
          else
            c2 = rsouin/dble(nobch)+resc2(m)
            c3 = dble(maxcas)+resc3(m)
          end if
        end if

*-----------------------------------------------------------------------

        if ( icntl == 13 .and.
     &       (isumtally(m)==2 .or. isumtally(m)==3) ) then ! Weighted average
           write(iot,*)                                       ! blank line
           goto 9999
        end if

        write(iot,*) ! brank line
        write(iot,'("# Information for Restart Calculation")')

        if( irestart .ne. 0 .and. lextrf(m) ) then
            write(iot,
     &        '("# This calculation was restarted using ",a)')
     &        crfln(m)(1:irfll(m))
        else
            write(iot,'("# This calculation was newly started")')
        end if

        write(iot,'("# istdev =",i3,
     &  " # 1:Batch variance, 2:History variance")') istdev
        write(iot,'("# resc2  =",1pe24.17,
     &  " # Total source weight or Total source weight / maxcas")') c2
        write(iot,'("# resc3  =",1pe24.17,
     &  " # Total history number or Total batch number")')  c3

        write(iot,'("# maxcas =",i12,
     &  " # History / Batch, only used for istdev=1")')  maxcas
        !FURUTA_2012/0809
       if ( nrandgen .eq. 0 ) then ! when LCG (2021.4.20)
        write(iot,'("# rseed = ",1p1e25.16e3,
     &" # Next initial random number")') rijklst
       else ! when xorshift
        write(iot,'("# bitrseed = ",b64.64," # bit data of rseed")')
     &         rijklst
       end if

 9999   continue

*-----------------------------------------------------------------------

      end subroutine


************************************************************************
*                                                                      *
      subroutine rrestart(jsn,jsi,dsin,idsi,ill,ilf,
     &                    istdev,resc2,resc3,maxcas,rijklst,ierr)
*                                                                      *
*                                                                      *
************************************************************************

        implicit double precision (a-h,o-z)

        include 'err.inc'

        character dsin(0:9)*200
        dimension idsi(0:9)
        dimension ill(0:9), ilf(0:9)

*-----------------------------------------------------------------------

        character chin*200, chlw*200, chcm*200

*-----------------------------------------------------------------------

      integer*8 :: ibitrseed ! S.H. xorshift (2020.2.6)
      common /paraj/  mstz(300), parz(300)
      integer nrandgen
      common /randn/ nrandgen

*-----------------------------------------------------------------------

        parameter(icsu = 6) ! S.H. xorshift (2020.2.6)
        dimension lschn(icsu), ischn(icsu)
        character schan(icsu)*10
        character stemp*20
        logical   lresprm(icsu)

        data ( schan(i), i = 1, icsu ) /
     &    'istdev ',
     &    'resc2  ',
     &    'resc3  ',
     &    'maxcas ',
     &    'rseed',
     &    'bitrseed' ! S.H. xorshift (2020.2.6)
     &  /

        data ( lschn(i), i = 1, icsu ) / 6, 5, 5, 6, 5, 8 / ! S.H. xorshift (2020.2.6)

*-----------------------------------------------------------------------

        jpn   = 0
        i1    = 0
        i2    = 0
        i3    = 0
        i4    = 0
        iskip = 0
        ierr  = 0
        prn   = 0.d0

        lresprm(:) = .false.
        lresprm(6) = .true. ! S.H. xorshift (2020.2.6)
        ierr  = 0

*-----------------------------------------------------------------------
* restart file format:
*-----------------------------------------------------------------------
* # Information for Successive Calculation
* # istdev = i3     # 1:Batch variance, 2:History variance
* # resc2  = e24.17 # Total source weight or Total source weight / maxcas
* # resc3  = e24.17 # Total history number or Total batch number
* # maxcas = i,     # History / Batch, only used for istdev=1
* # rseed  = f17.1  # Next initial random number
* # bitrseed=b64    # bit data of rseed
*-----------------------------------------------------------------------

  100 continue

        call readl(jsn,jsi,dsin,idsi,ill,ilf,'%!$',
     &             jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

        if ( ierr  .ne. 0 ) goto 900
        if ( jpn   .eq. 3 ) goto 900
        if ( iskip .ne. 0 ) goto 100

        do i = 1, icsu

          stemp = '#' // schan(i)(1:lschn(i)) // '='

          il = i1 + lschn(i) + 1
          if ( chcm(i1:il) .eq. stemp(1:lschn(i)+2) ) then

            ivl = i4
            do k = il+1, i4
              if ( chcm(k:k) .eq. '#' ) then
                ivl = k - 1
                exit
              end if
            end do
            call onum(chcm,il+1,ivl,prn,ierr)

            lresprm(i) = .true.
*-----------------------------------------------------------------------
* read istdev ( i = 1 )
*-----------------------------------------------------------------------
            if ( i .eq. 1 ) then

              istdev = int(prn)

*-----------------------------------------------------------------------
* read resc2 ( i = 2 )
*-----------------------------------------------------------------------
            else if ( i .eq. 2 ) then

              resc2 = prn

*-----------------------------------------------------------------------
* read resc3 ( i = 3 )
*-----------------------------------------------------------------------
            else if ( i .eq. 3 ) then

              resc3 = prn

*-----------------------------------------------------------------------
* read maxcas ( i = 4 )
*-----------------------------------------------------------------------
            else if ( i .eq. 4 ) then

              maxcas = int(prn)

*-----------------------------------------------------------------------
* read rseed ( i = 5 )
*-----------------------------------------------------------------------
            else if ( i .eq. 5 ) then

              if ( mstz(139) .eq. 1 ) then
               write(ErrCha,'("Error: nrandgen should be set to the ",
     &               "same as that in the past calculation. ",
     &               "Current nrandgen = ",i3 )')
     &                mstz(139)
                 ErrID = 'L:261/R:rrestart/F:restart.f'
                 call ErrWrite(ErrID,ErrCha)
                 ierr = 1
                 return
              end if

              rijklst = prn

c S.H. xorshift (2020.2.6)
*-----------------------------------------------------------------------
* read rijklst ( i = 6 )
*-----------------------------------------------------------------------
            else if ( i .eq. 6 ) then

              if ( mstz(139) .ne. 1 ) then
               write(ErrCha,'("Error: nrandgen should be set to the ",
     &               "same as that in the past calculation. ",
     &               "Current nrandgen = ",i3 )')
     &                mstz(139)
                 ErrID = 'L:280/R:rrestart/F:restart.f'
                 call ErrWrite(ErrID,ErrCha)
                 ierr = 1
                 return
              end if

              read(chcm(il+1:ivl),'(b64)') ibitrseed
              rijklst = transfer(ibitrseed,rijklst)
              lresprm(5) = .true.

*-----------------------------------------------------------------------
            endif
            exit

          endif

        enddo

        goto 100

*-----------------------------------------------------------------------
* end of file
*-----------------------------------------------------------------------
  900 continue

*-----------------------------------------------------------------------
* check found all parameter
*-----------------------------------------------------------------------
      ierr = 0
      if( .not. all( lresprm(:) )) then
        write(ErrCha,'("Error: No restart parameters in ", a )')
     &        dsin(0)(1:idsi(0))
       ErrID = 'L:312/R:rrestart/F:restart.f' !E83_001_001
       call ErrWrite(ErrID,ErrCha)

        ierr = 1
      end if

      end subroutine


************************************************************************
*                                                                      *
      subroutine read_resfiles(io,jo,ivers,ierr)
*                                                                      *
*                                                                      *
************************************************************************

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

*-----------------------------------------------------------------------

        character dsin(0:9)*200
        dimension idsi(0:9)
        dimension ill(0:9), ilf(0:9)

*-----------------------------------------------------------------------
        common /mpi00/ npe, me
        common /restbch/ rescb2(itlmax,isrc),rescb3(itlmax,isrc)
        common /itl4res/ itl4flg
        integer :: itl4flg = 0
        character*100 fname
        character*3 fnume(itlmax)
        logical l_exists1

*-----------------------------------------------------------------------

        iristdev(:) = -1
        irmaxcas(:) = maxcas
        rrijklst(:) = 0.d0

        lextrf(:)    = .false.

        inoprm = 0

        irestartable = 0

*-----------------------------------------------------------------------
* check icntl
*-----------------------------------------------------------------------
        if ( mstz(1) .eq. 8 ) then  ! mstz(1) is 'icntl' parameter.
          ErrCha = ''
          MsgID = 'L:398/R:read_resfiles/F:restart.f'
          call ErrWrite(MsgID, ErrCha)
          write(jo,'("istdev parameter is ignored in icntl = 8")')
          irestart = 0 ! cancel restart, it be newly started at all talls.
          return
        end if

*-----------------------------------------------------------------------
* read each restart files
*-----------------------------------------------------------------------
        do 100 m = 1, itnm

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
     &         18, 19, 21 /) ) ) then
            goto 100

          end if

          irestartable = irestartable + 1

*-----------------------------------------------------------------------
* read restart prameters in read resfile
*-----------------------------------------------------------------------
          if ( irfll(m) .eq. 0 ) then
            irfll(m) = itfll(m,1)
            crfln(m)(1:) = ctfln(m,1)(1:itfll(m,1))
          end if

          if( mstz(26).eq.4 ) itl4flg = 1
          mxbcht = mstz(4)
          if ( itl4flg.eq.1 ) then
            do itt = 1, mxbcht
              call mk_2dnumfn(crfln(m),fname,irfll(m),
     &             itt,mxbcht,npe)
              inquire(file=fname,exist=l_exists1)
              if ( l_exists1 ) then
                idsi(jsn) = irfll(m)+4
                dsin(jsn)(1:idsi(jsn)) = fname
                open( jsi, file=dsin(jsn)(1:idsi(jsn)),
     &                status='old', action='read',
     &                form='formatted', iostat=ierr )
                if ( ierr.eq.0 ) then
                  call rrestart( jsn,jsi,dsin,idsi,ill,ilf,
     &                           iristdev(m),rescb2(m,itt),
     &                           rescb3(m,itt),irmaxcas(m),
     &                           rrijklst(m),ierr )
                  jsn = jsn + 1
                  close(jsi)
                  if ( ierr.eq.0 ) then
                    lextrf(m) = .true.
                  else
                    inoprm = inoprm + 1
                  end if
                else
                  rescb2(m,itt) = 0.d0
                  rescb3(m,itt) = 0.d0
                end if
              else
                rescb2(m,itt) = 0.d0
                rescb3(m,itt) = 0.d0
              end if
            end do
          resc2(m) = rescb2(m,1)
          resc3(m) = rescb3(m,1)
          goto 100
          end if

          idsi(jsn) = irfll(m)
          dsin(jsn)(1:idsi(jsn)) = crfln(m)(1:idsi(jsn))

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
          MsgID = 'L:516/R:read_resfiles/F:restart.f'
          call ErrWrite(MsgID, ErrCha)
          write(jo,
     &      '("Error: Some restart file have no restart parameter.")')
          ierr = 1
          return
        end if

*-----------------------------------------------------------------------
* check exit restartable tally
*-----------------------------------------------------------------------
        if ( irestartable .eq. 0 ) then
          irestart = 0 ! cancel restart, it be newly started at all talls.
          return
        end if

*-----------------------------------------------------------------------
* check exit restart files
*-----------------------------------------------------------------------
        if ( .not. any( lextrf(1:itnm) ) ) then
          !! all restart files are not found
          ErrCha = ''
          MsgID = 'L:538/R:read_resfiles/F:restart.f'
          call ErrWrite(MsgID, ErrCha)
          write(jo,'("Warning: All resfile(s) do not exist, ",
     &              "so new calculation started with istdev = ",
     &              i1 )') istdev
          irestart = 0 ! cancel restart, it be newly started at all talls.
          return
        end if

*-----------------------------------------------------------------------
* tall number of primery restart file
*-----------------------------------------------------------------------
        irdrf = 0
        do m = 1, itnm
          if ( lextrf(m) ) then
            irdrf = m
            exit
          end if
        end do

        if ( .false. ) then !! debug
          write (0,*) iristdev(1:itnm)
          write (0,*) irmaxcas(1:itnm)
          write (0,*) rrijklst(1:itnm)
          write (0,*) lextrf(1:itnm)
        end if

*-----------------------------------------------------------------------
* left align values read from restart file
*-----------------------------------------------------------------------
        im = 0
        do m = 1, itnm
          if ( lextrf(m) ) then
            im = im + 1
            iristdev(im) = iristdev(m)
            irmaxcas(im) = irmaxcas(m)
            rrijklst(im) = rrijklst(m)
          end if
        end do

*-----------------------------------------------------------------------
* check all read istdev values are same
*-----------------------------------------------------------------------
        if( itnm .gt. 1 .and. im .gt. 1 ) then
          if( any(iristdev(1) .ne. iristdev(2:im)) ) then
            !! invalid istdev
            ErrCha = ''
            MsgID = 'L:585/R:read_resfiles/F:restart.f'
            call ErrWrite(MsgID, ErrCha)
            write(jo,
     &        '("Error: istdev should be the same for all resfiles.",
     &          " Check your restart files.")')
            ierr = 1
            return
          end if
        end if

*-----------------------------------------------------------------------
* check all read maxcas values are same
*-----------------------------------------------------------------------
        if( ireschk .eq. 0 .and. iristdev(1) .eq. 1) then ! T.Sato 2015/05/08
         if( itnm .gt. 1 .and. im .gt. 1 ) then
           if( any(irmaxcas(1) .ne. irmaxcas(2:im)) ) then
             !! invalid maxcas
             ErrCha = ''
             MsgID = 'L:603/R:read_resfiles/F:restart.f'
             call ErrWrite(MsgID, ErrCha)
             write(jo,
     &       '("Error: maxcas should be the same for all resfiles.",
     &       " Check your restart files.")')
             ierr = 1
             return
           end if
         end if
        end if
*-----------------------------------------------------------------------
* check all read rijklst values are same
*-----------------------------------------------------------------------
        lrijkeqrf = all(rrijklst(1) .eq. rrijklst(2:im))

*-----------------------------------------------------------------------
* restore values read from restart file
*-----------------------------------------------------------------------
        istdevres  = iristdev(1)
        maxcasres  = irmaxcas(1)
        rijklstres = rrijklst(1)

        irdrfll = irfll(irdrf)
        crdrfln(1:irdrfll) = crfln(irdrf)(1:irdrfll)
        if ( itl4flg.eq.1 ) then
          irdrfll = irfll(irdrf)+4
          mxbcht = mstz(4)
          call mk_2dnumfn(crfln(irdrf),fname,irfll(irdrf),
     &         mxbcht,mxbcht,npe)
          crdrfln(1:irdrfll) = fname
        end if

        if ( .false. ) then !! debug
          write (0,*) irdrf
          write (0,*) crdrfln(1:irdrfll)
          write (0,*) istdevres
          write (0,*) maxcasres
          write (0,*) resc2(1:itnm)
          write (0,*) resc3(1:itnm)
          write (0,*) rijklstres
          write (0,*) lrijkeqrf
        end if

      end subroutine


************************************************************************
*                                                                      *
      subroutine setrnd_restart()

*   initialize random for restart.                                     *
************************************************************************

        implicit real*8 (a-h,o-z)

        include 'param.inc'

*-----------------------------------------------------------------------

        parameter (p=2d0**24,q=2d0**(-24),rm=5d0**19)

*-----------------------------------------------------------------------

        integer nrandgen ! S.H. xorshift (2020.2.6)
        common /randn/ nrandgen
        integer*8 :: iranji64 ! S.H. xorshift (2020.2.6)
        common /randm4/ rnfb,rnfs,rngb,rngs,rnmult,ranj,rani,
     &                  rnrtc,nstrid,inif, iranji64
        common /randm5/ rijklst, rijkinit
        integer*8 :: iransb64 ! S.H. xorshift (2020.2.6)
        common /randtp/ rijk,rans,ranb, iransb64
!$OMP THREADPRIVATE(/randtp/)

        common /iradkk/ randkk,irskip
        common /randsv/ srijk,rrijk

        common /stat / istdev, irestart, ireschk

        common /res01/ istdevres,maxcasres,rijklstres,irdrf

*-----------------------------------------------------------------------

        if( irestart .ne. 0 ) then
          rijklst  = rijklstres
          rijkinit = rijklstres
          rijk = rijkinit
c S.H. xorshift (2020.2.6)
         if ( nrandgen .eq. 0 ) then
          rani=aint(rijk*q)
          ranj=rijk-rani*p
         else
          iranji64 = transfer(rijk,iranji64)
         end if
        else
          rijkinit = rijk
        end if

      end subroutine


************************************************************************
*                                                                      *
      subroutine setresval(io,jo,ivers,ierr)
*                                                                      *
*       programed by OBINATA on 2012/06/15                             *
*                                                                      *
************************************************************************
!$    use omp_lib

      implicit real*8 (a-h,o-z)

      include 'param.inc'

*-----------------------------------------------------------------------

      common /mpi00/ npe, me
      common /stat / istdev, irestart, ireschk
      common /cparm/ maxbch,maxcas

*-----------------------------------------------------------------------

      common /res01/ istdevres,maxcasres,rijklstres,irdrf
      integer italsh
      common /talsh/ italsh

*-----------------------------------------------------------------------

      if( irestart .ne. 0 ) then
        istdev = istdevres
c  Changed by T.Sato 2012/11/22
        if (istdev .eq. 1) maxcas = maxcasres
      end if
      ierr = 0
      if( me .eq. 0 ) then
!$      if( istdev .ne. 1 .and. omp_get_max_threads() .gt. 1) then
!$      if( italsh .eq. 1 )then !FURUTA20200907
!$        write(jo,
!$   &       '(''Error: istdev=1 or italsh=0 should be chosen'',
!$   &         '' on multi-thread execution.'')')
!$        ierr = 1
!$      end if
!$      endif
      end if

      end subroutine



************************************************************************
* read_talls:                                                          *
*                                                                      *
************************************************************************
      subroutine read_talls(io,ierr)
        use GGBANKMOD

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
        common /mpi00/ npe, me
        common /jcomon/ nabov,nobch,nocas,nomax
!$OMP   THREADPRIVATE(/jcomon/)
        common /restbch/ rescb2(itlmax,isrc),rescb3(itlmax,isrc)
        common /itl4res/ itl4flg

*-----------------------------------------------------------------------

        ierr    = 0

        if ( itl4flg.eq.2 ) goto 101 !frtati 2021/03/06
        call ALLOCATE_GGBANK
        call INIT_GGBANK
 101    continue

*-----------------------------------------------------------------------
* read each tally files
*-----------------------------------------------------------------------
        do m = 1, itnm

        if ( itl4flg.eq.1 .or. itl4flg.eq.2 ) then
          if ( nobch.eq.maxbch ) return
          nbcht = nobch
          if ( npe.gt.1 ) nbcht = nobch/(npe-1)
          if ( rescb2(m,nbcht+1).eq.0.d0 ) cycle
          resc2(m) = rescb2(m,nbcht+1)
          resc3(m) = rescb3(m,nbcht+1)
        end if

cOBINATA(2012.9.21): only head axis
        do 800 iax = 1, 1

*-----------------------------------------------------------------------
*         [t-track]
*-----------------------------------------------------------------------
          if ( ital(m) .eq. 1 ) then

            call read_ttrack(m,iax,ierr)
            if( ierr .ne. 0 ) then
              goto 900
            end if

*-----------------------------------------------------------------------
*         [t-adjoint]
*-----------------------------------------------------------------------

          else if ( ital(m) .eq. 19 ) then

            call read_tadjnt(m,iax,ierr)
            if( ierr .ne. 0 ) then
              goto 900
            end if

*-----------------------------------------------------------------------
*         [t-cross]
*-----------------------------------------------------------------------
          else if ( ital(m) .eq. 2 ) then

            call read_tcross(m,iax,ierr)
            if( ierr .ne. 0 ) then
              goto 900
            end if

*-----------------------------------------------------------------------
*         [t-yield]
*-----------------------------------------------------------------------
          else if ( ital(m) .eq. 3 ) then

            call read_tyield(m,iax,ierr)
            if( ierr .ne. 0 ) then
              goto 900
            end if

*-----------------------------------------------------------------------
*         [t-heat]
*-----------------------------------------------------------------------
          else if ( ital(m) .eq. 4 ) then

            call read_theat(m,iax,ierr)
            if( ierr .ne. 0 ) then
              goto 900
            end if

*-----------------------------------------------------------------------
*         [t-star]
*-----------------------------------------------------------------------
          else if ( ital(m) .eq. 5 ) then

            call read_tstar(m,iax,ierr)
            if( ierr .ne. 0 ) then
              goto 900
            end if

*-----------------------------------------------------------------------
*         [t-time]
*-----------------------------------------------------------------------
          else if ( ital(m) .eq. 6 ) then

            call read_ttime(m,iax,ierr)
            if( ierr .ne. 0 ) then
              goto 900
            end if

*-----------------------------------------------------------------------
*         [t-dpa]
*-----------------------------------------------------------------------
          else if ( ital(m) .eq. 7 ) then

            call read_tdpa(m,iax,ierr)
            if( ierr .ne. 0 ) then
              goto 900
            end if

*-----------------------------------------------------------------------
*         [t-product]
*-----------------------------------------------------------------------
          else if ( ital(m) .eq. 8 ) then

            call read_tproduct(m,iax,ierr)
            if( ierr .ne. 0 ) then
              goto 900
            end if

*-----------------------------------------------------------------------
*         [t-let]
*-----------------------------------------------------------------------
          else if ( ital(m) .eq. 12 ) then

            call read_tlet(m,iax,ierr)
            if( ierr .ne. 0 ) then
              goto 900
            end if

*-----------------------------------------------------------------------
*         [t-deposit]
*-----------------------------------------------------------------------
          else if ( ital(m) .eq. 13 ) then

            call read_tdeposit(m,iax,ierr)
            if( ierr .ne. 0 ) then
              goto 900
            end if

*-----------------------------------------------------------------------
*         [t-deposit2]
*-----------------------------------------------------------------------
          else if ( ital(m) .eq. 14 ) then

            call read_tdeposit2(m,iax,ierr)
            if( ierr .ne. 0 ) then
              goto 900
            end if

*-----------------------------------------------------------------------
*         [t-sed]
*-----------------------------------------------------------------------
          else if ( ital(m) .eq. 15 ) then

            call read_tsed(m,iax,ierr)
            if( ierr .ne. 0 ) then
              goto 900
            end if

*-----------------------------------------------------------------------
*         [t-point]
*-----------------------------------------------------------------------
          else if ( ital(m) .eq. 17 ) then

            call read_tpoint(m,iax,ierr)
            if( ierr .ne. 0 ) then
              goto 900
            end if

*-----------------------------------------------------------------------
*         [t-wwg]
*-----------------------------------------------------------------------
          else if ( ital(m) .eq. 18 ) then

            call read_twwg(m,iax,ierr)
            if( ierr .ne. 0 ) then
              goto 900
            end if

*-----------------------------------------------------------------------
*         [t-volume]
*-----------------------------------------------------------------------
          else if ( ital(m) .eq. 21 ) then

            call read_tvlm(m,iax,ierr)
            if( ierr .ne. 0 ) then
              goto 900
            end if

*-----------------------------------------------------------------------
          end if

 800    continue
        end do

        call trans_trRES2tr0()

 900    continue
        if ( itl4flg.eq.2 ) goto 102 !frtati 2021/03/06
        call DEALLOCATE_GGBANK
 102    continue

      end subroutine


************************************************************************
*                                                                      *
      subroutine trans_trRES2tr0()
*                                                                      *
*                                                                      *
************************************************************************

        use TALMOD
        use RESTALMOD

        implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

        include 'param.inc'

*-----------------------------------------------------------------------

        common /mpi00/ npe, me
        common /tall00/ itnm, ital(itlmax), itals(itlmax), italm(itlmax)

*-----------------------------------------------------------------------

        if ( npe .le. 1 ) then
          do m = 1,itnm
            tr0(italhead(m):italhead(m)+lrestalm(m)-1)
     &        = trRES(irestalm(m):irestalm(m)+lrestalm(m)-1)
            trRES(irestalm(m):irestalm(m)+lrestalm(m)-1) = 0.d0
          enddo

! sumover
          do m = 1,itnm
            if(mrestalm_sum(m) > 0) then
              tr0_sum(italhead_sum(m,1):
     &              italhead_sum(m,1)+mrestalm_sum(m)-1)
     &        = trRES_sum(irestalm_sum(m,1):
     &                    irestalm_sum(m,1)+mrestalm_sum(m)-1)
              trRES_sum(irestalm_sum(m,1):
     &                irestalm_sum(m,1)+mrestalm_sum(m)-1) = 0.d0
            endif
          enddo

        end if

      end subroutine

************************************************************************
*                                                                      *
      subroutine read_dmpinfo(ierr)
*                                                                      *
*       programmed by T.Furuta on 2015/05/15                           *
*                                                                      *
************************************************************************
      implicit real*8 (a-h,o-z)
      include 'param.inc'
      include 'err.inc'

      integer,intent(out) :: ierr

      integer i
      integer maxbchdmp,maxcasdmp
      real(8) rsouindmp,rsouinbch
      common /dmpinfo/ rsouinbch,maxbchdmp,maxcasdmp
      integer isorf,lsfile
      character sfile*100
      common /isorsf/ isorf(isrc),lsfile(isrc), sfile(isrc)

      integer jsn,jsi,idsi(0:9),ill(0:9),ilf(0:9)
      character dsin(0:9)*200
      integer iristdev,irmaxcas
      real(8) resc2,resc3,rrijklst

      jsn=0
      jsi=150
      idsi(:)=0
      ill(:)=1
      ilf(:)=10000000
      ierr=0

      do i=lsfile(1),1,-1
       if(sfile(1)(i:i).eq.'.')exit
      enddo
      if(i.eq.0)i=lsfile(1)+1

      dsin(jsn)(1:i-5) = sfile(1)(1:i-5)
      dsin(jsn)(i-4:lsfile(1)-4) = sfile(1)(i:lsfile(1))
      idsi(jsn) = lsfile(1)-4

      open(jsi,
     &     file=dsin(jsn)(1:idsi(jsn)),status='old',action='read',
     &     form='formatted',iostat=ierr)

      if(ierr.eq.0) then
       call rrestart(jsn,jsi,dsin,idsi,ill,ilf,
     &      iristdev,resc2,resc3,irmaxcas,rrijklst,ierr)
       close(jsi)
      endif

      if(ierr.eq.0)then
       maxcasdmp=irmaxcas
       if(iristdev.eq.1)then
        maxbchdmp=int(resc3)
        rsouindmp=resc2*maxcasdmp
       elseif(iristdev.eq.2)then
        maxbchdmp=int(resc3/maxcasdmp)
        rsouindmp=resc2
       else
        write(ErrCha,'("ERROR in read_dmpinfo: irstdev =",i5)')iristdev
        ErrID = 'L:1102/R:read_dmpinfo/F:restart.f' !E83_002_001
        call ErrWrite(ErrID,ErrCha)
        ierr=1
       endif
       rsouinbch=rsouindmp/dble(maxbchdmp)
      endif

      return
      end subroutine

************************************************************************
*                                                                      *
      subroutine setdmpinfo
*                                                                      *
*       programmed by T.Furuta on 2015/05/15                           *
*                                                                      *
************************************************************************
!$    use omp_lib

      implicit real*8 (a-h,o-z)

      include 'param.inc'

*-----------------------------------------------------------------------

      common /cparm/ maxbch,maxcas

*-----------------------------------------------------------------------

      real(8) dmpmulti
      integer idmpmode,ibchjmp,idmpjmp
      common /stat2/ dmpmulti,idmpmode,ibchjmp,idmpjmp(2)
      real(8) rsouinbch
      integer maxbchdmp,maxcasdmp
      common /dmpinfo/ rsouinbch,maxbchdmp,maxcasdmp

*-----------------------------------------------------------------------
      if(idmpmode.gt.0)then
       maxcas=maxcasdmp
       maxbch=maxbchdmp
      endif

      return
      end subroutine

