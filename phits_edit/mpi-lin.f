************************************************************************
*                                                                      *
      subroutine paraint(mn)
*                                                                      *
*        Initialization of MPI                                         *
*        Open input filename file                                      *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'mpif.h'
      include 'err.inc'

      common /mpi00/ npe, me
      common /mpi03/ filhd, nfihd
      common /mpi04/ stim1, stim2
      character filhd*100

*-----------------------------------------------------------------------

      common /inout/  in, io

*-----------------------------------------------------------------------

      logical   exex
      character filnm*100

      character user*20

*-----------------------------------------------------------------------
*        MPI initialization
*-----------------------------------------------------------------------

            call mpi_init(ierr)

*-----------------------------------------------------------------------

               if( ierr .ne. mpi_success ) then

                 write(ErrCha,*) ' *** Something wrong in MPI init. '
                 ErrID = 'L:44/R:paraint/F:mpi-lin.f' !E84_001_001
                 call ErrWrite(ErrID,ErrCha)

                  mn = -1
                  return

               end if

               mn = 0

*-----------------------------------------------------------------------

            call mpi_comm_size(mpi_comm_world, npe, ierr)
            call mpi_comm_rank(mpi_comm_world, me, ierr)

*-----------------------------------------------------------------------
*     open input file ( unit = 5 )
*     file name is 'phits.in', which is fixed in this version.
*-----------------------------------------------------------------------

         filnm = 'phits.in'
         nfiln =  8

         inquire( file = filnm, exist = exex )

            if( exex .eqv. .false. ) then

               if( me .eq. 0 ) then

               write(ErrCha,*) 'Input File Name Error. File not exist'//
     &                       ' ->> '//filnm(1:nfiln)
                 ErrID = 'L:75/R:paraint/F:mpi-lin.f' !E84_002_001
                 call ErrWrite(ErrID,ErrCha)

               end if

               mn = -1
               return

            end if

         open(in, file = filnm, status = 'old' )

*-----------------------------------------------------------------------
*     get user name and send it to all PE
*-----------------------------------------------------------------------

         if( me .eq. 0 ) then

               call getenv("LOGNAME",user)

               do k = 20, 1, -1

                  if( user(k:k) .ne. ' ' ) goto 500

               end do

                  write(*,*) '*** User name cannot be obtained'

                  mn = -1
                  return

  500          ku = k

                  filhd = '/wk/'//user(1:ku)//'/'
                  nfihd = ku + 5

         end if

                  call mpi_bcast(filhd,100,mpi_character,0,
     &                           mpi_comm_world,ierr)

                  call mpi_bcast(nfihd,1,mpi_integer,0,
     &                           mpi_comm_world,ierr)

*-----------------------------------------------------------------------
*     starting time
*-----------------------------------------------------------------------

            stim1 = paratim()

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      function paratim()
*                                                                      *
*        get elasptime                                                 *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'mpif.h'

*-----------------------------------------------------------------------
*     get time
*-----------------------------------------------------------------------

            paratim = mpi_wtime()

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine parafin
*                                                                      *
*        Finalization of MPI                                           *
*                                                                      *
************************************************************************
      use moddas_bends, only: moddas_bends_deallocate
      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'mpif.h'

      common /mpi00/ npe, me

cFURUTA20220106---------------------------------------------------------
      if(me.eq.0.and.npe.gt.1)call moddas_bends_deallocate
*-----------------------------------------------------------------------
*     finalization
*-----------------------------------------------------------------------

            call mpi_finalize(ierr)

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine paraiset(icc)
*                                                                      *
*        initial set for irskip, maxcas and maxbch                     *
*        modified by K.Niita on 2003/10/14                             *
*                                                                      *
************************************************************************
      use moddas_bends, only: moddas_bends_allocate
      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'mpif.h'
      include 'param.inc'
      include 'err.inc'

      common /mpi00/ npe, me
      common /mpi02/ nmbch0

*-----------------------------------------------------------------------

      common /inout/  in,io
      common /paraj/ mstz(300), parz(300)

*-----------------------------------------------------------------------

      common /iradkk/ randkk,irskip
      common /cparm/  maxbch,maxcas
      common /sumbnd/ nbnd(nbchmax) !FURUTA20220106

*-----------------------------------------------------------------------

      common /stat / istdev, irestart, ireschk
      common /res01/ istdevres,maxcasres,rijklstres,irdrf

*-----------------------------------------------------------------------

      if( me .eq. 0 .and. npe .gt. 1 ) then

         if( irskip .gt. 0 ) then

               write(io,'(/''*** Error in Parallel ***''/
     &                     ''    irskip should be =< 0'',
     &                     ''  : irskip ='',i6/)') irskip

               icc = 1

               goto 100

         end if

         if( istdev .eq. 1 .and. irestart .eq. 1 ) then

               maxca0 = maxcasres

               nmbch0 = max( 1, maxbch / ( npe - 1 ) )
               maxbc0 = nmbch0 * ( npe - 1 )

               if(maxbc0 .ne. maxbch) then
                write(io,'(/''*** Warning in Restart Parallel ***''/
     &          ''    maxbch should be a multiple of ( npe - 1 ):''
     &          , i4/''    maxbch was changed from'', i5,'' to'',i5/
     &          ''    while maxcas ='',i8,'' was read from resfile'')')
     &          npe-1,maxbch,maxbc0,maxca0
               endif

         elseif( mod( maxbch, npe-1 ) .eq. 0 ) then

               nmbch0 = maxbch / ( npe - 1 )
               maxca0 = maxcas
               maxbc0 = maxbch

         else

               nmbch0 = max( 1, maxbch / ( npe - 1 ) )

               rmtote = dble( nmbch0 * ( npe - 1 ) ) * dble( maxcas )
               rmdiff = dble( maxcas ) * dble( maxbch ) - rmtote
               nmcadd = nint( rmdiff / ( nmbch0 * ( npe - 1 ) ) )

               maxca0 = max( 1, maxcas + nmcadd )
               maxbc0 = nmbch0 * ( npe - 1 )

               write(io,'(/''*** Warning in Parallel ***''/
     &         ''    maxbch should be a multiple of ( npe - 1 ).''/
     &         ''    maxbch ='',i8,''  : npe-1 ='',i4)')
     &               maxbch, npe-1

               write(io,'(''    We have changed maxbch and maxcas as:''/
     &         ''    maxbch ='',i8,''  : maxcas ='',i8)')
     &               maxbc0, maxca0

         end if

           if( nmbch0 .gt. nbchmax ) then !FURUTA2016/03/10 kvlmax->nbchmax



             ErrCha = ''
             ErrID = 'L:287/R:paraiset/F:mpi-lin.f' !E84_003_001
             call ErrWriteIO(ErrID,ErrCha,io)
             call ErrWrite(ErrID,ErrCha)

             write(*,*) '/*** Error in Parallel ***/'

             write(*,'(''npe set to ('',i2,'' -1 ) = '',i2)')
     &       npe , npe-1

             write(*,'(''maximum maxbch is,'',i10,''*(npe-1) = '',i10)')
     &       nbchmax,nbchmax*(npe-1)

             write(*,'(''Please change maxbch less than'', i10,
     &       '' in this case. ( current maxbch = '',i10,'')'')')
     &       nbchmax*(npe-1),maxbch

c ----------------------------------------------------------------------

               icc = 1

               goto 100

            end if


            do i = 1, nmbch0

               nbnd(i) = ( i - 1 ) * ( npe - 1 )

            end do


      end if

  100    continue

               call parabcsti(icc)

               if( icc .ne. 0 ) return

*-----------------------------------------------------------------------

      if( me .eq. 0 .and. npe .gt. 1 ) then

         do i = 1, npe - 1

               irski0 = irskip !FURUTA

               call parasi(irski0,1,i)
               call parasi(maxca0,1,i)
               call parasi(nmbch0,1,i)

         end do

               maxbch = maxbc0
               maxcas = maxca0

      end if

      if( me .gt. 0 .and. npe .gt. 1 ) then

               call parari(mstz(2),1,0)
               call parari(mstz(3),1,0)
               call parari(mstz(4),1,0)

      end if

cFURUTA20220106---------------------------------------------------------
      if( me .eq. 0 .and. npe .gt. 1 )then
       call moddas_bends_allocate(nmbch0*(npe-1))
      endif
*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine paratal(tr,br,nr)
*                                                                      *
*        send and receive the tally results                            *
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'mpif.h'

      common /mpi00/ npe, me
      common /mpi01/ iccp(20000)

      common /talsav/ iptall

*-----------------------------------------------------------------------

      dimension tr(1), br(1)
      dimension ioc(0:1)

      data ioc / 16, 17 /

*-----------------------------------------------------------------------

         if( me .gt. 0 ) then

                        call parasr(tr(1),nr,0)

         else

*-----------------------------------------------------------------------

            if( iptall .eq. 0 ) then

                  do j = 1, nr

                     tr(j) = 0.0

                  end do

               do i = 1, npe - 1

                  if( iccp(i) .eq. 0 ) then

                           call pararr(br(1),nr,i)

                        do j = 1, nr

                           tr(j) = tr(j) + br(j)

                        end do

                  end if

               end do

*-----------------------------------------------------------------------

            else

                  open(ioc(0),form='unformatted',status='scratch')
                  open(ioc(1),form='unformatted',status='scratch')

                           jc = 1
                           ic = mod(jc,2)
                           io = ioc(ic)

                        do j = 1, nr

                           write(io) 0.0d0

                        end do

               do i = 1, npe - 1

                  if( iccp(i) .eq. 0 ) then

                           call pararr(tr(1),nr,i)

                           is = io
                           jc = jc + 1
                           ic = mod(jc,2)
                           io = ioc(ic)
                           rewind is
                           rewind io

                        do j = 1, nr

                           read(is) trt
                           trt = trt + tr(j)
                           write(io) trt

                        end do

                  end if

               end do

                           rewind io

                        do j = 1, nr

                           read(io) tr(j)

                        end do

                  close( io )
                  close( is )

            end if

*-----------------------------------------------------------------------

         end if

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      subroutine paratal_sumover(br,m)
*                                                                      *
*        send and receive the tally results                            *
*                                                                      *
*                                                                      *
************************************************************************

      use TALMOD
!$      use TALMOD0

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

C for nonshared_tally option
      integer italsh
      common /talsh/ italsh

      dimension br(*)

      real(8),pointer :: tr_sum(:)

      iax = 1

C for nonshared_tally option
!$       if(italsh .eq. 0) then
!$          call GET_TR_HEAD_POINTER0_SUM(tr_sum,m,iax)
!$       else
            call GET_TR_HEAD_POINTER_SUM(tr_sum,m,iax)
C for nonshared_tally option
!$       end if

      nr = mtalsize_sum(m)
      call paratal_sumover_sub(tr_sum,br,nr)

      return
      end

************************************************************************
*                                                                      *
      subroutine paratal_sumover_sub(tr,br,nr)
*                                                                      *
*        send and receive the tally results                            *
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'mpif.h'

      common /mpi00/ npe, me
      common /mpi01/ iccp(20000)

      common /talsav/ iptall

*-----------------------------------------------------------------------

      dimension tr(1), br(1)
      dimension ioc(0:1)

      data ioc / 151, 152 /

*-----------------------------------------------------------------------

         if( me .gt. 0 ) then

                        call parasr(tr(1),nr,0)

         else

*-----------------------------------------------------------------------

            if( iptall .eq. 0 ) then

                  do j = 1, nr

                     tr(j) = 0.0

                  end do

               do i = 1, npe - 1

                  if( iccp(i) .eq. 0 ) then

                           call pararr(br(1),nr,i)

                        do j = 1, nr

                           tr(j) = tr(j) + br(j)

                        end do

                  end if

               end do

*-----------------------------------------------------------------------

            else

                  open(ioc(0),form='unformatted',status='scratch')
                  open(ioc(1),form='unformatted',status='scratch')

                           jc = 1
                           ic = mod(jc,2)
                           io = ioc(ic)

                        do j = 1, nr

                           write(io) 0.0d0

                        end do

               do i = 1, npe - 1

                  if( iccp(i) .eq. 0 ) then

                           call pararr(tr(1),nr,i)

                           is = io
                           jc = jc + 1
                           ic = mod(jc,2)
                           io = ioc(ic)
                           rewind is
                           rewind io

                        do j = 1, nr

                           read(is) trt
                           trt = trt + tr(j)
                           write(io) trt

                        end do

                  end if

               end do

                           rewind io

                        do j = 1, nr

                           read(io) tr(j)

                        end do

                  close( io )
                  close( is )

            end if

*-----------------------------------------------------------------------

         end if

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine parabcsti(icc)
*                                                                      *
*        broadcast icc                                                 *
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'mpif.h'

      common /mpi00/ npe, me

*-----------------------------------------------------------------------

         if( npe .gt. 1 ) then

               call mpi_bcast(icc,1,mpi_integer,
     &                        0,mpi_comm_world,ierr)

         end if

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine parasi(icc,nc,ip)
*                                                                      *
*        send integer icc to ip                                        *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'mpif.h'

      common /mpi00/ npe, me

*-----------------------------------------------------------------------

      dimension icc(1)

*-----------------------------------------------------------------------

           call mpi_send(icc(1),nc,mpi_integer,ip,me,
     &                   mpi_comm_world,ierr)

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      subroutine parari(icc,nc,ip)
*                                                                      *
*        receive integer icc from ip                                   *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'mpif.h'

      common /mpi00/ npe, me

*-----------------------------------------------------------------------

      integer istatus(mpi_status_size)

*-----------------------------------------------------------------------

      dimension icc(1)

*-----------------------------------------------------------------------

           call mpi_recv(icc(1),nc,mpi_integer,ip,ip,
     &                   mpi_comm_world,istatus,ierr)

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      subroutine parasr(rcc,nc,ip)
*                                                                      *
*        send double rcc to ip                                         *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'mpif.h'

      common /mpi00/ npe, me

*-----------------------------------------------------------------------

      dimension rcc(1)

*-----------------------------------------------------------------------

           call mpi_send(rcc(1),nc,mpi_double_precision,ip,me,
     &                   mpi_comm_world,ierr)

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      subroutine pararr(rcc,nc,ip)
*                                                                      *
*        receive double rcc from ip                                    *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'mpif.h'

      common /mpi00/ npe, me

*-----------------------------------------------------------------------

      integer istatus(mpi_status_size)

*-----------------------------------------------------------------------

      dimension rcc(1)

*-----------------------------------------------------------------------

           call mpi_recv(rcc(1),nc,mpi_double_precision,ip,ip,
     &                   mpi_comm_world,istatus,ierr)

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine paraiccp
*                                                                      *
*        receive integer icc from all IP                               *
*                                                                      *
************************************************************************


      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'mpif.h'

      common /mpi00/ npe, me
      common /mpi01/ iccp(20000)

*-----------------------------------------------------------------------

      integer istatus(mpi_status_size)

*-----------------------------------------------------------------------

         do i = 1, npe - 1

            if( iccp(i) .eq. 0 ) then

               call mpi_recv(iccp(i),1,mpi_integer,i,i,
     &                       mpi_comm_world,istatus,ierr)

            end if

         end do

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine parastop(idnum)
*                                                                      *
*        abnormal stop and send message to main PE                     *
*                                                                      *
*           idnum : stop code                                          *
*                                                                      *
************************************************************************


      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'mpif.h'
      include 'err.inc'

      common /mpi00/ npe, me

*-----------------------------------------------------------------------
! T.Sato 2015/03/09, to avoid forced stop due to parastop
      common /errorcom/ncascerr
!$OMP THREADPRIVATE(/errorcom/)

      integer*8 :: iransb64 ! S.H. xorshift (2020.2.6)
      common /randtp/ rijk,rans,ranb, iransb64
!$OMP THREADPRIVATE(/randtp/)

      if(idnum.eq.111.or.idnum.eq.115.or.idnum.eq.116.or.idnum.eq.119
     &  .or.idnum.eq.120)then ! error from JAM
       if(ncascerr.eq.0) then ! first time called this routine
        write(ErrCha,*) 'parastop is called by idnum =',idnum
        ErrID = 'L:890/R:parastop/F:mpi-lin.f' !E84_004_001
        call ErrWrite(ErrID,ErrCha)
        Select Case(idnum)
        case(111)
        write(*,*) 'Failure in jet production. Resampling is required'
        case(115)
        write(*,*) 'Unknown parton distribution. Resampling is required'
        case(116)
        write(*,*) 'Too low CM energy. Resampling is required'
        case(119)
        write(*,*) '(g,N)pi reaction skipped'
        case(120)
         write(*,*) 'ATIMA calculates stopping power of particle with
     &   charge <= 0 '
        end select
        write(*,*) 'me=',me
       endif
       ncascerr = -1
       return
      endif
! End of revision

         write(6,'('' *** Program is stopped, stop number ='',i6)')
     &                 idnum
         write(6,'('' Stop IP ='',i3)') me

*-----------------------------------------------------------------------

               icc = 1

               call parasi(icc,1,0)

               call mpi_finalize(ierr)

*-----------------------------------------------------------------------

      stop
      end

