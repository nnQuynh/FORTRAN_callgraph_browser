!$    use omp_lib  !OBINATA(2012.5.23)
      use RESTALMOD
C for noshared_trally option
!$    use TALMOD0
      use TALMOD
      use dedx_file
      use CHARVARMOD, only: close_rwtfile
      use sumtallymod
      use sangelmod
      use LAFDATAMOD, only: DEALLOCATE_laf !FURUTA20201127
      use partmod ! frtati 2021/10/05
! y.sakaki 2023/07/19, For user defined model -->
      use udm_Parameter
      use udm_Utility
      use udm_Manager

      implicit real*8 (a-h,o-z)

************************************************************************
*                                                                      *
*             _/_/_/_/                                                 *
*            _/      _/                  _/_/_/_/_/_/                  *
*           _/      _/  _/      _/   _/      _/      _/_/_/_/_/        *
*          _/_/_/_/    _/      _/   _/      _/      _/                 *
*         _/          _/_/_/_/_/   _/      _/       _/_/_/_/           *
*        _/          _/      _/   _/      _/              _/           *
*       _/          _/      _/   _/      _/      _/_/_/_/_/            *
*                                                                      *
*______________________________________________________________________*
*                                                                      *
*                              P H I T S                               *
*                                                                      *
*              Particle and Heavy Ion Transport code System            *
*                                                                      *
                      parameter ( Version = 3.350 )
*                                                                      *
                 parameter ( LastRevised = 20250328 )
*                                                                      *
************************************************************************

      include 'param.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      common /mpi00/ npe, me
      common /mpi01/ iccp(20000)

*-----------------------------------------------------------------------

      common /inout/  in,io
      common /qparm/  ielas,icasc,iqstep,lvlopt,igamma
      common /cparm/  maxbch,maxcas

      common /inpec/  ititl, ipara, ibody, iregn, llarr, itby, itar
      common /ccggg/  icgg

*-----------------------------------------------------------------------

      common /tcntl/  icntl, inucr

*-----------------------------------------------------------------------


      common /startf/ iday0,imon0,iyer0,ihor0,imin0,isec0
      common /verjam/ versn, lastr, iyeav, imonv, idayv

      common /emode/  emodem, ge1, ge2, iemode

      common /egsemi/ iegsemi, iegsout
C for nonshared_tally option
      integer italsh
      common /talsh/ italsh

*-----------------------------------------------------------------------


      data iccp /20000*0/

      common /stat / istdev, irestart, ireschk
      data istdev / 0 /
      data irestart   / 0 /

      common /restart/ resc2(itlmax), resc3(itlmax)

*-----------------------------------------------------------------------
      common /stat2/ dmpmulti,idmpmode,ibchjmp,idmpjmp(2)
      data dmpmulti /0.0d0/
      data idmpmode /0/
      data ibchjmp /0/
      data idmpjmp /0,0/
*-----------------------------------------------------------------------
cFURUTA20150714 TETRA
      integer nlat3,itetcl(10),itetvol,itetauto,itgchk
      common /tetf1/ nlat3,itetcl,itetvol,itetauto,itgchk
*-----------------------------------------------------------------------

      call INIT_TALMOD
      resc2(:) = 0.0
      resc3(:) = 0.0

      call init_sumtal
      call init_anatal


*-----------------------------------------------------------------------
*        io configuration
*-----------------------------------------------------------------------

*           io =  5      is input file
*           io =  6      is output file

*           io =  7      is reserved for xsdir file(7)
*           io =  8      is reserved for getxst
*
*           io =  9      is reserved for dumpall
*
*           io = 10      is reserved for pcut file(10)
*           io = 11      is reserved for nuclear reaction file(11)
*           io = 12      is reserved for ncut file(12)
*           io = 13      is reserved for gcut file(13)
*           io = 14      is reserved for photon data file(14)
*
*-----------------------------------------------------------------------
*           temporary use ( 15 - 26, 71 - 79, 80 )
*-----------------------------------------------------------------------

*           io = 15      is temporary used for read01 (isc)

*           io = 15      is temporary used for CG temp (iot)
*           io = 16      is temporary used for CG setup (n16)
*           io = 17      is temporary used for CG setup (n17)
*           io = 18      is temporary used for CG body (itby)
*           io = 19      is temporary used for CG array (itar)
*           io = 20      is temporary used for CG,GG info (iog)
*
*           io = 21      is temporary used for CGVIEW file(2) (iot)
*           io = 22      is temporary used for CGVIEW file(3) (ios)

*           io = 15      is temporary used for GG fill      (ioc)
*           io = 16      is temporary used for GG surface   (ioa)
*           io = 71      is temporary used for GG transform (iob)
*
*           io = 26      is temporary used for region/cell info.  (iod)
*           io = 22      is temporary used for setpar of region (ioe)
*
*-----------------------------------------------------------------------

*           io = 15      is temporary used for tallsave & multiplier
*           io = 16      is temporary used for paratal
*           io = 17      is temporary used for paratal

*-----------------------------------------------------------------------

*           io = 18      is temporary used for gshow
*           io = 19      is temporary used for 3dshow
*           io = 20      is temporary used for tally output.err

*           io = 23      is temporary used for material (iom1)
*           io = 24      is temporary used for material (iom2)
*           io = 25      is temporary used for material (iom3) in setpar

*           io = 28      is temporary used for cosmic-ray source
*           io = 29      is temporary used for natural abundance
*           io = 71      is temporary used for magnet field

*-----------------------------------------------------------------------
*           ( 18 - 26, 71 - 79 )
*-----------------------------------------------------------------------

*           io = 18-26, 71-79  is temporary used for ANGEL

*-----------------------------------------------------------------------

*           io = 80      is temporary used for read00

*-----------------------------------------------------------------------
*           ( 81 - 100 )
*-----------------------------------------------------------------------

*           io = 81-100  is reserved for dump files

*-----------------------------------------------------------------------
*           ( 151 - 200 )
*-----------------------------------------------------------------------

*           io = 151-200  is reserved for user defined tally (usrtally)

*-----------------------------------------------------------------------
*           ( 201 )
*-----------------------------------------------------------------------

*           io = 201  is used in subroutine tdchain (tallsm4.f)

*-----------------------------------------------------------------------
*           ( 202 )
*-----------------------------------------------------------------------

*           io = 202  is reserved for sangel parameter in tally (sangel.f)

*-----------------------------------------------------------------------
*           ( 201 - 299 )
*-----------------------------------------------------------------------

*           io = 201-299 is used for multiplier

*-----------------------------------------------------------------------
*           ( 27 - 65 )
*-----------------------------------------------------------------------

*           io = 27      is reserved for batch.out file
*           io = 28      is reserved for file(6)
*           io = 29      is used for ECHOI

*           io = 30-39   are reserved for input files
*           io = 40-54   are reserved for bank temporary files

*           io = 55 - 60 are reserved for JAM
*           io = 61 - 65 are reserved for ovly13 (file(11),....)

*-----------------------------------------------------------------------
*           ( 66 - 67 )
*-----------------------------------------------------------------------

*           io = 66      is reserved for getGDRxsec (temp)
*           io = 66      is reserved for voxel.bin (use for long time) ! T.Sato 2018/09/10
*           io = 67      is reserved for Angel file name (only Windows Openmp)

*-----------------------------------------------------------------------
*           ( 68 - 69 )
*-----------------------------------------------------------------------

*           io = 68      is reserved for node file in tetrainit (temp)
*           io = 69      is reserved for element file in tetrainit (temp)

*-----------------------------------------------------------------------
*           ( 301 -- 399 )
*-----------------------------------------------------------------------

*           io = 301--  are reserved for RIPL data read
*                     (gammod.f)

*-----------------------------------------------------------------------
*           ( 790 - 792 )
*-----------------------------------------------------------------------

*           io = 790-792 are reserved for magnetic field map ! T.Sato 2018/12/19

*-----------------------------------------------------------------------
*           ( 801 -- 899 )
*-----------------------------------------------------------------------

*           io = 801--  are reserved for dumpall or source with dump data
*                     (ovly15.f, ompdump.f)

*-----------------------------------------------------------------------
*           ( 901 -- 909 )
*-----------------------------------------------------------------------

*            io = 901       : Re Write File
*            io = 902 - 909 : infl File

*-----------------------------------------------------------------------
*           ( 2000 -- 2200 )
*-----------------------------------------------------------------------

*           io = 2000 - 2200 are reserved for [t-4Dtrack] (talls_4dtrack.f)

*-----------------------------------------------------------------------
*           ( 2201 -- 2300 )
*-----------------------------------------------------------------------

*           io = 2201 - 2300 are reserved for [use defined interaction]

*-----------------------------------------------------------------------
*           ( 3000 -- 3010 )
*-----------------------------------------------------------------------

*           io = 3000-3003 -  are reserved for atomic deexcitation data
*                     (elemdatamod.f)

*-----------------------------------------------------------------------
*           ( 4001 - 4999 )
*-----------------------------------------------------------------------

*           io = 4001-4999 is reserved for source files (isorf)
*                (changed from 100-150)

*-----------------------------------------------------------------------
*        input and out put file unit (temporary)
*-----------------------------------------------------------------------

               in = 5
               io = 6
               jo = 6

*-----------------------------------------------------------------------
*        error flag
*-----------------------------------------------------------------------

               ierr = 0
               icc  = 0

*-----------------------------------------------------------------------
*        parallel flag
*-----------------------------------------------------------------------

               npe = 0
               me  = 0
*-----------------------------------------------------------------------
*        Zenkaku Space flag
*-----------------------------------------------------------------------
               zspc_count       = 0
               ZSPC_Eflag       = 0
               ZenkakuSpaceLine = ''

*-----------------------------------------------------------------------
*        MPI initialization
*-----------------------------------------------------------------------

               call paraint(ierr)

               if( ierr .eq. -1 ) goto 999

*-----------------------------------------------------------------------
*        Version and Last Reveiced
*-----------------------------------------------------------------------

               versn = Version
               lastr = LastRevised

               iyeav = lastr / 10000
               imonv = ( lastr - iyeav * 10000 ) / 100
               idayv = lastr - iyeav * 10000 - imonv * 100

*-----------------------------------------------------------------------
*        get starting time and date
*-----------------------------------------------------------------------

               call date_a_time(iyer0,imon0,iday0,
     &                          ihor0,imin0,isec0)

*-----------------------------------------------------------------------
*        total CPU time
*-----------------------------------------------------------------------

               call cputime(1)

*-----------------------------------------------------------------------
*        read input data and initialization of control PE
*-----------------------------------------------------------------------

         if( me .eq. 0 ) then

               call read00(in,ivers,ierr)

                  if( ierr .ne. 0 ) goto 200
                  if ( icntl .eq. 16 ) goto 999

*-----------------------------------------------------------------------
*        set up geometry of tetrahedrons
*-----------------------------------------------------------------------
             if( nlat3 .gt. 0 ) call read_tetra(ierr)
                  if( ierr .ne. 0 ) goto 200
*-----------------------------------------------------------------------
*        read input [source] section check
*-----------------------------------------------------------------------
               call SourceCheck(ierr)
                  if( ierr .ne. 0 ) goto 200
*-----------------------------------------------------------------------
*        read input [paramater] section zspc mode check
*-----------------------------------------------------------------------
               if(ZSPC_Eflag == 1) call ZspcLog

*-----------------------------------------------------------------------
*        event dump source.
* FURUTA20150515
*-----------------------------------------------------------------------
               if(idmpmode.eq.1)then
                call read_dmpinfo(ierr)
                if(ierr.ne.0)goto 200
               endif
*-----------------------------------------------------------------------
*        read restart file.
* OBINATA(2012.6.4)
*-----------------------------------------------------------------------

               if( irestart .ne. 0 ) then
                 call read_resfiles(io,jo,ivers,ierr)
                 if( ierr .ne. 0 ) goto 200
               end if

               call setpar(io,jo,ivers,ierr)
                  if( ierr .ne. 0 ) goto 200


            if( iegsemi .ne. 0 ) then

               call setegs(io,jo,ierr)
                  if( ierr .ne. 0 ) goto 200

               call egs5init

            end if

               call setresval(io,jo,ivers,ierr)
                 if( ierr .ne. 0 ) goto 200
               if(idmpmode.eq.1)call setdmpinfo

               call setgg(io,jo,ierrg)
               call setcg(io,jo,ierrg)

            if( ierrg .eq. 0 ) then

               call setmd(io,jo,ierr)
                  if( ierr .ne. 0 ) goto 200

               call setpag(io,jo,ierr)
                  if( ierr .ne. 0 ) goto 200

               call set_anatal(io,jo,ierr)
                  if( ierr .ne. 0 ) goto 200

               if ( icntl .eq. 17 ) then
                  call anatal_set_tal(io,jo,ivers,ierr)
                  if( ierr .ne. 0 ) goto 200
               end if

               call settal(io,jo,ierr)
                  if( ierr .ne. 0 ) goto 200

               call setbnk(io,jo,ierr)
                  if( ierr .ne. 0 ) goto 200

            end if
!$           if(.true.) then
!$            if(italsh .eq. 0 ) then
!$             call ALLOCATE_TAL
!$             call INIT_TALMOD0       ! duplicate from TALMOD0
!$            else
!$             call ALLOCATE_TAL
!$            end if
!$           else
               call ALLOCATE_TAL
!$           endif

            call ALLOCATE_SANGEL
            call set_sangel

            call ALLOCATE_PART ! frtati 2021/10/05

! T.Sato 2016/06/24 Output icntl in console
      if(icntl.ne.0) then
       if(icntl.eq.1) then
        write(*,*) 'icntl = 1: Nuclear reaction mode'
       elseif(icntl.eq.5) then
        write(*,*) 'icntl = 5: No reaction, no ionization mode'
       elseif(icntl.eq.6) then
        write(*,*) 'icntl = 6: Source check mode'
       elseif(icntl.eq.7) then
        write(*,*) 'icntl = 7: Draw geometry by [t-gshow]'
       elseif(icntl.eq.8) then
        write(*,*) 'icntl = 8: Draw geometry by gshow option'
       elseif(icntl.eq.9) then
        write(*,*) 'icntl = 9: Draw geometry by [t-rshow]'
       elseif(icntl.eq.10) then
        write(*,*) 'icntl = 10: Draw geometry by rshow option'
       elseif(icntl.eq.11) then
        write(*,*) 'icntl = 11: Draw 3D geometry by [t-3dshow]'
       elseif(icntl.eq.12) then
        write(*,*) 'icntl = 12: Read data from dumpall'
       elseif(icntl.eq.13) then
        write(*,*) 'icntl = 13: Sumtally mode'
       elseif(icntl.eq.14) then
        write(*,*) 'icntl = 14: T-Volume mode'
       elseif(icntl.eq.15) then
        write(*,*) 'icntl = 15: T-WWBG mode'
       elseif(icntl.eq.17) then
        write(*,*) 'icntl = 17: Anatally mode'
       endif
      endif


                call set_sumtal(io,jo,ierr)
                if( ierr .ne. 0 ) goto 200


               call echoi(io,jo,ivers,ierr,ierrg)


                  if( ierr  .ne. 0 ) goto 200
                  if( ivers .eq. 0 ) goto 200

               call setrnd(io,jo,ierr)

                  if( ierr .ne. 0 ) goto 200

                  if( icntl .eq.  2 ) then
                     call cgview
                  else if( icntl .eq.  7 ) then
                     call gshows
                  else if( icntl .eq.  8 ) then
                     call gshowt
                  else if( icntl .eq.  9 ) then
                     call rshows
                  else if( icntl .eq. 10 ) then
                     call rshowt
                  else if( icntl .eq. 11 ) then
                     call dshows
                  end if

            if( icntl .eq. 0 .or. icntl .eq. 1 .or.
     &          icntl .eq. 5 .or. icntl .eq. 6 .or.
     &          icntl .eq. 12 .or.
     &          icntl .eq. 13 .or.
     &          icntl .eq. 14 .or.
     &          icntl .eq. 15 .or.
     &          icntl .eq. 17 ) goto 100

  200             icc = 1

         end if

*-----------------------------------------------------------------------
*        input error or end of simple job
*-----------------------------------------------------------------------

  100          continue

               call parabcsti(icc)

               if( icc .eq. 1 ) goto 900

*-----------------------------------------------------------------------
*     initialization for parallel
*-----------------------------------------------------------------------

      if( npe .ge. 1 ) then

*-----------------------------------------------------------------------
*        read input data in each PE
*-----------------------------------------------------------------------

         if( me .gt. 0 ) then

               call read00(in,ivers,ierr)

*-----------------------------------------------------------------------
*        set up geometry of tetrahedrons
*-----------------------------------------------------------------------
             if( nlat3 .gt. 0 ) call read_tetra(ierr)

                  if( ierr .ne. 0 ) goto 900
*-----------------------------------------------------------------------
*        read restart file.
* OBINATA(2012.6.4)
*-----------------------------------------------------------------------
            if( irestart .ne. 0 ) then
              call read_resfiles(io,jo,ivers,ierr)

              if( ierr .ne. 0 ) goto 900
            end if

         end if

*-----------------------------------------------------------------------
*        set new irskip, maxcas and maxbch for each PE
*-----------------------------------------------------------------------

               call paraiset(icc)

               if( icc .eq. 1 ) goto 900

*-----------------------------------------------------------------------
*        initialization of execution PE
*-----------------------------------------------------------------------

         if( me .gt. 0 ) then

               call setpar(io,jo,ivers,ierr)

            if( iegsemi .ne. 0 ) then
               call egs5init

            end if

               call setresval(io,jo,ivers,ierr)  !OBINATA(2012.6.18)
               call setgg(io,jo,ierr)
               call setcg(io,jo,ierr)
               call setmd(io,jo,ierr)
               call setpag(io,jo,ierr)
               call settal(io,jo,ierr)
               call setbnk(io,jo,ierr)
               call setrnd(io,jo,ierr)

!$             if(.true.) then
!$              if(italsh .eq. 0 ) then
!$               call ALLOCATE_TAL
!$               call INIT_TALMOD0       ! duplicate from TALMOD0
!$              else
!$               call ALLOCATE_TAL
!$              end if
!$             else
               call ALLOCATE_TAL

!$             endif


               call ALLOCATE_PART ! frtati 2021/10/05

         end if

      end if

      if ( me .eq. 0 )
     &  call ALLOCATE_RESTAL


*-----------------------------------------------------------------------
*     initialization for restart, when istdev is less than 0.
* OBINATA(2012.6.4)
*-----------------------------------------------------------------------
       call setrnd_restart()

      if ( irestart .ne. 0 ) then
        if ( me .eq. 0 )
     &    call read_talls(io,ierr)


        call parabcsti(ierr)

        if ( ierr .ne. 0 ) goto 890

        call paraiznm_restart ! frtati 2022/05/02

      end if

C -- set char parameter rwt file delete flag ---------------------------
      call close_rwtfile(0)

*-----------------------------------------------------------------------
*     initialization For user defined model
*-----------------------------------------------------------------------
      if(iudmodel .ge. 1) then
*       ----------------------------------------------------------------
        call udm_initialize_Utility
C       -------------
        do i=1,udm_int_num
          call user_defined_interaction(1,i) ! Fill Parameters
        enddo
        call user_defined_particle(1) ! Fill Parameters
C       -------------
C       Check whether the entered modules exist
        call udm_check(1)
        call udm_check(2)
C       ----
        do i=1,udm_int_num
          udm_logical=.false.
          call user_defined_interaction(2,i) ! check whether name match
          if(.not. udm_logical) then
            print*,
     &"In [ User defined interaction ], there is no module corresponding
     & to the entered Name -> ",trim(udm_int_name(i))
            stop
          endif
        enddo
C       ----
        udm_integer=0
        call user_defined_particle(2) ! Count number of existing modules (for particle)
        if(udm_integer .ne. udm_part_num) then
          print*,
     &   "In [ User defined particle ], there is no module corresponding
     & to the entered Name"
          stop
        endif
C       -------------
C       'Initialize' should be after 'Fill Parameters'.
        do i=1,udm_int_num
          call user_defined_interaction(3,i) ! Initialize
        enddo
        call user_defined_particle(3) ! Initialize
C       -------------
        do i=1,udm_int_num
          call user_defined_interaction(4,i) ! Print Comments
        enddo
        call user_defined_particle(4) ! Print Comments
*       ----------------------------------------------------------------
      endif

*-----------------------------------------------------------------------
*     ovly12, ovly13 or ovly14
*-----------------------------------------------------------------------

            if( icntl .eq. 0 .or. icntl .eq. 5 .or.
     &          icntl .eq. 6 .or.
     &          icntl .eq. 14 .or.
     &          icntl .eq. 15 ) then

               call ovly12

            else if( icntl .eq. 1 .and. inucr .eq. 20 ) then

               call ovly14

            else if( icntl .eq. 1 ) then

               call ovly13

            else if( icntl .eq. 12 ) then

               call ovly15

            else if( icntl .eq. 13 ) then    ! sum tally function

               call sumtally(io,jo,ivers,ierr)

            else if( icntl .eq. 17 ) then    ! anatally function

               call anatally(io,jo,ivers,ierr)

            end if

  890 continue

CCSE change for sangel parameter (2017.09.30) >>>>>
      if ( me .eq. 0 ) then
        call DEALLOCATE_RESTAL
        call DEALLOCATE_SANGEL
        call DEALLOCATE_PART ! frtati 2021/10/05
      end if
CCSE change for sangel parameter (2017.09.30) <<<<<
      call DEALLOCATE_TAL
      call dedx_file_deallocate(mxmat)
      call DEALLOCATE_laf !FURUTA20201127

*-----------------------------------------------------------------------
*     summary of cpu time
*-----------------------------------------------------------------------

  900 continue

               call sumcput(ierr,icc)

  999 continue
! T.Sato 2015/8/24, deallocate EGS memory only when EGS is used
            if( iegsemi .ne. 0 .and. icc .eq. 0 ) then

               CALL DEALLOCATES_EGS5 !<-
            endif

C -- set char parameter rwt file delete flag ---------------------------
               call close_rwtfile(1)

               call parafin

      stop
      end
