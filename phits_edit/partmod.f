************************************************************************
*
      module partmod
*
*     created by N.Furutachi on 2021/10/05
*
************************************************************************
      implicit none

      integer,save :: itptfirst = 0
      integer,save :: itmxpt = 1, itmxptt = 1
      integer,save :: itmxgp = 1, itmxgpt = 1
      integer,allocatable,save :: itmxang(:), itpan(:)
      integer,allocatable,save :: itpat(:,:,:), jtpat(:,:,:,:)
      integer,allocatable :: temp_itmxang(:), temp_itpan(:)
      integer,allocatable :: temp_itpat(:,:,:), temp_jtpat(:,:,:,:)
      integer,allocatable :: ipan(:), ipat(:,:,:)
      character*3,allocatable,save :: chl(:)
      character*5,allocatable,save :: chm(:)
      character*1,allocatable,save :: chb(:)
      character*8,allocatable,save :: group(:)
      integer,allocatable,save :: iznmmx(:), iznmstore(:,:)
      logical,allocatable,save :: firstwarning(:)
      integer,allocatable,save :: temp_iznmmx(:), temp_iznmstore(:,:)
      integer,allocatable,save :: iznmturn(:,:)

      contains

************************************************************************
        subroutine ALLOCATE_PART_TEMP
************************************************************************
        include 'param.inc'

        integer :: iflag_realloc
        integer :: itmxptold,itmxgpold,i1,i2

        if ( itptfirst.eq.0 ) then
          if ( itmxptt.gt.0 ) then
            itmxpt = itmxptt
          else if ( itmxptt.eq.0 ) then
            itmxpt = 1
          end if
          itmxgp = itmxgpt
          allocate( itmxang(itlmax), itpan(itlmax) )
          allocate( itpat(itlmax,itmxpt,3) )       ! kitamura22/03/31
          allocate( jtpat(itlmax,itmxpt,itmxgp,2) )
          itptfirst = 1
          return
        end if

        iflag_realloc = 0
        itmxptold = itmxpt
        itmxgpold = itmxgp
        if ( itmxptt.gt.itmxpt ) then
          itmxpt = itmxptt
          iflag_realloc = 1
        end if
        if ( itmxgpt.gt.itmxgp ) then
          itmxgp = itmxgpt
          iflag_realloc = 1
        end if

        if ( iflag_realloc.eq.1 ) then
          allocate( temp_itpat(itlmax,itmxpt,3) )  ! kitamura22/03/31
          allocate( temp_jtpat(itlmax,itmxpt,itmxgp,2) )
          temp_itpat(:,:,:)=0.0
          temp_jtpat(:,:,:,:)=0.0
          do i1=1,itmxptold
           temp_itpat(:,i1,:) = itpat(:,i1,:)
           do i2=1,itmxgpold
            temp_jtpat(:,i1,i2,:) = jtpat(:,i1,i2,:)
           enddo
          enddo
          deallocate( itpat )
          deallocate( jtpat )
          do i1=1,100
          enddo
          allocate( itpat(itlmax,itmxpt,3) )       ! kitamura22/03/31
          allocate( jtpat(itlmax,itmxpt,itmxgp,2) )
          itpat(:,:,:) = temp_itpat(:,:,:)
          jtpat(:,:,:,:) = temp_jtpat(:,:,:,:)
          deallocate( temp_itpat )
          deallocate( temp_jtpat )
        end if

        return
        end subroutine ALLOCATE_PART_TEMP

************************************************************************
        subroutine ALLOCATE_PART
************************************************************************
        include 'param.inc'

        integer :: i
        integer :: itnm, ital, itals, italm, itnzn, itndm ! frtati 2022/02/28
        common /tall00/ itnm, ital(itlmax), itals(itlmax), italm(itlmax)
        common /tall83/ itnzn(itlmax), itndm(itlmax) ! frtati 2022/02/28
        integer :: npe, me
        common /mpi00/ npe, me

        if ( itptfirst.eq.0 ) then
          call ALLOCATE_PART_TEMP
        end if

        allocate( temp_itmxang(itnm), temp_itpan(itnm) )
        allocate( temp_itpat(itnm,itmxpt,3) )      ! kitamura22/03/31
        allocate( temp_jtpat(itnm,itmxpt,itmxgp,2) )
        do i = 1, itnm
          temp_itmxang(i) = itmxang(i)
          temp_itpan(i) = itpan(i)
          temp_itpat(i,:,:) = itpat(i,:,:)
          temp_jtpat(i,:,:,:) = jtpat(i,:,:,:)
        end do
        deallocate( itmxang, itpan )
        deallocate( itpat )
        deallocate( jtpat )
        allocate( itmxang(itnm), itpan(itnm) )
        allocate( itpat(itnm,itmxpt,3) )           ! kitamura22/03/31
        allocate( jtpat(itnm,itmxpt,itmxgp,2) )
        itmxang(:) = temp_itmxang(:)
        itpan(:) = temp_itpan(:)
        itpat(:,:,:) = temp_itpat(:,:,:)
        jtpat(:,:,:,:) = temp_jtpat(:,:,:,:)
        deallocate( temp_itmxang, temp_itpan )
        deallocate( temp_itpat )
        deallocate( temp_jtpat )

        allocate( ipan(itnm) )
        allocate( ipat(itnm,itmxpt,3) )            ! kitamura22/03/31

        ipat(:,:,:)=0  ! T.Sato 2022/09/06 initialization is necessary

        allocate( chl(itmxpt), chm(itmxpt), chb(itmxpt) )
        allocate( group(itmxpt) )

        call SET_ANGELLINES

        allocate(iznmmx(itnm))
        allocate(iznmstore(maxval(itnzn),itnm))
        allocate(firstwarning(itnm))
        iznmmx(:) = 1
        iznmstore(:,:) = 0
        firstwarning(:) = .true.
        if( npe.gt.1 ) then
          allocate(temp_iznmstore(maxval(itnzn),itnm))
          allocate(iznmturn(maxval(itnzn),itnm))
          temp_iznmstore(:,:) = 0
          iznmturn(:,:) = 1
        end if

        return
        end subroutine ALLOCATE_PART

************************************************************************
        subroutine DEALLOCATE_PART
************************************************************************
        integer :: npe, me
        common /mpi00/ npe, me

        deallocate( itmxang )
        deallocate( itpan )
        deallocate( itpat )
        deallocate( jtpat )
        deallocate( ipan )
        deallocate( ipat )
        deallocate( chl, chm, chb )
        deallocate( group )
        deallocate( iznmmx, iznmstore, firstwarning ) ! frtati 2022/02/28
        if( npe.gt.1 ) then ! frtati 2022/05/02
          deallocate( temp_iznmstore, iznmturn )
        end if

        return
        end subroutine DEALLOCATE_PART

************************************************************************
        subroutine SET_ANGELLINES
************************************************************************

        integer :: i, ilim, datamax
        parameter ( datamax = 36 )

        character*3 :: chld(datamax)
        character*5 :: chmd(datamax)
        character*6 :: chc
        character*20 :: chi

        chld(1:6)  = (/'l  ','dr ','ub ','mg ','qrr','vbb'/)
        chld(7:12) = (/'lr ','db ','ug ','mrr','qbb','v  '/)
        chld(13:18)= (/'lb ','dg ','urr','mbb','q  ','vr '/)
        chld(19:24)= (/'lg ','drr','ubb','m  ','qr ','vb '/)
        chld(25:30)= (/'lrr','dbb','u  ','mr ','qb ','vg '/)
        chld(31:36)= (/'lbb','d  ','ur ','mb ','qg ','vrr'/)
        chmd(1:6)  = (/'l3   ','d4r  ','u5b  ','m6g  ','q7rr ','v8bb '/)
        chmd(7:12) = (/'l9r  ','d10b ','u11g ','m12rr','q13bb','v14  '/)
        chmd(13:18)= (/'l3b  ','d4g  ','u5rr ','m6bb ','q7   ','v8r  '/)
        chmd(19:24)= (/'l9g  ','d10rr','u11bb','m12  ','q13r ','v14b '/)
        chmd(25:30)= (/'l3rr ','d4bb ','u5   ','m6r  ','q7b  ','v8g  '/)
        chmd(31:36)= (/'l9bb ','d10  ','u11r ','m12b ','q13g ','v14rr'/)

        chb(:) = ' '
        chl(:) = 'l  '
        chm(:) = 'l3  '
        if ( itmxpt.le.datamax ) then
          ilim = itmxpt
        else
          ilim = datamax
        end if
        do i = 1, ilim
          chl(i) = chld(i)
          chm(i) = chmd(i)
        end do

        if ( itmxpt.lt.10 ) then
          do i = 1, itmxpt
            chc = '-group'
            write(chi,*) i
            group(i) = 'p' // trim(adjustl(chi)) // chc
          end do
        else
          do i = 1, 9
            chc = '-group'
            write(chi,*) i
            group(i) = 'p' // trim(adjustl(chi)) // chc
          end do
          do i = 10, itmxpt
            chc = 'group '
            write(chi,*) i
            group(i) = 'p' // trim(adjustl(chi)) // trim(chc)
          end do
        end if

        return
        end subroutine SET_ANGELLINES

************************************************************************
      end module partmod
