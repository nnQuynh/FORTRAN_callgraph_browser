!***********************************************************************
      module sumtallymod
!                                                                      *
!     The update to multiple tally inputs (T.Miura 2015/07/31)         *
!                                                                      *
!***********************************************************************
      implicit double precision (a-h,o-z)

      integer, allocatable :: nsumtalRead(:)  ! number of sumtally data section

      integer, allocatable :: isumtally(:)   ! sumtally option

      integer, allocatable :: nfile(:)       ! number of tally files

      character(len=200), allocatable ::  tallyfname(:,:) ! tally file name(s)
      integer           , allocatable :: ltallyfname(:,:) ! Number of characters

      double precision  , allocatable :: weightRate(:,:) ! weighting rate
      double precision  , allocatable :: sumWR(:) ! Sum of the weight rates

      character(len=200), allocatable ::  sfile(:) ! file name of output
      integer           , allocatable :: lsfile(:) ! Number of characters

      double precision  , allocatable :: sumfactor(:) ! normalization factor

      character(len=200), allocatable ::  sumang(:)  ! angel parameter
      integer           , allocatable :: lsumang(:)

!------------------------------------------------------------------------

      double precision, allocatable,save:: trSUMTAL(:)

      double precision,save,allocatable :: resc2SUMTAL(:)
      double precision,save,allocatable :: resc3SUMTAL(:)

      double precision, allocatable,save,target:: trSUMTAL_SUM(:)

      contains

!------------------------------------------------------------------------
       subroutine ALLOCATE_SUMTAL
!
!      called form sumtally subroutine
!........................................................................
       use RESTALMOD, only : lrestalm,mrestalm_sum
       include 'param.inc'
       common /tall00/ itnm, ital(itlmax), itals(itlmax), italm(itlmax)

       integer imax,m,imax_sum
!------------------------------------------------------------------------

       imax=1
       do m=1, itnm
        imax = imax + lrestalm(m)
       enddo

       imax_sum = 1
       imax_sum = imax_sum + sum(mrestalm_sum(:))

       allocate( trSUMTAL(imax) )
       trSUMTAL(1:imax) = 0.0d0

       resc2SUMTAL(:) = 0.0d0
       resc3SUMTAL(:) = 0.0d0

       allocate( trSUMTAL_SUM(imax_sum) )
       trSUMTAL_SUM(1:imax_sum) = 0.0d0

       end subroutine ALLOCATE_SUMTAL
!------------------------------------------------------------------------


!------------------------------------------------------------------------
       subroutine COPY_SUMTAL
!
!      called form sumtally subroutine
!........................................................................
       use TALMOD,   only : tr0, tr0_sum

       tr0(:)   = trSUMTAL(:)
       tr0_sum(:)   = trSUMTAL_SUM(:)

       end subroutine COPY_SUMTAL
!------------------------------------------------------------------------


!------------------------------------------------------------------------
       subroutine DEALLOCATE_SUMTAL
!
!      called form sumtally subroutine
!........................................................................
       deallocate( trSUMTAL )

       deallocate( trSUMTAL_SUM )

       deallocate( nsumtalRead )  ! set_sumtal subroutine
       deallocate( isumtally   )
       deallocate( nfile       )
       deallocate(  tallyfname )
       deallocate( ltallyfname )
       deallocate( weightRate  )
       deallocate( sumWR       )
       deallocate( sumfactor   )

       deallocate(  sfile      )
       deallocate( lsfile      )
       deallocate(  sumang     )
       deallocate( lsumang     )


       end subroutine DEALLOCATE_SUMTAL
!------------------------------------------------------------------------


!***********************************************************************
      end module sumtallymod
!***********************************************************************
