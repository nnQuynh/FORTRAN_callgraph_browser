!***********************************************************************
      module sangelmod
!                                                                      *
!     create by T.Miura on 2017/09/30                                  *
!     input data for special angel parameters of tally section         *
!                                                                      *
!***********************************************************************
      implicit none

      integer                            :: itsanf    ! sangel read flag
      integer                       ,save:: itsanm=0  ! total of sangels

      integer           ,allocatable,save:: itsans(:) ! Num of sangel
      integer           ,allocatable,save:: ltsans(:) ! location itsang

      character(len=200),allocatable,save:: itsang(:) ! sangel param
      integer           ,allocatable,save:: ltsang(:) ! itsang's length

      contains
!------------------------------------------------------------------------
        subroutine ALLOCATE_SANGEL
!........................................................................
        include 'param.inc'

        integer :: itnm
        integer :: ital
        integer :: itals
        integer :: italm
        common /tall00/ itnm, ital(itlmax), itals(itlmax), italm(itlmax)
!------------------------------------------------------------------------

        if ( itnm == 0 ) return
           allocate( itsans(itnm) )
           allocate( ltsans(itnm) )
           itsans(1:itnm) = 0
           ltsans(1:itnm) = 0

        if ( itsanm == 0 ) return
           allocate( itsang(itsanm) )
           allocate( ltsang(itsanm) )
           itsang(1:itsanm) = ' '
           ltsang(1:itsanm) = 0

        return
        end subroutine ALLOCATE_SANGEL

!------------------------------------------------------------------------
        subroutine DEALLOCATE_SANGEL
!........................................................................
        include 'param.inc'

        integer :: itnm
        integer :: ital
        integer :: itals
        integer :: italm
        common /tall00/ itnm, ital(itlmax), itals(itlmax), italm(itlmax)
!------------------------------------------------------------------------

        if ( itnm == 0 ) return
        if ( itsanm == 0 ) return

           deallocate ( itsans )
           deallocate ( ltsans )
           deallocate ( itsang )
           deallocate ( ltsang )

        return
        end subroutine DEALLOCATE_SANGEL

!------------------------------------------------------------------------
      end module sangelmod
