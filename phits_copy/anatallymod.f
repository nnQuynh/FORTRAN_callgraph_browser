!***********************************************************************
      module anatallymod
!                                                                      *
!     The update to multiple tally inputs (T.Miura 2015/07/31)         *
!                                                                      *
!***********************************************************************
      implicit double precision (a-h,o-z)

      integer, allocatable :: nanatalRead(:)  ! number of anatally data section

      integer, allocatable :: manatally(:)   ! anatally option
      integer, allocatable :: ianataldchain(:)   ! anatally for dchain

      character(len=200), allocatable ::  tallyfolder(:) ! folder name of output
      integer           , allocatable :: ltallyfolder(:) ! Number of characters
      character(len=200), allocatable ::  tallyfile(:) ! file name of output
      integer           , allocatable :: ltallyfile(:) ! Number of characters

      character(len=200), allocatable ::  cx_txt(:)  ! x-txt
      integer           , allocatable :: lcx_txt(:)
      character(len=200), allocatable ::  cy_txt(:)  ! y-txt
      integer           , allocatable :: lcy_txt(:)
      character(len=200), allocatable ::  cz_txt(:)  ! z-txt
      integer           , allocatable :: lcz_txt(:)

!------------------------------------------------------------------------

      double precision, allocatable,save:: trANATAL(:)
      integer         , allocatable,save:: lanatalm(:) ! S.H. 2020.6.17
!sumover
      integer         , allocatable,save:: lanatalm_sum(:)
      integer         , allocatable,save:: ianatalm_sum(:,:,:)
      integer         , allocatable,save:: manatalm_sum(:,:,:)

      double precision, allocatable,save, target:: trANATAL_SUM(:)

      double precision,save,allocatable :: resc2ANATAL(:)
      double precision,save,allocatable :: resc3ANATAL(:)

!------------------------------------------------------------------------
      double precision,allocatable:: anataldata(:,:,:,:)
      double precision,allocatable:: fg(:), fgaxs(:),fgaxs2(:),fgaxs3(:)
      double precision,allocatable:: delvol(:,:)

!------------------------------------------------------------------------

      contains

!------------------------------------------------------------------------
       subroutine allocate_anatal
!
!      called form anatally subroutine
!........................................................................
       use RESTALMOD, only : lrestalm, mrestalm_sum, irestalm_sum,
     &                       lrestalm_sum
       use sumtallymod, only: nfile
       include 'param.inc'
       common /tall00/ itnm, ital(itlmax), itals(itlmax), italm(itlmax)

       integer imax,m,imax_sum
!------------------------------------------------------------------------

       allocate( lanatalm(itnm) ) ! S.H. 2020.6.17
!sumover
       allocate( lanatalm_sum(itnm) )

       nf_sum = 0
       do m=1, itnm
         nf_sum = nf_sum + nfile(m)
       enddo

! mrestalm_sum(m)     : one tally length                 module RESTALMOD
! irestalm_sum(m,iax) : start address of one sumover     module RESTALMOD
! lrestalm_sum(m,iax) : one sumover length               module RESTALMOD

! lanatalm_sum(m)     :  start address of one tally at nfile(m)=1
! ianatalm_sum(m,iax,nfile(m)) :  start address of one sumover
! manatalm_sum(m,iax,nfile(m)) :  one sumover length
! nfile(m)            :  file numbers of same tally

        nfilemax = 0
       imax=1
       imax_sum=1
       do m=1, itnm
        lanatalm(m) = imax           ! start address
        lanatalm_sum(m) = imax_sum   ! start address
        imax=imax + lrestalm(m) * nfile(m)
        imax_sum=imax_sum + mrestalm_sum(m) * nfile(m)
        nfilemax = max(nfilemax,nfile(m))
       enddo

       allocate( ianatalm_sum(itnm,6,nfilemax) )
       allocate( manatalm_sum(itnm,6,nfilemax) )

       imax_sum_file=1
       do m=1, itnm
         do n=1, nfile(m)
           do iax=1,6
             manatalm_sum(m,iax,n)  = lrestalm_sum(m,iax)
             ianatalm_sum(m,iax,n)  = imax_sum_file
             imax_sum_file = imax_sum_file + lrestalm_sum(m,iax)
           enddo
         enddo
       enddo

       allocate( tranatal(imax) )
       tranatal(1:imax) = 0.0d0

       allocate( tranatal_sum(imax_sum) )
       tranatal_sum(1:imax_sum) = 0.0d0

       resc2ANATAL(:) = 0.0d0
       resc3ANATAL(:) = 0.0d0

       end subroutine allocate_anatal
!------------------------------------------------------------------------

!------------------------------------------------------------------------
       subroutine copy_anatal
!
!      called form anatally subroutine
!........................................................................
       use RESTALMOD, only : lrestalm, mrestalm_sum
       use TALMOD,   only : tr0, tr0_sum
       use sumtallymod, only: nfile

       include 'param.inc'
       common /tall00/ itnm, ital(itlmax), itals(itlmax), italm(itlmax)

       deallocate(tr0)
       deallocate(tr0_sum)
       imax = 1
       imax_sum = 1
       do m=1, itnm
        imax=imax+lrestalm(m) * nfile(m)
        imax_sum = imax_sum + mrestalm_sum(m) * nfile(m)
       enddo

       allocate( tr0(imax) )
       allocate( tr0_sum(imax_sum) )

       tr0(:)   = tranatal(:)

       tr0_sum(:)   = tranatal_sum(:)

       end subroutine copy_anatal
!------------------------------------------------------------------------


!------------------------------------------------------------------------
       subroutine DEALLOCATE_ANATAL
!
!      called form anatally subroutine
!........................................................................
       deallocate( trANATAL )
       deallocate( trANATAL_SUM )
       deallocate( lanatalm ) ! S.H. 2020.6.17
       deallocate( lanatalm_sum )
       deallocate( ianatalm_sum )
       deallocate( manatalm_sum )

       deallocate( nanatalRead )  ! set_anatal subroutine
       deallocate( manatally   )
       deallocate(ianataldchain)
       deallocate( tallyfolder )
       deallocate(ltallyfolder )
       deallocate( tallyfile   )
       deallocate(ltallyfile   )


       end subroutine DEALLOCATE_ANATAL
!------------------------------------------------------------------------


!***********************************************************************
      end module anatallymod
!***********************************************************************
