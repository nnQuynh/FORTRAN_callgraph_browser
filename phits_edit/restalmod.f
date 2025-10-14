************************************************************************
      module RESTALMOD
*                                                                      *
*     - Buffer module to hold tally value of restart initial,          *
*       on master process when MPI parallel execution.                 *
*                                                                      *
*     by D.OBINATA on 2012.6.6                                         *
*                                                                      *
************************************************************************
      use partmod, only: itpan, itpat, jtpat ! frtati 2021/10/05
      use talmod, only: italsize_sum, mtalsize_sum, italhead_sum


!------------------------------------------------------------------------

      implicit double precision (a-h,o-z)
      real(8),allocatable,save:: trRES(:)
      integer,allocatable,save:: irestalm(:)
      integer,allocatable,save,target:: lrestalm(:)

      real(8),allocatable,save,target:: trRES_sum(:)
      integer,allocatable,save,target:: irestalm_sum(:,:)
      integer,allocatable,save,target:: mrestalm_sum(:)
      integer,allocatable,save,target:: lrestalm_sum(:,:)

      contains
!------------------------------------------------------------------------
      subroutine ALLOCATE_RESTAL
      include 'param.inc'
      common /tall00/ itnm, ital(itlmax), itals(itlmax), italm(itlmax)
      common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                itout(itlmax), ittwo(itlmax)
      common /tall02/ itrgn(itlmax), itrgm(itlmax), itreg(itlmax)
      common /tall03/ itrty(itlmax), itrnm(itlmax), itrrg(itlmax),
     &                rtrmi(itlmax), rtrma(itlmax), rtrdl(itlmax)
      common /tall04/ itxty(itlmax), itxnm(itlmax), itxrg(itlmax),
     &                rtxmi(itlmax), rtxma(itlmax), rtxdl(itlmax)
      common /tall05/ ityty(itlmax), itynm(itlmax), ityrg(itlmax),
     &                rtymi(itlmax), rtyma(itlmax), rtydl(itlmax)
      common /tall06/ itzty(itlmax), itznm(itlmax), itzrg(itlmax),
     &                rtzmi(itlmax), rtzma(itlmax), rtzdl(itlmax)
      common /tall07/ itety(itlmax), itenm(itlmax), iterg(itlmax),
     &                rtemi(itlmax), rtema(itlmax), rtedl(itlmax)
      common /tall12/ itman(itlmax), itmat(itlmax), itmct(itlmax),
     &                itnun(itlmax), itnuc(itlmax),
     &                itndz(itlmax), itndn(itlmax),
     &                itnkz(itlmax), itnkn(itlmax)
      common /tall13/ itrcn(itlmax), itrcr(itlmax), itrca(itlmax)
      common /tall17/ itndy(itlmax)
      common /tall20/ ittty(itlmax), ittnm(itlmax), ittrg(itlmax),
     &                rttmi(itlmax), rttma(itlmax), rttdl(itlmax)
      common /tall23/ itaty(itlmax), itanm(itlmax), itarg(itlmax),
     &                rtami(itlmax), rtama(itlmax), rtadl(itlmax)
      common /tall45/ itmlp(itlmax), itmln(itlmax,6), itmst(itlmax),
     &                itmli(itlmax,6), rtmme(itlmax,6), itmnt(itlmax,6),
     &                itmpn(itlmax,6), itmpt(itlmax,6,6,2)
      common /tall56/ itety2(itlmax), itenm2(itlmax), iterg2(itlmax),
     &                rtemi2(itlmax), rtema2(itlmax), rtedl2(itlmax)
      common /tall72/ itactnm(itlmax), itactrg(itlmax), itactmax(itlmax)      ! S.Abe 2018/02/15
      common /tall83/ itnzn(itlmax), itndm(itlmax) ! frtati 2022/02/18

      common /tall91/ itextstat(itlmax), mftal(itlmax) ! T.Sato 2024/05/17


      integer imax,m
      integer imax_sum,iax,mnmax_sum

      allocate( irestalm(max(itnm,1)) )
      allocate( lrestalm(max(itnm,1)) )
      irestalm(:) = 0
      lrestalm(:) = 0

      allocate( irestalm_sum(max(itnm,1),6) )
      allocate( mrestalm_sum(max(itnm,1)) )
      allocate( lrestalm_sum(max(itnm,1),6) )
      irestalm_sum(:,:) = 0
      lrestalm_sum(:,:) = 0
      mrestalm_sum(:) = 0

      imax=1

      do m=1, itnm

*-----------------------------------------------------------------------
*        t-track tally
*-----------------------------------------------------------------------

        if( ital(m) .eq. 1 ) then
          if ( itmsh(m) .eq. 1 ) then   ! ttracreg
            lrestalm(m) = itenm(m) * itpan(m) * itrgn(m)
     &                  * itmst(m) * ittnm(m) * mftal(m)
            ! trRES(np,ne,nt,nr,nm)

          else if ( itmsh(m) .eq. 2 ) then   ! ttracrz
            lrestalm(m) = itenm(m) * itpan(m) * mftal(m)
     &                  * itrnm(m) * itznm(m) * itanm(m) * itmst(m)
     &                  * ittnm(m)
            ! trRES(np,ne,nt,nr*nz,na,nm)


          else if ( itmsh(m) .eq. 3 ) then   ! ttracxyz
            lrestalm(m) = itenm(m) * itpan(m) * mftal(m)
     &                  * itxnm(m) * itynm(m) * itznm(m) * itmst(m)
     &                  * ittnm(m)
            ! trRES(np,ne,nt,nx*ny*nz,nm)
          else if ( itmsh(m) .eq. 4 ) then ! ttractet
            lrestalm(m) = itenm(m) * itpan(m) * itrgn(m)
     &                  * itmst(m) * ittnm(m) * mftal(m)
            ! trRES(np,ne,nt,nr,nm)
          endif

*-----------------------------------------------------------------------
*        t-adjoint tally
*-----------------------------------------------------------------------

        else if( ital(m) .eq. 19 ) then
          if ( itmsh(m) .eq. 1 ) then   ! tadjreg
            lrestalm(m) = itenm(m) * itpan(m) * itrgn(m)
     &                  * itmst(m) * ittnm(m) * 2
            ! trRES(np,ne,nt,nr,nm)

          else if ( itmsh(m) .eq. 2 ) then   ! tadjrz
            lrestalm(m) = itenm(m) * itpan(m) * 2
     &                  * itrnm(m) * itznm(m) * itanm(m) * itmst(m)
     &                  * ittnm(m)
            ! trRES(np,ne,nt,nr*nz,na,nm)

          else if ( itmsh(m) .eq. 3 ) then   ! tadjxyz
            lrestalm(m) = itenm(m) * itpan(m) * 2
     &                  * itxnm(m) * itynm(m) * itznm(m) * itmst(m)
     &                  * ittnm(m)
            ! trRES(np,ne,nt,nx*ny*nz,nm)
          endif

*-----------------------------------------------------------------------
*        t-cross tally
*-----------------------------------------------------------------------

        else if( ital(m) .eq. 2 ) then
          if ( itmsh(m) .eq. 1 ) then   ! tcrsreg
            lrestalm(m) = itenm(m) * itpan(m) * 2
     &                  * itrcn(m) * itanm(m) * itmst(m)
     &                  * ittnm(m)
          else if ( itmsh(m) .eq. 2 ) then   ! tcrsrz
            lrestalm(m) = itenm(m) * itpan(m) * 2
     &                  * ( itrnm(m) + 1 ) * itznm(m)
     &                  * itanm(m) * itmst(m)
     &                  * ittnm(m)
     &                  + itenm(m) * itpan(m) * 2
     &                  * ( itznm(m) + 1 ) * itrnm(m)
*    &                  * itanm(m)
*nais
     &                  * itanm(m) * itmst(m)
     &                  * ittnm(m)
          else if ( itmsh(m) .eq. 3 ) then   ! tcrsxyz
            lrestalm(m) = itenm(m) * itpan(m) * 2
     &                  * itxnm(m) * itynm(m) * ( itznm(m) + 1 )
     &                  * itanm(m) * itmst(m)
     &                  * ittnm(m)
          endif

*-----------------------------------------------------------------------
*         [t-yield]
*-----------------------------------------------------------------------
        else if ( ital(m) .eq. 3 ) then
cfrtati 2022/02/18 3 -> itndm+1 dimension for isomer states
          if ( itmsh(m) .eq. 1 ) then   ! tyildreg
            lrestalm(m) = itndz(m) * itndn(m) * 2
     &                  * itrgn(m) * (itndm(m)+1)
          else if ( itmsh(m) .eq. 2 ) then   ! tyildrz
            lrestalm(m) = itndz(m) * itndn(m) * 2
     &                  * itrnm(m) * itznm(m) * (itndm(m)+1)
          else if ( itmsh(m) .eq. 3 ) then   ! tyildxyz
            lrestalm(m) = itndz(m) * itndn(m) * 2
     &                  * itxnm(m) * itynm(m) * itznm(m) * (itndm(m)+1)
          else if ( itmsh(m) .eq. 4 ) then   ! tyildtet
            lrestalm(m) = itndz(m) * itndn(m) * 2
     &                  * itrgn(m) * (itndm(m)+1)
          endif

*-----------------------------------------------------------------------
*        t-heat tally
*-----------------------------------------------------------------------
        else if( ital(m) .eq. 4 ) then
          if( itmsh(m) .eq. 1 ) then    ! thetreg
            lrestalm(m) = itndy(m) * 2 * ( itenm(m) + 1 )
     &                  * itrgn(m) + 5
            ! trRES(nd,nr)
          else if(itmsh(m) .eq. 2) then ! thetrz
            lrestalm(m) = itndy(m) * 2 * ( itenm(m) + 1 )
     &                  * itrnm(m) * itznm(m) + 5
            ! trRES(nd,nr,nz)
          else if(itmsh(m) .eq. 3) then ! thetxyz
            lrestalm(m) = itndy(m) * 2 * ( itenm(m) + 1 )
     &                  * itxnm(m) * itynm(m) * itznm(m) + 5
            ! trRES(nd,nx,ny,nz)
          endif

*-----------------------------------------------------------------------
*         [t-star]
*-----------------------------------------------------------------------
        else if ( ital(m) .eq. 5 ) then
          if( itmsh(m) .eq. 1 ) then ! tstarreg
            lrestalm(m) = itenm(m) * ittnm(m) * itpan(m) * 2
     &                  * itrgn(m)
     &                  * itactnm(m)      ! S.Abe 2018/02/15
          ! trRES(np,ne,nt,nr,2)
          else if( itmsh(m) .eq. 2 ) then ! tstarrz
            lrestalm(m) = itenm(m) * itpan(m) * 2
     &                  * ittnm(m)
     &                  * itrnm(m) * itznm(m)
          ! trRES(np,ne,nt,nr,nz,2)
          else if( itmsh(m) .eq. 3 ) then ! tstarxyz
            lrestalm(m) = itenm(m) * itpan(m) * 2
     &                  * ittnm(m)
     &                  * itxnm(m) * itynm(m) * itznm(m)
          ! trRES(np,ne,nt,nx,ny,nz,2)
          endif

*-----------------------------------------------------------------------
*         [t-time]
*-----------------------------------------------------------------------
        else if ( ital(m) .eq. 6 ) then
          if( itmsh(m) .eq. 1 ) then ! ttimereg
            lrestalm(m) = itenm(m) * itpan(m) * 2
     &                  * ittnm(m)
     &                  * itrgn(m)
          ! trRES(np,ne,nt,nr,2)
          else if( itmsh(m) .eq. 2 ) then ! ttimedrz
            lrestalm(m) = itenm(m) * itpan(m) * 2
     &                  * ittnm(m)
     &                  * itrnm(m) * itznm(m)
          ! trRES(np,ne,nt,nr,nz,2)
          else if( itmsh(m) .eq. 3 ) then ! ttimexyz
            lrestalm(m) = itenm(m) * itpan(m) * 2
     &                  * ittnm(m)
     &                  * itxnm(m) * itynm(m) * itznm(m)
          ! trRES(np,ne,nt,nx,ny,nz,2)
          endif

*-----------------------------------------------------------------------
*         [t-dpa]
*-----------------------------------------------------------------------
        else if ( ital(m) .eq. 7 ) then
          if( itmsh(m) .eq. 1 ) then ! tdpareg
            lrestalm(m) = itndy(m) * itpan(m) * 2
     &                  * itrgn(m)
          ! trRES(nd,np,nr,2)
          else if( itmsh(m) .eq. 2 ) then ! tdpadrz
            lrestalm(m) = itndy(m) * itpan(m) * 2
     &                  * itrnm(m) * itznm(m)
          ! trRES(nd,np,nr,nz,2)
          else if( itmsh(m) .eq. 3 ) then ! tdpaxyz
            lrestalm(m) = itndy(m) * itpan(m) * 2
     &                  * itxnm(m) * itynm(m) * itznm(m)
          ! trRES(nd,np,nx,ny,nz,2)
          else if( itmsh(m) .eq. 4 ) then ! tdpatet
            lrestalm(m) = itndy(m) * itpan(m) * 2
     &                  * itrgn(m)
          ! trRES(nd,np,nr,2)
          endif

*-----------------------------------------------------------------------
*         [t-product]
*-----------------------------------------------------------------------
        else if ( ital(m) .eq. 8 ) then
          if( itmsh(m) .eq. 1 ) then ! tpdctreg
            lrestalm(m) = itenm(m) * ittnm(m) * itpan(m) * 2
     &                  * itrgn(m) * itanm(m)
          ! trRES(np,ne,nt,na,nr,2)
          else if( itmsh(m) .eq. 2 ) then ! tpdctdrz
            lrestalm(m) = itenm(m) * itpan(m) * 2
     &                  * ittnm(m) * itanm(m)
     &                  * itrnm(m) * itznm(m)
          ! trRES(np,ne,nt,na,nr,nz,2)
          else if( itmsh(m) .eq. 3 ) then ! tpdctxyz
            lrestalm(m) = itenm(m) * itpan(m) * 2
     &                  * ittnm(m) * itanm(m)
     &                  * itxnm(m) * itynm(m) * itznm(m)
          ! trRES(np,ne,nt,na,nx*ny*nz,2)
          else if( itmsh(m) .eq. 4 ) then ! tpdcttet
            lrestalm(m) = itenm(m) * ittnm(m) * itpan(m) * 2
     &                  * itrgn(m) * itanm(m)
          ! trRES(np,ne,nt,na,nr,2)
          endif

*-----------------------------------------------------------------------
*         [t-let]
*-----------------------------------------------------------------------
          else if ( ital(m) .eq. 12 ) then
            if( itmsh(m) .eq. 1 ) then ! tletreg
              lrestalm(m) = itenm(m) * itpan(m) * 2
     &                    * itrgn(m)
            ! trRES(np,ne,nr)
            else if( itmsh(m) .eq. 2 ) then ! tletdrz
              lrestalm(m) = itenm(m) * itpan(m) * 2
     &                    * itrnm(m) * itznm(m)
            ! trRES(np,ne,nr,nz)
            else if( itmsh(m) .eq. 3 ) then ! tletxyz
              lrestalm(m) = itenm(m) * itpan(m) * 2
     &                    * itxnm(m) * itynm(m) * itznm(m)
            ! trRES(np,ne,nx*ny*nz)
            endif


*-----------------------------------------------------------------------
*        t-deposit tally
*-----------------------------------------------------------------------
        else if( ital(m) .eq. 13 ) then
          if( itmsh(m) .eq. 1 ) then    ! tdepstreg
            lrestalm(m) = ( itenm(m) + 1 ) * itpan(m) * 2
     &                  * itrgn(m) * ittnm(m)
            ! trRES(np,ne,nr,nt)
          else if(itmsh(m) .eq. 2) then ! tdepstrz
            lrestalm(m) = ( itenm(m) + 1 ) * itpan(m) * 2
     &                  * ittnm(m)
     &                  * itrnm(m) * itznm(m)
            ! trRES(np,ne,nt,nr,nz)
          else if(itmsh(m) .eq. 3) then ! tdepstxyz
            lrestalm(m) = ( itenm(m) + 1 ) * itpan(m) * 2
     &                  * ittnm(m)
     &                  * itxnm(m) * itynm(m) * itznm(m)
            ! trRES(np,ne,nt,nx*ny*nz)
          else if( itmsh(m) .eq. 4 ) then ! tdepsttet
            lrestalm(m) = ( itenm(m) + 1 ) * itpan(m) * 2
     &                  * itrgn(m) * ittnm(m)
            ! trRES(np,ne,nr,nt)
          endif

*-----------------------------------------------------------------------
*        t-deposit2 tally
*-----------------------------------------------------------------------
        else if( ital(m) .eq. 14 ) then
          if( itmsh(m) .eq. 1 ) then    ! tdpst2reg
            lrestalm(m) = ( itenm(m) + 1 ) * ( itenm2(m) + 1 )
     &                  * itpan(m) * ittnm(m) * 2
            ! trRES(np,0:e1,0:e2,nt,2)
          endif

*-----------------------------------------------------------------------
*         [t-sed]
*-----------------------------------------------------------------------
        else if ( ital(m) .eq. 15 ) then
          if( itmsh(m) .eq. 1 ) then ! tletreg
            lrestalm(m) = itenm(m) * itpan(m) * 2
     &                  * itrgn(m)
            ! trRES(np,ne,nr)
          else if( itmsh(m) .eq. 2 ) then ! tletdrz
            lrestalm(m) = itenm(m) * itpan(m) * 2
     &                  * itrnm(m) * itznm(m)
            ! trRES(np,ne,nr,nz)
          else if( itmsh(m) .eq. 3 ) then ! tletxyz
            lrestalm(m) = itenm(m) * itpan(m) * 2
     &                  * itxnm(m) * itynm(m) * itznm(m)
            ! trRES(np,ne,nx*ny*nz)
          end if

*-----------------------------------------------------------------------
*         [t-point]
*-----------------------------------------------------------------------

        else if ( ital(m) .eq. 17 ) then
            lrestalm(m) = itenm(m) * itpan(m) * itmsh(m)
     &                  * itmst(m) * ittnm(m) * mftal(m)
            ! trRES(np,ne,nt,nr,nm)

*-----------------------------------------------------------------------
*         [t-wwg]
*-----------------------------------------------------------------------

        else if ( ital(m) .eq. 18 ) then
          if( itmsh(m) .eq. 1 ) then
            lrestalm(m) = itenm(m) * itpan(m) * itrgn(m)
     &                  * itmst(m) * ittnm(m) * mftal(m)  ! T.Sato 2024/05/17
            ! trRES(np,ne,nt,nr,nm)
          else if ( itmsh(m) .eq. 3 ) then
            lrestalm(m) = itenm(m) * itpan(m) * mftal(m)  ! T.Sato 2024/05/17
     &                  * itxnm(m) * itynm(m) * itznm(m) * itmst(m)
     &                  * ittnm(m)
            ! trRES(np,ne,nt,nx*ny*nz,nm)
cFURUTA20240110
          else if ( itmsh(m) .eq. 4 ) then ! ttractet
            lrestalm(m) = itenm(m) * itpan(m) * itrgn(m)
     &                  * itmst(m) * ittnm(m) * mftal(m) ! T.Sato 2024/05/17
            ! trRES(np,ne,nt,nr,nm)
          endif

*-----------------------------------------------------------------------
*         [t-volume]
*-----------------------------------------------------------------------

        else if ( ital(m) .eq. 21 ) then
            lrestalm(m) = itrgn(m) * 2
            ! trRES(nr)

*-----------------------------------------------------------------------
*         [t-wwbg]
*-----------------------------------------------------------------------

        else if ( ital(m) .eq. 22 ) then
          if ( itmsh(m) .eq. 1 ) then
            lrestalm(m) = itrgn(m) * 2 * 2
            ! trRES(nr)
          else if ( itmsh(m) .eq. 3 ) then
            lrestalm(m) = itxnm(m) * itynm(m) * itznm(m) * 2 * 2
            ! trRES(nx*ny*nz,2)
          else if ( itmsh(m) .eq. 4 ) then
            lrestalm(m) = itrgn(m) * 2 * 2
            ! trRES(nr,2)
          endif

*-----------------------------------------------------------------------
        endif

        irestalm(m)=imax
        imax=imax+lrestalm(m)

      enddo


       imax_sum=1

       if(itnm > 0) then

         lrestalm_sum(:,:) = italsize_sum(:,:)
         irestalm_sum(:,:) = italhead_sum(:,:)
         mrestalm_sum(:)  = mtalsize_sum(:)

         imax_sum= imax_sum + sum(mrestalm_sum(:))

       endif

      allocate( trRES(imax) )
      trRES(1:imax)=0.0d0

      allocate( trRES_sum(imax_sum) )
      trRES_sum(1:imax_sum)=0.0d0

      end subroutine ALLOCATE_RESTAL
!------------------------------------------------------------------------
      subroutine DEALLOCATE_RESTAL

      deallocate( trRES )
      deallocate( irestalm )
      deallocate( lrestalm )

      deallocate( trRES_sum )
      deallocate( irestalm_sum )
      deallocate( lrestalm_sum )
      deallocate( mrestalm_sum )

      end subroutine DEALLOCATE_RESTAL



      subroutine CALC_TALSIZE_SUM(m,iax,mnmax)



      end subroutine

      subroutine GET_TRRES_HEAD_POINTER_SUM(p_sum,m,iax)

        implicit none
        real(8),pointer,intent(out):: p_sum(:)
        integer,intent(in):: m,iax

        p_sum => trRES_sum(irestalm_sum(m,iax):)

      end subroutine

!------------------------------------------------------------------------
      end module RESTALMOD

