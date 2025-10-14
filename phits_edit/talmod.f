      module TALMOD


      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      real(8),allocatable,save,target:: tr0(:)
      integer,allocatable,save,target:: italsize(:)
      integer,allocatable,save,target:: italhead(:)
!$    real(8),allocatable,save,target:: tr0ref(:)
      real(8),allocatable,save,target:: deist(:,:) ! S.H. 2022.12.23
      real(8),allocatable,save,target:: prodenmnmx(:,:) ! S.H. 2024.4.3
      integer,allocatable,save,target:: iprodenhead(:)

* sumover
      real(8),allocatable,save,target:: tr0_sum(:)
      integer,allocatable,save,target:: italsize_sum(:,:)
      integer,allocatable,save,target:: italsize_2_sum(:,:)
      integer,allocatable,save,target:: mtalsize_sum(:)
      integer,allocatable,save,target:: italhead_sum(:,:)
!$    real(8),allocatable,save,target:: tr0ref_sum(:)

*-----------------------------------------------------------------------
      contains

      subroutine INIT_TALMOD

      end subroutine

      subroutine GET_TR_HEAD_POINTER(p,m)

        implicit none
        real(8),pointer,intent(out):: p(:)
        integer,intent(in):: m

        p => tr0(italhead(m):)

      end subroutine


      subroutine GET_TR_HEAD_POINTER_SUM(p_sum,m,iax)

        implicit none
        real(8),pointer,intent(out):: p_sum(:)
        integer,intent(in):: m,iax
        real(8),pointer :: pp(:)

        p_sum => tr0_sum(italhead_sum(m,iax):)

      end subroutine

      subroutine GET_TZ_HEAD_POINTER_SUM(p_sum,m,iax)

        implicit none
        real(8),pointer,intent(out):: p_sum(:)
        integer,intent(in):: m,iax
        real(8),pointer :: pp(:)

        p_sum => tr0_sum(italhead_sum(m,iax)+italsize_2_sum(m,iax):)

      end subroutine

      subroutine GET_TR_HEAD_POINTER_SUM_NTF(p_sum,m,ntf,iax)

        implicit none
        real(8),pointer,intent(out):: p_sum(:)
        integer,intent(in):: m,iax,ntf
        integer :: ntfbase
        real(8),pointer :: pp(:)

        ntfbase = mtalsize_sum(m) * (ntf-1)
        p_sum => tr0_sum(italhead_sum(m,iax) + ntfbase:)

      end subroutine

      subroutine GET_TZ_HEAD_POINTER_SUM_NTF(p_sum,m,ntf,iax)

        implicit none
        real(8),pointer,intent(out):: p_sum(:)
        integer,intent(in):: m,iax,ntf
        integer :: ntfbase
        real(8),pointer :: pp(:)

        ntfbase = mtalsize_sum(m) * (ntf-1)
        p_sum => 
     &  tr0_sum(italhead_sum(m,iax) + italsize_2_sum(m,iax) + ntfbase:)

      end subroutine

      subroutine CALC_TALSIZE(m,mnmax)

        use partmod, only: itpan, itpat, jtpat ! frtati 2021/10/05

        integer,intent(in):: m
        integer,intent(out):: mnmax

        include 'param.inc'
        common /tall00/ itnm, ital(itlmax), itals(itlmax), italm(itlmax)
        common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                  itout(itlmax), ittwo(itlmax)
        common /tall02/ itrgn(itlmax), itrgm(itlmax), itreg(itlmax)
        common /tall03/ itrty(itlmax), itrnm(itlmax), itrrg(itlmax),
     &                  rtrmi(itlmax), rtrma(itlmax), rtrdl(itlmax)
        common /tall04/ itxty(itlmax), itxnm(itlmax), itxrg(itlmax),
     &                  rtxmi(itlmax), rtxma(itlmax), rtxdl(itlmax)
        common /tall05/ ityty(itlmax), itynm(itlmax), ityrg(itlmax),
     &                  rtymi(itlmax), rtyma(itlmax), rtydl(itlmax)
        common /tall06/ itzty(itlmax), itznm(itlmax), itzrg(itlmax),
     &                  rtzmi(itlmax), rtzma(itlmax), rtzdl(itlmax)
        common /tall07/ itety(itlmax), itenm(itlmax), iterg(itlmax),
     &                  rtemi(itlmax), rtema(itlmax), rtedl(itlmax)
        common /tall12/ itman(itlmax), itmat(itlmax), itmct(itlmax),
     &                  itnun(itlmax), itnuc(itlmax),
     &                  itndz(itlmax), itndn(itlmax),
     &                  itnkz(itlmax), itnkn(itlmax)
        common /tall13/ itrcn(itlmax), itrcr(itlmax), itrca(itlmax)
        common /tall17/ itndy(itlmax)
        common /tall20/ ittty(itlmax), ittnm(itlmax), ittrg(itlmax),
     &                  rttmi(itlmax), rttma(itlmax), rttdl(itlmax)
        common /tall23/ itaty(itlmax), itanm(itlmax), itarg(itlmax),
     &                  rtami(itlmax), rtama(itlmax), rtadl(itlmax)
        common /tall45/ itmlp(itlmax), itmln(itlmax,6), itmst(itlmax),
     &                  itmli(itlmax,6),rtmme(itlmax,6),itmnt(itlmax,6),
     &                  itmpn(itlmax,6), itmpt(itlmax,6,6,2)
        common /tall56/ itety2(itlmax), itenm2(itlmax), iterg2(itlmax),
     &                  rtemi2(itlmax), rtema2(itlmax), rtedl2(itlmax)
        common /tall72/ itactnm(itlmax),itactrg(itlmax),itactmax(itlmax)      ! S.Abe 2018/02/15
        common /tall83/ itnzn(itlmax), itndm(itlmax) ! frtati 2022/02/18

      common /tall91/ itextstat(itlmax), mftal(itlmax) ! S.H. extstat 2024.3.18

*-----------------------------------------------------------------------

        mnmax = 0

*-----------------------------------------------------------------------
*        t-track tally
*-----------------------------------------------------------------------
        if( ital(m) .eq. 1 ) then
          if ( itmsh(m) .eq. 1 ) then   ! ttracreg
            mnmax = itenm(m) * itpan(m) * itrgn(m)
     &            * itmst(m) * ittnm(m) * mftal(m)
            ! tr(np,ne,nt,nr,nm)


          else if ( itmsh(m) .eq. 2 ) then   ! ttracrz
            mnmax = itenm(m) * itpan(m) * mftal(m)
     &            * itrnm(m) * itznm(m) * itanm(m) * itmst(m)
     &            * ittnm(m)
            ! tr(np,ne,nt,nr,nz,na,nm)


          else if ( itmsh(m) .eq. 3 ) then   ! ttracxyz
            mnmax = itenm(m) * itpan(m) * mftal(m)
     &            * itxnm(m) * itynm(m) * itznm(m) * itmst(m)
     &            * ittnm(m)
            ! tr(np,ne,nt,nx*ny*nz,nm)
          else if ( itmsh(m) .eq. 4 ) then ! ttractet
           mnmax = itenm(m) * itpan(m) * itrgn(m)
     &            * itmst(m) * ittnm(m) * mftal(m)
            ! tr(np,ne,nt,nr,nm)
          endif

*-----------------------------------------------------------------------
*        t-adjoint tally
*-----------------------------------------------------------------------

        else if( ital(m) .eq. 19 ) then
          if ( itmsh(m) .eq. 1 ) then   ! tadjreg
            mnmax = itenm(m) * itpan(m) * itrgn(m)
     &            * itmst(m) * ittnm(m) * 2
            ! tr(np,ne,nt,nr,nm)

          else if ( itmsh(m) .eq. 2 ) then   ! tadjrz
            mnmax = itenm(m) * itpan(m) * 2
     &            * itrnm(m) * itznm(m) * itanm(m) * itmst(m)
     &            * ittnm(m)
            ! tr(np,ne,nt,nr,nz,na,nm)

          else if ( itmsh(m) .eq. 3 ) then   ! tadjxyz
            mnmax = itenm(m) * itpan(m) * 2
     &            * itxnm(m) * itynm(m) * itznm(m) * itmst(m)
     &            * ittnm(m)
            ! tr(np,ne,nt,nx*ny*nz,nm)
          endif

*-----------------------------------------------------------------------
*        t-cross tally
*-----------------------------------------------------------------------

        else if( ital(m) .eq. 2 ) then
          if ( itmsh(m) .eq. 1 ) then   ! tcrsreg
            mnmax = itenm(m) * itpan(m) * 2
     &            * itrcn(m) * itanm(m) * itmst(m)
     &            * ittnm(m)
          else if ( itmsh(m) .eq. 2 ) then   ! tcrsrz
            mnmax = itenm(m) * itpan(m) * 2
     &            * ( itrnm(m) + 1 ) * itznm(m)
     &            * itanm(m) * itmst(m)
     &            * ittnm(m)
     &            + itenm(m) * itpan(m) * 2
     &            * ( itznm(m) + 1 ) * itrnm(m)
     &            * itanm(m) * itmst(m)
     &            * ittnm(m)
          else if ( itmsh(m) .eq. 3 ) then   ! tcrsxyz
            mnmax = itenm(m) * itpan(m) * 2
     &            * itxnm(m) * itynm(m) * ( itznm(m) + 1 )
     &            * itanm(m) * itmst(m)
     &            * ittnm(m)
          endif

*-----------------------------------------------------------------------
*         [t-yield]
*-----------------------------------------------------------------------
        else if ( ital(m) .eq. 3 ) then
cfrtati 2022/02/18 3 -> itndm+1 dimension for isomer states
          if ( itmsh(m) .eq. 1 ) then   ! tyildreg
            mnmax = itndz(m) * itndn(m) * 2
     &            * itrgn(m) * (itndm(m)+1)
          else if ( itmsh(m) .eq. 2 ) then   ! tyildrz
            mnmax = itndz(m) * itndn(m) * 2
     &            * itrnm(m) * itznm(m) * (itndm(m)+1)
          else if ( itmsh(m) .eq. 3 ) then   ! tyildxyz
            mnmax = itndz(m) * itndn(m) * 2
     &            * itxnm(m) * itynm(m) * itznm(m) * (itndm(m)+1)
          else if ( itmsh(m) .eq. 4 ) then   ! tyildtet
            mnmax = itndz(m) * itndn(m) * 2
     &            * itrgn(m) * (itndm(m)+1)
          endif

*-----------------------------------------------------------------------
*        t-heat tally
*-----------------------------------------------------------------------
        else if( ital(m) .eq. 4 ) then
          if( itmsh(m) .eq. 1 ) then    ! thetreg
            mnmax = itndy(m) * 2 * ( itenm(m) + 1 )
     &            * itrgn(m) + 5
            ! tr(nd,nr)
          else if(itmsh(m) .eq. 2) then ! thetrz
            mnmax = itndy(m) * 2 * ( itenm(m) + 1 )
     &            * itrnm(m) * itznm(m) + 5
            ! tr(nd,nr,nz)
          else if(itmsh(m) .eq. 3) then ! thetxyz
            mnmax = itndy(m) * 2 * ( itenm(m) + 1 )
     &            * itxnm(m) * itynm(m) * itznm(m) + 5
            ! tr(nd,nx,ny,nz)
          endif

*-----------------------------------------------------------------------
*         [t-star]
*-----------------------------------------------------------------------
        else if ( ital(m) .eq. 5 ) then
          if( itmsh(m) .eq. 1 ) then ! tstarreg
            mnmax = itenm(m) * ittnm(m) * itpan(m) * 2
     &            * itrgn(m)
     &            * itactnm(m)      ! S.Abe 2018/02/15
          ! tr(np,ne,nt,nr,2)
          else if( itmsh(m) .eq. 2 ) then ! tstarrz
            mnmax = itenm(m) * itpan(m) * 2
     &            * ittnm(m)
     &            * itrnm(m) * itznm(m)
          ! tr(np,ne,nt,nr,nz,2)
          else if( itmsh(m) .eq. 3 ) then ! tstarxyz
            mnmax = itenm(m) * itpan(m) * 2
     &            * ittnm(m)
     &            * itxnm(m) * itynm(m) * itznm(m)
          ! tr(np,ne,nt,nx,ny,nz,2)
          endif

*-----------------------------------------------------------------------
*         [t-time]
*-----------------------------------------------------------------------
        else if ( ital(m) .eq. 6 ) then
          if( itmsh(m) .eq. 1 ) then ! ttimereg
            mnmax = itenm(m) * itpan(m) * 2
     &            * ittnm(m)
     &            * itrgn(m)
          ! tr(np,ne,nt,nr,2)
          else if( itmsh(m) .eq. 2 ) then ! ttimedrz
            mnmax = itenm(m) * itpan(m) * 2
     &            * ittnm(m)
     &            * itrnm(m) * itznm(m)
          ! tr(np,ne,nt,nr,nz,2)
          else if( itmsh(m) .eq. 3 ) then ! ttimexyz
            mnmax = itenm(m) * itpan(m) * 2
     &            * ittnm(m)
     &            * itxnm(m) * itynm(m) * itznm(m)
          ! tr(np,ne,nt,nx,ny,nz,2)
          endif

*-----------------------------------------------------------------------
*         [t-dpa]
*-----------------------------------------------------------------------
        else if ( ital(m) .eq. 7 ) then
          if( itmsh(m) .eq. 1 ) then ! tdpareg
            mnmax = itndy(m) * itpan(m) * 2
     &            * itrgn(m)
          ! tr(nd,np,nr,2)
          else if( itmsh(m) .eq. 2 ) then ! tdpadrz
            mnmax = itndy(m) * itpan(m) * 2
     &            * itrnm(m) * itznm(m)
          ! tr(nd,np,nr,nz,2)
          else if( itmsh(m) .eq. 3 ) then ! tdpaxyz
            mnmax = itndy(m) * itpan(m) * 2
     &            * itxnm(m) * itynm(m) * itznm(m)
          ! tr(nd,np,nx,ny,nz,2)
          else if( itmsh(m) .eq. 4 ) then ! tdpatet
            mnmax = itndy(m) * itpan(m) * 2
     &            * itrgn(m)
          ! tr(nd,np,nr,2)
          endif

*-----------------------------------------------------------------------
*         [t-product]
*-----------------------------------------------------------------------
        else if ( ital(m) .eq. 8 ) then
          if( itmsh(m) .eq. 1 ) then ! tprodreg
            mnmax = itenm(m) * ittnm(m) * itpan(m) * 2
     &            * itrgn(m) * itanm(m)
          ! tr(np,ne,nt,na,nr,2)
          else if( itmsh(m) .eq. 2 ) then ! tproddrz
            mnmax = itenm(m) * itpan(m) * 2
     &            * ittnm(m) * itanm(m)
     &            * itrnm(m) * itznm(m)
          ! tr(np,ne,nt,na,nr,nz,2)
          else if( itmsh(m) .eq. 3 ) then ! tprodxyz
            mnmax = itenm(m) * itpan(m) * 2
     &            * ittnm(m) * itanm(m)
     &            * itxnm(m) * itynm(m) * itznm(m)
          ! tr(np,ne,nt,na,nx*ny*nz,2)
          else if( itmsh(m) .eq. 4 ) then ! tprodtet
            mnmax = itenm(m) * ittnm(m) * itpan(m) * 2
     &            * itrgn(m) * itanm(m)
          ! tr(np,ne,nt,na,nr,2)
          endif

*-----------------------------------------------------------------------
*         [g-show]
*-----------------------------------------------------------------------
        else if ( ital(m) .eq. 9 ) then

*-----------------------------------------------------------------------
*         [r-show]
*-----------------------------------------------------------------------
        else if ( ital(m) .eq. 10 ) then

*-----------------------------------------------------------------------
*         [3d-show]
*-----------------------------------------------------------------------
        else if ( ital(m) .eq. 11 ) then

*-----------------------------------------------------------------------
*         [t-let]
*-----------------------------------------------------------------------
        else if ( ital(m) .eq. 12 ) then
          if( itmsh(m) .eq. 1 ) then ! tletreg
            mnmax = itenm(m) * itpan(m) * 2
     &            * itrgn(m)
          ! tr(np,ne,nr)
          else if( itmsh(m) .eq. 2 ) then ! tletdrz
            mnmax = itenm(m) * itpan(m) * 2
     &            * itrnm(m) * itznm(m)
          ! tr(np,ne,nr,nz)
          else if( itmsh(m) .eq. 3 ) then ! tletxyz
            mnmax = itenm(m) * itpan(m) * 2
     &            * itxnm(m) * itynm(m) * itznm(m)
          ! tr(np,ne,nx*ny*nz)
          endif


*-----------------------------------------------------------------------
*        t-deposit tally
*-----------------------------------------------------------------------
        else if( ital(m) .eq. 13 ) then
          if( itmsh(m) .eq. 1 ) then    ! tdepstreg
            mnmax = ( itenm(m) + 1 ) * itpan(m) * 2
     &            * itrgn(m) * ittnm(m)
            ! tr(np,ne,nr,nt)
          else if(itmsh(m) .eq. 2) then ! tdepstrz
            mnmax = ( itenm(m) + 1 ) * itpan(m) * 2
     &            * ittnm(m)
     &            * itrnm(m) * itznm(m)
            ! tr(np,ne,nt,nr,nz)
          else if(itmsh(m) .eq. 3) then ! tdepstxyz
            mnmax = ( itenm(m) + 1 ) * itpan(m) * 2
     &            * ittnm(m)
     &            * itxnm(m) * itynm(m) * itznm(m)
            ! tr(np,ne,nt,nx*ny*nz)
          else if( itmsh(m) .eq. 4 ) then ! tdepsttet
            mnmax = ( itenm(m) + 1 ) * itpan(m) * 2
     &            * itrgn(m) * ittnm(m)
            ! tr(np,ne,nr,nt)
          endif

*-----------------------------------------------------------------------
*        t-deposit2 tally
*-----------------------------------------------------------------------
        else if( ital(m) .eq. 14 ) then
          if( itmsh(m) .eq. 1 ) then    ! tdpst2reg
            mnmax = ( itenm(m) + 1 ) * ( itenm2(m) + 1 )
     &            * itpan(m) * ittnm(m) * 2
            ! tr(np,0:e1,0:e2,nt,2)
          endif

*-----------------------------------------------------------------------
*         [t-sed]
*-----------------------------------------------------------------------
        else if ( ital(m) .eq. 15 ) then
          if( itmsh(m) .eq. 1 ) then ! tsedreg
            mnmax = itenm(m) * itpan(m) * 2
     &            * itrgn(m)
            ! tr(np,ne,nr)
          else if( itmsh(m) .eq. 2 ) then ! tseddrz
            mnmax = itenm(m) * itpan(m) * 2
     &            * itrnm(m) * itznm(m)
            ! tr(np,ne,nr,nz)
          else if( itmsh(m) .eq. 3 ) then ! tsedxyz
            mnmax = itenm(m) * itpan(m) * 2
     &            * itxnm(m) * itynm(m) * itznm(m)
            ! tr(np,ne,nx*ny*nz)
          end if

*-----------------------------------------------------------------------
*         [t-point]
*-----------------------------------------------------------------------

        else if ( ital(m) .eq. 17 ) then
            mnmax = itenm(m) * itpan(m) * itmsh(m)
     &            * itmst(m) * ittnm(m) * mftal(m)
            ! tr(np,ne,nt,nr,nm)

*-----------------------------------------------------------------------
*         [t-wwg]
*-----------------------------------------------------------------------

        else if ( ital(m) .eq. 18 ) then
          if ( itmsh(m) .eq. 1 ) then
            mnmax = itenm(m) * itpan(m) * itrgn(m)
     &            * itmst(m) * ittnm(m) * mftal(m)
            ! tr(np,ne,nt,nr,nm)

          else if ( itmsh(m) .eq. 3 ) then
            mnmax = itenm(m) * itpan(m) * mftal(m)
     &            * itxnm(m) * itynm(m) * itznm(m) * itmst(m)
     &            * ittnm(m)
            ! tr(np,ne,nt,nx*ny*nz,nm)
          else if ( itmsh(m) .eq. 4 ) then ! twwgtet

            mnmax = itenm(m) * itpan(m) * itrgn(m)
     &            * itmst(m) * ittnm(m) * mftal(m)
            ! tr(np,ne,nt,nr,nm)
          endif

*-----------------------------------------------------------------------
*        [t-volume]
*-----------------------------------------------------------------------

        else if( ital(m) .eq. 21 ) then
            mnmax = itrgn(m) * 2
            ! tr(nr)

*-----------------------------------------------------------------------
*        [t-wwbg]
*-----------------------------------------------------------------------

        else if( ital(m) .eq. 22 ) then
          if ( itmsh(m) .eq. 1 ) then
            mnmax = itrgn(m) * 2 * 2
            ! tr(nr,2)
          else if ( itmsh(m) .eq. 3 ) then
            mnmax = itxnm(m) * itynm(m) * itznm(m) * 2 * 2
            ! tr(nx*ny*nz,2)
          else if ( itmsh(m) .eq. 4 ) then
            mnmax = itrgn(m) * 2 * 2
            ! tr(nr,2)
          endif

*-----------------------------------------------------------------------
        endif

      end subroutine

      subroutine CALC_TALSIZE_SUM(m,iax,mnmax,mnmax2)

        use partmod, only: itpan, itpat, jtpat ! frtati 2021/10/05

        integer,intent(in):: m,iax
        integer,intent(out):: mnmax

        include 'param.inc'
        common /tall00/ itnm, ital(itlmax), itals(itlmax), italm(itlmax)
        common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                  itout(itlmax), ittwo(itlmax)
        common /tall02/ itrgn(itlmax), itrgm(itlmax), itreg(itlmax)
        common /tall03/ itrty(itlmax), itrnm(itlmax), itrrg(itlmax),
     &                  rtrmi(itlmax), rtrma(itlmax), rtrdl(itlmax)
        common /tall04/ itxty(itlmax), itxnm(itlmax), itxrg(itlmax),
     &                  rtxmi(itlmax), rtxma(itlmax), rtxdl(itlmax)
        common /tall05/ ityty(itlmax), itynm(itlmax), ityrg(itlmax),
     &                  rtymi(itlmax), rtyma(itlmax), rtydl(itlmax)
        common /tall06/ itzty(itlmax), itznm(itlmax), itzrg(itlmax),
     &                  rtzmi(itlmax), rtzma(itlmax), rtzdl(itlmax)
        common /tall07/ itety(itlmax), itenm(itlmax), iterg(itlmax),
     &                  rtemi(itlmax), rtema(itlmax), rtedl(itlmax)
*        common /tall09/ itpan(itlmax), itpat(itlmax,6,2),
*    &                  jtpat(itlmax,6,6,2) ! frtati 2021/10/05
        common /tall10/ itaxn(itlmax), itaxs(itlmax,6)
        common /tall12/ itman(itlmax), itmat(itlmax), itmct(itlmax),
     &                  itnun(itlmax), itnuc(itlmax),
     &                  itndz(itlmax), itndn(itlmax),
     &                  itnkz(itlmax), itnkn(itlmax)
        common /tall13/ itrcn(itlmax), itrcr(itlmax), itrca(itlmax)
        common /tall17/ itndy(itlmax)
        common /tall20/ ittty(itlmax), ittnm(itlmax), ittrg(itlmax),
     &                  rttmi(itlmax), rttma(itlmax), rttdl(itlmax)
        common /tall23/ itaty(itlmax), itanm(itlmax), itarg(itlmax),
     &                  rtami(itlmax), rtama(itlmax), rtadl(itlmax)
        common /tall45/ itmlp(itlmax), itmln(itlmax,6), itmst(itlmax),
     &                  itmli(itlmax,6),rtmme(itlmax,6),itmnt(itlmax,6),
     &                  itmpn(itlmax,6), itmpt(itlmax,6,6,2)
        common /tall56/ itety2(itlmax), itenm2(itlmax), iterg2(itlmax),
     &                  rtemi2(itlmax), rtema2(itlmax), rtedl2(itlmax)
        common /tall72/ itactnm(itlmax),itactrg(itlmax),itactmax(itlmax)      ! S.Abe 2018/02/15
        common /tall83/ itnzn(itlmax), itndm(itlmax) ! frtati 2022/02/18

        common /tall91/ itextstat(itlmax), mftal(itlmax) ! S.H. extstat 2024.3.18

        common /tallsum/itpan_sum(itlmax,6),itenm_sum(itlmax,6),
     &                  ittnm_sum(itlmax,6),itrgn_sum(itlmax,6),
     &                  itmst_sum(itlmax,6),mftal_sum(itlmax,6),
     &                  itrnm_sum(itlmax,6),itanm_sum(itlmax,6),
     &                  itxnm_sum(itlmax,6),itynm_sum(itlmax,6),
     &                  itznm_sum(itlmax,6),itndy_sum(itlmax,6),
     &                  itndz_sum(itlmax,6),itndn_sum(itlmax,6),
     &                  itndm_sum(itlmax,6),itenm2_sum(itlmax,6),
     &                  itrcn_sum(itlmax,6),itactnm_sum(itlmax,6),
     &                  itmsh_sum(itlmax,6)

*-----------------------------------------------------------------------

      mnmax = 0
      mnmax2 = 0

*-----------------------------------------------------------------------
*         [g-show]
*-----------------------------------------------------------------------
        if ( ital(m) .eq. 9 ) then

*-----------------------------------------------------------------------
*         [r-show]
*-----------------------------------------------------------------------
        else if ( ital(m) .eq. 10 ) then

*-----------------------------------------------------------------------
*         [3d-show]
*-----------------------------------------------------------------------
        else if ( ital(m) .eq. 11 ) then


*-----------------------------------------------------------------------
*        [user defined tally]
*-----------------------------------------------------------------------
        else if ( ital(m) .eq. 20 ) then

*-----------------------------------------------------------------------
*         [t-volume]
*-----------------------------------------------------------------------
        else if ( ital(m) .eq. 21 ) then
          itrgn_sum(m,iax) = itrgn(m)

*-----------------------------------------------------------------------
*         other
*-----------------------------------------------------------------------
        else if(  ital(m) .ge. 1 .and. ital(m) .le. 22 ) then 

          itpan_sum(m,iax) = itpan(m)
          itenm_sum(m,iax) = itenm(m)
          ittnm_sum(m,iax) = ittnm(m)
          itrgn_sum(m,iax) = itrgn(m)
          itmst_sum(m,iax) = itmst(m)
*         mftal_sum(m,iax) = mftal(m)
          mftal_sum(m,iax) = 2
          itrnm_sum(m,iax) = itrnm(m)
          itanm_sum(m,iax) = itanm(m)
          itxnm_sum(m,iax) = itxnm(m)
          itynm_sum(m,iax) = itynm(m)
          itznm_sum(m,iax) = itznm(m)
          itndy_sum(m,iax) = itndy(m)
          itndz_sum(m,iax) = itndz(m)
          itndn_sum(m,iax) = itndn(m)
          itndm_sum(m,iax) = itndm(m)
          itenm2_sum(m,iax) = itenm2(m)
          itrcn_sum(m,iax) = itrcn(m)
          itactnm_sum(m,iax) = itactnm(m)
          itmsh_sum(m,iax) = itmsh(m)

        endif


*-----------------------------------------------------------------------
*        t-track tally
*-----------------------------------------------------------------------
* itaxs(m,iax)
*     1:eng , 2:reg,   3:x.     4:y,      5:z,
*     6:r,    7:xy yx, 8:yz zy. 9:zx xz, 10:rz zr,
*    11:t,   12:rad,  13: deg, 14:tet,   15 dchain

        if( ital(m) .eq. 1 ) then
          if ( itmsh(m) .eq. 1 ) then   ! ttracreg

            if(iax <= itaxn(m)) then
              iaxnum = 0
              select case(itaxs(m,iax))
              case (1)   ! eng
                itenm_sum(m,iax) = 1
                iaxnum = 1
              case (2)   ! reg
                itrgn_sum(m,iax) = 1
                iaxnum = 1
              case (11)  ! t
                ittnm_sum(m,iax)  = 1
                iaxnum = 1
              case default
                
              end select

              mnmax = itenm_sum(m,iax) * itpan_sum(m,iax)
     &              * itrgn_sum(m,iax) * itmst_sum(m,iax)
     &              * ittnm_sum(m,iax) * mftal_sum(m,iax)
              mnmax = mnmax * iaxnum
              ! sum of tr(np,ne,nt,nr,nm)
            endif

          else if ( itmsh(m) .eq. 2 ) then   ! ttracrz

            if(iax <= itaxn(m)) then
              iaxnum = 0
              select case(itaxs(m,iax))
              case (1)   ! eng
                itenm_sum(m,iax) = 1
                iaxnum = 1
              case (5)   ! z
                itznm_sum(m,iax) = 1
                iaxnum = 1
              case (6)   ! r
                itrnm_sum(m,iax) = 1
                iaxnum = 1
              case (11)  ! t
                ittnm_sum(m,iax)  = 1
                iaxnum = 1
              case (12)  ! rad
                itanm_sum(m,iax)  = 1
                iaxnum = 1
              case (13)  ! deg
                itanm_sum(m,iax)  = 1
                iaxnum = 1
              case default

              end select 

              mnmax = itenm_sum(m,iax) * itpan_sum(m,iax)
     &              * mftal_sum(m,iax) * itrnm_sum(m,iax)
     &              * itznm_sum(m,iax) * itanm_sum(m,iax)
     &              * itmst_sum(m,iax) * ittnm_sum(m,iax)
              mnmax = mnmax * iaxnum
            ! tr(np,ne,nt,nr,nz,na,nm)
            endif

          else if ( itmsh(m) .eq. 3 ) then   ! ttracxyz

            if(iax <= itaxn(m)) then
              iaxnum = 0
              select case(itaxs(m,iax))
              case (1)   ! eng
                itenm_sum(m,iax) = 1
                iaxnum = 1
              case (3)   ! x
                itxnm_sum(m,iax) = 1
                iaxnum = 1
              case (4)   ! y
                itynm_sum(m,iax) = 1
                iaxnum = 1
              case (5)   ! z
                itznm_sum(m,iax) = 1
                iaxnum = 1
              case (11)  ! t
                ittnm_sum(m,iax)  = 1
                iaxnum = 1
              case default

              end select 

              mnmax = itenm_sum(m,iax) * itpan_sum(m,iax)
     &              * mftal_sum(m,iax) * itxnm_sum(m,iax)
     &              * itynm_sum(m,iax) * itznm_sum(m,iax)
     &              * itmst_sum(m,iax) * ittnm_sum(m,iax)
              mnmax = mnmax * iaxnum
              ! tr(np,ne,nt,nx*ny*nz,nm)
            endif

          else if ( itmsh(m) .eq. 4 ) then ! ttractet

            if(iax <= itaxn(m)) then
              iaxnum = 0
              select case(itaxs(m,iax))
              case (1)   ! eng
                itenm_sum(m,iax) = 1
                iaxnum = 1
              case (2)   ! reg
                itrgn_sum(m,iax) = 1
                iaxnum = 1
              case (11)  ! t
                ittnm_sum(m,iax)  = 1
                iaxnum = 1
              case default
                
              end select

              mnmax = itenm_sum(m,iax) * itpan_sum(m,iax)
     &              * itrgn_sum(m,iax) * itmst_sum(m,iax)
     &              * ittnm_sum(m,iax) * mftal_sum(m,iax)
              mnmax = mnmax * iaxnum
              ! sum of tr(np,ne,nt,nr,nm)

            endif

          endif

*-----------------------------------------------------------------------
*        t-adjoint tally
*-----------------------------------------------------------------------
*     1:eng , 2:reg,   3:x.     4:y,      5:z,
*     6:r,    7:xy yx, 8:yz zy. 9:zx xz, 10:rz zr,
*    11:t,   12:rad,  13: deg

        else if( ital(m) .eq. 19 ) then
          if ( itmsh(m) .eq. 1 ) then   ! tadjreg

            if(iax <= itaxn(m)) then
              iaxnum = 0
              select case(itaxs(m,iax))
              case (1)   ! eng
                itenm_sum(m,iax) = 1
                iaxnum = 1
              case (2)   ! reg
                itrgn_sum(m,iax) = 1
                iaxnum = 1
              case (11)  ! t
                ittnm_sum(m,iax)  = 1
                iaxnum = 1
              case default
                
              end select

              mnmax = itenm_sum(m,iax) * itpan_sum(m,iax)
     &              * itrgn_sum(m,iax) * itmst_sum(m,iax)
     &              * ittnm_sum(m,iax) * mftal_sum(m,iax)
              mnmax = mnmax * iaxnum
              ! sum of tr(np,ne,nt,nr,nm)

            endif

          else if ( itmsh(m) .eq. 2 ) then   ! tadjrz

            if(iax <= itaxn(m)) then
              iaxnum = 0
              select case(itaxs(m,iax))
              case (1)   ! eng
                itenm_sum(m,iax) = 1
                iaxnum = 1
              case (5)   ! z
                itznm_sum(m,iax) = 1
                iaxnum = 1
              case (6)   ! r
                itrnm_sum(m,iax) = 1
                iaxnum = 1
              case (11)  ! t
                ittnm_sum(m,iax)  = 1
                iaxnum = 1
              case (12)  ! rad
                itanm_sum(m,iax)  = 1
                iaxnum = 1
              case (13)  ! deg
                itanm_sum(m,iax)  = 1
                iaxnum = 1
              case default

              end select 

              mnmax = itenm_sum(m,iax) * itpan_sum(m,iax)
     &              * mftal_sum(m,iax) * itrnm_sum(m,iax)
     &              * itznm_sum(m,iax) * itanm_sum(m,iax)
     &              * itmst_sum(m,iax) * ittnm_sum(m,iax)
              mnmax = mnmax * iaxnum
            ! tr(np,ne,nt,nr,nz,na,nm)
            endif

          else if ( itmsh(m) .eq. 3 ) then   ! tadjxyz

            if(iax <= itaxn(m)) then
              iaxnum = 0
              select case(itaxs(m,iax))
              case (1)   ! eng
                itenm_sum(m,iax) = 1
                iaxnum = 1
              case (3)   ! x
                itxnm_sum(m,iax) = 1
                iaxnum = 1
              case (4)   ! y
                itynm_sum(m,iax) = 1
                iaxnum = 1
              case (5)   ! z
                itznm_sum(m,iax) = 1
                iaxnum = 1
              case (11)  ! t
                ittnm_sum(m,iax)  = 1
                iaxnum = 1
              case default

              end select 

              mnmax = itenm_sum(m,iax) * itpan_sum(m,iax)
     &              * mftal_sum(m,iax) * itxnm_sum(m,iax)
     &              * itynm_sum(m,iax) * itznm_sum(m,iax)
     &              * itmst_sum(m,iax) * ittnm_sum(m,iax)
              mnmax = mnmax * iaxnum
              ! tr(np,ne,nt,nx*ny*nz,nm)
            endif

          endif

*-----------------------------------------------------------------------
*        t-cross tally
*-----------------------------------------------------------------------
*     1:eng ,  2:reg,    3:x.     4:y,    5:z,
*     6:r,     7:xy yx,  8:cos,   9:t,   10:the
*    11:yz zy.12:zx xz, 13:rz zr,14:let

        else if( ital(m) .eq. 2 ) then
          if ( itmsh(m) .eq. 1 ) then   ! tcrsreg

            if(iax <= itaxn(m)) then
              iaxnum = 0
              select case(itaxs(m,iax))
              case (1)   ! eng
                itenm_sum(m,iax) = 1
                iaxnum = 1
              case (2)   ! r
                itrcn_sum(m,iax) = 1
                iaxnum = 1
              case (8)  ! cos
                itanm_sum(m,iax)  = 1
                iaxnum = 1
              case (9)  ! t
                ittnm_sum(m,iax)  = 1
                iaxnum = 1
              case (10)  ! the
                itanm_sum(m,iax)  = 1
                iaxnum = 1
              case (14)   ! let
                itenm_sum(m,iax) = 1
                iaxnum = 1
              case default
                
              end select

              mnmax = itenm_sum(m,iax) * itpan_sum(m,iax)
     &              * mftal_sum(m,iax) * itrcn_sum(m,iax)
     &              * itanm_sum(m,iax) * itmst_sum(m,iax)
     &              * ittnm_sum(m,iax)
              mnmax = mnmax * iaxnum
            ! tr(np,ne,na,nt,nr,nm)

            endif

          else if ( itmsh(m) .eq. 2 ) then   ! tcrsrz

            if(iax <= itaxn(m)) then
              iaxnum = 0
              select case(itaxs(m,iax))
              case (1)   ! eng
                itenm_sum(m,iax) = 1
                iaxnum = 1
              case (5)   ! z
                itznm_sum(m,iax) = 1
                iaxnum = 1
              case (6)   ! r
                itrnm_sum(m,iax) = 1
                iaxnum = 1
              case (8)  ! cos
                itanm_sum(m,iax)  = 1
                iaxnum = 1
              case (9)  ! t
                ittnm_sum(m,iax)  = 1
                iaxnum = 1
              case (10)  ! the
                itanm_sum(m,iax)  = 1
                iaxnum = 1
              case (14)   ! let
                itenm_sum(m,iax) = 1
                iaxnum = 1
              case default

              end select 

            mnmax = itenm_sum(m,iax) * itpan_sum(m,iax)
     &            * mftal_sum(m,iax) * ( itrnm_sum(m,iax) + 1 )
     &            * itznm_sum(m,iax) * itanm_sum(m,iax)
     &            * itmst_sum(m,iax) * ittnm_sum(m,iax)
     &            + itenm_sum(m,iax) * itpan_sum(m,iax)
     &            * mftal_sum(m,iax) * ( itznm_sum(m,iax) + 1 )
     &            * itrnm_sum(m,iax) * itanm_sum(m,iax)
     &            * itmst_sum(m,iax) * ittnm_sum(m,iax)
            mnmax2 = itenm_sum(m,iax) * itpan_sum(m,iax)
     &            * mftal_sum(m,iax) * ( itrnm_sum(m,iax) + 1 )
     &            * itznm_sum(m,iax) * itanm_sum(m,iax)
     &            * itmst_sum(m,iax) * ittnm_sum(m,iax)
              mnmax = mnmax * iaxnum
              mnmax2 = mnmax2 * iaxnum
            ! tr(np,ne,na,nt,nr+1,nz,nm) + tz(np,ne,na,nt,nz+1,nr,nm)
            endif


          else if ( itmsh(m) .eq. 3 ) then   ! tcrsxyz

            if(iax <= itaxn(m)) then
              iaxnum = 0
              select case(itaxs(m,iax))
              case (1)   ! eng
                itenm_sum(m,iax) = 1
                iaxnum = 1
              case (3)   ! x
                itxnm_sum(m,iax) = 1
                iaxnum = 1
              case (4)   ! y
                itynm_sum(m,iax) = 1
                iaxnum = 1
              case (5)   ! z
                itznm_sum(m,iax) = 1
                iaxnum = 1
              case (8)  ! cos
                itanm_sum(m,iax)  = 1
                iaxnum = 1
              case (9)  ! t
                ittnm_sum(m,iax)  = 1
                iaxnum = 1
              case (10)  ! the
                itanm_sum(m,iax)  = 1
                iaxnum = 1
              case (14)   ! let
                itenm_sum(m,iax) = 1
                iaxnum = 1
              case default

              end select 

            mnmax = itenm_sum(m,iax) * itpan_sum(m,iax)
     &            * mftal_sum(m,iax) * itxnm_sum(m,iax)
     &            * itynm_sum(m,iax) * ( itznm_sum(m,iax) + 1 )
     &            * itanm_sum(m,iax) * itmst_sum(m,iax)
     &            * ittnm_sum(m,iax)
              mnmax = mnmax * iaxnum
            ! tr(np,ne,na,nt,nx,ny,nz+1,nm)
            endif

          endif

*-----------------------------------------------------------------------
*         [t-yield]
*-----------------------------------------------------------------------
*     1:mass,   2:reg,    3:x.     4:y,     5:z,
*     6:r,      7:charge, 8:chart, 9:xy yx, 10:yz zy.
*    11:zx xz, 12:rz zr, 13:dchain,14:tet

        else if ( ital(m) .eq. 3 ) then
          if ( itmsh(m) .eq. 1 ) then   ! tyildreg

            if(iax <= itaxn(m)) then
              iaxnum = 0
              select case(itaxs(m,iax))
              case (2)   ! reg
                itrgn_sum(m,iax) = 1
                iaxnum = 1
              case default
                
              end select

              mnmax = itndz_sum(m,iax) * itndn_sum(m,iax)
     &              * mftal_sum(m,iax) * itrgn_sum(m,iax)
     &              * (itndm_sum(m,iax)+1)
              mnmax = mnmax * iaxnum

            endif

          else if ( itmsh(m) .eq. 2 ) then   ! tyildrz

            if(iax <= itaxn(m)) then
              iaxnum = 0
              select case(itaxs(m,iax))
              case (5)   ! z
                itznm_sum(m,iax) = 1
                iaxnum = 1
              case (6)   ! r
                itrnm_sum(m,iax) = 1
                iaxnum = 1
              case default
                
              end select

              mnmax = itndz_sum(m,iax) * itndn_sum(m,iax)
     &              * mftal_sum(m,iax) * itrnm_sum(m,iax)
     &              * itznm_sum(m,iax) * (itndm_sum(m,iax)+1)
              mnmax = mnmax * iaxnum

            endif


          else if ( itmsh(m) .eq. 3 ) then   ! tyildxyz

            if(iax <= itaxn(m)) then
              iaxnum = 0
              select case(itaxs(m,iax))
              case (3)   ! x
                itxnm_sum(m,iax) = 1
                iaxnum = 1
              case (4)   ! y
                itynm_sum(m,iax) = 1
                iaxnum = 1
              case (5)   ! z
                itznm_sum(m,iax) = 1
                iaxnum = 1
              case default
                
              end select

              mnmax = itndz_sum(m,iax) * itndn_sum(m,iax)
     &              * mftal_sum(m,iax) * itxnm_sum(m,iax)
     &              * itynm_sum(m,iax) * itznm_sum(m,iax)
     &              * (itndm_sum(m,iax)+1)
              mnmax = mnmax * iaxnum

            endif
          else if ( itmsh(m) .eq. 4 ) then   ! tyildtet

            if(iax <= itaxn(m)) then
              iaxnum = 0
              select case(itaxs(m,iax))
              case (2)   ! reg
                itrgn_sum(m,iax) = 1
                iaxnum = 1
              case default
                
              end select

              mnmax = itndz_sum(m,iax) * itndn_sum(m,iax)
     &              * mftal_sum(m,iax) * itrgn_sum(m,iax)
     &              * (itndm_sum(m,iax)+1)
              mnmax = mnmax * iaxnum

            endif
          endif

*-----------------------------------------------------------------------
*        t-heat tally
*-----------------------------------------------------------------------
*     1:eng,    2:reg,    3:x.     4:y,     5:z,
*     6:r,      7:xy yx,  8:yz zy. 9:zx xz, 10:rz zr
        else if( ital(m) .eq. 4 ) then
          if( itmsh(m) .eq. 1 ) then    ! thetreg

            if(iax <= itaxn(m)) then
              iaxnum = 0
              select case(itaxs(m,iax))
              case (1)   ! eng
                itenm_sum(m,iax) = 1
                iaxnum = 1
              case (2)   ! reg
                itrgn_sum(m,iax) = 1
                iaxnum = 1
              case default
                
              end select

              mnmax = itndy_sum(m,iax) * mftal_sum(m,iax)
     &              * ( itenm_sum(m,iax) + 1 ) * itrgn_sum(m,iax) + 5
              mnmax = mnmax * iaxnum
            ! sum of tr(nd,nr)
            endif

          else if(itmsh(m) .eq. 2) then ! thetrz

            if(iax <= itaxn(m)) then
              iaxnum = 0
              select case(itaxs(m,iax))
              case (1)   ! eng
                itenm_sum(m,iax) = 1
                iaxnum = 1
              case (5)   ! z
                itznm_sum(m,iax) = 1
                iaxnum = 1
              case (6)   ! r
                itrnm_sum(m,iax) = 1
                iaxnum = 1
              case default
                
              end select

              mnmax = itndy_sum(m,iax) * mftal_sum(m,iax)
     &              * ( itenm_sum(m,iax) + 1 ) * itrnm_sum(m,iax)
     &              * itznm_sum(m,iax) + 5
              mnmax = mnmax * iaxnum
            ! sum of tr(nd,nr,nz)
            endif

          else if(itmsh(m) .eq. 3) then ! thetxyz

            if(iax <= itaxn(m)) then
              iaxnum = 0
              select case(itaxs(m,iax))
              case (1)   ! eng
                itenm_sum(m,iax) = 1
                iaxnum = 1
              case (3)   ! x
                itxnm_sum(m,iax) = 1
                iaxnum = 1
              case (4)   ! y
                itynm_sum(m,iax) = 1
                iaxnum = 1
              case (5)   ! z
                itznm_sum(m,iax) = 1
                iaxnum = 1
              case default
                
              end select

              mnmax = itndy_sum(m,iax) * mftal_sum(m,iax)
     &              * ( itenm_sum(m,iax) + 1 ) * itxnm_sum(m,iax)
     &              * itynm_sum(m,iax) * itznm_sum(m,iax) + 5
              mnmax = mnmax * iaxnum
            ! sum of tr(nd,nx,ny,nz)
            endif

          endif

*-----------------------------------------------------------------------
*         [t-star]
*-----------------------------------------------------------------------
*     1:eng,    2:reg,    3:x.     4:y,     5:z,
*     6:r,      7:xy yx,  8:yz zy. 9:zx xz, 10:rz zr
*    11:t      12:act
        else if ( ital(m) .eq. 5 ) then
          if( itmsh(m) .eq. 1 ) then ! tstarreg

            if(iax <= itaxn(m)) then
              iaxnum = 0
              select case(itaxs(m,iax))
              case (1)   ! eng
                itenm_sum(m,iax) = 1
                iaxnum = 1
              case (2)   ! reg
                itrgn_sum(m,iax) = 1
                iaxnum = 1
              case (11)   ! t
                ittnm_sum(m,iax) = 1
                iaxnum = 1
              case (12)   ! act
                itactnm_sum(m,iax) = 1
                iaxnum = 1
              case default
                
              end select

              mnmax = itenm_sum(m,iax) * ittnm_sum(m,iax)
     &              * itpan_sum(m,iax) * mftal_sum(m,iax)
     &              * itrgn_sum(m,iax) * itactnm_sum(m,iax)
              mnmax = mnmax * iaxnum
          ! sum of tr(np,ne,nt,nr,2)
            endif

          else if( itmsh(m) .eq. 2 ) then ! tstarrz

            if(iax <= itaxn(m)) then
              iaxnum = 0
              select case(itaxs(m,iax))
              case (1)   ! eng
                itenm_sum(m,iax) = 1
                iaxnum = 1
              case (5)   ! z
                itznm_sum(m,iax) = 1
                iaxnum = 1
              case (6)   ! r
                itrnm_sum(m,iax) = 1
                iaxnum = 1
              case (11)   ! t
                ittnm_sum(m,iax) = 1
                iaxnum = 1
              case default
                
              end select

              mnmax = itenm_sum(m,iax) * itpan_sum(m,iax)
     &              * mftal_sum(m,iax) * ittnm_sum(m,iax)
     &              * itrnm_sum(m,iax) * itznm_sum(m,iax)
              mnmax = mnmax * iaxnum
          ! sum of tr(np,ne,nt,nr,nz,2)
            endif

          else if( itmsh(m) .eq. 3 ) then ! tstarxyz

            if(iax <= itaxn(m)) then
              iaxnum = 0
              select case(itaxs(m,iax))
              case (1)   ! eng
                itenm_sum(m,iax) = 1
                iaxnum = 1
              case (3)   ! x
                itxnm_sum(m,iax) = 1
                iaxnum = 1
              case (4)   ! y
                itynm_sum(m,iax) = 1
                iaxnum = 1
              case (5)   ! z
                itznm_sum(m,iax) = 1
                iaxnum = 1
              case (11)   ! t
                ittnm_sum(m,iax) = 1
                iaxnum = 1
              case default
                
              end select

              mnmax = itenm_sum(m,iax) * itpan_sum(m,iax)
     &              * mftal_sum(m,iax) * ittnm_sum(m,iax)
     &              * itxnm_sum(m,iax) * itynm_sum(m,iax)
     &              * itznm_sum(m,iax)
              mnmax = mnmax * iaxnum
          ! sum of tr(np,ne,nt,nx,ny,nz,2)
            endif
          endif

*-----------------------------------------------------------------------
*         [t-time]
*-----------------------------------------------------------------------
*     1:eng,    2:reg,    3:x.     4:y,     5:z,
*     6:r,      7:xy yx,  8:yz zy. 9:zx xz, 10:rz zr
*    11:t
        else if ( ital(m) .eq. 6 ) then
          if( itmsh(m) .eq. 1 ) then ! ttimereg

            if(iax <= itaxn(m)) then
              iaxnum = 0
              select case(itaxs(m,iax))
              case (1)   ! eng
                itenm_sum(m,iax) = 1
                iaxnum = 1
              case (2)   ! reg
                itrgn_sum(m,iax) = 1
                iaxnum = 1
              case (11)   ! t
                ittnm_sum(m,iax) = 1
                iaxnum = 1
              case default
                
              end select

              mnmax = itenm_sum(m,iax) * itpan_sum(m,iax)
     &              * mftal_sum(m,iax) * ittnm_sum(m,iax)
     &              * itrgn_sum(m,iax)
              mnmax = mnmax * iaxnum
          ! sum of tr(np,ne,nt,nr,2)
            endif

          else if( itmsh(m) .eq. 2 ) then ! ttimedrz

            if(iax <= itaxn(m)) then
              iaxnum = 0
              select case(itaxs(m,iax))
              case (1)   ! eng
                itenm_sum(m,iax) = 1
                iaxnum = 1
              case (5)   ! z
                itznm_sum(m,iax) = 1
                iaxnum = 1
              case (6)   ! r
                itrnm_sum(m,iax) = 1
                iaxnum = 1
              case (11)   ! t
                ittnm_sum(m,iax) = 1
                iaxnum = 1
              case default
                
              end select

              mnmax = itenm_sum(m,iax) * itpan_sum(m,iax)
     &              * mftal_sum(m,iax) * ittnm_sum(m,iax)
     &              * itrnm_sum(m,iax) * itznm_sum(m,iax)
              mnmax = mnmax * iaxnum
          ! sum of tr(np,ne,nt,nr,nz,2)
            endif

          else if( itmsh(m) .eq. 3 ) then ! ttimexyz

            if(iax <= itaxn(m)) then
              iaxnum = 0
              select case(itaxs(m,iax))
              case (1)   ! eng
                itenm_sum(m,iax) = 1
                iaxnum = 1
              case (3)   ! x
                itxnm_sum(m,iax) = 1
                iaxnum = 1
              case (4)   ! y
                itynm_sum(m,iax) = 1
                iaxnum = 1
              case (5)   ! z
                itznm_sum(m,iax) = 1
                iaxnum = 1
              case (11)   ! t
                ittnm_sum(m,iax) = 1
                iaxnum = 1
              case default
                
              end select

              mnmax = itenm_sum(m,iax) * itpan_sum(m,iax)
     &              * mftal_sum(m,iax) * ittnm_sum(m,iax)
     &              * itxnm_sum(m,iax) * itynm_sum(m,iax)
     &              * itznm_sum(m,iax)
              mnmax = mnmax * iaxnum
          ! sum of tr(np,ne,nt,nx,ny,nz,2)
            endif
          endif

*-----------------------------------------------------------------------
*         [t-dpa]
*-----------------------------------------------------------------------
*     1:reg,    2:x.     3:y,     4:z,     5:r,
*     6:xy yx,  7:yz zy. 8:zx xz, 9:rz zr, 10:tet
        else if ( ital(m) .eq. 7 ) then
          if( itmsh(m) .eq. 1 ) then ! tdpareg

            if(iax <= itaxn(m)) then
              iaxnum = 0
              select case(itaxs(m,iax))
              case (1)   ! reg
                itrgn_sum(m,iax) = 1
                iaxnum = 1
              case default
                
              end select

              mnmax = itndy_sum(m,iax) * itpan_sum(m,iax)
     &              * mftal_sum(m,iax) * itrgn_sum(m,iax)
              mnmax = mnmax * iaxnum
          ! sum of tr(nd,np,nr,2)
            endif

          else if( itmsh(m) .eq. 2 ) then ! tdpadrz

            if(iax <= itaxn(m)) then
              iaxnum = 0
              select case(itaxs(m,iax))
              case (4)   ! z
                itznm_sum(m,iax) = 1
                iaxnum = 1
              case (5)   ! r
                itrnm_sum(m,iax) = 1
                iaxnum = 1
              case default
                
              end select

            mnmax = itndy_sum(m,iax) * itpan_sum(m,iax)
     &            * mftal_sum(m,iax) * itrnm_sum(m,iax)
     &            * itznm_sum(m,iax)
              mnmax = mnmax * iaxnum
          ! sum of tr(nd,np,nr,nz,2)
            endif

          else if( itmsh(m) .eq. 3 ) then ! tdpaxyz

            if(iax <= itaxn(m)) then
              iaxnum = 0
              select case(itaxs(m,iax))
              case (4)   ! z
                itznm_sum(m,iax) = 1
                iaxnum = 1
              case (5)   ! r
                itrnm_sum(m,iax) = 1
                iaxnum = 1
              case default
                
              end select

              mnmax = itndy_sum(m,iax) * itpan_sum(m,iax)
     &              * mftal_sum(m,iax) * itxnm_sum(m,iax)
     &              * itynm_sum(m,iax) * itznm_sum(m,iax)
              mnmax = mnmax * iaxnum
          ! sum of tr(nd,np,nx,ny,nz,2)
            endif
          else if( itmsh(m) .eq. 4 ) then ! tdpatet

            if(iax <= itaxn(m)) then
              iaxnum = 0
              select case(itaxs(m,iax))
              case (1)   ! reg
                itrgn_sum(m,iax) = 1
                iaxnum = 1
              case default
                
              end select

              mnmax = itndy_sum(m,iax) * itpan_sum(m,iax)
     &              * mftal_sum(m,iax) * itrgn_sum(m,iax)
              mnmax = mnmax * iaxnum
          ! sum of tr(nd,np,nr,2)
            endif
          endif

*-----------------------------------------------------------------------
*         [t-product]
*-----------------------------------------------------------------------
*     1:eng,  2:reg,    3:x.     4:y,      5:z,
*     6:r,    7:xy yx,  8:yz zy. 9:zx xz, 10:rz zr,
*    11:t,   12:cos,   13:the,   14:tet,  15:let

        else if ( ital(m) .eq. 8 ) then
          if( itmsh(m) .eq. 1 ) then ! tprodreg

            if(iax <= itaxn(m)) then
              iaxnum = 0
              select case(itaxs(m,iax))
              case (1)   ! eng
                itenm_sum(m,iax) = 1
                iaxnum = 1
              case (2)   ! reg
                itrgn_sum(m,iax) = 1
                iaxnum = 1
              case (11)   ! t
                ittnm_sum(m,iax) = 1
                iaxnum = 1
              case (12)   ! cos
                itanm_sum(m,iax) = 1
                iaxnum = 1
              case (13)   ! the
                itanm_sum(m,iax) = 1
                iaxnum = 1
              case (15)   ! let
                itenm_sum(m,iax) = 1
                iaxnum = 1
              case default
                
              end select

              mnmax = itenm_sum(m,iax) * ittnm_sum(m,iax)
     &              * itpan_sum(m,iax) * mftal_sum(m,iax)
     &              * itrgn_sum(m,iax) * itanm_sum(m,iax)
              mnmax = mnmax * iaxnum
          ! sum of tr(np,ne,nt,na,nr,2)
            endif

          else if( itmsh(m) .eq. 2 ) then ! tproddrz

            if(iax <= itaxn(m)) then
              iaxnum = 0
              select case(itaxs(m,iax))
              case (1)   ! eng
                itenm_sum(m,iax) = 1
                iaxnum = 1
              case (5)   ! z
                itznm_sum(m,iax) = 1
                iaxnum = 1
              case (6)   ! r
                itrnm_sum(m,iax) = 1
                iaxnum = 1
              case (11)   ! t
                ittnm_sum(m,iax) = 1
                iaxnum = 1
              case (12)   ! cos
                itanm_sum(m,iax) = 1
                iaxnum = 1
              case (13)   ! the
                itanm_sum(m,iax) = 1
                iaxnum = 1
              case (15)   ! let
                itenm_sum(m,iax) = 1
                iaxnum = 1
              case default
                
              end select

              mnmax = itenm_sum(m,iax) * itpan_sum(m,iax)
     &              * mftal_sum(m,iax) * ittnm_sum(m,iax)
     &              * itanm_sum(m,iax) * itrnm_sum(m,iax)
     &              * itznm_sum(m,iax)
              mnmax = mnmax * iaxnum
          ! tr(sum of np,ne,nt,na,nr,nz,2)
            endif


          else if( itmsh(m) .eq. 3 ) then ! tprodxyz

            if(iax <= itaxn(m)) then
              iaxnum = 0
              select case(itaxs(m,iax))
              case (1)   ! eng
                itenm_sum(m,iax) = 1
                iaxnum = 1
              case (3)   ! x
                itxnm_sum(m,iax) = 1
                iaxnum = 1
              case (4)   ! y
                itynm_sum(m,iax) = 1
                iaxnum = 1
              case (5)   ! z
                itznm_sum(m,iax) = 1
                iaxnum = 1
              case (11)   ! t
                ittnm_sum(m,iax) = 1
                iaxnum = 1
              case (12)   ! cos
                itanm_sum(m,iax) = 1
                iaxnum = 1
              case (13)   ! the
                itanm_sum(m,iax) = 1
                iaxnum = 1
              case (15)   ! let
                itenm_sum(m,iax) = 1
                iaxnum = 1
              case default
                
              end select

              mnmax = itenm_sum(m,iax) * itpan_sum(m,iax)
     &              * mftal_sum(m,iax) * ittnm_sum(m,iax)
     &              * itanm_sum(m,iax) * itxnm_sum(m,iax)
     &              * itynm_sum(m,iax) * itznm_sum(m,iax)
              mnmax = mnmax * iaxnum
          ! sum of tr(sum of np,ne,nt,na,nx*ny*nz,2)
            endif


          else if( itmsh(m) .eq. 4 ) then ! tprodtet

            if(iax <= itaxn(m)) then
              iaxnum = 0
              select case(itaxs(m,iax))
              case (1)   ! eng
                itenm_sum(m,iax) = 1
                iaxnum = 1
              case (2)   ! reg
                itrgn_sum(m,iax) = 1
                iaxnum = 1
              case (11)   ! t
                ittnm_sum(m,iax) = 1
                iaxnum = 1
              case (12)   ! cos
                itanm_sum(m,iax) = 1
                iaxnum = 1
              case (13)   ! the
                itanm_sum(m,iax) = 1
                iaxnum = 1
              case (15)   ! let
                itenm_sum(m,iax) = 1
                iaxnum = 1
              case default
                
              end select

              mnmax = itenm_sum(m,iax) * ittnm_sum(m,iax)
     &              * itpan_sum(m,iax) * mftal_sum(m,iax)
     &              * itrgn_sum(m,iax) * itanm_sum(m,iax)
              mnmax = mnmax * iaxnum
          ! sum of tr(np,ne,nt,na,nr,2)
            endif
          endif

*-----------------------------------------------------------------------
*         [g-show]
*-----------------------------------------------------------------------
        else if ( ital(m) .eq. 9 ) then

*-----------------------------------------------------------------------
*         [r-show]
*-----------------------------------------------------------------------
        else if ( ital(m) .eq. 10 ) then

*-----------------------------------------------------------------------
*         [3d-show]
*-----------------------------------------------------------------------
        else if ( ital(m) .eq. 11 ) then

*-----------------------------------------------------------------------
*         [t-let]
*-----------------------------------------------------------------------
*     1:let,  2:reg,    3:x.     4:y,      5:z,
*     6:r,    7:xy yx,  8:yz zy. 9:zx xz, 10:rz zr,

        else if ( ital(m) .eq. 12 ) then
          if( itmsh(m) .eq. 1 ) then ! tletreg

            if(iax <= itaxn(m)) then
              iaxnum = 0
              select case(itaxs(m,iax))
              case (1)   ! let
                itenm_sum(m,iax) = 1
                iaxnum = 1
              case (2)   ! reg
                itrgn_sum(m,iax) = 1
                iaxnum = 1
              case default
                
              end select

              mnmax = itenm_sum(m,iax) * itpan_sum(m,iax)
     &              * mftal_sum(m,iax) * itrgn_sum(m,iax)
              mnmax = mnmax * iaxnum
          ! sum of tr(np,ne,nr)
            endif

          else if( itmsh(m) .eq. 2 ) then ! tletdrz

            if(iax <= itaxn(m)) then
              iaxnum = 0
              select case(itaxs(m,iax))
              case (1)   ! let
                itenm_sum(m,iax) = 1
                iaxnum = 1
              case (5)   ! z
                itznm_sum(m,iax) = 1
                iaxnum = 1
              case (6)   ! r
                itrnm_sum(m,iax) = 1
                iaxnum = 1
              case default
                
              end select

              mnmax = itenm_sum(m,iax) * itpan_sum(m,iax)
     &              * mftal_sum(m,iax) * itrnm_sum(m,iax)
     &              * itznm_sum(m,iax)
              mnmax = mnmax * iaxnum
          ! sum of tr(np,ne,nr,nz)
            endif

          else if( itmsh(m) .eq. 3 ) then ! tletxyz

            if(iax <= itaxn(m)) then
              iaxnum = 0
              select case(itaxs(m,iax))
              case (1)   ! let
                itenm_sum(m,iax) = 1
                iaxnum = 1
              case (3)   ! x
                itxnm_sum(m,iax) = 1
                iaxnum = 1
              case (4)   ! y
                itynm_sum(m,iax) = 1
                iaxnum = 1
              case (5)   ! z
                itznm_sum(m,iax) = 1
                iaxnum = 1
              case default
                
              end select

              mnmax = itenm_sum(m,iax) * itpan_sum(m,iax)
     &              * mftal_sum(m,iax) * itxnm_sum(m,iax)
     &              * itynm_sum(m,iax) * itznm_sum(m,iax)
              mnmax = mnmax * iaxnum
          ! sum of tr(np,ne,nx*ny*nz)
            endif
          endif


*-----------------------------------------------------------------------
*        t-deposit tally
*-----------------------------------------------------------------------
*     1:eng,  2:reg,    3:x.      4:y,      5:z,
*     6:r,    7:xy yx,  8:yz zy.  9:zx xz, 10:rz zr,
*    11:t,   12:t-eng, 13;eng-t, 14:tet

        else if( ital(m) .eq. 13 ) then
          if( itmsh(m) .eq. 1 ) then    ! tdepstreg

            if(iax <= itaxn(m)) then
              iaxnum = 0
              select case(itaxs(m,iax))
              case (1)   ! eng
                itenm_sum(m,iax) = 1
                iaxnum = 1
              case (2)   ! reg
                itrgn_sum(m,iax) = 1
                iaxnum = 1
              case (11)   ! t
                ittnm_sum(m,iax) = 1
                iaxnum = 1
              case (12)   ! t-eng

              case (13)   ! eng-t

              case default
                
              end select

              mnmax = ( itenm_sum(m,iax) + 1 ) * itpan_sum(m,iax)
     &              * mftal_sum(m,iax) * itrgn_sum(m,iax)
     &              * ittnm_sum(m,iax)
              mnmax = mnmax * iaxnum
            ! sum of tr(np,ne,nr,nt)
            endif

          else if(itmsh(m) .eq. 2) then ! tdepstrz

            if(iax <= itaxn(m)) then
              iaxnum = 0
              select case(itaxs(m,iax))
              case (1)   ! eng
                itenm_sum(m,iax) = 1
                iaxnum = 1
              case (5)   ! z
                itznm_sum(m,iax) = 1
                iaxnum = 1
              case (6)   ! r
                itrnm_sum(m,iax) = 1
                iaxnum = 1
              case (11)   ! t
                ittnm_sum(m,iax) = 1
                iaxnum = 1
              case (12)   ! t-eng

              case (13)   ! eng-t

              case default
                
              end select

              mnmax = ( itenm_sum(m,iax) + 1 ) * itpan_sum(m,iax)
     &              * mftal_sum(m,iax) * ittnm_sum(m,iax)
     &              * itrnm_sum(m,iax) * itznm_sum(m,iax)
              mnmax = mnmax * iaxnum
            ! sum of tr(np,ne,nt,nr,nz)
            endif

          else if(itmsh(m) .eq. 3) then ! tdepstxyz

            if(iax <= itaxn(m)) then
              iaxnum = 0
              select case(itaxs(m,iax))
              case (1)   ! eng
                itenm_sum(m,iax) = 1
                iaxnum = 1
              case (3)   ! x
                itxnm_sum(m,iax) = 1
                iaxnum = 1
              case (4)   ! y
                itynm_sum(m,iax) = 1
                iaxnum = 1
              case (5)   ! z
                itznm_sum(m,iax) = 1
                iaxnum = 1
              case (11)   ! t
                ittnm_sum(m,iax) = 1
                iaxnum = 1
              case (12)   ! t-eng

              case (13)   ! eng-t

              case default
                
              end select

              mnmax = ( itenm_sum(m,iax) + 1 ) * itpan_sum(m,iax)
     &              * mftal_sum(m,iax) * ittnm_sum(m,iax)
     &              * itxnm_sum(m,iax) * itynm_sum(m,iax)
     &              * itznm_sum(m,iax)
              mnmax = mnmax * iaxnum
            ! sum of tr(np,ne,nt,nx*ny*nz)
            endif

          else if( itmsh(m) .eq. 4 ) then ! tdepsttet

            if(iax <= itaxn(m)) then
              iaxnum = 0
              select case(itaxs(m,iax))
              case (1)   ! eng
                itenm_sum(m,iax) = 1
                iaxnum = 1
              case (2)   ! reg
                itrgn_sum(m,iax) = 1
                iaxnum = 1
              case (11)   ! t
                ittnm_sum(m,iax) = 1
                iaxnum = 1
              case (12)   ! t-eng

              case (13)   ! eng-t

              case default
                
              end select

              mnmax = ( itenm_sum(m,iax) + 1 ) * itpan_sum(m,iax)
     &              * mftal_sum(m,iax) * itrgn_sum(m,iax)
     &              * ittnm_sum(m,iax)
              mnmax = mnmax * iaxnum
            ! sum of tr(np,ne,nr,nt)
            endif

          endif

*-----------------------------------------------------------------------
*        t-deposit2 tally
*-----------------------------------------------------------------------
*     1:eng1,  2:eng2,  3:e12.   4:e21,      5:t-e1,
*     6:e1-t,  7:t-e2,  8:e2-t.  9:t
        else if( ital(m) .eq. 14 ) then
          if( itmsh(m) .eq. 1 ) then    ! tdpst2reg

            if(iax <= itaxn(m)) then
              iaxnum = 0
              select case(itaxs(m,iax))
              case (1)   ! eng1
                itenm_sum(m,iax) = 1
                iaxnum = 1
              case (2)   ! eng2
                itenm2_sum(m,iax) = 1
                iaxnum = 1
              case (3)   ! e12

              case (4)   ! e21

              case (5)   ! t-e1

              case (6)   ! e1-t

              case (7)   ! t-e2

              case (8)   ! e2-t

              case (9)   ! t
                ittnm_sum(m,iax) = 1
                iaxnum = 1
              case default
                
              end select

              mnmax = ( itenm_sum(m,iax) + 1 )
     &              * ( itenm2_sum(m,iax) + 1 ) * itpan_sum(m,iax)
     &              * ittnm_sum(m,iax) * mftal_sum(m,iax)
              mnmax = mnmax * iaxnum
            ! sum of tr(np,0:e1,0:e2,nt,2)
            endif

          endif

*-----------------------------------------------------------------------
*         [t-sed]
*-----------------------------------------------------------------------
*     1:sed,  2:reg,    3:x.      4:y,      5:z,
*     6:r,    7:xy yx,  8:yz zy.  9:zx xz, 10:rz zr,
        else if ( ital(m) .eq. 15 ) then
          if( itmsh(m) .eq. 1 ) then ! tsedreg

            if(iax <= itaxn(m)) then
              iaxnum = 0
              select case(itaxs(m,iax))
              case (1)   ! sed
                itenm_sum(m,iax) = 1
                iaxnum = 1
              case (2)   ! reg
                itrgn_sum(m,iax) = 1
                iaxnum = 1
              case default
                
              end select

              mnmax = itenm_sum(m,iax) * itpan_sum(m,iax)
     &              * mftal_sum(m,iax) * itrgn_sum(m,iax)
              mnmax = mnmax * iaxnum
            ! sum of tr(np,ne,nr)
            endif

          else if( itmsh(m) .eq. 2 ) then ! tseddrz

            if(iax <= itaxn(m)) then
              iaxnum = 0
              select case(itaxs(m,iax))
              case (1)   ! sed
                itenm_sum(m,iax) = 1
                iaxnum = 1
              case (5)   ! z
                itznm_sum(m,iax) = 1
                iaxnum = 1
              case (6)   ! r
                itrnm_sum(m,iax) = 1
                iaxnum = 1
              case default
                
              end select

              mnmax = itenm_sum(m,iax) * itpan_sum(m,iax)
     &              * mftal_sum(m,iax) * itrnm_sum(m,iax)
     &              * itznm_sum(m,iax)
              mnmax = mnmax * iaxnum
            ! sum of tr(np,ne,nr,nz)
            endif

          else if( itmsh(m) .eq. 3 ) then ! tsedxyz

            if(iax <= itaxn(m)) then
              iaxnum = 0
              select case(itaxs(m,iax))
              case (1)   ! sed
                itenm_sum(m,iax) = 1
                iaxnum = 1
              case (3)   ! x
                itxnm_sum(m,iax) = 1
                iaxnum = 1
              case (4)   ! y
                itynm_sum(m,iax) = 1
                iaxnum = 1
              case (5)   ! z
                itznm_sum(m,iax) = 1
                iaxnum = 1
              case default
                
              end select

              mnmax = itenm_sum(m,iax) * itpan_sum(m,iax)
     &              * mftal_sum(m,iax) * itxnm_sum(m,iax)
     &              * itynm_sum(m,iax) * itznm_sum(m,iax)
              mnmax = mnmax * iaxnum
            ! sum of tr(np,ne,nx*ny*nz)
            endif
          end if

*-----------------------------------------------------------------------
*         [t-point]
*-----------------------------------------------------------------------
*     1:eng,  2:t

        else if ( ital(m) .eq. 17 ) then

            if(iax <= itaxn(m)) then
              iaxnum = 0
              select case(itaxs(m,iax))
              case (1)   ! eng
                itenm_sum(m,iax) = 1
                iaxnum = 1
              case (2)   ! t
                ittnm_sum(m,iax) = 1
                iaxnum = 1
              case default
                
              end select
              mnmax = itenm_sum(m,iax) * itpan_sum(m,iax)
     &              * itmsh_sum(m,iax) * itmst_sum(m,iax)
     &              * ittnm_sum(m,iax) * mftal_sum(m,iax)
              mnmax = mnmax * iaxnum
            ! sum of tr(np,ne,nt,nr,nm)
            endif
*-----------------------------------------------------------------------
*         [t-wwg]
*-----------------------------------------------------------------------
*     1:eng,  2:reg,  3:xy yx,  4:yz zy.  5:zx xz,
*     6:t,    7:wwg   8:x.      9:y,     10:z,
*    11:tet

        else if ( ital(m) .eq. 18 ) then
          if ( itmsh(m) .eq. 1 ) then

            if(iax <= itaxn(m)) then
              iaxnum = 0
              select case(itaxs(m,iax))
              case (1)   ! eng
                itenm_sum(m,iax) = 1
                iaxnum = 1
              case (2)   ! reg
                itrgn_sum(m,iax) = 1
                iaxnum = 1
              case (6)   ! t
                ittnm_sum(m,iax) = 1
                iaxnum = 1
              case default
                
              end select
              mnmax = itenm_sum(m,iax) * itpan_sum(m,iax)
     &              * itrgn_sum(m,iax) * itmst_sum(m,iax)
     &              * ittnm_sum(m,iax) * mftal_sum(m,iax)
              mnmax = mnmax * iaxnum
            ! sum of tr(np,ne,nt,nr,nm)
            endif

          else if ( itmsh(m) .eq. 3 ) then

            if(iax <= itaxn(m)) then
              iaxnum = 0
              select case(itaxs(m,iax))
              case (1)   ! eng
                itenm_sum(m,iax) = 1
                iaxnum = 1
              case (6)   ! t
                ittnm_sum(m,iax) = 1
                iaxnum = 1
              case (8)   ! x
                itxnm_sum(m,iax) = 1
                iaxnum = 1
              case (9)   ! y
                itynm_sum(m,iax) = 1
                iaxnum = 1
              case (10)   ! z
                itznm_sum(m,iax) = 1
                iaxnum = 1
              case default
                
              end select
              mnmax = itenm_sum(m,iax) * itpan_sum(m,iax)
     &              * mftal_sum(m,iax) * itxnm_sum(m,iax)
     &              * itynm_sum(m,iax) * itznm_sum(m,iax)
     &              * itmst_sum(m,iax) * ittnm_sum(m,iax)
              mnmax = mnmax * iaxnum
            ! sum of tr(np,ne,nt,nx*ny*nz,nm)
            endif


          else if ( itmsh(m) .eq. 4 ) then ! twwgtet

            if(iax <= itaxn(m)) then
              iaxnum = 0
              select case(itaxs(m,iax))
              case (1)   ! eng
                itenm_sum(m,iax) = 1
                iaxnum = 1
              case (2)   ! reg
                itrgn_sum(m,iax) = 1
                iaxnum = 1
              case (6)   ! t
                ittnm_sum(m,iax) = 1
                iaxnum = 1
              case default
                
              end select
              mnmax = itenm_sum(m,iax) * itpan_sum(m,iax)
     &              * itrgn_sum(m,iax) * itmst_sum(m,iax)
     &               * ittnm_sum(m,iax) * mftal_sum(m,iax)
              mnmax = mnmax * iaxnum
            ! sum of tr(np,ne,nt,nr,nm)
            endif
          endif

*-----------------------------------------------------------------------
*        [t-volume]
*-----------------------------------------------------------------------

        else if( ital(m) .eq. 21 ) then
            mnmax = itrgn_sum(m,iax) * 2
            ! tr(nr)

*-----------------------------------------------------------------------
*        [t-wwbg]
*-----------------------------------------------------------------------
*     1:wwbg,  2:xy yx,  3:yz zy.  4:zx xz,

        else if( ital(m) .eq. 22 ) then
          if ( itmsh(m) .eq. 1 ) then

            if(iax <= itaxn(m)) then
              iaxnum = 0
              select case(itaxs(m,iax))
              case (0)   !

              case default
                
              end select

              mnmax = itrgn_sum(m,iax) * 2 * mftal_sum(m,iax)
              mnmax = mnmax * iaxnum
            ! sum of tr(nr,2)
            endif

          else if ( itmsh(m) .eq. 3 ) then

            if(iax <= itaxn(m)) then
              iaxnum = 0
              select case(itaxs(m,iax))
              case (0)   !

              case default
                
              end select

              mnmax = itxnm_sum(m,iax) * itynm_sum(m,iax)
     &              * itznm_sum(m,iax) * 2 * mftal_sum(m,iax)
              mnmax = mnmax * iaxnum
            ! sum of tr(nx*ny*nz,2)
            endif

          else if ( itmsh(m) .eq. 4 ) then

            if(iax <= itaxn(m)) then
              iaxnum = 0
              select case(itaxs(m,iax))
              case (0)   !

              case default
                
              end select

              mnmax = itrgn_sum(m,iax) * 2 * mftal_sum(m,iax)
              mnmax = mnmax * iaxnum
            ! sum of tr(nr,2)
            endif
          endif

*-----------------------------------------------------------------------
        endif

      end subroutine

      subroutine ALLOCATE_TAL

        implicit real*8 (a-h,o-z)
        include 'param.inc'
        common /tall00/ itnm, ital(itlmax), itals(itlmax), italm(itlmax)

        integer m,mnmax,imax

        integer iax,mnmax_sum,imax_sum

      integer maxbch,maxcas
      common /cparm/  maxbch,maxcas
      common /mpi00/ npe, me
      common /paraj/  mstz(300), parz(300) ! S.H. 2022.12.22
      common /tall78/ itism(itlmax), itist(10,itlmax), itjst(10,itlmax),
     &                itkst(10,itlmax), itstt(itlmax), itsdd(itlmax)
      common /tall91/ itextstat(itlmax), mftal(itlmax) ! S.H. extstat 2024.3.18

      common /tall10/ itaxn(itlmax), itaxs(itlmax,6)

        if(itnm .eq. 0) return  ! T.Sato 2013/11/06
*-----------------------------------------------------------------------

        allocate(italsize(itnm))
        allocate(italhead(itnm))

        allocate(italsize_sum(itnm,6))
        allocate(italsize_2_sum(itnm,6))
        allocate(mtalsize_sum(itnm))
        allocate(italhead_sum(itnm,6))

*-----------------------------------------------------------------------

        mtalsize_sum(:) = 0

        do m=1, itnm

          call CALC_TALSIZE(m,mnmax)
          italsize(m) = mnmax

          do iax=1,6
            call CALC_TALSIZE_SUM(m,iax,mnmax_sum,mnmax2_sum)
            italsize_sum(m,iax)   = mnmax_sum
            italsize_2_sum(m,iax) = mnmax2_sum
            mtalsize_sum(m) = mtalsize_sum(m) + mnmax_sum
          enddo

        enddo

*-----------------------------------------------------------------------

        italhead(1) = 1
        do m = 1, itnm-1
          italhead(m+1) = italhead(m) + italsize(m)
        end do

        do m = 1, itnm
          if(m == 1 ) then
            italhead_sum(m,1) = 1
           else
            italhead_sum(m,1) = italhead_sum(m-1,6)
     &                        + italsize_sum(m-1,6)
           endif
           do iax=1,5
            italhead_sum(m,iax+1) = italhead_sum(m,iax)
     &                            + italsize_sum(m,iax)
          enddo
        enddo

*-----------------------------------------------------------------------

        imax = sum(italsize(1:itnm)) + 1
        allocate(tr0(imax))
        tr0(1:imax) = 0.0d0
!$      allocate(tr0ref(imax))
!$      tr0ref(1:imax) = 0.0d0

        imax_sum = sum(italsize_sum(1:itnm,1:6)) +1
        imax_sum = max(imax_sum,1)
        allocate(tr0_sum(imax_sum))
        tr0_sum(1:imax_sum) = 0.0d0
!$      allocate(tr0ref_sum(imax_sum))
!$      tr0ref_sum(1:imax_sum) = 0.0d0

*-----------------------------------------------------------------------

        if ( mstz(26) .eq. 3 ) then ! S.H. 2022.12.22
           itmpanatal = sum(itism(1:itnm)) ! total mesh points in anatally subsection
         if ( npe .le. 1 ) then
           allocate(deist(4,itmpanatal*maxbch))
           deist(1:4,1:itmpanatal*maxbch) = 0.0d0
         else ! when MPI
           maxbchmpi = max( 1, maxbch / ( npe - 1 ) )
           allocate(deist(4,itmpanatal*maxbchmpi))
           deist(1:4,1:itmpanatal*maxbchmpi) = 0.0d0
         end if
        end if

*-----------------------------------------------------------------------

c S.H. 2024.3.4 extstat
        allocate(iprodenhead(itnm))
        iprodenhead(1) = 1
        nprodensize = 0
        nprodentally = 0
        do m = 1, itnm-1
           if ( mftal(m) .gt. 5 ) then
              iprodenhead(m+1) = iprodenhead(m) + italsize(m)/mftal(m)
              nprodentally = nprodentally + 1
              nprodensize = nprodensize + italsize(m)/mftal(m)
           end if
        end do
        if ( mftal(itnm) .gt. 5 ) then
              nprodentally = nprodentally + 1
              nprodensize = nprodensize + italsize(itnm)/mftal(itnm)
        end if
        if ( nprodentally .gt. 0 ) then
           allocate(prodenmnmx(2,nprodensize))
           prodenmnmx(1:2,1:nprodensize) = 0.0d0
        end if

*-----------------------------------------------------------------------
c S.H. 2024.3.4 extstat

      end subroutine

      subroutine DEALLOCATE_TAL

        implicit real*8 (a-h,o-z)
        include 'param.inc'
        common /tall00/ itnm, ital(itlmax), itals(itlmax), italm(itlmax)
        common /paraj/  mstz(300), parz(300) ! S.H. 2022.12.22

        if(itnm .eq. 0) return  ! T.Sato 2013/11/06

        deallocate(tr0)
        deallocate(italsize)
        deallocate(italhead)
!$      deallocate(tr0ref)
        if ( mstz(26) .eq. 3 ) then ! S.H. 2022.12.22
           deallocate(deist)    ! S.H. 2022.12.19
        end if

        deallocate(tr0_sum)
        deallocate(italsize_sum)
        deallocate(italsize_2_sum)
        deallocate(mtalsize_sum)
        deallocate(italhead_sum)
!$      deallocate(tr0ref_sum)

      end subroutine

      end module TALMOD

!$    module TALMOD0

!$    use TALMOD
!$    implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

!$    real(8),allocatable,save,target:: tr00(:)
!$    integer,allocatable,save,target:: italsize0(:)
!$    integer,allocatable,save,target:: italhead0(:)
!$OMP THREADPRIVATE( tr00 )
!$OMP THREADPRIVATE( italsize0 )
!$OMP THREADPRIVATE( italhead0 )

* sumover
 
!$    real(8),allocatable,save,target:: tr00_sum(:)
!$    integer,allocatable,save,target:: italsize0_sum(:,:)
!$    integer,allocatable,save,target:: italsize0_2_sum(:,:)
!$    integer,allocatable,save,target:: mtalsize0_sum(:)
!$    integer,allocatable,save,target:: italhead0_sum(:,:)
!$OMP THREADPRIVATE( tr00_sum )
!$OMP THREADPRIVATE( italsize0_sum )
!$OMP THREADPRIVATE( italsize0_2_sum )
!$OMP THREADPRIVATE( mtalsize0_sum )
!$OMP THREADPRIVATE( italhead0_sum )

*-----------------------------------------------------------------------
!$    contains

!$    subroutine INIT_TALMOD0
!$      implicit real*8 (a-h,o-z)
!$      include 'param.inc'
!$      common /tall00/ itnm, ital(itlmax), itals(itlmax), italm(itlmax)
!$      integer italsh
!$      common /talsh/ italsh

!$      integer imax,imax_sum

!$      if(itnm .eq. 0) return  ! T.Sato 2013/11/06

!$      if (italsh .eq. 0 ) then
!$        allocate(italsize0(itnm))
!$        allocate(italhead0(itnm))

!$        allocate(italsize0_sum(itnm,6))
!$        allocate(italsize0_2_sum(itnm,6))
!$        allocate(mtalsize0_sum(itnm))
!$        allocate(italhead0_sum(itnm,6))

!$        italsize0(:) = italsize(:)
!$        italhead0(:) = italhead(:)

!$        imax = sum(italsize0(1:itnm)) + 1
!$        allocate(tr00(imax))
!$        tr00(:) = tr0(:)

!$        italsize0_sum(:,:) = italsize_sum(:,:)
!$        italsize0_2_sum(:,:) = italsize_2_sum(:,:)
!$        mtalsize0_sum(:) = mtalsize_sum(:)
!$        italhead0_sum(:,:) = italhead_sum(:,:)

!$        imax_sum = sum(italsize0_sum(1:itnm,1:6)) +1
!$        imax_sum = max(imax_sum,1)
!$        allocate(tr00_sum(imax_sum))

!$        tr00_sum(:) = tr0_sum(:)

!$      end if
!$    end subroutine

!$    subroutine GET_TR_HEAD_POINTER0(p,m)
!$      implicit none
!$      real(8),pointer,intent(out):: p(:)
!$      integer,intent(in):: m
!$      p => tr00(italhead0(m):)
!$    end subroutine

!$    subroutine GET_TR_HEAD_POINTER0_SUM(p,m,iax)
!$      implicit none
!$      real(8),pointer,intent(out):: p(:)
!$      integer,intent(in):: m,iax
!$      p => tr00_sum(italhead0_sum(m,iax):)
!$    end subroutine

!$    subroutine GET_TZ_HEAD_POINTER0_SUM(p,m,iax)
!$      implicit none
!$      real(8),pointer,intent(out):: p(:)
!$      integer,intent(in):: m,iax
!$      p => tr00_sum(italhead0_sum(m,iax)+italsize0_2_sum(m,iax):)
!$    end subroutine

!$    subroutine GET_TR_HEAD_POINTER0_SUM_NTF(p,m,ntf,iax)
!$      implicit none
!$      real(8),pointer,intent(out):: p(:)
!$      integer,intent(in):: m,iax,ntf
!$      integer :: ntfbase
!$      ntfbase = mtalsize0_sum(m) * (ntf - 1)
!$      p => tr00_sum(italhead0_sum(m,iax)+ntfbase:)
!$    end subroutine

!$    subroutine GET_TZ_HEAD_POINTER0_SUM_NTF(p,m,ntf,iax)
!$      implicit none
!$      real(8),pointer,intent(out):: p(:)
!$      integer,intent(in):: m,iax,ntf
!$      integer :: ntfbase
!$      ntfbase = mtalsize0_sum(m) * (ntf - 1)
!$      p =>
!$   &  tr00_sum(italhead0_sum(m,iax)+italsize0_2_sum(m,iax)+ntfbase:)
!$    end subroutine

!$    subroutine DEALLOCATE_TAL0

!$      implicit real*8 (a-h,o-z)
!$      include 'param.inc'
!$      common /tall00/ itnm, ital(itlmax), itals(itlmax), italm(itlmax)

!$      if(itnm .eq. 0) return  ! T.Sato 2013/11/06

!$      deallocate(tr00)
!$      deallocate(italsize0)
!$      deallocate(italhead0)

!$      deallocate(tr00_sum)
!$      deallocate(italsize0_sum)
!$      deallocate(italsize0_2_sum)
!$      deallocate(mtalsize0_sum)
!$      deallocate(italhead0_sum)

!$    end subroutine

!$    subroutine UPDATE_TAL00REF
!$      implicit real*8 (a-h,o-z)  ! T.Sato 2025/03/19
!$      include 'param.inc'        ! T.Sato 2025/03/19
!$      common /tall00/ itnm, ital(itlmax), itals(itlmax), italm(itlmax)  ! T.Sato 2025/03/19

!$      if(itnm .eq. 0) return     ! T.Sato 2025/03/19

!$      tr0ref(:) = tr0(:)
!$      tr0ref_sum(:) = tr0_sum(:)
!$    end subroutine

!$    subroutine SHOW_TAL
!$      implicit real*8 (a-h,o-z)
!$      include 'param.inc'
!$      common /tall00/ itnm, ital(itlmax), itals(itlmax), italm(itlmax)
!$      integer i
!$      do i=1, sum(italsize0(1:itnm))
!$      enddo
!$      do i=1, sum(italsize0_sum(1:itnm,1:6))
!$      enddo
!$    end subroutine

!$    subroutine COMPOSE_TALMOD
!$      implicit real*8 (a-h,o-z)
!$      include 'param.inc'
!$      common /tall00/ itnm, ital(itlmax), itals(itlmax), italm(itlmax)
!$      integer i
C for nonshared_tally option
!$      integer italsh
!$      common /talsh/ italsh
!$      if (italsh .eq. 0 ) then
!$        do i=1, sum(italsize0(1:itnm))
!$            tr0(i) = tr0(i) + (tr00(i) - tr0ref(i))
!$        enddo

!$        do i=1,  sum(italsize0_sum(1:itnm,1:6))
!$            tr0_sum(i) = tr0_sum(i) + (tr00_sum(i) - tr0ref_sum(i))
!$        enddo
!$      end if

!$    end subroutine

!$    subroutine SYNC_TALMOD
C for nonshared_tally option
!$      implicit real*8 (a-h,o-z) ! T.Sato 2025/03/19
!$      include 'param.inc'       ! T.Sato 2025/03/19
!$      common /tall00/ itnm, ital(itlmax), itals(itlmax), italm(itlmax)  ! T.Sato 2025/03/19
!$      integer italsh
!$      common /talsh/ italsh

!$      if(itnm .eq. 0) return

!$      if (italsh .eq. 0 ) then
!$        tr00(:) = tr0(:)
!$        tr00_sum(:) = tr0_sum(:)
!$      end if
!$    end subroutine

!$    end module TALMOD0
C=======================================================================

