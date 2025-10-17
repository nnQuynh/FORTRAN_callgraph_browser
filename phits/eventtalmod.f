************************************************************************
      module EVENTTALMOD
*                                                                      *
*       Preparation of tally in each thread                            *
*       for event by event analysis                                    *
*       [T-Deposit2]                                                   *
*       [T-Deposit] output=deposit                                     *
*       [T-Heat] output=deposit-*                                      *
************************************************************************
      use partmod, only: itpan, itpat, jtpat ! frtati 2021/10/05

*-----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      real(8),allocatable,save:: trEVENT0(:)
!$OMP THREADPRIVATE( trEVENT0 )
      integer,allocatable,save:: iomptalm(:)
!$OMP THREADPRIVATE( iomptalm )
      logical,save:: FIRSTsrc=.true.
!$OMP THREADPRIVATE( FIRSTsrc )

      integer,allocatable,save:: itrmaxm(:)
!$OMP THREADPRIVATE( itrmaxm )
      integer,allocatable,save:: itrminm(:)
!$OMP THREADPRIVATE( itrminm )
      integer,allocatable,save:: itridxm(:)
!$OMP THREADPRIVATE( itridxm )
      integer,allocatable,save:: ltridxm(:)
!$OMP THREADPRIVATE( ltridxm )

      contains
!------------------------------------------------------------------------
      subroutine ALLOCATE_EVENTTAL
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
      common /tall75/ itnlatcel(itlmax), itnlatmem(itlmax)
      common /tall83/ itnzn(itlmax), itndm(itlmax) ! frtati 2022/02/18

! T.Sato 2024/03/18 weighted history counter ID
      common /tall92/ ichnum(itlmax),chbias(itlmax),pedest(itlmax)
     &               ,ictnum(itlmax),ctbias(itlmax),ictidx(itlmax)

      integer iompmax,m
      integer iim

      allocate( iomptalm(itnm) )
      iompmax=1

*-----------------------------------------------------------------------

      do m=1, itnm

        iomptalm(m)=iompmax

*-----------------------------------------------------------------------
        if( ital(m) .eq. 1 ) then   !OBINATA(2012.5.29):Ct in t-track
*-----------------------------------------------------------------------

          if ( itmsh(m) .eq. 1 ) then   ! ttracreg

            iompmax=iompmax
     &           + itenm(m) * itpan(m) * itrgn(m) * itmst(m) * ittnm(m)
            ! trEVENT(np,ne,nt,nr,nm)

          else if( itmsh(m) .eq. 2 ) then ! ttracrz

            iompmax=iompmax
     &           + itenm(m) * itpan(m)
     &           * itrnm(m) * itznm(m) * itanm(m) * itmst(m)
     &           * ittnm(m)
            ! trEVENT(np,ne,nt,nr*nz,na,nm)


          else if( itmsh(m) .eq. 3 ) then ! ttraxyz

            iompmax=iompmax
     &           + itenm(m) * itpan(m)
     &           * itxnm(m) * itynm(m) * itznm(m) * itmst(m)
     &           * ittnm(m)
            ! trEVENT(np,ne,nt,nx*ny*nz,nm)

          else if ( itmsh(m) .eq. 4 ) then   ! ttractet
            iompmax=iompmax
     &           + itenm(m) * itpan(m) * itrgn(m) * itmst(m) * ittnm(m)
            ! trEVENT(np,ne,nt,nr,nm)

          endif

*-----------------------------------------------------------------------
        else if( ital(m) .eq. 19 ) then  ! K.Niita adjoint
*-----------------------------------------------------------------------

          if ( itmsh(m) .eq. 1 ) then   ! tadjreg

            iompmax=iompmax
     &           + itenm(m) * itpan(m) * itrgn(m) * itmst(m) * ittnm(m)
            ! trEVENT(np,ne,nt,nr,nm)

          else if( itmsh(m) .eq. 2 ) then ! tadjrz

            iompmax=iompmax
     &           + itenm(m) * itpan(m)
     &           * itrnm(m) * itznm(m) * itanm(m) * itmst(m)
     &           * ittnm(m)
            ! trEVENT(np,ne,nt,nr*nz,na,nm)

          else if( itmsh(m) .eq. 3 ) then ! tadjxyz

            iompmax=iompmax
     &           + itenm(m) * itpan(m)
     &           * itxnm(m) * itynm(m) * itznm(m) * itmst(m)
     &           * ittnm(m)
            ! trEVENT(np,ne,nt,nx*ny*nz,nm)

          endif

*-----------------------------------------------------------------------
        else if( ital(m) .eq. 2 ) then
*-----------------------------------------------------------------------

          if ( itmsh(m) .eq. 1 ) then   ! tcrsreg

            iompmax=iompmax
     &                  + itenm(m) * itpan(m)
     &                  * itrcn(m) * itanm(m)
     &                  * ittnm(m) * itmst(m)

          else if ( itmsh(m) .eq. 2 ) then   ! tcrsrz

            iompmax=iompmax
     &                  + itenm(m) * itpan(m)
     &                  * ( itrnm(m) + 1 ) * itznm(m)
     &                  * itanm(m)
     &                  * ittnm(m) * itmst(m)
     &                  + itenm(m) * itpan(m)
     &                  * ( itznm(m) + 1 ) * itrnm(m)
     &                  * itanm(m)
     &                  * ittnm(m) * itmst(m)

          else if ( itmsh(m) .eq. 3 ) then   ! tcrsxyz

            iompmax=iompmax
     &                  + itenm(m) * itpan(m)
     &                  * itxnm(m) * itynm(m) * ( itznm(m) + 1 )
     &                  * itanm(m)
     &                  * ittnm(m) * itmst(m)
          endif

*-----------------------------------------------------------------------
        else if( ital(m) .eq. 4 ) then
*-----------------------------------------------------------------------

          if( itmsh(m) .eq. 1 ) then    ! thetreg

            if ( itout(m) .ge. 4 ) then
              iompmax=iompmax
     &             +itndy(m)*itrgn(m)
              ! trEVENT(nd,nr)
            else
              iompmax=iompmax
     &               + itndy(m) * ( itenm(m) + 1 )
     &               * itrgn(m) + 5
              ! trEVENT(nd,nr,0:ne)
            end if

          else if(itmsh(m) .eq. 2) then ! thetrz

            if ( itout(m) .ge. 4 ) then
              iompmax=iompmax
     &             +itndy(m)*itrnm(m)*itznm(m)
              ! trEVENT(nd,nr,nz)
            else
              iompmax=iompmax
     &               + itndy(m) * ( itenm(m) + 1 )
     &               * itrnm(m) * itznm(m) + 5
              ! trEVENT(nd,nr,nz,0:ne)
            end if

          else if(itmsh(m) .eq. 3) then ! thetxyz

            if ( itout(m) .ge. 4 ) then
              ! trEVENT(nd,nx,ny,nz)
              iompmax=iompmax
     &             +itndy(m)*itxnm(m)*itynm(m)*itznm(m)
            else
              iompmax=iompmax
     &               + itndy(m) * ( itenm(m) + 1 )
     &               * itxnm(m) * itynm(m) * itznm(m) + 5
              ! trEVENT(nd,nx,ny,nz,0:ne)
            end if

          endif

*-----------------------------------------------------------------------
        else if( ital(m) .eq. 13 .and. itout(m) .le. 1 ) then
*-----------------------------------------------------------------------

          if( itmsh(m) .eq. 1 ) then    ! tdepstreg

            iompmax=iompmax
     &             + ( itenm(m) + 1 ) * itpan(m) * 2
     &             * itrgn(m) * ittnm(m)
            ! trEVENT(np,nr,nt)

          else if(itmsh(m) .eq. 2) then ! tdepstrz

            iompmax=iompmax
     &             + ( itenm(m) + 1 ) * itpan(m) * 2
     &             * ittnm(m)
     &             * itrnm(m) * itznm(m)
            ! trEVENT(np,nt,nr,nz)

          else if(itmsh(m) .eq. 3) then ! tdepstxyz

            iompmax=iompmax
     &             + ( itenm(m) + 1 ) * itpan(m) * 2
     &             * ittnm(m)
     &             * itxnm(m) * itynm(m) * itznm(m)
            ! trEVENT(np,nt,nx*ny*nz)
          else if( itmsh(m) .eq. 4 ) then    ! tdepsttet
            iompmax=iompmax
     &             + ( itenm(m) + 1 ) * itpan(m) * 2
     &             * itrgn(m) * ittnm(m)
            ! trEVENT(np,nr,nt)

          endif

*-----------------------------------------------------------------------
        else if( ital(m) .eq. 13 .and. itout(m) .ge. 2 ) then
*-----------------------------------------------------------------------

          if( itmsh(m) .eq. 1 ) then    ! tdepstreg

cABE 2018/10/25, add nl for identify lattice cell, nl=1 for normal case
            iompmax=iompmax
     &           +itpan(m)*itrgn(m)*ittnm(m)*itnlatmem(m)
            ! trEVENT(np,nr,nt,nl)

          else if(itmsh(m) .eq. 2) then ! tdepstrz

            iompmax=iompmax
     &           +itpan(m)*ittnm(m)*itrnm(m)*itznm(m)
            ! trEVENT(np,nt,nr,nz)

          else if(itmsh(m) .eq. 3) then ! tdepstxyz

            iompmax=iompmax
     &           +itpan(m)*ittnm(m)*itxnm(m)*itynm(m)*itznm(m)
            ! trEVENT(np,nt,nx*ny*nz)
          else if( itmsh(m) .eq. 4 ) then    ! tdepsttet
            iompmax=iompmax
     &           +itpan(m)*itrgn(m)*ittnm(m)
            ! trEVENT(np,nr,nt)

          endif

*-----------------------------------------------------------------------
        else if( ital(m) .eq. 14 ) then
*-----------------------------------------------------------------------

          if( itmsh(m) .eq. 1 ) then    ! tdpst2reg

            iompmax=iompmax +itpan(m)*2*ittnm(m)
            ! trEVENT(np,2,nt)

          endif

*-----------------------------------------------------------------------
        else if( itals(m) .eq. 6 ) then ! tyildreg
*-----------------------------------------------------------------------

          iompmax=iompmax
     &           + itndz(m) * itndn(m)
     &           * itrgn(m) * (itndm(m)+1) ! frtati 2022/02/18 3 -> itndm+1

*-----------------------------------------------------------------------
        else if( itals(m) .eq. 7 ) then ! tyildrz
*-----------------------------------------------------------------------

          iompmax=iompmax
     &           + itndz(m) * itndn(m)
     &           * itrnm(m) * itznm(m) * (itndm(m)+1) ! frtati 2022/02/18 3 -> itndm+1

*-----------------------------------------------------------------------
        else if( itals(m) .eq. 8 ) then ! tyildxyz
*-----------------------------------------------------------------------

          iompmax=iompmax
     &           + itndz(m) * itndn(m)
     &           * itxnm(m) * itynm(m) * itznm(m) * (itndm(m)+1) ! frtati 2022/02/18 3 -> itndm+1

*-----------------------------------------------------------------------
        else if( itals(m) .eq. 51 ) then ! tyildtet
*-----------------------------------------------------------------------
          iompmax=iompmax
     &           + itndz(m) * itndn(m)
     &           * itrgn(m) * (itndm(m)+1) ! frtati 2022/02/18 3 -> itndm+1

*-----------------------------------------------------------------------
        else if( itals(m) .eq. 13 ) then ! tstarreg
*-----------------------------------------------------------------------

          iompmax=iompmax
     &           + itenm(m) * ittnm(m) * itpan(m)
     &           * itrgn(m)*itnlatmem(m) ! hirata 20230422 for nlatcell in t-interact

*-----------------------------------------------------------------------
        else if( itals(m) .eq. 14 ) then ! tstarrz
*-----------------------------------------------------------------------

          iompmax=iompmax
     &           + itenm(m) * itpan(m)
     &           * ittnm(m)
     &           * itrnm(m) * itznm(m)

*-----------------------------------------------------------------------
        else if( itals(m) .eq. 15 ) then ! tstarxyz
*-----------------------------------------------------------------------

          iompmax=iompmax
     &           + itenm(m) * itpan(m)
     &           * ittnm(m)
     &           * itxnm(m) * itynm(m) * itznm(m)

*-----------------------------------------------------------------------
        else if( itals(m) .eq. 16 ) then ! ttimereg
*-----------------------------------------------------------------------

          iompmax=iompmax
     &           + itenm(m) * itpan(m)
     &           * ittnm(m)
     &                     * itrgn(m)

*-----------------------------------------------------------------------
        else if( itals(m) .eq. 17 ) then ! ttimerz
*-----------------------------------------------------------------------

          iompmax=iompmax
     &           + itenm(m) * itpan(m)
     &           * ittnm(m)
     &           * itrnm(m) * itznm(m)

*-----------------------------------------------------------------------
        else if( itals(m) .eq. 18 ) then ! ttimexyz
*-----------------------------------------------------------------------

          iompmax=iompmax
     &           + itenm(m) * itpan(m)
     &           * ittnm(m)
     &           * itxnm(m) * itynm(m) * itznm(m)

*-----------------------------------------------------------------------
        else if( itals(m) .eq. 19 ) then ! tdpareg
*-----------------------------------------------------------------------

          iompmax=iompmax
     &           + itndy(m) * itpan(m)
     &           * itrgn(m)

*-----------------------------------------------------------------------
        else if( itals(m) .eq. 20 ) then ! tdparz
*-----------------------------------------------------------------------

          iompmax=iompmax
     &           + itndy(m) * itpan(m)
     &           * itrnm(m) * itznm(m)

*-----------------------------------------------------------------------
        else if( itals(m) .eq. 21 ) then ! tdpaxyz
*-----------------------------------------------------------------------

          iompmax=iompmax
     &           + itndy(m) * itpan(m)
     &           * itxnm(m) * itynm(m) * itznm(m)

*-----------------------------------------------------------------------
        else if( itals(m) .eq. 53 ) then ! tdpatet
*-----------------------------------------------------------------------
          iompmax=iompmax
     &           + itndy(m) * itpan(m)
     &           * itrgn(m)

*-----------------------------------------------------------------------
        else if( itals(m) .eq. 22 ) then ! tpdctreg
*-----------------------------------------------------------------------

          iompmax=iompmax
     &           + itenm(m) * ittnm(m) * itpan(m)
     &           * itrgn(m) * itanm(m)

*-----------------------------------------------------------------------
        else if( itals(m) .eq. 23 ) then ! tpdctrz
*-----------------------------------------------------------------------

          iompmax=iompmax
     &           + itenm(m) * itpan(m)
     &           * ittnm(m) * itanm(m)
     &           * itrnm(m) * itznm(m)

*-----------------------------------------------------------------------
        else if( itals(m) .eq. 24 ) then ! tpdctxyz
*-----------------------------------------------------------------------

          iompmax=iompmax
     &           + itenm(m) * itpan(m)
     &           * ittnm(m) * itanm(m)
     &           * itxnm(m) * itynm(m) * itznm(m)

*-----------------------------------------------------------------------
        else if( itals(m) .eq. 52 ) then ! tpdcttet
*-----------------------------------------------------------------------
          iompmax=iompmax
     &           + itenm(m) * ittnm(m) * itpan(m)
     &           * itrgn(m) * itanm(m)

*-----------------------------------------------------------------------
        else if( itals(m) .eq. 28 ) then ! tletreg
*-----------------------------------------------------------------------

          iompmax=iompmax
     &           + itenm(m) * itpan(m)
     &           * itrgn(m)

*-----------------------------------------------------------------------
        else if( itals(m) .eq. 29 ) then ! tletdrz
*-----------------------------------------------------------------------

          iompmax=iompmax
     &           + itenm(m) * itpan(m)
     &           * itrnm(m) * itznm(m)

*-----------------------------------------------------------------------
        else if( itals(m) .eq. 30 ) then ! tletxyz
*-----------------------------------------------------------------------

          iompmax=iompmax
     &           + itenm(m) * itpan(m)
     &           * itxnm(m) * itynm(m) * itznm(m)

*-----------------------------------------------------------------------
        else if( itals(m) .eq. 35 ) then ! tsedreg
*-----------------------------------------------------------------------

          iompmax=iompmax
     &           + itenm(m) * itpan(m)
     &           * itrgn(m)

*-----------------------------------------------------------------------
        else if( itals(m) .eq. 36 ) then ! tsedrz
*-----------------------------------------------------------------------

          iompmax=iompmax
     &           + itenm(m) * itpan(m)
     &           * itrnm(m) * itznm(m)

*-----------------------------------------------------------------------
        else if( itals(m) .eq. 37 ) then ! tsedxyz
*-----------------------------------------------------------------------

          iompmax=iompmax
     &           + itenm(m) * itpan(m)
     &           * itxnm(m) * itynm(m) * itznm(m)

*-----------------------------------------------------------------------
        else if( itals(m) .eq. 41 ) then ! tpoint
*-----------------------------------------------------------------------

          iompmax=iompmax
     &           + itenm(m) * itpan(m) * itmsh(m)
     &           * itmst(m) * ittnm(m)

*-----------------------------------------------------------------------
        else if( itals(m) .eq. 42 ) then ! twwg
*-----------------------------------------------------------------------

            iompmax=iompmax
     &           + itenm(m) * itpan(m) * itrgn(m) * itmst(m) * ittnm(m)
     &           * ictidx(m)  ! T.Sato 2024/05/18
            ! trEVENT(np,ne,nt,nr,nm)

*-----------------------------------------------------------------------
        else if( itals(m) .eq. 48 ) then ! twwg
*-----------------------------------------------------------------------

            iompmax=iompmax
     &           + itenm(m) * itpan(m) * itmst(m) * ittnm(m)
     &           * itxnm(m) * itynm(m) * itznm(m)
     &           * ictidx(m)  ! T.Sato 2024/05/18
            ! trEVENT(np,ne,nt,nx*ny*nz,nm)

*-----------------------------------------------------------------------
        else if( itals(m) .eq. 55 ) then ! twwg
*-----------------------------------------------------------------------
cKN 2015/12/27

            iompmax=iompmax
     &           + itenm(m) * itpan(m) * itrgn(m) * itmst(m) * ittnm(m)
     &           * ictidx(m)  ! T.Sato 2024/05/18
            ! trEVENT(np,ne,nt,nr,nm)

*-----------------------------------------------------------------------
        else if( itals(m) .eq. 46 ) then ! tvlm
*-----------------------------------------------------------------------

            iompmax=iompmax
     &           + itrgn(m)
            ! trEVENT(nr)

*-----------------------------------------------------------------------
        else if( itals(m) .eq. 47 ) then ! twwbg
*-----------------------------------------------------------------------

            iompmax=iompmax
     &           + itrgn(m) * 2
            ! trEVENT(nr,2)

        else if( itals(m) .eq. 56 ) then ! twwbg
*-----------------------------------------------------------------------
cKN 2017/02/18

            iompmax=iompmax
     &           + itxnm(m) * itynm(m) * itznm(m) * 2
            ! trEVENT(nx*ny*nz,2)

        else if( itals(m) .eq. 57 ) then ! twwbg
*-----------------------------------------------------------------------
cKN 2017/02/18

            iompmax=iompmax
     &           + itrgn(m) * 2
            ! trEVENT(nr,2)

*-----------------------------------------------------------------------

        endif

      enddo

      allocate( trEVENT0(iompmax-1) )
      trEVENT0(1:iompmax-1)=0.0d0

      call ALLOCATE_EVENTTAL_ITR

      end subroutine ALLOCATE_EVENTTAL
!-----------------------------------------------------------------------
      subroutine ALLOCATE_EVENTTAL_ITR
      include 'param.inc'
      common /tall00/ itnm, ital(itlmax), itals(itlmax), italm(itlmax)
      common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                itout(itlmax), ittwo(itlmax)

      allocate( itridxm(itnm) )
      allocate( ltridxm(itnm) )

      itridxm(:)  = 0
      ltridxm(:)  = 0

      iim =1

      do m=1, itnm

*-----------------------------------------------------------------------
*        t-track tally
*-----------------------------------------------------------------------
        if( ital(m) .eq. 1 ) then
          if ( itmsh(m) .eq. 1 ) then   ! ttracreg
            ltridxm(m)  = 5
            ! trRES(np,ne,nt,nr,nm)

          else if ( itmsh(m) .eq. 2 ) then   ! ttracrz
            ltridxm(m)  = 7
            ! trRES(np,ne,nt,nr*nz,na,nm)


          else if ( itmsh(m) .eq. 3 ) then   ! ttracxyz
            ltridxm(m)  = 7
            ! trRES(np,ne,nt,nx*ny*nz,nm)
cFURUTA20190121 bugfix
          else if ( itmsh(m) .eq. 4 ) then   ! ttractet
            ltridxm(m)  = 5
            ! trRES(np,ne,nt,nr,nm)
          endif

*-----------------------------------------------------------------------
*        t-adjoint tally
*-----------------------------------------------------------------------

        else if( ital(m) .eq. 19 ) then

          if ( itmsh(m) .eq. 1 ) then   ! tadjreg
            ltridxm(m)  = 5
            ! trRES(np,ne,nt,nr,nm)

          else if ( itmsh(m) .eq. 2 ) then   ! tadjrz
            ltridxm(m)  = 7
            ! trRES(np,ne,nt,nr*nz,na,nm)

          else if ( itmsh(m) .eq. 3 ) then   ! tadjxyz
            ltridxm(m)  = 7
            ! trRES(np,ne,nt,nx*ny*nz,nm)
          endif

*-----------------------------------------------------------------------
*        t-cross tally
*-----------------------------------------------------------------------

        else if( ital(m) .eq. 2 ) then
          if ( itmsh(m) .eq. 1 ) then   ! tcrsreg
            ltridxm(m)  = 6

          else if ( itmsh(m) .eq. 2 ) then   ! tcrsrz
            ltridxm(m)  = 7 * 2

          else if ( itmsh(m) .eq. 3 ) then   ! tcrsxyz
            ltridxm(m)  = 8

          endif

*-----------------------------------------------------------------------
*         [t-yield]
*-----------------------------------------------------------------------
        else if ( ital(m) .eq. 3 ) then
          if ( itmsh(m) .eq. 1 ) then   ! tyildreg
            ltridxm(m)  = 3
          else if ( itmsh(m) .eq. 2 ) then   ! tyildrz
            ltridxm(m)  = 4
          else if ( itmsh(m) .eq. 3 ) then   ! tyildxyz
            ltridxm(m)  = 5
          else if ( itmsh(m) .eq. 4 ) then   ! tyildtet
            ltridxm(m)  = 3
          endif

*-----------------------------------------------------------------------
*        t-heat tally
*-----------------------------------------------------------------------
        else if( ital(m) .eq. 4 ) then
          if( itmsh(m) .eq. 1 ) then    ! thetreg
            ltridxm(m)  = 3
            ! trRES(nd,nr)
          else if(itmsh(m) .eq. 2) then ! thetrz
            ltridxm(m)  = 4
            ! trRES(nd,nr,nz)
          else if(itmsh(m) .eq. 3) then ! thetxyz
            ltridxm(m)  = 5
            ! trRES(nd,nx,ny,nz)
          endif

*-----------------------------------------------------------------------
*         [t-star]
*-----------------------------------------------------------------------
        else if ( ital(m) .eq. 5 ) then
          if( itmsh(m) .eq. 1 ) then ! tstarreg
            ltridxm(m)  = 4
          ! trRES(np,ne,nt,nr,2)
          else if( itmsh(m) .eq. 2 ) then ! tstarrz
            ltridxm(m)  = 5
          ! trRES(np,ne,nt,nr,nz,2)
          else if( itmsh(m) .eq. 3 ) then ! tstarxyz
            ltridxm(m)  = 6
          ! trRES(np,ne,nt,nx,ny,nz,2)
          endif

*-----------------------------------------------------------------------
*         [t-time]
*-----------------------------------------------------------------------
        else if ( ital(m) .eq. 6 ) then
          if( itmsh(m) .eq. 1 ) then ! ttimereg
            ltridxm(m)  = 4
          ! trRES(np,ne,nt,nr,2)
          else if( itmsh(m) .eq. 2 ) then ! ttimedrz
            ltridxm(m)  = 5
          ! trRES(np,ne,nt,nr,nz,2)
          else if( itmsh(m) .eq. 3 ) then ! ttimexyz
            ltridxm(m)  = 6
          ! trRES(np,ne,nt,nx,ny,nz,2)
          endif

*-----------------------------------------------------------------------
*         [t-dpa]
*-----------------------------------------------------------------------
        else if ( ital(m) .eq. 7 ) then
          if( itmsh(m) .eq. 1 ) then ! tdpareg
            ltridxm(m)  = 3
          ! trRES(nd,np,nr,2)
          else if( itmsh(m) .eq. 2 ) then ! tdpadrz
            ltridxm(m)  = 4
          ! trRES(nd,np,nr,nz,2)
          else if( itmsh(m) .eq. 3 ) then ! tdpaxyz
            ltridxm(m)  = 5
          ! trRES(nd,np,nx,ny,nz,2)
          else if( itmsh(m) .eq. 4 ) then ! tdpatet
            ltridxm(m)  = 3
          ! trRES(nd,np,nr,2)
          endif

*-----------------------------------------------------------------------
*         [t-product]
*-----------------------------------------------------------------------
        else if ( ital(m) .eq. 8 ) then
          if( itmsh(m) .eq. 1 ) then ! tpdctreg
            ltridxm(m)  = 5
          ! trRES(np,ne,nt,na,nr,2)
          else if( itmsh(m) .eq. 2 ) then ! tpdctdrz
            ltridxm(m)  = 6
          ! trRES(np,ne,nt,na,nr,nz,2)
          else if( itmsh(m) .eq. 3 ) then ! tpdctxyz
            ltridxm(m)  = 7
          ! trRES(np,ne,nt,na,nx*ny*nz,2)
          else if( itmsh(m) .eq. 4 ) then ! tpdcttet
            ltridxm(m)  = 5
          ! trRES(np,ne,nt,na,nr,2)
          endif

*-----------------------------------------------------------------------
*         [t-let]
*-----------------------------------------------------------------------
          else if ( ital(m) .eq. 12 ) then
            if( itmsh(m) .eq. 1 ) then ! tletreg
              ltridxm(m)  = 3
            ! trRES(np,ne,nr)
            else if( itmsh(m) .eq. 2 ) then ! tletdrz
              ltridxm(m)  = 4
            ! trRES(np,ne,nr,nz)
            else if( itmsh(m) .eq. 3 ) then ! tletxyz
              ltridxm(m)  = 5
            ! trRES(np,ne,nx*ny*nz)
            endif


*-----------------------------------------------------------------------
*        t-deposit tally
*-----------------------------------------------------------------------
        else if( ital(m) .eq. 13 ) then
          if( itmsh(m) .eq. 1 ) then    ! tdepstreg
            ltridxm(m)  = 4
            ! trRES(np,ne,nr,nt)
          else if(itmsh(m) .eq. 2) then ! tdepstrz
            ltridxm(m)  = 5
            ! trRES(np,ne,nt,nr,nz)
          else if(itmsh(m) .eq. 3) then ! tdepstxyz
            ltridxm(m)  = 6
            ! trRES(np,ne,nt,nx*ny*nz)
          else if( itmsh(m) .eq. 4 ) then ! tdepsttet
            ltridxm(m)  = 4
            ! trRES(np,ne,nr,nt)
          endif

*-----------------------------------------------------------------------
*        t-deposit2 tally
*-----------------------------------------------------------------------
        else if( ital(m) .eq. 14 ) then
          if( itmsh(m) .eq. 1 ) then    ! tdpst2reg
            ltridxm(m)  = 4
            ! trRES(np,0:e1,0:e2,nt,2)
          endif

*-----------------------------------------------------------------------
*         [t-let]
*-----------------------------------------------------------------------
        else if ( ital(m) .eq. 15 ) then
          if( itmsh(m) .eq. 1 ) then ! tletreg
            ltridxm(m)  = 3
            ! trRES(np,ne,nr)
          else if( itmsh(m) .eq. 2 ) then ! tletdrz
            ltridxm(m)  = 4
            ! trRES(np,ne,nr,nz)
          else if( itmsh(m) .eq. 3 ) then ! tletxyz
            ltridxm(m)  = 5
            ! trRES(np,ne,nx*ny*nz)
          end if

*-----------------------------------------------------------------------
*        t-point tally
*-----------------------------------------------------------------------

        else if( ital(m) .eq. 17 ) then
            ltridxm(m)  = 5
            ! trRES(np,ne,nt,nr,nm)

*-----------------------------------------------------------------------
*        t-wwg tally
*-----------------------------------------------------------------------

        else if( ital(m) .eq. 18 ) then
          if( itmsh(m) .eq. 1 ) then
            ltridxm(m)  = 5
            ! trRES(np,ne,nt,nr,nm)
          else if( itmsh(m) .eq. 3 ) then
            ltridxm(m)  = 7
            ! trRES(np,ne,nt,nx*ny*nz,nm)
          end if

*-----------------------------------------------------------------------
*        t-volume tally
*-----------------------------------------------------------------------

        else if( ital(m) .eq. 21 ) then
            ltridxm(m)  = 1
            ! trRES(nr)

*-----------------------------------------------------------------------
*        t-wwbg tally
*-----------------------------------------------------------------------

        else if( ital(m) .eq. 22 ) then
            ltridxm(m)  = 2
            ! trRES(nr,2)

*-----------------------------------------------------------------------
        endif

        itridxm(m)=iim
        iim=iim+ltridxm(m)

      enddo

      allocate( itrmaxm(iim) )
      allocate( itrminm(iim) )
      itrmaxm(:) = 0
      itrminm(:) = 0

      end subroutine ALLOCATE_EVENTTAL_ITR
!------------------------------------------------------------------------
      subroutine DEALLOCATE_EVENTTAL

      deallocate( trEVENT0 )
      call DEALLOCATE_EVENTTAL_ITR

      end subroutine DEALLOCATE_EVENTTAL
!------------------------------------------------------------------------
      subroutine DEALLOCATE_EVENTTAL_ITR

      deallocate( itridxm )
      deallocate( ltridxm )
      deallocate( itrmaxm )
      deallocate( itrminm )

      end subroutine DEALLOCATE_EVENTTAL_ITR
!------------------------------------------------------------------------
      end module EVENTTALMOD
