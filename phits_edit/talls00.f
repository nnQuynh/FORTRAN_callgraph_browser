! Tally source program location
! [t-heat]: talls01.f
! [t-track]: talls01.f
! [t-cross]: talls02.f
! [t-yield]: talls03.f
! [t-time]: talls04.f
! [t-star]: talls04.f
! [t-dpa]: talls04.f
! [t-product]: talls05.f
! [t-gshow]: talls05.f
! [t-let]: talls06.f
! [t-deposit]:talls06.f
! [t-deposit2]: talls06.f
! [t-sed]: talls07.f
! [t-heat] event-by-event: talls08.f
! [t-deposit] event-by-event: talls08.f
! [t-deposit2] event-by-event: talls08.f
! [t-dchain]: tallsm4.f

************************************************************************
*                                                                      *
      subroutine talls01(ncol)
*                                                                      *
*       tally the status of particle transpot.                         *
*       modified by K.Niita on 2005/11/14                              *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*       ncol = 0 : end of batch for deposit of heat tally              *
*                                                                      *
*              1 : start of calculation                                *
*              2 : end of calculation                                  *
*              3 : end of a batch                                      *
*              4 : source                                              *
*              5 : detection of geometry error                         *
*              6 : recovery of geometry error                          *
*              7 : termination by geometry error                       *
*              8 : termination by weight cut-off                       *
*              9 : termination by time cut-off                         *
*             10 : geometry boundary crossing                          *
*             11 : termination by energy cut-off                       *
*             12 : termination by escape or leakage                    *
*             13 : nuclear reaction (n,x)                              *
*             14 : nuclear reaction (n,n'x)                            *
*             15 : sequential transport only for tally                 *
*             16 : surface cross for WW of xyz mesh
*             101: detector resolution,temporary used in this routine  *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
************************************************************************
      use EVENTTALMOD !FURUTA
      use RESTALMOD   !OBINATA
      use TALMOD
!$       use TALMOD0
      use sumtallymod !Hashimoto (2015.1.26)
      use MMBANKMOD !S.H. for tally unit of MeV/n (2016.12.25)
      use partmod, only: itpan, itpat, jtpat ! frtati 2021/10/05
      use moddas_mesh
      use moddas_region
      use moddas_tally

*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'param01.inc'

*-----------------------------------------------------------------------

      common /mpi00/ npe, me

      common /talout/ itall
      common /talmm/  nmmax, lmmax, itlmx

*-----------------------------------------------------------------------

      common /tcntl/  icntl, inucr

      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
!$OMP THREADPRIVATE(/icomon/)
      integer           iMeVperu
      common /cMeVperu/ iMeVperu

*-----------------------------------------------------------------------

      common /tall00/ itnm, ital(itlmax), itals(itlmax), italm(itlmax)
      common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                itout(itlmax), ittwo(itlmax)

      common /tall02/ itrgn(itlmax), itrgm(itlmax), itreg(itlmax)
      common /tall25/ itrnv(itlmax), itriv(itlmax), itrrv(itlmax)

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
      common /tall10/ itaxn(itlmax), itaxs(itlmax,6)


      common /tall12/ itman(itlmax), itmat(itlmax), itmct(itlmax),
     &                itnun(itlmax), itnuc(itlmax),
     &                itndz(itlmax), itndn(itlmax),
     &                itnkz(itlmax), itnkn(itlmax)

      common /tall13/ itrcn(itlmax), itrcr(itlmax), itrca(itlmax)
      common /tall26/ itrcs(itlmax), itrcc(itlmax), itrss(itlmax)

      common /tall16/ itnfr(itlmax), itnfz(itlmax), itnfn(itlmax)
      common /tall17/ itndy(itlmax)
      common /tall18/ ithet(itlmax)

      common /tall20/ ittty(itlmax), ittnm(itlmax), ittrg(itlmax),
     &                rttmi(itlmax), rttma(itlmax), rttdl(itlmax)

      common /tall22/ itmcn(itlmax), itmtn(itlmax), itmtt(itlmax)

      common /tall23/ itaty(itlmax), itanm(itlmax), itarg(itlmax),
     &                rtami(itlmax), rtama(itlmax), rtadl(itlmax)
      common /tall24/ itrcm(itlmax), itrcg(itlmax)

      common /tall45/ itmlp(itlmax), itmln(itlmax,6), itmst(itlmax),
     &                itmli(itlmax,6), rtmme(itlmax,6), itmnt(itlmax,6),
     &                itmpn(itlmax,6), itmpt(itlmax,6,6,2)
      common /tall48/ itmdp(itlmax,0:30), itmdf(itlmax)

      common /tall56/ itety2(itlmax), itenm2(itlmax), iterg2(itlmax),
     &                rtemi2(itlmax), rtema2(itlmax), rtedl2(itlmax)

      common /tall67/ itstd(itlmax), rtstd(itlmax)
      common / stdstop / idstop, instop

      common /tall72/ itactnm(itlmax), itactrg(itlmax), itactmax(itlmax)

      common /tall73/ itenclo(itlmax), itangform(itlmax)

      common /tall83/ itnzn(itlmax), itndm(itlmax)

      common /tall91/ itextstat(itlmax), mftal(itlmax) ! S.H. extstat 2024.3.18

! T.Sato 2024/03/18 weighted history counter ID
      common /tall92/ ichnum(itlmax),chbias(itlmax),pedest(itlmax)
     &               ,ictnum(itlmax),ctbias(itlmax),ictidx(itlmax)

! T.Sato 2017/6/27
      common /egs5cmn10/denstepold,denstepnew,deinit
      real*8            denstepold,denstepnew,deinit
!$OMP THREADPRIVATE(/egs5cmn10/)

      dimension     idas(mdas*2) !20220906frtati temporary for fbounds-check
      equivalence ( das, idas )

      real(8),pointer:: tr(:)
      real(8),allocatable:: trTEMP(:)
      real(8),allocatable:: trTEMP_sum(:)
      real(8),allocatable:: trEVENT_sum(:)

*-----------------------------------------------------------------------
      common /istcut/ ist_cut, ist_bat
C for nonshared_tally option
      integer italsh
      common /talsh/ italsh

*-----------------------------------------------------------------------
      common /stat / istdev, irestart, ireschk
      common /itl4res/ itl4flg

*-----------------------------------------------------------------------
      common /multipl/ imltp,iimlt(multmax),inmlt(multmax),
C MATSUDA 2024.11.25 (lagrange: ilmlt)
     &  ilmlt(multmax),
     &  idmlt(multmax),ismlt(multmax),impan(multmax),impat(multmax,6,3),
     &                 jmpat(multmax,6,6,2), imdfl(multmax)

C MATSUDA 2024.11.25 (multplf: file name)
      common /multplf/ imltf(multmax),lmltfile(multmax),mltfile(multmax)
      character mltfile*100
*-----------------------------------------------------------------------

! T.Sato 2014/8/21, for considering detector resolution
      if( ncol .eq. 101) then
       do m = 1, itnm
C for nonshared_tally option
!$       if(italsh .eq. 0) then
!$          call GET_TR_HEAD_POINTER0(tr,m)
!$       else
        call GET_TR_HEAD_POINTER(tr,m)
C for nonshared_tally option
!$       end if
        if( itals(m).eq.31.and.itout(m).ge.2 ) then
                allocate(trEVENT_sum(
     &             itpan(m)*(itenm(m)+1)*itrgn(m)*ittnm(m) ))
                call tdepstregEVENT(ncol,m,
     &               itmtn(m),ismte(itmtt(m)),
     &               itpan(m),itrgn(m),itrgm(m),itenm(m),
     &               ittnm(m),
     &               idas_itreg(itreg(m)),
     &               das_iterg(iterg(m)),das_ittrg(ittrg(m)),tr,
     &               trEVENT0(iomptalm(m)),trEVENT_sum)
                deallocate(trEVENT_sum)
        else if( itals(m).eq.32.and.itout(m).ge.2 ) then
                allocate(trEVENT_sum(
     &             itpan(m)*(itenm(m)+1)*ittnm(m)*itrnm(m)*itznm(m) ))
                call tdepstrzEVENT(ncol,m,
     &               itmtn(m),ismte(itmtt(m)),
     &               itpan(m),itrnm(m),itznm(m),itenm(m),
     &               ittnm(m),
     &               das_itrrg(itrrg(m)),das_itzrg(itzrg(m)),
     &               das_iterg(iterg(m)),das_ittrg(ittrg(m)),tr,
     &               trEVENT0(iomptalm(m)),trEVENT_sum)
                deallocate(trEVENT_sum)
        else if( itals(m).eq.33.and.itout(m).ge.2 ) then
                allocate(trEVENT_sum(
     &    itpan(m)*(itenm(m)+1)*ittnm(m)*itxnm(m)*itynm(m)*itznm(m) ))
                call tdepstxyzEVENT(ncol,m,
     &               itmtn(m),ismte(itmtt(m)),
     &               itpan(m),itxnm(m),itynm(m),itznm(m),itenm(m),
     &               ittnm(m),
     &               das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &               das_itzrg(itzrg(m)),
     &               das_iterg(iterg(m)),das_ittrg(ittrg(m)),tr,
     &               trEVENT0(iomptalm(m)),trEVENT_sum)
                deallocate(trEVENT_sum)
        else if( itals(m).eq.50.and.itout(m).ge.2 ) then
                allocate(trEVENT_sum(
     &             itpan(m)*(itenm(m)+1)*itrgn(m)*ittnm(m) ))
                call tdepsttetEVENT(ncol,m,
     &               itmtn(m),ismte(itmtt(m)),
     &               itpan(m),itrgn(m),itrgm(m),itreg(m),
     &               itenm(m),ittnm(m),
     &               das_iterg(iterg(m)),das_ittrg(ittrg(m)),tr,
     &               trEVENT0(iomptalm(m)),trEVENT_sum)
                deallocate(trEVENT_sum)
        else if( itals(m).eq.34) then ! T.Sato 2023/08/19
                allocate(trEVENT_sum(
     &             itpan(m)*(itenm(m)+1)*(itenm2(m)+1)*ittnm(m) ))
               call tdpst2regEVENT(ncol,m,
     &              itpan(m),itrgn(m),itrgm(m),
     &              itenm(m),itenm2(m),ittnm(m),
     &              idas_itreg(itreg(m)),
     &              das_iterg(iterg(m)),das_iterg2(iterg2(m)),
     &              das_ittrg(ittrg(m)),tr,
     &              trEVENT0(iomptalm(m)),trEVENT_sum)
                deallocate(trEVENT_sum)


        endif
       enddo
       goto 999
      endif
*------ End of revision, 2014/8/21 -------------------------------------

      if( ncol .eq. 0 .or. ncol .eq. 4 .or. ncol .ge. 9 ) then

*-----------------------------------------------------------------------

         do m = 1, itnm

C for nonshared_tally option
!$          if(italsh .eq. 0) then
!$            call GET_TR_HEAD_POINTER0(tr,m)
!$          else
            call GET_TR_HEAD_POINTER(tr,m)
!$          endif
*-----------------------------------------------------------------------

            if( itals(m) .eq. 1 ) then

               if( ist_cut .eq. 0 .or. itstd(m) .eq. 1 ) then

               call ttracreg(ncol,m,
     &              itmtn(m),ismte(itmtt(m)),
     &              itpan(m),itrgn(m),itrgm(m),itenm(m),itmst(m),
     &              ittnm(m),mftal(m),
     &              idas_itreg(itreg(m)),
     &              das_iterg(iterg(m)),das_ittrg(ittrg(m)),tr,
     &              trEVENT0(iomptalm(m)))

               end if

*-----------------------------------------------------------------------

            else if( itals(m) .eq. 2 ) then

               if( ist_cut .eq. 0 .or. itstd(m) .eq. 1 ) then

               call ttracrz(ncol,m,
     &              itmtn(m),ismte(itmtt(m)),
     &              itpan(m),itrnm(m),itznm(m),itenm(m),itmst(m),
     &              ittnm(m),itanm(m),mftal(m),
     &              das_itrrg(itrrg(m)),das_itzrg(itzrg(m)),
     &              das_itarg(itarg(m)),
     &              das_iterg(iterg(m)),das_ittrg(ittrg(m)),tr,
     &              trEVENT0(iomptalm(m)),
     &              itrmaxm(itridxm(m)),itrminm(itridxm(m)))

               end if

*-----------------------------------------------------------------------

            else if( itals(m) .eq. 3 ) then

               if( ist_cut .eq. 0 .or. itstd(m) .eq. 1 ) then
               call ttracxyz(ncol,m,
     &              itmtn(m),ismte(itmtt(m)),
     &              itpan(m),itxnm(m),itynm(m),itznm(m),itenm(m),
     &              itmst(m),ittnm(m),mftal(m),
     &              das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &              das_itzrg(itzrg(m)),
     &              das_iterg(iterg(m)),das_ittrg(ittrg(m)),tr,
     &              trEVENT0(iomptalm(m)),
     &              itrmaxm(itridxm(m)),itrminm(itridxm(m)))

               end if

*-----------------------------------------------------------------------

            else if( itals(m) .eq. 49 ) then

               call ttractet(ncol,m,
     &              itmtn(m),ismte(itmtt(m)),
     &              itpan(m),itrgn(m),itrgm(m),
     &              itenm(m),itmst(m),ittnm(m),mftal(m),
     &              idas_itreg(itreg(m)),
     &              das_iterg(iterg(m)),das_ittrg(ittrg(m)),tr,
     &              trEVENT0(iomptalm(m)))
               
*-----------------------------------------------------------------------

            else if( itals(m) .eq. 43 ) then

               call tadjntreg(ncol,m,
     &              itmtn(m),ismte(itmtt(m)),
     &              itpan(m),itrgn(m),itrgm(m),itenm(m),itmst(m),
     &              ittnm(m),
     &              idas_itreg(itreg(m)),
     &              das_iterg(iterg(m)),das_ittrg(ittrg(m)),tr,
     &              trEVENT0(iomptalm(m)))

*-----------------------------------------------------------------------

            else if( itals(m) .eq. 44 ) then

               call tadjntrz(ncol,m,
     &              itmtn(m),ismte(itmtt(m)),
     &              itpan(m),itrnm(m),itznm(m),itenm(m),itmst(m),
     &              ittnm(m),itanm(m),
     &              das_itrrg(itrrg(m)),das_itzrg(itzrg(m)),
     &              das_itarg(itarg(m)),
     &              das_iterg(iterg(m)),das_ittrg(ittrg(m)),tr,
     &              trEVENT0(iomptalm(m)),
     &              itrmaxm(itridxm(m)),itrminm(itridxm(m)))

*-----------------------------------------------------------------------

            else if( itals(m) .eq. 45 ) then

               call tadjntxyz(ncol,m,
     &              itmtn(m),ismte(itmtt(m)),
     &              itpan(m),itxnm(m),itynm(m),itznm(m),itenm(m),
     &              itmst(m),ittnm(m),
     &              das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &              das_itzrg(itzrg(m)),
     &              das_iterg(iterg(m)),das_ittrg(ittrg(m)),tr,
     &              trEVENT0(iomptalm(m)),
     &              itrmaxm(itridxm(m)),itrminm(itridxm(m)))

*-----------------------------------------------------------------------

            else if( itals(m) .eq. 4 ) then

               call tsufreg(ncol,m,
     &              itpan(m),itrcn(m),itenm(m),itanm(m),ittnm(m),
     &              itmst(m),itrcs(m),idas_itrcr(itrcr(m)),
     &              das_iterg(iterg(m)),das_itarg(itarg(m)),
     &              das_ittrg(ittrg(m)),
     &              tr,
     &              trEVENT0(iomptalm(m)))

*-----------------------------------------------------------------------

            else if( itals(m) .eq. 5 ) then

               idas0 = itenm(m) * itpan(m)
     &               * ( itrnm(m) + 1 ) * itznm(m)
     &               * itanm(m)
     &               * ittnm(m)
     &               * itmst(m)
               idasz = idas0 * 2
               idas1 = italm(m) + idas0 * 2

               call tsufrz(ncol,m,
     &              itpan(m),itrnm(m),itznm(m),itenm(m),itanm(m),
     &              ittnm(m),itmst(m),
     &              das_itrrg(itrrg(m)),das_itzrg(itzrg(m)),
     &              das_iterg(iterg(m)),das_itarg(itarg(m)),
     &              das_ittrg(ittrg(m)),
     &              tr,tr(1+idasz),
     &              trEVENT0(iomptalm(m)),
     &              trEVENT0(iomptalm(m)+idas0),
     &              itrmaxm(itridxm(m)),itrminm(itridxm(m)),
     &              itrmaxm(itridxm(m)+7),itrminm(itridxm(m)+7))


*-----------------------------------------------------------------------

            else if( itals(m) .eq. 12 ) then

               call tsufxyz(ncol,m,
     &              itpan(m),itxnm(m),itynm(m),itznm(m),
     &              itenm(m),itanm(m),
     &              ittnm(m),itmst(m),
     &              das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &              das_itzrg(itzrg(m)),
     &              das_iterg(iterg(m)),das_itarg(itarg(m)),
     &              das_ittrg(ittrg(m)),
     &              tr,
     &              trEVENT0(iomptalm(m)),
     &              itrmaxm(itridxm(m)),itrminm(itridxm(m)))

*-----------------------------------------------------------------------

            else if( itals(m) .eq. 6 ) then

               call tyilreg(ncol,m,itndz(m),itndn(m),itndm(m),
     &              itmtn(m),ismte(itmtt(m)),
     &              itrgn(m),itrgm(m),itman(m),
     &              idas_itreg(itreg(m)),ismat(itmat(m)),
     &              idas(itnkz(m)),idas(itnkn(m)),tr,
     &              trEVENT0(iomptalm(m)))

*-----------------------------------------------------------------------

            else if( itals(m) .eq. 7 ) then

               call tyilrz(ncol,m,itndz(m),itndn(m),itndm(m),
     &              itmtn(m),ismte(itmtt(m)),
     &              itrnm(m),itznm(m),itman(m),
     &              das_itrrg(itrrg(m)),das_itzrg(itzrg(m)),
     &              ismat(itmat(m)),
     &              idas(itnkz(m)),idas(itnkn(m)),tr,
     &              trEVENT0(iomptalm(m)),
     &              itrmaxm(itridxm(m)),itrminm(itridxm(m)))

*-----------------------------------------------------------------------

            else if( itals(m) .eq. 8 ) then

               call tyilxyz(ncol,m,itndz(m),itndn(m),itndm(m),
     &              itmtn(m),ismte(itmtt(m)),
     &              itxnm(m),itynm(m),itznm(m),
     &              itman(m),
     &              das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &              das_itzrg(itzrg(m)),
     &              ismat(itmat(m)),
     &              idas(itnkz(m)),idas(itnkn(m)),tr,
     &              trEVENT0(iomptalm(m)),
     &              itrmaxm(itridxm(m)),itrminm(itridxm(m)))

*-----------------------------------------------------------------------
            else if( itals(m) .eq. 51 ) then

               call tyiltet(ncol,m,itndz(m),itndn(m),itndm(m),
     &              itmtn(m),ismte(itmtt(m)),
     &              itrgn(m),itrgm(m),itman(m),
     &              ismat(itmat(m)),idas_itreg(itreg(m)),
     &              idas(itnkz(m)),idas(itnkn(m)),tr,
     &              trEVENT0(iomptalm(m)))

*-----------------------------------------------------------------------

            else if( itals(m) .eq. 9 ) then

              if( itout(m) .le. 3 )then
                call thetreg(ncol,m,itndy(m),itout(m),
     &               itmtn(m),ismte(itmtt(m)),
     &               itpan(m),itrgn(m),itrgm(m),
     &               idas_itreg(itreg(m)),itenm(m),das_iterg(iterg(m)),
     &               tr,das(ithet(m)),
     &               trEVENT0(iomptalm(m)))
              else if( itout(m) .ge. 4 )then
                allocate(trEVENT_sum(
     &             itndy(m)*itrgn(m)*(itenm(m)+1) ))
                call thetregEVENT(ncol,m,itndy(m),itout(m),
     &               itmtn(m),ismte(itmtt(m)),
     &               itpan(m),itrgn(m),itrgm(m),
     &               idas_itreg(itreg(m)),itenm(m),das_iterg(iterg(m)),
     &               tr,das(ithet(m)),
     &               trEVENT0(iomptalm(m)),trEVENT_sum)
                deallocate(trEVENT_sum)

              endif

*-----------------------------------------------------------------------

            else if( itals(m) .eq. 10 ) then

              if( itout(m) .le. 3 )then
                call thetrz(ncol,m,itndy(m),itout(m),
     &               itmtn(m),ismte(itmtt(m)),
     &               itpan(m),itrnm(m),itznm(m),
     &               das_itrrg(itrrg(m)),das_itzrg(itzrg(m)),
     &               itenm(m),das_iterg(iterg(m)),
     &               tr,das(ithet(m)),
     &               trEVENT0(iomptalm(m)),
     &               itrmaxm(itridxm(m)),itrminm(itridxm(m)))
              else if( itout(m) .ge. 4 )then
                allocate(trEVENT_sum(
     &             itndy(m)*itrnm(m)*itznm(m)*(itenm(m)+1) ))
                call thetrzEVENT(ncol,m,itndy(m),itout(m),
     &               itmtn(m),ismte(itmtt(m)),
     &               itpan(m),itrnm(m),itznm(m),
     &               das_itrrg(itrrg(m)),das_itzrg(itzrg(m)),
     &               itenm(m),das_iterg(iterg(m)),
     &               tr,das(ithet(m)),
     &               trEVENT0(iomptalm(m)),trEVENT_sum)
                deallocate(trEVENT_sum)

              endif

*-----------------------------------------------------------------------

            else if( itals(m) .eq. 11 ) then

              if( itout(m) .le. 3 )then
                call thetxyz(ncol,m,itndy(m),itout(m),
     &               itmtn(m),ismte(itmtt(m)),
     &               itpan(m),itxnm(m),itynm(m),itznm(m),
     &               das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &               das_itzrg(itzrg(m)),
     &               itenm(m),das_iterg(iterg(m)),
     &               tr,das(ithet(m)),
     &               trEVENT0(iomptalm(m)),
     &              itrmaxm(itridxm(m)),itrminm(itridxm(m)))
              else if( itout(m) .ge. 4 )then
                allocate(trEVENT_sum(
     &             itndy(m)*itxnm(m)*itynm(m)*itznm(m)*(itenm(m)+1) ))
                call thetxyzEVENT(ncol,m,itndy(m),itout(m),
     &               itmtn(m),ismte(itmtt(m)),
     &               itpan(m),itxnm(m),itynm(m),itznm(m),
     &               das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &               das_itzrg(itzrg(m)),
     &               itenm(m),das_iterg(iterg(m)),
     &               tr,das(ithet(m)),
     &               trEVENT0(iomptalm(m)),trEVENT_sum)
                deallocate(trEVENT_sum)

              endif

*-----------------------------------------------------------------------

            else if( itals(m) .eq. 13 ) then

               call tstarreg(ncol,m,
     &              itpan(m),itrgn(m),itrgm(m),itenm(m),ittnm(m),
     &              itactnm(m),                                         ! S.Abe 2018/02/15
     &              itman(m),ismat(itmat(m)),
     &              itmtn(m),ismte(itmtt(m)),
     &              idas_itreg(itreg(m)),
     &              das_iterg(iterg(m)),das_ittrg(ittrg(m)),
     &              das(itactrg(m)),tr,     ! S.Abe 2018/02/15, add das(itactrg(m))
     &              trEVENT0(iomptalm(m)))

*-----------------------------------------------------------------------

            else if( itals(m) .eq. 14 ) then

               call tstarrz(ncol,m,
     &              itpan(m),itrnm(m),itznm(m),itenm(m),ittnm(m),
     &              itman(m),ismat(itmat(m)),
     &              itmtn(m),ismte(itmtt(m)),
     &              das_itrrg(itrrg(m)),das_itzrg(itzrg(m)),
     &              das_iterg(iterg(m)),das_ittrg(ittrg(m)),tr,
     &              trEVENT0(iomptalm(m)),
     &              itrmaxm(itridxm(m)),itrminm(itridxm(m)))

*-----------------------------------------------------------------------

            else if( itals(m) .eq. 15 ) then

               call tstarxyz(ncol,m,
     &              itpan(m),itxnm(m),itynm(m),itznm(m),itenm(m),
     &              ittnm(m),
     &              itman(m),ismat(itmat(m)),
     &              itmtn(m),ismte(itmtt(m)),
     &              das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &              das_itzrg(itzrg(m)),
     &              das_iterg(iterg(m)),das_ittrg(ittrg(m)),tr,
     &              trEVENT0(iomptalm(m)),
     &              itrmaxm(itridxm(m)),itrminm(itridxm(m)))

*-----------------------------------------------------------------------

            else if( itals(m) .eq. 16 ) then

               call ttimereg(ncol,m,
     &              itmtn(m),ismte(itmtt(m)),
     &              itpan(m),itrgn(m),itrgm(m),itenm(m),ittnm(m),
     &              idas_itreg(itreg(m)),
     &              das_iterg(iterg(m)),das_ittrg(ittrg(m)),tr,
     &              trEVENT0(iomptalm(m)))

*-----------------------------------------------------------------------

            else if( itals(m) .eq. 17 ) then

               call ttimerz(ncol,m,
     &              itmtn(m),ismte(itmtt(m)),
     &              itpan(m),itrnm(m),itznm(m),itenm(m),ittnm(m),
     &              das_itrrg(itrrg(m)),das_itzrg(itzrg(m)),
     &              das_iterg(iterg(m)),das_ittrg(ittrg(m)),tr,
     &              trEVENT0(iomptalm(m)),
     &              itrmaxm(itridxm(m)),itrminm(itridxm(m)))

*-----------------------------------------------------------------------

            else if( itals(m) .eq. 18 ) then

               call ttimexyz(ncol,m,
     &              itmtn(m),ismte(itmtt(m)),
     &              itpan(m),itxnm(m),itynm(m),itznm(m),itenm(m),
     &              ittnm(m),
     &              das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &              das_itzrg(itzrg(m)),
     &              das_iterg(iterg(m)),das_ittrg(ittrg(m)),tr,
     &              trEVENT0(iomptalm(m)),
     &              itrmaxm(itridxm(m)),itrminm(itridxm(m)))

*-----------------------------------------------------------------------

            else if( itals(m) .eq. 19 ) then

               call tdpareg(ncol,m,itndy(m),itout(m),
     &              itpan(m),itrgn(m),itrgm(m),
     &              itman(m),ismat(itmat(m)),
     &              itmtn(m),ismte(itmtt(m)),
     &              idas_itreg(itreg(m)),tr,
     &              trEVENT0(iomptalm(m)))

*-----------------------------------------------------------------------

            else if( itals(m) .eq. 20 ) then

               call tdparz(ncol,m,itndy(m),itout(m),
     &              itpan(m),itrnm(m),itznm(m),
     &              itman(m),ismat(itmat(m)),
     &              itmtn(m),ismte(itmtt(m)),
     &              das_itrrg(itrrg(m)),das_itzrg(itzrg(m)),
     &              tr,
     &              trEVENT0(iomptalm(m)),
     &              itrmaxm(itridxm(m)),itrminm(itridxm(m)))

*-----------------------------------------------------------------------

            else if( itals(m) .eq. 21 ) then

               call tdpaxyz(ncol,m,itndy(m),itout(m),
     &              itpan(m),itxnm(m),itynm(m),itznm(m),
     &              itman(m),ismat(itmat(m)),
     &              itmtn(m),ismte(itmtt(m)),
     &              das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &              das_itzrg(itzrg(m)),
     &              tr,
     &              trEVENT0(iomptalm(m)),
     &              itrmaxm(itridxm(m)),itrminm(itridxm(m)))

*-----------------------------------------------------------------------
            else if( itals(m) .eq. 53 ) then

               call tdpatet(ncol,m,itndy(m),itout(m),
     &              itpan(m),itrgn(m),itrgm(m),
     &              itman(m),ismat(itmat(m)),
     &              itmtn(m),ismte(itmtt(m)),
     &              tr,
     &              idas_itreg(itreg(m)),
     &              trEVENT0(iomptalm(m)))

*-----------------------------------------------------------------------

            else if( itals(m) .eq. 22 ) then

               call tpdctreg(ncol,m,
     &              itpan(m),itrgn(m),itrgm(m),itenm(m),ittnm(m),
     &              itanm(m),itman(m),ismat(itmat(m)),
     &              itmtn(m),ismte(itmtt(m)),
     &              idas_itreg(itreg(m)),
     &              das_iterg(iterg(m)),das_ittrg(ittrg(m)),
     &              das_itarg(itarg(m)),tr,
     &              trEVENT0(iomptalm(m)))

*-----------------------------------------------------------------------

            else if( itals(m) .eq. 23 ) then

               call tpdctrz(ncol,m,
     &              itpan(m),itrnm(m),itznm(m),itenm(m),ittnm(m),
     &              itanm(m),itman(m),ismat(itmat(m)),
     &              itmtn(m),ismte(itmtt(m)),
     &              das_itrrg(itrrg(m)),das_itzrg(itzrg(m)),
     &              das_iterg(iterg(m)),das_ittrg(ittrg(m)),
     &              das_itarg(itarg(m)),tr,
     &              trEVENT0(iomptalm(m)),
     &              itrmaxm(itridxm(m)),itrminm(itridxm(m)))

*-----------------------------------------------------------------------

            else if( itals(m) .eq. 24 ) then

               call tpdctxyz(ncol,m,
     &              itpan(m),itxnm(m),itynm(m),itznm(m),itenm(m),
     &              ittnm(m),itanm(m),
     &              itman(m),ismat(itmat(m)),
     &              itmtn(m),ismte(itmtt(m)),
     &              das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &              das_itzrg(itzrg(m)),
     &              das_iterg(iterg(m)),das_ittrg(ittrg(m)),
     &              das_itarg(itarg(m)),tr,
     &              trEVENT0(iomptalm(m)),
     &              itrmaxm(itridxm(m)),itrminm(itridxm(m)))

*-----------------------------------------------------------------------

            else if( itals(m) .eq. 52 ) then
               call tpdcttet(ncol,m,
     &              itpan(m),itrgn(m),itrgm(m),
     &              itenm(m),ittnm(m),
     &              itanm(m),itman(m),ismat(itmat(m)),
     &              itmtn(m),ismte(itmtt(m)),
     &              idas_itreg(itreg(m)),
     &              das_iterg(iterg(m)),das_ittrg(ittrg(m)),
     &              das_itarg(itarg(m)),tr,
     &              trEVENT0(iomptalm(m)))

*-----------------------------------------------------------------------

            else if( itals(m) .eq. 28 ) then

               call tletreg(ncol,m,
     &              itmtn(m),ismte(itmtt(m)),
     &              itpan(m),itrgn(m),itrgm(m),itenm(m),
     &              idas_itreg(itreg(m)),
     &              das_iterg(iterg(m)),tr,
     &              trEVENT0(iomptalm(m)))

*-----------------------------------------------------------------------

            else if( itals(m) .eq. 29 ) then

               call tletrz(ncol,m,
     &              itmtn(m),ismte(itmtt(m)),
     &              itpan(m),itrnm(m),itznm(m),itenm(m),
     &              das_itrrg(itrrg(m)),das_itzrg(itzrg(m)),
     &              das_iterg(iterg(m)),tr,
     &              trEVENT0(iomptalm(m)),
     &              itrmaxm(itridxm(m)),itrminm(itridxm(m)))

*-----------------------------------------------------------------------

            else if( itals(m) .eq. 30 ) then

               call tletxyz(ncol,m,
     &              itmtn(m),ismte(itmtt(m)),
     &              itpan(m),itxnm(m),itynm(m),itznm(m),itenm(m),
     &              das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &              das_itzrg(itzrg(m)),
     &              das_iterg(iterg(m)),tr,
     &              trEVENT0(iomptalm(m)),
     &              itrmaxm(itridxm(m)),itrminm(itridxm(m)))

*-----------------------------------------------------------------------

            else if( itals(m) .eq. 31 ) then

              if( itout(m) .le. 1 )then
                call tdepstreg(ncol,m,
     &               itmtn(m),ismte(itmtt(m)),
     &               itpan(m),itrgn(m),itrgm(m),itenm(m),
     &               ittnm(m),itman(m),ismat(itmat(m)),
     &               idas_itreg(itreg(m)),
     &               das_iterg(iterg(m)),das_ittrg(ittrg(m)),tr,
     &               trEVENT0(iomptalm(m)))
              else if( itout(m) .ge. 2 )then
                allocate(trEVENT_sum(
     &             itpan(m)*(itenm(m)+1)*itrgn(m)*ittnm(m) ))
                call tdepstregEVENT(ncol,m,
     &               itmtn(m),ismte(itmtt(m)),
     &               itpan(m),itrgn(m),itrgm(m),itenm(m),
     &               ittnm(m),
     &               idas_itreg(itreg(m)),
     &               das_iterg(iterg(m)),das_ittrg(ittrg(m)),tr,
     &               trEVENT0(iomptalm(m)),trEVENT_sum)
                deallocate(trEVENT_sum)

              endif

*-----------------------------------------------------------------------

            else if( itals(m) .eq. 32 ) then

              if( itout(m) .le. 1 )then
                call tdepstrz(ncol,m,
     &               itmtn(m),ismte(itmtt(m)),
     &               itpan(m),itrnm(m),itznm(m),itenm(m),
     &               ittnm(m),itman(m),ismat(itmat(m)),
     &               das_itrrg(itrrg(m)),das_itzrg(itzrg(m)),
     &               das_iterg(iterg(m)),das_ittrg(ittrg(m)),tr,
     &               trEVENT0(iomptalm(m)),
     &               itrmaxm(itridxm(m)),itrminm(itridxm(m)))
              else if( itout(m) .ge. 2)then
                allocate(trEVENT_sum(
     &             itpan(m)*(itenm(m)+1)*ittnm(m)*itrnm(m)*itznm(m) ))
                call tdepstrzEVENT(ncol,m,
     &               itmtn(m),ismte(itmtt(m)),
     &               itpan(m),itrnm(m),itznm(m),itenm(m),
     &               ittnm(m),
     &               das_itrrg(itrrg(m)),das_itzrg(itzrg(m)),
     &               das_iterg(iterg(m)),das_ittrg(ittrg(m)),tr,
     &               trEVENT0(iomptalm(m)),trEVENT_sum)
                deallocate(trEVENT_sum)

              endif

*-----------------------------------------------------------------------

            else if( itals(m) .eq. 33 ) then

              if( itout(m) .le. 1 )then
                call tdepstxyz(ncol,m,
     &               itmtn(m),ismte(itmtt(m)),
     &               itpan(m),itxnm(m),itynm(m),itznm(m),itenm(m),
     &               ittnm(m),itman(m),ismat(itmat(m)),
     &               das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &               das_itzrg(itzrg(m)),
     &               das_iterg(iterg(m)),das_ittrg(ittrg(m)),tr,
     &               trEVENT0(iomptalm(m)),
     &               itrmaxm(itridxm(m)),itrminm(itridxm(m)))
              else if( itout(m) .ge. 2)then
                allocate(trEVENT_sum(
     &    itpan(m)*(itenm(m)+1)*ittnm(m)*itxnm(m)*itynm(m)*itznm(m) ))
                call tdepstxyzEVENT(ncol,m,
     &               itmtn(m),ismte(itmtt(m)),
     &               itpan(m),itxnm(m),itynm(m),itznm(m),itenm(m),
     &               ittnm(m),
     &               das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &               das_itzrg(itzrg(m)),
     &               das_iterg(iterg(m)),das_ittrg(ittrg(m)),tr,
     &               trEVENT0(iomptalm(m)),trEVENT_sum)
                deallocate(trEVENT_sum)

              endif

*-----------------------------------------------------------------------
            else if( itals(m) .eq. 50 ) then

              if( itout(m) .le. 1 )then
                call tdepsttet(ncol,m,
     &               itmtn(m),ismte(itmtt(m)),
     &               itpan(m),itrgn(m),itrgm(m),
     &               itenm(m),ittnm(m),itman(m),ismat(itmat(m)),
     &               idas_itreg(itreg(m)),
     &               das_iterg(iterg(m)),das_ittrg(ittrg(m)),tr,
     &               trEVENT0(iomptalm(m)))
              else if( itout(m) .ge. 2 )then
                allocate(trEVENT_sum(
     &             itpan(m)*(itenm(m)+1)*itrgn(m)*ittnm(m) ))
                call tdepsttetEVENT(ncol,m,
     &               itmtn(m),ismte(itmtt(m)),
     &               itpan(m),itrgn(m),itrgm(m),
     &               itenm(m),ittnm(m),
     &               idas_itreg(itreg(m)),
     &               das_iterg(iterg(m)),das_ittrg(ittrg(m)),tr,
     &               trEVENT0(iomptalm(m)),trEVENT_sum)
                deallocate(trEVENT_sum)


              endif

*-----------------------------------------------------------------------

            else if( itals(m) .eq. 34 ) then
                allocate(trEVENT_sum(
     &             itpan(m)*(itenm(m)+1)*(itenm2(m)+1)*ittnm(m) ))
               call tdpst2regEVENT(ncol,m,
     &              itpan(m),itrgn(m),itrgm(m),
     &              itenm(m),itenm2(m),ittnm(m),
     &              idas_itreg(itreg(m)),
     &              das_iterg(iterg(m)),das_iterg2(iterg2(m)),
     &              das_ittrg(ittrg(m)),tr,
     &              trEVENT0(iomptalm(m)),trEVENT_sum)
                deallocate(trEVENT_sum)

*-----------------------------------------------------------------------

            else if( itals(m) .eq. 35 ) then

               call tsedreg(ncol,m,
     &              itmtn(m),ismte(itmtt(m)),
     &              itpan(m),itrgn(m),itrgm(m),itenm(m),
     &              idas_itreg(itreg(m)),
     &              das_iterg(iterg(m)),tr,
     &              trEVENT0(iomptalm(m)))

*-----------------------------------------------------------------------

            else if( itals(m) .eq. 36 ) then

               call tsedrz(ncol,m,
     &              itmtn(m),ismte(itmtt(m)),
     &              itpan(m),itrnm(m),itznm(m),itenm(m),
     &              das_itrrg(itrrg(m)),das_itzrg(itzrg(m)),
     &              das_iterg(iterg(m)),tr,
     &              trEVENT0(iomptalm(m)),
     &              itrmaxm(itridxm(m)),itrminm(itridxm(m)))

*-----------------------------------------------------------------------

            else if( itals(m) .eq. 37 ) then

               call tsedxyz(ncol,m,
     &              itmtn(m),ismte(itmtt(m)),
     &              itpan(m),itxnm(m),itynm(m),itznm(m),itenm(m),
     &              das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &              das_itzrg(itzrg(m)),
     &              das_iterg(iterg(m)),tr,
     &              trEVENT0(iomptalm(m)),
     &              itrmaxm(itridxm(m)),itrminm(itridxm(m)))

*-----------------------------------------------------------------------

            else if( itals(m) .eq. 41 ) then

               if( ist_cut .eq. 0 .or. itstd(m) .eq. 1 ) then

               call tpointp(ncol,m,
     &              itpan(m),itmsh(m),itenm(m),itmst(m),ittnm(m),
     &              mftal(m),
     &              das_iterg(iterg(m)),das_ittrg(ittrg(m)),tr,
     &              trEVENT0(iomptalm(m)))

               end if
*-----------------------------------------------------------------------

            else if( itals(m) .eq. 42 ) then

               call twwgreg(ncol,m,
     &              itmtn(m),ismte(itmtt(m)),
     &              itpan(m),itrgn(m),itrgm(m),itenm(m),itmst(m),
     &              ittnm(m),mftal(m),ictidx(m), ! S.H. 2024.3.21, T.Sato 2024/05/18
     &              idas_itreg(itreg(m)),
     &              das_iterg(iterg(m)),das_ittrg(ittrg(m)),tr,
     &              trEVENT0(iomptalm(m)))

*-----------------------------------------------------------------------

            else if( itals(m) .eq. 48 ) then

               call twwgxyz(ncol,m,
     &              itmtn(m),ismte(itmtt(m)),
     &              itpan(m),itxnm(m),itynm(m),itznm(m),itenm(m),
     &              itmst(m),ittnm(m),mftal(m),ictidx(m), ! S.H. 2024.3.21, T.Sato 2024/05/18
     &              das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &              das_itzrg(itzrg(m)),
     &              das_iterg(iterg(m)),das_ittrg(ittrg(m)),tr,
     &              trEVENT0(iomptalm(m)),
     &              itrmaxm(itridxm(m)),itrminm(itridxm(m)))

*-----------------------------------------------------------------------
cFURUTA20231227

            else if( itals(m) .eq. 55 ) then

               call twwgtet(ncol,m,
     &              itmtn(m),ismte(itmtt(m)),
     &              itpan(m),itrgn(m),itrgm(m),
     &              itenm(m),itmst(m),ittnm(m),mftal(m),ictidx(m), ! S.H. 2024.3.21, T.Sato 2024/05/18
     &              idas_itreg(itreg(m)),
     &              das_iterg(iterg(m)),das_ittrg(ittrg(m)),tr,
     &              trEVENT0(iomptalm(m)))
!--->

*-----------------------------------------------------------------------
cKN 2017/02/16

            else if( itals(m) .eq. 46 .and.
     &             ( icntl .eq. 13 .or. icntl .eq. 14 ) ) then

               call tvlmreg(ncol,m,
     &              itrgn(m),itrgm(m),idas_itreg(itreg(m)),tr,
     &              trEVENT0(iomptalm(m)))

*-----------------------------------------------------------------------

            else if( itals(m) .eq. 47 .and.
     &             ( icntl .eq. 13 .or. icntl .eq. 15 ) ) then

               call twwbgreg(ncol,m,
     &              itrgn(m),itrgm(m),idas_itreg(itreg(m)),tr,
     &              trEVENT0(iomptalm(m)))

*-----------------------------------------------------------------------
cKN 2018/01/08

            else if( itals(m) .eq. 56 .and.
     &             ( icntl .eq. 13 .or. icntl .eq. 15 ) ) then

               call twwbgxyz(ncol,m,
     &              itxnm(m),itynm(m),itznm(m),
     &              das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &              das_itzrg(itzrg(m)),tr,
     &              trEVENT0(iomptalm(m)))

*-----------------------------------------------------------------------

            else if( itals(m) .eq. 57 .and.
     &             ( icntl .eq. 13 .or. icntl .eq. 15 ) ) then

               call twwbgtet(ncol,m,
     &              itrgn(m),itrgm(m),idas_itreg(itreg(m)),tr,
     &              trEVENT0(iomptalm(m)))
!--->

*-----------------------------------------------------------------------

            end if

         end do

*-----------------------------------------------------------------------

         if(ncol.eq.0)then
           FIRSTsrc=.true.
         elseif(ncol.eq.4)then
           FIRSTsrc=.false.
         endif

      end if

*-----------------------------------------------------------------------
*     ncol = 2, 3 : output of results
*-----------------------------------------------------------------------

      if( ncol .eq. 2 .or. ncol .eq. 3 ) then
C for nonshared_tally option
!$          if(italsh .eq. 0) then
!$            allocate(trTEMP(sum(italsize0(1:itnm)+1)))
!$            allocate(trTEMP_sum(sum(mtalsize0_sum(1:itnm)+1)))
!$          else
            allocate(trTEMP(sum(italsize(1:itnm)+1)))
            allocate(trTEMP_sum(sum(mtalsize_sum(1:itnm)+1)))
!$          end if

            igsh = 0

*-----------------------------------------------------------------------

         do m = 1, itnm
*-----------------------------------------------------------------------

          if( icntl.eq.13 .and. nsumtalRead(m).lt.1 ) cycle

*-----------------------------------------------------------------------

C for nonshared_tally option
!$          if(italsh .eq. 0) then
!$            call GET_TR_HEAD_POINTER0(tr,m)
!$            mnmax = italsize0(m)
!$          else
            call GET_TR_HEAD_POINTER(tr,m)
            mnmax = italsize(m)
!$          end if

*-----------------------------------------------------------------------

            if( itals(m) .eq. 1 ) then

               if( ist_cut .eq. 0 .or. itstd(m) .eq. 1 ) then

                     idas0 = nmmax
                     idas1 = lmmax
                     idas2 = idas1 + itenm(m)
                     idas3 = idas2 + ittnm(m)
                     idas4 = ( idas3 + itrgn(m) - 1 ) * 2 + 1
                     idas5 = idas3 + itrgn(m)
     &                     + ( itrgn(m) + mod(itrgn(m),2) ) / 2
                     idas6 = ( idas5 + itrgn(m) - 1 ) * 2 + 1
                     idasa = idas5

                  if( npe .gt. 1 ) then
                     call paratal(tr,trTEMP,mnmax)
                     call paratal_sumover(trTEMP_sum,m)
                  endif

               if( me .eq. 0 ) then

                  !OBINATA(2012.6.6)
                  if( npe .gt. 1 ) then
                    call sumtr(tr,trRES(irestalm(m)),mnmax,ncol)
                    call sumtr_sumover(m,ncol)   !  tr_sum -> trRES_sum
                   endif

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(0,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(0,trTEMP_sum,ncol,m)
                  endif


                  call ptracreg(m,
     &              itpan(m),itrgn(m),itrgm(m),itenm(m),itmst(m),
     &              ittnm(m),mftal(m),
     &              idas_itreg(itreg(m)),
     &              das_iterg(iterg(m)),das_ittrg(ittrg(m)),tr,
     &              itrnv(m),idas(itriv(m)),das(itrrv(m)),
     &              itxnm(m),itynm(m),itznm(m),
     &              das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &              das_itzrg(itzrg(m)),
     &              igsh,idasa)

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(1,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(1,trTEMP_sum,ncol,m)
                  endif

               end if

               end if

*-----------------------------------------------------------------------

            else if( itals(m) .eq. 2 ) then

               if( ist_cut .eq. 0 .or. itstd(m) .eq. 1 ) then
                     idas0 = nmmax
                     idas1 = lmmax
                     idas2 = idas1 + itenm(m)
                     idas3 = idas2 + ittnm(m)
                     idasa = idas3 + itanm(m)

                  if( npe .gt. 1 ) then
                     call paratal(tr,trTEMP,mnmax)
                     call paratal_sumover(trTEMP_sum,m)
                  endif

               if( me .eq. 0 ) then

                  !OBINATA(2012.6.26)
                  if( npe .gt. 1 ) then
                    call sumtr(tr,trRES(irestalm(m)),mnmax,ncol)
                    call sumtr_sumover(m,ncol)   !  tr_sum -> trRES_sum
                   endif
                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(0,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(0,trTEMP_sum,ncol,m)
                  endif

                  call ptracrz(m,
     &              itpan(m),itrnm(m),itznm(m),itanm(m),
     &              itenm(m),itmst(m),ittnm(m),mftal(m),
     &              das_itrrg(itrrg(m)),das_itzrg(itzrg(m)),
     &              das_itarg(itarg(m)),
     &              das_iterg(iterg(m)),das_ittrg(ittrg(m)),tr,
     &              idasa)


                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(1,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(1,trTEMP_sum,ncol,m)
                  endif

               end if

               end if
*-----------------------------------------------------------------------

            else if( itals(m) .eq. 3 ) then

               if( ist_cut .eq. 0 .or. itstd(m) .eq. 1 ) then

                     idas0 = nmmax
                     idas1 = lmmax
                     idas2 = idas1 + itenm(m)
                     idas3 = ( idas2 + ittnm(m) - 1 ) * 2 + 1
                     idasa = idas2 + ittnm(m)

                  if( npe .gt. 1 ) then
                     call paratal(tr,trTEMP,mnmax)
                     call paratal_sumover(trTEMP_sum,m)
                  endif

               if( me .eq. 0 ) then

                  !OBINATA(2012.6.26)
                  if( npe .gt. 1 ) then
                    call sumtr(tr,trRES(irestalm(m)),mnmax,ncol)
                    call sumtr_sumover(m,ncol)   !  tr_sum -> trRES_sum
                   endif

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(0,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(0,trTEMP_sum,ncol,m)
                  endif

                  call ptracxyz(m,
     &              itpan(m),itmtn(m),ismte(itmtt(m)),
     &              itxnm(m),itynm(m),itznm(m),itenm(m),itmst(m),
     &              ittnm(m),mftal(m),
     &              das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &              das_itzrg(itzrg(m)),
     &              das_iterg(iterg(m)),das_ittrg(ittrg(m)),tr,
     &              igsh,idasa)

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(1,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(1,trTEMP_sum,ncol,m)
                  endif

               end if

               end if
*-----------------------------------------------------------------------
            else if( itals(m) .eq. 49 ) then

                     idas0 = nmmax
                     idas1 = lmmax
                     idas2 = idas1 + itenm(m)
                     idas3 = idas2 + ittnm(m)
                     idas4 = ( idas3 + itrgn(m) - 1 ) * 2 + 1
                     idas5 = idas3 + itrgn(m)
     &                     + ( itrgn(m) + mod(itrgn(m),2) ) / 2
                     idas6 = ( idas5 + itrgn(m) - 1 ) * 2 + 1
                     idasa = idas5

                  if( npe .gt. 1 ) then
                     call paratal(tr,trTEMP,mnmax)
                     call paratal_sumover(trTEMP_sum,m)
                  endif

               if( me .eq. 0 ) then

                  !OBINATA(2012.6.6)
                  if( npe .gt. 1 ) then
                    call sumtr(tr,trRES(irestalm(m)),mnmax,ncol)
                    call sumtr_sumover(m,ncol)   !  tr_sum -> trRES_sum
                   endif

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(0,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(0,trTEMP_sum,ncol,m)
                  endif

!<-20211222murofushi update
!--                  call ptractet(m,
!--     &              itpan(m),itrgn(m),itenm(m),itmst(m),ittnm(m),
!--     &              das(iterg(m)),das(ittrg(m)),tr,
!--     &              itreg(m),
!--     &              itxnm(m),itynm(m),itznm(m),
!--     &              das(itxrg(m)),das(ityrg(m)),das(itzrg(m)),
!--     &              igsh,idasa)
                  call ptractet(m,
     &              itpan(m),itrgn(m),itrgm(m),
     &              itenm(m),itmst(m),ittnm(m),
     &              mftal(m),
     &              das_iterg(iterg(m)),das_ittrg(ittrg(m)),tr,
     &              idas_itreg(itreg(m)),
     &              itxnm(m),itynm(m),itznm(m),
     &              das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &              das_itzrg(itzrg(m)),
     &              igsh,idasa)
                  
                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(1,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(1,trTEMP_sum,ncol,m)
                  endif

               end if

*-----------------------------------------------------------------------

            else if( itals(m) .eq. 43 ) then

                     idas0 = nmmax
                     idas1 = lmmax
                     idas2 = idas1 + itenm(m)
                     idas3 = idas2 + ittnm(m)
                     idas4 = ( idas3 + itrgn(m) - 1 ) * 2 + 1
                     idas5 = idas3 + itrgn(m)
     &                     + ( itrgn(m) + mod(itrgn(m),2) ) / 2
                     idas6 = ( idas5 + itrgn(m) - 1 ) * 2 + 1
                     idasa = idas5

                  if( npe .gt. 1 ) then
                     call paratal(tr,trTEMP,mnmax)
                     call paratal_sumover(trTEMP_sum,m)
                  endif

               if( me .eq. 0 ) then

                  !OBINATA(2012.6.6)
                  if( npe .gt. 1 ) then
                    call sumtr(tr,trRES(irestalm(m)),mnmax,ncol)
                    call sumtr_sumover(m,ncol)   !  tr_sum -> trRES_sum
                   endif

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(0,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(0,trTEMP_sum,ncol,m)
                  endif

                  call padjntreg(m,
     &              itpan(m),itrgn(m),itrgm(m),itenm(m),itmst(m),
     &              ittnm(m),
     &              idas_itreg(itreg(m)),
     &              das_iterg(iterg(m)),das_ittrg(ittrg(m)),tr,
     &              itrnv(m),idas(itriv(m)),das(itrrv(m)),
     &              itxnm(m),itynm(m),itznm(m),
     &              das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &              das_itzrg(itzrg(m)),
     &              igsh,idasa)

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(1,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(1,trTEMP_sum,ncol,m)
                  endif

               end if

*-----------------------------------------------------------------------

            else if( itals(m) .eq. 44 ) then

                     idas0 = nmmax
                     idas1 = lmmax
                     idas2 = idas1 + itenm(m)
                     idas3 = idas2 + ittnm(m)
                     idasa = idas3 + itanm(m)

                  if( npe .gt. 1 ) then
                     call paratal(tr,trTEMP,mnmax)
                     call paratal_sumover(trTEMP_sum,m)
                  endif

               if( me .eq. 0 ) then

                  !OBINATA(2012.6.26)
                  if( npe .gt. 1 ) then
                    call sumtr(tr,trRES(irestalm(m)),mnmax,ncol)
                    call sumtr_sumover(m,ncol)   !  tr_sum -> trRES_sum
                   endif

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(0,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(0,trTEMP_sum,ncol,m)
                  endif

                  call padjntrz(m,
     &              itpan(m),itrnm(m),itznm(m),itanm(m),
     &              itenm(m),itmst(m),ittnm(m),
     &              das_itrrg(itrrg(m)),das_itzrg(itzrg(m)),
     &              das_itarg(itarg(m)),
     &              das_iterg(iterg(m)),das_ittrg(ittrg(m)),tr,
     &              idasa)

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(1,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(1,trTEMP_sum,ncol,m)
                  endif

               end if

*-----------------------------------------------------------------------

            else if( itals(m) .eq. 45 ) then

                     idas0 = nmmax
                     idas1 = lmmax
                     idas2 = idas1 + itenm(m)
                     idas3 = ( idas2 + ittnm(m) - 1 ) * 2 + 1
                     idasa = idas2 + ittnm(m)

                  if( npe .gt. 1 ) then
                     call paratal(tr,trTEMP,mnmax)
                     call paratal_sumover(trTEMP_sum,m)
                  endif

               if( me .eq. 0 ) then

                  !OBINATA(2012.6.26)
                  if( npe .gt. 1 ) then
                    call sumtr(tr,trRES(irestalm(m)),mnmax,ncol)
                    call sumtr_sumover(m,ncol)   !  tr_sum -> trRES_sum
                   endif

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(0,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(0,trTEMP_sum,ncol,m)
                  endif

                  call padjntxyz(m,
     &              itpan(m),itmtn(m),ismte(itmtt(m)),
     &              itxnm(m),itynm(m),itznm(m),itenm(m),itmst(m),
     &              ittnm(m),
     &              das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &              das_itzrg(itzrg(m)),
     &              das_iterg(iterg(m)),das_ittrg(ittrg(m)),tr,
     &              igsh,idasa)

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(1,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(1,trTEMP_sum,ncol,m)
                  endif

               end if

*-----------------------------------------------------------------------

            else if( itals(m) .eq. 4 ) then

                     idas0 = nmmax
                     idas1 = lmmax
                     idas2 = idas1 + itenm(m)
                     idas3 = idas2 + itanm(m)
                     idasa = idas3 + ittnm(m)

                  if( npe .gt. 1 ) then
                     call paratal(tr,trTEMP,mnmax)
                     call paratal_sumover(trTEMP_sum,m)
                  endif

               if( me .eq. 0 ) then

                  !OBINATA(2012.7.13)
                  if( npe .gt. 1 ) then
                    call sumtr(tr,trRES(irestalm(m)),mnmax,ncol)
                    call sumtr_sumover(m,ncol)   !  tr_sum -> trRES_sum
                   endif

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(0,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(0,trTEMP_sum,ncol,m)
                  endif

                  call psufreg(m,
     &              itpan(m),itrcn(m),itenm(m),itanm(m),ittnm(m),
     &              itmst(m),itrss(m),
     &              idas_itrcc(itrcc(m)),das_itrca(itrca(m)),
     &              das_iterg(iterg(m)),das_itarg(itarg(m)),
     &              das_ittrg(ittrg(m)),
     &              tr,idasa)

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(1,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(1,trTEMP_sum,ncol,m)
                  endif

                  if( ( me .gt. 0 .or. npe .le. 1 ) .and.
     &               ncol .eq. 2 .and. itmdp(m,0) .ne. 0 )
     &               close(itmdf(m))

               end if

*-----------------------------------------------------------------------

            else if( itals(m) .eq. 5 ) then

                     idas0 = nmmax
                     idas1 = lmmax
                     idas2 = idas1 + itenm(m)
                     idas3 = idas2 + itanm(m)
                     idasa = idas3 + ittnm(m)
                     idasz = itenm(m) * itpan(m) * 2
     &                     * ( itrnm(m) + 1 ) * itznm(m)
     &                     * itanm(m)
     &                     * ittnm(m) * itmst(m)

                  if( npe .gt. 1 ) then
                     call paratal(tr,trTEMP,mnmax)
                     call paratal_sumover(trTEMP_sum,m)
                  endif

               if( me .eq. 0 ) then

                  !OBINATA(2012.7.13)
                  if( npe .gt. 1 ) then
                    call sumtr(tr,trRES(irestalm(m)),mnmax,ncol)
                    call sumtr_sumover(m,ncol)   !  tr_sum -> trRES_sum
                   endif

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(0,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(0,trTEMP_sum,ncol,m)
                  endif

                  if( itenclo(m) .eq. 1 ) then

                     call psufrz_rcc(m,
     &                 itpan(m),itrnm(m),itznm(m),itenm(m),itanm(m),
     &                 ittnm(m),itmst(m),
     &                 das_itrrg(itrrg(m)),das_itzrg(itzrg(m)),
     &                 das_iterg(iterg(m)),das_itarg(itarg(m)),
     &                 das_ittrg(ittrg(m)),
     &                 tr,tr(1+idasz),idasa)

                  else
                     call psufrz(m,
     &                 itpan(m),itrnm(m),itznm(m),itenm(m),itanm(m),
     &                 ittnm(m),itmst(m),
     &                 das_itrrg(itrrg(m)),das_itzrg(itzrg(m)),
     &                 das_iterg(iterg(m)),das_itarg(itarg(m)),
     &                 das_ittrg(ittrg(m)),
     &                 tr,tr(1+idasz),idasa)

                  end if

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(1,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(1,trTEMP_sum,ncol,m)
                  endif

               endif

*-----------------------------------------------------------------------

            else if( itals(m) .eq. 12 ) then

                     idas0 = nmmax
                     idas1 = lmmax
                     idas2 = idas1 + itenm(m)
                     idas3 = idas2 + itanm(m)
                     idas4 = ( idas3 + ittnm(m) - 1 ) * 2 + 1
                     idasa = idas3 + ittnm(m)

                  if( npe .gt. 1 ) then
                     call paratal(tr,trTEMP,mnmax)
                     call paratal_sumover(trTEMP_sum,m)
                  endif

               if( me .eq. 0 ) then

                  !OBINATA(2012.7.13)
                  if( npe .gt. 1 ) then
                    call sumtr(tr,trRES(irestalm(m)),mnmax,ncol)
                    call sumtr_sumover(m,ncol)   !  tr_sum -> trRES_sum
                   endif

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(0,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(0,trTEMP_sum,ncol,m)
                  endif

                  if( itenclo(m) .eq. 1 ) then

                     call psufxyz_rpp(m,
     &                 itpan(m),itxnm(m),itynm(m),itznm(m),
     &                 itenm(m),itanm(m),ittnm(m),itmst(m),
     &                 das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &                 das_itzrg(itzrg(m)),
     &                 das_iterg(iterg(m)),das_itarg(itarg(m)),
     &                 das_ittrg(ittrg(m)),
     &                 tr,igsh,idasa)

                  else
                     call psufxyz(m,
     &                 itpan(m),itxnm(m),itynm(m),itznm(m),
     &                 itenm(m),itanm(m),ittnm(m),itmst(m),
     &                 das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &                 das_itzrg(itzrg(m)),
     &                 das_iterg(iterg(m)),das_itarg(itarg(m)),
     &                 das_ittrg(ittrg(m)),
     &                 tr,igsh,idasa)

                  end if

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(1,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(1,trTEMP_sum,ncol,m)
                  endif

               endif

*-----------------------------------------------------------------------

            else if( itals(m) .eq. 6 ) then

                     idas0 = nmmax
                     idas1 = lmmax
                     idas2 = idas1 + ( maxnt + maxpt ) * 2
                     idas3 = ( idas2 + itrgn(m) - 1 ) * 2 + 1
                     idas4 = idas2 + itrgn(m)
     &                     + ( itrgn(m) + mod(itrgn(m),2) ) / 2
                     idas5 = ( idas4 + itrgn(m) - 1 ) * 2 + 1
                     idasa = idas4

                  if( npe .gt. 1 ) then
                     call paratal(tr,trTEMP,mnmax)
                     call paratal_sumover(trTEMP_sum,m)
                  endif
                  if( itnzn(m).ne.0 .and. npe.gt.1 ) call paraiznm(m,2)

               if( me .eq. 0 ) then

                  !OBINATA(2012.8.23)
                  if( npe .gt. 1 ) then
                    call sumtr(tr,trRES(irestalm(m)),mnmax,ncol)
                    call sumtr_sumover(m,ncol)   !  tr_sum -> trRES_sum
                   endif

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(0,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(0,trTEMP_sum,ncol,m)
                  endif

                  call pyildreg(m,itndz(m),itndn(m),itndm(m),
     &              itrgn(m),itrgm(m),itnun(m),
     &              idas_itreg(itreg(m)),
     &              isnuc(itnuc(m)),idas(itnkz(m)),idas(itnkn(m)),
     &              tr,itrnv(m),idas(itriv(m)),das(itrrv(m)),
     &              itxnm(m),itynm(m),itznm(m),
     &              das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &              das_itzrg(itzrg(m)),
     &              igsh,idasa)

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(1,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(1,trTEMP_sum,ncol,m)
                  endif

               end if

*-----------------------------------------------------------------------

            else if( itals(m) .eq. 7 ) then

                     idas0 = nmmax
                     idas1 = lmmax
                     idas2 = idas1 + ( maxpt + maxnt ) * 2
                     idasa = idas2 + itnfr(m) * itnfz(m)

                  if( npe .gt. 1 ) then
                     call paratal(tr,trTEMP,mnmax)
                     call paratal_sumover(trTEMP_sum,m)
                  endif
                  if( itnzn(m).ne.0 .and. npe.gt.1 ) call paraiznm(m,2)

               if( me .eq. 0 ) then

                  !OBINATA(2012.8.23)
                  if( npe .gt. 1 ) then
                    call sumtr(tr,trRES(irestalm(m)),mnmax,ncol)
                    call sumtr_sumover(m,ncol)   !  tr_sum -> trRES_sum
                   endif

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(0,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(0,trTEMP_sum,ncol,m)
                  endif

                  call pyildrz(m,itndz(m),itndn(m),itndm(m),
     &              itnfr(m),itnfz(m),
     &              itrnm(m),itznm(m),itnun(m),
     &              das_itrrg(itrrg(m)),das_itzrg(itzrg(m)),
     &              isnuc(itnuc(m)),idas(itnkz(m)),idas(itnkn(m)),
     &              tr,idasa)

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(1,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(1,trTEMP_sum,ncol,m)
                  endif

               end if

*-----------------------------------------------------------------------

            else if( itals(m) .eq. 8 ) then

                     idas0 = nmmax
                     idas1 = lmmax
                     idas2 = idas1 + ( maxpt + maxnt ) * 2
                     idas3 = ( idas2 + itnfn(m) * itnfn(m) - 1 ) * 2 + 1
                     idasa = idas2 + itnfn(m) * itnfn(m)

                  if( npe .gt. 1 ) then
                     call paratal(tr,trTEMP,mnmax)
                     call paratal_sumover(trTEMP_sum,m)
                  endif
                  if( itnzn(m).ne.0 .and. npe.gt.1 ) call paraiznm(m,2)

               if( me .eq. 0 ) then

                  !OBINATA(2012.8.23)
                  if( npe .gt. 1 ) then
                    call sumtr(tr,trRES(irestalm(m)),mnmax,ncol)
                    call sumtr_sumover(m,ncol)   !  tr_sum -> trRES_sum
                   endif

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(0,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(0,trTEMP_sum,ncol,m)
                  endif

                  call pyildxyz(m,itndz(m),itndn(m),itndm(m),itnfn(m),
     &              itmtn(m),ismte(itmtt(m)),
     &              itxnm(m),itynm(m),itznm(m),itnun(m),
     &              das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &              das_itzrg(itzrg(m)),
     &              isnuc(itnuc(m)),idas(itnkz(m)),idas(itnkn(m)),
     &              tr,igsh,idasa)

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(1,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(1,trTEMP_sum,ncol,m)
                  endif

               end if

*-----------------------------------------------------------------------
            else if( itals(m) .eq. 51 ) then

                     idas0 = nmmax
                     idas1 = lmmax
                     idas2 = idas1 + ( maxnt + maxpt ) * 2
                     idas3 = ( idas2 + itrgn(m) - 1 ) * 2 + 1
                     idas4 = idas2 + itrgn(m)
     &                     + ( itrgn(m) + mod(itrgn(m),2) ) / 2
                     idas5 = ( idas4 + itrgn(m) - 1 ) * 2 + 1
                     idasa = idas4

                  if( npe .gt. 1 ) then
                     call paratal(tr,trTEMP,mnmax)
                     call paratal_sumover(trTEMP_sum,m)
                  endif
                  if( itnzn(m).ne.0 .and. npe.gt.1 ) call paraiznm(m,2)

               if( me .eq. 0 ) then

                  !OBINATA(2012.8.23)
                  if( npe .gt. 1 ) then
                    call sumtr(tr,trRES(irestalm(m)),mnmax,ncol)
                    call sumtr_sumover(m,ncol)   !  tr_sum -> trRES_sum
                   endif

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(0,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(0,trTEMP_sum,ncol,m)
                  endif

                  call pyildtet(m,itndz(m),itndn(m),itndm(m),
     &              itrgn(m),itrgm(m),itnun(m),
     &              isnuc(itnuc(m)),idas(itnkz(m)),idas(itnkn(m)),
     &              tr,
     &              itxnm(m),itynm(m),itznm(m),
     &              idas_itreg(itreg(m)),
     &              das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &              das_itzrg(itzrg(m)),
     &              igsh,idasa)

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(1,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(1,trTEMP_sum,ncol,m)
                  endif

               end if

*-----------------------------------------------------------------------

            else if( itals(m) .eq. 9 ) then

                     idas0 = nmmax
                     idas1 = lmmax
                     idas2 = ( idas1 + itrgn(m) - 1 ) * 2 + 1
                     idas3 = idas1 + itrgn(m)
     &                     + ( itrgn(m) + mod(itrgn(m),2) ) / 2
                     idas4 = ( idas3 + itrgn(m) - 1 ) * 2 + 1
                     idasa = idas3

                  if( npe .gt. 1 ) then
                     call paratal(tr,trTEMP,mnmax)
                     call paratal_sumover(trTEMP_sum,m)
                  endif

               if( me .eq. 0 ) then

                  !OBINATA(2012.7.23)
                  if( npe .gt. 1 ) then
                    call sumtr(tr,trRES(irestalm(m)),mnmax,ncol)
                    call sumtr_sumover(m,ncol)   !  tr_sum -> trRES_sum
                   endif

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(0,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(0,trTEMP_sum,ncol,m)
                  endif

                  call phetreg(m,itndy(m),itout(m),
     &              itpan(m),itrgn(m),itrgm(m),
     &              idas_itreg(itreg(m)),itenm(m),das_iterg(iterg(m)),
     &              tr,das(ithet(m)),
     &              itrnv(m),idas(itriv(m)),das(itrrv(m)),
     &              itxnm(m),itynm(m),itznm(m),
     &              das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &              das_itzrg(itzrg(m)),
     &              igsh,idasa)

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(1,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(1,trTEMP_sum,ncol,m)
                  endif

               end if

*-----------------------------------------------------------------------

            else if( itals(m) .eq. 10 ) then

                     idas0 = nmmax
                     idasa = lmmax

                  if( npe .gt. 1 ) then
                     call paratal(tr,trTEMP,mnmax)
                     call paratal_sumover(trTEMP_sum,m)
                  endif

               if( me .eq. 0 ) then

                  !OBINATA(2012.7.23)
                  if( npe .gt. 1 ) then
                    call sumtr(tr,trRES(irestalm(m)),mnmax,ncol)
                    call sumtr_sumover(m,ncol)   !  tr_sum -> trRES_sum
                   endif

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(0,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(0,trTEMP_sum,ncol,m)
                  endif

                  call phetrz(m,itndy(m),itout(m),
     &              itpan(m),itrnm(m),itznm(m),
     &              das_itrrg(itrrg(m)),das_itzrg(itzrg(m)),
     &              itenm(m),das_iterg(iterg(m)),
     &              tr,das(ithet(m)),idasa)

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(1,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(1,trTEMP_sum,ncol,m)
                  endif

               end if

*-----------------------------------------------------------------------

            else if( itals(m) .eq. 11 ) then

                     idas0 = nmmax
                     idas1 = ( lmmax - 1 ) * 2 + 1
                     idasa = lmmax

                  if( npe .gt. 1 ) then
                     call paratal(tr,trTEMP,mnmax)
                     call paratal_sumover(trTEMP_sum,m)
                  endif

               if( me .eq. 0 ) then

                  !OBINATA(2012.7.23)
                  if( npe .gt. 1 ) then
                    call sumtr(tr,trRES(irestalm(m)),mnmax,ncol)
                    call sumtr_sumover(m,ncol)   !  tr_sum -> trRES_sum
                   endif

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(0,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(0,trTEMP_sum,ncol,m)
                  endif

                  call phetxyz(m,itndy(m),itout(m),
     &              itpan(m),itmtn(m),ismte(itmtt(m)),
     &              itxnm(m),itynm(m),itznm(m),
     &              das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &              das_itzrg(itzrg(m)),
     &              itenm(m),das_iterg(iterg(m)),
     &              tr,das(ithet(m)),
     &              igsh,idasa)

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(1,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(1,trTEMP_sum,ncol,m)
                  endif

               end if

*-----------------------------------------------------------------------

            else if( itals(m) .eq. 13 ) then

                     idas0 = nmmax
                     idas1 = lmmax
                     idas2 = idas1 + itenm(m)
                     idas3 = idas2 + ittnm(m)
                     idas4 = idas3 + itactnm(m)   ! S.Abe 2018/02/15
                     idas5 = ( idas4 + itrgn(m) - 1 ) * 2 + 1
                     idas6 = idas4 + itrgn(m)
     &                     + ( itrgn(m) + mod(itrgn(m),2) ) / 2
                     idas7 = ( idas6 + itrgn(m) - 1 ) * 2 + 1
                     idasa = idas6

                  if( npe .gt. 1 ) then
                     call paratal(tr,trTEMP,mnmax)
                     call paratal_sumover(trTEMP_sum,m)
                  endif

               if( me .eq. 0 ) then

                  !OBINATA(2012.7.25)
                  if( npe .gt. 1 ) then
                    call sumtr(tr,trRES(irestalm(m)),mnmax,ncol)
                    call sumtr_sumover(m,ncol)   !  tr_sum -> trRES_sum
                   endif

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(0,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(0,trTEMP_sum,ncol,m)
                  endif

                  call pstarreg(m,
     &              itpan(m),itrgn(m),itrgm(m),itenm(m),ittnm(m),
     &              itactnm(m),                                         ! S.Abe 2018/02/15
     &              idas_itreg(itreg(m)),
     &              das_iterg(iterg(m)),das_ittrg(ittrg(m)),
     &              das(itactrg(m)),tr,     ! S.Abe 2018/02/15, add das(itactrg(m))
     &              itrnv(m),idas(itriv(m)),das(itrrv(m)),
     &              itxnm(m),itynm(m),itznm(m),
     &              das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &              das_itzrg(itzrg(m)),
     &              igsh,idasa)

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(1,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(1,trTEMP_sum,ncol,m)
                  endif

               end if

*-----------------------------------------------------------------------

            else if( itals(m) .eq. 14 ) then

                     idas0 = nmmax
                     idas1 = lmmax
                     idas2 = idas1 + itenm(m)
                     idasa = idas2 + ittnm(m)

                  if( npe .gt. 1 ) then
                     call paratal(tr,trTEMP,mnmax)
                     call paratal_sumover(trTEMP_sum,m)
                  endif

               if( me .eq. 0 ) then

                  !OBINATA(2012.7.25)
                  if( npe .gt. 1 ) then
                    call sumtr(tr,trRES(irestalm(m)),mnmax,ncol)
                    call sumtr_sumover(m,ncol)   !  tr_sum -> trRES_sum
                   endif

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(0,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(0,trTEMP_sum,ncol,m)
                  endif

                  call pstarrz(m,
     &              itpan(m),itrnm(m),itznm(m),itenm(m),ittnm(m),
     &              das_itrrg(itrrg(m)),das_itzrg(itzrg(m)),
     &              das_iterg(iterg(m)),das_ittrg(ittrg(m)),tr,
     &              idasa)

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(1,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(1,trTEMP_sum,ncol,m)
                  endif

               end if

*-----------------------------------------------------------------------

            else if( itals(m) .eq. 15 ) then

                     idas0 = nmmax
                     idas1 = lmmax
                     idas2 = idas1 + itenm(m)
                     idas3 = ( idas1 + ittnm(m) - 1 ) * 2 + 1
                     idasa = idas2 + ittnm(m)

                  if( npe .gt. 1 ) then
                     call paratal(tr,trTEMP,mnmax)
                     call paratal_sumover(trTEMP_sum,m)
                  endif

               if( me .eq. 0 ) then

                  !OBINATA(2012.7.25)
                  if( npe .gt. 1 ) then
                    call sumtr(tr,trRES(irestalm(m)),mnmax,ncol)
                    call sumtr_sumover(m,ncol)   !  tr_sum -> trRES_sum
                   endif

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(0,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(0,trTEMP_sum,ncol,m)
                  endif

                  call pstarxyz(m,
     &              itpan(m),itmtn(m),ismte(itmtt(m)),
     &              itxnm(m),itynm(m),itznm(m),itenm(m),
     &              ittnm(m),
     &              das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &              das_itzrg(itzrg(m)),
     &              das_iterg(iterg(m)),das_ittrg(ittrg(m)),tr,
     &              igsh,idasa)

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(1,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(1,trTEMP_sum,ncol,m)
                  endif

               end if

*-----------------------------------------------------------------------

            else if( itals(m) .eq. 16 ) then

                     idas0 = nmmax
                     idas1 = lmmax
                     idas2 = idas1 + itenm(m)
                     idas3 = idas2 + ittnm(m)
                     idas4 = ( idas3 + itrgn(m) - 1 ) * 2 + 1
                     idas5 = idas3 + itrgn(m)
     &                     + ( itrgn(m) + mod(itrgn(m),2) ) / 2
                     idas6 = ( idas5 + itrgn(m) - 1 ) * 2 + 1
                     idasa = idas5

                  if( npe .gt. 1 ) then
                     call paratal(tr,trTEMP,mnmax)
                     call paratal_sumover(trTEMP_sum,m)
                  endif

               if( me .eq. 0 ) then

                  !OBINATA(2012.9.19)
                  if( npe .gt. 1 ) then
                    call sumtr(tr,trRES(irestalm(m)),mnmax,ncol)
                    call sumtr_sumover(m,ncol)   !  tr_sum -> trRES_sum
                   endif

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(0,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(0,trTEMP_sum,ncol,m)
                  endif

                  call ptimereg(m,
     &              itpan(m),itrgn(m),itrgm(m),itenm(m),ittnm(m),
     &              idas_itreg(itreg(m)),
     &              das_iterg(iterg(m)),das_ittrg(ittrg(m)),tr,
     &              itrnv(m),idas(itriv(m)),das(itrrv(m)),
     &              itxnm(m),itynm(m),itznm(m),
     &              das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &              das_itzrg(itzrg(m)),
     &              igsh,idasa)

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(1,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(1,trTEMP_sum,ncol,m)
                  endif

                  if( ( me .gt. 0 .or. npe .le. 1 ) .and.
     &               ncol .eq. 2 .and. itmdp(m,0) .ne. 0 )
     &               close(itmdf(m))

               end if

*-----------------------------------------------------------------------

            else if( itals(m) .eq. 17 ) then

                     idas0 = nmmax
                     idas1 = lmmax
                     idas2 = idas1 + itenm(m)
                     idasa = idas2 + ittnm(m)

                  if( npe .gt. 1 ) then
                     call paratal(tr,trTEMP,mnmax)
                     call paratal_sumover(trTEMP_sum,m)
                  endif

               if( me .eq. 0 ) then

                  !OBINATA(2012.9.19)
                  if( npe .gt. 1 ) then
                    call sumtr(tr,trRES(irestalm(m)),mnmax,ncol)
                    call sumtr_sumover(m,ncol)   !  tr_sum -> trRES_sum
                   endif

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(0,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(0,trTEMP_sum,ncol,m)
                  endif

                  call ptimerz(m,
     &              itpan(m),itrnm(m),itznm(m),itenm(m),ittnm(m),
     &              das_itrrg(itrrg(m)),das_itzrg(itzrg(m)),
     &              das_iterg(iterg(m)),das_ittrg(ittrg(m)),tr,
     &              idasa)

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(1,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(1,trTEMP_sum,ncol,m)
                  endif

               end if

*-----------------------------------------------------------------------

            else if( itals(m) .eq. 18 ) then

                     idas0 = nmmax
                     idas1 = lmmax
                     idas2 = idas1 + itenm(m)
                     idas3 = ( idas2 + ittnm(m) - 1 ) * 2 + 1
                     idasa = idas2 + ittnm(m)

                  if( npe .gt. 1 ) then
                     call paratal(tr,trTEMP,mnmax)
                     call paratal_sumover(trTEMP_sum,m)
                  endif

               if( me .eq. 0 ) then

                  !OBINATA(2012.9.19)
                  if( npe .gt. 1 ) then
                    call sumtr(tr,trRES(irestalm(m)),mnmax,ncol)
                    call sumtr_sumover(m,ncol)   !  tr_sum -> trRES_sum
                   endif

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(0,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(0,trTEMP_sum,ncol,m)
                  endif

                  call ptimexyz(m,
     &              itmtn(m),ismte(itmtt(m)),
     &              itpan(m),itxnm(m),itynm(m),itznm(m),itenm(m),
     &              ittnm(m),
     &              das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &              das_itzrg(itzrg(m)),
     &              das_iterg(iterg(m)),das_ittrg(ittrg(m)),tr,
     &              igsh,idasa)

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(1,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(1,trTEMP_sum,ncol,m)
                  endif

               end if

*-----------------------------------------------------------------------

            else if( itals(m) .eq. 19 ) then

                     idas0 = nmmax
                     idas1 = lmmax
                     idas2 = ( idas1 + itrgn(m) - 1 ) * 2 + 1
                     idas3 = idas1 + itrgn(m)
     &                     + ( itrgn(m) + mod(itrgn(m),2) ) / 2
                     idas4 = ( idas3 + itrgn(m) - 1 ) * 2 + 1
                     idasa = idas3

                  if( npe .gt. 1 ) then
                     call paratal(tr,trTEMP,mnmax)
                     call paratal_sumover(trTEMP_sum,m)
                  endif

               if( me .eq. 0 ) then

                  !OBINATA(2012.9.14)
                  if( npe .gt. 1 ) then
                    call sumtr(tr,trRES(irestalm(m)),mnmax,ncol)
                    call sumtr_sumover(m,ncol)   !  tr_sum -> trRES_sum
                   endif

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(0,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(0,trTEMP_sum,ncol,m)
                  endif

                  call pdpareg(m,itndy(m),itout(m),
     &              itpan(m),itrgn(m),itrgm(m),idas_itreg(itreg(m)),
     &              tr,
     &              itrnv(m),idas(itriv(m)),das(itrrv(m)),
     &              itxnm(m),itynm(m),itznm(m),
     &              das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &              das_itzrg(itzrg(m)),
     &              igsh,idasa)

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(1,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(1,trTEMP_sum,ncol,m)
                  endif

               end if

*-----------------------------------------------------------------------

            else if( itals(m) .eq. 20 ) then

                     idas0 = nmmax
                     idasa = lmmax

                  if( npe .gt. 1 ) then
                     call paratal(tr,trTEMP,mnmax)
                     call paratal_sumover(trTEMP_sum,m)
                  endif

               if( me .eq. 0 ) then

                  !OBINATA(2012.9.14)
                  if( npe .gt. 1 ) then
                    call sumtr(tr,trRES(irestalm(m)),mnmax,ncol)
                    call sumtr_sumover(m,ncol)   !  tr_sum -> trRES_sum
                   endif

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(0,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(0,trTEMP_sum,ncol,m)
                  endif

                  call pdparz(m,itndy(m),itout(m),
     &              itpan(m),itrnm(m),itznm(m),
     &              das_itrrg(itrrg(m)),das_itzrg(itzrg(m)),
     &              tr,idasa)

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(1,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(1,trTEMP_sum,ncol,m)
                  endif

               end if

*-----------------------------------------------------------------------

            else if( itals(m) .eq. 21 ) then

                     idas0 = nmmax
                     idas1 = ( lmmax - 1 ) * 2 + 1
                     idasa = lmmax

                  if( npe .gt. 1 ) then
                     call paratal(tr,trTEMP,mnmax)
                     call paratal_sumover(trTEMP_sum,m)
                  endif

               if( me .eq. 0 ) then

                  !OBINATA(2012.9.14)
                  if( npe .gt. 1 ) then
                    call sumtr(tr,trRES(irestalm(m)),mnmax,ncol)
                    call sumtr_sumover(m,ncol)   !  tr_sum -> trRES_sum
                   endif

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(0,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(0,trTEMP_sum,ncol,m)
                  endif


                  call pdpaxyz(m,itndy(m),itout(m),
     &              itpan(m),itmtn(m),ismte(itmtt(m)),
     &              itxnm(m),itynm(m),itznm(m),
     &              das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &              das_itzrg(itzrg(m)),
     &              tr,igsh,idasa)

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(1,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(1,trTEMP_sum,ncol,m)
                  endif

               end if

*-----------------------------------------------------------------------
            else if( itals(m) .eq. 53 ) then

                     idas0 = nmmax
                     idas1 = lmmax
                     idas2 = ( idas1 + itrgn(m) - 1 ) * 2 + 1
                     idas3 = idas1 + itrgn(m)
     &                     + ( itrgn(m) + mod(itrgn(m),2) ) / 2
                     idas4 = ( idas3 + itrgn(m) - 1 ) * 2 + 1
                     idasa = idas3

                  if( npe .gt. 1 ) then
                     call paratal(tr,trTEMP,mnmax)
                     call paratal_sumover(trTEMP_sum,m)
                  endif

               if( me .eq. 0 ) then

                  !OBINATA(2012.9.14)
                  if( npe .gt. 1 ) then
                    call sumtr(tr,trRES(irestalm(m)),mnmax,ncol)
                    call sumtr_sumover(m,ncol)   !  tr_sum -> trRES_sum
                   endif

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(0,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(0,trTEMP_sum,ncol,m)
                  endif

                  call pdpatet(m,itndy(m),itout(m),
     &              itpan(m),itrgn(m),itrgm(m),
     &              tr,
     &              itxnm(m),itynm(m),itznm(m),
     &              idas_itreg(itreg(m)),   
     &              das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &              das_itzrg(itzrg(m)),
     &              igsh,idasa)

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(1,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(1,trTEMP_sum,ncol,m)
                  endif

               end if

*-----------------------------------------------------------------------

            else if( itals(m) .eq. 22 ) then

                     idas0 = nmmax
                     idas1 = lmmax
                     idas2 = idas1 + itenm(m)
                     idas3 = idas2 + ittnm(m)
                     idas4 = idas3 + itanm(m)
                     idas5 = ( idas4 + itrgn(m) - 1 ) * 2 + 1
                     idas6 = idas4 + itrgn(m)
     &                     + ( itrgn(m) + mod(itrgn(m),2) ) / 2
                     idas7 = ( idas6 + itrgn(m) - 1 ) * 2 + 1
                     idasa = idas6

                  if( npe .gt. 1 ) then
                     call paratal(tr,trTEMP,mnmax)
                     call paratal_sumover(trTEMP_sum,m)
                  endif

               if( me .eq. 0 ) then

                  !OBINATA(2012.9.11)
                  if( npe .gt. 1 ) then
                    call sumtr(tr,trRES(irestalm(m)),mnmax,ncol)
                    call sumtr_sumover(m,ncol)   !  tr_sum -> trRES_sum
                   endif

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(0,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(0,trTEMP_sum,ncol,m)
                  endif

                  call ppdctreg(m,
     &              itpan(m),itrgn(m),itrgm(m),itenm(m),ittnm(m),
     &              itanm(m),idas_itreg(itreg(m)),
     &              das_iterg(iterg(m)),das_ittrg(ittrg(m)),
     &              das_itarg(itarg(m)),
     &              tr,
     &              itrnv(m),idas(itriv(m)),das(itrrv(m)),
     &              itxnm(m),itynm(m),itznm(m),
     &              das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &              das_itzrg(itzrg(m)),
     &              igsh,idasa)

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(1,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(1,trTEMP_sum,ncol,m)
                  endif

                  if( ( me .gt. 0 .or. npe .le. 1 ) .and.
     &               ncol .eq. 2 .and. itmdp(m,0) .ne. 0 )
     &               close(itmdf(m))

               end if

*-----------------------------------------------------------------------

            else if( itals(m) .eq. 23 ) then

                     idas0 = nmmax
                     idas1 = lmmax
                     idas2 = idas1 + itenm(m)
                     idas3 = idas2 + ittnm(m)
                     idasa = idas3 + itanm(m)

                  if( npe .gt. 1 ) then
                     call paratal(tr,trTEMP,mnmax)
                     call paratal_sumover(trTEMP_sum,m)
                  endif

               if( me .eq. 0 ) then

                  !OBINATA(2012.9.11)
                  if( npe .gt. 1 ) then
                    call sumtr(tr,trRES(irestalm(m)),mnmax,ncol)
                    call sumtr_sumover(m,ncol)   !  tr_sum -> trRES_sum
                   endif

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(0,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(0,trTEMP_sum,ncol,m)
                  endif

                  call ppdctrz(m,
     &              itpan(m),itrnm(m),itznm(m),itenm(m),ittnm(m),
     &              itanm(m),das_itrrg(itrrg(m)),das_itzrg(itzrg(m)),
     &              das_iterg(iterg(m)),das_ittrg(ittrg(m)),
     &              das_itarg(itarg(m)),
     &              tr,idasa)

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(1,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(1,trTEMP_sum,ncol,m)
                  endif

               end if

*-----------------------------------------------------------------------

            else if( itals(m) .eq. 24 ) then

                     idas0 = nmmax
                     idas1 = lmmax
                     idas2 = idas1 + itenm(m)
                     idas3 = idas2 + ittnm(m)
                     idas4 = ( idas3 + itanm(m) - 1 ) * 2 + 1
                     idasa = idas3 + itanm(m)

                  if( npe .gt. 1 ) then
                     call paratal(tr,trTEMP,mnmax)
                     call paratal_sumover(trTEMP_sum,m)
                  endif

               if( me .eq. 0 ) then

                  !OBINATA(2012.9.11)
                  if( npe .gt. 1 ) then
                    call sumtr(tr,trRES(irestalm(m)),mnmax,ncol)
                    call sumtr_sumover(m,ncol)   !  tr_sum -> trRES_sum
                   endif

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(0,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(0,trTEMP_sum,ncol,m)
                  endif

                  call ppdctxyz(m,
     &              itpan(m),itmtn(m),ismte(itmtt(m)),
     &              itxnm(m),itynm(m),itznm(m),itenm(m),
     &              ittnm(m),itanm(m),
     &              das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &              das_itzrg(itzrg(m)),
     &              das_iterg(iterg(m)),das_ittrg(ittrg(m)),
     &              das_itarg(itarg(m)),
     &              tr,igsh,idasa)

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(1,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(1,trTEMP_sum,ncol,m)
                  endif

               end if

*-----------------------------------------------------------------------
            else if( itals(m) .eq. 52 ) then

                     idas0 = nmmax
                     idas1 = lmmax
                     idas2 = idas1 + itenm(m)
                     idas3 = idas2 + ittnm(m)
                     idas4 = idas3 + itanm(m)
                     idas5 = ( idas4 + itrgn(m) - 1 ) * 2 + 1
                     idas6 = idas4 + itrgn(m)
     &                     + ( itrgn(m) + mod(itrgn(m),2) ) / 2
                     idas7 = ( idas6 + itrgn(m) - 1 ) * 2 + 1
                     idasa = idas6

                  if( npe .gt. 1 ) then
                     call paratal(tr,trTEMP,mnmax)
                     call paratal_sumover(trTEMP_sum,m)
                  endif

               if( me .eq. 0 ) then

                  !OBINATA(2012.9.11)
                  if( npe .gt. 1 ) then
                    call sumtr(tr,trRES(irestalm(m)),mnmax,ncol)
                    call sumtr_sumover(m,ncol)   !  tr_sum -> trRES_sum
                   endif

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(0,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(0,trTEMP_sum,ncol,m)
                  endif

                  call ppdcttet(m,
     &              itpan(m),itrgn(m),itrgm(m),
     &              itenm(m),ittnm(m),itanm(m),
     &              das_iterg(iterg(m)),das_ittrg(ittrg(m)),
     &              das_itarg(itarg(m)),
     &              tr,
     &              itxnm(m),itynm(m),itznm(m),
     &              idas_itreg(itreg(m)),
     &              das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &              das_itzrg(itzrg(m)),
     &              igsh,idasa)

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(1,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(1,trTEMP_sum,ncol,m)
                  endif

                  if( ( me .gt. 0 .or. npe .le. 1 ) .and.
     &               ncol .eq. 2 .and. itmdp(m,0) .ne. 0 )
     &               close(itmdf(m))

               end if

*-----------------------------------------------------------------------

            else if( itals(m) .eq. 25 .and.
     &                ncol .eq. 2 ) then

               if( me .eq. 0 ) then

                     idas1 = ( lmmax - 1 ) * 2 + 1
                     idasa = lmmax

                  call pgshxyz(m,
     &              itxnm(m),itynm(m),itznm(m),
     &              das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &              das_itzrg(itzrg(m)),
     &              idasa)

               end if

*-----------------------------------------------------------------------

            else if( itals(m) .eq. 28 ) then

                     idas0 = nmmax
                     idas1 = lmmax
                     idas2 = idas1 + itenm(m)
                     idas3 = ( idas2 + itrgn(m) - 1 ) * 2 + 1
                     idas4 = idas2 + itrgn(m)
     &                     + ( itrgn(m) + mod(itrgn(m),2) ) / 2
                     idas5 = ( idas4 + itrgn(m) - 1 ) * 2 + 1
                     idasa = idas4

                  if( npe .gt. 1 ) then
                     call paratal(tr,trTEMP,mnmax)
                     call paratal_sumover(trTEMP_sum,m)
                  endif

               if( me .eq. 0 ) then

                  !OBINATA(2012.9.5)
                  if( npe .gt. 1 ) then
                    call sumtr(tr,trRES(irestalm(m)),mnmax,ncol)
                    call sumtr_sumover(m,ncol)   !  tr_sum -> trRES_sum
                   endif

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(0,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(0,trTEMP_sum,ncol,m)
                  endif

                  call pletreg(m,
     &              itpan(m),itrgn(m),itrgm(m),itenm(m),
     &              idas_itreg(itreg(m)),
     &              das_iterg(iterg(m)),tr,
     &              itrnv(m),idas(itriv(m)),das(itrrv(m)),
     &              itxnm(m),itynm(m),itznm(m),
     &              das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &              das_itzrg(itzrg(m)),
     &              igsh,idasa)

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(1,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(1,trTEMP_sum,ncol,m)
                  endif

               end if

*-----------------------------------------------------------------------

            else if( itals(m) .eq. 29 ) then

                     idas0 = nmmax
                     idas1 = lmmax
                     idasa = idas1 + itenm(m)

                  if( npe .gt. 1 ) then
                     call paratal(tr,trTEMP,mnmax)
                     call paratal_sumover(trTEMP_sum,m)
                  endif

               if( me .eq. 0 ) then

                  !OBINATA(2012.9.5)
                  if( npe .gt. 1 ) then
                    call sumtr(tr,trRES(irestalm(m)),mnmax,ncol)
                    call sumtr_sumover(m,ncol)   !  tr_sum -> trRES_sum
                   endif

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(0,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(0,trTEMP_sum,ncol,m)
                  endif

                  call pletrz(m,
     &              itpan(m),itrnm(m),itznm(m),itenm(m),
     &              das_itrrg(itrrg(m)),das_itzrg(itzrg(m)),
     &              das_iterg(iterg(m)),tr,idasa)

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(1,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(1,trTEMP_sum,ncol,m)
                  endif

               end if

*-----------------------------------------------------------------------

            else if( itals(m) .eq. 30 ) then

                     idas0 = nmmax
                     idas1 = lmmax
                     idas2 = ( idas1 + itenm(m) - 1 ) * 2 + 1
                     idasa = idas1 + itenm(m)

                  if( npe .gt. 1 ) then
                     call paratal(tr,trTEMP,mnmax)
                     call paratal_sumover(trTEMP_sum,m)
                  endif

               if( me .eq. 0 ) then

                  !OBINATA(2012.9.5)
                  if( npe .gt. 1 ) then
                    call sumtr(tr,trRES(irestalm(m)),mnmax,ncol)
                    call sumtr_sumover(m,ncol)   !  tr_sum -> trRES_sum
                   endif

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(0,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(0,trTEMP_sum,ncol,m)
                  endif

                  call pletxyz(m,
     &              itpan(m),itmtn(m),ismte(itmtt(m)),
     &              itxnm(m),itynm(m),itznm(m),itenm(m),
     &              das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &              das_itzrg(itzrg(m)),
     &              das_iterg(iterg(m)),tr,igsh,idasa)

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(1,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(1,trTEMP_sum,ncol,m)
                  endif

               end if

*-----------------------------------------------------------------------

            else if( itals(m) .eq. 31 ) then

                     idas0 = nmmax
                     idas1 = lmmax
                     idas2 = idas1 + ittnm(m)

                     idas3 = ( idas2 + itrgn(m) - 1 ) * 2 + 1
                     idas4 = idas2 + itrgn(m)
     &                     + ( itrgn(m) + mod(itrgn(m),2) ) / 2
                     idas5 = ( idas4 + itrgn(m) - 1 ) * 2 + 1
                     idasa = idas3

                  if( npe .gt. 1 ) then
                     call paratal(tr,trTEMP,mnmax)
                     call paratal_sumover(trTEMP_sum,m)
                  endif

               if( me .eq. 0 ) then

                  !OBINATA(2012.8.6)
                  if( npe .gt. 1 ) then
                    call sumtr(tr,trRES(irestalm(m)),mnmax,ncol)
                    call sumtr_sumover(m,ncol)   !  tr_sum -> trRES_sum
                   endif

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(0,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(0,trTEMP_sum,ncol,m)
                  endif

                  call pdepstreg(m,
     &              itpan(m),itrgn(m),itrgm(m),itenm(m),
     &              ittnm(m),
     &              idas_itreg(itreg(m)),
     &              das_iterg(iterg(m)),das_ittrg(ittrg(m)),tr,
     &              itrnv(m),idas(itriv(m)),das(itrrv(m)),
     &              itxnm(m),itynm(m),itznm(m),
     &              das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &              das_itzrg(itzrg(m)),
     &              igsh,idasa)

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(1,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(1,trTEMP_sum,ncol,m)
                  endif

               end if

*-----------------------------------------------------------------------

            else if( itals(m) .eq. 32 ) then

                     idas0 = nmmax
                     idas1 = lmmax
                     idasa = idas1 + ittnm(m)

                  if( npe .gt. 1 ) then
                     call paratal(tr,trTEMP,mnmax)
                     call paratal_sumover(trTEMP_sum,m)
                  endif

               if( me .eq. 0 ) then

                  !OBINATA(2012.8.6)
                  if( npe .gt. 1 ) then
                    call sumtr(tr,trRES(irestalm(m)),mnmax,ncol)
                    call sumtr_sumover(m,ncol)   !  tr_sum -> trRES_sum
                   endif

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(0,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(0,trTEMP_sum,ncol,m)
                  endif

                  call pdepstrz(m,
     &              itpan(m),itrnm(m),itznm(m),itenm(m),
     &              ittnm(m),
     &              das_itrrg(itrrg(m)),das_itzrg(itzrg(m)),
     &              das_iterg(iterg(m)),das_ittrg(ittrg(m)),tr,
     &              idasa)

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(1,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(1,trTEMP_sum,ncol,m)
                  endif

               end if

*-----------------------------------------------------------------------

            else if( itals(m) .eq. 33 ) then

                     idas0 = nmmax
                     idas1 = lmmax
                     idas2 = ( idas1 + ittnm(m) - 1 ) * 2 + 1
                     idasa = idas1 + ittnm(m)

                  if( npe .gt. 1 ) then
                     call paratal(tr,trTEMP,mnmax)
                     call paratal_sumover(trTEMP_sum,m)
                  endif

               if( me .eq. 0 ) then

                  !OBINATA(2012.8.6)
                  if( npe .gt. 1 ) then
                    call sumtr(tr,trRES(irestalm(m)),mnmax,ncol)
                    call sumtr_sumover(m,ncol)   !  tr_sum -> trRES_sum
                   endif

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(0,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(0,trTEMP_sum,ncol,m)
                  endif

                  call pdepstxyz(m,
     &              itpan(m),itmtn(m),ismte(itmtt(m)),
     &              itxnm(m),itynm(m),itznm(m),itenm(m),
     &              ittnm(m),
     &              das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &              das_itzrg(itzrg(m)),
     &              das_iterg(iterg(m)),das_ittrg(ittrg(m)),tr,
     &              igsh,idasa)

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(1,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(1,trTEMP_sum,ncol,m)
                  endif

               end if

*-----------------------------------------------------------------------

            else if( itals(m) .eq. 50 ) then

                     idas0 = nmmax
                     idas1 = lmmax
                     idas2 = idas1 + ittnm(m)

                     idas3 = ( idas2 + itrgn(m) - 1 ) * 2 + 1
                     idas4 = idas2 + itrgn(m)
     &                     + ( itrgn(m) + mod(itrgn(m),2) ) / 2
                     idas5 = ( idas4 + itrgn(m) - 1 ) * 2 + 1
                     idasa = idas3

                  if( npe .gt. 1 ) then
                     call paratal(tr,trTEMP,mnmax)
                     call paratal_sumover(trTEMP_sum,m)
                  endif

               if( me .eq. 0 ) then

                  !OBINATA(2012.8.6)
                  if( npe .gt. 1 ) then
                    call sumtr(tr,trRES(irestalm(m)),mnmax,ncol)
                    call sumtr_sumover(m,ncol)   !  tr_sum -> trRES_sum
                   endif

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(0,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(0,trTEMP_sum,ncol,m)
                  endif

                  call pdepsttet(m,
     &              itpan(m),itrgn(m),itrgm(m),itenm(m),ittnm(m),
     &              das_iterg(iterg(m)),das_ittrg(ittrg(m)),tr,
     &              itxnm(m),itynm(m),itznm(m),
     &              idas_itreg(itreg(m)),
     &              das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &              das_itzrg(itzrg(m)),
     &              igsh,idasa)

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(1,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(1,trTEMP_sum,ncol,m)
                  endif

               end if

*-----------------------------------------------------------------------

            else if( itals(m) .eq. 34 ) then

                     idas0 = nmmax
                     idas1 = lmmax
                     idas2 = idas1 + ittnm(m)
                     idas3 = ( idas2 + itrgn(m) - 1 ) * 2 + 1
                     idasa = idas2 + itrgn(m)
     &                     + ( itrgn(m) + mod(itrgn(m),2) ) / 2

                  if( npe .gt. 1 ) then
                     call paratal(tr,trTEMP,mnmax)
                     call paratal_sumover(trTEMP_sum,m)
                  endif

               if( me .eq. 0 ) then

                  !OBINATA(2012.9.20)
                  if( npe .gt. 1 ) then
                    call sumtr(tr,trRES(irestalm(m)),mnmax,ncol)
                    call sumtr_sumover(m,ncol)   !  tr_sum -> trRES_sum
                   endif

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(0,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(0,trTEMP_sum,ncol,m)
                  endif

                  call pdpst2reg(m,
     &              itpan(m),itrgn(m),itrgm(m),itenm(m),itenm2(m),
     &              ittnm(m),idas_itreg(itreg(m)),
     &              das_iterg(iterg(m)),das_iterg2(iterg2(m)),
     &              das_ittrg(ittrg(m)),
     &              tr,
     &              itrnv(m),idas(itriv(m)),das(itrrv(m)),idasa)

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(1,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(1,trTEMP_sum,ncol,m)
                  endif

               end if

*-----------------------------------------------------------------------

            else if( itals(m) .eq. 35 ) then

                     idas0 = nmmax
                     idas1 = lmmax
                     idas2 = idas1 + itenm(m)
                     idas3 = ( idas2 + itrgn(m) - 1 ) * 2 + 1
                     idas4 = idas2 + itrgn(m)
     &                     + ( itrgn(m) + mod(itrgn(m),2) ) / 2
                     idas5 = ( idas4 + itrgn(m) - 1 ) * 2 + 1
                     idasa = idas4

                  if( npe .gt. 1 ) then
                     call paratal(tr,trTEMP,mnmax)
                     call paratal_sumover(trTEMP_sum,m)
                  endif

               if( me .eq. 0 ) then

                  !OBINATA(2012.9.7)
                  if( npe .gt. 1 ) then
                    call sumtr(tr,trRES(irestalm(m)),mnmax,ncol)
                    call sumtr_sumover(m,ncol)   !  tr_sum -> trRES_sum
                   endif

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(0,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(0,trTEMP_sum,ncol,m)
                  endif

                  call psedreg(m,
     &              itpan(m),itrgn(m),itrgm(m),itenm(m),
     &              idas_itreg(itreg(m)),
     &              das_iterg(iterg(m)),tr,
     &              itrnv(m),idas(itriv(m)),das(itrrv(m)),
     &              itxnm(m),itynm(m),itznm(m),
     &              das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &              das_itzrg(itzrg(m)),
     &              igsh,idasa)

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(1,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(1,trTEMP_sum,ncol,m)
                  endif

               end if

*-----------------------------------------------------------------------

            else if( itals(m) .eq. 36 ) then

                     idas0 = nmmax
                     idas1 = lmmax
                     idasa = idas1 + itenm(m)

                  if( npe .gt. 1 ) then
                     call paratal(tr,trTEMP,mnmax)
                     call paratal_sumover(trTEMP_sum,m)
                  endif

               if( me .eq. 0 ) then

                  !OBINATA(2012.9.7)
                  if( npe .gt. 1 ) then
                    call sumtr(tr,trRES(irestalm(m)),mnmax,ncol)
                    call sumtr_sumover(m,ncol)   !  tr_sum -> trRES_sum
                   endif

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(0,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(0,trTEMP_sum,ncol,m)
                  endif

                  call psedrz(m,
     &              itpan(m),itrnm(m),itznm(m),itenm(m),
     &              das_itrrg(itrrg(m)),das_itzrg(itzrg(m)),
     &              das_iterg(iterg(m)),tr,
     &              idasa)

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(1,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(1,trTEMP_sum,ncol,m)
                  endif

               end if

*-----------------------------------------------------------------------

            else if( itals(m) .eq. 37 ) then

                     idas0 = nmmax
                     idas1 = lmmax
                     idas2 = ( idas1 + itenm(m) - 1 ) * 2 + 1
                     idasa = idas1 + itenm(m)

                  if( npe .gt. 1 ) then
                     call paratal(tr,trTEMP,mnmax)
                     call paratal_sumover(trTEMP_sum,m)
                  endif

               if( me .eq. 0 ) then

                  !OBINATA(2012.9.7)
                  if( npe .gt. 1 ) then
                    call sumtr(tr,trRES(irestalm(m)),mnmax,ncol)
                    call sumtr_sumover(m,ncol)   !  tr_sum -> trRES_sum
                   endif

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(0,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(0,trTEMP_sum,ncol,m)
                  endif

                  call psedxyz(m,
     &              itpan(m),
     &              itxnm(m),itynm(m),itznm(m),itenm(m),
     &              das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &              das_itzrg(itzrg(m)),
     &              das_iterg(iterg(m)),tr,
     &              igsh,idasa)

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(1,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(1,trTEMP_sum,ncol,m)
                  endif

               end if

*-----------------------------------------------------------------------

            else if( itals(m) .eq. 41 ) then

               if( ist_cut .eq. 0 .or. itstd(m) .eq. 1 ) then

                     idas0 = nmmax
                     idas1 = lmmax
                     idas2 = idas1 + itenm(m)
                     idas3 = idas2 + ittnm(m)
                     idasa = idas3

                  if( npe .gt. 1 ) then
                     call paratal(tr,trTEMP,mnmax)
                     call paratal_sumover(trTEMP_sum,m)
                  endif

               if( me .eq. 0 ) then

                  !OBINATA(2012.6.6)
                  if( npe .gt. 1 ) then
                    call sumtr(tr,trRES(irestalm(m)),mnmax,ncol)
                    call sumtr_sumover(m,ncol)   !  tr_sum -> trRES_sum
                   endif

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(0,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(0,trTEMP_sum,ncol,m)
                  endif


                  call ptpointp(m,
     &              itpan(m),itmsh(m),itenm(m),itmst(m),ittnm(m),
     &              mftal(m),
     &              das_iterg(iterg(m)),das_ittrg(ittrg(m)),tr,idasa)

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(1,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(1,trTEMP_sum,ncol,m)
                  endif

               end if

               end if

*-----------------------------------------------------------------------

            else if( itals(m) .eq. 42 ) then

                     idas0 = nmmax
                     idas1 = lmmax
                     idas2 = idas1 + itenm(m)
                     idas3 = idas2 + ittnm(m)
                     idas4 = ( idas3 + itrgn(m) - 1 ) * 2 + 1
                     idas5 = idas3 + itrgn(m)
     &                     + ( itrgn(m) + mod(itrgn(m),2) ) / 2
                     idas6 = ( idas5 + itrgn(m) - 1 ) * 2 + 1
                     idasa = idas5

                  if( npe .gt. 1 ) then
                     call paratal(tr,trTEMP,mnmax)
                     call paratal_sumover(trTEMP_sum,m)
                  endif

               if( me .eq. 0 ) then

                  !OBINATA(2012.6.6)
                  if( npe .gt. 1 ) then
                    call sumtr(tr,trRES(irestalm(m)),mnmax,ncol)
                    call sumtr_sumover(m,ncol)   !  tr_sum -> trRES_sum
                   endif

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(0,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(0,trTEMP_sum,ncol,m)
                  endif

                  call pwwgreg(m,
     &              itpan(m),itrgn(m),itrgm(m),itenm(m),itmst(m),
     &              ittnm(m),
     &              idas_itreg(itreg(m)),
     &              das_iterg(iterg(m)),das_ittrg(ittrg(m)),tr,
     &              itrnv(m),idas(itriv(m)),das(itrrv(m)),
     &              itxnm(m),itynm(m),itznm(m),mftal(m), ! S.H. 2024.3.21
     &              das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &              das_itzrg(itzrg(m)),
     &              igsh,idasa)

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(1,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(1,trTEMP_sum,ncol,m)
                  endif

               end if

*-----------------------------------------------------------------------

            else if( itals(m) .eq. 48 ) then

                     idas0 = nmmax
                     idas1 = lmmax
                     idas2 = idas1 + itenm(m)
                     idas3 = ( idas2 + ittnm(m) - 1 ) * 2 + 1
                     idasa = idas2 + ittnm(m)

                  if( npe .gt. 1 ) then
                     call paratal(tr,trTEMP,mnmax)
                     call paratal_sumover(trTEMP_sum,m)
                  endif

               if( me .eq. 0 ) then

                  !OBINATA(2012.6.26)
                  if( npe .gt. 1 ) then
                    call sumtr(tr,trRES(irestalm(m)),mnmax,ncol)
                    call sumtr_sumover(m,ncol)   !  tr_sum -> trRES_sum
                   endif

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(0,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(0,trTEMP_sum,ncol,m)
                  endif

                  call pwwgxyz(m,
     &              itpan(m),itmtn(m),ismte(itmtt(m)),
     &              itxnm(m),itynm(m),itznm(m),itenm(m),itmst(m),
     &              ittnm(m),mftal(m), ! S.H. 2024.3.21
     &              das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &              das_itzrg(itzrg(m)),
     &              das_iterg(iterg(m)),das_ittrg(ittrg(m)),tr,
     &              igsh,idasa)

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(1,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(1,trTEMP_sum,ncol,m)
                  endif

               end if

cKN 2015/03/17
*-----------------------------------------------------------------------
cKN 2016/08/13

            else if( itals(m) .eq. 55 ) then

                     idas0 = nmmax
                     idas1 = lmmax
                     idas2 = idas1 + itenm(m)
                     idas3 = idas2 + ittnm(m)
                     idas4 = ( idas3 + itrgn(m) - 1 ) * 2 + 1
                     idas5 = idas3 + itrgn(m)
     &                     + ( itrgn(m) + mod(itrgn(m),2) ) / 2
                     idas6 = ( idas5 + itrgn(m) - 1 ) * 2 + 1
                     idasa = idas5

                  if( npe .gt. 1 ) then
                     call paratal(tr,trTEMP,mnmax)
                     call paratal_sumover(trTEMP_sum,m)
                  endif

               if( me .eq. 0 ) then

                  !OBINATA(2012.6.6)
                  if( npe .gt. 1 ) then
                    call sumtr(tr,trRES(irestalm(m)),mnmax,ncol)
                    call sumtr_sumover(m,ncol)   !  tr_sum -> trRES_sum
                   endif

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(0,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(0,trTEMP_sum,ncol,m)
                  endif
                  
!<-20220509murofushi update
!--                  call pwwgreg(m,
!--     &              itpan(m),itrgn(m),itrgm(m),itenm(m),itmst(m),
!--     &              ittnm(m),
!--     &              idas(itreg(m)),
!--     &              das(iterg(m)),das(ittrg(m)),tr,
!--     &              itrnv(m),idas(itriv(m)),das(itrrv(m)),
!--     &              itxnm(m),itynm(m),itznm(m),
!--     &              das(itxrg(m)),das(ityrg(m)),das(itzrg(m)),
!--     &              igsh,idasa)
                  call pwwgtet(m,
     &              itpan(m),itrgn(m),itrgm(m),itenm(m),
     &              itmst(m),ittnm(m),
     &              das_iterg(iterg(m)),das_ittrg(ittrg(m)),tr,
     &              itxnm(m),itynm(m),itznm(m),mftal(m), ! S.H. 2024.3.21
     &              idas_itreg(itreg(m)),
     &              das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &              das_itzrg(itzrg(m)),
     &              igsh,idasa)
!--->

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(1,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(1,trTEMP_sum,ncol,m)
                  endif

               end if

*-----------------------------------------------------------------------
c iwamoto 2012/12/4
c        : production of dchain files at the end of job

*-----------------------------------------------------------------------

            else if( itals(m) .eq. 46 .and.
     &             ( icntl .eq. 13 .or. icntl .eq. 14 ) ) then

                     idas1 = lmmax
                     idas2 = ( idas1 + itrgn(m) - 1 ) * 2 + 1
                     idas3 = idas2 + itrgn(m)
     &                     + ( itrgn(m) + mod(itrgn(m),2) ) / 2

                  if( npe .gt. 1 ) then
                     call paratal(tr,trTEMP,mnmax)
                     call paratal_sumover(trTEMP_sum,m)
                  endif

               if( me .eq. 0 ) then

                  !OBINATA(2012.6.6)
                  if( npe .gt. 1 ) then
                    call sumtr(tr,trRES(irestalm(m)),mnmax,ncol)
                    call sumtr_sumover(m,ncol)   !  tr_sum -> trRES_sum
                   endif

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(0,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(0,trTEMP_sum,ncol,m)
                  endif

                  call pvlmreg(m,itrcm(m),idas_itrcg(itrcg(m)),
     &                         itrgn(m),itrgm(m),
     &                         idas_itreg(itreg(m)),tr,
     &                         das(idas1),idas(idas2),
     &                         itrnv(m),idas(itriv(m)),das(itrrv(m)),
     &                         idas3)

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(1,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(1,trTEMP_sum,ncol,m)
                  endif

               end if

*-----------------------------------------------------------------------

            else if( itals(m) .eq. 47 .and.
     &             ( icntl .eq. 13 .or. icntl .eq. 15 ) ) then

                     idas1 = lmmax
                     idas2 = ( idas1 + itrgn(m) - 1 ) * 2 + 1
                     idas3 = idas2 + itrgn(m)
     &                     + ( itrgn(m) + mod(itrgn(m),2) ) / 2
                     idasa = idas3

                  if( npe .gt. 1 ) then
                     call paratal(tr,trTEMP,mnmax)
                     call paratal_sumover(trTEMP_sum,m)
                  endif

               if( me .eq. 0 ) then

                  !OBINATA(2012.6.6)
                  if( npe .gt. 1 ) then
                    call sumtr(tr,trRES(irestalm(m)),mnmax,ncol)
                    call sumtr_sumover(m,ncol)   !  tr_sum -> trRES_sum
                   endif

cABE 2020/08/05 change "itall .ne. 0" to "itall .ne. -1"
                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(0,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(0,trTEMP_sum,ncol,m)
                  endif

                  call pwwbgreg(m,itrcm(m),idas_itrcg(itrcg(m)),
     &                      itrgn(m),itrgm(m),
     &                      idas_itreg(itreg(m)),tr,
     &                      das(idas1),idas(idas2),
     &                      itrnv(m),idas(itriv(m)),das(itrrv(m)),
     &                      idas3,
     &                      itxnm(m),itynm(m),itznm(m),
     &                      das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &                      das_itzrg(itzrg(m)),
     &                      igsh,idasa)

cABE 2020/08/05 change "itall .ne. 0" to "itall .ne. -1"
                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(1,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(1,trTEMP_sum,ncol,m)
                  endif

               end if

*-----------------------------------------------------------------------
cFURUTA20240118               
            else if( itals(m) .eq. 56 .and.
     &             ( icntl .eq. 13 .or. icntl .eq. 15 ) ) then

                     idas0 = nmmax
                     idas1 = lmmax
                     idas2 = idas1 + itenm(m)
                     idas3 = ( idas2 + ittnm(m) - 1 ) * 2 + 1
                     idasa = idas2 + ittnm(m)

                  if( npe .gt. 1 ) then
                     call paratal(tr,trTEMP,mnmax)
                     call paratal_sumover(trTEMP_sum,m)
                  endif

               if( me .eq. 0 ) then

                  !OBINATA(2012.6.26)
                  if( npe .gt. 1 ) then
                    call sumtr(tr,trRES(irestalm(m)),mnmax,ncol)
                    call sumtr_sumover(m,ncol)   !  tr_sum -> trRES_sum
                   endif

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(0,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(0,trTEMP_sum,ncol,m)
                  endif

                  call pwwbgxyz(m,tr,
     &                      itxnm(m),itynm(m),itznm(m),
     &                      das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &                      das_itzrg(itzrg(m)),
     &                      igsh,idasa)
!--->

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(1,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(1,trTEMP_sum,ncol,m)
                  endif

               end if

*-----------------------------------------------------------------------
cFURUTA20240118
            else if( itals(m) .eq. 57 .and.
     &             ( icntl .eq. 13 .or. icntl .eq. 15 ) ) then

                     idas0 = nmmax
                     idas1 = lmmax
                     idas2 = idas1 + itenm(m)
                     idas3 = idas2 + ittnm(m)
                     idas4 = ( idas3 + itrgn(m) - 1 ) * 2 + 1
                     idas5 = idas3 + itrgn(m)
     &                     + ( itrgn(m) + mod(itrgn(m),2) ) / 2
                     idas6 = ( idas5 + itrgn(m) - 1 ) * 2 + 1
                     idasa = idas5

                  if( npe .gt. 1 ) then
                     call paratal(tr,trTEMP,mnmax)
                     call paratal_sumover(trTEMP_sum,m)
                  endif

               if( me .eq. 0 ) then

                  !OBINATA(2012.6.6)
                  if( npe .gt. 1 ) then
                    call sumtr(tr,trRES(irestalm(m)),mnmax,ncol)
                    call sumtr_sumover(m,ncol)   !  tr_sum -> trRES_sum
                   endif

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(0,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(0,trTEMP_sum,ncol,m)
                  endif
                  
                  call pwwbgtet(m,itrgn(m),itrgm(m),tr,
     &              itxnm(m),itynm(m),itznm(m),
     &              idas_itreg(itreg(m)),
     &              das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &              das_itzrg(itzrg(m)),
     &              igsh,idasa)
!--->

                  if( itall .ne.-1 .and. npe .le. 1 .and. ncol .ne. 2 )
     &                                                             then
                     call tallsave(1,tr,trTEMP,mnmax,ncol)
                     call tallsave_sumover(1,trTEMP_sum,ncol,m)
                  endif

               end if

*-----------------------------------------------------------------------
            end if

*-----------------------------------------------------------------------7

            if( me .eq. 0 ) then

               if( rtstd(m) .gt. 0.0 ) then

                  instop = instop + 1
                  idstop = idstop + itstd(m)

               end if

            end if

*-----------------------------------------------------------------------

         end do

*-----------------------------------------------------------------------
                  if ( itall.eq.4 ) then
                    tr0 = 0.d0
                    if ( me.eq.0 .and. irestart.ne.0 ) then
                      itl4flg = 2
                      call read_talls(28,ierr)
                    end if
                  end if

*-----------------------------------------------------------------------
         deallocate(trTEMP)
         deallocate(trTEMP_sum)

      end if

*-----------------------------------------------------------------------
C MATSUDA 2024.12.09 (pmultip)

      if( ncol .eq. 2 ) then

         do m = 1, imltp

            if( imltf(m) .gt. 0 ) then
               call pmultip(m,1)
            end if

         end do

      end if

*-----------------------------------------------------------------------

  999 denstepold = 0.0d0 ! Reset EGS5 energy hinge related parameter used for tally, T.Sato 2017/06/27
      denstepnew = 0.0d0
      deinit = 0.0d0

      return
      end


************************************************************************
*                                                                      *
      subroutine gshows
*                                                                      *
*       output gshow tally only                                        *
*       last modified by K.Niita on 2002/01/04                         *
*                                                                      *
************************************************************************
      use GGBANKMOD !FURUTA
      use moddas_mesh
*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

      include 'param.inc'

*-----------------------------------------------------------------------

      common /talmm/  nmmax, lmmax, itlmx
      common /tall00/ itnm, ital(itlmax), itals(itlmax), italm(itlmax)

      common /tall04/ itxty(itlmax), itxnm(itlmax), itxrg(itlmax),
     &                rtxmi(itlmax), rtxma(itlmax), rtxdl(itlmax)
      common /tall05/ ityty(itlmax), itynm(itlmax), ityrg(itlmax),
     &                rtymi(itlmax), rtyma(itlmax), rtydl(itlmax)
      common /tall06/ itzty(itlmax), itznm(itlmax), itzrg(itlmax),
     &                rtzmi(itlmax), rtzma(itlmax), rtzdl(itlmax)

      dimension     idas(1)
      equivalence ( das, idas )

      common /ccggg/  icgg

*-----------------------------------------------------------------------

         do m = 1, itnm

            if( itals(m) .eq. 25 ) then

              if(icgg.ne.0)then
                  call ALLOCATE_GGBANK !FURUTA
                  call INIT_GGBANK     !FURUTA
              endif

                     idas1 = ( lmmax - 1 ) * 2 + 1
                     idasa = lmmax

                  call pgshxyz(m,
     &              itxnm(m),itynm(m),itznm(m),
     &              das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &              das_itzrg(itzrg(m)),
     &              idasa)

                  if(icgg.ne.0) call DEALLOCATE_GGBANK !FURUTA

            end if

         end do

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine gshowt
*                                                                      *
*       output gshow in xyz tally only                                 *
*       last modified by K.Niita on 2002/07/18                         *
*                                                                      *
************************************************************************
      use GGBANKMOD !FURUTA
      use TALMOD
!$       use TALMOD0
      use partmod, only: itpan, itpat, jtpat ! frtati 2021/10/05
      use moddas_mesh
      use moddas_tally
*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

      include 'param.inc'
      include 'param01.inc'

*-----------------------------------------------------------------------

      common /talout/ itall
      common /talmm/  nmmax, lmmax, itlmx

*-----------------------------------------------------------------------

      common /tall00/ itnm, ital(itlmax), itals(itlmax), italm(itlmax)
      common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                itout(itlmax), ittwo(itlmax)
      common /tall02/ itrgn(itlmax), itrgm(itlmax), itreg(itlmax)
      common /tall25/ itrnv(itlmax), itriv(itlmax), itrrv(itlmax)
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
      common /tall08/ rtrx0(itlmax), rtry0(itlmax)
      common /tall10/ itaxn(itlmax), itaxs(itlmax,6)
      common /tall11/ itfln(itlmax), itfll(itlmax,6), ctfln(itlmax,6)
      character ctfln*100

      common /tall12/ itman(itlmax), itmat(itlmax), itmct(itlmax),
     &                itnun(itlmax), itnuc(itlmax),
     &                itndz(itlmax), itndn(itlmax),
     &                itnkz(itlmax), itnkn(itlmax)

      common /tall13/ itrcn(itlmax), itrcr(itlmax), itrca(itlmax)
      common /tall26/ itrcs(itlmax), itrcc(itlmax), itrss(itlmax)

      common /tall14/ ittll(itlmax), ittle(itlmax)
      character ittle*80

      common /tall15/ itanl(itlmax), itang(itlmax)
      character itang*200

      common /tall27/ itaxl(itlmax), itaxt(itlmax),
     &                itayl(itlmax), itayt(itlmax),
     &                itazl(itlmax), itazt(itlmax)
      character itaxt*200, itayt*200, itazt*200

      common /tall16/ itnfr(itlmax), itnfz(itlmax), itnfn(itlmax)
      common /tall17/ itndy(itlmax)
      common /tall18/ ithet(itlmax)
      common /tall19/ itsmn(itlmax), itstm(itlmax)
      common /tall20/ ittty(itlmax), ittnm(itlmax), ittrg(itlmax),
     &                rttmi(itlmax), rttma(itlmax), rttdl(itlmax)
      common /tall21/ rtfac(itlmax)
      common /tall22/ itmcn(itlmax), itmtn(itlmax), itmtt(itlmax)
      common /tall23/ itaty(itlmax), itanm(itlmax), itarg(itlmax),
     &                rtami(itlmax), rtama(itlmax), rtadl(itlmax)
      common /tall28/ itgsh(itlmax)
      common /tall38/ itres(itlmax)
      common /tall39/ rtwid(itlmax)
      common /tall45/ itmlp(itlmax), itmln(itlmax,6), itmst(itlmax),
     &                itmli(itlmax,6), rtmme(itlmax,6), itmnt(itlmax,6),
     &                itmpn(itlmax,6), itmpt(itlmax,6,6,2)
      common /tall73/ itenclo(itlmax), itangform(itlmax)

      common /tall83/ itnzn(itlmax), itndm(itlmax)

      common /tall91/ itextstat(itlmax), mftal(itlmax) ! S.H. extstat 2024.3.18

C --- add NS 2020.04 ---
C for nonshared_tally option
      integer italsh
      common /talsh/ italsh

      dimension     idas(1)
      equivalence ( das, idas )

      common /ccggg/  icgg

      real(8),pointer:: tr(:)

*-----------------------------------------------------------------------

            igsh = 1

         do m = 1, itnm

         if( itgsh(m) .ne. 0 .and.
     &     ( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) ) then

C for nonshared_tally option
!$          if(italsh .eq. 0) then
!$            call GET_TR_HEAD_POINTER0(tr,m)
!$          else
              call GET_TR_HEAD_POINTER(tr,m)
!$          end if
              if(icgg.ne.0)then
                  call ALLOCATE_GGBANK !FURUTA
                  call INIT_GGBANK     !FURUTA
              endif

               if( itals(m) .eq. 3 ) then

                     idas1 = lmmax
                     idas2 = idas1 + itenm(m)
                     idas3 = ( idas2 + ittnm(m) - 1 ) * 2 + 1
                     idasa = idas2 + ittnm(m)

                  call ptracxyz(m,
     &              itpan(m),itmtn(m),ismte(itmtt(m)),
     &              itxnm(m),itynm(m),itznm(m),itenm(m),itmst(m),
     &              ittnm(m),mftal(m),
     &              das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &              das_itzrg(itzrg(m)),
     &              das_iterg(iterg(m)),das_ittrg(ittrg(m)),tr,
     &              igsh,idasa)

               else if( itals(m) .eq. 12 ) then

                     idas1 = lmmax
                     idas2 = idas1 + itenm(m)
                     idas3 = idas2 + itanm(m)
                     idas4 = ( idas3 + ittnm(m) - 1 ) * 2 + 1
                     idasa = idas3 + ittnm(m)

                  if( itenclo(m) .eq. 1 ) then

                     call psufxyz_rpp(m,
     &                 itpan(m),itxnm(m),itynm(m),itznm(m),
     &                 itenm(m),itanm(m),ittnm(m),itmst(m),
     &                 das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &                 das_itzrg(itzrg(m)),
     &                 das_iterg(iterg(m)),das_itarg(itarg(m)),
     &                 das_ittrg(ittrg(m)),
     &                 tr,igsh,idasa)

                  else
cKN 2015/10/30 missed
                     call psufxyz(m,
     &                 itpan(m),itxnm(m),itynm(m),itznm(m),
     &                 itenm(m),itanm(m),ittnm(m),itmst(m),
     &                 das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &                 das_itzrg(itzrg(m)),
     &                 das_iterg(iterg(m)),das_itarg(itarg(m)),
     &                 das_ittrg(ittrg(m)),
     &                 tr,igsh,idasa)

                  endif

               else if( itals(m) .eq. 8 ) then

                     idas1 = lmmax
                     idas2 = idas1 + ( maxpt + maxnt ) * 2
                     idas3 = ( idas2 + itnfn(m) * itnfn(m) - 1 ) * 2 + 1
                     idasa = idas2 + itnfn(m) * itnfn(m)

C S.H. revised to fix a bug of t-yield with gshow option (2022.9.14)
                  call pyildxyz(m,itndz(m),itndn(m),itndm(m),itnfn(m),
     &              itmtn(m),ismte(itmtt(m)),
     &              itxnm(m),itynm(m),itznm(m),itnun(m),
     &              das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &              das_itzrg(itzrg(m)),
     &              isnuc(itnuc(m)),idas(itnkz(m)),idas(itnkn(m)),
     &              tr,igsh,idasa)

               else if( itals(m) .eq. 11 ) then

                     idas1 = ( lmmax - 1 ) * 2 + 1
                     idasa = lmmax

                  call phetxyz(m,itndy(m),itout(m),
     &              itpan(m),itmtn(m),ismte(itmtt(m)),
     &              itxnm(m),itynm(m),itznm(m),
     &              das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &              das_itzrg(itzrg(m)),
     &              itenm(m),das_iterg(iterg(m)),
     &              tr,das(ithet(m)),
     &              igsh,idasa)

               else if( itals(m) .eq. 15 ) then

                     idas1 = lmmax
                     idas2 = idas1 + itenm(m)
                     idas3 = ( idas1 + ittnm(m) - 1 ) * 2 + 1
                     idasa = idas2 + ittnm(m)

                  call pstarxyz(m,
     &              itpan(m),itmtn(m),ismte(itmtt(m)),
     &              itxnm(m),itynm(m),itznm(m),itenm(m),
     &              ittnm(m),
     &              das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &              das_itzrg(itzrg(m)),
     &              das_iterg(iterg(m)),das_ittrg(ittrg(m)),tr,
     &              igsh,idasa)

               else if( itals(m) .eq. 18 ) then

                     idas1 = lmmax
                     idas2 = idas1 + itenm(m)
                     idas3 = ( idas2 + ittnm(m) - 1 ) * 2 + 1
                     idasa = idas2 + ittnm(m)

                  call ptimexyz(m,
     &              itmtn(m),ismte(itmtt(m)),
     &              itpan(m),itxnm(m),itynm(m),itznm(m),itenm(m),
     &              ittnm(m),
     &              das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &              das_itzrg(itzrg(m)),
     &              das_iterg(iterg(m)),das_ittrg(ittrg(m)),tr,
     &              igsh,idasa)

               else if( itals(m) .eq. 21 ) then

                     idas1 = ( lmmax - 1 ) * 2 + 1
                     idasa = lmmax

                  call pdpaxyz(m,itndy(m),itout(m),
     &              itpan(m),itmtn(m),ismte(itmtt(m)),
     &              itxnm(m),itynm(m),itznm(m),
     &              das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &              das_itzrg(itzrg(m)),
     &              tr,igsh,idasa)

               else if( itals(m) .eq. 24 ) then

                     idas1 = lmmax
                     idas2 = idas1 + itenm(m)
                     idas3 = idas2 + ittnm(m)
                     idas4 = ( idas3 + itanm(m) - 1 ) * 2 + 1
                     idasa = idas3 + itanm(m)

                  call ppdctxyz(m,
     &              itpan(m),itmtn(m),ismte(itmtt(m)),
     &              itxnm(m),itynm(m),itznm(m),itenm(m),
     &              ittnm(m),itanm(m),
     &              das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &              das_itzrg(itzrg(m)),
     &              das_iterg(iterg(m)),das_ittrg(ittrg(m)),
     &              das_itarg(itarg(m)),
     &              tr,igsh,idasa)

               else if( itals(m) .eq. 30 ) then

                     idas1 = lmmax
                     idas2 = ( idas1 + itenm(m) - 1 ) * 2 + 1
                     idasa = idas1 + itenm(m)

                  call pletxyz(m,
     &              itpan(m),itmtn(m),ismte(itmtt(m)),
     &              itxnm(m),itynm(m),itznm(m),itenm(m),
     &              das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &              das_itzrg(itzrg(m)),
     &              das_iterg(iterg(m)),tr,
     &              igsh,idasa)

               else if( itals(m) .eq. 33 ) then

                     idas1 = lmmax
                     idas2 = ( idas1 + ittnm(m) - 1 ) * 2 + 1
                     idasa = idas1 + ittnm(m)

                  call pdepstxyz(m,
     &              itpan(m),itmtn(m),ismte(itmtt(m)),
     &              itxnm(m),itynm(m),itznm(m),itenm(m),
     &              ittnm(m),
     &              das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &              das_itzrg(itzrg(m)),
     &              das_iterg(iterg(m)),das_ittrg(ittrg(m)),tr,
     &              igsh,idasa)

               else if( itals(m) .eq. 37 ) then

                     idas1 = lmmax
                     idas2 = ( idas1 + itenm(m) - 1 ) * 2 + 1
                     idasa = idas1 + itenm(m)

                  call psedxyz(m,
     &              itpan(m),
     &              itxnm(m),itynm(m),itznm(m),itenm(m),
     &              das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &              das_itzrg(itzrg(m)),
     &              das_iterg(iterg(m)),tr,
     &              igsh,idasa)

               end if

                  if(icgg.ne.0) call DEALLOCATE_GGBANK !FURUTA

            end if

         end do

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine rshows
*                                                                      *
*       output rshow tally only                                        *
*       last modified by K.Niita on 2002/01/09                         *
*                                                                      *
************************************************************************
      use GGBANKMOD !FURUTA
      use moddas_mesh
      use moddas_region
*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

      include 'param.inc'

*-----------------------------------------------------------------------

      common /talmm/  nmmax, lmmax, itlmx
      common /tall00/ itnm, ital(itlmax), itals(itlmax), italm(itlmax)

      common /tall02/ itrgn(itlmax), itrgm(itlmax), itreg(itlmax)
      common /tall25/ itrnv(itlmax), itriv(itlmax), itrrv(itlmax)

      common /tall04/ itxty(itlmax), itxnm(itlmax), itxrg(itlmax),
     &                rtxmi(itlmax), rtxma(itlmax), rtxdl(itlmax)
      common /tall05/ ityty(itlmax), itynm(itlmax), ityrg(itlmax),
     &                rtymi(itlmax), rtyma(itlmax), rtydl(itlmax)
      common /tall06/ itzty(itlmax), itznm(itlmax), itzrg(itlmax),
     &                rtzmi(itlmax), rtzma(itlmax), rtzdl(itlmax)

      dimension     idas(1)
      equivalence ( das, idas )

      common /ccggg/  icgg

*-----------------------------------------------------------------------

         do m = 1, itnm

            if( itals(m) .eq. 26 ) then

              if(icgg.ne.0)then
                  call ALLOCATE_GGBANK !FURUTA
                  call INIT_GGBANK     !FURUTA
              endif

                     idas1 = lmmax
                     idas2 = ( idas1 + itrgn(m) - 1 ) * 2 + 1
                     idas3 = idas2 + itrgn(m)
                     idasa = idas1 + itrgn(m)
     &                     + ( itrgn(m) + mod(itrgn(m),2) ) / 2

                  call prshxyz(m,
     &              itxnm(m),itynm(m),itznm(m),
     &              das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &              das_itzrg(itzrg(m)),
     &              itrgn(m),itrgm(m),idas_itreg(itreg(m)),
     &              itrnv(m),idas(itriv(m)),das(itrrv(m)),
     &              idasa)

                  if(icgg.ne.0) call DEALLOCATE_GGBANK !FURUTA

            end if

         end do

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine rshowt
*                                                                      *
*       output rshow in reg tally only                                 *
*       last modified by K.Niita on 2002/01/18                         *
*                                                                      *
************************************************************************
      use GGBANKMOD !FURUTA
      use TALMOD
!$       use TALMOD0
      use partmod, only: itpan, itpat, jtpat ! frtati 2021/10/05
      use moddas_mesh
      use moddas_region
      use moddas_tally
*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

      include 'param.inc'
      include 'param01.inc'

*-----------------------------------------------------------------------

      common /talout/ itall
      common /talmm/  nmmax, lmmax, itlmx

*-----------------------------------------------------------------------

      common /tall00/ itnm, ital(itlmax), itals(itlmax), italm(itlmax)
      common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                itout(itlmax), ittwo(itlmax)
      common /tall02/ itrgn(itlmax), itrgm(itlmax), itreg(itlmax)
      common /tall25/ itrnv(itlmax), itriv(itlmax), itrrv(itlmax)
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
      common /tall08/ rtrx0(itlmax), rtry0(itlmax)
      common /tall10/ itaxn(itlmax), itaxs(itlmax,6)
      common /tall11/ itfln(itlmax), itfll(itlmax,6), ctfln(itlmax,6)
      character ctfln*100

      common /tall12/ itman(itlmax), itmat(itlmax), itmct(itlmax),
     &                itnun(itlmax), itnuc(itlmax),
     &                itndz(itlmax), itndn(itlmax),
     &                itnkz(itlmax), itnkn(itlmax)

      common /tall13/ itrcn(itlmax), itrcr(itlmax), itrca(itlmax)
      common /tall26/ itrcs(itlmax), itrcc(itlmax), itrss(itlmax)

      common /tall14/ ittll(itlmax), ittle(itlmax)
      character ittle*80

      common /tall15/ itanl(itlmax), itang(itlmax)
      character itang*200

      common /tall27/ itaxl(itlmax), itaxt(itlmax),
     &                itayl(itlmax), itayt(itlmax),
     &                itazl(itlmax), itazt(itlmax)
      character itaxt*200, itayt*200, itazt*200

      common /tall16/ itnfr(itlmax), itnfz(itlmax), itnfn(itlmax)
      common /tall17/ itndy(itlmax)
      common /tall18/ ithet(itlmax)
      common /tall19/ itsmn(itlmax), itstm(itlmax)
      common /tall20/ ittty(itlmax), ittnm(itlmax), ittrg(itlmax),
     &                rttmi(itlmax), rttma(itlmax), rttdl(itlmax)
      common /tall21/ rtfac(itlmax)
      common /tall22/ itmcn(itlmax), itmtn(itlmax), itmtt(itlmax)
      common /tall23/ itaty(itlmax), itanm(itlmax), itarg(itlmax),
     &                rtami(itlmax), rtama(itlmax), rtadl(itlmax)
      common /tall28/ itgsh(itlmax)
      common /tall29/ itrsh(itlmax)
      common /tall38/ itres(itlmax)
      common /tall39/ rtwid(itlmax)
      common /tall45/ itmlp(itlmax), itmln(itlmax,6), itmst(itlmax),
     &                itmli(itlmax,6), rtmme(itlmax,6), itmnt(itlmax,6),
     &                itmpn(itlmax,6), itmpt(itlmax,6,6,2)
      common /tall72/ itactnm(itlmax), itactrg(itlmax), itactmax(itlmax)      ! S.Abe 2018/02/15
      common /tall83/ itnzn(itlmax), itndm(itlmax) ! frtati 2022/09/15

      common /tall91/ itextstat(itlmax), mftal(itlmax) ! S.H. extstat 2024.3.18

      dimension     idas(1)
      equivalence ( das, idas )

      common /ccggg/  icgg

      real(8),pointer:: tr(:)

*-----------------------------------------------------------------------

            igsh = 1

         do m = 1, itnm

            if( itrsh(m) .ne. 0 ) then

C for nonshared_tally option
!$          if(italsh .eq. 0) then
!$            call GET_TR_HEAD_POINTER0(tr,m)
!$          else
              call GET_TR_HEAD_POINTER(tr,m)
!$          end if
              if(icgg.ne.0)then
                  call ALLOCATE_GGBANK !FURUTA
                  call INIT_GGBANK     !FURUTA
              endif

               if( itals(m) .eq. 1 ) then

                     idas1 = lmmax
                     idas2 = idas1 + itenm(m)
                     idas3 = idas2 + ittnm(m)
                     idas4 = ( idas3 + itrgn(m) - 1 ) * 2 + 1
                     idas5 = idas3 + itrgn(m)
     &                     + ( itrgn(m) + mod(itrgn(m),2) ) / 2
                     idas6 = ( idas5 + itrgn(m) - 1 ) * 2 + 1
                     idasa = idas5

                  call ptracreg(m,
     &              itpan(m),itrgn(m),itrgm(m),itenm(m),itmst(m),
     &              ittnm(m),mftal(m),
     &              idas_itreg(itreg(m)),
     &              das_iterg(iterg(m)),das_ittrg(ittrg(m)),tr,
     &              itrnv(m),idas(itriv(m)),das(itrrv(m)),
     &              itxnm(m),itynm(m),itznm(m),
     &              das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &              das_itzrg(itzrg(m)),
     &              igsh,idasa)

               else if( itals(m) .eq. 6 ) then

                     idas1 = lmmax
                     idas2 = idas1 + ( maxnt + maxpt ) * 2
                     idas3 = ( idas2 + itrgn(m) - 1 ) * 2 + 1
                     idas4 = idas2 + itrgn(m)
     &                     + ( itrgn(m) + mod(itrgn(m),2) ) / 2
                     idas5 = ( idas4 + itrgn(m) - 1 ) * 2 + 1
                     idasa = idas4

                  call pyildreg(m,itndz(m),itndn(m),itndm(m),
     &              itrgn(m),itrgm(m),itnun(m),
     &              idas_itreg(itreg(m)),
     &              isnuc(itnuc(m)),idas(itnkz(m)),idas(itnkn(m)),
     &              tr,itrnv(m),idas(itriv(m)),das(itrrv(m)),
     &              itxnm(m),itynm(m),itznm(m),
     &              das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &              das_itzrg(itzrg(m)),
     &              igsh,idasa)

               else if( itals(m) .eq. 9 ) then

                     idas1 = lmmax
                     idas2 = ( idas1 + itrgn(m) - 1 ) * 2 + 1
                     idas3 = idas1 + itrgn(m)
     &                     + ( itrgn(m) + mod(itrgn(m),2) ) / 2
                     idas4 = ( idas3 + itrgn(m) - 1 ) * 2 + 1
                     idasa = idas3

                  call phetreg(m,itndy(m),itout(m),
     &              itpan(m),itrgn(m),itrgm(m),
     &              idas_itreg(itreg(m)),itenm(m),das_iterg(iterg(m)),
     &              tr,das(ithet(m)),
     &              itrnv(m),idas(itriv(m)),das(itrrv(m)),
     &              itxnm(m),itynm(m),itznm(m),
     &              das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &              das_itzrg(itzrg(m)),
     &              igsh,idasa)

               else if( itals(m) .eq. 13 ) then

                     idas1 = lmmax
                     idas2 = idas1 + itenm(m)
                     idas3 = idas2 + ittnm(m)
                     idas4 = idas3 + itactnm(m)   ! S.Abe 2018/02/15
                     idas5 = ( idas4 + itrgn(m) - 1 ) * 2 + 1
                     idas6 = idas4 + itrgn(m)
     &                     + ( itrgn(m) + mod(itrgn(m),2) ) / 2
                     idas7 = ( idas6 + itrgn(m) - 1 ) * 2 + 1
                     idasa = idas6

                  call pstarreg(m,
     &              itpan(m),itrgn(m),itrgm(m),itenm(m),ittnm(m),
     &              itactnm(m),                                         ! S.Abe 2018/02/15, add itactnm(m)
     &              idas_itreg(itreg(m)),
     &              das_iterg(iterg(m)),
     &              das_ittrg(ittrg(m)),das(itactrg(m)),tr,     ! S.Abe 2018/02/15, add das(itactrg(m))
     &              itrnv(m),idas(itriv(m)),das(itrrv(m)),
     &              itxnm(m),itynm(m),itznm(m),
     &              das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &              das_itzrg(itzrg(m)),
     &              igsh,idasa)

               else if( itals(m) .eq. 16 ) then

                     idas1 = lmmax
                     idas2 = idas1 + itenm(m)
                     idas3 = idas2 + ittnm(m)
                     idas4 = ( idas3 + itrgn(m) - 1 ) * 2 + 1
                     idas5 = idas3 + itrgn(m)
     &                     + ( itrgn(m) + mod(itrgn(m),2) ) / 2
                     idas6 = ( idas5 + itrgn(m) - 1 ) * 2 + 1
                     idasa = idas5

                  call ptimereg(m,
     &              itpan(m),itrgn(m),itrgm(m),itenm(m),ittnm(m),
     &              idas_itreg(itreg(m)),
     &              das_iterg(iterg(m)),das_ittrg(ittrg(m)),tr,
     &              itrnv(m),idas(itriv(m)),das(itrrv(m)),
     &              itxnm(m),itynm(m),itznm(m),
     &              das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &              das_itzrg(itzrg(m)),
     &              igsh,idasa)

               else if( itals(m) .eq. 19 ) then

                     idas1 = lmmax
                     idas2 = ( idas1 + itrgn(m) - 1 ) * 2 + 1
                     idas3 = idas1 + itrgn(m)
     &                     + ( itrgn(m) + mod(itrgn(m),2) ) / 2
                     idas4 = ( idas3 + itrgn(m) - 1 ) * 2 + 1
                     idasa = idas3

                  call pdpareg(m,itndy(m),itout(m),
     &              itpan(m),itrgn(m),itrgm(m),idas_itreg(itreg(m)),
     &              tr,
     &              itrnv(m),idas(itriv(m)),das(itrrv(m)),
     &              itxnm(m),itynm(m),itznm(m),
     &              das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &              das_itzrg(itzrg(m)),
     &              igsh,idasa)

               else if( itals(m) .eq. 22 ) then

                     idas1 = lmmax
                     idas2 = idas1 + itenm(m)
                     idas3 = idas2 + ittnm(m)
                     idas4 = idas3 + itanm(m)
                     idas5 = ( idas4 + itrgn(m) - 1 ) * 2 + 1
                     idas6 = idas4 + itrgn(m)
     &                     + ( itrgn(m) + mod(itrgn(m),2) ) / 2
                     idas7 = ( idas6 + itrgn(m) - 1 ) * 2 + 1
                     idasa = idas6

                  call ppdctreg(m,
     &              itpan(m),itrgn(m),itrgm(m),itenm(m),ittnm(m),
     &              itanm(m),idas_itreg(itreg(m)),
     &              das_iterg(iterg(m)),das_ittrg(ittrg(m)),
     &              das_itarg(itarg(m)),
     &              tr,
     &              itrnv(m),idas(itriv(m)),das(itrrv(m)),
     &              itxnm(m),itynm(m),itznm(m),
     &              das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &              das_itzrg(itzrg(m)),
     &              igsh,idasa)

               else if( itals(m) .eq. 28 ) then

                     idas1 = lmmax
                     idas2 = idas1 + itenm(m)
                     idas3 = ( idas2 + itrgn(m) - 1 ) * 2 + 1
                     idas4 = idas2 + itrgn(m)
     &                     + ( itrgn(m) + mod(itrgn(m),2) ) / 2
                     idas5 = ( idas4 + itrgn(m) - 1 ) * 2 + 1
                     idasa = idas4

                  call pletreg(m,
     &              itpan(m),itrgn(m),itrgm(m),itenm(m),
     &              idas_itreg(itreg(m)),
     &              das_iterg(iterg(m)),tr,
     &              itrnv(m),idas(itriv(m)),das(itrrv(m)),
     &              itxnm(m),itynm(m),itznm(m),
     &              das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &              das_itzrg(itzrg(m)),
     &              igsh,idasa)

               else if( itals(m) .eq. 31 ) then

                     idas1 = lmmax
                     idas2 = idas1 + ittnm(m)
                     idas3 = ( idas2 + itrgn(m) - 1 ) * 2 + 1
                     idas4 = idas2 + itrgn(m)
     &                     + ( itrgn(m) + mod(itrgn(m),2) ) / 2
                     idas5 = ( idas4 + itrgn(m) - 1 ) * 2 + 1
                     idasa = idas3

                  call pdepstreg(m,
     &              itpan(m),itrgn(m),itrgm(m),itenm(m),
     &              ittnm(m),
     &              idas_itreg(itreg(m)),
     &              das_iterg(iterg(m)),das_ittrg(ittrg(m)),tr,
     &              itrnv(m),idas(itriv(m)),das(itrrv(m)),
     &              itxnm(m),itynm(m),itznm(m),
     &              das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &              das_itzrg(itzrg(m)),
     &              igsh,idasa)

               else if( itals(m) .eq. 35 ) then

                     idas1 = lmmax
                     idas2 = idas1 + itenm(m)
                     idas3 = ( idas2 + itrgn(m) - 1 ) * 2 + 1
                     idas4 = idas2 + itrgn(m)
     &                     + ( itrgn(m) + mod(itrgn(m),2) ) / 2
                     idas5 = ( idas4 + itrgn(m) - 1 ) * 2 + 1
                     idasa = idas4

                  call psedreg(m,
     &              itpan(m),itrgn(m),itrgm(m),itenm(m),
     &              idas_itreg(itreg(m)),
     &              das_iterg(iterg(m)),tr,
     &              itrnv(m),idas(itriv(m)),das(itrrv(m)),
     &              itxnm(m),itynm(m),itznm(m),
     &              das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &              das_itzrg(itzrg(m)),
     &              igsh,idasa)

               else if( itals(m) .eq. 42 ) then

                     idas0 = nmmax
                     idas1 = lmmax
                     idas2 = idas1 + itenm(m)
                     idas3 = idas2 + ittnm(m)
                     idas4 = ( idas3 + itrgn(m) - 1 ) * 2 + 1
                     idas5 = idas3 + itrgn(m)
     &                     + ( itrgn(m) + mod(itrgn(m),2) ) / 2
                     idas6 = ( idas5 + itrgn(m) - 1 ) * 2 + 1
                     idasa = idas5

                  call pwwgreg(m,
     &              itpan(m),itrgn(m),itrgm(m),itenm(m),itmst(m),
     &              ittnm(m),
     &              idas_itreg(itreg(m)),
     &              das_iterg(iterg(m)),das_ittrg(ittrg(m)),tr,
     &              itrnv(m),idas(itriv(m)),das(itrrv(m)),
     &              itxnm(m),itynm(m),itznm(m),mftal(m), ! S.H. 2024.3.21
     &              das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &              das_itzrg(itzrg(m)),
     &              igsh,idasa)


               else if( itals(m) .eq. 48 ) then

                     idas0 = nmmax
                     idas1 = lmmax
                     idas2 = idas1 + itenm(m)
                     idas3 = ( idas2 + ittnm(m) - 1 ) * 2 + 1
                     idasa = idas2 + ittnm(m)

                  call pwwgxyz(m,
     &              itpan(m),itmtn(m),ismte(itmtt(m)),
     &              itxnm(m),itynm(m),itznm(m),itenm(m),itmst(m),
     &              ittnm(m),mftal(m), ! S.H. 2024.3.21
     &              das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &              das_itzrg(itzrg(m)),
     &              das_iterg(iterg(m)),das_ittrg(ittrg(m)),tr,
     &              igsh,idasa)

               end if

                  if(icgg.ne.0) call DEALLOCATE_GGBANK !FURUTA

            end if

         end do

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine dshows
*                                                                      *
*       output 3dshow tally only                                       *
*       last modified by K.Niita on 2002/10/24                         *
*                                                                      *
************************************************************************
      use GGBANKMOD !FURUTA
      use moddas_region
      use moddas_tally
*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

      include 'param.inc'

*-----------------------------------------------------------------------

      common /talmm/  nmmax, lmmax, itlmx
      common /tall00/ itnm, ital(itlmax), itals(itlmax), italm(itlmax)
      common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                itout(itlmax), ittwo(itlmax)

      common /tall22/ itmcn(itlmax), itmtn(itlmax), itmtt(itlmax)

      common /tall02/ itrgn(itlmax), itrgm(itlmax), itreg(itlmax)
      common /tall25/ itrnv(itlmax), itriv(itlmax), itrrv(itlmax)

      common /tall40/ rtorg(itlmax,3), rteye(itlmax,3), rtlit(itlmax,3),
     &                rtwin(itlmax,3), rtbrt(itlmax,2), ithvn(itlmax),
     &                itwin(itlmax,2), itbox(itlmax), itmir(itlmax),
     &                rtbox(itlmax,5,10), rtout(itlmax), rthet(itlmax),
     &                itlin(itlmax), itshd(itlmax), itgxs(itlmax)

      common /tall42/ itmbn(itlmax), itmbt(itlmax)
      common /tall43/ itmgn(itlmax), itmgm(itlmax), itmeg(itlmax)

      dimension     idas(1)
      equivalence ( das, idas )

      common /ccggg/  icgg

*-----------------------------------------------------------------------

         do m = 1, itnm

            if( itals(m) .eq. 27 ) then

                     nx = itwin(m,1) + 1
                     ny = itwin(m,2) + 1

              if(icgg.ne.0)then
                  call ALLOCATE_GGBANK !FURUTA
                  call INIT_GGBANK     !FURUTA
              endif

                     idas1 = ( lmmax - 1 ) * 2 + 1
                     idas2 = idas1 + nx * ny * 2
                     idas3 = idas2 + nx * ny
                     idas4 = idas3 + nx * ny * 2
                     idasa = lmmax

                  call pdshow(m,
     &               itout(m),itmtn(m),ismte(itmtt(m)),
     &               itmbn(m),ismte_itmbt(itmbt(m)),
     &               itrgn(m),itrgm(m),idas_itreg(itreg(m)),
     &               itmgn(m),itmgm(m),idas_itmeg(itmeg(m)),
     &               nx,ny,idasa)

                  if(icgg.ne.0) call DEALLOCATE_GGBANK !FURUTA

            end if

         end do

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      function vls(nl,lt,ilt,nsm,itrns,xxs0,xxf0,yys0,yyf0,zzs0,zzf0)
*                                                                      *
*       calculate the volume of xyz elements with certain material     *
*       last modified by K.Niita on 2004/12/27                         *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'

*-----------------------------------------------------------------------

      parameter( eps  = 1.0d-08 )
      parameter( epss = 1.0d-04 )

*-----------------------------------------------------------------------

      common /kmat1d/ idmn(0:kvlmax), idnm(kvmmax)

*-----------------------------------------------------------------------

      dimension   lt(nl)

      dimension   xm(1000)
      dimension   ym(1000)
      dimension   zm(1000)

      dimension   xx(2)
      dimension   yy(2)
      dimension   zz(2)

      dimension   nmedn(8), isn(2), uvw(3,3)

      data isn / 1, -1 /

*-----------------------------------------------------------------------

            dx = xxf0 - xxs0
            dy = yyf0 - yys0
            dz = zzf0 - zzs0

            vls0 = dx * dy * dz

*-----------------------------------------------------------------------
*     no material selection
*-----------------------------------------------------------------------

         if( nl .eq. 0 .or. nsm .eq. 0 ) goto 1000

*-----------------------------------------------------------------------
*     transform
*-----------------------------------------------------------------------

         call trnsxv(xxs0,yys0,zzs0,xxs,yys,zzs,itrns)
         call trnsxv(xxf0,yyf0,zzf0,xxf,yyf,zzf,itrns)

         call trnsuv(1.d0,0.d0,0.d0,uvw(1,1),uvw(1,2),uvw(1,3),itrns)
         call trnsuv(0.d0,1.d0,0.d0,uvw(2,1),uvw(2,2),uvw(2,3),itrns)
         call trnsuv(0.d0,0.d0,1.d0,uvw(3,1),uvw(3,2),uvw(3,3),itrns)

*-----------------------------------------------------------------------

            nsn = abs( nsm )

            xx(1) = xxs
            xx(2) = xxf
            yy(1) = yys
            yy(2) = yyf
            zz(1) = zzs
            zz(2) = zzf

*-----------------------------------------------------------------------
*     check only of the eight corners only for nsm > 0
*-----------------------------------------------------------------------

      if( nsm .gt. 0 ) then

                  in = 0

            do ix = 1, 2
            do iy = 1, 2
            do iz = 1, 2

                  in = in + 1

                  ici   = -1
                  mark  =  1
                  markp =  0

                  xi = xx(ix) + uvw(1,1) * isn(ix) * eps
                  yi = yy(iy) + uvw(2,2) * isn(iy) * eps
                  zi = zz(iz) + uvw(3,3) * isn(iz) * eps

                  u = uvw(1,1) * isn(ix)
                  v = uvw(2,2) * isn(iy)
                  w = uvw(3,3) * isn(iz)

                  rnm = sqrt( u**2 + v**2 + w**2 )
                  u = u / rnm
                  v = v / rnm
                  w = w / rnm

                  call gomsor(xi,yi,zi,u,v,w,
     &                        nmedn(in),iblz,mark,markp,ici)

                  if( mark .lt. -1 ) goto 1000

            end do
            end do
            end do

                  ichk = 0

               do i = 1, 8

                     do k = 1, nl

                        if( ( ilt .gt. 0 .and.
     &                        idmn(nmedn(i)) .eq. lt(k) ) .or.
     &                      ( ilt .lt. 0 .and.
     &                        idmn(nmedn(i)) .ne. lt(k) ) ) then

                           ichk = ichk + 1

                           goto 2000

                        end if

                     end do

 2000                continue

               end do

               if( ichk .eq. 8 .or. ichk .eq. 0 ) goto 1000

      end if

*-----------------------------------------------------------------------
*     with material selection
*-----------------------------------------------------------------------

               ddx = dx / nsn
               ddy = dy / nsn
               ddz = dz / nsn

            do i = 1, nsn

               xm(i) = xx(1) + ddx / 2.0 + dble(i-1) * ddx
               ym(i) = yy(1) + ddy / 2.0 + dble(i-1) * ddy
               zm(i) = zz(1) + ddz / 2.0 + dble(i-1) * ddz

            end do

               xs = xx(1) + uvw(1,1) * ddx * eps
               xf = xx(2) - uvw(1,1) * ddx * eps
               ys = yy(1) + uvw(2,2) * ddy * eps
               yf = yy(2) - uvw(2,2) * ddy * eps
               zs = zz(1) + uvw(3,3) * ddz * eps
               zf = zz(2) - uvw(3,3) * ddz * eps

               voly = 0.0
               volx = 0.0
               volz = 0.0

*-----------------------------------------------------------------------
*     scan z direction
*-----------------------------------------------------------------------

                  u = uvw(3,1)
                  v = uvw(3,2)
                  w = uvw(3,3)

         do 100 i = 1, nsn
         do 200 j = 1, nsn

                  xi = xm(i)
                  yi = ym(j)
                  zi = zs

                  ici   = -1
                  mark  =  1
                  markp =  0

                  call gomsor(xi,yi,zi,u,v,w,
     &                        nmed,iblz,mark,markp,ici)

               if( mark .ge. -1 ) then

                  nmed0 = nmed
                  mark  = 1
                  markp = 1

               else

                  goto 1000

               end if

*-----------------------------------------------------------------------

  300       continue

                  xc = xm(i)
                  yc = ym(j)
                  zc = zf

                  call gomprp(0,xi,yi,zi,xc,yc,zc,u,v,w,
     &                        nmed,iblz,mark,markp)

                        if( mark .lt. -1 ) goto 1000

                  if( nmed0 .gt. 0 ) then

                     do k = 1, nl

                        if( ( ilt .gt. 0 .and.
     &                        idmn(nmed0) .eq. lt(k) ) .or.
     &                      ( ilt .lt. 0 .and.
     &                        idmn(nmed0) .ne. lt(k) ) ) then

                           volz = volz + zc - zi

                           goto 400

                        end if

                     end do

  400                continue

                  end if

                  if( mark .eq. 1 .or. mark .eq. -1 ) then

                           goto 200

                  else if( mark .eq. 0 .or. mark .eq. 2 ) then

                           nmed0 = nmed

                           zi = zc

                           goto 300

                  end if

*-----------------------------------------------------------------------

  200    continue
  100    continue

               volz = volz * ddx * ddy

*-----------------------------------------------------------------------
*        no boundary
*-----------------------------------------------------------------------

         if( volz .eq. 0. .or. abs( volz - vls0 ) .lt. epss ) goto 1000

*-----------------------------------------------------------------------
*     scan x direction
*-----------------------------------------------------------------------

                  u = uvw(1,1)
                  v = uvw(1,2)
                  w = uvw(1,3)

         do 110 i = 1, nsn
         do 210 j = 1, nsn

                  xi = xs
                  yi = ym(i)
                  zi = zm(j)

                  ici   = -1
                  mark  =  1
                  markp =  0

                  call gomsor(xi,yi,zi,u,v,w,
     &                        nmed,iblz,mark,markp,ici)

               if( mark .ge. -1 ) then

                  nmed0 = nmed
                  mark  = 1
                  markp = 1

               else

                  goto 1000

               end if

*-----------------------------------------------------------------------

  310       continue

                  xc = xf
                  yc = ym(i)
                  zc = zm(j)

                  call gomprp(0,xi,yi,zi,xc,yc,zc,u,v,w,
     &                        nmed,iblz,mark,markp)

                        if( mark .lt. -1 ) goto 1000

                  if( nmed0 .gt. 0 ) then

                     do k = 1, nl

                        if( ( ilt .gt. 0 .and.
     &                        idmn(nmed0) .eq. lt(k) ) .or.
     &                      ( ilt .lt. 0 .and.
     &                        idmn(nmed0) .ne. lt(k) ) ) then

                           volx = volx + xc - xi

                           goto 410

                        end if

                     end do

  410                continue

                  end if

                  if( mark .eq. 1 .or. mark .eq. -1 ) then

                           goto 210

                  else if( mark .eq. 0 .or. mark .eq. 2 ) then

                           nmed0 = nmed

                           xi = xc

                           goto 310

                  end if

*-----------------------------------------------------------------------

  210    continue
  110    continue

               volx = volx * ddy * ddz

*-----------------------------------------------------------------------
*        no boundary
*-----------------------------------------------------------------------

         if( volx .eq. 0. .or. abs( volx - vls0 ) .lt. epss ) goto 1000

*-----------------------------------------------------------------------
*     scan y direction
*-----------------------------------------------------------------------

                  u = uvw(2,1)
                  v = uvw(2,2)
                  w = uvw(2,3)

         do 120 i = 1, nsn
         do 220 j = 1, nsn

                  xi = xm(i)
                  yi = ys
                  zi = zm(j)

                  ici   = -1
                  mark  =  1
                  markp =  0

                  call gomsor(xi,yi,zi,u,v,w,
     &                        nmed,iblz,mark,markp,ici)

               if( mark .ge. -1 ) then

                  nmed0 = nmed
                  mark  = 1
                  markp = 1

               else

                  goto 1000

               end if

*-----------------------------------------------------------------------

  320       continue

                  xc = xm(i)
                  yc = yf
                  zc = zm(j)

                  call gomprp(0,xi,yi,zi,xc,yc,zc,u,v,w,
     &                        nmed,iblz,mark,markp)

                        if( mark .lt. -1 ) goto 1000

                  if( nmed0 .gt. 0 ) then

                     do k = 1, nl

                        if( ( ilt .gt. 0 .and.
     &                        idmn(nmed0) .eq. lt(k) ) .or.
     &                      ( ilt .lt. 0 .and.
     &                        idmn(nmed0) .ne. lt(k) ) ) then

                           voly = voly + yc - yi

                           goto 420

                        end if

                     end do

  420                continue

                  end if

                  if( mark .eq. 1 .or. mark .eq. -1 ) then

                           goto 220

                  else if( mark .eq. 0 .or. mark .eq. 2 ) then

                           nmed0 = nmed

                           yi = yc

                           goto 320

                  end if

*-----------------------------------------------------------------------

  220    continue
  120    continue

               voly = voly * ddx * ddz

*-----------------------------------------------------------------------
*        no boundary
*-----------------------------------------------------------------------

         if( voly .eq. 0. .or. abs( voly - vls0 ) .lt. epss ) goto 1000

*-----------------------------------------------------------------------
*        average of three scans
*-----------------------------------------------------------------------

            vls = ( volx + voly + volz ) / 3.0

            if( vls .lt. vls0 / dble(nsn)**2 ) goto 1000

*-----------------------------------------------------------------------

            return

 1000 continue

            vls = vls0

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine sumtr(tr,br,nr,ncol)
*                                                                      *
*        sum tally values ( tr <- br + tr )                            *
*                                                                      *
*        by D.OBINATA on 2012.6.6                                      *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      common /talout/ itall ! T.Sato 2017/08/05
      common /imagecom/ imageout ! T.Sato 2017/08/05
      common /tcntl/  icntl, inucr ! T.Sato 2017/09/27

*-----------------------------------------------------------------------

      dimension tr(nr), br(nr)

*-----------------------------------------------------------------------

      tr(1:nr) = br(1:nr) + tr(1:nr)

! T.Sato 2019/01/06, activate itall for MPI
      if(itall.eq.0.and.ncol.ne.2) then ! output image (eps,vtk,bmp) or not
       imageout=1 ! do not output image
      else
       imageout=0
      endif
      if(icntl.eq.1.or.icntl.eq.13) imageout=0 ! for icntl = 1 or 13, always output eps, 2017/09/27

      return
      end subroutine


************************************************************************
*                                                                      *
      module tdepwgtsum_local
*                                                                      *
*     Modules for reading parameters of weighted summation             *
*     in [T-Deposit]                                                   *
*                                                                      *
*     Following subroutines are called externally                      *
*        allocate_depwgtsum01                                          *
*        allocate_depwgtsum02                                          *
*        reallocate_depwgtsum01                                        *
*        deallocate_depwgtsum                                          *
*                                                                      *
*     Following subroutines use this module                            *
*        tdeposit      tallsm2.f                                       *
*         tdepreg      tallsm1.f                                       *
*        tdepreg0      result.f                                        *
*                                                                      *
*     made by S.Abe on 2016/11/24                                      *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     common arrays                                                    *
*                                                                      *
*    icond(1,i,j)      condition number of i-th line and j-th add      *
*    icond(2,i,j)      operator for i-th linen and j-th add            *
*    icond(3,i,j)      list number for i-th line in condition          *
*     ethres(i,j)      threshold energy of condition                   *
*                                                                      *
*      numlist(l)      list number for l-th index in efficiency        *
*     depeff(k,l)      efficiency for k-th cell in l-th list           *
*                                                                      *
************************************************************************

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      integer,allocatable:: icond(:,:,:)
      real(8),allocatable:: ethres(:,:)
      integer,allocatable:: numlist(:)
      real(8),allocatable:: depeff(:,:)

*-----------------------------------------------------------------------
c temporal storage for reallocation
      integer,allocatable:: idummy(:,:,:)
      real(8),allocatable:: rdummy(:,:)

*-----------------------------------------------------------------------

      contains

************************************************************************
*                                                                      *
      subroutine allocate_depwgtsum01(ncond,nadd)
*                                                                      *
************************************************************************

      allocate( icond(3,ncond,nadd) )
      allocate( ethres(ncond,nadd) )
      icond(:,:,:) = 0
      ethres(:,:) = 0.d0

*-----------------------------------------------------------------------
      return
      end subroutine

************************************************************************
*                                                                      *
      subroutine allocate_depwgtsum02(ncell,nefflist)
*                                                                      *
************************************************************************

      allocate( numlist(nefflist) )
      allocate( depeff(ncell,0:nefflist) )
      numlist(:) = -1
      depeff(:,:) = 0.d0

*-----------------------------------------------------------------------
      return
      end subroutine

************************************************************************
*                                                                      *
      subroutine reallocate_depwgtsum01(ncond,nadd,madd)
*                                                                      *
************************************************************************

      allocate( idummy(3,ncond,nadd) )
      allocate( rdummy(ncond,nadd) )
      idummy(:,:,:) = icond(:,:,:)
      rdummy(:,:) = ethres(:,:)

      deallocate( icond )
      deallocate( ethres )
      allocate( icond(3,ncond,madd) )
      allocate( ethres(ncond,madd) )
      icond(:,:,:) = 0
      ethres(:,:) = 0.d0

      do j = 1, ncond
       do k = 1, nadd
        do i = 1, 3
         icond(i,j,k) = idummy(i,j,k)
        enddo
        ethres(j,k) = rdummy(j,k)
       enddo
      enddo

      deallocate( idummy )
      deallocate( rdummy )

*-----------------------------------------------------------------------
      return
      end subroutine

************************************************************************
*                                                                      *
      subroutine deallocate_depwgtsum
*                                                                      *
************************************************************************

      deallocate( icond )
      deallocate( ethres )
      deallocate( numlist )
      deallocate( depeff )

*-----------------------------------------------------------------------
      return
      end subroutine


      endmodule tdepwgtsum_local
************************************************************************



************************************************************************
*                                                                      *
      module tdepwgtsum_global
*                                                                      *
*     Modules for storaging parameters of weighted summation           *
*     in [T-Deposit]                                                   *
*                                                                      *
*     Following subroutines are called externally                      *
*        allocate_depwgtsum03                                          *
*        reallocate_depwgtsum03                                        *
*        reallocate_depwgtsum04                                        *
*                                                                      *
*     Following subroutines use this module                            *
*       tdepstreg      talls06.f                                       *
*  tdepstregEVENT      talls08.f                                       *
*        tdeposit      tallsm2.f                                       *
*    echrg_wgtsum      tallsm3.f                                       *
*                                                                      *
*     made by S.Abe on 2016/11/24                                      *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     common arrays                                                    *
*                                                                      *
* itcond(m,1,i,j)      condition number of i-th line and j-th add      *
* itcond(m,2,i,j)      operator for i-th linen and j-th add            *
* itcond(m,3,i,j)      list number for i-th line in condition          *
*    rteth(m,i,j)      threshold energy of condition                   *
*                                                                      *
*     itlist(m,l)      list number for l-th index in efficiency        *
*    rteff(m,k,l)      efficiency for k-th cell in l-th list           *
*                                                                      *
************************************************************************

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      integer,save:: allocate_first=0
      integer,save:: ncond_old=1
      integer,save:: nadd_old=1
      integer,save:: ncell_old=1
      integer,save:: nlist_old=1

*-----------------------------------------------------------------------

      integer,allocatable:: itcond(:,:,:,:)
      real(8),allocatable:: rteth(:,:,:)
      integer,allocatable:: itlist(:,:)
      real(8),allocatable:: rteff(:,:,:)

*-----------------------------------------------------------------------
c temporal storage for reallocation
      integer,allocatable:: idummy1(:,:,:,:)
      real(8),allocatable:: rdummy1(:,:,:)
      integer,allocatable:: idummy2(:,:)
      real(8),allocatable:: rdummy2(:,:,:)

*-----------------------------------------------------------------------
      contains

************************************************************************
*                                                                      *
      subroutine allocate_depwgtsum03(ncond_old,nadd_old,
     &                              ncell_old,nlist_old)
*                                                                      *
************************************************************************

      include 'param.inc'

      allocate( itcond(itlmax,3,ncond_old,nadd_old) )
      allocate( rteth(itlmax,ncond_old,nadd_old) )
      allocate( itlist(itlmax,nlist_old) )
      allocate( rteff(itlmax,ncell_old,0:nlist_old) )
      itcond(:,:,:,:) = 0
      rteth(:,:,:) = 0.d0
      itlist(:,:) = -1
      rteff(:,:,:) = 0.d0

*-----------------------------------------------------------------------
      return
      end subroutine

************************************************************************
*                                                                      *
      subroutine reallocate_depwgtsum03(ncond,nadd,ncond_old,nadd_old)
*                                                                      *
************************************************************************

      include 'param.inc'

      allocate( idummy1(itlmax,3,ncond_old,nadd_old) )
      allocate( rdummy1(itlmax,ncond_old,nadd_old) )
      idummy1(:,:,:,:) = itcond(:,:,:,:)
      rdummy1(:,:,:) = rteth(:,:,:)

      deallocate( itcond )
      deallocate( rteth )
      allocate( itcond(itlmax,3,ncond,nadd) )
      allocate( rteth(itlmax,ncond,nadd) )
      itcond(:,:,:,:) = 0
      rteth(:,:,:) = 0.d0

      do m = 1, itlmax
       do j = 1, ncond_old
        do k = 1, nadd_old
         do i = 1, 3
          itcond(m,i,j,k) = idummy1(m,i,j,k)
         enddo
         rteth(m,j,k) = rdummy1(m,j,k)
        enddo
       enddo
      enddo

      deallocate( idummy1 )
      deallocate( rdummy1 )

*-----------------------------------------------------------------------
      return
      end subroutine

************************************************************************
*                                                                      *
      subroutine reallocate_depwgtsum04(ncell,nefflist,
     &                                ncell_old,nlist_old)
*                                                                      *
************************************************************************

      include 'param.inc'

      allocate( idummy2(itlmax,nlist_old) )
      allocate( rdummy2(itlmax,ncell_old,0:nlist_old) )
      idummy2(:,:) = itlist(:,:)
      rdummy2(:,:,:) = rteff(:,:,:)

      deallocate( itlist )
      deallocate( rteff )
      allocate( itlist(itlmax,nefflist) )
      allocate( rteff(itlmax,ncell,0:nefflist) )
      itlist(:,:) = -1
      rteff(:,:,:) = 0.d0

      ijk = itlmax
      do m = 1, itlmax
       do j = 1, nlist_old
        itlist(m,j) = idummy2(m,j)
       enddo
      enddo
      do m = 1, itlmax
       do i = 1, ncell_old
        do j = 0, nlist_old
         rteff(m,i,j) = rdummy2(m,i,j)
        enddo
       enddo
      enddo

      deallocate( idummy2 )
      deallocate( rdummy2 )

*-----------------------------------------------------------------------
      return
      end subroutine

      endmodule tdepwgtsum_global


************************************************************************
*                                                                      *
      subroutine inittal_extstat
*                                                                      *
*       initialization of tr for extended statistical indicators       *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
************************************************************************
      use EVENTTALMOD !FURUTA
      use TALMOD
!$       use TALMOD0
      use sumtallymod !Hashimoto (2015.1.26)
      use MMBANKMOD !S.H. for tally unit of MeV/n (2016.12.25)
      use partmod, only: itpan, itpat, jtpat ! frtati 2021/10/05
      use moddas_mesh
      use moddas_region
      use moddas_tally

*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'param01.inc'

*-----------------------------------------------------------------------

      common /mpi00/ npe, me

      common /talout/ itall
      common /talmm/  nmmax, lmmax, itlmx

*-----------------------------------------------------------------------

      common /tall00/ itnm, ital(itlmax), itals(itlmax), italm(itlmax)
      common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                itout(itlmax), ittwo(itlmax)

      common /tall02/ itrgn(itlmax), itrgm(itlmax), itreg(itlmax)
      common /tall25/ itrnv(itlmax), itriv(itlmax), itrrv(itlmax)

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

c      common /tall09/ itpan(itlmax), itpat(itlmax,6,2),
c     &                jtpat(itlmax,6,6,2) ! frtati 2021/10/05

      common /tall12/ itman(itlmax), itmat(itlmax), itmct(itlmax),
     &                itnun(itlmax), itnuc(itlmax),
     &                itndz(itlmax), itndn(itlmax),
     &                itnkz(itlmax), itnkn(itlmax)

      common /tall13/ itrcn(itlmax), itrcr(itlmax), itrca(itlmax)
      common /tall26/ itrcs(itlmax), itrcc(itlmax), itrss(itlmax)

      common /tall16/ itnfr(itlmax), itnfz(itlmax), itnfn(itlmax)
      common /tall17/ itndy(itlmax)
      common /tall18/ ithet(itlmax)

      common /tall20/ ittty(itlmax), ittnm(itlmax), ittrg(itlmax),
     &                rttmi(itlmax), rttma(itlmax), rttdl(itlmax)

      common /tall22/ itmcn(itlmax), itmtn(itlmax), itmtt(itlmax)

      common /tall23/ itaty(itlmax), itanm(itlmax), itarg(itlmax),
     &                rtami(itlmax), rtama(itlmax), rtadl(itlmax)
cKN 2016/12/31
      common /tall24/ itrcm(itlmax), itrcg(itlmax)

      common /tall45/ itmlp(itlmax), itmln(itlmax,6), itmst(itlmax),
     &                itmli(itlmax,6), rtmme(itlmax,6), itmnt(itlmax,6),
     &                itmpn(itlmax,6), itmpt(itlmax,6,6,2)
      common /tall48/ itmdp(itlmax,0:30), itmdf(itlmax)

      common /tall56/ itety2(itlmax), itenm2(itlmax), iterg2(itlmax),
     &                rtemi2(itlmax), rtema2(itlmax), rtedl2(itlmax)

cKN 2016/12/31
      common /tall67/ itstd(itlmax), rtstd(itlmax)
      common / stdstop / idstop, instop

cABE 2018/02/08
      common /tall72/ itactnm(itlmax), itactrg(itlmax), itactmax(itlmax)

cABE 2018/10/18, 2019/11/26
      common /tall73/ itenclo(itlmax), itangform(itlmax)

cfrtati 2022/02/18
      common /tall83/ itnzn(itlmax), itndm(itlmax)

      common /tall91/ itextstat(itlmax), mftal(itlmax) ! S.H. extstat 2024.3.18
      common /tall93/ tprodenmn(itlmax), tprodenmx(itlmax),
     &     nbtproden(itlmax), itprodenchk(itlmax) !S.H. extstat 2024.4.28

      dimension     idas(mdas*2) !20220906frtati temporary for fbounds-check
      equivalence ( das, idas )

      real(8),pointer:: tr(:)

*-----------------------------------------------------------------------

C --- add NS 2020.04 ---
C for nonshared_tally option
      integer italsh
      common /talsh/ italsh
C --- end add NS 2020.04 ---

*-----------------------------------------------------------------------

      do m = 1, itnm

       if( mftal(m) .gt. 5 ) then

C --- add NS 2020.04 ---
C for nonshared_tally option
!$          if(italsh .eq. 0) then
!$            call GET_TR_HEAD_POINTER0(tr,m)
!$          else
C --- end add NS 2020.04 ---
            call GET_TR_HEAD_POINTER(tr,m)
C --- add NS 2020.04 ---
!$          endif
C --- end add NS 2020.04
*-----------------------------------------------------------------------

            if( itals(m) .eq. 1 ) then

               call inittal_extstat_ttracreg(m,
     &              itpan(m),itrgn(m),itrgm(m),itenm(m),itmst(m),
     &              ittnm(m),mftal(m),
     &              idas_itreg(itreg(m)),
     &              das_iterg(iterg(m)),das_ittrg(ittrg(m)),tr,
     &              itrnv(m),idas(itriv(m)),das(itrrv(m)),
     &              tprodenmn(m),tprodenmx(m),nbtproden(m))

*-----------------------------------------------------------------------

            else if( itals(m) .eq. 2 ) then

               call inittal_extstat_ttracrz(m,
     &              itpan(m),itrnm(m),itznm(m),itenm(m),itmst(m),
     &              ittnm(m),itanm(m),mftal(m),
     &              das_itrrg(itrrg(m)),das_itzrg(itzrg(m)),
     &              das_itarg(itarg(m)),
     &              das_iterg(iterg(m)),das_ittrg(ittrg(m)),tr,
     &              tprodenmn(m),tprodenmx(m),nbtproden(m))

*-----------------------------------------------------------------------

            else if( itals(m) .eq. 3 ) then

               call inittal_extstat_ttracxyz(m,
     &              itmtn(m),ismte(itmtt(m)),
     &              itpan(m),itxnm(m),itynm(m),itznm(m),itenm(m),
     &              itmst(m),ittnm(m),mftal(m),
     &              das_itxrg(itxrg(m)),das_ityrg(ityrg(m)),
     &              das_itzrg(itzrg(m)),
     &              das_iterg(iterg(m)),das_ittrg(ittrg(m)),tr,
     &              tprodenmn(m),tprodenmx(m),nbtproden(m))

*-----------------------------------------------------------------------

            else if( itals(m) .eq. 41 ) then

               call inittal_extstat_tpointp(m,
     &              itpan(m),itmsh(m),itenm(m),itmst(m),ittnm(m),
     &              mftal(m),
     &              das_iterg(iterg(m)),das_ittrg(ittrg(m)),tr,
     &              tprodenmn(m),tprodenmx(m),nbtproden(m))

*-----------------------------------------------------------------------

            end if

         end if

      end do

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine sumtr_sumover(m,ncol)
*                                                                      *
************************************************************************

      use TALMOD
!$      use TALMOD0
      use RESTALMOD

      implicit double precision (a-h,o-z)

C for nonshared_tally option
      integer italsh
      common /talsh/ italsh

      real(8),pointer :: tr_sum(:)
      real(8),pointer :: br_sum(:)

      iax = 1

C for nonshared_tally option
!$       if(italsh .eq. 0) then
!$          call GET_TR_HEAD_POINTER0_SUM(tr_sum,m,iax)
!$       else
            call GET_TR_HEAD_POINTER_SUM(tr_sum,m,iax)
C for nonshared_tally option
!$       end if

      call GET_TRRES_HEAD_POINTER_SUM(br_sum,m,iax)

      n_tr_sum = mtalsize_sum(m)

      call sumtr_sumover_sub(tr_sum,br_sum,n_tr_sum,ncol)

      return
      end

************************************************************************
*                                                                      *
      subroutine sumtr_sumover_sub(tr,br,nr,ncol)
*                                                                      *
*        sum tally values ( tr <- br + tr )                            *
*                                                                      *
*        by D.OBINATA on 2012.6.6                                      *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)
      
      dimension tr(nr), br(nr)

      call sumtr(tr,br,nr,ncol)

      return
      end subroutine

      subroutine print_tally_sum(m,iax)

      use TALMOD
!$      use TALMOD0

      implicit double precision (a-h,o-z)

C for nonshared_tally option
      integer italsh
      common /talsh/ italsh

      real(8),pointer :: p_sum(:)

C for nonshared_tally option
!$       if(italsh .eq. 0) then
!$          n_tr_sum = italsize0_sum(m,iax)
!$          call GET_TR_HEAD_POINTER0_SUM(p_sum,m,iax)
!$       else
            n_tr_sum = italsize_sum(m,iax)
            call GET_TR_HEAD_POINTER_SUM(p_sum,m,iax)
C for nonshared_tally option
!$       end if

      call print_tally_sum_sub(n_tr_sum,p_sum,m,iax)

      return
      end

      subroutine print_tally_sum_sub(n_tr_sum,tr_sum,m,iax)

      implicit double precision (a-h,o-z)

      real(8) :: tr_sum(n_tr_sum)

      write(*,*) 'm,iax,tr=',m,iax,(tr_sum(i),i=1,n_tr_sum)


      end

      subroutine print_tallyres_sum(m,iax)

      use RESTALMOD

      implicit double precision (a-h,o-z)

      real(8),pointer :: p_sum(:)

      n_tr_sum = lrestalm_sum(m,iax)

      p_sum =>  trRES_sum(irestalm_sum(m,iax):)

      call print_tallyres_sum_sub(n_tr_sum,p_sum,m,iax)

      return
      end

      subroutine print_tallyres_sum_sub(n_tr_sum,tr_sum,m,iax)

      implicit double precision (a-h,o-z)

      real(8) :: tr_sum(n_tr_sum)

      write(*,*) 'm,iax,trRES=',m,iax,(tr_sum(i),i=1,n_tr_sum)


      end

      subroutine print_tallysum_sum(m,iax)

      use RESTALMOD
      use sumtallymod


      implicit double precision (a-h,o-z)

      real(8),pointer :: p_sum(:)

      n_tr_sum = lrestalm_sum(m,iax)

      p_sum =>  trSUMTAL_SUM(irestalm_sum(m,iax):)

      call print_tallyres_sum_sub(n_tr_sum,p_sum,m,iax)

      return
      end

      subroutine print_tallysum_sum_sub(n_tr_sum,tr_sum,m,iax)

      implicit double precision (a-h,o-z)

      real(8) :: tr_sum(n_tr_sum)

      write(*,*) 'm,iax,trRES=',m,iax,(tr_sum(i),i=1,n_tr_sum)


      end


************************************************************************
*                                                                      *
      subroutine tr_sum_to_tott(ntf, ndim, nstep, istart, ibase,
     &                          n_sum_l,n_sum, maxtott, tott, tr_sum)
*
* ntf          : file number in nfiles
* ndim         : counts of parameters
* nstep(ndim)  : do step number
* istart(ndim) : start of do loop
* ibase(ndime) : base of do loop
* n_sum        : length of each parameter
* maxtott      : tott(maxtott,2)
* tott         : output data
* tr_sum(n_sum(1),-,n_sum(ndim),2*nfles) : database
*                                                                     *
************************************************************************

      implicit double precision (a-h,o-z)

      integer :: ntf, ndim, nstep(*), istart(*), ibase(*),
     &           n_sum_l(*),n_sum(*)
      integer :: maxtott
      real(8) :: tott(maxtott,2)
      real(8) :: tr_sum(*)

      integer :: n_sum_b(11), ist(11), ien(11), iloop(11), kbase,
     &           lb(11), ib(11), ia(11)

      kbase = 0
      n_sum_b(1) = 1
      do i=2,ndim+1
        n_sum_b(i) = n_sum_b(i-1)*(n_sum(i-1) + 1 - n_sum_l(i-1))
      enddo
      kbase = n_sum_b(ndim+1)
      n_sum_b(1) = 0

      lb(:) = 1
      ib(:) = 1
      ist(:) = 1
      do i=1,ndim
        lb(i) = n_sum_l(i)
        ib(i) = ibase(i)
        ist(i) = istart(i)
      enddo
      ien(:) = ist(:)
      do i=1,ndim
        ien(i) = nstep(i)
      enddo


      i = 0
      do ia10=ist(10),ien(10)
         ia(10) = ia10
       do ia9=ist(9),ien(9)
          ia(9) = ia9
        do ia8=ist(8),ien(8)
           ia(8) = ia8
         do ia7=ist(7),ien(7)
            ia(7) = ia7
          do ia6=ist(6),ien(6)
             ia(6) = ia6
           do ia5=ist(5),ien(5)
              ia(5) = ia5
            do ia4=ist(4),ien(4)
               ia(4) = ia4
             do ia3=ist(3),ien(3)
                ia(3) = ia3
              do ia2=ist(2),ien(2)
                 ia(2) = ia2
               do ia1=ist(1),ien(1)
                  ia(1) = ia1
                  i = i + 1
                  do k=1,2
                    call get_tr_sum_address
     &                   (ntf,k,ndim,lb,ib,ia,n_sum_b,kbase,ig)
                    tott(i,k) = tr_sum(ig)
                  enddo
               enddo
              enddo
             enddo
            enddo
           enddo
          enddo
         enddo
        enddo
       enddo
      enddo

      return
      end

************************************************************************
*                                                                      *
      subroutine tott_to_tr_sum(ntf, ndim, nstep, istart, ibase,
     &                          n_sum_l, n_sum,
     &                          maxtott, tott, tr_sum)
*
* ntf          : file number in nfiles
* ndim         : counts of parameters
* nstep(ndim)  : do step number
* istart(ndim) : start of do loop
* ibase(ndime) : base of do loop
* n_sum        : length of each parameter
* maxtott      : tott(maxtott,2)
* tott         : set data
* tr_sum(n_sum(1),-,n_sum(ndim),2*nfles) : database
*                                                                     *
************************************************************************

      implicit double precision (a-h,o-z)

      integer :: ntf, ndim, nstep(*), istart(*), ibase(*),
     &          n_sum_l(*), n_sum(*)
      integer :: maxtott
      real(8) :: tott(maxtott,2)
      real(8) :: tr_sum(*)

      integer :: n_sum_b(11), ist(11), ien(11), iloop(11), kbase,
     &           lb(11), ib(11), ia(11)

      kbase = 0
      n_sum_b(1) = 1
      do i=2,ndim+1
        n_sum_b(i) = n_sum_b(i-1)*(n_sum(i-1) + 1 - n_sum_l(i-1))
      enddo
      kbase = n_sum_b(ndim+1)
      n_sum_b(1) = 0

      lb(:) = 1
      ib(:) = 1
      ist(:) = 1
      do i=1,ndim
        lb(i) = n_sum_l(i)
        ib(i) = ibase(i)
        ist(i) = istart(i)
      enddo
      ien(:) = ist(:)
      do i=1,ndim
        ien(i) = nstep(i)
      enddo


      i = 0
      do ia10=ist(10),ien(10)
         ia(10) = ia10
       do ia9=ist(9),ien(9)
          ia(9) = ia9
        do ia8=ist(8),ien(8)
           ia(8) = ia8
         do ia7=ist(7),ien(7)
            ia(7) = ia7
          do ia6=ist(6),ien(6)
             ia(6) = ia6
           do ia5=ist(5),ien(5)
              ia(5) = ia5
            do ia4=ist(4),ien(4)
               ia(4) = ia4
             do ia3=ist(3),ien(3)
                ia(3) = ia3
              do ia2=ist(2),ien(2)
                 ia(2) = ia2
               do ia1=ist(1),ien(1)
                  ia(1) = ia1
                  i = i + 1
                  do k=1,2
                    call get_tr_sum_address
     &                   (ntf,k,ndim,lb,ib,ia,n_sum_b,kbase,ig)
                      tr_sum(ig) = tott(i,k)
                  enddo
               enddo
              enddo
             enddo
            enddo
           enddo
          enddo
         enddo
        enddo
       enddo
      enddo

      return
      end

************************************************************************
*                                                                      *
      subroutine trEVENT_add_tr_sum(
     &      ndim, nstep, istart, ibase, n_sum_l,n_sum,
     &      ndim0,nstep0,istart0,ibase0,n0_l,   n0,
     &                              tr, tr_sum, maxcas)
*
* ntf          : file number in nfiles
* ndim         : counts of parameters
* nstep(ndim)  : do step number
* istart(ndim) : start of do loop
* ibase(ndime) : base of do loop
* n_sum_l      : lower of each parameter
* n_sum        : max of each parameter
* maxtott      : tott(maxtott,2)
* tott         : set data
* tr(n_sum(1),-,n_sum(ndim))       : database
* tr_sum(n_sum(1),-,n_sum(ndim),2) : database
* maxcas       : occurrence number
*                                                                     *
************************************************************************

      implicit double precision (a-h,o-z)

      integer :: ndim, nstep(*), istart(*), ibase(*),
     &           n_sum_l(*),n_sum(*)
      integer :: ndim0,nstep0(*),istart0(*),ibase0(*),
     &           n0_l(*),n0(*)
      real(8) :: tr(*)
      real(8) :: tr_sum(*)
      integer :: maxcas

      integer :: ntf,ntf0
      integer :: n_sum_b(11), ist(11), ien(11), iloop(11), kbase,
     &           lb(11),ib(11), ia(11)
      integer :: n0_b(11), ist0(11), ien0(11), iloop0(11), kbase0,
     &           lb0(11),ib0(11), ia0(11)

      real(8),allocatable :: tr_add(:)

      ntf0 = 1
      ntf = 1

      kbase = 0
      n_sum_b(1) = 1
      do i=2,ndim+1
        n_sum_b(i) = n_sum_b(i-1)*(n_sum(i-1) + 1- n_sum_l(i-1))
      enddo
      kbase = n_sum_b(ndim+1)
      n_sum_b(1) = 0

      allocate (tr_add(kbase))
      tr_add(:) = 0.0d0

      lb(:) = 1
      ib(:) = 1
      ist(:) = 1
      do i=1,ndim
        lb(i) = n_sum_l(i)
        ib(i) = ibase(i)
        ist(i) = istart(i)
      enddo
      ien(:) = ist(:)
      do i=1,ndim
        ien(i) = nstep(i)
      enddo

      kbase0 = 0
      n0_b(1) = 1
      do i=2,ndim0+1
        n0_b(i) = n0_b(i-1)*(n0(i-1) + 1 - n0_l(i-1))
      enddo
      kbase0 = n0_b(ndim0+1)
      n0_b(1) = 0

      lb0(:) = 1
      ib0(:) = 1
      ist0(:) = 1
      do i=1,ndim0
        lb0(i) = n0_l(i)
        ib0(i) = ibase0(i)
        ist0(i) = istart0(i)
      enddo
      ien0(:) = ist0(:)
      do i=1,ndim0
        ien0(i) = nstep0(i)
      enddo

        do ia10=ist0(10),ien0(1)
          ia0(10) = ia10
          ia(10)  = ia10
          if(ist(10) == ien(10)) then
            ia(10)=ist(10)
          endif
          do ia9=ist0(9),ien0(9)
            ia0(9) = ia9
            ia(9)  = ia9
            if(ist(9) == ien(9)) then
              ia(9)=ist(9)
            endif
            do ia8=ist0(8),ien0(8)
              ia0(8) = ia8
              ia(8)  = ia8
              if(ist(8) == ien(8)) then
                ia(8)=ist(8)
              endif
              do ia7=ist0(7),ien0(7)
                ia0(7) = ia7
                ia(7)  = ia7
                if(ist(7) == ien(7)) then
                  ia(7)=ist(7)
                endif
                do ia6=ist0(6),ien0(6)
                  ia0(6) = ia6
                  ia(6)  = ia6
                  if(ist(6) == ien(6)) then
                    ia(6)=ist(6)
                  endif
                  do ia5=ist0(5),ien0(5)
                    ia0(5) = ia5
                    ia(5)  = ia5
                    if(ist(5) == ien(5)) then
                      ia(5)=ist(5)
                    endif
                    do ia4=ist0(4),ien0(4)
                      ia0(4) = ia4
                      ia(4)  = ia4
                      if(ist(4) == ien(4)) then
                        ia(4)=ist(4)
                      endif
                      do ia3=ist0(3),ien0(3)
                        ia0(3) = ia3
                        ia(3)  = ia3
                        if(ist(3) == ien(3)) then
                          ia(3)=ist(3)
                        endif
                        do ia2=ist0(2),ien0(2)
                          ia0(2) = ia2
                          ia(2)  = ia2
                          if(ist(2) == ien(2)) then
                            ia(2)=ist(2)
                          endif
                          do ia1=ist0(1),ien0(1)
                            ia0(1) = ia1
                            ia(1)  = ia1
                            if(ist(1) == ien(1)) then
                              ia(1)=ist(1)
                            endif
                              k = 1
                              call get_tr_sum_address
     &                        (ntf0,k,ndim0,lb0,ib0,ia0,n0_b,kbase0,ig0)
                              call get_tr_sum_address
     &                        (ntf,k,ndim,lb,ib,ia,n_sum_b,kbase,ig)
                              tr_add(ig) = tr_add(ig) + tr(ig0)/maxcas
                          enddo
                        enddo
                      enddo
                    enddo
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo

        do ia10=ist(10),ien(1)
          ia(10)  = ia10
          do ia9=ist(9),ien(9)
            ia(9)  = ia9
            do ia8=ist(8),ien(8)
              ia(8)  = ia8
              do ia7=ist(7),ien(7)
                ia(7)  = ia7
                do ia6=ist(6),ien(6)
                  ia(6)  = ia6
                  do ia5=ist(5),ien(5)
                    ia(5)  = ia5
                    do ia4=ist(4),ien(4)
                      ia(4)  = ia4
                      do ia3=ist(3),ien(3)
                        ia(3)  = ia3
                        do ia2=ist(2),ien(2)
                          ia(2)  = ia2
                          do ia1=ist(1),ien(1)
                            ia(1)  = ia1
                              k = 1
                              call get_tr_sum_address
     &                         (ntf,k,ndim,lb,ib,ia,n_sum_b,kbase,ig)
                              tr_sum(ig)=tr_sum(ig)+tr_add(ig)
                              k = 2
                              call get_tr_sum_address
     &                         (ntf,k,ndim,lb,ib,ia,n_sum_b,kbase,ig2)
                              tr_sum(ig2)=tr_sum(ig2)+ tr_add(ig)**2
                          enddo
                        enddo
                      enddo
                    enddo
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
!
      deallocate(tr_add) 

      return
      end

************************************************************************
*                                                                      *
      subroutine get_tr_sum_address
     &          (ntf, k, ndim, lb,  ib, ia, n_sum_b, kbase,ig)
*
* ig : address of data
*                                                                      *
************************************************************************

      implicit double precision (a-h,o-z)

      integer :: ndim

      integer :: n_sum_b(11),lb(11),ib(11),ia(11)

      ig = ia(1) + ib(1) - lb(1) + kbase * (k-1) + kbase * 2 * (ntf-1)
      do i=1,ndim
        ig = ig + (ia(i) + ib(i) - lb(i) -1) * n_sum_b(i)
      enddo

      return
      end

************************************************************************
*                                                                      *
      subroutine set_integer_2(ndata,n1,n2)
*                                                                      *
************************************************************************

      implicit double precision (a-h,o-z)
      
      dimension :: ndata(*)

      ndata(1) = n1
      ndata(2) = n2

      return
      end
************************************************************************
*                                                                      *
      subroutine set_integer_3(ndata,n1,n2,n3)
*                                                                      *
************************************************************************

      implicit double precision (a-h,o-z)
      
      dimension :: ndata(*)

      ndata(1) = n1
      ndata(2) = n2
      ndata(3) = n3

      return
      end
************************************************************************
*                                                                      *
      subroutine set_integer_4(ndata,n1,n2,n3,n4)
*                                                                      *
************************************************************************

      implicit double precision (a-h,o-z)
      
      dimension :: ndata(*)

      ndata(1) = n1
      ndata(2) = n2
      ndata(3) = n3
      ndata(4) = n4

      return
      end
************************************************************************
*                                                                      *
      subroutine set_integer_5(ndata,n1,n2,n3,n4,n5)
*                                                                      *
************************************************************************

      implicit double precision (a-h,o-z)
      
      dimension :: ndata(*)

      ndata(1) = n1
      ndata(2) = n2
      ndata(3) = n3
      ndata(4) = n4
      ndata(5) = n5

      return
      end
************************************************************************
*                                                                      *
      subroutine set_integer_6(ndata,n1,n2,n3,n4,n5,n6)
*                                                                      *
************************************************************************

      implicit double precision (a-h,o-z)
      
      dimension :: ndata(*)

      ndata(1) = n1
      ndata(2) = n2
      ndata(3) = n3
      ndata(4) = n4
      ndata(5) = n5
      ndata(6) = n6
      return
      end
************************************************************************
*                                                                      *
      subroutine set_integer_7(ndata,n1,n2,n3,n4,n5,n6,n7)
*                                                                      *
************************************************************************

      implicit double precision (a-h,o-z)
      
      dimension :: ndata(*)

      ndata(1) = n1
      ndata(2) = n2
      ndata(3) = n3
      ndata(4) = n4
      ndata(5) = n5
      ndata(6) = n6
      ndata(7) = n7

      return
      end
************************************************************************
*                                                                      *
      subroutine set_integer_8(ndata,n1,n2,n3,n4,n5,n6,n7,n8)
*                                                                      *
************************************************************************

      implicit double precision (a-h,o-z)
      
      dimension :: ndata(*)

      ndata(1) = n1
      ndata(2) = n2
      ndata(3) = n3
      ndata(4) = n4
      ndata(5) = n5
      ndata(6) = n6
      ndata(7) = n7
      ndata(8) = n8

      return
      end

************************************************************************
*                                                                      *
      subroutine set_integer_9(ndata,n1,n2,n3,n4,n5,n6,n7,n8,n9)
*                                                                      *
************************************************************************

      implicit double precision (a-h,o-z)
      
      dimension :: ndata(*)

      ndata(1) = n1
      ndata(2) = n2
      ndata(3) = n3
      ndata(4) = n4
      ndata(5) = n5
      ndata(6) = n6
      ndata(7) = n7
      ndata(8) = n8
      ndata(9) = n9

      return
      end

************************************************************************
*                                                                      *
      subroutine set_integer_10(ndata,n1,n2,n3,n4,n5,n6,n7,n8,n9,n10)
*                                                                      *
************************************************************************

      implicit double precision (a-h,o-z)
      
      dimension :: ndata(*)

      ndata(1)  = n1
      ndata(2)  = n2
      ndata(3)  = n3
      ndata(4)  = n4
      ndata(5)  = n5
      ndata(6)  = n6
      ndata(7)  = n7
      ndata(8)  = n8
      ndata(9)  = n9
      ndata(10) = n10

      return
      end

