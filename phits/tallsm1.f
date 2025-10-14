************************************************************************
*                                                                      *
      subroutine settal(io,jo,ierr)
*                                                                      *
*       set tally input                                                *
*       modified by S. Hashimoto and K.Niita on 2011/06/20             *
*                                                                      *
************************************************************************
      use partmod, only: itpan, itpat, jtpat ! frtati 2021/10/05
      use moddas
      use moddas_region
      use moddas_tally
      use talmod, only: deist ! S.H. 2022.12.16
      use tetramod, only: nelemtot

      implicit real*8 (a-h,o-z)

      include 'param.inc'
      include 'param01.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      common /mpi00/  npe, me
      common /tcntl/  icntl, inucr
      common /talmm/  nmmax, lmmax, itlmx
      common /dmpfil/ idmpf
      common /eparm/  esmax, esmin, emin(20)

      common /talout/ itall
      common /talsav/ iptall

      common /tall00/ itnm, ital(itlmax), itals(itlmax), italm(itlmax)
      common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                itout(itlmax), ittwo(itlmax)

      common /tall02/ itrgn(itlmax), itrgm(itlmax), itreg(itlmax)
      common /tall24/ itrcm(itlmax), itrcg(itlmax)

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
      common /udtall1/ iudtfll(itlmax,50), cudtfln(itlmax,50)
      character cudtfln*100

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
      common /tall31/ iterl(itlmax)

      common /tall32/ itelc(itlmax)

      common /tall33/ itln(itlmax,2), itli(itlmax,2), itlr(itlmax,2),
     &                rtdm(itlmax,2)
      common /tall34/ itvm(itlmax)
      common /tall35/ iteps(itlmax)

      common /tall40/ rtorg(itlmax,3), rteye(itlmax,3), rtlit(itlmax,3),
     &                rtwin(itlmax,3), rtbrt(itlmax,2), ithvn(itlmax),
     &                itwin(itlmax,2), itbox(itlmax), itmir(itlmax),
     &                rtbox(itlmax,5,10), rtout(itlmax), rthet(itlmax),
     &                itlin(itlmax), itshd(itlmax), itgxs(itlmax)

      common /tall42/ itmbn(itlmax), itmbt(itlmax)
      common /tall43/ itmgn(itlmax), itmgm(itlmax), itmeg(itlmax)
      common /tall44/ itmcm(itlmax), itmcg(itlmax)

      common /tall45/ itmlp(itlmax), itmln(itlmax,6), itmst(itlmax),
     &                itmli(itlmax,6), rtmme(itlmax,6), itmnt(itlmax,6),
     &                itmpn(itlmax,6), itmpt(itlmax,6,6,2)

      common /tall46/ itmtr(itlmax,4), rtmtr(itlmax,13)
      common /tall47/ itbtr(itlmax,5,4), rtbtr(itlmax,5,13)
      common /tall48/ itmdp(itlmax,0:30), itmdf(itlmax)

      common /tall49/ itglt(itlmax)
      common /tall50/ itlmt(itlmax), ite2l(itlmax)   ! CCSE 2022.08.31
      common /tall51/ eletb(2000), nlete

      common /tall55/ itlmt2(itlmax)

      common /tall56/ itety2(itlmax), itenm2(itlmax), iterg2(itlmax),
     &                rtemi2(itlmax), rtema2(itlmax), rtedl2(itlmax)
      common /tall57/ itsun(itlmax),rtdim(itlmax),rtucv(itlmax),
     &rtrho(itlmax),itmodel(itlmax) ! T.Sato 2022/08/14

      common /inggs/  iog, igcel, ioa, igsuf, iob, igtrs
      common /celdb/  idsn(kvlmax), idtn(kvlmax)

      common /cusrtally/ iusrtally, iudtf(50)


      common /gsline/ nowgshow, igsline

      dimension idas(1)
      equivalence ( das, idas )

*-----------------------------------------------------------------------

      common /regcrs/ icrsflx

*-----------------------------------------------------------------------

      common /inpec/ ititl, ipara, ibody, iregn, llarr, itby, itar
      common /regdc/ idrg(kvlmax), idgr(kvmmax)
      common /kmat1d/ idmn(0:kvlmax), idnm(kvmmax)
      common /ndemax/ dnmax(20)

      common /pointt/ itpont

*-----------------------------------------------------------------------

C S.H. added for Dump-Restart on 2014/5/7
      common /stat / istdev, irestart, ireschk

*-----------------------------------------------------------------------
      common /cparm/  maxbch, maxcas
      common /tall78/ itism(itlmax), itist(10,itlmax), itjst(10,itlmax),
     &                itkst(10,itlmax), itstt(itlmax), itsdd(itlmax)

*-----------------------------------------------------------------------

      character chme*5, filnm*100

*-----------------------------------------------------------------------

      integer,parameter:: mnmax = 0

*-----------------------------------------------------------------------

      common /tall67/ itstd(itlmax), rtstd(itlmax)
      data itstd / itlmax*1 /

*-----------------------------------------------------------------------

      common /tall91/ itextstat(itlmax), mftal(itlmax) ! S.H. extstat 2024.3.18
      common /tall93/ tprodenmn(itlmax), tprodenmx(itlmax),
     &     nbtproden(itlmax), itprodenchk(itlmax) !S.H. extstat 2024.4.28

*-----------------------------------------------------------------------

      integer,allocatable :: ielem2ir(:)

*-----------------------------------------------------------------------

c           if( icntl .eq. 1 ) return
            if( icntl .eq. 2 ) return
            if( icntl .eq. 4 ) return

            ierr = 0

*-----------------------------------------------------------------------
*     flag for point tally
*-----------------------------------------------------------------------
               itpont = 0

*-----------------------------------------------------------------------
*     set tally
*-----------------------------------------------------------------------

            ntstar = mmmax
            nmmax  = mmmax
            lmmax  = mmmax
            itlmx  = 0

*-----------------------------------------------------------------------
*           mesh for LET energy bin
*           esmin to msmax by nlete is increase by T.Sato
*           esmin to esmax by nlete (=400, 1900)
*-----------------------------------------------------------------------

                  nlete = 1900
                  emins = min( esmin, 1.d-6)
                  emaxs = max( esmax, 3.d+5)

                  xlinc = log( emaxs / emins ) /  nlete

               do i = 1, nlete + 1

                  eletb(i) = exp( ( i - 1.0d0 ) * xlinc ) * emins

               end do

*-----------------------------------------------------------------------

      if( itnm .gt. 0 ) then

*-----------------------------------------------------------------------
*        check file name in tally
*-----------------------------------------------------------------------

            do m = 1, itnm
              do k = 1, itfln(m)
                 do l = m + 1, itnm
                    do i = 1, itfln(l)

                      if( ital(m).ne.20 .and. ital(l).ne.20 ) then

                       do j = 1, itfll(l,i)
                          if( ctfln(l,i)(j:j) .ne. ctfln(m,k)(j:j) )
     &                        goto 88
                       end do

                      else

                       if( ital(m).eq.20 ) then

                        do j = 1, itfll(l,i)
                          if( ctfln(l,i)(j:j) .ne. cudtfln(m,k)(j:j) )
     &                        goto 88
                        end do

                       else if( ital(l).eq.20 ) then

                        do j = 1, iudtfll(l,i)
                          if( cudtfln(l,i)(j:j) .ne. ctfln(m,k)(j:j) )
     &                        goto 88
                        end do

                       else

                        do j = 1, iudtfll(l,i)
                          if( cudtfln(l,i)(j:j) .ne. cudtfln(m,k)(j:j) )
     &                        goto 88
                        end do

                       end if

                      end if

                        ierr = 1
                        write(io,'(''*Error: in '',i2,
     &                  ''-th tally, file name is the same as that of ''
     &                  ,i2,''-th tally'')') m, l
                      if( ital(m).ne.20 .and. ital(l).ne.20 ) then

                        write(io,'(''  file = '',100a1)')
     &                        (ctfln(l,i)(j:j),j=1,itfll(l,i))

                      else

                       if( ital(m).eq.20 ) then

                        write(io,'(''  file = '',100a1)')
     &                        (ctfln(l,i)(j:j),j=1,itfll(l,i))

                       else if( ital(l).eq.20 ) then

                        write(io,'(''  file = '',100a1)')
     &                        (cudtfln(l,i)(j:j),j=1,iudtfll(l,i))

                       else

                        write(io,'(''  file = '',100a1)')
     &                        (cudtfln(l,i)(j:j),j=1,iudtfll(l,i))

                       end if

                      end if

                        ErrCha = ''
                        MsgID = 'L:309/R:settal/F:tallsm1.f'
                        call ErrWrite(MsgID, ErrCha)
                        write(jo,'(''*Error: in '',i2,
     &                  ''-th tally, file name is the same as that of ''
     &                  ,i2,''-th tally'')') m, l
                      if( ital(m).ne.20 .and. ital(l).ne.20 ) then

                        write(jo,'(''  file = '',100a1)')
     &                        (ctfln(l,i)(j:j),j=1,itfll(l,i))

                      else

                       if( ital(m).eq.20 ) then

                        write(jo,'(''  file = '',100a1)')
     &                        (ctfln(l,i)(j:j),j=1,itfll(l,i))

                       else if( ital(l).eq.20 ) then

                        write(jo,'(''  file = '',100a1)')
     &                        (cudtfln(l,i)(j:j),j=1,iudtfll(l,i))

                       else

                        write(jo,'(''  file = '',100a1)')
     &                        (cudtfln(l,i)(j:j),j=1,iudtfll(l,i))

                       end if

                      end if

                        goto 999

   88                  continue

                    end do
                 end do
              end do
            end do

*-----------------------------------------------------------------------
*        check transform id in tally and in box
*-----------------------------------------------------------------------

            do m = 1, itnm

*-----------------------------------------------------------------------

               if( itmtr(m,1) .gt. 0 ) then

                     l = 0

                  do k = 1, igtrs

                     if( itmtr(m,3) .eq. idtn(k) ) l = k

                  end do

                  if( l .gt. 0 ) then

                     itmtr(m,4) = l

                  else

                     ierr = 1
                     write(io,'(''*Error: in '',i2,
     &               ''-th tally, trcl ='',i4,'' is not defined '',
     &               ''in [transform]'')') m, itmtr(m,3)

                     ErrCha = ''
                     MsgID = 'L:379/R:settal/F:tallsm1.f'
                     call ErrWrite(MsgID, ErrCha)
                     write(jo,'(''*Error: in '',i2,
     &               ''-th tally, trcl ='',i4,'' is not defined '',
     &               ''in [transform]'')') m, itmtr(m,3)

                     goto 999

                  end if

               end if

*-----------------------------------------------------------------------

               if( itbox(m) .gt. 0 ) then

                  do i = 1, itbox(m)

                     if( itbtr(m,i,1) .gt. 0 ) then

                           l = 0

                        do k = 1, igtrs

                           if( itbtr(m,i,3) .eq. idtn(k) ) l = k

                        end do

                        if( l .gt. 0 ) then

                           itbtr(m,i,4) = l

                        else

                           ierr = 1
                           write(io,'(''*Error: in '',i2,
     &                     ''-th tally, trcl ='',i4,'' is not '',
     &                     ''definend in [transform]'')') m, itmtr(m,3)

                           ErrCha = ''
                           MsgID = 'L:419/R:settal/F:tallsm1.f'
                           call ErrWrite(MsgID, ErrCha)
                           write(jo,'(''*Error: in '',i2,
     &                     ''-th tally, trcl ='',i4,'' is not '',
     &                     ''definend in [transform]'')') m, itmtr(m,3)

                           goto 999

                        end if

                     end if

                  end do

               end if

*-----------------------------------------------------------------------

            end do

*-----------------------------------------------------------------------
*     check rshow with mesh=tet
*-----------------------------------------------------------------------

            do m = 1, itnm

*-----------------------------------------------------------------------

             if( itrsh(m) .gt. 0 )then ! rshow > 0

              if( itmsh(m) .eq. 4 )then ! mesh = tet

               if( itglt(m) .ge. 2 ! gslat >= 2
     &              .or. (itglt(m) .eq. -1 .and. igsline .ge. 2) )then

                write(*,'(3a,i2,a)')'*** Warning: ',
     &               'gslat=1 is recommended for rshow>0 ',
     &               'with mesh=tet in ', m,'-th tally'
                write(*,'(2a)')
     &               '    Otherwise, results for each tetrahedron ',
     &               'are not shown'

               endif

              endif

             endif

*-----------------------------------------------------------------------

            end do

*-----------------------------------------------------------------------
*     set memory and the others
*-----------------------------------------------------------------------

            lmmax = 0
            nmmax = 0

            immax = 0
            jmmax = 0

      itrcg(1:itlmax) = itreg(1:itlmax)
      if( allocated(idas_itreg) ) then
         call moddas_allocate_int(itreg(itlmax), idas_itrcg)
         idas_itrcg(1:itrcg(itlmax)) = idas_itreg(1:itreg(itlmax))
         deallocate( idas_itreg )
      end if

      call moddas_allocate_int(MAX_NUM_ITREG, idas_itreg_temporary)

      itrcc(1:itlmax) = itrcr(1:itlmax)
      if( allocated(idas_itrcr) ) then
         call moddas_allocate_int(itrcr(itlmax), idas_itrcc)
         idas_itrcc(1:itrcc(itlmax)) = idas_itrcr(1:itrcr(itlmax))
      end if

      itmcg(1:itlmax) = itmeg(1:itlmax)
      if( allocated(idas_itmeg) ) then
         call moddas_allocate_int(itmeg(itlmax), idas_itmcg)
         idas_itmcg(1:itmcg(itlmax)) = idas_itmeg(1:itmeg(itlmax))
         deallocate( idas_itmeg )
      end if
      do m = 1, itnm


*-----------------------------------------------------------------------
*        t-track tally
*-----------------------------------------------------------------------

         if( ital(m) .eq. 1 ) then

*-----------------------------------------------------------------------
*              check region
*-----------------------------------------------------------------------

            if( itmsh(m) .eq. 1 ) then

                  igm  = 1

                  itrcm(m) = itrgm(m)

                  ntrn = itrgn(m)
                  mtrn = itrgm(m)

                  ign  = itrcg(m)
                  call tregion2(io,jo,ierr,ntrn,mtrn,idas_itrcg(ign)
     &                          ,igm,MAX_NUM_ITREG,idas_itreg_temporary)

                  if( mtrn > MAX_NUM_ITREG ) then
                     write(ErrCha,'(a,a,i5,a,a,i5,a)')
     &                    'sub.settal t-track@tallsm1.f'
     &                       //' ?dimension over idas_itreg_temporary?'
     &                       //' mtrn > MAX_NUM_ITREG'
     &                    ,' (mtrn=',mtrn,')'
     &                    ,' (MAX_NUM_ITREG@moddas.f=',MAX_NUM_ITREG,')'
                     ErrID = 'L:535/R:settal/F:tallsm1.f'
                     call ErrWrite(ErrID,ErrCha)
                  endif

                     if( ierr .ne. 0 ) then

                        write(io,'('' **** above error in '',i2,
     &                  ''-th tally'')') m

                        ErrCha = ''
                        MsgID = 'L:545/R:settal/F:tallsm1.f'
                        call ErrWrite(MsgID, ErrCha)
                        write(jo,'('' **** above error in '',i2,
     &                  ''-th tally'')') m

                        goto 999

                     end if

                  call moddas_reallocate_int(
     &                    itlmax, m, mtrn, itreg, idas_itreg)
                  idas_itreg(itreg(m):itreg(m)+mtrn-1)
     &                      = idas_itreg_temporary(1:mtrn)
                  itrgn(m) = ntrn
                  itrgm(m) = mtrn

                  itsmn(m) = itsmn(m) + ( mtrn + mod(mtrn,2) ) / 2


            else if( itmsh(m) .eq. 4 ) then

                  itrcm(m) = itrgm(m)

                  ntrn = itrgn(m)
                  mtrn = itrgm(m)

                  ign  = itrcg(m)

                  call moddas_allocate_int(nelemtot,ielem2ir)

                  call ttetmesh2(io,jo,ierr,ntrn,mtrn,idas_itrcg(ign),
     &                 MAX_NUM_ITREG,idas_itreg_temporary,ielem2ir)

                  if( mtrn > MAX_NUM_ITREG ) then
                     write(ErrCha,'(a,a,i5,a,a,i5,a)')
     &                    'sub.settal t-track@tallsm1.f'
     &                       //' ?dimension over idas_itreg_temporary?'
     &                       //' mtrn > MAX_NUM_ITREG'
     &                    ,' (mtrn=',mtrn,')'
     &                    ,' (MAX_NUM_ITREG@moddas.f=',MAX_NUM_ITREG,')'
                     ErrID = 'L:585/R:settal/F:tallsm1.f'
                     call ErrWrite(ErrID,ErrCha)
                  endif

                     if( ierr .ne. 0 ) then

                        write(io,'('' **** above error in '',i2,
     &                  ''-th tally'')') m

                        ErrCha = ''
                        MsgID = 'L:595/R:settal/F:tallsm1.f'
                        call ErrWrite(MsgID, ErrCha)
                        write(jo,'('' **** above error in '',i2,
     &                  ''-th tally'')') m

                        goto 999

                     end if

                  if(ntrn.gt.0)then

                   call moddas_reallocate_int(
     &                  itlmax, m, mtrn+ntrn,
     &                  itreg, idas_itreg)
                   idas_itreg(itreg(m):itreg(m)+mtrn-1)
     &                  = idas_itreg_temporary(1:mtrn)
                   idas_itreg(itreg(m)+mtrn:itreg(m)+mtrn+ntrn-1)
     &                  = ielem2ir(1:ntrn)

                   itrgn(m) = ntrn
                   itrgm(m) = mtrn+ntrn

                  else

                   ntrn0=idas_itreg_temporary(4)
                   idas_itreg_temporary(4)=0

                   call moddas_reallocate_int(
     &                  itlmax, m, mtrn, itreg, idas_itreg)
                   idas_itreg(itreg(m):itreg(m)+mtrn-1)
     &                  = idas_itreg_temporary(1:mtrn)
                   itrgn(m) = ntrn0
                   itrgm(m) = mtrn

                  endif

                  call moddas_deallocate_int(ielem2ir)

            end if

*-----------------------------------------------------------------------
*              check material
*-----------------------------------------------------------------------

               if( itmtn(m) .gt. 0 ) then

                  do i = 1, itmtn(m)

                        imat = ismte( itmtt(m) + i - 1 )

                     if( imat .gt. 0 ) then
                     if( idnm( imat ) .eq. 0 ) then

                        write(io,'('' **** Error in '',i2,
     &                  ''-th tally : material ['',i4,
     &                  ''] is not defined '',
     &                  ''in material section.'')') m, imat

                        ErrCha = ''
                        MsgID = 'L:654/R:settal/F:tallsm1.f'
                        call ErrWrite(MsgID, ErrCha)
                        write(jo,'('' **** Error in '',i2,
     &                  ''-th tally : material ['',i4,
     &                  ''] is not defined '',
     &                  ''in material section.'')') m, imat

                        goto 999

                     end if
                     end if

                  end do

               end if

*-----------------------------------------------------------------------
*           multiplier
*-----------------------------------------------------------------------

            if( itmst(m) .eq. 0 ) itmst(m) = 1

*-----------------------------------------------------------------------
*           specify the subroutine number and storage position
*-----------------------------------------------------------------------

               italm(m) = mmmax

            if( itmsh(m) .eq. 1 ) then

               itals(m) = 1

               if( itrsh(m) .ne. 0 .and. itrgn(m) .le. 1 ) goto 940

            else if( itmsh(m) .eq. 2 ) then

               itals(m) = 2

            else if( itmsh(m) .eq. 3 ) then

               itals(m) = 3

            else if( itmsh(m) .eq. 4 ) then
               itals(m) = 49

            end if

*-----------------------------------------------------------------------
*           stop OpenMP or MPI execution without setting proden prameters
*           (prodenmn, prodenmx, nbproden)
*-----------------------------------------------------------------------

!$          if( mftal(m).gt.5 .and. itprodenchk(m).eq.0 ) then
!$             write(ErrCha,'("Error: OpenMP execution cannot ",
!$   &           "be used in test calculation to estimate scale of ",
!$   &           "probability density function.")')
!$             ErrID = 'L:800/R:settal/F:tallsm1.f'
!$             call ErrWrite(ErrID,ErrCha)
!$
!$             goto 999
!$
!$          endif

            if( npe.gt.1 .and. mftal(m).gt.5 .and.
     &         itprodenchk(m).eq.0 ) then
               write(ErrCha,'("Error: MPI execution cannot ",
     &           "be used in test calculation to estimate scale of ",
     &           "probability density function.")')
               ErrID = 'L:722/R:settal/F:tallsm1.f'
               call ErrWrite(ErrID,ErrCha)

               goto 999

            endif

*-----------------------------------------------------------------------
*        t-adjoint tally
*-----------------------------------------------------------------------

         else if( ital(m) .eq. 19 ) then

*-----------------------------------------------------------------------
*              check region
*-----------------------------------------------------------------------

            if( itmsh(m) .eq. 1 ) then

                  igm  = 1

                  itrcm(m) = itrgm(m)

                  ntrn = itrgn(m)
                  mtrn = itrgm(m)
                  ign  = itrcg(m)
                  call tregion2(io,jo,ierr,ntrn,mtrn,idas_itrcg(ign)
     &                          ,igm,MAX_NUM_ITREG,idas_itreg_temporary)

                  if( mtrn > MAX_NUM_ITREG ) then
                     write(ErrCha,'(a,a,i5,a,a,i5,a)')
     &                    'sub.settal t-adjoint@tallsm1.f'
     &                       //' ?dimension over idas_itreg_temporary?'
     &                       //' mtrn > MAX_NUM_ITREG'
     &                    ,' (mtrn=',mtrn,')'
     &                    ,' (MAX_NUM_ITREG@moddas.f=',MAX_NUM_ITREG,')'
                     ErrID = 'L:758/R:settal/F:tallsm1.f'
                     call ErrWrite(ErrID,ErrCha)
                  endif

                     if( ierr .ne. 0 ) then

                        write(io,'('' **** above error in '',i2,
     &                  ''-th tally'')') m

                        ErrCha = ''
                        MsgID = 'L:768/R:settal/F:tallsm1.f'
                        call ErrWrite(MsgID, ErrCha)
                        write(jo,'('' **** above error in '',i2,
     &                  ''-th tally'')') m

                        goto 999

                     end if

                  call moddas_reallocate_int(
     &                    itlmax, m, mtrn, itreg, idas_itreg)
                  idas_itreg(itreg(m):itreg(m)+mtrn-1)
     &                      = idas_itreg_temporary(1:mtrn)
                  itrgn(m) = ntrn
                  itrgm(m) = mtrn

                  itsmn(m) = itsmn(m) + ( mtrn + mod(mtrn,2) ) / 2


            end if

*-----------------------------------------------------------------------
*              check material
*-----------------------------------------------------------------------

               if( itmtn(m) .gt. 0 ) then

                  do i = 1, itmtn(m)

                        imat = ismte( itmtt(m) + i - 1 )

                     if( imat .gt. 0 ) then
                     if( idnm( imat ) .eq. 0 ) then

                        write(io,'('' **** Error in '',i2,
     &                  ''-th tally : material ['',i4,
     &                  ''] is not defined '',
     &                  ''in material section.'')') m, imat

                        ErrCha = ''
                        MsgID = 'L:808/R:settal/F:tallsm1.f'
                        call ErrWrite(MsgID, ErrCha)
                        write(jo,'('' **** Error in '',i2,
     &                  ''-th tally : material ['',i4,
     &                  ''] is not defined '',
     &                  ''in material section.'')') m, imat

                        goto 999

                     end if
                     end if

                  end do

               end if

*-----------------------------------------------------------------------
*           multiplier
*-----------------------------------------------------------------------

            if( itmst(m) .eq. 0 ) itmst(m) = 1

*-----------------------------------------------------------------------
*           specify the subroutine number and storage position
*-----------------------------------------------------------------------

               italm(m) = mmmax

            if( itmsh(m) .eq. 1 ) then

               itals(m) = 43

               if( itrsh(m) .ne. 0 .and. itrgn(m) .le. 1 ) goto 940

            else if( itmsh(m) .eq. 2 ) then

               itals(m) = 44

            else if( itmsh(m) .eq. 3 ) then

               itals(m) = 45

            end if

*-----------------------------------------------------------------------
*        t-cross tally
*-----------------------------------------------------------------------

         else if( ital(m) .eq. 2 ) then

*-----------------------------------------------------------------------
*           open file for dump data
*-----------------------------------------------------------------------

            if( itmdp(m,0) .ne. 0 ) then

                  idmpf = idmpf + 1

                  if( idmpf .gt. 20 ) then

                        write(io,'('' **** Error in '',i2,
     &                  ''-th tally : number of dump file''
     &                  '' exceeds 20 files'')') m

                        ErrCha = ''
                        MsgID = 'L:873/R:settal/F:tallsm1.f'
                        call ErrWrite(MsgID, ErrCha)
                        write(jo,'('' **** Error in '',i2,
     &                  ''-th tally : number of dump file''
     &                  '' exceeds 20 files'')') m

                        goto 999

                  end if

                  itmdf(m) = idmpf + 80

C S.H. revised for Dump-Restart on 2014/5/7
                  filnm = ctfln(m,1)

                  do i = itfll(m,1), 1, -1
                     if ( ctfln(m,1)(i:i) .eq. '.' ) exit
                  end do
                  if ( i .eq. 0 ) i = itfll(m,1)+1

                  filnm(1:i-1) = ctfln(m,1)(1:i-1)
                  filnm(i:i+3) = '_dmp'
                  if ( i .lt. itfll(m,1) ) then
                     filnm(i+4:itfll(m,1)+4) = ctfln(m,1)(i:itfll(m,1))
                  end if

               if( npe .le. 1 ) then

                  if( itmdp(m,0) .gt. 0 ) then

                     if ( irestart .eq. 1 ) then

                        open(itmdf(m), file = filnm,
     &                       form='unformatted', status = 'unknown',
     &                       access = 'append' )

                     else

                        open(itmdf(m), file = filnm,
     &                       form='unformatted', status = 'unknown' )

                     end if

                  else

                     if ( irestart .eq. 1 ) then

                        open(itmdf(m), file = filnm,
     &                       form='formatted', status = 'unknown',
     &                       access = 'append' )
                     else

                        open(itmdf(m), file = filnm,
     &                       form='formatted', status = 'unknown' )

                     end if

                  end if

               else if( me .gt. 0 ) then

                     iorder = aint(log10(real(npe)))+1
                     if ( iorder .lt. 3) iorder = 3

                     write(chme,'(i5.5)') me

                     lengfilnm = len_trim(filnm)
                     filnm(lengfilnm+1:lengfilnm+1+iorder)
     &                    = '.' // chme(6-iorder:5)

                  if( itmdp(m,0) .gt. 0 ) then

                     if ( irestart .eq. 1 ) then

                        open(itmdf(m), file = filnm,
     &                       form='unformatted', status = 'unknown',
     &                       access = 'append' )

                     else

                        open(itmdf(m), file = filnm,
     &                       form='unformatted', status = 'unknown' )

                     end if

                  else

                     if ( irestart .eq. 1 ) then

                        open(itmdf(m), file = filnm,
     &                       form='formatted', status = 'unknown',
     &                       access = 'append' )

                     else

                        open(itmdf(m), file = filnm,
     &                       form='formatted', status = 'unknown' )

                     end if

                  end if

               end if

               if(npe.le.1.or.me.gt.0)then           !FURUTA20150427
                call chkdmpomp(io,m,itmdp(m,0),ierr) !FURUTA20150427
                if(ierr.ne.0)goto 999                !FURUTA20150427
               endif                                 !FURUTA20150427

            end if

*-----------------------------------------------------------------------
*              cg/gg flag for angle information
*-----------------------------------------------------------------------

               if( itmsh(m) .eq. 1 .and.
     &           ( itout(m) .eq. 1 .or. itout(m) .ge. 8 )) then  ! T.Sato 2021/05/05

                  icrsflx = 1

               end if

*-----------------------------------------------------------------------
*              check region
*-----------------------------------------------------------------------

            if( itmsh(m) .eq. 1 ) then

*-----------------------------------------------------------------------
*              exchange reg1 and reg2
*-----------------------------------------------------------------------

               if( itrcn(m) .lt. 0 ) then

                        itrcn(m) = -itrcn(m)

                        kdsm = itrcr(m)
                        ldsm = -1
                        idsm = -1

                  do i = 1, abs( itrcn(m) )

                        ldsm = ldsm + 1
                        ntr1 = idas_itrcr(kdsm+ldsm)
                        ldsm = ldsm + 1
                        mtr1 = idas_itrcr(kdsm+ldsm)

                        mdsm = 1
                        call moddas_allocate_int(
     &                          mtr1, idas_itrcr_temporary1)

                     do k = 1, mtr1

                        ldsm = ldsm + 1
                        idas_itrcr_temporary1( mdsm + k - 1 )
     &                     = idas_itrcr(kdsm+ldsm)

                     end do

                        ldsm = ldsm + 1
                        ntr2 = idas_itrcr(kdsm+ldsm)
                        ldsm = ldsm + 1
                        mtr2 = idas_itrcr(kdsm+ldsm)

                        ndsm = 1
                        call moddas_allocate_int(
     &                          mtr2, idas_itrcr_temporary2)

                     do k = 1, mtr2

                        ldsm = ldsm + 1
                        idas_itrcr_temporary2( ndsm + k - 1 )
     &                     = idas_itrcr(kdsm+ldsm)

                     end do

                        idsm = idsm + 1
                        idas_itrcr(kdsm+idsm) = ntr2
                        idsm = idsm + 1
                        idas_itrcr(kdsm+idsm) = mtr2

                     do k = 1, mtr2

                        idsm = idsm + 1
                        idas_itrcr(kdsm+idsm)
     &                     = idas_itrcr_temporary2( ndsm + k - 1 )

                     end do

                        idsm = idsm + 1
                        idas_itrcr(kdsm+idsm) = ntr1
                        idsm = idsm + 1
                        idas_itrcr(kdsm+idsm) = mtr1

                     do k = 1, mtr1

                        idsm = idsm + 1
                        idas_itrcr(kdsm+idsm)
     &                     = idas_itrcr_temporary1( mdsm + k - 1 )

                     end do
                     call moddas_deallocate_int(idas_itrcr_temporary1)
                     call moddas_deallocate_int(idas_itrcr_temporary2)

                  end do

               end if

*-----------------------------------------------------------------------

                     itrss(m) = itrcs(m)

                     kdsm = itrcc(m)
                     ldsm = -1

                  call moddas_reallocate_int(
     &                    itlmax, m, MAX_NUM_ITRCR, itrcr, idas_itrcr)
                     idsm = itrcr(m)
                     jdsm = -1

               do i = 1, itrcn(m) * 2

                     ldsm = ldsm + 1
                     ntrn = idas_itrcc(kdsm+ldsm)
                     ldsm = ldsm + 1
                     mtrn = idas_itrcc(kdsm+ldsm)
                     ign  = kdsm + ldsm + 1
                     ldsm = ldsm + mtrn

                     igm  = idsm + jdsm + 3

                     call tregion2(io,jo,ierr,ntrn,mtrn,idas_itrcc(ign)
     &                             ,igm,MAX_NUM_ITRCR,idas_itrcr)

                     if( ierr .ne. 0 ) then

                        write(io,'('' **** above error in '',i2,
     &                  ''-th tally'')') m

                        ErrCha = ''
                        MsgID = 'L:1113/R:settal/F:tallsm1.f'
                        call ErrWrite(MsgID, ErrCha)
                        write(jo,'('' **** above error in '',i2,
     &                  ''-th tally'')') m

                        goto 999

                     end if

                     jdsm = jdsm + 1
                     idas_itrcr(idsm+jdsm) = ntrn

                     jdsm = jdsm + 1
                     idas_itrcr(idsm+jdsm) = mtrn

                     jdsm = jdsm + mtrn

               end do

                  itrcs(m) = jdsm

               if( jdsm > MAX_NUM_ITRCR ) then
                  write(ErrCha,'(a,a,i5,a,a,i5,a)')
     &                 'sub.settal@tallsm1.f'
     &                    //' ?dimension over idas_itrcr?'
     &                    //' jdsm > MAX_NUM_ITRCR'
     &                 ,' (jdsm=',jdsm,')'
     &                 ,' (MAX_NUM_ITRCR@moddas.f=',MAX_NUM_ITRCR,')'
                  ErrID = 'L:1141/R:settal/F:tallsm1.f'
                  call ErrWrite(ErrID,ErrCha)
               endif

               call moddas_reduce_int(
     &                 itlmax, m, jdsm+1, itrcr, idas_itrcr)

                  itsmn(m) = itsmn(m) + ( jdsm + mod(jdsm,2) ) / 2 + 1

            end if

*-----------------------------------------------------------------------
*           specify the subroutine number and storage position
*-----------------------------------------------------------------------

               italm(m) = mmmax

            if( itmsh(m) .eq. 1 ) then

               itals(m) = 4

            else if( itmsh(m) .eq. 2 ) then

               itals(m) = 5

            else if( itmsh(m) .eq. 3 ) then

               itals(m) = 12

            end if

*-----------------------------------------------------------------------
*           multiplier
*-----------------------------------------------------------------------

            if( itmst(m) .eq. 0 ) itmst(m) = 1

*-----------------------------------------------------------------------
*        yield tally
*-----------------------------------------------------------------------

         else if( ital(m) .eq. 3 ) then

*-----------------------------------------------------------------------
*              check region
*-----------------------------------------------------------------------

            if( itmsh(m) .eq. 1 ) then

                  igm  = 1

                  itrcm(m) = itrgm(m)

                  ntrn = itrgn(m)
                  mtrn = itrgm(m)
                  ign  = itrcg(m)
                  call tregion2(io,jo,ierr,ntrn,mtrn,idas_itrcg(ign)
     &                          ,igm,MAX_NUM_ITREG,idas_itreg_temporary)

                  if( mtrn > MAX_NUM_ITREG ) then
                     write(ErrCha,'(a,a,i5,a,a,i5,a)')
     &                    'sub.settal t-yield@tallsm1.f'
     &                       //' ?dimension over idas_itreg_temporary?'
     &                       //' mtrn > MAX_NUM_ITREG'
     &                    ,' (mtrn=',mtrn,')'
     &                    ,' (MAX_NUM_ITREG@moddas.f=',MAX_NUM_ITREG,')'
                     ErrID = 'L:1207/R:settal/F:tallsm1.f'
                     call ErrWrite(ErrID,ErrCha)
                  endif

                     if( ierr .ne. 0 ) then

                        write(io,'('' **** above error in '',i2,
     &                  ''-th tally'')') m

                        ErrCha = ''
                        MsgID = 'L:1217/R:settal/F:tallsm1.f'
                        call ErrWrite(MsgID, ErrCha)
                        write(jo,'('' **** above error in '',i2,
     &                  ''-th tally'')') m

                        goto 999

                     end if

                  call moddas_reallocate_int(
     &                    itlmax, m, mtrn, itreg, idas_itreg)
                  idas_itreg(itreg(m):itreg(m)+mtrn-1)
     &                      = idas_itreg_temporary(1:mtrn)
                  itrgn(m) = ntrn
                  itrgm(m) = mtrn

                  itsmn(m) = itsmn(m) + ( mtrn + mod(mtrn,2) ) / 2


            else if( itmsh(m) .eq. 4 ) then

                  itrcm(m) = itrgm(m)

                  ntrn = itrgn(m)
                  mtrn = itrgm(m)

                  ign  = itrcg(m)

                  call moddas_allocate_int(nelemtot,ielem2ir)

                  call ttetmesh2(io,jo,ierr,ntrn,mtrn,idas_itrcg(ign),
     &                 MAX_NUM_ITREG,idas_itreg_temporary,ielem2ir)

                  if( mtrn > MAX_NUM_ITREG ) then
                     write(ErrCha,'(a,a,i5,a,a,i5,a)')
     &                    'sub.settal t-track@tallsm1.f'
     &                       //' ?dimension over idas_itreg_temporary?'
     &                       //' mtrn > MAX_NUM_ITREG'
     &                    ,' (mtrn=',mtrn,')'
     &                    ,' (MAX_NUM_ITREG@moddas.f=',MAX_NUM_ITREG,')'
                     ErrID = 'L:1257/R:settal/F:tallsm1.f'
                     call ErrWrite(ErrID,ErrCha)
                  endif

                     if( ierr .ne. 0 ) then

                        write(io,'('' **** above error in '',i2,
     &                  ''-th tally'')') m

                        ErrCha = ''
                        MsgID = 'L:1267/R:settal/F:tallsm1.f'
                        call ErrWrite(MsgID, ErrCha)
                        write(jo,'('' **** above error in '',i2,
     &                  ''-th tally'')') m

                        goto 999

                     end if

                  if(ntrn.gt.0)then

                   call moddas_reallocate_int(
     &                  itlmax, m, mtrn+ntrn,
     &                  itreg, idas_itreg)
                   idas_itreg(itreg(m):itreg(m)+mtrn-1)
     &                  = idas_itreg_temporary(1:mtrn)
                   idas_itreg(itreg(m)+mtrn:itreg(m)+mtrn+ntrn-1)
     &                  = ielem2ir(1:ntrn)

                   itrgn(m) = ntrn
                   itrgm(m) = mtrn+ntrn

                  else

                   ntrn0=idas_itreg_temporary(4)
                   idas_itreg_temporary(4)=0

                   call moddas_reallocate_int(
     &                  itlmax, m, mtrn, itreg, idas_itreg)
                   idas_itreg(itreg(m):itreg(m)+mtrn-1)
     &                  = idas_itreg_temporary(1:mtrn)
                   itrgn(m) = ntrn0
                   itrgm(m) = mtrn

                  endif

                  call moddas_deallocate_int(ielem2ir)

            end if

*-----------------------------------------------------------------------
*              check material
*-----------------------------------------------------------------------

               if( itmtn(m) .gt. 0 ) then

                  do i = 1, itmtn(m)

                        imat = ismte( itmtt(m) + i - 1 )

                     if( imat .gt. 0 ) then
                     if( idnm( imat ) .eq. 0 ) then

                        write(io,'('' **** Error in '',i2,
     &                  ''-th tally : material ['',i4,
     &                  ''] is not defined '',
     &                  ''in material section.'')') m, imat

                        ErrCha = ''
                        MsgID = 'L:1326/R:settal/F:tallsm1.f'
                        call ErrWrite(MsgID, ErrCha)
                        write(jo,'('' **** Error in '',i2,
     &                  ''-th tally : material ['',i4,
     &                  ''] is not defined '',
     &                  ''in material section.'')') m, imat

                        goto 999

                     end if
                     end if

                  end do

               end if

*-----------------------------------------------------------------------
*           specify the subroutine number and storage position
*-----------------------------------------------------------------------

               italm(m) = mmmax

               maxpn    = itndz(m)
               maxnn    = itndn(m)

            if( itmsh(m) .eq. 1 ) then

               itals(m) = 6

               if( itrsh(m) .ne. 0 .and. itrgn(m) .le. 1 ) goto 940

            else if( itmsh(m) .eq. 2 ) then

               itals(m) = 7

            else if( itmsh(m) .eq. 3 ) then

               itals(m) = 8

            else if( itmsh(m) .eq. 4 ) then
               itals(m) = 51

            end if

*-----------------------------------------------------------------------
*        dchain tally
*-----------------------------------------------------------------------

         else if( ital(m) .eq. 16 ) then

*            write(*,'(5x,''tallsm1.f_settal_line1117'')')

*-----------------------------------------------------------------------
*           specify the subroutine number and storage position
*-----------------------------------------------------------------------

            if( itmsh(m-2) .eq. 1 ) then

               itals(m) = 38

            else if( itmsh(m-2) .eq. 2 ) then

CCSE change for mesh=r-z parameter (2017.11.30) >>>>>
               itals(m) = 39   ! reserved for dchain (r-z mesh)
CCSE change for mesh=r-z parameter (2017.11.30) <<<<<
            else if( itmsh(m-2) .eq. 3 ) then

CCSE change for mesh=xyz parameter (2017.11.30) >>>>>
               itals(m) = 40   ! reserved for dchain (xyz mesh)
CCSE change for mesh=xyz parameter (2017.11.30) <<<<<
            else if( itmsh(m-2) .eq. 4 ) then
               itals(m) = 54   ! reserved for dchain (tet mesh)

            end if

*-----------------------------------------------------------------------
*        heat tally
*-----------------------------------------------------------------------

         else if( ital(m) .eq. 4 ) then

*-----------------------------------------------------------------------
*              check region
*-----------------------------------------------------------------------

            if( itmsh(m) .eq. 1 ) then

                  igm  = 1

                  itrcm(m) = itrgm(m)

                  ntrn = itrgn(m)
                  mtrn = itrgm(m)
                  ign  = itrcg(m)
                  call tregion2(io,jo,ierr,ntrn,mtrn,idas_itrcg(ign)
     &                          ,igm,MAX_NUM_ITREG,idas_itreg_temporary)

                  if( mtrn > MAX_NUM_ITREG ) then
                     write(ErrCha,'(a,a,i5,a,a,i5,a)')
     &                    'sub.settal t-heat@tallsm1.f'
     &                       //' ?dimension over idas_itreg_temporary?'
     &                       //' mtrn > MAX_NUM_ITREG'
     &                    ,' (mtrn=',mtrn,')'
     &                    ,' (MAX_NUM_ITREG@moddas.f=',MAX_NUM_ITREG,')'
                     ErrID = 'L:1430/R:settal/F:tallsm1.f'
                     call ErrWrite(ErrID,ErrCha)
                  endif

                     if( ierr .ne. 0 ) then

                        write(io,'('' **** above error in '',i2,
     &                  ''-th tally'')') m

                        ErrCha = ''
                        MsgID = 'L:1440/R:settal/F:tallsm1.f'
                        call ErrWrite(MsgID, ErrCha)
                        write(jo,'('' **** above error in '',i2,
     &                  ''-th tally'')') m

                        goto 999

                     end if

                  call moddas_reallocate_int(
     &                    itlmax, m, mtrn, itreg, idas_itreg)
                  idas_itreg(itreg(m):itreg(m)+mtrn-1)
     &                      = idas_itreg_temporary(1:mtrn)
                  itrgn(m) = ntrn
                  itrgm(m) = mtrn

                  itsmn(m) = itsmn(m) + ( mtrn + mod(mtrn,2) ) / 2


            end if

*-----------------------------------------------------------------------
*              check material
*-----------------------------------------------------------------------

               if( itmtn(m) .gt. 0 ) then

                  do i = 1, itmtn(m)

                        imat = ismte( itmtt(m) + i - 1 )

                     if( imat .gt. 0 ) then
                     if( idnm( imat ) .eq. 0 ) then

                        write(io,'('' **** Error in '',i2,
     &                  ''-th tally : material ['',i4,
     &                  ''] is not defined '',
     &                  ''in material section.'')') m, imat

                        ErrCha = ''
                        MsgID = 'L:1480/R:settal/F:tallsm1.f'
                        call ErrWrite(MsgID, ErrCha)
                        write(jo,'('' **** Error in '',i2,
     &                  ''-th tally : material ['',i4,
     &                  ''] is not defined '',
     &                  ''in material section.'')') m, imat

                        goto 999

                     end if
                     end if

                  end do

               end if

*-----------------------------------------------------------------------
*           specify the subroutine number and storage position
*-----------------------------------------------------------------------

               italm(m) = mmmax

            if( itmsh(m) .eq. 1 ) then

               itals(m) = 9

               if( itrsh(m) .ne. 0 .and. itrgn(m) .le. 1 ) goto 940

            else if( itmsh(m) .eq. 2 ) then

               itals(m) = 10

            else if( itmsh(m) .eq. 3 ) then

               itals(m) = 11

            end if

               ithet(m) = mmmax  !FURUTA20131001
               mmmax = mmmax + 5 !FURUTA20131001

               if( mmmax .gt. mdas ) goto 950

            do i = italm(m), mmmax - 1

               das(i) = 0.0d0

            end do

*-----------------------------------------------------------------------
*        star tally
*-----------------------------------------------------------------------

         else if( ital(m) .eq. 5 ) then

*-----------------------------------------------------------------------
*           open file for dump data
*-----------------------------------------------------------------------

            if( itmdp(m,0) .ne. 0 ) then

                  idmpf = idmpf + 1

                  if( idmpf .gt. 20 ) then

                        write(io,'('' **** Error in '',i2,
     &                  ''-th tally : number of dump file''
     &                  '' exceeds 20 files'')') m

                        ErrCha = ''
                        MsgID = 'L:1550/R:settal/F:tallsm1.f'
                        call ErrWrite(MsgID, ErrCha)
                        write(jo,'('' **** Error in '',i2,
     &                  ''-th tally : number of dump file''
     &                  '' exceeds 20 files'')') m

                        goto 999

                  end if

                  itmdf(m) = idmpf + 80

C S.H. revised for Dump-Restart on 2014/5/7
                  filnm = ctfln(m,1)

                  do i = itfll(m,1), 1, -1
                     if ( ctfln(m,1)(i:i) .eq. '.' ) exit
                  end do
                  if ( i .eq. 0 ) i = itfll(m,1)+1

                  filnm(1:i-1) = ctfln(m,1)(1:i-1)
                  filnm(i:i+3) = '_dmp'
                  if ( i .lt. itfll(m,1) ) then
                     filnm(i+4:itfll(m,1)+4) = ctfln(m,1)(i:itfll(m,1))
                  end if

               if( npe .le. 1 ) then

                  if( itmdp(m,0) .gt. 0 ) then

                     if ( irestart .eq. 1 ) then

                        open(itmdf(m), file = filnm,
     &                       form='unformatted', status = 'unknown',
     &                       access = 'append' )

                     else

                        open(itmdf(m), file = filnm,
     &                       form='unformatted', status = 'unknown' )

                     end if

                  else

                     if ( irestart .eq. 1 ) then

                        open(itmdf(m), file = filnm,
     &                       form='formatted', status = 'unknown',
     &                       access = 'append' )
                     else

                        open(itmdf(m), file = filnm,
     &                       form='formatted', status = 'unknown' )

                     end if

                  end if

               else if( me .gt. 0 ) then

                     iorder = aint(log10(real(npe)))+1
                     if ( iorder .lt. 3) iorder = 3

                     write(chme,'(i5.5)') me

                     lengfilnm = len_trim(filnm)
                     filnm(lengfilnm+1:lengfilnm+1+iorder)
     &                    = '.' // chme(6-iorder:5)

                  if( itmdp(m,0) .gt. 0 ) then

                     if ( irestart .eq. 1 ) then

                        open(itmdf(m), file = filnm,
     &                       form='unformatted', status = 'unknown',
     &                       access = 'append' )

                     else

                        open(itmdf(m), file = filnm,
     &                       form='unformatted', status = 'unknown' )

                     end if

                  else

                     if ( irestart .eq. 1 ) then

                        open(itmdf(m), file = filnm,
     &                       form='formatted', status = 'unknown',
     &                       access = 'append' )

                     else

                        open(itmdf(m), file = filnm,
     &                       form='formatted', status = 'unknown' )

                     end if

                  end if

               end if

               if(npe.le.1.or.me.gt.0)then           !FURUTA20150427
                call chkdmpomp(io,m,itmdp(m,0),ierr) !FURUTA20150427
                if(ierr.ne.0)goto 999                !FURUTA20150427
               endif                                 !FURUTA20150427

            end if

*-----------------------------------------------------------------------
*              check region
*-----------------------------------------------------------------------

            if( itmsh(m) .eq. 1 ) then

                  igm  = 1

                  itrcm(m) = itrgm(m)

                  ntrn = itrgn(m)
                  mtrn = itrgm(m)
                  ign  = itrcg(m)
                  call tregion2(io,jo,ierr,ntrn,mtrn,idas_itrcg(ign)
     &                          ,igm,MAX_NUM_ITREG,idas_itreg_temporary)

                  if( mtrn > MAX_NUM_ITREG ) then
                     write(ErrCha,'(a,a,i5,a,a,i5,a)')
     &                    'sub.settal t-star@tallsm1.f'
     &                       //' ?dimension over idas_itreg_temporary?'
     &                       //' mtrn > MAX_NUM_ITREG'
     &                    ,' (mtrn=',mtrn,')'
     &                    ,' (MAX_NUM_ITREG@moddas.f=',MAX_NUM_ITREG,')'
                     ErrID = 'L:1684/R:settal/F:tallsm1.f'
                     call ErrWrite(ErrID,ErrCha)
                  endif

                     if( ierr .ne. 0 ) then

                        write(io,'('' **** above error in '',i2,
     &                  ''-th tally'')') m

                        ErrCha = ''
                        MsgID = 'L:1694/R:settal/F:tallsm1.f'
                        call ErrWrite(MsgID, ErrCha)
                        write(jo,'('' **** above error in '',i2,
     &                  ''-th tally'')') m

                        goto 999

                     end if

                  call moddas_reallocate_int(
     &                    itlmax, m, mtrn, itreg, idas_itreg)
                  idas_itreg(itreg(m):itreg(m)+mtrn-1)
     &                      = idas_itreg_temporary(1:mtrn)
                  itrgn(m) = ntrn
                  itrgm(m) = mtrn

                  itsmn(m) = itsmn(m) + ( mtrn + mod(mtrn,2) ) / 2


            end if

*-----------------------------------------------------------------------
*              check material
*-----------------------------------------------------------------------

               if( itmtn(m) .gt. 0 ) then

                  do i = 1, itmtn(m)

                        imat = ismte( itmtt(m) + i - 1 )

                     if( imat .gt. 0 ) then
                     if( idnm( imat ) .eq. 0 ) then

                        write(io,'('' **** Error in '',i2,
     &                  ''-th tally : material ['',i4,
     &                  ''] is not defined '',
     &                  ''in material section.'')') m, imat

                        ErrCha = ''
                        MsgID = 'L:1734/R:settal/F:tallsm1.f'
                        call ErrWrite(MsgID, ErrCha)
                        write(jo,'('' **** Error in '',i2,
     &                  ''-th tally : material ['',i4,
     &                  ''] is not defined '',
     &                  ''in material section.'')') m, imat

                        goto 999

                     end if
                     end if

                  end do

               end if

*-----------------------------------------------------------------------
*           specify the subroutine number and storage position
*-----------------------------------------------------------------------

               italm(m) = mmmax

            if( itmsh(m) .eq. 1 ) then

               itals(m) = 13

               if( itrsh(m) .ne. 0 .and. itrgn(m) .le. 1 ) goto 940

            else if( itmsh(m) .eq. 2 ) then

               itals(m) = 14

            else if( itmsh(m) .eq. 3 ) then

               itals(m) = 15

            end if

*-----------------------------------------------------------------------
*        time tally
*-----------------------------------------------------------------------

         else if( ital(m) .eq. 6 ) then

*-----------------------------------------------------------------------
*           open file for dump data
*-----------------------------------------------------------------------

            if( itmdp(m,0) .ne. 0 ) then

                  idmpf = idmpf + 1

                  if( idmpf .gt. 20 ) then

                        write(io,'('' **** Error in '',i2,
     &                  ''-th tally : number of dump file''
     &                  '' exceeds 20 files'')') m

                        ErrCha = ''
                        MsgID = 'L:1793/R:settal/F:tallsm1.f'
                        call ErrWrite(MsgID, ErrCha)
                        write(jo,'('' **** Error in '',i2,
     &                  ''-th tally : number of dump file''
     &                  '' exceeds 20 files'')') m

                        goto 999

                  end if

                  itmdf(m) = idmpf + 80

C S.H. revised for Dump-Restart on 2014/5/7
                  filnm = ctfln(m,1)

                  do i = itfll(m,1), 1, -1
                     if ( ctfln(m,1)(i:i) .eq. '.' ) exit
                  end do
                  if ( i .eq. 0 ) i = itfll(m,1)+1

                  filnm(1:i-1) = ctfln(m,1)(1:i-1)
                  filnm(i:i+3) = '_dmp'
                  if ( i .lt. itfll(m,1) ) then
                     filnm(i+4:itfll(m,1)+4) = ctfln(m,1)(i:itfll(m,1))
                  end if

               if( npe .le. 1 ) then

                  if( itmdp(m,0) .gt. 0 ) then

                     if ( irestart .eq. 1 ) then

                        open(itmdf(m), file = filnm,
     &                       form='unformatted', status = 'unknown',
     &                       access = 'append' )

                     else

                        open(itmdf(m), file = filnm,
     &                       form='unformatted', status = 'unknown' )

                     end if

                  else

                     if ( irestart .eq. 1 ) then

                        open(itmdf(m), file = filnm,
     &                       form='formatted', status = 'unknown',
     &                       access = 'append' )
                     else

                        open(itmdf(m), file = filnm,
     &                       form='formatted', status = 'unknown' )

                     end if

                  end if

               else if( me .gt. 0 ) then

                     iorder = aint(log10(real(npe)))+1
                     if ( iorder .lt. 3) iorder = 3

                     write(chme,'(i5.5)') me

                     lengfilnm = len_trim(filnm)
                     filnm(lengfilnm+1:lengfilnm+1+iorder)
     &                    = '.' // chme(6-iorder:5)

                  if( itmdp(m,0) .gt. 0 ) then

                     if ( irestart .eq. 1 ) then

                        open(itmdf(m), file = filnm,
     &                       form='unformatted', status = 'unknown',
     &                       access = 'append' )

                     else

                        open(itmdf(m), file = filnm,
     &                       form='unformatted', status = 'unknown' )

                     end if

                  else

                     if ( irestart .eq. 1 ) then

                        open(itmdf(m), file = filnm,
     &                       form='formatted', status = 'unknown',
     &                       access = 'append' )

                     else

                        open(itmdf(m), file = filnm,
     &                       form='formatted', status = 'unknown' )

                     end if

                  end if

               end if

               if(npe.le.1.or.me.gt.0)then           !FURUTA20150427
                call chkdmpomp(io,m,itmdp(m,0),ierr) !FURUTA20150427
                if(ierr.ne.0)goto 999                !FURUTA20150427
               endif                                 !FURUTA20150427

            end if

*-----------------------------------------------------------------------
*              check region
*-----------------------------------------------------------------------

            if( itmsh(m) .eq. 1 ) then

                  igm  = 1

                  itrcm(m) = itrgm(m)

                  ntrn = itrgn(m)
                  mtrn = itrgm(m)
                  ign  = itrcg(m)
                  call tregion2(io,jo,ierr,ntrn,mtrn,idas_itrcg(ign)
     &                          ,igm,MAX_NUM_ITREG,idas_itreg_temporary)

                  if( mtrn > MAX_NUM_ITREG ) then
                     write(ErrCha,'(a,a,i5,a,a,i5,a)')
     &                    'sub.settal t-time@tallsm1.f'
     &                       //' ?dimension over idas_itreg_temporary?'
     &                       //' mtrn > MAX_NUM_ITREG'
     &                    ,' (mtrn=',mtrn,')'
     &                    ,' (MAX_NUM_ITREG@moddas.f=',MAX_NUM_ITREG,')'
                     ErrID = 'L:1927/R:settal/F:tallsm1.f'
                     call ErrWrite(ErrID,ErrCha)
                  endif

                     if( ierr .ne. 0 ) then

                        write(io,'('' **** above error in '',i2,
     &                  ''-th tally'')') m

                        ErrCha = ''
                        MsgID = 'L:1937/R:settal/F:tallsm1.f'
                        call ErrWrite(MsgID, ErrCha)
                        write(jo,'('' **** above error in '',i2,
     &                  ''-th tally'')') m

                        goto 999

                     end if

                  call moddas_reallocate_int(
     &                    itlmax, m, mtrn, itreg, idas_itreg)
                  idas_itreg(itreg(m):itreg(m)+mtrn-1)
     &                      = idas_itreg_temporary(1:mtrn)
                  itrgn(m) = ntrn
                  itrgm(m) = mtrn

                  itsmn(m) = itsmn(m) + ( mtrn + mod(mtrn,2) ) / 2


            end if

*-----------------------------------------------------------------------
*              check material
*-----------------------------------------------------------------------

               if( itmtn(m) .gt. 0 ) then

                  do i = 1, itmtn(m)

                        imat = ismte( itmtt(m) + i - 1 )

                     if( imat .gt. 0 ) then
                     if( idnm( imat ) .eq. 0 ) then

                        write(io,'('' **** Error in '',i2,
     &                  ''-th tally : material ['',i4,
     &                  ''] is not defined '',
     &                  ''in material section.'')') m, imat

                        ErrCha = ''
                        MsgID = 'L:1977/R:settal/F:tallsm1.f'
                        call ErrWrite(MsgID, ErrCha)
                        write(jo,'('' **** Error in '',i2,
     &                  ''-th tally : material ['',i4,
     &                  ''] is not defined '',
     &                  ''in material section.'')') m, imat

                        goto 999

                     end if
                     end if

                  end do

               end if

*-----------------------------------------------------------------------
*           specify the subroutine number and storage position
*-----------------------------------------------------------------------

               italm(m) = mmmax

            if( itmsh(m) .eq. 1 ) then

               itals(m) = 16

               if( itrsh(m) .ne. 0 .and. itrgn(m) .le. 1 ) goto 940

            else if( itmsh(m) .eq. 2 ) then

               itals(m) = 17

            else if( itmsh(m) .eq. 3 ) then

               itals(m) = 18

            end if

*-----------------------------------------------------------------------
*        dpa tally
*-----------------------------------------------------------------------

         else if( ital(m) .eq. 7 ) then

*-----------------------------------------------------------------------
*              check region
*-----------------------------------------------------------------------

            if( itmsh(m) .eq. 1 ) then

                  igm  = 1

                  itrcm(m) = itrgm(m)

                  ntrn = itrgn(m)
                  mtrn = itrgm(m)
                  ign  = itrcg(m)
                  call tregion2(io,jo,ierr,ntrn,mtrn,idas_itrcg(ign)
     &                          ,igm,MAX_NUM_ITREG,idas_itreg_temporary)

                  if( mtrn > MAX_NUM_ITREG ) then
                     write(ErrCha,'(a,a,i5,a,a,i5,a)')
     &                    'sub.settal t-dpa@tallsm1.f'
     &                       //' ?dimension over idas_itreg_temporary?'
     &                       //' mtrn > MAX_NUM_ITREG'
     &                    ,' (mtrn=',mtrn,')'
     &                    ,' (MAX_NUM_ITREG@moddas.f=',MAX_NUM_ITREG,')'
                     ErrID = 'L:2044/R:settal/F:tallsm1.f'
                     call ErrWrite(ErrID,ErrCha)
                  endif

                     if( ierr .ne. 0 ) then

                        write(io,'('' **** above error in '',i2,
     &                  ''-th tally'')') m

                        ErrCha = ''
                        MsgID = 'L:2054/R:settal/F:tallsm1.f'
                        call ErrWrite(MsgID, ErrCha)
                        write(jo,'('' **** above error in '',i2,
     &                  ''-th tally'')') m

                        goto 999

                     end if

                  call moddas_reallocate_int(
     &                    itlmax, m, mtrn, itreg, idas_itreg)
                  idas_itreg(itreg(m):itreg(m)+mtrn-1)
     &                      = idas_itreg_temporary(1:mtrn)
                  itrgn(m) = ntrn
                  itrgm(m) = mtrn

                  itsmn(m) = itsmn(m) + ( mtrn + mod(mtrn,2) ) / 2


            else if( itmsh(m) .eq. 4 ) then

                  itrcm(m) = itrgm(m)

                  ntrn = itrgn(m)
                  mtrn = itrgm(m)

                  ign  = itrcg(m)

                  call moddas_allocate_int(nelemtot,ielem2ir)

                  call ttetmesh2(io,jo,ierr,ntrn,mtrn,idas_itrcg(ign),
     &                 MAX_NUM_ITREG,idas_itreg_temporary,ielem2ir)

                  if( mtrn > MAX_NUM_ITREG ) then
                     write(ErrCha,'(a,a,i5,a,a,i5,a)')
     &                    'sub.settal t-track@tallsm1.f'
     &                       //' ?dimension over idas_itreg_temporary?'
     &                       //' mtrn > MAX_NUM_ITREG'
     &                    ,' (mtrn=',mtrn,')'
     &                    ,' (MAX_NUM_ITREG@moddas.f=',MAX_NUM_ITREG,')'
                     ErrID = 'L:2094/R:settal/F:tallsm1.f'
                     call ErrWrite(ErrID,ErrCha)
                  endif

                     if( ierr .ne. 0 ) then

                        write(io,'('' **** above error in '',i2,
     &                  ''-th tally'')') m

                        ErrCha = ''
                        MsgID = 'L:2104/R:settal/F:tallsm1.f'
                        call ErrWrite(MsgID, ErrCha)
                        write(jo,'('' **** above error in '',i2,
     &                  ''-th tally'')') m

                        goto 999

                     end if

                  if(ntrn.gt.0)then

                   call moddas_reallocate_int(
     &                  itlmax, m, mtrn+ntrn,
     &                  itreg, idas_itreg)
                   idas_itreg(itreg(m):itreg(m)+mtrn-1)
     &                  = idas_itreg_temporary(1:mtrn)
                   idas_itreg(itreg(m)+mtrn:itreg(m)+mtrn+ntrn-1)
     &                  = ielem2ir(1:ntrn)

                   itrgn(m) = ntrn
                   itrgm(m) = mtrn+ntrn

                  else

                   ntrn0=idas_itreg_temporary(4)
                   idas_itreg_temporary(4)=0

                   call moddas_reallocate_int(
     &                  itlmax, m, mtrn, itreg, idas_itreg)
                   idas_itreg(itreg(m):itreg(m)+mtrn-1)
     &                  = idas_itreg_temporary(1:mtrn)
                   itrgn(m) = ntrn0
                   itrgm(m) = mtrn

                  endif

                  call moddas_deallocate_int(ielem2ir)

            end if

*-----------------------------------------------------------------------
*              check material
*-----------------------------------------------------------------------

               if( itmtn(m) .gt. 0 ) then

                  do i = 1, itmtn(m)

                        imat = ismte( itmtt(m) + i - 1 )

                     if( imat .gt. 0 ) then
                     if( idnm( imat ) .eq. 0 ) then

                        write(io,'('' **** Error in '',i2,
     &                  ''-th tally : material ['',i4,
     &                  ''] is not defined '',
     &                  ''in material section.'')') m, imat

                        ErrCha = ''
                        MsgID = 'L:2163/R:settal/F:tallsm1.f'
                        call ErrWrite(MsgID, ErrCha)
                        write(jo,'('' **** Error in '',i2,
     &                  ''-th tally : material ['',i4,
     &                  ''] is not defined '',
     &                  ''in material section.'')') m, imat

                        goto 999

                     end if
                     end if

                  end do

               end if

*-----------------------------------------------------------------------
*              check proton library
*-----------------------------------------------------------------------

               if( itln(m,1) .gt. 0 ) then

                  do i = 1, itln(m,1)

                        imat = mlib(itli(m,1),  i)
                        jmat = mlib(itli(m,1)+1,i)

                     if( idnm( imat ) .eq. 0 ) then

                        write(io,'('' **** Error in '',i2,
     &                  ''-th tally : material ['',i4,
     &                  ''] is not defined '',
     &                  ''in material section.'')') m, imat

                        ErrCha = ''
                        MsgID = 'L:2198/R:settal/F:tallsm1.f'
                        call ErrWrite(MsgID, ErrCha)
                        write(jo,'('' **** Error in '',i2,
     &                  ''-th tally : material ['',i4,
     &                  ''] is not defined '',
     &                  ''in material section.'')') m, imat

                        goto 999

                     end if

                     if( idnm( jmat ) .eq. 0 ) then

                        write(io,'('' **** Error in '',i2,
     &                  ''-th tally : library material ['',i4,
     &                  ''] is not defined '',
     &                  ''in material section.'')') m, jmat

                        ErrCha = ''
                        MsgID = 'L:2217/R:settal/F:tallsm1.f'
                        call ErrWrite(MsgID, ErrCha)
                        write(jo,'('' **** Error in '',i2,
     &                  ''-th tally : library material ['',i4,
     &                  ''] is not defined '',
     &                  ''in material section.'')') m, jmat

                        goto 999

                     end if

                  end do

               end if

*-----------------------------------------------------------------------
*              check neutron library
*-----------------------------------------------------------------------

               if( itln(m,2) .gt. 0 ) then

                  do i = 1, itln(m,2)

                        imat = mlib(itli(m,2),  i)
                        jmat = mlib(itli(m,2)+1,i)

                     if( idnm( imat ) .eq. 0 ) then

                        write(io,'('' **** Error in '',i2,
     &                  ''-th tally : material ['',i4,
     &                  ''] is not defined '',
     &                  ''in material section.'')') m, imat

                        ErrCha = ''
                        MsgID = 'L:2251/R:settal/F:tallsm1.f'
                        call ErrWrite(MsgID, ErrCha)
                        write(jo,'('' **** Error in '',i2,
     &                  ''-th tally : material ['',i4,
     &                  ''] is not defined '',
     &                  ''in material section.'')') m, imat

                        goto 999

                     end if

                     if( idnm( jmat ) .eq. 0 ) then

                        write(io,'('' **** Error in '',i2,
     &                  ''-th tally : library material ['',i4,
     &                  ''] is not defined '',
     &                  ''in material section.'')') m, jmat

                        ErrCha = ''
                        MsgID = 'L:2270/R:settal/F:tallsm1.f'
                        call ErrWrite(MsgID, ErrCha)
                        write(jo,'('' **** Error in '',i2,
     &                  ''-th tally : library material ['',i4,
     &                  ''] is not defined '',
     &                  ''in material section.'')') m, jmat

                        goto 999

                     end if

                  end do

               end if

*-----------------------------------------------------------------------
*           specify the subroutine number and storage position
*-----------------------------------------------------------------------

               italm(m) = mmmax

            if( itmsh(m) .eq. 1 ) then

               itals(m) = 19

               if( itrsh(m) .ne. 0 .and. itrgn(m) .le. 1 ) goto 940

            else if( itmsh(m) .eq. 2 ) then

               itals(m) = 20

            else if( itmsh(m) .eq. 3 ) then

               itals(m) = 21

            else if( itmsh(m) .eq. 4 ) then
               itals(m) = 53

            end if

*-----------------------------------------------------------------------
*        product tally
*-----------------------------------------------------------------------

         else if( ital(m) .eq. 8 ) then

*-----------------------------------------------------------------------
*           open file for dump data
*-----------------------------------------------------------------------

            if( itmdp(m,0) .ne. 0 ) then

                  idmpf = idmpf + 1

                  if( idmpf .gt. 20 ) then

                        write(io,'('' **** Error in '',i2,
     &                  ''-th tally : number of dump file''
     &                  '' exceeds 20 files'')') m

                        ErrCha = ''
                        MsgID = 'L:2331/R:settal/F:tallsm1.f'
                        call ErrWrite(MsgID, ErrCha)
                        write(jo,'('' **** Error in '',i2,
     &                  ''-th tally : number of dump file''
     &                  '' exceeds 20 files'')') m

                        goto 999

                  end if

                  itmdf(m) = idmpf + 80

C S.H. revised for Dump-Restart on 2014/5/7
                  filnm = ctfln(m,1)

                  do i = itfll(m,1), 1, -1
                     if ( ctfln(m,1)(i:i) .eq. '.' ) exit
                  end do
                  if ( i .eq. 0 ) i = itfll(m,1)+1

                  filnm(1:i-1) = ctfln(m,1)(1:i-1)
                  filnm(i:i+3) = '_dmp'
                  if ( i .lt. itfll(m,1) ) then
                     filnm(i+4:itfll(m,1)+4) = ctfln(m,1)(i:itfll(m,1))
                  end if

               if( npe .le. 1 ) then

                  if( itmdp(m,0) .gt. 0 ) then

                     if ( irestart .eq. 1 ) then

                        open(itmdf(m), file = filnm,
     &                       form='unformatted', status = 'unknown',
     &                       access = 'append' )

                     else

                        open(itmdf(m), file = filnm,
     &                       form='unformatted', status = 'unknown' )

                     end if

                  else

                     if ( irestart .eq. 1 ) then

                        open(itmdf(m), file = filnm,
     &                       form='formatted', status = 'unknown',
     &                       access = 'append' )
                     else

                        open(itmdf(m), file = filnm,
     &                       form='formatted', status = 'unknown' )

                     end if

                  end if

               else if( me .gt. 0 ) then

                     iorder = aint(log10(real(npe)))+1
                     if ( iorder .lt. 3) iorder = 3

                     write(chme,'(i5.5)') me

                     lengfilnm = len_trim(filnm)
                     filnm(lengfilnm+1:lengfilnm+1+iorder)
     &                    = '.' // chme(6-iorder:5)

                  if( itmdp(m,0) .gt. 0 ) then

                     if ( irestart .eq. 1 ) then

                        open(itmdf(m), file = filnm,
     &                       form='unformatted', status = 'unknown',
     &                       access = 'append' )

                     else

                        open(itmdf(m), file = filnm,
     &                       form='unformatted', status = 'unknown' )

                     end if

                  else

                     if ( irestart .eq. 1 ) then

                        open(itmdf(m), file = filnm,
     &                       form='formatted', status = 'unknown',
     &                       access = 'append' )

                     else

                        open(itmdf(m), file = filnm,
     &                       form='formatted', status = 'unknown' )

                     end if

                  end if

               end if

               if(npe.le.1.or.me.gt.0)then           !FURUTA20150427
                call chkdmpomp(io,m,itmdp(m,0),ierr) !FURUTA20150427
                if(ierr.ne.0)goto 999                !FURUTA20150427
               endif                                 !FURUTA20150427

            end if

*-----------------------------------------------------------------------
*              check region
*-----------------------------------------------------------------------

            if( itmsh(m) .eq. 1 ) then

                  igm  = 1

                  itrcm(m) = itrgm(m)

                  ntrn = itrgn(m)
                  mtrn = itrgm(m)
                  ign  = itrcg(m)
                  call tregion2(io,jo,ierr,ntrn,mtrn,idas_itrcg(ign)
     &                          ,igm,MAX_NUM_ITREG,idas_itreg_temporary)

                  if( mtrn > MAX_NUM_ITREG ) then
                     write(ErrCha,'(a,a,i5,a,a,i5,a)')
     &                    'sub.settal t-product@tallsm1.f'
     &                       //' ?dimension over idas_itreg_temporary?'
     &                       //' mtrn > MAX_NUM_ITREG'
     &                    ,' (mtrn=',mtrn,')'
     &                    ,' (MAX_NUM_ITREG@moddas.f=',MAX_NUM_ITREG,')'
                     ErrID = 'L:2465/R:settal/F:tallsm1.f'
                     call ErrWrite(ErrID,ErrCha)
                  endif

                     if( ierr .ne. 0 ) then

                        write(io,'('' **** above error in '',i2,
     &                  ''-th tally'')') m

                        ErrCha = ''
                        MsgID = 'L:2475/R:settal/F:tallsm1.f'
                        call ErrWrite(MsgID, ErrCha)
                        write(jo,'('' **** above error in '',i2,
     &                  ''-th tally'')') m

                        goto 999

                     end if

                  call moddas_reallocate_int(
     &                    itlmax, m, mtrn, itreg, idas_itreg)
                  idas_itreg(itreg(m):itreg(m)+mtrn-1)
     &                      = idas_itreg_temporary(1:mtrn)
                  itrgn(m) = ntrn
                  itrgm(m) = mtrn

                  itsmn(m) = itsmn(m) + ( mtrn + mod(mtrn,2) ) / 2


            else if( itmsh(m) .eq. 4 ) then

                  itrcm(m) = itrgm(m)

                  ntrn = itrgn(m)
                  mtrn = itrgm(m)

                  ign  = itrcg(m)

                  call moddas_allocate_int(nelemtot,ielem2ir)

                  call ttetmesh2(io,jo,ierr,ntrn,mtrn,idas_itrcg(ign),
     &                 MAX_NUM_ITREG,idas_itreg_temporary,ielem2ir)

                  if( mtrn > MAX_NUM_ITREG ) then
                     write(ErrCha,'(a,a,i5,a,a,i5,a)')
     &                    'sub.settal t-track@tallsm1.f'
     &                       //' ?dimension over idas_itreg_temporary?'
     &                       //' mtrn > MAX_NUM_ITREG'
     &                    ,' (mtrn=',mtrn,')'
     &                    ,' (MAX_NUM_ITREG@moddas.f=',MAX_NUM_ITREG,')'
                     ErrID = 'L:2515/R:settal/F:tallsm1.f'
                     call ErrWrite(ErrID,ErrCha)
                  endif

                     if( ierr .ne. 0 ) then

                        write(io,'('' **** above error in '',i2,
     &                  ''-th tally'')') m

                        ErrCha = ''
                        MsgID = 'L:2525/R:settal/F:tallsm1.f'
                        call ErrWrite(MsgID, ErrCha)
                        write(jo,'('' **** above error in '',i2,
     &                  ''-th tally'')') m

                        goto 999

                     end if

                  if(ntrn.gt.0)then

                   call moddas_reallocate_int(
     &                  itlmax, m, mtrn+ntrn,
     &                  itreg, idas_itreg)
                   idas_itreg(itreg(m):itreg(m)+mtrn-1)
     &                  = idas_itreg_temporary(1:mtrn)
                   idas_itreg(itreg(m)+mtrn:itreg(m)+mtrn+ntrn-1)
     &                  = ielem2ir(1:ntrn)

                   itrgn(m) = ntrn
                   itrgm(m) = mtrn+ntrn

                  else

                   ntrn0=idas_itreg_temporary(4)
                   idas_itreg_temporary(4)=0

                   call moddas_reallocate_int(
     &                  itlmax, m, mtrn, itreg, idas_itreg)
                   idas_itreg(itreg(m):itreg(m)+mtrn-1)
     &                  = idas_itreg_temporary(1:mtrn)
                   itrgn(m) = ntrn0
                   itrgm(m) = mtrn

                  endif

                  call moddas_deallocate_int(ielem2ir)

            end if

*-----------------------------------------------------------------------
*              check material
*-----------------------------------------------------------------------

               if( itmtn(m) .gt. 0 ) then

                  do i = 1, itmtn(m)

                        imat = ismte( itmtt(m) + i - 1 )

                     if( imat .gt. 0 ) then
                     if( idnm( imat ) .eq. 0 ) then

                        write(io,'('' **** Error in '',i2,
     &                  ''-th tally : material ['',i4,
     &                  ''] is not defined '',
     &                  ''in material section.'')') m, imat

                        ErrCha = ''
                        MsgID = 'L:2584/R:settal/F:tallsm1.f'
                        call ErrWrite(MsgID, ErrCha)
                        write(jo,'('' **** Error in '',i2,
     &                  ''-th tally : material ['',i4,
     &                  ''] is not defined '',
     &                  ''in material section.'')') m, imat

                        goto 999

                     end if
                     end if

                  end do

               end if

*-----------------------------------------------------------------------
*           specify the subroutine number and storage position
*-----------------------------------------------------------------------

               italm(m) = mmmax

            if( itmsh(m) .eq. 1 ) then

               itals(m) = 22

               if( itrsh(m) .ne. 0 .and. itrgn(m) .le. 1 ) goto 940

            else if( itmsh(m) .eq. 2 ) then

               itals(m) = 23

            else if( itmsh(m) .eq. 3 ) then

               itals(m) = 24

            else if( itmsh(m) .eq. 4 ) then
               itals(m) = 52

            end if

*-----------------------------------------------------------------------
*        g-show tally
*-----------------------------------------------------------------------

         else if( ital(m) .eq. 9 ) then

*-----------------------------------------------------------------------

               itals(m) = 25

               italm(m) = 0

            if( itaxs(m,1) .eq. 1 ) then

               immax = ( itxnm(m) + 1 ) * ( itynm(m) + 1 )

            else if( itaxs(m,1) .eq. 2 ) then

               immax = ( itynm(m) + 1 ) * ( itznm(m) + 1 )

            else if( itaxs(m,1) .eq. 3 ) then

               immax = ( itxnm(m) + 1 ) * ( itznm(m) + 1 )

            end if

               immax    = ( immax + mod(immax,2) ) / 2
               lmmax    = max( lmmax, immax )

               itstm(m) = itstm(m) + immax

*-----------------------------------------------------------------------
*        r-show tally
*-----------------------------------------------------------------------

         else if( ital(m) .eq. 10 ) then

*-----------------------------------------------------------------------

               itals(m) = 26

               italm(m) = 0

*-----------------------------------------------------------------------
*              check region
*-----------------------------------------------------------------------

                  igm  = 1

                  itrcm(m) = itrgm(m)

                  ntrn = itrgn(m)
                  mtrn = itrgm(m)
                  ign  = itrcg(m)

               do k = 1, mtrn

                  ntrkg = idas_itrcg( ign + k - 1 )

                  if( ntrkg .eq. 60000 ) then

                        write(io,'('' *** Error at '',i2,
     &                  ''-th [t-rshow] tally.''/
     &                  ''  (all) cannot be used'')') m

                        ErrCha = ''
                        MsgID = 'L:2691/R:settal/F:tallsm1.f'
                        call ErrWrite(MsgID, ErrCha)
                        write(jo,'('' *** Error at '',i2,
     &                  ''-th [t-rshow] tally.''/
     &                  ''  (all) cannot be used'')') m

                        goto 999

                  end if

               end do

                  call tregion2(io,jo,ierr,ntrn,mtrn,idas_itrcg(ign)
     &                          ,igm,MAX_NUM_ITREG,idas_itreg_temporary)

                  if( mtrn > MAX_NUM_ITREG ) then
                     write(ErrCha,'(a,a,i5,a,a,i5,a)')
     &                    'sub.settal t-rshow@tallsm1.f'
     &                       //' ?dimension over idas_itreg_temporary?'
     &                       //' mtrn > MAX_NUM_ITREG'
     &                    ,' (mtrn=',mtrn,')'
     &                    ,' (MAX_NUM_ITREG@moddas.f=',MAX_NUM_ITREG,')'
                     ErrID = 'L:2713/R:settal/F:tallsm1.f'
                     call ErrWrite(ErrID,ErrCha)
                  endif

                     if( ierr .ne. 0 ) then

                        write(io,'('' **** above error in '',i2,
     &                  ''-th tally'')') m

                        ErrCha = ''
                        MsgID = 'L:2723/R:settal/F:tallsm1.f'
                        call ErrWrite(MsgID, ErrCha)
                        write(jo,'('' **** above error in '',i2,
     &                  ''-th tally'')') m

                        goto 999

                     end if

                  call moddas_reallocate_int(
     &                    itlmax, m, mtrn, itreg, idas_itreg)
                  idas_itreg(itreg(m):itreg(m)+mtrn-1)
     &                      = idas_itreg_temporary(1:mtrn)
                  itrgn(m) = ntrn
                  itrgm(m) = mtrn

                  itsmn(m) = itsmn(m) + ( mtrn + mod(mtrn,2) ) / 2


*-----------------------------------------------------------------------

            if( itaxs(m,1) .eq. 1 ) then

               immax = ( itxnm(m) + 1 ) * ( itynm(m) + 1 )

            else if( itaxs(m,1) .eq. 2 ) then

               immax = ( itynm(m) + 1 ) * ( itznm(m) + 1 )

            else if( itaxs(m,1) .eq. 3 ) then

               immax = ( itxnm(m) + 1 ) * ( itznm(m) + 1 )

            end if

               immax = immax + itrgn(m)

               immax    = itrgn(m) + ( immax + mod(immax,2) ) / 2
               lmmax    = max( lmmax, immax )

               itstm(m) = itstm(m) + immax

*-----------------------------------------------------------------------
*        3d-show tally
*-----------------------------------------------------------------------

         else if( ital(m) .eq. 11 ) then

*-----------------------------------------------------------------------

               itals(m) = 27

               italm(m) = 0

               icrsflx = 1

*-----------------------------------------------------------------------
*              check region
*-----------------------------------------------------------------------

            if( itrgn(m) .gt. 0 ) then

                  igm  = 1

                  itrcm(m) = itrgm(m)

                  ntrn = itrgn(m)
                  mtrn = itrgm(m)
                  ign  = itrcg(m)

               do k = 1, mtrn

                  ntrkg = idas_itrcg( ign + k - 1 )

                  if( ntrkg .eq. 60000 ) then

                        write(io,'('' *** Error at '',i2,
     &                  ''-th [t-3dshow] tally.''/
     &                  ''  (all) cannot be used'')') m

                        ErrCha = ''
                        MsgID = 'L:2804/R:settal/F:tallsm1.f'
                        call ErrWrite(MsgID, ErrCha)
                        write(jo,'('' *** Error at '',i2,
     &                  ''-th [t-3dshow] tally.''/
     &                  ''  (all) cannot be used'')') m

                        goto 999

                  end if

               end do

                  call tregion2(io,jo,ierr,ntrn,mtrn,idas_itrcg(ign)
     &                          ,igm,MAX_NUM_ITREG,idas_itreg_temporary)

                  if( mtrn > MAX_NUM_ITREG ) then
                     write(ErrCha,'(a,a,i5,a,a,i5,a)')
     &                    'sub.settal t-3dshow@tallsm1.f'
     &                       //' ?dimension over idas_itreg_temporary?'
     &                       //' mtrn > MAX_NUM_ITREG'
     &                    ,' (mtrn=',mtrn,')'
     &                    ,' (MAX_NUM_ITREG@moddas.f=',MAX_NUM_ITREG,')'
                     ErrID = 'L:2826/R:settal/F:tallsm1.f'
                     call ErrWrite(ErrID,ErrCha)
                  endif

                     if( ierr .ne. 0 ) then

                        write(io,'('' **** above error in '',i2,
     &                  ''-th tally'')') m

                        ErrCha = ''
                        MsgID = 'L:2836/R:settal/F:tallsm1.f'
                        call ErrWrite(MsgID, ErrCha)
                        write(jo,'('' **** above error in '',i2,
     &                  ''-th tally'')') m

                        goto 999

                     end if

                  call moddas_reallocate_int(
     &                    itlmax, m, mtrn, itreg, idas_itreg)
                  idas_itreg(itreg(m):itreg(m)+mtrn-1)
     &                      = idas_itreg_temporary(1:mtrn)
                  itrgn(m) = ntrn
                  itrgm(m) = mtrn

                  itsmn(m) = itsmn(m) + ( mtrn + mod(mtrn,2) ) / 2


            end if

*-----------------------------------------------------------------------
*              check region in box
*-----------------------------------------------------------------------

            if( itmgn(m) .gt. 0 ) then


                  itmcm(m) = itmgm(m)

                  ntrn = itmgn(m)
                  mtrn = itmgm(m)
                  ign  = itmcg(m)

               do k = 1, mtrn

                  ntrkg = idas_itmcg( ign + k - 1 )

                  if( ntrkg .eq. 60000 ) then

                        write(io,'('' *** Error at '',i2,
     &                  ''-th [t-3dshow] tally.''/
     &                  ''  (all) cannot be used'')') m

                        ErrCha = ''
                        MsgID = 'L:2881/R:settal/F:tallsm1.f'
                        call ErrWrite(MsgID, ErrCha)
                        write(jo,'('' *** Error at '',i2,
     &                  ''-th [t-3dshow] tally.''/
     &                  ''  (all) cannot be used'')') m

                        goto 999

                  end if

               end do

                  call moddas_reallocate_int(
     &                    itlmax, m, MAX_NUM_ITMEG, itmeg, idas_itmeg)
                  igm = itmeg(m)
                  call tregion2(io,jo,ierr,ntrn,mtrn,idas_itmcg(ign)
     &                          ,igm,MAX_NUM_ITMEG,idas_itmeg)

                     if( ierr .ne. 0 ) then

                        write(io,'('' **** above error in '',i2,
     &                  ''-th tally'')') m

                        ErrCha = ''
                        MsgID = 'L:2905/R:settal/F:tallsm1.f'
                        call ErrWrite(MsgID, ErrCha)
                        write(jo,'('' **** above error in '',i2,
     &                  ''-th tally'')') m

                        goto 999

                     end if

                  itmeg(m) = igm
                  itmgn(m) = ntrn
                  itmgm(m) = mtrn

               if( mtrn > MAX_NUM_ITMEG ) then
                  write(ErrCha,'(a,a,i5,a,a,i5,a)')
     &                 'sub.settal@tallsm1.f'
     &                    //' ?dimension over idas_itmeg?'
     &                    //' mtrn > MAX_NUM_ITMEG'
     &                 ,' (mtrn=',mtrn,')'
     &                 ,' (MAX_NUM_ITMEG@moddas.f=',MAX_NUM_ITMEG,')'
                  ErrID = 'L:2925/R:settal/F:tallsm1.f'
                  call ErrWrite(ErrID,ErrCha)
               endif

               call moddas_reduce_int(
     &                 itlmax, m, mtrn+1, itmeg, idas_itmeg)
                  itsmn(m) = itsmn(m) + ( mtrn + mod(mtrn,2) ) / 2


            end if

*-----------------------------------------------------------------------
*              check material
*-----------------------------------------------------------------------

               if( itmtn(m) .gt. 0 ) then

                  do i = 1, itmtn(m)

                        imat = ismte( itmtt(m) + i - 1 )

                     if( imat .gt. 0 ) then
                     if( idnm( imat ) .eq. 0 ) then

                        write(io,'('' **** Error in '',i2,
     &                  ''-th tally : material ['',i4,
     &                  ''] is not defined '',
     &                  ''in material section.'')') m, imat

                        ErrCha = ''
                        MsgID = 'L:2955/R:settal/F:tallsm1.f'
                        call ErrWrite(MsgID, ErrCha)
                        write(jo,'('' **** Error in '',i2,
     &                  ''-th tally : material ['',i4,
     &                  ''] is not defined '',
     &                  ''in material section.'')') m, imat

                        goto 999

                     end if
                     end if

                  end do

               end if

*-----------------------------------------------------------------------
*              check material in box
*-----------------------------------------------------------------------

               if( itmbn(m) .gt. 0 ) then

                  do i = 1, itmbn(m)

                        imat = ismte_itmbt( itmbt(m) + i - 1 )

                     if( imat .gt. 0 ) then
                     if( idnm( imat ) .eq. 0 ) then

                        write(io,'('' **** Error in '',i2,
     &                  ''-th tally : material ['',i4,
     &                  ''] is not defined '',
     &                  ''in material section.'')') m, imat

                        ErrCha = ''
                        MsgID = 'L:2990/R:settal/F:tallsm1.f'
                        call ErrWrite(MsgID, ErrCha)
                        write(jo,'('' **** Error in '',i2,
     &                  ''-th tally : material ['',i4,
     &                  ''] is not defined '',
     &                  ''in material section.'')') m, imat

                        goto 999

                     end if
                     end if

                  end do

               end if

*-----------------------------------------------------------------------

               immax = ( itwin(m,1) + 1 ) * ( itwin(m,2) + 1 ) * 9
               immax = ( immax + mod(immax,2) ) / 2

               lmmax = max( lmmax, immax )

               itstm(m) = itstm(m) + immax

*-----------------------------------------------------------------------
*        t-let tally
*-----------------------------------------------------------------------

         else if( ital(m) .eq. 12 ) then

*-----------------------------------------------------------------------
*              check region
*-----------------------------------------------------------------------

            if( itmsh(m) .eq. 1 ) then

                  igm  = 1

                  itrcm(m) = itrgm(m)

                  ntrn = itrgn(m)
                  mtrn = itrgm(m)
                  ign  = itrcg(m)
                  call tregion2(io,jo,ierr,ntrn,mtrn,idas_itrcg(ign)
     &                          ,igm,MAX_NUM_ITREG,idas_itreg_temporary)

                  if( mtrn > MAX_NUM_ITREG ) then
                     write(ErrCha,'(a,a,i5,a,a,i5,a)')
     &                    'sub.settal t-let@tallsm1.f'
     &                       //' ?dimension over idas_itreg_temporary?'
     &                       //' mtrn > MAX_NUM_ITREG'
     &                    ,' (mtrn=',mtrn,')'
     &                    ,' (MAX_NUM_ITREG@moddas.f=',MAX_NUM_ITREG,')'
                     ErrID = 'L:3044/R:settal/F:tallsm1.f'
                     call ErrWrite(ErrID,ErrCha)
                  endif

                     if( ierr .ne. 0 ) then

                        write(io,'('' **** above error in '',i2,
     &                  ''-th tally'')') m

                        ErrCha = ''
                        MsgID = 'L:3054/R:settal/F:tallsm1.f'
                        call ErrWrite(MsgID, ErrCha)
                        write(jo,'('' **** above error in '',i2,
     &                  ''-th tally'')') m

                        goto 999

                     end if

                  call moddas_reallocate_int(
     &                    itlmax, m, mtrn, itreg, idas_itreg)
                  idas_itreg(itreg(m):itreg(m)+mtrn-1)
     &                      = idas_itreg_temporary(1:mtrn)
                  itrgn(m) = ntrn
                  itrgm(m) = mtrn

                  itsmn(m) = itsmn(m) + ( mtrn + mod(mtrn,2) ) / 2


            end if

*-----------------------------------------------------------------------
*              check LET material
*-----------------------------------------------------------------------

               if( itlmt(m) .gt. 0 ) then

                     imat = itlmt(m)

                     if( imat .ne. 0 ) then
                     if( idnm( abs(imat) ) .eq. 0 ) then

                        write(io,'('' **** Error in '',i2,
     &                  ''-th tally : let material ['',i4,
     &                  ''] is not defined '',
     &                  ''in material section.'')') m, imat

                        ErrCha = ''
                        MsgID = 'L:3092/R:settal/F:tallsm1.f'
                        call ErrWrite(MsgID, ErrCha)
                        write(jo,'('' **** Error in '',i2,
     &                  ''-th tally : let material ['',i4,
     &                  ''] is not defined '',
     &                  ''in material section.'')') m, imat

                        goto 999

                     end if
                     end if

               end if

*-----------------------------------------------------------------------
*              check material
*-----------------------------------------------------------------------

               if( itmtn(m) .gt. 0 ) then

                  do i = 1, itmtn(m)

                        imat = ismte( itmtt(m) + i - 1 )

                     if( imat .gt. 0 ) then
                     if( idnm( imat ) .eq. 0 ) then

                        write(io,'('' **** Error in '',i2,
     &                  ''-th tally : material ['',i4,
     &                  ''] is not defined '',
     &                  ''in material section.'')') m, imat

                        ErrCha = ''
                        MsgID = 'L:3125/R:settal/F:tallsm1.f'
                        call ErrWrite(MsgID, ErrCha)
                        write(jo,'('' **** Error in '',i2,
     &                  ''-th tally : material ['',i4,
     &                  ''] is not defined '',
     &                  ''in material section.'')') m, imat

                        goto 999

                     end if
                     end if

                  end do

               end if

*-----------------------------------------------------------------------
*           specify the subroutine number and storage position
*-----------------------------------------------------------------------

               italm(m) = mmmax

            if( itmsh(m) .eq. 1 ) then

               itals(m) = 28

               if( itrsh(m) .ne. 0 .and. itrgn(m) .le. 1 ) goto 940

            else if( itmsh(m) .eq. 2 ) then

               itals(m) = 29

            else if( itmsh(m) .eq. 3 ) then

               itals(m) = 30

            end if

*-----------------------------------------------------------------------
*        t-deposit tally
*-----------------------------------------------------------------------

         else if( ital(m) .eq. 13 ) then

*-----------------------------------------------------------------------
*              check region
*-----------------------------------------------------------------------

            if( itmsh(m) .eq. 1 ) then

                  igm  = 1

                  itrcm(m) = itrgm(m)

                  ntrn = itrgn(m)
                  mtrn = itrgm(m)
                  ign  = itrcg(m)
                  call tregion2(io,jo,ierr,ntrn,mtrn,idas_itrcg(ign)
     &                          ,igm,MAX_NUM_ITREG,idas_itreg_temporary)

                  if( mtrn > MAX_NUM_ITREG ) then
                     write(ErrCha,'(a,a,i5,a,a,i5,a)')
     &                    'sub.settal t-deposit@tallsm1.f'
     &                       //' ?dimension over idas_itreg_temporary?'
     &                       //' mtrn > MAX_NUM_ITREG'
     &                    ,' (mtrn=',mtrn,')'
     &                    ,' (MAX_NUM_ITREG@moddas.f=',MAX_NUM_ITREG,')'
                     ErrID = 'L:3192/R:settal/F:tallsm1.f'
                     call ErrWrite(ErrID,ErrCha)
                  endif

                     if( ierr .ne. 0 ) then

                        write(io,'('' **** above error in '',i2,
     &                  ''-th tally'')') m

                        ErrCha = ''
                        MsgID = 'L:3202/R:settal/F:tallsm1.f'
                        call ErrWrite(MsgID, ErrCha)
                        write(jo,'('' **** above error in '',i2,
     &                  ''-th tally'')') m

                        goto 999

                     end if

                  call moddas_reallocate_int(
     &                    itlmax, m, mtrn, itreg, idas_itreg)
                  idas_itreg(itreg(m):itreg(m)+mtrn-1)
     &                      = idas_itreg_temporary(1:mtrn)
                  itrgn(m) = ntrn
                  itrgm(m) = mtrn

                  itsmn(m) = itsmn(m) + ( mtrn + mod(mtrn,2) ) / 2


            else if( itmsh(m) .eq. 4 ) then

                  itrcm(m) = itrgm(m)

                  ntrn = itrgn(m)
                  mtrn = itrgm(m)

                  ign  = itrcg(m)

                  call moddas_allocate_int(nelemtot,ielem2ir)

                  call ttetmesh2(io,jo,ierr,ntrn,mtrn,idas_itrcg(ign),
     &                 MAX_NUM_ITREG,idas_itreg_temporary,ielem2ir)

                  if( mtrn > MAX_NUM_ITREG ) then
                     write(ErrCha,'(a,a,i5,a,a,i5,a)')
     &                    'sub.settal t-track@tallsm1.f'
     &                       //' ?dimension over idas_itreg_temporary?'
     &                       //' mtrn > MAX_NUM_ITREG'
     &                    ,' (mtrn=',mtrn,')'
     &                    ,' (MAX_NUM_ITREG@moddas.f=',MAX_NUM_ITREG,')'
                     ErrID = 'L:3242/R:settal/F:tallsm1.f'
                     call ErrWrite(ErrID,ErrCha)
                  endif

                     if( ierr .ne. 0 ) then

                        write(io,'('' **** above error in '',i2,
     &                  ''-th tally'')') m

                        ErrCha = ''
                        MsgID = 'L:3252/R:settal/F:tallsm1.f'
                        call ErrWrite(MsgID, ErrCha)
                        write(jo,'('' **** above error in '',i2,
     &                  ''-th tally'')') m

                        goto 999

                     end if

                  if(ntrn.gt.0)then

                   call moddas_reallocate_int(
     &                  itlmax, m, mtrn+ntrn,
     &                  itreg, idas_itreg)
                   idas_itreg(itreg(m):itreg(m)+mtrn-1)
     &                  = idas_itreg_temporary(1:mtrn)
                   idas_itreg(itreg(m)+mtrn:itreg(m)+mtrn+ntrn-1)
     &                  = ielem2ir(1:ntrn)

                   itrgn(m) = ntrn
                   itrgm(m) = mtrn+ntrn

                  else

                   ntrn0=idas_itreg_temporary(4)
                   idas_itreg_temporary(4)=0

                   call moddas_reallocate_int(
     &                  itlmax, m, mtrn, itreg, idas_itreg)
                   idas_itreg(itreg(m):itreg(m)+mtrn-1)
     &                  = idas_itreg_temporary(1:mtrn)
                   itrgn(m) = ntrn0
                   itrgm(m) = mtrn

                  endif

                  call moddas_deallocate_int(ielem2ir)

            end if

*-----------------------------------------------------------------------
*              check LET material
*-----------------------------------------------------------------------

               if( itlmt(m) .gt. 0 ) then

                     imat = itlmt(m)

                     if( imat .ne. 0 ) then
                     if( idnm( abs(imat) ) .eq. 0 ) then

                        write(io,'('' **** Error in '',i2,
     &                  ''-th tally : let material ['',i4,
     &                  ''] is not defined '',
     &                  ''in material section.'')') m, imat

                        ErrCha = ''
                        MsgID = 'L:3309/R:settal/F:tallsm1.f'
                        call ErrWrite(MsgID, ErrCha)
                        write(jo,'('' **** Error in '',i2,
     &                  ''-th tally : let material ['',i4,
     &                  ''] is not defined '',
     &                  ''in material section.'')') m, imat

                        goto 999

                     end if
                     end if

               end if

*-----------------------------------------------------------------------
*              check material
*-----------------------------------------------------------------------

               if( itmtn(m) .gt. 0 ) then

                  do i = 1, itmtn(m)

                        imat = ismte( itmtt(m) + i - 1 )

                     if( imat .gt. 0 ) then
                     if( idnm( imat ) .eq. 0 ) then

                        write(io,'('' **** Error in '',i2,
     &                  ''-th tally : material ['',i4,
     &                  ''] is not defined '',
     &                  ''in material section.'')') m, imat

                        ErrCha = ''
                        MsgID = 'L:3342/R:settal/F:tallsm1.f'
                        call ErrWrite(MsgID, ErrCha)
                        write(jo,'('' **** Error in '',i2,
     &                  ''-th tally : material ['',i4,
     &                  ''] is not defined '',
     &                  ''in material section.'')') m, imat

                        goto 999

                     end if
                     end if

                  end do

               end if

*-----------------------------------------------------------------------
*           specify the subroutine number and storage position
*-----------------------------------------------------------------------

               italm(m) = mmmax

            if( itmsh(m) .eq. 1 ) then

               itals(m) = 31

               if( itrsh(m) .ne. 0 .and. itrgn(m) .le. 1 ) goto 940

            else if( itmsh(m) .eq. 2 ) then

               itals(m) = 32

            else if( itmsh(m) .eq. 3 ) then

               itals(m) = 33

            else if( itmsh(m) .eq. 4 ) then
               itals(m) = 50

            end if

*-----------------------------------------------------------------------
*        t-deposit2 tally
*-----------------------------------------------------------------------

         else if( ital(m) .eq. 14 ) then

*-----------------------------------------------------------------------
*              check region
*-----------------------------------------------------------------------

            if( itmsh(m) .eq. 1 ) then

                  igm  = 1

                  itrcm(m) = itrgm(m)

                  ntrn = itrgn(m)
                  mtrn = itrgm(m)
                  ign  = itrcg(m)
                  call tregion2(io,jo,ierr,ntrn,mtrn,idas_itrcg(ign)
     &                          ,igm,MAX_NUM_ITREG,idas_itreg_temporary)

                  if( mtrn > MAX_NUM_ITREG ) then
                     write(ErrCha,'(a,a,i5,a,a,i5,a)')
     &                    'sub.settal t-deposit2@tallsm1.f'
     &                       //' ?dimension over idas_itreg_temporary?'
     &                       //' mtrn > MAX_NUM_ITREG'
     &                    ,' (mtrn=',mtrn,')'
     &                    ,' (MAX_NUM_ITREG@moddas.f=',MAX_NUM_ITREG,')'
                     ErrID = 'L:3412/R:settal/F:tallsm1.f'
                     call ErrWrite(ErrID,ErrCha)
                  endif

                     if( ierr .ne. 0 ) then

                        write(io,'('' **** above error in '',i2,
     &                  ''-th tally'')') m

                        ErrCha = ''
                        MsgID = 'L:3422/R:settal/F:tallsm1.f'
                        call ErrWrite(MsgID, ErrCha)
                        write(jo,'('' **** above error in '',i2,
     &                  ''-th tally'')') m

                        goto 999

                     end if

                  call moddas_reallocate_int(
     &                    itlmax, m, mtrn, itreg, idas_itreg)
                  idas_itreg(itreg(m):itreg(m)+mtrn-1)
     &                      = idas_itreg_temporary(1:mtrn)
                  itrgn(m) = ntrn
                  itrgm(m) = mtrn

                  itsmn(m) = itsmn(m) + ( mtrn + mod(mtrn,2) ) / 2


            end if

*-----------------------------------------------------------------------
*              check LET material
*-----------------------------------------------------------------------

               if( itlmt(m) .gt. 0 ) then

                     imat = itlmt(m)

                     if( imat .ne. 0 ) then
                     if( idnm( abs(imat) ) .eq. 0 ) then

                        write(io,'('' **** Error in '',i2,
     &                  ''-th tally : let material ['',i4,
     &                  ''] is not defined '',
     &                  ''in material section.'')') m, imat

                        ErrCha = ''
                        MsgID = 'L:3460/R:settal/F:tallsm1.f'
                        call ErrWrite(MsgID, ErrCha)
                        write(jo,'('' **** Error in '',i2,
     &                  ''-th tally : let material ['',i4,
     &                  ''] is not defined '',
     &                  ''in material section.'')') m, imat

                        goto 999

                     end if
                     end if

               end if

               if( itlmt2(m) .gt. 0 ) then

                     imat = itlmt2(m)

                     if( imat .ne. 0 ) then
                     if( idnm( abs(imat) ) .eq. 0 ) then

                        write(io,'('' **** Error in '',i2,
     &                  ''-th tally : let material ['',i4,
     &                  ''] is not defined '',
     &                  ''in material section.'')') m, imat

                        ErrCha = ''
                        MsgID = 'L:3487/R:settal/F:tallsm1.f'
                        call ErrWrite(MsgID, ErrCha)
                        write(jo,'('' **** Error in '',i2,
     &                  ''-th tally : let material ['',i4,
     &                  ''] is not defined '',
     &                  ''in material section.'')') m, imat

                        goto 999

                     end if
                     end if

               end if

*-----------------------------------------------------------------------
*           specify the subroutine number and storage position
*-----------------------------------------------------------------------

               italm(m) = mmmax

            if( itmsh(m) .eq. 1 ) then

               itals(m) = 34

            end if

*-----------------------------------------------------------------------
*        t-sed tally
*-----------------------------------------------------------------------

         else if( ital(m) .eq. 15 ) then

*-----------------------------------------------------------------------
*              check region
*-----------------------------------------------------------------------

            if( itmsh(m) .eq. 1 ) then

                  igm  = 1

                  itrcm(m) = itrgm(m)

                  ntrn = itrgn(m)
                  mtrn = itrgm(m)
                  ign  = itrcg(m)
                  call tregion2(io,jo,ierr,ntrn,mtrn,idas_itrcg(ign)
     &                          ,igm,MAX_NUM_ITREG,idas_itreg_temporary)

                  if( mtrn > MAX_NUM_ITREG ) then
                     write(ErrCha,'(a,a,i5,a,a,i5,a)')
     &                    'sub.settal t-sed@tallsm1.f'
     &                       //' ?dimension over idas_itreg_temporary?'
     &                       //' mtrn > MAX_NUM_ITREG'
     &                    ,' (mtrn=',mtrn,')'
     &                    ,' (MAX_NUM_ITREG@moddas.f=',MAX_NUM_ITREG,')'
                     ErrID = 'L:3542/R:settal/F:tallsm1.f'
                     call ErrWrite(ErrID,ErrCha)
                  endif

                     if( ierr .ne. 0 ) then

                        write(io,'('' **** above error in '',i2,
     &                  ''-th tally'')') m

                        ErrCha = ''
                        MsgID = 'L:3552/R:settal/F:tallsm1.f'
                        call ErrWrite(MsgID, ErrCha)
                        write(jo,'('' **** above error in '',i2,
     &                  ''-th tally'')') m

                        goto 999

                     end if

                  call moddas_reallocate_int(
     &                    itlmax, m, mtrn, itreg, idas_itreg)
                  idas_itreg(itreg(m):itreg(m)+mtrn-1)
     &                      = idas_itreg_temporary(1:mtrn)
                  itrgn(m) = ntrn
                  itrgm(m) = mtrn

                  itsmn(m) = itsmn(m) + ( mtrn + mod(mtrn,2) ) / 2


            end if

*-----------------------------------------------------------------------
*              check LET material
*-----------------------------------------------------------------------

               if( itlmt(m) .gt. 0 ) then

                     imat = itlmt(m)

                     if( imat .ne. 0 ) then
                     if( idnm( abs(imat) ) .eq. 0 ) then

                        write(io,'('' **** Error in '',i2,
     &                  ''-th tally : let material ['',i4,
     &                  ''] is not defined '',
     &                  ''in material section.'')') m, imat

                        ErrCha = ''
                        MsgID = 'L:3590/R:settal/F:tallsm1.f'
                        call ErrWrite(MsgID, ErrCha)
                        write(jo,'('' **** Error in '',i2,
     &                  ''-th tally : let material ['',i4,
     &                  ''] is not defined '',
     &                  ''in material section.'')') m, imat

                        goto 999

                     end if
                     end if

               end if

*-----------------------------------------------------------------------
*              check material
*-----------------------------------------------------------------------

               if( itmtn(m) .gt. 0 ) then

                  do i = 1, itmtn(m)

                        imat = ismte( itmtt(m) + i - 1 )

                     if( imat .gt. 0 ) then
                     if( idnm( imat ) .eq. 0 ) then

                        write(io,'('' **** Error in '',i2,
     &                  ''-th tally : material ['',i4,
     &                  ''] is not defined '',
     &                  ''in material section.'')') m, imat

                        ErrCha = ''
                        MsgID = 'L:3623/R:settal/F:tallsm1.f'
                        call ErrWrite(MsgID, ErrCha)
                        write(jo,'('' **** Error in '',i2,
     &                  ''-th tally : material ['',i4,
     &                  ''] is not defined '',
     &                  ''in material section.'')') m, imat

                        goto 999

                     end if
                     end if

                  end do

               end if

*-----------------------------------------------------------------------
*           specify the subroutine number and storage position
*-----------------------------------------------------------------------

               italm(m) = mmmax

            if( itmsh(m) .eq. 1 ) then

               itals(m) = 35

               if( itrsh(m) .ne. 0 .and. itrgn(m) .le. 1 ) goto 940

            else if( itmsh(m) .eq. 2 ) then

               itals(m) = 36

            else if( itmsh(m) .eq. 3 ) then

               itals(m) = 37

            end if

*-----------------------------------------------------------------------
*           stop OpenMP or MPI execution without setting proden prameters
*           (prodenmn, prodenmx, nbproden)
*-----------------------------------------------------------------------

!$          if( mftal(m).gt.5 .and. itprodenchk(m).eq.0 ) then
!$             write(ErrCha,'("Error: OpenMP execution cannot ",
!$   &           "be used in test calculation to estimate scale of ",
!$   &           "probability density function.")')
!$             ErrID = 'L:800/R:settal/F:tallsm1.f'
!$             call ErrWrite(ErrID,ErrCha)
!$
!$             goto 999
!$
!$          endif

            if( npe.gt.1 .and. mftal(m).gt.5 .and.
     &         itprodenchk(m).eq.0 ) then
               write(ErrCha,'("Error: MPI execution cannot ",
     &           "be used in test calculation to estimate scale of ",
     &           "probability density function.")')
               ErrID = 'L:3682/R:settal/F:tallsm1.f'
               call ErrWrite(ErrID,ErrCha)

               goto 999

            endif

*-----------------------------------------------------------------------
*        t-point tally
*-----------------------------------------------------------------------

         else if( ital(m) .eq. 17 ) then

               itals(m) = 41
               itpont = 1

*-----------------------------------------------------------------------
*           multiplier
*-----------------------------------------------------------------------

            if( itmst(m) .eq. 0 ) itmst(m) = 1

*-----------------------------------------------------------------------
*        t-wwg tally
*-----------------------------------------------------------------------

         else if( ital(m) .eq. 18 ) then

*-----------------------------------------------------------------------
*              check region
*-----------------------------------------------------------------------

            if( itmsh(m) .eq. 1 ) then

                  igm  = 1

                  itrcm(m) = itrgm(m)

                  ntrn = itrgn(m)
                  mtrn = itrgm(m)
                  ign  = itrcg(m)
                  call tregion2(io,jo,ierr,ntrn,mtrn,idas_itrcg(ign)
     &                          ,igm,MAX_NUM_ITREG,idas_itreg_temporary)

                  if( mtrn > MAX_NUM_ITREG ) then
                     write(ErrCha,'(a,a,i5,a,a,i5,a)')
     &                    'sub.settal t-wwg@tallsm1.f'
     &                       //' ?dimension over idas_itreg_temporary?'
     &                       //' mtrn > MAX_NUM_ITREG'
     &                    ,' (mtrn=',mtrn,')'
     &                    ,' (MAX_NUM_ITREG@moddas.f=',MAX_NUM_ITREG,')'
                     ErrID = 'L:3733/R:settal/F:tallsm1.f'
                     call ErrWrite(ErrID,ErrCha)
                  endif

                     if( ierr .ne. 0 ) then

                        write(io,'('' **** above error in '',i2,
     &                  ''-th tally'')') m

                        ErrCha = ''
                        MsgID = 'L:3743/R:settal/F:tallsm1.f'
                        call ErrWrite(MsgID, ErrCha)
                        write(jo,'('' **** above error in '',i2,
     &                  ''-th tally'')') m

                        goto 999

                     end if

                  call moddas_reallocate_int(
     &                    itlmax, m, mtrn, itreg, idas_itreg)
                  idas_itreg(itreg(m):itreg(m)+mtrn-1)
     &                      = idas_itreg_temporary(1:mtrn)
                  itrgn(m) = ntrn
                  itrgm(m) = mtrn

                  itsmn(m) = itsmn(m) + ( mtrn + mod(mtrn,2) ) / 2


            else if( itmsh(m) .eq. 4 ) then

                  itrcm(m) = itrgm(m)

                  ntrn = itrgn(m)
                  mtrn = itrgm(m)

                  ign  = itrcg(m)

                  call moddas_allocate_int(nelemtot,ielem2ir)

                  call ttetmesh2(io,jo,ierr,ntrn,mtrn,idas_itrcg(ign),
     &                 MAX_NUM_ITREG,idas_itreg_temporary,ielem2ir)

                  if( mtrn > MAX_NUM_ITREG ) then
                     write(ErrCha,'(a,a,i5,a,a,i5,a)')
     &                    'sub.settal t-track@tallsm1.f'
     &                       //' ?dimension over idas_itreg_temporary?'
     &                       //' mtrn > MAX_NUM_ITREG'
     &                    ,' (mtrn=',mtrn,')'
     &                    ,' (MAX_NUM_ITREG@moddas.f=',MAX_NUM_ITREG,')'
                     ErrID = 'L:3783/R:settal/F:tallsm1.f'
                     call ErrWrite(ErrID,ErrCha)
                  endif

                     if( ierr .ne. 0 ) then

                        write(io,'('' **** above error in '',i2,
     &                  ''-th tally'')') m

                        ErrCha = ''
                        MsgID = 'L:3793/R:settal/F:tallsm1.f'
                        call ErrWrite(MsgID, ErrCha)
                        write(jo,'('' **** above error in '',i2,
     &                  ''-th tally'')') m

                        goto 999

                     end if

                  if(ntrn.gt.0)then

                   call moddas_reallocate_int(
     &                  itlmax, m, mtrn+ntrn,
     &                  itreg, idas_itreg)
                   idas_itreg(itreg(m):itreg(m)+mtrn-1)
     &                  = idas_itreg_temporary(1:mtrn)
                   idas_itreg(itreg(m)+mtrn:itreg(m)+mtrn+ntrn-1)
     &                  = ielem2ir(1:ntrn)

                   itrgn(m) = ntrn
                   itrgm(m) = mtrn+ntrn

                  else

                   ntrn0=idas_itreg_temporary(4)
                   idas_itreg_temporary(4)=0

                   call moddas_reallocate_int(
     &                  itlmax, m, mtrn, itreg, idas_itreg)
                   idas_itreg(itreg(m):itreg(m)+mtrn-1)
     &                  = idas_itreg_temporary(1:mtrn)
                   itrgn(m) = ntrn0
                   itrgm(m) = mtrn

                  endif

                  call moddas_deallocate_int(ielem2ir)

            end if

*-----------------------------------------------------------------------
*              check material
*-----------------------------------------------------------------------

               if( itmtn(m) .gt. 0 ) then

                  do i = 1, itmtn(m)

                        imat = ismte( itmtt(m) + i - 1 )

                     if( imat .gt. 0 ) then
                     if( idnm( imat ) .eq. 0 ) then

                        write(io,'('' **** Error in '',i2,
     &                  ''-th tally : material ['',i4,
     &                  ''] is not defined '',
     &                  ''in material section.'')') m, imat

                        ErrCha = ''
                        MsgID = 'L:3852/R:settal/F:tallsm1.f'
                        call ErrWrite(MsgID, ErrCha)
                        write(jo,'('' **** Error in '',i2,
     &                  ''-th tally : material ['',i4,
     &                  ''] is not defined '',
     &                  ''in material section.'')') m, imat

                        goto 999

                     end if
                     end if

                  end do

               end if

*-----------------------------------------------------------------------
*           multiplier
*-----------------------------------------------------------------------

            if( itmst(m) .eq. 0 ) itmst(m) = 1

*-----------------------------------------------------------------------
*           specify the subroutine number and storage position
*-----------------------------------------------------------------------

               italm(m) = mmmax

            if( itmsh(m) .eq. 1 ) then

               itals(m) = 42

               if( itrsh(m) .ne. 0 .and. itrgn(m) .le. 1 ) goto 940

            else if( itmsh(m) .eq. 3 ) then

               itals(m) = 48

            else if( itmsh(m) .eq. 4 ) then

               itals(m) = 55

            end if

*-----------------------------------------------------------------------
*        t-volume tally
*-----------------------------------------------------------------------

         else if( ital(m) .eq. 21 ) then

*-----------------------------------------------------------------------
*              check region
*-----------------------------------------------------------------------

                  igm  = 1

                  itrcm(m) = itrgm(m)

                  ntrn = itrgn(m)
                  mtrn = itrgm(m)
                  ign  = itrcg(m)
                  call tregion2(io,jo,ierr,ntrn,mtrn,idas_itrcg(ign)
     &                          ,igm,MAX_NUM_ITREG,idas_itreg_temporary)

                  if( mtrn > MAX_NUM_ITREG ) then
                     write(ErrCha,'(a,a,i5,a,a,i5,a)')
     &                    'sub.settal t-volume@tallsm1.f'
     &                       //' ?dimension over idas_itreg_temporary?'
     &                       //' mtrn > MAX_NUM_ITREG'
     &                    ,' (mtrn=',mtrn,')'
     &                    ,' (MAX_NUM_ITREG@moddas.f=',MAX_NUM_ITREG,')'
                     ErrID = 'L:3923/R:settal/F:tallsm1.f'
                     call ErrWrite(ErrID,ErrCha)
                  endif

                     if( ierr .ne. 0 ) then

                        write(io,'('' **** above error in '',i2,
     &                  ''-th tally'')') m

                        ErrCha = ''
                        MsgID = 'L:3933/R:settal/F:tallsm1.f'
                        call ErrWrite(MsgID, ErrCha)
                        write(jo,'('' **** above error in '',i2,
     &                  ''-th tally'')') m

                        goto 999

                     end if

                  call moddas_reallocate_int(
     &                    itlmax, m, mtrn, itreg, idas_itreg)
                  idas_itreg(itreg(m):itreg(m)+mtrn-1)
     &                      = idas_itreg_temporary(1:mtrn)
                  itrgn(m) = ntrn
                  itrgm(m) = mtrn

                  itsmn(m) = itsmn(m) + ( mtrn + mod(mtrn,2) ) / 2


                  italm(m) = mmmax

                  itals(m) = 46

*-----------------------------------------------------------------------
*        t-wwbg tally
*-----------------------------------------------------------------------

         else if( ital(m) .eq. 22 ) then

*-----------------------------------------------------------------------
*              check region
*-----------------------------------------------------------------------

            if( itmsh(m) .eq. 1 ) then

!<-20220407murofushi update
!--                  igm  = ( mmmax - 1 ) * 2 + 1
                  igm  = 1

                  itrcm(m) = itrgm(m)

                  ntrn = itrgn(m)
                  mtrn = itrgm(m)
                  ign  = itrcg(m)
                  call tregion2(io,jo,ierr,ntrn,mtrn,idas_itrcg(ign)
     &                          ,igm,MAX_NUM_ITREG,idas_itreg_temporary)

                  if( mtrn > MAX_NUM_ITREG ) then
                     write(ErrCha,'(a,a,i5,a,a,i5,a)')
     &                    'sub.settal t-wwbg@tallsm1.f'
     &                       //' ?dimension over idas_itreg_temporary?'
     &                       //' mtrn > MAX_NUM_ITREG'
     &                    ,' (mtrn=',mtrn,')'
     &                    ,' (MAX_NUM_ITREG@moddas.f=',MAX_NUM_ITREG,')'
                     ErrID = 'L:3987/R:settal/F:tallsm1.f'
                     call ErrWrite(ErrID,ErrCha)
                  endif

                     if( ierr .ne. 0 ) then

                        write(io,'('' **** above error in '',i2,
     &                  ''-th tally'')') m

                        ErrCha = ''
                        MsgID = 'L:3997/R:settal/F:tallsm1.f'
                        call ErrWrite(MsgID, ErrCha)
                        write(jo,'('' **** above error in '',i2,
     &                  ''-th tally'')') m

                        goto 999

                     end if

                  call moddas_reallocate_int(
     &                    itlmax, m, mtrn, itreg, idas_itreg)
                  idas_itreg(itreg(m):itreg(m)+mtrn-1)
     &                      = idas_itreg_temporary(1:mtrn)
                  itrgn(m) = ntrn
                  itrgm(m) = mtrn

                  itsmn(m) = itsmn(m) + ( mtrn + mod(mtrn,2) ) / 2


            else if( itmsh(m) .eq. 4 ) then

                  itrcm(m) = itrgm(m)

                  ntrn = itrgn(m)
                  mtrn = itrgm(m)

                  ign  = itrcg(m)

                  call moddas_allocate_int(nelemtot,ielem2ir)

                  call ttetmesh2(io,jo,ierr,ntrn,mtrn,idas_itrcg(ign),
     &                 MAX_NUM_ITREG,idas_itreg_temporary,ielem2ir)

                  if( mtrn > MAX_NUM_ITREG ) then
                     write(ErrCha,'(a,a,i5,a,a,i5,a)')
     &                    'sub.settal t-track@tallsm1.f'
     &                       //' ?dimension over idas_itreg_temporary?'
     &                       //' mtrn > MAX_NUM_ITREG'
     &                    ,' (mtrn=',mtrn,')'
     &                    ,' (MAX_NUM_ITREG@moddas.f=',MAX_NUM_ITREG,')'
                     ErrID = 'L:4037/R:settal/F:tallsm1.f'
                     call ErrWrite(ErrID,ErrCha)
                  endif

                     if( ierr .ne. 0 ) then

                        write(io,'('' **** above error in '',i2,
     &                  ''-th tally'')') m

                        ErrCha = ''
                        MsgID = 'L:4047/R:settal/F:tallsm1.f'
                        call ErrWrite(MsgID, ErrCha)
                        write(jo,'('' **** above error in '',i2,
     &                  ''-th tally'')') m

                        goto 999

                     end if

                  if(ntrn.gt.0)then

                   call moddas_reallocate_int(
     &                  itlmax, m, mtrn+ntrn,
     &                  itreg, idas_itreg)
                   idas_itreg(itreg(m):itreg(m)+mtrn-1)
     &                  = idas_itreg_temporary(1:mtrn)
                   idas_itreg(itreg(m)+mtrn:itreg(m)+mtrn+ntrn-1)
     &                  = ielem2ir(1:ntrn)

                   itrgn(m) = ntrn
                   itrgm(m) = mtrn+ntrn

                  else

                   ntrn0=idas_itreg_temporary(4)
                   idas_itreg_temporary(4)=0

                   call moddas_reallocate_int(
     &                  itlmax, m, mtrn, itreg, idas_itreg)
                   idas_itreg(itreg(m):itreg(m)+mtrn-1)
     &                  = idas_itreg_temporary(1:mtrn)
                   itrgn(m) = ntrn0
                   itrgm(m) = mtrn

                  endif

                  call moddas_deallocate_int(ielem2ir)

            end if

*-----------------------------------------------------------------------
*           specify the subroutine number and storage position
*-----------------------------------------------------------------------

                  italm(m) = mmmax

            if( itmsh(m) .eq. 1 ) then

                  itals(m) = 47

            else if( itmsh(m) .eq. 3 ) then

                  itals(m) = 56

            else if( itmsh(m) .eq. 4 ) then

                  itals(m) = 57

            end if

cKN 2017/01/03

*-----------------------------------------------------------------------
*        t-userdefined tally
*-----------------------------------------------------------------------

         else if( ital(m) .eq. 20 ) then

            iusrtally = 1

*-----------------------------------------------------------------------
*           open file for writing outputs
*-----------------------------------------------------------------------

                  infil = itfln(m)

                  if( infil .gt. 50 ) then

                        write(io,'('' **** Error in '',i2,
     &                  ''-th tally : number of file for user defined''
     &                  '' tally exceeds 50 files'')') m

                        ErrCha = ''
                        MsgID = 'L:4130/R:settal/F:tallsm1.f'
                        call ErrWrite(MsgID, ErrCha)
                        write(jo,'('' **** Error in '',i2,
     &                  ''-th tally : number of file for user defined''
     &                  '' tally exceeds 50 files'')') m

                        goto 999

                  end if

            do i = 1, infil

                  iudtf(i) = i + 150

               if( npe .le. 1 ) then

                     open(iudtf(i), file = cudtfln(m,i),
     &                    form='formatted', status = 'unknown' )

               else

                  if( me .eq. 0 ) then

                     open(iudtf(i), file = cudtfln(m,i),
     &                    form='formatted', status = 'unknown' )

                  else if( me .gt. 0 ) then

                     iorder = aint(log10(real(npe)))+1
                     if ( iorder .lt. 3) iorder = 3

                     write(chme,'(i5.5)') me

                     filnm = cudtfln(m,i)(1:iudtfll(m,i))
     &                    //'.' // chme(6-iorder:5)

                     open(iudtf(i), file = filnm,
     &                    form='formatted', status = 'unknown' )

                  end if

               end if

            end do

*-----------------------------------------------------------------------

         end if

*-----------------------------------------------------------------------

      end do

               m = itnm
               kmmax = mmmax + 1

               mmmax = kmmax + jmmax + lmmax

               if( mmmax .gt. mdas ) goto 950

               itlmx = jmmax + lmmax
               lmmax = kmmax + jmmax
               nmmax = kmmax

      end if

      goto 998

*-----------------------------------------------------------------------

  940    write(io,'(/
     &         ''<<< ERROR : at the '',i2,''-th tally''/
     &         ''    rshow is not zero but # of reg is 1'')')
     &                m

         ErrCha = ''
         MsgID = 'L:4206/R:settal/F:tallsm1.f'
         call ErrWrite(MsgID, ErrCha)
         write(jo,'(/
     &         ''<<< ERROR : at the '',i2,''-th tally''/
     &         ''    rshow is not zero but # of reg is 1'')')
     &                m

         goto 999

*-----------------------------------------------------------------------

  950    write(io,'(/
     &         ''<<< Memory ERROR : at the '',i2,''-th tally'',
     &         '' in '',i2,'' total tally >>>''/
     &         ''*  memory exceeds mdas'',/
     &         ''*  start of tally memory ='',i9,/
     &         ''*    end of tally memory ='',i9,/
     &         ''*     total tally memory ='',i9,/
     &         ''*     total memory: mdas ='',i9,/
     &         ''<<<  Please extend mdas in param.inc >>>'')')
     &                m, itnm, ntstar, mmmax, mmmax-ntstar, mdas

         ErrCha = ''
         MsgID = 'L:4229/R:settal/F:tallsm1.f'
         call ErrWrite(MsgID, ErrCha)
         write(jo,'(/
     &         ''<<< Memory ERROR : at the '',i2,''-th tally'',
     &         '' in '',i2,'' total tally >>>''/
     &         ''*  memory exceeds mdas'',/
     &         ''*  start of tally memory ='',i9,/
     &         ''*    end of tally memory ='',i9,/
     &         ''*     total tally memory ='',i9,/
     &         ''*     total memory: mdas ='',i9,/
     &         ''<<<  Please extend mdas in param.inc >>>'')')
     &                m, itnm, ntstar, mmmax, mmmax-ntstar, mdas

         goto 999

*-----------------------------------------------------------------------

  999 continue

         ierr = 1

  998 continue

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      subroutine tallsave(ic,tr,br,nr,ncol)
*                                                                      *
*        temporally save the tally data                                *
*                                                                      *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      parameter ( io = 15 )

*-----------------------------------------------------------------------

      common /talsav/ iptall
      common /talout/ itall ! T.Sato 2017/08/05
      common /imagecom/ imageout ! T.Sato 2017/08/05
      common /tcntl/  icntl, inucr ! T.Sato 2017/09/27

      dimension tr(*), br(*)

      if(itall.eq.0.and.ncol.ne.2.and.ic.eq.0) then ! output image (eps,vtk,bmp) or not
       imageout=1 ! do not output image
      else
       imageout=0
      endif
      if(icntl.eq.1.or.icntl.eq.13.or.icntl.eq.17) imageout=0 ! for icntl = 1, 13, or 17 always output eps, 2017/09/27
*-----------------------------------------------------------------------

         if( ic .eq. 0 ) then

            if( iptall .eq. 0 ) then

               do i = 1, nr
                  br(i) = tr(i)
               end do

            else

                  open(io,form='unformatted',status='scratch')
               do i = 1, nr
                  write(io) tr(i)
               end do

            end if

         else

            if( iptall .eq. 0 ) then

               do i = 1, nr
                  tr(i) = br(i)
               end do

            else

                  rewind io
               do i = 1, nr
                  read(io) tr(i)
               end do
                  close(io)

            end if

         end if

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      subroutine tallsave_sumover(ic,br,ncol,m)
*                                                                      *
*        temporally save the tally data                                *
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
      call tallsave_sumover_sub(ic,br,ncol,m,tr_sum,nr)

      return
      end

************************************************************************
*                                                                      *
      subroutine tallsave_sumover_sub(ic,br,ncol,m,tr_sum,nr)
*                                                                      *
************************************************************************

      use TALMOD
!$      use TALMOD0

      implicit double precision (a-h,o-z)

      parameter ( io = 150 )

*-----------------------------------------------------------------------

      common /talsav/ iptall
      common /talout/ itall ! T.Sato 2017/08/05
      common /imagecom/ imageout ! T.Sato 2017/08/05
      common /tcntl/  icntl, inucr ! T.Sato 2017/09/27

      dimension br(*)

      real(8) :: tr_sum(nr)

      if(itall.eq.0.and.ncol.ne.2.and.ic.eq.0) then ! output image (eps,vtk,bmp) or not
       imageout=1 ! do not output image
      else
       imageout=0
      endif
      if(icntl.eq.1.or.icntl.eq.13.or.icntl.eq.17) imageout=0 ! for icntl = 1, 13, or 17 always output eps, 2017/09/27
*-----------------------------------------------------------------------

         if( ic .eq. 0 ) then

            if( iptall .eq. 0 ) then

               do i = 1, nr
                  br(i) = tr_sum(i)
               end do

            else

                  open(io,form='unformatted',status='scratch')
               do i = 1, nr
                  write(io) tr_sum(i)
               end do

            end if

         else

            if( iptall .eq. 0 ) then

               do i = 1, nr
                  tr_sum(i) = br(i)
               end do

            else

                  rewind io
               do i = 1, nr
                  read(io) tr_sum(i)
               end do
                  close(io)

            end if

         end if

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      subroutine rdpname(ic,icl,chlw,istyp,inkf0,jstyp,jnkf0,ierr)
*                                                                      *
*       read particle name and nucleus                                 *
*       last modified by K.Niita on 2005/11/30                         *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      include 'param01.inc'

*-----------------------------------------------------------------------

      common /ptname/ pname(20), ipln(20)
      character       pname*8

      character chlw*200
      character element(104)*3,cnuc*3

      data element/
     & 'h  ','he ','li ','be ','b  ','c  ','n  ','o  ',
     & 'f  ','ne ','na ','mg ','al ','si ','p  ','s  ',
     & 'cl ','ar ','k  ','ca ','sc ','ti ','v  ','cr ',
     & 'mn ','fe ','co ','ni ','cu ','zn ','ga ','ge ',
     & 'as ','se ','br ','kr ','rb ','sr ','y  ','zr ',
     & 'nb ','mo ','tc ','ru ','rh ','pd ','ag ','cd ',
     & 'in ','sn ','sb ','te ','i  ','xe ','cs ','ba ',
     & 'la ','ce ','pr ','nd ','pm ','sm ','eu ','gd ',
     & 'tb ','dy ','ho ','er ','tm ','yb ','lu ','hf ',
     & 'ta ','w  ','re ','os ','ir ','pt ','au ','hg ',
     & 'tl ','pb ','bi ','po ','at ','rn ','fr ','ra ',
     & 'ac ','th ','pa ','u  ','np ','pu ','am ','cm ',
     & 'bk ','cf ','es ','fm ','md ','no ','lr ','ku '/

      logical deqn4
      logical dcom2

      dimension jstyp(mxpart), jnkf0(mxpart) ! frtati 2021/10/05 6 -> mxpart

      common /subtra/ isubt, ipsub(mxpart)   ! kitamura22/03/31
      dimension jsubt(mxpart)                ! kitamura22/03/31

*-----------------------------------------------------------------------

            ierr = 0
            ipar = 0
            ipnm = 0
            ipsu = 0  ! kitamura22/03/31

*-----------------------------------------------------------------------

  100 continue

               if( ipar .eq. 0 .and. ipnm .eq. 1 ) then

                  istyp = jstyp(1)
                  inkf0 = jnkf0(1)
                  isubt = jsubt(1)  ! kitamura22/03/31

                  return

               else if(( ipar .eq. 1 .and. ipnm .eq. 0 ) .or.   ! kitamura22/03/31
     &                 ( ipsu .eq. 1 .and. ipnm .eq. 0 )) then  ! kitamura22/03/31

                  goto 998

              else if( chlw(ic:ic) .eq. '-' .and.             ! kitamura22/03/31
     &                  chlw(ic+1:ic+1) .ge. 'a' .and.        ! kitamura22/03/31
     &                  chlw(ic+1:ic+1) .le. 'z' ) then       ! kitamura22/03/31

                  if( ipar .ne. 0 ) goto 998                  ! kitamura22/03/31

                  ipsu = ipsu + 1                             ! kitamura22/03/31
                  if( ipsu .gt. 2 ) goto 998                  ! kitamura22/03/31

                  ic = ic + 1                                 ! kitamura22/03/31
                  ic = jnumc(chlw,ic,icl)                     ! kitamura22/03/31
                  if( ic .gt. icl ) goto 998                  ! kitamura22/03/31

               else if( chlw(ic:ic+1) .eq. '-(' ) then        ! kitamura22/03/31

                  ipar = ipar + 1                             ! kitamura22/03/31
                  ipsu = ipsu + 2                             ! kitamura22/03/31
                  if( ipar .ne. 1 ) goto 998                  ! kitamura22/03/31
                  if( ipsu .ne. 2 ) goto 998                  ! kitamura22/03/31

                  ic = ic + 2                                 ! kitamura22/03/31
                  ic = jnumc(chlw,ic,icl)                     ! kitamura22/03/31
                  if( ic .gt. icl ) goto 998                  ! kitamura22/03/31

               else if( chlw(ic:ic) .eq. '(' ) then

                  ipar = ipar + 1
                  if( ipar .ne. 1 ) goto 998

                  ic = ic + 1
                  ic = jnumc(chlw,ic,icl)
                  if( ic .gt. icl ) goto 998

               else if( chlw(ic:ic) .eq. ')' ) then

                  ipar = ipar - 1
                  if( ipsu .eq. 2 ) ipsu = ipsu - 1           ! kitamura22/03/31
                  if( ipar .ne. 0 ) goto 998
                  if( ipnm .eq. 0 ) goto 998

                  ic = ic + 1
                  ic = jnumc(chlw,ic,icl)

                  if( ipnm .gt. 1 ) then                      ! kitamura22/03/31
                   istyp = -ipnm
                   inkf0 = 0
                   isubt = ipsu                               ! kitamura22/03/31
                  else if( ipnm .eq. 1 ) then                 ! kitamura22/03/31
                   istyp = jstyp(ipnm)                        ! kitamura22/03/31
                   inkf0 = jnkf0(1)                           ! kitamura22/03/31
                   isubt = ipsu                               ! kitamura22/03/31
                  end if                                      ! kitamura22/03/31

                  return

               end if

*-----------------------------------------------------------------------
*           check pname
*-----------------------------------------------------------------------

            do i = 1, 20

               if( chlw(ic:ic+ipln(i)-1) .eq. pname(i)(1:ipln(i)) ) then

                  ic = ic + ipln(i)
                  ic = jnumc(chlw,ic,icl)

                  ipnm = ipnm + 1

                  jstyp(ipnm) = i
                  jnkf0(ipnm) = kfft( jstyp(ipnm) )

                  if( ipsu .eq. 1) then                       ! kitamura22/03/31
                     jsubt(ipnm) = ipsu                       ! kitamura22/03/31
                  else                                        ! kitamura22/03/31
                     jsubt(ipnm) = 0                          ! kitamura22/03/31
                  end if                                      ! kitamura22/03/31

                  goto 100

               end if

            end do

*-----------------------------------------------------------------------
*        check nucleus or kf code
*-----------------------------------------------------------------------

                  ipnm = ipnm + 1

                        isa = 0
                        isn = 0

                        ica = 0
                        icb = 0
                        icm = 0
                        icn = 0
                        ism = 0 ! 2025/03/10 Ogawa. Isomer flag. Accept Na24m. kfcode or 24Nam are not accepted.

                     do i = ic, icl

                        if( i .eq. icl .and. isa .gt. 0 .and. isn .gt. 0
     &                    .and. chlw(i:i) .eq. 'm' ) then
                           ism = 1
                        elseif(i .eq. icl .and. isa .gt. 0 .and.
     &                   isn .gt. 0 .and. chlw(i:i) .eq. 'n') then
                           ism = 2
                        else if( chlw(i:i) .ge. 'a' .and.
     &                      chlw(i:i) .le. 'z' ) then
                           isa = isa + 1
                           if( isa .eq. 1 ) ica = i
                           icb = i
                        else if( deqn4( chlw(i:i) ) ) then
                           isn = isn + 1
                           if( isn .eq. 1 ) icm = i
                           icn = i
                        else if( chlw(i:i) .eq. ')' ) then
                           icd = i-1
                           goto 502
                        else if( dcom2( chlw(i:i) ) .or.
     &                           i .eq. icl ) then
                           icd = i
                           goto 502
                        else
                           goto 994
                        end if

                     end do

                        icd = icl

  502                continue

*-----------------------------------------------------------------------

               if( ica .gt. 0 ) then

                    if( icb-ica .lt. 0 .or. icb-ica .gt. 1 ) goto 994

                     cnuc = chlw(ica:icb)//'  '

                     do j = 1, 104
                        if( cnuc(1:3) .eq. element(j)(1:3) ) then
                           icha = j
                           goto 452
                        end if
                     end do

                           goto 994

  452                continue

                     if( icha .gt. 104 )  goto 994

                  if( icm .eq. 0 .or. icn .eq. 0 ) then

                     jnkf0(ipnm) = - icha * 1000000
                     jstyp(ipnm) = 19
                     jsubt(ipnm) = ipsu                       ! kitamura22/03/31

                  else

                     if( isn .gt. 3 ) goto 994

                     read(chlw(icm:icn),'(i5)') masi

                     if( ica .gt. icm .and. masi .lt. 0 ) then !FURUTA20240813
                        masi = -masi                          ! kitamura22/03/31
                        ipsu = ipsu + 1                       ! kitamura22/03/31
                     elseif( masi .lt. 0 ) then               ! kitamura22/03/31
                        masi = -masi                          ! kitamura22/03/31
                     end if                                   ! kitamura22/03/31

                     if( masi .lt. icha ) goto 994
                     if( masi-icha .gt. maxnt ) goto 994

                     jnkf0(ipnm) = icha * 1000000 + masi + 1000000000 *
     &                ism
                     jstyp(ipnm) = 19
                     jsubt(ipnm) = ipsu                       ! kitamura22/03/31

                  end if

                     ic = icd + 1

               else

                  if( icm .eq. 0 .or. icn .eq. 0 ) goto 994

                     call snum(chlw,icm,icn,ic2,cvvv,ierr) !FURUTA20240813

                     if( ierr .ne. 0 ) goto 998

                     jnkf0(ipnm) = nint( cvvv )
                     jstyp(ipnm) = kftp( jnkf0(ipnm) )
                     jsubt(ipnm) = ipsu                       ! kitamura22/03/31

                     if(chlw(icn+1:icn+1).eq.')')then !FURUTA20240813
                      ic = icn+1
                     else
                      ic = ic2
                     endif

               end if

                     ic = jnumc(chlw,ic,icl)

*-----------------------------------------------------------------------
*           1H --> proton
*-----------------------------------------------------------------------

               if( jnkf0(ipnm) .eq. 1000001 ) then

                     jnkf0(ipnm) = 2212
                     jstyp(ipnm) = 1

*-----------------------------------------------------------------------
*           2H --> deuteron
*-----------------------------------------------------------------------

               else if( jnkf0(ipnm) .eq. 1000002 ) then

                     jstyp(ipnm) = 15

*-----------------------------------------------------------------------
*           3H --> triton
*-----------------------------------------------------------------------

               else if( jnkf0(ipnm) .eq. 1000003 ) then

                     jstyp(ipnm) = 16

*-----------------------------------------------------------------------
*           3He --> 3He (ityp: 19 -> 17)
*-----------------------------------------------------------------------

               else if( jnkf0(ipnm) .eq. 2000003 ) then

                     jstyp(ipnm) = 17

*-----------------------------------------------------------------------
*           4He --> alpha
*-----------------------------------------------------------------------

               else if( jnkf0(ipnm) .eq. 2000004 ) then

                     jstyp(ipnm) = 18

               end if

*-----------------------------------------------------------------------

         goto 100

*-----------------------------------------------------------------------

  994    continue
         ierr = 994
         return

  998    continue
         ierr = 998
         return

      end


************************************************************************
*                                                                      *
      subroutine pchmlt(m,ityp,ktyp,jtyp,ipn,ips,
     &                  itpan,itpat,jtpat,itldd)
*                                                                      *
*       check of particle in Multiplier                                *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      include 'param.inc'

*-----------------------------------------------------------------------

      dimension itpan(itldd), itpat(itldd,6,2),
     &          jtpat(itldd,6,6,2)

      dimension ips(6)

*-----------------------------------------------------------------------
*        check of particles
*-----------------------------------------------------------------------

                     ipn = 0

         do i = 1, itpan(m)

                     ipp = 0

*-----------------------------------------------------------------------

            if( itpat(m,i,1) .gt. 0 ) then

*-----------------------------------------------------------------------

                  if( itpat(m,i,1) .eq. 20 ) then

                        ipp = 1

                  else if( itpat(m,i,1) .ne. 19 ) then

                        if( ktyp .eq. itpat(m,i,2) ) ipp = 1

                  else if( itpat(m,i,1) .eq. 19 .and.
     &                     itpat(m,i,2) .eq. -2000000 .and.
     &                   ( ityp .eq. 17 .or. ityp .eq. 18 ) ) then

                        ipp = 1

                  else if( itpat(m,i,1) .eq. 19 .and.
     &                     itpat(m,i,2) .eq. -1000000 .and.
     &                   ( ityp .eq. 1 .or. ityp .eq. 15 .or.
     &                     ityp .eq. 16 ) ) then

                        ipp = 1

                  else if( itpat(m,i,1) .eq. 19 .and.
     &                     ityp .ge. 15 ) then

                     if( itpat(m,i,2) .eq. 0 .and.
     &                   ityp .eq. 19 ) then

                        ipp = 1

                     else if( itpat(m,i,2) .gt. 0 ) then

                        if( ktyp .eq. itpat(m,i,2) ) ipp = 1

                     else if( itpat(m,i,2) .lt. 0 ) then

                        if( jtyp .eq. -itpat(m,i,2)/1000000 )
     &                      ipp = 1

                     end if

                  end if

*-----------------------------------------------------------------------

            else if( itpat(m,i,1) .lt. 0 ) then

*-----------------------------------------------------------------------

               do k = 1, -itpat(m,i,1)

                  if( jtpat(m,i,k,1) .eq. 20 ) then

                        ipp = 1

                  else if( jtpat(m,i,k,1) .ne. 19 ) then

                        if( ktyp .eq. jtpat(m,i,k,2) ) ipp = 1

                  else if( itpat(m,i,1) .eq. 19 .and.
     &                     itpat(m,i,2) .eq. -2000000 .and.
     &                   ( ityp .eq. 17 .or. ityp .eq. 18 ) ) then

                        ipp = 1

                  else if( itpat(m,i,1) .eq. 19 .and.
     &                     itpat(m,i,2) .eq. -1000000 .and.
     &                   ( ityp .eq. 1 .or. ityp .eq. 15 .or.
     &                     ityp .eq. 16 ) ) then

                        ipp = 1

                  else if( jtpat(m,i,k,1) .eq. 19 .and.
     &                     ityp .ge. 15 ) then

                     if( jtpat(m,i,k,2) .eq. 0 .and.
     &                   ityp .eq. 19 ) then

                        ipp = 1

                     else if( jtpat(m,i,k,2) .gt. 0 ) then

                        if( ktyp .eq. jtpat(m,i,k,2) ) ipp = 1

                     else if( jtpat(m,i,k,2) .lt. 0 ) then

                        if( jtyp .eq. -jtpat(m,i,k,2)/1000000 )
     &                      ipp = 1

                     end if

                  end if

                  if( ipp .ne. 0 ) goto 100

               end do

  100          continue

*-----------------------------------------------------------------------

            end if

*-----------------------------------------------------------------------

               if( ipp .ne. 0 ) then

                  ipn = ipn + 1
                  ips(ipn) = i

               end if

*-----------------------------------------------------------------------

         end do

      return
      end


************************************************************************
*                                                                      *
      subroutine pcheck(m,np,ityp,ktyp,jtyp,ipn,ips)
*                                                                      *
*       check of particle in itpat                                     *
*                                                                      *
************************************************************************
      use partmod, only: itmxpt, itpan, itpat, jtpat ! frtati 2021/10/05
*-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)

      include 'param.inc'
      include 'param01.inc'

*-----------------------------------------------------------------------


      dimension ips(itmxpt) ! frtati 2021/10/05 6 -> itmxpt

*-----------------------------------------------------------------------
*        check of particles
*-----------------------------------------------------------------------

                     ipn = 0

         do i = 1, np

                     ipp = 0

*-----------------------------------------------------------------------

            if( itpat(m,i,1) .gt. 0 ) then

*-----------------------------------------------------------------------

                  if( itpat(m,i,1) .eq. 20 ) then

                        ipp = 1

                  else if( itpat(m,i,1) .ne. 19 ) then

                        if( ktyp .eq. itpat(m,i,2) ) ipp = 1

                  else if( itpat(m,i,1) .eq. 19 .and.
     &                     itpat(m,i,2) .eq. -2000000 .and.
     &                   ( ityp .eq. 17 .or. ityp .eq. 18 ) ) then

                        ipp = 1

                  else if( itpat(m,i,1) .eq. 19 .and.
     &                     itpat(m,i,2) .eq. -1000000 .and.
     &                   ( ityp .eq. 1 .or. ityp .eq. 15 .or.
     &                     ityp .eq. 16 ) ) then

                        ipp = 1

                  else if( itpat(m,i,1) .eq. 19 .and.
     &                     ityp .ge. 15 ) then

                     if( itpat(m,i,2) .eq. 0 .and.
     &                   ityp .eq. 19 ) then

                        ipp = 1

                     else if( itpat(m,i,2) .gt. 0 ) then

                        if( ktyp .eq. itpat(m,i,2) ) ipp = 1

                     else if( itpat(m,i,2) .lt. 0 ) then

                        if( jtyp .eq. -itpat(m,i,2)/1000000 )
     &                      ipp = 1

                     end if

                  end if

*-----------------------------------------------------------------------

            else if( itpat(m,i,1) .lt. 0 ) then

*-----------------------------------------------------------------------

               do k = 1, -itpat(m,i,1)

                  if( jtpat(m,i,k,1) .eq. 20 ) then

                        ipp = 1

                  else if( jtpat(m,i,k,1) .ne. 19 ) then

                        if( ktyp .eq. jtpat(m,i,k,2) ) ipp = 1

                  else if( itpat(m,i,1) .eq. 19 .and.
     &                     itpat(m,i,2) .eq. -2000000 .and.
     &                   ( ityp .eq. 17 .or. ityp .eq. 18 ) ) then

                        ipp = 1

                  else if( itpat(m,i,1) .eq. 19 .and.
     &                     itpat(m,i,2) .eq. -1000000 .and.
     &                   ( ityp .eq. 1 .or. ityp .eq. 15 .or.
     &                     ityp .eq. 16 ) ) then

                        ipp = 1

                  else if( jtpat(m,i,k,1) .eq. 19 .and.
     &                     ityp .ge. 15 ) then

                     if( jtpat(m,i,k,2) .eq. 0 .and.
     &                   ityp .eq. 19 ) then

                        ipp = 1

                     else if( jtpat(m,i,k,2) .gt. 0 ) then

                        if( ktyp .eq. jtpat(m,i,k,2) ) ipp = 1

                     else if( jtpat(m,i,k,2) .lt. 0 ) then

                        if( jtyp .eq. -jtpat(m,i,k,2)/1000000 )
     &                      ipp = 1

                     end if

                  end if

                  if( ipp .ne. 0 ) goto 100

               end do

  100          continue

*-----------------------------------------------------------------------

            end if

*-----------------------------------------------------------------------

               if( ipp .ne. itpat(m,i,3) ) then  ! kitamura22/03/31

                  ipn = ipn + 1
                  ips(ipn) = i

               end if

*-----------------------------------------------------------------------

         end do

      return
      end


************************************************************************
*                                                                      *
      subroutine pcchck(m,ityp,ktyp,jtyp,icpan,icpat,icc)
*                                                                      *
*       check of particle in counter                                   *
*       last modified by K.Niita on 2006/01/16                         *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      include 'param.inc'

*-----------------------------------------------------------------------

      dimension icpan(3), icpat(3,20,2)

*-----------------------------------------------------------------------
*        check of particles
*-----------------------------------------------------------------------

                  icc = 0
                  if( icpan(m) .gt. 0 ) icc = 1

            do i = 1, abs( icpan(m) )

                     ipp = 0

                  if( icpat(m,i,1) .eq. 20 ) then

                        ipp = 1

                  else if( icpat(m,i,1) .ne. 19 ) then

                        if( ktyp .eq. icpat(m,i,2) ) ipp = 1

                  else if( icpat(m,i,1) .eq. 19 .and.
     &                     ityp .ge. 15 ) then

                     if( icpat(m,i,2) .eq. 0 .and.
     &                   ityp .eq. 19 ) then

                        ipp = 1

                     else if( icpat(m,i,2) .gt. 0 ) then

                        if( ktyp .eq. icpat(m,i,2) ) ipp = 1

                     else if( icpat(m,i,2) .lt. 0 ) then

                        if( jtyp .eq. -icpat(m,i,2)/1000000 )
     &                      ipp = 1

                     end if

                  else if( icpat(m,i,1) .eq. 19 .and.
     &                     icpat(m,i,2) .eq. -1000000 .and.
     &                     ityp .eq. 1 ) then

                        ipp = 1

                  end if

                  if( ipp .eq. 1 ) return

            end do

                  icc = 1
                  if( icpan(m) .gt. 0 ) icc = 0

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine tregck(kreg,ilev,ilat,mr,kr,jj,icc)
*                                                                      *
*       region check in tally                                          *
*       last modified by K.Niita on 2011/02/03                         *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      dimension kr(mr)

      dimension ilat(5,10)
      dimension lats(6)

      dimension ipar(0:20)
      dimension jpar(0:20)
      dimension klev(0:20)
      dimension iccp(0:20)

*-----------------------------------------------------------------------

                     icc = 0

                     jj = jj + 1
                     nreg = kr(jj)

*-----------------------------------------------------------------------
*              single entry or ( all )
*-----------------------------------------------------------------------

                  if( nreg .eq. 6000000 ) then

                     icc = 1
                     return

                  else if( nreg .gt. 0 .and. nreg .lt. 1000000 ) then

                     if( nreg .eq. kreg ) icc = 1
                     return

                  end if

*-----------------------------------------------------------------------
*           combine, lattice or level structure entries
*-----------------------------------------------------------------------

                     kpar = 0
                     jlev = 0

                  do i = 0, 20

                     ipar(i) = 0
                     jpar(i) = 0
                     klev(i) = 0
                     iccp(i) = 0

                  end do

*-----------------------------------------------------------------------

  600       continue

               if( nreg .gt. 0 .and. nreg .lt. 1000000 ) then

                        ipar(kpar) = ipar(kpar) + 1

                     if( jlev .eq. 0 ) then

                        if( nreg .eq. kreg ) then

                              iccp(kpar) = iccp(kpar) + 1

                        end if

                     else if( jlev .le. ilev ) then

                        do l = jlev, ilev

                           if( nreg .eq. ilat(1,l) ) then

                              iccp(kpar) = iccp(kpar) + 1

                           end if

                        end do

                     end if

               else if( nreg .lt. 0 ) then

                        kpar = kpar + 1
                        jpar(kpar) = -nreg

               else if( nreg .gt. 3000000 .and. nreg .lt. 4000000 ) then

                        knum = nreg - 3000000

                        kpar = kpar + 1
                        jpar(kpar) = knum
                        klev(kpar) = 1

               else if( nreg .gt. 4000000 .and. nreg .lt. 5000000 ) then

                        kpar = kpar + 1
                        jpar(kpar) = -1
                        ipar(kpar) = -1
                        klev(kpar) = 0

                        nreg = nreg - 4000000

                        jj = jj + 1
                        jlat = kr(jj)

                  do ll = 1, jlat

                     do m = 1, 6

                           jj = jj + 1
                           lats(m) = kr(jj)

                     end do

                     if( jlev .gt. 0 .and. jlev .le. ilev ) then

                        do l = jlev, ilev

                           if( nreg .eq. ilat(1,l) .and.
     &                         ilat(2,l) .ne. 0 ) then

                              if( ilat(3,l) .ge. lats(1) .and.
     &                            ilat(3,l) .le. lats(2) .and.
     &                            ilat(4,l) .ge. lats(3) .and.
     &                            ilat(4,l) .le. lats(4) .and.
     &                            ilat(5,l) .ge. lats(5) .and.
     &                            ilat(5,l) .le. lats(6) ) then

                                 iccp(kpar) = iccp(kpar) + 1

                              end if

                           end if

                        end do

                     end if

                  end do

               end if

*-----------------------------------------------------------------------

  660          continue

               if( ipar(kpar) .eq. jpar(kpar) ) then

                  if( klev(kpar) .ne. 0 .and.
     &                iccp(kpar) .ne. jpar(kpar) ) iccp(kpar) = 0

                  if( iccp(kpar) .ne. 0 )
     &                iccp(kpar-1) = iccp(kpar-1) + 1

                     ipar(kpar) = 0
                     jpar(kpar) = 0
                     iccp(kpar) = 0

                  if( klev(kpar) .gt. 0 ) then

                     klev(kpar) = 0
                     jlev       = 0

                  end if

                     kpar = kpar - 1
                     ipar(kpar) = ipar(kpar) + 1

                  if( kpar .eq. 0 ) goto 560

                     goto 660

               end if

               if( klev(kpar) .gt. 0 .and. ipar(kpar) .gt. 0 ) then

                  if( iccp(kpar) .eq. 0 ) then

                     jlev = ilev + 1

                  else

                     jlev = jlev + 1

                  end if

               end if

                     jj = jj + 1
                     nreg = kr(jj)

         goto 600

  560    continue

         icc = iccp(0)

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine ttetck(ic,mr,kr,itet,ihelem,ir,icc)
*                                                                      *
*       tetra check in tally                                           *
*       Modified by T.Furuta on 2025/01/16                             *
*                                                                      *
************************************************************************
      use TETRAMOD, only:nelem
      use moddas_ggs !frtati20220905
      implicit real*8 (a-h,o-z)
      include 'param.inc'
      include 'ggsparam.inc'
*-----------------------------------------------------------------------

      integer,intent(in) :: ic,itet
      integer,intent(in) :: ihelem
      integer,intent(out) :: icc
      integer,intent(in) :: mr,kr(mr)
      integer :: jhelem,numelem
      logical :: iflag

      common /regdc/  idrg(kvlmax), idgr(kvmmax)
*-----------------------------------------------------------------------

      icc = 0

      if(lat(llat+1,ic).eq.3)then ! in box but not in Tetra elements
       icc=-1
       return
      endif
      if(lat(llat+1,ic).ne.-3)return ! Not in Tetra regions icc=0
      itet0=kr(2)
      nr=kr(3)
      iflag=.false.
      do ir=1,nr
       if(ic.eq.idgr(kr(3+ir)))then
        iflag=.true.
        exit
       endif
      enddo
      if(nr.gt.0.and..not.iflag)then ! No in REG
       icc=-2
       return
      endif
      numelem=kr(4+nr)
      if(itet.eq.itet0)then
       jhelem=ihelem-nelem(itet-1)
       if(jhelem.gt.0.and.jhelem.le.nelem(itet)-nelem(itet-1))then
        if(numelem.gt.0)then
         ir=kr(4+nr+jhelem)
        else
         ir=jhelem
        endif
        icc=1 ! ir is always >0
       endif
      endif
      if(ir.eq.0)icc=-2

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      subroutine tregvl(nr,mr,kr,vl,lr,nvl,ivl,rvl)
*                                                                      *
*       calculate region volume in tally                               *
*       last modified by K.Niita on 2011/02/03                         *
*                                                                      *
************************************************************************
      use moddas_ggs !frtati20220905

      implicit real*8 (a-h,o-z)

      include 'param.inc'
      include 'ggsparam.inc'

*-----------------------------------------------------------------------

      common /inpec/ ititl, ipara, ibody, iregn, llarr, itby, itar

      common /volreg/ dvol(kvlmax)
      common /regdc/  idrg(kvlmax), idgr(kvmmax)
      common /regdm/  idmg(kvlmax)

      dimension kr(mr)
      dimension vl(nr)
      dimension lr(nr)
      dimension ivl(nvl)
      dimension rvl(nvl)

      dimension ipar(0:20)
      dimension jpar(0:20)
      dimension klev(0:20)
      dimension vccp(0:20)
      dimension jlat(6)

*-----------------------------------------------------------------------

                     jj = 0
                     jc = 0

      do 100 ir = 1, nr

                     jj = jj + 1
                     nreg = kr(jj)

*-----------------------------------------------------------------------
*              single entry or ( all )
*-----------------------------------------------------------------------

               if( nreg .eq. 6000000 ) then

                     jc = jc + 1
                     lr(ir) = jc + 1000000

                     vl(ir) = 0.0

                  do l = 1, iregn

                     if( ( junf .eq. 0 .or. jun(l) .eq. 0 ) .and.
     &                     idmg(l) .ge. 0 )
     &                  vl(ir) = vl(ir) + dvol(l)

                  end do

                     goto 150

               else if( nreg .gt. 0 .and. nreg .lt. 1000000 ) then

                     lr(ir) = nreg
                     vl(ir) = dvol( idgr( nreg ) )

                     goto 150

               end if

*-----------------------------------------------------------------------
*           combine, lattice or level structure entries
*-----------------------------------------------------------------------

                     jc = jc + 1
                     lr(ir) = jc + 1000000

                     kpar = 0
                     jlev = 0

                  do i = 0, 20

                     ipar(i) = 0
                     jpar(i) = 0
                     klev(i) = 0
                     vccp(i) = 0.0

                  end do

*-----------------------------------------------------------------------

  600       continue

               if( nreg .gt. 0 .and. nreg .lt. 1000000 ) then

                        ipar(kpar) = ipar(kpar) + 1

                     if( jlev .eq. 0 ) then

                        vccp(kpar) = vccp(kpar) + dvol( idgr( nreg ) )

                     end if

               else if( nreg .lt. 0 ) then

                        kpar = kpar + 1
                        jpar(kpar) = -nreg

               else if( nreg .gt. 3000000 .and. nreg .lt. 4000000 ) then

                        knum = nreg - 3000000

                        kpar = kpar + 1
                        jpar(kpar) = knum
                        klev(kpar) = 1

               else if( nreg .gt. 4000000 .and. nreg .lt. 5000000 ) then

                        kpar = kpar + 1
                        jpar(kpar) = -1
                        ipar(kpar) = -1
                        klev(kpar) = 0

                        nreg = nreg - 4000000

                        jj = jj + 1
                        ilat = kr(jj)

                        lats = 0

                  do ll = 1, ilat

                     do m = 1, 6

                        jj = jj + 1
                        jlat(m) = kr(jj)

                     end do

                        lats = lats + ( 1 + jlat(2) - jlat(1) )
     &                              * ( 1 + jlat(4) - jlat(3) )
     &                              * ( 1 + jlat(6) - jlat(5) )

                  end do

                     if( jlev .eq. 0 ) then

                        vccp(kpar) = vccp(kpar)
     &                             + dvol( idgr( nreg ) ) * lats

                     end if

               end if

*-----------------------------------------------------------------------

  660          continue

               if( ipar(kpar) .eq. jpar(kpar) ) then

                     vccp(kpar-1) = vccp(kpar-1) + vccp(kpar)

                     ipar(kpar) = 0
                     jpar(kpar) = 0
                     vccp(kpar) = 0.0

                  if( klev(kpar) .gt. 0 ) then

                     klev(kpar) = 0
                     jlev       = 0

                  end if

                     kpar = kpar - 1
                     ipar(kpar) = ipar(kpar) + 1

                  if( kpar .eq. 0 ) goto 560

                     goto 660

               end if

               if( klev(kpar) .gt. 0 .and. ipar(kpar) .gt. 0 ) then

                     jlev = jlev + 1

               end if

                     jj = jj + 1
                     nreg = kr(jj)

         goto 600

  560    continue

         vl(ir) = vccp(0)

*-----------------------------------------------------------------------

  150 continue

         if( nvl .gt. 0 ) then

            do k = 1, nvl

               if( ivl(k) .eq. lr(ir) ) then

                  vl(ir) = rvl(k)
                  goto 100

               end if

            end do

         end if

  100 continue

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine ttetvl(mr,kr,nr,vl,lr)
*                                                                      *
*       calculate tetra volume in tally                                *
*       Modified by T.Furuta on 2025/01/16                             *
*                                                                      *
************************************************************************
      use TETRAMOD, only:nelem,ielem2id,volelems
      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      dimension kr(mr)
      dimension vl(nr),lr(nr)

*-----------------------------------------------------------------------

      itet0=kr(2)
      nr0=kr(3)
      numelem=kr(4+nr0)
      if(nr0.gt.0)then
       do ielem=1,numelem
        ir=kr(4+nr0+ielem)
        if(ir.gt.0)then
         lr(ir)=ielem2id(ielem+nelem(itet0-1))
         vl(ir)=volelems(ielem+nelem(itet0-1))
        endif
       enddo
      else
       do ir = 1,nr
        lr(ir)=ielem2id(ir+nelem(itet0-1))
        vl(ir)=volelems(ir+nelem(itet0-1))
       enddo
      endif

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      subroutine ttetcm(mr,kr,nr,xcm)
*                                                                      *
*       calculate tetra CM coordinates in tally                        *
*       Created by T.Furuta on 2019/10/28                              *
*                                                                      *
************************************************************************
      use TETRAMOD, only:nelem,ielem2point,pointxyz
      implicit real*8 (a-h,o-z)
*-----------------------------------------------------------------------

      dimension kr(mr)
      dimension xcm(3,nr)

*-----------------------------------------------------------------------
      integer ip(4)
      real(8) xsum(3)

      itet0=kr(2)
      nr0=kr(3)
      numelem=kr(4+nr0)
      if(nr0.gt.0)then
       do ielem=1,numelem
        ir=kr(4+nr0+ielem)
        if(ir.gt.0)then
         ip(1:4)=ielem2point(1:4,ir+nelem(itet0-1))
         xsum(1:3)=0.0d0
         do j=1,4
          xsum(1:3)=xsum(1:3)+pointxyz(1:3,ip(j))
         enddo
         xcm(1:3,ir)=0.25d0*xsum(1:3)
        endif
       enddo
      else
       do ir=1,nr
        ip(1:4)=ielem2point(1:4,ir+nelem(itet0-1))
        xsum(1:3)=0.0d0
        do j=1,4
         xsum(1:3)=xsum(1:3)+pointxyz(1:3,ip(j))
        enddo
        xcm(1:3,ir)=0.25d0*xsum(1:3)
       enddo
      endif

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      subroutine remove_spaces(buf)
*                                                                      *
*       remove spaces in buffer characters                             *
*       Created by T.Furuta on 2019/10/28                              *
*                                                                      *
************************************************************************
      implicit none
      character(*),intent(inout) :: buf
      character(len=len(buf)) ctmp
      integer i,j
*-----------------------------------------------------------------------
      j=0
      do i=1,len(buf)
       if(buf(i:i).eq.' ')cycle
       j=j+1
       ctmp(j:j)=buf(i:i)
      enddo
      buf=ctmp(1:j)
*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      subroutine tregval(nr,mr,kr,vl,lr,nvl,ivl,rvl)
*                                                                      *
*       calculate region value in rshow tally                          *
*       last modified by K.Niita on 2011/02/03                         *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      include 'param.inc'
      include 'ggsparam.inc'

*-----------------------------------------------------------------------

      common /inpec/ ititl, ipara, ibody, iregn, llarr, itby, itar

      common /volreg/ dvol(kvlmax)
      common /regdc/  idrg(kvlmax), idgr(kvmmax)
      common /regdm/  idmg(kvlmax)

      dimension kr(mr)
      dimension vl(nr)
      dimension lr(nr)
      dimension ivl(nvl)
      dimension rvl(nvl)

      dimension ipar(0:20)
      dimension jpar(0:20)
      dimension klev(0:20)
      dimension jlat(6)

      data dnon /-1.d+40/

*-----------------------------------------------------------------------

                     jj = 0
                     jc = 0

      do 100 ir = 1, nr

                     jj = jj + 1
                     nreg = kr(jj)

*-----------------------------------------------------------------------
*              single entry,  ( all ) is not available
*-----------------------------------------------------------------------

               if( nreg .gt. 0 .and. nreg .lt. 1000000 ) then

                     lr(ir) = nreg
                     vl(ir) = dnon

                     goto 150

               end if

*-----------------------------------------------------------------------
*           combine, lattice or level structure entries
*-----------------------------------------------------------------------

                     jc = jc + 1
                     lr(ir) = jc + 1000000

                     kpar = 0
                     jlev = 0

                  do i = 0, 20

                     ipar(i) = 0
                     jpar(i) = 0
                     klev(i) = 0

                  end do

*-----------------------------------------------------------------------

  600       continue

               if( nreg .gt. 0 .and. nreg .lt. 1000000 ) then

                        ipar(kpar) = ipar(kpar) + 1

               else if( nreg .lt. 0 ) then

                        kpar = kpar + 1
                        jpar(kpar) = -nreg

               else if( nreg .gt. 3000000 .and. nreg .lt. 4000000 ) then

                        knum = nreg - 3000000

                        kpar = kpar + 1
                        jpar(kpar) = knum
                        klev(kpar) = 1

               else if( nreg .gt. 4000000 .and. nreg .lt. 5000000 ) then

                        kpar = kpar + 1
                        jpar(kpar) = -1
                        ipar(kpar) = -1
                        klev(kpar) = 0

                        nreg = nreg - 4000000

                        jj = jj + 1
                        ilat = kr(jj)

                        lats = 0

                  do ll = 1, ilat

                     do m = 1, 6

                        jj = jj + 1
                        jlat(m) = kr(jj)

                     end do

                        lats = lats + ( 1 + jlat(2) - jlat(1) )
     &                              * ( 1 + jlat(4) - jlat(3) )
     &                              * ( 1 + jlat(6) - jlat(5) )

                  end do

               end if

*-----------------------------------------------------------------------

  660          continue

               if( ipar(kpar) .eq. jpar(kpar) ) then

                     ipar(kpar) = 0
                     jpar(kpar) = 0

                  if( klev(kpar) .gt. 0 ) then

                     klev(kpar) = 0
                     jlev       = 0

                  end if

                     kpar = kpar - 1
                     ipar(kpar) = ipar(kpar) + 1

                  if( kpar .eq. 0 ) goto 560

                     goto 660

               end if

               if( klev(kpar) .gt. 0 .and. ipar(kpar) .gt. 0 ) then

                     jlev = jlev + 1

               end if

                     jj = jj + 1
                     nreg = kr(jj)

         goto 600

  560    continue

         vl(ir) = dnon

*-----------------------------------------------------------------------

  150 continue

         if( nvl .gt. 0 ) then

            do k = 1, nvl

               if( ivl(k) .eq. lr(ir) ) then

                  vl(ir) = rvl(k)
                  goto 100

               end if

            end do

         end if

  100 continue

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine tregion(icc,jsn,jsi,dsin,idsi,ill,ilf,
     &                   jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr,
     &                   ntrn,mtrn,ndsm,nvol,ivl,irvl,icfl
     &                   ,ndim_region,idas_region)
*                                                                      *
*       read 'reg =' sub-section of input tally section                *
*       modified by K.Niita on 2011/02/03                              *
*                                                                      *
************************************************************************
      use moddas
      use moddas_region_mtrg

      implicit real*8 (a-h,o-z)

      include 'param.inc'
      include 'risrcparam.inc' ! T.Sato 2022/03/25
      include 'err.inc'

*-----------------------------------------------------------------------

      character m_err*200
      common /error/ m_err, l_err, k_err

*-----------------------------------------------------------------------

      character chin*200, chlw*200, chcm*200
      character chlc*200

      character dsin(0:9)*200
      dimension idsi(0:9)

      dimension ill(0:9), ilf(0:9)

*-----------------------------------------------------------------------

      dimension ipar(0:20)
      dimension jpar(0:20)
      dimension lpar(0:20)

      data klnmax /1000000/

      integer, intent(in) :: ndim_region
      integer, intent(inout) :: idas_region(ndim_region)

      logical iprojall   ! multi-source subsection, proj=all
*-----------------------------------------------------------------------

            ierr = 0

      if( icfl .eq. 1 ) goto 150

*-----------------------------------------------------------------------
*     read one line from jsi
*-----------------------------------------------------------------------

  140 continue

            call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

               if( ierr .ne. 0 ) return
               if( jpn  .eq. 3 ) return

               if( iskip .ne. 0 ) goto 140

  150 continue

               if( ierr .ne. 0 ) return
               if( jpn  .eq. 3 ) return

ccse multi-source subsection, proj=all, write scrach file
            icl=i1
            chlc = chlw
            call chcomp(chlc,icl,i3,i5)
            if( chlc(icl:icl+7) .ne. '<source>' .and.
     &          chlc(icl:icl) .ne. '[' ) then
               inquire(niws,opened=iprojall)
               if(iprojall) write(niws,'(a)') chin(1:i2)
            endif
*-----------------------------------------------------------------------
*        error not for reg =
*-----------------------------------------------------------------------

         if( icc .eq. 1 ) then

            if( chcm(i1:i1+3) .ne. 'reg=' ) goto 998

         else if( icc .eq. 2 ) then

            if( chcm(i1:i1+8) .ne. 'reginbox=' ) goto 998

         end if

               ic = inumc(chlw,i1,i3,'=') + 1

*-----------------------------------------------------------------------
*        read region number
*-----------------------------------------------------------------------

               ic = jnumc(chlw,ic,i3)

               ntrn = 0
               mtrn = 0
               ibra = 0
               ilat = 0
               klat = 0
               ifis = 0
               kpar = 0
               iuni = 0

               nvol = 0

               do i = 0, 20
                  ipar(i) = 0
                  jpar(i) = 0
                  lpar(i) = 0
               end do

               igm = 0
               call moddas_allocate_int2(0, 20, MAX_NUM_MTRG, mtrg)
               ipmax = 0

            do i = 1, klnmax

                  ic = jnumc(chlw,ic,i3)

                  call tregion3(icn,chlw,i1,i3,ic,
     &                          igm,ipmax,ipar,jpar,
     &                          ibra,ilat,klat,ifis,kpar,iuni
     &                          ,MAX_NUM_MTRG,mtrg)

                     if( icn .gt. 900 ) goto 900
                     if( icn .eq. 500 ) goto 500

               if( i .lt. klnmax .and. ic .gt. i3 ) then

  147             call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

                     if( ierr .ne. 0 ) return
                     if( jpn  .eq. 3 ) goto 500

                     if( iskip .ne. 0 ) goto 147

                     ifis = ifis + 1

                     ic = i1
ccse multi-source subsection, proj=all, write scrach file
            icl=i1
            chlc = chlw
            call chcomp(chlc,icl,i3,i5)
            if( chlc(icl:icl+7) .ne. '<source>' .and.
     &          chlc(icl:icl) .ne. '[' ) then
               inquire(niws,opened=iprojall)
               if(iprojall) write(niws,'(a)') chin(1:i2)
            endif

               end if

            end do

                  goto 996

*-----------------------------------------------------------------------
*        modify the input
*-----------------------------------------------------------------------

  500    continue

               ntrn = ipar(0)
               mtrn = ipar(0)

            do k = 1, mtrn
               idas_region( ndsm + k - 1 ) = mtrg(igm+0,k)
            end do

            call tregion4(ierr,ntrn,mtrn,idas_region(ndsm),ndsm,igm
     &                    ,ndim_region,idas_region)

               if( ierr .ne. 0 ) goto 983

*-----------------------------------------------------------------------


*-----------------------------------------------------------------------
*        volume
*-----------------------------------------------------------------------

            if( chcm(i1:i1+5) .eq. 'volume' .or.
     &          chcm(i1:i1+4) .eq. 'value' ) then

                  call tregvol(jsn,jsi,dsin,idsi,ill,ilf,
     &                         jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr,
     &                         nvol,ivl,irvl)

            end if

*-----------------------------------------------------------------------
*     errors
*-----------------------------------------------------------------------

         return

  900 continue

         if( icn .eq. 980 ) goto 980
         if( icn .eq. 984 ) goto 984
         if( icn .eq. 985 ) goto 985
         if( icn .eq. 986 ) goto 986
         if( icn .eq. 987 ) goto 987
         if( icn .eq. 988 ) goto 988
         if( icn .eq. 989 ) goto 989
         if( icn .eq. 990 ) goto 990
         if( icn .eq. 991 ) goto 991
         if( icn .eq. 992 ) goto 992
         if( icn .eq. 993 ) goto 993
         if( icn .eq. 995 ) goto 995
         if( icn .eq. 997 ) goto 997

*-----------------------------------------------------------------------

  982 continue
         m_err = 'volume of value section is wrong'
         ErrCha = ''
         ErrID = 'L:6209/R:tregion/F:tallsm1.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  983 continue
         m_err = 'region is too many or memory is lack'
         ErrCha = ''
         ErrID = 'L:6220/R:tregion/F:tallsm1.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  980 continue

         m_err = 'n1-n2 should be used like {n1-n2}'
         ErrCha = ''
         ErrID = 'L:6232/R:tregion/F:tallsm1.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  984 continue

         m_err = 'maximum lattice elements by commas is 1000.'
         ErrCha = ''
         ErrID = 'L:6244/R:tregion/F:tallsm1.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  985 continue

         m_err = 'maximum level of < is 10.'
         ErrCha = ''
         ErrID = 'L:6256/R:tregion/F:tallsm1.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  986 continue

         m_err = 'usage of u=# is wrong.'
         ErrCha = ''
         ErrID = 'L:6268/R:tregion/F:tallsm1.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  987 continue

         m_err = '< should be inside ( ).'
         ErrCha = ''
         ErrID = 'L:6280/R:tregion/F:tallsm1.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  988 continue

         m_err = 'maximun level of ( ) is 10.'
         ErrCha = ''
         ErrID = 'L:6292/R:tregion/F:tallsm1.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  989 continue

         m_err = 'usage of latice [i1:i2 i3:i4 i5:i6] is wrong.'
         ErrCha = ''
         ErrID = 'L:6304/R:tregion/F:tallsm1.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  990 continue

         m_err = 'usage of latice [i1 i2 i3] is wrong.'
         ErrCha = ''
         ErrID = 'L:6316/R:tregion/F:tallsm1.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  991 continue

         m_err = 'all or (all) is available. (all<4), (all 3) are not.'
         ErrCha = ''
         ErrID = 'L:6328/R:tregion/F:tallsm1.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  992 continue

         m_err = 'Usage of Parenthesis { n1 - n2 } is wrong'
         ErrCha = ''
         ErrID = 'L:6340/R:tregion/F:tallsm1.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  993 continue

         m_err = 'Usage of Parenthesis ( ) is wrong'
         ErrCha = ''
         ErrID = 'L:6352/R:tregion/F:tallsm1.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  995 continue

         m_err = 'Description of region number '//
     &           '{n1-n2} (n1<n2) is wrong.'
         ErrCha = ''
         ErrID = 'L:6365/R:tregion/F:tallsm1.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  996 continue

         m_err = 'Number of region is too large. max(klnmax)=1000000'
         ErrCha = ''
         ErrID = 'L:6377/R:tregion/F:tallsm1.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  997 continue

         m_err = 'Description of region number is wrong.'
         ErrCha = ''
         ErrID = 'L:6389/R:tregion/F:tallsm1.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  998 continue

         m_err = 'After { mesh = reg } line should be { reg = }.'
         ErrCha = ''
         ErrID = 'L:6401/R:tregion/F:tallsm1.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

      end


************************************************************************
*                                                                      *
      subroutine tregion2(io,jo,ierr,ntrn,mtrn,ntrg,igm_argument
     &                    ,ndim_region,idas_region)
*                                                                      *
*       modify the region mesh                                         *
*       modified by K.Niita on 2011/02/03                              *
*                                                                      *
************************************************************************
      use moddas
      use moddas_region_mtrg
      use moddas_ggs

      implicit real*8 (a-h,o-z)

      include 'param.inc'
      include 'ggsparam.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      common /ccggg/  icgg

      common /tall02/ itrgn(itlmax), itrgm(itlmax), itreg(itlmax)
      common /tall24/ itrcm(itlmax), itrcg(itlmax)


*-----------------------------------------------------------------------

      common /inpec/ ititl, ipara, ibody, iregn, llarr, itby, itar
      common /regdc/ idrg(kvlmax), idgr(kvmmax)
      common /regdm/ idmg(kvlmax)

*-----------------------------------------------------------------------

      dimension ipar(0:20)
      dimension jpar(0:20)
      dimension lpar(0:20)
      dimension klev(0:20)
      dimension nlev(0:20)
      dimension nsav(0:20)
      dimension latj(6)

      dimension ntrg(mtrn)

      integer, intent(in) :: igm_argument
      integer, intent(in) :: ndim_region
      integer, intent(inout) :: idas_region(ndim_region)

      igm = 0
      call moddas_allocate_int2(0, 20, MAX_NUM_MTRG, mtrg)

*-----------------------------------------------------------------------

            ierr = 0

*-----------------------------------------------------------------------


*-----------------------------------------------------------------------

                  ibra = 0
                  kpar = 0
                  mlev = 0
                  ipmax = 0

               do i = 0, 20
                  ipar(i) = 0
                  jpar(i) = 0
                  lpar(i) = 0
                  klev(i) = 0
                  nlev(i) = 0
               end do

         k = 0
  510    k = k + 1

         if( k .gt. mtrn ) goto 520

            if( ntrg(k) .eq. 5000000 ) then

                  do l = 1, iregn

                     if( idmg(l) .ge. 0)then
                      if( junf .eq. 0 .or.
     &                     ( junf .ne. 0 .and. mfl(1,l) .eq. 0 )) then

                        lpar(kpar) = lpar(kpar) + 1
                        ipar(kpar) = ipar(kpar) + 1
                        if( ipar(kpar) .gt. ipmax ) ipmax = ipar(kpar)
                        mtrg(igm+kpar,ipar(kpar)) = idrg(l)

                      end if
                     endif

                  end do

            else if( ntrg(k) .gt. 7000000 ) then

                        if( icgg .eq. 0 ) goto 996

                        muni = ntrg(k) - 7000000

                        iuni = 0

                  do l = 1, mxa

                     if( junf .ne. 0 .and.
     &                   muni .eq. abs(jun(l)) ) then

                        iuni = iuni + 1

                        lpar(kpar) = lpar(kpar) + 1
                        ipar(kpar) = ipar(kpar) + 1
                        if( ipar(kpar) .gt. ipmax ) ipmax = ipar(kpar)
                        mtrg(igm+kpar,ipar(kpar)) = idrg(l)

                        ireg = idrg(l)

                     end if

                  end do

                     if( iuni .eq. 0 ) then

                        write(io,'(/'' **** Error'',
     &                  '' : universe '',i6,'' is not'',
     &                  '' defined in geometry.'')') muni

                        ErrCha = ''
                        MsgID = 'L:6542/R:tregion2/F:tallsm1.f'
                        call ErrWrite(MsgID, ErrCha)
                        write(jo,'(/'' **** Error'',
     &                  '' : universe '',i6,'' is not'',
     &                  '' defined in geometry.'')') muni

                        goto 999

                     end if

            else if( ntrg(k) .eq. 1000000 ) then

                        ibra = ibra + 1

                  if( ibra .eq. 1 ) then

                        lpar(kpar) = lpar(kpar) + 1

                        kpar = kpar + 1
                        ipar(kpar) = 0
                        lpar(kpar) = 0

                  end if

                        k = k + 1

            else if( ntrg(k) .eq. 1000001 ) then

                        ibra = ibra - 1

                  if( ibra .eq. 0 ) then

                        kpar = kpar - 1

                        ipar(kpar) = ipar(kpar) + 1
                        if( ipar(kpar) .gt. ipmax ) ipmax = ipar(kpar)
                        mtrg(igm+kpar,ipar(kpar)) = -lpar(kpar+1)

                        j = 0
  600                   j = j + 1

                     if( j .gt. ipar(kpar+1) ) goto 610

                        ipar(kpar) = ipar(kpar) + 1
                        if( ipar(kpar) .gt. ipmax ) ipmax = ipar(kpar)
                        mtrg(igm+kpar,ipar(kpar)) = mtrg(igm+kpar+1,j)

                        ireg = mtrg(igm+kpar+1,j)

                        if( ireg .gt. 0 .and. ireg .lt. 1000000) then
                            if(idgr(ireg) .eq. 0 ) goto 998
                        endif

                     if( ireg .gt. 4000000 .and.
     &                   ireg .lt. 5000000 ) then

                           if( icgg .eq. 0 ) goto 995

                           ireg = ireg - 4000000

                           if( idgr(ireg) .eq. 0 ) goto 998
                           if( lat(1,idgr(ireg)) .eq. 0 ) goto 997

                           j = j + 1

                           ipar(kpar) = ipar(kpar) + 1
                           if( ipar(kpar) .gt. ipmax )
     &                                         ipmax = ipar(kpar)
                           mtrg(igm+kpar,ipar(kpar)) =
     &                     mtrg(igm+kpar+1,j)

                           nlat = mtrg(igm+kpar+1,j)

                        do l = 1, nlat * 6

                           j = j + 1
                           ipar(kpar) = ipar(kpar) + 1
                           if( ipar(kpar) .gt. ipmax )
     &                                         ipmax = ipar(kpar)
                           mtrg(igm+kpar,ipar(kpar)) =
     &                     mtrg(igm+kpar+1,j)

                        end do

                     end if

                        goto 600

  610                continue

                  end if

            else if( ntrg(k) .eq. 3000000 ) then

                        if( icgg .eq. 0 ) goto 994

                        kpar = kpar + 1
                        ipar(kpar) = 0
                        lpar(kpar) = 0
                        jpar(kpar) = jpar(kpar-1) + 1

                        mlev = mlev + 1

            else if( ntrg(k) .eq. 2000000 ) then

                        if( icgg .eq. 0 ) goto 994

                        kpar = kpar + 1
                        ipar(kpar) = 0
                        lpar(kpar) = 0
                        jpar(kpar) = 1

                        k = k + 1
                        ilev = ntrg(k)
                        k = k + ilev

                        ibrs = ibra
                        ibra = 0

            else if( ntrg(k) .eq. 2000001 ) then

                        nelm = 1
                        kint = kpar - jpar(kpar)

                     do j = kint + 1, kpar

                        nelm = nelm * lpar(j)

                     end do

                        lpar(kint) = lpar(kint) + nelm

                           jlev = jpar(kpar)

                        do j = 1, jlev

                           klev(j) = 1
                           nlev(j) = 0
                           nsav(j) = 1

                        end do

                           ii = 0
                           jj = 1
                           klev(jj) = 0

  100                klev(jj) = klev(jj) + 1

                     if( klev(jj) .gt. lpar(jj+kint) ) then

                           klev(jj) = 1
                           nsav(jj) = 1

                           ii = jj + 1

  400                      continue

                        if( klev(ii) + 1 .le. lpar(ii+kint) ) then

                           klev(ii) = klev(ii) + 1
                           goto 200

                        else

                           if( ii .eq. jlev ) goto 300

                           klev(ii) = 1
                           nsav(ii) = 1
                           ii = ii + 1
                           goto 400

                        end if

                     else

                           ii = jj

                     end if

  200                continue

                           ipar(kint) = ipar(kint) + 1
                           if( ipar(kint) .gt. ipmax )
     &                                         ipmax = ipar(kint)
                           mtrg(igm+kint,ipar(kint)) = jlev + 3000000

               do j = 1, jlev

                        if( j .eq. ii ) then

                           nlev(j) = nlev(j) + 1
                           nsav(j) = nlev(j)

                        else

                           nlev(j) = nsav(j)

                        end if

                           ipar(kint) = ipar(kint) + 1
                           if( ipar(kint) .gt. ipmax )
     &                                         ipmax = ipar(kint)
                           mtrg(igm+kint,ipar(kint)) =
     &                     mtrg(igm+j+kint,nlev(j))
                           mtemp = mtrg(igm+j+kint,nlev(j))

                           ireg = mtemp

                  if( mtemp .le. 0 ) then

                    do l = 1, -mtemp

                           nlev(j) = nlev(j) + 1

                           ipar(kint) = ipar(kint) + 1
                           if( ipar(kint) .gt. ipmax )
     &                                         ipmax = ipar(kint)
                           mtrg(igm+kint,ipar(kint)) =
     &                     mtrg(igm+j+kint,nlev(j))
                           ntemp = mtrg(igm+j+kint,nlev(j))

                           ireg = ntemp

                      if( ntemp .gt. 4000000 .and.
     &                    ntemp .lt. 5000000 ) then

                           if( icgg .eq. 0 ) goto 995

                           ireg = ireg - 4000000

                           if( idgr(ireg) .eq. 0 ) goto 998
                           if( lat(1,idgr(ireg)) .eq. 0 ) goto 997

                           nlev(j) = nlev(j) + 1

                           ipar(kint) = ipar(kint) + 1
                           if( ipar(kint) .gt. ipmax )
     &                                         ipmax = ipar(kint)
                           mtrg(igm+kint,ipar(kint)) =
     &                     mtrg(igm+j+kint,nlev(j))
                           ltemp = mtrg(igm+j+kint,nlev(j))

                        do ll = 1, ltemp * 6

                           nlev(j) = nlev(j) + 1

                           ipar(kint) = ipar(kint) + 1
                           if( ipar(kint) .gt. ipmax )
     &                                         ipmax = ipar(kint)
                           mtrg(igm+kint,ipar(kint)) =
     &                     mtrg(igm+j+kint,nlev(j))

                        end do

                      else

                           if( idgr(ireg) .eq. 0 ) goto 998

                      end if

                    end do

                  else if( mtemp .gt. 4000000 .and.
     &                     mtemp .lt. 5000000 ) then

                           if( icgg .eq. 0 ) goto 995

                           ireg = ireg - 4000000

                           if( idgr(ireg) .eq. 0 ) goto 998
                           if( lat(1,idgr(ireg)) .eq. 0 ) goto 997

                           nlev(j) = nlev(j) + 1

                           ipar(kint) = ipar(kint) + 1
                           if( ipar(kint) .gt. ipmax )
     &                                         ipmax = ipar(kint)
                           mtrg(igm+kint,ipar(kint)) =
     &                     mtrg(igm+j+kint,nlev(j))
                           ltemp = mtrg(igm+j+kint,nlev(j))

                        do ll = 1, ltemp * 6

                           nlev(j) = nlev(j) + 1

                           ipar(kint) = ipar(kint) + 1
                           if( ipar(kint) .gt. ipmax )
     &                                         ipmax = ipar(kint)
                           mtrg(igm+kint,ipar(kint)) =
     &                     mtrg(igm+j+kint,nlev(j))

                        end do

                  else

                           if( idgr(ireg) .eq. 0 ) goto 998

                  end if

               end do

                        goto 100

  300          continue

                        kpar = kint
                        ibra = ibrs
                        mlev = 0

            else

                        lpar(kpar) = lpar(kpar) + 1
                        ipar(kpar) = ipar(kpar) + 1
                        if( ipar(kpar) .gt. ipmax ) ipmax = ipar(kpar)
                        mtrg(igm+kpar,ipar(kpar)) = ntrg(k)

                        ireg = ntrg(k)

                  if( ntrg(k) .gt. 4000000 .and.
     &                ntrg(k) .lt. 5000000 ) then

                        if( icgg .eq. 0 ) goto 995

                        ireg = ireg - 4000000

                        if( idgr(ireg) .eq. 0 ) goto 998
                        if( lat(1,idgr(ireg)) .eq. 0 ) goto 997

                        k = k + 1
                        nlat = ntrg(k)

                     if( ibra .gt. 0 ) then

                           ipar(kpar) = ipar(kpar) + 1
                           if( ipar(kpar) .gt. ipmax )
     &                                         ipmax = ipar(kpar)
                           mtrg(igm+kpar,ipar(kpar)) = nlat

                        do j = 1, nlat * 6

                           k = k + 1
                           ipar(kpar) = ipar(kpar) + 1
                           if( ipar(kpar) .gt. ipmax )
     &                                         ipmax = ipar(kpar)
                           mtrg(igm+kpar,ipar(kpar)) = ntrg(k)

                        end do

                     else

                                 lnum = 0

                        do j = 1, nlat

                           do m = 1, 6

                              k = k + 1
                              latj(m) = ntrg(k)

                           end do

                           do m3 = latj(5), latj(6)
                           do m2 = latj(3), latj(4)
                           do m1 = latj(1), latj(2)

                                 lnum = lnum + 1

                              if( lnum .ne. 1 ) then

                                 lpar(kpar) = lpar(kpar) + 1
                                 ipar(kpar) = ipar(kpar) + 1
                                 if( ipar(kpar) .gt. ipmax )
     &                                               ipmax = ipar(kpar)
                                 mtrg(igm+kpar,ipar(kpar)) =
     &                                               ireg + 4000000

                              end if

                                 ipar(kpar) = ipar(kpar) + 1
                                 mtrg(igm+kpar,ipar(kpar)) = 1

                                 ipar(kpar) = ipar(kpar) + 1
                                 mtrg(igm+kpar,ipar(kpar)) = m1
                                 ipar(kpar) = ipar(kpar) + 1
                                 mtrg(igm+kpar,ipar(kpar)) = m1
                                 ipar(kpar) = ipar(kpar) + 1
                                 mtrg(igm+kpar,ipar(kpar)) = m2
                                 ipar(kpar) = ipar(kpar) + 1
                                 mtrg(igm+kpar,ipar(kpar)) = m2
                                 ipar(kpar) = ipar(kpar) + 1
                                 mtrg(igm+kpar,ipar(kpar)) = m3
                                 ipar(kpar) = ipar(kpar) + 1
                                 mtrg(igm+kpar,ipar(kpar)) = m3

                                 if( ipar(kpar) .gt. ipmax )
     &                               ipmax = ipar(kpar) + 7
                           end do
                           end do
                           end do

                        end do

                     end if

                  else

                        if( ireg .gt. 0 .and. ireg .lt. 1000000) then
                            if(idgr(ireg) .eq. 0 ) goto 998
                        endif

                  end if

            end if

            goto 510

  520       continue

*-----------------------------------------------------------------------

               ntrn = lpar(0)
               mtrn = ipar(0)

            if( ipmax > MAX_NUM_MTRG ) then
               write(ErrCha,'(a,a,i5,a,a,i5,a)')
     &              'sub.tregion2@tallsm1.f ?dimension over mtrg?'
     &                 //' ipmax > MAX_NUM_MTRG'
     &              ,' (ipmax=',ipmax,')'
     &              ,' (MAX_NUM_MTRG@moddas.f=',MAX_NUM_MTRG,')'
               ErrID = 'L:6971/R:tregion2/F:tallsm1.f'
               call ErrWrite(ErrID,ErrCha)
            endif

            do k = 1, mtrn
               idas_region( igm_argument + k - 1) = mtrg(igm+0,k)
            end do

*-----------------------------------------------------------------------


*-----------------------------------------------------------------------

      return

*-----------------------------------------------------------------------

  993 continue
                  write(io,'(/'' **** Error'',
     &            '' : region is too many or memory is lack'')')
                  ErrCha = ''
                  MsgID = 'L:6992/R:tregion2/F:tallsm1.f'
                  call ErrWrite(MsgID, ErrCha)
                  write(jo,'(/'' **** Error'',
     &            '' : region is too many or memory is lack'')')
      ierr = 2
      return

*-----------------------------------------------------------------------

  994 continue
                  write(io,'(/'' **** Error'',
     &            '' : fill is used but giometry is CG'')')
                  ErrCha = ''
                  MsgID = 'L:7005/R:tregion2/F:tallsm1.f'
                  call ErrWrite(MsgID, ErrCha)
                  write(jo,'(/'' **** Error'',
     &            '' : fill is used but giometry is CG'')')
      ierr = 1
      return

*-----------------------------------------------------------------------

  995 continue
                  write(io,'(/'' **** Error'',
     &            '' : lattice is used but geometry is CG'')')
                  ErrCha = ''
                  MsgID = 'L:7018/R:tregion2/F:tallsm1.f'
                  call ErrWrite(MsgID, ErrCha)
                  write(jo,'(/'' **** Error'',
     &            '' : lattice is used but geometry is CG'')')
      ierr = 1
      return

*-----------------------------------------------------------------------

  996 continue
                  write(io,'(/'' **** Error'',
     &            '' : universe is used but giometry is CG'')')
                  ErrCha = ''
                  MsgID = 'L:7031/R:tregion2/F:tallsm1.f'
                  call ErrWrite(MsgID, ErrCha)
                  write(jo,'(/'' **** Error'',
     &            '' : universe is used but giometry is CG'')')
      ierr = 1
      return

*-----------------------------------------------------------------------

  997 continue
                  write(io,'(/'' **** Error'',
     &            '' : region '',i6,'' is not'',
     &            '' lattice region.'')') ireg
                  ErrCha = ''
                  MsgID = 'L:7045/R:tregion2/F:tallsm1.f'
                  call ErrWrite(MsgID, ErrCha)
                  write(jo,'(/'' **** Error'',
     &            '' : region '',i6,'' is not'',
     &            '' lattice region.'')') ireg
      ierr = 1
      return

*-----------------------------------------------------------------------

  998 continue
                  write(io,'(/'' **** Error'',
     &            '' : region '',i6,'' is not'',
     &            '' defined in geometry.'')') ireg
                  ErrCha = ''
                  MsgID = 'L:7060/R:tregion2/F:tallsm1.f'
                  call ErrWrite(MsgID, ErrCha)
                  write(jo,'(/'' **** Error'',
     &            '' : region '',i6,'' is not'',
     &            '' defined in geometry.'')') ireg

*-----------------------------------------------------------------------

  999 continue

      ierr = 1

      return
      end


************************************************************************
*                                                                      *
      subroutine tregion3(icn,chlw,i1,i3,ic,
     &                    igm,ipmax,ipar,jpar,
     &                    ibra,ilat,klat,ifis,kpar,iuni
     &                    ,ndim_mtrg,mtrg)
*                                                                      *
*       region part 2                                                  *
*       modified by K.Niita on 2011/02/03                              *
*                                                                      *
************************************************************************
      use moddas

      implicit real*8 (a-h,o-z)

      include 'param.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      character chlw*200

      dimension lats(1000,6)

      dimension ipar(0:20)
      dimension jpar(0:20)

      logical deqn4
      logical dnen1

      save lats, ntrf, ntri, mlat

      integer, intent(in) :: ndim_mtrg
      integer, intent(inout) :: mtrg(0:20,ndim_mtrg)

*-----------------------------------------------------------------------

            icn = 0

*-----------------------------------------------------------------------

               if( chlw(ic:ic) .eq. '{' ) then !FURUTA20230313

                     if( iuni .ne. 0 ) goto 986
                     if( ibra .ne. 0 ) goto 992

                           ibra = 1

                           ic = ic + 1

               else if( ibra .eq. 2 .and. chlw(ic:ic) .eq. '-' ) then

                     if( iuni .ne. 0 ) goto 986

                           ibra = 3

                           ic = ic + 1

               else if( chlw(ic:ic) .eq. '-' .and. ilat .eq. 0 ) then

                     goto 980

               else if( chlw(ic:ic) .eq. '}' ) then

                     if( iuni .ne. 0 ) goto 986
                     if( ibra .ne. 4 ) goto 992

                           ibra = 0

                           ic = ic + 1

                     if( ntrf .le. ntri ) goto 995

                     do j = 0, ntrf - ntri

                           ipar(kpar) = ipar(kpar) + 1
                           if( ipar(kpar) .gt. ipmax )
     &                                         ipmax = ipar(kpar)
                           mtrg(igm+kpar,ipar(kpar)) = ntri + j

                     end do

               else if( chlw(ic:ic) .eq. ':' ) then

                           if( iuni .ne. 0 ) goto 986
                           if( ibra .ne. 0 ) goto 992
                           if( ilat .eq. 0 ) goto 989
                           if( mlat .ne. 0 ) goto 989

                           mlat = 1

                           ic = ic + 1

               else if( chlw(ic:ic) .eq. ',' ) then

                           if( iuni .ne. 0 ) goto 986
                           if( ibra .ne. 0 ) goto 992
                           if( ilat .eq. 0 ) goto 989
                           if( mlat .ne. 0 ) goto 990
                           if( klat .le. 4 ) goto 990

                        if( klat .eq. 5 ) then

                           klat = klat + 1
                           lats(ilat,klat) = lats(ilat,klat-1)

                        end if

                           ilat = ilat + 1
                           if( ilat .gt. 1000 ) goto 984
                           klat = 0
                           mlat = 0

                           ic = ic + 1

               else if( chlw(ic:ic) .eq. ']' ) then

                           if( iuni .ne. 0 ) goto 986
                           if( ibra .ne. 0 ) goto 992
                           if( ilat .eq. 0 ) goto 989
                           if( mlat .ne. 0 ) goto 990
                           if( klat .le. 4 ) goto 990

                        if( klat .eq. 5 ) then

                           klat = klat + 1
                           lats(ilat,klat) = lats(ilat,klat-1)

                        end if

                           ipar(kpar) = ipar(kpar) + 1
                           if( ipar(kpar) .gt. ipmax )
     &                                         ipmax = ipar(kpar)
                           mtrg(igm+kpar,ipar(kpar)) = ilat

                     do k = 1, ilat

                        do l = 1, 6

                           ipar(kpar) = ipar(kpar) + 1
                           if( ipar(kpar) .gt. ipmax )
     &                                         ipmax = ipar(kpar)
                           mtrg(igm+kpar,ipar(kpar)) = lats(k,l)

                        end do

                     end do

                           ilat = 0
                           klat = 0
                           mlat = 0

                           ic = ic + 1

                else if( chlw(ic:ic) .eq. '[' .and. ic .gt. 5 ) then ! T.Sato 2023/04/14 back to original definition because [ is used for section name

                           if( iuni .ne. 0 ) goto 986
                           if( ibra .ne. 0 ) goto 992
                           if( ilat .ne. 0 ) goto 990

                           mtrg(igm+kpar,ipar(kpar)) =
     &                     mtrg(igm+kpar,ipar(kpar)) + 4000000

                           ilat = 1
                           klat = 0
                           mlat = 0

                           ic = ic + 1

               else if( chlw(ic:ic) .eq. 'u' ) then

                           ic0 = ic
                           ic = jnumc(chlw,ic+1,i3)

                           if( ic0 .eq. i1 .and.
     &                         chlw(ic:ic) .ne. '=' ) goto 500

                           if( chlw(ic:ic) .ne. '=' ) goto 986

                           if( iuni .ne. 0 ) goto 986
                           if( ibra .ne. 0 ) goto 992
                           if( ilat .ne. 0 ) goto 990

                           ic = ic + 1

                           iuni = 1

               else if( chlw(ic:ic) .eq. '<' ) then

                           if( iuni .ne. 0 ) goto 986
                           if( ilat .ne. 0 ) goto 990
                           if( ibra .ne. 0 ) goto 992
                           if( kpar .eq. 0 ) goto 987

                           ipar(kpar) = ipar(kpar) + 1
                           if( ipar(kpar) .gt. ipmax )
     &                                         ipmax = ipar(kpar)
                           mtrg(igm+kpar,ipar(kpar)) = 3000000

                           jpar(kpar-1) = jpar(kpar-1) + 1
                           if( jpar(kpar-1) .gt. 10 ) goto 985

                           ic = ic + 1

               else if( chlw(ic:ic) .eq. '(' ) then !FURUTA202300313

                           if( iuni .ne. 0 ) goto 986
                           if( ibra .ne. 0 ) goto 992
                           if( ilat .ne. 0 ) goto 990

                           jpar(kpar) = 0
                           ipar(kpar) = ipar(kpar) + 1
                           if( ipar(kpar) .gt. ipmax )
     &                                         ipmax = ipar(kpar)
                           mtrg(igm+kpar,ipar(kpar)) = 1000000

                           kpar = kpar + 1
                           if( kpar .gt. 10 ) goto 988

                           ipar(kpar) = 0

                           ic = ic + 1

               else if( chlw(ic:ic) .eq. ')' ) then

                           if( iuni .ne. 0 ) goto 986
                           if( ipar(kpar) .eq. 0 ) goto 993

                           kpar = kpar - 1

                     if( kpar .eq. 0 .and.
     &                   mtrg(igm+1,1) .eq. 5000000 .and.
     &                   ipar(1) .eq. 1 ) then

                              mtrg(igm+kpar,ipar(kpar)) = 6000000

                     else

                        if( jpar(kpar) .gt. 0 ) then

                              mtrg(igm+kpar,ipar(kpar)) = 2000000

                        end if

                        do k = 1, ipar(kpar+1)

                              if( mtrg(igm+kpar+1,k) .eq. 5000000 )
     &                            goto 991

                              ipar(kpar) = ipar(kpar) + 1
                              if( ipar(kpar) .gt. ipmax )
     &                                            ipmax = ipar(kpar)
                              mtrg(igm+kpar,ipar(kpar)) =
     &                        mtrg(igm+kpar+1,k)

                        end do

                              ipar(kpar) = ipar(kpar) + 1
                              if( ipar(kpar) .gt. ipmax )
     &                                            ipmax = ipar(kpar)

                           if( jpar(kpar) .gt. 0 ) then

                              mtrg(igm+kpar,ipar(kpar)) = 2000001

                           else

                              mtrg(igm+kpar,ipar(kpar)) = 1000001

                           end if

                        end if

                           ic = ic + 1
                           jpar(kpar) = 0

               else if( chlw(ic:ic+2) .eq. 'all' ) then

                        if( iuni .ne. 0 ) goto 986

                           ipar(kpar) = ipar(kpar) + 1
                           if( ipar(kpar) .gt. ipmax )
     &                                         ipmax = ipar(kpar)
                           mtrg(igm+kpar,ipar(kpar)) = 5000000

                           ic = ic + 3

*-----------------------------------------------------------------------

               else if( deqn4( chlw(ic:ic) ) ) then

                           ici = ic

                     do j = ic + 1, i3

                        if( dnen1( chlw(j:j) ) ) goto 167

                     end do

  167                      icf = j - 1

                           if( icf .lt. ici ) goto 997

                     call onum(chlw,ici,icf,cvvv,ierr)

                        if( ierr .ne. 0 ) goto 997

*-----------------------------------------------------------------------

                  if( ibra .eq. 1 ) then

                           ntri = nint( cvvv )

                           ibra = 2

                  else if( ibra .eq. 3 ) then

                           ntrf = nint( cvvv )

                           ibra = 4

                  else if( ibra .eq. 2 .or. ibra .eq. 4 ) then

                           goto 992

                  else if( ilat .gt. 0 ) then

                           klat = klat + 1

                     if( klat/2*2 .eq. klat ) then

                        if( mlat .eq. 0 ) then

                           lats(ilat,klat) = lats(ilat,klat-1)
                           klat = klat + 1

                        else

                           mlat = 0

                        end if

                     end if

                           lats(ilat,klat) = nint( cvvv )

                  else

                           iadd = 0
                           if( iuni .ne. 0 ) iadd = 7000000
                           iuni = 0

                           ipar(kpar) = ipar(kpar) + 1
                           if( ipar(kpar) .gt. ipmax )
     &                                         ipmax = ipar(kpar)
                           mtrg(igm+kpar,ipar(kpar)) =
     &                     nint( cvvv ) + iadd

                  end if

                           ic = icf + 1

               else if( ifis .eq. 0 .or.
     &                ( ifis .gt. 0 .and. ic .gt. i1 ) ) then

                        goto 997

               else if( ifis .gt. 0 .and. ic .eq. i1 ) then

                        goto 500

               end if
               if( ipmax > ndim_mtrg ) then
                  write(ErrCha,'(a,a,i5,a,a,i5,a)')
     &                 'sub.tregion3@tallsm1.f ?dimension over mtrg?'
     &                    //' ipmax > ndim_mtrg'
     &                 ,' (ipmax=',ipmax,')'
     &                 ,' ndim_mtrg=',ndim_mtrg,')'
                  ErrID = 'L:7454/R:tregion3/F:tallsm1.f'
                  call ErrWrite(ErrID,ErrCha)
               endif

      return

*-----------------------------------------------------------------------

  500 continue

               if( ibra .ne. 0 ) goto 997
               if( kpar .ne. 0 ) goto 993

               icn = 500
               if( ipmax > ndim_mtrg ) then
                  write(ErrCha,'(a,a,i5,a,a,i5,a)')
     &                 'sub.tregion3@tallsm1.f ?dimension over mtrg?'
     &                    //' ipmax > ndim_mtrg'
     &                 ,' (ipmax=',ipmax,')'
     &                 ,' (ndim_mtrg=',ndim_mtrg,')'
                  ErrID = 'L:7474/R:tregion3/F:tallsm1.f'
                  call ErrWrite(ErrID,ErrCha)
               endif

      return

*-----------------------------------------------------------------------

  980 continue
      icn = 980
      return
  984 continue
      icn = 984
      return
  985 continue
      icn = 985
      return
  986 continue
      icn = 986
      return
  987 continue
      icn = 987
      return
  988 continue
      icn = 988
      return
  989 continue
      icn = 989
      return
  990 continue
      icn = 990
      return
  991 continue
      icn = 991
      return
  992 continue
      icn = 992
      return
  993 continue
      icn = 993
      return
  995 continue
      icn = 995
      return
  997 continue
      icn = 997
      return

*-----------------------------------------------------------------------

      end


************************************************************************
*                                                                      *
      subroutine tregion4(ierr,ntrn,mtrn,ntrg,ign,igm_argument
     &                    ,ndim_region,idas_region)
*                                                                      *
*       region part 3                                                  *
*       modified by K.Niita on 2011/02/03                              *
*                                                                      *
************************************************************************
      use moddas
      use moddas_region_mtrg

      implicit real*8 (a-h,o-z)

      include 'param.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      dimension ntrg(mtrn)

      dimension ipar(0:20)
      dimension jpar(0:20)
      dimension lpar(0:20)

      integer, intent(in) :: igm_argument
      integer, intent(in) :: ndim_region
      integer, intent(inout) :: idas_region(ndim_region)

      igm = 0
      call moddas_allocate_int2(0, 20, MAX_NUM_MTRG, mtrg)

*-----------------------------------------------------------------------


*-----------------------------------------------------------------------

                  ibra  = 0
                  kpar  = 0
                  ipmax = 0
                  ierr = 0

               do i = 0, 20
                  ipar(i) = 0
                  jpar(i) = 0
                  lpar(i) = 0
               end do

         k = 0
  510    k = k + 1

         if( k .gt. mtrn ) goto 520

            if( ntrg(k) .eq. 1000000 ) then

                        ibra = ibra + 1

               if( ibra .eq. 1 ) then

                        lpar(kpar) = lpar(kpar) + 1

                        ipar(kpar) = ipar(kpar) + 1
                        if( ipar(kpar) .gt. ipmax ) ipmax = ipar(kpar)
                        mtrg(igm+kpar,ipar(kpar)) = 1000000

                        kpar = kpar + 1
                        ipar(kpar) = 0
                        lpar(kpar) = 0

               end if

            else if( ntrg(k) .eq. 1000001 ) then

                        ibra = ibra - 1

               if( ibra .eq. 0 ) then

                        kpar = kpar - 1

                        ipar(kpar) = ipar(kpar) + 1
                        if( ipar(kpar) .gt. ipmax ) ipmax = ipar(kpar)
                        mtrg(igm+kpar,ipar(kpar)) = lpar(kpar+1)

                     do j = 1, ipar(kpar+1)

                        ipar(kpar) = ipar(kpar) + 1
                        if( ipar(kpar) .gt. ipmax ) ipmax = ipar(kpar)
                        mtrg(igm+kpar,ipar(kpar)) = mtrg(igm+kpar+1,j)

                     end do

                        ipar(kpar) = ipar(kpar) + 1
                        if( ipar(kpar) .gt. ipmax ) ipmax = ipar(kpar)
                        mtrg(igm+kpar,ipar(kpar)) = 1000001

               end if

            else if( ntrg(k) .eq. 3000000 ) then

                        kpar = kpar + 1
                        ipar(kpar) = 0
                        lpar(kpar) = 0
                        jpar(kpar) = jpar(kpar-1) + 1

                        ipar(kpar) = ipar(kpar) + 1
                        if( ipar(kpar) .gt. ipmax ) ipmax = ipar(kpar)
                        mtrg(igm+kpar,ipar(kpar)) = 3000000

            else if( ntrg(k) .eq. 2000000 ) then

                        kpar = kpar + 1
                        ipar(kpar) = 0
                        lpar(kpar) = 0
                        jpar(kpar) = 1

                        ibrs = ibra
                        ibra = 0

            else if( ntrg(k) .eq. 2000001 ) then

                        nelm = 1
                        kint = kpar - jpar(kpar)

                     do j = kint + 1, kpar

                        nelm = nelm * lpar(j)

                     end do

                        lpar(kint) = lpar(kint) + nelm

                        ipar(kint) = ipar(kint) + 1
                        if( ipar(kint) .gt. ipmax ) ipmax = ipar(kint)
                        mtrg(igm+kint,ipar(kint)) = 2000000
                        ipar(kint) = ipar(kint) + 1
                        if( ipar(kint) .gt. ipmax ) ipmax = ipar(kint)
                        mtrg(igm+kint,ipar(kint)) = jpar(kpar)

                     do j = kint + 1, kpar

                        ipar(kint) = ipar(kint) + 1
                        if( ipar(kint) .gt. ipmax ) ipmax = ipar(kint)
                        mtrg(igm+kint,ipar(kint)) = lpar(j)

                     end do

                  do j = kint + 1, kpar

                     do l = 1, ipar(j)

                        ipar(kint) = ipar(kint) + 1
                        if( ipar(kint) .gt. ipmax ) ipmax = ipar(kint)
                        mtrg(igm+kint,ipar(kint)) = mtrg(igm+j,l)

                     end do

                  end do

                        ipar(kint) = ipar(kint) + 1
                        if( ipar(kint) .gt. ipmax ) ipmax = ipar(kint)
                        mtrg(igm+kint,ipar(kint)) = 2000001

                        kpar = kint
                        ibra = ibrs

            else

                        lpar(kpar) = lpar(kpar) + 1
                        ipar(kpar) = ipar(kpar) + 1
                        if( ipar(kpar) .gt. ipmax ) ipmax = ipar(kpar)
                        mtrg(igm+kpar,ipar(kpar)) = ntrg(k)

                  if( ntrg(k) .gt. 4000000 .and.
     &                ntrg(k) .lt. 5000000 ) then

                        k = k + 1
                        nlat = ntrg(k)
                        ipar(kpar) = ipar(kpar) + 1
                        if( ipar(kpar) .gt. ipmax ) ipmax = ipar(kpar)
                        mtrg(igm+kpar,ipar(kpar)) = nlat

                     do j = 1, nlat * 6

                        k = k + 1
                        ipar(kpar) = ipar(kpar) + 1
                        if( ipar(kpar) .gt. ipmax ) ipmax = ipar(kpar)
                        mtrg(igm+kpar,ipar(kpar)) = ntrg(k)

                     end do

                  end if

            end if

            goto 510

*-----------------------------------------------------------------------

  520       continue

               ntrn = lpar(0)
               mtrn = ipar(0)

            if( ipmax > MAX_NUM_MTRG ) then
               write(ErrCha,'(a,a,i5,a,a,i5,a)')
     &              'sub.tregion4@tallsm1.f ?dimension over mtrg?'
     &                 //' ipmax > MAX_NUM_MTRG'
     &              ,' (ipmax=',ipmax,')'
     &              ,' (MAX_NUM_MTRG@moddas.f=',MAX_NUM_MTRG,')'
               ErrID = 'L:7736/R:tregion4/F:tallsm1.f'
               call ErrWrite(ErrID,ErrCha)
            endif

            do k = 1, mtrn
               idas_region( ign + k - 1 ) = mtrg(igm+0,k)
            end do

*-----------------------------------------------------------------------


*-----------------------------------------------------------------------

      return

  993 ierr = 1

      return
      end

************************************************************************
*                                                                      *
      subroutine tregion5(chin,ic,i3,ic2,ntrn,mtrn,igm_argument,ierr
     &                    ,ndim_region,idas_region)
*                                                                      *
*       read region in colomn                                          *
*       modified by K.Niita on 2011/02/03                              *
*                                                                      *
************************************************************************
      use moddas
      use moddas_region_mtrg

      implicit real*8 (a-h,o-z)

      include 'param.inc'

*-----------------------------------------------------------------------

      dimension ipar(0:20)
      dimension jpar(0:20)
      dimension lpar(0:20)

      character chin*200
      character chlw*200

      logical deqn4

      data klnmax /1000000/

      integer, intent(in) :: igm_argument
      integer, intent(in) :: ndim_region
      integer, intent(inout) :: idas_region(ndim_region)

*-----------------------------------------------------------------------

         ierr = 0
         ipmax = 0 !FURUTA

*-----------------------------------------------------------------------

         if( chin(ic:ic) .eq. '(' ) then

                  ibra = 1

                  chlw(ic+5:ic+5) = chin(ic:ic)

            do i = ic+1, i3

                  chlw(i+5:i+5) = chin(i:i)

               if( chin(i:i) .eq. ')' ) then

                  ibra = ibra - 1

                  if( ibra .eq. 0 ) goto 100

               else if( chin(i:i) .eq. '(' ) then

                  ibra = ibra + 1

               end if

            end do

               goto 999

  100       continue

               ici = ic+5
               icf = i+5
               ic2 = i+1

         else if( chin(ic:ic) .eq. '{' ) then

                  chlw(ic+5:ic+5) = '('
                  chlw(ic+6:ic+6) = chin(ic:ic)

            do i = ic+1, i3

                  chlw(i+6:i+6) = chin(i:i)

               if( chin(i:i) .eq. '}' ) goto 200

            end do

               goto 999

  200       continue

                  chlw(i+7:i+7) = ')'

               ici = ic+5
               icf = i+7
               ic2 = i+1

         else if( deqn4( chin(ic:ic) ) ) then

               call snum(chin,ic,i3,ic2,cvvv,ierr)

                  if( ierr .ne. 0 ) goto 999

               igm = 0
               call moddas_allocate_int2(0, 20, 1, mtrg)
               ipar(0) = 1
               mtrg(igm+0,1) = nint( cvvv )

               goto 500

         else

               goto 999

         end if

*-----------------------------------------------------------------------

               ic = ici
               i1 = ici
               j3 = icf

               ntrn = 0
               mtrn = 0
               ibra = 0
               ilat = 0
               klat = 0
               ifis = 0
               kpar = 0
               iuni = 0

               nvol = 0

               do i = 0, 20
                  ipar(i) = 0
                  jpar(i) = 0
                  lpar(i) = 0
               end do

               igm = 0
               call moddas_allocate_int2(0, 20, MAX_NUM_MTRG, mtrg)
               ipmax = 0

            do i = 1, klnmax

                  ic = jnumc(chlw,ic,j3)

                  call tregion3(icn,chlw,i1,j3,ic,
     &                          igm,ipmax,ipar,jpar,
     &                          ibra,ilat,klat,ifis,kpar,iuni
     &                          ,MAX_NUM_MTRG,mtrg)


                     if( icn .gt. 900 ) goto 999
                     if( icn .eq. 500 ) goto 500
                     if( ic  .gt.  j3 ) goto 500

            end do

                  goto 999

*-----------------------------------------------------------------------

  500 continue

               ntrn = ipar(0)
               mtrn = ipar(0)

            do k = 1, mtrn
               idas_region( igm_argument + k - 1) = mtrg(igm+0,k)
            end do

            call tregion4(ierr,ntrn,mtrn,idas_region(igm_argument)
     &                    ,igm_argument,ign,ndim_region,idas_region)


               if( ierr .ne. 0 ) goto 999


*-----------------------------------------------------------------------

      return

 999  continue
      ierr = 1

      return
      end

************************************************************************
*                                                                      *
      subroutine ttetmesh(jsn,jsi,dsin,idsi,ill,ilf,
     &     jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr,
     &     mtrn,ntrn,
     &     ndim_region,idas_region)
*                                                                      *
*       read 'tet =' sub-section of input tally section                *
*       Modified by T.Furuta on 2025/01/16                             *
*                                                                      *
************************************************************************
      use moddas

      implicit real*8 (a-h,o-z)

      include 'param.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      character m_err*200
      common /error/ m_err, l_err, k_err

*-----------------------------------------------------------------------

      character(200),intent(out) :: chin, chlw, chcm

      character dsin(0:9)*200
      dimension idsi(0:9)

      dimension ill(0:9), ilf(0:9)

*-----------------------------------------------------------------------

      integer ipar

      data klnmax /1000000/

      integer, intent(in) :: ndim_region
      integer, intent(inout) :: idas_region(ndim_region)
      integer,allocatable :: mtetreg(:)

*-----------------------------------------------------------------------

            ierr = 0

*-----------------------------------------------------------------------
*     read one line from jsi
*-----------------------------------------------------------------------

  140 continue

            call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

               if( ierr .ne. 0 ) return
               if( jpn  .eq. 3 ) return

               if( iskip .ne. 0 ) goto 140

  150 continue

               if( ierr .ne. 0 ) return
               if( jpn  .eq. 3 ) return

*-----------------------------------------------------------------------
*        error not for reg =
*-----------------------------------------------------------------------

               if( chcm(i1:i1+3) .ne. 'reg=' ) goto 998

               ic = inumc(chlw,i1,i3,'=') + 1

*-----------------------------------------------------------------------
*        read region number
*-----------------------------------------------------------------------

               ic = jnumc(chlw,ic,i3)

               ifis = 0
               ipar = 0
               jpar = 0
               kpar = 0

               call moddas_allocate_int(MAX_NUM_MTRG,mtetreg)

            do i = 1, klnmax

                  ic = jnumc(chlw,ic,i3)

                  call ttetmesh3(icn,chlw,i1,i3,ic,
     &                 ipar,jpar,kpar,ifis,
     &                 max_num_MTRG,mtetreg)

                     if( icn .gt. 900 ) goto 900
                     if( icn .eq. 500 ) goto 500

               if( i .lt. klnmax .and. ic .gt. i3 ) then

  147             call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

                     if( ierr .ne. 0 ) return
                     if( jpn  .eq. 3 ) goto 500

                     if( iskip .ne. 0 ) goto 147

                     ifis = ifis + 1

                     ic = i1

               end if

            end do

                  goto 996

*-----------------------------------------------------------------------
*        modify the input
*-----------------------------------------------------------------------

  500    continue

         if(ipar.gt.1)then
          icn=996
          goto 900
         endif

         ntrn = kpar+1
         mtrn = 3+kpar

         ndsm=1
         do k=1,mtrn
          idas_region( ndsm + k - 1) = mtetreg(k)
         enddo

         call moddas_deallocate_int(mtetreg)

*-----------------------------------------------------------------------
*     errors
*-----------------------------------------------------------------------

         return

  900 continue

         if( icn .eq. 996 ) goto 996
         if( icn .eq. 997 ) goto 997

*-----------------------------------------------------------------------

  996 continue

         m_err = 'Only 1 region is allowed with mesh = tet. '
     &   // 'Use (*** ***) for more than one TETRA regions.'
         ErrCha = ''
         ErrID = 'L:8099/R:ttetmesh/F:tallsm1.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  997 continue

         m_err = 'Description of region number is wrong.'
         ErrCha = ''
         ErrID = 'L:8111/R:ttetmesh/F:tallsm1.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  998 continue
         m_err = 'After { mesh = tet } line should be { reg = }.'
         ErrCha = ''
         ErrID = 'L:8122/R:ttetmesh/F:tallsm1.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1

         return

*-----------------------------------------------------------------------

      end


************************************************************************
*                                                                      *
      subroutine ttetmesh2(io,jo,ierr,ntrn,mtrn,ntrg,
     &     ndim_region,idas_region,ielem2ir)
*                                                                      *
*       identify cell of tet and compute element volume                *
*       Modified by T.Furuta on 2025/01/16                             *
*                                                                      *
************************************************************************
      use TETRAMOD, only:nelem,tetratvol,ntetcl,itetst,itetlist,
     &     itetelem,nelemtot
      implicit real*8 (a-h,o-z)

      include 'param.inc'
      include 'err.inc'

*-----------------------------------------------------------------------
      integer,intent(in) :: io,jo,ndim_region,ntrg(mtrn)
      integer,intent(inout) :: ntrn,mtrn,idas_region(ndim_region)
      integer,intent(out) :: ierr,ielem2ir(nelemtot)
*-----------------------------------------------------------------------
      common /regdc/ idrg(kvlmax), idgr(kvmmax)

      integer nlat3,itetcl(10),itetvol,itetauto,itgchk
      common /tetf1/ nlat3,itetcl,itetvol,itetauto,itgchk

*-----------------------------------------------------------------------

      integer :: itettal
      common /itettal1/ itettal
      data itettal /0/

*-----------------------------------------------------------------------

      logical iflag,jflag

*-----------------------------------------------------------------------

      call ttetmesh5(io,jo,ierr,ntrn,mtrn,ntrg,
     &     ndim_region,idas_region,ielem2ir)

      if(ierr.eq.0)then
*-----------------------------------------------------------------------
*        volume
*-----------------------------------------------------------------------
       if(itettal.eq.0)then     ! only once
        call tetratvol(nlat3)
       endif
       itettal=itettal+1
*-----------------------------------------------------------------------
      endif

      return
      end


************************************************************************
*                                                                      *
      subroutine ttetmesh3(icn,chlw,i1,i3,ic,
     &     ipar,jpar,kpar,ifis,
     &     ndim_mtetreg,mtetreg)
*                                                                      *
*       tetmesh part 2                                                 *
*       read cell number                                               *
*       Modified by T.FURUTA on 2025/01/16                             *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      include 'param.inc'

*-----------------------------------------------------------------------

      integer,intent(in) :: i1,i3
      integer,intent(out) :: icn
      integer,intent(inout) :: ic,ipar,ifis
      character(200),intent(in) :: chlw

      logical deqn4
      logical dnen1

      integer, intent(in) :: ndim_mtetreg
      integer, intent(out) :: mtetreg(ndim_mtetreg)

*-----------------------------------------------------------------------

            icn = 0

*-----------------------------------------------------------------------

               if( chlw(ic:ic) .eq. '(' ) then

                jpar = 1
                ic=ic+1

               else if( chlw(ic:ic) .eq. ')' ) then

                if( jpar .eq. 0) goto 997
                if( kpar .eq. 1 .and. ipar .eq. 0) goto 997

                mtetreg(3) = kpar

                jpar = 0
                ic=ic+1

               else if( chlw(ic:ic) .eq. '<' ) then

                if( jpar .ne. 1) goto 997

                jpar = 2
                ic=ic+1

               elseif( deqn4( chlw(ic:ic) ) ) then

                           ici = ic

                     do j = ic + 1, i3

                        if( dnen1( chlw(j:j) ) ) goto 167

                     end do

  167                      icf = j - 1

                           if( icf .lt. ici ) goto 997

                     call onum(chlw,ici,icf,cvvv,ierr)

                        if( ierr .ne. 0 ) goto 997

*-----------------------------------------------------------------------

                        if(jpar .eq. 1)then

                         kpar = kpar + 1

                         mtetreg(3+kpar) = nint( cvvv )

                        elseif(jpar .eq. 0 .or. jpar .eq. 2)then

                         ipar = ipar + 1

                         mtetreg(1) = nint( cvvv )

                        endif

                           ic = icf + 1

               else if( ifis .eq. 0 .or.
     &                ( ifis .gt. 0 .and. ic .gt. i1 ) ) then

                        goto 997

               else if( ifis .gt. 0 .and. ic .eq. i1 ) then

                        goto 500

               end if

      return

*-----------------------------------------------------------------------

  500 continue

               icn = 500

      return

*-----------------------------------------------------------------------

  997 continue
      icn = 997
      return

*-----------------------------------------------------------------------

      end


************************************************************************
*                                                                      *
      subroutine ttetmesh4(io,jo,ntrn,mtrn,ntrg,
     &     ndim_region,idas_region,ndata,idata,iw,nw,
     &     jdata,ierr)
*                                                                      *
*       identify cell of tet and check element IDs                     *
*       Created by T.Furuta on 2019/01/10                              *
*                                                                      *
************************************************************************
      use TETRAMOD, only:nelem,ielem2id,nelemtot
      implicit real*8 (a-h,o-z)

      include 'param.inc'
      include 'err.inc'

*-----------------------------------------------------------------------
      integer,intent(in) :: io,jo,ndim_region,ndata,idata(ndata),iw,nw
      integer,intent(inout) :: mtrn,ntrg(mtrn)
      integer,intent(out) :: ntrn,idas_region(ndim_region)
      integer,intent(out) :: jdata(nelemtot)
      integer,intent(out) :: ierr
*-----------------------------------------------------------------------

      logical iflag,jflag
      integer i,itet,numelem
      integer,allocatable :: iidata(:,:),ielem2ir(:)
      integer,allocatable,save :: ielem2idarray(:,:)

*-----------------------------------------------------------------------

      allocate( ielem2ir(nelemtot) )

      call ttetmesh5(io,jo,ierr,ntrn,mtrn,ntrg,
     &     ndim_region,idas_region,ielem2ir)

      deallocate( ielem2ir )

      if(ierr.eq.1)return

      itet=idas_region(2)
      numelem=nelem(itet)-nelem(itet-1)

      if(ndata.gt.numelem)goto 998

*-----------------------------------------------------------------------

      if(iw.eq.1)then
       allocate( ielem2idarray(numelem,2) )
       do i=1,numelem
        ielem2idarray(i,1)=ielem2id(nelem(itet-1)+i)
        ielem2idarray(i,2)=i
       enddo
       call insertionsort(ielem2idarray,numelem)
      endif

      allocate( iidata(ndata,2) )
      do i=1,ndata
       iidata(i,1)=idata(i)
       iidata(i,2)=i
      enddo
      call insertionsort(iidata,ndata)

*-----------------------------------------------------------------------
*        check element IDs
*-----------------------------------------------------------------------
      jflag=.true.
      k=1
      do i=1,ndata
       iflag=.false.
       do j=k,numelem
        if(iidata(i,1).eq.ielem2idarray(j,1))then
         iflag=.true.
         k=j+1
         exit
        endif
       enddo
       if(.not.iflag)then
        jflag=.false.
        exit
       else
        jdata(ielem2idarray(j,2))=iidata(i,2)
       endif
      enddo
      if(.not.jflag)goto 997

      if(iw.eq.nw)deallocate( ielem2idarray )
      deallocate( iidata )

*-----------------------------------------------------------------------

      return

*-----------------------------------------------------------------------

 997  continue
                  write(io,'(/'' **** Error'',
     &            '' : tetra ID '',i6,'' of '',i2,
     &            ''-th WWG or WWBG is not found in tetmesh '',i2)')
     &             iidata(i,1),iw,itet
                  ErrCha = ''
                  MsgID = 'L:8416/R:ttetmesh4/F:tallsm1.f'
                  call ErrWrite(MsgID, ErrCha)
                  write(jo,'(/'' **** Error'',
     &            '' : tetra ID '',i6,'' of '',i2,
     &            ''-th WWG or WWBG is not found in tetmesh '',i2)')
     &             iidata(i,1),iw,itet
                  goto 999

 998  continue
                  write(io,'(/'' **** Error'',
     &            '' : ndata '',i6,'' is '',
     &            ''more than numelem '',i6)') ndata,numelem
                  ErrCha = ''
                  MsgID = 'L:8429/R:ttetmesh4/F:tallsm1.f'
                  call ErrWrite(MsgID, ErrCha)
                  write(jo,'(/'' **** Error'',
     &            '' : ndata '',i6,'' is '',
     &            ''more than numelem '',i6)') ndata,numelem
                  goto 999

*-----------------------------------------------------------------------

  999 continue

      ierr = 1

      return
      end

************************************************************************
*                                                                      *
      subroutine ttetmesh5(io,jo,ierr,ntrn,mtrn,ntrg,
     &     ndim_region,idas_region,ielem2ir)
*                                                                      *
*       identify cell of tet and compute element volume                *
*       Modified by T.Furuta on 2025/01/16                             *
*                                                                      *
************************************************************************
      use TETRAMOD, only:nelem,tetratvol,ntetcl,itetst,itetlist,
     &     itetelem,nelemtot
      implicit real*8 (a-h,o-z)

      include 'param.inc'
      include 'err.inc'

*-----------------------------------------------------------------------
      integer,intent(in) :: io,jo,ndim_region,ntrg(mtrn)
      integer,intent(inout) :: mtrn
      integer,intent(out) :: ntrn,idas_region(ndim_region)
      integer,intent(out) :: ierr,ielem2ir(nelemtot)
*-----------------------------------------------------------------------
      common /regdc/ idrg(kvlmax), idgr(kvmmax)

      integer nlat3,itetcl(10),itetvol,itetauto,itgchk
      common /tetf1/ nlat3,itetcl,itetvol,itetauto,itgchk

*-----------------------------------------------------------------------

      logical iflag,jflag

*-----------------------------------------------------------------------

      ierr = 0

*-----------------------------------------------------------------------

      iflag=.false.
      ireg=ntrg(1)
      if(ireg.gt.0)then
       idas_region(1)=ireg
       do itet=1,nlat3
        if(ireg.eq.idrg(itetcl(itet)))then
         iflag=.true.
         numelem=nelem(itet)-nelem(itet-1)
         idas_region(2)=itet
         exit
        endif
       enddo
      endif
      if(.not.iflag)then
       if(nlat3.eq.1)then
        iflag=.true.
        itet=1
        ireg=-idrg(itetcl(itet))
        idas_region(1)=ireg
        idas_region(2)=1
        numelem=nelem(itet)-nelem(itet-1)
       endif
      endif
      if(.not.iflag)goto 998
      ielem2ir(1:numelem)=0
      ntrn=ntrg(3)
      if(ntrn.eq.0.and.ireg.lt.0)ntrn=1
      idas_region(3)=ntrn
      jflag=.true.
      n=0
      do i=1,ntrn
       if(ntrn.eq.1.and.ireg.lt.0)then
        icl=idgr(ntrg(1))
        idas_region(3+i)=ntrg(1)
       else
        icl=idgr(ntrg(3+i))
        idas_region(3+i)=ntrg(3+i)
       endif
       do j=1,ntetcl
        if(icl.eq.itetlist(j))then
         iflag=.true.
         exit
        endif
       enddo
       jflag=.false.
       if(iflag)then
        do k=itetst(j)+1,itetst(j+1)
         ielem=itetelem(k)
         if(ielem.gt.nelem(itet-1).and.ielem.le.nelem(itet))then
          n=n+1
          ielem2ir(ielem)=n
          jflag=.true.
         endif ! elems belong to other itet are excluded
        enddo
       endif
       if(.not.jflag)exit
      enddo
      if(.not.jflag)goto 997

      idas_region(4+ntrn)=numelem
      mtrn=4+ntrn
      ntrn=n

      return

*-----------------------------------------------------------------------

 997  continue

      if(ntrn.eq.1.and.ireg.lt.0)then

                  write(io,'(/'' **** Error'',
     &            '' : region '',i6,'' is not'',
     &            '' in TETRA geometry of '',i6,''.'')')
     &            ntrg(1),-ireg
                  ErrCha = ''
                  MsgID = 'L:8558/R:ttetmesh5/F:tallsm1.f'
                  call ErrWrite(MsgID, ErrCha)
                  write(jo,'(/'' **** Error'',
     &            '' : region '',i6,'' is not'',
     &            '' in TETRA geometry of '',i6,''.'')')
     &            ntrg(1),-ireg

      else

                  write(io,'(/'' **** Error'',
     &            '' : region '',i6,'' is not'',
     &            '' in TETRA geometry of '',i6,''.'')')
     &            ntrg(3+i),abs(ireg)
                  ErrCha = ''
                  MsgID = 'L:8572/R:ttetmesh5/F:tallsm1.f'
                  call ErrWrite(MsgID, ErrCha)
                  write(jo,'(/'' **** Error'',
     &            '' : region '',i6,'' is not'',
     &            '' in TETRA geometry of '',i6,''.'')')
     &            ntrg(3+i),abs(ireg)

      endif

                  goto 999
*-----------------------------------------------------------------------

  998 continue
                  write(io,'(/'' **** Error'',
     &            '' : region '',i6,'' is not'',
     &            '' with TETRA geometry.'')') ireg
                  ErrCha = ''
                  MsgID = 'L:8589/R:ttetmesh5/F:tallsm1.f'
                  call ErrWrite(MsgID, ErrCha)
                  write(jo,'(/'' **** Error'',
     &            '' : region '',i6,'' is not'',
     &            '' with TETRA geometry.'')') ireg

                  if(nlat3.gt.1)then

                   write(io,'(/'' **** Error'',
     &            '' : direct TETRA region specification is not'',
     &            '' valid with more than one TETRA geometry. '',
     &            '' Use ( *** *** < *** ) format instead.'')')
                   ErrCha = ''
                   MsgID = 'L:8602/R:ttetmesh5/F:tallsm1.f'
                   call ErrWrite(MsgID, ErrCha)
                   write(jo,'(/'' **** Error'',
     &            '' : direct TETRA region specfication is not'',
     &            '' valid with more than one TETRA geometry. '',
     &            '' Use ( *** *** < *** ) format instead.'')')

                  endif

                  goto 999
*-----------------------------------------------------------------------

  999 continue

      ierr = 1

      return
      end


************************************************************************
*                                                                      *
      subroutine insertionsort(iarray,ndata)
*                                                                      *
*       insertion sort                                                 *
*                                                                      *
************************************************************************
      implicit none
      integer,intent(in) :: ndata
      integer,intent(inout) :: iarray(ndata,2)
      integer :: i,j
      integer :: itmp(2)

      do i=2,ndata
       j=i
       do
        if(j.le.1)exit
        if(iarray(j-1,1).le.iarray(j,1))exit
        itmp(:)=iarray(j-1,:)
        iarray(j-1,:)=iarray(j,:)
        iarray(j,:)=itmp(:)
        j=j-1
       enddo
      enddo
      return
      end

************************************************************************
*                                                                      *
      subroutine txymesh(jsn,jsi,dsin,idsi,ill,ilf,
     &                  jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr,
     &                  ixtp,inx,xmin,xmax,xdel,istxg,
     &                  iytp,iny,ymin,ymax,ydel,istyg,
     &                  iztp,inz,zmin,zmax,zdel,istzg)
*                                                                      *
*       read xyz mesh sub-section of input tally section               *
*       modified by K.Niita on 2000/07/06                              *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      include 'param.inc'
      include 'risrcparam.inc' ! T.Sato 2022/03/25
      include 'err.inc'

*-----------------------------------------------------------------------

      character m_err*200
      common /error/ m_err, l_err, k_err

*-----------------------------------------------------------------------

      character chin*200, chlw*200, chcm*200
      character chlc*200

      character dsin(0:9)*200
      dimension idsi(0:9)

      dimension ill(0:9), ilf(0:9)

      logical iprojall   ! multi-source subsection, proj=all
*-----------------------------------------------------------------------
*     initial values
*-----------------------------------------------------------------------

            ierr  = 0

            ixtp = -1
            iytp = -1
            iztp = -1

*-----------------------------------------------------------------------
*     read one line from jsi
*-----------------------------------------------------------------------

  140 continue

            call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

               if( ierr .ne. 0 ) return
               if( jpn  .eq. 3 ) return

               if( iskip .ne. 0 ) goto 140

  150 continue

               if( ierr .ne. 0 ) return
               if( jpn  .eq. 3 ) return

ccse multi-source subsection, proj=all, write scrach file
            icl=i1
            chlc = chlw
            call chcomp(chlc,icl,i3,i5)
            if( chlc(icl:icl+7) .ne. '<source>' .and.
     &          chlc(icl:icl) .ne. '[' ) then
               inquire(niws,opened=iprojall)
               if(iprojall) write(niws,'(a)') chin(1:i2)
            endif
*-----------------------------------------------------------------------

         if( chcm(i1:i1+6) .eq. 'x-type=' ) then

               ic = inumc(chlw,i1,i3,'=') + 1

               call snum(chlw,ic,i3,ic2,cvvv,ierr)

               if( ierr .ne. 0 ) goto 999

               ixtp = nint( cvvv )

               call getmsh(jsn,jsi,dsin,idsi,ill,ilf,
     &                     jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr,
     &                     'x',ixtp,inx,xmin,xmax,xdel,istxg)

               if( ierr .ne. 0 ) return

               goto 150

         else if( chcm(i1:i1+6) .eq. 'y-type=' ) then

               ic = inumc(chlw,i1,i3,'=') + 1

               call snum(chlw,ic,i3,ic2,cvvv,ierr)

               if( ierr .ne. 0 ) goto 999

               iytp = nint( cvvv )

               call getmsh(jsn,jsi,dsin,idsi,ill,ilf,
     &                     jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr,
     &                     'y',iytp,iny,ymin,ymax,ydel,istyg)

               if( ierr .ne. 0 ) return

               goto 150

         else if( chcm(i1:i1+6) .eq. 'z-type=' ) then

               ic = inumc(chlw,i1,i3,'=') + 1

               call snum(chlw,ic,i3,ic2,cvvv,ierr)

               if( ierr .ne. 0 ) goto 999

               iztp = nint( cvvv )

               call getmsh(jsn,jsi,dsin,idsi,ill,ilf,
     &                     jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr,
     &                     'z',iztp,inz,zmin,zmax,zdel,istzg)

               if( ierr .ne. 0 ) return

               goto 150

         else

               if( ixtp .eq. -1 ) then

                  m_err = 'x-type is missing in xyz mesh'
                  ErrCha = ''
                  ErrID = 'L:8784/R:txymesh/F:tallsm1.f'
                  goto 998

               else if( iytp .eq. -1 ) then

                  m_err = 'y-type is missing in xyz mesh'
                  ErrCha = ''
                  ErrID = 'L:8791/R:txymesh/F:tallsm1.f'
                  goto 998

               else if( iztp .eq. -1 ) then

                  m_err = 'z-type is missing in xyz mesh'
                  ErrCha = ''
                  ErrID = 'L:8798/R:txymesh/F:tallsm1.f'
                  goto 998

               end if

            return

         end if

*-----------------------------------------------------------------------
*     errors
*-----------------------------------------------------------------------

  999 continue

         m_err = 'Description of xyz mesh sub-section is wrong.'
         ErrCha = ''
         ErrID = 'L:8815/R:txymesh/F:tallsm1.f'

  998 continue

         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

      end


************************************************************************
*                                                                      *
      subroutine trzmesh(jsn,jsi,dsin,idsi,ill,ilf,
     &                  jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr,
     &                  rzx0,rzy0,
     &                  irtp,inr,rmin,rmax,rdel,istrg,
     &                  iztp,inz,zmin,zmax,zdel,istzg)
*                                                                      *
*       read r-z mesh sub-section of input tally section               *
*       modified by K.Niita on 2000/07/06                              *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      include 'param.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      character m_err*200
      common /error/ m_err, l_err, k_err

*-----------------------------------------------------------------------

      character chin*200, chlw*200, chcm*200

      character dsin(0:9)*200
      dimension idsi(0:9)

      dimension ill(0:9), ilf(0:9)

*-----------------------------------------------------------------------
*     initial values
*-----------------------------------------------------------------------

            ierr  = 0

            rzx0 = 0.0
            rzy0 = 0.0

            irtp = -1
            iztp = -1

*-----------------------------------------------------------------------
*     read one line from jsi
*-----------------------------------------------------------------------

  140 continue

            call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

               if( ierr .ne. 0 ) return
               if( jpn  .eq. 3 ) return

               if( iskip .ne. 0 ) goto 140

  150 continue

               if( ierr .ne. 0 ) return
               if( jpn  .eq. 3 ) return

*-----------------------------------------------------------------------

         if( chcm(i1:i1+2) .eq. 'x0=' ) then

               ic = inumc(chlw,i1,i3,'=') + 1

               call snum(chlw,ic,i3,ic2,cvvv,ierr)

               if( ierr .ne. 0 ) goto 999

               rzx0 = cvvv

               goto 140

         else if( chcm(i1:i1+2) .eq. 'y0=' ) then

               ic = inumc(chlw,i1,i3,'=') + 1

               call snum(chlw,ic,i3,ic2,cvvv,ierr)

               if( ierr .ne. 0 ) goto 999

               rzy0 = cvvv

               goto 140

         else if( chcm(i1:i1+6) .eq. 'r-type=' ) then

               ic = inumc(chlw,i1,i3,'=') + 1

               call snum(chlw,ic,i3,ic2,cvvv,ierr)

               if( ierr .ne. 0 ) goto 999

               irtp = nint( cvvv )

               call getmsh(jsn,jsi,dsin,idsi,ill,ilf,
     &                     jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr,
     &                     'r',irtp,inr,rmin,rmax,rdel,istrg)

               if( ierr .ne. 0 ) return

               goto 150

         else if( chcm(i1:i1+6) .eq. 'z-type=' ) then

               ic = inumc(chlw,i1,i3,'=') + 1

               call snum(chlw,ic,i3,ic2,cvvv,ierr)

               if( ierr .ne. 0 ) goto 999

               iztp = nint( cvvv )

               call getmsh(jsn,jsi,dsin,idsi,ill,ilf,
     &                     jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr,
     &                     'z',iztp,inz,zmin,zmax,zdel,istzg)

               if( ierr .ne. 0 ) return

               goto 150

         else

               if( irtp .eq. -1 ) then

                  m_err = 'r-type is missing in r-z mesh'
                  ErrCha = ''
                  ErrID = 'L:8960/R:trzmesh/F:tallsm1.f'
                  goto 998

               else if( iztp .eq. -1 ) then

                  m_err = 'z-type is missing in r-z mesh'
                  ErrCha = ''
                  ErrID = 'L:8967/R:trzmesh/F:tallsm1.f'
                  goto 998

               end if

            return

         end if

*-----------------------------------------------------------------------
*     errors
*-----------------------------------------------------------------------

  999 continue

         m_err = 'Description of r-z mesh sub-section is wrong.'
         ErrCha = ''
         ErrID = 'L:8984/R:trzmesh/F:tallsm1.f'

  998 continue

         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

      end


************************************************************************
*                                                                      *
      subroutine tregvol(jsn,jsi,dsin,idsi,ill,ilf,
     &                   jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr,
     &                   nvol,ivl,irvl)
*                                                                      *
*       read 'volume' sub-section of input region tally section        *
*       modified by K.Niita on 2001/11/22                              *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      include 'param.inc'
      include 'risrcparam.inc' ! T.Sato 2022/03/25
      include 'err.inc'

*-----------------------------------------------------------------------

      character m_err*200
      common /error/ m_err, l_err, k_err

*-----------------------------------------------------------------------

      character chin*200, chlw*200, chcm*200
      character chlc*200

      character dsin(0:9)*200
      dimension idsi(0:9)

      dimension ill(0:9), ilf(0:9)
      dimension irsq(10)

      logical dnen2

      dimension idas(mdas*2) !20220906frtati temporary for fbounds-check
      equivalence ( das, idas )

      logical iprojall   ! multi-source subsection, proj=all
*-----------------------------------------------------------------------

            ierr  = 0
            nrsq  = 0
            nvol  = 0

            irsq(1) = 1
            irsq(2) = 2

         do i = 3, 10

            irsq(i) = 0

         end do

            igm = mmmax

*-----------------------------------------------------------------------
*     read one line from jsi
*-----------------------------------------------------------------------

  140 continue

            call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

               if( ierr .ne. 0 ) goto 1000
               if( jpn  .eq. 3 ) goto 1000

               if( iskip .ne. 0 ) goto 140

ccse multi-source subsection, proj=all, write scrach file
ccse multi-source subsection, proj=all, write scrach file
            icl=i1
            chlc = chlw
            call chcomp(chlc,icl,i3,i5)
            if( chlc(icl:icl+7) .ne. '<source>' .and.
     &          chlc(icl:icl) .ne. '[' ) then
               inquire(niws,opened=iprojall)
               if(iprojall) write(niws,'(a)') chin(1:i2)
            endif
*-----------------------------------------------------------------------

         if( nrsq .eq. 0 ) then

               mrsq = 0

               ic = i1

  100       if( ic .gt. i3 ) goto 200

            if(      chlw(ic:ic+2) .eq. 'non' ) then

               irsq( mrsq + 1 ) = 0
               ic = ic + 3

            else if( chlw(ic:ic+2) .eq. 'reg' ) then

               irsq( mrsq + 1 ) = 1
               ic = ic + 3

            else if( chlw(ic:ic+2) .eq. 'vol' .or.
     &               chlw(ic:ic+2) .eq. 'val' ) then

               irsq( mrsq + 1 ) = 2
               ic = ic + 3

            else

               goto 998

            end if

               mrsq = mrsq + 1

               ic = jnumc(chlw,ic,i3)
               goto 100

  200       continue

            if( mrsq .gt. 0 ) then

                  ireg = 0
                  ivol = 0

               do k = 1, mrsq
                  if( irsq(k) .eq. 1 ) ireg = ireg + 1
                  if( irsq(k) .eq. 2 ) ivol = ivol + 1
               end do

                  if( ireg .ne. 1 .or. ivol .ne. 1 ) goto 998

                  nrsq = mrsq

                  goto 140

            else

                  nrsq = 2

            end if

         end if

*-----------------------------------------------------------------------
*        end of this sub-section
*-----------------------------------------------------------------------

         if( dnen2( chlw(i1:i1) ) ) then

            if( ( chlw(i1:i1) .eq. '[' .and. i1 .gt. 5 ) .or.
     &          ( chlw(i1:i1) .eq. '{' ) ) goto 500

            goto 1000

         end if

  500    continue

*-----------------------------------------------------------------------

               ic2  = i1

               nvol = nvol + 1

         do k = 1, nrsq

               ic = jnumc(chlw,ic2,i3)

            if( irsq(k) .eq. 0 ) then

               if( chlw(ic:ic) .eq. '[' ) then

                  ic2 = inumc(chlw,ic,i3,']') + 1

               else if( chlw(ic:ic) .eq. '{' ) then

                  ic2 = inumc(chlw,ic,i3,'}') + 1

               else

                  ic2 = inumc(chlw,ic,i3,' ')

               end if

            else if( irsq(k) .eq. 1 ) then

               call snum(chlw,ic,i3,ic2,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               das( igm + 2 * ( nvol - 1 ) ) = cvvv

            else if( irsq(k) .eq. 2 ) then

               call snum(chlw,ic,i3,ic2,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               das( igm + 2 * ( nvol - 1 ) + 1 ) = cvvv

            end if

         end do

            goto 140

*-----------------------------------------------------------------------
*        summary
*-----------------------------------------------------------------------

 1000 continue

         if( nvol .gt. 0 ) then

            if( igm + nvol * 4 .gt. mdas ) goto 996

               mtmax  = igm + 2 * nvol

            do k = 1, nvol

               das( mtmax + 2 * ( k - 1 ) ) =
     &         das( igm   + 2 * ( k - 1 ) )

               das( mtmax + 2 * ( k - 1 ) + 1 ) =
     &         das( igm   + 2 * ( k - 1 ) + 1 )

            end do

               ivl = ( igm - 1 ) * 2 + 1

            do k = 1, nvol

               idas( ivl + k - 1 ) =
     &         nint( das( mtmax + 2 * ( k - 1 ) ) )

            end do

               mmmax = mmmax + ( nvol + mod(nvol,2) ) / 2

               irvl = mmmax

            do k = 1, nvol

               das( irvl + k - 1 ) = das( mtmax + 2 * ( k - 1 ) + 1 )

            end do

               mmmax = mmmax + nvol

         end if

      return

*-----------------------------------------------------------------------
*     errors
*-----------------------------------------------------------------------

  996 continue

         m_err = 'Total tally storage number exceeds mdas'
         ErrCha = ''
         ErrID = 'L:9259/R:tregvol/F:tallsm1.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  997 continue

         m_err = 'Description of region number or vol/val is wrong.'
         ErrCha = ''
         ErrID = 'L:9271/R:tregvol/F:tallsm1.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  998 continue

         m_err = 'Index of vol/val in region description is wrong.'
         ErrCha = ''
         ErrID = 'L:9283/R:tregvol/F:tallsm1.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

      end


************************************************************************
*                                                                      *
      subroutine tcrsreg(jsn,jsi,dsin,idsi,ill,ilf,
     &                   jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr,
     &                   ntcn,ntcr,ktcr,iseq,mtcr,mmcr
     &                   ,ndim_ntcr,ndim_ktcr,idas_ntcr,das_ktcr)
*                                                                      *
*       read 'r-in r-out area' sub-section of input tally section      *
*       modified by K.Niita on 2010/04/28                              *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      include 'param.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      character m_err*200
      common /error/ m_err, l_err, k_err

*-----------------------------------------------------------------------

      character chin*200, chlw*200, chcm*200

      character dsin(0:9)*200
      dimension idsi(0:9)

      dimension ill(0:9), ilf(0:9)

      integer, intent(in) :: ndim_ntcr
      integer, intent(in) :: ndim_ktcr
      integer, intent(inout) :: idas_ntcr(ndim_ntcr)
      double precision, intent(inout) :: das_ktcr(ndim_ktcr)

*-----------------------------------------------------------------------

      dimension irsq(10)

*-----------------------------------------------------------------------

            ierr  = 0
            nrsq  = 0

            ntcn  = -1
            mtcn  = 0

            irsq(1) = 1
            irsq(2) = 2
            irsq(3) = 3

         do i = 4, 10

            irsq(i) = 0

         end do

*-----------------------------------------------------------------------
*     read one line from jsi
*-----------------------------------------------------------------------

  140 continue

            call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

               if( ierr .ne. 0 ) return
               if( jpn  .eq. 3 ) goto 1000

               if( iskip .ne. 0 ) goto 140

               if( mtcn .eq. ntcn ) goto 1000

*-----------------------------------------------------------------------
*        error not for reg =
*-----------------------------------------------------------------------

         if( ntcn .eq. -1 ) then

            if( chcm(i1:i1+3) .ne. 'reg=' ) goto 999

               ic = inumc(chlw,i1,i3,'=') + 1
               ic = jnumc(chlw,ic,i3)

               call snum(chlw,ic,i3,ic2,cvvv,ierr)

               ntcn = nint( cvvv )

               ktcr  = 1

               goto 140

         end if

*-----------------------------------------------------------------------

         if( nrsq .eq. 0 ) then

               mrsq = 0

               ic = i1

  100       if( ic .gt. i3 ) goto 200

            if(      chlw(ic:ic+2) .eq. 'non' ) then

               irsq( mrsq + 1 ) = 0
               ic = ic + 3

            else if( chlw(ic:ic+3) .eq. 'r-in' ) then

               irsq( mrsq + 1 ) = 1
               ic = ic + 4

            else if( chlw(ic:ic+4) .eq. 'r-out' ) then

               irsq( mrsq + 1 ) = 2
               ic = ic + 5

            else if( chlw(ic:ic+3) .eq. 'area' ) then

               irsq( mrsq + 1 ) = 3
               ic = ic + 4

C S.H. revised for adding r-from and r-to in [t-cross] (2017.9.12)
            else if( chlw(ic:ic+5) .eq. 'r-from' ) then

               irsq( mrsq + 1 ) = 1
               ic = ic + 6

            else if( chlw(ic:ic+3) .eq. 'r-to' ) then

               irsq( mrsq + 1 ) = 2
               ic = ic + 4

C S.H. revised for adding r-from and r-to in [t-cross] (2017.9.12)
            else

               goto 200

            end if

               mrsq = mrsq + 1

               ic = jnumc(chlw,ic,i3)
               goto 100

  200       continue

            if( mrsq .gt. 0 ) then

                  irin  = 0
                  irout = 0
                  iarea = 0

               do k = 1, mrsq
                  if( irsq(k) .eq. 1 ) irin  = irin  + 1
                  if( irsq(k) .eq. 2 ) irout = irout + 1
                  if( irsq(k) .eq. 3 ) iarea = iarea + 1
               end do

                  if( irin .ne. 1 .or. irout .ne. 1 .or.
     &                iarea .ne. 1 ) goto 998

                  nrsq = mrsq

                  goto 140

            else

                  nrsq = 3

            end if

         end if

            iseq = 3

         do k = 1, nrsq

            if( iseq .eq. 3 .and. irsq(k) .eq. 1 ) iseq = 0
            if( iseq .eq. 3 .and. irsq(k) .eq. 2 ) iseq = 1

         end do

*-----------------------------------------------------------------------

               mtcn = mtcn + 1

            if( mtcn .eq. 1 ) then

               ntcr  = 1
                  idsm = ntcr
                  jdsm = -1

            end if

               ic2  = i1
               ntrn = 0

         do k = 1, nrsq

               ic = jnumc(chlw,ic2,i3)

            if( irsq(k) .eq. 0 ) then

               if( chlw(ic:ic) .eq. '[' ) then

                  ic2 = inumc(chlw,ic,i3,']') + 1

               else if( chlw(ic:ic) .eq. '{' ) then

                  ic2 = inumc(chlw,ic,i3,'}') + 1

               else if( chlw(ic:ic) .eq. '(' ) then

                  ic2 = knump(chlw,ic,i3) + 1

               else

                  ic2 = inumc(chlw,ic,i3,' ')

               end if

            else if( irsq(k) .eq. 1 .or. irsq(k) .eq. 2 ) then

                     ndsm = idsm + jdsm + 3

                  call tregion5(chlw,ic,i3,ic2,ntrn,mtrn,ndsm,ierr
     &                          ,ndim_ntcr,idas_ntcr)

                  if( ierr .ne. 0 ) goto 997

                     jdsm = jdsm + 1
                     idas_ntcr(idsm+jdsm) = ntrn

                     jdsm = jdsm + 1
                     idas_ntcr(idsm+jdsm) = mtrn

                     jdsm = jdsm + mtrn

            else if( irsq(k) .eq. 3 ) then

               call snum(chlw,ic,i3,ic2,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               das_ktcr(ktcr-1+mtcn) = cvvv

            end if

         end do

            goto 140

*-----------------------------------------------------------------------
*     summary
*-----------------------------------------------------------------------

 1000 continue

            if( mtcn .ne. ntcn ) goto 995

               mtcr = jdsm


               mmcr = mmmax - ktcr

         return

*-----------------------------------------------------------------------
*     errors
*-----------------------------------------------------------------------

  970 continue

         m_err = 'Memory error: mmmax exceeds mdas '//
     &           ': Please extend mdas in param.inc'
         ErrCha = ''
         ErrID = 'L:9574/R:tcrsreg/F:tallsm1.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  995 continue

         m_err = 'Number of region ( reg = number ) is inconsistent'
         ErrCha = ''
         ErrID = 'L:9586/R:tcrsreg/F:tallsm1.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  997 continue

         m_err = 'Description of region number is wrong.'
         ErrCha = ''
         ErrID = 'L:9598/R:tcrsreg/F:tallsm1.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  998 continue

         m_err = 'Index of crossing region is wrong.'
         ErrCha = ''
         ErrID = 'L:9610/R:tcrsreg/F:tallsm1.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  999 continue

         m_err = 'Number of region ( reg = number ) is missing'
         ErrCha = ''
         ErrID = 'L:9622/R:tcrsreg/F:tallsm1.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

      end


************************************************************************
*                                                                      *
      subroutine tdepreg(icc,jsn,jsi,dsin,idsi,ill,ilf,
     &                   jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr,
     &                   ntrn,mtrn,ndsm,nvol,ivl,irvl,icfl,
     &                   iwgtsum,ntrn1,mtrn1,
     &                   ncond,nadd,ncell,nefflist
     &                   ,ndim_region,idas_region)
*                                                                      *
*       read 'reg =' sub-section of input tally section                *
*                                                                      *
*       modified subroutine 'tregion' to read                          *
*        "no cell operator ethres" sub-section and                     *
*        "cell cond0 cond1 ..." sub-section                            *
*       for 'reg=weightsum' in [T-Deposit]                             *
*                                                                      *
*       modified by S.Abe on 2016/11/24                                *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*       arguments for read input file:                                 *
*           jsn,jsi,dsin,idsi,ill,ilf,                                 *
*           jpn,chin,chlw,chcm,i1,i2,i3,i4                             *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*    add output :                                                      *
*          ntrn      total number of cell                              *
*          mtrn      total number of cell information x ... kr(x)      *
*         ntrn1      number of cell for condition                      *
*         mtrn1      number of cell information for condition          *
*                                                                      *
*         ncond      number of condition                               *
*          nadd      counter for number of additional condition        *
*         ncell      number of cell in efficiency list                 *
*      nefflist      number of "list" in index                         *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     local variable :                                                 *
*           ign      start position of das memory                      *
*         mcond      counter for number of condition                   *
*          nrsq      counter for number of reading index               *
*       irsq(5)      array for  reading index                          *
*                                                                      *
*         mcell      counter for number of efficiency list             *
*         nrsq2      counter for number of reading index               *
*        noflag      position of "cell" in index                       *
*       ncount1      counter for number of reading efficiency line     *
*       ncount2      counter for number of reading efficiency value    *
*         ntrn2      number of cell for efficiency list                *
*         mtrn2      number of cell information for efficiency list    *
*                                                                      *
************************************************************************
      use moddas
      use moddas_region
      use moddas_region_mtrg

      use tdepwgtsum_local

      implicit real*8 (a-h,o-z)

      include 'param.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      character m_err*200
      common /error/ m_err, l_err, k_err

*-----------------------------------------------------------------------

      character chin*200, chlw*200, chcm*200

      character dsin(0:9)*200
      dimension idsi(0:9)

      dimension ill(0:9), ilf(0:9)

*-----------------------------------------------------------------------

      dimension ipar(0:20)
      dimension jpar(0:20)
      dimension lpar(0:20)

      data klnmax /1000000/

      integer, intent(in) :: ndim_region
      integer, intent(inout) :: idas_region(ndim_region)

      dimension irsq(5)

*-----------------------------------------------------------------------

      ierr = 0

      if( icfl .eq. 1 ) goto 150

*-----------------------------------------------------------------------
*     read one line from jsi
*-----------------------------------------------------------------------

  140 continue

      call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &           jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

      if( ierr .ne. 0 ) return
      if( jpn  .eq. 3 ) return

      if( iskip .ne. 0 ) goto 140

  150 continue

      if( ierr .ne. 0 ) return
      if( jpn  .eq. 3 ) return

*-----------------------------------------------------------------------
*     error not for reg =
*-----------------------------------------------------------------------

      if( icc .eq. 1 ) then

         if( chcm(i1:i1+3) .ne. 'reg=' ) goto 998

      else if( icc .eq. 2 ) then

         if( chcm(i1:i1+8) .ne. 'reginbox=' ) goto 998

      end if

      ic = inumc(chlw,i1,i3,'=') + 1

*-----------------------------------------------------------------------
*     read region number for weighted summation
*-----------------------------------------------------------------------

      ic = jnumc(chlw,ic,i3)

      if( chlw(ic:ic+8) .eq. 'weightsum' ) then

         iwgtsum = 1
         igm = 1

         nvol = 0   ! S.Abe 2018/11/15

*-----------------------------------------------------------------------
*        initialize
*-----------------------------------------------------------------------

         ntrn = 0
         mtrn = 0
         ntrn1 = 0
         mtrn1 = 0
         ntrn2 = 0
         mtrn2 = 0

         ign = igm

         ncond = -1
         mcond = 0
         nadd = 1
         madd = 1
         nrsq = 0
         do i = 1, 5
            irsq(i) = i
         enddo

         ncell = -1
         mcell = 0
         nrsq2 = 0
         noflag = 0
         nefflist = 0
         ncount1 = 0

*-----------------------------------------------------------------------
*        read one line for conditon section
*-----------------------------------------------------------------------

  100    continue
         call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &              jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)
         if( iskip .ne. 0 ) goto 100

*-----------------------------------------------------------------------
*        read "ncond=" value
*-----------------------------------------------------------------------

         if( ncond .lt. 0 ) then

            if( ierr .ne. 0 ) goto 961
            if( jpn  .eq. 3 ) goto 961

            if( chcm(i1:i1+5) .ne. 'ncond=' ) goto 961

            ic = inumc(chlw,i1,i3,'=') + 1
            ic = jnumc(chlw,ic,i3)

            call snum(chlw,ic,i3,ic2,cvvv,ierr)
            if( ierr .ne. 0 ) goto 962

            ncond = nint( cvvv )
            if( ncond .lt. 0 ) goto 962
            goto 100

         endif

*-----------------------------------------------------------------------
*        read index of condition
*-----------------------------------------------------------------------

         if( nrsq .eq. 0 ) then

            if( ierr .ne. 0 ) goto 963
            if( jpn  .eq. 3 ) goto 963

            ic = i1

  110       if( ic .gt. i3 ) goto 120

            if(     chlw(ic:ic+1) .eq. 'no' ) then
               irsq(nrsq+1) = 1
               ic = ic + 2
            elseif( chlw(ic:ic+3) .eq. 'cell' ) then
               irsq(nrsq+1) = 2
               ic = ic + 4
            elseif( chlw(ic:ic+7) .eq. 'operator' ) then
               irsq(nrsq+1) = 3
               ic = ic + 8
            elseif( chlw(ic:ic+6) .eq. 'ethres' ) then
               irsq(nrsq+1) = 4
               ic = ic + 6
            elseif( chlw(ic:ic+6) .eq. 'list' ) then
               irsq(nrsq+1) = 5
               ic = ic + 4
            else
               goto 120
            endif

            nrsq = nrsq + 1
            if( nrsq .gt. 5 ) goto 963

            ic = jnumc(chlw,ic,i3)
            goto 110

*-----------------------------------------------------------------------
*        check index of condition
*-----------------------------------------------------------------------

  120       continue

            if( nrsq .gt. 0 ) then

               irno = 0
               ircel = 0
               irope = 0
               ireth = 0
               irlist = 0

               do k = 1, nrsq
                  if( irsq(k) .eq. 1 ) irno   = irno + 1
                  if( irsq(k) .eq. 2 ) ircel  = ircel + 1
                  if( irsq(k) .eq. 3 ) irope  = irope + 1
                  if( irsq(k) .eq. 4 ) ireth  = ireth + 1
                  if( irsq(k) .eq. 5 ) irlist = irlist + 1
               enddo

               if( irno .ne. 1 .or. ircel .ne. 1 .or.
     &             irope .ne. 1 .or. ireth .ne. 1 .or.
     &             irlist .ne. 1 ) goto 963

               call allocate_depwgtsum01(ncond,nadd)

               goto 100

            else

               goto 963

            endif

         endif

*-----------------------------------------------------------------------
*        read condition parameter
*-----------------------------------------------------------------------

         if( ierr .ne. 0 ) goto 964
         if( jpn  .eq. 3 ) goto 964

         ic = i1
         ic2 = i1
         ntrn = 0

         do k = 1, nrsq

            ic = jnumc(chlw,ic2,i3)

            if( irsq(k) .eq. 1 ) then   ! no

               if( chlw(ic:ic+2) .eq. 'and' ) then

                  if( mcond .eq. 0 ) goto 963

                  itmp1 = 1
                  madd = madd + 1
                  if( madd .gt. nadd ) then
                     call reallocate_depwgtsum01(ncond,nadd,madd)
                     nadd = madd
                  endif

                  ic = ic + 3
                  ic2 = ic + 1

               elseif( chlw(ic:ic+4) .eq. 'ncell' ) then

                  if( mcond .ne. ncond ) goto 964
                  goto 201   ! goto next index section

               else

                  call snum(chlw,ic,i3,ic2,cvvv,ierr)
                  if( ierr .ne. 0 ) goto 964

                  itmp1 = nint( cvvv )
                  if( itmp1 .eq. 0 ) goto 965

                  do im = 1, mcond
                     if( itmp1 .eq. icond(1,im,1) ) goto 966
                  enddo

                  mcond = mcond + 1
                  if( mcond .gt. ncond ) goto 964

                  madd = 1

               endif

            elseif( irsq(k) .eq. 2 ) then   ! cell

               call tregion5(chlw,ic,i3,ic2,ntrn,mtrn,ign,ierr
     &                       ,ndim_region,idas_region)
               if( ierr .ne. 0 ) goto 964

               ign = ign + mtrn
               ntrn1 = ntrn1 + ntrn
               mtrn1 = mtrn1 + mtrn

            elseif( irsq(k) .eq. 3 ) then   ! operator

               if(     chlw(ic:ic+1) .eq. 'le' .or.
     &                 chlw(ic:ic+1) .eq. '<=' .or.
     &                 chlw(ic:ic+1) .eq. '=<' ) then
                  itmp2 = 2
                  ic = ic + 2
                  ic2 = ic + 1
               elseif( chlw(ic:ic+1) .eq. 'ge' .or.
     &                 chlw(ic:ic+1) .eq. '>=' .or.
     &                 chlw(ic:ic+1) .eq. '=>' ) then
                  itmp2 = 4
                  ic = ic + 2
                  ic2 = ic + 1
               elseif( chlw(ic:ic+1) .eq. 'lt' .or.
     &                 chlw(ic:ic)   .eq. '<' ) then
                  itmp2 = 1
                  ic = ic + 2
                  ic2 = ic + 1
               elseif( chlw(ic:ic+1) .eq. 'gt' .or.
     &                 chlw(ic:ic)   .eq. '>' ) then
                  itmp2 = 5
                  ic = ic + 2
                  ic2 = ic + 1
               elseif( chlw(ic:ic+1) .eq. 'eq' .or.
     &                 chlw(ic:ic)   .eq. '=' ) then
                  itmp2 = 3
                  ic = ic + 2
                  ic2 = ic + 1
               else
                  goto 964
               endif

            elseif( irsq(k) .eq. 4 ) then   ! ethres

               call snum(chlw,ic,i3,ic2,cvvv,ierr)
               if( ierr .ne. 0 ) goto 964

               tmp3 = dble( cvvv )
               if( tmp3 .lt. 0.d0 ) goto 967

            elseif( irsq(k) .eq. 5 ) then   ! list

               call snum(chlw,ic,i3,ic2,cvvv,ierr)
               if( ierr .ne. 0 ) goto 964

               itmp4 = nint( cvvv )
               if( itmp4 .eq. 0 ) goto 968

            else

               goto 963

            endif

         enddo

         icond(1,mcond,madd) = itmp1
         icond(2,mcond,madd) = itmp2
         icond(3,mcond,madd) = itmp4
         ethres(mcond,madd) = tmp3

         goto 100

*-----------------------------------------------------------------------
*        read one line for efficiency section
*-----------------------------------------------------------------------

  200    continue
         call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &              jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)
         if( iskip .ne. 0 ) goto 200

*-----------------------------------------------------------------------
*        read "ncell=" value
*-----------------------------------------------------------------------

  201    continue

         if( ncell .lt. 0 ) then

            if( ierr .ne. 0 ) goto 964
            if( jpn  .eq. 3 ) goto 964

            if( chcm(i1:i1+5) .ne. 'ncell=' ) goto 971

            ic = inumc(chlw,i1,i3,'=') + 1
            ic = jnumc(chlw,ic,i3)

            call snum(chlw,ic,i3,ic2,cvvv,ierr)
            if( ierr .ne. 0 ) goto 972

            ncell = nint( cvvv )
            if( ncell .le. 0 ) goto 972

            goto 200

         endif

*-----------------------------------------------------------------------
*        read index of efficiency
*-----------------------------------------------------------------------

         if( nefflist .eq. 0 ) then

            if( ierr .ne. 0 ) goto 973
            if( jpn  .eq. 3 ) goto 973

            ic = i1

  210       if( ic .gt. i3 ) goto 220

            if(     chlw(ic:ic+3) .eq. 'cell' ) then
               ic = ic + 4
               ic2 = ic + 1
               if( noflag .ne. 0 ) goto 973
               noflag = nrsq2 + 1
            elseif( chlw(ic:ic+3) .eq. 'list' ) then
               ic = ic + 4
               ic2 = ic + 1
               call snum(chlw,ic,i3,ic2,cvvv,ierr)
               if( ierr .ne. 0 ) goto 973
               ic = ic2 - 1
               nefflist = nefflist + 1
            else
               goto 973
            endif

            nrsq2 = nrsq2 + 1

            ic = jnumc(chlw,ic,i3)
            goto 210

*-----------------------------------------------------------------------

  220       continue

            call allocate_depwgtsum02(ncell,nefflist)

            nefflist = 0
            ic = i1

  230       if( ic .gt. i3 ) goto 240

            if(     chlw(ic:ic+3) .eq. 'cell' ) then
               ic = ic + 4
               ic2 = ic + 1
            elseif( chlw(ic:ic+3) .eq. 'list' ) then
               ic = ic + 4
               ic2 = ic + 1
               call snum(chlw,ic,i3,ic2,cvvv,ierr)
               if( ierr .ne. 0 ) goto 973
               ic = ic2 - 1
               nefflist = nefflist + 1
               numlist(nefflist) = nint( cvvv )
               do i = 1, nefflist-1
                  if( numlist(nefflist) .eq. numlist(i) ) goto 973
               enddo
            else
               goto 973
            endif

            ic = jnumc(chlw,ic,i3)
            goto 230

*-----------------------------------------------------------------------
*        check correspondence of list number in condition
*-----------------------------------------------------------------------

  240       continue

            if( noflag .eq. 0 ) goto 973
            if( nefflist .le. 0 ) goto 973

            do i = 1, mcond

               nagree = 0

               do j = 1, nefflist
                  if( numlist(j) .eq. 0 ) cycle
                  if( icond(3,i,1) .eq. numlist(j) ) then
                     nagree = 1
                     exit
                  endif
               enddo

               if( nagree .eq. 0 ) goto 974

            enddo

            goto 200

         endif

*-----------------------------------------------------------------------
*        read efficiency parameter
*-----------------------------------------------------------------------

         if( ncount1 .eq. ncell ) goto 300

         if( ierr .ne. 0 ) goto 975
         if( jpn  .eq. 3 ) goto 975

         ic = i1
         ic2 = i1
         ntrn = 0

         ncount1 = ncount1 + 1
         ncount2 = 0

         do k = 1, nefflist + 1

            ic = jnumc(chlw,ic2,i3)

            if( k .eq. noflag ) then

               call tregion5(chlw,ic,i3,ic2,ntrn,mtrn,ign,ierr
     &                       ,ndim_region,idas_region)
               if( ierr .ne. 0 ) goto 975

               ign = ign + mtrn
               ntrn2 = ntrn2 + ntrn
               mtrn2 = mtrn2 + mtrn

            else

               call snum(chlw,ic,i3,ic2,cvvv,ierr)
               if( ierr .ne. 0 ) goto 975

               ncount2 = ncount2 + 1
               depeff(ncount1,ncount2) = cvvv

            endif

         enddo

         goto 200

*-----------------------------------------------------------------------

  300    continue

         ntrn = ntrn1 + ntrn2
         mtrn = mtrn1 + mtrn2

*-----------------------------------------------------------------------
*     normal case
*-----------------------------------------------------------------------

      else

         iwgtsum = 0

         ntrn = 0
         mtrn = 0
         ibra = 0
         ilat = 0
         klat = 0
         ifis = 0
         kpar = 0
         iuni = 0

         nvol = 0

         do i = 0, 20
            ipar(i) = 0
            jpar(i) = 0
            lpar(i) = 0
         end do

         igm = 0
         call moddas_allocate_int2(0, 20, MAX_NUM_MTRG, mtrg)
         ipmax = 0

         do i = 1, klnmax

            ic = jnumc(chlw,ic,i3)

            call tregion3(icn,chlw,i1,i3,ic,
     &                    igm,ipmax,ipar,jpar,
     &                    ibra,ilat,klat,ifis,kpar,iuni
     &                    ,MAX_NUM_MTRG,mtrg)
            if( icn .gt. 900 ) goto 900
            if( icn .eq. 500 ) goto 500
            if( jpn .eq. 3 ) return

            if( i .lt. klnmax .and. ic .gt. i3 ) then

  147          call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)
               if( ierr .ne. 0 ) return
               if( jpn  .eq. 3 ) goto 500

               if( iskip .ne. 0 ) goto 147

               ifis = ifis + 1

               ic = i1

            end if

         end do

         goto 996

*-----------------------------------------------------------------------
*        modify the input
*-----------------------------------------------------------------------

  500    continue

         ntrn = ipar(0)
         mtrn = ipar(0)

         do k = 1, mtrn
            idas_region( ndsm + k - 1 ) = mtrg(igm+0,k)
         end do

         call tregion4(ierr,ntrn,mtrn,idas_region(ndsm),ndsm,igm
     &                 ,ndim_region,idas_region)
         if( ierr .ne. 0 ) goto 951

      endif

*-----------------------------------------------------------------------


*-----------------------------------------------------------------------
*     volume
*-----------------------------------------------------------------------

      if( chcm(i1:i1+5) .eq. 'volume' .or.
     &    chcm(i1:i1+4) .eq. 'value' ) then

         call tregvol(jsn,jsi,dsin,idsi,ill,ilf,
     &                jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr,
     &                nvol,ivl,irvl)

      end if

*-----------------------------------------------------------------------

      return

*-----------------------------------------------------------------------
*     error messages
*-----------------------------------------------------------------------

  950 continue
      m_err = 'Memory error: mmmax exceeds mdas '//
     &        ': Please extend mdas in param.inc'
      ErrCha = ''
      ErrID = 'L:10333/R:tdepreg/F:tallsm1.f'
      l_err = ill(jsn)
      k_err = jsn
      ierr = 1
      return
  951 continue
      m_err = 'region is too many or memory is lack'
      ErrCha = ''
      ErrID = 'L:10341/R:tdepreg/F:tallsm1.f'
      l_err = ill(jsn)
      k_err = jsn
      ierr  = 1
      return

  961 continue
      m_err = 'After {reg = weightsum} line should be {ncond =}.'
      ErrCha = ''
      ErrID = 'L:10350/R:tdepreg/F:tallsm1.f'
      l_err = ill(jsn)
      k_err = jsn
      ierr = 1
      return
  962 continue
      m_err = 'Number of condition is wrong.'
      ErrCha = ''
      ErrID = 'L:10358/R:tdepreg/F:tallsm1.f'
      l_err = ill(jsn)
      k_err = jsn
      ierr = 1
      return
  963 continue
      m_err = 'Index of condition is wrong.'
      ErrCha = ''
      ErrID = 'L:10366/R:tdepreg/F:tallsm1.f'
      l_err = ill(jsn)
      k_err = jsn
      ierr = 1
      return
  964 continue
      m_err = 'Input value in condition is wrong.'
      ErrCha = ''
      ErrID = 'L:10374/R:tdepreg/F:tallsm1.f'
      l_err = ill(jsn)
      k_err = jsn
      ierr = 1
      return
  965 continue
      m_err = 'Zero is not allowed for condition number.'
      ErrCha = ''
      ErrID = 'L:10382/R:tdepreg/F:tallsm1.f'
      l_err = ill(jsn)
      k_err = jsn
      ierr = 1
      return
  966 continue
      m_err = 'Condition number is duplicated.'
      ErrCha = ''
      ErrID = 'L:10390/R:tdepreg/F:tallsm1.f'
      l_err = ill(jsn)
      k_err = jsn
      ierr = 1
      return
  967 continue
      m_err = 'Threshold energy is less than zero.'
      ErrCha = ''
      ErrID = 'L:10398/R:tdepreg/F:tallsm1.f'
      l_err = ill(jsn)
      k_err = jsn
      ierr = 1
      return
  968 continue
      m_err = 'Zero is not allowded for list number.'
      ErrCha = ''
      ErrID = 'L:10406/R:tdepreg/F:tallsm1.f'
      l_err = ill(jsn)
      k_err = jsn
      ierr = 1
      return

  971 continue
      m_err = 'After condition of weightsum should be {ncell =}.'
      ErrCha = ''
      ErrID = 'L:10415/R:tdepreg/F:tallsm1.f'
      l_err = ill(jsn)
      k_err = jsn
      ierr = 1
      return
  972 continue
      m_err = 'Number of ncell is wrong.'
      ErrCha = ''
      ErrID = 'L:10423/R:tdepreg/F:tallsm1.f'
      l_err = ill(jsn)
      k_err = jsn
      ierr = 1
      return
  973 continue
      m_err = 'Index of ncell is wrong.'
      ErrCha = ''
      ErrID = 'L:10431/R:tdepreg/F:tallsm1.f'
      l_err = ill(jsn)
      k_err = jsn
      ierr = 1
      return
  974 continue
      m_err = 'List number defined in condition does not exist.'
      ErrCha = ''
      ErrID = 'L:10439/R:tdepreg/F:tallsm1.f'
      l_err = ill(jsn)
      k_err = jsn
      ierr = 1
      return
  975 continue
      m_err = 'Input value in ncell is wrong.'
      ErrCha = ''
      ErrID = 'L:10447/R:tdepreg/F:tallsm1.f'
      l_err = ill(jsn)
      k_err = jsn
      ierr = 1
      return

*-----------------------------------------------------------------------
*     errors for tregion3
*-----------------------------------------------------------------------

  900 continue

         if( icn .eq. 980 ) goto 980
         if( icn .eq. 984 ) goto 984
         if( icn .eq. 985 ) goto 985
         if( icn .eq. 986 ) goto 986
         if( icn .eq. 987 ) goto 987
         if( icn .eq. 988 ) goto 988
         if( icn .eq. 989 ) goto 989
         if( icn .eq. 990 ) goto 990
         if( icn .eq. 991 ) goto 991
         if( icn .eq. 992 ) goto 992
         if( icn .eq. 993 ) goto 993
         if( icn .eq. 995 ) goto 995
         if( icn .eq. 997 ) goto 997
         if( icn .eq. 999 ) goto 999

*-----------------------------------------------------------------------

  980 continue

         m_err = 'n1-n2 should be used like {n1-n2}'
         ErrCha = ''
         ErrID = 'L:10480/R:tdepreg/F:tallsm1.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  984 continue

         m_err = 'maximum lattice elements by commas is 1000.'
         ErrCha = ''
         ErrID = 'L:10492/R:tdepreg/F:tallsm1.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  985 continue

         m_err = 'maximum level of < is 10.'
         ErrCha = ''
         ErrID = 'L:10504/R:tdepreg/F:tallsm1.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  986 continue

         m_err = 'usage of u=# is wrong.'
         ErrCha = ''
         ErrID = 'L:10516/R:tdepreg/F:tallsm1.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  987 continue

         m_err = '< should be inside ( ).'
         ErrCha = ''
         ErrID = 'L:10528/R:tdepreg/F:tallsm1.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  988 continue

         m_err = 'maximun level of ( ) is 10.'
         ErrCha = ''
         ErrID = 'L:10540/R:tdepreg/F:tallsm1.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  989 continue

         m_err = 'usage of latice [i1:i2 i3:i4 i5:i6] is wrong.'
         ErrCha = ''
         ErrID = 'L:10552/R:tdepreg/F:tallsm1.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  990 continue

         m_err = 'usage of latice [i1 i2 i3] is wrong.'
         ErrCha = ''
         ErrID = 'L:10564/R:tdepreg/F:tallsm1.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  991 continue

         m_err = 'all or (all) is available. (all<4), (all 3) are not.'
         ErrCha = ''
         ErrID = 'L:10576/R:tdepreg/F:tallsm1.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  992 continue

         m_err = 'Usage of Parenthesis { n1 - n2 } is wrong'
         ErrCha = ''
         ErrID = 'L:10588/R:tdepreg/F:tallsm1.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  993 continue

         m_err = 'Usage of Parenthesis ( ) is wrong'
         ErrCha = ''
         ErrID = 'L:10600/R:tdepreg/F:tallsm1.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  995 continue

         m_err = 'Description of region number '//
     &           '{n1-n2} (n1<n2) is wrong.'
         ErrCha = ''
         ErrID = 'L:10613/R:tdepreg/F:tallsm1.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  996 continue

         m_err = 'Number of region is too large. max(klnmax)=1000000'
         ErrCha = ''
         ErrID = 'L:10625/R:tdepreg/F:tallsm1.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  997 continue

         m_err = 'Description of region number is wrong.'
         ErrCha = ''
         ErrID = 'L:10637/R:tdepreg/F:tallsm1.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  998 continue

         m_err = 'After { mesh = reg } line should be { reg = }.'
         ErrCha = ''
         ErrID = 'L:10649/R:tdepreg/F:tallsm1.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  999 continue
         l_err = ill(jsn)
         k_err = jsn
         ierr = 1
         return

*-----------------------------------------------------------------------

      end


************************************************************************
*                                                                      *
      subroutine tdpalib(jsn,jsi,dsin,idsi,ill,ilf,
     &                   jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr,
     &                   icnp,ienp,dlmax,jtln,jtli,jtlr)
*                                                                      *
*       read dpa library information                                   *
*       modified by K.Niita on 2002/12/13                              *
*                                                                      *
************************************************************************
      use moddas
      use moddas_tally

      implicit real*8 (a-h,o-z)

      include 'param.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      character m_err*200
      common /error/ m_err, l_err, k_err

*-----------------------------------------------------------------------

      character chin*200, chlw*200, chcm*200

      character dsin(0:9)*200
      dimension idsi(0:9)

      dimension ill(0:9), ilf(0:9)


*-----------------------------------------------------------------------

      dimension irsq(10)

*-----------------------------------------------------------------------

            ierr  = 0
            nrsq  = 0
            icnp  = 0
            ienp  = 0
            dlmax = -1.0d0

            ntcn  = jtln
            mtcn  = 0
            mrsq  = 0


            irsq(1) = 1
            irsq(2) = 4
            irsq(3) = 2
            irsq(4) = 3

         do i = 5, 10

            irsq(i) = 0

         end do

*-----------------------------------------------------------------------

            jtli = (iaddress_mlib(icurrent_mlib) - 1)*3 + 1
            call moddas_reallocate_int2(
     &              3, icurrent_mlib, 1, 3, jtln, iaddress_mlib, mlib)

            jtlr = iaddress_flib(icurrent_flib)
            call moddas_reallocate_dbl(
     &              3, icurrent_flib, jtln, iaddress_flib, flib)

            iaddress_mlib(2) = iaddress_mlib(3)
            icurrent_mlib = 2
            iaddress_flib(2) = iaddress_flib(3)
            icurrent_flib = 2


            if( mmmax .gt. mdas ) goto 970

*-----------------------------------------------------------------------
*     read one line from jsi
*-----------------------------------------------------------------------

  140 continue

            call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

               if( ierr .ne. 0 ) return
               if( jpn  .eq. 3 ) goto 1000

               if( iskip .ne. 0 ) goto 140

               if( mtcn .eq. ntcn ) goto 1000

*-----------------------------------------------------------------------
*        read particle and max energy
*-----------------------------------------------------------------------

         if( nrsq .eq. 0 .and. mrsq .eq. 0 ) then

               ic = i1

            if( chlw(ic:ic+3) .eq. 'part' ) then

               if( ienp .ne. 0 ) goto 993
               if( icnp .ne. 0 ) goto 993

               ic = inumc(chlw,ic+3,i3,'=') + 1
               ic = jnumc(chlw,ic,i3)

               if( ic .gt. i3 ) goto 994

               if( chlw(ic:ic+5) .eq. 'proton' ) then

                  icnp = 1

               else if( chlw(ic:ic+6) .eq. 'neutron' ) then

                  icnp = 2

               else

                  goto 994

               end if

                  goto 140

            end if

            if( chlw(ic:ic+3) .eq. 'emax' ) then

               if( ienp .ne. 0 ) goto 993
               if( icnp .eq. 0 ) icnp = 2
               ienp = ienp + 1

               ic = inumc(chlw,ic+3,i3,'=') + 1
               ic = jnumc(chlw,ic,i3)

               if( ic .gt. i3 ) goto 992

               call snum(chlw,ic,i3,ic2,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               dlmax = cvvv

                  goto 140

            end if

         end if

*-----------------------------------------------------------------------
*        read data definition
*-----------------------------------------------------------------------

         if( nrsq .eq. 0 ) then

               ic = i1

  100       if( ic .gt. i3 ) goto 200

            if(      chlw(ic:ic+2) .eq. 'non' ) then

               irsq( mrsq + 1 ) = 0
               ic = ic + 3

            else if( chlw(ic:ic+2) .eq. 'mat' ) then

               irsq( mrsq + 1 ) = 1
               ic = ic + 3

            else if( chlw(ic:ic+2) .eq. 'fac' ) then

               irsq( mrsq + 1 ) = 4
               ic = ic + 3

            else if( chlw(ic:ic+2) .eq. 'lib' ) then

               irsq( mrsq + 1 ) = 2
               ic = ic + 3

            else if( chlw(ic:ic+1) .eq. 'mt' ) then

               irsq( mrsq + 1 ) = 3
               ic = ic + 2

            else

               goto 998

            end if

               mrsq = mrsq + 1

               ic = jnumc(chlw,ic,i3)
               goto 100

  200       continue

            if( mrsq .gt. 0 ) then

                  imat = 0
                  ifac = 0
                  ilib = 0
                  imtt = 0

               do k = 1, mrsq
                  if( irsq(k) .eq. 1 ) imat = imat + 1
                  if( irsq(k) .eq. 4 ) ifac = ifac + 1
                  if( irsq(k) .eq. 2 ) ilib = ilib + 1
                  if( irsq(k) .eq. 3 ) imtt = imtt + 1
               end do

                  if( imat .ne. 1 .or. ilib .ne. 1 ) goto 998

                  nrsq = mrsq

                  goto 140

            else

                  nrsq = 4
                  mrsq = 4

            end if

         end if

*-----------------------------------------------------------------------

               mtcn = mtcn + 1
               ic2  = i1
               ntrn = 0

               mlib(jtli+2,mtcn) = 444
               flib(jtlr-1+mtcn) = 1.0

         do k = 1, nrsq

               ic = jnumc(chlw,ic2,i3)

            if( irsq(k) .eq. 0 ) then

               if( chlw(ic:ic) .eq. '[' ) then

                  ic2 = inumc(chlw,ic,i3,']') + 1

               else if( chlw(ic:ic) .eq. '{' ) then

                  ic2 = inumc(chlw,ic,i3,'}') + 1

               else if( chlw(ic:ic) .eq. '(' ) then

                  ic2 = knump(chlw,ic,i3) + 1

               else

                  ic2 = inumc(chlw,ic,i3,' ')

               end if

            else if( irsq(k) .eq. 4 ) then

               call snum(chlw,ic,i3,ic2,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               flib(jtlr-1+mtcn) = cvvv

            else if( irsq(k) .le. 3 ) then

               call snum(chlw,ic,i3,ic2,cvvv,ierr)

               if( ierr .ne. 0  ) goto 997
               if( cvvv .lt. 0. ) goto 997

               mlib(jtli-1+irsq(k),mtcn) = nint( cvvv )

            end if

         end do

            goto 140

*-----------------------------------------------------------------------
*     summary
*-----------------------------------------------------------------------

 1000 continue

            if( mtcn .ne. ntcn ) goto 995

         return

*-----------------------------------------------------------------------
*     errors
*-----------------------------------------------------------------------

  970 continue

         m_err = 'Memory error: mmmax exceeds mdas '//
     &           ': Please extend mdas in param.inc'
         ErrCha = ''
         ErrID = 'L:10974/R:tdpalib/F:tallsm1.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  992 continue

         m_err = 'emax =  description is wrong'
         ErrCha = ''
         ErrID = 'L:10986/R:tdpalib/F:tallsm1.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  993 continue

         m_err = 'part = should be before emax =.'
         ErrCha = ''
         ErrID = 'L:10998/R:tdpalib/F:tallsm1.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  994 continue

         m_err = 'part = neutron or proton : something wrong'
         ErrCha = ''
         ErrID = 'L:11010/R:tdpalib/F:tallsm1.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  995 continue

         m_err = 'Number of library is inconsistent with the # of data'
         ErrCha = ''
         ErrID = 'L:11022/R:tdpalib/F:tallsm1.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  997 continue

         m_err = 'Description of library data is wrong.'
         ErrCha = ''
         ErrID = 'L:11034/R:tdpalib/F:tallsm1.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  998 continue

         m_err = 'Index of dpa library is wrong.'
         ErrCha = ''
         ErrID = 'L:11046/R:tdpalib/F:tallsm1.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

      end


************************************************************************
*                                                                      *
      subroutine tmultipl(jsn,jsi,dsin,idsi,ill,ilf,
     &                   jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr,
     &                   jtln,jtal,inpat,iptyp,ipnkf,dlmax,
     &                   jtli,itls,mmst,imst,imtinf) ! frtati 2023/12/07
*                                                                      *
*       read multiplier information                                    *
*       modified by K.Niita on 2002/12/23                              *
*                                                                      *
************************************************************************
      use moddas
      use moddas_tally

      implicit real*8 (a-h,o-z)

      include 'param.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      character m_err*200
      common /error/ m_err, l_err, k_err

*-----------------------------------------------------------------------

      character chin*200, chlw*200, chcm*200

      character dsin(0:9)*200
      dimension idsi(0:9)

      dimension ill(0:9), ilf(0:9)


*-----------------------------------------------------------------------

      dimension iptyp(6), ipnkf(6)
      dimension jstyp(6), jnkf0(6)

      dimension irsq(100)
      dimension imst(6)

      dimension imtinf(4) ! frtati 2023/12/07

      logical dnen3
      logical dcom3 ! T.Sato 2021/02/12

*-----------------------------------------------------------------------

            ierr  = 0
            inpat = 0
            ienp  = 0
            icnp  = 0

            nrsq  = 0
            mrsq  = 0
            mmst  = 0

            dlmax = 1.0d+10

            ntcn  = jtln
            mtcn  = 0

            irsq(1) = 1
            irsq(2) = 2

            imtinf = 0 ! frtati 2023/12/07
            iemtp = 0  ! frtati 2023/12/07

*-----------------------------------------------------------------------

            jtli = (iaddress_mltp(icurrent_mltp) - 1)*13 + 1
            call moddas_reallocate_int2(
     &             3, icurrent_mltp, 1, 13, jtln, iaddress_mltp, mltp)

            if( mmmax .gt. mdas ) goto 970

            jtls = iaddress_slib(icurrent_slib)
            call moddas_reallocate_dbl(
     &              3, icurrent_slib, MAX_NUM_SLIB, iaddress_slib, slib)
            itls = 0

*-----------------------------------------------------------------------
*     read one line from jsi
*-----------------------------------------------------------------------

  140 continue

            call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

               if( ierr .ne. 0 ) return
               if( jpn  .eq. 3 ) goto 1000
               if( iskip .ne. 0 ) goto 140

               if( mtcn .eq. ntcn ) goto 1000

*-----------------------------------------------------------------------
*        read particle and max energy
*-----------------------------------------------------------------------

         if( nrsq .eq. 0 .and. mrsq .eq. 0 ) then

               ic = i1

*-----------------------------------------------------------------------

            if( chlw(ic:ic+3) .eq. 'part' ) then

                  icnp = icnp + 1

                  if( icnp .gt. 1 ) goto 980

                  ic  = inumc(chlw,ic+3,i3,'=') + 1
                  ic  = jnumc(chlw,ic,i3)
                  icl = inumc(chlw,ic,i3,';') - 1

  400          continue

                  ic = jnumc(chlw,ic,icl)

                  if( ic .gt. icl ) then

                     icl = jnumc(chlw,icl+2,i3)

                     if( inpat .eq. 0 ) goto 994

                     goto 140

                  end if

*-----------------------------------------------------------------------

               call rdpname(ic,icl,chlw,istyp,inkf0,jstyp,jnkf0,ierr)

                  if( ierr .eq. 994 ) goto 994
                  if( ierr .eq. 998 ) goto 994

*-----------------------------------------------------------------------

                  if( istyp .gt. 0 ) then

                        inpat = inpat + 1
                        if( inpat .gt. 19 ) goto 992

                        iptyp(inpat) = istyp
                        ipnkf(inpat) = inkf0

                  else if( istyp .lt. 0 ) then

                     do i = 1, -istyp

                        inpat = inpat + 1
                        if( inpat .gt. 19 ) goto 992

                        iptyp(inpat) = jstyp(i)
                        ipnkf(inpat) = jnkf0(i)

                     end do

                  end if

                  goto 400

            end if

*-----------------------------------------------------------------------

            if( chlw(ic:ic+3) .eq. 'emax' ) then

               ic = inumc(chlw,ic+3,i3,'=') + 1
               ic = jnumc(chlw,ic,i3)

               if( ic .gt. i3 ) goto 992

               call snum(chlw,ic,i3,ic2,cvvv,ierr)

               if( ierr .ne. 0 ) goto 997

               ienp  = ienp + 1

               if( ienp .gt. 1 ) goto 981

               dlmax = cvvv

                  goto 140

            end if

*-----------------------------------------------------------------------
cfrtati 2023/12/07
            if( chlw(ic:ic+5).eq.'mtinfo' ) then
              ic = inumc(chlw,ic+3,i3,'=') + 1
              ic = jnumc(chlw,ic,i3)
              if( ic .gt. i3 ) goto 992
              call snum(chlw,ic,i3,ic2,cvvv,ierr)
              if( ierr .ne. 0 ) goto 997
              imtinf(1) = nint(cvvv)
              goto 140
            end if
            if( chlw(ic:ic+5).eq.'e-type' ) then
              ic = inumc(chlw,ic+3,i3,'=') + 1
              ic = jnumc(chlw,ic,i3)
              if( ic .gt. i3 ) goto 992
              call onum(chlw,ic,i3,cvvv,ierr)
              if( ierr .ne. 0 ) goto 997
              iemtp = nint(cvvv)
              call getmsh(jsn,jsi,dsin,idsi,ill,ilf,
     &                    jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr,
     &                    'e',iemtp,inem,emin,emax,edel,istegm)
              if( ierr .ne. 0 ) goto 997
            end if

*-----------------------------------------------------------------------

         end if

*-----------------------------------------------------------------------
*        read data definition
*-----------------------------------------------------------------------

         if( nrsq .eq. 0 ) then

               ic = i1

  100       if( ic .gt. i3 ) goto 200

            if(      chlw(ic:ic+2) .eq. 'non' ) then

               irsq( mrsq + 1 ) = 0
               ic = ic + 3

            else if( chlw(ic:ic+2) .eq. 'mat' ) then

               irsq( mrsq + 1 ) = 1
               ic = ic + 3

            else if( chlw(ic:ic+3) .eq. 'mset' ) then

               irsq( mrsq + 1 ) = 2
               ic = ic + 4

               call snum(chlw,ic,i3,ic2,cvvv,ierr)
               ic = ic2

               if( ierr .ne. 0 ) goto 982

               mmst = mmst + 1
               if( mmst .gt. 6 ) goto 983

               imst( mmst ) = nint( cvvv )
               if( imst(mmst) .le. 0 .or. imst(mmst) .gt. 6 ) goto 983

            else

               goto 998

            end if

               mrsq = mrsq + 1

               ic = jnumc(chlw,ic,i3)
               goto 100

  200       continue

            if( mrsq .gt. 0 ) then

                  imat = 0
                  jmst = 0
                  inon = 0

               do k = 1, mrsq
                  if( irsq(k) .eq. 1 ) imat = imat + 1
                  if( irsq(k) .eq. 2 ) jmst = jmst + 1
                  if( irsq(k) .eq. 0 ) inon = inon + 1
               end do

                  if( imat .ne. 1 .or. jmst .eq. 0 ) goto 998

                  nrsq = mrsq

                  goto 140

            else

                  nrsq = 2
                  mrsq = 2
                  mmst = 1
                  imst(1) = 1

            end if

         end if

*-----------------------------------------------------------------------

               mtcn = mtcn + 1
               ic2  = i1
               nmst = 0

         do k = 1, nrsq

               ic = jnumc(chlw,ic2,i3)

            if( irsq(k) .eq. 0 ) then

               if( chlw(ic:ic) .eq. '[' ) then

                  ic2 = inumc(chlw,ic,i3,']') + 1

               else if( chlw(ic:ic) .eq. '{' ) then

                  ic2 = inumc(chlw,ic,i3,'}') + 1

               else if( chlw(ic:ic) .eq. '(' ) then

                  ic2 = knump(chlw,ic,i3) + 1

               else

                  ic2 = inumc(chlw,ic,i3,' ')

               end if

            else if( irsq(k) .eq. 1 ) then

               if( chlw(ic:ic+2) .eq. 'all' ) then

                  mltp(1,jtli/13+1) = -1

                  ic2 = inumc(chlw,ic,i3,' ')

               else

                  if( jtal .ne. 0 ) goto 985

                  call snum(chlw,ic,i3,ic2,cvvv,ierr)

                  if( ierr .ne. 0 ) goto 997

                  mltp(1,jtli/13+mtcn) = nint( cvvv )

               end if

*-----------------------------------------------------------------------
*           read mset definitions
*                 ( :  1
*                 ) :  2
*                 : :  3
*                 C :  4
*                 m :  5
*                 R :  6
*              -1 m :  7
*           -1 m px :  8
*-----------------------------------------------------------------------

            else if( irsq(k) .eq. 2 ) then

                     ipar = 0
                     nmst = nmst + 1
                     if( nmst .gt. mmst ) goto 988

                     nnst = nmst * 2 - 1
                     mltp(1+nnst,jtli/13+mtcn) = jtls + itls

  500          continue

                  ic = jnumc(chlw,ic,i3)

               if( ic .gt. i3 .or. ipar .eq. 2 ) then

                  if( mmst .eq. 1 .and. ione .eq. 1 ) then

                     slib(jtls+itls) = 2.d0
                     itls = itls + 1

                     ibra = ibra - 1

                  end if

                     goto 600

               end if

               if( ipar .eq. 0 ) then

                     ibra = 0
                     icnt = 0
                     imtn = 0
                     ione = 0
                     ipar = 1

                  if( mmst .eq. 1 .and. chlw(ic:ic) .ne. '(' ) then

                     ibra = ibra + 1
                     ione = 1

                     slib(jtls+itls) = 1.d0
                     itls = itls + 1
                     slib(jtls+itls) = 0.d0
                     itls = itls + 1

                  end if

               end if

               if( chlw(ic:ic) .eq. '(' ) then

                     slib(jtls+itls) = 1.d0
                     itls = itls + 1
                     slib(jtls+itls) = 0.d0
                     itls = itls + 1

                     ibra = ibra + 1
                     ic = ic + 1

                     icnt = 0
                     imtn = 0

               else if( chlw(ic:ic) .eq. ')' ) then

                     slib(jtls+itls) = 2.d0
                     itls = itls + 1
                     slib(jtls+itls) = 0.d0
                     itls = itls + 1

                     ibra = ibra - 1
                     ic = ic + 1

                     if( ibra .eq. 0 ) ipar = 2

               else if( chlw(ic:ic) .eq. ':' ) then

                     if( icnt .eq. 0 .or. imtn .eq. 0 ) goto 987

                     slib(jtls+itls) = 3.d0
                     itls = itls + 1
                     slib(jtls+itls) = 0.d0
                     itls = itls + 1

                     ic = ic + 1

               else

                  do jj = ic, i3

                     if( dcom3( chlw(jj:jj) ) ) goto 510 ! T.Sato 2021/02/12 to use parameter and equation except for brancket

                  end do

  510                id = jj - 1


                     if( id .lt. ic ) goto 997

                     call onum(chlw,ic,id,cvvv,ierr)

            if( ierr .ne. 0 ) then
             ErrID = 'L:11516/R:tmultipl/F:tallsm1.f'
             write(ErrCha,'(''Wrong multiplier set definition.'',
     &       '' You cannot use branket in the equation here'')')
             call ErrWrite(ErrID,ErrCha)
             goto 997
            endif

                     ic = id + 1

                  if( icnt .eq. 0 ) then

                     cnst = cvvv
                     icnt = 1

                     slib(jtls+itls) = 4.d0
                     itls = itls + 1
                     slib(jtls+itls) = cvvv
                     itls = itls + 1

                  else if( imtn .eq. 0 ) then

                     mtnm = nint( cvvv )
                     imtn = 1

                     if( mtnm .eq. -1 ) imtn =  2
                     if( mtnm .le. -2 ) imtn = -2

                     slib(jtls+itls) = 5.d0
                     itls = itls + 1
                     slib(jtls+itls) = cvvv
                     itls = itls + 1

                  else if( imtn .eq. 1 ) then

                     slib(jtls+itls) = 6.d0
                     itls = itls + 1
                     slib(jtls+itls) = cvvv
                     itls = itls + 1

                  else if( imtn .ge. 2 .and. mod(imtn,2) .eq. 0 ) then

                     slib(jtls+itls) = 7.d0
                     itls = itls + 1
                     slib(jtls+itls) = cvvv
                     itls = itls + 1

                     imtn = imtn + 1

                  else if( imtn .ge. 2 .and. mod(imtn,2) .eq. 1 ) then

                     if( cvvv .le. 0.0 ) goto 978

                     slib(jtls+itls) = 8.d0
                     itls = itls + 1
                     slib(jtls+itls) = cvvv
                     itls = itls + 1

                     imtn = imtn + 1

                  else

                     goto 997

                  end if

               end if

                  goto 500

  600          continue

               if( ibra .ne. 0 ) goto 987

               if( ( slib(jtls+itls-3) .eq. 5 .and.
     &               slib(jtls+itls-2) .gt. 0.0 ) .or.
     &               slib(jtls+itls-3) .eq. 7 ) goto 979

                  mltp(1+nnst+1,jtli/13+mtcn) = jtls + itls - 1

                  ic2 = ic

            end if

*-----------------------------------------------------------------------

         end do

            goto 140

*-----------------------------------------------------------------------
*     summary
*-----------------------------------------------------------------------

 1000 continue

            if( inpat .eq. 0 ) then

                  inpat = 1
                  iptyp(1) = 20
                  ipnkf(1) = 0

            end if

            if( mtcn .ne. ntcn ) goto 995

            if( itls > MAX_NUM_SLIB ) then
               write(ErrCha,'(a,a,i5,a,a,i5,a)')
     &              'sub.tmultipl@tallsm1.f'
     &                 //' ?dimension over slib?'
     &                 //' itls > MAX_NUM_SLIB'
     &              ,' (itls=',itls,')'
     &              ,' (MAX_NUM_SLIB@moddas.f=',MAX_NUM_SLIB,')'
               ErrID = 'L:11628/R:tmultipl/F:tallsm1.f'
               call ErrWrite(ErrID,ErrCha)
            endif

            iaddress_mltp(2) = iaddress_mltp(3)
            icurrent_mltp = 2

            iaddress_slib(2) = iaddress_slib(3)
            icurrent_slib = 2

            if( mmmax .gt. mdas ) goto 970

cfrtati 2023/12/07 energy mesh for mtinfo
            if( imtinf(1).eq.1 ) then
              if( iemtp.eq.0 ) goto 968
              imtinf(2) = iemtp
              imtinf(3) = inem
              imtinf(4) = istegm
            end if

         return

*-----------------------------------------------------------------------
*     errors
*-----------------------------------------------------------------------

  968 continue

         m_err = 'e-type should be defined in multiplier subsection'//
     &           ' when mtinfo = 1'
         ErrCha = ''
         ErrID = 'L:11659/R:tmultipl/F:tallsm1.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

  970 continue

         m_err = 'Memory error: mmmax exceeds mdas '//
     &           ': Please extend mdas in param.inc'
         ErrCha = ''
         ErrID = 'L:11670/R:tmultipl/F:tallsm1.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  978 continue

         m_err = 'attenuation factor should be positive.'
         ErrCha = ''
         ErrID = 'L:11682/R:tmultipl/F:tallsm1.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  979 continue

         m_err = 'mset definition is wrong'
         ErrCha = ''
         ErrID = 'L:11694/R:tmultipl/F:tallsm1.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  980 continue

         m_err = 'part =  definition is appeared twice'
         ErrCha = ''
         ErrID = 'L:11706/R:tmultipl/F:tallsm1.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  981 continue

         m_err = 'emax =  definition is appeared twice'
         ErrCha = ''
         ErrID = 'L:11718/R:tmultipl/F:tallsm1.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  982 continue

         m_err = 'number after mset is wrong'
         ErrCha = ''
         ErrID = 'L:11730/R:tmultipl/F:tallsm1.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  983 continue

         m_err = 'number of mset should be 1 - 6'
         ErrCha = ''
         ErrID = 'L:11742/R:tmultipl/F:tallsm1.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  984 continue

         m_err = 'all is specified in mat but multiplier is not all.'
         ErrCha = ''
         ErrID = 'L:11754/R:tmultipl/F:tallsm1.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  985 continue

         m_err = 'all is specified in multiplier but mat is not all.'
         ErrCha = ''
         ErrID = 'L:11766/R:tmultipl/F:tallsm1.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  986 continue

         m_err = '{n1-n2} description is wrong.'
         ErrCha = ''
         ErrID = 'L:11778/R:tmultipl/F:tallsm1.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  987 continue

         m_err = 'mset discription is wrong.'
         ErrCha = ''
         ErrID = 'L:11790/R:tmultipl/F:tallsm1.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  988 continue

         m_err = 'number of mset definition exceeds mset index.'
         ErrCha = ''
         ErrID = 'L:11802/R:tmultipl/F:tallsm1.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  992 continue

         m_err = 'emax =  description is wrong'
         ErrCha = ''
         ErrID = 'L:11814/R:tmultipl/F:tallsm1.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  993 continue

         m_err = 'part = should be before emax =.'
         ErrCha = ''
         ErrID = 'L:11826/R:tmultipl/F:tallsm1.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  994 continue

         m_err = 'part = something wrong'
         ErrCha = ''
         ErrID = 'L:11838/R:tmultipl/F:tallsm1.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  995 continue

         m_err = 'Number of multiplier is inconsistent '//
     &           'with the # of data'
         ErrCha = ''
         ErrID = 'L:11851/R:tmultipl/F:tallsm1.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  996 continue

         m_err = 'number of part should be less than 6'
         ErrCha = ''
         ErrID = 'L:11863/R:tmultipl/F:tallsm1.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  997 continue

         m_err = 'Description of multiplier is wrong.'
         ErrCha = ''
         ErrID = 'L:11875/R:tmultipl/F:tallsm1.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  998 continue

         m_err = 'Index of multiplier is wrong.'
         ErrCha = ''
         ErrID = 'L:11887/R:tmultipl/F:tallsm1.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

      end


************************************************************************
*                                                                      *
      subroutine getmsh(jsn,jsi,dsin,idsi,ill,ilf,
     &                  jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr,
     &                  gdum,igtp,igr,gmin,gmax,gdel,igm)
*                                                                      *
*       get mesh information                                           *
*       modified by K.Niita on 2000/07/08                              *
*                                                                      *
************************************************************************
      use moddas
      use moddas_mesh
      use CHARVARMOD, only: ErrLine_Adjust

      implicit real*8 (a-h,o-z)

      include 'param.inc'
      include 'risrcparam.inc' ! T.Sato 2022/03/25
      include 'err.inc'

*-----------------------------------------------------------------------

      character m_err*200
      common /error/ m_err, l_err, k_err

*-----------------------------------------------------------------------

      common/nzerocom/nzerocheck ! T.Sato 2018/3/6, allow nz = 0 for [t-cross]

      character gdum*1

      dimension lschn(10), ischn(10)
      character schan(10)*8


      logical iprojall   ! multi-source subsection, proj=all
*-----------------------------------------------------------------------

      data icsu / 4 /

      data ( schan(i), i = 1, 4 ) /
     &    'n       ',' min    ',' max    ',' del    '/

      data ( lschn(i), i = 1, 4 ) /
     &     2,        4,          4,         4/

*-----------------------------------------------------------------------

      character chin*200, chlw*200, chcm*200
      character chlc*200
      character cwarn*200

      character dsin(0:9)*200
      dimension idsi(0:9)

      dimension ill(0:9), ilf(0:9)

*-----------------------------------------------------------------------
      ierrMSG=0 ! T.Sato 2017/07/14

         do i = 1, icsu

            ischn(i) = 0

         end do

*-----------------------------------------------------------------------
*     initial memory point
*-----------------------------------------------------------------------

         igm = iaddress_gmsh(icurrent_gmsh)

*-----------------------------------------------------------------------
*     check mesh type
*-----------------------------------------------------------------------

         if( igtp .lt. 1 .or. igtp .gt. 5 ) goto 999

*-----------------------------------------------------------------------
*     mesh type
*-----------------------------------------------------------------------

            schan(1)(2:2) = gdum

         do i = 2, icsu

            schan(i)(1:1) = gdum

         end do

*-----------------------------------------------------------------------

            ierr  = 0

            icmp = 0
            irct = 0

*-----------------------------------------------------------------------
*     read one line from jsi
*-----------------------------------------------------------------------

  140 continue

            call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

               if( ierr .ne. 0 ) return
               if( jpn  .eq. 3 ) goto 800

               if( iskip .ne. 0 ) goto 140

*-----------------------------------------------------------------------
*        end of the section
*-----------------------------------------------------------------------

            if( chlw(i1:i1) .eq. '[' ) then

               jpn = 1
               goto 800

            end if

*-----------------------------------------------------------------------
*        identify the parameters
*-----------------------------------------------------------------------

            icl = i1

  200    continue

            chlc = chlw
            call chcomp(chlc,icl,i3,i5)

         do i = 1, icsu

            il = icl + lschn(i) - 1

            if( chlc(icl:il) .eq. schan(i)(1:lschn(i)) ) goto 100

         end do

               goto 800

*-----------------------------------------------------------------------
*        read value of parameters
*-----------------------------------------------------------------------

  100    continue
ccse multi-source subsection, proj=all, write scrach file
            icl=i1
            chlc = chlw
            call chcomp(chlc,icl,i3,i5)
            if( chlc(icl:icl+7) .ne. '<source>' .and.
     &          chlc(icl:icl) .ne. '[' ) then
               inquire(niws,opened=iprojall)
               if(iprojall) write(niws,'(a)') chin(1:i2)
            endif

               ipm = i

               ischn( ipm ) = 1

               ic = inumc(chlw,il+1,i3,'=') + 1
               ic = jnumc(chlw,ic,i3)

               if( ic .gt. i3 ) goto 998

               icl = inumc(chlw,ic,i3,';') - 1

*-----------------------------------------------------------------------
*        mesh informations
*-----------------------------------------------------------------------

               call onum(chlw,ic,icl,cvvv,ierr)

               if( ierr .ne. 0 ) goto 998

         if( ipm .eq. 1 ) then

               igrr = nint( cvvv )
               igr  = iabs( igrr )

               if( igr.le.0.and.nzerocheck.eq.0 ) goto 997 ! T.Sato 2018/3/6

               call moddas_reallocate_dbl(
     &                 3, icurrent_gmsh, igr+1, iaddress_gmsh, gmsh)

               if( igtp .eq. 1 ) goto 500

         else if( ipm .eq. 2 ) then

               gmin = cvvv

         else if( ipm .eq. 3 ) then

               gmax = cvvv

         else if( ipm .eq. 4 ) then

               gdel = cvvv

         end if

               icl = jnumc(chlw,icl+2,i3)

               if( icl .le. i3 ) goto 200

               goto 140

*-----------------------------------------------------------------------
*        read mesh points
*-----------------------------------------------------------------------

  500    continue

                  if( chlw(icl+1:icl+1) .eq. ';' ) goto 996

  148             call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

                  if( ierr .ne. 0 ) return
                  if( jpn  .eq. 3 ) goto 997

                  if( iskip .ne. 0 ) goto 148

ccse multi-source subsection, proj=all, write scrach file
            icl=i1
            chlc = chlw
            call chcomp(chlc,icl,i3,i5)
            if( chlc(icl:icl+7) .ne. '<source>' .and.
     &          chlc(icl:icl) .ne. '[' ) then
               inquire(niws,opened=iprojall)
               if(iprojall) write(niws,'(a)') chin(1:i2)
            endif
                  ic = i1

            do i = 1, igr + 1

               call snum(chlw,ic,i3,ic2,cvvv,ierr)

               if( ierr .ne. 0 ) goto 994 ! T.Sato 2024/09/23

                  gmsh(igm-1+i) = cvvv

               ic = ic2

               if( i .le. igr .and. ic .gt. i3 ) then

  149             call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

                  if( ierr .ne. 0 ) return
                  if( jpn  .eq. 3 ) goto 997

                  if( iskip .ne. 0 ) goto 149

ccse multi-source subsection, proj=all, write scrach file
            icl=i1
            chlc = chlw
            call chcomp(chlc,icl,i3,i5)
            if( chlc(icl:icl+7) .ne. '<source>' .and.
     &          chlc(icl:icl) .ne. '[' ) then
               inquire(niws,opened=iprojall)
               if(iprojall) write(niws,'(a)') chin(1:i2)
            endif

                  ic = i1

               end if

            end do

! T.Sato 2024/09/23, add warning
            call snum(chlw,ic,i3,ic2,cvvv,ierr) ! try to read additional data
            if(ierr.eq.0) then ! normal data
             if(cvvv.gt.gmsh(igm-1+igr+1)) then ! greater than the last data, so may be number of mesh is wrong
              l_warn=ill(jsn)
              call ErrLine_Adjust(cwarn,l_warn)
         write(ErrCha,'("Warning in ",a," at line ",i8,
     &   ": number of data is more than that of meshs")')
     &   trim(cwarn),l_warn
         ErrID = 'L:12179/R:getmsh/F:tallsm1.f'
         call ErrWrite(ErrID,ErrCha)

!         write(*,*) 'Warning in ',trim(cwarn),' at line ',
!     &        l_warn,': number of data is more than that of meshs'
             endif
            else
             ierr=0
            endif

                  icmp = 1

               goto 140

*-----------------------------------------------------------------------
*     summary and check
*-----------------------------------------------------------------------

  800 continue

            if( igtp .eq. 1 ) then

               if( ischn(1) .eq. 0 .or.
     &             (igr.le.0.and. nzerocheck.eq.0)  .or.  ! T.Sato 2018/3/6
     &             icmp .eq. 0 ) then
                ierrMSG=1 ! T.Sato 2017/07/14
                goto 998
               endif

               do i = 1, igr

                     gmins = gmsh(igm-1+i)

                  do j = i + 1, igr + 1

                     if( gmsh(igm-1+j) .lt. gmins ) then

                        gmins         = gmsh(igm-1+j)
                        gmsh(igm-1+j) = gmsh(igm-1+i)
                        gmsh(igm-1+i) = gmins

                     end if

                  end do

               end do

            else if( igtp .eq. 2 ) then

               if( ischn(1) .eq. 0 .or.
     &             ischn(2) .eq. 0 .or.
     &             ischn(3) .eq. 0 .or.
     &             gmax .le. gmin  .or.
     &             igr  .le. 0 ) then
                  ierrMSG=2 ! T.Sato 2017/07/14
                  goto 998
               endif

                  gdel = ( gmax - gmin ) / dble( igr )

               do i = 1, igr

                  gmsh(igm-1+i) = gmin + gdel * dble( i - 1 )

               end do

                  gmsh(igm-1+igr+1) = gmax

            else if( igtp .eq. 3 ) then

               if( ischn(1) .eq. 0 .or.
     &             ischn(2) .eq. 0 .or.
     &             ischn(3) .eq. 0 .or.
     &             gmin .le. 0.0   .or.
     &             gmax .le. gmin  .or.
     &             igr  .le. 0 ) then
                ierrMSG=3 ! T.Sato 2017/07/14
                goto 998
               endif

                  gdel = log( gmax / gmin ) / dble( igr )

               do i = 1, igr

                  gmsh(igm-1+i) = gmin * exp( gdel * dble( i - 1 ) )

               end do

                  gmsh(igm-1+igr+1) = gmax

            else if( igtp .eq. 4 ) then

               if( ischn(4) .eq. 0 .or.
     &             ischn(2) .eq. 0 .or.
     &             ischn(3) .eq. 0 .or.
     &             gmax .le. gmin  .or.
     &             gdel .le. 0.0 ) then
                ierrMSG=4 ! T.Sato 2017/07/14
                goto 998
               endif

                  igr = int( ( gmax - gmin ) / gdel ) + 1

               call moddas_reallocate_dbl(
     &                 3, icurrent_gmsh, igr+1, iaddress_gmsh, gmsh)

               do i = 1, igr + 1

                  gmsh(igm-1+i) = gmin + gdel * dble( i - 1 )

               end do

                  if( gmsh(igm-1+igr) .eq. gmax ) igr = igr - 1

                  igrr = igr

            else if( igtp .eq. 5 ) then

               if( ischn(4) .eq. 0 .or.
     &             ischn(2) .eq. 0 .or.
     &             ischn(3) .eq. 0 .or.
     &             gmin .le. 0.0   .or.
     &             gmax .le. gmin  .or.
     &             gdel .le. 0.0 ) then
                ierrMSG=5 ! T.Sato 2017/07/14
                goto 998
               endif

                  igr = int( log( gmax / gmin ) / gdel ) + 1

               call moddas_reallocate_dbl(
     &                 3, icurrent_gmsh, igr+1, iaddress_gmsh, gmsh)

               do i = 1, igr + 1

                  gmsh(igm-1+i) = gmin * exp( gdel * dble( i - 1 ) )

               end do

                  if( gmsh(igm-1+igr) .eq. gmax ) igr = igr - 1

                  igrr = igr

            end if

                  iaddress_gmsh(2) = iaddress_gmsh(3)
                  icurrent_gmsh = 2

                  if( mmmax .ge. mdas ) goto 995

               igr = igrr

      return

*-----------------------------------------------------------------------
*     errors
*-----------------------------------------------------------------------

  994 continue

         m_err ='Number of data is insufficient or mesh value is'//
     &   ' wrong in this or above line(s)'
         ErrCha = ''
         ErrID = 'L:12342/R:getmsh/F:tallsm1.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

  995 continue

         m_err = 'Memory error: mmmax exceeds mdas '//
     &           ': Please extend mdas in param.inc'
         ErrCha = ''
         ErrID = 'L:12353/R:getmsh/F:tallsm1.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  996 continue

         m_err = 'In this line, [ ; ] cannot be used.'
         ErrCha = ''
         ErrID = 'L:12365/R:getmsh/F:tallsm1.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  997 continue

         m_err = 'Description of mesh points is wrong'
         ErrCha = ''
         ErrID = 'L:12377/R:getmsh/F:tallsm1.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  998 continue
      if(ierrMSG.le.1.and.ischn(1).eq.0) then ! T.Sato 2017/07/14
       m_err='Number of mesh points should be defined for *-type =< 3'
       ErrCha = ''
       ErrID = 'L:12389/R:getmsh/F:tallsm1.f'
      elseif(ierrMSG.eq.1) then
       m_err='Each mesh point must be defined for *-type = 1'
       ErrCha = ''
       ErrID = 'L:12393/R:getmsh/F:tallsm1.f'
      elseif(ierrMSG.ge.2.and.(ischn(2).eq.0.or.ischn(3).eq.0)) then
       m_err='Both *min= & *max= should be defined for *-type >= 2'
       ErrCha = ''
       ErrID = 'L:12397/R:getmsh/F:tallsm1.f'
      elseif(ierrMSG.ge.4.and.ischn(4).eq.0) then
       m_err='*del= should be defined for *-type = 4 or 5'
       ErrCha = ''
       ErrID = 'L:12401/R:getmsh/F:tallsm1.f'
      elseif((ierrMSG.eq.3.or.ierrMSG.eq.5).and.gmin.le.0.0) then
       m_err='Minimum value should be positive for *-type = 3 or 5'
       ErrCha = ''
       ErrID = 'L:12405/R:getmsh/F:tallsm1.f'
      elseif(ierrMSG.ge.2.and.ierrMSG.le.5.and.gmax.le.gmin) then
       m_err='max value is less than min value in mesh definition'
       ErrCha = ''
       ErrID = 'L:12409/R:getmsh/F:tallsm1.f'
      else
       m_err = 'Description of mesh parameter is wrong'
       ErrCha = ''
       ErrID = 'L:12413/R:getmsh/F:tallsm1.f'
      endif
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

  999 continue

         m_err = 'Mesh type should be 1 - 5'
         ErrCha = ''
         ErrID = 'L:12426/R:getmsh/F:tallsm1.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

      end


************************************************************************
*                                                                      *
      subroutine ttrans(jsn,jsi,dsin,idsi,ill,ilf,
     &                  jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr,
     &                  ic,itrs,igkst,idtt,vtrs)
*                                                                      *
*       read trcl = subsection of tallies                              *
*       modified by K.Niita on 2010/01/14                              *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      include 'param.inc'
      include 'risrcparam.inc' ! T.Sato 2022/03/25
      include 'err.inc'

*-----------------------------------------------------------------------

      character m_err*200
      common /error/ m_err, l_err, k_err

*-----------------------------------------------------------------------

      common /inggs/ iog, igcel, ioa, igsuf, iob, igtrs

      common /celdb/ idsn(kvlmax), idtn(kvlmax)

*-----------------------------------------------------------------------

      character chin*200, chlw*200, chcm*200
      character chlc*200

      character dsin(0:9)*200
      dimension idsi(0:9)

      dimension ill(0:9), ilf(0:9)

      logical iprojall   ! multi-source subsection, proj=all
*-----------------------------------------------------------------------

      dimension vtrs(13)

*-----------------------------------------------------------------------

            ierr  = 0

*-----------------------------------------------------------------------

               vtrs( 1) = 0.0
               vtrs( 2) = 0.0
               vtrs( 3) = 0.0
               vtrs( 4) = 1.0
               vtrs( 5) = 0.0
               vtrs( 6) = 0.0
               vtrs( 7) = 0.0
               vtrs( 8) = 1.0
               vtrs( 9) = 0.0
               vtrs(10) = 0.0
               vtrs(11) = 0.0
               vtrs(12) = 1.0
               vtrs(13) = 1.0

            if( itrs .eq. 1 ) then

               vtrs( 4) =  0.0
               vtrs( 5) = 90.0
               vtrs( 6) = 90.0
               vtrs( 7) = 90.0
               vtrs( 8) =  0.0
               vtrs( 9) = 90.0
               vtrs(10) = 90.0
               vtrs(11) = 90.0
               vtrs(12) =  0.0

            end if

            goto 100

*-----------------------------------------------------------------------
*     read one more line from jsi
*-----------------------------------------------------------------------

  140 continue

            call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

               if( ierr .ne. 0 ) goto 1000
               if( jpn  .eq. 3 ) goto 1000

               if( iskip .ne. 0 ) goto 140

ccse multi-source subsection, proj=all, write scrach file
            icl=i1
            chlc = chlw
            call chcomp(chlc,icl,i3,i5)
            if( chlc(icl:icl+7) .ne. '<source>' .and.
     &          chlc(icl:icl) .ne. '[' ) then
               inquire(niws,opened=iprojall)
               if(iprojall) write(niws,'(a)') chin(1:i2)
            endif
*-----------------------------------------------------------------------
*        end of this subsection
*-----------------------------------------------------------------------

            if( i1 .le. 5 .and. chlw(i1:i1) .eq. '[' ) then

               jpn = 1
               goto 1000

            end if

            if( i1 .le. 11 ) goto 1000

               ic = i1

*-----------------------------------------------------------------------
*        rest part or sequential line
*-----------------------------------------------------------------------

  100    continue

            if( ic .gt. i3 ) goto 140

                  call snum(chlw,ic,i3,ic2,cvvv,ierr)

                  if( ierr .ne. 0 ) goto 999

                  igkst = igkst + 1

                  if( igkst .gt. 13 ) goto 999

                  vtrs(igkst) = cvvv

                  ic = jnumc(chlw,ic2,i3)

                  goto 100

*-----------------------------------------------------------------------
*        write information on temporary file unit = 71
*-----------------------------------------------------------------------

 1000 continue

            if( igkst .le. 0 ) goto 999

            if( igkst .eq. 1 ) then

               idtt = nint( vtrs(1) )

*-----------------------------------------------------------------------
*              open temporary file : unit 71
*-----------------------------------------------------------------------

            else if( igkst .gt. 1 ) then

               if( igtrs .eq. 0 ) then

                  iob = 71

                  open(iob,status='scratch',form='unformatted')

               end if

*-----------------------------------------------------------------------

               igtrs = igtrs + 1

               idtn(igtrs) = 1000000 + igtrs

               write(iob) itrs, ( vtrs(i), i = 1, 13 )

               idtt = 1000000 + igtrs

            end if

            return

*-----------------------------------------------------------------------
*     errors
*-----------------------------------------------------------------------

  999 continue

         m_err = 'Description of trcl is wrong.'
         ErrCha = ''
         ErrID = 'L:12624/R:ttrans/F:tallsm1.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

      end

************************************************************************
*                                                                      *
      subroutine stbatm(jsn,jsi,dsin,idsi,ill,ilf,
     &                  jpn,chin,chlw,chcm,i1,i2,i3,i4,ierr,
     &                  ist_nm,isd_nm,kst_nm)
*                                                                      *
*       read anatally subsection of tallies                            *
*       modified by K.Niita on 2019/07/31                              *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      include 'param.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      character m_err*200
      common /error/ m_err, l_err, k_err

*-----------------------------------------------------------------------

      parameter ( npanatal = 12 ) ! S.H. 2020.3.27

      dimension ist_nm(npanatal), isd_nm(npanatal,kvlmax),
     &     kst_nm(npanatal)

*-----------------------------------------------------------------------

      dimension lschn(npanatal)
      character schan(npanatal)*9

*-----------------------------------------------------------------------

      data icsu / npanatal /

      data ( schan(i), i = 1, npanatal ) /
     &    'ireg     ','ix       ','iy       ','iz       ','ir       ',
     &    'ie       ','it       ','ia       ','ipart    ','imul     ',
     &    'manatally','sfile    '/

      data ( lschn(i), i = 1, npanatal ) /
     &     4,      2,      2,      2,      2,
     &     2,      2,      2,      5,      4,
     &     9,      5/


*-----------------------------------------------------------------------

      character chin*200, chlw*200, chcm*200
      character chlc*200

      character dsin(0:9)*200
      dimension idsi(0:9)

      dimension ill(0:9), ilf(0:9)

      dimension ntrg(kvlmax)

      logical deqn4
      logical dnen1

*-----------------------------------------------------------------------

            ierr  = 0
            jpn = 0

*-----------------------------------------------------------------------
*     read one more line from jsi
*-----------------------------------------------------------------------

  140 continue

            call readl(jsn,jsi,dsin,idsi,ill,ilf,'#!$',
     &                 jpn,chin,chlw,chcm,i1,i2,i3,i4,iskip,ierr)

               if( ierr .ne. 0 ) goto 1000
               if( jpn  .eq. 3 ) goto 1000

               if( iskip .ne. 0 ) goto 140

*-----------------------------------------------------------------------
*        end of this subsection
*-----------------------------------------------------------------------

            if( i1 .le. 5 .and. chlw(i1:i1) .eq. '[' ) then

               jpn = 1
               goto 1000

            end if

            if ( chlw(i1:i3-5) .eq. 'anatally end' ) goto 1000 ! "anatally end  nnn"

*-----------------------------------------------------------------------
*        identify the parameters
*-----------------------------------------------------------------------

            icl = i1

  200    continue

            chlc = chlw
            call chcomp(chlc,icl,i3,i5)

         do i = 1, icsu

            il = icl + lschn(i) - 1

            if( chlc(icl:il) .eq. schan(i)(1:lschn(i)) ) then

               ic = inumc(chlw,il+1,i3,'=') + 1
               ic = jnumc(chlw,ic,i3)

               goto 100

            end if

         end do

               goto 1000

*-----------------------------------------------------------------------
*        rest part or sequential line
*-----------------------------------------------------------------------

  100    continue

                  ntrn = 0

            if( ic .gt. i3 ) goto 140

                  ipm = i

               if( ist_nm(ipm) .eq. -1 ) goto 140

               if ( ipm.eq.12 ) goto 140 ! to avoid error when reading sfile (S.H. 2020.3.27)

               if( chlw(ic:ic+2) .eq. 'all' ) then

                  ist_nm(ipm) = -1
                  goto 140

               end if

*-----------------------------------------------------------------------

                  call snum(chlw,ic,i3,ic2,cvvv,ierr)

                  if( ierr .ne. 0 ) goto 999

                  ist_nm(ipm) = ist_nm(ipm) + 1

                  if( ist_nm(ipm) .gt. kvlmax ) goto 998

                  isd_nm(ipm,ist_nm(ipm)) = nint(cvvv)

                  ic = jnumc(chlw,ic2,i3)

                  goto 100

*-----------------------------------------------------------------------
*        summary
*-----------------------------------------------------------------------

 1000 continue

            do i = 1, npanatal
               kst_nm(i) = ist_nm(i)
            end do

         return

*-----------------------------------------------------------------------
*     errors
*-----------------------------------------------------------------------
  995 continue

         m_err = 'Description of region number '//
     &           '{n1-n2} (n1<n2) is wrong.'
         ErrCha = ''
         ErrID = 'L:12816/R:stbatm/F:tallsm1.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------
  996 continue

         m_err = 'Description of region number '//
     &           '( {n1-n2} n3 n4 ) is wrong.'
         ErrCha = ''
         ErrID = 'L:12828/R:stbatm/F:tallsm1.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------
  998 continue

         m_err = 'Max number of one variable is kvlmax.'
         ErrCha = ''
         ErrID = 'L:12839/R:stbatm/F:tallsm1.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------
  999 continue

         m_err = 'Description of anatally subsection is wrong.'
         ErrCha = ''
         ErrID = 'L:12850/R:stbatm/F:tallsm1.f'
         l_err = ill(jsn)
         k_err = jsn
         ierr  = 1
         return

*-----------------------------------------------------------------------

      end


************************************************************************
*                                                                      *
      subroutine anset(chlw,icl,
     &                 form,xfac,afac,cmin,cmax,izlog,inocm,inolg)
*                                                                      *
*        purpose : get angel para of form,xfac,afac,cmin,cmax,ylog     *
*                  for 2d-type = 3 legend                              *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      character chlw*200
      character chin*200

*-----------------------------------------------------------------------

         do i = 1, 200

            if( i .le. icl ) then

               chin(i:i) = chlw(i:i)

            else

               chin(i:i) = ' '

            end if

         end do

            call chcaps(chin,1,icl,i3,'#!$')

*-----------------------------------------------------------------------

         ic = 0

  100    ic = ic + 1

         if( ic .gt. icl ) return

         if( chin(ic:ic+3) .eq. 'form' ) then

               ic = ic + 4

               call snum(chin,ic,icl,ic2,cvvv,ierr)

               if( ierr .ne. 0 ) goto 100

               form = cvvv

               ic = ic2 - 1

         else if( chin(ic:ic+3) .eq. 'xfac' ) then

               ic = ic + 4

               call snum(chin,ic,icl,ic2,cvvv,ierr)

               if( ierr .ne. 0 ) goto 100

               xfac = cvvv

               ic = ic2 - 1

         else if( chin(ic:ic+3) .eq. 'afac' ) then

               ic = ic + 4

               call snum(chin,ic,icl,ic2,cvvv,ierr)

               if( ierr .ne. 0 ) goto 100

               afac = cvvv

               ic = ic2 - 1

         else if( chin(ic:ic+3) .eq. 'cmin' ) then

               ic = ic + 4

               call snum(chin,ic,icl,ic2,cvvv,ierr)

               if( ierr .ne. 0 ) goto 100

               cmin = cvvv

               ic = ic2 - 1

         else if( chin(ic:ic+3) .eq. 'cmax' ) then

               ic = ic + 4

               call snum(chin,ic,icl,ic2,cvvv,ierr)

               if( ierr .ne. 0 ) goto 100

               cmax = cvvv

               ic = ic2 - 1

         else if( chin(ic:ic+3) .eq. 'zlin' ) then

               izlog = 0

         else if( chin(ic:ic+3) .eq. 'nocm' ) then

               inocm = 0

         else if( chin(ic:ic+3) .eq. 'nolg' ) then

               inolg = 0

         end if

            goto 100

*-----------------------------------------------------------------------

      return
      end

************************************************************************
      subroutine avoidint(ixtp,inx,xmin,xmax,xdel) ! T.Sato 2024/03/23
! avoid integer value for mesh, originally introduced for [t-wwg]
! use x-mesh as an example, but can be applied to other parameter
************************************************************************
      implicit real*8 (a-h,o-z)
      common /paraj/ mstz(300), parz(300)

!      write(*,*) 'original',xmin,xmax,xdel
      if(ixtp.eq.2.or.ixtp.eq.4) then ! necessary only for type 2 and 4
       xmin=xmin-parz(204)       ! slightly decrease
       xmax=xmax+parz(204)*2.0d0 ! slihgtly increase
       xdel=(xmax-xmin)/inx
      endif
!      write(*,*) 'adjusted',xmin,xmax,xdel
      return
      end

