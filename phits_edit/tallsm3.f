************************************************************************
*                                                                      *
      subroutine twbgech(iot,m,iaxs,icax)
*                                                                      *
*       input echo of the t-wwbg tally                                 *
*       last modified by K.Niita on 2017/02/16                         *
*                                                                      *
************************************************************************

      use sumtallymod
      use partmod, only: itmxpt, itmxgp, itpan, itpat, jtpat ! frtati 2021/10/05
      use moddas_mesh
      use moddas_region

*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'

*-----------------------------------------------------------------------

      common /tall00/ itnm, ital(itlmax), itals(itlmax), italm(itlmax)
      common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                itout(itlmax), ittwo(itlmax)

      common /tall02/ itrgn(itlmax), itrgm(itlmax), itreg(itlmax)
      common /tall04/ itxty(itlmax), itxnm(itlmax), itxrg(itlmax),
     &                rtxmi(itlmax), rtxma(itlmax), rtxdl(itlmax)
      common /tall05/ ityty(itlmax), itynm(itlmax), ityrg(itlmax),
     &                rtymi(itlmax), rtyma(itlmax), rtydl(itlmax)
      common /tall06/ itzty(itlmax), itznm(itlmax), itzrg(itlmax),
     &                rtzmi(itlmax), rtzma(itlmax), rtzdl(itlmax)
      common /tall07/ itety(itlmax), itenm(itlmax), iterg(itlmax),
     &                rtemi(itlmax), rtema(itlmax), rtedl(itlmax)
      common /tall10/ itaxn(itlmax), itaxs(itlmax,6)

      common /tall20/ ittty(itlmax), ittnm(itlmax), ittrg(itlmax),
     &                rttmi(itlmax), rttma(itlmax), rttdl(itlmax)
      common /tall24/ itrcm(itlmax), itrcg(itlmax)
      common /tall25/ itrnv(itlmax), itriv(itlmax), itrrv(itlmax)

      common /tall14/ ittll(itlmax), ittle(itlmax)
      character ittle*80
      common /tall31/ iterl(itlmax)

      common /tall11/ itfln(itlmax), itfll(itlmax,6), ctfln(itlmax,6)
      character ctfln*100

      common /tall19/ itsmn(itlmax), itstm(itlmax)
      common /tall46/ itmtr(itlmax,4), rtmtr(itlmax,13)

      common /tall40/ rtorg(itlmax,3), rteye(itlmax,3), rtlit(itlmax,3),
     &                rtwin(itlmax,3), rtbrt(itlmax,2), ithvn(itlmax),
     &                itwin(itlmax,2), itbox(itlmax), itmir(itlmax),
     &                rtbox(itlmax,5,10), rtout(itlmax), rthet(itlmax),
     &                itlin(itlmax), itshd(itlmax), itgxs(itlmax)

      common /tall68/ itstp(itlmax), itmth(itlmax), rtvr0(itlmax),
     &                rtvx0(itlmax), rtvy0(itlmax), rtvz0(itlmax),
     &                rtvx1(itlmax), rtvy1(itlmax), rtvz1(itlmax)

      common /tall67/ itstd(itlmax), rtstd(itlmax)

      common /tall69/ itnms(itlmax), trmsh(itlmax,10),
     &                tzmsh(itlmax,10), tfmsh(itlmax,10),
     &                trmpo(itlmax,10), tzmpo(itlmax,10)

      dimension idas(mdas*2) !20220906frtati temporary for fbounds-check
      equivalence ( das, idas )

*-----------------------------------------------------------------------

      common /inpec/ ititl, ipara, ibody, iregn, llarr, itby, itar

*-----------------------------------------------------------------------
       integer   icntl
       integer   inucr

       common /tcntl/  icntl, inucr
!-----------------------------------------------------------------------

      character asfil*100

      dimension ian(6)

      character aname(4)*4
      data aname / 'wwbg',' xy ',' yz ', ' xz '/

*-----------------------------------------------------------------------
      common /tall79/ itallech
      if( itallech .eq. 0 ) return

*-----------------------------------------------------------------------
*     input echo
*-----------------------------------------------------------------------

               if( icax .eq. 0 ) then
                  write(iot,'("[ T-WWBG ]")')
               else
                  write(iot,'("[ T-WWBG ] off")')
               end if

            if( itmsh(m) .eq. 1 ) then

                  write(iot,'("     mesh =  reg",11x,
     &            " # mesh type is region-wise")')

                     idas1 = mmmax
                     idas2 = ( idas1 + itrgn(m) - 1 ) * 2 + 1
                     idas3 = idas2 + itrgn(m)

             call echrg_wwvm(1,iterl(m),0,1,iot,itrcm(m),
     &                       idas_itrcg(itrcg(m)),
     &                       itrgn(m),itrgm(m),idas_itreg(itreg(m)),
     &                       das(idas1),idas(idas2),
     &                       itrnv(m),idas(itriv(m)),das(itrrv(m)),
     &                       idas3)

            else if( itmsh(m) .eq. 3 ) then

                  write(iot,'("     mesh =  xyz",11x,
     &            " # mesh type is xyz scoring mesh")')

!<-20220912murofushi update
!--                  call echmty(0,iot,'x',itxty(m),itxnm(m),
!--     &                        rtxdl(m),rtxmi(m),rtxma(m),itxrg(m))
!--
!--                  call echmty(0,iot,'y',ityty(m),itynm(m),
!--     &                        rtydl(m),rtymi(m),rtyma(m),ityrg(m))
!--
!--                  call echmty(0,iot,'z',itzty(m),itznm(m),
!--     &                        rtzdl(m),rtzmi(m),rtzma(m),itzrg(m))
                  call echmty(0,iot,'x',itxty(m),itxnm(m),
     &                        rtxdl(m),rtxmi(m),rtxma(m),itxrg(m)
     &                        ,abs(itxnm(m)+1),das_itxrg)
                  call echmty(0,iot,'y',ityty(m),itynm(m),
     &                        rtydl(m),rtymi(m),rtyma(m),ityrg(m)
     &                        ,abs(itynm(m)+1),das_ityrg)
                  call echmty(0,iot,'z',itzty(m),itznm(m),
     &                        rtzdl(m),rtzmi(m),rtzma(m),itzrg(m)
     &                        ,abs(itznm(m)+1),das_itzrg)
!--->


            else if( itmsh(m) .eq. 4 ) then

                  write(iot,'("     mesh =  tet",11x,
     &            " # mesh type is tetra geometry")')

                  call echtet(iot,itrgm(m),idas_itreg(itreg(m)))

            end if


*-----------------------------------------------------------------------

                  asfil = '  # file name of output for the above axis'

            if( icax .eq. 0 ) then

                  icn = itaxn(m)

               do i = 1, icn

                  ian(i) = i

               end  do

            else

                  icn = 1
                  ian(1) = iaxs

            end if

*-----------------------------------------------------------------------

            do i = 1, icn

                  iax = ian(i)

                  write(iot,'("     axis = ",a4,11x,
     &            " # axis of output")')
     &            aname( itaxs(m,iax) )

            if ( itfll(m,iax) .gt. 0 ) then

                  msfile = max( 14, itfll(m,iax) )

                  write(iot,'("     file = ",100a1)')
     &                  ( ctfln(m,iax)(j:j), j = 1, msfile ),
     &                  ( asfil(j:j), j = 1, 42 )

            end if

            end do

*-----------------------------------------------------------------------
               call echprt(1,iot,m,itpan,itpat,jtpat,itnm,itmxpt,itmxgp)

                  call echmty(0,iot,'e',itety(m),itenm(m),
     &                        rtedl(m),rtemi(m),rtema(m),iterg(m)
     &                        ,abs(itenm(m)+1),das_iterg)

*-----------------------------------------------------------------------

               if( ittty(m) .ne. 0 ) then

                  call echmty(0,iot,'t',ittty(m),ittnm(m),
     &                        rttdl(m),rttmi(m),rttma(m),ittrg(m)
     &                        ,abs(ittnm(m)+1),das_ittrg)

               end if

*-----------------------------------------------------------------------

               if( rtstd(m) .gt. 0.d0 ) then

                  write(iot,'("   stdcut =",1p1g14.7,1x,1x,
     &            " # (D=-1) threshold value of STD cutoff")')
     &            rtstd(m)

               end if

*-----------------------------------------------------------------------

                  write(iot,'("       x0 =",1p1g14.7,1x,1x,
     &            " #  x of initial point")')
     &            rtvx0(m)

                  write(iot,'("       y0 =",1p1g14.7,1x,1x,
     &            " #  y of initial point")')
     &            rtvy0(m)

                  write(iot,'("       z0 =",1p1g14.7,1x,1x,
     &            " #  z of initial point")')
     &            rtvz0(m)

                  write(iot,'("       x1 =",1p1g14.7,1x,1x,
     &            " #  x of final point")')
     &            rtvx1(m)

                  write(iot,'("       y1 =",1p1g14.7,1x,1x,
     &            " #  y of final point")')
     &            rtvy1(m)

                  write(iot,'("       z1 =",1p1g14.7,1x,1x,
     &            " #  z of final point")')
     &            rtvz1(m)

*-----------------------------------------------------------------------

                  write(iot,'("   n-mesh = ",i6,9x,
     &            " # number of region")') itnms(m)

                  write(iot,'("   r-mesh = ",
     &                        5(1p1e13.5),
     &                        (/12x,5(1p1e13.5)))')
     &                        (trmsh(m,j),j=1,itnms(m))

                  write(iot,'("   z-mesh = ",
     &                        5(1p1e13.5),
     &                        (/12x,5(1p1e13.5)))')
     &                        (tzmsh(m,j),j=1,itnms(m))

                  write(iot,'("   f-mesh = ",
     &                        5(1p1e13.5),
     &                        (/12x,5(1p1e13.5)))')
     &                        (tfmsh(m,j),j=1,itnms(m)+1)

*-----------------------------------------------------------------------

              if( rtout(m) .gt. 1.0 ) then

                  write(iot,'("    r-out =",1p1g14.7,1x,1x,
     &            " # (D=0:no set) radius of outer sphere")')
     &            rtout(m)

              end if

*-----------------------------------------------------------------------
*           no write transform
*-----------------------------------------------------------------------

            if( itmtr(m,1) .gt. 10000 ) then

               if( itmtr(m,1) .eq. 1 ) then

                     write(iot,'("     trcl = ",i6,9x,
     &               " # transform ID")') itmtr(m,3)

               else

                  if( itmtr(m,2) .eq. 0 ) then

                     write(iot,'("     trcl = ",
     &                           3(1p1e15.7),
     &                           3(/12x,3(1p1e15.7)),i5)')
     &               ( rtmtr(m,j), j = 1, 12 ),
     &               nint( rtmtr(m,13) )

                  else

                     write(iot,'("    *trcl = ",
     &                           3(1p1e15.7),
     &                           3(/12x,3(1p1e15.7)),i5)')
     &               ( rtmtr(m,j), j = 1, 12 ),
     &               nint( rtmtr(m,13) )

                  end if

               end if

            end if

*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*        sumtally subsection
*-----------------------------------------------------------------------
            if ( icntl==13 .and. nsumtalRead(m)==1 ) then

               call sumtal_echo(iot,m)

            end if

*-----------------------------------------------------------------------

            write(iot,'("#    used :",
     &                  "        main (  %)",
     &                  "        temp (  %)",
     &                  "       total (  %)"/
     &                  "#  memory :",
     &                  3(i12," (",i3,")"))')
     &                  itsmn(m), nint(dble(itsmn(m))/dble(mdas)*100),
     &                  itstm(m), nint(dble(itstm(m))/dble(mdas)*100),
     &                  itsmn(m)+itstm(m),
     &                  nint(dble(itsmn(m)+itstm(m))/dble(mdas)*100)

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine tvolech(iot,m,iaxs,icax)
*                                                                      *
*       input echo of the t-volume tally                               *
*       last modified by K.Niita on 2017/01/02                         *
*                                                                      *
************************************************************************

      use sumtallymod
      use moddas_region

*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'

*-----------------------------------------------------------------------

      common /tall00/ itnm, ital(itlmax), itals(itlmax), italm(itlmax)
      common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                itout(itlmax), ittwo(itlmax)

      common /tall02/ itrgn(itlmax), itrgm(itlmax), itreg(itlmax)
      common /tall24/ itrcm(itlmax), itrcg(itlmax)
      common /tall25/ itrnv(itlmax), itriv(itlmax), itrrv(itlmax)

      common /tall14/ ittll(itlmax), ittle(itlmax)
      character ittle*80
      common /tall31/ iterl(itlmax)

      common /tall11/ itfln(itlmax), itfll(itlmax,6), ctfln(itlmax,6)
      character ctfln*100

      common /tall19/ itsmn(itlmax), itstm(itlmax)


      common /tall68/ itstp(itlmax), itmth(itlmax), rtvr0(itlmax),
     &                rtvx0(itlmax), rtvy0(itlmax), rtvz0(itlmax),
     &                rtvx1(itlmax), rtvy1(itlmax), rtvz1(itlmax)
      common /tall60/ itger(itlmax)

      common /tall67/ itstd(itlmax), rtstd(itlmax)

      common /tall40/ rtorg(itlmax,3), rteye(itlmax,3), rtlit(itlmax,3),
     &                rtwin(itlmax,3), rtbrt(itlmax,2), ithvn(itlmax),
     &                itwin(itlmax,2), itbox(itlmax), itmir(itlmax),
     &                rtbox(itlmax,5,10), rtout(itlmax), rthet(itlmax),
     &                itlin(itlmax), itshd(itlmax), itgxs(itlmax)

      dimension idas(mdas*2) !20220906frtati temporary for fbounds-check
      equivalence ( das, idas )

*-----------------------------------------------------------------------

      common /inpec/ ititl, ipara, ibody, iregn, llarr, itby, itar

*-----------------------------------------------------------------------
       integer   icntl
       integer   inucr

       common /tcntl/  icntl, inucr
!-----------------------------------------------------------------------

      character asfil*100

      dimension ian(6)

*-----------------------------------------------------------------------
      common /tall79/ itallech
      if( itallech .eq. 0 ) return

*-----------------------------------------------------------------------
*     input echo
*-----------------------------------------------------------------------

               if( icax .eq. 0 ) then
                  write(iot,'("[ T-Volume ]")')
               else
                  write(iot,'("[ T-Volume ] off")')
               end if

                  write(iot,'("     mesh =  reg",11x,
     &            " # mesh type is region-wise")')

                     idas1 = mmmax
                     idas2 = ( idas1 + itrgn(m) - 1 ) * 2 + 1
                     idas3 = idas2 + itrgn(m)

             call echrg_wwvm(1,iterl(m),0,1,iot,itrcm(m),
     &                       idas_itrcg(itrcg(m)),
     &                       itrgn(m),itrgm(m),idas_itreg(itreg(m)),
     &                       das(idas1),idas(idas2),
     &                       itrnv(m),idas(itriv(m)),das(itrrv(m)),
     &                       idas3)

*-----------------------------------------------------------------------

                  icn = 1
                  asfil = '  # file name of output for [volume]'

            do i = 1, icn

                  iax = i

            if ( itfll(m,iax) .gt. 0 ) then

                  msfile = max( 14, itfll(m,iax) )

                  write(iot,'("     file = ",100a1)')
     &                  ( ctfln(m,iax)(j:j), j = 1, msfile ),
     &                  ( asfil(j:j), j = 1, 36 )

            end if

            end do

*-----------------------------------------------------------------------

               if( rtstd(m) .gt. 0.d0 ) then

                  write(iot,'("   stdcut =",1p1g14.7,1x,1x,
     &            " # (D=-1) threshold value of STD cutoff")')
     &            rtstd(m)

               end if

*-----------------------------------------------------------------------

                  write(iot,'("   s-type = ",i4,10x,
     &            "  # 1: Sphere source, 2: Rectangular source")')
     &            itstp(m)

               if( itstp(m) .eq. 1 ) then

                  write(iot,'("       x0 =",1p1g14.7,1x,1x,
     &            " # (D=0.0) x of sphere center")')
     &            rtvx0(m)

                  write(iot,'("       y0 =",1p1g14.7,1x,1x,
     &            " # (D=0.0) y of sphere center")')
     &            rtvy0(m)

                  write(iot,'("       z0 =",1p1g14.7,1x,1x,
     &            " # (D=0.0) z of sphere center")')
     &            rtvz0(m)

                  write(iot,'("       r0 =",1p1g14.7,1x,1x,
     &            " # radius of sphere")')
     &            rtvr0(m)

               else if( itstp(m) .eq. 2 ) then

                  write(iot,'("       x0 =",1p1g14.7,1x,1x,
     &            " #  x-min of rectangular")')
     &            rtvx0(m)

                  write(iot,'("       y0 =",1p1g14.7,1x,1x,
     &            " #  y-min of rectangular")')
     &            rtvy0(m)

                  write(iot,'("       z0 =",1p1g14.7,1x,1x,
     &            " #  z-min of rectangular")')
     &            rtvz0(m)

                  write(iot,'("       x1 =",1p1g14.7,1x,1x,
     &            " #  x-max of rectangular")')
     &            rtvx1(m)

                  write(iot,'("       y1 =",1p1g14.7,1x,1x,
     &            " #  y-max of rectangular")')
     &            rtvy1(m)

                  write(iot,'("       z1 =",1p1g14.7,1x,1x,
     &            " #  z-max of rectangular")')
     &            rtvz1(m)

               end if

*-----------------------------------------------------------------------

              if( rtout(m) .gt. 1.0 ) then

                  write(iot,'("    r-out =",1p1g14.7,1x,1x,
     &            " # (D=0:no set) radius of outer sphere")')
     &            rtout(m)

              end if

*-----------------------------------------------------------------------
              if( itmth(m) .gt. 0 ) then

                  write(iot,'("   method = ",i4,10x,
     &            "  # (D=0:track) =1:point for int. method")')
     &            itmth(m)

              end if

*-----------------------------------------------------------------------
*        sumtally subsection
*-----------------------------------------------------------------------
            if ( icntl==13 .and. nsumtalRead(m)==1 ) then

               call sumtal_echo(iot,m)

            end if

*-----------------------------------------------------------------------

            write(iot,'("#    used :",
     &                  "        main (  %)",
     &                  "        temp (  %)",
     &                  "       total (  %)"/
     &                  "#  memory :",
     &                  3(i12," (",i3,")"))')
     &                  itsmn(m), nint(dble(itsmn(m))/dble(mdas)*100),
     &                  itstm(m), nint(dble(itstm(m))/dble(mdas)*100),
     &                  itsmn(m)+itstm(m),
     &                  nint(dble(itsmn(m)+itstm(m))/dble(mdas)*100)

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine twwgech(iot,m,iaxs,icax)
*                                                                      *
*       input echo of the wwg tally                                    *
*       last modified by K.Niita on 2015/12/27                         *
*                                                                      *
************************************************************************

      use sumtallymod
      use sangelmod, only: itsans
      use partmod, only: itmxpt, itmxgp, itpan, itpat, jtpat, ipat, ipan ! frtati 2021/10/05
      use moddas_mesh
      use moddas_region
      use moddas_tally

*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'

*-----------------------------------------------------------------------

      common /tall00/ itnm, ital(itlmax), itals(itlmax), italm(itlmax)
      common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                itout(itlmax), ittwo(itlmax)

      common /tall02/ itrgn(itlmax), itrgm(itlmax), itreg(itlmax)
      common /tall24/ itrcm(itlmax), itrcg(itlmax)
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

      common /tall16/ itnfr(itlmax), itnfz(itlmax), itnfn(itlmax)

      common /tall14/ ittll(itlmax), ittle(itlmax)
      character ittle*80

      common /tall15/ itanl(itlmax), itang(itlmax)
      character itang*200

      common /tall27/ itaxl(itlmax), itaxt(itlmax),
     &                itayl(itlmax), itayt(itlmax),
     &                itazl(itlmax), itazt(itlmax)
      character itaxt*200, itayt*200, itazt*200

      common /tall19/ itsmn(itlmax), itstm(itlmax)

      common /tall23/ itaty(itlmax), itanm(itlmax), itarg(itlmax),
     &                rtami(itlmax), rtama(itlmax), rtadl(itlmax)


      common /tall20/ ittty(itlmax), ittnm(itlmax), ittrg(itlmax),
     &                rttmi(itlmax), rttma(itlmax), rttdl(itlmax)
      common /tall21/ rtfac(itlmax)
      common /tall28/ itgsh(itlmax)
      common /tall29/ itrsh(itlmax)
      common /tall31/ iterl(itlmax)
      common /tall22/ itmcn(itlmax), itmtn(itlmax), itmtt(itlmax)
      common /tall34/ itvm(itlmax)
      common /tall35/ iteps(itlmax)
      common /tall37/ itcnt(9,itlmax)
      common /cntreg/ ncntc(3), ncreg(3), incrc(3), incrt(3), kcont(3)
      common /cntpat/ icpan(3), icpat(3,20,2)

      common /tall38/ itres(itlmax)
      common /tall39/ rtwid(itlmax)

      common /tall45/ itmlp(itlmax), itmln(itlmax,6), itmst(itlmax),
     &                itmli(itlmax,6), rtmme(itlmax,6), itmnt(itlmax,6),
     &                itmpn(itlmax,6), itmpt(itlmax,6,6,2)
      common /tall46/ itmtr(itlmax,4), rtmtr(itlmax,13)
      common /tall49/ itglt(itlmax)

      common /tall60/ itger(itlmax)

      common /tall67/ itstd(itlmax), rtstd(itlmax)

      common /tall82/ itcnth(9,itlmax)
cfrtati 2023/12/07
      common /tall84/ itmto(itlmax,6,4)

! T.Sato 2024/03/18 weighted history counter ID
      common /tall92/ ichnum(itlmax),chbias(itlmax),pedest(itlmax)
     &               ,ictnum(itlmax),ctbias(itlmax),ictidx(itlmax)

      common /tall94/ normwwg(itlmax),wwmax(itlmax),ichplane(itlmax) ! twwg normalization parameter
     &               ,elowthre(itlmax),elowbias(itlmax)  ! T.Sato 2025/01/04

      dimension idas(mdas*2) !20220906frtati temporary for fbounds-check
      equivalence ( das, idas )


cfrtati 2021/10/05 moved to partmod

      dimension ilgt(6), jlgt(6)
      character ch13*13, ch6*6
      character ch200*400
      character ch100*400

*-----------------------------------------------------------------------

      common /inpec/ ititl, ipara, ibody, iregn, llarr, itby, itar

*-----------------------------------------------------------------------
       integer   icntl
       integer   inucr

       common /tcntl/  icntl, inucr
!-----------------------------------------------------------------------

      character asfil*100

      character aname(11)*3
      data aname / 'eng','reg',' xy',' yz',' xz','  t','wwg',
     &             '  x','  y','  z','tet'/

      dimension ian(6)
      dimension jmat(20)

*-----------------------------------------------------------------------
      common /tall79/ itallech
      if( itallech .eq. 0 ) return

*-----------------------------------------------------------------------
*     input echo
*-----------------------------------------------------------------------

                  write(iot,'("[ T-WWG ]")')

                  write(iot,'("    title = ",80a1)')
     &            ( ittle(m)(i:i), i = 1, ittll(m) )

            if( itmsh(m) .eq. 1 ) then

                  write(iot,'("     mesh =  reg",11x,
     &            " # mesh type is region-wise")')

                     idas1 = mmmax
                     idas2 = ( idas1 + itrgn(m) - 1 ) * 2 + 1
                     idas3 = idas2 + itrgn(m)

                  call echrg(1,iterl(m),0,1,iot,itrcm(m),
     &                       idas_itrcg(itrcg(m)),
     &                       itrgn(m),itrgm(m),idas_itreg(itreg(m)),
     &                       das(idas1),idas(idas2),
     &                       itrnv(m),idas(itriv(m)),das(itrrv(m)),
     &                       idas3)

            else if( itmsh(m) .eq. 3 ) then

                  write(iot,'("     mesh =  xyz",11x,
     &            " # mesh type is xyz scoring mesh")')

                  call echmty(0,iot,'x',itxty(m),itxnm(m),
     &                        rtxdl(m),rtxmi(m),rtxma(m),itxrg(m)
     &                        ,abs(itxnm(m)+1),das_itxrg)
                  call echmty(0,iot,'y',ityty(m),itynm(m),
     &                        rtydl(m),rtymi(m),rtyma(m),ityrg(m)
     &                        ,abs(itynm(m)+1),das_ityrg)
                  call echmty(0,iot,'z',itzty(m),itznm(m),
     &                        rtzdl(m),rtzmi(m),rtzma(m),itzrg(m)
     &                        ,abs(itznm(m)+1),das_itzrg)

            else if( itmsh(m) .eq. 4 ) then

                  write(iot,'("     mesh =  tet",11x,
     &            " # mesh type is tetra geometry")')

                  call echtet(iot,itrgm(m),idas_itreg(itreg(m)))

            end if

*-----------------------------------------------------------------------

            if( itmtr(m,1) .gt. 0 ) then

               if( itmtr(m,1) .eq. 1 ) then
                     write(iot,'("     trcl = ",i6,9x,
     &               " # transform ID")') itmtr(m,3)

               else

                  if( itmtr(m,2) .eq. 0 ) then

                     write(iot,'("     trcl = ",
     &                           3(1p1e15.7),
     &                           3(/12x,3(1p1e15.7)),i5)')
     &               ( rtmtr(m,j), j = 1, 12 ),
     &               nint( rtmtr(m,13) )

                  else

                     write(iot,'("    *trcl = ",
     &                           3(1p1e15.7),
     &                           3(/12x,3(1p1e15.7)),i5)')
     &               ( rtmtr(m,j), j = 1, 12 ),
     &               nint( rtmtr(m,13) )

                  end if

               end if

            end if

*-----------------------------------------------------------------------

                  call echmty(0,iot,'e',itety(m),itenm(m),
     &                        rtedl(m),rtemi(m),rtema(m),iterg(m)
     &                        ,abs(itenm(m)+1),das_iterg)

*-----------------------------------------------------------------------

               if( ittty(m) .ne. 0 ) then

                  call echmty(0,iot,'t',ittty(m),ittnm(m),
     &                        rttdl(m),rttmi(m),rttma(m),ittrg(m)
     &                        ,abs(ittnm(m)+1),das_ittrg)

               end if

*-----------------------------------------------------------------------

                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [1/cm^2/source]")') itunt(m)

*-----------------------------------------------------------------------

                  asfil = '  # file name of output for the above axis'

            if( icax .eq. 0 ) then

                  icn = itaxn(m)

               do i = 1, icn

                  ian(i) = i

               end  do

            else

                  icn = 1
                  ian(1) = iaxs

            end if

*-----------------------------------------------------------------------

               if( itmtn(m) .eq. 0 ) then

                  write(iot,'(" material =  all"11x,
     &            " # (D=all) number of specific material")')

               else if( itmtn(m) .gt. 0 ) then

                  write(iot,'(" material =",i5,11x,
     &            " # (D=all) number of specific material")')
     &            itmtn(m) * itmcn(m)

                     k = 0

                  do j = 1, itmtn(m)

                     k = k + 1

                     jmat(k) = ismte( itmtt(m) + j - 1 )

                     if( k / 10 * 10 .eq. k .or. j .eq. itmtn(m) ) then

                        write(iot,'(11x,10(i6))') ( jmat(i),i = 1, k )

                        k = 0

                     end if

                  end do

               end if

*-----------------------------------------------------------------------

               idtyp = 0

            do i = 1, icn

               if( itaxs(m,ian(i)) .ge. 7 .and.
     &             itaxs(m,ian(i)) .le. 10 )  idtyp = ittwo(m)

            end do

            if( idtyp .ne. 0 ) then

               write(iot,'("  2D-type = ",i4,11x,
     &         " # 1:Cont, 2:Clust, 3:Color, 4:xyz, 5:mat,",
     &         " 6:Clust+Cont, 7:Col+Cont")') idtyp

            end if

*-----------------------------------------------------------------------

            do i = 1, icn

                  iax = ian(i)

                  write(iot,'("     axis =  ",a3,11x,
     &            " # axis of output")')
     &            aname( itaxs(m,iax) )

            if ( itfll(m,iax) .gt. 0 ) then

                  msfile = max( 14, itfll(m,iax) )

                  write(iot,'("     file = ",100a1)')
     &                  ( ctfln(m,iax)(j:j), j = 1, msfile ),
     &                  ( asfil(j:j), j = 1, 42 )

            end if

            end do
               call echprt(1,iot,m,itpan,itpat,jtpat,itnm,itmxpt,itmxgp)

*-----------------------------------------------------------------------
*        multiplier
*-----------------------------------------------------------------------

         if( itmlp(m) .gt. 0 ) then

            do k = 1, itmlp(m)

*-----------------------------------------------------------------------

               if( itmln(m,k) .eq. 1 .and.
     &             mltp(1,itmli(m,k)/13+1) .eq. -1 ) then

                  imall = 0

                  write(iot,'("  multiplier = all",9x,
     &            " # number of material group")')

               else

                  imall = itmln(m,k)

                  write(iot,'("  multiplier = ",i3,9x,
     &            " # number of material group")') imall

               end if

*-----------------------------------------------------------------------

                     ipan(m) = itmpn(m,k)

                  do j = 1, itmpn(m,k)

                     ipat(m,j,1) = itmpt(m,k,j,1)
                     ipat(m,j,2) = itmpt(m,k,j,2)

                  end do
                 call echprt(1,iot,m,ipan,ipat,jtpat,itnm,itmxpt,itmxgp)

               if( rtmme(m,k) .gt. 0.0 ) then

                  write(iot,'("     emax = ",1p1g13.5)') rtmme(m,k)

               end if

*-----------------------------------------------------------------------

                     do i = 1, itmst(m)

                        ilgt(i) = 0

                     end do

               do j = 1, itmln(m,k)

                  do im = 1, itmst(m)

                        nn = im * 2 - 1
                        jlgt(im) = 0

                  do jj = mltp(1+nn,itmli(m,k)/13+j),
     &                    mltp(1+nn+1,itmli(m,k)/13+j), 2

                        if( nint(slib(jj)) .eq. 1 .or.
     &                      nint(slib(jj)) .eq. 2 .or.
     &                      nint(slib(jj)) .eq. 3 ) then

                           jlgt(im) = jlgt(im) + 2

                        else if( nint(slib(jj)) .eq. 4 .or.
     &                           nint(slib(jj)) .eq. 8 ) then

                           write(ch13,'(1p1g13.5)') slib(jj+1)

                           do l = 1, 13
                              if( ch13(l:l) .ne. ' ' ) goto 410
                           end do
  410                      isi = l - 1
                           do l = 13, 1, -1
                              if( ch13(l:l) .ne. ' ' ) goto 420
                           end do
  420                      isf = 13 - l

                           jlgt(im) = jlgt(im) + 13 - isi - isf + 1

                        else if( nint(slib(jj)) .ge. 5 .and.
     &                           nint(slib(jj)) .le. 7 ) then

                           write(ch6,'(i6)') nint( slib(jj+1) )

                           do l = 1, 6
                              if( ch6(l:l) .ne. ' ' ) goto 430
                           end do
  430                      isi = l - 1

                           jlgt(im) = jlgt(im) + 6 - isi + 1

                        end if

                  end do
                  end do

                     do i = 1, itmst(m)

                        if( jlgt(i) .gt. ilgt(i) ) ilgt(i) = jlgt(i)

                     end do

               end do

                     do i = 1, itmst(m)

                        if( ilgt(i) .lt. 7 ) ilgt(i) = 7

                     end do

                  ch200(1:9) = '      mat'
                  ic = 9

               do i = 1, itmst(m)

                  isp = ilgt(i) - 7

                  do j = 1, isp
                     ch200(ic+j:ic+j) = ' '
                  end do
                     ic = ic + isp
                     ch200(ic+1:ic+4) = 'mset'
                     ic = ic + 4

                  write(ch200(ic+1:ic+1),'(i1)') itmnt(m,i)
                     ic = ic + 1

                  ch200(ic+1:ic+2) = '  '
                     ic = ic + 2

               end do

                  write(iot,'(200a1)') (ch200(j:j),j=1,ic)

*-----------------------------------------------------------------------

               do j = 1, itmln(m,k)

                     do ik = 1, 400
                        ch100(ik:ik) = ' '
                     end do

                     ch200(1:3) = '   '
                     ic = 3

                  if( mltp(1,itmli(m,k)/13+j) .gt. 0 ) then
                     write(ch200(ic+1:ic+6),'(i6)')
     &                  mltp(1,itmli(m,k)/13+j)

                   else

                     ch200(ic+1:ic+6) = '   all'

                   end if

                     ic = ic + 6

                  do im = 1, itmst(m)

                        nn = im * 2 - 1
                        id = 0

                  do jj = mltp(1+nn,itmli(m,k)/13+j),
     &                    mltp(1+nn+1,itmli(m,k)/13+j), 2

                        if( nint(slib(jj)) .eq. 1 ) then

                           ch100(id+1:id+2) = ' ('
                           id = id + 2

                        else if( nint(slib(jj)) .eq. 2 ) then

                           ch100(id+1:id+2) = ' )'
                           id = id + 2

                        else if( nint(slib(jj)) .eq. 3 ) then

                           ch100(id+1:id+2) = ' :'
                           id = id + 2

                        else if( nint(slib(jj)) .eq. 4 .or.
     &                           nint(slib(jj)) .eq. 8 ) then

                           write(ch13,'(1p1g13.5)') slib(jj+1)

                           do l = 1, 13
                              if( ch13(l:l) .ne. ' ' ) goto 440
                           end do
  440                      isi = l
                           do l = 13, 1, -1
                              if( ch13(l:l) .ne. ' ' ) goto 450
                           end do
  450                      isf = l

                           ch100(id+1:id+1) = ' '
                           id = id + 1
                           ch100(id+1:id+isf-isi+1) = ch13(isi:isf)
                           id = id+isf-isi+1

                        else if( nint(slib(jj)) .ge. 5 .and.
     &                           nint(slib(jj)) .le. 7 ) then

                           write(ch6,'(i6)') nint( slib(jj+1) )

                           do l = 1, 6
                              if( ch6(l:l) .ne. ' ' ) goto 460
                           end do
  460                      isi = l

                           ch100(id+1:id+1) = ' '
                           id = id + 1
                           ch100(id+1:id+6-isi+1) = ch6(isi:6)
                           id = id+6-isi+1

                        end if

                  end do

                           isd = ilgt(im) - id

                           do kk = 1, isd
                              ch200(ic+kk:ic+kk) = ' '
                           end do
                           ic = ic + isd

                           ch200(ic+1:ic+id) = ch100(1:id)
                           ic = ic + id

                  end do

                  write(iot,'(200a1)') (ch200(kk:kk),kk=1,ic)

               end do

*-----------------------------------------------------------------------

            end do

         end if

*-----------------------------------------------------------------------

               if( rtstd(m) .gt. 0.d0 ) then

                  write(iot,'("   stdcut =",1p1g14.7,1x,1x,
     &            " # (D=-1) threshold value of STD cutoff")')
     &            rtstd(m)

               end if

*-----------------------------------------------------------------------

               if( rtfac(m) .ne. 1.0d0 ) then

                  write(iot,'("   factor =",1p1g14.7,1x,1x,
     &            " # (D=1.0) normalization factor")')
     &            rtfac(m)

               end if

*-----------------------------------------------------------------------

               if( itmtn(m) .gt. 0 .and. itmsh(m) .eq. 3 ) then

                  write(iot,'("   volmat = ",i4,11x,
     &            " # (D=9) vol.corr., 0:no, number of scan / dir.")')
     &            itvm(m)

               end if

*-----------------------------------------------------------------------

               if( itanl(m) .gt. 0 ) then

                  write(iot,'("    angel = ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               end if

*-----------------------------------------------------------------------

               if( itsans(m) .gt. 0 ) then

                  write(iot,'("   sangel = ",i2)') itsans(m)
                  call write_sangel(iot,m,1)

               end if

*-----------------------------------------------------------------------

               if( itaxl(m) .gt. 0 ) then

                  write(iot,'("    x-txt = ",200a1)')
     &            ( itaxt(m)(i:i),i = 1, itaxl(m) )

               end if

*-----------------------------------------------------------------------

               if( itayl(m) .gt. 0 ) then

                  write(iot,'("    y-txt = ",200a1)')
     &            ( itayt(m)(i:i),i = 1, itayl(m) )

               end if

*-----------------------------------------------------------------------

               if( itazl(m) .gt. 0 ) then

                  write(iot,'("    z-txt = ",200a1)')
     &            ( itazt(m)(i:i),i = 1, itazl(m) )

               end if

*-----------------------------------------------------------------------

               if( itgsh(m) .gt. 0 ) then

                  write(iot,'("    gshow = ",i4,10x,
     &            "  # 0: no 1:bnd, 2:bnd+mat, ",
     &     "3:bnd+reg 4:bnd+lat 5:bmp style")')
     &            itgsh(m)

               end if

*-----------------------------------------------------------------------
               if( itglt(m) .ge. 0 ) then

                  write(iot,'("    gslat = ",i4,10x,
     &            "  # 0: no line for lat, 1: with,",
     &               " 2: no line between the same cell")')
     &            itglt(m)

               end if

*-----------------------------------------------------------------------

               if( itgsh(m) .gt. 0 .and. itger(m) .gt. 0 ) then

                  write(iot,'("    ginfo = ",i4,11x,
     &            " # (D=0) 0:no, 1:show in graph, 2:add file")')
     &            itger(m)

               end if

*-----------------------------------------------------------------------

               if( itrsh(m) .gt. 0 ) then

                  write(iot,'("    rshow = ",i4,10x,
     &            "  # 0: no 1:bnd, 2:bnd+mat, ",
     &     "3:bnd+reg 4:bnd+lat 5:bmp style")')
     &            itrsh(m)

                  call echmty(1,iot,'x',itxty(m),itxnm(m),
     &                        rtxdl(m),rtxmi(m),rtxma(m),itxrg(m)
     &                        ,abs(itxnm(m)+1),das_itxrg)
                  call echmty(1,iot,'y',ityty(m),itynm(m),
     &                        rtydl(m),rtymi(m),rtyma(m),ityrg(m)
     &                        ,abs(itynm(m)+1),das_ityrg)
                  call echmty(1,iot,'z',itzty(m),itznm(m),
     &                        rtzdl(m),rtzmi(m),rtzma(m),itzrg(m)
     &                        ,abs(itznm(m)+1),das_itzrg)

               end if

*-----------------------------------------------------------------------

               if( itgsh(m) .gt. 0 .or. itrsh(m) .gt. 0 ) then

                  write(iot,'("    resol = ",i4,10x,
     &            "  # (D=1) resolution of gshow or rshow")')
     &            itres(m)

                  write(iot,'("    width =",1p1g14.7,1x,1x,
     &            " # (D=0.5) width of lines for gshow or rshow")')
     &            rtwid(m)

               end if

*-----------------------------------------------------------------------

               if( iterl(m) .ne. 72 ) then

                  write(iot,'("   iechrl =",i5,10x,
     &            "  # (D=72) length reg of val or vol")')
     &            iterl(m)

               end if

*-----------------------------------------------------------------------

               if( iteps(m) .ne. 0 ) then

                  write(iot,'("   epsout =",i5,10x,
     &            "  # (D=0) generate eps file by ANGEL")')
     &            iteps(m)

               end if

*-----------------------------------------------------------------------
*        check of counter
*-----------------------------------------------------------------------

            do i = 1, 3

               if( itcnt(i,m) .ne. 0 ) then

                  write(iot,'(" ctmin(",i1,") =",i5,10x,
     &            "  # (D=-9999) min. of counter ",i1)')
     &            i, itcnt(i*2+2,m), i

                  write(iot,'(" ctmax(",i1,") =",i5,10x,
     &            "  # (D= 9999) max. of counter ",i1)')
     &            i, itcnt(i*2+3,m), i

                  if( ncntc(i) .eq. 0 )
     &            write(iot,'("# Warning: there is no counter = ",
     &                        i1," definition in [counter].")') i

               end if

            end do

            do i = 1, 3
               if( itcnth(i,m) .ne. 0 ) then
                  write(iot,'(" chmin(",i1,") =",i5,10x,
     &            "  # (D=-9999) min. of history counter ",i1)')
     &            i, itcnth(i*2+2,m), i
                  write(iot,'(" chmax(",i1,") =",i5,10x,
     &            "  # (D= 9999) max. of history counter ",i1)')
     &            i, itcnth(i*2+3,m), i
               end if
            end do

! T.Sato 2024/04/04, chwei parameter
           if(ichnum(m).ne.0) then ! chwei is specified
             write(iot,'(" chwei(",i1,") =",f8.2,7x,
     &       "  # (D=0) history counter bias factor")')
     &       ichnum(m), chbias(m)
           end if

! T.Sato 2024/05/18, ctwei parameter
           if(ictnum(m).ne.0) then ! ctwei is specified
             write(iot,'(" ctwei(",i1,") =",f8.2,7x,
     &       "  # (D=0) counter bias factor")')
     &       ictnum(m), ctbias(m)
           end if

! T.Sato 2024/05/28, pedestal parameter
           if(pedest(m).ne.0.1d0) then ! pedestal is specified
             write(iot,'(" pedestal =",f8.2,7x,
     &       "  # (D=0.1) pedestal weight from the bottom")')
     &       pedest(m)
           end if

! T.Sato 2024/06/09, normwwg parameter
           if(normwwg(m).ne.0) then ! normwwg is specified
             write(iot,'(" normww   =",i8,7x,
     &       "  # (D=0) normalization method")')
     &       normwwg(m)
           end if

! T.Sato 2024/06/09, wwmax
           if(wwmax(m).ne.1.0d0) then ! normwwg is specified
             write(iot,'(" wwmax    =",es15.4,
     &       "  # (D=1.0) maximum lower weight window")')
     &       wwmax(m)
           end if

! T.Sato 2024/06/13, chplane
           if(ichplane(m).eq.1) then ! xy
             write(iot,'(" chplane   =             xy",
     &       "  # (D=all) xy, yz, or xz plane for chwei")')
           elseif(ichplane(m).eq.2) then ! yz
             write(iot,'(" chplane   =             yz",
     &       "  # (D=all) xy, yz, or xz plane for chwei")')
           elseif(ichplane(m).eq.3) then ! xz
             write(iot,'(" chplane   =             xz",
     &       "  # (D=all) xy, yz, or xz plane for chwei")')
           endif

C S.Hashimoto, 2014.11.27
*-----------------------------------------------------------------------
*        sumtally subsection
*-----------------------------------------------------------------------
            if ( icntl==13 .and. nsumtalRead(m)==1 ) then

               call sumtal_echo(iot,m)

            end if

*-----------------------------------------------------------------------

            write(iot,'("#    used :",
     &                  "        main (  %)",
     &                  "        temp (  %)",
     &                  "       total (  %)"/
     &                  "#  memory :",
     &                  3(i12," (",i3,")"))')
     &                  itsmn(m), nint(dble(itsmn(m))/dble(mdas)*100),
     &                  itstm(m), nint(dble(itstm(m))/dble(mdas)*100),
     &                  itsmn(m)+itstm(m),
     &                  nint(dble(itsmn(m)+itstm(m))/dble(mdas)*100)

*-----------------------------------------------------------------------
cfrtati 2024/02/13
            if( itmlp(m).gt.0 ) then
              call displaymtinfo(iot,m)
            end if

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine tponech(iot,m,iaxs,icax)
*                                                                      *
*       input echo of the point tally                                  *
*       last modified by K.Niita on 2015/03/11                         *
*                                                                      *
************************************************************************

      use sumtallymod
      use sangelmod, only: itsans
      use partmod, only: itmxpt, itmxgp, itpan, itpat, jtpat, ipat, ipan ! frtati 2021/10/05
      use moddas_mesh
      use moddas_tally
*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'

*-----------------------------------------------------------------------

      common /tall00/ itnm, ital(itlmax), itals(itlmax), italm(itlmax) ! frtati 2021/10/05
      common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                itout(itlmax), ittwo(itlmax)
      common /tall07/ itety(itlmax), itenm(itlmax), iterg(itlmax),
     &                rtemi(itlmax), rtema(itlmax), rtedl(itlmax)

      common /tall10/ itaxn(itlmax), itaxs(itlmax,6)
      common /tall11/ itfln(itlmax), itfll(itlmax,6), ctfln(itlmax,6)
      character ctfln*100

      common /tall14/ ittll(itlmax), ittle(itlmax)
      character ittle*80

      common /tall15/ itanl(itlmax), itang(itlmax)
      character itang*200

      common /tall27/ itaxl(itlmax), itaxt(itlmax),
     &                itayl(itlmax), itayt(itlmax),
     &                itazl(itlmax), itazt(itlmax)
      character itaxt*200, itayt*200, itazt*200

      common /tall19/ itsmn(itlmax), itstm(itlmax)

      common /tall20/ ittty(itlmax), ittnm(itlmax), ittrg(itlmax),
     &                rttmi(itlmax), rttma(itlmax), rttdl(itlmax)
      common /tall21/ rtfac(itlmax)
      common /tall35/ iteps(itlmax)
      common /tall37/ itcnt(9,itlmax)
      common /cntreg/ ncntc(3), ncreg(3), incrc(3), incrt(3), kcont(3)

      common /tall45/ itmlp(itlmax), itmln(itlmax,6), itmst(itlmax),
     &                itmli(itlmax,6), rtmme(itlmax,6), itmnt(itlmax,6),
     &                itmpn(itlmax,6), itmpt(itlmax,6,6,2)

      common /tall62/ itpon(itlmax), rtpon(itlmax,20,4)

      common /tall67/ itstd(itlmax), rtstd(itlmax)

      common /tall82/ itcnth(9,itlmax)
cfrtati 2023/12/07
      common /tall84/ itmto(itlmax,6,4)

      common /tall91/ itextstat(itlmax), mftal(itlmax) ! S.H. extstat 2024.3.18
      common /tall93/ tprodenmn(itlmax), tprodenmx(itlmax),
     &     nbtproden(itlmax), itprodenchk(itlmax) !S.H. extstat 2024.4.28

*-----------------------------------------------------------------------
      common /tall78/ itism(itlmax), itist(10,itlmax), itjst(10,itlmax),
     &                itkst(10,itlmax), itstt(itlmax), itsdd(itlmax)

      character schan(10)*5
      data (schan(i), i = 1, 10 ) /
     &    'ireg ','ix   ','iy   ','iz   ','ir   ',
     &    'ie   ','it   ','ia   ','ipart','imul '/
      character ch73*200

      character yen*1

      dimension idas(mdas*2) ! S.H. 2022.12.16 temporary for fbounds-check
      equivalence ( das, idas )

*-----------------------------------------------------------------------


cfrtati 2021/10/05 moved to partmod

      dimension ilgt(6), jlgt(6)
      character ch13*13, ch6*6
      character ch200*400
      character ch100*400

*-----------------------------------------------------------------------
       integer   icntl
       integer   inucr

       common /tcntl/  icntl, inucr
!-----------------------------------------------------------------------

      character asfil*100

      character aname(2)*3
      data aname / 'eng','  t'/

      dimension ian(6)

      integer           iMeVperu
      common /cMeVperu/ iMeVperu

*-----------------------------------------------------------------------
      common /tall79/ itallech
      if( itallech .eq. 0 ) return

*-----------------------------------------------------------------------
*     input echo
*-----------------------------------------------------------------------

            write(iot,'("[ T-Point ]")')

                  write(iot,'("    title = ",80a1)')
     &            ( ittle(m)(i:i), i = 1, ittll(m) )

*-----------------------------------------------------------------------
*           point
*-----------------------------------------------------------------------

            if( itpon(m) .eq. 1 ) then

                  write(iot,'("    point = ",i4,11x,
     &            " # number of point estimators")') itmsh(m)

                  write(iot,'("  non      x             y",
     &             "             z            r0")')

               do i = 1, itmsh(m)

                  write(iot,'(i5,3x,4(1p1g14.7))') i,
     &                       (rtpon(m,i,j),j=1,4)

               end do

*-----------------------------------------------------------------------
*           ring
*-----------------------------------------------------------------------

            else if( itpon(m) .eq. 2 ) then

                  write(iot,'("     ring = ",i4,11x,
     &            " # number of ring estimators")') itmsh(m)

                  write(iot,'("  non   axis     ar",
     &              "           rr            r0")')

               do i = 1, itmsh(m)

                  if( nint(rtpon(m,i,1)) .eq. 1 ) then
                     write(iot,'(i5,4x,"x",4x,3(1p1g14.7))') i,
     &                       (rtpon(m,i,j),j=2,4)
                  else if( nint(rtpon(m,i,1)) .eq. 2 ) then
                     write(iot,'(i5,4x,"y",4x,3(1p1g14.7))') i,
     &                       (rtpon(m,i,j),j=2,4)
                  else if( nint(rtpon(m,i,1)) .eq. 3 ) then
                     write(iot,'(i5,4x,"z",4x,3(1p1g14.7))') i,
     &                       (rtpon(m,i,j),j=2,4)
                  end if

               end do

            end if

*-----------------------------------------------------------------------

                  call echmty(0,iot,'e',itety(m),itenm(m),
     &                        rtedl(m),rtemi(m),rtema(m),iterg(m)
     &                        ,abs(itenm(m)+1),das_iterg)

*-----------------------------------------------------------------------

               if( ittty(m) .ne. 0 ) then

                  call echmty(0,iot,'t',ittty(m),ittnm(m),
     &                        rttdl(m),rttmi(m),rttma(m),ittrg(m)
     &                        ,abs(ittnm(m)+1),das_ittrg)

               end if

*-----------------------------------------------------------------------

               if( itunt(m) .eq. 1 ) then

                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [1/cm^2/source]")') itunt(m)

               else if( itunt(m) .eq. 2 ) then

                if ( iMeVperu.eq. 1 ) then
                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [1/cm^2/(MeV/n)/source]")') itunt(m)
                else
                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [1/cm^2/MeV/source]")') itunt(m)
                end if

               else if( itunt(m) .eq. 3 ) then

                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [1/cm^2/Lethargy/source]")') itunt(m)

               else if( itunt(m) .eq. 11 ) then

                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [1/cm^2/nsec/source]")') itunt(m)

               else if( itunt(m) .eq. 12 ) then

                if ( iMeVperu.eq. 1 ) then
                  write(iot,'("     unit = ",i4,11x,
     &           " # unit is [1/cm^2/(MeV/n)/nsec/source]")') itunt(m)
                else
                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [1/cm^2/MeV/nsec/source]")') itunt(m)
                end if

               else if( itunt(m) .eq. 13 ) then

                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [1/cm^2/Lethargy/nsec/source]")')
     &               itunt(m)

               end if

*-----------------------------------------------------------------------

                  asfil = '  # file name of output for the above axis'

            if( icax .eq. 0 ) then

                  icn = itaxn(m)

               do i = 1, icn

                  ian(i) = i

               end  do

            else

                  icn = 1
                  ian(1) = iaxs

            end if

*-----------------------------------------------------------------------

            do i = 1, icn

                  iax = ian(i)

                  write(iot,'("     axis =  ",a3,11x,
     &            " # axis of output")')
     &            aname( itaxs(m,iax) )

             if ( itfll(m,iax) .gt. 0 ) then

                  msfile = max( 14, itfll(m,iax) )

                  write(iot,'("     file = ",100a1)')
     &                  ( ctfln(m,iax)(j:j), j = 1, msfile ),
     &                  ( asfil(j:j), j = 1, 42 )

            end if

            end do

*-----------------------------------------------------------------------
               call echprt(1,iot,m,itpan,itpat,jtpat,itnm,itmxpt,itmxgp)

*-----------------------------------------------------------------------
*        multiplier
*-----------------------------------------------------------------------

         if( itmlp(m) .gt. 0 ) then

            do k = 1, itmlp(m)

*-----------------------------------------------------------------------

               if( itmln(m,k) .eq. 1 .and.
     &             mltp(1,itmli(m,k)/13+1) .eq. -1 ) then

                  imall = 0

                  write(iot,'("  multiplier = all",9x,
     &            " # number of material group")')

               else

                  imall = itmln(m,k)

                  write(iot,'("  multiplier = ",i3,9x,
     &            " # number of material group")') imall

               end if

*-----------------------------------------------------------------------

                     ipan(m) = itmpn(m,k)

                  do j = 1, itmpn(m,k)

                     ipat(m,j,1) = itmpt(m,k,j,1)
                     ipat(m,j,2) = itmpt(m,k,j,2)

                  end do

                 call echprt(1,iot,m,ipan,ipat,jtpat,itnm,itmxpt,itmxgp)

               if( rtmme(m,k) .gt. 0.0 ) then

                  write(iot,'("     emax = ",1p1g13.5)') rtmme(m,k)

               end if

*-----------------------------------------------------------------------

                     do i = 1, itmst(m)

                        ilgt(i) = 0

                     end do

               do j = 1, itmln(m,k)

                  do im = 1, itmst(m)

                        nn = im * 2 - 1
                        jlgt(im) = 0

                  do jj = mltp(1+nn,itmli(m,k)/13+j),
     &                    mltp(1+nn+1,itmli(m,k)/13+j), 2

                        if( nint(slib(jj)) .eq. 1 .or.
     &                      nint(slib(jj)) .eq. 2 .or.
     &                      nint(slib(jj)) .eq. 3 ) then

                           jlgt(im) = jlgt(im) + 2

                        else if( nint(slib(jj)) .eq. 4 .or.
     &                           nint(slib(jj)) .eq. 8 ) then

                           write(ch13,'(1p1g13.5)') slib(jj+1)

                           do l = 1, 13
                              if( ch13(l:l) .ne. ' ' ) goto 410
                           end do
  410                      isi = l - 1
                           do l = 13, 1, -1
                              if( ch13(l:l) .ne. ' ' ) goto 420
                           end do
  420                      isf = 13 - l

                           jlgt(im) = jlgt(im) + 13 - isi - isf + 1

                        else if( nint(slib(jj)) .ge. 5 .and.
     &                           nint(slib(jj)) .le. 7 ) then

                           write(ch6,'(i6)') nint( slib(jj+1) )

                           do l = 1, 6
                              if( ch6(l:l) .ne. ' ' ) goto 430
                           end do
  430                      isi = l - 1

                           jlgt(im) = jlgt(im) + 6 - isi + 1

                        end if

                  end do
                  end do

                     do i = 1, itmst(m)

                        if( jlgt(i) .gt. ilgt(i) ) ilgt(i) = jlgt(i)

                     end do

               end do

                     do i = 1, itmst(m)

                        if( ilgt(i) .lt. 7 ) ilgt(i) = 7

                     end do

                  ch200(1:9) = '      mat'
                  ic = 9

               do i = 1, itmst(m)

                  isp = ilgt(i) - 7

                  do j = 1, isp
                     ch200(ic+j:ic+j) = ' '
                  end do
                     ic = ic + isp
                     ch200(ic+1:ic+4) = 'mset'
                     ic = ic + 4

                  write(ch200(ic+1:ic+1),'(i1)') itmnt(m,i)
                     ic = ic + 1

                  ch200(ic+1:ic+2) = '  '
                     ic = ic + 2

               end do

                  write(iot,'(200a1)') (ch200(j:j),j=1,ic)

*-----------------------------------------------------------------------

               do j = 1, itmln(m,k)

                     do ik = 1, 400
                        ch100(ik:ik) = ' '
                     end do

                     ch200(1:3) = '   '
                     ic = 3

                  if( mltp(1,itmli(m,k)/13+j) .gt. 0 ) then
                     write(ch200(ic+1:ic+6),'(i6)')
     &                  mltp(1,itmli(m,k)/13+j)

                   else

                     ch200(ic+1:ic+6) = '   all'

                   end if

                     ic = ic + 6

                  do im = 1, itmst(m)

                        nn = im * 2 - 1
                        id = 0

                  do jj = mltp(1+nn,itmli(m,k)/13+j),
     &                    mltp(1+nn+1,itmli(m,k)/13+j), 2

                        if( nint(slib(jj)) .eq. 1 ) then

                           ch100(id+1:id+2) = ' ('
                           id = id + 2

                        else if( nint(slib(jj)) .eq. 2 ) then

                           ch100(id+1:id+2) = ' )'
                           id = id + 2

                        else if( nint(slib(jj)) .eq. 3 ) then

                           ch100(id+1:id+2) = ' :'
                           id = id + 2

                        else if( nint(slib(jj)) .eq. 4 .or.
     &                           nint(slib(jj)) .eq. 8 ) then

                           write(ch13,'(1p1g13.5)') slib(jj+1)

                           do l = 1, 13
                              if( ch13(l:l) .ne. ' ' ) goto 440
                           end do
  440                      isi = l
                           do l = 13, 1, -1
                              if( ch13(l:l) .ne. ' ' ) goto 450
                           end do
  450                      isf = l

                           ch100(id+1:id+1) = ' '
                           id = id + 1
                           ch100(id+1:id+isf-isi+1) = ch13(isi:isf)
                           id = id+isf-isi+1

                        else if( nint(slib(jj)) .ge. 5 .and.
     &                           nint(slib(jj)) .le. 7 ) then

                           write(ch6,'(i6)') nint( slib(jj+1) )

                           do l = 1, 6
                              if( ch6(l:l) .ne. ' ' ) goto 460
                           end do
  460                      isi = l

                           ch100(id+1:id+1) = ' '
                           id = id + 1
                           ch100(id+1:id+6-isi+1) = ch6(isi:6)
                           id = id+6-isi+1

                        end if

                  end do

                           isd = ilgt(im) - id

                           do kk = 1, isd
                              ch200(ic+kk:ic+kk) = ' '
                           end do
                           ic = ic + isd

                           ch200(ic+1:ic+id) = ch100(1:id)
                           ic = ic + id

                  end do

                  write(iot,'(200a1)') (ch200(kk:kk),kk=1,ic)

               end do

*-----------------------------------------------------------------------

            end do

         end if

*-----------------------------------------------------------------------

               if( rtfac(m) .ne. 1.0d0 ) then

                  write(iot,'("   factor =",1p1g14.7,1x,1x,
     &            " # (D=1.0) normalization factor")')
     &            rtfac(m)

               end if

*-----------------------------------------------------------------------

               if( itanl(m) .gt. 0 ) then

                  write(iot,'("    angel = ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               end if

*-----------------------------------------------------------------------

               if( itsans(m) .gt. 0 ) then

                  write(iot,'("   sangel = ",i2)') itsans(m)
                  call write_sangel(iot,m,1)

               end if

*-----------------------------------------------------------------------

               if( itaxl(m) .gt. 0 ) then

                  write(iot,'("    x-txt = ",200a1)')
     &            ( itaxt(m)(i:i),i = 1, itaxl(m) )

               end if

*-----------------------------------------------------------------------

               if( itayl(m) .gt. 0 ) then

                  write(iot,'("    y-txt = ",200a1)')
     &            ( itayt(m)(i:i),i = 1, itayl(m) )

               end if

*-----------------------------------------------------------------------

               if( iteps(m) .ne. 0 ) then

                  write(iot,'("   epsout =",i5,10x,
     &            "  # (D=0) generate eps file by ANGEL")')
     &            iteps(m)

               end if

*-----------------------------------------------------------------------
*        check of counter
*-----------------------------------------------------------------------

            do i = 1, 3

               if( itcnt(i,m) .ne. 0 ) then

                  write(iot,'(" ctmin(",i1,") =",i5,10x,
     &            "  # (D=-9999) min. of counter ",i1)')
     &            i, itcnt(i*2+2,m), i

                  write(iot,'(" ctmax(",i1,") =",i5,10x,
     &            "  # (D= 9999) max. of counter ",i1)')
     &            i, itcnt(i*2+3,m), i

                  if( ncntc(i) .eq. 0 )
     &            write(iot,'("# Warning: there is no counter = ",
     &                        i1," definition in [counter].")') i

               end if

            end do

            do i = 1, 3
               if( itcnth(i,m) .ne. 0 ) then
                  write(iot,'(" chmin(",i1,") =",i5,10x,
     &            "  # (D=-9999) min. of history counter ",i1)')
     &            i, itcnth(i*2+2,m), i
                  write(iot,'(" chmax(",i1,") =",i5,10x,
     &            "  # (D= 9999) max. of history counter ",i1)')
     &            i, itcnth(i*2+3,m), i
               end if
            end do

*-----------------------------------------------------------------------

               if( rtstd(m) .gt. 0.d0 ) then

                  write(iot,'("   stdcut =",1p1g14.7,1x,1x,
     &            " # (D=-1) threshold value of STD cutoff")')
     &            rtstd(m)

               end if

*-----------------------------------------------------------------------

         if( itextstat(m) .gt. 0 ) then

            write(iot,'(" iextstat =",i5,10x,
     &      "  # (D=0) option for extended statistical indicators")')
     &           itextstat(m)

         end if

*-----------------------------------------------------------------------

         if( itprodenchk(m) .gt. 0 ) then

            write(iot,'(" prodenmn =",1p1g14.7,1x,
     &     " # (D=1e-2) lower limit of probability density function")')
     &           tprodenmn(m)

            write(iot,'(" prodenmx =",1p1g14.7,1x,
     &     " # (D=1e+2) upper limit of probability density function")')
     &           tprodenmx(m)

            write(iot,'(" nbproden =",i5,10x,
     &   "  # (D=100) number of bins of probability density function")')
     &           nbtproden(m)

         end if

*-----------------------------------------------------------------------
*     anatally
*-----------------------------------------------------------------------
                  yen  = char(92)

         if( itism(m) .gt. 0 ) then

                  write(iot,'(/"anatally start")')

            do i = 1, 10
            if( itkst(i,m) .ge. -1 ) then
                     ch73(1:11) = '   '//schan(i)//' = '
                     ilst = 11
               if( itkst(i,m) .eq. -1 ) then
                     ch73(ilst+1:ilst+4) = ' all'
                     ilst = ilst+4
                     write(iot,'(73(a1))') (ch73(j:j),j=1,ilst)
               else if( itkst(i,m) .ne. 0 ) then
                  do k = 1, itist(i,m)
                        isdnm = idas(itstt(m)+itjst(i,m)+k)
                     if( isdnm .lt. 10 ) then
                        write(ch73(ilst+1:ilst+3),'(i3)') isdnm
                        ilst = ilst + 3
                     else if( isdnm .lt. 100 ) then
                        write(ch73(ilst+1:ilst+4),'(i4)') isdnm
                        ilst = ilst + 4
                     else if( isdnm .lt. 1000 ) then
                        write(ch73(ilst+1:ilst+5),'(i5)') isdnm
                        ilst = ilst + 5
                     else if( isdnm .lt. 10000 ) then
                        write(ch73(ilst+1:ilst+6),'(i6)') isdnm
                        ilst = ilst + 6
                     end if
                  end do

                  if( ilst .le. 70 ) then
                           write(iot,'(73(a1))') (ch73(j:j),j=1,ilst)
                  else if( ilst .le. 130 ) then
                     do l = 65, 73
                        if( ch73(l:l) .eq. ' ' ) goto 111
                     end do
  111                   jlst = l
                        write(iot,'(74(a1))') (ch73(j:j),j=1,jlst), yen
                        write(iot,'(14x,74(a1))')
     &                                        (ch73(j:j),j=jlst+1,ilst)
                  else if( ilst .le. 190 ) then
                     do l = 65, 73
                        if( ch73(l:l) .eq. ' ' ) goto 112
                     end do
  112                   jlst = l
                        write(iot,'(74(a1))') (ch73(j:j),j=1,jlst), yen
                     do l = 125, 133
                        if( ch73(l:l) .eq. ' ' ) goto 113
                     end do
  113                   klst = l
                        write(iot,'(14x,74(a1))')
     &                                   (ch73(j:j),j=jlst+1,klst), yen
                        write(iot,'(14x,74(a1))')
     &                                        (ch73(j:j),j=klst+1,ilst)
                  end if
               end if

            end if
            end do

                  write(iot,'("anatally end")')

         end if


*-----------------------------------------------------------------------
*        sumtally subsection
*-----------------------------------------------------------------------
            if ( icntl==13 .and. nsumtalRead(m)==17 ) then

               call sumtal_echo(iot,m)

            end if
*-----------------------------------------------------------------------

            write(iot,'("#    used :",
     &                  "        main (  %)",
     &                  "        temp (  %)",
     &                  "       total (  %)"/
     &                  "#  memory :",
     &                  3(i12," (",i3,")"))')
     &                  itsmn(m), nint(dble(itsmn(m))/dble(mdas)*100),
     &                  itstm(m), nint(dble(itstm(m))/dble(mdas)*100),
     &                  itsmn(m)+itstm(m),
     &                  nint(dble(itsmn(m)+itstm(m))/dble(mdas)*100)

*-----------------------------------------------------------------------
cfrtati 2024/02/13
            if( itmlp(m).gt.0 ) then
              call displaymtinfo(iot,m)
            end if

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine sedech(iot,m,iaxs,icax)
*                                                                      *
*       input echo of the sed tally                                    *
*       last modified by K.Niita on 2006/02/09                         *
*                                                                      *
************************************************************************

      use sumtallymod
      use sangelmod, only: itsans
      use partmod ! frtati 2021/10/05
      use moddas_mesh
      use moddas_region
      use moddas_tally
*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'

*-----------------------------------------------------------------------

      common /tall00/ itnm, ital(itlmax), itals(itlmax), italm(itlmax)
      common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                itout(itlmax), ittwo(itlmax)

      common /tall02/ itrgn(itlmax), itrgm(itlmax), itreg(itlmax)
      common /tall24/ itrcm(itlmax), itrcg(itlmax)
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

      common /tall16/ itnfr(itlmax), itnfz(itlmax), itnfn(itlmax)

      common /tall14/ ittll(itlmax), ittle(itlmax)
      character ittle*80

      common /tall15/ itanl(itlmax), itang(itlmax)
      character itang*200

      common /tall27/ itaxl(itlmax), itaxt(itlmax),
     &                itayl(itlmax), itayt(itlmax),
     &                itazl(itlmax), itazt(itlmax)
      character itaxt*200, itayt*200, itazt*200

      common /tall19/ itsmn(itlmax), itstm(itlmax)

      common /tall21/ rtfac(itlmax)
      common /tall28/ itgsh(itlmax)
      common /tall29/ itrsh(itlmax)
      common /tall31/ iterl(itlmax)
      common /tall22/ itmcn(itlmax), itmtn(itlmax), itmtt(itlmax)
      common /tall34/ itvm(itlmax)
      common /tall35/ iteps(itlmax)
      common /tall37/ itcnt(9,itlmax)
      common /cntreg/ ncntc(3), ncreg(3), incrc(3), incrt(3), kcont(3)
      common /cntpat/ icpan(3), icpat(3,20,2)

      common /tall38/ itres(itlmax)
      common /tall39/ rtwid(itlmax)

      common /tall46/ itmtr(itlmax,4), rtmtr(itlmax,13)
      common /tall49/ itglt(itlmax)
      common /tall50/ itlmt(itlmax), ite2l(itlmax)   ! CCSE 2022.08.31

      common /tall57/ itsun(itlmax),rtdim(itlmax),rtucv(itlmax),
     &rtrho(itlmax),itmodel(itlmax) ! T.Sato 2022/08/14

      common /tall60/ itger(itlmax)

      common /tall67/ itstd(itlmax), rtstd(itlmax)

      common /tall82/ itcnth(9,itlmax)

      dimension idas(mdas*2) !20220906frtati temporary for fbounds-check
      equivalence ( das, idas )

*-----------------------------------------------------------------------

      common /inpec/ ititl, ipara, ibody, iregn, llarr, itby, itar

*-----------------------------------------------------------------------
       integer   icntl
       integer   inucr

       common /tcntl/  icntl, inucr
!-----------------------------------------------------------------------

      character asfil*100

      character aname(10)*3
      data aname / 'sed','reg','  x','  y','  z','  r',
     &             ' xy',' yz',' xz',' rz'/

      dimension ian(6)
      dimension jmat(20)

*-----------------------------------------------------------------------
      common /tall79/ itallech
      if( itallech .eq. 0 ) return

*-----------------------------------------------------------------------
*     input echo
*-----------------------------------------------------------------------

            write(iot,'("[ T-SED ]")')

                  write(iot,'("    title = ",80a1)')
     &            ( ittle(m)(i:i), i = 1, ittll(m) )

            if( itmsh(m) .eq. 1 ) then

                  write(iot,'("     mesh =  reg",11x,
     &            " # mesh type is region-wise")')

                     idas1 = mmmax
                     idas2 = ( idas1 + itrgn(m) - 1 ) * 2 + 1
                     idas3 = idas2 + itrgn(m)

                  call echrg(1,iterl(m),0,1,iot,itrcm(m),
     &                       idas_itrcg(itrcg(m)),
     &                       itrgn(m),itrgm(m),idas_itreg(itreg(m)),
     &                       das(idas1),idas(idas2),
     &                       itrnv(m),idas(itriv(m)),das(itrrv(m)),
     &                       idas3)

            else if( itmsh(m) .eq. 2 ) then

                  write(iot,'("     mesh =  r-z",11x,
     &            " # mesh type is r-z scoring mesh")')

                  write(iot,'("       x0 = ",1p1g14.7,1x,
     &            " # center x-position of r-z mesh")') rtrx0(m)

                  write(iot,'("       y0 = ",1p1g14.7,1x,
     &            " # center y-position of r-z mesh")') rtry0(m)

                  call echmty(0,iot,'r',itrty(m),itrnm(m),
     &                        rtrdl(m),rtrmi(m),rtrma(m),itrrg(m)
     &                        ,abs(itrnm(m)+1),das_itrrg)
                  call echmty(0,iot,'z',itzty(m),itznm(m),
     &                        rtzdl(m),rtzmi(m),rtzma(m),itzrg(m)
     &                        ,abs(itznm(m)+1),das_itzrg)

            else if( itmsh(m) .eq. 3 ) then

                  write(iot,'("     mesh =  xyz",11x,
     &            " # mesh type is xyz scoring mesh")')

                  call echmty(0,iot,'x',itxty(m),itxnm(m),
     &                        rtxdl(m),rtxmi(m),rtxma(m),itxrg(m)
     &                        ,abs(itxnm(m)+1),das_itxrg)
                  call echmty(0,iot,'y',ityty(m),itynm(m),
     &                        rtydl(m),rtymi(m),rtyma(m),ityrg(m)
     &                        ,abs(itynm(m)+1),das_ityrg)
                  call echmty(0,iot,'z',itzty(m),itznm(m),
     &                        rtzdl(m),rtzmi(m),rtzma(m),itzrg(m)
     &                        ,abs(itznm(m)+1),das_itzrg)

            end if

*-----------------------------------------------------------------------

            if( itmtr(m,1) .gt. 0 ) then

               if( itmtr(m,1) .eq. 1 ) then

                     write(iot,'("     trcl = ",i6,9x,
     &               " # transform ID")') itmtr(m,3)

               else

                  if( itmtr(m,2) .eq. 0 ) then

                     write(iot,'("     trcl = ",
     &                           3(1p1e15.7),
     &                           3(/12x,3(1p1e15.7)),i5)')
     &               ( rtmtr(m,j), j = 1, 12 ),
     &               nint( rtmtr(m,13) )

                  else

                     write(iot,'("    *trcl = ",
     &                           3(1p1e15.7),
     &                           3(/12x,3(1p1e15.7)),i5)')
     &               ( rtmtr(m,j), j = 1, 12 ),
     &               nint( rtmtr(m,13) )

                  end if

               end if

            end if

*-----------------------------------------------------------------------

                  call echmty2(0,iot,'e','se',itety(m),itenm(m),
     &                         rtedl(m),rtemi(m),rtema(m),iterg(m)
     &                        ,abs(itenm(m)+1),das_iterg)

*-----------------------------------------------------------------------

               if( itsun(m) .eq. 0 ) then

                  write(iot,'("  se-unit = ",i4,11x,
     &            " # se-unit is [event]")') itsun(m)

               elseif( itsun(m) .eq. 1 ) then

                  write(iot,'("  se-unit = ",i4,11x,
     &            " # se-unit is [MeV]")') itsun(m)

               else if( itsun(m) .eq. 2 ) then

                  write(iot,'("  se-unit = ",i4,11x,
     &            " # se-unit is [keV/um]")') itsun(m)

               else if( itsun(m) .eq. 3 ) then

                  write(iot,'("  se-unit = ",i4,11x,
     &            " # se-unit is [Gy]")') itsun(m)

               end if

*-----------------------------------------------------------------------

               if( itunt(m) .eq. 1 ) then

                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [cm/z/source)]")') itunt(m)

               else if( itunt(m) .eq. 2 ) then

                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [MeV/z/source)]")') itunt(m)

               else if( itunt(m) .eq. 3 ) then

                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [cm/ln(z)/source)]")') itunt(m)

               else if( itunt(m) .eq. 4 ) then

                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [MeV/ln(z)/source)]")') itunt(m)

               else if( itunt(m) .eq. 5 ) then

                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [cm/source]")') itunt(m)

               else if( itunt(m) .eq. 6 ) then

                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [MeV/source]")') itunt(m)

               else if( itunt(m) .eq. 7 ) then  ! T.Sato 2017/11/21

                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is z*f(z) [dimensionless]")') itunt(m)

               else if( itunt(m) .eq. 8 ) then ! T.Sato 2017/11/21

                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is z*d(z) [Gy]")') itunt(m)

               end if

*-----------------------------------------------------------------------

                  asfil = '  # file name of output for the above axis'

            if( icax .eq. 0 ) then

                  icn = itaxn(m)

               do i = 1, icn

                  ian(i) = i

               end  do

            else

                  icn = 1
                  ian(1) = iaxs

            end if

*-----------------------------------------------------------------------

                  write(iot,'("    model =",i5,11x,
     &            " # (D=0) Model type 0:old, 1:new")')
     &            itmodel(m)

*-----------------------------------------------------------------------

                  write(iot,'("    cdiam =",1p1g14.7,1x,1x,
     &            " # (D=1.0) diameter of target sphere [um]")')
     &            rtdim(m)

*-----------------------------------------------------------------------

                  write(iot,'("   letmat =",i5,11x,
     &            " # (D=0) mat ID for LET, 0:real mat, ",
     &                                    "<0: electron for H2O")')
     &            itlmt(m)

*-----------------------------------------------------------------------

                if( rtrho(m) .ne. 1.0d0) then
                  write(iot,'("   rhomat =",1p1g14.7,1x,1x,
     &            " # (D=1.0) density of letmat [g/cm^3]")')
     &            rtrho(m)
                end if

*-----------------------------------------------------------------------

               if( itmtn(m) .eq. 0 ) then

                  write(iot,'(" material =  all"11x,
     &            " # (D=all) number of specific material")')

               else if( itmtn(m) .gt. 0 ) then

                  write(iot,'(" material =",i5,11x,
     &            " # (D=all) number of specific material")')
     &            itmtn(m) * itmcn(m)

                     k = 0

                  do j = 1, itmtn(m)

                     k = k + 1

                     jmat(k) = ismte( itmtt(m) + j - 1 )

                     if( k / 10 * 10 .eq. k .or. j .eq. itmtn(m) ) then

                        write(iot,'(11x,10(i6))') ( jmat(i),i = 1, k )

                        k = 0

                     end if

                  end do

               end if

*-----------------------------------------------------------------------

               idtyp = 0

            do i = 1, icn

               if( itaxs(m,ian(i)) .ge. 7 .and.
     &             itaxs(m,ian(i)) .le. 10 )  idtyp = ittwo(m)

            end do

            if( idtyp .ne. 0 ) then

               write(iot,'("  2D-type = ",i4,11x,
     &         " # 1:Cont, 2:Clust, 3:Color, 4:xyz, 5:mat,",
     &         " 6:Clust+Cont, 7:Col+Cont")') idtyp

            end if

*-----------------------------------------------------------------------

            do i = 1, icn

                  iax = ian(i)

                  write(iot,'("     axis =  ",a3,11x,
     &            " # axis of output")')
     &            aname( itaxs(m,iax) )

             if ( itfll(m,iax) .gt. 0 ) then

                  msfile = max( 14, itfll(m,iax) )

                  write(iot,'("     file = ",100a1)')
     &                  ( ctfln(m,iax)(j:j), j = 1, msfile ),
     &                  ( asfil(j:j), j = 1, 42 )

            end if

            end do

               call echprt(1,iot,m,itpan,itpat,jtpat,itnm,itmxpt,itmxgp)

*-----------------------------------------------------------------------

               if( rtfac(m) .ne. 1.0d0 ) then

                  write(iot,'("   factor =",1p1g14.7,1x,1x,
     &            " # (D=1.0) normalization factor")')
     &            rtfac(m)

               end if

*-----------------------------------------------------------------------

               if( itmtn(m) .gt. 0 .and. itmsh(m) .eq. 3 ) then

                  write(iot,'("   volmat = ",i4,11x,
     &            " # (D=9) vol.corr., 0:no, number of scan / dir.")')
     &            itvm(m)

               end if

*-----------------------------------------------------------------------

               if( itanl(m) .gt. 0 ) then

                  write(iot,'("    angel = ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               end if

*-----------------------------------------------------------------------

               if( itsans(m) .gt. 0 ) then

                  write(iot,'("   sangel = ",i2)') itsans(m)
                  call write_sangel(iot,m,1)

               end if

*-----------------------------------------------------------------------

               if( itaxl(m) .gt. 0 ) then

                  write(iot,'("    x-txt = ",200a1)')
     &            ( itaxt(m)(i:i),i = 1, itaxl(m) )

               end if

*-----------------------------------------------------------------------

               if( itayl(m) .gt. 0 ) then

                  write(iot,'("    y-txt = ",200a1)')
     &            ( itayt(m)(i:i),i = 1, itayl(m) )

               end if

*-----------------------------------------------------------------------

               if( itazl(m) .gt. 0 ) then

                  write(iot,'("    z-txt = ",200a1)')
     &            ( itazt(m)(i:i),i = 1, itazl(m) )

               end if

*-----------------------------------------------------------------------

               if( itgsh(m) .gt. 0 ) then

                  write(iot,'("    gshow = ",i4,10x,
     &            "  # 0: no 1:bnd, 2:bnd+mat, ",
     &     "3:bnd+reg 4:bnd+lat 5:bmp style")')
     &            itgsh(m)

               end if

*-----------------------------------------------------------------------

               if( itglt(m) .ge. 0 ) then

                  write(iot,'("    gslat = ",i4,10x,
     &            "  # 0: no line for lat, 1: with,",
     &               " 2: no line between the same cell")')
     &            itglt(m)

               end if

*-----------------------------------------------------------------------

               if( itgsh(m) .gt. 0 .and. itger(m) .gt. 0 ) then

                  write(iot,'("    ginfo = ",i4,11x,
     &            " # (D=0) 0:no, 1:show in graph, 2:add file")')
     &            itger(m)

               end if

*-----------------------------------------------------------------------

               if( itrsh(m) .gt. 0 ) then

                  write(iot,'("    rshow = ",i4,10x,
     &            "  # 0: no 1:bnd, 2:bnd+mat, ",
     &     "3:bnd+reg 4:bnd+lat 5:bmp style")')
     &            itrsh(m)

                  call echmty(1,iot,'x',itxty(m),itxnm(m),
     &                        rtxdl(m),rtxmi(m),rtxma(m),itxrg(m)
     &                        ,abs(itxnm(m)+1),das_itxrg)
                  call echmty(1,iot,'y',ityty(m),itynm(m),
     &                        rtydl(m),rtymi(m),rtyma(m),ityrg(m)
     &                        ,abs(itynm(m)+1),das_ityrg)
                  call echmty(1,iot,'z',itzty(m),itznm(m),
     &                        rtzdl(m),rtzmi(m),rtzma(m),itzrg(m)
     &                        ,abs(itznm(m)+1),das_itzrg)

               end if

*-----------------------------------------------------------------------

               if( itgsh(m) .gt. 0 .or. itrsh(m) .gt. 0 ) then

                  write(iot,'("    resol = ",i4,10x,
     &            "  # (D=1) resolution of gshow or rshow")')
     &            itres(m)

                  write(iot,'("    width =",1p1g14.7,1x,1x,
     &            " # (D=0.5) width of lines for gshow or rshow")')
     &            rtwid(m)

               end if

*-----------------------------------------------------------------------

               if( iterl(m) .ne. 72 ) then

                  write(iot,'("   iechrl =",i5,10x,
     &            "  # (D=72) length reg of val or vol")')
     &            iterl(m)

               end if

*-----------------------------------------------------------------------

               if( iteps(m) .ne. 0 ) then

                  write(iot,'("   epsout =",i5,10x,
     &            "  # (D=0) generate eps file by ANGEL")')
     &            iteps(m)

               end if

*-----------------------------------------------------------------------
*        check of counter
*-----------------------------------------------------------------------

            do i = 1, 3

               if( itcnt(i,m) .ne. 0 ) then

                  write(iot,'(" ctmin(",i1,") =",i5,10x,
     &            "  # (D=-9999) min. of counter ",i1)')
     &            i, itcnt(i*2+2,m), i

                  write(iot,'(" ctmax(",i1,") =",i5,10x,
     &            "  # (D= 9999) max. of counter ",i1)')
     &            i, itcnt(i*2+3,m), i

                  if( ncntc(i) .eq. 0 )
     &            write(iot,'("# Warning: there is no counter = ",
     &                        i1," definition in [counter].")') i

               end if

            end do

            do i = 1, 3
               if( itcnth(i,m) .ne. 0 ) then
                  write(iot,'(" chmin(",i1,") =",i5,10x,
     &            "  # (D=-9999) min. of history counter ",i1)')
     &            i, itcnth(i*2+2,m), i
                  write(iot,'(" chmax(",i1,") =",i5,10x,
     &            "  # (D= 9999) max. of history counter ",i1)')
     &            i, itcnth(i*2+3,m), i
               end if
            end do

*-----------------------------------------------------------------------

               if( rtstd(m) .gt. 0.d0 ) then

                  write(iot,'("   stdcut =",1p1g14.7,1x,1x,
     &            " # (D=-1) threshold value of STD cutoff")')
     &            rtstd(m)

               end if

*-----------------------------------------------------------------------
               if( itmxang(m) .ne. 0 ) then

                  write(iot,'(" maxangel =",i5,10x,
     &            "  # (D=0) number of output particles")')
     &            itmxang(m)

               end if

*-----------------------------------------------------------------------
*        sumtally subsection
*-----------------------------------------------------------------------
            if ( icntl==13 .and. nsumtalRead(m)==15 ) then

               call sumtal_echo(iot,m)

            end if
*-----------------------------------------------------------------------

            write(iot,'("#    used :",
     &                  "        main (  %)",
     &                  "        temp (  %)",
     &                  "       total (  %)"/
     &                  "#  memory :",
     &                  3(i12," (",i3,")"))')
     &                  itsmn(m), nint(dble(itsmn(m))/dble(mdas)*100),
     &                  itstm(m), nint(dble(itstm(m))/dble(mdas)*100),
     &                  itsmn(m)+itstm(m),
     &                  nint(dble(itsmn(m)+itstm(m))/dble(mdas)*100)

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine dps2tech(iot,m,iaxs,icax)
*                                                                      *
*       input echo of the deposit2 tally                               *
*       last modified by K.Niita on 2005/12/12                         *
*                                                                      *
************************************************************************

      use sumtallymod
      use sangelmod, only: itsans
      use partmod ! frtati 2021/10/05
      use moddas_mesh
      use moddas_region
*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'

*-----------------------------------------------------------------------

      common /tall00/ itnm, ital(itlmax), itals(itlmax), italm(itlmax)
      common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                itout(itlmax), ittwo(itlmax)

      common /tall02/ itrgn(itlmax), itrgm(itlmax), itreg(itlmax)
      common /tall24/ itrcm(itlmax), itrcg(itlmax)
      common /tall25/ itrnv(itlmax), itriv(itlmax), itrrv(itlmax)

      common /tall07/ itety(itlmax), itenm(itlmax), iterg(itlmax),
     &                rtemi(itlmax), rtema(itlmax), rtedl(itlmax)


      common /tall10/ itaxn(itlmax), itaxs(itlmax,6)
      common /tall11/ itfln(itlmax), itfll(itlmax,6), ctfln(itlmax,6)
      character ctfln*100

      common /tall14/ ittll(itlmax), ittle(itlmax)
      character ittle*80

      common /tall15/ itanl(itlmax), itang(itlmax)
      character itang*200

      common /tall27/ itaxl(itlmax), itaxt(itlmax),
     &                itayl(itlmax), itayt(itlmax),
     &                itazl(itlmax), itazt(itlmax)
      character itaxt*200, itayt*200, itazt*200

      common /tall19/ itsmn(itlmax), itstm(itlmax)
      common /tall20/ ittty(itlmax), ittnm(itlmax), ittrg(itlmax),
     &                rttmi(itlmax), rttma(itlmax), rttdl(itlmax)
      common /tall21/ rtfac(itlmax)

      common /tall31/ iterl(itlmax)
      common /tall35/ iteps(itlmax)
      common /tall37/ itcnt(9,itlmax)
      common /cntreg/ ncntc(3), ncreg(3), incrc(3), incrt(3), kcont(3)
      common /cntpat/ icpan(3), icpat(3,20,2)

      common /tall50/ itlmt(itlmax), ite2l(itlmax)   ! CCSE 2022.08.31
      common /tall53/ itdfn(itlmax,2)

      common /tall54/ itdfn2(itlmax,2)
      common /tall55/ itlmt2(itlmax)

      common /tall56/ itety2(itlmax), itenm2(itlmax), iterg2(itlmax),
     &                rtemi2(itlmax), rtema2(itlmax), rtedl2(itlmax)

      common /tall60/ itger(itlmax)

      common /tall67/ itstd(itlmax), rtstd(itlmax)

      common /tall82/ itcnth(9,itlmax)

      dimension idas(mdas*2) !20220906frtati temporary for fbounds-check
      equivalence ( das, idas )

*-----------------------------------------------------------------------

      common /inpec/ ititl, ipara, ibody, iregn, llarr, itby, itar

*-----------------------------------------------------------------------
       integer   icntl
       integer   inucr

       common /tcntl/  icntl, inucr
!-----------------------------------------------------------------------

      character asfil*100

      character aname(9)*4
      data aname / 'eng1','eng2',' e12',' e21','t-e1',
     &             'e1-t','t-e2','e2-t','   t'/

      dimension ian(6)
      dimension jmat(20)

*-----------------------------------------------------------------------
      common /tall79/ itallech
      if( itallech .eq. 0 ) return

*-----------------------------------------------------------------------
*     input echo
*-----------------------------------------------------------------------

            write(iot,'("[ T-Deposit2 ]")')

                  write(iot,'("    title = ",80a1)')
     &            ( ittle(m)(i:i), i = 1, ittll(m) )

            if( itmsh(m) .eq. 1 ) then

                  write(iot,'("     mesh =  reg",11x,
     &            " # mesh type is only region-wise, two region")')

                     idas1 = mmmax
                     idas2 = ( idas1 + itrgn(m) - 1 ) * 2 + 1
                     idas3 = idas2 + itrgn(m)

                  call echrg(1,iterl(m),0,1,iot,itrcm(m),
     &                       idas_itrcg(itrcg(m)),
     &                       itrgn(m),itrgm(m),idas_itreg(itreg(m)),
     &                       das(idas1),idas(idas2),
     &                       itrnv(m),idas(itriv(m)),das(itrrv(m)),
     &                       idas3)

            end if

*-----------------------------------------------------------------------

               if( itunt(m) .eq. 1 ) then

                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is Number [1/source]")') itunt(m)

               else if( itunt(m) .eq. 2 ) then

                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is Number [1/nsec/source]")') itunt(m)

               end if

*-----------------------------------------------------------------------

                  asfil = '  # file name of output for the above axis'

            if( icax .eq. 0 ) then

                  icn = itaxn(m)

               do i = 1, icn

                  ian(i) = i

               end  do

            else

                  icn = 1
                  ian(1) = iaxs

            end if

*-----------------------------------------------------------------------

                  write(iot,'("  letmat1 =",i5,11x,
     &            " # (D=0) mat ID for LET, 0:real mat, ",
     &                                    "<0: electron for H2O")')
     &            itlmt(m)

                  write(iot,'("  letmat2 =",i5,11x,
     &            " # (D=0) mat ID for LET, 0:real mat, ",
     &                                    "<0: electron for H2O")')
     &            itlmt2(m)

*-----------------------------------------------------------------------

                  write(iot,'(" dedxfnc1 =",i5,11x,
     &            " # (D=0) user defined multiplier, 0(no), 1, 2")')
     &            itdfn(m,1)

                  write(iot,'(" dedxfnc2 =",i5,11x,
     &            " # (D=0) user defined multiplier, 0(no), 1, 2")')
     &            itdfn2(m,1)

*-----------------------------------------------------------------------

                  call echmty2(0,iot,'e','e1',itety(m),itenm(m),
     &                         rtedl(m),rtemi(m),rtema(m),iterg(m)
     &                        ,abs(itenm(m)+1),das_iterg)
                  call echmty2(0,iot,'e','e2',itety2(m),itenm2(m),
     &                         rtedl2(m),rtemi2(m),rtema2(m),iterg2(m)
     &                        ,abs(itenm2(m)+1),das_iterg2)

*-----------------------------------------------------------------------

               if( ittty(m) .gt. 0 ) then

                  call echmty(0,iot,'t',ittty(m),ittnm(m),
     &                        rttdl(m),rttmi(m),rttma(m),ittrg(m)
     &                        ,abs(ittnm(m)+1),das_ittrg)

               end if

*-----------------------------------------------------------------------

               idtyp = 0

            do i = 1, icn

               if( itaxs(m,ian(i)) .ge. 3 .and.
     &             itaxs(m,ian(i)) .le. 8 )  idtyp = ittwo(m)

            end do

            if( idtyp .ne. 0 ) then

               write(iot,'("  2D-type = ",i4,11x,
     &         " # 1:Cont, 2:Clust, 3:Color, 4:xyz, 5:mat,",
     &         " 6:Clust+Cont, 7:Col+Cont")') idtyp

            end if

*-----------------------------------------------------------------------

            do i = 1, icn

                  iax = ian(i)

                  write(iot,'("     axis =  ",a4,10x,
     &            " # axis of output")')
     &            aname( itaxs(m,iax) )

             if ( itfll(m,iax) .gt. 0 ) then

                  msfile = max( 14, itfll(m,iax) )

                  write(iot,'("     file = ",100a1)')
     &                  ( ctfln(m,iax)(j:j), j = 1, msfile ),
     &                  ( asfil(j:j), j = 1, 42 )

             end if

            end do

*-----------------------------------------------------------------------
               call echprt(1,iot,m,itpan,itpat,jtpat,itnm,itmxpt,itmxgp)

*-----------------------------------------------------------------------

               if( rtfac(m) .ne. 1.0d0 ) then

                  write(iot,'("   factor =",1p1g14.7,1x,1x,
     &            " # (D=1.0) normalization factor")')
     &            rtfac(m)

               end if

*-----------------------------------------------------------------------

               if( itanl(m) .gt. 0 ) then

                  write(iot,'("    angel = ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               end if

*-----------------------------------------------------------------------

               if( itsans(m) .gt. 0 ) then

                  write(iot,'("   sangel = ",i2)') itsans(m)
                  call write_sangel(iot,m,1)

               end if

*-----------------------------------------------------------------------

               if( itaxl(m) .gt. 0 ) then

                  write(iot,'("    x-txt = ",200a1)')
     &            ( itaxt(m)(i:i),i = 1, itaxl(m) )

               end if

*-----------------------------------------------------------------------

               if( itayl(m) .gt. 0 ) then

                  write(iot,'("    y-txt = ",200a1)')
     &            ( itayt(m)(i:i),i = 1, itayl(m) )

               end if

*-----------------------------------------------------------------------

               if( itazl(m) .gt. 0 ) then

                  write(iot,'("    z-txt = ",200a1)')
     &            ( itazt(m)(i:i),i = 1, itazl(m) )

               end if

*-----------------------------------------------------------------------

               if( iterl(m) .ne. 72 ) then

                  write(iot,'("   iechrl =",i5,10x,
     &            "  # (D=72) length reg of val or vol")')
     &            iterl(m)

               end if

*-----------------------------------------------------------------------

               if( iteps(m) .ne. 0 ) then

                  write(iot,'("   epsout =",i5,10x,
     &            "  # (D=0) generate eps file by ANGEL")')
     &            iteps(m)

               end if

*-----------------------------------------------------------------------
*        check of counter
*-----------------------------------------------------------------------

            do i = 1, 3

               if( itcnt(i,m) .ne. 0 ) then

                  write(iot,'(" ctmin(",i1,") =",i5,10x,
     &            "  # (D=-9999) min. of counter ",i1)')
     &            i, itcnt(i*2+2,m), i

                  write(iot,'(" ctmax(",i1,") =",i5,10x,
     &            "  # (D= 9999) max. of counter ",i1)')
     &            i, itcnt(i*2+3,m), i

                  if( ncntc(i) .eq. 0 )
     &            write(iot,'("# Warning: there is no counter = ",
     &                        i1," definition in [counter].")') i

               end if

            end do

            do i = 1, 3
               if( itcnth(i,m) .ne. 0 ) then
                  write(iot,'(" chmin(",i1,") =",i5,10x,
     &            "  # (D=-9999) min. of history counter ",i1)')
     &            i, itcnth(i*2+2,m), i
                  write(iot,'(" chmax(",i1,") =",i5,10x,
     &            "  # (D= 9999) max. of history counter ",i1)')
     &            i, itcnth(i*2+3,m), i
               end if
            end do

*-----------------------------------------------------------------------

               if( rtstd(m) .gt. 0.d0 ) then

                  write(iot,'("   stdcut =",1p1g14.7,1x,1x,
     &            " # (D=-1) threshold value of STD cutoff")')
     &            rtstd(m)

               end if

*-----------------------------------------------------------------------
               if( itmxang(m) .ne. 0 ) then

                  write(iot,'(" maxangel =",i5,10x,
     &            "  # (D=0) number of output particles")')
     &            itmxang(m)

               end if

*-----------------------------------------------------------------------

            write(iot,'("#    used :",
     &                  "        main (  %)",
     &                  "        temp (  %)",
     &                  "       total (  %)"/
     &                  "#  memory :",
     &                  3(i12," (",i3,")"))')
     &                  itsmn(m), nint(dble(itsmn(m))/dble(mdas)*100),
     &                  itstm(m), nint(dble(itstm(m))/dble(mdas)*100),
     &                  itsmn(m)+itstm(m),
     &                  nint(dble(itsmn(m)+itstm(m))/dble(mdas)*100)

*-----------------------------------------------------------------------
*        sumtally subsection
*-----------------------------------------------------------------------
            if ( icntl==13 .and. nsumtalRead(m)==14 ) then

               call sumtal_echo(iot,m)

            end if
*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine depstech(iot,m,iaxs,icax)
*                                                                      *
*       input echo of the deposit tally                                *
*       last modified by K.Niita on 2005/12/16                         *
*                                                                      *
************************************************************************

      use sumtallymod
      use sangelmod, only: itsans
      use partmod ! frtati 2021/10/05
      use moddas_mesh
      use moddas_region
      use moddas_tally

*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'

*-----------------------------------------------------------------------

      common /tall00/ itnm, ital(itlmax), itals(itlmax), italm(itlmax)
      common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                itout(itlmax), ittwo(itlmax)

      common /tall02/ itrgn(itlmax), itrgm(itlmax), itreg(itlmax)
      common /tall24/ itrcm(itlmax), itrcg(itlmax)
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

      common /tall16/ itnfr(itlmax), itnfz(itlmax), itnfn(itlmax)

      common /tall14/ ittll(itlmax), ittle(itlmax)
      character ittle*80

      common /tall15/ itanl(itlmax), itang(itlmax)
      character itang*200

      common /tall27/ itaxl(itlmax), itaxt(itlmax),
     &                itayl(itlmax), itayt(itlmax),
     &                itazl(itlmax), itazt(itlmax)
      character itaxt*200, itayt*200, itazt*200

      common /tall19/ itsmn(itlmax), itstm(itlmax)
      common /tall20/ ittty(itlmax), ittnm(itlmax), ittrg(itlmax),
     &                rttmi(itlmax), rttma(itlmax), rttdl(itlmax)

      common /tall21/ rtfac(itlmax)
      common /tall28/ itgsh(itlmax)
      common /tall29/ itrsh(itlmax)
      common /tall31/ iterl(itlmax)
      common /tall22/ itmcn(itlmax), itmtn(itlmax), itmtt(itlmax)
      common /tall34/ itvm(itlmax)
      common /tall35/ iteps(itlmax)
      common /tall36/ itdpo(itlmax)   ! S.Abe 2015/12/03
      common /tall37/ itcnt(9,itlmax)
      common /cntreg/ ncntc(3), ncreg(3), incrc(3), incrt(3), kcont(3)
      common /cntpat/ icpan(3), icpat(3,20,2)

      common /tall38/ itres(itlmax)
      common /tall39/ rtwid(itlmax)

      common /tall46/ itmtr(itlmax,4), rtmtr(itlmax,13)
      common /tall49/ itglt(itlmax)
      common /tall50/ itlmt(itlmax), ite2l(itlmax)   ! CCSE 2022.08.31
      common /tall53/ itdfn(itlmax,2)

      common /tall60/ itger(itlmax)

c T.Sato 2014/8/19 for detector resolution
      common /tall61/ rtdre(itlmax),rtdfa(itlmax)
      common /tall63/ itbmp(itlmax)
      common /tall64/ itvtk(itlmax),itvtkfmt(itlmax)

      common /tall70/ itrwgtsum(itlmax), itrgn1(itlmax),
     &                itrgm1(itlmax), itrncd(itlmax)

      common /tall67/ itstd(itlmax), rtstd(itlmax)

      common /tall75/ itnlatcel(itlmax), itnlatmem(itlmax)

      integer :: itfoam
      common /tall76/ itfoam(itlmax)

      common /tall82/ itcnth(9,itlmax)

      dimension idas(mdas*2) !20220906frtati temporary for fbounds-check
      equivalence ( das, idas )

*-----------------------------------------------------------------------

      common /inpec/ ititl, ipara, ibody, iregn, llarr, itby, itar

*-----------------------------------------------------------------------
       integer   icntl
       integer   inucr

       common /tcntl/  icntl, inucr
!-----------------------------------------------------------------------

      character asfil*100

      character aname(14)*5
      data aname / '  eng','  reg','    x','    y','    z','    r',
     &             '   xy','   yz','   xz','   rz','    t','t-eng',
     &             'eng-t','  tet'/

      dimension ian(6)
      dimension jmat(20)

! T.Sato 2025/03/05 for mother parameter -------------------------------

      common /tall12/ itman(itlmax), itmat(itlmax), itmct(itlmax),
     &                itnun(itlmax), itnuc(itlmax),
     &                itndz(itlmax), itndn(itlmax),
     &                itnkz(itlmax), itnkn(itlmax)

      character chau*8
      character chaus(7)*8
      character elmnt(104)*3

      data elmnt /
     &    'H  ','He ','Li ','Be ','B  ','C  ','N  ','O  ','F  ','Ne ',
     &    'Na ','Mg ','Al ','Si ','P  ','S  ','Cl ','Ar ','K  ','Ca ',
     &    'Sc ','Ti ','V  ','Cr ','Mn ','Fe ','Co ','Ni ','Cu ','Zn ',
     &    'Ga ','Ge ','As ','Se ','Br ','Kr ','Rb ','Sr ','Y  ','Zr ',
     &    'Nb ','Mo ','Tc ','Ru ','Rh ','Pd ','Ag ','Cd ','In ','Sn ',
     &    'Sb ','Te ','I  ','Xe ','Cs ','Ba ','La ','Ce ','Pr ','Nd ',
     &    'Pm ','Sm ','Eu ','Gd ','Tb ','Dy ','Ho ','Er ','Tm ','Yb ',
     &    'Lu ','Hf ','Ta ','W  ','Re ','Os ','Ir ','Pt ','Au ','Hg ',
     &    'Tl ','Pb ','Bi ','Po ','At ','Rn ','Fr ','Ra ','Ac ','Th ',
     &    'Pa ','U  ','Np ','Pu ','Am ','Cm ','Bk ','Cf ','Es ','Fm ',
     &    'Md ','No ','Lr ','Ku ' /


*-----------------------------------------------------------------------
      common /tall79/ itallech
      if( itallech .eq. 0 ) return

*-----------------------------------------------------------------------
*     input echo
*-----------------------------------------------------------------------

            write(iot,'("[ T-Deposit ]")')

                  write(iot,'("    title = ",80a1)')
     &            ( ittle(m)(i:i), i = 1, ittll(m) )

            if( itmsh(m) .eq. 1 ) then

                  write(iot,'("     mesh =  reg",11x,
     &            " # mesh type is region-wise")')

                     idas1 = mmmax
                     idas2 = ( idas1 + itrgn(m) - 1 ) * 2 + 1
                     idas3 = idas2 + itrgn(m)

               if( itrwgtsum(m) .eq. 1 ) then

                  call echrg_wgtsum(m,1,
     &                       iterl(m),0,1,iot,itrcm(m),
     &                       idas_itrcg(itrcg(m)),
     &                       itrgn(m),itrgm(m),idas_itreg(itreg(m)),
     &                       das(idas1),idas(idas2),
     &                       itrnv(m),idas(itriv(m)),das(itrrv(m)),
     &                       idas3)

               else

                  call echrg(1,iterl(m),0,1,iot,itrcm(m),
     &                       idas_itrcg(itrcg(m)),
     &                       itrgn(m),itrgm(m),idas_itreg(itreg(m)),
     &                       das(idas1),idas(idas2),
     &                       itrnv(m),idas(itriv(m)),das(itrrv(m)),
     &                       idas3)

               endif

            else if( itmsh(m) .eq. 2 ) then

                  write(iot,'("     mesh =  r-z",11x,
     &            " # mesh type is r-z scoring mesh")')

                  write(iot,'("       x0 = ",1p1g14.7,1x,
     &            " # center x-position of r-z mesh")') rtrx0(m)

                  write(iot,'("       y0 = ",1p1g14.7,1x,
     &            " # center y-position of r-z mesh")') rtry0(m)

                  call echmty(0,iot,'r',itrty(m),itrnm(m),
     &                        rtrdl(m),rtrmi(m),rtrma(m),itrrg(m)
     &                        ,abs(itrnm(m)+1),das_itrrg)
                  call echmty(0,iot,'z',itzty(m),itznm(m),
     &                        rtzdl(m),rtzmi(m),rtzma(m),itzrg(m)
     &                        ,abs(itznm(m)+1),das_itzrg)

            else if( itmsh(m) .eq. 3 ) then

                  write(iot,'("     mesh =  xyz",11x,
     &            " # mesh type is xyz scoring mesh")')

                  call echmty(0,iot,'x',itxty(m),itxnm(m),
     &                        rtxdl(m),rtxmi(m),rtxma(m),itxrg(m)
     &                        ,abs(itxnm(m)+1),das_itxrg)
                  call echmty(0,iot,'y',ityty(m),itynm(m),
     &                        rtydl(m),rtymi(m),rtyma(m),ityrg(m)
     &                        ,abs(itynm(m)+1),das_ityrg)
                  call echmty(0,iot,'z',itzty(m),itznm(m),
     &                        rtzdl(m),rtzmi(m),rtzma(m),itzrg(m)
     &                        ,abs(itznm(m)+1),das_itzrg)

            else if( itmsh(m) .eq. 4 ) then

                  write(iot,'("     mesh =  tet",11x,
     &            " # mesh type is tetra geoemtry")')

                  call echtet(iot,itrgm(m),idas_itreg(itreg(m)))

            end if

*-----------------------------------------------------------------------

            if( itmtr(m,1) .gt. 0 ) then

               if( itmtr(m,1) .eq. 1 ) then

                     write(iot,'("     trcl = ",i6,9x,
     &               " # transform ID")') itmtr(m,3)

               else

                  if( itmtr(m,2) .eq. 0 ) then

                     write(iot,'("     trcl = ",
     &                           3(1p1e15.7),
     &                           3(/12x,3(1p1e15.7)),i5)')
     &               ( rtmtr(m,j), j = 1, 12 ),
     &               nint( rtmtr(m,13) )

                  else

                     write(iot,'("    *trcl = ",
     &                           3(1p1e15.7),
     &                           3(/12x,3(1p1e15.7)),i5)')
     &               ( rtmtr(m,j), j = 1, 12 ),
     &               nint( rtmtr(m,13) )

                  end if

               end if

            end if

*-----------------------------------------------------------------------
               if( itunt(m) .eq. 0 ) then

                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [Gy/source]")') itunt(m)

               else if( itunt(m) .eq. 1 ) then

                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [MeV/cm^3/source]")') itunt(m)

               else if( itunt(m) .eq. 2 ) then

                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [MeV/source)]")') itunt(m)

               else if( itunt(m) .eq. 3 ) then

                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [1/source]",
     &            " : only for output=deposit")') itunt(m)

               else if( itunt(m) .eq. 4 ) then

                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [1/nsec/source]",
     &            " : only for output=deposit")') itunt(m)

               end if

*-----------------------------------------------------------------------

                  asfil = '  # file name of output for the above axis'

            if( icax .eq. 0 ) then

                  icn = itaxn(m)

               do i = 1, icn

                  ian(i) = i

               end  do

            else

                  icn = 1
                  ian(1) = iaxs

            end if

*-----------------------------------------------------------------------

                  write(iot,'("   letmat =",i5,11x,
     &            " # (D=0) mat ID for LET, 0:real mat, ",
     &                                    "<0: electron for H2O")')
     &            itlmt(m)

*-----------------------------------------------------------------------

                  write(iot,'("  dedxfnc =",i5,11x,
     &            " # (D=0) user defined multiplier, 0(no), 1, 2")')
     &            itdfn(m,1)

*-----------------------------------------------------------------------

            if (rtdre(m).gt.0.0) then

                  write(iot,'("   dresol =",es15.7,1x,
     &            " # (D=0) width = sqrt(dresol**2 + dfano*E)")')
     &            rtdre(m)

            end if

*-----------------------------------------------------------------------

            if (rtdfa(m).gt.0.0) then

                  write(iot,'("    dfano =",es15.7,1x,
     &            " # (D=0) width = sqrt(dresol**2 + dfano*E)")')
     &            rtdfa(m)

            end if

*-----------------------------------------------------------------------


               if( itmtn(m) .eq. 0 ) then

                  write(iot,'(" material =  all"11x,
     &            " # (D=all) number of specific material")')

               else if( itmtn(m) .gt. 0 ) then

                  write(iot,'(" material =",i5,11x,
     &            " # (D=all) number of specific material")')
     &            itmtn(m) * itmcn(m)

                     k = 0

                  do j = 1, itmtn(m)

                     k = k + 1

                     jmat(k) = ismte( itmtt(m) + j - 1 )

                     if( k / 10 * 10 .eq. k .or. j .eq. itmtn(m) ) then

                        write(iot,'(11x,10(i6))') ( jmat(i),i = 1, k )

                        k = 0

                     end if

                  end do

               end if

*-----------------------------------------------------------------------

               if( itout(m) .eq. 1 ) then

                  write(iot,'("   output =  dose",9x,
     &            "  # total deposit energy")')

               else if( itout(m) .eq. 2 ) then

                  write(iot,'("   output = deposit",7x,
     &            "  # deposit enerygy distribution")')

               end if

               write(iot,'("  deposit = ",i4,11x,
     &         " # (D=0) 0-> total deposit dist,",
     &                 " 1-> each process")')
     &         itdpo(m)

*-----------------------------------------------------------------------

               if( itout(m) .ge. 2 ) then

                  call echmty(0,iot,'e',itety(m),itenm(m),
     &                        rtedl(m),rtemi(m),rtema(m),iterg(m)
     &                        ,abs(itenm(m)+1),das_iterg)

               end if

*-----------------------------------------------------------------------

               if( ittty(m) .gt. 0 ) then

                  call echmty(0,iot,'t',ittty(m),ittnm(m),
     &                        rttdl(m),rttmi(m),rttma(m),ittrg(m)
     &                        ,abs(ittnm(m)+1),das_ittrg)

               end if

*-----------------------------------------------------------------------

               idtyp = 0

            do i = 1, icn

               if( itaxs(m,ian(i)) .ge. 7 .and.
     &             itaxs(m,ian(i)) .le. 10 )  idtyp = ittwo(m)

            end do

            if( idtyp .ne. 0 ) then

               write(iot,'("  2D-type = ",i4,11x,
     &         " # 1:Cont, 2:Clust, 3:Color, 4:xyz, 5:mat,",
     &         " 6:Clust+Cont, 7:Col+Cont")') idtyp

            end if

*-----------------------------------------------------------------------

            do i = 1, icn

                  iax = ian(i)

                  write(iot,'("     axis = ",a5,10x,
     &            " # axis of output")')
     &            aname( itaxs(m,iax) )

            if ( itfll(m,iax) .gt. 0 ) then

                  msfile = max( 14, itfll(m,iax) )

                  write(iot,'("     file = ",100a1)')
     &                  ( ctfln(m,iax)(j:j), j = 1, msfile ),
     &                  ( asfil(j:j), j = 1, 42 )

            end if

            end do

               call echprt(1,iot,m,itpan,itpat,jtpat,itnm,itmxpt,itmxgp)

*-----------------------------------------------------------------------

               if( rtfac(m) .ne. 1.0d0 ) then

                  write(iot,'("   factor =",1p1g14.7,1x,1x,
     &            " # (D=1.0) normalization factor")')
     &            rtfac(m)

               end if

*-----------------------------------------------------------------------

               if( itmtn(m) .gt. 0 .and. itmsh(m) .eq. 3 ) then

                  write(iot,'("   volmat = ",i4,11x,
     &            " # (D=9) vol.corr., 0:no, number of scan / dir.")')
     &            itvm(m)

               end if

*-----------------------------------------------------------------------

               if( itanl(m) .gt. 0 ) then

                  write(iot,'("    angel = ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               end if

*-----------------------------------------------------------------------

               if( itsans(m) .gt. 0 ) then

                  write(iot,'("   sangel = ",i2)') itsans(m)
                  call write_sangel(iot,m,1)

               end if

*-----------------------------------------------------------------------

               if( itaxl(m) .gt. 0 ) then

                  write(iot,'("    x-txt = ",200a1)')
     &            ( itaxt(m)(i:i),i = 1, itaxl(m) )

               end if

*-----------------------------------------------------------------------

               if( itayl(m) .gt. 0 ) then

                  write(iot,'("    y-txt = ",200a1)')
     &            ( itayt(m)(i:i),i = 1, itayl(m) )

               end if

*-----------------------------------------------------------------------

               if( itazl(m) .gt. 0 ) then

                  write(iot,'("    z-txt = ",200a1)')
     &            ( itazt(m)(i:i),i = 1, itazl(m) )

               end if

*-----------------------------------------------------------------------

               if( itgsh(m) .gt. 0 ) then

                  write(iot,'("    gshow = ",i4,10x,
     &            "  # 0: no 1:bnd, 2:bnd+mat, ",
     &     "3:bnd+reg 4:bnd+lat 5:bmp style")')
     &            itgsh(m)

               end if

*-----------------------------------------------------------------------

               if( itglt(m) .ge. 0 ) then

                  write(iot,'("    gslat = ",i4,10x,
     &            "  # 0: no line for lat, 1: with,",
     &               " 2: no line between the same cell")')
     &            itglt(m)

               end if

*-----------------------------------------------------------------------

               if( itgsh(m) .gt. 0 .and. itger(m) .gt. 0 ) then

                  write(iot,'("    ginfo = ",i4,11x,
     &            " # (D=0) 0:no, 1:show in graph, 2:add file")')
     &            itger(m)

               end if

*-----------------------------------------------------------------------

               if( itrsh(m) .gt. 0 ) then

                  write(iot,'("    rshow = ",i4,10x,
     &            "  # 0: no 1:bnd, 2:bnd+mat, ",
     &     "3:bnd+reg 4:bnd+lat 5:bmp style")')
     &            itrsh(m)

                  call echmty(1,iot,'x',itxty(m),itxnm(m),
     &                        rtxdl(m),rtxmi(m),rtxma(m),itxrg(m)
     &                        ,abs(itxnm(m)+1),das_itxrg)
                  call echmty(1,iot,'y',ityty(m),itynm(m),
     &                        rtydl(m),rtymi(m),rtyma(m),ityrg(m)
     &                        ,abs(itynm(m)+1),das_ityrg)
                  call echmty(1,iot,'z',itzty(m),itznm(m),
     &                        rtzdl(m),rtzmi(m),rtzma(m),itzrg(m)
     &                        ,abs(itznm(m)+1),das_itzrg)

               end if

*-----------------------------------------------------------------------

               if( itgsh(m) .gt. 0 .or. itrsh(m) .gt. 0 ) then

                  write(iot,'("    resol = ",i4,10x,
     &            "  # (D=1) resolution of gshow or rshow")')
     &            itres(m)

                  write(iot,'("    width =",1p1g14.7,1x,1x,
     &            " # (D=0.5) width of lines for gshow or rshow")')
     &            rtwid(m)

               end if

*-----------------------------------------------------------------------

               if( iterl(m) .ne. 72 ) then

                  write(iot,'("   iechrl =",i5,10x,
     &            "  # (D=72) length reg of val or vol")')
     &            iterl(m)

               end if

*-----------------------------------------------------------------------

               if( iteps(m) .ne. 0 ) then

                  write(iot,'("   epsout =",i5,10x,
     &            "  # (D=0) generate eps file by ANGEL")')
     &            iteps(m)

               end if

*-----------------------------------------------------------------------

               if( itbmp(m) .ne. 0 ) then

                  write(iot,'("   bmpout =",i5,10x,
     &            "  # (D=0) generate bitmap file")')
     &            itbmp(m)

               end if

*-----------------------------------------------------------------------

               if( itvtk(m) .ne. 0 ) then

                  write(iot,'("   vtkout =",i5,10x,
     &            "  # (D=0) generate vtk file")')
     &            itvtk(m)

                  if ( itvtkfmt(m) .eq. 0 ) then

                     write(iot,'("   vtkfmt =",i5,10x,
     &               "  # (D=0) generate vtk ascii file")')
     &               itvtkfmt(m)

                  else

                     write(iot,'("   vtkfmt =",i5,10x,
     &               "  # (D=0) generate vtk binary file")')
     &               itvtkfmt(m)

                  end if

               end if

*-----------------------------------------------------------------------

               if( itnlatcel(m) .ne. 0 ) then

                  write(iot,'("  nlatcel =",i5,10x,
     &            "  # (D=0) lattice cell ID for",
     &            " individual dealing in lattice")')
     &            itnlatcel(m)

                  write(iot,'("  nlatmem =",i5,10x,
     &            "  # (D=0) number of memory for",
     &            " individual dealing in lattice")')
     &            itnlatmem(m)

               end if

*-----------------------------------------------------------------------
               if( itfoam(m) .ne. 0 ) then

                  write(iot,'("  foamout =",i5,10x,
     &            "  # (D=0) generate OpenFOAM data")')
     &            itfoam(m)

               end if

*-----------------------------------------------------------------------
*        check of counter
*-----------------------------------------------------------------------

            do i = 1, 3

               if( itcnt(i,m) .ne. 0 ) then

                  write(iot,'(" ctmin(",i1,") =",i5,10x,
     &            "  # (D=-9999) min. of counter ",i1)')
     &            i, itcnt(i*2+2,m), i

                  write(iot,'(" ctmax(",i1,") =",i5,10x,
     &            "  # (D= 9999) max. of counter ",i1)')
     &            i, itcnt(i*2+3,m), i

                  if( ncntc(i) .eq. 0 )
     &            write(iot,'("# Warning: there is no counter = ",
     &                        i1," definition in [counter].")') i

               end if

            end do

            do i = 1, 3
               if( itcnth(i,m) .ne. 0 ) then
                  write(iot,'(" chmin(",i1,") =",i5,10x,
     &            "  # (D=-9999) min. of history counter ",i1)')
     &            i, itcnth(i*2+2,m), i
                  write(iot,'(" chmax(",i1,") =",i5,10x,
     &            "  # (D= 9999) max. of history counter ",i1)')
     &            i, itcnth(i*2+3,m), i
               end if
            end do

*-----------------------------------------------------------------------

               if( rtstd(m) .gt. 0.d0 ) then

                  write(iot,'("   stdcut =",1p1g14.7,1x,1x,
     &            " # (D=-1) threshold value of STD cutoff")')
     &            rtstd(m)

               end if

*-----------------------------------------------------------------------
               if( itmxang(m) .ne. 0 ) then

                  write(iot,'(" maxangel =",i5,10x,
     &            "  # (D=0) number of output particles")')
     &            itmxang(m)

               end if

*-----------------------------------------------------------------------

               if( itman(m) .eq. 0 ) then

                  write(iot,'("   mother =  all"11x,
     &            " # (D=all) number of specific mother")')

               else if( itman(m) .gt. 0 ) then

                  write(iot,'("   mother = ",i4,11x,
     &            " # (D=all) number of specific mother")')
     &            itman(m) * itmct(m)

                     k = 0

                  do j = 1, itman(m)

                     k = k + 1

                     iaz = ismat( itmat(m) + j - 1 )

                     iz = iaz / 1000
                     ia = iaz - iz * 1000

                     if( ia .eq. 0 ) then

                        chaus(k) = elmnt(iz)

                     else

                        call chname(idum,ia,iz,chau)

                        chaus(k) = chau

                     end if

                     if( k / 7 * 7 .eq. k .or. j .eq. itman(m) ) then

                        write(iot,'(12x,7(a8))') ( chaus(i),i = 1, k )

                        k = 0

                     end if

                  end do

               end if

*-----------------------------------------------------------------------
*        sumtally subsection
*-----------------------------------------------------------------------
            if ( icntl==13 .and. nsumtalRead(m)==13 ) then

               call sumtal_echo(iot,m)

            end if

*-----------------------------------------------------------------------

            write(iot,'("#    used :",
     &                  "        main (  %)",
     &                  "        temp (  %)",
     &                  "       total (  %)"/
     &                  "#  memory :",
     &                  3(i12," (",i3,")"))')
     &                  itsmn(m), nint(dble(itsmn(m))/dble(mdas)*100),
     &                  itstm(m), nint(dble(itstm(m))/dble(mdas)*100),
     &                  itsmn(m)+itstm(m),
     &                  nint(dble(itsmn(m)+itstm(m))/dble(mdas)*100)

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine letech(iot,m,iaxs,icax)
*                                                                      *
*       input echo of the let tally                                    *
*       last modified by K.Niita on 2005/08/31                         *
*                                                                      *
************************************************************************

      use sumtallymod
      use sangelmod, only: itsans
      use partmod ! frtati 2021/10/05
      use moddas_mesh
      use moddas_region
      use moddas_tally
*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'

*-----------------------------------------------------------------------

      common /tall00/ itnm, ital(itlmax), itals(itlmax), italm(itlmax)
      common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                itout(itlmax), ittwo(itlmax)

      common /tall02/ itrgn(itlmax), itrgm(itlmax), itreg(itlmax)
      common /tall24/ itrcm(itlmax), itrcg(itlmax)
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

      common /tall16/ itnfr(itlmax), itnfz(itlmax), itnfn(itlmax)

      common /tall14/ ittll(itlmax), ittle(itlmax)
      character ittle*80

      common /tall15/ itanl(itlmax), itang(itlmax)
      character itang*200

      common /tall27/ itaxl(itlmax), itaxt(itlmax),
     &                itayl(itlmax), itayt(itlmax),
     &                itazl(itlmax), itazt(itlmax)
      character itaxt*200, itayt*200, itazt*200

      common /tall19/ itsmn(itlmax), itstm(itlmax)

      common /tall21/ rtfac(itlmax)
      common /tall28/ itgsh(itlmax)
      common /tall29/ itrsh(itlmax)
      common /tall31/ iterl(itlmax)
      common /tall22/ itmcn(itlmax), itmtn(itlmax), itmtt(itlmax)
      common /tall34/ itvm(itlmax)
      common /tall35/ iteps(itlmax)
      common /tall37/ itcnt(9,itlmax)
      common /cntreg/ ncntc(3), ncreg(3), incrc(3), incrt(3), kcont(3)
      common /cntpat/ icpan(3), icpat(3,20,2)

      common /tall38/ itres(itlmax)
      common /tall39/ rtwid(itlmax)

      common /tall46/ itmtr(itlmax,4), rtmtr(itlmax,13)
      common /tall49/ itglt(itlmax)
      common /tall50/ itlmt(itlmax), ite2l(itlmax)   ! CCSE 2022.08.31

      common /tall60/ itger(itlmax)

      common /tall67/ itstd(itlmax), rtstd(itlmax)

      common /tall82/ itcnth(9,itlmax)

      dimension idas(mdas*2) !20220906frtati temporary for fbounds-check
      equivalence ( das, idas )

*-----------------------------------------------------------------------

      common /inpec/ ititl, ipara, ibody, iregn, llarr, itby, itar

*-----------------------------------------------------------------------
       integer   icntl
       integer   inucr

       common /tcntl/  icntl, inucr
!-----------------------------------------------------------------------

      character asfil*100

      character aname(10)*3
      data aname / 'let','reg','  x','  y','  z','  r',
     &             ' xy',' yz',' xz',' rz'/

      dimension ian(6)
      dimension jmat(20)

*-----------------------------------------------------------------------
      common /tall79/ itallech
      if( itallech .eq. 0 ) return

*-----------------------------------------------------------------------
*     input echo
*-----------------------------------------------------------------------

            write(iot,'("[ T-LET ]")')

                  write(iot,'("    title = ",80a1)')
     &            ( ittle(m)(i:i), i = 1, ittll(m) )

            if( itmsh(m) .eq. 1 ) then

                  write(iot,'("     mesh =  reg",11x,
     &            " # mesh type is region-wise")')

                     idas1 = mmmax
                     idas2 = ( idas1 + itrgn(m) - 1 ) * 2 + 1
                     idas3 = idas2 + itrgn(m)

                  call echrg(1,iterl(m),0,1,iot,itrcm(m),
     &                       idas_itrcg(itrcg(m)),
     &                       itrgn(m),itrgm(m),idas_itreg(itreg(m)),
     &                       das(idas1),idas(idas2),
     &                       itrnv(m),idas(itriv(m)),das(itrrv(m)),
     &                       idas3)

            else if( itmsh(m) .eq. 2 ) then

                  write(iot,'("     mesh =  r-z",11x,
     &            " # mesh type is r-z scoring mesh")')

                  write(iot,'("       x0 = ",1p1g14.7,1x,
     &            " # center x-position of r-z mesh")') rtrx0(m)

                  write(iot,'("       y0 = ",1p1g14.7,1x,
     &            " # center y-position of r-z mesh")') rtry0(m)

                  call echmty(0,iot,'r',itrty(m),itrnm(m),
     &                        rtrdl(m),rtrmi(m),rtrma(m),itrrg(m)
     &                        ,abs(itrnm(m)+1),das_itrrg)
                  call echmty(0,iot,'z',itzty(m),itznm(m),
     &                        rtzdl(m),rtzmi(m),rtzma(m),itzrg(m)
     &                        ,abs(itznm(m)+1),das_itzrg)

            else if( itmsh(m) .eq. 3 ) then

                  write(iot,'("     mesh =  xyz",11x,
     &            " # mesh type is xyz scoring mesh")')

                  call echmty(0,iot,'x',itxty(m),itxnm(m),
     &                        rtxdl(m),rtxmi(m),rtxma(m),itxrg(m)
     &                        ,abs(itxnm(m)+1),das_itxrg)
                  call echmty(0,iot,'y',ityty(m),itynm(m),
     &                        rtydl(m),rtymi(m),rtyma(m),ityrg(m)
     &                        ,abs(itynm(m)+1),das_ityrg)
                  call echmty(0,iot,'z',itzty(m),itznm(m),
     &                        rtzdl(m),rtzmi(m),rtzma(m),itzrg(m)
     &                        ,abs(itznm(m)+1),das_itzrg)

            end if

*-----------------------------------------------------------------------

            if( itmtr(m,1) .gt. 0 ) then

               if( itmtr(m,1) .eq. 1 ) then

                     write(iot,'("     trcl = ",i6,9x,
     &               " # transform ID")') itmtr(m,3)

               else

                  if( itmtr(m,2) .eq. 0 ) then

                     write(iot,'("     trcl = ",
     &                           3(1p1e15.7),
     &                           3(/12x,3(1p1e15.7)),i5)')
     &               ( rtmtr(m,j), j = 1, 12 ),
     &               nint( rtmtr(m,13) )

                  else

                     write(iot,'("    *trcl = ",
     &                           3(1p1e15.7),
     &                           3(/12x,3(1p1e15.7)),i5)')
     &               ( rtmtr(m,j), j = 1, 12 ),
     &               nint( rtmtr(m,13) )

                  end if

               end if

            end if

*-----------------------------------------------------------------------

                  call echmty(0,iot,'l',itety(m),itenm(m),
     &                        rtedl(m),rtemi(m),rtema(m),iterg(m)
     &                        ,abs(itenm(m)+1),das_iterg)

*-----------------------------------------------------------------------

               if( itunt(m) .eq. 1 ) then

                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [cm/(keV/um)/source]")') itunt(m)

               else if( itunt(m) .eq. 2 ) then

                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [MeV/(keV/um)/source]")') itunt(m)

               else if( itunt(m) .eq. 3 ) then

                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [cm/ln(keV/um)/source]")') itunt(m)

               else if( itunt(m) .eq. 4 ) then

                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [MeV/ln(keV/um)/source]")') itunt(m)

               else if( itunt(m) .eq. 5 ) then

                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [cm/source]")') itunt(m)

               else if( itunt(m) .eq. 6 ) then

                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [MeV/source]")') itunt(m)

               else if( itunt(m) .eq. 7 ) then

                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [1/cm^2/(keV/um)/source]")') itunt(m)

               else if( itunt(m) .eq. 8 ) then

                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [MeV/cm^3/(keV/um)/source]")') itunt(m)

               else if( itunt(m) .eq. 9 ) then

                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [1/cm^2/ln(keV/um)/source]")') itunt(m)

               else if( itunt(m) .eq. 10 ) then

                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [MeV/cm^3/ln(keV/um)/source]")')
     &               itunt(m)

               else if( itunt(m) .eq. 11 ) then

                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [1/cm^2/source]")') itunt(m)

               else if( itunt(m) .eq. 12 ) then

                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [MeV/cm^3/source]")') itunt(m)

               else if( itunt(m) .eq. 13 ) then ! T.Sato 2017/11/21

                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is L*f(L) [dimensionless]")') itunt(m)

               else if( itunt(m) .eq. 14 ) then ! T.Sato 2017/11/21

                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is L*d(L) [keV/um]")') itunt(m)

               end if

*-----------------------------------------------------------------------

                  asfil = '  # file name of output for the above axis'

            if( icax .eq. 0 ) then

                  icn = itaxn(m)

               do i = 1, icn

                  ian(i) = i

               end  do

            else

                  icn = 1
                  ian(1) = iaxs

            end if

*-----------------------------------------------------------------------

                  write(iot,'("   letmat =",i5,11x,
     &            " # (D=0) mat ID for LET, 0:real mat, ",
     &                                    "<0: electron for H2O")')
     &            itlmt(m)

*-----------------------------------------------------------------------

               if( itmtn(m) .eq. 0 ) then

                  write(iot,'(" material =  all"11x,
     &            " # (D=all) number of specific material")')

               else if( itmtn(m) .gt. 0 ) then

                  write(iot,'(" material =",i5,11x,
     &            " # (D=all) number of specific material")')
     &            itmtn(m) * itmcn(m)

                     k = 0

                  do j = 1, itmtn(m)

                     k = k + 1

                     jmat(k) = ismte( itmtt(m) + j - 1 )

                     if( k / 10 * 10 .eq. k .or. j .eq. itmtn(m) ) then

                        write(iot,'(11x,10(i6))') ( jmat(i),i = 1, k )

                        k = 0

                     end if

                  end do

               end if

*-----------------------------------------------------------------------

               idtyp = 0

            do i = 1, icn

               if( itaxs(m,ian(i)) .ge. 7 .and.
     &             itaxs(m,ian(i)) .le. 10 )  idtyp = ittwo(m)

            end do

            if( idtyp .ne. 0 ) then

               write(iot,'("  2D-type = ",i4,11x,
     &         " # 1:Cont, 2:Clust, 3:Color, 4:xyz, 5:mat,",
     &         " 6:Clust+Cont, 7:Col+Cont")') idtyp

            end if

*-----------------------------------------------------------------------

            do i = 1, icn

                  iax = ian(i)

                  write(iot,'("     axis =  ",a3,11x,
     &            " # axis of output")')
     &            aname( itaxs(m,iax) )

             if ( itfll(m,iax) .gt. 0 ) then

                  msfile = max( 14, itfll(m,iax) )

                  write(iot,'("     file = ",100a1)')
     &                  ( ctfln(m,iax)(j:j), j = 1, msfile ),
     &                  ( asfil(j:j), j = 1, 42 )

             end if

            end do

               call echprt(1,iot,m,itpan,itpat,jtpat,itnm,itmxpt,itmxgp)

*-----------------------------------------------------------------------

               if( rtfac(m) .ne. 1.0d0 ) then

                  write(iot,'("   factor =",1p1g14.7,1x,1x,
     &            " # (D=1.0) normalization factor")')
     &            rtfac(m)

               end if

*-----------------------------------------------------------------------

               if( itmtn(m) .gt. 0 .and. itmsh(m) .eq. 3 ) then

                  write(iot,'("   volmat = ",i4,11x,
     &            " # (D=9) vol.corr., 0:no, number of scan / dir.")')
     &            itvm(m)

               end if

*-----------------------------------------------------------------------

               if( itanl(m) .gt. 0 ) then

                  write(iot,'("    angel = ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               end if

*-----------------------------------------------------------------------

               if( itsans(m) .gt. 0 ) then

                  write(iot,'("   sangel = ",i2)') itsans(m)
                  call write_sangel(iot,m,1)

               end if

*-----------------------------------------------------------------------

               if( itaxl(m) .gt. 0 ) then

                  write(iot,'("    x-txt = ",200a1)')
     &            ( itaxt(m)(i:i),i = 1, itaxl(m) )

               end if

*-----------------------------------------------------------------------

               if( itayl(m) .gt. 0 ) then

                  write(iot,'("    y-txt = ",200a1)')
     &            ( itayt(m)(i:i),i = 1, itayl(m) )

               end if

*-----------------------------------------------------------------------

               if( itazl(m) .gt. 0 ) then

                  write(iot,'("    z-txt = ",200a1)')
     &            ( itazt(m)(i:i),i = 1, itazl(m) )

               end if

*-----------------------------------------------------------------------

               if( itgsh(m) .gt. 0 ) then

                  write(iot,'("    gshow = ",i4,10x,
     &            "  # 0: no 1:bnd, 2:bnd+mat, ",
     &     "3:bnd+reg 4:bnd+lat 5:bmp style")')
     &            itgsh(m)

               end if

*-----------------------------------------------------------------------

               if( itglt(m) .ge. 0 ) then

                  write(iot,'("    gslat = ",i4,10x,
     &            "  # 0: no line for lat, 1: with,",
     &               " 2: no line between the same cell")')
     &            itglt(m)

               end if

*-----------------------------------------------------------------------

               if( itgsh(m) .gt. 0 .and. itger(m) .gt. 0 ) then

                  write(iot,'("    ginfo = ",i4,11x,
     &            " # (D=0) 0:no, 1:show in graph, 2:add file")')
     &            itger(m)

               end if

*-----------------------------------------------------------------------

               if( itrsh(m) .gt. 0 ) then

                  write(iot,'("    rshow = ",i4,10x,
     &            "  # 0: no 1:bnd, 2:bnd+mat, ",
     &     "3:bnd+reg 4:bnd+lat 5:bmp style")')
     &            itrsh(m)

                  call echmty(1,iot,'x',itxty(m),itxnm(m),
     &                        rtxdl(m),rtxmi(m),rtxma(m),itxrg(m)
     &                        ,abs(itxnm(m)+1),das_itxrg)
                  call echmty(1,iot,'y',ityty(m),itynm(m),
     &                        rtydl(m),rtymi(m),rtyma(m),ityrg(m)
     &                        ,abs(itynm(m)+1),das_ityrg)
                  call echmty(1,iot,'z',itzty(m),itznm(m),
     &                        rtzdl(m),rtzmi(m),rtzma(m),itzrg(m)
     &                        ,abs(itznm(m)+1),das_itzrg)

               end if

*-----------------------------------------------------------------------

               if( itgsh(m) .gt. 0 .or. itrsh(m) .gt. 0 ) then

                  write(iot,'("    resol = ",i4,10x,
     &            "  # (D=1) resolution of gshow or rshow")')
     &            itres(m)

                  write(iot,'("    width =",1p1g14.7,1x,1x,
     &            " # (D=0.5) width of lines for gshow or rshow")')
     &            rtwid(m)

               end if

*-----------------------------------------------------------------------

               if( iterl(m) .ne. 72 ) then

                  write(iot,'("   iechrl =",i5,10x,
     &            "  # (D=72) length reg of val or vol")')
     &            iterl(m)

               end if

*-----------------------------------------------------------------------

               if( iteps(m) .ne. 0 ) then

                  write(iot,'("   epsout =",i5,10x,
     &            "  # (D=0) generate eps file by ANGEL")')
     &            iteps(m)

               end if

*-----------------------------------------------------------------------
*        check of counter
*-----------------------------------------------------------------------

            do i = 1, 3

               if( itcnt(i,m) .ne. 0 ) then

                  write(iot,'(" ctmin(",i1,") =",i5,10x,
     &            "  # (D=-9999) min. of counter ",i1)')
     &            i, itcnt(i*2+2,m), i

                  write(iot,'(" ctmax(",i1,") =",i5,10x,
     &            "  # (D= 9999) max. of counter ",i1)')
     &            i, itcnt(i*2+3,m), i

                  if( ncntc(i) .eq. 0 )
     &            write(iot,'("# Warning: there is no counter = ",
     &                        i1," definition in [counter].")') i

               end if

            end do

            do i = 1, 3
               if( itcnth(i,m) .ne. 0 ) then
                  write(iot,'(" chmin(",i1,") =",i5,10x,
     &            "  # (D=-9999) min. of history counter ",i1)')
     &            i, itcnth(i*2+2,m), i
                  write(iot,'(" chmax(",i1,") =",i5,10x,
     &            "  # (D= 9999) max. of history counter ",i1)')
     &            i, itcnth(i*2+3,m), i
               end if
            end do

*-----------------------------------------------------------------------

               if( rtstd(m) .gt. 0.d0 ) then

                  write(iot,'("   stdcut =",1p1g14.7,1x,1x,
     &            " # (D=-1) threshold value of STD cutoff")')
     &            rtstd(m)

               end if

*-----------------------------------------------------------------------
               if( itmxang(m) .ne. 0 ) then

                  write(iot,'(" maxangel =",i5,10x,
     &            "  # (D=0) number of output particles")')
     &            itmxang(m)

               end if

*-----------------------------------------------------------------------
*        sumtally subsection
*-----------------------------------------------------------------------
            if ( icntl==13 .and. nsumtalRead(m)==12 ) then

               call sumtal_echo(iot,m)

            end if
*-----------------------------------------------------------------------

            write(iot,'("#    used :",
     &                  "        main (  %)",
     &                  "        temp (  %)",
     &                  "       total (  %)"/
     &                  "#  memory :",
     &                  3(i12," (",i3,")"))')
     &                  itsmn(m), nint(dble(itsmn(m))/dble(mdas)*100),
     &                  itstm(m), nint(dble(itstm(m))/dble(mdas)*100),
     &                  itsmn(m)+itstm(m),
     &                  nint(dble(itsmn(m)+itstm(m))/dble(mdas)*100)

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine trckech(iot,m,iaxs,icax)
*                                                                      *
*       input echo of the track length tally                           *
*       last modified by K.Niita on 2014/12/04                         *
*                                                                      *
************************************************************************

      use sumtallymod
      use sangelmod, only: itsans
      use partmod ! frtati 2021/10/05
      use moddas_mesh
      use moddas_region
      use moddas_tally

*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'

*-----------------------------------------------------------------------

      common /tall00/ itnm, ital(itlmax), itals(itlmax), italm(itlmax)
      common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                itout(itlmax), ittwo(itlmax)

      common /tall02/ itrgn(itlmax), itrgm(itlmax), itreg(itlmax)
      common /tall24/ itrcm(itlmax), itrcg(itlmax)
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

      common /tall16/ itnfr(itlmax), itnfz(itlmax), itnfn(itlmax)

      common /tall14/ ittll(itlmax), ittle(itlmax)
      character ittle*80

      common /tall15/ itanl(itlmax), itang(itlmax)
      character itang*200

      common /tall27/ itaxl(itlmax), itaxt(itlmax),
     &                itayl(itlmax), itayt(itlmax),
     &                itazl(itlmax), itazt(itlmax)
      character itaxt*200, itayt*200, itazt*200

      common /tall19/ itsmn(itlmax), itstm(itlmax)

      common /tall23/ itaty(itlmax), itanm(itlmax), itarg(itlmax),
     &                rtami(itlmax), rtama(itlmax), rtadl(itlmax)


      common /tall20/ ittty(itlmax), ittnm(itlmax), ittrg(itlmax),
     &                rttmi(itlmax), rttma(itlmax), rttdl(itlmax)
      common /tall21/ rtfac(itlmax)
      common /tall28/ itgsh(itlmax)
      common /tall29/ itrsh(itlmax)
      common /tall31/ iterl(itlmax)
      common /tall22/ itmcn(itlmax), itmtn(itlmax), itmtt(itlmax)
      common /tall34/ itvm(itlmax)
      common /tall35/ iteps(itlmax)
      common /tall37/ itcnt(9,itlmax)
      common /cntreg/ ncntc(3), ncreg(3), incrc(3), incrt(3), kcont(3)
      common /cntpat/ icpan(3), icpat(3,20,2)

      common /tall38/ itres(itlmax)
      common /tall39/ rtwid(itlmax)

      common /tall45/ itmlp(itlmax), itmln(itlmax,6), itmst(itlmax),
     &                itmli(itlmax,6), rtmme(itlmax,6), itmnt(itlmax,6),
     &                itmpn(itlmax,6), itmpt(itlmax,6,6,2)
      common /tall46/ itmtr(itlmax,4), rtmtr(itlmax,13)
      common /tall49/ itglt(itlmax)

      common /tall60/ itger(itlmax)
      common /tall63/ itbmp(itlmax)
      common /tall64/ itvtk(itlmax),itvtkfmt(itlmax)

      common /tall67/ itstd(itlmax), rtstd(itlmax)

*-----------------------------------------------------------------------
      common /tall78/ itism(itlmax), itist(10,itlmax), itjst(10,itlmax),
     &                itkst(10,itlmax), itstt(itlmax), itsdd(itlmax)

      character schan(10)*5
      data ( schan(i), i = 1, 10 ) /
     &    'ireg ','ix   ','iy   ','iz   ','ir   ',
     &    'ie   ','it   ','ia   ','ipart','imul '/
      character ch73*200

*-----------------------------------------------------------------------

      character yen*1

*-----------------------------------------------------------------------
      integer :: itfoam
      common /tall76/ itfoam(itlmax)

      common /tall82/ itcnth(9,itlmax)
cfrtati 2023/12/07
      common /tall84/ itmto(itlmax,6,4)

      common /tall91/ itextstat(itlmax), mftal(itlmax) ! S.H. extstat 2024.3.18
      common /tall93/ tprodenmn(itlmax), tprodenmx(itlmax),
     &     nbtproden(itlmax), itprodenchk(itlmax) !S.H. extstat 2024.4.28

      dimension idas(mdas*2) !20220906frtati temporary for fbounds-check
      equivalence ( das, idas )


cfrtati 2021/10/05 moved to partmod

      dimension ilgt(6), jlgt(6)
      character ch13*13, ch6*6
      character ch200*400
      character ch100*400

*-----------------------------------------------------------------------

      common /inpec/ ititl, ipara, ibody, iregn, llarr, itby, itar

*-----------------------------------------------------------------------
       integer   icntl
       integer   inucr

       common /tcntl/  icntl, inucr
!-----------------------------------------------------------------------

      character asfil*100

cFURUTA20190110 !FURUTA20200601
      character aname(15)*6
      data aname / '   eng','   reg','     x','     y','     z',
     &     '     r','    xy','    yz','    xz','    rz','     t',
     &     '   rad','   deg','   tet','dchain'/

      dimension ian(6)
      dimension jmat(20)

      integer           iMeVperu
      common /cMeVperu/ iMeVperu

*-----------------------------------------------------------------------
      common /tall79/ itallech
      if( itallech .eq. 0 ) return

*-----------------------------------------------------------------------
*     input echo
*-----------------------------------------------------------------------
            if( ital(m+1) .eq. 16 .and. icax .eq. 0 ) then
                  write(iot,'("[ T-Track ] off")')
            else
                  write(iot,'("[ T-Track ]")')
            end if

                  write(iot,'("    title = ",80a1)')
     &            ( ittle(m)(i:i), i = 1, ittll(m) )

            if( itmsh(m) .eq. 1 ) then

                  write(iot,'("     mesh =  reg",11x,
     &            " # mesh type is region-wise")')

                     idas1 = mmmax
                     idas2 = ( idas1 + itrgn(m) - 1 ) * 2 + 1
                     idas3 = idas2 + itrgn(m)

                  call echrg(1,iterl(m),0,1,iot,itrcm(m),
     &                       idas_itrcg(itrcg(m)),
     &                       itrgn(m),itrgm(m),idas_itreg(itreg(m)),
     &                       das(idas1),idas(idas2),
     &                       itrnv(m),idas(itriv(m)),das(itrrv(m)),
     &                       idas3)

            else if( itmsh(m) .eq. 2 ) then

                  write(iot,'("     mesh =  r-z",11x,
     &            " # mesh type is r-z scoring mesh")')

                  write(iot,'("       x0 = ",1p1g14.7,1x,
     &            " # center x-position of r-z mesh")') rtrx0(m)

                  write(iot,'("       y0 = ",1p1g14.7,1x,
     &            " # center y-position of r-z mesh")') rtry0(m)

                  call echmty(0,iot,'r',itrty(m),itrnm(m),
     &                        rtrdl(m),rtrmi(m),rtrma(m),itrrg(m)
     &                        ,abs(itrnm(m)+1),das_itrrg)
                  call echmty(0,iot,'z',itzty(m),itznm(m),
     &                        rtzdl(m),rtzmi(m),rtzma(m),itzrg(m)
     &                        ,abs(itznm(m)+1),das_itzrg)

*-----------------------------------------------------------------------
               if( itaty(m) .ne. 0 ) then

                  call echmty(0,iot,'a',itaty(m),itanm(m),
     &                        rtadl(m),rtami(m),rtama(m),itarg(m)
     &                        ,abs(itanm(m)+1),das_itarg)

               end if
*-----------------------------------------------------------------------

            else if( itmsh(m) .eq. 3 ) then

                  write(iot,'("     mesh =  xyz",11x,
     &            " # mesh type is xyz scoring mesh")')

                  call echmty(0,iot,'x',itxty(m),itxnm(m),
     &                        rtxdl(m),rtxmi(m),rtxma(m),itxrg(m)
     &                        ,abs(itxnm(m)+1),das_itxrg)
                  call echmty(0,iot,'y',ityty(m),itynm(m),
     &                        rtydl(m),rtymi(m),rtyma(m),ityrg(m)
     &                        ,abs(itynm(m)+1),das_ityrg)
                  call echmty(0,iot,'z',itzty(m),itznm(m),
     &                        rtzdl(m),rtzmi(m),rtzma(m),itzrg(m)
     &                        ,abs(itznm(m)+1),das_itzrg)

*-----------------------------------------------------------------------
            else if( itmsh(m) .eq. 4 ) then

                  write(iot,'("     mesh =  tet",11x,
     &            " # mesh type is tetra geometry")')

                  call echtet(iot,itrgm(m),idas_itreg(itreg(m)))

            end if

*-----------------------------------------------------------------------

            if( itmtr(m,1) .gt. 0 ) then

               if( itmtr(m,1) .eq. 1 ) then

                     write(iot,'("     trcl = ",i6,9x,
     &               " # transform ID")') itmtr(m,3)

               else

                  if( itmtr(m,2) .eq. 0 ) then

                     write(iot,'("     trcl = ",
     &                           3(1p1e15.7),
     &                           3(/12x,3(1p1e15.7)),i5)')
     &               ( rtmtr(m,j), j = 1, 12 ),
     &               nint( rtmtr(m,13) )

                  else

                     write(iot,'("    *trcl = ",
     &                           3(1p1e15.7),
     &                           3(/12x,3(1p1e15.7)),i5)')
     &               ( rtmtr(m,j), j = 1, 12 ),
     &               nint( rtmtr(m,13) )

                  end if

               end if

            end if

*-----------------------------------------------------------------------

                  call echmty(0,iot,'e',itety(m),itenm(m),
     &                        rtedl(m),rtemi(m),rtema(m),iterg(m)
     &                        ,abs(itenm(m)+1),das_iterg)

*-----------------------------------------------------------------------

               if( ittty(m) .ne. 0 ) then

                  call echmty(0,iot,'t',ittty(m),ittnm(m),
     &                        rttdl(m),rttmi(m),rttma(m),ittrg(m)
     &                        ,abs(ittnm(m)+1),das_ittrg)

               end if

*-----------------------------------------------------------------------

               if( itunt(m) .eq. 1 ) then

                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [1/cm^2/source]")') itunt(m)

               else if( itunt(m) .eq. 2 ) then

                if ( iMeVperu.eq. 1 ) then
                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [1/cm^2/(MeV/n)/source]")') itunt(m)
                else
                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [1/cm^2/MeV/source]")') itunt(m)
                end if

               else if( itunt(m) .eq. 3 ) then

                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [1/cm^2/Lethargy/source]")') itunt(m)

               else if( itunt(m) .eq. 4 ) then

                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [1/cm^2/source] with vol=1cm^3")')
     &               itunt(m)

               else if( itunt(m) .eq. 11 ) then

                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [1/cm^2/nsec/source]")') itunt(m)

               else if( itunt(m) .eq. 12 ) then

                if ( iMeVperu.eq. 1 ) then
                  write(iot,'("     unit = ",i4,11x,
     &           " # unit is [1/cm^2/(MeV/n)/nsec/source]")') itunt(m)
                else
                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [1/cm^2/MeV/nsec/source]")') itunt(m)
                end if

               else if( itunt(m) .eq. 13 ) then

                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [1/cm^2/Lethargy/nsec/source]")')
     &               itunt(m)

               else if( itunt(m) .eq. 14 ) then

                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [1/cm^2/nsec/source] with vol=1cm^3")')
     &               itunt(m)

               end if

*-----------------------------------------------------------------------

                  asfil = '  # file name of output for the above axis'

            if( icax .eq. 0 ) then

                  icn = itaxn(m)

               do i = 1, icn

                  ian(i) = i

               end  do

            else

                  icn = 1
                  ian(1) = iaxs

            end if

*-----------------------------------------------------------------------

               if( itmtn(m) .eq. 0 ) then

                  write(iot,'(" material =  all"11x,
     &            " # (D=all) number of specific material")')

               else if( itmtn(m) .gt. 0 ) then

                  write(iot,'(" material =",i5,11x,
     &            " # (D=all) number of specific material")')
     &            itmtn(m) * itmcn(m)

                     k = 0

                  do j = 1, itmtn(m)

                     k = k + 1

                     jmat(k) = ismte( itmtt(m) + j - 1 )

                     if( k / 10 * 10 .eq. k .or. j .eq. itmtn(m) ) then

                        write(iot,'(11x,10(i6))') ( jmat(i),i = 1, k )

                        k = 0

                     end if

                  end do

               end if

*-----------------------------------------------------------------------

               idtyp = 0

            do i = 1, icn

               if( itaxs(m,ian(i)) .ge. 7 .and.
     &             itaxs(m,ian(i)) .le. 10 )  idtyp = ittwo(m)

            end do

            if( idtyp .ne. 0 ) then

               write(iot,'("  2D-type = ",i4,11x,
     &         " # 1:Cont, 2:Clust, 3:Color, 4:xyz, 5:mat,",
     &         " 6:Clust+Cont, 7:Col+Cont")') idtyp

            end if

*-----------------------------------------------------------------------

            do i = 1, icn

                  iax = ian(i)

                  write(iot,'("     axis =  ",a6,11x,
     &            " # axis of output")')
     &            aname( itaxs(m,iax) )

            if ( itfll(m,iax) .gt. 0 ) then

                  msfile = max( 14, itfll(m,iax) )

                  write(iot,'("     file = ",100a1)')
     &                  ( ctfln(m,iax)(j:j), j = 1, msfile ),
     &                  ( asfil(j:j), j = 1, 42 )

            end if

            end do

               call echprt(1,iot,m,itpan,itpat,jtpat,itnm,itmxpt,itmxgp)

*-----------------------------------------------------------------------
*        multiplier
*-----------------------------------------------------------------------

         if( itmlp(m) .gt. 0 ) then

            do k = 1, itmlp(m)

*-----------------------------------------------------------------------

               if( itmln(m,k) .eq. 1 .and.
     &             mltp(1,itmli(m,k)/13+1) .eq. -1 ) then

                  imall = 0

                  write(iot,'("  multiplier = all",9x,
     &            " # number of material group")')

               else

                  imall = itmln(m,k)

                  write(iot,'("  multiplier = ",i3,9x,
     &            " # number of material group")') imall

               end if

*-----------------------------------------------------------------------

                     ipan(m) = itmpn(m,k)

                  do j = 1, itmpn(m,k)

                     ipat(m,j,1) = itmpt(m,k,j,1)
                     ipat(m,j,2) = itmpt(m,k,j,2)

                  end do

                 call echprt(1,iot,m,ipan,ipat,jtpat,itnm,itmxpt,itmxgp)

               if( rtmme(m,k) .gt. 0.0 ) then

                  write(iot,'("     emax = ",1p1g13.5)') rtmme(m,k)

               end if

*-----------------------------------------------------------------------

                     do i = 1, itmst(m)

                        ilgt(i) = 0

                     end do

               do j = 1, itmln(m,k)

                  do im = 1, itmst(m)

                        nn = im * 2 - 1
                        jlgt(im) = 0

                  do jj = mltp(1+nn,itmli(m,k)/13+j),
     &                    mltp(1+nn+1,itmli(m,k)/13+j), 2

                        if( nint(slib(jj)) .eq. 1 .or.
     &                      nint(slib(jj)) .eq. 2 .or.
     &                      nint(slib(jj)) .eq. 3 ) then

                           jlgt(im) = jlgt(im) + 2

                        else if( nint(slib(jj)) .eq. 4 .or.
     &                           nint(slib(jj)) .eq. 8 ) then

                           write(ch13,'(1p1g13.5)') slib(jj+1)

                           do l = 1, 13
                              if( ch13(l:l) .ne. ' ' ) goto 410
                           end do
  410                      isi = l - 1
                           do l = 13, 1, -1
                              if( ch13(l:l) .ne. ' ' ) goto 420
                           end do
  420                      isf = 13 - l

                           jlgt(im) = jlgt(im) + 13 - isi - isf + 1

                        else if( nint(slib(jj)) .ge. 5 .and.
     &                           nint(slib(jj)) .le. 7 ) then

                           write(ch6,'(i6)') nint( slib(jj+1) )

                           do l = 1, 6
                              if( ch6(l:l) .ne. ' ' ) goto 430
                           end do
  430                      isi = l - 1

                           jlgt(im) = jlgt(im) + 6 - isi + 1

                        end if

                  end do
                  end do

                     do i = 1, itmst(m)

                        if( jlgt(i) .gt. ilgt(i) ) ilgt(i) = jlgt(i)

                     end do

               end do

                     do i = 1, itmst(m)

                        if( ilgt(i) .lt. 7 ) ilgt(i) = 7

                     end do

                  ch200(1:9) = '      mat'
                  ic = 9

               do i = 1, itmst(m)

                  isp = ilgt(i) - 7

                  do j = 1, isp
                     ch200(ic+j:ic+j) = ' '
                  end do
                     ic = ic + isp
                     ch200(ic+1:ic+4) = 'mset'
                     ic = ic + 4

                  write(ch200(ic+1:ic+1),'(i1)') itmnt(m,i)
                     ic = ic + 1

                  ch200(ic+1:ic+2) = '  '
                     ic = ic + 2

               end do

                  write(iot,'(200a1)') (ch200(j:j),j=1,ic)

*-----------------------------------------------------------------------

               do j = 1, itmln(m,k)

                     do ik = 1, 400
                        ch100(ik:ik) = ' '
                     end do

                     ch200(1:3) = '   '
                     ic = 3

                  if( mltp(1,itmli(m,k)/13+j) .gt. 0 ) then
                     write(ch200(ic+1:ic+6),'(i6)')
     &                  mltp(1,itmli(m,k)/13+j)

                   else

                     ch200(ic+1:ic+6) = '   all'

                   end if

                     ic = ic + 6

                  do im = 1, itmst(m)

                        nn = im * 2 - 1
                        id = 0

                  do jj = mltp(1+nn,itmli(m,k)/13+j),
     &                    mltp(1+nn+1,itmli(m,k)/13+j), 2

                        if( nint(slib(jj)) .eq. 1 ) then

                           ch100(id+1:id+2) = ' ('
                           id = id + 2

                        else if( nint(slib(jj)) .eq. 2 ) then

                           ch100(id+1:id+2) = ' )'
                           id = id + 2

                        else if( nint(slib(jj)) .eq. 3 ) then

                           ch100(id+1:id+2) = ' :'
                           id = id + 2

                        else if( nint(slib(jj)) .eq. 4 .or.
     &                           nint(slib(jj)) .eq. 8 ) then

                           write(ch13,'(1p1g13.5)') slib(jj+1)

                           do l = 1, 13
                              if( ch13(l:l) .ne. ' ' ) goto 440
                           end do
  440                      isi = l
                           do l = 13, 1, -1
                              if( ch13(l:l) .ne. ' ' ) goto 450
                           end do
  450                      isf = l

                           ch100(id+1:id+1) = ' '
                           id = id + 1
                           ch100(id+1:id+isf-isi+1) = ch13(isi:isf)
                           id = id+isf-isi+1

                        else if( nint(slib(jj)) .ge. 5 .and.
     &                           nint(slib(jj)) .le. 7 ) then

                           write(ch6,'(i6)') nint( slib(jj+1) )

                           do l = 1, 6
                              if( ch6(l:l) .ne. ' ' ) goto 460
                           end do
  460                      isi = l

                           ch100(id+1:id+1) = ' '
                           id = id + 1
                           ch100(id+1:id+6-isi+1) = ch6(isi:6)
                           id = id+6-isi+1

                        end if

                  end do

                           isd = ilgt(im) - id

                           do kk = 1, isd
                              ch200(ic+kk:ic+kk) = ' '
                           end do
                           ic = ic + isd

                           ch200(ic+1:ic+id) = ch100(1:id)
                           ic = ic + id

                  end do

                  write(iot,'(200a1)') (ch200(kk:kk),kk=1,ic)

               end do

*-----------------------------------------------------------------------

            end do

         end if

*-----------------------------------------------------------------------

               if( rtfac(m) .ne. 1.0d0 ) then

                  write(iot,'("   factor =",1p1g14.7,1x,1x,
     &            " # (D=1.0) normalization factor")')
     &            rtfac(m)

               end if

*-----------------------------------------------------------------------

               if( itmtn(m) .gt. 0 .and. itmsh(m) .eq. 3 ) then

                  write(iot,'("   volmat = ",i4,11x,
     &            " # (D=9) vol.corr., 0:no, number of scan / dir.")')
     &            itvm(m)

               end if

*-----------------------------------------------------------------------

               if( itanl(m) .gt. 0 ) then

                  write(iot,'("    angel = ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               end if

*-----------------------------------------------------------------------

               if( itsans(m) .gt. 0 ) then

                  write(iot,'("   sangel = ",i2)') itsans(m)
                  call write_sangel(iot,m,1)

               end if

*-----------------------------------------------------------------------

               if( itaxl(m) .gt. 0 ) then

                  write(iot,'("    x-txt = ",200a1)')
     &            ( itaxt(m)(i:i),i = 1, itaxl(m) )

               end if

*-----------------------------------------------------------------------

               if( itayl(m) .gt. 0 ) then

                  write(iot,'("    y-txt = ",200a1)')
     &            ( itayt(m)(i:i),i = 1, itayl(m) )

               end if

*-----------------------------------------------------------------------

               if( itazl(m) .gt. 0 ) then

                  write(iot,'("    z-txt = ",200a1)')
     &            ( itazt(m)(i:i),i = 1, itazl(m) )

               end if

*-----------------------------------------------------------------------

               if( itgsh(m) .gt. 0 ) then

                  write(iot,'("    gshow = ",i4,10x,
     &            "  # 0: no 1:bnd, 2:bnd+mat, ",
     &     "3:bnd+reg 4:bnd+lat 5:bmp style")')
     &            itgsh(m)

               end if

*-----------------------------------------------------------------------

               if( itglt(m) .ge. 0 ) then

                  write(iot,'("    gslat = ",i4,10x,
     &            "  # 0: no line for lat, 1: with,",
     &               " 2: no line between the same cell")')
     &            itglt(m)

               end if

*-----------------------------------------------------------------------

               if( itgsh(m) .gt. 0 .and. itger(m) .gt. 0 ) then

                  write(iot,'("    ginfo = ",i4,11x,
     &            " # (D=0) 0:no, 1:show in graph, 2:add file")')
     &            itger(m)

               end if

*-----------------------------------------------------------------------

               if( itrsh(m) .gt. 0 ) then

                  write(iot,'("    rshow = ",i4,10x,
     &            "  # 0: no 1:bnd, 2:bnd+mat, ",
     &     "3:bnd+reg 4:bnd+lat 5:bmp style")')
     &            itrsh(m)

                  call echmty(1,iot,'x',itxty(m),itxnm(m),
     &                        rtxdl(m),rtxmi(m),rtxma(m),itxrg(m)
     &                        ,abs(itxnm(m)+1),das_itxrg)
                  call echmty(1,iot,'y',ityty(m),itynm(m),
     &                        rtydl(m),rtymi(m),rtyma(m),ityrg(m)
     &                        ,abs(itynm(m)+1),das_ityrg)
                  call echmty(1,iot,'z',itzty(m),itznm(m),
     &                        rtzdl(m),rtzmi(m),rtzma(m),itzrg(m)
     &                        ,abs(itznm(m)+1),das_itzrg)

               end if

*-----------------------------------------------------------------------

               if( itgsh(m) .gt. 0 .or. itrsh(m) .gt. 0 ) then

                  write(iot,'("    resol = ",i4,10x,
     &            "  # (D=1) resolution of gshow or rshow")')
     &            itres(m)

                  write(iot,'("    width =",1p1g14.7,1x,1x,
     &            " # (D=0.5) width of lines for gshow or rshow")')
     &            rtwid(m)

               end if

*-----------------------------------------------------------------------

               if( iterl(m) .ne. 72 ) then

                  write(iot,'("   iechrl =",i5,10x,
     &            "  # (D=72) length reg of val or vol")')
     &            iterl(m)

               end if

*-----------------------------------------------------------------------

               if( iteps(m) .ne. 0 ) then

                  write(iot,'("   epsout =",i5,10x,
     &            "  # (D=0) generate eps file by ANGEL")')
     &            iteps(m)

               end if

*-----------------------------------------------------------------------

               if( itbmp(m) .ne. 0 ) then

                  write(iot,'("   bmpout =",i5,10x,
     &            "  # (D=0) generate bitmap file")')
     &            itbmp(m)

               end if

*-----------------------------------------------------------------------

               if( itvtk(m) .ne. 0 ) then

                  write(iot,'("   vtkout =",i5,10x,
     &            "  # (D=0) generate vtk file")')
     &            itvtk(m)

                  if ( itvtkfmt(m) .eq. 0 ) then

                     write(iot,'("   vtkfmt =",i5,10x,
     &               "  # (D=0) generate vtk ascii file")')
     &               itvtkfmt(m)

                  else

                     write(iot,'("   vtkfmt =",i5,10x,
     &               "  # (D=0) generate vtk binary file")')
     &               itvtkfmt(m)

                  end if

               end if

*-----------------------------------------------------------------------
*        check of counter
*-----------------------------------------------------------------------

            do i = 1, 3

               if( itcnt(i,m) .ne. 0 ) then

                  write(iot,'(" ctmin(",i1,") =",i5,10x,
     &            "  # (D=-9999) min. of counter ",i1)')
     &            i, itcnt(i*2+2,m), i

                  write(iot,'(" ctmax(",i1,") =",i5,10x,
     &            "  # (D= 9999) max. of counter ",i1)')
     &            i, itcnt(i*2+3,m), i

                  if( ncntc(i) .eq. 0 )
     &            write(iot,'("# Warning: there is no counter = ",
     &                        i1," definition in [counter].")') i

               end if

            end do

            do i = 1, 3
               if( itcnth(i,m) .ne. 0 ) then
                  write(iot,'(" chmin(",i1,") =",i5,10x,
     &            "  # (D=-9999) min. of history counter ",i1)')
     &            i, itcnth(i*2+2,m), i
                  write(iot,'(" chmax(",i1,") =",i5,10x,
     &            "  # (D= 9999) max. of history counter ",i1)')
     &            i, itcnth(i*2+3,m), i
               end if
            end do

*-----------------------------------------------------------------------

               if( rtstd(m) .gt. 0.d0 ) then

                  write(iot,'("   stdcut =",1p1g14.7,1x,1x,
     &            " # (D=-1) threshold value of STD cutoff")')
     &            rtstd(m)

               end if

*-----------------------------------------------------------------------

         if( itextstat(m) .gt. 0 ) then

            write(iot,'(" iextstat =",i5,10x,
     &      "  # (D=0) option for extended statistical indicators")')
     &           itextstat(m)

         end if

*-----------------------------------------------------------------------

         if( itprodenchk(m) .gt. 0 ) then

            write(iot,'(" prodenmn =",1p1g14.7,1x,1x,
     &     " # (D=1e-2) lower limit of probability density function")')
     &           tprodenmn(m)

            write(iot,'(" prodenmx =",1p1g14.7,1x,1x,
     &     " # (D=1e+2) upper limit of probability density function")')
     &           tprodenmx(m)

            write(iot,'(" nbproden =",i5,10x,
     &   "  # (D=100) number of bins of probability density function")')
     &           nbtproden(m)

         end if

*-----------------------------------------------------------------------
*     anatally
*-----------------------------------------------------------------------
                  yen  = char(92)

         if( itism(m) .gt. 0 ) then

                  write(iot,'(/"anatally start")')

            do i = 1, 10
            if( itkst(i,m) .ge. -1 ) then
                     ch73(1:11) = '   '//schan(i)//' = '
                     ilst = 11
               if( itkst(i,m) .eq. -1 ) then
                     ch73(ilst+1:ilst+4) = ' all'
                     ilst = ilst+4
                     write(iot,'(73(a1))') (ch73(j:j),j=1,ilst)
               else if( itkst(i,m) .ne. 0 ) then
                  do k = 1, itist(i,m)
                        isdnm = idas(itstt(m)+itjst(i,m)+k)
                     if( isdnm .lt. 10 ) then
                        write(ch73(ilst+1:ilst+3),'(i3)') isdnm
                        ilst = ilst + 3
                     else if( isdnm .lt. 100 ) then
                        write(ch73(ilst+1:ilst+4),'(i4)') isdnm
                        ilst = ilst + 4
                     else if( isdnm .lt. 1000 ) then
                        write(ch73(ilst+1:ilst+5),'(i5)') isdnm
                        ilst = ilst + 5
                     else if( isdnm .lt. 10000 ) then
                        write(ch73(ilst+1:ilst+6),'(i6)') isdnm
                        ilst = ilst + 6
                     end if
                  end do

                  if( ilst .le. 70 ) then
                           write(iot,'(73(a1))') (ch73(j:j),j=1,ilst)
                  else if( ilst .le. 130 ) then
                     do l = 65, 73
                        if( ch73(l:l) .eq. ' ' ) goto 111
                     end do
  111                   jlst = l
                        write(iot,'(74(a1))') (ch73(j:j),j=1,jlst), yen
                        write(iot,'(14x,74(a1))')
     &                                        (ch73(j:j),j=jlst+1,ilst)
                  else if( ilst .le. 190 ) then
                     do l = 65, 73
                        if( ch73(l:l) .eq. ' ' ) goto 112
                     end do
  112                   jlst = l
                        write(iot,'(74(a1))') (ch73(j:j),j=1,jlst), yen
                     do l = 125, 133
                        if( ch73(l:l) .eq. ' ' ) goto 113
                     end do
  113                   klst = l
                        write(iot,'(14x,74(a1))')
     &                                   (ch73(j:j),j=jlst+1,klst), yen
                        write(iot,'(14x,74(a1))')
     &                                        (ch73(j:j),j=klst+1,ilst)
                  end if
               end if

            end if
            end do

                  write(iot,'("anatally end")')

         end if

*-----------------------------------------------------------------------
               if( itfoam(m) .ne. 0 ) then

                  write(iot,'("  foamout =",i5,10x,
     &            "  # (D=0) generate OpenFOAM data")')
     &            itfoam(m)

               end if

*-----------------------------------------------------------------------
               if( itmxang(m) .ne. 0 ) then

                  write(iot,'(" maxangel =",i5,10x,
     &            "  # (D=0) number of output particles")')
     &            itmxang(m)

               end if

*-----------------------------------------------------------------------

*-----------------------------------------------------------------------
*        sumtally subsection
*-----------------------------------------------------------------------
            if ( icntl==13 .and. nsumtalRead(m)==1 ) then

               call sumtal_echo(iot,m)

            end if

*-----------------------------------------------------------------------

            write(iot,'("#    used :",
     &                  "        main (  %)",
     &                  "        temp (  %)",
     &                  "       total (  %)"/
     &                  "#  memory :",
     &                  3(i12," (",i3,")"))')
     &                  itsmn(m), nint(dble(itsmn(m))/dble(mdas)*100),
     &                  itstm(m), nint(dble(itstm(m))/dble(mdas)*100),
     &                  itsmn(m)+itstm(m),
     &                  nint(dble(itsmn(m)+itstm(m))/dble(mdas)*100)

*-----------------------------------------------------------------------
cfrtati 2024/02/13
            if( itmlp(m).gt.0 ) then
              call displaymtinfo(iot,m)
            end if

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine tadjech(iot,m,iaxs,icax)
*                                                                      *
*       input echo of the adjoint tally                                *
*       last modified by K.Niita on 2016/01/22                         *
*                                                                      *
************************************************************************

      use sumtallymod
      use sangelmod, only: itsans
      use partmod, only: itmxpt, itmxgp, itpan, itpat, jtpat, ipat, ipan ! frtati 2021/10/05
      use moddas_mesh
      use moddas_region
      use moddas_tally

*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'

*-----------------------------------------------------------------------

      common /tall00/ itnm, ital(itlmax), itals(itlmax), italm(itlmax)
      common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                itout(itlmax), ittwo(itlmax)

      common /tall02/ itrgn(itlmax), itrgm(itlmax), itreg(itlmax)
      common /tall24/ itrcm(itlmax), itrcg(itlmax)
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

      common /tall16/ itnfr(itlmax), itnfz(itlmax), itnfn(itlmax)

      common /tall14/ ittll(itlmax), ittle(itlmax)
      character ittle*80

      common /tall15/ itanl(itlmax), itang(itlmax)
      character itang*200

      common /tall27/ itaxl(itlmax), itaxt(itlmax),
     &                itayl(itlmax), itayt(itlmax),
     &                itazl(itlmax), itazt(itlmax)
      character itaxt*200, itayt*200, itazt*200

      common /tall19/ itsmn(itlmax), itstm(itlmax)

      common /tall23/ itaty(itlmax), itanm(itlmax), itarg(itlmax),
     &                rtami(itlmax), rtama(itlmax), rtadl(itlmax)


      common /tall20/ ittty(itlmax), ittnm(itlmax), ittrg(itlmax),
     &                rttmi(itlmax), rttma(itlmax), rttdl(itlmax)
      common /tall21/ rtfac(itlmax)
      common /tall28/ itgsh(itlmax)
      common /tall29/ itrsh(itlmax)
      common /tall31/ iterl(itlmax)
      common /tall22/ itmcn(itlmax), itmtn(itlmax), itmtt(itlmax)
      common /tall34/ itvm(itlmax)
      common /tall35/ iteps(itlmax)
      common /tall37/ itcnt(9,itlmax)
      common /cntreg/ ncntc(3), ncreg(3), incrc(3), incrt(3), kcont(3)
      common /cntpat/ icpan(3), icpat(3,20,2)

      common /tall38/ itres(itlmax)
      common /tall39/ rtwid(itlmax)

      common /tall45/ itmlp(itlmax), itmln(itlmax,6), itmst(itlmax),
     &                itmli(itlmax,6), rtmme(itlmax,6), itmnt(itlmax,6),
     &                itmpn(itlmax,6), itmpt(itlmax,6,6,2)
      common /tall46/ itmtr(itlmax,4), rtmtr(itlmax,13)
      common /tall49/ itglt(itlmax)

      common /tall60/ itger(itlmax)

      common /tall66/ itadm(itlmax), rtade(itlmax), rtadw(itlmax)

      common /tall67/ itstd(itlmax), rtstd(itlmax)

      common /tall82/ itcnth(9,itlmax)
cfrtati 2023/12/07
      common /tall84/ itmto(itlmax,6,4)

      dimension idas(1)
      equivalence ( das, idas )


cfrtati 2021/10/05 moved to partmod

      dimension ilgt(6), jlgt(6)
      character ch13*13, ch6*6
      character ch200*400
      character ch100*400

*-----------------------------------------------------------------------

      common /inpec/ ititl, ipara, ibody, iregn, llarr, itby, itar

*-----------------------------------------------------------------------
       integer   icntl
       integer   inucr

       common /tcntl/  icntl, inucr
!-----------------------------------------------------------------------

      character asfil*100

      character aname(13)*3
      data aname / 'eng','reg','  x','  y','  z','  r',
     &             ' xy',' yz',' xz',' rz','  t','rad','deg'/

      dimension ian(6)
      dimension jmat(20)

*-----------------------------------------------------------------------
      common /tall79/ itallech
      if( itallech .eq. 0 ) return

*-----------------------------------------------------------------------
*     input echo
*-----------------------------------------------------------------------

                  write(iot,'("[ T-Adjoint ]")')

                  write(iot,'("    title = ",80a1)')
     &            ( ittle(m)(i:i), i = 1, ittll(m) )

            if( itmsh(m) .eq. 1 ) then

                  write(iot,'("     mesh =  reg",11x,
     &            " # mesh type is region-wise")')

                     idas1 = mmmax
                     idas2 = ( idas1 + itrgn(m) - 1 ) * 2 + 1
                     idas3 = idas2 + itrgn(m)

                  call echrg(1,iterl(m),0,1,iot,itrcm(m),
     &                       idas_itrcg(itrcg(m)),
     &                       itrgn(m),itrgm(m),idas_itreg(itreg(m)),
     &                       das(idas1),idas(idas2),
     &                       itrnv(m),idas(itriv(m)),das(itrrv(m)),
     &                       idas3)

            else if( itmsh(m) .eq. 2 ) then

                  write(iot,'("     mesh =  r-z",11x,
     &            " # mesh type is r-z scoring mesh")')

                  write(iot,'("       x0 = ",1p1g14.7,1x,
     &            " # center x-position of r-z mesh")') rtrx0(m)

                  write(iot,'("       y0 = ",1p1g14.7,1x,
     &            " # center y-position of r-z mesh")') rtry0(m)

                  call echmty(0,iot,'r',itrty(m),itrnm(m),
     &                        rtrdl(m),rtrmi(m),rtrma(m),itrrg(m)
     &                        ,abs(itrnm(m)+1),das_itrrg)
                  call echmty(0,iot,'z',itzty(m),itznm(m),
     &                        rtzdl(m),rtzmi(m),rtzma(m),itzrg(m)
     &                        ,abs(itznm(m)+1),das_itzrg)

*-----------------------------------------------------------------------
               if( itaty(m) .ne. 0 ) then

                  call echmty(0,iot,'a',itaty(m),itanm(m),
     &                        rtadl(m),rtami(m),rtama(m),itarg(m)
     &                        ,abs(itanm(m)+1),das_itarg)

               end if
*-----------------------------------------------------------------------

            else if( itmsh(m) .eq. 3 ) then

                  write(iot,'("     mesh =  xyz",11x,
     &            " # mesh type is xyz scoring mesh")')

                  call echmty(0,iot,'x',itxty(m),itxnm(m),
     &                        rtxdl(m),rtxmi(m),rtxma(m),itxrg(m)
     &                        ,abs(itxnm(m)+1),das_itxrg)
                  call echmty(0,iot,'y',ityty(m),itynm(m),
     &                        rtydl(m),rtymi(m),rtyma(m),ityrg(m)
     &                        ,abs(itynm(m)+1),das_ityrg)
                  call echmty(0,iot,'z',itzty(m),itznm(m),
     &                        rtzdl(m),rtzmi(m),rtzma(m),itzrg(m)
     &                        ,abs(itznm(m)+1),das_itzrg)

            end if

*-----------------------------------------------------------------------

            if( itmtr(m,1) .gt. 0 ) then

               if( itmtr(m,1) .eq. 1 ) then

                     write(iot,'("     trcl = ",i6,9x,
     &               " # transform ID")') itmtr(m,3)

               else

                  if( itmtr(m,2) .eq. 0 ) then

                     write(iot,'("     trcl = ",
     &                           3(1p1e15.7),
     &                           3(/12x,3(1p1e15.7)),i5)')
     &               ( rtmtr(m,j), j = 1, 12 ),
     &               nint( rtmtr(m,13) )

                  else

                     write(iot,'("    *trcl = ",
     &                           3(1p1e15.7),
     &                           3(/12x,3(1p1e15.7)),i5)')
     &               ( rtmtr(m,j), j = 1, 12 ),
     &               nint( rtmtr(m,13) )

                  end if

               end if

            end if

*-----------------------------------------------------------------------

                  call echmty(0,iot,'e',itety(m),itenm(m),
     &                        rtedl(m),rtemi(m),rtema(m),iterg(m)
     &                        ,abs(itenm(m)+1),das_iterg)

*-----------------------------------------------------------------------

               if( ittty(m) .ne. 0 ) then

                  call echmty(0,iot,'t',ittty(m),ittnm(m),
     &                        rttdl(m),rttmi(m),rttma(m),ittrg(m)
     &                        ,abs(ittnm(m)+1),das_ittrg)

               end if

*-----------------------------------------------------------------------

               if( itunt(m) .eq. 1 ) then

                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [1/cm^2/source]")') itunt(m)

               else if( itunt(m) .eq. 2 ) then

                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [1/cm^2/MeV/source]")') itunt(m)

               else if( itunt(m) .eq. 3 ) then

                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [1/cm^2/Lethargy/source]")') itunt(m)

               else if( itunt(m) .eq. 4 ) then

                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [1/cm^2/source] with vol=1cm^3")')
     &               itunt(m)

               else if( itunt(m) .eq. 11 ) then

                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [1/cm^2/nsec/source]")') itunt(m)

               else if( itunt(m) .eq. 12 ) then

                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [1/cm^2/MeV/nsec/source]")') itunt(m)

               else if( itunt(m) .eq. 13 ) then

                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [1/cm^2/Lethargy/nsec/source]")')
     &               itunt(m)

               else if( itunt(m) .eq. 14 ) then

                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [1/cm^2/nsec/source] with vol=1cm^3")')
     &               itunt(m)

               end if

*-----------------------------------------------------------------------

                  asfil = '  # file name of output for the above axis'

            if( icax .eq. 0 ) then

                  icn = itaxn(m)

               do i = 1, icn

                  ian(i) = i

               end  do

            else

                  icn = 1
                  ian(1) = iaxs

            end if

*-----------------------------------------------------------------------

               if( itmtn(m) .eq. 0 ) then

                  write(iot,'(" material =  all"11x,
     &            " # (D=all) number of specific material")')

               else if( itmtn(m) .gt. 0 ) then

                  write(iot,'(" material =",i5,11x,
     &            " # (D=all) number of specific material")')
     &            itmtn(m) * itmcn(m)

                     k = 0

                  do j = 1, itmtn(m)

                     k = k + 1

                     jmat(k) = ismte( itmtt(m) + j - 1 )

                     if( k / 10 * 10 .eq. k .or. j .eq. itmtn(m) ) then

                        write(iot,'(11x,10(i6))') ( jmat(i),i = 1, k )

                        k = 0

                     end if

                  end do

               end if

*-----------------------------------------------------------------------

               idtyp = 0

            do i = 1, icn

               if( itaxs(m,ian(i)) .ge. 7 .and.
     &             itaxs(m,ian(i)) .le. 10 )  idtyp = ittwo(m)

            end do

            if( idtyp .ne. 0 ) then

               write(iot,'("  2D-type = ",i4,11x,
     &         " # 1:Cont, 2:Clust, 3:Color, 4:xyz, 5:mat,",
     &         " 6:Clust+Cont, 7:Col+Cont")') idtyp

            end if

*-----------------------------------------------------------------------

            do i = 1, icn

                  iax = ian(i)

                  write(iot,'("     axis =  ",a3,11x,
     &            " # axis of output")')
     &            aname( itaxs(m,iax) )

            if ( itfll(m,iax) .gt. 0 ) then

                  msfile = max( 14, itfll(m,iax) )

                  write(iot,'("     file = ",100a1)')
     &                  ( ctfln(m,iax)(j:j), j = 1, msfile ),
     &                  ( asfil(j:j), j = 1, 42 )

            end if

            end do

*-----------------------------------------------------------------------
*        multiplier
*-----------------------------------------------------------------------

         if( itmlp(m) .gt. 0 ) then

            do k = 1, itmlp(m)

*-----------------------------------------------------------------------

               if( itmln(m,k) .eq. 1 .and.
     &             mltp(1,itmli(m,k)/13+1) .eq. -1 ) then

                  imall = 0

                  write(iot,'("  multiplier = all",9x,
     &            " # number of material group")')

               else

                  imall = itmln(m,k)

                  write(iot,'("  multiplier = ",i3,9x,
     &            " # number of material group")') imall

               end if

*-----------------------------------------------------------------------

                     ipan(m) = itmpn(m,k)

                  do j = 1, itmpn(m,k)

                     ipat(m,j,1) = itmpt(m,k,j,1)
                     ipat(m,j,2) = itmpt(m,k,j,2)

                  end do

                 call echprt(1,iot,m,ipan,ipat,jtpat,itnm,itmxpt,itmxgp)

               if( rtmme(m,k) .gt. 0.0 ) then

                  write(iot,'("     emax = ",1p1g13.5)') rtmme(m,k)

               end if

*-----------------------------------------------------------------------

                     do i = 1, itmst(m)

                        ilgt(i) = 0

                     end do

               do j = 1, itmln(m,k)

                  do im = 1, itmst(m)

                        nn = im * 2 - 1
                        jlgt(im) = 0

                  do jj = mltp(1+nn,itmli(m,k)/13+j),
     &                    mltp(1+nn+1,itmli(m,k)/13+j), 2

                        if( nint(slib(jj)) .eq. 1 .or.
     &                      nint(slib(jj)) .eq. 2 .or.
     &                      nint(slib(jj)) .eq. 3 ) then

                           jlgt(im) = jlgt(im) + 2

                        else if( nint(slib(jj)) .eq. 4 .or.
     &                           nint(slib(jj)) .eq. 8 ) then

                           write(ch13,'(1p1g13.5)') slib(jj+1)

                           do l = 1, 13
                              if( ch13(l:l) .ne. ' ' ) goto 410
                           end do
  410                      isi = l - 1
                           do l = 13, 1, -1
                              if( ch13(l:l) .ne. ' ' ) goto 420
                           end do
  420                      isf = 13 - l

                           jlgt(im) = jlgt(im) + 13 - isi - isf + 1

                        else if( nint(slib(jj)) .ge. 5 .and.
     &                           nint(slib(jj)) .le. 7 ) then

                           write(ch6,'(i6)') nint( slib(jj+1) )

                           do l = 1, 6
                              if( ch6(l:l) .ne. ' ' ) goto 430
                           end do
  430                      isi = l - 1

                           jlgt(im) = jlgt(im) + 6 - isi + 1

                        end if

                  end do
                  end do

                     do i = 1, itmst(m)

                        if( jlgt(i) .gt. ilgt(i) ) ilgt(i) = jlgt(i)

                     end do

               end do

                     do i = 1, itmst(m)

                        if( ilgt(i) .lt. 7 ) ilgt(i) = 7

                     end do

                  ch200(1:9) = '      mat'
                  ic = 9

               do i = 1, itmst(m)

                  isp = ilgt(i) - 7

                  do j = 1, isp
                     ch200(ic+j:ic+j) = ' '
                  end do
                     ic = ic + isp
                     ch200(ic+1:ic+4) = 'mset'
                     ic = ic + 4

                  write(ch200(ic+1:ic+1),'(i1)') itmnt(m,i)
                     ic = ic + 1

                  ch200(ic+1:ic+2) = '  '
                     ic = ic + 2

               end do

                  write(iot,'(200a1)') (ch200(j:j),j=1,ic)

*-----------------------------------------------------------------------

               do j = 1, itmln(m,k)

                     do ik = 1, 400
                        ch100(ik:ik) = ' '
                     end do

                     ch200(1:3) = '   '
                     ic = 3

                  if( mltp(1,itmli(m,k)/13+j) .gt. 0 ) then
                     write(ch200(ic+1:ic+6),'(i6)')
     &                  mltp(1,itmli(m,k)/13+j)

                   else

                     ch200(ic+1:ic+6) = '   all'

                   end if

                     ic = ic + 6

                  do im = 1, itmst(m)

                        nn = im * 2 - 1
                        id = 0

                  do jj = mltp(1+nn,itmli(m,k)/13+j),
     &                    mltp(1+nn+1,itmli(m,k)/13+j), 2

                        if( nint(slib(jj)) .eq. 1 ) then

                           ch100(id+1:id+2) = ' ('
                           id = id + 2

                        else if( nint(slib(jj)) .eq. 2 ) then

                           ch100(id+1:id+2) = ' )'
                           id = id + 2

                        else if( nint(slib(jj)) .eq. 3 ) then

                           ch100(id+1:id+2) = ' :'
                           id = id + 2

                        else if( nint(slib(jj)) .eq. 4 .or.
     &                           nint(slib(jj)) .eq. 8 ) then

                           write(ch13,'(1p1g13.5)') slib(jj+1)

                           do l = 1, 13
                              if( ch13(l:l) .ne. ' ' ) goto 440
                           end do
  440                      isi = l
                           do l = 13, 1, -1
                              if( ch13(l:l) .ne. ' ' ) goto 450
                           end do
  450                      isf = l

                           ch100(id+1:id+1) = ' '
                           id = id + 1
                           ch100(id+1:id+isf-isi+1) = ch13(isi:isf)
                           id = id+isf-isi+1

                        else if( nint(slib(jj)) .ge. 5 .and.
     &                           nint(slib(jj)) .le. 7 ) then

                           write(ch6,'(i6)') nint( slib(jj+1) )

                           do l = 1, 6
                              if( ch6(l:l) .ne. ' ' ) goto 460
                           end do
  460                      isi = l

                           ch100(id+1:id+1) = ' '
                           id = id + 1
                           ch100(id+1:id+6-isi+1) = ch6(isi:6)
                           id = id+6-isi+1

                        end if

                  end do

                           isd = ilgt(im) - id

                           do kk = 1, isd
                              ch200(ic+kk:ic+kk) = ' '
                           end do
                           ic = ic + isd

                           ch200(ic+1:ic+id) = ch100(1:id)
                           ic = ic + id

                  end do

                  write(iot,'(200a1)') (ch200(kk:kk),kk=1,ic)

               end do

*-----------------------------------------------------------------------

            end do

         end if

*-----------------------------------------------------------------------

               if( rtfac(m) .ne. 1.0d0 ) then

                  write(iot,'("   factor =",1p1g14.7,1x,1x,
     &            " # (D=1.0) normalization factor")')
     &            rtfac(m)

               end if

*-----------------------------------------------------------------------

               if( itmtn(m) .gt. 0 .and. itmsh(m) .eq. 3 ) then

                  write(iot,'("   volmat = ",i4,11x,
     &            " # (D=9) vol.corr., 0:no, number of scan / dir.")')
     &            itvm(m)

               end if

*-----------------------------------------------------------------------

               if( itanl(m) .gt. 0 ) then

                  write(iot,'("    angel = ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               end if

*-----------------------------------------------------------------------

               if( itsans(m) .gt. 0 ) then

                  write(iot,'("   sangel = ",i2)') itsans(m)
                  call write_sangel(iot,m,1)

               end if

*-----------------------------------------------------------------------

               if( itaxl(m) .gt. 0 ) then

                  write(iot,'("    x-txt = ",200a1)')
     &            ( itaxt(m)(i:i),i = 1, itaxl(m) )

               end if

*-----------------------------------------------------------------------

               if( itayl(m) .gt. 0 ) then

                  write(iot,'("    y-txt = ",200a1)')
     &            ( itayt(m)(i:i),i = 1, itayl(m) )

               end if

*-----------------------------------------------------------------------

               if( itazl(m) .gt. 0 ) then

                  write(iot,'("    z-txt = ",200a1)')
     &            ( itazt(m)(i:i),i = 1, itazl(m) )

               end if

*-----------------------------------------------------------------------

               if( itgsh(m) .gt. 0 ) then

                  write(iot,'("    gshow = ",i4,10x,
     &            "  # 0: no 1:bnd, 2:bnd+mat, ",
     &     "3:bnd+reg 4:bnd+lat 5:bmp style")')
     &            itgsh(m)

               end if

*-----------------------------------------------------------------------

               if( itglt(m) .ge. 0 ) then

                  write(iot,'("    gslat = ",i4,10x,
     &            "  # 0: no line for lat, 1: with,",
     &               " 2: no line between the same cell")')
     &            itglt(m)

               end if

*-----------------------------------------------------------------------

               if( itgsh(m) .gt. 0 .and. itger(m) .gt. 0 ) then

                  write(iot,'("    ginfo = ",i4,11x,
     &            " # (D=0) 0:no, 1:show in graph, 2:add file")')
     &            itger(m)

               end if

*-----------------------------------------------------------------------

               if( itrsh(m) .gt. 0 ) then

                  write(iot,'("    rshow = ",i4,10x,
     &            "  # 0: no 1:bnd, 2:bnd+mat, ",
     &     "3:bnd+reg 4:bnd+lat 5:bmp style")')
     &            itrsh(m)

                  call echmty(1,iot,'x',itxty(m),itxnm(m),
     &                        rtxdl(m),rtxmi(m),rtxma(m),itxrg(m)
     &                        ,abs(itxnm(m)+1),das_itxrg)
                  call echmty(1,iot,'y',ityty(m),itynm(m),
     &                        rtydl(m),rtymi(m),rtyma(m),ityrg(m)
     &                        ,abs(itynm(m)+1),das_ityrg)
                  call echmty(1,iot,'z',itzty(m),itznm(m),
     &                        rtzdl(m),rtzmi(m),rtzma(m),itzrg(m)
     &                        ,abs(itznm(m)+1),das_itzrg)

               end if

*-----------------------------------------------------------------------

               if( itgsh(m) .gt. 0 .or. itrsh(m) .gt. 0 ) then

                  write(iot,'("    resol = ",i4,10x,
     &            "  # (D=1) resolution of gshow or rshow")')
     &            itres(m)

                  write(iot,'("    width =",1p1g14.7,1x,1x,
     &            " # (D=0.5) width of lines for gshow or rshow")')
     &            rtwid(m)

               end if

*-----------------------------------------------------------------------

               if( iterl(m) .ne. 72 ) then

                  write(iot,'("   iechrl =",i5,10x,
     &            "  # (D=72) length reg of val or vol")')
     &            iterl(m)

               end if

*-----------------------------------------------------------------------

               if( iteps(m) .ne. 0 ) then

                  write(iot,'("   epsout =",i5,10x,
     &            "  # (D=0) generate eps file by ANGEL")')
     &            iteps(m)

               end if

*-----------------------------------------------------------------------
*        check of counter
*-----------------------------------------------------------------------

            do i = 1, 3

               if( itcnt(i,m) .ne. 0 ) then

                  write(iot,'(" ctmin(",i1,") =",i5,10x,
     &            "  # (D=-9999) min. of counter ",i1)')
     &            i, itcnt(i*2+2,m), i

                  write(iot,'(" ctmax(",i1,") =",i5,10x,
     &            "  # (D= 9999) max. of counter ",i1)')
     &            i, itcnt(i*2+3,m), i

                  if( ncntc(i) .eq. 0 )
     &            write(iot,'("# Warning: there is no counter = ",
     &                        i1," definition in [counter].")') i

               end if

            end do

            do i = 1, 3
               if( itcnth(i,m) .ne. 0 ) then
                  write(iot,'(" chmin(",i1,") =",i5,10x,
     &            "  # (D=-9999) min. of history counter ",i1)')
     &            i, itcnth(i*2+2,m), i
                  write(iot,'(" chmax(",i1,") =",i5,10x,
     &            "  # (D= 9999) max. of history counter ",i1)')
     &            i, itcnth(i*2+3,m), i
               end if
            end do

*-----------------------------------------------------------------------
*      e-source and e-width
*-----------------------------------------------------------------------

                  write(iot,'("   e-smin =",1p1g14.7,1x,1x,
     &            " # (D=0.58) adjoint min source energy")')
     &            rtade(m)

                  write(iot,'("   e-smax =",1p1g14.7,1x,1x,
     &            " # (D=0.62) adjoint max source energy")')
     &            rtadw(m)

                  write(iot,'(" m-source =",i5,10x,
     &            "  # (D=1) number of multi-source")')
     &            itadm(m)

*-----------------------------------------------------------------------

               if( rtstd(m) .gt. 0.d0 ) then

                  write(iot,'("   stdcut =",1p1g14.7,1x,1x,
     &            " # (D=-1) threshold value of STD cutoff")')
     &            rtstd(m)

               end if

*-----------------------------------------------------------------------
*        sumtally subsection
*-----------------------------------------------------------------------
            if ( icntl==13 .and. nsumtalRead(m)==1 ) then

               call sumtal_echo(iot,m)

            end if

*-----------------------------------------------------------------------

            write(iot,'("#    used :",
     &                  "        main (  %)",
     &                  "        temp (  %)",
     &                  "       total (  %)"/
     &                  "#  memory :",
     &                  3(i12," (",i3,")"))')
     &                  itsmn(m), nint(dble(itsmn(m))/dble(mdas)*100),
     &                  itstm(m), nint(dble(itstm(m))/dble(mdas)*100),
     &                  itsmn(m)+itstm(m),
     &                  nint(dble(itsmn(m)+itstm(m))/dble(mdas)*100)

*-----------------------------------------------------------------------
cfrtati 2024/02/13
            if( itmlp(m).gt.0 ) then
              call displaymtinfo(iot,m)
            end if

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine tcrsech(iot,m,iaxs,icax)
*                                                                      *
*       input echo of the surface crossing flux tally                  *
*       last modified by K.Niita on 2005/11/23                         *
*                                                                      *
************************************************************************

      use sumtallymod
      use sangelmod, only: itsans
      use partmod ! frtati 2021/10/05
      use moddas_mesh
      use moddas_region
      use moddas_tally
*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'

*-----------------------------------------------------------------------

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


      common /tall23/ itaty(itlmax), itanm(itlmax), itarg(itlmax),
     &                rtami(itlmax), rtama(itlmax), rtadl(itlmax)

      common /tall08/ rtrx0(itlmax), rtry0(itlmax)

      common /tall10/ itaxn(itlmax), itaxs(itlmax,6)
      common /tall11/ itfln(itlmax), itfll(itlmax,6), ctfln(itlmax,6)
      character ctfln*100

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

      common /tall19/ itsmn(itlmax), itstm(itlmax)
      common /tall20/ ittty(itlmax), ittnm(itlmax), ittrg(itlmax),
     &                rttmi(itlmax), rttma(itlmax), rttdl(itlmax)
      common /tall21/ rtfac(itlmax)
      common /tall28/ itgsh(itlmax)
      common /tall35/ iteps(itlmax)
      common /tall37/ itcnt(9,itlmax)
      common /cntreg/ ncntc(3), ncreg(3), incrc(3), incrt(3), kcont(3)
      common /cntpat/ icpan(3), icpat(3,20,2)

      common /tall38/ itres(itlmax)
      common /tall39/ rtwid(itlmax)
      common /tall46/ itmtr(itlmax,4), rtmtr(itlmax,13)
      common /tall48/ itmdp(itlmax,0:30), itmdf(itlmax)
      common /tall49/ itglt(itlmax)

      common /tall60/ itger(itlmax)

*-----------------------------------------------------------------------
      common /tall45/ itmlp(itlmax), itmln(itlmax,6), itmst(itlmax),
     &                itmli(itlmax,6), rtmme(itlmax,6), itmnt(itlmax,6),
     &                itmpn(itlmax,6), itmpt(itlmax,6,6,2)

      common /tall67/ itstd(itlmax), rtstd(itlmax)

      common /tall73/ itenclo(itlmax), itangform(itlmax)

      common /tall82/ itcnth(9,itlmax)
cfrtati 2023/12/07
      common /tall84/ itmto(itlmax,6,4)

ccse 2022/08/31 LET parameter
      common /tall50/ itlmt(itlmax), ite2l(itlmax)   ! CCSE 2022.08.31


cfrtati 2021/10/05 moved to partmod

      dimension ilgt(6), jlgt(6)
      character ch13*13, ch6*6
      character ch200*400
      character ch100*400

*-----------------------------------------------------------------------

      dimension idas(100)
      equivalence ( das, idas )

*-----------------------------------------------------------------------

      common /inpec/ ititl, ipara, ibody, iregn, llarr, itby, itar

*-----------------------------------------------------------------------
       integer   icntl
       integer   inucr

       common /tcntl/  icntl, inucr
!-----------------------------------------------------------------------

      character asfil*100
      character bsfil*100
      character csfil*100
      character dsfil*100

ccse 2022/08/31, add LET axis, 13 -> 14
cABE 2018/10/18, add yz, xz, rz axis, 10 -> 13
      character aname(14)*3
      data aname / 'eng','reg','  x','  y','  z','  r',
     &             ' xy','cos','  t','the',' yz',' xz',
     &             ' rz','let'/

      dimension ian(6)

      character dum1*10000
      character dum2*10000

      character dmpc(30)*4
      data dmpc / '  kf','   x','   y','   z','   u','   v','   w',
     &            '   e','  wt','  tm','  c1','  c2','  c3',
     &            '  sx','  sy','  sz','  n0','  nc','  nb','  no',
     &            '    ','    ','    ','    ','    ','    ',
     &            '    ','    ','    ','    '/

      integer           iMeVperu
      common /cMeVperu/ iMeVperu

*-----------------------------------------------------------------------
      common /tall79/ itallech
      if( itallech .eq. 0 ) return

*-----------------------------------------------------------------------
*     input echo
*-----------------------------------------------------------------------

            write(iot,'("[ T-Cross ]")')

                  write(iot,'("    title = ",80a1)')
     &            ( ittle(m)(i:i), i = 1, ittll(m) )

*-----------------------------------------------------------------------

            if( itmsh(m) .eq. 1 ) then

                  write(iot,'("     mesh =  reg",11x,
     &            " # mesh type is region-wise")')

                  write(iot,'("      reg = ",i4,11x,
     &            " # number of crossing regions")') itrcn(m)

*-----------------------------------------------------------------------

                  igm = ( mmmax - 1 ) * 2 + 2

                  lgmax1 = 0
                  lgmax2 = 0

                  idsm = itrcc(m)
                  jdsm = -1

            do j = 1, itrcn(m)

                  jdsm = jdsm + 1
                  ntrn = idas_itrcc(idsm+jdsm)
                  jdsm = jdsm + 1
                  mtrn = idas_itrcc(idsm+jdsm)
                  kssm = idsm + jdsm + 1
                  jdsm = jdsm + mtrn

                  call echrg2(mtrn,idas_itrcc(kssm),dum1,lng1,icmb,igm)

                  if( lng1 .gt. lgmax1 ) lgmax1 = lng1

                  jdsm = jdsm + 1
                  ntrn = idas_itrcc(idsm+jdsm)
                  jdsm = jdsm + 1
                  mtrn = idas_itrcc(idsm+jdsm)
                  kssm = idsm + jdsm + 1
                  jdsm = jdsm + mtrn

                  call echrg2(mtrn,idas_itrcc(kssm),dum1,lng1,icmb,igm)

                  if( lng1 .gt. lgmax2 ) lgmax2 = lng1

            end do

                  dum2(1:9) = '      non'
                  l = 9

C S.H. revised for adding r-from and r-to in [t-cross] (2017.9.12)
                  dum2(l+1:l+11) = '     r-from'
                  l = l + 11
                  lblk = lgmax1
                do i = 1, lblk
                  dum2(i+l:i+l) = ' '
                end do
                  l = l + lblk

C S.H. revised for adding r-from and r-to in [t-cross] (2017.9.12)
                  dum2(l+1:l+8) =
     &                       ' r-to   '
                  l = l + 8
                  lblk = lgmax2
                do i = 1, lblk
                  dum2(i+l:i+l) = ' '
                end do
                  l = l + lblk

C S.H. revised for adding r-from and r-to in [t-cross] (2017.9.27)
                  dum2(l+1:l+5) =
     &                       ' area'
                  l = l + 5

                  write(iot,'(600a1)') (dum2(i:i),i=1,l)

*-----------------------------------------------------------------------

                  idsm = itrcc(m)
                  jdsm = -1

                  kdsm = itrca(m)

            do j = 1, itrcn(m)

                  jdsm = jdsm + 1
                  ntrn = idas_itrcc(idsm+jdsm)
                  jdsm = jdsm + 1
                  mtrn = idas_itrcc(idsm+jdsm)
                  kssm = idsm + jdsm + 1
                  jdsm = jdsm + mtrn

                  area = das_itrca(kdsm+j-1)

                  call echrg2(mtrn,idas_itrcc(kssm),dum1,lng1,icmb,igm)


C S.H. revised for adding r-from and r-to in [t-cross] (2017.9.27)
                     l = 4
                     dum2(1:l) = '    '
                     write(dum2(l+1:l+3),'(i3)') j
                     l = l + 3

                     dum2(l+1:l+7) = '       '
                     l = l + 7
                  do i = 1, lng1
                     dum2(i+l:i+l) = dum1(i:i)
                  end do
                     l = l + lng1
                     lngb = lgmax1 - lng1
                  do i = 1, lngb
                     dum2(i+l:i+l) = ' '
                  end do
                     l = l + lngb

                  jdsm = jdsm + 1
                  ntrn = idas_itrcc(idsm+jdsm)
                  jdsm = jdsm + 1
                  mtrn = idas_itrcc(idsm+jdsm)
                  kssm = idsm + jdsm + 1
                  jdsm = jdsm + mtrn

                  call echrg2(mtrn,idas_itrcc(kssm),dum1,lng1,icmb,igm)


C S.H. revised for adding r-from and r-to in [t-cross] (2017.9.27)
                     dum2(l+1:l+7) = '       '
                     l = l + 7
                  do i = 1, lng1
                     dum2(i+l:i+l) = dum1(i:i)
                  end do
                     l = l + lng1
                     lngb = lgmax2 - lng1
                  do i = 1, lngb
                     dum2(i+l:i+l) = ' '
                  end do
                     l = l + lngb


C S.H. revised for adding r-from and r-to in [t-cross] (2017.9.27)
                  write(dum2(l+1:l+20),
     &                      '(5x,1p1g15.7)') area
                     l = l + 20

                     write(iot,'(600a1)') (dum2(i:i),i=1,l)

            end do

*-----------------------------------------------------------------------

            else if( itmsh(m) .eq. 2 ) then

                  write(iot,'("     mesh =  r-z",11x,
     &            " # mesh type is r-z scoring mesh")')

                  write(iot,'("       x0 = ",1p1g14.7,1x,
     &            " # center x-position of r-z mesh")') rtrx0(m)

                  write(iot,'("       y0 = ",1p1g14.7,1x,
     &            " # center y-position of r-z mesh")') rtry0(m)

                  call echmty(0,iot,'r',itrty(m),itrnm(m),
     &                        rtrdl(m),rtrmi(m),rtrma(m),itrrg(m)
     &                        ,abs(itrnm(m)+1),das_itrrg)
                  call echmty(0,iot,'z',itzty(m),itznm(m),
     &                        rtzdl(m),rtzmi(m),rtzma(m),itzrg(m)
     &                        ,abs(itznm(m)+1),das_itzrg)

            else if( itmsh(m) .eq. 3 ) then

                  write(iot,'("     mesh =  xyz",11x,
     &            " # mesh type is xyz scoring mesh")')

                  call echmty(0,iot,'x',itxty(m),itxnm(m),
     &                        rtxdl(m),rtxmi(m),rtxma(m),itxrg(m)
     &                        ,abs(itxnm(m)+1),das_itxrg)
                  call echmty(0,iot,'y',ityty(m),itynm(m),
     &                        rtydl(m),rtymi(m),rtyma(m),ityrg(m)
     &                        ,abs(itynm(m)+1),das_ityrg)
                  call echmty(0,iot,'z',itzty(m),itznm(m),
     &                        rtzdl(m),rtzmi(m),rtzma(m),itzrg(m)
     &                        ,abs(itznm(m)+1),das_itzrg)

            end if

*-----------------------------------------------------------------------

            if( itmtr(m,1) .gt. 0 ) then

               if( itmtr(m,1) .eq. 1 ) then

                     write(iot,'("     trcl = ",i6,9x,
     &               " # transform ID")') itmtr(m,3)

               else

                  if( itmtr(m,2) .eq. 0 ) then

                     write(iot,'("     trcl = ",
     &                           3(1p1e15.7),
     &                           3(/12x,3(1p1e15.7)),i5)')
     &               ( rtmtr(m,j), j = 1, 12 ),
     &               nint( rtmtr(m,13) )

                  else

                     write(iot,'("    *trcl = ",
     &                           3(1p1e15.7),
     &                           3(/12x,3(1p1e15.7)),i5)')
     &               ( rtmtr(m,j), j = 1, 12 ),
     &               nint( rtmtr(m,13) )

                  end if

               end if

            end if

*-----------------------------------------------------------------------

               if( itout(m) .le. 4 .or. itout(m) .eq. 8 .or.
     &             itout(m) .eq. 10 ) then ! T.Sato 2023/04/13, add a-flux case

                  call echmty(0,iot,'e',itety(m),itenm(m),
     &                        rtedl(m),rtemi(m),rtema(m),iterg(m)
     &                        ,abs(itenm(m)+1),das_iterg)

               end if

*-----------------------------------------------------------------------

               if( itout(m) .ge. 8 .and. itout(m) .le. 11 ) then ! T.Sato 2023/04/13

                  call echmty(0,iot,'a',itaty(m),itanm(m),
     &                        rtadl(m),rtami(m),rtama(m),itarg(m)
     &                        ,abs(itanm(m)+1),das_itarg)

               end if

*-----------------------------------------------------------------------

               if( ittty(m) .ne. 0 ) then

                  call echmty(0,iot,'t',ittty(m),ittnm(m),
     &                        rttdl(m),rttmi(m),rttma(m),ittrg(m)
     &                        ,abs(ittnm(m)+1),das_ittrg)

               end if

*-----------------------------------------------------------------------

               write(iot,'("  eng2let =",i5,11x,
     &         " # (D=0) Flag to convert energy to LET, 0: not convert",
     &                                               ", 1: convert")')
     &         ite2l(m)

               if ( ite2l(m) .ne. 0 ) then
                  write(iot,'("   letmat =",i5,11x,
     &            " # (D=0) mat ID for LET, 0:real mat, ",
     &                                    "<0: electron for H2O")')
     &            itlmt(m)
               end if

*-----------------------------------------------------------------------

               if( itunt(m) .eq. 1 ) then

                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [1/cm^2/source]")') itunt(m)

               else if( itunt(m) .eq. 2 ) then

                if (ite2l(m) .eq. 0) then   ! ccse 2022/08/31
                if ( iMeVperu.eq. 1 ) then
                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [1/cm^2/(MeV/n)/source]")') itunt(m)
                else
                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [1/cm^2/MeV/source]")') itunt(m)
                end if
                else if (ite2l(m) .eq. 1) then   ! convert to energy to LET
                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [1/cm^2/(keV/um)/source]")') itunt(m)
                end if

               else if( itunt(m) .eq. 3 ) then

                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [1/cm^2/Lethargy/source]")') itunt(m)

               else if( itunt(m) .eq. 4 ) then

                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [1/cm^2/sr/source]")') itunt(m)

               else if( itunt(m) .eq. 5 ) then

                if (ite2l(m) .eq. 0) then   ! ccse 2022/08/31
                if ( iMeVperu.eq. 1 ) then
                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [1/cm^2/(MeV/n)/sr/source]")') itunt(m)
                else
                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [1/cm^2/MeV/sr/source]")') itunt(m)
                end if
                else if (ite2l(m) .eq. 1) then   ! convert to energy to LET
                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [1/cm^2/(keV/um)/sr/source]")') itunt(m)
                end if

               else if( itunt(m) .eq. 6 ) then

                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [1/cm^2/Lethargy/sr/source]")') itunt(m)

               else if( itunt(m) .eq. 11 ) then

                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [1/cm^2/nsec/source]")') itunt(m)

               else if( itunt(m) .eq. 12 ) then

                if (ite2l(m) .eq. 0) then   ! ccse 2022/08/31
                if ( iMeVperu.eq. 1 ) then
                  write(iot,'("     unit = ",i4,11x,
     &           " # unit is [1/cm^2/(MeV/n)/nsec/source]")') itunt(m)
                else
                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [1/cm^2/MeV/nsec/source]")') itunt(m)
                end if
                else if (ite2l(m) .eq. 1) then   ! convert to energy to LET
                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [1/cm^2/(keV/um)/nsec/source]")') itunt(m)
                end if

               else if( itunt(m) .eq. 13 ) then

                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [1/cm^2/Lethargy/nsec/source]")')
     &            itunt(m)

               else if( itunt(m) .eq. 14 ) then

                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [1/cm^2/nsec/sr/source]")') itunt(m)

               else if( itunt(m) .eq. 15 ) then

                if (ite2l(m) .eq. 0) then   ! ccse 2022/08/31
                if ( iMeVperu.eq. 1 ) then
                  write(iot,'("     unit = ",i4,11x,
     &        " # unit is [1/cm^2/(MeV/n)/nsec/sr/source]")') itunt(m)
                else
                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [1/cm^2/MeV/nsec/sr/source]")') itunt(m)
                end if
                else if (ite2l(m) .eq. 1) then   ! convert to energy to LET
                  write(iot,'("     unit = ",i4,11x,
     &         " # unit is [1/cm^2/(keV/um)/nsec/sr/source]")') itunt(m)
                end if

               else if( itunt(m) .eq. 16 ) then

                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [1/cm^2/Lethargy/nsec/sr/source]")')
     &            itunt(m)

               end if

*-----------------------------------------------------------------------

            if( icax .eq. 0 ) then

                  icn = itaxn(m)

               do i = 1, icn

                  ian(i) = i

               end  do

            else

                  icn = 1
                  ian(1) = iaxs

            end if

*-----------------------------------------------------------------------

               idtyp = 0

            do i = 1, icn

               if( itaxs(m,ian(i)) .eq. 7 .or.
     &             itaxs(m,ian(i)) .eq. 11 .or.
     &             itaxs(m,ian(i)) .eq. 12 .or.
     &             itaxs(m,ian(i)) .eq. 13 ) idtyp = ittwo(m)

            end do

            if( idtyp .ne. 0 ) then

               write(iot,'("  2D-type = ",i4,11x,
     &         " # 1:Cont, 2:Clust, 3:Color, 4:xyz, 5:mat,",
     &         " 6:Clust+Cont, 7:Col+Cont")') idtyp

            end if

*-----------------------------------------------------------------------

                  asfil = '  # file name of output for the above axis'
                  bsfil = '  # file name of dumped data'
                  csfil = '  # file name of dumped data summary'

            do i = 1, icn

                  iax = ian(i)

                  write(iot,'("     axis =  ",a3,11x,
     &            " # axis of output")')
     &            aname( itaxs(m,iax) )

             if ( itfll(m,iax) .gt. 0 ) then

               if( itmdp(m,0) .eq. 0 ) then

                  msfile = max( 14, itfll(m,iax) )

                  write(iot,'("     file = ",100a1)')
     &                  ( ctfln(m,iax)(j:j), j = 1, msfile ),
     &                  ( asfil(j:j), j = 1, 42 )

               else

C S.H. added for Dump-Restart on 2014/5/7
                  msfile = max( 14, itfll(m,iax) )

                  write(iot,'("     file = ",100a1)')
     &                  ( ctfln(m,iax)(j:j), j = 1, msfile ),
     &                  ( csfil(j:j), j = 1, 36 )

                  msfile = max( 14, itfll(m,iax)+4 )

                  do k = 1, msfile
                     dsfil(k:k) = ' '
                  end do

                  do is = itfll(m,iax), 1, -1
                     if ( ctfln(m,iax)(is:is) .eq. '.' ) exit
                  end do
                  if ( is .eq. 0 ) is = itfll(m,iax)+1

                  dsfil(1:is-1) = ctfln(m,iax)(1:is-1)
                  dsfil(is:is+3) = '_dmp'
                  if ( is .lt. itfll(m,iax) ) then
                     dsfil(is+4:itfll(m,iax)+4)
     &                    = ctfln(m,iax)(is:itfll(m,iax))
                  end if

                  write(iot,'("#    file = ",100a1)')
     &                  ( dsfil(j:j), j = 1, msfile ),
     &                  ( bsfil(j:j), j = 1, 28 )

               end if

             end if

            end do

*-----------------------------------------------------------------------

               if( itout(m) .eq. 1 ) then

                  write(iot,'("   output = flux",10x,
     &            "  # surface crossing flux")')

               else if( itout(m) .eq. 2 ) then

                  write(iot,'("   output = current",7x,
     &            "  # surface crossing current spectrum")')

               else if( itout(m) .eq. 3 ) then

                  write(iot,'("   output = f-curr",8x
     &            "  # surface crossing forward current spectrum")')

               else if( itout(m) .eq. 4 ) then

                  write(iot,'("   output = b-curr",8x
     &            "  # surface crossing backward current spectrum")')

               else if( itout(m) .eq. 5 ) then

                  write(iot,'("   output = o-curr",8x
     &            "  # surface crossing omni current")')

               else if( itout(m) .eq. 6 ) then

                  write(iot,'("   output = of-curr",7x,
     &            "  # surface crossing omni forward current")')

               else if( itout(m) .eq. 7 ) then

                  write(iot,'("   output = ob-curr",7x,
     &            "  # surface crossing omni backward current")')

               else if( itout(m) .eq. 8 ) then

                  write(iot,'("   output = a-curr",7x,
     &           "  # surface crossing current spectrum with angle")')

               else if( itout(m) .eq. 9 ) then

                  write(iot,'("   output = oa-curr",7x,
     &            "  # surface crossing omni current with angle")')

               else if( itout(m) .eq. 10 ) then

                  write(iot,'("   output = a-flux",8x,
     &            "  # surface crossing flux with angle")')

               else if( itout(m) .eq. 11 ) then

                  write(iot,'("   output = oa-flux",7x,
     &            "  # surface crossing omni flux with angle")')

               end if

*-----------------------------------------------------------------------
               call echprt(1,iot,m,itpan,itpat,jtpat,itnm,itmxpt,itmxgp)

*-----------------------------------------------------------------------
*        multiplier
*-----------------------------------------------------------------------

         if( itmlp(m) .gt. 0 ) then

            do k = 1, itmlp(m)

*-----------------------------------------------------------------------

               if( itmln(m,k) .eq. 1 .and.
     &             mltp(1,itmli(m,k)/13+1) .eq. -1 ) then

                  imall = 0

                  write(iot,'("  multiplier = all",9x,
     &            " # number of material group")')

               else

                  imall = itmln(m,k)

                  write(iot,'("  multiplier = ",i3,9x,
     &            " # number of material group")') imall

               end if

*-----------------------------------------------------------------------

                     ipan(m) = itmpn(m,k)

                  do j = 1, itmpn(m,k)

                     ipat(m,j,1) = itmpt(m,k,j,1)
                     ipat(m,j,2) = itmpt(m,k,j,2)

                  end do

                 call echprt(1,iot,m,ipan,ipat,jtpat,itnm,itmxpt,itmxgp)

               if( rtmme(m,k) .gt. 0.0 ) then

                  write(iot,'("     emax = ",1p1g13.5)') rtmme(m,k)

               end if

*-----------------------------------------------------------------------

                     do i = 1, itmst(m)

                        ilgt(i) = 0

                     end do

               do j = 1, itmln(m,k)

                  do im = 1, itmst(m)

                        nn = im * 2 - 1
                        jlgt(im) = 0

                  do jj = mltp(1+nn,itmli(m,k)/13+j),
     &                    mltp(1+nn+1,itmli(m,k)/13+j), 2

                        if( nint(slib(jj)) .eq. 1 .or.
     &                      nint(slib(jj)) .eq. 2 .or.
     &                      nint(slib(jj)) .eq. 3 ) then

                           jlgt(im) = jlgt(im) + 2

                        else if( nint(slib(jj)) .eq. 4 .or.
     &                           nint(slib(jj)) .eq. 8 ) then

                           write(ch13,'(1p1g13.5)') slib(jj+1)

                           do l = 1, 13
                              if( ch13(l:l) .ne. ' ' ) goto 410
                           end do
  410                      isi = l - 1
                           do l = 13, 1, -1
                              if( ch13(l:l) .ne. ' ' ) goto 420
                           end do
  420                      isf = 13 - l

                           jlgt(im) = jlgt(im) + 13 - isi - isf + 1

                        else if( nint(slib(jj)) .ge. 5 .and.
     &                           nint(slib(jj)) .le. 7 ) then

                           write(ch6,'(i6)') nint( slib(jj+1) )

                           do l = 1, 6
                              if( ch6(l:l) .ne. ' ' ) goto 430
                           end do
  430                      isi = l - 1

                           jlgt(im) = jlgt(im) + 6 - isi + 1

                        end if

                  end do
                  end do

                     do i = 1, itmst(m)

                        if( jlgt(i) .gt. ilgt(i) ) ilgt(i) = jlgt(i)

                     end do

               end do

                     do i = 1, itmst(m)

                        if( ilgt(i) .lt. 7 ) ilgt(i) = 7

                     end do

                  ch200(1:9) = '      mat'
                  ic = 9

               do i = 1, itmst(m)

                  isp = ilgt(i) - 7

                  do j = 1, isp
                     ch200(ic+j:ic+j) = ' '
                  end do
                     ic = ic + isp
                     ch200(ic+1:ic+4) = 'mset'
                     ic = ic + 4

                  write(ch200(ic+1:ic+1),'(i1)') itmnt(m,i)
                     ic = ic + 1

                  ch200(ic+1:ic+2) = '  '
                     ic = ic + 2

               end do

                  write(iot,'(200a1)') (ch200(j:j),j=1,ic)

*-----------------------------------------------------------------------

               do j = 1, itmln(m,k)

                     do ik = 1, 400
                        ch100(ik:ik) = ' '
                     end do

                     ch200(1:3) = '   '
                     ic = 3

                  if( mltp(1,itmli(m,k)/13+j) .gt. 0 ) then
                     write(ch200(ic+1:ic+6),'(i6)')
     &                  mltp(1,itmli(m,k)/13+j)

                   else

                     ch200(ic+1:ic+6) = '   all'

                   end if

                     ic = ic + 6

                  do im = 1, itmst(m)

                        nn = im * 2 - 1
                        id = 0

                  do jj = mltp(1+nn,itmli(m,k)/13+j),
     &                    mltp(1+nn+1,itmli(m,k)/13+j), 2

                        if( nint(slib(jj)) .eq. 1 ) then

                           ch100(id+1:id+2) = ' ('
                           id = id + 2

                        else if( nint(slib(jj)) .eq. 2 ) then

                           ch100(id+1:id+2) = ' )'
                           id = id + 2

                        else if( nint(slib(jj)) .eq. 3 ) then

                           ch100(id+1:id+2) = ' :'
                           id = id + 2

                        else if( nint(slib(jj)) .eq. 4 .or.
     &                           nint(slib(jj)) .eq. 8 ) then

                           write(ch13,'(1p1g13.5)') slib(jj+1)

                           do l = 1, 13
                              if( ch13(l:l) .ne. ' ' ) goto 440
                           end do
  440                      isi = l
                           do l = 13, 1, -1
                              if( ch13(l:l) .ne. ' ' ) goto 450
                           end do
  450                      isf = l

                           ch100(id+1:id+1) = ' '
                           id = id + 1
                           ch100(id+1:id+isf-isi+1) = ch13(isi:isf)
                           id = id+isf-isi+1

                        else if( nint(slib(jj)) .ge. 5 .and.
     &                           nint(slib(jj)) .le. 7 ) then

                           write(ch6,'(i6)') nint( slib(jj+1) )

                           do l = 1, 6
                              if( ch6(l:l) .ne. ' ' ) goto 460
                           end do
  460                      isi = l

                           ch100(id+1:id+1) = ' '
                           id = id + 1
                           ch100(id+1:id+6-isi+1) = ch6(isi:6)
                           id = id+6-isi+1

                        end if

                  end do

                           isd = ilgt(im) - id

                           do kk = 1, isd
                              ch200(ic+kk:ic+kk) = ' '
                           end do
                           ic = ic + isd

                           ch200(ic+1:ic+id) = ch100(1:id)
                           ic = ic + id

                  end do

                  write(iot,'(200a1)') (ch200(kk:kk),kk=1,ic)

               end do

*-----------------------------------------------------------------------

            end do

         end if

*-----------------------------------------------------------------------

               if( rtfac(m) .ne. 1.0d0 ) then

                  write(iot,'("   factor =",1p1g14.7,1x,1x,
     &            " # (D=1.0) normalization factor")')
     &            rtfac(m)

               end if

*-----------------------------------------------------------------------

               if( itanl(m) .gt. 0 ) then

                  write(iot,'("    angel = ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               end if

*-----------------------------------------------------------------------

               if( itsans(m) .gt. 0 ) then

                  write(iot,'("   sangel = ",i2)') itsans(m)
                  call write_sangel(iot,m,1)

               end if

*-----------------------------------------------------------------------

               if( itaxl(m) .gt. 0 ) then

                  write(iot,'("    x-txt = ",200a1)')
     &            ( itaxt(m)(i:i),i = 1, itaxl(m) )

               end if

*-----------------------------------------------------------------------

               if( itayl(m) .gt. 0 ) then

                  write(iot,'("    y-txt = ",200a1)')
     &            ( itayt(m)(i:i),i = 1, itayl(m) )

               end if

*-----------------------------------------------------------------------

               if( itazl(m) .gt. 0 ) then

                  write(iot,'("    z-txt = ",200a1)')
     &            ( itazt(m)(i:i),i = 1, itazl(m) )

               end if

*-----------------------------------------------------------------------

               if( itgsh(m) .gt. 0 ) then

                  write(iot,'("    gshow = ",i4,10x,
     &            "  # 0: no 1:bnd, 2:bnd+mat, ",
     &     "3:bnd+reg 4:bnd+lat 5:bmp style")')
     &            itgsh(m)

               end if

*-----------------------------------------------------------------------

               if( itglt(m) .ge. 0 ) then

                  write(iot,'("    gslat = ",i4,10x,
     &            "  # 0: no line for lat, 1: with,",
     &               " 2: no line between the same cell")')
     &            itglt(m)

               end if

*-----------------------------------------------------------------------

               if( itgsh(m) .gt. 0 .and. itger(m) .gt. 0 ) then

                  write(iot,'("    ginfo = ",i4,11x,
     &            " # (D=0) 0:no, 1:show in graph, 2:add file")')
     &            itger(m)

               end if

*-----------------------------------------------------------------------

               if( itgsh(m) .gt. 0 ) then

                  write(iot,'("    resol = ",i4,10x,
     &            "  # (D=1) resolution of gshow")')
     &            itres(m)

                  write(iot,'("    width =",1p1g14.7,1x,1x,
     &            " # (D=0.5) width of lines for gshow")')
     &            rtwid(m)

               end if

*-----------------------------------------------------------------------

               if( iteps(m) .ne. 0 ) then

                  write(iot,'("   epsout =",i5,10x,
     &            "  # (D=0) generate eps file by ANGEL")')
     &            iteps(m)

               end if

*-----------------------------------------------------------------------
*        check of counter
*-----------------------------------------------------------------------

            do i = 1, 3

               if( itcnt(i,m) .ne. 0 ) then

                  write(iot,'(" ctmin(",i1,") =",i5,10x,
     &            "  # (D=-9999) min. of counter ",i1)')
     &            i, itcnt(i*2+2,m), i

                  write(iot,'(" ctmax(",i1,") =",i5,10x,
     &            "  # (D= 9999) max. of counter ",i1)')
     &            i, itcnt(i*2+3,m), i

                  if( ncntc(i) .eq. 0 )
     &            write(iot,'("# Warning: there is no counter = ",
     &                        i1," definition in [counter].")') i

               end if

            end do

            do i = 1, 3
               if( itcnth(i,m) .ne. 0 ) then
                  write(iot,'(" chmin(",i1,") =",i5,10x,
     &            "  # (D=-9999) min. of history counter ",i1)')
     &            i, itcnth(i*2+2,m), i
                  write(iot,'(" chmax(",i1,") =",i5,10x,
     &            "  # (D= 9999) max. of history counter ",i1)')
     &            i, itcnth(i*2+3,m), i
               end if
            end do

*-----------------------------------------------------------------------

               if( rtstd(m) .gt. 0.d0 ) then

                  write(iot,'("   stdcut =",1p1g14.7,1x,1x,
     &            " # (D=-1) threshold value of STD cutoff")')
     &            rtstd(m)

               end if

*-----------------------------------------------------------------------

               if( itenclo(m) .ne. 0 ) then

                  write(iot,'("   enclos =",i5,10x,
     &            "  # (D=0) 0-> normal mode,",
     &                     " 1-> enclosure mode")')
     &            itenclo(m)

               end if

*-----------------------------------------------------------------------

               if( itangform(m) .eq. 1 ) then

                  write(iot,'(" iangform =",i5,10x,
     &            "  # 1-> angle formed from the x-axis")')
     &            itangform(m)

               elseif( itangform(m) .eq. 2 ) then

                  write(iot,'(" iangform =",i5,10x,
     &            "  # 2-> angle formed from the y-axis")')
     &            itangform(m)

               elseif( itangform(m) .eq. 3 ) then

                  write(iot,'(" iangform =",i5,10x,
     &            "  # 3-> angle formed from the z-axis")')
     &            itangform(m)

               end if

*-----------------------------------------------------------------------
               if( itmxang(m) .ne. 0 ) then

                  write(iot,'(" maxangel =",i5,10x,
     &            "  # (D=0) number of output particles")')
     &            itmxang(m)

               end if

*-----------------------------------------------------------------------
*        dump the data
*-----------------------------------------------------------------------

         if( itmdp(m,0) .ne. 0 ) then

                  write(iot,'("     dump =",i5,11x,
     &            " # (D=0) number of dumped data,",
     &            " <0: ascii, >0: binary")') itmdp(m,0)

                  write(iot,'(13x,30(i4))')
     &               ( itmdp(m,j), j = 1, abs( itmdp(m,0) ) )

                  write(iot,'("# dump data  ",30(a4))')
     &               ( dmpc(itmdp(m,j)), j = 1, abs( itmdp(m,0) ) )

         end if

*-----------------------------------------------------------------------
*        sumtally subsection
*-----------------------------------------------------------------------
            if ( icntl==13 .and. nsumtalRead(m)==2 ) then

               call sumtal_echo(iot,m)

            end if
*-----------------------------------------------------------------------

            write(iot,'("#    used :",
     &                  "        main (  %)",
     &                  "        temp (  %)",
     &                  "       total (  %)"/
     &                  "#  memory :",
     &                  3(i12," (",i3,")"))')
     &                  itsmn(m), nint(dble(itsmn(m))/dble(mdas)*100),
     &                  itstm(m), nint(dble(itstm(m))/dble(mdas)*100),
     &                  itsmn(m)+itstm(m),
     &                  nint(dble(itsmn(m)+itstm(m))/dble(mdas)*100)

*-----------------------------------------------------------------------
cfrtati 2024/02/13
            if( itmlp(m).gt.0 ) then
              call displaymtinfo(iot,m)
            end if

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine tyilech(iot,m,iaxs,icax)
*                                                                      *
*       input echo of the yield tally                                  *
*       last modified by K.Niita on 2011/11/11                         *
*                                                                      *
************************************************************************

      use sumtallymod
      use sangelmod, only: itsans
      use partmod, only: itmxpt, itmxgp, itpan, itpat, jtpat ! frtati 2021/10/05
      use moddas_mesh
      use moddas_region
      use moddas_tally
*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'

*-----------------------------------------------------------------------

      common /tall00/ itnm, ital(itlmax), itals(itlmax), italm(itlmax)
      common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                itout(itlmax), ittwo(itlmax)

      common /tall02/ itrgn(itlmax), itrgm(itlmax), itreg(itlmax)
      common /tall24/ itrcm(itlmax), itrcg(itlmax)
      common /tall25/ itrnv(itlmax), itriv(itlmax), itrrv(itlmax)

      common /tall03/ itrty(itlmax), itrnm(itlmax), itrrg(itlmax),
     &                rtrmi(itlmax), rtrma(itlmax), rtrdl(itlmax)
      common /tall04/ itxty(itlmax), itxnm(itlmax), itxrg(itlmax),
     &                rtxmi(itlmax), rtxma(itlmax), rtxdl(itlmax)
      common /tall05/ ityty(itlmax), itynm(itlmax), ityrg(itlmax),
     &                rtymi(itlmax), rtyma(itlmax), rtydl(itlmax)
      common /tall06/ itzty(itlmax), itznm(itlmax), itzrg(itlmax),
     &                rtzmi(itlmax), rtzma(itlmax), rtzdl(itlmax)

      common /tall08/ rtrx0(itlmax), rtry0(itlmax)

      common /tall10/ itaxn(itlmax), itaxs(itlmax,6)
      common /tall11/ itfln(itlmax), itfll(itlmax,6), ctfln(itlmax,6)
      character ctfln*100

      common /tall12/ itman(itlmax), itmat(itlmax), itmct(itlmax),
     &                itnun(itlmax), itnuc(itlmax),
     &                itndz(itlmax), itndn(itlmax),
     &                itnkz(itlmax), itnkn(itlmax)

      common /tall16/ itnfr(itlmax), itnfz(itlmax), itnfn(itlmax)

      common /tall14/ ittll(itlmax), ittle(itlmax)
      character ittle*80

      common /tall15/ itanl(itlmax), itang(itlmax)
      character itang*200

      common /tall27/ itaxl(itlmax), itaxt(itlmax),
     &                itayl(itlmax), itayt(itlmax),
     &                itazl(itlmax), itazt(itlmax)
      character itaxt*200, itayt*200, itazt*200

      common /tall19/ itsmn(itlmax), itstm(itlmax)

      common /tall21/ rtfac(itlmax)
      common /tall22/ itmcn(itlmax), itmtn(itlmax), itmtt(itlmax)
      common /tall28/ itgsh(itlmax)
      common /tall29/ itrsh(itlmax)
      common /tall30/ itnda(itlmax)
      common /tall31/ iterl(itlmax)
      common /tall34/ itvm(itlmax)
      common /tall35/ iteps(itlmax)
      common /tall36/ itdpo(itlmax)

      common /tall37/ itcnt(9,itlmax)

      common /cntreg/ ncntc(3), ncreg(3), incrc(3), incrt(3), kcont(3)
      common /cntpat/ icpan(3), icpat(3,20,2)

      common /tall38/ itres(itlmax)
      common /tall39/ rtwid(itlmax)
      common /tall46/ itmtr(itlmax,4), rtmtr(itlmax,13)
      common /tall49/ itglt(itlmax)
      common /tall52/ itprd(itlmax)

      common /tall60/ itger(itlmax)
      common /tall63/ itbmp(itlmax)
      common /tall64/ itvtk(itlmax),itvtkfmt(itlmax)

      common /tall67/ itstd(itlmax), rtstd(itlmax)

      integer :: itfoam
      common /tall76/ itfoam(itlmax)

      common /tall82/ itcnth(9,itlmax)

      common /tall83/ itnzn(itlmax), itndm(itlmax)

      dimension idas(mdas*2) !20220906frtati temporary for fbounds-check
      equivalence ( das, idas )

*-----------------------------------------------------------------------

      common /inpec/ ititl, ipara, ibody, iregn, llarr, itby, itar
      common /kmat1d/ idmn(0:kvlmax), idnm(kvmmax)

*-----------------------------------------------------------------------
       integer   icntl
       integer   inucr

       common /tcntl/  icntl, inucr
!-----------------------------------------------------------------------

      character asfil*100
      character esfil*100
      character erfnm*100
      character aname(14)*7
      data aname / 'mass   ',' reg   ','   x   ','   y   ',
     &             '   z   ','   r   ','charge ','chart  ',
     &             '  xy   ','  yz   ','  xz   ','  rz   ',
     &             'dchain ',' tet   '/

      dimension ian(6)
      dimension jmat(20)
      character chau*8
      character chaus(7)*8
      character elmnt(104)*3

*-----------------------------------------------------------------------

      data elmnt /
     &    'H  ','He ','Li ','Be ','B  ','C  ','N  ','O  ','F  ','Ne ',
     &    'Na ','Mg ','Al ','Si ','P  ','S  ','Cl ','Ar ','K  ','Ca ',
     &    'Sc ','Ti ','V  ','Cr ','Mn ','Fe ','Co ','Ni ','Cu ','Zn ',
     &    'Ga ','Ge ','As ','Se ','Br ','Kr ','Rb ','Sr ','Y  ','Zr ',
     &    'Nb ','Mo ','Tc ','Ru ','Rh ','Pd ','Ag ','Cd ','In ','Sn ',
     &    'Sb ','Te ','I  ','Xe ','Cs ','Ba ','La ','Ce ','Pr ','Nd ',
     &    'Pm ','Sm ','Eu ','Gd ','Tb ','Dy ','Ho ','Er ','Tm ','Yb ',
     &    'Lu ','Hf ','Ta ','W  ','Re ','Os ','Ir ','Pt ','Au ','Hg ',
     &    'Tl ','Pb ','Bi ','Po ','At ','Rn ','Fr ','Ra ','Ac ','Th ',
     &    'Pa ','U  ','Np ','Pu ','Am ','Cm ','Bk ','Cf ','Es ','Fm ',
     &    'Md ','No ','Lr ','Ku ' /

*-----------------------------------------------------------------------
      common /tall79/ itallech
      if( itallech .eq. 0 ) return

*-----------------------------------------------------------------------
*     input echo
*-----------------------------------------------------------------------
            if( ital(m+2) .eq. 16 .and. icax .eq. 0 ) then
                  write(iot,'("[ T-Yield ] off")')
            else
                  write(iot,'("[ T-Yield ]")')
            end if

                  write(iot,'("    title = ",80a1)')
     &            ( ittle(m)(i:i), i = 1, ittll(m) )

*-----------------------------------------------------------------------

            if( itmsh(m) .eq. 1 ) then

                  write(iot,'("     mesh =  reg",11x,
     &            " # mesh type is region-wise")')

                     idas1 = mmmax
                     idas2 = ( idas1 + itrgn(m) - 1 ) * 2 + 1
                     idas3 = idas2 + itrgn(m)

                  call echrg(1,iterl(m),0,1,iot,itrcm(m),
     &                       idas_itrcg(itrcg(m)),
     &                       itrgn(m),itrgm(m),idas_itreg(itreg(m)),
     &                       das(idas1),idas(idas2),
     &                       itrnv(m),idas(itriv(m)),das(itrrv(m)),
     &                       idas3)

            else if( itmsh(m) .eq. 2 ) then

                  write(iot,'("     mesh =  r-z",11x,
     &            " # mesh type is r-z scoring mesh")')

                  write(iot,'("       x0 = ",1p1g14.7,1x,
     &            " # center x-position of r-z mesh")') rtrx0(m)

                  write(iot,'("       y0 = ",1p1g14.7,1x,
     &            " # center y-position of r-z mesh")') rtry0(m)

                  call echmty(0,iot,'r',itrty(m),itrnm(m),
     &                        rtrdl(m),rtrmi(m),rtrma(m),itrrg(m)
     &                        ,abs(itrnm(m)+1),das_itrrg)
                  call echmty(0,iot,'z',itzty(m),itznm(m),
     &                        rtzdl(m),rtzmi(m),rtzma(m),itzrg(m)
     &                        ,abs(itznm(m)+1),das_itzrg)

            else if( itmsh(m) .eq. 3 ) then

                  write(iot,'("     mesh =  xyz",11x,
     &            " # mesh type is xyz scoring mesh")')

                  call echmty(0,iot,'x',itxty(m),itxnm(m),
     &                        rtxdl(m),rtxmi(m),rtxma(m),itxrg(m)
     &                        ,abs(itxnm(m)+1),das_itxrg)
                  call echmty(0,iot,'y',ityty(m),itynm(m),
     &                        rtydl(m),rtymi(m),rtyma(m),ityrg(m)
     &                        ,abs(itynm(m)+1),das_ityrg)
                  call echmty(0,iot,'z',itzty(m),itznm(m),
     &                        rtzdl(m),rtzmi(m),rtzma(m),itzrg(m)
     &                        ,abs(itznm(m)+1),das_itzrg)

            else if( itmsh(m) .eq. 4 ) then

                  write(iot,'("     mesh =  tet",11x,
     &            " # mesh type is tetra geometry")')

                  call echtet(iot,itrgm(m),idas_itreg(itreg(m)))

            end if

*-----------------------------------------------------------------------

            if( itmtr(m,1) .gt. 0 ) then

               if( itmtr(m,1) .eq. 1 ) then

                     write(iot,'("     trcl = ",i6,9x,
     &               " # transform ID")') itmtr(m,3)

               else

                  if( itmtr(m,2) .eq. 0 ) then

                     write(iot,'("     trcl = ",
     &                           3(1p1e15.7),
     &                           3(/12x,3(1p1e15.7)),i5)')
     &               ( rtmtr(m,j), j = 1, 12 ),
     &               nint( rtmtr(m,13) )

                  else

                     write(iot,'("    *trcl = ",
     &                           3(1p1e15.7),
     &                           3(/12x,3(1p1e15.7)),i5)')
     &               ( rtmtr(m,j), j = 1, 12 ),
     &               nint( rtmtr(m,13) )

                  end if

               end if

            end if

*-----------------------------------------------------------------------

                  write(iot,'("  special = ",i8,7x,
     &            " # (D=0) number of repeated nuclear reactions")')
     &            itspc(m)

*-----------------------------------------------------------------------

               if( itmtn(m) .eq. 0 ) then

                  write(iot,'(" material =  all"11x,
     &            " # (D=all) number of specific material")')

               else if( itmtn(m) .gt. 0 ) then

                  write(iot,'(" material =",i5,11x,
     &            " # (D=all) number of specific material")')
     &            itmtn(m) * itmcn(m)

                     k = 0

                  do j = 1, itmtn(m)

                     k = k + 1

                     jmat(k) = ismte( itmtt(m) + j - 1 )

                     if( k / 10 * 10 .eq. k .or. j .eq. itmtn(m) ) then

                        write(iot,'(11x,10(i6))') ( jmat(i),i = 1, k )

                        k = 0

                     end if

                  end do

               end if

*-----------------------------------------------------------------------

               if( itman(m) .eq. 0 ) then

                  write(iot,'("   mother =  all"11x,
     &            " # (D=all) number of specific mother")')

               else if( itman(m) .gt. 0 ) then

                  write(iot,'("   mother = ",i4,11x,
     &            " # (D=all) number of specific mother")')
     &            itman(m) * itmct(m)

                     k = 0

                  do j = 1, itman(m)

                     k = k + 1

                     iaz = ismat( itmat(m) + j - 1 )

                     iz = iaz / 1000
                     ia = iaz - iz * 1000

                     if( ia .eq. 0 ) then

                        chaus(k) = elmnt(iz)

                     else

                        call chname(idum,ia,iz,chau)

                        chaus(k) = chau

                     end if

                     if( k / 7 * 7 .eq. k .or. j .eq. itman(m) ) then

                        write(iot,'(12x,7(a8))') ( chaus(i),i = 1, k )

                        k = 0

                     end if

                  end do

               end if

*-----------------------------------------------------------------------

               if( itnun(m) .eq. 0 ) then

                  write(iot,'("  nucleus =  all"11x,
     &            " # (D=all) number of detected nucleus")')

               else if( itnun(m) .gt. 0 ) then

                  write(iot,'("  nucleus = ",i4,11x,
     &            " # (D=all) number of detected nucleus")')
     &            itnun(m)

                     k = 0

                  do j = 1, itnun(m)

                     k = k + 1

                     iaz = isnuc( itnuc(m) + j - 1 )

                     iz = iaz / 1000
                     ia = iaz - iz * 1000

                     if( ia .eq. 0 ) then

                        chaus(k) = elmnt(iz)

                     else

                        call chname(idum,ia,iz,chau)

                        chaus(k) = chau

                     end if


                     if( k / 7 * 7 .eq. k .or. j .eq. itnun(m) ) then

                        write(iot,'(12x,7(a8))') ( chaus(i),i = 1, k )

                        k = 0

                     end if

                  end do

               end if

*-----------------------------------------------------------------------

               if( itunt(m) .eq. 1 ) then

                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [1/source]")') itunt(m)

               else if( itunt(m) .eq. 2 ) then

                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [1/cm^3/source]")') itunt(m)

               end if

*-----------------------------------------------------------------------

               if( itnda(m) .eq. 0 ) then

                  write(iot,'("    ndata = ",i4,11x,
     &            " # 0:no, this option is not used")') itnda(m)

               else if( itnda(m) .eq. 1 ) then

                  write(iot,'("    ndata = ",i4,11x,
     &            " # 0:no, 1:use ndata for 4He, 14N, 16O")') itnda(m)

               else if( itnda(m) .eq. 2 ) then

                  write(iot,'("    ndata = ",i4,11x,
     &            " # 0:no, 2:use yield data by High Energy File")')
     &               itnda(m)

               else if( itnda(m) .eq. 3 ) then

                  write(iot,'("    ndata = ",i4,11x,
     &            " # 0:no, 3:use yield data by User defined File")')
     &               itnda(m)

               end if

*-----------------------------------------------------------------------

                  asfil = '  # file name of output for the above axis'
                  esfil = '  # file name of output for the error'

            if( icax .eq. 0 ) then

                  icn = itaxn(m)

               do i = 1, icn

                  ian(i) = i

               end  do

            else

                  icn = 1
                  ian(1) = iaxs

            end if

*-----------------------------------------------------------------------

               idtyp = 0

            do i = 1, icn

               if( itaxs(m,ian(i)) .ge. 8 .and.
     &             itaxs(m,ian(i)) .le. 12 )  idtyp = ittwo(m)

            end do

            if( idtyp .ne. 0 ) then

               write(iot,'("  2D-type = ",i4,11x,
     &         " # 1:Cont, 2:Clust, 3:Color, 4:xyz, 5:mat,",
     &         " 6:Clust+Cont, 7:Col+Cont")') idtyp

            end if

*-----------------------------------------------------------------------

            do i = 1, icn

                  iax = ian(i)

                  write(iot,'("     axis = ",a7,8x,
     &            " # axis of output")')
     &            aname( itaxs(m,iax) )

             if ( itfll(m,iax) .gt. 0 ) then

                  msfile = max( 14, itfll(m,iax) )

                  write(iot,'("     file = ",100a1)')
     &                  ( ctfln(m,iax)(j:j), j = 1, msfile ),
     &                  ( asfil(j:j), j = 1, 42 )

               if( i .eq. 1 .and. itaxs(m,iax) .eq. 8 ) then

                  write(iot,'("     info = ",i4,11x,
     &            " # 0-> no,",
     &            " 1-> write magic number and stable nuclei")')
     &            itout(m)

               end if

               if( i .eq. 1 .and. itaxs(m,iax) .eq. 13 ) then

                  write(iot,'("     info = ",i4,11x,
     &            " # 0-> no, 1-> write error")') itout(m)

                  if( itout(m) .ne. 0 ) then

                     do j = itfll(m,iax), 1, -1

                        if( ctfln(m,iax)(j:j) .eq. '.' ) goto 50

                     end do

                        j = itfll(m,iax)

   50                   itfp = j - 1

                     do j = 1, itfp

                        erfnm(j:j) = ctfln(m,iax)(j:j)

                     end do

                        erfnm(itfp+1:itfp+4) = '.err'

                     do j = itfp + 5, 100

                        erfnm(j:j) = ' '

                     end do

                        msfile = max( 14, itfp+4 )

                        write(iot,'("#   efile = ",100a1)')
     &                  ( erfnm(j:j), j = 1, msfile ),
     &                  ( esfil(j:j), j = 1, 37 )

                  end if

               end if

             end if

            end do

*-----------------------------------------------------------------------

               if( itpan(m) .gt. 0 ) then

               call echprt(1,iot,m,itpan,itpat,jtpat,itnm,itmxpt,itmxgp)

               end if

*-----------------------------------------------------------------------

               if( itprd(m) .eq. 0 ) then

                  write(iot,'("   output = ","product        ",
     &            " # (D=product) normal yeild")')

               else if( itprd(m) .eq. 1 ) then

                  write(iot,'("   output = ","cutoff         ",
     &            " # (D=product) cutoff particles")')

               end if

*-----------------------------------------------------------------------

               if( itdpo(m) .eq. 0 ) then

                  write(iot,'("  elastic =    0","           ",
     &            " # (D=1) without elastic collision")')

               else if( itdpo(m) .eq. -1 ) then

                  write(iot,'("  elastic =   -1","           ",
     &            " # (D=1) without the same target")')

               else

                  write(iot,'("  elastic =    1","           ",
     &            " # (D=1) with elastic collision")')

               end if

*-----------------------------------------------------------------------

               if( rtfac(m) .ne. 1.0d0 ) then

                  write(iot,'("   factor =",1p1g14.7,1x,1x,
     &            " # (D=1.0) normalization factor")')
     &            rtfac(m)

               end if

*-----------------------------------------------------------------------

               if( itmtn(m) .gt. 0 .and. itmsh(m) .eq. 3 ) then

                  write(iot,'("   volmat = ",i4,11x,
     &            " # (D=9) vol.corr., 0:no, number of scan / dir.")')
     &            itvm(m)

               end if

*-----------------------------------------------------------------------

               if( itanl(m) .gt. 0 ) then

                  write(iot,'("    angel = ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               end if

*-----------------------------------------------------------------------

               if( itsans(m) .gt. 0 ) then

                  write(iot,'("   sangel = ",i2)') itsans(m)
                  call write_sangel(iot,m,1)

               end if

*-----------------------------------------------------------------------

               if( itaxl(m) .gt. 0 ) then

                  write(iot,'("    x-txt = ",200a1)')
     &            ( itaxt(m)(i:i),i = 1, itaxl(m) )

               end if

*-----------------------------------------------------------------------

               if( itayl(m) .gt. 0 ) then

                  write(iot,'("    y-txt = ",200a1)')
     &            ( itayt(m)(i:i),i = 1, itayl(m) )

               end if

*-----------------------------------------------------------------------

               if( itazl(m) .gt. 0 ) then

                  write(iot,'("    z-txt = ",200a1)')
     &            ( itazt(m)(i:i),i = 1, itazl(m) )

               end if

*-----------------------------------------------------------------------

               if( itgsh(m) .gt. 0 ) then

                  write(iot,'("    gshow = ",i4,10x,
     &            "  # 0: no 1:bnd, 2:bnd+mat, ",
     &     "3:bnd+reg 4:bnd+lat 5:bmp style")')
     &            itgsh(m)

               end if

*-----------------------------------------------------------------------

               if( itglt(m) .ge. 0 ) then

                  write(iot,'("    gslat = ",i4,10x,
     &            "  # 0: no line for lat, 1: with,",
     &               " 2: no line between the same cell")')
     &            itglt(m)

               end if

*-----------------------------------------------------------------------

               if( itgsh(m) .gt. 0 .and. itger(m) .gt. 0 ) then

                  write(iot,'("    ginfo = ",i4,11x,
     &            " # (D=0) 0:no, 1:show in graph, 2:add file")')
     &            itger(m)

               end if

*-----------------------------------------------------------------------

               if( itrsh(m) .gt. 0 ) then

                  write(iot,'("    rshow = ",i4,10x,
     &            "  # 0: no 1:bnd, 2:bnd+mat, ",
     &     "3:bnd+reg 4:bnd+lat 5:bmp style")')
     &            itrsh(m)

                  call echmty(1,iot,'x',itxty(m),itxnm(m),
     &                        rtxdl(m),rtxmi(m),rtxma(m),itxrg(m)
     &                        ,abs(itxnm(m)+1),das_itxrg)
                  call echmty(1,iot,'y',ityty(m),itynm(m),
     &                        rtydl(m),rtymi(m),rtyma(m),ityrg(m)
     &                        ,abs(itynm(m)+1),das_ityrg)
                  call echmty(1,iot,'z',itzty(m),itznm(m),
     &                        rtzdl(m),rtzmi(m),rtzma(m),itzrg(m)
     &                        ,abs(itznm(m)+1),das_itzrg)

               end if

*-----------------------------------------------------------------------

               if( itgsh(m) .gt. 0 .or. itrsh(m) .gt. 0 ) then

                  write(iot,'("    resol = ",i4,10x,
     &            "  # (D=1) resolution of gshow or rshow")')
     &            itres(m)

                  write(iot,'("    width =",1p1g14.7,1x,1x,
     &            " # (D=0.5) width of lines for gshow or rshow")')
     &            rtwid(m)

               end if

*-----------------------------------------------------------------------

               if( iterl(m) .ne. 72 ) then

                  write(iot,'("   iechrl =",i5,10x,
     &            "  # (D=72) length reg of val or vol")')
     &            iterl(m)

               end if

*-----------------------------------------------------------------------

               if( iteps(m) .ne. 0 ) then

                  write(iot,'("   epsout =",i5,10x,
     &            "  # (D=0) generate eps file by ANGEL")')
     &            iteps(m)

               end if

*-----------------------------------------------------------------------

               if( itbmp(m) .ne. 0 ) then

                  write(iot,'("   bmpout =",i5,10x,
     &            "  # (D=0) generate bitmap file")')
     &            itbmp(m)

               end if

*-----------------------------------------------------------------------

               if( itvtk(m) .ne. 0 ) then

                  write(iot,'("   vtkout =",i5,10x,
     &            "  # (D=0) generate vtk file")')
     &            itvtk(m)

                  if ( itvtkfmt(m) .eq. 0 ) then

                     write(iot,'("   vtkfmt =",i5,10x,
     &               "  # (D=0) generate vtk ascii file")')
     &               itvtkfmt(m)

                  else

                     write(iot,'("   vtkfmt =",i5,10x,
     &               "  # (D=0) generate vtk binary file")')
     &               itvtkfmt(m)

                  end if

               end if

*-----------------------------------------------------------------------
*        check of counter
*-----------------------------------------------------------------------

            do i = 1, 3

               if( itcnt(i,m) .ne. 0 ) then

                  write(iot,'(" ctmin(",i1,") =",i5,10x,
     &            "  # (D=-9999) min. of counter ",i1)')
     &            i, itcnt(i*2+2,m), i

                  write(iot,'(" ctmax(",i1,") =",i5,10x,
     &            "  # (D= 9999) max. of counter ",i1)')
     &            i, itcnt(i*2+3,m), i

                  if( ncntc(i) .eq. 0 )
     &            write(iot,'("# Warning: there is no counter = ",
     &                        i1," definition in [counter].")') i

               end if

            end do

            do i = 1, 3
               if( itcnth(i,m) .ne. 0 ) then
                  write(iot,'(" chmin(",i1,") =",i5,10x,
     &            "  # (D=-9999) min. of history counter ",i1)')
     &            i, itcnth(i*2+2,m), i
                  write(iot,'(" chmax(",i1,") =",i5,10x,
     &            "  # (D= 9999) max. of history counter ",i1)')
     &            i, itcnth(i*2+3,m), i
               end if
            end do

*-----------------------------------------------------------------------

               if( rtstd(m) .gt. 0.d0 ) then

                  write(iot,'("   stdcut =",1p1g14.7,1x,1x,
     &            " # (D=-1) threshold value of STD cutoff")')
     &            rtstd(m)

               end if

*-----------------------------------------------------------------------
               if( itfoam(m) .ne. 0 ) then

                  write(iot,'("  foamout =",i5,10x,
     &            "  # (D=0) generate OpenFOAM data")')
     &            itfoam(m)

               end if

*-----------------------------------------------------------------------
               if( itnzn(m) .ne. 0 ) then

                  write(iot,'(" mxnuclei =",i5,10x,
     &            "  # (D=3000) maximum of product nuclei.",
     &            " 0: unlimited")') itnzn(m)-1

               end if

*-----------------------------------------------------------------------

*-----------------------------------------------------------------------
*        sumtally subsection
*-----------------------------------------------------------------------
            if ( icntl==13 .and. nsumtalRead(m)==3 ) then

               call sumtal_echo(iot,m)

            end if
*-----------------------------------------------------------------------

            write(iot,'("#    used :",
     &                  "        main (  %)",
     &                  "        temp (  %)",
     &                  "       total (  %)"/
     &                  "#  memory :",
     &                  3(i12," (",i3,")"))')
     &                  itsmn(m), nint(dble(itsmn(m))/dble(mdas)*100),
     &                  itstm(m), nint(dble(itstm(m))/dble(mdas)*100),
     &                  itsmn(m)+itstm(m),
     &                  nint(dble(itsmn(m)+itstm(m))/dble(mdas)*100)

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine thetech(iot,m,iaxs,icax)
*                                                                      *
*       input echo of the heat tally                                   *
*       last modified by K.Niita on 2004/09/15                         *
*                                                                      *
************************************************************************

      use sumtallymod
      use sangelmod, only: itsans
      use partmod ! frtati 2021/10/05
      use moddas_mesh
      use moddas_region
      use moddas_tally
*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'

*-----------------------------------------------------------------------

      common /tall00/ itnm, ital(itlmax), itals(itlmax), italm(itlmax)
      common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                itout(itlmax), ittwo(itlmax)

      common /tall02/ itrgn(itlmax), itrgm(itlmax), itreg(itlmax)
      common /tall24/ itrcm(itlmax), itrcg(itlmax)
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

      common /tall14/ ittll(itlmax), ittle(itlmax)
      character ittle*80

      common /tall15/ itanl(itlmax), itang(itlmax)
      character itang*200

      common /tall27/ itaxl(itlmax), itaxt(itlmax),
     &                itayl(itlmax), itayt(itlmax),
     &                itazl(itlmax), itazt(itlmax)
      character itaxt*200, itayt*200, itazt*200

      common /tall17/ itndy(itlmax)

      common /tall19/ itsmn(itlmax), itstm(itlmax)

      common /tall21/ rtfac(itlmax)
      common /tall22/ itmcn(itlmax), itmtn(itlmax), itmtt(itlmax)

      common /tall28/ itgsh(itlmax)
      common /tall29/ itrsh(itlmax)
      common /tall31/ iterl(itlmax)

      common /tall32/ itelc(itlmax)
      common /tall34/ itvm(itlmax)
      common /tall35/ iteps(itlmax)
      common /tall36/ itdpo(itlmax)
      common /tall37/ itcnt(9,itlmax)
      common /cntreg/ ncntc(3), ncreg(3), incrc(3), incrt(3), kcont(3)
      common /cntpat/ icpan(3), icpat(3,20,2)

      common /tall38/ itres(itlmax)
      common /tall39/ rtwid(itlmax)
      common /tall46/ itmtr(itlmax,4), rtmtr(itlmax,13)
      common /tall49/ itglt(itlmax)

      common /tall58/ iterr(itlmax)

      common /tall60/ itger(itlmax)
      common /tall63/ itbmp(itlmax)
      common /tall64/ itvtk(itlmax),itvtkfmt(itlmax)

      common /tall67/ itstd(itlmax), rtstd(itlmax)

      common /tall82/ itcnth(9,itlmax)

      dimension idas(mdas*2) !20220906frtati temporary for fbounds-check
      equivalence ( das, idas )

*-----------------------------------------------------------------------

      common /inpec/ ititl, ipara, ibody, iregn, llarr, itby, itar
      common /kmat1d/ idmn(0:kvlmax), idnm(kvmmax)

*-----------------------------------------------------------------------
       integer   icntl
       integer   inucr

       common /tcntl/  icntl, inucr
!-----------------------------------------------------------------------

      character asfil*100

      character aname(10)*3
      data aname / 'eng','reg','  x','  y','  z','  r',
     &             ' xy',' yz',' xz',' rz'/

      dimension ian(6)
      dimension jmat(20)

*-----------------------------------------------------------------------
      common /tall79/ itallech
      if( itallech .eq. 0 ) return

*-----------------------------------------------------------------------
*     input echo
*-----------------------------------------------------------------------

            write(iot,'("[ T-Heat ]")')

                  write(iot,'("    title = ",80a1)')
     &            ( ittle(m)(i:i), i = 1, ittll(m) )

            if( itmsh(m) .eq. 1 ) then

                  write(iot,'("     mesh =  reg",11x,
     &            " # mesh type is region-wise")')

                     idas1 = mmmax
                     idas2 = ( idas1 + itrgn(m) - 1 ) * 2 + 1
                     idas3 = idas2 + itrgn(m)

                  call echrg(1,iterl(m),0,1,iot,itrcm(m),
     &                       idas_itrcg(itrcg(m)),
     &                       itrgn(m),itrgm(m),idas_itreg(itreg(m)),
     &                       das(idas1),idas(idas2),
     &                       itrnv(m),idas(itriv(m)),das(itrrv(m)),
     &                       idas3)

            else if( itmsh(m) .eq. 2 ) then

                  write(iot,'("     mesh =  r-z",11x,
     &            " # mesh type is r-z scoring mesh")')

                  write(iot,'("       x0 = ",1p1g14.7,1x,
     &            " # center x-position of r-z mesh")') rtrx0(m)

                  write(iot,'("       y0 = ",1p1g14.7,1x,
     &            " # center y-position of r-z mesh")') rtry0(m)

                  call echmty(0,iot,'r',itrty(m),itrnm(m),
     &                        rtrdl(m),rtrmi(m),rtrma(m),itrrg(m)
     &                        ,abs(itrnm(m)+1),das_itrrg)
                  call echmty(0,iot,'z',itzty(m),itznm(m),
     &                        rtzdl(m),rtzmi(m),rtzma(m),itzrg(m)
     &                        ,abs(itznm(m)+1),das_itzrg)

            else if( itmsh(m) .eq. 3 ) then

                  write(iot,'("     mesh =  xyz",11x,
     &            " # mesh type is xyz scoring mesh")')

                  call echmty(0,iot,'x',itxty(m),itxnm(m),
     &                        rtxdl(m),rtxmi(m),rtxma(m),itxrg(m)
     &                        ,abs(itxnm(m)+1),das_itxrg)
                  call echmty(0,iot,'y',ityty(m),itynm(m),
     &                        rtydl(m),rtymi(m),rtyma(m),ityrg(m)
     &                        ,abs(itynm(m)+1),das_ityrg)
                  call echmty(0,iot,'z',itzty(m),itznm(m),
     &                        rtzdl(m),rtzmi(m),rtzma(m),itzrg(m)
     &                        ,abs(itznm(m)+1),das_itzrg)

            end if

*-----------------------------------------------------------------------

            if( itmtr(m,1) .gt. 0 ) then

               if( itmtr(m,1) .eq. 1 ) then

                     write(iot,'("     trcl = ",i6,9x,
     &               " # transform ID")') itmtr(m,3)

               else

                  if( itmtr(m,2) .eq. 0 ) then

                     write(iot,'("     trcl = ",
     &                           3(1p1e15.7),
     &                           3(/12x,3(1p1e15.7)),i5)')
     &               ( rtmtr(m,j), j = 1, 12 ),
     &               nint( rtmtr(m,13) )

                  else

                     write(iot,'("    *trcl = ",
     &                           3(1p1e15.7),
     &                           3(/12x,3(1p1e15.7)),i5)')
     &               ( rtmtr(m,j), j = 1, 12 ),
     &               nint( rtmtr(m,13) )

                  end if

               end if

            end if

*-----------------------------------------------------------------------
               if( itunt(m) .eq. 0 ) then

                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [Gy/source]")') itunt(m)

               else if( itunt(m) .eq. 1 ) then

                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [MeV/cm^3/source]")') itunt(m)

               else if( itunt(m) .eq. 2 ) then

                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [MeV/source]")') itunt(m)

               else if( itunt(m) .eq. 3 ) then

                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [1/source]",
     &            " : only for output=deposit")') itunt(m)

               end if

                  asfil = '  # file name of output for the above axis'

            if( icax .eq. 0 ) then

                  icn = itaxn(m)

               do i = 1, icn

                  ian(i) = i

               end  do

            else

                  icn = 1
                  ian(1) = iaxs

            end if

*-----------------------------------------------------------------------

               idtyp = 0

            do i = 1, icn

               if( itaxs(m,ian(i)) .ge. 7 .and.
     &             itaxs(m,ian(i)) .le. 10 )  idtyp = ittwo(m)

            end do

            if( idtyp .ne. 0 ) then

               write(iot,'("  2D-type = ",i4,11x,
     &         " # 1:Cont, 2:Clust, 3:Color, 4:xyz, 5:mat,",
     &         " 6:Clust+Cont, 7:Col+Cont")') idtyp

            end if

*-----------------------------------------------------------------------

            do i = 1, icn

                  iax = ian(i)

                  write(iot,'("     axis =  ",a3,11x,
     &            " # axis of output")')
     &            aname( itaxs(m,iax) )

             if ( itfll(m,iax) .gt. 0 ) then

                  msfile = max( 14, itfll(m,iax) )

                  write(iot,'("     file = ",100a1)')
     &                  ( ctfln(m,iax)(j:j), j = 1, msfile ),
     &                  ( asfil(j:j), j = 1, 42 )

             end if

            end do

*-----------------------------------------------------------------------

               if( itmtn(m) .eq. 0 ) then

                  write(iot,'(" material =  all"11x,
     &            " # (D=all) number of specific material")')

               else if( itmtn(m) .gt. 0 ) then

                  write(iot,'(" material =",i5,11x,
     &            " # (D=all) number of specific material")')
     &            itmtn(m) * itmcn(m)

                     k = 0

                  do j = 1, itmtn(m)

                     k = k + 1

                     jmat(k) = ismte( itmtt(m) + j - 1 )

                     if( k / 10 * 10 .eq. k .or. j .eq. itmtn(m) ) then

                        write(iot,'(11x,10(i6))') ( jmat(i),i = 1, k )

                        k = 0

                     end if

                  end do

               end if

*-----------------------------------------------------------------------

               if( itout(m) .eq. 1 ) then

                  write(iot,'("   output =  all",10x,
     &            "  # all information is written")')

               else if( itout(m) .eq. 2 ) then

                  write(iot,'("   output = simple",8x,
     &            "  # simple essential information is written")')

               else if( itout(m) .eq. 3 ) then

                  write(iot,'("   output = heat",10x,
     &            "  # only heat is written")')

               else if( itout(m) .eq. 4 ) then

                  write(iot,'("   output = deposit_all",3x,
     &            "  # all information is written")')

               else if( itout(m) .eq. 5 ) then

                  write(iot,'("   output = deposit_simple",
     &            "  # simple essential information is written")')

               else if( itout(m) .eq. 6 ) then

                  write(iot,'("   output = deposit_heat",2x,
     &            "  # only deposit energy is written")')

               end if

*-----------------------------------------------------------------------

               if( itout(m) .ge. 4 ) then

                  write(iot,'("  deposit = ",i4,11x,
     &            " # (D=0) 0-> total deposit dist,",
     &                    " 1-> each process")')
     &            itdpo(m)

                  call echmty(0,iot,'e',itety(m),itenm(m),
     &                        rtedl(m),rtemi(m),rtema(m),iterg(m)
     &                        ,abs(itenm(m)+1),das_iterg)

               end if

*-----------------------------------------------------------------------

               if( itpan(m) .gt. 0 ) then

               call echprt(1,iot,m,itpan,itpat,jtpat,itnm,itmxpt,itmxgp)

               end if

*-----------------------------------------------------------------------

               if( rtfac(m) .ne. 1.0d0 ) then

                  write(iot,'("   factor =",1p1g14.7,1x,1x,
     &            " # (D=1.0) normalization factor")')
     &            rtfac(m)

               end if

*-----------------------------------------------------------------------

                  write(iot,'(" electron = ",i4,11x,
     &            " # (D=0) 0-> photon library,",
     &                    " 1-> electron ionization")')
     &            itelc(m)

*-----------------------------------------------------------------------

               if( itmtn(m) .gt. 0 .and. itmsh(m) .eq. 3 ) then

                  write(iot,'("   volmat = ",i4,11x,
     &            " # (D=9) vol.corr., 0:no, number of scan / dir.")')
     &            itvm(m)

               end if

*-----------------------------------------------------------------------

               if( itanl(m) .gt. 0 ) then

                  write(iot,'("    angel = ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               end if

*-----------------------------------------------------------------------

               if( itsans(m) .gt. 0 ) then

                  write(iot,'("   sangel = ",i2)') itsans(m)
                  call write_sangel(iot,m,1)

               end if

*-----------------------------------------------------------------------

               if( itaxl(m) .gt. 0 ) then

                  write(iot,'("    x-txt = ",200a1)')
     &            ( itaxt(m)(i:i),i = 1, itaxl(m) )

               end if

*-----------------------------------------------------------------------

               if( itayl(m) .gt. 0 ) then

                  write(iot,'("    y-txt = ",200a1)')
     &            ( itayt(m)(i:i),i = 1, itayl(m) )

               end if

*-----------------------------------------------------------------------

               if( itazl(m) .gt. 0 ) then

                  write(iot,'("    z-txt = ",200a1)')
     &            ( itazt(m)(i:i),i = 1, itazl(m) )

               end if

*-----------------------------------------------------------------------

               if( itgsh(m) .gt. 0 ) then

                  write(iot,'("    gshow = ",i4,10x,
     &            "  # 0: no 1:bnd, 2:bnd+mat, ",
     &     "3:bnd+reg 4:bnd+lat 5:bmp style")')
     &            itgsh(m)

               end if

*-----------------------------------------------------------------------

               if( itglt(m) .ge. 0 ) then

                  write(iot,'("    gslat = ",i4,10x,
     &            "  # 0: no line for lat, 1: with,",
     &               " 2: no line between the same cell")')
     &            itglt(m)

               end if

*-----------------------------------------------------------------------

               if( itgsh(m) .gt. 0 .and. itger(m) .gt. 0 ) then

                  write(iot,'("    ginfo = ",i4,11x,
     &            " # (D=0) 0:no, 1:show in graph, 2:add file")')
     &            itger(m)

               end if

*-----------------------------------------------------------------------

               if( itrsh(m) .gt. 0 ) then

                  write(iot,'("    rshow = ",i4,10x,
     &            "  # 0: no 1:bnd, 2:bnd+mat, ",
     &     "3:bnd+reg 4:bnd+lat 5:bmp style")')
     &            itrsh(m)

                  call echmty(1,iot,'x',itxty(m),itxnm(m),
     &                        rtxdl(m),rtxmi(m),rtxma(m),itxrg(m)
     &                        ,abs(itxnm(m)+1),das_itxrg)
                  call echmty(1,iot,'y',ityty(m),itynm(m),
     &                        rtydl(m),rtymi(m),rtyma(m),ityrg(m)
     &                        ,abs(itynm(m)+1),das_ityrg)
                  call echmty(1,iot,'z',itzty(m),itznm(m),
     &                        rtzdl(m),rtzmi(m),rtzma(m),itzrg(m)
     &                        ,abs(itznm(m)+1),das_itzrg)

               end if

*-----------------------------------------------------------------------

               if( itgsh(m) .gt. 0 .or. itrsh(m) .gt. 0 ) then

                  write(iot,'("    resol = ",i4,10x,
     &            "  # (D=1) resolution of gshow or rshow")')
     &            itres(m)

                  write(iot,'("    width =",1p1g14.7,1x,1x,
     &            " # (D=0.5) width of lines for gshow or rshow")')
     &            rtwid(m)

               end if

*-----------------------------------------------------------------------

               if( iterl(m) .ne. 72 ) then

                  write(iot,'("   iechrl =",i5,10x,
     &            "  # (D=72) length reg of val or vol")')
     &            iterl(m)

               end if

*-----------------------------------------------------------------------

               if( iteps(m) .ne. 0 ) then

                  write(iot,'("   epsout =",i5,10x,
     &            "  # (D=0) generate eps file by ANGEL")')
     &            iteps(m)

               end if

*-----------------------------------------------------------------------

               if( itbmp(m) .ne. 0 ) then

                  write(iot,'("   bmpout =",i5,10x,
     &            "  # (D=0) generate bitmap file")')
     &            itbmp(m)

               end if

*-----------------------------------------------------------------------

               if( itvtk(m) .ne. 0 ) then

                  write(iot,'("   vtkout =",i5,10x,
     &            "  # (D=0) generate vtk file")')
     &            itvtk(m)

                  if ( itvtkfmt(m) .eq. 0 ) then

                     write(iot,'("   vtkfmt =",i5,10x,
     &               "  # (D=0) generate vtk ascii file")')
     &               itvtkfmt(m)

                  else

                     write(iot,'("   vtkfmt =",i5,10x,
     &               "  # (D=0) generate vtk binary file")')
     &               itvtkfmt(m)

                  end if

               end if

*-----------------------------------------------------------------------

               if( iterr(m) .ne. 0 ) then

                  write(iot,'("    err2d =",i5,10x,
     &            "  # (D=0) print errors in 2d")')
     &            iterr(m)

               end if

*-----------------------------------------------------------------------
*        check of counter
*-----------------------------------------------------------------------

            do i = 1, 3

               if( itcnt(i,m) .ne. 0 ) then

                  write(iot,'(" ctmin(",i1,") =",i5,10x,
     &            "  # (D=-9999) min. of counter ",i1)')
     &            i, itcnt(i*2+2,m), i

                  write(iot,'(" ctmax(",i1,") =",i5,10x,
     &            "  # (D= 9999) max. of counter ",i1)')
     &            i, itcnt(i*2+3,m), i

                  if( ncntc(i) .eq. 0 )
     &            write(iot,'("# Warning: there is no counter = ",
     &                        i1," definition in [counter].")') i

               end if

            end do

            do i = 1, 3
               if( itcnth(i,m) .ne. 0 ) then
                  write(iot,'(" chmin(",i1,") =",i5,10x,
     &            "  # (D=-9999) min. of history counter ",i1)')
     &            i, itcnth(i*2+2,m), i
                  write(iot,'(" chmax('',i1,'') =",i5,10x,
     &            "  # (D= 9999) max. of history counter ",i1)')
     &            i, itcnth(i*2+3,m), i
               end if
            end do

*-----------------------------------------------------------------------

               if( rtstd(m) .gt. 0.d0 ) then

                  write(iot,'("   stdcut =",1p1g14.7,1x,1x,
     &            " # (D=-1) threshold value of STD cutoff")')
     &            rtstd(m)

               end if

*-----------------------------------------------------------------------
               if( itmxang(m) .ne. 0 ) then

                  write(iot,'(" maxangel =",i5,10x,
     &            "  # (D=0) number of output particles")')
     &            itmxang(m)

               end if

*-----------------------------------------------------------------------
*        sumtally subsection
*-----------------------------------------------------------------------
            if ( icntl==13 .and. nsumtalRead(m)==4 ) then

               call sumtal_echo(iot,m)

            end if
*-----------------------------------------------------------------------

            write(iot,'("#    used :",
     &                  "        main (  %)",
     &                  "        temp (  %)",
     &                  "       total (  %)"/
     &                  "#  memory :",
     &                  3(i12," (",i3,")"))')
     &                  itsmn(m), nint(dble(itsmn(m))/dble(mdas)*100),
     &                  itstm(m), nint(dble(itstm(m))/dble(mdas)*100),
     &                  itsmn(m)+itstm(m),
     &                  nint(dble(itsmn(m)+itstm(m))/dble(mdas)*100)

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      subroutine tstaech(iot,m,iaxs,icax)
*                                                                      *
*       input echo of the star tally                                   *
*       last modified by K.Niita on 2008/01/29                         *
*                                                                      *
************************************************************************

      use sumtallymod
      use sangelmod, only: itsans
      use partmod ! frtati 2021/10/05
      use moddas_mesh
      use moddas_region
      use moddas_tally
*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'

*-----------------------------------------------------------------------

      common /tall00/ itnm, ital(itlmax), itals(itlmax), italm(itlmax)
      common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                itout(itlmax), ittwo(itlmax)

      common /tall02/ itrgn(itlmax), itrgm(itlmax), itreg(itlmax)
      common /tall24/ itrcm(itlmax), itrcg(itlmax)
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

      common /tall16/ itnfr(itlmax), itnfz(itlmax), itnfn(itlmax)

      common /tall14/ ittll(itlmax), ittle(itlmax)
      character ittle*80

      common /tall15/ itanl(itlmax), itang(itlmax)
      character itang*200

      common /tall27/ itaxl(itlmax), itaxt(itlmax),
     &                itayl(itlmax), itayt(itlmax),
     &                itazl(itlmax), itazt(itlmax)
      character itaxt*200, itayt*200, itazt*200

      common /tall19/ itsmn(itlmax), itstm(itlmax)

      common /tall20/ ittty(itlmax), ittnm(itlmax), ittrg(itlmax),
     &                rttmi(itlmax), rttma(itlmax), rttdl(itlmax)
      common /tall21/ rtfac(itlmax)
      common /tall22/ itmcn(itlmax), itmtn(itlmax), itmtt(itlmax)
      common /tall28/ itgsh(itlmax)
      common /tall29/ itrsh(itlmax)
      common /tall31/ iterl(itlmax)
      common /tall34/ itvm(itlmax)
      common /tall35/ iteps(itlmax)
      common /tall37/ itcnt(9,itlmax)

      common /cntreg/ ncntc(3), ncreg(3), incrc(3), incrt(3), kcont(3)
      common /cntpat/ icpan(3), icpat(3,20,2)

      common /tall38/ itres(itlmax)
      common /tall39/ rtwid(itlmax)
      common /tall46/ itmtr(itlmax,4), rtmtr(itlmax,13)
      common /tall48/ itmdp(itlmax,0:30), itmdf(itlmax)
      common /tall49/ itglt(itlmax)

      common /tall60/ itger(itlmax)
      common /tall63/ itbmp(itlmax)
      common /tall64/ itvtk(itlmax),itvtkfmt(itlmax)

      common /tall67/ itstd(itlmax), rtstd(itlmax)
      common /tall75/ itnlatcel(itlmax), itnlatmem(itlmax)
      common /tall82/ itcnth(9,itlmax)

      dimension idas(mdas*2) !20220906frtati temporary for fbounds-check
      equivalence ( das, idas )

*-----------------------------------------------------------------------

      common /inpec/ ititl, ipara, ibody, iregn, llarr, itby, itar
      common /kmat1d/ idmn(0:kvlmax), idnm(kvmmax)

*-----------------------------------------------------------------------
       integer   icntl
       integer   inucr

       common /tcntl/  icntl, inucr
!-----------------------------------------------------------------------

      character asfil*100
      character bsfil*100
      character csfil*100
      character dsfil*100

      character aname(12)*3     ! S.Abe 2018/02/15, 11->12
      data aname / 'eng','reg','  x','  y','  z','  r',
     &             ' xy',' yz',' xz',' rz','  t','act'/

      character dmpc(30)*4
      data dmpc / '  kf','   x','   y','   z','   u','   v','   w',
     &            '   e','  wt','  tm','  c1','  c2','  c3',
     &            '  sx','  sy','  sz','  n0','  nc','  nb','  no',
     &            '    ','    ','    ','    ','    ','    ',
     &            '    ','    ','    ','    '/

      dimension ian(6)
      dimension jmat(20)
      character chau*8
      character chaus(7)*8
      character elmnt(104)*3

*-----------------------------------------------------------------------

      data elmnt /
     &    'H  ','He ','Li ','Be ','B  ','C  ','N  ','O  ','F  ','Ne ',
     &    'Na ','Mg ','Al ','Si ','P  ','S  ','Cl ','Ar ','K  ','Ca ',
     &    'Sc ','Ti ','V  ','Cr ','Mn ','Fe ','Co ','Ni ','Cu ','Zn ',
     &    'Ga ','Ge ','As ','Se ','Br ','Kr ','Rb ','Sr ','Y  ','Zr ',
     &    'Nb ','Mo ','Tc ','Ru ','Rh ','Pd ','Ag ','Cd ','In ','Sn ',
     &    'Sb ','Te ','I  ','Xe ','Cs ','Ba ','La ','Ce ','Pr ','Nd ',
     &    'Pm ','Sm ','Eu ','Gd ','Tb ','Dy ','Ho ','Er ','Tm ','Yb ',
     &    'Lu ','Hf ','Ta ','W  ','Re ','Os ','Ir ','Pt ','Au ','Hg ',
     &    'Tl ','Pb ','Bi ','Po ','At ','Rn ','Fr ','Ra ','Ac ','Th ',
     &    'Pa ','U  ','Np ','Pu ','Am ','Cm ','Bk ','Cf ','Es ','Fm ',
     &    'Md ','No ','Lr ','Ku ' /

      integer           iMeVperu
      common /cMeVperu/ iMeVperu

      common /tall71/ itmorp(itlmax),itname(itlmax)
      common /tall72/ itactnm(itlmax), itactrg(itlmax), itactmax(itlmax)

*-----------------------------------------------------------------------
      common /tall79/ itallech
      if( itallech .eq. 0 ) return

*-----------------------------------------------------------------------
*     input echo
*-----------------------------------------------------------------------

      if(itname(m).eq.1) then
            write(iot,'("[ T-Star ]")')
      else
            write(iot,'("[ T-Interact ]")')
      endif
                  write(iot,'("    title = ",80a1)')
     &            ( ittle(m)(i:i), i = 1, ittll(m) )

*-----------------------------------------------------------------------

            if( itmsh(m) .eq. 1 ) then

                  write(iot,'("     mesh =  reg",11x,
     &            " # mesh type is region-wise")')

                     idas1 = mmmax
                     idas2 = ( idas1 + itrgn(m) - 1 ) * 2 + 1
                     idas3 = idas2 + itrgn(m)

                  call echrg(1,iterl(m),0,1,iot,itrcm(m),
     &                       idas_itrcg(itrcg(m)),
     &                       itrgn(m),itrgm(m),idas_itreg(itreg(m)),
     &                       das(idas1),idas(idas2),
     &                       itrnv(m),idas(itriv(m)),das(itrrv(m)),
     &                       idas3)

            else if( itmsh(m) .eq. 2 ) then

                  write(iot,'("     mesh =  r-z",11x,
     &            " # mesh type is r-z scoring mesh")')

                  write(iot,'("       x0 = ",1p1g14.7,1x,
     &            " # center x-position of r-z mesh")') rtrx0(m)

                  write(iot,'("       y0 = ",1p1g14.7,1x,
     &            " # center y-position of r-z mesh")') rtry0(m)

                  call echmty(0,iot,'r',itrty(m),itrnm(m),
     &                        rtrdl(m),rtrmi(m),rtrma(m),itrrg(m)
     &                        ,abs(itrnm(m)+1),das_itrrg)
                  call echmty(0,iot,'z',itzty(m),itznm(m),
     &                        rtzdl(m),rtzmi(m),rtzma(m),itzrg(m)
     &                        ,abs(itznm(m)+1),das_itzrg)

            else if( itmsh(m) .eq. 3 ) then

                  write(iot,'("     mesh =  xyz",11x,
     &            " # mesh type is xyz scoring mesh")')

                  call echmty(0,iot,'x',itxty(m),itxnm(m),
     &                        rtxdl(m),rtxmi(m),rtxma(m),itxrg(m)
     &                        ,abs(itxnm(m)+1),das_itxrg)
                  call echmty(0,iot,'y',ityty(m),itynm(m),
     &                        rtydl(m),rtymi(m),rtyma(m),ityrg(m)
     &                        ,abs(itynm(m)+1),das_ityrg)
                  call echmty(0,iot,'z',itzty(m),itznm(m),
     &                        rtzdl(m),rtzmi(m),rtzma(m),itzrg(m)
     &                        ,abs(itznm(m)+1),das_itzrg)

            end if

*-----------------------------------------------------------------------

            if( itmtr(m,1) .gt. 0 ) then

               if( itmtr(m,1) .eq. 1 ) then

                     write(iot,'("     trcl = ",i6,9x,
     &               " # transform ID")') itmtr(m,3)

               else

                  if( itmtr(m,2) .eq. 0 ) then

                     write(iot,'("     trcl = ",
     &                           3(1p1e15.7),
     &                           3(/12x,3(1p1e15.7)),i5)')
     &               ( rtmtr(m,j), j = 1, 12 ),
     &               nint( rtmtr(m,13) )

                  else

                     write(iot,'("    *trcl = ",
     &                           3(1p1e15.7),
     &                           3(/12x,3(1p1e15.7)),i5)')
     &               ( rtmtr(m,j), j = 1, 12 ),
     &               nint( rtmtr(m,13) )

                  end if

               end if

            end if

*-----------------------------------------------------------------------

                  call echmty(0,iot,'e',itety(m),itenm(m),
     &                        rtedl(m),rtemi(m),rtema(m),iterg(m)
     &                        ,abs(itenm(m)+1),das_iterg)

*-----------------------------------------------------------------------

               if( ittty(m) .ne. 0 ) then

                  call echmty(0,iot,'t',ittty(m),ittnm(m),
     &                        rttdl(m),rttmi(m),rttma(m),ittrg(m)
     &                        ,abs(ittnm(m)+1),das_ittrg)

               end if

*-----------------------------------------------------------------------

               if( itmorp(m) .eq. 0 ) then

                  write(iot,'("     morp =  mean"10x,
     &            " # (D=mean) output mode")')

               elseif( itmorp(m) .ne. 0 ) then

                  write(iot,'("     morp =  prob"10x,
     &            " # (D=mean) output mode")')

               endif

*-----------------------------------------------------------------------

               if( itactnm(m) .gt. 1 ) then

                  write(iot,'("   maxact = ",i8,7x,
     &            " # number of act-mesh points")') itactnm(m)-1

               endif

*-----------------------------------------------------------------------

               if( itmtn(m) .eq. 0 ) then

                  write(iot,'(" material =  all"11x,
     &            " # (D=all) number of specific material")')

               else if( itmtn(m) .gt. 0 ) then

                  write(iot,'(" material =",i5,11x,
     &            " # (D=all) number of specific material")')
     &            itmtn(m) * itmcn(m)

                     k = 0

                  do j = 1, itmtn(m)

                     k = k + 1

                     jmat(k) = ismte( itmtt(m) + j - 1 )

                     if( k / 10 * 10 .eq. k .or. j .eq. itmtn(m) ) then

                        write(iot,'(11x,10(i6))') ( jmat(i),i = 1, k )

                        k = 0

                     end if

                  end do

               end if

*-----------------------------------------------------------------------

               if( itman(m) .eq. 0 ) then

                  write(iot,'("   mother =  all"11x,
     &            " # (D=all) number of specific mother")')

               else if( itman(m) .gt. 0 ) then

                  write(iot,'("   mother = ",i4,11x,
     &            " # (D=all) number of specific mother")')
     &            itman(m) * itmct(m)

                     k = 0

                  do j = 1, itman(m)

                     k = k + 1

                     iaz = ismat( itmat(m) + j - 1 )

                     iz = iaz / 1000
                     ia = iaz - iz * 1000

                     if( ia .eq. 0 ) then

                        chaus(k) = elmnt(iz)

                     else

                        call chname(idum,ia,iz,chau)

                        chaus(k) = chau

                     end if

                     if( k / 7 * 7 .eq. k .or. j .eq. itman(m) ) then

                        write(iot,'(12x,7(a8))') ( chaus(i),i = 1, k )

                        k = 0

                     end if

                  end do

               end if

*-----------------------------------------------------------------------

               if( itunt(m) .eq. 1 ) then

                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [1/cm^3/source]")') itunt(m)

               else if( itunt(m) .eq. 2 ) then

                if ( iMeVperu.eq. 1 ) then
                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [1/cm^3/(MeV/n)/source]")') itunt(m)
                else
                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [1/cm^3/MeV/source]")') itunt(m)
                end if

               else if( itunt(m) .eq. 11 ) then

                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [1/cm^3/nsec/source]")') itunt(m)

               else if( itunt(m) .eq. 12 ) then

                if ( iMeVperu.eq. 1 ) then
                  write(iot,'("     unit = ",i4,11x,
     &           " # unit is [1/cm^3/(MeV/n)/nsec/source]")') itunt(m)
                else
                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [1/cm^3/MeV/nsec/source]")') itunt(m)
                end if

               end if

*-----------------------------------------------------------------------

                  asfil = '  # file name of output for the above axis'
                  bsfil = '  # file name of dumped data'
                  csfil = '  # file name of dumped data summary'

            if( icax .eq. 0 ) then

                  icn = itaxn(m)

               do i = 1, icn

                  ian(i) = i

               end  do

            else

                  icn = 1
                  ian(1) = iaxs

            end if

*-----------------------------------------------------------------------

               idtyp = 0

            do i = 1, icn

               if( itaxs(m,ian(i)) .ge. 7 .and.
     &             itaxs(m,ian(i)) .le. 10 )  idtyp = ittwo(m)

            end do

            if( idtyp .ne. 0 ) then

               write(iot,'("  2D-type = ",i4,11x,
     &         " # 1:Cont, 2:Clust, 3:Color, 4:xyz, 5:mat,",
     &         " 6:Clust+Cont, 7:Col+Cont")') idtyp

            end if

*-----------------------------------------------------------------------

            do i = 1, icn

                  iax = ian(i)

                  write(iot,'("     axis =  ",a3,11x,
     &            " # axis of output")')
     &            aname( itaxs(m,iax) )

             if ( itfll(m,iax) .gt. 0 ) then

               if( itmdp(m,0) .eq. 0 ) then

                  msfile = max( 14, itfll(m,iax) )

                  write(iot,'("     file = ",100a1)')
     &                  ( ctfln(m,iax)(j:j), j = 1, msfile ),
     &                  ( asfil(j:j), j = 1, 42 )

               else

C S.H. added for Dump-Restart on 2014/5/7
                  msfile = max( 14, itfll(m,iax) )

                  write(iot,'("     file = ",100a1)')
     &                  ( ctfln(m,iax)(j:j), j = 1, msfile ),
     &                  ( csfil(j:j), j = 1, 36 )

                  msfile = max( 14, itfll(m,iax)+4 )

                  do k = 1, msfile
                     dsfil(k:k) = ' '
                  end do

                  do is = itfll(m,iax), 1, -1
                     if ( ctfln(m,iax)(is:is) .eq. '.' ) exit
                  end do
                  if ( is .eq. 0 ) is = itfll(m,iax)+1

                  dsfil(1:is-1) = ctfln(m,iax)(1:is-1)
                  dsfil(is:is+3) = '_dmp'
                  if ( is .lt. itfll(m,iax) ) then
                     dsfil(is+4:itfll(m,iax)+4)
     &                    = ctfln(m,iax)(is:itfll(m,iax))
                  end if

                  write(iot,'("#    file = ",100a1)')
     &                  ( dsfil(j:j), j = 1, msfile ),
     &                  ( bsfil(j:j), j = 1, 28 )

               end if

             end if

            end do

*-----------------------------------------------------------------------

               if( itout(m) .eq. 1 ) then

                  write(iot,'("   output =  all",10x,
     &            "  # (D=nuclear) decay, elastic and nuclear")')

               else if( itout(m) .eq. 2 ) then

                  write(iot,'("   output = decay",9x,
     &            "  # (D=nuclear) star density of decay reaction")')

               else if( itout(m) .eq. 3 ) then

                  write(iot,'("   output = elastic",7x,
     &           "  # (D=nuclear) star density of elastic reaction")')

               else if( itout(m) .eq. 4 ) then

                  write(iot,'("   output = nuclear",7x,
     &           "  # (D=nuclear) star density of nuclear reaction")')

               else if( itout(m) .eq. 5 ) then

                  write(iot,'("   output = fission",7x,
     &           "  # (D=nuclear) star density of fission reaction")')

               else if( itout(m) .eq. 6 ) then

                  write(iot,'("   output = absorption",4x,
     &            "  # (D=nuclear) absorption reaction")')

               else if( itout(m) .eq. 7 ) then

                  write(iot,'("   output = heavyion",6x,
     &            "  # (D=nuclear) heavy ion reaction")')

               else if( itout(m) .eq. 8 ) then

                  write(iot,'("   output = transmut",6x,
     &            "  # (D=nuclear) transmutation reaction")')

               else if( itout(m) .eq. 9 ) then

                  write(iot,'("   output = atomic",8x,
     &            "  # (D=nuclear) atomic reaction")')

               else if( itout(m) .eq. 19 ) then

                  write(iot,'("   output = deltaray",6x,
     &            "  # (D=nuclear) delta ray production")')

               else if( itout(m) .eq. 20 ) then

                  write(iot,'("   output = knockelec",5x,
     &            "  # (D=nuclear) knock-on electron production")')

               else if( itout(m) .eq. 21 ) then

                  write(iot,'("   output = atmflu",8x,
     &            "  # (D=nuclear) atomic fluorescence x-ray prod.")')

               else if( itout(m) .eq. 22 ) then

                  write(iot,'("   output = auger",9x,
     &            "  # (D=nuclear) auger electron prod.")')

               else if( itout(m) .eq. 23 ) then

                  write(iot,'("   output = brems",9x,
     &            "  # (D=nuclear) bremsstrahlung")')

               else if( itout(m) .eq. 24 ) then

                  write(iot,'("   output = photoelec",5x,
     &            "  # (D=nuclear) photoelectric effect")')

               else if( itout(m) .eq. 25 ) then

                  write(iot,'("   output = compton",7x,
     &            "  # (D=nuclear) compton scattering")')

               else if( itout(m) .eq. 26 ) then

                  write(iot,'("   output = pairprod",6x,
     &            "  # (D=nuclear) pair production")')

               else if( itout(m) .eq. 27 ) then

                  write(iot,'("   output = annih",9x,
     &            "  # (D=nuclear) positron annihilation")')
               else if( itout(m) .eq. 30 ) then

                  write(iot,'("   output = ets_elas",7x,
     &            "  # (D=nuclear) elastic scattering in ETS mode")')

               else if( itout(m) .eq. 31 ) then

                  write(iot,'("   output = ets_ioniz",7x,
     &            "  # (D=nuclear) ionization in ETS mode")')

               else if( itout(m) .eq. 32 ) then

                  write(iot,'("   output = ets_e-exc",7x,
     &            "  # (D=nuclear) electronic excitation",
     &            " in ETS mode")')

               else if( itout(m) .eq. 33 ) then

                  write(iot,'("   output = ets_dea",7x,
     &            "  # (D=nuclear) dissociative electron attachment",
     &            " in ETS mode")')

               else if( itout(m) .eq. 34 ) then

                  write(iot,'("   output = ets_v-exc",7x,
     &            "  # (D=nuclear) vibration exciation in ETS mode")')

               else if( itout(m) .eq. 35 ) then

                  write(iot,'("   output = ets_p-exc",7x,
     &            "  # (D=nuclear) phonon excitation in ETS mode")')

               else if( itout(m) .eq. 36 ) then

                  write(iot,'("   output = ets_r-exc",7x,
     &            "  # (D=nuclear) rotation excitation in ETS mode")')

               else if( itout(m) .eq. 37 ) then

                  write(iot,'("   output = ets_plasmon",7x,
     &            "  # (D=nuclear) plasmon excitation in ETS mode")')

               else if( itout(m) .eq. 38 ) then

                  write(iot,'("   output =",1x,
     &            "ets_ioniz_e-exc",7x,
     &            "  # (D=nuclear) ionization and",
     &            " electronic excitation in ETS mode")')

               else if( itout(m) .eq. 39 ) then

                  write(iot,'("   output =",1x,
     &            "ets_hit",7x,
     &            "  # (D=nuclear) all reactions in ETS mode")')


               else if( itout(m) .eq. 40 ) then

                  write(iot,'("   output = kurbuc_elas",7x,
     &            "  # (D=nuclear) elastic scattering in KURBUC mode")')

               else if( itout(m) .eq. 41 ) then

                  write(iot,'("   output = kurbuc_ioniz",7x,
     &            "  # (D=nuclear) ionization in KURBUC mode")')

               else if( itout(m) .eq. 42 ) then

                  write(iot,'("   output = kurbuc_e-exc",7x,
     &            "  # (D=nuclear) electronic excitation in KURBUC mode"
     &            )')

               else if( itout(m) .eq. 43 ) then

                  write(iot,'("   output = kurbuc_e-cap",7x,
     &            "  # (D=nuclear) electron capture in KURBUC mode")')

               else if( itout(m) .eq. 44 ) then

                  write(iot,'("   output = kurbuc_e-stp",7x,
     &            "  # (D=nuclear) electron stripping in KURBUC mode")')

               else if( itout(m) .eq. 48 ) then

                  write(iot,'("   output = kurbuc_ioniz_e-exc_e-cap",7x,
     &            "  # (D=nuclear) ionization + capture + strip in
     &               KURBUC mode")')

               else if( itout(m) .eq. 49 ) then

                  write(iot,'("   output = kurbuc_hit",7x,
     &            "  # (D=nuclear) all reactions in KURBUC mode")')

               else if( itout(m) .eq. 50 ) then

                  write(iot,'("   output = its_elas",7x,
     &            "  # (D=nuclear) elastic scattering in ITSART mode")')

               else if( itout(m) .eq. 51 ) then

                  write(iot,'("   output = its_ioniz",7x,
     &            "  # (D=nuclear) ionization in ITSART mode")')

               else if( itout(m) .eq. 52 ) then

                  write(iot,'("   output = its_e-exc",7x,
     &            "  # (D=nuclear) electronic excitation in ITSART mode"
     &               )')

               else if( itout(m) .eq. 58 ) then

                  write(iot,'("   output = its_ioniz_e-exc_e-cap",7x,
     &            "  # (D=nuclear) ionization, excitation and electron
     &               capture in ITSART mode")')

               else if( itout(m) .eq. 59 ) then

                  write(iot,'("   output = its_hit",7x,
     &            "  # (D=nuclear) all reaction in ITSART mode")')

               else if( itout(m) .eq. 96 ) then

                  write(iot,'("   output = ts_w",7x,
     &            "  # (D=nuclear) sum of the events that generate",
     &            " an electron in all track structure modes")')

               else if( itout(m) .eq. 97 ) then

                  write(iot,'("   output = ts_ioniz",7x,
     &            "  # (D=nuclear) ionization in all track structure",
     &            " modes")')

               else if( itout(m) .eq. 98 ) then

                  write(iot,'("   output = ts_bio",7x,
     &           "   # (D=nuclear) sum of the biological events (ioniz",
     &           ", e-exc, e-cap, dea) in all track structure modes")')

               else if( itout(m) .eq. 99 ) then

                  write(iot,'("   output = ts_hit",7x,
     &            "  # (D=nuclear) all reactions in track structure",
     &            " modes")')

               end if

*-----------------------------------------------------------------------
               call echprt(1,iot,m,itpan,itpat,jtpat,itnm,itmxpt,itmxgp)

*-----------------------------------------------------------------------

               if( rtfac(m) .ne. 1.0d0 ) then

                  write(iot,'("   factor =",1p1g14.7,1x,1x,
     &            " # (D=1.0) normalization factor")')
     &            rtfac(m)

               end if

*-----------------------------------------------------------------------

               if( itmtn(m) .gt. 0 .and. itmsh(m) .eq. 3 ) then

                  write(iot,'("   volmat = ",i4,11x,
     &            " # (D=9) vol.corr., 0:no, number of scan / dir.")')
     &            itvm(m)

               end if

*-----------------------------------------------------------------------

               if( itanl(m) .gt. 0 ) then

                  write(iot,'("    angel = ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               end if

*-----------------------------------------------------------------------

               if( itsans(m) .gt. 0 ) then

                  write(iot,'("   sangel = ",i2)') itsans(m)
                  call write_sangel(iot,m,1)

               end if

*-----------------------------------------------------------------------

               if( itaxl(m) .gt. 0 ) then

                  write(iot,'("    x-txt = ",200a1)')
     &            ( itaxt(m)(i:i),i = 1, itaxl(m) )

               end if

*-----------------------------------------------------------------------

               if( itayl(m) .gt. 0 ) then

                  write(iot,'("    y-txt = ",200a1)')
     &            ( itayt(m)(i:i),i = 1, itayl(m) )

               end if

*-----------------------------------------------------------------------

               if( itazl(m) .gt. 0 ) then

                  write(iot,'("    z-txt = ",200a1)')
     &            ( itazt(m)(i:i),i = 1, itazl(m) )

               end if

*-----------------------------------------------------------------------

               if( itgsh(m) .gt. 0 ) then

                  write(iot,'("    gshow = ",i4,10x,
     &            "  # 0: no 1:bnd, 2:bnd+mat, ",
     &     "3:bnd+reg 4:bnd+lat 5:bmp style")')
     &            itgsh(m)

               end if

*-----------------------------------------------------------------------

               if( itglt(m) .ge. 0 ) then

                  write(iot,'("    gslat = ",i4,10x,
     &            "  # 0: no line for lat, 1: with,",
     &               " 2: no line between the same cell")')
     &            itglt(m)

               end if

*-----------------------------------------------------------------------

               if( itgsh(m) .gt. 0 .and. itger(m) .gt. 0 ) then

                  write(iot,'("    ginfo = ",i4,11x,
     &            " # (D=0) 0:no, 1:show in graph, 2:add file")')
     &            itger(m)

               end if

*-----------------------------------------------------------------------

               if( itrsh(m) .gt. 0 ) then

                  write(iot,'("    rshow = ",i4,10x,
     &            "  # 0: no 1:bnd, 2:bnd+mat, ",
     &     "3:bnd+reg 4:bnd+lat 5:bmp style")')
     &            itrsh(m)

                  call echmty(1,iot,'x',itxty(m),itxnm(m),
     &                        rtxdl(m),rtxmi(m),rtxma(m),itxrg(m)
     &                        ,abs(itxnm(m)+1),das_itxrg)
                  call echmty(1,iot,'y',ityty(m),itynm(m),
     &                        rtydl(m),rtymi(m),rtyma(m),ityrg(m)
     &                        ,abs(itynm(m)+1),das_ityrg)
                  call echmty(1,iot,'z',itzty(m),itznm(m),
     &                        rtzdl(m),rtzmi(m),rtzma(m),itzrg(m)
     &                        ,abs(itznm(m)+1),das_itzrg)

               end if

*-----------------------------------------------------------------------

               if( itgsh(m) .gt. 0 .or. itrsh(m) .gt. 0 ) then

                  write(iot,'("    resol = ",i4,10x,
     &            "  # (D=1) resolution of gshow or rshow")')
     &            itres(m)

                  write(iot,'("    width =",1p1g14.7,1x,1x,
     &            " # (D=0.5) width of lines for gshow or rshow")')
     &            rtwid(m)

               end if

*-----------------------------------------------------------------------

               if( iterl(m) .ne. 72 ) then

                  write(iot,'("   iechrl =",i5,10x,
     &            "  # (D=72) length reg of val or vol")')
     &            iterl(m)

               end if

*-----------------------------------------------------------------------

               if( iteps(m) .ne. 0 ) then

                  write(iot,'("   epsout =",i5,10x,
     &            "  # (D=0) generate eps file by ANGEL")')
     &            iteps(m)

               end if

*-----------------------------------------------------------------------

               if( itbmp(m) .ne. 0 ) then

                  write(iot,'("   bmpout =",i5,10x,
     &            "  # (D=0) generate bitmap file")')
     &            itbmp(m)

               end if

*-----------------------------------------------------------------------

               if( itvtk(m) .ne. 0 ) then

                  write(iot,'("   vtkout =",i5,10x,
     &            "  # (D=0) generate vtk file")')
     &            itvtk(m)

                  if ( itvtkfmt(m) .eq. 0 ) then

                     write(iot,'("   vtkfmt =",i5,10x,
     &               "  # (D=0) generate vtk ascii file")')
     &               itvtkfmt(m)

                  else

                     write(iot,'("   vtkfmt =",i5,10x,
     &               "  # (D=0) generate vtk binary file")')
     &               itvtkfmt(m)

                  end if

               end if

*-----------------------------------------------------------------------

               if( itnlatcel(m) .ne. 0 ) then

                  write(iot,'("  nlatcel =",i5,10x,
     &            "  # (D=0) lattice cell ID for",
     &            " individual dealing in lattice")')
     &            itnlatcel(m)

                  write(iot,'("  nlatmem =",i5,10x,
     &            "  # (D=0) number of memory for",
     &            " individual dealing in lattice")')
     &            itnlatmem(m)

               end if


*-----------------------------------------------------------------------
*        check of counter
*-----------------------------------------------------------------------

            do i = 1, 3

               if( itcnt(i,m) .ne. 0 ) then

                  write(iot,'(" ctmin(",i1,") =",i5,10x,
     &            "  # (D=-9999) min. of counter ",i1)')
     &            i, itcnt(i*2+2,m), i

                  write(iot,'(" ctmax(",i1,") =",i5,10x,
     &            "  # (D= 9999) max. of counter ",i1)')
     &            i, itcnt(i*2+3,m), i

                  if( ncntc(i) .eq. 0 )
     &            write(iot,'("# Warning: there is no counter = ",
     &                        i1," definition in [counter].")') i

               end if

            end do

            do i = 1, 3
               if( itcnth(i,m) .ne. 0 ) then
                  write(iot,'(" chmin(",i1,") =",i5,10x,
     &            "  # (D=-9999) min. of history counter ",i1)')
     &            i, itcnth(i*2+2,m), i
                  write(iot,'(" chmax(",i1,") =",i5,10x,
     &            "  # (D= 9999) max. of history counter ",i1)')
     &            i, itcnth(i*2+3,m), i
               end if
            end do

*-----------------------------------------------------------------------

               if( rtstd(m) .gt. 0.d0 ) then

                  write(iot,'("   stdcut =",1p1g14.7,1x,1x,
     &            " # (D=-1) threshold value of STD cutoff")')
     &            rtstd(m)

               end if

*-----------------------------------------------------------------------
               if( itmxang(m) .ne. 0 ) then

                  write(iot,'(" maxangel =",i5,10x,
     &            "  # (D=0) number of output particles")')
     &            itmxang(m)

               end if

*-----------------------------------------------------------------------
*        dump the data
*-----------------------------------------------------------------------

         if( itmdp(m,0) .ne. 0 ) then

                  write(iot,'("     dump =",i5,11x,
     &            " # (D=0) number of dumped data,",
     &            " <0: ascii, >0: binary")') itmdp(m,0)

                  write(iot,'(13x,30(i4))')
     &               ( itmdp(m,j), j = 1, abs( itmdp(m,0) ) )

                  write(iot,'("# dump data  ",30(a4))')
     &               ( dmpc(itmdp(m,j)), j = 1, abs( itmdp(m,0) ) )

         end if

*-----------------------------------------------------------------------
*        sumtally subsection
*-----------------------------------------------------------------------
            if ( icntl==13 .and. nsumtalRead(m)==5 ) then

               call sumtal_echo(iot,m)

            end if
*-----------------------------------------------------------------------

            write(iot,'("#    used :",
     &                  "        main (  %)",
     &                  "        temp (  %)",
     &                  "       total (  %)"/
     &                  "#  memory :",
     &                  3(i12," (",i3,")"))')
     &                  itsmn(m), nint(dble(itsmn(m))/dble(mdas)*100),
     &                  itstm(m), nint(dble(itstm(m))/dble(mdas)*100),
     &                  itsmn(m)+itstm(m),
     &                  nint(dble(itsmn(m)+itstm(m))/dble(mdas)*100)

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine tproech(iot,m,iaxs,icax)
*                                                                      *
*       input echo of the product tally                                *
*       last modified by K.Niita on 2005/11/24                         *
*                                                                      *
************************************************************************

      use sumtallymod
      use sangelmod, only: itsans
      use partmod ! frtati 2021/10/05
      use moddas_mesh
      use moddas_region
      use moddas_tally
*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'

*-----------------------------------------------------------------------

      common /tall00/ itnm, ital(itlmax), itals(itlmax), italm(itlmax)
      common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                itout(itlmax), ittwo(itlmax)

      common /tall02/ itrgn(itlmax), itrgm(itlmax), itreg(itlmax)
      common /tall24/ itrcm(itlmax), itrcg(itlmax)
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

      common /tall16/ itnfr(itlmax), itnfz(itlmax), itnfn(itlmax)

      common /tall14/ ittll(itlmax), ittle(itlmax)
      character ittle*80

      common /tall15/ itanl(itlmax), itang(itlmax)
      character itang*200

      common /tall27/ itaxl(itlmax), itaxt(itlmax),
     &                itayl(itlmax), itayt(itlmax),
     &                itazl(itlmax), itazt(itlmax)
      character itaxt*200, itayt*200, itazt*200

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
      common /tall34/ itvm(itlmax)
      common /tall35/ iteps(itlmax)
      common /tall37/ itcnt(9,itlmax)

      common /cntreg/ ncntc(3), ncreg(3), incrc(3), incrt(3), kcont(3)
      common /cntpat/ icpan(3), icpat(3,20,2)

      common /tall38/ itres(itlmax)
      common /tall39/ rtwid(itlmax)
      common /tall46/ itmtr(itlmax,4), rtmtr(itlmax,13)
      common /tall48/ itmdp(itlmax,0:30), itmdf(itlmax)
      common /tall49/ itglt(itlmax)

      common /tall60/ itger(itlmax)
      common /tall63/ itbmp(itlmax)
      common /tall64/ itvtk(itlmax),itvtkfmt(itlmax)

      common /tall67/ itstd(itlmax), rtstd(itlmax)

      common /tall74/ iprim(itlmax)
      integer :: itfoam
      common /tall76/ itfoam(itlmax)

      common /tall82/ itcnth(9,itlmax)

ccse 2022/08/31 LET parameter
      common /tall50/ itlmt(itlmax), ite2l(itlmax)   ! CCSE 2022.08.31

      dimension idas(mdas*2) !20220906frtati temporary for fbounds-check
      equivalence ( das, idas )

*-----------------------------------------------------------------------

      common /inpec/ ititl, ipara, ibody, iregn, llarr, itby, itar
      common /kmat1d/ idmn(0:kvlmax), idnm(kvmmax)

*-----------------------------------------------------------------------
       integer   icntl
       integer   inucr

       common /tcntl/  icntl, inucr
!-----------------------------------------------------------------------

      character asfil*100
      character bsfil*100
      character csfil*100
      character dsfil*100
ccse 2022/08/31, add LET axis, 14 -> 15
      character aname(15)*3
      data aname / 'eng','reg','  x','  y','  z','  r',
     &             ' xy',' yz',' xz',' rz','  t','cos','the','tet',
     &             'let'/

      character dmpc(30)*4
      data dmpc / '  kf','   x','   y','   z','   u','   v','   w',
     &            '   e','  wt','  tm','  c1','  c2','  c3',
     &            '  sx','  sy','  sz','  n0','  nc','  nb','  no',
     &            '    ','    ','    ','    ','    ','    ',
     &            '    ','    ','    ','    '/

      dimension ian(6)
      dimension jmat(20)
      character chau*8
      character chaus(7)*8
      character elmnt(104)*3

*-----------------------------------------------------------------------

      data elmnt /
     &    'H  ','He ','Li ','Be ','B  ','C  ','N  ','O  ','F  ','Ne ',
     &    'Na ','Mg ','Al ','Si ','P  ','S  ','Cl ','Ar ','K  ','Ca ',
     &    'Sc ','Ti ','V  ','Cr ','Mn ','Fe ','Co ','Ni ','Cu ','Zn ',
     &    'Ga ','Ge ','As ','Se ','Br ','Kr ','Rb ','Sr ','Y  ','Zr ',
     &    'Nb ','Mo ','Tc ','Ru ','Rh ','Pd ','Ag ','Cd ','In ','Sn ',
     &    'Sb ','Te ','I  ','Xe ','Cs ','Ba ','La ','Ce ','Pr ','Nd ',
     &    'Pm ','Sm ','Eu ','Gd ','Tb ','Dy ','Ho ','Er ','Tm ','Yb ',
     &    'Lu ','Hf ','Ta ','W  ','Re ','Os ','Ir ','Pt ','Au ','Hg ',
     &    'Tl ','Pb ','Bi ','Po ','At ','Rn ','Fr ','Ra ','Ac ','Th ',
     &    'Pa ','U  ','Np ','Pu ','Am ','Cm ','Bk ','Cf ','Es ','Fm ',
     &    'Md ','No ','Lr ','Ku ' /

      integer           iMeVperu
      common /cMeVperu/ iMeVperu

*-----------------------------------------------------------------------
      common /tall79/ itallech
      if( itallech .eq. 0 ) return

*-----------------------------------------------------------------------
*     input echo
*-----------------------------------------------------------------------

            write(iot,'("[ T-Product ]")')

                  write(iot,'("    title = ",80a1)')
     &            ( ittle(m)(i:i), i = 1, ittll(m) )

*-----------------------------------------------------------------------

            if( itmsh(m) .eq. 1 ) then

                  write(iot,'("     mesh =  reg",11x,
     &            " # mesh type is region-wise")')

                     idas1 = mmmax
                     idas2 = ( idas1 + itrgn(m) - 1 ) * 2 + 1
                     idas3 = idas2 + itrgn(m)

                  call echrg(1,iterl(m),0,1,iot,itrcm(m),
     &                       idas_itrcg(itrcg(m)),
     &                       itrgn(m),itrgm(m),idas_itreg(itreg(m)),
     &                       das(idas1),idas(idas2),
     &                       itrnv(m),idas(itriv(m)),das(itrrv(m)),
     &                       idas3)

            else if( itmsh(m) .eq. 2 ) then

                  write(iot,'("     mesh =  r-z",11x,
     &            " # mesh type is r-z scoring mesh")')

                  write(iot,'("       x0 = ",1p1g14.7,1x,
     &            " # center x-position of r-z mesh")') rtrx0(m)

                  write(iot,'("       y0 = ",1p1g14.7,1x,
     &            " # center y-position of r-z mesh")') rtry0(m)

                  call echmty(0,iot,'r',itrty(m),itrnm(m),
     &                        rtrdl(m),rtrmi(m),rtrma(m),itrrg(m)
     &                        ,abs(itrnm(m)+1),das_itrrg)
                  call echmty(0,iot,'z',itzty(m),itznm(m),
     &                        rtzdl(m),rtzmi(m),rtzma(m),itzrg(m)
     &                        ,abs(itznm(m)+1),das_itzrg)

            else if( itmsh(m) .eq. 3 ) then

                  write(iot,'("     mesh =  xyz",11x,
     &            " # mesh type is xyz scoring mesh")')

                  call echmty(0,iot,'x',itxty(m),itxnm(m),
     &                        rtxdl(m),rtxmi(m),rtxma(m),itxrg(m)
     &                        ,abs(itxnm(m)+1),das_itxrg)
                  call echmty(0,iot,'y',ityty(m),itynm(m),
     &                        rtydl(m),rtymi(m),rtyma(m),ityrg(m)
     &                        ,abs(itynm(m)+1),das_ityrg)
                  call echmty(0,iot,'z',itzty(m),itznm(m),
     &                        rtzdl(m),rtzmi(m),rtzma(m),itzrg(m)
     &                        ,abs(itznm(m)+1),das_itzrg)

            else if( itmsh(m) .eq. 4 ) then

                  write(iot,'("     mesh =  tet",11x,
     &            " # mesh type is tetra geometry")')

                  call echtet(iot,itrgm(m),idas_itreg(itreg(m)))

            end if

*-----------------------------------------------------------------------

            if( itmtr(m,1) .gt. 0 ) then

               if( itmtr(m,1) .eq. 1 ) then

                     write(iot,'("     trcl = ",i6,9x,
     &               " # transform ID")') itmtr(m,3)

               else

                  if( itmtr(m,2) .eq. 0 ) then

                     write(iot,'("     trcl = ",
     &                           3(1p1e15.7),
     &                           3(/12x,3(1p1e15.7)),i5)')
     &               ( rtmtr(m,j), j = 1, 12 ),
     &               nint( rtmtr(m,13) )

                  else

                     write(iot,'("    *trcl = ",
     &                           3(1p1e15.7),
     &                           3(/12x,3(1p1e15.7)),i5)')
     &               ( rtmtr(m,j), j = 1, 12 ),
     &               nint( rtmtr(m,13) )

                  end if

               end if

            end if

*-----------------------------------------------------------------------

                  call echmty(0,iot,'e',itety(m),itenm(m),
     &                        rtedl(m),rtemi(m),rtema(m),iterg(m)
     &                        ,abs(itenm(m)+1),das_iterg)

*-----------------------------------------------------------------------

               if( ittty(m) .ne. 0 ) then

                  call echmty(0,iot,'t',ittty(m),ittnm(m),
     &                        rttdl(m),rttmi(m),rttma(m),ittrg(m)
     &                        ,abs(ittnm(m)+1),das_ittrg)

               end if

*-----------------------------------------------------------------------

               if( itaty(m) .ne. 0 ) then

                  call echmty(0,iot,'a',itaty(m),itanm(m),
     &                        rtadl(m),rtami(m),rtama(m),itarg(m)
     &                        ,abs(itanm(m)+1),das_itarg)

               end if

*-----------------------------------------------------------------------
               write(iot,'("  eng2let =",i5,11x,
     &         " # (D=0) Flag to convert energy to LET, 0: not convert",
     &                                               ", 1: convert")')
     &         ite2l(m)

               if ( ite2l(m) .ne. 0 ) then
                  write(iot,'("   letmat =",i5,11x,
     &            " # (D=0) mat ID for LET, 0:real mat, ",
     &                                    "<0: electron for H2O")')
     &            itlmt(m)
               end if

*-----------------------------------------------------------------------

               if( itmtn(m) .eq. 0 ) then

                  write(iot,'(" material =  all"11x,
     &            " # (D=all) number of specific material")')

               else if( itmtn(m) .gt. 0 ) then

                  write(iot,'(" material =",i5,11x,
     &            " # (D=all) number of specific material")')
     &            itmtn(m) * itmcn(m)

                     k = 0

                  do j = 1, itmtn(m)

                     k = k + 1

                     jmat(k) = ismte( itmtt(m) + j - 1 )

                     if( k / 10 * 10 .eq. k .or. j .eq. itmtn(m) ) then

                        write(iot,'(11x,10(i6))') ( jmat(i),i = 1, k )

                        k = 0

                     end if

                  end do

               end if

*-----------------------------------------------------------------------

               if( itman(m) .eq. 0 ) then

                  write(iot,'("   mother =  all"11x,
     &            " # (D=all) number of specific mother")')

               else if( itman(m) .gt. 0 ) then

                  write(iot,'("   mother = ",i4,11x,
     &            " # (D=all) number of specific mother")')
     &            itman(m) * itmct(m)

                     k = 0

                  do j = 1, itman(m)

                     k = k + 1

                     iaz = ismat( itmat(m) + j - 1 )

                     iz = iaz / 1000
                     ia = iaz - iz * 1000

                     if( ia .eq. 0 ) then

                        chaus(k) = elmnt(iz)

                     else

                        call chname(idum,ia,iz,chau)

                        chaus(k) = chau

                     end if

                     if( k / 7 * 7 .eq. k .or. j .eq. itman(m) ) then

                        write(iot,'(12x,7(a8))') ( chaus(i),i = 1, k )

                        k = 0

                     end if

                  end do

               end if

*-----------------------------------------------------------------------

               if( itunt(m) .eq. 1 ) then

                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [1/source]")') itunt(m)

               else if( itunt(m) .eq. 2 ) then

                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [1/cm^3/source]")') itunt(m)

               else if( itunt(m) .eq. 3 ) then

                if( ite2l(m).eq.0 ) then  ! ccse 2022/09/30
                if ( iMeVperu.eq. 1 ) then
                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [1/(MeV/n)/source]")') itunt(m)
                else
                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [1/MeV/source]")') itunt(m)
                end if
                else if( ite2l(m).eq.1 ) then
                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [1/(keV/um)/source]")') itunt(m)
                end if

               else if( itunt(m) .eq. 4 ) then

                if( ite2l(m).eq.0 ) then  ! ccse 2022/09/30
                if ( iMeVperu.eq. 1 ) then
                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [1/cm^3/(MeV/n)/source]")') itunt(m)
                else
                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [1/cm^3/MeV/source]")') itunt(m)
                end if
                else if( ite2l(m).eq.1 ) then
                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [1/cm^3/(keV/um)/source]")') itunt(m)
                end if

               else if( itunt(m) .eq. 5 ) then

                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [1/Lethargy/source]")') itunt(m)

               else if( itunt(m) .eq. 6 ) then

                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [1/cm^3/Lethargy/source]")') itunt(m)

               else if( itunt(m) .eq. 11 ) then

                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [1/nsec/source]")') itunt(m)

               else if( itunt(m) .eq. 12 ) then

                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [1/cm^3/nsec/source]")') itunt(m)

               else if( itunt(m) .eq. 13 ) then

                if( ite2l(m).eq.0 ) then  ! ccse 2022/09/30
                if ( iMeVperu.eq. 1 ) then
                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [1/(MeV/n)/nsec/source]")') itunt(m)
                else
                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [1/MeV/nsec/source]")') itunt(m)
                end if
                else if( ite2l(m).eq.1 ) then
                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [1/(keV/um)/nsec/source]")') itunt(m)
                end if

               else if( itunt(m) .eq. 14 ) then

                if( ite2l(m).eq.0 ) then  ! ccse 2022/09/30
                if ( iMeVperu.eq. 1 ) then
                  write(iot,'("     unit = ",i4,11x,

     &           " # unit is [1/cm^3/(MeV/n)/nsec/source]")') itunt(m)
                else
                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [1/cm^3/MeV/nsec/source]")') itunt(m)
                end if
                else if( ite2l(m).eq.1 ) then
                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [1/cm^3/(keV/um)/nsec/source]")') itunt(m)
                end if

               else if( itunt(m) .eq. 15 ) then

                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [1/Lethargy/nsec/source]")') itunt(m)

               else if( itunt(m) .eq. 16 ) then

                  write(iot,'("     unit = ",i4,11x,
     &          " # unit is [1/cm^3/Lethargy/nsec/source]")') itunt(m)

               else if( itunt(m) .eq. 21 ) then

                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [1/sr/source]")') itunt(m)

               else if( itunt(m) .eq. 22 ) then

                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [1/cm^3/sr/source]")') itunt(m)

               else if( itunt(m) .eq. 23 ) then

                if( ite2l(m).eq.0 ) then  ! ccse 2022/09/30
                if ( iMeVperu.eq. 1 ) then
                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [1/(MeV/n)/sr/source]")') itunt(m)
                else
                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [1/MeV/sr/source]")') itunt(m)
                end if
                else if( ite2l(m).eq.1 ) then
                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [1/(keV/um)/sr/source]")') itunt(m)
                end if

               else if( itunt(m) .eq. 24 ) then

                if( ite2l(m).eq.0 ) then  ! ccse 2022/09/30
                if ( iMeVperu.eq. 1 ) then
                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [1/cm^3/(MeV/n)/sr/source]")') itunt(m)
                else
                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [1/cm^3/MeV/sr/source]")') itunt(m)
                end if
                else if( ite2l(m).eq.1 ) then
                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [1/cm^3/(keV/um)/sr/source]")') itunt(m)
                end if

               else if( itunt(m) .eq. 25 ) then

                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [1/Lethargy/sr/source]")') itunt(m)

               else if( itunt(m) .eq. 26 ) then

                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [1/cm^3/Lethargy/sr/source]")') itunt(m)

               else if( itunt(m) .eq. 31 ) then

                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [1/nsec/sr/source]")') itunt(m)

               else if( itunt(m) .eq. 32 ) then

                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [1/cm^3/nsec/sr/source]")') itunt(m)

               else if( itunt(m) .eq. 33 ) then

                if( ite2l(m).eq.0 ) then  ! ccse 2022/09/30
                if ( iMeVperu.eq. 1 ) then
                  write(iot,'("     unit = ",i4,11x,

     &            " # unit is [1/(MeV/n)/nsec/sr/source]")') itunt(m)
                else
                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [1/MeV/nsec/sr/source]")') itunt(m)
                end if
                else if( ite2l(m).eq.1 ) then
                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [1/(keV/um)/nsec/sr/source]")') itunt(m)
                end if

               else if( itunt(m) .eq. 34 ) then

                if( ite2l(m).eq.0 ) then  ! ccse 2022/09/30
                if ( iMeVperu.eq. 1 ) then
                  write(iot,'("     unit = ",i4,11x,

     &        " # unit is [1/cm^3/(MeV/n)/nsec/sr/source]")') itunt(m)
                else
                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [1/cm^3/MeV/nsec/sr/source]")') itunt(m)
                end if
                else if( ite2l(m).eq.1 ) then
                  write(iot,'("     unit = ",i4,11x,
     &       " # unit is [1/cm^3/(keV/um)/nsec/sr/source]")') itunt(m)
                end if

               else if( itunt(m) .eq. 35 ) then

                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [1/Lethargy/nsec/sr/source]")') itunt(m)

               else if( itunt(m) .eq. 36 ) then

                  write(iot,'("     unit = ",i4,11x,
     &          " # unit is [1/cm^3/Lethargy/nsec/sr/source]")')
     &             itunt(m)

               end if

*-----------------------------------------------------------------------

                  asfil = '  # file name of output for the above axis'
                  bsfil = '  # file name of dumped data'
                  csfil = '  # file name of dumped data summary'

            if( icax .eq. 0 ) then

                  icn = itaxn(m)

               do i = 1, icn

                  ian(i) = i

               end  do

            else

                  icn = 1
                  ian(1) = iaxs

            end if

*-----------------------------------------------------------------------

               idtyp = 0

            do i = 1, icn

               if( itaxs(m,ian(i)) .ge. 7 .and.
     &             itaxs(m,ian(i)) .le. 10 )  idtyp = ittwo(m)

            end do

            if( idtyp .ne. 0 ) then

               write(iot,'("  2D-type = ",i4,11x,
     &         " # 1:Cont, 2:Clust, 3:Color, 4:xyz, 5:mat,",
     &         " 6:Clust+Cont, 7:Col+Cont")') idtyp

            end if

*-----------------------------------------------------------------------

            do i = 1, icn

                  iax = ian(i)

                  write(iot,'("     axis =  ",a3,11x,
     &            " # axis of output")')
     &            aname( itaxs(m,iax) )

             if ( itfll(m,iax) .gt. 0 ) then

               if( itmdp(m,0) .eq. 0 ) then

                  msfile = max( 14, itfll(m,iax) )

                  write(iot,'("     file = ",100a1)')
     &                  ( ctfln(m,iax)(j:j), j = 1, msfile ),
     &                  ( asfil(j:j), j = 1, 42 )

               else

C S.H. added for Dump-Restart on 2014/5/7
                  msfile = max( 14, itfll(m,iax) )

                  write(iot,'("     file = ",100a1)')
     &                  ( ctfln(m,iax)(j:j), j = 1, msfile ),
     &                  ( csfil(j:j), j = 1, 36 )

                  msfile = max( 14, itfll(m,iax)+4 )

                  do k = 1, msfile
                     dsfil(k:k) = ' '
                  end do

                  do is = itfll(m,iax), 1, -1
                     if ( ctfln(m,iax)(is:is) .eq. '.' ) exit
                  end do
                  if ( is .eq. 0 ) is = itfll(m,iax)+1

                  dsfil(1:is-1) = ctfln(m,iax)(1:is-1)
                  dsfil(is:is+3) = '_dmp'
                  if ( is .lt. itfll(m,iax) ) then
                     dsfil(is+4:itfll(m,iax)+4)
     &                    = ctfln(m,iax)(is:itfll(m,iax))
                  end if

                  write(iot,'("#    file = ",100a1)')
     &                  ( dsfil(j:j), j = 1, msfile ),
     &                  ( bsfil(j:j), j = 1, 28 )

               end if

             end if

            end do

*-----------------------------------------------------------------------

               if( itout(m) .eq. 1 ) then

                  write(iot,'("   output = source",8x,
     &            "  # (D=nuclear) source distribution")')

               else if( itout(m) .eq. 2 ) then

                  write(iot,'("   output = nuclear",7x,
     &            "  # (D=nuclear) products from all reactions")')

               else if( itout(m) .eq. 3 ) then

                  write(iot,'("   output = decay",9x,
     &           "  # (D=nuclear) products from decay")')

               else if( itout(m) .eq. 4 ) then

                  write(iot,'("   output = fission",7x,
     &           "  # (D=nuclear) products from fission")')

               else if( itout(m) .eq. 5 ) then

                  write(iot,'("   output = elastic",7x,
     &           "  # (D=nuclear) products from elastic")')

               else if( itout(m) .eq. 6 ) then

                  write(iot,'("   output = nonela",8x,
     &           "  # (D=nuclear) products from nonelastic")')

               else if( itout(m) .eq. 7 ) then

                  write(iot,'("   output = atomic",8x,
     &           "  # (D=nuclear) atomic reaction")')

               end if

*-----------------------------------------------------------------------

               if( iprim(m) .eq. 1 ) then

                  write(iot,'("  primary = 1",13x,
     &            "  # (D=1) score primary particles")')

               else if( iprim(m) .eq. 0 ) then

                  write(iot,'("    primary = 0",7x,
     &            "  # (D=1) do not score primary particles")')

               end if

*-----------------------------------------------------------------------
               call echprt(1,iot,m,itpan,itpat,jtpat,itnm,itmxpt,itmxgp)

*-----------------------------------------------------------------------

               if( rtfac(m) .ne. 1.0d0 ) then

                  write(iot,'("   factor =",1p1g14.7,1x,1x,
     &            " # (D=1.0) normalization factor")')
     &            rtfac(m)

               end if

*-----------------------------------------------------------------------

               if( itmtn(m) .gt. 0 .and. itmsh(m) .eq. 3 ) then

                  write(iot,'("   volmat = ",i4,11x,
     &            " # (D=9) vol.corr., 0:no, number of scan / dir.")')
     &            itvm(m)

               end if

*-----------------------------------------------------------------------

               if( itanl(m) .gt. 0 ) then

                  write(iot,'("    angel = ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               end if

*-----------------------------------------------------------------------

               if( itsans(m) .gt. 0 ) then

                  write(iot,'("   sangel = ",i2)') itsans(m)
                  call write_sangel(iot,m,1)

               end if

*-----------------------------------------------------------------------

               if( itaxl(m) .gt. 0 ) then

                  write(iot,'("    x-txt = ",200a1)')
     &            ( itaxt(m)(i:i),i = 1, itaxl(m) )

               end if

*-----------------------------------------------------------------------

               if( itayl(m) .gt. 0 ) then

                  write(iot,'("    y-txt = ",200a1)')
     &            ( itayt(m)(i:i),i = 1, itayl(m) )

               end if

*-----------------------------------------------------------------------

               if( itazl(m) .gt. 0 ) then

                  write(iot,'("    z-txt = ",200a1)')
     &            ( itazt(m)(i:i),i = 1, itazl(m) )

               end if

*-----------------------------------------------------------------------

               if( itgsh(m) .gt. 0 ) then

                  write(iot,'("    gshow = ",i4,10x,
     &            "  # 0: no 1:bnd, 2:bnd+mat, ",
     &     "3:bnd+reg 4:bnd+lat 5:bmp style")')
     &            itgsh(m)

               end if

*-----------------------------------------------------------------------

               if( itglt(m) .ge. 0 ) then

                  write(iot,'("    gslat = ",i4,10x,
     &            "  # 0: no line for lat, 1: with,",
     &               " 2: no line between the same cell")')
     &            itglt(m)

               end if

*-----------------------------------------------------------------------

               if( itgsh(m) .gt. 0 .and. itger(m) .gt. 0 ) then

                  write(iot,'("    ginfo = ",i4,11x,
     &            " # (D=0) 0:no, 1:show in graph, 2:add file")')
     &            itger(m)

               end if

*-----------------------------------------------------------------------

               if( itrsh(m) .gt. 0 ) then

                  write(iot,'("    rshow = ",i4,10x,
     &            "  # 0: no 1:bnd, 2:bnd+mat, ",
     &     "3:bnd+reg 4:bnd+lat 5:bmp style")')
     &            itrsh(m)

                  call echmty(1,iot,'x',itxty(m),itxnm(m),
     &                        rtxdl(m),rtxmi(m),rtxma(m),itxrg(m)
     &                        ,abs(itxnm(m)+1),das_itxrg)
                  call echmty(1,iot,'y',ityty(m),itynm(m),
     &                        rtydl(m),rtymi(m),rtyma(m),ityrg(m)
     &                        ,abs(itynm(m)+1),das_ityrg)
                  call echmty(1,iot,'z',itzty(m),itznm(m),
     &                        rtzdl(m),rtzmi(m),rtzma(m),itzrg(m)
     &                        ,abs(itznm(m)+1),das_itzrg)

               end if

*-----------------------------------------------------------------------

               if( itgsh(m) .gt. 0 .or. itrsh(m) .gt. 0 ) then

                  write(iot,'("    resol = ",i4,10x,
     &            "  # (D=1) resolution of gshow or rshow")')
     &            itres(m)

                  write(iot,'("    width =",1p1g14.7,1x,1x,
     &            " # (D=0.5) width of lines for gshow or rshow")')
     &            rtwid(m)

               end if

*-----------------------------------------------------------------------

               if( iterl(m) .ne. 72 ) then

                  write(iot,'("   iechrl =",i5,10x,
     &            "  # (D=72) length reg of val or vol")')
     &            iterl(m)

               end if

*-----------------------------------------------------------------------

               if( iteps(m) .ne. 0 ) then

                  write(iot,'("   epsout =",i5,10x,
     &            "  # (D=0) generate eps file by ANGEL")')
     &            iteps(m)

               end if

*-----------------------------------------------------------------------

               if( itbmp(m) .ne. 0 ) then

                  write(iot,'("   bmpout =",i5,10x,
     &            "  # (D=0) generate bitmap file")')
     &            itbmp(m)

               end if

*-----------------------------------------------------------------------

               if( itvtk(m) .ne. 0 ) then

                  write(iot,'("   vtkout =",i5,10x,
     &            "  # (D=0) generate vtk file")')
     &            itvtk(m)

                  if ( itvtkfmt(m) .eq. 0 ) then

                     write(iot,'("   vtkfmt =",i5,10x,
     &               "  # (D=0) generate vtk ascii file")')
     &               itvtkfmt(m)

                  else

                     write(iot,'("   vtkfmt =",i5,10x,
     &               "  # (D=0) generate vtk binary file")')
     &               itvtkfmt(m)

                  end if

               end if

*-----------------------------------------------------------------------
*        check of counter
*-----------------------------------------------------------------------

            do i = 1, 3

               if( itcnt(i,m) .ne. 0 ) then

                  write(iot,'(" ctmin(",i1,") =",i5,10x,
     &            "  # (D=-9999) min. of counter ",i1)')
     &            i, itcnt(i*2+2,m), i

                  write(iot,'(" ctmax(",i1,") =",i5,10x,
     &            "  # (D= 9999) max. of counter ",i1)')
     &            i, itcnt(i*2+3,m), i

                  if( ncntc(i) .eq. 0 )
     &            write(iot,'("# Warning: there is no counter = ",
     &                        i1," definition in [counter].")') i

               end if

            end do

            do i = 1, 3
               if( itcnth(i,m) .ne. 0 ) then
                  write(iot,'(" chmin(",i1,") =",i5,10x,
     &            "  # (D=-9999) min. of history counter ",i1)')
     &            i, itcnth(i*2+2,m), i
                  write(iot,'(" chmax(",i1,") =",i5,10x,
     &            "  # (D= 9999) max. of history counter ",i1)')
     &            i, itcnth(i*2+3,m), i
               end if
            end do

*-----------------------------------------------------------------------

               if( rtstd(m) .gt. 0.d0 ) then

                  write(iot,'("   stdcut =",1p1g14.7,1x,1x,
     &            " # (D=-1) threshold value of STD cutoff")')
     &            rtstd(m)

               end if

*-----------------------------------------------------------------------
               if( itfoam(m) .ne. 0 ) then

                  write(iot,'("  foamout =",i5,10x,
     &            "  # (D=0) generate OpenFOAM data")')
     &            itfoam(m)

               end if

*-----------------------------------------------------------------------
               if( itmxang(m) .ne. 0 ) then

                  write(iot,'(" maxangel =",i5,10x,
     &            "  # (D=0) number of output particles")')
     &            itmxang(m)

               end if

*-----------------------------------------------------------------------
*        dump the data
*-----------------------------------------------------------------------

         if( itmdp(m,0) .ne. 0 ) then

                  write(iot,'("     dump =",i5,11x,
     &            " # (D=0) number of dumped data,",
     &            " <0: ascii, >0: binary")') itmdp(m,0)

                  write(iot,'(13x,30(i4))')
     &               ( itmdp(m,j), j = 1, abs( itmdp(m,0) ) )

                  write(iot,'("# dump data  ",30(a4))')
     &               ( dmpc(itmdp(m,j)), j = 1, abs( itmdp(m,0) ) )

         end if

*-----------------------------------------------------------------------
*        sumtally subsection
*-----------------------------------------------------------------------
            if ( icntl==13 .and. nsumtalRead(m)==8 ) then

               call sumtal_echo(iot,m)

            end if
*-----------------------------------------------------------------------

            write(iot,'("#    used :",
     &                  "        main (  %)",
     &                  "        temp (  %)",
     &                  "       total (  %)"/
     &                  "#  memory :",
     &                  3(i12," (",i3,")"))')
     &                  itsmn(m), nint(dble(itsmn(m))/dble(mdas)*100),
     &                  itstm(m), nint(dble(itstm(m))/dble(mdas)*100),
     &                  itsmn(m)+itstm(m),
     &                  nint(dble(itsmn(m)+itstm(m))/dble(mdas)*100)

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine ttimech(iot,m,iaxs,icax)
*                                                                      *
*       input echo of the time tally                                   *
*       last modified by K.Niita on 2005/08/15                         *
*                                                                      *
************************************************************************

      use sumtallymod
      use sangelmod, only: itsans
      use partmod ! frtati 2021/10/05
      use moddas_mesh
      use moddas_region
      use moddas_tally
*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'

*-----------------------------------------------------------------------

      common /tall00/ itnm, ital(itlmax), itals(itlmax), italm(itlmax)
      common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                itout(itlmax), ittwo(itlmax)

      common /tall02/ itrgn(itlmax), itrgm(itlmax), itreg(itlmax)
      common /tall24/ itrcm(itlmax), itrcg(itlmax)
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


      common /tall20/ ittty(itlmax), ittnm(itlmax), ittrg(itlmax),
     &                rttmi(itlmax), rttma(itlmax), rttdl(itlmax)

      common /tall08/ rtrx0(itlmax), rtry0(itlmax)

      common /tall10/ itaxn(itlmax), itaxs(itlmax,6)
      common /tall11/ itfln(itlmax), itfll(itlmax,6), ctfln(itlmax,6)
      character ctfln*100

      common /tall12/ itman(itlmax), itmat(itlmax), itmct(itlmax),
     &                itnun(itlmax), itnuc(itlmax),
     &                itndz(itlmax), itndn(itlmax),
     &                itnkz(itlmax), itnkn(itlmax)

      common /tall16/ itnfr(itlmax), itnfz(itlmax), itnfn(itlmax)

      common /tall14/ ittll(itlmax), ittle(itlmax)
      character ittle*80

      common /tall15/ itanl(itlmax), itang(itlmax)
      character itang*200

      common /tall27/ itaxl(itlmax), itaxt(itlmax),
     &                itayl(itlmax), itayt(itlmax),
     &                itazl(itlmax), itazt(itlmax)
      character itaxt*200, itayt*200, itazt*200

      common /tall19/ itsmn(itlmax), itstm(itlmax)

      common /tall21/ rtfac(itlmax)
      common /tall28/ itgsh(itlmax)
      common /tall29/ itrsh(itlmax)
      common /tall31/ iterl(itlmax)
      common /tall22/ itmcn(itlmax), itmtn(itlmax), itmtt(itlmax)
      common /tall34/ itvm(itlmax)
      common /tall35/ iteps(itlmax)
      common /tall37/ itcnt(9,itlmax)

      common /cntreg/ ncntc(3), ncreg(3), incrc(3), incrt(3), kcont(3)
      common /cntpat/ icpan(3), icpat(3,20,2)

      common /tall38/ itres(itlmax)
      common /tall39/ rtwid(itlmax)
      common /tall46/ itmtr(itlmax,4), rtmtr(itlmax,13)
      common /tall48/ itmdp(itlmax,0:30), itmdf(itlmax)
      common /tall49/ itglt(itlmax)

      common /tall60/ itger(itlmax)
      common /tall63/ itbmp(itlmax)
      common /tall64/ itvtk(itlmax),itvtkfmt(itlmax)

      common /tall67/ itstd(itlmax), rtstd(itlmax)

      common /tall82/ itcnth(9,itlmax)

      dimension idas(mdas*2) !20220906frtati temporary for fbounds-check
      equivalence ( das, idas )

*-----------------------------------------------------------------------

      common /inpec/ ititl, ipara, ibody, iregn, llarr, itby, itar

*-----------------------------------------------------------------------
       integer   icntl
       integer   inucr

       common /tcntl/  icntl, inucr
!-----------------------------------------------------------------------

      character asfil*100
      character bsfil*100
      character csfil*100
      character dsfil*100

      character aname(11)*3
      data aname / 'eng','reg','  x','  y','  z','  r',
     &             ' xy',' yz',' xz',' rz','  t'/

      dimension ian(6)
      dimension jmat(20)

      character dmpc(30)*4
      data dmpc / '  kf','   x','   y','   z','   u','   v','   w',
     &            '   e','  wt','  tm','  c1','  c2','  c3',
     &            '  sx','  sy','  sz','  n0','  nc','  nb','  no',
     &            '    ','    ','    ','    ','    ','    ',
     &            '    ','    ','    ','    '/

      integer           iMeVperu
      common /cMeVperu/ iMeVperu

*-----------------------------------------------------------------------
      common /tall79/ itallech
      if( itallech .eq. 0 ) return

*-----------------------------------------------------------------------
*     input echo
*-----------------------------------------------------------------------

            write(iot,'("[ T-Time ]")')

                  write(iot,'("    title = ",80a1)')
     &            ( ittle(m)(i:i), i = 1, ittll(m) )

*-----------------------------------------------------------------------

            if( itmsh(m) .eq. 1 ) then

                  write(iot,'("     mesh =  reg",11x,
     &            " # mesh type is region-wise")')

                     idas1 = mmmax
                     idas2 = ( idas1 + itrgn(m) - 1 ) * 2 + 1
                     idas3 = idas2 + itrgn(m)

                  call echrg(1,iterl(m),0,1,iot,itrcm(m),
     &                       idas_itrcg(itrcg(m)),
     &                       itrgn(m),itrgm(m),idas_itreg(itreg(m)),
     &                       das(idas1),idas(idas2),
     &                       itrnv(m),idas(itriv(m)),das(itrrv(m)),
     &                       idas3)

            else if( itmsh(m) .eq. 2 ) then

                  write(iot,'("     mesh =  r-z",11x,
     &            " # mesh type is r-z scoring mesh")')

                  write(iot,'("       x0 = ",1p1g14.7,1x,
     &            " # center x-position of r-z mesh")') rtrx0(m)

                  write(iot,'("       y0 = ",1p1g14.7,1x,
     &            " # center y-position of r-z mesh")') rtry0(m)

                  call echmty(0,iot,'r',itrty(m),itrnm(m),
     &                        rtrdl(m),rtrmi(m),rtrma(m),itrrg(m)
     &                        ,abs(itrnm(m)+1),das_itrrg)
                  call echmty(0,iot,'z',itzty(m),itznm(m),
     &                        rtzdl(m),rtzmi(m),rtzma(m),itzrg(m)
     &                        ,abs(itznm(m)+1),das_itzrg)

            else if( itmsh(m) .eq. 3 ) then

                  write(iot,'("     mesh =  xyz",11x,
     &            " # mesh type is xyz scoring mesh")')

                  call echmty(0,iot,'x',itxty(m),itxnm(m),
     &                        rtxdl(m),rtxmi(m),rtxma(m),itxrg(m)
     &                        ,abs(itxnm(m)+1),das_itxrg)
                  call echmty(0,iot,'y',ityty(m),itynm(m),
     &                        rtydl(m),rtymi(m),rtyma(m),ityrg(m)
     &                        ,abs(itynm(m)+1),das_ityrg)
                  call echmty(0,iot,'z',itzty(m),itznm(m),
     &                        rtzdl(m),rtzmi(m),rtzma(m),itzrg(m)
     &                        ,abs(itznm(m)+1),das_itzrg)

            end if

*-----------------------------------------------------------------------

            if( itmtr(m,1) .gt. 0 ) then

               if( itmtr(m,1) .eq. 1 ) then

                     write(iot,'("     trcl = ",i6,9x,
     &               " # transform ID")') itmtr(m,3)

               else

                  if( itmtr(m,2) .eq. 0 ) then

                     write(iot,'("     trcl = ",
     &                           3(1p1e15.7),
     &                           3(/12x,3(1p1e15.7)),i5)')
     &               ( rtmtr(m,j), j = 1, 12 ),
     &               nint( rtmtr(m,13) )

                  else

                     write(iot,'("    *trcl = ",
     &                           3(1p1e15.7),
     &                           3(/12x,3(1p1e15.7)),i5)')
     &               ( rtmtr(m,j), j = 1, 12 ),
     &               nint( rtmtr(m,13) )

                  end if

               end if

            end if

*-----------------------------------------------------------------------

                  call echmty(0,iot,'e',itety(m),itenm(m),
     &                        rtedl(m),rtemi(m),rtema(m),iterg(m)
     &                        ,abs(itenm(m)+1),das_iterg)

*-----------------------------------------------------------------------

                  call echmty(0,iot,'t',ittty(m),ittnm(m),
     &                        rttdl(m),rttmi(m),rttma(m),ittrg(m)
     &                        ,abs(ittnm(m)+1),das_ittrg)

*-----------------------------------------------------------------------

               if( itunt(m) .eq. 1 ) then

                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [1/source]")') itunt(m)

               else if( itunt(m) .eq. 2 ) then

                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [1/nsec/source]")') itunt(m)


               else if( itunt(m) .eq. 3 ) then

                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [1/nsec/cm^3/source]")') itunt(m)


               else if( itunt(m) .eq. 4 ) then

                if ( iMeVperu.eq. 1 ) then
                  write(iot,'("     unit = ",i4,11x,
     &           " # unit is [1/nsec/cm^3/(MeV/n)/source]")') itunt(m)
                else
                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [1/nsec/cm^3/MeV/source]")') itunt(m)
                end if

               end if

*-----------------------------------------------------------------------

                  asfil = '  # file name of output for the above axis'
                  bsfil = '  # file name of dumped data'
                  csfil = '  # file name of dumped data summary'

            if( icax .eq. 0 ) then

                  icn = itaxn(m)

               do i = 1, icn

                  ian(i) = i

               end  do

            else

                  icn = 1
                  ian(1) = iaxs

            end if

*-----------------------------------------------------------------------

               if( itmtn(m) .eq. 0 ) then

                  write(iot,'(" material =  all"11x,
     &            " # (D=all) number of specific material")')

               else if( itmtn(m) .gt. 0 ) then

                  write(iot,'(" material =",i5,11x,
     &            " # (D=all) number of specific material")')
     &            itmtn(m) * itmcn(m)

                     k = 0

                  do j = 1, itmtn(m)

                     k = k + 1

                     jmat(k) = ismte( itmtt(m) + j - 1 )

                     if( k / 10 * 10 .eq. k .or. j .eq. itmtn(m) ) then

                        write(iot,'(11x,10(i6))') ( jmat(i),i = 1, k )

                        k = 0

                     end if

                  end do

               end if

*-----------------------------------------------------------------------

               idtyp = 0

            do i = 1, icn

               if( itaxs(m,ian(i)) .ge. 7 .and.
     &             itaxs(m,ian(i)) .le. 10 )  idtyp = ittwo(m)

            end do

            if( idtyp .ne. 0 ) then

               write(iot,'("  2D-type = ",i4,11x,
     &         " # 1:Cont, 2:Clust, 3:Color, 4:xyz, 5:mat,",
     &         " 6:Clust+Cont, 7:Col+Cont")') idtyp

            end if

*-----------------------------------------------------------------------

            do i = 1, icn

                  iax = ian(i)

                  write(iot,'("     axis =  ",a3,11x,
     &            " # axis of output")')
     &            aname( itaxs(m,iax) )

             if ( itfll(m,iax) .gt. 0 ) then

               if( itmdp(m,0) .eq. 0 ) then

                  msfile = max( 14, itfll(m,iax) )

                  write(iot,'("     file = ",100a1)')
     &                  ( ctfln(m,iax)(j:j), j = 1, msfile ),
     &                  ( asfil(j:j), j = 1, 42 )

               else

C S.H. added for Dump-Restart on 2014/5/7
                  msfile = max( 14, itfll(m,iax) )

                  write(iot,'("     file = ",100a1)')
     &                  ( ctfln(m,iax)(j:j), j = 1, msfile ),
     &                  ( csfil(j:j), j = 1, 36 )

                  msfile = max( 14, itfll(m,iax)+4 )

                  do k = 1, msfile
                     dsfil(k:k) = ' '
                  end do

                  do is = itfll(m,iax), 1, -1
                     if ( ctfln(m,iax)(is:is) .eq. '.' ) exit
                  end do
                  if ( is .eq. 0 ) is = itfll(m,iax)+1

                  dsfil(1:is-1) = ctfln(m,iax)(1:is-1)
                  dsfil(is:is+3) = '_dmp'
                  if ( is .lt. itfll(m,iax) ) then
                     dsfil(is+4:itfll(m,iax)+4)
     &                    = ctfln(m,iax)(is:itfll(m,iax))
                  end if

                  write(iot,'("#    file = ",100a1)')
     &                  ( dsfil(j:j), j = 1, msfile ),
     &                  ( bsfil(j:j), j = 1, 28 )

               end if

             end if

            end do

*-----------------------------------------------------------------------

               if( itout(m) .eq. 1 ) then

                  write(iot,'("   output =  all",11x,
     &           "  # (D=all) time for cutoff and escape particles")')

               else if( itout(m) .eq. 2 ) then

                  write(iot,'("   output = cutoff",8x,
     &            "  # (D=all) time for cutoff particles")')

               else if( itout(m) .eq. 3 ) then

                  write(iot,'("   output = escape",8x,
     &            "  # (D=all) time for escape particles")')

               else if( itout(m) .eq. 4 ) then

                  write(iot,'("   output = decay",9x,
     &            "  # (D=all) time for decay particles")')

               end if

*-----------------------------------------------------------------------
               call echprt(1,iot,m,itpan,itpat,jtpat,itnm,itmxpt,itmxgp)

*-----------------------------------------------------------------------

               if( rtfac(m) .ne. 1.0d0 ) then

                  write(iot,'("   factor =",1p1g14.7,1x,1x,
     &            " # (D=1.0) normalization factor")')
     &            rtfac(m)

               end if

*-----------------------------------------------------------------------

               if( itmtn(m) .gt. 0 .and. itmsh(m) .eq. 3 ) then

                  write(iot,'("   volmat = ",i4,11x,
     &            " # (D=9) vol.corr., 0:no, number of scan / dir.")')
     &            itvm(m)

               end if

*-----------------------------------------------------------------------

               if( itanl(m) .gt. 0 ) then

                  write(iot,'("    angel = ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               end if

*-----------------------------------------------------------------------

               if( itsans(m) .gt. 0 ) then

                  write(iot,'("   sangel = ",i2)') itsans(m)
                  call write_sangel(iot,m,1)

               end if

*-----------------------------------------------------------------------

               if( itaxl(m) .gt. 0 ) then

                  write(iot,'("    x-txt = ",200a1)')
     &            ( itaxt(m)(i:i),i = 1, itaxl(m) )

               end if

*-----------------------------------------------------------------------

               if( itayl(m) .gt. 0 ) then

                  write(iot,'("    y-txt = ",200a1)')
     &            ( itayt(m)(i:i),i = 1, itayl(m) )

               end if

*-----------------------------------------------------------------------

               if( itazl(m) .gt. 0 ) then

                  write(iot,'("    z-txt = ",200a1)')
     &            ( itazt(m)(i:i),i = 1, itazl(m) )

               end if

*-----------------------------------------------------------------------

               if( itgsh(m) .gt. 0 ) then

                  write(iot,'("    gshow = ",i4,10x,
     &            "  # 0: no 1:bnd, 2:bnd+mat, ",
     &     "3:bnd+reg 4:bnd+lat 5:bmp style")')
     &            itgsh(m)

               end if

*-----------------------------------------------------------------------

               if( itglt(m) .ge. 0 ) then

                  write(iot,'("    gslat = ",i4,10x,
     &            "  # 0: no line for lat, 1: with,",
     &               " 2: no line between the same cell")')
     &            itglt(m)

               end if

*-----------------------------------------------------------------------

               if( itgsh(m) .gt. 0 .and. itger(m) .gt. 0 ) then

                  write(iot,'("    ginfo = ",i4,11x,
     &            " # (D=0) 0:no, 1:show in graph, 2:add file")')
     &            itger(m)

               end if

*-----------------------------------------------------------------------

               if( itrsh(m) .gt. 0 ) then

                  write(iot,'("    rshow = ",i4,10x,
     &            "  # 0: no 1:bnd, 2:bnd+mat, ",
     &     "3:bnd+reg 4:bnd+lat 5:bmp style")')
     &            itrsh(m)

                  call echmty(1,iot,'x',itxty(m),itxnm(m),
     &                        rtxdl(m),rtxmi(m),rtxma(m),itxrg(m)
     &                        ,abs(itxnm(m)+1),das_itxrg)
                  call echmty(1,iot,'y',ityty(m),itynm(m),
     &                        rtydl(m),rtymi(m),rtyma(m),ityrg(m)
     &                        ,abs(itynm(m)+1),das_ityrg)
                  call echmty(1,iot,'z',itzty(m),itznm(m),
     &                        rtzdl(m),rtzmi(m),rtzma(m),itzrg(m)
     &                        ,abs(itznm(m)+1),das_itzrg)

               end if

*-----------------------------------------------------------------------

               if( itgsh(m) .gt. 0 .or. itrsh(m) .gt. 0 ) then

                  write(iot,'("    resol = ",i4,10x,
     &            "  # (D=1) resolution of gshow or rshow")')
     &            itres(m)

                  write(iot,'("    width =",1p1g14.7,1x,1x,
     &            " # (D=0.5) width of lines for gshow or rshow")')
     &            rtwid(m)

               end if

*-----------------------------------------------------------------------

               if( iterl(m) .ne. 72 ) then

                  write(iot,'("   iechrl =",i5,10x,
     &            "  # (D=72) length reg of val or vol")')
     &            iterl(m)

               end if

*-----------------------------------------------------------------------

               if( iteps(m) .ne. 0 ) then

                  write(iot,'("   epsout =",i5,10x,
     &            "  # (D=0) generate eps file by ANGEL")')
     &            iteps(m)

               end if

*-----------------------------------------------------------------------

               if( itbmp(m) .ne. 0 ) then

                  write(iot,'("   bmpout =",i5,10x,
     &            "  # (D=0) generate bitmap file")')
     &            itbmp(m)

               end if

*-----------------------------------------------------------------------

               if( itvtk(m) .ne. 0 ) then

                  write(iot,'("   vtkout =",i5,10x,
     &            "  # (D=0) generate vtk file")')
     &            itvtk(m)

                  if ( itvtkfmt(m) .eq. 0 ) then

                     write(iot,'("   vtkfmt =",i5,10x,
     &               "  # (D=0) generate vtk ascii file")')
     &               itvtkfmt(m)

                  else

                     write(iot,'("   vtkfmt =",i5,10x,
     &               "  # (D=0) generate vtk binary file")')
     &               itvtkfmt(m)

                  end if

               end if

*-----------------------------------------------------------------------
*        check of counter
*-----------------------------------------------------------------------

            do i = 1, 3

               if( itcnt(i,m) .ne. 0 ) then

                  write(iot,'(" ctmin(",i1,") =",i5,10x,
     &            "  # (D=-9999) min. of counter ",i1)')
     &            i, itcnt(i*2+2,m), i

                  write(iot,'(" ctmax(",i1,") =",i5,10x,
     &            "  # (D= 9999) max. of counter ",i1)')
     &            i, itcnt(i*2+3,m), i

                  if( ncntc(i) .eq. 0 )
     &            write(iot,'("# Warning: there is no counter = ",
     &                        i1," definition in [counter].")') i

               end if

            end do

            do i = 1, 3
               if( itcnth(i,m) .ne. 0 ) then
                  write(iot,'(" chmin(",i1,") =",i5,10x,
     &            "  # (D=-9999) min. of history counter ",i1)')
     &            i, itcnth(i*2+2,m), i
                  write(iot,'(" chmax(",i1,") =",i5,10x,
     &            "  # (D= 9999) max. of history counter ",i1)')
     &            i, itcnth(i*2+3,m), i
               end if
            end do

*-----------------------------------------------------------------------

               if( rtstd(m) .gt. 0.d0 ) then

                  write(iot,'("   stdcut =",1p1g14.7,1x,1x,
     &            " # (D=-1) threshold value of STD cutoff")')
     &            rtstd(m)

               end if

*-----------------------------------------------------------------------
               if( itmxang(m) .ne. 0 ) then

                  write(iot,'(" maxangel =",i5,10x,
     &            "  # (D=0) number of output particles")')
     &            itmxang(m)

               end if

*-----------------------------------------------------------------------
*        dump the data
*-----------------------------------------------------------------------

         if( itmdp(m,0) .ne. 0 ) then

                  write(iot,'("     dump =",i5,11x,
     &            " # (D=0) number of dumped data,",
     &            " <0: ascii, >0: binary")') itmdp(m,0)

                  write(iot,'(13x,30(i4))')
     &               ( itmdp(m,j), j = 1, abs( itmdp(m,0) ) )

                  write(iot,'("# dump data  ",30(a4))')
     &               ( dmpc(itmdp(m,j)), j = 1, abs( itmdp(m,0) ) )

         end if

*-----------------------------------------------------------------------
*        sumtally subsection
*-----------------------------------------------------------------------
            if ( icntl==13 .and. nsumtalRead(m)==6 ) then

               call sumtal_echo(iot,m)

            end if
*-----------------------------------------------------------------------

            write(iot,'("#    used :",
     &                  "        main (  %)",
     &                  "        temp (  %)",
     &                  "       total (  %)"/
     &                  "#  memory :",
     &                  3(i12," (",i3,")"))')
     &                  itsmn(m), nint(dble(itsmn(m))/dble(mdas)*100),
     &                  itstm(m), nint(dble(itstm(m))/dble(mdas)*100),
     &                  itsmn(m)+itstm(m),
     &                  nint(dble(itsmn(m)+itstm(m))/dble(mdas)*100)

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine tdpaech(iot,m,iaxs,icax)
*                                                                      *
*       input echo of the dpa tally                                    *
*       last modified by K.Niita on 2004/09/15                         *
*                                                                      *
************************************************************************

      use sumtallymod
      use sangelmod, only: itsans
      use partmod ! frtati 2021/10/05
      use moddas_mesh
      use moddas_region
      use moddas_tally
*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'

*-----------------------------------------------------------------------

      common /tall00/ itnm, ital(itlmax), itals(itlmax), italm(itlmax)
      common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                itout(itlmax), ittwo(itlmax)


      common /tall02/ itrgn(itlmax), itrgm(itlmax), itreg(itlmax)
      common /tall24/ itrcm(itlmax), itrcg(itlmax)
      common /tall25/ itrnv(itlmax), itriv(itlmax), itrrv(itlmax)

      common /tall03/ itrty(itlmax), itrnm(itlmax), itrrg(itlmax),
     &                rtrmi(itlmax), rtrma(itlmax), rtrdl(itlmax)
      common /tall04/ itxty(itlmax), itxnm(itlmax), itxrg(itlmax),
     &                rtxmi(itlmax), rtxma(itlmax), rtxdl(itlmax)
      common /tall05/ ityty(itlmax), itynm(itlmax), ityrg(itlmax),
     &                rtymi(itlmax), rtyma(itlmax), rtydl(itlmax)
      common /tall06/ itzty(itlmax), itznm(itlmax), itzrg(itlmax),
     &                rtzmi(itlmax), rtzma(itlmax), rtzdl(itlmax)

      common /tall08/ rtrx0(itlmax), rtry0(itlmax)

      common /tall10/ itaxn(itlmax), itaxs(itlmax,6)
      common /tall11/ itfln(itlmax), itfll(itlmax,6), ctfln(itlmax,6)
      character ctfln*100

      common /tall12/ itman(itlmax), itmat(itlmax), itmct(itlmax),
     &                itnun(itlmax), itnuc(itlmax),
     &                itndz(itlmax), itndn(itlmax),
     &                itnkz(itlmax), itnkn(itlmax)

      common /tall16/ itnfr(itlmax), itnfz(itlmax), itnfn(itlmax)

      common /tall14/ ittll(itlmax), ittle(itlmax)
      character ittle*80

      common /tall15/ itanl(itlmax), itang(itlmax)
      character itang*200

      common /tall27/ itaxl(itlmax), itaxt(itlmax),
     &                itayl(itlmax), itayt(itlmax),
     &                itazl(itlmax), itazt(itlmax)
      character itaxt*200, itayt*200, itazt*200

      common /tall19/ itsmn(itlmax), itstm(itlmax)

      common /tall21/ rtfac(itlmax)
      common /tall22/ itmcn(itlmax), itmtn(itlmax), itmtt(itlmax)
      common /tall28/ itgsh(itlmax)
      common /tall29/ itrsh(itlmax)
      common /tall31/ iterl(itlmax)

      common /tall33/ itln(itlmax,2), itli(itlmax,2), itlr(itlmax,2),
     &                rtdm(itlmax,2)
      common /tall34/ itvm(itlmax)
      common /tall35/ iteps(itlmax)
      common /tall37/ itcnt(9,itlmax)

      common /cntreg/ ncntc(3), ncreg(3), incrc(3), incrt(3), kcont(3)
      common /cntpat/ icpan(3), icpat(3,20,2)

      common /tall38/ itres(itlmax)
      common /tall39/ rtwid(itlmax)
      common /tall46/ itmtr(itlmax,4), rtmtr(itlmax,13)
      common /tall49/ itglt(itlmax)

      common /tall60/ itger(itlmax)
      common /tall63/ itbmp(itlmax)
      common /tall64/ itvtk(itlmax),itvtkfmt(itlmax)

      common /tall67/ itstd(itlmax), rtstd(itlmax)

      integer :: itfoam
      common /tall76/ itfoam(itlmax)

      integer :: itdpa
      integer :: iteth
      common /tall80/ itdpa(itlmax)
      common /tall81/ iteth(itlmax)

      common /tall82/ itcnth(9,itlmax)

      dimension idas(mdas*2) !20220906frtati temporary for fbounds-check
      equivalence ( das, idas )


*-----------------------------------------------------------------------

      common /inpec/ ititl, ipara, ibody, iregn, llarr, itby, itar

      common /kmat1d/ idmn(0:kvlmax), idnm(kvmmax)

*-----------------------------------------------------------------------
       integer   icntl
       integer   inucr

       common /tcntl/  icntl, inucr
!-----------------------------------------------------------------------

      character asfil*100
      character aname(10)*3
      data aname / 'reg','  x','  y','  z','  r',
     &             ' xy',' yz',' xz',' rz','tet'/

      dimension ian(6)
      dimension jmat(20)
      character chau*8
      character chaus(7)*8
      character elmnt(104)*3

*-----------------------------------------------------------------------

      data elmnt /
     &    'H  ','He ','Li ','Be ','B  ','C  ','N  ','O  ','F  ','Ne ',
     &    'Na ','Mg ','Al ','Si ','P  ','S  ','Cl ','Ar ','K  ','Ca ',
     &    'Sc ','Ti ','V  ','Cr ','Mn ','Fe ','Co ','Ni ','Cu ','Zn ',
     &    'Ga ','Ge ','As ','Se ','Br ','Kr ','Rb ','Sr ','Y  ','Zr ',
     &    'Nb ','Mo ','Tc ','Ru ','Rh ','Pd ','Ag ','Cd ','In ','Sn ',
     &    'Sb ','Te ','I  ','Xe ','Cs ','Ba ','La ','Ce ','Pr ','Nd ',
     &    'Pm ','Sm ','Eu ','Gd ','Tb ','Dy ','Ho ','Er ','Tm ','Yb ',
     &    'Lu ','Hf ','Ta ','W  ','Re ','Os ','Ir ','Pt ','Au ','Hg ',
     &    'Tl ','Pb ','Bi ','Po ','At ','Rn ','Fr ','Ra ','Ac ','Th ',
     &    'Pa ','U  ','Np ','Pu ','Am ','Cm ','Bk ','Cf ','Es ','Fm ',
     &    'Md ','No ','Lr ','Ku ' /

*-----------------------------------------------------------------------
      common /tall79/ itallech
      if( itallech .eq. 0 ) return

*-----------------------------------------------------------------------
*     input echo
*-----------------------------------------------------------------------

            write(iot,'("[ T-DPA ]")')

                  write(iot,'("    title = ",80a1)')
     &            ( ittle(m)(i:i), i = 1, ittll(m) )

*-----------------------------------------------------------------------

            if( itmsh(m) .eq. 1 ) then

                  write(iot,'("     mesh =  reg",11x,
     &            " # mesh type is region-wise")')

                     idas1 = mmmax
                     idas2 = ( idas1 + itrgn(m) - 1 ) * 2 + 1
                     idas3 = idas2 + itrgn(m)

                  call echrg(1,iterl(m),0,1,iot,itrcm(m),
     &                       idas_itrcg(itrcg(m)),
     &                       itrgn(m),itrgm(m),idas_itreg(itreg(m)),
     &                       das(idas1),idas(idas2),
     &                       itrnv(m),idas(itriv(m)),das(itrrv(m)),
     &                       idas3)

            else if( itmsh(m) .eq. 2 ) then

                  write(iot,'("     mesh =  r-z",11x,
     &            " # mesh type is r-z scoring mesh")')

                  write(iot,'("       x0 = ",1p1g14.7,1x,
     &            " # center x-position of r-z mesh")') rtrx0(m)

                  write(iot,'("       y0 = ",1p1g14.7,1x,
     &            " # center y-position of r-z mesh")') rtry0(m)

                  call echmty(0,iot,'r',itrty(m),itrnm(m),
     &                        rtrdl(m),rtrmi(m),rtrma(m),itrrg(m)
     &                        ,abs(itrnm(m)+1),das_itrrg)
                  call echmty(0,iot,'z',itzty(m),itznm(m),
     &                        rtzdl(m),rtzmi(m),rtzma(m),itzrg(m)
     &                        ,abs(itznm(m)+1),das_itzrg)

            else if( itmsh(m) .eq. 3 ) then

                  write(iot,'("     mesh =  xyz",11x,
     &            " # mesh type is xyz scoring mesh")')

                  call echmty(0,iot,'x',itxty(m),itxnm(m),
     &                        rtxdl(m),rtxmi(m),rtxma(m),itxrg(m)
     &                        ,abs(itxnm(m)+1),das_itxrg)
                  call echmty(0,iot,'y',ityty(m),itynm(m),
     &                        rtydl(m),rtymi(m),rtyma(m),ityrg(m)
     &                        ,abs(itynm(m)+1),das_ityrg)
                  call echmty(0,iot,'z',itzty(m),itznm(m),
     &                        rtzdl(m),rtzmi(m),rtzma(m),itzrg(m)
     &                        ,abs(itznm(m)+1),das_itzrg)

            else if( itmsh(m) .eq. 4 ) then

                  write(iot,'("     mesh =  tet",11x,
     &            " # mesh type is tetra geometry")')

                  call echtet(iot,itrgm(m),idas_itreg(itreg(m)))

            end if

*-----------------------------------------------------------------------

            if( itmtr(m,1) .gt. 0 ) then

               if( itmtr(m,1) .eq. 1 ) then

                     write(iot,'("     trcl = ",i6,9x,
     &               " # transform ID")') itmtr(m,3)

               else

                  if( itmtr(m,2) .eq. 0 ) then

                     write(iot,'("     trcl = ",
     &                           3(1p1e15.7),
     &                           3(/12x,3(1p1e15.7)),i5)')
     &               ( rtmtr(m,j), j = 1, 12 ),
     &               nint( rtmtr(m,13) )

                  else

                     write(iot,'("    *trcl = ",
     &                           3(1p1e15.7),
     &                           3(/12x,3(1p1e15.7)),i5)')
     &               ( rtmtr(m,j), j = 1, 12 ),
     &               nint( rtmtr(m,13) )

                  end if

               end if

            end if

*-----------------------------------------------------------------------

               if( itmtn(m) .eq. 0 ) then

                  write(iot,'(" material =  all"11x,
     &            " # (D=all) number of specific material")')

               else if( itmtn(m) .gt. 0 ) then

                  write(iot,'(" material =",i5,11x,
     &            " # (D=all) number of specific material")')
     &            itmtn(m) * itmcn(m)

                     k = 0

                  do j = 1, itmtn(m)

                     k = k + 1

                     jmat(k) = ismte( itmtt(m) + j - 1 )

                     if( k / 10 * 10 .eq. k .or. j .eq. itmtn(m) ) then

                        write(iot,'(11x,10(i6))') ( jmat(i),i = 1, k )

                        k = 0

                     end if

                  end do

               end if

*-----------------------------------------------------------------------

               if( itln(m,1) .gt. 0 ) then

                  write(iot,'("  library = ",i4,11x,
     &            " # (D=0) number DPA library for proton")')
     &            itln(m,1)

                  write(iot,'("     part = proton")')
                  write(iot,'("     emax = ",1p1g13.5)') rtdm(m,1)

                     jtli = itli(m,1)
                     jtlr = itlr(m,1)

                  write(iot,'(3x,"   mat",4x,"fac"6x,
     &                        3x,"lib",4x,"mt")')

                  do k = 1, itln(m,1)

                     write(iot,'(3x,i6,1p1g13.5,2i6)')
     &               mlib(jtli-1+1,k),
     &               flib(jtlr-1+k),
     &               mlib(jtli-1+2,k),
     &               mlib(jtli-1+3,k)

                  end do

               end if

*-----------------------------------------------------------------------

               if( itln(m,2) .gt. 0 ) then

                  write(iot,'("  library = ",i4,11x,
     &            " # (D=0) number DPA library for neutron")')
     &            itln(m,2)

                  write(iot,'("     part = neutron")')
                  write(iot,'("     emax = ",1p1g13.5)') rtdm(m,2)

                     jtli = itli(m,2)
                     jtlr = itlr(m,2)

                  write(iot,'(3x,"   mat",4x,"fac"6x,
     &                        3x,"lib",4x,"mt")')

                  do k = 1, itln(m,2)

                     write(iot,'(3x,i6,1p1g13.5,2i6)')
     &               mlib(jtli-1+1,k),
     &               flib(jtlr-1+k),
     &               mlib(jtli-1+2,k),
     &               mlib(jtli-1+3,k)

                  end do

               end if

*-----------------------------------------------------------------------

               if( itman(m) .eq. 0 ) then

                  write(iot,'("   mother =  all"11x,
     &            " # (D=all) number of specific mother")')

               else if( itman(m) .gt. 0 ) then

                  write(iot,'("   mother = ",i4,11x,
     &            " # (D=all) number of specific mother")')
     &            itman(m) * itmct(m)

                     k = 0

                  do j = 1, itman(m)

                     k = k + 1

                     iaz = ismat( itmat(m) + j - 1 )

                     iz = iaz / 1000
                     ia = iaz - iz * 1000

                     if( ia .eq. 0 ) then

                        chaus(k) = elmnt(iz)

                     else

                        call chname(idum,ia,iz,chau)

                        chaus(k) = chau

                     end if

                     if( k / 7 * 7 .eq. k .or. j .eq. itman(m) ) then

                        write(iot,'(12x,7(a8))') ( chaus(i),i = 1, k )

                        k = 0

                     end if

                  end do

               end if

*-----------------------------------------------------------------------

               if( itunt(m) .eq. 1 ) then

                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [DPA*1.E-24/source]")') itunt(m)

               else if( itunt(m) .eq. 2 ) then

                  write(iot,'("     unit = ",i4,11x,
     &            " # unit is [DPA/source]")') itunt(m)

               end if

*-----------------------------------------------------------------------

                  asfil = '  # file name of output for the above axis'

            if( icax .eq. 0 ) then

                  icn = itaxn(m)

               do i = 1, icn

                  ian(i) = i

               end  do

            else

                  icn = 1
                  ian(1) = iaxs

            end if

*-----------------------------------------------------------------------

               idtyp = 0

            do i = 1, icn

               if( itaxs(m,ian(i)) .ge. 6 .and.
     &             itaxs(m,ian(i)) .le. 9 )  idtyp = ittwo(m)

            end do

            if( idtyp .ne. 0 ) then

               write(iot,'("  2D-type = ",i4,11x,
     &         " # 1:Cont, 2:Clust, 3:Color, 4:xyz, 5:mat,",
     &         " 6:Clust+Cont, 7:Col+Cont")') idtyp

            end if

*-----------------------------------------------------------------------

            do i = 1, icn

                  iax = ian(i)

                  write(iot,'("     axis =  ",a3,11x,
     &            " # axis of output")')
     &            aname( itaxs(m,iax) )

             if ( itfll(m,iax) .gt. 0 ) then

                  msfile = max( 14, itfll(m,iax) )

                  write(iot,'("     file = ",100a1)')
     &                  ( ctfln(m,iax)(j:j), j = 1, msfile ),
     &                  ( asfil(j:j), j = 1, 42 )

             end if

            end do

*-----------------------------------------------------------------------

               if( itout(m) .eq. 1 ) then

                  write(iot,'("   output =  dpa",10x,
     &            "  # (D=dpa) toatl DPA")')

               else if( itout(m) .eq. 2 ) then

                  write(iot,'("   output = simple",8x,
     &           "  # (D=dpa) total, elastic, inelastic, stopped")')

               else if( itout(m) .eq. 3 ) then

                  write(iot,'("   output = all",10x,
     &            "  # (D=dpa) DPA of simple plus fragments")')

               end if

*-----------------------------------------------------------------------
               call echprt(1,iot,m,itpan,itpat,jtpat,itnm,itmxpt,itmxgp)

*-----------------------------------------------------------------------

               if( rtfac(m) .ne. 1.0d0 ) then

                  write(iot,'("   factor =",1p1g14.7,1x,1x,
     &            " # (D=1.0) normalization factor")')
     &            rtfac(m)

               end if

*-----------------------------------------------------------------------

               if( itmtn(m) .gt. 0 .and. itmsh(m) .eq. 3 ) then

                  write(iot,'("   volmat = ",i4,11x,
     &            " # (D=9) vol.corr., 0:no, number of scan / dir.")')
     &            itvm(m)

               end if

*-----------------------------------------------------------------------

               if( itanl(m) .gt. 0 ) then

                  write(iot,'("    angel = ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               end if

*-----------------------------------------------------------------------

               if( itsans(m) .gt. 0 ) then

                  write(iot,'("   sangel = ",i2)') itsans(m)
                  call write_sangel(iot,m,1)

               end if

*-----------------------------------------------------------------------

               if( itaxl(m) .gt. 0 ) then

                  write(iot,'("    x-txt = ",200a1)')
     &            ( itaxt(m)(i:i),i = 1, itaxl(m) )

               end if

*-----------------------------------------------------------------------

               if( itayl(m) .gt. 0 ) then

                  write(iot,'("    y-txt = ",200a1)')
     &            ( itayt(m)(i:i),i = 1, itayl(m) )

               end if

*-----------------------------------------------------------------------

               if( itazl(m) .gt. 0 ) then

                  write(iot,'("    z-txt = ",200a1)')
     &            ( itazt(m)(i:i),i = 1, itazl(m) )

               end if

*-----------------------------------------------------------------------

               if( itgsh(m) .gt. 0 ) then

                  write(iot,'("    gshow = ",i4,10x,
     &            "  # 0: no 1:bnd, 2:bnd+mat, ",
     &     "3:bnd+reg 4:bnd+lat 5:bmp style")')
     &            itgsh(m)

               end if

*-----------------------------------------------------------------------

               if( itglt(m) .ge. 0 ) then

                  write(iot,'("    gslat = ",i4,10x,
     &            "  # 0: no line for lat, 1: with,",
     &               " 2: no line between the same cell")')
     &            itglt(m)

               end if

*-----------------------------------------------------------------------

               if( itgsh(m) .gt. 0 .and. itger(m) .gt. 0 ) then

                  write(iot,'("    ginfo = ",i4,11x,
     &            " # (D=0) 0:no, 1:show in graph, 2:add file")')
     &            itger(m)

               end if

*-----------------------------------------------------------------------

               if( itrsh(m) .gt. 0 ) then

                  write(iot,'("    rshow = ",i4,10x,
     &            "  # 0: no 1:bnd, 2:bnd+mat, ",
     &     "3:bnd+reg 4:bnd+lat 5:bmp style")')
     &            itrsh(m)

                  call echmty(1,iot,'x',itxty(m),itxnm(m),
     &                        rtxdl(m),rtxmi(m),rtxma(m),itxrg(m)
     &                        ,abs(itxnm(m)+1),das_itxrg)
                  call echmty(1,iot,'y',ityty(m),itynm(m),
     &                        rtydl(m),rtymi(m),rtyma(m),ityrg(m)
     &                        ,abs(itynm(m)+1),das_ityrg)
                  call echmty(1,iot,'z',itzty(m),itznm(m),
     &                        rtzdl(m),rtzmi(m),rtzma(m),itzrg(m)
     &                        ,abs(itznm(m)+1),das_itzrg)

               end if

*-----------------------------------------------------------------------

               if( itgsh(m) .gt. 0 .or. itrsh(m) .gt. 0 ) then

                  write(iot,'("    resol = ",i4,10x,
     &            "  # (D=1) resolution of gshow or rshow")')
     &            itres(m)

                  write(iot,'("    width =",1p1g14.7,1x,1x,
     &            " # (D=0.5) width of lines for gshow or rshow")')
     &            rtwid(m)

               end if

*-----------------------------------------------------------------------

               if( iterl(m) .ne. 72 ) then

                  write(iot,'("   iechrl =",i5,10x,
     &            "  # (D=72) length reg of val or vol")')
     &            iterl(m)

               end if

*-----------------------------------------------------------------------

               if( iteps(m) .ne. 0 ) then

                  write(iot,'("   epsout =",i5,10x,
     &            "  # (D=0) generate eps file by ANGEL")')
     &            iteps(m)

               end if

*-----------------------------------------------------------------------

               if( itbmp(m) .ne. 0 ) then

                  write(iot,'("   bmpout =",i5,10x,
     &            "  # (D=0) generate bitmap file")')
     &            itbmp(m)

               end if

*-----------------------------------------------------------------------

               if( itvtk(m) .ne. 0 ) then

                  write(iot,'("   vtkout =",i5,10x,
     &            "  # (D=0) generate vtk file")')
     &            itvtk(m)

                  if ( itvtkfmt(m) .eq. 0 ) then

                     write(iot,'("   vtkfmt =",i5,10x,
     &               "  # (D=0) generate vtk ascii file")')
     &               itvtkfmt(m)

                  else

                     write(iot,'("   vtkfmt =",i5,10x,
     &               "  # (D=0) generate vtk binary file")')
     &               itvtkfmt(m)

                  end if

               end if

*-----------------------------------------------------------------------
*        check of counter
*-----------------------------------------------------------------------

            do i = 1, 3

               if( itcnt(i,m) .ne. 0 ) then

                  write(iot,'(" ctmin(",i1,") =",i5,10x,
     &            "  # (D=-9999) min. of counter ",i1)')
     &            i, itcnt(i*2+2,m), i

                  write(iot,'(" ctmax(",i1,") =",i5,10x,
     &            "  # (D= 9999) max. of counter ",i1)')
     &            i, itcnt(i*2+3,m), i

                  if( ncntc(i) .eq. 0 )
     &            write(iot,'("# Warning: there is no counter = ",
     &                        i1," definition in [counter].")') i

               end if

            end do

            do i = 1, 3
               if( itcnth(i,m) .ne. 0 ) then
                  write(iot,'(" chmin(",i1,") =",i5,10x,
     &            "  # (D=-9999) min. of history counter ",i1)')
     &            i, itcnth(i*2+2,m), i
                  write(iot,'(" chmax(",i1,") =",i5,10x,
     &            "  # (D= 9999) max. of history counter ",i1)')
     &            i, itcnth(i*2+3,m), i
               end if
            end do

*-----------------------------------------------------------------------

               if( rtstd(m) .gt. 0.d0 ) then

                  write(iot,'("   stdcut =",1p1g14.7,1x,1x,
     &            " # (D=-1) threshold value of STD cutoff")')
     &            rtstd(m)

               end if

*-----------------------------------------------------------------------
               if( itfoam(m) .ne. 0 ) then

                  write(iot,'("  foamout =",i5,10x,
     &            "  # (D=0) generate OpenFOAM data")')
     &            itfoam(m)

               end if

*-----------------------------------------------------------------------
               if( itdpa(m) .ne. 0 ) then

                  write(iot,'("     idpa =",i5,10x,
     &            "  # dpa with defect efficiency, arc-dpa")')
     &            itdpa(m)

               end if
*-----------------------------------------------------------------------
               if( itdpa(m) .eq. 0 ) then

                  write(iot,'("     idpa =",i5,10x,
     &          "  # (D=0) dpa without defect efficiency, NRT-dpa")')
     &            itdpa(m)

               end if
*-----------------------------------------------------------------------
               if( iteth(m) .ne. 0 ) then

                  write(iot,'("     ieth =",i5,10x,
     &            "  # new displacement threshold energy set")')
     &            iteth(m)

               end if
*-----------------------------------------------------------------------
               if( iteth(m) .eq. 0 ) then

                  write(iot,'("     ieth =",i5,10x,
     &       "  # (D=0) original displacement threshold energy set")')
     &            iteth(m)

               end if
*-----------------------------------------------------------------------
               if( itmxang(m) .ne. 0 ) then

                  write(iot,'(" maxangel =",i5,10x,
     &            "  # (D=0) number of output particles")')
     &            itmxang(m)

               end if

*-----------------------------------------------------------------------

*-----------------------------------------------------------------------
*        sumtally subsection
*-----------------------------------------------------------------------
            if ( icntl==13 .and. nsumtalRead(m)==7 ) then

               call sumtal_echo(iot,m)

            end if
*-----------------------------------------------------------------------

            write(iot,'("#    used :",
     &                  "        main (  %)",
     &                  "        temp (  %)",
     &                  "       total (  %)"/
     &                  "#  memory :",
     &                  3(i12," (",i3,")"))')
     &                  itsmn(m), nint(dble(itsmn(m))/dble(mdas)*100),
     &                  itstm(m), nint(dble(itstm(m))/dble(mdas)*100),
     &                  itsmn(m)+itstm(m),
     &                  nint(dble(itsmn(m)+itstm(m))/dble(mdas)*100)

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine tgshech(iot,m,iaxs,icax)
*                                                                      *
*       input echo of the g-show tally                                 *
*       last modified by K.Niita on 2003/09/13                         *
*                                                                      *
************************************************************************
      use sangelmod, only: itsans
      use moddas_mesh

*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'

*-----------------------------------------------------------------------

      common /tall00/ itnm, ital(itlmax), itals(itlmax), italm(itlmax)
      common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                itout(itlmax), ittwo(itlmax)

      common /tall04/ itxty(itlmax), itxnm(itlmax), itxrg(itlmax),
     &                rtxmi(itlmax), rtxma(itlmax), rtxdl(itlmax)
      common /tall05/ ityty(itlmax), itynm(itlmax), ityrg(itlmax),
     &                rtymi(itlmax), rtyma(itlmax), rtydl(itlmax)
      common /tall06/ itzty(itlmax), itznm(itlmax), itzrg(itlmax),
     &                rtzmi(itlmax), rtzma(itlmax), rtzdl(itlmax)

      common /tall10/ itaxn(itlmax), itaxs(itlmax,6)
      common /tall11/ itfln(itlmax), itfll(itlmax,6), ctfln(itlmax,6)
      character ctfln*100

      common /tall14/ ittll(itlmax), ittle(itlmax)
      character ittle*80

      common /tall15/ itanl(itlmax), itang(itlmax)
      character itang*200

      common /tall27/ itaxl(itlmax), itaxt(itlmax),
     &                itayl(itlmax), itayt(itlmax),
     &                itazl(itlmax), itazt(itlmax)
      character itaxt*200, itayt*200, itazt*200

      common /tall19/ itsmn(itlmax), itstm(itlmax)

      common /tall21/ rtfac(itlmax)
      common /tall35/ iteps(itlmax)

      common /tall38/ itres(itlmax)
      common /tall39/ rtwid(itlmax)
      common /tall46/ itmtr(itlmax,4), rtmtr(itlmax,13)
      common /tall49/ itglt(itlmax)
      common /tall52/ itprd(itlmax)

      common /tall60/ itger(itlmax)
      common /tall63/ itbmp(itlmax)
      common /tall64/ itvtk(itlmax),itvtkfmt(itlmax)

*-----------------------------------------------------------------------

      character asfil*100

      character aname(3)*3
      data aname / ' xy',' yz',' xz'/

      dimension ian(6)

*-----------------------------------------------------------------------
      common /tall79/ itallech
      if( itallech .eq. 0 ) return

*-----------------------------------------------------------------------
*     input echo
*-----------------------------------------------------------------------

            write(iot,'("[ T-Gshow ]")')

                  write(iot,'("    title = ",80a1)')
     &            ( ittle(m)(i:i), i = 1, ittll(m) )

*-----------------------------------------------------------------------

                  write(iot,'("     mesh =  xyz",11x,
     &            " # mesh type is xyz scoring mesh")')

                  call echmty(1,iot,'x',itxty(m),itxnm(m),
     &                        rtxdl(m),rtxmi(m),rtxma(m),itxrg(m)
     &                        ,abs(itxnm(m)+1),das_itxrg)
                  call echmty(1,iot,'y',ityty(m),itynm(m),
     &                        rtydl(m),rtymi(m),rtyma(m),ityrg(m)
     &                        ,abs(itynm(m)+1),das_ityrg)
                  call echmty(1,iot,'z',itzty(m),itznm(m),
     &                        rtzdl(m),rtzmi(m),rtzma(m),itzrg(m)
     &                        ,abs(itznm(m)+1),das_itzrg)

*-----------------------------------------------------------------------

            if( itmtr(m,1) .gt. 0 ) then

               if( itmtr(m,1) .eq. 1 ) then

                     write(iot,'("     trcl = ",i6,9x,
     &               " # transformation ID")') itmtr(m,3)

               else

                  if( itmtr(m,2) .eq. 0 ) then

                     write(iot,'("     trcl = ",
     &                           3(1p1e15.7),
     &                           3(/12x,3(1p1e15.7)),i5)')
     &               ( rtmtr(m,j), j = 1, 12 ),
     &               nint( rtmtr(m,13) )

                  else

                     write(iot,'("    *trcl = ",
     &                           3(1p1e15.7),
     &                           3(/12x,3(1p1e15.7)),i5)')
     &               ( rtmtr(m,j), j = 1, 12 ),
     &               nint( rtmtr(m,13) )

                  end if

               end if

            end if

*-----------------------------------------------------------------------

                  asfil = '  # file name of output for the above axis'

            if( icax .eq. 0 ) then

                  icn = itaxn(m)

               do i = 1, icn

                  ian(i) = i

               end  do

            else

                  icn = 1
                  ian(1) = iaxs

            end if

*-----------------------------------------------------------------------

            do i = 1, icn

                  iax = ian(i)

                  write(iot,'("     axis =  ",a3,11x,
     &            " # axis of output")')
     &            aname( itaxs(m,iax) )

                  msfile = max( 14, itfll(m,iax) )

                  write(iot,'("     file = ",100a1)')
     &                  ( ctfln(m,iax)(j:j), j = 1, msfile ),
     &                  ( asfil(j:j), j = 1, 42 )

            end do

*-----------------------------------------------------------------------

                  write(iot,'("   output = ",i4,11x,
     &            " # (D=2) 1:bnd, 2:bnd+mat, ",
     &                     "3:bnd+num 4:bnd+mat+num")')
     &            itout(m)

*-----------------------------------------------------------------------
                  write(iot,'("    ginfo = ",i4,11x,
     &            " # (D=0) 0:no, 1:show in graph, 2:add file")')
     &            itger(m)

*-----------------------------------------------------------------------

                  write(iot,'("    resol = ",i4,10x,
     &            "  # (D=1) resolution of gshow or rshow")')
     &            itres(m)

                  write(iot,'("    width =",1p1g14.7,1x,1x,
     &            " # (D=0.5) width of lines for gshow or rshow")')
     &            rtwid(m)

*-----------------------------------------------------------------------

               if( itglt(m) .ge. 0 ) then

                  write(iot,'("    gslat = ",i4,10x,
     &            "  # 0: no line for lat, 1: with,",
     &               " 2: no line between the same cell")')
     &            itglt(m)

               end if

*-----------------------------------------------------------------------

               if( itanl(m) .gt. 0 ) then

                  write(iot,'("    angel = ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               end if

*-----------------------------------------------------------------------

               if( itsans(m) .gt. 0 ) then

                  write(iot,'("   sangel = ",i2)') itsans(m)
                  call write_sangel(iot,m,1)

               end if

*-----------------------------------------------------------------------

               if( itaxl(m) .gt. 0 ) then

                  write(iot,'("    x-txt = ",200a1)')
     &            ( itaxt(m)(i:i),i = 1, itaxl(m) )

               end if

*-----------------------------------------------------------------------

               if( itayl(m) .gt. 0 ) then

                  write(iot,'("    y-txt = ",200a1)')
     &            ( itayt(m)(i:i),i = 1, itayl(m) )

               end if

*-----------------------------------------------------------------------

               if( itazl(m) .gt. 0 ) then

                  write(iot,'("    z-txt = ",200a1)')
     &            ( itazt(m)(i:i),i = 1, itazl(m) )

               end if

*-----------------------------------------------------------------------

               if( iteps(m) .ne. 0 ) then

                  write(iot,'("   epsout =",i5,10x,
     &            "  # (D=0) generate eps file by ANGEL")')
     &            iteps(m)

               end if

*-----------------------------------------------------------------------

               if( itbmp(m) .ne. 0 ) then

                  write(iot,'("   bmpout =",i5,10x,
     &            "  # (D=0) generate bitmap file")')
     &            itbmp(m)

               end if

*-----------------------------------------------------------------------

               if( itvtk(m) .ne. 0 ) then

                  write(iot,'("   vtkout =",i5,10x,
     &            "  # (D=0) generate vtk file")')
     &            itvtk(m)

                  if ( itvtkfmt(m) .eq. 0 ) then

                     write(iot,'("   vtkfmt =",i5,10x,
     &               "  # (D=0) generate vtk ascii file")')
     &               itvtkfmt(m)

                  else

                     write(iot,'("   vtkfmt =",i5,10x,
     &               "  # (D=0) generate vtk binary file")')
     &               itvtkfmt(m)

                  end if

               end if

*-----------------------------------------------------------------------

            write(iot,'("#    used :",
     &                  "        main (  %)",
     &                  "        temp (  %)",
     &                  "       total (  %)"/
     &                  "#  memory :",
     &                  3(i12," (",i3,")"))')
     &                  itsmn(m), nint(dble(itsmn(m))/dble(mdas)*100),
     &                  itstm(m), nint(dble(itstm(m))/dble(mdas)*100),
     &                  itsmn(m)+itstm(m),
     &                  nint(dble(itsmn(m)+itstm(m))/dble(mdas)*100)

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine trshech(iot,m,iaxs,icax)
*                                                                      *
*       input echo of the r-show tally                                 *
*       last modified by K.Niita on 2003/08/13                         *
*                                                                      *
************************************************************************
      use sangelmod, only: itsans
      use moddas_mesh
      use moddas_region

*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'

*-----------------------------------------------------------------------

      common /tall00/ itnm, ital(itlmax), itals(itlmax), italm(itlmax)
      common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                itout(itlmax), ittwo(itlmax)

      common /tall04/ itxty(itlmax), itxnm(itlmax), itxrg(itlmax),
     &                rtxmi(itlmax), rtxma(itlmax), rtxdl(itlmax)
      common /tall05/ ityty(itlmax), itynm(itlmax), ityrg(itlmax),
     &                rtymi(itlmax), rtyma(itlmax), rtydl(itlmax)
      common /tall06/ itzty(itlmax), itznm(itlmax), itzrg(itlmax),
     &                rtzmi(itlmax), rtzma(itlmax), rtzdl(itlmax)

      common /tall10/ itaxn(itlmax), itaxs(itlmax,6)
      common /tall11/ itfln(itlmax), itfll(itlmax,6), ctfln(itlmax,6)
      character ctfln*100

      common /tall14/ ittll(itlmax), ittle(itlmax)
      character ittle*80

      common /tall15/ itanl(itlmax), itang(itlmax)
      character itang*200

      common /tall27/ itaxl(itlmax), itaxt(itlmax),
     &                itayl(itlmax), itayt(itlmax),
     &                itazl(itlmax), itazt(itlmax)
      character itaxt*200, itayt*200, itazt*200

      common /tall19/ itsmn(itlmax), itstm(itlmax)

      common /tall21/ rtfac(itlmax)

      common /tall02/ itrgn(itlmax), itrgm(itlmax), itreg(itlmax)
      common /tall24/ itrcm(itlmax), itrcg(itlmax)
      common /tall25/ itrnv(itlmax), itriv(itlmax), itrrv(itlmax)
      common /tall31/ iterl(itlmax)
      common /tall35/ iteps(itlmax)

      common /tall38/ itres(itlmax)
      common /tall39/ rtwid(itlmax)
      common /tall46/ itmtr(itlmax,4), rtmtr(itlmax,13)
      common /tall49/ itglt(itlmax)

      common /tall60/ itger(itlmax)
      common /tall63/ itbmp(itlmax)

      dimension idas(mdas*2) !20220906frtati temporary for fbounds-check
      equivalence ( das, idas )

*-----------------------------------------------------------------------

      character asfil*100

      character aname(3)*3
      data aname / ' xy',' yz',' xz'/

      dimension ian(6)

*-----------------------------------------------------------------------
      common /tall79/ itallech
      if( itallech .eq. 0 ) return

*-----------------------------------------------------------------------
*     input echo
*-----------------------------------------------------------------------

            write(iot,'("[ T-Rshow ]")')

                  write(iot,'("    title = ",80a1)')
     &            ( ittle(m)(i:i), i = 1, ittll(m) )

*-----------------------------------------------------------------------

                  write(iot,'("     mesh =  xyz",11x,
     &            " # mesh type is xyz scoring mesh")')

                  call echmty(1,iot,'x',itxty(m),itxnm(m),
     &                        rtxdl(m),rtxmi(m),rtxma(m),itxrg(m)
     &                        ,abs(itxnm(m)+1),das_itxrg)
                  call echmty(1,iot,'y',ityty(m),itynm(m),
     &                        rtydl(m),rtymi(m),rtyma(m),ityrg(m)
     &                        ,abs(itynm(m)+1),das_ityrg)
                  call echmty(1,iot,'z',itzty(m),itznm(m),
     &                        rtzdl(m),rtzmi(m),rtzma(m),itzrg(m)
     &                        ,abs(itznm(m)+1),das_itzrg)

*-----------------------------------------------------------------------

            if( itmtr(m,1) .gt. 0 ) then

               if( itmtr(m,1) .eq. 1 ) then

                     write(iot,'("     trcl = ",i6,9x,
     &               " # transformation ID")') itmtr(m,3)

               else

                  if( itmtr(m,2) .eq. 0 ) then

                     write(iot,'("     trcl = ",
     &                           3(1p1e15.7),
     &                           3(/12x,3(1p1e15.7)),i5)')
     &               ( rtmtr(m,j), j = 1, 12 ),
     &               nint( rtmtr(m,13) )

                  else

                     write(iot,'("    *trcl = ",
     &                           3(1p1e15.7),
     &                           3(/12x,3(1p1e15.7)),i5)')
     &               ( rtmtr(m,j), j = 1, 12 ),
     &               nint( rtmtr(m,13) )

                  end if

               end if

            end if

*-----------------------------------------------------------------------

                  asfil = '  # file name of output for the above axis'

            if( icax .eq. 0 ) then

                  icn = itaxn(m)

               do i = 1, icn

                  ian(i) = i

               end  do

            else

                  icn = 1
                  ian(1) = iaxs

            end if

*-----------------------------------------------------------------------

            do i = 1, icn

                  iax = ian(i)

                  write(iot,'("     axis =  ",a3,11x,
     &            " # axis of output")')
     &            aname( itaxs(m,iax) )

                  msfile = max( 14, itfll(m,iax) )

                  write(iot,'("     file = ",100a1)')
     &                  ( ctfln(m,iax)(j:j), j = 1, msfile ),
     &                  ( asfil(j:j), j = 1, 42 )

            end do

*-----------------------------------------------------------------------

                     idas1 = mmmax
                     idas2 = ( idas1 + itrgn(m) - 1 ) * 2 + 1
                     idas3 = idas2 + itrgn(m)

                  call echrg(1,iterl(m),1,1,iot,itrcm(m),
     &                       idas_itrcg(itrcg(m)),
     &                       itrgn(m),itrgm(m),idas_itreg(itreg(m)),
     &                       das(idas1),idas(idas2),
     &                       itrnv(m),idas(itriv(m)),das(itrrv(m)),
     &                       idas3)

*-----------------------------------------------------------------------

                  write(iot,'("   output = ",i4,11x,
     &            " # (D=1) 1:bnd, 2:bnd+num")')
     &            itout(m)

*-----------------------------------------------------------------------

                  write(iot,'("    resol = ",i4,10x,
     &            "  # (D=1) resolution of gshow or rshow")')
     &            itres(m)

                  write(iot,'("    width =",1p1g14.7,1x,1x,
     &            " # (D=0.5) width of lines for gshow or rshow")')
     &            rtwid(m)

*-----------------------------------------------------------------------

               if( itglt(m) .ge. 0 ) then

                  write(iot,'("    gslat = ",i4,10x,
     &            "  # 0: no for lat, 1: with,",
     &               " 2: no between the same cell")')
     &            itglt(m)

               end if

*-----------------------------------------------------------------------

               if( itanl(m) .gt. 0 ) then

                  write(iot,'("    angel = ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               end if

*-----------------------------------------------------------------------

               if( itsans(m) .gt. 0 ) then

                  write(iot,'("   sangel = ",i2)') itsans(m)
                  call write_sangel(iot,m,1)

               end if

*-----------------------------------------------------------------------

               if( itaxl(m) .gt. 0 ) then

                  write(iot,'("    x-txt = ",200a1)')
     &            ( itaxt(m)(i:i),i = 1, itaxl(m) )

               end if

*-----------------------------------------------------------------------

               if( itayl(m) .gt. 0 ) then

                  write(iot,'("    y-txt = ",200a1)')
     &            ( itayt(m)(i:i),i = 1, itayl(m) )

               end if

*-----------------------------------------------------------------------

               if( itazl(m) .gt. 0 ) then

                  write(iot,'("    z-txt = ",200a1)')
     &            ( itazt(m)(i:i),i = 1, itazl(m) )

               end if

*-----------------------------------------------------------------------

               if( iterl(m) .ne. 72 ) then

                  write(iot,'("   iechrl =",i5,10x,
     &            "  # (D=72) length reg of val or vol")')
     &            iterl(m)

               end if

*-----------------------------------------------------------------------

               if( iteps(m) .ne. 0 ) then

                  write(iot,'("   epsout =",i5,10x,
     &            "  # (D=0) generate eps file by ANGEL")')
     &            iteps(m)

               end if

*-----------------------------------------------------------------------

               if( itbmp(m) .ne. 0 ) then

                  write(iot,'("   bmpout =",i5,10x,
     &            "  # (D=0) generate bitmap file")')
     &            itbmp(m)

               end if

*-----------------------------------------------------------------------
                  write(iot,'("    ginfo = ",i4,11x,
     &            " # (D=0) 0:no, 1:show in graph, 2:add file")')
     &            itger(m)

*-----------------------------------------------------------------------

            write(iot,'("#    used :",
     &                  "        main (  %)",
     &                  "        temp (  %)",
     &                  "       total (  %)"/
     &                  "#  memory :",
     &                  3(i12," (",i3,")"))')
     &                  itsmn(m), nint(dble(itsmn(m))/dble(mdas)*100),
     &                  itstm(m), nint(dble(itstm(m))/dble(mdas)*100),
     &                  itsmn(m)+itstm(m),
     &                  nint(dble(itsmn(m)+itstm(m))/dble(mdas)*100)

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine tdshech(iot,m,iaxs,icax)
*                                                                      *
*       input echo of the 3d-show tally                                *
*       last modified by K.Niita on 2002/10/23                         *
*                                                                      *
************************************************************************
      use sangelmod, only: itsans
      use moddas_region
      use moddas_tally

*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'

*-----------------------------------------------------------------------

      common /tall00/ itnm, ital(itlmax), itals(itlmax), italm(itlmax)
      common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                itout(itlmax), ittwo(itlmax)

      common /tall11/ itfln(itlmax), itfll(itlmax,6), ctfln(itlmax,6)
      character ctfln*100

      common /tall14/ ittll(itlmax), ittle(itlmax)
      character ittle*80

      common /tall15/ itanl(itlmax), itang(itlmax)
      character itang*200

      common /tall02/ itrgn(itlmax), itrgm(itlmax), itreg(itlmax)
      common /tall25/ itrnv(itlmax), itriv(itlmax), itrrv(itlmax)

      common /tall22/ itmcn(itlmax), itmtn(itlmax), itmtt(itlmax)
      common /tall24/ itrcm(itlmax), itrcg(itlmax)
      common /tall27/ itaxl(itlmax), itaxt(itlmax),
     &                itayl(itlmax), itayt(itlmax),
     &                itazl(itlmax), itazt(itlmax)
      character itaxt*200, itayt*200, itazt*200

      common /tall19/ itsmn(itlmax), itstm(itlmax)
      common /tall35/ iteps(itlmax)

      common /tall38/ itres(itlmax)
      common /tall39/ rtwid(itlmax)

      common /tall40/ rtorg(itlmax,3), rteye(itlmax,3), rtlit(itlmax,3),
     &                rtwin(itlmax,3), rtbrt(itlmax,2), ithvn(itlmax),
     &                itwin(itlmax,2), itbox(itlmax), itmir(itlmax),
     &                rtbox(itlmax,5,10), rtout(itlmax), rthet(itlmax),
     &                itlin(itlmax), itshd(itlmax), itgxs(itlmax)

      common /tall42/ itmbn(itlmax), itmbt(itlmax)
      common /tall43/ itmgn(itlmax), itmgm(itlmax), itmeg(itlmax)
      common /tall44/ itmcm(itlmax), itmcg(itlmax)
      common /tall47/ itbtr(itlmax,5,4), rtbtr(itlmax,5,13)

*-----------------------------------------------------------------------

      character asfil*100

      dimension jmat(20)
      character aname(6)*3
      data aname / '  x',' -x','  y',' -y','  z',' -z'/

      dimension idas(1)
      equivalence ( das, idas )

*-----------------------------------------------------------------------
      common /tall79/ itallech
      if( itallech .eq. 0 ) return

*-----------------------------------------------------------------------
*     input echo
*-----------------------------------------------------------------------

            write(iot,'("[ T-3Dshow ]")')

                  write(iot,'("    title = ",80a1)')
     &            ( ittle(m)(i:i), i = 1, ittll(m) )

                  iax = 1

                  msfile = max( 14, itfll(m,iax) )
                  asfil = '  # file name of output'

                  write(iot,'("     file = ",100a1)')
     &                  ( ctfln(m,iax)(j:j), j = 1, msfile ),
     &                  ( asfil(j:j), j = 1, 23 )

*-----------------------------------------------------------------------

               if( itmtn(m) .eq. 0 ) then

                  write(iot,'(" material =  all"11x,
     &            " # (D=all) all materials are shown")')

               else if( itmtn(m) .gt. 0 ) then

                  write(iot,'(" material =",i5,11x,
     &            " # (D=all) show(+) or transparent(-) material")')
     &            itmtn(m) * itmcn(m)

                     k = 0

                  do j = 1, itmtn(m)

                     k = k + 1

                     jmat(k) = ismte( itmtt(m) + j - 1 )

                     if( k / 10 * 10 .eq. k .or. j .eq. itmtn(m) ) then

                        write(iot,'(11x,10(i6))') ( jmat(i),i = 1, k )

                        k = 0

                     end if

                  end do

               end if

*-----------------------------------------------------------------------

               if( itrgn(m) .gt. 0 ) then

                 write(iot,'(28x,"# transparent or shown regions.")')

                     idas1 = mmmax
                     idas2 = ( idas1 + itrgn(m) - 1 ) * 2 + 1
                     idas3 = idas2 + itrgn(m)

                  call echrg(1,72,0,0,iot,itrcm(m),
     &                       idas_itrcg(itrcg(m)),
     &                       itrgn(m),itrgm(m),idas_itreg(itreg(m)),
     &                       das(idas1),idas(idas2),
     &                       itrnv(m),idas(itriv(m)),das(itrrv(m)),
     &                       idas3)

               end if

*-----------------------------------------------------------------------

                  write(iot,'("       x0 =",1p1g14.7,1x,1x,
     &            " # (D=0.0) x-coordinate of the origin")')
     &            rtorg(m,1)
                  write(iot,'("       y0 =",1p1g14.7,1x,1x,
     &            " # (D=0.0) y-coordinate of the origin")')
     &            rtorg(m,2)
                  write(iot,'("       z0 =",1p1g14.7,1x,1x,
     &            " # (D=0.0) z-coordinate of the origin")')
     &            rtorg(m,3)

*-----------------------------------------------------------------------

                  write(iot,'("    e-the =",1p1g14.7,1x,1x,
     &            " # (D=80.0) eye point theta(degree) from z-axis")')
     &            rteye(m,1)
                  write(iot,'("    e-phi =",1p1g14.7,1x,1x,
     &            " # (D=140.0) eye point phi(degree) from x-axis")')
     &            rteye(m,2)
                  write(iot,'("    e-dst =",1p1g14.7,1x,1x,
     &            " # (D=w-dst*10) eye point distance from origin")')
     &            rteye(m,3)

*-----------------------------------------------------------------------

                  write(iot,'("    l-the =",1p1g14.7,1x,1x,
     &            " # (D=e-the) light point theta from z-axis")')
     &            rtlit(m,1)
                  write(iot,'("    l-phi =",1p1g14.7,1x,1x,
     &            " # (D=e-phi) light point phi from x-axis")')
     &            rtlit(m,2)
                  write(iot,'("    l-dst =",1p1g14.7,1x,1x,
     &            " # (D=e-dst) light point distance from origin")')
     &            rtlit(m,3)

*-----------------------------------------------------------------------

                  write(iot,'("    w-wdt =",1p1g14.7,1x,1x,
     &            " # (D=100) width of window (cm)")')
     &            rtwin(m,1)
                  write(iot,'("    w-hgt =",1p1g14.7,1x,1x,
     &            " # (D=100) hight of window (cm)")')
     &            rtwin(m,2)
                  write(iot,'("    w-dst =",1p1g14.7,1x,1x,
     &            " # (D=200) window distance from origin")')
     &            rtwin(m,3)

*-----------------------------------------------------------------------

                  write(iot,'("    w-mnw = ",i4,11x,
     &            " # (D=100) mesh number of window width")')
     &            itwin(m,1)
                  write(iot,'("    w-mnh = ",i4,11x
     &            " # (D=100) mesh number of window hight")')
     &            itwin(m,2)

*-----------------------------------------------------------------------

                  write(iot,'("   heaven = ",a3,12x,
     &            " # (D=y) direction to heaven")')
     &            aname(ithvn(m))

*-----------------------------------------------------------------------

                  write(iot,'("    w-ang =",1p1g14.7,1x,1x,
     &            " # (D=0.0) angle of frame (degree)")')
     &            rthet(m)

*-----------------------------------------------------------------------

               if( itmir(m) .eq. -1 ) then

                  write(iot,'("   mirror = ",i4,11x,
     &            " # (D=0) =-1: mirror display")')
     &            itmir(m)

               end if

*-----------------------------------------------------------------------

                  write(iot,'("   bright =",1p1g14.7,1x,1x,
     &            " # (D=0.8) top brightness")')
     &            rtbrt(m,1)
                  write(iot,'("     dark =",1p1g14.7,1x,1x,
     &            " # (D=0.2) bottom darkness")')
     &            rtbrt(m,2)

*-----------------------------------------------------------------------

                  write(iot,'("  axishow = ",i4,11x,
     &            " # (D=1) =1: small axis is shown,",
     &                    " =2: large, =0: no")')
     &            itgxs(m)

*-----------------------------------------------------------------------

            if( itbox(m) .gt. 0 ) then

                  write(iot,'("      box = ",i4,11x,
     &            " # (D=0) number of boxes")')
     &            itbox(m)

               do i = 1, itbox(m)

                  if( itbtr(m,i,1) .eq. 0 ) then

                        write(iot,'("      box   ",1p3e15.7/
     &                              12x,1p3e15.7/12x,1p4e15.7)')
     &                  (rtbox(m,i,j),j=1,10)

                  else if( itbtr(m,i,1) .eq. 1 ) then
                        write(iot,'("      box   trcl = ",i6/
     &                              12x,1p3e15.7/
     &                              12x,1p3e15.7/12x,1p4e15.7)')
     &                  itbtr(m,i,3), (rtbox(m,i,j),j=1,10)

                  else

                     if( itbtr(m,i,2) .eq. 0 ) then

                        write(iot,'("      box   trcl = (",
     &                                  1p3e15.7/
     &                              20x,1p3e15.7/
     &                              20x,1p3e15.7/
     &                              20x,1p3e15.7,i5," )"/
     &                              12x,1p3e15.7/
     &                              12x,1p3e15.7/12x,1p4e15.7)')
     &                  ( rtbtr(m,i,j), j = 1, 12 ),
     &                    nint( rtbtr(m,i,13) ),
     &                  (rtbox(m,i,j),j=1,10)

                     else

                        write(iot,'("      box  *trcl = (",
     &                                  1p3e15.7/
     &                              20x,1p3e15.7/
     &                              20x,1p3e15.7/
     &                              20x,1p3e15.7,i5," )"/
     &                              12x,1p3e15.7/
     &                              12x,1p3e15.7/12x,1p4e15.7)')
     &                  ( rtbtr(m,i,j), j = 1, 12 ),
     &                    nint( rtbtr(m,i,13) ),
     &                  (rtbox(m,i,j),j=1,10)

                     end if

                  end if

               end do

*-----------------------------------------------------------------------

               if( itmbn(m) .gt. 0 ) then

                  write(iot,'(" matinbox =",i5,11x,
     &            " # shown materials in box")')
     &            itmbn(m)

                     k = 0

                  do j = 1, itmbn(m)

                     k = k + 1

                     jmat(k) = ismte_itmbt( itmbt(m) + j - 1 )

                     if( k / 10 * 10 .eq. k .or. j .eq. itmbn(m) ) then

                        write(iot,'(11x,10(i6))') ( jmat(i),i = 1, k )

                        k = 0

                     end if

                  end do

               end if

*-----------------------------------------------------------------------

               if( itmgn(m) .gt. 0 ) then

                 write(iot,'(28x,
     &                "# transparent or shown regions in box.")')

                     idas1 = mmmax
                     idas2 = ( idas1 + itmgn(m) - 1 ) * 2 + 1
                     idas3 = idas2 + itmgn(m)

                  call echrg(2,72,0,0,iot,itmcm(m),idas_itmcg(itmcg(m)),
     &                       itmgn(m),itmgm(m),idas_itmeg(itmeg(m)),
     &                       das(idas1),idas(idas2),
     &                       itrnv(m),idas(itriv(m)),das(itrrv(m)),
     &                       idas3)

               end if

*-----------------------------------------------------------------------

            end if

*-----------------------------------------------------------------------

                  write(iot,'("    resol = ",i4,10x,
     &            "  # (D=1) resolution of 3dshow")')
     &            itres(m)

                  write(iot,'("    width =",1p1g14.7,1x,1x,
     &            " # (D=0.5) width of lines for 3dshow")')
     &            rtwid(m)

*-----------------------------------------------------------------------

                  write(iot,'("    r-out =",1p1g14.7,1x,1x,
     &            " # (D=50000) radius of outer sphere")')
     &            rtout(m)

*-----------------------------------------------------------------------

                  write(iot,'("   output = ",i4,11x,
     &            " # (D=3) 0:draft, 1:line, 2:col, 3:line+col")')
     &            itout(m)

*-----------------------------------------------------------------------

                  write(iot,'("     line = ",i4,11x,
     &            " # (D=0) 0:surface+mat, 1:+region")')
     &            itlin(m)

*-----------------------------------------------------------------------

                  write(iot,'("   shadow = ",i4,11x,
     &            " # (D=0) 0:no, 1:shadow")')
     &            itshd(m)

*-----------------------------------------------------------------------

               if( itanl(m) .gt. 0 ) then

                  write(iot,'("    angel = ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               end if

*-----------------------------------------------------------------------

               if( itsans(m) .gt. 0 ) then

                  write(iot,'("   sangel = ",i2)') itsans(m)
                  call write_sangel(iot,m,1)

               end if

*-----------------------------------------------------------------------

               if( itaxl(m) .gt. 0 ) then

                  write(iot,'("    x-txt = ",200a1)')
     &            ( itaxt(m)(i:i),i = 1, itaxl(m) )

               end if

*-----------------------------------------------------------------------

               if( itayl(m) .gt. 0 ) then

                  write(iot,'("    y-txt = ",200a1)')
     &            ( itayt(m)(i:i),i = 1, itayl(m) )

               end if

*-----------------------------------------------------------------------

               if( itazl(m) .gt. 0 ) then

                  write(iot,'("    z-txt = ",200a1)')
     &            ( itazt(m)(i:i),i = 1, itazl(m) )

               end if

*-----------------------------------------------------------------------

               if( iteps(m) .ne. 0 ) then

                  write(iot,'("   epsout =",i5,10x,
     &            "  # (D=0) generate eps file by ANGEL")')
     &            iteps(m)

               end if

*-----------------------------------------------------------------------

            write(iot,'("#    used :",
     &                  "        main (  %)",
     &                  "        temp (  %)",
     &                  "       total (  %)"/
     &                  "#  memory :",
     &                  3(i12," (",i3,")"))')
     &                  itsmn(m), nint(dble(itsmn(m))/dble(mdas)*100),
     &                  itstm(m), nint(dble(itstm(m))/dble(mdas)*100),
     &                  itsmn(m)+itstm(m),
     &                  nint(dble(itsmn(m)+itstm(m))/dble(mdas)*100)

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine wtpname(m,np,chp,chq)
*                                                                      *
*       write particle name for tally                                  *
*                                                                      *
************************************************************************
      use partmod, only: itmxpt, itpan, itpat, jtpat, group ! frtati 2021/10/05
*-----------------------------------------------------------------------

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'

*-----------------------------------------------------------------------


*-----------------------------------------------------------------------

      common /ptname/ pname(20), ipln(20)
      character       pname*8

      character qname*16

cfrtati 2021/10/05 6 -> itmxpt
      character chp(itmxpt)*11  ! kitamura22/03/31
      character chq(itmxpt)*9   ! kitamura22/03/31

      character chau*8
      character elmnt(104)*3

cfrtati 2021/10/05 moved to partmod


*-----------------------------------------------------------------------

      data elmnt /
     &    'H  ','He ','Li ','Be ','B  ','C  ','N  ','O  ','F  ','Ne ',
     &    'Na ','Mg ','Al ','Si ','P  ','S  ','Cl ','Ar ','K  ','Ca ',
     &    'Sc ','Ti ','V  ','Cr ','Mn ','Fe ','Co ','Ni ','Cu ','Zn ',
     &    'Ga ','Ge ','As ','Se ','Br ','Kr ','Rb ','Sr ','Y  ','Zr ',
     &    'Nb ','Mo ','Tc ','Ru ','Rh ','Pd ','Ag ','Cd ','In ','Sn ',
     &    'Sb ','Te ','I  ','Xe ','Cs ','Ba ','La ','Ce ','Pr ','Nd ',
     &    'Pm ','Sm ','Eu ','Gd ','Tb ','Dy ','Ho ','Er ','Tm ','Yb ',
     &    'Lu ','Hf ','Ta ','W  ','Re ','Os ','Ir ','Pt ','Au ','Hg ',
     &    'Tl ','Pb ','Bi ','Po ','At ','Rn ','Fr ','Ra ','Ac ','Th ',
     &    'Pa ','U  ','Np ','Pu ','Am ','Cm ','Bk ','Cf ','Es ','Fm ',
     &    'Md ','No ','Lr ','Ku ' /

*-----------------------------------------------------------------------
*        particle name
*-----------------------------------------------------------------------

               do i = 1, np

                        istyp = itpat(m,i,1)
                        inkf0 = mod(itpat(m,i,2),1000000000)
                        isubt = itpat(m,i,3)                      ! kitamura22/03/31

                  if( istyp .lt. 0 ) then

                     mk = -istyp                                  ! kitamura22/03/31
                     if( isubt .eq. 0 .or. isubt .eq. 1 .and.     ! kitamura22/03/31
     &                   mk .gt. 1 ) then                         ! kitamura22/03/31

                        chp(i) = 'y('//group(i)(1:8)
                        chq(i) = group(i)(1:8)

                     else if( isubt .eq. 0 .and. mk .eq. 1 ) then ! kitamura22/03/31

                        istyp  = jtpat(m,i,mk,1)                  ! kitamura22/03/31
                        chp(i) = 'y('//pname(istyp)(1:8)          ! kitamura22/03/31
                        chq(i) = pname(istyp)(1:8)                ! kitamura22/03/31

                     else if( isubt .eq. 1 .and. mk .eq. 1 ) then ! kitamura22/03/31

                        istyp  = jtpat(m,i,mk,1)                  ! kitamura22/03/31
                        chp(i) = 'y(-'//pname(istyp)(1:8)         ! kitamura22/03/31
                        chq(i) = '-'//pname(istyp)(1:8)           ! kitamura22/03/31

                     end if                                       ! kitamura22/03/31

                  else if( istyp .eq. 19 .and. inkf0 .ne. 0 ) then

                     if( inkf0 .gt. 0 ) then

                        iz = ichgf(istyp,inkf0)
                        ia = ibryf(istyp,inkf0)

                        call chname(idum,ia,iz,chau)

                        ism = itpat(m,i,2)/1000000000 ! isomer
                        if(ism .gt. 0) then
                            do ii = 1, 8
                                if(chau(ii:ii) .eq. ' ') exit
                            enddo
                            if(ism .eq. 1) chau(ii:ii+1) = '-m'
                            if(ism .eq. 2) chau(ii:ii+1) = '-n'
                        endif

                        if( isubt .eq. 0 ) then                   ! kitamura22/03/31

                           chp(i) = 'y('//chau(1:8)
                           chq(i) = chau(1:8)

                        else                                      ! kitamura22/03/31

                           chp(i) = 'y(-'//chau(1:8)              ! kitamura22/03/31
                           chq(i) = '-'//chau(1:8)                ! kitamura22/03/31

                        end if                                    ! kitamura22/03/31

                     else

                        iz = ichgf(istyp,inkf0)
                        chau = elmnt(iz)//'     '

                        if( isubt .eq. 0 ) then                   ! kitamura22/03/31

                          chp(i) = 'y('//chau(1:8)
                          chq(i) = chau(1:8)

                        else                                      ! kitamura22/03/31

                           chp(i) = 'y(-'//chau(1:8)              ! kitamura22/03/31
                           chq(i) = '-'//chau(1:8)                ! kitamura22/03/31

                        end if                                    ! kitamura22/03/31

                     end if

                  else if( istyp .ne. 11 ) then

                     if( isubt .eq. 0 ) then                      ! kitamura22/03/31

                        chp(i) = 'y('//pname(istyp)(1:8)
                        chq(i) = pname(istyp)(1:8)

                     else                                         ! kitamura22/03/31

                        chp(i) = 'y(-'//pname(istyp)(1:8)         ! kitamura22/03/31
                        chq(i) = '-'//pname(istyp)(1:8)           ! kitamura22/03/31

                     end if                                       ! kitamura22/03/31

                  else

                        call jamname(inkf0,0,0,qname)
                        chp(i) = 'y('//qname(1:8)
                        chq(i) = qname(1:8)

                  end if

               end do

      return
      end

************************************************************************
*                                                                      *
      subroutine wtpname2(m,i,chqfull)
*                                                                      *
*       write particle name for tally without using p*-group           *
*                                                                      *
************************************************************************
      use partmod, only: itpat,jtpat
*-----------------------------------------------------------------------
      parameter(maxpart=18) ! maximum particle number specified in one group
      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'

      character chatmp180*180,chqfull*200
      dimension idimtmp(maxpart) ! temporary used dimension
*-----------------------------------------------------------------------

      if(-itpat(m,i,1).gt.maxpart) then
       chqfull='more than 18 part is specified, '//
     &         'please specify part by yourself'
      else
       do mk = 1, -itpat(m,i,1)
        idimtmp(mk) = jtpat(m,i,mk,2)
       enddo
       write(chatmp180,'(18i10)') (idimtmp(mk),mk=1,-itpat(m,i,1))
       chqfull='('//trim(chatmp180)//' )'
      endif

      return
      end



************************************************************************
cfrtati 2021/10/05 added argument imxpt, imxgp
*                                                                      *
      subroutine echprt(icc,iot,m,itpan,itpat,jtpat,itldd,imxpt,imxgp)
*                                                                      *
*       input echo for multiplier particle                                        *
*       last modified by K.Niita on 2005/11/30                         *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      include 'param.inc'

*-----------------------------------------------------------------------

      common /ptname/ pname(20), ipln(20)
      character       pname*8

cfrtati 2021/10/05 6 -> imxpt
      dimension itpan(itldd), itpat(itldd,imxpt,3),  ! kitamura22/03/31
     &          jtpat(itldd,imxpt,imxgp,2)           ! kitamura22/03/31

      character kname*9
      character tname*9
      character klin*2048
      character tlin*2048
      character qname*16

*-----------------------------------------------------------------------

      character chau*8
      character elmnt(104)*3

*-----------------------------------------------------------------------

      data elmnt /
     &    'H  ','He ','Li ','Be ','B  ','C  ','N  ','O  ','F  ','Ne ',
     &    'Na ','Mg ','Al ','Si ','P  ','S  ','Cl ','Ar ','K  ','Ca ',
     &    'Sc ','Ti ','V  ','Cr ','Mn ','Fe ','Co ','Ni ','Cu ','Zn ',
     &    'Ga ','Ge ','As ','Se ','Br ','Kr ','Rb ','Sr ','Y  ','Zr ',
     &    'Nb ','Mo ','Tc ','Ru ','Rh ','Pd ','Ag ','Cd ','In ','Sn ',
     &    'Sb ','Te ','I  ','Xe ','Cs ','Ba ','La ','Ce ','Pr ','Nd ',
     &    'Pm ','Sm ','Eu ','Gd ','Tb ','Dy ','Ho ','Er ','Tm ','Yb ',
     &    'Lu ','Hf ','Ta ','W  ','Re ','Os ','Ir ','Pt ','Au ','Hg ',
     &    'Tl ','Pb ','Bi ','Po ','At ','Rn ','Fr ','Ra ','Ac ','Th ',
     &    'Pa ','U  ','Np ','Pu ','Am ','Cm ','Bk ','Cf ','Es ','Fm ',
     &    'Md ','No ','Lr ','Ku ' /

*-----------------------------------------------------------------------

                        inm = 0
                        ic  = 0

                        tlin(ic+1:ic+11) = '     part ='
                        klin(ic+1:ic+11) = '# kf/name :'
                        ic = ic + 11

            do i = 1, itpan(m)

                        kname = '         '
                        tname = '         '

                        istyp = itpat(m,i,1)
                        inkf0 = mod(itpat(m,i,2),1000000000) ! 2025/3/10 Ogawa 1000000000 is isomer index
                        isubt = itpat(m,i,3)  ! kitamura22/03/31

                        ik = 1
                        kk = 0
                        is = 0                   ! kitamura22/03/31
                        ks = 0                   ! kitamura22/03/31

                     if( isubt .ne. 0 ) then     ! kitamura22/03/31
                        is = isubt               ! kitamura22/03/31
                     end if                      ! kitamura22/03/31

                     if( istyp .lt. 0 ) then

                        ik = -istyp
                        kk =  1
                        ks = -istyp              ! kitamura22/03/31

                     end if

                     if( is .ne. 0 ) then        ! kitamura22/03/31
                        tlin(ic+1:ic+3) = '  -'  ! kitamura22/03/31
                        klin(ic+1:ic+3) = '  -'  ! kitamura22/03/31
                        ic = ic + 3              ! kitamura22/03/31
                     end if                      ! kitamura22/03/31

               do k = 1, ik

                        ict = 0                  ! kitamura22/03/31
                        inm = inm + 1

                     if( kk .eq. 1 ) then

                        istyp = jtpat(m,i,k,1)
                        inkf0 = jtpat(m,i,k,2)

                     end if

                  if( istyp .eq. 19 .and. inkf0 .ne. 0 ) then

                        iz = ichgf(istyp,abs(inkf0))
                        ia = ibryf(istyp,abs(inkf0))

                        call chname(idum,ia,iz,chau)

                        ism = itpat(m,i,2)/1000000000 ! isomer
                        if(ism .gt. 0) then
                            do ii = 1, 8
                                if(chau(ii:ii) .eq. ' ') exit
                            enddo
                            if(ism .eq. 1) chau(ii:ii+1) = '-m'
                            if(ism .eq. 2) chau(ii:ii+1) = '-n'
                        endif

                        tname = chau
                        call kfcname(abs(inkf0),9,kname)

                  else if( istyp .ne. 11 ) then

                        tname = pname(istyp)(1:8)
                        call kfcname(inkf0,9,kname)

                  else

                        call jamname(inkf0,0,0,qname)
                        kname(1:8) = qname(1:8)

                        call kfcname(inkf0,9,tname)

                  end if

                  if( kk .eq. 1 .and. k .eq. 1 ) then

                     if( is .eq. 0 .and. ks .gt. 1 ) then       ! kitamura22/03/31
                       tlin(ic+1:ic+3) = ' ( '
                        klin(ic+1:ic+3) = ' ( '
                        ic = ic + 3
                     else if( is .eq. 1 .and. ks .gt. 1 ) then  ! kitamura22/03/31
                        tlin(ic+1:ic+2) = '( '                  ! kitamura22/03/31
                        klin(ic+1:ic+2) = '( '                  ! kitamura22/03/31
                        ic = ic + 2                             ! kitamura22/03/31
                     else if( is .eq. 1 ) then                  ! kitamura22/03/31
                        tlin(ic+1:ic+2) = '  '                  ! kitamura22/03/31
                        klin(ic+1:ic+2) = '( '                  ! kitamura22/03/31
                        ic = ic + 2                             ! kitamura22/03/31
                        ict = 1                                 ! kitamura22/03/31
                     end if                                     ! kitamura22/03/31

                  else

                     if( is .eq. 0 ) then                       ! kitamura22/03/31
                        tlin(ic+1:ic+2) = ' '
                        klin(ic+1:ic+2) = ' '
                        ic = ic + 2
                     end if                                     ! kitamura22/03/31

                  end if

                     if( is .eq. 1 .and. ks .eq. 0) then        ! kitamura22/03/31
                        tlin(ic+1:ic+2) = '  '                  ! kitamura22/03/31
                        klin(ic+1:ic+2) = '( '                  ! kitamura22/03/31
                        ic = ic + 2                             ! kitamura22/03/31
                        ict = 1                                 ! kitamura22/03/31
                     end if                                     ! kitamura22/03/31

                     if( ict .eq. 0 ) then                      ! kitamura22/03/31
                        tlin(ic+1:ic+9) = tname(1:9)
                        klin(ic+1:ic+9) = kname(1:9)
                        ic = ic + 9
                     else                                       ! kitamura22/03/31
                       ict = ic - 2                             ! kitamura22/03/31
                       tlin(ict+1:ict+11) = tname(1:9)          ! kitamura22/03/31
                       klin(ic+1:ic+9) = kname(1:9)             ! kitamura22/03/31
                       ic = ic + 9                              ! kitamura22/03/31
                     end if                                     ! kitamura22/03/31

                     if( is .eq. 1 .and. ks .eq. 0) then        ! kitamura22/03/31
                        tlin(ic+1:ic+2) = ' '                   ! kitamura22/03/31
                        klin(ic+1:ic+2) = ')'                   ! kitamura22/03/31
                        ic = ic + 2                             ! kitamura22/03/31
                     end if                                     ! kitamura22/03/31

                     if( kk .eq. 1 .and. k .eq. ik ) then

                        if( ks .gt. 1 ) then                    ! kitamura22/03/31
                         tlin(ic+1:ic+2) = ')'
                         klin(ic+1:ic+2) = ')'
                         ic = ic + 2
                        else                                    ! kitamura22/03/31
                         tlin(ic+1:ic+2) = ' '                  ! kitamura22/03/31
                         klin(ic+1:ic+2) = ')'                  ! kitamura22/03/31
                         ic = ic + 2                            ! kitamura22/03/31
                        end if                                  ! kitamura22/03/31
                     end if

               end do

            end do

                  write(iot,'(500a1)') ( tlin(i:i), i = 1, ic )

               if( icc .eq. 1 ) then

                  write(iot,'(500a1)') ( klin(i:i), i = 1, ic )

               end if

*-----------------------------------------------------------------------

      return
      end


*****************************************************************************************
*                                                                      *
      subroutine echcprt(icc,iot,m,icpan,icpat)
*                                                                      *
*       input echo for particle of counter                             *
*       last modified by K.Niita on 2006/01/16                         *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      include 'param.inc'

*-----------------------------------------------------------------------

      common /ptname/ pname(20), ipln(20)
      character       pname*8

      dimension icpan(3), icpat(3,20,2)

      character kname*9
      character tname*9
      character klin*2048
      character tlin*2048
      character qname*16

*-----------------------------------------------------------------------

      character chau*8
      character elmnt(104)*3

*-----------------------------------------------------------------------

      data elmnt /
     &    'H  ','He ','Li ','Be ','B  ','C  ','N  ','O  ','F  ','Ne ',
     &    'Na ','Mg ','Al ','Si ','P  ','S  ','Cl ','Ar ','K  ','Ca ',
     &    'Sc ','Ti ','V  ','Cr ','Mn ','Fe ','Co ','Ni ','Cu ','Zn ',
     &    'Ga ','Ge ','As ','Se ','Br ','Kr ','Rb ','Sr ','Y  ','Zr ',
     &    'Nb ','Mo ','Tc ','Ru ','Rh ','Pd ','Ag ','Cd ','In ','Sn ',
     &    'Sb ','Te ','I  ','Xe ','Cs ','Ba ','La ','Ce ','Pr ','Nd ',
     &    'Pm ','Sm ','Eu ','Gd ','Tb ','Dy ','Ho ','Er ','Tm ','Yb ',
     &    'Lu ','Hf ','Ta ','W  ','Re ','Os ','Ir ','Pt ','Au ','Hg ',
     &    'Tl ','Pb ','Bi ','Po ','At ','Rn ','Fr ','Ra ','Ac ','Th ',
     &    'Pa ','U  ','Np ','Pu ','Am ','Cm ','Bk ','Cf ','Es ','Fm ',
     &    'Md ','No ','Lr ','Ku ' /

*-----------------------------------------------------------------------

                        ic  = 0

                     if( icpan(m) .gt. 0 ) then

                        tlin(ic+1:ic+11) = '     part ='

                     else

                        tlin(ic+1:ic+11) = '    *part ='

                     end if

                        klin(ic+1:ic+11) = '# kf/name :'
                        ic = ic + 11

            do i = 1, abs( icpan(m) )

                        kname = '         '
                        tname = '         '

                        istyp = icpat(m,i,1)
                        inkf0 = icpat(m,i,2)

                  if( istyp .eq. 19 .and. inkf0 .ne. 0 ) then

                        iz = ichgf(istyp,abs(inkf0))
                        ia = ibryf(istyp,abs(inkf0))

                        call chname(idum,ia,iz,chau)

                        tname = chau
                        call kfcname(abs(inkf0),9,kname)

                  else if( istyp .ne. 11 ) then

                        tname = pname(istyp)(1:8)
                        call kfcname(inkf0,9,kname)

                  else

                        call jamname(inkf0,0,0,qname)
                        kname(1:8) = qname(1:8)

                        call kfcname(inkf0,9,tname)

                  end if

                        tlin(ic+1:ic+2) = ' '
                        klin(ic+1:ic+2) = ' '
                        ic = ic + 2
                        tlin(ic+1:ic+9) = tname(1:9)
                        klin(ic+1:ic+9) = kname(1:9)
                        ic = ic + 9

            end do

                  write(iot,'(500a1)') ( tlin(i:i), i = 1, ic )

               if( icc .eq. 1 ) then

                  write(iot,'(500a1)') ( klin(i:i), i = 1, ic )

               end if

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine echmty(iccs,iot,dm,ityp,itnm,del,rmin,rmax,n
     &                  ,ndim_mesh,das_mesh)
*                                                                      *
*       input echo for mesh                                            *
*       modified by K.Niita on 2000/07/12                              *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      include 'param.inc'

*-----------------------------------------------------------------------

      character dm*1
      character ds*10
      double precision, intent(in) :: das_mesh(n+ndim_mesh)

*-----------------------------------------------------------------------

                  ds = '#         '

                  jtyp = abs( ityp )

            if( jtyp .eq. 1 ) then

                  write(iot,'("   ",a1,"-type =   ",i2,11x,
     &            " # ",a1,"-mesh is given by the below data")')
     &            dm, ityp, dm

            else if( jtyp .eq. 2 ) then

                  write(iot,'("   ",a1,"-type =   ",i2,11x,
     &            " # ",a1,"-mesh is linear given by ",a1,
     &            "min, ",a1,"max and n",a1)')
     &            dm, ityp, dm, dm, dm, dm

            else if( jtyp .eq. 3 ) then

                  write(iot,'("   ",a1,"-type =   ",i2,11x,
     &            " # ",a1,"-mesh is log given by ",a1,
     &            "min, ",a1,"max and n",a1)')
     &            dm, ityp, dm, dm, dm, dm

            else if( jtyp .eq. 4 ) then

                  write(iot,'("   ",a1,"-type =   ",i2,11x,
     &            " # ",a1,"-mesh is linear given by ",a1,
     &            "min, ",a1,"max and ",a1,"del")')
     &            dm, ityp, dm, dm, dm, dm

            else if( jtyp .eq. 5 ) then

                  write(iot,'("   ",a1,"-type =   ",i2,11x,
     &            " # ",a1,"-mesh is log given by ",a1,
     &            "min, ",a1,"max and ",a1,"del")')
     &            dm, ityp, dm, dm, dm, dm

            end if

         if( jtyp .eq. 1 ) then

                  write(iot,'("       n",a1," =",i7,9x,
     &            " # number of ",a1,"-mesh points")') dm, itnm, dm

                  ds = '          '

         else if( jtyp .eq. 2 .or. jtyp .eq. 3 ) then

                  write(iot,'("     ",a1,"min = ",1p1g14.7,1x,
     &            " # minimum value of ",a1,"-mesh points")')
     &            dm, rmin, dm

                  write(iot,'("     ",a1,"max = ",1p1g14.7,1x,
     &            " # maximum value of ",a1,"-mesh points")')
     &            dm, rmax, dm

                  write(iot,'("#    ",a1,"del = ",1p1g14.7,1x,
     &            " # mesh width of ",a1,"-mesh points")')
     &            dm, del, dm

                  write(iot,'("       n",a1," =",i7,9x,
     &            " # number of ",a1,"-mesh points")') dm, itnm, dm

         else if( jtyp .eq. 4 .or. jtyp .eq. 5 ) then

                  write(iot,'("     ",a1,"min = ",1p1g14.7,1x,
     &            " # minimum value of ",a1,"-mesh points")')
     &            dm, rmin, dm

                  write(iot,'("     ",a1,"max = ",1p1g14.7,1x,
     &            " # maximum value of ",a1,"-mesh points")')
     &            dm, rmax, dm

                  write(iot,'("     ",a1,"del = ",1p1g14.7,1x,
     &            " # mesh width of ",a1,"-mesh points")')
     &            dm, del, dm

                  write(iot,'("#      n",a1," =",i5,11x,
     &            " # number of ",a1,"-mesh points")') dm, itnm, dm

         end if

*-----------------------------------------------------------------------

         if( jtyp .eq. 1 .or. iccs .eq. 0 ) then

                  write(iot,'("#    data = ( ",a1,"(i), i = 1,",
     &            " n",a1," + 1 )")') dm, dm

                  ism = ( iabs(itnm) + 1 ) / 5
                  isa = ( iabs(itnm) + 1 ) - ( iabs(itnm) + 1 ) / 5 * 5

            if( ism .gt. 0 ) then

               do i = 1, ism

                  k = ( i - 1 ) * 5 + 1

                  write(iot,'(a10,1p5e13.5)')
     &                  ds, (das_mesh(n+j-1),j=k,k+4)

               end do

            end if

            if( isa .gt. 0 ) then

                  i = ( iabs(itnm) + 1 ) / 5 * 5 + 1

                  write(iot,'(a10,1p5e13.5)')
     &                  ds, (das_mesh(n+j-1),j=i,iabs(itnm)+1)

            end if

         end if

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine echmty2(iccs,iot,dm,dm2,ityp,itnm,del,rmin,rmax,n
     &                   ,ndim_mesh,das_mesh)
*                                                                      *
*       input echo for mesh with 2 characters                          *
*       modified by K.Niita on 2005/11/14                              *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      include 'param.inc'

*-----------------------------------------------------------------------

      character dm*1
      character dm2*2
      character ds*10
      double precision, intent(in) :: das_mesh(n+ndim_mesh)

*-----------------------------------------------------------------------

                  ds = '#         '

                  jtyp = abs( ityp )

            if( jtyp .eq. 1 ) then

                  write(iot,'("  ",a2,"-type =   ",i2,11x,
     &            " # ",a2,"-mesh is given by the below data")')
     &            dm2, ityp, dm2

            else if( jtyp .eq. 2 ) then

                  write(iot,'("  ",a2,"-type =   ",i2,11x,
     &            " # ",a2,"-mesh is linear given by ",a1,
     &            "min, ",a1,"max and n",a1)')
     &            dm2, ityp, dm2, dm, dm, dm

            else if( jtyp .eq. 3 ) then

                  write(iot,'("  ",a2,"-type =   ",i2,11x,
     &            " # ",a2,"-mesh is log given by ",a1,
     &            "min, ",a1,"max and n",a1)')
     &            dm2, ityp, dm2, dm, dm, dm

            else if( jtyp .eq. 4 ) then

                  write(iot,'("  ",a2,"-type =   ",i2,11x,
     &            " # ",a2,"-mesh is linear given by ",a1,
     &            "min, ",a1,"max and ",a1,"del")')
     &            dm2, ityp, dm2, dm, dm, dm

            else if( jtyp .eq. 5 ) then

                  write(iot,'("  ",a2,"-type =   ",i2,11x,
     &            " # ",a2,"-mesh is log given by ",a1,
     &            "min, ",a1,"max and ",a1,"del")')
     &            dm2, ityp, dm2, dm, dm, dm

            end if

         if( jtyp .eq. 1 ) then

                  write(iot,'("       n",a1," =",i7,9x,
     &            " # number of ",a1,"-mesh points")') dm, itnm, dm

                  ds = '          '

         else if( jtyp .eq. 2 .or. jtyp .eq. 3 ) then

                  write(iot,'("     ",a1,"min = ",1p1g14.7,1x,
     &            " # minimum value of ",a1,"-mesh points")')
     &            dm, rmin, dm

                  write(iot,'("     ",a1,"max = ",1p1g14.7,1x,
     &            " # maximum value of ",a1,"-mesh points")')
     &            dm, rmax, dm

                  write(iot,'("#    ",a1,"del = ",1p1g14.7,1x,
     &            " # mesh width of ",a1,"-mesh points")')
     &            dm, del, dm

                  write(iot,'("       n",a1," =",i5,11x,
     &            " # number of ",a1,"-mesh points")') dm, itnm, dm

         else if( jtyp .eq. 4 .or. jtyp .eq. 5 ) then

                  write(iot,'("     ",a1,"min = ",1p1g14.7,1x,
     &            " # minimum value of ",a1,"-mesh points")')
     &            dm, rmin, dm

                  write(iot,'("     ",a1,"max = ",1p1g14.7,1x,
     &            " # maximum value of ",a1,"-mesh points")')
     &            dm, rmax, dm

                  write(iot,'("     ",a1,"del = ",1p1g14.7,1x,
     &            " # mesh width of ",a1,"-mesh points")')
     &            dm, del, dm

                  write(iot,'("#      n",a1," =",i5,11x,
     &            " # number of ",a1,"-mesh points")') dm, itnm, dm

         end if

*-----------------------------------------------------------------------

         if( jtyp .eq. 1 .or. iccs .eq. 0 ) then

                  write(iot,'("#    data = ( ",a1,"(i), i = 1,",
     &            " n",a1," + 1 )")') dm, dm

                  ism = ( iabs(itnm) + 1 ) / 5
                  isa = ( iabs(itnm) + 1 ) - ( iabs(itnm) + 1 ) / 5 * 5

            if( ism .gt. 0 ) then

               do i = 1, ism

                  k = ( i - 1 ) * 5 + 1

                  write(iot,'(a10,1p5e13.5)')
     &                  ds, (das_mesh(n+j-1),j=k,k+4)

               end do

            end if

            if( isa .gt. 0 ) then

                  i = ( iabs(itnm) + 1 ) / 5 * 5 + 1

                  write(iot,'(a10,1p5e13.5)')
     &                  ds, (das_mesh(n+j-1),j=i,iabs(itnm)+1)

            end if

         end if

*-----------------------------------------------------------------------

      return
      end


************************************************************************

*                                                                      *
      subroutine echrg_wwvm(icc,ierl,iva,ivs,iot,mc,kc,nr,mr,kr,
     &                      vl,lr,nvl,ivl,rvl,igm)
*                                                                      *

*       input echo for region mesh of t-volume and t-wwbg              *
*       modified by K.Niita on 2017/02/17                              *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      include 'param.inc'
      include 'ggsparam.inc'

      common /regdc/ idrg(kvlmax), idgr(kvmmax)

*-----------------------------------------------------------------------

      character dum1*10000
      character dum2*10000
      character chdf(1000)*100
      dimension ildf(1000)
      character cblan*100

      dimension kc(mc)
      dimension kr(mr)

      dimension vl(nr)
      dimension lr(nr)
      dimension ivl(nvl)
      dimension rvl(nvl)

*-----------------------------------------------------------------------
*     read from kc(mc)
*-----------------------------------------------------------------------

            call echrg2(mc,kc,dum1,lng1,icmb,igm)

*-----------------------------------------------------------------------
*     input echo for tally region mesh
*-----------------------------------------------------------------------

                  cblan = ' '

                  ild1 = 12
                  ild0 = 72 - ild1

                  ilrm = lng1
                  isqd = 0
                  isrm = 0

  430             isqd = isqd + 1

                  if( ilrm .le. ild0 ) then

                        ildf(isqd) = ilrm

                     do k = 1, ilrm

                        chdf(isqd)(k:k) = dum1(k+isrm:k+isrm)

                     end do

                  else

                     do k = ild0, 1, -1

                        if( dum1(k+isrm:k+isrm) .eq. ' ' ) goto 420

                     end do

  420                k2 = k

                        ildf(isqd) = k2 - 1

                     do k = 1, k2 - 1

                        chdf(isqd)(k:k) = dum1(k+isrm:k+isrm)

                     end do

                     isrm = isrm + k2
                     ilrm = ilrm - k2

                     goto 430

                  end if

                  write(iot,'("      reg = ",100a1)')
     &                    (chdf(1)(j:j),j=1,ildf(1))

               if( isqd .gt. 1 ) then

                  do k = 2, isqd

                     write(iot,'(100a1)') (cblan(j:j),j=1,ild1),
     &                    (chdf(k)(j:j),j=1,ildf(k))

                  end do

               end if

*-----------------------------------------------------------------------
*        echo for combined regions
*-----------------------------------------------------------------------

         if( junf .eq. 0 .and. icmb .eq. 0 ) return


            write(iot,'("#   non     non   reg")')

            call tregvl(nr,mr,kr,vl,lr,nvl,ivl,rvl)

*-----------------------------------------------------------------------
*     read from kr(mr)
*-----------------------------------------------------------------------

                  j = 0

         do 500 ir = 1, nr

               call echrg3(j,mr,kr,dum2,lng1,icmb,igm)

                  ild1 = 31
                  ild0 = max(10,ierl-ild1)

                  ilrm = lng1
                  isqd = 0
                  isrm = 0

  530             isqd = isqd + 1

                  if( ilrm .le. ild0 ) then

                        ildf(isqd) = ilrm

                     do k = 1, ilrm

                        chdf(isqd)(k:k) = dum2(k+isrm:k+isrm)

                     end do

                  else

                     do k = ild0, 1, -1

                        if( dum2(k+isrm:k+isrm) .eq. ' ' ) goto 520

                     end do

  520                k2 = k

                        ildf(isqd) = k2 - 1

                     do k = 1, k2 - 1

                        chdf(isqd)(k:k) = dum2(k+isrm:k+isrm)

                     end do

                     isrm = isrm + k2
                     ilrm = ilrm - k2

                     goto 530

                  end if


                  write(iot,'("#",i5,2x,i7,"   ",100a1)')
     &                 ir, lr(ir), (chdf(1)(l:l),l=1,ildf(1))

               if( isqd .gt. 1 ) then

                  do k = 2, isqd

                     write(iot,'(100a1)') (cblan(l:l),'#',l=2,ild1-2),
     &                    '  ',(chdf(k)(l:l),l=1,ildf(k))

                  end do

               end if

*-----------------------------------------------------------------------

  500    continue

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine echrg(icc,ierl,iva,ivs,iot,mc,kc,nr,mr,kr,
     &                 vl,lr,nvl,ivl,rvl,igm)
*                                                                      *
*       input echo for tally region mesh                               *
*       modified by K.Niita on 2011/02/03                              *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

      include 'param.inc'
      include 'ggsparam.inc'

      common /regdc/ idrg(kvlmax), idgr(kvmmax)

*-----------------------------------------------------------------------

      character dum1*10000
      character dum2*10000
      character chdf(1000)*100
      dimension ildf(1000)
      character cblan*100

      dimension kc(mc)
      dimension kr(mr)

      dimension vl(nr)
      dimension lr(nr)
      dimension ivl(nvl)
      dimension rvl(nvl)

*-----------------------------------------------------------------------
*     read from kc(mc)
*-----------------------------------------------------------------------

            call echrg2(mc,kc,dum1,lng1,icmb,igm)

*-----------------------------------------------------------------------
*     input echo for tally region mesh
*-----------------------------------------------------------------------

                  cblan = ' '

                  ild1 = 12
                  ild0 = 72 - ild1

                  ilrm = lng1
                  isqd = 0
                  isrm = 0

  430             isqd = isqd + 1

                  if( ilrm .le. ild0 ) then

                        ildf(isqd) = ilrm

                     do k = 1, ilrm

                        chdf(isqd)(k:k) = dum1(k+isrm:k+isrm)

                     end do

                  else

                     do k = ild0, 1, -1

                        if( dum1(k+isrm:k+isrm) .eq. ' ' ) goto 420

                     end do

  420                k2 = k

                        ildf(isqd) = k2 - 1

                     do k = 1, k2 - 1

                        chdf(isqd)(k:k) = dum1(k+isrm:k+isrm)

                     end do

                     isrm = isrm + k2
                     ilrm = ilrm - k2

                     goto 430

                  end if

               if( icc .eq. 1 ) then

                  write(iot,'("      reg = ",100a1)')
     &                    (chdf(1)(j:j),j=1,ildf(1))

               else if( icc .eq. 2 ) then

                  write(iot,'(" reginbox = ",100a1)')
     &                    (chdf(1)(j:j),j=1,ildf(1))

               end if

               if( isqd .gt. 1 ) then

                  do k = 2, isqd

                     write(iot,'(100a1)') (cblan(j:j),j=1,ild1),
     &                    (chdf(k)(j:j),j=1,ildf(k))

                  end do

               end if

*-----------------------------------------------------------------------
*        echo for combined regions
*-----------------------------------------------------------------------

         if( ivs .eq. 0 .or. ( iva .eq. 0 .and.
     &       junf .eq. 0 .and. icmb .eq. 0 .and. nvl .eq. 0 ) ) return

         if( iva .eq. 0 ) then

            write(iot,'("   volume  ",17x,
     &                  "# combined, lattice or level structure ")')
            write(iot,'("   non     reg      vol     #",
     &                  " reg definition")')

         else

            write(iot,'("   value  ",18x,
     &                  "# values for each region")')
            write(iot,'("   non     reg      val     #",
     &                  " reg definition")')

         end if

               call tregvl(nr,mr,kr,vl,lr,nvl,ivl,rvl)

*-----------------------------------------------------------------------
*     read from kr(mr)
*-----------------------------------------------------------------------

                  j = 0

         do 500 ir = 1, nr

               call echrg3(j,mr,kr,dum2,lng1,icmb,igm)

                  ild1 = 30
                  ild0 = max(10,ierl-ild1)

                  ilrm = lng1
                  isqd = 0
                  isrm = 0

  530             isqd = isqd + 1

                  if( ilrm .le. ild0 ) then

                        ildf(isqd) = ilrm

                     do k = 1, ilrm

                        chdf(isqd)(k:k) = dum2(k+isrm:k+isrm)

                     end do

                  else

                     do k = ild0, 1, -1

                        if( dum2(k+isrm:k+isrm) .eq. ' ' ) goto 520

                     end do

  520                k2 = k

                        ildf(isqd) = k2 - 1

                     do k = 1, k2 - 1

                        chdf(isqd)(k:k) = dum2(k+isrm:k+isrm)

                     end do

                     isrm = isrm + k2
                     ilrm = ilrm - k2

                     goto 530

                  end if

                  write(iot,'(i5,2x,i7,1p1e13.4," # ",100a1)')
     &                 ir, lr(ir), vl(ir), (chdf(1)(l:l),l=1,ildf(1))

               if( isqd .gt. 1 ) then

                  do k = 2, isqd

                     write(iot,'(100a1)') (cblan(l:l),l=1,ild1-2),
     &                    '#',' ',(chdf(k)(l:l),l=1,ildf(k))

                  end do

               end if

*-----------------------------------------------------------------------

  500    continue

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine echrg2(mc,kc,dum1,lng1,icmb,igm_argument)
*                                                                      *
*       write region decription from kc(mc)                            *
*       modified by K.Niita on 2001/11/28                              *
*                                                                      *
************************************************************************
      use moddas
      use moddas_region

      implicit real*8 (a-h,o-z)

      include 'param.inc'

*-----------------------------------------------------------------------

      dimension kc(mc)
      dimension lats(6)

      character dum1*10000
      character dum2*100000

      include 'err.inc'
      integer, intent(in) :: igm_argument
      igm = 1
      call moddas_allocate_int(MAX_NUM_NRST, nrst)

*-----------------------------------------------------------------------

               j = 0
               l = 1
               icmb = 0
               dum2 = ' '
               nnm = 0

  100       j = j + 1

               if( j .gt. mc ) goto 200

               nreg = kc(j)

            if( nreg .eq. 5000000 ) then

                     dum2(l:l+4) = ' all '
                     l = l + 5

            else if( nreg .eq. 6000000 ) then

                     dum2(l:l+8) = ' ( all ) '
                     l = l + 9

                     icmb = icmb + 1

            else if( nreg .eq. 1000000 .or. nreg .eq. 2000000 ) then

                     dum2(l:l+2) = ' ( '
                     l = l + 3

               if( nreg .eq. 1000000 ) then

                     j = j + 1

               else if( nreg .eq. 2000000 ) then

                     j = j + 1
                     jlev = kc(j)

                  do k = 1, jlev

                     j = j + 1

                  end do

               end if

                     icmb = icmb + 1

            else if( nreg .eq. 1000001 .or. nreg .eq. 2000001 ) then

                     dum2(l:l+2) = ' ) '
                     l = l + 3

            else if( nreg .eq. 3000000 ) then

                     dum2(l:l+2) = ' < '
                     l = l + 3

            else if( nreg .gt. 7000000 ) then

                     ireg = nreg - 7000000

                     dum2(l:l+3) = ' u ='
                     l = l + 43

                     write(dum2(l:l+6),'(i7)') ireg
                     l = l + 7

            else if( nreg .gt. 4000000 .and. nreg .lt. 5000000 ) then

                     icmb = icmb + 1

                     ireg = nreg - 4000000

                     write(dum2(l:l+6),'(i7)') ireg
                     l = l + 7

                     j = j + 1
                     ilat = kc(j)

                     dum2(l:l+1) = '[ '
                     l = l + 2

                     do ll = 1, ilat

                        do m = 1, 6

                           j = j + 1
                           lats(m) = kc(j)

                        end do

                        if( lats(1) .eq. lats(2) .and.
     &                      lats(3) .eq. lats(4) .and.
     &                      lats(5) .eq. lats(6) ) then

                           do m = 1, 5, 2

                              write(dum2(l:l+6),'(i7)') lats(m)
                              l = l + 7

                           end do

                        else

                           do m = 1, 6

                                 write(dum2(l:l+6),'(i7)') lats(m)
                                 l = l + 7

                              if( m .eq. 1 .or. m .eq. 3 .or.
     &                            m .eq. 5 ) then

                                 dum2(l:l) = ':'
                                 l = l + 1

                              end if

                           end do

                        end if

                        if( ll .ne. ilat ) then

                           dum2(l:l+1) = ', '
                           l = l + 2

                        end if

                     end do

                        dum2(l:l+2) = ' ] '
                        l = l + 3

            else if( nreg .gt. 0 .and. nreg .lt. 1000000 ) then

                           nnm = nnm + 1
                           nrst(nnm+igm) = nreg

                     if( nnm .eq. 2 .and.
     &                   nrst(nnm+igm) .ne. nrst(nnm-1+igm)+1 ) then

                           write(dum2(l:l+6),'(i7)') nrst(nnm-1+igm)
                           l = l + 7

                           nnm = 1
                           nrst(nnm+igm) = nreg

                     else if( nnm .gt. 2 .and.
     &                    nrst(nnm+igm) .ne. nrst(nnm-1+igm)+1 ) then

                        if( nnm .eq. 3 ) then

                           write(dum2(l:l+6),'(i7)') nrst(1+igm)
                           l = l + 7

                           write(dum2(l:l+6),'(i7)') nrst(2+igm)
                           l = l + 7

                        else

                           dum2(l:l+2) = ' { '
                           l = l + 3

                           write(dum2(l:l+6),'(i7)') nrst(1+igm)
                           l = l + 7

                           dum2(l:l+2) = ' - '
                           l = l + 3

                           write(dum2(l:l+6),'(i7)') nrst(nnm-1+igm)
                           l = l + 7

                           dum2(l:l+2) = ' } '
                           l = l + 3

                        end if

                           nnm = 1
                           nrst(nnm+igm) = nreg

                     end if

                  idoprocess = 0
                  if( j+1 .gt. mc ) then
                    idoprocess = 1
                  else if( kc(j+1) .ge. 1000000 .or.
     &                     kc(j+1) .lt. 0 ) then
                    idoprocess = 1
                  end if
                  if( idoprocess .eq. 1 ) then

                        if( nnm .eq. 1 ) then

                           write(dum2(l:l+6),'(i7)') nrst(1+igm)
                           l = l + 7

                        else if( nnm .eq. 2 ) then

                           write(dum2(l:l+6),'(i7)') nrst(1+igm)
                           l = l + 7

                           write(dum2(l:l+6),'(i7)') nrst(2+igm)
                           l = l + 7

                        else

                           dum2(l:l+2) = ' { '
                           l = l + 3

                           write(dum2(l:l+6),'(i7)') nrst(1+igm)
                           l = l + 7

                           dum2(l:l+2) = ' - '
                           l = l + 3

                           write(dum2(l:l+6),'(i7)') nrst(nnm+igm)
                           l = l + 7

                           dum2(l:l+2) = ' } '
                           l = l + 3

                        end if

                           nnm = 0

                  end if

            end if

            goto 100

  200    continue

               call chlngt(dum2,10000,i1,i2)

                  l = 0

               do k = i1, i2

                 if( dum2(k:k) .ne. ' ' .or.
     &               dum2(k+1:k+1) .ne. ' ' ) then
                   idoprocess = 0
                   if( l.eq.0 ) then
                     idoprocess = 1
                   else if ( dum1(l:l) .ne. ':' .or.
     &                       dum2(k:k) .ne. ' ' ) then
                     idoprocess = 1
                   end if
                   if( idoprocess .eq. 1 ) then

                     l = l + 1
                     dum1(l:l) = dum2(k:k)

                   end if
                 end if

               end do

                     lng1 = l

*-----------------------------------------------------------------------
      if( nnm > MAX_NUM_NRST ) then
         write(ErrCha,'(a,a,i5,a,a,i5,a)')
     &        'sub.echrg2@tallsm3.f'
     &           //' ?dimension over nrst?'
     &           //' nnm > MAX_NUM_NRST'
     &        ,' (nnm=',nnm,')'
     &        ,' (MAX_NUM_NRST@moddas.f=',MAX_NUM_NRST,')'
         ErrID = 'L:15637/R:echrg2/F:tallsm3.f'
         call ErrWrite(ErrID,ErrCha)
      endif
      call moddas_deallocate_int(nrst)

      return
      end

************************************************************************
*                                                                      *
      subroutine echrg3(j,mr,kr,dum2,lng1,icmb,igm_argument)
*                                                                      *
*       write region decription from kr(mr)                            *
*       modified by K.Niita on 2011/02/03                              *
*                                                                      *
************************************************************************
      use moddas
      use moddas_region

      implicit real*8 (a-h,o-z)

      include 'param.inc'

*-----------------------------------------------------------------------

      dimension kr(mr)
      dimension lats(6)

      dimension ipar(0:20)
      dimension jpar(0:20)
      dimension klev(0:20)

      character dum1*10000
      character dum2*10000

      include 'err.inc'
      integer, intent(in) :: igm_argument
      igm = 1
      call moddas_allocate_int(MAX_NUM_NRST, nrst)

*-----------------------------------------------------------------------

               j = j + 1
               nreg = kr(j)

               icmb = 0
               dum1 = ' '
               l    = 1
               nnm = 0

*-----------------------------------------------------------------------

            if( nreg .eq. 6000000 ) then

                     dum1(l:l+8) = ' ( all ) '
                     l = l + 9

                     icmb = icmb + 1

                     goto 560

            else if( nreg .gt. 0 .and. nreg .lt. 1000000 ) then

                     write(dum1(l:l+6),'(i7)') nreg
                     l = l + 7

                     goto 560

            end if

*-----------------------------------------------------------------------
*           combine, lattice or level structure entries
*-----------------------------------------------------------------------

                     kpar = 0

                  do i = 0, 20

                     ipar(i) = 0
                     jpar(i) = 0
                     klev(i) = 0

                  end do

                     icmb = icmb + 1

*-----------------------------------------------------------------------

  600       continue

               if( nreg .lt. 0 ) then

                        dum1(l:l+2) = ' ( '
                        l = l + 3

                        kpar = kpar + 1
                        jpar(kpar) = -nreg

               else if( nreg .gt. 3000000 .and. nreg .lt. 4000000 ) then

                        dum1(l:l+2) = ' ( '
                        l = l + 3

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

                        write(dum1(l:l+6),'(i7)') nreg
                        l = l + 7

                        j = j + 1
                        ilat = kr(j)

                        dum1(l:l+1) = '[ '
                        l = l + 2

                     do ll = 1, ilat

                        do m = 1, 6

                           j = j + 1
                           lats(m) = kr(j)

                        end do

                        if( lats(1) .eq. lats(2) .and.
     &                      lats(3) .eq. lats(4) .and.
     &                      lats(5) .eq. lats(6) ) then

                           do m = 1, 5, 2

                              write(dum1(l:l+6),'(i7)') lats(m)
                              l = l + 7

                           end do

                        else

                           do m = 1, 6

                                 write(dum1(l:l+6),'(i7)') lats(m)
                                 l = l + 7

                              if( m .eq. 1 .or. m .eq. 3 .or.
     &                            m .eq. 5 ) then

                                 dum1(l:l) = ':'
                                 l = l + 1

                              end if

                           end do

                        end if

                        if( ll .ne. ilat ) then

                           dum1(l:l+1) = ', '
                           l = l + 2

                        end if

                     end do

               else if( nreg .gt. 0 .and. nreg .lt. 1000000 ) then

                           ipar(kpar) = ipar(kpar) + 1

                           nnm = nnm + 1
                           nrst(nnm+igm) = nreg

                     if( nnm .eq. 2 .and.
     &                   nrst(nnm+igm) .ne. nrst(nnm-1+igm)+1 ) then

                           write(dum1(l:l+6),'(i7)') nrst(nnm-1+igm)
                           l = l + 7

                           nnm = 1
                           nrst(nnm+igm) = nreg

                     else if( nnm .gt. 2 .and.
     &                    nrst(nnm+igm) .ne. nrst(nnm-1+igm)+1 ) then

                        if( nnm .eq. 3 ) then

                           write(dum1(l:l+6),'(i7)') nrst(1+igm)
                           l = l + 7

                           write(dum1(l:l+6),'(i7)') nrst(2+igm)
                           l = l + 7

                        else

                           dum1(l:l+2) = ' { '
                           l = l + 3

                           write(dum1(l:l+6),'(i7)') nrst(1+igm)
                           l = l + 7

                           dum1(l:l+2) = ' - '
                           l = l + 3

                           write(dum1(l:l+6),'(i7)') nrst(nnm-1+igm)
                           l = l + 7

                           dum1(l:l+2) = ' } '
                           l = l + 3

                        end if

                           nnm = 1
                           nrst(nnm+igm) = nreg

                     end if

                  idoprocess = 0
                  if( ipar(kpar) .eq. jpar(kpar) .or.
     &              ( klev(kpar) .gt. 0 .and. ipar(kpar) .gt. 0 ) ) then
                    idoprocess = 1
                  else if( j+1 .le. mr )  then
                    if( kr(j+1) .ge. 1000000 .or.
     &                  kr(j+1) .lt. 0 ) then
                      idoprocess = 1
                    end if
                  end if
                  if( idoprocess .eq. 1 ) then

                        if( nnm .eq. 1 ) then

                           write(dum1(l:l+6),'(i7)') nrst(1+igm)
                           l = l + 7

                        else if( nnm .eq. 2 ) then

                           write(dum1(l:l+6),'(i7)') nrst(1+igm)
                           l = l + 7

                           write(dum1(l:l+6),'(i7)') nrst(2+igm)
                           l = l + 7

                        else

                           dum1(l:l+2) = ' { '
                           l = l + 3

                           write(dum1(l:l+5),'(i6)') nrst(1+igm)
                           l = l + 6

                           dum1(l:l+2) = ' - '
                           l = l + 3

                           write(dum1(l:l+5),'(i6)') nrst(nnm+igm)
                           l = l + 6

                           dum1(l:l+2) = ' } '
                           l = l + 3

                        end if

                           nnm = 0

                  end if

               end if

*-----------------------------------------------------------------------

  660          continue

               if( ipar(kpar) .eq. jpar(kpar) ) then

                  if( jpar(kpar) .gt. 0 ) then

                     dum1(l:l+2) = ' ) '
                     l = l + 3

                  else

                     dum1(l:l+2) = ' ] '
                     l = l + 3

                  end if

                     ipar(kpar) = 0
                     jpar(kpar) = 0

                  if( klev(kpar) .gt. 0 ) then

                     klev(kpar) = 0

                  end if

                     kpar = kpar - 1
                     ipar(kpar) = ipar(kpar) + 1

                  if( kpar .eq. 0 ) goto 560

                     goto 660

               end if


               if( klev(kpar) .gt. 0 .and. ipar(kpar) .gt. 0 ) then

                     dum1(l:l+2) = ' < '
                     l = l + 3

               end if

                     j = j + 1
                     nreg = kr(j)

                     goto 600

*-----------------------------------------------------------------------

  560    continue

               call chlngt(dum1,10000,i1,i2)

                  l = 0

               do k = i1, i2

                 if( dum1(k:k) .ne. ' ' .or.
     &               dum1(k+1:k+1) .ne. ' ' ) then
                   idoprocess = 0
                   if( l.eq.0 ) then
                     idoprocess = 1
                   else if ( dum2(l:l) .ne. ':' .or.
     &                       dum1(k:k) .ne. ' ' ) then
                     idoprocess = 1
                   end if
                   if( idoprocess .eq. 1 ) then

                     l = l + 1
                     dum2(l:l) = dum1(k:k)

                   end if
                 end if

               end do

                     lng1 = l

*-----------------------------------------------------------------------
      if( nnm > MAX_NUM_NRST ) then
         write(ErrCha,'(a,a,i5,a,a,i5,a)')
     &        'sub.echrg3tallsm3.f'
     &           //' ?dimension over nrst?'
     &           //' nnm > MAX_NUM_NRST'
     &        ,' (nnm=',nnm,')'
     &        ,' (MAX_NUM_NRST@moddas.f=',MAX_NUM_NRST,')'
         ErrID = 'L:16001/R:echrg3/F:tallsm3.f'
         call ErrWrite(ErrID,ErrCha)
      endif
      call moddas_deallocate_int(nrst)

      return
      end


************************************************************************
*                                                                      *
      subroutine echtet(iot,mr,kr)
*                                                                      *
*       input echo for tally tetra mesh                                *
*       Last Modified by T.Furuta on 2025/01/16                        *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      integer,intent(in) :: iot,mr
      integer,intent(in) :: kr(mr)
      character dum*100

      ireg=kr(1)
      nr=kr(3)
      if(nr.gt.1.or.(ireg.gt.0.and.nr.gt.0))then
       dum(1:2)='( '
       l=3
       do ir=1,nr
        write(dum(l:l+6),'(i7)') kr(3+ir)
        l=l+7
        dum(l:l)=' '
        l=l+1
       enddo
       if(ireg.gt.0)then
        dum(l:l+1)='< '
        l=l+2
        write(dum(l:l+6),'(i7)') ireg
        l=l+7
       endif
       dum(l:l+1)=' )'
       l=l+1
      else
       l=1
       if(ireg.gt.0)then
        write(dum(l:l+6),'(i7)') ireg
       else
        write(dum(l:l+6),'(i7)') kr(4)
       endif
       l=l+6
      endif
      write(iot,'("      reg = ",100a1)')
     &     (dum(j:j),j=1,l)

*-----------------------------------------------------------------------
      return
      end


************************************************************************
*                                                                      *
      subroutine echrg_wgtsum(m,icc,ierl,iva,ivs,iot,mc,kc,nr,mr,kr,
     &                 vl,lr,nvl,ivl,rvl,igm)
*                                                                      *
*       input echo for tally region mesh                               *
*                                                                      *
*       modified subroutine 'echrg' to write                           *
*        "no cell operator ethres" sub-section and                     *
*        "cell cond0 cond1 ..." sub-section                            *
*       for 'reg=weightsum' in [T-Deposit]                             *
*                                                                      *
*       modified by S.Abe on 2016/11/24                                *
*                                                                      *
************************************************************************

      use tdepwgtsum_global

*-----------------------------------------------------------------------

      implicit real*8 (a-h,o-z)

      include 'param.inc'
      include 'ggsparam.inc'

      common /regdc/ idrg(kvlmax), idgr(kvmmax)

*-----------------------------------------------------------------------

      character dum1*10000
      character dum2*10000
      character chdf(1000)*100
      dimension ildf(1000)
      character cblan*100

      dimension kc(mc)
      dimension kr(mr)

      dimension vl(nr)
      dimension lr(nr)
      dimension ivl(nvl)
      dimension rvl(nvl)

      common /tall70/ itrwgtsum(itlmax), itrgn1(itlmax),
     &                itrgm1(itlmax), itrncd(itlmax)

*-----------------------------------------------------------------------
*     write reg = weightsum
*-----------------------------------------------------------------------

      write(iot,'("      reg = weightsum")')

      write(iot,'("    ncond = ",i4,11x,
     &            " # number of condition")') itrncd(m)

*-----------------------------------------------------------------------
*     check maximum length of cell and write index of condition
*-----------------------------------------------------------------------

      j = 0
      lgmax = 0
      do ir = 1, itrgn1(m)
       call echrg3(j,mr,kr,dum2,lng1,icmb,igm)
       if( lng1 .gt. lgmax ) lgmax = lng1 + 1
      enddo
      if( lgmax .lt. 7 ) lgmax = 7

      dum2(1:8) = '      no'
      l = 8
      if( lgmax .gt. 7 ) then
       do i = 1, lgmax-7
        dum2(l+i:l+i) = ' '
       end do
       l = l + lgmax - 7
      endif
      dum2(l+1:l+7) = '   cell'
      l = l + 7
      dum2(l+1:l+10) = '  operator'
      l = l + 10
      dum2(l+1:l+16) = '          ethres'
      l = l + 16
      dum2(l+1:l+8) = '    list'
      l = l + 8

      write(iot,'(600a1)') (dum2(i:i),i=1,l)

*-----------------------------------------------------------------------
*     write each condition
*-----------------------------------------------------------------------

      j = 0

      do ir = 1, ncond_old

       do k = 1, nadd_old

        if( itcond(m,1,ir,k) .ne. 0 ) then

         l = 0

         if( k .eq. 1 ) then
          write(dum2(l+1:l+8),'(1x,i7)') itcond(m,1,ir,k)
         else
          dum2(l+1:l+8) = '     and'
         endif
         l = l + 8

         call echrg3(j,mr,kr,dum1,lng1,icmb,igm)
         do i = 1, lgmax-lng1
          dum2(l+i:l+i) = ' '
         end do
         l = l + lgmax-lng1
         do i = 1, lng1
          dum2(l+i:l+i) = dum1(i:i)
         end do
         l = l + lng1

         if(     itcond(m,2,ir,k) .eq. 1 ) then
          dum2(l+1:l+10) = '     lt   '
         elseif( itcond(m,2,ir,k) .eq. 2 ) then
          dum2(l+1:l+10) = '     le   '
         elseif( itcond(m,2,ir,k) .eq. 3 ) then
          dum2(l+1:l+10) = '     eq   '
         elseif( itcond(m,2,ir,k) .eq. 4 ) then
          dum2(l+1:l+10) = '     ge   '
         elseif( itcond(m,2,ir,k) .eq. 5 ) then
          dum2(l+1:l+10) = '     gt   '
         endif
         l = l + 10

         write(dum2(l+1:l+16),'(1x,1p1e15.7)') rteth(m,ir,k)
         l = l + 16

         write(dum2(l+1:l+8),'(1x,i7)') itcond(m,3,ir,k)
         l = l + 8

         write(iot,'(600a1)') (dum2(i:i),i=1,l)

        endif

       enddo

      enddo

*-----------------------------------------------------------------------
*     write weightsum cell and efficiency list
*-----------------------------------------------------------------------

      j2 = j
      nr2 = nr - itrgn1(m)

      write(iot,'("    ncell = ",i4,11x,
     &            " # number of weightsum cell")') nr2

*-----------------------------------------------------------------------

      j = j2
      lgmax = 0
      do ir = 1, nr2
       call echrg3(j,mr,kr,dum2,lng1,icmb,igm)
       if( lng1 .gt. lgmax ) lgmax = lng1 + 1
      enddo
      if( lgmax .lt. 9 ) lgmax = 9

      l = 0
      if( lgmax .gt. 9 ) then
       do i = 1, lgmax-9
        dum2(l+i:l+i) = ' '
       end do
       l = l + lgmax - 9
      endif
      dum2(l+1:l+9) = '     cell'
      l = l + 9

      do k = 1, nlist_old
       if( itlist(m,k) .eq. -1 ) exit
       dum2(l+1:l+12) = '        list'
       l = l + 12
       write(dum2(l+1:l+4),'(i4.4)') itlist(m,k)
       l = l + 4
      enddo

      write(iot,'(600a1)') (dum2(i:i),i=1,l)

*-----------------------------------------------------------------------

      j = j2

      do ir = 1, nr2

       l = 0

       call echrg3(j,mr,kr,dum1,lng1,icmb,igm)
       do i = 1, lgmax-lng1
        dum2(l+i:l+i) = ' '
       end do
       l = l + lgmax-lng1
       do i = 1, lng1
        dum2(l+i:l+i) = dum1(i:i)
       end do
       l = l + lng1

       do k = 1, nlist_old
        if( itlist(m,k) .eq. -1 ) exit
        write(dum2(l+1:l+16),'(1x,1p1e15.7)') rteff(m,ir,k)
        l = l + 16
       enddo

       write(iot,'(600a1)') (dum2(i:i),i=1,l)

      enddo

*-----------------------------------------------------------------------
*        echo for combined regions
*-----------------------------------------------------------------------

         if( ivs .eq. 0 .or. ( iva .eq. 0 .and.
     &       junf .eq. 0 .and. icmb .eq. 0 .and. nvl .eq. 0 ) ) return

         if( iva .eq. 0 ) then

            write(iot,'("   volume  ",17x,
     &                  "# combined, lattice or level structure ")')
            write(iot,'("   non     reg      vol     #",
     &                  " reg definition")')

         else

            write(iot,'("   value  ",18x,
     &                  "# values for each region")')
            write(iot,'("   non     reg      val     #",
     &                  " reg definition")')

         end if

               call tregvl(nr,mr,kr,vl,lr,nvl,ivl,rvl)

*-----------------------------------------------------------------------
*     read from kr(mr)
*-----------------------------------------------------------------------

                  j = j2

         do 500 ir = itrgn1(m)+1, nr

               call echrg3(j,mr,kr,dum2,lng1,icmb,igm)

                  ild1 = 30
                  ild0 = max(10,ierl-ild1)

                  ilrm = lng1
                  isqd = 0
                  isrm = 0

  530             isqd = isqd + 1

                  if( ilrm .le. ild0 ) then

                        ildf(isqd) = ilrm

                     do k = 1, ilrm

                        chdf(isqd)(k:k) = dum2(k+isrm:k+isrm)

                     end do

                  else

                     do k = ild0, 1, -1

                        if( dum2(k+isrm:k+isrm) .eq. ' ' ) goto 520

                     end do

  520                k2 = k

                        ildf(isqd) = k2 - 1

                     do k = 1, k2 - 1

                        chdf(isqd)(k:k) = dum2(k+isrm:k+isrm)

                     end do

                     isrm = isrm + k2
                     ilrm = ilrm - k2

                     goto 530

                  end if

                  write(iot,'(i5,2x,i7,1p1e13.4," # ",100a1)')
     &                 ir-itrgn1(m), lr(ir), vl(ir),
     &                 (chdf(1)(l:l),l=1,ildf(1))

               if( isqd .gt. 1 ) then

                  do k = 2, isqd

                     write(iot,'(100a1)') (cblan(l:l),l=1,ild1-2),
     &                    '#',' ',(chdf(k)(l:l),l=1,ildf(k))

                  end do

               end if

*-----------------------------------------------------------------------

  500    continue

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine usrdfech(iot,m,iaxs,icax)
*                                                                      *
*       input echo of the user defined tally                           *
*       last modified by S.Hashimoto on 2017/11/01                     *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------

      use usrtalmod
      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'

*-----------------------------------------------------------------------

      common /tall11/ itfln(itlmax), itfll(itlmax,6), ctfln(itlmax,6)
      character ctfln*100
      common /udtall1/ iudtfll(itlmax,50), cudtfln(itlmax,50)
      character cudtfln*100

      common /cudtpara/ udtpara(0:9)
      common /ciudtpara/ iudtpara(0:9)

      integer nudtvar
      common /cnudtvar/ nudtvar

*-----------------------------------------------------------------------

      character asfil*100

*-----------------------------------------------------------------------
      common /tall79/ itallech
      if( itallech .eq. 0 ) return

*-----------------------------------------------------------------------
*     input echo
*-----------------------------------------------------------------------

            write(iot,'("[ T-Userdefined ]")')

*-----------------------------------------------------------------------

                  infil = itfln(m)

            do i = 1, infil

                  msfile = max( 14, iudtfll(m,i) ) ! S.H. revised (2016.12.21)
                  asfil = '  # file name of output'

                  write(iot,'("     file = ",100a1)')
     &                  ( cudtfln(m,i)(j:j), j = 1, msfile ), ! S.H. revised (2016.12.21)
     &                  ( asfil(j:j), j = 1, 23 )

            end do

*-----------------------------------------------------------------------

            do i = 0, 9

               if( iudtpara(i) .eq. 1 ) then

                  write(iot,'(" udtpara",i1," =",1p1g14.7,1x,1x,
     &            " # (D=0.0) parameters for user defined tally")')
     &            i,udtpara(i)

               end if

            end do

*-----------------------------------------------------------------------

            if( nudtvar .gt. 0 ) then

               write(iot,'("  nudtvar =",1x,i5,10x,
     &        " # (D=0) number of variables for user defined tally")')
     &         nudtvar

               do i = 1, nudtvar

                  if( iudtvar(i) .eq. 1 ) then

                     write(iot,'(5x,"udtvar(",i3,") =",1p1g14.7)')
     &                    i,udtvar(i)

               end if

            end do
            end if

*-----------------------------------------------------------------------

      return
      end


************************************************************************
      subroutine terrang(iot,m)
*                                                                      *
*       write angel para. for error file                               *
*       last modified by K. Niita on 2017/08/16                        *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'

      common /tall15/ itanl(itlmax), itang(itlmax)
      character itang*200
      character jtang*200
      character ktang*200

*-----------------------------------------------------------------------

         do i = 1, 200
            if( i .le. itanl(m) ) then
               jtang(i:i) = itang(m)(i:i)
            else
               jtang(i:i) = ' '
            end if
         end do

            call chcaps(jtang,1,itanl(m),i3,'#!$')

            l = 0
            i = 0

  200       i = i + 1
            if( i .gt. itanl(m) ) goto 300

            if( jtang(i:i+3) .eq. 'cmin' .or.
     &          jtang(i:i+3) .eq. 'cmax' .or.
     &          jtang(i:i+3) .eq. 'dmin' .or.
     &          jtang(i:i+3) .eq. 'dmax' .or.
     &          jtang(i:i+3) .eq. 'icut' .or.
     &          jtang(i:i+3) .eq. 'cols' .or.
     &          jtang(i:i+3) .eq. 'cuts' ) then

               do k = i + 4, itanl(m)
                  if( jtang(k:k) .eq. ' ' ) goto 100
               end do
  100             i = k - 1

            else

               l = l + 1
               ktang(l:l) = itang(m)(i:i)

            end if

         goto 200
  300    continue

                  write(iot,'( "p: ",200a1)')
     &            ( ktang(i:i),i = 1, l )


*-----------------------------------------------------------------------

      return
      end

