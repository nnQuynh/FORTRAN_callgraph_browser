************************************************************************
*                                                                      *
      subroutine whone(icw,jof,jhf,noned,
     &                 ixlog,iylog,xmax,xmin,ymax,ymin,
     &                 ierr,clal,clmo,idbg,xfac,rlptl,sybw,erwd)
cKN 2024/01/24
*                                                                      *
*                                                                      *
*        PURPOSE    :  PRE PROCEDURE OF LINES READ FROM JOL            *
*                                                                      *
*             ICW   :  0-> NORMAL, 1-> FOR INTERIA REGION              *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

*     dimension xyp(6,1)
      real(8),allocatable:: xyp(:,:)

      character m_err*200
      common /error/ m_err, l_err, k_err

      common /frm/  xal, yal

      common /fdfn/ ifd(0:13),ikan1,ikan2,ikan3,ikan4,
     &              idsm(9), idln(6), idxe, idye, ibd1, ibd2, ibd3,
     &              ibd4, ibd5, ibd6, ibd7, ibd8, ibd9, ibd10, ibd11

      dimension clal(3), rcol(3), rcob(3)

*-----------------------------------------------------------------------


      do 300 ione = 1, noned

         read(jof) ild, idyr, ixdc, iydc

         read(jof) ityp,iwi,ipol,imar,
     &             rcol,rcob,icoma,facsiz

         allocate(xyp(6,max(ild,1)))

         if( ityp .ne. 11 ) then

               if( clal(1) .gt. -r0max ) then
                  rcol(1) = clal(1)
                  rcol(2) = clal(2)
                  rcol(3) = clal(3)
               end if

               if( clmo .gt. -r0max .and. rcol(1) .gt. 0.0 )
     &             rcol(1) = -2.0
               if( rcol(1) .lt. -2.5 ) rcol(1) = -2.0

            if( imar .eq.  3 .or. imar .eq.  5 .or.
     &          imar .eq.  7 .or. imar .eq.  9 .or.
     &          imar .eq. 11 .or. imar .eq. 13 ) then

               if( clal(1) .gt. -r0max ) rcob(1) = -1.0

               if( clmo .gt. -r0max .and. rcob(1) .gt. 0.0 )
     &             rcob(1) = -1.0

               if( clmo .gt. -r0max .and. rcob(1) .lt. -r0max )
     &             rcob(1) = -1.0

               if( rcob(1) .lt. -2.5 ) rcob(1) = -1.0

            end if

         else if( ityp .eq. 11 ) then

               if( clmo .gt. -r0max .and. rcol(1) .gt. 0.0 )
     &         call ctomo(rcol,clmo)

               if( rcol(1) .lt. -2.5 ) rcol(1) = -1.0

         end if


*-----------------------------------------------------------------------
*        WRITE DEFINITION OF LINES
*-----------------------------------------------------------------------

            if( icoma .ne. 1 .and. ityp .ge. 1 .and. ityp .le. 6 ) then

               if( idln(ityp) .eq. 0 ) then

                  call lndfn(jhf,ityp,idbg)

                  idln(ityp) = 1

               end if

            end if

*-----------------------------------------------------------------------
*        WRITE DEFINITION OF SYMBOLS
*-----------------------------------------------------------------------

            if( imar .gt. 0 ) then

                  call smdfn(jhf,imar,idbg)

            end if

*-----------------------------------------------------------------------

                     rleng = 0.0

         do 310 ils = 1, ild

               read(jof) xpo,ypo,barxl,barxr,baryl,baryu

               if( barxl .lt. 0.0 .or.
     &             barxr .lt. 0.0 .or.
     &             baryl .lt. 0.0 .or.
     &             baryu .lt. 0.0 ) then

                     m_err = 'Erorr Numbers Should be Positive'
                     ErrCha = ''
                     ErrID = 'L:129/R:whone/F:a-line.f'
                     goto 999

               end if

               if( ixlog .eq. 1 ) then

                  if( xpo .lt. 0.0 ) then


                     xpo = r9min

                  else if( xpo .eq. 0.0 ) then

                     xpo = r9min

                  end if

                  if( xpo - barxl .le. 0.0 ) barxl = xpo - xmin * 0.1

               end if

               if( iylog .eq. 1 ) then

                  if( ypo .lt. 0.0 ) then


                     ypo = r9min

                  else if( ypo .eq. 0.0 ) then

                     ypo = r9min

                  end if

                  if( ypo - baryl .le. 0.0 ) baryl = ypo - ymin * 0.1

               end if


                     xyp(1,ils) = zonep(xmin,xmax,xpo,ixlog)
                     xyp(2,ils) = zonep(ymin,ymax,ypo,iylog)


                  if( ixdc .eq. 1 ) then

                     xyp(3,ils) = zonep(xmin,xmax, xpo - barxl ,ixlog)

                  else if( ixdc .eq. 2 ) then

                     xyp(4,ils) = zonep(xmin,xmax, xpo + barxr ,ixlog)

                  else if( ixdc .eq. 3 ) then

                     xyp(3,ils) = zonep(xmin,xmax, xpo - barxl ,ixlog)
                     xyp(4,ils) = zonep(xmin,xmax, xpo + barxr ,ixlog)

                  end if


                  if( iydc .eq. 1 ) then

                     xyp(5,ils) = zonep(ymin,ymax, ypo - baryl ,iylog)

                  else if( iydc .eq. 2 ) then

                     xyp(6,ils) = zonep(ymin,ymax, ypo + baryu ,iylog)

                  else if( iydc .eq. 3 ) then

                     xyp(5,ils) = zonep(ymin,ymax, ypo - baryl ,iylog)
                     xyp(6,ils) = zonep(ymin,ymax, ypo + baryu ,iylog)

                  end if

                  if( ils .gt. 1 .and.
     &                icoma .eq. 0 .and.
     &                ityp .ne. 0 ) then

                     rleng = rleng + sqrt(
     &                       ( xyp(1,ils) - xyp(1,ils-1) )**2 * xal**2
     &                     + ( xyp(2,ils) - xyp(2,ils-1) )**2 * yal**2 )

                  end if

  310    continue


*-----------------------------------------------------------------------

      if( icw .eq. 0 ) then

         if( ityp .ne. 11 ) then

            call preline(jhf,ild,ixdc,iydc,rleng,idbg,
     &                   ityp,iwi,ipol,imar,rcol,rcob,icoma,facsiz,
     &                   xyp,xfac,sybw,rlptl,erwd)
cKN 2024/01/24

         else

            call wline01(2,jhf,ild,ixdc,iydc,rleng,idbg,
     &                   ityp,iwi,ipol,imar,rcol,rcob,icoma,facsiz,
     &                   xyp,xfac,rlptl)

         end if

      else

            call wline01(0,jhf,ild,ixdc,iydc,rleng,idbg,
     &                   ityp,iwi,ipol,imar,rcol,rcob,icoma,facsiz,
     &                   xyp,xfac,rlptl)

      end if

*-----------------------------------------------------------------------

      deallocate(xyp)

  300 continue

      return
  999 ierr=1
      deallocate(xyp)
      return
      end

************************************************************************
*                                                                      *
      subroutine preline(jhf,ild,ixdc,iydc,rleng,idbg,
     &                   ityp,iwi,ipol,imar,rcol,rcob,icoma,facsiz,
     &                   xyp,xfac,sybw,rlptl,erwd)
cKN 2024/01/24
*                                                                      *
*                                                                      *
*        PURPOSE    :  PRE PROCESS FOR WRITE LINES ON jhf              *
*                                                                      *
*            IMAR  = 0 NO SYMBOL                                       *
*            ICOMA = 1 NO LINE                                         *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'

*-----------------------------------------------------------------------

      dimension xyp(6,1)
      dimension rcol(3), rcob(3)

*-----------------------------------------------------------------------
*        MARKER
*-----------------------------------------------------------------------

         if( imar .gt. 0 ) then

            if( imar .eq. 17 .or. imar .eq. 18 .or.
     &          imar .eq. 19 .or. imar .eq. 20 .or.
     &          imar .eq. 21 .or. imar .eq. 22 ) then

               if( icoma .ne. 1 ) then

*                 ---- CLIP ALL MARKERS AND DRAW LINE ----

                     call wsymb01(1,jhf,ild,ixdc,iydc,rleng,idbg,
     &                    ityp,iwi,ipol,imar,rcol,rcob,icoma,facsiz,
     &                    xyp,xfac,sybw)

                     call wline01(1,jhf,ild,ixdc,iydc,rleng,idbg,
     &                    ityp,iwi,ipol,imar,rcol,rcob,icoma,facsiz,
     &                    xyp,xfac,rlptl)

               end if


               if( ixdc .ne. 0 .or. iydc .ne. 0 ) then

*                 ---- CLIP ONE MARKER AND DRAW ERROR BARS ----

                     call weror01(1,jhf,ild,ixdc,iydc,rleng,idbg,
     &                    ityp,iwi,ipol,imar,rcol,rcob,icoma,facsiz,
     &                    xyp,xfac,erwd)
cKN 2024/01/24

               end if

            else

               if( icoma .ne. 1 ) then

*                 ---- DRAW LINE ----

                     call wline01(0,jhf,ild,ixdc,iydc,rleng,idbg,
     &                    ityp,iwi,ipol,imar,rcol,rcob,icoma,facsiz,
     &                    xyp,xfac,rlptl)

               end if

               if( ixdc .ne. 0 .or. iydc .ne. 0 ) then

*                 ---- DRAW ERROR BARS ----

                     call weror01(0,jhf,ild,ixdc,iydc,rleng,idbg,
     &                    ityp,iwi,ipol,imar,rcol,rcob,icoma,facsiz,
     &                    xyp,xfac,erwd)
cKN 2024/01/24

               end if

            end if

            if( imar .gt. 0 ) then

*                 ---- DRAW MAKERS ----

                     call wsymb01(0,jhf,ild,ixdc,iydc,rleng,idbg,
     &                    ityp,iwi,ipol,imar,rcol,rcob,icoma,facsiz,
     &                    xyp,xfac,sybw)

            end if


         else

               if( icoma .ne. 1 ) then

*                 ---- DRAW LINE ----

                     call wline01(0,jhf,ild,ixdc,iydc,rleng,idbg,
     &                    ityp,iwi,ipol,imar,rcol,rcob,icoma,facsiz,
     &                    xyp,xfac,rlptl)

               end if

               if( ixdc .ne. 0 .or. iydc .ne. 0 ) then

*                 ---- DRAW ERROR BARS ----

                     call weror01(0,jhf,ild,ixdc,iydc,rleng,idbg,
     &                    ityp,iwi,ipol,imar,rcol,rcob,icoma,facsiz,
     &                    xyp,xfac,erwd)
cKN 2024/01/24

               end if

         end if


*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      subroutine wline01(icw,jhf,ild,ixdc,iydc,rleng,idbg,
     &                   ityp,iwi,ipol,imar,rcol,rcob,icoma,facsiz,
     &                   xyp,xfac,rlptl)
*                                                                      *
*                                                                      *
*        PURPOSE    :  WRITE LINES ON jhf                              *
*                                                                      *
*              ICW  : 0 -> FIRST, CLIP THE AXIS, DRAW AND GRESAVE      *
*                   : 1 -> AFTER CLIP OF SYMBOLS AND AXIS,             *
*                          DRAW AND GRESAVE                            *
*                   : 2 -> WITHOUT CLIP                                *
*                   : 3 -> AFTER CLIP OF SYMBOLS, WITHOUT AXIS CLIP    *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'

*-----------------------------------------------------------------------

      parameter ( slu0 = 1.70 )
      parameter ( smin = 0.65 )
      parameter ( smmm = 3.5277773d-05 )
      parameter ( ddcm = 8.4666667d-3 )

      dimension xyp(6,*)
      dimension sls0(6)
      dimension rcol(3), rcob(3), rrcol(3)

      data sls0/ 0.375, 0.360, 0.440, 0.001, 0.3125, 0.25/

*-----------------------------------------------------------------------

               if( icw .eq. 0 ) then

                  write(jhf,'(''gs ax cl np'')')

               end if

*-----------------------------------------------------------------------

            iwi = max(1,iwi)

            if( ityp .ne. 0 .and. ityp .ne. 11 ) then

                     slu = slu0

               if( rlptl .gt. -r0max ) then

                  if( rlptl .lt. 0.0 ) then

                     slu = slu * sqrt( dble(iwi) / 4.0 )

                  else if( rlptl .gt. 0.0 ) then

                     slu = slu * rlptl

                  end if

               end if

                  if( ityp .eq. 2 ) then

                     slu = slu / 2.0

                  else if( ityp .eq. 4 ) then

                     slu = slu / 4.0

                  end if


                  if( rleng .lt. smin * slu ) then

                     slu = smin * slu

                  else if( rleng .lt. slu ) then

                     slu = rleng

                  end if


                  if( ityp .eq. 4 ) then

                        rleng = rleng - smmm * 1000.0
                        sls   = smmm
                        slu   = slu + sls

                  else

                     if( ityp .eq. 2 .or. ityp .eq. 3 ) then

                        sdw0 = dble( 4 ) * ddcm / 2.0
                        sdw  = dble(iwi) * ddcm / 2.0

                     else

                        sdw0 = dble( 4 ) * ddcm
                        sdw  = dble(iwi) * ddcm

                     end if

                        sls = sls0( ityp ) * slu - sdw0


                     if( ityp .eq. 1 ) then

                        slu = slu - smmm

                     else if( ityp .eq. 5 ) then

                        slu = slu - 2.0 * smmm

                     else if( ityp .eq. 6 ) then

                        slu = slu - 3.0 * smmm

                     end if

                  end if


               if( rleng .gt. slu ) then

                  if( xyp(1,1) .eq. xyp(1,ild) .and.
     &                xyp(2,1) .eq. xyp(2,ild) .and.
     &                ityp .ne. 4 ) then

                     nlu =  nint( rleng / ( slu - sls ) )

                  else

                     nlu =  nint( ( rleng - sls ) / ( slu - sls ) )

                  end if


                  if( ityp .eq. 4 ) then

                     rlu = ( rleng - sls ) / dble( nlu ) + sls

                  else


                     if( xyp(1,1) .eq. xyp(1,ild) .and.
     &                   xyp(2,1) .eq. xyp(2,ild) ) then

                        rlu = ( rleng - dble( nlu ) * sdw )
     &                      / ( dble( nlu ) * ( 1.0 - sls0( ityp ) ) )

                     else

                        rlu = ( rleng - dble( nlu - 1 ) * sdw )
     &                      / ( dble( nlu ) * ( 1.0 - sls0( ityp ) )
     &                          + sls0( ityp ) )

                     end if

                  end if


               else

                     rlu = slu

               end if

                  write(jhf,'(''/lu {'',e16.8,'' cm mul } N '',
     &                        ''/dw '',I3,'' dd N'')') rlu, iwi

            end if

*-----------------------------------------------------------------------

               if( ityp .ne. 11 ) then

                  itsd = ityp

               else

                  itsd = 0

               end if


                  if( rcol(1) .lt. -r0max ) then

                     rrcol(1) = -2.0
                     rrcol(2) =  1.0
                     rrcol(3) =  1.0

                  else

                     rrcol(1) = rcol(1)
                     rrcol(2) = rcol(2)
                     rrcol(3) = rcol(3)

                  end if


                  if( itsd .ne. 0 ) then

                     write(jhf,'(3f7.3,'' sc sd'',i1,'' dw lw'')')
     &                          rrcol, itsd

                  else

                     write(jhf,'(3f7.3,'' sc sd'',i1,i4,'' dd lw'')')
     &                          rrcol, itsd, iwi

                  end if


*-----------------------------------------------------------------------

            do 100 i = 1, ild

               if( i .eq. 1 ) then

                     write(jhf,'(2g14.5,'' m'')')
     &                           xyp(1,i),xyp(2,i)

               else if( i .eq. ild .and. ityp .eq. 0 .and.
     &                  xyp(1,1) .eq. xyp(1,ild) .and.
     &                  xyp(2,1) .eq. xyp(2,ild) ) then

                     write(jhf,'(2g14.5,'' l cp'')')
     &                           xyp(1,i),xyp(2,i)

               else

                     write(jhf,'(2g14.5,'' l'')')
     &                           xyp(1,i),xyp(2,i)

               end if

                        xp1 = xyp(1,i)
                        yp1 = xyp(2,i)

                     if( icw .eq. 0 .or. icw .eq. 1 ) then

                        xpm = max(0.0d0,xp1)
                        xp1 = min(1.0d0,xpm)
                        ypm = max(0.0d0,yp1)
                        yp1 = min(1.0d0,ypm)

                     end if

                        call bbox(3,0,xp1,yp1,0.d0)

  100       continue

*-----------------------------------------------------------------------


            if( icw .eq. 0 .or. icw .eq. 1 .or. icw .eq. 3 ) then

                  if( ityp .ne. 11 ) then

                     write(jhf,'(''st gr np'')')

                  else if( ityp .eq. 11 ) then

                     write(jhf,'(''cp fl gr np'')')

                  end if

            else

                  if( ityp .ne. 11 ) then

                     write(jhf,'(''st'')')

                  else if( ityp .eq. 11 ) then

                     write(jhf,'(''cp fl'')')

                  end if

            end if

*-----------------------------------------------------------------------

            if( idbg .eq. 1 )  write(jhf,'()')

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      subroutine weror01(icw,jhf,ild,ixdc,iydc,rleng,idbg,
     &                   ityp,iwi,ipol,imar,rcol,rcob,icoma,facsiz,
     &                   xyp,xfac,erwd)
cKN 2024/01/24
*                                                                      *
*                                                                      *
*        PURPOSE    :  WRITE ERROR BARS ON jhf                         *
*                                                                      *
*              ICW  : 0 -> DRAW ERROR BARS INSIDE THE AXIS             *
*                   : 1 -> DRAW ERROR BARS AFTER CLIP THE SYMBOL       *
*                          INSIDE THE AXIS                             *
*                   : 2 -> DRAW ERROR BARS                             *
*                   : 3 -> DRAW ERROR BARS AFTER CLIP THE SYMBOL       *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'

*-----------------------------------------------------------------------

      parameter (ddcm = 8.4666667d-3)

      parameter ( ss0 = 0.254 )
      parameter ( il0 = 4 )

*-----------------------------------------------------------------------

      common /frm/  xal, yal

      dimension xyp(6,*)

      dimension sfac(22), sfxm(22), sfxp(22), sfym(22), sfyp(22)
      character*3 sdif(22)

      data sfac/ 0.5,  1.0, 1.0, -1.0,  0.9,    -0.9,    1.2, -1.2,
     &           1.2, -1.2, 0.9, -0.9,  0.6364, -0.6364, 0.7,  2.5,
     &           1.0,  0.9, 1.2,  1.2,  0.9,     0.6364/

      data sfxm/ 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
     &           0.5, 0.5, 0.5, 0.5, 1.0, 1.0, 0.5, 0.2,
     &           0.5, 0.5, 0.5, 0.5, 0.5, 1.0/

      data sfxp/ 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
     &           0.5, 0.5, 0.5, 0.5, 1.0, 1.0, 0.5, 0.2,
     &           0.5, 0.5, 0.5, 0.5, 0.5, 1.0/

      data sfym/ 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.288675, 0.288675,
     &           0.57735, 0.57735, 0.866025, 0.866025,
     &           1.0, 1.0, 0.5, 0.2,
     &           0.5, 0.5, 0.288675, 0.57735, 0.866015, 1.0/

      data sfyp/ 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.57735, 0.57735,
     &           0.288675, 0.288675, 0.866025, 0.866025,
     &           1.0, 1.0, 0.5, 0.2,
     &           0.5, 0.5, 0.288675, 0.57735, 0.866015, 1.0/

      data sdif/ 'cir','pls','cir','cir','squ','squ','tr1','tr1',
     &           'tr2','tr2','di1','di1','di2','di2','crs','ast',
     &           'cir','squ','tr1','tr2','di1','di2'/


      common /fdfn/ ifd(0:13),ikan1,ikan2,ikan3,ikan4,
     &              idsm(9), idln(6), idxe, idye, ibd1, ibd2, ibd3,
     &              ibd4, ibd5, ibd6, ibd7, ibd8, ibd9, ibd10, ibd11

      dimension rcol(3), rcob(3), rrcol(3)

*-----------------------------------------------------------------------
*        FIRST DEFINITION OF ERROR OPERATORS DX AND DY
*        AND CLIPPING PATH BX AND BY
*-----------------------------------------------------------------------

*                 bx : x inside of the axis and y is free
*                 by : y inside of the axis and x is free

*-----------------------------------------------------------------------

            if( idxe .eq. 0 .and. idye .eq. 0 .and.
     &        ( ixdc .ne. 0 .or. iydc .ne. 0 ) ) then

                  write(jhf,'(
     &            ''/dx {/yps S N yps m yps l} N ''
     &            ''/dy {/xps S N xps S m xps S l} N''
     &            )')

            end if

            if( idxe .eq. 0 .and. ixdc .ne. 0 ) then

                  write(jhf,'(
     &            ''/bx { 0 yneg M 0 ypst L xal ypst '',
     &            ''L xal yneg L cp} N''
     &            )')

                  if( idbg .eq. 1 )  write(jhf,'()')

                  idxe = 1

            end if

            if( idye .eq. 0 .and. iydc .ne. 0 ) then

                  write(jhf,'(
     &            ''/by { xneg 0 M xneg yal L xpst '',
     &            ''yal L xpst 0 L cp} N''
     &            )')

                  if( idbg .eq. 1 )  write(jhf,'()')

                  idye = 1

            end if

*-----------------------------------------------------------------------

               il = il0

            if( facsiz .gt. 1.0 ) then

               il = nint( facsiz * dble( il0 ) ) - 1

            else if( facsiz .le. 0.51 ) then

               il = 1

            else if( facsiz .le. 0.61 ) then

               il = 2

            else if( facsiz .le. 0.81 ) then

               il = 3

            end if

            if( imar .eq. -1 ) il = iwi

*-----------------------------------------------------------------------

                  ss  = ss0 * facsiz

cKN 2024/01/24
                  ss = ss * erwd

                  sxl = ss * 0.75 / 2.0 / xal
                  syl = ss * 0.75 / 2.0 / yal


            if( imar .gt. 0 ) then

               if( sfac(imar) .gt. 0.0 ) then

                  ss = ss * sfac(imar)

               else if( sfac(imar) .gt. 0.0 ) then

                  ss = - ss * sfac(imar) + dble(il) * ddcm

               end if

               if( imar .eq. 14 ) then

                  ss = ss / 1.41421

               end if

                  sxm = ( ss * sfxm(imar) + dble(il) * ddcm ) / xal
                  sxp = ( ss * sfxp(imar) + dble(il) * ddcm ) / xal
                  sym = ( ss * sfym(imar) + dble(il) * ddcm ) / yal
                  syp = ( ss * sfyp(imar) + dble(il) * ddcm ) / yal

            else

                  sxm = 0.0
                  sxp = 0.0
                  sym = 0.0
                  syp = 0.0

            end if

*-----------------------------------------------------------------------

               if( rcol(1) .lt. -r0max ) then

                  rrcol(1) = -2.0
                  rrcol(2) =  1.0
                  rrcol(3) =  1.0

               else

                  rrcol(1) = rcol(1)
                  rrcol(2) = rcol(2)
                  rrcol(3) = rcol(3)

               end if


            if( icw .eq. 1 .or. icw .eq. 3 ) then

                  write(jhf,'(3f7.3,'' sc /sm '',g14.5,'' cm N sd0 '',
     &                        i2,'' dd lw'')')  rrcol, ss, il


                  if( imar .eq. 16 ) then

                     write(jhf,'(''/sm '',e14.6,'' cm N smc'')') ss

                  end if


            else if( icw .eq. 0 .or. icw .eq. 2 ) then


               write(jhf,'(3f7.3,'' sc sd0 '',i2,'' dd lw'')') rrcol, il


            end if

*-----------------------------------------------------------------------

         do 100 i = 1, ild

            if( ixdc .eq. 3 ) then

               if( ( ( xyp(4,i) .ge. 0.0 .and.
     &                 xyp(3,i) .le. 1.0 .and.
     &                 xyp(2,i) .ge. 0.0 .and. xyp(2,i) .le. 1.0 .and.
     &               ( icw .eq. 0 .or. icw .eq. 1 ) ) .or.
     &               ( icw .eq. 2 .or. icw .eq. 3 ) ) .and.
     &               ( xyp(3,i) .lt. xyp(1,i) - sxm .or.
     &                 xyp(4,i) .gt. xyp(1,i) + sxp ) ) then

                  if( idbg .eq. 1 )  write(jhf,'()')

                  if( icw .eq. 1 ) then

                        write(jhf,'(''gs bx cl '',2g14.5,
     &                              '' '',a3,'' cl np'')')
     &                              xyp(1,i),xyp(2,i),sdif(imar)

                  else if( icw .eq. 3 ) then

                        write(jhf,'(''gs ba cl '',2g14.5,
     &                              '' '',a3,'' cl np'')')
     &                              xyp(1,i),xyp(2,i),sdif(imar)

                  else if( icw .eq. 0 ) then

                        write(jhf,'(''gs bx cl np'')')

                  end if


                  if( xyp(3,i) .lt. xyp(1,i) - sxm .and.
     &                xyp(4,i) .gt. xyp(1,i) + sxp ) then

                        write(jhf,'(3g14.5,'' dx'')')
     &                              xyp(3,i),xyp(4,i),xyp(2,i)
                        write(jhf,'(3g14.5,'' dy'')')
     &                  xyp(2,i) + syl, xyp(2,i) - syl, xyp(3,i)

                     if( icw .eq. 2 ) then

                        write(jhf,'(3g14.5,'' dy st'')')
     &                  xyp(2,i) + syl, xyp(2,i) - syl, xyp(4,i)

                     else

                        write(jhf,'(3g14.5,'' dy st gr np'')')
     &                  xyp(2,i) + syl, xyp(2,i) - syl, xyp(4,i)

                     end if


                        xp1 = xyp(3,i)
                        xp2 = xyp(4,i)
                        yp1 = xyp(2,i) + syl
                        yp2 = xyp(2,i) - syl

                     if( icw .eq. 0 .or. icw .eq. 1 ) then

                        xpm = max(0.0d0,xp1)
                        xp1 = min(1.0d0,xpm)
                        xpm = max(0.0d0,xp2)
                        xp2 = min(1.0d0,xpm)

                     end if

                        call bbox(3,0,xp1,yp1,0.d0)
                        call bbox(3,0,xp1,yp2,0.d0)
                        call bbox(3,0,xp2,yp1,0.d0)
                        call bbox(3,0,xp2,yp2,0.d0)


                  else if( xyp(3,i) .lt. xyp(1,i) - sxm ) then

                        write(jhf,'(3g14.5,'' dx'')')
     &                              xyp(3,i),xyp(1,i),xyp(2,i)

                     if( icw .eq. 2 ) then

                        write(jhf,'(3g14.5,'' dy st'')')
     &                  xyp(2,i) + syl, xyp(2,i) - syl, xyp(3,i)

                     else

                        write(jhf,'(3g14.5,'' dy st gr np'')')
     &                  xyp(2,i) + syl, xyp(2,i) - syl, xyp(3,i)

                     end if


                        xp1 = xyp(3,i)
                        yp1 = xyp(2,i) + syl
                        yp2 = xyp(2,i) - syl

                     if( icw .eq. 0 .or. icw .eq. 1 ) then

                        xpm = max(0.0d0,xp1)
                        xp1 = min(1.0d0,xpm)

                     end if

                        call bbox(3,0,xp1,yp1,0.d0)
                        call bbox(3,0,xp1,yp2,0.d0)


                  else if( xyp(4,i) .gt. xyp(1,i) + sxp ) then

                        write(jhf,'(3g14.5,'' dx'')')
     &                              xyp(1,i),xyp(4,i),xyp(2,i)

                     if( icw .eq. 2 ) then

                        write(jhf,'(3g14.5,'' dy st'')')
     &                  xyp(2,i) + syl, xyp(2,i) - syl, xyp(4,i)

                     else

                        write(jhf,'(3g14.5,'' dy st gr np'')')
     &                  xyp(2,i) + syl, xyp(2,i) - syl, xyp(4,i)

                     end if


                        xp2 = xyp(4,i)
                        yp1 = xyp(2,i) + syl
                        yp2 = xyp(2,i) - syl

                     if( icw .eq. 0 .or. icw .eq. 1 ) then

                        xpm = max(0.0d0,xp2)
                        xp2 = min(1.0d0,xpm)

                     end if

                        call bbox(3,0,xp2,yp1,0.d0)
                        call bbox(3,0,xp2,yp2,0.d0)


                  end if

               end if

            else if( ixdc .eq. 2 ) then

               if( ( ( xyp(4,i) .ge. 0.0 .and.
     &                 xyp(1,i) .le. 1.0 .and.
     &                 xyp(2,i) .ge. 0.0 .and. xyp(2,i) .le. 1.0 .and.
     &               ( icw .eq. 0 .or. icw .eq. 1 ) ) .or.
     &               ( icw .eq. 2 .or. icw .eq. 3 ) ) .and.
     &                 xyp(4,i) .gt. xyp(1,i) + sxp ) then

                  if( idbg .eq. 1 )  write(jhf,'()')

                  if( icw .eq. 1 ) then

                        write(jhf,'(''gs bx cl '',2g14.5,
     &                              '' '',a3,'' cl np'')')
     &                              xyp(1,i),xyp(2,i),sdif(imar)

                  else if( icw .eq. 3 ) then

                        write(jhf,'(''gs ba cl '',2g14.5,
     &                              '' '',a3,'' cl np'')')
     &                              xyp(1,i),xyp(2,i),sdif(imar)

                  else if( icw .eq. 0 ) then

                        write(jhf,'(''gs bx cl np'')')

                  end if

                        write(jhf,'(3g14.5,'' dx'')')
     &                              xyp(1,i),xyp(4,i),xyp(2,i)

                     if( icw .eq. 2 ) then

                        write(jhf,'(3g14.5,'' dy st'')')
     &                  xyp(2,i) + syl, xyp(2,i) - syl, xyp(4,i)

                     else

                        write(jhf,'(3g14.5,'' dy st gr np'')')
     &                  xyp(2,i) + syl, xyp(2,i) - syl, xyp(4,i)

                     end if


                        xp2 = xyp(4,i)
                        yp1 = xyp(2,i) + syl
                        yp2 = xyp(2,i) - syl

                     if( icw .eq. 0 .or. icw .eq. 1 ) then

                        xpm = max(0.0d0,xp2)
                        xp2 = min(1.0d0,xpm)

                     end if

                        call bbox(3,0,xp2,yp1,0.d0)
                        call bbox(3,0,xp2,yp2,0.d0)


               end if

            else if( ixdc .eq. 1 ) then

               if( ( ( xyp(1,i) .ge. 0.0 .and.
     &                 xyp(3,i) .le. 1.0 .and.
     &                 xyp(2,i) .ge. 0.0 .and. xyp(2,i) .le. 1.0 .and.
     &               ( icw .eq. 0 .or. icw .eq. 1 ) ) .or.
     &               ( icw .eq. 2 .or. icw .eq. 3 ) ) .and.
     &                 xyp(3,i) .lt. xyp(1,i) - sxm ) then

                  if( idbg .eq. 1 )  write(jhf,'()')

                  if( icw .eq. 1 ) then

                        write(jhf,'(''gs bx cl '',2g14.5,
     &                              '' '',a3,'' cl np'')')
     &                              xyp(1,i),xyp(2,i),sdif(imar)

                  else if( icw .eq. 3 ) then

                        write(jhf,'(''gs ba cl '',2g14.5,
     &                              '' '',a3,'' cl np'')')
     &                              xyp(1,i),xyp(2,i),sdif(imar)

                  else if( icw .eq. 0 ) then

                        write(jhf,'(''gs bx cl np'')')

                  end if

                        write(jhf,'(3g14.5,'' dx'')')
     &                              xyp(3,i),xyp(1,i),xyp(2,i)

                     if( icw .eq. 2 ) then

                        write(jhf,'(3g14.5,'' dy st'')')
     &                  xyp(2,i) + syl, xyp(2,i) - syl, xyp(3,i)

                     else

                        write(jhf,'(3g14.5,'' dy st gr np'')')
     &                  xyp(2,i) + syl, xyp(2,i) - syl, xyp(3,i)

                     end if


                        xp1 = xyp(3,i)
                        yp1 = xyp(2,i) + syl
                        yp2 = xyp(2,i) - syl

                     if( icw .eq. 0 .or. icw .eq. 1 ) then

                        xpm = max(0.0d0,xp1)
                        xp1 = min(1.0d0,xpm)

                     end if

                        call bbox(3,0,xp1,yp1,0.d0)
                        call bbox(3,0,xp1,yp2,0.d0)


               end if

            end if

*-----------------------------------------------------------------------

                  indd = 0

            if( iydc .eq. 3 ) then

               if( ( ( xyp(6,i) .ge. 0.0 .and.
     &                 xyp(5,i) .le. 1.0 .and.
     &                 xyp(1,i) .ge. 0.0 .and. xyp(1,i) .le. 1.0 .and.
     &               ( icw .eq. 0 .or. icw .eq. 1 ) ) .or.
     &               ( icw .eq. 2 .or. icw .eq. 3 ) ) .and.
     &               ( xyp(5,i) .lt. xyp(2,i) - sym .or.
     &                 xyp(6,i) .gt. xyp(2,i) + syp ) ) then

                  if( idbg .eq. 1 )  write(jhf,'()')

                  if( icw .eq. 1 ) then

                        write(jhf,'(''gs by cl '',2g14.5,
     &                              '' '',a3,'' cl np'')')
     &                              xyp(1,i),xyp(2,i),sdif(imar)

                  else if( icw .eq. 3 ) then

                        write(jhf,'(''gs ba cl '',2g14.5,
     &                              '' '',a3,'' cl np'')')
     &                              xyp(1,i),xyp(2,i),sdif(imar)

                  else if( icw .eq. 0 ) then

                        write(jhf,'(''gs by cl np'')')

                  end if

                  if( xyp(5,i) .lt. xyp(2,i) - sym .and.
     &                xyp(6,i) .gt. xyp(2,i) + syp ) then

                        write(jhf,'(3g14.5,'' dy'')')
     &                              xyp(5,i), xyp(6,i), xyp(1,i)
                        write(jhf,'(3g14.5,'' dx'')')
     &                  xyp(1,i) + sxl, xyp(1,i) - sxl, xyp(5,i)

                     if( icw .eq. 2 ) then

                        write(jhf,'(3g14.5,'' dx st'')')
     &                  xyp(1,i) + sxl, xyp(1,i) - sxl, xyp(6,i)

                     else

                        write(jhf,'(3g14.5,'' dx st gr np'')')
     &                  xyp(1,i) + sxl, xyp(1,i) - sxl, xyp(6,i)

                     end if


                        xp1 = xyp(1,i) + syl
                        xp2 = xyp(1,i) - syl
                        yp1 = xyp(5,i)
                        yp2 = xyp(6,i)

                     if( icw .eq. 0 .or. icw .eq. 1 ) then

                        ypm = max(0.0d0,yp1)
                        yp1 = min(1.0d0,ypm)
                        ypm = max(0.0d0,yp2)
                        yp2 = min(1.0d0,ypm)

                     end if

                        call bbox(3,0,xp1,yp1,0.d0)
                        call bbox(3,0,xp1,yp2,0.d0)
                        call bbox(3,0,xp2,yp1,0.d0)
                        call bbox(3,0,xp2,yp2,0.d0)


                  else if( xyp(5,i) .lt. xyp(2,i) - sym ) then

                        write(jhf,'(3g14.5,'' dy'')')
     &                              xyp(5,i), xyp(2,i), xyp(1,i)

                     if( icw .eq. 2 ) then

                        write(jhf,'(3g14.5,'' dx st'')')
     &                  xyp(1,i) + sxl, xyp(1,i) - sxl, xyp(5,i)

                     else

                        write(jhf,'(3g14.5,'' dx st gr np'')')
     &                  xyp(1,i) + sxl, xyp(1,i) - sxl, xyp(5,i)

                     end if


                        xp1 = xyp(1,i) + syl
                        xp2 = xyp(1,i) - syl
                        yp1 = xyp(5,i)

                     if( icw .eq. 0 .or. icw .eq. 1 ) then

                        ypm = max(0.0d0,yp1)
                        yp1 = min(1.0d0,ypm)

                     end if

                        call bbox(3,0,xp1,yp1,0.d0)
                        call bbox(3,0,xp2,yp1,0.d0)


                  else if( xyp(6,i) .gt. xyp(2,i) + syp ) then

                        write(jhf,'(3g14.5,'' dy'')')
     &                              xyp(2,i), xyp(6,i), xyp(1,i)

                     if( icw .eq. 2 ) then

                        write(jhf,'(3g14.5,'' dx st'')')
     &                  xyp(1,i) + sxl, xyp(1,i) - sxl, xyp(6,i)

                     else

                        write(jhf,'(3g14.5,'' dx st gr np'')')
     &                  xyp(1,i) + sxl, xyp(1,i) - sxl, xyp(6,i)

                     end if


                        xp1 = xyp(1,i) + syl
                        xp2 = xyp(1,i) - syl
                        yp2 = xyp(6,i)

                     if( icw .eq. 0 .or. icw .eq. 1 ) then

                        ypm = max(0.0d0,yp2)
                        yp2 = min(1.0d0,ypm)

                     end if

                        call bbox(3,0,xp1,yp2,0.d0)
                        call bbox(3,0,xp2,yp2,0.d0)


                  end if

               end if

            else if( iydc .eq. 2 ) then

               if( ( ( xyp(6,i) .ge. 0.0 .and.
     &                 xyp(2,i) .le. 1.0 .and.
     &                 xyp(1,i) .ge. 0.0 .and. xyp(1,i) .le. 1.0 .and.
     &               ( icw .eq. 0 .or. icw .eq. 1 ) ) .or.
     &               ( icw .eq. 2 .or. icw .eq. 3 ) ) .and.
     &                 xyp(6,i) .gt. xyp(2,i) + syp ) then

                  if( idbg .eq. 1 )  write(jhf,'()')

                  if( icw .eq. 1 ) then

                        write(jhf,'(''gs by cl '',2g14.5,
     &                              '' '',a3,'' cl np'')')
     &                              xyp(1,i),xyp(2,i),sdif(imar)

                  else if( icw .eq. 3 ) then

                        write(jhf,'(''gs ba cl '',2g14.5,
     &                              '' '',a3,'' cl np'')')
     &                              xyp(1,i),xyp(2,i),sdif(imar)

                  else if( icw .eq. 0 ) then

                        write(jhf,'(''gs by cl np'')')

                  end if

                        write(jhf,'(3g14.5,'' dy'')')
     &                              xyp(2,i), xyp(6,i), xyp(1,i)

                     if( icw .eq. 2 ) then

                        write(jhf,'(3g14.5,'' dx st'')')
     &                  xyp(1,i) + sxl, xyp(1,i) - sxl, xyp(6,i)

                     else

                        write(jhf,'(3g14.5,'' dx st gr np'')')
     &                  xyp(1,i) + sxl, xyp(1,i) - sxl, xyp(6,i)

                     end if


                        xp1 = xyp(1,i) + syl
                        xp2 = xyp(1,i) - syl
                        yp2 = xyp(6,i)

                     if( icw .eq. 0 .or. icw .eq. 1 ) then

                        ypm = max(0.0d0,yp2)
                        yp2 = min(1.0d0,ypm)

                     end if

                        call bbox(3,0,xp1,yp2,0.d0)
                        call bbox(3,0,xp2,yp2,0.d0)


               end if

            else if( iydc .eq. 1 ) then

               if( ( ( xyp(2,i) .ge. 0.0 .and.
     &                 xyp(5,i) .le. 1.0 .and.
     &                 xyp(1,i) .ge. 0.0 .and. xyp(1,i) .le. 1.0 .and.
     &               ( icw .eq. 0 .or. icw .eq. 1 ) ) .or.
     &               ( icw .eq. 2 .or. icw .eq. 3 ) ) .and.
     &                 xyp(5,i) .lt. xyp(2,i) - sym ) then

                  if( idbg .eq. 1 )  write(jhf,'()')

                  if( icw .eq. 1 ) then

                        write(jhf,'(''gs by cl '',2g14.5,
     &                              '' '',a3,'' cl np'')')
     &                              xyp(1,i),xyp(2,i),sdif(imar)

                  else if( icw .eq. 3 ) then

                        write(jhf,'(''gs ba cl '',2g14.5,
     &                              '' '',a3,'' cl np'')')
     &                              xyp(1,i),xyp(2,i),sdif(imar)

                  else if( icw .eq. 0 ) then

                        write(jhf,'(''gs by cl np'')')

                  end if

                        write(jhf,'(3g14.5,'' dy'')')
     &                              xyp(5,i), xyp(2,i), xyp(1,i)

                        write(jhf,'(3g14.5,'' dy'')')
     &                              xyp(5,i), xyp(2,i), xyp(1,i)

                     if( icw .eq. 2 ) then

                        write(jhf,'(3g14.5,'' dx st'')')
     &                  xyp(1,i) + sxl, xyp(1,i) - sxl, xyp(5,i)

                     else

                        write(jhf,'(3g14.5,'' dx st gr np'')')
     &                  xyp(1,i) + sxl, xyp(1,i) - sxl, xyp(5,i)

                     end if


                        xp1 = xyp(1,i) + syl
                        xp2 = xyp(1,i) - syl
                        yp1 = xyp(5,i)

                     if( icw .eq. 0 .or. icw .eq. 1 ) then

                        ypm = max(0.0d0,yp1)
                        yp1 = min(1.0d0,ypm)

                     end if

                        call bbox(3,0,xp1,yp1,0.d0)
                        call bbox(3,0,xp2,yp1,0.d0)


               end if

            end if

  100    continue


*-----------------------------------------------------------------------

            if( idbg .eq. 1 )  write(jhf,'()')

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      subroutine wsymb01(icw,jhf,ild,ixdc,iydc,rleng,idbg,
     &                   ityp,iwi,ipol,imar,rcol,rcob,icoma,facsiz,
     &                   xyp,xfac,sybw)
*                                                                      *
*                                                                      *
*        PURPOSE    :  WRITE SYMBOLS ON jhf                            *
*                                                                      *
*              ICW  : 0 -> DRAW SYMBOL INSIDE THE AXIS                 *
*                   : 1 -> CLIP INSIDE THE AXIS AND SYMBOLS            *
*                   : 2 -> DRAW SYMBOLS WITHOUT CLIP                   *
*                   : 3 -> CLIP INSIDE THE SYMBOLS                     *
*                                                                      *
************************************************************************

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'angel01.inc'

*-----------------------------------------------------------------------

      parameter (ddcm = 8.4666667d-3)

      parameter ( ss0 = 0.254 )
      parameter ( il0 = 4 )

*-----------------------------------------------------------------------

      dimension xyp(6,*)

      dimension sfac(22), sfxm(22), sfxp(22), sfym(22), sfyp(22)
      character*3 sdif(22)

      common /frm/  xal, yal

      data sfac/ 0.5,  1.0, 1.0, -1.0,  0.9,    -0.9,    1.2, -1.2,
     &           1.2, -1.2, 0.9, -0.9,  0.6364, -0.6364, 0.7,  2.5,
     &           1.0,  0.9, 1.2,  1.2,  0.9,     0.6364/

      data sfxm/ 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
     &           0.5, 0.5, 0.5, 0.5, 1.0, 1.0, 0.5, 0.2,
     &           0.5, 0.5, 0.5, 0.5, 0.5, 1.0/

      data sfxp/ 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
     &           0.5, 0.5, 0.5, 0.5, 1.0, 1.0, 0.5, 0.2,
     &           0.5, 0.5, 0.5, 0.5, 0.5, 1.0/

      data sfym/ 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.288675, 0.288675,
     &           0.57735, 0.57735, 0.866025, 0.866025,
     &           1.0, 1.0, 0.5, 0.2,
     &           0.5, 0.5, 0.288675, 0.57735, 0.866015, 1.0/

      data sfyp/ 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.57735, 0.57735,
     &           0.288675, 0.288675, 0.866025, 0.866025,
     &           1.0, 1.0, 0.5, 0.2,
     &           0.5, 0.5, 0.288675, 0.57735, 0.866015, 1.0/

      data sdif/ 'cir','pls','cir','cir','squ','squ','tr1','tr1',
     &           'tr2','tr2','di1','di1','di2','di2','crs','ast',
     &           'cir','squ','tr1','tr2','di1','di2'/

      dimension rcol(3), rcob(3)

*-----------------------------------------------------------------------

               if( icw .eq. 1 ) then

                  write(jhf,'(''gs ax cl'')')

               else if( icw .eq. 3 ) then

                  write(jhf,'(''gs ba cl'')')

               end if

*-----------------------------------------------------------------------

                  il = il0

         if( sybw .gt. -r0max ) then

            if( sybw .lt. 0.0 ) then

               if( facsiz .gt. 1.0 ) then

                  il = nint( facsiz * dble( il0 ) ) - 1

               else if( facsiz .le. 0.51 ) then

                  il = 1

               else if( facsiz .le. 0.61 ) then

                  il = 2

               else if( facsiz .le. 0.81 ) then

                  il = 3

               end if

            else if( sybw .gt. 0.0 ) then

                  il = nint( sybw )

            end if

         end if

*-----------------------------------------------------------------------

                     ss = ss0 * facsiz

            if( imar .gt. 0 ) then

               if( sfac(imar) .gt. 0.0 ) then

                     ss = ss * sfac(imar)

               else if( sfac(imar) .lt. 0.0 ) then

                  if( imar .eq.  8 .or.
     &                imar .eq. 10 ) then

                     ss = - ss * sfac(imar) + dble(il) * ddcm
     &                         * ( 1.0 / sqrt(3.0) + 0.5 )

                  else if( imar .eq. 12 ) then

                     ss = - ss * sfac(imar) + dble(il) * ddcm
     &                         * 2.0 / sqrt(3.0)

                  else if( imar .eq. 14 ) then

                     ss = - ss * sfac(imar) + dble(il) * ddcm
     &                         / sqrt(2.0)

                  else

                     ss = - ss * sfac(imar) + dble(il) * ddcm

                  end if

               end if

            end if

*-----------------------------------------------------------------------

            if( imar .eq.  3 .or. imar .eq.  5 .or.
     &          imar .eq.  7 .or. imar .eq.  9 .or.
     &          imar .eq. 11 .or. imar .eq. 13 ) then


               write(jhf,'(3f7.3,'' sc /sm '',e14.6,'' cm N sd0 '',
     &                     i2,'' dd lw'')') rcob, ss, il

            else if( imar .eq. 16 ) then

                  write(jhf,'(''/sm '',e14.6,'' cm N smc'')') ss

            else

               write(jhf,'(3f7.3,'' sc /sm '',e14.6,'' cm N sd0 '',
     &                     i2,'' dd lw'')') rcol, ss, il

            end if

*-----------------------------------------------------------------------

            do 100 i = 1, ild

               if( ( ( icw .eq. 0 .or. icw .eq. 1 ) .and.
     &             xyp(1,i) .ge. 0.0 .and. xyp(1,i) .le. 1.0 .and.
     &             xyp(2,i) .ge. 0.0 .and. xyp(2,i) .le. 1.0 ) .or.
     &             icw .eq. 2 .or. icw .eq. 3 ) then

                  if( icw .eq. 0 .or. icw .eq. 2 ) then


                     if( imar .eq.  3 .or. imar .eq.  5 .or.
     &                   imar .eq.  7 .or. imar .eq.  9 .or.
     &                   imar .eq. 11 .or. imar .eq. 13 ) then


                        write(jhf,'(2g14.5,'' '',a3,'' gs fl gr '',/
     &                              3f7.3,'' sc st '',3f7.3,'' sc'')')
     &                  xyp(1,i),xyp(2,i),sdif(imar),rcol,rcob


                     else


                        write(jhf,'(2g14.5,'' '',a3)')
     &                              xyp(1,i),xyp(2,i),sdif(imar)

                     end if


                        xp1 = xyp(1,i) * xal - sfxm(imar) * ss
                        xp2 = xyp(1,i) * xal + sfxp(imar) * ss
                        yp1 = xyp(2,i) * yal - sfym(imar) * ss
                        yp2 = xyp(2,i) * yal + sfyp(imar) * ss

                        call bbox(2,0,xp1,yp1,0.d0)
                        call bbox(2,0,xp1,yp2,0.d0)
                        call bbox(2,0,xp2,yp1,0.d0)
                        call bbox(2,0,xp2,yp2,0.d0)


                  else if( icw .eq. 1 .or. icw .eq. 3 ) then

                        write(jhf,'(2g14.5,'' '',a3,'' cl'')')
     &                              xyp(1,i),xyp(2,i),sdif(imar)

                  end if

               end if

  100       continue

*-----------------------------------------------------------------------

            if( icw .eq. 0 .or. icw .eq. 2 ) then

                  if( imar .eq.  2 .or.
     &                imar .eq. 17 .or.
     &                imar .eq. 18 .or.
     &                imar .eq. 19 .or.
     &                imar .eq. 20 .or.
     &                imar .eq. 21 .or.
     &                imar .eq. 22 .or.
     &                imar .eq. 15 ) then

                     write(jhf,'(''st'')')

                  else if( imar .le. 16 ) then

                     write(jhf,'(''fl'')')

                  end if

            else if( icw .eq. 1 .or. icw .eq. 3 ) then

                  write(jhf,'(''np'')')

            end if

*-----------------------------------------------------------------------

            if( idbg .eq. 1 )  write(jhf,'()')

*-----------------------------------------------------------------------

      return
      end

