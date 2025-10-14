************************************************************************
*                                                                      *
      subroutine setrisrc(j33,nk0,nk1,ierr)
*                                                                      *
*       read decay data file based on DECDC2                           *
*             + Index file: RIsource.dat                               *
*             + Radiation data file: RIsource.rad                      *
*             + Beta spectra file: RIsource.bet                        *
*             + Auger-CK electron spectra file: RIsource.ack           *
*       , construct decay chain (linear chain)                         *
*       , calculate the activity by using the Bateman solution         *
*       and set in [Source] the energy and emission rate               *
*       made by N.Matsuda on 2016/07/31 for Gamma                      *
*       revised by N.Matsuda on 2017/05/29 for Alpha, Beta and Gamma   *
*                                                                      *
*       input :                                                        *
*          chalct    nuclide  ex) 55137.0 ( = Cs-137 )                 *
*          act000    initial activity (Bq)                             *
*          activy    activity (atoms) after cooling time               *
*          decayt    cooling time (sec) from reference date            *
*       output :                                                       *
*          ierr      error flag                                        *
*                                                                      *
************************************************************************
      use moddas
      use moddas_source

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'risrcparam.inc'
      include 'err.inc'
*     risrcparam.inc: include parameters
*     integer nuclinmax, nuclpumax, maxchain, chainmax
*     parameter ( nuclinmax =    100 ) nuclides including daughter
*     parameter ( nuclpumax =  1,000 ) Not used
*     parameter ( maxchain  =     23 ) from Number of chain (.NDX)
*     parameter ( chainmax  =  4,050 ) from Number of linear chain (.NDX)
*
*     integer rimax1, rimax2
*     parameter ( rimax1    =  4,000 ) from Number of auger electron (.ACK)
*     parameter ( rimax2    = 20,000 ) rimax1 * 5

*-----------------------------------------------------------------------

      common /error/  m_err, l_err, k_err
      character       m_err*200

*     for PHITS ( e-type = 8, or 9 ) -> ( = 22, or 23 )  2017.05.29
      common /isorst/ jstyp(isrc), istyp(isrc), inkf0(isrc), lstyp(isrc)
!$OMP THREADPRIVATE(/isorst/)
      common /ptname/ pname(20), ipln(20)
      character       pname*8

      common /isorse/ ngrp(isrc), ngei(isrc), ngea(isrc), ngfe(isrc),
     &                ngft(isrc), ngll(isrc), ngpi(isrc), ngpw(isrc)

      integer j33, nk0, nk1, nk2, nk3, nk4, nk5, nk6  ! nk5, 6 for beta
*     NEW parameters ===
      common /risource/ chalct(nuclinmax,isrc), act000(nuclinmax,isrc),
     &                  activy(nuclinmax,isrc), thf(nuclinmax,isrc),
     &                  decayt(isrc)
      real  chalct  ! ex) 55137.0 ( = Cs-137 )
      double precision  act000, activy, thf, decayt

C MATSUDA 2018.08.15 (icharacterx)
C MATSUDA 2019.05.08 (iannih)
      common /risrc00/ niorg(isrc), nicur(isrc), norm(isrc),
     &                 iaugers(isrc), icharacterx(isrc), iannih(isrc),
     &                 aclow(isrc)
      integer  niorg, nicur, norm, iaugers, icharacterx, iannih
      double precision  aclow
      common /risrc01/ normfact(isrc), asfsum(isrc)
      double precision  normfact, asfsum

*     IO number for this subroutine ===
      integer  iot, iodecdc2, ioddc2r, ioddc2b, ioddc2a, ioddc2s
      integer  ierr

      common /paran/  icfn(100), ilfn(100), chfn(100)
      character chfn*200
*     file name for READ
      character fname*200 ! T.Sato 2022/11/30, increase from 100 to 200 because inumc requests 200 character
*     file type (0: ORIGINAL of DCHAIN format, 1: DECDC2 format)
      integer  ndtype

*     INTERNAL variables
      real  nucld(4)
      real, allocatable :: chalctd(:,:)  ! chalctd(4,nuclinmax)
      double precision  thfa, bratio(4)
      double precision, allocatable :: br(:,:), r(:), nuclo(:)
*                       br(4,nuclinmax), r(nuclinmax), nuclo(nuclinmax)

      character  chaa*7
      real  chza
      logical  lex
      integer  nsf  ! Number of spontaneous fission

      character rcha*1
      double precision  asum00, actsum, abqsum, a0, rnum, timehf, anorm,
     &                  ratsum, kariv

      integer nm, kkk
      double precision, allocatable :: eng(:), rat(:)
      double precision, allocatable :: engall(:), ratall(:)
*     parameter ( rimax1    =  4,000 )
*     parameter ( rimax2    = 20,000 )
*     eng(rimax1), rat(rimax1), engall(rimax2), ratall(rimax2)

      character chin*200

C     for BETA -/+/both (.BET)  MATSUDA 2017.05.29
      parameter ( nddc2b_minu = 480 )  ! beta minus
      parameter ( nddc2b_plus = 429 )  ! beta plus
      parameter ( nddc2b_both =  22 )  ! both
      common /ddc2beta0/ tbetas00(3,nddc2b_minu)  ! from (.BET)
      parameter ( nddc2bbemax = 120 )  ! Maximum number of energy mesh
*     parameter ( nddc2b_both =  22 )  ! both
      common /ddc2beta3/ ddc2beng(nddc2bbemax), ddc2bend(2*nddc2b_both),
     &                   ddc2byld(2*nddc2b_both,nddc2bbemax),
     &                   ddc2bnuc(nddc2b_both)
      common /ddc2beta4/ nnddc2bn(2*nddc2b_both)

      character ddc2bnuc*7

      double precision smallv
      data smallv / 1.0d-36 /

*-----------------------------------------------------------------------
*     file number of scratch file for CHECKING decay calculation
           iot      = 820
*     file number of RIsource.dat or ICRP-07.NDX based on DECDC2
           iodecdc2 = 811
*     file number of RIsource.rad (radiation data) based on DECDC2
           ioddc2r  = 812
*     file number of RIsource.bet (beta spectra) based on DECDC2
           ioddc2b  = 813
*     file number of RIsource.ack (Auger-CK spectra) based on DECDC2
           ioddc2a  = 814
*     file number of RIsource.dat (decay data) formatted with DCHAIN
           ioddc2s  = 815

c     err flag for this subroutine
           ierr     =   0
*-----------------------------------------------------------------------
      allocate(chalctd(4,nuclinmax))
      allocate(br(4,nuclinmax), r(nuclinmax), nuclo(nuclinmax))
      allocate(eng(rimax1), rat(rimax1))
      rat = 0.d0
      allocate(engall(rimax2), ratall(rimax2))
*-----------------------------------------------------------------------
*     DATA CHECK for this subroutine (scratch: default)
*-----------------------------------------------------------------------
*     SCRATCH
        open( iot, status = 'scratch' )

*     OUTPUT
*     inquire( file = 'out_RIsource.txt', exist = lex )
*     if( .not. lex ) then
*       open( iot, status = 'replace', form = 'formatted',
*    &         file = 'out_RIsource.txt' )
*         rewind iot
*     else
*       open( iot, status = 'old',     form = 'formatted',
*    &         file = 'out_RIsource.txt', position = 'append' )
*     end if

        write(iot,'(  "# ",78a)') ( '=', i = 1, 78 )

*-----------------------------------------------------------------------
*     subroutine: set RI source Start
*-----------------------------------------------------------------------
  200 continue
      write(iot,'("SUBROUTINE START!! --- Set RI source ---")')

*     dummy entry of RI source (chalct and act000) for check  ###
*     ====================
*     j33 = 1
*     nk0 = 1
*     chalct(1,j33) = 55134.0  ! 'Cs-134 '
*     chalct(2,j33) = 55137.0  ! 'Cs-137 '
*     act000(1,j33) =   100.0  ! (Bq)
*     act000(2,j33) =    20.0  ! (Bq)
C ex) Pb-214 : 82214.0, At-217 : 85217.0, Ac-225 : 89225.0,
C     Cm-250 : 96250.0 <- Spontaneous fission (SF)
*                                                   ====================
        niorg(j33) = nk0

      write(iot,'("## RI nuclide entry")') ! CHECK the entry data
      do n = 1, nk0
        call listoutw(iot,j33,1,n,ierr)
      end do

*-----------------------------------------------------------------------
*     address check of decay data file (Decay dara lib. or ICRP-07.NDX)
*-----------------------------------------------------------------------
      write(iot,'("## address check of Decay Data File (DECDC2)")')

*     dummy entry of data file address (chfn(24))  ###
*     ====================
*     chfn(24) = 'C:\phits\data\RIsource.dat        '
*                                                   ====================

        write(iot,'(2x,"file(24) = ",a)') chfn(24)

        call fchkd2(24,ndtype,fname,ierr)

          if( ierr .eq. 1 ) then
*           l_err = ill(jsn)
*           k_err = jsn
            go to 998
          end if

        icl = inumc(fname,1,100,' ') - 1
*       do n = 100, 1, -1
*         if( fname(n:n) .ne. ' ' ) go to 201
*       end do
*
* 201     icl = n

        write(iot,'(2x,"Address:   ",a)') fname(1:icl)

*     CHECK THE ADDRESS OF THE OTHER EXTRA DATA FILE  2017.05.29
        write(iot,'(2x,"Extra data:")')
        lex = .false.

      if( inkf0(j33) .eq. 22 ) then  ! GAMMA-RAY
        if( ndtype .eq.  0 ) then
          fname(icl-3:icl) = '.rad'
          write(iot,'(2x,"fname = ",a)') fname
          inquire( file = fname, exist = lex )
          if( lex ) then
*           O.K. ... RIsource.rad exists.
                write(iot,'(13x,"RIsource.rad for gamma-rays")')
          else
            if( icl .gt. 12 ) then
              inquire( file = fname(1:icl-12)//'ICRP-07.RAD ',
     &               exist = lex )
              if( lex ) then
                write(iot,'(13x,"(ICRP-07.RAD) for gamma-rays")')
              else
*           O.K. ... RIsource.dat exists.
*               go to 946  ! but, RIsource.rad and ICRP-07.RAD does NOT exist.
                if( icharacterx(j33) .eq. 0 .or.
     &              icharacterx(j33) .eq. 2 ) go to 946
                if( iannih(j33) .eq. 1 ) go to 946
                fname(icl-3:icl) = '.dat'
              end if
            else if( icl .eq. 12 ) then
              inquire( file = 'ICRP-07.RAD', exist = lex )
              if( lex ) then
                write(iot,'(13x,"(ICRP-07.RAD) for gamma-rays")')
              else
*           O.K. ... RIsource.dat exists.
*               go to 946  ! but, RIsource.rad and ICRP-07.RAD does NOT exist.
                if( icharacterx(j33) .eq. 0 .or.
     &              icharacterx(j33) .eq. 2 ) go to 946
                if( iannih(j33) .eq. 1 ) go to 946
                fname(icl-3:icl) = '.dat'
              end if
            end if
          end if
        else
          fname(icl-2:icl) = 'RAD'
          inquire( file = fname, exist = lex )
          if( .not. lex ) go to 956  ! ICRP-07.RAD does NOT exist.
                write(iot,'(13x,"ICRP-07.RAD for gamma-rays")')
        end if

      else if( inkf0(j33) .eq. 2000004 ) then  ! ALPHA-PARTICLE
        if( ndtype .eq.  0 ) then
          fname(icl-2:icl) = 'rad'
          write(iot,'(2x,"fname = ",a)') fname
          inquire( file = fname, exist = lex )
          if( lex ) then
*           O.K. ... RIsource.rad exists.
                write(iot,'(13x,"RIsource.rad for alpha-particles")')
          else
            if( icl .gt. 12 ) then
              inquire( file = fname(1:icl-12)//'ICRP-07.RAD ',
     &                 exist = lex )
              if( lex ) then
                write(iot,'(13x,"(ICRP-07.RAD) for alpha-particles")')
              else
                go to 942  ! RIsource.rad does NOT exist.
              end if
            else if( icl .eq. 12 ) then
              inquire( file = 'ICRP-07.RAD', exist = lex )
              if( lex ) then
                write(iot,'(13x,"(ICRP-07.RAD) for alpha-particles")')
              else
                go to 942  ! RIsource.rad does NOT exist.
              end if
            end if
          end if
        else
          fname(icl-2:icl)  = 'RAD'
          inquire( file = fname, exist = lex )
          if( .not. lex ) go to 952  ! ICRP-07.RAD does NOT exist.
                write(iot,'(13x,"ICRP-07.RAD for alpha-particles")')
        end if

      else if( inkf0(j33) .eq. 11 ) then  ! ELECTRON
        if( iaugers(j33) .eq. 2 ) go to 206
        if( ndtype .eq.  0 ) then
          fname(icl-2:icl) = 'bet'
          write(iot,'(2x,"fname = ",a)') fname
          inquire( file = fname, exist = lex )
          if( lex ) then
*           O.K. ... RIsource.bet exists.
                write(iot,'(13x,"RIsource.bet for electrons")')
          else
            if( icl .gt. 12 ) then
              inquire( file = fname(1:icl-12)//'ICRP-07.BET ',
     &                 exist = lex )
              if( lex ) then
                write(iot,'(13x,"(ICRP-07.BET) for electrons")')
              else
                go to 943  ! RIsource.bet does NOT exist.
              end if
            else if( icl .eq. 12 ) then
              inquire( file = 'ICRP-07.BET', exist = lex )
              if( lex ) then
                write(iot,'(13x,"(ICRP-07.BET) for electrons")')
              else
                go to 943  ! RIsource.bet does NOT exist.
              end if
            end if
          end if
        else
          fname(icl-2:icl)  = 'BET'
          inquire( file = fname, exist = lex )
          if( .not. lex ) go to 953  ! ICRP-07.BET does NOT exist.
*           O.K. ... ICRP-07.BET exists.
                write(iot,'(13x,"ICRP-07.BET for electrons")')
        end if

  206 continue
        if( iaugers(j33) .eq. 1 ) go to 207
        if( ndtype .eq.  0 ) then
          fname(icl-2:icl)  = 'ack'
          write(iot,'(2x,"fname = ",a)') fname
          inquire( file = fname, exist = lex )
          if( lex ) then
*           O.K. ... RIsource.ack exists.
                write(iot,'(13x,"RIsource.ack for auger electrons")')
          else
            if( icl .gt. 12 ) then
              inquire( file = fname(1:icl-12)//'ICRP-07.ACK ',
     &                 exist = lex )
              if( lex ) then
                write(iot,'(13x,"(ICRP-07.ACK) for auger electrons")')
              else
                go to 944  ! RIsource.ack does NOT exist.
              end if
            else if( icl .eq. 12 ) then
              inquire( file = 'ICRP-07.ACK', exist = lex )
              if( lex ) then
                write(iot,'(13x,"(ICRP-07.ACK) for auger electrons")')
              else
                go to 944  ! RIsource.ack does NOT exist.
              end if
            end if
          end if
        else
          fname(icl-2:icl)  = 'ACK'
          inquire( file = fname, exist = lex )
          if( .not. lex ) go to 954  ! ICRP-07.ACK does NOT exist.
                write(iot,'(13x,"ICRP-07.ACK for auger electrons")')
        end if

C       for internal conversion electron  ! MATSUDA 2018.08.15
        if( ndtype .eq.  0 ) then
          fname(icl-2:icl) = 'rad'
          write(iot,'(2x,"fname = ",a)') fname
          inquire( file = fname, exist = lex )
          if( lex ) then
*           O.K. ... RIsource.rad exists.
                write(iot,'(13x,"RIsource.rad for IE")')
          else
            if( icl .gt. 12 ) then
              inquire( file = fname(1:icl-12)//'ICRP-07.RAD ',
     &                 exist = lex )
              if( lex ) then
                write(iot,'(13x,"(ICRP-07.RAD) for IE")')
              else
                go to 947  ! RIsource.rad does NOT exist.
              end if
            else if( icl .eq. 12 ) then
              inquire( file = 'ICRP-07.RAD', exist = lex )
              if( lex ) then
                write(iot,'(13x,"(ICRP-07.RAD) for IE")')
              else
                go to 947  ! RIsource.rad does NOT exist.
              end if
            end if
          end if
        else
          fname(icl-2:icl)  = 'RAD'
          inquire( file = fname, exist = lex )
          if( .not. lex ) go to 957  ! ICRP-07.RAD does NOT exist.
*           O.K. ... ICRP-07.RAD exists.
                write(iot,'(13x,"ICRP-07.RAD for IE")')
        end if
  207 continue

      else if( inkf0(j33) .eq. -11 ) then  ! POSITRON
        if( ndtype .eq.  0 ) then
          fname(icl-2:icl) = 'bet'
          write(iot,'(2x,"fname = ",a)') fname
          inquire( file = fname, exist = lex )
          if( lex ) then
*           O.K. ... RIsource.bet exists.
                write(iot,'(13x,"RIsource.ack for positrons")')
          else
            if( icl .gt. 12 ) then
              inquire( file = fname(1:icl-12)//'ICRP-07.BET ',
     &                 exist = lex )
              if( lex ) then
                write(iot,'(13x,"(ICRP-07.BET) for positrons")')
              else
                go to 945  ! RIsource.bet does NOT exist.
              end if
            else if( icl .eq. 12 ) then
              inquire( file = 'ICRP-07.BET', exist = lex )
              if( lex ) then
                write(iot,'(13x,"(ICRP-07.BET) for positrons")')
              else
                go to 945  ! RIsource.bet does NOT exist.
              end if
            end if
          end if
        else
          fname(icl-2:icl)  = 'BET'
          inquire( file = fname, exist = lex )
          if( .not. lex ) go to 955  ! ICRP-07.BET does NOT exist.
                write(iot,'(13x,"ICRP-07.RAD for positrons")')
        end if

      else if( inkf0(j33) .eq. 2112 ) then
        if( ndtype .eq. 0 ) then
          fname(icl-2:icl) = 'nsf'
          write(iot,'(2x,"fname = ",a)') fname
          inquire( file = fname, exist = lex )
          if( lex ) then
*           O.K. ... RIsource.nsf exists.
                write(iot,'(13x,"RIsource.nsf for SF neutron")')
          else
            if( icl .gt. 12 ) then
              inquire( file = fname(1:icl-12)//'ICRP-07.NSF ',
     &                 exist = lex )
              if( lex ) then
                write(iot,'(13x,"(ICRP-07.NSF) for SF neutron")')
              else
                go to 958  ! RIsource.nsf does NOT exist.
              end if
            else if( icl .eq. 12 ) then
              inquire( file = 'ICRP-07.NSF', exist = lex )
              if( lex ) then
                write(iot,'(13x,"(ICRP-07.NSF) for positrons")')
              else
                go to 958  ! RIsource.nsf does NOT exist.
              end if
            end if
          end if
        else
          fname(icl-2:icl)  = 'NSF'
          inquire( file = fname, exist = lex )
          if( .not. lex ) go to 958  ! ICRP-07.NSF does NOT exist.
          write(iot,'(13x,"ICRP-07.NSF for spontaeneous fission")')
        end if

      else

              go to 969

      end if


        if( ndtype .eq.  0 ) then
          fname(icl-2:icl) = 'dat'
        else
          fname(icl-2:icl)  = 'NDX'
        end if

      open( iodecdc2, file = fname(1:icl),
     &      form = 'formatted', status = 'unknown' )
        rewind iodecdc2

*-----------------------------------------------------------------------
*     read and search decay data file (Decay data lib.(0) of DECDC2(1))
*-----------------------------------------------------------------------
  210 continue

          nm = 0
        call rddecdc2(999,nm,chza,thfa,nucld,bratio,ierr)

          if( ierr .eq. 1 ) goto 998

        write(iot,'("# ",78a)') ( '=', i = 1, 78 )
        write(iot,'("# ENTRY DATA ...")')

      do n = 1, nk0
        chza = chalct(n,j33)

        call rddecdc2(ndtype,nm,chza,thfa,nucld,bratio,ierr)

          if( ierr .eq. 1 ) goto 998

          if( chza .ne. chalct(n,j33) ) go to 923
          if( thfa .eq. -1.0d+00 ) then
            call getelm(11,chaa,chza,ierr)

            if( ierr .eq. 1 ) then
*             l_err = ill(jsn)
*             k_err = jsn
              go to 998
            end if

              go to 924
          end if

*         chalct(n,j33) = chza
          thf(n,j33) = thfa
          do m = 1, 4
            chalctd(m,n) = nucld(m)
            br(m,n) = bratio(m)
          end do
*         a0(n) = act000(n,j33)  ! (Bq), act000(n,j33): initial activity
          activy(n,j33) = 0.0d+00  ! (atoms)
          nuclo(n) = n  ! Parent nuclide

      end do

*-----------------------------------------------------------------------
*     searching daughter
*-----------------------------------------------------------------------
      write(iot,'("# SEARCHING DAUGHTER ...")')
        nk1 = nk0
      if( decayt(j33) .eq. 0.0 ) go to 300  ! No cooling time
        ni  = 1  ! at least 1 nuclide
        nsf = 0
      write(iot,'("  START:",i3," ... ",f7.1)')
     &                         ni,chalct(ni,j33)

  260 continue   ! Create nuclides list
      do i = 1, 4
        if( chalctd(i,ni) .gt. 0.0 ) then
          chza = chalctd(i,ni)

          if( chza .eq. 990000.0 ) then  ! for Spontaneous fission nuclide
              nk1 = nk1 + 1  ! NEW entry (nuclide)
              nsf = nsf + 1
            chalctd(i,ni) = chza + nsf
            chalct(nk1,j33) = chza + nsf  ! = chalctd(i,ni)
            thf(nk1,j33) = -1.0d+00
*           write(iot,'(2x,''Half life:'',1pe13.5)') thfa
            do m = 1, 4
              chalctd(m,nk1) = 0.0
              br(m,nk1) =  0.0d+00
            end do
              br(1,nk1) =  br(i,ni)
            act000(nk1,j33) = 0.0d+00  ! (Bq)
            activy(nk1,j33) = 0.0d+00  ! (atoms)
            nuclo(nk1) = ni  ! Spontaneous fission nuclide

            call listoutw(iot,j33,99,ni,ierr)
            write(iot,'(4x,"Spontaneous fission (",i3," )")') nsf

          else  ! for Daughter nuclide
C           DATA CHECK
            do j = 1, nk1
              if( chza .eq. chalct(j,j33) ) then

                call listoutw(iot,j33,90,j,ierr)  ! Already registered
                  go to 280
              end if
            end do

  270     continue

            call rddecdc2(ndtype,nm,chza,thfa,nucld,bratio,ierr)

              if( ierr .eq. 1 ) goto 998

              if( chza .ne. chalctd(i,ni) ) go to 923

                nk1 = nk1 + 1  ! NEW entry (nuclide)

              chalct(nk1,j33) = chza

              thf(nk1,j33) = thfa
*             write(iot,'(2x,''Half life:'',1pe13.5)') thfa
              do m = 1, 4
                chalctd(m,nk1) = nucld(m)
                br(m,nk1) = bratio(m)
              end do

              act000(nk1,j33) = 0.0d+00  ! (Bq)
              activy(nk1,j33) = 0.0d+00  ! (atoms)
              if( ni .le. nk0 ) then  ! Parent nuclide
                nuclo(nk1) = ni
              else
                nuclo(nk1) = nuclo(ni)
              end if

              if( thfa .eq. -1.0d+00 ) then  ! Stable nuclide
                call listoutw(iot,j33,89,nk1,ierr)
              end if

  280     continue
          end if
        else
          if( chalctd(i,ni) .eq. 0.0 ) then
            if( i .eq. 1 ) write(iot,'("    pass")')
              exit  ! 2017.06.30
          else
            go to 928  ! negative number
          end if
        end if
      end do

      if( (ni+1) .gt. nk1 ) then
        write(iot,'("  END ( nuclides list )")')
      else
        ni = ni + 1
        write(iot,'("  NEXT: ",i3," ... ",f7.1)')
     &                           ni,chalct(ni,j33)
          go to 260
      end if

      close(iodecdc2)

  290 continue
        write(iot,'("# Finish to check the registered data ...")')
      do n = 1, nk1
        call listoutw(iot,j33,2,n,ierr)
      end do
        write(iot,'("* Total entry: ",i4," nuclides")') nk1

        nicur(j33) = nk1

*-----------------------------------------------------------------------
*     set initial activity and activity conversion
*         from Becquere (Bq) to number of atoms (atoms): a0
*-----------------------------------------------------------------------
  300 continue
        write(iot,'("## READ RESULT - nuclides list -")')
        asum00 = 0.0d+00  ! (atoms)

      do n = 1, nk1
        call listoutw(iot,j33,11,n,ierr)

        if( act000(n,j33) .ge. 0.0d+00 ) then
          if( thf(n,j33) .ne. -1.0d+00 ) then
c           act000(n,j33) = .......  ! (Bq) initial value
            r(n) = dble( log(2.) / thf(n,j33) )  ! decay constant
            asum00 = asum00 + act000(n,j33) / r(n)  ! (atoms)
            if( decayt(j33) .eq. 0.0d+00 ) then
              activy(n,j33) = act000(n,j33) / r(n)
*             because, decayt(j33) = 0.0d+00 -> go to 600
            end if
          else
c           act000(n,j33) = 0.0d+00  ! (Bq)
            r(n) = 0.0d+00  ! log(2.)/infini.
          end if

        end if

      end do
        write(iot,'("* Number of atoms (initial): ",1pe13.6
     &           ," (atoms)")') asum00

*-----------------------------------------------------------------------
*     construction of linear chain
*-----------------------------------------------------------------------
  400 continue
      write(iot,'("## DECAY CHAIN on nuclides list")')
      nk2 = nk1
*     dummy entry of Cooling time (decayt) for check  ###
*     ====================
*         decayt(j33) = 10.*24.*60.*60.  ! 10 days (positive entry)
*         decayt(j33) = -1.0             ! x Half life (negative entry)
*                                                   ====================
      if( decayt(j33) .eq. 0.0d+00 ) go to 600

*       Parameters for lchain
C         nk0: number of initial entry (nuclides)
C         nk1: number of registered nuclide
C         -> nk2: number of linear chain
C         chalctd: daughter
C         br: branching ratio
        call lchain(j33,nk0,nk2,chalctd,br,ierr)

          if( ierr .eq. 1 ) goto 998

*-----------------------------------------------------------------------
*     Lower limit value
*-----------------------------------------------------------------------
  600 continue
        write(iot,'("## Lower limit (original) for underflow")')
        write(iot,'("   After a cooling time: ",1pe13.6," (sec)")')
     &            decayt(j33)

            actsum = 0.0d+00  ! (atoms)
            abqsum = 0.0d+00  ! (Bq)
        do n = 1, nk1
            actsum = actsum + activy(n,j33)  ! incl. SF and stable
            abqsum = abqsum + activy(n,j33) * r(n)

          if( activy(n,j33) .lt. smallv ) then
c           underflow: data smallv / 1.0d-36 /
            activy(n,j33) = 0.0d+00  ! 1.0d-50 -> 0.0d+00
          end if

          call listoutw(iot,j33,21,n,ierr)

        end do

          write(iot,'("* Number of atoms (atoms): initial =",1pe13.6
     &             , ", after =",1pe13.6)') asum00, actsum
          write(iot,'("* Ratio: ",1pe16.9)') actsum/asum00
          write(iot,'("* Total: ",1pe13.6," (Bq)")') abqsum

          go to 610

        if( decayt(j33) .lt. 0.0d+00 ) then
          write(iot,'("  Normalize")')

*           actsum = 0.0d+00  ! (atoms)
*           abqsum = 0.0d+00  ! (Bq)

            activy(1,j33) = act000(1,j33) / r(1)
            actsum = activy(1,j33)  ! (atoms)
            abqsum = act000(1,j33)  ! (Bq)
            anorm  = act000(1,j33) / ( activy(1,j33) * r(1) )  ! (Bq/Bq)

          if( nk1 .gt. 1 ) then

          do n = 2, nk1
            activy(n,j33) = anorm * activy(n,j33)
            actsum = actsum + activy(n,j33)  ! incl. SF and stable
            abqsum = abqsum + activy(n,j33) * r(n)

            if( activy(n,j33) .lt. smallv ) then
c             underflow: data smallv / 1.0d-36 /
              activy(n,j33) = 0.0d+00  ! 1.0d-50 -> 0.0d+00
            end if

            call listoutw(iot,j33,21,n,ierr)

          end do

          end if
          write(iot,'(2x,"Ratio: ",1pe16.9)') actsum/asum00
          write(iot,'(2x,"Total: ",1pe13.6," (Bq)")') abqsum

        else
          write(iot,'("  No normalize")')

        end if

  610 continue
*       Spontaneous fission (Bq)
          asfsum(j33) = 0.0d+00

        do n = 1, nk1
          if( chalct(n,j33) .gt. 990000.0 ) then
            asfsum(j33) = asfsum(j33) +
     &                    activy(nuclo(n),j33) * r(nuclo(n)) * br(1,n)
          end if
        end do

        if ( asfsum(j33) .gt. 0.0d+00 ) then
            write(iot,'(2x,"Total (Spontaneous fission): ",1pe13.6,
     &               " (fission/sec)")') asfsum(j33)
        end if

  620 continue

*-----------------------------------------------------------------------
*     Energies and absolute yields (each particle)
*-----------------------------------------------------------------------

      icl = inumc(fname,1,100,' ') - 1
*     do n = 100, 1, -1
*       if( fname(n:n) .ne. ' ' ) go to 641
*     end do
*
* 641   icl = n

      write(iot,'(/,"## LIST UP (energy spectrum)")')

        nk3 = 0
        nk5 = rimax1
        normfact(j33) = 0.0d+00
*     parameter ( rimax1    =  4,000 ) from Number of auger electron (.ACK)
*     parameter ( rimax2    = 20,000 ) rimax1 * 5
        do k = 1, rimax1
          engall(k) = 0.0d+00
          ratall(k) = 0.0d+00
        end do

      do n = 1, nk1
        if( chalct(n,j33) .ge. 990000.0 ) then
          ! Spontaneous fission nuclide
          cycle ! frtati 2022/12/26
        else if( thf(n,j33) .eq. -1.0d+00 ) then
          ! Stable nuclide
          chza = chalct(n,j33)

          call getelm(11,chaa,chza,ierr)

            if( ierr .eq. 1 ) then
              go to 999
            end if

          write(iot,'(2x,a7,"(stable)")') chaa

        else
          ! Radioactive isotope
          chza = chalct(n,j33)

          call getelm(11,chaa,chza,ierr)

            if( ierr .eq. 1 ) then
              go to 999
            end if

          a0 = activy(n,j33) * r(n)  ! (Bq)

          if( a0 .lt. aclow(j33) ) then

            write(iot,'(2x,a7," <- under the lower limit (Bq)")') chaa
            kkk = 0

          else

            write(iot,'(2x,a7,1pe13.6," (Bq)  $ ")') chaa, a0

            if( inkf0(j33) .eq. 2000004 ) then
*             alpha ... OK
*             write(iot,'(''ALPHA'')')
              call rdalpha(ndtype,chza,a0,kkk,eng,rat,ierr)

            else if( inkf0(j33) .eq. 11 .or. inkf0(j33) .eq. -11 ) then
*             electron or positron ... OK
*             write(iot,'(''BETA - ELECTRON or POSITRON -'')')
              if( inkf0(j33) .eq. 11 .and. iaugers(j33) .eq. 2 ) then
                    kkk = 0
                    go to 652
              end if
*             do  nn = 1, nnbetas(3)   ! both electron and positron
              do  nn = 1, nddc2b_both  ! both electron and positron
*               if( tbetas00(3,nn) .eq. chza ) then
                if( abs(tbetas00(3,nn)-chza) .lt. 0.01 ) then
                  if( inkf0(j33) .eq. 11 ) then
                    kkk = nnddc2bn(nn) + 1
                  else
                    kkk = nnddc2bn(22+nn) + 1
                  end if

                    eng(1) = ddc2beng(1)
                  if( inkf0(j33) .eq. 11 ) then
                    rat(1) = ddc2byld(nn,1)
                  else
                    rat(1) = ddc2byld(22+nn,1)
                  end if

                  do nnn = 2, (kkk-1)
                    eng(nnn) = ddc2beng(nnn)
                    if( inkf0(j33) .eq. 11 ) then
                      rat(nnn) = ddc2byld(nn,nnn)
                    else
                      rat(nnn) = ddc2byld(22+nn,nnn)
                    end if
                    rat(nnn-1) = ( rat(nnn-1) + rat(nnn) ) / 2.0d+00
     &                         * ( eng(nnn) - eng(nnn-1) )
                    if( nnn .eq. 2 ) then
                      write(iot,'(2x,1pe15.7,1pe15.7,2x,
     &                          "$ Energy (MeV) and Absolute yield")')
     &                            eng(1), rat(1)

                    else
                      write(iot,'(2x,1pe15.7,1pe15.7)')
     &                                          eng(nnn-1), rat(nnn-1)

                    end if
                  end do

                  if( inkf0(j33) .eq. 11 ) then
                    eng(kkk) = ddc2bend(nn)
                  else
                    eng(kkk) = ddc2bend(22+nn)
                  end if
                    rat(kkk) = 0.0d+00
                    rat(kkk-1) = ( rat(kkk-1) + rat(kkk) ) / 2.0d+00
     &                         * ( eng(kkk) - eng(kkk-1) )
                      write(iot,'(2x,1pe15.7,1pe15.7/2x,1pe15.7)')
     &                                eng(kkk-1), rat(kkk-1), eng(kkk)

                    go to 652
                end if
              end do

              if( inkf0(j33) .eq. 11 ) then
*             electron
*               do  nn = 1, nnbetas(1)   ! electron (beta minus)
                do  nn = 1, nddc2b_minu  ! electron (beta minus)
*                 if( tbetas00(1,nn) .eq. chza ) then
                  if( abs(tbetas00(1,nn)-chza) .lt. 0.01 ) then
                    call rdbetaS(ndtype,chza,a0,kkk,eng,rat,ierr)

                      go to 652

                  end if
                end do
              else if( inkf0(j33) .eq. -11 ) then
*             positron
*               do  nn = 1, nnbetas(2)   ! positron (beta plus)
                do  nn = 1, nddc2b_plus  ! positron (beta plus)
*                 if( tbetas00(2,nn) .eq. chza ) then
                  if( abs(tbetas00(2,nn)-chza) .lt. 0.01 ) then
                    call rdbetaS(ndtype,chza,a0,kkk,eng,rat,ierr)

                      go to 652

                  end if
                end do

              end if

                    kkk = 0

  652 continue  ! for beta minus and plus

            else if( inkf0(j33) .eq. 22 ) then
*             gamma ... OK
*             write(iot,'(''GAMMA and ANNIHILATION'')')
              if( ( icharacterx(j33) .eq. 0 ) .and.
     &            ( iannih(j33) .eq. 0 ) ) then
                  kkk =  0
                call rdgamma(ndtype,chza,a0,kkk,eng,rat,ierr)
                call rdgamma2(ndtype,chza,a0,kkk,eng,rat,ierr)
                call rdgamma3(ndtype,chza,a0,kkk,eng,rat,ierr)
              else if( ( icharacterx(j33) .eq. 0 ) .and.
     &                 ( iannih(j33) .eq. 1 ) ) then
                  kkk =  0
                call rdgamma(ndtype,chza,a0,kkk,eng,rat,ierr)
                call rdgamma2(ndtype,chza,a0,kkk,eng,rat,ierr)
              else if( ( icharacterx(j33) .eq. 1 ) .and.
     &                 ( iannih(j33) .eq. 0 ) ) then
                  kkk =  0
                call rdgamma(ndtype,chza,a0,kkk,eng,rat,ierr)
                call rdgamma3(ndtype,chza,a0,kkk,eng,rat,ierr)
              else if( ( icharacterx(j33) .eq. 1 ) .and.
     &                 ( iannih(j33) .eq. 1 ) ) then
                  kkk =  0
                call rdgamma(ndtype,chza,a0,kkk,eng,rat,ierr)
              else if( icharacterx(j33) .eq. 2 ) then
                  kkk =  0
                call rdgamma2(ndtype,chza,a0,kkk,eng,rat,ierr)
              end if

            else if( inkf0(j33) .eq. 2112 ) then
*             neutron ... OK
*             write(iot,'(''NEUTRON'')')
              call rdsfneutron(ndtype,chza,a0,kkk,eng,rat,ierr)

            else
                go to 969

            end if

            if( ierr .eq. 1 ) go to 998

          end if

*     SUM of data
            ratsum = 0.0d+00  ! for this loop

          if( kkk .gt. 0 .and. a0 .gt. 0.0d+00 ) then
            if( inkf0(j33) .eq. 11 .or. inkf0(j33) .eq. -11 ) then
              if( kkk .gt. rimax1 ) go to 965
              if( kkk .gt. rimax2 ) go to 966

              if( nk3 .lt. kkk ) then

                do k = 1, kkk
                  engall(k) = eng(k)
                  ratall(k) = ratall(k) + a0*rat(k)
                  ratsum = ratsum + rat(k)
                end do
                  nk3 = kkk

              else if( nk3 .ge. kkk ) then

                do k = 1, kkk-1
*                 engall(k) = eng(k)
                  ratall(k) = ratall(k) + a0*rat(k)
                  ratsum = ratsum + rat(k)
                end do
                if( engall(nk3) .lt. eng(kkk) ) then
                  engall(nk3) = eng(kkk)
                else ! ( engall(kkk) .ge. eng(kkk) )
*                 engall(nk3) = engall(nk3)
                end if
*                 ratall(nk3) = 0.0d+00
*                 nk3 = nk3

              end if

            else  ! alpha, gamma-ray, SF-neutron
              if( kkk .gt. rimax1 ) go to 965
              if( nk3+kkk .gt. rimax2 ) go to 966

              do k = 1, kkk
C MATSUDA 2017.05.29 ( MeV -> MeV/n for alpha )
                if( inkf0(j33) .eq. 2000004 ) then
                  engall(nk3+k) = eng(k)/4.0d+00
                else  ! inkf0(j33) = 22 (gamma)
                  engall(nk3+k) = eng(k)
                end if
                  ratall(nk3+k) = a0*rat(k)
                  ratsum = ratsum + rat(k)
              end do
                  nk3 = nk3 + kkk

            end if

            if( inkf0(j33) .eq. 2000004 .or.
     &          inkf0(j33) .eq. 11 .or. inkf0(j33) .eq. -11 ) then
              write(iot,'("* Total =",1pe15.7," (particles/sec):"
     &                      ,1pe13.5," (Bq) X",1pe13.5," (%)",/)')
     &                     a0*ratsum, a0, ratsum*100
                normfact(j33) = normfact(j33) + a0 * ratsum

            else if( inkf0(j33) .eq. 22 ) then
              write(iot,'("* Total =",1pe15.7," (photon/sec):"
     &                      ,1pe13.5," (Bq) X",1pe13.5," (%)",/)')
     &                     a0*ratsum, a0, ratsum*100
              normfact(j33) = normfact(j33) + a0 * ratsum

            else if( inkf0(j33) .eq. 2112 ) then
              write(iot,'("* Total =",1pe15.7," (neutron/sec):"
     &                      ,1pe13.5," (Bq) X",1pe13.5," (%)",/)')
     &                     a0*ratsum, a0, ratsum*100
              normfact(j33) = normfact(j33) + a0 * ratsum

            else

              go to 969

            end if

          end if

    ! for Auger-CK electron spectrum
          if( a0 .gt. aclow(j33) .and. inkf0(j33) .eq. 11 ) then
*           write(iot,'(/,''## LIST UP (energy spectrum for Auger)'')')
*             electron (Auger-CK) ... OK
*             write(iot,'(''BETA - Auger-CK -'')')
            if( iaugers(j33) .eq. 1 ) then
                kkk = 0
                go to 658
            else
              call rdbetaD(ndtype,chza,a0,kkk,eng,rat,ierr)

              if( ierr .eq. 1 ) go to 998

              if( kkk .eq. 0 )
     &        call rdbetaD2(ndtype,chza,a0,kkk,eng,rat,ierr)

            end if

            if( ierr .eq. 1 ) go to 998

*     SUM of data
            ratsum = 0.0d+00  ! for this loop

            if( kkk .gt. 0 .and. a0 .gt. 0.0d+00 ) then
              if( kkk .gt. rimax1 ) go to 965
              if( nk5+kkk .gt. rimax2 ) go to 966

              do k = 1, kkk
                  engall(nk5+k) = eng(k)
                  ratall(nk5+k) = a0*rat(k)
                  ratsum = ratsum + rat(k)
              end do
                  nk5 = nk5 + kkk

              write(iot,'("* Total =",1pe15.7," (electron/sec):"
     &                      ,1pe13.5," (Bq) X",1pe13.5," (%)",/)')
     &                     a0*ratsum, a0, ratsum*100
              normfact(j33) = normfact(j33) + a0 * ratsum

            end if

          end if

    ! for internal conversion electron
          if( a0 .gt. aclow(j33) .and. inkf0(j33) .eq. 11 ) then
*           write(iot,'(/,''## LIST UP (energy spectrum for IE)'')')
*             electron (internal conversion electron) ... OK
*             write(iot,'(''BETA - internal conversion electron -'')')
            if( iaugers(j33) .eq. 1 ) then
                kkk = 0
                go to 658
            else
              call rdbetaD3(ndtype,chza,a0,kkk,eng,rat,ierr)

              if( ierr .eq. 1 ) go to 998

            end if

*     SUM of data
            ratsum = 0.0d+00  ! for this loop

            if( kkk .gt. 0 .and. a0 .gt. 0.0d+00 ) then
              if( kkk .gt. rimax1 ) go to 965
              if( nk5+kkk .gt. rimax2 ) go to 966

              do k = 1, kkk
                  engall(nk5+k) = eng(k)
                  ratall(nk5+k) = a0*rat(k)
                  ratsum = ratsum + rat(k)
              end do
                  nk5 = nk5 + kkk

              write(iot,'("* Total =",1pe15.7," (electron/sec):"
     &                      ,1pe13.5," (Bq) X",1pe13.5," (%)",/)')
     &                     a0*ratsum, a0, ratsum*100
              normfact(j33) = normfact(j33) + a0 * ratsum

            end if

          end if

        end if

  658 continue
      end do

      if( nk3 .eq. 0 .and. nk5 .eq. rimax1 ) go to 967

*-----------------------------------------------------------------------
*     Energy sorting
*-----------------------------------------------------------------------
  660 continue

        write(iot,'(/,"## Energy sorting (ascending order)")')

      if( inkf0(j33) .eq. 11 .or. inkf0(j33) .eq. -11 .or.
     &    inkf0(j33) .eq. 2112 ) then  ! for *.BET, frtati2023/04/27 and SF

          nk4 = nk3

*         write(iot,'(2x,''ne = '',i5)') nk4-1
*       do n = 1, nk4-1
*         write(iot,'(i4,''. '',1pe15.7,1pe15.7)')
*    &                 n, engall(n), ratall(n)
*       end do
*         write(iot,'(i4,''. '',1pe15.7)') nk4, engall(nk4)

          nk6 = 0
        if( nk5 .gt. rimax1 ) then
          write(iot,'(/,"## Energy sorting for Auger and IE")')
          do n = rimax1+1, nk5-1
            kariv = engall(n)
              nn = n
            do m = n, nk5
              if( engall(m) .lt. kariv ) then
                kariv = engall(m)
                  nn = m
              else
              end if
            end do

            if( n .ne. nn ) then
              kariv = engall(n)
              engall(n)  = engall(nn)
              engall(nn) = kariv
              kariv = ratall(n)
              ratall(n)  = ratall(nn)
              ratall(nn) = kariv
            end if
          end do

*         write(iot,'(2x,''ne = '',i5)') (nk5-rimax1)
*       do n = 1, nk5-rimax1
*         write(iot,'(i4,''. '',1pe15.7,1pe15.7)')
*    &                 nk4+n, engall(rimax1+n), ratall(rimax1+n)
*       end do

*         SUM UP the yield of the same energy
*           nk6 = 0
          do n = rimax1+1, nk5-1
            if( ratall(n) .gt. 0.0d+00 ) then
              if( engall(n) .eq. engall(n+1) ) then
                ratall(n+1) = ratall(n) + ratall(n+1)
                ratall(n) = 0.0d+00
              else
                nk6 = nk6 + 1
              end if
            end if
          end do

            nk6 = nk6 + 1
*           nk6 = nk6

        end if
            nk6 = rimax1 + nk6

      else
          do n = 1, nk3-1
            kariv = engall(n)
              nn = n
            do m = n, nk3
              if( engall(m) .lt. kariv ) then
                kariv = engall(m)
                  nn = m
              else
              end if
            end do

            if( n .ne. nn ) then
              kariv = engall(n)
              engall(n)  = engall(nn)
              engall(nn) = kariv
              kariv = ratall(n)
              ratall(n)  = ratall(nn)
              ratall(nn) = kariv
            end if
          end do

*         write(iot,'(2x,''ne = '',i5)') nk3
*       do n = 1, nk3
*         write(iot,'(i4,''. '',1pe15.7,1pe15.7)')
*    &                 n, engall(n), ratall(n)
*       end do

*         SUM UP the yield of the same energy
            nk4 = 0
          do n = 1, nk3-1
            if( ratall(n) .gt. 0.0d+00 ) then
              if( engall(n) .eq. engall(n+1) ) then
                ratall(n+1) = ratall(n) + ratall(n+1)
                ratall(n) = 0.0d+00
              else
                nk4 = nk4 + 1
              end if
            end if
          end do
            nk4 = nk4 + 1

      end if

      write(iot,'(/,"## Energy spectrum (final)")')
      if( inkf0(j33) .eq.  11 .or. inkf0(j33) .eq. -11 ) then
        if( nk4 .ne. 0 ) then
          write(iot,'(2x,"ne = ",i5)') (nk4-1)+(nk6-rimax1)
        else
          ! iauger = 2
          write(iot,'(2x,"ne = ",i5)')         (nk6-rimax1)
        end if

      else
        write(iot,'(2x,"ne = ",i5)')  nk4
      end if

C****
        nn = 0

      if( inkf0(j33) .eq.  11 .or. inkf0(j33) .eq. -11 ) then
        if( nk4 .ne. 0 ) then
          ngrp(j33) = (nk4-1)+(nk6-rimax1)
        else
          ! iauger = 2
          ngrp(j33) =         (nk6-rimax1)
        end if
          ngll(j33) =  1  ! positive @ e-type = 22, 23
c                         ! if negative, emin(1) must change from Zero.
      else
          ngrp(j33) = nk4
          ngll(j33) =  1  ! positive @ e-type = 22, 23
      end if

        call moddas_reallocate_dbl(isrc, j33, ngrp(j33), ngei, egmin)
        call moddas_reallocate_dbl(isrc, j33, ngrp(j33), ngea, egmax)
        call moddas_reallocate_dbl(isrc, j33, ngrp(j33), ngfe, fegrp)
        call moddas_reallocate_dbl(isrc, j33, ngrp(j33), ngft, rfe)
        call moddas_reallocate_dbl(isrc, j33, ngrp(j33), ngpi, prw)
        call moddas_reallocate_dbl(isrc, j33, ngrp(j33), ngpw, pwt)


        write(iot,'(3x,"$    Energy -lower, -upper  (MeV),"
     &                ,"  Activity X yield")')
      if( inkf0(j33) .eq. 11 .or. inkf0(j33) .eq. -11 .or.
     &    inkf0(j33) .eq. 2112 ) then  ! for *.BET, frtati2023/04/27 and SF
        if( nk4 .eq. 0 ) go to 710
        do n = 1, nk4-1
            egmin(ngei(j33)+n) = engall(n)
            egmax(ngea(j33)+n) = engall(n+1)
            fegrp(ngfe(j33)+n) = ratall(n)
            write(iot,'(i4,". ",1pe15.7,1pe15.7,x,1pe15.7)')
     &                   n, egmin(ngei(j33)+n), egmax(ngea(j33)+n),
     &                      fegrp(ngfe(j33)+n)
        end do

  710 continue
        if( nk4 .ne. 0 ) then
          nn = nk4-1 ! nn = 0 -> nk4-1  for Auger-CK and IE
        else
          nn = 0
        end if

        if( inkf0(j33) .eq. 11 .and. nk6 .gt. rimax1 ) then
          do n = 1, nk5-rimax1
            if( ratall(rimax1+n) .eq. 0.0d+00 ) then
            else
              nn = nn + 1
              egmin(ngei(j33)+nn) = engall(rimax1+n)
              egmax(ngea(j33)+nn) = engall(rimax1+n)
              fegrp(ngfe(j33)+nn) = ratall(rimax1+n)
            write(iot,'(i4,". ",1pe15.7,1pe15.7,x,1pe15.7)')
     &                  nn, egmin(ngei(j33)+nn), egmax(ngea(j33)+nn),
     &                      fegrp(ngfe(j33)+nn)
            end if
          end do

        end if

      else

          do n = 1, nk3
            if( ratall(n) .eq. 0.0d+00 ) then
            else
              nn = nn + 1
              egmin(ngei(j33)+nn) = engall(n)
              egmax(ngea(j33)+nn) = engall(n)
              fegrp(ngfe(j33)+nn) = ratall(n)
            write(iot,'(i4,". ",1pe15.7,1pe15.7,x,1pe15.7)')
     &                  nn, egmin(ngei(j33)+nn), egmax(ngea(j33)+nn),
     &                      fegrp(ngfe(j33)+nn)
            end if
          end do

      end if


C****

      if( norm(j33) .eq. 0 ) then  ! per Bq(0)
        write(iot,'("* Total:",1pe15.7," (Bq)"/"*"7x,1pe15.7
     &             ," (n/sec)")') abqsum, normfact(j33)
      else                         ! per particle
        write(iot,'("* Total:",1pe15.7," (/sec)")') normfact(j33)
          normfact(j33) = 1.0d+00
      end if

      write(iot,'(/,"# Output of sub. ddc2echo")')
*     CHECK
      do n = 1, nk1
        call ddc2echo(iot,n,j33,ierr)
      end do

*-----------------------------------------------------------------------

       goto 1000
*-----------------------------------------------------------------------
  920    m_err = 'RIsource.dat does not exist. Please check file(24).'
         ErrCha = ''
         ErrID = 'L:1375/R:setrisrc/F:risource.f'
           goto 998
  921    m_err = 'RIsource.dat does not exist. Please check file(24).'
         ErrCha = ''
         ErrID = 'L:1379/R:setrisrc/F:risource.f'
           goto 998
  922    m_err = 'Entry nuclide (stable?) can not be found in DECDC2.'
         ErrCha = ''
         ErrID = 'L:1383/R:setrisrc/F:risource.f'
           goto 998
  923    m_err = 'Read and search in RIsource.dat was failed.'
         ErrCha = ''
         ErrID = 'L:1387/R:setrisrc/F:risource.f'
           goto 998
  924    m_err = 'No decay database is available for nuclide: '//chaa ! T.Sato 2016/08/18
         ErrCha = ''
         ErrID = 'L:1391/R:setrisrc/F:risource.f'
           goto 998
  926    m_err = 'Atomic number of entry nuclide was exceeded 70.'
         ErrCha = ''
         ErrID = 'L:1395/R:setrisrc/F:risource.f'
           goto 998
  928    m_err = 'Data error in RIsource.dat (negative number)'
         ErrCha = ''
         ErrID = 'L:1399/R:setrisrc/F:risource.f'
           goto 998
  931    m_err = 'Invalid value (Bq) for Stable nuclide: '//chaa
         ErrCha = ''
         ErrID = 'L:1403/R:setrisrc/F:risource.f'
           goto 998

  942    m_err = 'RI source for alpha need RIsource.rad file.'
         ErrCha = ''
         ErrID = 'L:1408/R:setrisrc/F:risource.f'
           goto 998
  943    m_err = 'RI source for electron need RIsource.bet file.'
         ErrCha = ''
         ErrID = 'L:1412/R:setrisrc/F:risource.f'
           goto 998
  944    m_err = 'RI source for electron need RIsource.ack file.'
         ErrCha = ''
         ErrID = 'L:1416/R:setrisrc/F:risource.f'
           goto 998
  945    m_err = 'RI source for positron need RIsource.bet file.'
         ErrCha = ''
         ErrID = 'L:1420/R:setrisrc/F:risource.f'
           goto 998
  946    m_err = 'RI source for photon need RIsource.rad file.'
         ErrCha = ''
         ErrID = 'L:1424/R:setrisrc/F:risource.f'
           goto 998
  947    m_err = 'RI source for electron need RIsource.rad file.'
         ErrCha = ''
         ErrID = 'L:1428/R:setrisrc/F:risource.f'
           goto 998
  948    m_err = 'RI source for spontaneous fission neutron
     &            need RIsource.nsf file.'
         ErrCha = ''
         ErrID = 'L:1433/R:setrisrc/F:risource.f'
           goto 998
  952    m_err = 'RI source for alpha need RIsource.RAD file.'
         ErrCha = ''
         ErrID = 'L:1437/R:setrisrc/F:risource.f'
           goto 998
  953    m_err = 'RI source for electron need RIsource.BET file.'
         ErrCha = ''
         ErrID = 'L:1441/R:setrisrc/F:risource.f'
           goto 998
  954    m_err = 'RI source for electron need RIsource.ACK file.'
         ErrCha = ''
         ErrID = 'L:1445/R:setrisrc/F:risource.f'
           goto 998
  955    m_err = 'RI source for positron need RIsource.BET file.'
         ErrCha = ''
         ErrID = 'L:1449/R:setrisrc/F:risource.f'
           goto 998
  956    m_err = 'RI source for photon need RIsource.RAD file.'
         ErrCha = ''
         ErrID = 'L:1453/R:setrisrc/F:risource.f'
           goto 998
  957    m_err = 'RI source for electron need RIsource.RAD file.'
         ErrCha = ''
         ErrID = 'L:1457/R:setrisrc/F:risource.f'
           goto 998
  958    m_err = 'RI source for spontaneous fission neutron need
     &            RIsource.RAD file.'
         ErrCha = ''
         ErrID = 'L:1462/R:setrisrc/F:risource.f'
           goto 998

  961    m_err = 'RI source is limited for gamma-rays (photon): proj = '
     &           //pname(istyp(j33))(1:8)
         ErrCha = ''
         ErrID = 'L:1468/R:setrisrc/F:risource.f'
           goto 998
  969    m_err = 'proj is wrong or specified RI has too long half-life'
         ErrCha = ''
         ErrID = 'L:1472/R:setrisrc/F:risource.f'
           goto 998

  965    m_err = 'Number of radiation (1 nuclide) was exceeded rimax1.'
         ErrCha = ''
         ErrID = 'L:1477/R:setrisrc/F:risource.f'
           goto 998
  966    m_err = 'Total number of radiation was exceeded rimax2.'
         ErrCha = ''
         ErrID = 'L:1481/R:setrisrc/F:risource.f'
           goto 998
  967    m_err = 'Activities of all nuclides are 0.0 due to too long'
     &           //' decay time or too low activity.'
         ErrCha = ''
         ErrID = 'L:1486/R:setrisrc/F:risource.f'
           goto 998

  998 continue

        write(iot,'(/"***** Error Message from Input File *****"/)')

            icf = 200
         do i = 200, 1, -1
            if( m_err(i:i) .ne. ' ' ) goto 999
         end do
  999       icf = i

         call ErrWrite(ErrID, ErrCha)
         write(iot,'("  error = ",200a1/)') ( m_err(i:i), i=1,icf )

            ierr = 1

*-----------------------------------------------------------------------
 1000 continue
*       write(*,'(/''SUBROUTINE END!!''/)')
      close(iot)

      deallocate(chalctd)
      deallocate(br,r,nuclo)
      deallocate(eng,rat,engall,ratall)

      return
      end

************************************************************************
*                                                                      *
      subroutine rddecdc2(ntyp,nm,nucl,thf,nucld,bratio,ierr)
*                                                                      *
*       read DECDC2 index file and create nuclide list                 *
*         made by N.Matsuda on 2016/07/31                              *
*                                                                      *
************************************************************************

      implicit double precision (a-h,o-z)
      include 'risrcparam.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      common /error/  m_err, l_err, k_err
      character       m_err*200

      character chin*200, chlw*200, chcm*200

      integer  iot, iodecdc2, ioddc2r, ioddc2b, ioddc2a, ioddc2s
      integer  ierr

      character chaa*7
      real nucl, nucld, chza
      double precision thf, bratio
      dimension nucld(4), bratio(4)

      integer ntyp, klmax, kariint
      data klmax /100000/

*-----------------------------------------------------------------------
*     default values
*-----------------------------------------------------------------------
           m_err = ' Error in rddecdc2 sub.(read RIsource.dat)'
           ErrCha = ''
           ErrID = 'L:1552/R:rddecdc2/F:risource.f'
           l_err = 1
           k_err = 0
           ierr  = 0
*-----------------------------------------------------------------------
*     file number of scratch file for CHECKING decay calculation
           iot      = 820
*     file number of RIsource.dat or ICRP-07.NDX based on DECDC2
           iodecdc2 = 811
*     file number of RIsource.rad (radiation data) based on DECDC2
*          ioddc2r  = 812
*     file number of RIsource.bet (beta spectra) based on DECDC2
*          ioddc2b  = 813
*     file number of RIsource.ack (Auger-CK spectra) based on DECDC2
*          ioddc2a  = 814
*     file number of RIsource.dat (decay data) formatted with DCHAIN
           ioddc2s  = 815

c     err flag for this subroutine
           ierr     =   0
*-----------------------------------------------------------------------
*     START
*-----------------------------------------------------------------------
          rewind iodecdc2
*-----------------------------------------------------------------------
*     read Decay data file
*-----------------------------------------------------------------------
      if( ntyp .eq. 999 ) then  ! Search mode
        write(iot,'("## Read Decay data file (first 10 lines)",
     &              " by sub.rddecdc2")')

        do n = 1, klmax

          if( n .le. 10 ) then
            read( iodecdc2,'(a200)', iostat = ios ) chin
            if( ios .eq. -1 .and. n .eq. 1 ) go to 902
            if( ios .eq. -1 ) goto 903

            write(iot,'(i4,". ",a)') n, chin  ! lines: 1 to 10

          else
            read( iodecdc2,'()', iostat = ios )

            if( ios .eq. -1 ) then  ! end of data
              write(iot,'(8x,"..... (to be continued)")')
*             write(iot,'(''# '',78a)') ( '=', i = 1, 78 )

                nm = n
              write(iot ,'("/* Total number of the data: ",i6)') nm

                nucl = 0.0  ! dummy data
                thf  = 0.0d+00
              do k = 1, 4  ! initialization
                nucld(k)  = 0.0
                bratio(k) = 0.0d+00
              end do

                return
*               go to 300

            end if
          end if

        end do

          go to 904

      else  ! ntyp = 0 (ORIGINAL), or 1 (extra)
  200 continue
        chza = nucl

        if( ntyp .eq. 0 ) then  ! ORIGINAL

          call getelm(14,chaa,chza,ierr)

C       Sub. information: getelm(int.,char.,real,int.)
*       you can choose INPUT or OUTPUT style
*         INPUT:
*           2: chaa(ANY)   -> chza,8digit(1001.0(Zaid+))
*         OUTPUT:
*          11: chza(Zaid+) -> chaa,7digit(H-1    )
*          12: chza(Zaid+) -> chaa,7digit(H1     )
*          13: chza(Zaid+) -> chaa,7digit(1H     )
*          14: chza(Zaid+) -> chaa,7digit(%H   1 )

        else if( ntyp .eq. 1 ) then  ! extra

          call getelm(11,chaa,chza,ierr)

        else
          go to 999
        end if

          if( ierr .eq. 1 ) then
*           l_err = ill(jsn)
*           k_err = jsn
            go to 998
          end if

            thf = 0.0d+00  ! initialization

          do k = 1, 4  ! initialization
            nucld(k)  = 0.0
            bratio(k) = 0.0d+00
          end do

*-----------------------------------------------------------------------
  210 continue

        if( ntyp .eq. 0 ) then  ! Decay data file (default)

          do n = 1, nm
            read( iodecdc2,'(a200)', iostat = ios ) chin

            if( ios .eq. -1 ) exit  ! unregistered data

            if( chin(1:1) .eq. '%' ) then
*     ===== Nuclide: chin(3:8)
*             write(iot,'(''Nuclide: _|'',a6,''|_'')') chin(3:3+6-1)

C       Sub. information: getelm(int.,char.,real,int.)
*       you can choose INPUT or OUTPUT style
*         OUTPUT:
*          14: chza(Zaid+) -> chaa,7digit(%H   1 )

*         COMPARISON
              if( chin(3:8) .eq. chaa(2:7) ) then

*     =====  ZAID: chin(9:20)
*             write(iot,'(''ZAID: _|'',a12,''|_'')') chin(9:9+12-1)
*             read(chin(9:9+12-1),*) chza
*                                                             ==========

                write(iot,'(i6,". ",a)') n, chin

                call pickupd2(chin,ntyp,thf,nucld,bratio)
                  go to 300  ! registered data
*                 exit
              end if

            end if
          end do

        else if( ntyp .eq. 1 ) then  ! Decay data file (extra: ICRP-07)

          do n = 1, nm
            read( iodecdc2,'(a200)', iostat = ios ) chin

            if( ios .eq. -1 ) exit  ! unregistered data

*         COMPARISON
            if( chin(1:1+7-1) .eq. chaa ) then

              write(iot,'(i6,". ",a)') n, chin

              call pickupd2(chin,ntyp,thf,nucld,bratio)
                go to 300  ! registered data
*               exit
            end if

          end do

        else
          go to 905
        end if

      end if

*     Unregistered nuclide
        thf = -1.0d+00

  300 continue

*-----------------------------------------------------------------------

        goto 999

*-----------------------------------------------------------------------
  902    m_err = 'Data in DECDC2 index file is not found.'
         ErrCha = ''
         ErrID = 'L:1732/R:rddecdc2/F:risource.f'
           goto 998
  903    m_err = 'Decay data file (file(24)) was not read (, or < 10).'
         ErrCha = ''
         ErrID = 'L:1736/R:rddecdc2/F:risource.f'
           goto 998
  904    m_err = 'Decay data file exceeded the upper limit (klmax).'
         ErrCha = ''
         ErrID = 'L:1740/R:rddecdc2/F:risource.f'
           goto 998
  905    m_err = 'Decay data file is unknown.'
         ErrCha = ''
         ErrID = 'L:1744/R:rddecdc2/F:risource.f'
           goto 998
*-----------------------------------------------------------------------
  998 continue
         ierr = 1

  999 continue
*-----------------------------------------------------------------------
      return
       end

************************************************************************
*                                                                      *
      subroutine pickupd2(chin,ntyp,timehf,nucld,bratio)
*                                                                      *
*       data pickup from DECDC2 index file                             *
*         made by N.Matsuda on 2016/07/31                              *
*                                                                      *
************************************************************************

      implicit double precision (a-h,o-z)
      include 'risrcparam.inc'

*-----------------------------------------------------------------------

*     common /error/  m_err, l_err, k_err
*     character       m_err*200

      character chin*200
      integer ntyp, ndk

      integer  iot, iodecdc2, ioddc2r, ioddc2b, ioddc2a, ioddc2s
*     integer ierr

      character chaa*7, tnucl*7, rcha*2, brnt*10
      real nucl, nucld, chza
      double precision rnum, timehf, bratio
      dimension nucld(4), bratio(4)

*-----------------------------------------------------------------------
*     file number of scratch file for CHECKING decay calculation
           iot      = 820
*     file number of RIsource.dat or ICRP-07.NDX based on DECDC2
           iodecdc2 = 811
*     file number of RIsource.rad (radiation data) based on DECDC2
*          ioddc2r  = 812
*     file number of RIsource.bet (beta spectra) based on DECDC2
*          ioddc2b  = 813
*     file number of RIsource.ack (Auger-CK spectra) based on DECDC2
*          ioddc2a  = 814
*     file number of RIsource.dat (decay data) formatted with DCHAIN
*          ioddc2s  = 815

c     err flag for this subroutine
*          ierr     =   0
*-----------------------------------------------------------------------
*     data pickup
*-----------------------------------------------------------------------
  100 continue
      if( ntyp .eq. 0 ) then  ! Decay data file (default)
C     Nuclide (mother): nucla
          chaa = chin(3:3+6-1)//' '
        call getelm(2,chaa,chza,ierr)

C       Sub. information: getelm(int.,char.,real,int.)
*         INPUT:
*           2: chaa(ANY)   -> chza,8digit(1001.0(Zaid+))

*         if( ierr .eq. 1 ) then
*           l_err = ill(jsn)
*           k_err = jsn
*           go to 998
*         end if

*         read(chin(9:9+12-1),*) chza

        call getelm(11,chaa,chza,ierr)

C       Sub. information: getelm(int.,char.,real,int.)
*         OUTPUT:
*          11: chza(Zaid+) -> chaa,7digit(H-1    )

*         if( ierr .eq. 1 ) then
*           l_err = ill(jsn)
*           k_err = jsn
*           go to 998
*         end if

        tnucl = chaa
        nucl  = chza

C     Half life: timehf (sec)
*       write(iot,'(''Half life: _|'',a14,''|_'')') chin(33:33+14-1)
        read(chin(33:33+14-1),*) timehf

C     Nuclide (daughters) and the branching ratio
        read( iodecdc2,'(a200)', iostat = ios ) chin  ! COMMENT
        read( iodecdc2,'(a200)', iostat = ios ) chin  ! 3rd record
        read(chin(1:1+6-1),*) ndk
        if( ndk .eq. 0 ) then
C         Half life: timehf (sec)
            timehf = -1.0d+00
          do k = 1, 4
c         Nuclide (daughters): nucld
            nucld(k)  = 0.0
c         Branting ratio: bratio
            bratio(k) = 0.0d+00
          end do

        else
          if( ndk .gt. 4 ) ndk = 4  ! go to XXX <- error!!
          do k = 1, ndk
            read( iodecdc2,'(a200)', iostat = ios ) chin  ! 4th record
*           write(iot,'(''Daughter: _|'',a6,''|_'')') chin(43:43+6-1)
            chaa = chin(43:43+6-1)//' '

            if( chaa .eq. '       ' ) then
*           if( chin(1:1+12-1) .eq '        6.00' ) then
c         Nuclide (daughters): nucld
              nucld(k)  = 990000.0
            else

              call getelm(2,chaa,chza,ierr)

*               if( ierr .eq. 1 ) then
*                 l_err = ill(jsn)
*                 k_err = jsn
*                 go to 998
*               end if

c         Nuclide (daughters): nucld
              nucld(k) = chza
            end if
c         Branting ratio: bratio
*           write(iot,'(''Fraction: _|'',a6,''|_'')') chin(25:25+12-1)
            read(chin(25:25+12-1),'(d12.5)') bratio(k)

          end do

          if( ndk .eq. 4 ) then
          else
            do k = ndk+1, 4
c         Nuclide (daughters): nucld
              nucld(k)  = 0.0
c         Branting ratio: bratio
              bratio(k) = 0.0d+00
            end do
          end if

        end if

      else                    ! Decay data file (extra: ICRP-07.NDX)
C     Nuclide (mother): nucla
          chaa = chin(1:1+7-1)
        call getelm(2,chaa,chza,ierr)

C       Sub. information: getelm(int.,char.,real,int.)
*         INPUT:
*           2: chaa(ANY)   -> chza,8digit(1001.0(Zaid+))

*         if( ierr .eq. 1 ) then
*           l_err = ill(jsn)
*           k_err = jsn
*           go to 998
*         end if

        tnucl = chaa
        nucl  = chza

C     Half life: timehf (sec)
        read(chin(8:8+8-1),*) rnum
        rcha = chin(8+9-1:8+9-1+2-1)

        if( rcha .eq. 'y ' ) then
          timehf = dble( rnum*365.25*24.*60.*60. )
        else if( rcha .eq. 'd ' ) then
          timehf = dble( rnum*24.*60.*60. )
        else if( rcha .eq. 'h ' ) then
          timehf = dble( rnum*60.*60. )
        else if( rcha .eq. 'm ' ) then
          timehf = dble( rnum*60. )
        else if( rcha .eq. 's ' ) then
          timehf = dble( rnum )
        else if( rcha .eq. 'ms' ) then
          timehf = dble( rnum*1.0e-03 )
        else if( rcha .eq. 'um' ) then
          timehf = dble( rnum*1.0e-06 )
        else
        end if

C     Nuclide (daughters) and the branching ratio
        do k = 1, 4
          chaa = chin(54+25*(k-1):54+25*(k-1)+7-1)
          if( chaa .eq. '       ' ) then
c         Nuclide (daughters): nucld
            nucld(k)  = 0.0
c         Branting ratio: bratio
            bratio(k) = 0.0d+00
          else if( chaa .eq. 'SF     ' ) then
c         Nuclide (daughters): nucld
            nucld(k)  = 990000.0
c         Branting ratio: bratio
            brnt = chin(68+25*(k-1):68+25*(k-1)+10-1)
            brnt(7:7) = 'D'
            read(brnt,'(d12.5)') bratio(k)
          else
            call getelm(2,chaa,chza,ierr)

*             if( ierr .eq. 1 ) then
*               l_err = ill(jsn)
*               k_err = jsn
*               go to 998
*             end if

c         Nuclide (daughters): nucld
            nucld(k) = chza
c         Branting ratio: bratio
            brnt = chin(68+25*(k-1):68+25*(k-1)+10-1)
            brnt(7:7) = 'D'
            read(brnt,'(d12.5)') bratio(k)

          end if
        end do

C     Pointers of the other files
      end if

*-----------------------------------------------------------------------
*     CHECK and WRITE
*-----------------------------------------------------------------------
  200 continue
*     write(iot,'(4x,a,'' (Half life): '',1pe13.6,'' (sec)'')')
*    &            tnucl, timehf

      if( nucld(1) .ne. 0.0 ) then
        write(iot,fmt='(4x,a,f8.1)',advance='no') 'Daughter: ',nucld(1)
        do n = 2, 4
          if( nucld(n) .ne. 0.0 ) then
            write(iot,fmt='(a,f8.1)',advance='no') ', ',nucld(n)
          else
            write(iot,fmt='(a)') ''
*           go to 800
            exit
          end if
        end do
      else
      end if

 800  continue
      return
      end

************************************************************************
*                                                                      *
      subroutine lchain(j33,nk0,mi,nucld,br,ierr)
*                                                                      *
*       construct decay chain (linear) for all nuclides in list        *
*         made by N.Matsuda on 2016/07/31                              *
*                                                                      *
************************************************************************

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'risrcparam.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      common /error/  m_err, l_err, k_err
      character       m_err*200

*     NEW parameters ===
      common /risource/ chalct(nuclinmax,isrc), act000(nuclinmax,isrc),
     &                  activy(nuclinmax,isrc), thf(nuclinmax,isrc),
     &                  decayt(isrc)
      real  chalct
      double precision  act000, activy, thf, decayt

      integer  iot, iodecdc2, ioddc2r, ioddc2b, ioddc2a, ioddc2s

      real nucld(4,mi), nucl(mi)
      character rcha*1
      double precision br(4,mi), r(mi), rnum

      integer nk0, mi, mk, mcheck, ierr  ! nk1 = mi <- initial value
! T.Sato 2016/08/20, change to allocatable for OpenMP
      integer, allocatable:: nfch(:,:),long(:),ll(:)
      dimension a(mi), a0(mi)
      double precision deltt, actsum, asum00

      real chza, chza1, chza2
      character chaa*7, chaa1*7, chaa2*7
! T.Sato 2016/08/20, chage to allocatable for OpenMP
      real*8, allocatable :: b(:,:)
*     =>
      integer nfch1(maxchain), l1, l2
      dimension b1(maxchain), a1(mi)

*-----------------------------------------------------------------------
*     file number of scratch file for CHECKING decay calculation
           iot      = 820
*     file number of RIsource.dat or ICRP-07.NDX based on DECDC2
*          iodecdc2 = 811
*     file number of RIsource.rad (radiation data) based on DECDC2
*          ioddc2r  = 812
*     file number of RIsource.bet (beta spectra) based on DECDC2
*          ioddc2b  = 813
*     file number of RIsource.ack (Auger-CK spectra) based on DECDC2
*          ioddc2a  = 814
*     file number of RIsource.dat (decay data) formatted with DCHAIN
*          ioddc2s  = 815

      allocate(b(chainmax,maxchain))
      allocate(nfch(chainmax,maxchain), long(chainmax), ll(chainmax))

c     err flag for this subroutine
           ierr     =   0
*-----------------------------------------------------------------------
        i = 1
        k = 1       ! length
       nk1 = mi     ! nk1: registered nuclides
        mi = nk0
*-----------------------------------------------------------------------
*     Initialization
*-----------------------------------------------------------------------
      do n = 1, mi
        long(n) = 0  ! length of linear chain
        ll(n)   = 0  ! branching point of the linear chain (=length)
        nfch(n,k) = n  ! list number,
*                        n: number of linear chain  k: length
      end do

        asum00 = 0.0d+00  ! initial activity
        actsum = 0.0d+00  ! activity after cooling time

      do n = 1, nk1
c     Decay constant
        if( thf(n,j33) .eq. -1.0d+00 ) then  ! stable of SF nuclide
            r(n)  = 0.0d+00  ! log(2.)/infini.
          if( act000(n,j33) .ne. 0.0d+00 ) then
            chza = chalct(n,j33)
            call getelm(11,chaa,chza,ierr)

              if( ierr .eq. 1 ) then
*               l_err = ill(jsn)
*               k_err = jsn
*               close(iot)
*                 return
                go to 998
              end if

              go to 901  ! error
          else
            a0(n) = 0.0d+00  ! initial value (atoms)
          end if
        else  ! radioactive nuclide

            r(n)  = dble( log(2.) / thf(n,j33) )
            a0(n) = act000(n,j33) / r(n)  ! initial value (atoms)
        end if

c     Activity
          asum00 = asum00 + a0(n)

      end do

*-----------------------------------------------------------------------
*     Construct linear chain
*-----------------------------------------------------------------------
  100 continue
        mk = mi  ! mk(fix) and mi(var.): Total number of linear chain
        mcheck = 0    ! check - Continue or End  -
        do n = 1, mk
*         write (iot,'(''  n= '',i2,'', mk= '',i2)') n, mk
          if( long(n) .eq. 0 ) then
            do nn = 1, 4
              if( nucld(nn,nfch(n,k)) .eq. 0.0 ) then
                if( nn .eq. 1 ) then
*               Stable of SF nuclide
                  b(n,k) = br(nn,nfch(n,k))
                  long(n) = k
                  mcheck = mcheck + 1
                    go to 300
                else
*               Data not avaiable
                    go to 300
                end if
              else
*               Data check
                do i = 1, nk1
                  if( chalct(i,j33) .eq. nucld(nn,nfch(n,k)) ) then
*                 Registered nuclide
                    if( nn .eq. 1 ) then
                      nfch(n,k+1) = i
                      if( nucld(nn,nfch(n,k)) .ge. 990000.0 ) then
*                       Spontaneous fission nuclide
                        b(n,k) = br(nn,nfch(n,k))
                        long(n) = k + 1
                        mcheck = mcheck + 1
                          go to 300
                      else
                        b(n,k) = br(nn,nfch(n,k))
*                       ll(n) = ...
*                       long(n) = 0
                          go to 200
                      end if
                    else
*                   construct NEW linear chain; mi = mi + 1
                      nfch(mi+nn-1,k+1) = i
                      if( nucld(nn,nfch(n,k)) .ge. 990000.0 ) then
*                       Spontaneous fission nuclide
                        b(mi+nn-1,k) = br(nn,nfch(n,k))
                        ll(mi+nn-1) = k
                        long(mi+nn-1) = k + 1
                        mcheck = mcheck + 1
                          go to 300
                      else
                        b(mi+nn-1,k) = br(nn,nfch(n,k))
                        ll(mi+nn-1) = k
                        long(mi+nn-1) = 0
                          go to 200
                      end if
                    end if
                  end if
                end do
*                 Unregistered nuclide
                  go to 920 ! error
  200           continue
              end if
            end do
  300     continue

            if( nn .gt. 1 ) then
              if( nucld(nn,nfch(n,k)) .eq. 0.0 ) nn = nn - 1  ! last
*             complement NEW linear chain
              do i = mi + 1, mi + nn - 1
                do kk = 1, k - 1
                  nfch(i,kk) = nfch(n,kk)
                     b(i,kk) =    b(n,kk)
                end do
                  nfch(i,k)  = nfch(n,k)
              end do
            end if
              mi = mi + nn - 1
          else
*         completed linear chain
              mcheck = mcheck + 1
          end if

        end do

          k = k + 1

        if( mi .gt. chainmax ) goto 921 ! error
        if( mi .gt. mcheck ) goto 100

        write(iot,'("* Total number of linear chain:",i4)') mi
        write(iot,'("* Maximum length of linear chain:",i3)') (k-1)

*-----------------------------------------------------------------------
  500 continue
      write(iot,'("## RESULT - linear chain -")')
      do i = 1, mi  ! Num. of linear chain

        write (iot,fmt='(i4,a)',advance='no') i, '. '
        do k = 1, long(i) - 1
*       Chains ( from 1 to long-1 )
          chza = chalct(nfch(i,k),j33)
          call getelm(11,chaa,chza,ierr)

            if( ierr .eq. 1 ) then
*             l_err = ill(jsn)
*             k_err = jsn
*             close(iot)
*               return
              go to 998
            end if

          write (iot,fmt='(a)',advance='no') chaa//'->'
        end do

*       Chain end
        if( chalct(nfch(i,long(i)),j33) .ge. 990000.0 ) then
          write (iot,fmt='(a)',advance='no') 'SF(end)'
        else
          chza = chalct(nfch(i,long(i)),j33)
          call getelm(11,chaa,chza,ierr)

            if( ierr .eq. 1 ) then
*             l_err = ill(jsn)
*             k_err = jsn
*             close(iot)
*               return
              go to 998
            end if

          write (iot,fmt='(a)',advance='no') chaa//'(end)'
        end if
*       Chain length
          write (iot,fmt='(a9,i2)',advance='no') ', length:', long(i)
          write (iot,fmt='(a9,i2)') ', branch:', ll(i)

*       half lives (arbitrary) and
        if( long(i) .gt. 1 ) then
          write (iot,fmt='(4x,a4)',advance='no') 'HL  '
          do k = 1, long(i) - 1
            write (iot,fmt='(1pe9.2)',advance='no') r(nfch(i,k))
          end do
            write (iot,fmt='(a)') ''

*       branching ratio
          write (iot,fmt='(4x,a5)',advance='no') 'b:   '
          do k = 1, long(i) - 1
            write (iot,fmt='(f7.4,a2)',advance='no') b(i,k), '  '
          end do
            write (iot,fmt='(a)') ''
        else
        end if

      end do

*-----------------------------------------------------------------------
*     avoidance of double count
*-----------------------------------------------------------------------
  800 continue
      if( decayt(j33) .ge. 0.0d+00 ) go to 810
*     write(iot,'(''## CHECK and DELETE (same linear chain)'')')
      do i = 1, mi  ! Num. of linear chain
        if( long(i) .gt. 1 ) then
          do k = 2, long(i)

            if( nfch(i,k) .gt. nk0 ) then
            else

              do nn = 1, mi
                if( nfch(nn,1) .eq. nfch(i,k) ) then
*                 write (iot,'(2x,''nn='',i2)') nn
*                 long(nn) = 0  ! skip this linear chain

                  chza1 = chalct(nfch(i,1),j33)
                  call getelm(11,chaa1,chza1,ierr)

                    if( ierr .eq. 1 ) then
*                     l_err = ill(jsn)
*                     k_err = jsn
*                     close(iot)
*                       return
                      go to 998
                    end if

                  chza2 = chalct(nfch(i,k),j33)
                  call getelm(11,chaa2,chza2,ierr)

                    if( ierr .eq. 1 ) then
*                     l_err = ill(jsn)
*                     k_err = jsn
*                     close(iot)
*                       return
                      go to 998
                    end if

                  go to 980  ! decayt(j33) < 0.0

                end if
              end do

            end if

          end do
        end if
      end do

        go to 810

      do n = 1, mi
        if( long(n) .eq. 0 ) then
          write (iot,'(2x,i4,". DELETE, length: ",i2)') n, long(n)
        else
          write (iot,'(2x,i4,". OK.   , length: ",i2)') n, long(n)
        end if
      end do

*-----------------------------------------------------------------------
*     Bateman calculation for each linear chain
*-----------------------------------------------------------------------
  810 continue
***   INITIAL VALUE
        actsum = 0.0d+00  ! (atoms)

***   SET PARAMETERS
      do n = 1, mi
        if( long(n) .eq. 0 ) then  ! skip
        else
          if( n .eq. 1 ) then
            write(iot,'(  "# ",78a)') ( '=', i = 1, 78 )
            write(iot,'("## CALCULATION START - Bateman solution -")')
            if( decayt(j33) .ge. 0.0 ) then
              write(iot,'(2x,"Cooling time: ",1pe12.5," (sec)")')
     &                    decayt(j33)
            else
              write(iot,'(2x,"Cooling time: ",1pe12.5," (times)")')
     &                    abs( decayt(j33) )
            end if

          end if
*         do nn = 1, nk1
*           a1(nn) = a0(nn)
*         end do
            l1 = ll(n)
            l2 = long(n)
          do nn = 1, l2
            nfch1(nn) = nfch(n,nn)
            if( nn .ne. l2 ) b1(nn) = b(n,nn)
          end do

*     Bateman calculation
          write(iot,'(2x,"Linear chain No.",i3)') n
          call batesol(j33,nk1,nfch1,l1,l2,a0,r,b1,ierr)
            if(ierr .eq. 1 ) go to 999

*         do nn = 1, nk1
*           a(nn) = a(nn) + a1(nn)
*         end do
        end if
      end do

***   RESULT CHECK (atoms)
      write (iot,'("*- Results  ----------------------------")')
*     do n = 1, nk1
*       write(iot,'(2x,i2,''. '',f8.1,'', '',1pe15.7,'' ->'',1pe15.7
*    &        ,'' (atoms)'')') n, chalct(n,j33), a0(n), a(n)
*       asumaf = asumaf + a(n)
*     end do

*-----------------------------------------------------------------------

        goto 999

*-----------------------------------------------------------------------
  901    m_err = 'Activity entry was wrong for stable nuclide.'
         ErrCha = ''
         ErrID = 'L:2388/R:lchain/F:risource.f'
         goto 998
  920    m_err = 'Target nuclide is not found in the listed DATA.'
         ErrCha = ''
         ErrID = 'L:2392/R:lchain/F:risource.f'
         goto 998
  921    m_err = 'Number of linear chain was exceeded chainmax.' ! T.Sato 2016/08/18
         ErrCha = ''
         ErrID = 'L:2396/R:lchain/F:risource.f'
         goto 998
  980    m_err = 'Half-life decay mode (dtime < 0.0) cannot be used '
     &   //'when both parent and daughter nuclides ('//chaa1//' and '
     &   //chaa2//') are specified.'
         ErrCha = ''
         ErrID = 'L:2402/R:lchain/F:risource.f'
         goto 998

*-----------------------------------------------------------------------
  998 continue
         ierr = 1

  999 continue
      deallocate(b)
      deallocate(nfch,long,ll)

*-----------------------------------------------------------------------
      return
      end

************************************************************************
*                                                                      *
      subroutine batesol(j33,nk1,nfch,l1,l2,a0,r,b,ierr)
*                                                                      *
*       bateman solution (decay mode) for linear chain                 *
*         made by N.Matsuda on 2016/07/31                              *
*                                                                      *
************************************************************************

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'risrcparam.inc'

*-----------------------------------------------------------------------

      common /error/  m_err, l_err, k_err
      character       m_err*200

*     NEW parameters ===
      common /risource/ chalct(nuclinmax,isrc), act000(nuclinmax,isrc),
     &                  activy(nuclinmax,isrc), thf(nuclinmax,isrc),
     &                  decayt(isrc)
      real  chalct
      double precision  act000, activy, thf, decayt

      integer  iot, iodecdc2, ioddc2r, ioddc2b, ioddc2a, ioddc2s

      integer j33, nk1, nfch, l1, l2, ierr
      double precision a, a0, r, b, t
      dimension nfch(l2), a(nk1), a0(nk1), r(nk1), b(l2-1)

      integer mi, mk, j, mcheck
      double precision ramdat, dex, bunsi, bunbo, s, ssame
      double precision qg, qs

*     Parameters
*       nfch: linear chain
*       l1:   Branching point of the linear chain
*       l2:   Maximum length of linear chain
*       r:    Decay constant (/sec): from 1 to nk1
*       b:    Branching ratio
*       a0:   Number of atoms

      data smallv / 1.0d-36 /

*-----------------------------------------------------------------------
*     file number of scratch file for CHECKING decay calculation
           iot      = 820
*     file number of RIsource.dat or ICRP-07.NDX based on DECDC2
*          iodecdc2 = 811
*     file number of RIsource.rad (radiation data) based on DECDC2
*          ioddc2r  = 812
*     file number of RIsource.bet (beta spectra) based on DECDC2
*          ioddc2b  = 813
*     file number of RIsource.ack (Auger-CK spectra) based on DECDC2
*          ioddc2a  = 814
*     file number of RIsource.dat (decay data) formatted with DCHAIN
*          ioddc2s  = 815

c     err flag for this subroutine
           ierr     =   0
*-----------------------------------------------------------------------
***   INITIALIZE
      do n = 1, nk1
        a(n) = 0.0d+00  ! (atoms)
      end do

      if( l2 .eq. 1 ) then  ! stable or SF nuclide
        a(nfch(1)) = a0(nfch(1))  ! may be 0.0d+00
          go to 300
      end if

C     Cooling time
      if( decayt(j33) .eq. 0.0d+00 ) then
        go to 500
      else if( decayt(j33) .gt. 0.0d+00 ) then
        t = decayt(j33)
      else  ! Half life (1st entry)
        t = abs( decayt(j33) * dble( log(2.) / r(nfch(1)) ) )
      end if
        write(iot,'(4x,"Time =",1pe12.5," (sec)")') t

***   CALCULATION START!! (Main routine)
      do k = 1, l2  ! Target nuclide: nucl(nfch(k))
        if( k .le. l1 ) go to 260  ! skip <- double count

  100 continue
          j = 1  ! Fixing

  130   continue
          a00 = a0(nfch(j))
          if( a00 .eq. 0.0d+00 ) go to 260  ! error
          qg = 0.0d+00
          qs = 0.0d+00
          ksi = 0
          ksm = 0
        do ki = j, k  ! bateman factor: loop from j to k
          if( j .eq. k ) then  ! j = k = 1
            if( r(nfch(k)) .eq. 0.0d+00 ) then  ! stable
*             write (iot,'(2x,''stable'')')
              a(nfch(k)) = a0(nfch(k))
            else
*             write (iot,'(2x,''decay - itself -'')')
              ramdat = r(nfch(ki)) * t
              if( ramdat .le. 1.0d-10 ) then
                dex = 1.0d+00 - ramdat
              else
                dex = exp(-ramdat)
              end if
                a(nfch(k)) = a(nfch(k)) + a00 * dex
            end if
          else
            bunsi = 1.0d+00
            bunbo = 1.0d+00
*             s   = 1.0d+00  ( = bunsi / bunbo )
            do km = j, k  ! bateman factor: loop from j to k
              if( km .ne. k ) then  ! km = j, k-1, 1
                bunsi = bunsi * r(nfch(km)) * b(km)
*               write(iot,'(2x,''bunsi='',1pe12.5)') bunsi
              end if
              if( km .eq. ki ) then  ! km .eq. ki -> skip
              else
                rmri = abs(r(nfch(km))-r(nfch(ki))) * t
                if( rmri .ge. 1.0d-20 ) then
                  bunbo = bunbo * (r(nfch(km))-r(nfch(ki)))
*                 write(iot,'(2x,''bunbo='',1pe12.5)') bunbo
                else
                  ksi = ki
                  ksm = km
                  go to 190  ! loop stop
                end if
              end if
            end do
*             s =  bunsi / bunbo
  162     continue
              ramdat = r(nfch(ki)) * t
            if( ramdat .le. 1.0d-10 ) then
              dex = 1.0d+00 - ramdat
            else
              dex = exp(-ramdat)
            end if

  180     continue

              qg = qg + ( bunsi / bunbo ) * dex
              qs = qs + ( bunsi / bunbo )

          end if
  190   continue
        end do

*     branch for same lambda
        if( ksm .ne. 0 ) then  ! different or same lambda
*         write (iot,'(2x,''same lambda'')')
          bunsi = 1.0d+00
          bunbo = 1.0d+00
          ssame = 1.0d+00  ! ( = bunsi / bunbo )
          do km = j, k  ! bateman factor
            if( km .ne. k ) then  ! km = j, k-1, 1
              bunsi = bunsi * r(nfch(km)) * b(km)
            end if
            if( km .eq. ksi .or. km .eq. ksm ) then
            else
              rmri = abs(r(nfch(km))-r(nfch(ksm))) * t
              if( rmri .ge. 1.0d-20 ) then
                bunbo = bunbo * (r(nfch(km))-r(nfch(ksm)))
              end if
            end if
          end do
            ssame = bunsi / bunbo
  192   continue
            ramdat = r(nfch(ksm)) * t
          if( ramdat .le. 1.0d-10 ) then
            dex = 1.0d+00 - ramdat
          else
            dex = exp(-ramdat)
          end if

  193   continue

            qg = qg - ( qs * dex ) + ( ssame * t * dex )

        end if

            a(nfch(k)) = a(nfch(k)) + a00 * qg

*         Normalization for auto-setting mode ( dtime < 0.0 )
          if( decayt(j33) .lt. 0.0d+00 ) then
            a(nfch(k)) = ( 2.0**abs(decayt(j33)) ) * a(nfch(k))
          end if

  210 continue
*       j = j + 1
*       if( j .gt. k ) go to 260
*         go to 130
  260 continue
      end do

  300 continue
*     if( a(nfch(k)) .le. smallv ) a(nfch(k)) = 0.0d+00

*-----------------------------------------------------------------------
*     Check the result
*-----------------------------------------------------------------------
  500 continue

      do n = 1, nk1

        activy(n,j33) = activy(n,j33) + a(n)

        write(iot,'(4x,i2,". ",1pe12.5," ->",1pe12.5,"  Total: ",
     &            1pe12.5," (atoms)")') n, a0(n), a(n), activy(n,j33)
      end do

*-----------------------------------------------------------------------

        goto 999

*-----------------------------------------------------------------------
  998 continue
         ierr = 1

  999 continue

*-----------------------------------------------------------------------
      return
      end

************************************************************************
*                                                                      *
      subroutine rdgamma(ntyp,chza,act,kkk,eng,rat,ierr)
*                                                                      *
*       read RIsource.dat file (energy and yield)                      *
*         made by N.Matsuda on 2016/07/31                              *
*         revised by N.Matsuda on 2017/05/29                           *
*                                                                      *
************************************************************************

      implicit double precision (a-h,o-z)
      include 'risrcparam.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      common /error/  m_err, l_err, k_err
      character       m_err*200

      character chin*200, chlw*200, chcm*200

      integer  iot, iodecdc2, ioddc2r, ioddc2b, ioddc2a, ioddc2s
      integer ierr

      common /paran/  icfn(100), ilfn(100), chfn(100)
      character chfn*200
      character fname*200 ! T.Sato 2022/11/30, increase from 100 to 200 because inumc requests 200 character

      logical lex

      character chaa*7
      real chza
*     parameter ( rimax1    =  4,000 )
*     parameter ( rimax2    = 20,000 )
      double precision act, eng(rimax1), rat(rimax1)
      integer ntyp, icl, kkk, ndk, nsp, nnn

      integer klmax
      data klmax /500000/  !  50,977 lines @ RIsource.dat
                           ! 455,625 lines @ ICRP-07.RAD
      double precision smallv
      data smallv / 1.0d-36 /

*-----------------------------------------------------------------------
*     default values
*-----------------------------------------------------------------------
           m_err = ' Error in rdgamma sub.(read RIsource.dat file)'
           ErrCha = ''
           ErrID = 'L:2696/R:rdgamma/F:risource.f'
           l_err = 1
           k_err = 0
           ierr  = 0
*-----------------------------------------------------------------------
*     file number of scratch file for CHECKING decay calculation
           iot      = 820
*     file number of RIsource.dat or ICRP-07.NDX based on DECDC2
*          iodecdc2 = 811
*     file number of RIsource.rad (radiation data) based on DECDC2
           ioddc2r  = 812
*     file number of RIsource.bet (beta spectra) based on DECDC2
*          ioddc2b  = 813
*     file number of RIsource.ack (Auger-CK spectra) based on DECDC2
*          ioddc2a  = 814
*     file number of RIsource.dat (decay data) formatted with DCHAIN
*          ioddc2s  = 815

c     err flag for this subroutine
           ierr     =   0
*-----------------------------------------------------------------------
*     address check of decay data file
*-----------------------------------------------------------------------

        call fchkd2(24,ntyp,fname,ierr)

          if( ierr .eq. 1 ) then
*           l_err = ill(jsn)
*           k_err = jsn
            go to 998
          end if

        icl = inumc(fname,1,100,' ') - 1
*       do n = 100, 1, -1
*         if( fname(n:n) .ne. ' ' ) go to 201
*       end do
*
* 201     icl = n
          lex = .false.

        if( ntyp .eq. 0 ) then

          fname = fname(1:icl-3)//'rad'
          inquire( file = fname, exist = lex )
          if( .not. lex ) then
            if( icl .gt. 12 ) then
              inquire( file = fname(1:icl-12)//'ICRP-07.RAD ',
     &                 exist = lex )
              if( lex ) then
                fname = fname(1:icl-12)//'ICRP-07.RAD '
              else
                fname = fname(1:icl-3)//'dat'
              end if
            else if( icl .eq. 12 ) then
              inquire( file = 'ICRP-07.RAD', exist = lex )
              if( lex ) then
                fname = 'ICRP-07.RAD '
              else
                fname = fname(1:icl-3)//'dat'
              end if
            end if
          end if

        else  ! ntyp = 1

          fname = fname(1:icl-3)//'RAD'

        end if

      open( ioddc2r, file = fname,
     &      form='formatted', status = 'unknown' )
        rewind ioddc2r

*-----------------------------------------------------------------------
*     Target nuclide
*-----------------------------------------------------------------------
*         chza = nucl
        call getelm(11,chaa,chza,ierr)

C       Sub. information: getelm(int.,char.,real,int.)
*         OUTPUT:
*          11: chza(Zaid+) -> chaa,7digit(H-1    )

          if( ierr .eq. 1 ) then
*           l_err = ill(jsn)
*           k_err = jsn
*             return
            go to 999
          end if

*-----------------------------------------------------------------------
*     Read the data on the energy and yield of gamma-ray
*-----------------------------------------------------------------------
  220 continue

        nnn = -1
        kkk =  0
      if( ntyp .eq. 0 .and. fname(icl-2:icl) .eq. 'dat' ) then  ! RIsource.dat

        call getelm(14,chaa,chza,ierr)

C       Sub. information: getelm(int.,char.,real,int.)
*         OUTPUT:
*          14: chza(Zaid+) -> chaa,7digit(%H   1 )

          if( ierr .eq. 1 ) then
*           l_err = ill(jsn)
*           k_err = jsn
*             return
            go to 999
          end if

        do n = 1, klmax
          read( ioddc2r,'(a200)', iostat = ios ) chin
*         write(iot,'(i6,'': '',a)') n, chin

          if( n .eq. 1 ) then
            if( ios .eq. -1 ) go to 903
          else
            if( ios .eq. -1 ) go to 280  ! read end
          end if

          if( nnn .eq. 0 ) then
            go to 280  ! successfully end
          else if( nnn .gt. 0 ) then
            do nn = 1, 3
              nnn = nnn - 1
              kkk = kkk + 1
              read(chin(1+12*(2*nn-2):12*(2*nn-1)),'(e15.7)') eng(kkk)
              read(chin(1+12*(2*nn-1):12*(2*nn)),'(e15.7)')   rat(kkk)
              if( kkk .eq. 1 ) then
                write(iot,'(2x,1pe15.7,1pe15.7,2x,
     &                    "$ Energy (MeV) and Absolute yield")')
     &                      eng(kkk), rat(kkk)

              else
                write(iot,'(2x,1pe15.7,1pe15.7)') eng(kkk), rat(kkk)

              end if

              if( nnn .eq. 0 ) go to 280

            end do

          end if

          if( chin(1:1) .eq. '%' ) then
            if( chin(3:3+6-1) .eq. chaa(2:7) ) then
              read( ioddc2r, '()' )
              read( ioddc2r,'(a200)', iostat = ios ) chin
              read(chin(7:7+6-1),'(i6)') nsp
              if( nsp .eq. 0 ) then
                  go to 280  ! no gamma-rays
              end if
              read(chin(1:1+6-1),*) ndk
              if( ndk .eq. 0 ) then
                  go to 280  ! no daughter
              end if
              do nn = 1, ndk
                read( ioddc2r, '()' )
              end do
              read( ioddc2r,'(a200)', iostat = ios ) chin
              read(chin(7:7+6-1),*) nnn
*               if( nnn .eq. 0 ) go to 280
                if( nnn .gt. rimax1 ) go to 965

            else
              read( ioddc2r, '()' )
              read( ioddc2r, '()' )

            end if
          end if
        end do

          go to 940

      else  ! ntyp = 1 (Extra file: ICRP-07.RAD)

        do n = 1, klmax
          read( ioddc2r,'(a200)', iostat = ios ) chin
*         write(iot,'(i6,'': '',a)') n, chin

          if( n .eq. 1 ) then
            if( ios .eq. -1 ) go to 903
          else
            if( ios .eq. -1 ) go to 280  ! read end
          end if

          if( nnn .eq. 0 ) then
            go to 280  ! successfully end
          else if( nnn .gt. 0 ) then
            nnn = nnn - 1
            if( chin(27:27+3-1) .eq. '  G' ) then
                kkk = kkk + 1
              if( kkk .gt. rimax1 ) go to 965
              read(chin(15:15+12-1),'(e15.7)') eng(kkk)
              read(chin(3:3+12-1),'(e15.7)')   rat(kkk)
              if( kkk .eq. 1 ) then
                write(iot,'(2x,1pe15.7,1pe15.7,2x,
     &                    "$ Energy (MeV) and Absolute yield")')
     &                      eng(kkk), rat(kkk)

              else
                write(iot,'(2x,1pe15.7,1pe15.7)') eng(kkk), rat(kkk)

              end if

            end if

          end if

          if( chin(1:1+7-1) .eq. chaa ) then
            read(chin(21:21+9-1),'(i9)') nnn

          end if
        end do

      end if

      if( act .lt. 0.0d+00 ) then
      else
        go to 940
      end if

  280 continue

*-----------------------------------------------------------------------

        goto 999

*-----------------------------------------------------------------------
  901    m_err = 'RIsource.rad does NOT exist.'
         ErrCha = ''
         ErrID = 'L:2929/R:rdgamma/F:risource.f'
           goto 998
  903    m_err = 'Data in RIsource.rad is not read.'
         ErrCha = ''
         ErrID = 'L:2933/R:rdgamma/F:risource.f'
           goto 998
  940    m_err = 'Klmax in sub. rdgamma is smaller than the data line.'
         ErrCha = ''
         ErrID = 'L:2937/R:rdgamma/F:risource.f'
           goto 998
  965    m_err = 'Number of gamma-rays was exceeded rimax1.'
         ErrCha = ''
         ErrID = 'L:2941/R:rdgamma/F:risource.f'
           goto 998
*-----------------------------------------------------------------------
  998 continue
         ierr = 1

  999 continue
*-----------------------------------------------------------------------

      close(ioddc2r)
      return
      end

************************************************************************
*                                                                      *
      subroutine rdgamma2(ntyp,chza,act,kkk,eng,rat,ierr)
*                                                                      *
*       read RIsource.rad file (energy and yield)                      *
*         made by N.Matsuda on 2018/08/15                              *
*                                                                      *
************************************************************************

      implicit double precision (a-h,o-z)
      include 'risrcparam.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      common /error/  m_err, l_err, k_err
      character       m_err*200

      character chin*200, chlw*200, chcm*200

      integer  iot, iodecdc2, ioddc2r, ioddc2b, ioddc2a, ioddc2s
      integer ierr

      common /paran/  icfn(100), ilfn(100), chfn(100)
      character chfn*200
      character fname*200 ! T.Sato 2022/11/30, increase from 100 to 200 because inumc requests 200 character

      logical lex

      character chaa*7
      real chza
*     parameter ( rimax1    =  4,000 )
*     parameter ( rimax2    = 20,000 )
      double precision act, eng(rimax1), rat(rimax1)
      integer ntyp, icl, kkk, ndk, nsp, nnn

      integer klmax
      data klmax /500000/  !  50,977 lines @ RIsource.dat
                           ! 455,625 lines @ ICRP-07.RAD
      double precision smallv
      data smallv / 1.0d-36 /

*-----------------------------------------------------------------------
*     default values
*-----------------------------------------------------------------------
           m_err = ' Error in rdgamma2 sub.(read RIsource.rad file)'
           ErrCha = ''
           ErrID = 'L:3001/R:rdgamma2/F:risource.f'
           l_err = 1
           k_err = 0
           ierr  = 0
*-----------------------------------------------------------------------
*     file number of scratch file for CHECKING decay calculation
           iot      = 820
*     file number of RIsource.dat or ICRP-07.NDX based on DECDC2
*          iodecdc2 = 811
*     file number of RIsource.rad (radiation data) based on DECDC2
           ioddc2r  = 812
*     file number of RIsource.bet (beta spectra) based on DECDC2
*          ioddc2b  = 813
*     file number of RIsource.ack (Auger-CK spectra) based on DECDC2
*          ioddc2a  = 814
*     file number of RIsource.dat (decay data) formatted with DCHAIN
*          ioddc2s  = 815

c     err flag for this subroutine
           ierr     =   0
*-----------------------------------------------------------------------
*     address check of decay data file
*-----------------------------------------------------------------------

        call fchkd2(24,ntyp,fname,ierr)

          if( ierr .eq. 1 ) then
*           l_err = ill(jsn)
*           k_err = jsn
            go to 998
          end if

        icl = inumc(fname,1,100,' ') - 1
*       do n = 100, 1, -1
*         if( fname(n:n) .ne. ' ' ) go to 201
*       end do
*
* 201     icl = n
          lex = .false.

        if( ntyp .eq. 0 ) then

          fname = fname(1:icl-3)//'rad'
          inquire( file = fname, exist = lex )
          if( .not. lex ) then
            if( icl .gt. 12 ) then
              inquire( file = fname(1:icl-12)//'ICRP-07.RAD ',
     &                 exist = lex )
              if( lex ) then
                fname = fname(1:icl-12)//'ICRP-07.RAD '
              else
                fname = fname(1:icl-3)//'dat'  ! Error
                go to 901
              end if
            else if( icl .eq. 12 ) then
              inquire( file = 'ICRP-07.RAD', exist = lex )
              if( lex ) then
                fname = 'ICRP-07.RAD '
              else
                fname = fname(1:icl-3)//'dat'  ! Error
                go to 901
              end if
            end if
          end if

        else  ! ntyp = 1

          fname = fname(1:icl-3)//'RAD'

        end if

      open( ioddc2r, file = fname,
     &      form='formatted', status = 'unknown' )
        rewind ioddc2r

*-----------------------------------------------------------------------
*     Target nuclide
*-----------------------------------------------------------------------
*         chza = nucl
        call getelm(11,chaa,chza,ierr)

C       Sub. information: getelm(int.,char.,real,int.)
*         OUTPUT:
*          11: chza(Zaid+) -> chaa,7digit(H-1    )

          if( ierr .eq. 1 ) then
*           l_err = ill(jsn)
*           k_err = jsn
*             return
            go to 999
          end if

*-----------------------------------------------------------------------
*     Read the data on the energy and yield of gamma-ray
*-----------------------------------------------------------------------
  220 continue

        nnn = -1

C     else  ! ntyp = 1 (Extra file: ICRP-07.RAD)
  250 continue

        do n = 1, klmax
          read( ioddc2r,'(a200)', iostat = ios ) chin
*         write(iot,'(i6,'': '',a)') n, chin

          if( n .eq. 1 ) then
            if( ios .eq. -1 ) go to 903
          else
            if( ios .eq. -1 ) go to 280  ! read end
          end if

          if( nnn .eq. 0 ) then
            go to 280  ! successfully end
          else if( nnn .gt. 0 ) then
            nnn = nnn - 1
            if( chin(27:27+3-1) .eq. '  X' ) then
                kkk = kkk + 1
              if( kkk .gt. rimax1 ) go to 965
              read(chin(15:15+12-1),'(e15.7)') eng(kkk)
              read(chin(3:3+12-1),'(e15.7)')   rat(kkk)
              if( kkk .eq. 1 ) then
                write(iot,'(2x,1pe15.7,1pe15.7,2x,
     &                    "$ Energy (MeV) and Absolute yield")')
     &                      eng(kkk), rat(kkk)

              else
                write(iot,'(2x,1pe15.7,1pe15.7)') eng(kkk), rat(kkk)

              end if

            end if

          end if

          if( chin(1:1+7-1) .eq. chaa ) then
            read(chin(21:21+9-1),'(i9)') nnn

          end if
        end do


      if( act .lt. 0.0d+00 ) then
      else
        go to 940
      end if

  280 continue

*-----------------------------------------------------------------------

        goto 999

*-----------------------------------------------------------------------
  901    m_err = 'RIsource.rad does NOT exist.'
         ErrCha = ''
         ErrID = 'L:3157/R:rdgamma2/F:risource.f'
           goto 998
  903    m_err = 'Data in RIsource.rad is not read.'
         ErrCha = ''
         ErrID = 'L:3161/R:rdgamma2/F:risource.f'
           goto 998
  940    m_err = 'Klmax in sub. rdgamma2 is smaller than the data line.'
         ErrCha = ''
         ErrID = 'L:3165/R:rdgamma2/F:risource.f'
           goto 998
  965    m_err = 'Number of gamma-rays was exceeded rimax1.'
         ErrCha = ''
         ErrID = 'L:3169/R:rdgamma2/F:risource.f'
           goto 998
*-----------------------------------------------------------------------
  998 continue
         ierr = 1

  999 continue
*-----------------------------------------------------------------------

      close(ioddc2r)
      return
      end

************************************************************************
*                                                                      *
      subroutine rdgamma3(ntyp,chza,act,kkk,eng,rat,ierr)
*                                                                      *
*       read RIsource.rad file (energy and yield)                      *
*         made by N.Matsuda on 2019/05/08  for annihilation photons    *
*                                                                      *
************************************************************************

      implicit double precision (a-h,o-z)
      include 'risrcparam.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      common /error/  m_err, l_err, k_err
      character       m_err*200

      character chin*200, chlw*200, chcm*200

      integer  iot, iodecdc2, ioddc2r, ioddc2b, ioddc2a, ioddc2s
      integer ierr

      common /paran/  icfn(100), ilfn(100), chfn(100)
      character chfn*200
      character fname*200 ! T.Sato 2022/11/30, increase from 100 to 200 because inumc requests 200 character

      logical lex

      character chaa*7
      real chza
*     parameter ( rimax1    =  4,000 )
*     parameter ( rimax2    = 20,000 )
      double precision act, eng(rimax1), rat(rimax1)
      integer ntyp, icl, kkk, ndk, nsp, nnn

      integer klmax
      data klmax /500000/  !  50,977 lines @ RIsource.dat
                           ! 455,625 lines @ ICRP-07.RAD
      double precision smallv
      data smallv / 1.0d-36 /

*-----------------------------------------------------------------------
*     default values
*-----------------------------------------------------------------------
           m_err = ' Error in rdgamma3 sub.(read RIsource.rad file)'
           ErrCha = ''
           ErrID = 'L:3229/R:rdgamma3/F:risource.f'
           l_err = 1
           k_err = 0
           ierr  = 0
*-----------------------------------------------------------------------
*     file number of scratch file for CHECKING decay calculation
           iot      = 820
*     file number of RIsource.dat or ICRP-07.NDX based on DECDC2
*          iodecdc2 = 811
*     file number of RIsource.rad (radiation data) based on DECDC2
           ioddc2r  = 812
*     file number of RIsource.bet (beta spectra) based on DECDC2
*          ioddc2b  = 813
*     file number of RIsource.ack (Auger-CK spectra) based on DECDC2
*          ioddc2a  = 814
*     file number of RIsource.dat (decay data) formatted with DCHAIN
*          ioddc2s  = 815

c     err flag for this subroutine
           ierr     =   0
*-----------------------------------------------------------------------
*     address check of decay data file
*-----------------------------------------------------------------------

        call fchkd2(24,ntyp,fname,ierr)

          if( ierr .eq. 1 ) then
*           l_err = ill(jsn)
*           k_err = jsn
            go to 998
          end if

        icl = inumc(fname,1,100,' ') - 1
*       do n = 100, 1, -1
*         if( fname(n:n) .ne. ' ' ) go to 201
*       end do
*
* 201     icl = n
          lex = .false.

        if( ntyp .eq. 0 ) then

          fname = fname(1:icl-3)//'rad'
          inquire( file = fname, exist = lex )
          if( .not. lex ) then
            if( icl .gt. 12 ) then
              inquire( file = fname(1:icl-12)//'ICRP-07.RAD ',
     &                 exist = lex )
              if( lex ) then
                fname = fname(1:icl-12)//'ICRP-07.RAD '
              else
                fname = fname(1:icl-3)//'dat'  ! Error
                go to 901
              end if
            else if( icl .eq. 12 ) then
              inquire( file = 'ICRP-07.RAD', exist = lex )
              if( lex ) then
                fname = 'ICRP-07.RAD '
              else
                fname = fname(1:icl-3)//'dat'  ! Error
                go to 901
              end if
            end if
          end if

        else  ! ntyp = 1

          fname = fname(1:icl-3)//'RAD'

        end if

      open( ioddc2r, file = fname,
     &      form='formatted', status = 'unknown' )
        rewind ioddc2r

*-----------------------------------------------------------------------
*     Target nuclide
*-----------------------------------------------------------------------
*         chza = nucl
        call getelm(11,chaa,chza,ierr)

C       Sub. information: getelm(int.,char.,real,int.)
*         OUTPUT:
*          11: chza(Zaid+) -> chaa,7digit(H-1    )

          if( ierr .eq. 1 ) then
*           l_err = ill(jsn)
*           k_err = jsn
*             return
            go to 999
          end if

*-----------------------------------------------------------------------
*     Read the data on the energy and yield of gamma-ray
*-----------------------------------------------------------------------
  220 continue

        nnn = -1

      if( ntyp .eq. 0 .and. fname(icl-2:icl) .eq. 'dat' ) then  ! RIsource.dat

            go to 280  ! pass

      else  ! ntyp = 1 (Extra file: ICRP-07.RAD)
  250 continue

        do n = 1, klmax
          read( ioddc2r,'(a200)', iostat = ios ) chin
*         write(iot,'(i6,'': '',a)') n, chin

          if( n .eq. 1 ) then
            if( ios .eq. -1 ) go to 903
          else
            if( ios .eq. -1 ) go to 280  ! read end
          end if

          if( nnn .eq. 0 ) then
            go to 280  ! successfully end
          else if( nnn .gt. 0 ) then
            nnn = nnn - 1
            if( chin(27:27+3-1) .eq. ' AQ' ) then
                kkk = kkk + 1
              if( kkk .gt. rimax1 ) go to 965
              read(chin(15:15+12-1),'(e15.7)') eng(kkk)
              read(chin(3:3+12-1),'(e15.7)')   rat(kkk)
              if( kkk .eq. 1 ) then
                write(iot,'(2x,1pe15.7,1pe15.7,2x,
     &                    "$ Energy (MeV) and Absolute yield")')
     &                      eng(kkk), rat(kkk)

              else
                write(iot,'(2x,1pe15.7,1pe15.7)') eng(kkk), rat(kkk)

              end if

            end if

          end if

          if( chin(1:1+7-1) .eq. chaa ) then
            read(chin(21:21+9-1),'(i9)') nnn

          end if
        end do

      end if

      if( act .lt. 0.0d+00 ) then
      else
        go to 940
      end if

  280 continue

*-----------------------------------------------------------------------

        goto 999

*-----------------------------------------------------------------------
  901    m_err = 'RIsource.rad does NOT exist.'
         ErrCha = ''
         ErrID = 'L:3390/R:rdgamma3/F:risource.f'
           goto 998
  903    m_err = 'Data in RIsource.rad is not read.'
         ErrCha = ''
         ErrID = 'L:3394/R:rdgamma3/F:risource.f'
           goto 998
  940    m_err = 'Klmax in sub. rdgamma3 is smaller than the data line.'
         ErrCha = ''
         ErrID = 'L:3398/R:rdgamma3/F:risource.f'
           goto 998
  965    m_err = 'Number of gamma-rays was exceeded rimax1.'
         ErrCha = ''
         ErrID = 'L:3402/R:rdgamma3/F:risource.f'
           goto 998
*-----------------------------------------------------------------------
  998 continue
         ierr = 1

  999 continue
*-----------------------------------------------------------------------

      close(ioddc2r)
      return
      end

************************************************************************
*                                                                      *
      subroutine ddc2echo(iot,i,j,ierr)
*                                                                      *
*       for input echo (phits.out)                                     *
*         made by N.Matsuda on 2016/07/31                              *
*         revised by N.Matsuda on 2017/05/29                           *
*                                                                      *
************************************************************************

      implicit double precision (a-h,o-z)
      include 'param.inc'
      include 'risrcparam.inc'

*-----------------------------------------------------------------------

*     common /error/  m_err, l_err, k_err
*     character       m_err*200

      common /isorst/ jstyp(isrc), istyp(isrc), inkf0(isrc), lstyp(isrc)
!$OMP THREADPRIVATE(/isorst/)
      common /ptname/ pname(20), ipln(20)
      character       pname*8

      character chin*200, chlw*200, chcm*200

*     NEW parameters ===
      common /risource/ chalct(nuclinmax,isrc), act000(nuclinmax,isrc),
     &                  activy(nuclinmax,isrc), thf(nuclinmax,isrc),
     &                  decayt(isrc)
      real  chalct  ! ex) 55137.0 ( = Cs-137 )
      double precision  act000, activy, thf, decayt

C MATSUDA 2018.08.15 (icharacterx)
C MATSUDA 2019.05.08 (iannih)
      common /risrc00/ niorg(isrc), nicur(isrc), norm(isrc),
     &                 iaugers(isrc), icharacterx(isrc), iannih(isrc),
     &                 aclow(isrc)
      integer  niorg, nicur, norm, iaugers, icharacterx, iannih
      double precision  aclow

*     IO number for this subroutine ===
      integer  iot, ioddc2
      integer  ierr

      common /paran/  icfn(100), ilfn(100), chfn(100)
      character chfn*200
*     file name for READ
      character fname*200 ! T.Sato 2022/11/30, increase from 100 to 200 because inumc requests 200 character , fname2*100

      character chaa*7
      real chza
      logical lex

      double precision act, rmb, rnum, eng, rat, sumsf, eng2, rat2
      integer ndtype, nnn, ndk, kkk, kk

C     for BETA -/+/both (.BET)  MATSUDA 2017.05.29
      parameter ( nddc2b_minu = 480 )  ! beta minus
      parameter ( nddc2b_plus = 429 )  ! beta plus
      parameter ( nddc2b_both =  22 )  ! both
      common /ddc2beta0/ tbetas00(3,nddc2b_minu)  ! from (.BET)
      parameter ( nddc2bbemax = 120 )  ! Maximum number of energy mesh
*     parameter ( nddc2b_both =  22 )  ! both
      common /ddc2beta3/ ddc2beng(nddc2bbemax), ddc2bend(2*nddc2b_both),
     &                   ddc2byld(2*nddc2b_both,nddc2bbemax),
     &                   ddc2bnuc(nddc2b_both)
      common /ddc2beta4/ nnddc2bn(2*nddc2b_both)

      character ddc2bnuc*7

      integer klmax
      data klmax /500000/

*-----------------------------------------------------------------------
*     file number of scratch file for CHECKING decay calculation
*          iot      = 820
*     file number of RIsource.dat or ICRP-07.NDX based on DECDC2
*           iodecdc2 = 811
*     file number of RIsource.rad (radiation data) based on DECDC2
*          ioddc2r  = 812
*     file number of RIsource.bet (beta spectra) based on DECDC2
*          ioddc2b  = 813
*     file number of RIsource.ack (Auger-CK spectra) based on DECDC2
*          ioddc2a  = 814
*     file number of RIsource.dat (decay data) formatted with DCHAIN
*          ioddc2s  = 815

c     err flag for this subroutine
           ierr     =   0
*-----------------------------------------------------------------------
           ioddc2   = 811
*-----------------------------------------------------------------------
*     Nuclide
*-----------------------------------------------------------------------
        chza = chalct(i,j)

      if( chza .gt. 990000.0 ) then  ! Spontaneous fission nuclide
        sumsf = activy(i,j)  ! (atoms)
          return

      end if

        call getelm(11,chaa,chza,ierr)

C       Sub. information: getelm(int.,char.,real,int.)
*       you can choose INPUT or OUTPUT style
*         OUTPUT:
*          11: chza(Zaid+) -> chaa,7digit(H-1    )

*         if( ierr .eq. 1 ) then
*           l_err = ill(jsn)
*           k_err = jsn
*           go to 998
*         end if

      if( thf(i,j) .eq. -1.0d+00 ) then  ! Stable nuclide
        if( i .le. niorg(j) ) then
          write(iot,'("    ",a7,16x,
     &                      " # Stable nuclide")') chaa
        else
          write(iot,'("$   ",a7,16x,
     &                      " # Stable nuclide")') chaa
        end if

      else

*-----------------------------------------------------------------------
*     Decay constant and Activity
*-----------------------------------------------------------------------
        rmb = dble( log(2.) / thf(i,j) )
        act = activy(i,j) * rmb  ! (Bq)

*-----------------------------------------------------------------------
        if( i .le. niorg(j) ) then
          if( act .gt. aclow(j) ) then
            write(iot,advance='no',fmt='(" ",3x,a7,1pe13.5,4x,
     &          "# (Bq)"/28x,"# ->",1pe13.5," (Bq)")')
     &            chaa, act000(i,j), act
          else
            write(iot,advance='no',fmt='(" ",3x,a7,1pe13.5,4x,
     &          "# (Bq)"/28x,"# -> --- L.L. --- (Bq)")')
     &            chaa, act000(i,j)
*             return
          end if
        else
          if( act .gt. aclow(j) ) then
            write(iot,advance='no',fmt='("$",3x,a7,1pe13.5,4x,
     &          "# daughter"/28x,"# ->",1pe13.5," (Bq)")')
     &            chaa, 0.0d+00,     act
          else
            write(iot,advance='no',fmt='("$",3x,a7,1pe13.5,4x,
     &          "# daughter"/28x,"# -> --- L.L. --- (Bq)")')
     &            chaa, 0.0d+00
*             return
          end if
        end if

C     Half life
            write(iot,advance='no',fmt='(", half life:")')
        if( thf(i,j) .lt. 100. ) then
            write(iot,fmt='(x,f12.8,x,"(sec)")')
     &            thf(i,j)
        else if( thf(i,j)/60. .lt. 100. ) then
            write(iot,fmt='(x,f12.8,x,"(min.)")')
     &            thf(i,j)/60.
        else if( thf(i,j)/(60.*60.) .lt. 50. ) then
            write(iot,fmt='(x,f12.8,x,"(hour)")')
     &            thf(i,j)/(60.*60.)
        else if( thf(i,j)/(60.*60.*24.) .lt. 300. ) then
            write(iot,fmt='(f13.8,x,"(day)")')
     &            thf(i,j)/(60.*60.*24.)
        else if( thf(i,j)/(60.*60.*24.*365.26) .lt. 100. ) then
            write(iot,fmt='(x,f12.8,x,"(year)")')
     &            thf(i,j)/(60.*60.*24.*365.26)
        else
            write(iot,fmt='(1pe13.4,x,"(year)")')
     &            thf(i,j)/(60.*60.*24.*365.26)
        end if

        if( act .lt. aclow(j) ) return

*-----------------------------------------------------------------------
*     address check of decay data file
*-----------------------------------------------------------------------
  200 continue

        call fchkd2(24,ndtype,fname,ierr)

*         if( ierr .eq. 1 ) then
*           l_err = ill(jsn)
*           k_err = jsn
*           go to 998
*         end if

        do n = 100, 1, -1
          if( fname(n:n) .ne. ' ' ) go to 201
        end do

  201     icl = n
          lex = .false.

        if( inkf0(j) .eq. 22 ) then  ! photon
          if( ndtype .eq. 0 ) then
            inquire( file = fname(1:icl-3)//'rad', exist = lex )
            if( lex ) then    ! RIsource.rad
                fname = fname(1:icl-3)//'rad'
                ioddc2 = 812
            else
              if( icl .gt. 12 ) then
                inquire( file = fname(1:icl-12)//'ICRP-07.RAD ',
     &                   exist = lex )
                if( lex ) then  ! ICRP-07.RAD
                  fname = fname(1:icl-12)//'ICRP-07.RAD '
                  ioddc2 = 812
                else            ! RIsource.dat
                  fname = fname(1:icl-3)//'dat'
                  ioddc2 = 815
                end if
              else if( icl .eq. 12 ) then
                inquire( file = 'ICRP-07.RAD', exist = lex )
                if( lex ) then  ! ICRP-07.RAD
                  fname = 'ICRP-07.RAD '
                  ioddc2 = 812
                else            ! RIsource.dat
                  fname = fname(1:icl-3)//'dat'
                  ioddc2 = 815
                end if
              end if
            end if

          else                ! ICRP-07.RAD ( ndtype = 1 )
                fname = fname(1:icl-3)//'RAD'
                ioddc2 = 812

          end if

        else if( inkf0(j) .eq. 2000004 ) then  ! alpha particle
          if( ndtype .eq. 0 ) then
            inquire( file = fname(1:icl-3)//'rad', exist = lex )
            if( lex ) then    ! RIsource.rad
                fname = fname(1:icl-3)//'rad'
                ioddc2 = 811
            else
              if( icl .gt. 12 ) then
                inquire( file = fname(1:icl-12)//'ICRP-07.RAD ',
     &                   exist = lex )
                if( lex ) then  ! ICRP-07.RAD
                  fname = fname(1:icl-12)//'ICRP-07.RAD '
                  ioddc2 = 811
                end if
              else if( icl .eq. 12 ) then
                inquire( file = 'ICRP-07.RAD', exist = lex )
                if( lex ) then  ! ICRP-07.RAD
                  fname = 'ICRP-07.RAD '
                  ioddc2 = 811
                end if
              end if
            end if

          else                ! ICRP-07.RAD ( ndtype = 1 )
                fname  = fname(1:icl-3)//'RAD'
                ioddc2 = 811

          end if

        else if( inkf0(j) .eq. 11 .or. inkf0(j) .eq. -11 ) then  ! electron or positron
          if( ndtype .eq. 0 ) then
            inquire( file = fname(1:icl-3)//'bet', exist = lex )
            if( lex ) then    ! RIsource.bet
                fname = fname(1:icl-3)//'bet'
                ioddc2 = 813
            else
              if( icl .gt. 12 ) then
                inquire( file = fname(1:icl-12)//'ICRP-07.BET ',
     &                   exist = lex )
                if( lex ) then  ! ICRP-07.BET
                  fname = fname(1:icl-12)//'ICRP-07.BET '
                  ioddc2 = 813
                end if
              else if( icl .eq. 12 ) then
                inquire( file = 'ICRP-07.BET', exist = lex )
                if( lex ) then  ! ICRP-07.BET
                  fname = 'ICRP-07.BET '
                  ioddc2 = 813
                end if
              end if
            end if

          else                ! ICRP-07.BET ( ndtype = 1 )
                fname  = fname(1:icl-3)//'BET'
                ioddc2 = 813

          end if

        else if( inkf0(j) .eq. 2112 ) then
          if( ndtype .eq. 0 ) then
            inquire( file = fname(1:icl-3)//'nsf', exist = lex )
            if( lex ) then
              fname = fname(1:icl-3)//'nsf'
              ioddc2 = 816
            else
              if( icl .gt. 12 ) then
                inquire( file = fname(1:icl-12)//'ICRP-07.NSF ',
     &                   exist = lex )
                if( lex ) then
                  fname = fname(1:icl-12)//'ICRP-07.NSF '
                  ioddc2 = 816
                end if
              else if( icl .eq. 12 ) then
                inquire( file = 'ICRP-07.NSF', exist = lex )
                if( lex ) then
                  fname = 'ICRP-07.NSF '
                  ioddc2 = 816
                end if
              end if
            end if
          else
            fname  = fname(1:icl-3)//'NSF'
            ioddc2 = 816

          end if

        end if
*       Alpha: 811, Beta: 813, Gamma: 812 or 815 (only RIsource.dat)

  210   continue
          open( ioddc2, file = fname,
     &                  form = 'formatted', status = 'unknown' )
          rewind ioddc2

*-----------------------------------------------------------------------
*     Read the data and write
*-----------------------------------------------------------------------
  220   continue
          nnn = -1


        if( ioddc2 .eq. 815 ) then  ! RIsource.dat
          call getelm(14,chaa,chza,ierr)

C       Sub. information: getelm(int.,char.,real,int.)
*         OUTPUT:
*          14: chza(Zaid+) -> chaa,7digit(%H   1 )

*           if( ierr .eq. 1 ) then
*             l_err = ill(jsn)
*             k_err = jsn
*             go to 999
*           end if

          do n = 1, klmax
            read( ioddc2,'(a200)', iostat = ios ) chin
*           write(iot,'(i6,''. '',a)') n, chin

            if( ios .eq. -1 ) go to 500  ! end sub.

            if( nnn .eq.  0 ) then
              go to 500  ! end sub. successfully

            else if( nnn .gt.  0 ) then
              do nn = 1, 3
                  nnn = nnn - 1
                read(chin(1+12*(2*nn-2):12*(2*nn-1)),'(e12.5)') eng
                read(chin(1+12*(2*nn-1):12*(2*nn)),'(e12.5)')   rat

                write(iot,'("$ ",1pe13.5,1pe13.5,x,1pe13.5,"*"
     &                            ,1pe11.5)') eng, eng, act, abs(rat)
                if( nnn .eq. 0 ) exit
              end do

            end if

            if( chin(1:1) .eq. '%' ) then
              if( chin(3:3+6-1) .eq. chaa(2:7) ) then
                read( ioddc2, '()' )
                read( ioddc2,'(a200)', iostat = ios ) chin
                read(chin(7:7+6-1),*) nnn

                if( nnn .eq. 0 ) go to 500  ! end sub.

                read(chin(1:1+6-1),*) ndk

                if( ndk .eq. 0 ) go to 500  ! end sub.

                do nn = 1, ndk
                  read( ioddc2, '()' )
                end do

                read( ioddc2,'(a200)', iostat = ios ) chin
                read(chin(7:7+6-1),*) nnn

                if( nnn .eq. 0 ) go to 500  ! end sub.

                  write(iot,'("$   Gamma-rays with annihilation",/
     &                        "$   Energy (MeV/n)",13x,
     &                        "Activity (Bq) X Yield (abs.)",/
     &                        "$      Lower        Upper")')

              else
                read( ioddc2, '()' )
                read( ioddc2, '()' )

              end if
            end if

          end do

        else if( ioddc2 .eq. 811 .or. ioddc2 .eq. 812 ) then  ! ***.rad

C MATSUDA 2018.08.15  for character X-rays
          IF( ioddc2 .eq. 812 .and.
     &      ( icharacterx(j) .eq. 0 .or. icharacterx(j) .eq. 2 ) ) then

          do n = 1, klmax
            read( ioddc2,'(a200)', iostat = ios ) chin
*           write(iot,'(i6,''. '',a)') n, chin

            if( ios .eq. -1 ) go to 240

            if( nnn .eq.  0 ) then
              go to 240

            else if( nnn .gt.  0 ) then
              nnn = nnn - 1

              if( ioddc2 .eq. 812 ) then  ! X-rays
                if( chin(27:27+3-1) .eq. '  X' ) then
                    kkk = kkk + 1
                  read(chin(15:15+12-1),'(e15.7)') eng
                  read(chin(3:3+12-1),'(e15.7)')   rat

                  if( kkk .eq. 1 ) then
                    write(iot,'("$   X-rays (X)",/
     &                          "$   Energy (MeV/n)",13x,
     &                          "Activity (Bq) X Yield (abs.)",/
     &                          "$      Lower        Upper")')
                  end if

                  write(iot,'("$ ",1pe13.5,1pe13.5,x,1pe13.5,"*"
     &                              ,1pe11.5)') eng, eng, act, abs(rat)

                end if

              end if

            end if

            if( chin(1:1+7-1) .eq. chaa ) then
                read(chin(21:21+9-1),*) nnn
                  kkk =  0

            end if

          end do
              rewind ioddc2

          END IF

  240 continue
          rewind ioddc2
C MATSUDA 2019.06.20  BUGs
          if( ioddc2 .eq. 812 .and. icharacterx(j) .eq. 2 ) go to 245
              nnn = -1

          do n = 1, klmax
            read( ioddc2,'(a200)', iostat = ios ) chin
*           write(iot,'(i6,''. '',a)') n, chin

            if( ios .eq. -1 ) go to 245

            if( nnn .eq.  0 ) then
              go to 245

            else if( nnn .gt.  0 ) then
              nnn = nnn - 1

              if( ioddc2 .eq. 812 ) then  ! gamma-rays
                if( chin(27:27+3-1) .eq. '  G' ) then
                    kkk = kkk + 1
                  read(chin(15:15+12-1),'(e15.7)') eng
                  read(chin(3:3+12-1),'(e15.7)')   rat

                  if( kkk .eq. 1 ) then
                    write(iot,'("$   Gamma-rays (G)",/
     &                          "$   Energy (MeV/n)",13x,
     &                          "Activity (Bq) X Yield (abs.)",/
     &                          "$      Lower        Upper")')
                  end if

                  write(iot,'("$ ",1pe13.5,1pe13.5,x,1pe13.5,"*"
     &                              ,1pe11.5)') eng, eng, act, abs(rat)

                end if

              else if( ioddc2 .eq. 811 ) then  ! alpha particles
                if( chin(27:27+3-1) .eq. '  A' ) then
                    kkk = kkk + 1
                  read(chin(15:15+12-1),'(e15.7)') eng
                  read(chin(3:3+12-1),'(e15.7)')   rat

                  if( kkk .eq. 1 ) then
                    write(iot,'("$   Alpha particles (A)",/
     &                          "$   Energy (MeV/n)",13x,
     &                          "Activity (Bq) X Yield (abs.)",/
     &                          "$      Lower        Upper")')
                  end if

C MATSUDA 2017.05.29 ( MeV -> MeV/n for alpha )
                  write(iot,'("$ ",1pe13.5,1pe13.5,x,1pe13.5,"*"
     &                              ,1pe11.5)') eng/4.0d+00, eng/4.0d+00
     &                              , act, abs(rat)

                end if

              end if

            end if

            if( chin(1:1+7-1) .eq. chaa ) then
                read(chin(21:21+9-1),*) nnn
                  kkk =  0

            end if

          end do

C MATSUDA 2019.05.08  for annihilation photons
  245 continue
          rewind ioddc2
              nnn = -1

          IF( ioddc2 .eq. 812 .and. iannih(j) .eq. 0 ) then

          do n = 1, klmax
            read( ioddc2,'(a200)', iostat = ios ) chin
*           write(iot,'(i6,''. '',a)') n, chin

            if( ios .eq. -1 ) go to 500  ! end sub.

            if( nnn .eq.  0 ) then
              go to 500  ! end sub. successfully

            else if( nnn .gt.  0 ) then
              nnn = nnn - 1

              if( ioddc2 .eq. 812 ) then  ! annihilation photons
                if( chin(27:27+3-1) .eq. ' AQ' ) then
                    kkk = kkk + 1
                  read(chin(15:15+12-1),'(e15.7)') eng
                  read(chin(3:3+12-1),'(e15.7)')   rat

                  if( kkk .eq. 1 ) then
                    write(iot,'("$   Annihilation photons (AQ)",/
     &                          "$   Energy (MeV/n)",13x,
     &                          "Activity (Bq) X Yield (abs.)",/
     &                          "$      Lower        Upper")')
                  end if

                  write(iot,'("$ ",1pe13.5,1pe13.5,x,1pe13.5,"*"
     &                              ,1pe11.5)') eng, eng, act, abs(rat)

                end if

              end if

            end if

            if( chin(1:1+7-1) .eq. chaa ) then
                read(chin(21:21+9-1),*) nnn
                  kkk =  0

            end if

          end do

          END IF

        else if( ioddc2 .eq. 813 ) then  ! electron or positron
          if( inkf0(j) .eq. 11 .and. iaugers(j) .eq. 2 ) go to 280

*         do  n = 1, nnbetas(3)   ! both
          do  n = 1, nddc2b_both  ! both
*           if( tbetas00(3,n) .eq. chza ) then
            if( abs(tbetas00(3,n)-chza) .lt. 0.01 ) then

              if( inkf0(j) .eq.  11 ) then
                kk =  n
                    write(iot,'("$   Electron spectrum",/
     &                          "$   Energy (MeV/n)",13x,
     &                          "Activity (Bq) X Yield (abs.)",/
     &                          "$      Lower        Upper")')
              else if( inkf0(j) .eq. -11 ) then
                kk =  nddc2b_both + n
                    write(iot,'("$   Positron spectrum",/
     &                          "$   Energy (MeV/n)",13x,
     &                          "Activity (Bq) X Yield (abs.)",/
     &                          "$      Lower        Upper")')
              end if

                eng2 = ddc2beng(1)
                rat2 = ddc2byld(kk,1)
              do nn = 2, nnddc2bn(kk)
                eng  = eng2
                rat  = rat2
                eng2 = ddc2beng(nn)
                rat2 = ddc2byld(kk,nn)

                  write(iot,'("$ ",1pe13.5,1pe13.5,x,1pe13.5,"*"
     &                              ,1pe11.5)') eng, eng2, act
     &                              , abs((rat+rat2)/2.*(eng2-eng))

              end do

                eng  = eng2
                rat  = rat2
                eng2 = ddc2bend(kk)
                rat2 = 0.0d+00

                  write(iot,'("$ ",1pe13.5,1pe13.5,x,1pe13.5,"*"
     &                              ,1pe11.5)') eng, eng2, act
     &                              , abs((rat+rat2)/2.*(eng2-eng))

              if( inkf0(j) .eq.  11 ) go to 280
*             if( inkf0(j) .eq. -11 ) go to 500
                  go to 500  ! end sub. successfully

            end if

          end do

  250 continue
          if( inkf0(j) .eq.  11 ) then

*           do  n = 1, nnbetas(1)   ! electron (beta minus)
            do  n = 1, nddc2b_minu  ! electron (beta minus)
*             if( tbetas00(1,n) .eq. chza ) then
              if( abs(tbetas00(1,n)-chza) .lt. 0.01 ) then
                  go to 260

              end if
            end do
                  go to 280

          else if( inkf0(j) .eq. -11 ) then

*           do  n = 1, nnbetas(2)   ! positron (beta plus)
            do  n = 1, nddc2b_plus  ! positron (beta plus)
*             if( tbetas00(2,n) .eq. chza ) then
              if( abs(tbetas00(2,n)-chza) .lt. 0.01 ) then
                  go to 260

              end if
            end do
                  go to 500  ! end sub. successfully

          end if

  260 continue
              nnn = -1

          do n = 1, klmax
            read( ioddc2,'(a200)', iostat = ios ) chin
*           write(iot,'(i6,''. '',a)') n, chin

            if( ios .eq. -1 ) then
              if( inkf0(j) .eq.  11 ) go to 280
              if( inkf0(j) .eq. -11 ) go to 500  ! end sub.

            end if
            if( nnn .eq.  0 ) then
              if( inkf0(j) .eq.  11 ) go to 280
              if( inkf0(j) .eq. -11 ) go to 500  ! end sub. successfully

            else if( nnn .gt.  0 ) then
              nnn = nnn - 1
              kkk = kkk + 1

                eng  = eng2
                rat  = rat2
                read(chin(1:1+7-1),'(e15.7)')  eng2
                read(chin(8:8+10-1),'(e15.7)') rat2
                if( kkk .eq. 1 ) then
                  if( inkf0(j) .eq.  11 ) then
                    write(iot,'("$   Electron spectrum",/
     &                          "$   Energy (MeV/n)",13x,
     &                          "Activity (Bq) X Yield (abs.)",/
     &                          "$      Lower        Upper")')
                  else if( inkf0(j) .eq. -11 ) then
                    write(iot,'("$   Positron spectrum",/
     &                          "$   Energy (MeV/n)",13x,
     &                          "Activity (Bq) X Yield (abs.)",/
     &                          "$      Lower        Upper")')
                  end if

                else if( kkk .gt. 2 ) then
                  write(iot,'("$ ",1pe13.5,1pe13.5,x,1pe13.5,"*"
     &                              ,1pe11.5)') eng, eng2, act
     &                              , abs((rat+rat2)/2.*(eng2-eng))

                end if

            end if

            if( chin(1:1+7-1) .eq. chaa ) then
                read(chin(8:8+10-1),*) nnn
                eng2 = 0.0d+00
                rat2 = 0.0d+00
              kkk =  0

            end if

          end do

  280 continue  ! Auger-CK electron spectrum
          if( iaugers(j) .eq. 1 ) go to 500  ! end sub. successfully

            close(ioddc2)

          if( ndtype .eq. 0 ) then
            inquire( file = fname(1:icl-3)//'ack', exist = lex )
            if( lex ) then    ! RIsource.ack
                fname = fname(1:icl-3)//'ack'
                ioddc2 = 814
            else
              if( icl .gt. 12 ) then
                inquire( file = fname(1:icl-12)//'ICRP-07.ACK ',
     &                   exist = lex )
                if( lex ) then  ! ICRP-07.ACK
                  fname = fname(1:icl-12)//'ICRP-07.ACK '
                  ioddc2 = 814
                end if
              else if( icl .eq. 12 ) then
                inquire( file = 'ICRP-07.ACK ', exist = lex )
                if( lex ) then  ! ICRP-07.ACK
                  fname = 'ICRP-07.ACK '
                  ioddc2 = 814
                end if
              end if
            end if

          else                ! ICRP-07.ACK ( ndtype = 1 )
                fname  = fname(1:icl-3)//'ACK'
                ioddc2 = 814

          end if

          open( ioddc2, file = fname,
     &                  form = 'formatted', status = 'unknown' )
          rewind ioddc2

              nnn = -1

          do n = 1, klmax
            read( ioddc2,'(a200)', iostat = ios ) chin
*           write(iot,'(i7,''. '',a)') n, chin

*           if( ios .eq. -1 ) go to 500  ! end sub.
            if( ios .eq. -1 ) go to 290  !

            if( nnn .eq.  0 ) then
              go to 300  ! to internal conversion electron

            else if( nnn .gt.  0 ) then
              nnn = nnn - 1
              kkk = kkk + 1

                read(chin(12:12+12-1),'(e15.7)') eng
                read(chin(1:1+11-1),'(e15.7)')   rat

                  write(iot,'("$ ",1pe13.5,1pe13.5,x,1pe13.5,"*"
     &                              ,1pe11.5)') eng/1.0d+06, eng/1.0d+06
     &                              , act, abs(rat)

            end if

            if( chin(1:1+7-1) .eq. chaa ) then
                read(chin(22:22+11-1),*) nnn
                  write(iot,'("$   Auger-CK electron spectrum",/
     &                        "$   Energy (MeV/n)",13x,
     &                        "Activity (Bq) X Yield (abs.)",/
     &                        "$      Lower        Upper")')
              kkk =  0

            end if

          end do

  290 continue  ! Auger electron in RIsource.rad

            close(ioddc2)

          if( ndtype .eq. 0 ) then
            inquire( file = fname(1:icl-3)//'rad', exist = lex )
            if( lex ) then    ! RIsource.rad
                fname = fname(1:icl-3)//'rad'
                ioddc2 = 812
            else
              if( icl .gt. 12 ) then
                inquire( file = fname(1:icl-12)//'ICRP-07.RAD ',
     &                   exist = lex )
                if( lex ) then  ! ICRP-07.RAD
                  fname = fname(1:icl-12)//'ICRP-07.RAD '
                  ioddc2 = 812
                end if
              else if( icl .eq. 12 ) then
                inquire( file = 'ICRP-07.RAD ', exist = lex )
                if( lex ) then  ! ICRP-07.RAD
                  fname = 'ICRP-07.RAD '
                  ioddc2 = 812
                end if
              end if
            end if

          else                ! ICRP-07.RAD ( ndtype = 1 )
                fname  = fname(1:icl-3)//'RAD'
                ioddc2 = 812

          end if

          open( ioddc2, file = fname,
     &                  form = 'formatted', status = 'unknown' )
          rewind ioddc2

              nnn = -1

          do n = 1, klmax
            read( ioddc2,'(a200)', iostat = ios ) chin
*           write(iot,'(i7,''. '',a)') n, chin

            if( ios .eq. -1 ) go to 300  ! to internal conversion electron

            if( nnn .eq.  0 ) then
              go to 300  ! to internal conversion electron

            else if( nnn .gt.  0 ) then
              nnn = nnn - 1
                if( chin(27:27+3-1) .eq. ' AE' ) then
                    kkk = kkk + 1
                  read(chin(15:15+12-1),'(e15.7)') eng
                  read(chin(3:3+12-1),'(e15.7)')   rat

                  if( kkk .eq. 1 ) then
                    write(iot,'("$   Auger electron spectrum (AE)",/
     &                          "$   Energy (MeV/n)",13x,
     &                          "Activity (Bq) X Yield (abs.)",/
     &                          "$      Lower        Upper")')
                  end if

                  write(iot,'("$ ",1pe13.5,1pe13.5,x,1pe13.5,"*"
     &                              ,1pe11.5)') eng, eng, act, abs(rat)

                end if

            end if

            if( chin(1:1+7-1) .eq. chaa ) then
                read(chin(21:21+9-1),*) nnn
              kkk =  0

            end if

          end do

  300 continue  ! Internal conversion electron in RIsource.rad

            close(ioddc2)

          if( ndtype .eq. 0 ) then
            inquire( file = fname(1:icl-3)//'rad', exist = lex )
            if( lex ) then    ! RIsource.rad
                fname = fname(1:icl-3)//'rad'
                ioddc2 = 812
            else
              if( icl .gt. 12 ) then
                inquire( file = fname(1:icl-12)//'ICRP-07.RAD ',
     &                   exist = lex )
                if( lex ) then  ! ICRP-07.RAD
                  fname = fname(1:icl-12)//'ICRP-07.RAD '
                  ioddc2 = 812
                end if
              else if( icl .eq. 12 ) then
                inquire( file = 'ICRP-07.RAD ', exist = lex )
                if( lex ) then  ! ICRP-07.RAD
                  fname = 'ICRP-07.RAD '
                  ioddc2 = 812
                end if
              end if
            end if

          else                ! ICRP-07.RAD ( ndtype = 1 )
                fname  = fname(1:icl-3)//'RAD'
                ioddc2 = 812

          end if

          open( ioddc2, file = fname,
     &                  form = 'formatted', status = 'unknown' )
          rewind ioddc2

              nnn = -1

          do n = 1, klmax
            read( ioddc2,'(a200)', iostat = ios ) chin
*           write(iot,'(i7,''. '',a)') n, chin

            if( ios .eq. -1 ) go to 500  ! end sub.

            if( nnn .eq.  0 ) then
              go to 500  ! end sub. successfully

            else if( nnn .gt.  0 ) then
              nnn = nnn - 1
                if( chin(27:27+3-1) .eq. ' IE' ) then
                    kkk = kkk + 1
                  read(chin(15:15+12-1),'(e15.7)') eng
                  read(chin(3:3+12-1),'(e15.7)')   rat

                  if( kkk .eq. 1 ) then
                    write(iot,'("$   Internal conversion electron"
     &                          " (IE)",/
     &                          "$   Energy (MeV/n)",13x,
     &                          "Activity (Bq) X Yield (abs.)",/
     &                          "$      Lower        Upper")')
                  end if

                  write(iot,'("$ ",1pe13.5,1pe13.5,x,1pe13.5,"*"
     &                              ,1pe11.5)') eng, eng, act, abs(rat)

                end if

            end if

            if( chin(1:1+7-1) .eq. chaa ) then
                read(chin(21:21+9-1),*) nnn
              kkk =  0

            end if

          end do

        else if( ioddc2 .eq. 816 ) then  ! SF neutron
          do n = 1, klmax
            read( ioddc2,'(a200)', iostat = ios ) chin

            if( ios .eq. -1 ) go to 500

            if( nnn .eq.  0 ) then
              go to 500

            else if( nnn .gt.  0 ) then
              nnn = nnn - 1
              kkk = kkk + 1
              read(chin(1:1+8-1),'(e8.2)') eng_lw
              read(chin(10:10+8-1),'(e8.2)') eng_up
              read(chin(19:19+11-1),'(e11.5)') rat
              if( kkk .eq. 1 ) then
                write(iot,'("$   Spontaneous fission neutron (SF)",/
     &                      "$   Energy (MeV/n)",13x,
     &                      "Activity (Bq) X Yield (abs.)",/
     &                      "$      Lower        Upper")')
              end if

                write(iot,'("$ ",1pe13.5,1pe13.5,x,1pe13.5,"*"
     &          ,1pe11.5)') eng_lw, eng_up, act, abs(rat)

            end if

            if( chin(1:1+7-1) .eq. chaa ) then
              read(chin(27:27+3-1),'(i3)') nnn
              kkk =  0

            end if

          end do
C
        end if

      end if

*-----------------------------------------------------------------------

  500 continue
*-----------------------------------------------------------------------

        close(ioddc2)
        go to 999

*-----------------------------------------------------------------------
  998 continue
         ierr = 1

  999 continue
*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      subroutine rdbetaS(ntyp,chza,act,kkk,eng,rat,ierr)
*                                                                      *
*       read RIsource.bet file (energy and yield)                      *
*         made by N.Matsuda on 2017/05/29                              *
*                                                                      *
************************************************************************

      implicit double precision (a-h,o-z)
      include 'risrcparam.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      common /error/  m_err, l_err, k_err
      character       m_err*200

      character chin*200, chlw*200, chcm*200

      integer  iot, iodecdc2, ioddc2r, ioddc2b, ioddc2a, ioddc2s
      integer ierr

      common /paran/  icfn(100), ilfn(100), chfn(100)
      character chfn*200
      character fname*200 ! T.Sato 2022/11/30, increase from 100 to 200 because inumc requests 200 character

      logical lex

      character chaa*7
      real chza
*     parameter ( rimax1    =  4,000 )
*     parameter ( rimax2    = 20,000 )
      double precision act, eng(rimax1), rat(rimax1)
      integer ntyp, icl, kkk

      integer klmax
      data klmax /150000/  ! 111,325 lines @ ICRP-07.BET
      double precision smallv
      data smallv / 1.0d-36 /

*-----------------------------------------------------------------------
*     default values
*-----------------------------------------------------------------------
           m_err = ' Error in rdbetaS sub.(read RIsource.bet file)'
           ErrCha = ''
           ErrID = 'L:4459/R:rdbetaS/F:risource.f'
           l_err = 1
           k_err = 0
           ierr  = 0
*-----------------------------------------------------------------------
*     file number of scratch file for CHECKING decay calculation
           iot      = 820
*     file number of RIsource.dat or ICRP-07.NDX based on DECDC2
*          iodecdc2 = 811
*     file number of RIsource.rad (radiation data) based on DECDC2
*          ioddc2r  = 812
*     file number of RIsource.bet (beta spectra) based on DECDC2
           ioddc2b  = 813
*     file number of RIsource.ack (Auger-CK spectra) based on DECDC2
*          ioddc2a  = 814
*     file number of RIsource.dat (decay data) formatted with DCHAIN
*          ioddc2s  = 815

c     err flag for this subroutine
           ierr     =   0
*-----------------------------------------------------------------------
*     address check of decay data file
*-----------------------------------------------------------------------

        call fchkd2(24,ntyp,fname,ierr)

          if( ierr .eq. 1 ) then
*           l_err = ill(jsn)
*           k_err = jsn
            go to 998
          end if

        icl = inumc(fname,1,100,' ') - 1
*       do n = 100, 1, -1
*         if( fname(n:n) .ne. ' ' ) go to 201
*       end do
*
* 201     icl = n
          lex = .false.

        if( ntyp .eq. 0 ) then

          fname = fname(1:icl-3)//'bet'
          inquire( file = fname, exist = lex )
          if( .not. lex ) then
            if( icl .gt. 12 ) then
              inquire( file = fname(1:icl-12)//'ICRP-07.BET ',
     &                 exist = lex )
              if( lex ) then
                fname = fname(1:icl-12)//'ICRP-07.BET '
              end if
            else if( icl .eq. 12 ) then
              inquire( file = 'ICRP-07.BET', exist = lex )
              if( lex ) then
                fname = 'ICRP-07.BET '
              end if
            end if
          end if

        else  ! ntyp = 1

          fname = fname(1:icl-3)//'BET'

        end if

      open( ioddc2b, file = fname,
     &      form='formatted', status = 'unknown' )
        rewind ioddc2b

*-----------------------------------------------------------------------
*     Target nuclide
*-----------------------------------------------------------------------
*         chza = nucl
        call getelm(11,chaa,chza,ierr)

C       Sub. information: getelm(int.,char.,real,int.)
*         OUTPUT:
*          11: chza(Zaid+) -> chaa,7digit(H-1    )

          if( ierr .eq. 1 ) then
*           l_err = ill(jsn)
*           k_err = jsn
*             return
            go to 999
          end if

*-----------------------------------------------------------------------
*     Read the data on the energy and yield of electron or positron
*-----------------------------------------------------------------------
  220 continue

        nnn = -1
        kkk =  0

c     else  ! ntyp = 1 (Extra file: ICRP-07.BET)

        do n = 1, klmax
          read( ioddc2b,'(a200)', iostat = ios ) chin
*         write(iot,'(i7,'': '',a)') n, chin

          if( n .eq. 1 ) then
            if( ios .eq. -1 ) go to 903
          else
            if( ios .eq. -1 ) go to 280  ! read end
          end if

          if( nnn .eq. 0 ) then
            write(iot, '(2x,1pe15.7)') eng(kkk)
            go to 280  ! successfully end
          else if( nnn .gt. 0 ) then
            nnn = nnn - 1
            kkk = kkk + 1
              read(chin(1:1+7-1),'(e15.7)')  eng(kkk)
              read(chin(8:8+10-1),'(e15.7)') rat(kkk)
              if( kkk .eq. 1 ) then

              else if( kkk .eq. 2 ) then
                rat(kkk-1) = ( rat(kkk-1) + rat(kkk) ) / 2.0d+00
     &                     * ( eng(kkk) - eng(kkk-1) )
                write(iot,'(2x,1pe15.7,1pe15.7,2x,
     &                    "$ Energy (MeV) and Absolute yield")')
     &                      eng(kkk-1), rat(kkk-1)

              else
                rat(kkk-1) = ( rat(kkk-1) + rat(kkk) ) / 2.0d+00
     &                     * ( eng(kkk) - eng(kkk-1) )
                write(iot,'(2x,1pe15.7,1pe15.7)') eng(kkk-1), rat(kkk-1)

              end if

          end if

          if( chin(1:1+7-1) .eq. chaa ) then
            read(chin(8:8+10-1),'(i10)') nnn
            if( nnn .gt. rimax1 ) go to 965

          end if
        end do


      if( act .lt. 0.0d+00 ) then
      else
        go to 940
      end if

  280 continue

*-----------------------------------------------------------------------

        goto 999

*-----------------------------------------------------------------------
  901    m_err = 'RIsource.bet does NOT exist.'
         ErrCha = ''
         ErrID = 'L:4613/R:rdbetaS/F:risource.f'
           goto 998
  903    m_err = 'Data in RIsource.bet is not read.'
         ErrCha = ''
         ErrID = 'L:4617/R:rdbetaS/F:risource.f'
           goto 998
  925    m_err = 'RI source for electron need DECDC2(ICRP-07.BET) file.'
         ErrCha = ''
         ErrID = 'L:4621/R:rdbetaS/F:risource.f'
           goto 998
         ierr = 1
           return
  940    m_err = 'Klmax in sub. rdbetaS is smaller than the data line.'
         ErrCha = ''
         ErrID = 'L:4627/R:rdbetaS/F:risource.f'
           goto 998

  965    m_err = 'Number of electron/positron was exceeded rimax1.'
         ErrCha = ''
         ErrID = 'L:4632/R:rdbetaS/F:risource.f'
           goto 998
*-----------------------------------------------------------------------
  998 continue
         ierr = 1

  999 continue
*-----------------------------------------------------------------------

      close(ioddc2b)
      return
      end

************************************************************************
*                                                                      *
      subroutine rdbetaD(ntyp,chza,act,kkk,eng,rat,ierr)
*                                                                      *
*       read RIsource.ack file (energy and yield)                      *
*                                               for Auger-CK electron  *
*         made by N.Matsuda on 2017/05/29                              *
*                                                                      *
************************************************************************

      implicit double precision (a-h,o-z)
      include 'risrcparam.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      common /error/  m_err, l_err, k_err
      character       m_err*200

      character chin*200, chlw*200, chcm*200

      integer  iot, iodecdc2, ioddc2r, ioddc2b, ioddc2a, ioddc2s
      integer ierr

      common /paran/  icfn(100), ilfn(100), chfn(100)
      character chfn*200
      character fname*200 ! T.Sato 2022/11/30, increase from 100 to 200 because inumc requests 200 character

      logical lex

      character chaa*7
      real chza
*     parameter ( rimax1    =  4,000 )
*     parameter ( rimax2    = 20,000 )
      double precision act, eng(rimax1), rat(rimax1)
      integer ntyp, icl, kkk

      integer klmax
      data klmax /150000/  ! 131,571 lines @ ICRP-07.ACK
      double precision smallv
      data smallv / 1.0d-36 /

*-----------------------------------------------------------------------
*     default values
*-----------------------------------------------------------------------
           m_err = ' Error in rdbetaD sub.(read RIsource.ack file)'
           ErrCha = ''
           ErrID = 'L:4692/R:rdbetaD/F:risource.f'
           l_err = 1
           k_err = 0
           ierr  = 0
*-----------------------------------------------------------------------
*     file number of scratch file for CHECKING decay calculation
           iot      = 820
*     file number of RIsource.dat or ICRP-07.NDX based on DECDC2
*          iodecdc2 = 811
*     file number of RIsource.rad (radiation data) based on DECDC2
*          ioddc2r  = 812
*     file number of RIsource.bet (beta spectra) based on DECDC2
*          ioddc2b  = 813
*     file number of RIsource.ack (Auger-CK spectra) based on DECDC2
           ioddc2a  = 814
*     file number of RIsource.dat (decay data) formatted with DCHAIN
*          ioddc2s  = 815

c     err flag for this subroutine
           ierr     =   0
*-----------------------------------------------------------------------
*     address check of decay data file
*-----------------------------------------------------------------------

        call fchkd2(24,ntyp,fname,ierr)

          if( ierr .eq. 1 ) then
*           l_err = ill(jsn)
*           k_err = jsn
            go to 998
          end if

        icl = inumc(fname,1,100,' ') - 1
*       do n = 100, 1, -1
*         if( fname(n:n) .ne. ' ' ) go to 201
*       end do
*
* 201     icl = n
          lex = .false.

        if( ntyp .eq. 0 ) then

          fname = fname(1:icl-3)//'ack'
          inquire( file = fname, exist = lex )
          if( .not. lex ) then
            if( icl .gt. 12 ) then
              inquire( file = fname(1:icl-12)//'ICRP-07.ACK ',
     &                 exist = lex )
              if( lex ) then
                fname = fname(1:icl-12)//'ICRP-07.ACK '
              end if
            else if( icl .eq. 12 ) then
              inquire( file = 'ICRP-07.ACK', exist = lex )
              if( lex ) then
                fname = 'ICRP-07.ACK '
              end if
            end if
          end if

        else  ! ntyp = 1

          fname = fname(1:icl-3)//'ACK'

        end if

      open( ioddc2a, file = fname,
     &      form='formatted', status = 'unknown' )
        rewind ioddc2a

*-----------------------------------------------------------------------
*     Target nuclide
*-----------------------------------------------------------------------
*         chza = nucl
        call getelm(11,chaa,chza,ierr)

C       Sub. information: getelm(int.,char.,real,int.)
*         OUTPUT:
*          11: chza(Zaid+) -> chaa,7digit(H-1    )

          if( ierr .eq. 1 ) then
*           l_err = ill(jsn)
*           k_err = jsn
*             return
            go to 999
          end if

*-----------------------------------------------------------------------
*     Read the data on the energy and yield of Auger-CK electron
*-----------------------------------------------------------------------
  220 continue

        nnn = -1
        kkk =  0

c     else  ! ntyp = 1 (Extra file: ICRP-07.ACK)

        do n = 1, klmax
          read( ioddc2a,'(a200)', iostat = ios ) chin
*         write(iot,'(i7,'': '',a)') n, chin

          if( n .eq. 1 ) then
            if( ios .eq. -1 ) go to 903
          else
            if( ios .eq. -1 ) go to 280  ! read end
          end if

          if( nnn .eq. 0 ) then
            go to 280  ! successfully end
          else if( nnn .gt. 0 ) then
            nnn = nnn - 1
            kkk = kkk + 1
              read(chin(12:12+12-1),'(e15.7)') eng(kkk)
              read(chin(1:1+11-1),'(e15.7)')   rat(kkk)
              eng(kkk) = eng(kkk) / 1.0d+06
              if( kkk .eq. 1 ) then
                write(iot,'("  $ Auger-CK electron spectrum")')

                write(iot,'(2x,1pe15.7,1pe15.7,2x,
     &                    "$ Energy (MeV) and Absolute yield")')
     &                      eng(kkk), rat(kkk)

              else
                write(iot,'(2x,1pe15.7,1pe15.7)') eng(kkk), rat(kkk)

              end if

          end if

          if( chin(1:1+7-1) .eq. chaa ) then
            read(chin(22:22+11-1),'(i11)') nnn
            if( nnn .gt. rimax1 ) go to 965

          end if
        end do


      if( act .lt. 0.0d+00 ) then
      else
        go to 940
      end if

  280 continue

*-----------------------------------------------------------------------

        goto 999

*-----------------------------------------------------------------------
  901    m_err = 'RIsource.ack does NOT exist.'
         ErrCha = ''
         ErrID = 'L:4842/R:rdbetaD/F:risource.f'
           goto 998
  903    m_err = 'Data in RIsource.ack is not read.'
         ErrCha = ''
         ErrID = 'L:4846/R:rdbetaD/F:risource.f'
           goto 998
  925    m_err = 'RI source for electron need DECDC2(ICRP-07.ACK) file.'
         ErrCha = ''
         ErrID = 'L:4850/R:rdbetaD/F:risource.f'
         ierr = 1
           return
  940    m_err = 'Klmax in sub. rdbetaD is smaller than the data line.'
         ErrCha = ''
         ErrID = 'L:4855/R:rdbetaD/F:risource.f'
           goto 998
  965    m_err = 'Number of electron/positron was exceeded rimax1.'
         ErrCha = ''
         ErrID = 'L:4859/R:rdbetaD/F:risource.f'
           goto 998
*-----------------------------------------------------------------------
  998 continue
         ierr = 1

  999 continue
*-----------------------------------------------------------------------

      close(ioddc2a)
      return
      end

************************************************************************
*                                                                      *
      subroutine rdbetaD2(ntyp,chza,act,kkk,eng,rat,ierr)
*                                                                      *
*       read RIsource.rad file (energy and yield) for Auger electron   *
*         made by N.Matsuda on 2018/02/23                              *
*                                                                      *
************************************************************************

      implicit double precision (a-h,o-z)
      include 'risrcparam.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      common /error/  m_err, l_err, k_err
      character       m_err*200

      character chin*200, chlw*200, chcm*200

      integer  iot, iodecdc2, ioddc2r, ioddc2b, ioddc2a, ioddc2s
      integer ierr

      common /paran/  icfn(100), ilfn(100), chfn(100)
      character chfn*200
      character fname*200 ! T.Sato 2022/11/30, increase from 100 to 200 because inumc requests 200 character

      logical lex

      character chaa*7
      real chza
*     parameter ( rimax1    =  4,000 )
*     parameter ( rimax2    = 20,000 )
      double precision act, eng(rimax1), rat(rimax1)
      integer ntyp, icl, kkk

      integer klmax
      data klmax /500000/  ! 455,625 lines @ ICRP-07.RAD
      double precision smallv
      data smallv / 1.0d-36 /

*-----------------------------------------------------------------------
*     default values
*-----------------------------------------------------------------------
           m_err = ' Error in rdbetaD2 sub.(read RIsource.rad file)'
           ErrCha = ''
           ErrID = 'L:4918/R:rdbetaD2/F:risource.f'
           l_err = 1
           k_err = 0
           ierr  = 0
*-----------------------------------------------------------------------
*     file number of scratch file for CHECKING decay calculation
           iot      = 820
*     file number of RIsource.dat or ICRP-07.NDX based on DECDC2
*          iodecdc2 = 811
*     file number of RIsource.rad (radiation data) based on DECDC2
           ioddc2r  = 812
*     file number of RIsource.bet (beta spectra) based on DECDC2
*          ioddc2b  = 813
*     file number of RIsource.ack (Auger-CK spectra) based on DECDC2
*          ioddc2a  = 814
*     file number of RIsource.dat (decay data) formatted with DCHAIN
*          ioddc2s  = 815

c     err flag for this subroutine
           ierr     =   0
*-----------------------------------------------------------------------
*     address check of decay data file
*-----------------------------------------------------------------------

        call fchkd2(24,ntyp,fname,ierr)

          if( ierr .eq. 1 ) then
*           l_err = ill(jsn)
*           k_err = jsn
            go to 998
          end if

        icl = inumc(fname,1,100,' ') - 1
*       do n = 100, 1, -1
*         if( fname(n:n) .ne. ' ' ) go to 201
*       end do
*
* 201     icl = n
          lex = .false.

        if( ntyp .eq. 0 ) then

          fname = fname(1:icl-3)//'rad'
          inquire( file = fname, exist = lex )
          if( .not. lex ) then
            if( icl .gt. 12 ) then
              inquire( file = fname(1:icl-12)//'ICRP-07.RAD ',
     &                 exist = lex )
              if( lex ) then
                fname = fname(1:icl-12)//'ICRP-07.RAD '
              end if
            else if( icl .eq. 12 ) then
              inquire( file = 'ICRP-07.RAD', exist = lex )
              if( lex ) then
                fname = 'ICRP-07.RAD '
              end if
            end if
          end if

        else  ! ntyp = 1

          fname = fname(1:icl-3)//'RAD'

        end if


      open( ioddc2r, file = fname,
     &      form='formatted', status = 'unknown' )
        rewind ioddc2r

*-----------------------------------------------------------------------
*     Target nuclide
*-----------------------------------------------------------------------
*         chza = nucl
        call getelm(11,chaa,chza,ierr)

C       Sub. information: getelm(int.,char.,real,int.)
*         OUTPUT:
*          11: chza(Zaid+) -> chaa,7digit(H-1    )

          if( ierr .eq. 1 ) then
*           l_err = ill(jsn)
*           k_err = jsn
*             return
            go to 999
          end if

*-----------------------------------------------------------------------
*     Read the data on the energy and yield of Auger electron
*-----------------------------------------------------------------------
  220 continue

        nnn = -1
        kkk =  0

c     else  ! ntyp = 1 (Extra file: ICRP-07.RAD)

        do n = 1, klmax
          read( ioddc2r,'(a200)', iostat = ios ) chin

          if( n .eq. 1 ) then
            if( ios .eq. -1 ) go to 903
          else
            if( ios .eq. -1 ) go to 280  ! read end
          end if

          if( nnn .eq. 0 ) then
            go to 280  ! successfully end
          else if( nnn .gt. 0 ) then
            nnn = nnn - 1
            if( chin(27:27+3-1) .eq. ' AE' ) then
                kkk = kkk + 1
              if( kkk .gt. rimax1 ) go to 965
              read(chin(15:15+12-1),'(e15.7)') eng(kkk)
              read(chin(3:3+12-1),'(e15.7)')   rat(kkk)
              if( kkk .eq. 1 ) then
                write(iot,'(2x,1pe15.7,1pe15.7,2x,
     &                    "$ Energy (MeV) and Absolute yield")')
     &                      eng(kkk), rat(kkk)

              else
                write(iot,'(2x,1pe15.7,1pe15.7)') eng(kkk), rat(kkk)

              end if

            end if

          end if

          if( chin(1:1+7-1) .eq. chaa ) then
            read(chin(21:21+9-1),'(i9)') nnn

          end if
        end do


      if( act .lt. 0.0d+00 ) then
      else
        go to 940
      end if

  280 continue

*-----------------------------------------------------------------------

        goto 999

*-----------------------------------------------------------------------
  901    m_err = 'RIsource.rad does NOT exist.'
         ErrCha = ''
         ErrID = 'L:5068/R:rdbetaD2/F:risource.f'
           goto 998
  903    m_err = 'Data in RIsource.rad is not read.'
         ErrCha = ''
         ErrID = 'L:5072/R:rdbetaD2/F:risource.f'
           goto 998
  925    m_err = 'RI source for electron need DECDC2(ICRP-07.RAD) file.'
         ErrCha = ''
         ErrID = 'L:5076/R:rdbetaD2/F:risource.f'
         ierr = 1
           return
  940    m_err = 'Klmax in sub. rdbetaD2 is smaller than the data line.'
         ErrCha = ''
         ErrID = 'L:5081/R:rdbetaD2/F:risource.f'
           goto 998
  965    m_err = 'Number of electron was exceeded rimax1.'
         ErrCha = ''
         ErrID = 'L:5085/R:rdbetaD2/F:risource.f'
           goto 998
*-----------------------------------------------------------------------
  998 continue
         ierr = 1

  999 continue
*-----------------------------------------------------------------------

      close(ioddc2r)
      return
      end

************************************************************************
*                                                                      *
      subroutine rdbetaD3(ntyp,chza,act,kkk,eng,rat,ierr)
*                                                                      *
*       read RIsource.rad file (energy and yield)                      *
*                                    for internal conversion electron  *
*         made by N.Matsuda on 2018/08/15                              *
*                                                                      *
************************************************************************

      implicit double precision (a-h,o-z)
      include 'risrcparam.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      common /error/  m_err, l_err, k_err
      character       m_err*200

      character chin*200, chlw*200, chcm*200

      integer  iot, iodecdc2, ioddc2r, ioddc2b, ioddc2a, ioddc2s
      integer ierr

      common /paran/  icfn(100), ilfn(100), chfn(100)
      character chfn*200
      character fname*200 ! T.Sato 2022/11/30, increase from 100 to 200 because inumc requests 200 character

      logical lex

      character chaa*7
      real chza
*     parameter ( rimax1    =  4,000 )
*     parameter ( rimax2    = 20,000 )
      double precision act, eng(rimax1), rat(rimax1)
      integer ntyp, icl, kkk

      integer klmax
      data klmax /500000/  ! 455,625 lines @ ICRP-07.RAD
      double precision smallv
      data smallv / 1.0d-36 /

*-----------------------------------------------------------------------
*     default values
*-----------------------------------------------------------------------
           m_err = ' Error in rdbetaD3 sub.(read RIsource.rad file)'
           ErrCha = ''
           ErrID = 'L:5145/R:rdbetaD3/F:risource.f'
           l_err = 1
           k_err = 0
           ierr  = 0
*-----------------------------------------------------------------------
*     file number of scratch file for CHECKING decay calculation
           iot      = 820
*     file number of RIsource.dat or ICRP-07.NDX based on DECDC2
*          iodecdc2 = 811
*     file number of RIsource.rad (radiation data) based on DECDC2
           ioddc2r  = 812
*     file number of RIsource.bet (beta spectra) based on DECDC2
*          ioddc2b  = 813
*     file number of RIsource.ack (Auger-CK spectra) based on DECDC2
*          ioddc2a  = 814
*     file number of RIsource.dat (decay data) formatted with DCHAIN
*          ioddc2s  = 815

c     err flag for this subroutine
           ierr     =   0
*-----------------------------------------------------------------------
*     address check of decay data file
*-----------------------------------------------------------------------

        call fchkd2(24,ntyp,fname,ierr)

          if( ierr .eq. 1 ) then
*           l_err = ill(jsn)
*           k_err = jsn
            go to 998
          end if

        icl = inumc(fname,1,100,' ') - 1
*       do n = 100, 1, -1
*         if( fname(n:n) .ne. ' ' ) go to 201
*       end do
*
* 201     icl = n
          lex = .false.

        if( ntyp .eq. 0 ) then

          fname = fname(1:icl-3)//'rad'
          inquire( file = fname, exist = lex )
          if( .not. lex ) then
            if( icl .gt. 12 ) then
              inquire( file = fname(1:icl-12)//'ICRP-07.RAD ',
     &                 exist = lex )
              if( lex ) then
                fname = fname(1:icl-12)//'ICRP-07.RAD '
              end if
            else if( icl .eq. 12 ) then
              inquire( file = 'ICRP-07.RAD', exist = lex )
              if( lex ) then
                fname = 'ICRP-07.RAD '
              end if
            end if
          end if

        else  ! ntyp = 1

          fname = fname(1:icl-3)//'RAD'

        end if


      open( ioddc2r, file = fname,
     &      form='formatted', status = 'unknown' )
        rewind ioddc2r

*-----------------------------------------------------------------------
*     Target nuclide
*-----------------------------------------------------------------------
*         chza = nucl
        call getelm(11,chaa,chza,ierr)

C       Sub. information: getelm(int.,char.,real,int.)
*         OUTPUT:
*          11: chza(Zaid+) -> chaa,7digit(H-1    )

          if( ierr .eq. 1 ) then
*           l_err = ill(jsn)
*           k_err = jsn
*             return
            go to 999
          end if

*-----------------------------------------------------------------------
*     Read the data on the energy and yield of internal conversion electron
*-----------------------------------------------------------------------
  220 continue

        nnn = -1
        kkk =  0

c     else  ! ntyp = 1 (Extra file: ICRP-07.RAD)

        do n = 1, klmax
          read( ioddc2r,'(a200)', iostat = ios ) chin

          if( n .eq. 1 ) then
            if( ios .eq. -1 ) go to 903
          else
            if( ios .eq. -1 ) go to 280  ! read end
          end if

          if( nnn .eq. 0 ) then
            go to 280  ! successfully end
          else if( nnn .gt. 0 ) then
            nnn = nnn - 1
            if( chin(27:27+3-1) .eq. ' IE' ) then
                kkk = kkk + 1
              if( kkk .gt. rimax1 ) go to 965
              read(chin(15:15+12-1),'(e15.7)') eng(kkk)
              read(chin(3:3+12-1),'(e15.7)')   rat(kkk)
              if( kkk .eq. 1 ) then
                write(iot,'(2x,1pe15.7,1pe15.7,2x,
     &                    "$ Energy (MeV) and Absolute yield")')
     &                      eng(kkk), rat(kkk)

              else
                write(iot,'(2x,1pe15.7,1pe15.7)') eng(kkk), rat(kkk)

              end if

            end if

          end if

          if( chin(1:1+7-1) .eq. chaa ) then
            read(chin(21:21+9-1),'(i9)') nnn

          end if
        end do


      if( act .lt. 0.0d+00 ) then
      else
        go to 940
      end if

  280 continue

*-----------------------------------------------------------------------

        goto 999

*-----------------------------------------------------------------------
  901    m_err = 'RIsource.rad does NOT exist.'
         ErrCha = ''
         ErrID = 'L:5295/R:rdbetaD3/F:risource.f'
           goto 998
  903    m_err = 'Data in RIsource.rad is not read.'
         ErrCha = ''
         ErrID = 'L:5299/R:rdbetaD3/F:risource.f'
           goto 998
  925    m_err = 'RI source for electron need DECDC2(ICRP-07.RAD) file.'
         ErrCha = ''
         ErrID = 'L:5303/R:rdbetaD3/F:risource.f'
         ierr = 1
           return
  940    m_err = 'Klmax in sub. rdbetaD3 is smaller than the data line.'
         ErrCha = ''
         ErrID = 'L:5308/R:rdbetaD3/F:risource.f'
           goto 998
  965    m_err = 'Number of electron was exceeded rimax1.'
         ErrCha = ''
         ErrID = 'L:5312/R:rdbetaD3/F:risource.f'
           goto 998
*-----------------------------------------------------------------------
  998 continue
         ierr = 1

  999 continue
*-----------------------------------------------------------------------

      close(ioddc2r)
      return
      end

************************************************************************
*                                                                      *
      subroutine rdalpha(ntyp,chza,act,kkk,eng,rat,ierr)
*                                                                      *
*       read RIsource.rad file (energy and yield)                      *
*         made by N.Matsuda on 2017/05/29                              *
*                                                                      *
************************************************************************

      implicit double precision (a-h,o-z)
      include 'risrcparam.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      common /error/  m_err, l_err, k_err
      character       m_err*200

      character chin*200, chlw*200, chcm*200

      integer  iot, iodecdc2, ioddc2r, ioddc2b, ioddc2a, ioddc2s
      integer ierr

      common /paran/  icfn(100), ilfn(100), chfn(100)
      character chfn*200
      character fname*200 ! T.Sato 2022/11/30, increase from 100 to 200 because inumc requests 200 character

      logical lex

      character chaa*7
      real chza
*     parameter ( rimax1    =  4,000 )
*     parameter ( rimax2    = 20,000 )
      double precision act, eng(rimax1), rat(rimax1)
      integer ntyp, icl, kkk

      integer klmax
      data klmax /500000/  ! 455,625 lines @ ICRP-07.RAD
      double precision smallv
      data smallv / 1.0d-36 /

*-----------------------------------------------------------------------
*     default values
*-----------------------------------------------------------------------
           m_err = ' Error in rdalpha sub.(read RIsource.rad file)'
           ErrCha = ''
           ErrID = 'L:5371/R:rdalpha/F:risource.f'
           l_err = 1
           k_err = 0
           ierr  = 0
*-----------------------------------------------------------------------
*     file number of scratch file for CHECKING decay calculation
           iot      = 820
*     file number of RIsource.dat or ICRP-07.NDX based on DECDC2
*          iodecdc2 = 811
*     file number of RIsource.rad (radiation data) based on DECDC2
           ioddc2r  = 812
*     file number of RIsource.bet (beta spectra) based on DECDC2
*          ioddc2b  = 813
*     file number of RIsource.ack (Auger-CK spectra) based on DECDC2
*          ioddc2a  = 814
*     file number of RIsource.dat (decay data) formatted with DCHAIN
*          ioddc2s  = 815

c     err flag for this subroutine
           ierr     =   0
*-----------------------------------------------------------------------
*     address check of decay data file
*-----------------------------------------------------------------------

        call fchkd2(24,ntyp,fname,ierr)

          if( ierr .eq. 1 ) then
*           l_err = ill(jsn)
*           k_err = jsn
            go to 998
          end if

        icl = inumc(fname,1,100,' ') - 1
*       do n = 100, 1, -1
*         if( fname(n:n) .ne. ' ' ) go to 201
*       end do
*
* 201     icl = n
          lex = .false.

        if( ntyp .eq. 0 ) then

          fname = fname(1:icl-3)//'rad'
          inquire( file = fname, exist = lex )
          if( .not. lex ) then
            if( icl .gt. 12 ) then
              inquire( file = fname(1:icl-12)//'ICRP-07.RAD ',
     &                 exist = lex )
              if( lex ) then
                fname = fname(1:icl-12)//'ICRP-07.RAD '
              end if
            else if( icl .eq. 12 ) then
              inquire( file = 'ICRP-07.RAD', exist = lex )
              if( lex ) then
                fname = 'ICRP-07.RAD '
              end if
            end if
          end if

        else  ! ntyp = 1

          fname = fname(1:icl-3)//'RAD'

        end if


      open( ioddc2r, file = fname,
     &      form='formatted', status = 'unknown' )
        rewind ioddc2r

*-----------------------------------------------------------------------
*     Target nuclide
*-----------------------------------------------------------------------
*         chza = nucl
        call getelm(11,chaa,chza,ierr)

C       Sub. information: getelm(int.,char.,real,int.)
*         OUTPUT:
*          11: chza(Zaid+) -> chaa,7digit(H-1    )

          if( ierr .eq. 1 ) then
*           l_err = ill(jsn)
*           k_err = jsn
*             return
            go to 999
          end if

*-----------------------------------------------------------------------
*     Read the data on the energy and yield of alpha particle
*-----------------------------------------------------------------------
  220 continue

        nnn = -1
        kkk =  0

c     else  ! ntyp = 1 (Extra file: ICRP-07.RAD)

        do n = 1, klmax
          read( ioddc2r,'(a200)', iostat = ios ) chin

          if( n .eq. 1 ) then
            if( ios .eq. -1 ) go to 903
          else
            if( ios .eq. -1 ) go to 280  ! read end
          end if

          if( nnn .eq. 0 ) then
            go to 280  ! successfully end
          else if( nnn .gt. 0 ) then
            nnn = nnn - 1
            if( chin(27:27+3-1) .eq. '  A' ) then
                kkk = kkk + 1
              if( kkk .gt. rimax1 ) go to 965
              read(chin(15:15+12-1),'(e15.7)') eng(kkk)
              read(chin(3:3+12-1),'(e15.7)')   rat(kkk)
              if( kkk .eq. 1 ) then
                write(iot,'(2x,1pe15.7,1pe15.7,2x,
     &                    "$ Energy (MeV) and Absolute yield")')
     &                      eng(kkk), rat(kkk)

              else
                write(iot,'(2x,1pe15.7,1pe15.7)') eng(kkk), rat(kkk)

              end if

            end if

          end if

          if( chin(1:1+7-1) .eq. chaa ) then
            read(chin(21:21+9-1),'(i9)') nnn

          end if
        end do


      if( act .lt. 0.0d+00 ) then
      else
        go to 940
      end if

  280 continue

*-----------------------------------------------------------------------

        goto 999

*-----------------------------------------------------------------------
  901    m_err = 'RIsource.rad does NOT exist.'
         ErrCha = ''
         ErrID = 'L:5521/R:rdalpha/F:risource.f'
           goto 998
  903    m_err = 'Data in RIsource.rad is not read.'
         ErrCha = ''
         ErrID = 'L:5525/R:rdalpha/F:risource.f'
           goto 998
  925    m_err = 'RI source for alpha need DECDC2(ICRP-07.RAD) file.'
         ErrCha = ''
         ErrID = 'L:5529/R:rdalpha/F:risource.f'
         ierr = 1
           return
  940    m_err = 'Klmax in sub. rdalpha is smaller than the data line.'
         ErrCha = ''
         ErrID = 'L:5534/R:rdalpha/F:risource.f'
           goto 998
  965    m_err = 'Number of alpha particle was exceeded rimax1.'
         ErrCha = ''
         ErrID = 'L:5538/R:rdalpha/F:risource.f'
           goto 998
*-----------------------------------------------------------------------
  998 continue
         ierr = 1

  999 continue
*-----------------------------------------------------------------------

      close(ioddc2r)
      return
      end

************************************************************************
*                                                                      *
      subroutine rdsfneutron(ntyp,chza,act,kkk,eng,rat,ierr)
*                                                                      *
*       read RIsource.rad file (energy and yield)                      *
*         added by N.Furutachi on 2022/12/26                           *
*                                                                      *
************************************************************************

      implicit double precision (a-h,o-z)
      include 'risrcparam.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      common /error/  m_err, l_err, k_err
      character       m_err*200

      character chin*200, chlw*200, chcm*200

      integer  iot, iodecdc2, ioddc2r, ioddc2b, ioddc2a, ioddc2s
      integer ierr

      common /paran/  icfn(100), ilfn(100), chfn(100)
      character chfn*200
      character fname*200 ! T.Sato 2022/11/30, increase from 100 to 200 because inumc requests 200 character

      logical lex

      character chaa*7
      real chza
*     parameter ( rimax1    =  4,000 )
*     parameter ( rimax2    = 20,000 )
      double precision act, eng(rimax1), rat(rimax1)
      integer ntyp, icl, kkk

      integer klmax
      data klmax /500000/  ! 455,625 lines @ ICRP-07.RAD
      double precision smallv
      data smallv / 1.0d-36 /

*-----------------------------------------------------------------------
*     default values
*-----------------------------------------------------------------------
           m_err = ' Error in rdalpha sub.(read RIsource.rad file)'
           ErrCha = ''
           ErrID = 'L:5597/R:rdsfneutron/F:risource.f'
           l_err = 1
           k_err = 0
           ierr  = 0
*-----------------------------------------------------------------------
*     file number of scratch file for CHECKING decay calculation
           iot      = 820
*     file number of RIsource.dat or ICRP-07.NDX based on DECDC2
*          iodecdc2 = 811
*     file number of RIsource.rad (radiation data) based on DECDC2
*          ioddc2r  = 812
*     file number of RIsource.bet (beta spectra) based on DECDC2
*          ioddc2b  = 813
*     file number of RIsource.ack (Auger-CK spectra) based on DECDC2
*          ioddc2a  = 814
*     file number of RIsource.dat (decay data) formatted with DCHAIN
*          ioddc2s  = 815
*     file number of RIsource.nsf (SF neutron spectra) based on DECDC2
           ioddc2n  = 816

c     err flag for this subroutine
           ierr     =   0
*-----------------------------------------------------------------------
*     address check of decay data file
*-----------------------------------------------------------------------

        call fchkd2(24,ntyp,fname,ierr)

          if( ierr .eq. 1 ) then
            go to 998
          end if

        icl = inumc(fname,1,100,' ') - 1

        lex = .false.

        if( ntyp .eq. 0 ) then

          fname = fname(1:icl-3)//'nsf'
          inquire( file = fname, exist = lex )
          if( .not. lex ) then
            if( icl .gt. 12 ) then
              inquire( file = fname(1:icl-12)//'ICRP-07.NSF ',
     &                 exist = lex )
              if( lex ) then
                fname = fname(1:icl-12)//'ICRP-07.NSF '
              end if
            else if( icl .eq. 12 ) then
              inquire( file = 'ICRP-07.NSF', exist = lex )
              if( lex ) then
                fname = 'ICRP-07.NSF '
              end if
            end if
          end if

        else  ! ntyp = 1

          fname = fname(1:icl-3)//'NSF'

        end if

      open( ioddc2n, file = fname,
     &      form='formatted', status = 'unknown' )
        rewind ioddc2n

*-----------------------------------------------------------------------
*     Target nuclide
*-----------------------------------------------------------------------
        call getelm(11,chaa,chza,ierr)

*         OUTPUT:
*          11: chza(Zaid+) -> chaa,7digit(H-1    )

          if( ierr .eq. 1 ) then
            go to 999
          end if

*-----------------------------------------------------------------------
*     Read the data on the energy and the ratio of SF neutron
*-----------------------------------------------------------------------
  220 continue

        nnn = -1
        kkk =  0

        do n = 1, klmax
          read( ioddc2n,'(a200)', iostat = ios ) chin

          if( n .eq. 1 ) then
            if( ios .eq. -1 ) go to 903
          else
            if( ios .eq. -1 ) go to 280  ! read end
          end if

          if( nnn .eq. 0 ) then
            go to 280  ! successfully end
          else if( nnn .gt. 0 ) then
            nnn = nnn - 1
            kkk = kkk + 1
            if( kkk .gt. rimax1 ) go to 965
            read(chin(1:1+8-1),'(e8.2)') eng(kkk)
            read(chin(19:19+11-1),'(e11.5)') rat(kkk)
            if( nnn.eq.0 ) then
              kkk = kkk + 1
              read(chin(10:10+8-1),'(e8.2)') eng(kkk)
            end if
            if( kkk .eq. 1 ) then
              write(iot,'(2x,1pe15.7,1pe15.7,2x,
     &                    "$ Energy (MeV) and ratio")')
     &                     eng(kkk), rat(kkk)

            else
              write(iot,'(2x,1pe15.7,1pe15.7)') eng(kkk), rat(kkk)

            end if

          end if

          if( chin(1:1+7-1) .eq. chaa ) then
            read(chin(27:27+3-1),'(i3)') nnn

          end if
        end do


      if( act .lt. 0.0d+00 ) then
      else
        go to 940
      end if

  280 continue

*-----------------------------------------------------------------------

        goto 999

*-----------------------------------------------------------------------
  901    m_err = 'RIsource.rad does NOT exist.'
         ErrCha = ''
         ErrID = 'L:5736/R:rdsfneutron/F:risource.f'
           goto 998
  903    m_err = 'Data in RIsource.rad is not read.'
         ErrCha = ''
         ErrID = 'L:5740/R:rdsfneutron/F:risource.f'
           goto 998
  925    m_err = 'RI source for alpha need DECDC2(ICRP-07.RAD) file.'
         ErrCha = ''
         ErrID = 'L:5744/R:rdsfneutron/F:risource.f'
         ierr = 1
           return
  940    m_err = 'Klmax in sub. rdsfneutron is smaller than data line.'
         ErrCha = ''
         ErrID = 'L:5749/R:rdsfneutron/F:risource.f'
           goto 998
  965    m_err = 'Number of SF neutron spectra was exceeded rimax1.'
         ErrCha = ''
         ErrID = 'L:5753/R:rdsfneutron/F:risource.f'
           goto 998
*-----------------------------------------------------------------------
  998 continue
         ierr = 1

  999 continue
*-----------------------------------------------------------------------

      close(ioddc2n)
      return
      end

************************************************************************
*                                                                      *
      subroutine fchkd2(id,ndtype,fname,ierr)
*                                                                      *
*       file or folder check for RI source                             *
*         made by N.Matsuda on 2016/07/31                              *
*         revised by N.Matsuda on 2017/05/29                           *
*                                                                      *
************************************************************************

      implicit double precision (a-h,o-z)

      include 'err.inc'

*-----------------------------------------------------------------------

      common /error/  m_err, l_err, k_err
      character       m_err*200

      integer iot, ierr

      common /paran/  icfn(100), ilfn(100), chfn(100)
      character chfn*200
      character fname*200 ! T.Sato 2022/11/30, increase from 100 to 200 because inumc requests 200 character

      logical lex
      integer id, ndtype

*-----------------------------------------------------------------------
*     file number of scratch file for CHECKING decay calculation
           iot      = 820

c     err flag for this subroutine
           ierr     =   0
*-----------------------------------------------------------------------
*     address check of decay data file (Decay dara lib. or ICRP-07.NDX)
*-----------------------------------------------------------------------

*       icl = inumc(chfn(id),1,100,' ') - 1
        do n = 100, 1, -1
          if( chfn(id)(n:n) .ne. ' ' ) go to 201
        end do

  201     icl = n
          lex = .false.

C     DATA TYPE: Decay data lib. (RIsource.dat)  ! File
      if( chfn(id)(icl-11:icl-4) .eq. 'RIsource' ) then
          ndtype = 0
          fname = chfn(id)(1:icl-4)//'.dat'
            icl = icl

      else if( chfn(id)(icl-8:icl-1) .eq. 'RIsource' ) then
          ndtype = 0
          fname = chfn(id)(1:icl-1)//'.dat'
            icl = icl + 3

      else if( chfn(id)(icl-7:icl) .eq. 'RIsource' ) then
          ndtype = 0
          fname = chfn(id)(1:icl)//'.dat'
            icl = icl + 4

C     DATA TYPE: DECDC2 (original data)          ! Extra file
      else if( chfn(id)(icl-10:icl-4) .eq. 'ICRP-07' ) then
          ndtype = 1
          fname = chfn(id)(1:icl-4)//'.NDX'
            icl = icl

      else if( chfn(id)(icl-7:icl-1) .eq. 'ICRP-07' ) then
          ndtype = 1
          fname = chfn(id)(1:icl-1)//'.NDX'
            icl = icl + 3

      else if( chfn(id)(icl-6:icl) .eq. 'ICRP-07' ) then
          ndtype = 1
          fname = chfn(id)(1:icl)//'.NDX'
            icl = icl + 4

C     DATA TYPE: Decay data lib. (RIsource.dat)  ! Folder
      else if( chfn(id)(icl:icl) .eq. '/' .or.
     &         chfn(id)(icl:icl) .eq. '\\' ) then
          ndtype = 0
          fname = chfn(id)(1:icl-1)//'/RIsource.dat'
            icl = icl + 13

      else
          ndtype = 0
          fname = chfn(id)(1:icl)//'/RIsource.dat'
            icl = icl + 13

      end if

*-----------------------------------------------------------------------
      inquire( file = fname, exist = lex )
      if( .not. lex ) then
        if( ndtype .eq. 0 ) then
*         write(iot,'(2x,''RIsource.dat does NOT exist.'')')
            go to 901
        else
*         write(iot,'(2x,''DECDC2 index file does NOT exist.'')')
            go to 911
        end if

      else
*       if( ndtype .eq. 0 ) then
*         O.K. ... RIsource.dat file exists.
*       else
*         O.K. ... DECDC2 index file exists.
*       end if
*
      end if

          go to 500

*-----------------------------------------------------------------------
*     the other files
*-----------------------------------------------------------------------

      if( ndtype .eq. 0 ) then
        inquire( file = fname(1:icl-3)//'rad', exist = lex )
        if( .not. lex ) go to 902

        inquire( file = fname(1:icl-3)//'bet', exist = lex )
        if( .not. lex ) go to 903

        inquire( file = fname(1:icl-3)//'ack', exist = lex )
        if( .not. lex ) go to 904
        inquire( file = fname(1:icl-3)//'nsf', exist = lex )
        if( .not. lex ) go to 905

      else
        inquire( file = fname(1:icl-3)//'RAD', exist = lex )
        if( .not. lex ) go to 912

        inquire( file = fname(1:icl-3)//'BET', exist = lex )
        if( .not. lex ) go to 913

        inquire( file = fname(1:icl-3)//'ACK', exist = lex )
        if( .not. lex ) go to 914
        inquire( file = fname(1:icl-3)//'NSF', exist = lex )
        if( .not. lex ) go to 915

      end if

*-----------------------------------------------------------------------
  500 continue

      if( icl .gt. 100 ) then
         go to 951
      else if( icl .eq. 100 ) then
      else
        do n = icl+1, 100
           fname(n:n) = ' '
        end do

      end if

*-----------------------------------------------------------------------

         goto 999

*-----------------------------------------------------------------------
  901    m_err = 'Decay data file does NOT exist: file(24) = '//chfn(id)
         ErrCha = ''
         ErrID = 'L:5930/R:fchkd2/F:risource.f'
           goto 998
  902    m_err = 'RIsource.rad does NOT exist.'
         ErrCha = ''
         ErrID = 'L:5934/R:fchkd2/F:risource.f'
           goto 998
  903    m_err = 'RIsource.bet does NOT exist.'
         ErrCha = ''
         ErrID = 'L:5938/R:fchkd2/F:risource.f'
           goto 998
  904    m_err = 'RIsource.ack does NOT exist.'
         ErrCha = ''
         ErrID = 'L:5942/R:fchkd2/F:risource.f'
           goto 998
  905    m_err = 'RIsource.nsf does NOT exist.'
         ErrCha = ''
         ErrID = 'L:5946/R:fchkd2/F:risource.f'
           goto 998

  911    m_err = 'ICRP-07.NDX does NOT exist: file(24) = '//chfn(id)
         ErrCha = ''
         ErrID = 'L:5951/R:fchkd2/F:risource.f'
           goto 998
  912    m_err = 'ICRP-07.RAD does NOT exist.'
         ErrCha = ''
         ErrID = 'L:5955/R:fchkd2/F:risource.f'
           goto 998
  913    m_err = 'ICRP-07.BET does NOT exist.'
         ErrCha = ''
         ErrID = 'L:5959/R:fchkd2/F:risource.f'
           goto 998
  914    m_err = 'ICRP-07.ACK does NOT exist.'
         ErrCha = ''
         ErrID = 'L:5963/R:fchkd2/F:risource.f'
           goto 998
  915    m_err = 'ICRP-07.NSF does NOT exist.'
         ErrCha = ''
         ErrID = 'L:5967/R:fchkd2/F:risource.f'
           goto 998

  951    m_err = 'Address of file(24) was exceeded 100 characters.'
         ErrCha = ''
         ErrID = 'L:5972/R:fchkd2/F:risource.f'
           goto 998
*-----------------------------------------------------------------------
  998 continue
         ierr = 1

  999 continue
*-----------------------------------------------------------------------
      return
      end

************************************************************************
*                                                                      *
      subroutine getelm(m,nucl,chal,ierr)
*                                                                      *
*     get and translation ( H-1, H1, 1H, 1001 <=> zaid )               *
*     made by N.Matsuda on 2016/07/31                                  *
*                                                                      *
*        m : 1 or 2 get for read input, 11-13 translation for write    *
*       nucl <- chlw (lower-case character), 7 characters              *
*       chal <=> zzzaaa.m (actual number), iz, ia and im               *
*                                                                      *
************************************************************************

      implicit double precision (a-h,o-z)

      include 'err.inc'

*-----------------------------------------------------------------------

      common /error/  m_err, l_err, k_err
      character       m_err*200

      integer m, iz, ia, im
      real chal
      character nucl*7, zcha*2, acha3*3, acha2*2, acha1*1, mcha*1
*     dimension

*-----------------------------------------------------------------------
      character element(104)*2,elmnt(104)*2

      data element /
     &    'h ','he','li','be','b ','c ','n ','o ','f ','ne',
     &    'na','mg','al','si','p ','s ','cl','ar','k ','ca',
     &    'sc','ti','v ','cr','mn','fe','co','ni','cu','zn',
     &    'ga','ge','as','se','br','kr','rb','sr','y ','zr',
     &    'nb','mo','tc','ru','rh','pd','ag','cd','in','sn',
     &    'sb','te','i ','xe','cs','ba','la','ce','pr','nd',
     &    'pm','sm','eu','gd','tb','dy','ho','er','tm','yb',
     &    'lu','hf','ta','w ','re','os','ir','pt','au','hg',
     &    'tl','pb','bi','po','at','rn','fr','ra','ac','th',
     &    'pa','u ','np','pu','am','cm','bk','cf','es','fm',
     &    'md','no','lr','ku'/

      data elmnt /
     &    'H ','He','Li','Be','B ','C ','N ','O ','F ','Ne',
     &    'Na','Mg','Al','Si','P ','S ','Cl','Ar','K ','Ca',
     &    'Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn',
     &    'Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y ','Zr',
     &    'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn',
     &    'Sb','Te','I ','Xe','Cs','Ba','La','Ce','Pr','Nd',
     &    'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb',
     &    'Lu','Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg',
     &    'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th',
     &    'Pa','U ','Np','Pu','Am','Cm','Bk','Cf','Es','Fm',
     &    'Md','No','Lr','Ku' /

*-----------------------------------------------------------------------
*     get and check
*-----------------------------------------------------------------------
  100 continue
      if( m .eq. 1 .or. m .eq. 2 ) then
        do n = 7, 3, -1
          if( nucl(n:n) .ne. ' ' ) then
            icl = n
            go to 101
          else
          end if
        end do
            icl = 2
 101    continue

        if( nucl(1:1) .ge. 'A' .and. nucl(1:1) .le. 'z' ) then
          if( nucl(2:2) .eq. '-' ) then
            zcha = nucl(1:1)//' '
            j1 = 1
          else if( nucl(2:2) .ge. '0' .and. nucl(2:2) .le. '9' ) then
            zcha = nucl(1:1)//' '
            j1 = 1
          else if( nucl(3:3) .ge. 'A' .and. nucl(3:3) .le. 'z' ) then
            go to 955
          else
            zcha = nucl(1:2)
            j1 = 2
          end if

*       search the element
          do iz = 1, 104
            if( zcha(1:1) .ge. 'A' .and. zcha(1:1) .le. 'Z' ) then
              if ( zcha(1:2) .eq. elmnt(iz)(1:2) ) goto 150
            else
              if ( zcha(1:2) .eq. element(iz)(1:2) ) goto 150
            end if
          end do
            go to 951      ! error -not found in the list
  150     continue

          j1 = j1 + 1
          if( j1 .gt. icl ) then
                   ia = 0
                     im = 0
            if( m .eq. 2 ) go to 952    ! natural abundance
          else
            if( nucl(j1:j1) .eq. '-' ) j1 = j1 + 1
            if( nucl(icl:icl) .eq. 'm' .or. nucl(icl:icl) .eq. 'n'
     &                               .or. nucl(icl:icl) .eq. 'g' ) then
*             call onum(nucl,j1,icl-1,cvvv,ierr)
              read(nucl(j1:icl-1),*) cvvv
                   ia = nint( cvvv )
                   if( nucl(icl:icl) .eq. 'm' ) then
                     im = 1
                   else if( nucl(icl:icl) .eq. 'n' ) then
                     im = 2
                   else if( nucl(icl:icl) .eq. 'g' ) then
                     im = 0
                   else
                     im = 0
                   end if
              if( m .eq. 2 .and. ia .eq. 0 ) go to 952
              if( ia .gt. 276 ) go to 953      ! Rf-276
            else
*             call onum(nucl,j1,icl,cvvv,ierr)
              read(nucl(j1:icl),*) cvvv
                   ia = nint( cvvv )
                     im = 0
              if( m .eq. 2 .and. ia .eq. 0 ) go to 952
              if( ia .gt. 276 ) go to 953      ! Rf-276
            end if
          end if

               chal = 1000 * iz + ia + 0.1 * im

        else
          if( nucl(icl:icl) .ge. 'A' .and. nucl(icl:icl) .le. 'z' ) then
            if( nucl(icl-1:icl-1) .ge. 'A'
     &                           .and. nucl(icl-1:icl-1) .le. 'z' ) then
              zcha = nucl(icl-1:icl)
*             call onum(nucl,1,icl-2,cvvv,ierr)
              read(nucl(1:icl-2),*) cvvv
                 ia = nint( cvvv )
                   im = 0
            else
              zcha = nucl(icl:icl)
*             call onum(nucl,1,icl-1,cvvv,ierr)
              read(nucl(1:icl-1),*) cvvv
                 ia = nint( cvvv )
                   im = 0
            end if

*       search the element
            do iz = 1, 104
              if( zcha(1:1) .ge. 'A' .and. zcha(1:1) .le. 'Z' ) then
                if ( zcha(1:2) .eq. elmnt(iz)(1:2) ) goto 170
              else
                if ( zcha(1:2) .eq. element(iz)(1:2) ) goto 170
              end if
            end do
              go to 951      ! error -not found in the list
  170       continue

            if( ia .gt. 276 ) go to 953      ! Rf-276

               chal = 1000 * iz + ia + 0.1 * im

          else
*             call onum(nucl,1,icl,cvvv,ierr)
              read(nucl(1:icl),*) cvvv
               chal = cvvv
               iz = int( chal / 1000 )
               if( iz .gt. 104 ) then
                 zcha = ''
                 go to 951
               else
                 zcha = elmnt(iz)
               end if

               ia = int( chal - iz * 1000 )
               if( m .eq. 2 .and. ia .eq. 0 ) go to 952
               if( ia .gt. 276 ) go to 953      ! Rf-276

               im = nint( chal * 10 - iz * 10000 -ia * 10 )
               if( im .eq. 0 ) then
                 mcha = ' '
               else if( im .eq. 1 ) then
                 mcha = 'm'
               else if( im . eq. 2 ) then
                 mcha = 'n'
               else
                 write( mcha, '(i1.1)' ) im
                 go to 954
               end if

          end if
        end if


      end if
*-----------------------------------------------------------------------
*     translation
*-----------------------------------------------------------------------
  600 continue
      if( m .eq. 11 .or. m .eq. 12 .or. m .eq. 13 .or. m .eq. 14 ) then
        nucl = '       '
        iz = int( chal / 1000 )
        if( iz .gt. 104 ) then
          zcha = ''
          go to 951
        else
          zcha = elmnt(iz)
        end if

        ia = int( chal - iz * 1000 )
          if( ia .gt. 276 ) go to 953      ! Rf-276

        im = nint( chal * 10 - iz * 10000 -ia * 10 )
        if( im .eq. 0 ) then
          mcha = ' '
        else if( im .eq. 1 ) then
          mcha = 'm'
        else if( im . eq. 2 ) then
          mcha = 'n'
        else
          write( mcha, '(i1.1)' ) im
          go to 954
        end if

        if( m .eq. 11 ) then
          if( zcha(2:2) .eq. ' ' ) then
            if( ia .gt. 99 ) then
              write( acha3, '(i3.3)' ) ia
                nucl = zcha(1:1)//'-'//acha3//mcha//'  '
            else if( ia .gt. 9 .and. ia .lt. 100 ) then
              write( acha2, '(i2.2)' ) ia
                nucl = zcha(1:1)//'-'//acha2//mcha//'   '
            else if( ia .lt. 10 ) then
              write( acha1, '(i1.1)' ) ia
                nucl = zcha(1:1)//'-'//acha1//mcha//'    '
            end if
          else
            if( ia .gt. 99 ) then
              write( acha3, '(i3.3)' ) ia
                nucl = zcha(1:2)//'-'//acha3//mcha//' '
            else if( ia .gt. 9 .and. ia .lt. 100 ) then
              write( acha2, '(i2.2)' ) ia
                nucl = zcha(1:2)//'-'//acha2//mcha//'  '
            else if( ia .lt. 10 ) then
              write( acha1, '(i1.1)' ) ia
                nucl = zcha(1:2)//'-'//acha1//mcha//'   '
            end if
          end if
        else if( m .eq. 12 ) then
          if( zcha(2:2) .eq. ' ' ) then
            if( ia .gt. 99 ) then
              write( acha3, '(i3.3)' ) ia
                nucl = zcha(1:1)//acha3//mcha//'   '
            else if( ia .gt. 9 .and. ia .lt. 100 ) then
              write( acha2, '(i2.2)' ) ia
                nucl = zcha(1:1)//acha2//mcha//'    '
            else if( ia .lt. 10 ) then
              write( acha1, '(i1.1)' ) ia
                nucl = zcha(1:1)//acha1//mcha//'     '
            end if
          else
            if( ia .gt. 99 ) then
              write( acha3, '(i3.3)' ) ia
                nucl = zcha(1:2)//acha3//mcha//'  '
            else if( ia .gt. 9 .and. ia .lt. 100 ) then
              write( acha2, '(i2.2)' ) ia
                nucl = zcha(1:2)//acha2//mcha//'   '
            else if( ia .lt. 10 ) then
              write( acha1, '(i1.1)' ) ia
                nucl = zcha(1:2)//acha1//mcha//'    '
            end if
          end if
        else if( m .eq. 13 ) then
          if( zcha(2:2) .eq. ' ' ) then
            if( ia .gt. 99 ) then
              write( acha3, '(i3.3)' ) ia
                nucl = acha3//zcha(1:1)//'    '
            else if( ia .gt. 9 .and. ia .lt. 100 ) then
              write( acha2, '(i2.2)' ) ia
                nucl = acha2//zcha(1:1)//'     '
            else if( ia .lt. 10 ) then
              write( acha1, '(i1.1)' ) ia
                nucl = acha1//zcha(1:1)//'      '
            end if
          else
            if( ia .gt. 99 ) then
              write( acha3, '(i3.3)' ) ia
                nucl = acha3//zcha(1:2)//'   '
            else if( ia .gt. 9 .and. ia .lt. 100 ) then
              write( acha2, '(i2.2)' ) ia
                nucl = acha2//zcha(1:2)//'    '
            else if( ia .lt. 10 ) then
              write( acha1, '(i1.1)' ) ia
                nucl = acha1//zcha(1:2)//'     '
            end if
          end if
        else if( m .eq. 14 ) then
          nucl(1:2) = zcha(1:2)
          write( acha3, '(i3)' ) ia
          nucl = '%'//zcha(1:2)//acha3//mcha
        end if


      end if
      return

*-----------------------------------------------------------------------
*     errors
*-----------------------------------------------------------------------

  951 continue
         m_err = 'Unknown or Z > 104 nuclide: '//zcha
         ErrCha = ''
         ErrID = 'L:6297/R:getelm/F:risource.f'
         ierr  = 1
         return

  952 continue
         m_err = 'natural abundance is not allowed: '//nucl
         ErrCha = ''
         ErrID = 'L:6304/R:getelm/F:risource.f'
         ierr  = 1
         return

  953 continue
         m_err = 'anomalous value (atomic number) was detected: '//nucl
         ErrCha = ''
         ErrID = 'L:6311/R:getelm/F:risource.f'
         ierr  = 1
         return

  954 continue
         m_err = 'metastable number ( 0, 1, or 2 ) was wrong: '//mcha
         ErrCha = ''
         ErrID = 'L:6318/R:getelm/F:risource.f'
         ierr  = 1
         return

  955 continue
         m_err = 'unexpected data: '//nucl
         ErrCha = ''
         ErrID = 'L:6325/R:getelm/F:risource.f'
         ierr  = 1
         return

*-----------------------------------------------------------------------

      end

************************************************************************
*                                                                      *
      subroutine caltot()
*                                                                      *
*     Calc. revised totfact and <source> for RI source                 *
*     made by N.Matsuda on 2016/12/01                                  *
*                                                                      *
*                                                                      *
************************************************************************

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'risrcparam.inc'

      common /error/  m_err, l_err, k_err
      character       m_err*200

      common /isomul/ smlwt(isrc), totfact, imsrc
      common /isorsg/ seg0(isrc),seg1(isrc),seg2(isrc),seg3(isrc),
     &                set0(isrc),set1(isrc),set2(isrc),set3(isrc),
     &                jetyp(isrc),jptyp(isrc)

      common /risrc01/ normfact(isrc), asfsum(isrc)
      double precision normfact, asfsum
      common /risrc02/ smlwt2(isrc), totfact2
      double precision karival, smlwt2

*-----------------------------------------------------------------------

      if( imsrc .eq. 1 ) then
        if( jetyp(1) .eq. 28 .or. jetyp(1) .eq. 29 ) then
*         totfact = smlwts
          totfact2 = totfact
          totfact = totfact * normfact(1)
        else
        end if
      else
        do n = 1, imsrc
          smlwt2(n) = smlwt(n)
        end do
          totfact2 = totfact
        do n = 1, imsrc
          if( jetyp(n) .eq. 28 .or. jetyp(n) .eq. 29 ) then
              karival = 0.0d+00
            do nn = 1, imsrc
              karival = karival + smlwt(nn)
            end do
              smlwt(n) = smlwt(n) * normfact(n)
              totfact = totfact
     &              * ( karival - smlwt2(n) + smlwt(n) ) / karival
          else
          end if
        end do

      end if

*-----------------------------------------------------------------------
      return
      end

************************************************************************
*                                                                      *
      subroutine listoutw(iot,j33,mm,nkx,ierr)
*                                                                      *
*     write RI list to OUTPUT file (dummy)                             *
*     made by N.Matsuda on 2016/07/31                                  *
*                                                                      *
*                                                                      *
************************************************************************

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'risrcparam.inc'

      common /error/  m_err, l_err, k_err
      character       m_err*200

      integer iot, j33, mm, nkx
*     NEW parameters ===
      common /risource/ chalct(nuclinmax,isrc), act000(nuclinmax,isrc),
     &                  activy(nuclinmax,isrc), thf(nuclinmax,isrc),
     &                  decayt(isrc)
      real  chalct
      double precision  act000, activy, thf, decayt

*     integer  iot, iodecdc2, ioddc2r, ioddc2b, ioddc2a, ioddc2s

      character chaa*7, rcha*1
      real chza
      double precision rnum, dconst, a0, abq
      integer kno, nmm

      real chal
      character nucl*7, zcha*2, acha3*3, acha2*2, acha1*1, mcha*1

*-----------------------------------------------------------------------
*     file number of scratch file for CHECKING decay calculation
           iot      = 820
*     file number of RIsource.dat or ICRP-07.NDX based on DECDC2
*          iodecdc2 = 811
*     file number of RIsource.rad (radiation data) based on DECDC2
*          ioddc2r  = 812
*     file number of RIsource.bet (beta spectra) based on DECDC2
*          ioddc2b  = 813
*     file number of RIsource.ack (Auger-CK spectra) based on DECDC2
*          ioddc2a  = 814
*     file number of RIsource.dat (decay data) formatted with DCHAIN
*          ioddc2s  = 815

c     err flag for this subroutine
           ierr     =   0
*-----------------------------------------------------------------------

          nmm = mm
      chza = chalct(nkx,j33)
      if( chza .ge. 990000.0 ) then
          chaa = '- SF - '
        if( mm .eq. 11 ) then
          kno = int( chza - 990000.0 )
          nmm = 13
        else if( mm .eq. 21 ) then
          kno = int( chza - 990000.0 )
          nmm = 23
        end if
          abq = 0.0d+00
*     else if( chza .gt. 70000.0 ) then
      ! limit to eliminate Spontaneous fission nuclide
*        go to 926
      else

        call getelm(11,chaa,chza,ierr)

C       Sub. information: getelm(int.,char.,real,int.)
*       you can choose INPUT or OUTPUT style
*         INPUT:
*           2: chaa(ANY)   -> chza,8digit(1001.0(Zaid+))
*         OUTPUT:
*          11: chza(Zaid+) -> chaa,7digit(H-1    )
*          12: chza(Zaid+) -> chaa,7digit(H1     )
*          13: chza(Zaid+) -> chaa,7digit(1H     )
*          14: chza(Zaid+) -> chaa,7digit(%H   1 )

          if( ierr .eq. 1 ) then
*           l_err = ill(jsn)
*           k_err = jsn
*           close(iot)
*             return
            go to 998
          end if

        if( mm .gt. 10 .and. mm .lt. 86 ) then
          if( thf(nkx,j33) .eq. -1.0d+00 ) then
*           dconst = 0.0d+00  ! log(2.)/infini.
            if( mm .eq. 11 ) then
              nmm = nmm + 1
            else if( mm .eq. 21 ) then
              nmm = nmm + 1
            end if
          else
            dconst = dble( log(2.) / thf(nkx,j33) )

            if( mm .eq. 11 ) then
              if( act000(nkx,j33) .eq. 0.0d+00 ) then
                abq = 0.0d+00
                a0  = 0.0d+00
              else
                abq = act000(nkx,j33)
                a0  = abq / dconst
              end if
            elseif( mm .eq. 21 ) then
              if( activy(nkx,j33) .eq. 0.0d+00 ) then
                abq = 0.0d+00
                a0  = 0.0d+00
              else
                abq = activy(nkx,j33) * dconst
                a0  = activy(nkx,j33)
              end if
            end if
          end if
        else
          if( mm .eq. 2 ) then
            abq = act000(nkx,j33)
          else if( mm .eq. 3 ) then
            abq = activy(nkx,j33)
          end if
        end if

      end if

        if( nmm .eq. 1 ) then
C     Ex. OUTPUT  1. Cs-137
          write(iot,'(2x,i2,". ",a7)') nkx, chaa
        else if( nmm .eq. 2 ) then
C     Ex. OUTPUT  1. Cs-137 :  1.000000E+02 (Bq)
          write(iot,'(2x,i2,". ",a7,": ",1pe13.6," (Bq)")')
     &                nkx, chaa, abq
        else if( nmm .eq. 3 ) then
C     Ex. OUTPUT  1. Cs-137 :  1.000000E+02 (Bq)
          write(iot,'(2x,i2,". ",a7,": ",1pe13.6," (Bq)")')
     &                nkx, chaa, abq

        else if( nmm .eq. 11 .or. nmm .eq. 21 ) then
C     Ex. OUTPUT  1. Cs-137
          write(iot,'(2x,i2,". ",a,": ",1pe11.4," (sec) ",
     &                1pe13.6," (Bq) -> ",1pe13.6," (atoms)")')
     &                nkx, chaa, thf(nkx,j33), abq, a0
        else if( nmm .eq. 12 .or. nmm .eq. 22 ) then
C     Ex. OUTPUT  1. Ba-137  ( stable )
          write(iot,'(2x,i2,". ",a7,"  (stable)")') nkx, chaa
        else if( nmm .eq. 13 .or. nmm .eq. 23 ) then
C     Ex. OUTPUT  1. Spontaneous fission nuclide
          write(iot,'(2x,i2,". Spontaneous fission nuclide ("
     &                  ,i2,")")') nkx, kno

        else if( nmm .eq. 89 ) then
C     Stable nuclide
          write(iot,'(4x,a," (stable)")') chaa
        else if( nmm .eq. 90 ) then
C     Registered nuclide
          write(iot,'(4x,"Already registered: ",a)') chaa
        else if( nmm .eq. 99 ) then
C     Warning!! for spontaneous fission nuclide
          write(iot,'("* Warning!! ",a," was SF nclide.")') chaa
        end if

        go to 999
*-----------------------------------------------------------------------
  998 continue
         ierr = 1

  999 continue


*-----------------------------------------------------------------------
      return
      end

************************************************************************
*                                                                      *
      block data decdc2beta0
*                                                                      *
*              type of beta (minus, plus, and  both)                   *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      parameter ( nddc2b_minu = 480 )  ! beta minus
      parameter ( nddc2b_plus = 429 )  ! beta plus
      parameter ( nddc2b_both =  22 )  ! both

      common /ddc2beta0/ tbetas00(3,nddc2b_minu)

*-----------------------------------------------------------------------

      data ( tbetas00(1,i), i = 1, nddc2b_minu ) / ! beta minus (ZAID+)
     &      1003.0,
     &      4010.0,
     &      6014.0,
     &      7016.0,
     &      8019.0,
     &     10024.0,
     &     11024.0,
     &     12027.0,  12028.0,
     &     13028.0,  13029.0,
     &     14031.0,  14032.0,
     &     15032.0,  15033.0,
     &     16035.0,  16037.0,  16038.0,
     &     17038.0,  17039.0,  17040.0,
     &     18039.0,  18041.0,  18042.0,  18043.0,  18044.0,
     &     19042.0,  19043.0,  19044.0,  19045.0,  19046.0,
     &     20045.0,  20047.0,  20049.0,
     &     21046.0,  21047.0,  21048.0,  21049.0,  21050.0,
     &     22051.0,  22052.0,
     &     23050.0,  23052.0,  23053.0,
     &     24055.0,  24056.0,
     &     25056.0,  25057.0,  25058.1,
     &     26059.0,  26060.0,  26061.0,  26062.0,
     &     27060.0,  27060.1,  27061.0,  27062.0,  27062.1,
     &     28063.0,  28065.0,  28066.0,
     &     29066.0,  29067.0,  29069.0,
     &     30069.0,  30069.1,  30071.0,  30071.1,  30072.0,
     &     31070.0,  31072.0,  31073.0,  31074.0,
     &     32075.0,  32077.0,  32078.0,
     &     33076.0,  33077.0,  33078.0,  33079.0,
     &     34079.0,  34079.1,  34081.0,  34081.1,  34083.0,  34083.1,
     &     34084.0,
     &     35082.0,  35082.1,  35083.0,  35084.0,  35084.1,  35085.0,
     &     36085.0,  36085.1,  36087.0,  36088.0,  36089.0,
     &     37086.0,  37087.0,  37088.0,  37089.0,  37090.0,  37090.1,
     &     38089.0,  38090.0,  38091.0,  38092.0,  38093.0,  38094.0,
     &     39090.0,  39090.1,  39091.0,  39092.0,  39093.0,  39094.0,
     &     39095.0,
     &     40093.0,  40095.0,  40097.0,
     &     41094.0,  41094.1,  41095.0,  41095.1,  41096.0,  41097.0,
     &     41098.1,  41099.0,  41099.1,
     &     42099.0,  42101.0,  42102.0,
     &     43098.0,  43099.0,  43099.1,  43101.0,  43102.0,  43102.1,
     &     43104.0,  43105.0,
     &     44103.0,  44105.0,  44106.0,  44107.0,  44108.0,
     &     45104.0,  45104.1,  45105.0,  45106.0,  45106.1,  45107.0,
     &     45108.0,  45109.0,
     &     46107.0,  46109.0,  46111.0,  46112.0,  46114.0,
     &     47110.0,  47110.1,  47111.0,  47111.1,  47112.0,  47113.0,
     &     47113.1,  47114.0,  47115.0,  47116.0,  47117.0,
     &     48113.0,  48113.1,  48115.0,  48115.1,  48117.0,  48117.1,
     &     48118.0,  48119.0,  48119.1,
     &     49115.0,  49115.1,  49116.1,  49117.0,  49117.1,  49118.0,
     &     49118.1,  49119.0,  49119.1,  49121.0,  49121.1,
     &     50121.0,  50121.1,  50123.0,  50123.1,  50125.0,  50125.1,
     &     50126.0,  50127.0,  50127.1,  50128.0,  50129.0,  50130.0,
     &     50130.1,
     &     51124.0,  51124.1,  51125.0,  51126.0,  51126.1,  51127.0,
     &     51128.0,  51128.1,  51129.0,  51130.0,  51130.1,  51131.0,
     &     51133.0,
     &     52127.0,  52127.1,  52129.0,  52129.1,  52131.0,  52131.1,
     &     52132.0,  52133.0,  52133.1,  52134.0,
     &     53129.0,  53130.0,  53130.1,  53131.0,  53132.0,  53132.1,
     &     53133.0,  53134.0,  53134.1,  53135.0,
     &     54133.0,  54135.0,  54135.1,  54137.0,  54138.0,
     &     55134.0,  55135.0,  55136.0,  55137.0,  55138.0,  55138.1,
     &     55139.0,  55140.0,
     &     56139.0,  56140.0,  56141.0,  56142.0,
     &     57138.0,  57140.0,  57141.0,  57142.0,  57143.0,
     &     58141.0,  58143.0,  58144.0,  58145.0,
     &     59142.0,  59143.0,  59144.0,  59144.1,  59145.0,  59146.0,
     &     59147.0,  59148.0,  59148.1,
     &     60147.0,  60149.0,  60151.0,  60152.0,  61146.0,
     &     61147.0,  61148.0,  61148.1,  61149.0,  61150.0,  61151.0,
     &     61152.0,  61152.1,  61153.0,  61154.0,  61154.1,
     &     62151.0,  62153.0,  62155.0,  62156.0,  62157.0,
     &     63154.0,  63155.0,  63156.0,  63157.0,  63158.0,  63159.0,
     &     64159.0,  64162.0,
     &     65158.0,  65160.0,  65161.0,  65162.0,  65163.0,  65164.0,
     &     65165.0,
     &     66165.0,  66165.1,  66166.0,  66167.0,  66168.0,
     &     67164.0,  67166.0,  67166.1,  67167.0,  67168.0,  67170.0,
     &     68169.0,  68171.0,  68172.0,  68173.0,
     &     69170.0,  69171.0,  69172.0,  69173.0,  69174.0,  69175.0,
     &     69176.0,
     &     70175.0,  70177.0,  70178.0,  70179.0,
     &     71176.0,  71176.1,  71177.0,  71177.1,  71178.0,  71178.1,
     &     71179.0,  71180.0,  71181.0,
     &     72180.1,  72181.0,  72182.0,  72182.1,  72183.0,  72184.0,
     &     73180.0,  73182.0,  73183.0,  73184.0,  73185.0,  73186.0,
     &     74185.0,  74187.0,  74188.0,  74190.0,
     &     75186.0,  75187.0,  75188.0,  75189.0,  75190.0,  75190.1,
     &     76191.0,  76193.0,  76194.0,  76196.0,
     &     77192.0,  77192.1,  77194.0,  77194.1,  77195.0,  77195.1,
     &     77196.0,  77196.1,
     &     78197.0,  78197.1,  78199.0,  78200.0,  78202.0,
     &     79196.0,  79198.0,  79199.0,  79200.0,  79200.1,  79201.0,
     &     79202.0,
     &     80203.0,  80205.0,  80206.0,  80207.0,
     &     81204.0,  81206.0,  81207.0,  81208.0,  81209.0,  81210.0,
     &     82209.0,  82210.0,  82211.0,  82212.0,  82214.0,
     &     83210.0,  83211.0,  83212.0,  83212.2,  83213.0,  83214.0,
     &     83215.0,  83216.0,
     &     84218.0,
     &     85218.0,  85220.0,
     &     86223.0,
     &     87220.0,  87222.0,  87223.0,  87224.0,  87227.0,
     &     88225.0,  88227.0,  88228.0,  88230.0,
     &     89226.0,  89227.0,  89228.0,  89230.0,  89231.0,  89232.0,
     &     89233.0,
     &     90231.0,  90233.0,  90234.0,  90235.0,  90236.0,
     &     91230.0,  91232.0,  91233.0,  91234.0,  91234.1,  91235.0,
     &     91236.0,  91237.0,
     &     92237.0,  92239.0,  92240.0,  92242.0,
     &     93236.0,  93236.1,  93238.0,  93239.0,  93240.0,  93240.1,
     &     93241.0,  93242.0,  93242.1,
     &     94241.0,  94243.0,  94245.0,  94246.0,
     &     95242.0,  95244.0,  95244.1,  95245.0,  95246.0,  95246.1,
     &     95247.0,
     &     96249.0,  96250.0,  96251.0,
     &     97248.1,  97249.0,  97250.0,  97251.0,
     &     98253.0,  98255.0,
     &     99254.0,  99254.1,  99255.0,  99256.0/

*-----------------------------------------------------------------------
      data ( tbetas00(2,i), i = 1, nddc2b_plus ) / ! beta plus (ZAID+)
     &      6010.0,   6011.0,
     &      7013.0,
     &      8014.0,   8015.0,
     &      9017.0,   9018.0,
     &     10019.0,
     &     11022.0,
     &     13026.0,
     &     15030.0,
     &     17034.0,  17034.1,
     &     19038.0,
     &     21042.1,  21043.0,  21044.0,
     &     22045.0,
     &     23047.0,  23048.0,
     &     24048.0,  24049.0,
     &     25050.1,  25051.0,  25052.0,  25052.1,
     &     26052.0,  26053.0,
     &     27054.1,  27055.0,  27056.0,  27058.0,
     &     28056.0,  28057.0,  28059.0,
     &     29057.0,  29059.0,  29060.0,  29061.0,  29062.0,
     &     30060.0,  30061.0,  30062.0,  30063.0,  30065.0,
     &     31064.0,  31065.0,  31066.0,  31068.0,
     &     32066.0,  32067.0,  32069.0,
     &     33068.0,  33069.0,  33070.0,  33071.0,  33072.0,
     &     34070.0,  34071.0,  34073.0,  34073.1,
     &     35072.0,  35073.0,  35074.0,  35074.1,  35075.0,  35076.0,
     &     35076.1,  35077.0,
     &     36074.0,  36075.0,  36077.0,  36079.0,
     &     37077.0,  37078.0,  37078.1,  37079.0,  37080.0,  37081.0,
     &     37081.1,  37082.0,  37082.1,
     &     38079.0,  38080.0,  38081.0,  38083.0,  38085.1,
     &     39081.0,  39083.0,  39083.1,  39084.1,  39085.0,  39085.1,
     &     39086.0,  39086.1,  39087.0,  39087.1,  39088.0,
     &     40085.0,  40086.0,  40087.0,  40089.0,  40089.1,
     &     41087.0,  41088.0,  41088.1,  41089.0,  41089.1,  41090.0,
     &     41091.0,  41091.1,  41092.1,
     &     42089.0,  42090.0,  42091.0,  42091.1,
     &     43091.0,  43091.1,  43092.0,  43093.0,  43093.1,  43094.0,
     &     43094.1,  43095.1,  43096.1,
     &     44092.0,  44094.0,  44095.0,
     &     45094.0,  45095.0,  45095.1,  45096.0,  45096.1,  45097.0,
     &     45097.1,  45098.0,  45099.0,  45099.1,
     &     45100.0,  45100.1,  45102.1,
     &     46096.0,  46097.0,  46098.0,  46099.0,  46101.0,
     &     47099.0,  47100.1,  47101.0,  47102.0,  47102.1,  47103.0,
     &     47104.0,  47104.1,  47105.1,
     &     48101.0,  48102.0,  48103.0,  48105.0,  48107.0,
     &     49103.0,  49105.0,  49106.0,  49106.1,  49107.0,  49108.0,
     &     49108.1,  49109.0,  49110.0,  49110.1,
     &     50106.0,  50108.0,  50109.0,  50111.0,
     &     51111.0,  51113.0,  51114.0,  51115.0,  51116.0,  51116.1,
     &     51117.0,  51118.0,  51118.1,  51120.0,
     &     52113.0,  52114.0,  52115.0,  52115.1,  52116.0,  52117.0,
     &     52119.0,  52119.1,
     &     53118.0,  53118.1,  53119.0,  53120.0,  53120.1,  53121.0,
     &     53122.0,  53124.0,
     &     54120.0,  54121.0,  54123.0,  54125.0,
     &     55121.0,  55121.1,  55123.0,  55124.0,  55125.0,  55126.0,
     &     55127.0,  55128.0,  55129.0,
     &     56124.0,  56126.0,  56127.0,  56129.0,  56129.1,
     &     57128.0,  57129.0,  57130.0,  57131.0,  57132.0,  57132.1,
     &     57133.0,  57134.0,  57135.0,  57136.0,
     &     58130.0,  58131.0,  58133.0,  58133.1,  58135.0,  58137.0,
     &     59134.0,  59134.1,  59135.0,  59136.0,  59137.0,  59138.0,
     &     59138.1,  59139.0,  59140.0,
     &     60134.0,  60135.0,  60136.0,  60137.0,  60139.0,  60139.1,
     &     60141.0,  60141.1,
     &     61136.0,  61137.1,  61139.0,  61140.0,  61140.1,  61141.0,
     &     61142.0,
     &     62139.0,  62140.0,  62141.0,  62141.1,  62142.0,  62143.0,
     &     62143.1,
     &     63142.0,  63142.1,  63143.0,  63144.0,  63145.0,  63146.0,
     &     63147.0,  63148.0,  63150.0,
     &     64142.0,  64143.1,  64144.0,  64145.0,  64145.1,  64147.0,
     &     64149.0,
     &     65146.0,  65147.0,  65147.1,  65148.0,  65148.1,  65149.0,
     &     65149.1,  65150.0,  65150.1,  65151.0,  65151.1,  65152.0,
     &     65152.1,  65153.0,  65154.0,
     &     66148.0,  66149.0,  66150.0,  66151.0,  66153.0,  66155.0,
     &     67150.0,  67153.0,  67153.1,  67154.0,  67154.1,  67155.0,
     &     67156.0,  67157.0,  67159.0,  67160.0,  67162.0,  67162.1,
     &     68154.0,  68159.0,  68161.0,  68163.0,  69161.0,
     &     69162.0,  69163.0,  69164.0,  69165.0,  69166.0,
     &     70162.0,  70163.0,  70165.0,  70167.0,
     &     71165.0,  71167.0,  71169.0,  71170.0,  71171.0,  71172.0,
     &     71174.0,
     &     72167.0,  72169.0,  72173.0,
     &     73170.0,  73172.0,  73173.0,  73174.0,  73175.0,  73176.0,
     &     73178.0,
     &     74177.0,
     &     75178.0,  75179.0,  75180.0,  75181.0,  75182.1,  75184.0,
     &     76180.0,  76181.0,  76183.0,  76183.1,
     &     77180.0,  77182.0,  77183.0,  77184.0,  77185.0,  77186.0,
     &     77186.1,  77187.0,  77188.0,
     &     78184.0,  78187.0,  78189.0,
     &     79186.0,  79187.0,  79190.0,  79191.0,  79192.0,  79194.0,
     &     80190.0,  80191.1,  80193.0,  80193.1,  80195.0,  80195.1,
     &     81190.0,  81190.1,  81194.0,  81194.1,  81195.0,  81196.0,
     &     81197.0,  81198.0,  81198.1,  81199.0,
     &     81200.0,  82194.0,  82195.1,  82196.0,  82197.0,  82197.1,
     &     82199.0,
     &     82201.0,
     &     83197.0,  83200.0,  83201.0,  83202.0,  83203.0,  83204.0,
     &     83205.0,  83206.0,  83207.0,
     &     84203.0,  84205.0,  84207.0,
     &     85204.0,  85205.0,  85206.0,  85207.0,  85208.0,  85209.0,
     &     85210.0,
     &     86207.0,  86209.0,  86211.0,
     &     87212.0,
     &     91228.0,
     &     93232.0,  93234.0,
     &     95238.0,
     &     96239.0,
     &     99249.0,  99250.1,
     &    100251.0/

*-----------------------------------------------------------------------
      data ( tbetas00(3,i), i = 1, nddc2b_both ) / ! both (ZAID+)
     &     17036.0,
     &     19040.0,
     &     25054.0,
     &     29064.0,
     &     33074.0,
     &     35078.0,  35080.0,
     &     37084.0,
     &     45102.0,
     &     47106.0,  47108.0,
     &     49112.0,  49114.0,
     &     51122.0,
     &     53126.0,  53128.0,
     &     55130.0,  55132.0,
     &     63150.1,  63152.0,  63152.1,
     &     69168.0/

*-----------------------------------------------------------------------

      end


************************************************************************
*                                                                      *
      block data decdc2beta3
*                                                                      *
*              picked up beta spectrum from DECDC2 (3: both)           *
*              calculated by A. Endo
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      parameter ( nddc2bbemax = 120 )  ! Maximum number of energy mesh
      parameter ( nddc2b_both =  22 )  ! both

      common /ddc2beta3/ ddc2beng(nddc2bbemax), ddc2bend(2*nddc2b_both),
     &                   ddc2byld(2*nddc2b_both,nddc2bbemax),
     &                   ddc2bnuc(nddc2b_both)
      common /ddc2beta4/ nnddc2bn(2*nddc2b_both)

      character ddc2bnuc*7
*     integer   nnddc2bn

*-----------------------------------------------------------------------

      data  ( ddc2beng(i), i = 1, nddc2bbemax ) /
     &  0.000000, 0.000100, 0.000110, 0.000120, 0.000130, 0.000140,
     &  0.000150, 0.000160, 0.000180, 0.000200, 0.000220, 0.000240,
     &  0.000260, 0.000280, 0.000300, 0.000320, 0.000360, 0.000400,
     &  0.000450, 0.000500, 0.000550, 0.000600, 0.000650, 0.000700,
     &  0.000750, 0.000800, 0.000850, 0.000900, 0.001000, 0.001100,
     &  0.001200, 0.001300, 0.001400, 0.001500, 0.001600, 0.001800,
     &  0.002000, 0.002200, 0.002400, 0.002600, 0.002800, 0.003000,
     &  0.003200, 0.003600, 0.004000, 0.004500, 0.005000, 0.005500,
     &  0.006000, 0.006500, 0.007000, 0.007500, 0.008000, 0.008500,
     &  0.009000, 0.010000, 0.011000, 0.012000, 0.013000, 0.014000,
     &  0.015000, 0.016000, 0.018000, 0.020000, 0.022000, 0.024000,
     &  0.026000, 0.028000, 0.030000, 0.032000, 0.036000, 0.040000,
     &  0.045000, 0.050000, 0.055000, 0.060000, 0.065000, 0.070000,
     &  0.075000, 0.080000, 0.085000, 0.090000, 0.100000, 0.110000,
     &  0.120000, 0.130000, 0.140000, 0.150000, 0.160000, 0.180000,
     &  0.200000, 0.220000, 0.240000, 0.260000, 0.280000, 0.300000,
     &  0.320000, 0.360000, 0.400000, 0.450000, 0.500000, 0.550000,
     &  0.600000, 0.650000, 0.700000, 0.750000, 0.800000, 0.850000,
     &  0.900000, 1.000000, 1.100000, 1.200000, 1.300000, 1.400000,
     &  1.500000, 1.600000, 1.800000, 2.000000, 2.200000, 2.400000/

      data  ( ddc2bend(i), i = 1, 2*nddc2b_both ) /
     &  0.708600, 1.311090, 0.697100, 0.578700, 1.353000, 0.272000,
     &  2.001000, 0.894000, 1.151000, 0.195000, 1.649000, 0.663000,
     &  1.988700, 1.978600, 1.258000, 2.119000, 0.361000, 0.815500,
     &  1.013100, 1.474520, 1.864400, 0.169270,
     &  0.120070, 0.482870, 0.355100, 0.653100, 1.540400, 2.552000,
     &  0.849100, 1.658900, 1.301000, 1.943300, 0.896000, 1.564000,
     &  0.430000, 0.597700, 1.133000, 0.230000, 1.957000, 0.440280,
     &  1.281100, 0.730520, 0.897900, 0.577200/

      data  ( nnddc2bn(i), i = 1, 2*nddc2b_both ) /
     &  105, 113, 104, 102, 113,  94, 118, 108, 111,  90,
     &  116, 104, 117, 117, 112, 118,  98, 107, 110, 114,
     &  117,  89,
     &   85, 100,  97, 104, 115, 120, 107, 116, 113, 117,
     &  108, 115,  99, 102, 111,  92, 117,  99, 112, 105,
     &  108, 102/

      data  ( ddc2bnuc(i), i = 1, nddc2b_both ) /
     &  'Cl-36  ', 'K-40   ', 'Mn-54  ', 'Cu-64  ', 'As-74  ',
     &  'Br-78  ', 'Br-80  ', 'Rb-84  ', 'Rh-102 ', 'Ag-106 ',
     &  'Ag-108 ', 'In-112 ', 'In-114 ', 'Sb-122 ', 'I-126  ',
     &  'I-128  ', 'Cs-130 ', 'Cs-132 ', 'Eu-150m', 'Eu-152 ',
     &  'Eu-152m', 'Tm-168'/

      data ( ddc2byld(1,i), i = 1, 105 ) / ! Cl-36 beta minus
     &  1.434E+00, 1.434E+00, 1.434E+00, 1.434E+00, 1.434E+00,
     &  1.434E+00, 1.433E+00, 1.433E+00, 1.433E+00, 1.433E+00,
     &  1.433E+00, 1.433E+00, 1.433E+00, 1.433E+00, 1.433E+00,
     &  1.433E+00, 1.432E+00, 1.432E+00, 1.432E+00, 1.432E+00,
     &  1.431E+00, 1.431E+00, 1.431E+00, 1.431E+00, 1.430E+00,
     &  1.430E+00, 1.430E+00, 1.430E+00, 1.429E+00, 1.428E+00,
     &  1.428E+00, 1.427E+00, 1.427E+00, 1.426E+00, 1.426E+00,
     &  1.425E+00, 1.424E+00, 1.424E+00, 1.424E+00, 1.424E+00,
     &  1.424E+00, 1.424E+00, 1.424E+00, 1.424E+00, 1.424E+00,
     &  1.424E+00, 1.424E+00, 1.424E+00, 1.425E+00, 1.426E+00,
     &  1.427E+00, 1.428E+00, 1.429E+00, 1.430E+00, 1.432E+00,
     &  1.435E+00, 1.438E+00, 1.442E+00, 1.447E+00, 1.451E+00,
     &  1.456E+00, 1.460E+00, 1.470E+00, 1.480E+00, 1.491E+00,
     &  1.501E+00, 1.512E+00, 1.522E+00, 1.532E+00, 1.542E+00,
     &  1.562E+00, 1.582E+00, 1.605E+00, 1.627E+00, 1.648E+00,
     &  1.669E+00, 1.688E+00, 1.706E+00, 1.724E+00, 1.741E+00,
     &  1.757E+00, 1.772E+00, 1.800E+00, 1.826E+00, 1.849E+00,
     &  1.870E+00, 1.889E+00, 1.905E+00, 1.920E+00, 1.945E+00,
     &  1.962E+00, 1.973E+00, 1.976E+00, 1.972E+00, 1.961E+00,
     &  1.941E+00, 1.912E+00, 1.825E+00, 1.696E+00, 1.471E+00,
     &  1.180E+00, 8.372E-01, 4.798E-01, 1.699E-01, 4.424E-03/

      data ( ddc2byld(2,i), i = 1, 113 ) / ! K-40 beta minus
     &  4.075E-01, 4.073E-01, 4.073E-01, 4.073E-01, 4.072E-01,
     &  4.072E-01, 4.072E-01, 4.072E-01, 4.072E-01, 4.071E-01,
     &  4.071E-01, 4.071E-01, 4.070E-01, 4.070E-01, 4.070E-01,
     &  4.069E-01, 4.069E-01, 4.068E-01, 4.067E-01, 4.066E-01,
     &  4.065E-01, 4.065E-01, 4.064E-01, 4.063E-01, 4.062E-01,
     &  4.061E-01, 4.060E-01, 4.060E-01, 4.058E-01, 4.056E-01,
     &  4.055E-01, 4.053E-01, 4.051E-01, 4.050E-01, 4.048E-01,
     &  4.045E-01, 4.041E-01, 4.038E-01, 4.038E-01, 4.039E-01,
     &  4.041E-01, 4.042E-01, 4.043E-01, 4.045E-01, 4.047E-01,
     &  4.050E-01, 4.053E-01, 4.056E-01, 4.060E-01, 4.064E-01,
     &  4.068E-01, 4.072E-01, 4.077E-01, 4.082E-01, 4.087E-01,
     &  4.099E-01, 4.111E-01, 4.125E-01, 4.139E-01, 4.154E-01,
     &  4.169E-01, 4.185E-01, 4.218E-01, 4.253E-01, 4.289E-01,
     &  4.325E-01, 4.361E-01, 4.399E-01, 4.436E-01, 4.473E-01,
     &  4.548E-01, 4.623E-01, 4.716E-01, 4.809E-01, 4.900E-01,
     &  4.990E-01, 5.078E-01, 5.166E-01, 5.252E-01, 5.338E-01,
     &  5.422E-01, 5.505E-01, 5.667E-01, 5.826E-01, 5.980E-01,
     &  6.131E-01, 6.279E-01, 6.422E-01, 6.563E-01, 6.834E-01,
     &  7.092E-01, 7.337E-01, 7.570E-01, 7.790E-01, 7.998E-01,
     &  8.193E-01, 8.376E-01, 8.702E-01, 8.976E-01, 9.244E-01,
     &  9.428E-01, 9.527E-01, 9.540E-01, 9.466E-01, 9.303E-01,
     &  9.051E-01, 8.706E-01, 8.265E-01, 7.719E-01, 6.266E-01,
     &  4.255E-01, 1.806E-01, 2.842E-03/

      data ( ddc2byld(3,i), i = 1, 104 ) / ! Mn-54 beta minus
     &  3.951E-06, 3.949E-06, 3.949E-06, 3.948E-06, 3.948E-06,
     &  3.948E-06, 3.948E-06, 3.947E-06, 3.947E-06, 3.946E-06,
     &  3.946E-06, 3.945E-06, 3.945E-06, 3.944E-06, 3.944E-06,
     &  3.943E-06, 3.942E-06, 3.941E-06, 3.940E-06, 3.939E-06,
     &  3.937E-06, 3.936E-06, 3.935E-06, 3.934E-06, 3.932E-06,
     &  3.931E-06, 3.930E-06, 3.928E-06, 3.926E-06, 3.923E-06,
     &  3.921E-06, 3.918E-06, 3.916E-06, 3.913E-06, 3.911E-06,
     &  3.905E-06, 3.900E-06, 3.895E-06, 3.890E-06, 3.885E-06,
     &  3.880E-06, 3.875E-06, 3.870E-06, 3.872E-06, 3.874E-06,
     &  3.877E-06, 3.879E-06, 3.881E-06, 3.884E-06, 3.886E-06,
     &  3.889E-06, 3.891E-06, 3.894E-06, 3.896E-06, 3.899E-06,
     &  3.905E-06, 3.911E-06, 3.917E-06, 3.924E-06, 3.931E-06,
     &  3.939E-06, 3.947E-06, 3.963E-06, 3.981E-06, 4.000E-06,
     &  4.020E-06, 4.041E-06, 4.062E-06, 4.083E-06, 4.105E-06,
     &  4.150E-06, 4.195E-06, 4.252E-06, 4.310E-06, 4.368E-06,
     &  4.425E-06, 4.482E-06, 4.538E-06, 4.593E-06, 4.647E-06,
     &  4.701E-06, 4.753E-06, 4.855E-06, 4.952E-06, 5.044E-06,
     &  5.132E-06, 5.214E-06, 5.291E-06, 5.363E-06, 5.489E-06,
     &  5.592E-06, 5.672E-06, 5.728E-06, 5.760E-06, 5.767E-06,
     &  5.749E-06, 5.705E-06, 5.536E-06, 5.249E-06, 4.700E-06,
     &  3.910E-06, 2.874E-06, 1.663E-06, 5.217E-07/

      data ( ddc2byld(4,i), i = 1, 102 ) / ! Cu-64 beta minus
     &  1.108E+00, 1.108E+00, 1.108E+00, 1.108E+00, 1.108E+00,
     &  1.108E+00, 1.107E+00, 1.107E+00, 1.107E+00, 1.107E+00,
     &  1.107E+00, 1.107E+00, 1.107E+00, 1.107E+00, 1.107E+00,
     &  1.107E+00, 1.107E+00, 1.107E+00, 1.107E+00, 1.107E+00,
     &  1.107E+00, 1.107E+00, 1.107E+00, 1.107E+00, 1.107E+00,
     &  1.106E+00, 1.106E+00, 1.106E+00, 1.106E+00, 1.106E+00,
     &  1.106E+00, 1.106E+00, 1.105E+00, 1.105E+00, 1.105E+00,
     &  1.105E+00, 1.104E+00, 1.104E+00, 1.104E+00, 1.103E+00,
     &  1.103E+00, 1.103E+00, 1.102E+00, 1.102E+00, 1.101E+00,
     &  1.101E+00, 1.101E+00, 1.101E+00, 1.101E+00, 1.101E+00,
     &  1.101E+00, 1.101E+00, 1.101E+00, 1.101E+00, 1.101E+00,
     &  1.101E+00, 1.101E+00, 1.101E+00, 1.101E+00, 1.102E+00,
     &  1.102E+00, 1.102E+00, 1.103E+00, 1.103E+00, 1.104E+00,
     &  1.105E+00, 1.106E+00, 1.107E+00, 1.109E+00, 1.110E+00,
     &  1.113E+00, 1.116E+00, 1.119E+00, 1.122E+00, 1.125E+00,
     &  1.128E+00, 1.130E+00, 1.132E+00, 1.133E+00, 1.134E+00,
     &  1.134E+00, 1.134E+00, 1.132E+00, 1.129E+00, 1.123E+00,
     &  1.115E+00, 1.106E+00, 1.094E+00, 1.081E+00, 1.048E+00,
     &  1.009E+00, 9.636E-01, 9.124E-01, 8.562E-01, 7.955E-01,
     &  7.313E-01, 6.642E-01, 5.251E-01, 3.857E-01, 2.238E-01,
     &  9.304E-02, 1.367E-02/

      data ( ddc2byld(5,i), i = 1, 113 ) / ! As-74 beta minus
     &  4.825E-01, 4.824E-01, 4.824E-01, 4.824E-01, 4.824E-01,
     &  4.824E-01, 4.824E-01, 4.824E-01, 4.824E-01, 4.824E-01,
     &  4.824E-01, 4.824E-01, 4.823E-01, 4.823E-01, 4.823E-01,
     &  4.823E-01, 4.823E-01, 4.822E-01, 4.822E-01, 4.822E-01,
     &  4.821E-01, 4.821E-01, 4.821E-01, 4.820E-01, 4.820E-01,
     &  4.820E-01, 4.819E-01, 4.819E-01, 4.818E-01, 4.817E-01,
     &  4.817E-01, 4.816E-01, 4.815E-01, 4.815E-01, 4.814E-01,
     &  4.813E-01, 4.811E-01, 4.810E-01, 4.808E-01, 4.807E-01,
     &  4.806E-01, 4.804E-01, 4.803E-01, 4.800E-01, 4.797E-01,
     &  4.794E-01, 4.794E-01, 4.794E-01, 4.795E-01, 4.795E-01,
     &  4.795E-01, 4.796E-01, 4.796E-01, 4.797E-01, 4.797E-01,
     &  4.798E-01, 4.799E-01, 4.800E-01, 4.801E-01, 4.802E-01,
     &  4.803E-01, 4.805E-01, 4.808E-01, 4.811E-01, 4.815E-01,
     &  4.820E-01, 4.825E-01, 4.830E-01, 4.836E-01, 4.843E-01,
     &  4.856E-01, 4.871E-01, 4.890E-01, 4.910E-01, 4.931E-01,
     &  4.953E-01, 4.974E-01, 4.997E-01, 5.019E-01, 5.042E-01,
     &  5.064E-01, 5.086E-01, 5.126E-01, 5.161E-01, 5.192E-01,
     &  5.216E-01, 5.236E-01, 5.250E-01, 5.259E-01, 5.261E-01,
     &  5.242E-01, 5.203E-01, 5.145E-01, 5.069E-01, 4.975E-01,
     &  4.865E-01, 4.741E-01, 4.452E-01, 4.121E-01, 3.666E-01,
     &  3.191E-01, 2.728E-01, 2.310E-01, 1.974E-01, 1.760E-01,
     &  1.671E-01, 1.578E-01, 1.468E-01, 1.340E-01, 1.032E-01,
     &  6.716E-02, 3.103E-02, 4.681E-03/

      data ( ddc2byld(6,i), i = 1, 94 ) / ! Br-78 beta minus
     &  8.531E-04, 8.527E-04, 8.526E-04, 8.526E-04, 8.525E-04,
     &  8.525E-04, 8.524E-04, 8.524E-04, 8.523E-04, 8.522E-04,
     &  8.521E-04, 8.520E-04, 8.519E-04, 8.519E-04, 8.518E-04,
     &  8.517E-04, 8.515E-04, 8.513E-04, 8.511E-04, 8.508E-04,
     &  8.506E-04, 8.504E-04, 8.502E-04, 8.499E-04, 8.497E-04,
     &  8.495E-04, 8.492E-04, 8.490E-04, 8.485E-04, 8.481E-04,
     &  8.476E-04, 8.472E-04, 8.467E-04, 8.463E-04, 8.458E-04,
     &  8.449E-04, 8.440E-04, 8.430E-04, 8.421E-04, 8.412E-04,
     &  8.403E-04, 8.394E-04, 8.384E-04, 8.366E-04, 8.348E-04,
     &  8.325E-04, 8.302E-04, 8.286E-04, 8.269E-04, 8.252E-04,
     &  8.235E-04, 8.218E-04, 8.201E-04, 8.184E-04, 8.167E-04,
     &  8.133E-04, 8.099E-04, 8.065E-04, 8.030E-04, 7.996E-04,
     &  7.961E-04, 7.927E-04, 7.858E-04, 7.789E-04, 7.720E-04,
     &  7.652E-04, 7.584E-04, 7.515E-04, 7.447E-04, 7.379E-04,
     &  7.243E-04, 7.105E-04, 6.933E-04, 6.758E-04, 6.582E-04,
     &  6.403E-04, 6.222E-04, 6.038E-04, 5.852E-04, 5.664E-04,
     &  5.473E-04, 5.281E-04, 4.892E-04, 4.498E-04, 4.101E-04,
     &  3.705E-04, 3.312E-04, 2.925E-04, 2.548E-04, 1.832E-04,
     &  1.193E-04, 6.605E-05, 2.649E-05, 3.938E-06/

      data ( ddc2byld(7,i), i = 1, 118 ) / ! Br-80 beta minus
     &  3.227E-01, 3.227E-01, 3.227E-01, 3.227E-01, 3.227E-01,
     &  3.227E-01, 3.227E-01, 3.228E-01, 3.228E-01, 3.228E-01,
     &  3.228E-01, 3.228E-01, 3.228E-01, 3.228E-01, 3.228E-01,
     &  3.228E-01, 3.228E-01, 3.228E-01, 3.228E-01, 3.228E-01,
     &  3.229E-01, 3.229E-01, 3.229E-01, 3.229E-01, 3.229E-01,
     &  3.229E-01, 3.229E-01, 3.230E-01, 3.230E-01, 3.230E-01,
     &  3.230E-01, 3.231E-01, 3.231E-01, 3.231E-01, 3.232E-01,
     &  3.232E-01, 3.233E-01, 3.233E-01, 3.234E-01, 3.234E-01,
     &  3.235E-01, 3.235E-01, 3.236E-01, 3.237E-01, 3.238E-01,
     &  3.239E-01, 3.241E-01, 3.245E-01, 3.249E-01, 3.252E-01,
     &  3.256E-01, 3.260E-01, 3.264E-01, 3.268E-01, 3.271E-01,
     &  3.279E-01, 3.286E-01, 3.294E-01, 3.302E-01, 3.309E-01,
     &  3.317E-01, 3.325E-01, 3.340E-01, 3.356E-01, 3.373E-01,
     &  3.389E-01, 3.406E-01, 3.423E-01, 3.441E-01, 3.458E-01,
     &  3.495E-01, 3.532E-01, 3.579E-01, 3.627E-01, 3.676E-01,
     &  3.725E-01, 3.774E-01, 3.824E-01, 3.875E-01, 3.925E-01,
     &  3.975E-01, 4.026E-01, 4.127E-01, 4.227E-01, 4.327E-01,
     &  4.426E-01, 4.525E-01, 4.622E-01, 4.718E-01, 4.907E-01,
     &  5.091E-01, 5.270E-01, 5.443E-01, 5.609E-01, 5.770E-01,
     &  5.923E-01, 6.070E-01, 6.343E-01, 6.586E-01, 6.847E-01,
     &  7.060E-01, 7.223E-01, 7.337E-01, 7.403E-01, 7.421E-01,
     &  7.392E-01, 7.313E-01, 7.187E-01, 7.014E-01, 6.540E-01,
     &  5.918E-01, 5.184E-01, 4.380E-01, 3.557E-01, 2.728E-01,
     &  1.920E-01, 5.743E-02, 1.668E-06/

      data ( ddc2byld(8,i), i = 1, 108 ) / ! Rb-84 beta minus
     &  5.947E-02, 5.946E-02, 5.946E-02, 5.945E-02, 5.945E-02,
     &  5.945E-02, 5.945E-02, 5.945E-02, 5.945E-02, 5.944E-02,
     &  5.944E-02, 5.944E-02, 5.943E-02, 5.943E-02, 5.943E-02,
     &  5.943E-02, 5.942E-02, 5.941E-02, 5.941E-02, 5.940E-02,
     &  5.939E-02, 5.938E-02, 5.938E-02, 5.937E-02, 5.936E-02,
     &  5.936E-02, 5.935E-02, 5.934E-02, 5.933E-02, 5.931E-02,
     &  5.930E-02, 5.928E-02, 5.927E-02, 5.925E-02, 5.924E-02,
     &  5.921E-02, 5.918E-02, 5.915E-02, 5.912E-02, 5.909E-02,
     &  5.906E-02, 5.903E-02, 5.900E-02, 5.895E-02, 5.889E-02,
     &  5.881E-02, 5.874E-02, 5.867E-02, 5.869E-02, 5.870E-02,
     &  5.870E-02, 5.870E-02, 5.870E-02, 5.870E-02, 5.870E-02,
     &  5.871E-02, 5.871E-02, 5.871E-02, 5.871E-02, 5.872E-02,
     &  5.872E-02, 5.872E-02, 5.873E-02, 5.874E-02, 5.875E-02,
     &  5.877E-02, 5.878E-02, 5.881E-02, 5.883E-02, 5.886E-02,
     &  5.893E-02, 5.900E-02, 5.910E-02, 5.921E-02, 5.932E-02,
     &  5.944E-02, 5.955E-02, 5.966E-02, 5.978E-02, 5.989E-02,
     &  5.999E-02, 6.009E-02, 6.028E-02, 6.044E-02, 6.058E-02,
     &  6.069E-02, 6.078E-02, 6.085E-02, 6.089E-02, 6.090E-02,
     &  6.082E-02, 6.066E-02, 6.041E-02, 6.008E-02, 5.967E-02,
     &  5.917E-02, 5.859E-02, 5.713E-02, 5.525E-02, 5.219E-02,
     &  4.823E-02, 4.324E-02, 3.720E-02, 3.019E-02, 2.247E-02,
     &  1.455E-02, 7.275E-03, 1.865E-03/

      data ( ddc2byld(9,i), i = 1, 111 ) / ! Rh-102 beta minus
     &  2.688E-01, 2.688E-01, 2.688E-01, 2.688E-01, 2.688E-01,
     &  2.688E-01, 2.688E-01, 2.688E-01, 2.688E-01, 2.688E-01,
     &  2.688E-01, 2.688E-01, 2.688E-01, 2.688E-01, 2.688E-01,
     &  2.687E-01, 2.687E-01, 2.687E-01, 2.687E-01, 2.687E-01,
     &  2.687E-01, 2.687E-01, 2.687E-01, 2.687E-01, 2.687E-01,
     &  2.687E-01, 2.687E-01, 2.687E-01, 2.687E-01, 2.687E-01,
     &  2.687E-01, 2.687E-01, 2.687E-01, 2.687E-01, 2.687E-01,
     &  2.687E-01, 2.687E-01, 2.687E-01, 2.687E-01, 2.687E-01,
     &  2.687E-01, 2.687E-01, 2.686E-01, 2.686E-01, 2.686E-01,
     &  2.686E-01, 2.686E-01, 2.686E-01, 2.685E-01, 2.685E-01,
     &  2.685E-01, 2.687E-01, 2.689E-01, 2.690E-01, 2.692E-01,
     &  2.695E-01, 2.699E-01, 2.702E-01, 2.705E-01, 2.709E-01,
     &  2.712E-01, 2.715E-01, 2.721E-01, 2.728E-01, 2.734E-01,
     &  2.740E-01, 2.747E-01, 2.753E-01, 2.759E-01, 2.766E-01,
     &  2.778E-01, 2.791E-01, 2.807E-01, 2.823E-01, 2.838E-01,
     &  2.854E-01, 2.870E-01, 2.885E-01, 2.900E-01, 2.915E-01,
     &  2.930E-01, 2.944E-01, 2.972E-01, 2.998E-01, 3.022E-01,
     &  3.044E-01, 3.064E-01, 3.083E-01, 3.099E-01, 3.126E-01,
     &  3.144E-01, 3.153E-01, 3.154E-01, 3.147E-01, 3.131E-01,
     &  3.108E-01, 3.076E-01, 2.993E-01, 2.884E-01, 2.718E-01,
     &  2.528E-01, 2.325E-01, 2.123E-01, 1.910E-01, 1.678E-01,
     &  1.433E-01, 1.183E-01, 9.341E-02, 6.956E-02, 2.867E-02,
     &  3.691E-03/

      data ( ddc2byld(10,i), i = 1, 90 ) / ! Ag-106 beta minus
     &  1.314E-01, 1.313E-01, 1.313E-01, 1.312E-01, 1.312E-01,
     &  1.312E-01, 1.312E-01, 1.312E-01, 1.312E-01, 1.312E-01,
     &  1.311E-01, 1.311E-01, 1.311E-01, 1.311E-01, 1.310E-01,
     &  1.310E-01, 1.310E-01, 1.309E-01, 1.309E-01, 1.308E-01,
     &  1.308E-01, 1.307E-01, 1.307E-01, 1.306E-01, 1.306E-01,
     &  1.305E-01, 1.304E-01, 1.304E-01, 1.303E-01, 1.302E-01,
     &  1.301E-01, 1.300E-01, 1.299E-01, 1.297E-01, 1.296E-01,
     &  1.294E-01, 1.292E-01, 1.290E-01, 1.288E-01, 1.286E-01,
     &  1.283E-01, 1.281E-01, 1.279E-01, 1.275E-01, 1.270E-01,
     &  1.265E-01, 1.260E-01, 1.254E-01, 1.249E-01, 1.243E-01,
     &  1.238E-01, 1.233E-01, 1.228E-01, 1.224E-01, 1.219E-01,
     &  1.210E-01, 1.201E-01, 1.192E-01, 1.183E-01, 1.174E-01,
     &  1.165E-01, 1.156E-01, 1.138E-01, 1.120E-01, 1.101E-01,
     &  1.083E-01, 1.065E-01, 1.047E-01, 1.029E-01, 1.010E-01,
     &  9.740E-02, 9.378E-02, 8.928E-02, 8.480E-02, 8.036E-02,
     &  7.595E-02, 7.159E-02, 6.727E-02, 6.301E-02, 5.882E-02,
     &  5.469E-02, 5.063E-02, 4.278E-02, 3.534E-02, 2.838E-02,
     &  2.198E-02, 1.622E-02, 1.119E-02, 6.970E-03, 1.356E-03/

      data ( ddc2byld(11,i), i = 1, 116 ) / ! Ag-108 beta minus
     &  5.381E-01, 5.381E-01, 5.381E-01, 5.382E-01, 5.382E-01,
     &  5.382E-01, 5.382E-01, 5.382E-01, 5.382E-01, 5.382E-01,
     &  5.382E-01, 5.382E-01, 5.382E-01, 5.382E-01, 5.382E-01,
     &  5.382E-01, 5.382E-01, 5.383E-01, 5.383E-01, 5.383E-01,
     &  5.383E-01, 5.383E-01, 5.384E-01, 5.384E-01, 5.384E-01,
     &  5.384E-01, 5.384E-01, 5.385E-01, 5.385E-01, 5.385E-01,
     &  5.386E-01, 5.386E-01, 5.387E-01, 5.387E-01, 5.387E-01,
     &  5.388E-01, 5.389E-01, 5.390E-01, 5.391E-01, 5.391E-01,
     &  5.392E-01, 5.393E-01, 5.394E-01, 5.395E-01, 5.397E-01,
     &  5.399E-01, 5.401E-01, 5.403E-01, 5.405E-01, 5.406E-01,
     &  5.408E-01, 5.410E-01, 5.417E-01, 5.423E-01, 5.429E-01,
     &  5.440E-01, 5.452E-01, 5.463E-01, 5.474E-01, 5.486E-01,
     &  5.497E-01, 5.508E-01, 5.531E-01, 5.554E-01, 5.576E-01,
     &  5.599E-01, 5.621E-01, 5.644E-01, 5.666E-01, 5.689E-01,
     &  5.735E-01, 5.781E-01, 5.839E-01, 5.897E-01, 5.956E-01,
     &  6.015E-01, 6.075E-01, 6.134E-01, 6.194E-01, 6.254E-01,
     &  6.314E-01, 6.374E-01, 6.493E-01, 6.612E-01, 6.729E-01,
     &  6.846E-01, 6.961E-01, 7.074E-01, 7.185E-01, 7.402E-01,
     &  7.609E-01, 7.807E-01, 7.995E-01, 8.172E-01, 8.338E-01,
     &  8.492E-01, 8.634E-01, 8.882E-01, 9.079E-01, 9.252E-01,
     &  9.342E-01, 9.351E-01, 9.278E-01, 9.126E-01, 8.898E-01,
     &  8.598E-01, 8.232E-01, 7.805E-01, 7.325E-01, 6.238E-01,
     &  5.034E-01, 3.769E-01, 2.532E-01, 1.425E-01, 5.608E-02,
     &  6.637E-03/

      data ( ddc2byld(12,i), i = 1, 104 ) / ! In-112 beta minus
     &  1.159E+00, 1.159E+00, 1.159E+00, 1.159E+00, 1.159E+00,
     &  1.159E+00, 1.159E+00, 1.159E+00, 1.158E+00, 1.158E+00,
     &  1.158E+00, 1.158E+00, 1.158E+00, 1.158E+00, 1.158E+00,
     &  1.158E+00, 1.158E+00, 1.158E+00, 1.158E+00, 1.158E+00,
     &  1.158E+00, 1.158E+00, 1.158E+00, 1.158E+00, 1.158E+00,
     &  1.158E+00, 1.158E+00, 1.158E+00, 1.158E+00, 1.157E+00,
     &  1.157E+00, 1.157E+00, 1.157E+00, 1.157E+00, 1.157E+00,
     &  1.157E+00, 1.156E+00, 1.156E+00, 1.156E+00, 1.156E+00,
     &  1.155E+00, 1.155E+00, 1.155E+00, 1.154E+00, 1.154E+00,
     &  1.153E+00, 1.153E+00, 1.152E+00, 1.151E+00, 1.151E+00,
     &  1.150E+00, 1.150E+00, 1.149E+00, 1.149E+00, 1.150E+00,
     &  1.150E+00, 1.150E+00, 1.151E+00, 1.151E+00, 1.151E+00,
     &  1.151E+00, 1.152E+00, 1.152E+00, 1.152E+00, 1.153E+00,
     &  1.153E+00, 1.153E+00, 1.154E+00, 1.154E+00, 1.154E+00,
     &  1.154E+00, 1.154E+00, 1.154E+00, 1.154E+00, 1.154E+00,
     &  1.153E+00, 1.153E+00, 1.152E+00, 1.151E+00, 1.149E+00,
     &  1.148E+00, 1.146E+00, 1.141E+00, 1.135E+00, 1.129E+00,
     &  1.121E+00, 1.112E+00, 1.101E+00, 1.090E+00, 1.063E+00,
     &  1.032E+00, 9.971E-01, 9.577E-01, 9.144E-01, 8.676E-01,
     &  8.177E-01, 7.650E-01, 6.533E-01, 5.363E-01, 3.895E-01,
     &  2.512E-01, 1.323E-01, 4.486E-02, 2.076E-03/

      data ( ddc2byld(13,i), i = 1, 117 ) / ! In-114 beta minus
     &  3.799E-01, 3.799E-01, 3.799E-01, 3.799E-01, 3.799E-01,
     &  3.799E-01, 3.799E-01, 3.799E-01, 3.799E-01, 3.800E-01,
     &  3.800E-01, 3.800E-01, 3.800E-01, 3.800E-01, 3.800E-01,
     &  3.800E-01, 3.800E-01, 3.800E-01, 3.800E-01, 3.801E-01,
     &  3.801E-01, 3.801E-01, 3.801E-01, 3.801E-01, 3.802E-01,
     &  3.802E-01, 3.802E-01, 3.802E-01, 3.802E-01, 3.803E-01,
     &  3.803E-01, 3.804E-01, 3.804E-01, 3.804E-01, 3.805E-01,
     &  3.805E-01, 3.806E-01, 3.807E-01, 3.808E-01, 3.808E-01,
     &  3.809E-01, 3.810E-01, 3.810E-01, 3.812E-01, 3.813E-01,
     &  3.815E-01, 3.817E-01, 3.819E-01, 3.821E-01, 3.822E-01,
     &  3.824E-01, 3.826E-01, 3.828E-01, 3.833E-01, 3.838E-01,
     &  3.847E-01, 3.856E-01, 3.865E-01, 3.874E-01, 3.882E-01,
     &  3.891E-01, 3.900E-01, 3.918E-01, 3.935E-01, 3.953E-01,
     &  3.971E-01, 3.988E-01, 4.006E-01, 4.024E-01, 4.041E-01,
     &  4.077E-01, 4.113E-01, 4.158E-01, 4.203E-01, 4.249E-01,
     &  4.295E-01, 4.342E-01, 4.388E-01, 4.435E-01, 4.482E-01,
     &  4.529E-01, 4.577E-01, 4.671E-01, 4.766E-01, 4.860E-01,
     &  4.955E-01, 5.048E-01, 5.142E-01, 5.234E-01, 5.417E-01,
     &  5.595E-01, 5.770E-01, 5.939E-01, 6.103E-01, 6.262E-01,
     &  6.415E-01, 6.562E-01, 6.836E-01, 7.084E-01, 7.352E-01,
     &  7.574E-01, 7.748E-01, 7.873E-01, 7.948E-01, 7.973E-01,
     &  7.948E-01, 7.874E-01, 7.750E-01, 7.579E-01, 7.102E-01,
     &  6.465E-01, 5.694E-01, 4.824E-01, 3.893E-01, 2.947E-01,
     &  2.039E-01, 5.675E-02/

      data ( ddc2byld(14,i), i = 1, 117 ) / ! Sb-122 beta minus
     &  7.459E-01, 7.459E-01, 7.459E-01, 7.459E-01, 7.459E-01,
     &  7.459E-01, 7.459E-01, 7.459E-01, 7.459E-01, 7.459E-01,
     &  7.459E-01, 7.459E-01, 7.459E-01, 7.459E-01, 7.459E-01,
     &  7.459E-01, 7.459E-01, 7.460E-01, 7.460E-01, 7.460E-01,
     &  7.460E-01, 7.460E-01, 7.460E-01, 7.460E-01, 7.460E-01,
     &  7.460E-01, 7.460E-01, 7.460E-01, 7.461E-01, 7.461E-01,
     &  7.461E-01, 7.461E-01, 7.461E-01, 7.462E-01, 7.462E-01,
     &  7.462E-01, 7.462E-01, 7.463E-01, 7.463E-01, 7.464E-01,
     &  7.464E-01, 7.464E-01, 7.465E-01, 7.465E-01, 7.466E-01,
     &  7.467E-01, 7.468E-01, 7.469E-01, 7.469E-01, 7.470E-01,
     &  7.471E-01, 7.472E-01, 7.473E-01, 7.473E-01, 7.483E-01,
     &  7.495E-01, 7.507E-01, 7.519E-01, 7.531E-01, 7.543E-01,
     &  7.555E-01, 7.567E-01, 7.590E-01, 7.614E-01, 7.637E-01,
     &  7.660E-01, 7.683E-01, 7.706E-01, 7.729E-01, 7.752E-01,
     &  7.797E-01, 7.843E-01, 7.899E-01, 7.956E-01, 8.013E-01,
     &  8.069E-01, 8.125E-01, 8.181E-01, 8.237E-01, 8.292E-01,
     &  8.347E-01, 8.401E-01, 8.508E-01, 8.612E-01, 8.714E-01,
     &  8.813E-01, 8.908E-01, 9.000E-01, 9.088E-01, 9.253E-01,
     &  9.402E-01, 9.535E-01, 9.651E-01, 9.750E-01, 9.831E-01,
     &  9.895E-01, 9.942E-01, 9.985E-01, 9.961E-01, 9.845E-01,
     &  9.639E-01, 9.356E-01, 9.009E-01, 8.612E-01, 8.181E-01,
     &  7.728E-01, 7.225E-01, 6.676E-01, 6.092E-01, 4.860E-01,
     &  3.625E-01, 2.498E-01, 1.606E-01, 1.089E-01, 8.758E-02,
     &  6.481E-02, 2.013E-02/

      data ( ddc2byld(15,i), i = 1, 112 ) / ! I-126 beta minus
     &  9.144E-01, 9.143E-01, 9.143E-01, 9.143E-01, 9.143E-01,
     &  9.143E-01, 9.143E-01, 9.143E-01, 9.142E-01, 9.142E-01,
     &  9.142E-01, 9.142E-01, 9.142E-01, 9.141E-01, 9.141E-01,
     &  9.141E-01, 9.140E-01, 9.140E-01, 9.140E-01, 9.139E-01,
     &  9.139E-01, 9.138E-01, 9.138E-01, 9.137E-01, 9.137E-01,
     &  9.136E-01, 9.135E-01, 9.135E-01, 9.134E-01, 9.133E-01,
     &  9.132E-01, 9.131E-01, 9.130E-01, 9.129E-01, 9.128E-01,
     &  9.126E-01, 9.124E-01, 9.122E-01, 9.120E-01, 9.118E-01,
     &  9.115E-01, 9.113E-01, 9.111E-01, 9.107E-01, 9.103E-01,
     &  9.098E-01, 9.093E-01, 9.087E-01, 9.082E-01, 9.077E-01,
     &  9.072E-01, 9.066E-01, 9.061E-01, 9.056E-01, 9.050E-01,
     &  9.057E-01, 9.059E-01, 9.061E-01, 9.063E-01, 9.065E-01,
     &  9.067E-01, 9.068E-01, 9.071E-01, 9.074E-01, 9.076E-01,
     &  9.077E-01, 9.079E-01, 9.080E-01, 9.081E-01, 9.081E-01,
     &  9.082E-01, 9.081E-01, 9.079E-01, 9.076E-01, 9.072E-01,
     &  9.066E-01, 9.059E-01, 9.050E-01, 9.040E-01, 9.029E-01,
     &  9.016E-01, 9.002E-01, 8.969E-01, 8.930E-01, 8.884E-01,
     &  8.832E-01, 8.774E-01, 8.710E-01, 8.639E-01, 8.482E-01,
     &  8.302E-01, 8.105E-01, 7.891E-01, 7.664E-01, 7.429E-01,
     &  7.188E-01, 6.946E-01, 6.477E-01, 6.040E-01, 5.434E-01,
     &  4.766E-01, 4.057E-01, 3.337E-01, 2.634E-01, 1.982E-01,
     &  1.416E-01, 9.743E-02, 7.002E-02, 5.923E-02, 3.913E-02,
     &  1.866E-02, 3.188E-03/

      data ( ddc2byld(16,i), i = 1, 118 ) / ! I-128 beta minus
     &  3.576E-01, 3.576E-01, 3.576E-01, 3.576E-01, 3.576E-01,
     &  3.576E-01, 3.576E-01, 3.576E-01, 3.576E-01, 3.576E-01,
     &  3.576E-01, 3.577E-01, 3.577E-01, 3.577E-01, 3.577E-01,
     &  3.577E-01, 3.577E-01, 3.577E-01, 3.577E-01, 3.577E-01,
     &  3.578E-01, 3.578E-01, 3.578E-01, 3.578E-01, 3.578E-01,
     &  3.578E-01, 3.579E-01, 3.579E-01, 3.579E-01, 3.579E-01,
     &  3.580E-01, 3.580E-01, 3.580E-01, 3.581E-01, 3.581E-01,
     &  3.582E-01, 3.583E-01, 3.583E-01, 3.584E-01, 3.585E-01,
     &  3.585E-01, 3.586E-01, 3.587E-01, 3.588E-01, 3.589E-01,
     &  3.591E-01, 3.593E-01, 3.594E-01, 3.596E-01, 3.598E-01,
     &  3.599E-01, 3.601E-01, 3.603E-01, 3.604E-01, 3.606E-01,
     &  3.616E-01, 3.624E-01, 3.633E-01, 3.641E-01, 3.649E-01,
     &  3.657E-01, 3.665E-01, 3.682E-01, 3.698E-01, 3.714E-01,
     &  3.731E-01, 3.747E-01, 3.763E-01, 3.779E-01, 3.795E-01,
     &  3.828E-01, 3.860E-01, 3.901E-01, 3.942E-01, 3.983E-01,
     &  4.024E-01, 4.066E-01, 4.107E-01, 4.149E-01, 4.191E-01,
     &  4.232E-01, 4.274E-01, 4.358E-01, 4.442E-01, 4.525E-01,
     &  4.608E-01, 4.691E-01, 4.773E-01, 4.854E-01, 5.015E-01,
     &  5.171E-01, 5.324E-01, 5.473E-01, 5.616E-01, 5.755E-01,
     &  5.888E-01, 6.016E-01, 6.253E-01, 6.466E-01, 6.696E-01,
     &  6.883E-01, 7.028E-01, 7.128E-01, 7.183E-01, 7.195E-01,
     &  7.163E-01, 7.088E-01, 6.973E-01, 6.819E-01, 6.404E-01,
     &  5.867E-01, 5.234E-01, 4.521E-01, 3.759E-01, 2.985E-01,
     &  2.243E-01, 9.899E-02, 1.600E-02/

      data ( ddc2byld(17,i), i = 1, 98 ) / ! Cs-130 beta minus
     &  1.000E-01, 1.000E-01, 1.000E-01, 1.000E-01, 1.000E-01,
     &  9.999E-02, 9.999E-02, 9.999E-02, 9.998E-02, 9.997E-02,
     &  9.997E-02, 9.996E-02, 9.995E-02, 9.994E-02, 9.994E-02,
     &  9.993E-02, 9.992E-02, 9.990E-02, 9.988E-02, 9.987E-02,
     &  9.985E-02, 9.983E-02, 9.981E-02, 9.979E-02, 9.978E-02,
     &  9.976E-02, 9.974E-02, 9.972E-02, 9.969E-02, 9.965E-02,
     &  9.962E-02, 9.958E-02, 9.955E-02, 9.951E-02, 9.947E-02,
     &  9.940E-02, 9.933E-02, 9.926E-02, 9.919E-02, 9.912E-02,
     &  9.905E-02, 9.898E-02, 9.890E-02, 9.876E-02, 9.862E-02,
     &  9.844E-02, 9.826E-02, 9.808E-02, 9.790E-02, 9.772E-02,
     &  9.754E-02, 9.736E-02, 9.718E-02, 9.700E-02, 9.682E-02,
     &  9.658E-02, 9.635E-02, 9.612E-02, 9.588E-02, 9.564E-02,
     &  9.540E-02, 9.516E-02, 9.467E-02, 9.418E-02, 9.368E-02,
     &  9.317E-02, 9.266E-02, 9.214E-02, 9.161E-02, 9.108E-02,
     &  9.001E-02, 8.891E-02, 8.752E-02, 8.610E-02, 8.465E-02,
     &  8.318E-02, 8.168E-02, 8.016E-02, 7.862E-02, 7.705E-02,
     &  7.546E-02, 7.385E-02, 7.057E-02, 6.721E-02, 6.379E-02,
     &  6.031E-02, 5.679E-02, 5.324E-02, 4.967E-02, 4.252E-02,
     &  3.546E-02, 2.863E-02, 2.217E-02, 1.622E-02, 1.094E-02,
     &  6.497E-03, 3.070E-03, 1.993E-06/

      data ( ddc2byld(18,i), i = 1, 107 ) / ! Cs-132 beta minus
     &  6.550E-02, 6.548E-02, 6.548E-02, 6.547E-02, 6.547E-02,
     &  6.547E-02, 6.547E-02, 6.547E-02, 6.546E-02, 6.546E-02,
     &  6.545E-02, 6.545E-02, 6.544E-02, 6.544E-02, 6.543E-02,
     &  6.543E-02, 6.542E-02, 6.541E-02, 6.540E-02, 6.539E-02,
     &  6.537E-02, 6.536E-02, 6.535E-02, 6.534E-02, 6.533E-02,
     &  6.532E-02, 6.531E-02, 6.529E-02, 6.527E-02, 6.525E-02,
     &  6.522E-02, 6.520E-02, 6.518E-02, 6.516E-02, 6.513E-02,
     &  6.509E-02, 6.504E-02, 6.499E-02, 6.495E-02, 6.490E-02,
     &  6.485E-02, 6.481E-02, 6.476E-02, 6.467E-02, 6.458E-02,
     &  6.446E-02, 6.435E-02, 6.423E-02, 6.411E-02, 6.400E-02,
     &  6.388E-02, 6.377E-02, 6.365E-02, 6.353E-02, 6.342E-02,
     &  6.326E-02, 6.312E-02, 6.297E-02, 6.282E-02, 6.267E-02,
     &  6.251E-02, 6.236E-02, 6.205E-02, 6.174E-02, 6.143E-02,
     &  6.111E-02, 6.080E-02, 6.047E-02, 6.015E-02, 5.982E-02,
     &  5.917E-02, 5.850E-02, 5.766E-02, 5.681E-02, 5.596E-02,
     &  5.510E-02, 5.423E-02, 5.335E-02, 5.248E-02, 5.160E-02,
     &  5.072E-02, 4.983E-02, 4.806E-02, 4.630E-02, 4.456E-02,
     &  4.284E-02, 4.115E-02, 3.952E-02, 3.794E-02, 3.503E-02,
     &  3.252E-02, 3.054E-02, 2.920E-02, 2.849E-02, 2.776E-02,
     &  2.695E-02, 2.605E-02, 2.401E-02, 2.171E-02, 1.854E-02,
     &  1.517E-02, 1.174E-02, 8.418E-03, 5.383E-03, 2.832E-03,
     &  9.807E-04, 5.898E-05/

      data ( ddc2byld(19,i), i = 1, 110 ) / ! Eu-150m beta minus
     &  1.281E+00, 1.281E+00, 1.281E+00, 1.281E+00, 1.281E+00,
     &  1.281E+00, 1.281E+00, 1.281E+00, 1.281E+00, 1.281E+00,
     &  1.281E+00, 1.281E+00, 1.281E+00, 1.281E+00, 1.281E+00,
     &  1.281E+00, 1.281E+00, 1.281E+00, 1.281E+00, 1.281E+00,
     &  1.281E+00, 1.281E+00, 1.281E+00, 1.281E+00, 1.281E+00,
     &  1.281E+00, 1.281E+00, 1.281E+00, 1.281E+00, 1.281E+00,
     &  1.281E+00, 1.281E+00, 1.281E+00, 1.281E+00, 1.281E+00,
     &  1.281E+00, 1.281E+00, 1.281E+00, 1.281E+00, 1.281E+00,
     &  1.281E+00, 1.281E+00, 1.281E+00, 1.281E+00, 1.281E+00,
     &  1.281E+00, 1.281E+00, 1.281E+00, 1.281E+00, 1.281E+00,
     &  1.281E+00, 1.281E+00, 1.281E+00, 1.281E+00, 1.281E+00,
     &  1.281E+00, 1.281E+00, 1.283E+00, 1.284E+00, 1.286E+00,
     &  1.287E+00, 1.289E+00, 1.292E+00, 1.294E+00, 1.297E+00,
     &  1.300E+00, 1.303E+00, 1.306E+00, 1.309E+00, 1.311E+00,
     &  1.317E+00, 1.322E+00, 1.328E+00, 1.334E+00, 1.340E+00,
     &  1.346E+00, 1.352E+00, 1.357E+00, 1.362E+00, 1.367E+00,
     &  1.372E+00, 1.377E+00, 1.386E+00, 1.394E+00, 1.401E+00,
     &  1.407E+00, 1.413E+00, 1.418E+00, 1.422E+00, 1.427E+00,
     &  1.429E+00, 1.428E+00, 1.423E+00, 1.414E+00, 1.401E+00,
     &  1.386E+00, 1.366E+00, 1.318E+00, 1.257E+00, 1.166E+00,
     &  1.058E+00, 9.387E-01, 8.100E-01, 6.761E-01, 5.414E-01,
     &  4.104E-01, 2.882E-01, 1.803E-01, 9.232E-02, 1.396E-03/

      data ( ddc2byld(20,i), i = 1, 114 ) / ! Eu-152 beta minus
     &  8.829E-01, 8.826E-01, 8.826E-01, 8.825E-01, 8.825E-01,
     &  8.824E-01, 8.824E-01, 8.824E-01, 8.823E-01, 8.822E-01,
     &  8.822E-01, 8.821E-01, 8.820E-01, 8.820E-01, 8.819E-01,
     &  8.818E-01, 8.817E-01, 8.815E-01, 8.814E-01, 8.812E-01,
     &  8.810E-01, 8.808E-01, 8.807E-01, 8.805E-01, 8.803E-01,
     &  8.801E-01, 8.800E-01, 8.798E-01, 8.794E-01, 8.791E-01,
     &  8.787E-01, 8.784E-01, 8.780E-01, 8.777E-01, 8.773E-01,
     &  8.766E-01, 8.759E-01, 8.752E-01, 8.745E-01, 8.738E-01,
     &  8.731E-01, 8.724E-01, 8.717E-01, 8.703E-01, 8.689E-01,
     &  8.672E-01, 8.654E-01, 8.637E-01, 8.620E-01, 8.602E-01,
     &  8.585E-01, 8.567E-01, 8.550E-01, 8.532E-01, 8.515E-01,
     &  8.480E-01, 8.445E-01, 8.423E-01, 8.398E-01, 8.374E-01,
     &  8.350E-01, 8.325E-01, 8.276E-01, 8.226E-01, 8.177E-01,
     &  8.127E-01, 8.077E-01, 8.027E-01, 7.977E-01, 7.926E-01,
     &  7.825E-01, 7.723E-01, 7.596E-01, 7.469E-01, 7.343E-01,
     &  7.216E-01, 7.091E-01, 6.967E-01, 6.843E-01, 6.722E-01,
     &  6.602E-01, 6.484E-01, 6.254E-01, 6.035E-01, 5.828E-01,
     &  5.636E-01, 5.461E-01, 5.303E-01, 5.165E-01, 4.958E-01,
     &  4.775E-01, 4.584E-01, 4.385E-01, 4.181E-01, 3.974E-01,
     &  3.765E-01, 3.559E-01, 3.163E-01, 2.807E-01, 2.366E-01,
     &  1.932E-01, 1.531E-01, 1.191E-01, 9.427E-02, 8.207E-02,
     &  7.582E-02, 6.918E-02, 6.232E-02, 5.539E-02, 4.160E-02,
     &  2.878E-02, 1.717E-02, 7.656E-03, 1.532E-03/

      data ( ddc2byld(21,i), i = 1, 117 ) / ! Eu-152m beta minus
     &  4.215E-01, 4.215E-01, 4.215E-01, 4.215E-01, 4.215E-01,
     &  4.215E-01, 4.215E-01, 4.215E-01, 4.215E-01, 4.215E-01,
     &  4.215E-01, 4.215E-01, 4.215E-01, 4.215E-01, 4.215E-01,
     &  4.215E-01, 4.215E-01, 4.215E-01, 4.215E-01, 4.215E-01,
     &  4.215E-01, 4.215E-01, 4.215E-01, 4.215E-01, 4.215E-01,
     &  4.215E-01, 4.215E-01, 4.215E-01, 4.215E-01, 4.215E-01,
     &  4.215E-01, 4.215E-01, 4.215E-01, 4.215E-01, 4.215E-01,
     &  4.215E-01, 4.215E-01, 4.215E-01, 4.215E-01, 4.215E-01,
     &  4.215E-01, 4.215E-01, 4.215E-01, 4.215E-01, 4.215E-01,
     &  4.215E-01, 4.215E-01, 4.214E-01, 4.214E-01, 4.214E-01,
     &  4.214E-01, 4.214E-01, 4.214E-01, 4.214E-01, 4.214E-01,
     &  4.214E-01, 4.214E-01, 4.220E-01, 4.226E-01, 4.231E-01,
     &  4.236E-01, 4.241E-01, 4.252E-01, 4.262E-01, 4.272E-01,
     &  4.283E-01, 4.293E-01, 4.303E-01, 4.313E-01, 4.323E-01,
     &  4.344E-01, 4.364E-01, 4.389E-01, 4.415E-01, 4.440E-01,
     &  4.466E-01, 4.492E-01, 4.518E-01, 4.545E-01, 4.571E-01,
     &  4.598E-01, 4.626E-01, 4.681E-01, 4.738E-01, 4.795E-01,
     &  4.851E-01, 4.906E-01, 4.959E-01, 5.010E-01, 5.110E-01,
     &  5.203E-01, 5.291E-01, 5.372E-01, 5.448E-01, 5.517E-01,
     &  5.581E-01, 5.639E-01, 5.739E-01, 5.819E-01, 5.894E-01,
     &  5.948E-01, 5.986E-01, 6.000E-01, 5.973E-01, 5.908E-01,
     &  5.805E-01, 5.665E-01, 5.491E-01, 5.283E-01, 4.777E-01,
     &  4.169E-01, 3.487E-01, 2.765E-01, 2.043E-01, 1.366E-01,
     &  7.818E-02, 5.418E-03/

      data ( ddc2byld(22,i), i = 1, 89 ) / ! Tm-168 beta minus
     &  1.574E-03, 1.573E-03, 1.573E-03, 1.572E-03, 1.572E-03,
     &  1.572E-03, 1.572E-03, 1.572E-03, 1.572E-03, 1.571E-03,
     &  1.571E-03, 1.571E-03, 1.570E-03, 1.570E-03, 1.570E-03,
     &  1.569E-03, 1.569E-03, 1.568E-03, 1.567E-03, 1.567E-03,
     &  1.566E-03, 1.565E-03, 1.564E-03, 1.564E-03, 1.563E-03,
     &  1.562E-03, 1.561E-03, 1.560E-03, 1.559E-03, 1.557E-03,
     &  1.556E-03, 1.554E-03, 1.553E-03, 1.551E-03, 1.550E-03,
     &  1.547E-03, 1.543E-03, 1.540E-03, 1.537E-03, 1.534E-03,
     &  1.531E-03, 1.528E-03, 1.525E-03, 1.519E-03, 1.513E-03,
     &  1.505E-03, 1.498E-03, 1.490E-03, 1.482E-03, 1.475E-03,
     &  1.467E-03, 1.459E-03, 1.452E-03, 1.444E-03, 1.437E-03,
     &  1.422E-03, 1.407E-03, 1.392E-03, 1.377E-03, 1.365E-03,
     &  1.351E-03, 1.338E-03, 1.311E-03, 1.285E-03, 1.258E-03,
     &  1.232E-03, 1.205E-03, 1.179E-03, 1.153E-03, 1.127E-03,
     &  1.075E-03, 1.023E-03, 9.598E-04, 8.971E-04, 8.354E-04,
     &  7.750E-04, 7.158E-04, 6.580E-04, 6.018E-04, 5.472E-04,
     &  4.944E-04, 4.435E-04, 3.480E-04, 2.617E-04, 1.857E-04,
     &  1.211E-04, 6.904E-05, 3.069E-05, 7.284E-06/

      data ( ddc2byld(23,i), i = 1, 85 ) / ! Cl-36 beta plus
     &  0.000E+00, 5.363E-19, 3.012E-18, 1.357E-17, 5.121E-17,
     &  1.672E-16, 4.839E-16, 1.266E-15, 6.768E-15, 2.796E-14,
     &  9.475E-14, 2.748E-13, 7.031E-13, 1.624E-12, 3.448E-12,
     &  6.818E-12, 2.247E-11, 6.200E-11, 1.835E-10, 4.661E-10,
     &  1.056E-09, 2.185E-09, 4.202E-09, 7.600E-09, 1.304E-08,
     &  2.137E-08, 3.363E-08, 5.104E-08, 1.070E-07, 2.016E-07,
     &  3.484E-07, 5.607E-07, 8.504E-07, 1.228E-06, 1.701E-06,
     &  2.960E-06, 4.653E-06, 6.787E-06, 9.352E-06, 1.233E-05,
     &  1.570E-05, 1.944E-05, 2.353E-05, 3.267E-05, 4.295E-05,
     &  5.724E-05, 7.293E-05, 8.986E-05, 1.079E-04, 1.269E-04,
     &  1.469E-04, 1.676E-04, 1.891E-04, 2.113E-04, 2.341E-04,
     &  2.814E-04, 3.306E-04, 3.814E-04, 4.336E-04, 4.868E-04,
     &  5.410E-04, 5.958E-04, 7.069E-04, 8.187E-04, 9.302E-04,
     &  1.040E-03, 1.148E-03, 1.253E-03, 1.355E-03, 1.451E-03,
     &  1.630E-03, 1.784E-03, 1.938E-03, 2.044E-03, 2.099E-03,
     &  2.101E-03, 2.050E-03, 1.949E-03, 1.801E-03, 1.611E-03,
     &  1.388E-03, 1.141E-03, 6.264E-04, 1.910E-04, 1.106E-08/

      data ( ddc2byld(24,i), i = 1, 100 ) / ! K-40 beta plus
     &  0.000E+00, 1.762E-21, 1.227E-20, 6.674E-20, 2.973E-19,
     &  1.125E-18, 3.717E-18, 1.096E-17, 7.228E-17, 3.564E-16,
     &  1.406E-15, 4.657E-15, 1.340E-14, 3.436E-14, 8.014E-14,
     &  1.726E-13, 6.608E-13, 2.075E-12, 7.063E-12, 2.029E-11,
     &  5.132E-11, 1.174E-10, 2.473E-10, 4.861E-10, 9.005E-10,
     &  1.583E-09, 2.655E-09, 4.273E-09, 9.922E-09, 2.038E-08,
     &  3.786E-08, 6.475E-08, 1.034E-07, 1.562E-07, 2.249E-07,
     &  4.164E-07, 6.872E-07, 1.042E-06, 1.481E-06, 2.004E-06,
     &  2.607E-06, 3.288E-06, 4.044E-06, 5.761E-06, 7.733E-06,
     &  1.052E-05, 1.362E-05, 1.701E-05, 2.065E-05, 2.453E-05,
     &  2.863E-05, 3.293E-05, 3.741E-05, 4.206E-05, 4.687E-05,
     &  5.693E-05, 6.752E-05, 7.858E-05, 9.006E-05, 1.019E-04,
     &  1.141E-04, 1.266E-04, 1.525E-04, 1.793E-04, 2.069E-04,
     &  2.353E-04, 2.643E-04, 2.939E-04, 3.239E-04, 3.544E-04,
     &  4.165E-04, 4.799E-04, 5.606E-04, 6.424E-04, 7.252E-04,
     &  8.085E-04, 8.922E-04, 9.761E-04, 1.060E-03, 1.143E-03,
     &  1.227E-03, 1.309E-03, 1.472E-03, 1.632E-03, 1.786E-03,
     &  1.936E-03, 2.079E-03, 2.215E-03, 2.344E-03, 2.576E-03,
     &  2.772E-03, 2.927E-03, 3.038E-03, 3.099E-03, 3.107E-03,
     &  3.056E-03, 2.939E-03, 2.479E-03, 1.694E-03, 4.504E-04/

      data ( ddc2byld(25,i), i = 1, 97 ) / ! Mn-54 beta plus
     &  0.000E+00, 4.815E-32, 6.406E-31, 6.124E-30, 4.488E-29,
     &  2.646E-28, 1.302E-27, 5.508E-27, 6.807E-26, 5.712E-25,
     &  3.562E-24, 1.758E-23, 7.199E-23, 2.530E-22, 7.840E-22,
     &  2.187E-21, 1.323E-20, 6.180E-20, 3.257E-19, 1.380E-18,
     &  4.962E-18, 1.564E-17, 4.421E-17, 1.137E-16, 2.692E-16,
     &  5.919E-16, 1.217E-15, 2.358E-15, 7.552E-15, 2.025E-14,
     &  4.700E-14, 9.690E-14, 1.814E-13, 3.135E-13, 5.073E-13,
     &  1.136E-12, 2.175E-12, 3.712E-12, 5.815E-12, 8.530E-12,
     &  1.189E-11, 1.590E-11, 2.059E-11, 3.195E-11, 4.592E-11,
     &  6.695E-11, 9.175E-11, 1.201E-10, 1.519E-10, 1.870E-10,
     &  2.251E-10, 2.662E-10, 3.101E-10, 3.567E-10, 4.059E-10,
     &  5.114E-10, 6.260E-10, 7.487E-10, 8.791E-10, 1.016E-09,
     &  1.160E-09, 1.310E-09, 1.625E-09, 1.958E-09, 2.307E-09,
     &  2.670E-09, 3.045E-09, 3.429E-09, 3.822E-09, 4.222E-09,
     &  5.040E-09, 5.876E-09, 6.936E-09, 8.005E-09, 9.077E-09,
     &  1.014E-08, 1.120E-08, 1.225E-08, 1.328E-08, 1.429E-08,
     &  1.527E-08, 1.623E-08, 1.807E-08, 1.978E-08, 2.134E-08,
     &  2.276E-08, 2.400E-08, 2.507E-08, 2.593E-08, 2.702E-08,
     &  2.714E-08, 2.616E-08, 2.399E-08, 2.059E-08, 1.607E-08,
     &  1.079E-08, 5.437E-09/

      data ( ddc2byld(26,i), i = 1, 104 ) / ! Cu-64 beta plus
     &  0.000E+00, 0.000E+00, 3.596E-27, 5.007E-26, 5.113E-25,
     &  4.052E-24, 2.601E-23, 1.398E-22, 2.628E-21, 3.143E-20,
     &  2.659E-19, 1.714E-18, 8.881E-18, 3.856E-17, 1.447E-16,
     &  4.808E-16, 3.976E-15, 2.448E-14, 1.758E-13, 9.852E-13,
     &  4.570E-12, 1.819E-11, 6.354E-11, 1.981E-10, 5.578E-10,
     &  1.434E-09, 3.392E-09, 7.449E-09, 2.956E-08, 9.426E-08,
     &  2.519E-07, 5.835E-07, 1.203E-06, 2.257E-06, 3.916E-06,
     &  9.813E-06, 2.046E-05, 3.736E-05, 6.179E-05, 9.482E-05,
     &  1.372E-04, 1.896E-04, 2.523E-04, 4.099E-04, 6.109E-04,
     &  9.234E-04, 1.303E-03, 1.747E-03, 2.255E-03, 2.824E-03,
     &  3.452E-03, 4.137E-03, 4.877E-03, 5.670E-03, 6.514E-03,
     &  8.347E-03, 1.036E-02, 1.254E-02, 1.487E-02, 1.734E-02,
     &  1.995E-02, 2.266E-02, 2.841E-02, 3.451E-02, 4.090E-02,
     &  4.753E-02, 5.436E-02, 6.134E-02, 6.843E-02, 7.562E-02,
     &  9.017E-02, 1.048E-01, 1.232E-01, 1.413E-01, 1.592E-01,
     &  1.767E-01, 1.937E-01, 2.103E-01, 2.264E-01, 2.420E-01,
     &  2.571E-01, 2.716E-01, 2.990E-01, 3.243E-01, 3.474E-01,
     &  3.685E-01, 3.875E-01, 4.044E-01, 4.194E-01, 4.438E-01,
     &  4.609E-01, 4.712E-01, 4.752E-01, 4.731E-01, 4.656E-01,
     &  4.530E-01, 4.357E-01, 3.894E-01, 3.306E-01, 2.467E-01,
     &  1.603E-01, 8.225E-02, 2.446E-02, 9.276E-05/

      data ( ddc2byld(27,i), i = 1, 115 ) / ! As-74 beta plus
     &  0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     &  1.475E-27, 1.235E-26, 8.441E-26, 2.412E-24, 4.113E-23,
     &  4.724E-22, 3.978E-21, 2.615E-20, 1.406E-19, 6.408E-19,
     &  2.549E-18, 2.918E-17, 2.406E-16, 2.408E-15, 1.829E-14,
     &  1.119E-13, 5.731E-13, 2.513E-12, 9.601E-12, 3.240E-11,
     &  9.777E-11, 2.667E-10, 6.642E-10, 3.254E-09, 1.227E-08,
     &  3.751E-08, 9.701E-08, 2.192E-07, 4.438E-07, 8.218E-07,
     &  2.287E-06, 5.168E-06, 1.006E-05, 1.754E-05, 2.811E-05,
     &  4.222E-05, 6.023E-05, 8.245E-05, 1.404E-04, 2.174E-04,
     &  3.416E-04, 4.978E-04, 6.863E-04, 9.071E-04, 1.160E-03,
     &  1.445E-03, 1.762E-03, 2.110E-03, 2.488E-03, 2.897E-03,
     &  3.802E-03, 4.821E-03, 5.947E-03, 7.176E-03, 8.502E-03,
     &  9.921E-03, 1.143E-02, 1.467E-02, 1.821E-02, 2.200E-02,
     &  2.601E-02, 3.022E-02, 3.459E-02, 3.911E-02, 4.375E-02,
     &  5.335E-02, 6.325E-02, 7.594E-02, 8.883E-02, 1.018E-01,
     &  1.148E-01, 1.278E-01, 1.407E-01, 1.535E-01, 1.661E-01,
     &  1.785E-01, 1.907E-01, 2.146E-01, 2.376E-01, 2.596E-01,
     &  2.807E-01, 3.009E-01, 3.200E-01, 3.382E-01, 3.717E-01,
     &  4.014E-01, 4.274E-01, 4.497E-01, 4.686E-01, 4.840E-01,
     &  4.961E-01, 5.050E-01, 5.135E-01, 5.104E-01, 4.916E-01,
     &  4.582E-01, 4.128E-01, 3.581E-01, 2.973E-01, 2.339E-01,
     &  1.718E-01, 1.152E-01, 6.881E-02, 3.751E-02, 2.491E-02,
     &  2.072E-02, 1.548E-02, 9.623E-03, 4.074E-03, 4.165E-04/

      data ( ddc2byld(28,i), i = 1, 120 ) / ! Br-78 beta plus
     &  0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     &  0.000E+00, 0.000E+00, 8.918E-28, 3.142E-26, 6.400E-25,
     &  8.566E-24, 8.252E-23, 6.114E-22, 3.663E-21, 1.844E-20,
     &  8.046E-20, 1.091E-18, 1.050E-17, 1.255E-16, 1.122E-15,
     &  7.964E-15, 4.659E-14, 2.300E-13, 9.759E-13, 3.611E-12,
     &  1.181E-11, 3.454E-11, 9.149E-11, 4.960E-10, 2.024E-09,
     &  6.588E-09, 1.792E-08, 4.222E-08, 8.853E-08, 1.688E-07,
     &  4.926E-07, 1.155E-06, 2.314E-06, 4.131E-06, 6.756E-06,
     &  1.032E-05, 1.494E-05, 2.072E-05, 3.606E-05, 5.684E-05,
     &  9.100E-05, 1.347E-04, 1.881E-04, 2.516E-04, 3.251E-04,
     &  4.087E-04, 5.025E-04, 6.064E-04, 7.204E-04, 8.443E-04,
     &  1.122E-03, 1.438E-03, 1.792E-03, 2.182E-03, 2.607E-03,
     &  3.066E-03, 3.557E-03, 4.630E-03, 5.816E-03, 7.104E-03,
     &  8.484E-03, 9.949E-03, 1.149E-02, 1.310E-02, 1.476E-02,
     &  1.826E-02, 2.194E-02, 2.675E-02, 3.173E-02, 3.684E-02,
     &  4.207E-02, 4.738E-02, 5.274E-02, 5.816E-02, 6.361E-02,
     &  6.908E-02, 7.457E-02, 8.557E-02, 9.657E-02, 1.075E-01,
     &  1.185E-01, 1.293E-01, 1.401E-01, 1.508E-01, 1.720E-01,
     &  1.929E-01, 2.134E-01, 2.335E-01, 2.532E-01, 2.726E-01,
     &  2.915E-01, 3.101E-01, 3.459E-01, 3.800E-01, 4.201E-01,
     &  4.573E-01, 4.914E-01, 5.223E-01, 5.499E-01, 5.741E-01,
     &  5.949E-01, 6.121E-01, 6.257E-01, 6.358E-01, 6.453E-01,
     &  6.410E-01, 6.233E-01, 5.935E-01, 5.527E-01, 5.028E-01,
     &  4.459E-01, 3.213E-01, 2.019E-01, 9.593E-02, 2.066E-02/

      data ( ddc2byld(29,i), i = 1, 107 ) / ! Br-80 beta plus
     &  0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     &  0.000E+00, 0.000E+00, 2.631E-28, 9.269E-27, 1.888E-25,
     &  2.527E-24, 2.434E-23, 1.803E-22, 1.080E-21, 5.439E-21,
     &  2.373E-20, 3.217E-19, 3.096E-18, 3.702E-17, 3.307E-16,
     &  2.348E-15, 1.373E-14, 6.780E-14, 2.876E-13, 1.064E-12,
     &  3.479E-12, 1.018E-11, 2.695E-11, 1.461E-10, 5.960E-10,
     &  1.940E-09, 5.277E-09, 1.243E-08, 2.605E-08, 4.969E-08,
     &  1.449E-07, 3.396E-07, 6.804E-07, 1.214E-06, 1.985E-06,
     &  3.031E-06, 4.387E-06, 6.081E-06, 1.057E-05, 1.666E-05,
     &  2.665E-05, 3.940E-05, 5.499E-05, 7.347E-05, 9.485E-05,
     &  1.192E-04, 1.464E-04, 1.765E-04, 2.095E-04, 2.453E-04,
     &  3.253E-04, 4.163E-04, 5.178E-04, 6.295E-04, 7.508E-04,
     &  8.814E-04, 1.021E-03, 1.324E-03, 1.657E-03, 2.017E-03,
     &  2.401E-03, 2.805E-03, 3.228E-03, 3.667E-03, 4.119E-03,
     &  5.059E-03, 6.035E-03, 7.290E-03, 8.569E-03, 9.861E-03,
     &  1.116E-02, 1.245E-02, 1.373E-02, 1.500E-02, 1.625E-02,
     &  1.748E-02, 1.869E-02, 2.104E-02, 2.329E-02, 2.544E-02,
     &  2.747E-02, 2.940E-02, 3.121E-02, 3.292E-02, 3.603E-02,
     &  3.872E-02, 4.098E-02, 4.283E-02, 4.428E-02, 4.533E-02,
     &  4.601E-02, 4.632E-02, 4.593E-02, 4.428E-02, 4.068E-02,
     &  3.570E-02, 2.972E-02, 2.315E-02, 1.648E-02, 1.023E-02,
     &  4.969E-03, 1.335E-03/

      data ( ddc2byld(30,i), i = 1, 116 ) / ! Rb-84 beta plus
     &  0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     &  0.000E+00, 0.000E+00, 0.000E+00, 3.054E-27, 7.430E-26,
     &  1.159E-24, 1.278E-23, 1.069E-22, 7.145E-22, 3.982E-21,
     &  1.910E-20, 3.086E-19, 3.492E-18, 5.034E-17, 5.338E-16,
     &  4.424E-15, 2.971E-14, 1.657E-13, 7.812E-13, 3.168E-12,
     &  1.121E-11, 3.514E-11, 9.876E-11, 5.898E-10, 2.592E-09,
     &  8.950E-09, 2.552E-08, 6.246E-08, 1.352E-07, 2.650E-07,
     &  8.071E-07, 1.956E-06, 4.024E-06, 7.339E-06, 1.222E-05,
     &  1.895E-05, 2.780E-05, 3.900E-05, 6.919E-05, 1.108E-04,
     &  1.803E-04, 2.703E-04, 3.818E-04, 5.155E-04, 6.718E-04,
     &  8.509E-04, 1.053E-03, 1.279E-03, 1.527E-03, 1.799E-03,
     &  2.411E-03, 3.115E-03, 3.907E-03, 4.786E-03, 5.748E-03,
     &  6.791E-03, 7.911E-03, 1.037E-02, 1.310E-02, 1.608E-02,
     &  1.928E-02, 2.268E-02, 2.626E-02, 3.000E-02, 3.387E-02,
     &  4.198E-02, 5.047E-02, 6.148E-02, 7.279E-02, 8.429E-02,
     &  9.589E-02, 1.075E-01, 1.191E-01, 1.306E-01, 1.420E-01,
     &  1.533E-01, 1.644E-01, 1.859E-01, 2.067E-01, 2.264E-01,
     &  2.453E-01, 2.631E-01, 2.798E-01, 2.956E-01, 3.240E-01,
     &  3.483E-01, 3.687E-01, 3.853E-01, 3.981E-01, 4.074E-01,
     &  4.133E-01, 4.160E-01, 4.123E-01, 3.981E-01, 3.679E-01,
     &  3.275E-01, 2.808E-01, 2.323E-01, 1.870E-01, 1.502E-01,
     &  1.278E-01, 1.232E-01, 1.213E-01, 1.185E-01, 1.102E-01,
     &  9.779E-02, 8.131E-02, 6.126E-02, 3.918E-02, 1.807E-02,
     &  3.029E-03/

      data ( ddc2byld(31,i), i = 1, 113 ) / ! Rh-102 beta plus
     &  0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     &  0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     &  0.000E+00, 0.000E+00, 0.000E+00, 9.166E-27, 7.965E-26,
     &  5.842E-25, 2.112E-23, 5.088E-22, 1.755E-20, 4.067E-19,
     &  6.674E-18, 8.043E-17, 7.348E-16, 5.242E-15, 3.004E-14,
     &  1.420E-13, 5.674E-13, 1.958E-12, 1.617E-11, 9.045E-11,
     &  3.757E-10, 1.240E-09, 3.417E-09, 8.161E-09, 1.737E-08,
     &  6.043E-08, 1.623E-07, 3.627E-07, 7.085E-07, 1.250E-06,
     &  2.039E-06, 3.126E-06, 4.560E-06, 8.648E-06, 1.464E-05,
     &  2.524E-05, 3.974E-05, 5.856E-05, 8.205E-05, 1.105E-04,
     &  1.442E-04, 1.834E-04, 2.283E-04, 2.790E-04, 3.358E-04,
     &  4.680E-04, 6.258E-04, 8.097E-04, 1.020E-03, 1.258E-03,
     &  1.522E-03, 1.812E-03, 2.473E-03, 3.237E-03, 4.100E-03,
     &  5.057E-03, 6.103E-03, 7.234E-03, 8.443E-03, 9.726E-03,
     &  1.249E-02, 1.549E-02, 1.952E-02, 2.380E-02, 2.828E-02,
     &  3.292E-02, 3.768E-02, 4.253E-02, 4.744E-02, 5.238E-02,
     &  5.735E-02, 6.232E-02, 7.222E-02, 8.200E-02, 9.159E-02,
     &  1.009E-01, 1.100E-01, 1.188E-01, 1.273E-01, 1.434E-01,
     &  1.580E-01, 1.714E-01, 1.833E-01, 1.938E-01, 2.029E-01,
     &  2.106E-01, 2.171E-01, 2.261E-01, 2.303E-01, 2.293E-01,
     &  2.220E-01, 2.096E-01, 1.931E-01, 1.740E-01, 1.537E-01,
     &  1.339E-01, 1.162E-01, 1.019E-01, 8.746E-02, 5.741E-02,
     &  2.945E-02, 8.462E-03, 9.363E-07/

      data ( ddc2byld(32,i), i = 1, 117 ) / ! Ag-106 beta plus
     &  0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     &  0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     &  0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 7.221E-27,
     &  6.001E-26, 2.760E-24, 8.331E-23, 3.707E-21, 1.068E-19,
     &  2.100E-18, 2.931E-17, 3.011E-16, 2.358E-15, 1.456E-14,
     &  7.305E-14, 3.063E-13, 1.099E-12, 9.633E-12, 5.623E-11,
     &  2.411E-10, 8.154E-10, 2.292E-09, 5.563E-09, 1.201E-08,
     &  4.272E-08, 1.168E-07, 2.649E-07, 5.240E-07, 9.347E-07,
     &  1.540E-06, 2.381E-06, 3.500E-06, 6.731E-06, 1.153E-05,
     &  2.013E-05, 3.205E-05, 4.768E-05, 6.738E-05, 9.147E-05,
     &  1.202E-04, 1.539E-04, 1.927E-04, 2.369E-04, 2.867E-04,
     &  4.035E-04, 5.444E-04, 7.102E-04, 9.016E-04, 1.119E-03,
     &  1.363E-03, 1.634E-03, 2.254E-03, 2.981E-03, 3.811E-03,
     &  4.742E-03, 5.770E-03, 6.891E-03, 8.101E-03, 9.395E-03,
     &  1.222E-02, 1.533E-02, 1.957E-02, 2.416E-02, 2.904E-02,
     &  3.416E-02, 3.950E-02, 4.501E-02, 5.068E-02, 5.647E-02,
     &  6.237E-02, 6.835E-02, 8.052E-02, 9.288E-02, 1.053E-01,
     &  1.179E-01, 1.304E-01, 1.428E-01, 1.552E-01, 1.798E-01,
     &  2.038E-01, 2.272E-01, 2.500E-01, 2.720E-01, 2.933E-01,
     &  3.138E-01, 3.335E-01, 3.704E-01, 4.039E-01, 4.408E-01,
     &  4.721E-01, 4.978E-01, 5.178E-01, 5.321E-01, 5.410E-01,
     &  5.445E-01, 5.428E-01, 5.361E-01, 5.247E-01, 4.890E-01,
     &  4.389E-01, 3.781E-01, 3.114E-01, 2.441E-01, 1.804E-01,
     &  1.196E-01, 2.506E-02/

      data ( ddc2byld(33,i), i = 1, 108 ) / ! Ag-108 beta plus
     &  0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     &  0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     &  0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 2.101E-28,
     &  1.745E-27, 8.027E-26, 2.423E-24, 1.078E-22, 3.107E-21,
     &  6.108E-20, 8.524E-19, 8.755E-18, 6.857E-17, 4.233E-16,
     &  2.124E-15, 8.904E-15, 3.196E-14, 2.800E-13, 1.634E-12,
     &  7.005E-12, 2.369E-11, 6.657E-11, 1.616E-10, 3.487E-10,
     &  1.240E-09, 3.391E-09, 7.689E-09, 1.520E-08, 2.711E-08,
     &  4.465E-08, 6.903E-08, 1.015E-07, 1.950E-07, 3.338E-07,
     &  5.825E-07, 9.269E-07, 1.378E-06, 1.946E-06, 2.641E-06,
     &  3.469E-06, 4.438E-06, 5.554E-06, 6.824E-06, 8.252E-06,
     &  1.160E-05, 1.563E-05, 2.036E-05, 2.582E-05, 3.201E-05,
     &  3.894E-05, 4.661E-05, 6.417E-05, 8.464E-05, 1.079E-04,
     &  1.340E-04, 1.626E-04, 1.937E-04, 2.271E-04, 2.627E-04,
     &  3.400E-04, 4.243E-04, 5.383E-04, 6.600E-04, 7.881E-04,
     &  9.210E-04, 1.058E-03, 1.197E-03, 1.339E-03, 1.481E-03,
     &  1.625E-03, 1.768E-03, 2.053E-03, 2.334E-03, 2.609E-03,
     &  2.875E-03, 3.132E-03, 3.378E-03, 3.614E-03, 4.049E-03,
     &  4.436E-03, 4.774E-03, 5.062E-03, 5.300E-03, 5.490E-03,
     &  5.632E-03, 5.729E-03, 5.794E-03, 5.700E-03, 5.391E-03,
     &  4.903E-03, 4.264E-03, 3.517E-03, 2.713E-03, 1.910E-03,
     &  1.168E-03, 5.531E-04, 1.384E-04/

      data ( ddc2byld(34,i), i = 1, 115 ) / ! In-112 beta plus
     &  0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     &  0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     &  0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     &  3.638E-27, 2.175E-25, 8.378E-24, 4.877E-22, 1.758E-20,
     &  4.138E-19, 6.652E-18, 7.628E-17, 6.506E-16, 4.292E-15,
     &  2.269E-14, 9.911E-14, 3.676E-13, 3.382E-12, 2.042E-11,
     &  8.977E-11, 3.095E-10, 8.832E-10, 2.171E-09, 4.738E-09,
     &  1.716E-08, 4.759E-08, 1.093E-07, 2.184E-07, 3.931E-07,
     &  6.528E-07, 1.017E-06, 1.506E-06, 2.930E-06, 5.070E-06,
     &  8.952E-06, 1.439E-05, 2.159E-05, 3.075E-05, 4.202E-05,
     &  5.558E-05, 7.155E-05, 9.008E-05, 1.113E-04, 1.353E-04,
     &  1.920E-04, 2.609E-04, 3.426E-04, 4.376E-04, 5.461E-04,
     &  6.686E-04, 8.051E-04, 1.121E-03, 1.493E-03, 1.921E-03,
     &  2.405E-03, 2.942E-03, 3.530E-03, 4.168E-03, 4.853E-03,
     &  6.358E-03, 8.025E-03, 1.031E-02, 1.280E-02, 1.545E-02,
     &  1.825E-02, 2.116E-02, 2.419E-02, 2.729E-02, 3.048E-02,
     &  3.371E-02, 3.700E-02, 4.369E-02, 5.047E-02, 5.729E-02,
     &  6.412E-02, 7.092E-02, 7.768E-02, 8.437E-02, 9.750E-02,
     &  1.102E-01, 1.225E-01, 1.342E-01, 1.454E-01, 1.560E-01,
     &  1.660E-01, 1.755E-01, 1.926E-01, 2.073E-01, 2.222E-01,
     &  2.333E-01, 2.407E-01, 2.443E-01, 2.446E-01, 2.415E-01,
     &  2.353E-01, 2.264E-01, 2.149E-01, 2.013E-01, 1.691E-01,
     &  1.316E-01, 9.218E-02, 5.468E-02, 2.361E-02, 3.997E-03/

      data ( ddc2byld(35,i), i = 1, 99 ) / ! In-114 beta plus
     &  0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     &  0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     &  0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     &  1.213E-29, 7.251E-28, 2.793E-26, 1.625E-24, 5.858E-23,
     &  1.379E-21, 2.216E-20, 2.541E-19, 2.167E-18, 1.429E-17,
     &  7.552E-17, 3.299E-16, 1.223E-15, 1.125E-14, 6.792E-14,
     &  2.984E-13, 1.028E-12, 2.934E-12, 7.211E-12, 1.573E-11,
     &  5.693E-11, 1.578E-10, 3.620E-10, 7.231E-10, 1.301E-09,
     &  2.159E-09, 3.361E-09, 4.972E-09, 9.663E-09, 1.670E-08,
     &  2.943E-08, 4.723E-08, 7.076E-08, 1.006E-07, 1.372E-07,
     &  1.812E-07, 2.329E-07, 2.927E-07, 3.610E-07, 4.381E-07,
     &  6.196E-07, 8.392E-07, 1.098E-06, 1.398E-06, 1.739E-06,
     &  2.121E-06, 2.546E-06, 3.519E-06, 4.655E-06, 5.949E-06,
     &  7.393E-06, 8.979E-06, 1.070E-05, 1.254E-05, 1.450E-05,
     &  1.871E-05, 2.327E-05, 2.935E-05, 3.573E-05, 4.230E-05,
     &  4.898E-05, 5.568E-05, 6.233E-05, 6.889E-05, 7.530E-05,
     &  8.152E-05, 8.752E-05, 9.875E-05, 1.088E-04, 1.176E-04,
     &  1.251E-04, 1.312E-04, 1.360E-04, 1.394E-04, 1.424E-04,
     &  1.405E-04, 1.342E-04, 1.242E-04, 1.112E-04, 9.597E-05,
     &  7.928E-05, 6.203E-05, 2.951E-05, 6.257E-06/

      data ( ddc2byld(36,i), i = 1, 102 ) / ! Sb-122 beta plus
     &  0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     &  0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     &  0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     &  0.000E+00, 8.103E-29, 4.054E-27, 3.125E-25, 1.411E-23,
     &  3.957E-22, 7.265E-21, 9.206E-20, 8.464E-19, 5.909E-18,
     &  3.261E-17, 1.473E-16, 5.608E-16, 5.356E-15, 3.316E-14,
     &  1.484E-13, 5.188E-13, 1.497E-12, 3.714E-12, 8.168E-12,
     &  2.998E-11, 8.408E-11, 1.949E-10, 3.928E-10, 7.126E-10,
     &  1.191E-09, 1.868E-09, 2.780E-09, 5.466E-09, 9.542E-09,
     &  1.701E-08, 2.756E-08, 4.164E-08, 5.966E-08, 8.198E-08,
     &  1.089E-07, 1.409E-07, 1.780E-07, 2.207E-07, 2.691E-07,
     &  3.841E-07, 5.244E-07, 6.914E-07, 8.860E-07, 1.109E-06,
     &  1.361E-06, 1.642E-06, 2.292E-06, 3.060E-06, 3.943E-06,
     &  4.938E-06, 6.040E-06, 7.245E-06, 8.547E-06, 9.940E-06,
     &  1.298E-05, 1.631E-05, 2.082E-05, 2.564E-05, 3.070E-05,
     &  3.594E-05, 4.131E-05, 4.675E-05, 5.224E-05, 5.774E-05,
     &  6.322E-05, 6.866E-05, 7.934E-05, 8.967E-05, 9.957E-05,
     &  1.090E-04, 1.179E-04, 1.262E-04, 1.340E-04, 1.478E-04,
     &  1.594E-04, 1.687E-04, 1.756E-04, 1.802E-04, 1.825E-04,
     &  1.824E-04, 1.798E-04, 1.672E-04, 1.452E-04, 1.062E-04,
     &  6.009E-05, 1.830E-05/

      data ( ddc2byld(37,i), i = 1, 111 ) / ! I-126 beta plus
     &  0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     &  0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     &  0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     &  0.000E+00, 1.477E-27, 9.766E-26, 1.004E-23, 5.662E-22,
     &  1.872E-20, 3.877E-19, 5.358E-18, 5.245E-17, 3.831E-16,
     &  2.186E-15, 1.012E-14, 3.924E-14, 3.843E-13, 2.417E-12,
     &  1.094E-11, 3.854E-11, 1.119E-10, 2.792E-10, 6.171E-10,
     &  2.285E-09, 6.456E-09, 1.507E-08, 3.057E-08, 5.580E-08,
     &  9.383E-08, 1.479E-07, 2.213E-07, 4.391E-07, 7.728E-07,
     &  1.390E-06, 2.272E-06, 3.458E-06, 4.989E-06, 6.898E-06,
     &  9.221E-06, 1.199E-05, 1.523E-05, 1.896E-05, 2.323E-05,
     &  3.342E-05, 4.597E-05, 6.102E-05, 7.868E-05, 9.905E-05,
     &  1.222E-04, 1.482E-04, 2.088E-04, 2.810E-04, 3.648E-04,
     &  4.600E-04, 5.662E-04, 6.830E-04, 8.100E-04, 9.466E-04,
     &  1.247E-03, 1.579E-03, 2.032E-03, 2.520E-03, 3.036E-03,
     &  3.573E-03, 4.124E-03, 4.685E-03, 5.251E-03, 5.818E-03,
     &  6.382E-03, 6.941E-03, 8.033E-03, 9.077E-03, 1.006E-02,
     &  1.098E-02, 1.182E-02, 1.259E-02, 1.328E-02, 1.442E-02,
     &  1.525E-02, 1.578E-02, 1.605E-02, 1.608E-02, 1.590E-02,
     &  1.555E-02, 1.507E-02, 1.389E-02, 1.271E-02, 1.180E-02,
     &  1.181E-02, 1.171E-02, 1.141E-02, 1.090E-02, 1.017E-02,
     &  9.234E-03, 8.090E-03, 6.762E-03, 5.297E-03, 2.292E-03,
     &  1.857E-04/

      data ( ddc2byld(38,i), i = 1, 92 ) / ! I-128 beta plus
     &  0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     &  0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     &  0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     &  0.000E+00, 7.968E-29, 5.267E-27, 5.414E-25, 3.052E-23,
     &  1.009E-21, 2.089E-20, 2.886E-19, 2.824E-18, 2.063E-17,
     &  1.176E-16, 5.445E-16, 2.111E-15, 2.066E-14, 1.299E-13,
     &  5.875E-13, 2.069E-12, 6.005E-12, 1.497E-11, 3.308E-11,
     &  1.223E-10, 3.453E-10, 8.052E-10, 1.632E-09, 2.975E-09,
     &  4.998E-09, 7.869E-09, 1.176E-08, 2.329E-08, 4.091E-08,
     &  7.341E-08, 1.196E-07, 1.817E-07, 2.614E-07, 3.605E-07,
     &  4.806E-07, 6.232E-07, 7.895E-07, 9.808E-07, 1.198E-06,
     &  1.715E-06, 2.347E-06, 3.098E-06, 3.974E-06, 4.976E-06,
     &  6.105E-06, 7.363E-06, 1.026E-05, 1.366E-05, 1.753E-05,
     &  2.185E-05, 2.659E-05, 3.170E-05, 3.714E-05, 4.288E-05,
     &  5.509E-05, 6.800E-05, 8.468E-05, 1.015E-04, 1.179E-04,
     &  1.337E-04, 1.485E-04, 1.620E-04, 1.740E-04, 1.845E-04,
     &  1.932E-04, 2.002E-04, 2.088E-04, 2.102E-04, 2.049E-04,
     &  1.935E-04, 1.769E-04, 1.562E-04, 1.324E-04, 8.103E-05,
     &  3.420E-05, 4.377E-06/

      data ( ddc2byld(39,i), i = 1, 117 ) / ! Cs-130 beta plus
     &  0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     &  0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     &  0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     &  0.000E+00, 0.000E+00, 1.170E-25, 1.605E-23, 1.119E-21,
     &  4.295E-20, 9.862E-19, 1.463E-17, 1.503E-16, 1.135E-15,
     &  6.622E-15, 3.114E-14, 1.221E-13, 1.211E-12, 7.675E-12,
     &  3.487E-11, 1.233E-10, 3.589E-10, 8.978E-10, 1.989E-09,
     &  7.403E-09, 2.103E-08, 4.936E-08, 1.007E-07, 1.848E-07,
     &  3.124E-07, 4.949E-07, 7.442E-07, 1.491E-06, 2.647E-06,
     &  4.809E-06, 7.931E-06, 1.218E-05, 1.770E-05, 2.466E-05,
     &  3.319E-05, 4.343E-05, 5.550E-05, 6.954E-05, 8.565E-05,
     &  1.245E-04, 1.730E-04, 2.317E-04, 3.013E-04, 3.824E-04,
     &  4.755E-04, 5.809E-04, 8.301E-04, 1.132E-03, 1.488E-03,
     &  1.898E-03, 2.363E-03, 2.881E-03, 3.453E-03, 4.076E-03,
     &  5.474E-03, 7.062E-03, 9.296E-03, 1.178E-02, 1.449E-02,
     &  1.740E-02, 2.049E-02, 2.374E-02, 2.713E-02, 3.064E-02,
     &  3.426E-02, 3.798E-02, 4.566E-02, 5.358E-02, 6.170E-02,
     &  6.995E-02, 7.829E-02, 8.670E-02, 9.512E-02, 1.120E-01,
     &  1.287E-01, 1.452E-01, 1.613E-01, 1.771E-01, 1.925E-01,
     &  2.075E-01, 2.219E-01, 2.493E-01, 2.745E-01, 3.029E-01,
     &  3.276E-01, 3.487E-01, 3.659E-01, 3.794E-01, 3.891E-01,
     &  3.952E-01, 3.976E-01, 3.964E-01, 3.919E-01, 3.732E-01,
     &  3.431E-01, 3.035E-01, 2.569E-01, 2.061E-01, 1.541E-01,
     &  1.038E-01, 2.411E-02/

      data ( ddc2byld(40,i), i = 1, 99 ) / ! Cs-132 beta plus
     &  0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     &  0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     &  0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     &  0.000E+00, 0.000E+00, 4.533E-26, 6.221E-24, 4.333E-22,
     &  1.663E-20, 3.819E-19, 5.665E-18, 5.818E-17, 4.392E-16,
     &  2.563E-15, 1.205E-14, 4.722E-14, 4.685E-13, 2.967E-12,
     &  1.348E-11, 4.762E-11, 1.386E-10, 3.466E-10, 7.676E-10,
     &  2.855E-09, 8.105E-09, 1.901E-08, 3.875E-08, 7.106E-08,
     &  1.200E-07, 1.900E-07, 2.855E-07, 5.712E-07, 1.013E-06,
     &  1.837E-06, 3.024E-06, 4.634E-06, 6.726E-06, 9.352E-06,
     &  1.256E-05, 1.641E-05, 2.094E-05, 2.618E-05, 3.219E-05,
     &  4.664E-05, 6.454E-05, 8.614E-05, 1.116E-04, 1.412E-04,
     &  1.749E-04, 2.128E-04, 3.019E-04, 4.087E-04, 5.332E-04,
     &  6.751E-04, 8.341E-04, 1.009E-03, 1.201E-03, 1.407E-03,
     &  1.860E-03, 2.362E-03, 3.048E-03, 3.785E-03, 4.562E-03,
     &  5.366E-03, 6.186E-03, 7.015E-03, 7.843E-03, 8.665E-03,
     &  9.473E-03, 1.026E-02, 1.177E-02, 1.316E-02, 1.440E-02,
     &  1.549E-02, 1.642E-02, 1.719E-02, 1.778E-02, 1.846E-02,
     &  1.850E-02, 1.795E-02, 1.689E-02, 1.538E-02, 1.354E-02,
     &  1.145E-02, 9.226E-03, 4.855E-03, 1.417E-03/

      data ( ddc2byld(41,i), i = 1, 112 ) / ! Eu-150m beta plus
     &  0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     &  0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     &  0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     &  0.000E+00, 0.000E+00, 3.924E-28, 1.363E-25, 1.660E-23,
     &  8.739E-22, 2.364E-20, 3.770E-19, 3.948E-18, 2.949E-17,
     &  1.676E-16, 7.614E-16, 2.878E-15, 2.660E-14, 1.585E-13,
     &  6.852E-13, 2.329E-12, 6.581E-12, 1.609E-11, 3.505E-11,
     &  1.278E-10, 3.598E-10, 8.434E-10, 1.727E-09, 3.190E-09,
     &  5.439E-09, 8.703E-09, 1.322E-08, 2.708E-08, 4.916E-08,
     &  9.176E-08, 1.552E-07, 2.439E-07, 3.626E-07, 5.155E-07,
     &  7.073E-07, 9.422E-07, 1.225E-06, 1.559E-06, 1.949E-06,
     &  2.913E-06, 4.148E-06, 5.683E-06, 7.546E-06, 9.762E-06,
     &  1.235E-05, 1.534E-05, 2.259E-05, 3.163E-05, 4.256E-05,
     &  5.545E-05, 7.036E-05, 8.731E-05, 1.063E-04, 1.274E-04,
     &  1.755E-04, 2.316E-04, 3.122E-04, 4.036E-04, 5.051E-04,
     &  6.155E-04, 7.341E-04, 8.598E-04, 9.918E-04, 1.129E-03,
     &  1.272E-03, 1.418E-03, 1.721E-03, 2.033E-03, 2.351E-03,
     &  2.672E-03, 2.993E-03, 3.313E-03, 3.629E-03, 4.246E-03,
     &  4.836E-03, 5.392E-03, 5.912E-03, 6.392E-03, 6.831E-03,
     &  7.229E-03, 7.586E-03, 8.177E-03, 8.614E-03, 8.961E-03,
     &  9.117E-03, 9.119E-03, 8.949E-03, 8.602E-03, 8.096E-03,
     &  7.452E-03, 6.695E-03, 5.854E-03, 4.959E-03, 3.139E-03,
     &  1.497E-03, 3.411E-04/

      data ( ddc2byld(42,i), i = 1, 105 ) / ! Eu-152 beta plus
     &  0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     &  0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     &  0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     &  0.000E+00, 0.000E+00, 4.028E-29, 1.399E-26, 1.704E-24,
     &  8.970E-23, 2.426E-21, 3.870E-20, 4.052E-19, 3.026E-18,
     &  1.719E-17, 7.812E-17, 2.952E-16, 2.728E-15, 1.625E-14,
     &  7.027E-14, 2.388E-13, 6.747E-13, 1.649E-12, 3.592E-12,
     &  1.309E-11, 3.686E-11, 8.637E-11, 1.768E-10, 3.265E-10,
     &  5.565E-10, 8.901E-10, 1.352E-09, 2.767E-09, 5.021E-09,
     &  9.365E-09, 1.583E-08, 2.486E-08, 3.692E-08, 5.246E-08,
     &  7.193E-08, 9.575E-08, 1.244E-07, 1.582E-07, 1.976E-07,
     &  2.949E-07, 4.194E-07, 5.737E-07, 7.607E-07, 9.826E-07,
     &  1.242E-06, 1.540E-06, 2.260E-06, 3.155E-06, 4.232E-06,
     &  5.498E-06, 6.955E-06, 8.604E-06, 1.044E-05, 1.247E-05,
     &  1.708E-05, 2.240E-05, 2.995E-05, 3.841E-05, 4.767E-05,
     &  5.761E-05, 6.813E-05, 7.911E-05, 9.046E-05, 1.021E-04,
     &  1.139E-04, 1.259E-04, 1.499E-04, 1.738E-04, 1.970E-04,
     &  2.194E-04, 2.406E-04, 2.606E-04, 2.791E-04, 3.115E-04,
     &  3.373E-04, 3.564E-04, 3.688E-04, 3.748E-04, 3.749E-04,
     &  3.694E-04, 3.590E-04, 3.262E-04, 2.824E-04, 2.226E-04,
     &  1.702E-04, 1.195E-04, 7.069E-05, 3.014E-05, 4.811E-06/

      data ( ddc2byld(43,i), i = 1, 108 ) / ! Eu-152m beta plus
     &  0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     &  0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     &  0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     &  0.000E+00, 0.000E+00, 1.004E-29, 3.486E-27, 4.247E-25,
     &  2.236E-23, 6.048E-22, 9.647E-21, 1.010E-19, 7.545E-19,
     &  4.287E-18, 1.948E-17, 7.362E-17, 6.805E-16, 4.055E-15,
     &  1.753E-14, 5.959E-14, 1.684E-13, 4.116E-13, 8.966E-13,
     &  3.268E-12, 9.203E-12, 2.157E-11, 4.416E-11, 8.157E-11,
     &  1.391E-10, 2.225E-10, 3.381E-10, 6.921E-10, 1.256E-09,
     &  2.344E-09, 3.963E-09, 6.228E-09, 9.254E-09, 1.315E-08,
     &  1.804E-08, 2.403E-08, 3.122E-08, 3.973E-08, 4.966E-08,
     &  7.417E-08, 1.056E-07, 1.445E-07, 1.918E-07, 2.480E-07,
     &  3.136E-07, 3.892E-07, 5.723E-07, 8.003E-07, 1.075E-06,
     &  1.399E-06, 1.773E-06, 2.198E-06, 2.673E-06, 3.198E-06,
     &  4.396E-06, 5.784E-06, 7.772E-06, 1.002E-05, 1.249E-05,
     &  1.517E-05, 1.804E-05, 2.106E-05, 2.421E-05, 2.747E-05,
     &  3.083E-05, 3.427E-05, 4.130E-05, 4.847E-05, 5.567E-05,
     &  6.283E-05, 6.990E-05, 7.681E-05, 8.355E-05, 9.632E-05,
     &  1.080E-04, 1.186E-04, 1.279E-04, 1.358E-04, 1.425E-04,
     &  1.479E-04, 1.520E-04, 1.564E-04, 1.560E-04, 1.492E-04,
     &  1.364E-04, 1.185E-04, 9.704E-05, 7.355E-05, 5.013E-05,
     &  2.935E-05, 1.388E-05, 3.624E-06/

      data ( ddc2byld(44,i), i = 1, 102 ) / ! Tm-168 beta plus
     &  0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     &  0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     &  0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00, 0.000E+00,
     &  0.000E+00, 0.000E+00, 5.912E-29, 2.181E-26, 2.376E-24,
     &  1.046E-22, 2.324E-21, 3.067E-20, 2.701E-19, 1.729E-18,
     &  8.583E-18, 3.467E-17, 1.183E-16, 9.266E-16, 4.870E-15,
     &  1.914E-14, 6.042E-14, 1.612E-13, 3.769E-13, 7.924E-13,
     &  2.747E-12, 7.497E-12, 1.723E-11, 3.487E-11, 6.400E-11,
     &  1.088E-10, 1.741E-10, 2.649E-10, 5.459E-10, 1.000E-09,
     &  1.890E-09, 3.240E-09, 5.161E-09, 7.771E-09, 1.119E-08,
     &  1.553E-08, 2.092E-08, 2.747E-08, 3.532E-08, 4.458E-08,
     &  6.782E-08, 9.811E-08, 1.364E-07, 1.835E-07, 2.403E-07,
     &  3.075E-07, 3.859E-07, 5.788E-07, 8.237E-07, 1.124E-06,
     &  1.484E-06, 1.904E-06, 2.386E-06, 2.932E-06, 3.542E-06,
     &  4.950E-06, 6.605E-06, 9.005E-06, 1.174E-05, 1.479E-05,
     &  1.810E-05, 2.165E-05, 2.540E-05, 2.931E-05, 3.335E-05,
     &  3.749E-05, 4.170E-05, 5.024E-05, 5.877E-05, 6.714E-05,
     &  7.522E-05, 8.291E-05, 9.012E-05, 9.680E-05, 1.084E-04,
     &  1.173E-04, 1.235E-04, 1.271E-04, 1.280E-04, 1.263E-04,
     &  1.224E-04, 1.164E-04, 9.927E-05, 7.731E-05, 4.737E-05,
     &  2.033E-05, 2.895E-06/

*-----------------------------------------------------------------------

      end


