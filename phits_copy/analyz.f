************************************************************************
*                                                                      *
      subroutine analyz(ncol,mark)
*                                                                      *
*       main control routine of analysis                               *
*       modified by K.Niita on 2011/06/14                              *
*                                                                      *
************************************************************************
C for USE_MOD_COUNTER
      use mod_counter, only: rncnt,rnint,rnintr,rnpnt,rnpntr
     &                      ,aevts,aevtr,bevts,bevtr
     &                      ,iaevt, ibevt, jaevt, jbevt
C for REDUCTION_COUNTER
!$   &                      ,rncnt2,rnint2,rnintr2,rnpnt2,rnpntr2
!$   &                      ,aevts2,aevtr2,bevts2,bevtr2
      use TETRAMOD, only: tetrafin
      use TDCHAINMOD, only: itdc,pdchreg2,talldcfin !FURUTA20200522
      use MMBANKMOD             !FURUTA

      use NDATA2MOD
      use moddas_region
      use moddas_bends, only: bends
      use t4dtrack_mod, only: it4dtrack, t4dtrack
*-----------------------------------------------------------------------

      implicit real*8(a-h,o-z)

*-----------------------------------------------------------------------

      include 'param00.inc'
      include 'param.inc'
      include 'err.inc'

*-----------------------------------------------------------------------

      common /mpi00/ npe, me
      common /mpi01/ iccp(20000)
      common /mpi02/ nmbch0
      integer nrmnbch0
      common /batchprocess/ nrmnbch0
      common /tall00/ itnm, ital(itlmax), itals(itlmax), italm(itlmax)

*-----------------------------------------------------------------------

      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
!$OMP THREADPRIVATE(/icomon/)
      common /randn/ nrandgen ! S.H. xorshift (2020.2.6)
      integer*8 :: iranji64 ! S.H. xorshift (2020.2.6)
      common /randm4/ rnfb,rnfs,rngb,rngs,rnmult,ranj,rani,
     &                rnrtc,nstrid,inif, iranji64
      common /irad/   irands,irandf,nseed !FURUTA
!$OMP THREADPRIVATE(/irad/)
      common /ncall0/ ncall               !FURUTA
      common /randsv/ srijk,rrijk         !FURUTA
      common /cparm/  maxbch,maxcas
      common /jcomon/ nabov,nobch,nocas,nomax
!$OMP THREADPRIVATE(/jcomon/)
      common /eparm/  esmax, esmin, emin(20)
      common /rcomon/ rcasc
      common /neulo/  reutn
      common /inout/  in,io

      common /isomul/ smlwt(isrc), totfact, imsrc
      common /isorsf/ isorf(isrc),lsfile(isrc), sfile(isrc)
      character sfile*100
      common /isorst/ jstyp(isrc), istyp(isrc), inkf0(isrc), lstyp(isrc)
!$OMP THREADPRIVATE(/isorst/)

      common /wtsave/ oldwt
!$OMP THREADPRIVATE(/wtsave/)
      common /ngcut/  incut, igcut, ipcut
      common /ncutne/ wtneut,rtneut
      common /ncutgm/ wtgamm,rtgamm
      common /ncuten/ wtelen,rtelen
      common /ncutep/ wtelep,rtelep
      common /ncutpr/ wtprot,rtprot
      common /regdc/  idrg(kvlmax), idgr(kvmmax)



      common /clustt/ nclsts, iclusts(nnn)
!$OMP THREADPRIVATE(/clustt/)
      common /clustw/ jclusts(0:8,nnn),  qclusts(0:12,nnn)
!$OMP THREADPRIVATE(/clustw/)
      common /clustl/ rumpal(0:20), numpal(0:20)
!$OMP THREADPRIVATE(/clustl/)

      common /prodp/  aprodp(20), bprodp(20)
      common /dcayp/  adcayp(20), bdcayp(20)
      common /stopp/  astopp(20), bstopp(20)
      common /timep/  atimep(20), btimep(20)
      common /leakp/  aleakp(20), bleakp(20)

      common /otherp/ nothp(2000,2), nothn
      common /othjmp/ nojmp(2000,2), nojmn
      common /othstp/ nostp(2000,2), nostn
      common /othtim/ notip(2000,2), notin
      common /othlkp/ nolkp(2000,2), nolkn
      common /othdcp/ nodcp(2000,6), nodcn

      common /kounts/ rount(20)

      common /inpec/ ititl, ipara, ibody, iregn, llarr, itby, itar
      common /kmat1a/ mxmat, mxmat0, mxnel
      common /talout/ itall
      common /taliin/ rsouin, nzztin, nrgnin
      common /kcomsi/ iibnk, jjbnk               !FURUTA
!$OMP THREADPRIVATE(/kcomsi/)
      common /kcomsm/ imbnk, jmbnk, itbnk, jtbnk !FURUTA

      common /tlgeom/ iblz1,iblz2
!$OMP THREADPRIVATE(/tlgeom/)
      common /tlglat/ ilev1,ilev2,ilat1(5,10),ilat2(5,10)
!$OMP THREADPRIVATE(/tlglat/)
      common /cntreg/ ncntc(3), ncreg(3), incrc(3), incrt(3), kcont(3)
      common /cntpat/ icpan(3), icpat(3,20,2)
      common /dumpall/ idumpall
      common /cusrtally/ iusrtally, iudtf(50)
      common /tmtreg/ ntmrg, intmc, intmt, ktime

      common /tcntl/  icntl, inucr

      common /sumbnd/ nbnd(nbchmax) !FURUTA20220106

      common /paran/ icfn(100), ilfn(100), chfn(100)
      character chfn*200



*-----------------------------------------------------------------------

      dimension data(11)

      dimension bend(10)
      common / stdstop / idstop, instop
      common /tall67/ itstd(itlmax), rtstd(itlmax)

*-----------------------------------------------------------------------

      data rount/20*0.0d0/
      data iibnk/0/,jjbnk/0/ !FURUTA

      common /randm5/ rijklst,rijkinit !OBINATA(2012.6.18)

      integer irndmode,idmprijk
      common /irndm/ irndmode,idmprijk

      common /infprint/ infout

      common /stat / istdev, irestart, ireschk
      common /res01/ istdevres,maxcasres,rijklstres,irdrf
      common /res02/ crdrfln, irdrfll
      character crdrfln*100

      common /res04/ lrijkeqrf
      logical lrijkeqrf

*-----------------------------------------------------------------------
      common /dmpinfo/ rsouinbch,maxbchdmp,maxcasdmp
      common /stat2/ dmpmulti,idmpmode,ibchjmp,idmpjmp(2)

*-----------------------------------------------------------------------
cFURUTA20150714 TETRA
      integer nlat3,itetcl(10),itetvol,itetauto,itgchk
      common /tetf1/ nlat3,itetcl,itetvol,itetauto,itgchk

*-----------------------------------------------------------------------
      common /timecut/ timeout

      common /  tscmsg   / ktsc(kvlmax), mntsc, ntsc(kvlmax)
      common /  tscreg   / ntscell(kvlmax)
      common / etsminmax / etsmin, etsmax

*-----------------------------------------------------------------------
      common /cntmx/ ncntmx(3), nctmxsr(3)
!$OMP THREADPRIVATE(/cntmx/)

*-----------------------------------------------------------------------
*     count ncol in analyz
*-----------------------------------------------------------------------

      if( ncol .ge. 1 .and. ncol .le. 16 ) then
!$OMP CRITICAL (rount_crit)
         rount( ncol ) = rount( ncol ) + 1.0
!$OMP END CRITICAL (rount_crit)

      end if

*-----------------------------------------------------------------------
*     dumpall
*-----------------------------------------------------------------------

         if( abs(idumpall) .eq. 1 .and. icntl .ne. 12 .and.
     &       ncol .ne. 0 .and. ncol .ne. 101 ) then

            call dumpdat(ncol)

         end if

*-----------------------------------------------------------------------
*     user defined tally
*-----------------------------------------------------------------------

         if( iusrtally .eq. 1 ) then

!$OMP CRITICAL (usrtally_crit)
            call usrtally(ncol)
!$OMP END CRITICAL (usrtally_crit)

         end if


         if( it4dtrack > 0 ) then
            call t4dtrack(ncol)
         endif


*-----------------------------------------------------------------------
*     geometry errors and weight cut-off
*-----------------------------------------------------------------------

      if( ncol .ge. 5 .and. ncol .le. 8 ) return


*-----------------------------------------------------------------------
*     ncol = 1 : start of calulation
*-----------------------------------------------------------------------

      if( ncol .eq. 1 ) then

               call timex(dummy,dummy,0)

               ncall  = 0

               wtneut = 0.0d0
               rtneut = 0.0d0

               wtgamm = 0.0d0
               rtgamm = 0.0d0

               wtelen = 0.0d0
               rtelen = 0.0d0

               wtelep = 0.0d0
               rtelep = 0.0d0

               wtprot = 0.0d0
               rtprot = 0.0d0


            do i = 1, 20

               bprodp(i) = 0.0
               aprodp(i) = 0.0

               bdcayp(i) = 0.0
               adcayp(i) = 0.0

               bstopp(i) = 0.0
               astopp(i) = 0.0

               btimep(i) = 0.0
               atimep(i) = 0.0

               bleakp(i) = 0.0
               aleakp(i) = 0.0

            end do

            do i = 1, 200
            do j = 1, 2

               nothp(i,j) = 0
               nojmp(i,j) = 0
               nostp(i,j) = 0
               notip(i,j) = 0
               nolkp(i,j) = 0

            end do
            do j = 1, 6

               nodcp(i,j) = 0

            end do
            end do

               nothn = 0
               nojmn = 0
               nostn = 0
               notin = 0
               nolkn = 0

               nodcn = 0

               reutn = 0.0

               imbnk = 0
               jmbnk = 0
               itbnk = 0
               jtbnk = 0

*-----------------------------------------------------------------------

      if( me .eq. 0 ) then

            if (irestart .ne. 0 ) then
                write(*,
     &           '('' This is restart calculation.'',
     &           '' Initial random seed is read from '', a )')
     &           crdrfln(1:irdrfll)
            end if

            if (irestart .ne. 0 .and. .not. lrijkeqrf ) then

                ErrCha = ''
                ErrID = 'L:316/R:analyz/F:analyz.f' !W00_001_001
                call ErrWrite(ErrID,ErrCha)

                write(*,
     &                  '('' Warning: This random seed is different'',
     &                    '' from those written in other resfiles.'')')
            end if

      end if

*-----------------------------------------------------------------------
*     ncol = 4 : start of new source
*                summary of bank access
*                ncut, pcut and gcut for source particles
*-----------------------------------------------------------------------

      else if( ncol .eq. 4 ) then
!$OMP CRITICAL (imbnk_crit)
                  if( iibnk .gt. imbnk ) imbnk = iibnk
                  if( jjbnk .gt. jmbnk ) jmbnk = jjbnk
!$OMP END CRITICAL (imbnk_crit)
                  iibnk = 0
                  jjbnk = 0

                  call cputime(3)

               call talls01(ncol)

               ncntmx(:) = nctmxsr(:) !frtati 2021/03/05
                  rncnt(3) = rncnt(3) + 1.0
                  call cputime(3)

*-----------------------------------------------------------------------
*     ncol = 9 : booking of stopped particles by time cutoff
*-----------------------------------------------------------------------

      else if( ncol .eq. 9 ) then
!$OMP CRITICAL (timep_crit)
                  btimep(ityp) = btimep(ityp) + 1.0
                  atimep(ityp) = atimep(ityp) + oldwt
!$OMP END CRITICAL (timep_crit)

               if( ityp .eq. 11 ) then

                  kf = nkf(ibknkf+no,ipomp+1)
!$OMP CRITICAL (tin_crit)
                  if( notin .eq. 0 ) then

                        notin = notin + 1
                        notip(notin,1) = 1
                        notip(notin,2) = kf

                  else

                     do m = 1, notin

                        if( kf .eq. notip(m,2) ) then

                           notip(m,1) = notip(m,1) + 1
                           goto 220

                        end if

                     end do

                        notin = notin + 1
                        notip(notin,1) = 1
                        notip(notin,2) = kf

  220                continue

                  end if
!$OMP END CRITICAL (tin_crit)
               end if

*-----------------------------------------------------------------------
*     ncol = 11 : booking of stopped particles by energy cutoff
*-----------------------------------------------------------------------

      else if( ncol .eq. 11 ) then
!$OMP CRITICAL (stopp_crit)
                  bstopp(ityp) = bstopp(ityp) + 1.0
                  astopp(ityp) = astopp(ityp) + oldwt
!$OMP END CRITICAL (stopp_crit)
               if( ityp .eq. 11 ) then

                  kf = nkf(ibknkf+no,ipomp+1)
!$OMP CRITICAL (stn_crit)
                  if( nostn .eq. 0 ) then

                        nostn = nostn + 1
                        nostp(nostn,1) = 1
                        nostp(nostn,2) = kf

                  else

                     do m = 1, nostn

                        if( kf .eq. nostp(m,2) ) then

                           nostp(m,1) = nostp(m,1) + 1
                           goto 210

                        end if

                     end do

                        nostn = nostn + 1
                        nostp(nostn,1) = 1
                        nostp(nostn,2) = kf

  210                continue

                  end if
!$OMP END CRITICAL (stn_crit)
               end if

*-----------------------------------------------------------------------
*        write ncut, pcut and gcut file
*-----------------------------------------------------------------------
!$OMP CRITICAL (ncut_crit)
            if( ityp .eq. 2 .and. incut .ne. 0 ) then

                  rtneut = rtneut + 1.0
                  wtneut = wtneut + wt(ibkwt+no,ipomp+1)

                  data(1)  = ec(ibkec+no,ipomp+1)
                  data(2)  = xc(ibkxc+no,ipomp+1)
                  data(3)  = yc(ibkyc+no,ipomp+1)
                  data(4)  = zc(ibkzc+no,ipomp+1)
                  data(5)  = u(ibku+no,ipomp+1)
                  data(6)  = v(ibkv+no,ipomp+1)
                  data(7)  = w(ibkw+no,ipomp+1)
                  data(8)  = wt(ibkwt+no,ipomp+1)
                  data(9)  = 1.0d0
                  data(10) = iblz1
                  data(11) = abs(tc(ibktc+no,ipomp+1))

                  call wrnt12(data)

            end if

            if( ( ityp .eq. 14 .and.
     &            igcut .ge. 1 .and. igcut .le. 2  ) .or.
     &          ( ityp .ge. 12 .and. ityp .le. 14 .and.
     &            igcut .eq. 3 ) ) then

               if( ityp .eq. 12 ) then

                  dii = 4.0d0
                  rtelen = rtelen + 1.0
                  wtelen = wtelen + wt(ibkwt+no,ipomp+1)

               else if( ityp .eq. 13 ) then

                  dii = 5.0d0
                  rtelep = rtelep + 1.0
                  wtelep = wtelep + wt(ibkwt+no,ipomp+1)

               else if( ityp .eq. 14 ) then

                  dii = 3.0d0
                  rtgamm = rtgamm + 1.0
                  wtgamm = wtgamm + wt(ibkwt+no,ipomp+1)

               end if

                  data(1)  = ec(ibkec+no,ipomp+1)
                  data(2)  = xc(ibkxc+no,ipomp+1)
                  data(3)  = yc(ibkyc+no,ipomp+1)
                  data(4)  = zc(ibkzc+no,ipomp+1)
                  data(5)  = u(ibku+no,ipomp+1)
                  data(6)  = v(ibkv+no,ipomp+1)
                  data(7)  = w(ibkw+no,ipomp+1)
                  data(8)  = wt(ibkwt+no,ipomp+1)
                  data(9)  = dii
                  data(10) = iblz1
                  data(11) = abs(tc(ibktc+no,ipomp+1))

                  call wrnt13(data)

            end if

            if( ityp .eq. 1 .and. ipcut .ne. 0 ) then

                  rtprot = rtprot + 1.0
                  wtprot = wtprot + wt(ibkwt+no,ipomp+1)

                  data(1)  = ec(ibkec+no,ipomp+1)
                  data(2)  = xc(ibkxc+no,ipomp+1)
                  data(3)  = yc(ibkyc+no,ipomp+1)
                  data(4)  = zc(ibkzc+no,ipomp+1)
                  data(5)  = u(ibku+no,ipomp+1)
                  data(6)  = v(ibkv+no,ipomp+1)
                  data(7)  = w(ibkw+no,ipomp+1)
                  data(8)  = wt(ibkwt+no,ipomp+1)
                  data(9)  = 2.0d0
                  data(10) = iblz1
                  data(11) = abs(tc(ibktc+no,ipomp+1))

                  call wrnt10(data)

            end if
!$OMP END CRITICAL (ncut_crit)
*-----------------------------------------------------------------------
*     ncol = 12 : booking of leakage particles
*-----------------------------------------------------------------------

      else if( ncol .eq. 12 ) then
!$OMP CRITICAL (leakp_crit)
                  bleakp(ityp) = bleakp(ityp) + 1.0
                  aleakp(ityp) = aleakp(ityp) + oldwt
!$OMP END CRITICAL (leakp_crit)

               if( ityp .eq. 11 ) then

                  kf = nkf(ibknkf+no,ipomp+1)
!$OMP CRITICAL (lkn_crit)
                  if( nolkn .eq. 0 ) then

                        nolkn = nolkn + 1
                        nolkp(nolkn,1) = 1
                        nolkp(nolkn,2) = kf

                  else

                     do m = 1, nolkn

                        if( kf .eq. nolkp(m,2) ) then

                           nolkp(m,1) = nolkp(m,1) + 1
                           goto 300

                        end if

                     end do

                        nolkn = nolkn + 1
                        nolkp(nolkn,1) = 1
                        nolkp(nolkn,2) = kf

  300                continue

                  end if
!$OMP END CRITICAL (lkn_crit)
               end if

*-----------------------------------------------------------------------
*     ncol = 13, 14 :  after nuclear reactions
*-----------------------------------------------------------------------

      else if( ncol .eq. 13 .or. ncol .eq. 14 ) then

c*-----------------------------------------------------------------------
c*        collision type check, apsorption or fission
c*-----------------------------------------------------------------------
*-----------------------------------------------------------------------
*        booking and data down of produced particles
*-----------------------------------------------------------------------

         if( nclsts .gt. 0 ) then

         do i = 1, nclsts

               knn = jclusts(3,i)
               kf  = jclusts(7,i)

*-----------------------------------------------------------------------

            if( mntsc .gt. 0 .and.
     &          ntscell( idgr(iblz(ibkblz+no,ipomp+1)) ) .ne. 0 .and.
     &        ( knn .eq. 12 .or. knn .eq. 13 ) .and.
     &        ( e(ibke+no,ipomp+1) .le. etsmax .and.
     &          e(ibke+no,ipomp+1) .ge. etsmin ) ) then

               emint = etsmin
            else
               emint = emin(knn) * dble(max(1,ibryf(knn,kf)))
            endif

*-----------------------------------------------------------------------
*                 booking of the other particles
*-----------------------------------------------------------------------

                  if( knn .eq. 11 ) then
!$OMP CRITICAL (thn_crit)
                     if( nothn .eq. 0 ) then

                           nothn = nothn + 1
                           nothp(nothn,1) = 1
                           nothp(nothn,2) = kf

                     else

                        do m = 1, nothn

                           if( kf .eq. nothp(m,2) ) then

                              nothp(m,1) = nothp(m,1) + 1
                              goto 310

                           end if

                        end do

                           nothn = nothn + 1
                           nothp(nothn,1) = 1
                           nothp(nothn,2) = kf

  310                   continue

                     end if
!$OMP END CRITICAL (thn_crit)
                  end if

*-----------------------------------------------------------------------
*           booking of stopped particles
*-----------------------------------------------------------------------

            if( jclusts(4,i) .eq. -1 ) then
!$OMP CRITICAL (stopp_crit)
                     bstopp(knn) = bstopp(knn) + 1.0
                     astopp(knn) = astopp(knn) + qclusts(8,i)
!$OMP END CRITICAL (stopp_crit)
                  if( knn .eq. 11 ) then
!$OMP CRITICAL (stn_crit)
                     if( nostn .eq. 0 ) then

                           nostn = nostn + 1
                           nostp(nostn,1) = 1
                           nostp(nostn,2) = kf

                     else

                        do m = 1, nostn

                           if( kf .eq. nostp(m,2) ) then

                              nostp(m,1) = nostp(m,1) + 1
                              goto 320

                           end if

                        end do

                           nostn = nostn + 1
                           nostp(nostn,1) = 1
                           nostp(nostn,2) = kf

  320                   continue

                     end if
!$OMP END CRITICAL (stn_crit)
                  end if

            end if

*-----------------------------------------------------------------------
*           write ncut, pcut and gcut file
*-----------------------------------------------------------------------

            if( jclusts(4,i) .eq. -1 .and.
     &          qclusts(8,i) .gt. 0.0 ) then

*-----------------------------------------------------------------------
*              write the information of low energy neutrons
*              for MCNP on 12
*-----------------------------------------------------------------------
!$OMP CRITICAL (ncut_crit)
               if( knn .eq. 2 .and. incut .ne. 0 .and.
     &             qclusts(7,i) .gt. 0.0 .and.
     &             qclusts(7,i) .le. emint ) then

                     jclusts(4,i) = -2

                     rtneut = rtneut + 1.0
                     wtneut = wtneut + qclusts(8,i)

                     data(1)  = qclusts(7,i)
                     data(2)  = qclusts(10,i)
                     data(3)  = qclusts(11,i)
                     data(4)  = qclusts(12,i)
                     data(5)  = qclusts(1,i)
                     data(6)  = qclusts(2,i)
                     data(7)  = qclusts(3,i)
                     data(8)  = qclusts(8,i)
                     data(9)  = 1.0d0
                     data(10) = iblz1
                     data(11) = abs(qclusts(9,i))

                     call wrnt12(data)

               end if

*-----------------------------------------------------------------------
*              write the information of low energy photons
*              for MCNP on 13
*-----------------------------------------------------------------------

               if( ( ( knn .eq. 14 .and.
     &                 igcut .ge. 1 .and. igcut .le. 2  ) .or.
     &               ( knn .ge. 12 .and. knn .le. 14 .and.
     &                 igcut .eq. 3 ) ) .and.
     &                 qclusts(7,i) .gt. 0.0 .and.
     &                 qclusts(7,i) .le. emint ) then

                     if( knn .eq. 12 ) then

                        dii = 4.0d0
                        rtelen = rtelen + 1.0
                        wtelen = wtelen + qclusts(8,i)

                     else if( knn .eq. 13 ) then

                        dii = 5.0d0
                        rtelep = rtelep + 1.0
                        wtelep = wtelep + qclusts(8,i)

                     else if( knn .eq. 14 ) then

                        dii = 3.0d0
                        rtgamm = rtgamm + 1.0
                        wtgamm = wtgamm + qclusts(8,i)

                     end if

                     jclusts(4,i) = -2

                     data(1)  = qclusts(7,i)
                     data(2)  = qclusts(10,i)
                     data(3)  = qclusts(11,i)
                     data(4)  = qclusts(12,i)
                     data(5)  = qclusts(1,i)
                     data(6)  = qclusts(2,i)
                     data(7)  = qclusts(3,i)
                     data(8)  = qclusts(8,i)
                     data(9)  = dii
                     data(10) = iblz1
                     data(11) = abs(qclusts(9,i))

                     call wrnt13(data)

               end if

*-----------------------------------------------------------------------
*              write the information of low energy proton
*              for MCNP on 10
*-----------------------------------------------------------------------

               if( knn .eq. 1 .and. ipcut .ne. 0 .and.
     &             qclusts(7,i) .gt. 0.0 .and.
     &             qclusts(7,i) .le. emint ) then

                     jclusts(4,i) = -2

                     rtprot = rtprot + 1.0
                     wtprot = wtprot + qclusts(8,i)

                     data(1)  = qclusts(7,i)
                     data(2)  = qclusts(10,i)
                     data(3)  = qclusts(11,i)
                     data(4)  = qclusts(12,i)
                     data(5)  = qclusts(1,i)
                     data(6)  = qclusts(2,i)
                     data(7)  = qclusts(3,i)
                     data(8)  = qclusts(8,i)
                     data(9)  = 2.0d0
                     data(10) = iblz1
                     data(11) = abs(qclusts(9,i))

                     call wrnt10(data)

               end if
!$OMP END CRITICAL (ncut_crit)
*-----------------------------------------------------------------------

            end if

*-----------------------------------------------------------------------

         end do

*-----------------------------------------------------------------------
*           booking of the produced particles and fission (i=20)
*-----------------------------------------------------------------------
!$OMP CRITICAL
               do i = 1, 18
                  bprodp(i) = bprodp(i) + numpal(i)
                  aprodp(i) = aprodp(i) + rumpal(i)
               end do
                  bprodp(19) = bprodp(19) + numpal(0)
     &                       - numpal(15) - numpal(16)
     &                       - numpal(17) - numpal(18)
                  aprodp(19) = aprodp(19) + rumpal(0)
     &                       - rumpal(15) - rumpal(16)
     &                       - rumpal(17) - rumpal(18)
                  bprodp(20) = bprodp(20) + numpal(20)
                  aprodp(20) = aprodp(20) + rumpal(20)
!$OMP END CRITICAL
*-----------------------------------------------------------------------

         end if

*-----------------------------------------------------------------------
*     ncol = 0  : end of batch for deposit of heat tally
*     ncol = 101: for detector resolution
*-----------------------------------------------------------------------

      else if( ncol .eq. 0 .or. ncol .eq. 101 ) then

         call talls01(ncol)
         return

      end if

*-----------------------------------------------------------------------
*     ncol >= 10 :  call tally subroutine
*-----------------------------------------------------------------------

      if( ncol .ge. 10 ) then

                  call cputime(3)

               call talls01(ncol)
                  rncnt(3) = rncnt(3) + 1.0
                  call cputime(3)

      end if

*-----------------------------------------------------------------------
*     ncol = 4, 10 :  timer for source and crossing
*-----------------------------------------------------------------------

      if( ( ncol .eq. 10 .or. ncol .eq. 4 ) .and.
     &      ntmrg .gt. 0 ) then

cABE 2024/11/27 ... check whether cells where particle are located before/after transport are the same
            ichkcell = 1
            if( iblz1 .ne. iblz2 ) ichkcell = 0
            if( ichkcell .eq. 1 .and.
     &          ilev1 .gt. 0 .and. ilev2 .gt. 0 .and.
     &          ilev1 .eq. ilev2 ) then
               do kkk = 1, ilev1
                  if( ichkcell .eq. 0 ) exit
                  do jjj = 1, 5
                     if( ilat1(jjj,kkk) .ne. ilat2(jjj,kkk) ) then
                        ichkcell = 0
                        exit
                     endif
                  enddo
               enddo
            else
               ichkcell = 0
            endif
cABE end
                     idsm = intmc
                     jdsm = 0

                     kdsm = ktime
                     ldsm = 0

               do m = 1, ntmrg

                     jdsm = jdsm + 1
                     ntrn = idas_intmc(idsm+jdsm)
                     jdsm = jdsm + 1
                     mtrn = idas_intmc(idsm+jdsm)

                     ldsm  = ldsm + 1
                     inin  = idas_ktime(kdsm+ldsm)
                     ldsm  = ldsm + 1
                     inout = idas_ktime(kdsm+ldsm)
                     ldsm  = ldsm + 1
                     incol = idas_ktime(kdsm+ldsm)
                     ldsm  = ldsm + 1
                     inref = idas_ktime(kdsm+ldsm)

                  if( ( ( inout .ne. 0 .and. mark .eq. 0 .and.
     &                    ichkcell .eq. 0) .or.
     &                  ( inref .ne. 0 .and. mark .eq. 2 ) ) .and.
     &                    ncol .eq. 10 ) then   ! S.Abe avoid misoperation

                        jj  = 0

                     do ii = 1, ntrn

                        call tregck(iblz1,ilev1,ilat1,
     &                              mtrn,idas_intmc(idsm+jdsm+1),jj,icc)

                        if( icc .ne. 0 ) then

                           if( inout .ne. 0 .and.
     &                              mark .eq. 0 ) then

                              if( inout .eq. 1 ) then

                                 tc(ibktc+no,ipomp+1) =
     &                                    -abs(tc(ibktc+no,ipomp+1))

                              else if( inout .eq. -1 ) then

                                 tc(ibktc+no,ipomp+1) = 0.d0

                              end if

                           else if( inref .ne. 0 .and.
     &                              mark .eq. 2 ) then

                              if( inref .eq. 1 ) then

                                 tc(ibktc+no,ipomp+1) =
     &                               -abs(tc(ibktc+no,ipomp+1))

                              else if( inref .eq. -1 ) then

                                 tc(ibktc+no,ipomp+1) = 0.d0

                              end if

                           end if

                        end if

                     end do

                  end if

                  if( inin .ne. 0 .and.
     &                ( ncol .eq. 4 .or.
     &                  ( ncol .eq. 10 .and. ichkcell .eq. 0 ) )) then   ! S.Abe avoid misoperation

                        jj  = 0

                     do ii = 1, ntrn

                        call tregck(iblz2,ilev2,ilat2,
     &                              mtrn,idas_intmc(idsm+jdsm+1),jj,icc)

                        if( icc .ne. 0 ) then

                              if( inin .eq. 1 ) then

                                 tc(ibktc+no,ipomp+1) =
     &                                    -abs(tc(ibktc+no,ipomp+1))

                              else if( inin .eq. -1 ) then

                                 tc(ibktc+no,ipomp+1) = 0.d0

                              end if

                        end if

                     end do

                  end if

                     jdsm = jdsm + mtrn

               end do

      end if

*-----------------------------------------------------------------------
*     ncol = 4, 10 :  counter for source and crossing
*-----------------------------------------------------------------------

      if( ( ncol .eq. 10 .or. ncol .eq. 4 ) .and.
     &    ( ncntc(1) .eq. 1 .or. ncntc(2) .eq. 1 .or.
     &      ncntc(3) .eq. 1 ) ) then

cABE 2024/11/27 ... check whether cells where particle are located before/after transport are the same
            ichkcell = 1
            if( iblz1 .ne. iblz2 ) ichkcell = 0
            if( ichkcell .eq. 1 .and.
     &          ilev1 .gt. 0 .and. ilev2 .gt. 0 .and.
     &          ilev1 .eq. ilev2 ) then
               do kkk = 1, ilev1
                  if( ichkcell .eq. 0 ) exit
                  do jjj = 1, 5
                     if( ilat1(jjj,kkk) .ne. ilat2(jjj,kkk) ) then
                        ichkcell = 0
                        exit
                     endif
                  enddo
               enddo
            else
               ichkcell = 0
            endif
cABE end

            do k = 1, 3

            if( ncntc(k) .eq. 1 ) then
                call pcchck(k,ityp,ktyp,jtyp,icpan,icpat,icc)
            if( icc .eq. 1 ) then

                     idsm = incrt(k)
                     jdsm = 0

                     kdsm = kcont(k)
                     ldsm = 0

               do m = 1, ncreg(k)

                     jdsm = jdsm + 1
                     ntrn = idas_incrt(idsm+jdsm)
                     jdsm = jdsm + 1
                     mtrn = idas_incrt(idsm+jdsm)

                     ldsm  = ldsm + 1
                     inin  = idas_kcont(kdsm+ldsm)
                     ldsm  = ldsm + 1
                     inout = idas_kcont(kdsm+ldsm)
                     ldsm  = ldsm + 1
                     incol = idas_kcont(kdsm+ldsm)
                     ldsm  = ldsm + 1
                     inref = idas_kcont(kdsm+ldsm)
                     ldsm  = ldsm + 1
                     infis = idas_kcont(kdsm+ldsm)
                     ldsm   = ldsm + 1
                     inels  = idas_kcont(kdsm+ldsm)
                     ldsm   = ldsm + 1
                     iniel  = idas_kcont(kdsm+ldsm)
                     ldsm   = ldsm + 1
                     inncr  = idas_kcont(kdsm+ldsm)
                     ldsm   = ldsm + 1
                     indcy  = idas_kcont(kdsm+ldsm)
                     ldsm   = ldsm + 1
                     inato  = idas_kcont(kdsm+ldsm)
                     ldsm   = ldsm + 1
                     indlr  = idas_kcont(kdsm+ldsm)
                     ldsm   = ldsm + 1
                     influ  = idas_kcont(kdsm+ldsm)
                     ldsm   = ldsm + 1
                     inaug  = idas_kcont(kdsm+ldsm)
                     ldsm   = ldsm + 1
                     inbrm  = idas_kcont(kdsm+ldsm)
                     ldsm   = ldsm + 1
                     inphe  = idas_kcont(kdsm+ldsm)
                     ldsm   = ldsm + 1
                     incmp  = idas_kcont(kdsm+ldsm)
                     ldsm   = ldsm + 1
                     inppd  = idas_kcont(kdsm+ldsm)
                     ldsm   = ldsm + 1
                     inanh  = idas_kcont(kdsm+ldsm)
                     ldsm   = ldsm + 1
                     inmst  = idas_kcont(kdsm+ldsm)
                     ldsm   = ldsm + 1
                     inray  = idas_kcont(kdsm+ldsm)
                     ldsm   = ldsm + 1
                     inkoe  = idas_kcont(kdsm+ldsm)
                     ldsm   = ldsm + 1
                     indat  = idas_kcont(kdsm+ldsm)


                  if( ( ( inout .ne. 0 .and. mark .eq. 0 ) .or.
     &                  ( inref .ne. 0 .and. mark .eq. 2 ) ) .and.
     &                    ncol .eq. 10 ) then

                        jj  = 0

                     do ii = 1, ntrn

                        call tregck(iblz1,ilev1,ilat1,
     &                              mtrn,idas_incrt(idsm+jdsm+1),jj,icc)

                        if( icc .ne. 0 ) then

                           if( inout .ne. 0 .and.
     &                              mark .eq. 0 .and.
     &                         ichkcell .eq. 0 ) then   ! S.Abe 2024/11/27 avoid double counting

                              if( inout .eq. 10000 ) then

                                 ncnt(ibknct+k,no,ipomp+1) = 0

                              else if( abs( ncnt(ibknct+k,no,ipomp+1)
     &                                 + inout ) .le. 9999 ) then

                                 ncnt(ibknct+k,no,ipomp+1) =
     &                           ncnt(ibknct+k,no,ipomp+1) + inout
                                 if ( ncnt(ibknct+k,no,ipomp+1)
     &                               .gt.ncntmx(k) )
     &                           ncntmx(k) = ncnt(ibknct+k,no,ipomp+1)

                              end if

                           else if( inref .ne. 0 .and.
     &                              mark .eq. 2 ) then

                              if( inref .eq. 10000 ) then

                                 ncnt(ibknct+k,no,ipomp+1) = 0

                              else if( abs( ncnt(ibknct+k,no,ipomp+1)
     &                                 + inref ) .le. 9999 ) then

                                 ncnt(ibknct+k,no,ipomp+1) =
     &                           ncnt(ibknct+k,no,ipomp+1) + inref
                                 if ( ncnt(ibknct+k,no,ipomp+1)
     &                               .gt.ncntmx(k) )
     &                           ncntmx(k) = ncnt(ibknct+k,no,ipomp+1)

                              end if

                           end if

                        end if

                     end do

                  end if

                  if( inin .ne. 0 ) then

                        jj  = 0

                     do ii = 1, ntrn

                        call tregck(iblz2,ilev2,ilat2,
     &                              mtrn,idas_incrt(idsm+jdsm+1),jj,icc)

                        if( icc .ne. 0 ) then

                           if( ncol. eq. 4 .or.
     &                         ( ncol. eq. 10 .and.
     &                           ichkcell .eq. 0 ) ) then   ! S.Abe 2024/11/27 avoid double counting

                              if( inin .eq. 10000 ) then

                                 ncnt(ibknct+k,no,ipomp+1) = 0

                              else if( abs( ncnt(ibknct+k,no,ipomp+1)
     &                                 + inin ) .le. 9999 ) then

                                 ncnt(ibknct+k,no,ipomp+1) =
     &                           ncnt(ibknct+k,no,ipomp+1) + inin
                                 if ( ncnt(ibknct+k,no,ipomp+1)
     &                               .gt.ncntmx(k) )
     &                           ncntmx(k) = ncnt(ibknct+k,no,ipomp+1)

                              end if

                           end if

                        end if

                     end do

                  end if

                     jdsm = jdsm + mtrn

               end do

            end if
            end if
            end do

      end if

*-----------------------------------------------------------------------
*     ncol = 3 : end of a batch
*-----------------------------------------------------------------------

      if( ncol .eq. 3 ) then

               itstop = 1
               instop = 0
               idstop = 0

*-----------------------------------------------------------------------
*        single case, write batch status
*-----------------------------------------------------------------------

         if( npe .le. 1 ) then

*-----------------------------------------------------------------------

            if( infout .eq. 2 .or. infout .eq. 4 .or.
     &          infout .eq. 6 .or. infout .eq. 7 .or.
     &          infout .eq. 8 .or. timeout .gt. 0.d0 ) then

*-----------------------------------------------------------------------

              if ( nrandgen .eq. 0 ) then ! when LCG (2021.4.20)
               write(io,'(/79("-")/
     &                    "bat[",i8,"] ncas =",f16.0)')
     &               nobch, rcasc
               write(io,'(" rseed = ",1p1e25.16e3)')
     &               srijk
              else ! when xorshift
               write(io,'(/79("-")/
     &                    "bat[",i8,"] ncas =",f16.0)')
     &               nobch, rcasc
               write(io,'(" bitrseed = ",b64.64)')
     &               srijk
              end if


               if( incut .ne. 0 )
     &         write(io,'(''              ncut ='',f16.0,
     &                    '' :  ncut/s = '',1pe16.9)')
     &             rtneut, rtneut / rcasc

               if( igcut .ne. 0 )
     &         write(io,'(''              gcut ='',f16.0,
     &                    '' :  gcut/s = '',1pe16.9)')
     &             rtgamm, rtgamm / rcasc

               if( ipcut .ne. 0 )
     &         write(io,'(''              pcut ='',f16.0,
     &                    '' :  pcut/s = '',1pe16.9)')
     &             rtprot, rtprot / rcasc

               if( imbnk .gt. 0 .or. jmbnk .gt. 0 )
     &         write(io,'(''     tmp bnk acces ='',f16.0,
     &                    '' : f acces = '',f16.0)')
     &             dble(jmbnk), dble(imbnk)

                  call timex(time3,time2,0)

               if( timeout .gt. 0.0d0 ) then

                  if( time2 + time3 .gt. timeout ) itstop = 0

               end if

            else

                  call timex(time3,time2,1)

            end if

*-----------------------------------------------------------------------
*           batch status on file 'batch.out'
*-----------------------------------------------------------------------

               iot = 27

               open(iot,file=chfn(22),status='unknown')
               read(iot,*) nrmnbch
               close(iot)

               if ( nrmnbch .gt. nrmnbch0 ) nrmnbch = nrmnbch0
               nrmnbch = nrmnbch - 1
               nrmnbch0 = nrmnbch0 - 1
               istop = 1
               if ( nrmnbch.ne.nrmnbch0 .and. nrmnbch.le.0 ) then
                  istop = 0
                  nrmnbch = 0
               end if

               open(iot,file=chfn(22),status='unknown')

               write(iot,'(i0,a)')
     &              nrmnbch,' <--- number of remaining batches ' !S.H.(2019.3.14)

              if ( nrandgen .eq. 0 ) then ! when LCG (2021.4.20)
               write(iot,'(/79("-")/
     &                    "bat[",i8,"] ncas =",f16.0)')
     &               nobch, rcasc
               write(iot,'(" rseed = ",1p1e25.16e3)')
     &               srijk
              else ! when xorshift
               write(iot,'(/79("-")/
     &                    "bat[",i8,"] ncas =",f16.0)')
     &               nobch, rcasc
               write(iot,'(" bitrseed = ",b64.64)')
     &               srijk
              end if

               if( incut .ne. 0 )
     &         write(iot,'(''              ncut ='',f16.0,
     &                     '' :  ncut/s = '',1pe16.9)')
     &             rtneut, rtneut / rcasc

               if( igcut .ne. 0 )
     &         write(iot,'(''              gcut ='',f16.0,
     &                     '' :  gcut/s = '',1pe16.9)')
     &             rtgamm, rtgamm / rcasc

               if( ipcut .ne. 0 )
     &         write(iot,'(''              pcut ='',f16.0,
     &                     '' :  pcut/s = '',1pe16.9)')
     &             rtprot, rtprot / rcasc

                  ih = int( time3 / 3600. )
                  rm = time3 - dble(ih) * 3600.
                  im = int( rm / 60. )
                  ts = rm - dble(im) * 60.

               if( ih .ge. 1 ) then

                  write(iot,'(''          cpu time = '',
     &                       i3,'' h.'',i3,'' m.'',f6.2,'' s.'')')
     &                       ih, im, ts

               else if( im .ge. 1 ) then

                  write(iot,'(''          cpu time = '',
     &                       i3,'' m.'',f6.2,'' s.'')')
     &                       im, ts

               else

                  write(iot,'(''          cpu time = '',
     &                       f7.3,'' s.'')')
     &                       ts

               end if

               call date_a_time(iyer0,imon0,iday0,
     &                          ihor0,imin0,isec0)

               write(iot,'(/'' date = '',
     &                            i4,''-'',i2.2,''-'',i2.2)')
     &                            iyer0,imon0,iday0
               write(iot,'( '' time = '',
     &                            i2.2,''h '',i2.2,''m '',i2.2,''s''/)')
     &                            ihor0,imin0,isec0

               write(*,'(''bat['',i8,''] ncas ='',f16.0,
     &                   '' : date = '',i4,''-'',i2.2,''-'',i2.2,
     &              '' : time = '',i2.2,''h '',i2.2,''m '',i2.2,''s'')')
     &               nobch, rcasc,
     &               iyer0,imon0,iday0,
     &               ihor0,imin0,isec0

               if ( nrandgen .eq. 0 ) then ! when LCG
                  write(iot,'(79("-")/,
     &         "next initial random seed:",/1x,"rseed = ",1p1e25.16e3)')
     &            rijklst
               else ! when xorshift
                  write(iot,'(79("-")/,
     &         "next initial random seed:",/1x,"bitrseed = ",b64.64)')
     &            rijklst
               end if

               close(iot)

*-----------------------------------------------------------------------

            if( istop .eq. 0 .or. itstop .eq. 0 ) then

! T.Sato 2019/8/16, iot=27 is not opened at this moment

                  close(iot)

                  if( istop .eq. 0 ) then
                     write(*,'(
     &               ''**** calculation is stopped by batch.out !!'')')
                  else
                     write(*,'(
     &               ''**** calculation is stopped by timeout !!'')')
                  end if

                  ncol = -10
                  return

            end if

         end if

*-----------------------------------------------------------------------
*        parallel case
*        send parallel status and batch information to root
*-----------------------------------------------------------------------

         if( me .gt. 0 ) then

                  call timex(time3,time2,1)

                  icc = 0

                  bend(1) = srijk
                  bend(2) = reutn
                  bend(3) = dble(ncall)/dble(maxcas) !FURUTA
                  bend(4) = time3

                  bend(10)= time2

                  bend(5) = rtneut
                  bend(6) = rtgamm
                  bend(7) = rtprot
                  bend(8) = dble(imbnk)
                  bend(9) = dble(jmbnk)

                  call parasi(icc,1,0)

                  call parasr(bend,10,0)

                  call parasr(rsouin,1,0)
                  call parasi(maxbch,1,0)
                  if(itall.ne.-1)call parasr(rijklst,1,0) ! T.Sato 2017/08/05

                  call parari(istop,1,0)
                  call parari(itstop,1,0)

                  call date_a_time(iyer0,imon0,iday0,
     &                             ihor0,imin0,isec0)

                  write(*,'(''bat['',i6,''] ncas='',f12.0,
     &              '' : date= '',i4,''-'',i2.2,''-'',i2.2,
     &              '' : time= '',i2.2,''h '',i2.2,''m '',i2.2,''s'',
     &              '': pe= '',i5.0)')
     &              nobch, rcasc,
     &              iyer0,imon0,iday0,
     &              ihor0,imin0,isec0,me

*-----------------------------------------------------------------------

                  if( istop .eq. 0 .or. itstop .eq. 0 ) then

                     ncol = -10
                     return

                  end if

*-----------------------------------------------------------------------

         end if

*-----------------------------------------------------------------------
*        reset of the variables at the end of batch
*        tally summary at the end of each batch
*-----------------------------------------------------------------------

         if( me .gt. 0 .or. npe .le. 1 ) then

               ncall = 0
               reutn = 0.0

            if( itall .ne. -1 .and. nobch .lt. maxbch ) then ! T.Sato 2017/08/05

                  if(npe.le.1.and.idmpmode.eq.1)rsouin=rsouinbch*nobch !FURUTA20150515

                  ncol = 3

                  call cputime(3)

                  call talls01(ncol)
                  rncnt(3) = rncnt(3) + 1.0
                  call cputime(3)

*-----------------------------------------------------------------------
               if( me .gt. 0 ) then
                  call parari(itstd,itlmax,0)
                  call parari(instop,1,0)
                  call parari(idstop,1,0)
               end if

*-----------------------------------------------------------------------

               if( instop .gt. 0 .and. idstop .eq. 0 ) then

                  if( npe .le. 1 ) then

                     iot = 27
                     open(iot,file=chfn(22),status='unknown',
     &                    access='append')

                        write(iot,'( 79(''-'')
     &                  /'' calculation is stopped by stdcut''
     &                                 /79(''-''))')
                     close(iot)

                     write(*,'(
     &               ''**** calculation is stopped by stdcut !!'')')

                  end if

                     ncol = -10
                     return

               end if

*-----------------------------------------------------------------------

            end if

         end if

*-----------------------------------------------------------------------
*        receive the batch data from each PE
*-----------------------------------------------------------------------

         if( me .eq. 0 .and. npe .gt. 1 ) then

                     stim4 = paratim()
                     nobch = 0

            do i = 1, nmbch0

                     instop = 0
                     idstop = 0
*-----------------------------------------------------------------------

                     cpusm  = 0.0
                     cpual  = 0.0

                     ipusm  = 0
                     totevt = 0
                     maxbch = 0
                     sumneu = 0.0
                     sumgam = 0.0
                     sumpro = 0.0
                     smibnk = 0.0
                     smjbnk = 0.0

                     call paraiccp

               do j = 1, npe - 1

                  if( iccp(j) .eq. 0 ) then

                     call pararr(bend,10,j)

                     call pararr(rsouir,1,j)
                     call parari(maxbcr,1,j)

                     do k = 1, 4

                        bends(k,nbnd(i)+j) = bend(k) !FURUTA20220106

                     end do

                     cpusm  = cpusm  + bend(4)
                     cpual  = cpual  + bend(10)

                     sumneu = sumneu + bend(5)
                     sumgam = sumgam + bend(6)
                     sumpro = sumpro + bend(7)
                     if( smibnk .lt. bend(8) ) smibnk = bend(8)
                     if( smjbnk .lt. bend(9) ) smjbnk = bend(9)
                     ipusm  = ipusm  + 1
                     nobch  = nobch  + 1
                     totevt = totevt + rsouir
                     maxbch = maxbch + maxbcr

                  end if

               end do

               if(itall.ne.-1)then  ! T.Sato 2017/08/05
                 if(irndmode.eq.0)then
                   do j=1,npe-1
                     if( iccp(j) .eq. 0)then
                       call pararr(rijklst,1,j)
                     endif
                   enddo
                 else
                   do j=npe-1,1,-1
                     if( iccp(j) .eq. 0)then
                       call pararr(rijklst,1,j)
                     endif
                   enddo
                 endif
               endif

               if(idmpmode.eq.0)then
                     rsouin = totevt
               else
                     rsouin=rsouinbch*nobch
               endif

                     stim6 = paratim()
                     stim5 = stim6 - stim4
                     stim4 = stim6

                     avcpu = cpusm / dble( ipusm )

                     totcpu = cpual + cpusm

               if( timeout .gt. 0.0d0 ) then

                  if( totcpu .gt. timeout ) itstop = 0

               end if

*-----------------------------------------------------------------------

               iot = 27

               open(iot,file=chfn(22),status='unknown')
               read(iot,*) nrmnbch
               close(iot)

               if ( nrmnbch .gt. nrmnbch0 ) nrmnbch = nrmnbch0
               nrmnbch = nrmnbch - 1
               nrmnbch0 = nrmnbch0 - 1
               istop = 1
               if ( nrmnbch.ne.nrmnbch0 .and. nrmnbch.le.0 ) then
                  istop = 0
                  nrmnbch = 0
               end if

               do j = 1, npe - 1

                  if( iccp(j) .eq. 0 ) call parasi(istop,1,j)
                  if( iccp(j) .eq. 0 ) call parasi(itstop,1,j)

               end do

               open(iot,file=chfn(22),status='unknown')

               write(iot,'(i0,a)')
     &              nrmnbch,' <--- number of remaining batches ' !S.H.(2019.3.14)

               write(iot,'(/79(''-''))')
               write(iot,'('' Local Batch = ['',i6,'' ] :'',
     &                    '' Parallel Status: '',
     &                    '' 0-> normal, 1-> abnormal stop'')') i
               write(iot,'( 79(''-''))')

               write(iot,'('' ip  status'')')
               write(iot,'(2i5)') (k,iccp(k),k=1,npe-1)

               write(iot,'(/'' total source  ='',f16.0)')
     &                  rsouin

               if( incut .ne. 0 )
     &         write(iot,'(''              ncut ='',f16.0,
     &                     '' :  ncut/s = '',1pe16.9)')
     &                  sumneu, sumneu / rsouin

               if( igcut .ne. 0 )
     &         write(iot,'(''              gcut ='',f16.0,
     &                     '' :  gcut/s = '',1pe16.9)')
     &                  sumgam, sumgam / rsouin

               if( ipcut .ne. 0 )
     &         write(iot,'(''              pcut ='',f16.0,
     &                     '' :  pcut/s = '',1pe16.9)')
     &                  sumpro, sumpro / rcasc


               if( smibnk .gt. 0.0 .or. smjbnk .gt. 0.0 )
     &         write(io,'(  '' tmp bnk acces ='',f16.0,
     &                      '' : f acces = '',f16.0)')
     &                  smjbnk, smibnk

               write(iot,'(/'' elapse time ='',f13.2,'' sec.''/
     &                      '' av cpu time ='',f13.2,'' sec.'')')
     &                    stim5, avcpu


               call date_a_time(iyer0,imon0,iday0,
     &                          ihor0,imin0,isec0)

               write(iot,'(/'' date = '',
     &                            i4,''-'',i2.2,''-'',i2.2)')
     &                            iyer0,imon0,iday0
               write(iot,'( '' time = '',
     &                            i2.2,''h '',i2.2,''m '',i2.2,''s'')')
     &                            ihor0,imin0,isec0

               write(iot,'(/'' Starting Random Number of each PE''/
     &                      ''  ip num      rijk'')')

               do j = 1, npe - 1

                 if ( nrandgen .eq. 0 ) then ! when LCG (2021.4.20)
                  write(iot,'(1x,i5,3x,1p1e25.16e3)') j,
     &                           bends(1,nbnd(i)+j) !FURUTA20220106
                 else ! when xorshift
                  write(iot,'(1x,i5,3x,b64.64)') j,
     &                           bends(1,nbnd(i)+j) !FURUTA20220106
                 end if

               end do

               if ( nrandgen .eq. 0 ) then ! when LCG
                  write(iot,'(79("-")/,
     &         "next initial random seed:",/1x,"rseed = ",1p1e25.16e3)')
     &            rijklst
               else ! when xorshift
                  write(iot,'(79("-")/,
     &         "next initial random seed:",/1x,"bitrseed = ",b64.64)')
     &            rijklst
               end if

*-----------------------------------------------------------------------

               if( istop .eq. 0 .or. itstop .eq. 0 ) then

                  if( istop .eq. 0 ) then
                     write(iot,'( 79(''-'')
     &               /'' calculation is stopped by batch.out''
     &                           /79(''-''))')
                  else
                     write(iot,'( 79(''-'')
     &               /'' calculation is stopped by timeout''
     &                           /79(''-''))')
                  end if

                  if( istop .eq. 0 ) then
                     write(*,'(
     &               ''**** calculation is stopped by batch.out !!'')')
                  else
                     write(*,'(
     &               ''**** calculation is stopped by timeout !!'')')
                  end if

                  nmbch0 = i
                  ncol = -10
                  return

               end if

                  close(iot)

*-----------------------------------------------------------------------
*           tally summary at the end of each batch
*-----------------------------------------------------------------------

            if( itall .ne. -1 .and. nobch .lt. maxbch ) then ! T.Sato 2017/08/05

                  ncol = 3

                  call cputime(3)

                  call talls01(ncol)
                  rncnt(3) = rncnt(3) + 1.0
                  call cputime(3)

*-----------------------------------------------------------------------

               do j = 1, npe - 1
                  if( iccp(j) .eq. 0 ) call parasi(itstd,itlmax,j)
                  if( iccp(j) .eq. 0 ) call parasi(instop,1,j)
                  if( iccp(j) .eq. 0 ) call parasi(idstop,1,j)
               end do

               if( instop .gt. 0 .and. idstop .eq. 0 ) then

                  iot = 27
                  open(iot,file=chfn(22),status='unknown',
     &            access='append')

                  write(iot,'( 79(''-'')
     &               /'' calculation is stopped by stdcut''
     &                              /79(''-''))')
                  close(iot)

                  write(*,'(
     &            ''**** calculation is stopped by stdcut !!'')')

                  nmbch0 = i
                  ncol = -10
                  return

               end if

*-----------------------------------------------------------------------

            end if

*-----------------------------------------------------------------------

            end do

         end if

*-----------------------------------------------------------------------

               if( imbnk .gt. itbnk ) itbnk = imbnk
               if( jmbnk .gt. jtbnk ) jtbnk = jmbnk

               imbnk = 0
               jmbnk = 0

*-----------------------------------------------------------------------

      end if

*-----------------------------------------------------------------------
*     ncol = 2 : end of job
*                 close files
*-----------------------------------------------------------------------

      if( ncol .eq. 2 ) then

*-----------------------------------------------------------------------

         if( me .gt. 0 .or. npe .le. 1 ) then

            if( ipcut .ne. 0 ) then

               data(1) = -1.0d0
               call wrnt10(data)
               close(10)

            end if

            if( incut .ne. 0 ) then

               data(1) = -1.0d0
               call wrnt12(data)
               close(12)

            end if

            if( igcut .ne. 0 ) then

               data(1) = -1.0d0
               call wrnt13(data)
               close(13)

            end if

               do k = 1, imsrc

                  if( jstyp(k) .eq. 12 )  close(isorf(k))

               end do

         end if

            call sumout



*-----------------------------------------------------------------------
c        : creating dchain files at the end of job
c             MPI, openMP, Windows
            if ( me .eq. 0 ) then
         do m = 1, itnm
CCSE chg for mesh=r-z/xyz (2018.05.31) >>>>>
               if( itals(m) .eq. 38 .or.
     &             itals(m) .eq. 39 .or.
     &             itals(m) .eq. 40 .or.
     &             itals(m) .eq. 54 ) then ! add 39,40,54 for dchain mesh
                call pdchreg2(m)
               end if
CCSE chg for mesh=r-z/xyz (2018.05.31) <<<<<
         end do
            end if
c ----------------------------------------------------------------------

         if( itallo .ne. 0 ) then
               call DEALLOCATE_NDATA2
         end if

         if(nlat3.gt.0) call tetrafin !FURUTA20200505

         if(itdc.gt.0) call talldcfin !FURUTA20200522
*-----------------------------------------------------------------------
      end if
*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine sumout
*                                                                      *
*       summary of PHITS                                               *
*       modified by K.Niita on 2010/12/21                              *
*                                                                      *
************************************************************************
C for USE_MOD_COUNTER
      use mod_counter, only: rncnt,rnint,rnintr,rnpnt,rnpntr
     &                      ,aevts,aevtr,bevts,bevtr
     &                      ,iaevt, ibevt, jaevt, jbevt
C for REDUCTION_COUNTER
!$   &                      ,rncnt2,rnint2,rnintr2,rnpnt2,rnpntr2
!$   &                      ,aevts2,aevtr2,bevts2,bevtr2
      use MEMBANKMOD,only:mdbatima,ndbatima !FURUTA20160128
      use moddas_region
      use moddas_variance_reduction
      use moddas_bends, only: bends
      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'

*-----------------------------------------------------------------------

      common /mpi00/ npe, me
      common /mpi01/ iccp(20000)
      common /mpi02/ nmbch0

*-----------------------------------------------------------------------

      common /randsv/ srijk,rrijk !FURUTA
      common /randm5/ rijklst, rijkinit !OBINATA

*-----------------------------------------------------------------------

      common /kmat1a/ mxmat, mxmat0, mxnel
      common /inout/  in,io
      common /cparm/  maxbch,maxcas
      common /randn/ nrandgen ! S.H. xorshift (2020.2.6)
      common /iradkk/ randkk,irskip
      common /taliin/ rsouin, nzztin, nrgnin
      common /rcomon/ rcasc
      common /ngcut/  incut, igcut, ipcut

*-----------------------------------------------------------------------

      common /kounts/ rount(20)
      common /prodp/  aprodp(20), bprodp(20)
      common /dcayp/  adcayp(20), bdcayp(20)
      common /stopp/  astopp(20), bstopp(20)
      common /timep/  atimep(20), btimep(20)
      common /leakp/  aleakp(20), bleakp(20)

      common /otherp/ nothp(2000,2), nothn
      common /othjmp/ nojmp(2000,2), nojmn
      common /othstp/ nostp(2000,2), nostn
      common /othtim/ notip(2000,2), notin
      common /othlkp/ nolkp(2000,2), nolkn
      common /othdcp/ nodcp(2000,6), nodcn

      common /otheid/ idpat(20), idoth(200), idono
      common /ncutne/ wtneut,rtneut
      common /ncutgm/ wtgamm,rtgamm
      common /ncuten/ wtelen,rtelen
      common /ncutep/ wtelep,rtelep
      common /ncutpr/ wtprot,rtprot

      common /wtcntl/ iwt, icimp(20), ifcls(20), iwwin(20), ircls(20)

      common /isorfs/ ispfs(isrc), rspfn, rspfz, ispfn

      common /impmsg/ iimpn, isimp, iswct, isstr, maxip,
     &                mnimp(7,0:20),
     &                kfimp(7), inimc(7), inimt(7)

      common /fclmsg/ ifcln, isfcl, maxrg,
     &                mnfcl(6,0:20), mrfcl(kvlmax),
     &                kfcls(6), inflc(6), inflt(6)

      common /wwindn/ eenww(6,100), iwwdp, mnwwp(6,0:20), kfwwp(6),
     &                inwwc(6), ienww(6), inwwt(6), iswwp, maxww

      common /wwind0/ dvals, iwmsh, kecho, kgwwp(6), idval
      common /wwxyz/ iwxty(6), iwxnm(6), iwxrg(6),
     &               rwxmi(6), rwxma(6), rwxdl(6),
     &               iwyty(6), iwynm(6), iwyrg(6),
     &               rwymi(6), rwyma(6), rwydl(6),
     &               iwzty(6), iwznm(6), iwzrg(6),
     &               rwzmi(6), rwzma(6), rwzdl(6)

      common /parai/  ipsq(400)
      common /paraj/  mstz(300), parz(300) !OBINATA
      common /inpec/  ititl, ipara, ibody, iregn, llarr, itby, itar
      common /wparm/  swtm(20), wc1(20), wc2(20) !FURUTA
      common /kmat1d/ idmn(0:kvlmax), idnm(kvmmax)
      common /regdc/  idrg(kvlmax), idgr(kvmmax)
      common /cgerr/  nlost, ilost, igerr, icger, ncger, nrecover
      common /cgstr/  novp, nrovp(3,1000)
      common /talout/ itall
      common /ptname/ pname(20), ipln(20)
      character       pname*8
      common /kcomsm/ imbnk, jmbnk, itbnk, jtbnk !FURUTA

      common /sumbnd/ nbnd(nbchmax) !FURUTA20220106


      common /infprint/ infout

*-----------------------------------------------------------------------
      common /stat / istdev, irestart, ireschk
      common /res01/ istdevres,maxcasres,rijklstres,irdrf
      common /res02/ crdrfln, irdrfll
      character crdrfln*100

      common /res04/ lrijkeqrf
      logical lrijkeqrf

      dimension ios(2)

*-----------------------------------------------------------------------

      dimension noths(2000)
      character uname(20)*8

      character dum1*10000
      character dum2*10000

*-----------------------------------------------------------------------



*-----------------------------------------------------------------------
*     for parallel
*-----------------------------------------------------------------------

      dimension rounr(20)
      character aount(20)*35

      dimension aprodr(20), bprodr(20)
      dimension adcayr(20), bdcayr(20)
      dimension astopr(20), bstopr(20)
      dimension atimer(20), btimer(20)
      dimension aleakr(20), bleakr(20)

      dimension nothr(2000,2)
      dimension nojmr(2000,2)
      dimension nostr(2000,2)
      dimension notir(2000,2)
      dimension nolkr(2000,2)

      dimension nodcr(2000,6)

      dimension nrovr(3,1000)

*-----------------------------------------------------------------------

      character qname*16
      character qnamd(0:3)*16

*-----------------------------------------------------------------------
*                                                        ncol
      data (aount(i),i=1,15)/
     &            'start of calculation',              !  1
     &            'end of calculation',                !  2
     &            'end of a batch',                    !  3
     &            'source',                            !  4
     &            'detection of geometry error/warn',  !  5
     &            'recovery from geometry warning',    !  6
     &            'termination by geometry error',     !  7
     &            'termination by weight cut-off',     !  8
     &            'termination by time cut-off',       !  9
     &            'geometry boundary crossing',        ! 10
     &            'termination by energy cut-off',     ! 11
     &            'termination by escape or leakage',  ! 12
     &            'reaction : (n,x) type',             ! 13
     &            "reaction : (n,n'x) type",           ! 14
     &            "sequential transport for tally"/    ! 15

*-----------------------------------------------------------------------
      common /dmpinfo/ rsouinbch,maxbchdmp,maxcasdmp
      common /stat2/ dmpmulti,idmpmode,ibchjmp,idmpjmp(2)
*-----------------------------------------------------------------------
      real(8),allocatable :: aevt1(:),bevt1(:)
*-----------------------------------------------------------------------

      ios = (/ io, 6 /) !OBINATA(2012.7.26)

*-----------------------------------------------------------------------
*     write batch information on output file
*-----------------------------------------------------------------------

      if( infout .eq. 2 .or. infout .eq. 4 .or. infout .eq. 6 .or.
     &    infout .eq. 7 .or. infout .eq. 8 ) then

*-----------------------------------------------------------------------

      if( me .eq. 0 .and. npe .gt. 1 ) then

                     ibch = 0

         do i = 1, nmbch0

            do j = 1, npe - 1

               if( iccp(j) .eq. 0 ) then

                     ibch  = ibch + 1
                     rcasc = dble( maxcas ) * dble( ibch )

                 if ( nrandgen .eq. 0 ) then ! when LCG (2021.4.20)
                  write(io,'(/79("-")/
     &                       "bat[",i8,"] ncas =",f16.0,
     &                       " : pe[",i5," ]")')
     &                  ibch, rcasc, j
                  write(io,'(" rijk = ",1p1e25.16e3)')
     &                  bends(1,nbnd(i)+j) !FURUTA20220106
                 else ! when xorshift
                  write(io,'(/79("-")/
     &                       "bat[",i8,"] ncas =",f16.0,
     &                       " : pe[",i5," ]")')
     &                  ibch, rcasc, j
                  write(io,'(" bitrseed = ",b64.64)')
     &                  bends(1,nbnd(i)+j) !FURUTA20220106
                 end if


                  time3 = bends(4,nbnd(i)+j) !FURUTA20220106

                  ih = int( time3 / 3600. )
                  rm = time3 - dble(ih) * 3600.
                  im = int( rm / 60. )
                  ts = rm - dble(im) * 60.

                  if( ih .ge. 1 ) then

                     write(io,'(''          cpu time = '',
     &                          i3,'' h.'',i3,'' m.'',f6.2,'' s.'')')
     &                          ih, im, ts

                  else if( im .ge. 1 ) then

                     write(io,'(''          cpu time = '',
     &                          i3,'' m.'',f6.2,'' s.'')')
     &                          im, ts

                  else

                     write(io,'(''          cpu time = '',
     &                          f7.3,'' s.'')')
     &                          ts

                  end if

               end if

            end do

         end do

      end if

*-----------------------------------------------------------------------

      end if

*-----------------------------------------------------------------------
*     final check of parallel staus
*-----------------------------------------------------------------------

         if( npe .gt. 1 ) then

            if( me .gt. 0 ) then

               icc = 0

               call parasi(icc,1,0)

            else

               call paraiccp

            end if

         end if

*-----------------------------------------------------------------------

         if( me .eq. 0 ) then

               write(io,'(/79(''=''))')
               write(io,'(/'' Summary for the end of job'')')

         end if

*-----------------------------------------------------------------------
*     summary of analyz call
*-----------------------------------------------------------------------

         if( npe .gt. 1 ) then

            if( me .gt. 0 ) then

               call parasr(rount,20,0)

            else

                     do k = 1, 20

                        rount(k) = 0.0

                     end do

               do j = 1, npe - 1

                  if( iccp(j) .eq. 0 ) then

                     call pararr(rounr,20,j)

                     do k = 1, 20

                        rount(k) = rount(k) + rounr(k)

                     end do

                  end if

               end do

            end if

         end if

*-----------------------------------------------------------------------

         if( me .eq. 0 ) then

               write(io,'(/79(''-''))')

               write(io,'('' number of analyz call vs ncol''/
     &                     79(''-'')/
     &                    ''     ncol          number'')')

               write(io,'(i9,f16.0,''   : '',a32)')
     &                          (i,rount(i),aount(i), i=1,15)

         end if

*-----------------------------------------------------------------------
*     total number of source
*-----------------------------------------------------------------------

            totevt = rsouin

*-----------------------------------------------------------------------

         if( npe .gt. 1 ) then

            if( me .gt. 0 ) then

                     call parasr(rsouin,1,0)

            else

                     totevt = 0

               do i = 1, npe - 1

                  if( iccp(i) .eq. 0 ) then

                     call pararr(rsouir,1,i)

                     totevt = totevt + rsouir

                  end if

               end do

                     rsouin = totevt

            end if

         end if

*-----------------------------------------------------------------------
*        total number of batch
*-----------------------------------------------------------------------

         if( npe .gt. 1 ) then

            if( me .gt. 0 ) then

                     call parasi(maxbch,1,0)

            else

                     maxbch = 0

               do i = 1, npe - 1

                  if( iccp(i) .eq. 0 ) then

                     call parari(maxbcr,1,i)

                     maxbch = maxbch + maxbcr

                  end if

               end do

            end if

         end if
         if(idmpmode.eq.1)then
          totevt = dble(rsouinbch*maxbch)
          rsouin = totevt
         endif

*-----------------------------------------------------------------------
*        final random number
*-----------------------------------------------------------------------
         if( npe .gt. 1 ) then

           if( me .gt. 0 ) then
             call parasr(rijklst,1,0)
           else
               do i = 1, npe - 1
                 if( iccp(i) .eq. 0 ) then
                   call pararr(rijklst,1,i)
                 endif
               enddo
           end if

         end if
         rijk=rijklst
c........

*-----------------------------------------------------------------------
*        call tally summary subroutine
*-----------------------------------------------------------------------

            ncol = 2

                  call cputime(3)

                  call talls01(ncol)
                  rncnt(3) = rncnt(3) + 1.0
                  call cputime(3)

*-----------------------------------------------------------------------

      if( infout .eq. 8 ) then

*-----------------------------------------------------------------------
*     summary of cell important function
*-----------------------------------------------------------------------

      if( iwt .ne. 0 ) then

         if( npe .gt. 1 ) then

            do kk = 1, 20

                  iini = ( kk - 1 ) * maxip * 7 !FURUTA20130820

               if( me .gt. 0 ) then

                        call parasr
     &                     (wtximp_pointer(isimp+iini+1), maxip*7, 0)

               else

                        do j = 1, maxip * 7

                           wtximp_pointer(isimp+j+iini) = 0.0d0

                        end do

                  do i = 1, npe - 1

                     if( iccp(i) .eq. 0 ) then

                        call pararr(das_isstr(isstr+1), maxip*7, i)

                        do j = 1, maxip * 7

                           wtximp_pointer(isimp+j+iini)
     &                        = wtximp_pointer(isimp+j+iini)
     &                          + das_isstr(isstr+j)

                        end do

                     end if

                  end do

               end if

            end do

         end if

*-----------------------------------------------------------------------

      if( me .eq. 0 ) then

            do ip = 1, 19
            do id = 1, 7
            do ir = 1, maxip

               wtximp(ip+isimp,id,ir) = wtximp(ip+isimp,id,ir) / totevt

            end do
            end do
            end do

*-----------------------------------------------------------------------
*        Splitting
*-----------------------------------------------------------------------

                     igm = ( mmmax - 1 ) * 2 + 1

         do kk = 1, iimpn + 1

                     kimps = 0

         do k = 1, mnimp(kk,20)

                     lngmax = 0
                     jsek   = 0
                     t1     = 0.0d0

                     ii = 0

                     idsm = inimt(kk)
                     jdsm = 0

               do m = 1, mnimp(kk,0) + 1

                     jdsm = jdsm + 1
                     ntrn = idas_inimt(idsm+jdsm)
                     jdsm = jdsm + 1
                     mtrn = idas_inimt(idsm+jdsm)
                     kssm = idsm + jdsm + 1
                     jdsm = jdsm + mtrn

                     jj = 0

                  do j = 1, ntrn

                     ii = ii + 1

                     call echrg3(jj,mtrn,idas_inimt(kssm)
     &                           ,dum1,lng1,icmb,igm)

                     if( wtximp(mnimp(kk,k)+isimp,3,ii) .gt. 0.d0 ) then

                         jsek = jsek + 1

                         if( lng1 .gt. lngmax ) lngmax = lng1

                         t1 = t1 + wtximp(mnimp(kk,k)+isimp,3,ii)

                     end if

                  end do

               end do

*-----------------------------------------------------------------------

            if( jsek .ne. 0 ) then

               if( kimps .eq. 0 ) then

                  write(io,'(/79(''-''))')
                  write(io,'('' Summary of importance per source'',
     &            '' (non zero) for '',i2,''-th [importance]''/
     &            '' === Splitting ==='')') kk
                     write(io,'(79(''-''))')

                  kimps = kimps + 1

               end if

                     uname(k) = pname(mnimp(kk,k))(1:8)

*-----------------------------------------------------------------------

                     lblk = lngmax
                     dum2(1:5) = '  reg'

                   do i = 1, lblk
                     dum2(i+5:i+5) = ' '
                   end do

                     dum2(lblk+5+1:lblk+5+59) =
     &              '[ '// uname(k) //
     &              ']    incoming    splitted    splitted'//
     &              '   splitting'

                     write(io,'(/600a1)') (dum2(i:i),i=1,lblk+5+59)

                     dum2(1:5) = '     '
                   do i = 1, lblk
                     dum2(i+5:i+5) = ' '
                   end do

                     dum2(lblk+5+1:lblk+5+59) =
     &              ' importance      weight      weight      number'//
     &              '      events'

                     write(io,'(600a1)') (dum2(i:i),i=1,lblk+5+59)

                     ii = 0

                     idsm = inimt(kk)
                     jdsm = 0

                     kdsm = kfimp(kk)

               do m = 1, mnimp(kk,0) + 1

                     jdsm = jdsm + 1
                     ntrn = idas_inimt(idsm+jdsm)
                     jdsm = jdsm + 1
                     mtrn = idas_inimt(idsm+jdsm)
                     kssm = idsm + jdsm + 1
                     jdsm = jdsm + mtrn

                     rimp = das_kfimp(kdsm-1+m)

                     jj = 0

                  do j = 1, ntrn

                     ii = ii + 1

                     call echrg3(jj,mtrn,idas_inimt(kssm)
     &                           ,dum1,lng1,icmb,igm)

                  if( wtximp(mnimp(kk,k)+isimp,3,ii) .gt. 0.0 ) then

                     dum2(1:3) = '   '

                     do i = 1, lng1
                        dum2(i+3:i+3) = dum1(i:i)
                     end do

                     lngb = lngmax - lng1

                     do i = 1, lngb
                        dum2(i+lng1+3:i+lng1+3) = ' '
                     end do

                     if( wtximp(mnimp(kk,k)+isimp,2,ii) .le. 0.d0 ) then
                        rn = 0.0
                     else
                        rn = wtximp(mnimp(kk,k)+isimp,1,ii)
     &                     / wtximp(mnimp(kk,k)+isimp,2,ii)
                     endif

                     write(dum2(lngb+lng1+4:lngb+lng1+65),
     &                         '(1x,1p5e12.4)') rimp,
     &                         wtximp(mnimp(kk,k)+isimp,1,ii),
     &                         wtximp(mnimp(kk,k)+isimp,2,ii),
     &                         rn,
     &                         wtximp(mnimp(kk,k)+isimp,3,ii)

                     write(io,'(600a1)') (dum2(i:i),i=1,lngb+lng1+65)

                  end if
                  end do

               end do

                     dum2(1:13) = ' total events'

                     lstl = 3+lngmax+49

                  do i = 14, lstl

                     dum2(i:i) = ' '

                  end do

                     write(dum2(lstl+1:lstl+13),
     &                         '(1p1e12.4)') t1

                     write(io,'(600a1)') (dum2(i:i),i=1,lstl+13)


            end if

         end do
         end do

*-----------------------------------------------------------------------
*        Russian Roulette
*-----------------------------------------------------------------------

                     igm = ( mmmax - 1 ) * 2 + 1

         do kk = 1, iimpn + 1

                     kimps = 0

         do k = 1, mnimp(kk,20)

                     lngmax = 0
                     jsek   = 0
                     t2     = 0.0d0
                     t3     = 0.0d0

                     ii = 0

                     idsm = inimt(kk)
                     jdsm = 0

               do m = 1, mnimp(kk,0) + 1

                     jdsm = jdsm + 1
                     ntrn = idas_inimt(idsm+jdsm)
                     jdsm = jdsm + 1
                     mtrn = idas_inimt(idsm+jdsm)
                     kssm = idsm + jdsm + 1
                     jdsm = jdsm + mtrn

                     jj = 0

                  do j = 1, ntrn

                     ii = ii + 1

                     call echrg3(jj,mtrn,idas_inimt(kssm)
     &                           ,dum1,lng1,icmb,igm)

                     if( wtximp(mnimp(kk,k)+isimp,5,ii) .gt. 0.d0 .or.
     &                   wtximp(mnimp(kk,k)+isimp,7,ii) .gt. 0.d0 ) then

                         jsek = jsek + 1

                         if( lng1 .gt. lngmax ) lngmax = lng1

                         t2 = t2 + wtximp(mnimp(kk,k)+isimp,5,ii)
                         t3 = t3 + wtximp(mnimp(kk,k)+isimp,7,ii)

                     end if

                  end do

               end do

*-----------------------------------------------------------------------

            if( jsek .ne. 0 ) then

               if( kimps .eq. 0 ) then

                  write(io,'(/79(''-''))')
                  write(io,'('' Summary of importance per source'',
     &            '' (non zero) for '',i2,''-th [importance]''/
     &            '' === Russian Roulette ==='')') kk
                     write(io,'(79(''-''))')

                  kimps = kimps + 1

               end if

                     uname(k) = pname(mnimp(kk,k))(1:8)

*-----------------------------------------------------------------------

                     lblk = lngmax
                     dum2(1:5) = '  reg'

                   do i = 1, lblk
                     dum2(i+5:i+5) = ' '
                   end do

                     dum2(lblk+5+1:lblk+5+59) =
     &              '[ '// uname(k) //
     &              ']    survived    survived      killed'//
     &              '      killed'

                     write(io,'(/600a1)') (dum2(i:i),i=1,lblk+5+59)

                     dum2(1:5) = '     '
                   do i = 1, lblk
                     dum2(i+5:i+5) = ' '
                   end do

                     dum2(lblk+5+1:lblk+5+59) =
     &              ' importance      weight      events      weight'//
     &              '      events'

                     write(io,'(600a1)') (dum2(i:i),i=1,lblk+5+59)

                     ii = 0

                     idsm = inimt(kk)
                     jdsm = 0

                     kdsm = kfimp(kk)

               do m = 1, mnimp(kk,0) + 1

                     jdsm = jdsm + 1
                     ntrn = idas_inimt(idsm+jdsm)
                     jdsm = jdsm + 1
                     mtrn = idas_inimt(idsm+jdsm)
                     kssm = idsm + jdsm + 1
                     jdsm = jdsm + mtrn

                     rimp = das_kfimp(kdsm-1+m)

                     jj = 0

                  do j = 1, ntrn

                     ii = ii + 1

                     call echrg3(jj,mtrn,idas_inimt(kssm)
     &                           ,dum1,lng1,icmb,igm)

                  if( wtximp(mnimp(kk,k)+isimp,5,ii) .gt. 0.0 .or.
     &                wtximp(mnimp(kk,k)+isimp,7,ii) .gt. 0.0 ) then

                     dum2(1:3) = '   '

                     do i = 1, lng1
                        dum2(i+3:i+3) = dum1(i:i)
                     end do

                     lngb = lngmax - lng1

                     do i = 1, lngb
                        dum2(i+lng1+3:i+lng1+3) = ' '
                     end do

                     write(dum2(lngb+lng1+4:lngb+lng1+65),
     &                         '(1x,1p5e12.4)') rimp,
     &                         wtximp(mnimp(kk,k)+isimp,4,ii),
     &                         wtximp(mnimp(kk,k)+isimp,5,ii),
     &                         wtximp(mnimp(kk,k)+isimp,6,ii),
     &                         wtximp(mnimp(kk,k)+isimp,7,ii)

                     write(io,'(600a1)') (dum2(i:i),i=1,lngb+lng1+65)

                  end if
                  end do

               end do

                     dum2(1:13) = ' total events'

                     lstl = 3+lngmax+25

                  do i = 14, lstl

                     dum2(i:i) = ' '

                  end do

                     write(dum2(lstl+1:lstl+13),
     &                         '(1p1e12.4)') t2

                     lstl = 3+lngmax+49

                  do i = 3+lngmax+38, lstl

                     dum2(i:i) = ' '

                  end do

                     write(dum2(lstl+1:lstl+13),
     &                         '(1p1e12.4)') t3

                     write(io,'(600a1)') (dum2(i:i),i=1,lstl+13)

            end if

         end do
         end do

*-----------------------------------------------------------------------

      end if
      end if

*-----------------------------------------------------------------------
*     summary of weight window
*-----------------------------------------------------------------------

      if( iwwdp .ne. 0 ) then

         if( npe .gt. 1 ) then

            do kk = 1, 20

                  iini = ( kk - 1 ) * maxww * 14 !FURUTA20130820

               if( me .gt. 0 ) then

                        call parasr
     &                     (wtxwwp_pointer(iswwp+iini+1), maxww*14, 0)

               else

                        do j = 1, maxww * 14

                           wtxwwp_pointer(iswwp+j+iini) = 0.0d0

                        end do

                  do i = 1, npe - 1

                     if( iccp(i) .eq. 0 ) then

                        call pararr(das_isstr(isstr+1), maxww*14, i)

                        do j = 1, maxww * 14

                           wtxwwp_pointer(iswwp+j+iini)
     &                        = wtxwwp_pointer(iswwp+j+iini)
     &                          + das_isstr(isstr+j)

                        end do

                     end if

                  end do

               end if

            end do

         end if

*-----------------------------------------------------------------------

      if( me .eq. 0 ) then

            do ip = 1, 19
            do id = 1, 14
            do ir = 1, maxww

               wtxwwp(ip+iswwp,id,ir) = wtxwwp(ip+iswwp,id,ir) / totevt

            end do
            end do
            end do

*-----------------------------------------------------------------------
*        Splitting at Boundary
*-----------------------------------------------------------------------

                     igm = ( mmmax - 1 ) * 2 + 1

         do kk = 1, iwwdp

                     kimps = 0

         do k = 1, mnwwp(kk,20)

                     lngmax = 0
                     jsek   = 0
                     t1     = 0.0d0

*-----------------------------------------------------------------------

            if( iwmsh .eq. 1 ) then

                     ii = 0

                     idsm = inwwt(kk)
                     jdsm = 0

               do m = 1, mnwwp(kk,0)

                     jdsm = jdsm + 1
                     ntrn = idas_inwwt(idsm+jdsm)
                     jdsm = jdsm + 1
                     mtrn = idas_inwwt(idsm+jdsm)
                     kssm = idsm + jdsm + 1
                     jdsm = jdsm + mtrn

                     jj = 0

                  do j = 1, ntrn

                     ii = ii + 1

                     call echrg3(jj,mtrn,idas_inwwt(kssm)
     &                           ,dum1,lng1,icmb,igm)

                     if( wtxwwp(mnwwp(kk,k)+iswwp,3,ii) .gt. 0.d0 ) then

                         jsek = jsek + 1

                         if( lng1 .gt. lngmax ) lngmax = lng1

                         t1 = t1 + wtxwwp(mnwwp(kk,k)+iswwp,3,ii)

                     end if

                  end do

               end do

*-----------------------------------------------------------------------

            else if( iwmsh .eq. 3 ) then

                     inx = iwxnm(kk)
                     iny = iwynm(kk)
                     inz = iwznm(kk)
                     ixyz = inx * iny * inz

                  do jx = 1, inx
                  do jy = 1, iny
                  do jz = 1, inz

                     ii = jz + inz * ( jy - 1 )+ iny * inz * ( jx - 1 )

                     if( wtxwwp(mnwwp(kk,k)+iswwp,3,ii) .gt. 0.d0 ) then

                         jsek = jsek + 1

                         t1 = t1 + wtxwwp(mnwwp(kk,k)+iswwp,3,ii)

                     end if

                  end do
                  end do
                  end do

*-----------------------------------------------------------------------

            else if( iwmsh .eq. 4 ) then

                  do j = 1, mnwwp(kk,0)

                     if( wtxwwp(mnwwp(kk,k)+iswwp,3,j) .gt. 0.d0 ) then

                         jsek = jsek + 1

                         t1 = t1 + wtxwwp(mnwwp(kk,k)+iswwp,3,j)

                     end if

                  end do

            end if

*-----------------------------------------------------------------------

            if( jsek .ne. 0 ) then

               if( kimps .eq. 0 ) then

                  write(io,'(/79(''-''))')
                  write(io,'('' Summary of weight window per source'',
     &            '' (non zero) for '',i2,''-th [weight window]''/
     &            '' === Splitting at Boundary ==='')') kk
                     write(io,'(79(''-''))')

                  kimps = kimps + 1

               end if

                     uname(k) = pname(mnwwp(kk,k))(1:8)

*-----------------------------------------------------------------------

            if( iwmsh .eq. 1 ) then

                     lblk = lngmax
                     dum2(1:5) = '  reg'

                   do i = 1, lblk
                     dum2(i+5:i+5) = ' '
                   end do

                     lbin = 5 + lblk

*-----------------------------------------------------------------------

            else if( iwmsh .eq. 3 ) then

                     dum2(1:2) = ' ('
                     lni = 2 + 1

                  if( iwxnm(kk) .lt. 10 ) then
                     laix = 1
                     dum2(lni:lni) = 'x'
                  else if( iwxnm(kk) .lt. 100 ) then
                     laix = 2
                     dum2(lni:lni+laix-1) = ' x'
                  else if( iwxnm(kk) .lt. 1000 ) then
                     laix = 3
                     dum2(lni:lni+laix-1) = '  x'
                  else if( iwxnm(kk) .lt. 10000 ) then
                     laix = 4
                     dum2(lni:lni+laix-1) = '   x'
                  else if( iwxnm(kk) .lt. 100000 ) then
                     laix = 5
                     dum2(lni:lni+laix-1) = '    x'
                  end if

                     lni = lni + laix
                     dum2(lni:lni) = ' '
                     lni = lni + 1

                  if( iwynm(kk) .lt. 10 ) then
                     laiy = 1
                     dum2(lni:lni) = 'y'
                  else if( iwynm(kk) .lt. 100 ) then
                     laiy = 2
                     dum2(lni:lni+laiy-1) = ' y'
                  else if( iwynm(kk) .lt. 1000 ) then
                     laiy = 3
                     dum2(lni:lni+laiy-1) = '  y'
                  else if( iwynm(kk) .lt. 10000 ) then
                     laiy = 4
                     dum2(lni:lni+laiy-1) = '   y'
                  else if( iwynm(kk) .lt. 100000 ) then
                     laiy = 5
                     dum2(lni:lni+laix-1) = '    y'
                  end if

                     lni = lni + laiy
                     dum2(lni:lni) = ' '
                     lni = lni + 1

                  if( iwznm(kk) .lt. 10 ) then
                     laiz = 1
                     dum2(lni:lni) = 'z'
                  else if( iwznm(kk) .lt. 100 ) then
                     laiz = 2
                     dum2(lni:lni+laiz-1) = ' z'
                  else if( iwznm(kk) .lt. 1000 ) then
                     laiz = 3
                     dum2(lni:lni+laiz-1) = '  z'
                  else if( iwznm(kk) .lt. 10000 ) then
                     laiz = 4
                     dum2(lni:lni+laiz-1) = '   z'
                  else if( iwznm(kk) .lt. 100000 ) then
                     laiz = 5
                     dum2(lni:lni+laiz-1) = '    z'
                  end if

                     lni = lni + laiz
                     dum2(lni:lni+1) = ') '

                     lbin = lni+1

*-----------------------------------------------------------------------

            else if( iwmsh .eq. 4 ) then

                     lngmax = 8
                     lblk = lngmax
                     dum2(1:6) = '   tet'

                   do i = 1, lblk
                     dum2(i+6:i+6) = ' '
                   end do

                     lbin = lngmax + 6

            end if

*-----------------------------------------------------------------------

                     dum2(lbin+1:lbin+59) =
     &              '[ '// uname(k) //
     &              ']    incoming    splitted    splitted'//
     &              '   splitting'

                     write(io,'(/600a1)') (dum2(i:i),i=1,lbin+59)

                   do i = 1, lbin
                     dum2(i:i) = ' '
                   end do

                     dum2(lbin+1:lbin+59) =
     &              ' ww1 bound       weight      events      weight'//
     &              '      events'

                     write(io,'(600a1)') (dum2(i:i),i=1,lbin+59)

*-----------------------------------------------------------------------

            if( iwmsh .eq. 1 ) then

                     ii = 0

                     idsm = inwwt(kk)
                     jdsm = 0

                     kdsm = kfwwp(kk)

               do m = 1, mnwwp(kk,0)

                     jdsm = jdsm + 1
                     ntrn = idas_inwwt(idsm+jdsm)
                     jdsm = jdsm + 1
                     mtrn = idas_inwwt(idsm+jdsm)
                     kssm = idsm + jdsm + 1
                     jdsm = jdsm + mtrn

                     rimp = das_kfwwp(kdsm+m-1)

                     jj = 0

                  do j = 1, ntrn

                     ii = ii + 1

                     call echrg3(jj,mtrn,idas_inwwt(kssm)
     &                           ,dum1,lng1,icmb,igm)

                  if( wtxwwp(mnwwp(kk,k)+iswwp,3,ii) .gt. 0.0d0 ) then

                     dum2(1:3) = '   '

                     do i = 1, lng1
                        dum2(i+3:i+3) = dum1(i:i)
                     end do

                     lngb = lngmax - lng1

                     do i = 1, lngb
                        dum2(i+lng1+3:i+lng1+3) = ' '
                     end do

                     if( wtxwwp(mnwwp(kk,k)+iswwp,2,ii) .le. 0.d0 ) then
                        rn = 0.0
                     else
                        rn = wtxwwp(mnwwp(kk,k)+iswwp,1,ii)
     &                     / wtxwwp(mnwwp(kk,k)+iswwp,2,ii)
                     endif

                     write(dum2(lngb+lng1+4:lngb+lng1+65),
     &                         '(1x,1p5e12.4)') rimp,
     &                         wtxwwp(mnwwp(kk,k)+iswwp,1,ii),
     &                         wtxwwp(mnwwp(kk,k)+iswwp,2,ii),
     &                         rn,
     &                         wtxwwp(mnwwp(kk,k)+iswwp,3,ii)

                     write(io,'(600a1)') (dum2(i:i),i=1,lngb+lng1+65)

                  end if
                  end do

               end do

                     dum2(1:13) = ' total events'

                     lstl = 3+lngmax+49

                  do i = 14, lstl

                     dum2(i:i) = ' '

                  end do

                     write(dum2(lstl+1:lstl+13),
     &                         '(1p1e12.4)') t1

                     write(io,'(600a1)') (dum2(i:i),i=1,lstl+13)

*-----------------------------------------------------------------------

            else if( iwmsh .eq. 3 ) then

                     inx = iwxnm(kk)
                     iny = iwynm(kk)
                     inz = iwznm(kk)
                     ixyz = inx * iny * inz

               do jx = 1, inx
               do jy = 1, iny
               do jz = 1, inz

                     dum2(1:2) = ' ('
                     lni = 2 + 1

                  if( laix .eq. 1 ) then
                     write(dum2(lni:lni),'(i1)') jx
                  else if( laix .eq. 2 ) then
                     write(dum2(lni:lni+1),'(i2)') jx
                  else if( laix .eq. 3 ) then
                     write(dum2(lni:lni+2),'(i3)') jx
                  else if( laix .eq. 4 ) then
                     write(dum2(lni:lni+3),'(i4)') jx
                  else if( laix .eq. 5 ) then
                     write(dum2(lni:lni+4),'(i5)') jx
                  end if

                     lni = lni + laix
                     dum2(lni:lni) = ' '
                     lni = lni + 1

                  if( laiy .eq. 1 ) then
                     write(dum2(lni:lni),'(i1)') jy
                  else if( laiy .eq. 2 ) then
                     write(dum2(lni:lni+1),'(i2)') jy
                  else if( laiy .eq. 3 ) then
                     write(dum2(lni:lni+2),'(i3)') jy
                  else if( laiy .eq. 4 ) then
                     write(dum2(lni:lni+3),'(i4)') jy
                  else if( laiy .eq. 5 ) then
                     write(dum2(lni:lni+4),'(i5)') jy
                  end if

                     lni = lni + laiy
                     dum2(lni:lni) = ' '
                     lni = lni + 1

                  if( laiz .eq. 1 ) then
                     write(dum2(lni:lni),'(i1)') jz
                  else if( laiz .eq. 2 ) then
                     write(dum2(lni:lni+1),'(i2)') jz
                  else if( laiz .eq. 3 ) then
                     write(dum2(lni:lni+2),'(i3)') jz
                  else if( laiz .eq. 4 ) then
                     write(dum2(lni:lni+3),'(i4)') jz
                  else if( laiz .eq. 5 ) then
                     write(dum2(lni:lni+4),'(i5)') jz
                  end if

                     lni = lni + laiz
                     dum2(lni:lni+1) = ') '
                     lni = lni + 1

                  ii = jz + inz * ( jy - 1 )+ iny * inz * ( jx - 1 )

                  if( wtxwwp(mnwwp(kk,k)+iswwp,3,ii) .gt. 0.0d0 ) then

                     if( wtxwwp(mnwwp(kk,k)+iswwp,2,ii) .le. 0.d0 ) then
                        rn = 0.0
                     else
                        rn = wtxwwp(mnwwp(kk,k)+iswwp,1,ii)
     &                     / wtxwwp(mnwwp(kk,k)+iswwp,2,ii)
                     endif

                     if( allocated(das_kgwwp) ) then
                        rimp = das_kgwwp(kgwwp(kk)+ii-1)
                     end if

                     write(dum2(lni:lni+60),
     &                         '(1p5e12.4)') rimp,
     &                         wtxwwp(mnwwp(kk,k)+iswwp,1,ii),
     &                         wtxwwp(mnwwp(kk,k)+iswwp,2,ii),
     &                         rn,
     &                         wtxwwp(mnwwp(kk,k)+iswwp,3,ii)

                     write(io,'(600a1)') (dum2(i:i),i=1,lni+60)

                  end if

               end do
               end do
               end do

                     dum2(1:13) = ' total events'

                  do i = 14, lni+47

                     dum2(i:i) = ' '

                  end do

                     write(dum2(lni+48:lni+60),
     &                         '(1p1e12.4)') t1

                     write(io,'(600a1)') (dum2(i:i),i=1,lni+60)

*-----------------------------------------------------------------------

            else if( iwmsh .eq. 4 ) then

                     idsm = inwwt(kk)
                     jdsm = inwwc(kk)
                     kdsm = kfwwp(kk)

                   do j = 1, mnwwp(kk,0)

                     rimp = das_kfwwp(kdsm+j-1)

                     lng1=8
                     write(dum1(1:lng1),'(i8)')
     &                    idas_inwwc(jdsm+1+idas_inwwt(idsm+1+j))

                     dum2(1:3) = '   '

                     do i = 1, lng1
                      dum2(i+3:i+3) = dum1(i:i)
                     end do

                     lngb = lngmax - lng1

                     do i = 1, lngb + 2
                      dum2(i+lng1+3:i+lng1+3) = ' '
                     end do

                     lni = lngb + lng1 + 3 + 2 + 1

                   if( wtxwwp(mnwwp(kk,k)+iswwp,3,j) .gt. 0.0d0 ) then

                     if( wtxwwp(mnwwp(kk,k)+iswwp,2,j) .le. 0.d0 ) then
                        rn = 0.0
                     else
                        rn = wtxwwp(mnwwp(kk,k)+iswwp,1,j)
     &                     / wtxwwp(mnwwp(kk,k)+iswwp,2,j)
                     endif

                     write(dum2(lni:lni+60),
     &                         '(1x,1p5e12.4)') rimp,
     &                         wtxwwp(mnwwp(kk,k)+iswwp,1,j),
     &                         wtxwwp(mnwwp(kk,k)+iswwp,2,j),
     &                         rn,
     &                         wtxwwp(mnwwp(kk,k)+iswwp,3,j)

                     write(io,'(600a1)') (dum2(i:i),i=1,lni+60)

                  end if

                  end do

                     dum2(1:13) = ' total events'

                  do i = 14, lni+47

                     dum2(i:i) = ' '

                  end do

                     write(dum2(lni+48:lni+60),
     &                         '(1p1e12.4)') t1

                     write(io,'(600a1)') (dum2(i:i),i=1,lni+60)

            end if

*-----------------------------------------------------------------------

            end if

         end do
         end do

*-----------------------------------------------------------------------
*        Splitting at Collision
*-----------------------------------------------------------------------

                     igm = ( mmmax - 1 ) * 2 + 1

         do kk = 1, iwwdp

                     kimps = 0

         do k = 1, mnwwp(kk,20)

                     lngmax = 0
                     jsek   = 0
                     t1     = 0.0d0

*-----------------------------------------------------------------------

            if( iwmsh .eq. 1 ) then

                     ii = 0

                     idsm = inwwt(kk)
                     jdsm = 0

               do m = 1, mnwwp(kk,0)

                     jdsm = jdsm + 1
                     ntrn = idas_inwwt(idsm+jdsm)
                     jdsm = jdsm + 1
                     mtrn = idas_inwwt(idsm+jdsm)
                     kssm = idsm + jdsm + 1
                     jdsm = jdsm + mtrn

                     jj = 0

                  do j = 1, ntrn

                     ii = ii + 1

                     call echrg3(jj,mtrn,idas_inwwt(kssm)
     &                           ,dum1,lng1,icmb,igm)

                    if( wtxwwp(mnwwp(kk,k)+iswwp,14,ii) .gt. 0.d0 ) then

                         jsek = jsek + 1

                         if( lng1 .gt. lngmax ) lngmax = lng1

                         t1 = t1 + wtxwwp(mnwwp(kk,k)+iswwp,14,ii)

                    end if

                  end do

               end do

*-----------------------------------------------------------------------

            else if( iwmsh .eq. 3 ) then

                     inx = iwxnm(kk)
                     iny = iwynm(kk)
                     inz = iwznm(kk)
                     ixyz = inx * iny * inz

                  do jx = 1, inx
                  do jy = 1, iny
                  do jz = 1, inz

                     ii = jz + inz * ( jy - 1 )+ iny * inz * ( jx - 1 )

                    if( wtxwwp(mnwwp(kk,k)+iswwp,14,ii) .gt. 0.d0 ) then

                         jsek = jsek + 1

                         t1 = t1 + wtxwwp(mnwwp(kk,k)+iswwp,14,ii)

                     end if

                  end do
                  end do
                  end do

*-----------------------------------------------------------------------

            else if( iwmsh .eq. 4 ) then

                  do j = 1, mnwwp(kk,0)

                     if( wtxwwp(mnwwp(kk,k)+iswwp,3,j) .gt. 0.d0 ) then

                         jsek = jsek + 1

                         t1 = t1 + wtxwwp(mnwwp(kk,k)+iswwp,14,j)

                     end if

                  end do

            end if

*-----------------------------------------------------------------------

            if( jsek .ne. 0 ) then

               if( kimps .eq. 0 ) then

                  write(io,'(/79(''-''))')
                  write(io,'('' Summary of weight window per source'',
     &            '' (non zero) for '',i2,''-th [weight window]''/
     &            '' === Splitting at Collision ==='')') kk
                     write(io,'(79(''-''))')

                  kimps = kimps + 1

               end if

                     uname(k) = pname(mnwwp(kk,k))(1:8)

*-----------------------------------------------------------------------

            if( iwmsh .eq. 1 ) then

                     lblk = lngmax
                     dum2(1:5) = '  reg'

                   do i = 1, lblk
                     dum2(i+5:i+5) = ' '
                   end do

                     lbin = 5 + lblk

*-----------------------------------------------------------------------

            else if( iwmsh .eq. 3 ) then

                     dum2(1:2) = ' ('
                     lni = 2 + 1

                  if( iwxnm(kk) .lt. 10 ) then
                     laix = 1
                     dum2(lni:lni) = 'x'
                  else if( iwxnm(kk) .lt. 100 ) then
                     laix = 2
                     dum2(lni:lni+laix-1) = ' x'
                  else if( iwxnm(kk) .lt. 1000 ) then
                     laix = 3
                     dum2(lni:lni+laix-1) = '  x'
                  else if( iwxnm(kk) .lt. 10000 ) then
                     laix = 4
                     dum2(lni:lni+laix-1) = '   x'
                  else if( iwxnm(kk) .lt. 100000 ) then
                     laix = 5
                     dum2(lni:lni+laix-1) = '    x'
                  end if

                     lni = lni + laix
                     dum2(lni:lni) = ' '
                     lni = lni + 1

                  if( iwynm(kk) .lt. 10 ) then
                     laiy = 1
                     dum2(lni:lni) = 'y'
                  else if( iwynm(kk) .lt. 100 ) then
                     laiy = 2
                     dum2(lni:lni+laiy-1) = ' y'
                  else if( iwynm(kk) .lt. 1000 ) then
                     laiy = 3
                     dum2(lni:lni+laiy-1) = '  y'
                  else if( iwynm(kk) .lt. 10000 ) then
                     laiy = 4
                     dum2(lni:lni+laiy-1) = '   y'
                  else if( iwynm(kk) .lt. 100000 ) then
                     laiy = 5
                     dum2(lni:lni+laix-1) = '    y'
                  end if

                     lni = lni + laiy
                     dum2(lni:lni) = ' '
                     lni = lni + 1

                  if( iwznm(kk) .lt. 10 ) then
                     laiz = 1
                     dum2(lni:lni) = 'z'
                  else if( iwznm(kk) .lt. 100 ) then
                     laiz = 2
                     dum2(lni:lni+laiz-1) = ' z'
                  else if( iwznm(kk) .lt. 1000 ) then
                     laiz = 3
                     dum2(lni:lni+laiz-1) = '  z'
                  else if( iwznm(kk) .lt. 10000 ) then
                     laiz = 4
                     dum2(lni:lni+laiz-1) = '   z'
                  else if( iwznm(kk) .lt. 100000 ) then
                     laiz = 5
                     dum2(lni:lni+laiz-1) = '    z'
                  end if

                     lni = lni + laiz
                     dum2(lni:lni+1) = ') '

                     lbin = lni+1

*-----------------------------------------------------------------------

            else if( iwmsh .eq. 4 ) then

                     lngmax = 8
                     lblk = lngmax
                     dum2(1:6) = '   tet'

                   do i = 1, lblk
                     dum2(i+6:i+6) = ' '
                   end do

                     lbin = lngmax + 6

            end if

*-----------------------------------------------------------------------

                     dum2(lbin+1:lbin+59) =
     &              '[ '// uname(k) //
     &              ']    incoming    splitted    splitted'//
     &              '   splitting'

                     write(io,'(/600a1)') (dum2(i:i),i=1,lbin+59)

                   do i = 1, lbin
                     dum2(i:i) = ' '
                   end do

                     dum2(lbin+1:lbin+59) =
     &              ' ww1 bound       weight      events      weight'//
     &              '      events'

                     write(io,'(600a1)') (dum2(i:i),i=1,lbin+59)

*-----------------------------------------------------------------------

            if( iwmsh .eq. 1 ) then

                     ii = 0

                     idsm = inwwt(kk)
                     jdsm = 0

                     kdsm = kfwwp(kk)

               do m = 1, mnwwp(kk,0)

                     jdsm = jdsm + 1
                     ntrn = idas_inwwt(idsm+jdsm)
                     jdsm = jdsm + 1
                     mtrn = idas_inwwt(idsm+jdsm)
                     kssm = idsm + jdsm + 1
                     jdsm = jdsm + mtrn

                     rimp = das_kfwwp(kdsm+m-1)

                     jj = 0

                  do j = 1, ntrn

                     ii = ii + 1

                     call echrg3(jj,mtrn,idas_inwwt(kssm)
     &                           ,dum1,lng1,icmb,igm)

                  if( wtxwwp(mnwwp(kk,k)+iswwp,14,ii) .gt. 0.0d0 ) then

                     dum2(1:3) = '   '

                     do i = 1, lng1
                        dum2(i+3:i+3) = dum1(i:i)
                     end do

                     lngb = lngmax - lng1

                     do i = 1, lngb
                        dum2(i+lng1+3:i+lng1+3) = ' '
                     end do

                    if( wtxwwp(mnwwp(kk,k)+iswwp,13,ii) .le. 0.d0 ) then
                        rn = 0.0
                    else
                        rn = wtxwwp(mnwwp(kk,k)+iswwp,12,ii)
     &                     / wtxwwp(mnwwp(kk,k)+iswwp,13,ii)
                    endif

                     if( allocated(das_kgwwp) ) then
                        rimp = das_kgwwp(kgwwp(kk)+ii-1)
                     end if

                     write(dum2(lngb+lng1+4:lngb+lng1+65),
     &                         '(1x,1p5e12.4)') rimp,
     &                         wtxwwp(mnwwp(kk,k)+iswwp,12,ii),
     &                         wtxwwp(mnwwp(kk,k)+iswwp,13,ii),
     &                         rn,
     &                         wtxwwp(mnwwp(kk,k)+iswwp,14,ii)

                     write(io,'(600a1)') (dum2(i:i),i=1,lngb+lng1+65)

                  end if
                  end do

               end do

                     dum2(1:13) = ' total events'

                     lstl = 3+lngmax+49

                  do i = 14, lstl

                     dum2(i:i) = ' '

                  end do

                     write(dum2(lstl+1:lstl+13),
     &                         '(1p1e12.4)') t1

                     write(io,'(600a1)') (dum2(i:i),i=1,lstl+13)

*-----------------------------------------------------------------------

            else if( iwmsh .eq. 3 ) then

                     inx = iwxnm(kk)
                     iny = iwynm(kk)
                     inz = iwznm(kk)
                     ixyz = inx * iny * inz

               do jx = 1, inx
               do jy = 1, iny
               do jz = 1, inz

                     dum2(1:2) = ' ('
                     lni = 2 + 1

                  if( laix .eq. 1 ) then
                     write(dum2(lni:lni),'(i1)') jx
                  else if( laix .eq. 2 ) then
                     write(dum2(lni:lni+1),'(i2)') jx
                  else if( laix .eq. 3 ) then
                     write(dum2(lni:lni+2),'(i3)') jx
                  else if( laix .eq. 4 ) then
                     write(dum2(lni:lni+3),'(i4)') jx
                  else if( laix .eq. 5 ) then
                     write(dum2(lni:lni+4),'(i5)') jx
                  end if

                     lni = lni + laix
                     dum2(lni:lni) = ' '
                     lni = lni + 1

                  if( laiy .eq. 1 ) then
                     write(dum2(lni:lni),'(i1)') jy
                  else if( laiy .eq. 2 ) then
                     write(dum2(lni:lni+1),'(i2)') jy
                  else if( laiy .eq. 3 ) then
                     write(dum2(lni:lni+2),'(i3)') jy
                  else if( laiy .eq. 4 ) then
                     write(dum2(lni:lni+3),'(i4)') jy
                  else if( laiy .eq. 5 ) then
                     write(dum2(lni:lni+4),'(i5)') jy
                  end if

                     lni = lni + laiy
                     dum2(lni:lni) = ' '
                     lni = lni + 1

                  if( laiz .eq. 1 ) then
                     write(dum2(lni:lni),'(i1)') jz
                  else if( laiz .eq. 2 ) then
                     write(dum2(lni:lni+1),'(i2)') jz
                  else if( laiz .eq. 3 ) then
                     write(dum2(lni:lni+2),'(i3)') jz
                  else if( laiz .eq. 4 ) then
                     write(dum2(lni:lni+3),'(i4)') jz
                  else if( laiz .eq. 5 ) then
                     write(dum2(lni:lni+4),'(i5)') jz
                  end if

                     lni = lni + laiz
                     dum2(lni:lni+1) = ') '
                     lni = lni + 1

                  ii = jz + inz * ( jy - 1 )+ iny * inz * ( jx - 1 )

                  if( wtxwwp(mnwwp(kk,k)+iswwp,14,ii) .gt. 0.0d0 ) then

                    if( wtxwwp(mnwwp(kk,k)+iswwp,13,ii) .le. 0.d0 ) then
                        rn = 0.0
                    else
                        rn = wtxwwp(mnwwp(kk,k)+iswwp,12,ii)
     &                     / wtxwwp(mnwwp(kk,k)+iswwp,13,ii)
                    endif

                     if( allocated(das_kgwwp) ) then
                        rimp = das_kgwwp(kgwwp(kk)+ii-1)
                     end if

                     write(dum2(lni:lni+60),
     &                         '(1p5e12.4)') rimp,
     &                         wtxwwp(mnwwp(kk,k)+iswwp,12,ii),
     &                         wtxwwp(mnwwp(kk,k)+iswwp,13,ii),
     &                         rn,
     &                         wtxwwp(mnwwp(kk,k)+iswwp,14,ii)

                     write(io,'(600a1)') (dum2(i:i),i=1,lni+60)

                  end if

               end do
               end do
               end do

                     dum2(1:13) = ' total events'

                  do i = 14, lni+47

                     dum2(i:i) = ' '

                  end do

                     write(dum2(lni+48:lni+60),
     &                         '(1p1e12.4)') t1

                     write(io,'(600a1)') (dum2(i:i),i=1,lni+60)

*-----------------------------------------------------------------------

            else if( iwmsh .eq. 4 ) then

                     idsm = inwwt(kk)
                     jdsm = inwwc(kk)
                     kdsm = kfwwp(kk)

                   do j = 1, mnwwp(kk,0)

                     rimp = das_kfwwp(kdsm+j-1)

                     lng1=8
                     write(dum1(1:lng1),'(i8)')
     &                    idas_inwwc(jdsm+1+idas_inwwt(idsm+1+j))

                     dum2(1:3) = '   '

                     do i = 1, lng1
                      dum2(i+3:i+3) = dum1(i:i)
                     end do

                     lngb = lngmax - lng1

                     do i = 1, lngb + 2
                      dum2(i+lng1+3:i+lng1+3) = ' '
                     end do

                     lni = lngb + lng1 + 3 + 2 + 1

                   if( wtxwwp(mnwwp(kk,k)+iswwp,14,j) .gt. 0.0d0 ) then

                     if( wtxwwp(mnwwp(kk,k)+iswwp,13,j) .le. 0.d0 ) then
                        rn = 0.0
                     else
                        rn = wtxwwp(mnwwp(kk,k)+iswwp,12,j)
     &                     / wtxwwp(mnwwp(kk,k)+iswwp,13,j)
                     endif

                     write(dum2(lni:lni+60),
     &                         '(1x,1p5e12.4)') rimp,
     &                         wtxwwp(mnwwp(kk,k)+iswwp,12,j),
     &                         wtxwwp(mnwwp(kk,k)+iswwp,13,j),
     &                         rn,
     &                         wtxwwp(mnwwp(kk,k)+iswwp,14,j)

                     write(io,'(600a1)') (dum2(i:i),i=1,lni+60)

                  end if

                  end do

                     dum2(1:13) = ' total events'

                  do i = 14, lni+47

                     dum2(i:i) = ' '

                  end do

                     write(dum2(lni+48:lni+60),
     &                         '(1p1e12.4)') t1

                     write(io,'(600a1)') (dum2(i:i),i=1,lni+60)

            end if

*-----------------------------------------------------------------------

            end if

         end do
         end do

*-----------------------------------------------------------------------
*        Russian Roulette at Boundary
*-----------------------------------------------------------------------

                     igm = ( mmmax - 1 ) * 2 + 1

         do kk = 1, iwwdp

                     kimps = 0

         do k = 1, mnwwp(kk,20)

                     lngmax = 0
                     jsek   = 0
                     t2     = 0.0d0
                     t3     = 0.0d0

*-----------------------------------------------------------------------

            if( iwmsh .eq. 1 ) then

                     ii = 0

                     idsm = inwwt(kk)
                     jdsm = 0

               do m = 1, mnwwp(kk,0)

                     jdsm = jdsm + 1
                     ntrn = idas_inwwt(idsm+jdsm)
                     jdsm = jdsm + 1
                     mtrn = idas_inwwt(idsm+jdsm)
                     kssm = idsm + jdsm + 1
                     jdsm = jdsm + mtrn

                     jj = 0

                  do j = 1, ntrn

                     ii = ii + 1

                     call echrg3(jj,mtrn,idas_inwwt(kssm)
     &                           ,dum1,lng1,icmb,igm)

                     if( wtxwwp(mnwwp(kk,k)+iswwp,5,ii) .gt. 0.d0 .or.
     &                   wtxwwp(mnwwp(kk,k)+iswwp,7,ii) .gt. 0.d0 ) then

                         jsek = jsek + 1

                         if( lng1 .gt. lngmax ) lngmax = lng1

                         t2 = t2 + wtxwwp(mnwwp(kk,k)+iswwp,5,ii)
                         t3 = t3 + wtxwwp(mnwwp(kk,k)+iswwp,7,ii)

                     end if

                  end do

               end do

*-----------------------------------------------------------------------

            else if( iwmsh .eq. 3 ) then

                     inx = iwxnm(kk)
                     iny = iwynm(kk)
                     inz = iwznm(kk)
                     ixyz = inx * iny * inz

                  do jx = 1, inx
                  do jy = 1, iny
                  do jz = 1, inz

                     ii = jz + inz * ( jy - 1 )+ iny * inz * ( jx - 1 )

                     if( wtxwwp(mnwwp(kk,k)+iswwp,5,ii) .gt. 0.d0 .or.
     &                   wtxwwp(mnwwp(kk,k)+iswwp,7,ii) .gt. 0.d0 ) then

                         jsek = jsek + 1

                         t2 = t2 + wtxwwp(mnwwp(kk,k)+iswwp,5,ii)
                         t3 = t3 + wtxwwp(mnwwp(kk,k)+iswwp,7,ii)

                     end if

                  end do
                  end do
                  end do

*-----------------------------------------------------------------------

            else if( iwmsh .eq. 4 ) then

                  do j = 1, mnwwp(kk,0)

                     if( wtxwwp(mnwwp(kk,k)+iswwp,5,j) .gt. 0.d0 .or.
     &                   wtxwwp(mnwwp(kk,k)+iswwp,7,j) .gt. 0.d0 ) then

                         jsek = jsek + 1

                         t2 = t2 + wtxwwp(mnwwp(kk,k)+iswwp,5,j)
                         t3 = t3 + wtxwwp(mnwwp(kk,k)+iswwp,7,j)

                     end if

                  end do

            end if

*-----------------------------------------------------------------------

            if( jsek .ne. 0 ) then

               if( kimps .eq. 0 ) then

                  write(io,'(/79(''-''))')
                  write(io,'('' Summary of weight window per source'',
     &            '' (non zero) for '',i2,''-th [weight window]''/
     &            '' === Russian Roulette at Boundary ==='')') kk
                     write(io,'(79(''-''))')

                  kimps = kimps + 1

               end if

                     uname(k) = pname(mnwwp(kk,k))(1:8)

*-----------------------------------------------------------------------

            if( iwmsh .eq. 1 ) then

                     lblk = lngmax
                     dum2(1:5) = '  reg'

                   do i = 1, lblk
                     dum2(i+5:i+5) = ' '
                   end do

                     lbin = 5 + lblk

*-----------------------------------------------------------------------

            else if( iwmsh .eq. 3 ) then

                     dum2(1:2) = ' ('
                     lni = 2 + 1

                  if( iwxnm(kk) .lt. 10 ) then
                     laix = 1
                     dum2(lni:lni) = 'x'
                  else if( iwxnm(kk) .lt. 100 ) then
                     laix = 2
                     dum2(lni:lni+laix-1) = ' x'
                  else if( iwxnm(kk) .lt. 1000 ) then
                     laix = 3
                     dum2(lni:lni+laix-1) = '  x'
                  else if( iwxnm(kk) .lt. 10000 ) then
                     laix = 4
                     dum2(lni:lni+laix-1) = '   x'
                  else if( iwxnm(kk) .lt. 100000 ) then
                     laix = 5
                     dum2(lni:lni+laix-1) = '    x'
                  end if

                     lni = lni + laix
                     dum2(lni:lni) = ' '
                     lni = lni + 1

                  if( iwynm(kk) .lt. 10 ) then
                     laiy = 1
                     dum2(lni:lni) = 'y'
                  else if( iwynm(kk) .lt. 100 ) then
                     laiy = 2
                     dum2(lni:lni+laiy-1) = ' y'
                  else if( iwynm(kk) .lt. 1000 ) then
                     laiy = 3
                     dum2(lni:lni+laiy-1) = '  y'
                  else if( iwynm(kk) .lt. 10000 ) then
                     laiy = 4
                     dum2(lni:lni+laiy-1) = '   y'
                  else if( iwynm(kk) .lt. 100000 ) then
                     laiy = 5
                     dum2(lni:lni+laix-1) = '    y'
                  end if

                     lni = lni + laiy
                     dum2(lni:lni) = ' '
                     lni = lni + 1

                  if( iwznm(kk) .lt. 10 ) then
                     laiz = 1
                     dum2(lni:lni) = 'z'
                  else if( iwznm(kk) .lt. 100 ) then
                     laiz = 2
                     dum2(lni:lni+laiz-1) = ' z'
                  else if( iwznm(kk) .lt. 1000 ) then
                     laiz = 3
                     dum2(lni:lni+laiz-1) = '  z'
                  else if( iwznm(kk) .lt. 10000 ) then
                     laiz = 4
                     dum2(lni:lni+laiz-1) = '   z'
                  else if( iwznm(kk) .lt. 100000 ) then
                     laiz = 5
                     dum2(lni:lni+laiz-1) = '    z'
                  end if

                     lni = lni + laiz
                     dum2(lni:lni+1) = ') '

                     lbin = lni+1

*-----------------------------------------------------------------------

            else if( iwmsh .eq. 4 ) then

                     lngmax = 8
                     lblk = lngmax
                     dum2(1:6) = '   tet'

                   do i = 1, lblk
                     dum2(i+6:i+6) = ' '
                   end do

                     lbin = lngmax + 6

            end if

*-----------------------------------------------------------------------

                     dum2(lbin+1:lbin+59) =
     &              '[ '// uname(k) //
     &              ']    survived    survived      killed'//
     &              '      killed'

                     write(io,'(/600a1)') (dum2(i:i),i=1,lbin+59)

                   do i = 1, lbin
                     dum2(i:i) = ' '
                   end do

                     dum2(lbin+1:lbin+59) =
     &              ' ww1 bound       weight      events      weight'//
     &              '      events'

                     write(io,'(600a1)') (dum2(i:i),i=1,lbin+59)

*-----------------------------------------------------------------------

            if( iwmsh .eq. 1 ) then

                     ii = 0

                     idsm = inwwt(kk)
                     jdsm = 0

                     kdsm = kfwwp(kk)

               do m = 1, mnwwp(kk,0)

                     jdsm = jdsm + 1
                     ntrn = idas_inwwt(idsm+jdsm)
                     jdsm = jdsm + 1
                     mtrn = idas_inwwt(idsm+jdsm)
                     kssm = idsm + jdsm + 1
                     jdsm = jdsm + mtrn

                     rimp = das_kfwwp(kdsm+m-1)

                     jj = 0

                  do j = 1, ntrn

                     ii = ii + 1

                     call echrg3(jj,mtrn,idas_inwwt(kssm)
     &                           ,dum1,lng1,icmb,igm)

                  if( wtxwwp(mnwwp(kk,k)+iswwp,5,ii) .gt. 0.0 .or.
     &                wtxwwp(mnwwp(kk,k)+iswwp,7,ii) .gt. 0.0 ) then

                     dum2(1:3) = '   '

                     do i = 1, lng1
                        dum2(i+3:i+3) = dum1(i:i)
                     end do

                     lngb = lngmax - lng1

                     do i = 1, lngb
                        dum2(i+lng1+3:i+lng1+3) = ' '
                     end do

                     write(dum2(lngb+lng1+4:lngb+lng1+65),
     &                         '(1x,1p5e12.4)') rimp,
     &                         wtxwwp(mnwwp(kk,k)+iswwp,4,ii),
     &                         wtxwwp(mnwwp(kk,k)+iswwp,5,ii),
     &                         wtxwwp(mnwwp(kk,k)+iswwp,6,ii),
     &                         wtxwwp(mnwwp(kk,k)+iswwp,7,ii)

                     write(io,'(600a1)') (dum2(i:i),i=1,lngb+lng1+65)

                  end if
                  end do

               end do

                     dum2(1:13) = ' total events'

                     lstl = 3+lngmax+25

                  do i = 14, lstl

                     dum2(i:i) = ' '

                  end do

                     write(dum2(lstl+1:lstl+13),
     &                         '(1p1e12.4)') t2

                     lstl = 3+lngmax+49

                  do i = 3+lngmax+38, lstl

                     dum2(i:i) = ' '

                  end do

                     write(dum2(lstl+1:lstl+13),
     &                         '(1p1e12.4)') t3

                     write(io,'(600a1)') (dum2(i:i),i=1,lstl+13)

*-----------------------------------------------------------------------

            else if( iwmsh .eq. 3 ) then

                     inx = iwxnm(kk)
                     iny = iwynm(kk)
                     inz = iwznm(kk)
                     ixyz = inx * iny * inz

               do jx = 1, inx
               do jy = 1, iny
               do jz = 1, inz

                     dum2(1:2) = ' ('
                     lni = 2 + 1

                  if( laix .eq. 1 ) then
                     write(dum2(lni:lni),'(i1)') jx
                  else if( laix .eq. 2 ) then
                     write(dum2(lni:lni+1),'(i2)') jx
                  else if( laix .eq. 3 ) then
                     write(dum2(lni:lni+2),'(i3)') jx
                  else if( laix .eq. 4 ) then
                     write(dum2(lni:lni+3),'(i4)') jx
                  else if( laix .eq. 5 ) then
                     write(dum2(lni:lni+4),'(i5)') jx
                  end if

                     lni = lni + laix
                     dum2(lni:lni) = ' '
                     lni = lni + 1

                  if( laiy .eq. 1 ) then
                     write(dum2(lni:lni),'(i1)') jy
                  else if( laiy .eq. 2 ) then
                     write(dum2(lni:lni+1),'(i2)') jy
                  else if( laiy .eq. 3 ) then
                     write(dum2(lni:lni+2),'(i3)') jy
                  else if( laiy .eq. 4 ) then
                     write(dum2(lni:lni+3),'(i4)') jy
                  else if( laiy .eq. 5 ) then
                     write(dum2(lni:lni+4),'(i5)') jy
                  end if

                     lni = lni + laiy
                     dum2(lni:lni) = ' '
                     lni = lni + 1

                  if( laiz .eq. 1 ) then
                     write(dum2(lni:lni),'(i1)') jz
                  else if( laiz .eq. 2 ) then
                     write(dum2(lni:lni+1),'(i2)') jz
                  else if( laiz .eq. 3 ) then
                     write(dum2(lni:lni+2),'(i3)') jz
                  else if( laiz .eq. 4 ) then
                     write(dum2(lni:lni+3),'(i4)') jz
                  else if( laiz .eq. 5 ) then
                     write(dum2(lni:lni+4),'(i5)') jz
                  end if

                     lni = lni + laiz
                     dum2(lni:lni+1) = ') '
                     lni = lni + 1

                  ii = jz + inz * ( jy - 1 )+ iny * inz * ( jx - 1 )

                  if( wtxwwp(mnwwp(kk,k)+iswwp,5,ii) .gt. 0.0 .or.
     &                wtxwwp(mnwwp(kk,k)+iswwp,7,ii) .gt. 0.0 ) then

                     if( allocated(das_kgwwp) ) then
                        rimp = das_kgwwp(kgwwp(kk)+ii-1)
                     end if

                     write(dum2(lni:lni+60),
     &                         '(1p5e12.4)') rimp,
     &                         wtxwwp(mnwwp(kk,k)+iswwp,4,ii),
     &                         wtxwwp(mnwwp(kk,k)+iswwp,5,ii),
     &                         wtxwwp(mnwwp(kk,k)+iswwp,6,ii),
     &                         wtxwwp(mnwwp(kk,k)+iswwp,7,ii)

                     write(io,'(600a1)') (dum2(i:i),i=1,lni+60)

                  end if

               end do
               end do
               end do

                     dum2(1:13) = ' total events'

                  do i = 14, lni+23

                     dum2(i:i) = ' '

                  end do

                     write(dum2(lni+24:lni+36),
     &                         '(1p1e12.4)') t2

                  do i = lni+37, lni+48

                     dum2(i:i) = ' '

                  end do

                     write(dum2(lni+48:lni+60),
     &                         '(1p1e12.4)') t3

                     write(io,'(600a1)') (dum2(i:i),i=1,lni+60)

*-----------------------------------------------------------------------

            else if( iwmsh .eq. 4 ) then

                     idsm = inwwt(kk)
                     jdsm = inwwc(kk)
                     kdsm = kfwwp(kk)

                   do j = 1, mnwwp(kk,0)

                     rimp = das_kfwwp(kdsm+j-1)

                     lng1=8
                     write(dum1(1:lng1),'(i8)')
     &                    idas_inwwc(jdsm+1+idas_inwwt(idsm+1+j))

                     dum2(1:3) = '   '

                     do i = 1, lng1
                      dum2(i+3:i+3) = dum1(i:i)
                     end do

                     lngb = lngmax - lng1

                     do i = 1, lngb + 2
                      dum2(i+lng1+3:i+lng1+3) = ' '
                     end do

                     lni = lngb + lng1 + 3 + 2 + 1

                  if( wtxwwp(mnwwp(kk,k)+iswwp,5,j) .gt. 0.0 .or.
     &                wtxwwp(mnwwp(kk,k)+iswwp,7,j) .gt. 0.0 ) then

                     if( allocated(das_kgwwp) ) then
                        rimp = das_kgwwp(kgwwp(kk)+j-1)
                     end if

                     write(dum2(lni:lni+60),
     &                         '(1x,1p5e12.4)') rimp,
     &                         wtxwwp(mnwwp(kk,k)+iswwp,4,j),
     &                         wtxwwp(mnwwp(kk,k)+iswwp,5,j),
     &                         wtxwwp(mnwwp(kk,k)+iswwp,6,j),
     &                         wtxwwp(mnwwp(kk,k)+iswwp,7,j)

                     write(io,'(600a1)') (dum2(i:i),i=1,lni+60)

                  end if

                  end do

                     dum2(1:13) = ' total events'

                  do i = 14, lni+23

                     dum2(i:i) = ' '

                  end do

                     write(dum2(lni+24:lni+36),
     &                         '(1p1e12.4)') t2

                  do i = lni+37, lni+48

                     dum2(i:i) = ' '

                  end do

                     write(dum2(lni+48:lni+60),
     &                         '(1p1e12.4)') t3

                     write(io,'(600a1)') (dum2(i:i),i=1,lni+60)

            end if

*-----------------------------------------------------------------------

            end if

         end do
         end do

*-----------------------------------------------------------------------
*        Russian Roulette at Collision
*-----------------------------------------------------------------------

                     igm = ( mmmax - 1 ) * 2 + 1

         do kk = 1, iwwdp

                     kimps = 0

         do k = 1, mnwwp(kk,20)

                     lngmax = 0
                     jsek   = 0
                     t2     = 0.0d0
                     t3     = 0.0d0

*-----------------------------------------------------------------------

            if( iwmsh .eq. 1 ) then

                     ii = 0

                     idsm = inwwt(kk)
                     jdsm = 0

               do m = 1, mnwwp(kk,0)

                     jdsm = jdsm + 1
                     ntrn = idas_inwwt(idsm+jdsm)
                     jdsm = jdsm + 1
                     mtrn = idas_inwwt(idsm+jdsm)
                     kssm = idsm + jdsm + 1
                     jdsm = jdsm + mtrn

                     jj = 0

                  do j = 1, ntrn

                     ii = ii + 1

                     call echrg3(jj,mtrn,idas_inwwt(kssm)
     &                           ,dum1,lng1,icmb,igm)

                    if( wtxwwp(mnwwp(kk,k)+iswwp, 9,ii) .gt. 0.d0 .or.
     &                  wtxwwp(mnwwp(kk,k)+iswwp,11,ii) .gt. 0.d0 ) then

                         jsek = jsek + 1

                         if( lng1 .gt. lngmax ) lngmax = lng1

                         t2 = t2 + wtxwwp(mnwwp(kk,k)+iswwp, 9,ii)
                         t3 = t3 + wtxwwp(mnwwp(kk,k)+iswwp,11,ii)

                     end if

                  end do

               end do

*-----------------------------------------------------------------------

            else if( iwmsh .eq. 3 ) then

                     inx = iwxnm(kk)
                     iny = iwynm(kk)
                     inz = iwznm(kk)
                     ixyz = inx * iny * inz

                  do jx = 1, inx
                  do jy = 1, iny
                  do jz = 1, inz

                     ii = jz + inz * ( jy - 1 )+ iny * inz * ( jx - 1 )

                    if( wtxwwp(mnwwp(kk,k)+iswwp, 9,ii) .gt. 0.d0 .or.
     &                  wtxwwp(mnwwp(kk,k)+iswwp,11,ii) .gt. 0.d0 ) then

                         jsek = jsek + 1

                         t2 = t2 + wtxwwp(mnwwp(kk,k)+iswwp, 9,ii)
                         t3 = t3 + wtxwwp(mnwwp(kk,k)+iswwp,11,ii)

                     end if

                  end do
                  end do
                  end do

*-----------------------------------------------------------------------

            else if( iwmsh .eq. 4 ) then

                  do j = 1, mnwwp(kk,0)

                     if( wtxwwp(mnwwp(kk,k)+iswwp,9,j) .gt. 0.d0 .or.
     &                   wtxwwp(mnwwp(kk,k)+iswwp,11,j) .gt. 0.d0 ) then

                         jsek = jsek + 1

                         t2 = t2 + wtxwwp(mnwwp(kk,k)+iswwp,9,j)
                         t3 = t3 + wtxwwp(mnwwp(kk,k)+iswwp,11,j)

                     end if

                  end do

            end if

*-----------------------------------------------------------------------

            if( jsek .ne. 0 ) then

               if( kimps .eq. 0 ) then

                  write(io,'(/79(''-''))')
                  write(io,'('' Summary of weight window per source'',
     &            '' (non zero) for '',i2,''-th [weight window]''/
     &            '' === Russian Roulette at Collision ==='')') kk
                     write(io,'(79(''-''))')

                  kimps = kimps + 1

               end if

                     uname(k) = pname(mnwwp(kk,k))(1:8)

*-----------------------------------------------------------------------

            if( iwmsh .eq. 1 ) then

                     lblk = lngmax
                     dum2(1:5) = '  reg'

                   do i = 1, lblk
                     dum2(i+5:i+5) = ' '
                   end do

                     lbin = 5 + lblk

*-----------------------------------------------------------------------

            else if( iwmsh .eq. 3 ) then

                     dum2(1:2) = ' ('
                     lni = 2 + 1

                  if( iwxnm(kk) .lt. 10 ) then
                     laix = 1
                     dum2(lni:lni) = 'x'
                  else if( iwxnm(kk) .lt. 100 ) then
                     laix = 2
                     dum2(lni:lni+laix-1) = ' x'
                  else if( iwxnm(kk) .lt. 1000 ) then
                     laix = 3
                     dum2(lni:lni+laix-1) = '  x'
                  else if( iwxnm(kk) .lt. 10000 ) then
                     laix = 4
                     dum2(lni:lni+laix-1) = '   x'
                  else if( iwxnm(kk) .lt. 100000 ) then
                     laix = 5
                     dum2(lni:lni+laix-1) = '    x'
                  end if

                     lni = lni + laix
                     dum2(lni:lni) = ' '
                     lni = lni + 1

                  if( iwynm(kk) .lt. 10 ) then
                     laiy = 1
                     dum2(lni:lni) = 'y'
                  else if( iwynm(kk) .lt. 100 ) then
                     laiy = 2
                     dum2(lni:lni+laiy-1) = ' y'
                  else if( iwynm(kk) .lt. 1000 ) then
                     laiy = 3
                     dum2(lni:lni+laiy-1) = '  y'
                  else if( iwynm(kk) .lt. 10000 ) then
                     laiy = 4
                     dum2(lni:lni+laiy-1) = '   y'
                  else if( iwynm(kk) .lt. 100000 ) then
                     laiy = 5
                     dum2(lni:lni+laix-1) = '    y'
                  end if

                     lni = lni + laiy
                     dum2(lni:lni) = ' '
                     lni = lni + 1

                  if( iwznm(kk) .lt. 10 ) then
                     laiz = 1
                     dum2(lni:lni) = 'z'
                  else if( iwznm(kk) .lt. 100 ) then
                     laiz = 2
                     dum2(lni:lni+laiz-1) = ' z'
                  else if( iwznm(kk) .lt. 1000 ) then
                     laiz = 3
                     dum2(lni:lni+laiz-1) = '  z'
                  else if( iwznm(kk) .lt. 10000 ) then
                     laiz = 4
                     dum2(lni:lni+laiz-1) = '   z'
                  else if( iwznm(kk) .lt. 100000 ) then
                     laiz = 5
                     dum2(lni:lni+laiz-1) = '    z'
                  end if

                     lni = lni + laiz
                     dum2(lni:lni+1) = ') '

                     lbin = lni+1

*-----------------------------------------------------------------------

            else if( iwmsh .eq. 4 ) then

                     lngmax = 8
                     lblk = lngmax
                     dum2(1:6) = '   tet'

                   do i = 1, lblk
                     dum2(i+6:i+6) = ' '
                   end do

                     lbin = lngmax + 6

            end if

*-----------------------------------------------------------------------

                     dum2(lbin+1:lbin+59) =
     &              '[ '// uname(k) //
     &              ']    survived    survived      killed'//
     &              '      killed'

                     write(io,'(/600a1)') (dum2(i:i),i=1,lbin+59)

                   do i = 1, lbin
                     dum2(i:i) = ' '
                   end do

                     dum2(lbin+1:lbin+59) =
     &              ' ww1 bound       weight      events      weight'//
     &              '      events'

                     write(io,'(600a1)') (dum2(i:i),i=1,lbin+59)

*-----------------------------------------------------------------------

            if( iwmsh .eq. 1 ) then

                     ii = 0

                     idsm = inwwt(kk)
                     jdsm = 0

                     kdsm = kfwwp(kk)

               do m = 1, mnwwp(kk,0)

                     jdsm = jdsm + 1
                     ntrn = idas_inwwt(idsm+jdsm)
                     jdsm = jdsm + 1
                     mtrn = idas_inwwt(idsm+jdsm)
                     kssm = idsm + jdsm + 1
                     jdsm = jdsm + mtrn

                     rimp = das_kfwwp(kdsm+m-1)

                     jj = 0

                  do j = 1, ntrn

                     ii = ii + 1

                     call echrg3(jj,mtrn,idas_inwwt(kssm)
     &                           ,dum1,lng1,icmb,igm)

                  if( wtxwwp(mnwwp(kk,k)+iswwp, 9,ii) .gt. 0.0 .or.
     &                wtxwwp(mnwwp(kk,k)+iswwp,11,ii) .gt. 0.0 ) then

                     dum2(1:3) = '   '

                     do i = 1, lng1
                        dum2(i+3:i+3) = dum1(i:i)
                     end do

                     lngb = lngmax - lng1

                     do i = 1, lngb
                        dum2(i+lng1+3:i+lng1+3) = ' '
                     end do

                     write(dum2(lngb+lng1+4:lngb+lng1+65),
     &                         '(1x,1p5e12.4)') rimp,
     &                         wtxwwp(mnwwp(kk,k)+iswwp, 8,ii),
     &                         wtxwwp(mnwwp(kk,k)+iswwp, 9,ii),
     &                         wtxwwp(mnwwp(kk,k)+iswwp,10,ii),
     &                         wtxwwp(mnwwp(kk,k)+iswwp,11,ii)

                     write(io,'(600a1)') (dum2(i:i),i=1,lngb+lng1+65)

                  end if
                  end do

               end do

                     dum2(1:13) = ' total events'

                     lstl = 3+lngmax+25

                  do i = 14, lstl

                     dum2(i:i) = ' '

                  end do

                     write(dum2(lstl+1:lstl+13),
     &                         '(1p1e12.4)') t2

                     lstl = 3+lngmax+49

                  do i = 3+lngmax+38, lstl

                     dum2(i:i) = ' '

                  end do

                     write(dum2(lstl+1:lstl+13),
     &                         '(1p1e12.4)') t3

                     write(io,'(600a1)') (dum2(i:i),i=1,lstl+13)

*-----------------------------------------------------------------------

            else if( iwmsh .eq. 3 ) then

                     inx = iwxnm(kk)
                     iny = iwynm(kk)
                     inz = iwznm(kk)
                     ixyz = inx * iny * inz

               do jx = 1, inx
               do jy = 1, iny
               do jz = 1, inz

                     dum2(1:2) = ' ('
                     lni = 2 + 1

                  if( laix .eq. 1 ) then
                     write(dum2(lni:lni),'(i1)') jx
                  else if( laix .eq. 2 ) then
                     write(dum2(lni:lni+1),'(i2)') jx
                  else if( laix .eq. 3 ) then
                     write(dum2(lni:lni+2),'(i3)') jx
                  else if( laix .eq. 4 ) then
                     write(dum2(lni:lni+3),'(i4)') jx
                  else if( laix .eq. 5 ) then
                     write(dum2(lni:lni+4),'(i5)') jx
                  end if

                     lni = lni + laix
                     dum2(lni:lni) = ' '
                     lni = lni + 1

                  if( laiy .eq. 1 ) then
                     write(dum2(lni:lni),'(i1)') jy
                  else if( laiy .eq. 2 ) then
                     write(dum2(lni:lni+1),'(i2)') jy
                  else if( laiy .eq. 3 ) then
                     write(dum2(lni:lni+2),'(i3)') jy
                  else if( laiy .eq. 4 ) then
                     write(dum2(lni:lni+3),'(i4)') jy
                  else if( laiy .eq. 5 ) then
                     write(dum2(lni:lni+4),'(i5)') jy
                  end if

                     lni = lni + laiy
                     dum2(lni:lni) = ' '
                     lni = lni + 1

                  if( laiz .eq. 1 ) then
                     write(dum2(lni:lni),'(i1)') jz
                  else if( laiz .eq. 2 ) then
                     write(dum2(lni:lni+1),'(i2)') jz
                  else if( laiz .eq. 3 ) then
                     write(dum2(lni:lni+2),'(i3)') jz
                  else if( laiz .eq. 4 ) then
                     write(dum2(lni:lni+3),'(i4)') jz
                  else if( laiz .eq. 5 ) then
                     write(dum2(lni:lni+4),'(i5)') jz
                  end if

                     lni = lni + laiz
                     dum2(lni:lni) = ')'
                     lni = lni + 1

                  ii = jz + inz * ( jy - 1 )+ iny * inz * ( jx - 1 )

                  if( wtxwwp(mnwwp(kk,k)+iswwp, 9,ii) .gt. 0.0 .or.
     &                wtxwwp(mnwwp(kk,k)+iswwp,11,ii) .gt. 0.0 ) then

                     if( allocated(das_kgwwp) ) then
                        rimp = das_kgwwp(kgwwp(kk)+ii-1)
                     end if

                     write(dum2(lni:lni+60),
     &                         '(1p5e12.4)') rimp,
     &                         wtxwwp(mnwwp(kk,k)+iswwp, 8,ii),
     &                         wtxwwp(mnwwp(kk,k)+iswwp, 9,ii),
     &                         wtxwwp(mnwwp(kk,k)+iswwp,10,ii),
     &                         wtxwwp(mnwwp(kk,k)+iswwp,11,ii)

                     write(io,'(600a1)') (dum2(i:i),i=1,lni+60)

                  end if

               end do
               end do
               end do

                     dum2(1:13) = ' total events'

                  do i = 14, lni+23

                     dum2(i:i) = ' '

                  end do

                     write(dum2(lni+24:lni+36),
     &                         '(1p1e12.4)') t2

                  do i = lni+37, lni+48

                     dum2(i:i) = ' '

                  end do

                     write(dum2(lni+48:lni+60),
     &                         '(1p1e12.4)') t3

                     write(io,'(600a1)') (dum2(i:i),i=1,lni+60)

*-----------------------------------------------------------------------

            else if( iwmsh .eq. 4 ) then

                     idsm = inwwt(kk)
                     jdsm = inwwc(kk)
                     kdsm = kfwwp(kk)

                   do j = 1, mnwwp(kk,0)

                     rimp = das_kfwwp(kdsm+j-1)

                     lng1=8
                     write(dum1(1:lng1),'(i8)')
     &                    idas_inwwc(jdsm+1+idas_inwwt(idsm+1+j))

                     dum2(1:3) = '   '

                     do i = 1, lng1
                      dum2(i+3:i+3) = dum1(i:i)
                     end do

                     lngb = lngmax - lng1

                     do i = 1, lngb + 2
                      dum2(i+lng1+3:i+lng1+3) = ' '
                     end do

                     lni = lngb + lng1 + 3 + 2 + 1

                  if( wtxwwp(mnwwp(kk,k)+iswwp,9,j) .gt. 0.0 .or.
     &                wtxwwp(mnwwp(kk,k)+iswwp,11,j) .gt. 0.0 ) then

                     if( allocated(das_kgwwp) ) then
                        rimp = das_kgwwp(kgwwp(kk)+j-1)
                     end if

                     write(dum2(lni:lni+60),
     &                         '(1x,1p5e12.4)') rimp,
     &                         wtxwwp(mnwwp(kk,k)+iswwp,8,j),
     &                         wtxwwp(mnwwp(kk,k)+iswwp,9,j),
     &                         wtxwwp(mnwwp(kk,k)+iswwp,10,j),
     &                         wtxwwp(mnwwp(kk,k)+iswwp,11,j)

                     write(io,'(600a1)') (dum2(i:i),i=1,lni+60)

                  end if

                  end do

                     dum2(1:13) = ' total events'

                     lstl = 3+lngmax+25

                  do i = 14, lstl

                     dum2(i:i) = ' '

                  end do

                     write(dum2(lstl+1:lstl+13),
     &                         '(1p1e12.4)') t2

                     lstl = 3+lngmax+49

                  do i = 3+lngmax+38, lstl

                     dum2(i:i) = ' '

                  end do

                     write(dum2(lstl+1:lstl+13),
     &                         '(1p1e12.4)') t3

                     write(io,'(600a1)') (dum2(i:i),i=1,lstl+13)

            end if

*-----------------------------------------------------------------------

            end if

         end do
         end do

*-----------------------------------------------------------------------

      end if
      end if

*-----------------------------------------------------------------------
*     summary of weight cutoff
*-----------------------------------------------------------------------

      if( iwt .ne. 0 ) then

         if( npe .gt. 1 ) then

            do kk = 1, 20

                  iini = ( kk - 1 ) * maxip * 8 !FURUTA20130820

               if( me .gt. 0 ) then

                        call parasr
     &                     (wtxcut_pointer(iswct+iini+1), maxip*8, 0)

               else

                        do j = 1, maxip * 8

                           wtxcut_pointer(iswct+j+iini) = 0.0d0

                        end do

                  do i = 1, npe - 1

                     if( iccp(i) .eq. 0 ) then

                        call pararr(das_isstr(isstr+1), maxip*8, i)

                        do j = 1, maxip * 8

                           wtxcut_pointer(iswct+j+iini)
     &                        = wtxcut_pointer(iswct+j+iini)
     &                          + das_isstr(isstr+j)

                        end do

                     end if

                  end do

               end if

            end do

         end if

*-----------------------------------------------------------------------

      if( me .eq. 0 ) then

            do ip = 1, 19
            do id = 1, 8
            do ir = 1, maxip

               wtxcut(ip+iswct,id,ir) = wtxcut(ip+iswct,id,ir) / totevt

            end do
            end do
            end do

*-----------------------------------------------------------------------

            do ip = 1, 20

                  k1 = 48 + ip
                  k2 = 68 + ip
                  k3 = 88 + ip

                  do l = 1, ipara

                     if( k1 .eq. -ipsq(l) .or.
     &                   k2 .eq. -ipsq(l) .or.
     &                   k3 .eq. -ipsq(l) ) goto 560

                  end do

            end do

            goto 570

  560       continue

            write(io,'(/79(''-''))')
            write(io,'('' Summary of weight cutoff factors'',
     &                 '' other than default values.'')')
            write(io,'( 79(''-''))')

            write(io,'(/'' particle     '',
     &                  ''   wc1         wc2         swtm'')')

            do 540 ip = 1, 20

                  k1 = 48 + ip
                  k2 = 68 + ip
                  k3 = 88 + ip

                  do l = 1, ipara

                     if( k1 .eq. -ipsq(l) .or.
     &                   k2 .eq. -ipsq(l) .or.
     &                   k3 .eq. -ipsq(l) ) goto 550

                  end do

                  goto 540

  550          continue

               write(io,'(2x,a8,3x,1p9e12.4)') pname(ip),
     &                              wc1(ip),wc2(ip),swtm(ip)

  540       continue

  570       continue

*-----------------------------------------------------------------------
*        transport and produced particles
*-----------------------------------------------------------------------

                     igm = ( mmmax - 1 ) * 2 + 1

         do ll = 0, 4, 4

         do kk = 1, iimpn + 1

                     kimps = 0

         do k = 1, mnimp(kk,20)

                     lngmax = 0
                     jsek   = 0
                     t1     = 0.0d0
                     t2     = 0.0d0

                     ii = 0

                     idsm = inimt(kk)
                     jdsm = 0

               do m = 1, mnimp(kk,0) + 1

                     jdsm = jdsm + 1
                     ntrn = idas_inimt(idsm+jdsm)
                     jdsm = jdsm + 1
                     mtrn = idas_inimt(idsm+jdsm)
                     kssm = idsm + jdsm + 1
                     jdsm = jdsm + mtrn

                     jj = 0

                  do j = 1, ntrn

                     ii = ii + 1

                     call echrg3(jj,mtrn,idas_inimt(kssm)
     &                           ,dum1,lng1,icmb,igm)

                     if( wtxcut(mnimp(kk,k)+iswct,2+ll,ii) .gt. 0.d0
     &                   .or.
     &                   wtxcut(mnimp(kk,k)+iswct,4+ll,ii) .gt. 0.d0 )
     &               then

                         jsek = jsek + 1

                         if( lng1 .gt. lngmax ) lngmax = lng1

                         t1 = t1 + wtxcut(mnimp(kk,k)+iswct,2+ll,ii)
                         t2 = t2 + wtxcut(mnimp(kk,k)+iswct,4+ll,ii)

                     end if

                  end do

               end do

*-----------------------------------------------------------------------

            if( jsek .ne. 0 ) then

               if( kimps .eq. 0 ) then

                  write(io,'(/79(''-''))')

                  if( ll .eq. 0 ) then

                  write(io,'('' Summary of weight cutoff per source'',
     &            '' (non zero)''/
     &         '' === Russian Roulette for Transport Particles ==='')')
                     write(io,'(79(''-''))')

                  else

                  write(io,'('' Summary of weight cutoff per source'',
     &            '' (non zero)''/
     &         '' === Russian Roulette for Produced Particles ==='')')
                     write(io,'(79(''-''))')

                  end if

                  kimps = kimps + 1

               end if

                     uname(k) = pname(mnimp(kk,k))(1:8)

*-----------------------------------------------------------------------

                     lblk = lngmax
                     dum2(1:5) = '  reg'

                   do i = 1, lblk
                     dum2(i+5:i+5) = ' '
                   end do

                     dum2(lblk+5+1:lblk+5+59) =
     &              '[ '// uname(k) //
     &              ']    survived    survived      killed'//
     &              '      killed'

                     write(io,'(/600a1)') (dum2(i:i),i=1,lblk+5+59)

                     dum2(1:5) = '     '
                   do i = 1, lblk
                     dum2(i+5:i+5) = ' '
                   end do

                     dum2(lblk+5+1:lblk+5+59) =
     &              ' importance      weight      events      weight'//
     &              '      events'

                     write(io,'(600a1)') (dum2(i:i),i=1,lblk+5+59)

                     ii = 0

                     idsm = inimt(kk)
                     jdsm = 0

                     kdsm = kfimp(kk)

               do m = 1, mnimp(kk,0) + 1

                     jdsm = jdsm + 1
                     ntrn = idas_inimt(idsm+jdsm)
                     jdsm = jdsm + 1
                     mtrn = idas_inimt(idsm+jdsm)
                     kssm = idsm + jdsm + 1
                     jdsm = jdsm + mtrn

                     rimp = das_kfimp(kdsm-1+m)

                     jj = 0

                  do j = 1, ntrn

                     ii = ii + 1

                     call echrg3(jj,mtrn,idas_inimt(kssm)
     &                           ,dum1,lng1,icmb,igm)

                  if( wtxcut(mnimp(kk,k)+iswct,2+ll,ii) .gt. 0.0 .or.
     &                wtxcut(mnimp(kk,k)+iswct,4+ll,ii) .gt. 0.0 ) then

                     dum2(1:3) = '   '

                     do i = 1, lng1
                        dum2(i+3:i+3) = dum1(i:i)
                     end do

                     lngb = lngmax - lng1

                     do i = 1, lngb
                        dum2(i+lng1+3:i+lng1+3) = ' '
                     end do

                     write(dum2(lngb+lng1+4:lngb+lng1+65),
     &                         '(1x,1p5e12.4)') rimp,
     &                         wtxcut(mnimp(kk,k)+iswct,1+ll,ii),
     &                         wtxcut(mnimp(kk,k)+iswct,2+ll,ii),
     &                         wtxcut(mnimp(kk,k)+iswct,3+ll,ii),
     &                         wtxcut(mnimp(kk,k)+iswct,4+ll,ii)

                     write(io,'(600a1)') (dum2(i:i),i=1,lngb+lng1+65)

                  end if
                  end do

               end do

                     dum2(1:13) = ' total events'

                     lstl = 3+lngmax+25

                  do i = 14, lstl

                     dum2(i:i) = ' '

                  end do

                     write(dum2(lstl+1:lstl+13),
     &                         '(1p1e12.4)') t1

                     lstl = 3+lngmax+49

                  do i = 3+lngmax+38, lstl

                     dum2(i:i) = ' '

                  end do

                     write(dum2(lstl+1:lstl+13),
     &                         '(1p1e12.4)') t2

                     write(io,'(600a1)') (dum2(i:i),i=1,lstl+13)

            end if

         end do
         end do
         end do

*-----------------------------------------------------------------------

      end if
      end if

*-----------------------------------------------------------------------
*     summary of forced collisions
*-----------------------------------------------------------------------

      if( ifcln .ne. 0 ) then

         if( npe .gt. 1 ) then

            do kk = 1, 20

                  iini = ( kk - 1 ) * maxrg * 2 !FURUTA20130820

               if( me .gt. 0 ) then

                        call parasr
     &                     (wtxfcl_pointer(isfcl+iini+1), maxrg*2, 0)

               else

                        do j = 1, maxrg * 2

                           wtxfcl_pointer(isfcl+j+iini) = 0.0d0

                        end do

                  do i = 1, npe - 1

                     if( iccp(i) .eq. 0 ) then

                        call pararr(das_isstr(isstr+1), maxrg*2, i)

                        do j = 1, maxrg * 2

                           wtxfcl_pointer(isfcl+j+iini)
     &                        = wtxfcl_pointer(isfcl+j+iini)
     &                          + das_isstr(isstr+j)

                        end do

                     end if

                  end do

               end if

            end do

         end if

*-----------------------------------------------------------------------

      if( me .eq. 0 ) then

            do ip = 1, 19
            do id = 1, 2
            do ir = 1, maxrg

               wtxfcl(ip+isfcl,id,ir) = wtxfcl(ip+isfcl,id,ir) / totevt

            end do
            end do
            end do

*-----------------------------------------------------------------------

                     igm = ( mmmax - 1 ) * 2 + 1

         do kk = 1, ifcln

                     kimps = 0

         do k = 1, mnfcl(kk,20)

                     lngmax = 0
                     jsek   = 0
                     t2     = 0.0d0

                     ii = 0

                     idsm = inflt(kk)
                     jdsm = 0

               do m = 1, mnfcl(kk,0)

                     jdsm = jdsm + 1
                     ntrn = idas_inflt(idsm+jdsm)
                     jdsm = jdsm + 1
                     mtrn = idas_inflt(idsm+jdsm)
                     kssm = idsm + jdsm + 1
                     jdsm = jdsm + mtrn

                     jj = 0

                  do j = 1, ntrn

                     ii = ii + 1

                     call echrg3(jj,mtrn,idas_inflt(kssm)
     &                           ,dum1,lng1,icmb,igm)

                     if( wtxfcl(mnfcl(kk,k)+isfcl,2,ii) .gt. 0.d0 ) then

                         jsek = jsek + 1

                         if( lng1 .gt. lngmax ) lngmax = lng1

                         t2 = t2 + wtxfcl(mnfcl(kk,k)+isfcl,2,ii)

                     end if

                  end do

               end do

*-----------------------------------------------------------------------

            if( jsek .ne. 0 ) then

               if( kimps .eq. 0 ) then

                  write(io,'(/79(''-''))')
                  write(io,'('' Summary of forced collisions'',
     &                    '' per source (non zero)''/
     &                    '' for '',i2,''-th [forced collisions]'')') kk
                  write(io,'( 79(''-''))')

                  kimps = kimps + 1

               end if

                  uname(k) = pname(mnfcl(kk,k))(1:8)

*-----------------------------------------------------------------------

                     lblk = lngmax
                     dum2(1:5) = '  reg'

                   do i = 1, lblk
                     dum2(i+5:i+5) = ' '
                   end do

                     dum2(lblk+5+1:lblk+5+35) =
     &              '[ '// uname(k) //
     &              ']  forced col  forced col'

                     write(io,'(/600a1)') (dum2(i:i),i=1,lblk+5+35)

                     dum2(1:5) = '     '
                   do i = 1, lblk
                     dum2(i+5:i+5) = ' '
                   end do

                     dum2(lblk+5+1:lblk+5+35) =
     &              ' fcl factor      weight      events'

                     write(io,'(600a1)') (dum2(i:i),i=1,lblk+5+35)

                     ii = 0

                     idsm = inflt(kk)
                     jdsm = 0

                     kdsm = kfcls(kk)

               do m = 1, mnfcl(kk,0)

                     jdsm = jdsm + 1
                     ntrn = idas_inflt(idsm+jdsm)
                     jdsm = jdsm + 1
                     mtrn = idas_inflt(idsm+jdsm)
                     kssm = idsm + jdsm + 1
                     jdsm = jdsm + mtrn

                     rfcls = das_kfcls(kdsm-1+m)

                     jj = 0

                  do j = 1, ntrn

                     ii = ii + 1

                     call echrg3(jj,mtrn,idas_inflt(kssm)
     &                           ,dum1,lng1,icmb,igm)

                  if( wtxfcl(mnfcl(kk,k)+isfcl,2,ii) .gt. 0.0 ) then

                     dum2(1:3) = '   '

                     do i = 1, lng1
                        dum2(i+3:i+3) = dum1(i:i)
                     end do

                     lngb = lngmax - lng1

                     do i = 1, lngb
                        dum2(i+lng1+3:i+lng1+3) = ' '
                     end do

                     write(dum2(lngb+lng1+4:lngb+lng1+41),
     &                         '(1x,1p3e12.4)') rfcls,
     &                         wtxfcl(mnfcl(kk,k)+isfcl,1,ii),
     &                         wtxfcl(mnfcl(kk,k)+isfcl,2,ii)

                     write(io,'(600a1)') (dum2(i:i),i=1,lngb+lng1+41)

                  end if
                  end do

               end do

                     dum2(1:13) = ' total events'

                     lstl = 3+lngmax+25

                  do i = 14, lstl

                     dum2(i:i) = ' '

                  end do

                     write(dum2(lstl+1:lstl+13),
     &                         '(1p1e12.4)') t2

                     write(io,'(600a1)') (dum2(i:i),i=1,lstl+13)

            end if

         end do
         end do

      end if
      end if

*-----------------------------------------------------------------------
*     region-wise total number of collision
*-----------------------------------------------------------------------

         if( npe .gt. 1 ) then

          allocate( aevt1(inevt*iregn) ) !FURUTA20210623

            if( me .gt. 0 ) then

!-----------------------------------------------------------------------
             do m=1,iregn
              do n=1,inevt
               aevt1(n+(m-1)*inevt)=aevts(iaevt+n,m)
              enddo
             enddo
             call parasr(aevt1(1),inevt*iregn,0)
!------------------------------------------------------------------------

            else

                     do n = 1, inevt
                     do m = 1, iregn

                        aevts(iaevt+n,m) = 0.0

                     end do
                     end do

               do i = 1, npe - 1

                  if( iccp(i) .eq. 0 ) then

!------------------------------------------------------------------------
                   call pararr(aevt1(1),inevt*iregn,i)

                   do m = 1, iregn
                    do n = 1, inevt
                     aevts(iaevt+n,m) = aevts(iaevt+n,m)
     &                    + aevt1(n+(m-1)*inevt)
                    end do
                   end do
!------------------------------------------------------------------------

                  end if

               end do

            end if

            deallocate( aevt1 ) !FURUTA20210623

         end if

      if( me .eq. 0 ) then

*-----------------------------------------------------------------------
*     Region-wise total number of collisions
*     for High energy part and  HI collisions
*-----------------------------------------------------------------------

               seky = 0.0
               sekc = 0.0
               seke = 0.0
               sekn = 0.0
               sekh = 0.0
               sekd = 0.0

            do m = 1, iregn

               seky = seky+ aevts(iaevt+1,m)
               sekc = sekc+ aevts(iaevt+2,m)
               seke = seke+ aevts(iaevt+3,m)
               sekn = sekn+ aevts(iaevt+4,m)

               sekh = sekh
     &              + aevts(iaevt+5,m) + aevts(iaevt+6,m)
     &              + aevts(iaevt+7,m) + aevts(iaevt+8,m)
     &              + aevts(iaevt+9,m)

            end do

         if( seky+sekc+seke+sekn+sekh .gt. 0.0 ) then

               write(io,'(/79(''-''))')
               write(io,'('' Region-wise total number of collisions'',
     &                    '' for High energy and HI per source''
     &                    /79(''-'')/
     &                    '' Region'',
     &                    ''    Hydrogen'',
     &                    ''       Decay'',
     &                    ''     Elastic'',
     &                    ''   Non-Elast'',
     &                    ''   Heavy Ion'',
     &                    ''    Hi-Total'')')

            do m = 1, iregn

               if( aevts(iaevt+ 1,m)+aevts(iaevt+ 2,m)
     &            +aevts(iaevt+ 3,m)+aevts(iaevt+ 4,m)
     &            +aevts(iaevt+ 5,m)+aevts(iaevt+ 6,m)
     &            +aevts(iaevt+ 7,m)+aevts(iaevt+ 8,m)
     &            +aevts(iaevt+ 9,m) .gt. 0.0 )
     &               write(io,'(i7,1p6e12.4)')
     &               idrg(m),aevts(iaevt+ 1,m) / totevt,
     &                       aevts(iaevt+ 2,m) / totevt,
     &                       aevts(iaevt+ 3,m) / totevt,
     &                       aevts(iaevt+ 4,m) / totevt,
     &                     ( aevts(iaevt+ 5,m) + aevts(iaevt+6,m)
     &                     + aevts(iaevt+ 7,m) + aevts(iaevt+8,m)
     &                     + aevts(iaevt+ 9,m) ) / totevt,
     &                     ( aevts(iaevt+ 1,m) + aevts(iaevt+2,m)
     &                     + aevts(iaevt+ 3,m) + aevts(iaevt+4,m)
     &                     + aevts(iaevt+ 5,m) + aevts(iaevt+6,m)
     &                     + aevts(iaevt+ 7,m) + aevts(iaevt+8,m)
     &                     + aevts(iaevt+ 9,m) ) / totevt

            end do

               write(io,'(''  Total'',1p6e12.4)')
     &         seky/totevt,sekc/totevt,seke/totevt,sekn/totevt,
     &         sekh/totevt,(seky+sekc+seke+sekn+sekh)/totevt

         end if

*-----------------------------------------------------------------------
*     Region-wise total number of collisions
*     for low energy part per source
*-----------------------------------------------------------------------

               seky = 0.0
               sekc = 0.0
               seke = 0.0
               sekf = 0.0
               sekg = 0.0
               sekd = 0.0

            do m = 1, iregn

               seky = seky+ aevts(iaevt+10,m)
               sekc = sekc+ aevts(iaevt+11,m)
               seke = seke+ aevts(iaevt+12,m)
               sekf = sekf+ aevts(iaevt+57,m)
               sekg = sekg+ aevts(iaevt+63,m)
               sekd = sekd
     &              + aevts(iaevt+10,m) + aevts(iaevt+11,m)
     &              + aevts(iaevt+12,m) + aevts(iaevt+57,m)
     &              + aevts(iaevt+63,m)

            end do

         if( sekd .gt. 0.0 ) then

               write(io,'(/79(''-''))')
               write(io,'('' Region-wise total number of collisions'',
     &                    '' for Nuclear Data part per source''
     &                    /79(''-'')/
     &                    '' Region'',
     &                    ''     Neutron'',
     &                    ''      Photon'',
     &                    ''    Electron'',
     &                    ''      Proton'',
     &                    ''     N-Emode'',
     &                    ''  Data-Total'')')

            do m = 1, iregn

               if( aevts(iaevt+10,m)+aevts(iaevt+11,m)
     &            +aevts(iaevt+12,m)+aevts(iaevt+57,m)
     &            +aevts(iaevt+63,m) .gt. 0.0 )
     &               write(io,'(i7,1p6e12.4)')
     &               idrg(m),aevts(iaevt+10,m) / totevt,
     &                       aevts(iaevt+11,m) / totevt,
     &                       aevts(iaevt+12,m) / totevt,
     &                       aevts(iaevt+57,m) / totevt,
     &                       aevts(iaevt+63,m) / totevt,
     &                     ( aevts(iaevt+10,m) + aevts(iaevt+11,m)
     &                     + aevts(iaevt+12,m) + aevts(iaevt+57,m)
     &                     + aevts(iaevt+63,m) ) / totevt

            end do

               write(io,'(''  Total'',1p6e12.4)')
     &         seky/totevt,sekc/totevt,seke/totevt,sekf/totevt,
     &         sekg/totevt,sekd/totevt

         end if

*-----------------------------------------------------------------------
*      Region-wise total number of n-coll
*      for low energy per source
*-----------------------------------------------------------------------

               sekg = 0.0
               seke = 0.0
               sekn = 0.0
               sekf = 0.0
               sekd = 0.0

            do m = 1, iregn

               sekg = sekg + aevts(iaevt+15,m)
               seke = seke + aevts(iaevt+17,m)
               sekn = sekn + aevts(iaevt+18,m)
               sekf = sekf + aevts(iaevt+19,m)
               sekd = sekd + aevts(iaevt+24,m)

            end do

         if( sekg+seke+sekn+sekf+sekd .gt. 0 ) then

               write(io,'(/79(''-''))')
               write(io,'('' Region-wise total number of n-coll'',
     &                    '' for low energy per source''
     &                    /79(''-'')/
     &                    '' Region'',
     &                    ''     Capture'',
     &                    ''     Elastic'',
     &                    ''      (N,N'',a1,'')'',
     &                    ''      (N,xN)'',
     &                    ''     Fission'',
     &                    ''  Total-Coll'')') "'"

            do m = 1, iregn

               if( aevts(iaevt+15,m)
     &            +aevts(iaevt+17,m)+aevts(iaevt+18,m)
     &            +aevts(iaevt+19,m)+aevts(iaevt+24,m) .gt. 0.0 )
     &               write(io,'(i7,1p6e12.4)')
     &               idrg(m),aevts(iaevt+15,m)/totevt,
     &                       aevts(iaevt+17,m)/totevt,
     &                       aevts(iaevt+24,m)/totevt,
     &                       aevts(iaevt+18,m)/totevt,
     &                       aevts(iaevt+19,m)/totevt,
     &                       aevts(iaevt+17,m)/totevt+
     &                       aevts(iaevt+18,m)/totevt+
     &                       aevts(iaevt+19,m)/totevt+
     &                       aevts(iaevt+24,m)/totevt

            end do

               write(io,'(''  Total'',1p6e12.4)')
     &         sekg/totevt,
     &         seke/totevt,sekd/totevt,sekn/totevt,sekf/totevt,
     &         sekg/totevt+
     &         seke/totevt+sekn/totevt+sekf/totevt+sekd/totevt

         end if

*-----------------------------------------------------------------------
*     Region-wise neutron flux gain and loss
*     for low energy per source
*-----------------------------------------------------------------------

               sekc = 0.0
               sekn = 0.0
               sekf = 0.0
               sekd = 0.0
               sekg = 0.0

            do m = 1, iregn

               sekc = sekc + aevts(iaevt+15,m)
               sekn = sekn + aevts(iaevt+18,m)
               sekf = sekf + aevts(iaevt+20,m)
               sekd = sekd + aevts(iaevt+19,m)
               sekg = sekg + aevts(iaevt+21,m)

            end do

         if( sekc+sekn+sekf+sekd+sekg .gt. 0 ) then

               write(io,'(/79(''-''))')
               write(io,'('' Region-wise neutron flux gain and loss'',
     &                    '' for low energy per source''
     &                    /79(''-'')/
     &                    '' Region'',
     &                    ''   Capt-Loss'',
     &                    ''  (N,xN)Loss'',
     &                    ''  (N,xN)Gain'',
     &                    ''   Fiss-Loss'',
     &                    ''   Fiss-Gain'',
     &                    ''       Total'')')

            do m = 1, iregn

               if( aevts(iaevt+15,m)+aevts(iaevt+18,m)
     &            +aevts(iaevt+20,m)+aevts(iaevt+19,m)
     &            +aevts(iaevt+21,m) .gt. 0.0 )
     &               write(io,'(i7,1p6e12.4)')
     &               idrg(m),-aevts(iaevt+15,m)/totevt,
     &                       -aevts(iaevt+18,m)/totevt,
     &                       +aevts(iaevt+20,m)/totevt,
     &                       -aevts(iaevt+19,m)/totevt,
     &                       +aevts(iaevt+21,m)/totevt,
     &                       -aevts(iaevt+15,m)/totevt
     &                       -aevts(iaevt+18,m)/totevt
     &                       +aevts(iaevt+20,m)/totevt
     &                       -aevts(iaevt+19,m)/totevt
     &                       +aevts(iaevt+21,m)/totevt

            end do

               write(io,'(''  Total'',1p6e12.4)')
     &         -sekc/totevt,-sekn/totevt,+sekf/totevt,
     &         -sekd/totevt,+sekg/totevt,
     &         -sekc/totevt -sekn/totevt +sekf/totevt
     &         -sekd/totevt +sekg/totevt

         end if

*-----------------------------------------------------------------------
*      Region-wise neutron energy gain and loss
*      for low energy per source
*-----------------------------------------------------------------------

               sekc = 0.0
               seku = 0.0
               sekd = 0.0

            do m = 1, iregn

               sekc = sekc + aevts(iaevt+16,m)
               seku = seku + aevts(iaevt+22,m)
               sekd = sekd + aevts(iaevt+23,m)

            end do

         if( sekg+sekc+seku+sekd .gt. 0 ) then

               write(io,'(/79(''-''))')
               write(io,'('' Region-wise neutron energy gain and loss'',
     &                    '' for low energy per source''
     &                    /79(''-'')/
     &                    '' Region'',
     &                    ''   Capt-Loss'',
     &                    ''   UpScatter'',
     &                    ''   DownScatt'')')

            do m = 1, iregn

               if( aevts(iaevt+16,m)
     &            +aevts(iaevt+22,m)+aevts(iaevt+23,m) .gt. 0.0 )
     &               write(io,'(i7,1p6e12.4)')
     &               idrg(m),aevts(iaevt+16,m)/totevt,
     &                       aevts(iaevt+22,m)/totevt,
     &                       aevts(iaevt+23,m)/totevt

            end do

               write(io,'(''  Total'',1p6e12.4)')
     &         sekc/totevt,seku/totevt,
     &         sekd/totevt

         end if

*-----------------------------------------------------------------------
*     Region-wise photon flux gain
*     for low energy per source
*-----------------------------------------------------------------------

               sekc = 0.0
               sekn = 0.0
               sekf = 0.0
               sekd = 0.0
               sekg = 0.0

            do m = 1, iregn

               sekc = sekc + aevts(iaevt+13,m)
               sekn = sekn + aevts(iaevt+26,m)
               sekf = sekf + aevts(iaevt+30,m)
               sekd = sekd + aevts(iaevt+32,m)
               sekg = sekg + aevts(iaevt+34,m)

            end do

         if( sekc+sekn+sekf+sekd+sekg .gt. 0 ) then

               write(io,'(/79(''-''))')
               write(io,'('' Region-wise photon flux gain'',
     &                    '' for low energy per source''
     &                    /79(''-'')/
     &                    '' Region'',
     &                    ''   from Neut'',
     &                    ''      bremss'',
     &                    ''    p-annihi'',
     &                    ''   el x-rays'',
     &                    ''   fluoresnc'',
     &                    ''       Total'')')

            do m = 1, iregn

               if( aevts(iaevt+13,m)+aevts(iaevt+26,m)
     &            +aevts(iaevt+30,m)+aevts(iaevt+32,m)
     &            +aevts(iaevt+34,m) .gt. 0.0 )
     &               write(io,'(i7,1p6e12.4)')
     &               idrg(m),aevts(iaevt+13,m)/totevt,
     &                       aevts(iaevt+26,m)/totevt,
     &                       aevts(iaevt+30,m)/totevt,
     &                       aevts(iaevt+32,m)/totevt,
     &                       aevts(iaevt+34,m)/totevt,
     &                       aevts(iaevt+13,m)/totevt
     &                      +aevts(iaevt+26,m)/totevt
     &                      +aevts(iaevt+30,m)/totevt
     &                      +aevts(iaevt+32,m)/totevt
     &                      +aevts(iaevt+34,m)/totevt

            end do

               write(io,'(''  Total'',1p6e12.4)')
     &          sekc/totevt, sekn/totevt, sekf/totevt,
     &          sekd/totevt, sekg/totevt,
     &         (sekc+sekn+sekf+sekd+sekg)/totevt

         end if

*-----------------------------------------------------------------------
*     Region-wise photon energy gain
*     for low energy per source
*-----------------------------------------------------------------------

               sekc = 0.0
               sekn = 0.0
               sekf = 0.0
               sekd = 0.0
               sekg = 0.0

            do m = 1, iregn

               sekc = sekc + aevts(iaevt+14,m)
               sekn = sekn + aevts(iaevt+27,m)
               sekf = sekf + aevts(iaevt+31,m)
               sekd = sekd + aevts(iaevt+33,m)
               sekg = sekg + aevts(iaevt+35,m)

            end do

         if( sekc+sekn+sekf+sekd+sekg .gt. 0 ) then

               write(io,'(/79(''-''))')
               write(io,'('' Region-wise photon energy gain'',
     &                    '' for low energy per source''
     &                    /79(''-'')/
     &                    '' Region'',
     &                    ''   from Neut'',
     &                    ''      bremss'',
     &                    ''    p-annihi'',
     &                    ''   el x-rays'',
     &                    ''   fluoresnc'',
     &                    ''       Total'')')

            do m = 1, iregn

               if( aevts(iaevt+14,m)+aevts(iaevt+27,m)
     &            +aevts(iaevt+31,m)+aevts(iaevt+33,m)
     &            +aevts(iaevt+35,m) .gt. 0.0 )
     &               write(io,'(i7,1p6e12.4)')
     &               idrg(m),aevts(iaevt+14,m)/totevt,
     &                       aevts(iaevt+27,m)/totevt,
     &                       aevts(iaevt+31,m)/totevt,
     &                       aevts(iaevt+33,m)/totevt,
     &                       aevts(iaevt+35,m)/totevt,
     &                       aevts(iaevt+14,m)/totevt
     &                      +aevts(iaevt+27,m)/totevt
     &                      +aevts(iaevt+31,m)/totevt
     &                      +aevts(iaevt+33,m)/totevt
     &                      +aevts(iaevt+35,m)/totevt

            end do

               write(io,'(''  Total'',1p6e12.4)')
     &          sekc/totevt, sekn/totevt, sekf/totevt,
     &          sekd/totevt, sekg/totevt,
     &         (sekc+sekn+sekf+sekd+sekg)/totevt

         end if

*-----------------------------------------------------------------------
*     Region-wise photon flux loss
*     for low energy per source
*-----------------------------------------------------------------------

               sekc = 0.0
               sekn = 0.0

            do m = 1, iregn

               sekc = sekc + aevts(iaevt+28,m)
               sekn = sekn + aevts(iaevt+36,m)

            end do

         if( sekc+sekn .gt. 0 ) then

               write(io,'(/79(''-''))')
               write(io,'('' Region-wise photon flux loss'',
     &                    '' for low energy per source''
     &                    /79(''-'')/
     &                    '' Region'',
     &                    ''     capture'',
     &                    ''   pair-prod'',
     &                    ''       Total'')')

            do m = 1, iregn

               if( aevts(iaevt+28,m)+aevts(iaevt+36,m) .gt. 0.0 )
     &               write(io,'(i7,1p6e12.4)')
     &               idrg(m),aevts(iaevt+28,m)/totevt,
     &                       aevts(iaevt+36,m)/totevt,
     &                       aevts(iaevt+28,m)/totevt
     &                      +aevts(iaevt+36,m)/totevt

            end do

               write(io,'(''  Total'',1p6e12.4)')
     &          sekc/totevt, sekn/totevt,
     &         (sekc+sekn)/totevt

         end if

*-----------------------------------------------------------------------
*     Region-wise photon energy loss
*     for low energy per source
*-----------------------------------------------------------------------

               sekc = 0.0
               sekn = 0.0
               sekf = 0.0

            do m = 1, iregn

               sekc = sekc + aevts(iaevt+29,m)
               sekn = sekn + aevts(iaevt+37,m)
               sekf = sekf + aevts(iaevt+62,m)

            end do

         if( sekc+sekn .gt. 0 ) then

               write(io,'(/79(''-''))')
               write(io,'('' Region-wise photon energy loss'',
     &                    '' for low energy per source''
     &                    /79(''-'')/
     &                    '' Region'',
     &                    ''     capture'',
     &                    ''   pair-prod'',
     &                    ''     compton'',
     &                    ''       Total'')')

            do m = 1, iregn

               if( aevts(iaevt+29,m)+aevts(iaevt+37,m)+
     &             aevts(iaevt+62,m) .gt. 0.0 )
     &               write(io,'(i7,1p6e12.4)')
     &               idrg(m),aevts(iaevt+29,m)/totevt,
     &                       aevts(iaevt+37,m)/totevt,
     &                       aevts(iaevt+62,m)/totevt,
     &                       aevts(iaevt+29,m)/totevt
     &                      +aevts(iaevt+37,m)/totevt
     &                      +aevts(iaevt+62,m)/totevt

            end do

               write(io,'(''  Total'',1p6e12.4)')
     &          sekc/totevt, sekn/totevt, sekf/totevt,
     &         (sekc+sekn+sekf)/totevt

         end if

*-----------------------------------------------------------------------
*     Region-wise electron flux gain
*     for low energy per source
*-----------------------------------------------------------------------

               sekc = 0.0
               sekn = 0.0
               sekf = 0.0
               sekd = 0.0
               sekg = 0.0

            do m = 1, iregn

               sekc = sekc + aevts(iaevt+38,m)
               sekn = sekn + aevts(iaevt+40,m)
               sekf = sekf + aevts(iaevt+42,m)
               sekd = sekd + aevts(iaevt+44,m)
               sekg = sekg + aevts(iaevt+46,m)

            end do

         if( sekc+sekn+sekf+sekd+sekg .gt. 0 ) then

               write(io,'(/79(''-''))')
               write(io,'('' Region-wise electron flux gain'',
     &                    '' for low energy per source''
     &                    /79(''-'')/
     &                    '' Region'',
     &                    ''   pair prod'',
     &                    ''     compton'',
     &                    ''   photo-ele'',
     &                    ''       auger'',
     &                    ''    knock-on'',
     &                    ''       Total'')')

            do m = 1, iregn

               if( aevts(iaevt+38,m)+aevts(iaevt+40,m)
     &            +aevts(iaevt+42,m)+aevts(iaevt+44,m)
     &            +aevts(iaevt+46,m) .gt. 0.0 )
     &               write(io,'(i7,1p6e12.4)')
     &               idrg(m),aevts(iaevt+38,m)/totevt,
     &                       aevts(iaevt+40,m)/totevt,
     &                       aevts(iaevt+42,m)/totevt,
     &                       aevts(iaevt+44,m)/totevt,
     &                       aevts(iaevt+46,m)/totevt,
     &                       aevts(iaevt+38,m)/totevt
     &                      +aevts(iaevt+40,m)/totevt
     &                      +aevts(iaevt+42,m)/totevt
     &                      +aevts(iaevt+44,m)/totevt
     &                      +aevts(iaevt+46,m)/totevt

            end do

               write(io,'(''  Total'',1p6e12.4)')
     &          sekc/totevt, sekn/totevt, sekf/totevt,
     &          sekd/totevt, sekg/totevt,
     &         (sekc+sekn+sekf+sekd+sekg)/totevt

         end if

*-----------------------------------------------------------------------
*     Region-wise electron energy gain
*     for low energy per source
*-----------------------------------------------------------------------

               sekc = 0.0
               sekn = 0.0
               sekf = 0.0
               sekd = 0.0
               sekg = 0.0

            do m = 1, iregn

               sekc = sekc + aevts(iaevt+39,m)
               sekn = sekn + aevts(iaevt+41,m)
               sekf = sekf + aevts(iaevt+43,m)
               sekd = sekd + aevts(iaevt+45,m)
               sekg = sekg + aevts(iaevt+47,m)

            end do

         if( sekc+sekn+sekf+sekd+sekg .gt. 0 ) then

               write(io,'(/79(''-''))')
               write(io,'('' Region-wise electron energy gain'',
     &                    '' for low energy per source''
     &                    /79(''-'')/
     &                    '' Region'',
     &                    ''   pair prod'',
     &                    ''     compton'',
     &                    ''   photo-ele'',
     &                    ''       auger'',
     &                    ''    knock-on'',
     &                    ''       Total'')')

            do m = 1, iregn

               if( aevts(iaevt+39,m)+aevts(iaevt+41,m)
     &            +aevts(iaevt+43,m)+aevts(iaevt+45,m)
     &            +aevts(iaevt+47,m) .gt. 0.0 )
     &               write(io,'(i7,1p6e12.4)')
     &               idrg(m),aevts(iaevt+39,m)/totevt,
     &                       aevts(iaevt+41,m)/totevt,
     &                       aevts(iaevt+43,m)/totevt,
     &                       aevts(iaevt+45,m)/totevt,
     &                       aevts(iaevt+47,m)/totevt,
     &                       aevts(iaevt+39,m)/totevt
     &                      +aevts(iaevt+41,m)/totevt
     &                      +aevts(iaevt+43,m)/totevt
     &                      +aevts(iaevt+45,m)/totevt
     &                      +aevts(iaevt+47,m)/totevt

            end do

               write(io,'(''  Total'',1p6e12.4)')
     &          sekc/totevt, sekn/totevt, sekf/totevt,
     &          sekd/totevt, sekg/totevt,
     &         (sekc+sekn+sekf+sekd+sekg)/totevt

         end if

*-----------------------------------------------------------------------
*     Region-wise electron energy loss
*     for low energy per source
*-----------------------------------------------------------------------

               sekc = 0.0
               sekn = 0.0

            do m = 1, iregn

               sekc = sekc + aevts(iaevt+48,m)
               sekn = sekn + aevts(iaevt+49,m)

            end do

         if( sekc+sekn .gt. 0 ) then

               write(io,'(/79(''-''))')
               write(io,'('' Region-wise electron energy loss'',
     &                    '' for low energy per source''
     &                    /79(''-'')/
     &                    '' Region'',
     &                    ''     scatter'',
     &                    ''      bremss'',
     &                    ''       Total'')')

            do m = 1, iregn

               if( aevts(iaevt+48,m)+
     &             aevts(iaevt+49,m) .gt. 0.0 )
     &               write(io,'(i7,1p6e12.4)')
     &               idrg(m),aevts(iaevt+48,m)/totevt,
     &                       aevts(iaevt+49,m)/totevt,
     &                       aevts(iaevt+48,m)/totevt
     &                      +aevts(iaevt+49,m)/totevt

            end do

               write(io,'(''  Total'',1p6e12.4)')
     &          sekc/totevt, sekn/totevt,
     &         (sekc+sekn)/totevt

         end if

*-----------------------------------------------------------------------
*      Region-wise total number of n-coll
*      for high energy library per source
*-----------------------------------------------------------------------

               sekg = 0.0
               seke = 0.0
               sekn = 0.0

            do m = 1, iregn

               sekg = sekg + aevts(iaevt+50,m)
               seke = seke + aevts(iaevt+51,m)
               sekn = sekn + aevts(iaevt+52,m)

            end do

         if( sekg+seke+sekn .gt. 0 ) then

               write(io,'(/79(''-''))')
               write(io,'('' Region-wise total number of n-coll'',
     &                    '' for high energy library per source''
     &                    /79(''-'')/
     &                    '' Region'',
     &                    ''      (N,xP)'',
     &                    ''     (N,xPI)'',
     &                    ''   (N,other)'',
     &                    ''  Total-Coll'')')

            do m = 1, iregn

               if( aevts(iaevt+50,m)
     &            +aevts(iaevt+51,m)+aevts(iaevt+52,m) .gt. 0.0 )
     &               write(io,'(i7,1p6e12.4)')
     &               idrg(m),aevts(iaevt+50,m)/totevt,
     &                       aevts(iaevt+51,m)/totevt,
     &                       aevts(iaevt+52,m)/totevt,
     &                       aevts(iaevt+50,m)/totevt +
     &                       aevts(iaevt+51,m)/totevt +
     &                       aevts(iaevt+52,m)/totevt

            end do

               write(io,'(''  Total'',1p6e12.4)')
     &         sekg/totevt,
     &         seke/totevt,sekn/totevt,
     &         sekg/totevt+seke/totevt+sekn/totevt

         end if

*-----------------------------------------------------------------------
*      Region-wise total number of p-coll
*      for high energy library per source
*-----------------------------------------------------------------------

               sekg = 0.0
               seke = 0.0
               sekn = 0.0

            do m = 1, iregn

               sekg = sekg + aevts(iaevt+59,m)
               seke = seke + aevts(iaevt+60,m)
               sekn = sekn + aevts(iaevt+61,m)

            end do

         if( sekg+seke+sekn .gt. 0 ) then

               write(io,'(/79(''-''))')
               write(io,'('' Region-wise total number of p-coll'',
     &                    '' for high energy library per source''
     &                    /79(''-'')/
     &                    '' Region'',
     &                    ''     Capture'',
     &                    ''     Elastic'',
     &                    ''   Non-Elast'',
     &                    ''  Total-Coll'')')

            do m = 1, iregn

               if( aevts(iaevt+59,m)+aevts(iaevt+60,m)
     &            +aevts(iaevt+61,m) .gt. 0.0 )
     &               write(io,'(i7,1p6e12.4)')
     &               idrg(m),aevts(iaevt+59,m)/totevt,
     &                       aevts(iaevt+60,m)/totevt,
     &                       aevts(iaevt+61,m)/totevt,
     &                       aevts(iaevt+59,m)/totevt+
     &                       aevts(iaevt+60,m)/totevt+
     &                       aevts(iaevt+61,m)/totevt

            end do

               write(io,'(''  Total'',1p6e12.4)')
     &         sekg/totevt,seke/totevt,sekn/totevt,
     &         sekg/totevt+seke/totevt+sekn/totevt

         end if

*-----------------------------------------------------------------------
*      Region-wise total number of p-coll
*      for high energy library per source
*-----------------------------------------------------------------------

               sekg = 0.0
               seke = 0.0
               sekn = 0.0
               sekf = 0.0
               sekd = 0.0

            do m = 1, iregn

               sekg = sekg + aevts(iaevt+53,m)
               seke = seke + aevts(iaevt+54,m)
               sekn = sekn + aevts(iaevt+55,m)
               sekf = sekf + aevts(iaevt+56,m)
               sekd = sekd + aevts(iaevt+58,m)

            end do

         if( sekg+seke+sekn+sekf+sekd .gt. 0 ) then

               write(io,'(/79(''-''))')
               write(io,'('' Region-wise total number of p-coll'',
     &                    '' for high energy library per source''
     &                    /79(''-'')/
     &                    '' Region'',
     &                    ''      (P,xP)'',
     &                    ''      (P,xN)'',
     &                    ''  (P,Photon)'',
     &                    ''     (P,xPI)'',
     &                    ''   (P,other)'',
     &                    ''  Total-Coll'')')

            do m = 1, iregn

               if( aevts(iaevt+53,m)+aevts(iaevt+54,m)
     &            +aevts(iaevt+55,m)+aevts(iaevt+56,m)
     &            +aevts(iaevt+58,m) .gt. 0.0 )
     &               write(io,'(i7,1p6e12.4)')
     &               idrg(m),aevts(iaevt+53,m)/totevt,
     &                       aevts(iaevt+54,m)/totevt,
     &                       aevts(iaevt+58,m)/totevt,
     &                       aevts(iaevt+55,m)/totevt,
     &                       aevts(iaevt+56,m)/totevt,
     &                       aevts(iaevt+53,m)/totevt+
     &                       aevts(iaevt+54,m)/totevt+
     &                       aevts(iaevt+55,m)/totevt+
     &                       aevts(iaevt+56,m)/totevt+
     &                       aevts(iaevt+58,m)/totevt

            end do

               write(io,'(''  Total'',1p6e12.4)')
     &         sekg/totevt,
     &         seke/totevt,sekd/totevt,sekn/totevt,sekf/totevt,
     &         sekg/totevt+seke/totevt+
     &         sekd/totevt+sekn/totevt+sekf/totevt

         end if

cABE add @2014/08/12, to output muon reaction
*-----------------------------------------------------------------------
*     Region-wise total number of collisions
*     for Muon Interaction
*-----------------------------------------------------------------------

               seky = 0.d0
               sekb = 0.d0
               sekp = 0.d0
               sekc = 0.d0

            do m = 1, iregn
               seky = seky + aevts(iaevt+64,m)
               sekb = sekb + aevts(iaevt+66,m)
               sekp = sekp + aevts(iaevt+67,m)
               sekc = sekc + aevts(iaevt+65,m)
            end do

         if( seky+sekb+sekp+sekc .gt. 0.d0 ) then

               write(io,'(/79(''-''))')
               write(io,'('' Region-wise total number of collisions'',
     &                    '' for Muon Interaction per source''
     &                    /79(''-'')/
     &                    '' Region'',
     &                    ''   photonucl'',
     &                    ''      bremss'',
     &                    ''   pair prod'',
     &                    ''     capture'',
     &                    ''       Total'')')

            do m = 1, iregn

               if( aevts(iaevt+64,m)+aevts(iaevt+66,m)+aevts(iaevt+67,m)
     &            +aevts(iaevt+65,m) .gt. 0.0 )
     &               write(io,'(i7,1p5e12.4)')
     &               idrg(m),aevts(iaevt+64,m) / totevt,
     &                       aevts(iaevt+66,m) / totevt,
     &                       aevts(iaevt+67,m) / totevt,
     &                       aevts(iaevt+65,m) / totevt,
     &                       aevts(iaevt+64,m) / totevt+
     &                       aevts(iaevt+66,m) / totevt+
     &                       aevts(iaevt+67,m) / totevt+
     &                       aevts(iaevt+65,m) / totevt

            end do

               write(io,'(''  Total'',1p5e12.4)')
     &         seky/totevt,sekb/totevt,sekp/totevt,sekc/totevt,
     &         seky/totevt+sekb/totevt+sekp/totevt+sekc/totevt

         end if
*-----------------------------------------------------------------------
*     Region-wise total number of collisions
*     for Track structure mode
*-----------------------------------------------------------------------

               sekk = 0.0
               sekp = 0.0
               sekc = 0.0
               seka = 0.0

            do m = 1, iregn
               sekk = sekk + aevts(iaevt+68,m)
               sekp = sekp + aevts(iaevt+69,m)
               sekc = sekc + aevts(iaevt+70,m)
               seka = seka + aevts(iaevt+71,m)
            end do

         if( sekk+sekp+sekc+seka .gt. 0.d0 ) then

               write(io,'(/79(''-''))')
               write(io,'('' Region-wise total number of collisions'',
     &                    '' for Track structure mode per source''
     &                    /79(''-'')/
     &                    '' Region'',
     &                    ''         ETS'',
     &                    ''    KURBUC p'',
     &                    ''    KURBUC C'',
     &                    ''      ITSART'',
     &                    ''       Total'')')

            do m = 1, iregn

               if( aevts(iaevt+68,m)+aevts(iaevt+69,m)+aevts(iaevt+70,m)
     &            +aevts(iaevt+71,m) .gt. 0.0 )
     &               write(io,'(i7,1p5e12.4)')
     &               idrg(m),aevts(iaevt+68,m) / totevt,
     &                       aevts(iaevt+69,m) / totevt,
     &                       aevts(iaevt+70,m) / totevt,
     &                       aevts(iaevt+71,m) / totevt,
     &                       aevts(iaevt+68,m) / totevt+
     &                       aevts(iaevt+69,m) / totevt+
     &                       aevts(iaevt+70,m) / totevt+
     &                       aevts(iaevt+71,m) / totevt

            end do

               write(io,'(''  Total'',1p5e12.4)')
     &         sekk/totevt,sekp/totevt,sekc/totevt,seka/totevt,
     &         sekk/totevt+sekp/totevt+sekc/totevt+seka/totevt

         end if
*-----------------------------------------------------------------------

      end if

*-----------------------------------------------------------------------
*     medium-wise total number of collision
*-----------------------------------------------------------------------

         if( npe .gt. 1 ) then

          allocate( bevt1(inevt*mxmat) ) !FURUTA20210623

            if( me .gt. 0 ) then

!-----------------------------------------------------------------------
             do m=1,mxmat
              do n=1,inevt
               bevt1(n+(m-1)*inevt)=bevts(ibevt+n,m)
              enddo
             enddo
             call parasr(bevt1(1),inevt*mxmat,0)
!-----------------------------------------------------------------------

            else

                     do n = 1, inevt
                     do m = 1, mxmat

                        bevts(ibevt+n,m) = 0.0

                     end do
                     end do

               do i = 1, npe - 1

                  if( iccp(i) .eq. 0 ) then

!------------------------------------------------------------------------
                   call pararr(bevt1(1),inevt*mxmat,i)

                   do m = 1, mxmat
                    do n = 1, inevt
                     bevts(ibevt+n,m) = bevts(ibevt+n,m)
     &                    + bevt1(n+(m-1)*inevt)
                    end do
                   end do
!------------------------------------------------------------------------

                  end if

               end do

            end if

            deallocate( bevt1 ) !FURUTA20210623

         end if

      if( me .eq. 0 ) then

*-----------------------------------------------------------------------
*     Medium-wise total number of collisions
*     for High energy part and  HI collisions
*-----------------------------------------------------------------------

               seky = 0.0
               sekc = 0.0
               seke = 0.0
               sekn = 0.0
               sekh = 0.0
               sekd = 0.0

            do m = 1, mxmat

               seky = seky+ bevts(ibevt+1,m)
               sekc = sekc+ bevts(ibevt+2,m)
               seke = seke+ bevts(ibevt+3,m)
               sekn = sekn+ bevts(ibevt+4,m)

               sekh = sekh
     &              + bevts(ibevt+5,m) + bevts(ibevt+6,m)
     &              + bevts(ibevt+7,m) + bevts(ibevt+8,m)
     &              + bevts(ibevt+9,m)

            end do

         if( seky+sekc+seke+sekn+sekh .gt. 0.0 ) then

               write(io,'(/79(''-''))')
               write(io,'('' Medium-wise total number of collisions'',
     &                    '' for High energy and HI per source''
     &                    /79(''-'')/
     &                    '' Medium'',
     &                    ''    Hydrogen'',
     &                    ''       Decay'',
     &                    ''     Elastic'',
     &                    ''   Non-Elast'',
     &                    ''   Heavy Ion'',
     &                    ''    Hi-Total'')')

            do m = 1, mxmat

               if( bevts(ibevt+ 1,m)+bevts(ibevt+ 2,m)
     &            +bevts(ibevt+ 3,m)+bevts(ibevt+ 4,m)
     &            +bevts(ibevt+ 5,m)+bevts(ibevt+ 6,m)
     &            +bevts(ibevt+ 7,m)+bevts(ibevt+ 8,m)
     &            +bevts(ibevt+ 9,m) .gt. 0.0 )
     &               write(io,'(i7,1p6e12.4)')
     &               idmn(m),bevts(ibevt+ 1,m) / totevt,
     &                       bevts(ibevt+ 2,m) / totevt,
     &                       bevts(ibevt+ 3,m) / totevt,
     &                       bevts(ibevt+ 4,m) / totevt,
     &                     ( bevts(ibevt+ 5,m) + bevts(ibevt+6,m)
     &                     + bevts(ibevt+ 7,m) + bevts(ibevt+8,m)
     &                     + bevts(ibevt+ 9,m) ) / totevt,
     &                     ( bevts(ibevt+ 1,m) + bevts(ibevt+2,m)
     &                     + bevts(ibevt+ 3,m) + bevts(ibevt+4,m)
     &                     + bevts(ibevt+ 5,m) + bevts(ibevt+6,m)
     &                     + bevts(ibevt+ 7,m) + bevts(ibevt+8,m)
     &                     + bevts(ibevt+ 9,m) ) / totevt

            end do

               write(io,'(''  Total'',1p6e12.4)')
     &         seky/totevt,sekc/totevt,seke/totevt,sekn/totevt,
     &         sekh/totevt,(seky+sekc+seke+sekn+sekh)/totevt

         end if

*-----------------------------------------------------------------------
*     Medium-wise total number of collisions
*     for low energy part per source
*-----------------------------------------------------------------------

               seky = 0.0
               sekc = 0.0
               seke = 0.0
               sekf = 0.0
               sekd = 0.0

            do m = 1, mxmat

               seky = seky+ bevts(ibevt+10,m)
               sekc = sekc+ bevts(ibevt+11,m)
               seke = seke+ bevts(ibevt+12,m)
               sekf = sekf+ bevts(ibevt+57,m)
               sekd = sekd
     &              + bevts(ibevt+10,m) + bevts(ibevt+11,m)
     &              + bevts(ibevt+12,m) + bevts(ibevt+57,m)

            end do

         if( sekd .gt. 0.0 ) then

               write(io,'(/79(''-''))')
               write(io,'('' Medium-wise total number of collisions'',
     &                    '' for Nuclear Data part per source''
     &                    /79(''-'')/
     &                    '' Medium'',
     &                    ''     Neutron'',
     &                    ''      Photon'',
     &                    ''    Electron'',
     &                    ''      Proton'',
     &                    ''  Data-Total'')')

            do m = 1, mxmat

               if( bevts(ibevt+10,m)
     &            +bevts(ibevt+11,m)+bevts(ibevt+12,m) .gt. 0.0 )
     &               write(io,'(i7,1p6e12.4)')
     &               idmn(m),bevts(ibevt+10,m) / totevt,
     &                       bevts(ibevt+11,m) / totevt,
     &                       bevts(ibevt+12,m) / totevt,
     &                       bevts(ibevt+57,m) / totevt,
     &                     ( bevts(ibevt+10,m) + bevts(ibevt+11,m)
     &                     + bevts(ibevt+12,m)
     &                     + bevts(ibevt+57,m) ) / totevt

            end do

               write(io,'(''  Total'',1p6e12.4)')
     &         seky/totevt,sekc/totevt,seke/totevt,sekf/totevt,
     &         sekd/totevt

         end if

*-----------------------------------------------------------------------
*      Medium-wise total number of n-coll
*      for low energy per source
*-----------------------------------------------------------------------

               sekg = 0.0
               seke = 0.0
               sekn = 0.0
               sekf = 0.0
               sekd = 0.0

            do m = 1, mxmat

               sekg = sekg + bevts(ibevt+15,m)
               seke = seke + bevts(ibevt+17,m)
               sekn = sekn + bevts(ibevt+18,m)
               sekf = sekf + bevts(ibevt+19,m)
               sekd = sekd + bevts(ibevt+24,m)

            end do

         if( sekg+seke+sekn+sekf+sekd .gt. 0 ) then

               write(io,'(/79(''-''))')
               write(io,'('' Medium-wise total number of n-coll'',
     &                    '' for low energy per source''
     &                    /79(''-'')/
     &                    '' Medium'',
     &                    ''     Capture'',
     &                    ''     Elastic'',
     &                    ''      (N,N'',a1,'')'',
     &                    ''      (N,xN)'',
     &                    ''     Fission'',
     &                    ''  Total-Coll'')') "'"

            do m = 1, mxmat

               if( bevts(ibevt+15,m)
     &            +bevts(ibevt+17,m)+bevts(ibevt+18,m)
     &            +bevts(ibevt+19,m)+bevts(ibevt+24,m) .gt. 0.0 )
     &               write(io,'(i7,1p6e12.4)')
     &               idmn(m),bevts(ibevt+15,m)/totevt,
     &                       bevts(ibevt+17,m)/totevt,
     &                       bevts(ibevt+24,m)/totevt,
     &                       bevts(ibevt+18,m)/totevt,
     &                       bevts(ibevt+19,m)/totevt,
     &                       bevts(ibevt+17,m)/totevt+
     &                       bevts(ibevt+18,m)/totevt+
     &                       bevts(ibevt+19,m)/totevt+
     &                       bevts(ibevt+24,m)/totevt

            end do

               write(io,'(''  Total'',1p6e12.4)')
     &         sekg/totevt,
     &         seke/totevt,sekd/totevt,sekn/totevt,sekf/totevt,
     &         sekg/totevt+
     &         seke/totevt+sekn/totevt+sekf/totevt+sekd/totevt

         end if

*-----------------------------------------------------------------------
*     Medium-wise neutron flux gain and loss
*     for low energy per source
*-----------------------------------------------------------------------

               sekc = 0.0
               sekn = 0.0
               sekf = 0.0
               sekd = 0.0
               sekg = 0.0

            do m = 1, mxmat

               sekc = sekc + bevts(ibevt+15,m)
               sekn = sekn + bevts(ibevt+18,m)
               sekf = sekf + bevts(ibevt+20,m)
               sekd = sekd + bevts(ibevt+19,m)
               sekg = sekg + bevts(ibevt+21,m)

            end do

         if( sekc+sekn+sekf+sekd+sekg .gt. 0 ) then

               write(io,'(/79(''-''))')
               write(io,'('' Medium-wise neutron flux gain and loss'',
     &                    '' for low energy per source''
     &                    /79(''-'')/
     &                    '' Medium'',
     &                    ''   Capt-Loss'',
     &                    ''  (N,xN)Loss'',
     &                    ''  (N,xN)Gain'',
     &                    ''   Fiss-Loss'',
     &                    ''   Fiss-Gain'',
     &                    ''       Total'')')

            do m = 1, mxmat

               if( bevts(ibevt+15,m)+bevts(ibevt+18,m)
     &            +bevts(ibevt+20,m)+bevts(ibevt+19,m)
     &            +bevts(ibevt+21,m) .gt. 0.0 )
     &               write(io,'(i7,1p6e12.4)')
     &               idmn(m),-bevts(ibevt+15,m)/totevt,
     &                       -bevts(ibevt+18,m)/totevt,
     &                       +bevts(ibevt+20,m)/totevt,
     &                       -bevts(ibevt+19,m)/totevt,
     &                       +bevts(ibevt+21,m)/totevt,
     &                       -bevts(ibevt+15,m)/totevt
     &                       -bevts(ibevt+18,m)/totevt
     &                       +bevts(ibevt+20,m)/totevt
     &                       -bevts(ibevt+19,m)/totevt
     &                       +bevts(ibevt+21,m)/totevt

            end do

               write(io,'(''  Total'',1p6e12.4)')
     &         -sekc/totevt,-sekn/totevt,+sekf/totevt,
     &         -sekd/totevt,+sekg/totevt,
     &         -sekc/totevt -sekn/totevt +sekf/totevt
     &         -sekd/totevt +sekg/totevt

         end if

*-----------------------------------------------------------------------
*      Medium-wise neutron energy gain and loss
*      for low energy per source
*-----------------------------------------------------------------------

               sekc = 0.0
               seku = 0.0
               sekd = 0.0

            do m = 1, mxmat

               sekc = sekc + bevts(ibevt+16,m)
               seku = seku + bevts(ibevt+22,m)
               sekd = sekd + bevts(ibevt+23,m)

            end do

         if( sekg+sekc+seku+sekd .gt. 0 ) then

               write(io,'(/79(''-''))')
               write(io,'('' Medium-wise neutron energy gain and loss'',
     &                    '' for low energy per source''
     &                    /79(''-'')/
     &                    '' Medium'',
     &                    ''   Capt-Loss'',
     &                    ''   UpScatter'',
     &                    ''   DownScatt'')')

            do m = 1, mxmat

               if( bevts(ibevt+16,m)
     &            +bevts(ibevt+24,m)+bevts(ibevt+25,m) .gt. 0.0 )
     &               write(io,'(i7,1p6e12.4)')
     &               idmn(m),bevts(ibevt+16,m)/totevt,
     &                       bevts(ibevt+22,m)/totevt,
     &                       bevts(ibevt+23,m)/totevt

            end do

               write(io,'(''  Total'',1p6e12.4)')
     &         sekc/totevt,seku/totevt,
     &         sekd/totevt

         end if

*-----------------------------------------------------------------------
*     Medium-wise photon flux gain
*     for low energy per source
*-----------------------------------------------------------------------

               sekc = 0.0
               sekn = 0.0
               sekf = 0.0
               sekd = 0.0
               sekg = 0.0

            do m = 1, mxmat

               sekc = sekc + bevts(ibevt+13,m)
               sekn = sekn + bevts(ibevt+26,m)
               sekf = sekf + bevts(ibevt+30,m)
               sekd = sekd + bevts(ibevt+32,m)
               sekg = sekg + bevts(ibevt+34,m)

            end do

         if( sekc+sekn+sekf+sekd+sekg .gt. 0 ) then

               write(io,'(/79(''-''))')
               write(io,'('' Medium-wise photon flux gain'',
     &                    '' for low energy per source''
     &                    /79(''-'')/
     &                    '' Medium'',
     &                    ''   from Neut'',
     &                    ''      bremss'',
     &                    ''    p-annihi'',
     &                    ''   el x-rays'',
     &                    ''   fluoresnc'',
     &                    ''       Total'')')

            do m = 1, mxmat

               if( bevts(ibevt+13,m)+bevts(ibevt+26,m)
     &            +bevts(ibevt+30,m)+bevts(ibevt+32,m)
     &            +bevts(ibevt+34,m) .gt. 0.0 )
     &               write(io,'(i7,1p6e12.4)')
     &               idmn(m),bevts(ibevt+13,m)/totevt,
     &                       bevts(ibevt+26,m)/totevt,
     &                       bevts(ibevt+30,m)/totevt,
     &                       bevts(ibevt+32,m)/totevt,
     &                       bevts(ibevt+34,m)/totevt,
     &                       bevts(ibevt+13,m)/totevt
     &                      +bevts(ibevt+26,m)/totevt
     &                      +bevts(ibevt+30,m)/totevt
     &                      +bevts(ibevt+32,m)/totevt
     &                      +bevts(ibevt+34,m)/totevt

            end do

               write(io,'(''  Total'',1p6e12.4)')
     &          sekc/totevt, sekn/totevt, sekf/totevt,
     &          sekd/totevt, sekg/totevt,
     &         (sekc+sekn+sekf+sekd+sekg)/totevt

         end if

*-----------------------------------------------------------------------
*     Medium-wise photon energy gain
*     for low energy per source
*-----------------------------------------------------------------------

               sekc = 0.0
               sekn = 0.0
               sekf = 0.0
               sekd = 0.0
               sekg = 0.0

            do m = 1, mxmat

               sekc = sekc + bevts(ibevt+14,m)
               sekn = sekn + bevts(ibevt+27,m)
               sekf = sekf + bevts(ibevt+31,m)
               sekd = sekd + bevts(ibevt+33,m)
               sekg = sekg + bevts(ibevt+35,m)

            end do

         if( sekc+sekn+sekf+sekd+sekg .gt. 0 ) then

               write(io,'(/79(''-''))')
               write(io,'('' Medium-wise photon energy gain'',
     &                    '' for low energy per source''
     &                    /79(''-'')/
     &                    '' Medium'',
     &                    ''   from Neut'',
     &                    ''      bremss'',
     &                    ''    p-annihi'',
     &                    ''   el x-rays'',
     &                    ''   fluoresnc'',
     &                    ''       Total'')')

            do m = 1, mxmat

               if( bevts(ibevt+14,m)+bevts(ibevt+27,m)
     &            +bevts(ibevt+31,m)+bevts(ibevt+33,m)
     &            +bevts(ibevt+35,m) .gt. 0.0 )
     &               write(io,'(i7,1p6e12.4)')
     &               idmn(m),bevts(ibevt+14,m)/totevt,
     &                       bevts(ibevt+27,m)/totevt,
     &                       bevts(ibevt+31,m)/totevt,
     &                       bevts(ibevt+33,m)/totevt,
     &                       bevts(ibevt+35,m)/totevt,
     &                       bevts(ibevt+14,m)/totevt
     &                      +bevts(ibevt+27,m)/totevt
     &                      +bevts(ibevt+31,m)/totevt
     &                      +bevts(ibevt+33,m)/totevt
     &                      +bevts(ibevt+35,m)/totevt

            end do

               write(io,'(''  Total'',1p6e12.4)')
     &          sekc/totevt, sekn/totevt, sekf/totevt,
     &          sekd/totevt, sekg/totevt,
     &         (sekc+sekn+sekf+sekd+sekg)/totevt

         end if

*-----------------------------------------------------------------------
*     Medium-wise photon flux loss
*     for low energy per source
*-----------------------------------------------------------------------

               sekc = 0.0
               sekn = 0.0

            do m = 1, mxmat

               sekc = sekc + bevts(ibevt+28,m)
               sekn = sekn + bevts(ibevt+36,m)

            end do

         if( sekc+sekn .gt. 0 ) then

               write(io,'(/79(''-''))')
               write(io,'('' Medium-wise photon flux loss'',
     &                    '' for low energy per source''
     &                    /79(''-'')/
     &                    '' Medium'',
     &                    ''     capture'',
     &                    ''   pair-prod'',
     &                    ''       Total'')')

            do m = 1, mxmat

               if( bevts(ibevt+28,m)+bevts(ibevt+36,m) .gt. 0.0 )
     &               write(io,'(i7,1p6e12.4)')
     &               idmn(m),bevts(ibevt+28,m)/totevt,
     &                       bevts(ibevt+36,m)/totevt,
     &                       bevts(ibevt+28,m)/totevt
     &                      +bevts(ibevt+36,m)/totevt

            end do

               write(io,'(''  Total'',1p6e12.4)')
     &          sekc/totevt, sekn/totevt,
     &         (sekc+sekn)/totevt

         end if

*-----------------------------------------------------------------------
*     Medium-wise photon energy loss
*     for low energy per source
*-----------------------------------------------------------------------

               sekc = 0.0
               sekn = 0.0
               sekf = 0.0

            do m = 1, mxmat

               sekc = sekc + bevts(ibevt+29,m)
               sekn = sekn + bevts(ibevt+37,m)
               sekf = sekf + bevts(ibevt+62,m)

            end do

         if( sekc+sekn .gt. 0 ) then

               write(io,'(/79(''-''))')
               write(io,'('' Medium-wise photon energy loss'',
     &                    '' for low energy per source''
     &                    /79(''-'')/
     &                    '' Medium'',
     &                    ''     capture'',
     &                    ''   pair-prod'',
     &                    ''     compton'',
     &                    ''       Total'')')

            do m = 1, mxmat

               if( bevts(ibevt+29,m)+bevts(ibevt+37,m)+
     &             bevts(ibevt+62,m) .gt. 0.0 )
     &               write(io,'(i7,1p6e12.4)')
     &               idmn(m),bevts(ibevt+29,m)/totevt,
     &                       bevts(ibevt+37,m)/totevt,
     &                       bevts(ibevt+62,m)/totevt,
     &                       bevts(ibevt+29,m)/totevt
     &                      +bevts(ibevt+37,m)/totevt
     &                      +bevts(ibevt+62,m)/totevt
            end do

               write(io,'(''  Total'',1p6e12.4)')
     &          sekc/totevt, sekn/totevt, sekf/totevt,
     &         (sekc+sekn+sekf)/totevt

         end if

*-----------------------------------------------------------------------
*     Medium-wise electron flux gain
*     for low energy per source
*-----------------------------------------------------------------------

               sekc = 0.0
               sekn = 0.0
               sekf = 0.0
               sekd = 0.0
               sekg = 0.0

            do m = 1, mxmat

               sekc = sekc + bevts(ibevt+38,m)
               sekn = sekn + bevts(ibevt+40,m)
               sekf = sekf + bevts(ibevt+42,m)
               sekd = sekd + bevts(ibevt+44,m)
               sekg = sekg + bevts(ibevt+46,m)

            end do

         if( sekc+sekn+sekf+sekd+sekg .gt. 0 ) then

               write(io,'(/79(''-''))')
               write(io,'('' Medium-wise electron flux gain'',
     &                    '' for low energy per source''
     &                    /79(''-'')/
     &                    '' Medium'',
     &                    ''   pair prod'',
     &                    ''     compton'',
     &                    ''   photo-ele'',
     &                    ''       auger'',
     &                    ''    knock-on'',
     &                    ''       Total'')')

            do m = 1, mxmat

               if( bevts(ibevt+38,m)+bevts(ibevt+40,m)
     &            +bevts(ibevt+42,m)+bevts(ibevt+44,m)
     &            +bevts(ibevt+46,m) .gt. 0.0 )
     &               write(io,'(i7,1p6e12.4)')
     &               idmn(m),bevts(ibevt+38,m)/totevt,
     &                       bevts(ibevt+40,m)/totevt,
     &                       bevts(ibevt+42,m)/totevt,
     &                       bevts(ibevt+44,m)/totevt,
     &                       bevts(ibevt+46,m)/totevt,
     &                       bevts(ibevt+38,m)/totevt
     &                      +bevts(ibevt+40,m)/totevt
     &                      +bevts(ibevt+42,m)/totevt
     &                      +bevts(ibevt+44,m)/totevt
     &                      +bevts(ibevt+46,m)/totevt

            end do

               write(io,'(''  Total'',1p6e12.4)')
     &          sekc/totevt, sekn/totevt, sekf/totevt,
     &          sekd/totevt, sekg/totevt,
     &         (sekc+sekn+sekf+sekd+sekg)/totevt

         end if

*-----------------------------------------------------------------------
*     Medium-wise electron energy gain
*     for low energy per source
*-----------------------------------------------------------------------

               sekc = 0.0
               sekn = 0.0
               sekf = 0.0
               sekd = 0.0
               sekg = 0.0

            do m = 1, mxmat

               sekc = sekc + bevts(ibevt+39,m)
               sekn = sekn + bevts(ibevt+41,m)
               sekf = sekf + bevts(ibevt+43,m)
               sekd = sekd + bevts(ibevt+45,m)
               sekg = sekg + bevts(ibevt+47,m)

            end do

         if( sekc+sekn+sekf+sekd+sekg .gt. 0 ) then

               write(io,'(/79(''-''))')
               write(io,'('' Medium-wise electron energy gain'',
     &                    '' for low energy per source''
     &                    /79(''-'')/
     &                    '' Medium'',
     &                    ''   pair prod'',
     &                    ''     compton'',
     &                    ''   photo-ele'',
     &                    ''       auger'',
     &                    ''    knock-on'',
     &                    ''       Total'')')

            do m = 1, mxmat

               if( bevts(ibevt+39,m)+bevts(ibevt+41,m)
     &            +bevts(ibevt+43,m)+bevts(ibevt+45,m)
     &            +bevts(ibevt+47,m) .gt. 0.0 )
     &               write(io,'(i7,1p6e12.4)')
     &               idmn(m),bevts(ibevt+39,m)/totevt,
     &                       bevts(ibevt+41,m)/totevt,
     &                       bevts(ibevt+43,m)/totevt,
     &                       bevts(ibevt+45,m)/totevt,
     &                       bevts(ibevt+47,m)/totevt,
     &                       bevts(ibevt+39,m)/totevt
     &                      +bevts(ibevt+41,m)/totevt
     &                      +bevts(ibevt+43,m)/totevt
     &                      +bevts(ibevt+45,m)/totevt
     &                      +bevts(ibevt+47,m)/totevt

            end do

               write(io,'(''  Total'',1p6e12.4)')
     &          sekc/totevt, sekn/totevt, sekf/totevt,
     &          sekd/totevt, sekg/totevt,
     &         (sekc+sekn+sekf+sekd+sekg)/totevt

         end if

*-----------------------------------------------------------------------
*     Medium-wise electron energy loss
*     for low energy per source
*-----------------------------------------------------------------------

               sekc = 0.0
               sekn = 0.0

            do m = 1, mxmat

               sekc = sekc + bevts(ibevt+48,m)
               sekn = sekn + bevts(ibevt+49,m)

            end do

         if( sekc+sekn .gt. 0 ) then

               write(io,'(/79(''-''))')
               write(io,'('' Medium-wise electron energy loss'',
     &                    '' for low energy per source''
     &                    /79(''-'')/
     &                    '' Medium'',
     &                    ''     scatter'',
     &                    ''      bremss'',
     &                    ''       Total'')')

            do m = 1, mxmat

               if( bevts(ibevt+48,m)+
     &             bevts(ibevt+49,m) .gt. 0.0 )
     &               write(io,'(i7,1p6e12.4)')
     &               idmn(m),bevts(ibevt+48,m)/totevt,
     &                       bevts(ibevt+49,m)/totevt,
     &                       bevts(ibevt+48,m)/totevt
     &                      +bevts(ibevt+49,m)/totevt

            end do

               write(io,'(''  Total'',1p6e12.4)')
     &          sekc/totevt, sekn/totevt,
     &         (sekc+sekn)/totevt

         end if

*-----------------------------------------------------------------------
*      Medium-wise total number of n-coll
*      for high energy library per source
*-----------------------------------------------------------------------

               sekg = 0.0
               seke = 0.0
               sekn = 0.0

            do m = 1, mxmat

               sekg = sekg + bevts(ibevt+50,m)
               seke = seke + bevts(ibevt+51,m)
               sekn = sekn + bevts(ibevt+52,m)

            end do

         if( sekg+seke+sekn .gt. 0 ) then

               write(io,'(/79(''-''))')
               write(io,'('' Medium-wise total number of n-coll'',
     &                    '' for high energy library per source''
     &                    /79(''-'')/
     &                    '' Medium'',
     &                    ''      (N,xP)'',
     &                    ''     (N,xPI)'',
     &                    ''   (N,other)'',
     &                    ''  Total-Coll'')')

            do m = 1, mxmat

               if( bevts(ibevt+50,m)
     &            +bevts(ibevt+51,m)+bevts(ibevt+52,m) .gt. 0.0 )
     &               write(io,'(i7,1p6e12.4)')
     &               idmn(m),bevts(ibevt+50,m)/totevt,
     &                       bevts(ibevt+51,m)/totevt,
     &                       bevts(ibevt+52,m)/totevt,
     &                       bevts(ibevt+50,m)/totevt+
     &                       bevts(ibevt+51,m)/totevt+
     &                       bevts(ibevt+52,m)/totevt

            end do

               write(io,'(''  Total'',1p6e12.4)')
     &         sekg/totevt,
     &         seke/totevt,sekn/totevt,
     &         sekg/totevt+
     &         seke/totevt+sekn/totevt

         end if

*-----------------------------------------------------------------------
*      Medium-wise total number of p-coll
*      for high energy library per source
*-----------------------------------------------------------------------

               sekg = 0.0
               seke = 0.0
               sekn = 0.0

            do m = 1, mxmat

               sekg = sekg + bevts(ibevt+59,m)
               seke = seke + bevts(ibevt+60,m)
               sekn = sekn + bevts(ibevt+61,m)

            end do

         if( sekg+seke+sekn .gt. 0 ) then

               write(io,'(/79(''-''))')
               write(io,'('' Medium-wise total number of p-coll'',
     &                    '' for high energy library per source''
     &                    /79(''-'')/
     &                    '' Medium'',
     &                    ''     Capture'',
     &                    ''     Elastic'',
     &                    ''   Non-Elast'',
     &                    ''  Total-Coll'')')

            do m = 1, mxmat

               if( bevts(ibevt+59,m)+bevts(ibevt+60,m)
     &            +bevts(ibevt+61,m) .gt. 0.0 )
     &               write(io,'(i7,1p6e12.4)')
     &               idmn(m),bevts(ibevt+59,m)/totevt,
     &                       bevts(ibevt+60,m)/totevt,
     &                       bevts(ibevt+61,m)/totevt,
     &                       bevts(ibevt+59,m)/totevt+
     &                       bevts(ibevt+60,m)/totevt+
     &                       bevts(ibevt+61,m)/totevt

            end do

               write(io,'(''  Total'',1p6e12.4)')
     &         sekg/totevt,seke/totevt,sekn/totevt,
     &         sekg/totevt+seke/totevt+sekn/totevt

         end if

*-----------------------------------------------------------------------
*      Medium-wise total number of p-coll
*      for high energy library per source
*-----------------------------------------------------------------------

               sekg = 0.0
               seke = 0.0
               sekn = 0.0
               sekf = 0.0
               sekd = 0.0

            do m = 1, mxmat

               sekg = sekg + bevts(ibevt+53,m)
               seke = seke + bevts(ibevt+54,m)
               sekn = sekn + bevts(ibevt+55,m)
               sekf = sekf + bevts(ibevt+56,m)
               sekd = sekd + bevts(ibevt+58,m)

            end do

         if( sekg+seke+sekn+sekf+sekd .gt. 0 ) then

               write(io,'(/79(''-''))')
               write(io,'('' Medium-wise total number of p-coll'',
     &                    '' for high energy library per source''
     &                    /79(''-'')/
     &                    '' Medium'',
     &                    ''      (P,xP)'',
     &                    ''      (P,xN)'',
     &                    ''  (P,Photon)'',
     &                    ''     (P,xPI)'',
     &                    ''   (P,other)'',
     &                    ''  Total-Coll'')')

            do m = 1, mxmat

               if( bevts(ibevt+53,m)+bevts(ibevt+54,m)
     &            +bevts(ibevt+55,m)+bevts(ibevt+56,m)
     &            +bevts(ibevt+58,m) .gt. 0.0 )
     &               write(io,'(i7,1p6e12.4)')
     &               idmn(m),bevts(ibevt+53,m)/totevt,
     &                       bevts(ibevt+54,m)/totevt,
     &                       bevts(ibevt+58,m)/totevt,
     &                       bevts(ibevt+55,m)/totevt,
     &                       bevts(ibevt+56,m)/totevt,
     &                       bevts(ibevt+53,m)/totevt+
     &                       bevts(ibevt+54,m)/totevt+
     &                       bevts(ibevt+55,m)/totevt+
     &                       bevts(ibevt+56,m)/totevt+
     &                       bevts(ibevt+58,m)/totevt

            end do

               write(io,'(''  Total'',1p6e12.4)')
     &         sekg/totevt,
     &         seke/totevt,sekd/totevt,sekn/totevt,sekf/totevt,
     &         sekg/totevt+
     &         seke/totevt+sekd/totevt+sekn/totevt+sekf/totevt

         end if

cABE add @2014/08/12, to output muon reaction
*-----------------------------------------------------------------------
*     Medium-wise total number of collisions
*     for Muon Photonuclear Interaction
*-----------------------------------------------------------------------

               seky = 0.0
               sekb = 0.0
               sekp = 0.0
               sekc = 0.0

            do m = 1, mxmat
               seky = seky+ bevts(ibevt+64,m)
               sekb = sekb+ bevts(ibevt+66,m)
               sekp = sekp+ bevts(ibevt+67,m)
               sekc = seky+ bevts(ibevt+65,m)
            end do

         if( seky+sekb+sekp+sekc .gt. 0.d0 ) then

               write(io,'(/79(''-''))')
               write(io,'('' Medium-wise total number of collisions'',
     &                    '' for Muon Interaction per source''
     &                    /79(''-'')/
     &                    '' Medium'',
     &                    ''   photonucl'',
     &                    ''      bremss'',
     &                    ''   pair prod'',
     &                    ''     capture'',
     &                    ''       Total'')')

            do m = 1, mxmat

               if( bevts(ibevt+64,m)+bevts(ibevt+66,m)+bevts(ibevt+67,m)
     &            +bevts(ibevt+65,m) .gt. 0.d0 )
     &               write(io,'(i7,1p5e12.4)')
     &               idmn(m),bevts(ibevt+64,m) / totevt,
     &                       bevts(ibevt+66,m) / totevt,
     &                       bevts(ibevt+67,m) / totevt,
     &                       bevts(ibevt+65,m) / totevt,
     &                       bevts(ibevt+64,m) / totevt+
     &                       bevts(ibevt+66,m) / totevt+
     &                       bevts(ibevt+67,m) / totevt+
     &                       bevts(ibevt+65,m) / totevt

            end do

               write(io,'(''  Total'',1p5e12.4)')
     &         seky/totevt,sekb/totevt,sekp/totevt,sekc/totevt,
     &         seky/totevt+sekb/totevt+sekp/totevt+sekc/totevt

         end if
*-----------------------------------------------------------------------
*     Medium-wise total number of collisions
*     for Track structure mode
*-----------------------------------------------------------------------

               sekk = 0.0
               sekp = 0.0
               sekc = 0.0
               seka = 0.0

            do m = 1, mxmat
               sekk = sekk+ bevts(ibevt+68,m)
               sekp = sekp+ bevts(ibevt+69,m)
               sekc = sekc+ bevts(ibevt+70,m)
               seka = seka+ bevts(ibevt+71,m)
            end do

         if( sekk+sekp+sekc+seka .gt. 0.d0 ) then

               write(io,'(/79(''-''))')
               write(io,'('' Medium-wise total number of collisions'',
     &                    '' for Track structure mode per source''
     &                    /79(''-'')/
     &                    '' Medium'',
     &                    ''         ETS'',
     &                    ''    KURBUC p'',
     &                    ''    KURBUC C'',
     &                    ''      ITSART'',
     &                    ''       Total'')')

            do m = 1, mxmat

               if( bevts(ibevt+68,m)+bevts(ibevt+69,m)+bevts(ibevt+70,m)
     &            +bevts(ibevt+71,m) .gt. 0.d0 )
     &               write(io,'(i7,1p5e12.4)')
     &               idmn(m),bevts(ibevt+68,m) / totevt,
     &                       bevts(ibevt+69,m) / totevt,
     &                       bevts(ibevt+70,m) / totevt,
     &                       bevts(ibevt+71,m) / totevt,
     &                       bevts(ibevt+68,m) / totevt+
     &                       bevts(ibevt+69,m) / totevt+
     &                       bevts(ibevt+70,m) / totevt+
     &                       bevts(ibevt+71,m) / totevt

            end do

               write(io,'(''  Total'',1p5e12.4)')
     &         sekk/totevt,sekp/totevt,sekc/totevt,seka/totevt,
     &         sekk/totevt+sekp/totevt+sekc/totevt+seka/totevt

         end if


*-----------------------------------------------------------------------

      end if

*-----------------------------------------------------------------------

      end if

      if( infout .eq. 3 .or. infout .eq. 5 .or. infout .eq. 6 .or.
     &    infout .eq. 7 .or. infout .eq. 8 ) then

*-----------------------------------------------------------------------
*     particle data transfer to root
*-----------------------------------------------------------------------

         if( npe .gt. 1 ) then

            if( me .gt. 0 ) then

                     call parasr(aprodp,20,0)
                     call parasr(adcayp,20,0)
                     call parasr(astopp,20,0)
                     call parasr(atimep,20,0)
                     call parasr(aleakp,20,0)

                     call parasr(bprodp,20,0)
                     call parasr(bdcayp,20,0)
                     call parasr(bstopp,20,0)
                     call parasr(btimep,20,0)
                     call parasr(bleakp,20,0)

                     call parasi(nothp,4000,0)
                     call parasi(nojmp,4000,0)
                     call parasi(nostp,4000,0)
                     call parasi(notip,4000,0)
                     call parasi(nolkp,4000,0)

                     call parasi(nodcp,12000,0)

                     call parasi(nothn,1,0)
                     call parasi(nojmn,1,0)
                     call parasi(nostn,1,0)
                     call parasi(notin,1,0)
                     call parasi(nolkn,1,0)

                     call parasi(nodcn,1,0)

            else

                     do m = 1, 20

                        aprodp(m) = 0.0
                        adcayp(m) = 0.0
                        astopp(m) = 0.0
                        atimep(m) = 0.0
                        aleakp(m) = 0.0

                        bprodp(m) = 0.0
                        bdcayp(m) = 0.0
                        bstopp(m) = 0.0
                        btimep(m) = 0.0
                        bleakp(m) = 0.0

                     end do

                     do m = 1, 2000
                     do n = 1, 2

                        nothp(m,n) = 0
                        nojmp(m,n) = 0
                        nostp(m,n) = 0
                        notip(m,n) = 0
                        nolkp(m,n) = 0

                     end do
                     end do

                     do m = 1, 2000
                     do n = 1, 6

                        nodcp(m,n) = 0

                     end do
                     end do

                        nothn = 0
                        nojmn = 0
                        nostn = 0
                        notin = 0
                        nolkn = 0

                        nodcn = 0

*-----------------------------------------------------------------------

               do i = 1, npe - 1

               if( iccp(i) .eq. 0 ) then

                     call pararr(aprodr,20,i)
                     call pararr(adcayr,20,i)
                     call pararr(astopr,20,i)
                     call pararr(atimer,20,i)
                     call pararr(aleakr,20,i)

                     call pararr(bprodr,20,i)
                     call pararr(bdcayr,20,i)
                     call pararr(bstopr,20,i)
                     call pararr(btimer,20,i)
                     call pararr(bleakr,20,i)

                     call parari(nothr,4000,i)
                     call parari(nojmr,4000,i)
                     call parari(nostr,4000,i)
                     call parari(notir,4000,i)
                     call parari(nolkr,4000,i)

                     call parari(nodcr,12000,i)

                     call parari(nothu,1,i)
                     call parari(nojmu,1,i)
                     call parari(nostu,1,i)
                     call parari(notiu,1,i)
                     call parari(nolku,1,i)

                     call parari(nodcu,1,i)

*-----------------------------------------------------------------------

                  do m = 1, 20

                     aprodp(m) = aprodp(m) + aprodr(m)
                     adcayp(m) = adcayp(m) + adcayr(m)
                     astopp(m) = astopp(m) + astopr(m)
                     atimep(m) = atimep(m) + atimer(m)
                     aleakp(m) = aleakp(m) + aleakr(m)

                     bprodp(m) = bprodp(m) + bprodr(m)
                     bdcayp(m) = bdcayp(m) + bdcayr(m)
                     bstopp(m) = bstopp(m) + bstopr(m)
                     btimep(m) = btimep(m) + btimer(m)
                     bleakp(m) = bleakp(m) + bleakr(m)

                  end do

*-----------------------------------------------------------------------

                  if( nothu .gt. 0 ) then

                     if( nothn .eq. 0 ) then

                           nothn = nothu

                        do m = 1, nothn

                           nothp(m,1) = nothp(m,1) + nothr(m,1)
                           nothp(m,2) = nothp(m,2) + nothr(m,2)

                        end do

                     else

                        do 310 n = 1, nothu

                           do m = 1, nothn

                              if( nothr(n,2) .eq. nothp(m,2) ) then

                                 nothp(m,1) = nothp(m,1) + nothr(n,1)
                                 goto 310

                              end if

                           end do

                              nothn = nothn + 1
                              nothp(nothn,1) = nothr(n,1)
                              nothp(nothn,2) = nothr(n,2)

  310                   continue

                     end if

                  end if

*-----------------------------------------------------------------------

                  if( notiu .gt. 0 ) then

                     if( notin .eq. 0 ) then

                           notin = notiu

                        do m = 1, notin

                           notip(m,1) = notip(m,1) + notir(m,1)
                           notip(m,2) = notip(m,2) + notir(m,2)

                        end do

                     else

                        do 290 n = 1, notiu

                           do m = 1, notin

                              if( notir(n,2) .eq. notip(m,2) ) then

                                 notip(m,1) = notip(m,1) + notir(n,1)
                                 goto 290

                              end if

                           end do

                              notin = notin + 1
                              notip(nostn,1) = notir(n,1)
                              notip(nostn,2) = notir(n,2)

  290                   continue

                     end if

                  end if

*-----------------------------------------------------------------------

                  if( nostu .gt. 0 ) then

                     if( nostn .eq. 0 ) then

                           nostn = nostu

                        do m = 1, nostn

                           nostp(m,1) = nostp(m,1) + nostr(m,1)
                           nostp(m,2) = nostp(m,2) + nostr(m,2)

                        end do

                     else

                        do 320 n = 1, nostu

                           do m = 1, nostn

                              if( nostr(n,2) .eq. nostp(m,2) ) then

                                 nostp(m,1) = nostp(m,1) + nostr(n,1)
                                 goto 320

                              end if

                           end do

                              nostn = nostn + 1
                              nostp(nostn,1) = nostr(n,1)
                              nostp(nostn,2) = nostr(n,2)

  320                   continue

                     end if

                  end if

*-----------------------------------------------------------------------

                  if( nojmu .gt. 0 ) then

                     if( nojmn .eq. 0 ) then

                           nojmn = nojmu

                        do m = 1, nojmn

                           nojmp(m,1) = nojmp(m,1) + nojmr(m,1)
                           nojmp(m,2) = nojmp(m,2) + nojmr(m,2)

                        end do

                     else

                        do 330 n = 1, nojmu

                           do m = 1, nojmn

                              if( nojmr(n,2) .eq. nojmp(m,2) ) then

                                 nojmp(m,1) = nojmp(m,1) + nojmr(n,1)
                                 goto 330

                              end if

                           end do

                              nojmn = nojmn + 1
                              nojmp(nojmn,1) = nojmr(n,1)
                              nojmp(nojmn,2) = nojmr(n,2)

  330                   continue

                     end if

                  end if

*-----------------------------------------------------------------------

                  if( nolku .gt. 0 ) then

                     if( nolkn .eq. 0 ) then

                           nolkn = nolku

                        do m = 1, nolkn

                           nolkp(m,1) = nolkp(m,1) + nolkr(m,1)
                           nolkp(m,2) = nolkp(m,2) + nolkr(m,2)

                        end do

                     else

                        do 340 n = 1, nolku

                           do m = 1, nolkn

                              if( nolkr(n,2) .eq. nolkp(m,2) ) then

                                 nolkp(m,1) = nolkp(m,1) + nolkr(n,1)
                                 goto 340

                              end if

                           end do

                              nolkn = nolkn + 1
                              nolkp(nolkn,1) = nolkr(n,1)
                              nolkp(nolkn,2) = nolkr(n,2)

  340                   continue

                     end if

                  end if

*-----------------------------------------------------------------------

                  if( nodcu .gt. 0 ) then

                     if( nodcn .eq. 0 ) then

                           nodcn = nodcu

                        do m = 1, nodcn

                           do k = 1, 6

                              nodcp(m,k) = nodcp(m,k) + nodcr(m,k)

                           end do

                        end do

                     else

                        do 350 n = 1, nodcu

                           do m = 1, nodcn

                              if( nodcr(n,2) .eq. 2 ) then

                                 if( nodcr(n,2) .eq. nodcp(m,2) .and.
     &                               nodcr(n,3) .eq. nodcp(m,3) .and.
     &                               nodcr(n,4) .eq. nodcp(m,4) .and.
     &                               nodcr(n,5) .eq. nodcp(m,5) ) then

                                    nodcp(m,1) = nodcp(m,1) + nodcr(n,1)
                                    goto 350

                                 end if

                              else if( nodcr(n,2) .eq. 3 ) then

                                 if( nodcr(n,2) .eq. nodcp(m,2) .and.
     &                               nodcr(n,3) .eq. nodcp(m,3) .and.
     &                               nodcr(n,4) .eq. nodcp(m,4) .and.
     &                               nodcr(n,5) .eq. nodcp(m,5) .and.
     &                               nodcr(n,6) .eq. nodcp(m,6) ) then

                                    nodcp(m,1) = nodcp(m,1) + nodcr(n,1)
                                    goto 350

                                 end if

                              end if

                           end do

                              nodcn = nodcn + 1
                              nodcp(nodcn,1) = nodcr(n,1)
                              nodcp(nodcn,2) = nodcr(n,2)

                           if( nodcr(n,2) .eq. 2 ) then

                              nodcp(nodcn,3) = nodcr(n,3)
                              nodcp(nodcn,4) = nodcr(n,4)
                              nodcp(nodcn,5) = nodcr(n,5)

                           else if( nodcr(n,2) .eq. 3 ) then

                              nodcp(nodcn,3) = nodcr(n,3)
                              nodcp(nodcn,4) = nodcr(n,4)
                              nodcp(nodcn,5) = nodcr(n,5)
                              nodcp(nodcn,6) = nodcr(n,6)

                           end if

  350                   continue

                     end if

                  end if

*-----------------------------------------------------------------------

               end if
               end do

            end if

         end if

*-----------------------------------------------------------------------
*     List of transport particles
*-----------------------------------------------------------------------

      if( me .eq. 0 ) then

            write(io,'(/79(''-''))')
            write(io,'('' List of transport particles'',
     &                 '' (not including source)'',
     &                 /79(''-'')/2x,
     &                 ''   Name      kf-code     '',
     &                 ''   mass       charge   baryon'')')

            do i = 1, 19

               if( i .ne. 11 ) then

               if( idpat(i) .gt. 0 .and. bprodp(i) .gt. 0.0 ) then

                  kf  = kfft(i)

               if( kf .ne. 0 .and. i .lt. 15 ) then

                  rmp = rmtyp(11,kf)
                  kch = ichgf(11,kf)
                  kby = ibryf(11,kf)

                  write(io,'(5x,a8,i9,3x,f10.1,4x,i6,2x,i6)')
     &                 pname(i), kf, rmp, kch, kby

               else if( i .ge. 15 .and. i .le. 18 ) then

                  rmp = rmtyp(i,kf)
                  kch = ichgf(i,kf)
                  kby = ibryf(i,kf)

                  write(io,'(5x,a8,i9,3x,f10.1,4x,i6,2x,i6)')
     &                 pname(i), kf, rmp, kch, kby

               else if( i .eq. 19 ) then

                  rmp = 0.0
                  kch = 0
                  kby = 0

                  write(io,'(5x,a8,i9,3x,f10.1,4x,i6,2x,i6)')
     &                 pname(i), kf, rmp, kch, kby

               end if

               end if
               end if

            end do

         if( idpat(11) .eq. 1 .and. idono .gt. 0 .and.
     &       nothn .gt. 0 ) then

                  write(io,'()')

            do 410 i = 1, idono

                  kf  = idoth(i)

                     do m = 1, nothn

                        if( kf .eq. nothp(m,2) ) goto 400

                     end do

                     goto 410

  400                continue

                  rmp = rmtyp(11,kf)
                  kch = ichgf(11,kf)
                  kby = ibryf(11,kf)

                  call jamname(kf,0,0,qname)

                  write(io,'(5x,a7,i10,3x,f10.1,4x,i6,2x,i6)')
     &                 qname(1:7), kf, rmp, kch, kby

  410       continue

         end if

*-----------------------------------------------------------------------
*        summary of produced particles and fission
*-----------------------------------------------------------------------

         do it = 1, 20

            if( aprodp(it) .gt. 0.0 ) goto 48

         end do

            goto 49

   48       write(io,'(/79(''-''))')
            write(io,'( '' prod. particles'',
     &                  ''       number        weight'',
     &                  ''        weight per source'')')
            write(io,'( 79(''-''))')

            do i = 1, 20
               noths(i) = i
            end do

            do i = 1, 19
                     bnmax = bprodp(noths(i))
               do j = i + 1, 20
                  if( bprodp(noths(j)) .gt. bnmax ) then
                     bnmax = bprodp(noths(j))
                     ii = noths(i)
                     noths(i) = noths(j)
                     noths(j) = ii
                  end if
               end do
            end do

         do is = 1, 20

            it = noths(is)

            if( aprodp(it) .gt. 0.0 ) then

               if( it .ne. 20 ) then

                  write(io,'(5x,a8,f16.0,4x,1pe14.7,4x,1pe14.7)')
     &                  pname(it), bprodp(it),
     &                  aprodp(it), aprodp(it) / totevt

               else

                  write(io,'(5x,''fission '',f16.0,4x,1pe14.7,4x,
     &                  1pe14.7)')
     &                  bprodp(it), aprodp(it), aprodp(it) / totevt

               end if

            end if

         end do

   49    continue

*-----------------------------------------------------------------------
*        summary of the other kind of produced particles
*-----------------------------------------------------------------------

         if( nothn .gt. 0 ) then

            do i = 1, nothn
               noths(i) = i
            end do

            do i = 1, nothn - 1
                     bnmax = nothp(noths(i),1)
               do j = i + 1, nothn
                  if( nothp(noths(j),1) .gt. bnmax ) then
                     bnmax = nothp(noths(j),1)
                     ii = noths(i)
                     noths(i) = noths(j)
                     noths(j) = ii
                  end if
               end do
            end do

            write(io,'(/79(''-''))')
            write(io,'('' the other kind of produced particles''
     &                 /79(''-'')/2x,
     &                 ''   Name      kf-code        number  '',
     &                 ''   mass       charge   baryon'')')

            do i = 1, nothn

               kf  = nothp(noths(i),2)
               rmp = rmtyp(11,kf)
               kch = ichgf(11,kf)
               kby = ibryf(11,kf)

               call jamname(kf,0,0,qname)

               write(io,'(5x,a7,i10,i14,f10.1,4x,i6,2x,i6)')
     &                        qname(1:7), kf, nothp(noths(i),1),
     &                        rmp, kch, kby

            end do

         end if

*-----------------------------------------------------------------------
*        summary of the other kind of particles into JAM
*-----------------------------------------------------------------------

         if( nojmn .gt. 0 ) then

            do i = 1, nojmn
               noths(i) = i
            end do

            do i = 1, nojmn - 1
                     bnmax = nojmp(noths(i),1)
               do j = i + 1, nojmn
                  if( nojmp(noths(j),1) .gt. bnmax ) then
                     bnmax = nojmp(noths(j),1)
                     ii = noths(i)
                     noths(i) = noths(j)
                     noths(j) = ii
                  end if
               end do
            end do

            write(io,'(/79(''-''))')
            write(io,'('' the other kind of particles into JAM''
     &                 /79(''-'')/2x,
     &                 ''   Name      kf-code        number  '',
     &                 ''   mass       charge   baryon'')')

            do i = 1, nojmn

               kf  = nojmp(noths(i),2)
               rmp = rmtyp(11,kf)
               kch = ichgf(11,kf)
               kby = ibryf(11,kf)

               call jamname(kf,0,0,qname)

               write(io,'(5x,a7,i10,i14,f10.1,4x,i6,2x,i6)')
     &                        qname(1:7), kf, nojmp(noths(i),1),
     &                        rmp, kch, kby

            end do

         end if

*-----------------------------------------------------------------------
*        summary of particle decay
*-----------------------------------------------------------------------

         do it = 1, 20

            if( adcayp(it) .gt. 0.0 ) goto 68

         end do

            goto 69

   68       write(io,'(/79(''-''))')
            write(io,'( '' particle decays'',
     &                  ''       number        weight'',
     &                  ''        weight per source'')')
            write(io,'( 79(''-''))')

            do i = 1, 20
               noths(i) = i
            end do

            do i = 1, 19
                     bnmax = bdcayp(noths(i))
               do j = i + 1, 20
                  if( bdcayp(noths(j)) .gt. bnmax ) then
                     bnmax = bdcayp(noths(j))
                     ii = noths(i)
                     noths(i) = noths(j)
                     noths(j) = ii
                  end if
               end do
            end do

         do is = 1, 20

            it = noths(is)

            if( adcayp(it) .gt. 0.0 ) then

               write(io,'(5x,a8,f16.0,4x,1pe14.7,4x,1pe14.7)')
     &               pname(it), bdcayp(it),
     &               adcayp(it), adcayp(it) / totevt

            end if

         end do

   69    continue

*-----------------------------------------------------------------------
*        summary of the other kind of decay particles
*-----------------------------------------------------------------------

         if( nodcn .gt. 0 ) then

            do i = 1, nodcn
               noths(i) = i
            end do

            do i = 1, nodcn - 1
                     bnmax = nodcp(noths(i),1)
               do j = i + 1, nodcn
                  if( nodcp(noths(j),1) .gt. bnmax ) then
                     bnmax = nodcp(noths(j),1)
                     ii = noths(i)
                     noths(i) = noths(j)
                     noths(j) = ii
                  end if
               end do
            end do

            write(io,'(/79(''-''))')
            write(io,'('' the other kind of decay particles''
     &                 /79(''-''))')

            do i = 1, nodcn

                  kf0 = nodcp(noths(i),2)
                  call jamname(kf0,0,0,qnamd(0))

               if( nodcp(noths(i),3) .eq. 2 ) then

                  kf1 = nodcp(noths(i),4)
                  kf2 = nodcp(noths(i),5)
                  call jamname(kf1,0,0,qnamd(1))
                  call jamname(kf2,0,0,qnamd(2))

                  write(io,'(5x,a7,'' ->  '',a7,
     &                       '' +  '',a7,11x,i14)')
     &            qnamd(0),qnamd(1),qnamd(2),nodcp(noths(i),1)

               else if( nodcp(noths(i),3) .eq. 3 ) then

                  kf1 = nodcp(noths(i),4)
                  kf2 = nodcp(noths(i),5)
                  kf3 = nodcp(noths(i),6)
                  call jamname(kf1,0,0,qnamd(1))
                  call jamname(kf2,0,0,qnamd(2))
                  call jamname(kf3,0,0,qnamd(3))

                  write(io,'(5x,a7,'' ->  '',a7,
     &                       '' +  '',a7,
     &                       '' +  '',a7,i14)')
     &            qnamd(0),qnamd(1),qnamd(2),qnamd(3),nodcp(noths(i),1)

               end if

            end do

         end if

*-----------------------------------------------------------------------
*        summary of time stopped particles
*-----------------------------------------------------------------------

         do it = 1, 20

            if( atimep(it) .gt. 0.0 ) goto 18

         end do

            goto 19

   18       write(io,'(/79(''-''))')
            write(io,'( '' time stop. part.'',
     &                  ''      number        weight'',
     &                  ''        weight per source'')')
            write(io,'( 79(''-''))')

            do i = 1, 20
               noths(i) = i
            end do

            do i = 1, 19
                     bnmax = btimep(noths(i))
               do j = i + 1, 20
                  if( btimep(noths(j)) .gt. bnmax ) then
                     bnmax = btimep(noths(j))
                     ii = noths(i)
                     noths(i) = noths(j)
                     noths(j) = ii
                  end if
               end do
            end do

         do is = 1, 20

            it = noths(is)

            if( atimep(it) .gt. 0.0 ) then

               write(io,'(5x,a8,f16.0,4x,1pe14.7,4x,1pe14.7)')
     &               pname(it), btimep(it),
     &               atimep(it), atimep(it) / totevt

            end if

         end do

   19    continue

*-----------------------------------------------------------------------
*        summary of the other kind of time stopped particles
*-----------------------------------------------------------------------

         if( notin .gt. 0 ) then

            do i = 1, notin
               noths(i) = i
            end do

            do i = 1, notin - 1
                     bnmax = notip(noths(i),1)
               do j = i + 1, notin
                  if( notip(noths(j),1) .gt. bnmax ) then
                     bnmax = notip(noths(j),1)
                     ii = noths(i)
                     noths(i) = noths(j)
                     noths(j) = ii
                  end if
               end do
            end do

            write(io,'(/79(''-''))')
            write(io,'('' the other kind of time stopped particles''
     &                 /79(''-'')/2x,
     &                 ''   Name      kf-code        number  '',
     &                 ''   mass       charge   baryon'')')

            do i = 1, notin

               kf  = notip(noths(i),2)

               rmp = rmtyp(11,kf)
               kch = ichgf(11,kf)
               kby = ibryf(11,kf)

               call jamname(kf,0,0,qname)

               write(io,'(5x,a7,i10,i14,f10.1,4x,i6,2x,i6)')
     &                        qname(1:7), kf, notip(noths(i),1),
     &                        rmp, kch, kby

            end do

         end if

*-----------------------------------------------------------------------
*        summary of energy stopped particles
*-----------------------------------------------------------------------

         do it = 1, 20

            if( astopp(it) .gt. 0.0 ) goto 38

         end do

            goto 39

   38       write(io,'(/79(''-''))')
            write(io,'( '' stop. particles.'',
     &                  ''    number        weight'',
     &                  ''        weight per source'')')
            write(io,'( 79(''-''))')

            do i = 1, 20
               noths(i) = i
            end do

            do i = 1, 19
                     bnmax = bstopp(noths(i))
               do j = i + 1, 20
                  if( bstopp(noths(j)) .gt. bnmax ) then
                     bnmax = bstopp(noths(j))
                     ii = noths(i)
                     noths(i) = noths(j)
                     noths(j) = ii
                  end if
               end do
            end do

         do is = 1, 20

            it = noths(is)

            if( astopp(it) .gt. 0.0 ) then

               write(io,'(5x,a8,f16.0,4x,1pe14.7,4x,1pe14.7)')
     &               pname(it), bstopp(it),
     &               astopp(it), astopp(it) / totevt

            end if

         end do

   39    continue

*-----------------------------------------------------------------------
*        summary of the other kind of energy stopped particles
*-----------------------------------------------------------------------

         if( nostn .gt. 0 ) then

            do i = 1, nostn
               noths(i) = i
            end do

            do i = 1, nostn - 1
                     bnmax = nostp(noths(i),1)
               do j = i + 1, nostn
                  if( nostp(noths(j),1) .gt. bnmax ) then
                     bnmax = nostp(noths(j),1)
                     ii = noths(i)
                     noths(i) = noths(j)
                     noths(j) = ii
                  end if
               end do
            end do

            write(io,'(/79(''-''))')
            write(io,'('' the other kind of stopped particles''
     &                 /79(''-'')/2x,
     &                 ''   Name      kf-code        number  '',
     &                 ''   mass       charge   baryon'')')

            do i = 1, nostn

               kf  = nostp(noths(i),2)

               rmp = rmtyp(11,kf)
               kch = ichgf(11,kf)
               kby = ibryf(11,kf)

               call jamname(kf,0,0,qname)

               write(io,'(5x,a7,i10,i14,f10.1,4x,i6,2x,i6)')
     &                        qname(1:7), kf, nostp(noths(i),1),
     &                        rmp, kch, kby

            end do

         end if

*-----------------------------------------------------------------------
*        summary of leakage particles
*-----------------------------------------------------------------------

         do it = 1, 20

            if( aleakp(it) .gt. 0.0 ) goto 28

         end do

            goto 29

   28       write(io,'(/79(''-''))')
            write(io,'( '' leak. particles'',
     &                  ''       number        weight'',
     &                  ''        weight per source'')')
            write(io,'( 79(''-''))')

            do i = 1, 20
               noths(i) = i
            end do

            do i = 1, 19
                     bnmax = bleakp(noths(i))
               do j = i + 1, 20
                  if( bleakp(noths(j)) .gt. bnmax ) then
                     bnmax = bleakp(noths(j))
                     ii = noths(i)
                     noths(i) = noths(j)
                     noths(j) = ii
                  end if
               end do
            end do

         do is = 1, 20

            it = noths(is)

            if( aleakp(it) .gt. 0.0 ) then

               write(io,'(5x,a8,f16.0,4x,1pe14.7,4x,1pe14.7)')
     &               pname(it), bleakp(it),
     &               aleakp(it), aleakp(it) / totevt

            end if

         end do

   29    continue

*-----------------------------------------------------------------------
*        summary of the other kind of leakage particles
*-----------------------------------------------------------------------

         if( nolkn .gt. 0 ) then

            do i = 1, nolkn
               noths(i) = i
            end do

            do i = 1, nolkn - 1
                     bnmax = nolkp(noths(i),1)
               do j = i + 1, nolkn
                  if( nolkp(noths(j),1) .gt. bnmax ) then
                     bnmax = nolkp(noths(j),1)
                     ii = noths(i)
                     noths(i) = noths(j)
                     noths(j) = ii
                  end if
               end do
            end do

            write(io,'(/79(''-''))')
            write(io,'('' the other kind of leakage particles''
     &                 /79(''-'')/2x,
     &                 ''   Name      kf-code        number  '',
     &                 ''   mass       charge   baryon'')')

            do i = 1, nolkn

               kf  = nolkp(noths(i),2)
               rmp = rmtyp(11,kf)
               kch = ichgf(11,kf)
               kby = ibryf(11,kf)

               call jamname(kf,0,0,qname)

               write(io,'(5x,a7,i10,i14,f10.1,4x,i6,2x,i6)')
     &                        qname(1:7), kf, nolkp(noths(i),1),
     &                        rmp, kch, kby

            end do

         end if

      end if

*-----------------------------------------------------------------------
*        cut-off neutron
*-----------------------------------------------------------------------

         if( npe .gt. 1 ) then

            if( me .gt. 0 ) then

                     call parasr(rtneut,1,0)
                     call parasr(wtneut,1,0)

            else

                     rtneut = 0.0
                     wtneut = 0.0

               do i = 1, npe - 1

                  if( iccp(i) .eq. 0 ) then

                     call pararr(rtneur,1,i)
                     call pararr(wtneur,1,i)

                     rtneut = rtneut + rtneur
                     wtneut = wtneut + wtneur

                  end if

               end do

            end if

         end if

*-----------------------------------------------------------------------

      if( me .eq. 0 .and. incut .ne. 0 ) then

            write(io,'(/79(''-''))')
            write(io,'( '' cut-off neutron'',
     &                  ''       number        weight'',
     &                  ''        number per source'')')
            write(io,'( 79(''-''))')

               write(io,'(5x,a8,f16.0,4x,1pe14.7,4x,1pe14.7)')
     &               pname(2), rtneut,
     &               wtneut, rtneut / totevt

      end if

*-----------------------------------------------------------------------
*        cut-off photon
*-----------------------------------------------------------------------

         if( npe .gt. 1 ) then

            if( me .gt. 0 ) then

                     call parasr(rtgamm,1,0)
                     call parasr(wtgamm,1,0)

            else

                     rtgamm = 0.0
                     wtgamm = 0.0

               do i = 1, npe - 1

                  if( iccp(i) .eq. 0 ) then

                     call pararr(rtgamr,1,i)
                     call pararr(wtgamr,1,i)

                     rtgamm = rtgamm + rtgamr
                     wtgamm = wtgamm + wtgamr

                  end if

               end do

            end if

         end if

*-----------------------------------------------------------------------

      if( me .eq. 0 .and. igcut .ne. 0 ) then

            write(io,'(/79(''-''))')
            write(io,'( '' cut-off photon '',
     &                  ''       number        weight'',
     &                  ''        number per source'')')
            write(io,'( 79(''-''))')

               write(io,'(5x,a8,f16.0,4x,1pe14.7,4x,1pe14.7)')
     &               pname(14), rtgamm,
     &               wtgamm, rtgamm / totevt

      end if

*-----------------------------------------------------------------------
*        cut-off electron
*-----------------------------------------------------------------------

         if( npe .gt. 1 ) then

            if( me .gt. 0 ) then

                     call parasr(rtelen,1,0)
                     call parasr(wtelen,1,0)

            else

                     rtelen = 0.0
                     wtelen = 0.0

               do i = 1, npe - 1

                  if( iccp(i) .eq. 0 ) then

                     call pararr(rteler,1,i)
                     call pararr(wteler,1,i)

                     rtelen = rtelen + rteler
                     wtelen = wtelen + wteler

                  end if

               end do

            end if

         end if

*-----------------------------------------------------------------------

      if( me .eq. 0 .and. igcut .eq. 3 ) then

            write(io,'(/79(''-''))')
            write(io,'( '' cut-off electron '',
     &                  ''       number        weight'',
     &                  ''        number per source'')')
            write(io,'( 79(''-''))')

               write(io,'(5x,a8,f16.0,4x,1pe14.7,4x,1pe14.7)')
     &               pname(12), rtelen,
     &               wtelen, rtelen / totevt

      end if

*-----------------------------------------------------------------------
*        cut-off positron
*-----------------------------------------------------------------------

         if( npe .gt. 1 ) then

            if( me .gt. 0 ) then

                     call parasr(rtelep,1,0)
                     call parasr(wtelep,1,0)

            else

                     rtelep = 0.0
                     wtelep = 0.0

               do i = 1, npe - 1

                  if( iccp(i) .eq. 0 ) then

                     call pararr(rteler,1,i)
                     call pararr(wteler,1,i)

                     rtelep = rtelep + rteler
                     wtelep = wtelep + wteler

                  end if

               end do

            end if

         end if

*-----------------------------------------------------------------------

      if( me .eq. 0 .and. igcut .eq. 3 ) then

            write(io,'(/79(''-''))')
            write(io,'( '' cut-off positron '',
     &                  ''       number        weight'',
     &                  ''        number per source'')')
            write(io,'( 79(''-''))')

               write(io,'(5x,a8,f16.0,4x,1pe14.7,4x,1pe14.7)')
     &               pname(12), rtelep,
     &               wtelep, rtelep / totevt

      end if

*-----------------------------------------------------------------------
*        cut-off proton
*-----------------------------------------------------------------------

         if( npe .gt. 1 ) then

            if( me .gt. 0 ) then

                     call parasr(rtprot,1,0)
                     call parasr(wtprot,1,0)

            else

                     rtprot = 0.0
                     wtprot = 0.0

               do i = 1, npe - 1

                  if( iccp(i) .eq. 0 ) then

                     call pararr(rtpror,1,i)
                     call pararr(wtpror,1,i)

                     rtprot = rtprot + rtpror
                     wtprot = wtprot + wtpror

                  end if

               end do

            end if

         end if

*-----------------------------------------------------------------------

      if( me .eq. 0 .and. ipcut .ne. 0 ) then

            write(io,'(/79(''-''))')
            write(io,'( '' cut-off proton '',
     &                  ''       number        weight'',
     &                  ''        number per source'')')
            write(io,'( 79(''-''))')

               write(io,'(5x,a8,f16.0,4x,1pe14.7,4x,1pe14.7)')
     &               pname(1), rtprot,
     &               wtprot, rtprot / totevt

      end if

*-----------------------------------------------------------------------

      end if

*-----------------------------------------------------------------------
*        extra bank access
*-----------------------------------------------------------------------

         if( npe .gt. 1 ) then

            if( me .gt. 0 ) then

                     call parasi(itbnk,1,0)
                     call parasi(jtbnk,1,0)

            else

                     itbnk = 0
                     jtbnk = 0

               do i = 1, npe - 1

                  if( iccp(i) .eq. 0 ) then

                     call parari(itbnkr,1,i)
                     call parari(jtbnkr,1,i)

                     if( itbnkr .gt. itbnk ) itbnk = itbnkr
                     if( jtbnkr .gt. jtbnk ) jtbnk = jtbnkr

                  end if

               end do

            end if

         end if

*-----------------------------------------------------------------------

      if( me .eq. 0 .and. ( jtbnk .gt. 0 .or. itbnk .gt. 0 ) ) then

            write(io,'(/79(''-''))')
            write(io,'( '' extra bank access'',
     &                  ''     memory        temp file'',
     &                  ''      max. number per source'')')
            write(io,'( 79(''-''))')

               write(io,'(13x,f16.0,f16.0)')
     &               dble(jtbnk), dble(itbnk)

      end if

*-----------------------------------------------------------------------

      if( me .eq. 0 ) then

            write(io,'(/79(''-''))')
            if (irestart .ne. 0 .and. mstz(3) .ne. maxcasres
     &      .and. istdev .eq. 1 ) then
              do i = 1, 2
                write(ios(i),'('' maxcas is adjusted from '',i8,
     &                   '' to '',i8,
     &                   '', which is written in restart file'')')
     &                 mstz(3), maxcasres
              end do
            end if

            if (idmpmode.eq.1)then
             do i=1,2
              write(ios(i),'('' idmpmode = 1'')')
              write(ios(i),'('' total source is adopted'',
     &              '' from dump source and totfact is ignored'')')
              if(mstz(3).ne.maxcasdmp)then
               write(ios(i),'('' maxcas is adjusted from '',i8,
     &                   '' to '',i8,
     &                   '', which is written in dump source'')')
     &                 mstz(3), maxcasdmp
              endif
              if(mstz(4).ne.maxbchdmp) then
               write(ios(i),'('' maxbch is adjusted from '',i8,
     &                   '' to '',i8,
     &                   '', which is written in dump source'')')
     &                 mstz(4), maxbchdmp
              endif
             enddo
             write(io,'( 79(''-''))')
            end if
c------------------
            write(io,'( '' source: maxcas'',
     &                  ''    maxbch'',
     &                  ''      irskip'',
     &                  ''   average weight          total source'')')
            write(io,'( 79(''-''))')

               write(io,'(3i12,3x,1pe14.7,3x,1pe20.13)')
     &               maxcas, maxbch,
     &               irskip,
     &               rsouin / rcasc, rsouin

         if( ispfn .ne. 0 ) then

            write(io,'(/79(''-''))')
            write(io,'( '' s-fiss source'',
     &                  ''      total source'',
     &                  ''        zero neutron'',
     &                  ''       total neutron'')')
            write(io,'( 79(''-''))')

               write(io,'(12x,1pe20.13,1pe20.13,1pe20.13)')
     &               rsouin, rspfz, rspfn

         end if

cFURUTA20160128--------------------------------------------------------
         if(ndbatima(1)+ndbatima(2)+ndbatima(3).gt.0)then
          write(io,'(/79(''-''))')
          write(io,'( '' used ATIMA database'')')
          write(io,'( 79(''-''))')
          if(ndbatima(1).gt.0)
     &         write(io,'( ''     for dedx  ='',i9,
     &                     '' / mdbatima ='',i9)') ndbatima(1), mdbatima
          if(ndbatima(2).gt.0)
     &         write(io,'( ''     for edisp ='',i9,
     &                     '' / mdbatima ='',i9)') ndbatima(2), mdbatima
          if(ndbatima(3).gt.0)
     &         write(io,'( ''     for nsprd ='',i9,
     &                     '' / mdbatima ='',i9)') ndbatima(3), mdbatima
          if(ndbatima(1).eq.mdbatima
     &         .or.ndbatima(2).eq.mdbatima
     &         .or.ndbatima(3).eq.mdbatima)then
           write(io,'( '' **** WARNING: maximum number of'',
     &          '' database reached *****************'')')
           write(io,'( ''  calculation may be speeded up'',
     &          '' by increasing MDBATIMA or DBCUTOFF'')')
           write(io,'( '' '',66(''*''))')
          endif
         endif
c----------------------------------------------------------------------

      end if

*-----------------------------------------------------------------------
*        summary CG error
*-----------------------------------------------------------------------

         if( npe .gt. 1 ) then

            if( me .gt. 0 ) then

                     call parasi(ilost,1,0)
                     call parasi(nlost,1,0)
                     call parasi(novp,1,0)
                     call parasi(icger,1,0)
                     call parasi(ncger,1,0)
                     call parasi(nrovp,3000,0)

            else

                     ilost = 0
                     nlost = 0
                     novp  = 0
                     icger  = 0
                     ncger  = 0

                  do m = 1, 3
                  do n = 1, 1000

                     nrovp(m,n) = 0.0

                  end do
                  end do

               do i = 1, npe - 1
               if( iccp(i) .eq. 0 ) then

                     call parari(ilosr,1,i)
                     call parari(nlosr,1,i)
                     call parari(novr,1,i)
                     call parari(icgeu,1,i)
                     call parari(ncgeu,1,i)
                     call parari(nrovr,3000,i)

                     ilost = ilost + ilosr
                     nlost = nlost + nlosr
                     ncger = ncger + ncgeu

                  if( icgeu .gt. 0 ) then

                     if( icger .eq. 0 ) then

                        icger = icgeu
                        novp  = novr

                        do m = 1, 3
                        do n = 1, novp

                           nrovp(m,n) = nrovr(m,n)

                        end do
                        end do

                     else

                        icger = icger + icgeu

                        do 200 j = 1, novr

                           do k = 1, novp

                              if( nrovp(1,k) .eq. nrovr(1,j) .and.
     &                            nrovp(2,k) .eq. nrovr(2,j) ) then

                                 nrovp(3,k) = nrovp(3,k) + nrovr(3,j)

                                 goto 200

                              end if

                           end do

                           novp = novp + 1

                           if( novp .le. 1000 ) then

                              nrovp(1,novp) = nrovr(1,j)
                              nrovp(2,novp) = nrovr(2,j)
                              nrovp(3,novp) = nrovr(3,j)

                           end if

  200                   continue

                     end if

                  end if

               end if
               end do

            end if

         end if

*-----------------------------------------------------------------------

      if( me .eq. 0 ) then

            write(io,'(/79(''-''))')
            write(io,'('' Geometry error summary'')')
            write(io,'( 79(''-''))')

               write(io,'( '' Number of lost particles     ='',i6,
     &                     '' / nlost ='',i9)') ilost, nlost

               write(io,'( '' Number of geometry recovering ='',i6)')
     &           icger

               write(io,'( '' Number of unrecovered errors ='',i6)')
     &           ncger

            if( icger .gt. 0 ) then

               do i = 1, novp
                  noths(i) = i
               end do

               do i = 1, novp - 1
                        bnmax = nrovp(3,noths(i))
                  do j = i + 1, novp
                     if( nrovp(3,noths(j)) .gt. bnmax ) then
                        bnmax = nrovp(3,noths(j))
                        ii = noths(i)
                        noths(i) = noths(j)
                        noths(j) = ii
                     end if
                  end do
               end do

               write(io,'(/'' List of errors/warnings of overlap '',
     &                     ''regions'')')
               write(io,'( ''   region1   region2    number'')')

               do i = 1, novp

               write(io,'(3i10)')
     &         nrovp(2,noths(i)), nrovp(1,noths(i)), nrovp(3,noths(i))

               end do

            end if

      end if


*-----------------------------------------------------------------------

      if( me .eq. 0 ) then

            write(io,'(/79(''-''))')
            if (irestart .ne. 0 ) then
              do i = 1, 1
                write(ios(i),
     &           '('' This is restart calculation.'',
     &           '' Initial random seed is read from '', a )')
     &           crdrfln(1:irdrfll)
              end do
            end if
            if (irestart .ne. 0 .and. .not. lrijkeqrf ) then
              do i = 1, 1
                write(ios(i),
     &                  '('' Warning: This random seed is different'',
     &                    '' from those written in other resfiles.'')')
              end do
            end if

           if ( nrandgen .eq. 0 ) then ! when LCG (2021.4.20)
            write(io,'(1x,
     &         "initial random seed:",/3x,"rseed = ",1p1e25.16e3)')
     &             rijkinit
            write(io,'(1x,
     &         "next initial random seed:",/3x,"rseed = ",1p1e25.16e3)')
     &           rijk
           else ! when xorshift
            write(io,'(1x,
     &         "initial random seed:",/3x,"bitrseed = ",b64.64)')
     &             rijkinit
            write(io,'(1x,
     &         "next initial random seed:",/3x,"bitrseed = ",b64.64)')
     &           rijk
           end if

      end if

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine sumcput(ierr,icc)
*                                                                      *
*       summary of cpu time                                            *
*       last modified by K.Niita and S. Hashimoto on 2011/08/02        *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     cputm(i) : time for CPU summary                                  *
*     rncnt(i) : counter for CPU summary                               *
*                                                                      *
*           i =  1, total cpu time                                     *
*                                                                      *
*             =  2, set data                                           *
*             =  3, analysis                                           *
*             =  4, nevap                                              *
*             =  5, dexgam                                             *
*             =  6, nreac                                              *
*                                                                      *
*             =  7, dklos                                              *
*             =  8, hydro                                              *
*             =  9, n-data                                             *
*             = 10, p-data                                             *
*             = 11, e-data                                             *
*             = 12, h-data                                             *
*             = 13, elast                                              *
*             = 14, ncacs                                              *
*                                                                      *
*             = 15, bertini                                            *
*             = 16, isobar                                             *
*             = 17, JAM                                                *
*             = 18, QMD                                                *
*             = 19, JAMQMD                                             *
*             = 20, INCL                                               *
*             = 21, INC-ELF                                            *
*             = 22, photonuclear                                       *
*             = 23, muon reaction                                      *
*             = 24, p-egs5                                             *
*             = 25, e-egs5                                             *
*             = 26, muon capture                                       *
*             = 27, muon bremsstrahlung                                *
*             = 28, pair production by muon                            *
*             = 29, electron track structure mode                      *
*             = 30, ion track structure mode                           *
*             = 31, frag data (user defined cross section)             *
*             = 32, SCINFUL mode                                       *
*             = 33, photonuclear library                               *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*      rnint(20*(irnmdl-1)+i) : counter for Model calculation          *
*     rnintr(20*(irnmdl-1)+i) : real counter for Model calculation     *
*                                                                      *
*      irnmdl =  1, bertini                                            *
*             =  2, isobar                                             *
*             =  3, JAM                                                *
*             =  4, QMD                                                *
*             =  5, JAMQMD                                             *
*             =  6, INCL                                               *
*             =  7, INCELF                                             *
*                                                                      *
*           i =  1, proton                                             *
*             =  2, neutron                                            *
*             =  3, pion+                                              *
*             =  4, pion0                                              *
*             =  5, pion-                                              *
*             =  6, muon+                                              *
*             =  7, muon-                                              *
*             =  8, kaon+                                              *
*             =  9, kaon0                                              *
*             = 10, kaon-                                              *
*             = 12, electron                                           *
*             = 13, positron                                           *
*             = 14, photon                                             *
*             = 15, deuteron                                           *
*             = 16, triton                                             *
*             = 17, 3He                                                *
*             = 18, Alpha                                              *
*             = 19, residual nucleus                                   *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*      rnpnt(20*(irnmdl-1)+i) : counter for Photonuclear reaction      *
*     rnpntr(20*(irnmdl-1)+i) : real counter for Photonuclear reaction *
*                                                                      *
*      irnpnr =  1, photon                                             *
*             =  2, muon                                               *
*                                                                      *
*           i =  1, Giant-dipole resonace                              *
*             =  2, Quasideuteron disintegration                       *
*             =  3, Nucleon resonance excitation                       *
*             =  4, Nucleon non-resonance excitation                   *
*                                                                      *
************************************************************************
C for USE_MOD_COUNTER
      use mod_counter, only: rncnt,rnint,rnintr,rnpnt,rnpntr
     &                      ,aevts,aevtr,bevts,bevtr
     &                      ,iaevt, ibevt, jaevt, jbevt
C for REDUCTION_COUNTER
!$   &                      ,rncnt2,rnint2,rnintr2,rnpnt2,rnpntr2
!$   &                      ,aevts2,aevtr2,bevts2,bevtr2

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      common /mpi00/ npe, me
      common /mpi01/ iccp(20000)
      common /mpi04/ stim1, stim2

*-----------------------------------------------------------------------

      common /inout/  in,io
      common /qparm/  ielas,icasc,iqstep,lvlopt,igamma

*-----------------------------------------------------------------------

      common /tcntl/  icntl, inucr

*-----------------------------------------------------------------------

      common /cputim/ stime(40), cputm(40)
      common /ptname/ pname(20), ipln(20)
      common /paraj/ mstz(300), parz(300)

*-----------------------------------------------------------------------

      common /crshi/  bplus, icrhi, ijudg, imadj, iqmax
      common /ccxsm/  icxsni, icxspi

      common /infprint/ infout
*-----------------------------------------------------------------------

      dimension stime0(40), cputm0(40)
      dimension rncnt0(40), rnint0(200), rnintr0(200)
      dimension rnpnt0(40), rnpntr0(40)
      dimension tint(19)

*-----------------------------------------------------------------------

      integer iegsemi, iegsout
      common /egsemi/ iegsemi, iegsout

! 2015/5/28 T.Sato, Datapath common
      common /paran/ icfn(100), ilfn(100), chfn(100)
      character chfn*200
      integer icfn, ilfn

*-----------------------------------------------------------------------
*     ending time
*-----------------------------------------------------------------------

      if( me .eq. 0 ) then

            call date_a_time(iyer1,imon1,iday1,
     &                       ihor1,imin1,isec1)

            write(io,'(/79(''-''))')
            write(io,'('' job termination date : '',
     &            i4,''/'',i2.2,''/'',i2.2)') iyer1, imon1, iday1

            write(io,'(''                 time :   '',
     &              i2.2,'':'',i2.2,'':'',i2.2)') ihor1, imin1, isec1

      end if

*-----------------------------------------------------------------------
*     EGS file delete
*-----------------------------------------------------------------------

      if( me .eq. 0 ) then

         if( iegsemi .ne. 0 ) then

            open(UNIT=567,file=chfn(23)(1:ilfn(23))//'.tmp',
     &       STATUS='unknown')
            close(UNIT=567,status='delete')

         end if

        if( iegsemi .ne. 0 .and. iegsout .eq. 0 ) then

            open(UNIT=512,FILE=chfn(23)(1:ilfn(23))//'.dat',
     &      STATUS='unknown')
            close(UNIT=512,status='delete')
            open(UNIT=525,FILE=chfn(23)(1:ilfn(23))//'.inp',
     &      STATUS='unknown')
            close(UNIT=525,STATUS='delete')
            open(UNIT=517,FILE=chfn(23)(1:ilfn(23))//'.msfit',
     &      STATUS='unknown')
            close(UNIT=517,STATUS='delete')
            open(UNIT=531,file=chfn(23)(1:ilfn(23))//'job.ssl',
     &      STATUS='unknown')
            close(UNIT=531,status='delete')

         end if

      end if

*-----------------------------------------------------------------------
*     error
*-----------------------------------------------------------------------

      if ( ierr .ne. 0 .or. icc .ne. 0 ) go to 900

*-----------------------------------------------------------------------
*     summary of cputime
*-----------------------------------------------------------------------

         call cputime(1)

*-----------------------------------------------------------------------
*        for parallel
*-----------------------------------------------------------------------

      if( npe .gt. 1 ) then

         if( me .gt. 0 ) then

                  call parasr(stime(1),40,0)
                  call parasr(cputm(1),40,0)
                  call parasr(rncnt(1),40,0)
                  call parasr(rnint(1),200,0)
                  call parasr(rnintr(1),200,0)
                  call parasr(rnpnt(1),40,0)
                  call parasr(rnpntr(1),40,0)

         end if

         if( me .eq. 0 ) then

               do i = 1, 40

                  stime(i) = 0.0
                  cputm(i) = 0.0
                  rncnt(i) = 0
                  rnpnt(i) = 0
                  rnpntr(i) = 0
               end do
               do i = 1, 200
                  rnint(i) = 0
                  rnintr(i) = 0
               end do

            do k = 1, npe - 1

               if( iccp(k) .eq. 0 ) then

                  call pararr(stime0(1),40,k)
                  call pararr(cputm0(1),40,k)
                  call pararr(rncnt0(1),40,k)
                  call pararr(rnint0(1),200,k)
                  call pararr(rnintr0(1),200,k)
                  call pararr(rnpnt0(1),40,k)
                  call pararr(rnpntr0(1),40,k)
                  do i = 1, 40

                     stime(i) = stime(i) + stime0(i)
                     cputm(i) = cputm(i) + cputm0(i)
                     rncnt(i) = rncnt(i) + rncnt0(i)
                     rnpnt(i) = rnpnt(i) + rnpnt0(i)
                     rnpntr(i) = rnpntr(i) + rnpntr0(i)
                  end do

                  do i = 1, 200
                     rnint(i) = rnint(i) + rnint0(i)
                     rnintr(i) = rnintr(i) + rnintr0(i)
                  end do

               end if

            end do

         end if

      end if

*-----------------------------------------------------------------------

      if( me .eq. 0 .and. cputm(1) .gt. 0.0d0 ) then

            tcpu1  = cputm( 1) / cputm(1) * 100.0
            tcpu2  = cputm( 2) / cputm(1) * 100.0
            tcpu3  = cputm( 3) / cputm(1) * 100.0
            tcpu4  = cputm( 4) / cputm(1) * 100.0
            tcpu5  = cputm( 5) / cputm(1) * 100.0
            tcpu6  = cputm( 6) / cputm(1) * 100.0
            tcpu7  = cputm( 7) / cputm(1) * 100.0
            tcpu8  = cputm( 8) / cputm(1) * 100.0
            tcpu9  = cputm( 9) / cputm(1) * 100.0
            tcpu10 = cputm(10) / cputm(1) * 100.0
            tcpu11 = cputm(11) / cputm(1) * 100.0
            tcpu12 = cputm(12) / cputm(1) * 100.0
            tcpu13 = cputm(13) / cputm(1) * 100.0
            tcpu14 = cputm(14) / cputm(1) * 100.0
            tcpu15 = cputm(15) / cputm(1) * 100.0
            tcpu16 = cputm(16) / cputm(1) * 100.0
            tcpu17 = cputm(17) / cputm(1) * 100.0
            tcpu18 = cputm(18) / cputm(1) * 100.0
            tcpu19 = cputm(19) / cputm(1) * 100.0
            tcpu20 = cputm(20) / cputm(1) * 100.0
            tcpu21 = cputm(21) / cputm(1) * 100.0
            tcpu22 = cputm(22) / cputm(1) * 100.0
            tcpu23 = cputm(23) / cputm(1) * 100.0
            tcpu24 = cputm(24) / cputm(1) * 100.0
            tcpu25 = cputm(25) / cputm(1) * 100.0
            tcpu26 = cputm(26) / cputm(1) * 100.0
            tcpu27 = cputm(27) / cputm(1) * 100.0
            tcpu28 = cputm(28) / cputm(1) * 100.0
            tcpu29 = cputm(29) / cputm(1) * 100.0
            tcpu30 = cputm(30) / cputm(1) * 100.0 ! T.Sato 2021/09/20 for avoiding uninitialization
!nais added  photonuclear library
            tcpu33 = cputm(33) / cputm(1) * 100.0

         if( ( icntl .eq. 0 .or. icntl .eq. 5 .or.
     &         icntl .eq. 6 .or.
     &         icntl .eq. 14 .or.
     &         icntl .eq. 15 ) .and. mstz(44) .ne. 0) then

            cputt  = cputm( 1) - cputm( 2) - cputm( 3) - cputm( 4)
     &             - cputm( 5) - cputm( 6)
            tcput  = cputt / cputm(1) * 100.0

            cpute  = cputm( 6) - cputm( 7) - cputm( 8) - cputm( 9)
     &             - cputm(10) - cputm(11) - cputm(12) - cputm(13)
     &             - cputm(14) - cputm(21) - cputm(22)
     &             - cputm(24) - cputm(25) - cputm(26) - cputm(27)
!
!nais added photonuclear library
     &             - cputm(28) -cputm(33)
            tcpue  = cpute / cputm(1) * 100.0

         else

            cputt  = 0.0
            tcput  = 0.0

            cpute  = 0.0
            tcpue  = 0.0

         end if

            write(io,'(/,79(''-''))')
            write(io,'(''    CPU time and number of event called in ''
     &               ''PHITS'')')
            write(io,'(79(''-''))')

         if( mstz(44) .ne. 0 ) then

            write(io,'( /16x,''         sec          %'',
     &                       ''             count'')')
            write(io,'( ''total cpu time ='',2f13.2)') cputm(1), tcpu1
            write(io,'(/''     transport ='',2f13.2)') cputt, tcput
            write(io,'( ''      set data ='',2f13.2)') cputm(2), tcpu2

            write(io,'( ''      analysis ='',2f13.2,f16.0)')
     &                            cputm(3), tcpu3, rncnt(3)

            write(io,'( ''         nevap ='',2f13.2,f16.0)')
     &                            cputm(4), tcpu4, rncnt(4)

            write(io,'( ''        dexgam ='',2f13.2,f16.0)')
     &                            cputm(5), tcpu5, rncnt(5)

            write(io,'( ''         nreac ='',2f13.2,f16.0)')
     &                            cputm(6), tcpu6, rncnt(6)

            write(io,'(/''         other ='',2f13.2)')
     &                            cpute, tcpue
            write(io,'( ''         dklos ='',2f13.2,f16.0)')
     &                            cputm(7), tcpu7, rncnt(7)
            write(io,'( ''         hydro ='',2f13.2,f16.0)')
     &                            cputm(8), tcpu8, rncnt(8)
            write(io,'( ''        n-data ='',2f13.2,f16.0)')
     &                            cputm(9), tcpu9, rncnt(9)
            write(io,'( ''        h-data ='',2f13.2,f16.0)')
     &                            cputm(12), tcpu12, rncnt(12)
            write(io,'( ''        p-data ='',2f13.2,f16.0)')
     &                            cputm(10), tcpu10, rncnt(10)
            write(io,'( ''        e-data ='',2f13.2,f16.0)')
     &                            cputm(11), tcpu11, rncnt(11)
            write(io,'( ''        p-egs5 ='',2f13.2,f16.0)')
     &                            cputm(24), tcpu24, rncnt(24)
            write(io,'( ''        e-egs5 ='',2f13.2,f16.0)')
     &                            cputm(25), tcpu25, rncnt(25)
            write(io,'( ''      e-tsmode ='',2f13.2,f16.0)')   ! Takeshi Kai
     &                            cputm(29), tcpu29, rncnt(29) ! Takeshi Kai
            write(io,'( ''    ion-tsmode ='',2f13.2,f16.0)')   ! Takeshi Kai
     &                            cputm(30), tcpu30, rncnt(30) ! Takeshi Kai
            write(io,'( ''     photonucl ='',2f13.2,f16.0)')
     &                            cputm(22), tcpu22, rncnt(22)
!nais added photonuclear library
            write(io,'( '' photonucl lib ='',2f13.2,f16.0)')
     &                            cputm(33), tcpu33, rncnt(33)
!
            write(io,'( ''   muon p-nucl ='',2f13.2,f16.0)')
     &                            cputm(23), tcpu23, rncnt(23)
            write(io,'( ''   muon bremss ='',2f13.2,f16.0)')
     &                            cputm(27), tcpu27, rncnt(27)
            write(io,'( ''    muon pprod ='',2f13.2,f16.0)')
     &                            cputm(28), tcpu28, rncnt(28)
            write(io,'( ''  muon capture ='',2f13.2,f16.0)')
     &                            cputm(26), tcpu26, rncnt(26)
            write(io,'( ''         elast ='',2f13.2,f16.0)')
     &                            cputm(13), tcpu13, rncnt(13)
            write(io,'( ''         ncasc ='',2f13.2,f16.0)')
     &                            cputm(14), tcpu14, rncnt(14)

            write(io,'(/''       bertini ='',2f13.2,f16.0)')
     &                            cputm(15), tcpu15, rncnt(15)
            write(io,'( ''        isobar ='',2f13.2,f16.0)')
     &                            cputm(16), tcpu16, rncnt(16)
            write(io,'( ''           JAM ='',2f13.2,f16.0)')
     &                            cputm(17), tcpu17, rncnt(17)
            write(io,'( ''           QMD ='',2f13.2,f16.0)')
     &                            cputm(18), tcpu18, rncnt(18)
            write(io,'( ''        JAMQMD ='',2f13.2,f16.0)')
     &                            cputm(19), tcpu19, rncnt(19)
            write(io,'( ''          INCL ='',2f13.2,f16.0)')
     &                            cputm(20), tcpu20, rncnt(20)
            write(io,'( ''        INCELF ='',2f13.2,f16.0)')
     &                            cputm(21), tcpu21, rncnt(21)

         else

            write(io,'( /16x,''             sec'')')
            write(io,'( ''total cpu time ='',f16.2)') cputm(1)

            write(io,'( //16x,''           count'')')

            write(io,'( ''      analysis ='',f16.0
     &              ,''  : data processing'')') rncnt(3)
            write(io,'( ''         nevap ='',f16.0
     &              ,''  : evaporation'')') rncnt(4)
            write(io,'( ''        dexgam ='',f16.0
     &              ,''  : de-excitation'')') rncnt(5)
            write(io,'( ''         nreac ='',f16.0
     &              ,''  : atomic and nuclear reactions'')') rncnt(6)
            write(io,'( ''         dklos ='',f16.0
     &              ,''  : particle decay'')') rncnt(7)
            write(io,'( ''         hydro ='',f16.0
     &              ,''  : nucleon-nucleon scattering'')') rncnt(8)
            write(io,'( ''        n-data ='',f16.0
     &              ,''  : neutron data library'')') rncnt(9)
            write(io,'( ''        h-data ='',f16.0
     &              ,''  : p, d, a data library'')') rncnt(12)
            write(io,'( ''        p-data ='',f16.0
     &              ,''  : photon data library'')') rncnt(10)
            write(io,'( ''        e-data ='',f16.0
     &              ,''  : electron data library'')') rncnt(11)
            write(io,'( ''        p-egs5 ='',f16.0
     &              ,''  : photon interaction with EGS5'')') rncnt(24)
            write(io,'( ''        e-egs5 ='',f16.0
     &              ,''  : electron interaction with EGS5'')') rncnt(25)
            write(io,'( ''      e-tsmode ='',f16.0                       ! Takeshi Kai
     &              ,''  : electron track structure mode'')') rncnt(29)  ! Takeshi Kai
            write(io,'( ''    ion-tsmode ='',f16.0                       ! Takeshi Kai
     &              ,''  : ion track structure mode'')') rncnt(30)       ! Takeshi Kai
            write(io,'( ''     photonucl ='',f16.0
     &              ,''  : photo-nuclear reaction'')') rncnt(22)
!nais added photonuclear library
            write(io,'( '' photonucl lib ='',f16.0
     &              ,''  : photo-nuclear reaction with library'')')
     &               rncnt(33)
            write(io,'( ''       mu-reac ='',f16.0
     &              ,''  : muon-induced nuclear reaction'')') rncnt(23)
            write(io,'( ''       mu-brem ='',f16.0
     &              ,''  : muon-induced bremsstrahlung'')') rncnt(27)
            write(io,'( ''       mu-pprd ='',f16.0
     &              ,''  : muon-induced pair production'')') rncnt(28)
            write(io,'( ''        mu-cap ='',f16.0
     &              ,''  : muon capture in nucleus'')') rncnt(26)
            write(io,'( ''         elast ='',f16.0
     &              ,''  : elastic scattering'')') rncnt(13)
            write(io,'( ''         ncasc ='',f16.0
     &              ,''  : nuclear reaction model'')') rncnt(14)

            write(io,'(/''       bertini ='',f16.0
     &              ,''  : Bertini model'')') rncnt(15)
            write(io,'( ''        isobar ='',f16.0
     &              ,''  : isobar model'')') rncnt(16)
            write(io,'( ''           JAM ='',f16.0
     &              ,''  : JAM model'')') rncnt(17)
            write(io,'( ''           QMD ='',f16.0
     &              ,''  : JQMD model'')') rncnt(18)
            write(io,'( ''        JAMQMD ='',f16.0
     &              ,''  : JAMQMD model'')') rncnt(19)
            write(io,'( ''          INCL ='',f16.0
     &              ,''  : INCL model'')') rncnt(20)
            write(io,'( ''        INCELF ='',f16.0
     &              ,''  : INCELF model'')') rncnt(21)
            write(io,'( ''     frag data ='',f16.0
     &              ,''  : user defined cross section'')')rncnt(31)
            write(io,'( ''       SCINFUL ='',f16.0
     &              ,''  : SCINFUL mode'')') rncnt(32)

         end if

*-----------------------------------------------------------------------

         if ( infout .eq. 8 ) then

*-----------------------------------------------------------------------
*        Model calculation summary
*-----------------------------------------------------------------------

            do i = 1, 7

               irnmdl = 20 * (i-1)
               rnina = 0.0d0
               rninra = 0.0d0

               do j = 1, 19
                  rnina  = rnina + rnint(irnmdl+j)
                  rninra  = rninra + rnintr(irnmdl+j)
               enddo

               if( rnina .gt. 0.0d0 ) then

                  tinta = rninra / rnina * 100.0d0

                  if ( i .eq. 1 )
     &            write(io,'(/''   === incident particle into'',
     &                        '' bert ====='')')
                  if ( i .eq. 2 )
     &            write(io,'(/''   === incident particle into'',
     &                        '' isobar ====='')')
                  if ( i .eq. 3 )
     &            write(io,'(/''   === incident particle into'',
     &                        '' JAM ====='')')
                  if ( i .eq. 4 )
     &            write(io,'(/''   === incident particle into'',
     &                        '' QMD ====='')')
                  if ( i .eq. 5 )
     &            write(io,'(/''   === incident particle into'',
     &                        '' JAMQMD ====='')')
                  if ( i .eq. 6 )
     &            write(io,'(/''   === incident particle into'',
     &                        '' INCL ====='')')
                  if ( i .eq. 7 )
     &            write(io,'(/''   === incident particle into'',
     &                        '' INC-ELF ====='')')

                  write(io,'(38x,''%'',13x,''count'',12x,''real'')')

                  do j = 1, 20
                     if( rnint(irnmdl+j) .gt. 0.0d0 ) then
                       tint(j) = rnintr(irnmdl+j) / rnint(irnmdl+j)
     &                         * 100.0d0

                       write(io,'(4x,a8,'' ='',15x,f13.2,2f16.0)')
     &                 pname(j),tint(j),rnint(irnmdl+j),rnintr(irnmdl+j)
                     endif
                  enddo

                  write(io,'(3x,71(''-''))')
                  write(io,'(20x,''total ='',2x,f13.2,2f16.0)')
     &                  tinta, rnina, rninra

               endif

            enddo

*-----------------------------------------------------------------------
*        Photonuclear reaction summary
*-----------------------------------------------------------------------

            do i = 1, 2

!!
               irnpnr = 5 * (i-1)
               rnpna = 0.0d0
               rnpnra = 0.0d0

               do j = 1, 5
                  rnpna  = rnpna + rnpnt(irnpnr+j)
                  rnpnra  = rnpnra + rnpntr(irnpnr+j)
               enddo

               if( rnpna .gt. 0.0d0 ) then

                  tinta = rnpnra / rnpna * 100.0d0

                  if ( i .eq. 1 )
     &            write(io,'(/''   === photonuclear interaction'',
     &                        '' caused by photon ====='')')
                  if ( i .eq. 2 )
     &            write(io,'(/''   === photonuclear interaction'',
     &                        '' caused by muon ====='')')

                  write(io,'(38x,''%'',13x,''count'',12x,''real'')')

                  do j = 1, 5
                     if( rnpnt(irnpnr+j) .gt. 0.0d0 ) then
                        tint(j) = rnpntr(irnpnr+j) / rnpnt(irnpnr+j)
     &                          * 100.0d0

                        if ( j .eq. 1 )
     &                  write(io,'(4x,''Giant-dipole resonace'',
     &                                '' ='',2x,f13.2,2f16.0)')
     &                  tint(j),rnpnt(irnpnr+j),rnpntr(irnpnr+j)
                        if ( j .eq. 2 )
     &                  write(io,'(4x,''        Quasideuteron'',
     &                                '' ='',2x,f13.2,2f16.0)')
     &                  tint(j),rnpnt(irnpnr+j),rnpntr(irnpnr+j)
                        if ( j .eq. 3 )
     &                  write(io,'(4x,''    Nucleon resonance'',
     &                                '' ='',2x,f13.2,2f16.0)')
     &                  tint(j),rnpnt(irnpnr+j),rnpntr(irnpnr+j)
                        if ( j .eq. 4 )
     &                  write(io,'(4x,''Nucleon non-resonance'',
     &                                '' ='',2x,f13.2,2f16.0)')
     &                  tint(j),rnpnt(irnpnr+j),rnpntr(irnpnr+j)
                        if ( j .eq. 5 )
     &                  write(io,'(4x,''         used library'',
     &                                '' ='',2x,f13.2,2f16.0)')
     &                  tint(j),rnpnt(irnpnr+j),rnpntr(irnpnr+j)

                     endif
                  enddo

                  write(io,'(3x,71(''-''))')
                  write(io,'(20x,''total ='',2x,f13.2,2f16.0)')
     &                  tinta, rnpna, rnpnra

               endif

            enddo

         endif

*-----------------------------------------------------------------------
*     print references when special models are used
*-----------------------------------------------------------------------

         write(io,'(/''>>> Citation Request >>>'',55(''=''))')

         if ( rncnt(20).ge.1 .or. rncnt(21).ge.1 .or.
     &        icrhi.eq.2 .or. icxsni.eq.1 .or. iegsemi .ne. 0 ) then

         write(io,'(
     &   /'' This execution uses model(s) that must be explicitly ''
     &   ,''cited in addition to'',/1x''the PHITS original document: ''
     &   ,''T.Sato et al., J.Nucl.Sci.Technol.61, 127-135 (2024).''
     &   ,/'' Please refer the following document(s) ''
     &   ,''in your publication using this result'')')

         if ( rncnt(20) .ge. 1 ) then

         write(io,'(/'' The INCL model:''
     &   ,/4x''A. Boudard et al., Phys. Rev C87, 014606 (2013).'')')

         end if

         if ( iegsemi .ne. 0 ) then

         write(io,'(/'' The EGS5 code:''
     &   ,/4x''H. Hirayama et al., SLAC-R-730 (2005) and KEK ''
     &   ,''Report 2005-8 (2005)'')')

         end if

         if ( rncnt(29) .gt. 0 ) then

         write(io,'(/'' e-tsmode:''
     &   ,/4x''T. Kai et al., Radiat. Phys. Chem., 115, 1-5 (2015).
     &   '')')

         end if

         if ( rncnt(21) .ge. 1 ) then

         write(io,'(/'' The INC-ELF model:''
     &   ,/4x''Y. Sawada et al., Nucl. Instr. & Meth. B 291, 38-44''
     &   ,'' (2012)'')')

         end if

         if ( icrhi .eq. 2 .or. icxsni .eq. 1 ) then

         write(io,'(/'' The KUROTAMA model:''
     &            ,/4x,''K. Iida, A. Kohama, and K. Oyamatsu,''
     &            ,'' J. Phys. Soc. Japan 76, 044201 (2007).'')')

         end if

         else

         write(io,'(
     &   /'' Please refer the following document(s) ''
     &   ,''in your publication using this result''
     &   ,/'' T.Sato et al., J.Nucl.Sci.Technol.61, 127-135 (2024).'')')

         end if

*-----------------------------------------------------------------------

      end if

*-----------------------------------------------------------------------
*     end of all job and  MPI finalization
*-----------------------------------------------------------------------

      if( npe .gt. 1  ) then

            stim2 = paratim()
            stim3 = stim2 - stim1

         if( me .gt. 0 ) then

                     call parasr(stim3,1,0)

         else

                     stim5 = 0.0

               do i = 1, npe - 1

                  if( iccp(i) .eq. 0 ) then

                     call pararr(stim4,1,i)

                     stim5 = stim5 + stim4

                  end if

               end do

               write(io,'(/79(''-''))')
               write(io,'('' Final Parallel Status:'',
     &                    '' 0-> normal, 1-> abnormal stop'')')
               write(io,'( 79(''-''))')

               write(io,'('' ip  status'')')
               write(io,'(2i5)') (i,iccp(i),i=1,npe-1)

               ntpe = 0

            do i = 1, npe - 1

               if( iccp(i) .eq. 0 ) ntpe = ntpe + 1

            end do

               write(io,'(/'' total executable PE ='',i13)') ntpe

               write(io,'( '' elapse time of root ='',f13.2)')
     &                        stim3

               write(io,'( '' average of cpu time ='',f13.2)')
     &                        cputm(1) / dble( ntpe )

               write(io,'( '' sum of elapse time  ='',f13.2)')
     &                        stim5

               write(io,'( '' sum of cpu time     ='',f13.2)')
     &                        cputm(1)

               actm = stim5 / stim3
               write(io,'( '' accel: ela-all/ela0 ='',f13.2)')
     &                        actm

               actm = cputm(1) / stim3
               write(io,'( '' accel: cpu-all/ela0 ='',f13.2)')
     &                        actm

               actm = cputm(1) / stim3 / dble( ntpe )
               write(io,'( '' acceleration rate   ='',f13.4)')
     &                        actm

         end if

      end if

*-----------------------------------------------------------------------

  900 continue

      if( me .eq. 0  ) then

               write(io,'(/'' END '')')

      end if

               close(io)

*-----------------------------------------------------------------------

      return
      end

