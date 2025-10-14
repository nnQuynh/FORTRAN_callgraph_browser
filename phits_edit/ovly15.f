************************************************************************
*                                                                      *
      subroutine ovly15
*                                                                      *
*       control routine for dump all calculation                       *
*       last modified by S. Hashimoto and K.Niita on 2011/07/06        *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*       ncol = 1 : start of calculation                                *
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
*             13 : (n,x) reaction                                      *
*             14 : (n,n'x) reaction                                    *
*             15 : sequential transport only for tally                 *
*                                                                      *
*----------------------------------------------------------------------*
*        jcoll reaction type identifier                                *
*----------------------------------------------------------------------*
*                                                                      *
*        jcoll : =  0, nothing happen                                  *
*                =  1, Hydrogen collisions                             *
*                =  2, Particle Decays                                 *
*                =  3, Elastic collisions                              *
*                =  4, High Energy Nuclear collisions                  *
*                =  5, Heavy Ion reactions                             *
*                =  6, Neutron reactions by data                       *
*                =  7, Photon reactions by data                        *
*                =  8, Electron reactions by data                      *
*                =  9, P,d,a, and photo-nuclear reactions by data      *
*                = 10, Neutron event mode                              *
*                = 11, Delta Ray production                            *
*                = 12, Muon atomic interaction                         *
*                = 13, Photon by EGS5                                  *
*                = 14, Electron by EGS5                                *
*                = 15, Photon photonuclear interaction                 *
*                = 16, Negative muon captured by nucleon               *
*                = 17, Muon photonuclear interaction                   *
*                = 18, Electron recoil by track strcuture mode         *
*                = 19, Muon pair production (photon -> mu+ mu-)        *
*                = 20, User defined interaction                        *
*                                                                      *
*----------------------------------------------------------------------*
*        kcoll reaction type identifier                                *
*----------------------------------------------------------------------*
*                                                                      *
*        kcoll : =  0, normal                                          *
*                =  1, high energy fission                             *
*                =  2, high energy absorption                          *
*                =  3, low energy n elastic                            *
*                =  4, low energy n non-elastic                        *
*                =  5, low energy n fission                            *
*                =  6, low energy n absorption                         *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*        nclsts   : total number of out going particles and nuclei     *
*                                                                      *
*        iclusts(nclsts)                                               *
*                                                                      *
*                i = 0, nucleus                                        *
*                  = 1, proton                                         *
*                  = 2, neutron                                        *
*                  = 3, pion                                           *
*                  = 4, photon                                         *
*                  = 5, kaon                                           *
*                  = 6, muon                                           *
*                  = 7, others                                         *
*                                                                      *
*        jclusts(i,nclsts)                                             *
*                                                                      *
*                i = 0, angular momentum                               *
*                  = 1, proton number                                  *
*                  = 2, neutron number                                 *
*                  = 3, ip, see below                                  *
*                  = 4, status of the particle 0: real, <0 : dead      *
*                  = 5, charge                                         *
*                  = 6, baryon number                                  *
*                  = 7, kf code                                        *
*                  = 8, isomer level (0: Ground, 1,2: 1st, 2nd isomer) *
*                                                                      *
*        qclusts(i,nclsts)                                             *
*                                                                      *
*                i = 0, impact parameter                               *
*                  = 1, x-component of unit vector of momentum         *
*                  = 2, y-component of unit vector of momentum         *
*                  = 3, z-component of unit vector of momentum         *
*                  = 4, etot = sqrt( p**2 + rm**2 ) (GeV)              *
*                  = 5, rest mass (GeV)                                *
*                  = 6, excitation energy (MeV)                        *
*                  = 7, kinetic energy (MeV)                           *
*                  = 8, weight                                         *
*                  = 9, time                                           *
*                  = 10, x                                             *
*                  = 11, y                                             *
*                  = 12, z                                             *
*                                                                      *
************************************************************************
C for USE_MOD_COUNTER !FURUTA20210506
      use mod_counter, only: ALLOCATE_EVTS,DEALLOCATE_EVTS,INIT_EVTS
      use MMBANKMOD  !FURUTA
      use MEMBANKMOD !FURUTA
      use GGBANKMOD  !FURUTA
      use GGMBANKMOD !FURUTA
      use EVENTTALMOD!FURUTA
      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param00.inc'
      include 'param.inc'

      include 'err.inc'

*-----------------------------------------------------------------------

      common /mpi00/ npe, me

*-----------------------------------------------------------------------

      common /jcomon/ nabov,nobch,nocas,nomax
!$OMP THREADPRIVATE(/jcomon/)
      common /rcomon/ rcasc

      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
!$OMP THREADPRIVATE(/icomon/)
      common /wtsave/ oldwt
!$OMP THREADPRIVATE(/wtsave/)
      common /tlgeom/ iblz1,iblz2
!$OMP THREADPRIVATE(/tlgeom/)
      common /tlglat/ ilev1,ilev2,ilat1(5,10),ilat2(5,10)
!$OMP THREADPRIVATE(/tlglat/)
      common /tlcost/ costha, uang(3), nsurf
!$OMP THREADPRIVATE(/tlcost/)

      common /elect/  uint(3), qs, qo, eint, delc, am, qsex,
     &                nq, ns, n1, noz, mtel
!$OMP THREADPRIVATE(/elect/)

      common /mathzn/ mathz, mathn, jcoll, kcoll
!$OMP THREADPRIVATE(/mathzn/)

      common /clustt/ nclsts, iclusts(nnn)
!$OMP THREADPRIVATE(/clustt/)
      common /clustw/ jclusts(0:8,nnn),  qclusts(0:12,nnn)
!$OMP THREADPRIVATE(/clustw/)
      common /cntcls/ jcount(3,nnn)
!$OMP THREADPRIVATE(/cntcls/)

      common /taliin/ rsouin, nzztin, nrgnin
      common /kmat1d/ idmn(0:kvlmax), idnm(kvmmax)

*-----------------------------------------------------------------------


      integer nrandgen
      common /randn/ nrandgen ! S.H. xorshift (2020.5.29)
      integer*8 :: iranji64 ! S.H. xorshift (2020.5.29)
      common /randm4/ rnfb,rnfs,rngb,rngs,rnmult,ranj,rani,
     &                rnrtc,nstrid,inif, iranji64
      integer*8 :: iransb64 ! S.H. xorshift (2020.2.6)
      common /randtp/ rijk,rans,ranb, iransb64
!$OMP THREADPRIVATE(/randtp/)

      common /dumpall/ idumpall
      common /ccggg/ icgg

! T.Sato 2016/05/28, for batch.out
      common /paran/ icfn(100), ilfn(100), chfn(100)
      character chfn*200
      integer icfn,ilfn

*-----------------------------------------------------------------------
*     write id is 9
*-----------------------------------------------------------------------

            io = 9

*-----------------------------------------------------------------------
*        initialize the 'batch.out' file
*-----------------------------------------------------------------------

         if( me .eq. 0 ) then

               iot = 27
               open(iot,file=chfn(22),status='unknown')

               write(iot,'(''1 <--- 1:continue, 0:stop'')')

               write(iot,'(/79(''-'')/'' start calculation''
     &                     /79(''-''))')

               call date_a_time(iyer0,imon0,iday0,
     &                          ihor0,imin0,isec0)
               write(iot,'(/'' date = '',
     &                            i4,''-'',i2.2,''-'',i2.2)')
     &                            iyer0,imon0,iday0
               write(iot,'( '' time = '',
     &                            i2.2,''h '',i2.2,''m '',i2.2/)')
     &                            ihor0,imin0,isec0

               close(iot)

         end if

*-----------------------------------------------------------------------
*        for control PE
*-----------------------------------------------------------------------

         if( me .eq. 0 .and. npe .gt. 1 ) then

               if(icgg.ne.0)then
                 call ALLOCATE_GGBANK !FURUTA
                 call INIT_GGBANK     !FURUTA
               endif
               call ALLOCATE_GGMBANK  !FURUTA
               call INIT_GGMBANK      !FURUTA

               ncol  = 3

               call analyz(ncol,mark)

               ncol  = 2

               call analyz(ncol,mark)

               if(icgg.ne.0) call DEALLOCATE_GGBANK !FURUTA
               call DEALLOCATE_GGMBANK              !FURUTA

               return

         end if

*-----------------------------------------------------------------------
*        random number generator
*-----------------------------------------------------------------------

               if( inif .eq. 0 ) call advijk

               inif = 0
            if ( nrandgen .eq. 0 ) then ! S.H. xorshift (2020.5.29)
               ranb = rani
               rans = ranj
            else
               iransb64 = iranji64
            end if

*-----------------------------------------------------------------------
      call ALLOCATE_MMBANK      !FURUTA
      call ALLOCATE_MEMBANK     !FURUTA
      if(icgg.ne.0)then
        call ALLOCATE_GGBANK    !FURUTA
        call INIT_GGBANK        !FURUTA
      endif
      call ALLOCATE_GGMBANK     !FURUTA
      call INIT_GGMBANK         !FURUTA
      call ALLOCATE_EVENTTAL    !FURUTA
      call ALLOCATE_EVTS        !FURUTA20210506
      call INIT_EVTS            !FURUTA20210506
      if( idumpall .eq. 1 ) then

*-----------------------------------------------------------------------
*     read dumpall (binary) file
*-----------------------------------------------------------------------

  100    continue

               read(io,iostat=ios,err=999) ncol
               if( ios .eq. -1 ) goto 900

*-----------------------------------------------------------------------

         if( ncol .eq. 1 .or. ncol .eq. 2 .or. ncol .eq. 3 ) then

               if( ncol.eq.3)then
                 call EVENTanalyz        !FURUTA
               endif

               call analyz(ncol,mark)

               if( ncol .eq. 2 )then                   !FURUTA
                 call DEALLOCATE_MMBANK                !FURUTA
                 call DEALLOCATE_MEMBANK               !FURUTA
                 if(icgg.ne.0) call DEALLOCATE_GGBANK  !FURUTA
                 call DEALLOCATE_GGMBANK               !FURUTA
                 call DEALLOCATE_EVENTTAL              !FURUTA
                 call DEALLOCATE_EVTS                  !FURUTA20210506
                 return                                !FURUTA
               endif                                   !FURUTA
               goto 100

*-----------------------------------------------------------------------
*        ncol = 4 - 15
*-----------------------------------------------------------------------

         else if( ncol .ge. 4 ) then

*-----------------------------------------------------------------------

            if( ncol .eq. 4 ) then

               read(io,iostat=ios,err=999) nocas, nobch, rcasc, rsouin
               if( ios .eq. -1 ) goto 900

            end if

*-----------------------------------------------------------------------

               read(io,iostat=ios,err=999)
     &                  no,mat,ityp,ktyp,jtyp,mtyp,rtyp, oldwt
               if( ios .eq. -1 ) goto 900

               mat = idnm(mat)

            if( ityp .eq. 12 .or. ityp .eq. 13 ) then

               read(io,iostat=ios,err=999) qs
               if( ios .eq. -1 ) goto 900

               mtel = mat

            end if

*-----------------------------------------------------------------------

               read(io,iostat=ios,err=999) iblz1,iblz2, ilev1,ilev2
               if( ios .eq. -1 ) goto 900

            if( ilev1 .gt. 0 ) then

               read(io,iostat=ios,err=999)
     &                  ( ( ilat1(i,j), i=1,5 ), j=1,ilev1 )
               if( ios .eq. -1 ) goto 900

            end if

            if( ilev2 .gt. 0 ) then

               read(io,iostat=ios,err=999)
     &                  ( ( ilat2(i,j), i=1,5 ), j=1,ilev2 )
               if( ios .eq. -1 ) goto 900

            end if

              read(io,iostat=ios,err=999)
     &             costha, uang(1), uang(2), uang(3), nsurf
              if( ios .eq. -1 ) goto 900

*-----------------------------------------------------------------------
               read(io,iostat=ios,err=999)
     &                  name(ibknam+no,ipomp+1),
     &                         ( ncnt(ibknct+i,no,ipomp+1), i=1,3 )
               if( ios .eq. -1 ) goto 900

               read(io,iostat=ios,err=999)
     &                  wt(ibkwt+no,ipomp+1),
     &                  u(ibku+no,ipomp+1),v(ibkv+no,ipomp+1),
     &                  w(ibkw+no,ipomp+1)
               if( ios .eq. -1 ) goto 900
               read(io,iostat=ios,err=999)
     &                  e(ibke+no,ipomp+1),t(ibkt+no,ipomp+1),
     &                  x(ibkx+no,ipomp+1),y(ibky+no,ipomp+1),
     &                  z(ibkz+no,ipomp+1)
               if( ios .eq. -1 ) goto 900
               read(io,iostat=ios,err=999)
     &                  ec(ibkec+no,ipomp+1),tc(ibktc+no,ipomp+1),
     &                  xc(ibkxc+no,ipomp+1),yc(ibkyc+no,ipomp+1),
     &                  zc(ibkzc+no,ipomp+1)
               if( ios .eq. -1 ) goto 900
               read(io,iostat=ios,err=999)
     &                  spx(ibkspx+no,ipomp+1),spy(ibkspy+no,ipomp+1),
     &                  spz(ibkspz+no,ipomp+1)
               if( ios .eq. -1 ) goto 900
               read(io,iostat=ios,err=999)
     &                  nzst(ibkzst+no,ipomp+1)
               if( ios .eq. -1 ) goto 900

*-----------------------------------------------------------------------

            if( ncol .eq. 13 .or. ncol .eq. 14 ) then

               read(io,iostat=ios,err=999) nclsts
               if( ios .eq. -1 ) goto 900

               if( nclsts .gt. 0 ) then

                  read(io,iostat=ios,err=999) mathz, mathn, jcoll, kcoll
                  if( ios .eq. -1 ) goto 900

               do i = 1, nclsts

                  read(io,iostat=ios,err=999) iclusts(i)
                  if( ios .eq. -1 ) goto 900
                  read(io,iostat=ios,err=999) ( jclusts(j,i), j=0,8)
                  if( ios .eq. -1 ) goto 900
                  read(io,iostat=ios,err=999) ( qclusts(j,i), j=0,12)
                  if( ios .eq. -1 ) goto 900
                  read(io,iostat=ios,err=999) ( jcount(j,i),  j=1,3)
                  if( ios .eq. -1 ) goto 900

               end do

               end if

            end if

*-----------------------------------------------------------------------

               call analyz(ncol,mark)

               goto 100

         end if

*-----------------------------------------------------------------------
*     read dumpall (ascii) file
*-----------------------------------------------------------------------

      else if( idumpall .eq. -1 ) then

  200    continue

               read(io,*,iostat=ios,err=999) ncol
               if( ios .eq. -1 ) goto 900

*-----------------------------------------------------------------------

         if( ncol .eq. 1 .or. ncol .eq. 2 .or. ncol .eq. 3 ) then

               if( ncol.eq.3)then
                 call EVENTanalyz        !FURUTA
               endif

               call analyz(ncol,mark)

               if( ncol .eq. 2 )then                   !FURUTA
                 call DEALLOCATE_MMBANK                !FURUTA
                 call DEALLOCATE_MEMBANK               !FURUTA
                 if(icgg.ne.0) call DEALLOCATE_GGBANK  !FURUTA
                 call DEALLOCATE_GGMBANK               !FURUTA
                 call DEALLOCATE_EVENTTAL              !FURUTA
                 call DEALLOCATE_EVTS                  !FURUTA20210506
                 return                                !FURUTA
               endif                                   !FURUTA
               goto 200

*-----------------------------------------------------------------------
*        ncol = 4 - 15
*-----------------------------------------------------------------------

         else if( ncol .ge. 4 ) then

*-----------------------------------------------------------------------

            if( ncol .eq. 4 ) then

               read(io,*,iostat=ios,err=999)
     &              nocas, nobch, rcasc, rsouin
               if( ios .eq. -1 ) goto 900

            end if

*-----------------------------------------------------------------------

               read(io,*,iostat=ios,err=999)
     &                  no,mat,ityp,ktyp,jtyp,mtyp,rtyp, oldwt
               if( ios .eq. -1 ) goto 900

               mat = idnm(mat)

            if( ityp .eq. 12 .or. ityp .eq. 13 ) then

               read(io,*,iostat=ios,err=999) qs
               if( ios .eq. -1 ) goto 900

               mtel = mat

            end if

*-----------------------------------------------------------------------

               read(io,*,iostat=ios,err=999) iblz1,iblz2, ilev1,ilev2
               if( ios .eq. -1 ) goto 900

            if( ilev1 .gt. 0 ) then

               read(io,*,iostat=ios,err=999)
     &                  ( ( ilat1(i,j), i=1,5 ), j=1,ilev1 )
               if( ios .eq. -1 ) goto 900

            end if

            if( ilev2 .gt. 0 ) then

               read(io,*,iostat=ios,err=999)
     &                  ( ( ilat2(i,j), i=1,5 ), j=1,ilev2 )
               if( ios .eq. -1 ) goto 900

            end if

              read(io,*,iostat=ios,err=999)
     &             costha, uang(1), uang(2), uang(3), nsurf
              if( ios .eq. -1 ) goto 900

*-----------------------------------------------------------------------
               read(io,*,iostat=ios,err=999)
     &                  name(ibknam+no,ipomp+1),
     &                       ( ncnt(ibknct+i,no,ipomp+1), i=1,3 )
               if( ios .eq. -1 ) goto 900

               read(io,*,iostat=ios,err=999)
     &                  wt(ibkwt+no,ipomp+1),
     &                  u(ibku+no,ipomp+1),v(ibkv+no,ipomp+1),
     &                  w(ibkw+no,ipomp+1)
               if( ios .eq. -1 ) goto 900
               read(io,*,iostat=ios,err=999)
     &                  e(ibke+no,ipomp+1),t(ibkt+no,ipomp+1),
     &                  x(ibkx+no,ipomp+1),y(ibky+no,ipomp+1),
     &                  z(ibkz+no,ipomp+1)
               if( ios .eq. -1 ) goto 900
               read(io,*,iostat=ios,err=999)
     &                  ec(ibkec+no,ipomp+1),tc(ibktc+no,ipomp+1),
     &                  xc(ibkxc+no,ipomp+1),yc(ibkyc+no,ipomp+1),
     &                  zc(ibkzc+no,ipomp+1)
               if( ios .eq. -1 ) goto 900
               read(io,*,iostat=ios,err=999)
     &                  spx(ibkspx+no,ipomp+1),spy(ibkspy+no,ipomp+1),
     &                  spz(ibkspz+no,ipomp+1)
               if( ios .eq. -1 ) goto 900
               read(io,*,iostat=ios,err=999)
     &                  nzst(ibkzst+no,ipomp+1)
               if( ios .eq. -1 ) goto 900

*-----------------------------------------------------------------------

            if( ncol .eq. 13 .or. ncol .eq. 14 ) then

               read(io,*,iostat=ios,err=999) nclsts
               if( ios .eq. -1 ) goto 900

               if( nclsts .gt. 0 ) then

                  read(io,*,iostat=ios,err=999)
     &                 mathz, mathn, jcoll, kcoll
                  if( ios .eq. -1 ) goto 900

               do i = 1, nclsts

                  read(io,*,iostat=ios,err=999) iclusts(i)
                  if( ios .eq. -1 ) goto 900
                  read(io,*,iostat=ios,err=999) ( jclusts(j,i), j=0,8)
                  if( ios .eq. -1 ) goto 900
                  read(io,*,iostat=ios,err=999) ( qclusts(j,i), j=0,12)
                  if( ios .eq. -1 ) goto 900
                  read(io,*,iostat=ios,err=999) ( jcount(j,i),  j=1,3)
                  if( ios .eq. -1 ) goto 900

               end do

               end if

            end if

*-----------------------------------------------------------------------

               call analyz(ncol,mark)

               goto 200

         end if

      end if

*-----------------------------------------------------------------------

  900 continue

         write(ErrCha,*) ' Error: dumpall file is ended abnormally'
         ErrID = 'L:588/R:ovly15/F:ovly15.f' !E00_011_001
         call ErrWrite(ErrID,ErrCha)
         call parastop( 700 )

  999 continue

      write(ErrCha,*)'Error: reading dumpall file is stopped abnormally'
         ErrID = 'L:595/R:ovly15/F:ovly15.f' !E00_011_002
         call ErrWrite(ErrID,ErrCha)
         call parastop( 701 )
*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine dumpdat(ncol)
*                                                                      *
*       dump all data on file 09                                       *
*       last modified by K.Niita on 2011/07/06                         *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*       ncol = 1 : start of calculation                                *
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
*             13 : (n,x) reaction                                      *
*             14 : (n,n'x) reaction                                    *
*             15 : sequential transport only for tally                 *
*                                                                      *
*----------------------------------------------------------------------*
*        jcoll reaction type identifier                                *
*----------------------------------------------------------------------*
*                                                                      *
*        jcoll : =  0, nothing happen                                  *
*                =  1, Hydrogen collisions                             *
*                =  2, Particle Decays                                 *
*                =  3, Elastic collisions                              *
*                =  4, High Energy Nuclear collisions                  *
*                =  5, Heavy Ion reactions                             *
*                =  6, Neutron reactions by data                       *
*                =  7, Photon reactions by data                        *
*                =  8, Electron reactions by data                      *
*                =  9, P,d,a, and photo-nuclear reactions by data      *
*                = 11, Delta Ray production                            *
*                = 12, Muon atomic interaction                         *
*                = 13, Photon by EGS5                                  *
*                = 14, Electron by EGS5                                *
*                = 15, Photon photonuclear interaction                 *
*                = 16, Negative muon captured by nucleon               *
*                = 17, Muon photonuclear interaction                   *
*                = 18, Electron recoil by track strcuture mode         *
*                = 19, Muon pair production (photon -> mu+ mu-)        *
*                = 20, User defined interaction                        *
*                                                                      *
*----------------------------------------------------------------------*
*        kcoll reaction type identifier                                *
*----------------------------------------------------------------------*
*                                                                      *
*        kcoll : =  0, normal                                          *
*                =  1, high energy fission                             *
*                =  2, high energy absorption                          *
*                =  3, low energy n elastic                            *
*                =  4, low energy n non-elastic                        *
*                =  5, low energy n fission                            *
*                =  6, low energy n absorption                         *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*        nclsts   : total number of out going particles and nuclei     *
*                                                                      *
*        iclusts(nclsts)                                               *
*                                                                      *
*                i = 0, nucleus                                        *
*                  = 1, proton                                         *
*                  = 2, neutron                                        *
*                  = 3, pion                                           *
*                  = 4, photon                                         *
*                  = 5, kaon                                           *
*                  = 6, muon                                           *
*                  = 7, others                                         *
*                                                                      *
*        jclusts(i,nclsts)                                             *
*                                                                      *
*                i = 0, angular momentum                               *
*                  = 1, proton number                                  *
*                  = 2, neutron number                                 *
*                  = 3, ip, see below                                  *
*                  = 4, status of the particle 0: real, <0 : dead      *
*                  = 5, charge                                         *
*                  = 6, baryon number                                  *
*                  = 7, kf code                                        *
*                  = 8, isomer level (0: Ground, 1,2: 1st, 2nd isomer) *
*                                                                      *
*        qclusts(i,nclsts)                                             *
*                                                                      *
*                i = 0, impact parameter                               *
*                  = 1, px (GeV/c)                                     *
*                  = 2, py (GeV/c)                                     *
*                  = 3, pz (GeV/c)                                     *
*                  = 4, etot = sqrt( p**2 + rm**2 ) (GeV)              *
*                  = 5, rest mass (GeV)                                *
*                  = 6, excitation energy (MeV)                        *
*                  = 7, kinetic energy (MeV)                           *
*                  = 8, weight                                         *
*                  = 9, time                                           *
*                  = 10, x                                             *
*                  = 11, y                                             *
*                  = 12, z                                             *
*                                                                      *
************************************************************************
      use MMBANKMOD !FURUTA
      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param00.inc'
      include 'param.inc'

*-----------------------------------------------------------------------

      common /mpi00/ npe, me

*-----------------------------------------------------------------------

      common /jcomon/ nabov,nobch,nocas,nomax
!$OMP THREADPRIVATE(/jcomon/)
      common /rcomon/ rcasc

      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
!$OMP THREADPRIVATE(/icomon/)
      common /wtsave/ oldwt
!$OMP THREADPRIVATE(/wtsave/)
      common /tlgeom/ iblz1,iblz2
!$OMP THREADPRIVATE(/tlgeom/)
      common /tlglat/ ilev1,ilev2,ilat1(5,10),ilat2(5,10)
!$OMP THREADPRIVATE(/tlglat/)
      common /tlcost/ costha, uang(3), nsurf
!$OMP THREADPRIVATE(/tlcost/)

      common /elect/  uint(3), qs, qo, eint, delc, am, qsex,
     &                nq, ns, n1, noz, mtel
!$OMP THREADPRIVATE(/elect/)

      common /mathzn/ mathz, mathn, jcoll, kcoll
!$OMP THREADPRIVATE(/mathzn/)

      common /clustt/ nclsts, iclusts(nnn)
!$OMP THREADPRIVATE(/clustt/)
      common /clustw/ jclusts(0:8,nnn),  qclusts(0:12,nnn)
!$OMP THREADPRIVATE(/clustw/)
      common /cntcls/ jcount(3,nnn)
!$OMP THREADPRIVATE(/cntcls/)

      common /taliin/ rsouin, nzztin, nrgnin
      common /kmat1d/ idmn(0:kvlmax), idnm(kvmmax)

      common /dumpall/ idumpall

      common /sors0/ rcasc00,rsouin00,initsor
!$OMP THREADPRIVATE(/sors0/)

      character(80) dummy

      logical,allocatable :: exex(:)
      integer io2,io2ini,imax,ncol0,nobch0,nocas0
      real(8) rcascMAX,rcasc0,rsouin0

*-----------------------------------------------------------------------
*     write id is 9
*-----------------------------------------------------------------------

            io = 9
            io2ini = 800
            io2= io2ini+ipomp

*-----------------------------------------------------------------------
*     for control PE
*-----------------------------------------------------------------------

      if( me .eq. 0 .and. npe .gt. 1 ) then

               return

      end if

*-----------------------------------------------------------------------
*     binary or ascii
*     binary : idumpall = 1
*-----------------------------------------------------------------------

      if( idumpall .eq. 1 ) then

*-----------------------------------------------------------------------
*        ncol = 1, 2, 3
*-----------------------------------------------------------------------

         if( ncol .eq. 1 .or. ncol .eq. 2 .or. ncol .eq. 3 ) then

           if(ncol.eq.3)then
             allocate(exex(0:npomp-1))
             imax=0
             rcascMAX=0.0d0
             do i=0,npomp-1
               inquire(unit=io2ini+i, opened=exex(i))
               if(exex(i))then
                 write(io2ini+i) 20
                 rewind(io2ini+i)
                 read(io2ini+i) ncol0
                 read(io2ini+i) nocas0, nobch0, rcasc0, rsouin0
                 if(rcasc0.ge.rcascMAX)then
                   imax=i
                   rcascMAX=rcasc0
                 endif
               endif
             enddo
             do i=0,npomp-1
               if(exex(i).and.i.ne.imax)then
                 call cpcasdump(io,io2ini+i)
               endif
             enddo
             call cpcasdump(io,io2ini+imax)
             deallocate(exex)
           endif

           write(io) ncol
           return

         end if

*-----------------------------------------------------------------------
*        ncol = 4 - 15
*-----------------------------------------------------------------------

         if( ncol .ge. 4 ) then

*-----------------------------------------------------------------------

            if( ncol .eq. 4 ) then

               if(initsor.eq.0)then
                 initsor=1
               else
                 write(io2) 20
!$OMP CRITICAL (cpcasdump_crit)
                 call cpcasdump(io,io2)
!$OMP END CRITICAL (cpcasdump_crit)
               endif

               call opencasdump(io2,nocas,nobch)

               write(io2) ncol
               write(io2) nocas, nobch, rcasc00, rsouin00

            else

               write(io2) ncol

            end if

*-----------------------------------------------------------------------

               write(io2) no,idmn(mat),ityp,ktyp,jtyp,mtyp,rtyp, oldwt

            if( ityp .eq. 12 .or. ityp .eq. 13 ) then

               write(io2) qs

            end if

*-----------------------------------------------------------------------

               write(io2) iblz1,iblz2, ilev1,ilev2

            if( ilev1 .gt. 0 ) then

               write(io2) ( ( ilat1(i,j), i=1,5 ), j=1,ilev1 )

            end if

            if( ilev2 .gt. 0 ) then

               write(io2) ( ( ilat2(i,j), i=1,5 ), j=1,ilev2 )

            end if

              write(io2) costha, uang(1), uang(2), uang(3), nsurf

*-----------------------------------------------------------------------
               write(io2) name(ibknam+no,ipomp+1),
     &                       ( ncnt(ibknct+i,no,ipomp+1), i=1,3 )

               write(io2) wt(ibkwt+no,ipomp+1),
     &                   u(ibku+no,ipomp+1),v(ibkv+no,ipomp+1),
     &                   w(ibkw+no,ipomp+1)
               write(io2) e(ibke+no,ipomp+1), t(ibkt+no,ipomp+1),
     &                   x(ibkx+no,ipomp+1),y(ibky+no,ipomp+1),
     &                   z(ibkz+no,ipomp+1)
               write(io2) ec(ibkec+no,ipomp+1), tc(ibktc+no,ipomp+1),
     &                   xc(ibkxc+no,ipomp+1),yc(ibkyc+no,ipomp+1),
     &                   zc(ibkzc+no,ipomp+1)
               write(io2) spx(ibkspx+no,ipomp+1),spy(ibkspy+no,ipomp+1),
     &                    spz(ibkspz+no,ipomp+1)
               write(io2) nzst(ibkzst+no,ipomp+1)

*-----------------------------------------------------------------------

            if( ncol .eq. 13 .or. ncol .eq. 14 ) then

                  write(io2) nclsts

               if( nclsts .gt. 0 ) then

                     write(io2) mathz, mathn, jcoll, kcoll

                  do i = 1, nclsts

                     write(io2) iclusts(i)
                     write(io2) ( jclusts(j,i), j=0,8)
                     write(io2) ( qclusts(j,i), j=0,12)
                     write(io2) ( jcount(j,i),  j=1,3)

                  end do

               end if

            end if

         end if

*-----------------------------------------------------------------------
*     ascii : idumpall = -1
*-----------------------------------------------------------------------

      else if( idumpall .eq. -1 ) then

*-----------------------------------------------------------------------
*        ncol = 1, 2, 3
*-----------------------------------------------------------------------

         if( ncol .eq. 1 .or. ncol .eq. 2 .or. ncol .eq. 3 ) then

           if(ncol.eq.3)then
             allocate(exex(0:npomp-1))
             imax=0
             rcascMAX=0.0d0
             do i=0,npomp-1
               inquire(unit=io2ini+i, opened=exex(i))
               if(exex(i))then
                 rewind(io2ini+i)
                 read(io2ini+i,*) ncol0
                 read(io2ini+i,*) nocas0, nobch0, rcasc0, rsouin0
                 if(rcasc0.ge.rcascMAX)then
                   imax=i
                   rcascMAX=rcasc0
                 endif
               endif
             enddo
             do i=0,npomp-1
               if(exex(i).and.i.ne.imax)then
                 call cpcasdump(io,io2ini+i)
               endif
             enddo
             call cpcasdump(io,io2ini+imax)
             deallocate(exex)
           endif

           write(io,*) ncol
           return

         end if

*-----------------------------------------------------------------------
*        ncol = 4 - 15
*-----------------------------------------------------------------------

         if( ncol .ge. 4 ) then

*-----------------------------------------------------------------------

            if( ncol .eq. 4 ) then

               if(initsor.eq.0)then
                 initsor=1
               else
!$OMP CRITICAL (cpcasdump_crit)
                 call cpcasdump(io,io2)
!$OMP END CRITICAL (cpcasdump_crit)
               endif

               call opencasdump(io2,nocas,nobch)

               write(io2,*) ncol
               write(io2,*) nocas, nobch, rcasc00, rsouin00

            else

               write(io2,*) ncol

            end if

*-----------------------------------------------------------------------

               write(io2,*) no,idmn(mat),ityp,ktyp,jtyp,mtyp,rtyp, oldwt

            if( ityp .eq. 12 .or. ityp .eq. 13 ) then

               write(io2,*) qs

            end if

*-----------------------------------------------------------------------

               write(io2,*) iblz1,iblz2, ilev1,ilev2

            if( ilev1 .gt. 0 ) then

               write(io2,*) ( ( ilat1(i,j), i=1,5 ), j=1,ilev1 )

            end if

            if( ilev2 .gt. 0 ) then

               write(io2,*) ( ( ilat2(i,j), i=1,5 ), j=1,ilev2 )

            end if

              write(io2,*) costha, uang(1), uang(2), uang(3), nsurf

*-----------------------------------------------------------------------
               write(io2,*) name(ibknam+no,ipomp+1),
     &                          ( ncnt(ibknct+i,no,ipomp+1), i=1,3 )

               write(io2,*) wt(ibkwt+no,ipomp+1),
     &                     u(ibku+no,ipomp+1),v(ibkv+no,ipomp+1),
     &                     w(ibkw+no,ipomp+1)
               write(io2,*) e(ibke+no,ipomp+1), t(ibkt+no,ipomp+1),
     &                     x(ibkx+no,ipomp+1),y(ibky+no,ipomp+1),
     &                      z(ibkz+no,ipomp+1)
               write(io2,*) ec(ibkec+no,ipomp+1), tc(ibktc+no,ipomp+1),
     &                     xc(ibkxc+no,ipomp+1),yc(ibkyc+no,ipomp+1),
     &                     zc(ibkzc+no,ipomp+1)
               write(io2,*) spx(ibkspx+no,ipomp+1),
     &                      spy(ibkspy+no,ipomp+1),
     &                      spz(ibkspz+no,ipomp+1)
               write(io2,*) nzst(ibkzst+no,ipomp+1)

*-----------------------------------------------------------------------

            if( ncol .eq. 13 .or. ncol .eq. 14 ) then

                  write(io2,*) nclsts

               if( nclsts .gt. 0 ) then

                     write(io2,*) mathz, mathn, jcoll, kcoll

                  do i = 1, nclsts

                     write(io2,*) iclusts(i)
                     write(io2,*) ( jclusts(j,i), j=0,8)
                     write(io2,*) ( qclusts(j,i), j=0,12)
                     write(io2,*) ( jcount(j,i),  j=1,3)

                  end do

               end if

            end if

*-----------------------------------------------------------------------

         end if

      end if

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine opencasdump(io2,nocas,nobch)
*                                                                      *
************************************************************************
      implicit real*8 (a-h,o-z)

      include 'param.inc'

*-----------------------------------------------------------------------

      common /mpi00/ npe, me
      common /mpi03/ filhd, nfihd
      character filhd*100

*-----------------------------------------------------------------------
      common /paraj/  mstz(300), parz(300)
      common /dumpall/ idumpall
      common /paran/  icfn(100), ilfn(100), chfn(100)
      character       chfn*200

      character chme*5,chcas*8,chbch*8
      character filnm*100

      iord1 = aint(log10(real(nocas))) + 1
      write(chcas,'(i8.8)') nocas
      iord2 = aint(log10(real(nobch))) + 1
      write(chbch,'(i8.8)') nobch

         if( npe .le. 1 ) then

           filnm = chfn(15)(1:ilfn(15))//'_'//chbch(9-iord2:8)//'_'
     &          //chcas(9-iord1:8)

            if( idumpall .eq. 1 ) then

               open(io2, file = filnm,
     &                 form='unformatted',status = 'unknown' )

            end if

            if( idumpall .eq. -1 ) then

               open(io2, file = filnm,
     &                 form='formatted',status = 'unknown' )

            end if

         elseif( me .gt. 0)then

            if( mstz(61) .eq. 0 ) then

                  filnm = filhd(1:nfihd) // chfn(15)(1:ilfn(15))
     &             //'_'//chbch(9-iord2:8)//'_'//chcas(9-iord1:8)
                  iend=index(filnm,' ')-1

               if( idumpall .eq. 1 ) then

                  open(io2, file = filnm(1:iend),
     &                    form='unformatted',status = 'unknown' )

               else if( idumpall .eq. -1 ) then

                  open(io2, file = filnm(1:iend),
     &                    form='formatted',status = 'unknown' )

               end if

            else if( mstz(61) .eq. 1 ) then

                  iorder = aint(log10(real(npe))) + 1
                  if ( iorder .lt. 3) iorder = 3

                  write(chme,'(i5.5)') me

                  filnm = filhd(1:nfihd) // chfn(15)(1:ilfn(15))
     &                 //'_'//chbch(9-iord2:8)//'_'//chcas(9-iord1:8)
     &                 //'.' // chme(6-iorder:5)
                  iend=index(filnm,' ')-1

               if( idumpall .eq. 1 ) then

                  open(io2, file = filnm(1:iend),
     &                    form='unformatted',status = 'unknown' )

               else if( idumpall .eq. -1 ) then

                  open(io2, file = filnm(1:iend),
     &                    form='formatted',status = 'unknown' )

               end if

            else if( mstz(61) .eq. 2 ) then

                  filnm = chfn(15)(1:ilfn(15))
     &             //'_'//chbch(9-iord2:8)//'_'//chcas(9-iord1:8)
                  iend=index(filnm,' ')-1

               if( idumpall .eq. 1 ) then

                  open(io2, file = filnm(1:iend),
     &                    form='unformatted',status = 'unknown' )

               else if( idumpall .eq. -1 ) then

                  open(io2, file = filnm(1:iend),
     &                    form='formatted',status = 'unknown' )

               end if

            else if( mstz(61) .eq. 3 ) then

                  iorder = aint(log10(real(npe))) + 1
                  if ( iorder .lt. 3) iorder = 3

                  write(chme,'(i5.5)') me

                  filnm = chfn(15)(1:ilfn(15))
     &                 //'_'//chbch(9-iord2:8)//'_'//chcas(9-iord1:8)
     &                 //'.' // chme(6-iorder:5)
                  iend=index(filnm,' ')-1

               if( idumpall .eq. 1 ) then

                  open(io2, file = filnm(1:iend),
     &                    form='unformatted',status = 'unknown' )

               else if( idumpall .eq. -1 ) then

                  open(io2, file = filnm(1:iend),
     &                    form='formatted',status = 'unknown' )

               end if

            end if


         end if

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      subroutine cpcasdump(io,io2)
*                                                                      *
************************************************************************
      implicit none
      integer nnn,nomp
      include 'param00.inc'
      include 'err.inc'
      integer idumpall
      common /dumpall/ idumpall

      character(80) dummy
      integer i,j,io,io2,ios

      integer ncol,nocas,nobch
      real(8) rcasc,rsouin
      integer no,mat,ityp,ktyp,jtyp,mtyp
      real(8) rtyp,oldwt,qs
      integer iblz1,iblz2,ilev1,ilev2
      integer ilat1(5,10),ilat2(5,10)
      real(8) costha,uang(3)
      integer nsurf
      integer name,ncnt(3)
      real(8) wt,u,v,w,e,t,x,y,z,ec,tc,xc,yc,zc,spx,spy,spz
      integer nzst,nclsts,mathz,mathn,jcoll,kcoll
      integer, allocatable, save :: iclusts(:),jclusts(:,:),jcount(:,:)
      double precision, allocatable, save :: qclusts(:,:)
!$OMP THREADPRIVATE(iclusts,jclusts,jcount,qclusts)

       if((.not. allocated(iclusts))) then
        allocate(iclusts(13000),jclusts(0:8,13000),jcount(3,13000),
     &   qclusts(0:12,13000))
       else
        iclusts = 0
        jclusts = 0
        jcount  = 0
        qclusts = 0.d0
       endif

      rewind(io2)

*-----------------------------------------------------------------------
*     read dumpall (binary) file
*-----------------------------------------------------------------------

      if( idumpall .eq. 1 ) then

        do                      ! DO

          read(io2,iostat=ios,err=999) ncol
          if( ios .eq. -1 ) goto 900
          if( ncol .eq. 20 )then
            close(io2,status='DELETE')
            exit                ! EXIT DO
          endif
          write(io) ncol

          if(ncol.eq.4)then
            read(io2,iostat=ios,err=999) nocas, nobch, rcasc, rsouin
            if( ios .eq. -1 ) goto 900
            write(io) nocas, nobch, rcasc, rsouin
          endif

          read(io2,iostat=ios,err=999)
     &         no,mat,ityp,ktyp,jtyp,mtyp,rtyp, oldwt
          if( ios .eq. -1 ) goto 900
          write(io)no,mat,ityp,ktyp,jtyp,mtyp,rtyp, oldwt
          if( ityp .eq. 12 .or. ityp .eq. 13 ) then
            read(io2,iostat=ios,err=999) qs
            if( ios .eq. -1 ) goto 900
            write(io) qs
          end if

*-----------------------------------------------------------------------

          read(io2,iostat=ios,err=999) iblz1,iblz2, ilev1,ilev2
          if( ios .eq. -1 ) goto 900
          write(io) iblz1,iblz2, ilev1,ilev2
          if( ilev1 .gt. 0 ) then
            read(io2,iostat=ios,err=999)
     &           ( ( ilat1(i,j), i=1,5 ), j=1,ilev1 )
            if( ios .eq. -1 ) goto 900
            write(io) ( ( ilat1(i,j), i=1,5 ), j=1,ilev1 )
          end if
          if( ilev2 .gt. 0 ) then
            read(io2,iostat=ios,err=999)
     &           ( ( ilat2(i,j), i=1,5 ), j=1,ilev2 )
            if( ios .eq. -1 ) goto 900
            write(io) ( ( ilat2(i,j), i=1,5 ), j=1,ilev2 )
          end if

          read(io2,iostat=ios,err=999)
     &         costha, uang(1), uang(2), uang(3), nsurf
          if( ios .eq. -1 ) goto 900
          write(io) costha, uang(1), uang(2), uang(3), nsurf

*-----------------------------------------------------------------------
          read(io2,iostat=ios,err=999) name, ncnt(1:3)
          if( ios .eq. -1 ) goto 900
          write(io) name, ncnt(1:3)
          read(io2,iostat=ios,err=999) wt,u,v,w
          if( ios .eq. -1 ) goto 900
          write(io) wt,u,v,w
          read(io2,iostat=ios,err=999) e,t,x,y,z
          if( ios .eq. -1 ) goto 900
          write(io) e,t,x,y,z
          read(io2,iostat=ios,err=999) ec,tc,xc,yc,zc
          if( ios .eq. -1 ) goto 900
          write(io) ec,tc,xc,yc,zc
          read(io2,iostat=ios,err=999) spx,spy,spz
          if( ios .eq. -1 ) goto 900
          write(io) spx,spy,spz
          read(io2,iostat=ios,err=999) nzst
          if( ios .eq. -1 ) goto 900
          write(io) nzst
*-----------------------------------------------------------------------

          if( ncol .eq. 13 .or. ncol .eq. 14 ) then
            read(io2,iostat=ios,err=999) nclsts
            if( ios .eq. -1 ) goto 900
            write(io) nclsts
            if( nclsts .gt. 0 ) then
              read(io2,iostat=ios,err=999)
     &             mathz, mathn, jcoll, kcoll
              if( ios .eq. -1 ) goto 900
              write(io) mathz, mathn, jcoll, kcoll
              do i = 1, nclsts
                read(io2,iostat=ios,err=999) iclusts(i)
                if( ios .eq. -1 ) goto 900
                write(io) iclusts(i)
                read(io2,iostat=ios,err=999) (jclusts(j,i),j=0,8)
                if( ios .eq. -1 ) goto 900
                write(io) (jclusts(j,i),j=0,8)
                read(io2,iostat=ios,err=999) (qclusts(j,i),j=0,12)
                if( ios .eq. -1 ) goto 900
                write(io) (qclusts(j,i),j=0,12)
                read(io2,iostat=ios,err=999) ( jcount(j,i),j=1,3)
                if( ios .eq. -1 ) goto 900
                write(io) ( jcount(j,i),j=1,3)
              end do
            end if
          end if

          cycle                 ! NEXT DO

*-----------------------------------------------------------------------
*        ERROR MESSAGE
*-----------------------------------------------------------------------
 900      continue

          write(ErrCha,*) ' Error: cascade file is ended abnormally'
          ErrID = 'L:1377/R:cpcasdump/F:ovly15.f' !E00_012_001
          call ErrWrite(ErrID,ErrCha)
          call parastop( 700 )

 999      continue

          write(ErrCha,*) ' Error: cascade file is stopped abnormally'
          ErrID = 'L:1384/R:cpcasdump/F:ovly15.f' !E00_012_002
          call ErrWrite(ErrID,ErrCha)
          call parastop( 701 )
*-----------------------------------------------------------------------

        enddo

*-----------------------------------------------------------------------
*     read dumpall (ascii) file
*-----------------------------------------------------------------------

      elseif( idumpall .eq. -1 ) then

        do
          read(io2,'(a)',end=200)dummy
          write(io,'(a)')dummy
        enddo
 200    close(io2,status='DELETE')

      end if
*-----------------------------------------------------------------------

      return
      end

