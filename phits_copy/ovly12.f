************************************************************************
*                                                                      *
      subroutine ovly12
*                                                                      *
*       main control routine of PHITS                                  *
*       modified by K.Niita on 2007/09/26                              *
*                                                                      *
************************************************************************
!$    use omp_lib
      use mod_ompparallel !--- add NS 2020.04 del THREADPRIVATE

      use udm_Parameter, only : udm_initialize, iudmodel

C for USE_MOD_COUNTER
      use mod_counter, only: ALLOCATE_EVTS,DEALLOCATE_EVTS,INIT_EVTS
     &     ,rncnt,rnint,rnintr,rnpnt,rnpntr
     &     ,aevts,aevtr,bevts,bevtr
C for REDUCTION_COUNTER
!$   &     ,ALLOCATE_EVTS2,DEALLOCATE_EVTS2,INIT_EVTS2
!$   &     ,rncnt2,rnint2,rnintr2,rnpnt2,rnpntr2
!$   &     ,aevts2,aevtr2,bevts2,bevtr2

C for NONSHARED_TALLY
!$    use TALMOD0

C for NONCRITICAL_OVLY12
C for STRICT_RSOUIN
!$    use mod_rsouin

      use MMBANKMOD  !FURUTA
      use MEMBANKMOD !FURUTA
      use GGBANKMOD  !FURUTA
      use GGMBANKMOD !FURUTA
      use EVENTTALMOD !FURUTA
      use QMD_COOD2_MOD, only : cooddealloc
      use TETRAMOD, only : deallocate_welem,deallocate_volelems

      use fragdatamod, only: deallocate_FragData

      use t4dtrack_mod, only: it4dtrack,t4dtrack,init_t4dtrack
      implicit none

*-----------------------------------------------------------------------
      integer nnn,nomp

      include 'param00.inc'
      include 'param.inc'

      include 'err.inc'

*-----------------------------------------------------------------------
*     !FURUTA
*-----------------------------------------------------------------------
      integer irndmode,idmprijk
      common /irndm/ irndmode,idmprijk
      integer initbch
      data initbch/0/

      integer inobch,inocas,nobch00,nocas00,ipe
      integer istop
      integer inocas0 !FURUTA20130713

      integer i,k
      integer iot,iyer0,imon0,iday0,ihor0,imin0,isec0
      integer ncol,nbeta,icge,nsav,itmak
      integer mark,markp
      integer jrskip
      integer ichgf,ibryf,idcyc
      real(8) rmtyp,emint
      external ichgf,ibryf,idcyc,rmtyp
*-----------------------------------------------------------------------
      integer npe,me,nmbch0

      common /mpi00/ npe, me
      common /mpi02/ nmbch0

      integer nrmnbch0
      common /batchprocess/ nrmnbch0

*-----------------------------------------------------------------------
      integer no,mat,ityp,ktyp,jtyp,mtyp
      real(8) rtyp,ctyp
      integer nabov,nobch,nocas,nomax
      real(8) rcasc,oldwt
      integer iblz1,iblz2,ilev1,ilev2,ilat1,ilat2,idrg,idgr
      integer idmn,idnm
      real(8) egs,uus,vvs,wws,wts,tms
      integer nms,nct

      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
!$OMP THREADPRIVATE(/icomon/)
      common /jcomon/ nabov,nobch,nocas,nomax
!$OMP THREADPRIVATE(/jcomon/)
      common /rcomon/ rcasc
      common /wtsave/ oldwt
!$OMP THREADPRIVATE(/wtsave/)
      common /tlgeom/ iblz1,iblz2
!$OMP THREADPRIVATE(/tlgeom/)
      common /tlglat/ ilev1,ilev2,ilat1(5,10),ilat2(5,10)
!$OMP THREADPRIVATE(/tlglat/)
      common /regdc/  idrg(kvlmax), idgr(kvmmax)
      common /kmat1d/ idmn(0:kvlmax), idnm(kvmmax)
      common /pnsave/ egs, uus, vvs, wws, wts, tms, nms, nct(3)
!$OMP THREADPRIVATE(/pnsave/)

*-----------------------------------------------------------------------
      real(8) rnfb,rnfs,rngb,rngs,rnmult,rijk,ranj,rani,rans,
     &                ranb,rnrtc
      real(8),parameter:: randp=2d0**24,randq=2d0**(-24)
      integer nstrid,inif
      integer irskip
      integer irands,irandf,nseed,ncall
      real(8) randkk,srijk,rrijk

      common/irad/irands,irandf,nseed
!$OMP THREADPRIVATE(/irad/)
      common/ncall0/ncall
      integer nrandgen
      common /randn/ nrandgen ! S.H. xorshift (2020.2.6)
      integer*8 :: iranji64 ! S.H. xorshift (2020.2.6)
      common /randm4/ rnfb,rnfs,rngb,rngs,rnmult,ranj,rani,
     &                rnrtc,nstrid,inif, iranji64
      integer*8 :: iransb64 ! S.H. xorshift (2020.2.6)
      common /randtp/ rijk,rans,ranb, iransb64
!$OMP THREADPRIVATE(/randtp/)
      common /iradkk/ randkk,irskip
      common /randsv/ srijk,rrijk !FURUTA

*-----------------------------------------------------------------------
      real(8) esmax,esmin,emin,dnmax
      integer maxbch,maxcas
      real(8) andt
      integer jevap,npidk,iabsms
      real(8) tmax
      integer icntl,inucr,mathz,mathn,jcoll,kcoll
      integer ircss

      common /eparm/  esmax, esmin, emin(20)
      common /ndemax/ dnmax(20)
      common /cparm/  maxbch,maxcas
      common /bparm/  andt,jevap,npidk
      common /poabs/  iabsms
!$OMP THREADPRIVATE(/poabs/)
      common /tparm/  tmax(20)
      common /tcntl/  icntl, inucr
      common /mathzn/ mathz, mathn, jcoll, kcoll
!$OMP THREADPRIVATE(/mathzn/)

*-----------------------------------------------------------------------
      integer idpat,idoth,idono
      integer mttcn,mttc1,mttc2
      real(8) smttc

      common /otheid/ idpat(20), idoth(200), idono
      common /mttmc/  smttc(kvlmax), mttcn, mttc1(kvlmax), mttc2(kvlmax)

*-----------------------------------------------------------------------
      integer nclst,iclust,jclust,nclsts,iclusts,jclusts
      real(8) qclust,qclusts
      integer ieleh

      common /clustf/ nclst, iclust(nnn)
!$OMP THREADPRIVATE(/clustf/)
      common /clustg/ jclust(0:8,nnn), qclust(0:12,nnn)
!$OMP THREADPRIVATE(/clustg/)
      common /clustt/ nclsts, iclusts(nnn)
!$OMP THREADPRIVATE(/clustt/)
      common /clustw/ jclusts(0:8,nnn),  qclusts(0:12,nnn)
!$OMP THREADPRIVATE(/clustw/)

      common /electh/ ieleh

*-----------------------------------------------------------------------
      real(8) uint, qs, qo, eint, delc, am, qsex
      integer nqe, ns, n1, noz, mtel

      common /elect/  uint(3), qs, qo, eint, delc, am, qsex,
     &                nqe, ns, n1, noz, mtel
!$OMP THREADPRIVATE(/elect/)
*-----------------------------------------------------------------------

      logical pseudpr
      external pseudpr

*-----------------------------------------------------------------------
      integer mbnk,ibnk,jbnk,iibnk,jjbnk,imbnk,jmbnk,itbnk,jtbnk

      common /kcomon/ mbnk, ibnk, jbnk
!$OMP THREADPRIVATE(/kcomon/)
      common /kcomsi/ iibnk, jjbnk               !FURUTA
!$OMP THREADPRIVATE(/kcomsi/)
      common /kcomsm/ imbnk, jmbnk, itbnk, jtbnk !FURUTA

      data mbnk  / 0 /
      data ibnk  / 0 /
      data jbnk  / 0 /

*-----------------------------------------------------------------------
      integer ispfs,ispfn
      integer jstyp, istyp, inkf0, lstyp
      real(8) rspfn, rspfz
      real(8) ssx, ssy, ssz

      common /isorfs/ ispfs(isrc), rspfn, rspfz, ispfn
      common /isorpn/ ssx(isrc), ssy(isrc), ssz(isrc)
!$OMP THREADPRIVATE(/isorpn/)
      common /isorst/ jstyp(isrc), istyp(isrc), inkf0(isrc), lstyp(isrc)
!$OMP THREADPRIVATE(/isorst/)

*-----------------------------------------------------------------------
      real(8) dedxfd
      common /dedxfac/ dedxfd
!$OMP THREADPRIVATE(/dedxfac/)

*-----------------------------------------------------------------------

      integer iwwxyz
      common /wwin00/ iwwxyz
!$OMP THREADPRIVATE(/wwin00/)

      integer iresample ! T.Sato 2023/07/17

*-----------------------------------------------------------------------
*     definition of transport particles
*-----------------------------------------------------------------------

*                  1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
      data idpat / 1,1,1,1,1,1,1,1,1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0 /

      data idono / 55 /
      data idoth / 12,     -12,    14,   -14,
     &             221,   -221,   331,  -311, -2112, -2212,
     &             3122,  3222,  3212,  3112,  3322,  3312,  3334,
     &            -3122, -3222, -3212, -3112, -3322, -3312, -3334,
     &             130 ,   310, ! y.sakaki 2021/8
     &               16,  !      nu_tau  (0.0)
     &              -16,  ! anti-nu_tau  (0.0)
     &               15,  !      tau- (1777.0)
     &              -15,  !      tau+ (1777.0)
     &              113,  !      rho0 (768.50)
     &              213,  !      rho+ (766.90)
     &             -213,  !      rho- (766.90)
     &              223,  !      omega(781.94)
     &            20213,  !      a_1+ (1230.0)
     &           -20213,  !      a_1- (1230.0)
     &             -313,  ! anti-K*0  (896.10)
     &              313,  !      K*0  (896.10)
     &              323,  !      K*+  (891.60)
     &             -323,  !      K*-  (891.60)
     &              333,  !      phi  (1019.4)
     &              411,  !      D+   (1869.3)
     &             -411,  !      D-   (1869.3)
     &              421,  !      D0   (1864.5)
     &             -421,  ! anti-D0   (1864.5)
     &              431,  !      D_s+ (1968.5)
     &             -431,  !      D_s- (1968.5)
     &              511,  !      B0   (5279.2)
     &             -511,  ! anti-B0   (5279.2)
     &              521,  !      B+   (5278.9)
     &             -521,  !      B-   (5278.9)
     &              531,  !      B_s0 (5369.3)
     &             -531,  ! anti-B_s0 (5369.3)
     &              541,  !      B_c+ (6594.0)
     &             -541,  !      B_c- (6594.0)
     &             145*0 /

*-----------------------------------------------------------------------
      integer ipln

      common /ptname/ pname(20), ipln(20)
      character       pname*8

      data pname /'proton  ',
     &            'neutron ',
     &            'pion+   ',
     &            'pion0   ',
     &            'pion-   ',
     &            'muon+   ',
     &            'muon-   ',
     &            'kaon+   ',
     &            'kaon0   ',
     &            'kaon-   ',
     &            'other   ',
     &            'electron',
     &            'positron',
     &            'photon  ',
     &            'deuteron',
     &            'triton  ',
     &            '3he     ',
     &            'alpha   ',
     &            'nucleus ',
     &            'all     '/

      data      ipln / 6, 7, 9*5, 8, 8, 6, 8, 6, 3, 5, 7, 3 /

*-----------------------------------------------------------------------
      real(8) stime,cputm
      integer nlost,ilost,igerr,icger,ncger,nrecover
      integer mkc,ngp,ntyn,mpan,iet,ipsc,ixre,ixcos,nter

      common /cputim/ stime(40), cputm(40)
      common /cgerr/  nlost, ilost, igerr, icger, ncger, nrecover
      common /gm011/ mkc,ngp,ntyn,mpan,iet,ipsc,ixre,ixcos,nter
!$OMP THREADPRIVATE(/gm011/)

      real(8) rsouin
      integer nzztin,nrgnin
      common /taliin/ rsouin, nzztin, nrgnin

      integer initsor
      real(8) rcasc00,rsouin00
      common /sors0/ rcasc00,rsouin00,initsor
!$OMP THREADPRIVATE(/sors0/)
      data initsor /0/


      integer irndmode0,inobch0
      character(5) cipe

      real(8),allocatable::  saveRIJK(:)
      real(8),allocatable:: RIJKtot(:,:)

      integer icgg
      common /ccggg/  icgg

      integer nby,nlev,nar,nq,iaw,iay,nf,nx1
      common/arar/nby,nlev,nar,nq,iaw,iay,nf,nx1(3)
!$OMP THREADPRIVATE(/arar/)
      integer ierror
      data ierror/0/
      common/flag/ierror
!$OMP THREADPRIVATE(/flag/)
      integer kma,kfpd,klcr,knbd,kior,kriz,krcz,kmiz,kmcz,kkr1,
     &     kkr2,knsr,kvoll,nadd,ldata,ltma,lfpd,numr,irtru,
     &     numb,nir,kbiz,kbcz
      common/gomloc/kma,kfpd,klcr,knbd,kior,kriz,krcz,kmiz,kmcz,kkr1,
     &              kkr2,knsr,kvoll,nadd,ldata,ltma,lfpd,numr,irtru,
     &              numb,nir,kbiz,kbcz
!$OMP THREADPRIVATE(/gomloc/)
      integer  mus,muz,ll,ipret,iflow,iect,nlo,igx
      common /mgomv/  mus,muz,ll,ipret,iflow,iect,nlo,igx
!$OMP THREADPRIVATE(/mgomv/)
      integer intt,iot2,iout,iou2,idm,ioe
      common/tape/intt,iot2,iout,iou2,idm(4),ioe
!$OMP THREADPRIVATE(/tape/)
      integer ljy
      common/parem/ljy(44)
!$OMP THREADPRIVATE(/parem/)

      real(8) rijklst,rijkinit         !FURUTA20130226
      common /randm5/ rijklst,rijkinit !FURUTA20130226

*-----------------------------------------------------------------------

      integer ircode
      common /egs5cmn5/ircode
!$OMP THREADPRIVATE(/egs5cmn5/)

      common /egs5cmn6/pstep,dpmfp,gmfp
!$OMP THREADPRIVATE(/egs5cmn6/)
      real*8 pstep,dpmfp,gmfp

      integer iegsemi, iegsout
      common /egsemi/ iegsemi, iegsout

      integer ifgsq
      common /egs5cmn9/ifgsq
!$OMP THREADPRIVATE(/egs5cmn9/)

      common /paraj/ mstz(300), parz(300)
      real*8  parz
      integer mstz
      integer iii
      real*8 eeee

      integer idmpomp !FURUTA20150427
      common /idmpomp0/idmpomp !FURUTA20150427

! T.Sato 2016/05/28, for batch.out
      common /paran/ icfn(100), ilfn(100), chfn(100)
      character chfn*200
      integer icfn,ilfn

      integer ntscell, ktsc, mntsc, ntsc, mat1
      real*8 etsmin,etsmax,tsmax,ptsmin,ptsmax,ctsmin,ctsmax
      common /  tscmsg   / ktsc(kvlmax), mntsc, ntsc(kvlmax)
      common /  tscreg   / ntscell(kvlmax)
      common / etsminmax / etsmin, etsmax
      common / tsminmax  / tsmax
      common / ptsminmax / ptsmin, ptsmax
      common / ctsminmax / ctsmin, ctsmax
*-----------------------------------------------------------------------
      real(8) dmpmulti
      integer idmpmode,ibchjmp,idmpjmp
      common /stat2/ dmpmulti,idmpmode,ibchjmp,idmpjmp(2)
      integer ksoutnode,itetreg
      common /itetsor/ ksoutnode,itetreg(isrc)
      integer itettal
      common /itettal1/ itettal
*-----------------------------------------------------------------------
      integer icnt14
      common /ivmth/  icnt14
*-----------------------------------------------------------------------
      integer matold
      common /matoldcom/ matold
!$OMP THREADPRIVATE(/matoldcom/)
*-----------------------------------------------------------------------
      integer nfdopt5
      common /nfdopt5/ nfdopt5
*-----------------------------------------------------------------------
      integer iscinful, iswitch
      common /scincom/iscinful  ! use SCINFUL or not
      common /scinswt/iswitch
!$OMP THREADPRIVATE(/scinswt/)
*-----------------------------------------------------------------------

C for NONCRITICAL_OVLY12
!$    integer inocas_pre, idx
!$    real(8) rani_local, ranj_local, a_local, b_local
!$    integer*8 :: iranji64_local,it_local
C for NONSHARED_TALLY
      integer italsh
      common /talsh/ italsh

      logical lbnkskp ! flag to prioritize secondary particle tracking
      data lbnkskp /.false./
      integer maxbnk, maxbn2
      real(8) rtrckflp
      common /bnkmem/ maxbnk, maxbn2, rtrckflp

*-----------------------------------------------------------------------
      real(8) rn
*-----------------------------------------------------------------------
      integer mcol, mark_dum

*-----------------------------------------------------------------------
*
*               ityp :  description
*
*                 1  :  proton
*                 2  :  neutron
*                 3  :  pion (+)
*                 4  :  pion (0)
*                 5  :  pion (-)
*                 6  :  muon (+)
*                 7  :  muon (-)
*                 8  :  kaon (+)
*                 9  :  kaon (0)
*                10  :  kaon (-)
*
*                11  :  other
*
*                12  :  electron
*                13  :  positron
*
*                14  :  photon
*
*                15  :  deuteron
*                16  :  triton
*                17  :  3He
*                18  :  Alpha
*                19  :  Neucleus
*
*                20  :  All
*
*-----------------------------------------------------------------------
*           the other transport particles ( 24 )
*
*                 +-   12 : nu_e
*                 +-   14 : nu_mu
*                 +-  221 : eta
*                     331 : eta'
*                    -311 : k0bar
*                   -2112 : nbar
*                   -2212 : pbar
*                 +- 3122 : Lambda0
*                 +- 3222 : Sigma+
*                 +- 3212 : Sigma0
*                 +- 3112 : Sigma-
*                 +- 3322 : Xi0
*                 +- 3312 : Xi-
*                 +- 3334 : Omega-
*
*-----------------------------------------------------------------------
*
*       ncol = 1 : start of calculation
*              2 : end of calculation
*              3 : end of a batch
*              4 : source
*              5 : detection of geometry error
*              6 : recovery of geometry error
*              7 : termination by geometry error
*              8 : termination by weight cut-off
*              9 : termination by time cut-off
*             10 : geometry boundary crossing
*             11 : termination by energy cut-off
*             12 : termination by escape or leakage
*             13 : (n,x) reaction
*             14 : (n,n'x) reaction
*             15 : sequential transport only for tally
*             16 : surface cross for WW of xyz mesh
*

C for STRICT_RSOUIN
!$    call ALLOCATE_MODRSOUIN

*-----------------------------------------------------------------------
*        initialize the 'batch.out' file
*-----------------------------------------------------------------------

         if( me .eq. 0 ) then

               iot = 27
               open(iot,file=chfn(22),status='unknown')

               if( npe .le. 1 ) then !S.H.(2019.3.1)
               write(iot,'(i0,a)') maxbch,' <--- remaining batch number'
               nrmnbch0 = maxbch
               else
               write(iot,'(i0,a)') nmbch0,' <--- remaining batch number'
               nrmnbch0 = nmbch0
               end if

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

           if(irndmode.lt.0)then

            ErrCha = ''
            ErrID = 'L:549/R:ovly12/F:ovly12.f' !E00_010_001
            call ErrWrite(ErrID,ErrCha)

             write(*,'(2a,i5)')
     &            'ERROR: Random # generator works only with ',
     &            'single calc ',
     &            npe
             call parafin
             stop
           elseif(irndmode.gt.0)then
             if(irndmode.ne.max(npe-1,1))then
               call parafin
               stop
             endif
           endif

               if(icgg.ne.0)then
                 call ALLOCATE_GGBANK !FURUTA
                 call INIT_GGBANK     !FURUTA
               endif
               call ALLOCATE_GGMBANK  !FURUTA
               call INIT_GGMBANK      !FURUTA

               call ALLOCATE_EVTS !FURUTA20210623

               ncol  = 3

               call analyz(ncol,mark)

               ncol  = 2

               call analyz(ncol,mark)

               call DEALLOCATE_EVTS !FURUTA20210623

               if(icgg.ne.0) call DEALLOCATE_GGBANK  !FURUTA
               call DEALLOCATE_GGMBANK               !FURUTA

               call inittal_extstat      ! S.H. extstat 2024.4.1

               return

         end if


*-----------------------------------------------------------------------
*        initial event
*-----------------------------------------------------------------------

         jrskip=0               !FURUTA
         istop=1

*-----------------------------------------------------------------------
*        control of random number skip ( irskip ) MPI
*-----------------------------------------------------------------------

         if( irskip .ne. 0 ) then
           jrskip = iabs( irskip )
           if(irndmode.ne.0 .and. irskip .gt. 0)then
             rrijk=rijk
           endif
           do i = 1, jrskip
             call advijk
           end do
           if(irskip .gt. 0)then
             if(irndmode.ne.0)then
               srijk=rijk
               rijk=rrijk
              if ( nrandgen .eq. 0 ) then ! S.H. xorshift (2020.2.6)
               rani=aint(rrijk*randq)
               ranj=rrijk-rani*randp
              else
               iranji64 = transfer(rrijk,iranji64)
              end if
             endif
           else
             jrskip = 0
           end if
         end if

         if(irndmode.eq.0)then ! irndmode = 0
           do ipe=1,(me-1)        ! idle spinning for MPI !FURUTA20130722
             do inocas=1,maxcas   ! idle spinning for MPI !FURUTA20130722
               call advijk        ! idle spinning for MPI
             enddo                ! idle spinning for MPI !FURUTA20130722
           enddo                  ! idle spinning for MPI
         else
*-----------------------------------------------------------------------
* Random seed reading from files (rseed*.inp)
*-----------------------------------------------------------------------
* irndmode < 0 : Random # generator for |irndmode| processes -> rseed*.inp
* irndmode > 0 : MC calculation reading random # from rseed*.inp
*-----------------------------------------------------------------------
           irndmode0=abs(irndmode)
           if(irndmode.lt.0)then ! irndmode <0
             if(npe.le.1)then
               if(mod(maxbch,irndmode0).ne.0)then

                 ErrCha = ''
                 ErrID = 'L:648/R:ovly12/F:ovly12.f' !E00_010_002
                 call ErrWrite(ErrID,ErrCha)

                 write(*,'(2a,2i7)')
     &                'ERROR: Random # generator works only with ',
     &                'maxbch = irndmode * integer',maxbch,irndmode0
               else
                 allocate( RIJKtot(maxbch/irndmode0,irndmode0) )
                 do inobch=1,maxbch/irndmode0
                   do ipe=1,irndmode0
                     RIJKtot(inobch,ipe)=rijk
                     do inocas=1,maxcas
                       call advijk
                     enddo
                   enddo
                 enddo
                 do ipe=1,irndmode0
                   write(cipe,'(i5.5)')ipe
                   open(200,file='rseed'//cipe//'.inp',status='unknown')
                   do inobch=1,maxbch/irndmode0
                    write(200,'(i10,x,b64.64)')
     &                   inobch,RIJKtot(inobch,ipe)
                   enddo
                   close(200)
                 enddo
                 write(*,'(a,i6,a)')
     &                'Random #s are successfully generated!! ',
     &                irndmode0,' files'
               endif
             endif
             call parafin
             stop
           else ! irndmode > 0
             if(irndmode0.ne.max(npe-1,1))then

              ErrCha = ''
              ErrID = 'L:684/R:ovly12/F:ovly12.f' !E00_010_003
              call ErrWrite(ErrID,ErrCha)

               write(*,'(a,2i6)')
     &              'ERROR: Random # irndmode .ne. (npe-1)',
     &              irndmode,max(npe-1,1)
               call parafin
               stop
             endif
             allocate(saveRIJK(maxbch))
             write(cipe,'(i5.5)')max(me-1,0)+1
             open(200,file='rseed'//cipe//'.inp',status='old')
             do inobch=1,maxbch
              read(200,'(i10,x,b64.64)')
     &             inobch0,saveRIJK(inobch)
             enddo
             close(200)
           endif
*-----------------------------------------------------------------------
         endif

*-----------------------------------------------------------------------
*     ncol = 1 : start of calculation
*-----------------------------------------------------------------------

               ncol = 1

               call analyz(ncol,mark)

*-----------------------------------------------------------------------
*           control of irskip
*-----------------------------------------------------------------------

      if(jrskip.gt.0)then
        nobch00=jrskip/maxcas+1
        nocas00=jrskip-(nobch00-1)*maxcas+1
        rcasc=dble(maxcas)*dble(nobch00-1)+dble(nocas00-1)
        if(irndmode.ne.0) saveRIJK(nobch00)=srijk
      else
        nobch00=1
        nocas00=1
        rcasc=0.0d0
      endif

!$      if(italsh .eq. 0) then
!$      call COMPOSE_TALMOD
!$      call UPDATE_TAL00REF
!$      call DEALLOCATE_TAL0
!$      end if

C for REDUCTION_COUNTER
!$    call ALLOCATE_EVTS2
!$    call INIT_EVTS2

*-----------------------------------------------------------------------
*     ncol = 4 : start of new source
*        rcasc : current event number ( maxbch * maxcas )
*        nobch : current batch number ( maxbch )
*        nocas : current event number in a batch ( maxcas )
*-----------------------------------------------------------------------
      CALL  COPYIN_EGS5
      call init_sors !FURUTA20210324

      ipomp=0
      npomp=1
!$OMP PARALLEL

!$OMP& COPYIN(/irad/,/randtp/)
!$OMP& COPYIN(/isorpn/,/isorst/)
!$OMP& COPYIN(/gm011/)
!$OMP& COPYIN(/arar/,/flag/,/gomloc/,/mgomv/,/tape/,/parem/)

!$OMP& PRIVATE(inobch,inocas)
!$OMP& SHARED(nobch00,nocas00,jrskip,initbch,inocas0) !FURUTA20130713

!$OMP& SHARED(/mdasa/,/mdasb/)

!$OMP& PRIVATE(i,k)
!$OMP& SHARED(iot,iyer0,imon0,iday0,ihor0,imin0,isec0)
!$OMP& PRIVATE(ncol,nbeta,icge,nsav,itmak)
!$OMP& PRIVATE(mark,markp)
!$OMP& PRIVATE(emint)
!$OMP& SHARED(npe,me)
!$OMP& SHARED(istop)

!$OMP& SHARED(/rcomon/,/randsv/,/randm4/)
!$OMP& SHARED(/isorfs/)
!$OMP& SHARED(/regdc/,/kmat1d/)
!$OMP& SHARED(/eparm/,/ndemax/,/cparm/,/bparm/)
!$OMP& SHARED(/tparm/,/tcntl/)
!$OMP& SHARED(/otheid/,/mttmc/)
!$OMP& SHARED(/ptname/)
!$OMP& SHARED(/electh/)
!$OMP& SHARED(/taliin/)
!$OMP& SHARED(/itetsor/)

!$OMP& REDUCTION(+:stime,cputm)
!$OMP& REDUCTION(+:ilost,icger,ncger)
C for REDUCTION_COUNTER
!$OMP& REDUCTION (+:rncnt2,rnint2,rnintr2,rnpnt2,rnpntr2)
!$OMP& REDUCTION (+:aevts2,aevtr2,bevts2,bevtr2)

C for NONCRITICAL_OVLY12
!$OMP& PRIVATE(inocas_pre, idx)
!$OMP& PRIVATE(rani_local, ranj_local, a_local, b_local)
!$OMP& PRIVATE(iranji64_local, it_local)

!$      npomp=OMP_GET_NUM_THREADS()
!$      ipomp=OMP_GET_THREAD_NUM()
!$      write(*,'(''OpenMP PARALLEL PROCESS'',
!$   &     i4,''/'',i4,''  @ IP(MPI)='', i5)') ipOMP+1,npomp,me
!$OMP BARRIER
!------------------------------------------------------------------------
      if(iudmodel .ge. 1) then
        call udm_initialize
      endif

      call ALLOCATE_MMBANK   !FURUTA
      call ALLOCATE_MEMBANK  !FURUTA
      if(icgg.ne.0)then
        call ALLOCATE_GGBANK !FURUTA
      endif
      call ALLOCATE_GGMBANK  !FURUTA
      call INIT_GGMBANK      !FURUTA
C for NONSHARED_TALLY
!$      if(italsh .eq. 0) then
!$      call INIT_TALMOD0       ! duplicate from TALMOD0
!$      end if
      call ALLOCATE_EVENTTAL !FURUTA
      call ALLOCATE_EVTS !FURUTA20210119
      call INIT_EVTS     !FURUTA20210119

      call inittal_extstat      ! S.H. extstat 2024.4.1

      if(it4dtrack.gt.0)call init_t4dtrack

*-----------------------------------------------------------------------
*           end of a batch ( ncol = 3 ) or final ( ncol = 2 )
*-----------------------------------------------------------------------
       do inobch=1,maxbch
        if(istop.eq.0)cycle ! batch stop option
        if(inobch.lt.nobch00)cycle
        nobch=inobch

C for REDUCTION_COUNTER
!$    rncnt  = 0.0d0
!$    rnint  = 0.0d0
!$    rnintr = 0.0d0
!$    rnpnt  = 0.0d0
!$    rnpntr = 0.0d0
!$    aevts  = 0.0d0
!$    bevts  = 0.0d0
!$    aevtr  = 0.0d0
!$    bevtr  = 0.0d0

!$OMP MASTER
        if(irndmode.eq.0)then
          if(initbch.eq.0)then
            initbch=1
          else
            do ipe=1,(npe-2)        ! idle spinning for MPI !FURUTA20130722
              do inocas=1,maxcas    ! idle spinning for MPI !FURUTA20130722
                call advijk         ! idle spinning for MPI
              enddo                 ! idle spinning for MPI !FURUTA20130722
            enddo                   ! idle spinning for MPI
            if(npe.le.2)then
             if ( nrandgen .eq. 0 ) then ! S.H. xorshift (2020.2.6)
              rijk=rani*randp+ranj
             else
              rijk = transfer(iranji64,rijk)
             end if
            endif
          endif
          srijk=rijk
        else
          srijk = saveRIJK(inobch)
        endif
       if ( nrandgen .eq. 0 ) then ! S.H. xorshift (2020.2.6)
        rani=aint(srijk*randq)
        ranj=srijk-rani*randp
       else
        iranji64 = transfer(srijk,iranji64)
       end if
        inocas0=0 !FURUTA20150427 BUGFIX
        ibchjmp=min(ibchjmp,0) !FURUTA20150515
!$OMP END MASTER
!$OMP BARRIER
        if(idmpomp.ne.0)call openompfile !FURUTA20150427
        if(it4dtrack.gt.0)call t4dtrack(-1)

C for NONCRITICAL_OVLY12
!$      inocas_pre = 1
!$      if ( nrandgen .eq. 0 ) then
!$          rani_local = rani
!$          ranj_local = ranj
!$      else
!$          iranji64_local = iranji64
!$      end if
C for STRICT_RSOUIN
!$OMP MASTER
!$      rsouin_pre = rsouin
!$OMP END MASTER
!$OMP BARRIER
C --- end add NS 2020.04 ---

!$OMP DO SCHEDULE(guided)
        do inocas=1,maxcas
          if(inocas.lt.nocas00)cycle

          if(icgg.ne.0) call INIT_GGBANK !FURUTA

          mcol = 101
          mark_dum = 0
          call analyz(mcol,mark_dum)

*-----------------------------------------------------------------------
*           generate random number of next event
*-----------------------------------------------------------------------

C for NONCRITICAL_OVLY12
!$    if (.true.) then
!$        if ( nrandgen .eq. 0 ) then
!$          do idx=inocas_pre, inocas
!$            nocas   = idx
!$            ranb = rani_local
!$            rans = ranj_local
!$            a_local =  rnfs*ranj_local
!$            b_local = (rnfb*ranj_local
!$   &                   - aint(rnfb*ranj_local*randq)*randp)
!$   &                + (rnfs*rani_local
!$   &                   - aint(rnfs*rani_local*randq)*randp)
!$   &                + aint(a_local*randq)
!$            ranj_local = a_local - aint(a_local*randq)*randp
!$            rani_local = b_local - aint(b_local*randq)*randp
!$          enddo
!$          rijk = rani_local*randp+ranj_local
!$          inocas_pre = inocas + 1
!$        else
!$          do idx=inocas_pre, inocas
!$            nocas   = idx
!$            iransb64 = iranji64_local
!$            it_local = xor(iranji64_local,ishft(iranji64_local,13))
!$            it_local = xor(it_local,ishft(it_local,-7))
!$            iranji64_local = xor(it_local,ishft(it_local,17))
!$          enddo
!$          rijk = transfer(iranji64_local,rijk)
!$          inocas_pre = inocas + 1
!$        endif
!$    else
            inocas0=inocas0+1 !FURUTA20130713
            nocas=inocas0     !FURUTA20130713
            if ( nrandgen .eq. 0 ) then ! S.H. xorshift (2020.2.6)
               ranb = rani
               rans = ranj
            else
               iransb64 = iranji64
            end if
            if(idmprijk.eq.1)call dmprijk
            call advijk
C for NONCRITICAL_OVLY12
!$    endif

*-----------------------------------------------------------------------
*        start of new source
*-----------------------------------------------------------------------
C for NONCRITICAL_OVLY12
!$        if (.true.) then
!$          rcasc=dble(maxcas)*dble(inobch-1)+dble(inocas)
!$        else
               rcasc=rcasc+1.0
C for NONCRITICAL_OVLY12
!$        endif
               jcoll = 0
               kcoll = 0
  90           call sors(iresample) ! T.Sato 2023/07/17, occasionally resampling is required
               if(iresample.eq.1) goto 90
               rcasc00=rcasc
               rsouin00=rsouin

      call ionts_setup


*-----------------------------------------------------------------------
*        EGS5 for the first electron step
*-----------------------------------------------------------------------

               ircode = -1
               dpmfp = 0 ! reset dpmfp, T.Sato 2019/7/6

*-----------------------------------------------------------------------
*           analyz ( ncol = 4 )
*-----------------------------------------------------------------------

               ncol = 4
               mark = 1

               call analyz(ncol,mark)

*-----------------------------------------------------------------------
cKN 2024/03/26

               call wwindw(ncol)

cKN 2024/03/26
*-----------------------------------------------------------------------
*  next cascading particle
*-----------------------------------------------------------------------

               no = no - 1

  600 continue

               no = no + 1

*-----------------------------------------------------------------------

               if( no .gt. nomax ) then

                  if( mbnk .ne. 0 ) call update(1)

                  if( mbnk .eq. 0  .or. no .eq. 0) then

                     if( lbnkskp ) then
                       lbnkskp = .false.
                       call update(3)
                       goto 600
                     endif

                     cycle

                  end if

               end if

               if( icntl .eq. 6 ) cycle

*-----------------------------------------------------------------------
               if( icnt14 .ne. 0 ) cycle

*-----------------------------------------------------------------------
*           electron and positron above dnmax(12)
*-----------------------------------------------------------------------

            if( ieleh .ne. 0 .and.
     &        ( nty(ibknty+no,ipomp+1) .eq. 12 .or.
     &          nty(ibknty+no,ipomp+1) .eq. 13 ) .and.
     &          e(ibke+no,ipomp+1) .gt. dnmax(12) ) then

               wt(ibkwt+no,ipomp+1) = wt(ibkwt+no,ipomp+1) *
     &                               e(ibke+no,ipomp+1) / dnmax(12)
               e(ibke+no,ipomp+1)   = dnmax(12) - 1.e-6

            end if

*-----------------------------------------------------------------------
            if(mstz(85).eq.99)then
               write(93,*)'------------ phitsdump -------------'
               do iii = no,nomax
                  if(nty(ibknty+iii,ipomp+1).eq.14)then
                     eeee = e(ibke+iii,ipomp+1)
                  else
                     eeee = e(ibke+iii,ipomp+1)+0.510998902
                  endif

                  write(93,'(i3,7f10.5)')nty(ibknty+iii,ipomp+1),
     $                 eeee, x(ibkx+iii,ipomp+1),y(ibky+iii,ipomp+1),
     &                       z(ibkz+iii,ipomp+1)
     $                      ,u(ibku+iii,ipomp+1),v(ibkv+iii,ipomp+1),
     &                       w(ibkw+iii,ipomp+1)
               enddo
               write(93,*)'------------------------------------'
            endif
               if(abs(nkf(ibknkf+no,ipomp+1)) .eq. 311) then
                 nty(ibknty+no,ipomp+1) = 11
                 if(rn(0) .gt. 0.5d0) then
                   nkf(ibknkf+no,ipomp+1) = 130
                 else
                   nkf(ibknkf+no,ipomp+1) = 310
                 endif
               endif
               ityp = nty(ibknty+no,ipomp+1)
               ktyp = nkf(ibknkf+no,ipomp+1)
               jtyp = ichgf(ityp,ktyp)
               mtyp = ibryf(ityp,ktyp)
               rtyp = rmtyp(ityp,ktyp)
               ctyp = nzst(ibknkf+no,ipomp+1) ! 2022/6/8 Ogawa. Store chage state as ctyp

               mark  = 1
               markp = 0
               nsav  = 0
               icge  = 0
               jcoll = 0

               ns = 0
               qo = 0.0d0

               nabov = 0

               xc(ibkxc+no,ipomp+1) = x(ibkx+no,ipomp+1)
               yc(ibkyc+no,ipomp+1) = y(ibky+no,ipomp+1)
               zc(ibkzc+no,ipomp+1) = z(ibkz+no,ipomp+1)
               ec(ibkec+no,ipomp+1) = e(ibke+no,ipomp+1)
               tc(ibktc+no,ipomp+1) = t(ibkt+no,ipomp+1)

               mat = nmed(ibknmd+no,ipomp+1)
               matold = mat !FURUTA20190904

               oldwt = wt(ibkwt+no,ipomp+1)

*-----------------------------------------------------------------------
*        reset EGS5 electron parameters
*-----------------------------------------------------------------------

            if( ( ityp .eq. 12 .or. ityp .eq. 13 ) .and.
     &            iegsemi .ne. 0 ) call egs5esteps0

*-----------------------------------------------------------------------
*           mat time change
*-----------------------------------------------------------------------

            if( mttcn .gt. 0 .and. mat .gt. 0 ) then

               do k = 1, mttcn

                  if( mttc1(k) .eq. idmn(mat) .and.
     &                smttc(k) .lt. abs(t(ibkt+no,ipomp+1)) ) then

                     if( mttc2(k) .gt. 0 ) then

                        nmed(ibknmd+no,ipomp+1) = idnm(mttc2(k))

                     else

                        nmed(ibknmd+no,ipomp+1) = mttc2(k)

                     end if

                        mat = nmed(ibknmd+no,ipomp+1)
                        matold = mat !FURUTA20190904

                  end if

               end do

            end if

*-----------------------------------------------------------------------

            if( ityp .lt. 15 ) then

                  emint = emin(ityp)

*-----------------------------------------------------------------------
               if( mntsc .gt. 0 .and.
     &             ntscell( idgr(iblz(ibkblz+no,ipomp+1)) ) .ne. 0 .and.
     &           ( ityp .eq. 12 .or. ityp .eq. 13 ) .and.
     &           ( e(ibke+no,ipomp+1) .le. etsmax .and.
     &             e(ibke+no,ipomp+1) .ge. etsmin ) )  emint = etsmin

*-----------------------------------------------------------------------

            else

               emint = emin(ityp) * dble( mtyp )

            end if

*-----------------------------------------------------------------------
*        cut off; low energy, time out or zero weight
*        stopped decay particle -> 400, iabsms = 2
*-----------------------------------------------------------------------

            if( mat .eq. 0 .and. ktyp .eq. -11) then ! positrons do not "decay" in vacuum

                  continue

            elseif( mat .eq. 0 .and. ktyp .eq. 11 .and.
     &              e(ibke+no,ipomp+1) .gt. etsmin) then ! do not kill for potential track structure calculation

                  continue

            elseif( e(ibke+no,ipomp+1) .le. emint ) then

               if( idcyc(ktyp) .eq. 1 .and. ktyp .ne. 2112) then ! 2016/9/21 Ogawa, Inhibit neutron decay when energy cut-off

                  iabsms = 2
                  goto 400

               end if

                  ncol = 11
                  goto 500

            end if

            if( wt(ibkwt+no,ipomp+1) .le. 0.0d0 ) then

                  ncol = 8
                  goto 500

            end if

            if( abs(t(ibkt+no,ipomp+1)) .gt. tmax(ityp) ) then

                  ncol = 9
                  goto 500

            end if

*-----------------------------------------------------------------------
*     continue the flight of particle
*-----------------------------------------------------------------------

  800 continue

      if( ncol .eq. 10 .and. mat .eq. matold ) then
         ifgsq = 1
      else
         ifgsq = 0
      end if

               mat   = nmed(ibknmd+no,ipomp+1)
               matold = mat !FURUTA20190904
               oldwt = wt(ibkwt+no,ipomp+1)

*-----------------------------------------------------------------------
*           reset EGS5 edep
*-----------------------------------------------------------------------
            if( iegsemi .ne. 0 ) call egs5init2

*-----------------------------------------------------------------------
*           mat time change
*-----------------------------------------------------------------------

            if( mttcn .gt. 0 .and. mat .gt. 0 ) then

               do k = 1, mttcn

                  if( mttc1(k) .eq. idmn(mat) .and.
     &                smttc(k) .lt. abs(t(ibkt+no,ipomp+1)) ) then

                     if( mttc2(k) .gt. 0 ) then

                        nmed(ibknmd+no,ipomp+1) = idnm(mttc2(k))

                     else

                        nmed(ibknmd+no,ipomp+1) = mttc2(k)

                     end if

                        mat = nmed(ibknmd+no,ipomp+1)
                        matold = mat !FURUTA20190904

                  end if

               end do

            end if

*-----------------------------------------------------------------------
*        space and time transport of the particle
*-----------------------------------------------------------------------
            dedxfd = 1.d0
            if( ityp .ne. 12 .and. ityp .ne. 13 ) ns = 0


            call partrs(mark,markp,nbeta,icge,itmak,emint) !(emint:Takeshi Kai)


*-----------------------------------------------------------------------
*        geometry error
*-----------------------------------------------------------------------
*            icge =  0 : no error
*                 >  0 : CG error try again : ncol = 5
*                 = -3 : recovered error    : ncol = 6
*                 = -2 : unrecovered error  : ncol = 7
*                 = -1 : lost particle      : ncol = 7
*-----------------------------------------------------------------------

         if( icge .ne. 0 ) then

               if( icge .eq. -1 .or. icge .eq. -2 ) then

                     ncol = 7

               else if( icge .eq. -3 ) then

                     icge = 0
                     ncol = 6

               else if( icge .ge. 1 ) then

                  if( ncol .eq. 6 ) then

                     ncol = 7

                  else

                     ncol = 5

                  end if

               end if

         else

*-----------------------------------------------------------------------
*        outgoing to the void region (mark = -1 ) : ncol = 12
*-----------------------------------------------------------------------

            if( mark .eq. -1 ) then

                  ncol = 12

*-----------------------------------------------------------------------
*        pass the forward boundary ( mark = 0 ) : ncol = 10
*        pass the reflect boundary ( mark = 2 ) : ncol = 10
*-----------------------------------------------------------------------

            else if( mark .eq. 0 .or. mark .eq. 2 ) then

                  ncol = 10
                  if( iwwxyz .eq. 2 ) ncol = 16

*-----------------------------------------------------------------------
*        time cutoff
*           mark = 1 and itmak = 1 : time cutoff : ncol = 9
*-----------------------------------------------------------------------

            else if( mark .eq. 1 .and. itmak .eq. 1 ) then

                  ncol = 9

*-----------------------------------------------------------------------
*        collision or stopped ( mark = 1 )
*-----------------------------------------------------------------------
*              nbeta < 3 ( nuclear collisions ) -> 400
*                 for proton : ( pseudo collision ) -> 300
*              nbeta = 3 ( stopped by energy loss ) : ncol = 11
*                 ec(ibkec+no) > emin ; ( pseudo energy dump ) -> 300
*                 for npidk = 0, charge<0, decay p -> 400, iabsms = 1
*                 for npidk = 1, decay particle    -> 400, iabsms = 2
*                 iabsms = 1, forced absorption
*                 iabsms = 2, forced decay
*-----------------------------------------------------------------------

            else if( mark .eq. 1 ) then

                     iabsms = 0

               if( nbeta .lt. 3 ) then

                  if( nfcs(ibknfc+no,ipomp+1) .ne. 2 ) then

                     if( pseudpr(ityp) ) goto 300

cfrtati 2023/12/22
                  else
                     call getrealflt(ityp)
                  end if

cKN 2018/01/08
                  if( iwwxyz .eq. 1 ) then

                     ncol = 16
                     goto 500

                  else

                     goto 400

                  end if

               else if( nbeta .eq. 3 ) then

                     ncol = 11

                  if( ec(ibkec+no,ipomp+1) .gt. emint ) then

                     goto 300

                  else if( npidk .eq. 0 .and.
     &                     idcyc(ktyp) .eq. 1 .and.
     &                     jtyp .lt. 0 .and.
     &                     ktyp .ne. -11 .and.
     &                     mat .gt. 0 ) then

                     iabsms = 1

                     goto 400

                  else if( idcyc(ktyp) .eq. 1 .and. ktyp .ne. 2112) then ! 2016/9/21 Ogawa, Inhibit neutron decay when energy cut-off

                     iabsms = 2

                     goto 400

                  end if

               end if

*-----------------------------------------------------------------------

            end if

         end if

*-----------------------------------------------------------------------
*     analysis and tally
*-----------------------------------------------------------------------

  500 continue

*-----------------------------------------------------------------------
*           analysis
*-----------------------------------------------------------------------

               call analyz(ncol,mark)

*-----------------------------------------------------------------------
*           ncol = 14 : sequential transport of scattered particle
*-----------------------------------------------------------------------


               if( ncol .eq. 14 ) then

                  ec(ibkec+no,ipomp+1)    = egs
                  tc(ibktc+no,ipomp+1)    = tms

                  u(ibku+no,ipomp+1)      = uus
                  v(ibkv+no,ipomp+1)      = vvs
                  w(ibkw+no,ipomp+1)      = wws
                  wt(ibkwt+no,ipomp+1)    = wts
                  name(ibknam+no,ipomp+1) = nms

                  ncnt(ibknct+1,no,ipomp+1) = nct(1)
                  ncnt(ibknct+2,no,ipomp+1) = nct(2)
                  ncnt(ibknct+3,no,ipomp+1) = nct(3)

                  call gomupp(mark,markp,uus,vvs,wws)

               end if

*-----------------------------------------------------------------------
*           icntl = 5 : no reaction, no ionization mode
*-----------------------------------------------------------------------

               if( icntl .eq.  5 ) then

                  x(ibkx+no,ipomp+1) = xc(ibkxc+no,ipomp+1)
                  y(ibky+no,ipomp+1) = yc(ibkyc+no,ipomp+1)
                  z(ibkz+no,ipomp+1) = zc(ibkzc+no,ipomp+1)

                  e(ibke+no,ipomp+1) = ec(ibkec+no,ipomp+1)
                  t(ibkt+no,ipomp+1) = tc(ibktc+no,ipomp+1)

               end if

*-----------------------------------------------------------------------
            if( icntl .ne.  5 .and.
     &          icntl .ne. 14 .and. icntl .ne. 15 ) then

*-----------------------------------------------------------------------
*           cell importance function for ncol = 10
*-----------------------------------------------------------------------

                  call celimp(ncol)

*-----------------------------------------------------------------------
*           forced collisions        for ncol = 10, 13, 14
*-----------------------------------------------------------------------

                  call fclsta(ncol,mark,markp)

*-----------------------------------------------------------------------
*           weight window            for ncol = 10, 13, 14, 16
*-----------------------------------------------------------------------

                  call wwindw(ncol)

*-----------------------------------------------------------------------
*           weight cutoff            for ncol = 10, 13, 14
*-----------------------------------------------------------------------

                  call wtcutof(ncol)

*-----------------------------------------------------------------------
*           splitting                for ncol = 10
*-----------------------------------------------------------------------

                  call splitng(ncol)

*-----------------------------------------------------------------------
*           update the position and energy
*-----------------------------------------------------------------------


                  x(ibkx+no,ipomp+1) = xc(ibkxc+no,ipomp+1)
                  y(ibky+no,ipomp+1) = yc(ibkyc+no,ipomp+1)
                  z(ibkz+no,ipomp+1) = zc(ibkzc+no,ipomp+1)

                  e(ibke+no,ipomp+1) = ec(ibkec+no,ipomp+1)
                  t(ibkt+no,ipomp+1) = tc(ibktc+no,ipomp+1)

*-----------------------------------------------------------------------

*-----------------------------------------------------------------------
*           tally for dead particles
*-----------------------------------------------------------------------

               if( wt(ibkwt+no,ipomp+1) .le. 0.0d0 ) then

                  ncol = 8
                  call analyz(ncol,mark)

               end if

*-----------------------------------------------------------------------
*           put generated particles into bank, ic = 2 : normal update
*-----------------------------------------------------------------------

                  call update(2)

*-----------------------------------------------------------------------
            end if

*-----------------------------------------------------------------------
*           next flight ( -> 800 ) or next particle ( -> 600 )
*-----------------------------------------------------------------------

            if( ncol .eq.  5 .or. ncol .eq.  6 .or.
     &          ncol .eq. 10 .or. ncol .eq. 14 .or.
     &          ncol .eq. 16 ) then

                  if(nomax .gt. int(dble(maxbnk)*rtrckflp) .and.
     &               nomax .gt. 100 * no      .and. .not. lbnkskp) then ! 2021/12/13 skip the current particle and track secondary particles
                      lbnkskp = .true.
                      call update(4)
                      goto 600
                  endif

                  goto 800

            else

                  goto 600

            end if

*-----------------------------------------------------------------------
*     pseudo collisions or pseudo energy dump, go back to 800
*-----------------------------------------------------------------------

  300 continue

            dedxfd = 1.d0

*-----------------------------------------------------------------------
*           call tally at ncol = 15
*-----------------------------------------------------------------------

                  ncol = 15

                  call analyz(ncol,mark)

*-----------------------------------------------------------------------

                  e(ibke+no,ipomp+1) = ec(ibkec+no,ipomp+1)
                  t(ibkt+no,ipomp+1) = tc(ibktc+no,ipomp+1)
                  x(ibkx+no,ipomp+1) = xc(ibkxc+no,ipomp+1)
                  y(ibky+no,ipomp+1) = yc(ibkyc+no,ipomp+1)
                  z(ibkz+no,ipomp+1) = zc(ibkzc+no,ipomp+1)

                  icge  = 0

                  goto 800

*-----------------------------------------------------------------------
*     nuclear reaction, decay and evaporation
*-----------------------------------------------------------------------
*           nclst =< 0 : pseudo collision            -> 300
*             fail of pi-, k- absorption : ncol = 11 -> 500
*           nclst > 0 : collision happened,   dataup -> 500
*             projectile disappeared    :  ncol = 13
*             projectile appeared       :  ncol = 14
*-----------------------------------------------------------------------

  400 continue

*-----------------------------------------------------------------------
*           iabsms = 0
*-----------------------------------------------------------------------

            if( iabsms .eq. 0 ) then

                     call nreac

               if( ( jcoll .ge. 1 .and. jcoll .le. 5 ) .or.
     &               jcoll .eq. 11 ) then

                     call nevap(0)

                  if( nclst .le. 0 ) goto 300

               end if

                  if( nclsts .lt. 0 ) goto 300

*-----------------------------------------------------------------------
*           repeated collsions
*-----------------------------------------------------------------------

                     call rcldol(ircss,ncol,mark)

                     if( ircss .ne. 0 ) goto 600

*-----------------------------------------------------------------------
*           iabsms = 1, 2
*-----------------------------------------------------------------------

            else

                     call nreac

               if( ( jcoll .ge. 1 .and. jcoll .le. 5 ) .or.
     &               jcoll .eq. 11 ) then

                     call nevap(0)

                  if( nclst .le. 0 .and. iabsms .eq. 1 ) then

                     iabsms = 2
                     goto 400

                  else if( nclst .le. 0 .and. iabsms .eq. 2 ) then

                     ncol = 11
                     goto 500

                  end if

                  if( nclst .le. 0 ) goto 300

               end if

                  if( nclsts .lt. 0 ) goto 300

            end if

                     call dataup(ncol,14)

                     goto 500

*-----------------------------------------------------------------------
        enddo ! inocas=1,maxcas
!$OMP END DO
C for NONCRITICAL_OVLY12
!$      if ((inocas_pre-1) .eq. maxcas) then
!$      if ( nrandgen .eq. 0 ) then
!$        rani = rani_local
!$        ranj = ranj_local
!$      else
!$        iranji64 = iranji64_local
!$      end if
!$        rcasc=dble(maxcas)*dble(inobch-1)+dble(maxcas)
!$      endif
!$OMP MASTER
!$      rsouin = rsouin_pre
!$      do idx=1, maxcas
!$        rsouin = rsouin + rsouin_arr(idx)
!$      enddo
!$OMP END MASTER

        mcol = 101
        mark_dum = 0
        call analyz(mcol,mark_dum)

!$OMP BARRIER
        call EVENTanalyz        !FURUTA
!$OMP CRITICAL (ompfile)
        if(idmpomp.ne.0)call cpompfile !FURUTA20150427
        if(it4dtrack.gt.0)call t4dtrack(-3)
!$OMP END CRITICAL (ompfile)
C for NONSHARED_TALLY
!$        if(italsh .eq. 0) then
!$OMP CRITICAL (TY1)
!$        call COMPOSE_TALMOD  ! update global
!$OMP END CRITICAL (TY1)
!$        end if

!$OMP BARRIER

C for NONSHARED_TALLY
!$        if(italsh .eq. 0) then
!$        call SYNC_TALMOD     ! copyback global to local
!$        end if

        ncol=3
        initsor=0
!$OMP MASTER
C for NONSHARED_TALLY
!$        if(italsh .eq. 0) then
!$        call UPDATE_TAL00REF
!$        end if
        if(irndmode.eq.0)then                    !FURUTA20130226
         if ( nrandgen .eq. 0 ) then ! S.H. xorshift (2020.2.6)
          rijklst=rani*randp+ranj                !FURUTA20130226
         else
          rijklst=transfer(iranji64,rijklst)
         end if
        else                                     !FURUTA20130226
          rijklst=saveRIJK(min(inobch+1,maxbch)) !FURUTA20130226
        endif                                    !FURUTA20130226
        if(icgg.ne.0) call INIT_GGBANK !FURUTA
        call analyz(ncol,mark)
        if(nocas00.ne.1)nocas00=1
        if(ncol.eq.-10)istop=0
        if(ibchjmp.eq.0.and.idmpmode.eq.1)then

         ErrCha = ''
         ErrID = 'L:1705/R:ovly12/F:ovly12.f' !E00_010_004
         call ErrWrite(ErrID,ErrCha)

         write(*,'(a)')'** ERROR: NOCAS of dump source exceeds MAXCAS'
         write(*,'(3a)')'MAXCAS should be equal to the number same as',
     &        ' the first process producing the dump source',
     &        ' for idmpmode=1'
         write(*,'(2(a,i5))')' ibchjmp = ',ibchjmp,' at NOBCH=',nobch
         call parastop( 712 )
        endif
c------------------
!$OMP END MASTER

C for REDUCTION_COUNTER
!$    rncnt2(:)  = rncnt(:)  + rncnt2(:)
!$    rnint2(:)  = rnint(:)  + rnint2(:)
!$    rnintr2(:) = rnintr(:) + rnintr2(:)
!$    rnpnt2(:)  = rnpnt(:)  + rnpnt2(:)
!$    rnpntr2(:) = rnpntr(:) + rnpntr2(:)
!$    aevts2 = aevts2 + aevts
!$    aevtr2 = aevtr2 + aevtr
!$    bevts2 = bevts2 + bevts
!$    bevtr2 = bevtr2 + bevtr

!$OMP BARRIER
      enddo ! inobch=1,maxbch

      call DEALLOCATE_MMBANK    !FURUTA
      call DEALLOCATE_MEMBANK   !FURUTA
      call DEALLOCATE_EVENTTAL  !FURUTA
      call cooddealloc          !Ogawa
!$      write(*,'(''OpenMP FINALIZE'',
!$   &     i4,''/'',i4,''  @ IP(MPI)='', i5)') ipOMP+1,npomp,me

C for REDUCTION_COUNTER
!$OMP END PARALLEL

C  restore rncnt...
!$    rncnt(:)  = rncnt2(:)
!$    rnint(:)  = rnint2(:)
!$    rnintr(:) = rnintr2(:)
!$    rnpnt(:)  = rnpnt2(:)
!$    rnpntr(:) = rnpntr2(:)

!$    aevts = aevts2
!$    aevtr = aevtr2
!$    bevts = bevts2
!$    bevtr = bevtr2


      if ( nrandgen .eq. 0 ) then ! S.H. xorshift (2020.2.6)
         rrijk=rani*randp+ranj
      else
         rrijk=transfer(iranji64,rrijk)
      end if
      if(irndmode.gt.0)then
        if(istop.eq.0)then
          maxbch=nobch
          rrijk=saveRIJK(nobch+1)
        endif
        deallocate( saveRIJK )
      endif
      rijklst=rrijk

      ncol=2
      call analyz(ncol,mark)

      call DEALLOCATE_EVTS      !FURUTA20210119
!$    call DEALLOCATE_EVTS2     !FURUTA20210119

      if(ksoutnode.ne.0)call deallocate_welem
      if(itettal.gt.0)call deallocate_volelems

      if ( nfdopt5 .gt. 0 ) call deallocate_FragData


      if(icgg.ne.0) call DEALLOCATE_GGBANK !FURUTA
      call DEALLOCATE_GGMBANK   !FURUTA

C for NONSHARED_TALLY
!$      if(italsh .eq. 0) then
!$      call DEALLOCATE_TAL0
!$      end if


C for NONCRITICAL_OVLY12
C for STRICT_RSOUIN
!$    call DEALLOCATE_MODRSOUIN
      return


      end

************************************************************************
*                                                                      *
      subroutine rcldol(ircss,ncol,mark)
*                                                                      *
*     do repeated collisions                                           *
*                                                                      *
*     2019/10/20 : last revised by K.Niita                             *
*                                                                      *
*        ircss =  0 : without repeated collisions                      *
*              =  1 : repeated collisions                              *
*                                                                      *
************************************************************************

      use MMBANKMOD  !FURUTA
      use MEMBANKMOD !FURUTA
      use moddas_region
      use moddas_repeated_collisions

*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'param00.inc'

      common /bnkmem/ maxbnk, maxbn2, rtrckflp
      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
!$OMP THREADPRIVATE(/icomon/)
      common /jcomon/ nabov,nobch,nocas,nomax
!$OMP THREADPRIVATE(/jcomon/)
      common /tlgeom/ iblz1,iblz2
!$OMP THREADPRIVATE(/tlgeom/)
      common /tlglat/ ilev1,ilev2,ilat1(5,10),ilat2(5,10)
!$OMP THREADPRIVATE(/tlglat/)

      common /poabs/  iabsms
!$OMP THREADPRIVATE(/poabs/)

      common /mathzn/ mathz, mathn, jcoll, kcoll
!$OMP THREADPRIVATE(/mathzn/)

      common /wtsave/ oldwt
!$OMP THREADPRIVATE(/wtsave/)

      common /clustf/ nclst, iclust(nnn)
!$OMP THREADPRIVATE(/clustf/)
      common /clustg/ jclust(0:8,nnn), qclust(0:12,nnn)
!$OMP THREADPRIVATE(/clustg/)
      common /clustt/ nclsts, iclusts(nnn)
!$OMP THREADPRIVATE(/clustt/)
      common /clustw/ jclusts(0:8,nnn),  qclusts(0:12,nnn)
!$OMP THREADPRIVATE(/clustw/)

      common /clustv/ kdecay(4)
!$OMP THREADPRIVATE(/clustv/)

      common /regdc/  idrg(kvlmax), idgr(kvmmax)

      dimension qcl(0:12,nnn)

*-----------------------------------------------------------------------

      common /wtcntl/ iwt, icimp(20), ifcls(20), iwwin(20), ircls(20)

      common /rclmsg/ ircln, isrcl, maxrr,
     &                mnrcl(6,0:20), mrrcl(kvlmax),
     &                krcls(6), lrcls(6), inrlc(6), inrlt(6),
     &                irpem(6,2), irman(6), irmct(6), jsmat(6)
      common /rclemm/ erpem(6,2)

      dimension     idas(1)
      equivalence ( das, idas )

*-----------------------------------------------------------------------

            ircss = 0

            if( iabsms .ne. 0 ) return
            if( ircln .eq. 0  ) return
            if( ircls(ityp) .eq. 0 ) return
            if( mrrcl(idgr(iblz1)) .eq. 0 ) return

            if( nclst .lt. 0 ) return

            if( mathz .le. 2 .or.
     &        (jcoll .ne. 4 .and. jcoll .ne. 5 .and. jcoll .ne. 15 .and.
     &         jcoll .ne. 16 .and. jcoll .ne. 17 .and. jcoll .ne. 19 ))
     &       return

*-----------------------------------------------------------------------

                     icd = 0

            do 100 j = 1, ircln

*-----------------------------------------------------------------------
*              check particle
*-----------------------------------------------------------------------

                     icp = 0

                  do l = 1, mnrcl(j,20)

                     if( mnrcl(j,l) .eq. ityp ) icp = 1

                  end do

                     if( icp .eq. 0 ) goto 100

*-----------------------------------------------------------------------
*              check energy
*-----------------------------------------------------------------------

                     eein  = ec(ibkec+no,ipomp+1)

                  if( irpem(j,1) .gt. 0 .and.
     &                eein .lt. erpem(j,1) ) return

                  if( irpem(j,2) .gt. 0 .and.
     &                eein .gt. erpem(j,2) ) return

*-----------------------------------------------------------------------
*              check mother
*-----------------------------------------------------------------------

                  if( irman(j) .gt. 0 ) then

                     do jj = 1, irman(j)

                        iaz = ismat_jsmat( jsmat(j) + jj - 1 )
                        iz = iaz / 1000
                        ia = iaz - iz * 1000

                        if( ia .eq. 0 ) then
                           if( iz .eq. mathz ) goto 400
                        else
                           if( iz .eq. mathz .and.
     &                         ia .eq. mathz + mathn ) goto 400
                        end if

                     end do

                        return

                  end if

  400                continue

*-----------------------------------------------------------------------
*              check region
*-----------------------------------------------------------------------

                        idsm = inrlc(j)
                        jdsm = 0

                        kdsm = krcls(j)
                        ldsm = lrcls(j)
                  iflag_evap = 0
                  if( allocated(idas_lrcls) ) then
                     if( j < 6 ) then
                        if( lrcls(j+1)-lrcls(j) > 0 ) then
                           iflag_evap = 1
                        end if
                     else
                        if( idas_lrcls(lrcls(j)) > 0 ) then
                           iflag_evap = 1
                        end if
                     end if
                  end if

                  do i = 1, mnrcl(j,0)

                        jj = 0

                        jdsm = jdsm + 1
                        ntrn = idas_inrlc(idsm+jdsm)
                        jdsm = jdsm + 1
                        mtrn = idas_inrlc(idsm+jdsm)

                        jrcls = idas_krcls(kdsm-1+i)
                        iecls = 1
                        if( iflag_evap == 1 )
     &                     iecls = idas_lrcls(ldsm-1+i)

                     do 200 k = 1, ntrn

                        call tregck(iblz2,ilev2,ilat2,
     &                              mtrn,idas_inrlc(idsm+jdsm+1),jj,icc)

                        if( icc .ne. 0 .and. jrcls .ne. 0 ) then

                              icd = icd + 1

                           do l = 1, mnrcl(j,20)

                              if( mnrcl(j,l) .eq. ityp ) goto 300

                           end do

                              goto 100

                        end if

  200                continue

                        jdsm = jdsm + mtrn

                  end do

  100       continue

                  if( icd .eq. 0 ) return

  300       continue

*-----------------------------------------------------------------------

               if( abs(jrcls) * iecls .eq. 1 ) return

               ircss = 1

*-----------------------------------------------------------------------
*           weight reduction
*-----------------------------------------------------------------------

            wt(ibkwt+no,ipomp+1) = wt(ibkwt+no,ipomp+1)  /
     &                                  dble( abs(jrcls) * iecls )

            oldwt = wt(ibkwt+no,ipomp+1)

*-----------------------------------------------------------------------

                        call dataup(ncol,0)
                        call analyz(ncol,mark)

                        call wwindw(ncol)
                        call wtcutof(ncol)
                        call update(2)

*-----------------------------------------------------------------------

                     do jj = 1, iecls - 1

                        call nevap(0)

                        call dataup(ncol,0)
                        call analyz(ncol,mark)

                        call wwindw(ncol)
                        call wtcutof(ncol)
                        call update(2)

                     end do

*-----------------------------------------------------------------------

               do ii = 1, abs(jrcls) - 1

                    ipim = 0

   22             continue

                     call ncasc(1,ityp,ktyp,eein,
     &                          mathz+mathn,mathz,bmax)

                  if( nclst .lt. 0 ) then

                     ipim = ipim + 1

                     if( ipim .le. 20 ) goto 22

                  end if

                     do jj = 1, iecls

                        call nevap(0)

                        call dataup(ncol,0)
                        call analyz(ncol,mark)

                        call wwindw(ncol)
                        call wtcutof(ncol)
                        call update(2)

                     end do

               end do

*-----------------------------------------------------------------------

                        x(ibkx+no,ipomp+1) = xc(ibkxc+no,ipomp+1)
                        t(ibkt+no,ipomp+1) = tc(ibktc+no,ipomp+1)
                        y(ibky+no,ipomp+1) = yc(ibkyc+no,ipomp+1)
                        z(ibkz+no,ipomp+1) = zc(ibkzc+no,ipomp+1)
                        e(ibke+no,ipomp+1) = ec(ibkec+no,ipomp+1)

*-----------------------------------------------------------------------

      return

      end

************************************************************************
*                                                                      *
      subroutine dmprijk(ircss,ncol,mark)
*                                                                      *
*     2024/12/06 : last revised by T.Furuta                            *
*                                                                      *
************************************************************************
      implicit real*8(a-h,o-z)
      common /mpi00/ npe, me
      common /randtp/ rijk,rans,ranb, iransb64
!$OMP THREADPRIVATE(/randtp/)
      character(100) :: filnm
      character(5) :: chme
      if(me.gt.0)then
         iorder=aint(log10(real(npe)))+1
         if(iorder.lt.3)iorder=3
         write(chme,'(i5.5)')me
         filnm='rijkdmp'//chme(6-iorder:5)//'.inp'
      else
         filnm='rijkdmp.inp'
      endif
      open(801,file=filnm,form='formatted',status='unknown')
      write(801,'(" bitrseed = ",b64.64)')rijk
      close(801)
      return
      end

