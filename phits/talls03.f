************************************************************************
*                                                                      *
      subroutine tyilreg(ncol,m,mz,mn,mm,nl,lt, ! frtati 2022/02/18 added mm
     &                   nr,mr,nm,kr,mt,ikzz,iknn,tr,
     &                   trEVENT)
*                                                                      *
*       nuclear yield (or production) tally in region mesh             *
*       last modified by K.Niita on 2011/11/11                         *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*      ncol  ..... reaction type                                       *
*             11 : termination by energy cut-off                       *
*             13 : nuclear reaction (n,x)                              *
*             14 : nuclear reaction (n,n'x)                            *
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
*        kcoll : =  0, normal                                          *
*                =  1, high energy fission                             *
*                =  2, high energy absorption                          *
*                =  3, low energy n elastic                            *
*                =  4, low energy n non-elastic                        *
*                =  5, low energy n fission                            *
*                =  6, low energy n absorption                         *
*                                                                      *
************************************************************************
      use QMD_COOD2_MOD, only : nqmdm
      use MMBANKMOD !FURUTA
      use partmod, only: itmxpt,itpan,itpat,jtpat,iznmmx,iznmturn ! frtati 2022/05/02
      use moddas_region ! S.H. (2022.11.9)
*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'param00.inc'
      include 'param01.inc'

*-----------------------------------------------------------------------

      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
!$OMP THREADPRIVATE(/icomon/)
      common /wtsave/ oldwt
!$OMP THREADPRIVATE(/wtsave/)
      common /tlgeom/ iblz1,iblz2
!$OMP THREADPRIVATE(/tlgeom/)
      common /tlglat/ ilev1,ilev2,ilat1(5,10),ilat2(5,10)
!$OMP THREADPRIVATE(/tlglat/)
      common /qparm/  ielas,icasc,iqstep,lvlopt,igamma
      common /bparm/  andt,jevap,npidk
      common /geosig/ geosig(250)
      common /cparm/  maxbch,maxcas
      common /mathzn/ mathz, mathn, jcoll, kcoll
!$OMP THREADPRIVATE(/mathzn/)
      common /eparm/  esmax, esmin, emin(20)

      common /regdc/  idrg(kvlmax), idgr(kvmmax)
*-----------------------------------------------------------------------

      common /clustf/ nclst, iclust(nnn)
!$OMP THREADPRIVATE(/clustf/)
      common /clustt/ nclsts, iclusts(nnn)
!$OMP THREADPRIVATE(/clustt/)
      common /clustw/ jclusts(0:8,nnn),  qclusts(0:12,nnn)
!$OMP THREADPRIVATE(/clustw/)
      common /clustl/ rumpal(0:20), numpal(0:20)
!$OMP THREADPRIVATE(/clustl/)
      common /cntcls/ jcount(3,nnn)
!$OMP THREADPRIVATE(/cntcls/)

*-----------------------------------------------------------------------

      dimension numsav(0:20), rumsav(0:20)

*-----------------------------------------------------------------------

      common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                itout(itlmax), ittwo(itlmax)
      common /tall12/ itman(itlmax), itmat(itlmax), itmct(itlmax),
     &                itnun(itlmax), itnuc(itlmax),
     &                itndz(itlmax), itndn(itlmax),
     &                itnkz(itlmax), itnkn(itlmax)
      common /tall22/ itmcn(itlmax), itmtn(itlmax), itmtt(itlmax)
      common /tall46/ itmtr(itlmax,4), rtmtr(itlmax,13)

      common /tall30/ itnda(itlmax)
      common /tall36/ itdpo(itlmax)
      common /tall37/ itcnt(9,itlmax)
      common /tall52/ itprd(itlmax)

      common /trstar/ itrstar
!$OMP THREADPRIVATE(/trstar/)

      common /tall82/ itcnth(9,itlmax)
      common /tall83/ itnzn(itlmax), itndm(itlmax)

      common /kmat1d/ idmn(0:kvlmax), idnm(kvmmax)
      common /clionprd/  lionprd

*-----------------------------------------------------------------------

      dimension   lt(nl)
      dimension   kr(mr)
      dimension   mt(nm)
      dimension   ikzz(maxpt,maxnt), iknn(maxpt,maxnt)
      dimension   tr(nr,mz,mn,0:mm,2)
      dimension   trEVENT(nr,mz,mn,0:mm)           !OBINATA(2012.8.20): as Ct
      real(8),allocatable,save:: tr0(:,:,:,:) !OBINATA(2012.8.20): as C

      dimension mnz(mxprodxs),mna(mxprodxs),xxn(mxprodxs),ildd(mxprodxs)

      dimension   ips(itmxpt) ! frtati 2021/10/05 6 -> itmxpt
      dimension   ncntt(3)

*-----------------------------------------------------------------------
      common /cntreg/ ncntc(3), ncreg(3), incrc(3), incrt(3), kcont(3)
      common /cntpat/ icpan(3), icpat(3,mxcntprt,2) ! S.H. set mxcntprt (2022.3.24)
      dimension     idas(1)
      equivalence ( das, idas )
      common /pnsave/ egs, uus, vvs, wws, wts, tms, nms, nct(3)
!$OMP THREADPRIVATE(/pnsave/)

*-----------------------------------------------------------------------
      common /cntmx/ ncntmx(3), nctmxsr(3)
!$OMP THREADPRIVATE(/cntmx/)

      common /mpi00/ npe, me
      real*8,allocatable :: tryld(:,:)

*-----------------------------------------------------------------------

      common /jcomon/ nabov,nobch,nocas,nomax
!$OMP THREADPRIVATE(/jcomon/)
      common /stat / istdev, irestart, ireschk

      integer, allocatable, save :: iclusav(:), jclusav(:,:)
      double precision, allocatable, save :: qclusav(:,:)
!$OMP THREADPRIVATE(iclusav, jclusav, qclusav)

       if( .not. allocated(iclusav) ) then ! initial allocation
        allocate(iclusav(nqmdm),jclusav(0:8,nqmdm),qclusav(0:12,nqmdm))
        iclusav = 0
        jclusav = 0
        qclusav = 0.d0
       elseif( ubound(iclusav,1) .lt. nqmdm ) then ! extend array
        deallocate(iclusav,jclusav,qclusav)
        allocate(iclusav(nqmdm),jclusav(0:8,nqmdm),qclusav(0:12,nqmdm))
        iclusav = 0
        jclusav = 0
        qclusav = 0.d0
       endif

*-----------------------------------------------------------------------
* check of history counter
*-----------------------------------------------------------------------
      ihistcount = 0 ! history counter index, T.Sato 2022/12/20 for speed up
         if ( ncol.eq.0 .or. ncol.eq.4 ) then
           do i = 1, 3
             if( itcnth(i,m) .eq. 1 ) then
               if( ncntmx(i) .lt. itcnth(i*2+2,m) .or.
     &             ncntmx(i) .gt. itcnth(i*2+3,m) ) then
                 ihistcount = 1
               end if
             end if
           end do
         end if

*-----------------------------------------------------------------------
*        source particle or end of batch ( in case of istdev = 2 )
*
* OBINATA(2012.8.20): change tr(,,,,3) to trEVENT(,,,,)
*-----------------------------------------------------------------------

         if (( ncol .eq. 0 .or. ncol .eq. 4 )
     &                              .and. istdev .eq. 2) then
           if ((nocas.gt.1.or.ncol.eq.0) .and. ihistcount.ne.1 ) then

             tr(:,:,:,:,1) = tr(:,:,:,:,1) + trEVENT(:,:,:,:)
             tr(:,:,:,:,2) = tr(:,:,:,:,2) + trEVENT(:,:,:,:) ** 2

           end if

           trEVENT(:,:,:,:) = 0

         end if

*-----------------------------------------------------------------------
*        end of batch ( in case of istdev = 1 )
*
* OBINATA(2012.8.20): modificate for thread parallel
*-----------------------------------------------------------------------

         if ( ncol .eq. 0 .and. istdev .eq. 1) then
!$OMP MASTER
             allocate( tr0(nr,mz,mn,0:mm) ) ! frtati 2022/02/18 2 -> mm
             tr0(:,:,:,:) = 0.d0
!$OMP END MASTER
!$OMP BARRIER
!$OMP CRITICAL (tyilreg_crit_ist1)
             tr0(:,:,:,:) = tr0(:,:,:,:) + trEVENT(:,:,:,:)
!$OMP END CRITICAL (tyilreg_crit_ist1)
!$OMP BARRIER
!$OMP MASTER
             tr(:,:,:,:,1) = tr(:,:,:,:,1) + tr0(:,:,:,:) / maxcas
             tr(:,:,:,:,2) = tr(:,:,:,:,2)
     &                   + ( tr0(:,:,:,:) / maxcas ) ** 2
             deallocate( tr0 )
!$OMP END MASTER

           trEVENT(:,:,:,:) = 0

         end if

*-----------------------------------------------------------------------
*        rearrangement of iznm for MPI run ! frtati 2022/05/02
*-----------------------------------------------------------------------

         if( npe.gt.1 .and. itnzn(m).ne.0 .and. ncol.eq.0 ) then
           call paraiznm(m,1)
           if( me.gt.1 ) then
             allocate( tryld(iznmmx(m),2) )
             do ir = 1, nr
               do iz = 1, iznmmx(m)
                 tryld(iz,:) = tr(ir,iz,1,0,:)
               end do
               do iz = 1, iznmmx(m)
                 tr(ir,iz,1,0,:) = tryld(iznmturn(iz,m),:)
               end do
             end do
             deallocate( tryld )
           end if
         end if

*-----------------------------------------------------------------------
*        check of ncol and nclsts ( outgoing particles )
*-----------------------------------------------------------------------

         if( ncol .eq. 11 ) then

            if( itprd(m) .ne. 1 ) return

               iccol = 0
               npart = 1

               if( ityp .lt. 15 .or. ityp .gt. 19 ) return

               jz    = jtyp
               jn    = ktyp - ktyp / 1000000 * 1000000 - jz

               ipart = ityp
               tlw   = oldwt
               il = 0  ! T.Sato 2016/3/29, isomer is not considered for output=cutoff
               ncntt(1) = ncnt(ibknct+1,no,ipomp+1)
               ncntt(2) = ncnt(ibknct+2,no,ipomp+1)
               ncntt(3) = ncnt(ibknct+3,no,ipomp+1)

         else if( ncol .eq. 13 .or. ncol .eq. 14 ) then

            if ( itnda(m) .ge. 2 .and. jcoll .eq. 9 ) then
               if( nclsts .lt. 0 ) return
            else
               if( nclsts .le. 0 ) return
            end if

            if( ( itnda(m) .eq. 2 .and.
     &          ( jcoll .eq. 6 .or. jcoll .eq. 9 ) ) .or.
     &          ( itnda(m) .eq. 3 .and.
     &          ( jcoll .eq. 6 .or. jcoll .eq. 9 .or.
     &            jcoll .eq. 4 .or. jcoll .eq. 5 .or.
     &            jcoll .eq. 15 .or. jcoll .eq. 10 ) ) ) then

            else
               if( itdpo(m) .eq. 0 .and.
     &           ( jcoll .eq. 3 .or. kcoll .eq. 3 ) ) return
               if( itdpo(m) .eq. -1 .and. itrstar .eq. 1 ) return
            end if

               iccol = 1

               npart = nclsts

               ncntt(1) = jcount(1,1)
               ncntt(2) = jcount(2,1)
               ncntt(3) = jcount(3,1)

         else

            return

         end if

*-----------------------------------------------------------------------
*        check of counter
*-----------------------------------------------------------------------

            do i = 1, 3

               if( itcnt(i,m) .ne. 0 ) then

                  if( ncntt(i) .lt. itcnt(i*2+2,m) .or.
     &                ncntt(i) .gt. itcnt(i*2+3,m) ) return

               end if

            end do

*-----------------------------------------------------------------------
*        check of mat
*-----------------------------------------------------------------------

            if( nl .gt. 0 ) then

                  do i = 1, nl

                     if( itmcn(m) .gt. 0 .and.
     &                   idmn(mat) .eq. lt(i) ) goto 502
                     if( itmcn(m) .lt. 0 .and.
     &                   idmn(mat) .eq. lt(i) ) return

                  end do

                     if( itmcn(m) .gt. 0 ) return

            end if

  502          continue

*-----------------------------------------------------------------------

         if( iccol .eq. 1 ) then

*-----------------------------------------------------------------------
*           check of particles
*-----------------------------------------------------------------------

            call pcheck(m,itpan(m),ityp,ktyp,jtyp,ipn,ips)

               if( ipn .eq. 0 ) return

*-----------------------------------------------------------------------
*           mother and energy
*-----------------------------------------------------------------------

               ata   = dble( mathz + mathn )
               atz   = dble( mathz )
               mmas  = nint( ata )
               mchg  = nint( atz )
               eein  = ec(ibkec+no,ipomp+1)

*-----------------------------------------------------------------------
*           for specific mother nuclei
*-----------------------------------------------------------------------

            if( nm .gt. 0 ) then

               do i = 1, nm

                  iz = mt(i) / 1000
                  ia = mt(i) - iz * 1000

                  if( ( ia .gt. 0 .and.
     &                  mmas .eq. ia .and. mchg .eq. iz ) .or.
     &                ( ia .eq. 0 .and. mchg .eq. iz ) ) then

                     if( itmct(m) .gt. 0 ) goto 30
                     if( itmct(m) .lt. 0 ) return

                  end if

               end do

                  if( itmct(m) .gt. 0 ) return

            end if

   30       continue

         end if

*-----------------------------------------------------------------------
*     check of region
*-----------------------------------------------------------------------

                  jj = 0

*-----------------------------------------------------------------------

      do 100 ii = 1, nr

*-----------------------------------------------------------------------

               call tregck(iblz1,ilev1,ilat1,mr,kr,jj,icc)

               if( icc .eq. 0 ) goto 100

               ir = ii

*-----------------------------------------------------------------------
*        normal weight
*-----------------------------------------------------------------------

               idoy = 1
               idev = 1
               wyld = oldwt

*-----------------------------------------------------------------------
*        for special : repeated nuclear reactions
*-----------------------------------------------------------------------

         if( iccol .eq. 1 .and. itspc(m) .gt. 0 .and.
     &       mathz .gt. 2 .and.
     &     ( jcoll .eq. 4 .or. jcoll .eq. 5 ) ) then

               idoy = itspc(m)
               idev = 10
               wyld = oldwt / dble( idoy * idev )

*-----------------------------------------------------------------------
*           save normal output and reaction parameters
*-----------------------------------------------------------------------

                  nclssav = nclsts

            do i = 1, nclsts

                  iclusav(i) = iclusts(i)

               do k = 0, 8

                  jclusav(k,i) = jclusts(k,i)

               end do

               do k = 0, 12

                  qclusav(k,i) = qclusts(k,i)

               end do

            end do

               do k = 0, 20

                  numsav(k) = numpal(k)
                  rumsav(k) = rumpal(k)

               end do

                  iprj  = ityp
                  kprj  = ktyp

                  if( ityp .lt. 15 ) then

                     bmax = sqrt( geosig(mmas) * 100.0 / 3.1415926 )

                  else if( ityp .ge. 15 ) then

                     ap = dble( ktyp - ktyp / 1000000 * 1000000 )
                     zp = dble( ktyp / 1000000 )

                     call sighi(ap,zp,eein,ata,zta,signe,sigel,bmax)

                  end if

         end if

*-----------------------------------------------------------------------
*     repeated do loop
*-----------------------------------------------------------------------

         do ireac = 1, idoy

*-----------------------------------------------------------------------
*     repeat calculation of nuclear reactions
*-----------------------------------------------------------------------

            if( idoy .gt. 1 .and. ireac .gt. 1 ) then

                  ipim = 0

   22          continue

                  call ncasc(1,iprj,kprj,eein,mmas,mchg,bmax)

               if( nclst .lt. 0 ) then

                     ipim = ipim + 1

                     if( ipim .le. 20 ) goto 22

               end if

            end if

*-----------------------------------------------------------------------
*     repeat calculation of evaporation
*-----------------------------------------------------------------------

            do ievap = 1, idev

               if( idev .gt. 1 ) then

                  call nevap(1)

                  do i = 1, nclsts

                     qclusts(8,i) = wyld * qclusts(8,i)

                  end do

                  npart = nclsts

               end if

*-----------------------------------------------------------------------
*        booking after nuclear reactions
*-----------------------------------------------------------------------

                  istat = 1
                  mnx   = 0

*-----------------------------------------------------------------------
*           NDATA=2,3 : Replace with yield data (activation cross section
*-----------------------------------------------------------------------

               if( iccol .eq. 1 .and.
     &           ( itnda(m) .eq. 2 .and.
     &           ( jcoll .eq. 6 .or. jcoll .eq. 9 ) ) .or.
     &           ( itnda(m) .eq. 3 .and.
     &           ( jcoll .eq. 6 .or. jcoll .eq. 9 .or.
     &             jcoll .eq. 4 .or. jcoll .eq. 5 .or.
     &             jcoll .eq. 15 .or. jcoll .eq. 10 ) ) ) then

                  mkk = mat
                  icl = idgr(iblz1)
                  jcl = jcoll

                  call prodxs2(ktyp,mchg,mmas,eein,
     &                 istat,mnx,mnz,mna,xxn,xxs,
     &                 ildd,mkk,icl,jcl,m)

                  iexclight = 0 ! Exclude light ions
                  if( itprd(m) .eq. 1 ) then ! output = cutoff
                     iexclight = 1
                     if ( jcoll.eq.6 .and.
     &                    ktyp.eq.2112 .and. eein.lt.20.0 )
     &                    iexclight = 0
                  end if

                  if( istat .eq. 0 .and. mnx .gt. 0 ) then

                   do nnx = 1, mnx

                      jz = mnz(nnx)
                      jn = mna(nnx) - mnz(nnx)
                      il = ildd(nnx)

                    if ( iexclight.ne.1 .or. jz.gt.2 ) then

                     if( jz .gt. 0 .and. jz .le. maxpt .and.
     &                      jn .gt. 0 .and. jn .le. maxnt ) then

                        kz = ikzz(jz,jn)
                        kn = iknn(jz,jn)

                      if( itnzn(m).eq.0 ) then ! frtati 2022/02/18
                       if( kz .le. mz .and. kn .le. mn ) then
                          trEVENT(ir,kz,kn,il) =
     &                         trEVENT(ir,kz,kn,il) + wyld * xxn(nnx)
                       end if

                      else
                         trEVENT(ir,igetiznm(kz,kn,il,m),1,0) =
     &                        trEVENT(ir,igetiznm(kz,kn,il,m),1,0) +
     &                        wyld * xxn(nnx)
                      end if

                     end if

                    end if

                   end do

                  end if

               end if

*-----------------------------------------------------------------------
*           nuclear data for 4He, 14N, 16O
*-----------------------------------------------------------------------

               if( iccol .eq. 1 .and.
     &             itnda(m) .eq. 1 .and. jcoll .eq. 4 ) then

                     call prodxs(ktyp,mchg,mmas,eein,
     &                           istat,mnx,mnz,mna,xxn,xxs,ildd)

                  if( istat .eq. 0 .and. mnx .gt. 0 ) then

                     do nnx = 1, mnx

                        jz = mnz(nnx)
                        jn = mna(nnx) - mnz(nnx)
                        il = ildd(nnx)

                        if( jz .gt. 0 .and. jz .le. maxpt .and.
     &                      jn .gt. 0 .and. jn .le. maxnt ) then

                           kz = ikzz(jz,jn)
                           kn = iknn(jz,jn)

                        if( itnzn(m).eq.0 ) then ! frtati 2022/02/18

                        if( kz .le. mz .and. kn .le. mn ) then

                           trEVENT(ir,kz,kn,il) =
     &                     trEVENT(ir,kz,kn,il) + wyld * xxn(nnx)

                        end if

                        else
                           trEVENT(ir,igetiznm(kz,kn,il,m),1,0) =
     &                     trEVENT(ir,igetiznm(kz,kn,il,m),1,0) +
     &                     wyld * xxn(nnx)
                        end if

                        end if

                     end do

                  end if

               end if


*-----------------------------------------------------------------------
*           normal case
*-----------------------------------------------------------------------

            if( npart .gt. 0 ) then

               do 500 j = 1, npart

*-----------------------------------------------------------------------

                  if( iccol .eq. 1 ) then

                     ipart = jclusts(3,j)
                     tlw   = qclusts(8,j)
                     jz    = jclusts(1,j)
                     jn    = jclusts(2,j)
                     il    = jclusts(8,j)

                  end if

*-----------------------------------------------------------------------
* When ndata=2 or 3, counter values change by ndata in [counter]
                 if ( (ipart.ge.3 .and. ipart.le.13) .or.
     &                 ipart.ge.19 ) then

                  if( ( itnda(m) .eq. 2 .or. itnda(m) .eq. 3 ) .and.
     &                 ( istat .eq. 0 .or. istat .eq. 3 ) ) then

                   if( ncntc(1) .eq. 1 .or. ncntc(2) .eq. 1 .or.
     &                    ncntc(3) .eq. 1 ) then
                    do k = 1, 3
                     ndatcount = 0
                     if( ncntc(k) .eq. 1 ) then
                      knn  = ityp
                      kcg  = jtyp
                      kkf  = ktyp
                      call pcchck(k,knn,kkf,kcg,icpan,icpat,icc)

                      if( icc .eq. 1 ) then
                       kdsm = kcont(k)
                       ldsm = 22
                       ndatcount = idas_kcont(kdsm+ldsm)
                      end if

                      if ( ndatcount .ne. 0 ) then

                       if ( ncol .eq. 13 ) then
                        ncnta(ibknct+k,j,ipomp+1)
     &                         = ncnta(ibknct+k,j,ipomp+1)
     &                         + ndatcount

                       else if ( ncol .eq. 14 ) then
                        if ( j .eq. 1 ) then
                         nct(k) = nct(k) + ndatcount
                        else if (j .gt. 1 ) then
                         ncnta(ibknct+k,j-1,ipomp+1)
     &                          = ncnta(ibknct+k,j-1,ipomp+1)
     &                          + ndatcount
                        end if

                       end if

                      end if

                     end if
                    end do
                   end if

                  end if

                 end if

*-----------------------------------------------------------------------

                  if( ipart .lt. 15 .or. ipart .gt. 19 ) goto 500

                  if( iccol .eq. 1 .and. itprd(m) .eq. 1 ) then

                     emint = emin(ipart) * dble( jz + jn )

                     if( qclusts(7,j) .gt. emint ) goto 500

                  end if

*-----------------------------------------------------------------------

                     if( ( itnda(m) .eq. 2 .or. itnda(m) .eq. 3 ) .and.
     &                 istat .eq. 0 .and.
     &                 ( iexclight.ne.1 .or. jz.gt.2 ) ) goto 500

                     if( ( itnda(m) .eq. 2 .or. itnda(m) .eq. 3 ) .and.
     &                    istat .eq. 3 .and. lionprd .ne. 1) goto 500

                     if( itnda(m) .eq. 1 .and.
     &                   istat .eq. 0 .and. mnx .gt. 0 ) then

                        do nnx = 1, mnx

                           if( jz .eq. mnz(nnx) .and.
     &                         jn .eq. mna(nnx) - mnz(nnx) ) goto 500

                        end do

                     end if

                        if( jz .gt. 0 .and. jz .le. maxpt .and.
     &                      jn .gt. 0 .and. jn .le. maxnt ) then

                           kz = ikzz(jz,jn)
                           kn = iknn(jz,jn)

                        if( itnzn(m).eq.0 ) then ! frtati 2022/02/18

                        if( kz .le. mz .and. kn .le. mn ) then
cOBINATA(2012.8.20): Ct = Ct + xi.wi
                           trEVENT(ir,kz,kn,il) = trEVENT(ir,kz,kn,il)
     &                                       + tlw
                        end if

                        else
                          trEVENT(ir,igetiznm(kz,kn,il,m),1,0) =
     &                    trEVENT(ir,igetiznm(kz,kn,il,m),1,0) + tlw
                        end if

                        end if

  500          continue

            end if

*-----------------------------------------------------------------------

            end do

         end do

*-----------------------------------------------------------------------
*        restore normal output
*-----------------------------------------------------------------------

         if( idoy .gt. 1 ) then

                  nclsts = nclssav

            do i = 1, nclsts

                  iclusts(i) = iclusav(i)
                  iclusav(i) = 0

               do k = 0, 8

                  jclusts(k,i) = jclusav(k,i)
                  jclusav(k,i) = 0

               end do

               do k = 0, 12

                  qclusts(k,i) = qclusav(k,i)
                  qclusav(k,i) = 0.d0

               end do

            end do

               do k = 0, 20

                  numpal(k) = numsav(k)
                  rumpal(k) = rumsav(k)

               end do

         end if

*-----------------------------------------------------------------------
*        store the decrease of mother by negative weight
*-----------------------------------------------------------------------

        if( ncol .ne. 11 ) then ! S.H. 2021.12.22
         if( itdpo(m) .eq. -1 .and. itrstar .eq. 0 ) then

            if( mathz .gt. 0 .and. mathz .le. maxpt .and.
     &          mathn .gt. 0 .and. mathn .le. maxnt ) then

               kz = ikzz(mathz,mathn)
               kn = iknn(mathz,mathn)

             if( itnzn(m).eq.0 ) then ! frtati 2022/02/18
              if( kz .le. mz .and. kn .le. mn ) then ! T.Sato 2024/05/02
               trEVENT(ir,kz,kn,0) = trEVENT(ir,kz,kn,0) - oldwt
              endif
             else
               trEVENT(ir,igetiznm(kz,kn,0,m),1,0) =
     &         trEVENT(ir,igetiznm(kz,kn,0,m),1,0) - oldwt
             end if

            end if

         end if
        end if

*-----------------------------------------------------------------------

  100 continue

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine pyildreg(m,mz,mn,mm, ! frtati 2022/02/18 added mm
     &                    nr,mr,nn,kr,nt,ikzz,iknn,
     &                    tr,nvl,ivl,rvl,
     &                    nx,ny,nz,xm,ym,zm,igsh,idasa)
*                                                                      *
*       output of nuclear yield (or production) tally in region mesh   *
*       last modified by K.Niita on 2004/09/15                         *
*                                                                      *
************************************************************************
      use sangelmod, only: itsans

*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'param01.inc'

*-----------------------------------------------------------------------

      common /taliin/ rsouin, nzztin, nrgnin

*-----------------------------------------------------------------------

      common /verjam/ versn, lastr, iyeav, imonv, idayv

      common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                itout(itlmax), ittwo(itlmax)
      common /tall02/ itrgn(itlmax), itrgm(itlmax), itreg(itlmax)
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

      common /tall21/ rtfac(itlmax)
      common /tall29/ itrsh(itlmax)
      common /tall35/ iteps(itlmax)
      common /tall38/ itres(itlmax)
      common /tall39/ rtwid(itlmax)
      common /tall46/ itmtr(itlmax,4), rtmtr(itlmax,13)
      common /tall49/ itglt(itlmax)
      common /tall83/ itnzn(itlmax), itndm(itlmax)

      common /istcut/ ist_cut, ist_bat

      common /fact01/ facmax(itlmax) ! kitamura23/03/31

*-----------------------------------------------------------------------

      common /volreg/ dvol(kvlmax)
      common /regdc/ idrg(kvlmax), idgr(kvmmax)

*-----------------------------------------------------------------------

      common /tcntl/  icntl, inucr

*-----------------------------------------------------------------------

      common /cparm/  maxbch,maxcas
      common /jcomon/ nabov,nobch,nocas,nomax
!$OMP THREADPRIVATE(/jcomon/)
* Use igamma
      common /qparm/  ielas,icasc,iqstep,lvlopt,igamma

      common /talout/ itall

      character fname*100, fnume*3

*-----------------------------------------------------------------------
      character gfnam*100
      common /tall60/ itger(itlmax)

      common /tall67/ itstd(itlmax), rtstd(itlmax)

      common /mpi00/ npe, me

*-----------------------------------------------------------------------

      dimension   kr(mr)
      dimension   nt(nn)
      dimension   ikzz(maxpt,maxnt), iknn(maxpt,maxnt)
      dimension   tm(maxnt+maxpt,2)
      dimension   vl(nr)
      dimension   lr(nr)
      dimension   tr(nr,mz,mn,0:mm,2) ! frtati 2022/02/18 2 -> mm
      dimension   ivl(nvl)
      dimension   rvl(nvl)

      dimension   xm(nx+1)
      dimension   ym(ny+1)
      dimension   zm(nz+1)
      dimension   val(nr)

      integer,allocatable :: ixyz(:)

*-----------------------------------------------------------------------

      character hsunit(2)*15

      data hsunit / '[1/source]     ',
     &              '[1/cm^3/source]'/

*-----------------------------------------------------------------------

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

      data ipstep / 12 /

      character erfnm*100

      character chau*8
      character cha*1
      data cha /"'"/

*-----------------------------------------------------------------------
      character(28) cfmt                !FURUTA20200615
      character(12) cir                 !FURUTA20200615
      common /redufmt/ iredufmt(itlmax) !FURUTA20200615
*-----------------------------------------------------------------------
      character rpa*1
      data rpa /'}'/
      character yen*1
      yen  = char(92)
      allocate( ixyz(max((nx+1)*(ny+1),(ny+1)*(nz+1),(nx+1)*(nz+1))) )

*-----------------------------------------------------------------------
*        set mesh volume and name
*-----------------------------------------------------------------------

               call tregvl(nr,mr,kr,vl,lr,nvl,ivl,rvl)

*-----------------------------------------------------------------------
*        relative error and  unit conversion
*        c1 : nomalization for source
*-----------------------------------------------------------------------

                  cmax = 0.0
                  cmin = 1.e+33
                  dnon = 1.e-33
                  isdz = 0
                  stdm = 0.0
                  itstd(m) = 1

         if( igsh .eq. 0 ) then

            if ( rsouin .gt. 0d0 ) then ! S.H. (2016.7.26)

               c1 = 1.0d+0 / rsouin

            else

               c1 = 0d0

            end if

            facmax(m) = 1.d0
            if( rtfac(m) .lt. 0.d0 ) then
               facmax(m) = 0.d0

               do 101 ir = 1, nr
               do 101 iz = 1, mz
               do 101 in = 1, mn
               do 101 il = 0, mm

                  if( tr(ir,iz,in,il,1) .ne. 0.d0 ) then

                     if( itunt(m) .eq. 1 ) then

                        fmaxfc = tr(ir,iz,in,il,1)

                     else if( itunt(m) .eq. 2 ) then

                        fmaxfc = tr(ir,iz,in,il,1) / vl(ir)

                     end if

                     if( facmax(m) .lt. fmaxfc) facmax(m) = fmaxfc

                  end if

  101          continue

            end if

            do 100 ir = 1, nr

               if( itunt(m) .eq. 1 ) then

                  cc = abs(rtfac(m)/facmax(m))

               else if( itunt(m) .eq. 2 ) then

                  cc = abs(rtfac(m)/facmax(m)) / vl(ir)

               end if

            do 100 iz = 1, mz
            do 100 in = 1, mn
            do 100 il = 0, mm ! frtati 2022/02/18

               if( tr(ir,iz,in,il,1) .ne. 0.d0 ) then

                  call calc_stdev(m,Xa,sigx,
     &                            tr(ir,iz,in,il,1),
     &                            tr(ir,iz,in,il,2),
     &                            cc)

                  tr(ir,iz,in,il,1) = Xa
                  tr(ir,iz,in,il,2) = sigx
                  if( sigx .gt. stdm ) stdm = sigx

                  if( tr(ir,iz,in,il,1) .gt. cmax )
     &                                   cmax = tr(ir,iz,in,il,1)

cABE 2022/03/02, avoid the negative value of cmin
                  if( tr(ir,iz,in,il,1) .gt. 0.d0 .and.
     &                tr(ir,iz,in,il,1) .lt. cmin )
     &                                   cmin = tr(ir,iz,in,il,1)

               else

                  if( il .eq. 0 ) isdz = 1
                  tr(ir,iz,in,il,2) = 0.0

               end if

  100       continue

            if( nobch .gt. ist_bat ) then
            if( rtstd(m) .gt. 0.d0 .and.
     &          stdm .gt. 0.d0 .and. isdz .eq. 0 .and.
     &          stdm .lt. rtstd(m) ) itstd(m) = 0
            end if

         end if

*-----------------------------------------------------------------------
*     out put unit = 31  : temporary number
*-----------------------------------------------------------------------

      do 900 iax = 1, itfln(m)

!OBINATA(2012.8.20): output *.err
        noe = 1
        if ( any( itaxs(m,iax) .eq. (/ 8, 13 /) )
     &      .and. ittwo(m) .ne. 4 .and. icntl .ne. 8 )
     &    noe = 2

        if(itaxs(m,iax).eq.13.and.iredufmt(m).eq.1) noe=1 !FURUTA20200615

        do 900 ioe = 1, noe

         if( igsh .ne. 0 .and.
     &     ( itaxs(m,iax) .lt. 7 .or. itrsh(m) .eq. 0 ) ) goto 900

         if( itall .eq. 2 .and. nobch .lt. maxbch .and.
     &       igsh .eq. 0 ) then

            write(fnume,'(i3.3)') nobch
!OBINATA(2012.8.20): output *.err
            if ( ioe .eq. 1 ) then
              fname = ctfln(m,iax)(1:itfll(m,iax))//'.'//fnume
            else
              call mk_2derrfn(ctfln(m,iax),fname,itfll(m,iax))
            end if

         else if ( itall.eq.4 .and. igsh.eq.0 ) then
            write(fnume,'(i3.3)') nobch
            if ( npe.gt.1 ) write (fnume,'(i3.3)') nobch/(npe-1)
            if ( ioe.eq.1 ) then
              call mk_2dnumfn(ctfln(m,iax),fname,itfll(m,iax),
     &             nobch,maxbch,npe)
            else
              call mk_2dnerfn(ctfln(m,iax),fname,itfll(m,iax),fnume)
            end if

         else

!OBINATA(2012.8.20): output *.err
            if ( ioe .eq. 1 ) then
              fname = ctfln(m,iax)
            else
              call mk_2derrfn(ctfln(m,iax),fname,itfll(m,iax))
            end if

         end if

            iot = 31
            open(iot, file = fname, status = 'unknown' )

*-----------------------------------------------------------------------
            gfnam = ctfln(m,iax)
            igfmn = itfll(m,iax)
            igser = itger(m)

*-----------------------------------------------------------------------
*     input echo
*-----------------------------------------------------------------------

               call tyilech(iot,m,iax,1)

*-----------------------------------------------------------------------
*        dchain axis
*-----------------------------------------------------------------------

         if( itaxs(m,iax) .eq. 13 ) then

*-----------------------------------------------------------------------
*           output error or not
*-----------------------------------------------------------------------

            if( itout(m) .ne. 0 ) then

               do i = itfll(m,iax), 1, -1

                  if( ctfln(m,iax)(i:i) .eq. '.' ) goto 50

               end do

                  i = itfll(m,iax)

   50             itfp = i - 1

               do i = 1, itfp

                  erfnm(i:i) = ctfln(m,iax)(i:i)

               end do

                  erfnm(itfp+1:itfp+4) = '.err'

                  iou = 15
                  open(iou, file = erfnm(1:itfp+4), status = 'unknown' )

                  call tyilech(iou,m,iax,1)

            end if

*-----------------------------------------------------------------------
*        output : do not change below expresion
*-----------------------------------------------------------------------

            write(iot,'(/"#",78("-"))')

            write(iot,'(/
     &           " regionwise nuclear yield (or production)"/
     &           " ----------------------------------------")')

            if( itout(m) .ne. 0 ) then

               write(iou,'(/"#",78("-"))')

               write(iou,'(/
     &         " Statistical Error(%) for regionwise nuclear ",
     &         "yield of above file."/
     &         " --------------------------------------------",
     &         "---------------------")')

            end if

            if(iredufmt(m).eq.0)then !FURUTA20200615
*-----------------------------------------------------------------------

             do 170 il = 0, 2 ! frtati 2022/03/11
            do 170 iz = 1, maxpt

               if( nn .gt. 0 ) then

                  do i = 1, nn

                     if( nt(i) / 1000 .eq. iz ) goto 150

                  end do

                     goto 170

               end if

  150          continue

                     do in = 1, maxnt
                     do ir = 1, nr

                        if( itnzn(m).eq.0 ) then
                          t0_tr = tr(ir,iz,in,il,1)
                        else
                          t0_tr = tr(ir,igetiznmp(iz,in,il,m),1,0,1)
                        end if
                        if( t0_tr .ne. 0.0 ) goto 120

                     end do
                     end do

                     goto 170

  120             im = in

                     do in = maxnt, im + 1, -1
                     do ir = 1, nr

                        if( itnzn(m).eq.0 ) then
                          t0_tr = tr(ir,iz,in,il,1)
                        else
                          t0_tr = tr(ir,igetiznmp(iz,in,il,m),1,0,1)
                        end if
                        if( t0_tr .ne. 0.0 ) goto 140

                     end do
                     end do

                  jm = im
  140             jm = in

                  km = jm - im + 1
                  lm = ( km - 1 ) / ipstep + 1

               do mmm = 1, lm

                  n1 = ipstep * ( mmm - 1 ) + 1
                  n2 = min( ipstep * mmm, km )
                  n3 = n1 + im - 1
                  n4 = n2 + im - 1

                  IF(il .eq. 0) then
                   write(iot,'(/1x,i4,"-",a2,
     &                       " isotope production",
     &                       " # n3 n4 = ",i3,1x,i3)')
     &                 iz, elmnt(iz), n3, n4

                  ELSEIF(il .eq. 1) then
                   write(iot,'(/1x,i4,"-",a2,
     &                       " 1st metastable isotope production",
     &                       " # n3 n4 = ",i3,1x,i3)')
     &                 iz, elmnt(iz), n3, n4

                  ELSEIF(il .eq. 2) then
                   write(iot,'(/1x,i4,"-",a2,
     &                       " 2nd metastable isotope production",
     &                       " # n3 n4 = ",i3,1x,i3)')
     &                 iz, elmnt(iz), n3, n4
                  ENDIF

                  IF(il .eq. 0) then
                     write(iot,'(" reg.",12i11)')
     &                      ( im - 1 + iz + n, n = n1, n2 )
                  ELSEIF(il .eq. 1) then
                     write(iot,'(" reg.",i10,"m",i10,"m",i10,"m"
     &               ,i10,"m",i10,"m",i10,"m",i10,"m",i10,"m"
     &               ,i10,"m",i10,"m",i10,"m",i10,"m")')
     &                      ( im - 1 + iz + n, n = n1, n2 )
                  ELSEIF(il .eq. 2) then
                     write(iot,'(" reg.",i10,"n",i10,"n",i10,"n"
     &               ,i10,"n",i10,"n",i10,"n",i10,"n",i10,"n"
     &               ,i10,"n",i10,"n",i10,"n",i10,"n")')
     &                      ( im - 1 + iz + n, n = n1, n2 )
                  ENDIF

                  do ir = 1, nr

c  *** Changed by T.Sato 2013/10/9, i5 -> i7
                   if( itnzn(m).eq.0 ) then ! frtati 2022/02/18
                     write(iot,'(i7,1p12e11.3)')
     &                     lr(ir), ( tr(ir,iz,i,il,ioe), i = n3, n4 )
                   else
                     write(iot,'(i7,1p12e11.3)') lr(ir),
     &               ( tr(ir,igetiznmp(iz,i,il,m),1,0,ioe),i = n3, n4 )
                   end if

                  end do

                  if( itout(m) .ne. 0 ) then

                  IF(il .eq. 0) then
                    write(iou,'(/1x,i4,"-",a2,
     &              " isotope production : ERROR(%)")')
     &              iz, elmnt(iz)

                  ELSEIF(il .eq. 1) then
                    write(iou,'(/1x,i4,"-",a2,
     &              " 1st metastable isotope production : ERROR(%)")')
     &              iz, elmnt(iz)

                  ELSEIF(il .eq. 2) then
                    write(iou,'(/1x,i4,"-",a2,
     &              " 2nd metastable isotope production : ERROR(%)")')
     &              iz, elmnt(iz)
                  ENDIF

                  IF(il .eq. 0) then
                     write(iou,'(" reg.",12i11)')
     &                      ( im - 1 + iz + n, n = n1, n2 )
                  ELSEIF(il .eq. 1) then
                     write(iou,'(" reg.",i10,"m",i10,"m",i10,"m"
     &               ,i10,"m",i10,"m",i10,"m",i10,"m",i10,"m"
     &               ,i10,"m",i10,"m",i10,"m",i10,"m")')
     &                      ( im - 1 + iz + n, n = n1, n2 )
                  ELSEIF(il .eq. 2) then
                     write(iou,'(" reg.",i10,"n",i10,"n",i10,"n"
     &               ,i10,"n",i10,"n",i10,"n",i10,"n",i10,"n"
     &               ,i10,"n",i10,"n",i10,"n",i10,"n")')
     &                      ( im - 1 + iz + n, n = n1, n2 )
                  ENDIF

                     do ir = 1, nr

                      if( itnzn(m).eq.0 ) then ! frtati 2022/02/18
                        write(iou,'(i7,1p12e11.3)')
     &                    lr(ir), ( tr(ir,iz,i,il,2)*100.0, i = n3, n4 )
                      else
                        write(iou,'(i7,1p12e11.3)') lr(ir),
     &                  ( tr(ir,igetiznmp(iz,i,il,m),1,0,2)*100.0,
     &                    i = n3, n4 )
                      end if

                     end do

                  end if

               end do

  170       continue

*-----------------------------------------------------------------------
            else
             write(iot,'(a)')
             write(iot,'("# num nucleusID yield r.err")')
             if(itout(m).ne.0)then
              write(iou,'(a)')
              write(iou,'("# num nucleusID yield r.err")')
             endif
             cfmt='(i#,x,i7,1p2e11.3)'
             do ir=1,nr
!NS 2021.04 change for NVIDIA TOOL KIT COMPILL ERR
              write(cir,"(I12)")ir
!NS 2021.04 end change for NVIDIA TOOL KIT COMPILL ERR
              nir=len_trim(adjustl(cir))
              write(cfmt(3:3),'(i1)')nir
              do iz=1,maxpt
               if(nn.gt.0) then
                do i = 1,nn
                 if(nt(i)/1000.eq.iz)exit
                enddo
                if(i.gt.nn)cycle
               endif
               do in=1,maxnt
                do il=0,2 ! frtati 2022/05/02
                 if( itnzn(m).eq.0 ) then ! frtati 2022/02/18
                  if(tr(ir,iz,in,il,ioe).ne.0d0)then
                   write(iot,cfmt)ir,iz*10000+(in+iz)*10+il,
     &                  (tr(ir,iz,in,il,ioee),ioee=1,2)
                   if(itout(m).ne.0)then
                    write(iou,cfmt)ir,iz*10000+(in+iz)*10+il,
     &                   tr(ir,iz,in,il,2)
                   endif
                  endif
               else
                  if(tr(ir,igetiznmp(iz,in,il,m),1,0,ioe).ne.0d0)then
                   write(iot,cfmt)ir,iz*10000+(in+iz)*10+il,
     &             (tr(ir,igetiznmp(iz,in,il,m),1,0,ioee),ioee=1,2)
                   if(itout(m).ne.0)then
                    write(iou,cfmt)ir,iz*10000+(in+iz)*10+il,
     &                   tr(ir,igetiznmp(iz,in,il,m),1,0,2)
                   endif
                  endif
                 endif
                enddo
               enddo
              enddo
             enddo
             write(iot,cfmt)0,0,0.0d0,0.0d0
             if(itout(m).ne.0)write(iou,cfmt)0,0,0.0d0,0.0d0
            endif
*-----------------------------------------------------------------------

            if( itout(m) .ne. 0 ) close(iou)

*-----------------------------------------------------------------------
*        reg axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 2 ) then

               inum = 0

                  if( nn .eq. 0 ) then

                     nc = 1

                  else

                     nc = nn

                  end if

*-----------------------------------------------------------------------

            do ic = 1, nc

               inum = inum + 1

*-----------------------------------------------------------------------

               if( nn .gt. 0 ) then

                  iz = nt(ic) / 1000

                  ia = nt(ic) - iz * 1000

                  if( ia .gt. 0 ) then

                     in = ia - iz

                     call chname(idum,ia,iz,chau)

                  end if

               end if

*-----------------------------------------------------------------------

                  write(iot,'(/"#",78("-"))')

               if( inum .eq. 1 ) then
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

               if( rtfac(m) .lt. 0.d0 )
     &           write(iot,
     &           '( "#restart maximum for normalization:",1p1e14.7)')
     &           facmax(m)

            if( nn .eq. 0 ) then

               write(iot,'("#   no. =",i3,3x,"Z = all")')
     &                     inum

            else if( ia .eq. 0 ) then

               write(iot,'("#   no. =",i3,3x,"Z = ",i3," : ",a3)')
     &                     inum, iz, elmnt(iz)

            else

               write(iot,'("#   no. =",i3,3x,a8)')
     &                     inum, chau

            end if

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: Serial Num. of Region")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: Number ",a15)') hsunit(itunt(m))


               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

                  write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')

               if( itanl(m) .gt. 0 ) then

                  write(iot,'( "p: ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               end if

               if( itsans(m) .gt. 0 ) then

                  call write_sangel(iot,m,0)

               end if

c S.H. added IF statement below for epsout=2 (2017.2.18)
              if ( iteps(m) .ne. 2 ) then

               write(iot,'( "h:   x      n",5x,"n",10x,
     &                    "y,l3        n")')

              else

               write(iot,'( "h:   x      n",5x,"n",10x,
     &                    "y1,l3       ny2 dy1=[y1*y2]")')

              end if

               write(iot,'( "#  num    reg     volume  ",
     &                    "   number      r.err")')

*-----------------------------------------------------------------------

                     seka = 0.0
                     sera = 0.0
                     seva = 0.0
                     voll = 0.0

               do ir = 1, nr

                     ireg = lr(ir)

                     sek = 0.0
                     ser = 0.0

                     if( itunt(m) .eq. 1 ) then

                        vm = 1.0

                     else if( itunt(m) .eq. 2 ) then

                        vm = vl(ir)

                     end if

                        voll = voll + vl(ir)

*-----------------------------------------------------------------------

                  if( nn .eq. 0 ) then

                     do kn = 1, mn
                     do kz = 1, mz
                     do il = 0, mm ! frtati 2022/02/18

                        lz  = ikzz( kz, kn )
                        ln  = iknn( kz, kn )

                        vn  = vm * tr(ir,lz,ln,il,1)
                        sek = sek + vn
                        ser = ser + ( vn * tr(ir,lz,ln,il,2) )**2

                     end do
                     end do
                     end do

                  else if( ia .eq. 0 ) then

                     do kn = 1, mn
                     do il = 0, mm ! frtati 2022/02/18

                        lz  = ikzz( iz, kn )
                        ln  = iknn( iz, kn )

                        vn  = vm * tr(ir,lz,ln,il,1)
                        sek = sek + vn
                        ser = ser + ( vn * tr(ir,lz,ln,il,2) )**2

                     end do
                     end do

                  else

                        lz  = ikzz( iz, in )
                        ln  = iknn( iz, in )

                     do il = 0, mm ! frtati 2022/02/18
                        vn  = vm * tr(ir,lz,ln,il,1)
                        sek = sek + vn
                        ser = ser
     &                      + ( vn * tr(ir,lz,ln,il,2) )**2
                     end do

                  end if

*-----------------------------------------------------------------------

                        seka = seka + sek
                        sera = sera + ser
                        seva = seva + vm

                     if( sek .ne. 0.0 ) then ! .gt. -> .ne. by T.Sato 2024/05/02 for negative sek

                        ser = sqrt( ser ) / sek
                        sek = sek / vm

                     end if

                  write(iot,'(i5,1x,i7,1pe13.4,1pe13.4,0pf8.4)')
     &                  ir, ireg, vl(ir), sek, ser

               end do

                     if( seka .gt. 0.0 ) then

                        sera = sqrt( sera ) / seka

                        if( itunt(m) .eq. 2 )
     &                  seka = seka / seva

                     end if

               write(iot,'(/"#   sum over ",1pe13.4,1pe13.4,0pf8.4)')
     &                     voll, seka, sera

            if( itunt(m) .eq. 2 ) then

               write(iot,'(
     &         "#   vol sum = ",1p1e13.4,"  [cm^3]")') seva

            end if

               if( nn .eq. 0 ) then

                  write(iot,'(/a1,"no. =",i3,
     &                        ",    Z  =  all", a1)')
     &                        cha, inum, cha

               else if( ia .eq. 0 ) then

                  write(iot,'(/a1,"no. =",i3,
     &                        ",    Z  =  ",i3, a1)')
     &                        cha, inum, iz, cha

               else

                  write(iot,'(/a1,"no. =",i3,
     &                        ",    ",a8, a1)')
     &                        cha, inum, chau, cha

               end if

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

               if( nn .eq. 0 ) then

                  write(iot,'("wt: s(0.7)",/
     &                        a1,"vspace{-3}"/
     &                        "      Z &=&   all"/
     &                        "e:")')
     &                        yen

               else if( ia .eq. 0 ) then

                  write(iot,'("wt: s(0.7)",/
     &                        a1,"vspace{-3}"/
     &                        "      Z &=& ",i3,"   :  ",a3/
     &                        "e:")')
     &                        yen, iz, elmnt(iz)

               else

                  write(iot,'("wt: s(0.7)",/
     &                        a1,"vspace{-3}"/
     &                        "  nucl. &=&   ",a8/
     &                        "e:")')
     &                        yen, chau

               end if

*-----------------------------------------------------------------------

            end do

*-----------------------------------------------------------------------
*        mass axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 1 ) then

               inum = 0

                  if( nn .eq. 0 ) then

                     nc = 1

                  else

                     nc = nn

                  end if

*-----------------------------------------------------------------------

            do ir = 1, nr

            do 260 ic = 1, nc

                  if( nn .gt. 0 ) iz = nt(ic) / 1000

                  do i = 1, maxnt + maxpt

                     tm(i,1) = 0.0d+0
                     tm(i,2) = 0.0d+0

                  end do

                  if( itunt(m) .eq. 1 ) then

                     vm = 1.0

                  else if( itunt(m) .eq. 2 ) then

                     vm = vl(ir)

                  end if

*-----------------------------------------------------------------------

                  if( nn .eq. 0 .and. itfln(m) .eq. 1 ) then

                     do ln = 1, mn
                     do il = 0, mm ! frtati 2022/02/18

                        lz = 1
                        im = ln

                        vn = vm * tr(ir,lz,ln,il,1)
                        tm(im,1) = tm(im,1) + vn
                        tm(im,2) = tm(im,2)
     &                           + ( vn * tr(ir,lz,ln,il,2) )**2

                     end do
                     end do

                  else if( nn .eq. 0 .and. itfln(m) .gt. 1 ) then

                     do kn = 1, mn
                     do kz = 1, mz
                     do il = 0, mm ! frtati 2022/02/18

                        lz = ikzz( kz, kn )
                        ln = iknn( kz, kn )

                        im = kz + kn

                        vn = vm * tr(ir,lz,ln,il,1)
                        tm(im,1) = tm(im,1) + vn
                        tm(im,2) = tm(im,2)
     &                           + ( vn * tr(ir,lz,ln,il,2) )**2

                     end do
                     end do
                     end do

                  else

                     do kn = 1, mn
                     do il = 0, mm ! frtati 2022/02/18

                        lz = ikzz( iz, kn )
                        ln = iknn( iz, kn )

                        im = iz + kn

                        vn = vm * tr(ir,lz,ln,il,1)
                        tm(im,1) = tm(im,1) + vn
                        tm(im,2) = tm(im,2)
     &                           + ( vn * tr(ir,lz,ln,il,2) )**2

                     end do
                     end do

                  end if

*-----------------------------------------------------------------------

               do i = 1, maxnt + maxpt

                  if( tm(i,1) .gt. 0.0d0 ) then

                     tm(i,2) = sqrt( tm(i,2) ) / tm(i,1)
                     tm(i,1) = tm(i,1) / vm

                  end if

               end do

               do i = 1, maxnt + maxpt

                  if( tm(i,1) .gt. 0.0d0 ) goto 200

               end do

               goto 260

  200          im = i

               do i = maxnt + maxpt, im + 1, -1

                  if( tm(i,1) .gt. 0.0d0 ) goto 210

               end do

               jm = im
  210          jm = i

               im = max( 1, im - 2 )
               jm = jm + 2

*-----------------------------------------------------------------------

               ireg = lr(ir)

                  write(iot,'(/"#",78("-"))')

               inum = inum + 1

               if( inum .eq. 1 ) then
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

               if( rtfac(m) .lt. 0.d0 )
     &           write(iot,
     &           '( "#restart maximum for normalization:",1p1e14.7)')
     &           facmax(m)

               if( nn .eq. 0 ) then

                  write(iot,'("#   no. =",i3,3x,
     &            "reg  =",i5,3x,"Z = all")')
     &                        inum, ireg

               else

                  write(iot,'("#   no. =",i3,3x,
     &            "reg  =",i5,3x,"Z = ",i3," : ",a3)')
     &                        inum, ireg, iz, elmnt(iz)

               end if

               write(iot,'("# im jm = ",i3,1x,i3)') im, jm

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: Mass")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: Number ",a15)') hsunit(itunt(m))

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

                  write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')

               if( itanl(m) .gt. 0 ) then

                  write(iot,'( "p: ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               end if

               if( itsans(m) .gt. 0 ) then

                  call write_sangel(iot,m,0)

               end if

c S.H. added IF statement below for epsout=2 (2017.2.18)
              if ( iteps(m) .ne. 2 ) then

               write(iot,'( "h: x-0.5    y,hl0       n")')

              else

               write(iot,'( "h: x-0.5    y1,hl0      ",
     &                      "ny2 dy1=[y1*y2]")')

              end if

               write(iot,'( "#  mass     number      r.err")')

*-----------------------------------------------------------------------

               if( im .eq. 1 ) then

                  write(iot,'(3x,f4.1,2x,1pe13.4,0pf8.4)')
     &                                   0.5, 0.0, 0.0

               end if

                  seka = 0.0
                  sera = 0.0

               do i = im, jm

                  vn  = vm * tm(i,1)
                  seka = seka + vn
                  sera = sera + ( vn * tm(i,2) )**2

                  write(iot,'(3x,i4,2x,1pe13.4,0pf8.4)')
     &                           i, tm(i,1), tm(i,2)

               end do

                  if( seka .gt. 0.0 ) then

                     sera = sqrt( sera ) / seka
                     seka = seka / vm

                  end if

               write(iot,'(/"#   sum  ",1pe13.4,0pf8.4)')
     &               seka, sera

               if( itunt(m) .eq. 2 ) then

                  write(iot,'(
     &            "#   vol sum = ",1p1e13.4,"  [cm^3]")') vm

               end if

               if( nn .eq. 0 ) then

                  write(iot,'(/a1,"no. =",i3,
     &                            ",    reg  =",i7,
     &                            ",    Z  =  all", a1)')
     &                        cha, inum, ireg, cha

               else

                  write(iot,'(/a1,"no. =",i3,
     &                            ",    reg  =",i7,
     &                            ",    Z  =  ",i3, a1)')
     &                        cha, inum, ireg, iz, cha

               end if

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

               if( itunt(m) .eq. 1 ) then

                  if( nn .eq. 0 ) then

                     write(iot,'("wt: s(0.7)",/
     &                           a1,"vspace{-3}"/
     &                           "  mass distribution"/
     &                           "  in region mesh"/
     &                           "    Z &=&   all"/
     &                           "e:")')
     &                           yen

                  else

                     write(iot,'("wt: s(0.7)",/
     &                           a1,"vspace{-3}"/
     &                           "  isotope distribution"/
     &                           "  in region mesh"/
     &                           "    Z &=& ",i3,"   :  ",a3/
     &                           "e:")')
     &                           yen, iz, elmnt(iz)

                  end if

               else if( itunt(m) .eq. 2 ) then

                  if( nn .eq. 0 ) then

                     write(iot,'("wt: s(0.7)",/
     &                           a1,"vspace{-3}"/
     &                           "  mass distribution"/
     &                           "  in region mesh"/
     &                           "  vol &=&",1pe13.4," [cm^3]"/
     &                           "    Z &=&   all"/
     &                           "e:")')
     &                           yen, vl(ir)

                  else

                     write(iot,'("wt: s(0.7)",/
     &                           a1,"vspace{-3}"/
     &                           "  isotope distribution"/
     &                           "  in region mesh"/
     &                           "  vol &=&",1pe13.4," [cm^3]"/
     &                           "    Z &=& ",i3,"   :  ",a3/
     &                           "e:")')
     &                           yen, vl(ir),
     &                           iz, elmnt(iz)

                  end if

               end if

  260       continue

            end do

*-----------------------------------------------------------------------
*        charge axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 7 ) then

               inum = 0

            do 370 ir = 1, nr

                  do i = 1, maxpt

                     tm(i,1) = 0.0d+0
                     tm(i,2) = 0.0d+0

                  end do

                  if( itunt(m) .eq. 1 ) then

                     vm = 1.0

                  else if( itunt(m) .eq. 2 ) then

                     vm = vl(ir)

                  end if

*-----------------------------------------------------------------------

                  do kz = 1, mz
                  do kn = 1, mn
                  do il = 0, mm ! frtati 2022/02/18

                        iz  = ikzz( kz, kn )
                        in  = iknn( kz, kn )

                        vn  = vm * tr(ir,iz,in,il,1)
                        tm(iz,1) = tm(iz,1) + vn
                        tm(iz,2) = tm(iz,2)
     &                           + ( vn * tr(ir,iz,in,il,2) )**2

                  end do
                  end do
                  end do

*-----------------------------------------------------------------------

               do i = 1, maxpt

                  if( tm(i,1) .gt. 0.0d0 ) then

                     tm(i,2) = sqrt( tm(i,2) ) / tm(i,1)
                     tm(i,1) = tm(i,1) / vm

                  end if

               end do

               do i = 1, maxpt

                  if( tm(i,1) .gt. 0.0d0 ) goto 300

               end do

               goto 370

  300          im = i

               do i = maxpt, im + 1, -1

                  if( tm(i,1) .gt. 0.0d0 ) goto 310

               end do

               jm = im
  310          jm = i

               im = max( 1, im - 2 )
               jm = jm + 2

*-----------------------------------------------------------------------

               ireg = lr(ir)

                  write(iot,'(/"#",78("-"))')

               inum = inum + 1

               if( inum .eq. 1 ) then
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

               if( rtfac(m) .lt. 0.d0 )
     &           write(iot,
     &           '( "#restart maximum for normalization:",1p1e14.7)')
     &           facmax(m)

               write(iot,'("#   no. =",i3,3x,"reg =",i7)')
     &                     inum, ireg

               write(iot,'("# im jm = ",i3,1x,i3)') im, jm

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: Charge")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: Number ",a15)') hsunit(itunt(m))

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

                  write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')

               if( itanl(m) .gt. 0 ) then

                  write(iot,'( "p: ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               end if

               if( itsans(m) .gt. 0 ) then

                  call write_sangel(iot,m,0)

               end if

c S.H. added IF statement below for epsout=2 (2017.2.18)
              if ( iteps(m) .ne. 2 ) then

               write(iot,'( "h: x-0.5    y,hl0       n")')

              else

               write(iot,'( "h: x-0.5    y1,hl0      ",
     &                      "ny2 dy1=[y1*y2]")')

              end if

               write(iot,'( "#  charge   number      r.err")')

*-----------------------------------------------------------------------

               if( im .eq. 1 ) then

                  write(iot,'(3x,f4.1,2x,1pe13.4,0pf8.4)')
     &                                   0.5, 0.0, 0.0

               end if

                  seka = 0.0
                  sera = 0.0

               do i = im, jm

                  vn  = vm * tm(i,1)
                  seka = seka + vn
                  sera = sera + ( vn * tm(i,2) )**2

                  write(iot,'(3x,i4,2x,1pe13.4,0pf8.4)')
     &                           i, tm(i,1), tm(i,2)

               end do

                  if( seka .gt. 0.0 ) then

                     sera = sqrt( sera ) / seka
                     seka = seka / vm

                  end if

               write(iot,'(/"#   sum  ",1pe13.4,0pf8.4)')
     &               seka, sera

               if( itunt(m) .eq. 2 ) then

                  write(iot,'(
     &            "#   vol sum = ",1p1e13.4,"  [cm^3]")') vm

               end if

                  write(iot,'(/a1,"no. =",i3,
     &                            ",    reg  =",i7,a1)')
     &                        cha, inum, ireg, cha

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

               if( itunt(m) .eq. 1 ) then

                     write(iot,'("wt: s(0.7)",/
     &                           a1,"vspace{-3}"/
     &                           "  charge distribution"/
     &                           "  in region mesh"/
     &                           "e:")')
     &                           yen

               else if( itunt(m) .eq. 2 ) then

                     write(iot,'("wt: s(0.7)",/
     &                           a1,"vspace{-3}"/
     &                           "  charge distribution"/
     &                           "  in region mesh"/
     &                           "  vol &=&",1pe13.4," [cm^3]"/
     &                           "e:")')
     &                           yen, vl(ir)

               end if

  370       continue

*-----------------------------------------------------------------------
*        chart axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 8 ) then

               inum = 0

            do 470 ir = 1, nr

            do 460 il = 0, 2 ! loop for isomeric level
*-----------------------------------------------------------------------

                  itmax  = 0
                  icmax  = 0
                  inmax  = 0

               do iz = 1, maxpt
               do in = 1, maxnt

                  if( itnzn(m).eq.0 ) then
                    t0_tr = tr(ir,iz,in,il,1)
                  else
                    t0_tr = tr(ir,igetiznmp(iz,in,il,m),1,0,1)
                  end if
                  if( t0_tr .gt. 0.d0 ) then

                     if( iz + in .gt. itmax ) itmax = iz + in
                     if( iz      .gt. icmax ) icmax = iz
                     if(      in .gt. inmax ) inmax =      in

                  end if

               end do
               end do

                     dxmax = dble(inmax+2)
                     dymax = dble(icmax+2)
                     dform = dymax / dxmax

                     inmag = inmax

                  if( dform .gt. 0.8 ) then

                     inmag = nint( dymax / 0.8 )
                     dxmax = dble( nint( dymax / 0.8 ) )
                     dform = dymax / dxmax

                  end if

               if( itmax .eq. 0 ) goto 470

*-----------------------------------------------------------------------

               ireg = lr(ir)

                  write(iot,'(/"#",78("-"))')

               inum = inum + 1

               if( inum .eq. 1 .and. il .eq. 0 ) then
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

               if( rtfac(m) .lt. 0.d0 )
     &           write(iot,
     &           '( "#restart maximum for normalization:",1p1e14.7)')
     &           facmax(m)

               write(iot,'("#   no. =",i3,3x,
     &         "ir  =",i3,3x,"il  =",i3)' )
     &                     inum, ir, il   ! tally data number, region number, isomeric level

               write(iot,'("# icmax inmax = ",i3,1x,i3)')
     &            icmax, inmax

*-----------------------------------------------------------------------

            if( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) then

                  write(iot,'(/a1,"no. =",i3,
     &                            ",    reg  =",i7,a1)')
     &                        cha, inum, ireg, cha

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

             if(il .eq. 0) then    ! specify isomeric level at the bottom of figure
               write(iot,'("msdc: {",a1,"huge Ground state}")') yen
             elseif(il .eq. 1) then
               write(iot,'("msdc: {",a1,"huge 1st isomer}")') yen
             else
               write(iot,'("msdc: {",a1,"huge 2nd isomer}")') yen
             endif

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: N Neutron Number")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: Z Proton Number")')

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

*-----------------------------------------------------------------------

               if( inum .eq. 1 ) then

                  form  = dform
                  xfac  = 1.1
                  afac  = 0.6
                  izlog = 1
                  inocm = 1
                  inolg = 1

                  if( itanl(m) .gt. 0 )
     &            call anset(itang(m),itanl(m),
     &                       form,xfac,afac,cmin,cmax,izlog,inocm,inolg)

                  if( form .le. 1.0 ) then

                     scal = form**0.35
                     xfac = xfac / form**0.5
                     xorg = 0.0
                     yorg = min( 20.d0, ( 1.0 / form - 1.0 ) / 2.5 )

                  else

                     scal = 1.0 / form**0.41
                     xfac = xfac / form**0.5
                     xorg = min( 20.d0, ( form - 1.0 ) / 1.5)
                     yorg = 0.0

                  end if

               end if

                  write(iot,'( "set: c1[",f6.3,
     &                            "] c2[",f6.3,
     &                            "] c5[",f6.3,
     &                            "] c6[",f6.3,
     &                            "] c7[",f6.3,
     &                            "] c8[",f6.3,"]")')
     &                    form, xfac, afac, scal, xorg, yorg
                  write(iot,'( "p: h2fs form[c1] xfac[c2]",
     &                 " afac[c5] scal[c6] xorg[c7] yorg[c8] nosp")')

               if( ( ( ittwo(m) .ge. 2 .and. ittwo(m) .le. 3 ) .or.
     &               ( ittwo(m) .ge. 6 .and. ittwo(m) .le. 7 ) ) .and.
     &             cmin .gt. 0.0 .and. cmax .gt. cmin ) then

                if (ioe .eq. 1 ) then
                  write(iot,'( "set: c3[",1p1e15.6,"] c4[",
     &                  1p1e15.6,"]")') cmin, cmax
                else
                  write(iot,'( "set: c3[1.0e-4] c4[1.0]")')
                end if
                  write(iot,'( "p: cmin[c3] cmax[c4]")')
                  write(iot,'( "p: dmin(1e-31)")')

                  if( ittwo(m) .eq. 3 .or. ittwo(m) .eq. 7 )
     &                 write(iot,'( "p: zlog")')

               end if

               write(iot,'( "p: xmin(0) xmax(",f5.1,")")') dxmax
               write(iot,'( "p: ymin(0) ymax(",f5.1,")")') dymax

               if( itanl(m) .gt. 0 .and. ioe .eq. 1 ) then

                  write(iot,'( "p: ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               else if( itanl(m) .gt. 0 .and. ioe .ne. 1 ) then

                   call terrang(iot,m)

               end if

               if( itsans(m) .gt. 0 ) then

                  call write_sangel(iot,m,0)

               end if

            end if

*-----------------------------------------------------------------------

            if( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) then

               write(iot,'(/"hc: y = ",i3," to 1 by -1 ;",
     &                         " x = 1 to ",i3," by 1 ;")')
     &                        icmax+2, inmax+2

               do i = icmax+2, 1, -1
                if( itnzn(m).eq.0 ) then ! frtati 2022/02/18
                  write(iot,'(1p10e11.3)')
     &                 ( tr(ir,i,l,il,ioe), l = 1, inmax+2 )
                else
                  write(iot,'(1p10e11.3)')
     &          ( tr(ir,igetiznmp(i,l,il,m),1,0,ioe), l = 1, inmax+2 )
                end if

               end do

*-----------------------------------------------------------------------

            else if( ittwo(m) .eq. 4 ) then

               write(iot,'(/"#   Z    N  Mass  ",
     &                      "  number    r.err")')

               do i = 1, icmax
               do j = 1, inmax

                if( itnzn(m).eq.0 ) then ! frtati 2022/02/18
                  if( tr(ir,i,j,il,1) .gt. 0.d0 ) then

                     write(iot,'(3i5,1pe13.4,0pf8.4)')
     &               i, j, i+j,
     &               tr(ir,i,j,il,1), tr(ir,i,j,il,2)

                  end if
                else
                  if( tr(ir,igetiznmp(i,j,il,m),1,0,1) .gt. 0.d0 ) then
                     write(iot,'(3i5,1pe13.4,0pf8.4)')
     &               i, j, i+j,
     &               tr(ir,igetiznmp(i,j,il,m),1,0,1),
     &               tr(ir,igetiznmp(i,j,il,m),1,0,2)
                  end if
                end if

               end do
               end do

*-----------------------------------------------------------------------

            else if( ittwo(m) .eq. 5 ) then

               write(iot,'(/"#   Z = 1 to ",i3/
     &                      "#   N = 1 to ",i3)')
     &                        icmax+2, inmax+2

                  write(iot,'(/4x,a3,4x,1p1000e11.3)')
     &            'Z/N', ( dble( in ), in = 1, inmax+2 )

               do i = icmax+2, 1, -1
                if( itnzn(m).eq.0 ) then ! frtati 2022/02/18
                  write(iot,'(1p1000e11.3)')
     &            dble( i ),
     &            ( tr(ir,i,l,il,ioe), l = 1, inmax+2 )
                else
                  write(iot,'(1p1000e11.3)')
     &            dble( i ),
     &            ( tr(ir,igetiznmp(i,l,il,m),1,0,ioe),l = 1,inmax+2 )
                end if

               end do

            end if

*-----------------------------------------------------------------------

            if( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) then

               if( itout(m) .ne. 0 ) then

                  call wmgcstb(iot,icmax,inmag)

               end if

            end if

*-----------------------------------------------------------------------

      if( inolg .eq. 1 .and. ( ittwo(m) .eq. 3 .or. ittwo(m) .eq. 7 )
     &    .and. cmin .gt. 0.0 .and. cmax .gt. cmin .and.
     &    igsh .eq. 0 ) then

      write(iot,'(/"#",78("-"))')
      write(iot,'(
     &"z: xorg(1.03)"/
     &"p: xfac[c2*0.03] form[c1/0.03] notc noxt noyt dmin(1)"/
     &"p: ymin(0) ymax(1) xmin(0) xmax(1) cmin(1) cmax(100) nosp"/
     &"hc: y= 0.005 to 0.995 by 0.01 ; x= 0.5 to 0.5 by 1 ;"/
     &" 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20"/
     &"21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40"/
     &"41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60"/
     &"61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80"/
     &"81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100"/
     &"z: xorg(0.0)"/
     &"p: xfac[c2*0.03] form[c1/0.03] ",
     &"nosp afac[c5*0.7] ymin[c3] ymax[c4]"/
     &"p: noxt noxn ytxt(-1) ynum(-1) xltd(-1) xstd(-1) itic(1) ylog"
     &)')

      if ( ioe .eq. 2 ) then

          write(iot,'("y: Relative Error")')

      else if( itazl(m) .eq. 0 ) then

         write(iot,'("y: Number ",a15)') hsunit(itunt(m))

      else

         write(iot,'("y: ",200a1)') (itazt(m)(i:i),i=1,itazl(m))

      end if

      if( izlog .eq. 0 ) write(iot,'( "p: ylin")')

      write(iot,'( "#",78("-"))')

      end if

*-----------------------------------------------------------------------

            if( inocm .eq. 1 ) then

               if( itunt(m) .eq. 1 ) then

                     write(iot,'(/"wt: s(0.7)",/
     &                           a1,"vspace{-3}"/
     &                           "  yield distribution"/
     &                           "  in region mesh"/
     &                           "e:")')
     &                           yen

               else if( itunt(m) .eq. 2 ) then

                     write(iot,'(/"wt: s(0.7)",/
     &                           a1,"vspace{-3}"/
     &                           "  yield distribution"/
     &                           "  in region mesh"/
     &                           "  vol &=&",1pe13.4," [cm^3]"/
     &                           "e:")')
     &                           yen, vl(ir)

               end if

            end if

*-----------------------------------------------------------------------

      if( inolg .eq. 1 .and. ( ittwo(m) .eq. 3 .or. ittwo(m) .eq. 7 )
     &    .and. cmin .gt. 0.0 .and. cmax .gt. cmin .and.
     &    igsh .eq. 0 ) then

         write(iot,'(/"#",78("-"))')
         write(iot,'("z: xorg[-1.03/0.05]"/
     &               "p: form[c1] xfac[c2] afac[c5] nosp notf")')

               write(iot,'( "p: xmin(0) xmax(",f5.1,")")') dxmax
               write(iot,'( "p: ymin(0) ymax(",f5.1,")")') dymax

      end if

*-----------------------------------------------------------------------

  460       continue

  470       continue

*-----------------------------------------------------------------------
*        xy axis (matrix)
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 9 ) then

*-----------------------------------------------------------------------

               inum = 0

                  if( nn .eq. 0 .or. igsh .ne. 0  ) then

                     nc = 1

                  else

                     nc = nn

                  end if

*-----------------------------------------------------------------------

            do jz = 1, nz
            do ic = 1, nc

               zval = ( zm(jz) + zm(jz+1) ) / 2.0d0

               inum = inum + 1

*-----------------------------------------------------------------------

               if( nn .gt. 0 ) then

                  iz = nt(ic) / 1000

                  ia = nt(ic) - iz * 1000

                  if( ia .gt. 0 ) then

                     in = ia - iz

                     call chname(idum,ia,iz,chau)

                  end if

               end if

*-----------------------------------------------------------------------

                  write(iot,'(/"#",78("-"))')

               if( inum .eq. 1 ) then
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

               if( rtfac(m) .lt. 0.d0 )
     &           write(iot,
     &           '( "#restart maximum for normalization:",1p1e14.7)')
     &           facmax(m)

            if( nn .eq. 0 ) then

               write(iot,'("#   no. =",i3,3x,
     &         "iz  =",i3,3x,"Z = all")')
     &                     inum, jz

            else if( ia .eq. 0 ) then

               write(iot,'("#   no. =",i3,3x,
     &         "iz  =",i3,3x,"Z = ",i3," : ",a3)')
     &                     inum, jz, iz, elmnt(iz)

            else

               write(iot,'("#   no. =",i3,3x,
     &         "iz  =",i3,3x,a8)')
     &                     inum, jz, chau

            end if

*-----------------------------------------------------------------------

               if( nn .eq. 0 ) then

                  write(iot,'(/a1,"no. =",i3,
     &                        ",    iz =",i3,
     &                        ",    Z  =  all", a1)')
     &                        cha, inum, jz, cha

               else if( ia .eq. 0 ) then

                  write(iot,'(/a1,"no. =",i3,
     &                        ",    iz =",i3,
     &                        ",    Z  =  ",i3, a1)')
     &                        cha, inum, jz, iz, cha

               else

                  write(iot,'(/a1,"no. =",i3,
     &                        ",    iz =",i3,
     &                        ",    ",a8, a1)')
     &                        cha, inum, jz, chau, cha

               end if

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: x [cm]")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: y [cm]")')

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

*-----------------------------------------------------------------------

         if( igsh .eq. 0 ) then

               do i = 1, nr

                     sek = 0.0

                  if( nn .eq. 0 ) then

                     do kn = 1, mn
                     do kz = 1, mz
                     do il = 0, mm ! frtati 2022/02/18

                        lz  = ikzz( kz, kn )
                        ln  = iknn( kz, kn )

                        sek = sek + tr(i,lz,ln,il,1)

                     end do
                     end do
                     end do

                  else if( ia .eq. 0 ) then

                     do kn = 1, mn
                     do il = 0, mm ! frtati 2022/02/18

                        lz  = ikzz( iz, kn )
                        ln  = iknn( iz, kn )

                        sek = sek + tr(i,lz,ln,il,1)

                     end do
                     end do

                  else

                        lz  = ikzz( iz, in )
                        ln  = iknn( iz, in )

                     do il = 0, mm ! frtati 2022/02/18
                        sek = sek + tr(i,lz,ln,il,1)
                     end do

                  end if

                  val(i) = sek

               end do

                  izlog = 1
                  if( itanl(m) .gt. 0 )
     &            call anset(itang(m),itanl(m),
     &                       form,xfac,afac,cmin,cmax,izlog,inocm,inolg)

            if( izlog .eq. 0 ) then

               do i = 1, nr

                  if( val(i) .gt. dnon .and. cmin .ne. cmax ) then

                     val(i) = ( val(i) - cmin ) / ( cmax - cmin )
                     val(i) = max(0.0d0,val(i))
                     val(i) = min(1.0d0,val(i))

                  else

                     val(i) = -1.0

                  end if

               end do

            else

               do i = 1, nr

                  if( val(i) .gt. dnon .and. cmin .ne. cmax ) then

                     val(i) = log10(val(i)/cmin) / log10(cmax/cmin)
                     val(i) = max(0.0d0,val(i))
                     val(i) = min(1.0d0,val(i))

                  else

                     val(i) = -1.0d0

                  end if

               end do

            end if

         end if

*-----------------------------------------------------------------------

               if( inum .eq. 1 ) then

                  form  = ( ym(ny+1) - ym(1) )
     &                  / ( xm(nx+1) - xm(1) )
                  xfac  = 0.9
                  afac  = 0.8
                  izlog = 1
                  inocm = 1
                  inolg = 1

                  if( itanl(m) .gt. 0 )
     &            call anset(itang(m),itanl(m),
     &                       form,xfac,afac,cmin,cmax,izlog,inocm,inolg)

                  if( form .le. 1.0 ) then

                     scal = form**0.35
                     xfac = xfac / form**0.5
                     xorg = 0.0
                     yorg = min( 20.d0, ( 1.0 / form - 1.0 ) / 2.5 )

                  else

                     scal = 1.0 / form**0.41
                     xfac = xfac / form**0.5
                     xorg = min( 20.d0, ( form - 1.0 ) / 1.5)
                     yorg = 0.0

                  end if

               end if

                  write(iot,'( "set: c1[",f6.3,
     &                            "] c2[",f6.3,
     &                            "] c5[",f6.3,
     &                            "] c6[",f6.3,
     &                            "] c7[",f6.3,
     &                            "] c8[",f6.3,"]")')
     &                    form, xfac, afac, scal, xorg, yorg
                  write(iot,'( "p: h2fs form[c1] xfac[c2]",
     &                 " afac[c5] scal[c6] xorg[c7] yorg[c8] nosp")')

               if( cmin .gt. 0.0 .and. cmax .gt. cmin ) then

                  write(iot,'( "set: c3[",1p1e15.6,"] c4[",
     &                  1p1e15.6,"]")') cmin, cmax

               end if

*-----------------------------------------------------------------------

               xmin = xm(1)
               xmax = xm(nx+1)
               ymin = ym(1)
               ymax = ym(ny+1)

               write(iot,'( "p: xmin(",1p1g14.7,") xmax(",1p1g14.7,
     &                       ")")') xmin, xmax

               write(iot,'( "p: ymin(",1p1g14.7,") ymax(",1p1g14.7,
     &                       ")")') ymin, ymax

               if( itanl(m) .gt. 0 ) then

                  write(iot,'( "p: ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               end if

               if( itsans(m) .gt. 0 ) then

                  call write_sangel(iot,m,0)

               end if

*-----------------------------------------------------------------------
*        rshow
*-----------------------------------------------------------------------

               write(iot,'(/"#",78("-"))')
               write(iot,'( "# rshow")')
               write(iot,'( "#",78("-"))')
               write(iot,'(/"p: legs[c5*0.875]")')

               none = 1
               iaxs = 1
               iuni = itrsh(m) * 2 - 1 + igsh
               ires = itres(m)
               widt = rtwid(m)

            if( igsh .eq. 0 ) then

               call gshow(1,iot,iaxs,iuni,widt,ires,
     &                    inum,igsh,igser,gfnam,igfmn,
     &                    nx+1,ny+1,none,xm,ym,zval,ixyz(1),
     &                    nr,mr,kr,val,itmtr(m,4),itglt(m))

            else

               call gshow(0,iot,iaxs,iuni,widt,ires,
     &                    inum,igsh,igser,gfnam,igfmn,
     &                    nx+1,ny+1,none,xm,ym,zval,ixyz(1),
     &                    1,1,krr,vll,itmtr(m,4),itglt(m))

            end if

*-----------------------------------------------------------------------

      if( inolg .eq. 1 .and. igsh .eq. 0 .and.
     &    cmin .gt. 0.0 .and. cmax .gt. cmin ) then

      write(iot,'(/"#",78("-"))')
      write(iot,'(
     &"z: xorg(1.03)"/
     &"p: xfac[c2*0.05] form[c1/0.05] notc noxt noyt dmin(1)"/
     &"p: ymin(0) ymax(1) xmin(0) xmax(1) cmin(1) cmax(100) nosp"/
     &"hc: y= 0.005 to 0.995 by 0.01 ; x= 0.5 to 0.5 by 1 ;"/
     &" 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20"/
     &"21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40"/
     &"41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60"/
     &"61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80"/
     &"81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100"/
     &"z: xorg(0.0)"/
     &"p: xfac[c2*0.05] form[c1/0.05] ",
     &"nosp afac[c5*0.625] ymin[c3] ymax[c4]"/
     &"p: noxt noxn ytxt(-1) ynum(-1) xltd(-1) xstd(-1) itic(1) ylog"
     &)')

      if( itazl(m) .eq. 0 ) then

         write(iot,'("y: Number ",a15)') hsunit(itunt(m))

      else

         write(iot,'("y: ",200a1)') (itazt(m)(i:i),i=1,itazl(m))

      end if

      if( izlog .eq. 0 ) write(iot,'( "p: ylin")')

      write(iot,'( "#",78("-"))')

      end if

*-----------------------------------------------------------------------

            if( inocm .eq. 1 ) then

               if( nn .eq. 0 ) then

                  write(iot,'(/"wt: s[c5]",/
     &                        a1,"vspace{-3}"/
     &                        "  z &=&",1pe13.4," [cm]"/
     &                        "  Z &=&   all"/
     &                        "e:")')
     &                        yen, zval

               else if( ia .eq. 0 ) then

                  write(iot,'(/"wt: s[c5]",/
     &                        a1,"vspace{-3}"/
     &                        "  z &=&",1pe13.4," [cm]"/
     &                        "  Z &=& ",i3,"   :  ",a3/
     &                        "e:")')
     &                        yen,
     &                        zval, iz, elmnt(iz)

               else

                  write(iot,'(/"wt: s[c5]",/
     &                        a1,"vspace{-3}"/
     &                        "  z     &=&",1pe13.4," [cm]"/
     &                        "  nucl. &=&   ",a8/
     &                        "e:")')
     &                        yen,
     &                        zval, chau

               end if

            end if

*-----------------------------------------------------------------------

      if( inolg .eq. 1 .and. igsh .eq. 0 .and.
     &    cmin .gt. 0.0 .and. cmax .gt. cmin ) then

         write(iot,'(/"#",78("-"))')
         write(iot,'("z: xorg[-1.03/0.05]"/
     &               "p: form[c1] xfac[c2] afac[c5] nosp notf")')

               write(iot,'( "p: xmin(",1p1g14.7,") xmax(",1p1g14.7,
     &                       ")")') xmin, xmax

               write(iot,'( "p: ymin(",1p1g14.7,") ymax(",1p1g14.7,
     &                       ")")') ymin, ymax

      end if

*-----------------------------------------------------------------------

            end do
            end do

*-----------------------------------------------------------------------
*        yz axis (matrix)
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 10 ) then

*-----------------------------------------------------------------------

               inum = 0

                  if( nn .eq. 0 .or. igsh .ne. 0  ) then

                     nc = 1

                  else

                     nc = nn

                  end if

*-----------------------------------------------------------------------

            do jx = 1, nx
            do ic = 1, nc

               xval = ( xm(jx) + xm(jx+1) ) / 2.0d0

               inum = inum + 1

*-----------------------------------------------------------------------

               if( nn .gt. 0 ) then

                  iz = nt(ic) / 1000

                  ia = nt(ic) - iz * 1000

                  if( ia .gt. 0 ) then

                     in = ia - iz

                     call chname(idum,ia,iz,chau)

                  end if

               end if

*-----------------------------------------------------------------------

                  write(iot,'(/"#",78("-"))')

               if( inum .eq. 1 ) then
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

               if( rtfac(m) .lt. 0.d0 )
     &           write(iot,
     &           '( "#restart maximum for normalization:",1p1e14.7)')
     &           facmax(m)

            if( nn .eq. 0 ) then

               write(iot,'("#   no. =",i3,3x,
     &         "ix  =",i3,3x,"Z = all")')
     &                     inum, jx

            else if( ia .eq. 0 ) then

               write(iot,'("#   no. =",i3,3x,
     &         "ix  =",i3,3x,"Z = ",i3," : ",a3)')
     &                     inum, jx, iz, elmnt(iz)

            else

               write(iot,'("#   no. =",i3,3x,
     &         "ix  =",i3,3x,a8)')
     &                     inum, jx, chau

            end if

*-----------------------------------------------------------------------

               if( nn .eq. 0 ) then

                  write(iot,'(/a1,"no. =",i3,
     &                        ",    ix =",i3,
     &                        ",    Z  =  all", a1)')
     &                        cha, inum, jx, cha

               else if( ia .eq. 0 ) then

                  write(iot,'(/a1,"no. =",i3,
     &                        ",    ix =",i3,
     &                        ",    Z  =  ",i3, a1)')
     &                        cha, inum, jx, iz, cha

               else

                  write(iot,'(/a1,"no. =",i3,
     &                        ",    ix =",i3,
     &                        ",    ",a8, a1)')
     &                        cha, inum, jx, chau, cha

               end if

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: z [cm]")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: y [cm]")')

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

*-----------------------------------------------------------------------

         if( igsh .eq. 0 ) then

               do i = 1, nr

                     sek = 0.0

                  if( nn .eq. 0 ) then

                     do kn = 1, mn
                     do kz = 1, mz
                     do il = 0, mm ! frtati 2022/02/18

                        lz  = ikzz( kz, kn )
                        ln  = iknn( kz, kn )

                        sek = sek + tr(i,lz,ln,il,1)

                     end do
                     end do
                     end do

                  else if( ia .eq. 0 ) then

                     do kn = 1, mn
                     do il = 0, mm ! frtati 2022/02/18

                        lz  = ikzz( iz, kn )
                        ln  = iknn( iz, kn )

                        sek = sek + tr(i,lz,ln,il,1)

                     end do
                     end do

                  else

                        lz  = ikzz( iz, in )
                        ln  = iknn( iz, in )

                     do il = 0, mm ! frtati 2022/02/18
                        sek = sek + tr(i,lz,ln,il,1)
                     end do

                  end if

                  val(i) = sek

               end do

                  izlog = 1
                  if( itanl(m) .gt. 0 )
     &            call anset(itang(m),itanl(m),
     &                       form,xfac,afac,cmin,cmax,izlog,inocm,inolg)

            if( izlog .eq. 0 ) then

               do i = 1, nr

                  if( val(i) .gt. dnon .and. cmin .ne. cmax ) then

                     val(i) = ( val(i) - cmin ) / ( cmax - cmin )
                     val(i) = max(0.0d0,val(i))
                     val(i) = min(1.0d0,val(i))

                  else

                     val(i) = -1.0

                  end if

               end do

            else

               do i = 1, nr

                  if( val(i) .gt. dnon .and. cmin .ne. cmax ) then

                     val(i) = log10(val(i)/cmin) / log10(cmax/cmin)
                     val(i) = max(0.0d0,val(i))
                     val(i) = min(1.0d0,val(i))

                  else

                     val(i) = -1.0d0

                  end if

               end do

            end if

         end if

*-----------------------------------------------------------------------

               if( inum .eq. 1 ) then

                  form  = ( ym(ny+1) - ym(1) )
     &                  / ( zm(nz+1) - zm(1) )
                  xfac  = 0.9
                  afac  = 0.8
                  izlog = 1
                  inocm = 1
                  inolg = 1

                  if( itanl(m) .gt. 0 )
     &            call anset(itang(m),itanl(m),
     &                       form,xfac,afac,cmin,cmax,izlog,inocm,inolg)

                  if( form .le. 1.0 ) then

                     scal = form**0.35
                     xfac = xfac / form**0.5
                     xorg = 0.0
                     yorg = min( 20.d0, ( 1.0 / form - 1.0 ) / 2.5 )

                  else

                     scal = 1.0 / form**0.41
                     xfac = xfac / form**0.5
                     xorg = min( 20.d0, ( form - 1.0 ) / 1.5)
                     yorg = 0.0

                  end if

               end if

                  write(iot,'( "set: c1[",f6.3,
     &                            "] c2[",f6.3,
     &                            "] c5[",f6.3,
     &                            "] c6[",f6.3,
     &                            "] c7[",f6.3,
     &                            "] c8[",f6.3,"]")')
     &                    form, xfac, afac, scal, xorg, yorg
                  write(iot,'( "p: h2fs form[c1] xfac[c2]",
     &                 " afac[c5] scal[c6] xorg[c7] yorg[c8] nosp")')

               if( cmin .gt. 0.0 .and. cmax .gt. cmin ) then

                  write(iot,'( "set: c3[",1p1e15.6,"] c4[",
     &                  1p1e15.6,"]")') cmin, cmax

               end if

*-----------------------------------------------------------------------

               zmin = zm(1)
               zmax = zm(nz+1)
               ymin = ym(1)
               ymax = ym(ny+1)

               write(iot,'( "p: xmin(",1p1g14.7,") xmax(",1p1g14.7,
     &                       ")")') zmin, zmax

               write(iot,'( "p: ymin(",1p1g14.7,") ymax(",1p1g14.7,
     &                       ")")') ymin, ymax

               if( itanl(m) .gt. 0 ) then

                  write(iot,'( "p: ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               end if

               if( itsans(m) .gt. 0 ) then

                  call write_sangel(iot,m,0)

               end if

*-----------------------------------------------------------------------
*        rshow
*-----------------------------------------------------------------------

               write(iot,'(/"#",78("-"))')
               write(iot,'( "# rshow")')
               write(iot,'( "#",78("-"))')
               write(iot,'(/"p: legs[c5*0.875]")')

               none = 1
               iaxs = 2
               iuni = itrsh(m) * 2 - 1 + igsh
               ires = itres(m)
               widt = rtwid(m)

            if( igsh .eq. 0 ) then

               call gshow(1,iot,iaxs,iuni,widt,ires,
     &                    inum,igsh,igser,gfnam,igfmn,
     &                    nz+1,ny+1,none,zm,ym,xval,ixyz(1),
     &                    nr,mr,kr,val,itmtr(m,4),itglt(m))

            else

               call gshow(0,iot,iaxs,iuni,widt,ires,
     &                    inum,igsh,igser,gfnam,igfmn,
     &                    nz+1,ny+1,none,zm,ym,xval,ixyz(1),
     &                    1,1,krr,vll,itmtr(m,4),itglt(m))

            end if

*-----------------------------------------------------------------------

      if( inolg .eq. 1 .and. igsh .eq. 0 .and.
     &    cmin .gt. 0.0 .and. cmax .gt. cmin ) then

      write(iot,'(/"#",78("-"))')
      write(iot,'(
     &"z: xorg(1.03)"/
     &"p: xfac[c2*0.05] form[c1/0.05] notc noxt noyt dmin(1)"/
     &"p: ymin(0) ymax(1) xmin(0) xmax(1) cmin(1) cmax(100) nosp"/
     &"hc: y= 0.005 to 0.995 by 0.01 ; x= 0.5 to 0.5 by 1 ;"/
     &" 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20"/
     &"21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40"/
     &"41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60"/
     &"61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80"/
     &"81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100"/
     &"z: xorg(0.0)"/
     &"p: xfac[c2*0.05] form[c1/0.05] ",
     &"nosp afac[c5*0.625] ymin[c3] ymax[c4]"/
     &"p: noxt noxn ytxt(-1) ynum(-1) xltd(-1) xstd(-1) itic(1) ylog"
     &)')

      if( itazl(m) .eq. 0 ) then

         write(iot,'("y: Number ",a15)') hsunit(itunt(m))

      else

         write(iot,'("y: ",200a1)') (itazt(m)(i:i),i=1,itazl(m))

      end if

      if( izlog .eq. 0 ) write(iot,'( "p: ylin")')

      write(iot,'( "#",78("-"))')

      end if

*-----------------------------------------------------------------------

            if( inocm .eq. 1 ) then

               if( nn .eq. 0 ) then

                  write(iot,'(/"wt: s[c5]",/
     &                        a1,"vspace{-3}"/
     &                        "  x &=&",1pe13.4," [cm]"/
     &                        "  Z &=&   all"/
     &                        "e:")')
     &                        yen, xval

               else if( ia .eq. 0 ) then

                  write(iot,'(/"wt: s[c5]",/
     &                        a1,"vspace{-3}"/
     &                        "  x &=&",1pe13.4," [cm]"/
     &                        "  Z &=& ",i3,"   :  ",a3/
     &                        "e:")')
     &                        yen,
     &                        xval, iz, elmnt(iz)

               else

                  write(iot,'(/"wt: s[c5]",/
     &                        a1,"vspace{-3}"/
     &                        "  x     &=&",1pe13.4," [cm]"/
     &                        "  nucl. &=&   ",a8/
     &                        "e:")')
     &                        yen,
     &                        xval, chau

               end if

            end if

*-----------------------------------------------------------------------

      if( inolg .eq. 1 .and. igsh .eq. 0 .and.
     &    cmin .gt. 0.0 .and. cmax .gt. cmin ) then

         write(iot,'(/"#",78("-"))')
         write(iot,'("z: xorg[-1.03/0.05]"/
     &               "p: form[c1] xfac[c2] afac[c5] nosp notf")')

               write(iot,'( "p: xmin(",1p1g14.7,") xmax(",1p1g14.7,
     &                       ")")') zmin, zmax

               write(iot,'( "p: ymin(",1p1g14.7,") ymax(",1p1g14.7,
     &                       ")")') ymin, ymax

      end if

*-----------------------------------------------------------------------

            end do
            end do

*-----------------------------------------------------------------------
*        xz axis (matrix)
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 11 ) then

*-----------------------------------------------------------------------

               inum = 0

                  if( nn .eq. 0 .or. igsh .ne. 0  ) then

                     nc = 1

                  else

                     nc = nn

                  end if

*-----------------------------------------------------------------------

            do jy = 1, ny
            do ic = 1, nc

               yval = ( ym(jy) + ym(jy+1) ) / 2.0d0

               inum = inum + 1

*-----------------------------------------------------------------------

               if( nn .gt. 0 ) then

                  iz = nt(ic) / 1000

                  ia = nt(ic) - iz * 1000

                  if( ia .gt. 0 ) then

                     in = ia - iz

                     call chname(idum,ia,iz,chau)

                  end if

               end if

*-----------------------------------------------------------------------

                  write(iot,'(/"#",78("-"))')

               if( inum .eq. 1 ) then
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

               if( rtfac(m) .lt. 0.d0 )
     &           write(iot,
     &           '( "#restart maximum for normalization:",1p1e14.7)')
     &           facmax(m)

            if( nn .eq. 0 ) then

               write(iot,'("#   no. =",i3,3x,
     &         "iy  =",i3,3x,"Z = all")')
     &                     inum, jy

            else if( ia .eq. 0 ) then

               write(iot,'("#   no. =",i3,3x,
     &         "iy  =",i3,3x,"Z = ",i3," : ",a3)')
     &                     inum, jy, iz, elmnt(iz)

            else

               write(iot,'("#   no. =",i3,3x,
     &         "iy  =",i3,3x,a8)')
     &                     inum, jy, chau

            end if

*-----------------------------------------------------------------------

               if( nn .eq. 0 ) then

                  write(iot,'(/a1,"no. =",i3,
     &                        ",    iy =",i3,
     &                        ",    Z  =  all", a1)')
     &                        cha, inum, jy, cha

               else if( ia .eq. 0 ) then

                  write(iot,'(/a1,"no. =",i3,
     &                        ",    iy =",i3,
     &                        ",    Z  =  ",i3, a1)')
     &                        cha, inum, jy, iz, cha

               else

                  write(iot,'(/a1,"no. =",i3,
     &                        ",    iy =",i3,
     &                        ",    ",a8, a1)')
     &                        cha, inum, jy, chau, cha

               end if

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: z [cm]")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: x [cm]")')

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

*-----------------------------------------------------------------------

         if( igsh .eq. 0 ) then

               do i = 1, nr

                     sek = 0.0

                  if( nn .eq. 0 ) then

                     do kn = 1, mn
                     do kz = 1, mz
                     do il = 0, mm ! frtati 2022/02/18

                        lz  = ikzz( kz, kn )
                        ln  = iknn( kz, kn )

                        sek = sek + tr(i,lz,ln,il,1)

                     end do
                     end do
                     end do

                  else if( ia .eq. 0 ) then

                     do kn = 1, mn
                     do il = 0, mm ! frtati 2022/02/18

                        lz  = ikzz( iz, kn )
                        ln  = iknn( iz, kn )

                        sek = sek + tr(i,lz,ln,il,1)

                     end do
                     end do

                  else

                        lz  = ikzz( iz, in )
                        ln  = iknn( iz, in )

                     do il = 0, mm ! frtati 2022/02/18
                        sek = sek + tr(i,lz,ln,il,1)
                     end do

                  end if

                  val(i) = sek

               end do

                  izlog = 1
                  if( itanl(m) .gt. 0 )
     &            call anset(itang(m),itanl(m),
     &                       form,xfac,afac,cmin,cmax,izlog,inocm,inolg)

            if( izlog .eq. 0 ) then

               do i = 1, nr

                  if( val(i) .gt. dnon .and. cmin .ne. cmax ) then

                     val(i) = ( val(i) - cmin ) / ( cmax - cmin )
                     val(i) = max(0.0d0,val(i))
                     val(i) = min(1.0d0,val(i))

                  else

                     val(i) = -1.0

                  end if

               end do

            else

               do i = 1, nr

                  if( val(i) .gt. dnon .and. cmin .ne. cmax ) then

                     val(i) = log10(val(i)/cmin) / log10(cmax/cmin)
                     val(i) = max(0.0d0,val(i))
                     val(i) = min(1.0d0,val(i))

                  else

                     val(i) = -1.0d0

                  end if

               end do

            end if

         end if

*-----------------------------------------------------------------------

               if( inum .eq. 1 ) then

                  form  = ( xm(nx+1) - xm(1) )
     &                  / ( zm(nz+1) - zm(1) )
                  xfac  = 0.9
                  afac  = 0.8
                  izlog = 1
                  inocm = 1
                  inolg = 1

                  if( itanl(m) .gt. 0 )
     &            call anset(itang(m),itanl(m),
     &                       form,xfac,afac,cmin,cmax,izlog,inocm,inolg)

                  if( form .le. 1.0 ) then

                     scal = form**0.35
                     xfac = xfac / form**0.5
                     xorg = 0.0
                     yorg = min( 20.d0, ( 1.0 / form - 1.0 ) / 2.5 )

                  else

                     scal = 1.0 / form**0.41
                     xfac = xfac / form**0.5
                     xorg = min( 20.d0, ( form - 1.0 ) / 1.5)
                     yorg = 0.0

                  end if

               end if

                  write(iot,'( "set: c1[",f6.3,
     &                            "] c2[",f6.3,
     &                            "] c5[",f6.3,
     &                            "] c6[",f6.3,
     &                            "] c7[",f6.3,
     &                            "] c8[",f6.3,"]")')
     &                    form, xfac, afac, scal, xorg, yorg
                  write(iot,'( "p: h2fs form[c1] xfac[c2]",
     &                 " afac[c5] scal[c6] xorg[c7] yorg[c8] nosp")')

               if( cmin .gt. 0.0 .and. cmax .gt. cmin ) then

                  write(iot,'( "set: c3[",1p1e15.6,"] c4[",
     &                  1p1e15.6,"]")') cmin, cmax

               end if

*-----------------------------------------------------------------------

               zmin = zm(1)
               zmax = zm(nz+1)
               xmin = xm(1)
               xmax = xm(nx+1)

               write(iot,'( "p: xmin(",1p1g14.7,") xmax(",1p1g14.7,
     &                       ")")') zmin, zmax

               write(iot,'( "p: ymin(",1p1g14.7,") ymax(",1p1g14.7,
     &                       ")")') xmin, xmax

               if( itanl(m) .gt. 0 ) then

                  write(iot,'( "p: ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               end if

               if( itsans(m) .gt. 0 ) then

                  call write_sangel(iot,m,0)

               end if

*-----------------------------------------------------------------------
*        rshow
*-----------------------------------------------------------------------

               write(iot,'(/"#",78("-"))')
               write(iot,'( "# rshow")')
               write(iot,'( "#",78("-"))')
               write(iot,'(/"p: legs[c5*0.875]")')

               none = 1
               iaxs = 3
               iuni = itrsh(m) * 2 - 1 + igsh
               ires = itres(m)
               widt = rtwid(m)

            if( igsh .eq. 0 ) then

               call gshow(1,iot,iaxs,iuni,widt,ires,
     &                    inum,igsh,igser,gfnam,igfmn,
     &                    nz+1,nx+1,none,zm,xm,yval,ixyz(1),
     &                    nr,mr,kr,val,itmtr(m,4),itglt(m))

            else

               call gshow(0,iot,iaxs,iuni,widt,ires,
     &                    inum,igsh,igser,gfnam,igfmn,
     &                    nz+1,nx+1,none,zm,xm,yval,ixyz(1),
     &                    1,1,krr,vll,itmtr(m,4),itglt(m))

            end if

*-----------------------------------------------------------------------

      if( inolg .eq. 1 .and. igsh .eq. 0 .and.
     &    cmin .gt. 0.0 .and. cmax .gt. cmin ) then

      write(iot,'(/"#",78("-"))')
      write(iot,'(
     &"z: xorg(1.03)"/
     &"p: xfac[c2*0.05] form[c1/0.05] notc noxt noyt dmin(1)"/
     &"p: ymin(0) ymax(1) xmin(0) xmax(1) cmin(1) cmax(100) nosp"/
     &"hc: y= 0.005 to 0.995 by 0.01 ; x= 0.5 to 0.5 by 1 ;"/
     &" 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20"/
     &"21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40"/
     &"41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60"/
     &"61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80"/
     &"81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100"/
     &"z: xorg(0.0)"/
     &"p: xfac[c2*0.05] form[c1/0.05] ",
     &"nosp afac[c5*0.625] ymin[c3] ymax[c4]"/
     &"p: noxt noxn ytxt(-1) ynum(-1) xltd(-1) xstd(-1) itic(1) ylog"
     &)')

      if( itazl(m) .eq. 0 ) then

         write(iot,'("y: Number ",a15)') hsunit(itunt(m))

      else

         write(iot,'("y: ",200a1)') (itazt(m)(i:i),i=1,itazl(m))

      end if

      if( izlog .eq. 0 ) write(iot,'( "p: ylin")')

      write(iot,'( "#",78("-"))')

      end if

*-----------------------------------------------------------------------

            if( inocm .eq. 1 ) then

               if( nn .eq. 0 ) then

                  write(iot,'(/"wt: s[c5]",/
     &                        a1,"vspace{-3}"/
     &                        "  y &=&",1pe13.4," [cm]"/
     &                        "  Z &=&   all"/
     &                        "e:")')
     &                        yen, yval

               else if( ia .eq. 0 ) then

                  write(iot,'(/"wt: s[c5]",/
     &                        a1,"vspace{-3}"/
     &                        "  y &=&",1pe13.4," [cm]"/
     &                        "  Z &=& ",i3,"   :  ",a3/
     &                        "e:")')
     &                        yen,
     &                        yval, iz, elmnt(iz)

               else

                  write(iot,'(/"wt: s[c5]",/
     &                        a1,"vspace{-3}"/
     &                        "  y     &=&",1pe13.4," [cm]"/
     &                        "  nucl. &=&   ",a8/
     &                        "e:")')
     &                        yen,
     &                        yval, chau

               end if

            end if

*-----------------------------------------------------------------------

      if( inolg .eq. 1 .and. igsh .eq. 0 .and.
     &    cmin .gt. 0.0 .and. cmax .gt. cmin ) then

         write(iot,'(/"#",78("-"))')
         write(iot,'("z: xorg[-1.03/0.05]"/
     &               "p: form[c1] xfac[c2] afac[c5] nosp notf")')

               write(iot,'( "p: xmin(",1p1g14.7,") xmax(",1p1g14.7,
     &                       ")")') zmin, zmax

               write(iot,'( "p: ymin(",1p1g14.7,") ymax(",1p1g14.7,
     &                       ")")') xmin, xmax

      end if

*-----------------------------------------------------------------------

            end do
            end do

*-----------------------------------------------------------------------

         end if

            call prestart(m,iot) !OBINATA(2012.8.20)

            close(iot)

         if( iteps(m) .ne. 0 .and. itaxs(m,iax) .ne. 13 ) then

            open(iot, file = fname, status = 'unknown' )
            call a_angel(idasa,fname)

         end if

  900 continue

*-----------------------------------------------------------------------
      deallocate( ixyz )
      return
      end


************************************************************************
*                                                                      *
      subroutine tyiltet(ncol,m,mz,mn,mm,nl,lt, ! frtati 2022/02/18 added mm
     &                   nr,mr,nm,mt,
     &                   kr,ikzz,iknn,tr,
     &                   trEVENT)
*                                                                      *
*       nuclear yield (or production) tally in tetra mesh              *
*       Last Modified by T.Furuta on 2025/01/16                        *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*      ncol  ..... reaction type                                       *
*             11 : termination by energy cut-off                       *
*             13 : nuclear reaction (n,x)                              *
*             14 : nuclear reaction (n,n'x)                            *
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
*        kcoll : =  0, normal                                          *
*                =  1, high energy fission                             *
*                =  2, high energy absorption                          *
*                =  3, low energy n elastic                            *
*                =  4, low energy n non-elastic                        *
*                =  5, low energy n fission                            *
*                =  6, low energy n absorption                         *
*                                                                      *
************************************************************************
      use QMD_COOD2_MOD, only : nqmdm
      use MMBANKMOD !FURUTA
      use partmod, only: itmxpt,itpan,itpat,jtpat,iznmmx,iznmturn ! frtati 2022/05/02
      use moddas_region ! S.H. (2022.11.9)
*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'param00.inc'
      include 'param01.inc'

*-----------------------------------------------------------------------

      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
!$OMP THREADPRIVATE(/icomon/)
      common /wtsave/ oldwt
!$OMP THREADPRIVATE(/wtsave/)
      common /tlgeom/ iblz1,iblz2
!$OMP THREADPRIVATE(/tlgeom/)
      common /tlglat/ ilev1,ilev2,ilat1(5,10),ilat2(5,10)
!$OMP THREADPRIVATE(/tlglat/)
      common /qparm/  ielas,icasc,iqstep,lvlopt,igamma
      common /bparm/  andt,jevap,npidk
      common /geosig/ geosig(250)
      common /cparm/  maxbch,maxcas
      common /mathzn/ mathz, mathn, jcoll, kcoll
!$OMP THREADPRIVATE(/mathzn/)
      common /eparm/  esmax, esmin, emin(20)

      common /regdc/  idrg(kvlmax), idgr(kvmmax)
*-----------------------------------------------------------------------

      common /clustf/ nclst, iclust(nnn)
!$OMP THREADPRIVATE(/clustf/)
      common /clustt/ nclsts, iclusts(nnn)
!$OMP THREADPRIVATE(/clustt/)
      common /clustw/ jclusts(0:8,nnn),  qclusts(0:12,nnn)
!$OMP THREADPRIVATE(/clustw/)
      common /clustl/ rumpal(0:20), numpal(0:20)
!$OMP THREADPRIVATE(/clustl/)
      common /cntcls/ jcount(3,nnn)
!$OMP THREADPRIVATE(/cntcls/)

*-----------------------------------------------------------------------

      dimension numsav(0:20), rumsav(0:20)

*-----------------------------------------------------------------------

      common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                itout(itlmax), ittwo(itlmax)
      common /tall12/ itman(itlmax), itmat(itlmax), itmct(itlmax),
     &                itnun(itlmax), itnuc(itlmax),
     &                itndz(itlmax), itndn(itlmax),
     &                itnkz(itlmax), itnkn(itlmax)
      common /tall22/ itmcn(itlmax), itmtn(itlmax), itmtt(itlmax)
      common /tall46/ itmtr(itlmax,4), rtmtr(itlmax,13)

      common /tall30/ itnda(itlmax)
      common /tall36/ itdpo(itlmax)
      common /tall37/ itcnt(9,itlmax)
      common /tall52/ itprd(itlmax)

      common /trstar/ itrstar
!$OMP THREADPRIVATE(/trstar/)

      common /tall82/ itcnth(9,itlmax)
      common /tall83/ itnzn(itlmax), itndm(itlmax)

      common /kmat1d/ idmn(0:kvlmax), idnm(kvmmax)
      common /clionprd/  lionprd

*-----------------------------------------------------------------------

      dimension   lt(nl)
      dimension   kr(mr)
      dimension   mt(nm)
      dimension   ikzz(maxpt,maxnt), iknn(maxpt,maxnt)
      dimension   tr(nr,mz,mn,0:mm,2)
      dimension   trEVENT(nr,mz,mn,0:mm)           !OBINATA(2012.8.20): as Ct
      real(8),allocatable,save:: tr0(:,:,:,:) !OBINATA(2012.8.20): as C

      dimension mnz(mxprodxs),mna(mxprodxs),xxn(mxprodxs),ildd(mxprodxs)

      dimension   ips(itmxpt) ! frtati 2021/10/05 6 -> itmxpt
      dimension   ncntt(3)

*-----------------------------------------------------------------------
      common /cntreg/ ncntc(3), ncreg(3), incrc(3), incrt(3), kcont(3)
      common /cntpat/ icpan(3), icpat(3,20,2)
      dimension     idas(1)
      equivalence ( das, idas )
      common /pnsave/ egs, uus, vvs, wws, wts, tms, nms, nct(3)
!$OMP THREADPRIVATE(/pnsave/)

*-----------------------------------------------------------------------
      common /cntmx/ ncntmx(3), nctmxsr(3)
!$OMP THREADPRIVATE(/cntmx/)

      common /mpi00/ npe, me
      real*8,allocatable :: tryld(:,:)

*-----------------------------------------------------------------------

      common /jcomon/ nabov,nobch,nocas,nomax
!$OMP THREADPRIVATE(/jcomon/)
      common /stat / istdev, irestart, ireschk

      integer, allocatable, save :: iclusav(:), jclusav(:,:)
      double precision, allocatable, save :: qclusav(:,:)
!$OMP THREADPRIVATE(iclusav, jclusav, qclusav)

      integer iii0,kkk0
      common /itettal2/ iii0,kkk0
!$OMP THREADPRIVATE(/itettal2/)

       if( .not. allocated(iclusav) ) then ! initial allocation
        allocate(iclusav(nqmdm),jclusav(0:8,nqmdm),qclusav(0:12,nqmdm))
        iclusav = 0
        jclusav = 0
        qclusav = 0.d0
       elseif( ubound(iclusav,1) .lt. nqmdm ) then ! extend array
        deallocate(iclusav,jclusav,qclusav)
        allocate(iclusav(nqmdm),jclusav(0:8,nqmdm),qclusav(0:12,nqmdm))
        iclusav = 0
        jclusav = 0
        qclusav = 0.d0
       endif

*-----------------------------------------------------------------------
* check of history counter
*-----------------------------------------------------------------------
      ihistcount = 0 ! history counter index, T.Sato 2022/12/20 for speed up
         if ( ncol.eq.0 .or. ncol.eq.4 ) then
           do i = 1, 3
             if( itcnth(i,m) .eq. 1 ) then
               if( ncntmx(i) .lt. itcnth(i*2+2,m) .or.
     &             ncntmx(i) .gt. itcnth(i*2+3,m) ) then
                 ihistcount = 1
               end if
             end if
           end do
         end if

*-----------------------------------------------------------------------
*        source particle or end of batch ( in case of istdev = 2 )
*
* OBINATA(2012.8.20): change tr(,,,,3) to trEVENT(,,,,)
*-----------------------------------------------------------------------

         if (( ncol .eq. 0 .or. ncol .eq. 4 )
     &                              .and. istdev .eq. 2) then
           if ((nocas.gt.1.or.ncol.eq.0) .and. ihistcount.ne.1 ) then

             tr(:,:,:,:,1) = tr(:,:,:,:,1) + trEVENT(:,:,:,:)
             tr(:,:,:,:,2) = tr(:,:,:,:,2) + trEVENT(:,:,:,:) ** 2

           end if

           trEVENT(:,:,:,:) = 0

         end if

*-----------------------------------------------------------------------
*        end of batch ( in case of istdev = 1 )
*
* OBINATA(2012.8.20): modificate for thread parallel
*-----------------------------------------------------------------------

         if ( ncol .eq. 0 .and. istdev .eq. 1) then
!$OMP MASTER
             allocate( tr0(nr,mz,mn,0:mm) ) ! frtati 2022/02/18 2 -> mm
             tr0(:,:,:,:) = 0.d0
!$OMP END MASTER
!$OMP BARRIER
!$OMP CRITICAL (tyiltet_crit_ist1)
             tr0(:,:,:,:) = tr0(:,:,:,:) + trEVENT(:,:,:,:)
!$OMP END CRITICAL (tyiltet_crit_ist1)
!$OMP BARRIER
!$OMP MASTER
             tr(:,:,:,:,1) = tr(:,:,:,:,1) + tr0(:,:,:,:) / maxcas
             tr(:,:,:,:,2) = tr(:,:,:,:,2)
     &                   + ( tr0(:,:,:,:) / maxcas ) ** 2
             deallocate( tr0 )
!$OMP END MASTER

           trEVENT(:,:,:,:) = 0

         end if

*-----------------------------------------------------------------------
*        rearrangement of iznm for MPI run ! frtati 2022/05/02
*-----------------------------------------------------------------------

         if( npe.gt.1 .and. itnzn(m).ne.0 .and. ncol.eq.0 ) then
           call paraiznm(m,1)
           if( me.gt.1 ) then
             allocate( tryld(iznmmx(m),2) )
             do ir = 1, nr
               do iz = 1, iznmmx(m)
                 tryld(iz,:) = tr(ir,iz,1,0,:)
               end do
               do iz = 1, iznmmx(m)
                 tr(ir,iz,1,0,:) = tryld(iznmturn(iz,m),:)
               end do
             end do
             deallocate( tryld )
           end if
         end if

*-----------------------------------------------------------------------
*        check of ncol and nclsts ( outgoing particles )
*-----------------------------------------------------------------------

         if( ncol .eq. 11 ) then

            if( itprd(m) .ne. 1 ) return

               iccol = 0
               npart = 1

               if( ityp .lt. 15 .or. ityp .gt. 19 ) return

               jz    = jtyp
               jn    = ktyp - ktyp / 1000000 * 1000000 - jz

               ipart = ityp
               tlw   = oldwt
               il = 0  ! T.Sato 2016/3/29, isomer is not considered for output=cutoff
               ncntt(1) = ncnt(ibknct+1,no,ipomp+1)
               ncntt(2) = ncnt(ibknct+2,no,ipomp+1)
               ncntt(3) = ncnt(ibknct+3,no,ipomp+1)

         else if( ncol .eq. 13 .or. ncol .eq. 14 ) then

            if ( itnda(m) .ge. 2 .and. jcoll .eq. 9 ) then
               if( nclsts .lt. 0 ) return
            else
               if( nclsts .le. 0 ) return
            end if

            if( ( itnda(m) .eq. 2 .and.
     &          ( jcoll .eq. 6 .or. jcoll .eq. 9 ) ) .or.
     &          ( itnda(m) .eq. 3 .and.
     &          ( jcoll .eq. 6 .or. jcoll .eq. 9 .or.
     &            jcoll .eq. 4 .or. jcoll .eq. 5 .or.
     &            jcoll .eq. 15 .or. jcoll .eq. 10 ) ) ) then

            else
               if( itdpo(m) .eq. 0 .and.
     &           ( jcoll .eq. 3 .or. kcoll .eq. 3 ) ) return
               if( itdpo(m) .eq. -1 .and. itrstar .eq. 1 ) return
            end if

               iccol = 1

               npart = nclsts

               ncntt(1) = jcount(1,1)
               ncntt(2) = jcount(2,1)
               ncntt(3) = jcount(3,1)

         else

            return

         end if

*-----------------------------------------------------------------------
*        check of counter
*-----------------------------------------------------------------------

            do i = 1, 3

               if( itcnt(i,m) .ne. 0 ) then

                  if( ncntt(i) .lt. itcnt(i*2+2,m) .or.
     &                ncntt(i) .gt. itcnt(i*2+3,m) ) return

               end if

            end do

*-----------------------------------------------------------------------
*        check of mat
*-----------------------------------------------------------------------

            if( nl .gt. 0 ) then

                  do i = 1, nl

                     if( itmcn(m) .gt. 0 .and.
     &                   idmn(mat) .eq. lt(i) ) goto 502
                     if( itmcn(m) .lt. 0 .and.
     &                   idmn(mat) .eq. lt(i) ) return

                  end do

                     if( itmcn(m) .gt. 0 ) return

            end if

  502          continue

*-----------------------------------------------------------------------

         if( iccol .eq. 1 ) then

*-----------------------------------------------------------------------
*           check of particles
*-----------------------------------------------------------------------

            call pcheck(m,itpan(m),ityp,ktyp,jtyp,ipn,ips)

               if( ipn .eq. 0 ) return

*-----------------------------------------------------------------------
*           mother and energy
*-----------------------------------------------------------------------

               ata   = dble( mathz + mathn )
               atz   = dble( mathz )
               mmas  = nint( ata )
               mchg  = nint( atz )
               eein  = ec(ibkec+no,ipomp+1)

*-----------------------------------------------------------------------
*           for specific mother nuclei
*-----------------------------------------------------------------------

            if( nm .gt. 0 ) then

               do i = 1, nm

                  iz = mt(i) / 1000
                  ia = mt(i) - iz * 1000

                  if( ( ia .gt. 0 .and.
     &                  mmas .eq. ia .and. mchg .eq. iz ) .or.
     &                ( ia .eq. 0 .and. mchg .eq. iz ) ) then

                     if( itmct(m) .gt. 0 ) goto 30
                     if( itmct(m) .lt. 0 ) return

                  end if

               end do

                  if( itmct(m) .gt. 0 ) return

            end if

   30       continue

         end if

*-----------------------------------------------------------------------
*     check of region and tetra
*-----------------------------------------------------------------------

               jj = 0

               ic=idgr(iblz1)
               itet=kkk0-10000
               ihelem=iii0

               call ttetck(ic,mr,kr,itet,ihelem,ir,icc)

               if( icc .ne. 1 ) goto 100

*-----------------------------------------------------------------------
*        normal weight
*-----------------------------------------------------------------------

               idoy = 1
               idev = 1
               wyld = oldwt

*-----------------------------------------------------------------------
*        for special : repeated nuclear reactions
*-----------------------------------------------------------------------

         if( iccol .eq. 1 .and. itspc(m) .gt. 0 .and.
     &       mathz .gt. 2 .and.
     &     ( jcoll .eq. 4 .or. jcoll .eq. 5 ) ) then

               idoy = itspc(m)
               idev = 10
               wyld = oldwt / dble( idoy * idev )

*-----------------------------------------------------------------------
*           save normal output and reaction parameters
*-----------------------------------------------------------------------

                  nclssav = nclsts

            do i = 1, nclsts

                  iclusav(i) = iclusts(i)

               do k = 0, 8

                  jclusav(k,i) = jclusts(k,i)

               end do

               do k = 0, 12

                  qclusav(k,i) = qclusts(k,i)

               end do

            end do

               do k = 0, 20

                  numsav(k) = numpal(k)
                  rumsav(k) = rumpal(k)

               end do

                  iprj  = ityp
                  kprj  = ktyp

                  if( ityp .lt. 15 ) then

                     bmax = sqrt( geosig(mmas) * 100.0 / 3.1415926 )

                  else if( ityp .ge. 15 ) then

                     ap = dble( ktyp - ktyp / 1000000 * 1000000 )
                     zp = dble( ktyp / 1000000 )

                     call sighi(ap,zp,eein,ata,zta,signe,sigel,bmax)

                  end if

         end if

*-----------------------------------------------------------------------
*     repeated do loop
*-----------------------------------------------------------------------

         do ireac = 1, idoy

*-----------------------------------------------------------------------
*     repeat calculation of nuclear reactions
*-----------------------------------------------------------------------

            if( idoy .gt. 1 .and. ireac .gt. 1 ) then

                  ipim = 0

   22          continue

                  call ncasc(1,iprj,kprj,eein,mmas,mchg,bmax)

               if( nclst .lt. 0 ) then

                     ipim = ipim + 1

                     if( ipim .le. 20 ) goto 22

               end if

            end if

*-----------------------------------------------------------------------
*     repeat calculation of evaporation
*-----------------------------------------------------------------------

            do ievap = 1, idev

               if( idev .gt. 1 ) then

                  call nevap(1)

                  do i = 1, nclsts

                     qclusts(8,i) = wyld * qclusts(8,i)

                  end do

                  npart = nclsts

               end if

*-----------------------------------------------------------------------
*        booking after nuclear reactions
*-----------------------------------------------------------------------

                  istat = 1
                  mnx   = 0

*-----------------------------------------------------------------------
*           NDATA=2,3 : Replace with yield data (activation cross section
*-----------------------------------------------------------------------

               if( iccol .eq. 1 .and.
     &           ( itnda(m) .eq. 2 .and.
     &           ( jcoll .eq. 6 .or. jcoll .eq. 9 ) ) .or.
     &           ( itnda(m) .eq. 3 .and.
     &           ( jcoll .eq. 6 .or. jcoll .eq. 9 .or.
     &             jcoll .eq. 4 .or. jcoll .eq. 5 .or.
     &             jcoll .eq. 15 .or. jcoll .eq. 10 ) ) ) then

                  mkk = mat
                  icl = idgr(iblz1)
                  jcl = jcoll

                  call prodxs2(ktyp,mchg,mmas,eein,
     &                 istat,mnx,mnz,mna,xxn,xxs,
     &                 ildd,mkk,icl,jcl,m)

                  iexclight = 0 ! Exclude light ions
                  if( itprd(m) .eq. 1 ) then ! output = cutoff
                     iexclight = 1
                     if ( jcoll.eq.6 .and.
     &                    ktyp.eq.2112 .and. eein.lt.20.0 )
     &                    iexclight = 0
                  end if

                  if( istat .eq. 0 .and. mnx .gt. 0 ) then

                   do nnx = 1, mnx

                      jz = mnz(nnx)
                      jn = mna(nnx) - mnz(nnx)
                      il = ildd(nnx)

                    if ( iexclight.ne.1 .or. jz.gt.2 ) then

                     if( jz .gt. 0 .and. jz .le. maxpt .and.
     &                      jn .gt. 0 .and. jn .le. maxnt ) then

                        kz = ikzz(jz,jn)
                        kn = iknn(jz,jn)

                      if( itnzn(m).eq.0 ) then ! frtati 2022/02/18
                       if( kz .le. mz .and. kn .le. mn ) then
                          trEVENT(ir,kz,kn,il) =
     &                         trEVENT(ir,kz,kn,il) + wyld * xxn(nnx)
                       end if

                      else
                         trEVENT(ir,igetiznm(kz,kn,il,m),1,0) =
     &                        trEVENT(ir,igetiznm(kz,kn,il,m),1,0) +
     &                        wyld * xxn(nnx)
                      end if

                     end if

                    end if

                   end do

                  end if

               end if

*-----------------------------------------------------------------------
*           nuclear data for 4He, 14N, 16O
*-----------------------------------------------------------------------

               if( iccol .eq. 1 .and.
     &             itnda(m) .eq. 1 .and. jcoll .eq. 4 ) then

                     call prodxs(ktyp,mchg,mmas,eein,
     &                           istat,mnx,mnz,mna,xxn,xxs,ildd)

                  if( istat .eq. 0 .and. mnx .gt. 0 ) then

                     do nnx = 1, mnx

                        jz = mnz(nnx)
                        jn = mna(nnx) - mnz(nnx)
                        il = ildd(nnx)

                        if( jz .gt. 0 .and. jz .le. maxpt .and.
     &                      jn .gt. 0 .and. jn .le. maxnt ) then

                           kz = ikzz(jz,jn)
                           kn = iknn(jz,jn)

                        if( itnzn(m).eq.0 ) then ! frtati 2022/02/18

                        if( kz .le. mz .and. kn .le. mn ) then

                           trEVENT(ir,kz,kn,il) =
     &                     trEVENT(ir,kz,kn,il) + wyld * xxn(nnx)

                        end if

                        else
                           trEVENT(ir,igetiznm(kz,kn,il,m),1,0) =
     &                     trEVENT(ir,igetiznm(kz,kn,il,m),1,0) +
     &                     wyld * xxn(nnx)
                        end if

                        end if

                     end do

                  end if

               end if


*-----------------------------------------------------------------------
*           normal case
*-----------------------------------------------------------------------

            if( npart .gt. 0 ) then

               do 500 j = 1, npart

*-----------------------------------------------------------------------

                  if( iccol .eq. 1 ) then

                     ipart = jclusts(3,j)
                     tlw   = qclusts(8,j)
                     jz    = jclusts(1,j)
                     jn    = jclusts(2,j)
                     il    = jclusts(8,j)

                  end if

*-----------------------------------------------------------------------
* When ndata=2 or 3, counter values change by ndata in [counter]
                 if ( (ipart.ge.3 .and. ipart.le.13) .or.
     &                 ipart.ge.19 ) then

                  if( ( itnda(m) .eq. 2 .or. itnda(m) .eq. 3 ) .and.
     &                 ( istat .eq. 0 .or. istat .eq. 3 ) ) then

                   if( ncntc(1) .eq. 1 .or. ncntc(2) .eq. 1 .or.
     &                    ncntc(3) .eq. 1 ) then
                    do k = 1, 3
                     ndatcount = 0
                     if( ncntc(k) .eq. 1 ) then
                      knn  = ityp
                      kcg  = jtyp
                      kkf  = ktyp
                      call pcchck(k,knn,kkf,kcg,icpan,icpat,icc)

                      if( icc .eq. 1 ) then
                       kdsm = kcont(k)
                       ldsm = 22
                       ndatcount = idas_kcont(kdsm+ldsm)
                      end if

                      if ( ndatcount .ne. 0 ) then

                       if ( ncol .eq. 13 ) then
                        ncnta(ibknct+k,j,ipomp+1)
     &                         = ncnta(ibknct+k,j,ipomp+1)
     &                         + ndatcount

                       else if ( ncol .eq. 14 ) then
                        if ( j .eq. 1 ) then
                         nct(k) = nct(k) + ndatcount
                        else if (j .gt. 1 ) then
                         ncnta(ibknct+k,j-1,ipomp+1)
     &                          = ncnta(ibknct+k,j-1,ipomp+1)
     &                          + ndatcount
                        end if

                       end if

                      end if

                     end if
                    end do
                   end if

                  end if

                 end if

*-----------------------------------------------------------------------

                  if( ipart .lt. 15 .or. ipart .gt. 19 ) goto 500

                  if( iccol .eq. 1 .and. itprd(m) .eq. 1 ) then

                     emint = emin(ipart) * dble( jz + jn )

                     if( qclusts(7,j) .gt. emint ) goto 500

                  end if

*-----------------------------------------------------------------------

                     if( ( itnda(m) .eq. 2 .or. itnda(m) .eq. 3 ) .and.
     &                 istat .eq. 0 .and.
     &                 ( iexclight.ne.1 .or. jz.gt.2 ) ) goto 500

                     if( ( itnda(m) .eq. 2 .or. itnda(m) .eq. 3 ) .and.
     &                    istat .eq. 3 .and. lionprd .ne. 1) goto 500

                     if( itnda(m) .eq. 1 .and.
     &                   istat .eq. 0 .and. mnx .gt. 0 ) then

                        do nnx = 1, mnx

                           if( jz .eq. mnz(nnx) .and.
     &                         jn .eq. mna(nnx) - mnz(nnx) ) goto 500

                        end do

                     end if

                        if( jz .gt. 0 .and. jz .le. maxpt .and.
     &                      jn .gt. 0 .and. jn .le. maxnt ) then

                           kz = ikzz(jz,jn)
                           kn = iknn(jz,jn)

                        if( itnzn(m).eq.0 ) then ! frtati 2022/02/18

                        if( kz .le. mz .and. kn .le. mn ) then
cOBINATA(2012.8.20): Ct = Ct + xi.wi
                           trEVENT(ir,kz,kn,il) = trEVENT(ir,kz,kn,il)
     &                                       + tlw
                        end if

                        else
                           trEVENT(ir,igetiznm(kz,kn,il,m),1,0) =
     &                     trEVENT(ir,igetiznm(kz,kn,il,m),1,0) + tlw
                        end if

                        end if

  500          continue

            end if

*-----------------------------------------------------------------------

            end do

         end do

*-----------------------------------------------------------------------
*        restore normal output
*-----------------------------------------------------------------------

         if( idoy .gt. 1 ) then

                  nclsts = nclssav

            do i = 1, nclsts

                  iclusts(i) = iclusav(i)
                  iclusav(i) = 0

               do k = 0, 8

                  jclusts(k,i) = jclusav(k,i)
                  jclusav(k,i) = 0

               end do

               do k = 0, 12

                  qclusts(k,i) = qclusav(k,i)
                  qclusav(k,i) = 0.d0

               end do

            end do

               do k = 0, 20

                  numpal(k) = numsav(k)
                  rumpal(k) = rumsav(k)

               end do

         end if

*-----------------------------------------------------------------------
*        store the decrease of mother by negative weight
*-----------------------------------------------------------------------

        if( ncol .ne. 11 ) then ! S.H. 2021.12.22
         if( itdpo(m) .eq. -1 .and. itrstar .eq. 0 ) then

            if( mathz .gt. 0 .and. mathz .le. maxpt .and.
     &          mathn .gt. 0 .and. mathn .le. maxnt ) then

               kz = ikzz(mathz,mathn)
               kn = iknn(mathz,mathn)

             if( itnzn(m).eq.0 ) then ! frtati 2022/02/18
              if( kz .le. mz .and. kn .le. mn ) then ! T.Sato 2024/05/02
               trEVENT(ir,kz,kn,0) = trEVENT(ir,kz,kn,0) - oldwt
              endif
             else
               trEVENT(ir,igetiznm(kz,kn,0,m),1,0) =
     &         trEVENT(ir,igetiznm(kz,kn,0,m),1,0) - oldwt
             end if

            end if

         end if
        end if

*-----------------------------------------------------------------------

  100 continue

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine pyildtet(m,mz,mn,mm, ! frtati 2022/02/18
     &                    nr,mr,nn,nt,ikzz,iknn,
     &                    tr,
     &                    nx,ny,nz,kr,xm,ym,zm,igsh,idasa)
*                                                                      *
*       output of nuclear yield (or production) tally in tetra mesh    *
*       Last modified by T.Furuta on 2025/01/16                        *
*                                                                      *
************************************************************************
      use sangelmod, only: itsans
*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'param01.inc'

*-----------------------------------------------------------------------

      common /taliin/ rsouin, nzztin, nrgnin

*-----------------------------------------------------------------------

      common /verjam/ versn, lastr, iyeav, imonv, idayv

      common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                itout(itlmax), ittwo(itlmax)
      common /tall02/ itrgn(itlmax), itrgm(itlmax), itreg(itlmax)
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

      common /tall21/ rtfac(itlmax)
      common /tall29/ itrsh(itlmax)
      common /tall35/ iteps(itlmax)
      common /tall38/ itres(itlmax)
      common /tall39/ rtwid(itlmax)
      common /tall46/ itmtr(itlmax,4), rtmtr(itlmax,13)
      common /tall49/ itglt(itlmax)
      common /tall83/ itnzn(itlmax), itndm(itlmax)

      common /istcut/ ist_cut, ist_bat

      common /fact01/ facmax(itlmax) ! kitamura23/03/31

*-----------------------------------------------------------------------

      common /volreg/ dvol(kvlmax)
      common /regdc/ idrg(kvlmax), idgr(kvmmax)

*-----------------------------------------------------------------------

      common /tcntl/  icntl, inucr

*-----------------------------------------------------------------------

      common /cparm/  maxbch,maxcas
      common /jcomon/ nabov,nobch,nocas,nomax
!$OMP THREADPRIVATE(/jcomon/)
* Use igamma
      common /qparm/  ielas,icasc,iqstep,lvlopt,igamma

      common /talout/ itall

      character fname*100, fnume*3

*-----------------------------------------------------------------------
      character gfnam*100
      common /tall60/ itger(itlmax)

      common /tall67/ itstd(itlmax), rtstd(itlmax)

      common /mpi00/ npe, me

*-----------------------------------------------------------------------

      dimension   kr(mr)
      dimension   nt(nn)
      dimension   ikzz(maxpt,maxnt), iknn(maxpt,maxnt)
      dimension   tm(maxnt+maxpt,2)
      dimension   tr(nr,mz,mn,0:mm,2) ! frtati 2022/02/18 2 -> mm

      dimension   xm(nx+1)
      dimension   ym(ny+1)
      dimension   zm(nz+1)

      integer,allocatable :: ixyz(:)
      integer,allocatable :: lr(:)        !FURUTA20190204
      real(8),allocatable :: vl(:),val(:) !FURUTA20190204
c Dont know why but necessary to avoid segmentation fault

*-----------------------------------------------------------------------

      character hsunit(2)*15

      data hsunit / '[1/source]     ',
     &              '[1/cm^3/source]'/

*-----------------------------------------------------------------------

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

      data ipstep / 12 /

      character erfnm*100

      character chau*8
      character cha*1
      data cha /"'"/

      character rpa*1
      data rpa /'}'/
*-----------------------------------------------------------------------
cFURUTA20190208 OpenFOAM output
      character(1),allocatable :: foamfIType(:)
      integer,allocatable :: foamfIndex(:)
      character(len=255) :: outFilename
      integer :: numIndex,ifilecount
      integer :: itfoam
      common /tall76/ itfoam(itlmax)
cFURUTA20191028 CSV output
      character(400) buf
      character(200) sbuf
      real(8),allocatable :: xcm(:,:)
*-----------------------------------------------------------------------
      character(28) cfmt                !FURUTA20200615
      character(12) cir                 !FURUTA20200615
      common /redufmt/ iredufmt(itlmax) !FURUTA20200615
*-----------------------------------------------------------------------
      character yen*1
      yen  = char(92)
      allocate( ixyz(max((nx+1)*(ny+1),(ny+1)*(nz+1),(nx+1)*(nz+1))) )

*-----------------------------------------------------------------------
*        set mesh volume and name
*-----------------------------------------------------------------------

               allocate ( lr(nr),vl(nr),val(nr) ) !FURUTA20190204
               call ttetvl(mr,kr,nr,vl,lr)

*-----------------------------------------------------------------------
*        relative error and  unit conversion
*        c1 : nomalization for source
*-----------------------------------------------------------------------

                  cmax = 0.0
                  cmin = 1.e+33
                  dnon = 1.e-33
                  isdz = 0
                  stdm = 0.0
                  itstd(m) = 1

         if( igsh .eq. 0 ) then

            if ( rsouin .gt. 0d0 ) then ! S.H. (2016.7.26)

               c1 = 1.0d+0 / rsouin

            else

               c1 = 0d0

            end if

            facmax(m) = 1.d0
            if( rtfac(m) .lt. 0.d0 ) then
               facmax(m) = 0.d0

               do 101 ir = 1, nr
               do 101 iz = 1, mz
               do 101 in = 1, mn
               do 101 il = 0, mm

                  if( tr(ir,iz,in,il,1) .ne. 0.d0 ) then

                     if( itunt(m) .eq. 1 ) then

                        fmaxfc = tr(ir,iz,in,il,1)

                     else if( itunt(m) .eq. 2 ) then

                        fmaxfc = tr(ir,iz,in,il,1) / vl(ir)

                     end if

                     if( facmax(m) .lt. fmaxfc) facmax(m) = fmaxfc

                  end if

  101          continue

            end if

            do 100 ir = 1, nr

               if( itunt(m) .eq. 1 ) then

                  cc = abs(rtfac(m)/facmax(m))

               else if( itunt(m) .eq. 2 ) then

                  cc = abs(rtfac(m)/facmax(m)) / vl(ir)

               end if

            do 100 iz = 1, mz
            do 100 in = 1, mn
            do 100 il = 0, mm ! frtati 2022/02/18

               if( tr(ir,iz,in,il,1) .ne. 0.d0 ) then

                  call calc_stdev(m,Xa,sigx,
     &                            tr(ir,iz,in,il,1),
     &                            tr(ir,iz,in,il,2),
     &                            cc)

                  tr(ir,iz,in,il,1) = Xa
                  tr(ir,iz,in,il,2) = sigx
                  if( sigx .gt. stdm ) stdm = sigx

                  if( tr(ir,iz,in,il,1) .gt. cmax )
     &                                   cmax = tr(ir,iz,in,il,1)

cABE 2022/03/02, avoid the negative value of cmin
                  if( tr(ir,iz,in,il,1) .gt. 0.d0 .and.
     &                tr(ir,iz,in,il,1) .lt. cmin )
     &                                   cmin = tr(ir,iz,in,il,1)

               else

                  if( il .eq. 0 ) isdz = 1
                  tr(ir,iz,in,il,2) = 0.0

               end if

  100       continue

            if( nobch .gt. ist_bat ) then
            if( rtstd(m) .gt. 0.d0 .and.
     &          stdm .gt. 0.d0 .and. isdz .eq. 0 .and.
     &          stdm .lt. rtstd(m) ) itstd(m) = 0
            end if

         end if

*-----------------------------------------------------------------------
*     out put unit = 31  : temporary number
*-----------------------------------------------------------------------

      do 900 iax = 1, itfln(m)

!OBINATA(2012.8.20): output *.err
        noe = 1
        if ( any( itaxs(m,iax) .eq. (/ 8, 13 /) )
     &      .and. ittwo(m) .ne. 4 .and. icntl .ne. 8 )
     &    noe = 2

        if(itaxs(m,iax).eq.13.and.iredufmt(m).eq.1) noe=1 !FURUTA20200615

        do 900 ioe = 1, noe

         if( igsh .ne. 0 .and.
     &     ( itaxs(m,iax) .lt. 7 .or. itrsh(m) .eq. 0 ) ) goto 900

         if( itall .eq. 2 .and. nobch .lt. maxbch .and.
     &       igsh .eq. 0 ) then

            write(fnume,'(i3.3)') nobch
!OBINATA(2012.8.20): output *.err
            if ( ioe .eq. 1 ) then
              fname = ctfln(m,iax)(1:itfll(m,iax))//'.'//fnume
            else
              call mk_2derrfn(ctfln(m,iax),fname,itfll(m,iax))
            end if

         else if ( itall.eq.4 .and. igsh.eq.0 ) then
            write(fnume,'(i3.3)') nobch
            if ( npe.gt.1 ) write (fnume,'(i3.3)') nobch/(npe-1)
            if ( ioe.eq.1 ) then
              call mk_2dnumfn(ctfln(m,iax),fname,itfll(m,iax),
     &             nobch,maxbch,npe)
            else
              call mk_2dnerfn(ctfln(m,iax),fname,itfll(m,iax),fnume)
            end if

         else

!OBINATA(2012.8.20): output *.err
            if ( ioe .eq. 1 ) then
              fname = ctfln(m,iax)
            else
              call mk_2derrfn(ctfln(m,iax),fname,itfll(m,iax))
            end if

         end if

            iot = 31
            open(iot, file = fname, status = 'unknown' )
            iot2 = 32 !FURUTA20190208 OpenFOAM output

*-----------------------------------------------------------------------
            gfnam = ctfln(m,iax)
            igfmn = itfll(m,iax)
            igser = itger(m)

*-----------------------------------------------------------------------
*     input echo
*-----------------------------------------------------------------------

               call tyilech(iot,m,iax,1)

*-----------------------------------------------------------------------
*        dchain axis
*-----------------------------------------------------------------------

         if( itaxs(m,iax) .eq. 13 ) then

*-----------------------------------------------------------------------
*           output error or not
*-----------------------------------------------------------------------

            if( itout(m) .ne. 0 ) then

               do i = itfll(m,iax), 1, -1

                  if( ctfln(m,iax)(i:i) .eq. '.' ) goto 50

               end do

                  i = itfll(m,iax)

   50             itfp = i - 1

               do i = 1, itfp

                  erfnm(i:i) = ctfln(m,iax)(i:i)

               end do

                  erfnm(itfp+1:itfp+4) = '.err'

                  iou = 15
                  open(iou, file = erfnm(1:itfp+4), status = 'unknown' )

                  call tyilech(iou,m,iax,1)

            end if

*-----------------------------------------------------------------------
*        output : do not change below expresion
*-----------------------------------------------------------------------

            write(iot,'(/"#",78("-"))')
            write(iot,'(/
     &           " tetra scoring mesh nuclear yield (or production)"/
     &           " ----------------------------------------")')

            if( itout(m) .ne. 0 ) then

               write(iou,'(/"#",78("-"))')
               write(iou,'(/
     &         " Statistical Error(%) for tetra scoring mesh nuclear ",
     &         "yield of above file."/
     &         " --------------------------------------------",
     &         "---------------------")')

            end if

            if(iredufmt(m).eq.0)then !FURUTA20200615
*-----------------------------------------------------------------------

             do 170 il = 0, 2 ! frtati 2022/03/11
            do 170 iz = 1, maxpt

               if( nn .gt. 0 ) then

                  do i = 1, nn

                     if( nt(i) / 1000 .eq. iz ) goto 150

                  end do

                     goto 170

               end if

  150          continue

                     do in = 1, maxnt
                     do ir = 1, nr

                        if( itnzn(m).eq.0 ) then
                          t0_tr = tr(ir,iz,in,il,1)
                        else
                          t0_tr = tr(ir,igetiznmp(iz,in,il,m),1,0,1)
                        end if
                        if( t0_tr .ne. 0.0 ) goto 120

                     end do
                     end do

                     goto 170

  120             im = in

                     do in = maxnt, im + 1, -1
                     do ir = 1, nr

                        if( itnzn(m).eq.0 ) then
                          t0_tr = tr(ir,iz,in,il,1)
                        else
                          t0_tr = tr(ir,igetiznmp(iz,in,il,m),1,0,1)
                        end if
                        if( t0_tr .ne. 0.0 ) goto 140

                     end do
                     end do

                  jm = im
  140             jm = in

                  km = jm - im + 1
                  lm = ( km - 1 ) / ipstep + 1

               do mmm = 1, lm

                  n1 = ipstep * ( mmm - 1 ) + 1
                  n2 = min( ipstep * mmm, km )
                  n3 = n1 + im - 1
                  n4 = n2 + im - 1

                  IF(il .eq. 0) then
                   write(iot,'(/1x,i4,"-",a2,
     &                       " isotope production",
     &                       " # n3 n4 = ",i3,1x,i3)')
     &                 iz, elmnt(iz), n3, n4

                  ELSEIF(il .eq. 1) then
                   write(iot,'(/1x,i4,"-",a2,
     &                       " 1st metastable isotope production",
     &                       " # n3 n4 = ",i3,1x,i3)')
     &                 iz, elmnt(iz), n3, n4

                  ELSEIF(il .eq. 2) then
                   write(iot,'(/1x,i4,"-",a2,
     &                       " 2nd metastable isotope production",
     &                       " # n3 n4 = ",i3,1x,i3)')
     &                 iz, elmnt(iz), n3, n4
                  ENDIF

                  IF(il .eq. 0) then
                     write(iot,'(" tetra",12i11)')
     &                      ( im - 1 + iz + n, n = n1, n2 )
                  ELSEIF(il .eq. 1) then
                   write(iot,
     &                  '(" tetra",i10,"m",i10,"m",i10,"m"
     &               ,i10,"m",i10,"m",i10,"m",i10,"m",i10,"m"
     &               ,i10,"m",i10,"m",i10,"m",i10,"m")')
     &                      ( im - 1 + iz + n, n = n1, n2 )
                  ELSEIF(il .eq. 2) then
                   write(iot,
     &                  '(" tetra",i10,"n",i10,"n",i10,"n"
     &               ,i10,"n",i10,"n",i10,"n",i10,"n",i10,"n"
     &               ,i10,"n",i10,"n",i10,"n",i10,"n")')
     &                      ( im - 1 + iz + n, n = n1, n2 )
                  ENDIF

                  do ir = 1, nr

c  *** Changed by T.Sato 2013/10/9, i5 -> i7
                   if( itnzn(m).eq.0 ) then ! frtati 2022/02/18
                     write(iot,'(i7,1p12e11.3)')
     &                     lr(ir), ( tr(ir,iz,i,il,ioe), i = n3, n4 )
                   else
                     write(iot,'(i7,1p12e11.3)') lr(ir),
     &               ( tr(ir,igetiznmp(iz,i,il,m),1,0,ioe),i = n3, n4 )
                   end if

                  end do


                  if( itout(m) .ne. 0 ) then

                  IF(il .eq. 0) then
                    write(iou,'(/1x,i4,"-",a2,
     &              " isotope production : ERROR(%)")')
     &              iz, elmnt(iz)

                  ELSEIF(il .eq. 1) then
                    write(iou,'(/1x,i4,"-",a2,
     &              " 1st metastable isotope production : ERROR(%)")')
     &              iz, elmnt(iz)

                  ELSEIF(il .eq. 2) then
                    write(iou,'(/1x,i4,"-",a2,
     &              " 2nd metastable isotope production : ERROR(%)")')
     &              iz, elmnt(iz)
                  ENDIF

                  IF(il .eq. 0) then
                     write(iou,'(" tetra",12i11)')
     &                      ( im - 1 + iz + n, n = n1, n2 )
                  ELSEIF(il .eq. 1) then
                   write(iou,
     &                  '(" tetra",i10,"m",i10,"m",i10,"m"
     &               ,i10,"m",i10,"m",i10,"m",i10,"m",i10,"m"
     &               ,i10,"m",i10,"m",i10,"m",i10,"m")')
     &                      ( im - 1 + iz + n, n = n1, n2 )
                  ELSEIF(il .eq. 2) then
                   write(iou,
     &                  '(" tetra",i10,"n",i10,"n",i10,"n"
     &               ,i10,"n",i10,"n",i10,"n",i10,"n",i10,"n"
     &               ,i10,"n",i10,"n",i10,"n",i10,"n")')
     &                      ( im - 1 + iz + n, n = n1, n2 )
                  ENDIF

                     do ir = 1, nr

                      if( itnzn(m).eq.0 ) then ! frtati 2022/02/18
                        write(iou,'(i7,1p12e11.3)')
     &                    lr(ir), ( tr(ir,iz,i,il,2)*100.0, i = n3, n4 )
                      else
                        write(iou,'(i7,1p12e11.3)') lr(ir),
     &                  ( tr(ir,igetiznmp(iz,i,il,m),1,0,2)*100.0,
     &                    i = n3, n4 )
                      end if

                     end do

                  end if

               end do

  170       continue

*-----------------------------------------------------------------------
           else
             write(iot,'(a)')
             write(iot,'("# num nucleusID yield r.err")')
             if(itout(m).ne.0)then
              write(iou,'(a)')
              write(iou,'("# num nucleusID yield r.err")')
             endif
             cfmt='(i#,x,i7,1p2e11.3)'
             do ir=1,nr
!NS 2021.04 change for NVIDIA TOOL KIT COMPILL ERR
              write(cir,"(I12)")ir
!NS 2021.04 end change for NVIDIA TOOL KIT COMPILL ERR
              nir=len_trim(adjustl(cir))
              write(cfmt(3:3),'(i1)')nir
              do iz=1,maxpt
               if(nn.gt.0) then
                do i = 1,nn
                 if(nt(i)/1000.eq.iz)exit
                enddo
                if(i.gt.nn)cycle
               endif
               do in=1,maxnt
                do il=0,2 ! frtati 2022/05/02
                 if( itnzn(m).eq.0 ) then ! frtati 2022/02/18
                  if(tr(ir,iz,in,il,ioe).ne.0d0)then
                   write(iot,cfmt)ir,iz*10000+(in+iz)*10+il,
     &                  (tr(ir,iz,in,il,ioee),ioee=1,2)
                   if(itout(m).ne.0)then
                    write(iou,cfmt)ir,iz*10000+(in+iz)*10+il,
     &                   tr(ir,iz,in,il,2)
                   endif
                  endif
                 else
                  if(tr(ir,igetiznmp(iz,in,il,m),1,0,ioe).ne.0d0)then
                   write(iot,cfmt)ir,iz*10000+(in+iz)*10+il,
     &             (tr(ir,igetiznmp(iz,in,il,m),1,0,ioee),ioee=1,2)
                   if(itout(m).ne.0)then
                    write(iou,cfmt)ir,iz*10000+(in+iz)*10+il,
     &                   tr(ir,igetiznmp(iz,in,il,m),1,0,2)
                   endif
                  end if
                 endif
                enddo
               enddo
              enddo
             enddo
             write(iot,cfmt)0,0,0.0d0,0.0d0
             if(itout(m).ne.0)write(iou,cfmt)0,0,0.0d0,0.0d0
            endif
*-----------------------------------------------------------------------

               if( itout(m) .ne. 0 ) close(iou)

*-----------------------------------------------------------------------
*        tet axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 14 ) then

               inum = 0

                  if( nn .eq. 0 ) then

                     nc = 1

                  else

                     nc = nn

                  end if

cFURUTA20190208 OpenFOAM output
                  if ( itfoam(m) .ne. 0) then
                   numIndex = 1
                   allocate( foamfIType(1:numIndex) )
                   allocate( foamfIndex(1:numIndex) )
                   foamfIType = (/ 'n' /)
                   foamfIndex = (/ nc /)
                   ifilecount=0
                  endif
cFURUTA20191028 CSV output
            if ( itfoam(m) .eq. 2) then
             allocate ( xcm(3,nr) ) !FURUTA20191028
             call ttetcm(mr,kr,nr,xcm)
            endif

*-----------------------------------------------------------------------

            do ic = 1, nc

               inum = inum + 1

*-----------------------------------------------------------------------
cFURUTA20190208 OpenFOAM output
               if ( itfoam(m) .ne. 0) then
                 ifilecount=ifilecount+1
                 call openfoam_create_filename(
     &                itfoam(m),fname, foamfIType, foamfIndex, numIndex,
     &                ifilecount, outFilename)
                 open(iot2, file = outFilename, status='unknown')
                 if(itfoam(m).eq.1)then
                  write(iot2,'(a)')'('
cFURUTA20191028 CSV output
                 elseif(itfoam(m).eq.2)then
                  if( itayl(m) .eq. 0 ) then
                   write(sbuf,'( "Number",a15)') hsunit(itunt(m))
                  else
                   write(sbuf,'(200a1)')
     &                  (itayt(m)(i:i),i=1,itayl(m))
                  end if
                  write(buf,'( "# tetra,xCM[cm],yCM[cm],zCM[cm],volume,
     &                  ",a200,",r.err")')sbuf
                  call remove_spaces(buf)
                  write(iot2,'(a)')trim(buf)

                 endif
               endif
*-----------------------------------------------------------------------

               if( nn .gt. 0 ) then

                  iz = nt(ic) / 1000

                  ia = nt(ic) - iz * 1000

                  if( ia .gt. 0 ) then

                     in = ia - iz

                     call chname(idum,ia,iz,chau)

                  end if

               end if

*-----------------------------------------------------------------------

                  write(iot,'(/"#",78("-"))')

               if( inum .eq. 1 ) then
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

               if( rtfac(m) .lt. 0.d0 )
     &           write(iot,
     &           '( "#restart maximum for normalization:",1p1e14.7)')
     &           facmax(m)

            if( nn .eq. 0 ) then

               write(iot,'("#   no. =",i3,3x,"Z = all")')
     &                     inum

            else if( ia .eq. 0 ) then

               write(iot,'("#   no. =",i3,3x,"Z = ",i3," : ",a3)')
     &                     inum, iz, elmnt(iz)

            else

               write(iot,'("#   no. =",i3,3x,a8)')
     &                     inum, chau

            end if

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: Serial Num. of Tetrahedron")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: Number ",a15)') hsunit(itunt(m))


               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

                  write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')

               if( itanl(m) .gt. 0 ) then

                  write(iot,'( "p: ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               end if

               if( itsans(m) .gt. 0 ) then

                  call write_sangel(iot,m,0)

               end if

c S.H. added IF statement below for epsout=2 (2017.2.18)
              if ( iteps(m) .ne. 2 ) then

               write(iot,'( "h:   x      n",5x,"n",10x,
     &                    "y,l3        n")')

              else

               write(iot,'( "h:   x      n",5x,"n",10x,
     &                    "y1,l3       ny2 dy1=[y1*y2]")')

              end if

               write(iot,'( "#    num    tetra   volume  ",
     &                    "   number      r.err")')

*-----------------------------------------------------------------------

                     seka = 0.0
                     sera = 0.0
                     seva = 0.0
                     voll = 0.0

               do ir = 1, nr

                     ireg = lr(ir)

                     sek = 0.0
                     ser = 0.0

                     if( itunt(m) .eq. 1 ) then

                        vm = 1.0

                     else if( itunt(m) .eq. 2 ) then

                        vm = vl(ir)

                     end if

                        voll = voll + vl(ir)

*-----------------------------------------------------------------------

                  if( nn .eq. 0 ) then

                     do kn = 1, mn
                     do kz = 1, mz
                     do il = 0, mm ! frtati 2022/02/18

                        lz  = ikzz( kz, kn )
                        ln  = iknn( kz, kn )

                        vn  = vm * tr(ir,lz,ln,il,1)
                        sek = sek + vn
                        ser = ser + ( vn * tr(ir,lz,ln,il,2) )**2

                     end do
                     end do
                     end do

                  else if( ia .eq. 0 ) then

                     do kn = 1, mn
                     do il = 0, mm ! frtati 2022/02/18

                        lz  = ikzz( iz, kn )
                        ln  = iknn( iz, kn )

                        vn  = vm * tr(ir,lz,ln,il,1)
                        sek = sek + vn
                        ser = ser + ( vn * tr(ir,lz,ln,il,2) )**2

                     end do
                     end do

                  else

                        lz  = ikzz( iz, in )
                        ln  = iknn( iz, in )

                     do il = 0, mm ! frtati 2022/02/18
                        vn  = vm * tr(ir,lz,ln,il,1)
                        sek = sek + vn
                        ser = ser
     &                      + ( vn * tr(ir,lz,ln,il,2) )**2
                     end do

                  end if

*-----------------------------------------------------------------------

                        seka = seka + sek
                        sera = sera + ser
                        seva = seva + vm

                     if( sek .ne. 0.0 ) then ! .gt. -> .ne. by T.Sato 2024/05/02 for negative sek

                        ser = sqrt( ser ) / sek
                        sek = sek / vm

                     end if

                  write(iot,'(i8,1x,i8,1pe13.4,1pe13.4,0pf8.4)')
     &                  ir, ireg, vl(ir), sek, ser
cFURUTA20190208 OpenFOAM output
                  if ( itfoam(m) .eq. 1) then
                   write(iot2,'(1pe13.4)') sek
cFURUTA20191028 CSV output
                  elseif(itfoam(m) .eq. 2) then
                   write(buf,'(i8,5(",",1pe13.4),",",0pf8.4)')
     &                  lr(ir), xcm(1:3,ir),vl(ir),
     &                  sek,ser
                   call remove_spaces(buf)
                   write(iot2,'(a)')trim(buf)
                  endif

               end do

                     if( seka .gt. 0.0 ) then

                        sera = sqrt( sera ) / seka

                        if( itunt(m) .eq. 2 )
     &                  seka = seka / seva

                     end if

              write(iot,'(/"#       sum over ",1pe13.4,1pe13.4,0pf8.4)')
     &                     voll, seka, sera

            if( itunt(m) .eq. 2 ) then

               write(iot,'(
     &         "#   vol sum = ",1p1e13.4,"  [cm^3]")') seva

            end if

               if( nn .eq. 0 ) then

                  write(iot,'(/a1,"no. =",i3,
     &                        ",    Z  =  all", a1)')
     &                        cha, inum, cha

               else if( ia .eq. 0 ) then

                  write(iot,'(/a1,"no. =",i3,
     &                        ",    Z  =  ",i3, a1)')
     &                        cha, inum, iz, cha

               else

                  write(iot,'(/a1,"no. =",i3,
     &                        ",    ",a8, a1)')
     &                        cha, inum, chau, cha

               end if

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

               if( nn .eq. 0 ) then

                  write(iot,'("wt: s(0.7)",/
     &                        a1,"vspace{-3}"/
     &                        "      Z &=&   all"/
     &                        "e:")')
     &                        yen

               else if( ia .eq. 0 ) then

                  write(iot,'("wt: s(0.7)",/
     &                        a1,"vspace{-3}"/
     &                        "      Z &=& ",i3,"   :  ",a3/
     &                        "e:")')
     &                        yen, iz, elmnt(iz)

               else

                  write(iot,'("wt: s(0.7)",/
     &                        a1,"vspace{-3}"/
     &                        "  nucl. &=&   ",a8/
     &                        "e:")')
     &                        yen, chau

               end if

*-----------------------------------------------------------------------
cFURUTA20190208 OpenFOAM output
             if ( itfoam(m) .eq. 1) then !FURUTA20191028
               write(iot2,'(a)')')'
               close(iot2)
             endif
*-----------------------------------------------------------------------

            end do

cFURUTA20190208 OpenFOAM output
            if ( itfoam(m) .ne. 0) then
              deallocate( foamfIType,foamfIndex )
            endif

*-----------------------------------------------------------------------
*        mass axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 1 ) then

               inum = 0

                  if( nn .eq. 0 ) then

                     nc = 1

                  else

                     nc = nn

                  end if

*-----------------------------------------------------------------------

            do ir = 1, nr

            do 260 ic = 1, nc

                  if( nn .gt. 0 ) iz = nt(ic) / 1000

                  do i = 1, maxnt + maxpt

                     tm(i,1) = 0.0d+0
                     tm(i,2) = 0.0d+0

                  end do

                  if( itunt(m) .eq. 1 ) then

                     vm = 1.0

                  else if( itunt(m) .eq. 2 ) then

                     vm = vl(ir)

                  end if

*-----------------------------------------------------------------------

                  if( nn .eq. 0 .and. itfln(m) .eq. 1 ) then

                     do ln = 1, mn
                     do il = 0, mm ! frtati 2022/02/18

                        lz = 1
                        im = ln

                        vn = vm * tr(ir,lz,ln,il,1)
                        tm(im,1) = tm(im,1) + vn
                        tm(im,2) = tm(im,2)
     &                           + ( vn * tr(ir,lz,ln,il,2) )**2

                     end do
                     end do

                  else if( nn .eq. 0 .and. itfln(m) .gt. 1 ) then

                     do kn = 1, mn
                     do kz = 1, mz
                     do il = 0, mm ! frtati 2022/02/18

                        lz = ikzz( kz, kn )
                        ln = iknn( kz, kn )

                        im = kz + kn

                        vn = vm * tr(ir,lz,ln,il,1)
                        tm(im,1) = tm(im,1) + vn
                        tm(im,2) = tm(im,2)
     &                           + ( vn * tr(ir,lz,ln,il,2) )**2

                     end do
                     end do
                     end do

                  else

                     do kn = 1, mn
                     do il = 0, mm ! frtati 2022/02/18

                        lz = ikzz( iz, kn )
                        ln = iknn( iz, kn )

                        im = iz + kn

                        vn = vm * tr(ir,lz,ln,il,1)
                        tm(im,1) = tm(im,1) + vn
                        tm(im,2) = tm(im,2)
     &                           + ( vn * tr(ir,lz,ln,il,2) )**2

                     end do
                     end do

                  end if

*-----------------------------------------------------------------------

               do i = 1, maxnt + maxpt

                  if( tm(i,1) .gt. 0.0d0 ) then

                     tm(i,2) = sqrt( tm(i,2) ) / tm(i,1)
                     tm(i,1) = tm(i,1) / vm

                  end if

               end do

               do i = 1, maxnt + maxpt

                  if( tm(i,1) .gt. 0.0d0 ) goto 200

               end do

               goto 260

  200          im = i

               do i = maxnt + maxpt, im + 1, -1

                  if( tm(i,1) .gt. 0.0d0 ) goto 210

               end do

               jm = im
  210          jm = i

               im = max( 1, im - 2 )
               jm = jm + 2

*-----------------------------------------------------------------------

               ireg = lr(ir)

                  write(iot,'(/"#",78("-"))')

               inum = inum + 1

               if( inum .eq. 1 ) then
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

               if( rtfac(m) .lt. 0.d0 )
     &           write(iot,
     &           '( "#restart maximum for normalization:",1p1e14.7)')
     &           facmax(m)

               if( nn .eq. 0 ) then

                  write(iot,'("#   no. =",i3,3x,
     &            "tetra=",i8,3x,"Z = all")')
     &                        inum, ireg

               else

                  write(iot,'("#   no. =",i3,3x,
     &            "tetra=",i8,3x,"Z = ",i3," : ",a3)')
     &                        inum, ireg, iz, elmnt(iz)

               end if

               write(iot,'("# im jm = ",i3,1x,i3)') im, jm

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: Mass")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: Number ",a15)') hsunit(itunt(m))

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

                  write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')

               if( itanl(m) .gt. 0 ) then

                  write(iot,'( "p: ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               end if

               if( itsans(m) .gt. 0 ) then

                  call write_sangel(iot,m,0)

               end if

c S.H. added IF statement below for epsout=2 (2017.2.18)
              if ( iteps(m) .ne. 2 ) then

               write(iot,'( "h: x-0.5    y,hl0       n")')

              else

               write(iot,'( "h: x-0.5    y1,hl0      ",
     &                      "ny2 dy1=[y1*y2]")')

              end if

               write(iot,'( "#  mass     number      r.err")')

*-----------------------------------------------------------------------

               if( im .eq. 1 ) then

                  write(iot,'(3x,f4.1,2x,1pe13.4,0pf8.4)')
     &                                   0.5, 0.0, 0.0

               end if

                  seka = 0.0
                  sera = 0.0

               do i = im, jm

                  vn  = vm * tm(i,1)
                  seka = seka + vn
                  sera = sera + ( vn * tm(i,2) )**2

                  write(iot,'(3x,i4,2x,1pe13.4,0pf8.4)')
     &                           i, tm(i,1), tm(i,2)

               end do

                  if( seka .gt. 0.0 ) then

                     sera = sqrt( sera ) / seka
                     seka = seka / vm

                  end if

               write(iot,'(/"#   sum  ",1pe13.4,0pf8.4)')
     &               seka, sera

               if( itunt(m) .eq. 2 ) then

                  write(iot,'(
     &            "#   vol sum = ",1p1e13.4,"  [cm^3]")') vm

               end if

               if( nn .eq. 0 ) then

                  write(iot,'(/a1,"no. =",i3,
     &                            ",  tetra  =",i8,
     &                            ",    Z  =  all", a1)')
     &                        cha, inum, ireg, cha

               else

                  write(iot,'(/a1,"no. =",i3,
     &                            ",  tetra  =",i8,
     &                            ",    Z  =  ",i3, a1)')
     &                        cha, inum, ireg, iz, cha

               end if

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

               if( itunt(m) .eq. 1 ) then

                  if( nn .eq. 0 ) then

                     write(iot,'("wt: s(0.7)",/
     &                           a1,"vspace{-3}"/
     &                           "  mass distribution"/
     &                           "  in tetra mesh"/
     &                           "    Z &=&   all"/
     &                           "e:")')
     &                           yen

                  else

                     write(iot,'("wt: s(0.7)",/
     &                           a1,"vspace{-3}"/
     &                           "  isotope distribution"/
     &                           "  in tetra mesh"/
     &                           "    Z &=& ",i3,"   :  ",a3/
     &                           "e:")')
     &                           yen, iz, elmnt(iz)

                  end if

               else if( itunt(m) .eq. 2 ) then

                  if( nn .eq. 0 ) then

                     write(iot,'("wt: s(0.7)",/
     &                           a1,"vspace{-3}"/
     &                           "  mass distribution"/
     &                           "  in tetra mesh"/
     &                           "  vol &=&",1pe13.4," [cm^3]"/
     &                           "    Z &=&   all"/
     &                           "e:")')
     &                           yen, vl(ir)

                  else

                     write(iot,'("wt: s(0.7)",/
     &                           a1,"vspace{-3}"/
     &                           "  isotope distribution"/
     &                           "  in tetra mesh"/
     &                           "  vol &=&",1pe13.4," [cm^3]"/
     &                           "    Z &=& ",i3,"   :  ",a3/
     &                           "e:")')
     &                           yen, vl(ir),
     &                           iz, elmnt(iz)

                  end if

               end if

  260       continue

            end do

*-----------------------------------------------------------------------
*        charge axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 7 ) then

               inum = 0

            do 370 ir = 1, nr

                  do i = 1, maxpt

                     tm(i,1) = 0.0d+0
                     tm(i,2) = 0.0d+0

                  end do

                  if( itunt(m) .eq. 1 ) then

                     vm = 1.0

                  else if( itunt(m) .eq. 2 ) then

                     vm = vl(ir)

                  end if

*-----------------------------------------------------------------------

                  do kz = 1, mz
                  do kn = 1, mn
                  do il = 0, mm ! frtati 2022/02/18

                        iz  = ikzz( kz, kn )
                        in  = iknn( kz, kn )

                        vn  = vm * tr(ir,iz,in,il,1)
                        tm(iz,1) = tm(iz,1) + vn
                        tm(iz,2) = tm(iz,2)
     &                           + ( vn * tr(ir,iz,in,il,2) )**2

                  end do
                  end do
                  end do

*-----------------------------------------------------------------------

               do i = 1, maxpt

                  if( tm(i,1) .gt. 0.0d0 ) then

                     tm(i,2) = sqrt( tm(i,2) ) / tm(i,1)
                     tm(i,1) = tm(i,1) / vm

                  end if

               end do

               do i = 1, maxpt

                  if( tm(i,1) .gt. 0.0d0 ) goto 300

               end do

               goto 370

  300          im = i

               do i = maxpt, im + 1, -1

                  if( tm(i,1) .gt. 0.0d0 ) goto 310

               end do

               jm = im
  310          jm = i

               im = max( 1, im - 2 )
               jm = jm + 2

*-----------------------------------------------------------------------

               ireg = lr(ir)

                  write(iot,'(/"#",78("-"))')

               inum = inum + 1

               if( inum .eq. 1 ) then
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

               if( rtfac(m) .lt. 0.d0 )
     &           write(iot,
     &           '( "#restart maximum for normalization:",1p1e14.7)')
     &           facmax(m)

               write(iot,'("#   no. =",i3,3x,"tetra=",i8)')
     &                     inum, ireg

               write(iot,'("# im jm = ",i3,1x,i3)') im, jm

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: Charge")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: Number ",a15)') hsunit(itunt(m))

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

                  write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')

               if( itanl(m) .gt. 0 ) then

                  write(iot,'( "p: ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               end if

               if( itsans(m) .gt. 0 ) then

                  call write_sangel(iot,m,0)

               end if

c S.H. added IF statement below for epsout=2 (2017.2.18)
              if ( iteps(m) .ne. 2 ) then

               write(iot,'( "h: x-0.5    y,hl0       n")')

              else

               write(iot,'( "h: x-0.5    y1,hl0      ",
     &                      "ny2 dy1=[y1*y2]")')

              end if

               write(iot,'( "#  charge   number      r.err")')

*-----------------------------------------------------------------------

               if( im .eq. 1 ) then

                  write(iot,'(3x,f4.1,2x,1pe13.4,0pf8.4)')
     &                                   0.5, 0.0, 0.0

               end if

                  seka = 0.0
                  sera = 0.0

               do i = im, jm

                  vn  = vm * tm(i,1)
                  seka = seka + vn
                  sera = sera + ( vn * tm(i,2) )**2

                  write(iot,'(3x,i4,2x,1pe13.4,0pf8.4)')
     &                           i, tm(i,1), tm(i,2)

               end do

                  if( seka .gt. 0.0 ) then

                     sera = sqrt( sera ) / seka
                     seka = seka / vm

                  end if

               write(iot,'(/"#   sum  ",1pe13.4,0pf8.4)')
     &               seka, sera

               if( itunt(m) .eq. 2 ) then

                  write(iot,'(
     &            "#   vol sum = ",1p1e13.4,"  [cm^3]")') vm

               end if

                  write(iot,'(/a1,"no. =",i3,
     &                            ",  tetra  =",i8,a1)')
     &                        cha, inum, ireg, cha

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

               if( itunt(m) .eq. 1 ) then

                     write(iot,'("wt: s(0.7)",/
     &                           a1,"vspace{-3}"/
     &                           "  charge distribution"/
     &                           "  in tetra mesh"/
     &                           "e:")')
     &                           yen

               else if( itunt(m) .eq. 2 ) then

                     write(iot,'("wt: s(0.7)",/
     &                           a1,"vspace{-3}"/
     &                           "  charge distribution"/
     &                           "  in tetra mesh"/
     &                           "  vol &=&",1pe13.4," [cm^3]"/
     &                           "e:")')
     &                           yen, vl(ir)

               end if

  370       continue

*-----------------------------------------------------------------------
*        chart axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 8 ) then

               inum = 0

            do 470 ir = 1, nr

            do 460 il = 0, 2 ! loop for isomeric level
*-----------------------------------------------------------------------

                  itmax  = 0
                  icmax  = 0
                  inmax  = 0

               do iz = 1, maxpt
               do in = 1, maxnt

                  if( itnzn(m).eq.0 ) then
                    t0_tr = tr(ir,iz,in,il,1)
                  else
                    t0_tr = tr(ir,igetiznmp(iz,in,il,m),1,0,1)
                  end if
                  if( t0_tr .gt. 0.d0 ) then

                     if( iz + in .gt. itmax ) itmax = iz + in
                     if( iz      .gt. icmax ) icmax = iz
                     if(      in .gt. inmax ) inmax =      in

                  end if

               end do
               end do

                     dxmax = dble(inmax+2)
                     dymax = dble(icmax+2)
                     dform = dymax / dxmax

                     inmag = inmax

                  if( dform .gt. 0.8 ) then

                     inmag = nint( dymax / 0.8 )
                     dxmax = dble( nint( dymax / 0.8 ) )
                     dform = dymax / dxmax

                  end if

               if( itmax .eq. 0 ) goto 470

*-----------------------------------------------------------------------

               ireg = lr(ir)

                  write(iot,'(/"#",78("-"))')

               inum = inum + 1

               if( inum .eq. 1 .and. il .eq. 0 ) then
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

               if( rtfac(m) .lt. 0.d0 )
     &           write(iot,
     &           '( "#restart maximum for normalization:",1p1e14.7)')
     &           facmax(m)

               write(iot,'("#   no. =",i3,3x,
     &         "ir  =",i3,3x,"il  =",i3)' )
     &                     inum, ir, il   ! tally data number, tetra number, isomeric level

               write(iot,'("# icmax inmax = ",i3,1x,i3)')
     &            icmax, inmax

*-----------------------------------------------------------------------

            if( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) then

                  write(iot,'(/a1,"no. =",i3,
     &                            ",  tetra  =",i8,a1)')
     &                        cha, inum, ireg, cha

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

             if(il .eq. 0) then    ! specify isomeric level at the bottom of figure
               write(iot,'("msdc: {",a1,"huge Ground state}")') yen
             elseif(il .eq. 1) then
               write(iot,'("msdc: {",a1,"huge 1st isomer}")') yen
             else
               write(iot,'("msdc: {",a1,"huge 2nd isomer}")') yen
             endif

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: N Neutron Number")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: Z Proton Number")')

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

*-----------------------------------------------------------------------

               if( inum .eq. 1 ) then

                  form  = dform
                  xfac  = 1.1
                  afac  = 0.6
                  izlog = 1
                  inocm = 1
                  inolg = 1

                  if( itanl(m) .gt. 0 )
     &            call anset(itang(m),itanl(m),
     &                       form,xfac,afac,cmin,cmax,izlog,inocm,inolg)

                  if( form .le. 1.0 ) then

                     scal = form**0.35
                     xfac = xfac / form**0.5
                     xorg = 0.0
                     yorg = min( 20.d0, ( 1.0 / form - 1.0 ) / 2.5 )

                  else

                     scal = 1.0 / form**0.41
                     xfac = xfac / form**0.5
                     xorg = min( 20.d0, ( form - 1.0 ) / 1.5)
                     yorg = 0.0

                  end if

               end if

                  write(iot,'( "set: c1[",f6.3,
     &                            "] c2[",f6.3,
     &                            "] c5[",f6.3,
     &                            "] c6[",f6.3,
     &                            "] c7[",f6.3,
     &                            "] c8[",f6.3,"]")')
     &                    form, xfac, afac, scal, xorg, yorg
                  write(iot,'( "p: h2fs form[c1] xfac[c2]",
     &                 " afac[c5] scal[c6] xorg[c7] yorg[c8] nosp")')

               if( ( ( ittwo(m) .ge. 2 .and. ittwo(m) .le. 3 ) .or.
     &               ( ittwo(m) .ge. 6 .and. ittwo(m) .le. 7 ) ) .and.
     &             cmin .gt. 0.0 .and. cmax .gt. cmin ) then

                if (ioe .eq. 1 ) then
                  write(iot,'( "set: c3[",1p1e15.6,"] c4[",
     &                  1p1e15.6,"]")') cmin, cmax
                else
                  write(iot,'( "set: c3[1.0e-4] c4[1.0]")')
                end if
                  write(iot,'( "p: cmin[c3] cmax[c4]")')
                  write(iot,'( "p: dmin(1e-31)")')

                  if( ittwo(m) .eq. 3 .or. ittwo(m) .eq. 7 )
     &                 write(iot,'( "p: zlog")')

               end if

               write(iot,'( "p: xmin(0) xmax(",f5.1,")")') dxmax
               write(iot,'( "p: ymin(0) ymax(",f5.1,")")') dymax

               if( itanl(m) .gt. 0 .and. ioe .eq. 1 ) then

                  write(iot,'( "p: ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               else if( itanl(m) .gt. 0 .and. ioe .ne. 1 ) then

                   call terrang(iot,m)

               end if

               if( itsans(m) .gt. 0 ) then

                  call write_sangel(iot,m,0)

               end if

            end if

*-----------------------------------------------------------------------

            if( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) then

               write(iot,'(/"hc: y = ",i3," to 1 by -1 ;",
     &                         " x = 1 to ",i3," by 1 ;")')
     &                        icmax+2, inmax+2

               do i = icmax+2, 1, -1
                if( itnzn(m).eq.0 ) then ! frtati 2022/02/18
                  write(iot,'(1p10e11.3)')
     &                 ( tr(ir,i,l,il,ioe), l = 1, inmax+2 )
                else
                  write(iot,'(1p10e11.3)')
     &          ( tr(ir,igetiznmp(i,l,il,m),1,0,ioe), l = 1, inmax+2 )
                end if

               end do

*-----------------------------------------------------------------------

            else if( ittwo(m) .eq. 4 ) then

               write(iot,'(/"#   Z    N  Mass  ",
     &                      "  number    r.err")')

               do i = 1, icmax
               do j = 1, inmax

                if( itnzn(m).eq.0 ) then ! frtati 2022/02/18
                  if( tr(ir,i,j,il,1) .gt. 0.d0 ) then

                     write(iot,'(3i5,1pe13.4,0pf8.4)')
     &               i, j, i+j,
     &               tr(ir,i,j,il,1), tr(ir,i,j,il,2)

                  end if
                else
                  if( tr(ir,igetiznmp(i,j,il,m),1,0,1) .gt. 0.d0 ) then
                     write(iot,'(3i5,1pe13.4,0pf8.4)')
     &               i, j, i+j,
     &               tr(ir,igetiznmp(i,j,il,m),1,0,1),
     &               tr(ir,igetiznmp(i,j,il,m),1,0,2)
                  end if
                end if

               end do
               end do

*-----------------------------------------------------------------------

            else if( ittwo(m) .eq. 5 ) then

               write(iot,'(/"#   Z = 1 to ",i3/
     &                      "#   N = 1 to ",i3)')
     &                        icmax+2, inmax+2

                  write(iot,'(/4x,a3,4x,1p1000e11.3)')
     &            'Z/N', ( dble( in ), in = 1, inmax+2 )

               do i = icmax+2, 1, -1
                if( itnzn(m).eq.0 ) then ! frtati 2022/02/18
                  write(iot,'(1p1000e11.3)')
     &            dble( i ),
     &            ( tr(ir,i,l,il,ioe), l = 1, inmax+2 )
                else
                  write(iot,'(1p1000e11.3)')
     &            dble( i ),
     &            ( tr(ir,igetiznmp(i,l,il,m),1,0,ioe),l = 1,inmax+2 )
                end if

               end do

            end if

*-----------------------------------------------------------------------

            if( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) then

               if( itout(m) .ne. 0 ) then

                  call wmgcstb(iot,icmax,inmag)

               end if

            end if

*-----------------------------------------------------------------------

      if( inolg .eq. 1 .and. ( ittwo(m) .eq. 3 .or. ittwo(m) .eq. 7 )
     &    .and. cmin .gt. 0.0 .and. cmax .gt. cmin .and.
     &    igsh .eq. 0 ) then

      write(iot,'(/"#",78("-"))')
      write(iot,'(
     &"z: xorg(1.03)"/
     &"p: xfac[c2*0.03] form[c1/0.03] notc noxt noyt dmin(1)"/
     &"p: ymin(0) ymax(1) xmin(0) xmax(1) cmin(1) cmax(100) nosp"/
     &"hc: y= 0.005 to 0.995 by 0.01 ; x= 0.5 to 0.5 by 1 ;"/
     &" 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20"/
     &"21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40"/
     &"41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60"/
     &"61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80"/
     &"81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100"/
     &"z: xorg(0.0)"/
     &"p: xfac[c2*0.03] form[c1/0.03] ",
     &"nosp afac[c5*0.7] ymin[c3] ymax[c4]"/
     &"p: noxt noxn ytxt(-1) ynum(-1) xltd(-1) xstd(-1) itic(1) ylog"
     &)')

      if ( ioe .eq. 2 ) then

          write(iot,'("y: Relative Error")')

      else if( itazl(m) .eq. 0 ) then

         write(iot,'("y: Number ",a15)') hsunit(itunt(m))

      else

         write(iot,'("y: ",200a1)') (itazt(m)(i:i),i=1,itazl(m))

      end if

      if( izlog .eq. 0 ) write(iot,'( "p: ylin")')

      write(iot,'( "#",78("-"))')

      end if

*-----------------------------------------------------------------------

            if( inocm .eq. 1 ) then

               if( itunt(m) .eq. 1 ) then

                     write(iot,'(/"wt: s(0.7)",/
     &                           a1,"vspace{-3}"/
     &                           "  yield distribution"/
     &                           "  in tetra mesh"/
     &                           "e:")')
     &                           yen

               else if( itunt(m) .eq. 2 ) then

                     write(iot,'(/"wt: s(0.7)",/
     &                           a1,"vspace{-3}"/
     &                           "  yield distribution"/
     &                           "  in tetra mesh"/
     &                           "  vol &=&",1pe13.4," [cm^3]"/
     &                           "e:")')
     &                           yen, vl(ir)

               end if

            end if

*-----------------------------------------------------------------------

      if( inolg .eq. 1 .and. ( ittwo(m) .eq. 3 .or. ittwo(m) .eq. 7 )
     &    .and. cmin .gt. 0.0 .and. cmax .gt. cmin .and.
     &    igsh .eq. 0 ) then

         write(iot,'(/"#",78("-"))')
         write(iot,'("z: xorg[-1.03/0.05]"/
     &               "p: form[c1] xfac[c2] afac[c5] nosp notf")')

               write(iot,'( "p: xmin(0) xmax(",f5.1,")")') dxmax
               write(iot,'( "p: ymin(0) ymax(",f5.1,")")') dymax

      end if

*-----------------------------------------------------------------------

  460       continue

  470       continue

*-----------------------------------------------------------------------
*        xy axis (matrix)
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 9 ) then

*-----------------------------------------------------------------------

               inum = 0

                  if( nn .eq. 0 .or. igsh .ne. 0  ) then

                     nc = 1

                  else

                     nc = nn

                  end if

*-----------------------------------------------------------------------

            do jz = 1, nz
            do ic = 1, nc

               zval = ( zm(jz) + zm(jz+1) ) / 2.0d0

               inum = inum + 1

*-----------------------------------------------------------------------

               if( nn .gt. 0 ) then

                  iz = nt(ic) / 1000

                  ia = nt(ic) - iz * 1000

                  if( ia .gt. 0 ) then

                     in = ia - iz

                     call chname(idum,ia,iz,chau)

                  end if

               end if

*-----------------------------------------------------------------------

                  write(iot,'(/"#",78("-"))')

               if( inum .eq. 1 ) then
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

            if( nn .eq. 0 ) then

               write(iot,'("#   no. =",i3,3x,
     &         "iz  =",i3,3x,"Z = all")')
     &                     inum, jz

            else if( ia .eq. 0 ) then

               write(iot,'("#   no. =",i3,3x,
     &         "iz  =",i3,3x,"Z = ",i3," : ",a3)')
     &                     inum, jz, iz, elmnt(iz)

            else

               write(iot,'("#   no. =",i3,3x,
     &         "iz  =",i3,3x,a8)')
     &                     inum, jz, chau

            end if

*-----------------------------------------------------------------------

               if( nn .eq. 0 ) then

                  write(iot,'(/a1,"no. =",i3,
     &                        ",    iz =",i3,
     &                        ",    Z  =  all", a1)')
     &                        cha, inum, jz, cha

               else if( ia .eq. 0 ) then

                  write(iot,'(/a1,"no. =",i3,
     &                        ",    iz =",i3,
     &                        ",    Z  =  ",i3, a1)')
     &                        cha, inum, jz, iz, cha

               else

                  write(iot,'(/a1,"no. =",i3,
     &                        ",    iz =",i3,
     &                        ",    ",a8, a1)')
     &                        cha, inum, jz, chau, cha

               end if

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: x [cm]")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: y [cm]")')

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

*-----------------------------------------------------------------------

         if( igsh .eq. 0 ) then

               do i = 1, nr

                     sek = 0.0

                  if( nn .eq. 0 ) then

                     do kn = 1, mn
                     do kz = 1, mz
                     do il = 0, mm ! frtati 2022/02/18

                        lz  = ikzz( kz, kn )
                        ln  = iknn( kz, kn )

                        sek = sek + tr(i,lz,ln,il,1)

                     end do
                     end do
                     end do

                  else if( ia .eq. 0 ) then

                     do kn = 1, mn
                     do il = 0, mm ! frtati 2022/02/18

                        lz  = ikzz( iz, kn )
                        ln  = iknn( iz, kn )

                        sek = sek + tr(i,lz,ln,il,1)

                     end do
                     end do

                  else

                        lz  = ikzz( iz, in )
                        ln  = iknn( iz, in )

                     do il = 0, mm ! frtati 2022/02/18
                        sek = sek + tr(i,lz,ln,il,1)
                     end do

                  end if

                  val(i) = sek

               end do

                  izlog = 1
                  if( itanl(m) .gt. 0 )
     &            call anset(itang(m),itanl(m),
     &                       form,xfac,afac,cmin,cmax,izlog,inocm,inolg)

            if( izlog .eq. 0 ) then

               do i = 1, nr

                  if( val(i) .gt. dnon .and. cmin .ne. cmax ) then

                     val(i) = ( val(i) - cmin ) / ( cmax - cmin )
                     val(i) = max(0.0d0,val(i))
                     val(i) = min(1.0d0,val(i))

                  else

                     val(i) = -1.0

                  end if

               end do

            else

               do i = 1, nr

                  if( val(i) .gt. dnon .and. cmin .ne. cmax ) then

                     val(i) = log10(val(i)/cmin) / log10(cmax/cmin)
                     val(i) = max(0.0d0,val(i))
                     val(i) = min(1.0d0,val(i))

                  else

                     val(i) = -1.0d0

                  end if

               end do

            end if

         end if

*-----------------------------------------------------------------------

               if( inum .eq. 1 ) then

                  form  = ( ym(ny+1) - ym(1) )
     &                  / ( xm(nx+1) - xm(1) )
                  xfac  = 0.9
                  afac  = 0.8
                  izlog = 1
                  inocm = 1
                  inolg = 1

                  if( itanl(m) .gt. 0 )
     &            call anset(itang(m),itanl(m),
     &                       form,xfac,afac,cmin,cmax,izlog,inocm,inolg)

                  if( form .le. 1.0 ) then

                     scal = form**0.35
                     xfac = xfac / form**0.5
                     xorg = 0.0
                     yorg = min( 20.d0, ( 1.0 / form - 1.0 ) / 2.5 )

                  else

                     scal = 1.0 / form**0.41
                     xfac = xfac / form**0.5
                     xorg = min( 20.d0, ( form - 1.0 ) / 1.5)
                     yorg = 0.0

                  end if

               end if

                  write(iot,'( "set: c1[",f6.3,
     &                            "] c2[",f6.3,
     &                            "] c5[",f6.3,
     &                            "] c6[",f6.3,
     &                            "] c7[",f6.3,
     &                            "] c8[",f6.3,"]")')
     &                    form, xfac, afac, scal, xorg, yorg
                  write(iot,'( "p: h2fs form[c1] xfac[c2]",
     &                 " afac[c5] scal[c6] xorg[c7] yorg[c8] nosp")')

               if( cmin .gt. 0.0 .and. cmax .gt. cmin ) then

                  write(iot,'( "set: c3[",1p1e15.6,"] c4[",
     &                  1p1e15.6,"]")') cmin, cmax

               end if

*-----------------------------------------------------------------------

               xmin = xm(1)
               xmax = xm(nx+1)
               ymin = ym(1)
               ymax = ym(ny+1)

               write(iot,'( "p: xmin(",1p1g14.7,") xmax(",1p1g14.7,
     &                       ")")') xmin, xmax

               write(iot,'( "p: ymin(",1p1g14.7,") ymax(",1p1g14.7,
     &                       ")")') ymin, ymax

               if( itanl(m) .gt. 0 ) then

                  write(iot,'( "p: ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               end if

               if( itsans(m) .gt. 0 ) then

                  call write_sangel(iot,m,0)

               end if

*-----------------------------------------------------------------------
*        rshow
*-----------------------------------------------------------------------

               write(iot,'(/"#",78("-"))')
               write(iot,'( "# rshow")')
               write(iot,'( "#",78("-"))')
               write(iot,'(/"p: legs[c5*0.875]")')

               none = 1
               iaxs = 1
               iuni = itrsh(m) * 2 - 1 + igsh
               ires = itres(m)
               widt = rtwid(m)

            if( igsh .eq. 0 ) then

               call gshow(2,iot,iaxs,iuni,widt,ires,
     &                    inum,igsh,igser,gfnam,igfmn,
     &                    nx+1,ny+1,none,xm,ym,zval,ixyz(1),
     &                    nr,1,krr,val,itmtr(m,4),itglt(m))

            else

               call gshow(0,iot,iaxs,iuni,widt,ires,
     &                    inum,igsh,igser,gfnam,igfmn,
     &                    nx+1,ny+1,none,xm,ym,zval,ixyz(1),
     &                    1,1,krr,vll,itmtr(m,4),itglt(m))

            end if

*-----------------------------------------------------------------------

      if( inolg .eq. 1 .and. igsh .eq. 0 .and.
     &    cmin .gt. 0.0 .and. cmax .gt. cmin ) then

      write(iot,'(/"#",78("-"))')
      write(iot,'(
     &"z: xorg(1.03)"/
     &"p: xfac[c2*0.05] form[c1/0.05] notc noxt noyt dmin(1)"/
     &"p: ymin(0) ymax(1) xmin(0) xmax(1) cmin(1) cmax(100) nosp"/
     &"hc: y= 0.005 to 0.995 by 0.01 ; x= 0.5 to 0.5 by 1 ;"/
     &" 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20"/
     &"21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40"/
     &"41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60"/
     &"61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80"/
     &"81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100"/
     &"z: xorg(0.0)"/
     &"p: xfac[c2*0.05] form[c1/0.05] ",
     &"nosp afac[c5*0.625] ymin[c3] ymax[c4]"/
     &"p: noxt noxn ytxt(-1) ynum(-1) xltd(-1) xstd(-1) itic(1) ylog"
     &)')

      if( itazl(m) .eq. 0 ) then

         write(iot,'("y: Number ",a15)') hsunit(itunt(m))

      else

         write(iot,'("y: ",200a1)') (itazt(m)(i:i),i=1,itazl(m))

      end if

      if( izlog .eq. 0 ) write(iot,'( "p: ylin")')

      write(iot,'( "#",78("-"))')

      end if

*-----------------------------------------------------------------------

            if( inocm .eq. 1 ) then

               if( nn .eq. 0 ) then

                  write(iot,'(/"wt: s[c5]",/
     &                        a1,"vspace{-3}"/
     &                        "  z &=&",1pe13.4," [cm]"/
     &                        "  Z &=&   all"/
     &                        "e:")')
     &                        yen, zval

               else if( ia .eq. 0 ) then

                  write(iot,'(/"wt: s[c5]",/
     &                        a1,"vspace{-3}"/
     &                        "  z &=&",1pe13.4," [cm]"/
     &                        "  Z &=& ",i3,"   :  ",a3/
     &                        "e:")')
     &                        yen,
     &                        zval, iz, elmnt(iz)

               else

                  write(iot,'(/"wt: s[c5]",/
     &                        a1,"vspace{-3}"/
     &                        "  z     &=&",1pe13.4," [cm]"/
     &                        "  nucl. &=&   ",a8/
     &                        "e:")')
     &                        yen,
     &                        zval, chau

               end if

            end if

*-----------------------------------------------------------------------

      if( inolg .eq. 1 .and. igsh .eq. 0 .and.
     &    cmin .gt. 0.0 .and. cmax .gt. cmin ) then

         write(iot,'(/"#",78("-"))')
         write(iot,'("z: xorg[-1.03/0.05]"/
     &               "p: form[c1] xfac[c2] afac[c5] nosp notf")')

               write(iot,'( "p: xmin(",1p1g14.7,") xmax(",1p1g14.7,
     &                       ")")') xmin, xmax

               write(iot,'( "p: ymin(",1p1g14.7,") ymax(",1p1g14.7,
     &                       ")")') ymin, ymax

      end if

*-----------------------------------------------------------------------

            end do
            end do

*-----------------------------------------------------------------------
*        yz axis (matrix)
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 10 ) then

*-----------------------------------------------------------------------

               inum = 0

                  if( nn .eq. 0 .or. igsh .ne. 0  ) then

                     nc = 1

                  else

                     nc = nn

                  end if

*-----------------------------------------------------------------------

            do jx = 1, nx
            do ic = 1, nc

               xval = ( xm(jx) + xm(jx+1) ) / 2.0d0

               inum = inum + 1

*-----------------------------------------------------------------------

               if( nn .gt. 0 ) then

                  iz = nt(ic) / 1000

                  ia = nt(ic) - iz * 1000

                  if( ia .gt. 0 ) then

                     in = ia - iz

                     call chname(idum,ia,iz,chau)

                  end if

               end if

*-----------------------------------------------------------------------

                  write(iot,'(/"#",78("-"))')

               if( inum .eq. 1 ) then
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

               if( rtfac(m) .lt. 0.d0 )
     &           write(iot,
     &           '( "#restart maximum for normalization:",1p1e14.7)')
     &           facmax(m)

            if( nn .eq. 0 ) then

               write(iot,'("#   no. =",i3,3x,
     &         "ix  =",i3,3x,"Z = all")')
     &                     inum, jx

            else if( ia .eq. 0 ) then

               write(iot,'("#   no. =",i3,3x,
     &         "ix  =",i3,3x,"Z = ",i3," : ",a3)')
     &                     inum, jx, iz, elmnt(iz)

            else

               write(iot,'("#   no. =",i3,3x,
     &         "ix  =",i3,3x,a8)')
     &                     inum, jx, chau

            end if

*-----------------------------------------------------------------------

               if( nn .eq. 0 ) then

                  write(iot,'(/a1,"no. =",i3,
     &                        ",    ix =",i3,
     &                        ",    Z  =  all", a1)')
     &                        cha, inum, jx, cha

               else if( ia .eq. 0 ) then

                  write(iot,'(/a1,"no. =",i3,
     &                        ",    ix =",i3,
     &                        ",    Z  =  ",i3, a1)')
     &                        cha, inum, jx, iz, cha

               else

                  write(iot,'(/a1,"no. =",i3,
     &                        ",    ix =",i3,
     &                        ",    ",a8, a1)')
     &                        cha, inum, jx, chau, cha

               end if

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: z [cm]")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: y [cm]")')

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

*-----------------------------------------------------------------------

         if( igsh .eq. 0 ) then

               do i = 1, nr

                     sek = 0.0

                  if( nn .eq. 0 ) then

                     do kn = 1, mn
                     do kz = 1, mz
                     do il = 0, mm ! frtati 2022/02/18

                        lz  = ikzz( kz, kn )
                        ln  = iknn( kz, kn )

                        sek = sek + tr(i,lz,ln,il,1)

                     end do
                     end do
                     end do

                  else if( ia .eq. 0 ) then

                     do kn = 1, mn
                     do il = 0, mm ! frtati 2022/02/18

                        lz  = ikzz( iz, kn )
                        ln  = iknn( iz, kn )

                        sek = sek + tr(i,lz,ln,il,1)

                     end do
                     end do

                  else

                        lz  = ikzz( iz, in )
                        ln  = iknn( iz, in )

                     do il = 0, mm ! frtati 2022/02/18
                        sek = sek + tr(i,lz,ln,il,1)
                     end do

                  end if

                  val(i) = sek

               end do

                  izlog = 1
                  if( itanl(m) .gt. 0 )
     &            call anset(itang(m),itanl(m),
     &                       form,xfac,afac,cmin,cmax,izlog,inocm,inolg)

            if( izlog .eq. 0 ) then

               do i = 1, nr

                  if( val(i) .gt. dnon .and. cmin .ne. cmax ) then

                     val(i) = ( val(i) - cmin ) / ( cmax - cmin )
                     val(i) = max(0.0d0,val(i))
                     val(i) = min(1.0d0,val(i))

                  else

                     val(i) = -1.0

                  end if

               end do

            else

               do i = 1, nr

                  if( val(i) .gt. dnon .and. cmin .ne. cmax ) then

                     val(i) = log10(val(i)/cmin) / log10(cmax/cmin)
                     val(i) = max(0.0d0,val(i))
                     val(i) = min(1.0d0,val(i))

                  else

                     val(i) = -1.0d0

                  end if

               end do

            end if

         end if

*-----------------------------------------------------------------------

               if( inum .eq. 1 ) then

                  form  = ( ym(ny+1) - ym(1) )
     &                  / ( zm(nz+1) - zm(1) )
                  xfac  = 0.9
                  afac  = 0.8
                  izlog = 1
                  inocm = 1
                  inolg = 1

                  if( itanl(m) .gt. 0 )
     &            call anset(itang(m),itanl(m),
     &                       form,xfac,afac,cmin,cmax,izlog,inocm,inolg)

                  if( form .le. 1.0 ) then

                     scal = form**0.35
                     xfac = xfac / form**0.5
                     xorg = 0.0
                     yorg = min( 20.d0, ( 1.0 / form - 1.0 ) / 2.5 )

                  else

                     scal = 1.0 / form**0.41
                     xfac = xfac / form**0.5
                     xorg = min( 20.d0, ( form - 1.0 ) / 1.5)
                     yorg = 0.0

                  end if

               end if

                  write(iot,'( "set: c1[",f6.3,
     &                            "] c2[",f6.3,
     &                            "] c5[",f6.3,
     &                            "] c6[",f6.3,
     &                            "] c7[",f6.3,
     &                            "] c8[",f6.3,"]")')
     &                    form, xfac, afac, scal, xorg, yorg
                  write(iot,'( "p: h2fs form[c1] xfac[c2]",
     &                 " afac[c5] scal[c6] xorg[c7] yorg[c8] nosp")')

               if( cmin .gt. 0.0 .and. cmax .gt. cmin ) then

                  write(iot,'( "set: c3[",1p1e15.6,"] c4[",
     &                  1p1e15.6,"]")') cmin, cmax

               end if

*-----------------------------------------------------------------------

               zmin = zm(1)
               zmax = zm(nz+1)
               ymin = ym(1)
               ymax = ym(ny+1)

               write(iot,'( "p: xmin(",1p1g14.7,") xmax(",1p1g14.7,
     &                       ")")') zmin, zmax

               write(iot,'( "p: ymin(",1p1g14.7,") ymax(",1p1g14.7,
     &                       ")")') ymin, ymax

               if( itanl(m) .gt. 0 ) then

                  write(iot,'( "p: ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               end if

               if( itsans(m) .gt. 0 ) then

                  call write_sangel(iot,m,0)

               end if

*-----------------------------------------------------------------------
*        rshow
*-----------------------------------------------------------------------

               write(iot,'(/"#",78("-"))')
               write(iot,'( "# rshow")')
               write(iot,'( "#",78("-"))')
               write(iot,'(/"p: legs[c5*0.875]")')

               none = 1
               iaxs = 2
               iuni = itrsh(m) * 2 - 1 + igsh
               ires = itres(m)
               widt = rtwid(m)

            if( igsh .eq. 0 ) then

               call gshow(2,iot,iaxs,iuni,widt,ires,
     &                    inum,igsh,igser,gfnam,igfmn,
     &                    nz+1,ny+1,none,zm,ym,xval,ixyz(1),
     &                    nr,1,krr,val,itmtr(m,4),itglt(m))

            else

               call gshow(0,iot,iaxs,iuni,widt,ires,
     &                    inum,igsh,igser,gfnam,igfmn,
     &                    nz+1,ny+1,none,zm,ym,xval,ixyz(1),
     &                    1,1,krr,vll,itmtr(m,4),itglt(m))

            end if

*-----------------------------------------------------------------------

      if( inolg .eq. 1 .and. igsh .eq. 0 .and.
     &    cmin .gt. 0.0 .and. cmax .gt. cmin ) then

      write(iot,'(/"#",78("-"))')
      write(iot,'(
     &"z: xorg(1.03)"/
     &"p: xfac[c2*0.05] form[c1/0.05] notc noxt noyt dmin(1)"/
     &"p: ymin(0) ymax(1) xmin(0) xmax(1) cmin(1) cmax(100) nosp"/
     &"hc: y= 0.005 to 0.995 by 0.01 ; x= 0.5 to 0.5 by 1 ;"/
     &" 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20"/
     &"21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40"/
     &"41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60"/
     &"61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80"/
     &"81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100"/
     &"z: xorg(0.0)"/
     &"p: xfac[c2*0.05] form[c1/0.05] ",
     &"nosp afac[c5*0.625] ymin[c3] ymax[c4]"/
     &"p: noxt noxn ytxt(-1) ynum(-1) xltd(-1) xstd(-1) itic(1) ylog"
     &)')

      if( itazl(m) .eq. 0 ) then

         write(iot,'("y: Number ",a15)') hsunit(itunt(m))

      else

         write(iot,'("y: ",200a1)') (itazt(m)(i:i),i=1,itazl(m))

      end if

      if( izlog .eq. 0 ) write(iot,'( "p: ylin")')

      write(iot,'( "#",78("-"))')

      end if

*-----------------------------------------------------------------------

            if( inocm .eq. 1 ) then

               if( nn .eq. 0 ) then

                  write(iot,'(/"wt: s[c5]",/
     &                        a1,"vspace{-3}"/
     &                        "  x &=&",1pe13.4," [cm]"/
     &                        "  Z &=&   all"/
     &                        "e:")')
     &                        yen, xval

               else if( ia .eq. 0 ) then

                  write(iot,'(/"wt: s[c5]",/
     &                        a1,"vspace{-3}"/
     &                        "  x &=&",1pe13.4," [cm]"/
     &                        "  Z &=& ",i3,"   :  ",a3/
     &                        "e:")')
     &                        yen,
     &                        xval, iz, elmnt(iz)

               else

                  write(iot,'(/"wt: s[c5]",/
     &                        a1,"vspace{-3}"/
     &                        "  x     &=&",1pe13.4," [cm]"/
     &                        "  nucl. &=&   ",a8/
     &                        "e:")')
     &                        yen,
     &                        xval, chau

               end if

            end if

*-----------------------------------------------------------------------

      if( inolg .eq. 1 .and. igsh .eq. 0 .and.
     &    cmin .gt. 0.0 .and. cmax .gt. cmin ) then

         write(iot,'(/"#",78("-"))')
         write(iot,'("z: xorg[-1.03/0.05]"/
     &               "p: form[c1] xfac[c2] afac[c5] nosp notf")')

               write(iot,'( "p: xmin(",1p1g14.7,") xmax(",1p1g14.7,
     &                       ")")') zmin, zmax

               write(iot,'( "p: ymin(",1p1g14.7,") ymax(",1p1g14.7,
     &                       ")")') ymin, ymax

      end if

*-----------------------------------------------------------------------

            end do
            end do

*-----------------------------------------------------------------------
*        xz axis (matrix)
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 11 ) then

*-----------------------------------------------------------------------

               inum = 0

                  if( nn .eq. 0 .or. igsh .ne. 0  ) then

                     nc = 1

                  else

                     nc = nn

                  end if

*-----------------------------------------------------------------------

            do jy = 1, ny
            do ic = 1, nc

               yval = ( ym(jy) + ym(jy+1) ) / 2.0d0

               inum = inum + 1

*-----------------------------------------------------------------------

               if( nn .gt. 0 ) then

                  iz = nt(ic) / 1000

                  ia = nt(ic) - iz * 1000

                  if( ia .gt. 0 ) then

                     in = ia - iz

                     call chname(idum,ia,iz,chau)

                  end if

               end if

*-----------------------------------------------------------------------

                  write(iot,'(/"#",78("-"))')

               if( inum .eq. 1 ) then
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

               if( rtfac(m) .lt. 0.d0 )
     &           write(iot,
     &           '( "#restart maximum for normalization:",1p1e14.7)')
     &           facmax(m)

            if( nn .eq. 0 ) then

               write(iot,'("#   no. =",i3,3x,
     &         "iy  =",i3,3x,"Z = all")')
     &                     inum, jy

            else if( ia .eq. 0 ) then

               write(iot,'("#   no. =",i3,3x,
     &         "iy  =",i3,3x,"Z = ",i3," : ",a3)')
     &                     inum, jy, iz, elmnt(iz)

            else

               write(iot,'("#   no. =",i3,3x,
     &         "iy  =",i3,3x,a8)')
     &                     inum, jy, chau

            end if

*-----------------------------------------------------------------------

               if( nn .eq. 0 ) then

                  write(iot,'(/a1,"no. =",i3,
     &                        ",    iy =",i3,
     &                        ",    Z  =  all", a1)')
     &                        cha, inum, jy, cha

               else if( ia .eq. 0 ) then

                  write(iot,'(/a1,"no. =",i3,
     &                        ",    iy =",i3,
     &                        ",    Z  =  ",i3, a1)')
     &                        cha, inum, jy, iz, cha

               else

                  write(iot,'(/a1,"no. =",i3,
     &                        ",    iy =",i3,
     &                        ",    ",a8, a1)')
     &                        cha, inum, jy, chau, cha

               end if

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: z [cm]")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: x [cm]")')

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

*-----------------------------------------------------------------------

         if( igsh .eq. 0 ) then

               do i = 1, nr

                     sek = 0.0

                  if( nn .eq. 0 ) then

                     do kn = 1, mn
                     do kz = 1, mz
                     do il = 0, mm ! frtati 2022/02/18

                        lz  = ikzz( kz, kn )
                        ln  = iknn( kz, kn )

                        sek = sek + tr(i,lz,ln,il,1)

                     end do
                     end do
                     end do

                  else if( ia .eq. 0 ) then

                     do kn = 1, mn
                     do il = 0, mm ! frtati 2022/02/18

                        lz  = ikzz( iz, kn )
                        ln  = iknn( iz, kn )

                        sek = sek + tr(i,lz,ln,il,1)

                     end do
                     end do

                  else

                        lz  = ikzz( iz, in )
                        ln  = iknn( iz, in )

                     do il = 0, mm ! frtati 2022/02/18
                        sek = sek + tr(i,lz,ln,il,1)
                     end do

                  end if

                  val(i) = sek

               end do

                  izlog = 1
                  if( itanl(m) .gt. 0 )
     &            call anset(itang(m),itanl(m),
     &                       form,xfac,afac,cmin,cmax,izlog,inocm,inolg)

            if( izlog .eq. 0 ) then

               do i = 1, nr

                  if( val(i) .gt. dnon .and. cmin .ne. cmax ) then

                     val(i) = ( val(i) - cmin ) / ( cmax - cmin )
                     val(i) = max(0.0d0,val(i))
                     val(i) = min(1.0d0,val(i))

                  else

                     val(i) = -1.0

                  end if

               end do

            else

               do i = 1, nr

                  if( val(i) .gt. dnon .and. cmin .ne. cmax ) then

                     val(i) = log10(val(i)/cmin) / log10(cmax/cmin)
                     val(i) = max(0.0d0,val(i))
                     val(i) = min(1.0d0,val(i))

                  else

                     val(i) = -1.0d0

                  end if

               end do

            end if

         end if

*-----------------------------------------------------------------------

               if( inum .eq. 1 ) then

                  form  = ( xm(nx+1) - xm(1) )
     &                  / ( zm(nz+1) - zm(1) )
                  xfac  = 0.9
                  afac  = 0.8
                  izlog = 1
                  inocm = 1
                  inolg = 1

                  if( itanl(m) .gt. 0 )
     &            call anset(itang(m),itanl(m),
     &                       form,xfac,afac,cmin,cmax,izlog,inocm,inolg)

                  if( form .le. 1.0 ) then

                     scal = form**0.35
                     xfac = xfac / form**0.5
                     xorg = 0.0
                     yorg = min( 20.d0, ( 1.0 / form - 1.0 ) / 2.5 )

                  else

                     scal = 1.0 / form**0.41
                     xfac = xfac / form**0.5
                     xorg = min( 20.d0, ( form - 1.0 ) / 1.5)
                     yorg = 0.0

                  end if

               end if

                  write(iot,'( "set: c1[",f6.3,
     &                            "] c2[",f6.3,
     &                            "] c5[",f6.3,
     &                            "] c6[",f6.3,
     &                            "] c7[",f6.3,
     &                            "] c8[",f6.3,"]")')
     &                    form, xfac, afac, scal, xorg, yorg
                  write(iot,'( "p: h2fs form[c1] xfac[c2]",
     &                 " afac[c5] scal[c6] xorg[c7] yorg[c8] nosp")')

               if( cmin .gt. 0.0 .and. cmax .gt. cmin ) then

                  write(iot,'( "set: c3[",1p1e15.6,"] c4[",
     &                  1p1e15.6,"]")') cmin, cmax

               end if

*-----------------------------------------------------------------------

               zmin = zm(1)
               zmax = zm(nz+1)
               xmin = xm(1)
               xmax = xm(nx+1)

               write(iot,'( "p: xmin(",1p1g14.7,") xmax(",1p1g14.7,
     &                       ")")') zmin, zmax

               write(iot,'( "p: ymin(",1p1g14.7,") ymax(",1p1g14.7,
     &                       ")")') xmin, xmax

               if( itanl(m) .gt. 0 ) then

                  write(iot,'( "p: ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               end if

               if( itsans(m) .gt. 0 ) then

                  call write_sangel(iot,m,0)

               end if

*-----------------------------------------------------------------------
*        rshow
*-----------------------------------------------------------------------

               write(iot,'(/"#",78("-"))')
               write(iot,'( "# rshow")')
               write(iot,'( "#",78("-"))')
               write(iot,'(/"p: legs[c5*0.875]")')

               none = 1
               iaxs = 3
               iuni = itrsh(m) * 2 - 1 + igsh
               ires = itres(m)
               widt = rtwid(m)

            if( igsh .eq. 0 ) then

               call gshow(2,iot,iaxs,iuni,widt,ires,
     &                    inum,igsh,igser,gfnam,igfmn,
     &                    nz+1,nx+1,none,zm,xm,yval,ixyz(1),
     &                    nr,1,krr,val,itmtr(m,4),itglt(m))

            else

               call gshow(0,iot,iaxs,iuni,widt,ires,
     &                    inum,igsh,igser,gfnam,igfmn,
     &                    nz+1,nx+1,none,zm,xm,yval,ixyz(1),
     &                    1,1,krr,vll,itmtr(m,4),itglt(m))

            end if

*-----------------------------------------------------------------------

      if( inolg .eq. 1 .and. igsh .eq. 0 .and.
     &    cmin .gt. 0.0 .and. cmax .gt. cmin ) then

      write(iot,'(/"#",78("-"))')
      write(iot,'(
     &"z: xorg(1.03)"/
     &"p: xfac[c2*0.05] form[c1/0.05] notc noxt noyt dmin(1)"/
     &"p: ymin(0) ymax(1) xmin(0) xmax(1) cmin(1) cmax(100) nosp"/
     &"hc: y= 0.005 to 0.995 by 0.01 ; x= 0.5 to 0.5 by 1 ;"/
     &" 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20"/
     &"21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40"/
     &"41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60"/
     &"61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80"/
     &"81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100"/
     &"z: xorg(0.0)"/
     &"p: xfac[c2*0.05] form[c1/0.05] ",
     &"nosp afac[c5*0.625] ymin[c3] ymax[c4]"/
     &"p: noxt noxn ytxt(-1) ynum(-1) xltd(-1) xstd(-1) itic(1) ylog"
     &)')

      if( itazl(m) .eq. 0 ) then

         write(iot,'("y: Number ",a15)') hsunit(itunt(m))

      else

         write(iot,'("y: ",200a1)') (itazt(m)(i:i),i=1,itazl(m))

      end if

      if( izlog .eq. 0 ) write(iot,'( "p: ylin")')

      write(iot,'( "#",78("-"))')

      end if

*-----------------------------------------------------------------------

            if( inocm .eq. 1 ) then

               if( nn .eq. 0 ) then

                  write(iot,'(/"wt: s[c5]",/
     &                        a1,"vspace{-3}"/
     &                        "  y &=&",1pe13.4," [cm]"/
     &                        "  Z &=&   all"/
     &                        "e:")')
     &                        yen, yval

               else if( ia .eq. 0 ) then

                  write(iot,'(/"wt: s[c5]",/
     &                        a1,"vspace{-3}"/
     &                        "  y &=&",1pe13.4," [cm]"/
     &                        "  Z &=& ",i3,"   :  ",a3/
     &                        "e:")')
     &                        yen,
     &                        yval, iz, elmnt(iz)

               else

                  write(iot,'(/"wt: s[c5]",/
     &                        a1,"vspace{-3}"/
     &                        "  y     &=&",1pe13.4," [cm]"/
     &                        "  nucl. &=&   ",a8/
     &                        "e:")')
     &                        yen,
     &                        yval, chau

               end if

            end if

*-----------------------------------------------------------------------

      if( inolg .eq. 1 .and. igsh .eq. 0 .and.
     &    cmin .gt. 0.0 .and. cmax .gt. cmin ) then

         write(iot,'(/"#",78("-"))')
         write(iot,'("z: xorg[-1.03/0.05]"/
     &               "p: form[c1] xfac[c2] afac[c5] nosp notf")')

               write(iot,'( "p: xmin(",1p1g14.7,") xmax(",1p1g14.7,
     &                       ")")') zmin, zmax

               write(iot,'( "p: ymin(",1p1g14.7,") ymax(",1p1g14.7,
     &                       ")")') xmin, xmax

      end if

*-----------------------------------------------------------------------

            end do
            end do

*-----------------------------------------------------------------------

         end if

            call prestart(m,iot) !OBINATA(2012.8.20)

            close(iot)

         if( iteps(m) .ne. 0 .and. itaxs(m,iax) .ne. 13 ) then

            open(iot, file = fname, status = 'unknown' )
            call a_angel(idasa,fname)

         end if

  900 continue

*-----------------------------------------------------------------------
      deallocate( ixyz )
      deallocate( lr,vl,val ) !FURUTA20190204
      if(itfoam(m).eq.2) deallocate(xcm) !FURUTA20191028

      return
      end


************************************************************************
*                                                                      *
      subroutine tyilrz(ncol,m,mz,mn,mm,nl,lt, ! frtati 2022/02/18 added mm
     &                  nr,nz,nm,rm,zm,mt,ikzz,iknn,tr,
     &                  trEVENT,
     &                  itrmax,itrmin)
*                                                                      *
*       nuclear yield (or production) tally in r-z mesh                *
*       last modified by K.Niita on 2014/11/11                         *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*      ncol  ..... reaction type                                       *
*             11 : termination by energy cut-off                       *
*             13 : nuclear reaction (n,x)                              *
*             14 : nuclear reaction (n,n'x)                            *
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
*        kcoll : =  0, normal                                          *
*                =  1, high energy fission                             *
*                =  2, high energy absorption                          *
*                =  3, low energy n elastic                            *
*                =  4, low energy n non-elastic                        *
*                =  5, low energy n fission                            *
*                =  6, low energy n absorption                         *
*                                                                      *
************************************************************************
      use QMD_COOD2_MOD, only : nqmdm
      use MMBANKMOD !FURUTA
      use partmod, only: itmxpt,itpan,itpat,jtpat,iznmmx,iznmturn ! frtati 2022/05/02
      use moddas_region ! S.H. (2022.11.9)
*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'param00.inc'
      include 'param01.inc'

*-----------------------------------------------------------------------

      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
!$OMP THREADPRIVATE(/icomon/)
      common /wtsave/ oldwt
!$OMP THREADPRIVATE(/wtsave/)
      common /qparm/  ielas,icasc,iqstep,lvlopt,igamma
      common /bparm/  andt,jevap,npidk
      common /geosig/ geosig(250)
      common /cparm/  maxbch,maxcas
      common /mathzn/ mathz, mathn, jcoll, kcoll
!$OMP THREADPRIVATE(/mathzn/)

      common /trstar/ itrstar
!$OMP THREADPRIVATE(/trstar/)

      common /kmat1d/ idmn(0:kvlmax), idnm(kvmmax)
      common /clionprd/  lionprd
      common /eparm/  esmax, esmin, emin(20)

      common /regdc/  idrg(kvlmax), idgr(kvmmax)
      common /tlgeom/ iblz1,iblz2
!$OMP THREADPRIVATE(/tlgeom/)
*-----------------------------------------------------------------------

      common /clustf/ nclst, iclust(nnn)
!$OMP THREADPRIVATE(/clustf/)
      common /clustt/ nclsts, iclusts(nnn)
!$OMP THREADPRIVATE(/clustt/)
      common /clustw/ jclusts(0:8,nnn),  qclusts(0:12,nnn)
!$OMP THREADPRIVATE(/clustw/)
      common /clustl/ rumpal(0:20), numpal(0:20)
!$OMP THREADPRIVATE(/clustl/)
      common /cntcls/ jcount(3,nnn)
!$OMP THREADPRIVATE(/cntcls/)

*-----------------------------------------------------------------------

      dimension numsav(0:20), rumsav(0:20)

*-----------------------------------------------------------------------

      common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                itout(itlmax), ittwo(itlmax)
      common /tall08/ rtrx0(itlmax), rtry0(itlmax)
      common /tall12/ itman(itlmax), itmat(itlmax), itmct(itlmax),
     &                itnun(itlmax), itnuc(itlmax),
     &                itndz(itlmax), itndn(itlmax),
     &                itnkz(itlmax), itnkn(itlmax)
      common /tall22/ itmcn(itlmax), itmtn(itlmax), itmtt(itlmax)
      common /tall30/ itnda(itlmax)
      common /tall36/ itdpo(itlmax)
      common /tall37/ itcnt(9,itlmax)
      common /tall46/ itmtr(itlmax,4), rtmtr(itlmax,13)
      common /tall52/ itprd(itlmax)

      common /tall82/ itcnth(9,itlmax)
      common /tall83/ itnzn(itlmax), itndm(itlmax)

*-----------------------------------------------------------------------

      dimension   lt(nl)
      dimension   rm(nr+1)
      dimension   zm(nz+1)
      dimension   mt(nm)
      dimension   ikzz(maxpt,maxnt), iknn(maxpt,maxnt)
      dimension   tr(nr,nz,mz,mn,0:mm,2)
      dimension   trEVENT(nr,nz,mz,mn,0:mm)        !OBINATA(2012.8.20): as Ct
      real(8),allocatable,save:: tr0(:,:,:,:,:) !OBINATA(2012.8.20): as C
      dimension   itrmax(4),itrmin(4)

      dimension mnz(mxprodxs),mna(mxprodxs),xxn(mxprodxs),ildd(mxprodxs)

      dimension   ips(itmxpt) ! frtati 2021/10/05 6 -> itmxpt
      dimension   ncntt(3)


*-----------------------------------------------------------------------
      common /cntreg/ ncntc(3), ncreg(3), incrc(3), incrt(3), kcont(3)
      common /cntpat/ icpan(3), icpat(3,20,2)
      dimension     idas(1)
      equivalence ( das, idas )
      common /pnsave/ egs, uus, vvs, wws, wts, tms, nms, nct(3)
!$OMP THREADPRIVATE(/pnsave/)

*-----------------------------------------------------------------------
      common /cntmx/ ncntmx(3), nctmxsr(3)
!$OMP THREADPRIVATE(/cntmx/)

      common /mpi00/ npe, me
      real*8,allocatable :: tryld(:,:)

*-----------------------------------------------------------------------

      common /jcomon/ nabov,nobch,nocas,nomax
!$OMP THREADPRIVATE(/jcomon/)
      common /stat / istdev, irestart, ireschk

      integer, allocatable, save :: iclusav(:), jclusav(:,:)
      double precision, allocatable, save :: qclusav(:,:)
!$OMP THREADPRIVATE(iclusav, jclusav, qclusav)

       if( .not. allocated(iclusav) ) then ! initial allocation
        allocate(iclusav(nqmdm),jclusav(0:8,nqmdm),qclusav(0:12,nqmdm))
        iclusav = 0
        jclusav = 0
        qclusav = 0.d0
       elseif( ubound(iclusav,1) .lt. nqmdm ) then ! extend array
        deallocate(iclusav,jclusav,qclusav)
        allocate(iclusav(nqmdm),jclusav(0:8,nqmdm),qclusav(0:12,nqmdm))
        iclusav = 0
        jclusav = 0
        qclusav = 0.d0
       endif

*-----------------------------------------------------------------------

      if (istdev .eq. 2) then
        call readitrminmax4(itrmin,itrmax,(/nr,nz,mz,mn/),
     &                      mir,miz,mimz,mimn,
     &                      mxr,mxz,mxmz,mxmn)
      endif

*-----------------------------------------------------------------------
* check of history counter
*-----------------------------------------------------------------------
      ihistcount = 0 ! history counter index, T.Sato 2022/12/20 for speed up
         if ( ncol.eq.0 .or. ncol.eq.4 ) then
           do i = 1, 3
             if( itcnth(i,m) .eq. 1 ) then
               if( ncntmx(i) .lt. itcnth(i*2+2,m) .or.
     &             ncntmx(i) .gt. itcnth(i*2+3,m) ) then
                 ihistcount = 1
               end if
             end if
           end do
         end if

*-----------------------------------------------------------------------
*        source particle or end of batch ( in case of istdev = 2 )
*
* OBINATA(2012.8.20): change tr(,,,,3) to trEVENT(,,,,)
*-----------------------------------------------------------------------

         if (( ncol .eq. 0 .or. ncol .eq. 4 )
     &                              .and. istdev .eq. 2) then
           if ((nocas.gt.1.or.ncol.eq.0) .and. ihistcount.ne.1 ) then

            if( itnzn(m).eq.0 ) then ! frtati 2022/02/18
             do imn = mimn,mxmn
             do imz = mimz,mxmz
             do iz = miz,mxz
             do ir = mir,mxr
               tr(ir,iz,imz,imn,:,1) = tr(ir,iz,imz,imn,:,1)
     &                               + trEVENT(ir,iz,imz,imn,:)
               tr(ir,iz,imz,imn,:,2) = tr(ir,iz,imz,imn,:,2)
     &                               + trEVENT(ir,iz,imz,imn,:)**2
             enddo
             enddo
             enddo
             enddo
            else
             do imz = 1,iznmmx(m)
             do iz = miz,mxz
             do ir = mir,mxr
               tr(ir,iz,imz,1,0,1) = tr(ir,iz,imz,1,0,1)
     &                               + trEVENT(ir,iz,imz,1,0)
               tr(ir,iz,imz,1,0,2) = tr(ir,iz,imz,1,0,2)
     &                               + trEVENT(ir,iz,imz,1,0)**2
             enddo
             enddo
             enddo
            end if

           end if

          if( itnzn(m).eq.0 ) then ! frtati 2022/02/18
           do imn = mimn,mxmn
           do imz = mimz,mxmz
           do iz = miz,mxz
           do ir = mir,mxr
             trEVENT(ir,iz,imz,imn,:) = 0
           enddo
           enddo
           enddo
           enddo
          else
           do imz = 1,iznmmx(m)
           do iz = miz,mxz
           do ir = mir,mxr
             trEVENT(ir,iz,imz,1,0) = 0
           enddo
           enddo
           enddo
          end if

           call resetitrminmax(itrmin,itrmax,4,(/nr,nz,mz,mn/))

         end if

*-----------------------------------------------------------------------
*        end of batch ( in case of istdev = 1 )
*
* OBINATA(2012.8.20): modificate for thread parallel
*-----------------------------------------------------------------------

         if ( ncol .eq. 0 .and. istdev .eq. 1) then
!$OMP MASTER
             allocate( tr0(nr,nz,mz,mn,0:mm) ) ! frtati 2022/02/18 2 -> mm
             tr0(:,:,:,:,:) = 0.d0
!$OMP END MASTER
!$OMP BARRIER
!$OMP CRITICAL (tyilrz_crit_ist1)
             tr0(:,:,:,:,:) = tr0(:,:,:,:,:) + trEVENT(:,:,:,:,:)
!$OMP END CRITICAL (tyilrz_crit_ist1)
!$OMP BARRIER
!$OMP MASTER
             tr(:,:,:,:,:,1) = tr(:,:,:,:,:,1) + tr0(:,:,:,:,:) / maxcas
             tr(:,:,:,:,:,2) = tr(:,:,:,:,:,2)
     &                   + ( tr0(:,:,:,:,:) / maxcas ) ** 2
             deallocate( tr0 )
!$OMP END MASTER

           trEVENT(:,:,:,:,:) = 0

         end if

*-----------------------------------------------------------------------
*        rearrangement of iznm for MPI run ! frtati 2022/05/02
*-----------------------------------------------------------------------

         if( npe.gt.1 .and. itnzn(m).ne.0 .and. ncol.eq.0 ) then
           call paraiznm(m,1)
           if( me.gt.1 ) then
             allocate( tryld(iznmmx(m),2) )
             do iz = 1, nz
             do ir = 1, nr
               do imz = 1, iznmmx(m)
                 tryld(imz,:) = tr(ir,iz,imz,1,0,:)
               end do
               do imz = 1, iznmmx(m)
                 tr(ir,iz,imz,1,0,:) = tryld(iznmturn(imz,m),:)
               end do
             end do
             end do
             deallocate( tryld )
           end if
         end if

*-----------------------------------------------------------------------
*        check of ncol and nclsts ( outgoing particles )
*-----------------------------------------------------------------------

         if( ncol .eq. 11 ) then

            if( itprd(m) .ne. 1 ) return

               iccol = 0
               npart = 1

               if( ityp .lt. 15 .or. ityp .gt. 19 ) return

               jz    = jtyp
               jn    = ktyp - ktyp / 1000000 * 1000000 - jz

               ipart = ityp
               tlw   = oldwt
               il = 0  ! T.Sato 2016/3/29, isomer is not considered for output=cutoff

               ncntt(1) = ncnt(ibknct+1,no,ipomp+1)
               ncntt(2) = ncnt(ibknct+2,no,ipomp+1)
               ncntt(3) = ncnt(ibknct+3,no,ipomp+1)

         else if( ncol .eq. 13 .or. ncol .eq. 14 ) then

            if ( itnda(m) .ge. 2 .and. jcoll .eq. 9 ) then
               if( nclsts .lt. 0 ) return
            else
               if( nclsts .le. 0 ) return
            end if

            if( ( itnda(m) .eq. 2 .and.
     &          ( jcoll .eq. 6 .or. jcoll .eq. 9 ) ) .or.
     &          ( itnda(m) .eq. 3 .and.
     &          ( jcoll .eq. 6 .or. jcoll .eq. 9 .or.
     &            jcoll .eq. 4 .or. jcoll .eq. 5 .or.
     &            jcoll .eq. 15 .or. jcoll .eq. 10 ) ) ) then

            else
               if( itdpo(m) .eq. 0 .and.
     &           ( jcoll .eq. 3 .or. kcoll .eq. 3 ) ) return
               if( itdpo(m) .eq. -1 .and. itrstar .eq. 1 ) return
            end if

               iccol = 1

               npart = nclsts

               ncntt(1) = jcount(1,1)
               ncntt(2) = jcount(2,1)
               ncntt(3) = jcount(3,1)

         else

            return

         end if

*-----------------------------------------------------------------------
*        check of counter
*-----------------------------------------------------------------------

            do i = 1, 3

               if( itcnt(i,m) .ne. 0 ) then

                  if( ncntt(i) .lt. itcnt(i*2+2,m) .or.
     &                ncntt(i) .gt. itcnt(i*2+3,m) ) return

               end if

            end do

*-----------------------------------------------------------------------
*        check of mat
*-----------------------------------------------------------------------

            if( nl .gt. 0 ) then

                  do i = 1, nl

                     if( itmcn(m) .gt. 0 .and.
     &                   idmn(mat) .eq. lt(i) ) goto 502
                     if( itmcn(m) .lt. 0 .and.
     &                   idmn(mat) .eq. lt(i) ) return

                  end do

                     if( itmcn(m) .gt. 0 ) return

            end if

  502          continue

*-----------------------------------------------------------------------

         if( iccol .eq. 1 ) then

*-----------------------------------------------------------------------
*           check of particles
*-----------------------------------------------------------------------

            call pcheck(m,itpan(m),ityp,ktyp,jtyp,ipn,ips)

               if( ipn .eq. 0 ) return

*-----------------------------------------------------------------------
*           mother and energy
*-----------------------------------------------------------------------

               ata   = dble( mathz + mathn )
               atz   = dble( mathz )
               mmas  = nint( ata )
               mchg  = nint( atz )
               eein  = ec(ibkec+no,ipomp+1)

*-----------------------------------------------------------------------
*           for specific mother nuclei
*-----------------------------------------------------------------------

            if( nm .gt. 0 ) then

               do i = 1, nm

                  iz = mt(i) / 1000
                  ia = mt(i) - iz * 1000

                  if( ( ia .gt. 0 .and.
     &                  mmas .eq. ia .and. mchg .eq. iz ) .or.
     &                ( ia .eq. 0 .and. mchg .eq. iz ) ) then

                     if( itmct(m) .gt. 0 ) goto 30
                     if( itmct(m) .lt. 0 ) return

                  end if

               end do

                  if( itmct(m) .gt. 0 ) return

            end if

   30       continue

         end if

*-----------------------------------------------------------------------
*        transform positions
cKN      yield position is assumed to be the reaction point.
*-----------------------------------------------------------------------

               xpart = xc(ibkxc+no,ipomp+1)
               ypart = yc(ibkyc+no,ipomp+1)
               zpart = zc(ibkzc+no,ipomp+1)

               call trnsxx(xpart,ypart,zpart,
     &                     xxc,yyc,zzc,itmtr(m,4))

*-----------------------------------------------------------------------
*              check z mesh and r mesh region
*-----------------------------------------------------------------------

                  xpp = xxc
                  ypp = yyc
                  zpp = zzc

                  x0 = rtrx0(m)
                  y0 = rtry0(m)

                  dis1 = sqrt( ( xpp - x0 )**2
     &                       + ( ypp - y0 )**2 )

               if(  zpp .lt. zm(1) .or. zpp  .ge. zm(nz+1) ) return
               if( dis1 .lt. rm(1) .or. dis1 .ge. rm(nr+1) ) return

                  do i = 2, nz + 1

                     if( zpp .lt. zm(i) ) goto 38

                  end do

   38          izc = i - 1

                  do i = 2, nr + 1

                     if( dis1 .lt. rm(i) ) goto 39

                  end do

   39          irc = i - 1

*-----------------------------------------------------------------------
*        normal weight
*-----------------------------------------------------------------------

               idoy = 1
               idev = 1
               wyld = oldwt

*-----------------------------------------------------------------------
*        for special : repeated nuclear reactions
*-----------------------------------------------------------------------

         if( iccol .eq. 1 .and. itspc(m) .gt. 0 .and.
     &       mathz .gt. 2 .and.
     &     ( jcoll .eq. 4 .or. jcoll .eq. 5 ) ) then

               idoy = itspc(m)
               idev = 10
               wyld = oldwt / dble( idoy * idev )

*-----------------------------------------------------------------------
*           save normal output and reaction parameters
*-----------------------------------------------------------------------

                  nclssav = nclsts

            do i = 1, nclsts

                  iclusav(i) = iclusts(i)

               do k = 0, 8

                  jclusav(k,i) = jclusts(k,i)

               end do

               do k = 0, 12

                  qclusav(k,i) = qclusts(k,i)

               end do

            end do

               do k = 0, 20

                  numsav(k) = numpal(k)
                  rumsav(k) = rumpal(k)

               end do

                  iprj  = ityp
                  kprj  = ktyp

                  if( ityp .lt. 15 ) then

                     bmax = sqrt( geosig(mmas) * 100.0 / 3.1415926 )

                  else if( ityp .ge. 15 ) then

                     ap = dble( ktyp - ktyp / 1000000 * 1000000 )
                     zp = dble( ktyp / 1000000 )

                     call sighi(ap,zp,eein,ata,zta,signe,sigel,bmax)

                  end if

         end if

*-----------------------------------------------------------------------
*     repeated do loop
*-----------------------------------------------------------------------

         do ireac = 1, idoy

*-----------------------------------------------------------------------
*     repeat calculation of nuclear reactions
*-----------------------------------------------------------------------

            if( idoy .gt. 1 .and. ireac .gt. 1 ) then

                  ipim = 0

   22          continue

                  call ncasc(1,iprj,kprj,eein,mmas,mchg,bmax)

               if( nclst .lt. 0 ) then

                     ipim = ipim + 1

                     if( ipim .le. 20 ) goto 22

               end if

            end if

*-----------------------------------------------------------------------
*     repeat calculation of evaporation
*-----------------------------------------------------------------------

            do ievap = 1, idev

               if( idev .gt. 1 ) then

                  call nevap(1)

                  do i = 1, nclsts

                     qclusts(8,i) = wyld * qclusts(8,i)

                  end do

                  npart = nclsts

               end if

*-----------------------------------------------------------------------
*           booking after nuclear reactions
*-----------------------------------------------------------------------

                  istat = 1
                  mnx   = 0

*-----------------------------------------------------------------------
*           NDATA=2,3 : Replace with yield data (activation cross section
*-----------------------------------------------------------------------

               if( iccol .eq. 1 .and.
     &           ( itnda(m) .eq. 2 .and.
     &           ( jcoll .eq. 6 .or. jcoll .eq. 9 ) ) .or.
     &           ( itnda(m) .eq. 3 .and.
     &           ( jcoll .eq. 6 .or. jcoll .eq. 9 .or.
     &             jcoll .eq. 4 .or. jcoll .eq. 5 .or.
     &             jcoll .eq. 15 .or. jcoll .eq. 10 ) ) ) then

                  mkk = mat
                  icl = idgr(iblz1)
                  jcl = jcoll

                  call prodxs2(ktyp,mchg,mmas,eein,
     &                 istat,mnx,mnz,mna,xxn,xxs,
     &                 ildd,mkk,icl,jcl,m)

                  iexclight = 0 ! Exclude light ions
                  if( itprd(m) .eq. 1 ) then ! output = cutoff
                     iexclight = 1
                     if ( jcoll.eq.6 .and.
     &                    ktyp.eq.2112 .and. eein.lt.20.0 )
     &                    iexclight = 0
                  end if

                  if( istat .eq. 0 .and. mnx .gt. 0 ) then

                   do nnx = 1, mnx

                      jz = mnz(nnx)
                      jn = mna(nnx) - mnz(nnx)
                      il = ildd(nnx)

                    if ( iexclight.ne.1 .or. jz.gt.2 ) then

                     if( jz .gt. 0 .and. jz .le. maxpt .and.
     &                      jn .gt. 0 .and. jn .le. maxnt ) then

                        kz = ikzz(jz,jn)
                        kn = iknn(jz,jn)

                      if( itnzn(m).eq.0 ) then ! frtati 2022/02/18
                       if( kz .le. mz .and. kn .le. mn ) then
                          trEVENT(irc,izc,kz,kn,il) =
     &                       trEVENT(irc,izc,kz,kn,il) + wyld * xxn(nnx)

                        if (istdev .eq. 2) then
                           call setitrminmax(itrmin,itrmax,1,4,
     &                          (/irc,izc,kz,kn/))
                        endif

                       end if

                      else
                         trEVENT(irc,izc,igetiznm(kz,kn,il,m),1,0) =
     &                       trEVENT(irc,izc,igetiznm(kz,kn,il,m),1,0) +
     &                        wyld * xxn(nnx)
                       if (istdev .eq. 2) then
                          call setitrminmax(itrmin,itrmax,1,4,
     &                         (/irc,izc,kz,kn/))
                       endif
                      end if

                     end if

                    end if

                   end do

                  end if

               end if

*-----------------------------------------------------------------------
*           NDATA1 : nuclear data for 4He, 14N, 16O
*-----------------------------------------------------------------------

               if( iccol .eq. 1 .and.
     &             itnda(m) .eq. 1 .and. jcoll .eq. 4 ) then

                     call prodxs(ktyp,mchg,mmas,eein,
     &                           istat,mnx,mnz,mna,xxn,xxs,ildd)

                  if( istat .eq. 0 .and. mnx .gt. 0 ) then

                     do nnx = 1, mnx

                        jz = mnz(nnx)
                        jn = mna(nnx) - mnz(nnx)
                        il = ildd(nnx)

                        if( jz .gt. 0 .and. jz .le. maxpt .and.
     &                      jn .gt. 0 .and. jn .le. maxnt ) then

                           kz = ikzz(jz,jn)
                           kn = iknn(jz,jn)

                        if( itnzn(m).eq.0 ) then ! frtati 2022/02/18

                        if( kz .le. mz .and. kn .le. mn ) then

                           trEVENT(irc,izc,kz,kn,il) =
     &                     trEVENT(irc,izc,kz,kn,il) + wyld * xxn(nnx)

                           if (istdev .eq. 2) then
                             call setitrminmax(itrmin,itrmax,1,4,
     &                                         (/irc,izc,kz,kn/))
                           endif

                        end if

                        else
                           trEVENT(irc,izc,igetiznm(kz,kn,il,m),1,0) =
     &                     trEVENT(irc,izc,igetiznm(kz,kn,il,m),1,0) +
     &                     wyld * xxn(nnx)
                           if (istdev .eq. 2) then
                             call setitrminmax(itrmin,itrmax,1,4,
     &                                         (/irc,izc,kz,kn/))
                           endif
                        end if

                        end if

                     end do

                  end if

               end if



*-----------------------------------------------------------------------
*           normal case
*-----------------------------------------------------------------------

            if( npart .gt. 0 ) then

            do 500 j = 1, npart

*-----------------------------------------------------------------------

                  if( iccol .eq. 1 ) then

                     ipart = jclusts(3,j)
                     tlw   = qclusts(8,j)
                     jz    = jclusts(1,j)
                     jn    = jclusts(2,j)
                     il    = jclusts(8,j)

                  end if

*-----------------------------------------------------------------------
* When ndata=2 or 3, counter values change by ndata in [counter]
                 if ( (ipart.ge.3 .and. ipart.le.13) .or.
     &                 ipart.ge.19 ) then

                  if( ( itnda(m) .eq. 2 .or. itnda(m) .eq. 3 ) .and.
     &                 ( istat .eq. 0 .or. istat .eq. 3 ) ) then

                   if( ncntc(1) .eq. 1 .or. ncntc(2) .eq. 1 .or.
     &                    ncntc(3) .eq. 1 ) then
                    do k = 1, 3
                     ndatcount = 0
                     if( ncntc(k) .eq. 1 ) then
                      knn  = ityp
                      kcg  = jtyp
                      kkf  = ktyp
                      call pcchck(k,knn,kkf,kcg,icpan,icpat,icc)

                      if( icc .eq. 1 ) then
                       kdsm = kcont(k)
                       ldsm = 22
                       ndatcount = idas_kcont(kdsm+ldsm)
                      end if

                      if ( ndatcount .ne. 0 ) then

                       if ( ncol .eq. 13 ) then
                        ncnta(ibknct+k,j,ipomp+1)
     &                         = ncnta(ibknct+k,j,ipomp+1)
     &                         + ndatcount

                       else if ( ncol .eq. 14 ) then
                        if ( j .eq. 1 ) then
                         nct(k) = nct(k) + ndatcount
                        else if (j .gt. 1 ) then
                         ncnta(ibknct+k,j-1,ipomp+1)
     &                          = ncnta(ibknct+k,j-1,ipomp+1)
     &                          + ndatcount
                        end if

                       end if

                      end if

                     end if
                    end do
                   end if

                  end if

                 end if

*-----------------------------------------------------------------------

                  if( ipart .lt. 15 .or. ipart .gt. 19 ) goto 500

                  if( iccol .eq. 1 .and. itprd(m) .eq. 1 ) then

                     emint = emin(ipart) * dble( jz + jn )

                     if( qclusts(7,j) .gt. emint ) goto 500

                  end if

*-----------------------------------------------------------------------

                     if( ( itnda(m) .eq. 2 .or. itnda(m) .eq. 3 ) .and.
     &                 istat .eq. 0 .and.
     &                 ( iexclight.ne.1 .or. jz.gt.2 ) ) goto 500

                     if( ( itnda(m) .eq. 2 .or. itnda(m) .eq. 3 ) .and.
     &                    istat .eq. 3 .and. lionprd .ne. 1) goto 500

                     if( itnda(m) .eq. 1 .and.
     &                   istat .eq. 0 .and. mnx .gt. 0 ) then

                        do nnx = 1, mnx

                           if( jz .eq. mnz(nnx) .and.
     &                         jn .eq. mna(nnx) - mnz(nnx) ) goto 500

                        end do

                     end if

                        if( jz .gt. 0 .and. jz .le. maxpt .and.
     &                      jn .gt. 0 .and. jn .le. maxnt ) then

                           kz = ikzz(jz,jn)
                           kn = iknn(jz,jn)

                        if( itnzn(m).eq.0 ) then ! frtati 2022/02/18

                        if( kz .le. mz .and. kn .le. mn ) then

                           trEVENT(irc,izc,kz,kn,il) =
     &                     trEVENT(irc,izc,kz,kn,il) + tlw

                           if (istdev .eq. 2) then
                             call setitrminmax(itrmin,itrmax,1,4,
     &                                         (/irc,izc,kz,kn/))
                           endif

                        end if

                        else
                           trEVENT(irc,izc,igetiznm(kz,kn,il,m),1,0) =
     &                     trEVENT(irc,izc,igetiznm(kz,kn,il,m),1,0) +
     &                     tlw
                           if (istdev .eq. 2) then
                             call setitrminmax(itrmin,itrmax,1,4,
     &                                         (/irc,izc,kz,kn/))
                           endif
                        end if

                        end if

  500       continue

            end if

*-----------------------------------------------------------------------

            end do

         end do

*-----------------------------------------------------------------------
*        store the decrease of mother by negative weight
*-----------------------------------------------------------------------

        if( ncol .ne. 11 ) then ! S.H. 2021.12.22
         if( itdpo(m) .eq. -1 .and. itrstar .eq. 0 ) then

            if( mathz .gt. 0 .and. mathz .le. maxpt .and.
     &          mathn .gt. 0 .and. mathn .le. maxnt ) then

               kz = ikzz(mathz,mathn)
               kn = iknn(mathz,mathn)
             if( itnzn(m).eq.0 ) then ! frtati 2022/02/18
              if( kz .le. mz .and. kn .le. mn ) then ! T.Sato 2024/05/02
               trEVENT(irc,izc,kz,kn,0) =
     &                    trEVENT(irc,izc,kz,kn,0) - oldwt

               if (istdev .eq. 2) then
                 call setitrminmax(itrmin,itrmax,1,4,
     &                                (/irc,izc,kz,kn/))
               endif
              endif
             else
               trEVENT(irc,izc,igetiznm(kz,kn,0,m),1,0) =
     &         trEVENT(irc,izc,igetiznm(kz,kn,0,m),1,0) - oldwt
                  if (istdev .eq. 2) then
                    call setitrminmax(itrmin,itrmax,1,4,
     &                                (/irc,izc,kz,kn/))
                  endif
             end if

            end if

         end if
        end if

*-----------------------------------------------------------------------
*        restore normal output
*-----------------------------------------------------------------------

         if( idoy .gt. 1 ) then

                  nclsts = nclssav

            do i = 1, nclsts

                  iclusts(i) = iclusav(i)
                  iclusav(i) = 0

               do k = 0, 8

                  jclusts(k,i) = jclusav(k,i)
                  jclusav(k,i) = 0

               end do

               do k = 0, 12

                  qclusts(k,i) = qclusav(k,i)
                  qclusav(k,i) = 0.d0

               end do

            end do

               do k = 0, 20

                  numpal(k) = numsav(k)
                  rumpal(k) = rumsav(k)

               end do

         end if

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine pyildrz(m,mz,mn,mm,nfr,nfz, ! frtati 2022/02/18 added mm
     &                   nr,nz,nn,rm,zm,nt,ikzz,iknn,
     &                   tr,idasa)
*                                                                      *
*       output of nuclear yield (or production) tally in r-z mesh      *
*       last modified by K.Niita on 2004/04/27                         *
*                                                                      *
************************************************************************
      use sangelmod, only: itsans

*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'param01.inc'

*-----------------------------------------------------------------------

      common /taliin/ rsouin, nzztin, nrgnin

*-----------------------------------------------------------------------

      common /verjam/ versn, lastr, iyeav, imonv, idayv

      common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                itout(itlmax), ittwo(itlmax)
      common /tall03/ itrty(itlmax), itrnm(itlmax), itrrg(itlmax),
     &                rtrmi(itlmax), rtrma(itlmax), rtrdl(itlmax)
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

      common /tall21/ rtfac(itlmax)
      common /tall35/ iteps(itlmax)
      common /tall46/ itmtr(itlmax,4), rtmtr(itlmax,13)
      common /tall83/ itnzn(itlmax), itndm(itlmax)

      common /istcut/ ist_cut, ist_bat

      common /fact01/ facmax(itlmax) ! kitamura23/03/31

*-----------------------------------------------------------------------

      common /tcntl/  icntl, inucr

*-----------------------------------------------------------------------

      common /cparm/  maxbch,maxcas
      common /jcomon/ nabov,nobch,nocas,nomax
!$OMP THREADPRIVATE(/jcomon/)

      common /talout/ itall

      character fname*100, fnume*3

*-----------------------------------------------------------------------
      character gfnam*100
      common /tall60/ itger(itlmax)

      common /tall67/ itstd(itlmax), rtstd(itlmax)

      common /mpi00/ npe, me

*-----------------------------------------------------------------------

      dimension   rm(nr+1)
      dimension   zm(nz+1)
      dimension   nt(nn)
      dimension   ikzz(maxpt,maxnt), iknn(maxpt,maxnt)
      dimension   tm(maxnt+maxpt,2)
      dimension   fm(nfz,nfr)
      dimension   tr(nr,nz,mz,mn,0:mm,2) ! frtati 2022/02/18 2 -> mm

*-----------------------------------------------------------------------

      character hsunit(2)*15

      data hsunit / '[1/source]     ',
     &              '[1/cm^3/source]'/

*-----------------------------------------------------------------------

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

      data pi/3.14159265d+0/

      character dc2*4

CCSE add to solve undefined parameter (2018.07.31) >>>>>
      data ipstep / 12 /
CCSE add to solve undefined parameter (2018.07.31) <<<<<
CCSE add for mesh=r-z parameter (2017.11.30) >>>>>
      character erfnm*100
CCSE add for mesh=r-z parameter (2017.11.30) <<<<<

      character chau*8
      character cha*1
      data cha /"'"/

      character rpa*1
      data rpa /'}'/
      character yen*1

*-----------------------------------------------------------------------
*        set mesh volume
*-----------------------------------------------------------------------

                 vl(jr,jz) = pi * ( rm(jr+1)**2 - rm(jr)**2 )
     &                          * ( zm(jz+1) - zm(jz) )

*-----------------------------------------------------------------------

      yen  = char(92)
      igsh = 0

*-----------------------------------------------------------------------
*        c1 : nomalization for source
*-----------------------------------------------------------------------

            if ( rsouin .gt. 0d0 ) then ! S.H. (2016.7.26)

               c1 = 1.0d+0 / rsouin

            else

               c1 = 0d0

            end if

*-----------------------------------------------------------------------

                  cmax = 0.0
                  cmin = 1.e+33
                  isdz = 0
                  stdm = 0.0
                  itstd(m) = 1

            facmax(m) = 1.d0
            if( rtfac(m) .lt. 0.d0 ) then
               facmax(m) = 0.d0

               do 101 jr = 1, nr
               do 101 jz = 1, nz
               do 101 iz = 1, mz
               do 101 in = 1, mn
               do 101 il = 0, mm

                  if( tr(jr,jz,iz,in,il,1) .ne. 0.d0 ) then

                     if( itunt(m) .eq. 1 ) then

                        fmaxfc = tr(jr,jz,iz,in,il,1)

                     else if( itunt(m) .eq. 2 ) then

                        fmaxfc = tr(jr,jz,iz,in,il,1) / vl(jr,jz)

                     end if

                     if( facmax(m) .lt. fmaxfc) facmax(m) = fmaxfc

                  end if

  101          continue

            end if

            do 100 jr = 1, nr
            do 100 jz = 1, nz

               if( itunt(m) .eq. 1 ) then

                  cc = abs(rtfac(m)/facmax(m))

               else if( itunt(m) .eq. 2 ) then

                  cc = abs(rtfac(m)/facmax(m)) / vl(jr,jz)

               end if

            do 100 iz = 1, mz
            do 100 in = 1, mn
            do 100 il = 0, mm ! frtati 2022/02/18

               if( tr(jr,jz,iz,in,il,1) .ne. 0.d0 ) then

                  call calc_stdev(m,Xa,sigx,
     &                            tr(jr,jz,iz,in,il,1),
     &                            tr(jr,jz,iz,in,il,2),
     &                            cc)

                  tr(jr,jz,iz,in,il,1) = Xa
                  tr(jr,jz,iz,in,il,2) = sigx
                  if( sigx .gt. stdm ) stdm = sigx

                  if( tr(jr,jz,iz,in,il,1) .gt. cmax )
     &                                   cmax = tr(jr,jz,iz,in,il,1)

cABE 2022/03/02, avoid the negative value of cmin
                  if( tr(jr,jz,iz,in,il,1) .gt. 0.d0 .and.
     &                tr(jr,jz,iz,in,il,1) .lt. cmin)
     &                                   cmin = tr(jr,jz,iz,in,il,1)

               else

                  if( il .eq. 0 ) isdz = 1
                  tr(jr,jz,iz,in,il,2) = 0.0

               end if

  100       continue

            if( nobch .gt. ist_bat ) then
            if( rtstd(m) .gt. 0.d0 .and.
     &          stdm .gt. 0.d0 .and. isdz .eq. 0 .and.
     &          stdm .lt. rtstd(m) ) itstd(m) = 0
            end if

*-----------------------------------------------------------------------
*     out put unit = 31  : temporary number
*-----------------------------------------------------------------------

      do iax = 1, itfln(m)

!OBINATA(2012.8.20): output *.err
        noe = 1
        if ( any( itaxs(m,iax) .eq. (/ 8, 12 /) )
     &      .and. ittwo(m) .ne. 4 .and. icntl .ne. 8 )
     &    noe = 2

        do 900 ioe = 1, noe

         if( itall .eq. 2 .and. nobch .lt. maxbch ) then

            write(fnume,'(i3.3)') nobch
!OBINATA(2012.8.20): output *.err
            if ( ioe .eq. 1 ) then
              fname = ctfln(m,iax)(1:itfll(m,iax))//'.'//fnume
            else
              call mk_2derrfn(ctfln(m,iax),fname,itfll(m,iax))
            end if

         else if ( itall.eq.4 ) then
            write(fnume,'(i3.3)') nobch
            if ( npe.gt.1 ) write (fnume,'(i3.3)') nobch/(npe-1)
            if ( ioe.eq.1 ) then
              call mk_2dnumfn(ctfln(m,iax),fname,itfll(m,iax),
     &             nobch,maxbch,npe)
            else
              call mk_2dnerfn(ctfln(m,iax),fname,itfll(m,iax),fnume)
            end if

         else

!OBINATA(2012.8.20): output *.err
            if ( ioe .eq. 1 ) then
              fname = ctfln(m,iax)
            else
              call mk_2derrfn(ctfln(m,iax),fname,itfll(m,iax))
            end if

         end if

            iot = 31
            open(iot, file = fname, status = 'unknown' )

*-----------------------------------------------------------------------
            gfnam = ctfln(m,iax)
            igfmn = itfll(m,iax)
            igser = itger(m)

*-----------------------------------------------------------------------
*     input echo
*-----------------------------------------------------------------------

               call tyilech(iot,m,iax,1)

CCSE add for mesh=r-z parameter (2017.11.30) >>>>>
*-----------------------------------------------------------------------
*        dchain axis
*-----------------------------------------------------------------------

         if( itaxs(m,iax) .eq. 13 ) then

*-----------------------------------------------------------------------
*           output error or not
*-----------------------------------------------------------------------

            if( itout(m) .ne. 0 ) then

               do i = itfll(m,iax), 1, -1

                  if( ctfln(m,iax)(i:i) .eq. '.' ) goto 50

               end do

               i = itfll(m,iax)

   50          itfp = i - 1

               do i = 1, itfp

                  erfnm(i:i) = ctfln(m,iax)(i:i)

               end do

               erfnm(itfp+1:itfp+4) = '.err'

               iou = 15
               open(iou, file = erfnm(1:itfp+4), status = 'unknown' )

               call tyilech(iou,m,iax,1)

            end if

*-----------------------------------------------------------------------
*        output : do not change below expresion
*-----------------------------------------------------------------------

            write(iot,'(/"#",78("-"))')

            write(iot,'(/
     &           " r-z scoring mesh nuclear yield (or production)"/
     &           " ----------------------------------------")')

            if( itout(m) .ne. 0 ) then

               write(iou,'(/"#",78("-"))')

               write(iou,'(/
     &         " Statistical Error(%) for r-z scoring mesh nuclear ",
     &         "yield of above file."/
     &         " --------------------------------------------",
     &         "---------------------")')

            end if

*-----------------------------------------------------------------------

            do 170 il = 0, 2 ! frtati 2022/03/11
               do 170 iz = 1, maxpt

                  if( nn .gt. 0 ) then

                     do i = 1, nn

                        if( nt(i) / 1000 .eq. iz ) goto 150

                     end do

                     goto 170

                  end if

  150             continue

                  do in = 1, maxnt
                     do jr = 1, nr
                     do jz = 1, nz

                        if( itnzn(m).eq.0 ) then
                          t0_tr = tr(jr,jz,iz,in,il,1)
                        else
                          t0_tr = tr(jr,jz,igetiznmp(iz,in,il,m),1,0,1)
                        end if
                        if( t0_tr .ne. 0.0 ) goto 120

                     end do
                     end do
                  end do

                  goto 170

  120             im = in

                  do in = maxnt, im + 1, -1
                     do jr = 1, nr
                     do jz = 1, nz

                        if( itnzn(m).eq.0 ) then
                          t0_tr = tr(jr,jz,iz,in,il,1)
                        else
                          t0_tr = tr(jr,jz,igetiznmp(iz,in,il,m),1,0,1)
                        end if
                        if( t0_tr .ne. 0.0 ) goto 140

                     end do
                     end do
                  end do

                  jm = im
  140             jm = in

                  km = jm - im + 1
                  lm = ( km - 1 ) / ipstep + 1

                  do mmm = 1, lm

                     n1 = ipstep * ( mmm - 1 ) + 1
                     n2 = min( ipstep * mmm, km )
                     n3 = n1 + im - 1
                     n4 = n2 + im - 1

                     IF(il .eq. 0) then
                        write(iot,'(/1x,i4,"-",a2,
     &                           " isotope production",
     &                           " # n3 n4 = ",i3,1x,i3)')
     &                  iz, elmnt(iz), n3, n4

                     ELSEIF(il .eq. 1) then
                        write(iot,'(/1x,i4,"-",a2,
     &                           " 1st metastable isotope production",
     &                           " # n3 n4 = ",i3,1x,i3)')
     &                 iz, elmnt(iz), n3, n4

                     ELSEIF(il .eq. 2) then
                        write(iot,'(/1x,i4,"-",a2,
     &                           " 2nd metastable isotope production",
     &                           " # n3 n4 = ",i3,1x,i3)')
     &                 iz, elmnt(iz), n3, n4
                     ENDIF

                     IF(il .eq. 0) then
                        write(iot,'("     jr","     jz",12i11)')
     &                       ( im - 1 + iz + n, n = n1, n2 )
                     ELSEIF(il .eq. 1) then
                      write(iot,'("     jr","     jz"
     &                ,i10,"m",i10,"m",i10,"m",i10,"m",i10,"m"
     &                ,i10,"m",i10,"m",i10,"m",i10,"m",i10,"m"
     &                ,i10,"m",i10,"m")')
     &                ( im - 1 + iz + n, n = n1, n2 )
                     ELSEIF(il .eq. 2) then
                      write(iot,'("     jr","     jz"
     &                ,i10,"n",i10,"n",i10,"n",i10,"n",i10,"n"
     &                ,i10,"n",i10,"n",i10,"n",i10,"n",i10,"n"
     &                ,i10,"n",i10,"n")')
     &                ( im - 1 + iz + n, n = n1, n2 )
                     ENDIF

                     do jr = 1, nr
                     do jz = 1, nz
                      if( itnzn(m).eq.0 ) then ! frtati 2022/02/18
                        write(iot,'(2i7,1p12e11.3)')
     &                       jr,jz, (tr(jr,jz,iz,i,il,ioe), i = n3, n4 )
                      else
                        write(iot,'(2i7,1p12e11.3)')
     &                  jr,jz, (tr(jr,jz,igetiznmp(iz,i,il,m),1,0,ioe),
     &                  i = n3, n4 )
                      end if
                     end do
                     end do

                     if( itout(m) .ne. 0 ) then

                        IF(il .eq. 0) then
                           write(iou,'(/1x,i4,"-",a2,
     &                             " isotope production : ERROR(%)")')
     &                         iz, elmnt(iz)

                        ELSEIF(il .eq. 1) then
                           write(iou,'(/1x,i4,"-",a2,
     &              " 1st metastable isotope production : ERROR(%)")')
     &              iz, elmnt(iz)

                        ELSEIF(il .eq. 2) then
                           write(iou,'(/1x,i4,"-",a2,
     &              " 2nd metastable isotope production : ERROR(%)")')
     &              iz, elmnt(iz)
                        ENDIF

                        IF(il .eq. 0) then
                       write(iou,'("     jr","     jz",12i11)')
     &                      ( im - 1 + iz + n, n = n1, n2 )
                        ELSEIF(il .eq. 1) then
                       write(iou,'("     jr","     jz"
     &                ,i10,"m",i10,"m",i10,"m",i10,"m",i10,"m"
     &                ,i10,"m",i10,"m",i10,"m",i10,"m",i10,"m"
     &                ,i10,"m",i10,"m")')
     &                ( im - 1 + iz + n, n = n1, n2 )
                        ELSEIF(il .eq. 2) then
                       write(iou,'("     jr","     jz"
     &                ,i10,"n",i10,"n",i10,"n",i10,"n",i10,"n"
     &                ,i10,"n",i10,"n",i10,"n",i10,"n",i10,"n"
     &                ,i10,"n",i10,"n")')
     &                 ( im - 1 + iz + n, n = n1, n2 )
                        ENDIF

                        do jr = 1, nr
                        do jz = 1, nz
                         if( itnzn(m).eq.0 ) then ! frtati 2022/02/18
                           write(iou,'(2i7,1p12e11.3)')
     &                   ir,iz, (tr(jr,jz,iz,i,il,2)*100.0, i = n3, n4 )
                         else
                           write(iou,'(2i7,1p12e11.3)') ir,iz,
     &                     (tr(jr,jz,igetiznmp(iz,i,il,m),1,0,2)*100.0,
     &                     i = n3, n4 )
                         end if
                        end do
                        end do

                     end if

                  end do

  170       continue

            if( itout(m) .ne. 0 ) close(iou)

CCSE add for mesh=r-z parameter (2017.11.30) <<<<<

*-----------------------------------------------------------------------
*        r axis
*-----------------------------------------------------------------------

CCSE change for mesh=r-z parameter (2017.11.30) >>>>>
         else if( itaxs(m,iax) .eq. 6 ) then
CCSE change for mesh=r-z parameter (2017.11.30) <<<<<

               inum = 0

                  if( nn .eq. 0 ) then

                     nc = 1

                  else

                     nc = nn

                  end if

*-----------------------------------------------------------------------

            do jz = 1, nz

            do ic = 1, nc

               inum = inum + 1

*-----------------------------------------------------------------------

               if( nn .gt. 0 ) then

                  iz = nt(ic) / 1000

                  ia = nt(ic) - iz * 1000

                  if( ia .gt. 0 ) then

                     in = ia - iz

                     call chname(idum,ia,iz,chau)

                  end if

               end if

*-----------------------------------------------------------------------

                  write(iot,'(/"#",78("-"))')

               if( inum .eq. 1 ) then
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

               if( rtfac(m) .lt. 0.d0 )
     &           write(iot,
     &           '( "#restart maximum for normalization:",1p1e14.7)')
     &           facmax(m)

            if( nn .eq. 0 ) then

               write(iot,'("#   no. =",i3,3x,
     &         "iz  =",i3,3x,"Z = all"/
     &         "#   z = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     inum, jz,
     &                     zm(jz), zm(jz+1)

            else if( ia .eq. 0 ) then

               write(iot,'("#   no. =",i3,3x,
     &         "iz  =",i3,3x,"Z = ",i3," : ",a3/
     &         "#   z = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     inum, jz, iz, elmnt(iz),
     &                     zm(jz), zm(jz+1)

            else

               write(iot,'("#   no. =",i3,3x,
     &         "iz  =",i3,3x,a8/
     &         "#   z = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     inum, jz, chau,
     &                     zm(jz), zm(jz+1)

            end if

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: r [cm]")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: Number ",a15)') hsunit(itunt(m))

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

                  write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')

               if( itanl(m) .gt. 0 ) then

                  write(iot,'( "p: ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               end if

               if( itsans(m) .gt. 0 ) then

                  call write_sangel(iot,m,0)

               end if

c S.H. added IF statement below for epsout=2 (2017.2.18)
              if ( iteps(m) .ne. 2 ) then

               write(iot,'( "h: n",12x,"x",12x,
     &                    "y,hhl0      n")')

              else

               write(iot,'( "h: n",12x,"x",12x,
     &                    "y1,hhl0     ny2 dy1=[y1*y2]")')

              end if

               write(iot,'( "#  r-lower      r-upper",3x,
     &                    "   number      r.err")')

*-----------------------------------------------------------------------

                     seka = 0.0
                     sera = 0.0
                     seva = 0.0

               do jr = 1, nr

                     sek = 0.0
                     ser = 0.0

                     if( itunt(m) .eq. 1 ) then

                        vm = 1.0

                     else if( itunt(m) .eq. 2 ) then

                        vm = vl(jr,jz)

                     end if

*-----------------------------------------------------------------------

                  if( nn .eq. 0 ) then

                     do kn = 1, mn
                     do kz = 1, mz
                     do il = 0, mm ! frtati 2022/02/18

                        lz  = ikzz( kz, kn )
                        ln  = iknn( kz, kn )

                        vn  = vm * tr(jr,jz,lz,ln,il,1)
                        sek = sek + vn
                        ser = ser + ( vn * tr(jr,jz,lz,ln,il,2) )**2

                     end do
                     end do
                     end do

                  else if( ia .eq. 0 ) then

                     do kn = 1, mn
                     do il = 0, mm ! frtati 2022/02/18

                        lz  = ikzz( iz, kn )
                        ln  = iknn( iz, kn )

                        vn  = vm * tr(jr,jz,lz,ln,il,1)
                        sek = sek + vn
                        ser = ser + ( vn * tr(jr,jz,lz,ln,il,2) )**2

                     end do
                     end do

                  else

                        lz  = ikzz( iz, in )
                        ln  = iknn( iz, in )

                     do il = 0, mm ! frtati 2022/02/18
                        vn  = vm * tr(jr,jz,lz,ln,il,1)
                        sek = sek + vn
                        ser = ser
     &                      + ( vn * tr(jr,jz,lz,ln,il,2) )**2
                     end do

                  end if

*-----------------------------------------------------------------------

                        seka = seka + sek
                        sera = sera + ser
                        seva = seva + vm

                     if( sek .ne. 0.0 ) then ! .gt. -> .ne. by T.Sato 2024/05/02 for negative sek

                        ser = sqrt( ser ) / sek
                        sek = sek / vm

                     end if

                  write(iot,'(1p2e13.4,1pe13.4,0pf8.4)')
     &                  rm(jr),rm(jr+1), sek, ser

               end do

                     if( seka .gt. 0.0 ) then

                        sera = sqrt( sera ) / seka

                        if( itunt(m) .eq. 2 )
     &                  seka = seka / seva

                     end if

               write(iot,'(/"#   sum over",14x,1pe13.4,0pf8.4)')
     &               seka, sera

            if( itunt(m) .eq. 2 ) then

               write(iot,'(
     &         "#   vol sum = ",1p1e13.4,"  [cm^3]")') seva

            end if

               if( nn .eq. 0 ) then

                  write(iot,'(/a1,"no. =",i3,
     &                        ",    iz =",i3,
     &                        ",    Z  =  all", a1)')
     &                        cha, inum, jz, cha

               else if( ia .eq. 0 ) then

                  write(iot,'(/a1,"no. =",i3,
     &                        ",    iz =",i3,
     &                        ",    Z  =  ",i3, a1)')
     &                        cha, inum, jz, iz, cha

               else

                  write(iot,'(/a1,"no. =",i3,
     &                        ",    iz =",i3,
     &                        ",    ",a8, a1)')
     &                        cha, inum, jz, chau, cha

               end if

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

               if( nn .eq. 0 ) then

                  write(iot,'("wt: s(0.7)",/
     &                        a1,"vspace{-3}"/
     &                        "  zmin  &=&",1pe13.4," [cm]"/
     &                        "  zmax  &=&",1pe13.4," [cm]"/
     &                        "      Z &=&   all"/
     &                        "e:")')
     &                        yen,
     &                        zm(jz), zm(jz+1)

               else if( ia .eq. 0 ) then

                  write(iot,'("wt: s(0.7)",/
     &                        a1,"vspace{-3}"/
     &                        "  zmin  &=&",1pe13.4," [cm]"/
     &                        "  zmax  &=&",1pe13.4," [cm]"/
     &                        "      Z &=& ",i3,"   :  ",a3/
     &                        "e:")')
     &                        yen,
     &                        zm(jz), zm(jz+1), iz, elmnt(iz)

               else

                  write(iot,'("wt: s(0.7)",/
     &                        a1,"vspace{-3}"/
     &                        "  zmin  &=&",1pe13.4," [cm]"/
     &                        "  zmax  &=&",1pe13.4," [cm]"/
     &                        "  nucl. &=&   ",a8/
     &                        "e:")')
     &                        yen,
     &                        zm(jz), zm(jz+1), chau

               end if

*-----------------------------------------------------------------------

            end do
            end do

*-----------------------------------------------------------------------
*        z axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 5 ) then

               inum = 0

                  if( nn .eq. 0 ) then

                     nc = 1

                  else

                     nc = nn

                  end if

*-----------------------------------------------------------------------

            do jr = 1, nr

            do ic = 1, nc

               inum = inum + 1

*-----------------------------------------------------------------------

               if( nn .gt. 0 ) then

                  iz = nt(ic) / 1000

                  ia = nt(ic) - iz * 1000

                  if( ia .gt. 0 ) then

                     in = ia - iz

                     call chname(idum,ia,iz,chau)

                  end if

               end if

*-----------------------------------------------------------------------

                  write(iot,'(/"#",78("-"))')

               if( inum .eq. 1 ) then
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

               if( rtfac(m) .lt. 0.d0 )
     &           write(iot,
     &           '( "#restart maximum for normalization:",1p1e14.7)')
     &           facmax(m)

            if( nn .eq. 0 ) then

               write(iot,'("#   no. =",i3,3x,
     &         "ir  =",i3,3x,"Z = all"/
     &         "#   r = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     inum, jr,
     &                     rm(jr), rm(jr+1)

            else if( ia .eq. 0 ) then

               write(iot,'("#   no. =",i3,3x,
     &         "ir  =",i3,3x,"Z = ",i3," : ",a3/
     &         "#   r = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     inum, jr, iz, elmnt(iz),
     &                     rm(jr), rm(jr+1)

            else

               write(iot,'("#   no. =",i3,3x,
     &         "ir  =",i3,3x,a8/
     &         "#   r = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     inum, jr, chau,
     &                     rm(jr), rm(jr+1)

            end if

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: z [cm]")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: Number ",a15)') hsunit(itunt(m))

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

                  write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')

               if( itanl(m) .gt. 0 ) then

                  write(iot,'( "p: ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               end if

               if( itsans(m) .gt. 0 ) then

                  call write_sangel(iot,m,0)

               end if

c S.H. added IF statement below for epsout=2 (2017.2.18)
              if ( iteps(m) .ne. 2 ) then

               write(iot,'( "h: n",12x,"x",12x,
     &                    "y,hhl0      n")')

              else

               write(iot,'( "h: n",12x,"x",12x,
     &                    "y1,hhl0     ny2 dy1=[y1*y2]")')

              end if

               write(iot,'( "#  z-lower      z-upper",3x,
     &                    "   number      r.err")')

*-----------------------------------------------------------------------

                     seka = 0.0
                     sera = 0.0
                     seva = 0.0

               do jz = 1, nz

                     sek = 0.0
                     ser = 0.0

                     if( itunt(m) .eq. 1 ) then

                        vm = 1.0

                     else if( itunt(m) .eq. 2 ) then

                        vm = vl(jr,jz)

                     end if

*-----------------------------------------------------------------------

                  if( nn .eq. 0 ) then

                     do kn = 1, mn
                     do kz = 1, mz
                     do il = 0, mm ! frtati 2022/02/18

                        lz  = ikzz( kz, kn )
                        ln  = iknn( kz, kn )

                        vn  = vm * tr(jr,jz,lz,ln,il,1)
                        sek = sek + vn
                        ser = ser + ( vn * tr(jr,jz,lz,ln,il,2) )**2

                     end do
                     end do
                     end do

                  else if( ia .eq. 0 ) then

                     do kn = 1, mn
                     do il = 0, mm ! frtati 2022/02/18

                        lz  = ikzz( iz, kn )
                        ln  = iknn( iz, kn )

                        vn  = vm * tr(jr,jz,lz,ln,il,1)
                        sek = sek + vn
                        ser = ser + ( vn * tr(jr,jz,lz,ln,il,2) )**2

                     end do
                     end do

                  else

                        lz  = ikzz( iz, in )
                        ln  = iknn( iz, in )

                     do il = 0, mm ! frtati 2022/02/18
                        vn  = vm * tr(jr,jz,lz,ln,il,1)
                        sek = sek + vn
                        ser = ser
     &                      + ( vn * tr(jr,jz,lz,ln,il,2) )**2
                     end do

                  end if

*-----------------------------------------------------------------------

                        seka = seka + sek
                        sera = sera + ser
                        seva = seva + vm

                     if( sek .ne. 0.0 ) then ! .gt. -> .ne. byT.Sato 2024/05/02 for negative sek

                        ser = sqrt( ser ) / sek
                        sek = sek / vm

                     end if

                  write(iot,'(1p2e13.4,1pe13.4,0pf8.4)')
     &                  zm(jz),zm(jz+1), sek, ser

               end do

                     if( seka .gt. 0.0 ) then

                        sera = sqrt( sera ) / seka

                        if( itunt(m) .eq. 2 )
     &                  seka = seka / seva

                     end if

               write(iot,'(/"#   sum over",14x,1pe13.4,0pf8.4)')
     &               seka, sera

            if( itunt(m) .eq. 2 ) then

               write(iot,'(
     &         "#   vol sum = ",1p1e13.4,"  [cm^3]")') seva

            end if

               if( nn .eq. 0 ) then

                  write(iot,'(/a1,"no. =",i3,
     &                        ",    ir =",i3,
     &                        ",    Z  =  all", a1)')
     &                        cha, inum, jr, cha

               else if( ia .eq. 0 ) then

                  write(iot,'(/a1,"no. =",i3,
     &                        ",    ir =",i3,
     &                        ",    Z  =  ",i3, a1)')
     &                        cha, inum, jr, iz, cha

               else

                  write(iot,'(/a1,"no. =",i3,
     &                        ",    ir =",i3,
     &                        ",    ",a8, a1)')
     &                        cha, inum, jr, chau, cha

               end if

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

               if( nn .eq. 0 ) then

                  write(iot,'("wt: s(0.7)",/
     &                        a1,"vspace{-3}"/
     &                        "  rmin  &=&",1pe13.4," [cm]"/
     &                        "  rmax  &=&",1pe13.4," [cm]"/
     &                        "      Z &=&   all"/
     &                        "e:")')
     &                        yen,
     &                        rm(jr), rm(jr+1)

               else if( ia .eq. 0 ) then

                  write(iot,'("wt: s(0.7)",/
     &                        a1,"vspace{-3}"/
     &                        "  rmin  &=&",1pe13.4," [cm]"/
     &                        "  rmax  &=&",1pe13.4," [cm]"/
     &                        "      Z &=& ",i3,"   :  ",a3/
     &                        "e:")')
     &                        yen,
     &                        rm(jr), rm(jr+1), iz, elmnt(iz)

               else

                  write(iot,'("wt: s(0.7)",/
     &                        a1,"vspace{-3}"/
     &                        "  rmin  &=&",1pe13.4," [cm]"/
     &                        "  rmax  &=&",1pe13.4," [cm]"/
     &                        "  nucl. &=&   ",a8/
     &                        "e:")')
     &                        yen,
     &                        rm(jr), rm(jr+1), chau

               end if

*-----------------------------------------------------------------------

            end do
            end do

*-----------------------------------------------------------------------
*        mass axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 1 ) then

               inum = 0

                  if( nn .eq. 0 ) then

                     nc = 1

                  else

                     nc = nn

                  end if

*-----------------------------------------------------------------------

            do jr = 1, nr
            do jz = 1, nz

            do 260 ic = 1, nc

                  if( nn .gt. 0 ) iz = nt(ic) / 1000

                  do i = 1, maxnt + maxpt

                     tm(i,1) = 0.0d+0
                     tm(i,2) = 0.0d+0

                  end do

                  if( itunt(m) .eq. 1 ) then

                     vm = 1.0

                  else if( itunt(m) .eq. 2 ) then

                     vm = vl(jr,jz)

                  end if

*-----------------------------------------------------------------------

                  if( nn .eq. 0 .and. itfln(m) .eq. 1 ) then

                     do ln = 1, mn
                     do il = 0, mm ! frtati 2022/02/18

                        lz = 1
                        im = ln

                        vn = vm * tr(jr,jz,lz,ln,il,1)
                        tm(im,1) = tm(im,1) + vn
                        tm(im,2) = tm(im,2)
     &                           + ( vn * tr(jr,jz,lz,ln,il,2) )**2

                     end do
                     end do

                  else if( nn .eq. 0 .and. itfln(m) .gt. 1 ) then

                     do kn = 1, mn
                     do kz = 1, mz
                     do il = 0, mm ! frtati 2022/02/18

                        lz = ikzz( kz, kn )
                        ln = iknn( kz, kn )

                        im = kz + kn

                        vn = vm * tr(jr,jz,lz,ln,il,1)
                        tm(im,1) = tm(im,1) + vn
                        tm(im,2) = tm(im,2)
     &                           + ( vn * tr(jr,jz,lz,ln,il,2) )**2

                     end do
                     end do
                     end do

                  else

                     do kn = 1, mn
                     do il = 0, mm ! frtati 2022/02/18

                        lz = ikzz( iz, kn )
                        ln = iknn( iz, kn )

                        im = iz + kn

                        vn = vm * tr(jr,jz,lz,ln,il,1)
                        tm(im,1) = tm(im,1) + vn
                        tm(im,2) = tm(im,2)
     &                           + ( vn * tr(jr,jz,lz,ln,il,2) )**2

                     end do
                     end do

                  end if

*-----------------------------------------------------------------------

               do i = 1, maxnt + maxpt

                  if( tm(i,1) .gt. 0.0d0 ) then

                     tm(i,2) = sqrt( tm(i,2) ) / tm(i,1)
                     tm(i,1) = tm(i,1) / vm

                  end if

               end do

               do i = 1, maxnt + maxpt

                  if( tm(i,1) .gt. 0.0d0 ) goto 200

               end do

                  goto 260

  200          im = i

               do i = maxnt + maxpt, im + 1, -1

                  if( tm(i,1) .gt. 0.0d0 ) goto 210

               end do

               jm = im
  210          jm = i

               im = max( 1, im - 2 )
               jm = jm + 2

*-----------------------------------------------------------------------

                  write(iot,'(/"#",78("-"))')

               inum = inum + 1

               if( inum .eq. 1 ) then
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

               if( rtfac(m) .lt. 0.d0 )
     &           write(iot,
     &           '( "#restart maximum for normalization:",1p1e14.7)')
     &           facmax(m)

               if( nn .eq. 0 ) then

                  write(iot,'("#   no. =",i3,3x,
     &            "ir =",i3,3x,"iz =",i3,3x,
     &            "Z = all"/
     &            "#   r = (",1p1e13.4,"  -",1p1e13.4,"  )"/
     &            "#   z = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        inum, jr, jz,
     &                        rm(jr), rm(jr+1),
     &                        zm(jz), zm(jz+1)

               else

                  write(iot,'("#   no. =",i3,3x,
     &            "ir =",i3,3x,"iz =",i3,3x,"Z = ",i3," : ",a3/
     &            "#   r = (",1p1e13.4,"  -",1p1e13.4,"  )"/
     &            "#   z = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        inum, jr, jz, iz, elmnt(iz),
     &                        rm(jr), rm(jr+1),
     &                        zm(jz), zm(jz+1)

               end if

               write(iot,'("# im jm = ",i3,1x,i3)') im, jm

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: Mass")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: Number ",a15)') hsunit(itunt(m))

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

                  write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')

               if( itanl(m) .gt. 0 ) then

                  write(iot,'( "p: ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               end if

               if( itsans(m) .gt. 0 ) then

                  call write_sangel(iot,m,0)

               end if

c S.H. added IF statement below for epsout=2 (2017.2.18)
              if ( iteps(m) .ne. 2 ) then

               write(iot,'( "h: x-0.5    y,hl0       n")')

              else

               write(iot,'( "h: x-0.5    y1,hl0      ",
     &                      "ny2 dy1=[y1*y2]")')

              end if

               write(iot,'( "#  mass     number      r.err")')

*-----------------------------------------------------------------------

               if( im .eq. 1 ) then

                  write(iot,'(3x,f4.1,2x,1pe13.4,0pf8.4)')
     &                                   0.5, 0.0, 0.0

               end if

                  seka = 0.0
                  sera = 0.0

               do i = im, jm

                  vn  = vm * tm(i,1)
                  seka = seka + vn
                  sera = sera + ( vn * tm(i,2) )**2

                  write(iot,'(3x,i4,2x,1pe13.4,0pf8.4)')
     &                           i, tm(i,1), tm(i,2)

               end do

                  if( seka .gt. 0.0 ) then

                     sera = sqrt( sera ) / seka
                     seka = seka / vm

                  end if

               write(iot,'(/"#   sum  ",1pe13.4,0pf8.4)')
     &               seka, sera

               if( itunt(m) .eq. 2 ) then

                  write(iot,'(
     &            "#   vol sum = ",1p1e13.4,"  [cm^3]")') vm

               end if

               if( nn .eq. 0 ) then

                  write(iot,'(/a1,"no. =",i3,
     &                            ",    ir  =  ",i3,
     &                            ",    iz  =  ",i3,
     &                            ",    Z  =  all", a1)')
     &                        cha, inum, jr, jz, cha

               else

                  write(iot,'(/a1,"no. =",i3,
     &                            ",    ir  =  ",i3,
     &                            ",    iz  =  ",i3,
     &                            ",    Z  =  ",i3, a1)')
     &                        cha, inum, jr, jz, iz, cha

               end if

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

               if( itunt(m) .eq. 1 ) then

                  if( nn .eq. 0 ) then

                     write(iot,'("wt: s(0.7)",/
     &                           a1,"vspace{-3}"/
     &                           "  mass distribution"/
     &                           "  in r-z mesh"/
     &                           "  rmin  &=&",1pe13.4," [cm]"/
     &                           "  rmax  &=&",1pe13.4," [cm]"/
     &                           "  zmin  &=&",1pe13.4," [cm]"/
     &                           "  zmax  &=&",1pe13.4," [cm]"/
     &                           "      Z &=&   all"/
     &                           "e:")')
     &                           yen,
     &                           rm(jr), rm(jr+1),
     &                           zm(jz), zm(jz+1)

                  else

                     write(iot,'("wt: s(0.7)",/
     &                           a1,"vspace{-3}"/
     &                           "  isotope distribution"/
     &                           "  in r-z mesh"/
     &                           "  rmin  &=&",1pe13.4," [cm]"/
     &                           "  rmax  &=&",1pe13.4," [cm]"/
     &                           "  zmin  &=&",1pe13.4," [cm]"/
     &                           "  zmax  &=&",1pe13.4," [cm]"/
     &                           "      Z &=& ",i3,"   :  ",a3/
     &                           "e:")')
     &                           yen,
     &                           rm(jr), rm(jr+1),
     &                           zm(jz), zm(jz+1),
     &                           iz, elmnt(iz)

                  end if

               else if( itunt(m) .eq. 2 ) then

                  if( nn .eq. 0 ) then

                     write(iot,'("wt: s(0.7)",/
     &                           a1,"vspace{-3}"/
     &                           "  mass distribution"/
     &                           "  in r-z mesh"/
     &                           "  rmin  &=&",1pe13.4," [cm]"/
     &                           "  rmax  &=&",1pe13.4," [cm]"/
     &                           "  zmin  &=&",1pe13.4," [cm]"/
     &                           "  zmax  &=&",1pe13.4," [cm]"/
     &                           "    vol &=&",1pe13.4," [cm^3]"/
     &                           "      Z &=&   all"/
     &                           "e:")')
     &                           yen,
     &                           rm(jr), rm(jr+1),
     &                           zm(jz), zm(jz+1),
     &                           vl(jr,jz)

                  else

                     write(iot,'("wt: s(0.7)",/
     &                           a1,"vspace{-3}"/
     &                           "  isotope distribution"/
     &                           "  in r-z mesh"/
     &                           "  rmin  &=&",1pe13.4," [cm]"/
     &                           "  rmax  &=&",1pe13.4," [cm]"/
     &                           "  zmin  &=&",1pe13.4," [cm]"/
     &                           "  zmax  &=&",1pe13.4," [cm]"/
     &                           "    vol &=&",1pe13.4," [cm^3]"/
     &                           "      Z &=& ",i3,"   :  ",a3/
     &                           "e:")')
     &                           yen,
     &                           rm(jr), rm(jr+1),
     &                           zm(jz), zm(jz+1),
     &                           vl(jr,jz),
     &                           iz, elmnt(iz)

                  end if

               end if

*-----------------------------------------------------------------------

  260       continue

            end do
            end do

*-----------------------------------------------------------------------
*        charge axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 7 ) then

               inum = 0

            do 380 jr = 1, nr
            do 370 jz = 1, nz

                  do i = 1, maxpt

                     tm(i,1) = 0.0d+0
                     tm(i,2) = 0.0d+0

                  end do

                  if( itunt(m) .eq. 1 ) then

                     vm = 1.0

                  else if( itunt(m) .eq. 2 ) then

                     vm = vl(jr,jz)

                  end if

*-----------------------------------------------------------------------

                  do kz = 1, mz
                  do kn = 1, mn
                  do il = 0, mm ! frtati 2022/02/18

                        iz  = ikzz( kz, kn )
                        in  = iknn( kz, kn )

                        vn  = vm * tr(jr,jz,iz,in,il,1)
                        tm(iz,1) = tm(iz,1) + vn
                        tm(iz,2) = tm(iz,2)
     &                           + ( vn * tr(jr,jz,iz,in,il,2) )**2

                  end do
                  end do
                  end do

*-----------------------------------------------------------------------

               do i = 1, maxpt

                  if( tm(i,1) .gt. 0.0d0 ) then

                     tm(i,2) = sqrt( tm(i,2) ) / tm(i,1)
                     tm(i,1) = tm(i,1) / vm

                  end if

               end do

               do i = 1, maxpt

                  if( tm(i,1) .gt. 0.0d0 ) goto 300

               end do

               goto 370

  300          im = i

               do i = maxpt, im + 1, -1

                  if( tm(i,1) .gt. 0.0d0 ) goto 310

               end do

               jm = im
  310          jm = i

               im = max( 1, im - 2 )
               jm = jm + 2

*-----------------------------------------------------------------------

                  write(iot,'(/"#",78("-"))')

               inum = inum + 1

               if( inum .eq. 1 ) then
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

               if( rtfac(m) .lt. 0.d0 )
     &           write(iot,
     &           '( "#restart maximum for normalization:",1p1e14.7)')
     &           facmax(m)

               write(iot,'("#   no. =",i3,3x,
     &         "ir =",i3,3x,"iz =",i3,3x/
     &         "#   r = (",1p1e13.4,"  -",1p1e13.4,"  )"/
     &         "#   z = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     inum, jr, jz,
     &                     rm(jr), rm(jr+1),
     &                     zm(jz), zm(jz+1)

               write(iot,'("# im jm = ",i3,1x,i3)') im, jm

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: Charge")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: Number ",a15)') hsunit(itunt(m))

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

                  write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')

               if( itanl(m) .gt. 0 ) then

                  write(iot,'( "p: ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               end if

               if( itsans(m) .gt. 0 ) then

                  call write_sangel(iot,m,0)

               end if

c S.H. added IF statement below for epsout=2 (2017.2.18)
              if ( iteps(m) .ne. 2 ) then

               write(iot,'( "h: x-0.5    y,hl0       n")')

              else

               write(iot,'( "h: x-0.5    y1,hl0      ",
     &                      "ny2 dy1=[y1*y2]")')

              end if

               write(iot,'( "#  charge   number      r.err")')

*-----------------------------------------------------------------------

               if( im .eq. 1 ) then

                  write(iot,'(3x,f4.1,2x,1pe13.4,0pf8.4)')
     &                                   0.5, 0.0, 0.0

               end if

                  seka = 0.0
                  sera = 0.0

               do i = im, jm

                  vn  = vm * tm(i,1)
                  seka = seka + vn
                  sera = sera + ( vn * tm(i,2) )**2

                  write(iot,'(3x,i4,2x,1pe13.4,0pf8.4)')
     &                           i, tm(i,1), tm(i,2)

               end do

                  if( seka .gt. 0.0 ) then

                     sera = sqrt( sera ) / seka
                     seka = seka / vm

                  end if

               write(iot,'(/"#   sum  ",1pe13.4,0pf8.4)')
     &               seka, sera

               if( itunt(m) .eq. 2 ) then

                  write(iot,'(
     &            "#   vol sum = ",1p1e13.4,"  [cm^3]")') vm

               end if

                  write(iot,'(/a1,"no. =",i3,
     &                            ",    ir  =  ",i3,
     &                            ",    iz  =  ",i3, a1)')
     &                        cha, inum, jr, jz, cha

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

               if( itunt(m) .eq. 1 ) then

                     write(iot,'("wt: s(0.7)",/
     &                           a1,"vspace{-3}"/
     &                           "  charge distribution"/
     &                           "  in r-z mesh"/
     &                           "  rmin  &=&",1pe13.4," [cm]"/
     &                           "  rmax  &=&",1pe13.4," [cm]"/
     &                           "  zmin  &=&",1pe13.4," [cm]"/
     &                           "  zmax  &=&",1pe13.4," [cm]"/
     &                           "e:")')
     &                           yen,
     &                           rm(jr), rm(jr+1),
     &                           zm(jz), zm(jz+1)

               else if( itunt(m) .eq. 2 ) then

                     write(iot,'("wt: s(0.7)",/
     &                           a1,"vspace{-3}"/
     &                           "  charge distribution"/
     &                           "  in r-z mesh"/
     &                           "  rmin  &=&",1pe13.4," [cm]"/
     &                           "  rmax  &=&",1pe13.4," [cm]"/
     &                           "  zmin  &=&",1pe13.4," [cm]"/
     &                           "  zmax  &=&",1pe13.4," [cm]"/
     &                           "    vol &=&",1pe13.4," [cm^3]"/
     &                           "e:")')
     &                           yen,
     &                           rm(jr), rm(jr+1),
     &                           zm(jz), zm(jz+1),
     &                           vl(jr,jz)

               end if

  370       continue
  380       continue

*-----------------------------------------------------------------------
*        chart axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 8 ) then

               inum = 0

            do 480 jr = 1, nr
            do 470 jz = 1, nz
            do 460 il = 0, 2 ! OGW loop for isomeric level

*-----------------------------------------------------------------------

                  itmax  = 0
                  icmax  = 0
                  inmax  = 0

               do iz = 1, maxpt
               do in = 1, maxnt

                  if( itnzn(m).eq.0 ) then
                    t0_tr = tr(jr,jz,iz,in,il,1)
                  else
                    t0_tr = tr(jr,jz,igetiznmp(iz,in,il,m),1,0,1)
                  end if
                  if( t0_tr .gt. 0.d0 ) then

                     if( iz + in .gt. itmax ) itmax = iz + in
                     if( iz      .gt. icmax ) icmax = iz
                     if(      in .gt. inmax ) inmax =      in

                  end if

               end do
               end do

                     dxmax = dble(inmax+2)
                     dymax = dble(icmax+2)
                     dform = dymax / dxmax

                     inmag = inmax

                  if( dform .gt. 0.8 ) then

                     inmag = nint( dymax / 0.8 )
                     dxmax = dble( nint( dymax / 0.8 ) )
                     dform = dymax / dxmax

                  end if

               if( itmax .eq. 0 ) goto 470

*-----------------------------------------------------------------------

                  write(iot,'(/"#",78("-"))')

               inum = inum + 1

               if( inum .eq. 1 .and. il .eq. 0 ) then
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

               if( rtfac(m) .lt. 0.d0 )
     &           write(iot,
     &           '( "#restart maximum for normalization:",1p1e14.7)')
     &           facmax(m)

               write(iot,'("#   no. =",i3,3x,
     &         "ir =",i3,3x,"iz  =",i3,3x,"il  =",i3)')
     &                     inum, jr, jz, il

               write(iot,'("# icmax inmax = ",i3,1x,i3)')
     &            icmax, inmax

*-----------------------------------------------------------------------

            if( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) then

                  write(iot,'(/a1,"no. =",i3,
     &                            ",    ir  =  ",i3,
     &                            ",    iz  =  ",i3, a1)')
     &                        cha, inum, jr, jz, cha

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

             if(il .eq. 0) then    ! specify isomeric level at the bottom of figure
               write(iot,'("msdc: {",a1,"huge Ground state}")') yen
             elseif(il .eq. 1) then
               write(iot,'("msdc: {",a1,"huge 1st isomer}")') yen
             else
               write(iot,'("msdc: {",a1,"huge 2nd isomer}")') yen
             endif

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: N Neutron Number")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: Z Proton Number")')

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

*-----------------------------------------------------------------------

               if( inum .eq. 1 ) then

                  form  = dform
                  xfac  = 1.1
                  afac  = 0.6
                  izlog = 1
                  inocm = 1
                  inolg = 1

                  if( itanl(m) .gt. 0 )
     &            call anset(itang(m),itanl(m),
     &                       form,xfac,afac,cmin,cmax,izlog,inocm,inolg)

                  if( form .le. 1.0 ) then

                     scal = form**0.35
                     xfac = xfac / form**0.5
                     xorg = 0.0
                     yorg = min( 20.d0, ( 1.0 / form - 1.0 ) / 2.5 )

                  else

                     scal = 1.0 / form**0.41
                     xfac = xfac / form**0.5
                     xorg = min( 20.d0, ( form - 1.0 ) / 1.5)
                     yorg = 0.0

                  end if

               end if

                  write(iot,'( "set: c1[",f6.3,
     &                            "] c2[",f6.3,
     &                            "] c5[",f6.3,
     &                            "] c6[",f6.3,
     &                            "] c7[",f6.3,
     &                            "] c8[",f6.3,"]")')
     &                    form, xfac, afac, scal, xorg, yorg
                  write(iot,'( "p: h2fs form[c1] xfac[c2]",
     &                 " afac[c5] scal[c6] xorg[c7] yorg[c8] nosp")')

               if( ( ( ittwo(m) .ge. 2 .and. ittwo(m) .le. 3 ) .or.
     &               ( ittwo(m) .ge. 6 .and. ittwo(m) .le. 7 ) ) .and.
     &             cmin .gt. 0.0 .and. cmax .gt. cmin ) then

                if (ioe .eq. 1 ) then
                  write(iot,'( "set: c3[",1p1e15.6,"] c4[",
     &                  1p1e15.6,"]")') cmin, cmax
                else
                  write(iot,'( "set: c3[1.0e-4] c4[1.0]")')
                end if
                  write(iot,'( "p: cmin[c3] cmax[c4]")')
                  write(iot,'( "p: dmin(1e-31)")')

                  if( ittwo(m) .eq. 3 .or. ittwo(m) .eq. 7 )
     &                 write(iot,'( "p: zlog")')

               end if

               write(iot,'( "p: xmin(0) xmax(",f5.1,")")') dxmax
               write(iot,'( "p: ymin(0) ymax(",f5.1,")")') dymax

               if( itanl(m) .gt. 0 .and. ioe .eq. 1 ) then

                  write(iot,'( "p: ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               else if( itanl(m) .gt. 0 .and. ioe .ne. 1 ) then

                   call terrang(iot,m)

               end if

               if( itsans(m) .gt. 0 ) then

                  call write_sangel(iot,m,0)

               end if

            end if

*-----------------------------------------------------------------------

            if( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) then

               write(iot,'(/"hc: y = ",i3," to 1 by -1 ;",
     &                         " x = 1 to ",i3," by 1 ;")')
     &                        icmax+2, inmax+2

               do i = icmax+2, 1, -1
                if( itnzn(m).eq.0 ) then ! frtati 2022/02/18
                  write(iot,'(1p10e11.3)')
     &                 ( tr(jr,jz,i,l,il,ioe), l = 1, inmax+2 )
                else
                  write(iot,'(1p10e11.3)')
     &          ( tr(jr,jz,igetiznmp(i,l,il,m),1,0,ioe),l = 1,inmax+2 )
                end if

               end do

*-----------------------------------------------------------------------

            else if( ittwo(m) .eq. 4 ) then

               write(iot,'(/"#   Z    N  Mass  ",
     &                      "  number    r.err")')

               do i = 1, icmax
               do j = 1, inmax

                if( itnzn(m).eq.0 ) then ! frtati 2022/02/18
                  if( tr(jr,jz,i,j,il,1) .gt. 0.d0 ) then

                     write(iot,'(3i5,1pe13.4,0pf8.4)')
     &               i, j, i+j,
     &               tr(jr,jz,i,j,il,1), tr(jr,jz,i,j,il,2)

                  end if
                else
                  if( tr(jr,jz,igetiznmp(i,j,il,m),1,0,1).gt.0.d0 ) then
                     write(iot,'(3i5,1pe13.4,0pf8.4)')
     &               i, j, i+j,
     &               tr(jr,jz,igetiznmp(i,j,il,m),1,0,1),
     &               tr(jr,jz,igetiznmp(i,j,il,m),1,0,2)
                  end if
                end if

               end do
               end do

*-----------------------------------------------------------------------

            else if( ittwo(m) .eq. 5 ) then

               write(iot,'(/"#   Z = 1 to ",i3/
     &                      "#   N = 1 to ",i3)')
     &                        icmax+2, inmax+2

                  write(iot,'(/4x,a3,4x,1p1000e11.3)')
     &            'Z/N', ( dble( in ), in = 1, inmax+2 )

               do i = icmax+2, 1, -1
                if( itnzn(m).eq.0 ) then ! frtati 2022/02/18
                  write(iot,'(1p1000e11.3)')
     &            dble( i ),
     &            ( tr(jr,jz,i,l,il,ioe), l = 1, inmax+2 )
                else
                  write(iot,'(1p1000e11.3)')
     &            dble( i ),
     &            ( tr(jz,jr,igetiznmp(i,l,il,m),1,0,ioe),
     &              l = 1,inmax+2 )
                end if

               end do

            end if

*-----------------------------------------------------------------------

            if( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) then

               if( itout(m) .ne. 0 ) then

                  call wmgcstb(iot,icmax,inmag)

               end if

            end if

*-----------------------------------------------------------------------

      if( inolg .eq. 1 .and. ( ittwo(m) .eq. 3 .or. ittwo(m) .eq. 7 )
     &    .and. cmin .gt. 0.0 .and. cmax .gt. cmin .and.
     &    igsh .eq. 0 ) then

      write(iot,'(/"#",78("-"))')
      write(iot,'(
     &"z: xorg(1.03)"/
     &"p: xfac[c2*0.03] form[c1/0.03] notc noxt noyt dmin(1)"/
     &"p: ymin(0) ymax(1) xmin(0) xmax(1) cmin(1) cmax(100) nosp"/
     &"hc: y= 0.005 to 0.995 by 0.01 ; x= 0.5 to 0.5 by 1 ;"/
     &" 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20"/
     &"21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40"/
     &"41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60"/
     &"61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80"/
     &"81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100"/
     &"z: xorg(0.0)"/
     &"p: xfac[c2*0.03] form[c1/0.03] ",
     &"nosp afac[c5*0.7] ymin[c3] ymax[c4]"/
     &"p: noxt noxn ytxt(-1) ynum(-1) xltd(-1) xstd(-1) itic(1) ylog"
     &)')

      if ( ioe .eq. 2 ) then

          write(iot,'("y: Relative Error")')

      else if( itazl(m) .eq. 0 ) then

         write(iot,'("y: Number ",a15)') hsunit(itunt(m))

      else

         write(iot,'("y: ",200a1)') (itazt(m)(i:i),i=1,itazl(m))

      end if

      if( izlog .eq. 0 ) write(iot,'( "p: ylin")')

      write(iot,'( "#",78("-"))')

      end if

*-----------------------------------------------------------------------

            if( inocm .eq. 1 ) then

               if( itunt(m) .eq. 1 ) then

                     write(iot,'(/"wt: s(0.7)",/
     &                           a1,"vspace{-3}"/
     &                           "  yield distribution"/
     &                           "  in r-z mesh"/
     &                           "  rmin  &=&",1pe13.4," [cm]"/
     &                           "  rmax  &=&",1pe13.4," [cm]"/
     &                           "  zmin  &=&",1pe13.4," [cm]"/
     &                           "  zmax  &=&",1pe13.4," [cm]"/
     &                           "e:")')
     &                           yen,
     &                           rm(jr), rm(jr+1),
     &                           zm(jz), zm(jz+1)

               else if( itunt(m) .eq. 2 ) then

                     write(iot,'(/"wt: s(0.7)",/
     &                           a1,"vspace{-3}"/
     &                           "  yield distribution"/
     &                           "  in r-z mesh"/
     &                           "  rmin  &=&",1pe13.4," [cm]"/
     &                           "  rmax  &=&",1pe13.4," [cm]"/
     &                           "  zmin  &=&",1pe13.4," [cm]"/
     &                           "  zmax  &=&",1pe13.4," [cm]"/
     &                           "    vol &=&",1pe13.4," [cm^3]"/
     &                           "e:")')
     &                           yen,
     &                           rm(jr), rm(jr+1),
     &                           zm(jz), zm(jz+1),
     &                           vl(jr,jz)

               end if

            end if

*-----------------------------------------------------------------------

      if( inolg .eq. 1 .and. ( ittwo(m) .eq. 3 .or. ittwo(m) .eq. 7 )
     &    .and. cmin .gt. 0.0 .and. cmax .gt. cmin .and.
     &    igsh .eq. 0 ) then

         write(iot,'(/"#",78("-"))')
         write(iot,'("z: xorg[-1.03/0.05]"/
     &               "p: form[c1] xfac[c2] afac[c5] nosp notf")')

               write(iot,'( "p: xmin(0) xmax(",f5.1,")")') dxmax
               write(iot,'( "p: ymin(0) ymax(",f5.1,")")') dymax

      end if

*-----------------------------------------------------------------------

  460       continue
  470       continue
  480       continue

*-----------------------------------------------------------------------
*        rz axis (matrix)
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 12 ) then

               if( ittwo(m) .eq. 1 ) then

                  dc2 = 'h2: '

               else if( ittwo(m) .eq. 2 ) then

                  dc2 = 'hd: '

               else if( ittwo(m) .eq. 3 ) then

                  dc2 = 'hc: '

               else if( ittwo(m) .eq. 6 ) then

                  dc2 = 'hd2:'

               else if( ittwo(m) .eq. 7 ) then

                  dc2 = 'hc2:'

               end if

*-----------------------------------------------------------------------

               inum = 0

                  if( nn .eq. 0 ) then

                     nc = 1

                  else

                     nc = nn

                  end if

*-----------------------------------------------------------------------
*              max and min
*-----------------------------------------------------------------------

                        fmax = 0.0
                        fmin = 1.e+33

                  do ic = 1, nc

                     if( nn .gt. 0 ) then

                        iz = nt(ic) / 1000

                        ia = nt(ic) - iz * 1000

                        if( ia .gt. 0 ) in = ia - iz

                     end if

                     do jz = 1, nz
                     do jr = 1, nr

                           sek = 0.0

                           if( itunt(m) .eq. 1 ) then

                              vm = 1.0

                           else if( itunt(m) .eq. 2 ) then

                              vm = vl(jr,jz)

                           end if

                        if( nn .eq. 0 ) then

                           do kn = 1, mn
                           do kz = 1, mz
                           do il = 0, mm ! frtati 2022/02/18

                              lz  = ikzz( kz, kn )
                              ln  = iknn( kz, kn )

                              vn  = vm * tr(jr,jz,lz,ln,il,1)
                              sek = sek + vn

                           end do
                           end do
                           end do

                        else if( ia .eq. 0 ) then

                           do kn = 1, mn
                           do il = 0, mm ! frtati 2022/02/18

                              lz  = ikzz( iz, kn )
                              ln  = iknn( iz, kn )

                              vn  = vm * tr(jr,jz,lz,ln,il,1)
                              sek = sek + vn

                           end do
                           end do

                        else

                              lz  = ikzz( iz, in )
                              ln  = iknn( iz, in )

                           do il = 0, mm ! frtati 2022/02/18
                              vn  = vm * tr(jr,jz,lz,ln,il,1)
                              sek = sek + vn
                           end do

                        end if

                              fmsv = sek / vm

                              if( fmsv .gt. fmax ) fmax = fmsv
                              if( fmsv .gt. 0.0 .and.
     &                            fmsv .lt. fmin ) fmin = fmsv

                     end do
                     end do

                  end do

*-----------------------------------------------------------------------

            do ic = 1, nc

               inum = inum + 1

*-----------------------------------------------------------------------

               if( nn .gt. 0 ) then

                  iz = nt(ic) / 1000

                  ia = nt(ic) - iz * 1000

                  if( ia .gt. 0 ) then

                     in = ia - iz

                     call chname(idum,ia,iz,chau)

                  end if

               end if

*-----------------------------------------------------------------------

                  write(iot,'(/"#",78("-"))')

               if( inum .eq. 1 ) then
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

               if( rtfac(m) .lt. 0.d0 )
     &           write(iot,
     &           '( "#restart maximum for normalization:",1p1e14.7)')
     &           facmax(m)

            if( nn .eq. 0 ) then

               write(iot,'("#   no. =",i3,3x,"Z = all")')
     &                     inum

            else if( ia .eq. 0 ) then

               write(iot,'("#   no. =",i3,3x,
     &         "Z = ",i3," : ",a3)')
     &                     inum, iz, elmnt(iz)

            else

               write(iot,'(/"#   no. =",i3,3x,a8)')
     &                     inum, chau

            end if

*-----------------------------------------------------------------------

            if( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) then

               if( nn .eq. 0 ) then

                  write(iot,'(/a1,"no. =",i3,
     &                        ",    Z  =  all", a1)')
     &                        cha, inum, cha

               else if( ia .eq. 0 ) then

                  write(iot,'(/a1,"no. =",i3,
     &                        ",    Z  =  ",i3, a1)')
     &                        cha, inum, iz, cha

               else

                  write(iot,'(/a1,"no. =",i3,
     &                        ",    ",a8, a1)')
     &                        cha, inum, chau, cha

               end if

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: z [cm]")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: r [cm]")')

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

*-----------------------------------------------------------------------

               if( inum .eq. 1 ) then

                  form  = ( rm(nr+1) - rm(1) )
     &                  / ( zm(nz+1) - zm(1) )
                  xfac  = 0.9
                  afac  = 0.8
                  izlog = 1
                  inocm = 1
                  inolg = 1

                  if( itanl(m) .gt. 0 )
     &            call anset(itang(m),itanl(m),
     &                       form,xfac,afac,fmin,fmax,izlog,inocm,inolg)

                  if( form .le. 1.0 ) then

                     scal = form**0.35
                     xfac = xfac / form**0.5
                     xorg = 0.0
                     yorg = min( 20.d0, ( 1.0 / form - 1.0 ) / 2.5 )

                  else

                     scal = 1.0 / form**0.41
                     xfac = xfac / form**0.5
                     xorg = min( 20.d0, ( form - 1.0 ) / 1.5)
                     yorg = 0.0

                  end if

               end if

                  write(iot,'( "set: c1[",f6.3,
     &                            "] c2[",f6.3,
     &                            "] c5[",f6.3,
     &                            "] c6[",f6.3,
     &                            "] c7[",f6.3,
     &                            "] c8[",f6.3,"]")')
     &                    form, xfac, afac, scal, xorg, yorg
                  write(iot,'( "p: h2fs form[c1] xfac[c2]",
     &                 " afac[c5] scal[c6] xorg[c7] yorg[c8] nosp")')

               if( ( ( ittwo(m) .ge. 2 .and. ittwo(m) .le. 3 ) .or.
     &               ( ittwo(m) .ge. 6 .and. ittwo(m) .le. 7 ) ) .and.
     &             fmin .gt. 0.0 .and. fmax .gt. cmin ) then

                if (ioe .eq. 1 ) then
                  write(iot,'( "set: c3[",1p1e15.6,"] c4[",
     &                  1p1e15.6,"]")') fmin, fmax
                else
                  write(iot,'( "set: c3[1.0e-4] c4[1.0]")')
                end if
                  write(iot,'( "p: cmin[c3] cmax[c4]")')
                  write(iot,'( "p: dmin(1e-31)")')

                  if( ittwo(m) .eq. 3 .or. ittwo(m) .eq. 7 )
     &                 write(iot,'( "p: zlog")')

               end if

*-----------------------------------------------------------------------

               zmin = zm(1)
               zmax = zm(nz+1)
               rmin = rm(1)
               rmax = rm(nr+1)

               write(iot,'( "p: xmin(",1p1g14.7,") xmax(",1p1g14.7,
     &                       ")")') zmin, zmax

               write(iot,'( "p: ymin(",1p1g14.7,") ymax(",1p1g14.7,
     &                       ")")') rmin, rmax

               if( itanl(m) .gt. 0 .and. ioe .eq. 1 ) then

                  write(iot,'( "p: ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               else if( itanl(m) .gt. 0 .and. ioe .ne. 1 ) then

                   call terrang(iot,m)

               end if

               if( itsans(m) .gt. 0 ) then

                  call write_sangel(iot,m,0)

               end if

            end if

*-----------------------------------------------------------------------

               do jz = 1, nz
               do jr = 1, nr

                     sek = 0.0

                     if( itunt(m) .eq. 1 ) then

                        vm = 1.0

                     else if( itunt(m) .eq. 2 ) then

                        vm = vl(jr,jz)

                     end if

*-----------------------------------------------------------------------

                  if( nn .eq. 0 ) then

                     do kn = 1, mn
                     do kz = 1, mz
                     do il = 0, mm ! frtati 2022/02/18

                        lz  = ikzz( kz, kn )
                        ln  = iknn( kz, kn )

                        vn  = vm * tr(jr,jz,lz,ln,il,ioe)
                        sek = sek + vn

                     end do
                     end do
                     end do

                  else if( ia .eq. 0 ) then

                     do kn = 1, mn
                     do il = 0, mm ! frtati 2022/02/18

                        lz  = ikzz( iz, kn )
                        ln  = iknn( iz, kn )

                        vn  = vm * tr(jr,jz,lz,ln,il,ioe)
                        sek = sek + vn

                     end do
                     end do

                  else

                        lz  = ikzz( iz, in )
                        ln  = iknn( iz, in )

                     do il = 0, mm ! frtati 2022/02/18
                        vn  = vm * tr(jr,jz,lz,ln,il,ioe)
                        sek = sek + vn
                     end do

                  end if

                        fm(jz,jr) = sek / vm

               end do
               end do

*-----------------------------------------------------------------------

               write(iot,'( "#  nr = ",i3,"   nz = ",i3)')
     &                       nr, nz

            if( ittwo(m) .ne. 4 ) then

               write(iot,'( "# ( ( data(z,r), z = 1, nz ),",
     &                      " r = nr, 1, -1 )")')

            end if

*-----------------------------------------------------------------------

            if( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) then

               write(iot,'(/a4," y = ",1p1g14.7," to ",
     &         1p1g14.7," by ",1p1g14.7," ; x = ",
     &         1p1g14.7," to ",1p1g14.7," by ",
     &         1p1g14.7," ;")') dc2,
     &         rm(nr) + rtrdl(m)/2.0, rm(1) + rtrdl(m)/2.0, rtrdl(m),
     &         zm(1) + rtzdl(m)/2.0, zm(nz) + rtzdl(m)/2.0, rtzdl(m)

               write(iot,'(1p10e11.3)')
     &         ( ( fm(jz,jr), jz = 1, nz ), jr = nr, 1, -1 )

*-----------------------------------------------------------------------

            else if( ittwo(m) .eq. 4 ) then

               write(iot,'(/"# r          z        ",
     &                      "  number")')

               do jz = 1, nz
               do jr = 1, nr

                  write(iot,'(1p10e11.3)')
     &               rm(jr)  + rtrdl(m)/2.0,
     &               zm(jz)  + rtzdl(m)/2.0,
     &               fm(jz,jr)

               end do
               end do

*-----------------------------------------------------------------------

            else if( ittwo(m) .eq. 5 ) then

               write(iot,'("#   r = ",1p1g14.7," to ",
     &         1p1g14.7," by ",1p1g14.7/
     &                     "#   z = ",1p1g14.7," to ",
     &         1p1g14.7," by ",1p1g14.7)')
     &         rm(1) + rtrdl(m)/2.0, rm(nr) + rtrdl(m)/2.0, rtrdl(m),
     &         zm(1) + rtzdl(m)/2.0, zm(nz) + rtzdl(m)/2.0, rtzdl(m)

                  write(iot,'(/4x,a3,4x,1p1000e11.3)')
     &            'r/z',( zm(jz) + rtzdl(m)/2.0, jz = 1, nz )

               do jr = nr, 1, -1

                  write(iot,'(1p1000e11.3)')
     &            rm(jr) + rtrdl(m)/2.0,
     &            ( fm(jz,jr), jz = 1, nz )

               end do

            end if

*-----------------------------------------------------------------------

      if( inolg .eq. 1 .and. ( ittwo(m) .eq. 3 .or. ittwo(m) .eq. 7 )
     &    .and. fmin .gt. 0.0 .and. fmax .gt. fmin .and.
     &    igsh .eq. 0 ) then

      write(iot,'(/"#",78("-"))')
      write(iot,'(
     &"z: xorg(1.03)"/
     &"p: xfac[c2*0.05] form[c1/0.05] notc noxt noyt dmin(1)"/
     &"p: ymin(0) ymax(1) xmin(0) xmax(1) cmin(1) cmax(100) nosp"/
     &"hc: y= 0.005 to 0.995 by 0.01 ; x= 0.5 to 0.5 by 1 ;"/
     &" 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20"/
     &"21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40"/
     &"41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60"/
     &"61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80"/
     &"81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100"/
     &"z: xorg(0.0)"/
     &"p: xfac[c2*0.05] form[c1/0.05] ",
     &"nosp afac[c5*0.625] ymin[c3] ymax[c4]"/
     &"p: noxt noxn ytxt(-1) ynum(-1) xltd(-1) xstd(-1) itic(1) ylog"
     &)')

      if ( ioe .eq. 2 ) then

          write(iot,'("y: Relative Error")')

      else if( itazl(m) .eq. 0 ) then

         write(iot,'("y: Number ",a15)') hsunit(itunt(m))

      else

         write(iot,'("y: ",200a1)') (itazt(m)(i:i),i=1,itazl(m))

      end if

      if( izlog .eq. 0 ) write(iot,'( "p: ylin")')

      write(iot,'( "#",78("-"))')

      end if

*-----------------------------------------------------------------------

            if( inocm .eq. 1 ) then

               if( nn .eq. 0 ) then

                  write(iot,'(/"wt: s(0.7)",/
     &                        a1,"vspace{-3}"/
     &                        "      Z &=&   all"/
     &                        "e:")')
     &                        yen

               else if( ia .eq. 0 ) then

                  write(iot,'(/"wt: s(0.7)",/
     &                        a1,"vspace{-3}"/
     &                        "      Z &=& ",i3,"   :  ",a3/
     &                        "e:")')
     &                        yen,
     &                        iz, elmnt(iz)

               else

                  write(iot,'(/"wt: s(0.7)",/
     &                        a1,"vspace{-3}"/
     &                        "  nucl. &=&   ",a8/
     &                        "e:")')
     &                        yen,
     &                        chau

               end if

            end if

*-----------------------------------------------------------------------

      if( inolg .eq. 1 .and. ( ittwo(m) .eq. 3 .or. ittwo(m) .eq. 7 )
     &    .and. cmin .gt. 0.0 .and. cmax .gt. cmin .and.
     &    igsh .eq. 0 ) then

         write(iot,'(/"#",78("-"))')
         write(iot,'("z: xorg[-1.03/0.05]"/
     &               "p: form[c1] xfac[c2] afac[c5] nosp notf")')

               write(iot,'( "p: xmin(",1p1g14.7,") xmax(",1p1g14.7,
     &                       ")")') zmin, zmax

               write(iot,'( "p: ymin(",1p1g14.7,") ymax(",1p1g14.7,
     &                       ")")') rmin, rmax

      end if

*-----------------------------------------------------------------------

            end do

*-----------------------------------------------------------------------

         end if

            call prestart(m,iot) !OBINATA(2012.8.20)

            close(iot)

         if( iteps(m) .ne. 0 ) then

            open(iot, file = fname, status = 'unknown' )
            call a_angel(idasa,fname)

         end if

  900    continue

      end do

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine tyilxyz(ncol,m,mz,mn,mm,nl,lt, ! frtati 2022/02/18 added mm
     &                   nx,ny,nz,nm,xm,ym,zm,
     &                   mt,ikzz,iknn,tr,
     &                   trEVENT,
     &                   itrmax,itrmin)
*                                                                      *
*       nuclear yield (or production) tally in xyz mesh                *
*       last modified by K.Niita on 2014/11/11                         *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*      ncol  ..... reaction type                                       *
*             11 : termination by energy cut-off                       *
*             13 : nuclear reaction (n,x)                              *
*             14 : nuclear reaction (n,n'x)                            *
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
*        kcoll : =  0, normal                                          *
*                =  1, high energy fission                             *
*                =  2, high energy absorption                          *
*                =  3, low energy n elastic                            *
*                =  4, low energy n non-elastic                        *
*                =  5, low energy n fission                            *
*                =  6, low energy n absorption                         *
*                                                                      *
************************************************************************
      use QMD_COOD2_MOD, only : nqmdm
      use MMBANKMOD !FURUTA
      use partmod, only: itmxpt,itpan,itpat,jtpat,iznmmx,iznmturn ! frtati 2022/05/02
      use moddas_region ! S.H. (2022.11.9)
*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'param00.inc'
      include 'param01.inc'

*-----------------------------------------------------------------------

      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
!$OMP THREADPRIVATE(/icomon/)
      common /wtsave/ oldwt
!$OMP THREADPRIVATE(/wtsave/)
      common /qparm/  ielas,icasc,iqstep,lvlopt,igamma
      common /bparm/  andt,jevap,npidk
      common /geosig/ geosig(250)
      common /cparm/  maxbch,maxcas
      common /mathzn/ mathz, mathn, jcoll, kcoll
!$OMP THREADPRIVATE(/mathzn/)

      common /trstar/ itrstar
!$OMP THREADPRIVATE(/trstar/)

      common /kmat1d/ idmn(0:kvlmax), idnm(kvmmax)
      common /clionprd/  lionprd
      common /eparm/  esmax, esmin, emin(20)

      common /regdc/  idrg(kvlmax), idgr(kvmmax)
      common /tlgeom/ iblz1,iblz2
!$OMP THREADPRIVATE(/tlgeom/)
*-----------------------------------------------------------------------

      common /clustf/ nclst, iclust(nnn)
!$OMP THREADPRIVATE(/clustf/)
      common /clustt/ nclsts, iclusts(nnn)
!$OMP THREADPRIVATE(/clustt/)
      common /clustw/ jclusts(0:8,nnn),  qclusts(0:12,nnn)
!$OMP THREADPRIVATE(/clustw/)
      common /clustl/ rumpal(0:20), numpal(0:20)
!$OMP THREADPRIVATE(/clustl/)
      common /cntcls/ jcount(3,nnn)
!$OMP THREADPRIVATE(/cntcls/)

*-----------------------------------------------------------------------

      dimension numsav(0:20), rumsav(0:20)

*-----------------------------------------------------------------------

      common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                itout(itlmax), ittwo(itlmax)
      common /tall12/ itman(itlmax), itmat(itlmax), itmct(itlmax),
     &                itnun(itlmax), itnuc(itlmax),
     &                itndz(itlmax), itndn(itlmax),
     &                itnkz(itlmax), itnkn(itlmax)
      common /tall22/ itmcn(itlmax), itmtn(itlmax), itmtt(itlmax)
      common /tall30/ itnda(itlmax)
      common /tall36/ itdpo(itlmax)
      common /tall37/ itcnt(9,itlmax)
      common /tall46/ itmtr(itlmax,4), rtmtr(itlmax,13)
      common /tall52/ itprd(itlmax)

      common /tall82/ itcnth(9,itlmax)
      common /tall83/ itnzn(itlmax), itndm(itlmax)

*-----------------------------------------------------------------------

      dimension   lt(nl)
      dimension   xm(nx+1)
      dimension   ym(ny+1)
      dimension   zm(nz+1)
      dimension   mt(nm)
      dimension   ikzz(maxpt,maxnt), iknn(maxpt,maxnt)
      dimension   tr(nx,ny,nz,mz,mn,0:mm,2)
      dimension   trEVENT(nx,ny,nz,mz,mn,0:mm)     !OBINATA(2012.8.20): as Ct
      real(8),allocatable,save:: tr0(:,:,:,:,:,:) !OBINATA(2012.8.20): as C
      dimension   itrmax(5),itrmin(5)

      dimension mnz(mxprodxs),mna(mxprodxs),xxn(mxprodxs),ildd(mxprodxs)

      dimension   ips(itmxpt) ! frtati 2021/10/05 6 -> itmxpt
      dimension   ncntt(3)


*-----------------------------------------------------------------------
      common /cntreg/ ncntc(3), ncreg(3), incrc(3), incrt(3), kcont(3)
      common /cntpat/ icpan(3), icpat(3,20,2)
      dimension     idas(1)
      equivalence ( das, idas )
      common /pnsave/ egs, uus, vvs, wws, wts, tms, nms, nct(3)
!$OMP THREADPRIVATE(/pnsave/)

*-----------------------------------------------------------------------
      common /cntmx/ ncntmx(3), nctmxsr(3)
!$OMP THREADPRIVATE(/cntmx/)

      common /mpi00/ npe, me
      real*8,allocatable :: tryld(:,:)

*-----------------------------------------------------------------------

      common /jcomon/ nabov,nobch,nocas,nomax
!$OMP THREADPRIVATE(/jcomon/)
      common /stat / istdev, irestart, ireschk

      integer, allocatable, save :: iclusav(:), jclusav(:,:)
      double precision, allocatable, save :: qclusav(:,:)
!$OMP THREADPRIVATE(iclusav, jclusav, qclusav)

       if( .not. allocated(iclusav) ) then ! initial allocation
        allocate(iclusav(nqmdm),jclusav(0:8,nqmdm),qclusav(0:12,nqmdm))
        iclusav = 0
        jclusav = 0
        qclusav = 0.d0
       elseif( ubound(iclusav,1) .lt. nqmdm ) then ! extend array
        deallocate(iclusav,jclusav,qclusav)
        allocate(iclusav(nqmdm),jclusav(0:8,nqmdm),qclusav(0:12,nqmdm))
        iclusav = 0
        jclusav = 0
        qclusav = 0.d0
       endif

*-----------------------------------------------------------------------

      if (istdev .eq. 2) then
        call readitrminmax5(itrmin,itrmax,(/nx,ny,nz,mz,mn/),
     &                      mix,miy,miz,mimz,mimn,
     &                      mxx,mxy,mxz,mxmz,mxmn)
      endif

*-----------------------------------------------------------------------
* check of history counter
*-----------------------------------------------------------------------
      ihistcount = 0 ! history counter index, T.Sato 2022/12/20 for speed up
         if ( ncol.eq.0 .or. ncol.eq.4 ) then
           do i = 1, 3
             if( itcnth(i,m) .eq. 1 ) then
               if( ncntmx(i) .lt. itcnth(i*2+2,m) .or.
     &             ncntmx(i) .gt. itcnth(i*2+3,m) ) then
                 ihistcount = 1
               end if
             end if
           end do
         end if

*-----------------------------------------------------------------------
*        source particle or end of batch ( in case of istdev = 2 )
*
* OBINATA(2012.8.20): change tr(,,,,3) to trEVENT(,,,,)
*-----------------------------------------------------------------------

      if (( ncol .eq. 0 .or. ncol .eq. 4 )
     &                            .and. istdev .eq. 2) then
        if ((nocas.gt.1.or.ncol.eq.0) .and. ihistcount.ne.1 ) then

            if( itnzn(m).eq.0 ) then ! frtati 2022/02/18
             do imn = mimn,mxmn
             do imz = mimz,mxmz
             do iz = miz,mxz
             do iy = miy,mxy
             do ix = mix,mxx
               tr(ix,iy,iz,imz,imn,:,1) = tr(ix,iy,iz,imz,imn,:,1)
     &                                + trEVENT(ix,iy,iz,imz,imn,:)
               tr(ix,iy,iz,imz,imn,:,2) = tr(ix,iy,iz,imz,imn,:,2)
     &                                + trEVENT(ix,iy,iz,imz,imn,:)**2
             enddo
             enddo
             enddo
             enddo
             enddo
            else
             do imz = 1,iznmmx(m)
             do iz = miz,mxz
             do iy = miy,mxy
             do ix = mix,mxx
               tr(ix,iy,iz,imz,1,0,1) = tr(ix,iy,iz,imz,1,0,1)
     &                                + trEVENT(ix,iy,iz,imz,1,0)
               tr(ix,iy,iz,imz,1,0,2) = tr(ix,iy,iz,imz,1,0,2)
     &                                + trEVENT(ix,iy,iz,imz,1,0)**2
             enddo
             enddo
             enddo
             enddo
            end if

       end if

          if( itnzn(m).eq.0 ) then ! frtati 2022/02/18
           do imn = mimn,mxmn
           do imz = mimz,mxmz
           do iz = miz,mxz
           do iy = miy,mxy
           do ix = mix,mxx
             trEVENT(ix,iy,iz,imz,imn,:) = 0
           enddo
           enddo
           enddo
           enddo
           enddo
          else
           do imz = 1,iznmmx(m)
           do iz = miz,mxz
           do iy = miy,mxy
           do ix = mix,mxx
             trEVENT(ix,iy,iz,imz,1,0) = 0
           enddo
           enddo
           enddo
           enddo
          end if

           call resetitrminmax(itrmin,itrmax,5,(/nx,ny,nz,mz,mn/))

      end if

*-----------------------------------------------------------------------
*        end of batch ( in case of istdev = 1 )
*
* OBINATA(2012.8.20): modificate for thread parallel
*-----------------------------------------------------------------------

         if ( ncol .eq. 0 .and. istdev .eq. 1) then
!$OMP MASTER
             allocate( tr0(nx,ny,nz,mz,mn,0:mm) ) ! frtati 2022/02/18
             tr0(:,:,:,:,:,:) = 0.d0
!$OMP END MASTER
!$OMP BARRIER
!$OMP CRITICAL (tyilxyz_crit_ist1)
             tr0(:,:,:,:,:,:) = tr0(:,:,:,:,:,:) + trEVENT(:,:,:,:,:,:)
!$OMP END CRITICAL (tyilxyz_crit_ist1)
!$OMP BARRIER
!$OMP MASTER
       tr(:,:,:,:,:,:,1) = tr(:,:,:,:,:,:,1) + tr0(:,:,:,:,:,:) / maxcas
       tr(:,:,:,:,:,:,2) = tr(:,:,:,:,:,:,2)
     &                   + ( tr0(:,:,:,:,:,:) / maxcas ) ** 2
             deallocate( tr0 )
!$OMP END MASTER

           trEVENT(:,:,:,:,:,:) = 0

         end if

*-----------------------------------------------------------------------
*        rearrangement of iznm for MPI run ! frtati 2022/05/02
*-----------------------------------------------------------------------

         if( npe.gt.1 .and. itnzn(m).ne.0 .and. ncol.eq.0 ) then
           call paraiznm(m,1)
           if( me.gt.1 ) then
             allocate( tryld(iznmmx(m),2) )
             do iz = 1,nz
             do iy = 1,ny
             do ix = 1,nx
               do imz = 1, iznmmx(m)
                 tryld(imz,:) = tr(ix,iy,iz,imz,1,0,:)
               end do
               do imz = 1, iznmmx(m)
                 tr(ix,iy,iz,imz,1,0,:) = tryld(iznmturn(imz,m),:)
               end do
             end do
             end do
             end do
             deallocate( tryld )
           end if
         end if

*-----------------------------------------------------------------------
*        check of ncol and nclsts ( outgoing particles )
*-----------------------------------------------------------------------

         if( ncol .eq. 11 ) then

            if( itprd(m) .ne. 1 ) return

               iccol = 0
               npart = 1

               if( ityp .lt. 15 .or. ityp .gt. 19 ) return

               jz    = jtyp
               jn    = ktyp - ktyp / 1000000 * 1000000 - jz

               ipart = ityp
               tlw   = oldwt
               il = 0  ! T.Sato 2016/3/29, isomer is not considered for output=cutoff

               ncntt(1) = ncnt(ibknct+1,no,ipomp+1)
               ncntt(2) = ncnt(ibknct+2,no,ipomp+1)
               ncntt(3) = ncnt(ibknct+3,no,ipomp+1)

         else if( ncol .eq. 13 .or. ncol .eq. 14 ) then

            if ( itnda(m) .ge. 2 .and. jcoll .eq. 9 ) then
               if( nclsts .lt. 0 ) return
            else
               if( nclsts .le. 0 ) return
            end if

            if( ( itnda(m) .eq. 2 .and.
     &          ( jcoll .eq. 6 .or. jcoll .eq. 9 ) ) .or.
     &          ( itnda(m) .eq. 3 .and.
     &          ( jcoll .eq. 6 .or. jcoll .eq. 9 .or.
     &            jcoll .eq. 4 .or. jcoll .eq. 5 .or.
     &            jcoll .eq. 15 .or. jcoll .eq. 10 ) ) ) then

            else
               if( itdpo(m) .eq. 0 .and.
     &           ( jcoll .eq. 3 .or. kcoll .eq. 3 ) ) return
               if( itdpo(m) .eq. -1 .and. itrstar .eq. 1 ) return
            end if

               iccol = 1

               npart = nclsts

               ncntt(1) = jcount(1,1)
               ncntt(2) = jcount(2,1)
               ncntt(3) = jcount(3,1)

         else

            return

         end if

*-----------------------------------------------------------------------
*        check of counter
*-----------------------------------------------------------------------

            do i = 1, 3

               if( itcnt(i,m) .ne. 0 ) then

                  if( ncntt(i) .lt. itcnt(i*2+2,m) .or.
     &                ncntt(i) .gt. itcnt(i*2+3,m) ) return

               end if

            end do

*-----------------------------------------------------------------------
*        check of mat
*-----------------------------------------------------------------------

            if( nl .gt. 0 ) then

                  do i = 1, nl

                     if( itmcn(m) .gt. 0 .and.
     &                   idmn(mat) .eq. lt(i) ) goto 502
                     if( itmcn(m) .lt. 0 .and.
     &                   idmn(mat) .eq. lt(i) ) return

                  end do

                     if( itmcn(m) .gt. 0 ) return

            end if

  502          continue

*-----------------------------------------------------------------------

         if( iccol .eq. 1 ) then

*-----------------------------------------------------------------------
*           check of particles
*-----------------------------------------------------------------------

            call pcheck(m,itpan(m),ityp,ktyp,jtyp,ipn,ips)

               if( ipn .eq. 0 ) return

*-----------------------------------------------------------------------
*           mother and energy
*-----------------------------------------------------------------------

               ata   = dble( mathz + mathn )
               atz   = dble( mathz )
               mmas  = nint( ata )
               mchg  = nint( atz )
               eein  = ec(ibkec+no,ipomp+1)

*-----------------------------------------------------------------------
*            for specific mother nuclei
*-----------------------------------------------------------------------

            if( nm .gt. 0 ) then

               do i = 1, nm

                  iz = mt(i) / 1000
                  ia = mt(i) - iz * 1000

                  if( ( ia .gt. 0 .and.
     &                  mmas .eq. ia .and. mchg .eq. iz ) .or.
     &                ( ia .eq. 0 .and. mchg .eq. iz ) ) then

                     if( itmct(m) .gt. 0 ) goto 30
                     if( itmct(m) .lt. 0 ) return

                  end if

               end do

                  if( itmct(m) .gt. 0 ) return

            end if

   30       continue

         end if


*-----------------------------------------------------------------------
*        transform positions
cKN      yield position is assumed to be the reaction point.
*-----------------------------------------------------------------------

            xpart = xc(ibkxc+no,ipomp+1)
            ypart = yc(ibkyc+no,ipomp+1)
            zpart = zc(ibkzc+no,ipomp+1)

            call trnsxx(xpart,ypart,zpart,
     &                  xxc,yyc,zzc,itmtr(m,4))

*-----------------------------------------------------------------------
*           check z mesh and r mesh region
*-----------------------------------------------------------------------

               xpp = xxc
               ypp = yyc
               zpp = zzc

            if(  xpp .lt. xm(1) .or. xpp  .ge. xm(nx+1) ) return
            if(  ypp .lt. ym(1) .or. ypp  .ge. ym(ny+1) ) return
            if(  zpp .lt. zm(1) .or. zpp  .ge. zm(nz+1) ) return

               do i = 2, nx + 1

                  if( xpp .lt. xm(i) ) goto 36

               end do

   36       ixc = i - 1

               do i = 2, ny + 1

                  if( ypp .lt. ym(i) ) goto 37

               end do

   37       iyc = i - 1

               do i = 2, nz + 1

                  if( zpp .lt. zm(i) ) goto 38

               end do

   38       izc = i - 1

*-----------------------------------------------------------------------
*        normal weight
*-----------------------------------------------------------------------

               idoy = 1
               idev = 1
               wyld = oldwt

*-----------------------------------------------------------------------
*        for special : repeated nuclear reactions
*-----------------------------------------------------------------------

         if( iccol .eq. 1 .and. itspc(m) .gt. 0 .and.
     &       mathz .gt. 2 .and.
     &     ( jcoll .eq. 4 .or. jcoll .eq. 5 ) ) then

               idoy = itspc(m)
               idev = 10
               wyld = oldwt / dble( idoy * idev )

*-----------------------------------------------------------------------
*           save normal output and reaction parameters
*-----------------------------------------------------------------------

                  nclssav = nclsts

            do i = 1, nclsts

                  iclusav(i) = iclusts(i)

               do k = 0, 8

                  jclusav(k,i) = jclusts(k,i)

               end do

               do k = 0, 12

                  qclusav(k,i) = qclusts(k,i)

               end do

            end do

               do k = 0, 20

                  numsav(k) = numpal(k)
                  rumsav(k) = rumpal(k)

               end do

                  iprj  = ityp
                  kprj  = ktyp

                  if( ityp .lt. 15 ) then

                     bmax = sqrt( geosig(mmas) * 100.0 / 3.1415926 )

                  else if( ityp .ge. 15 ) then

                     ap = dble( ktyp - ktyp / 1000000 * 1000000 )
                     zp = dble( ktyp / 1000000 )

                     call sighi(ap,zp,eein,ata,zta,signe,sigel,bmax)

                  end if

         end if

*-----------------------------------------------------------------------
*     repeated do loop
*-----------------------------------------------------------------------

         do ireac = 1, idoy

*-----------------------------------------------------------------------
*     repeat calculation of nuclear reactions
*-----------------------------------------------------------------------

            if( idoy .gt. 1 .and. ireac .gt. 1 ) then

                  ipim = 0

   22          continue

                  call ncasc(1,iprj,kprj,eein,mmas,mchg,bmax)

               if( nclst .lt. 0 ) then

                     ipim = ipim + 1

                     if( ipim .le. 20 ) goto 22

               end if

            end if

*-----------------------------------------------------------------------
*     repeat calculation of evaporation
*-----------------------------------------------------------------------

            do ievap = 1, idev

               if( idev .gt. 1 ) then

                  call nevap(1)

                  do i = 1, nclsts

                     qclusts(8,i) = wyld * qclusts(8,i)

                  end do

                  npart = nclsts

               end if

*-----------------------------------------------------------------------
*        booking after nuclear reactions
*-----------------------------------------------------------------------

                  istat = 1
                  mnx   = 0

*-----------------------------------------------------------------------
*           NDATA=2,3 : Replace with yield data (activation cross section
*-----------------------------------------------------------------------

               if( iccol .eq. 1 .and.
     &           ( itnda(m) .eq. 2 .and.
     &           ( jcoll .eq. 6 .or. jcoll .eq. 9 ) ) .or.
     &           ( itnda(m) .eq. 3 .and.
     &           ( jcoll .eq. 6 .or. jcoll .eq. 9 .or.
     &             jcoll .eq. 4 .or. jcoll .eq. 5 .or.
     &             jcoll .eq. 15 .or. jcoll .eq. 10 ) ) ) then

                  mkk = mat
                  icl = idgr(iblz1)
                  jcl = jcoll

                  call prodxs2(ktyp,mchg,mmas,eein,
     &                 istat,mnx,mnz,mna,xxn,xxs,
     &                 ildd,mkk,icl,jcl,m)

                  iexclight = 0 ! Exclude light ions
                  if( itprd(m) .eq. 1 ) then ! output = cutoff
                     iexclight = 1
                     if ( jcoll.eq.6 .and.
     &                    ktyp.eq.2112 .and. eein.lt.20.0 )
     &                    iexclight = 0
                  end if

                  if( istat .eq. 0 .and. mnx .gt. 0 ) then

                   do nnx = 1, mnx

                      jz = mnz(nnx)
                      jn = mna(nnx) - mnz(nnx)
                      il = ildd(nnx)

                    if ( iexclight.ne.1 .or. jz.gt.2 ) then

                     if( jz .gt. 0 .and. jz .le. maxpt .and.
     &                      jn .gt. 0 .and. jn .le. maxnt ) then

                        kz = ikzz(jz,jn)
                        kn = iknn(jz,jn)

                      if( itnzn(m).eq.0 ) then ! frtati 2022/02/18
                       if( kz .le. mz .and. kn .le. mn ) then
                          trEVENT(ixc,iyc,izc,kz,kn,il) =
     &                   trEVENT(ixc,iyc,izc,kz,kn,il) + wyld * xxn(nnx)

                        if (istdev .eq. 2) then
                           call setitrminmax(itrmin,itrmax,1,5,
     &                          (/ixc,iyc,izc,kz,kn/))
                        endif

                       end if

                      else
                         trEVENT(ixc,iyc,izc,
     &                        igetiznm(kz,kn,il,m),1,0) =
     &                        trEVENT(ixc,iyc,izc,
     &                       igetiznm(kz,kn,il,m),1,0) + wyld * xxn(nnx)
                       if (istdev .eq. 2) then
                          call setitrminmax(itrmin,itrmax,1,5,
     &                         (/ixc,iyc,izc,kz,kn/))
                       endif
                      end if

                     end if

                    end if

                   end do

                  end if

               end if

*-----------------------------------------------------------------------
*           nuclear data for 4He, 14N, 16O
*-----------------------------------------------------------------------

               if( iccol .eq. 1 .and.
     &             itnda(m) .eq. 1 .and. jcoll .eq. 4 ) then

                     call prodxs(ktyp,mchg,mmas,eein,
     &                           istat,mnx,mnz,mna,xxn,xxs,ildd)

                  if( istat .eq. 0 .and. mnx .gt. 0 ) then

                     do nnx = 1, mnx

                        jz = mnz(nnx)
                        jn = mna(nnx) - mnz(nnx)
                        il = ildd(nnx)

                        if( jz .gt. 0 .and. jz .le. maxpt .and.
     &                      jn .gt. 0 .and. jn .le. maxnt ) then

                           kz = ikzz(jz,jn)
                           kn = iknn(jz,jn)

                        if( itnzn(m).eq.0 ) then ! frtati 2022/02/18

                        if( kz .le. mz .and. kn .le. mn ) then

                         trEVENT(ixc,iyc,izc,kz,kn,il) =
     &                   trEVENT(ixc,iyc,izc,kz,kn,il) + wyld * xxn(nnx)

                           if (istdev .eq. 2) then
                             call setitrminmax(itrmin,itrmax,1,5,
     &                                         (/ixc,iyc,izc,kz,kn/))
                           endif

                        end if

                        else
                           trEVENT(ixc,iyc,izc,
     &                     igetiznm(kz,kn,il,m),1,0) =
     &                     trEVENT(ixc,iyc,izc,
     &                     igetiznm(kz,kn,il,m),1,0) + wyld * xxn(nnx)
                           if (istdev .eq. 2) then
                             call setitrminmax(itrmin,itrmax,1,5,
     &                                         (/ixc,iyc,izc,kz,kn/))
                           endif
                        end if

                        end if

                     end do

                  end if

               end if


*-----------------------------------------------------------------------
*           normal case
*-----------------------------------------------------------------------

            if( npart .gt. 0 ) then

               do 500 j = 1, npart

*-----------------------------------------------------------------------

                  if( iccol .eq. 1 ) then

                     ipart = jclusts(3,j)
                     tlw   = qclusts(8,j)
                     jz    = jclusts(1,j)
                     jn    = jclusts(2,j)
                     il    = jclusts(8,j)

                  end if

*-----------------------------------------------------------------------
* When ndata=2 or 3, counter values change by ndata in [counter]
                 if ( (ipart.ge.3 .and. ipart.le.13) .or.
     &                 ipart.ge.19 ) then

                  if( ( itnda(m) .eq. 2 .or. itnda(m) .eq. 3 ) .and.
     &                 ( istat .eq. 0 .or. istat .eq. 3 ) ) then

                   if( ncntc(1) .eq. 1 .or. ncntc(2) .eq. 1 .or.
     &                    ncntc(3) .eq. 1 ) then
                    do k = 1, 3
                     ndatcount = 0
                     if( ncntc(k) .eq. 1 ) then
                      knn  = ityp
                      kcg  = jtyp
                      kkf  = ktyp
                      call pcchck(k,knn,kkf,kcg,icpan,icpat,icc)

                      if( icc .eq. 1 ) then
                       kdsm = kcont(k)
                       ldsm = 22
                       ndatcount = idas_kcont(kdsm+ldsm)
                      end if

                      if ( ndatcount .ne. 0 ) then

                       if ( ncol .eq. 13 ) then
                        ncnta(ibknct+k,j,ipomp+1)
     &                         = ncnta(ibknct+k,j,ipomp+1)
     &                         + ndatcount

                       else if ( ncol .eq. 14 ) then
                        if ( j .eq. 1 ) then
                         nct(k) = nct(k) + ndatcount
                        else if (j .gt. 1 ) then
                         ncnta(ibknct+k,j-1,ipomp+1)
     &                          = ncnta(ibknct+k,j-1,ipomp+1)
     &                          + ndatcount
                        end if

                       end if

                      end if

                     end if
                    end do
                   end if

                  end if

                 end if

*-----------------------------------------------------------------------

                  if( ipart .lt. 15 .or. ipart .gt. 19 ) goto 500

                  if( iccol .eq. 1 .and. itprd(m) .eq. 1 ) then

                     emint = emin(ipart) * dble( jz + jn )

                     if( qclusts(7,j) .gt. emint ) goto 500

                  end if

*-----------------------------------------------------------------------

                     if( ( itnda(m) .eq. 2 .or. itnda(m) .eq. 3 ) .and.
     &                 istat .eq. 0 .and.
     &                 ( iexclight.ne.1 .or. jz.gt.2 ) ) goto 500

                     if( ( itnda(m) .eq. 2 .or. itnda(m) .eq. 3 ) .and.
     &                    istat .eq. 3 .and. lionprd .ne. 1) goto 500

                     if( itnda(m) .eq. 1 .and.
     &                   istat .eq. 0 .and. mnx .gt. 0 ) then

                        do nnx = 1, mnx

                           if( jz .eq. mnz(nnx) .and.
     &                         jn .eq. mna(nnx) - mnz(nnx) ) goto 500

                        end do

                     end if

                        if( jz .gt. 0 .and. jz .le. maxpt .and.
     &                      jn .gt. 0 .and. jn .le. maxnt ) then

                           kz = ikzz(jz,jn)
                           kn = iknn(jz,jn)

                        if( itnzn(m).eq.0 ) then ! frtati 2022/02/18

                        if( kz .le. mz .and. kn .le. mn ) then
                           trEVENT(ixc,iyc,izc,kz,kn,il) =
     &                     trEVENT(ixc,iyc,izc,kz,kn,il) + tlw

                           if (istdev .eq. 2) then
                             call setitrminmax(itrmin,itrmax,1,5,
     &                                         (/ixc,iyc,izc,kz,kn/))
                           endif

                        end if

                        else
                           trEVENT(ixc,iyc,izc,
     &                     igetiznm(kz,kn,il,m),1,0) =
     &                     trEVENT(ixc,iyc,izc,
     &                     igetiznm(kz,kn,il,m),1,0) + tlw
                           if (istdev .eq. 2) then
                             call setitrminmax(itrmin,itrmax,1,5,
     &                                         (/ixc,iyc,izc,kz,kn/))
                           endif
                        end if

                        end if

  500          continue

            end if

*-----------------------------------------------------------------------

            end do

         end do

*-----------------------------------------------------------------------
*        store the decrease of mother by negative weight
*-----------------------------------------------------------------------

        if( ncol .ne. 11 ) then ! S.H. 2021.12.22
         if( itdpo(m) .eq. -1 .and. itrstar .eq. 0 ) then

            if( mathz .gt. 0 .and. mathz .le. maxpt .and.
     &          mathn .gt. 0 .and. mathn .le. maxnt ) then

              kz = ikzz(mathz,mathn)
              kn = iknn(mathz,mathn)

              if( itnzn(m).eq.0 ) then ! frtati 2022/02/18
               if( kz .le. mz .and. kn .le. mn ) then ! T.Sato 2024/05/02
                trEVENT(ixc,iyc,izc,kz,kn,0) =
     &                     trEVENT(ixc,iyc,izc,kz,kn,0) - oldwt

                if(istdev .eq. 2) then
                 call setitrminmax(itrmin,itrmax,1,5,
     &                                (/ixc,iyc,izc,kz,kn/))
                endif
               endif
              else
               trEVENT(ixc,iyc,izc,
     &         igetiznm(kz,kn,0,m),1,0) =
     &         trEVENT(ixc,iyc,izc
     &         ,igetiznm(kz,kn,0,m),1,0) - oldwt
                  if (istdev .eq. 2) then
                    call setitrminmax(itrmin,itrmax,1,5,
     &                                (/ixc,iyc,izc,kz,kn/))
                  endif
              end if

            end if

         end if
        end if

*-----------------------------------------------------------------------
*        restore normal output
*-----------------------------------------------------------------------

         if( idoy .gt. 1 ) then

                  nclsts = nclssav

            do i = 1, nclsts

                  iclusts(i) = iclusav(i)
                  iclusav(i) = 0

               do k = 0, 8

                  jclusts(k,i) = jclusav(k,i)
                  jclusav(k,i) = 0

               end do

               do k = 0, 12

                  qclusts(k,i) = qclusav(k,i)
                  qclusav(k,i) = 0.d0

               end do

            end do

               do k = 0, 20

                  numpal(k) = numsav(k)
                  rumpal(k) = rumsav(k)

               end do

         end if

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine pyildxyz(m,mz,mn,mm,nf,nl,lt, ! frtati 2022/02/18 added mm
     &                    nx,ny,nz,nn,xm,ym,zm,nt,ikzz,iknn,
     &                    tr,igsh,idasa)
*                                                                      *
*       output of nuclear yield (or production) tally in xyz mesh      *
*       last modified by K.Niita on 2004/12/27                         *
*                                                                      *
************************************************************************
      use sangelmod, only: itsans

*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'
      include 'param01.inc'

*-----------------------------------------------------------------------

      common /taliin/ rsouin, nzztin, nrgnin

*-----------------------------------------------------------------------

      common /verjam/ versn, lastr, iyeav, imonv, idayv

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

      common /tall21/ rtfac(itlmax)
      common /tall28/ itgsh(itlmax)
      common /tall22/ itmcn(itlmax), itmtn(itlmax), itmtt(itlmax)
      common /tall34/ itvm(itlmax)
      common /tall35/ iteps(itlmax)
      common /tall38/ itres(itlmax)
      common /tall39/ rtwid(itlmax)
      common /tall46/ itmtr(itlmax,4), rtmtr(itlmax,13)
      common /tall49/ itglt(itlmax)
      common /tall83/ itnzn(itlmax), itndm(itlmax)

      common /istcut/ ist_cut, ist_bat

      common /fact01/ facmax(itlmax) ! kitamura23/03/31

*-----------------------------------------------------------------------

      common /tcntl/  icntl, inucr

*-----------------------------------------------------------------------

      common /cparm/  maxbch,maxcas
      common /jcomon/ nabov,nobch,nocas,nomax
!$OMP THREADPRIVATE(/jcomon/)

      common /talout/ itall

      character fname*100, fnume*3

*-----------------------------------------------------------------------
      character gfnam*100
      common /tall60/ itger(itlmax)

      common /tall67/ itstd(itlmax), rtstd(itlmax)

      common /mpi00/ npe, me

*-----------------------------------------------------------------------

      dimension   lt(nl)
      dimension   xm(nx+1)
      dimension   ym(ny+1)
      dimension   zm(nz+1)
      dimension   nt(nn)
      dimension   ikzz(maxpt,maxnt), iknn(maxpt,maxnt)
      dimension   tm(maxpt+maxnt,2)
      dimension   fm(nf,nf)
      dimension   tr(nx,ny,nz,mz,mn,0:mm,2) ! frtati 2022/02/18 added mm

      integer,allocatable :: ixyz(:)

*-----------------------------------------------------------------------

      character hsunit(2)*15

      data hsunit / '[1/source]     ',
     &              '[1/cm^3/source]'/

*-----------------------------------------------------------------------

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

      character dc2*4

CCSE add to solve undefined parameter (2018.07.31) >>>>>
      data ipstep / 12 /
CCSE add to solve undefined parameter (2018.07.31) <<<<<
CCSE add for mesh=xyz parameter (2017.11.30) >>>>>
      character erfnm*100
CCSE add for mesh=xyz parameter (2017.11.30) <<<<<

      character chau*8
      character cha*1
      data cha /"'"/

      character rpa*1
      data rpa /'}'/
      character yen*1

*-----------------------------------------------------------------------

      integer :: itbmp
      common /tall63/ itbmp(itlmax)

      integer :: bmpWidth, bmpHeight
      character(1), allocatable :: bmpfIType(:)
      integer, allocatable :: bmpfIndex(:)
      integer :: numIndex

      integer :: itvtk, itvtkfmt
      common /tall64/ itvtk(itlmax),itvtkfmt(itlmax)

      integer :: isunit_vtk_default = 91
      integer :: isunit_vtk_meta_default = 92
      integer :: isunit_vtk_rm_default = 93
      integer :: isunit_vtk_geom_default = 94
      integer :: isunit_vtk_geom_meta_default = 95
      integer :: iunit_vtk_g_default = 96
      integer :: isunit_vtk = 0
      integer :: isunit_vtk_meta = 0
      integer :: isunit_vtk_rm = 0
      integer :: isunit_vtk_geom = 0
      integer :: isunit_vtk_geom_meta = 0
      integer :: iunit_vtk_g = 0
      integer :: ios

      character(len=255) :: outFilename
      logical :: isText

      double precision, allocatable :: fmval(:,:,:,:)
*-----------------------------------------------------------------------
      character(28) cfmt                !FURUTA20200615
      character(12) cir                 !FURUTA20200615
      common /redufmt/ iredufmt(itlmax) !FURUTA20200615
*-----------------------------------------------------------------------

*-----------------------------------------------------------------------
*        set mesh volume
*-----------------------------------------------------------------------

         vl(ix,iy,iz) = vls(nl,lt,itmcn(m),itvm(m),itmtr(m,4),
     &                      xm(ix),xm(ix+1),
     &                      ym(iy),ym(iy+1),
     &                      zm(iz),zm(iz+1))

*-----------------------------------------------------------------------

      yen  = char(92)
      allocate( ixyz(max((nx+1)*(ny+1),(ny+1)*(nz+1),(nx+1)*(nz+1))) )

*-----------------------------------------------------------------------
*        c1 : nomalization for source
*-----------------------------------------------------------------------

                  cmax = 0.0
                  cmin = 1.e+33
                  isdz = 0
                  stdm = 0.0
                  itstd(m) = 1

         if( igsh .eq. 0 ) then

            if ( rsouin .gt. 0d0 ) then ! S.H. (2016.7.26)

               c1 = 1.0d+0 / rsouin

            else

               c1 = 0d0

            end if

            facmax(m) = 1.d0
            if( rtfac(m) .lt. 0.d0 ) then
               facmax(m) = 0.d0

               do 101 jx = 1, nx
               do 101 jy = 1, ny
               do 101 jz = 1, nz
               do 101 iz = 1, mz
               do 101 in = 1, mn
               do 101 il = 0, mm

                  if( tr(jx,jy,jz,iz,in,il,1) .ne. 0.d0 ) then

                     if( itunt(m) .eq. 1 ) then

                        fmaxfc = tr(jx,jy,jz,iz,in,il,1)

                     else if( itunt(m) .eq. 2 ) then

                        fmaxfc = tr(jx,jy,jz,iz,in,il,1) / vl(jx,jy,jz)

                     end if

                     if( facmax(m) .lt. fmaxfc) facmax(m) = fmaxfc

                  end if

  101          continue

            end if

            do 100 jx = 1, nx
            do 100 jy = 1, ny
            do 100 jz = 1, nz

               if( itunt(m) .eq. 1 ) then

                  cc = abs(rtfac(m)/facmax(m))

               else if( itunt(m) .eq. 2 ) then

                  cc = abs(rtfac(m)/facmax(m)) / vl(jx,jy,jz)

               end if

            do 100 iz = 1, mz
            do 100 in = 1, mn
            do 100 il = 0, mm ! frtati 2022/02/18

               if( tr(jx,jy,jz,iz,in,il,1) .ne. 0.d0 ) then

                  call calc_stdev(m,Xa,sigx,
     &                            tr(jx,jy,jz,iz,in,il,1),
     &                            tr(jx,jy,jz,iz,in,il,2),
     &                            cc)

                  tr(jx,jy,jz,iz,in,il,1) = Xa
                  tr(jx,jy,jz,iz,in,il,2) = sigx
                  if( sigx .gt. stdm ) stdm = sigx

                  if( tr(jx,jy,jz,iz,in,il,1) .gt. cmax )
     &                                   cmax = tr(jx,jy,jz,iz,in,il,1)

cABE 2022/03/02, avoid the negative value of cmin
                  if( tr(jx,jy,jz,iz,in,il,1) .gt. 0.d0 .and.
     &                tr(jx,jy,jz,iz,in,il,1) .lt. cmin )
     &                                   cmin = tr(jx,jy,jz,iz,in,il,1)

               else

                  if( il .eq. 0 ) isdz = 1
                  tr(jx,jy,jz,iz,in,il,2) = 0.0

               end if

  100       continue

            if( nobch .gt. ist_bat ) then
            if( rtstd(m) .gt. 0.d0 .and.
     &          stdm .gt. 0.d0 .and. isdz .eq. 0 .and.
     &          stdm .lt. rtstd(m) ) itstd(m) = 0
            end if

         end if

*-----------------------------------------------------------------------
*     out put unit = 31  : temporary number
*-----------------------------------------------------------------------

      do 900 iax = 1, itfln(m)

!OBINATA(2012.8.20): output *.err
        noe = 1
        if ( any( itaxs(m,iax) .eq. (/ 8, 9, 10, 11, 13 /) ) ! H.Ratliff 2020.04.09
     &      .and. ittwo(m) .ne. 4 .and. icntl .ne. 8 )
     &    noe = 2

        if(itaxs(m,iax).eq.13.and.iredufmt(m).eq.1) noe=1 !FURUTA20200615

        do 900 ioe = 1, noe

         if( igsh .ne. 0 .and.
     &     ( itaxs(m,iax) .lt. 7 .or.
     &       ittwo(m) .eq. 4 .or. ittwo(m) .eq. 5 ) ) goto 900

         if( itall .eq. 2 .and. nobch .lt. maxbch .and.
     &       igsh .eq. 0 ) then

            write(fnume,'(i3.3)') nobch
!OBINATA(2012.8.20): output *.err
            if ( ioe .eq. 1 ) then
              fname = ctfln(m,iax)(1:itfll(m,iax))//'.'//fnume
            else
              call mk_2derrfn(ctfln(m,iax),fname,itfll(m,iax))
            end if

         else if ( itall.eq.4 .and. igsh.eq.0 ) then
            write(fnume,'(i3.3)') nobch
            if ( npe.gt.1 ) write (fnume,'(i3.3)') nobch/(npe-1)
            if ( ioe.eq.1 ) then
              call mk_2dnumfn(ctfln(m,iax),fname,itfll(m,iax),
     &             nobch,maxbch,npe)
            else
              call mk_2dnerfn(ctfln(m,iax),fname,itfll(m,iax),fnume)
            end if

         else

!OBINATA(2012.8.20): output *.err
            if ( ioe .eq. 1 ) then
              fname = ctfln(m,iax)
            else
              call mk_2derrfn(ctfln(m,iax),fname,itfll(m,iax))
            end if

         end if

            iot = 31
            open(iot, file = fname, status = 'unknown' )

*-----------------------------------------------------------------------
            gfnam = ctfln(m,iax)
            igfmn = itfll(m,iax)
            igser = itger(m)

*-----------------------------------------------------------------------
*     input echo
*-----------------------------------------------------------------------

               call tyilech(iot,m,iax,1)

CCSE add for mesh=xyz parameter (2017.11.30) >>>>>
*-----------------------------------------------------------------------
*        dchain axis
*-----------------------------------------------------------------------

         if( itaxs(m,iax) .eq. 13 ) then

*-----------------------------------------------------------------------
*           output error or not
*-----------------------------------------------------------------------

            if( itout(m) .ne. 0 ) then

               do i = itfll(m,iax), 1, -1

                  if( ctfln(m,iax)(i:i) .eq. '.' ) goto 50

               end do

               i = itfll(m,iax)

   50          itfp = i - 1

               do i = 1, itfp

                  erfnm(i:i) = ctfln(m,iax)(i:i)

               end do

               erfnm(itfp+1:itfp+4) = '.err'

               iou = 15
               open(iou, file = erfnm(1:itfp+4), status = 'unknown' )

               call tyilech(iou,m,iax,1)

            end if

*-----------------------------------------------------------------------
*        output : do not change below expresion
*-----------------------------------------------------------------------

            write(iot,'(/"#",78("-"))')

            write(iot,'(/
     &           " xyz scoring mesh nuclear yield (or production)"/
     &           " ----------------------------------------")')

            if( itout(m) .ne. 0 ) then

               write(iou,'(/"#",78("-"))')

               write(iou,'(/
     &         " Statistical Error(%) for xyz scoring mesh nuclear ",
     &         "yield of above file."/
     &         " --------------------------------------------",
     &         "---------------------")')

            end if

           if(iredufmt(m).eq.0)then !FURUTA20200615
*-----------------------------------------------------------------------

            do 170 il = 0, 2 ! frtati 2022/03/11
               do 170 iz = 1, maxpt

                  if( nn .gt. 0 ) then

                     do i = 1, nn

                        if( nt(i) / 1000 .eq. iz ) goto 150

                     end do

                     goto 170

                  end if

  150             continue

                  do in = 1, maxnt
                     do jx = 1, nx
                     do jy = 1, ny
                     do jz = 1, nz

                        if( itnzn(m).eq.0 ) then
                          t0_tr = tr(jx,jy,jz,iz,in,il,1)
                        else
                          t0_tr =
     &                    tr(jx,jy,jz,igetiznmp(iz,in,il,m),1,0,1)
                        end if
                        if( t0_tr .ne. 0.0 ) goto 120

                     end do
                     end do
                     end do
                  end do

                  goto 170

  120             im = in

                  do in = maxnt, im + 1, -1
                     do jx = 1, nx
                     do jy = 1, ny
                     do jz = 1, nz

                        if( itnzn(m).eq.0 ) then
                          t0_tr = tr(jx,jy,jz,iz,in,il,1)
                        else
                          t0_tr =
     &                    tr(jx,jy,jz,igetiznmp(iz,in,il,m),1,0,1)
                        end if
                        if( t0_tr .ne. 0.0 ) goto 140

                     end do
                     end do
                     end do
                  end do

                  jm = im
  140             jm = in

                  km = jm - im + 1
                  lm = ( km - 1 ) / ipstep + 1

                  do mmm = 1, lm

                     n1 = ipstep * ( mmm - 1 ) + 1
                     n2 = min( ipstep * mmm, km )
                     n3 = n1 + im - 1
                     n4 = n2 + im - 1

                     IF(il .eq. 0) then
                        write(iot,'(/1x,i4,"-",a2,
     &                           " isotope production",
     &                           " # n3 n4 = ",i3,1x,i3)')
     &                  iz, elmnt(iz), n3, n4

                     ELSEIF(il .eq. 1) then
                        write(iot,'(/1x,i4,"-",a2,
     &                           " 1st metastable isotope production",
     &                           " # n3 n4 = ",i3,1x,i3)')
     &                 iz, elmnt(iz), n3, n4

                     ELSEIF(il .eq. 2) then
                        write(iot,'(/1x,i4,"-",a2,
     &                           " 2nd metastable isotope production",
     &                           " # n3 n4 = ",i3,1x,i3)')
     &                 iz, elmnt(iz), n3, n4
                     ENDIF

                     IF(il .eq. 0) then
                        write(iot,'("     jx","     jy","     jz"
     &                       ,12i11)')
     &                       ( im - 1 + iz + n, n = n1, n2 )
                     ELSEIF(il .eq. 1) then
                        write(iot,'("     jx","     jy","     jz"
     &                ,i10,"m",i10,"m",i10,"m",i10,"m",i10,"m"
     &                ,i10,"m",i10,"m",i10,"m",i10,"m",i10,"m"
     &                ,i10,"m",i10,"m")')
     &                ( im - 1 + iz + n, n = n1, n2 )
                     ELSEIF(il .eq. 2) then
                        write(iot,'("     jx","     jy","     jz"
     &                ,i10,"n",i10,"n",i10,"n",i10,"n",i10,"n"
     &                ,i10,"n",i10,"n",i10,"n",i10,"n",i10,"n"
     &                ,i10,"n",i10,"n")')
     &                ( im - 1 + iz + n, n = n1, n2 )
                     ENDIF

                     do jx = 1, nx
                     do jy = 1, ny
                     do jz = 1, nz
                      if( itnzn(m).eq.0 ) then ! frtati 2022/02/18
                        write(iot,'(3i7,1p12e11.3)')
     &                 jx,jy,jz, (tr(jx,jy,jz,iz,i,il,ioe), i = n3, n4 )
                      else
                        write(iot,'(3i7,1p12e11.3)') jx,jy,jz,
     &                  (tr(jx,jy,jz,igetiznmp(iz,i,il,m),1,0,ioe),
     &                   i = n3, n4 )
                      end if
                     end do
                     end do
                     end do

                     if( itout(m) .ne. 0 ) then

                        IF(il .eq. 0) then
                           write(iou,'(/1x,i4,"-",a2,
     &                             " isotope production : ERROR(%)")')
     &                         iz, elmnt(iz)

                        ELSEIF(il .eq. 1) then
                           write(iou,'(/1x,i4,"-",a2,
     &              " 1st metastable isotope production : ERROR(%)")')
     &              iz, elmnt(iz)

                        ELSEIF(il .eq. 2) then
                           write(iou,'(/1x,i4,"-",a2,
     &              " 2nd metastable isotope production : ERROR(%)")')
     &              iz, elmnt(iz)
                        ENDIF

                        IF(il .eq. 0) then
                       write(iou,'("     jx","     jy","     jz"
     &                     ,12i11)')
     &                      ( im - 1 + iz + n, n = n1, n2 )
                        ELSEIF(il .eq. 1) then
                       write(iou,'("     jx","     jy","     jz"
     &                ,i10,"m",i10,"m",i10,"m",i10,"m",i10,"m"
     &                ,i10,"m",i10,"m",i10,"m",i10,"m",i10,"m"
     &                ,i10,"m",i10,"m")')
     &                      ( im - 1 + iz + n, n = n1, n2 )
                        ELSEIF(il .eq. 2) then
                       write(iou,'("     jx","     jy","     jz"
     &                ,i10,"n",i10,"n",i10,"n",i10,"n",i10,"n"
     &                ,i10,"n",i10,"n",i10,"n",i10,"n",i10,"n"
     &                ,i10,"n",i10,"n")')
     &                      ( im - 1 + iz + n, n = n1, n2 )
                        ENDIF

                        do jx = 1, nx
                        do jy = 1, ny
                        do jz = 1, nz
                         if( itnzn(m).eq.0 ) then ! frtati 2022/02/18
                           write(iou,'(3i7,1p12e11.3)')
     &             jx,jy,jz, (tr(jx,jy,jz,iz,i,il,2)*100.0, i = n3, n4 )
                         else
                           write(iou,'(3i7,1p12e11.3)') jx,jy,jz,
     &                     (tr(jx,jy,jz,igetiznmp(iz,i,il,m),1,0,2)
     &                      *100.0,i = n3, n4 )
                         end if
                        end do
                        end do
                        end do

                     end if

                  end do

  170       continue

*-----------------------------------------------------------------------
           else
             write(iot,'(a)')
             write(iot,'("# num nucleusID yield r.err")')
             if(itout(m).ne.0)then
              write(iou,'(a)')
              write(iou,'("# num nucleusID yield r.err")')
             endif
             cfmt='(i#,x,i7,1p2e11.3)'
             do jz=1,nz
             do jy=1,ny
             do jx=1,nx
              ir=jx+(jy-1)*nx+(jz-1)*nx*ny
!NS 2021.04 change for NVIDIA TOOL KIT COMPILL ERR
              write(cir,"(I12)")ir
!NS 2021.04 end change for NVIDIA TOOL KIT COMPILL ERR
              nir=len_trim(adjustl(cir))
              write(cfmt(3:3),'(i1)')nir
              do iz=1,maxpt
               if(nn.gt.0) then
                do i = 1,nn
                 if(nt(i)/1000.eq.iz)exit
                enddo
                if(i.gt.nn)cycle
               endif
               do in=1,maxnt
                do il=0,2 ! frtati 2022/05/02
                 if( itnzn(m).eq.0 ) then ! frtati 2022/02/18
                  if(tr(jx,jy,jz,iz,in,il,ioe).ne.0d0)then
                   write(iot,cfmt)ir,iz*10000+(in+iz)*10+il,
     &                  (tr(jx,jy,jz,iz,in,il,ioee),ioee=1,2)
                   if(itout(m).ne.0)then
                    write(iou,cfmt)ir,iz*10000+(in+iz)*10+il,
     &                   tr(jx,jy,jz,iz,in,il,2)
                   endif
                  endif
                 else
                  if(tr(jx,jy,jz,igetiznmp(iz,in,il,m),1,0,ioe)
     &               .ne.0d0)then
                   write(iot,cfmt)ir,iz*10000+(in+iz)*10+il,
     &             (tr(jx,jy,jz,igetiznmp(iz,in,il,m),1,0,ioee),
     &              ioee=1,2)
                   if(itout(m).ne.0)then
                    write(iou,cfmt)ir,iz*10000+(in+iz)*10+il,
     &              tr(jx,jy,jz,igetiznmp(iz,in,il,m),1,0,2)
                   endif
                  endif
                 endif
                enddo
               enddo
              enddo
             enddo
             enddo
             enddo
             write(iot,cfmt)0,0,0.0d0,0.0d0
             if(itout(m).ne.0)write(iou,cfmt)0,0,0.0d0,0.0d0
           endif
*-----------------------------------------------------------------------

            if( itout(m) .ne. 0 ) close(iou)

CCSE add for mesh=xyz parameter (2017.11.30) <<<<<

*-----------------------------------------------------------------------
*        x axis
*-----------------------------------------------------------------------

CCSE change for mesh=xyz parameter (2017.11.30) >>>>>
         else if( itaxs(m,iax) .eq. 3 ) then
CCSE change for mesh=xyz parameter (2017.11.30) <<<<<

               inum = 0

                  if( nn .eq. 0 ) then

                     nc = 1

                  else

                     nc = nn

                  end if

*-----------------------------------------------------------------------

            do jy = 1, ny
            do jz = 1, nz

            do ic = 1, nc

               inum = inum + 1

*-----------------------------------------------------------------------

               if( nn .gt. 0 ) then

                  iz = nt(ic) / 1000

                  ia = nt(ic) - iz * 1000

                  if( ia .gt. 0 ) then

                     in = ia - iz

                     call chname(idum,ia,iz,chau)

                  end if

               end if

*-----------------------------------------------------------------------

                  write(iot,'(/"#",78("-"))')

               if( inum .eq. 1 ) then
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

               if( rtfac(m) .lt. 0.d0 )
     &           write(iot,
     &           '( "#restart maximum for normalization:",1p1e14.7)')
     &           facmax(m)

            if( nn .eq. 0 ) then

               write(iot,'("#   no. =",i3,3x,
     &         "iy  =",i3,3x,
     &         "iz  =",i3,3x,"Z = all"/
     &         "#   y = (",1p1e13.4,"  -",1p1e13.4,"  )"/
     &         "#   z = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     inum, jy, jz,
     &                     ym(jy), ym(jy+1),
     &                     zm(jz), zm(jz+1)

            else if( ia .eq. 0 ) then

               write(iot,'("#   no. =",i3,3x,
     &         "iy  =",i3,3x,
     &         "iz  =",i3,3x,"Z = ",i3," : ",a3/
     &         "#   y = (",1p1e13.4,"  -",1p1e13.4,"  )"/
     &         "#   z = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     inum, jy, jz, iz, elmnt(iz),
     &                     ym(jy), ym(jy+1),
     &                     zm(jz), zm(jz+1)

            else

               write(iot,'("#   no. =",i3,3x,
     &         "iy  =",i3,3x,
     &         "iz  =",i3,3x,a8/
     &         "#   y = (",1p1e13.4,"  -",1p1e13.4,"  )"/
     &         "#   z = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     inum, jy, jz, chau,
     &                     ym(jy), ym(jy+1),
     &                     zm(jz), zm(jz+1)

            end if

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: x [cm]")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: Number ",a15)') hsunit(itunt(m))

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

                  write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')

               if( itanl(m) .gt. 0 ) then

                  write(iot,'( "p: ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               end if

               if( itsans(m) .gt. 0 ) then

                  call write_sangel(iot,m,0)

               end if

c S.H. added IF statement below for epsout=2 (2017.2.18)
              if ( iteps(m) .ne. 2 ) then

               write(iot,'( "h: n",12x,"x",12x,
     &                    "y,hhl0      n")')

              else

               write(iot,'( "h: n",12x,"x",12x,
     &                    "y1,hhl0     ny2 dy1=[y1*y2]")')

              end if

               write(iot,'( "#  x-lower      x-upper",3x,
     &                    "   number      r.err")')

*-----------------------------------------------------------------------

                     seka = 0.0
                     sera = 0.0
                     seva = 0.0

               do jx = 1, nx

                     sek = 0.0
                     ser = 0.0

                     if( itunt(m) .eq. 1 ) then

                        vm = 1.0

                     else if( itunt(m) .eq. 2 ) then

                        vm = vl(jx,jy,jz)

                     end if

*-----------------------------------------------------------------------

                  if( nn .eq. 0 ) then

                     do kn = 1, mn
                     do kz = 1, mz
                     do il = 0, mm ! frtati 2022/02/18

                        lz  = ikzz( kz, kn )
                        ln  = iknn( kz, kn )

                        vn  = vm * tr(jx,jy,jz,lz,ln,il,1)
                        sek = sek + vn
                        ser = ser + ( vn * tr(jx,jy,jz,lz,ln,il,2) )**2

                     end do
                     end do
                     end do

                  else if( ia .eq. 0 ) then

                     do kn = 1, mn
                     do il = 0, mm ! frtati 2022/02/18

                        lz  = ikzz( iz, kn )
                        ln  = iknn( iz, kn )

                        vn  = vm * tr(jx,jy,jz,lz,ln,il,1)
                        sek = sek + vn
                        ser = ser + ( vn * tr(jx,jy,jz,lz,ln,il,2) )**2

                     end do
                     end do

                  else

                        lz  = ikzz( iz, in )
                        ln  = iknn( iz, in )

                     do il = 0, mm ! frtati 2022/02/18
                        vn  = vm * tr(jx,jy,jz,lz,ln,il,1)
                        sek = sek + vn
                        ser = ser
     &                      + ( vn * tr(jx,jy,jz,lz,ln,il,2) )**2
                     end do

                  end if

*-----------------------------------------------------------------------

                        seka = seka + sek
                        sera = sera + ser
                        seva = seva + vm

                     if( sek .ne. 0.0 ) then ! .gt. -> .ne. by T.Sato 2024/05/02 for negative sek

                        ser = sqrt( ser ) / sek
                        sek = sek / vm

                     end if

                  write(iot,'(1p2e13.4,1pe13.4,0pf8.4)')
     &                  xm(jx), xm(jx+1), sek, ser

               end do

                     if( seka .gt. 0.0 ) then

                        sera = sqrt( sera ) / seka

                        if( itunt(m) .eq. 2 )
     &                  seka = seka / seva

                     end if

               write(iot,'(/"#   sum over",14x,1pe13.4,0pf8.4)')
     &               seka, sera

            if( itunt(m) .eq. 2 ) then

               write(iot,'(
     &         "#   vol sum = ",1p1e13.4,"  [cm^3]")') seva

            end if

               if( nn .eq. 0 ) then

                  write(iot,'(/a1,"no. =",i3,
     &                        ",    iy =",i3,
     &                        ",    iz =",i3,
     &                        ",    Z  =  all", a1)')
     &                        cha, inum, jy, jz, cha

               else if( ia .eq. 0 ) then

                  write(iot,'(/a1,"no. =",i3,
     &                        ",    iy =",i3,
     &                        ",    iz =",i3,
     &                        ",    Z  =  ",i3, a1)')
     &                        cha, inum, jy, jz, iz, cha

               else

                  write(iot,'(/a1,"no. =",i3,
     &                        ",    iy =",i3,
     &                        ",    iz =",i3,
     &                        ",    ",a8, a1)')
     &                        cha, inum, jy, jz, chau, cha

               end if

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

               if( nn .eq. 0 ) then

                  write(iot,'("wt: s(0.7)",/
     &                        a1,"vspace{-3}"/
     &                        "  ymin  &=&",1pe13.4," [cm]"/
     &                        "  ymax  &=&",1pe13.4," [cm]"/
     &                        "  zmin  &=&",1pe13.4," [cm]"/
     &                        "  zmax  &=&",1pe13.4," [cm]"/
     &                        "      Z &=&   all"/
     &                        "e:")')
     &                        yen,
     &                        ym(jy), ym(jy+1),
     &                        zm(jz), zm(jz+1)

               else if( ia .eq. 0 ) then

                  write(iot,'("wt: s(0.7)",/
     &                        a1,"vspace{-3}"/
     &                        "  ymin  &=&",1pe13.4," [cm]"/
     &                        "  ymax  &=&",1pe13.4," [cm]"/
     &                        "  zmin  &=&",1pe13.4," [cm]"/
     &                        "  zmax  &=&",1pe13.4," [cm]"/
     &                        "      Z &=& ",i3,"   :  ",a3/
     &                        "e:")')
     &                        yen,
     &                        ym(jy), ym(jy+1),
     &                        zm(jz), zm(jz+1), iz, elmnt(iz)

               else

                  write(iot,'("wt: s(0.7)",/
     &                        a1,"vspace{-3}"/
     &                        "  ymin  &=&",1pe13.4," [cm]"/
     &                        "  ymax  &=&",1pe13.4," [cm]"/
     &                        "  zmin  &=&",1pe13.4," [cm]"/
     &                        "  zmax  &=&",1pe13.4," [cm]"/
     &                        "  nucl. &=&   ",a8/
     &                        "e:")')
     &                        yen,
     &                        ym(jy), ym(jy+1),
     &                        zm(jz), zm(jz+1), chau

               end if

*-----------------------------------------------------------------------

            end do
            end do
            end do

*-----------------------------------------------------------------------
*        y axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 4 ) then

               inum = 0

                  if( nn .eq. 0 ) then

                     nc = 1

                  else

                     nc = nn

                  end if

*-----------------------------------------------------------------------

            do jx = 1, nx
            do jz = 1, nz

            do ic = 1, nc

               inum = inum + 1

*-----------------------------------------------------------------------

               if( nn .gt. 0 ) then

                  iz = nt(ic) / 1000

                  ia = nt(ic) - iz * 1000

                  if( ia .gt. 0 ) then

                     in = ia - iz

                     call chname(idum,ia,iz,chau)

                  end if

               end if

*-----------------------------------------------------------------------

                  write(iot,'(/"#",78("-"))')

               if( inum .eq. 1 ) then
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

               if( rtfac(m) .lt. 0.d0 )
     &           write(iot,
     &           '( "#restart maximum for normalization:",1p1e14.7)')
     &           facmax(m)

            if( nn .eq. 0 ) then

               write(iot,'("#   no. =",i3,3x,
     &         "ix  =",i3,3x,
     &         "iz  =",i3,3x,"Z = all"/
     &         "#   x = (",1p1e13.4,"  -",1p1e13.4,"  )"/
     &         "#   z = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     inum, jx, jz,
     &                     xm(jx), xm(jx+1),
     &                     zm(jz), zm(jz+1)

            else if( ia .eq. 0 ) then

               write(iot,'("#   no. =",i3,3x,
     &         "ix  =",i3,3x,
     &         "iz  =",i3,3x,"Z = ",i3," : ",a3/
     &         "#   x = (",1p1e13.4,"  -",1p1e13.4,"  )"/
     &         "#   z = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     inum, jx, jz, iz, elmnt(iz),
     &                     xm(jx), xm(jx+1),
     &                     zm(jz), zm(jz+1)

            else

               write(iot,'("#   no. =",i3,3x,
     &         "ix  =",i3,3x,
     &         "iz  =",i3,3x,a8/
     &         "#   x = (",1p1e13.4,"  -",1p1e13.4,"  )"/
     &         "#   z = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     inum, jx, jz, chau,
     &                     xm(jx), xm(jx+1),
     &                     zm(jz), zm(jz+1)

            end if

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: y [cm]")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: Number ",a15)') hsunit(itunt(m))

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

                  write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')

               if( itanl(m) .gt. 0 ) then

                  write(iot,'( "p: ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               end if

               if( itsans(m) .gt. 0 ) then

                  call write_sangel(iot,m,0)

               end if

c S.H. added IF statement below for epsout=2 (2017.2.18)
              if ( iteps(m) .ne. 2 ) then

               write(iot,'( "h: n",12x,"x",12x,
     &                    "y,hhl0      n")')

              else

               write(iot,'( "h: n",12x,"x",12x,
     &                    "y1,hhl0     ny2 dy1=[y1*y2]")')

              end if

               write(iot,'( "#  y-lower      y-upper",3x,
     &                    "   number      r.err")')

*-----------------------------------------------------------------------

                     seka = 0.0
                     sera = 0.0
                     seva = 0.0

               do jy = 1, ny

                     sek = 0.0
                     ser = 0.0

                     if( itunt(m) .eq. 1 ) then

                        vm = 1.0

                     else if( itunt(m) .eq. 2 ) then

                        vm = vl(jx,jy,jz)

                     end if

*-----------------------------------------------------------------------

                  if( nn .eq. 0 ) then

                     do kn = 1, mn
                     do kz = 1, mz
                     do il = 0, mm ! frtati 2022/02/18

                        lz  = ikzz( kz, kn )
                        ln  = iknn( kz, kn )

                        vn  = vm * tr(jx,jy,jz,lz,ln,il,1)
                        sek = sek + vn
                        ser = ser + ( vn * tr(jx,jy,jz,lz,ln,il,2) )**2

                     end do
                     end do
                     end do

                  else if( ia .eq. 0 ) then

                     do kn = 1, mn
                     do il = 0, mm ! frtati 2022/02/18

                        lz  = ikzz( iz, kn )
                        ln  = iknn( iz, kn )

                        vn  = vm * tr(jx,jy,jz,lz,ln,il,1)
                        sek = sek + vn
                        ser = ser + ( vn * tr(jx,jy,jz,lz,ln,il,2) )**2

                     end do
                     end do

                  else

                        lz  = ikzz( iz, in )
                        ln  = iknn( iz, in )

                     do il = 0, mm ! frtati 2022/02/18
                        vn  = vm * tr(jx,jy,jz,lz,ln,il,1)
                        sek = sek + vn
                        ser = ser
     &                      + ( vn * tr(jx,jy,jz,lz,ln,il,2) )**2
                     end do

                  end if

*-----------------------------------------------------------------------

                        seka = seka + sek
                        sera = sera + ser
                        seva = seva + vm

                     if( sek .ne. 0.0 ) then ! .gt. -> .ne. by T.Sato 2024/05/02 for negative sek

                        ser = sqrt( ser ) / sek
                        sek = sek / vm

                     end if

                  write(iot,'(1p2e13.4,1pe13.4,0pf8.4)')
     &                  ym(jy),ym(jy+1), sek, ser

               end do

                     if( seka .gt. 0.0 ) then

                        sera = sqrt( sera ) / seka

                        if( itunt(m) .eq. 2 )
     &                  seka = seka / seva

                     end if

               write(iot,'(/"#   sum over",14x,1pe13.4,0pf8.4)')
     &               seka, sera

            if( itunt(m) .eq. 2 ) then

               write(iot,'(
     &         "#   vol sum = ",1p1e13.4,"  [cm^3]")') seva

            end if

               if( nn .eq. 0 ) then

                  write(iot,'(/a1,"no. =",i3,
     &                        ",    ix =",i3,
     &                        ",    iz =",i3,
     &                        ",    Z  =  all", a1)')
     &                        cha, inum, jx, jz, cha

               else if( ia .eq. 0 ) then

                  write(iot,'(/a1,"no. =",i3,
     &                        ",    ix =",i3,
     &                        ",    iz =",i3,
     &                        ",    Z  =  ",i3, a1)')
     &                        cha, inum, jx, jz, iz, cha

               else

                  write(iot,'(/a1,"no. =",i3,
     &                        ",    ix =",i3,
     &                        ",    iz =",i3,
     &                        ",    ",a8, a1)')
     &                        cha, inum, jx, jz, chau, cha

               end if

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

               if( nn .eq. 0 ) then

                  write(iot,'("wt: s(0.7)",/
     &                        a1,"vspace{-3}"/
     &                        "  xmin  &=&",1pe13.4," [cm]"/
     &                        "  xmax  &=&",1pe13.4," [cm]"/
     &                        "  zmin  &=&",1pe13.4," [cm]"/
     &                        "  zmax  &=&",1pe13.4," [cm]"/
     &                        "      Z &=&   all"/
     &                        "e:")')
     &                        yen,
     &                        xm(jx), xm(jx+1),
     &                        zm(jz), zm(jz+1)

               else if( ia .eq. 0 ) then

                  write(iot,'("wt: s(0.7)",/
     &                        a1,"vspace{-3}"/
     &                        "  xmin  &=&",1pe13.4," [cm]"/
     &                        "  xmax  &=&",1pe13.4," [cm]"/
     &                        "  zmin  &=&",1pe13.4," [cm]"/
     &                        "  zmax  &=&",1pe13.4," [cm]"/
     &                        "      Z &=& ",i3,"   :  ",a3/
     &                        "e:")')
     &                        yen,
     &                        xm(jx), xm(jx+1),
     &                        zm(jz), zm(jz+1), iz, elmnt(iz)

               else

                  write(iot,'("wt: s(0.7)",/
     &                        a1,"vspace{-3}"/
     &                        "  xmin  &=&",1pe13.4," [cm]"/
     &                        "  xmax  &=&",1pe13.4," [cm]"/
     &                        "  zmin  &=&",1pe13.4," [cm]"/
     &                        "  zmax  &=&",1pe13.4," [cm]"/
     &                        "  nucl. &=&   ",a8/
     &                        "e:")')
     &                        yen,
     &                        xm(jx), xm(jx+1),
     &                        zm(jz), zm(jz+1), chau

               end if

*-----------------------------------------------------------------------

            end do
            end do
            end do

*-----------------------------------------------------------------------
*        z axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 5 ) then

               inum = 0

                  if( nn .eq. 0 ) then

                     nc = 1

                  else

                     nc = nn

                  end if

*-----------------------------------------------------------------------

            do jx = 1, nx
            do jy = 1, ny

            do ic = 1, nc

               inum = inum + 1

*-----------------------------------------------------------------------

               if( nn .gt. 0 ) then

                  iz = nt(ic) / 1000

                  ia = nt(ic) - iz * 1000

                  if( ia .gt. 0 ) then

                     in = ia - iz

                     call chname(idum,ia,iz,chau)

                  end if

               end if

*-----------------------------------------------------------------------

                  write(iot,'(/"#",78("-"))')

               if( inum .eq. 1 ) then
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

               if( rtfac(m) .lt. 0.d0 )
     &           write(iot,
     &           '( "#restart maximum for normalization:",1p1e14.7)')
     &           facmax(m)

            if( nn .eq. 0 ) then

               write(iot,'("#   no. =",i3,3x,
     &         "ix  =",i3,3x,
     &         "iy  =",i3,3x,"Z = all"/
     &         "#   x = (",1p1e13.4,"  -",1p1e13.4,"  )"/
     &         "#   y = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     inum, jx, jy,
     &                     xm(jx), xm(jx+1),
     &                     ym(jy), ym(jy+1)

            else if( ia .eq. 0 ) then

               write(iot,'("#   no. =",i3,3x,
     &         "ix  =",i3,3x,
     &         "iy  =",i3,3x,"Z = ",i3," : ",a3/
     &         "#   x = (",1p1e13.4,"  -",1p1e13.4,"  )"/
     &         "#   y = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     inum, jx, jy, iz, elmnt(iz),
     &                     xm(jx), xm(jx+1),
     &                     ym(jy), ym(jy+1)

            else

               write(iot,'("#   no. =",i3,3x,
     &         "ix  =",i3,3x,
     &         "iy  =",i3,3x,a8/
     &         "#   x = (",1p1e13.4,"  -",1p1e13.4,"  )"/
     &         "#   y = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     inum, jx, jy, chau,
     &                     xm(jx), xm(jx+1),
     &                     ym(jy), ym(jy+1)

            end if

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: z [cm]")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: Number ",a15)') hsunit(itunt(m))

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

                  write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')

               if( itanl(m) .gt. 0 ) then

                  write(iot,'( "p: ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               end if

               if( itsans(m) .gt. 0 ) then

                  call write_sangel(iot,m,0)

               end if

c S.H. added IF statement below for epsout=2 (2017.2.18)
              if ( iteps(m) .ne. 2 ) then

               write(iot,'( "h: n",12x,"x",12x,
     &                    "y,hhl0      n")')

              else

               write(iot,'( "h: n",12x,"x",12x,
     &                    "y1,hhl0     ny2 dy1=[y1*y2]")')

              end if

               write(iot,'( "#  z-lower      z-upper",3x,
     &                    "   number      r.err")')

*-----------------------------------------------------------------------

                     seka = 0.0
                     sera = 0.0
                     seva = 0.0

               do jz = 1, nz

                     sek = 0.0
                     ser = 0.0

                     if( itunt(m) .eq. 1 ) then

                        vm = 1.0

                     else if( itunt(m) .eq. 2 ) then

                        vm = vl(jx,jy,jz)

                     end if

*-----------------------------------------------------------------------

                  if( nn .eq. 0 ) then

                     do kn = 1, mn
                     do kz = 1, mz
                     do il = 0, mm ! frtati 2022/02/18

                        lz  = ikzz( kz, kn )
                        ln  = iknn( kz, kn )

                        vn  = vm * tr(jx,jy,jz,lz,ln,il,1)
                        sek = sek + vn
                        ser = ser + ( vn * tr(jx,jy,jz,lz,ln,il,2) )**2

                     end do
                     end do
                     end do

                  else if( ia .eq. 0 ) then

                     do kn = 1, mn
                     do il = 0, mm ! frtati 2022/02/18

                        lz  = ikzz( iz, kn )
                        ln  = iknn( iz, kn )

                        vn  = vm * tr(jx,jy,jz,lz,ln,il,1)
                        sek = sek + vn
                        ser = ser + ( vn * tr(jx,jy,jz,lz,ln,il,2) )**2

                     end do
                     end do

                  else

                        lz  = ikzz( iz, in )
                        ln  = iknn( iz, in )

                     do il = 0, mm ! frtati 2022/02/18
                        vn  = vm * tr(jx,jy,jz,lz,ln,il,1)
                        sek = sek + vn
                        ser = ser
     &                      + ( vn * tr(jx,jy,jz,lz,ln,il,2) )**2
                     end do

                  end if

*-----------------------------------------------------------------------

                        seka = seka + sek
                        sera = sera + ser
                        seva = seva + vm

                     if( sek .ne. 0.0 ) then ! .gt. -> .ne. by T.Sato 2024/05/02 for negative sek

                        ser = sqrt( ser ) / sek
                        sek = sek / vm

                     end if

                  write(iot,'(1p2e13.4,1pe13.4,0pf8.4)')
     &                  zm(jz),zm(jz+1), sek, ser

               end do

                     if( seka .gt. 0.0 ) then

                        sera = sqrt( sera ) / seka

                        if( itunt(m) .eq. 2 )
     &                  seka = seka / seva

                     end if

               write(iot,'(/"#   sum over",14x,1pe13.4,0pf8.4)')
     &               seka, sera

            if( itunt(m) .eq. 2 ) then

               write(iot,'(
     &         "#   vol sum = ",1p1e13.4,"  [cm^3]")') seva

            end if

               if( nn .eq. 0 ) then

                  write(iot,'(/a1,"no. =",i3,
     &                        ",    ix =",i3,
     &                        ",    iy =",i3,
     &                        ",    Z  =  all", a1)')
     &                        cha, inum, jx, jy, cha

               else if( ia .eq. 0 ) then

                  write(iot,'(/a1,"no. =",i3,
     &                        ",    ix =",i3,
     &                        ",    iy =",i3,
     &                        ",    Z  =  ",i3, a1)')
     &                        cha, inum, jx, jy, iz, cha

               else

                  write(iot,'(/a1,"no. =",i3,
     &                        ",    ix =",i3,
     &                        ",    iy =",i3,
     &                        ",    ",a8, a1)')
     &                        cha, inum, jx, jy, chau, cha

               end if

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

               if( nn .eq. 0 ) then

                  write(iot,'("wt: s(0.7)",/
     &                        a1,"vspace{-3}"/
     &                        "  xmin  &=&",1pe13.4," [cm]"/
     &                        "  xmax  &=&",1pe13.4," [cm]"/
     &                        "  ymin  &=&",1pe13.4," [cm]"/
     &                        "  ymax  &=&",1pe13.4," [cm]"/
     &                        "      Z &=&   all"/
     &                        "e:")')
     &                        yen,
     &                        xm(jx), xm(jx+1),
     &                        ym(jy), ym(jy+1)

               else if( ia .eq. 0 ) then

                  write(iot,'("wt: s(0.7)",/
     &                        a1,"vspace{-3}"/
     &                        "  xmin  &=&",1pe13.4," [cm]"/
     &                        "  xmax  &=&",1pe13.4," [cm]"/
     &                        "  ymin  &=&",1pe13.4," [cm]"/
     &                        "  ymax  &=&",1pe13.4," [cm]"/
     &                        "      Z &=& ",i3,"   :  ",a3/
     &                        "e:")')
     &                        yen,
     &                        xm(jx), xm(jx+1),
     &                        ym(jy), ym(jy+1), iz, elmnt(iz)

               else

                  write(iot,'("wt: s(0.7)",/
     &                        a1,"vspace{-3}"/
     &                        "  xmin  &=&",1pe13.4," [cm]"/
     &                        "  xmax  &=&",1pe13.4," [cm]"/
     &                        "  ymin  &=&",1pe13.4," [cm]"/
     &                        "  ymax  &=&",1pe13.4," [cm]"/
     &                        "  nucl. &=&   ",a8/
     &                        "e:")')
     &                        yen,
     &                        xm(jx), xm(jx+1),
     &                        ym(jy), ym(jy+1), chau

               end if

*-----------------------------------------------------------------------

            end do
            end do
            end do

*-----------------------------------------------------------------------
*        mass axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 1 ) then

               inum = 0

                  if( nn .eq. 0 ) then

                     nc = 1

                  else

                     nc = nn

                  end if

*-----------------------------------------------------------------------

            do jx = 1, nx
            do jy = 1, ny
            do jz = 1, nz

            do 260 ic = 1, nc

                  if( nn .gt. 0 ) iz = nt(ic) / 1000

                  do i = 1, maxnt + maxpt

                     tm(i,1) = 0.0d+0
                     tm(i,2) = 0.0d+0

                  end do

                  if( itunt(m) .eq. 1 ) then

                     vm = 1.0

                  else if( itunt(m) .eq. 2 ) then

                     vm = vl(jx,jy,jz)

                  end if

*-----------------------------------------------------------------------

                  if( nn .eq. 0 .and. itfln(m) .eq. 1 ) then

                     do ln = 1, mn
                     do il = 0, mm ! frtati 2022/02/18

                        lz = 1
                        im = ln

                        vn = vm * tr(jx,jy,jz,lz,ln,il,1)
                        tm(im,1) = tm(im,1) + vn
                        tm(im,2) = tm(im,2)
     &                           + ( vn * tr(jx,jy,jz,lz,ln,il,2) )**2

                     end do
                     end do

                  else if( nn .eq. 0 .and. itfln(m) .gt. 1 ) then

                     do kn = 1, mn
                     do kz = 1, mz
                     do il = 0, mm ! frtati 2022/02/18

                        lz = ikzz( kz, kn )
                        ln = iknn( kz, kn )

                        im = kz + kn

                        vn = vm * tr(jx,jy,jz,lz,ln,il,1)
                        tm(im,1) = tm(im,1) + vn
                        tm(im,2) = tm(im,2)
     &                           + ( vn * tr(jx,jy,jz,lz,ln,il,2) )**2

                     end do
                     end do
                     end do

                  else

                     do kn = 1, mn
                     do il = 0, mm ! frtati 2022/02/18

                        lz = ikzz( iz, kn )
                        ln = iknn( iz, kn )

                        im = iz + kn

                        vn = vm * tr(jx,jy,jz,lz,ln,il,1)
                        tm(im,1) = tm(im,1) + vn
                        tm(im,2) = tm(im,2)
     &                           + ( vn * tr(jx,jy,jz,lz,ln,il,2) )**2

                     end do
                     end do

                  end if

*-----------------------------------------------------------------------

               do i = 1, maxnt + maxpt

                  if( tm(i,1) .gt. 0.0d0 ) then

                     tm(i,2) = sqrt( tm(i,2) ) / tm(i,1)
                     tm(i,1) = tm(i,1) / vm

                  end if

               end do

               do i = 1, maxnt + maxpt

                  if( tm(i,1) .gt. 0.0d0 ) goto 200

               end do

                  goto 260

  200          im = i

               do i = maxnt + maxpt, im + 1, -1

                  if( tm(i,1) .gt. 0.0d0 ) goto 210

               end do

               jm = im
  210          jm = i

               im = max( 1, im - 2 )
               jm = jm + 2

*-----------------------------------------------------------------------

                  write(iot,'(/"#",78("-"))')

               inum = inum + 1

               if( inum .eq. 1 ) then
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

               if( rtfac(m) .lt. 0.d0 )
     &           write(iot,
     &           '( "#restart maximum for normalization:",1p1e14.7)')
     &           facmax(m)

               if( nn .eq. 0 ) then

                  write(iot,'("#   no. =",i3,3x,
     &            "ix =",i3,3x,"iy =",i3,3x,"iz =",i3,3x,
     &            "Z = all"/
     &            "#   x = (",1p1e13.4,"  -",1p1e13.4,"  )"/
     &            "#   y = (",1p1e13.4,"  -",1p1e13.4,"  )"/
     &            "#   z = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        inum, jx, jy, jz,
     &                        xm(jx), xm(jx+1),
     &                        ym(jy), ym(jy+1),
     &                        zm(jz), zm(jz+1)

               else

                  write(iot,'("#   no. =",i3,3x,
     &            "ix =",i3,3x,"iy =",i3,3x,"iz =",i3,3x,
     &            "Z = ",i3," : ",a3/
     &            "#   x = (",1p1e13.4,"  -",1p1e13.4,"  )"/
     &            "#   y = (",1p1e13.4,"  -",1p1e13.4,"  )"/
     &            "#   z = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        inum, jx, jy, jz, iz, elmnt(iz),
     &                        xm(jx), xm(jx+1),
     &                        ym(jy), ym(jy+1),
     &                        zm(jz), zm(jz+1)

               end if

               write(iot,'("# im jm = ",i3,1x,i3)') im, jm

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: Mass")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: Number ",a15)') hsunit(itunt(m))

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

                  write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')

               if( itanl(m) .gt. 0 ) then

                  write(iot,'( "p: ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               end if

               if( itsans(m) .gt. 0 ) then

                  call write_sangel(iot,m,0)

               end if

c S.H. added IF statement below for epsout=2 (2017.2.18)
              if ( iteps(m) .ne. 2 ) then

               write(iot,'( "h: x-0.5    y,hl0       n")')

              else

               write(iot,'( "h: x-0.5    y1,hl0      ",
     &                      "ny2 dy1=[y1*y2]")')

              end if

               write(iot,'( "#  mass     number      r.err")')

*-----------------------------------------------------------------------

               if( im .eq. 1 ) then

                  write(iot,'(3x,f4.1,2x,1pe13.4,0pf8.4)')
     &                                   0.5, 0.0, 0.0

               end if

                  seka = 0.0
                  sera = 0.0

               do i = im, jm

                  vn  = vm * tm(i,1)
                  seka = seka + vn
                  sera = sera + ( vn * tm(i,2) )**2

                  write(iot,'(3x,i4,2x,1pe13.4,0pf8.4)')
     &                           i, tm(i,1), tm(i,2)

               end do

                  if( seka .gt. 0.0 ) then

                     sera = sqrt( sera ) / seka
                     seka = seka / vm

                  end if

               write(iot,'(/"#   sum  ",1pe13.4,0pf8.4)')
     &               seka, sera

               if( itunt(m) .eq. 2 ) then

                  write(iot,'(
     &            "#   vol sum = ",1p1e13.4,"  [cm^3]")') vm

               end if

               if( nn .eq. 0 ) then

                  write(iot,'(/a1,"no. =",i3,
     &                            ",   (ix,iy,iz)  =  (",i3,
     &                            ","i3,",",i3,")",
     &                            ",    Z  =  all", a1)')
     &                        cha, inum, jx, jy, jz, cha

               else

                  write(iot,'(/a1,"no. =",i3,
     &                            ",   (ix,iy,iz)  =  (",i3,
     &                            ","i3,",",i3,")",
     &                            ",    Z  =  ",i3, a1)')
     &                        cha, inum, jx, jy, jz, iz, cha

               end if

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

               if( itunt(m) .eq. 1 ) then

                  if( nn .eq. 0 ) then

                     write(iot,'("wt: s(0.7)",/
     &                           a1,"vspace{-3}"/
     &                           "  mass distribution"/
     &                           "  in xyz mesh"/
     &                           "  xmin  &=&",1pe13.4," [cm]"/
     &                           "  ymax  &=&",1pe13.4," [cm]"/
     &                           "  ymin  &=&",1pe13.4," [cm]"/
     &                           "  ymax  &=&",1pe13.4," [cm]"/
     &                           "  zmin  &=&",1pe13.4," [cm]"/
     &                           "  zmax  &=&",1pe13.4," [cm]"/
     &                           "      Z &=&   all"/
     &                           "e:")')
     &                           yen,
     &                           xm(jx), xm(jx+1),
     &                           ym(jy), ym(jy+1),
     &                           zm(jz), zm(jz+1)

                  else

                     write(iot,'("wt: s(0.7)",/
     &                           a1,"vspace{-3}"/
     &                           "  isotope distribution"/
     &                           "  in xyz mesh"/
     &                           "  xmin  &=&",1pe13.4," [cm]"/
     &                           "  xmax  &=&",1pe13.4," [cm]"/
     &                           "  ymin  &=&",1pe13.4," [cm]"/
     &                           "  ymax  &=&",1pe13.4," [cm]"/
     &                           "  zmin  &=&",1pe13.4," [cm]"/
     &                           "  zmax  &=&",1pe13.4," [cm]"/
     &                           "      Z &=& ",i3,"   :  ",a3/
     &                           "e:")')
     &                           yen,
     &                           xm(jx), xm(jx+1),
     &                           ym(jy), ym(jy+1),
     &                           zm(jz), zm(jz+1),
     &                           iz, elmnt(iz)

                  end if

               else if( itunt(m) .eq. 2 ) then

                  if( nn .eq. 0 ) then

                     write(iot,'("wt: s(0.7)",/
     &                           a1,"vspace{-3}"/
     &                           "  mass distribution"/
     &                           "  in xyz mesh"/
     &                           "  xmin  &=&",1pe13.4," [cm]"/
     &                           "  xmax  &=&",1pe13.4," [cm]"/
     &                           "  ymin  &=&",1pe13.4," [cm]"/
     &                           "  ymax  &=&",1pe13.4," [cm]"/
     &                           "  zmin  &=&",1pe13.4," [cm]"/
     &                           "  zmax  &=&",1pe13.4," [cm]"/
     &                           "    vol &=&",1pe13.4," [cm^3]"/
     &                           "      Z &=&   all"/
     &                           "e:")')
     &                           yen,
     &                           xm(jx), xm(jx+1),
     &                           ym(jy), ym(jy+1),
     &                           zm(jz), zm(jz+1),
     &                           vl(jx,jy,jz)

                  else

                     write(iot,'("wt: s(0.7)",/
     &                           a1,"vspace{-3}"/
     &                           "  isotope distribution"/
     &                           "  in xyz mesh"/
     &                           "  xmin  &=&",1pe13.4," [cm]"/
     &                           "  xmax  &=&",1pe13.4," [cm]"/
     &                           "  ymin  &=&",1pe13.4," [cm]"/
     &                           "  ymax  &=&",1pe13.4," [cm]"/
     &                           "  zmin  &=&",1pe13.4," [cm]"/
     &                           "  zmax  &=&",1pe13.4," [cm]"/
     &                           "    vol &=&",1pe13.4," [cm^3]"/
     &                           "      Z &=& ",i3,"   :  ",a3/
     &                           "e:")')
     &                           yen,
     &                           xm(jx), xm(jx+1),
     &                           ym(jy), ym(jy+1),
     &                           zm(jz), zm(jz+1),
     &                           vl(jx,jy,jz),
     &                           iz, elmnt(iz)

                  end if

               end if

*-----------------------------------------------------------------------

  260       continue

            end do
            end do
            end do

*-----------------------------------------------------------------------
*        charge axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 7 ) then

               inum = 0

            do 390 jx = 1, nx
            do 380 jy = 1, ny
            do 370 jz = 1, nz

                  do i = 1, maxpt

                     tm(i,1) = 0.0d+0
                     tm(i,2) = 0.0d+0

                  end do

                  if( itunt(m) .eq. 1 ) then

                     vm = 1.0

                  else if( itunt(m) .eq. 2 ) then

                     vm = vl(jx,jy,jz)

                  end if

*-----------------------------------------------------------------------

                  do kz = 1, mz
                  do kn = 1, mn
                  do il = 0, mm ! frtati 2022/02/18

                        iz  = ikzz( kz, kn )
                        in  = iknn( kz, kn )

                        vn  = vm * tr(jx,jy,jz,iz,in,il,1)
                        tm(iz,1) = tm(iz,1) + vn
                        tm(iz,2) = tm(iz,2)
     &                           + ( vn * tr(jx,jy,jz,iz,in,il,2) )**2

                  end do
                  end do
                  end do

*-----------------------------------------------------------------------

               do i = 1, maxpt

                  if( tm(i,1) .gt. 0.0d0 ) then

                     tm(i,2) = sqrt( tm(i,2) ) / tm(i,1)
                     tm(i,1) = tm(i,1) / vm

                  end if

               end do

               do i = 1, maxpt

                  if( tm(i,1) .gt. 0.0d0 ) goto 300

               end do

               goto 370

  300          im = i

               do i = maxpt, im + 1, -1

                  if( tm(i,1) .gt. 0.0d0 ) goto 310

               end do

               jm = im
  310          jm = i

               im = max( 1, im - 2 )
               jm = jm + 2

*-----------------------------------------------------------------------

                  write(iot,'(/"#",78("-"))')

               inum = inum + 1

               if( inum .eq. 1 ) then
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

               if( rtfac(m) .lt. 0.d0 )
     &           write(iot,
     &           '( "#restart maximum for normalization:",1p1e14.7)')
     &           facmax(m)

               write(iot,'("#   no. =",i3,3x,
     &         "ix =",i3,3x,"iy =",i3,3x,"iz =",i3,3x/
     &         "#   x = (",1p1e13.4,"  -",1p1e13.4,"  )"/
     &         "#   y = (",1p1e13.4,"  -",1p1e13.4,"  )"/
     &         "#   z = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     inum, jx, jy, jz,
     &                     xm(jx), xm(jx+1),
     &                     ym(jy), ym(jy+1),
     &                     zm(jz), zm(jz+1)

               write(iot,'("# im jm = ",i3,1x,i3)') im, jm

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: Charge")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: Number ",a15)') hsunit(itunt(m))

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

                  write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')

               if( itanl(m) .gt. 0 ) then

                  write(iot,'( "p: ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               end if

               if( itsans(m) .gt. 0 ) then

                  call write_sangel(iot,m,0)

               end if

c S.H. added IF statement below for epsout=2 (2017.2.18)
              if ( iteps(m) .ne. 2 ) then

               write(iot,'( "h: x-0.5    y,hl0       n")')

              else

               write(iot,'( "h: x-0.5    y1,hl0      ",
     &                      "ny2 dy1=[y1*y2]")')

              end if

               write(iot,'( "#  charge   number      r.err")')

*-----------------------------------------------------------------------

               if( im .eq. 1 ) then

                  write(iot,'(3x,f4.1,2x,1pe13.4,0pf8.4)')
     &                                   0.5, 0.0, 0.0

               end if

                  seka = 0.0
                  sera = 0.0

               do i = im, jm

                  vn  = vm * tm(i,1)
                  seka = seka + vn
                  sera = sera + ( vn * tm(i,2) )**2

                  write(iot,'(3x,i4,2x,1pe13.4,0pf8.4)')
     &                           i, tm(i,1), tm(i,2)

               end do

                  if( seka .gt. 0.0 ) then

                     sera = sqrt( sera ) / seka
                     seka = seka / vm

                  end if

               write(iot,'(/"#   sum  ",1pe13.4,0pf8.4)')
     &               seka, sera

               if( itunt(m) .eq. 2 ) then

                  write(iot,'(
     &            "#   vol sum = ",1p1e13.4,"  [cm^3]")') vm

               end if

                  write(iot,'(/a1,"no. =",i3,
     &                            ",   (ix,iy,iz)  =  (",i3,
     &                            ","i3,",",i3,")",
     &                             a1)')
     &                        cha, inum, jx, jy, jz, cha

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

               if( itunt(m) .eq. 1 ) then

                     write(iot,'("wt: s(0.7)",/
     &                           a1,"vspace{-3}"/
     &                           "  charge distribution"/
     &                           "  in xyz mesh"/
     &                           "  xmin  &=&",1pe13.4," [cm]"/
     &                           "  xmax  &=&",1pe13.4," [cm]"/
     &                           "  ymin  &=&",1pe13.4," [cm]"/
     &                           "  ymax  &=&",1pe13.4," [cm]"/
     &                           "  zmin  &=&",1pe13.4," [cm]"/
     &                           "  zmax  &=&",1pe13.4," [cm]"/
     &                           "e:")')
     &                           yen,
     &                           xm(jx), xm(jx+1),
     &                           ym(jy), ym(jy+1),
     &                           zm(jz), zm(jz+1)

               else if( itunt(m) .eq. 2 ) then

                     write(iot,'("wt: s(0.7)",/
     &                           a1,"vspace{-3}"/
     &                           "  charge distribution"/
     &                           "  in xyz mesh"/
     &                           "  xmin  &=&",1pe13.4," [cm]"/
     &                           "  xmax  &=&",1pe13.4," [cm]"/
     &                           "  ymin  &=&",1pe13.4," [cm]"/
     &                           "  ymax  &=&",1pe13.4," [cm]"/
     &                           "  zmin  &=&",1pe13.4," [cm]"/
     &                           "  zmax  &=&",1pe13.4," [cm]"/
     &                           "    vol &=&",1pe13.4," [cm^3]"/
     &                           "e:")')
     &                           yen,
     &                           xm(jx), xm(jx+1),
     &                           ym(jy), ym(jy+1),
     &                           zm(jz), zm(jz+1),
     &                           vl(jx,jy,jz)

               end if

  370       continue
  380       continue
  390       continue

*-----------------------------------------------------------------------
*        chart axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 8 ) then

               inum = 0

            do 490 jx = 1, nx
            do 480 jy = 1, ny
            do 470 jz = 1, nz
            do 460 il = 0, 2 ! loop for isomeric level

*-----------------------------------------------------------------------

                  itmax  = 0
                  icmax  = 0
                  inmax  = 0

               do iz = 1, maxpt
               do in = 1, maxnt

                  if( itnzn(m).eq.0 ) then
                    t0_tr = tr(jx,jy,jz,iz,in,il,1)
                  else
                    t0_tr = tr(jx,jy,jz,igetiznmp(iz,in,il,m),1,0,1)
                  end if
                  if( t0_tr .gt. 0.d0 ) then

                     if( iz + in .gt. itmax ) itmax = iz + in
                     if( iz      .gt. icmax ) icmax = iz
                     if(      in .gt. inmax ) inmax =      in

                  end if

               end do
               end do

                     dxmax = dble(inmax+2)
                     dymax = dble(icmax+2)
                     dform = dymax / dxmax

                     inmag = inmax

                  if( dform .gt. 0.8 ) then

                     inmag = nint( dymax / 0.8 )
                     dxmax = dble( nint( dymax / 0.8 ) )
                     dform = dymax / dxmax

                  end if

               if( itmax .eq. 0 ) goto 470

*-----------------------------------------------------------------------

                  write(iot,'(/"#",78("-"))')

               inum = inum + 1

               if( inum .eq. 1 .and. il .eq. 0 ) then
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

               if( rtfac(m) .lt. 0.d0 )
     &           write(iot,
     &           '( "#restart maximum for normalization:",1p1e14.7)')
     &           facmax(m)

               write(iot,'("#   no. =",i3,3x,"ix =",i3,3x,
     &         "iy =",i3,3x,"iz  =",i3,3x,"il  =",i3)')
     &                     inum, jx, jy, jz, il

               write(iot,'("# icmax inmax = ",i3,1x,i3)') icmax, inmax

*-----------------------------------------------------------------------

            if( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) then

                  write(iot,'(/a1,"no. =",i3,
     &                            ",   (ix,iy,iz)  =  (",i3,
     &                            ","i3,",",i3,")",
     &                            a1)')
     &                        cha, inum, jx, jy, jz, cha

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

             if(il .eq. 0) then    ! specify isomeric level at the bottom of figure
               write(iot,'("msdc: {",a1,"huge Ground state}")') yen
             elseif(il .eq. 1) then
               write(iot,'("msdc: {",a1,"huge 1st isomer}")') yen
             else
               write(iot,'("msdc: {",a1,"huge 2nd isomer}")') yen
             endif

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: N Neutron Number")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: Z Proton Number")')

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

*-----------------------------------------------------------------------

               if( inum .eq. 1 ) then

                  form  = dform
                  xfac  = 1.1
                  afac  = 0.6
                  izlog = 1
                  inocm = 1
                  inolg = 1

                  if( itanl(m) .gt. 0 )
     &            call anset(itang(m),itanl(m),
     &                       form,xfac,afac,cmin,cmax,izlog,inocm,inolg)

                  if( form .le. 1.0 ) then

                     scal = form**0.35
                     xfac = xfac / form**0.5
                     xorg = 0.0
                     yorg = min( 20.d0, ( 1.0 / form - 1.0 ) / 2.5 )

                  else

                     scal = 1.0 / form**0.41
                     xfac = xfac / form**0.5
                     xorg = min( 20.d0, ( form - 1.0 ) / 1.5)
                     yorg = 0.0

                  end if

               end if

                  write(iot,'( "set: c1[",f6.3,
     &                            "] c2[",f6.3,
     &                            "] c5[",f6.3,
     &                            "] c6[",f6.3,
     &                            "] c7[",f6.3,
     &                            "] c8[",f6.3,"]")')
     &                    form, xfac, afac, scal, xorg, yorg
                  write(iot,'( "p: h2fs form[c1] xfac[c2]",
     &                 " afac[c5] scal[c6] xorg[c7] yorg[c8] nosp")')

               if( ( ( ittwo(m) .ge. 2 .and. ittwo(m) .le. 3 ) .or.
     &               ( ittwo(m) .ge. 6 .and. ittwo(m) .le. 7 ) ) .and.
     &             cmin .gt. 0.0 .and. cmax .gt. cmin ) then

                if (ioe .eq. 1 ) then
                  write(iot,'( "set: c3[",1p1e15.6,"] c4[",
     &                  1p1e15.6,"]")') cmin, cmax
                else
                  write(iot,'( "set: c3[1.0e-4] c4[1.0]")')
                end if
                  write(iot,'( "p: cmin[c3] cmax[c4]")')
                  write(iot,'( "p: dmin(1e-31)")')

                  if( ittwo(m) .eq. 3 .or. ittwo(m) .eq. 7 )
     &                 write(iot,'( "p: zlog")')

               end if

               write(iot,'( "p: xmin(0) xmax(",f5.1,")")') dxmax
               write(iot,'( "p: ymin(0) ymax(",f5.1,")")') dymax

               if( itanl(m) .gt. 0 .and. ioe .eq. 1 ) then

                  write(iot,'( "p: ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               else if( itanl(m) .gt. 0 .and. ioe .ne. 1 ) then

                   call terrang(iot,m)

               end if

               if( itsans(m) .gt. 0 ) then

                  call write_sangel(iot,m,0)

               end if

            end if

*-----------------------------------------------------------------------

               write(iot,'("#  ny = ",i3,"   nx = ",i3)')
     &                      icmax+1, inmax+2

            if( ittwo(m) .ne. 4 ) then

               write(iot,'( "# ( ( data(x,y), x = 1, nx ),",
     &                      " y = ny, 1, -1 )")')

            end if

*-----------------------------------------------------------------------

            if( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) then

               write(iot,'(/"hc: y = ",i3," to 1 by -1 ;",
     &                         " x = 1 to ",i3," by 1 ;")')
     &                        icmax+2, inmax+2

               do i = icmax+2, 1, -1
                if( itnzn(m).eq.0 ) then ! frtati 2022/02/18
                  write(iot,'(1p10e11.3)')
     &                 ( tr(jx,jy,jz,i,l,il,ioe), l = 1, inmax+2 )
                else
                  write(iot,'(1p10e11.3)')
     &          ( tr(jx,jy,jz,igetiznmp(i,l,il,m),1,0,ioe),
     &            l = 1, inmax+2 )
                end if

               end do

*-----------------------------------------------------------------------

            else if( ittwo(m) .eq. 4 ) then

               write(iot,'(/"#   Z    N  Mass  ",
     &                      "  number    r.err")')

               do i = 1, icmax
               do j = 1, inmax

                if( itnzn(m).eq.0 ) then ! frtati 2022/02/18
                  if( tr(jx,jy,jz,i,j,il,1) .gt. 0.d0 ) then

                     write(iot,'(3i5,1pe13.4,0pf8.4)')
     &               i, j, i+j,
     &               tr(jx,jy,jz,i,j,il,1), tr(jx,jy,jz,i,j,il,2)

                  end if
                else
                  if( tr(jx,jy,jz,igetiznmp(i,j,il,m),1,0,1) .gt.
     &                0.d0 ) then
                     write(iot,'(3i5,1pe13.4,0pf8.4)')
     &               i, j, i+j,
     &               tr(jx,jy,jz,igetiznmp(i,j,il,m),1,0,1),
     &               tr(jx,jy,jz,igetiznmp(i,j,il,m),1,0,2)
                  end if
                end if

               end do
               end do

*-----------------------------------------------------------------------

            else if( ittwo(m) .eq. 5 ) then

               write(iot,'(/"#   Z = 1 to ",i3/
     &                      "#   N = 1 to ",i3)')
     &                        icmax+2, inmax+2
                  write(iot,'(/4x,a3,4x,1p1000e11.3)')
     &            'Z/N', ( dble( in ), in = 1, inmax+2 )

               do i = icmax+2, 1, -1

                if( itnzn(m).eq.0 ) then ! frtati 2022/02/18
                  write(iot,'(1p1000e11.3)')
     &            dble( i ),
     &            ( tr(jx,jy,jz,i,l,il,ioe), l = 1, inmax+2 )
                else
                  write(iot,'(1p1000e11.3)')
     &            dble( i ),
     &            ( tr(jx,jy,jz,igetiznmp(i,l,il,m),1,0,ioe),
     &              l = 1,inmax+2 )
                end if

               end do

            end if

*-----------------------------------------------------------------------

            if( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) then

               if( itout(m) .ne. 0 ) then

                  call wmgcstb(iot,icmax,inmag)

               end if

            end if

*-----------------------------------------------------------------------

      if( inolg .eq. 1 .and. ( ittwo(m) .eq. 3 .or. ittwo(m) .eq. 7 )
     &    .and. cmin .gt. 0.0 .and. cmax .gt. cmin .and.
     &    igsh .eq. 0 ) then

      write(iot,'(/"#",78("-"))')
      write(iot,'(
     &"z: xorg(1.03)"/
     &"p: xfac[c2*0.03] form[c1/0.03] notc noxt noyt dmin(1)"/
     &"p: ymin(0) ymax(1) xmin(0) xmax(1) cmin(1) cmax(100) nosp"/
     &"hc: y= 0.005 to 0.995 by 0.01 ; x= 0.5 to 0.5 by 1 ;"/
     &" 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20"/
     &"21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40"/
     &"41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60"/
     &"61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80"/
     &"81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100"/
     &"z: xorg(0.0)"/
     &"p: xfac[c2*0.03] form[c1/0.03] ",
     &"nosp afac[c5*0.7] ymin[c3] ymax[c4]"/
     &"p: noxt noxn ytxt(-1) ynum(-1) xltd(-1) xstd(-1) itic(1) ylog"
     &)')

      if ( ioe .eq. 2 ) then

          write(iot,'("y: Relative Error")')

      else if( itazl(m) .eq. 0 ) then

         write(iot,'("y: Number ",a15)') hsunit(itunt(m))

      else

         write(iot,'("y: ",200a1)') (itazt(m)(i:i),i=1,itazl(m))

      end if

      if( izlog .eq. 0 ) write(iot,'( "p: ylin")')

      write(iot,'( "#",78("-"))')

      end if

*-----------------------------------------------------------------------

            if( inocm .eq. 1 ) then

               if( itunt(m) .eq. 1 ) then

                     write(iot,'(/"wt: s(0.7)",/
     &                           a1,"vspace{-3}"/
     &                           "  yield distribution"/
     &                           "  in xyz mesh"/
     &                           "  xmin  &=&",1pe13.4," [cm]"/
     &                           "  xmax  &=&",1pe13.4," [cm]"/
     &                           "  ymin  &=&",1pe13.4," [cm]"/
     &                           "  ymax  &=&",1pe13.4," [cm]"/
     &                           "  zmin  &=&",1pe13.4," [cm]"/
     &                           "  zmax  &=&",1pe13.4," [cm]"/
     &                           "e:")')
     &                           yen,
     &                           xm(jx), xm(jx+1),
     &                           ym(jy), ym(jy+1),
     &                           zm(jz), zm(jz+1)

               else if( itunt(m) .eq. 2 ) then

                     write(iot,'(/"wt: s(0.7)",/
     &                           a1,"vspace{-3}"/
     &                           "  yield distribution"/
     &                           "  in xyz mesh"/
     &                           "  xmin  &=&",1pe13.4," [cm]"/
     &                           "  xmax  &=&",1pe13.4," [cm]"/
     &                           "  ymin  &=&",1pe13.4," [cm]"/
     &                           "  ymax  &=&",1pe13.4," [cm]"/
     &                           "  zmin  &=&",1pe13.4," [cm]"/
     &                           "  zmax  &=&",1pe13.4," [cm]"/
     &                           "    vol &=&",1pe13.4," [cm^3]"/
     &                           "e:")')
     &                           yen,
     &                           xm(jx), xm(jx+1),
     &                           ym(jy), ym(jy+1),
     &                           zm(jz), zm(jz+1),
     &                           vl(jx,jy,jz)

               end if

            end if

*-----------------------------------------------------------------------

      if( inolg .eq. 1 .and. ( ittwo(m) .eq. 3 .or. ittwo(m) .eq. 7 )
     &    .and. cmin .gt. 0.0 .and. cmax .gt. cmin .and.
     &    igsh .eq. 0 ) then

         write(iot,'(/"#",78("-"))')
         write(iot,'("z: xorg[-1.03/0.05]"/
     &               "p: form[c1] xfac[c2] afac[c5] nosp notf")')

               write(iot,'( "p: xmin(0) xmax(",f5.1,")")') dxmax
               write(iot,'( "p: ymin(0) ymax(",f5.1,")")') dymax

      end if

*-----------------------------------------------------------------------

  460       continue
  470       continue
  480       continue
  490       continue

*-----------------------------------------------------------------------
*        xy axis (matrix)
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 9 ) then

               if( ittwo(m) .eq. 1 ) then

                  dc2 = 'h2: '

               else if( ittwo(m) .eq. 2 ) then

                  dc2 = 'hd: '

               else if( ittwo(m) .eq. 3 ) then

                  dc2 = 'hc: '

               else if( ittwo(m) .eq. 6 ) then

                  dc2 = 'hd2:'

               else if( ittwo(m) .eq. 7 ) then

                  dc2 = 'hc2:'

               end if

*-----------------------------------------------------------------------

               inum = 0

                  if( nn .eq. 0 .or. igsh .ne. 0  ) then

                     nc = 1

                  else

                     nc = nn

                  end if

*-----------------------------------------------------------------------
*              max and min
*-----------------------------------------------------------------------

                        fmax = 0.0
                        fmin = 1.e+33

                  do jz = 1, nz
                  do ic = 1, nc

                     if( nn .gt. 0 ) then

                        iz = nt(ic) / 1000

                        ia = nt(ic) - iz * 1000

                        if( ia .gt. 0 ) in = ia - iz

                     end if

                     do jx = 1, nx
                     do jy = 1, ny

                           sek = 0.0

                           if( itunt(m) .eq. 1 ) then

                              vm = 1.0

                           else if( itunt(m) .eq. 2 ) then

                              vm = vl(jx,jy,jz)

                           end if

                        if( nn .eq. 0 ) then

                           do kn = 1, mn
                           do kz = 1, mz
                           do il = 0, mm ! frtati 2022/02/18

                              lz  = ikzz( kz, kn )
                              ln  = iknn( kz, kn )

                              vn  = vm * tr(jx,jy,jz,lz,ln,il,1)
                              sek = sek + vn

                           end do
                           end do
                           end do

                        else if( ia .eq. 0 ) then

                           do kn = 1, mn
                           do il = 0, mm ! frtati 2022/02/18

                              lz  = ikzz( iz, kn )
                              ln  = iknn( iz, kn )

                              vn  = vm * tr(jx,jy,jz,lz,ln,il,1)
                              sek = sek + vn

                           end do
                           end do

                        else

                              lz  = ikzz( iz, in )
                              ln  = iknn( iz, in )

                           do il = 0, mm ! frtati 2022/02/18
                              vn  = vm * tr(jx,jy,jz,lz,ln,il,1)
                              sek = sek + vn
                           end do

                        end if

                              fmsv = sek / vm

                              if( fmsv .gt. fmax ) fmax = fmsv
                              if( fmsv .gt. 0.0 .and.
     &                            fmsv .lt. fmin ) fmin = fmsv

                     end do
                     end do

                  end do
                  end do

*-----------------------------------------------------------------------

               allocate( fmval(nx,ny,nz,nc) )

            do jz = 1, nz

            do ic = 1, nc

               inum = inum + 1

*-----------------------------------------------------------------------

               if( nn .gt. 0 ) then

                  iz = nt(ic) / 1000

                  ia = nt(ic) - iz * 1000

                  if( ia .gt. 0 ) then

                     in = ia - iz

                     call chname(idum,ia,iz,chau)

                  end if

               end if

*-----------------------------------------------------------------------

                  write(iot,'(/"#",78("-"))')

               if( inum .eq. 1 ) then
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

               if( rtfac(m) .lt. 0.d0 )
     &           write(iot,
     &           '( "#restart maximum for normalization:",1p1e14.7)')
     &           facmax(m)

            if( nn .eq. 0 ) then

               write(iot,'("#   no. =",i3,3x,
     &         "iz  =",i3,3x,"Z = all")')
     &                     inum, jz

            else if( ia .eq. 0 ) then

               write(iot,'("#   no. =",i3,3x,
     &         "iz  =",i3,3x,"Z = ",i3," : ",a3)')
     &                     inum, jz, iz, elmnt(iz)

            else

               write(iot,'("#   no. =",i3,3x,
     &         "iz  =",i3,3x,a8)')
     &                     inum, jz, chau

            end if

*-----------------------------------------------------------------------

            if( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) then

               if( nn .eq. 0 ) then

                  write(iot,'(/a1,"no. =",i3,
     &                        ",    iz =",i3,
     &                        ",    Z  =  all", a1)')
     &                        cha, inum, jz, cha

               else if( ia .eq. 0 ) then

                  write(iot,'(/a1,"no. =",i3,
     &                        ",    iz =",i3,
     &                        ",    Z  =  ",i3, a1)')
     &                        cha, inum, jz, iz, cha

               else

                  write(iot,'(/a1,"no. =",i3,
     &                        ",    iz =",i3,
     &                        ",    ",a8, a1)')
     &                        cha, inum, jz, chau, cha

               end if

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: x [cm]")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: y [cm]")')

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

*-----------------------------------------------------------------------

               if( inum .eq. 1 ) then

                  form  = ( ym(ny+1) - ym(1) )
     &                  / ( xm(nx+1) - xm(1) )
                  xfac  = 0.9
                  afac  = 0.8
                  izlog = 1
                  inocm = 1
                  inolg = 1

                  if( itanl(m) .gt. 0 )
     &            call anset(itang(m),itanl(m),
     &                       form,xfac,afac,fmin,fmax,izlog,inocm,inolg)

                  if( form .le. 1.0 ) then

                     scal = form**0.35
                     xfac = xfac / form**0.5
                     xorg = 0.0
                     yorg = min( 20.d0, ( 1.0 / form - 1.0 ) / 2.5 )

                  else

                     scal = 1.0 / form**0.41
                     xfac = xfac / form**0.5
                     xorg = min( 20.d0, ( form - 1.0 ) / 1.5)
                     yorg = 0.0

                  end if

               end if

                  write(iot,'( "set: c1[",f6.3,
     &                            "] c2[",f6.3,
     &                            "] c5[",f6.3,
     &                            "] c6[",f6.3,
     &                            "] c7[",f6.3,
     &                            "] c8[",f6.3,"]")')
     &                    form, xfac, afac, scal, xorg, yorg
                  write(iot,'( "p: h2fs form[c1] xfac[c2]",
     &                 " afac[c5] scal[c6] xorg[c7] yorg[c8] nosp")')

               if( ( ( ittwo(m) .ge. 2 .and. ittwo(m) .le. 3 ) .or.
     &               ( ittwo(m) .ge. 6 .and. ittwo(m) .le. 7 ) ) .and.
     &             fmin .gt. 0.0 .and. fmax .gt. fmin ) then

                if (ioe .eq. 1 ) then
                  write(iot,'( "set: c3[",1p1e15.6,"] c4[",
     &                  1p1e15.6,"]")') fmin, fmax
                else
                  write(iot,'( "set: c3[1.0e-4] c4[1.0]")')
                end if
                  write(iot,'( "p: cmin[c3] cmax[c4]")')
                  write(iot,'( "p: dmin(1e-31)")')

                  if( ittwo(m) .eq. 3 .or. ittwo(m) .eq. 7 )
     &                 write(iot,'( "p: zlog")')

               end if

               xmin = xm(1)
               xmax = xm(nx+1)
               ymin = ym(1)
               ymax = ym(ny+1)

               write(iot,'( "p: xmin(",1p1g14.7,") xmax(",1p1g14.7,
     &                       ")")') xmin, xmax

               write(iot,'( "p: ymin(",1p1g14.7,") ymax(",1p1g14.7,
     &                       ")")') ymin, ymax

               if( itanl(m) .gt. 0 .and. ioe .eq. 1 ) then

                  write(iot,'( "p: ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               else if( itanl(m) .gt. 0 .and. ioe .ne. 1 ) then

                   call terrang(iot,m)

               end if

               if( itsans(m) .gt. 0 ) then

                  call write_sangel(iot,m,0)

               end if

            end if

*-----------------------------------------------------------------------

               do jx = 1, nx
               do jy = 1, ny

                     sek = 0.0

                     if( itunt(m) .eq. 1 ) then

                        vm = 1.0

                     else if( itunt(m) .eq. 2 ) then

                        vm = vl(jx,jy,jz)

                     end if

*-----------------------------------------------------------------------

                  if( nn .eq. 0 ) then

                     do kn = 1, mn
                     do kz = 1, mz
                     do il = 0, mm ! frtati 2022/02/18

                        lz  = ikzz( kz, kn )
                        ln  = iknn( kz, kn )

                        vn  = vm * tr(jx,jy,jz,lz,ln,il,ioe)
                        sek = sek + vn

                     end do
                     end do
                     end do

                  else if( ia .eq. 0 ) then

                     do kn = 1, mn
                     do il = 0, mm ! frtati 2022/02/18

                        lz  = ikzz( iz, kn )
                        ln  = iknn( iz, kn )

                        vn  = vm * tr(jx,jy,jz,lz,ln,il,ioe)
                        sek = sek + vn

                     end do
                     end do

                  else

                        lz  = ikzz( iz, in )
                        ln  = iknn( iz, in )

                     do il = 0, mm ! frtati 2022/02/18
                        vn  = vm * tr(jx,jy,jz,lz,ln,il,ioe)
                        sek = sek + vn
                     end do

                  end if


                        fm(jx,jy) = sek / vm

               end do
               end do

*-----------------------------------------------------------------------

               write(iot,'("#  ny = ",i3,"   nx = ",i3)')
     &                      ny, nx

            if( ittwo(m) .ne. 4 ) then

               write(iot,'( "# ( ( data(x,y), x = 1, nx ),",
     &                      " y = ny, 1, -1 )")')

            end if

*-----------------------------------------------------------------------

            if( ( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) .and.
     &            igsh .eq. 0 ) then

               write(iot,'(/a4," y = ",1p1g14.7," to ",
     &         1p1g14.7," by ",1p1g14.7," ; x = ",
     &         1p1g14.7," to ",1p1g14.7," by ",
     &         1p1g14.7," ;")') dc2,
     &         ym(ny) + rtydl(m)/2.0, ym(1) + rtydl(m)/2.0, rtydl(m),
     &         xm(1) + rtxdl(m)/2.0, xm(nx) + rtxdl(m)/2.0, rtxdl(m)

               write(iot,'(1p10e11.3)')
     &         ( ( fm(jx,jy), jx = 1, nx ), jy = ny, 1, -1 )

*-----------------------------------------------------------------------

            else if( ittwo(m) .eq. 4 ) then

               write(iot,'(/"# x          y        ",
     &                      "  number")')

               do jy = 1, ny
               do jx = 1, nx

                  write(iot,'(1p10e11.3)')
     &               xm(jx)  + rtxdl(m)/2.0,
     &               ym(jy)  + rtydl(m)/2.0,
     &               fm(jx,jy)

               end do
               end do

*-----------------------------------------------------------------------

            else if( ittwo(m) .eq. 5 ) then

               write(iot,'("#   x = ",1p1g14.7," to ",
     &         1p1g14.7," by ",1p1g14.7/
     &                     "#   y = ",1p1g14.7," to ",
     &         1p1g14.7," by ",1p1g14.7)')
     &         xm(1) + rtxdl(m)/2.0, xm(nx) + rtxdl(m)/2.0, rtxdl(m),
     &         ym(1) + rtydl(m)/2.0, ym(ny) + rtydl(m)/2.0, rtydl(m)

                  write(iot,'(/4x,a3,4x,1p1000e11.3)')
     &            'y/x',( xm(jx) + rtxdl(m)/2.0, jx = 1, nx )

               do jy = ny, 1, -1

                  write(iot,'(1p1000e11.3)')
     &            ym(jy) + rtydl(m)/2.0,
     &            ( fm(jx,jy), jx = 1, nx )

               end do

            end if

*-----------------------------------------------------------------------
*        gshow
*-----------------------------------------------------------------------

         if( itgsh(m) .ne. 0 .and.
     &     ( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) ) then

               write(iot,'(/"#",78("-"))')
               write(iot,'( "# gshow")')
               write(iot,'( "#",78("-"))')
               write(iot,'(/"p: legs[c5*0.875]")')

               zval = ( zm(jz) + zm(jz+1) ) / 2.0d0
               none = 1
               iaxs = 1
               iuni = itgsh(m) * 2 - 1 + igsh
               ires = itres(m)
               widt = rtwid(m)

               call gshow(0,iot,iaxs,iuni,widt,ires,
     &                    inum,igsh,igser,gfnam,igfmn,
     &                    nx+1,ny+1,none,xm,ym,zval,ixyz(1),
     &                    1,1,krr,vll,itmtr(m,4),itglt(m))

         end if

*-----------------------------------------------------------------------

      if( inolg .eq. 1 .and. ( ittwo(m) .eq. 3 .or. ittwo(m) .eq. 7 )
     &    .and. fmin .gt. 0.0 .and. fmax .gt. fmin .and.
     &    igsh .eq. 0 ) then

      write(iot,'(/"#",78("-"))')
      write(iot,'(
     &"z: xorg(1.03)"/
     &"p: xfac[c2*0.05] form[c1/0.05] notc noxt noyt dmin(1)"/
     &"p: ymin(0) ymax(1) xmin(0) xmax(1) cmin(1) cmax(100) nosp"/
     &"hc: y= 0.005 to 0.995 by 0.01 ; x= 0.5 to 0.5 by 1 ;"/
     &" 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20"/
     &"21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40"/
     &"41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60"/
     &"61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80"/
     &"81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100"/
     &"z: xorg(0.0)"/
     &"p: xfac[c2*0.05] form[c1/0.05] ",
     &"nosp afac[c5*0.625] ymin[c3] ymax[c4]"/
     &"p: noxt noxn ytxt(-1) ynum(-1) xltd(-1) xstd(-1) itic(1) ylog"
     &)')

      if ( ioe .eq. 2 ) then

          write(iot,'("y: Relative Error")')

      else if( itazl(m) .eq. 0 ) then

         write(iot,'("y: Number ",a15)') hsunit(itunt(m))

      else

         write(iot,'("y: ",200a1)') (itazt(m)(i:i),i=1,itazl(m))

      end if

      if( izlog .eq. 0 ) write(iot,'( "p: ylin")')

      write(iot,'( "#",78("-"))')

      end if

*-----------------------------------------------------------------------

            if( inocm .eq. 1 ) then

               if( nn .eq. 0 ) then

                  write(iot,'(/"wt: s[c5]",/
     &                        a1,"vspace{-3}"/
     &                        "  zmin  &=&",1pe13.4," [cm]"/
     &                        "  zmax  &=&",1pe13.4," [cm]"/
     &                        "      Z &=&   all"/
     &                        "e:")')
     &                        yen,
     &                        zm(jz), zm(jz+1)

               else if( ia .eq. 0 ) then

                  write(iot,'(/"wt: s[c5]",/
     &                        a1,"vspace{-3}"/
     &                        "  zmin  &=&",1pe13.4," [cm]"/
     &                        "  zmax  &=&",1pe13.4," [cm]"/
     &                        "      Z &=& ",i3,"   :  ",a3/
     &                        "e:")')
     &                        yen,
     &                        zm(jz), zm(jz+1), iz, elmnt(iz)

               else

                  write(iot,'(/"wt: s[c5]",/
     &                        a1,"vspace{-3}"/
     &                        "  zmin  &=&",1pe13.4," [cm]"/
     &                        "  zmax  &=&",1pe13.4," [cm]"/
     &                        "  nucl. &=&   ",a8/
     &                        "e:")')
     &                        yen,
     &                        zm(jz), zm(jz+1), chau

               end if

            end if

*-----------------------------------------------------------------------

      if( inolg .eq. 1 .and. ( ittwo(m) .eq. 3 .or. ittwo(m) .eq. 7 )
     &    .and. cmin .gt. 0.0 .and. cmax .gt. cmin .and.
     &    igsh .eq. 0 ) then

         write(iot,'(/"#",78("-"))')
         write(iot,'("z: xorg[-1.03/0.05]"/
     &               "p: form[c1] xfac[c2] afac[c5] nosp notf")')

               write(iot,'( "p: xmin(",1p1g14.7,") xmax(",1p1g14.7,
     &                       ")")') xmin, xmax

               write(iot,'( "p: ymin(",1p1g14.7,") ymax(",1p1g14.7,
     &                       ")")') ymin, ymax

      end if

*-----------------------------------------------------------------------

               if ( itvtk(m) .ne. 0 ) then

                  do jy = 1, ny
                     do jx = 1, nx
                        fmval(jx,jy,jz,ic) = fm(jx,jy)
                     end do
                  end do

               end if

            end do
            end do

*-----------------------------------------------------------------------

            if ( itbmp(m) .ne. 0 ) then

               numIndex = 2
               allocate( bmpfIType(1:numIndex) )
               allocate( bmpfIndex(1:numIndex) )
               bmpfIType = (/ 'n', 'z' /)
               bmpfIndex = (/ nc, nz /)
               bmpWidth  = nx
               bmpHeight = ny

            end if

            if ( itvtk(m) .ne. 0 ) then

               call open_file(isunit_vtk_meta_default,
     &                 "", isunit_vtk_meta, ios, .true.)
               call open_file(isunit_vtk_default,
     &                 "", isunit_vtk, ios, .true.)

               iaxs = 1
               nparam = 1
               ntg = 1

               write(isunit_vtk_meta) ntg

               write(isunit_vtk_meta) iaxs
               write(isunit_vtk_meta) nparam
               write(isunit_vtk_meta) 'n'
               write(isunit_vtk_meta) nc

               write(isunit_vtk_meta) nx,ny,nz
               write(isunit_vtk_meta)
     &                 ( xm(ix), ix=1,nx+1 )
               write(isunit_vtk_meta)
     &                 ( ym(iy), iy=1,ny+1 )
               write(isunit_vtk_meta)
     &                 ( zm(iz), iz=1,nz+1 )

               do ic = 1, nc

                  write(isunit_vtk_meta) ic

                  write(isunit_vtk)
     &                    ( ( ( fmval(ix,iy,iz,ic),
     &                          ix=1,nx ),
     &                          iy=1,ny ),
     &                          iz=1,nz )
               end do

               iuni = itgsh(m) * 2 - 1 + igsh
               ires = itres(m)

               call open_file(isunit_vtk_rm_default,
     &                 "", isunit_vtk_rm, ios, .true.)

               if ( itgsh(m) .ne. 0 .and. ioe .eq. 1 ) then

                  call open_file(isunit_vtk_geom_default,
     &                    "", isunit_vtk_geom, ios, .true.)
                  call open_file(isunit_vtk_geom_meta_default,
     &                    "", isunit_vtk_geom_meta, ios, .true.)

               else
                  isunit_vtk_geom = 0
                  isunit_vtk_geom_meta = 0
               end if

               call vtk_set_cell_regmat(
     &                 isunit_vtk_rm,
     &                 isunit_vtk_geom, isunit_vtk_geom_meta,
     &                 nx, ny, nz, xm, ym, zm,
     &                 iaxs,iuni,ires, igser,
     &                 1,1,krr,vll,itmtr(m,4))

            end if

*-----------------------------------------------------------------------
*        yz axis (matrix)
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 10 ) then

               if( ittwo(m) .eq. 1 ) then

                  dc2 = 'h2: '

               else if( ittwo(m) .eq. 2 ) then

                  dc2 = 'hd: '

               else if( ittwo(m) .eq. 3 ) then

                  dc2 = 'hc: '

               else if( ittwo(m) .eq. 6 ) then

                  dc2 = 'hd2:'

               else if( ittwo(m) .eq. 7 ) then

                  dc2 = 'hc2:'

               end if

*-----------------------------------------------------------------------

               inum = 0

                  if( nn .eq. 0 .or. igsh .ne. 0  ) then

                     nc = 1

                  else

                     nc = nn

                  end if

*-----------------------------------------------------------------------
*              max and min
*-----------------------------------------------------------------------

                        fmax = 0.0
                        fmin = 1.e+33

                  do jx = 1, nx
                  do ic = 1, nc

                     if( nn .gt. 0 ) then

                        iz = nt(ic) / 1000

                        ia = nt(ic) - iz * 1000

                        if( ia .gt. 0 ) in = ia - iz

                     end if

                     do jz = 1, nz
                     do jy = 1, ny

                           sek = 0.0

                           if( itunt(m) .eq. 1 ) then

                              vm = 1.0

                           else if( itunt(m) .eq. 2 ) then

                              vm = vl(jx,jy,jz)

                           end if

                        if( nn .eq. 0 ) then

                           do kn = 1, mn
                           do kz = 1, mz
                           do il = 0, mm ! frtati 2022/02/18

                              lz  = ikzz( kz, kn )
                              ln  = iknn( kz, kn )

                              vn  = vm * tr(jx,jy,jz,lz,ln,il,1)
                              sek = sek + vn

                           end do
                           end do
                           end do

                        else if( ia .eq. 0 ) then

                           do kn = 1, mn
                           do il = 0, mm ! frtati 2022/02/18

                              lz  = ikzz( iz, kn )
                              ln  = iknn( iz, kn )

                              vn  = vm * tr(jx,jy,jz,lz,ln,il,1)
                              sek = sek + vn

                           end do
                           end do

                        else

                              lz  = ikzz( iz, in )
                              ln  = iknn( iz, in )

                           do il = 0, mm ! frtati 2022/02/18
                              vn  = vm * tr(jx,jy,jz,lz,ln,il,1)
                              sek = sek + vn
                           end do

                        end if

                              fmsv = sek / vm

                              if( fmsv .gt. fmax ) fmax = fmsv
                              if( fmsv .gt. 0.0 .and.
     &                            fmsv .lt. fmin ) fmin = fmsv

                     end do
                     end do

                  end do
                  end do

*-----------------------------------------------------------------------

               allocate( fmval(nx,ny,nz,nc) )

            do jx = 1, nx

            do ic = 1, nc

               inum = inum + 1

*-----------------------------------------------------------------------

               if( nn .gt. 0 ) then

                  iz = nt(ic) / 1000

                  ia = nt(ic) - iz * 1000

                  if( ia .gt. 0 ) then

                     in = ia - iz

                     call chname(idum,ia,iz,chau)

                  end if

               end if

*-----------------------------------------------------------------------

                  write(iot,'(/"#",78("-"))')

               if( inum .eq. 1 ) then
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

               if( rtfac(m) .lt. 0.d0 )
     &           write(iot,
     &           '( "#restart maximum for normalization:",1p1e14.7)')
     &           facmax(m)

            if( nn .eq. 0 ) then

               write(iot,'("#   no. =",i3,3x,
     &         "ix  =",i3,3x,"Z = all")')
     &                     inum, jx

            else if( ia .eq. 0 ) then

               write(iot,'("#   no. =",i3,3x,
     &         "ix  =",i3,3x,"Z = ",i3," : ",a3)')
     &                     inum, jx, iz, elmnt(iz)

            else

               write(iot,'("#   no. =",i3,3x,
     &         "ix  =",i3,3x,a8)')
     &                     inum, jx, chau

            end if

*-----------------------------------------------------------------------

            if( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) then

               if( nn .eq. 0 ) then

                  write(iot,'(/a1,"no. =",i3,
     &                        ",    ix =",i3,
     &                        ",    Z  =  all", a1)')
     &                        cha, inum, jx, cha

               else if( ia .eq. 0 ) then

                  write(iot,'(/a1,"no. =",i3,
     &                        ",    ix =",i3,
     &                        ",    Z  =  ",i3, a1)')
     &                        cha, inum, jx, iz, cha

               else

                  write(iot,'(/a1,"no. =",i3,
     &                        ",    ix =",i3,
     &                        ",    ",a8, a1)')
     &                        cha, inum, jx, chau, cha

               end if

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: z [cm]")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: y [cm]")')

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

*-----------------------------------------------------------------------

               if( inum .eq. 1 ) then

                  form  = ( ym(ny+1) - ym(1) )
     &                  / ( zm(nz+1) - zm(1) )
                  xfac  = 0.9
                  afac  = 0.8
                  izlog = 1
                  inocm = 1
                  inolg = 1

                  if( itanl(m) .gt. 0 )
     &            call anset(itang(m),itanl(m),
     &                       form,xfac,afac,fmin,fmax,izlog,inocm,inolg)

                  if( form .le. 1.0 ) then

                     scal = form**0.35
                     xfac = xfac / form**0.5
                     xorg = 0.0
                     yorg = min( 20.d0, ( 1.0 / form - 1.0 ) / 2.5 )

                  else

                     scal = 1.0 / form**0.41
                     xfac = xfac / form**0.5
                     xorg = min( 20.d0, ( form - 1.0 ) / 1.5)
                     yorg = 0.0

                  end if

               end if

                  write(iot,'( "set: c1[",f6.3,
     &                            "] c2[",f6.3,
     &                            "] c5[",f6.3,
     &                            "] c6[",f6.3,
     &                            "] c7[",f6.3,
     &                            "] c8[",f6.3,"]")')
     &                    form, xfac, afac, scal, xorg, yorg
                  write(iot,'( "p: h2fs form[c1] xfac[c2]",
     &                 " afac[c5] scal[c6] xorg[c7] yorg[c8] nosp")')

               if( ( ( ittwo(m) .ge. 2 .and. ittwo(m) .le. 3 ) .or.
     &               ( ittwo(m) .ge. 6 .and. ittwo(m) .le. 7 ) ) .and.
     &             fmin .gt. 0.0 .and. fmax .gt. fmin ) then

                if (ioe .eq. 1 ) then
                  write(iot,'( "set: c3[",1p1e15.6,"] c4[",
     &                  1p1e15.6,"]")') fmin, fmax
                else
                  write(iot,'( "set: c3[1.0e-4] c4[1.0]")')
                end if
                  write(iot,'( "p: cmin[c3] cmax[c4]")')
                  write(iot,'( "p: dmin(1e-31)")')

                  if( ittwo(m) .eq. 3 .or. ittwo(m) .eq. 7 )
     &                 write(iot,'( "p: zlog")')

               end if

*-----------------------------------------------------------------------

               zmin = zm(1)
               zmax = zm(nz+1)
               ymin = ym(1)
               ymax = ym(ny+1)

               write(iot,'( "p: xmin(",1p1g14.7,") xmax(",1p1g14.7,
     &                       ")")') zmin, zmax

               write(iot,'( "p: ymin(",1p1g14.7,") ymax(",1p1g14.7,
     &                       ")")') ymin, ymax

               if( itanl(m) .gt. 0 .and. ioe .eq. 1 ) then

                  write(iot,'( "p: ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               else if( itanl(m) .gt. 0 .and. ioe .ne. 1 ) then

                   call terrang(iot,m)

               end if

               if( itsans(m) .gt. 0 ) then

                  call write_sangel(iot,m,0)

               end if

            end if

*-----------------------------------------------------------------------

               do jz = 1, nz
               do jy = 1, ny

                     sek = 0.0

                     if( itunt(m) .eq. 1 ) then

                        vm = 1.0

                     else if( itunt(m) .eq. 2 ) then

                        vm = vl(jx,jy,jz)

                     end if

*-----------------------------------------------------------------------

                  if( nn .eq. 0 ) then

                     do kn = 1, mn
                     do kz = 1, mz
                     do il = 0, mm ! frtati 2022/02/18

                        lz  = ikzz( kz, kn )
                        ln  = iknn( kz, kn )

                        vn  = vm * tr(jx,jy,jz,lz,ln,il,ioe)
                        sek = sek + vn

                     end do
                     end do
                     end do

                  else if( ia .eq. 0 ) then

                     do kn = 1, mn
                     do il = 0, mm ! frtati 2022/02/18

                        lz  = ikzz( iz, kn )
                        ln  = iknn( iz, kn )

                        vn  = vm * tr(jx,jy,jz,lz,ln,il,ioe)
                        sek = sek + vn

                     end do
                     end do

                  else

                        lz  = ikzz( iz, in )
                        ln  = iknn( iz, in )

                     do il = 0, mm ! frtati 2022/02/18
                        vn  = vm * tr(jx,jy,jz,lz,ln,il,ioe)
                        sek = sek + vn
                     end do

                  end if

                        fm(jz,jy) = sek / vm

               end do
               end do

*-----------------------------------------------------------------------

               write(iot,'("#  ny = ",i3,"   nz = ",i3)')
     &                      ny, nz

            if( ittwo(m) .ne. 4 ) then

               write(iot,'( "# ( ( data(z,y), z = 1, nz ),",
     &                      " y = ny, 1, -1 )")')

            end if

*-----------------------------------------------------------------------

            if( ( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) .and.
     &            igsh .eq. 0 ) then

               write(iot,'(/a4," y = ",1p1g14.7," to ",
     &         1p1g14.7," by ",1p1g14.7," ; x = ",
     &         1p1g14.7," to ",1p1g14.7," by ",
     &         1p1g14.7," ;")') dc2,
     &         ym(ny) + rtydl(m)/2.0, ym(1) + rtydl(m)/2.0, rtydl(m),
     &         zm(1) + rtzdl(m)/2.0, zm(nz) + rtzdl(m)/2.0, rtzdl(m)

               write(iot,'(1p10e11.3)')
     &         ( ( fm(jz,jy), jz = 1, nz ), jy = ny, 1, -1 )

*-----------------------------------------------------------------------

            else if( ittwo(m) .eq. 4 ) then

               write(iot,'(/"# y          z        ",
     &                      "  number")')

               do jz = 1, nz
               do jy = 1, ny

                  write(iot,'(1p10e11.3)')
     &               ym(jy)  + rtydl(m)/2.0,
     &               zm(jz)  + rtzdl(m)/2.0,
     &               fm(jz,jy)

               end do
               end do

*-----------------------------------------------------------------------

            else if( ittwo(m) .eq. 5 ) then

               write(iot,'("#   y = ",1p1g14.7," to ",
     &         1p1g14.7," by ",1p1g14.7/
     &                     "#   z = ",1p1g14.7," to ",
     &         1p1g14.7," by ",1p1g14.7)')
     &         ym(1) + rtydl(m)/2.0, ym(ny) + rtydl(m)/2.0, rtydl(m),
     &         zm(1) + rtzdl(m)/2.0, zm(nz) + rtzdl(m)/2.0, rtzdl(m)

                  write(iot,'(/4x,a3,4x,1p1000e11.3)')
     &            'y/z',( zm(jz) + rtzdl(m)/2.0, jz = 1, nz )

               do jy = ny, 1, -1

                  write(iot,'(1p1000e11.3)')
     &            ym(jy) + rtydl(m)/2.0,
     &            ( fm(jz,jy), jz = 1, nz )

               end do

            end if

*-----------------------------------------------------------------------
*        gshow
*-----------------------------------------------------------------------

         if( itgsh(m) .ne. 0 .and.
     &     ( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) ) then

               write(iot,'(/"#",78("-"))')
               write(iot,'( "# gshow")')
               write(iot,'( "#",78("-"))')
               write(iot,'(/"p: legs[c5*0.875]")')

               xval = ( xm(jx) + xm(jx+1) ) / 2.0d0
               none = 1
               iaxs = 2
               iuni = itgsh(m) * 2 - 1 + igsh
               ires = itres(m)
               widt = rtwid(m)

               call gshow(0,iot,iaxs,iuni,widt,ires,
     &                    inum,igsh,igser,gfnam,igfmn,
     &                    nz+1,ny+1,none,zm,ym,xval,ixyz(1),
     &                    1,1,krr,vll,itmtr(m,4),itglt(m))

         end if

*-----------------------------------------------------------------------

      if( inolg .eq. 1 .and. ( ittwo(m) .eq. 3 .or. ittwo(m) .eq. 7 )
     &    .and. fmin .gt. 0.0 .and. fmax .gt. fmin .and.
     &    igsh .eq. 0 ) then

      write(iot,'(/"#",78("-"))')
      write(iot,'(
     &"z: xorg(1.03)"/
     &"p: xfac[c2*0.05] form[c1/0.05] notc noxt noyt dmin(1)"/
     &"p: ymin(0) ymax(1) xmin(0) xmax(1) cmin(1) cmax(100) nosp"/
     &"hc: y= 0.005 to 0.995 by 0.01 ; x= 0.5 to 0.5 by 1 ;"/
     &" 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20"/
     &"21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40"/
     &"41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60"/
     &"61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80"/
     &"81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100"/
     &"z: xorg(0.0)"/
     &"p: xfac[c2*0.05] form[c1/0.05] ",
     &"nosp afac[c5*0.625] ymin[c3] ymax[c4]"/
     &"p: noxt noxn ytxt(-1) ynum(-1) xltd(-1) xstd(-1) itic(1) ylog"
     &)')

      if ( ioe .eq. 2 ) then

          write(iot,'("y: Relative Error")')

      else if( itazl(m) .eq. 0 ) then

         write(iot,'("y: Number ",a15)') hsunit(itunt(m))

      else

         write(iot,'("y: ",200a1)') (itazt(m)(i:i),i=1,itazl(m))

      end if

      if( izlog .eq. 0 ) write(iot,'( "p: ylin")')

      write(iot,'( "#",78("-"))')

      end if

*-----------------------------------------------------------------------

            if( inocm .eq. 1 ) then

               if( nn .eq. 0 ) then

                  write(iot,'(/"wt: s[c5]",/
     &                        a1,"vspace{-3}"/
     &                        "  xmin  &=&",1pe13.4," [cm]"/
     &                        "  xmax  &=&",1pe13.4," [cm]"/
     &                        "      Z &=&   all"/
     &                        "e:")')
     &                        yen,
     &                        xm(jx), xm(jx+1)

               else if( ia .eq. 0 ) then

                  write(iot,'(/"wt: s[c5]",/
     &                        a1,"vspace{-3}"/
     &                        "  xmin  &=&",1pe13.4," [cm]"/
     &                        "  xmax  &=&",1pe13.4," [cm]"/
     &                        "      Z &=& ",i3,"   :  ",a3/
     &                        "e:")')
     &                        yen,
     &                        xm(jx), xm(jx+1), iz, elmnt(iz)

               else

                  write(iot,'(/"wt: s[c5]",/
     &                        a1,"vspace{-3}"/
     &                        "  xmin  &=&",1pe13.4," [cm]"/
     &                        "  xmax  &=&",1pe13.4," [cm]"/
     &                        "  nucl. &=&   ",a8/
     &                        "e:")')
     &                        yen,
     &                        xm(jx), xm(jx+1), chau

               end if

            end if

*-----------------------------------------------------------------------

      if( inolg .eq. 1 .and. ( ittwo(m) .eq. 3 .or. ittwo(m) .eq. 7 )
     &    .and. cmin .gt. 0.0 .and. cmax .gt. cmin .and.
     &    igsh .eq. 0 ) then

         write(iot,'(/"#",78("-"))')
         write(iot,'("z: xorg[-1.03/0.05]"/
     &               "p: form[c1] xfac[c2] afac[c5] nosp notf")')

               write(iot,'( "p: xmin(",1p1g14.7,") xmax(",1p1g14.7,
     &                       ")")') zmin, zmax

               write(iot,'( "p: ymin(",1p1g14.7,") ymax(",1p1g14.7,
     &                       ")")') ymin, ymax

      end if

*-----------------------------------------------------------------------

               if ( itvtk(m) .ne. 0 ) then

                  do jy = 1, ny
                     do jz = 1, nz
                        fmval(jx,jy,jz,ic) = fm(jz,jy)
                     end do
                  end do

               end if

            end do
            end do

*-----------------------------------------------------------------------

            if ( itbmp(m) .ne. 0 ) then

               numIndex = 2
               allocate( bmpfIType(1:numIndex) )
               allocate( bmpfIndex(1:numIndex) )
               bmpfIType = (/ 'n', 'x' /)
               bmpfIndex = (/ nc, nx /)
               bmpWidth  = nz
               bmpHeight = ny

            end if

            if ( itvtk(m) .ne. 0 ) then

               call open_file(isunit_vtk_meta_default,
     &                 "", isunit_vtk_meta, ios, .true.)
               call open_file(isunit_vtk_default,
     &                 "", isunit_vtk, ios, .true.)

               iaxs = 2
               nparam = 1
               ntg = 1

               write(isunit_vtk_meta) ntg

               write(isunit_vtk_meta) iaxs
               write(isunit_vtk_meta) nparam
               write(isunit_vtk_meta) 'n'
               write(isunit_vtk_meta) nc

               write(isunit_vtk_meta) nx,ny,nz
               write(isunit_vtk_meta)
     &                 ( xm(ix), ix=1,nx+1 )
               write(isunit_vtk_meta)
     &                 ( ym(iy), iy=1,ny+1 )
               write(isunit_vtk_meta)
     &                 ( zm(iz), iz=1,nz+1 )

               do ic = 1, nc

                  write(isunit_vtk_meta) ic

                  write(isunit_vtk)
     &                    ( ( ( fmval(ix,iy,iz,ic),
     &                          ix=1,nx ),
     &                          iy=1,ny ),
     &                          iz=1,nz )
               end do

               iuni = itgsh(m) * 2 - 1 + igsh
               ires = itres(m)

               call open_file(isunit_vtk_rm_default,
     &                 "", isunit_vtk_rm, ios, .true.)

               if ( itgsh(m) .ne. 0 .and. ioe .eq. 1 ) then

                  call open_file(isunit_vtk_geom_default,
     &                    "", isunit_vtk_geom, ios, .true.)
                  call open_file(isunit_vtk_geom_meta_default,
     &                    "", isunit_vtk_geom_meta, ios, .true.)

               else
                  isunit_vtk_geom = 0
                  isunit_vtk_geom_meta = 0
               end if

               call vtk_set_cell_regmat(
     &                 isunit_vtk_rm,
     &                 isunit_vtk_geom, isunit_vtk_geom_meta,
     &                 nx, ny, nz, xm, ym, zm,
     &                 iaxs,iuni,ires, igser,
     &                 1,1,krr,vll,itmtr(m,4))

            end if

*-----------------------------------------------------------------------
*        xz axis (matrix)
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 11 ) then

               if( ittwo(m) .eq. 1 ) then

                  dc2 = 'h2: '

               else if( ittwo(m) .eq. 2 ) then

                  dc2 = 'hd: '

               else if( ittwo(m) .eq. 3 ) then

                  dc2 = 'hc: '

               else if( ittwo(m) .eq. 6 ) then

                  dc2 = 'hd2:'

               else if( ittwo(m) .eq. 7 ) then

                  dc2 = 'hc2:'

               end if

*-----------------------------------------------------------------------

               inum = 0

                  if( nn .eq. 0 .or. igsh .ne. 0  ) then

                     nc = 1

                  else

                     nc = nn

                  end if

*-----------------------------------------------------------------------
*              max and min
*-----------------------------------------------------------------------

                        fmax = 0.0
                        fmin = 1.e+33

                  do jy = 1, ny
                  do ic = 1, nc

                     if( nn .gt. 0 ) then

                        iz = nt(ic) / 1000

                        ia = nt(ic) - iz * 1000

                        if( ia .gt. 0 ) in = ia - iz

                     end if

                     do jz = 1, nz
                     do jx = 1, nx

                           sek = 0.0

                           if( itunt(m) .eq. 1 ) then

                              vm = 1.0

                           else if( itunt(m) .eq. 2 ) then

                              vm = vl(jx,jy,jz)

                           end if

                        if( nn .eq. 0 ) then

                           do kn = 1, mn
                           do kz = 1, mz
                           do il = 0, mm ! frtati 2022/02/18

                              lz  = ikzz( kz, kn )
                              ln  = iknn( kz, kn )

                              vn  = vm * tr(jx,jy,jz,lz,ln,il,1)
                              sek = sek + vn

                           end do
                           end do
                           end do

                        else if( ia .eq. 0 ) then

                           do kn = 1, mn
                           do il = 0, mm ! frtati 2022/02/18

                              lz  = ikzz( iz, kn )
                              ln  = iknn( iz, kn )

                              vn  = vm * tr(jx,jy,jz,lz,ln,il,1)
                              sek = sek + vn

                           end do
                           end do

                        else

                              lz  = ikzz( iz, in )
                              ln  = iknn( iz, in )

                           do il = 0, mm ! frtati 2022/02/18
                              vn  = vm * tr(jx,jy,jz,lz,ln,il,1)
                              sek = sek + vn
                           end do

                        end if

                              fmsv = sek / vm

                              if( fmsv .gt. fmax ) fmax = fmsv
                              if( fmsv .gt. 0.0 .and.
     &                            fmsv .lt. fmin ) fmin = fmsv

                     end do
                     end do

                  end do
                  end do

*-----------------------------------------------------------------------

               allocate( fmval(nx,ny,nz,nc) )

            do jy = 1, ny

            do ic = 1, nc

               inum = inum + 1

*-----------------------------------------------------------------------

               if( nn .gt. 0 ) then

                  iz = nt(ic) / 1000

                  ia = nt(ic) - iz * 1000

                  if( ia .gt. 0 ) then

                     in = ia - iz

                     call chname(idum,ia,iz,chau)

                  end if

               end if

*-----------------------------------------------------------------------

                  write(iot,'(/"#",78("-"))')

               if( inum .eq. 1 ) then
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

               if( rtfac(m) .lt. 0.d0 )
     &           write(iot,
     &           '( "#restart maximum for normalization:",1p1e14.7)')
     &           facmax(m)

            if( nn .eq. 0 ) then

               write(iot,'("#   no. =",i3,3x,
     &         "iy  =",i3,3x,"Z = all")')
     &                     inum, jy

            else if( ia .eq. 0 ) then

               write(iot,'("#   no. =",i3,3x,
     &         "iy  =",i3,3x,"Z = ",i3," : ",a3)')
     &                     inum, jy, iz, elmnt(iz)

            else

               write(iot,'("#   no. =",i3,3x,
     &         "iy  =",i3,3x,a8)')
     &                     inum, jy, chau

            end if

*-----------------------------------------------------------------------

            if( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) then

               if( nn .eq. 0 ) then

                  write(iot,'(/a1,"no. =",i3,
     &                        ",    iy =",i3,
     &                        ",    Z  =  all", a1)')
     &                        cha, inum, jy, cha

               else if( ia .eq. 0 ) then

                  write(iot,'(/a1,"no. =",i3,
     &                        ",    iy =",i3,
     &                        ",    Z  =  ",i3, a1)')
     &                        cha, inum, jy, iz, cha

               else

                  write(iot,'(/a1,"no. =",i3,
     &                        ",    iy =",i3,
     &                        ",    ",a8, a1)')
     &                        cha, inum, jy, chau, cha

               end if

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: z [cm]")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: x [cm]")')

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

*-----------------------------------------------------------------------

               if( inum .eq. 1 ) then

                  form  = ( xm(nx+1) - xm(1) )
     &                  / ( zm(nz+1) - zm(1) )
                  xfac  = 0.9
                  afac  = 0.8
                  izlog = 1
                  inocm = 1
                  inolg = 1

                  if( itanl(m) .gt. 0 )
     &            call anset(itang(m),itanl(m),
     &                       form,xfac,afac,fmin,fmax,izlog,inocm,inolg)

                  if( form .le. 1.0 ) then

                     scal = form**0.35
                     xfac = xfac / form**0.5
                     xorg = 0.0
                     yorg = min( 20.d0, ( 1.0 / form - 1.0 ) / 2.5 )

                  else

                     scal = 1.0 / form**0.41
                     xfac = xfac / form**0.5
                     xorg = min( 20.d0, ( form - 1.0 ) / 1.5)
                     yorg = 0.0

                  end if

               end if

                  write(iot,'( "set: c1[",f6.3,
     &                            "] c2[",f6.3,
     &                            "] c5[",f6.3,
     &                            "] c6[",f6.3,
     &                            "] c7[",f6.3,
     &                            "] c8[",f6.3,"]")')
     &                    form, xfac, afac, scal, xorg, yorg
                  write(iot,'( "p: h2fs form[c1] xfac[c2]",
     &                 " afac[c5] scal[c6] xorg[c7] yorg[c8] nosp")')

               if( ( ( ittwo(m) .ge. 2 .and. ittwo(m) .le. 3 ) .or.
     &               ( ittwo(m) .ge. 6 .and. ittwo(m) .le. 7 ) ) .and.
     &             fmin .gt. 0.0 .and. fmax .gt. fmin ) then

                if (ioe .eq. 1 ) then
                  write(iot,'( "set: c3[",1p1e15.6,"] c4[",
     &                  1p1e15.6,"]")') fmin, fmax
                else
                  write(iot,'( "set: c3[1.0e-4] c4[1.0]")')
                end if
                  write(iot,'( "p: cmin[c3] cmax[c4]")')
                  write(iot,'( "p: dmin(1e-31)")')

                  if( ittwo(m) .eq. 3 .or. ittwo(m) .eq. 7 )
     &                 write(iot,'( "p: zlog")')

               end if

*-----------------------------------------------------------------------

               zmin = zm(1)
               zmax = zm(nz+1)
               xmin = xm(1)
               xmax = xm(nx+1)

               write(iot,'( "p: xmin(",1p1g14.7,") xmax(",1p1g14.7,
     &                       ")")') zmin, zmax

               write(iot,'( "p: ymin(",1p1g14.7,") ymax(",1p1g14.7,
     &                       ")")') xmin, xmax

               if( itanl(m) .gt. 0 .and. ioe .eq. 1 ) then

                  write(iot,'( "p: ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               else if( itanl(m) .gt. 0 .and. ioe .ne. 1 ) then

                   call terrang(iot,m)

               end if

               if( itsans(m) .gt. 0 ) then

                  call write_sangel(iot,m,0)

               end if

            end if

*-----------------------------------------------------------------------

               do jz = 1, nz
               do jx = 1, nx

                     sek = 0.0

                     if( itunt(m) .eq. 1 ) then

                        vm = 1.0

                     else if( itunt(m) .eq. 2 ) then

                        vm = vl(jx,jy,jz)

                     end if

*-----------------------------------------------------------------------

                  if( nn .eq. 0 ) then

                     do kn = 1, mn
                     do kz = 1, mz
                     do il = 0, mm ! frtati 2022/02/18

                        lz  = ikzz( kz, kn )
                        ln  = iknn( kz, kn )

                        vn  = vm * tr(jx,jy,jz,lz,ln,il,ioe)
                        sek = sek + vn

                     end do
                     end do
                     end do

                  else if( ia .eq. 0 ) then

                     do kn = 1, mn
                     do il = 0, mm ! frtati 2022/02/18

                        lz  = ikzz( iz, kn )
                        ln  = iknn( iz, kn )

                        vn  = vm * tr(jx,jy,jz,lz,ln,il,ioe)
                        sek = sek + vn

                     end do
                     end do

                  else

                        lz  = ikzz( iz, in )
                        ln  = iknn( iz, in )

                     do il = 0, mm ! frtati 2022/02/18
                        vn  = vm * tr(jx,jy,jz,lz,ln,il,ioe)
                        sek = sek + vn
                     end do

                  end if


                        fm(jz,jx) = sek / vm

               end do
               end do

*-----------------------------------------------------------------------

               write(iot,'("#  nx = ",i3,"   nz = ",i3)')
     &                      nx, nz

            if( ittwo(m) .ne. 4 ) then

               write(iot,'( "# ( ( data(z,x), z = 1, nz ),",
     &                      " x = nx, 1, -1 )")')

            end if

*-----------------------------------------------------------------------

            if( ( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) .and.
     &            igsh .eq. 0 ) then

               write(iot,'(/a4," y = ",1p1g14.7," to ",
     &         1p1g14.7," by ",1p1g14.7," ; x = ",
     &         1p1g14.7," to ",1p1g14.7," by ",
     &         1p1g14.7," ;")') dc2,
     &         xm(nx) + rtxdl(m)/2.0, xm(1) + rtxdl(m)/2.0, rtxdl(m),
     &         zm(1) + rtzdl(m)/2.0, zm(nz) + rtzdl(m)/2.0, rtzdl(m)

               write(iot,'(1p10e11.3)')
     &         ( ( fm(jz,jx), jz = 1, nz ), jx = nx, 1, -1 )

*-----------------------------------------------------------------------

            else if( ittwo(m) .eq. 4 ) then

               write(iot,'(/"# x          z        ",
     &                      "  number")')

               do jz = 1, nz
               do jx = 1, nx

                  write(iot,'(1p10e11.3)')
     &               xm(jx)  + rtxdl(m)/2.0,
     &               zm(jz)  + rtzdl(m)/2.0,
     &               fm(jz,jx)

               end do
               end do

*-----------------------------------------------------------------------

            else if( ittwo(m) .eq. 5 ) then

               write(iot,'("#   x = ",1p1g14.7," to ",
     &         1p1g14.7," by ",1p1g14.7/
     &                     "#   z = ",1p1g14.7," to ",
     &         1p1g14.7," by ",1p1g14.7)')
     &         xm(1) + rtxdl(m)/2.0, xm(nx) + rtxdl(m)/2.0, rtxdl(m),
     &         zm(1) + rtzdl(m)/2.0, zm(nz) + rtzdl(m)/2.0, rtzdl(m)

                  write(iot,'(/4x,a3,4x,1p1000e11.3)')
     &            'x/z',( zm(jz) + rtzdl(m)/2.0, jz = 1, nz )

               do jx = nx, 1, -1

                  write(iot,'(1p1000e11.3)')
     &            xm(jx) + rtxdl(m)/2.0,
     &            ( fm(jz,jx), jz = 1, nz )

               end do

            end if

*-----------------------------------------------------------------------
*        gshow
*-----------------------------------------------------------------------

         if( itgsh(m) .ne. 0 .and.
     &     ( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) ) then

               write(iot,'(/"#",78("-"))')
               write(iot,'( "# gshow")')
               write(iot,'( "#",78("-"))')
               write(iot,'(/"p: legs[c5*0.875]")')

               yval = ( ym(jy) + ym(jy+1) ) / 2.0d0
               none = 1
               iaxs = 3
               iuni = itgsh(m) * 2 - 1 + igsh
               ires = itres(m)
               widt = rtwid(m)

               call gshow(0,iot,iaxs,iuni,widt,ires,
     &                    inum,igsh,igser,gfnam,igfmn,
     &                    nz+1,nx+1,none,zm,xm,yval,ixyz(1),
     &                    1,1,krr,vll,itmtr(m,4),itglt(m))

         end if

*-----------------------------------------------------------------------

      if( inolg .eq. 1 .and. ( ittwo(m) .eq. 3 .or. ittwo(m) .eq. 7 )
     &    .and. fmin .gt. 0.0 .and. fmax .gt. fmin .and.
     &    igsh .eq. 0 ) then

      write(iot,'(/"#",78("-"))')
      write(iot,'(
     &"z: xorg(1.03)"/
     &"p: xfac[c2*0.05] form[c1/0.05] notc noxt noyt dmin(1)"/
     &"p: ymin(0) ymax(1) xmin(0) xmax(1) cmin(1) cmax(100) nosp"/
     &"hc: y= 0.005 to 0.995 by 0.01 ; x= 0.5 to 0.5 by 1 ;"/
     &" 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20"/
     &"21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40"/
     &"41 42 43 44 45 46 47 48 49 50 51 52 53 54 55 56 57 58 59 60"/
     &"61 62 63 64 65 66 67 68 69 70 71 72 73 74 75 76 77 78 79 80"/
     &"81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96 97 98 99 100"/
     &"z: xorg(0.0)"/
     &"p: xfac[c2*0.05] form[c1/0.05] ",
     &"nosp afac[c5*0.625] ymin[c3] ymax[c4]"/
     &"p: noxt noxn ytxt(-1) ynum(-1) xltd(-1) xstd(-1) itic(1) ylog"
     &)')

      if ( ioe .eq. 2 ) then

          write(iot,'("y: Relative Error")')

      else if( itazl(m) .eq. 0 ) then

         write(iot,'("y: Number ",a15)') hsunit(itunt(m))

      else

         write(iot,'("y: ",200a1)') (itazt(m)(i:i),i=1,itazl(m))

      end if

      if( izlog .eq. 0 ) write(iot,'( "p: ylin")')

      write(iot,'( "#",78("-"))')

      end if

*-----------------------------------------------------------------------

            if( inocm .eq. 1 ) then

               if( nn .eq. 0 ) then

                  write(iot,'(/"wt: s[c5]",/
     &                        a1,"vspace{-3}"/
     &                        "  ymin  &=&",1pe13.4," [cm]"/
     &                        "  ymax  &=&",1pe13.4," [cm]"/
     &                        "      Z &=&   all"/
     &                        "e:")')
     &                        yen,
     &                        ym(jy), ym(jy+1)

               else if( ia .eq. 0 ) then

                  write(iot,'(/"wt: s[c5]",/
     &                        a1,"vspace{-3}"/
     &                        "  ymin  &=&",1pe13.4," [cm]"/
     &                        "  ymax  &=&",1pe13.4," [cm]"/
     &                        "      Z &=& ",i3,"   :  ",a3/
     &                        "e:")')
     &                        yen,
     &                        ym(jy), ym(jy+1), iz, elmnt(iz)

               else

                  write(iot,'(/"wt: s[c5]",/
     &                        a1,"vspace{-3}"/
     &                        "  ymin  &=&",1pe13.4," [cm]"/
     &                        "  ymax  &=&",1pe13.4," [cm]"/
     &                        "  nucl. &=&   ",a8/
     &                        "e:")')
     &                        yen,
     &                        ym(jy), ym(jy+1), chau

               end if

            end if

*-----------------------------------------------------------------------

      if( inolg .eq. 1 .and. ( ittwo(m) .eq. 3 .or. ittwo(m) .eq. 7 )
     &    .and. cmin .gt. 0.0 .and. cmax .gt. cmin .and.
     &    igsh .eq. 0 ) then

         write(iot,'(/"#",78("-"))')
         write(iot,'("z: xorg[-1.03/0.05]"/
     &               "p: form[c1] xfac[c2] afac[c5] nosp notf")')

               write(iot,'( "p: xmin(",1p1g14.7,") xmax(",1p1g14.7,
     &                       ")")') zmin, zmax

               write(iot,'( "p: ymin(",1p1g14.7,") ymax(",1p1g14.7,
     &                       ")")') xmin, xmax

      end if

*-----------------------------------------------------------------------

               if ( itvtk(m) .ne. 0 ) then

                  do jx = 1, nx
                     do jz = 1, nz
                        fmval(jx,jy,jz,ic) = fm(jz,jx)
                     end do
                  end do

               end if

            end do
            end do

*-----------------------------------------------------------------------

            if ( itbmp(m) .ne. 0 ) then

               numIndex = 2
               allocate( bmpfIType(1:numIndex) )
               allocate( bmpfIndex(1:numIndex) )
               bmpfIType = (/ 'n', 'y' /)
               bmpfIndex = (/ nc, ny /)
               bmpWidth  = nz
               bmpHeight = nx

            end if

            if ( itvtk(m) .ne. 0 ) then

               call open_file(isunit_vtk_meta_default,
     &                 "", isunit_vtk_meta, ios, .true.)
               call open_file(isunit_vtk_default,
     &                 "", isunit_vtk, ios, .true.)

               iaxs = 3
               nparam = 1
               ntg = 1

               write(isunit_vtk_meta) ntg

               write(isunit_vtk_meta) iaxs
               write(isunit_vtk_meta) nparam
               write(isunit_vtk_meta) 'n'
               write(isunit_vtk_meta) nc

               write(isunit_vtk_meta) nx,ny,nz
               write(isunit_vtk_meta)
     &                 ( xm(ix), ix=1,nx+1 )
               write(isunit_vtk_meta)
     &                 ( ym(iy), iy=1,ny+1 )
               write(isunit_vtk_meta)
     &                 ( zm(iz), iz=1,nz+1 )

               do ic = 1, nc

                  write(isunit_vtk_meta) ic

                  write(isunit_vtk)
     &                    ( ( ( fmval(ix,iy,iz,ic),
     &                          ix=1,nx ),
     &                          iy=1,ny ),
     &                          iz=1,nz )
               end do

               iuni = itgsh(m) * 2 - 1 + igsh
               ires = itres(m)

               call open_file(isunit_vtk_rm_default,
     &                 "", isunit_vtk_rm, ios, .true.)

               if ( itgsh(m) .ne. 0 .and. ioe .eq. 1 ) then

                  call open_file(isunit_vtk_geom_default,
     &                    "", isunit_vtk_geom, ios, .true.)
                  call open_file(isunit_vtk_geom_meta_default,
     &                    "", isunit_vtk_geom_meta, ios, .true.)

               else
                  isunit_vtk_geom = 0
                  isunit_vtk_geom_meta = 0
               end if

               call vtk_set_cell_regmat(
     &                 isunit_vtk_rm,
     &                 isunit_vtk_geom, isunit_vtk_geom_meta,
     &                 nx, ny, nz, xm, ym, zm,
     &                 iaxs,iuni,ires, igser,
     &                 1,1,krr,vll,itmtr(m,4))

            end if

*-----------------------------------------------------------------------

         end if

            call prestart(m,iot) !OBINATA(2012.8.20)

            close(iot)

*-----------------------------------------------------------------------

         if( iteps(m) .ne. 0 ) then

            open(iot, file = fname, status = 'unknown' )
            call a_angel(idasa,fname)

         end if

         if ( itbmp(m) .ne. 0 ) then

            open(iot, file = fname, status = 'unknown' )
            call bitmap_a_angel(
     &              idasa, fname,
     &              bmpWidth, bmpHeight,
     &              bmpfIType, bmpfIndex, numIndex)
            close(iot)

         end if

         if ( itvtk(m) .ne. 0 ) then

            rewind(isunit_vtk_meta)
            rewind(isunit_vtk)

            call vtk_write_tally(
     &              itvtkfmt(m),
     &              isunit_vtk_meta, isunit_vtk,
     &              isunit_vtk_rm,
     &              fname)

            call close_file(isunit_vtk_meta)
            call close_file(isunit_vtk)

            if ( itgsh(m) .ne. 0 .and. ioe .eq. 1 .and.
     &           isunit_vtk_geom.gt.0 .and.
     &           isunit_vtk_geom_meta.gt.0 ) then

               rewind(isunit_vtk_geom)
               rewind(isunit_vtk_geom_meta)

               isText = (itvtkfmt(m).eq.0)

               call vtk_create_filename(
     &                 fname, outFilename, .true., 0, 0)
               call open_file(iunit_vtk_g_default,
     &                 outFilename, iunit_vtk_g, ios, isText)
               call vtk_write_polydata(
     &                 itvtkfmt(m), iunit_vtk_g,
     &                 isunit_vtk_geom, isunit_vtk_geom_meta)
               call close_file(iunit_vtk_g)

            end if

         end if

         if ( allocated(fmval)     ) deallocate( fmval )
         if ( allocated(bmpfIType) ) deallocate( bmpfIType )
         if ( allocated(bmpfIndex) ) deallocate( bmpfIndex )

         call close_file(isunit_vtk_rm)
         call close_file(isunit_vtk_geom)
         call close_file(isunit_vtk_geom_meta)

*-----------------------------------------------------------------------

  900 continue

*-----------------------------------------------------------------------
      deallocate( ixyz )
      return
      end


************************************************************************
*                                                                      *
      subroutine prodxs2(kf,iz,im,ein,istat,mx,mz,ma,xx,sigpt,
     &                   il,mk,icl,jcl,mtallyID)
*                                                                      *
*     last modified by K.Niita on 2018/02/22                           *
*                                                                      *
*     This subroutine calculate nuclide production cross section       *
*     with an evaluated cross section data set.                        *
*     Data for neutron and proton incidence on LLFP                    *
*     are currently available.                                         *
*                                                                      *
*     input:  kf    - the kf code of incident particle                 *
*             iz    - atomic number of target nucleus                  *
*             im    - mass number of target nucleus                    *
*             ein   - incident energy in MeV                           *
*             mk    - material number                                  *
*             icl   - cell number                                      *
*             jcl   - jcoll ID                                         *
*             mtallID - tally ID                                       *
*                                                                      *
*     output  istat - status                                           *
*                   0  normal end                                      *
*                   1  the incident particle is not proton nor neutron,*
*                      deuteron, alpha, or gamma(photon)               *
*                   2  no cross section data are given for the target  *
*                      nucleus                                         *
*                   3  the incident energy is out of range(axis=dchain)*
*                   4  cross section given by model or nuclear data    *
*                      in getflt is not found (when sigpt=0)           *
*                                                                      *
*     output                                                           *
*          mx      - number of produced nuclides                       *
*          mxprodxs- maximum number of the following arrays            *
*          mz(mxprodxs) - atomic number of the m-th nuclide            *
*          ma(mxprodxs) - mass number of the m-th nuclide              *
*          xx(mxprodxs) - multiplicity for the m-th nuclide            *
*          il(mxprodxs) - 0,1,2 : excitation level                     *
*          sigpt   - reaction cross section (b)                        *
*                                                                      *
************************************************************************

      use NDATA2MOD
      use MEMBANKMOD
      use GGMARRAYMOD, only: fme, jmd
      use moddas_material
      implicit real*8(a-h,o-z)

      include 'param.inc'
      include 'ggmparam.inc'
      include 'err.inc'

      common /paraj/ mstz(300), parz(300)

      common /tall10/ itaxn(itlmax), itaxs(itlmax,6) ! T.Sato 2024/09/03

      common /pnint/ ipnint
      common /kmat1k/ kmatd(kvlmax), kmate(kvlmax)
      common /kmat1ka/ kmathd(kvlmax), kmathe(kvlmax)

      common /emode/  emodem, ge1, ge2, iemode

*-----------------------------------------------------------------------

      dimension mz(mxprodxs),ma(mxprodxs),xx(mxprodxs),il(mxprodxs)

! T.Sato 2024/09/03, ignore neutron below 20 MeV data for axis = dchain
      if(kf.eq.2112.and.ein.lt.20.0001.and.itaxs(mtallyID,1).eq.13) then
       istat = 3
       return
      endif

*-----------------------------------------------------------------------
*     check data
*-----------------------------------------------------------------------

      if( itallo .eq. 0 ) then

         istat = 2
         return

      end if

*-----------------------------------------------------------------------
*     calculate a cross section value
*     check incident particle type (p, n, d, a, g)
*-----------------------------------------------------------------------

      izam = iz * 1000 + im
      mx = 0

      if ( kf .eq. 2212 ) then  ! proton
         ityp = 1

         sigpt = 0d0
         if( jcl .eq. 4 ) then

          dmaxn = das_kmatd(kmatd(mk)+1)
          dminn = das_kmatd(kmatd(mk)+2)

          if ( ein .gt. dmaxn ) then

            call sigrc(ityp,ein,im,iz,sigt,signe,sigel)
            sigpt = signe

          else

            call sig_tot(ityp,sigt,sigaa,ein,mk)
            mii = 0
            do ii = jmd(1+mk), jmd(1+mk+1) - 1
               mii = mii + 1
               if( izam .eq. isigza(mii) ) then
                  sigpt =  siggcn(mii) /fme(ii)
                  exit
               end if
            end do
            dsmax = das_kmate(kmate(mk)+(mii-1)*5+11)

            if ( ein .gt. dsmax ) then ! high-energy emode
               call sigrc(ityp,ein,im,iz,sigt,signe,sigel)
               sigpt = signe
            end if

          end if

         else if( jcl .eq. 9 ) then
            call sig_tot(ityp,sigt,sigaa,ein,mk)
            mii = 0
            do ii = jmd(1+mk), jmd(1+mk+1) - 1
               mii = mii + 1
               if( izam .eq. isigza(mii) ) then
                  sigpt = siggc(mii) /fme(ii)
                  exit
               end if
            end do

         else
            write(ErrCha,*) '*** Error in NDATA=2,3, for proton'
            ErrID = 'L:18086/R:prodxs2/F:talls03.f' !E03_108_001
            call ErrWrite(ErrID,ErrCha)
            write(*,*)
     &           'Warning : Unexpected jcoll is found for proton'//
     &           ' in [t-yield] with ndata=2,3'

         end if
         if( sigpt .le. 0d0 ) then
            istat = 4
            return
         end if

         do i = 1, ndatp
          if( indatp(i) .eq. izam ) then
           do j = 1, jdatp(i)

            if( ein .gt. cxproton(i,j,1,1) .and.
     &             ein .le. cxproton(i,j,1,kedatp(i,j)) ) then

             do k = 2, kedatp(i,j)
                if( ein .le. cxproton(i,j,1,k) ) goto 100
             end do
             goto 110

 100         continue
             mx = mx + 1
             mz(mx) = kndatp(i,j) / 1000
             ma(mx) = kndatp(i,j)
     &            - kndatp(i,j) / 1000 * 1000
             intdat = kidatp(i,j)

             if( intdat .eq. 1 ) then

                xx(mx) = cxproton(i,j,2,k)

             else if( intdat .eq. 2 ) then

                xx(mx) = cxproton(i,j,2,k-1)
     &               + ( cxproton(i,j,2,k)-cxproton(i,j,2,k-1) )
     &               * ( ein - cxproton(i,j,1,k-1) )
     &               / ( cxproton(i,j,1,k)-cxproton(i,j,1,k-1) )

             else if( intdat .eq. 3 ) then

                xx(mx) = cxproton(i,j,2,k-1)
     &               + ( cxproton(i,j,2,k)-cxproton(i,j,2,k-1) )
     &               * log( ein
     &               / cxproton(i,j,1,k-1) )
     &               / log( cxproton(i,j,1,k)
     &               / cxproton(i,j,1,k-1) )

             else if( intdat .eq. 4 ) then

                xx(mx) = cxproton(i,j,2,k-1)
     &               * exp( log( cxproton(i,j,2,k)
     &               / cxproton(i,j,2,k-1) )
     &               * ( ein - cxproton(i,j,1,k-1) )
     &               / ( cxproton(i,j,1,k)-cxproton(i,j,1,k-1)))

             else if( intdat .eq. 5 ) then

                xx(mx) = cxproton(i,j,2,k-1)
     &               * exp( log( cxproton(i,j,2,k)
     &               / cxproton(i,j,2,k-1) )
     &               * log( ein
     &               / cxproton(i,j,1,k-1) )
     &               / log( cxproton(i,j,1,k)
     &               / cxproton(i,j,1,k-1) ) )

             end if

             xx(mx) = xx(mx) / sigpt
             il(mx) = kldatp(i,j)

            end if
 110        continue

           end do
          end if
         end do

*-----------------------------------------------------------------------

      else if( kf .eq. 2112 ) then ! neutron

         sigpt = 0d0
         if( jcl .eq. 4 ) then

          dmaxn = das_kmatd(kmatd(mk)+3)
          dminn = das_kmatd(kmatd(mk)+4)

          if ( ein .gt. dmaxn ) then

            ityp = 2
            call sigrc(ityp,ein,im,iz,sigt,signe,sigel)
            sigpt = signe

          else

            tme = 0.0
            call xstneu(0,sigt,sigaa,icl,ein,tme,mk)
            mii = 0
            do ii = jmd(1+mk), jmd(1+mk+1) - 1
               mii = mii + 1
               if( izam .eq. isigza(mii) ) then
                  sigpt = siggcn(mii) /fme(ii)
                  exit
               end if
            end do
            dsmax = das_kmate(kmate(mk)+(mii-1)*5+12)

            if ( ein .gt. dsmax ) then ! high-energy emode
               ityp = 2
               call sigrc(ityp,ein,im,iz,sigt,signe,sigel)
               sigpt = signe
            end if

          end if

         else if( jcl .eq. 6 .or. jcl .eq. 10 ) then
            tme = 0.0
            call xstneu(0,sigt,sigaa,icl,ein,tme,mk)
            mii = 0
            do ii = jmd(1+mk), jmd(1+mk+1) - 1
               mii = mii + 1
               if( izam .eq. isigza(mii) ) then
                  sigpt = siggc(mii) /fme(ii)
                  exit
               end if
            end do

         else
            write(ErrCha,*) '*** Error in NDATA=2,3, for neutron'
            ErrID = 'L:18219/R:prodxs2/F:talls03.f' !E03_108_002
            call ErrWrite(ErrID,ErrCha)
            write(*,*)
     &           'Warning : Unexpected jcoll is found for neutron'//
     &           ' in [t-yield] with ndata=2,3'

         end if
         if( sigpt .le. 0d0 ) then
            istat = 4
            return
         end if

         do i = 1, ndatn
          if( indatn(i) .eq. izam ) then
           do j = 1, jdatn(i)

            if( ein .gt. cxneutron(i,j,1,1) .and.
     &             ein .le. cxneutron(i,j,1,kedatn(i,j)) ) then

             do k = 2, kedatn(i,j)
                if( ein .le. cxneutron(i,j,1,k) ) goto 200
             end do
             goto 210

 200         continue
             mx = mx + 1
             mz(mx) = kndatn(i,j) / 1000
             ma(mx) = kndatn(i,j)
     &            - kndatn(i,j) / 1000 * 1000
             intdat = kidatn(i,j)

             if( intdat .eq. 1 ) then

                xx(mx) = cxneutron(i,j,2,k)

             else if( intdat .eq. 2 ) then

                xx(mx) = cxneutron(i,j,2,k-1)
     &               + ( cxneutron(i,j,2,k)-cxneutron(i,j,2,k-1) )
     &               * ( ein - cxneutron(i,j,1,k-1) )
     &               / ( cxneutron(i,j,1,k)-cxneutron(i,j,1,k-1) )

             else if( intdat .eq. 3 ) then

                xx(mx) = cxneutron(i,j,2,k-1)
     &               + ( cxneutron(i,j,2,k)-cxneutron(i,j,2,k-1) )
     &               * log( ein
     &               / cxneutron(i,j,1,k-1) )
     &               / log( cxneutron(i,j,1,k)
     &               / cxneutron(i,j,1,k-1) )

             else if( intdat .eq. 4 ) then

                xx(mx) = cxneutron(i,j,2,k-1)
     &               * exp( log( cxneutron(i,j,2,k)
     &               / cxneutron(i,j,2,k-1) )
     &               * ( ein - cxneutron(i,j,1,k-1) )
     &               / ( cxneutron(i,j,1,k)-cxneutron(i,j,1,k-1)))

             else if( intdat .eq. 5 ) then

                xx(mx) = cxneutron(i,j,2,k-1)
     &               * exp( log( cxneutron(i,j,2,k)
     &               / cxneutron(i,j,2,k-1) )
     &               * log( ein
     &               / cxneutron(i,j,1,k-1) )
     &               / log( cxneutron(i,j,1,k)
     &               / cxneutron(i,j,1,k-1) ) )

             end if

             xx(mx) = xx(mx) / sigpt
             il(mx) = kldatn(i,j)

            end if
 210        continue

           end do
          end if
         end do

*-----------------------------------------------------------------------

      else if( kf .eq. 1000002 ) then ! deuteron
         ityp = 15

         sigpt = 0d0
         if( jcl .eq. 5 ) then

          einmevpern = ein/2.0
          dmaxn = das_kmathd(kmathd(mk)+3)
          dminn = das_kmathd(kmathd(mk)+4)

          if ( einmevpern .gt. dmaxn ) then

            ap = 2.0
            zp = 1.0
            at = im
            zt = iz
            call sighi(ap,zp,ein,at,zt,signe,sigel,bmax)
            sigpt = signe

          else

            call sig_tot(ityp,sigt,sigaa,ein,mk)
            mii = 0
            do ii = jmd(1+mk), jmd(1+mk+1) - 1
               mii = mii + 1
               if( izam .eq. isigza(mii) ) then
                  sigpt = siggc(mii) /fme(ii)
                  exit
               end if
            end do
            dsmax = das_kmate(kmate(mk)+(mii-1)*5+14)

            if ( einmevpern .gt. dsmax ) then ! high-energy emode
               ap = 2.0
               zp = 1.0
               at = im
               zt = iz
               call sighi(ap,zp,ein,at,zt,signe,sigel,bmax)
               sigpt = signe
            end if

          end if

         else if( jcl .eq. 9 ) then
            call sig_tot(ityp,sigt,sigaa,ein,mk)
            mii = 0
            do ii = jmd(1+mk), jmd(1+mk+1) - 1
               mii = mii + 1
               if( izam .eq. isigza(mii) ) then
                  sigpt = siggc(mii) /fme(ii)
                  exit
               end if
            end do

         else
            write(ErrCha,*) '*** Error in NDATA=2,3, for deuteron'
            ErrID = 'L:18358/R:prodxs2/F:talls03.f' !E03_108_001
            call ErrWrite(ErrID,ErrCha)
            write(*,*)
     &           'Warning : Unexpected jcoll is found for deuteron'//
     &           ' in [t-yield] with ndata=2,3'

         end if
         if( sigpt .le. 0d0 ) then
            istat = 4
            return
         end if

         do i = 1, ndatd
          if( indatd(i) .eq. izam ) then
           do j = 1, jdatd(i)

            if( ein .gt. cxdeuteron(i,j,1,1) .and.
     &             ein .le. cxdeuteron(i,j,1,kedatd(i,j)) ) then

             do k = 2, kedatd(i,j)
                if( ein .le. cxdeuteron(i,j,1,k) ) goto 300
             end do
             goto 310

 300         continue
             mx = mx + 1
             mz(mx) = kndatd(i,j) / 1000
             ma(mx) = kndatd(i,j)
     &            - kndatd(i,j) / 1000 * 1000
             intdat = kidatd(i,j)

             if( intdat .eq. 1 ) then

                xx(mx) = cxdeuteron(i,j,2,k)

             else if( intdat .eq. 2 ) then

                xx(mx) = cxdeuteron(i,j,2,k-1)
     &               + ( cxdeuteron(i,j,2,k)-cxdeuteron(i,j,2,k-1) )
     &               * ( ein - cxdeuteron(i,j,1,k-1) )
     &               / ( cxdeuteron(i,j,1,k)-cxdeuteron(i,j,1,k-1) )

             else if( intdat .eq. 3 ) then

                xx(mx) = cxdeuteron(i,j,2,k-1)
     &               + ( cxdeuteron(i,j,2,k)-cxdeuteron(i,j,2,k-1) )
     &               * log( ein
     &               / cxdeuteron(i,j,1,k-1) )
     &               / log( cxdeuteron(i,j,1,k)
     &               / cxdeuteron(i,j,1,k-1) )

             else if( intdat .eq. 4 ) then

                xx(mx) = cxdeuteron(i,j,2,k-1)
     &               * exp( log( cxdeuteron(i,j,2,k)
     &               / cxdeuteron(i,j,2,k-1) )
     &               * ( ein - cxdeuteron(i,j,1,k-1) )
     &               / ( cxdeuteron(i,j,1,k)-cxdeuteron(i,j,1,k-1)))

             else if( intdat .eq. 5 ) then

                xx(mx) = cxdeuteron(i,j,2,k-1)
     &               * exp( log( cxdeuteron(i,j,2,k)
     &               / cxdeuteron(i,j,2,k-1) )
     &               * log( ein
     &               / cxdeuteron(i,j,1,k-1) )
     &               / log( cxdeuteron(i,j,1,k)
     &               / cxdeuteron(i,j,1,k-1) ) )

             end if

             xx(mx) = xx(mx) / sigpt
             il(mx) = kldatd(i,j)

            end if
 310        continue

           end do
          end if
         end do

*-----------------------------------------------------------------------

      else if( kf .eq. 2000004 ) then ! alpha
         ityp = 18

         sigpt = 0d0
         if( jcl .eq. 5 ) then

          einmevpern = ein/4.0
          dmaxn = das_kmathd(kmathd(mk)+5)
          dminn = das_kmathd(kmathd(mk)+6)

          if ( einmevpern .gt. dmaxn ) then

            ap = 4.0
            zp = 2.0
            at = im
            zt = iz
            call sighi(ap,zp,ein,at,zt,signe,sigel,bmax)
            sigpt = signe

          else

            call sig_tot(ityp,sigt,sigaa,ein,mk)
            mii = 0
            do ii = jmd(1+mk), jmd(1+mk+1) - 1
               mii = mii + 1
               if( izam .eq. isigza(mii) ) then
                  sigpt = siggc(mii) /fme(ii)
                  exit
               end if
            end do
            dsmax = das_kmate(kmate(mk)+(mii-1)*5+15)

            if ( einmevpern .gt. dsmax ) then ! high-energy emode
               ap = 4.0
               zp = 2.0
               at = im
               zt = iz
               call sighi(ap,zp,ein,at,zt,signe,sigel,bmax)
               sigpt = signe
            end if

          end if

         else if( jcl .eq. 9 ) then
            call sig_tot(ityp,sigt,sigaa,ein,mk)
            mii = 0
            do ii = jmd(1+mk), jmd(1+mk+1) - 1
               mii = mii + 1
               if( izam .eq. isigza(mii) ) then
                  sigpt = siggc(mii) /fme(ii)
                  exit
               end if
            end do

         else
            write(ErrCha,*) '*** Error in NDATA=2,3, for alpha'
            ErrID = 'L:18497/R:prodxs2/F:talls03.f' !E03_108_001
            call ErrWrite(ErrID,ErrCha)
            write(*,*)
     &           'Warning : Unexpected jcoll is found for alpha'//
     &           ' in [t-yield] with ndata=2,3'

         end if
         if( sigpt .le. 0d0 ) then
            istat = 4
            return
         end if

         do i = 1, ndata
          if( indata(i) .eq. izam ) then
           do j = 1, jdata(i)

            if( ein .gt. cxalpha(i,j,1,1) .and.
     &             ein .le. cxalpha(i,j,1,kedata(i,j)) ) then

             do k = 2, kedata(i,j)
                if( ein .le. cxalpha(i,j,1,k) ) goto 400
             end do
             goto 410

 400         continue
             mx = mx + 1
             mz(mx) = kndata(i,j) / 1000
             ma(mx) = kndata(i,j)
     &            - kndata(i,j) / 1000 * 1000
             intdat = kidata(i,j)

             if( intdat .eq. 1 ) then

                xx(mx) = cxalpha(i,j,2,k)

             else if( intdat .eq. 2 ) then

                xx(mx) = cxalpha(i,j,2,k-1)
     &               + ( cxalpha(i,j,2,k)-cxalpha(i,j,2,k-1) )
     &               * ( ein - cxalpha(i,j,1,k-1) )
     &               / ( cxalpha(i,j,1,k)-cxalpha(i,j,1,k-1) )

             else if( intdat .eq. 3 ) then

                xx(mx) = cxalpha(i,j,2,k-1)
     &               + ( cxalpha(i,j,2,k)-cxalpha(i,j,2,k-1) )
     &               * log( ein
     &               / cxalpha(i,j,1,k-1) )
     &               / log( cxalpha(i,j,1,k)
     &               / cxalpha(i,j,1,k-1) )

             else if( intdat .eq. 4 ) then

                xx(mx) = cxalpha(i,j,2,k-1)
     &               * exp( log( cxalpha(i,j,2,k)
     &               / cxalpha(i,j,2,k-1) )
     &               * ( ein - cxalpha(i,j,1,k-1) )
     &               / ( cxalpha(i,j,1,k)-cxalpha(i,j,1,k-1)))

             else if( intdat .eq. 5 ) then

                xx(mx) = cxalpha(i,j,2,k-1)
     &               * exp( log( cxalpha(i,j,2,k)
     &               / cxalpha(i,j,2,k-1) )
     &               * log( ein
     &               / cxalpha(i,j,1,k-1) )
     &               / log( cxalpha(i,j,1,k)
     &               / cxalpha(i,j,1,k-1) ) )

             end if

             xx(mx) = xx(mx) / sigpt
             il(mx) = kldata(i,j)

            end if
 410        continue

           end do
          end if
         end do

*-----------------------------------------------------------------------

      else if( kf .eq. 22 ) then ! gamma (photon)

         sigpt = 0d0
         if( ipnint .ge. 1 ) then

          if( jcl .eq. 9 ) then
           call pnctot(mk,totmpni,ein)
           mii = 0
           do ii = jmd(1+mk), jmd(1+mk+1) - 1
              mii = mii + 1
              if( izam .eq. isigza(mii) ) then
                 sigpt = siggc(mii) /fme(ii)
                 exit
              end if
           end do

          else if( jcl .eq. 15 ) then

           dmaxn = das_kmathd(kmathd(mk)+1)
           dminn = das_kmathd(kmathd(mk)+2)

           if ( ein .gt. dmaxn ) then

            call xstpni(0,sigt,ein,mk)
            do ii = jmd(1+mk), jmd(1+mk+1) - 1
               icntm = ii - jmd(1+mk) + 1
               if( izam .eq. isigza(icntm) ) then
                  sigpt = siggpn(icntm) /fme(ii)
                  exit
               end if
            end do

           else

            call pnctot(mk,totmpni,ein)
            mii = 0
            do ii = jmd(1+mk), jmd(1+mk+1) - 1
               mii = mii + 1
               if( izam .eq. isigza(mii) ) then
                  sigpt = siggc(mii) /fme(ii)
                  exit
               end if
            end do
            dsmax = das_kmate(kmate(mk)+(mii-1)*5+13)

            if ( ein .gt. dsmax ) then ! high-energy emode
             call xstpni(0,sigt,ein,mk)
             do ii = jmd(1+mk), jmd(1+mk+1) - 1
                icntm = ii - jmd(1+mk) + 1
                if( izam .eq. isigza(icntm) ) then
                   sigpt = siggpn(icntm) /fme(ii)
                   exit
                end if
             end do
            end if

           end if

          else
             write(ErrCha,*) '*** Error in NDATA=2,3, for gamma'
             ErrID = 'L:18640/R:prodxs2/F:talls03.f' !E03_108_003
             call ErrWrite(ErrID,ErrCha)
             write(*,*)
     &            'Warning : Unexpected jcoll is found for gamma'//
     &            ' in [t-yield] with ndata=2,3'

          end if

         endif


         if ( sigpt .gt. 0d0 ) then
          do i = 1, ndatg
           if( indatg(i) .eq. izam ) then
            do j = 1, jdatg(i)

             if( ein .gt. cxgamma(i,j,1,1) .and.
     &              ein .le. cxgamma(i,j,1,kedatg(i,j)) ) then

              do k = 2, kedatg(i,j)
                 if( ein .le. cxgamma(i,j,1,k) ) goto 500
              end do
              goto 510

 500          continue
              mx = mx + 1
              mz(mx) = kndatg(i,j) / 1000
              ma(mx) = kndatg(i,j)
     &             - kndatg(i,j) / 1000 * 1000
              intdat = kidatg(i,j)

              if( intdat .eq. 1 ) then

                 xx(mx) = cxgamma(i,j,2,k)

              else if( intdat .eq. 2 ) then

                 xx(mx) = cxgamma(i,j,2,k-1)
     &                + ( cxgamma(i,j,2,k)-cxgamma(i,j,2,k-1) )
     &                * ( ein - cxgamma(i,j,1,k-1) )
     &                / ( cxgamma(i,j,1,k)-cxgamma(i,j,1,k-1) )

              else if( intdat .eq. 3 ) then

                 xx(mx) = cxgamma(i,j,2,k-1)
     &                + ( cxgamma(i,j,2,k)-cxgamma(i,j,2,k-1) )
     &                * log( ein
     &                / cxgamma(i,j,1,k-1) )
     &                / log( cxgamma(i,j,1,k)
     &                / cxgamma(i,j,1,k-1) )

              else if( intdat .eq. 4 ) then

                 xx(mx) = cxgamma(i,j,2,k-1)
     &                * exp( log( cxgamma(i,j,2,k)
     &                / cxgamma(i,j,2,k-1) )
     &                * ( ein - cxgamma(i,j,1,k-1) )
     &                / ( cxgamma(i,j,1,k)-cxgamma(i,j,1,k-1)))

              else if( intdat .eq. 5 ) then

                 xx(mx) = cxgamma(i,j,2,k-1)
     &                * exp( log( cxgamma(i,j,2,k)
     &                / cxgamma(i,j,2,k-1) )
     &                * log( ein
     &                / cxgamma(i,j,1,k-1) )
     &                / log( cxgamma(i,j,1,k)
     &                / cxgamma(i,j,1,k-1) ) )

              end if

              xx(mx) = xx(mx) / sigpt
              il(mx) = kldatg(i,j)

             end if
 510         continue

            end do
           end if
          end do
         end if

*-----------------------------------------------------------------------

      else

         istat = 1
         return

      end if

*-----------------------------------------------------------------------

      if( mx .eq. 0 ) then
         istat = 2
      else
         istat = 0
      end if

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine prodxs(kf,iz,im,ein,istat,mx,mz,ma,xx,signe,il)
*                                                                      *
*     coded by F.Maekawa on 2001/12/12                                 *
*     last modified by K.Niita on 2005/08/10                           *
*                                                                      *
*     This subroutine calculate nuclide production cross section       *
*     with an evaluated cross section data set.                        *
*     Data for neutron and proton incidence on He-4, N-14 and O-16     *
*     are currently available.                                         *
*                                                                      *
*     input:  kf    - the kf code of incident particle                 *
*             iz    - atomic number of target nucleus                  *
*             im    - mass number of target nucleus                    *
*             ein   - incident energy in MeV                           *
*                                                                      *
*     output  istat - status                                           *
*                   0  normal end                                      *
*                   1  the incident particle is not ptoron nor neutron *
*                   2  no cross section data are given for the target  *
*                      nucleus                                         *
*                   3  the incident energy is out of range             *
*                                                                      *
*     output                                                           *
*          mx      - number of produced nuclides                       *
*          mxprodxs- maximum number of the following arrays            *
*          mz(mxprodxs) - atomic number of the m-th nuclide            *
*          ma(mxprodxs) - mass number of the m-th nuclide              *
*          xx(mxprodxs) - multiplicity for the m-th nuclide            *
*          il(mxprodxs) - 0,1,2 : excitation level                     *
*          signe   - reaction cross section (mb)                       *
*                                                                      *
*     ifc: file assignment number                                      *
*     ix:  maximum number of energy points                             *
*     jx:  maximum number of particle types                            *
*     kx:  maximum number of reactions                                 *
*     lx:  maximum number of interpolation intervals                   *
*                                                                      *
************************************************************************

      implicit real*8(a-h,o-z)

      include 'param.inc'

*-----------------------------------------------------------------------

      parameter (ix=150,jx=2,kx=20,lx=1)

*-----------------------------------------------------------------------

      common /xs_com/ nmax(jx),kmax(jx),en(jx,ix),id(jx,7,kx),
     &                xn(jx,kx,ix),nbt(lx),inx(lx)
!$OMP THREADPRIVATE(/xs_com/)

      dimension mz(mxprodxs),ma(mxprodxs),xx(mxprodxs),il(mxprodxs)

      data ifst /0/
      save ifst !FURUTA
!$OMP THREADPRIVATE(ifst)
*-----------------------------------------------------------------------
*     initial zero set
*-----------------------------------------------------------------------

      if( ifst .eq. 0 ) then

         ifst = ifst + 1

         do i = 1, ix
         do j = 1, kx
         do k = 1, jx

            if( xn(k,j,i) .lt. 1.d-3 ) xn(k,j,i) = 1.d-30

         end do
         end do
         end do

      end if

*-----------------------------------------------------------------------
*     check incident particle type (neutron or proton)
*-----------------------------------------------------------------------

      if (kf.eq.2112) then

         j=1
         ityp = 2

      else if (kf.eq.2212) then

         j=2
         ityp = 1

      else

         istat=1
         return

      endif

*-----------------------------------------------------------------------
*     check target nucleus (He-4, N-14, O-16)
*-----------------------------------------------------------------------

      if (iz.eq.2 .and. im.eq. 4) goto 21
      if (iz.eq.7 .and. im.eq.14) goto 21
      if (iz.eq.8 .and. im.eq.16) goto 21

         istat=2
         return

   21 continue

*-----------------------------------------------------------------------
*     check incident energy
*-----------------------------------------------------------------------

      if (ein.lt.en(j,1) .or. ein.gt.en(j,nmax(j))) then

        istat=3
        return

      endif

*-----------------------------------------------------------------------
*     calculate a cross section value
*-----------------------------------------------------------------------

      np=nmax(j)
      nr=1
      nbt(1)=nmax(j)
      inx(1)=5
      mx=0

          call sigrc(ityp,ein,im,iz,sigt,signe,sigel)
            signe = signe * 1000.0

      do 31 k=1,kmax(j)

        if (iz.eq.id(j,1,k) .and. im.eq.id(j,2,k)) then

          mx=mx+1
          mz(mx)=id(j,5,k)
          ma(mx)=id(j,6,k)
          call interpol(j,k,np,nr,ein,xsec)
          xx(mx)=xsec-1.0d-30
          if( xx(mx) .lt. 1.0d-30 ) xx(mx) = 0.0d0

          xx(mx) = xx(mx) / signe

          il(mx) = 0

        endif

   31 continue

         istat = 0

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine interpol(jjj,kkk,np,nr,ein,xsec)
*                                                                      *
************************************************************************

      implicit real*8(a-h,o-z)

      parameter (ix=150,jx=2,kx=20,lx=1)

      common /xs_com/ nmax(jx),kmax(jx),en(jx,ix),id(jx,7,kx),
     &                xn(jx,kx,ix),nbt(lx),inx(lx)
!$OMP THREADPRIVATE(/xs_com/)

        do 2 j=np-1,1,-1
          if (ein.gt.en(jjj,j)) then
            d=en(jjj,j+1)-en(jjj,j)
            if (d.le.0.0) then
              xsec=0.0
              return
            endif
            inte=0
            do 8 n=1,nr
              if (j+1.le.nbt(n)) then
                inte=inx(n)
                goto 9
              endif
    8       continue
    9       continue
            if (inte.eq.1) then
              xsec=xn(jjj,kkk,j+1)
              return
            endif
            if (inte.lt.2 .or. inte.gt.5) stop ' inte error'
            if (inte.eq.2 .or. inte.eq.4) then
              e1=en(jjj,j)
              e2=en(jjj,j+1)
              ei=ein
            else
              e1=log(en(jjj,j))
              e2=log(en(jjj,j+1))
              ei=log(ein)
            endif
            if (inte.eq.2 .or. inte.eq.3) then
              x1=xn(jjj,kkk,j)
              x2=xn(jjj,kkk,j+1)
            else
              x1=log(xn(jjj,kkk,j))
              x2=log(xn(jjj,kkk,j+1))
            endif
            d=e2-e1
            xsec=((ei-e1)*x2+(e2-ei)*x1)/d
            if (inte.eq.2 .or. inte.eq.3) then
              xsec=xsec
            else
              xsec=exp(xsec)
            endif
            goto 3
          endif
    2   continue
    3   continue

      return
      end

************************************************************************
*                                                                      *
      block data heon00
*                                                                      *
*      Data for neutron and proton incidence on He-4, N-14 and O-16    *
*      are currently available.                                        *
*                                                                      *
*      ix:  maximum number of energy points                            *
*      jx:  maximum number of particle types                           *
*      kx:  maximum number of reactions                                *
*      lx:  maximum number of interpolation intervals                  *
*                                                                      *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      parameter (ix=150,jx=2,kx=20,lx=1)
      parameter (iall = ix*jx*kx )

      common /xs_com/ nmax(jx),kmax(jx),en(jx,ix),id(jx,7,kx),
     &                xn(jx,kx,ix),nbt(lx),inx(lx)
!$OMP THREADPRIVATE(/xs_com/)

*-----------------------------------------------------------------------


*-----------------------------------------------------------------------
*    neutron
*-----------------------------------------------------------------------

       data kmax(1) / 19 /

       data nmax(1) / 142 /
       data (en(1,n),n=1,142) /
     & 1.000E+00, 1.110E+00, 1.111E+00, 1.170E+00, 1.171E+00, 1.230E+00,
     & 1.231E+00, 1.290E+00, 1.291E+00, 1.350E+00, 1.351E+00, 1.420E+00,
     & 1.421E+00, 1.500E+00, 1.501E+00, 1.570E+00, 1.571E+00, 1.650E+00,
     & 1.651E+00, 1.740E+00, 1.741E+00, 1.830E+00, 1.831E+00, 1.920E+00,
     & 1.921E+00, 2.020E+00, 2.021E+00, 2.120E+00, 2.121E+00, 2.230E+00,
     & 2.231E+00, 2.310E+00, 2.311E+00, 2.350E+00, 2.351E+00, 2.370E+00,
     & 2.371E+00, 2.390E+00, 2.391E+00, 2.470E+00, 2.471E+00, 2.590E+00,
     & 2.591E+00, 2.730E+00, 2.731E+00, 2.870E+00, 2.871E+00, 3.010E+00,
     & 3.011E+00, 3.170E+00, 3.171E+00, 3.330E+00, 3.331E+00, 3.680E+00,
     & 3.681E+00, 4.070E+00, 4.071E+00, 4.490E+00, 4.491E+00, 4.720E+00,
     & 4.721E+00, 4.970E+00, 4.971E+00, 5.220E+00, 5.221E+00, 5.490E+00,
     & 5.491E+00, 5.770E+00, 5.771E+00, 6.070E+00, 6.071E+00, 6.380E+00,
     & 6.381E+00, 6.590E+00, 6.591E+00, 6.700E+00, 6.701E+00, 7.050E+00,
     & 7.051E+00, 7.410E+00, 7.411E+00, 7.790E+00, 7.791E+00, 8.190E+00,
     & 8.191E+00, 8.610E+00, 8.611E+00, 9.050E+00, 9.051E+00, 9.510E+00,
     & 9.511E+00, 1.000E+01, 1.001E+01, 1.050E+01, 1.051E+01, 1.110E+01,
     & 1.111E+01, 1.160E+01, 1.161E+01, 1.220E+01, 1.221E+01, 1.250E+01,
     & 1.251E+01, 1.280E+01, 1.281E+01, 1.350E+01, 1.351E+01, 1.380E+01,
     & 1.381E+01, 1.420E+01, 1.421E+01, 1.460E+01, 1.461E+01, 1.490E+01,
     & 1.491E+01, 1.570E+01, 1.571E+01, 1.650E+01, 1.651E+01, 1.690E+01,
     & 1.691E+01, 1.730E+01, 1.731E+01, 1.960E+01, 2.000E+01, 2.500E+01,
     & 3.000E+01, 3.500E+01, 4.000E+01, 5.000E+01, 6.000E+01, 7.000E+01,
     & 8.000E+01, 9.000E+01, 1.000E+02, 1.500E+02, 3.000E+02, 6.000E+02,
     & 1.000E+03, 3.000E+03, 1.000E+04, 1.000E+10/

*-----------------------------------------------------------------------
*    He-4(n,x)H-3
*-----------------------------------------------------------------------
       data (id(1,n,1),n=1,7) /
     &     2,    4,    0,    1,    1,    3,   19/
       data (xn(1,1,i),i=124,142) /
     & 1.000E-03, 5.000E+01, 8.000E+01, 6.700E+01, 5.700E+01, 5.000E+01,
     & 4.500E+01, 4.000E+01, 3.700E+01, 3.500E+01, 3.200E+01, 3.000E+01,
     & 2.700E+01, 2.500E+01, 2.000E+01, 1.600E+01, 1.600E+01, 1.600E+01,
     & 1.600E+01/

*-----------------------------------------------------------------------
*    N-14(n,x)H-3
*-----------------------------------------------------------------------
       data (id(1,n,2),n=1,7) /
     &     7,   14,    0,    1,    1,    3,   87/
       data (xn(1,2,i),i=56,142) /
     & 1.000E-03, 2.141E-02, 2.141E-02, 2.319E-01, 2.319E-01, 8.758E-01,
     & 8.758E-01, 2.342E+00, 2.342E+00, 6.048E+00, 6.048E+00, 1.859E+01,
     & 1.859E+01, 1.940E+01, 1.940E+01, 2.005E+01, 2.005E+01, 2.383E+01,
     & 2.383E+01, 2.755E+01, 2.755E+01, 2.688E+01, 2.688E+01, 2.196E+01,
     & 2.196E+01, 1.781E+01, 1.781E+01, 1.607E+01, 1.607E+01, 1.720E+01,
     & 1.720E+01, 1.807E+01, 1.807E+01, 1.703E+01, 1.703E+01, 1.922E+01,
     & 1.922E+01, 2.529E+01, 2.529E+01, 2.839E+01, 2.839E+01, 3.006E+01,
     & 3.006E+01, 3.103E+01, 3.103E+01, 3.132E+01, 3.132E+01, 3.126E+01,
     & 3.126E+01, 3.047E+01, 3.047E+01, 2.972E+01, 2.972E+01, 2.882E+01,
     & 2.882E+01, 2.774E+01, 2.774E+01, 2.660E+01, 2.660E+01, 2.532E+01,
     & 2.532E+01, 2.006E+01, 2.006E+01, 1.672E+01, 1.672E+01, 1.477E+01,
     & 1.477E+01, 1.024E+01, 1.024E+01, 1.024E+01, 1.024E+01, 1.024E+01,
     & 1.024E+01, 1.040E+01, 1.560E+01, 1.790E+01, 1.910E+01, 1.950E+01,
     & 2.190E+01, 2.240E+01, 2.610E+01, 2.640E+01, 4.820E+01, 6.570E+01,
     & 6.900E+01, 6.900E+01, 6.900E+01/

*-----------------------------------------------------------------------
*    N-14(n,x)Be-7
*-----------------------------------------------------------------------
       data (id(1,n,3),n=1,7) /
     &     7,   14,    0,    1,    4,    7,   16/
       data (xn(1,3,i),i=127,142) /
     & 1.000E-03, 2.770E-01, 8.640E-01, 1.730E+00, 2.850E+00, 3.540E+00,
     & 3.670E+00, 3.620E+00, 3.670E+00, 4.010E+00, 4.690E+00, 7.900E+00,
     & 9.300E+00, 7.140E+00, 7.140E+00, 7.140E+00/

*-----------------------------------------------------------------------
*    N-14(n,x)Be-11
*-----------------------------------------------------------------------
       data (id(1,n,4),n=1,7) /
     &     7,   14,    0,    1,    4,   11,   15/
       data (xn(1,4,i),i=128,142) /
     & 1.000E-03, 6.520E-02, 1.630E-01, 1.630E-01, 1.630E-01, 2.120E-01,
     & 2.930E-01, 1.470E-01, 1.960E-01, 2.440E-01, 1.790E-01, 8.150E-02,
     & 4.890E-02, 4.890E-02, 4.890E-02/

*-----------------------------------------------------------------------
*    N-14(n,x)C-10
*-----------------------------------------------------------------------
       data (id(1,n,5),n=1,7) /
     &     7,   14,    0,    1,    6,   10,   14/
       data (xn(1,5,i),i=129,142) /
     & 1.000E-03, 1.630E-02, 8.150E-02, 1.470E-01, 6.520E-02, 2.770E-01,
     & 3.100E-01, 2.440E-01, 4.400E-01, 7.820E-01, 1.260E+00, 9.450E-01,
     & 9.450E-01, 9.450E-01/

*-----------------------------------------------------------------------
*    N-14(n,x)C-11
*-----------------------------------------------------------------------
       data (id(1,n,6),n=1,7) /
     &     7,   14,    0,    1,    6,   11,   18/
       data (xn(1,6,i),i=125,142) /
     & 1.000E-03, 1.630E-02, 1.430E+00, 1.960E+00, 3.620E+00, 4.400E+00,
     & 4.250E+00, 4.420E+00, 3.980E+00, 4.040E+00, 3.850E+00, 4.040E+00,
     & 4.870E+00, 7.580E+00, 9.040E+00, 6.470E+00, 6.470E+00, 6.470E+00/

*-----------------------------------------------------------------------
*    N-14(n,x)C-14
*-----------------------------------------------------------------------
       data (id(1,n,7),n=1,7) /
     &     7,   14,    0,    1,    6,   14,  142/
       data (xn(1,7,i),i=1,142) /
     & 7.280E+00, 7.280E+00, 6.992E+00, 6.992E+00, 1.738E+01, 1.738E+01,
     & 1.901E+01, 1.901E+01, 7.358E+01, 7.358E+01, 2.340E+02, 2.340E+02,
     & 9.283E+01, 9.283E+01, 2.808E+01, 2.808E+01, 1.990E+01, 1.990E+01,
     & 1.521E+01, 1.521E+01, 1.274E+01, 1.274E+01, 1.231E+01, 1.231E+01,
     & 1.176E+01, 1.176E+01, 1.213E+01, 1.213E+01, 3.584E+01, 3.584E+01,
     & 4.147E+01, 4.147E+01, 1.889E+01, 1.889E+01, 1.557E+01, 1.557E+01,
     & 1.422E+01, 1.422E+01, 1.354E+01, 1.354E+01, 1.827E+01, 1.827E+01,
     & 5.108E+01, 5.108E+01, 7.231E+01, 7.231E+01, 3.155E+01, 3.155E+01,
     & 3.059E+01, 3.059E+01, 3.603E+01, 3.603E+01, 6.565E+01, 6.565E+01,
     & 9.979E+01, 9.979E+01, 9.325E+01, 9.325E+01, 5.090E+01, 5.090E+01,
     & 3.275E+01, 3.275E+01, 2.857E+01, 2.857E+01, 1.980E+01, 1.980E+01,
     & 1.701E+01, 1.701E+01, 1.759E+01, 1.759E+01, 1.776E+01, 1.776E+01,
     & 1.773E+01, 1.773E+01, 1.442E+01, 1.442E+01, 1.620E+01, 1.620E+01,
     & 1.427E+01, 1.427E+01, 1.460E+01, 1.460E+01, 2.264E+01, 2.264E+01,
     & 2.756E+01, 2.756E+01, 3.220E+01, 3.220E+01, 4.516E+01, 4.516E+01,
     & 5.505E+01, 5.505E+01, 5.821E+01, 5.821E+01, 6.022E+01, 6.022E+01,
     & 6.117E+01, 6.117E+01, 5.894E+01, 5.894E+01, 5.690E+01, 5.690E+01,
     & 5.419E+01, 5.419E+01, 4.886E+01, 4.886E+01, 4.640E+01, 4.640E+01,
     & 4.438E+01, 4.438E+01, 4.252E+01, 4.252E+01, 4.064E+01, 4.064E+01,
     & 3.863E+01, 3.863E+01, 3.384E+01, 3.384E+01, 3.140E+01, 3.140E+01,
     & 2.987E+01, 2.987E+01, 2.571E+01, 2.571E+01, 2.570E+01, 3.970E+01,
     & 3.300E+01, 2.820E+01, 2.480E+01, 1.920E+01, 1.560E+01, 1.330E+01,
     & 1.210E+01, 1.080E+01, 9.700E+00, 6.310E+00, 3.370E+00, 1.810E+00,
     & 9.780E-01, 8.470E-01, 8.470E-01, 8.470E-01/

*-----------------------------------------------------------------------
*    N-14(n,x)N-13
*-----------------------------------------------------------------------
       data (id(1,n,8),n=1,7) /
     &     7,   14,    0,    1,    7,   13,   47/
       data (xn(1,8,i),i=96,142) /
     & 1.000E-03, 2.832E-02, 2.832E-02, 6.528E-01, 6.528E-01, 1.583E+00,
     & 1.583E+00, 2.495E+00, 2.495E+00, 4.107E+00, 4.107E+00, 5.087E+00,
     & 5.087E+00, 5.963E+00, 5.963E+00, 6.863E+00, 6.863E+00, 7.666E+00,
     & 7.666E+00, 8.358E+00, 8.358E+00, 9.574E+00, 9.574E+00, 1.011E+01,
     & 1.011E+01, 1.035E+01, 1.035E+01, 1.026E+01, 1.026E+01, 1.200E+01,
     & 1.800E+01, 1.600E+01, 1.470E+01, 1.470E+01, 1.470E+01, 1.460E+01,
     & 1.450E+01, 1.460E+01, 1.430E+01, 1.520E+01, 1.760E+01, 2.110E+01,
     & 2.790E+01, 2.930E+01, 2.310E+01, 2.310E+01, 2.310E+01/

*-----------------------------------------------------------------------
*    O-16(n,x)H-3
*-----------------------------------------------------------------------
       data (id(1,n,9),n=1,7) /
     &     8,   16,    0,    1,    1,    3,   29/
       data (xn(1,9,i),i=114,142) /
     & 1.000E-03, 1.567E-02, 1.567E-02, 2.194E+00, 2.194E+00, 2.551E+00,
     & 2.551E+00, 2.689E+00, 2.689E+00, 2.977E+00, 2.977E+00, 5.000E+00,
     & 1.000E+01, 1.500E+01, 2.000E+01, 2.500E+01, 2.770E+01, 2.770E+01,
     & 2.580E+01, 2.330E+01, 2.320E+01, 2.270E+01, 2.160E+01, 2.340E+01,
     & 4.980E+01, 7.190E+01, 8.430E+01, 8.430E+01, 8.430E+01/

*-----------------------------------------------------------------------
*    O-16(n,x)Be-7
*-----------------------------------------------------------------------
       data (id(1,n,10),n=1,7) /
     &     8,   16,    0,    1,    4,    7,   15/
       data (xn(1,10,i),i=128,142) /
     & 1.000E-03, 1.880E-01, 4.690E+00, 3.770E+00, 3.070E+00, 3.140E+00,
     & 3.140E+00, 3.090E+00, 3.330E+00, 4.040E+00, 6.970E+00, 8.470E+00,
     & 6.940E+00, 6.940E+00, 6.940E+00/

*-----------------------------------------------------------------------
*    O-16(n,x)Be-11
*-----------------------------------------------------------------------
       data (id(1,n,11),n=1,7) /
     &     8,   16,    0,    1,    4,   11,   16/
       data (xn(1,11,i),i=127,142) /
     & 1.000E-03, 6.820E-02, 1.530E-01, 3.070E-01, 4.600E-01, 4.770E-01,
     & 3.240E-01, 3.920E-01, 3.070E-01, 2.900E-01, 1.710E-01, 2.560E-01,
     & 1.530E-01, 1.710E-01, 1.710E-01, 1.710E-01/

*-----------------------------------------------------------------------
*    O-16(n,x)C-10
*-----------------------------------------------------------------------
       data (id(1,n,12),n=1,7) /
     &     8,   16,    0,    1,    6,   10,   14/
       data (xn(1,12,i),i=129,142) /
     & 1.000E-03, 6.820E-02, 6.820E-02, 1.530E-01, 2.900E-01, 3.240E-01,
     & 2.560E-01, 3.240E-01, 5.970E-01, 1.070E+00, 1.160E+00, 1.350E+00,
     & 1.350E+00, 1.350E+00/

*-----------------------------------------------------------------------
*    O-16(n,x)C-11
*-----------------------------------------------------------------------
       data (id(1,n,13),n=1,7) /
     &     8,   16,    0,    1,    6,   11,   17/
       data (xn(1,13,i),i=126,142) /
     & 1.000E-03, 2.220E-01, 3.890E+00, 8.130E+00, 1.360E+01, 1.000E+01,
     & 8.780E+00, 8.830E+00, 7.380E+00, 8.370E+00, 7.350E+00, 6.750E+00,
     & 7.880E+00, 8.870E+00, 5.580E+00, 5.580E+00, 5.580E+00/

*-----------------------------------------------------------------------
*    O-16(n,x)C-14
*-----------------------------------------------------------------------
       data (id(1,n,14),n=1,7) /
     &     8,   16,    0,    1,    6,   14,   29/
       data (xn(1,14,i),i=114,142) /
     & 1.000E-03, 3.381E-05, 3.381E-05, 7.741E-02, 7.741E-02, 1.510E-01,
     & 1.510E-01, 1.900E-01, 1.900E-01, 3.790E-01, 3.790E-01, 7.670E-01,
     & 3.290E+00, 6.460E+00, 7.400E+00, 8.180E+00, 8.050E+00, 7.470E+00,
     & 6.750E+00, 6.220E+00, 6.020E+00, 5.560E+00, 5.390E+00, 4.660E+00,
     & 4.760E+00, 4.010E+00, 3.340E+00, 3.340E+00, 3.340E+00/

*-----------------------------------------------------------------------
*    O-16(n,x)C-15
*-----------------------------------------------------------------------
       data (id(1,n,15),n=1,7) /
     &     8,   16,    0,    1,    6,   15,   18/
       data (xn(1,15,i),i=125,142) /
     & 1.000E-03, 5.290E-01, 6.310E-01, 9.380E-01, 1.160E+00, 9.890E-01,
     & 6.820E-01, 6.650E-01, 6.990E-01, 6.480E-01, 4.260E-01, 3.580E-01,
     & 2.560E-01, 5.120E-02, 4.500E-02, 3.410E-02, 3.410E-02, 3.410E-02/

*-----------------------------------------------------------------------
*    O-16(n,x)N-13
*-----------------------------------------------------------------------
       data (id(1,n,16),n=1,7) /
     &     8,   16,    0,    1,    7,   13,   17/
       data (xn(1,16,i),i=126,142) /
     & 1.000E-03, 2.220E-01, 2.300E+00, 3.000E+00, 4.500E+00, 4.540E+00,
     & 4.330E+00, 3.560E+00, 3.510E+00, 3.510E+00, 3.070E+00, 4.110E+00,
     & 5.370E+00, 5.420E+00, 3.720E+00, 3.720E+00, 3.720E+00/

*-----------------------------------------------------------------------
*    O-16(n,x)N-16
*-----------------------------------------------------------------------
       data (id(1,n,17),n=1,7) /
     &     8,   16,    0,    1,    7,   16,   51/
       data (xn(1,17,i),i=92,142) /
     & 1.000E-03, 1.186E-01, 1.186E-01, 2.910E+00, 2.910E+00, 2.463E+01,
     & 2.463E+01, 5.944E+01, 5.944E+01, 2.410E+01, 2.410E+01, 2.610E+01,
     & 2.610E+01, 4.352E+01, 4.352E+01, 4.781E+01, 4.781E+01, 4.674E+01,
     & 4.674E+01, 4.320E+01, 4.320E+01, 3.675E+01, 3.675E+01, 3.788E+01,
     & 3.788E+01, 3.300E+01, 3.300E+01, 2.918E+01, 2.918E+01, 2.703E+01,
     & 2.703E+01, 2.381E+01, 2.381E+01, 1.854E+01, 1.333E+01, 1.061E+01,
     & 8.000E+00, 6.569E+00, 4.527E+00, 3.201E+00, 2.547E+00, 2.000E+00,
     & 1.600E+00, 1.349E+00, 1.000E+00, 7.438E-01, 5.303E-01, 5.000E-01,
     & 5.000E-01, 5.000E-01, 5.000E-01/

*-----------------------------------------------------------------------
*    O-16(n,x)O-14
*-----------------------------------------------------------------------
       data (id(1,n,18),n=1,7) /
     &     8,   16,    0,    1,    8,   14,   18/
       data (xn(1,18,i),i=125,142) /
     & 1.000E-03, 1.000E-30, 1.710E-02, 5.120E-02, 2.220E-01, 3.750E-01,
     & 3.580E-01, 4.600E-01, 4.600E-01, 4.940E-01, 5.630E-01, 6.140E-01,
     & 1.280E+00, 2.570E+00, 3.000E+00, 2.370E+00, 2.370E+00, 2.370E+00/

*-----------------------------------------------------------------------
*    O-16(n,x)O-15
*-----------------------------------------------------------------------
       data (id(1,n,19),n=1,7) /
     &     8,   16,    0,    1,    8,   15,   25/
       data (xn(1,19,i),i=118,142) /
     & 1.000E-03, 5.677E-02, 5.677E-02, 3.462E-01, 3.462E-01, 1.489E+00,
     & 1.489E+00, 3.070E+00, 1.110E+01, 1.930E+01, 2.050E+01, 1.570E+01,
     & 1.440E+01, 1.500E+01, 1.550E+01, 1.630E+01, 1.760E+01, 1.860E+01,
     & 2.340E+01, 3.200E+01, 4.070E+01, 4.130E+01, 3.200E+01, 3.200E+01,
     & 3.200E+01/


*-----------------------------------------------------------------------
*    proton
*-----------------------------------------------------------------------

       data kmax(2) / 17 /

       data nmax(2) / 33 /
       data (en(2,n),n=1,33) /
     & 3.000E+00, 4.000E+00, 5.000E+00, 6.000E+00, 7.000E+00, 8.000E+00,
     & 9.000E+00, 1.000E+01, 1.100E+01, 1.200E+01, 1.300E+01, 1.400E+01,
     & 1.500E+01, 1.600E+01, 1.700E+01, 2.000E+01, 2.500E+01, 3.000E+01,
     & 3.500E+01, 4.000E+01, 5.000E+01, 6.000E+01, 7.000E+01, 8.000E+01,
     & 9.000E+01, 1.000E+02, 1.500E+02, 3.000E+02, 6.000E+02, 1.000E+03,
     & 3.000E+03, 1.000E+04, 1.000E+10/

*-----------------------------------------------------------------------
*    He-4(p,x)H-3
*-----------------------------------------------------------------------
       data (id(2,n,1),n=1,7) /
     &     2,    4,    1,    1,    1,    3,   19/
       data (xn(2,1,i),i=15,33) /
     & 1.000E-03, 1.000E-01, 2.375E+00, 5.313E+00, 7.500E+00, 9.000E+00,
     & 9.250E+00, 9.188E+00, 9.125E+00, 9.000E+00, 8.750E+00, 8.750E+00,
     & 9.125E+00, 1.100E+01, 1.363E+01, 1.850E+01, 3.000E+01, 3.000E+01,
     & 3.000E+01/

*-----------------------------------------------------------------------
*    N-14(p,x)H-3
*-----------------------------------------------------------------------
       data (id(2,n,2),n=1,7) /
     &     7,   14,    1,    1,    1,    3,   18/
       data (xn(2,2,i),i=16,33) /
     & 1.000E-03, 8.150E-02, 3.020E+00, 3.780E+00, 4.120E+00, 7.330E+00,
     & 7.950E+00, 8.590E+00, 1.100E+01, 1.090E+01, 1.200E+01, 1.350E+01,
     & 1.800E+01, 2.500E+01, 2.720E+01, 2.850E+01, 2.850E+01, 2.850E+01/

*-----------------------------------------------------------------------
*    N-14(p,x)Be-7
*-----------------------------------------------------------------------
       data (id(2,n,3),n=1,7) /
     &     7,   14,    1,    1,    4,    7,   23/
       data (xn(2,3,i),i=11,33) /
     & 1.000E-03, 1.000E+01, 1.700E+01, 2.500E+01, 3.500E+01, 4.500E+01,
     & 3.200E+01, 2.000E+01, 1.300E+01, 1.050E+01, 1.050E+01, 1.150E+01,
     & 1.150E+01, 1.150E+01, 1.100E+01, 1.000E+01, 1.000E+01, 1.000E+01,
     & 1.000E+01, 1.000E+01, 1.000E+01, 1.000E+01, 1.000E+01/

*-----------------------------------------------------------------------
*    N-14(p,x)Be-11
*-----------------------------------------------------------------------
       data (id(2,n,4),n=1,7) /
     &     7,   14,    1,    1,    4,   11,   14/
       data (xn(2,4,i),i=20,33) /
     & 1.000E-03, 3.000E-01, 3.000E-01, 3.000E-01, 3.000E-01, 3.000E-01,
     & 3.000E-01, 2.500E-01, 2.000E-01, 1.790E-01, 9.780E-02, 8.150E-02,
     & 8.150E-02, 8.150E-02/

*-----------------------------------------------------------------------
*    N-14(p,x)C-10
*-----------------------------------------------------------------------
       data (id(2,n,5),n=1,7) /
     &     7,   14,    1,    1,    6,   10,   20/
       data (xn(2,5,i),i=14,33) /
     & 1.000E-03, 3.260E-02, 5.870E-01, 2.920E+00, 6.570E+00, 6.940E+00,
     & 6.760E+00, 4.950E+00, 4.140E+00, 3.880E+00, 3.830E+00, 3.420E+00,
     & 2.900E+00, 2.530E+00, 2.150E+00, 1.920E+00, 1.350E+00, 7.820E-01,
     & 7.820E-01, 7.820E-01/

*-----------------------------------------------------------------------
*    N-14(p,x)C-11
*-----------------------------------------------------------------------
       data (id(2,n,6),n=1,7) /
     &     7,   14,    1,    1,    6,   11,   33/
       data (xn(2,6,i),i=1,33) /
     & 1.000E-03, 1.000E+00, 4.600E+01, 1.370E+02, 1.510E+02, 1.360E+02,
     & 1.250E+02, 1.230E+02, 1.180E+02, 1.160E+02, 1.250E+02, 1.250E+02,
     & 1.160E+02, 1.070E+02, 1.000E+02, 8.000E+01, 6.000E+01, 5.500E+01,
     & 5.000E+01, 4.500E+01, 4.000E+01, 3.300E+01, 2.800E+01, 2.500E+01,
     & 2.200E+01, 2.000E+01, 2.000E+01, 2.000E+01, 2.000E+01, 2.000E+01,
     & 2.000E+01, 2.000E+01, 2.000E+01/

*-----------------------------------------------------------------------
*    N-14(p,x)N-13
*-----------------------------------------------------------------------
       data (id(2,n,7),n=1,7) /
     &     7,   14,    1,    1,    7,   13,   26/
       data (xn(2,7,i),i=8,33) /
     & 1.000E-03, 1.000E+00, 3.000E+00, 6.000E+00, 1.000E+01, 1.500E+01,
     & 2.000E+01, 2.500E+01, 3.200E+01, 3.500E+01, 3.500E+01, 3.500E+01,
     & 3.000E+01, 2.300E+01, 1.800E+01, 1.600E+01, 1.400E+01, 1.300E+01,
     & 1.200E+01, 1.000E+01, 1.000E+01, 1.000E+01, 1.000E+01, 5.000E+00,
     & 5.000E+00, 5.000E+00/

*-----------------------------------------------------------------------
*    N-14(p,x)O-14
*-----------------------------------------------------------------------
       data (id(2,n,8),n=1,7) /
     &     7,   14,    1,    1,    8,   14,   31/
       data (xn(2,8,i),i=3,33) /
     & 1.000E-03, 1.000E+01, 5.000E+01, 1.000E+02, 1.000E+02, 1.000E+02,
     & 1.000E+02, 9.000E+01, 7.500E+01, 6.000E+01, 5.000E+01, 4.000E+01,
     & 3.500E+01, 2.000E+01, 1.000E+01, 6.000E+00, 4.000E+00, 3.000E+00,
     & 1.700E+00, 1.000E+00, 7.000E-01, 5.000E-01, 4.000E-01, 3.000E-01,
     & 1.000E-01, 9.780E-02, 1.140E-01, 2.120E-01, 3.590E-01, 3.590E-01,
     & 3.590E-01/

*-----------------------------------------------------------------------
*    O-16(p,x)H-3
*-----------------------------------------------------------------------
       data (id(2,n,9),n=1,7) /
     &     8,   16,    1,    1,    1,    3,   18/
       data (xn(2,9,i),i=16,33) /
     & 1.000E-03, 1.010E+00, 2.290E+00, 7.660E+00, 8.920E+00, 1.010E+01,
     & 8.610E+00, 9.020E+00, 9.600E+00, 9.330E+00, 9.600E+00, 1.210E+01,
     & 1.470E+01, 4.000E+01, 4.000E+01, 4.000E+01, 4.000E+01, 4.000E+01/

*-----------------------------------------------------------------------
*    O-16(p,x)Be-7
*-----------------------------------------------------------------------
       data (id(2,n,10),n=1,7) /
     &     8,   16,    1,    1,    4,    7,   17/
       data (xn(2,10,i),i=17,33) /
     & 1.000E-03, 1.000E-01, 1.000E+00, 3.000E+00, 7.000E+00, 9.000E+00,
     & 8.500E+00, 8.000E+00, 8.000E+00, 8.000E+00, 8.000E+00, 1.000E+01,
     & 1.000E+01, 1.000E+01, 1.000E+01, 1.000E+01, 1.000E+01/

*-----------------------------------------------------------------------
*    O-16(p,x)Be-11
*-----------------------------------------------------------------------
       data (id(2,n,11),n=1,7) /
     &     8,   16,    1,    1,    4,   11,    8/
       data (xn(2,11,i),i=26,33) /
     & 1.000E-03, 3.410E-02, 5.120E-02, 8.530E-02, 1.190E-01, 8.530E-02,
     & 8.530E-02, 8.530E-02/

*-----------------------------------------------------------------------
*    O-16(p,x)C-10
*-----------------------------------------------------------------------
       data (id(2,n,12),n=1,7) /
     &     8,   16,    1,    1,    6,   10,   16/
       data (xn(2,12,i),i=18,33) /
     & 1.000E-03, 5.000E-01, 1.431E+00, 3.390E+00, 4.920E+00, 6.090E+00,
     & 6.540E+00, 6.540E+00, 7.020E+00, 6.810E+00, 5.880E+00, 5.130E+00,
     & 4.260E+00, 4.000E+00, 4.000E+00, 4.000E+00/

*-----------------------------------------------------------------------
*    O-16(p,x)C-11
*-----------------------------------------------------------------------
       data (id(2,n,13),n=1,7) /
     &     8,   16,    1,    1,    6,   11,   18/
       data (xn(2,13,i),i=16,33) /
     & 1.000E-03, 1.000E+00, 2.000E+00, 3.000E+00, 1.000E+01, 1.800E+01,
     & 1.900E+01, 1.900E+01, 1.800E+01, 1.700E+01, 1.600E+01, 1.500E+01,
     & 1.500E+01, 1.000E+01, 1.000E+01, 1.000E+01, 1.000E+01, 1.000E+01/

*-----------------------------------------------------------------------
*    O-16(p,x)C-14
*-----------------------------------------------------------------------
       data (id(2,n,14),n=1,7) /
     &     8,   16,    1,    1,    6,   14,   17/
       data (xn(2,14,i),i=17,33) /
     & 1.000E-03, 1.000E-01, 5.000E-01, 2.000E+00, 3.000E+00, 3.000E+00,
     & 3.000E+00, 2.800E+00, 2.600E+00, 2.400E+00, 2.000E+00, 2.000E+00,
     & 2.000E+00, 2.000E+00, 2.000E+00, 2.000E+00, 2.000E+00/

*-----------------------------------------------------------------------
*    O-16(p,x)N-13
*-----------------------------------------------------------------------
       data (id(2,n,15),n=1,7) /
     &     8,   16,    1,    1,    7,   13,   30/
       data (xn(2,15,i),i=4,33) /
     & 1.000E-03, 5.270E+00, 8.000E+01, 2.500E+01, 2.000E+01, 7.000E+01,
     & 3.000E+01, 4.000E+01, 5.000E+01, 4.000E+01, 2.400E+01, 1.800E+01,
     & 1.000E+01, 8.000E+00, 8.000E+00, 8.000E+00, 8.000E+00, 8.000E+00,
     & 8.000E+00, 8.000E+00, 8.000E+00, 8.000E+00, 8.000E+00, 8.000E+00,
     & 8.000E+00, 7.000E+00, 7.000E+00, 7.000E+00, 7.000E+00, 7.000E+00/

*-----------------------------------------------------------------------
*    O-16(p,x)O-14
*-----------------------------------------------------------------------
       data (id(2,n,16),n=1,7) /
     &     8,   16,    1,    1,    8,   14,   18/
       data (xn(2,16,i),i=16,33) /
     & 1.000E-03, 4.000E+01, 5.500E+01, 6.500E+01, 7.000E+01, 7.000E+01,
     & 7.000E+01, 7.000E+01, 6.500E+01, 6.300E+01, 6.000E+01, 5.000E+01,
     & 4.200E+01, 3.700E+01, 3.400E+01, 3.000E+01, 3.000E+01, 3.000E+01/

*-----------------------------------------------------------------------
*    O-16(p,x)O-15
*-----------------------------------------------------------------------
       data (id(2,n,17),n=1,7) /
     &     8,   16,    1,    1,    8,   15,   22/
       data (xn(2,17,i),i=12,33) /
     & 1.000E-03, 1.000E+00, 2.000E+00, 4.000E+00, 2.000E+01, 5.500E+01,
     & 7.000E+01, 7.000E+01, 7.000E+01, 7.000E+01, 7.000E+01, 7.000E+01,
     & 7.000E+01, 7.000E+01, 6.000E+01, 5.000E+01, 5.000E+01, 4.000E+01,
     & 4.000E+01, 4.000E+01, 4.000E+01, 4.000E+01/

*-----------------------------------------------------------------------

      end

************************************************************************
      function igetiznm(iz,in,im,m) ! frtati 2022/02/18
*     get or assign numbers for product nuclei in t-yield tally
************************************************************************
      use partmod,only : iznmmx, iznmstore, firstwarning
      include 'param.inc'
      common /tall83/ itnzn(itlmax), itndm(itlmax)

      igetiznm = 1
      it_iznm = 10000*iz + 10*in + im
      do i = 2, iznmmx(m)
        if( it_iznm.eq.iznmstore(i,m) ) then
          igetiznm = i
          return
        end if
      end do
      if( iznmmx(m).eq.itnzn(m) ) then
        igetiznm = itnzn(m)
        if( firstwarning(m) ) then
!$OMP CRITICAL (mxnuclei)
          write(*,'("Warning: Number of product nuclei exceeded",
     &              " mxnuclei in tally ",i2)') m
!$OMP END CRITICAL (mxnuclei)
          firstwarning(m) = .false.
        end if
      else
        iznmmx(m) = iznmmx(m) + 1
        igetiznm = iznmmx(m)
        iznmstore(iznmmx(m),m) = it_iznm
      end if

      return
      end

************************************************************************
      function igetiznmp(iz,in,im,m) ! frtati 2022/02/18
*     get numbers for product nuclei in t-yield tally
************************************************************************
      use partmod,only : iznmmx,  iznmstore
      include 'param.inc'
      common /tall83/ itnzn(itlmax), itndm(itlmax)

      it_iznm = 10000*iz + 10*in + im
      do i = 2, iznmmx(m)
        if( it_iznm.eq.iznmstore(i,m) ) then
          igetiznmp = i
          return
        end if
      end do
      igetiznmp = 1

      return
      end

************************************************************************
      subroutine paraiznm(m,mode) ! frtati 2022/05/02
*     update iznm of each PE for MPI run
************************************************************************
      use partmod
      include 'param.inc'
      common /mpi00/ npe, me
      common /mpi01/ iccp(20000)
      common /tall83/ itnzn(itlmax), itndm(itlmax)

      do ip = 1, npe - 1
        if( iccp(ip).eq.0 ) then
          imast = ip
          exit
        end if
      end do

      if( mode.eq.1 ) then
        if( me.eq.imast ) then
          do ip = 1, npe - 1
            if( ip.ne.imast .and. iccp(ip).eq.0 ) then
              call parari(temp_iznmstore(1,m),itnzn(m),ip)
              iznmadd = 0
              do i = 2, itnzn(m)
                if( temp_iznmstore(i,m).eq.0 ) exit
                do j = 2, iznmmx(m)
                  if( iznmstore(j,m).eq.temp_iznmstore(i,m)  ) exit
                end do
                if( j.eq.iznmmx(m)+1 .and.
     &              iznmmx(m)+iznmadd.lt.itnzn(m)-1 ) then
                  iznmstore(j+iznmadd,m) = temp_iznmstore(i,m)
                  iznmadd = iznmadd + 1
                else if( iznmmx(m)+iznmadd.eq.itnzn(m)-1 .and.
     &                   firstwarning(m) ) then
                  write(*,'("Warning: Number of product nuclei exceeded",
     &                      " mxnuclei in tally ",i2)') m
                  firstwarning(m) = .false.
                end if
              end do
              iznmmx(m) = iznmmx(m) + iznmadd
            end if
          end do
          temp_iznmstore(:,m) = iznmstore(:,m)
          do ip = 1, npe - 1
            if( ip.ne.imast .and. iccp(ip).eq.0 ) then
              call parasi(temp_iznmstore(1,m),itnzn(m),ip)
            end if
          end do
        else
          temp_iznmstore(:,m) = iznmstore(:,m)
          call parasi(temp_iznmstore(1,m),itnzn(m),imast)
          call parari(temp_iznmstore(1,m),itnzn(m),imast)
          do i = 2, itnzn(m)
            if( temp_iznmstore(i,m).eq.0 ) exit
            do j = 2, iznmmx(m)
              if( temp_iznmstore(i,m).eq.iznmstore(j,m) ) then
                iznmturn(i,m) = j
                exit
              end if
            end do
          end do
        end if
        if( me.eq.imast ) then
          do ip = 1, npe - 1
            if( ip.ne.imast .and. iccp(ip).eq.0 ) then
              call parasi(iznmmx(m),1,ip)
            end if
          end do
        else
          call parari(iznmmx(m),1,imast)
          iznmstore(:,m) = temp_iznmstore(:,m)
        end if
      else if( mode.eq.2 ) then
        if( me.eq.0 ) then
          call parari(iznmmx(m),1,imast)
          call parari(iznmstore(1,m),itnzn(m),imast)
        else if( me.eq.imast ) then
          call parasi(iznmmx(m),1,0)
          call parasi(iznmstore(1,m),itnzn(m),0)
        end if
      end if

      return
      end

************************************************************************
      subroutine paraiznm_restart ! frtati 2022/05/02
*     set iznm of each PE to restart MPI run
************************************************************************
      use partmod
      include 'param.inc'
      common /mpi00/ npe, me
      common /tall00/ itnm, ital(itlmax), itals(itlmax), italm(itlmax)
      common /tall83/ itnzn(itlmax), itndm(itlmax)

      do m = 1, itnm
        if( (itals(m).eq.6 .or. itals(m).eq.7 .or. itals(m).eq.8 .or.
     &       itals(m).eq.51) .and. itnzn(m).ne.0 ) then
          if( me.eq.0 ) then
            do ip = 1, npe - 1
              call parasi(iznmmx(m),1,ip)
              call parasi(iznmstore(1,m),itnzn(m),ip)
            end do
          else
            call parari(iznmmx(m),1,0)
            call parari(iznmstore(1,m),itnzn(m),0)
          end if
        end if
      end do

      return
      end
