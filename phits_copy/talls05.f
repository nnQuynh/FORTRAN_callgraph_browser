************************************************************************
*                                                                      *
      subroutine tpdctreg(ncol,m,np,nr,mr,ne,nt,na,nm,mt,nl,lt,
     &                    kr,eb,tb,ab,tr,trEVENT)
*                                                                      *
*       product tally in region mesh                                   *
*       last modified by K.Niita on 2011/05/17                         *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*      ncol  ..... reaction type                                       *
*              4 : source                                              *
*             13 : nuclear reaction (n,x)                              *
*             14 : nuclear reaction (n,n'x)                            *
*                                                                      *
*        jcoll : =  0, nothing happen                                  *
*                =  1, Hydrogen collisions                             *
*                =  2, Particle Decays                                 *
*                =  3, Elastic collisions                              *
*                =  4, Nuclear collisions                              *
*                =  5, Heavy Ion collisions                            *
*                =  6, Neutron reactions by data                       *
*                =  7, Photon reactions by data                        *
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
      use MMBANKMOD !FURUTA
      use partmod, only: itmxpt, itpan, itpat, jtpat ! frtati 2021/10/05
*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param00.inc'
      include 'param.inc'

      parameter ( pi  = 3.141592653589793d0 )

*-----------------------------------------------------------------------

      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
!$OMP THREADPRIVATE(/icomon/)
      common /jcomon/ nabov,nobch,nocas,nomax
!$OMP THREADPRIVATE(/jcomon/)
      common /wtsave/ oldwt
!$OMP THREADPRIVATE(/wtsave/)
      common /tlgeom/ iblz1,iblz2
!$OMP THREADPRIVATE(/tlgeom/)
      common /tlglat/ ilev1,ilev2,ilat1(5,10),ilat2(5,10)
!$OMP THREADPRIVATE(/tlglat/)
      common /mathzn/ mathz, mathn, jcoll, kcoll
!$OMP THREADPRIVATE(/mathzn/)
      common /kmat1d/ idmn(0:kvlmax), idnm(kvmmax)
      common /pnint/  ipnint

*-----------------------------------------------------------------------

      common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                itout(itlmax), ittwo(itlmax)

      common /tall12/ itman(itlmax), itmat(itlmax), itmct(itlmax),
     &                itnun(itlmax), itnuc(itlmax),
     &                itndz(itlmax), itndn(itlmax),
     &                itnkz(itlmax), itnkn(itlmax)
      common /tall22/ itmcn(itlmax), itmtn(itlmax), itmtt(itlmax)
      common /tall37/ itcnt(9,itlmax)
      common /tall48/ itmdp(itlmax,0:30), itmdf(itlmax)

      common /tall23/ itaty(itlmax), itanm(itlmax), itarg(itlmax),
     &                rtami(itlmax), rtama(itlmax), rtadl(itlmax)
      common /tall74/ iprim(itlmax)

      common /tall82/ itcnth(9,itlmax)

      common /tall50/ itlmt(itlmax), ite2l(itlmax)

*-----------------------------------------------------------------------

      common /clustt/ nclsts, iclusts(nnn)
!$OMP THREADPRIVATE(/clustt/)
      common /clustw/ jclusts(0:8,nnn),  qclusts(0:12,nnn)
!$OMP THREADPRIVATE(/clustw/)
      common /cntcls/ jcount(3,nnn)
!$OMP THREADPRIVATE(/cntcls/)

*-----------------------------------------------------------------------

      dimension   kr(mr)
      dimension   mt(nm)
      dimension   lt(nl)
      dimension   eb(ne+1)
      dimension   tb(nt+1)
      dimension   ab(na+1)
      dimension   tr(np,ne,nt,na,nr,2)
      dimension   trEVENT(np,ne,nt,na,nr)       !OBINATA(2012.9.11): as Ct
      real(8),allocatable,save:: tr0(:,:,:,:,:) !OBINATA(2012.9.11): as C

      dimension ips(itmxpt) ! frtati 2021/10/05 6 -> itmxpt
      dimension ncntt(3)
      dimension dmpd(30)

*-----------------------------------------------------------------------

      common /stat / istdev, irestart, ireschk
      common /cparm/ maxbch,maxcas

*-----------------------------------------------------------------------
      integer idmpomp !FURUTA20150427
      common /idmpomp0/idmpomp !FURUTA20150427
*-----------------------------------------------------------------------
      common /tstara/ atmrc(10,8)      ! Ogawa 2023/7/18, (10,7) -> (10,8) for plasmon
!$OMP THREADPRIVATE(/tstara/)
      common /tpdcta/ nflumu
!$OMP THREADPRIVATE(/tpdcta/)
      common /trskip/ ntrskip
!$OMP THREADPRIVATE(/trskip/)

      integer           iMeVperu
      common /cMeVperu/ iMeVperu

      common /cntmx/ ncntmx(3), nctmxsr(3)
!$OMP THREADPRIVATE(/cntmx/)
      common / tsxcl / itsxcl

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
* OBINATA(2012.9.11): change tr(,,,,3) to trEVENT(,,,,)
*-----------------------------------------------------------------------

         if (( ncol .eq. 0 .or. ncol .eq. 4 )
     &                              .and. istdev .eq. 2) then
           if ((nocas.gt.1.or.ncol.eq.0) .and. ihistcount.ne.1 ) then

             tr(:,:,:,:,:,1) = tr(:,:,:,:,:,1) + trEVENT(:,:,:,:,:)
             tr(:,:,:,:,:,2) = tr(:,:,:,:,:,2) + trEVENT(:,:,:,:,:) ** 2

! sumover
             call tpdctreg_sumover(m,1,
     &                   np, ne, nt, na, nr, trEVENT)

           end if

           trEVENT(:,:,:,:,:) = 0

         end if

*-----------------------------------------------------------------------
*        end of batch ( in case of istdev = 1 )
*
* OBINATA(2012.9.11): modificate for thread parallel
*-----------------------------------------------------------------------

         if ( ncol .eq. 0 .and. istdev .eq. 1) then
!$OMP MASTER
             allocate( tr0(np,ne,nt,na,nr) )
             tr0(:,:,:,:,:) = 0.d0
!$OMP END MASTER
!$OMP BARRIER
!$OMP CRITICAL (tpdctreg_crit_ist1)
             tr0(:,:,:,:,:) = tr0(:,:,:,:,:) + trEVENT(:,:,:,:,:)
!$OMP END CRITICAL (tpdctreg_crit_ist1)
!$OMP BARRIER
!$OMP MASTER
             tr(:,:,:,:,:,1) = tr(:,:,:,:,:,1) + tr0(:,:,:,:,:) / maxcas
             tr(:,:,:,:,:,2) = tr(:,:,:,:,:,2)
     &                     + ( tr0(:,:,:,:,:) / maxcas ) ** 2

! sumover
             call tpdctreg_sumover(m,maxcas,
     &                   np, ne, nt, na, nr, tr0)

             deallocate( tr0 )
!$OMP END MASTER

           trEVENT(:,:,:,:,:) = 0

         end if

*-----------------------------------------------------------------------
*        out put : check of ncol, jcoll and kcoll
*-----------------------------------------------------------------------

         if( itout(m) .eq. 1 .and. ncol .ne. 4 ) return
         if( itout(m) .gt. 1 .and.
     &       ncol .ne. 13 .and. ncol .ne. 14 ) return


*-----------------------------------------------------------------------

         if( itout(m) .eq. 2 .and. ! output = nuclear (See subroutine tprodct)
     &       jcoll .ne. 1  .and. jcoll .ne. 3  .and.
     &       jcoll .ne. 4  .and. jcoll .ne. 5  .and.
     &       jcoll .ne. 6  .and. jcoll .ne. 9  .and.      ! ccse 2023/03/29
     &       jcoll .ne. 10 .and. jcoll .ne. 15 .and.
     &       jcoll .ne. 16 .and. jcoll .ne. 17 .and.      ! S.Abe 2017/02/08
     &       jcoll .ne. 19 .and. jcoll .ne. 20 ) return   ! y.sakaki 2024/01

         if( itout(m) .eq. 3 .and. jcoll .ne. 2 ) return
         if( itout(m) .eq. 3 .and.
     &       jcoll .eq. 2 .and. ityp .eq. 13 ) return     ! S.Abe 2017/02/08

         if( itout(m) .eq. 4 .and.
     &       kcoll .ne. 1 .and. kcoll .ne. 5 ) return

         if( itout(m) .eq. 5 .and.
     &       jcoll .ne. 3 .and. kcoll .ne. 3 ) return

         if( itout(m) .eq. 6 .and. ! output = nonela (See subroutine tprodct)
     &       jcoll .ne. 1  .and.
     &       jcoll .ne. 4  .and. jcoll .ne. 5  .and.
     &       jcoll .ne. 6  .and. jcoll .ne. 9  .and.      ! ccse 2023/03/29
     &       jcoll .ne. 10 .and. jcoll .ne. 15 .and.
     &       jcoll .ne. 16 .and. jcoll .ne. 17 .and.
     &       jcoll .ne. 19 .and. jcoll .ne. 20 ) return   ! y.sakaki 2024/01
         if( itout(m) .eq. 6 .and. kcoll .eq. 3 ) return

         sumatmrc = 0.d0
         do i = 1, 5
         do j = 1, 5    ! except mscat
            sumatmrc = sumatmrc + atmrc(i,j)
         enddo
         enddo
         if( itout(m) .eq. 7 .and.
     &       jcoll    .ne. 18 .and.
     &       ( jcoll .ne. 11 .and. sumatmrc .le. 0.d0 )
     &     ) return

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
*        for source
*-----------------------------------------------------------------------

         if( ncol .eq. 4 ) then

               iccol = 0
               npart = nomax

               ncntt(1) = ncnt(ibknct+1,no,ipomp+1)
               ncntt(2) = ncnt(ibknct+2,no,ipomp+1)
               ncntt(3) = ncnt(ibknct+3,no,ipomp+1)

*-----------------------------------------------------------------------
*        for reaction, check of reaction and mother
*-----------------------------------------------------------------------

         else

            if( nclsts .le. 0 ) return

               iccol = 1
               npart = nclsts

               ata   = dble( mathz + mathn )
               atz   = dble( mathz )
               mmas  = nint( ata )
               mchg  = nint( atz )

               ncntt(1) = jcount(1,1)
               ncntt(2) = jcount(2,1)
               ncntt(3) = jcount(3,1)

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

                     if( itmct(m) .gt. 0 ) goto 140
                     if( itmct(m) .lt. 0 ) return

                  end if

               end do

                  if( itmct(m) .gt. 0 ) return

            end if

  140       continue

*-----------------------------------------------------------------------
*           check elastic collision
*-----------------------------------------------------------------------

            if( ( jcoll .eq. 1 .or.
     &            jcoll .eq. 4 .or.
     &            jcoll .eq. 5 .or.
     &            jcoll .eq. 10 ) .and.
     &            npart .eq. 2 ) then

               if( jcoll .eq. 1 ) then

                  if( ( iclusts(1) .eq. 1 .and.
     &                  jclusts(7,2) .eq. ktyp .and.
     &                  qclusts(6,2) .eq. 0.d0 ) .or.
     &                ( iclusts(2) .eq. 1 .and.
     &                  jclusts(7,1) .eq. ktyp .and.
     &                  qclusts(6,1) .eq. 0.d0 ) ) then

                     if( itout(m) .ne. 2 .and. itout(m) .ne. 5 ) return

                  end if

               else

                  if( ( iclusts(1) .eq. 0 .and.
     &                  jclusts(1,1) .eq. mchg .and.
     &                  jclusts(2,1) .eq. mmas - mchg .and.
     &                  jclusts(7,2) .eq. ktyp .and.
     &                  qclusts(6,2) .eq. 0.d0 .and.
     &                  qclusts(6,1) .eq. 0.d0 ) .or.
     &                ( iclusts(2) .eq. 0 .and.
     &                  jclusts(1,2) .eq. mchg .and.
     &                  jclusts(2,2) .eq. mmas - mchg .and.
     &                  jclusts(7,1) .eq. ktyp .and.
     &                  qclusts(6,1) .eq. 0.d0 .and.
     &                  qclusts(6,2) .eq. 0.d0 ) ) then

                     if( itout(m) .ne. 2 .and. itout(m) .ne. 5 ) return

                  end if

               end if

            end if

*-----------------------------------------------------------------------

         end if

*-----------------------------------------------------------------------
*     check of counter
*-----------------------------------------------------------------------

            do i = 1, 3

               if( itcnt(i,m) .ne. 0 ) then

                  if( ncntt(i) .lt. itcnt(i*2+2,m) .or.
     &                ncntt(i) .gt. itcnt(i*2+3,m) ) return

               end if

            end do

*-----------------------------------------------------------------------
*        material for LET
*-----------------------------------------------------------------------

               if( itlmt(m) .gt. 0 ) then

                  lmat = idnm( itlmt(m) )

               else if( itlmt(m) .eq. 0 ) then

                  lmat = mat

               else

                  lmat = -idnm( -itlmt(m) )

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

               ir1 = ii

*-----------------------------------------------------------------------
*        do loop for particles
*-----------------------------------------------------------------------

         do 200 j = 1, npart

*-----------------------------------------------------------------------
            if( ityp .eq. 7 .and. (jcoll.eq.2 .or. jcoll.eq.16) ) then
               if( itout(m) .eq. 7 .and. j .gt. nflumu ) goto 200
               if( itout(m) .ne. 7 .and. j .le. nflumu ) goto 200
            endif
            if( j .eq. ntrskip ) goto 200

            if( iccol .eq. 1 ) then

                  ipart = jclusts(3,j)
                  kpart = jclusts(7,j)
                  lpart = jclusts(8,j)
                  jpart = ichgf(ipart,kpart)
                  rpart = qclusts(5,j) * 1000.d0   ! ccse 2022/08/31

                  epart = qclusts(7,j)
                  tlw   = qclusts(8,j)
                  tpart = abs(qclusts(9,j))
                  xpart = qclusts(10,j)
                  ypart = qclusts(11,j)
                  zpart = qclusts(12,j)

                  upart = qclusts(1,j)
                  vpart = qclusts(2,j)
                  wpart = qclusts(3,j)

                  sxpat = spx(ibkspx+no,ipomp+1)
                  sypat = spy(ibkspy+no,ipomp+1)
                  szpat = spz(ibkspz+no,ipomp+1)
                  nzpat = jclusts(5,j)
                  nname = name(ibknam+no,ipomp+1)

            else

                  nj = no + j - 1

                  ipart = nty(ibknty+nj,ipomp+1)
                  kpart = nkf(ibknkf+nj,ipomp+1)
                  jpart = ichgf(ityp,ktyp)
                  rpart = rtyp   ! ccse 2022/08/31

                  epart = e(ibke+nj,ipomp+1)
                  tlw   = wt(ibkwt+nj,ipomp+1)
                  tpart = abs(t(ibkt+nj,ipomp+1))

                  xpart = x(ibkx+nj,ipomp+1)
                  ypart = y(ibky+nj,ipomp+1)
                  zpart = z(ibkz+nj,ipomp+1)

                  upart = u(ibku+nj,ipomp+1)
                  vpart = v(ibkv+nj,ipomp+1)
                  wpart = w(ibkw+nj,ipomp+1)

                  sxpat = spx(ibkspx+nj,ipomp+1)
                  sypat = spy(ibkspy+nj,ipomp+1)
                  szpat = spz(ibkspz+nj,ipomp+1)
                  nzpat = ctyp
                  nname = name(ibknam+nj,ipomp+1)

            end if

             if( kpart .eq. ktyp .and. iprim(m) .eq. 0 .and.   !exclude primary particle in 
     &      ( (jcoll .eq. 18 .or. jcoll .eq. 14 .or. jcoll .eq. 13 .or. ! Track-structure and EGS
     &         jcoll .eq. 3) .and. j .eq. 1 .or.! elastic scattering
     &         jcoll .eq. 10 .and. j .eq. 1 .and. npart .eq. 2 .or.! elastic scattering
     &         jcoll .ge. 6 .AND. jcoll .le. 9   .and. j .eq. npart !exclude primary sampled from X-section data
     &      ) ) goto 200 ! rejection of remaining projectile
             
            if( epart .le. 0.d0 ) goto 200   ! S.Abe 2017/02/08

*-----------------------------------------------------------------------
*           MeV -> MeV/n energy conversion ! 2018/3/7 Ogawa
*-----------------------------------------------------------------------

               if(iMeVperu.eq.1 .and. ipart.ge.15 .and. ipart.le.19)then
                    ebm = dble(kpart - kpart / 1000000 * 1000000)
               else
                    ebm = 1.d0
               end if

*-----------------------------------------------------------------------
*           check of particles
*-----------------------------------------------------------------------

            call pcheck(m,np,ipart,kpart+1000000000*lpart,jpart,ipn,ips)

               if( ipn .eq. 0 ) goto 200

*-----------------------------------------------------------------------
*           check of energy
*-----------------------------------------------------------------------

            if (ite2l(m) .eq. 0) then   ! ccse 2022/08/31

               if( epart .lt. eb(1)*ebm ) goto 200
               if( epart .ge. eb(ne+1)*ebm ) goto 200

*-----------------------------------------------------------------------
*           energy
*-----------------------------------------------------------------------

            do i = 1, ne

               if( epart .ge. eb(i)*ebm .and.
     &             epart .lt. eb(i+1)*ebm ) goto 30

            end do

            else if (ite2l(m) .eq. 1) then   ! convert to energy to LET
               call dedxas(epart,dedx,lmat,ipart,kpart,jpart,rpart)
               dedx = dedx /10.0d0

               ! check of energy
               if( dedx .lt. eb(1) ) goto 200
               if( dedx .ge. eb(ne+1) ) goto 200

            ! energy
            do i = 1, ne

               if( dedx .ge. eb(i) .and.
     &             dedx .lt. eb(i+1) ) goto 30

            end do

            end if

   30          ie1 = i

*-----------------------------------------------------------------------
*           check of time
*-----------------------------------------------------------------------

               if( tpart .lt. tb(1) ) goto 200
               if( tpart .ge. tb(nt+1) ) goto 200

*-----------------------------------------------------------------------
*           time
*-----------------------------------------------------------------------

            do i = 1, nt

               if( tpart .ge. tb(i) .and.
     &             tpart .lt. tb(i+1) ) goto 40

            end do

   40          it1 = i

*-----------------------------------------------------------------------
*        check angle
*-----------------------------------------------------------------------

         if( itaty(m) .ne. 0 ) then

                  cst = 1.d0
                  pab = sqrt( upart**2 + vpart**2 + wpart**2 )
                  if( pab .gt. 0.0d0 ) cst = wpart / pab

               if( itaty(m) .gt. 0 ) then

                  if( cst .lt. ab(1) ) goto 200
                  if( cst .gt. ab(na+1) ) goto 200

               else

                  if( cst .lt. cos( ab(na+1) / 180.d0 * pi ) ) goto 200
                  if( cst .gt. cos( ab(1) / 180.d0 * pi ) ) goto 200

               end if

            do i = 1, na

               if( itaty(m) .gt. 0 ) then

                  if( cst .ge. ab(i) .and.
     &                cst .le. ab(i+1) ) goto 50

               else

                  if( cst .ge. cos( ab(i+1) / 180.d0 * pi ) .and.
     &                cst .le. cos( ab(i) / 180.d0 * pi ) ) goto 50

               end if

            end do

   50          ia1 = i

         else

               ia1 = 1

         end if

*-----------------------------------------------------------------------
*        tally
*-----------------------------------------------------------------------
            do ip = 1, ipn

               trEVENT(ips(ip),ie1,it1,ia1,ir1) =
     &         trEVENT(ips(ip),ie1,it1,ia1,ir1) + tlw

            end do
*-----------------------------------------------------------------------
*        dump data on file
*-----------------------------------------------------------------------

         if( itmdp(m,0) .ne. 0 ) then
               dmpd(1)  = dble( kpart )
               dmpd(2)  = xpart
               dmpd(3)  = ypart
               dmpd(4)  = zpart
               dmpd(5)  = upart
               dmpd(6)  = vpart
               dmpd(7)  = wpart
               dmpd(8)  = epart
               dmpd(9)  = tlw
               dmpd(10) = tpart
               dmpd(11) = ncntt(1)
               dmpd(12) = ncntt(2)
               dmpd(13) = ncntt(3)
               dmpd(14) = sxpat
               dmpd(15) = sypat
               dmpd(16) = szpat
               dmpd(17) = nname
               dmpd(18) = nocas
               dmpd(19) = nobch
               dmpd(20) = no
               dmpd(21) = nzpat

               if( ipart .ge. 15 )
     &         dmpd(8) = dmpd(8) / ibryf(ipart,kpart)

            if(idmpomp.ne.0)then       !FURUTA20150427
             call writeompfile(m,dmpd) !FURUTA20150427
            else                       !FURUTA20150427

             if( itmdp(m,0) .gt. 0 ) then

               write(itmdf(m))
     &         ( dmpd(itmdp(m,k)), k = 1, abs( itmdp(m,0) ) )

              else

               write(itmdf(m),'(30(1p1d24.15))')
     &         ( dmpd(itmdp(m,k)), k = 1, abs( itmdp(m,0) ) )

              end if

            endif !FURUTA20150427

         end if

*-----------------------------------------------------------------------

  200    continue

*-----------------------------------------------------------------------

  100 continue

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine tpdcttet(ncol,m,np,nr,mr,ne,nt,na,nm,mt,nl,lt,
     &                    kr,eb,tb,ab,tr,trEVENT)
*                                                                      *
*       product tally in tetra mesh                                    *
*       Last Modified by T.Furuta on 2025/01/16                        *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*      ncol  ..... reaction type                                       *
*              4 : source                                              *
*             13 : nuclear reaction (n,x)                              *
*             14 : nuclear reaction (n,n'x)                            *
*                                                                      *
*        jcoll : =  0, nothing happen                                  *
*                =  1, Hydrogen collisions                             *
*                =  2, Particle Decays                                 *
*                =  3, Elastic collisions                              *
*                =  4, Nuclear collisions                              *
*                =  5, Heavy Ion collisions                            *
*                =  6, Neutron reactions by data                       *
*                =  7, Photon reactions by data                        *
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
      use MMBANKMOD !FURUTA
      use partmod, only: itmxpt, itpan, itpat, jtpat ! frtati 2021/10/05
*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param00.inc'
      include 'param.inc'

      parameter ( pi  = 3.141592653589793d0 )

*-----------------------------------------------------------------------

      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
!$OMP THREADPRIVATE(/icomon/)
      common /jcomon/ nabov,nobch,nocas,nomax
!$OMP THREADPRIVATE(/jcomon/)
      common /wtsave/ oldwt
!$OMP THREADPRIVATE(/wtsave/)
      common /tlgeom/ iblz1,iblz2
!$OMP THREADPRIVATE(/tlgeom/)
      common /tlglat/ ilev1,ilev2,ilat1(5,10),ilat2(5,10)
!$OMP THREADPRIVATE(/tlglat/)
      common /mathzn/ mathz, mathn, jcoll, kcoll
!$OMP THREADPRIVATE(/mathzn/)
      common /kmat1d/ idmn(0:kvlmax), idnm(kvmmax)
      common /pnint/  ipnint

*-----------------------------------------------------------------------

      common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                itout(itlmax), ittwo(itlmax)

      common /tall12/ itman(itlmax), itmat(itlmax), itmct(itlmax),
     &                itnun(itlmax), itnuc(itlmax),
     &                itndz(itlmax), itndn(itlmax),
     &                itnkz(itlmax), itnkn(itlmax)
      common /tall22/ itmcn(itlmax), itmtn(itlmax), itmtt(itlmax)
      common /tall37/ itcnt(9,itlmax)
      common /tall48/ itmdp(itlmax,0:30), itmdf(itlmax)

      common /tall23/ itaty(itlmax), itanm(itlmax), itarg(itlmax),
     &                rtami(itlmax), rtama(itlmax), rtadl(itlmax)

      common /tall82/ itcnth(9,itlmax)
      common /tall74/ iprim(itlmax)

      common /tall50/ itlmt(itlmax), ite2l(itlmax)

*-----------------------------------------------------------------------

      common /clustt/ nclsts, iclusts(nnn)
!$OMP THREADPRIVATE(/clustt/)
      common /clustw/ jclusts(0:8,nnn),  qclusts(0:12,nnn)
!$OMP THREADPRIVATE(/clustw/)
      common /cntcls/ jcount(3,nnn)
!$OMP THREADPRIVATE(/cntcls/)

*-----------------------------------------------------------------------

      dimension   kr(mr)
      dimension   mt(nm)
      dimension   lt(nl)
      dimension   eb(ne+1)
      dimension   tb(nt+1)
      dimension   ab(na+1)
      dimension   tr(np,ne,nt,na,nr,2)
      dimension   trEVENT(np,ne,nt,na,nr)       !OBINATA(2012.9.11): as Ct
      real(8),allocatable,save:: tr0(:,:,:,:,:) !OBINATA(2012.9.11): as C

      dimension ips(itmxpt) ! frtati 2021/10/05 6 -> itmxpt
      dimension ncntt(3)
      dimension dmpd(30)

*-----------------------------------------------------------------------

      common /stat / istdev, irestart, ireschk
      common /cparm/ maxbch,maxcas

*-----------------------------------------------------------------------
      integer idmpomp !FURUTA20150427
      common /idmpomp0/idmpomp !FURUTA20150427
*-----------------------------------------------------------------------
      common /tstara/ atmrc(10,8)      ! Ogawa 2023/7/18, (10,7) -> (10,8) for plasmon
!$OMP THREADPRIVATE(/tstara/)
      common /tpdcta/ nflumu
!$OMP THREADPRIVATE(/tpdcta/)
      common /trskip/ ntrskip
!$OMP THREADPRIVATE(/trskip/)

      integer           iMeVperu
      common /cMeVperu/ iMeVperu

      common /regdc/  idrg(kvlmax), idgr(kvmmax)
      integer iii0,kkk0
      common /itettal2/ iii0,kkk0
!$OMP THREADPRIVATE(/itettal2/)

      common /cntmx/ ncntmx(3), nctmxsr(3)
!$OMP THREADPRIVATE(/cntmx/)

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
* OBINATA(2012.9.11): change tr(,,,,3) to trEVENT(,,,,)
*-----------------------------------------------------------------------

         if (( ncol .eq. 0 .or. ncol .eq. 4 )
     &                              .and. istdev .eq. 2) then
           if ((nocas.gt.1.or.ncol.eq.0) .and. ihistcount.ne.1 ) then

             tr(:,:,:,:,:,1) = tr(:,:,:,:,:,1) + trEVENT(:,:,:,:,:)
             tr(:,:,:,:,:,2) = tr(:,:,:,:,:,2) + trEVENT(:,:,:,:,:) ** 2

! sumover
             call tpdctreg_sumover(m,1,
     &                   np,  ne, nt, na, nr, trEVENT)

           end if

           trEVENT(:,:,:,:,:) = 0

         end if

*-----------------------------------------------------------------------
*        end of batch ( in case of istdev = 1 )
*
* OBINATA(2012.9.11): modificate for thread parallel
*-----------------------------------------------------------------------

         if ( ncol .eq. 0 .and. istdev .eq. 1) then
!$OMP MASTER
             allocate( tr0(np,ne,nt,na,nr) )
             tr0(:,:,:,:,:) = 0.d0
!$OMP END MASTER
!$OMP BARRIER
!$OMP CRITICAL (tpdcttet_crit_ist1)
             tr0(:,:,:,:,:) = tr0(:,:,:,:,:) + trEVENT(:,:,:,:,:)
!$OMP END CRITICAL (tpdcttet_crit_ist1)
!$OMP BARRIER
!$OMP MASTER
             tr(:,:,:,:,:,1) = tr(:,:,:,:,:,1) + tr0(:,:,:,:,:) / maxcas
             tr(:,:,:,:,:,2) = tr(:,:,:,:,:,2)
     &                     + ( tr0(:,:,:,:,:) / maxcas ) ** 2

! sumover
             call tpdctreg_sumover(m,maxcas,
     &                   np,  ne, nt, na, nr, tr0)

             deallocate( tr0 )
!$OMP END MASTER

           trEVENT(:,:,:,:,:) = 0

         end if

*-----------------------------------------------------------------------
*        out put : check of ncol, jcoll and kcoll
*-----------------------------------------------------------------------

         if( itout(m) .eq. 1 .and. ncol .ne. 4 ) return
         if( itout(m) .gt. 1 .and.
     &       ncol .ne. 13 .and. ncol .ne. 14 ) return


*-----------------------------------------------------------------------

         if( itout(m) .eq. 2 .and. ! output = nuclear (See subroutine tprodct)
     &       jcoll .ne. 1  .and. jcoll .ne. 3  .and.
     &       jcoll .ne. 4  .and. jcoll .ne. 5  .and.
     &       jcoll .ne. 6  .and. jcoll .ne. 9  .and.      ! ccse 2023/03/29
     &       jcoll .ne. 10 .and. jcoll .ne. 15 .and.
     &       jcoll .ne. 16 .and. jcoll .ne. 17 .and.      ! S.Abe 2017/02/08
     &       jcoll .ne. 19 .and. jcoll .ne. 20 ) return   ! y.sakaki 2024/01

         if( itout(m) .eq. 3 .and. jcoll .ne. 2 ) return
         if( itout(m) .eq. 3 .and.
     &       jcoll .eq. 2 .and. ityp .eq. 13 ) return     ! S.Abe 2017/02/08

         if( itout(m) .eq. 4 .and.
     &       kcoll .ne. 1 .and. kcoll .ne. 5 ) return

         if( itout(m) .eq. 5 .and.
     &       jcoll .ne. 3 .and. kcoll .ne. 3 ) return

         if( itout(m) .eq. 6 .and. ! output = nonela (See subroutine tprodct)
     &       jcoll .ne. 1  .and.
     &       jcoll .ne. 4  .and. jcoll .ne. 5  .and.
     &       jcoll .ne. 6  .and. jcoll .ne. 9  .and.      ! ccse 2023/03/29
     &       jcoll .ne. 10 .and. jcoll .ne. 15 .and.
     &       jcoll .ne. 16 .and. jcoll .ne. 17 .and.
     &       jcoll .ne. 19 .and. jcoll .ne. 20 ) return   ! y.sakaki 2024/01
         if( itout(m) .eq. 6 .and. kcoll .eq. 3 ) return

         sumatmrc = 0.d0
         do i = 1, 5
         do j = 1, 5    ! except mscat
            sumatmrc = sumatmrc + atmrc(i,j)
         enddo
         enddo
         if( itout(m) .eq. 7 .and.
     &       jcoll    .ne. 18 .and.
     &       ( jcoll .ne. 11 .and. sumatmrc .le. 0.d0 )
     &     ) return

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
*        for source
*-----------------------------------------------------------------------

         if( ncol .eq. 4 ) then

               iccol = 0
               npart = nomax

               ncntt(1) = ncnt(ibknct+1,no,ipomp+1)
               ncntt(2) = ncnt(ibknct+2,no,ipomp+1)
               ncntt(3) = ncnt(ibknct+3,no,ipomp+1)

*-----------------------------------------------------------------------
*        for reaction, check of reaction and mother
*-----------------------------------------------------------------------

         else

            if( nclsts .le. 0 ) return

               iccol = 1
               npart = nclsts

               ata   = dble( mathz + mathn )
               atz   = dble( mathz )
               mmas  = nint( ata )
               mchg  = nint( atz )

               ncntt(1) = jcount(1,1)
               ncntt(2) = jcount(2,1)
               ncntt(3) = jcount(3,1)

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

                     if( itmct(m) .gt. 0 ) goto 140
                     if( itmct(m) .lt. 0 ) return

                  end if

               end do

                  if( itmct(m) .gt. 0 ) return

            end if

  140       continue

*-----------------------------------------------------------------------
*           check elastic collision
*-----------------------------------------------------------------------

            if( ( jcoll .eq. 1 .or.
     &            jcoll .eq. 4 .or.
     &            jcoll .eq. 5 .or.
     &            jcoll .eq. 10 ) .and.
     &            npart .eq. 2 ) then

               if( jcoll .eq. 1 ) then

                  if( ( iclusts(1) .eq. 1 .and.
     &                  jclusts(7,2) .eq. ktyp .and.
     &                  qclusts(6,2) .eq. 0.d0 ) .or.
     &                ( iclusts(2) .eq. 1 .and.
     &                  jclusts(7,1) .eq. ktyp .and.
     &                  qclusts(6,1) .eq. 0.d0 ) ) then

                     if( itout(m) .ne. 2 .and. itout(m) .ne. 5 ) return

                  end if

               else

                  if( ( iclusts(1) .eq. 0 .and.
     &                  jclusts(1,1) .eq. mchg .and.
     &                  jclusts(2,1) .eq. mmas - mchg .and.
     &                  jclusts(7,2) .eq. ktyp .and.
     &                  qclusts(6,2) .eq. 0.d0 .and.
     &                  qclusts(6,1) .eq. 0.d0 ) .or.
     &                ( iclusts(2) .eq. 0 .and.
     &                  jclusts(1,2) .eq. mchg .and.
     &                  jclusts(2,2) .eq. mmas - mchg .and.
     &                  jclusts(7,1) .eq. ktyp .and.
     &                  qclusts(6,1) .eq. 0.d0 .and.
     &                  qclusts(6,2) .eq. 0.d0 ) ) then

                     if( itout(m) .ne. 2 .and. itout(m) .ne. 5 ) return

                  end if

               end if

            end if

*-----------------------------------------------------------------------

         end if

*-----------------------------------------------------------------------
*     check of counter
*-----------------------------------------------------------------------

            do i = 1, 3

               if( itcnt(i,m) .ne. 0 ) then

                  if( ncntt(i) .lt. itcnt(i*2+2,m) .or.
     &                ncntt(i) .gt. itcnt(i*2+3,m) ) return

               end if

            end do

*-----------------------------------------------------------------------
*        material for LET
*-----------------------------------------------------------------------

               if( itlmt(m) .gt. 0 ) then

                  lmat = idnm( itlmt(m) )

               else if( itlmt(m) .eq. 0 ) then

                  lmat = mat

               else

                  lmat = -idnm( -itlmt(m) )

               end if

*-----------------------------------------------------------------------
*     check of region
*-----------------------------------------------------------------------

               jj = 0

               ic=idgr(iblz1)
               itet=kkk0-10000
               ihelem=iii0

               call ttetck(ic,mr,kr,itet,ihelem,ir1,icc)

               if( icc .ne. 1 ) goto 100

*-----------------------------------------------------------------------
*        do loop for particles
*-----------------------------------------------------------------------

         do 200 j = 1, npart

*-----------------------------------------------------------------------
            if( ityp .eq. 7 .and. (jcoll.eq.2 .or. jcoll.eq.16) ) then
               if( itout(m) .eq. 7 .and. j .gt. nflumu ) goto 200
               if( itout(m) .ne. 7 .and. j .le. nflumu ) goto 200
            endif
            if( j .eq. ntrskip ) goto 200

            if( iccol .eq. 1 ) then

                  ipart = jclusts(3,j)
                  kpart = jclusts(7,j)
                  lpart = jclusts(8,j)
                  jpart = ichgf(ipart,kpart)
                  rpart = qclusts(5,j) * 1000.d0   ! ccse 2022/08/31

                  epart = qclusts(7,j)
                  tlw   = qclusts(8,j)
                  tpart = abs(qclusts(9,j))
                  xpart = qclusts(10,j)
                  ypart = qclusts(11,j)
                  zpart = qclusts(12,j)

                  upart = qclusts(1,j)
                  vpart = qclusts(2,j)
                  wpart = qclusts(3,j)

                  sxpat = spx(ibkspx+no,ipomp+1)
                  sypat = spy(ibkspy+no,ipomp+1)
                  szpat = spz(ibkspz+no,ipomp+1)
                  nzpat = jclusts(5,j)
                  nname = name(ibknam+no,ipomp+1)

            else

                  nj = no + j - 1

                  ipart = nty(ibknty+nj,ipomp+1)
                  kpart = nkf(ibknkf+nj,ipomp+1)
                  jpart = ichgf(ityp,ktyp)
                  rpart = rtyp   ! ccse 2022/08/31

                  epart = e(ibke+nj,ipomp+1)
                  tlw   = wt(ibkwt+nj,ipomp+1)
                  tpart = abs(t(ibkt+nj,ipomp+1))

                  xpart = x(ibkx+nj,ipomp+1)
                  ypart = y(ibky+nj,ipomp+1)
                  zpart = z(ibkz+nj,ipomp+1)

                  upart = u(ibku+nj,ipomp+1)
                  vpart = v(ibkv+nj,ipomp+1)
                  wpart = w(ibkw+nj,ipomp+1)

                  sxpat = spx(ibkspx+nj,ipomp+1)
                  sypat = spy(ibkspy+nj,ipomp+1)
                  szpat = spz(ibkspz+nj,ipomp+1)
                  nzpat = ctyp
                  nname = name(ibknam+nj,ipomp+1)

            end if

             if( kpart .eq. ktyp .and. iprim(m) .eq. 0 .and.   !exclude primary particle in 
     &      ( (jcoll .eq. 18 .or. jcoll .eq. 14 .or. jcoll .eq. 13 .or. ! Track-structure and EGS
     &         jcoll .eq. 3) .and. j .eq. 1 .or.! elastic scattering
     &         jcoll .eq. 10 .and. j .eq. 1 .and. npart .eq. 2 .or.! elastic scattering
     &         jcoll .ge. 6 .AND. jcoll .le. 9   .and. j .eq. npart !exclude primary sampled from X-section data
     &      ) ) goto 200 ! rejection of remaining projectile
             
            if( epart .le. 0.d0 ) goto 200   ! S.Abe 2017/02/08

*-----------------------------------------------------------------------
*           MeV -> MeV/n energy conversion ! 2018/3/7 Ogawa
*-----------------------------------------------------------------------

               if(iMeVperu.eq.1 .and. ipart.ge.15 .and. ipart.le.19)then
                    ebm = dble(kpart - kpart / 1000000 * 1000000)
               else
                    ebm = 1.d0
               end if

*-----------------------------------------------------------------------
*           check of particles
*-----------------------------------------------------------------------

            call pcheck(m,np,ipart,kpart+1000000000*lpart,jpart,ipn,ips)

               if( ipn .eq. 0 ) goto 200

*-----------------------------------------------------------------------
*           check of energy
*-----------------------------------------------------------------------

            if (ite2l(m) .eq. 0) then   ! ccse 2022/08/31

               if( epart .lt. eb(1)*ebm ) goto 200
               if( epart .ge. eb(ne+1)*ebm ) goto 200

*-----------------------------------------------------------------------
*           energy
*-----------------------------------------------------------------------

            do i = 1, ne

               if( epart .ge. eb(i)*ebm .and.
     &             epart .lt. eb(i+1)*ebm ) goto 30

            end do

            else if (ite2l(m) .eq. 1) then   ! convert to energy to LET
               call dedxas(epart,dedx,lmat,ipart,kpart,jpart,rpart)
               dedx = dedx /10.0d0

               ! check of energy
               if( dedx .lt. eb(1) ) goto 200
               if( dedx .ge. eb(ne+1) ) goto 200

            ! energy
            do i = 1, ne

               if( dedx .ge. eb(i) .and.
     &             dedx .lt. eb(i+1) ) goto 30

            end do

            end if

   30          ie1 = i

*-----------------------------------------------------------------------
*           check of time
*-----------------------------------------------------------------------

               if( tpart .lt. tb(1) ) goto 200
               if( tpart .ge. tb(nt+1) ) goto 200

*-----------------------------------------------------------------------
*           time
*-----------------------------------------------------------------------

            do i = 1, nt

               if( tpart .ge. tb(i) .and.
     &             tpart .lt. tb(i+1) ) goto 40

            end do

   40          it1 = i

*-----------------------------------------------------------------------
*        check angle
*-----------------------------------------------------------------------

         if( itaty(m) .ne. 0 ) then

                  cst = 1.d0
                  pab = sqrt( upart**2 + vpart**2 + wpart**2 )
                  if( pab .gt. 0.0d0 ) cst = wpart / pab

               if( itaty(m) .gt. 0 ) then

                  if( cst .lt. ab(1) ) goto 200
                  if( cst .gt. ab(na+1) ) goto 200

               else

                  if( cst .lt. cos( ab(na+1) / 180.d0 * pi ) ) goto 200
                  if( cst .gt. cos( ab(1) / 180.d0 * pi ) ) goto 200

               end if

            do i = 1, na

               if( itaty(m) .gt. 0 ) then

                  if( cst .ge. ab(i) .and.
     &                cst .le. ab(i+1) ) goto 50

               else

                  if( cst .ge. cos( ab(i+1) / 180.d0 * pi ) .and.
     &                cst .le. cos( ab(i) / 180.d0 * pi ) ) goto 50

               end if

            end do

   50          ia1 = i

         else

               ia1 = 1

         end if

*-----------------------------------------------------------------------
*        tally
*-----------------------------------------------------------------------
            do ip = 1, ipn

               trEVENT(ips(ip),ie1,it1,ia1,ir1) =
     &         trEVENT(ips(ip),ie1,it1,ia1,ir1) + tlw

            end do
*-----------------------------------------------------------------------
*        dump data on file
*-----------------------------------------------------------------------

         if( itmdp(m,0) .ne. 0 ) then
               dmpd(1)  = dble( kpart )
               dmpd(2)  = xpart
               dmpd(3)  = ypart
               dmpd(4)  = zpart
               dmpd(5)  = upart
               dmpd(6)  = vpart
               dmpd(7)  = wpart
               dmpd(8)  = epart
               dmpd(9)  = tlw
               dmpd(10) = tpart
               dmpd(11) = ncntt(1)
               dmpd(12) = ncntt(2)
               dmpd(13) = ncntt(3)
               dmpd(14) = sxpat
               dmpd(15) = sypat
               dmpd(16) = szpat
               dmpd(17) = nname
               dmpd(18) = nocas
               dmpd(19) = nobch
               dmpd(20) = no
               dmpd(21) = nzpat

               if( ipart .ge. 15 )
     &         dmpd(8) = dmpd(8) / ibryf(ipart,kpart)

            if(idmpomp.ne.0)then       !FURUTA20150427
             call writeompfile(m,dmpd) !FURUTA20150427
            else                       !FURUTA20150427

             if( itmdp(m,0) .gt. 0 ) then

               write(itmdf(m))
     &         ( dmpd(itmdp(m,k)), k = 1, abs( itmdp(m,0) ) )

              else

               write(itmdf(m),'(30(1p1d24.15))')
     &         ( dmpd(itmdp(m,k)), k = 1, abs( itmdp(m,0) ) )

              end if

            endif !FURUTA20150427

         end if

*-----------------------------------------------------------------------

  200    continue

*-----------------------------------------------------------------------

  100 continue

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine ppdctreg(m,np,nr,mr,ne,nt,na,kr,eb,tb,ab,tr,
     &                    nvl,ivl,rvl,
     &                    nx,ny,nz,xm,ym,zm,igsh,idasa)
*                                                                      *
*       output the product tally in region mesh                        *
*       last modified by K.Niita on 2005/11/24                         *
*                                                                      *
************************************************************************
      use sangelmod, only: itsans
      use partmod ! frtati 2021/10/05

*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'

      parameter ( pi  = 3.141592653589793d0 )

*-----------------------------------------------------------------------

      common /taliin/ rsouin, nzztin, nrgnin

*-----------------------------------------------------------------------

      common /verjam/ versn, lastr, iyeav, imonv, idayv

      common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                itout(itlmax), ittwo(itlmax)

      common /tall07/ itety(itlmax), itenm(itlmax), iterg(itlmax),
     &                rtemi(itlmax), rtema(itlmax), rtedl(itlmax)

      common /tall20/ ittty(itlmax), ittnm(itlmax), ittrg(itlmax),
     &                rttmi(itlmax), rttma(itlmax), rttdl(itlmax)


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
      common /tall48/ itmdp(itlmax,0:30), itmdf(itlmax)
      common /tall49/ itglt(itlmax)

      common /tall23/ itaty(itlmax), itanm(itlmax), itarg(itlmax),
     &                rtami(itlmax), rtama(itlmax), rtadl(itlmax)

      common /istcut/ ist_cut, ist_bat

      common /fact01/ facmax(itlmax) ! kitamura23/03/31

*-----------------------------------------------------------------------

      common /volreg/ dvol(kvlmax)
      common /regdc/ idrg(kvlmax), idgr(kvmmax)

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

      common /tall50/ itlmt(itlmax), ite2l(itlmax)
*-----------------------------------------------------------------------

      dimension   kr(mr)
      dimension   eb(ne+1)
      dimension   tb(nt+1)
      dimension   ab(na+1)
      dimension   vl(nr)
      dimension   lr(nr)
      dimension   ew(ne)
      dimension   tw(nt)
      dimension   aw(na)
      dimension   tr(np,ne,nt,na,nr,2)
      dimension   ivl(nvl)
      dimension   rvl(nvl)

      dimension   xm(nx+1)
      dimension   ym(ny+1)
      dimension   zm(nz+1)
      dimension   val(nr)

      integer,allocatable :: ixyz(:)

*-----------------------------------------------------------------------

      dimension tott(np,2) ! frtati 2021/10/05 6 -> np

      character hsunit(36)*32

      data hsunit( 1) / '[1/source]                      '/
      data hsunit( 2) / '[1/cm^3/source]                 '/
      data hsunit( 3) / '[1/MeV/source]                  '/
      data hsunit( 4) / '[1/cm^3/MeV/source]             '/
      data hsunit( 5) / '[1/Lethargy/source]             '/
      data hsunit( 6) / '[1/cm^3/Lethargy/source]        '/
      data hsunit(11) / '[1/nsec/source]                 '/
      data hsunit(12) / '[1/cm^3/nsec/source]            '/
      data hsunit(13) / '[1/MeV/nsec/source]             '/
      data hsunit(14) / '[1/cm^3/MeV/nsec/source]        '/
      data hsunit(15) / '[1/Lethargy/nsec/source]        '/
      data hsunit(16) / '[1/cm^3/Lethargy/nsec/source]   '/
      data hsunit(21) / '[1/sr/source]                   '/
      data hsunit(22) / '[1/cm^3/sr/source]              '/
      data hsunit(23) / '[1/MeV/sr/source]               '/
      data hsunit(24) / '[1/cm^3/MeV/sr/source]          '/
      data hsunit(25) / '[1/Lethargy/sr/source]          '/
      data hsunit(26) / '[1/cm^3/Lethargy/sr/source]     '/
      data hsunit(31) / '[1/nsec/sr/source]              '/
      data hsunit(32) / '[1/cm^3/nsec/sr/source]         '/
      data hsunit(33) / '[1/MeV/nsec/sr/source]          '/
      data hsunit(34) / '[1/cm^3/MeV/nsec/sr/source]     '/
      data hsunit(35) / '[1/Lethargy/nsec/sr/source]     '/
      data hsunit(36) / '[1/cm^3/Lethargy/nsec/sr/source]'/

      integer           iMeVperu
      common /cMeVperu/ iMeVperu

      character cha*1
      data cha /"'"/

      character chp(itmxpt)*11  ! kitamura22/03/31
      character chq(itmxpt)*9   ! kitamura22/03/31

      character aname*3

      character rpa*1
      data rpa /'}'/
      character yen*1

      include 'samepage_include/samepage000.inc'

      yen  = char(92)
      allocate( ixyz(max((nx+1)*(ny+1),(ny+1)*(nz+1),(nx+1)*(nz+1))) )

      include 'samepage_include/samepage001.inc'

*-----------------------------------------------------------------------

      if ( iMeVperu.eq. 1 ) then
         hsunit( 3) = '[1/(MeV/n)/source]              '
         hsunit( 4) = '[1/cm^3/(MeV/n)/source]         '
         hsunit(13) = '[1/(MeV/n)/nsec/source]         '
         hsunit(14) = '[1/cm^3/(MeV/n)/nsec/source]    '
         hsunit(23) = '[1/(MeV/n)/sr/source]           '
         hsunit(24) = '[1/cm^3/(MeV/n)/sr/source]      '
         hsunit(33) = '[1/(MeV/n)/nsec/sr/source]      '
         hsunit(34) = '[1/cm^3/(MeV/n)/nsec/sr/source] '
      end if

      if ( ite2l(m) .eq. 1 ) then
         hsunit( 3) = '[1/(keV/um)/source]              '
         hsunit( 4) = '[1/cm^3/(keV/um)/source]         '
         hsunit(13) = '[1/(keV/um)/nsec/source]         '
         hsunit(14) = '[1/cm^3/(keV/um)/nsec/source]    '
         hsunit(23) = '[1/(keV/um)/sr/source]           '
         hsunit(24) = '[1/cm^3/(keV/um)/sr/source]      '
         hsunit(33) = '[1/(keV/um)/nsec/sr/source]      '
         hsunit(34) = '[1/cm^3/(keV/um)/nsec/sr/source] '
      end if

*-----------------------------------------------------------------------
      np_mxang = np
      if ( itmxang(m).ne.0 .and. iabs(itmxang(m)).lt.np ) then
        np_mxang = iabs(itmxang(m))
      end if

*-----------------------------------------------------------------------
*     for igshow
*-----------------------------------------------------------------------

         if( igsh .eq. 0 ) then

            npg = np
            neg = ne
            ntg = nt
            nag = na

         else

            npg = 1
            neg = 1
            ntg = 1
            nag = 1

         end if

*-----------------------------------------------------------------------

            if( itaty(m) .gt. 0 ) then
               aname = 'cos'
            else
               aname = 'the'
            end if

*-----------------------------------------------------------------------
*        particle name
*-----------------------------------------------------------------------

               call wtpname(m,np,chp,chq)

*-----------------------------------------------------------------------
*           itunt(m) = 3 13 23 33 4 14 24 34: /MeV
*-----------------------------------------------------------------------

            if( itunt(m) .eq.  3 .or. itunt(m) .eq.  4 .or.
     &          itunt(m) .eq. 13 .or. itunt(m) .eq. 14 .or.
     &          itunt(m) .eq. 23 .or. itunt(m) .eq. 24 .or.
     &          itunt(m) .eq. 33 .or. itunt(m) .eq. 34 ) then

               do i = 1, ne

                  ew(i) = eb(i+1) - eb(i)

               end do
               ew_sum = eb(ne+1) - eb(1)

            else if( itunt(m) .eq.  5 .or. itunt(m) .eq.  6 .or.
     &               itunt(m) .eq. 15 .or. itunt(m) .eq. 16 .or.
     &               itunt(m) .eq. 25 .or. itunt(m) .eq. 26 .or.
     &               itunt(m) .eq. 35 .or. itunt(m) .eq. 36 ) then

               do i = 1, ne

                  ew(i) = log( eb(i+1) / eb(i) )

               end do
               ew_sum = log( eb(ne+1) / eb(1))

            else

               do i = 1, ne

                  ew(i) = 1.d+0

               end do
               ew_sum = 1.0d0

            end if

*-----------------------------------------------------------------------
*           itunt(m) = 11, 12, 13, 14, 31, 32, 33, 34 : /nsec
*-----------------------------------------------------------------------

            if( ( itunt(m) .ge. 11 .and. itunt(m) .le. 14 ) .or.
     &          ( itunt(m) .ge. 31 .and. itunt(m) .le. 34 ) ) then

               do i = 1, nt

                  tw(i) = tb(i+1) - tb(i)

               end do
               tw_sum = tb(nt+1) - tb(1)

            else

               do i = 1, nt

                  tw(i) = 1.d+0

               end do
               tw_sum = 1.0d0

            end if

*-----------------------------------------------------------------------
*           itunt(m) = 20+, 30+ : /SR
*-----------------------------------------------------------------------

            if( itunt(m) .lt. 20 ) then

               do i = 1, na

                  aw(i) = 1.d+0

               end do
               aw_sum = 1.0d0

            else

               aw_sum = 0.0d0

               do i = 1, na

                  if( itaty(m) .gt. 0 ) then

                     aw(i) = ( ab(i+1) - ab(i) ) * 2.0 * pi
                     aw_sum = aw_sum + aw(i)

                  else

                     aw(i) = ( cos( ab(i) / 180.d0 * pi )
     &                       - cos( ab(i+1) / 180.d0 * pi ) ) * 2.0 * pi
                     aw_sum = aw_sum + aw(i)

                  end if

               end do

            end if

*-----------------------------------------------------------------------
*           itunt(m) = 2, 4, 6, +10,20,30  : /cm^3
*-----------------------------------------------------------------------
*        set mesh volume and name
*-----------------------------------------------------------------------

               call tregvl(nr,mr,kr,vl,lr,nvl,ivl,rvl)

*-----------------------------------------------------------------------
*        relative error and  unit conversion
*        c1 : nomalization for source
*-----------------------------------------------------------------------

         if( igsh .eq. 0 ) then

                  cmax = 0.0
                  cmin = 1.e+33
                  dnon = 1.e-33
                  isdz = 0
                  stdm = 0.0
                  itstd(m) = 1

               if ( rsouin .gt. 0d0 ) then ! S.H. (2016.7.26)

                  c1 = 1.0d+0 / rsouin

               else

                  c1 = 0d0

               end if

            facmax(m) = 1.d0
            if( rtfac(m) .lt. 0.d0 ) then
               facmax(m) = 0.d0

               do 101 ir = 1, nr
               do 101 ie = 1, ne
               do 101 it = 1, nt
               do 101 ia = 1, na
               do 101 ip = 1, np

                  if( tr(ip,ie,it,ia,ir,1) .gt. 0.d0 ) then

                     fmaxfc = tr(ip,ie,it,ia,ir,1)
     &                                   / ew(ie) / tw(it) / aw(ia)

                     if( itunt(m) .eq.  2 .or. itunt(m) .eq.  4 .or.
     &                   itunt(m) .eq.  6 .or. itunt(m) .eq. 12 .or.
     &                   itunt(m) .eq. 14 .or. itunt(m) .eq. 16 .or.
     &                   itunt(m) .eq. 22 .or. itunt(m) .eq. 24 .or.
     &                   itunt(m) .eq. 26 .or. itunt(m) .eq. 32 .or.
     &                   itunt(m) .eq. 34 .or. itunt(m) .eq. 36 )
     &                   fmaxfc = fmaxfc / vl(ir)

                     if( facmax(m) .lt. fmaxfc) facmax(m) = fmaxfc

                  end if

  101          continue

            end if

            do 100 ir = 1, nr
            do 100 ie = 1, ne
            do 100 it = 1, nt
            do 100 ia = 1, na
            do 100 ip = 1, np

               if( tr(ip,ie,it,ia,ir,1) .gt. 0.d0 ) then

                  cc = abs(rtfac(m)/facmax(m))
     &                          / ew(ie) / tw(it) / aw(ia)

                  if( itunt(m) .eq.  2 .or. itunt(m) .eq.  4 .or.
     &                itunt(m) .eq.  6 .or. itunt(m) .eq. 12 .or.
     &                itunt(m) .eq. 14 .or. itunt(m) .eq. 16 .or.
     &                itunt(m) .eq. 22 .or. itunt(m) .eq. 24 .or.
     &                itunt(m) .eq. 26 .or. itunt(m) .eq. 32 .or.
     &                itunt(m) .eq. 34 .or. itunt(m) .eq. 36 )
     &                cc = cc / vl(ir)

                  call calc_stdev(m,Xa,sigx,
     &                            tr(ip,ie,it,ia,ir,1),
     &                            tr(ip,ie,it,ia,ir,2),
     &                            cc)

                  tr(ip,ie,it,ia,ir,1) = Xa
                  tr(ip,ie,it,ia,ir,2) = sigx
                  if( sigx .gt. stdm ) stdm = sigx

                  if( tr(ip,ie,it,ia,ir,1) .gt. cmax )
     &                             cmax = tr(ip,ie,it,ia,ir,1)

                  if( tr(ip,ie,it,ia,ir,1) .lt. cmin )
     &                             cmin = tr(ip,ie,it,ia,ir,1)

               else

                  isdz = 1
                  tr(ip,ie,it,ia,ir,2) = 0.0

               end if

! sumover
               fact_in = abs(rtfac(m)/facmax(m))
               call ppdctreg_sumover_stdev(0,m,ip,ie,it,ia,ir,
     &                fact_in,ew(ie),tw(it),aw(ia),1.0d0,
     &                ew_sum,tw_sum,aw_sum,1.0d0)


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

         if( igsh .ne. 0 .and.
     &     ( itaxs(m,iax) .lt. 7 .or. itrsh(m) .eq. 0 ) ) goto 900

         if( itmdp(m,0) .eq. 0 ) then

            if( (itall .eq. 2 .and. nobch .lt. maxbch .and. igsh .eq. 0)
     &           .or. ( itall .eq. 4 .and. igsh .eq. 0 ) ) then

               call mk_2dnumfn(ctfln(m,iax),fname,itfll(m,iax),
     &              nobch,maxbch,npe)

            else

               fname = ctfln(m,iax)

            end if

         else

            if( (itall .eq. 2 .and. nobch .lt. maxbch .and. igsh .eq. 0)
     &           .or. ( itall .eq. 4 .and. igsh .eq. 0 ) ) then

               call mk_2dnumfn(ctfln(m,iax),fname,itfll(m,iax),
     &              nobch,maxbch,npe)

            else

               fname = ctfln(m,iax)(1:itfll(m,iax))

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

               call tproech(iot,m,iax,1)

*-----------------------------------------------------------------------
*        energy axis  ( LET axis )
*-----------------------------------------------------------------------
      include 'samepage_include/samepage002_petar.inc'
      include 'samepage_include/samepagechp_petar.inc'
      include 'samepage_include/samepageseti.inc'
               inum = 0

         if( itaxs(m,iax) .eq. 1 .or.
     &       itaxs(m,iax) .eq. 15 ) then

            do iri = 1, nr, nrstepi
            do iti = 1, nt, ntstepi
            do iai = 1, na, nastepi
            do ipi = 1, np, npstepi
               ir = iri
               it = iti
               ia = iai
               ip = ipi

               inum = inum + 1

               ireg = lr(ir)

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

               if(iloopmode .ge. 1) then
                  write(iot,'("#   ",a)') chp(ipi)(3:10)
               endif

                  write(iot,'("#   no. =",i3,3x,"reg =",i7)')
     &                     inum, ireg

               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  it =",i3/
     &            "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                  it, tb(it), tb(it+ntstepi)
               end if

               if( itaty(m) .ne. 0 ) then
                  write(iot,'(
     &            "# ia  =",i3,/
     &            "# ",a3," = (",
     &                           1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     ia, aname, ab(ia), ab(ia+nastepi)
               end if

               if( itaxl(m) .eq. 0 ) then

                if ( itaxs(m,iax) .eq. 1 ) then   ! ccse 2022/09/30

                if ( iMeVperu.eq. 1 ) then
                  write(iot,'(/"x: Energy [MeV/n]")')
                else
                  write(iot,'(/"x: Energy [MeV]")')
                end if

                else if( itaxs(m,iax) .eq. 15 ) then
                 write(iot,'(/"x: LET [keV/um]")')
                end if

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: Number ",a32)') hsunit(itunt(m))

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

               if( itety(m) .eq.  3 .or. itety(m) .eq.  5 .or.
     &             itunt(m) .eq.  5 .or. itunt(m) .eq.  6 .or.
     &             itunt(m) .eq. 15 .or. itunt(m) .eq. 16 .or.
     &             itunt(m) .eq. 25 .or. itunt(m) .eq. 26 .or.
     &             itunt(m) .eq. 35 .or. itunt(m) .eq. 36 ) then

                  write(iot,'( "p: xlog ylog afac(0.8) form(0.9)")')

               else

                  write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')

               end if

               if( itanl(m) .gt. 0 ) then

                  write(iot,'( "p: ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               end if

               if( itsans(m) .gt. 0 ) then

                  call write_sangel(iot,m,0)

               end if

       include 'samepage_include/petar_e.inc'

*-----------------------------------------------------------------------

             if(iloopmode .eq. 0) then
               if( ittty(m) .ne. 0 .and. itaty(m) .ne. 0 ) then
                  write(iot,'(/a1,"no. =",i3,
     &                         ",    reg  =",i7,",   t =",
     &                     i3,"   ang =",i3,a1)')
     &                     cha, inum, ireg, it, ia, cha
               else if( ittty(m) .eq. 0 .and. itaty(m) .ne. 0 ) then
                  write(iot,'(/a1,"no. =",i3,
     &                         ",    reg  =",i7,",  ang =",
     &                     i3,a1)')
     &                     cha, inum, ireg, ia, cha
               else if( ittty(m) .ne. 0 .and. itaty(m) .eq. 0 ) then
                  write(iot,'(/a1,"no. =",i3,
     &                         ",    reg  =",i7,",   t =",
     &                     i3,a1)')
     &                     cha, inum, ireg, it, cha
               else
                  write(iot,'(/a1,"no. =",i3,
     &                         ",    reg  =",i7,a1)')
     &                     cha, inum, ireg, cha
               end if
             else
               write(changelsub(1),'(a1,"no. =",i3)') cha,inum
               write(changelsub(2),'(",  part = ",a)') chp(ipi)(3:10)
               write(changelsub(3),'(",  reg =",i7)') ireg
               write(changelsub(4),'(",  e =",i3)') ie
               write(changelsub(5),'(",  t =",i3)') it
               write(changelsub(6),'(",  ang =",i3)') ia
               write(changelsub(7),'(a1)') cha

               if( ne .eq. 0 ) then
                 changelsub(4) = " "
               end if
               if( ittty(m) .eq. 0 ) then
                 changelsub(5) = " "
               end if
               if( itaty(m) .eq. 0 ) then
                 changelsub(6) = " "
               end if

               changelsub(4) = " "

               if(iloopmode .eq. 2) then
                 changelsub(3) = " "
               end if
               if(iloopmode .eq. 1) then
                 changelsub(4) = " "
               end if
               if(iloopmode .eq. 10) then
                 changelsub(5) = " "
               end if
               if(iloopmode .eq. 8 .or. iloopmode .eq. 9) then
                 changelsub(6) = " "
               end if
                angeltitle = trim(changelsub(1))//trim(changelsub(2))//
     &                       trim(changelsub(3))//trim(changelsub(4))//
     &                       trim(changelsub(5))//trim(changelsub(6))//
     &                       trim(changelsub(7))
                write(iot,'(/a)') trim(angeltitle)
          end if

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn


                  write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "   vol  &=&",1pe13.4," [cm^3]")')
     &                     yen, vl(ir)
               if(iloopmode .ge. 1) then
                     write(iot,'("  part  &=& ",a)') chp(ipi)(3:10)
               endif

               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &                     "  tmin  &=&",1pe13.4," [nsec]"/
     &                     "  tmax  &=&",1pe13.4," [nsec]")')
     &                     tb(it), tb(it+ntstepi)
               end if
               if( itaty(m) .ne. 0 ) then
                  write(iot,'(
     &                     "  amin  &=&",1pe13.4,/
     &                     "  amax  &=&",1pe13.4)')
     &                     ab(ia), ab(ia+nastepi)
               end if

                  write(iot,'("e:")')

            end do

            end do
            end do
            end do

*-----------------------------------------------------------------------
*        time axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 11 ) then

            do iri = 1, nr, nrstepi
            do iei = 1, ne, nestepi
            do iai = 1, na, nastepi
            do ipi = 1, np, npstepi
               ir = iri
               ie = iei
               ia = iai
               ip = ipi

               inum = inum + 1

               ireg = lr(ir)

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

               if(iloopmode .ge. 1) then
                  write(iot,'("#   ",a)') chp(ipi)(3:10)
               endif


                  write(iot,'("#   no. =",i3,3x,"reg =",i7/
     &            "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     inum, ireg,
     &                     eb(ie), eb(ie+nestepi)

               if( itaty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  ia =",i3,/
     &            "# ",a3," = (",
     &                           1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     ia, aname, ab(ia), ab(ia+nastepi)
               end if

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: Time [nsec]")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: Number ",a32)') hsunit(itunt(m))

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

               if( ittty(m) .eq. 3 .or. ittty(m) .eq. 5 ) then

                  write(iot,'( "p: xlog ylog afac(0.8) form(0.9)")')

               else

                  write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')

               end if

               if( itanl(m) .gt. 0 ) then

                  write(iot,'( "p: ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               end if

               if( itsans(m) .gt. 0 ) then

                  call write_sangel(iot,m,0)

               end if

       include 'samepage_include/petar_t.inc'

*-----------------------------------------------------------------------

             if(iloopmode .eq. 0) then
               if( itaty(m) .ne. 0 ) then
                  write(iot,'(/a1,"no. =",i3,
     &                         ",    reg  =",i7,",  e =",
     &                     i3,",  ang =",i3,a1)')
     &                     cha, inum, ireg, ie, ia, cha
               else
                  write(iot,'(/a1,"no. =",i3,
     &                         ",    reg  =",i7,",  e =",
     &                     i3,a1)')
     &                     cha, inum, ireg, ie, cha
               end if
             else
               write(changelsub(1),'(a1,"no. =",i3)') cha,inum
               write(changelsub(2),'(",  part = ",a)') chp(ipi)(3:10)
               write(changelsub(3),'(",  reg =",i7)') ireg
               write(changelsub(4),'(",  e =",i3)') ie
               write(changelsub(5),'(",  t =",i3)') it
               write(changelsub(6),'(",  ang =",i3)') ia
               write(changelsub(7),'(a1)') cha

               if( ne .eq. 0 ) then
                 changelsub(4) = " "
               end if
               if( ittty(m) .eq. 0 ) then
                 changelsub(5) = " "
               end if
               if( itaty(m) .eq. 0 ) then
                 changelsub(6) = " "
               end if

               changelsub(5) = " "

               if(iloopmode .eq. 2) then
                 changelsub(3) = " "
               end if
               if(iloopmode .eq. 1) then
                 changelsub(4) = " "
               end if
               if(iloopmode .eq. 10) then
                 changelsub(5) = " "
               end if
               if(iloopmode .eq. 8 .or. iloopmode .eq. 9) then
                 changelsub(6) = " "
               end if
                angeltitle = trim(changelsub(1))//trim(changelsub(2))//
     &                       trim(changelsub(3))//trim(changelsub(4))//
     &                       trim(changelsub(5))//trim(changelsub(6))//
     &                       trim(changelsub(7))
                write(iot,'(/a)') trim(angeltitle)
          end if

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn


                if( ite2l(m).eq.0 ) then  ! ccse 2022/09/30
                  if( iMeVperu.eq.0 ) then
                     write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "   vol  &=&",1pe13.4," [cm^3]"/
     &                     "  emin  &=&",1pe13.4," [MeV]"/
     &                     "  emax  &=&",1pe13.4," [MeV]")')
     &                     yen, vl(ir),
     &                     eb(ie), eb(ie+nestepi)
                  else
                     write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "   vol  &=&",1pe13.4," [cm^3]"/
     &                     "  emin  &=&",1pe13.4," [MeV/n]"/
     &                     "  emax  &=&",1pe13.4," [MeV/n]")')
     &                     yen, vl(ir),
     &                     eb(ie), eb(ie+nestepi)
                  end if
                else if( ite2l(m).eq.1 ) then
                     write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "   vol  &=&",1pe13.4," [cm^3]"/
     &                     "  emin  &=&",1pe13.4," [keV/um]"/
     &                     "  emax  &=&",1pe13.4," [keV/um]")')
     &                     yen, vl(ir),
     &                     eb(ie), eb(ie+nestepi)
                end if
               if(iloopmode .ge. 1) then
                     write(iot,'("  part  &=& ",a)') chp(ipi)(3:10)
               endif
               if( itaty(m) .ne. 0 ) then
                  write(iot,'(
     &                     "  amin  &=&",1pe13.4,/
     &                     "  amax  &=&",1pe13.4)')
     &                     ab(ia), ab(ia+nastepi)
               end if

                  write(iot,'("e:")')

            end do

            end do
            end do
            end do

*-----------------------------------------------------------------------
*        angle axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 12 .or.
     &            itaxs(m,iax) .eq. 13 ) then

            do iri = 1, nr, nrstepi
            do iei = 1, ne, nestepi
            do iti = 1, nt, ntstepi
            do ipi = 1, np, npstepi
               ir = iri
               ie = iei
               it = iti
               ip = ipi

               inum = inum + 1

               ireg = lr(ir)

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

               if(iloopmode .ge. 1) then
                  write(iot,'("#   ",a)') chp(ipi)(3:10)
               endif

                  write(iot,'("#   no. =",i3,3x,"reg =",i7/
     &            "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     inum, ireg,
     &                     eb(ie), eb(ie+nestepi)

               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  it =",i3/
     &            "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                  it, tb(it), tb(it+ntstepi)
               end if

               if( itaxl(m) .eq. 0 ) then

                  if( itaxs(m,iax) .eq. 12 ) then
                     write(iot,'(/"x: cos(",a1,"theta)")') yen
                  else
                     write(iot,'(/"x: ",a1,"theta  [deg]")') yen
                  end if

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: Number ",a32)') hsunit(itunt(m))

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

       include 'samepage_include/petar_a.inc'

*-----------------------------------------------------------------------

             if(iloopmode .eq. 0) then
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(/a1,"no. =",i3,
     &                         ",    reg  =",i7,",  e =",
     &                     i3,",  t =",i3,a1)')
     &                     cha, inum, ireg, ie, it, cha
               else
                  write(iot,'(/a1,"no. =",i3,
     &                         ",    reg  =",i7,",  e =",
     &                     i3,a1)')
     &                     cha, inum, ireg, ie, cha
               end if
             else
               write(changelsub(1),'(a1,"no. =",i3)') cha,inum
               write(changelsub(2),'(",  part = ",a)') chp(ipi)(3:10)
               write(changelsub(3),'(",  reg =",i7)') ireg
               write(changelsub(4),'(",  e =",i3)') ie
               write(changelsub(5),'(",  t =",i3)') it
               write(changelsub(6),'(",  ang =",i3)') ia
               write(changelsub(7),'(a1)') cha

               if( ne .eq. 0 ) then
                 changelsub(4) = " "
               end if
               if( ittty(m) .eq. 0 ) then
                 changelsub(5) = " "
               end if
               if( itaty(m) .eq. 0 ) then
                 changelsub(6) = " "
               end if

               changelsub(6) = " "

               if(iloopmode .eq. 2) then
                 changelsub(3) = " "
               end if
               if(iloopmode .eq. 1) then
                 changelsub(4) = " "
               end if
               if(iloopmode .eq. 10) then
                 changelsub(5) = " "
               end if
               if(iloopmode .eq. 8 .or. iloopmode .eq. 9) then
                 changelsub(6) = " "
               end if
                angeltitle = trim(changelsub(1))//trim(changelsub(2))//
     &                       trim(changelsub(3))//trim(changelsub(4))//
     &                       trim(changelsub(5))//trim(changelsub(6))//
     &                       trim(changelsub(7))
                write(iot,'(/a)') trim(angeltitle)
          end if

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

                if( ite2l(m).eq.0 ) then  ! ccse 2022/09/30
                  if( iMeVperu.eq.0 ) then
                     write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "   vol  &=&",1pe13.4," [cm^3]"/
     &                     "  emin  &=&",1pe13.4," [MeV]"/
     &                     "  emax  &=&",1pe13.4," [MeV]")')
     &                     yen, vl(ir),
     &                     eb(ie), eb(ie+nestepi)
                  else
                     write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "   vol  &=&",1pe13.4," [cm^3]"/
     &                     "  emin  &=&",1pe13.4," [MeV/n]"/
     &                     "  emax  &=&",1pe13.4," [MeV/n]")')
     &                     yen, vl(ir),
     &                     eb(ie), eb(ie+nestepi)
                  end if
                else if( ite2l(m).eq.1 ) then
                     write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "   vol  &=&",1pe13.4," [cm^3]"/
     &                     "  emin  &=&",1pe13.4," [keV/um]"/
     &                     "  emax  &=&",1pe13.4," [keV/um]")')
     &                     yen, vl(ir),
     &                     eb(ie), eb(ie+nestepi)
                end if
               if(iloopmode .ge. 1) then
                     write(iot,'("  part  &=& ",a)') chp(ipi)(3:10)
               endif

               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &                     "  tmin  &=&",1pe13.4," [nsec]"/
     &                     "  tmax  &=&",1pe13.4," [nsec]")')
     &                     tb(it), tb(it+ntstepi)
               end if

                  write(iot,'("e:")')

            end do

            end do
            end do
            end do

*-----------------------------------------------------------------------
*        reg axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 2 ) then

            do iei = 1, ne, nestepi
            do iti = 1, nt, ntstepi
            do iai = 1, na, nastepi
            do ipi = 1, np, npstepi
               ie = iei
               it = iti
               ia = iai
               ip = ipi

               inum = inum + 1

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

               if(iloopmode .ge. 1) then
                  write(iot,'("#   ",a)') chp(ipi)(3:10)
               endif

                  write(iot,'("#   no. =",i3,3x,
     &            "ie  =",i3/
     &            "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        inum, ie,
     &                        eb(ie), eb(ie+nestepi)

               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  it =",i3/
     &            "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                  it, tb(it), tb(it+ntstepi)
               end if

               if( itaty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  ia =",i3,/
     &            "# ",a3," = (",
     &                           1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     ia, aname, ab(ia), ab(ia+nastepi)
               end if

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: Serial Num. of Region")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: Number ",a32)') hsunit(itunt(m))

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

       include 'samepage_include/petar_r.inc'

*-----------------------------------------------------------------------

             if(iloopmode .eq. 0) then
               if( ittty(m) .ne. 0 .and. itaty(m) .ne. 0 ) then
                  write(iot,'(/a1,"no. =",i3,
     &                         ",   e =",i3,
     &                         ",   t =",i3,"   ang =",i3,a1)')
     &                     cha, inum, ie, it, ia, cha
               else if( ittty(m) .eq. 0 .and. itaty(m) .ne. 0 ) then
                  write(iot,'(/a1,"no. =",i3,
     &                         ",   e =",i3,
     &                         ",  ang =",i3,a1)')
     &                     cha, inum, ie, ia, cha
               else if( ittty(m) .ne. 0 .and. itaty(m) .eq. 0 ) then
                  write(iot,'(/a1,"no. =",i3,
     &                         ",   e =",i3,
     &                         ",   t =",i3,a1)')
     &                     cha, inum, ie, it, cha
               else
                  write(iot,'(/a1,"no. =",i3,",  e =",
     &                     i3,a1)')
     &                     cha, inum, ie, cha
               end if
             else
               write(changelsub(1),'(a1,"no. =",i3)') cha,inum
               write(changelsub(2),'(",  part = ",a)') chp(ipi)(3:10)
               write(changelsub(3),'(",  reg =",i7)') ireg
               write(changelsub(4),'(",  e =",i3)') ie
               write(changelsub(5),'(",  t =",i3)') it
               write(changelsub(6),'(",  ang =",i3)') ia
               write(changelsub(7),'(a1)') cha

               if( ne .eq. 0 ) then
                 changelsub(4) = " "
               end if
               if( ittty(m) .eq. 0 ) then
                 changelsub(5) = " "
               end if
               if( itaty(m) .eq. 0 ) then
                 changelsub(6) = " "
               end if

               changelsub(3) = " "

               if(iloopmode .eq. 2) then
                 changelsub(3) = " "
               end if
               if(iloopmode .eq. 1) then
                 changelsub(4) = " "
               end if
               if(iloopmode .eq. 10) then
                 changelsub(5) = " "
               end if
               if(iloopmode .eq. 8 .or. iloopmode .eq. 9) then
                 changelsub(6) = " "
               end if
                angeltitle = trim(changelsub(1))//trim(changelsub(2))//
     &                       trim(changelsub(3))//trim(changelsub(4))//
     &                       trim(changelsub(5))//trim(changelsub(6))//
     &                       trim(changelsub(7))
                write(iot,'(/a)') trim(angeltitle)
          end if


               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

                if( ite2l(m).eq.0 ) then  ! ccse 2022/09/30
                  if( iMeVperu.eq.0 ) then
                     write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "  emin  &=&",1pe13.4," [MeV]"/
     &                     "  emax  &=&",1pe13.4," [MeV]")')
     &                     yen, eb(ie), eb(ie+nestepi)
                  else
                     write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "  emin  &=&",1pe13.4," [MeV/n]"/
     &                     "  emax  &=&",1pe13.4," [MeV/n]")')
     &                     yen, eb(ie), eb(ie+nestepi)
                  end if
                else if( ite2l(m).eq.1 ) then
                     write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "  emin  &=&",1pe13.4," [keV/um]"/
     &                     "  emax  &=&",1pe13.4," [keV/um]")')
     &                     yen, eb(ie), eb(ie+nestepi)
                end if
               if(iloopmode .ge. 1) then
                     write(iot,'("  part  &=& ",a)') chp(ipi)(3:10)
               endif

               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &                     "  tmin  &=&",1pe13.4," [nsec]"/
     &                     "  tmax  &=&",1pe13.4," [nsec]")')
     &                     tb(it), tb(it+ntstepi)
               end if
               if( itaty(m) .ne. 0 ) then
                  write(iot,'(
     &                     "  amin  &=&",1pe13.4,/
     &                     "  amax  &=&",1pe13.4)')
     &                     ab(ia), ab(ia+nastepi)
               end if

                  write(iot,'("e:")')
            end do

            end do
            end do
            end do

*-----------------------------------------------------------------------
*        xy axis (matrix)
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 7 ) then

*-----------------------------------------------------------------------

               inum = 0

            do iz = 1, nz
            do ip = 1, npg
            do ie = 1, neg
            do it = 1, ntg
            do ia = 1, nag

               zval = ( zm(iz) + zm(iz+1) ) / 2.0d0

               inum = inum + 1

                  write(iot,'(/"#",78("-"))')

               if( inum .eq. 1 .or. ip.gt.np_mxang ) then ! frtati 2021/10/05
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

               if( rtfac(m) .lt. 0.d0 )
     &           write(iot,
     &           '( "#restart maximum for normalization:",1p1e14.7)')
     &           facmax(m)

               if ( ip.gt.np_mxang ) then
                 write(iot,'( " SKIPPAGE:")')
                 inum = inum - 1
               end if

                  write(iot,'("#   no. =",i3,3x,
     &            "ie  =",i3,3x,
     &            "iz  =",i3,3x,
     &            "part. = ",a8,/
     &            "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )"/
     &            "#   z = ",1p1e13.4)')
     &                        inum, ie, iz, chq(ip),
     &                        eb(ie), eb(ie+1),
     &                        zval
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  it =",i3/
     &            "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                  it, tb(it), tb(it+1)
               end if

               if( itaty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  ia =",i3,/
     &            "# ",a3," = (",
     &                           1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     ia, aname, ab(ia), ab(ia+1)
               end if

*-----------------------------------------------------------------------

               if( ittty(m) .ne. 0 .and. itaty(m) .ne. 0 ) then
                  write(iot,'(/a1,"no. =",i3,
     &                     ",   e =",i3,",   z =",i3,
     &                     ",   t =",i3,"   ang =",i3,a1)')
     &                     cha, inum, ie, iz, it, ia, cha
               else if( ittty(m) .eq. 0 .and. itaty(m) .ne. 0 ) then
                  write(iot,'(/a1,"no. =",i3,
     &                     ",   e =",i3,",   z =",i3,
     &                     ",  ang =",i3,a1)')
     &                     cha, inum, ie, iz, ia, cha
               else if( ittty(m) .ne. 0 .and. itaty(m) .eq. 0 ) then
                  write(iot,'(/a1,"no. =",i3,
     &                     ",   e =",i3,",   z =",i3,
     &                     ",   t =",i3,a1)')
     &                     cha, inum, ie, iz, it, cha
               else
                  write(iot,'(/a1,"no. =",i3,
     &                     ",    e =",i3,
     &                     ",    z =",i3,a1)')
     &                     cha, inum, ie, iz, cha
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

                  val(i) = tr(ip,ie,it,ia,i,1)

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

         write(iot,'("y: Number ",a32)') hsunit(itunt(m))

      else

         write(iot,'("y: ",200a1)') (itazt(m)(i:i),i=1,itazl(m))

      end if

      if( izlog .eq. 0 ) write(iot,'( "p: ylin")')

      write(iot,'( "#",78("-"))')

      end if

*-----------------------------------------------------------------------

            if( inocm .eq. 1 ) then

                if( ite2l(m).eq.0 ) then  ! ccse 2022/09/30
                  if( iMeVperu.eq.0 ) then
                     write(iot,'(/"wt: s[c5]",/
     &                     a1,"vspace{-3}"/
     &                     "  emin  &=&",1pe13.4," [MeV]"/
     &                     "  emax  &=&",1pe13.4," [MeV]"/
     &                     "  z     &=&",1pe13.4," [cm]"/
     &                     "  part. &=&   ",a8)')
     &                     yen, eb(ie), eb(ie+1),
     &                     zval, chq(ip)
                  else
                     write(iot,'(/"wt: s[c5]",/
     &                     a1,"vspace{-3}"/
     &                     "  emin  &=&",1pe13.4," [MeV/n]"/
     &                     "  emax  &=&",1pe13.4," [MeV/n]"/
     &                     "  z     &=&",1pe13.4," [cm]"/
     &                     "  part. &=&   ",a8)')
     &                     yen, eb(ie), eb(ie+1),
     &                     zval, chq(ip)
                  end if
                else if( ite2l(m).eq.1 ) then
                     write(iot,'(/"wt: s[c5]",/
     &                     a1,"vspace{-3}"/
     &                     "  emin  &=&",1pe13.4," [keV/um]"/
     &                     "  emax  &=&",1pe13.4," [keV/um]"/
     &                     "  z     &=&",1pe13.4," [cm]"/
     &                     "  part. &=&   ",a8)')
     &                     yen, eb(ie), eb(ie+1),
     &                     zval, chq(ip)
                end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &                     "  tmin  &=&",1pe13.4," [nsec]"/
     &                     "  tmax  &=&",1pe13.4," [nsec]")')
     &                     tb(it), tb(it+1)
               end if
               if( itaty(m) .ne. 0 ) then
                  write(iot,'(
     &                     "  amin  &=&",1pe13.4,/
     &                     "  amax  &=&",1pe13.4)')
     &                     ab(ia), ab(ia+1)
               end if

                  write(iot,'("e:")')

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
            end do
            end do
            end do

*-----------------------------------------------------------------------
*        yz axis (matrix)
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 8 ) then

*-----------------------------------------------------------------------

               inum = 0

            do ix = 1, nx
            do ip = 1, npg
            do ie = 1, neg
            do it = 1, ntg
            do ia = 1, nag

               xval = ( xm(ix) + xm(ix+1) ) / 2.0d0

               inum = inum + 1

                  write(iot,'(/"#",78("-"))')

               if( inum .eq. 1 .or. ip.gt.np_mxang ) then ! frtati 2021/10/05
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

               if( rtfac(m) .lt. 0.d0 )
     &           write(iot,
     &           '( "#restart maximum for normalization:",1p1e14.7)')
     &           facmax(m)

               if ( ip.gt.np_mxang ) then
                 write(iot,'( " SKIPPAGE:")')
                 inum = inum - 1
               end if

                  write(iot,'( "#   no. =",i3,3x,
     &            "ie  =",i3,3x,
     &            "ix  =",i3,3x,
     &            "part. = ",a8,/
     &            "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )"/
     &            "#   x = ",1p1e13.4)')
     &                        inum, ie, ix, chq(ip),
     &                        eb(ie), eb(ie+1),
     &                        xval
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  it =",i3/
     &            "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                  it, tb(it), tb(it+1)
               end if

               if( itaty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  ia =",i3,/
     &            "# ",a3," = (",
     &                           1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     ia, aname, ab(ia), ab(ia+1)
               end if

*-----------------------------------------------------------------------

               if( ittty(m) .ne. 0 .and. itaty(m) .ne. 0 ) then
                  write(iot,'(/a1,"no. =",i3,
     &                     ",   e =",i3,",   x =",i3,
     &                     ",   t =",i3,"   ang =",i3,a1)')
     &                     cha, inum, ie, ix, it, ia, cha
               else if( ittty(m) .eq. 0 .and. itaty(m) .ne. 0 ) then
                  write(iot,'(/a1,"no. =",i3,
     &                     ",   e =",i3,",   x =",i3,
     &                     ",  ang =",i3,a1)')
     &                     cha, inum, ie, ix, ia, cha
               else if( ittty(m) .ne. 0 .and. itaty(m) .eq. 0 ) then
                  write(iot,'(/a1,"no. =",i3,
     &                     ",   e =",i3,",   x =",i3,
     &                     ",   t =",i3,a1)')
     &                     cha, inum, ie, ix, it, cha
               else
                  write(iot,'(/a1,"no. =",i3,
     &                     ",    e =",i3,
     &                     ",    x =",i3,a1)')
     &                     cha, inum, ie, ix, cha
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

                  val(i) = tr(ip,ie,it,ia,i,1)

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

         write(iot,'("y: Number ",a32)') hsunit(itunt(m))

      else

         write(iot,'("y: ",200a1)') (itazt(m)(i:i),i=1,itazl(m))

      end if

      if( izlog .eq. 0 ) write(iot,'( "p: ylin")')

      write(iot,'( "#",78("-"))')

      end if

*-----------------------------------------------------------------------

            if( inocm .eq. 1 ) then

                if( ite2l(m).eq.0 ) then  ! ccse 2022/09/30
                  if( iMeVperu.eq.0 ) then
                     write(iot,'(/"wt: s[c5]",/
     &                     a1,"vspace{-3}"/
     &                     "  emin  &=&",1pe13.4," [MeV]"/
     &                     "  emax  &=&",1pe13.4," [MeV]"/
     &                     "  x     &=&",1pe13.4," [cm]"/
     &                     "  part. &=&   ",a8)')
     &                     yen, eb(ie), eb(ie+1),
     &                     xval, chq(ip)
                  else
                     write(iot,'(/"wt: s[c5]",/
     &                     a1,"vspace{-3}"/
     &                     "  emin  &=&",1pe13.4," [MeV/n]"/
     &                     "  emax  &=&",1pe13.4," [MeV/n]"/
     &                     "  x     &=&",1pe13.4," [cm]"/
     &                     "  part. &=&   ",a8)')
     &                     yen, eb(ie), eb(ie+1),
     &                     xval, chq(ip)
                  end if
                else if( ite2l(m).eq.1 ) then
                     write(iot,'(/"wt: s[c5]",/
     &                     a1,"vspace{-3}"/
     &                     "  emin  &=&",1pe13.4," [keV/um]"/
     &                     "  emax  &=&",1pe13.4," [keV/um]"/
     &                     "  x     &=&",1pe13.4," [cm]"/
     &                     "  part. &=&   ",a8)')
     &                     yen, eb(ie), eb(ie+1),
     &                     xval, chq(ip)
                end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &                     "  tmin  &=&",1pe13.4," [nsec]"/
     &                     "  tmax  &=&",1pe13.4," [nsec]")')
     &                     tb(it), tb(it+1)
               end if
               if( itaty(m) .ne. 0 ) then
                  write(iot,'(
     &                     "  amin  &=&",1pe13.4,/
     &                     "  amax  &=&",1pe13.4)')
     &                     ab(ia), ab(ia+1)
               end if

                  write(iot,'("e:")')

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
            end do
            end do
            end do

*-----------------------------------------------------------------------
*        xz axis (matrix)
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 9 ) then

*-----------------------------------------------------------------------

               inum = 0

            do iy = 1, ny
            do ip = 1, npg
            do ie = 1, neg
            do it = 1, ntg
            do ia = 1, nag

               yval = ( ym(iy) + ym(iy+1) ) / 2.0d0

               inum = inum + 1

                  write(iot,'(/"#",78("-"))')

               if( inum .eq. 1 .or. ip.gt.np_mxang ) then ! frtati 2021/10/05
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

               if( rtfac(m) .lt. 0.d0 )
     &           write(iot,
     &           '( "#restart maximum for normalization:",1p1e14.7)')
     &           facmax(m)

               if ( ip.gt.np_mxang ) then
                 write(iot,'( " SKIPPAGE:")')
                 inum = inum - 1
               end if

                  write(iot,'("#   no. =",i3,3x,
     &            "ie  =",i3,3x,
     &            "iy  =",i3,3x,
     &            "part. = ",a8,/
     &            "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )"/
     &            "#   y = ",1p1e13.4)')
     &                        inum, ie, iy, chq(ip),
     &                        eb(ie), eb(ie+1),
     &                        yval
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  it =",i3/
     &            "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                  it, tb(it), tb(it+1)
               end if

               if( itaty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  ia =",i3,/
     &            "# ",a3," = (",
     &                           1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     ia, aname, ab(ia), ab(ia+1)
               end if

*-----------------------------------------------------------------------

               if( ittty(m) .ne. 0 .and. itaty(m) .ne. 0 ) then
                  write(iot,'(/a1,"no. =",i3,
     &                     ",   e =",i3,",   y =",i3,
     &                     ",   t =",i3,"   ang =",i3,a1)')
     &                     cha, inum, ie, iy, it, ia, cha
               else if( ittty(m) .eq. 0 .and. itaty(m) .ne. 0 ) then
                  write(iot,'(/a1,"no. =",i3,
     &                     ",   e =",i3,",   y =",i3,
     &                     ",  ang =",i3,a1)')
     &                     cha, inum, ie, iy, ia, cha
               else if( ittty(m) .ne. 0 .and. itaty(m) .eq. 0 ) then
                  write(iot,'(/a1,"no. =",i3,
     &                     ",   e =",i3,",   y =",i3,
     &                     ",   t =",i3,a1)')
     &                     cha, inum, ie, iy, it, cha
               else
                  write(iot,'(/a1,"no. =",i3,
     &                     ",    e =",i3,
     &                     ",    y =",i3,a1)')
     &                     cha, inum, ie, iy, cha
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

                  val(i) = tr(ip,ie,it,ia,i,1)

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

         write(iot,'("y: Number ",a32)') hsunit(itunt(m))

      else

         write(iot,'("y: ",200a1)') (itazt(m)(i:i),i=1,itazl(m))

      end if

      if( izlog .eq. 0 ) write(iot,'( "p: ylin")')

      write(iot,'( "#",78("-"))')

      end if

*-----------------------------------------------------------------------

            if( inocm .eq. 1 ) then

                if( ite2l(m).eq.0 ) then  ! ccse 2022/09/30
                  if( iMeVperu.eq.0 ) then
                     write(iot,'(/"wt: s[c5]",/
     &                     a1,"vspace{-3}"/
     &                     "  emin  &=&",1pe13.4," [MeV]"/
     &                     "  emax  &=&",1pe13.4," [MeV]"/
     &                     "  y     &=&",1pe13.4," [cm]"/
     &                     "  part. &=&   ",a8/
     &                     "e:")')
     &                     yen, eb(ie), eb(ie+1),
     &                     yval, chq(ip)
                  else
                     write(iot,'(/"wt: s[c5]",/
     &                     a1,"vspace{-3}"/
     &                     "  emin  &=&",1pe13.4," [MeV/n]"/
     &                     "  emax  &=&",1pe13.4," [MeV/n]"/
     &                     "  y     &=&",1pe13.4," [cm]"/
     &                     "  part. &=&   ",a8/
     &                     "e:")')
     &                     yen, eb(ie), eb(ie+1),
     &                     yval, chq(ip)
                  end if
                else if( ite2l(m).eq.1 ) then
                     write(iot,'(/"wt: s[c5]",/
     &                     a1,"vspace{-3}"/
     &                     "  emin  &=&",1pe13.4," [keV/um]"/
     &                     "  emax  &=&",1pe13.4," [keV/um]"/
     &                     "  y     &=&",1pe13.4," [cm]"/
     &                     "  part. &=&   ",a8/
     &                     "e:")')
     &                     yen, eb(ie), eb(ie+1),
     &                     yval, chq(ip)
                end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &                     "  tmin  &=&",1pe13.4," [nsec]"/
     &                     "  tmax  &=&",1pe13.4," [nsec]")')
     &                     tb(it), tb(it+1)
               end if
               if( itaty(m) .ne. 0 ) then
                  write(iot,'(
     &                     "  amin  &=&",1pe13.4,/
     &                     "  amax  &=&",1pe13.4)')
     &                     ab(ia), ab(ia+1)
               end if

                  write(iot,'("e:")')

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
            end do
            end do
            end do

*-----------------------------------------------------------------------

         end if

            call prestart(m,iot) !OBINATA(2012.9.11)

            close(iot)

         if( iteps(m) .ne. 0 ) then

            open(iot, file = fname, status = 'unknown' )
            call a_angel(idasa,fname)

         end if

  900 continue

*-----------------------------------------------------------------------
      deallocate( ixyz )
cnais 2023/01/31
      include 'samepage_include/samepage999.inc'

      return
      end


************************************************************************
*                                                                      *
      subroutine ppdcttet(m,np,nr,mr,ne,nt,na,eb,tb,ab,tr,
     &                    nx,ny,nz,kr,xm,ym,zm,igsh,idasa)
*                                                                      *
*       output the product tally in tetra mesh                         *
*       Last Modified by T.Furuta on 2025/01/16                        *
*                                                                      *
************************************************************************
      use sangelmod, only: itsans
      use partmod ! frtati 2021/10/05

*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'

      parameter ( pi  = 3.141592653589793d0 )

*-----------------------------------------------------------------------

      common /taliin/ rsouin, nzztin, nrgnin

*-----------------------------------------------------------------------

      common /verjam/ versn, lastr, iyeav, imonv, idayv

      common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                itout(itlmax), ittwo(itlmax)

      common /tall07/ itety(itlmax), itenm(itlmax), iterg(itlmax),
     &                rtemi(itlmax), rtema(itlmax), rtedl(itlmax)

      common /tall20/ ittty(itlmax), ittnm(itlmax), ittrg(itlmax),
     &                rttmi(itlmax), rttma(itlmax), rttdl(itlmax)


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
      common /tall48/ itmdp(itlmax,0:30), itmdf(itlmax)
      common /tall49/ itglt(itlmax)

      common /tall23/ itaty(itlmax), itanm(itlmax), itarg(itlmax),
     &                rtami(itlmax), rtama(itlmax), rtadl(itlmax)

      common /istcut/ ist_cut, ist_bat

      common /fact01/ facmax(itlmax) ! kitamura23/03/31

*-----------------------------------------------------------------------

      common /volreg/ dvol(kvlmax)
      common /regdc/ idrg(kvlmax), idgr(kvmmax)

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

      common /tall50/ itlmt(itlmax), ite2l(itlmax)
*-----------------------------------------------------------------------

      dimension   kr(mr)
      dimension   eb(ne+1)
      dimension   tb(nt+1)
      dimension   ab(na+1)
      dimension   ew(ne)
      dimension   tw(nt)
      dimension   aw(na)
      dimension   tr(np,ne,nt,na,nr,2)

      dimension   xm(nx+1)
      dimension   ym(ny+1)
      dimension   zm(nz+1)

      integer,allocatable :: ixyz(:)
      integer,allocatable :: lr(:)        !FURUTA20190204
      real(8),allocatable :: vl(:),val(:) !FURUTA20190204
c Dont know why but necessary to avoid segmentation fault

*-----------------------------------------------------------------------

      dimension tott(np,2) ! frtati 2021/10/05 6 -> np

      character hsunit(36)*32

      data hsunit( 1) / '[1/source]                      '/
      data hsunit( 2) / '[1/cm^3/source]                 '/
      data hsunit( 3) / '[1/MeV/source]                  '/
      data hsunit( 4) / '[1/cm^3/MeV/source]             '/
      data hsunit( 5) / '[1/Lethargy/source]             '/
      data hsunit( 6) / '[1/cm^3/Lethargy/source]        '/
      data hsunit(11) / '[1/nsec/source]                 '/
      data hsunit(12) / '[1/cm^3/nsec/source]            '/
      data hsunit(13) / '[1/MeV/nsec/source]             '/
      data hsunit(14) / '[1/cm^3/MeV/nsec/source]        '/
      data hsunit(15) / '[1/Lethargy/nsec/source]        '/
      data hsunit(16) / '[1/cm^3/Lethargy/nsec/source]   '/
      data hsunit(21) / '[1/sr/source]                   '/
      data hsunit(22) / '[1/cm^3/sr/source]              '/
      data hsunit(23) / '[1/MeV/sr/source]               '/
      data hsunit(24) / '[1/cm^3/MeV/sr/source]          '/
      data hsunit(25) / '[1/Lethargy/sr/source]          '/
      data hsunit(26) / '[1/cm^3/Lethargy/sr/source]     '/
      data hsunit(31) / '[1/nsec/sr/source]              '/
      data hsunit(32) / '[1/cm^3/nsec/sr/source]         '/
      data hsunit(33) / '[1/MeV/nsec/sr/source]          '/
      data hsunit(34) / '[1/cm^3/MeV/nsec/sr/source]     '/
      data hsunit(35) / '[1/Lethargy/nsec/sr/source]     '/
      data hsunit(36) / '[1/cm^3/Lethargy/nsec/sr/source]'/

      integer           iMeVperu
      common /cMeVperu/ iMeVperu

      character cha*1
      data cha /"'"/

      character chp(itmxpt)*11  ! kitamura22/03/31
      character chq(itmxpt)*9   ! kitamura22/03/31

      character aname*3

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
      character yen*1

      include 'samepage_include/samepage000.inc'

      yen  = char(92)
      allocate( ixyz(max((nx+1)*(ny+1),(ny+1)*(nz+1),(nx+1)*(nz+1))) )

      include 'samepage_include/samepage001.inc'



*-----------------------------------------------------------------------

      if ( iMeVperu.eq. 1 ) then
         hsunit( 3) = '[1/(MeV/n)/source]              '
         hsunit( 4) = '[1/cm^3/(MeV/n)/source]         '
         hsunit(13) = '[1/(MeV/n)/nsec/source]         '
         hsunit(14) = '[1/cm^3/(MeV/n)/nsec/source]    '
         hsunit(23) = '[1/(MeV/n)/sr/source]           '
         hsunit(24) = '[1/cm^3/(MeV/n)/sr/source]      '
         hsunit(33) = '[1/(MeV/n)/nsec/sr/source]      '
         hsunit(34) = '[1/cm^3/(MeV/n)/nsec/sr/source] '
      end if

      if ( ite2l(m) .eq. 1 ) then   ! convert to energy to LET
         hsunit( 3) = '[1/(keV/um)/source]              '
         hsunit( 4) = '[1/cm^3/(keV/um)/source]         '
         hsunit(13) = '[1/(keV/um)/nsec/source]         '
         hsunit(14) = '[1/cm^3/(keV/um)/nsec/source]    '
         hsunit(23) = '[1/(keV/um)/sr/source]           '
         hsunit(24) = '[1/cm^3/(keV/um)/sr/source]      '
         hsunit(33) = '[1/(keV/um)/nsec/sr/source]      '
         hsunit(34) = '[1/cm^3/(keV/um)/nsec/sr/source] '
      end if
*-----------------------------------------------------------------------
      np_mxang = np
      if ( itmxang(m).ne.0 .and. iabs(itmxang(m)).lt.np ) then
        np_mxang = iabs(itmxang(m))
      end if

*-----------------------------------------------------------------------
*     for igshow
*-----------------------------------------------------------------------

         if( igsh .eq. 0 ) then

            npg = np
            neg = ne
            ntg = nt
            nag = na

         else

            npg = 1
            neg = 1
            ntg = 1
            nag = 1

         end if

*-----------------------------------------------------------------------

            if( itaty(m) .gt. 0 ) then
               aname = 'cos'
            else
               aname = 'the'
            end if

*-----------------------------------------------------------------------
*        particle name
*-----------------------------------------------------------------------

               call wtpname(m,np,chp,chq)

*-----------------------------------------------------------------------
*           itunt(m) = 3 13 23 33 4 14 24 34: /MeV
*-----------------------------------------------------------------------

            if( itunt(m) .eq.  3 .or. itunt(m) .eq.  4 .or.
     &          itunt(m) .eq. 13 .or. itunt(m) .eq. 14 .or.
     &          itunt(m) .eq. 23 .or. itunt(m) .eq. 24 .or.
     &          itunt(m) .eq. 33 .or. itunt(m) .eq. 34 ) then

               do i = 1, ne

                  ew(i) = eb(i+1) - eb(i)

               end do

            else if( itunt(m) .eq.  5 .or. itunt(m) .eq.  6 .or.
     &               itunt(m) .eq. 15 .or. itunt(m) .eq. 16 .or.
     &               itunt(m) .eq. 25 .or. itunt(m) .eq. 26 .or.
     &               itunt(m) .eq. 35 .or. itunt(m) .eq. 36 ) then

               do i = 1, ne

                  ew(i) = log( eb(i+1) / eb(i) )

               end do

            else

               do i = 1, ne

                  ew(i) = 1.d+0

               end do

            end if

*-----------------------------------------------------------------------
*           itunt(m) = 11, 12, 13, 14, 31, 32, 33, 34 : /nsec
*-----------------------------------------------------------------------

            if( ( itunt(m) .ge. 11 .and. itunt(m) .le. 14 ) .or.
     &          ( itunt(m) .ge. 31 .and. itunt(m) .le. 34 ) ) then

               do i = 1, nt

                  tw(i) = tb(i+1) - tb(i)

               end do

            else

               do i = 1, nt

                  tw(i) = 1.d+0

               end do

            end if

*-----------------------------------------------------------------------
*           itunt(m) = 20+, 30+ : /SR
*-----------------------------------------------------------------------

            if( itunt(m) .lt. 20 ) then

               do i = 1, na

                  aw(i) = 1.d+0

               end do

            else

               do i = 1, na

                  if( itaty(m) .gt. 0 ) then

                     aw(i) = ( ab(i+1) - ab(i) ) * 2.0 * pi

                  else

                     aw(i) = ( cos( ab(i) / 180.d0 * pi )
     &                       - cos( ab(i+1) / 180.d0 * pi ) ) * 2.0 * pi

                  end if

               end do

            end if

*-----------------------------------------------------------------------
*           itunt(m) = 2, 4, 6, +10,20,30  : /cm^3
*-----------------------------------------------------------------------
*        set mesh volume and name
*-----------------------------------------------------------------------

               allocate ( lr(nr),vl(nr),val(nr) ) !FURUTA20190204
               call ttetvl(mr,kr,nr,vl,lr)

*-----------------------------------------------------------------------
*        relative error and  unit conversion
*        c1 : nomalization for source
*-----------------------------------------------------------------------

         if( igsh .eq. 0 ) then

                  cmax = 0.0
                  cmin = 1.e+33
                  dnon = 1.e-33
                  isdz = 0
                  stdm = 0.0
                  itstd(m) = 1

               if ( rsouin .gt. 0d0 ) then ! S.H. (2016.7.26)

                  c1 = 1.0d+0 / rsouin

               else

                  c1 = 0d0

               end if

            facmax(m) = 1.d0
            if( rtfac(m) .lt. 0.d0 ) then
               facmax(m) = 0.d0

               do 101 ir = 1, nr
               do 101 ie = 1, ne
               do 101 it = 1, nt
               do 101 ia = 1, na
               do 101 ip = 1, np

                  if( tr(ip,ie,it,ia,ir,1) .gt. 0.d0 ) then

                     fmaxfc = tr(ip,ie,it,ia,ir,1)
     &                                   / ew(ie) / tw(it) / aw(ia)

                     if( itunt(m) .eq.  2 .or. itunt(m) .eq.  4 .or.
     &                   itunt(m) .eq.  6 .or. itunt(m) .eq. 12 .or.
     &                   itunt(m) .eq. 14 .or. itunt(m) .eq. 16 .or.
     &                   itunt(m) .eq. 22 .or. itunt(m) .eq. 24 .or.
     &                   itunt(m) .eq. 26 .or. itunt(m) .eq. 32 .or.
     &                   itunt(m) .eq. 34 .or. itunt(m) .eq. 36 )
     &                   fmaxfc = fmaxfc / vl(ir)

                     if( facmax(m) .lt. fmaxfc) facmax(m) = fmaxfc

                  end if

  101          continue

            end if

            do 100 ir = 1, nr
            do 100 ie = 1, ne
            do 100 it = 1, nt
            do 100 ia = 1, na
            do 100 ip = 1, np


                  if( itunt(m) .eq.  2 .or. itunt(m) .eq.  4 .or.
     &                itunt(m) .eq.  6 .or. itunt(m) .eq. 12 .or.
     &                itunt(m) .eq. 14 .or. itunt(m) .eq. 16 .or.
     &                itunt(m) .eq. 22 .or. itunt(m) .eq. 24 .or.
     &                itunt(m) .eq. 26 .or. itunt(m) .eq. 32 .or.
     &                itunt(m) .eq. 34 .or. itunt(m) .eq. 36 ) then
                      vw =  vl(ir)
                      vw_sum = sum(vl(:))
                  else
                      vw = 1.0d0
                      vw_sum = 1.0d0
                  endif

               if( tr(ip,ie,it,ia,ir,1) .gt. 0.d0 ) then

                  cc = abs(rtfac(m)/facmax(m))
     &                          / ew(ie) / tw(it) / aw(ia)

                  if( itunt(m) .eq.  2 .or. itunt(m) .eq.  4 .or.
     &                itunt(m) .eq.  6 .or. itunt(m) .eq. 12 .or.
     &                itunt(m) .eq. 14 .or. itunt(m) .eq. 16 .or.
     &                itunt(m) .eq. 22 .or. itunt(m) .eq. 24 .or.
     &                itunt(m) .eq. 26 .or. itunt(m) .eq. 32 .or.
     &                itunt(m) .eq. 34 .or. itunt(m) .eq. 36 )
     &                cc = cc / vl(ir)

                  call calc_stdev(m,Xa,sigx,
     &                            tr(ip,ie,it,ia,ir,1),
     &                            tr(ip,ie,it,ia,ir,2),
     &                            cc)

                  tr(ip,ie,it,ia,ir,1) = Xa
                  tr(ip,ie,it,ia,ir,2) = sigx
                  if( sigx .gt. stdm ) stdm = sigx

                  if( tr(ip,ie,it,ia,ir,1) .gt. cmax )
     &                             cmax = tr(ip,ie,it,ia,ir,1)

                  if( tr(ip,ie,it,ia,ir,1) .lt. cmin )
     &                             cmin = tr(ip,ie,it,ia,ir,1)

               else

                  isdz = 1
                  tr(ip,ie,it,ia,ir,2) = 0.0

               end if

! sumover
               fact_in = abs(rtfac(m)/facmax(m))
               call ppdctreg_sumover_stdev(0,m,ip,ie,it,ia,ir,
     &                fact_in,ew(ie),tw(it),aw(ia),vw,
     &                ew_sum,tw_sum,aw_sum,vw_sum)


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

         if( igsh .ne. 0 .and.
     &     ( itaxs(m,iax) .lt. 7 .or. itrsh(m) .eq. 0 ) ) goto 900

         if( itmdp(m,0) .eq. 0 ) then

            if( (itall .eq. 2 .and. nobch .lt. maxbch .and. igsh .eq. 0)
     &           .or. ( itall .eq. 4 .and. igsh .eq. 0 ) ) then

               call mk_2dnumfn(ctfln(m,iax),fname,itfll(m,iax),
     &              nobch,maxbch,npe)

            else

               fname = ctfln(m,iax)

            end if

         else

            if( (itall .eq. 2 .and. nobch .lt. maxbch .and. igsh .eq. 0)
     &           .or. ( itall .eq. 4 .and. igsh .eq. 0 ) ) then

               call mk_2dnumfn(ctfln(m,iax),fname,itfll(m,iax),
     &              nobch,maxbch,npe)

            else

               fname = ctfln(m,iax)(1:itfll(m,iax))

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

               call tproech(iot,m,iax,1)

*-----------------------------------------------------------------------
*        energy axis  ( LET axis )
*-----------------------------------------------------------------------

      include 'samepage_include/samepage002_petar.inc'
      include 'samepage_include/samepagechp_petar_tet.inc'
      include 'samepage_include/samepageseti.inc'

               inum = 0

         if( itaxs(m,iax) .eq. 1 .or.
     &       itaxs(m,iax) .eq. 15 ) then

            do iri = 1, nr, nrstepi
            do iti = 1, nt, ntstepi
            do iai = 1, na, nastepi
            do ipi = 1, np, npstepi
               ir = iri
               it = iti
               ia = iai
               ip = ipi

               inum = inum + 1

               ireg = lr(ir)

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

               if(iloopmode .ge. 1) then
                  write(iot,'("#   ",a)') chp(ipi)(3:10)
               endif

                  write(iot,'("#   no. =",i3,3x,"tetra=",i8)')
     &                     inum, ireg

               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  it =",i3/
     &            "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                  it, tb(it), tb(it+ntstepi)
               end if

               if( itaty(m) .ne. 0 ) then
                  write(iot,'(
     &            "# ia  =",i3,/
     &            "# ",a3," = (",
     &                           1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     ia, aname, ab(ia), ab(ia+nastepi)
               end if

               if( itaxl(m) .eq. 0 ) then

                if ( itaxs(m,iax) .eq. 1 ) then   ! ccse 2022/09/30

                if ( iMeVperu.eq. 1 ) then
                  write(iot,'(/"x: Energy [MeV/n]")')
                else
                  write(iot,'(/"x: Energy [MeV]")')
                end if

                else if ( itaxs(m,iax) .eq. 15 ) then
                 write(iot,'(/"x: LET [keV/um]")')
                end if

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: Number ",a32)') hsunit(itunt(m))

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

               if( itety(m) .eq.  3 .or. itety(m) .eq.  5 .or.
     &             itunt(m) .eq.  5 .or. itunt(m) .eq.  6 .or.
     &             itunt(m) .eq. 15 .or. itunt(m) .eq. 16 .or.
     &             itunt(m) .eq. 25 .or. itunt(m) .eq. 26 .or.
     &             itunt(m) .eq. 35 .or. itunt(m) .eq. 36 ) then

                  write(iot,'( "p: xlog ylog afac(0.8) form(0.9)")')

               else

                  write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')

               end if

               if( itanl(m) .gt. 0 ) then

                  write(iot,'( "p: ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               end if

               if( itsans(m) .gt. 0 ) then

                  call write_sangel(iot,m,0)

               end if

       include 'samepage_include/petar_e.inc'
*-----------------------------------------------------------------------

             if(iloopmode .eq. 0) then

               if( ittty(m) .ne. 0 .and. itaty(m) .ne. 0 ) then
                  write(iot,'(/a1,"no. =",i3,
     &                         ",   tetra =",i8,",   t =",
     &                     i3,"   ang =",i3,a1)')
     &                     cha, inum, ireg, it, ia, cha
               else if( ittty(m) .eq. 0 .and. itaty(m) .ne. 0 ) then
                  write(iot,'(/a1,"no. =",i3,
     &                         ",   tetra =",i8,",  ang =",
     &                     i3,a1)')
     &                     cha, inum, ireg, ia, cha
               else if( ittty(m) .ne. 0 .and. itaty(m) .eq. 0 ) then
                  write(iot,'(/a1,"no. =",i3,
     &                         ",   tetra =",i8,",   t =",
     &                     i3,a1)')
     &                     cha, inum, ireg, it, cha
               else
                  write(iot,'(/a1,"no. =",i3,
     &                         ",   tetra =",i8,a1)')
     &                     cha, inum, ireg, cha
               end if
             else
               write(changelsub(1),'(a1,"no. =",i3)') cha,inum
               write(changelsub(2),'(",  part = ",a)') chp(ipi)(3:10)
               write(changelsub(3),'(",   tetra =",i8)') ireg
               write(changelsub(4),'(",  e =",i3)') ie
               write(changelsub(5),'(",  t =",i3)') it
               write(changelsub(6),'(",  ang =",i3)') ia
               write(changelsub(7),'(a1)') cha

               if( ne .eq. 0 ) then
                 changelsub(4) = " "
               end if
               if( ittty(m) .eq. 0 ) then
                 changelsub(5) = " "
               end if
               if( itaty(m) .eq. 0 ) then
                 changelsub(6) = " "
               end if

               changelsub(4) = " "

               if(iloopmode .eq. 2) then
                 changelsub(3) = " "
               end if
               if(iloopmode .eq. 1) then
                 changelsub(4) = " "
               end if
               if(iloopmode .eq. 10) then
                 changelsub(5) = " "
               end if
               if(iloopmode .eq. 8 .or. iloopmode .eq. 9) then
                 changelsub(6) = " "
               end if
                angeltitle = trim(changelsub(1))//trim(changelsub(2))//
     &                       trim(changelsub(3))//trim(changelsub(4))//
     &                       trim(changelsub(5))//trim(changelsub(6))//
     &                       trim(changelsub(7))
                write(iot,'(/a)') trim(angeltitle)
          end if

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn


                  write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "   vol  &=&",1pe13.4," [cm^3]")')
     &                     yen, vl(ir)
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &                     "  tmin  &=&",1pe13.4," [nsec]"/
     &                     "  tmax  &=&",1pe13.4," [nsec]")')
     &                     tb(it), tb(it+ntstepi)
               end if
               if(iloopmode .ge. 1) then
                     write(iot,'("  part  &=& ",a)') chp(ipi)(3:10)
               endif
               if( itaty(m) .ne. 0 ) then
                  write(iot,'(
     &                     "  amin  &=&",1pe13.4,/
     &                     "  amax  &=&",1pe13.4)')
     &                     ab(ia), ab(ia+nastepi)
               end if

                  write(iot,'("e:")')

            end do

            end do
            end do
            end do

*-----------------------------------------------------------------------
*        time axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 11 ) then

            do iri = 1, nr, nrstepi
            do iei = 1, ne, nestepi
            do iai = 1, na, nastepi
            do ipi = 1, np, npstepi
               ir = iri
               ie = iei
               ia = iai
               ip = ipi

               inum = inum + 1

               ireg = lr(ir)

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

               if(iloopmode .ge. 1) then
                  write(iot,'("#   ",a)') chp(ipi)(3:10)
               endif


                  write(iot,'("#   no. =",i3,3x,"tetra=",i8/
     &            "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     inum, ireg,
     &                     eb(ie), eb(ie+nestepi)

               if( itaty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  ia =",i3,/
     &            "# ",a3," = (",
     &                           1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     ia, aname, ab(ia), ab(ia+nastepi)
               end if

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: Time [nsec]")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: Number ",a32)') hsunit(itunt(m))

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

               if( ittty(m) .eq. 3 .or. ittty(m) .eq. 5 ) then

                  write(iot,'( "p: xlog ylog afac(0.8) form(0.9)")')

               else

                  write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')

               end if

               if( itanl(m) .gt. 0 ) then

                  write(iot,'( "p: ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               end if

               if( itsans(m) .gt. 0 ) then

                  call write_sangel(iot,m,0)

               end if

       include 'samepage_include/petar_t.inc'

*-----------------------------------------------------------------------

             if(iloopmode .eq. 0) then
               if( itaty(m) .ne. 0 ) then
                  write(iot,'(/a1,"no. =",i3,
     &                         ",   tetra =",i8,",  e =",
     &                     i3,",  ang =",i3,a1)')
     &                     cha, inum, ireg, ie, ia, cha
               else
                  write(iot,'(/a1,"no. =",i3,
     &                         ",   tetra =",i8,",  e =",
     &                     i3,a1)')
     &                     cha, inum, ireg, ie, cha
               end if
             else
               write(changelsub(1),'(a1,"no. =",i3)') cha,inum
               write(changelsub(2),'(",  part = ",a)') chp(ipi)(3:10)
               write(changelsub(3),'(",   tetra =",i8)') ireg
               write(changelsub(4),'(",  e =",i3)') ie
               write(changelsub(5),'(",  t =",i3)') it
               write(changelsub(6),'(",  ang =",i3)') ia
               write(changelsub(7),'(a1)') cha

               if( ne .eq. 0 ) then
                 changelsub(4) = " "
               end if
               if( ittty(m) .eq. 0 ) then
                 changelsub(5) = " "
               end if
               if( itaty(m) .eq. 0 ) then
                 changelsub(6) = " "
               end if

               changelsub(5) = " "

               if(iloopmode .eq. 2) then
                 changelsub(3) = " "
               end if
               if(iloopmode .eq. 1) then
                 changelsub(4) = " "
               end if
               if(iloopmode .eq. 10) then
                 changelsub(5) = " "
               end if
               if(iloopmode .eq. 8 .or. iloopmode .eq. 9) then
                 changelsub(6) = " "
               end if
                angeltitle = trim(changelsub(1))//trim(changelsub(2))//
     &                       trim(changelsub(3))//trim(changelsub(4))//
     &                       trim(changelsub(5))//trim(changelsub(6))//
     &                       trim(changelsub(7))
                write(iot,'(/a)') trim(angeltitle)
          end if

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

                if( ite2l(m).eq.0 ) then  ! ccse 2022/09/30
                  if( iMeVperu.eq.0 ) then
                     write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "   vol  &=&",1pe13.4," [cm^3]"/
     &                     "  emin  &=&",1pe13.4," [MeV]"/
     &                     "  emax  &=&",1pe13.4," [MeV]")')
     &                     yen, vl(ir),
     &                     eb(ie), eb(ie+nestepi)
                  else
                     write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "   vol  &=&",1pe13.4," [cm^3]"/
     &                     "  emin  &=&",1pe13.4," [MeV/n]"/
     &                     "  emax  &=&",1pe13.4," [MeV/n]")')
     &                     yen, vl(ir),
     &                     eb(ie), eb(ie+nestepi)
                  end if
                else if( ite2l(m).eq.1 ) then
                     write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "   vol  &=&",1pe13.4," [cm^3]"/
     &                     "  emin  &=&",1pe13.4," [keV/um]"/
     &                     "  emax  &=&",1pe13.4," [keV/um]")')
     &                     yen, vl(ir),
     &                     eb(ie), eb(ie+nestepi)
                end if
               if(iloopmode .ge. 1) then
                     write(iot,'("  part  &=& ",a)') chp(ipi)(3:10)
               endif

               if( itaty(m) .ne. 0 ) then
                  write(iot,'(
     &                     "  amin  &=&",1pe13.4,/
     &                     "  amax  &=&",1pe13.4)')
     &                     ab(ia), ab(ia+nastepi)
               end if

                  write(iot,'("e:")')

            end do

            end do
            end do
            end do

*-----------------------------------------------------------------------
*        angle axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 12 .or.
     &            itaxs(m,iax) .eq. 13 ) then

            do iri = 1, nr, nrstepi
            do iei = 1, ne, nestepi
            do iti = 1, nt, ntstepi
            do ipi = 1, np, npstepi
               ir = iri
               ie = iei
               it = iti
               ip = ipi

               inum = inum + 1

               ireg = lr(ir)

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

               if(iloopmode .ge. 1) then
                  write(iot,'("#   ",a)') chp(ipi)(3:10)
               endif

                  write(iot,'("#   no. =",i3,3x,"tetra=",i8/
     &            "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     inum, ireg,
     &                     eb(ie), eb(ie+nestepi)

               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  it =",i3/
     &            "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                  it, tb(it), tb(it+ntstepi)
               end if

               if( itaxl(m) .eq. 0 ) then

                  if( itaxs(m,iax) .eq. 12 ) then
                     write(iot,'(/"x: cos(",a1,"theta)")') yen
                  else
                     write(iot,'(/"x: ",a1,"theta  [deg]")') yen
                  end if

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: Number ",a32)') hsunit(itunt(m))

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

       include 'samepage_include/petar_a.inc'

*-----------------------------------------------------------------------

             if(iloopmode .eq. 0) then
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(/a1,"no. =",i3,
     &                         ",   tetra =",i8,",  e =",
     &                     i3,",  t =",i3,a1)')
     &                     cha, inum, ireg, ie, it, cha
               else
                  write(iot,'(/a1,"no. =",i3,
     &                         ",   tetra =",i8,",  e =",
     &                     i3,a1)')
     &                     cha, inum, ireg, ie, cha
               end if
             else
               write(changelsub(1),'(a1,"no. =",i3)') cha,inum
               write(changelsub(2),'(",  part = ",a)') chp(ipi)(3:10)
               write(changelsub(3),'(",   tetra =",i8)') ireg
               write(changelsub(4),'(",  e =",i3)') ie
               write(changelsub(5),'(",  t =",i3)') it
               write(changelsub(6),'(",  ang =",i3)') ia
               write(changelsub(7),'(a1)') cha

               if( ne .eq. 0 ) then
                 changelsub(4) = " "
               end if
               if( ittty(m) .eq. 0 ) then
                 changelsub(5) = " "
               end if
               if( itaty(m) .eq. 0 ) then
                 changelsub(6) = " "
               end if

               changelsub(6) = " "

               if(iloopmode .eq. 2) then
                 changelsub(3) = " "
               end if
               if(iloopmode .eq. 1) then
                 changelsub(4) = " "
               end if
               if(iloopmode .eq. 10) then
                 changelsub(5) = " "
               end if
               if(iloopmode .eq. 8 .or. iloopmode .eq. 9) then
                 changelsub(6) = " "
               end if
                angeltitle = trim(changelsub(1))//trim(changelsub(2))//
     &                       trim(changelsub(3))//trim(changelsub(4))//
     &                       trim(changelsub(5))//trim(changelsub(6))//
     &                       trim(changelsub(7))
                write(iot,'(/a)') trim(angeltitle)
          end if

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn


                if( ite2l(m).eq.0 ) then  ! ccse 2022/09/30
                  if( iMeVperu.eq.0 ) then
                     write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "   vol  &=&",1pe13.4," [cm^3]"/
     &                     "  emin  &=&",1pe13.4," [MeV]"/
     &                     "  emax  &=&",1pe13.4," [MeV]")')
     &                     yen, vl(ir),
     &                     eb(ie), eb(ie+nestepi)
                  else
                     write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "   vol  &=&",1pe13.4," [cm^3]"/
     &                     "  emin  &=&",1pe13.4," [MeV/n]"/
     &                     "  emax  &=&",1pe13.4," [MeV/n]")')
     &                     yen, vl(ir),
     &                     eb(ie), eb(ie+nestepi)
                  end if
                else if( ite2l(m).eq.1 ) then
                     write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "   vol  &=&",1pe13.4," [cm^3]"/
     &                     "  emin  &=&",1pe13.4," [keV/um]"/
     &                     "  emax  &=&",1pe13.4," [keV/um]")')
     &                     yen, vl(ir),
     &                     eb(ie), eb(ie+nestepi)
                end if
               if(iloopmode .ge. 1) then
                     write(iot,'("  part  &=& ",a)') chp(ipi)(3:10)
               endif
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &                     "  tmin  &=&",1pe13.4," [nsec]"/
     &                     "  tmax  &=&",1pe13.4," [nsec]")')
     &                     tb(it), tb(it+ntstepi)
               end if

                  write(iot,'("e:")')

            end do

            end do
            end do
            end do

*-----------------------------------------------------------------------
*        tet axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 14 ) then

            do iei = 1, ne, nestepi
            do iti = 1, nt, ntstepi
            do iai = 1, na, nastepi
            do ipi = 1, np, npstepi
               ie = iei
               it = iti
               ia = iai
               ip = ipi

               inum = inum + 1

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

               if(iloopmode .ge. 1) then
                  write(iot,'("#   ",a)') chp(ipi)(3:10)
               endif

                  write(iot,'("#   no. =",i3,3x,
     &            "ie  =",i3/
     &            "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        inum, ie,
     &                        eb(ie), eb(ie+nestepi)

               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  it =",i3/
     &            "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                  it, tb(it), tb(it+ntstepi)
               end if

               if( itaty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  ia =",i3,/
     &            "# ",a3," = (",
     &                           1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     ia, aname, ab(ia), ab(ia+nastepi)
               end if

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: Serial Num. of Tetrahedron")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: Number ",a32)') hsunit(itunt(m))

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

       include 'samepage_include/petar_r_tetra.inc'

*-----------------------------------------------------------------------

             if(iloopmode .eq. 0) then
               if( ittty(m) .ne. 0 .and. itaty(m) .ne. 0 ) then
                  write(iot,'(/a1,"no. =",i3,
     &                         ",   e =",i3,
     &                         ",   t =",i3,"   ang =",i3,a1)')
     &                     cha, inum, ie, it, ia, cha
               else if( ittty(m) .eq. 0 .and. itaty(m) .ne. 0 ) then
                  write(iot,'(/a1,"no. =",i3,
     &                         ",   e =",i3,
     &                         ",  ang =",i3,a1)')
     &                     cha, inum, ie, ia, cha
               else if( ittty(m) .ne. 0 .and. itaty(m) .eq. 0 ) then
                  write(iot,'(/a1,"no. =",i3,
     &                         ",   e =",i3,
     &                         ",   t =",i3,a1)')
     &                     cha, inum, ie, it, cha
               else
                  write(iot,'(/a1,"no. =",i3,",  e =",
     &                     i3,a1)')
     &                     cha, inum, ie, cha
               end if
             else
               write(changelsub(1),'(a1,"no. =",i3)') cha,inum
               write(changelsub(2),'(",  part = ",a)') chp(ipi)(3:10)
               write(changelsub(3),'(",   tetra =",i8)') ireg
               write(changelsub(4),'(",  e =",i3)') ie
               write(changelsub(5),'(",  t =",i3)') it
               write(changelsub(6),'(",  ang =",i3)') ia
               write(changelsub(7),'(a1)') cha

               if( ne .eq. 0 ) then
                 changelsub(4) = " "
               end if
               if( ittty(m) .eq. 0 ) then
                 changelsub(5) = " "
               end if
               if( itaty(m) .eq. 0 ) then
                 changelsub(6) = " "
               end if

               changelsub(3) = " "

               if(iloopmode .eq. 2) then
                 changelsub(3) = " "
               end if
               if(iloopmode .eq. 1) then
                 changelsub(4) = " "
               end if
               if(iloopmode .eq. 10) then
                 changelsub(5) = " "
               end if
               if(iloopmode .eq. 8 .or. iloopmode .eq. 9) then
                 changelsub(6) = " "
               end if
                angeltitle = trim(changelsub(1))//trim(changelsub(2))//
     &                       trim(changelsub(3))//trim(changelsub(4))//
     &                       trim(changelsub(5))//trim(changelsub(6))//
     &                       trim(changelsub(7))
                write(iot,'(/a)') trim(angeltitle)
          end if


               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

                if( ite2l(m).eq.0 ) then  ! ccse 2022/09/30
                  if( iMeVperu.eq.0 ) then
                     write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "  emin  &=&",1pe13.4," [MeV]"/
     &                     "  emax  &=&",1pe13.4," [MeV]")')
     &                     yen, eb(ie), eb(ie+nestepi)
                  else
                     write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "  emin  &=&",1pe13.4," [MeV/n]"/
     &                     "  emax  &=&",1pe13.4," [MeV/n]")')
     &                     yen, eb(ie), eb(ie+nestepi)
                  end if
                else if( ite2l(m).eq.1 ) then
                     write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "  emin  &=&",1pe13.4," [keV/um]"/
     &                     "  emax  &=&",1pe13.4," [keV/um]")')
     &                     yen, eb(ie), eb(ie+nestepi)
                end if
               if(iloopmode .ge. 1) then
                     write(iot,'("  part  &=& ",a)') chp(ipi)(3:10)
               endif

               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &                     "  tmin  &=&",1pe13.4," [nsec]"/
     &                     "  tmax  &=&",1pe13.4," [nsec]")')
     &                     tb(it), tb(it+ntstepi)
               end if
               if( itaty(m) .ne. 0 ) then
                  write(iot,'(
     &                     "  amin  &=&",1pe13.4,/
     &                     "  amax  &=&",1pe13.4)')
     &                     ab(ia), ab(ia+nastepi)
               end if

                  write(iot,'("e:")')

            end do

            end do
            end do
            end do

cFURUTA20190208 OpenFOAM output
            if ( itfoam(m) .ne. 0) then
             numIndex = 4
             allocate( foamfIType(1:numIndex) )
             allocate( foamfIndex(1:numIndex) )
             foamfIType = (/ 'a', 't', 'e', 'p' /)
             foamfIndex = (/ na, nt, ne, np/) !FURUTA20191028 bugfix
            endif
cFURUTA20191028 CSV output
            if ( itfoam(m) .eq. 2) then
             allocate ( xcm(3,nr) ) !FURUTA20191028
             call ttetcm(mr,kr,nr,xcm)
            endif

*-----------------------------------------------------------------------
*        xy axis (matrix)
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 7 ) then

*-----------------------------------------------------------------------

               inum = 0

            do iz = 1, nz
            do ip = 1, npg
            do ie = 1, neg
            do it = 1, ntg
            do ia = 1, nag

               zval = ( zm(iz) + zm(iz+1) ) / 2.0d0

               inum = inum + 1

                  write(iot,'(/"#",78("-"))')

               if( inum .eq. 1 .or. ip.gt.np_mxang ) then ! frtati 2021/10/05
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

               if( rtfac(m) .lt. 0.d0 )
     &           write(iot,
     &           '( "#restart maximum for normalization:",1p1e14.7)')
     &           facmax(m)

               if ( ip.gt.np_mxang ) then
                 write(iot,'( " SKIPPAGE:")')
                 inum = inum - 1
               end if

                  write(iot,'("#   no. =",i3,3x,
     &            "ie  =",i3,3x,
     &            "iz  =",i3,3x,
     &            "part. = ",a8,/
     &            "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )"/
     &            "#   z = ",1p1e13.4)')
     &                        inum, ie, iz, chq(ip),
     &                        eb(ie), eb(ie+1),
     &                        zval
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  it =",i3/
     &            "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                  it, tb(it), tb(it+1)
               end if

               if( itaty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  ia =",i3,/
     &            "# ",a3," = (",
     &                           1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     ia, aname, ab(ia), ab(ia+1)
               end if

*-----------------------------------------------------------------------

               if( ittty(m) .ne. 0 .and. itaty(m) .ne. 0 ) then
                  write(iot,'(/a1,"no. =",i3,
     &                     ",   e =",i3,",   z =",i3,
     &                     ",   t =",i3,"   ang =",i3,a1)')
     &                     cha, inum, ie, iz, it, ia, cha
               else if( ittty(m) .eq. 0 .and. itaty(m) .ne. 0 ) then
                  write(iot,'(/a1,"no. =",i3,
     &                     ",   e =",i3,",   z =",i3,
     &                     ",  ang =",i3,a1)')
     &                     cha, inum, ie, iz, ia, cha
               else if( ittty(m) .ne. 0 .and. itaty(m) .eq. 0 ) then
                  write(iot,'(/a1,"no. =",i3,
     &                     ",   e =",i3,",   z =",i3,
     &                     ",   t =",i3,a1)')
     &                     cha, inum, ie, iz, it, cha
               else
                  write(iot,'(/a1,"no. =",i3,
     &                     ",    e =",i3,
     &                     ",    z =",i3,a1)')
     &                     cha, inum, ie, iz, cha
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

                  val(i) = tr(ip,ie,it,ia,i,1)

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

         write(iot,'("y: Number ",a32)') hsunit(itunt(m))

      else

         write(iot,'("y: ",200a1)') (itazt(m)(i:i),i=1,itazl(m))

      end if

      if( izlog .eq. 0 ) write(iot,'( "p: ylin")')

      write(iot,'( "#",78("-"))')

      end if

*-----------------------------------------------------------------------

            if( inocm .eq. 1 ) then

                if( ite2l(m).eq.0 ) then  ! ccse 2022/09/30
                  if( iMeVperu.eq.0 ) then
                     write(iot,'(/"wt: s[c5]",/
     &                     a1,"vspace{-3}"/
     &                     "  emin  &=&",1pe13.4," [MeV]"/
     &                     "  emax  &=&",1pe13.4," [MeV]"/
     &                     "  z     &=&",1pe13.4," [cm]"/
     &                     "  part. &=&   ",a8)')
     &                     yen, eb(ie), eb(ie+1),
     &                     zval, chq(ip)
                  else
                     write(iot,'(/"wt: s[c5]",/
     &                     a1,"vspace{-3}"/
     &                     "  emin  &=&",1pe13.4," [MeV/n]"/
     &                     "  emax  &=&",1pe13.4," [MeV/n]"/
     &                     "  z     &=&",1pe13.4," [cm]"/
     &                     "  part. &=&   ",a8)')
     &                     yen, eb(ie), eb(ie+1),
     &                     zval, chq(ip)
                  end if
                else if( ite2l(m).eq.1 ) then
                     write(iot,'(/"wt: s[c5]",/
     &                     a1,"vspace{-3}"/
     &                     "  emin  &=&",1pe13.4," [keV/um]"/
     &                     "  emax  &=&",1pe13.4," [keV/um]"/
     &                     "  z     &=&",1pe13.4," [cm]"/
     &                     "  part. &=&   ",a8)')
     &                     yen, eb(ie), eb(ie+1),
     &                     zval, chq(ip)
                end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &                     "  tmin  &=&",1pe13.4," [nsec]"/
     &                     "  tmax  &=&",1pe13.4," [nsec]")')
     &                     tb(it), tb(it+1)
               end if
               if( itaty(m) .ne. 0 ) then
                  write(iot,'(
     &                     "  amin  &=&",1pe13.4,/
     &                     "  amax  &=&",1pe13.4)')
     &                     ab(ia), ab(ia+1)
               end if

                  write(iot,'("e:")')

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
            end do
            end do
            end do

*-----------------------------------------------------------------------
*        yz axis (matrix)
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 8 ) then

*-----------------------------------------------------------------------

               inum = 0

            do ix = 1, nx
            do ip = 1, npg
            do ie = 1, neg
            do it = 1, ntg
            do ia = 1, nag

               xval = ( xm(ix) + xm(ix+1) ) / 2.0d0

               inum = inum + 1

                  write(iot,'(/"#",78("-"))')

               if( inum .eq. 1 .or. ip.gt.np_mxang ) then ! frtati 2021/10/05
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

               if( rtfac(m) .lt. 0.d0 )
     &           write(iot,
     &           '( "#restart maximum for normalization:",1p1e14.7)')
     &           facmax(m)

               if ( ip.gt.np_mxang ) then
                 write(iot,'( " SKIPPAGE:")')
                 inum = inum - 1
               end if

                  write(iot,'( "#   no. =",i3,3x,
     &            "ie  =",i3,3x,
     &            "ix  =",i3,3x,
     &            "part. = ",a8,/
     &            "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )"/
     &            "#   x = ",1p1e13.4)')
     &                        inum, ie, ix, chq(ip),
     &                        eb(ie), eb(ie+1),
     &                        xval
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  it =",i3/
     &            "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                  it, tb(it), tb(it+1)
               end if

               if( itaty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  ia =",i3,/
     &            "# ",a3," = (",
     &                           1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     ia, aname, ab(ia), ab(ia+1)
               end if

*-----------------------------------------------------------------------

               if( ittty(m) .ne. 0 .and. itaty(m) .ne. 0 ) then
                  write(iot,'(/a1,"no. =",i3,
     &                     ",   e =",i3,",   x =",i3,
     &                     ",   t =",i3,"   ang =",i3,a1)')
     &                     cha, inum, ie, ix, it, ia, cha
               else if( ittty(m) .eq. 0 .and. itaty(m) .ne. 0 ) then
                  write(iot,'(/a1,"no. =",i3,
     &                     ",   e =",i3,",   x =",i3,
     &                     ",  ang =",i3,a1)')
     &                     cha, inum, ie, ix, ia, cha
               else if( ittty(m) .ne. 0 .and. itaty(m) .eq. 0 ) then
                  write(iot,'(/a1,"no. =",i3,
     &                     ",   e =",i3,",   x =",i3,
     &                     ",   t =",i3,a1)')
     &                     cha, inum, ie, ix, it, cha
               else
                  write(iot,'(/a1,"no. =",i3,
     &                     ",    e =",i3,
     &                     ",    x =",i3,a1)')
     &                     cha, inum, ie, ix, cha
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

                  val(i) = tr(ip,ie,it,ia,i,1)

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

         write(iot,'("y: Number ",a32)') hsunit(itunt(m))

      else

         write(iot,'("y: ",200a1)') (itazt(m)(i:i),i=1,itazl(m))

      end if

      if( izlog .eq. 0 ) write(iot,'( "p: ylin")')

      write(iot,'( "#",78("-"))')

      end if

*-----------------------------------------------------------------------

            if( inocm .eq. 1 ) then

                if( ite2l(m).eq.0 ) then  ! ccse 2022/09/30
                  if( iMeVperu.eq.0 ) then
                     write(iot,'(/"wt: s[c5]",/
     &                     a1,"vspace{-3}"/
     &                     "  emin  &=&",1pe13.4," [MeV]"/
     &                     "  emax  &=&",1pe13.4," [MeV]"/
     &                     "  x     &=&",1pe13.4," [cm]"/
     &                     "  part. &=&   ",a8)')
     &                     yen, eb(ie), eb(ie+1),
     &                     xval, chq(ip)
                  else
                     write(iot,'(/"wt: s[c5]",/
     &                     a1,"vspace{-3}"/
     &                     "  emin  &=&",1pe13.4," [MeV/n]"/
     &                     "  emax  &=&",1pe13.4," [MeV/n]"/
     &                     "  x     &=&",1pe13.4," [cm]"/
     &                     "  part. &=&   ",a8)')
     &                     yen, eb(ie), eb(ie+1),
     &                     xval, chq(ip)
                  end if
                else if( ite2l(m).eq.1 ) then
                     write(iot,'(/"wt: s[c5]",/
     &                     a1,"vspace{-3}"/
     &                     "  emin  &=&",1pe13.4," [keV/um]"/
     &                     "  emax  &=&",1pe13.4," [keV/um]"/
     &                     "  x     &=&",1pe13.4," [cm]"/
     &                     "  part. &=&   ",a8)')
     &                     yen, eb(ie), eb(ie+1),
     &                     xval, chq(ip)
                end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &                     "  tmin  &=&",1pe13.4," [nsec]"/
     &                     "  tmax  &=&",1pe13.4," [nsec]")')
     &                     tb(it), tb(it+1)
               end if
               if( itaty(m) .ne. 0 ) then
                  write(iot,'(
     &                     "  amin  &=&",1pe13.4,/
     &                     "  amax  &=&",1pe13.4)')
     &                     ab(ia), ab(ia+1)
               end if

                  write(iot,'("e:")')

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
            end do
            end do
            end do

*-----------------------------------------------------------------------
*        xz axis (matrix)
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 9 ) then

*-----------------------------------------------------------------------

               inum = 0

            do iy = 1, ny
            do ip = 1, npg
            do ie = 1, neg
            do it = 1, ntg
            do ia = 1, nag

               yval = ( ym(iy) + ym(iy+1) ) / 2.0d0

               inum = inum + 1

                  write(iot,'(/"#",78("-"))')

               if( inum .eq. 1 .or. ip.gt.np_mxang ) then ! frtati 2021/10/05
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

               if( rtfac(m) .lt. 0.d0 )
     &           write(iot,
     &           '( "#restart maximum for normalization:",1p1e14.7)')
     &           facmax(m)

               if ( ip.gt.np_mxang ) then
                 write(iot,'( " SKIPPAGE:")')
                 inum = inum - 1
               end if

                  write(iot,'("#   no. =",i3,3x,
     &            "ie  =",i3,3x,
     &            "iy  =",i3,3x,
     &            "part. = ",a8,/
     &            "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )"/
     &            "#   y = ",1p1e13.4)')
     &                        inum, ie, iy, chq(ip),
     &                        eb(ie), eb(ie+1),
     &                        yval
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  it =",i3/
     &            "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                  it, tb(it), tb(it+1)
               end if

               if( itaty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  ia =",i3,/
     &            "# ",a3," = (",
     &                           1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     ia, aname, ab(ia), ab(ia+1)
               end if

*-----------------------------------------------------------------------

               if( ittty(m) .ne. 0 .and. itaty(m) .ne. 0 ) then
                  write(iot,'(/a1,"no. =",i3,
     &                     ",   e =",i3,",   y =",i3,
     &                     ",   t =",i3,"   ang =",i3,a1)')
     &                     cha, inum, ie, iy, it, ia, cha
               else if( ittty(m) .eq. 0 .and. itaty(m) .ne. 0 ) then
                  write(iot,'(/a1,"no. =",i3,
     &                     ",   e =",i3,",   y =",i3,
     &                     ",  ang =",i3,a1)')
     &                     cha, inum, ie, iy, ia, cha
               else if( ittty(m) .ne. 0 .and. itaty(m) .eq. 0 ) then
                  write(iot,'(/a1,"no. =",i3,
     &                     ",   e =",i3,",   y =",i3,
     &                     ",   t =",i3,a1)')
     &                     cha, inum, ie, iy, it, cha
               else
                  write(iot,'(/a1,"no. =",i3,
     &                     ",    e =",i3,
     &                     ",    y =",i3,a1)')
     &                     cha, inum, ie, iy, cha
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

                  val(i) = tr(ip,ie,it,ia,i,1)

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

         write(iot,'("y: Number ",a32)') hsunit(itunt(m))

      else

         write(iot,'("y: ",200a1)') (itazt(m)(i:i),i=1,itazl(m))

      end if

      if( izlog .eq. 0 ) write(iot,'( "p: ylin")')

      write(iot,'( "#",78("-"))')

      end if

*-----------------------------------------------------------------------

            if( inocm .eq. 1 ) then

                if( ite2l(m).eq.0 ) then  ! ccse 2022/09/30
                  if( iMeVperu.eq.0 ) then
                     write(iot,'(/"wt: s[c5]",/
     &                     a1,"vspace{-3}"/
     &                     "  emin  &=&",1pe13.4," [MeV]"/
     &                     "  emax  &=&",1pe13.4," [MeV]"/
     &                     "  y     &=&",1pe13.4," [cm]"/
     &                     "  part. &=&   ",a8/
     &                     "e:")')
     &                     yen, eb(ie), eb(ie+1),
     &                     yval, chq(ip)
                  else
                     write(iot,'(/"wt: s[c5]",/
     &                     a1,"vspace{-3}"/
     &                     "  emin  &=&",1pe13.4," [MeV/n]"/
     &                     "  emax  &=&",1pe13.4," [MeV/n]"/
     &                     "  y     &=&",1pe13.4," [cm]"/
     &                     "  part. &=&   ",a8/
     &                     "e:")')
     &                     yen, eb(ie), eb(ie+1),
     &                     yval, chq(ip)
                  end if
                else if( ite2l(m).eq.1 ) then
                     write(iot,'(/"wt: s[c5]",/
     &                     a1,"vspace{-3}"/
     &                     "  emin  &=&",1pe13.4," [keV/um]"/
     &                     "  emax  &=&",1pe13.4," [keV/um]"/
     &                     "  y     &=&",1pe13.4," [cm]"/
     &                     "  part. &=&   ",a8/
     &                     "e:")')
     &                     yen, eb(ie), eb(ie+1),
     &                     yval, chq(ip)
                end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &                     "  tmin  &=&",1pe13.4," [nsec]"/
     &                     "  tmax  &=&",1pe13.4," [nsec]")')
     &                     tb(it), tb(it+1)
               end if
               if( itaty(m) .ne. 0 ) then
                  write(iot,'(
     &                     "  amin  &=&",1pe13.4,/
     &                     "  amax  &=&",1pe13.4)')
     &                     ab(ia), ab(ia+1)
               end if

                  write(iot,'("e:")')

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
            end do
            end do
            end do

*-----------------------------------------------------------------------

         end if

            call prestart(m,iot) !OBINATA(2012.9.11)

            close(iot)

         if( iteps(m) .ne. 0 ) then

            open(iot, file = fname, status = 'unknown' )
            call a_angel(idasa,fname)

         end if
cFURUTA20190208 OpenFOAM output
         if ( itfoam(m) .ne. 0) then
          ifilecount=0
          do ip = 1, np
          do ie = 1, ne
          do it = 1, nt
          do ia = 1, na
           ifilecount=ifilecount+1
           call openfoam_create_filename(
     &          itfoam(m),fname, foamfIType, foamfIndex, numIndex,
     &          ifilecount, outFilename)
           open(iot, file = outFilename, status='unknown')
           if( itfoam(m).eq.1)then !FURUTA20191028
            write(iot,'(a)')'('
            do ir=1,nr
             write(iot,'(1pe13.4)') tr(ip,ie,it,ia,ir,1)
            enddo
            write(iot,'(a)')')'
CFURUTA20191028 CSV output
           elseif( itfoam(m) .eq. 2)then
            if( itayl(m) .eq. 0 ) then
             write(sbuf,'( "Number",a32)') hsunit(itunt(m))
            else
             write(sbuf,'(200a1)')
     &            (itayt(m)(i:i),i=1,itayl(m))
            end if
            write(buf,'( "# tetra,xCM[cm],yCM[cm],zCM[cm],volume,
     &           ",a200,",r.err")')sbuf
            call remove_spaces(buf)
            write(iot,'(a)')trim(buf)
            do ir=1,nr
             write(buf,'(i8,5(",",1pe13.4),",",0pf8.4)')
     &            lr(ir), xcm(1:3,ir),vl(ir),
     &            (tr(ip,ie,it,ia,ir,k),k=1,2)
             call remove_spaces(buf)
             write(iot,'(a)')trim(buf)
            enddo
           endif
           close(iot)
          end do
          end do
          end do
          end do
          deallocate( foamfIType,foamfIndex )
         endif

  900 continue

*-----------------------------------------------------------------------
      deallocate( ixyz )
      deallocate( lr,vl,val ) !FURUTA20190204
      if(itfoam(m).eq.2) deallocate(xcm) !FURUTA20191028

      include 'samepage_include/samepage999.inc'

      return
      end


************************************************************************
*                                                                      *
      subroutine tpdctrz(ncol,m,np,nr,nz,ne,nt,na,nm,mt,nl,lt,
     &                   rm,zm,eb,tb,ab,tr,trEVENT,
     &                   itrmax,itrmin)
*                                                                      *
*       product tally in r-z scoring mesh                              *
*       last modified by K.Niita on 2011/05/17                         *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*      ncol  ..... reaction type                                       *
*              4 : source                                              *
*             13 : nuclear reaction (n,x)                              *
*             14 : nuclear reaction (n,n'x)                            *
*                                                                      *
*        jcoll : =  0, nothing happen                                  *
*                =  1, Hydrogen collisions                             *
*                =  2, Particle Decays                                 *
*                =  3, Elastic collisions                              *
*                =  4, Nuclear collisions                              *
*                =  5, Heavy Ion collisions                            *
*                =  6, Neutron reactions by data                       *
*                =  7, Photon reactions by data                        *
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
      use MMBANKMOD !FURUTA
      use partmod, only: itmxpt, itpan, itpat, jtpat ! frtati 2021/10/05
*-----------------------------------------------------------------------

      implicit double precision( a-h, o-z )

*-----------------------------------------------------------------------

      include 'param00.inc'
      include 'param.inc'

      parameter ( pi  = 3.141592653589793d0 )

*-----------------------------------------------------------------------

      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
!$OMP THREADPRIVATE(/icomon/)
      common /jcomon/ nabov,nobch,nocas,nomax
!$OMP THREADPRIVATE(/jcomon/)
      common /wtsave/ oldwt
!$OMP THREADPRIVATE(/wtsave/)
      common /mathzn/ mathz, mathn, jcoll, kcoll
!$OMP THREADPRIVATE(/mathzn/)
      common /kmat1d/ idmn(0:kvlmax), idnm(kvmmax)
      common /pnint/  ipnint

*-----------------------------------------------------------------------

      common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                itout(itlmax), ittwo(itlmax)

      common /tall08/ rtrx0(itlmax), rtry0(itlmax)
      common /tall12/ itman(itlmax), itmat(itlmax), itmct(itlmax),
     &                itnun(itlmax), itnuc(itlmax),
     &                itndz(itlmax), itndn(itlmax),
     &                itnkz(itlmax), itnkn(itlmax)
      common /tall22/ itmcn(itlmax), itmtn(itlmax), itmtt(itlmax)
      common /tall37/ itcnt(9,itlmax)
      common /tall46/ itmtr(itlmax,4), rtmtr(itlmax,13)
      common /tall48/ itmdp(itlmax,0:30), itmdf(itlmax)

      common /tall23/ itaty(itlmax), itanm(itlmax), itarg(itlmax),
     &                rtami(itlmax), rtama(itlmax), rtadl(itlmax)

      common /tall82/ itcnth(9,itlmax)
      common /tall74/ iprim(itlmax)

      common /tall50/ itlmt(itlmax), ite2l(itlmax)

*-----------------------------------------------------------------------

      common /clustt/ nclsts, iclusts(nnn)
!$OMP THREADPRIVATE(/clustt/)
      common /clustw/ jclusts(0:8,nnn),  qclusts(0:12,nnn)
!$OMP THREADPRIVATE(/clustw/)
      common /cntcls/ jcount(3,nnn)
!$OMP THREADPRIVATE(/cntcls/)

*-----------------------------------------------------------------------

      dimension   rm(nr+1)
      dimension   zm(nz+1)
      dimension   mt(nm)
      dimension   lt(nl)
      dimension   eb(ne+1)
      dimension   tb(nt+1)
      dimension   ab(na+1)
      dimension   tr(np,ne,nt,na,nr,nz,2)
      dimension   trEVENT(np,ne,nt,na,nr,nz)      !OBINATA(2012.9.11): as Ct
      real(8),allocatable,save:: tr0(:,:,:,:,:,:) !OBINATA(2012.9.11): as C
      dimension   itrmax(6),itrmin(6)

      dimension ips(itmxpt) ! frtati 2021/10/05 6 -> itmxpt
      dimension ncntt(3)
      dimension dmpd(30)

*-----------------------------------------------------------------------

      common /stat / istdev, irestart, ireschk
      common /cparm/ maxbch,maxcas

*-----------------------------------------------------------------------
      integer idmpomp !FURUTA20150427
      common /idmpomp0/idmpomp !FURUTA20150427
*-----------------------------------------------------------------------
      common /tstara/ atmrc(10,8)      ! Ogawa 2023/7/18, (10,7) -> (10,8) for plasmon
!$OMP THREADPRIVATE(/tstara/)
      common /tpdcta/ nflumu
!$OMP THREADPRIVATE(/tpdcta/)
      common /trskip/ ntrskip
!$OMP THREADPRIVATE(/trskip/)

      integer           iMeVperu
      common /cMeVperu/ iMeVperu

      common /cntmx/ ncntmx(3), nctmxsr(3)
!$OMP THREADPRIVATE(/cntmx/)

! sumover
      dimension ldo(2,6)

*-----------------------------------------------------------------------

      if (istdev .eq. 2) then

        call readitrminmax6(itrmin,itrmax,(/ np,ne,nt,na,nr,nz /),
     &                      mnp,mne,mnt,mna,mnr,mnz,
     &                      mxp,mxe,mxt,mxa,mxr,mxz)

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
* OBINATA(2012.9.11): change tr(,,,,3) to trEVENT(,,,,)
*-----------------------------------------------------------------------

         if (( ncol .eq. 0 .or. ncol .eq. 4 )
     &                              .and. istdev .eq. 2) then
           if ((nocas.gt.1.or.ncol.eq.0) .and. ihistcount.ne.1 ) then

             do iz = mnz,mxz
             do ir = mnr,mxr
             do ia = mna,mxa
             do it = mnt,mxt
             do ie = mne,mxe
               tr(:,ie,it,ia,ir,iz,1) = tr(:,ie,it,ia,ir,iz,1)
     &                                + trEVENT(:,ie,it,ia,ir,iz)
               tr(:,ie,it,ia,ir,iz,2) = tr(:,ie,it,ia,ir,iz,2)
     &                                + trEVENT(:,ie,it,ia,ir,iz) ** 2
             enddo
             enddo
             enddo
             enddo
             enddo

! sumover
             ldo(1,1) = 1
             ldo(2,1) = np
             ldo(1,2) = mne
             ldo(2,2) = mxe
             ldo(1,3) = mnt
             ldo(2,3) = mxt
             ldo(1,4) = mna
             ldo(2,4) = mxa
             ldo(1,5) = mnr
             ldo(2,5) = mxr
             ldo(1,6) = mnz
             ldo(2,6) = mxz
             call tpdctrz_sumover(m,1,
     &                   np,  ne, nt, na, nr, nz, trEVENT,ldo)


           end if

           do iz = mnz,mxz
           do ir = mnr,mxr
           do ia = mna,mxa
           do it = mnt,mxt
           do ie = mne,mxe
             trEVENT(:,ie,it,ia,ir,iz) = 0
           enddo
           enddo
           enddo
           enddo
           enddo

           call resetitrminmax(itrmin,itrmax,6,(/np,ne,nt,na,nr,nz/))

         end if

*-----------------------------------------------------------------------
*        end of batch ( in case of istdev = 1 )
*
* OBINATA(2012.9.11): modificate for thread parallel
*-----------------------------------------------------------------------

         if ( ncol .eq. 0 .and. istdev .eq. 1) then
!$OMP MASTER
             allocate( tr0(np,ne,nt,na,nr,nz) )
             tr0(:,:,:,:,:,:) = 0.d0
!$OMP END MASTER
!$OMP BARRIER
!$OMP CRITICAL (tpdctrz_crit_ist1)
             tr0(:,:,:,:,:,:) = tr0(:,:,:,:,:,:) + trEVENT(:,:,:,:,:,:)
!$OMP END CRITICAL (tpdctrz_crit_ist1)
!$OMP BARRIER
!$OMP MASTER
             tr(:,:,:,:,:,:,1) = tr(:,:,:,:,:,:,1)
     &                         + tr0(:,:,:,:,:,:) / maxcas
             tr(:,:,:,:,:,:,2) = tr(:,:,:,:,:,:,2)
     &                         + ( tr0(:,:,:,:,:,:) / maxcas ) ** 2

! sumover
             ldo(1,1) = 1
             ldo(2,1) = np
             ldo(1,2) = 1
             ldo(2,2) = ne
             ldo(1,3) = 1
             ldo(2,3) = nt
             ldo(1,4) = 1
             ldo(2,4) = na
             ldo(1,5) = 1
             ldo(2,5) = nr
             ldo(1,6) = 1
             ldo(2,6) = nz
             call tpdctrz_sumover(m,maxcas,
     &                   np,  ne, nt, na, nr, nz, tr0, ldo)

             deallocate( tr0 )
!$OMP END MASTER

           trEVENT(:,:,:,:,:,:) = 0

         end if

*-----------------------------------------------------------------------
*        out put : check of ncol, jcoll and kcoll
*-----------------------------------------------------------------------


         if( itout(m) .eq. 1 .and. ncol .ne. 4 ) return
         if( itout(m) .gt. 1 .and.
     &       ncol .ne. 13 .and. ncol .ne. 14 ) return


*-----------------------------------------------------------------------

         if( itout(m) .eq. 2 .and. ! output = nuclear (See subroutine tprodct)
     &       jcoll .ne. 1  .and. jcoll .ne. 3  .and.
     &       jcoll .ne. 4  .and. jcoll .ne. 5  .and.
     &       jcoll .ne. 6  .and. jcoll .ne. 9  .and.      ! ccse 2023/03/29
     &       jcoll .ne. 10 .and. jcoll .ne. 15 .and.
     &       jcoll .ne. 16 .and. jcoll .ne. 17 .and.      ! S.Abe 2017/02/08
     &       jcoll .ne. 19 .and. jcoll .ne. 20 ) return   ! y.sakaki 2024/01

         if( itout(m) .eq. 3 .and. jcoll .ne. 2 ) return
         if( itout(m) .eq. 3 .and.
     &       jcoll .eq. 2 .and. ityp .eq. 13 ) return     ! S.Abe 2017/02/08

         if( itout(m) .eq. 4 .and.
     &       kcoll .ne. 1 .and. kcoll .ne. 5 ) return

         if( itout(m) .eq. 5 .and.
     &       jcoll .ne. 3 .and. kcoll .ne. 3 ) return

         if( itout(m) .eq. 6 .and. ! output = nonela (See subroutine tprodct)
     &       jcoll .ne. 1  .and.
     &       jcoll .ne. 4  .and. jcoll .ne. 5  .and.
     &       jcoll .ne. 6  .and. jcoll .ne. 9  .and.      ! ccse 2023/03/29
     &       jcoll .ne. 10 .and. jcoll .ne. 15 .and.
     &       jcoll .ne. 16 .and. jcoll .ne. 17 .and.
     &       jcoll .ne. 19 .and. jcoll .ne. 20 ) return   ! y.sakaki 2024/01
         if( itout(m) .eq. 6 .and. kcoll .eq. 3 ) return

         sumatmrc = 0.d0
         do i = 1, 5
         do j = 1, 5    ! except mscat
            sumatmrc = sumatmrc + atmrc(i,j)
         enddo
         enddo
         if( itout(m) .eq. 7 .and.
     &       jcoll    .ne. 18 .and.
     &       ( jcoll .ne. 11 .and. sumatmrc .le. 0.d0 )
     &     ) return

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
*        for source
*-----------------------------------------------------------------------

         if( ncol .eq. 4 ) then

               iccol = 0
               npart = nomax

               ncntt(1) = ncnt(ibknct+1,no,ipomp+1)
               ncntt(2) = ncnt(ibknct+2,no,ipomp+1)
               ncntt(3) = ncnt(ibknct+3,no,ipomp+1)

*-----------------------------------------------------------------------
*        for reaction, check of reaction and mother
*-----------------------------------------------------------------------

         else

            if( nclsts .le. 0 ) return

               iccol = 1
               npart = nclsts

               ata   = dble( mathz + mathn )
               atz   = dble( mathz )
               mmas  = nint( ata )
               mchg  = nint( atz )

               ncntt(1) = jcount(1,1)
               ncntt(2) = jcount(2,1)
               ncntt(3) = jcount(3,1)

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

                     if( itmct(m) .gt. 0 ) goto 140
                     if( itmct(m) .lt. 0 ) return

                  end if

               end do

                  if( itmct(m) .gt. 0 ) return

            end if

  140       continue

*-----------------------------------------------------------------------
*           check elastic collision
*-----------------------------------------------------------------------

            if( ( jcoll .eq. 1 .or.
     &            jcoll .eq. 4 .or.
     &            jcoll .eq. 5 .or.
     &            jcoll .eq. 10 ) .and.
     &            npart .eq. 2 ) then

               if( jcoll .eq. 1 ) then

                  if( ( iclusts(1) .eq. 1 .and.
     &                  jclusts(7,2) .eq. ktyp .and.
     &                  qclusts(6,2) .eq. 0.d0 ) .or.
     &                ( iclusts(2) .eq. 1 .and.
     &                  jclusts(7,1) .eq. ktyp .and.
     &                  qclusts(6,1) .eq. 0.d0 ) ) then

                     if( itout(m) .ne. 2 .and. itout(m) .ne. 5 ) return

                  end if

               else

                  if( ( iclusts(1) .eq. 0 .and.
     &                  jclusts(1,1) .eq. mchg .and.
     &                  jclusts(2,1) .eq. mmas - mchg .and.
     &                  jclusts(7,2) .eq. ktyp .and.
     &                  qclusts(6,2) .eq. 0.d0 .and.
     &                  qclusts(6,1) .eq. 0.d0 ) .or.
     &                ( iclusts(2) .eq. 0 .and.
     &                  jclusts(1,2) .eq. mchg .and.
     &                  jclusts(2,2) .eq. mmas - mchg .and.
     &                  jclusts(7,1) .eq. ktyp .and.
     &                  qclusts(6,1) .eq. 0.d0 .and.
     &                  qclusts(6,2) .eq. 0.d0 ) ) then

                     if( itout(m) .ne. 2 .and. itout(m) .ne. 5 ) return

                  end if

               end if

            end if

*-----------------------------------------------------------------------

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
*        material for LET
*-----------------------------------------------------------------------

               if( itlmt(m) .gt. 0 ) then

                  lmat = idnm( itlmt(m) )

               else if( itlmt(m) .eq. 0 ) then

                  lmat = mat

               else

                  lmat = -idnm( -itlmt(m) )

               end if

*-----------------------------------------------------------------------
*        do loop for particles
*-----------------------------------------------------------------------

         do 100 j = 1, npart

*-----------------------------------------------------------------------
            if( ityp .eq. 7 .and. (jcoll.eq.2 .or. jcoll.eq.16) ) then
               if( itout(m) .eq. 7 .and. j .gt. nflumu ) goto 100
               if( itout(m) .ne. 7 .and. j .le. nflumu ) goto 100
            endif
            if( j .eq. ntrskip ) goto 100

            if( iccol .eq. 1 ) then

                  ipart = jclusts(3,j)
                  kpart = jclusts(7,j)
                  lpart = jclusts(8,j)
                  jpart = ichgf(ipart,kpart)
                  rpart = qclusts(5,j) * 1000.d0   ! ccse 2022/08/31

                  epart = qclusts(7,j)
                  tlw   = qclusts(8,j)
                  tpart = abs(qclusts(9,j))
                  xpart = qclusts(10,j)
                  ypart = qclusts(11,j)
                  zpart = qclusts(12,j)

                  upart = qclusts(1,j)
                  vpart = qclusts(2,j)
                  wpart = qclusts(3,j)

                  sxpat = spx(ibkspx+no,ipomp+1)
                  sypat = spy(ibkspy+no,ipomp+1)
                  szpat = spz(ibkspz+no,ipomp+1)
                  nzpat = jclusts(5,j)
                  nname = name(ibknam+no,ipomp+1)

            else

                  nj = no + j - 1

                  ipart = nty(ibknty+nj,ipomp+1)
                  kpart = nkf(ibknkf+nj,ipomp+1)
                  jpart = ichgf(ityp,ktyp)
                  rpart = rtyp   ! ccse 2022/08/31

                  epart = e(ibke+nj,ipomp+1)
                  tlw   = wt(ibkwt+nj,ipomp+1)
                  tpart = abs(t(ibkt+nj,ipomp+1))

                  xpart = x(ibkx+nj,ipomp+1)
                  ypart = y(ibky+nj,ipomp+1)
                  zpart = z(ibkz+nj,ipomp+1)

                  upart = u(ibku+nj,ipomp+1)
                  vpart = v(ibkv+nj,ipomp+1)
                  wpart = w(ibkw+nj,ipomp+1)

                  sxpat = spx(ibkspx+nj,ipomp+1)
                  sypat = spy(ibkspy+nj,ipomp+1)
                  szpat = spz(ibkspz+nj,ipomp+1)
                  nzpat = ctyp
                  nname = name(ibknam+nj,ipomp+1)

            end if

             if( kpart .eq. ktyp .and. iprim(m) .eq. 0 .and.   !exclude primary particle in 
     &      ( (jcoll .eq. 18 .or. jcoll .eq. 14 .or. jcoll .eq. 13 .or. ! Track-structure and EGS
     &         jcoll .eq. 3) .and. j .eq. 1 .or.! elastic scattering
     &         jcoll .eq. 10 .and. j .eq. 1 .and. npart .eq. 2 .or.! elastic scattering
     &         jcoll .ge. 6 .AND. jcoll .le. 9   .and. j .eq. npart !exclude primary sampled from X-section data
     &      ) ) goto 100 ! rejection of remaining projectile
             
            
            if( epart .le. 0.d0 ) goto 100   ! S.Abe 2017/02/08

*-----------------------------------------------------------------------
*           transform positions
*-----------------------------------------------------------------------

               call trnsxx(xpart,ypart,zpart,
     &                     xcc,ycc,zcc,itmtr(m,4))

*-----------------------------------------------------------------------
*           check z mesh and r mesh region
*-----------------------------------------------------------------------

               if(  zcc .lt. zm(1) ) goto 100
               if(  zcc .ge. zm(nz+1) ) goto 100

                  x0 = rtrx0(m)
                  y0 = rtry0(m)

                  dis0 = sqrt( (  xcc - x0 )**2
     &                       + (  ycc - y0 )**2 )

               if( dis0 .lt. rm(1) ) goto 100
               if( dis0 .ge. rm(nr + 1) ) goto 100

*-----------------------------------------------------------------------
*           z-position
*-----------------------------------------------------------------------

               do i = 1, nz

                  if( zcc .ge. zm(i) .and.
     &                zcc .lt. zm(i+1) ) goto 38

               end do

   38          iz1 = i

*-----------------------------------------------------------------------
*           r-position
*-----------------------------------------------------------------------

               do i = 1, nr

                  if( dis0 .ge. rm(i) .and.
     &                dis0 .lt. rm(i+1) ) goto 47

               end do

   47          ir1 = i

*-----------------------------------------------------------------------
*           MeV -> MeV/n energy conversion ! 2018/3/7 Ogawa
*-----------------------------------------------------------------------

               if(iMeVperu.eq.1 .and. ipart.ge.15 .and. ipart.le.19)then
                    ebm = dble(kpart - kpart / 1000000 * 1000000)
               else
                    ebm = 1.d0
               end if

*-----------------------------------------------------------------------
*           check of particles
*-----------------------------------------------------------------------

            call pcheck(m,np,ipart,kpart+1000000000*lpart,jpart,ipn,ips)

               if( ipn .eq. 0 ) goto 100

*-----------------------------------------------------------------------
*           check of energy
*-----------------------------------------------------------------------

            if (ite2l(m) .eq. 0) then   ! ccse 2022/08/31

               if( epart .lt. eb(1)*ebm ) goto 100
               if( epart .ge. eb(ne+1)*ebm ) goto 100

*-----------------------------------------------------------------------
*           energy
*-----------------------------------------------------------------------

            do i = 1, ne

               if( epart .ge. eb(i)*ebm .and.
     &             epart .lt. eb(i+1)*ebm ) goto 30

            end do

            else if (ite2l(m) .eq. 1) then   ! convert to energy to LET

               call dedxas(epart,dedx,lmat,ipart,kpart,jpart,rpart)
               dedx = dedx /10.0d0

               ! check of energy
               if( dedx .lt. eb(1) ) goto 100
               if( dedx .ge. eb(ne+1) ) goto 100

            ! energy
            do i = 1, ne

               if( dedx .ge. eb(i) .and.
     &             dedx .lt. eb(i+1) ) goto 30

            end do

            end if

   30          ie1 = i

*-----------------------------------------------------------------------
*           check of time
*-----------------------------------------------------------------------

               if( tpart .lt. tb(1) ) goto 100
               if( tpart .ge. tb(nt+1) ) goto 100

*-----------------------------------------------------------------------
*           time
*-----------------------------------------------------------------------

            do i = 1, nt

               if( tpart .ge. tb(i) .and.
     &             tpart .lt. tb(i+1) ) goto 40

            end do

   40          it1 = i

*-----------------------------------------------------------------------
*        check angle
*-----------------------------------------------------------------------

         if( itaty(m) .ne. 0 ) then

               call trnsuu(upart,vpart,wpart,
     &                     ucc,vcc,wcc,itmtr(m,4))

                  cst = 1.d0
                  pab = sqrt( ucc**2 + vcc**2 + wcc**2 )
                  if( pab .gt. 0.0d0 ) cst = wcc / pab

               if( itaty(m) .gt. 0 ) then

                  if( cst .lt. ab(1) ) goto 100
                  if( cst .gt. ab(na+1) ) goto 100

               else

                  if( cst .lt. cos( ab(na+1) / 180.d0 * pi ) ) goto 100
                  if( cst .gt. cos( ab(1) / 180.d0 * pi ) ) goto 100

               end if

            do i = 1, na

               if( itaty(m) .gt. 0 ) then

                  if( cst .ge. ab(i) .and.
     &                cst .le. ab(i+1) ) goto 50

               else

                  if( cst .ge. cos( ab(i+1) / 180.d0 * pi ) .and.
     &                cst .le. cos( ab(i) / 180.d0 * pi ) ) goto 50

               end if

            end do

   50          ia1 = i

         else

               ia1 = 1

         end if

*-----------------------------------------------------------------------
*        tally
*-----------------------------------------------------------------------
               do ip = 1, ipn

                  trEVENT(ips(ip),ie1,it1,ia1,ir1,iz1) =
     &            trEVENT(ips(ip),ie1,it1,ia1,ir1,iz1) + tlw

               end do

               if (istdev .eq. 2) then

                 call setitrminmax(itrmin,itrmax,2,6,
     &                             (/ie1,it1,ia1,ir1,iz1/))

               endif
*-----------------------------------------------------------------------
*        dump data on file
*-----------------------------------------------------------------------

         if( itmdp(m,0) .ne. 0 ) then
               dmpd(1)  = dble( kpart )
               dmpd(2)  = xpart
               dmpd(3)  = ypart
               dmpd(4)  = zpart
               dmpd(5)  = upart
               dmpd(6)  = vpart
               dmpd(7)  = wpart
               dmpd(8)  = epart
               dmpd(9)  = tlw
               dmpd(10) = tpart
               dmpd(11) = ncntt(1)
               dmpd(12) = ncntt(2)
               dmpd(13) = ncntt(3)
               dmpd(14) = sxpat
               dmpd(15) = sypat
               dmpd(16) = szpat
               dmpd(17) = nname
               dmpd(18) = nocas
               dmpd(19) = nobch
               dmpd(20) = no
               dmpd(21) = nzpat

               if( ipart .ge. 15 )
     &         dmpd(8) = dmpd(8) / ibryf(ipart,kpart)

            if(idmpomp.ne.0)then       !FURUTA20150427
             call writeompfile(m,dmpd) !FURUTA20150427
            else                       !FURUTA20150427

             if( itmdp(m,0) .gt. 0 ) then

               write(itmdf(m))
     &         ( dmpd(itmdp(m,k)), k = 1, abs( itmdp(m,0) ) )

              else

               write(itmdf(m),'(30(1p1d24.15))')
     &         ( dmpd(itmdp(m,k)), k = 1, abs( itmdp(m,0) ) )

              end if

            endif !FURUTA20150427

         end if

*-----------------------------------------------------------------------

  100    continue

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine ppdctrz(m,np,nr,nz,ne,nt,na,rm,zm,eb,tb,ab,tr,
     &                   idasa)
*                                                                      *
*       output r-z scoring mesh product tally                          *
*       last modified by K.Niita on 2005/11/24                         *
*                                                                      *
************************************************************************
      use sangelmod, only: itsans
      use partmod ! frtati 2021/10/05

*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'

      parameter ( pi  = 3.141592653589793d0 )

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
      common /tall07/ itety(itlmax), itenm(itlmax), iterg(itlmax),
     &                rtemi(itlmax), rtema(itlmax), rtedl(itlmax)

      common /tall20/ ittty(itlmax), ittnm(itlmax), ittrg(itlmax),
     &                rttmi(itlmax), rttma(itlmax), rttdl(itlmax)


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
      common /tall48/ itmdp(itlmax,0:30), itmdf(itlmax)

      common /tall23/ itaty(itlmax), itanm(itlmax), itarg(itlmax),
     &                rtami(itlmax), rtama(itlmax), rtadl(itlmax)

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

      common /tall50/ itlmt(itlmax), ite2l(itlmax)
*-----------------------------------------------------------------------

      dimension   rm(nr+1)
      dimension   zm(nz+1)
      dimension   eb(ne+1)
      dimension   tb(nt+1)
      dimension   ab(na+1)
      dimension   ew(ne)
      dimension   tw(nt)
      dimension   aw(na)
      dimension   tr(np,ne,nt,na,nr,nz,2)

*-----------------------------------------------------------------------

      dimension tott(np,2) ! frtati 2021/10/05 6 -> np

      character hsunit(36)*32

      data hsunit( 1) / '[1/source]                      '/
      data hsunit( 2) / '[1/cm^3/source]                 '/
      data hsunit( 3) / '[1/MeV/source]                  '/
      data hsunit( 4) / '[1/cm^3/MeV/source]             '/
      data hsunit( 5) / '[1/Lethargy/source]             '/
      data hsunit( 6) / '[1/cm^3/Lethargy/source]        '/
      data hsunit(11) / '[1/nsec/source]                 '/
      data hsunit(12) / '[1/cm^3/nsec/source]            '/
      data hsunit(13) / '[1/MeV/nsec/source]             '/
      data hsunit(14) / '[1/cm^3/MeV/nsec/source]        '/
      data hsunit(15) / '[1/Lethargy/nsec/source]        '/
      data hsunit(16) / '[1/cm^3/Lethargy/nsec/source]   '/
      data hsunit(21) / '[1/sr/source]                   '/
      data hsunit(22) / '[1/cm^3/sr/source]              '/
      data hsunit(23) / '[1/MeV/sr/source]               '/
      data hsunit(24) / '[1/cm^3/MeV/sr/source]          '/
      data hsunit(25) / '[1/Lethargy/sr/source]          '/
      data hsunit(26) / '[1/cm^3/Lethargy/sr/source]     '/
      data hsunit(31) / '[1/nsec/sr/source]              '/
      data hsunit(32) / '[1/cm^3/nsec/sr/source]         '/
      data hsunit(33) / '[1/MeV/nsec/sr/source]          '/
      data hsunit(34) / '[1/cm^3/MeV/nsec/sr/source]     '/
      data hsunit(35) / '[1/Lethargy/nsec/sr/source]     '/
      data hsunit(36) / '[1/cm^3/Lethargy/nsec/sr/source]'/

      integer           iMeVperu
      common /cMeVperu/ iMeVperu

      character cha*1
      data cha /"'"/

      character dc2*4
      character chp(itmxpt)*11  ! kitamura22/03/31
      character chq(itmxpt)*9   ! kitamura22/03/31

      character aname*3

*-----------------------------------------------------------------------

      character rpa*1
      data rpa /'}'/
      character yen*1

! sumover
      real(8),allocatable :: vl_r(:),vl_z(:)

      include 'samepage_include/samepage000.inc'

*-----------------------------------------------------------------------
*        set mesh volume
*-----------------------------------------------------------------------

                  vl(ir,iz) = pi * ( rm(ir+1)**2 - rm(ir)**2 )
     &                           * ( zm(iz+1) - zm(iz) )

*-----------------------------------------------------------------------

      include 'samepage_include/samepage001.inc'

      allocate (vl_r(nz),vl_z(nr))

      if( itunt(m) .eq.  2 .or. itunt(m) .eq.  4 .or.
     &    itunt(m) .eq.  6 .or. itunt(m) .eq. 12 .or.
     &    itunt(m) .eq. 14 .or. itunt(m) .eq. 16 .or.
     &    itunt(m) .eq. 22 .or. itunt(m) .eq. 24 .or.
     &    itunt(m) .eq. 26 .or. itunt(m) .eq. 32 .or.
     &    itunt(m) .eq. 34 .or. itunt(m) .eq. 36 ) then
        vl_r(:) = 0.0d0
        vl_z(:) = 0.0d0
        do iz = i, nz
          do ir = 1,nr
            vl_r(iz) = vl_r(iz) + vl(ir,iz)
            vl_z(ir) = vl_z(ir) + vl(ir,iz)
          enddo
        enddo
      else
        vl_r(:) = 1.0d0
        vl_z(:) = 1.0d0
      endif

      yen  = char(92)
      igsh = 0

*-----------------------------------------------------------------------

      if ( iMeVperu.eq. 1 ) then
         hsunit( 3) = '[1/(MeV/n)/source]              '
         hsunit( 4) = '[1/cm^3/(MeV/n)/source]         '
         hsunit(13) = '[1/(MeV/n)/nsec/source]         '
         hsunit(14) = '[1/cm^3/(MeV/n)/nsec/source]    '
         hsunit(23) = '[1/(MeV/n)/sr/source]           '
         hsunit(24) = '[1/cm^3/(MeV/n)/sr/source]      '
         hsunit(33) = '[1/(MeV/n)/nsec/sr/source]      '
         hsunit(34) = '[1/cm^3/(MeV/n)/nsec/sr/source] '
      end if

      if ( ite2l(m) .eq. 1 ) then   ! convert to energy to LET
         hsunit( 3) = '[1/(keV/um)/source]             '
         hsunit( 4) = '[1/cm^3/(keV/um)/source]        '
         hsunit(13) = '[1/(keV/um)/nsec/source]        '
         hsunit(14) = '[1/cm^3/(keV/um)/nsec/source]   '
         hsunit(23) = '[1/(keV/um)/sr/source]          '
         hsunit(24) = '[1/cm^3/(keV/um)/sr/source]     '
         hsunit(33) = '[1/(keV/um)/nsec/sr/source]     '
         hsunit(34) = '[1/cm^3/(keV/um)/nsec/sr/source]'
      end if
*-----------------------------------------------------------------------
      np_mxang = np
      if ( itmxang(m).ne.0 .and. iabs(itmxang(m)).lt.np ) then
        np_mxang = iabs(itmxang(m))
      end if

*-----------------------------------------------------------------------

            if( itaty(m) .gt. 0 ) then
               aname = 'cos'
            else
               aname = 'the'
            end if

*-----------------------------------------------------------------------
*        particle name
*-----------------------------------------------------------------------

               call wtpname(m,np,chp,chq)

*-----------------------------------------------------------------------
*        c1 : nomalization for source
*-----------------------------------------------------------------------

            if ( rsouin .gt. 0d0 ) then ! S.H. (2016.7.26)

               c1 = 1.0d+0 / rsouin

            else

               c1 = 0d0

            end if

*-----------------------------------------------------------------------
*           itunt(m) = 3 13 23 33 4 14 24 34: /MeV
*-----------------------------------------------------------------------

            if( itunt(m) .eq.  3 .or. itunt(m) .eq.  4 .or.
     &          itunt(m) .eq. 13 .or. itunt(m) .eq. 14 .or.
     &          itunt(m) .eq. 23 .or. itunt(m) .eq. 24 .or.
     &          itunt(m) .eq. 33 .or. itunt(m) .eq. 34 ) then

               do i = 1, ne

                  ew(i) = eb(i+1) - eb(i)

               end do
               ew_sum = eb(ne+1) - eb(1)

            else if( itunt(m) .eq.  5 .or. itunt(m) .eq.  6 .or.
     &               itunt(m) .eq. 15 .or. itunt(m) .eq. 16 .or.
     &               itunt(m) .eq. 25 .or. itunt(m) .eq. 26 .or.
     &               itunt(m) .eq. 35 .or. itunt(m) .eq. 36 ) then

               do i = 1, ne

                  ew(i) = log( eb(i+1) / eb(i) )

               end do
               ew_sum = log( eb(ne+1) / eb(1) )

            else

               do i = 1, ne

                  ew(i) = 1.d+0

               end do
               ew_sum = 1.0d0

            end if

*-----------------------------------------------------------------------
*           itunt(m) = 11, 12, 13, 14, 31, 32, 33, 34 : /nsec
*-----------------------------------------------------------------------

            if( ( itunt(m) .ge. 11 .and. itunt(m) .le. 14 ) .or.
     &          ( itunt(m) .ge. 31 .and. itunt(m) .le. 34 ) ) then

               do i = 1, nt

                  tw(i) = tb(i+1) - tb(i)

               end do
               tw_sum = tb(nt+1) - tb(1)

            else

               do i = 1, nt

                  tw(i) = 1.d+0

               end do
               tw_sum = 1.0d0

            end if

*-----------------------------------------------------------------------
*           itunt(m) = 20+, 30+ : /SR
*-----------------------------------------------------------------------

            if( itunt(m) .lt. 20 ) then

               do i = 1, na

                  aw(i) = 1.d+0

               end do
               aw_sum = 1.0d0

            else
               aw_sum = 0.0d0
               do i = 1, na

                  if( itaty(m) .gt. 0 ) then

                     aw(i) = ( ab(i+1) - ab(i) ) * 2.0 * pi

                  else

                     aw(i) = ( cos( ab(i) / 180.d0 * pi )
     &                       - cos( ab(i+1) / 180.d0 * pi ) ) * 2.0 * pi

                  end if
                  aw_sum = aw_sum + aw(i)

               end do

            end if

*-----------------------------------------------------------------------
*           itunt(m) = 2, 4, 6, +10,20,30  : /cm^3
*-----------------------------------------------------------------------
*        relative error and  unit conversion
*-----------------------------------------------------------------------

                  cmax = 0.0
                  cmin = 1.e+33
                  isdz = 0
                  stdm = 0.0
                  itstd(m) = 1

            facmax(m) = 1.d0
            if( rtfac(m) .lt. 0.d0 ) then
               facmax(m) = 0.d0

               do 101 iz = 1, nz
               do 101 ir = 1, nr
               do 101 ie = 1, ne
               do 101 it = 1, nt
               do 101 ia = 1, na
               do 101 ip = 1, np

                  if( tr(ip,ie,it,ia,ir,iz,1) .gt. 0.d0 ) then

                     if( itunt(m) .eq.  2 .or. itunt(m) .eq.  4 .or.
     &                   itunt(m) .eq.  6 .or. itunt(m) .eq. 12 .or.
     &                   itunt(m) .eq. 14 .or. itunt(m) .eq. 16 .or.
     &                   itunt(m) .eq. 22 .or. itunt(m) .eq. 24 .or.
     &                   itunt(m) .eq. 26 .or. itunt(m) .eq. 32 .or.
     &                   itunt(m) .eq. 34 .or. itunt(m) .eq. 36 ) then

                        vm= vl(ir,iz)

                     else

                        vm = 1.d+0

                     end if

                     fmaxfc = tr(ip,ie,it,ia,ir,iz,1)
     &                                / vm / ew(ie) / tw(it) / aw(ia)

                     if( facmax(m) .lt. fmaxfc) facmax(m) = fmaxfc

                  end if

  101          continue

            end if

            do 100 iz = 1, nz
            do 100 ir = 1, nr
            do 100 ie = 1, ne
            do 100 it = 1, nt
            do 100 ia = 1, na
            do 100 ip = 1, np


                  if( itunt(m) .eq.  2 .or. itunt(m) .eq.  4 .or.
     &                itunt(m) .eq.  6 .or. itunt(m) .eq. 12 .or.
     &                itunt(m) .eq. 14 .or. itunt(m) .eq. 16 .or.
     &                itunt(m) .eq. 22 .or. itunt(m) .eq. 24 .or.
     &                itunt(m) .eq. 26 .or. itunt(m) .eq. 32 .or.
     &                itunt(m) .eq. 34 .or. itunt(m) .eq. 36 ) then

                     vm= vl(ir,iz)

                  else

                     vm = 1.d+0

                  end if

               if( tr(ip,ie,it,ia,ir,iz,1) .gt. 0.d0 ) then

                  cc = abs(rtfac(m)/facmax(m))
     &                          / vm / ew(ie) / tw(it) / aw(ia)

                  call calc_stdev(m,Xa,sigx,
     &                            tr(ip,ie,it,ia,ir,iz,1),
     &                            tr(ip,ie,it,ia,ir,iz,2),
     &                            cc)

                  tr(ip,ie,it,ia,ir,iz,1) = Xa
                  tr(ip,ie,it,ia,ir,iz,2) = sigx
                  if( sigx .gt. stdm ) stdm = sigx

                  if( tr(ip,ie,it,ia,ir,iz,1) .gt. cmax )
     &                                   cmax = tr(ip,ie,it,ia,ir,iz,1)

                  if( tr(ip,ie,it,ia,ir,iz,1) .lt. cmin )
     &                                   cmin = tr(ip,ie,it,ia,ir,iz,1)

               else

                  isdz = 1
                  tr(ip,ie,it,ia,ir,iz,2) = 0.0

               end if

! sumover
               fact_in = abs(rtfac(m)/facmax(m))
               call ppdctrz_sumover_stdev(0,m,ip,ie,it,ia,ir,iz,
     &                fact_in,ew(ie),tw(it),aw(ia),vm,
     &                ew_sum,tw_sum,aw_sum,vl_r(iz),vl_z(ir))


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

!OBINATA(2012.9.11): output *.err
        noe = 1
        if ( any( itaxs(m,iax) .eq. (/ 7, 8, 9, 10 /) )
     &      .and. ittwo(m) .ne. 4 .and. icntl .ne. 8 )
     &    noe = 2

        do ioe = 1, noe

         if( itmdp(m,0) .eq. 0 ) then

            if( itall .eq. 2 .and. nobch .lt. maxbch ) then

               write(fnume,'(i3.3)') nobch
!OBINATA(2012.9.11): output *.err
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
     &               nobch,maxbch,npe)
              else
                call mk_2dnerfn(ctfln(m,iax),fname,itfll(m,iax),fnume)
              end if

            else

!OBINATA(2012.9.11): output *.err
              if ( ioe .eq. 1 ) then
                fname = ctfln(m,iax)
              else
                call mk_2derrfn(ctfln(m,iax),fname,itfll(m,iax))
              end if

            end if

         else

            if( ( itall .eq. 2 .and. nobch .lt. maxbch )
     &           .or. itall .eq. 4 ) then

               call mk_2dnumfn(ctfln(m,iax),fname,itfll(m,iax),
     &              nobch,maxbch,npe)

            else

               fname = ctfln(m,iax)(1:itfll(m,iax))

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

               call tproech(iot,m,iax,1)

*-----------------------------------------------------------------------
*        energy axis  ( LET axis )
*-----------------------------------------------------------------------

      include 'samepage_include/samepage002_petarz.inc'
      include 'samepage_include/samepagechp_ptarz.inc'
      include 'samepage_include/samepageseti.inc'

         if( itaxs(m,iax) .eq. 1 .or.
     &       itaxs(m,iax) .eq. 15 ) then

               inum = 0

            do iri = 1, nr, nrstepi
            do izi = 1, nz, nzstepi
            do iti = 1, nt, ntstepi
            do iai = 1, na, nastepi
            do ipi = 1, np, npstepi
               ir = iri
               iz = izi
               it = iti
               ia = iai
               ip = ipi

               inum = inum + 1

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

               if(iloopmode .ge. 1) then
                  write(iot,'("#   ",a)') chp(ipi)(3:10)
               endif

                  write(iot,'("#   no. =",i3,3x,
     &            "ir  =",i3,3x,
     &            "iz  =",i3/
     &            "#   r = (",1p1e13.4,"  -",1p1e13.4,"  )"/
     &            "#   z = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        inum, ir, iz,
     &                        rm(ir), rm(ir+nrstepi),
     &                        zm(iz), zm(iz+nzstepi)
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  it =",i3/
     &            "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                  it, tb(it), tb(it+ntstepi)
               end if

               if( itaty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  ia =",i3,/
     &            "# ",a3," = (",
     &                           1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     ia, aname, ab(ia), ab(ia+nastepi)
               end if

               if( itaxl(m) .eq. 0 ) then

                if ( itaxs(m,iax) .eq. 1 ) then   ! ccse 2022/09/30

                if ( iMeVperu.eq. 1 ) then
                  write(iot,'(/"x: Energy [MeV/n]")')
                else
                  write(iot,'(/"x: Energy [MeV]")')
                end if

                else if ( itaxs(m,iax) .eq. 15 ) then
                 write(iot,'(/"x: LET [keV/um]")')
                end if

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: Number ",a32)') hsunit(itunt(m))

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

               if( itety(m) .eq.  3 .or. itety(m) .eq.  5 .or.
     &             itunt(m) .eq.  5 .or. itunt(m) .eq.  6 .or.
     &             itunt(m) .eq. 15 .or. itunt(m) .eq. 16 .or.
     &             itunt(m) .eq. 25 .or. itunt(m) .eq. 26 .or.
     &             itunt(m) .eq. 35 .or. itunt(m) .eq. 36 ) then

                  write(iot,'( "p: xlog ylog afac(0.8) form(0.9)")')

               else

                  write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')

               end if

               if( itanl(m) .gt. 0 ) then

                  write(iot,'( "p: ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               end if

               if( itsans(m) .gt. 0 ) then

                  call write_sangel(iot,m,0)

               end if

       include 'samepage_include/petarz_e.inc'


             if(iloopmode .eq. 0) then
               if( ittty(m) .ne. 0 .and. itaty(m) .ne. 0 ) then
                  write(iot,'(/a1,"no. =",i3,
     &                     ",    (ir,iz) = (",i3,",",i3,")",
     &                     ",   t =",i3,"   ang =",i3,a1)')
     &                     cha, inum, ir, iz, it, ia, cha
               else if( ittty(m) .eq. 0 .and. itaty(m) .ne. 0 ) then
                  write(iot,'(/a1,"no. =",i3,
     &                     ",    (ir,iz) = (",i3,",",i3,")",
     &                     ",  ang =",i3,a1)')
     &                     cha, inum, ir, iz, ia, cha
               else if( ittty(m) .ne. 0 .and. itaty(m) .eq. 0 ) then
                  write(iot,'(/a1,"no. =",i3,
     &                     ",    (ir,iz) = (",i3,",",i3,")",
     &                     ",   t =",i3,a1)')
     &                     cha, inum, ir, iz, it, cha
               else
                  write(iot,'(/a1,"no. =",i3,
     &                     ",    (ir,iz) = (",i3,",",i3,")",
     &                     a1)')
     &                     cha, inum, ir, iz, cha
               end if
             else
               write(changelsub(1),'(a1,"no. =",i3)') cha,inum
               write(changelsub(2),'(",  part = ",a)') chp(ipi)(3:10)
               write(changelsub(3),'(",  e =",i3)') ie
               write(changelsub(4),'(",  t =",i3)') it
               write(changelsub(5),'(",  ang =",i3)') ia
               write(changelsub(6),'(",  r =",i3)') ir
               write(changelsub(7),'(",  z =",i3)') iz
               write(changelsub(8),'(a1)') cha

               if( ne .eq. 0 ) then
                 changelsub(3) = " "
               end if
               if( ittty(m) .eq. 0 ) then
                 changelsub(4) = " "
               end if
               if( itaty(m) .eq. 0 ) then
                 changelsub(5) = " "
               end if

               changelsub(3) = " "

               if(iloopmode .eq. 1) then
                 changelsub(3) = " "
               end if
               if(iloopmode .eq. 10) then
                 changelsub(4) = " "
               end if
               if(iloopmode .eq. 8 .or. iloopmode .eq. 9) then
                 changelsub(5) = " "
               end if
               if(iloopmode .eq. 6) then
                 changelsub(6) = " "
               end if
               if(iloopmode .eq. 5) then
                 changelsub(7) = " "
               end if
                angeltitle = trim(changelsub(1))//trim(changelsub(2))//
     &                       trim(changelsub(3))//trim(changelsub(4))//
     &                       trim(changelsub(5))//trim(changelsub(6))//
     &                       trim(changelsub(7))//trim(changelsub(8))
                write(iot,'(/a)') trim(angeltitle)
          end if

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn


                  write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "   vol  &=&",1pe13.4," [cm^3]"/
     &                     "  rmin  &=&",1pe13.4," [cm]"/
     &                     "  rmax  &=&",1pe13.4," [cm]"/
     &                     "  zmin  &=&",1pe13.4," [cm]"/
     &                     "  zmax  &=&",1pe13.4," [cm]")')
     &                     yen,
     &                     vm, rm(ir), rm(ir+nrstepi),
     &                     zm(iz), zm(iz+nzstepi)
               if(iloopmode .ge. 1) then
                     write(iot,'("  part  &=& ",a)') chp(ipi)(3:10)
               endif
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &                     "  tmin  &=&",1pe13.4," [nsec]"/
     &                     "  tmax  &=&",1pe13.4," [nsec]")')
     &                     tb(it), tb(it+ntstepi)
               end if
               if( itaty(m) .ne. 0 ) then
                  write(iot,'(
     &                     "  amin  &=&",1pe13.4,/
     &                     "  amax  &=&",1pe13.4)')
     &                     ab(ia), ab(ia+nastepi)
               end if

                  write(iot,'("e:")')

            end do

            end do
            end do
            end do
            end do

*-----------------------------------------------------------------------
*        time axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 11 ) then

               inum = 0

            do iri = 1, nr, nrstepi
            do izi = 1, nz, nzstepi
            do iei = 1, ne, nestepi
            do iai = 1, na, nastepi
            do ipi = 1, np, npstepi
               ir = iri
               iz = izi
               ie = iei
               ia = iai
               ip = ipi

               inum = inum + 1

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

               if(iloopmode .ge. 1) then
                  write(iot,'("#   ",a)') chp(ipi)(3:10)
               endif

               write(iot,'("#   no. =",i3,3x,
     &         "ir  =",i3,3x,
     &         "iz  =",i3,3x,
     &         "ie  =",i3,3x,/
     &         "#   r = (",1p1e13.4,"  -",1p1e13.4,"  )"/
     &         "#   z = (",1p1e13.4,"  -",1p1e13.4,"  )"/
     &         "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     inum, ir, iz, ie,
     &                     rm(ir), rm(ir+nrstepi),
     &                     zm(iz), zm(iz+nzstepi),
     &                     eb(ie), eb(ie+nestepi)

               if( itaty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  ia =",i3,/
     &            "# ",a3," = (",
     &                           1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     ia, aname, ab(ia), ab(ia+nastepi)
               end if

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: Time [nsec]")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: Number ",a32)') hsunit(itunt(m))

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

               if( ittty(m) .eq. 3 .or. ittty(m) .eq. 5 ) then

                  write(iot,'( "p: xlog ylog afac(0.8) form(0.9)")')

               else

                  write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')

               end if

               if( itanl(m) .gt. 0 ) then

                  write(iot,'( "p: ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               end if

               if( itsans(m) .gt. 0 ) then

                  call write_sangel(iot,m,0)

               end if

       include 'samepage_include/petarz_t.inc'

             if(iloopmode .eq. 0) then
               if( itaty(m) .ne. 0 ) then
                  write(iot,'(/a1,"no. =",i3,
     &                     ",    (ir,iz) = (",i3,",",i3,")",
     &                     ",    e =",i3,
     &                     ",  ang =",i3,a1)')
     &                     cha, inum, ir, iz, ie, ia, cha
               else
                  write(iot,'(/a1,"no. =",i3,
     &                     ",    (ir,iz) = (",i3,",",i3,")",
     &                     ",    e =",i3,a1)')
     &                     cha, inum, ir, iz, ie, cha
               end if

             else
               write(changelsub(1),'(a1,"no. =",i3)') cha,inum
               write(changelsub(2),'(",  part = ",a)') chp(ipi)(3:10)
               write(changelsub(3),'(",  e =",i3)') ie
               write(changelsub(4),'(",  t =",i3)') it
               write(changelsub(5),'(",  ang =",i3)') ia
               write(changelsub(6),'(",  r =",i3)') ir
               write(changelsub(7),'(",  z =",i3)') iz
               write(changelsub(8),'(a1)') cha

               if( ne .eq. 0 ) then
                 changelsub(3) = " "
               end if
               if( ittty(m) .eq. 0 ) then
                 changelsub(4) = " "
               end if
               if( itaty(m) .eq. 0 ) then
                 changelsub(5) = " "
               end if

               changelsub(4) = " "

               if(iloopmode .eq. 1) then
                 changelsub(3) = " "
               end if
               if(iloopmode .eq. 10) then
                 changelsub(4) = " "
               end if
               if(iloopmode .eq. 8 .or. iloopmode .eq. 9) then
                 changelsub(5) = " "
               end if
               if(iloopmode .eq. 6) then
                 changelsub(6) = " "
               end if
               if(iloopmode .eq. 5) then
                 changelsub(7) = " "
               end if
                angeltitle = trim(changelsub(1))//trim(changelsub(2))//
     &                       trim(changelsub(3))//trim(changelsub(4))//
     &                       trim(changelsub(5))//trim(changelsub(6))//
     &                       trim(changelsub(7))//trim(changelsub(8))
                write(iot,'(/a)') trim(angeltitle)
          end if

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn


                if( ite2l(m).eq.0 ) then  ! ccse 2022/09/30
                  if( iMeVperu.eq.0 ) then
                     write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "   vol  &=&",1pe13.4," [cm^3]"/
     &                     "  rmin  &=&",1pe13.4," [cm]"/
     &                     "  rmax  &=&",1pe13.4," [cm]"/
     &                     "  zmin  &=&",1pe13.4," [cm]"/
     &                     "  zmax  &=&",1pe13.4," [cm]"/
     &                     "  emin  &=&",1pe13.4," [MeV]"/
     &                     "  emax  &=&",1pe13.4," [MeV]")')
     &                     yen,
     &                     vm, rm(ir), rm(ir+nrstepi),
     &                     zm(iz), zm(iz+nzstepi),
     &                     eb(ie), eb(ie+nestepi)
                  else
                     write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "   vol  &=&",1pe13.4," [cm^3]"/
     &                     "  rmin  &=&",1pe13.4," [cm]"/
     &                     "  rmax  &=&",1pe13.4," [cm]"/
     &                     "  zmin  &=&",1pe13.4," [cm]"/
     &                     "  zmax  &=&",1pe13.4," [cm]"/
     &                     "  emin  &=&",1pe13.4," [MeV/n]"/
     &                     "  emax  &=&",1pe13.4," [MeV/n]")')
     &                     yen,
     &                     vm, rm(ir), rm(ir+nrstepi),
     &                     zm(iz), zm(iz+nzstepi),
     &                     eb(ie), eb(ie+nestepi)
                  end if
                else if( ite2l(m).eq.1 ) then
                     write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "   vol  &=&",1pe13.4," [cm^3]"/
     &                     "  rmin  &=&",1pe13.4," [cm]"/
     &                     "  rmax  &=&",1pe13.4," [cm]"/
     &                     "  zmin  &=&",1pe13.4," [cm]"/
     &                     "  zmax  &=&",1pe13.4," [cm]"/
     &                     "  emin  &=&",1pe13.4," [keV/um]"/
     &                     "  emax  &=&",1pe13.4," [keV/um]")')
     &                     yen,
     &                     vm, rm(ir), rm(ir+nrstepi),
     &                     zm(iz), zm(iz+nzstepi),
     &                     eb(ie), eb(ie+nestepi)
                end if
               if(iloopmode .ge. 1) then
                     write(iot,'("  part  &=& ",a)') chp(ipi)(3:10)
               endif
               if( itaty(m) .ne. 0 ) then
                  write(iot,'(
     &                     "  amin  &=&",1pe13.4,/
     &                     "  amax  &=&",1pe13.4)')
     &                     ab(ia), ab(ia+nastepi)
               end if

                  write(iot,'("e:")')

            end do

            end do
            end do
            end do
            end do

*-----------------------------------------------------------------------
*        angle axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 12 .or.
     &            itaxs(m,iax) .eq. 13 ) then

               inum = 0

            do iri = 1, nr, nrstepi
            do izi = 1, nz, nzstepi
            do iei = 1, ne, nestepi
            do iti = 1, nt, ntstepi
            do ipi = 1, np, npstepi
               ir = iri
               iz = izi
               ie = iei
               it = iti
               ip = ipi

               inum = inum + 1

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

               if(iloopmode .ge. 1) then
                  write(iot,'("#   ",a)') chp(ipi)(3:10)
               endif

               write(iot,'("#   no. =",i3,3x,
     &         "ir  =",i3,3x,
     &         "iz  =",i3,3x,
     &         "ie  =",i3,3x,/
     &         "#   r = (",1p1e13.4,"  -",1p1e13.4,"  )"/
     &         "#   z = (",1p1e13.4,"  -",1p1e13.4,"  )"/
     &         "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     inum, ir, iz, ie,
     &                     rm(ir), rm(ir+nrstepi),
     &                     zm(iz), zm(iz+nzstepi),
     &                     eb(ie), eb(ie+nestepi)

               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  it =",i3/
     &            "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                  it, tb(it), tb(it+ntstepi)
               end if

               if( itaxl(m) .eq. 0 ) then

                  if( itaxs(m,iax) .eq. 12 ) then
                     write(iot,'(/"x: cos(",a1,"theta)")') yen
                  else
                     write(iot,'(/"x: ",a1,"theta  [deg]")') yen
                  end if

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: Number ",a32)') hsunit(itunt(m))

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

       include 'samepage_include/petarz_a.inc'

             if(iloopmode .eq. 0) then
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(/a1,"no. =",i3,
     &                     ",    (ir,iz) = (",i3,",",i3,")",
     &                     ",    e =",i3,
     &                     ",    t =",i3,a1)')
     &                     cha, inum, ir, iz, ie, it, cha
               else
                  write(iot,'(/a1,"no. =",i3,
     &                     ",    (ir,iz) = (",i3,",",i3,")",
     &                     ",    e =",i3,a1)')
     &                     cha, inum, ir, iz, ie, cha
               end if

             else
               write(changelsub(1),'(a1,"no. =",i3)') cha,inum
               write(changelsub(2),'(",  part = ",a)') chp(ipi)(3:10)
               write(changelsub(3),'(",  e =",i3)') ie
               write(changelsub(4),'(",  t =",i3)') it
               write(changelsub(5),'(",  ang =",i3)') ia
               write(changelsub(6),'(",  r =",i3)') ir
               write(changelsub(7),'(",  z =",i3)') iz
               write(changelsub(8),'(a1)') cha

               if( ne .eq. 0 ) then
                 changelsub(3) = " "
               end if
               if( ittty(m) .eq. 0 ) then
                 changelsub(4) = " "
               end if
               if( itaty(m) .eq. 0 ) then
                 changelsub(5) = " "
               end if

               changelsub(5) = " "

               if(iloopmode .eq. 1) then
                 changelsub(3) = " "
               end if
               if(iloopmode .eq. 10) then
                 changelsub(4) = " "
               end if
               if(iloopmode .eq. 8 .or. iloopmode .eq. 9) then
                 changelsub(5) = " "
               end if
               if(iloopmode .eq. 6) then
                 changelsub(6) = " "
               end if
               if(iloopmode .eq. 5) then
                 changelsub(7) = " "
               end if
                angeltitle = trim(changelsub(1))//trim(changelsub(2))//
     &                       trim(changelsub(3))//trim(changelsub(4))//
     &                       trim(changelsub(5))//trim(changelsub(6))//
     &                       trim(changelsub(7))//trim(changelsub(8))
                write(iot,'(/a)') trim(angeltitle)
          end if

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

                if( ite2l(m).eq.0 ) then  ! ccse 2022/09/30
                  if( iMeVperu.eq.0 ) then
                     write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "   vol  &=&",1pe13.4," [cm^3]"/
     &                     "  rmin  &=&",1pe13.4," [cm]"/
     &                     "  rmax  &=&",1pe13.4," [cm]"/
     &                     "  zmin  &=&",1pe13.4," [cm]"/
     &                     "  zmax  &=&",1pe13.4," [cm]"/
     &                     "  emin  &=&",1pe13.4," [MeV]"/
     &                     "  emax  &=&",1pe13.4," [MeV]")')
     &                     yen,
     &                     vm, rm(ir), rm(ir+nrstepi),
     &                     zm(iz), zm(iz+nzstepi),
     &                     eb(ie), eb(ie+nestepi)
                  else
                     write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "   vol  &=&",1pe13.4," [cm^3]"/
     &                     "  rmin  &=&",1pe13.4," [cm]"/
     &                     "  rmax  &=&",1pe13.4," [cm]"/
     &                     "  zmin  &=&",1pe13.4," [cm]"/
     &                     "  zmax  &=&",1pe13.4," [cm]"/
     &                     "  emin  &=&",1pe13.4," [MeV/n]"/
     &                     "  emax  &=&",1pe13.4," [MeV/n]")')
     &                     yen,
     &                     vm, rm(ir), rm(ir+nrstepi),
     &                     zm(iz), zm(iz+nzstepi),
     &                     eb(ie), eb(ie+nestepi)
                  end if
                else if( ite2l(m).eq.1 ) then
                     write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "   vol  &=&",1pe13.4," [cm^3]"/
     &                     "  rmin  &=&",1pe13.4," [cm]"/
     &                     "  rmax  &=&",1pe13.4," [cm]"/
     &                     "  zmin  &=&",1pe13.4," [cm]"/
     &                     "  zmax  &=&",1pe13.4," [cm]"/
     &                     "  emin  &=&",1pe13.4," [keV/um]"/
     &                     "  emax  &=&",1pe13.4," [keV/um]")')
     &                     yen,
     &                     vm, rm(ir), rm(ir+nrstepi),
     &                     zm(iz), zm(iz+nzstepi),
     &                     eb(ie), eb(ie+nestepi)
                end if
               if(iloopmode .ge. 1) then
                     write(iot,'("  part  &=& ",a)') chp(ipi)(3:10)
               endif

               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &                     "  tmin  &=&",1pe13.4," [nsec]"/
     &                     "  tmax  &=&",1pe13.4," [nsec]")')
     &                     tb(it), tb(it+ntstepi)
               end if

                  write(iot,'("e:")')

            end do

            end do
            end do
            end do
            end do

*-----------------------------------------------------------------------
*        r axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 6 ) then

               inum = 0

            do iei = 1, ne, nestepi
            do iti = 1, nt, ntstepi
            do iai = 1, na, nastepi
            do izi = 1, nz, nzstepi
            do ipi = 1, np, npstepi
               ie = iei
               it = iti
               ia = iai
               iz = izi
               ip = ipi

               inum = inum + 1

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

               if(iloopmode .ge. 1) then
                  write(iot,'("#   ",a)') chp(ipi)(3:10)
               endif

                  write(iot,'("#   no. =",i3,3x,
     &            "ie  =",i3,3x,
     &            "iz  =",i3,/
     &            "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )"/
     &            "#   z = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        inum, ie, iz,
     &                        eb(ie), eb(ie+nestepi),
     &                        zm(iz), zm(iz+nzstepi)
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  it =",i3/
     &            "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                  it, tb(it), tb(it+ntstepi)
               end if

               if( itaty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  ia =",i3,/
     &            "# ",a3," = (",
     &                           1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     ia, aname, ab(ia), ab(ia+nastepi)
               end if

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: r [cm]")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: Number ",a32)') hsunit(itunt(m))

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

               if( itrty(m) .eq. 3 .or. itrty(m) .eq. 5 ) then

                  write(iot,'( "p: xlog ylog afac(0.8) form(0.9)")')

               else

                  write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')

               end if

               if( itanl(m) .gt. 0 ) then

                  write(iot,'( "p: ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               end if

               if( itsans(m) .gt. 0 ) then

                  call write_sangel(iot,m,0)

               end if

       include 'samepage_include/petarz_r.inc'

*-----------------------------------------------------------------------

             if(iloopmode .eq. 0) then
               if( ittty(m) .ne. 0 .and. itaty(m) .ne. 0 ) then
                  write(iot,'(/a1,"no. =",i3,
     &                     ",  e =",i3,
     &                     ",  z =",i3,
     &                     ",  t =",i3,"   ang =",i3,a1)')
     &                     cha, inum, ie, iz, it, ia, cha
               else if( ittty(m) .eq. 0 .and. itaty(m) .ne. 0 ) then
                  write(iot,'(/a1,"no. =",i3,
     &                     ",  e =",i3,
     &                     ",  z =",i3,
     &                     ",  ang =",i3,a1)')
     &                     cha, inum, ie, iz, ia, cha
               else if( ittty(m) .ne. 0 .and. itaty(m) .eq. 0 ) then
                  write(iot,'(/a1,"no. =",i3,
     &                     ",  e =",i3,
     &                     ",  z =",i3,
     &                     ",  t =",i3,a1)')
     &                     cha, inum, ie, iz, it, cha
               else
                  write(iot,'(/a1,"no. =",i3,
     &                     ",   e =",i3,
     &                     ",   z =",i3,a1)')
     &                     cha, inum, ie, iz, cha
               end if
             else
               write(changelsub(1),'(a1,"no. =",i3)') cha,inum
               write(changelsub(2),'(",  part = ",a)') chp(ipi)(3:10)
               write(changelsub(3),'(",  e =",i3)') ie
               write(changelsub(4),'(",  t =",i3)') it
               write(changelsub(5),'(",  ang =",i3)') ia
               write(changelsub(6),'(",  r =",i3)') ir
               write(changelsub(7),'(",  z =",i3)') iz
               write(changelsub(8),'(a1)') cha

               if( ne .eq. 0 ) then
                 changelsub(3) = " "
               end if
               if( ittty(m) .eq. 0 ) then
                 changelsub(4) = " "
               end if
               if( itaty(m) .eq. 0 ) then
                 changelsub(5) = " "
               end if

               changelsub(6) = " "

               if(iloopmode .eq. 1) then
                 changelsub(3) = " "
               end if
               if(iloopmode .eq. 10) then
                 changelsub(4) = " "
               end if
               if(iloopmode .eq. 8 .or. iloopmode .eq. 9) then
                 changelsub(5) = " "
               end if
               if(iloopmode .eq. 6) then
                 changelsub(6) = " "
               end if
               if(iloopmode .eq. 5) then
                 changelsub(7) = " "
               end if
                angeltitle = trim(changelsub(1))//trim(changelsub(2))//
     &                       trim(changelsub(3))//trim(changelsub(4))//
     &                       trim(changelsub(5))//trim(changelsub(6))//
     &                       trim(changelsub(7))//trim(changelsub(8))
                write(iot,'(/a)') trim(angeltitle)
          end if

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn


                if( ite2l(m).eq.0 ) then  ! ccse 2022/09/30
                  if( iMeVperu.eq.0 ) then
                     write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "  emin  &=&",1pe13.4," [MeV]"/
     &                     "  emax  &=&",1pe13.4," [MeV]"/
     &                     "  zmin  &=&",1pe13.4," [cm]"/
     &                     "  zmax  &=&",1pe13.4," [cm]")')
     &                     yen,
     &                     eb(ie), eb(ie+nestepi),
     &                     zm(iz), zm(iz+nzstepi)
                  else
                     write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "  emin  &=&",1pe13.4," [MeV/n]"/
     &                     "  emax  &=&",1pe13.4," [MeV/n]"/
     &                     "  zmin  &=&",1pe13.4," [cm]"/
     &                     "  zmax  &=&",1pe13.4," [cm]")')
     &                     yen,
     &                     eb(ie), eb(ie+nestepi),
     &                     zm(iz), zm(iz+nzstepi)
                  end if
                else if( ite2l(m).eq.1 ) then
                     write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "  emin  &=&",1pe13.4," [keV/um]"/
     &                     "  emax  &=&",1pe13.4," [keV/um]"/
     &                     "  zmin  &=&",1pe13.4," [cm]"/
     &                     "  zmax  &=&",1pe13.4," [cm]")')
     &                     yen,
     &                     eb(ie), eb(ie+nestepi),
     &                     zm(iz), zm(iz+nzstepi)
                end if
               if(iloopmode .ge. 1) then
                     write(iot,'("  part  &=& ",a)') chp(ipi)(3:10)
               endif
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &                     "  tmin  &=&",1pe13.4," [nsec]"/
     &                     "  tmax  &=&",1pe13.4," [nsec]")')
     &                     tb(it), tb(it+ntstepi)
               end if
               if( itaty(m) .ne. 0 ) then
                  write(iot,'(
     &                     "  amin  &=&",1pe13.4,/
     &                     "  amax  &=&",1pe13.4)')
     &                     ab(ia), ab(ia+nastepi)
               end if

                  write(iot,'("e:")')

            end do

            end do
            end do
            end do
            end do

*-----------------------------------------------------------------------
*        z axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 5 ) then

               inum = 0

            do iei = 1, ne, nestepi
            do iti = 1, nt, ntstepi
            do iai = 1, na, nastepi
            do iri = 1, nr, nrstepi
            do ipi = 1, np, npstepi
               ie = iei
               it = iti
               ia = iai
               ir = iri
               ip = ipi

               inum = inum + 1

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

               if(iloopmode .ge. 1) then
                  write(iot,'("#   ",a)') chp(ipi)(3:10)
               endif


                  write(iot,'("#   no. =",i3,3x,
     &            "ie  =",i3,3x,
     &            "ir  =",i3,/
     &            "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )"/
     &            "#   r = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        inum, ie, ir,
     &                        eb(ie), eb(ie+nestepi),
     &                        rm(ir), rm(ir+nrstepi)
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  it =",i3/
     &            "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                  it, tb(it), tb(it+ntstepi)
               end if

               if( itaty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  ia =",i3,/
     &            "# ",a3," = (",
     &                           1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     ia, aname, ab(ia), ab(ia+nastepi)
               end if

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: z [cm]")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: Number ",a32)') hsunit(itunt(m))

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

               if( itzty(m) .eq. 3 .or. itzty(m) .eq. 5 ) then

                  write(iot,'( "p: xlog ylog afac(0.8) form(0.9)")')

               else

                  write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')

               end if

               if( itanl(m) .gt. 0 ) then

                  write(iot,'( "p: ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               end if

               if( itsans(m) .gt. 0 ) then

                  call write_sangel(iot,m,0)

               end if

       include 'samepage_include/petarz_z.inc'

*-----------------------------------------------------------------------

             if(iloopmode .eq. 0) then
               if( ittty(m) .ne. 0 .and. itaty(m) .ne. 0 ) then
                  write(iot,'(/a1,"no. =",i3,
     &                     ",  e =",i3,
     &                     ",  r =",i3,
     &                     ",  t =",i3,"   ang =",i3,a1)')
     &                     cha, inum, ie, ir, it, ia, cha
               else if( ittty(m) .eq. 0 .and. itaty(m) .ne. 0 ) then
                  write(iot,'(/a1,"no. =",i3,
     &                     ",  e =",i3,
     &                     ",  r =",i3,
     &                     ",  ang =",i3,a1)')
     &                     cha, inum, ie, ir, ia, cha
               else if( ittty(m) .ne. 0 .and. itaty(m) .eq. 0 ) then
                  write(iot,'(/a1,"no. =",i3,
     &                     ",  e =",i3,
     &                     ",  r =",i3,
     &                     ",  t =",i3,a1)')
     &                     cha, inum, ie, ir, it, cha
               else
                  write(iot,'(/a1,"no. =",i3,
     &                     ",  e =",i3,
     &                     ",  r =",i3,a1)')
     &                     cha, inum, ie, ir, cha
               end if
             else
               write(changelsub(1),'(a1,"no. =",i3)') cha,inum
               write(changelsub(2),'(",  part = ",a)') chp(ipi)(3:10)
               write(changelsub(3),'(",  e =",i3)') ie
               write(changelsub(4),'(",  t =",i3)') it
               write(changelsub(5),'(",  ang =",i3)') ia
               write(changelsub(6),'(",  r =",i3)') ir
               write(changelsub(7),'(",  z =",i3)') iz
               write(changelsub(8),'(a1)') cha

               if( ne .eq. 0 ) then
                 changelsub(3) = " "
               end if
               if( ittty(m) .eq. 0 ) then
                 changelsub(4) = " "
               end if
               if( itaty(m) .eq. 0 ) then
                 changelsub(5) = " "
               end if

               changelsub(7) = " "

               if(iloopmode .eq. 1) then
                 changelsub(3) = " "
               end if
               if(iloopmode .eq. 10) then
                 changelsub(4) = " "
               end if
               if(iloopmode .eq. 8 .or. iloopmode .eq. 9) then
                 changelsub(5) = " "
               end if
               if(iloopmode .eq. 6) then
                 changelsub(6) = " "
               end if
               if(iloopmode .eq. 5) then
                 changelsub(7) = " "
               end if
                angeltitle = trim(changelsub(1))//trim(changelsub(2))//
     &                       trim(changelsub(3))//trim(changelsub(4))//
     &                       trim(changelsub(5))//trim(changelsub(6))//
     &                       trim(changelsub(7))//trim(changelsub(8))
                write(iot,'(/a)') trim(angeltitle)
          end if

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn


                if( ite2l(m).eq.0 ) then  ! ccse 2022/09/30
                  if( iMeVperu.eq.0 ) then
                     write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "  emin  &=&",1pe13.4," [MeV]"/
     &                     "  emax  &=&",1pe13.4," [MeV]"/
     &                     "  rmin  &=&",1pe13.4," [cm]"/
     &                     "  rmax  &=&",1pe13.4," [cm]")')
     &                     yen,
     &                     eb(ie), eb(ie+nestepi),
     &                     rm(ir), rm(ir+nrstepi)
                  else
                     write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "  emin  &=&",1pe13.4," [MeV/n]"/
     &                     "  emax  &=&",1pe13.4," [MeV/n]"/
     &                     "  rmin  &=&",1pe13.4," [cm]"/
     &                     "  rmax  &=&",1pe13.4," [cm]")')
     &                     yen,
     &                     eb(ie), eb(ie+nestepi),
     &                     rm(ir), rm(ir+nrstepi)
                  end if
                else if( ite2l(m).eq.1 ) then
                     write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "  emin  &=&",1pe13.4," [keV/um]"/
     &                     "  emax  &=&",1pe13.4," [keV/um]"/
     &                     "  rmin  &=&",1pe13.4," [cm]"/
     &                     "  rmax  &=&",1pe13.4," [cm]")')
     &                     yen,
     &                     eb(ie), eb(ie+nestepi),
     &                     rm(ir), rm(ir+nrstepi)
                end if
               if(iloopmode .ge. 1) then
                     write(iot,'("  part  &=& ",a)') chp(ipi)(3:10)
               endif
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &                     "  tmin  &=&",1pe13.4," [nsec]"/
     &                     "  tmax  &=&",1pe13.4," [nsec]")')
     &                     tb(it), tb(it+ntstepi)
               end if
               if( itaty(m) .ne. 0 ) then
                  write(iot,'(
     &                     "  amin  &=&",1pe13.4,/
     &                     "  amax  &=&",1pe13.4)')
     &                     ab(ia), ab(ia+nastepi)
               end if

                  write(iot,'("e:")')

            end do

            end do
            end do
            end do
            end do

*-----------------------------------------------------------------------
*        rz axis (matrix)
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

            do ip = 1, np
            do ie = 1, ne
            do it = 1, nt
            do ia = 1, na

               inum = inum + 1

                  write(iot,'(/"#",78("-"))')

               if( inum .eq. 1 .or. ip.gt.np_mxang ) then ! frtati 2021/10/05
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

               if( rtfac(m) .lt. 0.d0 )
     &           write(iot,
     &           '( "#restart maximum for normalization:",1p1e14.7)')
     &           facmax(m)

               if ( ip.gt.np_mxang ) then
                 write(iot,'( " SKIPPAGE:")')
                 inum = inum - 1
               end if

                  write(iot,'("#   no. =",i3,3x,
     &            "ie  =",i3,3x,
     &            "part. = ",a8,/
     &            "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        inum, ie, chq(ip),
     &                        eb(ie), eb(ie+1)
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  it =",i3/
     &            "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                  it, tb(it), tb(it+1)
               end if

               if( itaty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  ia =",i3,/
     &            "# ",a3," = (",
     &                           1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     ia, aname, ab(ia), ab(ia+1)
               end if

*-----------------------------------------------------------------------

               if( ittty(m) .ne. 0 .and. itaty(m) .ne. 0 ) then
                  write(iot,'(/a1,"no. =",i3,
     &                     ",  e =",i3,
     &                     ",  t =",i3,"   ang =",i3,a1)')
     &                     cha, inum, ie, it, ia, cha
               else if( ittty(m) .eq. 0 .and. itaty(m) .ne. 0 ) then
                  write(iot,'(/a1,"no. =",i3,
     &                     ",  e =",i3,
     &                     ",  ang =",i3,a1)')
     &                     cha, inum, ie, ia, cha
               else if( ittty(m) .ne. 0 .and. itaty(m) .eq. 0 ) then
                  write(iot,'(/a1,"no. =",i3,
     &                     ",  e =",i3,
     &                     ",  t =",i3,a1)')
     &                     cha, inum, ie, it, cha
               else
                  write(iot,'(/a1,"no. =",i3,
     &                     ",  e =",i3,a1)')
     &                     cha, inum, ie, cha
               end if

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

*-----------------------------------------------------------------------

            if( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) then

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

                  write(iot,'("#  nr = ",i3,"   nz = ",i3)')
     &                         nr, nz

            if( ittwo(m) .ne. 4 ) then

                  write(iot,'( "# ( ( data(z,r), z = 1, nz ),",
     &                         " r = nr, 1, -1 )")')

            end if

*-----------------------------------------------------------------------

            if( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) then

                  write(iot,'(/a4," y = ",1p1g14.7," to ",
     &            1p1g14.7," by ",1p1g14.7," ; x = ",
     &            1p1g14.7," to ",1p1g14.7," by ",
     &            1p1g14.7," ;")') dc2,
     &            rm(nr) + rtrdl(m)/2.0, rm(1) + rtrdl(m)/2.0, rtrdl(m),
     &            zm(1) + rtzdl(m)/2.0, zm(nz) + rtzdl(m)/2.0, rtzdl(m)

!OBINATA(2012.9.5): output *.err
               write(iot,'(1p10e11.3)')
     &         ( ( tr(ip,ie,it,ia,ir,iz,ioe),
     &             iz = 1, nz ), ir = nr, 1, -1 )

*-----------------------------------------------------------------------

            else if( ittwo(m) .eq. 4 ) then

               write(iot,'(/"# r          z        ",
     &                      "  number     r.err")')

               do iz = 1, nz
               do ir = 1, nr

                  write(iot,'(1p3e11.3,0pf8.4)')
     &               rm(ir)  + rtrdl(m)/2.0,
     &               zm(iz)  + rtzdl(m)/2.0,
     &               tr(ip,ie,it,ia,ir,iz,1), tr(ip,ie,it,ia,ir,iz,2)

               end do
               end do

*-----------------------------------------------------------------------

            else if( ittwo(m) .eq. 5 ) then

                  write(iot,'("#   r = ",1p1g14.7," to ",
     &            1p1g14.7," by ",1p1g14.7/
     &                        "#   z = ",1p1g14.7," to ",
     &            1p1g14.7," by ",1p1g14.7)')
     &            rm(1) + rtrdl(m)/2.0, rm(nr) + rtrdl(m)/2.0, rtrdl(m),
     &            zm(1) + rtzdl(m)/2.0, zm(nz) + rtzdl(m)/2.0, rtzdl(m)

                  write(iot,'(/4x,a3,4x,1p1000e11.3)')
     &            'r/z',( zm(iz) + rtzdl(m)/2.0, iz = 1, nz )

               do ir = nr, 1, -1

!OBINATA(2012.9.5): output *.err
                  write(iot,'(1p1000e11.3)')
     &            rm(ir) + rtrdl(m)/2.0,
     &            ( tr(ip,ie,it,ia,ir,iz,ioe), iz = 1, nz )

               end do

            end if

*-----------------------------------------------------------------------

      if( inolg .eq. 1 .and. ( ittwo(m) .eq. 3 .or. ittwo(m) .eq. 7 )
     &    .and. cmin .gt. 0.0 .and. cmax .gt. cmin .and.
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

         write(iot,'("y: Number ",a32)') hsunit(itunt(m))

      else

         write(iot,'("y: ",200a1)') (itazt(m)(i:i),i=1,itazl(m))

      end if

      if( izlog .eq. 0 ) write(iot,'( "p: ylin")')

      write(iot,'( "#",78("-"))')

      end if

*-----------------------------------------------------------------------

            if( inocm .eq. 1 ) then

                if( ite2l(m).eq.0 ) then  ! ccse 2022/09/30
                  if( iMeVperu.eq.0 ) then
                     write(iot,'(/"wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "  emin  &=&",1pe13.4," [MeV]"/
     &                     "  emax  &=&",1pe13.4," [MeV]"/
     &                     "  part. &=&   ",a8)')
     &                     yen, eb(ie), eb(ie+1),
     &                          chq(ip)
                  else
                     write(iot,'(/"wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "  emin  &=&",1pe13.4," [MeV/n]"/
     &                     "  emax  &=&",1pe13.4," [MeV/n]"/
     &                     "  part. &=&   ",a8)')
     &                     yen, eb(ie), eb(ie+1),
     &                          chq(ip)
                  end if
                else if( ite2l(m).eq.1 ) then
                     write(iot,'(/"wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "  emin  &=&",1pe13.4," [keV/um]"/
     &                     "  emax  &=&",1pe13.4," [keV/um]"/
     &                     "  part. &=&   ",a8)')
     &                     yen, eb(ie), eb(ie+1),
     &                          chq(ip)
                end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &                     "  tmin  &=&",1pe13.4," [nsec]"/
     &                     "  tmax  &=&",1pe13.4," [nsec]")')
     &                     tb(it), tb(it+1)
               end if
               if( itaty(m) .ne. 0 ) then
                  write(iot,'(
     &                     "  amin  &=&",1pe13.4,/
     &                     "  amax  &=&",1pe13.4)')
     &                     ab(ia), ab(ia+1)
               end if

                  write(iot,'("e:")')

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
            end do
            end do
            end do

*-----------------------------------------------------------------------

         end if

*-----------------------------------------------------------------------

            call prestart(m,iot) !OBINATA(2012.9.11)

            close(iot)

         if( iteps(m) .ne. 0 ) then

            open(iot, file = fname, status = 'unknown' )
            call a_angel(idasa,fname)

         end if

        end do

      end do

      include 'samepage_include/samepage999.inc'

*-----------------------------------------------------------------------

      deallocate (vl_r,vl_z)

      return
      end


************************************************************************
*                                                                      *
      subroutine tpdctxyz(ncol,m,np,nx,ny,nz,ne,nt,na,nm,mt,nl,lt,
     &                    xm,ym,zm,eb,tb,ab,tr,trEVENT,
     &                    itrmax,itrmin)
*                                                                      *
*       product tally in xyz scoring mesh                              *
*       last modified by K.Niita on 2011/05/17                         *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*      ncol  ..... reaction type                                       *
*              4 : source                                              *
*             13 : nuclear reaction (n,x)                              *
*             14 : nuclear reaction (n,n'x)                            *
*                                                                      *
*        jcoll : =  0, nothing happen                                  *
*                =  1, Hydrogen collisions                             *
*                =  2, Particle Decays                                 *
*                =  3, Elastic collisions                              *
*                =  4, Nuclear collisions                              *
*                =  5, Heavy Ion collisions                            *
*                =  6, Neutron reactions by data                       *
*                =  7, Photon reactions by data                        *
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
      use MMBANKMOD !FURUTA
      use partmod, only: itmxpt, itpan, itpat, jtpat ! frtati 2021/10/05
*-----------------------------------------------------------------------

      implicit double precision( a-h, o-z )

*-----------------------------------------------------------------------

      include 'param00.inc'
      include 'param.inc'

      parameter ( pi  = 3.141592653589793d0 )

*-----------------------------------------------------------------------

      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
!$OMP THREADPRIVATE(/icomon/)
      common /jcomon/ nabov,nobch,nocas,nomax
!$OMP THREADPRIVATE(/jcomon/)
      common /wtsave/ oldwt
!$OMP THREADPRIVATE(/wtsave/)
      common /mathzn/ mathz, mathn, jcoll, kcoll
!$OMP THREADPRIVATE(/mathzn/)
      common /kmat1d/ idmn(0:kvlmax), idnm(kvmmax)
      common /pnint/  ipnint

*-----------------------------------------------------------------------

      common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                itout(itlmax), ittwo(itlmax)

      common /tall12/ itman(itlmax), itmat(itlmax), itmct(itlmax),
     &                itnun(itlmax), itnuc(itlmax),
     &                itndz(itlmax), itndn(itlmax),
     &                itnkz(itlmax), itnkn(itlmax)
      common /tall22/ itmcn(itlmax), itmtn(itlmax), itmtt(itlmax)
      common /tall37/ itcnt(9,itlmax)
      common /tall46/ itmtr(itlmax,4), rtmtr(itlmax,13)
      common /tall48/ itmdp(itlmax,0:30), itmdf(itlmax)

      common /tall23/ itaty(itlmax), itanm(itlmax), itarg(itlmax),
     &                rtami(itlmax), rtama(itlmax), rtadl(itlmax)

      common /tall82/ itcnth(9,itlmax)
      common /tall74/ iprim(itlmax)

      common /tall50/ itlmt(itlmax), ite2l(itlmax)

*-----------------------------------------------------------------------

      common /clustt/ nclsts, iclusts(nnn)
!$OMP THREADPRIVATE(/clustt/)
      common /clustw/ jclusts(0:8,nnn),  qclusts(0:12,nnn)
!$OMP THREADPRIVATE(/clustw/)
      common /cntcls/ jcount(3,nnn)
!$OMP THREADPRIVATE(/cntcls/)

*-----------------------------------------------------------------------

      dimension   xm(nx+1)
      dimension   ym(ny+1)
      dimension   zm(nz+1)
      dimension   mt(nm)
      dimension   lt(nl)
      dimension   eb(ne+1)
      dimension   tb(nt+1)
      dimension   ab(na+1)
      dimension   tr(np,ne,nt,na,nx*ny*nz,2)
      dimension   trEVENT(np,ne,nt,na,nx*ny*nz) !OBINATA(2012.9.11): as Ct
      real(8),allocatable,save:: tr0(:,:,:,:,:) !OBINATA(2012.9.11): as C
      dimension   itrmax(7),itrmin(7)

      dimension ips(itmxpt) ! frtati 2021/10/05 6 -> itmxpt
      dimension ncntt(3)
      dimension dmpd(30)

*-----------------------------------------------------------------------

      common /stat / istdev, irestart, ireschk
      common /cparm/ maxbch,maxcas

*-----------------------------------------------------------------------
      integer idmpomp !FURUTA20150427
      common /idmpomp0/idmpomp !FURUTA20150427
*-----------------------------------------------------------------------
      common /tstara/ atmrc(10,8)      ! Ogawa 2023/7/18, (10,7) -> (10,8) for plasmon
!$OMP THREADPRIVATE(/tstara/)
      common /tpdcta/ nflumu
!$OMP THREADPRIVATE(/tpdcta/)
      common /trskip/ ntrskip
!$OMP THREADPRIVATE(/trskip/)

      integer           iMeVperu
      common /cMeVperu/ iMeVperu

      common /cntmx/ ncntmx(3), nctmxsr(3)
!$OMP THREADPRIVATE(/cntmx/)

! sumover
      dimension ldo(2,7)

*-----------------------------------------------------------------------

      icf(ix,iy,iz) = ix + ( iy - 1 ) * nx + ( iz - 1 ) * nx * ny

*-----------------------------------------------------------------------

      if (istdev .eq. 2) then

        call readitrminmax7(itrmin,itrmax,(/ np,ne,nt,na,nx,ny,nz /),
     &                      mnp,mne,mnt,mna,mnx,mny,mnz,
     &                      mxp,mxe,mxt,mxa,mxx,mxy,mxz)

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
* OBINATA(2012.9.11): change tr(,,,,3) to trEVENT(,,,,)
*-----------------------------------------------------------------------

         if (( ncol .eq. 0 .or. ncol .eq. 4 )
     &                              .and. istdev .eq. 2) then
           if ((nocas.gt.1.or.ncol.eq.0) .and. ihistcount.ne.1 ) then

             do iz = mnz,mxz
             do iy = mny,mxy
             do ix = mnx,mxx
             do ia = mna,mxa
             do it = mnt,mxt
             do ie = mne,mxe
               ixyz = icf(ix,iy,iz)
               tr(:,ie,it,ia,ixyz,1) = tr(:,ie,it,ia,ixyz,1)
     &                               + trEVENT(:,ie,it,ia,ixyz)
               tr(:,ie,it,ia,ixyz,2) = tr(:,ie,it,ia,ixyz,2)
     &                               + trEVENT(:,ie,it,ia,ixyz) ** 2
             enddo
             enddo
             enddo
             enddo
             enddo
             enddo

! sumover
             ldo(1,1) = 1
             ldo(2,1) = np
             ldo(1,2) = mne
             ldo(2,2) = mxe
             ldo(1,3) = mnt
             ldo(2,3) = mxt
             ldo(1,4) = mna
             ldo(2,4) = mxa
             ldo(1,5) = mnx
             ldo(2,5) = mxx
             ldo(1,6) = mny
             ldo(2,6) = mxy
             ldo(1,7) = mnz
             ldo(2,7) = mxz
             call tpdctxyz_sumover(m,1,
     &                   np,  ne, nt, na, nx, ny, nz, trEVENT,ldo)

           end if

           do iz = mnz,mxz
           do iy = mny,mxy
           do ix = mnx,mxx
           do ia = mna,mxa
           do it = mnt,mxt
           do ie = mne,mxe
             ixyz = icf(ix,iy,iz)
             trEVENT(:,ie,it,ia,ixyz) = 0
           enddo
           enddo
           enddo
           enddo
           enddo
           enddo

           call resetitrminmax(itrmin,itrmax,7,(/np,ne,nt,na,nx,ny,nz/))

         end if

*-----------------------------------------------------------------------
*        end of batch ( in case of istdev = 1 )
*
* OBINATA(2012.9.11): modificate for thread parallel
*-----------------------------------------------------------------------

         if ( ncol .eq. 0 .and. istdev .eq. 1) then
!$OMP MASTER
             allocate( tr0(np,ne,nt,na,nx*ny*nz) )
             tr0(:,:,:,:,:) = 0.d0
!$OMP END MASTER
!$OMP BARRIER
!$OMP CRITICAL (tpdctxyz_crit_ist1)
             tr0(:,:,:,:,:) = tr0(:,:,:,:,:) + trEVENT(:,:,:,:,:)
!$OMP END CRITICAL (tpdctxyz_crit_ist1)
!$OMP BARRIER
!$OMP MASTER
             tr(:,:,:,:,:,1) = tr(:,:,:,:,:,1) + tr0(:,:,:,:,:) / maxcas
             tr(:,:,:,:,:,2) = tr(:,:,:,:,:,2)
     &                     + ( tr0(:,:,:,:,:) / maxcas ) ** 2

! sumover
           ldo(1,1) = 1
           ldo(2,1) = np
           ldo(1,2) = 1
           ldo(2,2) = ne
           ldo(1,3) = 1
           ldo(2,3) = nt
           ldo(1,4) = 1
           ldo(2,4) = na
           ldo(1,5) = 1
           ldo(2,5) = nx
           ldo(1,6) = 1
           ldo(2,6) = ny
           ldo(1,7) = 1
           ldo(2,7) = nz
             call tpdctxyz_sumover(m,maxcas,
     &                   np,  ne, nt, na, nx, ny, nz, tr0, ldo)


             deallocate( tr0 )
!$OMP END MASTER

           trEVENT(:,:,:,:,:) = 0

         end if

*-----------------------------------------------------------------------
*        out put : check of ncol, jcoll and kcoll
*-----------------------------------------------------------------------

         if( itout(m) .eq. 1 .and. ncol .ne. 4 ) return
         if( itout(m) .gt. 1 .and.
     &        ncol .ne. 13 .and. ncol .ne. 14 ) return


*-----------------------------------------------------------------------

         if( itout(m) .eq. 2 .and. ! output = nuclear (See subroutine tprodct)
     &       jcoll .ne. 1  .and. jcoll .ne. 3  .and.
     &       jcoll .ne. 4  .and. jcoll .ne. 5  .and.
     &       jcoll .ne. 6  .and. jcoll .ne. 9  .and.      ! ccse 2023/03/29
     &       jcoll .ne. 10 .and. jcoll .ne. 15 .and.
     &       jcoll .ne. 16 .and. jcoll .ne. 17 .and.      ! S.Abe 2017/02/08
     &       jcoll .ne. 19 .and. jcoll .ne. 20 ) return   ! y.sakaki 2024/01

         if( itout(m) .eq. 3 .and. jcoll .ne. 2 ) return
         if( itout(m) .eq. 3 .and.
     &       jcoll .eq. 2 .and. ityp .eq. 13 ) return     ! S.Abe 2017/02/08

         if( itout(m) .eq. 4 .and.
     &        kcoll .ne. 1 .and. kcoll .ne. 5 ) return

         if( itout(m) .eq. 5 .and.
     &       jcoll .ne. 3 .and. kcoll .ne. 3 ) return

         if( itout(m) .eq. 6 .and. ! output = nonela (See subroutine tprodct)
     &       jcoll .ne. 1  .and.
     &       jcoll .ne. 4  .and. jcoll .ne. 5  .and.
     &       jcoll .ne. 6  .and. jcoll .ne. 9  .and.      ! ccse 2023/03/29
     &       jcoll .ne. 10 .and. jcoll .ne. 15 .and.
     &       jcoll .ne. 16 .and. jcoll .ne. 17 .and.
     &       jcoll .ne. 19 .and. jcoll .ne. 20 ) return   ! y.sakaki 2024/01
         if( itout(m) .eq. 6 .and. kcoll .eq. 3 ) return

         sumatmrc = 0.d0
         do i = 1, 5
         do j = 1, 5    ! except mscat
            sumatmrc = sumatmrc + atmrc(i,j)
         enddo
         enddo
         if( itout(m) .eq. 7 .and.
     &       jcoll    .ne. 18 .and.
     &       ( jcoll .ne. 11 .and. sumatmrc .le. 0.d0 )
     &     ) return

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
*        for source
*-----------------------------------------------------------------------

         if( ncol .eq. 4 ) then

               iccol = 0
               npart = nomax

               ncntt(1) = ncnt(ibknct+1,no,ipomp+1)
               ncntt(2) = ncnt(ibknct+2,no,ipomp+1)
               ncntt(3) = ncnt(ibknct+3,no,ipomp+1)

*-----------------------------------------------------------------------
*        for reaction, check of reaction and mother
*-----------------------------------------------------------------------

         else

            if( nclsts .le. 0 ) return

               iccol = 1
               npart = nclsts

               ata   = dble( mathz + mathn )
               atz   = dble( mathz )
               mmas  = nint( ata )
               mchg  = nint( atz )

               ncntt(1) = jcount(1,1)
               ncntt(2) = jcount(2,1)
               ncntt(3) = jcount(3,1)

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

                     if( itmct(m) .gt. 0 ) goto 140
                     if( itmct(m) .lt. 0 ) return

                  end if

               end do

                  if( itmct(m) .gt. 0 ) return

            end if

  140       continue

*-----------------------------------------------------------------------
*           check elastic collision
*-----------------------------------------------------------------------

            if( ( jcoll .eq. 1 .or.
     &            jcoll .eq. 4 .or.
     &            jcoll .eq. 5 .or.
     &            jcoll .eq. 10 ) .and.
     &            npart .eq. 2 ) then

               if( jcoll .eq. 1 ) then

                  if( ( iclusts(1) .eq. 1 .and.
     &                  jclusts(7,2) .eq. ktyp .and.
     &                  qclusts(6,2) .eq. 0.d0 ) .or.
     &                ( iclusts(2) .eq. 1 .and.
     &                  jclusts(7,1) .eq. ktyp .and.
     &                  qclusts(6,1) .eq. 0.d0 ) ) then

                     if( itout(m) .ne. 2 .and. itout(m) .ne. 5 ) return

                  end if

               else

                  if( ( iclusts(1) .eq. 0 .and.
     &                  jclusts(1,1) .eq. mchg .and.
     &                  jclusts(2,1) .eq. mmas - mchg .and.
     &                  jclusts(7,2) .eq. ktyp .and.
     &                  qclusts(6,2) .eq. 0.d0 .and.
     &                  qclusts(6,1) .eq. 0.d0 ) .or.
     &                ( iclusts(2) .eq. 0 .and.
     &                  jclusts(1,2) .eq. mchg .and.
     &                  jclusts(2,2) .eq. mmas - mchg .and.
     &                  jclusts(7,1) .eq. ktyp .and.
     &                  qclusts(6,1) .eq. 0.d0 .and.
     &                  qclusts(6,2) .eq. 0.d0 ) ) then

                     if( itout(m) .ne. 2 .and. itout(m) .ne. 5 ) return

                  end if

               end if

            end if

*-----------------------------------------------------------------------

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
*        material for LET
*-----------------------------------------------------------------------

               if( itlmt(m) .gt. 0 ) then

                  lmat = idnm( itlmt(m) )

               else if( itlmt(m) .eq. 0 ) then

                  lmat = mat

               else

                  lmat = -idnm( -itlmt(m) )

               end if

*-----------------------------------------------------------------------
*        do loop for particles
*-----------------------------------------------------------------------

         do 100 j = 1, npart

*-----------------------------------------------------------------------
            if( ityp .eq. 7 .and. (jcoll.eq.2 .or. jcoll.eq.16) ) then
               if( itout(m) .eq. 7 .and. j .gt. nflumu ) goto 100
               if( itout(m) .ne. 7 .and. j .le. nflumu ) goto 100
            endif
            if( j .eq. ntrskip ) goto 100

            if( iccol .eq. 1 ) then

                  ipart = jclusts(3,j)
                  kpart = jclusts(7,j)
                  lpart = jclusts(8,j)
                  jpart = ichgf(ipart,kpart)
                  rpart = qclusts(5,j) * 1000.d0   ! ccse 2022/08/31

                  epart = qclusts(7,j)
                  tlw   = qclusts(8,j)
                  tpart = abs(qclusts(9,j))
                  xpart = qclusts(10,j)
                  ypart = qclusts(11,j)
                  zpart = qclusts(12,j)

                  upart = qclusts(1,j)
                  vpart = qclusts(2,j)
                  wpart = qclusts(3,j)

                  sxpat = spx(ibkspx+no,ipomp+1)
                  sypat = spy(ibkspy+no,ipomp+1)
                  szpat = spz(ibkspz+no,ipomp+1)
                  nzpat = jclusts(5,j)
                  nname = name(ibknam+no,ipomp+1)

            else

                  nj = no + j - 1

                  ipart = nty(ibknty+nj,ipomp+1)
                  kpart = nkf(ibknkf+nj,ipomp+1)
                  jpart = ichgf(ityp,ktyp)
                  rpart = rtyp   ! ccse 2022/08/31

                  epart = e(ibke+nj,ipomp+1)
                  tlw   = wt(ibkwt+nj,ipomp+1)
                  tpart = abs(t(ibkt+nj,ipomp+1))

                  xpart = x(ibkx+nj,ipomp+1)
                  ypart = y(ibky+nj,ipomp+1)
                  zpart = z(ibkz+nj,ipomp+1)

                  upart = u(ibku+nj,ipomp+1)
                  vpart = v(ibkv+nj,ipomp+1)
                  wpart = w(ibkw+nj,ipomp+1)

                  sxpat = spx(ibkspx+nj,ipomp+1)
                  sypat = spy(ibkspy+nj,ipomp+1)
                  szpat = spz(ibkspz+nj,ipomp+1)
                  nzpat = ctyp
                  nname = name(ibknam+nj,ipomp+1)

            end if

             if( kpart .eq. ktyp .and. iprim(m) .eq. 0 .and.   !exclude primary particle in 
     &      ( (jcoll .eq. 18 .or. jcoll .eq. 14 .or. jcoll .eq. 13 .or. ! Track-structure and EGS
     &         jcoll .eq. 3) .and. j .eq. 1 .or.! elastic scattering
     &         jcoll .eq. 10 .and. j .eq. 1 .and. npart .eq. 2 .or.! elastic scattering
     &         jcoll .ge. 6 .AND. jcoll .le. 9   .and. j .eq. npart !exclude primary sampled from X-section data
     &      ) ) goto 100 ! rejection of remaining projectile
             
            
            if( epart .le. 0.d0 ) goto 100   ! S.Abe 2017/02/08

*-----------------------------------------------------------------------
*           transform positions
*-----------------------------------------------------------------------

               call trnsxx(xpart,ypart,zpart,
     &                     xcc,ycc,zcc,itmtr(m,4))

*-----------------------------------------------------------------------
*           check position : out of rainge
*-----------------------------------------------------------------------

               if(  xcc .lt. xm(1) ) goto 100
               if(  xcc .ge. xm(nx+1) ) goto 100

               if(  ycc .lt. ym(1) ) goto 100
               if(  ycc .ge. ym(ny+1) ) goto 100

               if(  zcc .lt. zm(1) ) goto 100
               if(  zcc .ge. zm(nz+1) ) goto 100

*-----------------------------------------------------------------------
*           x-position
*-----------------------------------------------------------------------

               do i = 1, nx

                  if( xcc .ge. xm(i) .and.
     &                xcc .lt. xm(i+1) ) goto 36

               end do

   36          ix1 = i

*-----------------------------------------------------------------------
*           y-position
*-----------------------------------------------------------------------

               do i = 1, ny

                  if( ycc .ge. ym(i) .and.
     &                ycc .lt. ym(i+1) ) goto 37

               end do

   37          iy1 = i

*-----------------------------------------------------------------------
*           z-position
*-----------------------------------------------------------------------

               do i = 1, nz

                  if( zcc .ge. zm(i) .and.
     &                zcc .lt. zm(i+1) ) goto 38

               end do

   38          iz1 = i

*-----------------------------------------------------------------------
*           MeV -> MeV/n energy conversion ! 2018/3/7 Ogawa
*-----------------------------------------------------------------------

               if(iMeVperu.eq.1 .and. ipart.ge.15 .and. ipart.le.19)then
                    ebm = dble(kpart - kpart / 1000000 * 1000000)
               else
                    ebm = 1.d0
               end if

*-----------------------------------------------------------------------
*           check of particles
*-----------------------------------------------------------------------

            call pcheck(m,np,ipart,kpart+1000000000*lpart,jpart,ipn,ips)

               if( ipn .eq. 0 ) goto 100

*-----------------------------------------------------------------------
*           check of energy
*-----------------------------------------------------------------------

            if (ite2l(m) .eq. 0) then   ! ccse 2022/08/31

            if( epart .lt. eb(1)*ebm ) goto 100
            if( epart .ge. eb(ne+1)*ebm ) goto 100 ! 2014/12/16  Previously, 'epart' was mistakenly 'ec(ibkec+no)'

*-----------------------------------------------------------------------
*           energy
*-----------------------------------------------------------------------

            do i = 1, ne

               if( epart .ge. eb(i)*ebm .and.
     &             epart .lt. eb(i+1)*ebm ) goto 30

            end do

            else if (ite2l(m) .eq. 1) then   ! convert to energy to LET

               call dedxas(epart,dedx,lmat,ipart,kpart,jpart,rpart)
               dedx = dedx /10.0d0

            ! check of energy
            if( dedx .lt. eb(1) ) goto 100
            if( dedx .ge. eb(ne+1) ) goto 100 ! 2014/12/16  Previously, 'epart' was mistakenly 'ec(ibkec+no)'

            ! energy
            do i = 1, ne

               if( dedx .ge. eb(i) .and.
     &             dedx .lt. eb(i+1) ) goto 30

            end do

            end if

   30          ie1 = i

*-----------------------------------------------------------------------
*           check of time
*-----------------------------------------------------------------------

               if( tpart .lt. tb(1) ) goto 100
               if( tpart .ge. tb(nt+1) ) goto 100

*-----------------------------------------------------------------------
*           time
*-----------------------------------------------------------------------

            do i = 1, nt

               if( tpart .ge. tb(i) .and.
     &             tpart .lt. tb(i+1) ) goto 40

            end do

   40          it1 = i

*-----------------------------------------------------------------------
*        check angle
*-----------------------------------------------------------------------

         if( itaty(m) .ne. 0 ) then

               call trnsuu(upart,vpart,wpart,
     &                     ucc,vcc,wcc,itmtr(m,4))

                  cst = 1.d0
                  pab = sqrt( ucc**2 + vcc**2 + wcc**2 )
                  if( pab .gt. 0.0d0 ) cst = wcc / pab

               if( itaty(m) .gt. 0 ) then

                  if( cst .lt. ab(1) ) goto 100
                  if( cst .gt. ab(na+1) ) goto 100

               else

                  if( cst .lt. cos( ab(na+1) / 180.d0 * pi ) ) goto 100
                  if( cst .gt. cos( ab(1) / 180.d0 * pi ) ) goto 100

               end if

            do i = 1, na

               if( itaty(m) .gt. 0 ) then

                  if( cst .ge. ab(i) .and.
     &                cst .le. ab(i+1) ) goto 50

               else

                  if( cst .ge. cos( ab(i+1) / 180.d0 * pi ) .and.
     &                cst .le. cos( ab(i) / 180.d0 * pi ) ) goto 50

               end if

            end do

   50          ia1 = i

         else

               ia1 = 1

         end if

*-----------------------------------------------------------------------
*        tally
*-----------------------------------------------------------------------
               do ip = 1, ipn

                  trEVENT(ips(ip),ie1,it1,ia1,icf(ix1,iy1,iz1)) =
     &            trEVENT(ips(ip),ie1,it1,ia1,icf(ix1,iy1,iz1)) + tlw

               end do

               if (istdev .eq. 2) then

                 call setitrminmax(itrmin,itrmax,2,7,
     &                             (/ie1,it1,ia1,ix1,iy1,iz1/))

               endif
*-----------------------------------------------------------------------
*        dump data on file
*-----------------------------------------------------------------------

         if( itmdp(m,0) .ne. 0 ) then
               dmpd(1)  = dble( kpart )
               dmpd(2)  = xpart
               dmpd(3)  = ypart
               dmpd(4)  = zpart
               dmpd(5)  = upart
               dmpd(6)  = vpart
               dmpd(7)  = wpart
               dmpd(8)  = epart
               dmpd(9)  = tlw
               dmpd(10) = tpart
               dmpd(11) = ncntt(1)
               dmpd(12) = ncntt(2)
               dmpd(13) = ncntt(3)
               dmpd(14) = sxpat
               dmpd(15) = sypat
               dmpd(16) = szpat
               dmpd(17) = nname
               dmpd(18) = nocas
               dmpd(19) = nobch
               dmpd(20) = no
               dmpd(21) = nzpat

               if( ipart .ge. 15 )
     &         dmpd(8) = dmpd(8) / ibryf(ipart,kpart)

            if(idmpomp.ne.0)then       !FURUTA20150427
             call writeompfile(m,dmpd) !FURUTA20150427
            else                       !FURUTA20150427

             if( itmdp(m,0) .gt. 0 ) then

               write(itmdf(m))
     &         ( dmpd(itmdp(m,k)), k = 1, abs( itmdp(m,0) ) )

              else

               write(itmdf(m),'(30(1p1d24.15))')
     &         ( dmpd(itmdp(m,k)), k = 1, abs( itmdp(m,0) ) )

              end if

            endif !FURUTA20150427

         end if

*-----------------------------------------------------------------------

  100    continue

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine ppdctxyz(m,np,nl,lt,
     &                    nx,ny,nz,ne,nt,na,xm,ym,zm,eb,tb,ab,tr,
     &                    igsh,idasa)
*                                                                      *
*       output xyz scoring mesh product tally                          *
*       last modified by K.Niita on 2005/11/24                         *
*                                                                      *
************************************************************************
      use sangelmod, only: itsans
      use partmod ! frtati 2021/10/05

*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'

      parameter ( pi  = 3.141592653589793d0 )

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
      common /tall07/ itety(itlmax), itenm(itlmax), iterg(itlmax),
     &                rtemi(itlmax), rtema(itlmax), rtedl(itlmax)

      common /tall20/ ittty(itlmax), ittnm(itlmax), ittrg(itlmax),
     &                rttmi(itlmax), rttma(itlmax), rttdl(itlmax)


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
      common /tall48/ itmdp(itlmax,0:30), itmdf(itlmax)
      common /tall49/ itglt(itlmax)

      common /tall23/ itaty(itlmax), itanm(itlmax), itarg(itlmax),
     &                rtami(itlmax), rtama(itlmax), rtadl(itlmax)

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

      common /tall50/ itlmt(itlmax), ite2l(itlmax)
*-----------------------------------------------------------------------

      dimension   lt(nl)
      dimension   xm(nx+1)
      dimension   ym(ny+1)
      dimension   zm(nz+1)
      dimension   eb(ne+1)
      dimension   tb(nt+1)
      dimension   ab(na+1)
      dimension   ew(ne)
      dimension   tw(nt)
      dimension   aw(na)
      dimension   tr(np,ne,nt,na,nx*ny*nz,2)

      integer,allocatable :: ixyz(:)

*-----------------------------------------------------------------------

      dimension tott(np,2) ! frtati 2021/10/05 6 -> np

      character hsunit(36)*32

      data hsunit( 1) / '[1/source]                      '/
      data hsunit( 2) / '[1/cm^3/source]                 '/
      data hsunit( 3) / '[1/MeV/source]                  '/
      data hsunit( 4) / '[1/cm^3/MeV/source]             '/
      data hsunit( 5) / '[1/Lethargy/source]             '/
      data hsunit( 6) / '[1/cm^3/Lethargy/source]        '/
      data hsunit(11) / '[1/nsec/source]                 '/
      data hsunit(12) / '[1/cm^3/nsec/source]            '/
      data hsunit(13) / '[1/MeV/nsec/source]             '/
      data hsunit(14) / '[1/cm^3/MeV/nsec/source]        '/
      data hsunit(15) / '[1/Lethargy/nsec/source]        '/
      data hsunit(16) / '[1/cm^3/Lethargy/nsec/source]   '/
      data hsunit(21) / '[1/sr/source]                   '/
      data hsunit(22) / '[1/cm^3/sr/source]              '/
      data hsunit(23) / '[1/MeV/sr/source]               '/
      data hsunit(24) / '[1/cm^3/MeV/sr/source]          '/
      data hsunit(25) / '[1/Lethargy/sr/source]          '/
      data hsunit(26) / '[1/cm^3/Lethargy/sr/source]     '/
      data hsunit(31) / '[1/nsec/sr/source]              '/
      data hsunit(32) / '[1/cm^3/nsec/sr/source]         '/
      data hsunit(33) / '[1/MeV/nsec/sr/source]          '/
      data hsunit(34) / '[1/cm^3/MeV/nsec/sr/source]     '/
      data hsunit(35) / '[1/Lethargy/nsec/sr/source]     '/
      data hsunit(36) / '[1/cm^3/Lethargy/nsec/sr/source]'/

      integer           iMeVperu
      common /cMeVperu/ iMeVperu

      character cha*1
      data cha /"'"/

      character dc2*4
      character chp(itmxpt)*11  ! kitamura22/03/31
      character chq(itmxpt)*9   ! kitamura22/03/31

      character aname*3


*-----------------------------------------------------------------------

      character rpa*1
      data rpa /'}'/
      character yen*1

! sumover
      real(8),allocatable :: vl_x(:,:),vl_y(:,:),vl_z(:,:)

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

      include 'samepage_include/samepage000.inc'

*-----------------------------------------------------------------------

      icf(ix,iy,iz) = ix + ( iy - 1 ) * nx + ( iz - 1 ) * nx * ny

*-----------------------------------------------------------------------
*        set mesh volume
*-----------------------------------------------------------------------

         vl(ix,iy,iz) = vls(nl,lt,itmcn(m),itvm(m),itmtr(m,4),
     &                      xm(ix),xm(ix+1),
     &                      ym(iy),ym(iy+1),
     &                      zm(iz),zm(iz+1))

*-----------------------------------------------------------------------

      include 'samepage_include/samepage001.inc'

! sumover
      allocate (vl_x(ny,nz),vl_y(nx,nz),vl_z(nx,ny))
      if( itunt(m) .eq.  2 .or. itunt(m) .eq.  4 .or.
     &    itunt(m) .eq.  6 .or. itunt(m) .eq. 12 .or.
     &    itunt(m) .eq. 14 .or. itunt(m) .eq. 16 .or.
     &    itunt(m) .eq. 22 .or. itunt(m) .eq. 24 .or.
     &    itunt(m) .eq. 26 .or. itunt(m) .eq. 32 .or.
     &    itunt(m) .eq. 34 .or. itunt(m) .eq. 36 ) then
        vl_x(:,:) = 0.0d0
        vl_y(:,:) = 0.0d0
        vl_z(:,:) = 0.0d0
        do iz = 1, nz
          do iy = 1, ny
            do ix = 1, nx
              vl_x(iy,iz) = vl_x(iy,iz) + vl(ix,iy,iz)
              vl_y(ix,iz) = vl_y(ix,iz) + vl(ix,iy,iz)
              vl_z(ix,iy) = vl_z(ix,iy) + vl(ix,iy,iz)
            enddo
          enddo
        enddo
      else
        vl_x(:,:) = 1.0d0
        vl_y(:,:) = 1.0d0
        vl_z(:,:) = 1.0d0
      endif

      yen  = char(92)
      allocate( ixyz(max((nx+1)*(ny+1),(ny+1)*(nz+1),(nx+1)*(nz+1))) )

*-----------------------------------------------------------------------

      if ( iMeVperu.eq. 1 ) then
         hsunit( 3) = '[1/(MeV/n)/source]              '
         hsunit( 4) = '[1/cm^3/(MeV/n)/source]         '
         hsunit(13) = '[1/(MeV/n)/nsec/source]         '
         hsunit(14) = '[1/cm^3/(MeV/n)/nsec/source]    '
         hsunit(23) = '[1/(MeV/n)/sr/source]           '
         hsunit(24) = '[1/cm^3/(MeV/n)/sr/source]      '
         hsunit(33) = '[1/(MeV/n)/nsec/sr/source]      '
         hsunit(34) = '[1/cm^3/(MeV/n)/nsec/sr/source] '
      end if

      if ( ite2l(m) .eq. 1 ) then   ! convert to energy to LET
         hsunit( 3) = '[1/(keV/um)/source]             '
         hsunit( 4) = '[1/cm^3/(keV/um)/source]        '
         hsunit(13) = '[1/(keV/um)/nsec/source]        '
         hsunit(14) = '[1/cm^3/(keV/um)/nsec/source]   '
         hsunit(23) = '[1/(keV/um)/sr/source]          '
         hsunit(24) = '[1/cm^3/(keV/um)/sr/source]     '
         hsunit(33) = '[1/(keV/um)/nsec/sr/source]     '
         hsunit(34) = '[1/cm^3/(keV/um)/nsec/sr/source]'
      end if
*-----------------------------------------------------------------------
      np_mxang = np
      if ( itmxang(m).ne.0 .and. iabs(itmxang(m)).lt.np ) then
        np_mxang = iabs(itmxang(m))
      end if

*-----------------------------------------------------------------------

            if( itaty(m) .gt. 0 ) then
               aname = 'cos'
            else
               aname = 'the'
            end if

*-----------------------------------------------------------------------
*     for igshow
*-----------------------------------------------------------------------

         if( igsh .eq. 0 ) then

            npg = np
            neg = ne
            ntg = nt
            nag = na

         else

            npg = 1
            neg = 1
            ntg = 1
            nag = 1

         end if

*-----------------------------------------------------------------------
*        particle name
*-----------------------------------------------------------------------

               call wtpname(m,np,chp,chq)

*-----------------------------------------------------------------------
*           itunt(m) = 3 13 23 33 4 14 24 34: /MeV
*-----------------------------------------------------------------------

            if( itunt(m) .eq.  3 .or. itunt(m) .eq.  4 .or.
     &          itunt(m) .eq. 13 .or. itunt(m) .eq. 14 .or.
     &          itunt(m) .eq. 23 .or. itunt(m) .eq. 24 .or.
     &          itunt(m) .eq. 33 .or. itunt(m) .eq. 34 ) then

               do i = 1, ne

                  ew(i) = eb(i+1) - eb(i)

               end do
               ew_sum = eb(ne+1) - eb(1)

            else if( itunt(m) .eq.  5 .or. itunt(m) .eq.  6 .or.
     &               itunt(m) .eq. 15 .or. itunt(m) .eq. 16 .or.
     &               itunt(m) .eq. 25 .or. itunt(m) .eq. 26 .or.
     &               itunt(m) .eq. 35 .or. itunt(m) .eq. 36 ) then

               do i = 1, ne

                  ew(i) = log( eb(i+1) / eb(i) )

               end do
               ew_sum = log( eb(ne+1) / eb(1))

            else

               do i = 1, ne

                  ew(i) = 1.d+0

               end do
               ew_sum = 1.0d0

            end if

*-----------------------------------------------------------------------
*           itunt(m) = 11, 12, 13, 14, 31, 32, 33, 34 : /nsec
*-----------------------------------------------------------------------

            if( ( itunt(m) .ge. 11 .and. itunt(m) .le. 14 ) .or.
     &          ( itunt(m) .ge. 31 .and. itunt(m) .le. 34 ) ) then

               do i = 1, nt

                  tw(i) = tb(i+1) - tb(i)

               end do
               tw_sum = tb(nt+1) - tb(1)

            else

               do i = 1, nt

                  tw(i) = 1.d+0

               end do
               tw_sum = 1.0d0

            end if

*-----------------------------------------------------------------------
*           itunt(m) = 20+, 30+ : /SR
*-----------------------------------------------------------------------

            if( itunt(m) .lt. 20 ) then

               do i = 1, na

                  aw(i) = 1.d+0

               end do
               aw_sum = 1.0d0

            else

               aw_sum = 0.0d0
               do i = 1, na

                  if( itaty(m) .gt. 0 ) then

                     aw(i) = ( ab(i+1) - ab(i) ) * 2.0 * pi

                  else

                     aw(i) = ( cos( ab(i) / 180.d0 * pi )
     &                       - cos( ab(i+1) / 180.d0 * pi ) ) * 2.0 * pi

                  end if
                  aw_sum = aw_sum + aw(i)

               end do

            end if

*-----------------------------------------------------------------------
*           itunt(m) = 2, 4, 6, +10,20,30  : /cm^3
*-----------------------------------------------------------------------
*        relative error and  unit conversion
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

               do 101 iz = 1, nz
               do 101 iy = 1, ny
               do 101 ix = 1, nx
               do 101 ie = 1, ne
               do 101 it = 1, nt
               do 101 ia = 1, na
               do 101 ip = 1, np

                  if( tr(ip,ie,it,ia,icf(ix,iy,iz),1) .gt. 0.d0 ) then

                     if( itunt(m) .eq.  2 .or. itunt(m) .eq.  4 .or.
     &                   itunt(m) .eq.  6 .or. itunt(m) .eq. 12 .or.
     &                   itunt(m) .eq. 14 .or. itunt(m) .eq. 16 .or.
     &                   itunt(m) .eq. 22 .or. itunt(m) .eq. 24 .or.
     &                   itunt(m) .eq. 26 .or. itunt(m) .eq. 32 .or.
     &                   itunt(m) .eq. 34 .or. itunt(m) .eq. 36 ) then

                        vm = vl(ix,iy,iz)

                     else

                        vm = 1.d+0

                     end if

                     fmaxfc = tr(ip,ie,it,ia,icf(ix,iy,iz),1)
     &                                      / vm / ew(ie) / aw(ia)

                     if( facmax(m) .lt. fmaxfc) facmax(m) = fmaxfc

                  end if

  101          continue

            end if

            do 100 iz = 1, nz
            do 100 iy = 1, ny
            do 100 ix = 1, nx
            do 100 ie = 1, ne
            do 100 it = 1, nt
            do 100 ia = 1, na
            do 100 ip = 1, np

                  if( itunt(m) .eq.  2 .or. itunt(m) .eq.  4 .or.
     &                itunt(m) .eq.  6 .or. itunt(m) .eq. 12 .or.
     &                itunt(m) .eq. 14 .or. itunt(m) .eq. 16 .or.
     &                itunt(m) .eq. 22 .or. itunt(m) .eq. 24 .or.
     &                itunt(m) .eq. 26 .or. itunt(m) .eq. 32 .or.
     &                itunt(m) .eq. 34 .or. itunt(m) .eq. 36 ) then

                     vm = vl(ix,iy,iz)

                  else

                     vm = 1.d+0

                  end if

               if( tr(ip,ie,it,ia,icf(ix,iy,iz),1) .gt. 0.d0 ) then


                  cc = abs(rtfac(m)/facmax(m))
     &                          / vm / ew(ie) / aw(ia)

                  call calc_stdev(m,Xa,sigx,
     &                            tr(ip,ie,it,ia,icf(ix,iy,iz),1),
     &                            tr(ip,ie,it,ia,icf(ix,iy,iz),2),
     &                            cc)

                  tr(ip,ie,it,ia,icf(ix,iy,iz),1) = Xa
                  tr(ip,ie,it,ia,icf(ix,iy,iz),2) = sigx
                  if( sigx .gt. stdm ) stdm = sigx

                  if( tr(ip,ie,it,ia,icf(ix,iy,iz),1) .gt. cmax )
     &                     cmax = tr(ip,ie,it,ia,icf(ix,iy,iz),1)

                  if( tr(ip,ie,it,ia,icf(ix,iy,iz),1) .lt. cmin )
     &                     cmin = tr(ip,ie,it,ia,icf(ix,iy,iz),1)

               else

                  isdz = 1
                  tr(ip,ie,it,ia,icf(ix,iy,iz),2) = 0.d+0

               end if

! sumover
               fact_in = abs(rtfac(m)/facmax(m))
               call ppdctxyz_sumover_stdev(0,m,ip,ie,it,ia,ix,iy,iz,
     &                fact_in,ew(ie),1.0d0,aw(ia),vm,
     &                ew_sum,1.0d0,aw_sum,
     &                vl_x(iy,iz),vl_y(ix,iz),vl_z(ix,iy))

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

!OBINATA(2012.9.11): output *.err
        noe = 1
        if ( any( itaxs(m,iax) .eq. (/ 7, 8, 9, 10 /) )
     &      .and. ittwo(m) .ne. 4 .and. icntl .ne. 8 )
     &    noe = 2

        do 900 ioe = 1, noe

         if( igsh .ne. 0 .and.
     &     ( itaxs(m,iax) .lt. 7 .or.
     &       ittwo(m) .eq. 4 .or. ittwo(m) .eq. 5 ) ) goto 900

         if( itmdp(m,0) .eq. 0 ) then

            if( itall .eq. 2 .and. nobch .lt. maxbch .and.
     &          igsh .eq. 0 ) then

               write(fnume,'(i3.3)') nobch
!OBINATA(2012.9.11): output *.err
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
     &               nobch,maxbch,npe)
              else
                call mk_2dnerfn(ctfln(m,iax),fname,itfll(m,iax),fnume)
              end if

            else

!OBINATA(2012.9.11): output *.err
              if ( ioe .eq. 1 ) then
                fname = ctfln(m,iax)
              else
                call mk_2derrfn(ctfln(m,iax),fname,itfll(m,iax))
              end if

            end if

         else

            if( (itall .eq. 2 .and. nobch .lt. maxbch .and. igsh .eq. 0)
     &           .or. ( itall .eq. 4 .and. igsh .eq. 0 ) ) then

               call mk_2dnumfn(ctfln(m,iax),fname,itfll(m,iax),
     &              nobch,maxbch,npe)

            else

               fname = ctfln(m,iax)(1:itfll(m,iax))

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

               call tproech(iot,m,iax,1)


      include 'samepage_include/samepage002_petaxyz.inc'
      include 'samepage_include/samepagechp_petaxyz.inc'
      include 'samepage_include/samepageseti.inc'

*-----------------------------------------------------------------------
*        energy axis  ( LET axis )
*-----------------------------------------------------------------------

         if( itaxs(m,iax) .eq. 1 .or.
     &       itaxs(m,iax) .eq. 15 ) then

               inum = 0

            do ixi = 1, nx, nxstepi
            do iyi = 1, ny, nystepi
            do izi = 1, nz, nzstepi
            do iti = 1, nt, ntstepi
            do iai = 1, na, nastepi
            do ipi = 1, np, npstepi
               ix = ixi
               iy = iyi
               iz = izi
               it = iti
               ia = iai
               ip = ip

               inum = inum + 1

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

               if(iloopmode .ge. 1) then
                  write(iot,'("#   ",a)') chp(ipi)(3:10)
               endif

                  write(iot,'("#   no. =",i3,3x,
     &            "ix  =",i3,3x,
     &            "iy  =",i3,3x,
     &            "iz  =",i3/
     &            "#   x = (",1p1e13.4,"  -",1p1e13.4,"  )"/
     &            "#   y = (",1p1e13.4,"  -",1p1e13.4,"  )"/
     &            "#   z = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        inum, ix, iy, iz,
     &                        xm(ix), xm(ix+nxstepi),
     &                        ym(iy), ym(iy+nystepi),
     &                        zm(iz), zm(iz+nzstepi)
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  it =",i3/
     &            "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                  it, tb(it), tb(it+ntstepi)
               end if

               if( itaty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  ia =",i3,/
     &            "# ",a3," = (",
     &                           1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     ia, aname, ab(ia), ab(ia+nastepi)
               end if

               if( itaxl(m) .eq. 0 ) then

                if ( itaxs(m,iax) .eq. 1 ) then   ! ccse 2022/09/30

                if ( iMeVperu.eq. 1 ) then
                  write(iot,'(/"x: Energy [MeV/n]")')
                else
                  write(iot,'(/"x: Energy [MeV]")')
                end if

                else if ( itaxs(m,iax) .eq. 15 ) then
                 write(iot,'(/"x: LET [keV/um]")')
                end if

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: Number ",a32)') hsunit(itunt(m))

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

               if( itety(m) .eq.  3 .or. itety(m) .eq.  5 .or.
     &             itunt(m) .eq.  5 .or. itunt(m) .eq.  6 .or.
     &             itunt(m) .eq. 15 .or. itunt(m) .eq. 16 .or.
     &             itunt(m) .eq. 25 .or. itunt(m) .eq. 26 .or.
     &             itunt(m) .eq. 35 .or. itunt(m) .eq. 36 ) then

                  write(iot,'( "p: xlog ylog afac(0.8) form(0.9)")')

               else

                  write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')

               end if

               if( itanl(m) .gt. 0 ) then

                  write(iot,'( "p: ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               end if

               if( itsans(m) .gt. 0 ) then

                  call write_sangel(iot,m,0)

               end if

       include 'samepage_include/petaxyz_e.inc'

             if(iloopmode .eq. 0) then
               if( ittty(m) .ne. 0 .and. itaty(m) .ne. 0 ) then
                  write(iot,'(/a1,"no. =",i3,
     &            ",  (ix,iy,iz) = (",i3,",",i3,",",i3,")",
     &                     ",   t =",i3,"   ang =",i3,a1)')
     &                     cha, inum, ix, iy, iz, it, ia, cha
               else if( ittty(m) .eq. 0 .and. itaty(m) .ne. 0 ) then
                  write(iot,'(/a1,"no. =",i3,
     &            ",  (ix,iy,iz) = (",i3,",",i3,",",i3,")",
     &                     ",  ang =",i3,a1)')
     &                     cha, inum, ix, iy, iz, ia, cha
               else if( ittty(m) .ne. 0 .and. itaty(m) .eq. 0 ) then
                  write(iot,'(/a1,"no. =",i3,
     &            ",  (ix,iy,iz) = (",i3,",",i3,",",i3,")",
     &                     ",   t =",i3,a1)')
     &                     cha, inum, ix, iy, iz, it, cha
               else
                  write(iot,'(/a1,"no. =",i3,
     &            ",  (ix,iy,iz) = (",i3,",",i3,",",i3,")",
     &            a1)')
     &            cha, inum, ix, iy, iz, cha
               end if
             else
               write(changelsub(1),'(a1,"no. =",i3)') cha,inum
               write(changelsub(2),'(",  part = ",a)') chp(ipi)(3:10)
               write(changelsub(3),'(",  e =",i3)') ie
               write(changelsub(4),'(",  t =",i3)') it
               write(changelsub(5),'(",  a =",i3)') ia
               write(changelsub(6),'(",  x =",i3)') ix
               write(changelsub(7),'(",  y =",i3)') iy
               write(changelsub(8),'(",  z =",i3)') iz
               write(changelsub(9),'(a1)') cha

               if( ne .eq. 0 ) then
                 changelsub(3) = " "
               end if
               if( ittty(m) .eq. 0 ) then
                 changelsub(4) = " "
               end if
               if( itaty(m) .eq. 0 ) then
                 changelsub(5) = " "
               end if

               changelsub(3) = " "

               if(iloopmode .eq. 1) then
                 changelsub(3) = " "
               end if
               if(iloopmode .eq. 10) then
                 changelsub(4) = " "
               end if
               if(iloopmode .eq. 8 .or. iloopmode .eq. 9) then
                 changelsub(5) = " "
               end if
               if(iloopmode .eq. 3) then
                 changelsub(6) = " "
               end if
               if(iloopmode .eq. 4) then
                 changelsub(7) = " "
               end if
               if(iloopmode .eq. 5) then
                 changelsub(8) = " "
               end if
                angeltitle = trim(changelsub(1))//trim(changelsub(2))//
     &                       trim(changelsub(3))//trim(changelsub(4))//
     &                       trim(changelsub(5))//trim(changelsub(6))//
     &                       trim(changelsub(7))//trim(changelsub(8))//
     &                       trim(changelsub(9))
                write(iot,'(/a)') trim(angeltitle)
          end if


               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

                  write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "   vol  &=&",1pe13.4," [cm^3]"/
     &                     "  xmin  &=&",1pe13.4," [cm]"/
     &                     "  xmax  &=&",1pe13.4," [cm]"/
     &                     "  ymin  &=&",1pe13.4," [cm]"/
     &                     "  ymax  &=&",1pe13.4," [cm]"/
     &                     "  zmin  &=&",1pe13.4," [cm]"/
     &                     "  zmax  &=&",1pe13.4," [cm]")')
     &                     yen, vm, xm(ix), xm(ix+nxstepi),
     &                  ym(iy), ym(iy+nystepi), zm(iz), zm(iz+nzstepi)
               if(iloopmode .ge. 1) then
                     write(iot,'("  part  &=& ",a)') chp(ipi)(3:10)
               endif

               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &                     "  tmin  &=&",1pe13.4," [nsec]"/
     &                     "  tmax  &=&",1pe13.4," [nsec]")')
     &                     tb(it), tb(it+ntstepi)
               end if
               if( itaty(m) .ne. 0 ) then
                  write(iot,'(
     &                     "  amin  &=&",1pe13.4,/
     &                     "  amax  &=&",1pe13.4)')
     &                     ab(ia), ab(ia+nastepi)
               end if

                  write(iot,'("e:")')

            end do

            end do
            end do
            end do
            end do
            end do

*-----------------------------------------------------------------------
*        time axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 11 ) then

               inum = 0

            do ixi = 1, nx, nxstepi
            do iyi = 1, ny, nystepi
            do izi = 1, nz, nzstepi
            do iei = 1, ne, nestepi
            do iai = 1, na, nastepi
            do ipi = 1, np, npstepi
               ix = ixi
               iy = iyi
               iz = izi
               ie = iei
               ia = iai
               ip = ipi

               inum = inum + 1

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

               if(iloopmode .ge. 1) then
                  write(iot,'("#   ",a)') chp(ipi)(3:10)
               endif

               write(iot,'("#   no. =",i3,3x,
     &         "ix  =",i3,3x,
     &         "iy  =",i3,3x,
     &         "iz  =",i3,3x,
     &         "ie  =",i3,3x,/
     &         "#   x = (",1p1e13.4,"  -",1p1e13.4,"  )"/
     &         "#   y = (",1p1e13.4,"  -",1p1e13.4,"  )"/
     &         "#   z = (",1p1e13.4,"  -",1p1e13.4,"  )"/
     &         "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     inum, ix, iy, iz, ie,
     &                     xm(ix), xm(ix+nxstepi),
     &                     ym(iy), ym(iy+nystepi),
     &                     zm(iz), zm(iz+nzstepi),
     &                     eb(ie), eb(ie+nestepi)

               if( itaty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  ia =",i3,/
     &            "# ",a3," = (",
     &                           1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     ia, aname, ab(ia), ab(ia+nastepi)
               end if

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: Time [nsec]")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: Number ",a32)') hsunit(itunt(m))

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

               if( ittty(m) .eq. 3 .or. ittty(m) .eq. 5 ) then

                  write(iot,'( "p: xlog ylog afac(0.8) form(0.9)")')

               else

                  write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')

               end if

               if( itanl(m) .gt. 0 ) then

                  write(iot,'( "p: ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               end if

               if( itsans(m) .gt. 0 ) then

                  call write_sangel(iot,m,0)

               end if

       include 'samepage_include/petaxyz_t.inc'

             if(iloopmode .eq. 0) then
               if( itaty(m) .ne. 0 ) then
                  write(iot,'(/a1,"no. =",i3,
     &            ",  (ix,iy,iz) = (",i3,",",i3,",",i3,")",
     &            ",    ie =",i3,
     &            ",  ang =",i3,a1)')
     &                     cha, inum, ix, iy, iz, ie, ia, cha
               else
                  write(iot,'(/a1,"no. =",i3,
     &            ",    (ix,iy,iz) = (",i3,",",i3,",",i3,")",
     &            ",    ie =",i3,a1)')
     &            cha, inum, ix, iy, iz, ie, cha
               end if
             else
               write(changelsub(1),'(a1,"no. =",i3)') cha,inum
               write(changelsub(2),'(",  part = ",a)') chp(ipi)(3:10)
               write(changelsub(3),'(",  e =",i3)') ie
               write(changelsub(4),'(",  t =",i3)') it
               write(changelsub(5),'(",  a =",i3)') ia
               write(changelsub(6),'(",  x =",i3)') ix
               write(changelsub(7),'(",  y =",i3)') iy
               write(changelsub(8),'(",  z =",i3)') iz
               write(changelsub(9),'(a1)') cha

               if( ne .eq. 0 ) then
                 changelsub(3) = " "
               end if
               if( ittty(m) .eq. 0 ) then
                 changelsub(4) = " "
               end if
               if( itaty(m) .eq. 0 ) then
                 changelsub(5) = " "
               end if

               changelsub(4) = " "

               if(iloopmode .eq. 1) then
                 changelsub(3) = " "
               end if
               if(iloopmode .eq. 10) then
                 changelsub(4) = " "
               end if
               if(iloopmode .eq. 8 .or. iloopmode .eq. 9) then
                 changelsub(5) = " "
               end if
               if(iloopmode .eq. 3) then
                 changelsub(6) = " "
               end if
               if(iloopmode .eq. 4) then
                 changelsub(7) = " "
               end if
               if(iloopmode .eq. 5) then
                 changelsub(8) = " "
               end if
                angeltitle = trim(changelsub(1))//trim(changelsub(2))//
     &                       trim(changelsub(3))//trim(changelsub(4))//
     &                       trim(changelsub(5))//trim(changelsub(6))//
     &                       trim(changelsub(7))//trim(changelsub(8))//
     &                       trim(changelsub(9))
                write(iot,'(/a)') trim(angeltitle)
          end if


               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

                if( ite2l(m).eq.0 ) then  ! ccse 2022/09/30
                  if( iMeVperu.eq.0 ) then
                     write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "   vol  &=&",1pe13.4," [cm^3]"/
     &                     "  xmin  &=&",1pe13.4," [cm]"/
     &                     "  xmax  &=&",1pe13.4," [cm]"/
     &                     "  ymin  &=&",1pe13.4," [cm]"/
     &                     "  ymax  &=&",1pe13.4," [cm]"/
     &                     "  zmin  &=&",1pe13.4," [cm]"/
     &                     "  zmax  &=&",1pe13.4," [cm]"/
     &                     "  emin  &=&",1pe13.4," [MeV]"/
     &                     "  emax  &=&",1pe13.4," [MeV]")')
     &                     yen, vm, xm(ix), xm(ix+nxstepi),
     &                  ym(iy), ym(iy+nystepi), zm(iz), zm(iz+nzstepi),
     &                     eb(ie), eb(ie+nestepi)
                  else
                     write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "   vol  &=&",1pe13.4," [cm^3]"/
     &                     "  xmin  &=&",1pe13.4," [cm]"/
     &                     "  xmax  &=&",1pe13.4," [cm]"/
     &                     "  ymin  &=&",1pe13.4," [cm]"/
     &                     "  ymax  &=&",1pe13.4," [cm]"/
     &                     "  zmin  &=&",1pe13.4," [cm]"/
     &                     "  zmax  &=&",1pe13.4," [cm]"/
     &                     "  emin  &=&",1pe13.4," [MeV/n]"/
     &                     "  emax  &=&",1pe13.4," [MeV/n]")')
     &                     yen, vm, xm(ix), xm(ix+nxstepi),
     &                  ym(iy), ym(iy+nystepi), zm(iz), zm(iz+nzstepi),
     &                     eb(ie), eb(ie+nestepi)
                  end if
                else if( ite2l(m).eq.1 ) then
                     write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "   vol  &=&",1pe13.4," [cm^3]"/
     &                     "  xmin  &=&",1pe13.4," [cm]"/
     &                     "  xmax  &=&",1pe13.4," [cm]"/
     &                     "  ymin  &=&",1pe13.4," [cm]"/
     &                     "  ymax  &=&",1pe13.4," [cm]"/
     &                     "  zmin  &=&",1pe13.4," [cm]"/
     &                     "  zmax  &=&",1pe13.4," [cm]"/
     &                     "  emin  &=&",1pe13.4," [keV/um]"/
     &                     "  emax  &=&",1pe13.4," [keV/um]")')
     &                     yen, vm, xm(ix), xm(ix+nxstepi),
     &                  ym(iy), ym(iy+nystepi), zm(iz), zm(iz+nzstepi),
     &                     eb(ie), eb(ie+nestepi)
                end if
               if(iloopmode .ge. 1) then
                     write(iot,'("  part  &=& ",a)') chp(ipi)(3:10)
               endif

               if( itaty(m) .ne. 0 ) then
                  write(iot,'(
     &                     "  amin  &=&",1pe13.4,/
     &                     "  amax  &=&",1pe13.4)')
     &                     ab(ia), ab(ia+nastepi)
               end if

                  write(iot,'("e:")')

            end do

            end do
            end do
            end do
            end do
            end do

*-----------------------------------------------------------------------
*        angle axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 12 .or.
     &            itaxs(m,iax) .eq. 13 ) then

               inum = 0

            do ixi = 1, nx, nxstepi
            do iyi = 1, ny, nystepi
            do izi = 1, nz, nzstepi
            do iei = 1, ne, nestepi
            do iti = 1, nt, ntstepi
            do ipi = 1, np, npstepi
               ix = ixi
               iy = iyi
               iz = izi
               ie = iei
               it = iti
               ip = ipi

               inum = inum + 1

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

               if(iloopmode .ge. 1) then
                  write(iot,'("#   ",a)') chp(ipi)(3:10)
               endif

               write(iot,'("#   no. =",i3,3x,
     &         "ix  =",i3,3x,
     &         "iy  =",i3,3x,
     &         "iz  =",i3,3x,
     &         "ie  =",i3,3x,/
     &         "#   x = (",1p1e13.4,"  -",1p1e13.4,"  )"/
     &         "#   y = (",1p1e13.4,"  -",1p1e13.4,"  )"/
     &         "#   z = (",1p1e13.4,"  -",1p1e13.4,"  )"/
     &         "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     inum, ix, iy, iz, ie,
     &                     xm(ix), xm(ix+nxstepi),
     &                     ym(iy), ym(iy+nystepi),
     &                     zm(iz), zm(iz+nzstepi),
     &                     eb(ie), eb(ie+nestepi)

               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  it =",i3/
     &            "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                  it, tb(it), tb(it+ntstepi)
               end if

               if( itaxl(m) .eq. 0 ) then

                  if( itaxs(m,iax) .eq. 12 ) then
                     write(iot,'(/"x: cos(",a1,"theta)")') yen
                  else
                     write(iot,'(/"x: ",a1,"theta  [deg]")') yen
                  end if

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: Number ",a32)') hsunit(itunt(m))

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

       include 'samepage_include/petaxyz_a.inc'

             if(iloopmode .eq. 0) then
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(/a1,"no. =",i3,
     &            ",  (ix,iy,iz) = (",i3,",",i3,",",i3,")",
     &            ",    e =",i3,
     &            ",    t =",i3,a1)')
     &                     cha, inum, ix, iy, iz, ie, it, cha
               else
                  write(iot,'(/a1,"no. =",i3,
     &            ",    (ix,iy,iz) = (",i3,",",i3,",",i3,")",
     &            ",   e =",i3,a1)')
     &            cha, inum, ix, iy, iz, ie, cha
               end if
             else
               write(changelsub(1),'(a1,"no. =",i3)') cha,inum
               write(changelsub(2),'(",  part = ",a)') chp(ipi)(3:10)
               write(changelsub(3),'(",  e =",i3)') ie
               write(changelsub(4),'(",  t =",i3)') it
               write(changelsub(5),'(",  a =",i3)') ia
               write(changelsub(6),'(",  x =",i3)') ix
               write(changelsub(7),'(",  y =",i3)') iy
               write(changelsub(8),'(",  z =",i3)') iz
               write(changelsub(9),'(a1)') cha

               if( ne .eq. 0 ) then
                 changelsub(3) = " "
               end if
               if( ittty(m) .eq. 0 ) then
                 changelsub(4) = " "
               end if
               if( itaty(m) .eq. 0 ) then
                 changelsub(5) = " "
               end if

               changelsub(5) = " "

               if(iloopmode .eq. 1) then
                 changelsub(3) = " "
               end if
               if(iloopmode .eq. 10) then
                 changelsub(4) = " "
               end if
               if(iloopmode .eq. 8 .or. iloopmode .eq. 9) then
                 changelsub(5) = " "
               end if
               if(iloopmode .eq. 3) then
                 changelsub(6) = " "
               end if
               if(iloopmode .eq. 4) then
                 changelsub(7) = " "
               end if
               if(iloopmode .eq. 5) then
                 changelsub(8) = " "
               end if
                angeltitle = trim(changelsub(1))//trim(changelsub(2))//
     &                       trim(changelsub(3))//trim(changelsub(4))//
     &                       trim(changelsub(5))//trim(changelsub(6))//
     &                       trim(changelsub(7))//trim(changelsub(8))//
     &                       trim(changelsub(9))
                write(iot,'(/a)') trim(angeltitle)
          end if


               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn



                if( ite2l(m).eq.0 ) then  ! ccse 2022/09/30
                  if( iMeVperu.eq.0 ) then
                     write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "   vol  &=&",1pe13.4," [cm^3]"/
     &                     "  xmin  &=&",1pe13.4," [cm]"/
     &                     "  xmax  &=&",1pe13.4," [cm]"/
     &                     "  ymin  &=&",1pe13.4," [cm]"/
     &                     "  ymax  &=&",1pe13.4," [cm]"/
     &                     "  zmin  &=&",1pe13.4," [cm]"/
     &                     "  zmax  &=&",1pe13.4," [cm]"/
     &                     "  emin  &=&",1pe13.4," [MeV]"/
     &                     "  emax  &=&",1pe13.4," [MeV]")')
     &                     yen, vm, xm(ix), xm(ix+nxstepi),
     &                  ym(iy), ym(iy+nystepi), zm(iz), zm(iz+nzstepi),
     &                     eb(ie), eb(ie+nestepi)
                  else
                     write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "   vol  &=&",1pe13.4," [cm^3]"/
     &                     "  xmin  &=&",1pe13.4," [cm]"/
     &                     "  xmax  &=&",1pe13.4," [cm]"/
     &                     "  ymin  &=&",1pe13.4," [cm]"/
     &                     "  ymax  &=&",1pe13.4," [cm]"/
     &                     "  zmin  &=&",1pe13.4," [cm]"/
     &                     "  zmax  &=&",1pe13.4," [cm]"/
     &                     "  emin  &=&",1pe13.4," [MeV/n]"/
     &                     "  emax  &=&",1pe13.4," [MeV/n]")')
     &                     yen, vm, xm(ix), xm(ix+nxstepi),
     &                  ym(iy), ym(iy+nystepi), zm(iz), zm(iz+nzstepi),
     &                     eb(ie), eb(ie+nestepi)
                  end if
                else if( ite2l(m).eq.1 ) then
                     write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "   vol  &=&",1pe13.4," [cm^3]"/
     &                     "  xmin  &=&",1pe13.4," [cm]"/
     &                     "  xmax  &=&",1pe13.4," [cm]"/
     &                     "  ymin  &=&",1pe13.4," [cm]"/
     &                     "  ymax  &=&",1pe13.4," [cm]"/
     &                     "  zmin  &=&",1pe13.4," [cm]"/
     &                     "  zmax  &=&",1pe13.4," [cm]"/
     &                     "  emin  &=&",1pe13.4," [keV/um]"/
     &                     "  emax  &=&",1pe13.4," [keV/um]")')
     &                     yen, vm, xm(ix), xm(ix+nxstepi),
     &                  ym(iy), ym(iy+nystepi), zm(iz), zm(iz+nzstepi),
     &                     eb(ie), eb(ie+nestepi)
                end if
               if(iloopmode .ge. 1) then
                     write(iot,'("  part  &=& ",a)') chp(ipi)(3:10)
               endif
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &                     "  tmin  &=&",1pe13.4," [nsec]"/
     &                     "  tmax  &=&",1pe13.4," [nsec]")')
     &                     tb(it), tb(it+ntstepi)
               end if

                  write(iot,'("e:")')

            end do

            end do
            end do
            end do
            end do
            end do

*-----------------------------------------------------------------------
*        x axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 3 ) then

               inum = 0

            do iei = 1, ne, nestepi
            do iti = 1, nt, ntstepi
            do iai = 1, na, nastepi
            do iyi = 1, ny, nystepi
            do izi = 1, nz, nzstepi
            do ipi = 1, np, npstepi
               ie = iei
               it = iti
               ia = iai
               iy = iyi
               iz = izi
               ip = ipi

               inum = inum + 1

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

               if(iloopmode .ge. 1) then
                  write(iot,'("#   ",a)') chp(ipi)(3:10)
               endif

                  write(iot,'("#   no. =",i3,3x,
     &            "ie  =",i3,3x,
     &            "iy  =",i3,3x,
     &            "iz  =",i3,/
     &            "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )"/
     &            "#   y = (",1p1e13.4,"  -",1p1e13.4,"  )"/
     &            "#   z = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        inum, ie, iy, iz,
     &                        eb(ie), eb(ie+nestepi),
     &                        ym(iy), ym(iy+nystepi),
     &                        zm(iz), zm(iz+nzstepi)
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  it =",i3/
     &            "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                  it, tb(it), tb(it+ntstepi)
               end if

               if( itaty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  ia =",i3,/
     &            "# ",a3," = (",
     &                           1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     ia, aname, ab(ia), ab(ia+nastepi)
               end if

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: x [cm]")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: Number ",a32)') hsunit(itunt(m))

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

               if( itxty(m) .eq. 3 .or. itxty(m) .eq. 5 ) then

                  write(iot,'( "p: xlog ylog afac(0.8) form(0.9)")')

               else

                  write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')

               end if

               if( itanl(m) .gt. 0 ) then

                  write(iot,'( "p: ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               end if

               if( itsans(m) .gt. 0 ) then

                  call write_sangel(iot,m,0)

               end if

       include 'samepage_include/petaxyz_x.inc'

*-----------------------------------------------------------------------

             if(iloopmode .eq. 0) then
               if( ittty(m) .ne. 0 .and. itaty(m) .ne. 0 ) then
                  write(iot,'(/a1,"no. =",i3,
     &                     ",   e =",i3,
     &                     ",   t =",i3,
     &                     "    a =",i3,
     &                     ",   y =",i3,
     &                     ",   z =",i3,a1)')
     &                     cha, inum, ie, it, ia, iy, iz, cha
               else if( ittty(m) .eq. 0 .and. itaty(m) .ne. 0 ) then
                  write(iot,'(/a1,"no. =",i3,
     &                     ",   e =",i3,
     &                     ",   a =",i3,
     &                     ",   y =",i3,
     &                     ",   z =",i3,a1)')
     &                     cha, inum, ie, ia, iy, iz, cha
               else if( ittty(m) .ne. 0 .and. itaty(m) .eq. 0 ) then
                  write(iot,'(/a1,"no. =",i3,
     &                     ",   e =",i3,
     &                     ",   t =",i3,
     &                     ",   y =",i3,
     &                     ",   z =",i3,a1)')
     &                     cha, inum, ie, it, iy, iz, cha
               else
                  write(iot,'(/a1,"no. =",i3,
     &                     ",   e =",i3,
     &                     ",   y =",i3,
     &                     ",   z =",i3,a1)')
     &                     cha, inum, ie, iy, iz, cha
               end if

             else
               write(changelsub(1),'(a1,"no. =",i3)') cha,inum
               write(changelsub(2),'(",  part = ",a)') chp(ipi)(3:10)
               write(changelsub(3),'(",  e =",i3)') ie
               write(changelsub(4),'(",  t =",i3)') it
               write(changelsub(5),'(",  a =",i3)') ia
               write(changelsub(6),'(",  x =",i3)') ix
               write(changelsub(7),'(",  y =",i3)') iy
               write(changelsub(8),'(",  z =",i3)') iz
               write(changelsub(9),'(a1)') cha

               if( ne .eq. 0 ) then
                 changelsub(3) = " "
               end if
               if( ittty(m) .eq. 0 ) then
                 changelsub(4) = " "
               end if
               if( itaty(m) .eq. 0 ) then
                 changelsub(5) = " "
               end if

               changelsub(6) = " "

               if(iloopmode .eq. 1) then
                 changelsub(3) = " "
               end if
               if(iloopmode .eq. 10) then
                 changelsub(4) = " "
               end if
               if(iloopmode .eq. 8 .or. iloopmode .eq. 9) then
                 changelsub(5) = " "
               end if
               if(iloopmode .eq. 3) then
                 changelsub(6) = " "
               end if
               if(iloopmode .eq. 4) then
                 changelsub(7) = " "
               end if
               if(iloopmode .eq. 5) then
                 changelsub(8) = " "
               end if
                angeltitle = trim(changelsub(1))//trim(changelsub(2))//
     &                       trim(changelsub(3))//trim(changelsub(4))//
     &                       trim(changelsub(5))//trim(changelsub(6))//
     &                       trim(changelsub(7))//trim(changelsub(8))//
     &                       trim(changelsub(9))
                write(iot,'(/a)') trim(angeltitle)
          end if

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

                if( ite2l(m).eq.0 ) then  ! ccse 2022/09/30
                  if( iMeVperu.eq.0 ) then
                     write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "  emin  &=&",1pe13.4," [MeV]"/
     &                     "  emax  &=&",1pe13.4," [MeV]"/
     &                     "  ymin  &=&",1pe13.4," [cm]"/
     &                     "  ymax  &=&",1pe13.4," [cm]"/
     &                     "  zmin  &=&",1pe13.4," [cm]"/
     &                     "  zmax  &=&",1pe13.4," [cm]")')
     &                     yen, eb(ie), eb(ie+nestepi),
     &                  ym(iy), ym(iy+nystepi), zm(iz), zm(iz+nzstepi)
                  else
                     write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "  emin  &=&",1pe13.4," [MeV/n]"/
     &                     "  emax  &=&",1pe13.4," [MeV/n]"/
     &                     "  ymin  &=&",1pe13.4," [cm]"/
     &                     "  ymax  &=&",1pe13.4," [cm]"/
     &                     "  zmin  &=&",1pe13.4," [cm]"/
     &                     "  zmax  &=&",1pe13.4," [cm]")')
     &                     yen, eb(ie), eb(ie+nestepi),
     &                  ym(iy), ym(iy+nystepi), zm(iz), zm(iz+nzstepi)
                  end if
                else if( ite2l(m).eq.1 ) then
                     write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "  emin  &=&",1pe13.4," [keV/um]"/
     &                     "  emax  &=&",1pe13.4," [keV/um]"/
     &                     "  ymin  &=&",1pe13.4," [cm]"/
     &                     "  ymax  &=&",1pe13.4," [cm]"/
     &                     "  zmin  &=&",1pe13.4," [cm]"/
     &                     "  zmax  &=&",1pe13.4," [cm]")')
     &                     yen, eb(ie), eb(ie+nestepi),
     &                  ym(iy), ym(iy+nystepi), zm(iz), zm(iz+nzstepi)
                end if
               if(iloopmode .ge. 1) then
                     write(iot,'("  part  &=& ",a)') chp(ipi)(3:10)
               endif

               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &                     "  tmin  &=&",1pe13.4," [nsec]"/
     &                     "  tmax  &=&",1pe13.4," [nsec]")')
     &                     tb(it), tb(it+ntstepi)
               end if
               if( itaty(m) .ne. 0 ) then
                  write(iot,'(
     &                     "  amin  &=&",1pe13.4,/
     &                     "  amax  &=&",1pe13.4)')
     &                     ab(ia), ab(ia+nastepi)
               end if

                  write(iot,'("e:")')

            end do

            end do
            end do
            end do
            end do
            end do

*-----------------------------------------------------------------------
*        y axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 4 ) then

               inum = 0

            do iei = 1, ne, nestepi
            do iti = 1, nt, ntstepi
            do iai = 1, na, nastepi
            do ixi = 1, nx, nxstepi
            do izi = 1, nz, nzstepi
            do ipi = 1, np, npstepi
               ie = iei
               it = iti
               ia = iai
               ix = ixi
               iz = izi
               ip = ipi

               inum = inum + 1

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

               if(iloopmode .ge. 1) then
                  write(iot,'("#   ",a)') chp(ipi)(3:10)
               endif

                  write(iot,'("#   no. =",i3,3x,
     &            "ie  =",i3,3x,
     &            "ix  =",i3,3x,
     &            "iz  =",i3,/
     &            "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )"/
     &            "#   x = (",1p1e13.4,"  -",1p1e13.4,"  )"/
     &            "#   z = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        inum, ie, ix, iz,
     &                        eb(ie), eb(ie+nestepi),
     &                        xm(ix), xm(ix+nxstepi),
     &                        zm(iz), zm(iz+nzstepi)
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  it =",i3/
     &            "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                  it, tb(it), tb(it+ntstepi)
               end if

               if( itaty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  ia =",i3,/
     &            "# ",a3," = (",
     &                           1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     ia, aname, ab(ia), ab(ia+nastepi)
               end if

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: y [cm]")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: Number ",a32)') hsunit(itunt(m))

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

               if( ityty(m) .eq. 3 .or. ityty(m) .eq. 5 ) then

                  write(iot,'( "p: xlog ylog afac(0.8) form(0.9)")')

               else

                  write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')

               end if

               if( itanl(m) .gt. 0 ) then

                  write(iot,'( "p: ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               end if

               if( itsans(m) .gt. 0 ) then

                  call write_sangel(iot,m,0)

               end if

       include 'samepage_include/petaxyz_y.inc'

*-----------------------------------------------------------------------

             if(iloopmode .eq. 0) then
               if( ittty(m) .ne. 0 .and. itaty(m) .ne. 0 ) then
                  write(iot,'(/a1,"no. =",i3,
     &                     ",   e =",i3,
     &                     ",   t =",i3,
     &                     "    a =",i3,
     &                     ",   x =",i3,
     &                     ",   z =",i3,a1)')
     &                     cha, inum, ie, it, ia, ix, iz, cha
               else if( ittty(m) .eq. 0 .and. itaty(m) .ne. 0 ) then
                  write(iot,'(/a1,"no. =",i3,
     &                     ",   e =",i3,
     &                     ",   a =",i3,
     &                     ",   x =",i3,
     &                     ",   z =",i3,a1)')
     &                     cha, inum, ie, ia, ix, iz, cha
               else if( ittty(m) .ne. 0 .and. itaty(m) .eq. 0 ) then
                  write(iot,'(/a1,"no. =",i3,
     &                     ",   e =",i3,
     &                     ",   t =",i3,
     &                     ",   x =",i3,
     &                     ",   z =",i3,a1)')
     &                     cha, inum, ie, it, ix, iz, cha
               else
                  write(iot,'(/a1,"no. =",i3,
     &                     ",   e =",i3,
     &                     ",   x =",i3,
     &                     ",   z =",i3,a1)')
     &                     cha, inum, ie, ix, iz, cha
               end if

             else
               write(changelsub(1),'(a1,"no. =",i3)') cha,inum
               write(changelsub(2),'(",  part = ",a)') chp(ipi)(3:10)
               write(changelsub(3),'(",  e =",i3)') ie
               write(changelsub(4),'(",  t =",i3)') it
               write(changelsub(5),'(",  a =",i3)') ia
               write(changelsub(6),'(",  x =",i3)') ix
               write(changelsub(7),'(",  y =",i3)') iy
               write(changelsub(8),'(",  z =",i3)') iz
               write(changelsub(9),'(a1)') cha

               if( ne .eq. 0 ) then
                 changelsub(3) = " "
               end if
               if( ittty(m) .eq. 0 ) then
                 changelsub(4) = " "
               end if
               if( itaty(m) .eq. 0 ) then
                 changelsub(5) = " "
               end if

               changelsub(7) = " "

               if(iloopmode .eq. 1) then
                 changelsub(3) = " "
               end if
               if(iloopmode .eq. 10) then
                 changelsub(4) = " "
               end if
               if(iloopmode .eq. 8 .or. iloopmode .eq. 9) then
                 changelsub(5) = " "
               end if
               if(iloopmode .eq. 3) then
                 changelsub(6) = " "
               end if
               if(iloopmode .eq. 4) then
                 changelsub(7) = " "
               end if
               if(iloopmode .eq. 5) then
                 changelsub(8) = " "
               end if
                angeltitle = trim(changelsub(1))//trim(changelsub(2))//
     &                       trim(changelsub(3))//trim(changelsub(4))//
     &                       trim(changelsub(5))//trim(changelsub(6))//
     &                       trim(changelsub(7))//trim(changelsub(8))//
     &                       trim(changelsub(9))
                write(iot,'(/a)') trim(angeltitle)
          end if

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

                if( ite2l(m).eq.0 ) then  ! ccse 2022/09/30
                  if( iMeVperu.eq.0 ) then
                     write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "  emin  &=&",1pe13.4," [MeV]"/
     &                     "  emax  &=&",1pe13.4," [MeV]"/
     &                     "  xmin  &=&",1pe13.4," [cm]"/
     &                     "  xmax  &=&",1pe13.4," [cm]"/
     &                     "  zmin  &=&",1pe13.4," [cm]"/
     &                     "  zmax  &=&",1pe13.4," [cm]")')
     &                     yen, eb(ie), eb(ie+nestepi),
     &                  xm(ix), xm(ix+nxstepi), zm(iz), zm(iz+nzstepi)
                  else
                     write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "  emin  &=&",1pe13.4," [MeV/n]"/
     &                     "  emax  &=&",1pe13.4," [MeV/n]"/
     &                     "  xmin  &=&",1pe13.4," [cm]"/
     &                     "  xmax  &=&",1pe13.4," [cm]"/
     &                     "  zmin  &=&",1pe13.4," [cm]"/
     &                     "  zmax  &=&",1pe13.4," [cm]")')
     &                     yen, eb(ie), eb(ie+nestepi),
     &                  xm(ix), xm(ix+nxstepi), zm(iz), zm(iz+nzstepi)
                  end if
                else if( ite2l(m).eq.1 ) then
                     write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "  emin  &=&",1pe13.4," [keV/um]"/
     &                     "  emax  &=&",1pe13.4," [keV/um]"/
     &                     "  xmin  &=&",1pe13.4," [cm]"/
     &                     "  xmax  &=&",1pe13.4," [cm]"/
     &                     "  zmin  &=&",1pe13.4," [cm]"/
     &                     "  zmax  &=&",1pe13.4," [cm]")')
     &                     yen, eb(ie), eb(ie+nestepi),
     &                  xm(ix), xm(ix+nxstepi), zm(iz), zm(iz+nzstepi)
                end if
               if(iloopmode .ge. 1) then
                     write(iot,'("  part  &=& ",a)') chp(ipi)(3:10)
               endif

               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &                     "  tmin  &=&",1pe13.4," [nsec]"/
     &                     "  tmax  &=&",1pe13.4," [nsec]")')
     &                     tb(it), tb(it+ntstepi)
               end if
               if( itaty(m) .ne. 0 ) then
                  write(iot,'(
     &                     "  amin  &=&",1pe13.4,/
     &                     "  amax  &=&",1pe13.4)')
     &                     ab(ia), ab(ia+nastepi)
               end if

                  write(iot,'("e:")')

            end do

            end do
            end do
            end do
            end do
            end do

*-----------------------------------------------------------------------
*        z axis
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 5 ) then

               inum = 0

            do iei = 1, ne, nestepi
            do iti = 1, nt, ntstepi
            do iai = 1, na, nastepi
            do ixi = 1, nx, nxstepi
            do iyi = 1, ny, nystepi
            do ipi = 1, np, npstepi
               ie = iei
               it = iti
               ia = iai
               ix = ixi
               iy = iyi
               ip = ipi

               inum = inum + 1

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

               if(iloopmode .ge. 1) then
                  write(iot,'("#   ",a)') chp(ipi)(3:10)
               endif

                  write(iot,'("#   no. =",i3,3x,
     &            "ie  =",i3,3x,
     &            "ix  =",i3,3x,
     &            "iy  =",i3,/
     &            "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )"/
     &            "#   x = (",1p1e13.4,"  -",1p1e13.4,"  )"/
     &            "#   y = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        inum, ie, ix, iy,
     &                        eb(ie), eb(ie+nestepi),
     &                        xm(ix), xm(ix+nxstepi),
     &                        ym(iy), ym(iy+nystepi)
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  it =",i3/
     &            "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                  it, tb(it), tb(it+ntstepi)
               end if

               if( itaty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  ia =",i3,/
     &            "# ",a3," = (",
     &                           1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     ia, aname, ab(ia), ab(ia+nastepi)
               end if

               if( itaxl(m) .eq. 0 ) then

                  write(iot,'(/"x: z[cm]")')

               else

                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .eq. 0 ) then

                  write(iot,'( "y: Number ",a32)') hsunit(itunt(m))

               else

                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

               if( itzty(m) .eq. 3 .or. itzty(m) .eq. 5 ) then

                  write(iot,'( "p: xlog ylog afac(0.8) form(0.9)")')

               else

                  write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')

               end if

               if( itanl(m) .gt. 0 ) then

                  write(iot,'( "p: ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               end if

               if( itsans(m) .gt. 0 ) then

                  call write_sangel(iot,m,0)

               end if

       include 'samepage_include/petaxyz_z.inc'

*-----------------------------------------------------------------------

             if(iloopmode .eq. 0) then
               if( ittty(m) .ne. 0 .and. itaty(m) .ne. 0 ) then
                  write(iot,'(/a1,"no. =",i3,
     &                     ",   e =",i3,
     &                     ",   t =",i3,
     &                     "    a =",i3,
     &                     ",   x =",i3,
     &                     ",   y =",i3,a1)')
     &                     cha, inum, ie, it, ia, ix, iy, cha
               else if( ittty(m) .eq. 0 .and. itaty(m) .ne. 0 ) then
                  write(iot,'(/a1,"no. =",i3,
     &                     ",   e =",i3,
     &                     ",   a =",i3,
     &                     ",   x =",i3,
     &                     ",   y =",i3,a1)')
     &                     cha, inum, ie, ia, ix, iy, cha
               else if( ittty(m) .ne. 0 .and. itaty(m) .eq. 0 ) then
                  write(iot,'(/a1,"no. =",i3,
     &                     ",   e =",i3,
     &                     ",   t =",i3,
     &                     ",   x =",i3,
     &                     ",   y =",i3,a1)')
     &                     cha, inum, ie, it, ix, iy, cha
               else
                  write(iot,'(/a1,"no. =",i3,
     &                     ",   e =",i3,
     &                     ",   x =",i3,
     &                     ",   y =",i3,a1)')
     &                     cha, inum, ie, ix, iy, cha
               end if

             else
               write(changelsub(1),'(a1,"no. =",i3)') cha,inum
               write(changelsub(2),'(",  part = ",a)') chp(ipi)(3:10)
               write(changelsub(3),'(",  e =",i3)') ie
               write(changelsub(4),'(",  t =",i3)') it
               write(changelsub(5),'(",  a =",i3)') ia
               write(changelsub(6),'(",  x =",i3)') ix
               write(changelsub(7),'(",  y =",i3)') iy
               write(changelsub(8),'(",  z =",i3)') iz
               write(changelsub(9),'(a1)') cha

               if( ne .eq. 0 ) then
                 changelsub(3) = " "
               end if
               if( ittty(m) .eq. 0 ) then
                 changelsub(4) = " "
               end if
               if( itaty(m) .eq. 0 ) then
                 changelsub(5) = " "
               end if

               changelsub(8) = " "

               if(iloopmode .eq. 1) then
                 changelsub(3) = " "
               end if
               if(iloopmode .eq. 10) then
                 changelsub(4) = " "
               end if
               if(iloopmode .eq. 8 .or. iloopmode .eq. 9) then
                 changelsub(5) = " "
               end if
               if(iloopmode .eq. 3) then
                 changelsub(6) = " "
               end if
               if(iloopmode .eq. 4) then
                 changelsub(7) = " "
               end if
               if(iloopmode .eq. 5) then
                 changelsub(8) = " "
               end if
                angeltitle = trim(changelsub(1))//trim(changelsub(2))//
     &                       trim(changelsub(3))//trim(changelsub(4))//
     &                       trim(changelsub(5))//trim(changelsub(6))//
     &                       trim(changelsub(7))//trim(changelsub(8))//
     &                       trim(changelsub(9))
                write(iot,'(/a)') trim(angeltitle)
          end if

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

                if( ite2l(m).eq.0 ) then  ! ccse 2022/09/30
                  if( iMeVperu.eq.0 ) then
                     write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "  emin  &=&",1pe13.4," [MeV]"/
     &                     "  emax  &=&",1pe13.4," [MeV]"/
     &                     "  xmin  &=&",1pe13.4," [cm]"/
     &                     "  xmax  &=&",1pe13.4," [cm]"/
     &                     "  ymin  &=&",1pe13.4," [cm]"/
     &                     "  ymax  &=&",1pe13.4," [cm]")')
     &                     yen, eb(ie), eb(ie+nestepi),
     &                  xm(ix), xm(ix+nxstepi), ym(iy), ym(iy+nystepi)
                  else
                     write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "  emin  &=&",1pe13.4," [MeV/n]"/
     &                     "  emax  &=&",1pe13.4," [MeV/n]"/
     &                     "  xmin  &=&",1pe13.4," [cm]"/
     &                     "  xmax  &=&",1pe13.4," [cm]"/
     &                     "  ymin  &=&",1pe13.4," [cm]"/
     &                     "  ymax  &=&",1pe13.4," [cm]")')
     &                     yen, eb(ie), eb(ie+nestepi),
     &                  xm(ix), xm(ix+nxstepi), ym(iy), ym(iy+nystepi)
                  end if
                else if( ite2l(m).eq.1 ) then
                     write(iot,'("wt: s(0.7)",/
     &                     a1,"vspace{-3}"/
     &                     "  emin  &=&",1pe13.4," [keV/um]"/
     &                     "  emax  &=&",1pe13.4," [keV/um]"/
     &                     "  xmin  &=&",1pe13.4," [cm]"/
     &                     "  xmax  &=&",1pe13.4," [cm]"/
     &                     "  ymin  &=&",1pe13.4," [cm]"/
     &                     "  ymax  &=&",1pe13.4," [cm]")')
     &                     yen, eb(ie), eb(ie+nestepi),
     &                  xm(ix), xm(ix+nxstepi), ym(iy), ym(iy+nystepi)
                end if
               if(iloopmode .ge. 1) then
                     write(iot,'("  part  &=& ",a)') chp(ipi)(3:10)
               endif

               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &                     "  tmin  &=&",1pe13.4," [nsec]"/
     &                     "  tmax  &=&",1pe13.4," [nsec]")')
     &                     tb(it), tb(it+ntstepi)
               end if
               if( itaty(m) .ne. 0 ) then
                  write(iot,'(
     &                     "  amin  &=&",1pe13.4,/
     &                     "  amax  &=&",1pe13.4)')
     &                     ab(ia), ab(ia+nastepi)
               end if

                  write(iot,'("e:")')

            end do

            end do
            end do
            end do
            end do
            end do

*-----------------------------------------------------------------------
*        xy axis (matrix)
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 7 ) then

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

            do iz = 1, nz
            do ip = 1, npg
            do ie = 1, neg
            do it = 1, ntg
            do ia = 1, nag

               inum = inum + 1

                  write(iot,'(/"#",78("-"))')

               if( inum .eq. 1 .or. ip.gt.np_mxang ) then ! frtati 2021/10/05
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

               if( rtfac(m) .lt. 0.d0 )
     &           write(iot,
     &           '( "#restart maximum for normalization:",1p1e14.7)')
     &           facmax(m)

               if ( ip.gt.np_mxang ) then
                 write(iot,'( " SKIPPAGE:")')
                 inum = inum - 1
               end if

                  write(iot,'("#   no. =",i3,3x,
     &            "ie  =",i3,3x,
     &            "iz  =",i3,3x,
     &            "part. = ",a8,/
     &            "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )"/
     &            "#   z = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        inum, ie, iz, chq(ip),
     &                        eb(ie), eb(ie+1),
     &                        zm(iz), zm(iz+1)
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  it =",i3/
     &            "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                  it, tb(it), tb(it+1)
               end if

               if( itaty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  ia =",i3,/
     &            "# ",a3," = (",
     &                           1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     ia, aname, ab(ia), ab(ia+1)
               end if

*-----------------------------------------------------------------------

            if( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) then

               if( ittty(m) .ne. 0 .and. itaty(m) .ne. 0 ) then
                  write(iot,'(/a1,"no. =",i3,
     &                     ",   e =",i3,
     &                     ",   t =",i3,
     &                     "    a =",i3,
     &                     ",   z =",i3,a1)')
     &                     cha, inum, ie, it, ia, iz, cha
               else if( ittty(m) .eq. 0 .and. itaty(m) .ne. 0 ) then
                  write(iot,'(/a1,"no. =",i3,
     &                     ",   e =",i3,
     &                     ",   a =",i3,
     &                     ",   z =",i3,a1)')
     &                     cha, inum, ie, ia, iz, cha
               else if( ittty(m) .ne. 0 .and. itaty(m) .eq. 0 ) then
                  write(iot,'(/a1,"no. =",i3,
     &                     ",   e =",i3,
     &                     ",   t =",i3,
     &                     ",   z =",i3,a1)')
     &                     cha, inum, ie, it, iz, cha
               else
                  write(iot,'(/a1,"no. =",i3,
     &                     ",   e =",i3,
     &                     ",   z =",i3,a1)')
     &                     cha, inum, ie, iz, cha
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

*-----------------------------------------------------------------------

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

                  write(iot,'("#  ny = ",i3,"   nx = ",i3)')
     &                         ny, nx

            if( ittwo(m) .ne. 4 ) then

                  write(iot,'( "# ( ( data(x,y), x = 1, nx ),",
     &                         " y = ny, 1, -1 )")')

            end if

*-----------------------------------------------------------------------

            if( ( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) .and.
     &            igsh .eq. 0 ) then

                  write(iot,'(/a4," y = ",1p1g14.7," to ",
     &            1p1g14.7," by ",1p1g14.7," ; x = ",
     &            1p1g14.7," to ",1p1g14.7," by ",
     &            1p1g14.7," ;")') dc2,
     &            ym(ny) + rtydl(m)/2.0, ym(1) + rtydl(m)/2.0, rtydl(m),
     &            xm(1) + rtxdl(m)/2.0, xm(nx) + rtxdl(m)/2.0, rtxdl(m)

!OBINATA(2012.9.11): output *.err
               write(iot,'(1p10e11.3)')
     &         ( ( tr(ip,ie,it,ia,icf(ix,iy,iz),ioe), ix = 1, nx ),
     &                                      iy = ny, 1, -1 )

*-----------------------------------------------------------------------

            else if( ittwo(m) .eq. 4 ) then

               write(iot,'(/"# x          y        ",
     &                      "  number     r.err")')

               do iy = 1, ny
               do ix = 1, nx

                  write(iot,'(1p3e11.3,0pf8.4)')
     &               xm(ix)  + rtxdl(m)/2.0,
     &               ym(iy)  + rtydl(m)/2.0,
     &               tr(ip,ie,it,ia,icf(ix,iy,iz),1),
     &               tr(ip,ie,it,ia,icf(ix,iy,iz),2)

               end do
               end do

*-----------------------------------------------------------------------

            else if( ittwo(m) .eq. 5 ) then

                  write(iot,'("#   x = ",1p1g14.7," to ",
     &            1p1g14.7," by ",1p1g14.7/
     &                        "#   y = ",1p1g14.7," to ",
     &            1p1g14.7," by ",1p1g14.7)')
     &            xm(1) + rtxdl(m)/2.0, xm(nx) + rtxdl(m)/2.0, rtxdl(m),
     &            ym(1) + rtydl(m)/2.0, ym(ny) + rtydl(m)/2.0, rtydl(m)

                  write(iot,'(/4x,a3,4x,1p1000e11.3)')
     &            'y/x',( xm(ix) + rtxdl(m)/2.0, ix = 1, nx )

               do iy = ny, 1, -1

!OBINATA(2012.9.11): output *.err
                  write(iot,'(1p1000e11.3)')
     &            ym(iy) + rtydl(m)/2.0,
     &            ( tr(ip,ie,it,ia,icf(ix,iy,iz),ioe), ix = 1, nx )

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

               zval = ( zm(iz) + zm(iz+1) ) / 2.0d0
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
     &    .and. cmin .gt. 0.0 .and. cmax .gt. cmin .and.
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

         write(iot,'("y: Number ",a32)') hsunit(itunt(m))

      else

         write(iot,'("y: ",200a1)') (itazt(m)(i:i),i=1,itazl(m))

      end if

      if( izlog .eq. 0 ) write(iot,'( "p: ylin")')

      write(iot,'( "#",78("-"))')

      end if

*-----------------------------------------------------------------------

            if( inocm .eq. 1 ) then

                if( ite2l(m).eq.0 ) then  ! ccse 2022/09/30
                  if( iMeVperu.eq.0 ) then
                     write(iot,'(/"wt: s[c5]",/
     &                     a1,"vspace{-3}"/
     &                     "  emin  &=&",1pe13.4," [MeV]"/
     &                     "  emax  &=&",1pe13.4," [MeV]"/
     &                     "  zmin  &=&",1pe13.4," [cm]"/
     &                     "  zmax  &=&",1pe13.4," [cm]"/
     &                     "  part. &=&   ",a8)')
     &                     yen, eb(ie), eb(ie+1),
     &                     zm(iz), zm(iz+1), chq(ip)
                  else
                     write(iot,'(/"wt: s[c5]",/
     &                     a1,"vspace{-3}"/
     &                     "  emin  &=&",1pe13.4," [MeV/n]"/
     &                     "  emax  &=&",1pe13.4," [MeV/n]"/
     &                     "  zmin  &=&",1pe13.4," [cm]"/
     &                     "  zmax  &=&",1pe13.4," [cm]"/
     &                     "  part. &=&   ",a8)')
     &                     yen, eb(ie), eb(ie+1),
     &                     zm(iz), zm(iz+1), chq(ip)
                  end if
                else if( ite2l(m).eq.1 ) then
                     write(iot,'(/"wt: s[c5]",/
     &                     a1,"vspace{-3}"/
     &                     "  emin  &=&",1pe13.4," [keV/um]"/
     &                     "  emax  &=&",1pe13.4," [keV/um]"/
     &                     "  zmin  &=&",1pe13.4," [cm]"/
     &                     "  zmax  &=&",1pe13.4," [cm]"/
     &                     "  part. &=&   ",a8)')
     &                     yen, eb(ie), eb(ie+1),
     &                     zm(iz), zm(iz+1), chq(ip)
                end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &                     "  tmin  &=&",1pe13.4," [nsec]"/
     &                     "  tmax  &=&",1pe13.4," [nsec]")')
     &                     tb(it), tb(it+1)
               end if
               if( itaty(m) .ne. 0 ) then
                  write(iot,'(
     &                     "  amin  &=&",1pe13.4,/
     &                     "  amax  &=&",1pe13.4)')
     &                     ab(ia), ab(ia+1)
               end if

                  write(iot,'("e:")')

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

            end do
            end do
            end do
            end do
            end do

*-----------------------------------------------------------------------

            if ( itbmp(m) .ne. 0 ) then

               numIndex = 5
               allocate( bmpfIType(1:numIndex) )
               allocate( bmpfIndex(1:numIndex) )
               bmpfIType = (/ 'a', 't', 'e', 'p', 'z' /)
               bmpfIndex = (/ nag, ntg, neg, npg, nz /)
               bmpWidth  = nx
               bmpHeight = ny

            end if

            if ( itvtk(m) .ne. 0 ) then

               call open_file(isunit_vtk_meta_default,
     &                 "", isunit_vtk_meta, ios, .true.)
               call open_file(isunit_vtk_default,
     &                 "", isunit_vtk, ios, .true.)

               iaxs = 1
               nparam = 3

               write(isunit_vtk_meta) ntg

               do it = 1, ntg

                  write(isunit_vtk_meta) iaxs
                  write(isunit_vtk_meta) nparam
                  write(isunit_vtk_meta) 'p', 'e', 'a'
                  write(isunit_vtk_meta) npg, neg, nag

                  write(isunit_vtk_meta) nx,ny,nz
                  write(isunit_vtk_meta)
     &                    ( xm(ix), ix=1,nx+1 )
                  write(isunit_vtk_meta)
     &                    ( ym(iy), iy=1,ny+1 )
                  write(isunit_vtk_meta)
     &                    ( zm(iz), iz=1,nz+1 )

                  do ip = 1, npg
                  do ie = 1, neg
                  do ia = 1, nag

                     write(isunit_vtk_meta) ip, ie, ia

                     write(isunit_vtk)
     &                       ( ( ( tr(ip,ie,it,ia,icf(ix,iy,iz),ioe),
     &                             ix=1,nx ),
     &                             iy=1,ny ),
     &                             iz=1,nz )

                  end do
                  end do
                  end do

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

         else if( itaxs(m,iax) .eq. 8 ) then

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

            do ix = 1, nx
            do ip = 1, npg
            do ie = 1, neg
            do it = 1, ntg
            do ia = 1, nag

               inum = inum + 1

                  write(iot,'(/"#",78("-"))')

               if( inum .eq. 1 .or. ip.gt.np_mxang ) then ! frtati 2021/10/05
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

               if( rtfac(m) .lt. 0.d0 )
     &           write(iot,
     &           '( "#restart maximum for normalization:",1p1e14.7)')
     &           facmax(m)

               if ( ip.gt.np_mxang ) then
                 write(iot,'( " SKIPPAGE:")')
                 inum = inum - 1
               end if

                  write(iot,'( "#   no. =",i3,3x,
     &            "ie  =",i3,3x,
     &            "ix  =",i3,3x,
     &            "part. = ",a8,/
     &            "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )"/
     &            "#   x = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        inum, ie, ix, chq(ip),
     &                        eb(ie), eb(ie+1),
     &                        xm(ix), xm(ix+1)
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  it =",i3/
     &            "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                  it, tb(it), tb(it+1)
               end if

               if( itaty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  ia =",i3,/
     &            "# ",a3," = (",
     &                           1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     ia, aname, ab(ia), ab(ia+1)
               end if

*-----------------------------------------------------------------------

            if( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) then

               if( ittty(m) .ne. 0 .and. itaty(m) .ne. 0 ) then
                  write(iot,'(/a1,"no. =",i3,
     &                     ",   e =",i3,
     &                     ",   t =",i3,
     &                     "    a =",i3,
     &                     ",   x =",i3,a1)')
     &                     cha, inum, ie, it, ia, ix, cha
               else if( ittty(m) .eq. 0 .and. itaty(m) .ne. 0 ) then
                  write(iot,'(/a1,"no. =",i3,
     &                     ",   e =",i3,
     &                     ",   a =",i3,
     &                     ",   x =",i3,a1)')
     &                     cha, inum, ie, ia, ix, cha
               else if( ittty(m) .ne. 0 .and. itaty(m) .eq. 0 ) then
                  write(iot,'(/a1,"no. =",i3,
     &                     ",   e =",i3,
     &                     ",   t =",i3,
     &                     ",   x =",i3,a1)')
     &                     cha, inum, ie, it, ix, cha
               else
                  write(iot,'(/a1,"no. =",i3,
     &                     ",   e =",i3,
     &                     ",   x =",i3,a1)')
     &                     cha, inum, ie, ix, cha
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

                  write(iot,'("#  ny = ",i3,"   nz = ",i3)')
     &                         ny, nz

            if( ittwo(m) .ne. 4 ) then

                  write(iot,'( "# ( ( data(z,y), z = 1, nz ),",
     &                         " y = ny, 1, -1 )")')

            end if

*-----------------------------------------------------------------------

            if( ( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) .and.
     &            igsh .eq. 0 ) then

                  write(iot,'(/a4," y = ",1p1g14.7," to ",
     &            1p1g14.7," by ",1p1g14.7," ; x = ",
     &            1p1g14.7," to ",1p1g14.7," by ",
     &            1p1g14.7," ;")') dc2,
     &            ym(ny) + rtydl(m)/2.0, ym(1) + rtydl(m)/2.0, rtydl(m),
     &            zm(1) + rtzdl(m)/2.0, zm(nz) + rtzdl(m)/2.0, rtzdl(m)

!OBINATA(2012.9.11): output *.err
               write(iot,'(1p10e11.3)')
     &         ( ( tr(ip,ie,it,ia,icf(ix,iy,iz),ioe), iz = 1, nz ),
     &                                      iy = ny, 1, -1 )

*-----------------------------------------------------------------------

            else if( ittwo(m) .eq. 4 ) then

               write(iot,'(/"# y          z        ",
     &                      "  number     r.err")')

               do iz = 1, nz
               do iy = 1, ny

                  write(iot,'(1p3e11.3,0pf8.4)')
     &               ym(iy)  + rtydl(m)/2.0,
     &               zm(iz)  + rtzdl(m)/2.0,
     &               tr(ip,ie,it,ia,icf(ix,iy,iz),1),
     &               tr(ip,ie,it,ia,icf(ix,iy,iz),2)

               end do
               end do

*-----------------------------------------------------------------------

            else if( ittwo(m) .eq. 5 ) then

                  write(iot,'("#   y = ",1p1g14.7," to ",
     &            1p1g14.7," by ",1p1g14.7/
     &                        "#   z = ",1p1g14.7," to ",
     &            1p1g14.7," by ",1p1g14.7)')
     &            ym(1) + rtydl(m)/2.0, ym(ny) + rtydl(m)/2.0, rtydl(m),
     &            zm(1) + rtzdl(m)/2.0, zm(nz) + rtzdl(m)/2.0, rtzdl(m)

                  write(iot,'(/4x,a3,4x,1p1000e11.3)')
     &            'y/z',( zm(iz) + rtzdl(m)/2.0, iz = 1, nz )

               do iy = ny, 1, -1

!OBINATA(2012.9.11): output *.err
                  write(iot,'(1p1000e11.3)')
     &            ym(iy) + rtydl(m)/2.0,
     &            ( tr(ip,ie,it,ia,icf(ix,iy,iz),ioe), iz = 1, nz )

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

               xval = ( xm(ix) + xm(ix+1) ) / 2.0d0
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
     &    .and. cmin .gt. 0.0 .and. cmax .gt. cmin .and.
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

         write(iot,'("y: Number ",a32)') hsunit(itunt(m))

      else

         write(iot,'("y: ",200a1)') (itazt(m)(i:i),i=1,itazl(m))

      end if

      if( izlog .eq. 0 ) write(iot,'( "p: ylin")')

      write(iot,'( "#",78("-"))')

      end if

*-----------------------------------------------------------------------

            if( inocm .eq. 1 ) then

                if( ite2l(m).eq.0 ) then  ! ccse 2022/09/30
                  if( iMeVperu.eq.0 ) then
                     write(iot,'(/"wt: s[c5]",/
     &                     a1,"vspace{-3}"/
     &                     "  emin  &=&",1pe13.4," [MeV]"/
     &                     "  emax  &=&",1pe13.4," [MeV]"/
     &                     "  xmin  &=&",1pe13.4," [cm]"/
     &                     "  xmax  &=&",1pe13.4," [cm]"/
     &                     "  part. &=&   ",a8)')
     &                     yen, eb(ie), eb(ie+1),
     &                     xm(ix), xm(ix+1), chq(ip)
                  else
                     write(iot,'(/"wt: s[c5]",/
     &                     a1,"vspace{-3}"/
     &                     "  emin  &=&",1pe13.4," [MeV/n]"/
     &                     "  emax  &=&",1pe13.4," [MeV/n]"/
     &                     "  xmin  &=&",1pe13.4," [cm]"/
     &                     "  xmax  &=&",1pe13.4," [cm]"/
     &                     "  part. &=&   ",a8)')
     &                     yen, eb(ie), eb(ie+1),
     &                     xm(ix), xm(ix+1), chq(ip)
                  end if
                else if( ite2l(m).eq.1 ) then
                     write(iot,'(/"wt: s[c5]",/
     &                     a1,"vspace{-3}"/
     &                     "  emin  &=&",1pe13.4," [keV/um]"/
     &                     "  emax  &=&",1pe13.4," [keV/um]"/
     &                     "  xmin  &=&",1pe13.4," [cm]"/
     &                     "  xmax  &=&",1pe13.4," [cm]"/
     &                     "  part. &=&   ",a8)')
     &                     yen, eb(ie), eb(ie+1),
     &                     xm(ix), xm(ix+1), chq(ip)
                end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &                     "  tmin  &=&",1pe13.4," [nsec]"/
     &                     "  tmax  &=&",1pe13.4," [nsec]")')
     &                     tb(it), tb(it+1)
               end if
               if( itaty(m) .ne. 0 ) then
                  write(iot,'(
     &                     "  amin  &=&",1pe13.4,/
     &                     "  amax  &=&",1pe13.4)')
     &                     ab(ia), ab(ia+1)
               end if

                  write(iot,'("e:")')

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

            end do
            end do
            end do
            end do
            end do

*-----------------------------------------------------------------------

            if ( itbmp(m) .ne. 0 ) then

               numIndex = 5
               allocate( bmpfIType(1:numIndex) )
               allocate( bmpfIndex(1:numIndex) )
               bmpfIType = (/ 'a', 't', 'e', 'p', 'x' /)
               bmpfIndex = (/ nag, ntg, neg, npg, nx /)
               bmpWidth  = nz
               bmpHeight = ny

            end if

            if ( itvtk(m) .ne. 0 ) then

               call open_file(isunit_vtk_meta_default,
     &                 "", isunit_vtk_meta, ios, .true.)
               call open_file(isunit_vtk_default,
     &                 "", isunit_vtk, ios, .true.)

               iaxs = 2
               nparam = 3

               write(isunit_vtk_meta) ntg

               do it = 1, ntg

                  write(isunit_vtk_meta) iaxs
                  write(isunit_vtk_meta) nparam
                  write(isunit_vtk_meta) 'p', 'e', 'a'
                  write(isunit_vtk_meta) npg, neg, nag

                  write(isunit_vtk_meta) nx,ny,nz
                  write(isunit_vtk_meta)
     &                    ( xm(ix), ix=1,nx+1 )
                  write(isunit_vtk_meta)
     &                    ( ym(iy), iy=1,ny+1 )
                  write(isunit_vtk_meta)
     &                    ( zm(iz), iz=1,nz+1 )

                  do ip = 1, npg
                  do ie = 1, neg
                  do ia = 1, nag

                     write(isunit_vtk_meta) ip, ie, ia

                     write(isunit_vtk)
     &                       ( ( ( tr(ip,ie,it,ia,icf(ix,iy,iz),ioe),
     &                             ix=1,nx ),
     &                             iy=1,ny ),
     &                             iz=1,nz )

                  end do
                  end do
                  end do

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

            do iy = 1, ny
            do ip = 1, npg
            do ie = 1, neg
            do it = 1, ntg
            do ia = 1, nag

               inum = inum + 1

                  write(iot,'(/"#",78("-"))')

               if( inum .eq. 1 .or. ip.gt.np_mxang ) then ! frtati 2021/10/05
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

               if( rtfac(m) .lt. 0.d0 )
     &           write(iot,
     &           '( "#restart maximum for normalization:",1p1e14.7)')
     &           facmax(m)

               if ( ip.gt.np_mxang ) then
                 write(iot,'( " SKIPPAGE:")')
                 inum = inum - 1
               end if

                  write(iot,'("#   no. =",i3,3x,
     &            "ie  =",i3,3x,
     &            "iy  =",i3,3x,
     &            "part. = ",a8,/
     &            "#   e = (",1p1e13.4,"  -",1p1e13.4,"  )"/
     &            "#   y = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                        inum, ie, iy, chq(ip),
     &                        eb(ie), eb(ie+1),
     &                        ym(iy), ym(iy+1)
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  it =",i3/
     &            "#   t = (",1p1e13.4,"  -",1p1e13.4,"  )")')
     &                  it, tb(it), tb(it+1)
               end if

               if( itaty(m) .ne. 0 ) then
                  write(iot,'(
     &            "#  ia =",i3,/
     &            "# ",a3," = (",
     &                           1p1e13.4,"  -",1p1e13.4,"  )")')
     &                     ia, aname, ab(ia), ab(ia+1)
               end if

*-----------------------------------------------------------------------

            if( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) then

               if( ittty(m) .ne. 0 .and. itaty(m) .ne. 0 ) then
                  write(iot,'(/a1,"no. =",i3,
     &                     ",   e =",i3,
     &                     ",   t =",i3,
     &                     "    a =",i3,
     &                     ",   y =",i3,a1)')
     &                     cha, inum, ie, it, ia, iy, cha
               else if( ittty(m) .eq. 0 .and. itaty(m) .ne. 0 ) then
                  write(iot,'(/a1,"no. =",i3,
     &                     ",   e =",i3,
     &                     ",   a =",i3,
     &                     ",   y =",i3,a1)')
     &                     cha, inum, ie, ia, iy, cha
               else if( ittty(m) .ne. 0 .and. itaty(m) .eq. 0 ) then
                  write(iot,'(/a1,"no. =",i3,
     &                     ",   e =",i3,
     &                     ",   t =",i3,
     &                     ",   y =",i3,a1)')
     &                     cha, inum, ie, it, iy, cha
               else
                  write(iot,'(/a1,"no. =",i3,
     &                     ",   e =",i3,
     &                     ",   y =",i3,a1)')
     &                     cha, inum, ie, iy, cha
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

                  write(iot,'("#  nx = ",i3,"   nz = ",i3)')
     &                         nx, nz

            if( ittwo(m) .ne. 4 ) then

                  write(iot,'( "# ( ( data(z,x), z = 1, nz ),",
     &                         " x = nx, 1, -1 )")')

            end if

*-----------------------------------------------------------------------

            if( ( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) .and.
     &            igsh .eq. 0 ) then

                  write(iot,'(/a4," y = ",1p1g14.7," to ",
     &            1p1g14.7," by ",1p1g14.7," ; x = ",
     &            1p1g14.7," to ",1p1g14.7," by ",
     &            1p1g14.7," ;")') dc2,
     &            xm(nx) + rtxdl(m)/2.0, xm(1) + rtxdl(m)/2.0, rtxdl(m),
     &            zm(1) + rtzdl(m)/2.0, zm(nz) + rtzdl(m)/2.0, rtzdl(m)

!OBINATA(2012.9.11): output *.err
               write(iot,'(1p10e11.3)')
     &         ( ( tr(ip,ie,it,ia,icf(ix,iy,iz),ioe), iz = 1, nz ),
     &                                      ix = nx, 1, -1 )

*-----------------------------------------------------------------------

            else if( ittwo(m) .eq. 4 ) then

               write(iot,'(/"# x          z        ",
     &                      "  number     r.err")')

               do iz = 1, nz
               do ix = 1, nx

                  write(iot,'(1p3e11.3,0pf8.4)')
     &               xm(ix)  + rtxdl(m)/2.0,
     &               zm(iz)  + rtzdl(m)/2.0,
     &               tr(ip,ie,it,ia,icf(ix,iy,iz),1),
     &               tr(ip,ie,it,ia,icf(ix,iy,iz),2)

               end do
               end do

*-----------------------------------------------------------------------

            else if( ittwo(m) .eq. 5 ) then

                  write(iot,'("#   x = ",1p1g14.7," to ",
     &            1p1g14.7," by ",1p1g14.7/
     &                        "#   z = ",1p1g14.7," to ",
     &            1p1g14.7," by ",1p1g14.7)')
     &            xm(1) + rtxdl(m)/2.0, xm(nx) + rtxdl(m)/2.0, rtxdl(m),
     &            zm(1) + rtzdl(m)/2.0, zm(nz) + rtzdl(m)/2.0, rtzdl(m)

                  write(iot,'(/4x,a3,4x,1p1000e11.3)')
     &            'x/z',( zm(iz) + rtzdl(m)/2.0, iz = 1, nz )

               do ix = nx, 1, -1

!OBINATA(2012.9.11): output *.err
                  write(iot,'(1p1000e11.3)')
     &            xm(ix) + rtxdl(m)/2.0,
     &            ( tr(ip,ie,it,ia,icf(ix,iy,iz),ioe), iz = 1, nz )

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

               yval = ( ym(iy) + ym(iy+1) ) / 2.0d0
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
     &    .and. cmin .gt. 0.0 .and. cmax .gt. cmin .and.
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

         write(iot,'("y: Number ",a32)') hsunit(itunt(m))

      else

         write(iot,'("y: ",200a1)') (itazt(m)(i:i),i=1,itazl(m))

      end if

      if( izlog .eq. 0 ) write(iot,'( "p: ylin")')

      write(iot,'( "#",78("-"))')

      end if

*-----------------------------------------------------------------------

            if( inocm .eq. 1 ) then

                if( ite2l(m).eq.0 ) then  ! ccse 2022/09/30
                  if( iMeVperu.eq.0 ) then
                     write(iot,'(/"wt: s[c5]",/
     &                     a1,"vspace{-3}"/
     &                     "  emin  &=&",1pe13.4," [MeV]"/
     &                     "  emax  &=&",1pe13.4," [MeV]"/
     &                     "  ymin  &=&",1pe13.4," [cm]"/
     &                     "  ymax  &=&",1pe13.4," [cm]"/
     &                     "  part. &=&   ",a8)')
     &                     yen, eb(ie), eb(ie+1),
     &                     ym(iy), ym(iy+1), chq(ip)
                  else
                     write(iot,'(/"wt: s[c5]",/
     &                     a1,"vspace{-3}"/
     &                     "  emin  &=&",1pe13.4," [MeV/n]"/
     &                     "  emax  &=&",1pe13.4," [MeV/n]"/
     &                     "  ymin  &=&",1pe13.4," [cm]"/
     &                     "  ymax  &=&",1pe13.4," [cm]"/
     &                     "  part. &=&   ",a8)')
     &                     yen, eb(ie), eb(ie+1),
     &                     ym(iy), ym(iy+1), chq(ip)
                  end if
                else if( ite2l(m).eq.1 ) then
                     write(iot,'(/"wt: s[c5]",/
     &                     a1,"vspace{-3}"/
     &                     "  emin  &=&",1pe13.4," [keV/um]"/
     &                     "  emax  &=&",1pe13.4," [keV/um]"/
     &                     "  ymin  &=&",1pe13.4," [cm]"/
     &                     "  ymax  &=&",1pe13.4," [cm]"/
     &                     "  part. &=&   ",a8)')
     &                     yen, eb(ie), eb(ie+1),
     &                     ym(iy), ym(iy+1), chq(ip)
                end if
               if( ittty(m) .ne. 0 ) then
                  write(iot,'(
     &                     "  tmin  &=&",1pe13.4," [nsec]"/
     &                     "  tmax  &=&",1pe13.4," [nsec]")')
     &                     tb(it), tb(it+1)
               end if
               if( itaty(m) .ne. 0 ) then
                  write(iot,'(
     &                     "  amin  &=&",1pe13.4,/
     &                     "  amax  &=&",1pe13.4)')
     &                     ab(ia), ab(ia+1)
               end if

                  write(iot,'("e:")')

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

            end do
            end do
            end do
            end do
            end do

*-----------------------------------------------------------------------

            if ( itbmp(m) .ne. 0 ) then

               numIndex = 5
               allocate( bmpfIType(1:numIndex) )
               allocate( bmpfIndex(1:numIndex) )
               bmpfIType = (/ 'a', 't', 'e', 'p', 'y' /)
               bmpfIndex = (/ nag, ntg, neg, npg, ny /)
               bmpWidth  = nz
               bmpHeight = nx

            end if

            if ( itvtk(m) .ne. 0 ) then

               call open_file(isunit_vtk_meta_default,
     &                 "", isunit_vtk_meta, ios, .true.)
               call open_file(isunit_vtk_default,
     &                 "", isunit_vtk, ios, .true.)

               iaxs = 3
               nparam = 3

               write(isunit_vtk_meta) ntg

               do it = 1, ntg

                  write(isunit_vtk_meta) iaxs
                  write(isunit_vtk_meta) nparam
                  write(isunit_vtk_meta) 'p', 'e', 'a'
                  write(isunit_vtk_meta) npg, neg, nag

                  write(isunit_vtk_meta) nx,ny,nz
                  write(isunit_vtk_meta)
     &                    ( xm(ix), ix=1,nx+1 )
                  write(isunit_vtk_meta)
     &                    ( ym(iy), iy=1,ny+1 )
                  write(isunit_vtk_meta)
     &                    ( zm(iz), iz=1,nz+1 )

                  do ip = 1, npg
                  do ie = 1, neg
                  do ia = 1, nag

                     write(isunit_vtk_meta) ip, ie, ia

                     write(isunit_vtk)
     &                       ( ( ( tr(ip,ie,it,ia,icf(ix,iy,iz),ioe),
     &                             ix=1,nx ),
     &                             iy=1,ny ),
     &                             iz=1,nz )

                  end do
                  end do
                  end do

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

         end if             ! itaxs

            call prestart(m,iot) !OBINATA(2012.9.11)

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

         if ( allocated(bmpfIType) ) deallocate( bmpfIType )
         if ( allocated(bmpfIndex) ) deallocate( bmpfIndex )

         call close_file(isunit_vtk_rm)
         call close_file(isunit_vtk_geom)
         call close_file(isunit_vtk_geom_meta)

*-----------------------------------------------------------------------

  900 continue

*-----------------------------------------------------------------------
      deallocate( ixyz )
cnais 2023/01/31
      include 'samepage_include/samepage999.inc'

      deallocate (vl_x,vl_y,vl_z)

      return
      end


************************************************************************
*                                                                      *
      subroutine pgshxyz(m,nx,ny,nz,xm,ym,zm,idasa)
*                                                                      *
*       output xyz scoring mesh gshow tally                            *
*       last modified by K.Niita on 2004/09/15                         *
*                                                                      *
************************************************************************
      use sangelmod, only: itsans

*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'

*-----------------------------------------------------------------------

      common /taliin/ rsouin, nzztin, nrgnin

      common /verjam/ versn, lastr, iyeav, imonv, idayv

*-----------------------------------------------------------------------

      common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                itout(itlmax), ittwo(itlmax)
      common /tall10/ itaxn(itlmax), itaxs(itlmax,6)
      common /tall11/ itfln(itlmax), itfll(itlmax,6), ctfln(itlmax,6)
      character ctfln*100

      common /tall14/ ittll(itlmax), ittle(itlmax)
      character ittle*80
      common /tall15/ itanl(itlmax), itang(itlmax)
      character itang*200

      common /tall21/ rtfac(itlmax)

      common /tall27/ itaxl(itlmax), itaxt(itlmax),
     &                itayl(itlmax), itayt(itlmax),
     &                itazl(itlmax), itazt(itlmax)
      character itaxt*200, itayt*200, itazt*200
      common /tall35/ iteps(itlmax)
      common /tall38/ itres(itlmax)
      common /tall39/ rtwid(itlmax)
      common /tall46/ itmtr(itlmax,4), rtmtr(itlmax,13)
      common /tall49/ itglt(itlmax)

*-----------------------------------------------------------------------

      common /cparm/  maxbch,maxcas
      common /jcomon/ nabov,nobch,nocas,nomax
!$OMP THREADPRIVATE(/jcomon/)

      common /talout/ itall

      character fname*100, fnume*3

*-----------------------------------------------------------------------
      character gfnam*100
      common /tall60/ itger(itlmax)

*-----------------------------------------------------------------------

      dimension   xm(nx+1)
      dimension   ym(ny+1)
      dimension   zm(nz+1)

      integer,allocatable :: ixyz(:)

*-----------------------------------------------------------------------

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

      integer :: isunit_vtk_geom_default = 94
      integer :: isunit_vtk_geom_meta_default = 95
      integer :: iunit_vtk_g_default = 96
      integer :: isunit_vtk_geom = 0
      integer :: isunit_vtk_geom_meta = 0
      integer :: iunit_vtk_g = 0
      integer :: ios

      character(len=255) :: outFilename
      logical :: isText

*-----------------------------------------------------------------------

      yen  = char(92)
      allocate( ixyz(max((nx+1)*(ny+1),(ny+1)*(nz+1),(nx+1)*(nz+1))) )

*-----------------------------------------------------------------------
*     out put unit = 31  : temporary number
*-----------------------------------------------------------------------

      do iax = 1, itfln(m)

         if( itall .eq. 2 .and. nobch .lt. maxbch ) then

            write(fnume,'(i3.3)') nobch
            fname = ctfln(m,iax)(1:itfll(m,iax))//'.'//fnume

         else

            fname = ctfln(m,iax)

         end if

            iot = 31
            open(iot, file = fname, status = 'unknown' )

*-----------------------------------------------------------------------
            igsh  = 1
            gfnam = ctfln(m,iax)
            igfmn = itfll(m,iax)
            igser = itger(m)

*-----------------------------------------------------------------------
*     input echo
*-----------------------------------------------------------------------

               call tgshech(iot,m,iax,1)

*-----------------------------------------------------------------------
*        xy axis (matrix)
*-----------------------------------------------------------------------

         if( itaxs(m,iax) .eq. 1 ) then

*-----------------------------------------------------------------------

               inum = 0

            do iz = 1, nz

               zval = ( zm(iz) + zm(iz+1) ) / 2.0d0

               inum = inum + 1

                  write(iot,'(/"#",78("-"))')

cfrtati 2022/03/22 restoring
               if( inum .eq. 1 ) then
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

               write(iot,'("#   no. =",i3,3x,"z = ",1p1e13.4)')
     &                     inum, zval

*-----------------------------------------------------------------------

               write(iot,'(/a1,"no. =",i3,
     &                     ",    z = ",1p1e13.4,a1)')
     &                     cha, inum, zval, cha

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

               write(iot,'(/"p: legs[c5*0.875]")')

*-----------------------------------------------------------------------

               none = 1

               call gshow(0,iot,itaxs(m,iax),itout(m),rtwid(m),itres(m),
     &                    inum,igsh,igser,gfnam,igfmn,
     &                    nx+1,ny+1,none,xm,ym,zval,ixyz(1),
     &                    1,1,krr,vll,itmtr(m,4),itglt(m))

*-----------------------------------------------------------------------

            end do

            if ( itbmp(m) .ne. 0 ) then

               numIndex = 1
               allocate( bmpfIType(1:numIndex) )
               allocate( bmpfIndex(1:numIndex) )
               bmpfIType = (/ 'z' /)
               bmpfIndex = (/ nz /)
               bmpWidth  = nx
               bmpHeight = ny

            end if

*-----------------------------------------------------------------------
*        yz axis (matrix)
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 2 ) then

*-----------------------------------------------------------------------

               inum = 0

            do ix = 1, nx

               xval = ( xm(ix) + xm(ix+1) ) / 2.0d0

               inum = inum + 1

                  write(iot,'(/"#",78("-"))')

cfrtati 2022/03/22 restoring
               if( inum .eq. 1 ) then
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

               write(iot,'("#   no. =",i3,3x,"x = ",1p1e13.4)')
     &                     inum, xval

*-----------------------------------------------------------------------

               write(iot,'(/a1,"no. =",i3,
     &                     ",    x = ",1p1e13.4,a1)')
     &                     cha, inum, xval, cha

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

*-----------------------------------------------------------------------

               xmin = zm(1)
               xmax = zm(nz+1)
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

               write(iot,'(/"p: legs[c5*0.875]")')

*-----------------------------------------------------------------------

               none = 1

               call gshow(0,iot,itaxs(m,iax),itout(m),rtwid(m),itres(m),
     &                    inum,igsh,igser,gfnam,igfmn,
     &                    nz+1,ny+1,none,zm,ym,xval,ixyz(1),
     &                    1,1,krr,vll,itmtr(m,4),itglt(m))

*-----------------------------------------------------------------------

            end do

            if ( itbmp(m) .ne. 0 ) then

               numIndex = 1
               allocate( bmpfIType(1:numIndex) )
               allocate( bmpfIndex(1:numIndex) )
               bmpfIType = (/ 'x' /)
               bmpfIndex = (/ nx /)
               bmpWidth  = nz
               bmpHeight = ny

            end if

*-----------------------------------------------------------------------
*        xz axis (matrix)
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 3 ) then

*-----------------------------------------------------------------------

               inum = 0

            do iy = 1, ny

               yval = ( ym(iy) + ym(iy+1) ) / 2.0d0

               inum = inum + 1

                  write(iot,'(/"#",78("-"))')

cfrtati 2022/03/22 restoring
               if( inum .eq. 1 ) then
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

               write(iot,'("#   no. =",i3,3x,"x = ",1p1e13.4)')
     &                     inum, yval

*-----------------------------------------------------------------------

               write(iot,'(/a1,"no. =",i3,
     &                     ",    y = ",1p1e13.4,a1)')
     &                     cha, inum, yval, cha

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

*-----------------------------------------------------------------------

               xmin = zm(1)
               xmax = zm(nz+1)
               ymin = xm(1)
               ymax = xm(nx+1)

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

               write(iot,'(/"p: legs[c5*0.875]")')

*-----------------------------------------------------------------------

               none = 1

               call gshow(0,iot,itaxs(m,iax),itout(m),rtwid(m),itres(m),
     &                    inum,igsh,igser,gfnam,igfmn,
     &                    nz+1,nx+1,none,zm,xm,yval,ixyz(1),
     &                    1,1,krr,vll,itmtr(m,4),itglt(m))

*-----------------------------------------------------------------------

            end do

            if ( itbmp(m) .ne. 0 ) then

               numIndex = 1
               allocate( bmpfIType(1:numIndex) )
               allocate( bmpfIndex(1:numIndex) )
               bmpfIType = (/ 'y' /)
               bmpfIndex = (/ ny /)
               bmpWidth  = nz
               bmpHeight = nx

            end if

*-----------------------------------------------------------------------

         end if

            close(iot)

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

            iaxs = itaxs(m,iax)
            iuni = itout(m)
            ires = itres(m)

            call open_file(isunit_vtk_geom_default,
     &              "", isunit_vtk_geom, ios, .true.)
            call open_file(isunit_vtk_geom_meta_default,
     &              "", isunit_vtk_geom_meta, ios, .true.)

            call vtk_set_gshow(
     &              isunit_vtk_geom, isunit_vtk_geom_meta,
     &              nx, ny, nz, xm, ym, zm,
     &              iaxs, iuni, ires, igser,
     &              1, 1, krr, vll, itmtr(m,4))

            isText = (itvtkfmt(m).eq.0)

            rewind(isunit_vtk_geom)
            rewind(isunit_vtk_geom_meta)

            call vtk_create_filename(
     &              fname, outFilename, .false., 0, 0)
            call open_file(iunit_vtk_g_default,
     &              outFilename, iunit_vtk_g, ios, isText)
            call vtk_write_polydata(
     &              itvtkfmt(m),
     &              iunit_vtk_g,
     &              isunit_vtk_geom, isunit_vtk_geom_meta)
            call close_file(iunit_vtk_g)

         end if

         if ( allocated(bmpfIType) ) deallocate( bmpfIType )
         if ( allocated(bmpfIndex) ) deallocate( bmpfIndex )

         call close_file(isunit_vtk_geom)
         call close_file(isunit_vtk_geom_meta)

*-----------------------------------------------------------------------

      end do
      deallocate( ixyz)
      
*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine prshxyz(m,nx,ny,nz,xm,ym,zm,
     &                   nr,mr,kr,nvl,ivl,rvl,idasa)
*                                                                      *
*       output xyz scoring mesh rshow tally                            *
*       last modified by K.Niita on 2004/09/15                         *
*                                                                      *
************************************************************************
      use sangelmod, only: itsans

*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'

*-----------------------------------------------------------------------

      common /talmm/  nmmax, lmmax, itlmx
      common /taliin/ rsouin, nzztin, nrgnin

*-----------------------------------------------------------------------

      common /verjam/ versn, lastr, iyeav, imonv, idayv

      common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                itout(itlmax), ittwo(itlmax)
      common /tall10/ itaxn(itlmax), itaxs(itlmax,6)
      common /tall11/ itfln(itlmax), itfll(itlmax,6), ctfln(itlmax,6)
      character ctfln*100

      common /tall14/ ittll(itlmax), ittle(itlmax)
      character ittle*80
      common /tall15/ itanl(itlmax), itang(itlmax)
      character itang*200

      common /tall21/ rtfac(itlmax)

      common /tall27/ itaxl(itlmax), itaxt(itlmax),
     &                itayl(itlmax), itayt(itlmax),
     &                itazl(itlmax), itazt(itlmax)
      character itaxt*200, itayt*200, itazt*200
      common /tall35/ iteps(itlmax)
      common /tall38/ itres(itlmax)
      common /tall39/ rtwid(itlmax)
      common /tall46/ itmtr(itlmax,4), rtmtr(itlmax,13)
      common /tall49/ itglt(itlmax)

*-----------------------------------------------------------------------

      common /cparm/  maxbch,maxcas
      common /jcomon/ nabov,nobch,nocas,nomax
!$OMP THREADPRIVATE(/jcomon/)

      common /talout/ itall

      character fname*100, fnume*3

*-----------------------------------------------------------------------
      character gfnam*100
      common /tall60/ itger(itlmax)

*-----------------------------------------------------------------------

      dimension   xm(nx+1)
      dimension   ym(ny+1)
      dimension   zm(nz+1)
      dimension   kr(mr)
      dimension   vl(nr)
      dimension   lr(nr)
      dimension   ivl(nvl)
      dimension   rvl(nvl)

      integer,allocatable :: ixyz(:)

*-----------------------------------------------------------------------

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

*-----------------------------------------------------------------------

      yen  = char(92)
      allocate( ixyz(max((nx+1)*(ny+1),(ny+1)*(nz+1),(nx+1)*(nz+1))) )

*-----------------------------------------------------------------------
*     out put unit = 31  : temporary number
*-----------------------------------------------------------------------

      do iax = 1, itfln(m)

         if( itall .eq. 2 .and. nobch .lt. maxbch ) then

            write(fnume,'(i3.3)') nobch
            fname = ctfln(m,iax)(1:itfll(m,iax))//'.'//fnume

         else

            fname = ctfln(m,iax)

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

               call trshech(iot,m,iax,1)

*-----------------------------------------------------------------------
*     set values and normalization between 0 and 1
*-----------------------------------------------------------------------

               call tregval(nr,mr,kr,vl,lr,nvl,ivl,rvl)

                  cmax = -1.e+33
                  cmin =  1.e+33
                  dnon =  1.e-33

               do i = 1, nr

                  if( vl(i) .gt. cmax .and. vl(i) .gt. dnon )
     &                           cmax = vl(i)
                  if( vl(i) .lt. cmin .and. vl(i) .gt. dnon )
     &                           cmin = vl(i)

               end do

                  izlog = 1
                  if( itanl(m) .gt. 0 )
     &            call anset(itang(m),itanl(m),
     &                       form,xfac,afac,cmin,cmax,izlog,inocm,inolg)

            if( izlog .eq. 0 ) then

               do i = 1, nr

                  if( vl(i) .gt. dnon .and. cmin .ne. cmax ) then

                     vl(i) = ( vl(i) - cmin ) / ( cmax - cmin )
                     vl(i) = max(0.0d0,vl(i))
                     vl(i) = min(1.0d0,vl(i))

                  else

                     vl(i) = -1.0

                  end if

               end do

            else

               do i = 1, nr

                  if( vl(i) .gt. dnon .and. cmin .ne. cmax ) then

                     vl(i) = log10(vl(i)/cmin) / log10(cmax/cmin)
                     vl(i) = max(0.0d0,vl(i))
                     vl(i) = min(1.0d0,vl(i))

                  else

                     vl(i) = -1.0d0

                  end if

               end do

            end if

*-----------------------------------------------------------------------
*        xy axis (matrix)
*-----------------------------------------------------------------------

         if( itaxs(m,iax) .eq. 1 ) then

*-----------------------------------------------------------------------

               inum = 0

            do iz = 1, nz

               zval = ( zm(iz) + zm(iz+1) ) / 2.0d0

               inum = inum + 1

                  write(iot,'(/"#",78("-"))')

cfrtati 2022/03/22 restoring
               if( inum .eq. 1 ) then
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

               write(iot,'("#   no. =",i3,3x,"z = ",1p1e13.4)')
     &                     inum, zval

*-----------------------------------------------------------------------

               write(iot,'(/a1,"no. =",i3,
     &                     ",    z = ",1p1e13.4,a1)')
     &                     cha, inum, zval, cha

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

                  write(iot,'( "set: c3[",1p1e15.6,"] c4[",
     &                  1p1e15.6,"]")') cmin, cmax

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

               write(iot,'(/"p: legs[c5*0.875]")')

*-----------------------------------------------------------------------

               none = 1

               iuni = itout(m) * 2 - 1

               call gshow(1,iot,itaxs(m,iax),iuni,rtwid(m),itres(m),
     &                    inum,igsh,igser,gfnam,igfmn,
     &                    nx+1,ny+1,none,xm,ym,zval,ixyz(1),
     &                    nr,mr,kr,vl,itmtr(m,4),itglt(m))

*-----------------------------------------------------------------------

      if( inolg .eq. 1 .and.
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

         write(iot,'("y: Values")')

      else

         write(iot,'("y: ",200a1)') (itazt(m)(i:i),i=1,itazl(m))

      end if

      if( izlog .eq. 0 ) write(iot,'( "p: ylin")')

      write(iot,'( "#",78("-"))')

      end if

*-----------------------------------------------------------------------

            if( inocm .eq. 1 ) then

               write(iot,'(/"wt: s[c5]",/
     &                     a1,"vspace{-3}"/
     &                     "  z &=&",1pe13.4," [cm]"/
     &                     "e:")')
     &                     yen, zval

            end if

*-----------------------------------------------------------------------

      if( inolg .eq. 1 .and.
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

            if ( itbmp(m) .ne. 0 ) then

               numIndex = 1
               allocate( bmpfIType(1:numIndex) )
               allocate( bmpfIndex(1:numIndex) )
               bmpfIType = (/ 'z' /)
               bmpfIndex = (/ nz /)
               bmpWidth  = nx
               bmpHeight = ny

            end if

*-----------------------------------------------------------------------
*        yz axis (matrix)
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 2 ) then

*-----------------------------------------------------------------------

               inum = 0

            do ix = 1, nx

               xval = ( xm(ix) + xm(ix+1) ) / 2.0d0

               inum = inum + 1

                  write(iot,'(/"#",78("-"))')

cfrtati 2022/03/22 restoring
               if( inum .eq. 1 ) then
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

               write(iot,'("#   no. =",i3,3x,"x = ",1p1e13.4)')
     &                     inum, xval

*-----------------------------------------------------------------------

               write(iot,'(/a1,"no. =",i3,
     &                     ",    x = ",1p1e13.4,a1)')
     &                     cha, inum, xval, cha

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

                  write(iot,'( "set: c3[",1p1e15.6,"] c4[",
     &                  1p1e15.6,"]")') cmin, cmax

*-----------------------------------------------------------------------

               xmin = zm(1)
               xmax = zm(nz+1)
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

               write(iot,'(/"p: legs[c5*0.875]")')

*-----------------------------------------------------------------------

               none = 1

               iuni = itout(m) * 2 - 1

               call gshow(1,iot,itaxs(m,iax),iuni,rtwid(m),itres(m),
     &                    inum,igsh,igser,gfnam,igfmn,
     &                    nz+1,ny+1,none,zm,ym,xval,ixyz(1),
     &                    nr,mr,kr,vl,itmtr(m,4),itglt(m))

*-----------------------------------------------------------------------

      if( inolg .eq. 1 .and.
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

         write(iot,'("y: Values")')

      else

         write(iot,'("y: ",200a1)') (itazt(m)(i:i),i=1,itazl(m))

      end if

      if( izlog .eq. 0 ) write(iot,'( "p: ylin")')

      write(iot,'( "#",78("-"))')

      end if

*-----------------------------------------------------------------------

            if( inocm .eq. 1 ) then

               write(iot,'(/"wt: s[c5]",/
     &                     a1,"vspace{-3}"/
     &                     "  x &=&",1pe13.4," [cm]"/
     &                     "e:")')
     &                     yen, xval

            end if

*-----------------------------------------------------------------------

      if( inolg .eq. 1 .and.
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

            if ( itbmp(m) .ne. 0 ) then

               numIndex = 1
               allocate( bmpfIType(1:numIndex) )
               allocate( bmpfIndex(1:numIndex) )
               bmpfIType = (/ 'x' /)
               bmpfIndex = (/ nx /)
               bmpWidth  = nz
               bmpHeight = ny

            end if

*-----------------------------------------------------------------------
*        xz axis (matrix)
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 3 ) then

*-----------------------------------------------------------------------

               inum = 0

            do iy = 1, ny

               yval = ( ym(iy) + ym(iy+1) ) / 2.0d0

               inum = inum + 1

                  write(iot,'(/"#",78("-"))')

cfrtati 2022/03/22 restoring
               if( inum .eq. 1 ) then
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

               write(iot,'("#   no. =",i3,3x,"x = ",1p1e13.4)')
     &                     inum, yval

*-----------------------------------------------------------------------

               write(iot,'(/a1,"no. =",i3,
     &                     ",    y = ",1p1e13.4,a1)')
     &                     cha, inum, yval, cha

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

                  write(iot,'( "set: c3[",1p1e15.6,"] c4[",
     &                  1p1e15.6,"]")') cmin, cmax

*-----------------------------------------------------------------------

               xmin = zm(1)
               xmax = zm(nz+1)
               ymin = xm(1)
               ymax = xm(nx+1)

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

               write(iot,'(/"p: legs[c5*0.875]")')

*-----------------------------------------------------------------------

               none = 1

               iuni = itout(m) * 2 - 1

               call gshow(1,iot,itaxs(m,iax),iuni,rtwid(m),itres(m),
     &                    inum,igsh,igser,gfnam,igfmn,
     &                    nz+1,nx+1,none,zm,xm,yval,ixyz(1),
     &                    nr,mr,kr,vl,itmtr(m,4),itglt(m))

*-----------------------------------------------------------------------

      if( inolg .eq. 1 .and.
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

         write(iot,'("y: Values")')

      else

         write(iot,'("y: ",200a1)') (itazt(m)(i:i),i=1,itazl(m))

      end if

      if( izlog .eq. 0 ) write(iot,'( "p: ylin")')

      write(iot,'( "#",78("-"))')

      end if

*-----------------------------------------------------------------------

            if( inocm .eq. 1 ) then

               write(iot,'(/"wt: s[c5]",/
     &                     a1,"vspace{-3}"/
     &                     "  y &=&",1pe13.4," [cm]"/
     &                     "e:")')
     &                     yen, yval


            end if

*-----------------------------------------------------------------------

      if( inolg .eq. 1 .and.
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

            if ( itbmp(m) .ne. 0 ) then

               numIndex = 1
               allocate( bmpfIType(1:numIndex) )
               allocate( bmpfIndex(1:numIndex) )
               bmpfIType = (/ 'y' /)
               bmpfIndex = (/ ny /)
               bmpWidth  = nz
               bmpHeight = nx

            end if

*-----------------------------------------------------------------------

         end if

            close(iot)

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

         if ( allocated(bmpfIType) ) deallocate( bmpfIType )
         if ( allocated(bmpfIndex) ) deallocate( bmpfIndex )

*-----------------------------------------------------------------------

      end do
      deallocate( ixyz )
      
*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine gshow(ior,iot,icod,iout0,widt,ires,
     &                 jnum,igsh,igser,gfnam,igfmn,
     &                 nx,ny,nz,vx,vy,vz,imat,
     &                 nr,mr,kr,vl,mtrns,iglat0)
*                                                                      *
*       output geometry boundary                                       *
*       last modified by K.Niita on 2020/04/23                         *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'

*-----------------------------------------------------------------------

      parameter ( smin = 0.2 ) ! min is blue

*-----------------------------------------------------------------------

      common /tlgeom/ iblz1,iblz2
!$OMP THREADPRIVATE(/tlgeom/)
      common /tlglat/ ilev1,ilev2,ilat1(5,10),ilat2(5,10)
!$OMP THREADPRIVATE(/tlglat/)

      common /kmat1a/ mxmat, mxmat0, mxnel

      common /verjam/ versn, lastr, iyeav, imonv, idayv
      common /regdc/ idrg(kvlmax), idgr(kvmmax)

*-----------------------------------------------------------------------
      common /igsherr/ igsher, icl01, icl02
      character i3iop*3
      character iopfile*100
      character gfnam*100

      common /tcntl/  icntl, inucr

*-----------------------------------------------------------------------

      dimension   vx(nx)
      dimension   vy(ny)
      dimension   vz(nz)
      dimension   imat(nx,ny,nz)
      dimension   kr(mr)
      dimension   vl(nr)

*-----------------------------------------------------------------------

      dimension itcl(kvlmax)

      dimension x(3)
      dimension xc(3)
      dimension xd(3)

      dimension rcd(8,3)
      dimension ics(8,3)

      dimension icdx(8)
      dimension icdy(8)

      dimension icd0(3)
      dimension icd1(3)
      dimension icd2(3)

      data icd0 / 3, 1, 2 /
      data icd1 / 1, 3, 3 /
      data icd2 / 2, 2, 1 /

      data icdx / 1,  1,  0, -1, -1, -1,  0,  1/
      data icdy / 0,  1,  1,  1,  0, -1, -1, -1/

*-----------------------------------------------------------------------

      data dlg / 1.0d+10 /
      data epsb / 1.0d-08 /

      data dps1 / 0.12345678d-5 / ! T.Sato 2024/04/23
      data dps2 / 0.23456789d-5 / ! T.Sato 2024/04/23

*-----------------------------------------------------------------------

      character huni7*10
      data huni7 /'h: x ny1,0'/

      character huni3*10
      data huni3 /'  y2=[y1]('/

      character huni4*80
      character huni6*30
      character huni5*3

      character huni11*16
      character huni8*16
      data huni8 /'hb:   line clip '/
      character huni9*16
      data huni9 /'hb: noline clip '/

      data huni5 /'),i'/


      character nlat*20
      character mlat*20

      character ct(7)*3
      data ct /'zzz','zz ','z  ','   ','t  ',
     &         'tt ','ttt'/

      character coln*8

      common /mtnmcl/ dmhsb(-1:kvlmax,4),
     &                dmtnm(-1:kvlmax), dmtcl(-1:kvlmax),
     &                nmtnm(-1:kvlmax), nmtcl(-1:kvlmax)
      character dmtnm*80, dmtcl*30

      common /mtnreg/ dmhsg(-1:kvlmax),
     &                nmtng(-1:kvlmax), dmtng(-1:kvlmax)
      character dmtng*80

c added  Nais 2019/11/19
      common /gsline/ nowgshow, igsline

      dimension ercol(3,3)
      data  ercol / 1.0, 0.133, 1.0, 1.0, 1.0, 1.0, 0.6, 1.0, 0.0/
      
*-----------------------------------------------------------------------
cFURUTA20150714 TETRA
      integer nlat3,itetcl(10),itetvol,itetauto,itgchk
      common /tetf1/ nlat3,itetcl,itetvol,itetauto,itgchk

      integer iii0,kkk0
      common /itettal2/ iii0,kkk0
!$OMP THREADPRIVATE(/itettal2/)

      if(iout0.eq.9) return ! need not draw geometry, T.Sato 2019/02/20

! 1st time call this routine, delete *.err, file T.Sato 2016/2/16
      data ifirst/0/

*-----------------------------------------------------------------------

         iout = iout0
         iglat = iglat0
         nowgshow = 0

      if( igsline .eq. 0 ) then
         if( iglat .eq. - 1) then
            iglat = 0
         end if
      else if( igsline .eq. 1 ) then
         if( iglat .eq. - 1) then
            iglat = 1
         end if
      else if( igsline .eq. 2 ) then
         if( iglat .eq. - 1) then
            iglat = 2
         end if
      else if( igsline .eq. 3 ) then
         if( iglat .eq. - 1) then
            iglat = 3
         end if
      end if

         if( iglat .eq. 0 ) nowgshow = 2
         if( iglat .eq. 2 ) nowgshow = 2
         if( iglat .eq. 3 ) nowgshow = 1

         if( nowgshow .eq. 1 ) then ! T.Sato 2020/10/07
            if( iout .eq. 5 ) iout = 3
            if( iout .eq. 6 ) iout = 4
            if( iout .eq. 7 ) iout = 3
            if( iout .eq. 8 ) iout = 4
         end if

*-----------------------------------------------------------------------

      if(ifirst.eq.0) then
       do ii = igfmn, 1, -1
        if( gfnam(ii:ii) .eq. '.' ) goto 153
       end do
       iffmn = igfmn
       goto 156
  153  if( ii .eq. 1 ) then
        iffmn = igfmn
        goto 156
       end if
       iffmn = ii - 1
  156  continue
       iop = 20
       iopfile = gfnam(1:iffmn)//'_geo.out'
       open(iop, file = iopfile, status = 'unknown')
       close(iop,status='delete')
       ifirst=1
      endif

      if(nlat3.gt.0) call tetranowgshow !FURUTA20221209

*-----------------------------------------------------------------------
*     open temporary file
*-----------------------------------------------------------------------

            ioh = 18
            open(ioh,status='scratch',form='unformatted')

*-----------------------------------------------------------------------

         if( igser .gt. 0 .or.
     &       icntl .eq. 7 .or. icntl .eq. 8 .or.
     &       icntl .eq. 9 .or. icntl .eq. 10 ) then

            if( igser .gt. 0 ) igsher = igser

            if( icntl .eq. 7 .or. icntl .eq. 8 .or.
     &          icntl .eq. 9 .or. icntl .eq. 10 ) igsher = 2

            ioe = 19
            open(ioe,status='scratch',form='unformatted')

         end if

            inne1 = 0
            inne2 = 0
            inne3 = 0

*-----------------------------------------------------------------------
*           color of material
*-----------------------------------------------------------------------

            k = 0

         do i = 1, mxmat

            itcl(i) = 0

         end do

            ivodi = 0
            ivodo = 0

*-----------------------------------------------------------------------
*           width of line
*-----------------------------------------------------------------------

            ifac = nint( widt * 10.0 )

            if( ifac .le. 1 ) then
               itt = 1
            else if( ifac .le. 3 ) then
               itt = 2
            else if( ifac .le. 5 ) then
               itt = 3
            else if( ifac .le. 10 ) then
               itt = 4
            else if( ifac .le. 20 ) then
               itt = 5
            else if( ifac .le. 30 ) then
               itt = 6
            else
               itt = 7
            end if

cFURUTA20221209 move to below

*-----------------------------------------------------------------------
*        save initial mesh
*-----------------------------------------------------------------------

               ires2 = ires**2

               xmin0 = vx(1)
               ymin0 = vy(1)
               zmin0 = vz(1)
               if(nz.eq.1) then
                vz(1) = vz(1) + dps1*0.01 ! T.Sato 2024/04/23 
               else
                vz(1) = vz(1) + (vz(nz)-vz(1)) * dps1 ! T.Sato 2024/04/23 
               endif

               xdel0 = ( vx(nx) - vx(1) ) / dble( nx - 1 )
               ydel0 = ( vy(ny) - vy(1) ) / dble( ny - 1 )

               xdel1 = ( vx(nx) - vx(1) ) / dble( ires )
               ydel1 = ( vy(ny) - vy(1) ) / dble( ires )

               xdeld = ( vx(nx) - vx(1) ) / dble( nx )
               ydeld = ( vy(ny) - vy(1) ) / dble( ny )

               xdel = ( xdel1 + 2.0 * xdeld ) / dble( nx - 1 )
               ydel = ( ydel1 + 2.0 * ydeld ) / dble( ny - 1 )

               ddel = min( xdel/dble(nx), ydel/dble(ny) )**2

!FURUTA20221209 initialization move to here
*-----------------------------------------------------------------------
*        initialization
*        icod = 1(xy), 2(yz), 3(xz)
*-----------------------------------------------------------------------
            do i = 1, 8

               srt = sqrt( (xdel*icdx(i))**2 + (ydel*icdy(i))**2 )

               rcd(i,icd0(icod)) = 0.0d0
               rcd(i,icd1(icod)) = xdel/srt*icdx(i)
               rcd(i,icd2(icod)) = ydel/srt*icdy(i)

               ics(i,1) = icdx(i)
               ics(i,2) = icdy(i)
               ics(i,3) = 0

            end do

*-----------------------------------------------------------------------
*     do loop 900 for fine resolution
*-----------------------------------------------------------------------

      do 900 ire = 1, ires2

               ire1 = ire - ( ire - 1 ) / ires * ires
               ire2 = ( ire - 1 ) / ires + 1

               xmin = xmin0 + xdel1 * dble(ire1-1) - xdeld + xdel1*dps1 ! T.Sato 2024/04/23
               ymin = ymin0 + ydel1 * dble(ire2-1) - ydeld + ydel1*dps2 ! T.Sato 2024/04/23

               do i = 1, nx
                  vx(i) = xmin + xdel * dble( i - 1 )
               end do

               do i = 1, ny
                  vy(i) = ymin + ydel * dble( i - 1 )
               end do

*-----------------------------------------------------------------------
*        min and max of x, y and frame
*-----------------------------------------------------------------------

               call gval(x,1,1,1,icod,nx,ny,nz,vx,vy,vz)

                  xvmin = x(icd1(icod))
                  yvmin = x(icd2(icod))

               call gval(x,nx,ny,1,icod,nx,ny,nz,vx,vy,vz)

                  xvmax = x(icd1(icod))
                  yvmax = x(icd2(icod))

*-----------------------------------------------------------------------

                  xminf = xmin + xdeld - xdeld*dps1 ! T.Sato 2024/04/23
                  yminf = ymin + ydeld - ydeld*dps2 ! T.Sato 2024/04/23
                  xmaxf = xminf + xdel1
                  ymaxf = yminf + ydel1

*-----------------------------------------------------------------------
*        initialization of imat
*-----------------------------------------------------------------------
                  
               do i = 1, nx
               do j = 1, ny
               do k = 1, nz

                  imat(i,j,k) = 0

               end do
               end do
               end do

*-----------------------------------------------------------------------
cKN 2014/01/22 write error positions
*-----------------------------------------------------------------------

            if( igsher .ne. 0 ) then

                     ic = 1

               do k = 1, nz
               do j = 1, ny
               do i = 1, nx

                     call gval(x,i,j,k,icod,nx,ny,nz,vx,vy,vz)

                     u  = rcd(ic,1)
                     v  = rcd(ic,2)
                     w  = rcd(ic,3)

                     mark  =  1
                     markp =  0
                     ici   = -1

                     call gomsort(x(1),x(2),x(3),u,v,w,
     &                            nmed,iblz,mark,markp,ici,mtrns)

                     if( mark .le. -2 ) then

                        if( mark .eq. -2 ) then
                           ieid = 2
                           inne2 = inne2 + 1
                           write(ioe) ieid, i, j, k, 0, 0
                        else if( mark .eq. -3 ) then
                           ieid = 3
                           inne3 = inne3 + 1
                           write(ioe) ieid, i, j, k, icl01, icl02
                        else
                           ieid = 1
                           inne1 = inne1 + 1
                           write(ioe) ieid, i, j, k, 0, 0
                        end if

                     end if

               end do
               end do
               end do

            end if

*-----------------------------------------------------------------------
*     start check each cell
*-----------------------------------------------------------------------

               inum = 0

  100 continue

                     ic = 1

               do k = 1, nz
               do j = 1, ny
               do i = 1, nx

                  if( imat(i,j,k) .eq. 0 ) then

                     call gval(x,i,j,k,icod,nx,ny,nz,vx,vy,vz)

                     u  = rcd(ic,1)
                     v  = rcd(ic,2)
                     w  = rcd(ic,3)

                     mark  =  1
                     markp =  0
                     ici   = -1

                     call gomsort(x(1),x(2),x(3),u,v,w,
     &                            nmed,iblz,mark,markp,ici,mtrns)

                     if( mark .gt. -2 ) then

                        goto 200

                     else

                        imat(i,j,k) = 90000

                     end if

                  end if

               end do
               end do
               end do

*-----------------------------------------------------------------------
cKN 2014/01/22 write error positions
*-----------------------------------------------------------------------

         if( igsher .ne. 0 ) then
         if( inne1 + inne2 + inne3 .gt. 0 ) then

            inne = inne1 + inne2 + inne3

            rewind(ioe)

            if( inne1 .gt. 0 ) then
               write(iot,'(/a)') huni7//
     &              huni3//'error points'//huni5//'c[darkred]'
               x1 = -1.0
               y1 = -1.0
               write(iot,'(1p2e14.6)') x1, y1
            end if

            if( inne2 .gt. 0 ) then
               write(iot,'(/a)') huni7//
     &              huni3//'undefined'//huni5//'c[violet]'
               x1 = -1.0
               y1 = -1.0
               write(iot,'(1p2e14.6)') x1, y1
            end if

            if( inne3 .gt. 0 ) then
               write(iot,'(/a)') huni7//
     &              huni3//'double defined'//huni5//'c[black]'
               x1 = -1.0
               y1 = -1.0
               write(iot,'(1p2e14.6)') x1, y1
            end if

               write(iot,'("hb:")')
               write(iot,'("bmap: width=",e14.6," hight=",e14.6)')
     &                       xdel, ydel

            if( igsher .gt. 1 ) then

                  do ii = igfmn, 1, -1
                     if( gfnam(ii:ii) .eq. '.' ) goto 53
                  end do
                     iffmn = igfmn
                     goto 56
   53             if( ii .eq. 1 ) then
                     iffmn = igfmn
                     goto 56
                  end if
                     iffmn = ii - 1
   56             continue

! T.Sato 2016/2/16, always open as append
               iopfile = gfnam(1:iffmn)//'_geo.out'
               iop = 20
               open(iop, file = iopfile, status = 'unknown',
     &                   position = 'append' )

               write(iop,'(/"Errors of cell definition",
     &                      " in EPS Page No. =",
     &               i4)') jnum
               write(iop,'(
     &         "Overlapped Cell IDs  x, y, z  coodinates"/
     &         "(Cells 0       0  indicate undefined region)")')

            end if

            do ii = 1, inne

               read(ioe) ieid, i, j, k, iicl01, iicl02

               call gval(xc,i,j,k,icod,nx,ny,nz,vx,vy,vz)
                  xv = xc(icd1(icod))
                  yv = xc(icd2(icod))

               write(iot,'(1p2e14.6,0p3f7.3)')
     &               xv, yv,
     &               ercol(ieid,1), ercol(ieid,2), ercol(ieid,3)

               if( igsher .gt. 1 ) write(iop,'(2i8,3x,1p3e14.6)')
     &                             iicl01, iicl02, xc(1), xc(2), xc(3)

            end do

            if( igsher .gt. 1 ) close(iop)

         end if
         end if

*-----------------------------------------------------------------------

               goto 900

*-----------------------------------------------------------------------
*        start new region
*-----------------------------------------------------------------------

  200    continue

               rewind ioh

               iblz0 = iblz
               nmed0 = nmed

               i0 = i
               j0 = j
               k0 = k

               ip = 0

            if( ilev1 .gt. 1 .and. ilat1(2,1) .ne. 0 ) then

               ilatc = 1

               ilatx = ilat1(3,1)
               ilaty = ilat1(4,1)
               ilatz = ilat1(5,1)

* added by nais 2019/11/19
cKN 2020/04/08 deleted  ???
            elseif( ilev1 .eq. 2 ) then
               ilatc = 2
* added by nais 2019/11/19

            else

               ilatc = 0

            end if

*-----------------------------------------------------------------------
*        region values : ior = 1
*-----------------------------------------------------------------------

         if( ior .eq. 1 ) then

               mm = 0

            do ir = 1, nr

               call tregck(iblz1,ilev1,ilat1,mr,kr,mm,icr)

               if( icr .ne. 0 ) then

                  ir1 = ir
                  goto 400

               end if

            end do

                  ir1 = 0

  400       continue

*-----------------------------------------------------------------------
*     region value is scaled as smin - 1.0 ( blue to red )
*-----------------------------------------------------------------------

               if( ir1 .ne. 0 ) then

                  rval = vl(ir1) * ( 1.0 - smin ) + smin

               else

                  rval = -1.0d0

               end if

*-----------------------------------------------------------------------
*        region values : ior = 2
*-----------------------------------------------------------------------

         elseif( ior .eq. 2 ) then

          ic=idgr(iblz1)
          itet=kkk0-10000
          itet0=itet
          ihelem=iii0

          call ttetck(ic,nr,itet,itet0,ihelem,icr)
          
          if(icr.gt.0.and.ihelem.gt.0)then

           ir1=ihelem

          else

           ir1=0

          endif

*-----------------------------------------------------------------------
*     region value is scaled as smin - 1.0 ( blue to red )
*-----------------------------------------------------------------------

               if( ir1 .ne. 0 ) then

                  rval = vl(ir1) * ( 1.0 - smin ) + smin

               else

                  rval = -1.0d0

               end if

         end if

*-----------------------------------------------------------------------
*           ic = 5 : first left check
*-----------------------------------------------------------------------

            ic = 5

                  ii = i + ics(ic,1)
                  jj = j + ics(ic,2)
                  kk = k + ics(ic,3)

            if( ii .gt. 0 ) then

                  call gval(xc,ii,jj,kk,icod,nx,ny,nz,vx,vy,vz)

                  u  = rcd(ic,1)
                  v  = rcd(ic,2)
                  w  = rcd(ic,3)

                  mark  = 1
                  markp = 0
                  iblz  = iblz0

                  call tetra0iii !FURUTA20160701 Bugfix

               if(iout.eq.10) then ! bitmap style, T.Sato 2019/02/20
                  mark = 0
               else

                  ici = -1
                if( nowgshow .eq. 1 )
     &          call gomsort(x(1),x(2),x(3),u,v,w,
     &                      nmed,iblz,mark,markp,ici,mtrns)
                call gomprpt(0,x(1),x(2),x(3),xc(1),xc(2),xc(3),u,v,w,
     &                      nmed,iblz,mark,markp,mtrns)
               endif

               if( mark .eq. 1 ) then

                  imat(i,j,k) = imat(ii,jj,kk)
                  goto 100

               else if( mark .eq. 0 .or. mark .eq. -1 .or.
     &             mark .eq. 2 ) then

                  x0 = xc(icd1(icod))
                  y0 = xc(icd2(icod))

               else

                  imat(i,j,k) = -90000
                  goto 100

               end if

            else

                  x0 = x(icd1(icod))
                  y0 = x(icd2(icod))

            end if

*-----------------------------------------------------------------------
*           ic = 6 : first left-down check
*-----------------------------------------------------------------------

            ic = 6

                  ii = i + ics(ic,1)
                  jj = j + ics(ic,2)
                  kk = k + ics(ic,3)

            if( ii .gt. 0 .and. jj .gt. 0 ) then

                  call gval(xc,ii,jj,kk,icod,nx,ny,nz,vx,vy,vz)

                  u  = rcd(ic,1)
                  v  = rcd(ic,2)
                  w  = rcd(ic,3)

                  mark  = 1
                  markp = 0
                  iblz  = iblz0

                  call tetra0iii !FURUTA20160701 Bugfix

               if(iout.eq.10) then ! bitmap style, T.Sato 2019/02/20
                  mark = 0
               else

                  ici = -1
                if( nowgshow .eq. 1 )
     &          call gomsort(x(1),x(2),x(3),u,v,w,
     &                      nmed,iblz,mark,markp,ici,mtrns)
                call gomprpt(0,x(1),x(2),x(3),xc(1),xc(2),xc(3),u,v,w,
     &                      nmed,iblz,mark,markp,mtrns)
               endif

               if( mark .eq. 0 .or. mark .eq. -1 .or.
     &             mark .eq. 2 ) then

                  x1 = xc(icd1(icod))
                  y1 = xc(icd2(icod))

                  ip = ip + 1
                  write(ioh) x1, y1

               end if

            end if

*-----------------------------------------------------------------------
*           ic = 7 : first down check
*-----------------------------------------------------------------------

            ic = 7

                  ii = i + ics(ic,1)
                  jj = j + ics(ic,2)
                  kk = k + ics(ic,3)

            if( jj .gt. 0 ) then

                  call gval(xc,ii,jj,kk,icod,nx,ny,nz,vx,vy,vz)

                  u  = rcd(ic,1)
                  v  = rcd(ic,2)
                  w  = rcd(ic,3)

                  mark  = 1
                  markp = 0
                  iblz  = iblz0

                  call tetra0iii !FURUTA20160701 Bugfix

               if(iout.eq.10) then ! bitmap style, T.Sato 2019/02/20
                  mark = 0
               else

                  ici = -1
                if( nowgshow .eq. 1 )
     &          call gomsort(x(1),x(2),x(3),u,v,w,
     &                      nmed,iblz,mark,markp,ici,mtrns)
                call gomprpt(0,x(1),x(2),x(3),xc(1),xc(2),xc(3),u,v,w,
     &                      nmed,iblz,mark,markp,mtrns)
               endif

               if( mark .eq. 1 ) then

                  imat(i,j,k) = imat(ii,jj,kk)
                  goto 100

               else if( mark .eq. 0 .or. mark .eq. -1 .or.
     &             mark .eq. 2 ) then

                  x2 = xc(icd1(icod))
                  y2 = xc(icd2(icod))

                  ip = ip + 1
                  write(ioh) x2, y2

               else

                  imat(i,j,k) = -90000
                  goto 100

               end if

            else

                  x2 = x(icd1(icod))
                  y2 = x(icd2(icod))

                  ip = ip + 1
                  write(ioh) x2, y2

            end if

*-----------------------------------------------------------------------
*        write new region
*-----------------------------------------------------------------------

                  inum = inum + 1

               if( nowgshow .ne. 1 ) then
                  inrg = inum * 10000 + iblz0
               else if( nowgshow .eq. 1 ) then
                  inrg = inum * 10000 + nmed0
               end if

               imat(i,j,k) = inrg

               ic = 8

*-----------------------------------------------------------------------
*        eight directions movement and final write information
*-----------------------------------------------------------------------

  300    continue

*-----------------------------------------------------------------------
*           iout = 1 : boundary
*                = 2 : boundary + material
*                = 3 : boundary + material number
*                = 4 : boundary + material + material number
*                = 5 : boundary + region number
*                = 6 : boundary + material + region number
*                = 7 : boundary + lattice number
*                = 8 : boundary + material + lattice number
*                = 9 : no geometry
*               = 10 : bitmap style

*-----------------------------------------------------------------------

         if( i .eq. i0 .and. j .eq. j0 .and. k .eq. k0 .and.
     &     ( ic .eq. 5 .or. ic .eq. 6 .or. ic .eq. 7 ) ) then

                  ip = ip + 1
                  write(ioh) x0, y0

            if( iout .eq. 3 .or. iout .eq. 4 .or.
     &          iout .eq. 5 .or. iout .eq. 6 .or.
     &          iout .eq. 7 .or. iout .eq. 8 ) then

                     rewind ioh

                     xpmn =  1.e33
                     xpma = -1.e33
                     ypmn =  1.e33
                     ypma = -1.e33

                  do l = 1, ip

                     read(ioh) xp, yp

                     if( xp .gt. xpma ) xpma = xp
                     if( xp .lt. xpmn ) xpmn = xp
                     if( yp .gt. ypma ) ypma = yp
                     if( yp .lt. ypmn ) ypmn = yp

                  end do

                     xp = ( xpma + xpmn ) / 2.0
                     yp = ( ypma + ypmn ) / 2.0

                     ixys = 2

*-----------------------------------------------------------------------

                     zp = vz(k0)

                     call gxval(xc,icod,xp,yp,zp)

                        ii = i0
                        jj = j0

                        u  = rcd(ic,1)
                        v  = rcd(ic,2)
                        w  = rcd(ic,3)

                        mark  =  1
                        markp =  0
                        ici   = -1

                        call gomsort(xc(1),xc(2),xc(3),u,v,w,
     &                              nmed,iblz,mark,markp,ici,mtrns)

                  if( mark .gt. -2 .and.
     &              ( nmed .ne. nmed0 .or. iblz .ne. iblz0 ) ) then

                        call gval(xc,ii,jj,k0,icod,nx,ny,nz,vx,vy,vz)

                        xp = xc(icd1(icod))
                        yp = xc(icd2(icod))
                        ixys = 1

  120                continue

                        ii = ii + 1
                        jj = jj + 1

                        mark  =  1
                        markp =  0
                        ici   = -1

                        call gomsort(xc(1),xc(2),xc(3),u,v,w,
     &                               nmed,iblz,mark,markp,ici,mtrns)

                     if( mark .gt. -2 .and.
     &                       nmed .eq. nmed0 .and. iblz .eq. iblz0 .and.
     &                       ii.le.nx .and. jj.le.ny) then

                        call gval(xc,ii,jj,k0,icod,nx,ny,nz,vx,vy,vz)

                        xp = xc(icd1(icod))
                        yp = xc(icd2(icod))

                        goto 120

                     else

                        ii = ( ii - 1 + i0 ) / 2
                        jj = ( jj - 1 + j0 ) / 2

                        call gval(xc,ii,jj,k0,icod,nx,ny,nz,vx,vy,vz)

                        xp = xc(icd1(icod))
                        yp = xc(icd2(icod))
                        ixys = 2

                     end if

                  end if

*-----------------------------------------------------------------------

                     sizm = dmhsb(nmed0,4) * 0.7
                     sizr = dmhsg(idgr(iblz0)) * 0.7

               if( iout .eq. 3 .or. iout .eq. 4 ) then

                  if( sizm .gt. 0.0 ) then

                     if4   = nmtnm(nmed0)
                     huni4 = dmtnm(nmed0)

                     write(iot,'(/"w:",a,"/ x(",e16.7,
     &               ") y(",e16.7,") ix(",i1,
     &               ") iy(",i1,") s(",
     &               e13.4,")")')
     &               huni4(1:if4), xp, yp, ixys, ixys, sizm

                  end if

               else if(   iout .eq. 5 .or. iout .eq. 6 .or.
     &                ( ( iout .eq. 7 .or. iout .eq. 8 ) .and.
     &                    ilatc .eq. 0 ) ) then

                  if( sizr .gt. 0.0 ) then

                     if4   = nmtng(idgr(iblz0))
                     huni4 = dmtng(idgr(iblz0))

                     write(iot,'(/"w:",a,"/ x(",e16.7,
     &               ") y(",e16.7,") ix(",i1,
     &               ") iy(",i1,") s(",
     &               e13.4,")")')
     &               huni4(1:if4), xp, yp, ixys, ixys, sizr

                  end if

               else if( ( iout .eq. 7 .or. iout .eq. 8 ) .and.
     &                    ilatc .ne. 0 ) then

                  if( sizr .gt. 0.0 ) then

                     write(nlat,'(a1,i4,a1,i4,a1,i4,a1)')
     &               '(', ilatx, ',', ilaty, ',', ilatz, ')'

                        jk = 0

                     do ik = 1, 20

                        if( nlat(ik:ik) .ne. ' ') then

                           jk = jk + 1
                           mlat(jk:jk) = nlat(ik:ik)

                        end if

                     end do

                     write(iot,'(/"w:",a,"/ x(",e13.4,
     &               ") y(",e16.7,") ix(",i1,
     &               ") iy(",i1,") s(",
     &               e13.4,")")')
     &                   mlat(1:jk), xp, yp, ixys, ixys, sizr

                  end if

               end if

            end if

*-----------------------------------------------------------------------
*           boundary is in it or not
*-----------------------------------------------------------------------

               ibond = 0

                  rewind ioh

               do l = 1, ip

                  read(ioh) x0, y0

                     if( abs( x0 - xvmin ) .le. epsb .or.
     &                   abs( x0 - xvmax ) .le. epsb .or.
     &                   abs( y0 - yvmin ) .le. epsb .or.
     &                   abs( y0 - yvmax ) .le. epsb ) then

                        ibond = 1
                        goto 230

                     end if

               end do

  230          continue

*-----------------------------------------------------------------------

               if( ilatc .eq. 0 .or. iglat .ne. 0 ) then

                  huni11 = huni8

               else

                  huni11 = huni9

               end if

               if(iout.eq.10) huni11 = huni9 ! Always no line for bitmap style, T.Sato 2019/02/20

*-----------------------------------------------------------------------

            if( ior .eq. 0 .and.
     &             ( iout .eq. 1 .or.
     &               iout .eq. 3 .or.
     &               iout .eq. 5 .or.
     &               iout .eq. 7 ) ) then

                  write(iot,'(/a)') huni11//ct(itt)

            else if( ior .eq. 0 .and.
     &             ( iout .eq. 2 .or.
     &               iout .eq. 4 .or.
     &               iout .eq. 6 .or.
     &               iout .eq. 8 .or.
     &               iout .eq.10) ) then  ! T.Sato 2019/02/20

                     if4   = nmtnm(nmed0)
                     huni4 = dmtnm(nmed0)
                     if6   = nmtcl(nmed0)+3
                     huni6(1:2) = 'c['
                     huni6(2+1:2+nmtcl(nmed0))
     &                   = dmtcl(nmed0)(1:nmtcl(nmed0))
                     huni6(3+nmtcl(nmed0):3+nmtcl(nmed0)) = ']'

               if( nmed0 .gt. 0 ) then

                  if( itcl(nmed0) .eq. 0 ) then

                     write(iot,'(/a)') huni7//
     &                    huni3//huni4(1:if4)//huni5//huni6(1:if6)

                     x1 = -1.0
                     y1 = -1.0

                     write(iot,'(1p2e14.6)') x1, y1

                     itcl(nmed0) = 1

                  end if

               else if( nmed0 .eq. 0 ) then

                  if( ivodi .eq. 0 ) then

                     write(iot,'(/a)') huni7//
     &                    huni3//huni4(1:if4)//huni5//huni6(1:if6)

                     x1 = -1.0
                     y1 = -1.0

                     write(iot,'(1p2e14.6)') x1, y1

                     ivodi = ivodi + 1

                  end if

               else if( nmed0 .eq. -100 ) then

                  if( ivodo .eq. 0 ) then

                     write(iot,'(/a)') huni7//
     &                    huni3//huni4(1:if4)//huni5//huni6(1:if6)

                     x1 = -1.0
                     y1 = -1.0

                     write(iot,'(1p2e14.6)') x1, y1

                     ivodo = ivodo + 1

                  end if

               end if

                     if6   = nmtcl(nmed0)+4
                     huni6(1:3) = 'cb('
                     huni6(3+1:3+nmtcl(nmed0))
     &                   = dmtcl(nmed0)(1:nmtcl(nmed0))
                     huni6(4+nmtcl(nmed0):4+nmtcl(nmed0)) = ')'

                  write(iot,'(/a)') huni11//ct(itt)//' '//
     &                              huni6(1:if6)

            else if( ior .ge. 1 .and.
     &             ( iout .eq. 1 .or.
     &               iout .eq. 3 .or.
     &               iout .eq. 5 .or.
     &               iout .eq. 7 ) ) then

                     write(coln,'(f8.5)') rval

                  write(iot,'(/a)') huni11//ct(itt)//' '//
     &                             'cb('//coln//')'

            end if

*-----------------------------------------------------------------------
*           frame:
*-----------------------------------------------------------------------

               if( ibond .eq. 1 ) then

                  write(iot,'("frame:")')
                  write(iot,'(1p2e15.7)') xminf, yminf
                  write(iot,'(1p2e15.7)') xmaxf, yminf
                  write(iot,'(1p2e15.7)') xmaxf, ymaxf
                  write(iot,'(1p2e15.7)') xminf, ymaxf

               end if

*-----------------------------------------------------------------------
*           clip:
*-----------------------------------------------------------------------

               write(iot,'("clip:")')

                  ipl = 2

               call ppline(ipl,ioh,iot,ip,ddel)

*-----------------------------------------------------------------------

               goto 100

         end if

*-----------------------------------------------------------------------
*        four directions movement
*-----------------------------------------------------------------------

  500    continue

               if( ic .eq. 1 ) goto 501
               if( ic .eq. 2 ) goto 502
               if( ic .eq. 3 ) goto 503
               if( ic .eq. 4 ) goto 504
               if( ic .eq. 5 ) goto 505
               if( ic .eq. 6 ) goto 506
               if( ic .eq. 7 ) goto 507
               if( ic .eq. 8 ) goto 508

*-----------------------------------------------------------------------
*        ic = 1
*-----------------------------------------------------------------------

  501    ic = 1

                        ii = i + ics(ic,1)
                        jj = j + ics(ic,2)
                        kk = k + ics(ic,3)

            if( ii .le. nx ) then

                     call gval(xc,ii,jj,kk,icod,nx,ny,nz,vx,vy,vz)

                        u  = rcd(ic,1)
                        v  = rcd(ic,2)
                        w  = rcd(ic,3)

                        mark  = 1
                        markp = 0
                        iblz  = iblz0

                        call tetra0iii !FURUTA20160701 Bugfix

                     if(iout.eq.10) then ! bitmap style, T.Sato 2019/02/20
                        mark = 0
                     else

                        ici = -1
                      if( nowgshow .eq. 1 )
     &                call gomsort(x(1),x(2),x(3),u,v,w,
     &                            nmed,iblz,mark,markp,ici,mtrns)
                      call gomprpt(0,x(1),x(2),x(3),xc(1),xc(2),xc(3),
     &                            u,v,w,nmed,iblz,mark,markp,mtrns)
                     endif

               if( mark .eq. 1 ) then

                  if( imat(ii,jj,kk) .lt. 0 ) goto 110

                  if( j + 1 .le. ny ) then

                        is = 3

                        xd(1) = x(1) + dlg * rcd(is,1)
                        xd(2) = x(2) + dlg * rcd(is,2)
                        xd(3) = x(3) + dlg * rcd(is,3)

                        u  = rcd(is,1)
                        v  = rcd(is,2)
                        w  = rcd(is,3)

                        mark  = 1
                        markp = 0
                        iblz  = iblz0

                        call tetra0iii !FURUTA20160701 Bugfix

                     if(iout.eq.10) then ! bitmap style, T.Sato 2019/02/20
                        mark = 0
                     else

                        ici = -1
                      if( nowgshow .eq. 1 )
     &                call gomsort(x(1),x(2),x(3),u,v,w,
     &                            nmed,iblz,mark,markp,ici,mtrns)
                      call gomprpt(0,x(1),x(2),x(3),xd(1),xd(2),xd(3),
     &                            u,v,w,nmed,iblz,mark,markp,mtrns)
                     endif

                     if( mark .eq. 0 .or. mark .eq. -1 .or.
     &                   mark .eq. 2 ) then

                        do jl = j + 1, ny

                           if( vy(jl) .gt. xd(icd2(icod)) ) goto 601

                           if( imat(i,jl,k) .gt. 0 .and.
     &                         imat(i,jl,k) .ne. inrg ) then

                              imat(i,j,k) = -90000
                              goto 110

                           else

                              imat(i,jl,k) = inrg

                           end if

                        end do

  601                   continue

                     else

                        imat(i,j,k) = -90000
                        goto 110

                     end if


                  end if

                        i = ii
                        j = jj
                        k = kk

                        imat(i,j,k) = inrg

                        x(1) = xc(1)
                        x(2) = xc(2)
                        x(3) = xc(3)

                        ic = 7
                        goto 300

               else if( mark .eq. 0 .or. mark .eq. -1 .or.
     &                  mark .eq. 2 ) then

                        ic = 2
                        goto 310

               else

                        imat(i,j,k) = -90000
                        goto 110

               end if

            else

                        xc(1) = x(1)
                        xc(2) = x(2)
                        xc(3) = x(3)

                        ic = 3
                        goto 310

            end if

*-----------------------------------------------------------------------
*        ic = 2
*-----------------------------------------------------------------------

  502    ic = 2

                        ii = i + ics(ic,1)
                        jj = j + ics(ic,2)
                        kk = k + ics(ic,3)

            if( ii .le. nx .and. jj .le. ny ) then

                     call gval(xc,ii,jj,kk,icod,nx,ny,nz,vx,vy,vz)

                        u  = rcd(ic,1)
                        v  = rcd(ic,2)
                        w  = rcd(ic,3)

                        mark  = 1
                        markp = 0
                        iblz  = iblz0

                        call tetra0iii !FURUTA20160701 Bugfix

                     if(iout.eq.10) then ! bitmap style, T.Sato 2019/02/20
                        mark = 0
                     else

                        ici = -1
                      if( nowgshow .eq. 1 )
     &                call gomsort(x(1),x(2),x(3),u,v,w,
     &                            nmed,iblz,mark,markp,ici,mtrns)
                      call gomprpt(0,x(1),x(2),x(3),xc(1),xc(2),xc(3),
     &                            u,v,w,nmed,iblz,mark,markp,mtrns)
                     endif

               if( mark .eq. 0 .or. mark .eq. -1 .or.
     &             mark .eq. 2 ) then

                        ic = 3
                        goto 310

               end if

            end if

*-----------------------------------------------------------------------
*        ic = 3
*-----------------------------------------------------------------------

  503    ic = 3

                        ii = i + ics(ic,1)
                        jj = j + ics(ic,2)
                        kk = k + ics(ic,3)

            if( jj .le. ny ) then

                     call gval(xc,ii,jj,kk,icod,nx,ny,nz,vx,vy,vz)

                        u  = rcd(ic,1)
                        v  = rcd(ic,2)
                        w  = rcd(ic,3)

                        mark  = 1
                        markp = 0
                        iblz  = iblz0

                        call tetra0iii !FURUTA20160701 Bugfix

                     if(iout.eq.10) then ! bitmap style, T.Sato 2019/02/20
                        mark = 0
                     else

                        ici = -1
                      if( nowgshow .eq. 1 )
     &                call gomsort(x(1),x(2),x(3),u,v,w,
     &                            nmed,iblz,mark,markp,ici,mtrns)
                      call gomprpt(0,x(1),x(2),x(3),xc(1),xc(2),xc(3),
     &                            u,v,w,nmed,iblz,mark,markp,mtrns)
                     endif

               if( mark .eq. 1 ) then

                  if( imat(ii,jj,kk) .lt. 0 ) goto 110

                  if( i - 1 .ge. 1 ) then

                        is = 5

                        xd(1) = x(1) + dlg * rcd(is,1)
                        xd(2) = x(2) + dlg * rcd(is,2)
                        xd(3) = x(3) + dlg * rcd(is,3)

                        u  = rcd(is,1)
                        v  = rcd(is,2)
                        w  = rcd(is,3)

                        mark  = 1
                        markp = 0
                        iblz  = iblz0

                        call tetra0iii !FURUTA20160701 Bugfix

                     if(iout.eq.10) then ! bitmap style, T.Sato 2019/02/20
                        mark = 0
                     else

                        ici = -1
                      if( nowgshow .eq. 1 )
     &                call gomsort(x(1),x(2),x(3),u,v,w,
     &                            nmed,iblz,mark,markp,ici,mtrns)
                      call gomprpt(0,x(1),x(2),x(3),xd(1),xd(2),xd(3),
     &                            u,v,w,nmed,iblz,mark,markp,mtrns)
                     endif

                     if( mark .eq. 0 .or. mark .eq. -1 .or.
     &                   mark .eq. 2 ) then

                        do il = i - 1, 1, -1

                           if( vx(il) .lt. xd(icd1(icod)) ) goto 603

                           if( imat(il,j,k) .gt. 0 .and.
     &                         imat(il,j,k) .ne. inrg ) then

                              imat(i,j,k) = -90000
                              goto 110

                           else

                              imat(il,j,k) = inrg

                           end if

                        end do

  603                   continue

                     else

                        imat(i,j,k) = -90000
                        goto 110

                     end if

                  end if

                        i = ii
                        j = jj
                        k = kk

                        imat(i,j,k) = inrg

                        x(1) = xc(1)
                        x(2) = xc(2)
                        x(3) = xc(3)

                        ic = 1
                        goto 300

               else if( mark .eq. 0 .or. mark .eq. -1 .or.
     &                  mark .eq. 2 ) then

                        ic = 4
                        goto 310

               else

                        imat(i,j,k) = -90000
                        goto 110

               end if

            else

                        xc(1) = x(1)
                        xc(2) = x(2)
                        xc(3) = x(3)

                        ic = 5
                        goto 310

            end if

*-----------------------------------------------------------------------
*        ic = 4
*-----------------------------------------------------------------------

  504    ic = 4

                        ii = i + ics(ic,1)
                        jj = j + ics(ic,2)
                        kk = k + ics(ic,3)

            if( ii .gt. 0 .and. jj .le. ny ) then

                     call gval(xc,ii,jj,kk,icod,nx,ny,nz,vx,vy,vz)

                        u  = rcd(ic,1)
                        v  = rcd(ic,2)
                        w  = rcd(ic,3)

                        mark  = 1
                        markp = 0
                        iblz  = iblz0

                        call tetra0iii !FURUTA20160701 Bugfix

                     if(iout.eq.10) then ! bitmap style, T.Sato 2019/02/20
                        mark = 0
                     else

                        ici = -1
                      if( nowgshow .eq. 1 )
     &                call gomsort(x(1),x(2),x(3),u,v,w,
     &                            nmed,iblz,mark,markp,ici,mtrns)
                      call gomprpt(0,x(1),x(2),x(3),xc(1),xc(2),xc(3),
     &                            u,v,w,nmed,iblz,mark,markp,mtrns)
                     endif

               if( mark .eq. 0 .or. mark .eq. -1 .or.
     &             mark .eq. 2 ) then

                        ic = 5
                        goto 310

               end if

            end if

*-----------------------------------------------------------------------
*        ic = 5
*-----------------------------------------------------------------------

  505    ic = 5

                        ii = i + ics(ic,1)
                        jj = j + ics(ic,2)
                        kk = k + ics(ic,3)

            if( ii .gt. 0 ) then

                     call gval(xc,ii,jj,kk,icod,nx,ny,nz,vx,vy,vz)

                        u  = rcd(ic,1)
                        v  = rcd(ic,2)
                        w  = rcd(ic,3)

                        mark  = 1
                        markp = 0
                        iblz  = iblz0

                        call tetra0iii !FURUTA20160701 Bugfix

                     if(iout.eq.10) then ! bitmap style, T.Sato 2019/02/20
                        mark = 0
                     else

                        ici = -1
                      if( nowgshow .eq. 1 )
     &                call gomsort(x(1),x(2),x(3),u,v,w,
     &                            nmed,iblz,mark,markp,ici,mtrns)
                      call gomprpt(0,x(1),x(2),x(3),xc(1),xc(2),xc(3),
     &                            u,v,w,nmed,iblz,mark,markp,mtrns)
                     endif

               if( mark .eq. 1 ) then

                  if( imat(ii,jj,kk) .lt. 0 ) goto 110

                  if( j - 1 .ge. 1 ) then

                        is = 7

                        xd(1) = x(1) + dlg * rcd(is,1)
                        xd(2) = x(2) + dlg * rcd(is,2)
                        xd(3) = x(3) + dlg * rcd(is,3)

                        u  = rcd(is,1)
                        v  = rcd(is,2)
                        w  = rcd(is,3)

                        mark  = 1
                        markp = 0
                        iblz  = iblz0

                        call tetra0iii !FURUTA20160701 Bugfix

                     if(iout.eq.10) then ! bitmap style, T.Sato 2019/02/20
                        mark = 0
                     else

                        ici = -1
                      if( nowgshow .eq. 1 )
     &                call gomsort(x(1),x(2),x(3),u,v,w,
     &                            nmed,iblz,mark,markp,ici,mtrns)
                      call gomprpt(0,x(1),x(2),x(3),xd(1),xd(2),xd(3),
     &                            u,v,w,nmed,iblz,mark,markp,mtrns)
                     endif

                     if( mark .eq. 0 .or. mark .eq. -1 .or.
     &                   mark .eq. 2 ) then

                        do jl = j - 1, 1, -1

                           if( vy(jl) .lt. xd(icd2(icod)) ) goto 605

                           if( imat(i,jl,k) .gt. 0 .and.
     &                         imat(i,jl,k) .ne. inrg ) then

                              imat(i,j,k) = -90000
                              goto 110

                           else

                              imat(i,jl,k) = inrg

                           end if

                        end do

  605                   continue

                     else

                        imat(i,j,k) = -90000
                        goto 110

                     end if

                  end if

                        i = ii
                        j = jj
                        k = kk

                        imat(i,j,k) = inrg

                        x(1) = xc(1)
                        x(2) = xc(2)
                        x(3) = xc(3)

                        ic = 3
                        goto 300

               else if( mark .eq. 0 .or. mark .eq. -1 .or.
     &                  mark .eq. 2 ) then

                        ic = 6
                        goto 310

               else

                        imat(i,j,k) = -90000
                        goto 110

               end if

            else

                        xc(1) = x(1)
                        xc(2) = x(2)
                        xc(3) = x(3)

                        ic = 7
                        goto 310

            end if

*-----------------------------------------------------------------------
*        ic = 6
*-----------------------------------------------------------------------

  506    ic = 6

                        ii = i + ics(ic,1)
                        jj = j + ics(ic,2)
                        kk = k + ics(ic,3)

            if( ii .gt. 0 .and. jj .gt. 0 ) then

                     call gval(xc,ii,jj,kk,icod,nx,ny,nz,vx,vy,vz)

                        u  = rcd(ic,1)
                        v  = rcd(ic,2)
                        w  = rcd(ic,3)

                        mark  = 1
                        markp = 0
                        iblz  = iblz0

                        call tetra0iii !FURUTA20160701 Bugfix

                     if(iout.eq.10) then ! bitmap style, T.Sato 2019/02/20
                        mark = 0
                     else

                        ici = -1
                      if( nowgshow .eq. 1 )
     &                call gomsort(x(1),x(2),x(3),u,v,w,
     &                            nmed,iblz,mark,markp,ici,mtrns)
                      call gomprpt(0,x(1),x(2),x(3),xc(1),xc(2),xc(3),
     &                            u,v,w,nmed,iblz,mark,markp,mtrns)
                     endif

               if( mark .eq. 0 .or. mark .eq. -1 .or.
     &             mark .eq. 2 ) then

                        ic = 7
                        goto 310

               end if

            end if

*-----------------------------------------------------------------------
*        ic = 7
*-----------------------------------------------------------------------

  507    ic = 7

                        ii = i + ics(ic,1)
                        jj = j + ics(ic,2)
                        kk = k + ics(ic,3)

            if( jj .gt. 0 ) then

                     call gval(xc,ii,jj,kk,icod,nx,ny,nz,vx,vy,vz)

                        u  = rcd(ic,1)
                        v  = rcd(ic,2)
                        w  = rcd(ic,3)

                        mark  = 1
                        markp = 0
                        iblz  = iblz0

                        call tetra0iii !FURUTA20160701 Bugfix

                     if(iout.eq.10) then ! bitmap style, T.Sato 2019/02/20
                        mark = 0
                     else

                        ici = -1
                      if( nowgshow .eq. 1 )
     &                call gomsort(x(1),x(2),x(3),u,v,w,
     &                            nmed,iblz,mark,markp,ici,mtrns)
                      call gomprpt(0,x(1),x(2),x(3),xc(1),xc(2),xc(3),
     &                            u,v,w,nmed,iblz,mark,markp,mtrns)
                     endif

               if( mark .eq. 1 ) then

                  if( imat(ii,jj,kk) .lt. 0 ) goto 110

                  if( i + 1 .le. nx ) then

                        is = 1

                        xd(1) = x(1) + dlg * rcd(is,1)
                        xd(2) = x(2) + dlg * rcd(is,2)
                        xd(3) = x(3) + dlg * rcd(is,3)

                        u  = rcd(is,1)
                        v  = rcd(is,2)
                        w  = rcd(is,3)

                        mark  = 1
                        markp = 0
                        iblz  = iblz0

                        call tetra0iii !FURUTA20160701 Bugfix

                     if(iout.eq.10) then ! bitmap style, T.Sato 2019/02/20
                        mark = 0
                     else

                        ici = -1
                      if( nowgshow .eq. 1 )
     &                call gomsort(x(1),x(2),x(3),u,v,w,
     &                            nmed,iblz,mark,markp,ici,mtrns)
                      call gomprpt(0,x(1),x(2),x(3),xd(1),xd(2),xd(3),
     &                            u,v,w,nmed,iblz,mark,markp,mtrns)
                     endif

                     if( mark .eq. 0 .or. mark .eq. -1 .or.
     &                   mark .eq. 2 ) then

                        do il = i + 1, nx

                           if( vx(il) .gt. xd(icd1(icod)) ) goto 607

                           if( imat(il,j,k) .gt. 0 .and.
     &                         imat(il,j,k) .ne. inrg ) then

                              imat(i,j,k) = -90000
                              goto 110

                           else

                              imat(il,j,k) = inrg

                           end if

                        end do

  607                   continue

                     else

                        imat(i,j,k) = -90000
                        goto 110

                     end if

                  end if

                        i = ii
                        j = jj
                        k = kk

                        imat(i,j,k) = inrg

                        x(1) = xc(1)
                        x(2) = xc(2)
                        x(3) = xc(3)

                        ic = 5
                        goto 300

               else if( mark .eq. 0 .or. mark .eq. -1 .or.
     &                  mark .eq. 2 ) then

                        ic = 8
                        goto 310

               else

                        imat(i,j,k) = -90000
                        goto 110

               end if

            else

                        xc(1) = x(1)
                        xc(2) = x(2)
                        xc(3) = x(3)

                        ic = 1
                        goto 310

            end if

*-----------------------------------------------------------------------
*        ic = 8
*-----------------------------------------------------------------------

  508    ic = 8

                        ii = i + ics(ic,1)
                        jj = j + ics(ic,2)
                        kk = k + ics(ic,3)

            if( ii .le. nx .and. jj .gt. 0 ) then

                     call gval(xc,ii,jj,kk,icod,nx,ny,nz,vx,vy,vz)

                        u  = rcd(ic,1)
                        v  = rcd(ic,2)
                        w  = rcd(ic,3)

                        mark  = 1
                        markp = 0
                        iblz  = iblz0

                        call tetra0iii !FURUTA20160701 Bugfix

                     if(iout.eq.10) then ! bitmap style, T.Sato 2019/02/20
                        mark = 0
                     else

                        ici = -1
                      if( nowgshow .eq. 1 )
     &                call gomsort(x(1),x(2),x(3),u,v,w,
     &                            nmed,iblz,mark,markp,ici,mtrns)
                      call gomprpt(0,x(1),x(2),x(3),xc(1),xc(2),xc(3),
     &                            u,v,w,nmed,iblz,mark,markp,mtrns)
                     endif

               if( mark .eq. 0 .or. mark .eq. -1 .or.
     &             mark .eq. 2 ) then

                        ic = 1
                        goto 310

               end if

            end if

                  goto 501

*-----------------------------------------------------------------------

  310       continue

                  x3 = xc(icd1(icod))
                  y3 = xc(icd2(icod))

                  ip = ip + 1
                  write(ioh) x3, y3

            goto 300

*-----------------------------------------------------------------------

  110 continue

               do k = 1, nz
               do j = 1, ny
               do i = 1, nx

                  if( imat(i,j,k) .eq. inrg ) then

                     imat(i,j,k) = -90000

                  end if

               end do
               end do
               end do

               goto 100

*-----------------------------------------------------------------------

  900 continue

*-----------------------------------------------------------------------
*        restore initial mesh
*-----------------------------------------------------------------------

               do i = 1, nx
                  vx(i) = xmin0 + xdel0 * dble( i - 1 )
               end do

               do i = 1, ny
                  vy(i) = ymin0 + ydel0 * dble( i - 1 )
               end do

                  vz(1) = zmin0

*-----------------------------------------------------------------------

      close( ioh )

      if( igsher .ne. 0 ) then
         close( ioe )
      end if
         igsher = 0

*-----------------------------------------------------------------------
      if(nlat3.gt.0) call tetrafingshow !FURUTA20221209

* added by NAIS 2019/11/19
      nowgshow = 0

      return
      end


************************************************************************
*                                                                      *
      subroutine gxval(x,icod,vx,vy,vz)
*                                                                      *
*       give x,y,z values from vx,vy,vz for gshow                      *
*       last modified by K.Niita on 2003/06/10                         *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

      dimension x(3)

*-----------------------------------------------------------------------

      if( icod .eq. 1 ) then

            x(1) = vx
            x(2) = vy
            x(3) = vz

      else if( icod .eq. 2 ) then

            x(1) = vz
            x(2) = vy
            x(3) = vx

      else if( icod .eq. 3 ) then

            x(1) = vy
            x(2) = vz
            x(3) = vx

      end if

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine gval(x,i,j,k,icod,nx,ny,nz,vx,vy,vz)
*                                                                      *
*       give x,y,z values from i,j,k for gshow                         *
*       last modified by K.Niita on 2001/12/28                         *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      dimension   x(3)
      dimension   vx(nx)
      dimension   vy(ny)
      dimension   vz(nz)

*-----------------------------------------------------------------------

      if( icod .eq. 1 ) then

            x(1) = vx(i)
            x(2) = vy(j)
            x(3) = vz(k)

      else if( icod .eq. 2 ) then

            x(1) = vz(k)
            x(2) = vy(j)
            x(3) = vx(i)

      else if( icod .eq. 3 ) then

            x(1) = vy(j)
            x(2) = vz(k)
            x(3) = vx(i)

      end if

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine pdshow(m,iout,nl,lt,nlb,ltb,
     &                  nr,mr,kr,ns,ms,ks,nx,ny,
     &                  idasa)
*                                                                      *
*       output 3d geometry                                             *
*       last modified by K.Niita on 2002/10/24                         *
*                                                                      *
************************************************************************
      use sangelmod, only: itsans

*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'

      parameter ( pi = 3.141592653589793d0 )

*-----------------------------------------------------------------------

      common /verjam/ versn, lastr, iyeav, imonv, idayv

      common /tall11/ itfln(itlmax), itfll(itlmax,6), ctfln(itlmax,6)
      character ctfln*100
      common /tall14/ ittll(itlmax), ittle(itlmax)
      character ittle*80
      common /tall15/ itanl(itlmax), itang(itlmax)
      character itang*200

      common /tall22/ itmcn(itlmax), itmtn(itlmax), itmtt(itlmax)
      common /tall27/ itaxl(itlmax), itaxt(itlmax),
     &                itayl(itlmax), itayt(itlmax),
     &                itazl(itlmax), itazt(itlmax)
      character itaxt*200, itayt*200, itazt*200

      common /tall35/ iteps(itlmax)
      common /tall38/ itres(itlmax)
      common /tall39/ rtwid(itlmax)

      common /tall40/ rtorg(itlmax,3), rteye(itlmax,3), rtlit(itlmax,3),
     &                rtwin(itlmax,3), rtbrt(itlmax,2), ithvn(itlmax),
     &                itwin(itlmax,2), itbox(itlmax), itmir(itlmax),
     &                rtbox(itlmax,5,10), rtout(itlmax), rthet(itlmax),
     &                itlin(itlmax), itshd(itlmax), itgxs(itlmax)

      common /tall47/ itbtr(itlmax,5,4), rtbtr(itlmax,5,13)

*-----------------------------------------------------------------------

      common /shadw/  shadd, ishdd
      common /cubdb/  nbnd

      data shadd / 1000.0 /
      data ishdd / 50 /
      data nbnd  / 4 /

*-----------------------------------------------------------------------

      common /tlgeom/ iblz1,iblz2
!$OMP THREADPRIVATE(/tlgeom/)
      common /tlglat/ ilev1,ilev2,ilat1(5,10),ilat2(5,10)
!$OMP THREADPRIVATE(/tlglat/)
      common /tlcost/ costha, uang(3), nsurf
!$OMP THREADPRIVATE(/tlcost/)

      common /kmat1a/ mxmat, mxmat0, mxnel
      common /kmat1d/ idmn(0:kvlmax), idnm(kvmmax)
      common /regdc/  idrg(kvlmax), idgr(kvmmax)

*-----------------------------------------------------------------------

      dimension   lt(nl)
      dimension   ltb(nlb)
      dimension   kr(mr)
      dimension   ks(ms)
      dimension   nmat(nx,ny,2)
      dimension   imat(nx,ny)
      dimension   jmat(nx,ny,2)
      dimension   ibon(nx,ny,4)

*-----------------------------------------------------------------------

      common /tboxs/ xwin1, ywin1, zwin1, xeye(3),
     &               unx, uny, unz, umx, umy, umz

      dimension itcl(kvlmax)
      dimension x(3), xc(3), angb(3), angd(3), ceye(3)
      dimension xbnd(16), ybnd(16)

      dimension arsh(3,3)
      character crsh(3)*1
      character ccsh(3)*4

*-----------------------------------------------------------------------

      data epsb / 1.0d-08 /

*-----------------------------------------------------------------------

      character huni7*10
      data huni7 /'h: x ny1,0'/

      character huni11*16
      character huni8*16
      data huni8 /'hb:   line clip '/
      character huni9*16
      data huni9 /'hb: noline clip '/

      character huni3*10
      data huni3 /'  y2=[y1]('/

      character huni4*80
      character huni6*30

      character huni5*3
      data huni5 /'),i'/

      character nlat*20
      character mlat*20

      character ct(7)*3
      data ct /'zzz','zz ','z  ','   ','t  ',
     &         'tt ','ttt'/

      common /mtnmcl/ dmhsb(-1:kvlmax,4),
     &                dmtnm(-1:kvlmax), dmtcl(-1:kvlmax),
     &                nmtnm(-1:kvlmax), nmtcl(-1:kvlmax)
      character dmtnm*80, dmtcl*30

      dimension org(3), eye(3), cli(3), win(3), nwn(2), brt(2)
      dimension box(5,10)

      dimension un1(3,6)
      data un1 / 1.d0, 0.d0, 0.d0,
     &          -1.d0, 0.d0, 0.d0,
     &           0.d0, 1.d0, 0.d0,
     &           0.d0,-1.d0, 0.d0,
     &           0.d0, 0.d0, 1.d0,
     &           0.d0, 0.d0,-1.d0/

*-----------------------------------------------------------------------

      character fname*100, fnume*3

      character cha*1
      data cha /"'"/

      character rpa*1
      data rpa /'}'/
      character yen*1

      data dps1 / 0.12345678d-05 /
      data dps2 / 0.23456789d-05 /

*-----------------------------------------------------------------------
cFURUTA20150714 TETRA

*-----------------------------------------------------------------------
*     internal function
*-----------------------------------------------------------------------

               vx(ik) = xmin + xdel * dble( ik - 1 )
               vy(ik) = ymin + ydel * dble( ik - 1 )

*-----------------------------------------------------------------------

      yen  = char(92)

*-----------------------------------------------------------------------
*     open output file
*-----------------------------------------------------------------------

            fname = ctfln(m,1)

            iot = 31
            open(iot, file = fname, status = 'unknown' )

*-----------------------------------------------------------------------
*     input echo
*-----------------------------------------------------------------------

               call tdshech(iot,m,0,0)

*-----------------------------------------------------------------------
*     set constants
*-----------------------------------------------------------------------

            org(1) = rtorg(m,1)
            org(2) = rtorg(m,2)
            org(3) = rtorg(m,3)

            win(1) = rtwin(m,1)
            win(2) = rtwin(m,2)
            win(3) = rtwin(m,3)

            eye(1) = rteye(m,1)
            eye(2) = rteye(m,2)
            eye(3) = rteye(m,3)

            cli(1) = rtlit(m,1)
            cli(2) = rtlit(m,2)
            cli(3) = rtlit(m,3)

            nwn(1) = itwin(m,1)
            nwn(2) = itwin(m,2)

            brt(1) = rtbrt(m,1)
            brt(2) = rtbrt(m,2)

            ihvn   = ithvn(m)
            ibox   = itbox(m)

            mirror = itmir(m)
            theta  = rthet(m) * mirror

            nsb    = 0
            if( nlb .gt. 0 .or. ns .gt. 0 ) nsb = 1

*-----------------------------------------------------------------------

            do i = 1, ibox

               do j = 1, 3

                  k = ( j - 1 ) * 3

                  if( itbtr(m,i,4) .eq. 0 ) then

                     box(i,k+1) = rtbox(m,i,k+1)
                     box(i,k+2) = rtbox(m,i,k+2)
                     box(i,k+3) = rtbox(m,i,k+3)

                  else

                     call trnsxv(rtbox(m,i,k+1),
     &                           rtbox(m,i,k+2),
     &                           rtbox(m,i,k+3),
     &                           xxa,yya,zza,itbtr(m,i,4))

                     box(i,k+1) = xxa
                     box(i,k+2) = yya
                     box(i,k+3) = zza

                  end if

               end do

                  box(i,10) = rtbox(m,i,10)

            end do

            call intbox(ibox,box,ierr)

*-----------------------------------------------------------------------
*     header and definitions for angel
*-----------------------------------------------------------------------

               write(iot,'(/"#",78("-"))')

               write(iot,'("msuc: {",80a1)')
     &                  ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

               write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
               write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

               write(iot,'(/"p: noxt noyt itic(-1)")')

               if( itaxl(m) .ne. 0 ) then

                  write(iot,'(/"p: xtxt(1)")')
                  write(iot,'(/"x: ",200a1)')
     &                       (itaxt(m)(i:i),i=1,itaxl(m))

               end if

               if( itayl(m) .ne. 0 ) then

                  write(iot,'(/"p: ytxt(1)")')
                  write(iot,'( "y: ",200a1)')
     &                       (itayt(m)(i:i),i=1,itayl(m))

               end if

                  form  = win(2) / win(1)
                  xfac  = 1.2
                  afac  = 0.5
                  izlog = 1
                  inocm = 1
                  inolg = 1

                  if( itanl(m) .gt. 0 )
     &            call anset(itang(m),itanl(m),
     &                       form,xfac,afac,cmin,cmax,izlog,inocm,inolg)

                  if( form .le. 1.0 ) then

                     scal = form**0.25
                     xfac = xfac / form**0.5
                     xorg = 0.0 - 0.07 * ( 1. - form )**0.25
                     yorg = min( 20.d0, (1.0/form-1.0)/2.5
     &                          - 0.1 / form**0.25 )

                  else

                     scal = 1.0 / form**0.6
                     xfac = xfac / form**0.5
                     xorg = min( 20.d0, ( form - 1.0 ) / 1.5)
                     yorg = 0.0 - 0.08 / form**0.6

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
               write(iot,'("p: legs[c5*1.3]")')

*-----------------------------------------------------------------------

               xmin = 0.0d0
               xmax = win(1)
               ymin = 0.0d0
               ymax = win(2)

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
*     open temporary file
*-----------------------------------------------------------------------

         if( iout .gt. 0 ) then

            ioh = 18
            open(ioh,status='scratch',form='unformatted')

         end if

         if( iout .gt. 1 ) then

            iog = 19
            open(iog,status='scratch',form='unformatted')

            iop = 20
            open(iop,status='scratch',form='unformatted')

            ioq = 22
            open(ioq,status='scratch',form='unformatted')

         end if

*-----------------------------------------------------------------------
*           color of material
*-----------------------------------------------------------------------

            k = 0

         do i = 1, mxmat

            itcl(i) = 0

         end do

            ivod = 0

*-----------------------------------------------------------------------
*           width of line
*-----------------------------------------------------------------------

            ifac = nint( rtwid(m) * 10.0 )

            if( ifac .le. 1 ) then
               itt = 1
            else if( ifac .le. 3 ) then
               itt = 2
            else if( ifac .le. 5 ) then
               itt = 3
            else if( ifac .le. 10 ) then
               itt = 4
            else if( ifac .le. 20 ) then
               itt = 5
            else if( ifac .le. 30 ) then
               itt = 6
            else
               itt = 7
            end if

*-----------------------------------------------------------------------
*           window dimension and resolution
*-----------------------------------------------------------------------

               ires  = itres(m)
               ires2 = ires**2

               xmin0 = 0.0d0
               ymin0 = 0.0d0

               xmax0 = win(1)
               ymax0 = win(2)

               xdel0 = win(1) / dble( nx - 1 )
               ydel0 = win(2) / dble( nx - 1 )

               xdel1 = win(1) / dble( ires )
               ydel1 = win(2) / dble( ires )

               xdeld = win(1) / dble( nx )
               ydeld = win(2) / dble( nx )

               xdel = ( xdel1 + 2.0 * xdeld ) / dble( nx - 1 )
               ydel = ( ydel1 + 2.0 * ydeld ) / dble( ny - 1 )

*-----------------------------------------------------------------------
*           bond resolution
*-----------------------------------------------------------------------

               mbnd = 2**nbnd

               xdelb = xdel / dble( mbnd )
               ydelb = ydel / dble( mbnd )

            do i = 1, mbnd

               xbnd(i) = xdelb / 2.0 + dble( i - 1 ) * xdelb
               ybnd(i) = ydelb / 2.0 + dble( i - 1 ) * ydelb

            end do

               ddel = ( min( xdel, ydel ) / dble( mbnd ) )**2

*-----------------------------------------------------------------------
*           light point ceye(3)
*-----------------------------------------------------------------------

               uwin1 = sin( cli(1)/180.d0*pi ) * cos( cli(2)/180.d0*pi )
               vwin1 = sin( cli(1)/180.d0*pi ) * sin( cli(2)/180.d0*pi )
               wwin1 = cos( cli(1)/180.d0*pi )

               ceye(1) = cli(3) * uwin1 + org(1)
               ceye(2) = cli(3) * vwin1 + org(2)
               ceye(3) = cli(3) * wwin1 + org(3)

*-----------------------------------------------------------------------
*           window direction(uwin0), center(xwin0) and eye(xeye) point
*-----------------------------------------------------------------------

               uwin0 = sin( eye(1)/180.d0*pi ) * cos( eye(2)/180.d0*pi )
               vwin0 = sin( eye(1)/180.d0*pi ) * sin( eye(2)/180.d0*pi )
               wwin0 = cos( eye(1)/180.d0*pi )

               xwin0 = win(3) * uwin0 + org(1)
               ywin0 = win(3) * vwin0 + org(2)
               zwin0 = win(3) * wwin0 + org(3)

               xeye(1) = eye(3) * uwin0 + org(1)
               xeye(2) = eye(3) * vwin0 + org(2)
               xeye(3) = eye(3) * wwin0 + org(3)

*-----------------------------------------------------------------------
*           window origin(xwin1) and vecter(unx, umx)
*-----------------------------------------------------------------------

               unx = un1(2,ihvn) * wwin0 - un1(3,ihvn) * vwin0
               uny = un1(3,ihvn) * uwin0 - un1(1,ihvn) * wwin0
               unz = un1(1,ihvn) * vwin0 - un1(2,ihvn) * uwin0

               unn = sqrt( unx**2 + uny**2 + unz**2 )

               unx =  unx / unn
               uny =  uny / unn
               unz =  unz / unn

               umx =  uny * wwin0 - unz * vwin0
               umy =  unz * uwin0 - unx * wwin0
               umz =  unx * vwin0 - uny * uwin0

               unn = sqrt( umx**2 + umy**2 + umz**2 )

               umx = -umx / unn
               umy = -umy / unn
               umz = -umz / unn

*-----------------------------------------------------------------------
*           rotation of window
*-----------------------------------------------------------------------

               csthe = cos( -theta / 180.d0 * pi )
               snthe = sin( -theta / 180.d0 * pi )

               xwin2 =  csthe * unx + snthe * umx
               ywin2 =  csthe * uny + snthe * umy
               zwin2 =  csthe * unz + snthe * umz

               xwin3 = -snthe * unx + csthe * umx
               ywin3 = -snthe * uny + csthe * umy
               zwin3 = -snthe * unz + csthe * umz

               unx = xwin2 * mirror
               uny = ywin2 * mirror
               unz = zwin2 * mirror

               umx = xwin3
               umy = ywin3
               umz = zwin3

*-----------------------------------------------------------------------

               xwin1 = xwin0 - win(1) / 2.d0 * unx
     &                       - win(2) / 2.d0 * umx
               ywin1 = ywin0 - win(1) / 2.d0 * uny
     &                       - win(2) / 2.d0 * umy
               zwin1 = zwin0 - win(1) / 2.d0 * unz
     &                       - win(2) / 2.d0 * umz

*-----------------------------------------------------------------------
*        for arrows
*-----------------------------------------------------------------------

               vnx = xeye(1) - org(1)
               vny = xeye(2) - org(2)
               vnz = xeye(3) - org(3)

               axm = 0.d0

         do mm = 1, 3

               v0x = 0.0d0
               v0y = 0.0d0
               v0z = 0.0d0

               if( mm .eq. 1 ) v0x = v0x + 10.0
               if( mm .eq. 2 ) v0y = v0y + 10.0
               if( mm .eq. 3 ) v0z = v0z + 10.0

               vpx = vnx - v0x
               vpy = vny - v0y
               vpz = vnz - v0z

               vqx = xwin0 - org(1)
               vqy = ywin0 - org(2)
               vqz = zwin0 - org(3)

               vnq = vnx * ( v0x - vqx )
     &             + vny * ( v0y - vqy )
     &             + vnz * ( v0z - vqz )
               vnp = vnx * vpx + vny * vpy + vnz * vpz

               if( abs(vnp) .gt. 1.d-10 ) then
                  vax = v0x - vnq / vnp * vpx
                  vay = v0y - vnq / vnp * vpy
                  vaz = v0z - vnq / vnp * vpz
               else
                  vax = v0x
                  vay = v0y
                  vaz = v0z
               end if

                  vax = vax + org(1)
                  vay = vay + org(2)
                  vaz = vaz + org(3)

                  bunb = umx * uny - umy * unx
                  bunc = umy * unz - umz * uny
                  buna = umz * unx - umx * unz

               if( abs(bunb) .gt. 1.d-10 ) then

                  xon = ( - umy * ( vax - xwin0 )
     &                    + umx * ( vay - ywin0 ) )
     &                / bunb
                  yon = (   uny * ( vax - xwin0 )
     &                    - unx * ( vay - ywin0 ) )
     &                / bunb

               else if( abs(bunc) .gt. 1.d-10 ) then

                  xon = ( - umz * ( vay - ywin0 )
     &                    + umy * ( vaz - zwin0 ) )
     &                / bunc
                  yon = (   unz * ( vay - ywin0 )
     &                    - uny * ( vaz - zwin0 ) )
     &                / bunc

               else if( abs(buna) .gt. 1.d-10 ) then

                  xon = ( - umx * ( vaz - zwin0 )
     &                    + umz * ( vax - xwin0 ) )
     &                / buna
                  yon = (   unx * ( vaz - zwin0 )
     &                    - unz * ( vax - xwin0 ) )
     &                / buna

               else

                  xon = 0.0d0
                  yon = 0.0d0

               end if

                  axy = sqrt( xon**2 + yon**2 )

                  if( axy .gt. axm ) axm = axy

                  arsh(mm,1) = xon
                  arsh(mm,2) = yon
                  arsh(mm,3) = axy

         end do

cKN 2016/09/12 default size
                  c0 = 10.d0
                  c1 = win(1) / 100.d0 * c0
                  c2 = c1

            if( itgxs(m) .le. 1 ) then
               write(iot,'(/"#----- Axis Arrow Here ------ ")')
               write(iot,'("# x-position")')
               write(iot,'("set: c1[",1p1g14.7,"]")') c1
               write(iot,'("# y-position")')
               write(iot,'("set: c2[",1p1g14.7,"]")') c2
               write(iot,'("# size of arrows")')
               write(iot,'("set: c3[",1p1g14.7,"]")') 1.d0
               write(iot,'("# size of character")')
               write(iot,'("set: c4[",1p1g14.7,"]")') 1.d0
            else if( itgxs(m) .eq. 2 ) then
                  c11 = win(1) / 100.d0 * 50.d0
                  c12 = c11
               write(iot,'(/"#----- Axis Arrow Here ------ ")')
               write(iot,'("# x-position")')
               write(iot,'("set: c1[",1p1g14.7,"]")') c11
               write(iot,'("# y-position")')
               write(iot,'("set: c2[",1p1g14.7,"]")') c12
               write(iot,'("# size of arrows")')
               write(iot,'("set: c3[",1p1g14.7,"]")') 6.d0
               write(iot,'("# size of character")')
               write(iot,'("set: c4[",1p1g14.7,"]")') 1.d0
            end if

                  crsh(1) = 'x'
                  crsh(2) = 'y'
                  crsh(3) = 'z'

                  ccsh(1) = 'c(g)'
                  ccsh(2) = 'c(b)'
                  ccsh(3) = 'c(r)'

         do mm = 1, 3

               ax = arsh(mm,1) / axm * c1 * 0.6667
               ay = arsh(mm,2) / axm * c1 * 0.6667

            if( arsh(mm,1) .ne. 0.0d0 .and. arsh(mm,2) .ne. 0.0d0 ) then

               bx = arsh(mm,1)
     &            / sqrt(arsh(mm,1)**2+arsh(mm,2)**2) * c1 *  0.15
               by = arsh(mm,2)
     &            / sqrt(arsh(mm,1)**2+arsh(mm,2)**2) * c1 *  0.15

            else

               bx = 0.0d0
               by = 0.0d0

            end if

            if( itgxs(m) .ne. 0 ) then

               if( abs(ax) .gt. 1.d-10 .or. abs(ay) .gt. 1.d-10)
     &         write(iot,'("a: x[c1] y[c2] ",
     &                     "ax[c1+(",1p1g14.7,")*c3] ",
     &                     "ay[c2+(",1p1g14.7,")*c3] ",
     &                       a4 )') ax, ay, ccsh(mm)

               write(iot,'("w:",a1,"/ ",
     &                     "x[c1+(",1p1g14.7,")*c3+(",
     &                                1p1g14.7,")*c4] "
     &                     "y[c2+(",1p1g14.7,")*c3+(",
     &                                1p1g14.7,")*c4] "
     &                     "ix(2) iy(2) s[c4] c(e)" )')
     &                     crsh(mm), ax, bx, ay, by

            else

               if( abs(ax) .gt. 1.d-10 .or. abs(ay) .gt. 1.d-10)
     &         write(iot,'("# a: x[c1] y[c2] ",
     &                     "ax[c1+(",1p1g14.7,")*c3] ",
     &                     "ay[c2+(",1p1g14.7,")*c3] ",
     &                       a4 )') ax, ay, ccsh(mm)

               write(iot,'("# w:",a1,"/ ",
     &                     "x[c1+(",1p1g14.7,")*c3+(",
     &                                1p1g14.7,")*c4] "
     &                     "y[c2+(",1p1g14.7,")*c3+(",
     &                                1p1g14.7,")*c4] "
     &                     "ix(2) iy(2) s[c4] c(e)" )')
     &                     crsh(mm), ax, bx, ay, by

            end if

         end do

               write(iot,'("#----- End of Axis Arrow ------ ")')


*-----------------------------------------------------------------------
*     do loop 900 for fine resolution
*-----------------------------------------------------------------------

      do 900 ire = 1, ires2

               ire1 = ire - ( ire - 1 ) / ires * ires
               ire2 = ( ire - 1 ) / ires + 1

               xmin = xmin0 + xdel1 * dble(ire1-1) - xdeld + xdel1*dps1 ! T.Sato 2024/04/23
               ymin = ymin0 + ydel1 * dble(ire2-1) - ydeld + ydel1*dps2 ! T.Sato 2024/04/23

*-----------------------------------------------------------------------
*           min and max of x, y and frame
*-----------------------------------------------------------------------

               xminf = xmin + xdeld - xdeld*dps1 ! T.Sato 2024/04/23
               yminf = ymin + ydeld - ydeld*dps2 ! T.Sato 2024/04/23
               xmaxf = xminf + xdel1
               ymaxf = yminf + ydel1

               xvmin = vx(1)
               yvmin = vy(1)
               xvmax = vx(nx)
               yvmax = vy(ny)

*-----------------------------------------------------------------------
*        initialization of imat and ibon
*-----------------------------------------------------------------------

            do i = 1, nx
            do j = 1, ny

                  ibon(i,j,1) = 0
                  ibon(i,j,2) = 0
                  ibon(i,j,3) = 0
                  ibon(i,j,4) = 0

               call xvwi(vx(i),vy(j),x(1),x(2),x(3),u,v,w)
               call smats(1,iout,m,nsb,nl,lt,nlb,ltb,
     &                    nr,mr,kr,ns,ms,ks,
     &                    x,u,v,w,ceye,nmed,isb,ishad)

               if( nmed .gt. 0 ) then

                  imat(i,j) = 0

               else

                  imat(i,j) = -1

               end if

                  nmat(i,j,1) = ishad
                  jmat(i,j,1) = isb
                  jmat(i,j,2) = nmed

            end do
            end do

*-----------------------------------------------------------------------
*        draft
*-----------------------------------------------------------------------

         if( iout .eq. 0 ) then

                  ncont = 0

               do ii = 1, nx
               do jj = 1, ny

                        nmed0 = jmat(ii,jj,2)

                  if( nmed0 .gt. 0 ) then

                        ncont = ncont + 1

                        if4   = nmtnm(nmed0)
                        huni4 = dmtnm(nmed0)
                        if6   = nmtcl(nmed0)+3
                        huni6(1:2) = 'c['
                        huni6(2+1:2+nmtcl(nmed0))
     &                      = dmtcl(nmed0)(1:nmtcl(nmed0))
                        huni6(3+nmtcl(nmed0):3+nmtcl(nmed0)) = ']'

                     if( itcl(nmed0) .eq. 0 ) then

                        write(iot,'(/a)') huni7//
     &                       huni3//huni4(1:if4)//huni5//huni6(1:if6)

                        x1 = -1.0
                        y1 = -1.0

                        write(iot,'(1p2e14.6)') x1, y1

                        itcl(nmed0) = 1

                     end if

                  end if

               end do
               end do

            if( ncont .gt. 0 ) then

                  write(iot,'(/"hb: ")')
                  write(iot,'( "bmap: width=",e16.7,
     &                         " hight=",e16.7)') xdel, ydel

               do ii = 1, nx
               do jj = 1, ny

                        nmed0 = jmat(ii,jj,2)

                  if( nmed0 .gt. 0 ) then

                        shad = dble( abs( nmat(ii,jj,1) ) ) / shadd

                        col1 = dmhsb(nmed0,1)
                        col2 = dmhsb(nmed0,2)
                        col3 = dmhsb(nmed0,3)

                     if( itcl(nmed0) .eq. 0 ) then

                        if4   = nmtnm(nmed0)
                        huni4 = dmtnm(nmed0)
                        if6   = nmtcl(nmed0)+3
                        huni6(1:2) = 'c['
                        huni6(2+1:2+nmtcl(nmed0))
     &                      = dmtcl(nmed0)(1:nmtcl(nmed0))
                        huni6(3+nmtcl(nmed0):3+nmtcl(nmed0)) = ']'

                        write(iot,'(/a)') huni7//
     &                       huni3//huni4(1:if4)//huni5//huni6(1:if6)

                        x1 = -1.0
                        y1 = -1.0

                        write(iot,'(1p2e14.6)') x1, y1

                        itcl(nmed0) = 1

                     end if

                        colx = col1
                        coly = col2
                        colz = col3

                     call givcol(shad,colx,coly,colz,brt)

                     write(iot,'(1p2e14.6,0p3f7.3)')
     &                  vx(ii), vy(jj), colx, coly, colz

                  end if

               end do
               end do

            end if

            goto 900

         end if

*-----------------------------------------------------------------------
*     start check each cell
*-----------------------------------------------------------------------

               inum = 2

  100 continue

               ic = 1

               do j = 1, ny
               do i = 1, nx

                  if( imat(i,j) .eq. 0 ) goto 200

               end do
               end do

               goto 900

*-----------------------------------------------------------------------
*        start new region
*-----------------------------------------------------------------------

  200    continue

               rewind ioh

               isb0   = jmat(i,j,1)
               nmed0  = jmat(i,j,2)
               ishad0 = nmat(i,j,1)

               ibtm  = 1
               if( isb0 .le. 0 .and. itshd(m) .eq. 0 ) ibtm = 0

               i0 = i
               j0 = j

               ig = 0
               ip = 0

               do jj = 1, ny
               do ii = 1, nx

                  nmat(ii,jj,2) = 0

               end do
               end do

               if( ibtm .gt. 0 .and. iout .ne. 1 ) then

                  ig = ig + 1
                  nmat(i,j,2) = ishad0

               end if

               if4   = nmtnm(nmed0)
               huni4 = dmtnm(nmed0)
               if6   = nmtcl(nmed0)+3
               huni6(1:2) = 'c['
               huni6(2+1:2+nmtcl(nmed0))
     &             = dmtcl(nmed0)(1:nmtcl(nmed0))
               huni6(3+nmtcl(nmed0):3+nmtcl(nmed0)) = ']'

*-----------------------------------------------------------------------
*           ic = 5 : first left check
*-----------------------------------------------------------------------

            ic = 5

                  ii = i - 1
                  jj = j

            if( ii .gt. 0 ) then

               if( isb0  .eq. jmat(ii,jj,1) .and.
     &             nmed0 .eq. jmat(ii,jj,2) ) then

                        imat(i,j) = imat(ii,jj)

                        goto 100

               else

                  if( ibon(ii,jj,1) .gt. 0 ) then

                        x0 = vx(ii) + xbnd( ibon(ii,jj,1) )
                        y0 = vy(jj)

                  else

                        delx = xdel / 2.0
                        ipsm = 1
                        vxx  = vx(ii) + delx
                        vyy  = vy(jj)

                     do k = 1, nbnd

                        call xvwi(vxx,vyy,x(1),x(2),x(3),u,v,w)
                        call smats(0,iout,m,nsb,nl,lt,nlb,ltb,
     &                             nr,mr,kr,ns,ms,ks,
     &                             x,u,v,w,ceye,nmed,isb,ishad)

                        ipsm = 1
                        if( isb .eq. isb0 .and.
     &                      nmed .eq. nmed0 ) ipsm = -1

                        delx = delx / 2.0
                        vxx  = vxx + ipsm * delx

                     end do

                        x0 = vxx
                        y0 = vyy

                        ibon(ii,jj,1) =
     &                  int( ( vxx - vx(ii) ) / delx / 2.0 ) + 1

                  end if

                     if( ibtm .gt. 0 .and. iout .ne. 1 ) then

                        ig = ig + 1
                        nmat(ii,jj,2) = ishad0

                     end if

               end if

            else

                        x0 = vx(i)
                        y0 = vy(j)

            end if

*-----------------------------------------------------------------------
*           ic = 6 : first left-down check
*-----------------------------------------------------------------------

            ic = 6

                  ii = i - 1
                  jj = j - 1

            if( ii .gt. 0 .and. jj .gt. 0 ) then

               if( isb0  .ne. jmat(ii,jj,1) .or.
     &             nmed0 .ne. jmat(ii,jj,2) ) then

                  if( ibon(ii,jj,2) .gt. 0 ) then

                        vxx = vx(ii) + xbnd( ibon(ii,jj,2) )
                        vyy = vy(jj) + ybnd( ibon(ii,jj,2) )

                        ip = ip + 1
                        write(ioh) vxx, vyy

                  else

                        delx = xdel / 2.0
                        dely = ydel / 2.0
                        ipsm = 1
                        vxx  = vx(ii) + delx
                        vyy  = vy(jj) + dely

                     do k = 1, nbnd

                        call xvwi(vxx,vyy,x(1),x(2),x(3),u,v,w)
                        call smats(0,iout,m,nsb,nl,lt,nlb,ltb,
     &                             nr,mr,kr,ns,ms,ks,
     &                             x,u,v,w,ceye,nmed,isb,ishad)

                        ipsm = 1
                        if( isb .eq. isb0 .and.
     &                      nmed .eq. nmed0 ) ipsm = -1

                        delx = delx / 2.0
                        dely = dely / 2.0
                        vxx  = vxx + ipsm * delx
                        vyy  = vyy + ipsm * dely

                     end do

                        ip = ip + 1
                        write(ioh) vxx, vyy

                        ibon(ii,jj,2) =
     &                  int( ( vxx - vx(ii) ) / delx / 2.0 ) + 1

                  end if

                     if( ibtm .gt. 0 .and. iout .ne. 1 ) then

                        ig = ig + 1
                        nmat(ii,jj,2) = ishad0

                     end if

               end if

            end if

*-----------------------------------------------------------------------
*           ic = 7 : first down check
*-----------------------------------------------------------------------

            ic = 7

                  ii = i
                  jj = j - 1

            if( jj .gt. 0 ) then

               if( isb0  .eq. jmat(ii,jj,1) .and.
     &             nmed0 .eq. jmat(ii,jj,2) ) then

                        imat(i,j) = imat(ii,jj)

                        goto 100

               else

                  if( ibon(ii,jj,3) .gt. 0 ) then

                        vxx = vx(ii)
                        vyy = vy(jj) + ybnd( ibon(ii,jj,3) )

                        ip = ip + 1
                        write(ioh) vxx, vyy

                  else

                        dely = ydel / 2.0
                        ipsm = 1
                        vxx  = vx(ii)
                        vyy  = vy(jj) + dely

                     do k = 1, nbnd

                        call xvwi(vxx,vyy,x(1),x(2),x(3),u,v,w)
                        call smats(0,iout,m,nsb,nl,lt,nlb,ltb,
     &                             nr,mr,kr,ns,ms,ks,
     &                             x,u,v,w,ceye,nmed,isb,ishad)

                        ipsm = 1
                        if( isb .eq. isb0 .and.
     &                      nmed .eq. nmed0 ) ipsm = -1

                        dely = dely / 2.0
                        vyy  = vyy + ipsm * dely

                     end do

                        ip = ip + 1
                        write(ioh) vxx, vyy

                        ibon(ii,jj,3) =
     &                  int( ( vyy - vy(jj) ) / dely / 2.0 ) + 1

                  end if

                     if( ibtm .gt. 0 .and. iout .ne. 1 ) then

                        ig = ig + 1
                        nmat(ii,jj,2) = ishad0

                     end if

               end if

            else

                        vxx = vx(i)
                        vyy = vy(j)

                        ip = ip + 1
                        write(ioh) vxx, vyy

            end if

*-----------------------------------------------------------------------
*        write new region
*-----------------------------------------------------------------------

               inum = inum + 2
               imat(i,j) = inum

               ic  = 8

*-----------------------------------------------------------------------
*        eight directions movement and final write information
*-----------------------------------------------------------------------

  300    continue

*-----------------------------------------------------------------------
*           iout = 1 : line
*                = 2 : color
*                = 3 : line + color
*-----------------------------------------------------------------------

      if( i .eq. i0 .and. j .eq. j0 .and.
     &  ( ic .eq. 5 .or. ic .eq. 6 .or. ic .eq. 7 ) ) then

                  ip = ip + 1
                  write(ioh) x0, y0

*-----------------------------------------------------------------------
*        color for legend
*-----------------------------------------------------------------------

            if( iout .ne. 1 ) then

                  col1 = dmhsb(nmed0,1)
                  col2 = dmhsb(nmed0,2)
                  col3 = dmhsb(nmed0,3)

               if( itcl(nmed0) .eq. 0 ) then

                  write(iot,'(/a)') huni7//
     &                 huni3//huni4(1:if4)//huni5//huni6(1:if6)

                  x1 = -1.0
                  y1 = -1.0

                  write(iot,'(1p2e14.6)') x1, y1

                  itcl(nmed0) = 1

               end if

            end if

*-----------------------------------------------------------------------
*        check real surfaces, plane or not
*-----------------------------------------------------------------------

                  iouts = iout

            if( ibtm .gt. 0 .and. iout .ne. 1 ) then

                  ishmin = int( shadd ) + 10
                  ishmax = 0
                  ishneg = 0
                  ishavr = 0
                  ishnm  = 0

               do jj = 1, ny
               do ii = 1, nx

               if( nmat(ii,jj,2) .ne. 0 ) then

                     ishad = nmat(ii,jj,2)
                     ishnm = ishnm + 1

                     if( ishnm .eq. 1 ) then

                        ishne0 = ishad

                     else

                        if( ishad * ishne0 .lt. 0 ) ishneg = 1

                     end if

                     ishad = abs( ishad )

                     ishavr = ishavr + ishad

                     if( ishad .gt. ishmax ) ishmax = ishad
                     if( ishad .lt. ishmin ) ishmin = ishad

                     nmat(ii,jj,2) = ishad

               end if
               end do
               end do

                  ishavr = nint( dble(ishavr) / dble(ishnm) )
                  ishdel = ishmax - ishmin

               if( ishneg .gt. 0 ) then

                  ibtm   = -1
                  if( iout .eq. 3 ) iouts = 1

               else if( ishmax - ishmin .le. ishdd ) then

                  ibtm   = 0
                  ishad0 = ishavr

               end if

            end if

*-----------------------------------------------------------------------
*        surfaces and lines
*-----------------------------------------------------------------------

               if( iouts .eq. 2 ) then

                  huni11 = huni9

               else

                  huni11 = huni8

               end if

*-----------------------------------------------------------------------
*           boundary is in it or not
*-----------------------------------------------------------------------

               ibond = 0

            if( iouts .ne. 2 ) then

                  rewind ioh

               do l = 1, ip

                  read(ioh) x0, y0

                     if( abs( x0 - xvmin ) .le. epsb .or.
     &                   abs( x0 - xvmax ) .le. epsb .or.
     &                   abs( y0 - yvmin ) .le. epsb .or.
     &                   abs( y0 - yvmax ) .le. epsb ) then

                        ibond = 1
                        goto 230

                     end if

               end do

  230          continue

            end if

*-----------------------------------------------------------------------
*        with shadow, ibtm = -1
*-----------------------------------------------------------------------

         if( ibtm .lt. 0 ) then

            call psshow(m,iout,nl,lt,nlb,ltb,
     &                  nr,mr,kr,ns,ms,ks,nx,ny,
     &                  nmat,imat,jmat,ibon,ibond,
     &                  inum,xmin,xdel,ymin,ydel,xbnd,ybnd,
     &                  ddel,xminf,xmaxf,yminf,ymaxf,
     &                  ioq,iog,iop,iot,nsb,itt,brt,ceye)

*-----------------------------------------------------------------------
*           iouts = 1
*-----------------------------------------------------------------------

            if( iouts .eq. 1 ) then

                  write(iot,'(/a)') huni11//ct(itt)

               if( ibond .eq. 1 ) then

                  write(iot,'("frame:")')
                  write(iot,'(1p2e14.6)') xminf, yminf
                  write(iot,'(1p2e14.6)') xmaxf, yminf
                  write(iot,'(1p2e14.6)') xmaxf, ymaxf
                  write(iot,'(1p2e14.6)') xminf, ymaxf

               end if

               write(iot,'("clip:")')

                  ipl = 1

               call ppline(ipl,ioh,iot,ip,ddel)

            end if

*-----------------------------------------------------------------------
*        no shadow
*-----------------------------------------------------------------------

         else if( ibtm .ge. 0 ) then

*-----------------------------------------------------------------------
*           line only
*-----------------------------------------------------------------------

            if( iouts .eq. 1 ) then

                  write(iot,'(/a)') huni11//ct(itt)

*-----------------------------------------------------------------------
*           real surfaces
*-----------------------------------------------------------------------

            else if( ibtm .gt. 0 ) then

                  write(iot,'(/a)') huni11//ct(itt)

*-----------------------------------------------------------------------
*           box surfaces
*-----------------------------------------------------------------------

            else if( ibtm .eq. 0 ) then

                     shad0 = dble( abs( ishad0 ) ) / shadd

                     colx = col1
                     coly = col2
                     colz = col3

                     call givcol(shad0,colx,coly,colz,brt)

                     write(huni6,'(a3,3f7.3,a1)')
     &               'cb(', colx, coly, colz, ')'

                     if6 = 25

                  write(iot,'(/a)') huni11//ct(itt)//' '//
     &                              huni6(1:if6)

            end if


*-----------------------------------------------------------------------
*           frame:
*-----------------------------------------------------------------------

               if( ibond .eq. 1 ) then

                  write(iot,'("frame:")')
                  write(iot,'(1p2e14.6)') xminf, yminf
                  write(iot,'(1p2e14.6)') xmaxf, yminf
                  write(iot,'(1p2e14.6)') xmaxf, ymaxf
                  write(iot,'(1p2e14.6)') xminf, ymaxf

               end if

*-----------------------------------------------------------------------
*           clip:
*-----------------------------------------------------------------------

               write(iot,'("clip:")')

                  ipl = 1

               call ppline(ipl,ioh,iot,ip,ddel)

*-----------------------------------------------------------------------
*           path: and bmap:
*-----------------------------------------------------------------------

               if( ibtm .gt. 0 .and. iouts .ne. 1 ) then

                  call ppbshow(iot,iog,iop,nx,ny,nmat,
     &                         xdel,ydel,xmin,ymin,
     &                         col1,col2,col3,brt,ishdel)

               end if

         end if

*-----------------------------------------------------------------------

               goto 100

      end if

*-----------------------------------------------------------------------
*        four directions movement
*-----------------------------------------------------------------------

  500    continue

               if( ic .eq. 1 ) goto 501
               if( ic .eq. 2 ) goto 502
               if( ic .eq. 3 ) goto 503
               if( ic .eq. 4 ) goto 504
               if( ic .eq. 5 ) goto 505
               if( ic .eq. 6 ) goto 506
               if( ic .eq. 7 ) goto 507
               if( ic .eq. 8 ) goto 508

*-----------------------------------------------------------------------
*        ic = 1
*-----------------------------------------------------------------------

  501    ic = 1

                        ii = i + 1
                        jj = j

            if( ii .le. nx ) then

               if( isb0  .eq. jmat(ii,jj,1) .and.
     &             nmed0 .eq. jmat(ii,jj,2) .and.
     &           ( imat(ii,jj) .eq. 0 .or.
     &             imat(ii,jj) .eq. inum ) ) then

                  if( j + 1 .le. ny ) then

                     do jl = j + 1, ny

                        if( isb0  .eq. jmat(ii,jl,1) .and.
     &                      nmed0 .eq. jmat(ii,jl,2) .and.
     &                    ( imat(ii,jl) .eq. 0 .or.
     &                      imat(ii,jl) .eq. inum ) ) then

                           if( imat(ii,jl) .eq. 0 ) then

                              imat(ii,jl) = inum

                              if( ibtm .gt. 0 .and. iout .ne. 1 ) then

                                 ig = ig + 1
                                 nmat(ii,jl,2) = nmat(ii,jl,1)

                              end if

                           end if

                        else

                              if( ibtm .gt. 0 .and. iout .ne. 1 ) then

                                 ig = ig + 1
                                 nmat(ii,jl,2) = nmat(ii,jl-1,1)

                              end if

                           goto 601

                        end if

                     end do

  601                continue

                  end if

                        i = ii
                        j = jj

                        imat(i,j) = inum

                        if( ibtm .gt. 0 .and. iout .ne. 1 ) then

                           ig = ig + 1
                           nmat(i,j,2) = nmat(i,j,1)

                        end if

                     if( j-1 .lt. 1 ) then

                        ic = 6

                     elseif( ibon(i,j-1,4) .gt. 0 ) then

                        ic = 7

                     else

                        ic = 6

                     end if

                        goto 300

               else

                  if( ibon(i,j,1) .gt. 0 ) then

                        xc(1) = vx(i) + xbnd( ibon(i,j,1) )
                        xc(2) = vy(j)

                  else

                        delx = xdel / 2.0
                        ipsm = 1
                        vxx  = vx(i) + delx
                        vyy  = vy(j)

                     do k = 1, nbnd

                        call xvwi(vxx,vyy,x(1),x(2),x(3),u,v,w)
                        call smats(0,iout,m,nsb,nl,lt,nlb,ltb,
     &                             nr,mr,kr,ns,ms,ks,
     &                             x,u,v,w,ceye,nmed,isb,ishad)

                        ipsm = -1
                        if( isb .eq. isb0 .and.
     &                      nmed .eq. nmed0 ) ipsm = 1

                        delx = delx / 2.0
                        vxx  = vxx + ipsm * delx

                     end do

                        xc(1) = vxx
                        xc(2) = vyy

                        ibon(i,j,1) =
     &                  int( ( vxx - vx(i) ) / delx / 2.0 ) + 1

                  end if

                     if( ibtm .gt. 0 .and. iout .ne. 1 ) then

                        ig = ig + 1
                        nmat(ii,jj,2) = nmat(i,j,1)

                     end if

                        ic = 2
                        goto 310

               end if

            else

                        xc(1) = vx(i)
                        xc(2) = vy(j)

                        ic = 3
                        goto 310

            end if

*-----------------------------------------------------------------------
*        ic = 2
*-----------------------------------------------------------------------

  502    ic = 2

                        ii = i + 1
                        jj = j + 1

            if( ii .le. nx .and. jj .le. ny ) then

               if( isb0  .ne. jmat(ii,jj,1) .or.
     &             nmed0 .ne. jmat(ii,jj,2) .or.
     &           ( imat(ii,jj) .ne. 0 .and.
     &             imat(ii,jj) .ne. inum ) ) then

                  if( ibon(i,j,2) .gt. 0 ) then

                        xc(1) = vx(i) + xbnd( ibon(i,j,2) )
                        xc(2) = vy(j) + ybnd( ibon(i,j,2) )

                  else

                        delx = xdel / 2.0
                        dely = ydel / 2.0
                        ipsm = 1
                        vxx  = vx(i) + delx
                        vyy  = vy(j) + dely

                     do k = 1, nbnd

                        call xvwi(vxx,vyy,x(1),x(2),x(3),u,v,w)
                        call smats(0,iout,m,nsb,nl,lt,nlb,ltb,
     &                             nr,mr,kr,ns,ms,ks,
     &                             x,u,v,w,ceye,nmed,isb,ishad)

                        ipsm = -1
                        if( isb .eq. isb0 .and.
     &                      nmed .eq. nmed0 ) ipsm = 1

                        delx = delx / 2.0
                        dely = dely / 2.0
                        vxx  = vxx + ipsm * delx
                        vyy  = vyy + ipsm * dely

                     end do

                        xc(1) = vxx
                        xc(2) = vyy

                        ibon(i,j,2) =
     &                  int( ( vxx - vx(i) ) / delx / 2.0 ) + 1

                  end if

                     if( ibtm .gt. 0 .and. iout .ne. 1 ) then

                        ig = ig + 1
                        nmat(ii,jj,2) = nmat(i,j,1)

                     end if

                        ic = 3
                        goto 310

               end if

            end if

*-----------------------------------------------------------------------
*        ic = 3
*-----------------------------------------------------------------------

  503    ic = 3

                        ii = i
                        jj = j + 1

            if( jj .le. ny ) then

               if( isb0  .eq. jmat(ii,jj,1) .and.
     &             nmed0 .eq. jmat(ii,jj,2) .and.
     &           ( imat(ii,jj) .eq. 0 .or.
     &             imat(ii,jj) .eq. inum ) ) then

                  if( i - 1 .ge. 1 ) then

                     do il = i - 1, 1, -1

                        if( isb0  .eq. jmat(il,jj,1) .and.
     &                      nmed0 .eq. jmat(il,jj,2) .and.
     &                    ( imat(il,jj) .eq. 0 .or.
     &                      imat(il,jj) .eq. inum ) ) then

                           if( imat(il,jj) .eq. 0 ) then

                              imat(il,jj) = inum

                              if( ibtm .gt. 0 .and. iout .ne. 1 ) then

                                 ig = ig + 1
                                 nmat(il,jj,2) = nmat(il,jj,1)

                              end if

                           end if

                        else

                              if( ibtm .gt. 0 .and. iout .ne. 1 ) then

                                 ig = ig + 1
                                 nmat(il,jj,2) = nmat(il+1,jj,1)

                              end if

                           goto 603

                        end if

                     end do

  603                continue

                  end if

                        i = ii
                        j = jj

                        imat(i,j) = inum

                        if( ibtm .gt. 0 .and. iout .ne. 1 ) then

                           ig = ig + 1
                           nmat(i,j,2) = nmat(i,j,1)

                        end if

                     if( ibon(i,j-1,2) .gt. 0 ) then

                        ic = 1

                     else

                        ic = 8

                     end if

                        goto 300

               else

                  if( ibon(i,j,3) .gt. 0 ) then

                        xc(1) = vx(i)
                        xc(2) = vy(j) + ybnd( ibon(i,j,3) )

                  else

                        dely = ydel / 2.0
                        ipsm = 1
                        vxx  = vx(i)
                        vyy  = vy(j) + dely

                     do k = 1, nbnd

                        call xvwi(vxx,vyy,x(1),x(2),x(3),u,v,w)
                        call smats(0,iout,m,nsb,nl,lt,nlb,ltb,
     &                             nr,mr,kr,ns,ms,ks,
     &                             x,u,v,w,ceye,nmed,isb,ishad)

                        ipsm = -1
                        if( isb .eq. isb0 .and.
     &                      nmed .eq. nmed0 ) ipsm = 1

                        dely = dely / 2.0
                        vyy  = vyy + ipsm * dely

                     end do

                        xc(1) = vxx
                        xc(2) = vyy

                        ibon(i,j,3) =
     &                  int( ( vyy - vy(j) ) / dely / 2.0 ) + 1

                  end if

                     if( ibtm .gt. 0 .and. iout .ne. 1 ) then

                        ig = ig + 1
                        nmat(ii,jj,2) = nmat(i,j,1)

                     end if

                        ic = 4
                        goto 310

               end if

            else

                        xc(1) = vx(i)
                        xc(2) = vy(j)

                        ic = 5
                        goto 310

            end if

*-----------------------------------------------------------------------
*        ic = 4
*-----------------------------------------------------------------------

  504    ic = 4

                        ii = i - 1
                        jj = j + 1

            if( ii .gt. 0 .and. jj .le. ny ) then

               if( isb0  .ne. jmat(ii,jj,1) .or.
     &             nmed0 .ne. jmat(ii,jj,2) .or.
     &           ( imat(ii,jj) .ne. 0 .and.
     &             imat(ii,jj) .ne. inum ) ) then

                  if( ibon(i,j,4) .gt. 0 ) then

                        xc(1) = vx(i) - xbnd( ibon(i,j,4) )
                        xc(2) = vy(j) + ybnd( ibon(i,j,4) )

                  else

                        delx = -xdel / 2.0
                        dely =  ydel / 2.0
                        ipsm = 1
                        vxx  = vx(i) + delx
                        vyy  = vy(j) + dely

                     do k = 1, nbnd

                        call xvwi(vxx,vyy,x(1),x(2),x(3),u,v,w)
                        call smats(0,iout,m,nsb,nl,lt,nlb,ltb,
     &                             nr,mr,kr,ns,ms,ks,
     &                             x,u,v,w,ceye,nmed,isb,ishad)

                        ipsm = -1
                        if( isb .eq. isb0 .and.
     &                      nmed .eq. nmed0 ) ipsm = 1

                        delx = delx / 2.0
                        dely = dely / 2.0
                        vxx  = vxx + ipsm * delx
                        vyy  = vyy + ipsm * dely

                     end do

                        xc(1) = vxx
                        xc(2) = vyy

                        ibon(i,j,4) =
     &                  int( ( vyy - vy(j) ) / dely / 2.0 ) + 1

                  end if

                     if( ibtm .gt. 0 .and. iout .ne. 1 ) then

                        ig = ig + 1
                        nmat(ii,jj,2) = nmat(i,j,1)

                     end if

                        ic = 5
                        goto 310

               end if

            end if

*-----------------------------------------------------------------------
*        ic = 5
*-----------------------------------------------------------------------

  505    ic = 5

                        ii = i - 1
                        jj = j

            if( ii .gt. 0 ) then

               if( isb0  .eq. jmat(ii,jj,1) .and.
     &             nmed0 .eq. jmat(ii,jj,2) .and.
     &           ( imat(ii,jj) .eq. 0 .or.
     &             imat(ii,jj) .eq. inum ) ) then

                  if( j - 1 .ge. 1 ) then

                     do jl = j - 1, 1, -1

                        if( isb0  .eq. jmat(ii,jl,1) .and.
     &                      nmed0 .eq. jmat(ii,jl,2) .and.
     &                    ( imat(ii,jl) .eq. 0 .or.
     &                      imat(ii,jl) .eq. inum ) ) then

                           if( imat(ii,jl) .eq. 0 ) then

                              imat(ii,jl) = inum

                              if( ibtm .gt. 0 .and. iout .ne. 1 ) then

                                 ig = ig + 1
                                 nmat(ii,jl,2) = nmat(ii,jl,1)

                              end if

                           end if

                        else

                              if( ibtm .gt. 0 .and. iout .ne. 1 ) then

                                 ig = ig + 1
                                 nmat(ii,jl,2) = nmat(ii,jl+1,1)

                              end if

                           goto 605

                        end if

                     end do

  605                continue

                  end if

                        i = ii
                        j = jj

                        imat(i,j) = inum

                        if( ibtm .gt. 0 .and. iout .ne. 1 ) then

                           ig = ig + 1
                           nmat(i,j,2) = nmat(i,j,1)

                        end if
                     if( i-1 .lt. 1 ) then

                        ic = 2

                     elseif( ibon(i-1,j,2) .gt. 0 ) then

                        ic = 3

                     else

                        ic = 2

                     end if

                        goto 300

               else

                  if( ibon(ii,jj,1) .gt. 0 ) then

                        xc(1) = vx(ii) + xbnd( ibon(ii,jj,1) )
                        xc(2) = vy(jj)

                  else

                        delx = xdel / 2.0
                        ipsm = 1
                        vxx  = vx(ii) + delx
                        vyy  = vy(jj)

                     do k = 1, nbnd

                        call xvwi(vxx,vyy,x(1),x(2),x(3),u,v,w)
                        call smats(0,iout,m,nsb,nl,lt,nlb,ltb,
     &                             nr,mr,kr,ns,ms,ks,
     &                             x,u,v,w,ceye,nmed,isb,ishad)

                        ipsm = 1
                        if( isb .eq. isb0 .and.
     &                      nmed .eq. nmed0 ) ipsm = -1

                        delx = delx / 2.0
                        vxx  = vxx + ipsm * delx

                     end do

                        xc(1) = vxx
                        xc(2) = vyy

                        ibon(ii,jj,1) =
     &                  int( ( vxx - vx(ii) ) / delx / 2.0 ) + 1

                  end if

                     if( ibtm .gt. 0 .and. iout .ne. 1 ) then

                        ig = ig + 1
                        nmat(ii,jj,2) = nmat(i,j,1)

                     end if

                        ic = 6
                        goto 310

               end if

            else

                        xc(1) = vx(i)
                        xc(2) = vy(j)

                        ic = 7
                        goto 310

            end if

*-----------------------------------------------------------------------
*        ic = 6
*-----------------------------------------------------------------------

  506    ic = 6

                        ii = i - 1
                        jj = j - 1

            if( ii .gt. 0 .and. jj .gt. 0 ) then

               if( isb0  .ne. jmat(ii,jj,1) .or.
     &             nmed0 .ne. jmat(ii,jj,2) .or.
     &           ( imat(ii,jj) .ne. 0 .and.
     &             imat(ii,jj) .ne. inum ) ) then

                  if( ibon(ii,jj,2) .gt. 0 ) then

                        xc(1) = vx(ii) + xbnd( ibon(ii,jj,2) )
                        xc(2) = vy(jj) + ybnd( ibon(ii,jj,2) )

                  else

                        delx = xdel / 2.0
                        dely = ydel / 2.0
                        ipsm = 1
                        vxx  = vx(ii) + delx
                        vyy  = vy(jj) + dely

                     do k = 1, nbnd

                        call xvwi(vxx,vyy,x(1),x(2),x(3),u,v,w)
                        call smats(0,iout,m,nsb,nl,lt,nlb,ltb,
     &                             nr,mr,kr,ns,ms,ks,
     &                             x,u,v,w,ceye,nmed,isb,ishad)

                        ipsm = 1
                        if( isb .eq. isb0 .and.
     &                      nmed .eq. nmed0 ) ipsm = -1

                        delx = delx / 2.0
                        dely = dely / 2.0
                        vxx  = vxx + ipsm * delx
                        vyy  = vyy + ipsm * dely

                     end do

                        xc(1) = vxx
                        xc(2) = vyy

                        ibon(ii,jj,2) =
     &                  int( ( vxx - vx(ii) ) / delx / 2.0 ) + 1

                  end if

                     if( ibtm .gt. 0 .and. iout .ne. 1 ) then

                        ig = ig + 1
                        nmat(ii,jj,2) = nmat(i,j,1)

                     end if

                        ic = 7
                        goto 310

               end if

            end if

*-----------------------------------------------------------------------
*        ic = 7
*-----------------------------------------------------------------------

  507    ic = 7

                        ii = i
                        jj = j - 1

            if( jj .gt. 0 ) then

               if( isb0  .eq. jmat(ii,jj,1) .and.
     &             nmed0 .eq. jmat(ii,jj,2) .and.
     &           ( imat(ii,jj) .eq. 0 .or.
     &             imat(ii,jj) .eq. inum ) ) then

                  if( i + 1 .le. nx ) then

                     do il = i + 1, nx

                        if( isb0  .eq. jmat(il,jj,1) .and.
     &                      nmed0 .eq. jmat(il,jj,2) .and.
     &                    ( imat(il,jj) .eq. 0 .or.
     &                      imat(il,jj) .eq. inum ) ) then

                           if( imat(il,jj) .eq. 0 ) then

                              imat(il,jj) = inum

                              if( ibtm .gt. 0 .and. iout .ne. 1 ) then

                                 ig = ig + 1
                                 nmat(il,jj,2) = nmat(il,jj,1)

                              end if

                           end if

                        else

                              if( ibtm .gt. 0 .and. iout .ne. 1 ) then

                                 ig = ig + 1
                                 nmat(il,jj,2) = nmat(il-1,jj,1)

                              end if

                           goto 607

                        end if

                     end do

  607                continue

                  end if

                        i = ii
                        j = jj

                        imat(i,j) = inum

                        if( ibtm .gt. 0 .and. iout .ne. 1 ) then

                           ig = ig + 1
                           nmat(i,j,2) = nmat(i,j,1)

                        end if

                     if( ibon(i,j,4) .gt. 0 ) then

                        ic = 5

                     else

                        ic = 4

                     end if

                        goto 300

               else

                  if( ibon(ii,jj,3) .gt. 0 ) then

                        xc(1) = vx(ii)
                        xc(2) = vy(jj) + ybnd( ibon(ii,jj,3) )

                  else

                        dely = ydel / 2.0
                        ipsm = 1
                        vxx  = vx(ii)
                        vyy  = vy(jj) + dely

                     do k = 1, nbnd

                        call xvwi(vxx,vyy,x(1),x(2),x(3),u,v,w)
                        call smats(0,iout,m,nsb,nl,lt,nlb,ltb,
     &                             nr,mr,kr,ns,ms,ks,
     &                             x,u,v,w,ceye,nmed,isb,ishad)

                        ipsm = 1
                        if( isb .eq. isb0 .and.
     &                      nmed .eq. nmed0 ) ipsm = -1

                        dely = dely / 2.0
                        vyy  = vyy + ipsm * dely

                     end do

                        xc(1) = vxx
                        xc(2) = vyy

                        ibon(ii,jj,3) =
     &                  int( ( vyy - vy(jj) ) / dely / 2.0 ) + 1

                  end if

                     if( ibtm .gt. 0 .and. iout .ne. 1 ) then

                        ig = ig + 1
                        nmat(ii,jj,2) = nmat(i,j,1)

                     end if

                        ic = 8
                        goto 310

               end if

            else

                        xc(1) = vx(i)
                        xc(2) = vy(j)

                        ic = 1
                        goto 310

            end if

*-----------------------------------------------------------------------
*        ic = 8
*-----------------------------------------------------------------------

  508    ic = 8

                     ii = i + 1
                     jj = j - 1

            if( ii .le. nx .and. jj .gt. 0 ) then

               if( isb0  .ne. jmat(ii,jj,1) .or.
     &             nmed0 .ne. jmat(ii,jj,2) .or.
     &           ( imat(ii,jj) .ne. 0 .and.
     &             imat(ii,jj) .ne. inum ) ) then

                  if( ibon(ii,jj,4) .gt. 0 ) then

                        xc(1) = vx(ii) - xbnd( ibon(ii,jj,4) )
                        xc(2) = vy(jj) + ybnd( ibon(ii,jj,4) )

                  else

                        delx = -xdel / 2.0
                        dely =  ydel / 2.0
                        ipsm = 1
                        vxx  = vx(ii) + delx
                        vyy  = vy(jj) + dely

                     do k = 1, nbnd

                        call xvwi(vxx,vyy,x(1),x(2),x(3),u,v,w)
                        call smats(0,iout,m,nsb,nl,lt,nlb,ltb,
     &                             nr,mr,kr,ns,ms,ks,
     &                             x,u,v,w,ceye,nmed,isb,ishad)

                        ipsm = 1
                        if( isb .eq. isb0 .and.
     &                      nmed .eq. nmed0 ) ipsm = -1

                        delx = delx / 2.0
                        dely = dely / 2.0
                        vxx  = vxx + ipsm * delx
                        vyy  = vyy + ipsm * dely

                     end do

                        xc(1) = vxx
                        xc(2) = vyy

                        ibon(ii,jj,4) =
     &                  int( ( vyy - vy(jj) ) / dely / 2.0 ) + 1

                  end if

                     if( ibtm .gt. 0 .and. iout .ne. 1 ) then

                        ig = ig + 1
                        nmat(ii,jj,2) = nmat(i,j,1)

                     end if

                        ic = 1
                        goto 310

               end if

            end if

                  goto 501

*-----------------------------------------------------------------------
*        write boundary points on ioh
*-----------------------------------------------------------------------

  310       continue

                  ip = ip + 1
                  write(ioh) xc(1), xc(2)

            goto 300

*-----------------------------------------------------------------------

  900 continue

*-----------------------------------------------------------------------

      if( iout .gt. 0 ) close( ioh )
      if( iout .gt. 1 ) close( iog )
      if( iout .gt. 1 ) close( iop )

*-----------------------------------------------------------------------

            close(iot)

         if( iteps(m) .ne. 0 ) then

            open(iot, file = fname, status = 'unknown' )
            call a_angel(idasa,fname)

         end if

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine psshow(m,iout,nl,lt,nlb,ltb,
     &                  nr,mr,kr,ns,ms,ks,nx,ny,
     &                  nmat,imat,jmat,ibon,ibond,
     &                  inum,xmin,xdel,ymin,ydel,xbnd,ybnd,
     &                  ddel,xminf,xmaxf,yminf,ymaxf,
     &                  ioh,iog,iop,iot,nsb,itt,brt,ceye)
*                                                                      *
*       sub program of output 3d geometry                              *
*       last modified by K.Niita on 2002/11/14                         *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'

      parameter ( pi = 3.141592653589793d0 )

*-----------------------------------------------------------------------

      common /shadw/  shadd, ishdd
      common /cubdb/  nbnd

*-----------------------------------------------------------------------

      dimension   lt(nl)
      dimension   ltb(nlb)
      dimension   kr(mr)
      dimension   ks(ms)
      dimension   nmat(nx,ny,2)
      dimension   imat(nx,ny)
      dimension   jmat(nx,ny,2)
      dimension   ibon(nx,ny,4)

*-----------------------------------------------------------------------

      dimension x(3), xc(3), xd(3), angb(3), angd(3), ceye(3)
      dimension xbnd(16), ybnd(16)

*-----------------------------------------------------------------------

      character huni11*16
      character huni9*16
      data huni9 /'hb: noline clip '/

      character huni4*80
      character huni6*30

      character ct(7)*3
      data ct /'zzz','zz ','z  ','   ','t  ',
     &         'tt ','ttt'/

      common /mtnmcl/ dmhsb(-1:kvlmax,4),
     &                dmtnm(-1:kvlmax), dmtcl(-1:kvlmax),
     &                nmtnm(-1:kvlmax), nmtcl(-1:kvlmax)
      character dmtnm*80, dmtcl*30

      dimension  brt(2)

*-----------------------------------------------------------------------
*     internal function
*-----------------------------------------------------------------------

               vx(ik) = xmin + xdel * dble( ik - 1 )
               vy(ik) = ymin + ydel * dble( ik - 1 )

*-----------------------------------------------------------------------

         goto 110

*-----------------------------------------------------------------------
*     start check each cell
*-----------------------------------------------------------------------

  100 continue

               do j = 1, ny
               do i = 1, nx

                  if( imat(i,j) .eq. -inum ) imat(i,j) = -inum - 1

               end do
               end do

  110 continue

               do j = 1, ny
               do i = 1, nx

                  if( imat(i,j) .eq. inum ) goto 200

               end do
               end do

               goto 900

*-----------------------------------------------------------------------
*        start new region
*-----------------------------------------------------------------------

  200    continue

               isb0   = jmat(i,j,1)
               nmed0  = jmat(i,j,2)
               ishad0 = nmat(i,j,1)

               ibtm  = 1
               if( isb0 .le. 0 ) ibtm = 0

               rewind ioh

               i0 = i
               j0 = j

               ig = 0
               ip = 0

               do jj = 1, ny
               do ii = 1, nx

                  nmat(ii,jj,2) = 0

               end do
               end do

               if( ibtm .gt. 0 ) then

                  ig = ig + 1
                  nmat(i,j,2) = ishad0

               end if

               if4   = nmtnm(nmed0)
               huni4 = dmtnm(nmed0)
               if6   = nmtcl(nmed0)+3
               huni6(1:2) = 'c['
               huni6(2+1:2+nmtcl(nmed0))
     &             = dmtcl(nmed0)(1:nmtcl(nmed0))
               huni6(3+nmtcl(nmed0):3+nmtcl(nmed0)) = ']'

*-----------------------------------------------------------------------
*           ic = 5 : first left check
*-----------------------------------------------------------------------

            ic = 5

                  ii = i - 1
                  jj = j

            if( ii .gt. 0 ) then

                  if( ibon(ii,jj,1) .gt. 0 ) then

                        x0 = vx(ii) + xbnd( ibon(ii,jj,1) )
                        y0 = vy(jj)

                  else

                        delx = xdel / 2.0
                        ipsm = 1
                        vxx  = vx(ii) + delx
                        vyy  = vy(jj)

                     do k = 1, nbnd

                        call xvwi(vxx,vyy,x(1),x(2),x(3),u,v,w)
                        call smats(1,iout,m,nsb,nl,lt,nlb,ltb,
     &                             nr,mr,kr,ns,ms,ks,
     &                             x,u,v,w,ceye,nmed,isb,ishad)

                        ipsm = 1
                        if( isb .eq. isb0 .and.
     &                      ishad * ishad0 .gt. 0 .and.
     &                      nmed .eq. nmed0 ) ipsm = -1

                        delx = delx / 2.0
                        vxx  = vxx + ipsm * delx

                     end do

                        x0 = vxx
                        y0 = vyy

                        ibon(ii,jj,1) =
     &                  int( ( vxx - vx(ii) ) / delx / 2.0 ) + 1

                  end if

                     if( ibtm .gt. 0 ) then

                        ig = ig + 1
                        nmat(ii,jj,2) = ishad0

                     end if

            else

                        x0 = vx(i)
                        y0 = vy(j)

            end if

*-----------------------------------------------------------------------
*           ic = 6 : first left-down check
*-----------------------------------------------------------------------

            ic = 6

                  icm = 0

                  ii = i - 1
                  jj = j - 1

            if( ii .gt. 0 .and. jj .gt. 0 ) then

                     icm = 1

                  if( ibon(ii,jj,2) .gt. 0 ) then

                        vxx = vx(ii) + xbnd( ibon(ii,jj,2) )
                        vyy = vy(jj) + ybnd( ibon(ii,jj,2) )

                        ip = ip + 1
                        write(ioh) vxx, vyy

                  else

                        delx = xdel / 2.0
                        dely = ydel / 2.0
                        ipsm = 1
                        vxx  = vx(ii) + delx
                        vyy  = vy(jj) + dely

                     do k = 1, nbnd

                        call xvwi(vxx,vyy,x(1),x(2),x(3),u,v,w)
                        call smats(1,iout,m,nsb,nl,lt,nlb,ltb,
     &                             nr,mr,kr,ns,ms,ks,
     &                             x,u,v,w,ceye,nmed,isb,ishad)

                        ipsm = 1
                        if( isb .eq. isb0 .and.
     &                      ishad * ishad0 .gt. 0 .and.
     &                      nmed .eq. nmed0 ) ipsm = -1

                        delx = delx / 2.0
                        dely = dely / 2.0
                        vxx  = vxx + ipsm * delx
                        vyy  = vyy + ipsm * dely

                     end do

                        ip = ip + 1
                        write(ioh) vxx, vyy

                        ibon(ii,jj,2) =
     &                  int( ( vxx - vx(ii) ) / delx / 2.0 ) + 1

                  end if

                     if( ibtm .gt. 0 ) then

                        ig = ig + 1
                        nmat(ii,jj,2) = ishad0

                     end if

            end if

*-----------------------------------------------------------------------
*           ic = 7 : first down check
*-----------------------------------------------------------------------

            ic = 7

                  ii = i
                  jj = j - 1

            if( jj .gt. 0 ) then

                  if( ibon(ii,jj,3) .gt. 0 ) then

                        vxx = vx(ii)
                        vyy = vy(jj) + ybnd( ibon(ii,jj,3) )

                        ip = ip + 1
                        write(ioh) vxx, vyy

                  else

                        dely = ydel / 2.0
                        ipsm = 1
                        vxx  = vx(ii)
                        vyy  = vy(jj) + dely

                     do k = 1, nbnd

                        call xvwi(vxx,vyy,x(1),x(2),x(3),u,v,w)
                        call smats(1,iout,m,nsb,nl,lt,nlb,ltb,
     &                             nr,mr,kr,ns,ms,ks,
     &                             x,u,v,w,ceye,nmed,isb,ishad)

                        ipsm = 1
                        if( isb .eq. isb0 .and.
     &                      ishad * ishad0 .gt. 0 .and.
     &                      nmed .eq. nmed0 ) ipsm = -1

                        dely = dely / 2.0
                        vyy  = vyy + ipsm * dely

                     end do

                        ip = ip + 1
                        write(ioh) vxx, vyy

                        ibon(ii,jj,3) =
     &                  int( ( vyy - vy(jj) ) / dely / 2.0 ) + 1

                  end if

                     if( ibtm .gt. 0 ) then

                        ig = ig + 1
                        nmat(ii,jj,2) = ishad0

                     end if

            else

                        vxx = vx(i)
                        vyy = vy(j)

                        ip = ip + 1
                        write(ioh) vxx, vyy

            end if

*-----------------------------------------------------------------------
*        write new region
*-----------------------------------------------------------------------

               imat(i,j) = -inum

               ic = 8

*-----------------------------------------------------------------------
*        eight directions movement and final write information
*-----------------------------------------------------------------------

  300    continue

*-----------------------------------------------------------------------
*           iouts = 2 : color
*-----------------------------------------------------------------------

      if( i .eq. i0 .and. j .eq. j0 .and.
     &  ( ic .eq. 5 .or. ic .eq. 6 .or. ic .eq. 7 ) ) then

                  ip = ip + 1
                  write(ioh) x0, y0

*-----------------------------------------------------------------------
*        check real surfaces, plane or not
*-----------------------------------------------------------------------

            if( ibtm .gt. 0 ) then

                  ishmin = int( shadd ) + 10
                  ishmax = 0
                  ishavr = 0
                  ishnm  = 0

               do jj = 1, ny
               do ii = 1, nx

               if( nmat(ii,jj,2) .ne. 0 ) then

                  ishad = abs( nmat(ii,jj,2) )
                  ishnm = ishnm + 1

                  ishavr = ishavr + ishad

                  if( ishad .gt. ishmax ) ishmax = ishad
                  if( ishad .lt. ishmin ) ishmin = ishad

                  nmat(ii,jj,2) = ishad

               end if
               end do
               end do

                  ishavr = nint( dble(ishavr) / dble(ishnm) )
                  ishdel = ishmax - ishmin

               if( ishmax - ishmin .le. ishdd ) then

                  ibtm   = 0
                  ishad0 = ishavr

               end if

            end if

*-----------------------------------------------------------------------
*        surfaces
*-----------------------------------------------------------------------

                  huni11 = huni9

                  col1 = dmhsb(nmed0,1)
                  col2 = dmhsb(nmed0,2)
                  col3 = dmhsb(nmed0,3)

*-----------------------------------------------------------------------
*           real surfaces
*-----------------------------------------------------------------------

            if( ibtm .gt. 0 ) then

                  write(iot,'(/a)') huni11//ct(itt)

*-----------------------------------------------------------------------
*           box surfaces
*-----------------------------------------------------------------------

            else if( ibtm .eq. 0 ) then

                     shad0 = dble( abs( ishad0 ) ) / shadd

                     colx = col1
                     coly = col2
                     colz = col3

                     call givcol(shad0,colx,coly,colz,brt)

                     write(huni6,'(a3,3f7.3,a1)')
     &               'cb(', colx, coly, colz, ')'

                     if6 = 25

                  write(iot,'(/a)') huni11//ct(itt)//' '//
     &                              huni6(1:if6)

            end if


*-----------------------------------------------------------------------
*           frame:
*-----------------------------------------------------------------------

               if( ibond .eq. 1 ) then

                  write(iot,'("frame:")')
                  write(iot,'(1p2e14.6)') xminf, yminf
                  write(iot,'(1p2e14.6)') xmaxf, yminf
                  write(iot,'(1p2e14.6)') xmaxf, ymaxf
                  write(iot,'(1p2e14.6)') xminf, ymaxf

               end if

*-----------------------------------------------------------------------

               write(iot,'("clip:")')

                  ipl = 1

               call ppline(ipl,ioh,iot,ip,ddel)

*-----------------------------------------------------------------------

            if( ibtm .gt. 0 ) then

               call ppbshow(iot,iog,iop,nx,ny,nmat,xdel,ydel,xmin,ymin,
     &                      col1,col2,col3,brt,ishdel)

            end if

*-----------------------------------------------------------------------

               goto 100

      end if

*-----------------------------------------------------------------------
*        four directions movement
*-----------------------------------------------------------------------

  500    continue

               if( ic .eq. 1 ) goto 501
               if( ic .eq. 2 ) goto 502
               if( ic .eq. 3 ) goto 503
               if( ic .eq. 4 ) goto 504
               if( ic .eq. 5 ) goto 505
               if( ic .eq. 6 ) goto 506
               if( ic .eq. 7 ) goto 507
               if( ic .eq. 8 ) goto 508

*-----------------------------------------------------------------------
*        ic = 1
*-----------------------------------------------------------------------

  501    ic = 1

                        ii = i + 1
                        jj = j

            if( ii .le. nx ) then

               if( inum .eq. abs(imat(ii,jj)) .and.
     &             ishad0 * nmat(ii,jj,1) .gt. 0 ) then

                  if( j + 1 .le. ny ) then

                     do jl = j + 1, ny

                        if( inum .eq. abs(imat(ii,jl)) .and.
     &                      ishad0 * nmat(ii,jl,1) .gt. 0 ) then

                           if( imat(ii,jl) .eq. inum ) then

                              imat(ii,jl) = -inum

                              if( ibtm .gt. 0 ) then

                                 ig = ig + 1
                                 nmat(ii,jl,2) = nmat(ii,jl,1)

                              end if

                           end if

                        else

                              if( ibtm .gt. 0 ) then

                                 ig = ig + 1
                                 nmat(ii,jl,2) = nmat(ii,jl-1,1)

                              end if

                           goto 601

                        end if

                     end do

  601                continue

                  end if

                        i = ii
                        j = jj

                        imat(i,j) = -inum

                        if( ibtm .gt. 0 ) then

                           ig = ig + 1
                           nmat(i,j,2) =  nmat(i,j,1)

                        end if

                     if( j-1 .lt. 1 ) then

                        ic = 6

                     elseif( ibon(i,j-1,4) .gt. 0 ) then

                        ic = 7

                     else

                        ic = 6

                     end if

                        goto 300

               else

                  if( ibon(i,j,1) .gt. 0 ) then

                        xc(1) = vx(i) + xbnd( ibon(i,j,1) )
                        xc(2) = vy(j)

                  else

                        delx = xdel / 2.0
                        ipsm = 1
                        vxx  = vx(i) + delx
                        vyy  = vy(j)

                     do k = 1, nbnd

                        call xvwi(vxx,vyy,x(1),x(2),x(3),u,v,w)
                        call smats(1,iout,m,nsb,nl,lt,nlb,ltb,
     &                             nr,mr,kr,ns,ms,ks,
     &                             x,u,v,w,ceye,nmed,isb,ishad)

                        ipsm = -1
                        if( isb .eq. isb0 .and.
     &                      ishad * ishad0 .gt. 0 .and.
     &                      nmed .eq. nmed0 ) ipsm = 1

                        delx = delx / 2.0
                        vxx  = vxx + ipsm * delx

                     end do

                        xc(1) = vxx
                        xc(2) = vyy

                        ibon(i,j,1) =
     &                  int( ( vxx - vx(i) ) / delx / 2.0 ) + 1

                  end if

                     if( ibtm .gt. 0 ) then

                        ig = ig + 1
                        nmat(ii,jj,2) =  nmat(i,j,1)

                     end if

                        ic = 2
                        goto 310

               end if

            else

                        xc(1) = vx(i)
                        xc(2) = vy(j)

                        ic = 3
                        goto 310

            end if

*-----------------------------------------------------------------------
*        ic = 2
*-----------------------------------------------------------------------

  502    ic = 2

                        ii = i + 1
                        jj = j + 1

            if( ii .le. nx .and. jj .le. ny ) then

               if( inum .ne. abs(imat(ii,jj)) .or.
     &             ishad0 * nmat(ii,jj,1) .lt. 0 ) then

                  if( ibon(i,j,2) .gt. 0 ) then

                        xc(1) = vx(i) + xbnd( ibon(i,j,2) )
                        xc(2) = vy(j) + ybnd( ibon(i,j,2) )

                  else

                        delx = xdel / 2.0
                        dely = ydel / 2.0
                        ipsm = 1
                        vxx  = vx(i) + delx
                        vyy  = vy(j) + dely

                     do k = 1, nbnd

                        call xvwi(vxx,vyy,x(1),x(2),x(3),u,v,w)
                        call smats(1,iout,m,nsb,nl,lt,nlb,ltb,
     &                             nr,mr,kr,ns,ms,ks,
     &                             x,u,v,w,ceye,nmed,isb,ishad)

                        ipsm = -1
                        if( isb .eq. isb0 .and.
     &                      ishad * ishad0 .gt. 0 .and.
     &                      nmed .eq. nmed0 ) ipsm = 1

                        delx = delx / 2.0
                        dely = dely / 2.0
                        vxx  = vxx + ipsm * delx
                        vyy  = vyy + ipsm * dely

                     end do

                        xc(1) = vxx
                        xc(2) = vyy

                        ibon(i,j,2) =
     &                  int( ( vxx - vx(i) ) / delx / 2.0 ) + 1

                  end if

                     if( ibtm .gt. 0 ) then

                        ig = ig + 1
                        nmat(ii,jj,2) = nmat(i,j,1)

                     end if

                        ic = 3
                        goto 310

               end if

            end if

*-----------------------------------------------------------------------
*        ic = 3
*-----------------------------------------------------------------------

  503    ic = 3

                        ii = i
                        jj = j + 1

            if( jj .le. ny ) then

               if( inum .eq. abs(imat(ii,jj)) .and.
     &             ishad0 * nmat(ii,jj,1) .gt. 0 ) then

                  if( i - 1 .ge. 1 ) then

                     do il = i - 1, 1, -1

                        if( inum .eq. abs(imat(il,jj)) .and.
     &                      ishad0 * nmat(il,jj,1) .gt. 0 ) then

                           if( imat(il,jj) .eq. inum ) then

                              imat(il,jj) = -inum

                              if( ibtm .gt. 0 ) then

                                 ig = ig + 1
                                 nmat(il,jj,2) = nmat(il,jj,1)

                              end if

                           end if

                        else

                              if( ibtm .gt. 0 ) then

                                 ig = ig + 1
                                 nmat(il,jj,2) = nmat(il+1,jj,1)

                              end if

                           goto 603

                        end if

                     end do

  603                continue

                  end if

                        i = ii
                        j = jj

                        imat(i,j) = -inum

                        if( ibtm .gt. 0 ) then

                           ig = ig + 1
                           nmat(i,j,2) = nmat(i,j,1)

                        end if

                     if( ibon(i,j-1,2) .gt. 0 ) then

                        ic = 1

                     else

                        ic = 8

                     end if

                        goto 300

               else

                  if( ibon(i,j,3) .gt. 0 ) then

                        xc(1) = vx(i)
                        xc(2) = vy(j) + ybnd( ibon(i,j,3) )

                  else

                        dely = ydel / 2.0
                        ipsm = 1
                        vxx  = vx(i)
                        vyy  = vy(j) + dely

                     do k = 1, nbnd

                        call xvwi(vxx,vyy,x(1),x(2),x(3),u,v,w)
                        call smats(1,iout,m,nsb,nl,lt,nlb,ltb,
     &                             nr,mr,kr,ns,ms,ks,
     &                             x,u,v,w,ceye,nmed,isb,ishad)

                        ipsm = -1
                        if( isb .eq. isb0 .and.
     &                      ishad * ishad0 .gt. 0 .and.
     &                      nmed .eq. nmed0 ) ipsm = 1

                        dely = dely / 2.0
                        vyy  = vyy + ipsm * dely

                     end do

                        xc(1) = vxx
                        xc(2) = vyy

                        ibon(i,j,3) =
     &                  int( ( vyy - vy(j) ) / dely / 2.0 ) + 1

                  end if

                     if( ibtm .gt. 0 ) then

                        ig = ig + 1
                        nmat(ii,jj,2) = nmat(i,j,1)

                     end if

                        ic = 4
                        goto 310

               end if

            else

                        xc(1) = vx(i)
                        xc(2) = vy(j)

                        ic = 5
                        goto 310

            end if

*-----------------------------------------------------------------------
*        ic = 4
*-----------------------------------------------------------------------

  504    ic = 4

                        ii = i - 1
                        jj = j + 1

            if( ii .gt. 0 .and. jj .le. ny ) then

               if( inum .ne. abs(imat(ii,jj)) .or.
     &             ishad0 * nmat(ii,jj,1) .lt. 0 ) then

                  if( ibon(i,j,4) .gt. 0 ) then

                        xc(1) = vx(i) - xbnd( ibon(i,j,4) )
                        xc(2) = vy(j) + ybnd( ibon(i,j,4) )

                  else

                        delx = -xdel / 2.0
                        dely =  ydel / 2.0
                        ipsm = 1
                        vxx  = vx(i) + delx
                        vyy  = vy(j) + dely

                     do k = 1, nbnd

                        call xvwi(vxx,vyy,x(1),x(2),x(3),u,v,w)
                        call smats(1,iout,m,nsb,nl,lt,nlb,ltb,
     &                             nr,mr,kr,ns,ms,ks,
     &                             x,u,v,w,ceye,nmed,isb,ishad)

                        ipsm = -1
                        if( isb .eq. isb0 .and.
     &                      ishad * ishad0 .gt. 0 .and.
     &                      nmed .eq. nmed0 ) ipsm = 1

                        delx = delx / 2.0
                        dely = dely / 2.0
                        vxx  = vxx + ipsm * delx
                        vyy  = vyy + ipsm * dely

                     end do

                        xc(1) = vxx
                        xc(2) = vyy

                        ibon(i,j,4) =
     &                  int( ( vyy - vy(j) ) / dely / 2.0 ) + 1

                  end if

                     if( ibtm .gt. 0 ) then

                        ig = ig + 1
                        nmat(ii,jj,2) = nmat(i,j,1)

                     end if

                        ic = 5
                        goto 310

               end if

            end if

*-----------------------------------------------------------------------
*        ic = 5
*-----------------------------------------------------------------------

  505    ic = 5

                        ii = i - 1
                        jj = j

            if( ii .gt. 0 ) then

               if( inum .eq. abs(imat(ii,jj)) .and.
     &             ishad0 * nmat(ii,jj,1) .gt. 0 ) then

                  if( j - 1 .ge. 1 ) then

                     do jl = j - 1, 1, -1

                        if( inum .eq. abs(imat(ii,jl)) .and.
     &                      ishad0 * nmat(ii,jl,1) .gt. 0 ) then

                           if( imat(ii,jl) .eq. inum ) then

                              imat(ii,jl) = -inum

                              if( ibtm .gt. 0 ) then

                                 ig = ig + 1
                                 nmat(ii,jl,2) = nmat(ii,jl,1)

                              end if

                           end if

                        else

                              if( ibtm .gt. 0 ) then

                                 ig = ig + 1
                                 nmat(ii,jl,2) = nmat(ii,jl+1,1)

                              end if

                           goto 605

                        end if

                     end do

  605                continue

                  end if

                        i = ii
                        j = jj

                        imat(i,j) = -inum

                        if( ibtm .gt. 0 ) then

                           ig = ig + 1
                           nmat(i,j,2) = nmat(i,j,1)

                        end if

                     if( i-1 .lt. 1 ) then

                        ic = 2

                     elseif( ibon(i-1,j,2) .gt. 0 ) then

                        ic = 3

                     else

                        ic = 2

                     end if

                        goto 300

               else

                  if( ibon(ii,jj,1) .gt. 0 ) then

                        xc(1) = vx(ii) + xbnd( ibon(ii,jj,1) )
                        xc(2) = vy(jj)

                  else

                        delx = xdel / 2.0
                        ipsm = 1
                        vxx  = vx(ii) + delx
                        vyy  = vy(jj)

                     do k = 1, nbnd

                        call xvwi(vxx,vyy,x(1),x(2),x(3),u,v,w)
                        call smats(1,iout,m,nsb,nl,lt,nlb,ltb,
     &                             nr,mr,kr,ns,ms,ks,
     &                             x,u,v,w,ceye,nmed,isb,ishad)

                        ipsm = 1
                        if( isb .eq. isb0 .and.
     &                      ishad * ishad0 .gt. 0 .and.
     &                      nmed .eq. nmed0 ) ipsm = -1

                        delx = delx / 2.0
                        vxx  = vxx + ipsm * delx

                     end do

                        xc(1) = vxx
                        xc(2) = vyy

                        ibon(ii,jj,1) =
     &                  int( ( vxx - vx(ii) ) / delx / 2.0 ) + 1

                  end if

                     if( ibtm .gt. 0 ) then

                        ig = ig + 1
                        nmat(ii,jj,2) = nmat(i,j,1)

                     end if

                        ic = 6
                        goto 310

               end if

            else

                        xc(1) = vx(i)
                        xc(2) = vy(j)

                        ic = 7
                        goto 310

            end if

*-----------------------------------------------------------------------
*        ic = 6
*-----------------------------------------------------------------------

  506    ic = 6

                        ii = i - 1
                        jj = j - 1

            if( ii .gt. 0 .and. jj .gt. 0 ) then

               if( inum .ne. abs(imat(ii,jj)) .or.
     &             ishad0 * nmat(ii,jj,1) .lt. 0 ) then

                  if( ibon(ii,jj,2) .gt. 0 ) then

                        xc(1) = vx(ii) + xbnd( ibon(ii,jj,2) )
                        xc(2) = vy(jj) + ybnd( ibon(ii,jj,2) )

                  else

                        delx = xdel / 2.0
                        dely = ydel / 2.0
                        ipsm = 1
                        vxx  = vx(ii) + delx
                        vyy  = vy(jj) + dely

                     do k = 1, nbnd

                        call xvwi(vxx,vyy,x(1),x(2),x(3),u,v,w)
                        call smats(1,iout,m,nsb,nl,lt,nlb,ltb,
     &                             nr,mr,kr,ns,ms,ks,
     &                             x,u,v,w,ceye,nmed,isb,ishad)

                        ipsm = 1
                        if( isb .eq. isb0 .and.
     &                      ishad * ishad0 .gt. 0 .and.
     &                      nmed .eq. nmed0 ) ipsm = -1

                        delx = delx / 2.0
                        dely = dely / 2.0
                        vxx  = vxx + ipsm * delx
                        vyy  = vyy + ipsm * dely

                     end do

                        xc(1) = vxx
                        xc(2) = vyy

                        ibon(ii,jj,2) =
     &                  int( ( vxx - vx(ii) ) / delx / 2.0 ) + 1

                  end if

                     if( ibtm .gt. 0 ) then

                        ig = ig + 1
                        nmat(ii,jj,2) = nmat(i,j,1)

                     end if

                        ic = 7
                        goto 310

               end if

            end if

*-----------------------------------------------------------------------
*        ic = 7
*-----------------------------------------------------------------------

  507    ic = 7

                        ii = i
                        jj = j - 1

            if( jj .gt. 0 ) then

               if( inum .eq. abs(imat(ii,jj)) .and.
     &             ishad0 * nmat(ii,jj,1) .gt. 0 ) then

                  if( i + 1 .le. nx ) then

                     do il = i + 1, nx

                        if( inum .eq. abs(imat(il,jj)) .and.
     &                      ishad0 * nmat(il,jj,1) .gt. 0 ) then

                           if( imat(il,jj) .eq. inum ) then

                              imat(il,jj) = -inum

                              if( ibtm .gt. 0 ) then

                                 ig = ig + 1
                                 nmat(il,jj,2) = nmat(il,jj,1)

                              end if

                           end if

                        else

                              if( ibtm .gt. 0 ) then

                                 ig = ig + 1
                                 nmat(il,jj,2) = nmat(il-1,jj,1)

                              end if

                           goto 607

                        end if

                     end do

  607                continue

                  end if

                        i = ii
                        j = jj

                        imat(i,j) = -inum

                        if( ibtm .gt. 0 ) then

                           ig = ig + 1
                           nmat(i,j,2) = nmat(i,j,1)

                        end if

                     if( ibon(i,j,4) .gt. 0 ) then

                        ic = 5

                     else

                        ic = 4

                     end if

                        goto 300

               else

                  if( ibon(ii,jj,3) .gt. 0 ) then

                        xc(1) = vx(ii)
                        xc(2) = vy(jj) + ybnd( ibon(ii,jj,3) )

                  else

                        dely = ydel / 2.0
                        ipsm = 1
                        vxx  = vx(ii)
                        vyy  = vy(jj) + dely

                     do k = 1, nbnd

                        call xvwi(vxx,vyy,x(1),x(2),x(3),u,v,w)
                        call smats(1,iout,m,nsb,nl,lt,nlb,ltb,
     &                             nr,mr,kr,ns,ms,ks,
     &                             x,u,v,w,ceye,nmed,isb,ishad)

                        ipsm = 1
                        if( isb .eq. isb0 .and.
     &                      ishad * ishad0 .gt. 0 .and.
     &                      nmed .eq. nmed0 ) ipsm = -1

                        dely = dely / 2.0
                        vyy  = vyy + ipsm * dely

                     end do

                        xc(1) = vxx
                        xc(2) = vyy

                        ibon(ii,jj,3) =
     &                  int( ( vyy - vy(jj) ) / dely / 2.0 ) + 1

                  end if

                     if( ibtm .gt. 0 ) then

                        ig = ig + 1
                        nmat(ii,jj,2) = nmat(i,j,1)

                     end if

                        ic = 8
                        goto 310

               end if

            else

                        xc(1) = vx(i)
                        xc(2) = vy(j)

                        ic = 1
                        goto 310

            end if

*-----------------------------------------------------------------------
*        ic = 8
*-----------------------------------------------------------------------

  508    ic = 8

                     ii = i + 1
                     jj = j - 1

            if( ii .le. nx .and. jj .gt. 0 ) then

               if( inum .ne. abs(imat(ii,jj)) .or.
     &             ishad0 * nmat(ii,jj,1) .lt. 0 ) then

                  if( ibon(ii,jj,4) .gt. 0 ) then

                        xc(1) = vx(ii) - xbnd( ibon(ii,jj,4) )
                        xc(2) = vy(jj) + ybnd( ibon(ii,jj,4) )

                  else

                        delx = -xdel / 2.0
                        dely =  ydel / 2.0
                        ipsm = 1
                        vxx  = vx(ii) + delx
                        vyy  = vy(jj) + dely

                     do k = 1, nbnd

                        call xvwi(vxx,vyy,x(1),x(2),x(3),u,v,w)
                        call smats(1,iout,m,nsb,nl,lt,nlb,ltb,
     &                             nr,mr,kr,ns,ms,ks,
     &                             x,u,v,w,ceye,nmed,isb,ishad)

                        ipsm = 1
                        if( isb .eq. isb0 .and.
     &                      ishad * ishad0 .gt. 0 .and.
     &                      nmed .eq. nmed0 ) ipsm = -1

                        delx = delx / 2.0
                        dely = dely / 2.0
                        vxx  = vxx + ipsm * delx
                        vyy  = vyy + ipsm * dely

                     end do

                        xc(1) = vxx
                        xc(2) = vyy

                        ibon(ii,jj,4) =
     &                  int( ( vyy - vy(jj) ) / dely / 2.0 ) + 1

                  end if

                     if( ibtm .gt. 0 ) then

                        ig = ig + 1
                        nmat(ii,jj,2) = nmat(i,j,1)

                     end if

                        ic = 1
                        goto 310

               end if

            end if

                  goto 501

*-----------------------------------------------------------------------

  310       continue

                  ip = ip + 1
                  write(ioh) xc(1), xc(2)

            goto 300

*-----------------------------------------------------------------------

  900 continue

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine ppbshow(iot,iog,ioh,nx,ny,nmat,xdel,ydel,xmin,ymin,
     &                   col1,col2,col3,brt,ishdel)
*                                                                      *
*       output bitmap                                                  *
*       last modified by K.Niita on 2002/11/15                         *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      common /shadw/  shadd, ishdd
      dimension   nmat(nx,ny,2)
      dimension   brt(2)

      data eps / 1.0d-08 /

      character huni12*7
      data huni12 /' width='/
      character huni13*7
      data huni13 /' hight='/
      character huni14*16
      character huni15*16

*-----------------------------------------------------------------------
*     internal function
*-----------------------------------------------------------------------

               vx(ik) = xmin + xdel * dble( ik - 1 )
               vy(ik) = ymin + ydel * dble( ik - 1 )

*-----------------------------------------------------------------------

               xdeh = xdel / 2.0
               ydeh = ydel / 2.0

               dst = min(xdeh,ydeh) * eps

               jshdd = nint( dble( ishdel * ishdd ) / shadd )
               if( jshdd .lt. 1 ) jshdd = 1

               rewind iog
               ig = 0

               inum = 0

               goto 110

*-----------------------------------------------------------------------
*     start check each cell
*-----------------------------------------------------------------------

  100    continue

               do j = 1, ny
               do i = 1, nx

                  if( nmat(i,j,2) .lt. 0 ) nmat(i,j,2) = 0

               end do
               end do

  110    continue

               do j = 1, ny
               do i = 1, nx

                  if( nmat(i,j,2) .gt. 0 ) goto 200

               end do
               end do

               goto 900

*-----------------------------------------------------------------------
*           alone cell : first right and up cell check
*-----------------------------------------------------------------------

  200    continue

                  icon = 0

            if( i+2 .le. nx ) then

               if( nmat(i+1,j,2) .gt. 0 .and.
     &             abs( nmat(i+1,j,2) - nmat(i,j,2) ) .lt. jshdd .and.
     &             nmat(i+2,j,2) .gt. 0 .and.
     &             abs( nmat(i+2,j,2) - nmat(i,j,2) ) .lt. jshdd )
     &             icon = 1

            end if

            if( icon .eq. 0 .and. j+2 .le. ny ) then

               if( nmat(i,j+1,2) .gt. 0 .and.
     &             abs( nmat(i,j+1,2) - nmat(i,j,2) ) .lt. jshdd .and.
     &             nmat(i,j+2,2) .gt. 0 .and.
     &             abs( nmat(i,j+2,2) - nmat(i,j,2) ) .lt. jshdd )
     &             icon = 1

            end if

            if( icon .eq. 0 ) then

                  write(iog) i, j, nmat(i,j,2)
                  ig = ig + 1

                  nmat(i,j,2) = 0
                  goto 110

            end if

*-----------------------------------------------------------------------
*        start new path
*-----------------------------------------------------------------------

               i0 = i
               j0 = j

               isavn = 1
               isavr = nmat(i,j,2)
               ismax = nmat(i,j,2) + jshdd / 2
               ismin = nmat(i,j,2) - jshdd / 2

               nmat(i,j,2) = -nmat(i,j,2)

               rewind ioh
               ip = 0

               inum = inum + 1

*-----------------------------------------------------------------------
*           first left and down point
*-----------------------------------------------------------------------

               x1 = vx(i) - xdeh
               y1 = vy(j) + ydeh

               x2 = vx(i) - xdeh
               y2 = vy(j) - ydeh

               x0 = x2
               y0 = y2

               icf = 0

               ic = 1

*-----------------------------------------------------------------------
*        four directions movement and final write information
*-----------------------------------------------------------------------

  300    continue

*-----------------------------------------------------------------------

      if( i .eq. i0 .and. j .eq. j0 .and. ic .eq. 7 ) then

            if( icf .eq. 0 ) then

               x3 = x0
               y3 = y0

               icf = 1

               goto 310

            end if

*-----------------------------------------------------------------------
*              write path: and c(color)
*-----------------------------------------------------------------------

               shad = dble(isavr) / dble(isavn) / shadd

               colx = col1
               coly = col2
               colz = col3

               call givcol(shad,colx,coly,colz,brt)

               write(iot,'("path: c(",3f7.3,")")') colx, coly, colz

               rewind ioh

               do l = 1, ip

                  read(ioh) x1, y1
                  write(iot,'(1p2e14.6)') x1, y1

               end do

*-----------------------------------------------------------------------

            goto 100

      end if

*-----------------------------------------------------------------------
*        four directions movement
*-----------------------------------------------------------------------

  500    continue

               if( ic .eq. 1 ) goto 501
               if( ic .eq. 3 ) goto 503
               if( ic .eq. 5 ) goto 505
               if( ic .eq. 7 ) goto 507

*-----------------------------------------------------------------------
*        ic = 1
*-----------------------------------------------------------------------

  501    ic = 1

                        ii = i + 1
                        jj = j

            if( ii .le. nx ) then

               if( nmat(ii,jj,2) .ne. 0 .and.
     &             abs(nmat(ii,jj,2)) - ismax .gt. -jshdd .and.
     &             abs(nmat(ii,jj,2)) - ismin .lt.  jshdd ) then

                     if( nmat(i,j,2) .gt. 0 ) then

                        isavn = isavn + 1
                        isavr = isavr + nmat(i,j,2)
                        nmat(i,j,2) = -nmat(i,j,2)

                     end if

                  if( j + 1 .le. ny ) then

                     do jl = j + 1, ny

                        if( nmat(ii,jl,2) .gt. 0 .and.
     &                      nmat(ii,jl,2) - ismax .gt. -jshdd .and.
     &                      nmat(ii,jl,2) - ismin .lt.  jshdd ) then

                              isavn = isavn + 1
                              isavr = isavr +  nmat(ii,jl,2)
                              nmat(ii,jl,2) = -nmat(ii,jl,2)

                        else

                           goto 601

                        end if

                     end do

  601                continue

                  end if

                        i = ii
                        j = jj

                        ic = 7
                        goto 300

               else

                        x3 = vx(i) + xdeh
                        y3 = vy(j) - ydeh

                        ic = 3
                        goto 310

               end if

            else

                        x3 = vx(i) + xdeh
                        y3 = vy(j) - ydeh

                        ic = 3
                        goto 310

            end if

*-----------------------------------------------------------------------
*        ic = 3
*-----------------------------------------------------------------------

  503    ic = 3

                        ii = i
                        jj = j + 1

            if( jj .le. ny ) then

               if( nmat(ii,jj,2) .ne. 0 .and.
     &             abs(nmat(ii,jj,2)) - ismax .gt. -jshdd .and.
     &             abs(nmat(ii,jj,2)) - ismin .lt.  jshdd ) then

                     if( nmat(i,j,2) .gt. 0 ) then

                        isavn = isavn + 1
                        isavr = isavr + nmat(i,j,2)
                        nmat(i,j,2) = -nmat(i,j,2)

                     end if

                  if( i - 1 .ge. 1 ) then

                     do il = i - 1, 1, -1

                        if( nmat(il,jj,2) .gt. 0 .and.
     &                      nmat(il,jj,2) - ismax .gt. -jshdd .and.
     &                      nmat(il,jj,2) - ismin .lt.  jshdd ) then

                              isavn = isavn + 1
                              isavr = isavr +  nmat(il,jj,2)
                              nmat(il,jj,2) = -nmat(il,jj,2)

                        else

                           goto 603

                        end if

                     end do

  603                continue

                  end if

                        i = ii
                        j = jj

                        ic = 1
                        goto 300

               else

                        x3 = vx(i) + xdeh
                        y3 = vy(j) + ydeh

                        ic = 5
                        goto 310

               end if

            else

                        x3 = vx(i) + xdeh
                        y3 = vy(j) + ydeh

                        ic = 5
                        goto 310

            end if

*-----------------------------------------------------------------------
*        ic = 5
*-----------------------------------------------------------------------

  505    ic = 5

                        ii = i - 1
                        jj = j

            if( ii .gt. 0 ) then

               if( nmat(ii,jj,2) .ne. 0 .and.
     &             abs(nmat(ii,jj,2)) - ismax .gt. -jshdd .and.
     &             abs(nmat(ii,jj,2)) - ismin .lt.  jshdd ) then

                     if( nmat(i,j,2) .gt. 0 ) then

                        isavn = isavn + 1
                        isavr = isavr + nmat(i,j,2)
                        nmat(i,j,2) = -nmat(i,j,2)

                     end if

                  if( j - 1 .ge. 1 ) then

                     do jl = j - 1, 1, -1

                        if( nmat(ii,jl,2) .gt. 0 .and.
     &                      nmat(ii,jl,2) - ismax .gt. -jshdd .and.
     &                      nmat(ii,jl,2) - ismin .lt.  jshdd ) then

                              isavn = isavn + 1
                              isavr = isavr +  nmat(ii,jl,2)
                              nmat(ii,jl,2) = -nmat(ii,jl,2)

                        else

                           goto 605

                        end if

                     end do

  605                continue

                  end if

                        i = ii
                        j = jj

                        ic = 3
                        goto 300

               else

                        x3 = vx(i) - xdeh
                        y3 = vy(j) + ydeh

                        ic = 7
                        goto 310

               end if

            else

                        x3 = vx(i) - xdeh
                        y3 = vy(j) + ydeh

                        ic = 7
                        goto 310

            end if

*-----------------------------------------------------------------------
*        ic = 7
*-----------------------------------------------------------------------

  507    ic = 7

                        ii = i
                        jj = j - 1

            if( jj .gt. 0 ) then

               if( nmat(ii,jj,2) .ne. 0 .and.
     &             abs(nmat(ii,jj,2)) - ismax .gt. -jshdd .and.
     &             abs(nmat(ii,jj,2)) - ismin .lt.  jshdd ) then

                     if( nmat(i,j,2) .gt. 0 ) then

                        isavn = isavn + 1
                        isavr = isavr + nmat(i,j,2)
                        nmat(i,j,2) = -nmat(i,j,2)

                     end if

                  if( i + 1 .le. nx ) then

                     do il = i + 1, nx

                        if( nmat(il,jj,2) .gt. 0 .and.
     &                      nmat(il,jj,2) - ismax .gt. -jshdd .and.
     &                      nmat(il,jj,2) - ismin .lt.  jshdd ) then

                              isavn = isavn + 1
                              isavr = isavr +  nmat(il,jj,2)
                              nmat(il,jj,2) = -nmat(il,jj,2)

                        else

                           goto 607

                        end if

                     end do

  607                continue

                  end if

                        i = ii
                        j = jj

                        ic = 5
                        goto 300

               else

                        x3 = vx(i) - xdeh
                        y3 = vy(j) - ydeh

                        ic = 1
                        goto 310

               end if

            else

                        x3 = vx(i) - xdeh
                        y3 = vy(j) - ydeh

                        ic = 1
                        goto 310

            end if

*-----------------------------------------------------------------------

  310       continue

            if( ( x2 - x1 )**2 + ( y2 - y1 )**2 .lt. dst ) then

                  x2 = x3
                  y2 = y3

            else if( abs( x2 - x1 ) .lt. dst .and.
     &               abs( x3 - x2 ) .lt. dst ) then

                  x2 = x3
                  y2 = y3

            else if( abs( y2 - y1 ) .lt. dst .and.
     &               abs( y3 - y2 ) .lt. dst ) then

                  x2 = x3
                  y2 = y3

            else

                  write(ioh) x2, y2
                  ip = ip + 1

                  x1 = x2
                  y1 = y2
                  x2 = x3
                  y2 = y3

            end if

            goto 300

*-----------------------------------------------------------------------

  900 continue

*-----------------------------------------------------------------------
*        write alone cell on bmap:
*-----------------------------------------------------------------------

            if( ig .gt. 0 ) then

               rewind iog

                  write(huni14,'(e16.7)') xdel
                  write(huni15,'(e16.7)') ydel

                  write(iot,'(a)') 'bmap: '//
     &                              huni12//huni14//huni13//huni15

               do l = 1, ig

                     read(iog) ii, jj, ishad

                     shad  = dble( abs( ishad ) ) / shadd

                     colx = col1
                     coly = col2
                     colz = col3

                     call givcol(shad,colx,coly,colz,brt)

                     write(iot,'(1p2e14.6,0p3f7.3)')
     &                  vx(ii), vy(jj), colx, coly, colz

               end do

            end if

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      subroutine ppline(ipl,ioh,iot,ip,ddel)
*                                                                      *
*       output boundary clip line                                      *
*       last modified by K.Niita on 2002/11/18                         *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

      parameter ( kc = 4 )

*-----------------------------------------------------------------------

      dimension x(kc), y(kc), xd(kc), yd(kc), drc(kc)

*-----------------------------------------------------------------------

               rewind ioh

*-----------------------------------------------------------------------

         if( ipl .eq. 0 ) then

               do l = 1, ip

                     read(ioh) x0, y0
                     write(iot,'(1p2e16.8)') x0, y0

               end do

*-----------------------------------------------------------------------

         else if( ipl .eq. 1 ) then

                        k = 1

                        x(1) = 0.d0
                        y(1) = 0.d0

               do 100 l = 1, ip

                     read(ioh) x(k), y(k)

                     if( k .gt. 1 ) then

                        dpr = ( x(k) - x(k-1) )**2
     &                      + ( y(k) - y(k-1) )**2

                        if( dpr .lt. ddel ) goto 100

                     end if

                           xd(k) = ( x(k) - x(1) )
                           yd(k) = ( y(k) - y(1) )

                  if( k .ge. 4 ) then

                        ddc  = xd(k)**2 + yd(k)**2
                        icc  = 0

                     if( ddc .gt. 0 ) then

                        dd = ( xd(k)*yd(2) - yd(k)*xd(2) )**2 / ddc
                        if( dd .lt. ddel ) icc = 1

                     end if

                     if( icc .eq. 1 ) then

                              x(2) = x(k)
                              y(2) = y(k)

                              xd(2) = xd(k)
                              yd(2) = yd(k)

                              k = 3

                     else

                        if( k .eq. kc ) then

                              write(iot,'(1p2e16.8)') x(1), y(1)

                           do i = 1, k - 1

                              x(i) = x(i+1)
                              y(i) = y(i+1)

                              xd(i) = x(i+1) - x(1)
                              yd(i) = y(i+1) - y(1)

                           end do

                        else

                              k = k + 1

                        end if

                     end if

                  else

                        k = k + 1

                  end if

  100             continue

                  if( k .gt. 1 ) then

                     do kk = 1, k - 1

                        write(iot,'(1p2e16.8)') x(kk), y(kk)

                     end do

                  end if

*-----------------------------------------------------------------------

         else if( ipl .eq. 2 ) then

                        sddel = sqrt( ddel )

                        k = 1

                        x(1) = 0.d0
                        y(1) = 0.d0

               do 200 l = 1, ip

                     read(ioh) x(k), y(k)

                  if( k .eq. 3 ) then

                           x23 = ( x(3) - x(2) )
                           y23 = ( y(3) - y(2) )
                           r23 = sqrt( x23**2 + y23**2 )

                           if( r23 .lt. sddel ) goto 200

                        rdcos = ( x12 * x23 + y12 * y23 ) / r12 / r23

                     if( rdcos .gt. 0.9999d0 ) then

                           x(2) = x(3)
                           y(2) = y(3)

                           x12 = ( x(2) - x(1) )
                           y12 = ( y(2) - y(1) )
                           r12 = sqrt( x12**2 + y12**2 )

                           if( r12 .lt. sddel ) then

                              k = 2
                              goto 200

                           end if

                     else

                           write(iot,'(1p2e16.8)') x(1), y(1)

                           x(1) = x(2)
                           y(1) = y(2)
                           x(2) = x(3)
                           y(2) = y(3)

                           r12  = r23

                     end if

                  else

                     if( k .eq. 2 ) then

                           x12 = ( x(2) - x(1) )
                           y12 = ( y(2) - y(1) )
                           r12 = sqrt( x12**2 + y12**2 )

                           if( r12 .lt. sddel ) goto 200

                     end if

                        k = k + 1

                  end if

  200             continue

                  if( k .gt. 1 ) then

                     do kk = 1, k - 1

                        write(iot,'(1p2e16.8)') x(kk), y(kk)

                     end do

                  end if

*-----------------------------------------------------------------------

         end if

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine smats(ics,iout,m,nsb,nl,lt,nlb,ltb,
     &                 nr,mr,kr,ns,ms,ks,
     &                 x,u,v,w,ceye,nmed,isb,ishad)
*                                                                      *
*       evaluate one point                                             *
*       last modified by K.Niita on 2002/10/31                         *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'

      parameter ( pi = 3.141592653589793d0 )

*-----------------------------------------------------------------------

      common /tlglat/ ilev1,ilev2,ilat1(5,10),ilat2(5,10)
!$OMP THREADPRIVATE(/tlglat/)
      common /tlcost/ costha, uang(3), nsurf
!$OMP THREADPRIVATE(/tlcost/)

      common /tall40/ rtorg(itlmax,3), rteye(itlmax,3), rtlit(itlmax,3),
     &                rtwin(itlmax,3), rtbrt(itlmax,2), ithvn(itlmax),
     &                itwin(itlmax,2), itbox(itlmax), itmir(itlmax),
     &                rtbox(itlmax,5,10), rtout(itlmax), rthet(itlmax),
     &                itlin(itlmax), itshd(itlmax), itgxs(itlmax)

      common /shadw/  shadd, ishdd

*-----------------------------------------------------------------------

      dimension   lt(nl)
      dimension   ltb(nlb)
      dimension   kr(mr)
      dimension   ks(ms)

      data dsg / 1.0d-08 /

*-----------------------------------------------------------------------

      dimension xo(3), x(3), xc(3), xd(3), angb(3), angd(3), ceye(3)

*-----------------------------------------------------------------------

            shad = 0.5d0

            xo(1) = x(1)
            xo(2) = x(2)
            xo(3) = x(3)

*-----------------------------------------------------------------------
*           check nmed and isb
*-----------------------------------------------------------------------

            call smatsp(iout,m,nsb,nl,lt,nlb,ltb,
     &                  nr,mr,kr,ns,ms,ks,
     &                  x,u,v,w,nmed,isb,angb,isbd)

            if( nmed .le. 0 .or. ics .eq. 0 ) return

*-----------------------------------------------------------------------
*           calculate the angle
*-----------------------------------------------------------------------

                  dx = x(1) - ceye(1)
                  dy = x(2) - ceye(2)
                  dz = x(3) - ceye(3)

                  uvw = sqrt( dx**2 + dy**2 + dz**2 )

                  ux = dx / uvw
                  uy = dy / uvw
                  uz = dz / uvw

            if( isb .gt. 0 ) then

                     angd(1) = uang(1)
                     angd(2) = uang(2)
                     angd(3) = uang(3)

               if( ilev2 .eq. 0 ) then

                     shad = abs( angd(1) * ux
     &                         + angd(2) * uy
     &                         + angd(3) * uz )

               else

                     xd(1) = x(1) + dsg * ux
                     xd(2) = x(2) + dsg * uy
                     xd(3) = x(3) + dsg * uz

                     xc(1) = x(1) - dsg * ux
                     xc(2) = x(2) - dsg * uy
                     xc(3) = x(3) - dsg * uz

                     mark  = 1
                     markp = 0
                     ici  = -1

                  call gomsor(xc(1),xc(2),xc(3),ux,uy,uz,
     &                        nmed0,iblz,mark,markp,ici)

                     markp = 1

                  call gomprp(0,xc(1),xc(2),xc(3),
     &                          xd(1),xd(2),xd(3),
     &                          ux,uy,uz,nmed0,iblz,mark,markp)

                  if( mark .eq. 0 .or. mark .eq. 2 .or.
     &                mark .eq. -1 ) then

                     shad = costha

                  else

                     shad = 0.0001

                  end if

               end if

            else

                     shad = abs( angb(1) * ux
     &                         + angb(2) * uy
     &                         + angb(3) * uz )

            end if

*-----------------------------------------------------------------------
*        calculate the shadow
*-----------------------------------------------------------------------

         if( itshd(m) .ne. 0 ) then

               dist = sqrt( ( x(1) - xo(1) )**2
     &                    + ( x(2) - xo(2) )**2
     &                    + ( x(3) - xo(3) )**2 )

               xc(1) = x(1) - ux * dist
               xc(2) = x(2) - uy * dist
               xc(3) = x(3) - uz * dist

            call smatsp(iout,m,nsb,nl,lt,nlb,ltb,
     &                  nr,mr,kr,ns,ms,ks,
     &                  xc,ux,uy,uz,nmed0,isb0,angb,isbd0)

            if( isb0 .ne. isb .or.
     &        ( itlin(m) .eq. 0 .and. isbd0 .ne. isbd ) .or.
     &          abs( xc(1) - x(1) ) .gt. dsg * 100.0 .or.
     &          abs( xc(2) - x(2) ) .gt. dsg * 100.0 .or.
     &          abs( xc(3) - x(3) ) .gt. dsg * 100.0 ) then

               if( itlin(m) .ne. 0 .and. isb0 .eq. isb ) then

                  shad =  shad / 2.0**itshd(m)

               else if( itlin(m) .eq. 0 .and. isbd0 .eq. isbd ) then

                  shad =  shad / 2.0**itshd(m)

               else

                  shad = -shad / 2.0**itshd(m)

               end if

            end if

         end if

               ishad = int( shad * shadd )
               if( ishad .eq. 0 ) ishad = 1

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine smatsp(iout,m,nsb,nl,lt,nlb,ltb,
     &                  nr,mr,kr,ns,ms,ks,
     &                  x,u,v,w,nmed,isb,angb,isbd0)
*                                                                      *
*       evaluate one point                                             *
*       last modified by K.Niita on 2002/11/08                         *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'

      parameter ( pi = 3.141592653589793d0 )

*-----------------------------------------------------------------------

      common /tall22/ itmcn(itlmax), itmtn(itlmax), itmtt(itlmax)

      common /kmat1d/ idmn(0:kvlmax), idnm(kvmmax)

      common /tlgeom/ iblz1,iblz2
!$OMP THREADPRIVATE(/tlgeom/)
      common /tlglat/ ilev1,ilev2,ilat1(5,10),ilat2(5,10)
!$OMP THREADPRIVATE(/tlglat/)
      common /tlcost/ costha, uang(3), nsurf
!$OMP THREADPRIVATE(/tlcost/)

      common /tall40/ rtorg(itlmax,3), rteye(itlmax,3), rtlit(itlmax,3),
     &                rtwin(itlmax,3), rtbrt(itlmax,2), ithvn(itlmax),
     &                itwin(itlmax,2), itbox(itlmax), itmir(itlmax),
     &                rtbox(itlmax,5,10), rtout(itlmax), rthet(itlmax),
     &                itlin(itlmax), itshd(itlmax), itgxs(itlmax)

*-----------------------------------------------------------------------

      dimension   lt(nl)
      dimension   ltb(nlb)
      dimension   kr(mr)
      dimension   ks(ms)

      data dlg / 1.0d+10 /

*-----------------------------------------------------------------------

      dimension x(3), xc(3), xd(3), angb(3)

*-----------------------------------------------------------------------

                  isb  = 0
                  nmed = -1

                  isbb = 0
                  isbs = 0
                  ijst = 0

*-----------------------------------------------------------------------
*              check position
*-----------------------------------------------------------------------

                  mark  =  1
                  markp =  0
                  ici   = -1

                  call gomsor(x(1),x(2),x(3),u,v,w,
     &                        nmed,iblz,mark,markp,ici)

                  if( nmed .eq. -1 .or. mark .le. -2 ) then

                     nmed = -1
                     goto 500

                  end if

*-----------------------------------------------------------------------

  210 continue

*-----------------------------------------------------------------------
*        check boxes
*-----------------------------------------------------------------------

            call dbox(icp,isbb,nsb,x(1),x(2),x(3),u,v,w,
     &                xd(1),xd(2),xd(3),dist,angb)

*-----------------------------------------------------------------------

         if( icp .eq. 1 ) then

*-----------------------------------------------------------------------

            if( nsb .gt. 0 ) then

*-----------------------------------------------------------------------
*                 check material
*-----------------------------------------------------------------------

  240             continue

                        icmat = 0

                  if( nlb .gt. 0 ) then

                     do k = 1, nlb

                        if( idmn(nmed) .eq. ltb(k) ) then

                           icmat = 1
                           goto 601

                        end if

                     end do

                  end if

*-----------------------------------------------------------------------
*                 check region
*-----------------------------------------------------------------------

  601             continue

                  if( ns .gt. 0 ) then

                     jj  = 0

                     do ii = 1, ns

                        call tregck(iblz2,ilev2,ilat2,ms,ks,jj,icc)

                        if( icc .ne. 0 ) then

                           if( icmat .eq. 1 ) then

                              goto 230

                           else

                              goto 500

                           end if

                        end if

                     end do

                  end if

                  if( icmat .eq. 1 ) goto 500

*-----------------------------------------------------------------------
*                 transparent material in box
*-----------------------------------------------------------------------

  230             continue

                        xc(1) = xd(1)
                        xc(2) = xd(2)
                        xc(3) = xd(3)

                        mark  = 1
                        markp = 1

                     call gomprp(0,x(1),x(2),x(3),xc(1),xc(2),xc(3),
     &                           u,v,w,nmed,iblz,mark,markp)

                        x(1) = xc(1)
                        x(2) = xc(2)
                        x(3) = xc(3)

               if( nmed .eq. -1 .or. mark .le. -1 ) then

                        nmed = -1
                        goto 500

               else if( mark .eq. 0 .or. mark .eq. 2 ) then

                  if( nmed .eq. 0 ) goto 230

                     ijst = 1
                     isbs = nsurf

                     goto 240

               end if

            end if

*-----------------------------------------------------------------------

                     x(1) = xd(1)
                     x(2) = xd(2)
                     x(3) = xd(3)

                     mark  =  1
                     markp =  0
                     ici   = -1

                  call gomsor(x(1),x(2),x(3),u,v,w,
     &                        nmed,iblz,mark,markp,ici)

                  if( nmed .eq. -1 .or. mark .le. -2 ) then

                     nmed = -1
                     goto 500

                  end if

                     ijst = 0

         end if

*-----------------------------------------------------------------------

  200 continue

*-----------------------------------------------------------------------
*              check material
*-----------------------------------------------------------------------

                  icr   = 1
                  icmat = 0

               if( nmed .eq. 0 ) then

                  icr = 0

               else

                  if( nl .gt. 0 ) then

                     do k = 1, nl

                        if( idmn(nmed) .eq. lt(k) ) then

                              icmat = 1

                           if( itmcn(m) .gt. 0 ) then

                              goto 501

                           else

                              icr = 0
                              goto 501

                           end if

                        end if

                     end do

                     if( itmcn(m) .gt. 0 ) icr = 0

  501                continue

                  else

                     icmat = 0

                  end if

               end if

*-----------------------------------------------------------------------
*              check region
*-----------------------------------------------------------------------

  220 continue

               if( ( icr .gt. 0 .or. icmat .eq. 1 ) .and.
     &               nr .gt. 0 ) then

                  jj  = 0

                  do ii = 1, nr

                     call tregck(iblz2,ilev2,ilat2,mr,kr,jj,icc)

                     if( icc .ne. 0 ) then

                        if( icr .eq. 0 .and. icmat .eq. 1 ) then

                           icr = 1
                           goto 503

                        else

                           icr = 0
                           goto 503

                        end if

                     end if

                  end do

  503                continue

               end if

*-----------------------------------------------------------------------
*           transparent material or region, go ahead
*-----------------------------------------------------------------------

            if( icr .eq. 0 ) then

                  if( icp .eq. 0 ) then

                           xc(1) = x(1) + dlg * u
                           xc(2) = x(2) + dlg * v
                           xc(3) = x(3) + dlg * w

                           mark  = 1
                           markp = 1

                        call gomprp(0,x(1),x(2),x(3),xc(1),xc(2),xc(3),
     &                              u,v,w,nmed,iblz,mark,markp)

                     if( mark .le. -1 .or. mark .eq. 1 ) then

                           nmed = -1

                     else

                           isbs = nsurf
                           ijst = 1

                           x(1) = xc(1)
                           x(2) = xc(2)
                           x(3) = xc(3)

                           goto 200

                     end if

*-----------------------------------------------------------------------

                  else if( icp .eq. 1 ) then

                           goto 210

*-----------------------------------------------------------------------

                  else if( icp .eq. 2 ) then

                           mark  = 1
                           markp = 1

                        call gomdis(dpr,x(1),x(2),x(3),
     &                              u,v,w,mark,markp,nmed,iblz)

                     if( dpr .le. dist ) then

                        call gomnew(0,dpr,
     &                              x(1),x(2),x(3),xc(1),xc(2),xc(3),
     &                              u,v,w,ec,nmed,iblz,mark,markp)

                        if( mark .le. -1 ) then

                           nmed = -1

                        else

                           isbs = nsurf
                           ijst = 1

                           x(1) = xc(1)
                           x(2) = xc(2)
                           x(3) = xc(3)

                           goto 210

                        end if

                     else

                           x(1) = xd(1)
                           x(2) = xd(2)
                           x(3) = xd(3)

                           mark  =  1
                           markp =  0
                           ici   = -1

                           call gomsor(x(1),x(2),x(3),u,v,w,
     &                                 nmed,iblz,mark,markp,ici)

                           if( nmed .eq. -1 .or. mark .le. -2 ) then

                              nmed = -1

                           else

                              ijst = 0
                              goto 210

                           end if

                     end if

                  end if

            end if

*-----------------------------------------------------------------------

  500    continue

*-----------------------------------------------------------------------

            if( nmed .gt. 0 ) then

                        isbd  = 0
                        isbd0 = iblz2 * 1000

                     do k = 1, ilev2

                        isbd0 = isbd0 + ilat2(1,k)*1000
     &                        + ilat2(2,k) * ( ilat2(3,k)*100
     &                                      + ilat2(4,k)*10
     &                                      + ilat2(5,k) )

                     end do

                        isbd0 = iabs( isbd0 )

                  if( itlin(m) .ne. 0 .and.
     &              ( iout .eq. 1 .or. iout .eq. 3 ) ) then

                        isbd = isbd0

                  end if

                  if( ijst .eq. 0 ) then

                        isb = isbb * 1000 - isbd

                  else

                        isb = isbs * 1000 + isbd

                  end if

            end if

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine givcol(shad,colx,coly,colz,brt)
*                                                                      *
*       give color from shad and brt                                   *
*       last modified by K.Niita on 2002/11/15                         *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      dimension brt(2)

*-----------------------------------------------------------------------

                  if( colx .gt. 0 ) then

                     if( shad .gt. 0.5 ) then

                        brii = coly
     &                       * min( 1.0d0, max( 0.d0, 1.0d0 - brt(1) ) )
                        brid = coly - brii
                        coly = coly - ( shad - 0.5 ) * 2.0 * brid

                     else

                        brii = colz * min( 1.0d0, max( 0.d0, brt(2) ) )
                        brid = colz - brii
                        colz = colz - ( 0.5 - shad ) * 2.0 * brid

                     end if

                  else

                        colc = -colx

                     if( shad .gt. 0.5 ) then

                        brii = ( 1.0 - colc )
     &                       * min( 1.0d0, max( 0.d0, brt(1) ) )
                        colc = colc + ( shad - 0.5 ) * 2.0 * brii

                     else

                        brii = colc * min( 1.0d0, max( 0.d0, brt(2) ) )
                        brid = colc - brii
                        colc = colc - ( 0.5 - shad ) * 2.0 * brid

                     end if

                        colx = -colc

                  end if

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine xvwi(xx,yy,x0,y0,z0,u,v,w)
*                                                                      *
*       give  the inital coordinate and vector                         *
*       last modified by K.Niita on 2002/10/30                         *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      common /tboxs/ xwin1, ywin1, zwin1, xeye(3),
     &               unx, uny, unz, umx, umy, umz

*-----------------------------------------------------------------------

*-----------------------------------------------------------------------
*     initial position and vector
*-----------------------------------------------------------------------

            x0 = xwin1 + xx * unx + yy * umx
            y0 = ywin1 + xx * uny + yy * umy
            z0 = zwin1 + xx * unz + yy * umz

            dx = x0 - xeye(1)
            dy = y0 - xeye(2)
            dz = z0 - xeye(3)

            uvw = sqrt( dx**2 + dy**2 + dz**2 )

            u = dx / uvw
            v = dy / uvw
            w = dz / uvw

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      subroutine dbox(ic,isb,nsb,x0,y0,z0,u,v,w,
     &                xc,yc,zc,dist,angb)
*                                                                      *
*       check trnsparent boxes                                         *
*       last modified by K.Niita on 2002/10/25                         *
*                                                                      *
*        ic = 0 : no box                                               *
*        ic = 1 : inside the box, xc,yc,zc, the boundary               *
*        ic = 2 : outside the box, dist : next boundary to the box     *
*             nsb=0  xc,yc,zc, the boundary of outside the box         *
*             nsb>0  xc,yc,zc, the boundary of the first box           *
*                                                                      *
*        isb = - ( 10 * ibox + isuf )                                  *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      common /tboxi/ suf(5,6,4), ibox

      dimension insd(5)
      dimension angb(3)

*-----------------------------------------------------------------------
*     no box, pass through ic = 0
*-----------------------------------------------------------------------

            ic  = 0

            xc = x0
            yc = y0
            zc = z0

         if( ibox .eq. 0 ) return

*-----------------------------------------------------------------------

               do k = 1, ibox

                  insd(k) = 1

               end do

*-----------------------------------------------------------------------
*     check initial point is inside of each box
*-----------------------------------------------------------------------

            inss = insdch(x0,y0,z0,u,v,w,ibox,insd,suf)

*-----------------------------------------------------------------------
*     initial point is inside the boxes
*-----------------------------------------------------------------------

         if( inss .gt. 0 ) then

            ic = 1

            dist = dstsuf(x0,y0,z0,u,v,w,xc,yc,zc,
     &                    inss,ibox,insd,suf,isb,angb)

*-----------------------------------------------------------------------
*     initial point is outside the boxes
*-----------------------------------------------------------------------

         else

            dist = dstbox(ic,nsb,x0,y0,z0,u,v,w,xc,yc,zc,
     &                    ibox,insd,suf,isb,angb)

         end if

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      subroutine intbox(jbox,box,ierr)
*                                                                      *
*       initialization of boxes                                        *
*       last modified by K.Niita on 2002/10/30                         *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      common /tboxi/ suf(5,6,4), ibox
      common /tboxj/ idss(5,6)

      dimension box(5,10)
      data dsg / 1.0d-08 /

*-----------------------------------------------------------------------

         ierr = 0
         ibox = jbox
         if( ibox .eq. 0 ) return

         ierr = 1

         do k = 1, ibox

*-----------------------------------------------------------------------

            x0 = box(k,1)
            y0 = box(k,2)
            z0 = box(k,3)

            x1 = box(k,4)
            y1 = box(k,5)
            z1 = box(k,6)

            x2 = box(k,7)
            y2 = box(k,8)
            z2 = box(k,9)

            d3 = box(k,10)

            u1 = ( y1 - y0 ) * ( z2 - z0 ) - ( z1 - z0 ) * ( y2 - y0 )
            v1 = ( z1 - z0 ) * ( x2 - x0 ) - ( x1 - x0 ) * ( z2 - z0 )
            w1 = ( x1 - x0 ) * ( y2 - y0 ) - ( y1 - y0 ) * ( x2 - x0 )

            srt = u1**2 + v1**2 + w1**2
            if( srt .eq. 0.0d0 ) return
            srt = sqrt( srt )

            u1 = u1 / srt
            v1 = v1 / srt
            w1 = w1 / srt

            x3 = x0 + d3 * u1
            y3 = y0 + d3 * v1
            z3 = z0 + d3 * w1

            u2 = ( y3 - y0 ) * ( z1 - z0 ) - ( z3 - z0 ) * ( y1 - y0 )
            v2 = ( z3 - z0 ) * ( x1 - x0 ) - ( x3 - x0 ) * ( z1 - z0 )
            w2 = ( x3 - x0 ) * ( y1 - y0 ) - ( y3 - y0 ) * ( x1 - x0 )

            srt = u2**2 + v2**2 + w2**2
            if( srt .eq. 0.0d0 ) return
            srt = sqrt( srt )

            u2 = u2 / srt
            v2 = v2 / srt
            w2 = w2 / srt

            u3 = ( y2 - y0 ) * ( z3 - z0 ) - ( z2 - z0 ) * ( y3 - y0 )
            v3 = ( z2 - z0 ) * ( x3 - x0 ) - ( x2 - x0 ) * ( z3 - z0 )
            w3 = ( x2 - x0 ) * ( y3 - y0 ) - ( y2 - y0 ) * ( x3 - x0 )

            srt = u3**2 + v3**2 + w3**2
            if( srt .eq. 0.0d0 ) return
            srt = sqrt( srt )

            u3 = u3 / srt
            v3 = v3 / srt
            w3 = w3 / srt

*-----------------------------------------------------------------------

            d1 =  ( u1 * x0 + v1 * y0 + w1 * z0 )
            d2 = -( u1 * x3 + v1 * y3 + w1 * z3 )
            d3 =  ( u2 * x0 + v2 * y0 + w2 * z0 )
            d4 = -( u2 * x2 + v2 * y2 + w2 * z2 )
            d5 =  ( u3 * x0 + v3 * y0 + w3 * z0 )
            d6 = -( u3 * x1 + v3 * y1 + w3 * z1 )

            suf(k,1,1) =  u1
            suf(k,1,2) =  v1
            suf(k,1,3) =  w1
            suf(k,1,4) =  d1

            suf(k,2,1) = -u1
            suf(k,2,2) = -v1
            suf(k,2,3) = -w1
            suf(k,2,4) =  d2

            suf(k,3,1) =  u2
            suf(k,3,2) =  v2
            suf(k,3,3) =  w2
            suf(k,3,4) =  d3

            suf(k,4,1) = -u2
            suf(k,4,2) = -v2
            suf(k,4,3) = -w2
            suf(k,4,4) =  d4

            suf(k,5,1) =  u3
            suf(k,5,2) =  v3
            suf(k,5,3) =  w3
            suf(k,5,4) =  d5

            suf(k,6,1) = -u3
            suf(k,6,2) = -v3
            suf(k,6,3) = -w3
            suf(k,6,4) =  d6

*-----------------------------------------------------------------------

         end do

*-----------------------------------------------------------------------
*     check of identical surface
*-----------------------------------------------------------------------

            do i = 1, ibox

               do j = 1, 6

                  idss(i,j) = -( 10 * i + j )

               end do

            end do

      if( ibox .gt. 1 ) then

         do i = 1, ibox - 1

            do j = 1, 6

               do ii = i + 1, ibox

                  do jj = 1, 6

                        jsuf = -( ii * 10 + jj )

                     if( idss(ii,jj) .eq. jsuf .and.
     &                   abs( abs( suf(i,j,1) + suf(i,j,2)
     &                           + suf(i,j,3) + suf(i,j,4) ) -
     &                        abs( suf(ii,jj,1) + suf(ii,jj,2)
     &                           + suf(ii,jj,3) + suf(ii,jj,4) ) )
     &                  .lt. dsg .and.
     &                   abs( abs(suf(i,j,1)) - abs(suf(ii,jj,1)) )
     &                  .lt. dsg .and.
     &                   abs( abs(suf(i,j,2)) - abs(suf(ii,jj,2)) )
     &                  .lt. dsg .and.
     &                   abs( abs(suf(i,j,3)) - abs(suf(ii,jj,3)) )
     &                  .lt. dsg .and.
     &                   abs( abs(suf(i,j,4)) - abs(suf(ii,jj,4)) )
     &                  .lt. dsg ) then

                           idss(ii,jj) = idss(i,j)

                     end if

                  end do

               end do

            end do

         end do

      end if

*-----------------------------------------------------------------------

      ierr = 0

      return
      end

************************************************************************
*                                                                      *
      function dstbox(ic,nsb,x,y,z,u,v,w,xc,yc,zc,
     &                ibox,insd,suf,isb,angb)
*                                                                      *
*       give the distance up to box or no box                          *
*       last modified by K.Niita on 2002/10/27                         *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

      parameter ( huge = 1d37 )
      parameter ( coincd = 1.d-10 )

      common /tboxj/ idss(5,6)

      dimension insd(5)
      dimension angb(3)
      dimension suf(5,6,4)

      data dsg / 1.0d-08 /

*-----------------------------------------------------------------------

               inss = 0
               dsti = huge
               dstf = huge

         do 200 i = 1, ibox

               dste = huge

            do 210 j = 1, 6

               tu = suf(i,j,1) * u + suf(i,j,2) * v + suf(i,j,3) * w
               if( tu .eq. 0. ) goto 210

               tt = suf(i,j,1) * x + suf(i,j,2) * y + suf(i,j,3) * z
     &            - suf(i,j,4)

               if( abs(tt) .le. coincd * abs(tu) ) tt = tu

               tt = - tt / tu

               if( tt .lt. 0. ) goto 210

                  dst = tt

                  xc = x + dst * u
                  yc = y + dst * v
                  zc = z + dst * w

               do 220 k = 1, 6

                  if( k .eq. j ) goto 220

                  tt = suf(i,k,1)*xc + suf(i,k,2)*yc + suf(i,k,3)*zc
     &               - suf(i,k,4)
                  tu = suf(i,k,1)*u + suf(i,k,2)*v + suf(i,k,3)*w

                  if( abs(tt) .le. coincd * abs(tu) ) tt = tu

                  if( tt .lt. 0.d0 ) goto 210

  220          continue

               if( dst .lt. dste ) then

                  dste = dst
                  ii = i
                  jj = j

               end if

  210       continue

               if( dste .lt. huge .and. dste .lt. dsti ) then

                  dsti = dste
                  inss = i

               end if

  200    continue

               dstbox = dsti

               ic = 0

*-----------------------------------------------------------------------
*        next boundary suface
*-----------------------------------------------------------------------

         if( inss .ne. 0 ) then

            ic = 2

            dsti = dsti + dsg

            xc = x + dsti * u
            yc = y + dsti * v
            zc = z + dsti * w

            isb = idss(ii,jj)

            angb(1) = suf(ii,jj,1)
            angb(2) = suf(ii,jj,2)
            angb(3) = suf(ii,jj,3)

         end if

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      function dstsuf(x,y,z,u,v,w,xc,yc,zc,
     &                in,ibox,insd,suf,isb,angb)
*                                                                      *
*       give the distance up to boundary start from box i              *
*       last modified by K.Niita on 2002/10/27                         *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

      parameter ( huge = 1d37 )
      parameter ( coincd = 1.d-10 )

      common /tboxj/ idss(5,6)

      dimension suf(5,6,4)
      dimension insd(5)
      dimension angb(3)

*-----------------------------------------------------------------------

               i = in

               xc = x
               yc = y
               zc = z

               dist = 0.d0

  100    continue

               ii = i

               dst = huge

            do 210 j = 1, 6

               tu = suf(i,j,1) * u + suf(i,j,2) * v + suf(i,j,3) * w
               if( tu .eq. 0. ) goto 210

               tt = suf(i,j,1) * xc + suf(i,j,2) * yc + suf(i,j,3) * zc
     &            - suf(i,j,4)

               if( abs(tt) .le. coincd * abs(tu) ) tt = tu

               tt = - tt / tu

               if( tt .gt. 0. and. tt .lt. dst ) then

                  dst = tt
                  jj  = j

               end if

  210       continue

               dist = dist + dst

               xc = xc + dst * u
               yc = yc + dst * v
               zc = zc + dst * w

               insd(i) = 0

            i = insdch(xc,yc,zc,u,v,w,ibox,insd,suf)

            if( i .gt. 0 ) goto 100

*-----------------------------------------------------------------------

            dstsuf = dist

            isb = idss(ii,jj)

            angb(1) = suf(ii,jj,1)
            angb(2) = suf(ii,jj,2)
            angb(3) = suf(ii,jj,3)

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
      function insdch(x,y,z,u,v,w,ibox,insd,suf)
*                                                                      *
*       give the first box which includes the point                    *
*       last modified by K.Niita on 2002/10/26                         *
*                                                                      *
************************************************************************

*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      parameter ( coincd = 1.d-10 )

      dimension suf(5,6,4)
      dimension insd(5)

*-----------------------------------------------------------------------

               inss = 0

         do 100 i = 1, ibox

            if( insd(i) .eq. 0 ) goto 100

            do j = 1, 6

               tt = suf(i,j,1) * x + suf(i,j,2) * y + suf(i,j,3) * z
     &            - suf(i,j,4)
               tu = suf(i,j,1) * u + suf(i,j,2) * v + suf(i,j,3) * w

               if( abs(tt) .le. coincd * abs(tu) ) tt = tu

               if( tt .lt. 0.d0 ) goto 100

            end do

               inss = i
               goto 110

  100    continue
  110    continue

*-----------------------------------------------------------------------

         insdch = inss

*-----------------------------------------------------------------------

      return
      end


************************************************************************
*                                                                      *
*  sumover subroutine                                                  *
*                                                                      *
************************************************************************

************************************************************************
*                                                                      *
      subroutine tpdctreg_sumover(m,maxcas,
     &                   np,  ne, nt,  na, nr, tr0x)
*                                                                      *
*     tr_sum(np,ne,nt,na,nr,2)                                           *
*                                                                      *
************************************************************************

      use TALMOD
!$      use TALMOD0

      implicit double precision (a-h,o-z)
      dimension   tr0x(np,ne,nt,na,nr)

      include 'param.inc'
      common /tall10/ itaxn(itlmax), itaxs(itlmax,6)
      common /tallsum/itpan_sum(itlmax,6),itenm_sum(itlmax,6),
     &                  ittnm_sum(itlmax,6),itrgn_sum(itlmax,6),
     &                  itmst_sum(itlmax,6),mftal_sum(itlmax,6),
     &                  itrnm_sum(itlmax,6),itanm_sum(itlmax,6),
     &                  itxnm_sum(itlmax,6),itynm_sum(itlmax,6),
     &                  itznm_sum(itlmax,6),itndy_sum(itlmax,6),
     &                  itndz_sum(itlmax,6),itndn_sum(itlmax,6),
     &                  itndm_sum(itlmax,6),itenm2_sum(itlmax,6),
     &                  itrcn_sum(itlmax,6),itactnm_sum(itlmax,6),
     &                  itmsh_sum(itlmax,6)


C for nonshared_tally option
      integer italsh
      common /talsh/ italsh

      real(8),pointer :: tr_sum(:)

      do iax=1,itaxn(m)

C for nonshared_tally option
!$       if(italsh .eq. 0) then
!$          n_tr_sum = italsize0_sum(m,iax)
!$          call GET_TR_HEAD_POINTER0_SUM(tr_sum,m,iax)
!$       else
            n_tr_sum = italsize_sum(m,iax)
            call GET_TR_HEAD_POINTER_SUM(tr_sum,m,iax)
C for nonshared_tally option
!$       end if
      if(n_tr_sum > 0) then

        call tpdctreg_sumover_sub(maxcas,itaxs(m,iax),
     &     np,  ne, nt, na, nr,    tr0x,
     &     itpan_sum(m,iax),itenm_sum(m,iax),ittnm_sum(m,iax),
     &     itanm_sum(m,iax),itrgn_sum(m,iax),
     &     tr_sum)

        endif

      enddo

      return
      end

************************************************************************
*                                                                      *
      subroutine tpdctreg_sumover_sub(maxcas,itaxs_in,
     &           np,  ne, nt, na, nr, tr0,
     &           np_sum,  ne_sum, nt_sum, na_sum, nr_sum,
     &           tr_sum)
*                                                                      *
************************************************************************

      implicit double precision (a-h,o-z)
      dimension   tr0(np,ne,nt,na,nr)
      dimension   tr_sum(np_sum,ne_sum,nt_sum,na_sum,nr_sum,2)

      if(itaxs_in == 1 .or. itaxs_in == 15) then  ! energ
        do ir = 1,nr
        do ia = 1,na
        do it = 1,nt
        do ip = 1,np
          tr0_sum = 0.0d0
          do ie = 1,ne
            tr0_sum = tr0_sum + tr0(ip,ie,it,ia,ir) / maxcas
          end do
          tr_sum(ip,1,it,ia,ir,1) =
     &           tr_sum(ip,1,it,ia,ir,1) + tr0_sum
          tr_sum(ip,1,it,ia,ir,2) =
     &          tr_sum(ip,1,it,ia,ir,2) + tr0_sum ** 2
        end do
        end do
        end do
        end do
!

      elseif(itaxs_in == 2) then    ! r
        do ia = 1,na
        do it = 1,nt
        do ie = 1,ne
        do ip = 1,np
          tr0_sum = 0.0d0
          do ir = 1,nr
            tr0_sum = tr0_sum + tr0(ip,ie,it,ia,ir) / maxcas
          end do
          tr_sum(ip,ie,it,ia,1,1) =
     &           tr_sum(ip,ie,it,ia,1,1) + tr0_sum
          tr_sum(ip,ie,it,ia,1,2) =
     &           tr_sum(ip,ie,it,ia,1,2) + tr0_sum ** 2
        end do
        end do
        end do
        end do

      elseif(itaxs_in == 11) then   ! time
        do ir = 1,nr
        do ia = 1,na
        do ie = 1,ne
        do ip = 1,np
          tr0_sum = 0.0d0
          do it = 1,nt
            tr0_sum = tr0_sum + tr0(ip,ie,it,ia,ir) / maxcas
          end do
          tr_sum(ip,ie,1,ia,ir,1) =
     &           tr_sum(ip,ie,1,ia,ir,1) + tr0_sum
          tr_sum(ip,ie,1,ia,ir,2) =
     &          tr_sum(ip,ie,1,ia,ir,2) + tr0_sum ** 2
        end do
        end do
        end do
        end do

      elseif(itaxs_in == 12 .or. itaxs_in == 13) then   ! ang
        do ir = 1,nr
        do it = 1,nt
        do ie = 1,ne
        do ip = 1,np
          tr0_sum = 0.0d0
          do ia = 1,na
            tr0_sum = tr0_sum + tr0(ip,ie,it,ia,ir) / maxcas
          end do
          tr_sum(ip,ie,it,1,ir,1) =
     &           tr_sum(ip,ie,it,1,ir,1) + tr0_sum
          tr_sum(ip,ie,it,1,ir,2) =
     &          tr_sum(ip,ie,it,1,ir,2) + tr0_sum ** 2
        end do
        end do
        end do
        end do

       end if
       return
       end



************************************************************************
*                                                                      *
      subroutine tpdctrz_sumover(m,maxcas,
     &                   np, ne, nt, na, nr, nz, tr0x, ldo)
*                                                                      *
*     tr_sum(np,ne,nt,na,nr,nz,2)                                         *
*                                                                      *
************************************************************************

      use TALMOD
!$      use TALMOD0

      implicit double precision (a-h,o-z)
      dimension   tr0x(np,ne,nt,na,nr,nz)

      include 'param.inc'
      common /tall10/ itaxn(itlmax), itaxs(itlmax,6)
      common /tallsum/itpan_sum(itlmax,6),itenm_sum(itlmax,6),
     &                  ittnm_sum(itlmax,6),itrgn_sum(itlmax,6),
     &                  itmst_sum(itlmax,6),mftal_sum(itlmax,6),
     &                  itrnm_sum(itlmax,6),itanm_sum(itlmax,6),
     &                  itxnm_sum(itlmax,6),itynm_sum(itlmax,6),
     &                  itznm_sum(itlmax,6),itndy_sum(itlmax,6),
     &                  itndz_sum(itlmax,6),itndn_sum(itlmax,6),
     &                  itndm_sum(itlmax,6),itenm2_sum(itlmax,6),
     &                  itrcn_sum(itlmax,6),itactnm_sum(itlmax,6),
     &                  itmsh_sum(itlmax,6)

      dimension ldo(2,*)

C for nonshared_tally option
      integer italsh
      common /talsh/ italsh

      real(8),pointer :: tr_sum(:)

      do iax=1,itaxn(m)

C for nonshared_tally option
!$       if(italsh .eq. 0) then
!$          n_tr_sum = italsize0_sum(m,iax)
!$          call GET_TR_HEAD_POINTER0_SUM(tr_sum,m,iax)
!$       else
            n_tr_sum = italsize_sum(m,iax)
            call GET_TR_HEAD_POINTER_SUM(tr_sum,m,iax)
C for nonshared_tally option
!$       end if
      if(n_tr_sum > 0) then

        call tpdctrz_sumover_sub(maxcas,itaxs(m,iax),
     &     np,  ne, nt, na, nr, nz,    tr0x,
     &     itpan_sum(m,iax),itenm_sum(m,iax),ittnm_sum(m,iax),
     &     itanm_sum(m,iax),itrnm_sum(m,iax),itznm_sum(m,iax),
     &     tr_sum,ldo)

        endif

      enddo

      return
      end

************************************************************************
*                                                                      *
      subroutine tpdctrz_sumover_sub(maxcas,itaxs_in,
     &           np,  ne, nt, na, nr, nz, tr0,
     &           np_sum,  ne_sum, nt_sum, na_sum, nr_sum, nz_sum,
     &           tr_sum,ldo)
*                                                                      *
************************************************************************

      implicit double precision (a-h,o-z)
      dimension   tr0(np,ne,nt,na,nr,nz)
      dimension   tr_sum(np_sum,ne_sum,nt_sum,na_sum,nr_sum,nz_sum,2)

      dimension ldo(2,*)

      if(itaxs_in == 1 .or. itaxs_in == 15) then  ! energ
        do iz = ldo(1,6),ldo(2,6)
        do ir = ldo(1,5),ldo(2,5)
        do ia = ldo(1,4),ldo(2,4)
        do it = ldo(1,3),ldo(2,3)
        do ip = ldo(1,1),ldo(2,1)
          tr0_sum = 0.0d0
          do ie = ldo(1,2),ldo(2,2)
            tr0_sum = tr0_sum + tr0(ip,ie,it,ia,ir,iz) / maxcas
          end do
          tr_sum(ip,1,it,ia,ir,iz,1) =
     &           tr_sum(ip,1,it,ia,ir,iz,1) + tr0_sum
          tr_sum(ip,1,it,ia,ir,iz,2) =
     &          tr_sum(ip,1,it,ia,ir,iz,2) + tr0_sum ** 2
        end do
        end do
        end do
        end do
        end do
!

      elseif(itaxs_in == 5) then    ! z
        do ir = ldo(1,5),ldo(2,5)
        do ia = ldo(1,4),ldo(2,4)
        do it = ldo(1,3),ldo(2,3)
        do ie = ldo(1,2),ldo(2,2)
        do ip = ldo(1,1),ldo(2,1)
          tr0_sum = 0.0d0
          do iz = ldo(1,6),ldo(2,6)
            tr0_sum = tr0_sum + tr0(ip,ie,it,ia,ir,iz) / maxcas
          end do
          tr_sum(ip,ie,it,ia,ir,1,1) =
     &           tr_sum(ip,ie,it,ia,ir,1,1) + tr0_sum
          tr_sum(ip,ie,it,ia,ir,1,2) =
     &           tr_sum(ip,ie,it,ia,ir,1,2) + tr0_sum ** 2
        end do
        end do
        end do
        end do
        end do

      elseif(itaxs_in == 6) then    ! r
        do iz = ldo(1,6),ldo(2,6)
        do ia = ldo(1,4),ldo(2,4)
        do it = ldo(1,3),ldo(2,3)
        do ie = ldo(1,2),ldo(2,2)
        do ip = ldo(1,1),ldo(2,1)
          tr0_sum = 0.0d0
          do ir = ldo(1,5),ldo(2,5)
            tr0_sum = tr0_sum + tr0(ip,ie,it,ia,ir,iz) / maxcas
          end do
          tr_sum(ip,ie,it,ia,1,iz,1) =
     &           tr_sum(ip,ie,it,ia,1,iz,1) + tr0_sum
          tr_sum(ip,ie,it,ia,1,iz,2) =
     &           tr_sum(ip,ie,it,ia,1,iz,2) + tr0_sum ** 2
        end do
        end do
        end do
        end do
        end do

      elseif(itaxs_in == 11) then   ! time
        do iz = ldo(1,6),ldo(2,6)
        do ir = ldo(1,5),ldo(2,5)
        do ia = ldo(1,4),ldo(2,4)
        do ie = ldo(1,2),ldo(2,2)
        do ip = ldo(1,1),ldo(2,1)
          tr0_sum = 0.0d0
          do it = ldo(1,3),ldo(2,3)
            tr0_sum = tr0_sum + tr0(ip,ie,it,ia,ir,iz) / maxcas
          end do
          tr_sum(ip,ie,1,ia,ir,iz,1) =
     &           tr_sum(ip,ie,1,ia,ir,iz,1) + tr0_sum
          tr_sum(ip,ie,1,ia,ir,iz,2) =
     &          tr_sum(ip,ie,1,ia,ir,iz,2) + tr0_sum ** 2
        end do
        end do
        end do
        end do
        end do

      elseif(itaxs_in == 12 .or. itaxs_in == 13) then   ! ang
        do iz = ldo(1,6),ldo(2,6)
        do ir = ldo(1,5),ldo(2,5)
        do it = ldo(1,3),ldo(2,3)
        do ie = ldo(1,2),ldo(2,2)
        do ip = ldo(1,1),ldo(2,1)
          tr0_sum = 0.0d0
          do ia = ldo(1,4),ldo(2,4)
            tr0_sum = tr0_sum + tr0(ip,ie,it,ia,ir,iz) / maxcas
          end do
          tr_sum(ip,ie,it,1,ir,iz,1) =
     &           tr_sum(ip,ie,it,1,ir,iz,1) + tr0_sum
          tr_sum(ip,ie,it,1,ir,iz,2) =
     &          tr_sum(ip,ie,it,1,ir,iz,2) + tr0_sum ** 2
        end do
        end do
        end do
        end do
        end do

       end if
       return
       end

************************************************************************
*                                                                      *
      subroutine tpdctxyz_sumover(m,maxcas,
     &                   np,  ne, nt,  na, nx, ny, nz, tr0x, ldo)
*                                                                      *
*     tr_sum(np,ne,nt,na,nx*ny*nz,2)                                   *
*                                                                      *
************************************************************************

      use TALMOD
!$      use TALMOD0

      implicit double precision (a-h,o-z)
      dimension   tr0x(np,ne,nt,na,nx*ny*nz)

      include 'param.inc'
      common /tall10/ itaxn(itlmax), itaxs(itlmax,6)
      common /tallsum/itpan_sum(itlmax,6),itenm_sum(itlmax,6),
     &                  ittnm_sum(itlmax,6),itrgn_sum(itlmax,6),
     &                  itmst_sum(itlmax,6),mftal_sum(itlmax,6),
     &                  itrnm_sum(itlmax,6),itanm_sum(itlmax,6),
     &                  itxnm_sum(itlmax,6),itynm_sum(itlmax,6),
     &                  itznm_sum(itlmax,6),itndy_sum(itlmax,6),
     &                  itndz_sum(itlmax,6),itndn_sum(itlmax,6),
     &                  itndm_sum(itlmax,6),itenm2_sum(itlmax,6),
     &                  itrcn_sum(itlmax,6),itactnm_sum(itlmax,6),
     &                  itmsh_sum(itlmax,6)

      dimension ldo(2,*)

C for nonshared_tally option
      integer italsh
      common /talsh/ italsh

      real(8),pointer :: tr_sum(:)

      do iax=1,itaxn(m)

C for nonshared_tally option
!$       if(italsh .eq. 0) then
!$          n_tr_sum = italsize0_sum(m,iax)
!$          call GET_TR_HEAD_POINTER0_SUM(tr_sum,m,iax)
!$       else
            n_tr_sum = italsize_sum(m,iax)
            call GET_TR_HEAD_POINTER_SUM(tr_sum,m,iax)
C for nonshared_tally option
!$       end if
      if(n_tr_sum > 0) then

        call tpdctxyz_sumover_sub(maxcas,itaxs(m,iax),
     &     np,  ne, nt, na, nx, ny, nz,    tr0x,
     &     itpan_sum(m,iax),itenm_sum(m,iax),ittnm_sum(m,iax),
     &     itanm_sum(m,iax),
     &     itxnm_sum(m,iax),itynm_sum(m,iax),itznm_sum(m,iax),
     &     tr_sum, ldo)

        endif

      enddo

      return
      end


************************************************************************
*                                                                      *
      subroutine tpdctxyz_sumover_sub(maxcas,itaxs_in,
     &         np,  ne, nt, na, nx, ny, nz, tr0,
     &         np_sum,  ne_sum, nt_sum, na_sum, nx_sum, ny_sum, nz_sum,
     &         tr_sum, ldo)
*                                                                      *
************************************************************************

      implicit double precision (a-h,o-z)
      dimension   tr0(np,ne,nt,na,nx*ny*nz)
      dimension   tr_sum(np_sum,ne_sum,nt_sum,na_sum,
     &                   nx_sum*ny_sum*nz_sum,2)
      dimension ldo(2,*)

      icf(ix,iy,iz) = ix + ( iy - 1 ) * nx + ( iz - 1 ) * nx * ny
      icf_sum(ix,iy,iz) = ix + ( iy - 1 ) * nx_sum
     &                  + ( iz - 1 ) * nx_sum * ny_sum

      if(itaxs_in == 1 .or. itaxs_in == 15 ) then  ! energ
        do iz = ldo(1,7),ldo(2,7)
        do iy = ldo(1,6),ldo(2,6)
        do ix = ldo(1,5),ldo(2,5)
        do ia = ldo(1,4),ldo(2,4)
        do it = ldo(1,3),ldo(2,3)
        do ip = ldo(1,1),ldo(2,1)
          tr0_sum = 0.0d0
          do ie = ldo(1,2),ldo(2,2)
            tr0_sum = tr0_sum 
     &              + tr0(ip,ie,it,ia,icf(ix,iy,iz)) / maxcas
          end do
          tr_sum(ip,1,it,ia,icf_sum(ix,iy,iz),1) =
     &           tr_sum(ip,1,it,ia,icf_sum(ix,iy,iz),1) + tr0_sum
          tr_sum(ip,1,it,ia,icf_sum(ix,iy,iz),2) =
     &          tr_sum(ip,1,it,ia,icf_sum(ix,iy,iz),2) + tr0_sum ** 2
        end do
        end do
        end do
        end do
        end do
        end do
!
      elseif(itaxs_in == 3) then    ! x
        do iz = ldo(1,7),ldo(2,7)
        do iy = ldo(1,6),ldo(2,6)
        do ia = ldo(1,4),ldo(2,4)
        do it = ldo(1,3),ldo(2,3)
        do ie = ldo(1,2),ldo(2,2)
        do ip = ldo(1,1),ldo(2,1)
          tr0_sum = 0.0d0
          do ix = ldo(1,5),ldo(2,5)
            tr0_sum = tr0_sum 
     &              + tr0(ip,ie,it,ia,icf(ix,iy,iz)) / maxcas
          end do
          tr_sum(ip,ie,it,ia,icf_sum(1,iy,iz),1) =
     &           tr_sum(ip,ie,it,ia,icf_sum(1,iy,iz),1) + tr0_sum
          tr_sum(ip,ie,it,ia,icf_sum(1,iy,iz),2) =
     &           tr_sum(ip,ie,it,ia,icf_sum(1,iy,iz),2) + tr0_sum ** 2
        end do
        end do
        end do
        end do
        end do
        end do

      elseif(itaxs_in == 4) then    ! y
        do iz = ldo(1,7),ldo(2,7)
        do ix = ldo(1,5),ldo(2,5)
        do ia = ldo(1,4),ldo(2,4)
        do it = ldo(1,3),ldo(2,3)
        do ie = ldo(1,2),ldo(2,2)
        do ip = ldo(1,1),ldo(2,1)
          tr0_sum = 0.0d0
          do iy = ldo(1,6),ldo(2,6)
            tr0_sum = tr0_sum 
     &              + tr0(ip,ie,it,ia,icf(ix,iy,iz)) / maxcas
          end do
          tr_sum(ip,ie,it,ia,icf_sum(ix,1,iz),1) =
     &           tr_sum(ip,ie,it,ia,icf_sum(ix,1,iz),1) + tr0_sum
          tr_sum(ip,ie,it,ia,icf_sum(ix,1,iz),2) =
     &           tr_sum(ip,ie,it,ia,icf_sum(ix,1,iz),2) + tr0_sum ** 2
        end do
        end do
        end do
        end do
        end do
        end do

      elseif(itaxs_in == 5) then    ! z
        do iy = ldo(1,6),ldo(2,6)
        do ix = ldo(1,5),ldo(2,5)
        do ia = ldo(1,4),ldo(2,4)
        do it = ldo(1,3),ldo(2,3)
        do ie = ldo(1,2),ldo(2,2)
        do ip = ldo(1,1),ldo(2,1)
          tr0_sum = 0.0d0
          do iz = ldo(1,7),ldo(2,7)
            tr0_sum = tr0_sum 
     &              + tr0(ip,ie,it,ia,icf(ix,iy,iz)) / maxcas
          end do
          tr_sum(ip,ie,it,ia,icf_sum(ix,iy,1),1) =
     &           tr_sum(ip,ie,it,ia,icf_sum(ix,iy,1),1) + tr0_sum
          tr_sum(ip,ie,it,ia,icf_sum(ix,iy,1),2) =
     &           tr_sum(ip,ie,it,ia,icf_sum(ix,iy,1),2) + tr0_sum ** 2
        end do
        end do
        end do
        end do
        end do
        end do

      elseif(itaxs_in == 11) then   ! time
        do iz = ldo(1,7),ldo(2,7)
        do iy = ldo(1,6),ldo(2,6)
        do ix = ldo(1,5),ldo(2,5)
        do ia = ldo(1,4),ldo(2,4)
        do ie = ldo(1,2),ldo(2,2)
        do ip = ldo(1,1),ldo(2,1)
          tr0_sum = 0.0d0
          do it = ldo(1,3),ldo(2,3)
            tr0_sum = tr0_sum 
     &              + tr0(ip,ie,it,ia,icf(ix,iy,iz)) / maxcas
          end do
          tr_sum(ip,ie,1,ia,icf_sum(ix,iy,iz),1) =
     &           tr_sum(ip,ie,1,ia,icf_sum(ix,iy,iz),1) + tr0_sum
          tr_sum(ip,ie,1,ia,icf_sum(ix,iy,iz),2) =
     &           tr_sum(ip,ie,1,ia,icf_sum(ix,iy,iz),2) + tr0_sum ** 2
        end do
        end do
        end do
        end do
        end do
        end do

      elseif(itaxs_in == 12 .or. itaxs_in == 13) then   ! ang
        do iz = ldo(1,7),ldo(2,7)
        do iy = ldo(1,6),ldo(2,6)
        do ix = ldo(1,5),ldo(2,5)
        do it = ldo(1,3),ldo(2,3)
        do ie = ldo(1,2),ldo(2,2)
        do ip = ldo(1,1),ldo(2,1)
          tr0_sum = 0.0d0
          do ia = ldo(1,4),ldo(2,4)
            tr0_sum = tr0_sum 
     &              + tr0(ip,ie,it,ia,icf(ix,iy,iz)) / maxcas
          end do
          tr_sum(ip,ie,it,1,icf_sum(ix,iy,iz),1) =
     &           tr_sum(ip,ie,it,1,icf_sum(ix,iy,iz),1) + tr0_sum
          tr_sum(ip,ie,it,1,icf_sum(ix,iy,iz),2) =
     &           tr_sum(ip,ie,it,1,icf_sum(ix,iy,iz),2) + tr0_sum ** 2
        end do
        end do
        end do
        end do
        end do
        end do

       end if
       return
       end

************************************************************************
*                                                                      *
      subroutine ppdctreg_sumover_stdev(mode,m,ip,ie,it,ia,ir,
     &                fact_in,ew,tw,aw,vw,ew_sum,tw_sum,aw_sum,vw_sum)
*                                                                      *
************************************************************************

      use TALMOD
!$      use TALMOD0
      use RESTALMOD

      implicit double precision (a-h,o-z)

      include 'param.inc'
      common /tall10/ itaxn(itlmax), itaxs(itlmax,6)
      common /tallsum/itpan_sum(itlmax,6),itenm_sum(itlmax,6),
     &                  ittnm_sum(itlmax,6),itrgn_sum(itlmax,6),
     &                  itmst_sum(itlmax,6),mftal_sum(itlmax,6),
     &                  itrnm_sum(itlmax,6),itanm_sum(itlmax,6),
     &                  itxnm_sum(itlmax,6),itynm_sum(itlmax,6),
     &                  itznm_sum(itlmax,6),itndy_sum(itlmax,6),
     &                  itndz_sum(itlmax,6),itndn_sum(itlmax,6),
     &                  itndm_sum(itlmax,6),itenm2_sum(itlmax,6),
     &                  itrcn_sum(itlmax,6),itactnm_sum(itlmax,6),
     &                  itmsh_sum(itlmax,6)

C for nonshared_tally option
      integer italsh
      common /talsh/ italsh

      real(8),pointer :: tr_sum(:)
     
      do iax=1,itaxn(m)

        if(mode == 0) then
C for nonshared_tally option
!$       if(italsh .eq. 0) then
!$          n_tr_sum = italsize0_sum(m,iax)
!$          call GET_TR_HEAD_POINTER0_SUM(tr_sum,m,iax)
!$       else
            n_tr_sum = italsize_sum(m,iax)
            call GET_TR_HEAD_POINTER_SUM(tr_sum,m,iax)
C for nonshared_tally option
!$       end if
        else
          n_tr_sum = italsize_sum(m,iax)
          tr_sum => trRES_sum(irestalm_sum(m,iax):)
        endif

        if(n_tr_sum > 0) then
          sum_fact = fact_in/ew/tw/aw/vw
          if((itaxs(m,iax) == 1 .or. itaxs(m,iax) == 15) .and.
     &        ie == 1) then
             sum_fact = sum_fact * ew / ew_sum
          else if(itaxs(m,iax) == 2 .and. ir == 1) then
             sum_fact = sum_fact * vw / vw_sum
          else if(itaxs(m,iax) == 11 .and. it == 1) then
             sum_fact = sum_fact * tw / tw_sum
          else if((itaxs(m,iax) == 12 .or. itaxs(m,iax) == 13)
     &            .and. ia == 1) then
             sum_fact = sum_fact * aw / aw_sum
          else
            sum_fact = 0.0d0
          endif
          if(sum_fact /= 0.0d0) then
            call ppdctreg_sumover_stdev_sub(mode,m,sum_fact,
     &         ip,ie,it,ia,ir,
     &         itpan_sum(m,iax),itenm_sum(m,iax),ittnm_sum(m,iax),
     &         itanm_sum(m,iax),itrgn_sum(m,iax),
     &         tr_sum)
          endif

        endif
      enddo

      return
      end

************************************************************************
*                                                                      *
      subroutine ppdctreg_sumover_stdev_ntf(mode,m,ntf,ip,ie,it,ia,ir,
     &                fact_in,ew,tw,aw,vw,ew_sum,tw_sum,aw_sum,vw_sum)
*                                                                      *
************************************************************************

      use TALMOD
!$      use TALMOD0
      use RESTALMOD

      implicit double precision (a-h,o-z)

      include 'param.inc'
      common /tall10/ itaxn(itlmax), itaxs(itlmax,6)
      common /tallsum/itpan_sum(itlmax,6),itenm_sum(itlmax,6),
     &                  ittnm_sum(itlmax,6),itrgn_sum(itlmax,6),
     &                  itmst_sum(itlmax,6),mftal_sum(itlmax,6),
     &                  itrnm_sum(itlmax,6),itanm_sum(itlmax,6),
     &                  itxnm_sum(itlmax,6),itynm_sum(itlmax,6),
     &                  itznm_sum(itlmax,6),itndy_sum(itlmax,6),
     &                  itndz_sum(itlmax,6),itndn_sum(itlmax,6),
     &                  itndm_sum(itlmax,6),itenm2_sum(itlmax,6),
     &                  itrcn_sum(itlmax,6),itactnm_sum(itlmax,6),
     &                  itmsh_sum(itlmax,6)

C for nonshared_tally option
      integer italsh
      common /talsh/ italsh

      real(8),pointer :: tr_sum(:)
     
      do iax=1,itaxn(m)

        if(mode == 0) then
C for nonshared_tally option
!$       if(italsh .eq. 0) then
!$          n_tr_sum = italsize0_sum(m,iax)
!$          call GET_TR_HEAD_POINTER0_SUM_NTF(tr_sum,m,ntf,iax)
!$       else
            n_tr_sum = italsize_sum(m,iax)
            call GET_TR_HEAD_POINTER_SUM_NTF(tr_sum,m,ntf,iax)
C for nonshared_tally option
!$       end if
        else
          ntfbase = mtalsize_sum(m) * (ntf-1)
          n_tr_sum = italsize_sum(m,iax)
          tr_sum => trRES_sum(irestalm_sum(m,iax)+ntfbase:)
        endif

        if(n_tr_sum > 0) then
          sum_fact = fact_in/ew/tw/aw/vw
          if((itaxs(m,iax) == 1 .or. itaxs(m,iax) == 15) .and.
     &        ie == 1) then
             sum_fact = sum_fact * ew / ew_sum
          else if(itaxs(m,iax) == 2 .and. ir == 1) then
             sum_fact = sum_fact * vw / vw_sum
          else if(itaxs(m,iax) == 11 .and. it == 1) then
             sum_fact = sum_fact * tw / tw_sum
          else if((itaxs(m,iax) == 12 .or. itaxs(m,iax) == 13)
     &            .and. ia == 1) then
             sum_fact = sum_fact * aw / aw_sum
          else
            sum_fact = 0.0d0
          endif
          if(sum_fact /= 0.0d0) then
            call ppdctreg_sumover_stdev_sub(mode,m,sum_fact,
     &         ip,ie,it,ia,ir,
     &         itpan_sum(m,iax),itenm_sum(m,iax),ittnm_sum(m,iax),
     &         itanm_sum(m,iax),itrgn_sum(m,iax),
     &         tr_sum)
          endif

        endif
      enddo

      return
      end

************************************************************************
*                                                                      *
      subroutine ppdctreg_sumover_stdev_sub(mode,m,sum_fact,
     &                   ip,ie,it,ia,ir,
     &                   np_sum,ne_sum,nt_sum,na_sum,nr_sum,
     &                   tr_sum)
*                                                                      *
************************************************************************

      implicit double precision (a-h,o-z)

      dimension   tr_sum(np_sum,ne_sum,nt_sum,na_sum,nr_sum,2)
      
      if(tr_sum(ip,ie,it,ia,ir,1) > 0.0) then
       if(mode == 0) then
        call calc_stdev(m,Xa,sigx,
     &                  tr_sum(ip,ie,it,ia,ir,1),
     &                  tr_sum(ip,ie,it,ia,ir,2),
     &                  sum_fact)
       else
        sum_fact_r = 1.0d0 /sum_fact
        call invert_stdev(m,Xa,sigx,
     &                  tr_sum(ip,ie,it,ia,ir,1),
     &                  tr_sum(ip,ie,it,ia,ir,2),
     &                  sum_fact_r)
       endif

        tr_sum(ip,ie,it,ia,ir,1) = Xa
        tr_sum(ip,ie,it,ia,ir,2) = sigx

      else

        tr_sum(ip,ie,it,ia,ir,2) = 0.0

      end if

      return
      end

************************************************************************
*                                                                      *
      subroutine  ppdctreg_sumover_getput(mode,m,iax,
     &    npstepi,nestepi,ntstepi,nastepi,nrstepi,
     &    ipi,iei,iti,iai,iri,maxtott,tott_sum)
*                                                                      *
************************************************************************

      use TALMOD
!$      use TALMOD0
      use RESTALMOD

      implicit double precision (a-h,o-z)

      real(8) :: tott_sum(maxtott,2)

C for nonshared_tally option
      integer italsh
      common /talsh/ italsh

      include 'param.inc'
      common /tall10/ itaxn(itlmax), itaxs(itlmax,6)
      common /tallsum/itpan_sum(itlmax,6),itenm_sum(itlmax,6),
     &                  ittnm_sum(itlmax,6),itrgn_sum(itlmax,6),
     &                  itmst_sum(itlmax,6),mftal_sum(itlmax,6),
     &                  itrnm_sum(itlmax,6),itanm_sum(itlmax,6),
     &                  itxnm_sum(itlmax,6),itynm_sum(itlmax,6),
     &                  itznm_sum(itlmax,6),itndy_sum(itlmax,6),
     &                  itndz_sum(itlmax,6),itndn_sum(itlmax,6),
     &                  itndm_sum(itlmax,6),itenm2_sum(itlmax,6),
     &                  itrcn_sum(itlmax,6),itactnm_sum(itlmax,6),
     &                  itmsh_sum(itlmax,6)

      real(8),pointer :: tr_sum(:)
      real(8),allocatable :: tott_in(:,:)
      character :: chin*200

      if(mode == 0) then
C for nonshared_tally option
!$       if(italsh .eq. 0) then
!$          call GET_TR_HEAD_POINTER0_SUM(tr_sum,m,iax)
!$       else
            call GET_TR_HEAD_POINTER_SUM(tr_sum,m,iax)
C for nonshared_tally option
!$       end if

         call ppdctreg_sumover_getput_sub(mode,
     &    npstepi,nestepi,ntstepi,nastepi,nrstepi,
     &    ipi,iei,iti,iai,iri,maxtott,tott_sum,
     &    itpan_sum(m,iax),itenm_sum(m,iax),ittnm_sum(m,iax),
     &    itanm_sum(m,iax),itrgn_sum(m,iax),
     &    tr_sum)

       else
          tr_sum => trRES_sum(irestalm_sum(m,iax):)

          nsame = npstepi * nestepi * ntstepi * nastepi * nrstepi
          allocate (tott_in(nsame,2))

          read(mode,'(a)') chin
          read(mode,'(26x,1000(1pe13.4,0pf8.4))')
     &    (tott_in(i,1),tott_in(i,2),i=1,nsame)

         call ppdctreg_sumover_getput_sub(mode,
     &    npstepi,nestepi,ntstepi,nastepi,nrstepi,
     &    ipi,iei,iti,iai,iri,nsame,tott_in,
     &    itpan_sum(m,iax),itenm_sum(m,iax),ittnm_sum(m,iax),
     &    itanm_sum(m,iax),itrgn_sum(m,iax),
     &    tr_sum)

          deallocate (tott_in)

      endif

      return
      end

************************************************************************
*                                                                      *
      subroutine ppdctreg_sumover_getput_sub(mode,
     &    npstepi,nestepi,ntstepi,nastepi,nrstepi,
     &    ipi,iei,iti,iai,iri,maxtott,tott_sum,
     &    np_sum,ne_sum,nt_sum,na_sum,nr_sum,tr_sum)
*                                                                      *
************************************************************************

      implicit double precision (a-h,o-z)

      real(8) :: tr_sum(np_sum,ne_sum,nt_sum,na_sum,nr_sum,2)

      real(8) :: tott_sum(maxtott,2)

        i = 0
         do irloop=1,nrstepi
          do ialoop=1,nastepi
            do itloop=1,ntstepi
              do ieloop=1,nestepi
                do iploop=1,npstepi
                  i = i + 1
                  do k=1,2
                   if(mode == 0) then
                    tott_sum(i,k) =
     &              tr_sum(ipi+iploop-1,iei+ieloop-1,iti+itloop-1,
     &              iai+ialoop-1,iri+irloop-1,k)
                   else
                    tr_sum(ipi+iploop-1,iei+ieloop-1,iti+itloop-1,
     &              iai+ialoop-1,iri+irloop-1,k)
     &              = tott_sum(i,k)
                   endif
                  enddo
                end do
              end do
            end do
          end do
        end do

      return
      end

************************************************************************
*                                                                      *
      subroutine ppdctrz_sumover_stdev(mode,m,ip,ie,it,ia,ir,iz,
     &                fact_in,ew,tw,aw,vm,ew_sum,tw_sum,aw_sum,
     &                vl_r,vl_z)
*                                                                      *
************************************************************************

      use TALMOD
!$      use TALMOD0
      use RESTALMOD

      implicit double precision (a-h,o-z)

      include 'param.inc'
      common /tall10/ itaxn(itlmax), itaxs(itlmax,6)
      common /tallsum/itpan_sum(itlmax,6),itenm_sum(itlmax,6),
     &                  ittnm_sum(itlmax,6),itrgn_sum(itlmax,6),
     &                  itmst_sum(itlmax,6),mftal_sum(itlmax,6),
     &                  itrnm_sum(itlmax,6),itanm_sum(itlmax,6),
     &                  itxnm_sum(itlmax,6),itynm_sum(itlmax,6),
     &                  itznm_sum(itlmax,6),itndy_sum(itlmax,6),
     &                  itndz_sum(itlmax,6),itndn_sum(itlmax,6),
     &                  itndm_sum(itlmax,6),itenm2_sum(itlmax,6),
     &                  itrcn_sum(itlmax,6),itactnm_sum(itlmax,6),
     &                  itmsh_sum(itlmax,6)

C for nonshared_tally option
      integer italsh
      common /talsh/ italsh

      real(8),pointer :: tr_sum(:)
     
      do iax=1,itaxn(m)

        if(mode == 0) then
C for nonshared_tally option
!$       if(italsh .eq. 0) then
!$          n_tr_sum = italsize0_sum(m,iax)
!$          call GET_TR_HEAD_POINTER0_SUM(tr_sum,m,iax)
!$       else
            n_tr_sum = italsize_sum(m,iax)
            call GET_TR_HEAD_POINTER_SUM(tr_sum,m,iax)
C for nonshared_tally option
!$       end if
        else
          n_tr_sum = italsize_sum(m,iax)
          tr_sum => trRES_sum(irestalm_sum(m,iax):)
        endif

        if(n_tr_sum > 0) then
          sum_fact = fact_in/vm/ew/tw/aw
          if((itaxs(m,iax) == 1 .or. itaxs(m,iax) == 15) .and.
     &        ie == 1) then
             sum_fact = sum_fact * ew / ew_sum
          else if(itaxs(m,iax) == 5 .and. iz == 1) then
             sum_fact = sum_fact * vm / vl_z
          else if(itaxs(m,iax) == 6 .and. ir == 1) then
             sum_fact = sum_fact * vm / vl_r
          else if(itaxs(m,iax) == 11 .and. it == 1) then
             sum_fact = sum_fact * tw / tw_sum
          else if((itaxs(m,iax) == 12 .or. itaxs(m,iax) == 13) .and.
     &             ia == 1) then
             sum_fact = sum_fact * aw / aw_sum
          else
            sum_fact = 0.0d0
          endif
          if(sum_fact /= 0.0d0) then
            call ppdctrz_sumover_stdev_sub(mode,m,sum_fact,
     &         ip,ie,it,ia,ir,iz,
     &         itpan_sum(m,iax),itenm_sum(m,iax),ittnm_sum(m,iax),
     &         itanm_sum(m,iax),itrnm_sum(m,iax),itznm_sum(m,iax),
     &         tr_sum)
          endif

        endif
      enddo

      return
      end

************************************************************************
*                                                                      *
      subroutine ppdctrz_sumover_stdev_ntf(mode,m,ntf,
     &                ip,ie,it,ia,ir,iz,
     &                fact_in,ew,tw,aw,vm,ew_sum,tw_sum,aw_sum,
     &                vl_r,vl_z)
*                                                                      *
************************************************************************

      use TALMOD
!$      use TALMOD0
      use RESTALMOD

      implicit double precision (a-h,o-z)

      include 'param.inc'
      common /tall10/ itaxn(itlmax), itaxs(itlmax,6)
      common /tallsum/itpan_sum(itlmax,6),itenm_sum(itlmax,6),
     &                  ittnm_sum(itlmax,6),itrgn_sum(itlmax,6),
     &                  itmst_sum(itlmax,6),mftal_sum(itlmax,6),
     &                  itrnm_sum(itlmax,6),itanm_sum(itlmax,6),
     &                  itxnm_sum(itlmax,6),itynm_sum(itlmax,6),
     &                  itznm_sum(itlmax,6),itndy_sum(itlmax,6),
     &                  itndz_sum(itlmax,6),itndn_sum(itlmax,6),
     &                  itndm_sum(itlmax,6),itenm2_sum(itlmax,6),
     &                  itrcn_sum(itlmax,6),itactnm_sum(itlmax,6),
     &                  itmsh_sum(itlmax,6)

C for nonshared_tally option
      integer italsh
      common /talsh/ italsh

      real(8),pointer :: tr_sum(:)
     
      do iax=1,itaxn(m)

        if(mode == 0) then
C for nonshared_tally option
!$       if(italsh .eq. 0) then
!$          n_tr_sum = italsize0_sum(m,iax)
!$          call GET_TR_HEAD_POINTER0_SUM_NTF(tr_sum,m,ntf,iax)
!$       else
            n_tr_sum = italsize_sum(m,iax)
            call GET_TR_HEAD_POINTER_SUM_NTF(tr_sum,m,ntf,iax)
C for nonshared_tally option
!$       end if
        else
          ntfbase = mtalsize_sum(m) * (ntf-1)
          n_tr_sum = italsize_sum(m,iax)
          tr_sum => trRES_sum(irestalm_sum(m,iax)+ntfbase:)
        endif

        if(n_tr_sum > 0) then
          sum_fact = fact_in/vm/ew/tw/aw
          if((itaxs(m,iax) == 1 .or. itaxs(m,iax) == 15) .and.
     &        ie == 1) then
             sum_fact = sum_fact * ew / ew_sum
          else if(itaxs(m,iax) == 5 .and. iz == 1) then
             sum_fact = sum_fact * vm / vl_z
          else if(itaxs(m,iax) == 6 .and. ir == 1) then
             sum_fact = sum_fact * vm / vl_r
          else if(itaxs(m,iax) == 11 .and. it == 1) then
             sum_fact = sum_fact * tw / tw_sum
          else if((itaxs(m,iax) == 12 .or. itaxs(m,iax) == 13) .and.
     &             ia == 1) then
             sum_fact = sum_fact * aw / aw_sum
          else
            sum_fact = 0.0d0
          endif
          if(sum_fact /= 0.0d0) then
            call ppdctrz_sumover_stdev_sub(mode,m,sum_fact,
     &         ip,ie,it,ia,ir,iz,
     &         itpan_sum(m,iax),itenm_sum(m,iax),ittnm_sum(m,iax),
     &         itanm_sum(m,iax),itrnm_sum(m,iax),itznm_sum(m,iax),
     &         tr_sum)
          endif

        endif
      enddo

      return
      end

************************************************************************
*                                                                      *
      subroutine ppdctrz_sumover_stdev_sub(mode,m,sum_fact,
     &                   ip,ie,it,ia,ir,iz,
     &                   np_sum,ne_sum,nt_sum,na_sum,nr_sum,nz_sum,
     &                   tr_sum)
*                                                                      *
************************************************************************

      implicit double precision (a-h,o-z)

      dimension   tr_sum(np_sum,ne_sum,nt_sum,na_sum,nr_sum,nz_sum,2)

      if(tr_sum(ip,ie,it,ia,ir,iz,1) > 0.0) then
        if(mode == 0) then
          call calc_stdev(m,Xa,sigx,
     &                  tr_sum(ip,ie,it,ia,ir,iz,1),
     &                  tr_sum(ip,ie,it,ia,ir,iz,2),
     &                  sum_fact)
        else
          sum_fact_r = 1.0d0 / sum_fact
          call invert_stdev(m,Xa,sigx,
     &                  tr_sum(ip,ie,it,ia,ir,iz,1),
     &                  tr_sum(ip,ie,it,ia,ir,iz,2),
     &                  sum_fact_r)
        endif
        tr_sum(ip,ie,it,ia,ir,iz,1) = Xa
        tr_sum(ip,ie,it,ia,ir,iz,2) = sigx

      else

        tr_sum(ip,ie,it,ia,ir,iz,2) = 0.0

      end if

      return
      end

************************************************************************
*                                                                      *
      subroutine  ppdctrz_sumover_getput(mode,m,iax,
     &    npstepi,nestepi,ntstepi,nastepi,nrstepi,nzstepi,
     &    ipi,iei,iti,iai,iri,izi,maxtott,tott_sum)
*                                                                      *
************************************************************************

      use TALMOD
!$      use TALMOD0
      use RESTALMOD

      implicit double precision (a-h,o-z)

      real(8) :: tott_sum(maxtott,2)

C for nonshared_tally option
      integer italsh
      common /talsh/ italsh

      include 'param.inc'
      common /tall10/ itaxn(itlmax), itaxs(itlmax,6)
      common /tallsum/itpan_sum(itlmax,6),itenm_sum(itlmax,6),
     &                  ittnm_sum(itlmax,6),itrgn_sum(itlmax,6),
     &                  itmst_sum(itlmax,6),mftal_sum(itlmax,6),
     &                  itrnm_sum(itlmax,6),itanm_sum(itlmax,6),
     &                  itxnm_sum(itlmax,6),itynm_sum(itlmax,6),
     &                  itznm_sum(itlmax,6),itndy_sum(itlmax,6),
     &                  itndz_sum(itlmax,6),itndn_sum(itlmax,6),
     &                  itndm_sum(itlmax,6),itenm2_sum(itlmax,6),
     &                  itrcn_sum(itlmax,6),itactnm_sum(itlmax,6),
     &                  itmsh_sum(itlmax,6)

      real(8),pointer :: tr_sum(:)
      real(8),allocatable :: tott_in(:,:)
      character :: chin*200

      if(mode == 0) then
C for nonshared_tally option
!$       if(italsh .eq. 0) then
!$          call GET_TR_HEAD_POINTER0_SUM(tr_sum,m,iax)
!$       else
            call GET_TR_HEAD_POINTER_SUM(tr_sum,m,iax)
C for nonshared_tally option
!$       end if

      call ppdctrz_sumover_getput_sub(mode,
     &    npstepi,nestepi,ntstepi,nastepi,nrstepi,nzstepi,
     &    ipi,iei,iti,iai,iri,izi,maxtott,tott_sum,
     &    itpan_sum(m,iax),itenm_sum(m,iax),ittnm_sum(m,iax),
     &    itanm_sum(m,iax),itrnm_sum(m,iax),itznm_sum(m,iax),
     &    tr_sum)

       else
          tr_sum => trRES_sum(irestalm_sum(m,iax):)

          nsame = npstepi * nestepi * ntstepi * nastepi
     &          * nrstepi * nzstepi
          allocate (tott_in(nsame,2))

          read(mode,'(a)') chin
          read(mode,'(26x,1000(1pe13.4,0pf8.4))')
     &    (tott_in(i,1),tott_in(i,2),i=1,nsame)

          call ppdctrz_sumover_getput_sub(mode,
     &    npstepi,nestepi,ntstepi,nastepi,nrstepi,nzstepi,
     &    ipi,iei,iti,iai,iri,izi,nsame,tott_in,
     &    itpan_sum(m,iax),itenm_sum(m,iax),ittnm_sum(m,iax),
     &    itanm_sum(m,iax),itrnm_sum(m,iax),itznm_sum(m,iax),
     &    tr_sum)

          deallocate (tott_in)

      endif

      return
      end

************************************************************************
*                                                                      *
      subroutine ppdctrz_sumover_getput_sub(mode,
     &    npstepi,nestepi,ntstepi,nastepi,nrstepi,nzstepi,
     &    ipi,iei,iti,iai,iri,izi,maxtott,tott_sum,
     &    np_sum,ne_sum,nt_sum,na_sum,nr_sum,nz_sum,tr_sum)
*                                                                      *
************************************************************************

      implicit double precision (a-h,o-z)

      real(8) :: tr_sum(np_sum,ne_sum,nt_sum,na_sum,nr_sum,nz_sum,2)

      real(8) :: tott_sum(maxtott,2)

        i = 0
         do izloop=1,nzstepi
           do irloop=1,nrstepi
            do ialoop=1,nastepi
             do itloop=1,ntstepi
              do ieloop=1,nestepi
                do iploop=1,npstepi
                  i = i + 1
                  do k=1,2
                   if(mode == 0) then
                    tott_sum(i,k) =
     &              tr_sum(ipi+iploop-1,iei+ieloop-1,iti+itloop-1,
     &              iai+ialoop-1,iri+irloop-1,izi+izloop-1,k)
                   else
                    tr_sum(ipi+iploop-1,iei+ieloop-1,iti+itloop-1,
     &              iai+ialoop-1,iri+irloop-1,izi+izloop-1,k)
     &              = tott_sum(i,k)
                   endif
                  enddo
                end do
              end do
            end do
           end do
          end do
        end do

      return
      end


************************************************************************
*                                                                      *
      subroutine ppdctxyz_sumover_stdev(mode,m,ip,ie,it,ia,ix,iy,iz,
     &                fact_in,ew,tw,aw,vw,ew_sum,tw_sum,aw_sum,
     &                vl_x,vl_y,vl_z)
*                                                                      *
************************************************************************

      use TALMOD
!$      use TALMOD0
      use RESTALMOD

      implicit double precision (a-h,o-z)

      include 'param.inc'
      common /tall10/ itaxn(itlmax), itaxs(itlmax,6)
      common /tallsum/itpan_sum(itlmax,6),itenm_sum(itlmax,6),
     &                  ittnm_sum(itlmax,6),itrgn_sum(itlmax,6),
     &                  itmst_sum(itlmax,6),mftal_sum(itlmax,6),
     &                  itrnm_sum(itlmax,6),itanm_sum(itlmax,6),
     &                  itxnm_sum(itlmax,6),itynm_sum(itlmax,6),
     &                  itznm_sum(itlmax,6),itndy_sum(itlmax,6),
     &                  itndz_sum(itlmax,6),itndn_sum(itlmax,6),
     &                  itndm_sum(itlmax,6),itenm2_sum(itlmax,6),
     &                  itrcn_sum(itlmax,6),itactnm_sum(itlmax,6),
     &                  itmsh_sum(itlmax,6)

C for nonshared_tally option
      integer italsh
      common /talsh/ italsh

      real(8),pointer :: tr_sum(:)
     
      do iax=1,itaxn(m)

        if(mode == 0) then
C for nonshared_tally option
!$       if(italsh .eq. 0) then
!$          n_tr_sum = italsize0_sum(m,iax)
!$          call GET_TR_HEAD_POINTER0_SUM(tr_sum,m,iax)
!$       else
            n_tr_sum = italsize_sum(m,iax)
            call GET_TR_HEAD_POINTER_SUM(tr_sum,m,iax)
C for nonshared_tally option
!$       end if
        else
          n_tr_sum = italsize_sum(m,iax)
          tr_sum => trRES_sum(irestalm_sum(m,iax):)
        endif

        if(n_tr_sum > 0) then
           sum_fact = fact_in/vw/ew/tw/aw
          if((itaxs(m,iax) == 1 .or. itaxs(m,iax) == 15) .and.
     &        ie == 1) then
             sum_fact = sum_fact * ew / ew_sum
          else if(itaxs(m,iax) == 3 .and. ix == 1) then
             sum_fact = sum_fact * vw / vl_x
          else if(itaxs(m,iax) == 4 .and. iy == 1) then
             sum_fact = sum_fact * vw / vl_y
          else if(itaxs(m,iax) == 5 .and. iz == 1) then
             sum_fact = sum_fact * vw / vl_z
          else if(itaxs(m,iax) == 11 .and. it == 1) then
             sum_fact = sum_fact * tw / tw_sum
          else if((itaxs(m,iax) == 12  .or. itaxs(m,iax) == 13) .and.
     &             ia == 1) then
             sum_fact = sum_fact * aw / aw_sum
          else
            sum_fact = 0.0d0
          endif
          if(sum_fact /= 0.0d0) then
            call ppdctxyz_sumover_stdev_sub(mode,m,sum_fact,
     &         ip,ie,it,ia,ix,iy,iz,
     &         itpan_sum(m,iax),itenm_sum(m,iax),ittnm_sum(m,iax),
     &         itanm_sum(m,iax),
     &         itxnm_sum(m,iax),itynm_sum(m,iax),itznm_sum(m,iax),
     &         tr_sum)
          endif

        endif
      enddo

      return
      end

************************************************************************
*                                                                      *
      subroutine ppdctxyz_sumover_stdev_ntf(mode,m,ntf,
     &                ip,ie,it,ia,ix,iy,iz,
     &                fact_in,ew,tw,aw,vw,ew_sum,tw_sum,aw_sum,
     &                vl_x,vl_y,vl_z)
*                                                                      *
************************************************************************

      use TALMOD
!$      use TALMOD0
      use RESTALMOD

      implicit double precision (a-h,o-z)

      include 'param.inc'
      common /tall10/ itaxn(itlmax), itaxs(itlmax,6)
      common /tallsum/itpan_sum(itlmax,6),itenm_sum(itlmax,6),
     &                  ittnm_sum(itlmax,6),itrgn_sum(itlmax,6),
     &                  itmst_sum(itlmax,6),mftal_sum(itlmax,6),
     &                  itrnm_sum(itlmax,6),itanm_sum(itlmax,6),
     &                  itxnm_sum(itlmax,6),itynm_sum(itlmax,6),
     &                  itznm_sum(itlmax,6),itndy_sum(itlmax,6),
     &                  itndz_sum(itlmax,6),itndn_sum(itlmax,6),
     &                  itndm_sum(itlmax,6),itenm2_sum(itlmax,6),
     &                  itrcn_sum(itlmax,6),itactnm_sum(itlmax,6),
     &                  itmsh_sum(itlmax,6)

C for nonshared_tally option
      integer italsh
      common /talsh/ italsh

      real(8),pointer :: tr_sum(:)
     
      do iax=1,itaxn(m)

        if(mode == 0) then
C for nonshared_tally option
!$       if(italsh .eq. 0) then
!$          n_tr_sum = italsize0_sum(m,iax)
!$          call GET_TR_HEAD_POINTER0_SUM_NTF(tr_sum,m,ntf,iax)
!$       else
            n_tr_sum = italsize_sum(m,iax)
            call GET_TR_HEAD_POINTER_SUM_NTF(tr_sum,m,ntf,iax)
C for nonshared_tally option
!$       end if
        else
          ntfbase = mtalsize_sum(m) * (ntf-1)
          n_tr_sum = italsize_sum(m,iax)
          tr_sum => trRES_sum(irestalm_sum(m,iax)+ntfbase:)
        endif

        if(n_tr_sum > 0) then
           sum_fact = fact_in/vw/ew/tw/aw
          if((itaxs(m,iax) == 1 .or. itaxs(m,iax) == 15) .and.
     &        ie == 1) then
             sum_fact = sum_fact * ew / ew_sum
          else if(itaxs(m,iax) == 3 .and. ix == 1) then
             sum_fact = sum_fact * vw / vl_x
          else if(itaxs(m,iax) == 4 .and. iy == 1) then
             sum_fact = sum_fact * vw / vl_y
          else if(itaxs(m,iax) == 5 .and. iz == 1) then
             sum_fact = sum_fact * vw / vl_z
          else if(itaxs(m,iax) == 11 .and. it == 1) then
             sum_fact = sum_fact * tw / tw_sum
          else if((itaxs(m,iax) == 12  .or. itaxs(m,iax) == 13) .and.
     &             ia == 1) then
             sum_fact = sum_fact * aw / aw_sum
          else
            sum_fact = 0.0d0
          endif
          if(sum_fact /= 0.0d0) then
            call ppdctxyz_sumover_stdev_sub(mode,m,sum_fact,
     &         ip,ie,it,ia,ix,iy,iz,
     &         itpan_sum(m,iax),itenm_sum(m,iax),ittnm_sum(m,iax),
     &         itanm_sum(m,iax),
     &         itxnm_sum(m,iax),itynm_sum(m,iax),itznm_sum(m,iax),
     &         tr_sum)
          endif

        endif
      enddo

      return
      end

************************************************************************
*                                                                      *
      subroutine ppdctxyz_sumover_stdev_sub(mode,m,sum_fact,
     &                   ip,ie,it,ia,ix,iy,iz,
     &                   np_sum,ne_sum,nt_sum,na_sum,
     &                   nx_sum,ny_sum,nz_sum,tr_sum)
*                                                                      *
************************************************************************

      implicit double precision (a-h,o-z)

      dimension   tr_sum(np_sum,ne_sum,nt_sum,na_sum,
     &                   nx_sum*ny_sum*nz_sum,2)
      
      icf_sum(ix,iy,iz) = ix + ( iy - 1 ) * nx_sum
     &                  + ( iz - 1 ) * nx_sum * ny_sum


      if(tr_sum(ip,ie,it,ia,icf_sum(ix,iy,iz),1) > 0.0) then
        if(mode == 0) then
          call calc_stdev(m,Xa,sigx,
     &                  tr_sum(ip,ie,it,ia,icf_sum(ix,iy,iz),1),
     &                  tr_sum(ip,ie,it,ia,icf_sum(ix,iy,iz),2),
     &                  sum_fact)

        else
          sum_fact_r = 1.0d0 / sum_fact
          call invert_stdev(m,Xa,sigx,
     &                  tr_sum(ip,ie,it,ia,icf_sum(ix,iy,iz),1),
     &                  tr_sum(ip,ie,it,ia,icf_sum(ix,iy,iz),2),
     &                  sum_fact_r)

        endif

        tr_sum(ip,ie,it,ia,icf_sum(ix,iy,iz),1) = Xa
        tr_sum(ip,ie,it,ia,icf_sum(ix,iy,iz),2) = sigx

      else

        tr_sum(ip,ie,it,ia,icf_sum(ix,iy,iz),2) = 0.0

      end if

      return
      end

************************************************************************
*                                                                      *
      subroutine  ppdctxyz_sumover_getput(mode,m,iax,
     &    npstepi,nestepi,ntstepi,nastepi,nxstepi,nystepi,nzstepi,
     &    ipi,iei,iti,iai,ixi,iyi,izi,maxtott,tott_sum)
*                                                                      *
************************************************************************

      use TALMOD
!$      use TALMOD0
      use RESTALMOD

      implicit double precision (a-h,o-z)

      real(8) :: tott_sum(maxtott,2)

C for nonshared_tally option
      integer italsh
      common /talsh/ italsh

      include 'param.inc'
      common /tall10/ itaxn(itlmax), itaxs(itlmax,6)
      common /tallsum/itpan_sum(itlmax,6),itenm_sum(itlmax,6),
     &                  ittnm_sum(itlmax,6),itrgn_sum(itlmax,6),
     &                  itmst_sum(itlmax,6),mftal_sum(itlmax,6),
     &                  itrnm_sum(itlmax,6),itanm_sum(itlmax,6),
     &                  itxnm_sum(itlmax,6),itynm_sum(itlmax,6),
     &                  itznm_sum(itlmax,6),itndy_sum(itlmax,6),
     &                  itndz_sum(itlmax,6),itndn_sum(itlmax,6),
     &                  itndm_sum(itlmax,6),itenm2_sum(itlmax,6),
     &                  itrcn_sum(itlmax,6),itactnm_sum(itlmax,6),
     &                  itmsh_sum(itlmax,6)

      real(8),pointer :: tr_sum(:)
      real(8),allocatable :: tott_in(:,:)
      character :: chin*200

      if(mode == 0) then
C for nonshared_tally option
!$       if(italsh .eq. 0) then
!$          call GET_TR_HEAD_POINTER0_SUM(tr_sum,m,iax)
!$       else
            call GET_TR_HEAD_POINTER_SUM(tr_sum,m,iax)
C for nonshared_tally option
!$       end if

          call ppdctxyz_sumover_getput_sub(mode,
     &    npstepi,nestepi,ntstepi,nastepi,nxstepi,nystepi,nzstepi,
     &    ipi,iei,iti,iai,ixi,iyi,izi,maxtott,tott_sum,
     &    itpan_sum(m,iax),itenm_sum(m,iax),ittnm_sum(m,iax),
     &    itanm_sum(m,iax),
     &    itxnm_sum(m,iax),itynm_sum(m,iax),itznm_sum(m,iax),
     &    tr_sum)

       else
          tr_sum => trRES_sum(irestalm_sum(m,iax):)

          nsame = npstepi * nestepi * ntstepi * nastepi
     &          * nxstepi * nystepi * nzstepi
          allocate (tott_in(nsame,2))

          read(mode,'(a)') chin
          read(mode,'(26x,1000(1pe13.4,0pf8.4))')
     &    (tott_in(i,1),tott_in(i,2),i=1,nsame)

          call ppdctxyz_sumover_getput_sub(mode,
     &    npstepi,nestepi,ntstepi,nastepi,nxstepi,nystepi,nzstepi,
     &    ipi,iei,iti,iai,ixi,iyi,izi,nsame,tott_in,
     &    itpan_sum(m,iax),itenm_sum(m,iax),ittnm_sum(m,iax),
     &    itanm_sum(m,iax),
     &    itxnm_sum(m,iax),itynm_sum(m,iax),itznm_sum(m,iax),
     &    tr_sum)

          deallocate (tott_in)

      endif

      return
      end

************************************************************************
*                                                                      *
      subroutine ppdctxyz_sumover_getput_sub(mode,
     &    npstepi,nestepi,ntstepi,nastepi,nxstepi,nystepi,nzstepi,
     &    ipi,iei,iti,iai,ixi,iyi,izi,maxtott,tott_sum,
     &    np_sum,ne_sum,nt_sum,na_sum,nx_sum,ny_sum,nz_sum,tr_sum)
*                                                                      *
************************************************************************

      implicit double precision (a-h,o-z)

      real(8) :: tr_sum(np_sum,ne_sum,nt_sum,na_sum,
     &                  nx_sum*ny_sum*nz_sum,2)

      real(8) :: tott_sum(maxtott,2)

      icf_sum(ix,iy,iz) = ix + ( iy - 1 ) * nx_sum
     &                  + ( iz - 1 ) * nx_sum * ny_sum

        i = 0
         do izloop=1,nzstepi
          do iyloop=1,nystepi
           do ixloop=1,nxstepi
            do ialoop=1,nastepi
             do itloop=1,ntstepi
              do ieloop=1,nestepi
                do iploop=1,npstepi
                  i = i + 1
                  do k=1,2
                   if(mode == 0) then
                    tott_sum(i,k) =
     &              tr_sum(ipi+iploop-1,iei+ieloop-1,iti+itloop-1,
     &              iai+ialoop-1,
     &              icf_sum(ixi+ixloop-1,iyi+iyloop-1,izi+izloop-1),k)
                   else
                    tr_sum(ipi+iploop-1,iei+ieloop-1,iti+itloop-1,
     &              iai+ialoop-1,
     &              icf_sum(ixi+ixloop-1,iyi+iyloop-1,izi+izloop-1),k)
     &              = tott_sum(i,k)
                   endif
                  enddo
                end do
              end do
            end do
           end do
          end do
         end do
        end do

      return
      end
