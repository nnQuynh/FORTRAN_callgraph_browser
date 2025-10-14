************************************************************************
*                                                                      *
      subroutine EVENTanalyz
*                                                                      *
*       analyz for deposit option                                      *
*----------------------------------------------------------------------*
*                                                                      *
************************************************************************
*-----------------------------------------------------------------------

      implicit real*8(a-h,o-z)

*-----------------------------------------------------------------------

      common /mpi00/ npe, me

*-----------------------------------------------------------------------
*        tally for deposit of heat tally : ncol = 0
*-----------------------------------------------------------------------

         if( me .gt. 0 .or. npe .le. 1 ) then

                  call cputime(3)

                  mcol = 0
                  mark_dum = 0
               call analyz(mcol,mark_dum)
                  call cputime(3)

         end if

*-----------------------------------------------------------------------
         return
         end

************************************************************************
*                                                                      *
      subroutine thetregEVENT(ncol,m,nd,ns0,nl,lt,np,nr,mr,kr,
     &                   ne,eb,tr,rabs,trEVENT,tr0)
*                                                                      *
*       nuclear heating tally in region mesh                           *
*       output=deposit-*                                               *
*       last modified by T.Furuta on 2012/05/08                        *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*      ncol  ..... reaction type                                       *
*              0 : end of batch                                        *
*              4 : source particle                                     *
*                                                                      *
*              9 : termination by time cut-off                         *
*             10 : geometry boundary crossing                          *
*             11 : termination by energy cut-off                       *
*             12 : termination by escape or leakage                    *
*             13 : nuclear reaction (n,x)                              *
*             14 : nuclear reaction (n,n'x)                            *
*             15 : sequential transport only for tally                 *
*                                                                      *
*     definition of tr( nd, ir, id )                                   *
*                                                                      *
*        ns = 1 : all, 2 : simple, 3 : heat                            *
*                                                                      *
*      nd ..... number of event type                                   *
*        ns = 3                                                        *
*          1 : total                                                   *
*                                                                      *
*        ns = 2                                                        *
*          2 : cutoff neutron                                          *
*          3 : cutoff photon                                           *
*          4 : cutoff proton                                           *
*          5 : leakage                                                 *
*                                                                      *
*          6 : recoil                                                  *
*          7 : ionization                                              *
*          8 : low neutron                                             *
*          9 : electron                                                *
*         10 : low proton                                              *
*         11 : other                                                   *
*                                                                      *
*        ns = 1                                                        *
*         12 : deuteron recoil                                         *
*         13 : triton   recoil                                         *
*         14 : 3He      recoil                                         *
*         15 : alpha    recoil                                         *
*         16 : residual recoil                                         *
*                                                                      *
*         17 : proton ionization                                       *
*         18 : pion+  ionization                                       *
*         19 : pion-  ionization                                       *
*         20 : other  ionization                                       *
*                                                                      *
*         21 : stopped proton                                          *
*         22 : stopped neutron                                         *
*         23 : stopped photon                                          *
*         24 : stopped pion+                                           *
*         25 : stopped pion-                                           *
*         26 : stopped other                                           *
*                                                                      *
*         27 : remaining excitation energy                             *
*         28 : fission high recoil                                     *
*         29 : fission low neutron                                     *
*                                                                      *
*         29 + np      : ionization of other paticle                   *
*         29 + np + np : stopped other paticle                         *
*                                                                      *
*       id ..... data type                                             *
*          1 : heating                                                 *
*          2 : relative error                                          *
*                                                                      *
*       unit = MeV                                                     *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
************************************************************************
      use MMBANKMOD   !FURUTA
      use EVENTTALMOD !FURUTA
      use NGSDATAMOD, only : bindeg
      use partmod, only: itmxpt, itpan, itpat, jtpat ! frtati 2021/10/05
*-----------------------------------------------------------------------

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param00.inc'
      include 'param.inc'

*-----------------------------------------------------------------------

      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
!$OMP THREADPRIVATE(/icomon/)
      common /wtsave/ oldwt
!$OMP THREADPRIVATE(/wtsave/)
      common /tlgeom/ iblz1,iblz2
!$OMP THREADPRIVATE(/tlgeom/)
      common /tlglat/ ilev1,ilev2,ilat1(5,10),ilat2(5,10)
!$OMP THREADPRIVATE(/tlglat/)
      common /mathzn/ mathz, mathn, jcoll, kcoll
!$OMP THREADPRIVATE(/mathzn/)
      common /regdc/  idrg(kvlmax), idgr(kvmmax)
      common /kmat1d/ idmn(0:kvlmax), idnm(kvmmax)
      common /eparm/  esmax, esmin, emin(20)
      common /celdn/  denr(kvlmax), denm(kvlmax), denc(kvlmax)

*-----------------------------------------------------------------------

      common /clustt/ nclsts, iclusts(nnn)
!$OMP THREADPRIVATE(/clustt/)
      common /clustw/ jclusts(0:8,nnn),  qclusts(0:12,nnn)
!$OMP THREADPRIVATE(/clustw/)
      common /cntcls/ jcount(3,nnn)
!$OMP THREADPRIVATE(/cntcls/)

*-----------------------------------------------------------------------

      common /tall22/ itmcn(itlmax), itmtn(itlmax), itmtt(itlmax)
      common /tall32/ itelc(itlmax)
      common /tall37/ itcnt(9,itlmax)

*-----------------------------------------------------------------------

      dimension   lt(nl)
      dimension   kr(mr)
      dimension   tr(nd,nr,0:ne,2)
      dimension   rabs(5)
      dimension   eb(ne+1)
      dimension   trEVENT(nd,nr)

      dimension ips(itmxpt) ! frtati 2021/10/05 6 -> itmxpt

*-----------------------------------------------------------------------

      data rmngv / 0.93895 /
      data rmnmv / 938.95 /

      save sweight
!$OMP THREADPRIVATE(sweight)
      common /jcomon/ nabov,nobch,nocas,nomax
!$OMP THREADPRIVATE(/jcomon/)
      common /tall36/ itdpo(itlmax)

*-----------------------------------------------------------------------
      common /dedxfac/ dedxfd
!$OMP THREADPRIVATE(/dedxfac/)

*-----------------------------------------------------------------------

      real*8 edep
      common/EPCONTedep/ edep
!$OMP THREADPRIVATE(/EPCONTedep/)

      integer iegsemi, iegsout
      common /egsemi/ iegsemi, iegsout

C for nonshared_tally option
      integer italsh
      common /talsh/ italsh

! sumover
      real(8) tr0(nd,nr,0:ne)
      tr0(:,:,:) = 0.0d0

*-----------------------------------------------------------------------
*     ns from ns0
*-----------------------------------------------------------------------

            ns = ns0
            if( ns .ge. 4 ) ns = ns - 3

*-----------------------------------------------------------------------
*     deposit : ncol = 0, 4 and ns = 4
*-----------------------------------------------------------------------

      if( ( ncol .eq. 0 .or. ncol .eq. 4 ) .and. ns0 .ge. 4 ) then

         if( .not. FIRSTsrc ) then

            do id = 1, nd
            do ir = 1, nr

               if( itdpo(m) .eq. 0
     &             .and. trEVENT(1,ir) .gt. 0.0 ) then

                  heats = trEVENT(1,ir)
                  ratio = trEVENT(id,ir) / trEVENT(1,ir)

               else

                  heats = trEVENT(id,ir)
                  ratio = 1.0d0

               end if

               if( heats .gt. 0.0d0 .and.
     &             heats .ge. eb(1) .and. heats .lt. eb(ne+1) ) then

                  do ie = 1, ne

                     if( heats .ge. eb(ie) .and.
     &                   heats .lt. eb(ie+1) ) then
! sumover
                        tr0(id,ir,ie) = tr0(id,ir,ie)
     &                                + sweight*ratio
                       if (italsh .eq. 0 ) then
                          tr(id,ir,ie,1) = tr(id,ir,ie,1)
     &                                   + sweight*ratio
                          tr(id,ir,ie,2) = tr(id,ir,ie,2)
     &                                   + (sweight*ratio)**2
                       else
!$OMP CRITICAL (tr_thetregEVENT)
                          tr(id,ir,ie,1) = tr(id,ir,ie,1)
     &                                   + sweight*ratio
                          tr(id,ir,ie,2) = tr(id,ir,ie,2)
     &                                   + (sweight*ratio)**2
!$OMP END CRITICAL (tr_thetregEVENT)
                       end if

                     end if

                  end do

               end if

            end do
            end do

! sumover
           call thetreg_sumover(m, 1,
     &                           nd, nr, ne, tr0)
           tr0(:,:,:) = 0.0d0

            do id = 1, nd
            do ir = 1, nr

               trEVENT(id,ir) = 0.0d0

            end do
            end do

         end if

*-----------------------------------------------------------------------
*        save initial weight and zero set
*-----------------------------------------------------------------------

         if( ncol .eq. 4 ) then

            sweight = oldwt

         end if

      end if

*-----------------------------------------------------------------------
*     check of ncol
*-----------------------------------------------------------------------

      if( ncol .lt. 9 ) return

*-----------------------------------------------------------------------
*        nothing for electron and positron with nuclear data at itelc=0
*                    and secondary photon
*-----------------------------------------------------------------------

         if( itelc(m) .eq. 0 .and. iegsemi .eq. 0 .and.   ! S.Abe 2015/08/31
     &     ( jcoll .eq. 8 .or.
     &     ( jcoll .eq. 7 .and. name(ibknam+no,ipomp+1) .gt. 1 ) ) )
     &                       return

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
*        incident particle energy loss
*-----------------------------------------------------------------------
*        enmin : final particle energy ( with mass for not nucleon )
*        enion : ionization loss energy for charged particles
*-----------------------------------------------------------------------

               if( ityp .le. 2 .or. ityp .ge. 12 ) then
                  enmin = 0.0
               else
                  enmin = rtyp - mtyp * rmnmv
               end if

               enmin = ( ec(ibkec+no,ipomp+1) + enmin ) * oldwt
             enion = abs( ec(ibkec+no,ipomp+1) - e(ibke+no,ipomp+1) ) *
     &                                      oldwt

               enion = enion * dedxfd

*-----------------------------------------------------------------------
cKN Iwase 2014/08/22 for EGS

               if( ( ityp.eq.12 .or. ityp.eq.13 .or. ityp.eq.14 ) .and.
     &               iegsemi .ne. 0 ) then

                     enion = edep * oldwt

                  if(ncol.eq.9 .or. ncol.eq.11) then
                     enmin = enmin + edep * oldwt
                     enion = 0d0
                  endif

               endif

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
*        check of counter
*-----------------------------------------------------------------------

            do i = 1, 3

               if( itcnt(i,m) .ne. 0 ) then

                  if( ncnt(ibknct+i,no,ipomp+1) .lt. itcnt(i*2+2,m) .or.
     &                ncnt(ibknct+i,no,ipomp+1) .gt. itcnt(i*2+3,m) )
     &                      goto 200

               end if

            end do

*-----------------------------------------------------------------------
*        ionization
*-----------------------------------------------------------------------

            if( enion .gt. 0.0 ) then

                        trEVENT(1,ir) = trEVENT(1,ir) + enion

               if( ns .le. 2 ) then

                        trEVENT(7,ir) = trEVENT(7,ir) + enion

               end if

               if( ns .eq. 1 ) then

                  if( ityp .eq. 1 ) then

                        k = 17

                  else if( ityp .eq. 3 ) then

                        k = 18

                  else if( ityp .eq. 5 ) then

                        k = 19

                  else if( ityp .eq. 12 .or. ityp .eq. 13 ) then

                        k = 20
                        rabs(4) = rabs(4) + enion

                  else

                        k = 20

                  end if

                        trEVENT(k,ir) = trEVENT(k,ir) + enion

                  if( np .gt. 0 ) then

                     call pcheck(m,np,ityp,ktyp,jtyp,ipn,ips)

                     if( ipn .gt. 0 ) then

                        do i = 1, ipn

                           k = 29 + ips(i)

                           trEVENT(k,ir) = trEVENT(k,ir) + enion

                        end do

                     end if

                  end if

               end if

            end if

*-----------------------------------------------------------------------
*        stopped energy cut-off or time cut-off particles
*-----------------------------------------------------------------------

            if( ncol .eq. 9 .or. ncol .eq. 11 ) then

                        trEVENT(1,ir) = trEVENT(1,ir) + enmin

               if( ns .le. 2 ) then

                  if( ityp .ne. 12 .and. ityp .ne. 13 ) then

                        trEVENT(11,ir)
     &                   = trEVENT(11,ir) + enmin

                  else

                        trEVENT(9,ir) = trEVENT(9,ir) + enmin

                  end if

               end if

               if( ns .eq. 1 ) then

                  if( ityp .eq. 1 ) then

                        k = 21

                  else if( ityp .eq. 3 ) then

                        k = 24

                  else if( ityp .eq. 5 ) then

                        k = 25

                  else if( ityp .eq. 12 .or. ityp .eq. 13 ) then

                        k = 26
                        rabs(4) = rabs(4) + enmin

                  else

                        k = 26

                  end if

                        trEVENT(k,ir) = trEVENT(k,ir) + enmin

                  if( np .gt. 0 ) then

                     call pcheck(m,np,ityp,ktyp,jtyp,ipn,ips)

                     if( ipn .gt. 0 ) then

                        do i = 1, ipn

                           k = 29 + np + ips(i)

                           trEVENT(k,ir) = trEVENT(k,ir) + enmin

                        end do

                     end if

                  end if

               end if

            end if

*-----------------------------------------------------------------------
*        termination by escape or leakage
*-----------------------------------------------------------------------

            if( ncol .eq. 12 ) then

               if( ns .le. 2 ) then

                        trEVENT(5,ir) = trEVENT(5,ir) + enmin

               end if

            end if

*-----------------------------------------------------------------------
*        neutron, proton and photon for nuclear data
*-----------------------------------------------------------------------

         if( mat .gt. 0 .and.
     &       ( jcoll .eq. 6 .or.
     &         jcoll .eq. 9 .or.
     &         ( ( jcoll .eq. 7 .or. jcoll .eq. 15 ) .and.   ! S.Abe 2017/01/17
     &           iegsemi .eq. 0 .and.   ! S.Abe 2015/08/31
     &           name(ibknam+no,ipomp+1) .eq. 1 .and.
     &                itelc(m) .eq. 0 ) ) ) then

*-----------------------------------------------------------------------

            icl = idgr(iblz1)
            erg = ( e(ibke+no,ipomp+1) + ec(ibkec+no,ipomp+1) ) / 2.0
            tgt = sqrt( ( xc(ibkxc+no,ipomp+1) - x(ibkx+no,ipomp+1) )**2
     &                + ( yc(ibkyc+no,ipomp+1) - y(ibky+no,ipomp+1) )**2
     &             + ( zc(ibkzc+no,ipomp+1) - z(ibkz+no,ipomp+1) )**2 )

               mk = mat
               rh = denm(mat)

            if( ityp .eq. 2 ) then

               call heatn(icl,erg,heatr,heatf,mk,rh,m,0,mtdum)

            else if( ityp .eq. 14 ) then

               call heatp(icl,erg,heatr,mk,rh,m,0,mtdum)

               heatf = 0.0

            else if( ityp .eq. 1 ) then

               call heath(icl,erg,heatr,mk,rh,m,0,mtdum)

               heatf = 0.0

            end if

            heatr = heatr * tgt * oldwt
            heatf = heatf * tgt * oldwt
            heat  = heatr + heatf

*-----------------------------------------------------------------------

                        trEVENT(1,ir) = trEVENT(1,ir) + heat

                  if( ns .le. 2 ) then

                     if( ityp .eq. 2 ) then

                        trEVENT(8,ir) = trEVENT(8,ir) + heat

                     else if( ityp .eq. 1 ) then

                        trEVENT(10,ir)
     &                      = trEVENT(10,ir) + heat

                     else

                        trEVENT(9,ir) = trEVENT(9,ir) + heat

                     end if

                  end if

                  if( ns .eq. 1 ) then

                    trEVENT(29,ir)
     &                   = trEVENT(29,ir) + heatf

                  end if

                  if( ityp .eq. 2 ) then

                        rabs(3) = rabs(3) + heat

                  else if( ityp .eq. 1 ) then

                        rabs(5) = rabs(5) + heat

                  else

                        rabs(4) = rabs(4) + heat

                  end if

*-----------------------------------------------------------------------

         end if

*-----------------------------------------------------------------------

  200 continue

*-----------------------------------------------------------------------
*     energy of produced particle and/or nucleus from nuclear reactions
*-----------------------------------------------------------------------
cc H.Iwase 2014/1/9  (add EGS particles)

      if( ( ncol .eq. 13 .or. ncol .eq. 14 ) .and.
     &    ( jcoll .lt. 6  .or. jcoll .eq. 10 .or.
     &      jcoll .eq. 12 .or. jcoll .eq. 16 .or. jcoll .eq. 17 .or.   ! S.Abe 2017/01/17
     &      jcoll .eq. 13 .or. jcoll .eq. 14 .or. jcoll .eq. 15 .or.   ! S.Abe 2017/01/17
     &      jcoll .eq. 11 .or.
     &      ( itelc(m) .eq. 1 .and.
     &        ( jcoll .eq. 7 .or. jcoll .eq. 8 ) ) ) ) then

                  totout = 0.0
                  totmas = 0.0

                  emathi = bindeg(mathz,mathn)

               if( ityp .ge. 15 .and. ityp .le. 19 ) then

                  ipnm = ktyp / 1000000
                  innm = ktyp - ktyp / 1000000 * 1000000 - ipnm

                  emathi = emathi + bindeg(ipnm,innm)

               end if

*-----------------------------------------------------------------------
*        from ordinary output
*-----------------------------------------------------------------------

         if( nclsts .gt. 0 ) then

*-----------------------------------------------------------------------
*        check of counter
*-----------------------------------------------------------------------

            do i = 1, 3

               if( itcnt(i,m) .ne. 0 ) then

                  if( jcount(i,1) .lt. itcnt(i*2+2,m) .or.
     &                jcount(i,1) .gt. itcnt(i*2+3,m) ) goto 300

               end if

            end do

*-----------------------------------------------------------------------

            do j = 1, nclsts

                  ipart = jclusts(3,j)
                  kpart = jclusts(7,j)
                  jpart = jclusts(5,j)

*-----------------------------------------------------------------------
*              nuclei
*-----------------------------------------------------------------------

               if( ipart .ge. 15 .and. ipart .le. 19 ) then

                     enrem  = qclusts(6,j) * qclusts(8,j)

                     totout = totout + enrem
                     totmas = totmas + bindeg(jclusts(1,j),jclusts(2,j))

*-----------------------------------------------------------------------
*                 enrem is remaining excitation energies,
*-----------------------------------------------------------------------

                        trEVENT(1,ir) = trEVENT(1,ir) + enrem

                  if( ns .le. 2 ) then

                    trEVENT(11,ir) = trEVENT(11,ir) + enrem

                  end if

                  if( ns .eq. 1 ) then

                    trEVENT(27,ir) = trEVENT(27,ir) + enrem

                     if( kcoll .eq. 1 ) then

                       trEVENT(28,ir)
     &                      = trEVENT(28,ir) + enrem

                     end if

                  end if

               end if

*-----------------------------------------------------------------------
*              out going particles and nuclei
*-----------------------------------------------------------------------

                  if( ipart .le. 2 .or. ipart .ge. 12 ) then
                     enoms = 0.0
                  else
                     enoms = ( qclusts(5,j) - jclusts(6,j) * rmngv )
                  end if

                     enout = qclusts(7,j)
                     enmas = enoms * 1000.0
                     enpar = ( enout + enmas ) * qclusts(8,j)

                     totout = totout + enpar

*-----------------------------------------------------------------------
*              dead particle and nucleus
*-----------------------------------------------------------------------

               if( jclusts(4,j) .lt. 0 ) then

                     if( jclusts(4,j) .eq. -1 ) then

                        trEVENT(1,ir) = trEVENT(1,ir) + enpar

                     end if

                  if( ns .le. 2 ) then

                     if( ipart .eq. 2 .and.
     &                   jclusts(4,j) .eq. -2 ) then

                        k = 2

                     else if( ipart .ge. 12 .and. ipart .le. 14 .and.
     &                        jclusts(4,j) .eq. -2 ) then

                        k = 3

                     else if( ipart .eq. 1 .and.
     &                        jclusts(4,j) .eq. -2 ) then

                        k = 4

                     else if( ipart .ge. 15 .and. ipart .le. 19 ) then

                        k = 6

                     else if( ipart .eq. 12 .or. ipart .eq. 13 ) then

                        k = 9
                        rabs(4) = rabs(4) + enpar

                     else if( ipart .eq. 14 .and. ns .eq. 2 ) then

                        k = 11

                     else if( ipart .eq. 14 .and. ns .eq. 1 ) then

                        k = 0

                     else

                        k = 11

                     end if

                     if( k .gt. 0 ) then

                        trEVENT(k,ir) = trEVENT(k,ir) + enpar

                     end if

                  end if

                  if( ns .eq. 1 .and. jclusts(4,j) .eq. -1 ) then

                     if( ipart .eq. 1 ) then

                        k = 21

                     else if( ipart .eq. 2 ) then

                        k = 22

                     else if( ipart .eq. 14 ) then

                        k = 23

                     else if( ipart .ge. 15 .and.
     &                        ipart .le. 19) then

                        k = ipart - 3

                     else if( ipart .eq. 3 ) then

                        k = 24

                     else if( ipart .eq. 5 ) then

                        k = 25

                     else if( ipart .ge. 12 .and. ipart .le. 13 ) then

                        k = 26

                     else

                        k = 26

                     end if

                     if( k .gt. 0 ) then

                        trEVENT(k,ir) = trEVENT(k,ir) + enpar

                     end if

                     if( np .gt. 0 ) then

                        call pcheck(m,np,ipart,kpart,jpart,ipn,ips)

                        if( ipn .gt. 0 ) then

                           do i = 1, ipn

                              k = 29 + np + ips(i)

                              trEVENT(k,ir) = trEVENT(k,ir) + enpar

                           end do

                        end if

                     end if

                     if( kcoll .eq. 1 ) then

                           trEVENT(28,ir) = trEVENT(28,ir) + enpar

                     end if

                  end if

               end if

*-----------------------------------------------------------------------

            end do

  300       continue

         end if

*-----------------------------------------------------------------------
*        total absorption energy
*-----------------------------------------------------------------------

         if( jcoll .lt. 6  .or. jcoll .eq. 10 .or.
     &       jcoll .eq. 12 .or. jcoll .eq. 16 .or. jcoll .eq. 17 .or.   ! S.Abe 2017/01/17
     &       jcoll .eq. 13 .or. jcoll .eq. 14 .or. jcoll .eq. 15 .or.   ! S.Abe 2017/01/17
     &       jcoll .eq. 11 ) then

            rabs(1) = rabs(1) + enmin - totout
            rabs(2) = rabs(2) + ( emathi - totmas ) * oldwt

*-----------------------------------------------------------------------
*           debug for energy conservation
*-----------------------------------------------------------------------




         end if

*-----------------------------------------------------------------------

      end if

*-----------------------------------------------------------------------

  100 continue

*-----------------------------------------------------------------------

      return
      end

************************************************************************
*                                                                      *
      subroutine thetrzEVENT(ncol,m,nd,ns0,nl,lt,np,nr,nz,rm,zm,
     &                  ne,eb,tr,rabs,trEVENT,tr0)
*                                                                      *
*       nuclear heating tally in r-z scoring mesh                      *
*       output = deposit-*                                             *
*       last modified by T.Furuta on 2012/05/08                        *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*      ncol  ..... reaction type                                       *
*              0 : end of batch                                        *
*              4 : source particle                                     *
*                                                                      *
*              9 : termination by time cut-off                         *
*             10 : geometry boundary crossing                          *
*             11 : termination by energy cut-off                       *
*             12 : termination by escape or leakage                    *
*             13 : nuclear reaction (n,x)                              *
*             14 : nuclear reaction (n,n'x)                            *
*             15 : sequential transport only for tally                 *
*                                                                      *
*     definition of tr( nd, ir, iz, id )                               *
*                                                                      *
*        ns = 1 : all, 2 : simple, 3 : heat                            *
*                                                                      *
*      nd ..... number of event type                                   *
*        ns = 3                                                        *
*          1 : total                                                   *
*                                                                      *
*        ns = 2                                                        *
*          2 : cutoff neutron                                          *
*          3 : cutoff photon                                           *
*          4 : cutoff proton                                           *
*          5 : leakage                                                 *
*                                                                      *
*          6 : recoil                                                  *
*          7 : ionization                                              *
*          8 : low neutron                                             *
*          9 : electron                                                *
*         10 : low proton                                              *
*         11 : other                                                   *
*                                                                      *
*        ns = 1                                                        *
*         12 : deuteron recoil                                         *
*         13 : triton   recoil                                         *
*         14 : 3He      recoil                                         *
*         15 : alpha    recoil                                         *
*         16 : residual recoil                                         *
*                                                                      *
*         17 : proton ionization                                       *
*         18 : pion+  ionization                                       *
*         19 : pion-  ionization                                       *
*         20 : other  ionization                                       *
*                                                                      *
*         21 : stopped proton                                          *
*         22 : stopped neutron                                         *
*         23 : stopped photon                                          *
*         24 : stopped pion+                                           *
*         25 : stopped pion-                                           *
*         26 : stopped other                                           *
*                                                                      *
*         27 : remaining excitation energy                             *
*         28 : fission high recoil                                     *
*         29 : fission low neutron                                     *
*                                                                      *
*         29 + np      : ionization of other paticle                   *
*         29 + np + np : stopped other paticle                         *
*                                                                      *
*       id ..... data type                                             *
*          1 : heating                                                 *
*          2 : relative error                                          *
*                                                                      *
*       unit = MeV                                                     *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
************************************************************************
      use MMBANKMOD   !FURUTA
      use EVENTTALMOD !FURUTA
      use NGSDATAMOD, only : bindeg
      use partmod, only: itmxpt, itpan, itpat, jtpat ! frtati 2021/10/05
*-----------------------------------------------------------------------

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param00.inc'
      include 'param.inc'

*-----------------------------------------------------------------------

      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
!$OMP THREADPRIVATE(/icomon/)
      common /wtsave/ oldwt
!$OMP THREADPRIVATE(/wtsave/)
      common /tlgeom/ iblz1,iblz2
!$OMP THREADPRIVATE(/tlgeom/)
      common /mathzn/ mathz, mathn, jcoll, kcoll
!$OMP THREADPRIVATE(/mathzn/)
      common /regdc/  idrg(kvlmax), idgr(kvmmax)
      common /kmat1d/ idmn(0:kvlmax), idnm(kvmmax)
      common /eparm/  esmax, esmin, emin(20)
      common /celdn/  denr(kvlmax), denm(kvlmax), denc(kvlmax)

*-----------------------------------------------------------------------

      common /clustt/ nclsts, iclusts(nnn)
!$OMP THREADPRIVATE(/clustt/)
      common /clustw/ jclusts(0:8,nnn),  qclusts(0:12,nnn)
!$OMP THREADPRIVATE(/clustw/)
      common /cntcls/ jcount(3,nnn)
!$OMP THREADPRIVATE(/cntcls/)

*-----------------------------------------------------------------------

      common /tall08/ rtrx0(itlmax), rtry0(itlmax)
      common /tall22/ itmcn(itlmax), itmtn(itlmax), itmtt(itlmax)
      common /tall32/ itelc(itlmax)
      common /tall37/ itcnt(9,itlmax)
      common /tall46/ itmtr(itlmax,4), rtmtr(itlmax,13)

*-----------------------------------------------------------------------

      dimension   lt(nl)
      dimension   rm(nr+1)
      dimension   zm(nz+1)
      dimension   tr(nd,nr,nz,0:ne,2)
      dimension   rabs(5)
      dimension   eb(ne+1)
      dimension   trEVENT(nd,nr,nz)

      dimension ips(itmxpt) ! frtati 2021/10/05 6 -> itmxpt

*-----------------------------------------------------------------------

      data rmngv / 0.93895 /
      data rmnmv / 938.95 /

*-----------------------------------------------------------------------

      dimension ud(3)

      data dmax  /1.0d+19/
      data dmax0 /1.0d+18/

*-----------------------------------------------------------------------

      save sweight
!$OMP THREADPRIVATE(sweight)
      common /jcomon/ nabov,nobch,nocas,nomax
!$OMP THREADPRIVATE(/jcomon/)
      common /tall36/ itdpo(itlmax)

*-----------------------------------------------------------------------

      common /paraj/ mstz(300), parz(300)

*-----------------------------------------------------------------------

      real*8 edep
      common/EPCONTedep/ edep
!$OMP THREADPRIVATE(/EPCONTedep/)

      integer iegsemi, iegsout
      common /egsemi/ iegsemi, iegsout

C for nonshared_tally option
      integer italsh
      common /talsh/ italsh

*------------------------------------------------------------------------

      common /spred/ nspred, nwsprd, nedisp, itstep, ndedx

*-----------------------------------------------------------------------
! sumover
      real(8) tr0(nd,nr,nz,0:ne)

      tr0(:,:,:,:) = 0.0d0

            small = parz(28) * 1.d-4 ! T.Sato 2023/11/01 avoid infinite loop

*-----------------------------------------------------------------------
*     ns from ns0
*-----------------------------------------------------------------------

            ns = ns0
            if( ns .ge. 4 ) ns = ns - 3

*-----------------------------------------------------------------------
*     deposit : ncol = 0, 4 and ns = 4
*-----------------------------------------------------------------------

      if( ( ncol .eq. 0 .or. ncol .eq. 4 ) .and. ns0 .ge. 4 ) then

         if( .not. FIRSTsrc ) then

            do id = 1, nd
            do ir = 1, nr
            do iz = 1, nz

               if( itdpo(m) .eq. 0
     &             .and. trEVENT(1,ir,iz) .gt. 0.0 ) then

                  heats = trEVENT(1,ir,iz)
                  ratio = trEVENT(id,ir,iz) / trEVENT(1,ir,iz)

               else

                  heats = trEVENT(id,ir,iz)
                  ratio = 1.0d0

               end if

               if( heats .gt. 0.0d0 .and.
     &             heats .ge. eb(1) .and. heats .lt. eb(ne+1) ) then

                  do ie = 1, ne

                     if( heats .ge. eb(ie) .and.
     &                   heats .lt. eb(ie+1) ) then

! sumover
                       tr0(id,ir,iz,ie) = tr0(id,ir,iz,ie)
     &                                     + sweight*ratio
                       if (italsh .eq. 0 ) then
                         tr(id,ir,iz,ie,1) = tr(id,ir,iz,ie,1)
     &                                     + sweight*ratio
                         tr(id,ir,iz,ie,2) = tr(id,ir,iz,ie,2)
     &                                     + (sweight*ratio)**2
                       else
!$OMP CRITICAL (tr_thetrzEVENT)
                        tr(id,ir,iz,ie,1) = tr(id,ir,iz,ie,1)
     &                                    + sweight*ratio
                        tr(id,ir,iz,ie,2) = tr(id,ir,iz,ie,2)
     &                                    + (sweight*ratio)**2
!$OMP END CRITICAL (tr_thetrzEVENT)
                      end if

                     end if

                  end do

               end if

            end do
            end do
            end do

! sumover
           call thetrz_sumover(m, 1,
     &                            nd, nr, nz, ne, tr0)
           tr0(:,:,:,:) = 0.0d0

            do id = 1, nd
            do ir = 1, nr
            do iz = 1, nz

               trEVENT(id,ir,iz) = 0.0d0

            end do
            end do
            end do

         end if

*-----------------------------------------------------------------------
*        save initial weight and zero set
*-----------------------------------------------------------------------

         if( ncol .eq. 4 ) then

            sweight = oldwt

         end if

      end if

*-----------------------------------------------------------------------
*        check of ncol
*-----------------------------------------------------------------------

         if( ncol .lt. 9 ) return

*-----------------------------------------------------------------------
*        nothing for electron and positron with nuclear data
*                    and secondary photon
*-----------------------------------------------------------------------

         if( itelc(m) .eq. 0 .and. iegsemi .eq. 0 .and.   ! S.Abe 2015/08/31
     &     ( jcoll .eq. 8 .or.
     &     ( jcoll .eq. 7 .and. name(ibknam+no,ipomp+1) .gt. 1 ) ) )
     &              return

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
*        transform positions
*-----------------------------------------------------------------------

            call trnsxx(x(ibkx+no,ipomp+1),y(ibky+no,ipomp+1),
     &                  z(ibkz+no,ipomp+1),
     &                  xxa,yya,zza,itmtr(m,4))

            call trnsxx(xc(ibkxc+no,ipomp+1),yc(ibkyc+no,ipomp+1),
     &                  zc(ibkzc+no,ipomp+1),
     &                  xxc,yyc,zzc,itmtr(m,4))

*-----------------------------------------------------------------------
* add dummy step for EGS5 electrons at ncol=11
*     (stopped but still have dE)
cc H.Iwase 2014/8/21 add the change of H.Iwase 2014/1/9
*-----------------------------------------------------------------------

         if( iegsemi .ne. 0 .and.
     &     ( ityp .eq. 12 .or. ityp .eq. 13 ) .and. ncol .eq. 11 ) then

               dum = ( xxc - xxa )**2
     &             + ( yyc - yya )**2
     &             + ( zzc - zza )**2

            if( dum .le. small**2 ) then

               dumstep = small * 10.d0

               xdum = x(ibkx+no,ipomp+1) + u(ibku+no,ipomp+1)*dumstep
               ydum = y(ibky+no,ipomp+1) + v(ibkv+no,ipomp+1)*dumstep
               zdum = z(ibkz+no,ipomp+1) + w(ibkw+no,ipomp+1)*dumstep

               call trnsxx(xdum,ydum,zdum,xxc,yyc,zzc,itmtr(m,4))

            end if

         end if
*-----------------------------------------------------------------------
*        check z mesh and r mesh region
*-----------------------------------------------------------------------

               x0 = rtrx0(m)
               y0 = rtry0(m)

            if( zza-small .lt. zm(1) .and.
     &          zzc-small .lt. zm(1) ) return

            if( zza+small .ge. zm(nz+1) .and.
     &          zzc+small .ge. zm(nz+1) ) return

               dis0 = sqrt( ( xxa - x0 )**2
     &                    + ( yya - y0 )**2 )

               dis1 = sqrt( ( xxc - x0 )**2
     &                    + ( yyc - y0 )**2 )

            if( dis0 .lt. rm(1) .and.
     &          dis1 .lt. rm(1) ) return

            if( dis0 .ge. rm(nr+1) .and.
     &          dis1 .ge. rm(nr+1) ) then

               aa =   yyc - yya
               bb = - xxc + xxa
               cc = - aa * xxa - bb * yya

               if( aa**2 + bb**2 .ne. 0.0d0 ) then

                  dd = abs( aa * x0 + bb * y0 + cc )
     &               / sqrt( aa**2 + bb**2 )

                  if( dd .ge. rm(nr+1) ) return

               end if

            end if

*-----------------------------------------------------------------------
*        enmin : final particle energy ( with mass for not nucleon )
*-----------------------------------------------------------------------

               if( ityp .le. 2 .or. ityp .ge. 12 ) then
                  enmin = 0.0
               else
                  enmin = rtyp - mtyp * rmnmv
               end if

              enmin = ( ec(ibkec+no,ipomp+1) + enmin ) * oldwt


*-----------------------------------------------------------------------
*        check of counter
*-----------------------------------------------------------------------

            do i = 1, 3

               if( itcnt(i,m) .ne. 0 ) then

                  if( ncnt(ibknct+i,no,ipomp+1) .lt. itcnt(i*2+2,m) .or.
     &                ncnt(ibknct+i,no,ipomp+1) .gt. itcnt(i*2+3,m) )
     &                goto 200

               end if

            end do

*-----------------------------------------------------------------------
*     ncol = 9, 11, 12 : time and enrgy cut-off, escape
*     reaction point is out of rainge -> 500
*-----------------------------------------------------------------------

      if( ncol .eq. 9 .or. ncol .eq. 11 .or. ncol .eq. 12 ) then

               if( zzc .lt. zm(1) .or.
     &             zzc .ge. zm(nz+1) ) goto 500

               if( dis1 .lt. rm(1) .or.
     &             dis1 .ge. rm(nr+1) ) goto 500

*-----------------------------------------------------------------------
*        get end point : iz and ir
*-----------------------------------------------------------------------

               do i = 1, nz

                  if( zzc .ge. zm(i) .and.
     &                zzc .lt. zm(i+1) ) goto 11

               end do

   11             iz = i

               do i = 1, nr

                  if( dis1 .ge. rm(i) .and.
     &                dis1 .lt. rm(i+1) ) goto 21

               end do

   21             ir = i

               if( iz .ge. nz + 1 .or. ir .ge. nr + 1 ) goto 500

*-----------------------------------------------------------------------
*           stopped energy cut-off or time cut-off particles
*-----------------------------------------------------------------------

            if( ncol .eq. 9 .or. ncol .eq. 11 ) then

              trEVENT(1,ir,iz) = trEVENT(1,ir,iz) + enmin

               if( ns .le. 2 ) then

                  if( ityp .ne. 12 .and. ityp .ne. 13 ) then

                    trEVENT(11,ir,iz)
     &                   = trEVENT(11,ir,iz) + enmin

                  else

                    trEVENT(9,ir,iz)
     &                   = trEVENT(9,ir,iz) + enmin

                  end if

               end if

               if( ns .eq. 1 ) then

                  if( ityp .eq. 1 ) then

                        k = 21

                  else if( ityp .eq. 3 ) then

                        k = 24

                  else if( ityp .eq. 5 ) then

                        k = 25

                  else if( ityp .eq. 12 .or. ityp .eq. 13 ) then

                        k = 26
                        rabs(4) = rabs(4) + enmin

                  else

                        k = 26

                  end if

                  trEVENT(k,ir,iz) = trEVENT(k,ir,iz) + enmin

                  if( np .gt. 0 ) then

                     call pcheck(m,np,ityp,ktyp,jtyp,ipn,ips)

                     if( ipn .gt. 0 ) then

                        do i = 1, ipn

                           k = 29 + np + ips(i)

                           trEVENT(k,ir,iz)
     &                          = trEVENT(k,ir,iz) + enmin

                        end do

                     end if

                  end if

               end if

            end if

*-----------------------------------------------------------------------
*           termination by escape or leakage
*-----------------------------------------------------------------------

            if( ncol .eq. 12 ) then

               if( ns .le. 2 ) then

                 trEVENT(5,ir,iz) = trEVENT(5,ir,iz) + enmin

               end if

            end if

      end if

*-----------------------------------------------------------------------
*     ionization or track length for low energy neutron and photon
*-----------------------------------------------------------------------

  500 if( mat .le. 0 ) goto 200

*-----------------------------------------------------------------------
*     neutron, proton and photon for nuclear data
*-----------------------------------------------------------------------

         if( jcoll .eq. 6 .or.
     &       jcoll .eq. 9 .or.   ! S.Abe 2017/01/17
     &       ( ( jcoll .eq. 7 .or. jcoll .eq. 15 ) .and.   ! S.Abe 2017/01/17
     &         iegsemi .eq. 0 .and.   ! S.Abe 2015/08/31
     &         name(ibknam+no,ipomp+1) .eq. 1 .and. itelc(m) .eq. 0 ) )
     &      then

*-----------------------------------------------------------------------

               icl = idgr(iblz1)
               erg = ( e(ibke+no,ipomp+1) + ec(ibkec+no,ipomp+1) ) / 2.0

               mk = mat
               rh = denm(mat)

            if( ityp .eq. 2 ) then

               call heatn(icl,erg,heatr,heatf,mk,rh,m,0,mtdum)

            else if( ityp .eq. 14 ) then

               call heatp(icl,erg,heatr,mk,rh,m,0,mtdum)

               heatf = 0.0

            end if

            heatr0 = heatr * oldwt
            heatf0 = heatf * oldwt
            heat0  = heatr0 + heatf0

            inpd = 1

         else if( jcoll .eq. 9 ) then

            inpd = 2

         else

            inpd = 0

         end if

*-----------------------------------------------------------------------
*        only for charged particles
*-----------------------------------------------------------------------
cc H.Iwase 2014/1/9  (EGS5)

         if( inpd .eq. 0 ) then

            if( iegsemi .eq. 0 ) then

               if( jtyp .eq. 0 ) goto 200
               if( abs( ec(ibkec+no,ipomp+1) - e(ibke+no,ipomp+1) )
     &             .eq. 0.0d0 )
     &             goto 200

             else

               if( jtyp .eq. 0 .and. ityp .ne. 14 ) goto 200
               if( ( abs( ec(ibkec+no,ipomp+1) - e(ibke+no,ipomp+1) )
     $                           .eq. 0.0d0 )
     $              .and.
     $          ( ityp .ne. 12 .and. ityp .ne. 13 .and. ityp .ne. 14))
     $            goto 200

             end if

         end if

*-----------------------------------------------------------------------
*        distance and unit vector ud(i)
*-----------------------------------------------------------------------

            dis = ( xxc - xxa )**2
     &          + ( yyc - yya )**2
     &          + ( zzc - zza )**2

            if( dis .le. small**2 ) goto 200

            dis = sqrt(dis)

            ud(1) = ( xxc - xxa ) / dis
            ud(2) = ( yyc - yya ) / dis
            ud(3) = ( zzc - zza ) / dis

*-----------------------------------------------------------------------
*     initial position, energy and initial range rng
*-----------------------------------------------------------------------

            tot = 0.0d0

            se  = e(ibke+no,ipomp+1)
            xpp = xxa
            ypp = yya
            zpp = zza

            if( inpd .ne. 1 )
     &      call rainge(se,rng,mat,ityp,ktyp,jtyp,rtyp)

*-----------------------------------------------------------------------
*     initial z-position
*-----------------------------------------------------------------------

               izm = 0
               izc = 0

            if( ud(3) .gt. 0.0d0 ) izk =  1
            if( ud(3) .eq. 0.0d0 ) izk =  0
            if( ud(3) .lt. 0.0d0 ) izk = -1

            if( izk .ge. 0 ) then

               do i = 1, nz + 1
                  if( zm(i) .gt. zpp + small ) then
                         izm = i - 1
                         izc = i - 1
                         goto 38
                  end if
               end do

            else

               izm = nz + 2
               izc = nz + 1

               do i = 1, nz + 1
                  if( zm(i) .ge. zpp - small ) then
                         izm = i
                         izc = i - 1
                         goto 38
                  end if
               end do

            end if

   38       continue

*-----------------------------------------------------------------------
*     initial r-position
*-----------------------------------------------------------------------

               uvec = ( xxa - x0 ) * ud(1) + ( yya - y0 ) * ud(2)

            if( uvec .ge. 0.0d0 ) then

               do i = 1, nr + 1
                  if( rm(i) .gt. dis0 + small ) then
                        irc = i - 1
                        goto 47
                  end if
               end do

            else

               do i = 1, nr + 1
                  if( rm(i) .ge. dis0 - small ) then
                        irc = i - 1
                        goto 47
                  end if
               end do

            end if

                  irc = nr + 1

   47          continue

*-----------------------------------------------------------------------
*     loop for finding mesh and booking upto total distance
*-----------------------------------------------------------------------

   50 continue

*-----------------------------------------------------------------------
*      calculation is finished
*-----------------------------------------------------------------------

         if( tot .ge. dis ) goto 200

*-----------------------------------------------------------------------
*        dr : distance to the nearest r mesh
*-----------------------------------------------------------------------

               dr = dmax

               a  = ud(1)**2 + ud(2)**2

         if( a .gt. 0.0d0 ) then

               b = 2.0d0 * ( ud(1) * ( xpp - x0 )
     &                     + ud(2) * ( ypp - y0 ) )

               ir1 = irc
               ir2 = irc + 1

              if( b .ge. 0.0 ) ir1 = irc + 1

            do 51 i = max( 1, ir1 ), min( nr + 1, ir2 )

               if( rm(i) .eq. 0.0d0 ) goto 51

               c = ( xpp - x0 )**2
     &           + ( ypp - y0 )**2 - rm(i)**2

               dchk = b**2 - 4.0d0 * a * c

               if( dchk .lt. 0.0d0 ) goto 51

               t1 = ( -b + sqrt(dchk) ) / ( a * 2.0 )
               t2 = ( -b - sqrt(dchk) ) / ( a * 2.0 )

               if( t1 .lt. small ) t1 = dmax
               if( t2 .lt. small ) t2 = dmax

               ttr = min(t1,t2)

               if( ttr .lt. dr ) then

                  dr  = ttr
                  idr = i

               end if

   51       continue

         end if

*-----------------------------------------------------------------------
*        dz : distance to the nearest z mesh
*-----------------------------------------------------------------------

               dz = dmax

         if( izk .ne. 0 ) then
            if( izm + izk .ge. 1 .and.
     &          izm + izk .le. nz + 1 ) then

               dz = abs( ( zm(izm+izk) - zpp ) / ud(3) )

            end if
         end if

*-----------------------------------------------------------------------
*        which boundary is the nearlist
*-----------------------------------------------------------------------

            if( dr .gt. dz ) then

               dd = dz
               jz = 1

            else

               dd = dr
               jz = 0

            end if

               if( dd .gt. dmax0 ) goto 200

*-----------------------------------------------------------------------
*        recalculate rng for ncol = 11 and nedisp
*-----------------------------------------------------------------------

            if( mat .gt. 0 .and. jtyp .ne. 0 .and. nedisp .ne. 0 .and.
     &           (  ityp .ne. 12 .and. ityp .ne. 13 ) ) then


               call rainge(se,rng,mat,ityp,ktyp,jtyp,rtyp)

               if( ncol .eq. 11 )  dis = tot + rng

            end if

*-----------------------------------------------------------------------
*        propagate position upto the boundary or the final point
*-----------------------------------------------------------------------

            if( tot + dd .ge. dis - small ) then

               dd  = dis - tot
               tot = dis
               xpp = xxc
               ypp = yyc
               zpp = zzc

            else

               tot = tot + dd
               xpp = xpp + dd * ud(1)
               ypp = ypp + dd * ud(2)
               zpp = zpp + dd * ud(3)

            end if

*-----------------------------------------------------------------------
*        final energy ee and range rng
*-----------------------------------------------------------------------

            if( inpd .ne. 1 ) then

               delt = dd

               call ecol(ee,delt,se,rng,
     &                   mat,ityp,ktyp,jtyp,rtyp)

               rng = max( 0.0d0, rng - delt )
               ee  = max( 0.0d0, ee )
               ee  = min( se, ee )

               enion = abs( se - ee ) * oldwt

               if( iegsemi .ne. 0 .and. ityp .eq. 14 ) then
                  enion = enion + edep * oldwt
               end if

            else

               ee = se

               heatr = heatr0 * dd
               heatf = heatf0 * dd
               heat  = heat0  * dd

                  if( ityp .eq. 2 ) then

                     rabs(3) = rabs(3) + heat

                  else

                     rabs(4) = rabs(4) + heat

                  end if

            end if

*-----------------------------------------------------------------------
*        booking the ionization energy
*-----------------------------------------------------------------------

            iz = izc
            ir = irc

         if( ir .ge. 1 .and. ir .lt. nr + 1 .and.
     &       iz .ge. 1 .and. iz .lt. nz + 1 ) then

            if( inpd .ne. 1 ) then

              trEVENT(1,ir,iz) = trEVENT(1,ir,iz) + enion

               if( ns .le. 2 ) then

                 trEVENT(7,ir,iz) = trEVENT(7,ir,iz) + enion

               end if

               if( ns .eq. 1 ) then

                  if( ityp .eq. 1 ) then

                        k = 17

                  else if( ityp .eq. 3 ) then

                        k = 18

                  else if( ityp .eq. 5 ) then

                        k = 19

                  else if( ityp .eq. 12 .or. ityp .eq. 13 ) then

                        k = 20
                        rabs(4) = rabs(4) + enion

                  else

                        k = 20

                  end if

                  trEVENT(k,ir,iz) = trEVENT(k,ir,iz) + enion

                  if( np .gt. 0 ) then

                     call pcheck(m,np,ityp,ktyp,jtyp,ipn,ips)

                     if( ipn .gt. 0 ) then

                        do i = 1, ipn

                           k = 29 + ips(i)

                           trEVENT(k,ir,iz)
     &                          = trEVENT(k,ir,iz) + enion

                        end do

                     end if

                  end if

               end if

            end if

*-----------------------------------------------------------------------
*           neutron and photon for nuclear data
*-----------------------------------------------------------------------

            if( inpd .eq. 1 ) then

              trEVENT(1,ir,iz) = trEVENT(1,ir,iz) + heat

                  if( ns .le. 2 ) then

                     if( ityp .eq. 2 ) then

                       trEVENT(8,ir,iz)
     &                      = trEVENT(8,ir,iz) + heat

                     else

                       trEVENT(9,ir,iz)
     &                      = trEVENT(9,ir,iz) + heat

                     end if

                  end if

                  if( ns .eq. 1 ) then

                    trEVENT(29,ir,iz)
     &                   = trEVENT(29,ir,iz) + heatf

                  end if

*-----------------------------------------------------------------------
*           proton for nuclear data
*-----------------------------------------------------------------------

            else if( inpd .eq. 2 ) then

                     erg = ( se + ee ) / 2.0
                     icl = idgr(iblz1)

                     mk = mat
                     rh = denm(mat)

                     call heath(icl,erg,heatr,mk,rh,m,0,mtdum)

                     heatf = 0.0

                     heatr0 = heatr * oldwt
                     heatf0 = heatf * oldwt
                     heat0  = heatr0 + heatf0

                     heatr = heatr0 * dd
                     heatf = heatf0 * dd
                     heat  = heat0  * dd

                     rabs(5) = rabs(5) + heat

                     trEVENT(1,ir,iz)
     &                    = trEVENT(1,ir,iz) + heat

                  if( ns .le. 2 ) then

                    trEVENT(10,ir,iz)
     &                   = trEVENT(10,ir,iz) + heat

                  end if

            end if

         end if

*-----------------------------------------------------------------------
*     next position
*     for z-crossing jz = 1, r-crossing jz = 0
*-----------------------------------------------------------------------

            if( jz .eq. 1 ) then

                  izm = izm + izk
                  izc = izc + izk

            else if( jz .eq. 0 ) then

               if( idr .gt. irc ) then

                  irc = idr

               else

                  irc = idr - 1

               end if

            end if

         if( inpd .ne. 1 ) then

            se  = ee

         end if

*-----------------------------------------------------------------------

      goto 50

*-----------------------------------------------------------------------
*     producced particles
*-----------------------------------------------------------------------

  200 continue

*-----------------------------------------------------------------------
*     energy of produced particle and/or nucleus from nuclear reactions
*-----------------------------------------------------------------------
cc H.Iwase 2014/1/9 (EGS5)

      if( ( ncol .eq. 13 .or. ncol .eq. 14 ) .and.
     &    ( jcoll .lt. 6  .or. jcoll .eq. 10 .or.
     &      jcoll .eq. 12 .or. jcoll .eq. 16 .or. jcoll .eq. 17 .or.   ! S.Abe 2017/01/17
     &      jcoll .eq. 13 .or. jcoll .eq. 14 .or. jcoll .eq. 15 .or.   ! S.Abe 2014/01/17
     &      jcoll .eq. 11 .or.
     &      ( itelc(m) .eq. 1 .and.
     &        ( jcoll .eq. 7 .or. jcoll .eq. 8 ) ) ) ) then

                  totout = 0.0
                  totmas = 0.0

                  emathi = bindeg(mathz,mathn)

               if( ityp .ge. 15 .and. ityp .le. 19 ) then

                  ipnm = ktyp / 1000000
                  innm = ktyp - ktyp / 1000000 * 1000000 - ipnm

                  emathi = emathi + bindeg(ipnm,innm)

               end if

*-----------------------------------------------------------------------
*        from ordinary output
*-----------------------------------------------------------------------

         if( nclsts .gt. 0 ) then

*-----------------------------------------------------------------------
*        check of counter
*-----------------------------------------------------------------------

            do i = 1, 3

               if( itcnt(i,m) .ne. 0 ) then

                  if( jcount(i,1) .lt. itcnt(i*2+2,m) .or.
     &                jcount(i,1) .gt. itcnt(i*2+3,m) ) goto 300

               end if

            end do

*-----------------------------------------------------------------------

         do 600 j = 1, nclsts

               ipart = jclusts(3,j)
               kpart = jclusts(7,j)
               jpart = jclusts(5,j)

*-----------------------------------------------------------------------
*           collision position
*-----------------------------------------------------------------------

               call trnsxx(qclusts(10,j),
     &                     qclusts(11,j),
     &                     qclusts(12,j),
     &                     xcc,ycc,zcc,itmtr(m,4))

               x0 = rtrx0(m)
               y0 = rtry0(m)

               if( zcc .lt. zm(1) .or.
     &             zcc .ge. zm(nz+1) ) goto 600

               dis1 = sqrt( (  xcc - x0 )**2
     &                    + (  ycc - y0 )**2 )

               if( dis1 .lt. rm(1) .or.
     &             dis1 .ge. rm(nr+1) ) goto 600

               do i = 1, nz

                  if( zcc .ge. zm(i) .and.
     &                zcc .lt. zm(i+1) ) goto 31

               end do

   31             iz = i

               do i = 1, nr

                  if( dis1 .ge. rm(i) .and.
     &                dis1 .lt. rm(i+1) ) goto 41

               end do

   41             ir = i

               if( iz .ge. nz + 1 .or. ir .ge. nr + 1 ) goto 600

*-----------------------------------------------------------------------
*              nuclei
*-----------------------------------------------------------------------

               if( ipart .ge. 15 .and. ipart .le. 19 ) then

                     enrem  = qclusts(6,j) * qclusts(8,j)

                     totout = totout + enrem
                     totmas = totmas + bindeg(jclusts(1,j),jclusts(2,j))

*-----------------------------------------------------------------------
*                 enrem is remaining excitation energies,
*-----------------------------------------------------------------------

                     trEVENT(1,ir,iz)
     &                    = trEVENT(1,ir,iz) + enrem

                  if( ns .le. 2 ) then

                    trEVENT(11,ir,iz)
     &                   = trEVENT(11,ir,iz) + enrem

                  end if

                  if( ns .eq. 1 ) then

                    trEVENT(27,ir,iz)
     &                   = trEVENT(27,ir,iz) + enrem

                     if( kcoll .eq. 1 ) then

                       trEVENT(28,ir,iz)
     &                      = trEVENT(28,ir,iz) + enrem

                     end if

                  end if

               end if

*-----------------------------------------------------------------------
*              out going particles and nuclei
*-----------------------------------------------------------------------

                  if( ipart .le. 2 .or. ipart .ge. 12 ) then
                     enoms = 0.0
                  else
                     enoms = ( qclusts(5,j) - jclusts(6,j) * rmngv )
                  end if

                     enout = qclusts(7,j)
                     enmas = enoms * 1000.0
                     enpar = ( enout + enmas ) * qclusts(8,j)

                     totout = totout + enpar

*-----------------------------------------------------------------------
*              dead particle and nucleus
*-----------------------------------------------------------------------

               if( jclusts(4,j) .lt. 0 ) then

                     if( jclusts(4,j) .eq. -1 ) then

                        trEVENT(1,ir,iz)
     &                      = trEVENT(1,ir,iz) + enpar

                     end if

                  if( ns .le. 2 ) then

                     if( ipart .eq. 2 .and.
     &                   jclusts(4,j) .eq. -2 ) then

                        k = 2

                     else if( ipart .ge. 12 .and. ipart .le. 14 .and.
     &                        jclusts(4,j) .eq. -2 ) then

                        k = 3

                     else if( ipart .eq. 1 .and.
     &                        jclusts(4,j) .eq. -2 ) then

                        k = 4

                     else if( ipart .ge. 15 .and. ipart .le. 19 ) then

                        k = 6

                     else if( ipart .eq. 12 .or. ipart .eq. 13 ) then

                        k = 9
                        rabs(4) = rabs(4) + enpar

                     else if( ipart .eq. 14 .and. ns .eq. 2 ) then

                        k = 11

                     else if( ipart .eq. 14 .and. ns .eq. 1 ) then

                        k = 0

                     else

                        k = 11

                     end if

                     if( k .gt. 0 ) then

                        trEVENT(k,ir,iz)
     &                      = trEVENT(k,ir,iz) + enpar

                     end if

                  end if

                  if( ns .eq. 1 .and. jclusts(4,j) .eq. -1 ) then

                     if( ipart .eq. 1 ) then

                        k = 21

                     else if( ipart .eq. 2 ) then

                        k = 22

                     else if( ipart .eq. 14 ) then

                        k = 23

                     else if( ipart .ge. 15 .and.
     &                        ipart .le. 19) then

                        k = ipart - 3

                     else if( ipart .eq. 3 ) then

                        k = 24

                     else if( ipart .eq. 5 ) then

                        k = 25

                     else if( ipart .ge. 12 .and. ipart .le. 13 ) then

                        k = 26

                     else

                        k = 26

                     end if

                     if( k .gt. 0 ) then

                       trEVENT(k,ir,iz)
     &                      = trEVENT(k,ir,iz) + enpar

                     end if

                     if( np .gt. 0 ) then

                        call pcheck(m,np,ipart,kpart,jpart,ipn,ips)

                        if( ipn .gt. 0 ) then

                           do i = 1, ipn

                              k = 29 + np + ips(i)

                              trEVENT(k,ir,iz)
     &                             = trEVENT(k,ir,iz)
     &                             + enpar

                           end do

                        end if

                     end if

                     if( kcoll .eq. 1 ) then

                       trEVENT(28,ir,iz) = trEVENT(28,ir,iz)
     &                      + enpar

                     end if

                  end if

               end if

*-----------------------------------------------------------------------

  600    continue

  300    continue

         end if

*-----------------------------------------------------------------------
*        total absorption energy
*-----------------------------------------------------------------------

         if( ( jcoll .lt. 6  .or. jcoll .eq. 10 .or.
     &         jcoll .eq. 12 .or. jcoll .eq. 16 .or. jcoll .eq. 17 .or.   ! S.Abe 2017/01/17
     &         jcoll .eq. 13 .or. jcoll .eq. 14 .or. jcoll .eq. 15 .or.   ! S.Abe 2017/01/17
     &         jcoll .eq. 11 ) .and.
     &       totout .gt. 0.0d0 ) then

            rabs(1) = rabs(1) + enmin - totout
            rabs(2) = rabs(2) + ( emathi - totmas ) * oldwt

         end if

*-----------------------------------------------------------------------

      end if

*-----------------------------------------------------------------------

      end

************************************************************************
*                                                                      *
      subroutine thetxyzEVENT(ncol,m,nd,ns0,nl,lt,np,nx,ny,nz,xm,ym,zm,
     &                   ne,eb,tr,rabs,trEVENT,tr0)
*                                                                      *
*       nuclear heating tally in xyz scoring mesh                      *
*       output = deposit-*                                             *
*       last modified by T.Furuta on 2012/05/08                        *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*      ncol  ..... reaction type                                       *
*              0 : end of batch                                        *
*              4 : source particle                                     *
*                                                                      *
*              9 : termination by time cut-off                         *
*             10 : geometry boundary crossing                          *
*             11 : termination by energy cut-off                       *
*             12 : termination by escape or leakage                    *
*             13 : nuclear reaction (n,x)                              *
*             14 : nuclear reaction (n,n'x)                            *
*             15 : sequential transport only for tally                 *
*                                                                      *
*     definition of tr( nd, ix, iy, iz, id )                           *
*                                                                      *
*        ns = 1 : all, 2 : simple, 3 : heat                            *
*                                                                      *
*      nd ..... number of event type                                   *
*        ns = 3                                                        *
*          1 : total                                                   *
*                                                                      *
*        ns = 2                                                        *
*          2 : cutoff neutron                                          *
*          3 : cutoff photon                                           *
*          4 : cutoff proton                                           *
*          5 : leakage                                                 *
*                                                                      *
*          6 : recoil                                                  *
*          7 : ionization                                              *
*          8 : low neutron                                             *
*          9 : electron                                                *
*         10 : low proton                                              *
*         11 : other                                                   *
*                                                                      *
*        ns = 1                                                        *
*         12 : deuteron recoil                                         *
*         13 : triton   recoil                                         *
*         14 : 3He      recoil                                         *
*         15 : alpha    recoil                                         *
*         16 : residual recoil                                         *
*                                                                      *
*         17 : proton ionization                                       *
*         18 : pion+  ionization                                       *
*         19 : pion-  ionization                                       *
*         20 : other  ionization                                       *
*                                                                      *
*         21 : stopped proton                                          *
*         22 : stopped neutron                                         *
*         23 : stopped photon                                          *
*         24 : stopped pion+                                           *
*         25 : stopped pion-                                           *
*         26 : stopped other                                           *
*                                                                      *
*         27 : remaining excitation energy                             *
*         28 : fission high recoil                                     *
*         29 : fission low neutron                                     *
*                                                                      *
*         29 + np      : ionization of other paticle                   *
*         29 + np + np : stopped other paticle                         *
*                                                                      *
*       id ..... data type                                             *
*          1 : heating                                                 *
*          2 : relative error                                          *
*                                                                      *
*       unit = MeV                                                     *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
************************************************************************
      use MMBANKMOD   !FURUTA
      use EVENTTALMOD !FURUTA
      use NGSDATAMOD, only : bindeg
      use partmod, only: itmxpt, itpan, itpat, jtpat ! frtati 2021/10/05
*-----------------------------------------------------------------------

      implicit real*8 (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param00.inc'
      include 'param.inc'

*-----------------------------------------------------------------------

      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
!$OMP THREADPRIVATE(/icomon/)
      common /wtsave/ oldwt
!$OMP THREADPRIVATE(/wtsave/)
      common /tlgeom/ iblz1,iblz2
!$OMP THREADPRIVATE(/tlgeom/)
      common /mathzn/ mathz, mathn, jcoll, kcoll
!$OMP THREADPRIVATE(/mathzn/)
      common /regdc/  idrg(kvlmax), idgr(kvmmax)
      common /kmat1d/ idmn(0:kvlmax), idnm(kvmmax)
      common /eparm/  esmax, esmin, emin(20)
      common /celdn/  denr(kvlmax), denm(kvlmax), denc(kvlmax)

*-----------------------------------------------------------------------

      common /clustt/ nclsts, iclusts(nnn)
!$OMP THREADPRIVATE(/clustt/)
      common /clustw/ jclusts(0:8,nnn),  qclusts(0:12,nnn)
!$OMP THREADPRIVATE(/clustw/)
      common /cntcls/ jcount(3,nnn)
!$OMP THREADPRIVATE(/cntcls/)

*-----------------------------------------------------------------------

      common /tall22/ itmcn(itlmax), itmtn(itlmax), itmtt(itlmax)
      common /tall32/ itelc(itlmax)
      common /tall37/ itcnt(9,itlmax)
      common /tall46/ itmtr(itlmax,4), rtmtr(itlmax,13)

*-----------------------------------------------------------------------

      dimension   lt(nl)
      dimension   xm(nx+1)
      dimension   ym(ny+1)
      dimension   zm(nz+1)
      dimension   tr(nd,nx,ny,nz,0:ne,2)
      dimension   rabs(5)
      dimension   eb(ne+1)
      dimension   trEVENT(nd,nx,ny,nz)

      dimension ips(itmxpt) ! frtati 2021/10/05 6 -> itmxpt

*-----------------------------------------------------------------------

      data rmngv / 0.93895 /
      data rmnmv / 938.95 /

*-----------------------------------------------------------------------

      dimension ud(3)

      data dmax  /1.0d+19/
      data dmax0 /1.0d+18/

*-----------------------------------------------------------------------

      save sweight
!$OMP THREADPRIVATE(sweight)
      common /jcomon/ nabov,nobch,nocas,nomax
!$OMP THREADPRIVATE(/jcomon/)
      common /tall36/ itdpo(itlmax)

*-----------------------------------------------------------------------

      common /paraj/ mstz(300), parz(300)

*-----------------------------------------------------------------------

      real*8 edep
      common/EPCONTedep/ edep
!$OMP THREADPRIVATE(/EPCONTedep/)

      integer iegsemi, iegsout
      common /egsemi/ iegsemi, iegsout

C for nonshared_tally option
      integer italsh
      common /talsh/ italsh

*------------------------------------------------------------------------

      common /spred/ nspred, nwsprd, nedisp, itstep, ndedx

*-----------------------------------------------------------------------
! sumover
      real(8) tr0(nd,nx,ny,nz,0:ne)

      tr0(:,:,:,:,:) = 0.0d0

            small = parz(28) * 1.d-4 ! T.Sato 2023/11/01 avoid infinite loop

*-----------------------------------------------------------------------
*     ns from ns0
*-----------------------------------------------------------------------

            ns = ns0
            if( ns .ge. 4 ) ns = ns - 3

*-----------------------------------------------------------------------
*     deposit : ncol = 0, 4 and ns = 4
*-----------------------------------------------------------------------

      if( ( ncol .eq. 0 .or. ncol .eq. 4 ) .and. ns0 .ge. 4 ) then

         if( .not. FIRSTsrc ) then

            do id = 1, nd
            do ix = 1, nx
            do iy = 1, ny
            do iz = 1, nz

               if( itdpo(m) .eq. 0 .and.
     &             trEVENT(1,ix,iy,iz) .gt. 0.0 ) then

                  heats = trEVENT(1,ix,iy,iz)
                  ratio = trEVENT(id,ix,iy,iz)
     &                 / trEVENT(1,ix,iy,iz)

               else

                  heats = trEVENT(id,ix,iy,iz)
                  ratio = 1.0d0

               end if

               if( heats .gt. 0.0d0 .and.
     &             heats .ge. eb(1) .and. heats .lt. eb(ne+1) ) then

                  do ie = 1, ne

                     if( heats .ge. eb(ie) .and.
     &                   heats .lt. eb(ie+1) ) then
! sumover
                       tr0(id,ix,iy,iz,ie) = tr0(id,ix,iy,iz,ie)
     &                                        + sweight*ratio
                       if (italsh .eq. 0 ) then
                         tr(id,ix,iy,iz,ie,1) = tr(id,ix,iy,iz,ie,1)
     &                                        + sweight*ratio
                         tr(id,ix,iy,iz,ie,2) = tr(id,ix,iy,iz,ie,2)
     &                                        + (sweight*ratio)**2
                       else
!$OMP CRITICAL (tr_thetxyzEVENT)
                        tr(id,ix,iy,iz,ie,1) = tr(id,ix,iy,iz,ie,1)
     &                                       + sweight*ratio
                        tr(id,ix,iy,iz,ie,2) = tr(id,ix,iy,iz,ie,2)
     &                                       + (sweight*ratio)**2
!$OMP END CRITICAL (tr_thetxyzEVENT)
                       end if

                     end if

                  end do

               end if

            end do
            end do
            end do
            end do

! sumover
           call thetxyz_sumover(m, 1,
     &                            nd, nx, ny, nz, ne, tr0)
           tr0(:,:,:,:,:) = 0.0d0

            do id = 1, nd
            do ix = 1, nx
            do iy = 1, ny
            do iz = 1, nz

               trEVENT(id,ix,iy,iz) = 0.0d0

            end do
            end do
            end do
            end do

         end if

*-----------------------------------------------------------------------
*        save initial weight and zero set
*-----------------------------------------------------------------------

         if( ncol .eq. 4 ) then

            sweight = oldwt

         end if

      end if

*-----------------------------------------------------------------------
*        check of ncol
*-----------------------------------------------------------------------

         if( ncol .lt. 9 ) return

*-----------------------------------------------------------------------
*        nothing for electron and positron with nuclear data
*                    and secondary photon
*-----------------------------------------------------------------------

         if( itelc(m) .eq. 0 .and. iegsemi .eq. 0 .and.   ! S.Abe 2015/08/31
     &     ( jcoll .eq. 8 .or.
     &     ( jcoll .eq. 7 .and. name(ibknam+no,ipomp+1) .gt. 1 ) ) )
     &                   return

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
*        transform positions
*-----------------------------------------------------------------------

            call trnsxx(x(ibkx+no,ipomp+1),y(ibky+no,ipomp+1),
     &                  z(ibkz+no,ipomp+1),
     &                  xxa,yya,zza,itmtr(m,4))

            call trnsxx(xc(ibkxc+no,ipomp+1),yc(ibkyc+no,ipomp+1),
     &                  zc(ibkzc+no,ipomp+1),
     &                  xxc,yyc,zzc,itmtr(m,4))

*-----------------------------------------------------------------------
* add dummy step for EGS5 electrons at ncol=11
*     (stopped but still have dE)
cc H.Iwase 2014/8/21 add the change of H.Iwase 2014/1/9
*-----------------------------------------------------------------------

         if( iegsemi .ne. 0 .and.
     &     ( ityp .eq. 12 .or. ityp .eq. 13 ) .and. ncol .eq. 11 ) then

               dum = ( xxc - xxa )**2
     &             + ( yyc - yya )**2
     &             + ( zzc - zza )**2

            if( dum .le. small**2 ) then

               dumstep = small * 10.d0

               xdum = x(ibkx+no,ipomp+1) + u(ibku+no,ipomp+1)*dumstep
               ydum = y(ibky+no,ipomp+1) + v(ibkv+no,ipomp+1)*dumstep
               zdum = z(ibkz+no,ipomp+1) + w(ibkw+no,ipomp+1)*dumstep

               call trnsxx(xdum,ydum,zdum,xxc,yyc,zzc,itmtr(m,4))

            end if

         end if
*-----------------------------------------------------------------------
*        check position : out of rainge
*-----------------------------------------------------------------------

            if( xxa-small .lt. xm(1) .and.
     &          xxc-small .lt. xm(1) ) return

            if( yya-small .lt. ym(1) .and.
     &          yyc-small .lt. ym(1) ) return

            if( zza-small .lt. zm(1) .and.
     &          zzc-small .lt. zm(1) ) return

            if( xxa+small .ge. xm(nx+1) .and.
     &          xxc+small .ge. xm(nx+1) ) return

            if( yya+small .ge. ym(ny+1) .and.
     &          yyc+small .ge. ym(ny+1) ) return

            if( zza+small .ge. zm(nz+1) .and.
     &          zzc+small .ge. zm(nz+1) ) return

*-----------------------------------------------------------------------
*        enmin : final particle energy ( with mass for not nucleon )
*-----------------------------------------------------------------------

               if( ityp .le. 2 .or. ityp .ge. 12 ) then
                  enmin = 0.0
               else
                  enmin = rtyp - mtyp * rmnmv
               end if

               enmin = ( ec(ibkec+no,ipomp+1) + enmin ) * oldwt


*-----------------------------------------------------------------------
*        check of counter
*-----------------------------------------------------------------------

            do i = 1, 3

               if( itcnt(i,m) .ne. 0 ) then

                  if( ncnt(ibknct+i,no,ipomp+1) .lt. itcnt(i*2+2,m) .or.
     &                ncnt(ibknct+i,no,ipomp+1) .gt. itcnt(i*2+3,m) )
     &                               goto 200

               end if

            end do

*-----------------------------------------------------------------------
*        ncol = 9, 11, 12 : time and enrgy cut-off, escape
*        reaction point is out of rainge -> 500
*-----------------------------------------------------------------------

      if( ncol .eq. 9 .or. ncol .eq. 11 .or. ncol .eq. 12 ) then

*-----------------------------------------------------------------------

               if( xxc .lt. xm(1) .or.
     &             xxc .ge. xm(nx+1) ) goto 500

               if( yyc .lt. ym(1) .or.
     &             yyc .ge. ym(ny+1) ) goto 500

               if( zzc .lt. zm(1) .or.
     &             zzc .ge. zm(nz+1) ) goto 500

*-----------------------------------------------------------------------
*           get end point : ix, iy, iz
*-----------------------------------------------------------------------

               do i = 1, nx

                  if( xxc .ge. xm(i) .and.
     &                xxc .lt. xm(i+1) ) goto 11

               end do

   11             ix = i

               do i = 1, ny

                  if( yyc .ge. ym(i) .and.
     &                yyc .lt. ym(i+1) ) goto 21

               end do

   21             iy = i

               do i = 1, nz

                  if( zzc .ge. zm(i) .and.
     &                zzc .lt. zm(i+1) ) goto 31

               end do

   31             iz = i

               if( ix .ge. nx + 1 .or.
     &             iy .ge. ny + 1 .or.
     &             iz .ge. nz + 1 ) goto 500

*-----------------------------------------------------------------------
*           stopped energy cut-off or time cut-off particles
*-----------------------------------------------------------------------

            if( ncol .eq. 9 .or. ncol .eq. 11 ) then

              trEVENT(1,ix,iy,iz) = trEVENT(1,ix,iy,iz)
     &             + enmin

               if( ns .le. 2 ) then

                  if( ityp .ne. 12 .and. ityp .ne. 13 ) then

                    trEVENT(11,ix,iy,iz) = trEVENT(11,ix,iy,iz)
     &                   + enmin

                  else

                    trEVENT(9,ix,iy,iz) = trEVENT(9,ix,iy,iz)
     &                   + enmin

                  end if

               end if

               if( ns .eq. 1 ) then

                  if( ityp .eq. 1 ) then

                        k = 21

                  else if( ityp .eq. 3 ) then

                        k = 24

                  else if( ityp .eq. 5 ) then

                        k = 25

                  else if( ityp .eq. 12 .or. ityp .eq. 13 ) then

                        k = 26
                        rabs(4) = rabs(4) + enmin

                  else

                        k = 26

                  end if

                  trEVENT(k,ix,iy,iz) = trEVENT(k,ix,iy,iz)
     &                 + enmin

                  if( np .gt. 0 ) then

                     call pcheck(m,np,ityp,ktyp,jtyp,ipn,ips)

                     if( ipn .gt. 0 ) then

                        do i = 1, ipn

                           k = 29 + np + ips(i)

                           trEVENT(k,ix,iy,iz)
     &                          = trEVENT(k,ix,iy,iz)
     &                          + enmin

                        end do

                     end if

                  end if

               end if

            end if

*-----------------------------------------------------------------------
*           termination by escape or leakage
*-----------------------------------------------------------------------

            if( ncol .eq. 12 ) then

               if( ns .le. 2 ) then

                 trEVENT(5,ix,iy,iz) = trEVENT(5,ix,iy,iz)
     &                + enmin

               end if

            end if

      end if

*-----------------------------------------------------------------------
*     ionization or track length for low energy neutron and photon
*-----------------------------------------------------------------------

  500 if( mat .le. 0 ) goto 200

*-----------------------------------------------------------------------
*     neutron, proton and photon for nuclear data
*-----------------------------------------------------------------------

      if( jcoll .eq. 6 .or.
     &    jcoll .eq. 9 .or.   ! S.Abe 2017/01/17
     &    ( ( jcoll .eq. 7 .or. jcoll .eq. 15 ) .and.   ! S.Abe 2017/01/17
     &      iegsemi .eq. 0 .and.   ! S.Abe 2015/08/31
     &  name(ibknam+no,ipomp+1) .eq. 1 .and. itelc(m) .eq. 0 ) ) then

*-----------------------------------------------------------------------

               icl = idgr(iblz1)
               erg = e(ibke+no,ipomp+1)

               mk = mat
               rh = denm(mat)

            if( ityp .eq. 2 ) then

               call heatn(icl,erg,heatr,heatf,mk,rh,m,0,mtdum)

            else if( ityp .eq. 14 ) then

               call heatp(icl,erg,heatr,mk,rh,m,0,mtdum)

               heatf = 0.0

            end if

            heatr0 = heatr * oldwt
            heatf0 = heatf * oldwt
            heat0  = heatr0 + heatf0

            inpd = 1

      else if( jcoll .eq. 9 ) then

            inpd = 2

      else

            inpd = 0

      end if

*-----------------------------------------------------------------------
*     only for charged particles
*-----------------------------------------------------------------------
cc H.Iwase 2014/1/9  (EGS5)

         if( inpd .eq. 0 ) then

            if( iegsemi .eq. 0 ) then

               if( jtyp .eq. 0 ) goto 200
               if( abs( ec(ibkec+no,ipomp+1) - e(ibke+no,ipomp+1) )
     &              .eq. 0.0d0 )
     &             goto 200

             else

               if( jtyp .eq. 0 .and. ityp .ne. 14 ) goto 200
               if( ( abs( ec(ibkec+no,ipomp+1) - e(ibke+no,ipomp+1) )
     &                          .eq. 0.0d0 )
     $              .and.
     $          ( ityp .ne. 12 .and. ityp .ne. 13 .and. ityp .ne. 14))
     $            goto 200

             end if

         end if

*-----------------------------------------------------------------------
*        distance and unit vector ud(i)
*-----------------------------------------------------------------------

            dis = ( xxc - xxa )**2
     &          + ( yyc - yya )**2
     &          + ( zzc - zza )**2

            if( dis .le. small**2 ) goto 200

            dis = sqrt(dis)

            ud(1) = ( xxc - xxa ) / dis
            ud(2) = ( yyc - yya ) / dis
            ud(3) = ( zzc - zza ) / dis

*-----------------------------------------------------------------------
*     initial position, energy and initial range
*-----------------------------------------------------------------------

            tot = 0.0d0

            se  = e(ibke+no,ipomp+1)
            xpp = xxa
            ypp = yya
            zpp = zza

            if( inpd .ne. 1 )
     &      call rainge(se,rng,mat,ityp,ktyp,jtyp,rtyp)

*-----------------------------------------------------------------------
*     initial xyz position
*-----------------------------------------------------------------------

               ixm = 0
               ixc = 0

            if( ud(1) .gt. 0.0d0 ) ixk =  1
            if( ud(1) .eq. 0.0d0 ) ixk =  0
            if( ud(1) .lt. 0.0d0 ) ixk = -1

            if( ixk .ge. 0 ) then

               do i = 1, nx + 1
                  if( xm(i) .gt. xpp + small ) then
                         ixm = i - 1
                         ixc = i - 1
                         goto 32
                  end if
               end do

            else

               ixm = nx + 2
               ixc = nx + 1

               do i = 1, nx + 1
                  if( xm(i) .ge. xpp - small ) then
                         ixm = i
                         ixc = i - 1
                         goto 32
                  end if
               end do

            end if

   32       continue

*-----------------------------------------------------------------------

               iym = 0
               iyc = 0

            if( ud(2) .gt. 0.0d0 ) iyk =  1
            if( ud(2) .eq. 0.0d0 ) iyk =  0
            if( ud(2) .lt. 0.0d0 ) iyk = -1

            if( iyk .ge. 0 ) then

               do i = 1, ny + 1
                  if( ym(i) .gt. ypp + small ) then
                         iym = i - 1
                         iyc = i - 1
                         goto 35
                  end if
               end do

            else

               iym = ny + 2
               iyc = ny + 1

               do i = 1, ny + 1
                  if( ym(i) .ge. ypp - small ) then
                         iym = i
                         iyc = i - 1
                         goto 35
                  end if
               end do

            end if

   35       continue

*-----------------------------------------------------------------------

               izm = 0
               izc = 0

            if( ud(3) .gt. 0.0d0 ) izk =  1
            if( ud(3) .eq. 0.0d0 ) izk =  0
            if( ud(3) .lt. 0.0d0 ) izk = -1

            if( izk .ge. 0 ) then

               do i = 1, nz + 1
                  if( zm(i) .gt. zpp + small ) then
                         izm = i - 1
                         izc = i - 1
                         goto 38
                  end if
               end do

            else

               izm = nz + 2
               izc = nz + 1

               do i = 1, nz + 1
                  if( zm(i) .ge. zpp - small ) then
                         izm = i
                         izc = i - 1
                         goto 38
                  end if
               end do

            end if

   38       continue

*-----------------------------------------------------------------------
*     loop for finding mesh and booking upto total distance
*-----------------------------------------------------------------------

   50 continue

*-----------------------------------------------------------------------
*     calculation is finished
*-----------------------------------------------------------------------

         if( tot .ge. dis ) goto 200

*-----------------------------------------------------------------------
*        dx,dy,dz : distance to the nearest xyz mesh
*-----------------------------------------------------------------------

            dx = dmax
            dy = dmax
            dz = dmax

         if( ixk .ne. 0 ) then
            if( ixm + ixk .ge. 1 .and.
     &         ixm + ixk .le. nx + 1 ) then

               dx = abs( ( xm(ixm+ixk) - xpp ) / ud(1) )

            end if
         end if

         if( iyk .ne. 0 ) then
            if( iym + iyk .ge. 1 .and.
     &          iym + iyk .le. ny + 1 ) then

               dy = abs( ( ym(iym+iyk) - ypp ) / ud(2) )

            end if
         end if

         if( izk .ne. 0 ) then
            if( izm + izk .ge. 1 .and.
     &          izm + izk .le. nz + 1 ) then

               dz = abs( ( zm(izm+izk) - zpp ) / ud(3) )

            end if
         end if

*-----------------------------------------------------------------------
*        which boundary is the nearlist
*-----------------------------------------------------------------------

               dd = dx
               jj = 1

            if( dy .lt. dd ) then

               dd = dy
               jj = 2

            end if

            if( dz .lt. dd ) then

               dd = dz
               jj = 3

            end if

            if( dd .gt. dmax0 ) goto 200

*-----------------------------------------------------------------------
*        recalculate rng for ncol = 11 and nedisp
*-----------------------------------------------------------------------

            if( mat .gt. 0 .and. jtyp .ne. 0 .and. nedisp .ne. 0 .and.
     &           (  ityp .ne. 12 .and. ityp .ne. 13 ) ) then


               call rainge(se,rng,mat,ityp,ktyp,jtyp,rtyp)

               if( ncol .eq. 11 )  dis = tot + rng

            end if

*-----------------------------------------------------------------------
*        propagate position upto the boundary or the final point
*-----------------------------------------------------------------------

            if( tot + dd .ge. dis - small ) then

               dd  = dis - tot
               tot = dis
               xpp = xxc
               ypp = yyc
               zpp = zzc

            else

               tot = tot + dd
               xpp = xpp + dd * ud(1)
               ypp = ypp + dd * ud(2)
               zpp = zpp + dd * ud(3)

            end if

*-----------------------------------------------------------------------
*        final energy ee and range rng
*-----------------------------------------------------------------------

            if( inpd .ne. 1 ) then

               delt = dd

               call ecol(ee,delt,se,rng,
     &                   mat,ityp,ktyp,jtyp,rtyp)

               rng = max( 0.0d0, rng - delt )
               ee  = max( 0.0d0, ee )
               ee  = min( se, ee )

               enion = abs( se - ee ) * oldwt

               if( iegsemi .ne. 0 .and. ityp .eq. 14 ) then
                  enion = enion + edep * oldwt
               end if

            else

               ee = se

               heatr = heatr0 * dd
               heatf = heatf0 * dd
               heat  = heat0  * dd

                  if( ityp .eq. 2 ) then

                     rabs(3) = rabs(3) + heat

                  else

                     rabs(4) = rabs(4) + heat

                  end if

            end if

*-----------------------------------------------------------------------
*        booking the ionization energy
*-----------------------------------------------------------------------

         if( ixc .ge. 1 .and. ixc .lt. nx + 1 .and.
     &       iyc .ge. 1 .and. iyc .lt. ny + 1 .and.
     &       izc .ge. 1 .and. izc .lt. nz + 1 ) then

            if( inpd .ne. 1 ) then

              trEVENT(1,ixc,iyc,izc) = trEVENT(1,ixc,iyc,izc)
     &             + enion

               if( ns .le. 2 ) then

                 trEVENT(7,ixc,iyc,izc) = trEVENT(7,ixc,iyc,izc)
     &                + enion

               end if

               if( ns .eq. 1 ) then

                  if( ityp .eq. 1 ) then

                        k = 17

                  else if( ityp .eq. 3 ) then

                        k = 18

                  else if( ityp .eq. 5 ) then

                        k = 19

                  else if( ityp .eq. 12 .or. ityp .eq. 13 ) then

                        k = 20
                        rabs(4) = rabs(4) + enion

                  else

                        k = 20

                  end if

                  trEVENT(k,ixc,iyc,izc)
     &                 = trEVENT(k,ixc,iyc,izc)
     &                 + enion

                  if( np .gt. 0 ) then

                     call pcheck(m,np,ityp,ktyp,jtyp,ipn,ips)

                     if( ipn .gt. 0 ) then

                        do i = 1, ipn

                           k = 29 + ips(i)

                           trEVENT(k,ixc,iyc,izc)
     &                          = trEVENT(k,ixc,iyc,izc)
     &                          + enion

                        end do

                     end if

                  end if

               end if

            end if

*-----------------------------------------------------------------------
*           neutron and photon for nuclear data
*-----------------------------------------------------------------------

            if( inpd .eq. 1 ) then

              trEVENT(1,ixc,iyc,izc) = trEVENT(1,ixc,iyc,izc)
     &             + heat

                if( ns .le. 2 ) then

                     if( ityp .eq. 2 ) then

                       trEVENT(8,ixc,iyc,izc)
     &                      = trEVENT(8,ixc,iyc,izc)
     &                      + heat

                     else

                       trEVENT(9,ixc,iyc,izc)
     &                      = trEVENT(9,ixc,iyc,izc)
     &                      + heat

                     end if

                end if

                if( ns .eq. 1 ) then

                  trEVENT(29,ixc,iyc,izc)
     &                 = trEVENT(29,ixc,iyc,izc)
     &                 + heatf

                end if

*-----------------------------------------------------------------------
*           proton for nuclear data
*-----------------------------------------------------------------------

            else if( inpd .eq. 2 ) then

                     erg = ( se + ee ) / 2.0
                     icl = idgr(iblz1)

                     mk = mat
                     rh = denm(mat)

                     call heath(icl,erg,heatr,mk,rh,m,0,mtdum)

                     heatf = 0.0

                     heatr0 = heatr * oldwt
                     heatf0 = heatf * oldwt
                     heat0  = heatr0 + heatf0

                     heatr = heatr0 * dd
                     heatf = heatf0 * dd
                     heat  = heat0  * dd

                     rabs(5) = rabs(5) + heat

                     trEVENT(1,ixc,iyc,izc)
     &                    = trEVENT(1,ixc,iyc,izc)
     &                    + heat

                if( ns .le. 2 ) then

                  trEVENT(10,ixc,iyc,izc)
     &                 = trEVENT(10,ixc,iyc,izc)
     &                 + heat

                end if

            end if

         end if

*-----------------------------------------------------------------------
*     next position
*-----------------------------------------------------------------------

            if( jj .eq. 1 ) then

               ixm = ixm + ixk
               ixc = ixc + ixk

            else if( jj .eq. 2 ) then

               iym = iym + iyk
               iyc = iyc + iyk

            else if( jj .eq. 3 ) then

               izm = izm + izk
               izc = izc + izk

            end if

            if( inpd .ne. 1 ) then

               se  = ee

            end if

*-----------------------------------------------------------------------

      goto 50

*-----------------------------------------------------------------------
*     produced particles
*-----------------------------------------------------------------------

  200 continue

*-----------------------------------------------------------------------
*     energy of produced particle and/or nucleus from nuclear reactions
*-----------------------------------------------------------------------
cc H.Iwase 2014/1/9  (EGS5)

      if( ( ncol .eq. 13 .or. ncol .eq. 14 ) .and.
     &    ( jcoll .lt. 6  .or. jcoll .eq. 10 .or.
     &      jcoll .eq. 12 .or. jcoll .eq. 16 .or. jcoll .eq. 17 .or.   ! S.Abe 2017/01/17
     &      jcoll .eq. 13 .or. jcoll .eq. 14 .or. jcoll .eq. 15 .or.   ! S.Abe 2017/01/17
     &      jcoll .eq. 11 .or.
     &      ( itelc(m) .eq. 1 .and.
     &        ( jcoll .eq. 7 .or. jcoll .eq. 8 ) ) ) ) then

                  totout = 0.0
                  totmas = 0.0

                  emathi = bindeg(mathz,mathn)

               if( ityp .ge. 15 .and. ityp .le. 19 ) then

                  ipnm = ktyp / 1000000
                  innm = ktyp - ktyp / 1000000 * 1000000 - ipnm

                  emathi = emathi + bindeg(ipnm,innm)

               end if

*-----------------------------------------------------------------------
*        from ordinary output
*-----------------------------------------------------------------------

         if( nclsts .gt. 0 ) then

*-----------------------------------------------------------------------
*        check of counter
*-----------------------------------------------------------------------

            do i = 1, 3

               if( itcnt(i,m) .ne. 0 ) then

                  if( jcount(i,1) .lt. itcnt(i*2+2,m) .or.
     &                jcount(i,1) .gt. itcnt(i*2+3,m) ) goto 300

               end if

            end do

*-----------------------------------------------------------------------

         do 600 j = 1, nclsts

               ipart = jclusts(3,j)
               kpart = jclusts(7,j)
               jpart = jclusts(5,j)

*-----------------------------------------------------------------------
*           collision position
*-----------------------------------------------------------------------

               call trnsxx(qclusts(10,j),
     &                     qclusts(11,j),
     &                     qclusts(12,j),
     &                     xcc,ycc,zcc,itmtr(m,4))

               if( xcc .lt. xm(1) .or.
     &             xcc .ge. xm(nx+1) ) goto 600

               if( ycc .lt. ym(1) .or.
     &             ycc .ge. ym(ny+1) ) goto 600

               if( zcc .lt. zm(1) .or.
     &             zcc .ge. zm(nz+1) ) goto 600

*-----------------------------------------------------------------------
*           get end point : ix, iy, iz
*-----------------------------------------------------------------------

               do i = 1, nx

                  if( xcc .ge. xm(i) .and.
     &                xcc .lt. xm(i+1) ) goto 41

               end do

   41             ix = i

               do i = 1, ny

                  if( ycc .ge. ym(i) .and.
     &                ycc .lt. ym(i+1) ) goto 51

               end do

   51             iy = i

               do i = 1, nz

                  if( zcc .ge. zm(i) .and.
     &                zcc .lt. zm(i+1) ) goto 61

               end do

   61             iz = i

               if( ix .ge. nx + 1 .or.
     &             iy .ge. ny + 1 .or.
     &             iz .ge. nz + 1 ) goto 600

*-----------------------------------------------------------------------
*              nuclei
*-----------------------------------------------------------------------

               if( ipart .ge. 15 .and. ipart .le. 19 ) then

                     enrem  = qclusts(6,j) * qclusts(8,j)

                     totout = totout + enrem
                     totmas = totmas + bindeg(jclusts(1,j),jclusts(2,j))

*-----------------------------------------------------------------------
*                 enrem is remaining excitation energies,
*-----------------------------------------------------------------------

                     trEVENT(1,ix,iy,iz) = trEVENT(1,ix,iy,iz)
     &                    + enrem

                  if( ns .le. 2 ) then

                    trEVENT(11,ix,iy,iz)
     &                   = trEVENT(11,ix,iy,iz)
     &                   + enrem

                  end if

                  if( ns .eq. 1 ) then

                    trEVENT(27,ix,iy,iz)
     &                   = trEVENT(27,ix,iy,iz)
     &                   + enrem

                     if( kcoll .eq. 1 ) then

                       trEVENT(28,ix,iy,iz)
     &                      = trEVENT(28,ix,iy,iz)
     &                      + enrem

                     end if

                  end if

               end if

*-----------------------------------------------------------------------
*              out going particles and nuclei
*-----------------------------------------------------------------------

                  if( ipart .le. 2 .or. ipart .ge. 12 ) then
                     enoms = 0.0
                  else
                     enoms = ( qclusts(5,j) - jclusts(6,j) * rmngv )
                  end if

                     enout = qclusts(7,j)
                     enmas = enoms * 1000.0
                     enpar = ( enout + enmas ) * qclusts(8,j)

                     totout = totout + enpar

*-----------------------------------------------------------------------
*              dead particle and nucleus
*-----------------------------------------------------------------------

               if( jclusts(4,j) .lt. 0 ) then

                     if( jclusts(4,j) .eq. -1 ) then

                       trEVENT(1,ix,iy,iz)
     &                      = trEVENT(1,ix,iy,iz)
     &                      + enpar

                     end if

                  if( ns .le. 2 ) then

                     if( ipart .eq. 2 .and.
     &                   jclusts(4,j) .eq. -2 ) then

                        k = 2

                     else if( ipart .ge. 12 .and. ipart .le. 14 .and.
     &                        jclusts(4,j) .eq. -2 ) then

                        k = 3

                     else if( ipart .eq. 1 .and.
     &                        jclusts(4,j) .eq. -2 ) then

                        k = 4

                     else if( ipart .ge. 15 .and. ipart .le. 19 ) then

                        k = 6

                     else if( ipart .eq. 12 .or. ipart .eq. 13 ) then

                        k = 9
                        rabs(4) = rabs(4) + enpar

                     else if( ipart .eq. 14 .and. ns .eq. 2 ) then

                        k = 11

                     else if( ipart .eq. 14 .and. ns .eq. 1 ) then

                        k = 0

                     else

                        k = 11

                     end if

                     if( k .gt. 0 ) then

                       trEVENT(k,ix,iy,iz)
     &                      = trEVENT(k,ix,iy,iz)
     &                      + enpar

                     end if

                  end if

                  if( ns .eq. 1 .and. jclusts(4,j) .eq. -1 ) then

                     if( ipart .eq. 1 ) then

                        k = 21

                     else if( ipart .eq. 2 ) then

                        k = 22

                     else if( ipart .eq. 14 ) then

                        k = 23

                     else if( ipart .ge. 15 .and.
     &                        ipart .le. 19) then

                        k = ipart - 3

                     else if( ipart .eq. 3 ) then

                        k = 24

                     else if( ipart .eq. 5 ) then

                        k = 25

                     else if( ipart .ge. 12 .and. ipart .le. 13 ) then

                        k = 26

                     else

                        k = 26

                     end if

                     if( k .gt. 0 ) then

                        trEVENT(k,ix,iy,iz)
     &                      = trEVENT(k,ix,iy,iz)
     &                      + enpar

                     end if

                     if( np .gt. 0 ) then

                        call pcheck(m,np,ipart,kpart,jpart,ipn,ips)

                        if( ipn .gt. 0 ) then

                           do i = 1, ipn

                              k = 29 + np + ips(i)

                              trEVENT(k,ix,iy,iz)
     &                             = trEVENT(k,ix,iy,iz)
     &                             + enpar

                           end do

                        end if

                     end if

                     if( kcoll .eq. 1 ) then


                       trEVENT(28,ix,iy,iz)
     &                      = trEVENT(28,ix,iy,iz)
     &                      + enpar

                     end if

                  end if

               end if

*-----------------------------------------------------------------------

  600    continue

  300    continue

         end if

*-----------------------------------------------------------------------
*        total absorption energy
*-----------------------------------------------------------------------

         if( ( jcoll .lt. 6  .or. jcoll .eq. 10 .or.
     &         jcoll .eq. 12 .or. jcoll .eq. 16 .or. jcoll .eq. 17 .or.   ! S.Abe 2017/01/17
     &         jcoll .eq. 13 .or. jcoll .eq. 14 .or. jcoll .eq. 15 .or.   ! S.Abe 2017/01/17
     &         jcoll .eq. 11) .and.
     &       totout .gt. 0.0d0 ) then

            rabs(1) = rabs(1) + enmin - totout
            rabs(2) = rabs(2) + ( emathi - totmas ) * oldwt

         end if

*-----------------------------------------------------------------------

      end if

*-----------------------------------------------------------------------

      end

************************************************************************
*                                                                      *
      subroutine tdepstregEVENT(ncol,m,nl,lt,np,nr,mr,ne,nt,kr,eb,tb,tr,
     &     trEVENT,tr0)
*                                                                      *
*       Deposit tally in region mesh                                   *
*       output = deposit                                               *
*       last modified by T.Furuta on 2012/05/08                        *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*      ncol  ..... reaction type                                       *
*              0 : end of batch                                        *
*              4 : source particle                                     *
*                                                                      *
*              9 : termination by time cut-off                         *
*             10 : geometry boundary crossing                          *
*             11 : termination by energy cut-off                       *
*             12 : termination by escape or leakage                    *
*             13 : nuclear reaction (n,x)                              *
*             14 : nuclear reaction (n,n'x)                            *
*             15 : sequential transport only for tally                 *
*                                                                      *
************************************************************************
      use MMBANKMOD   !FURUTA
      use EVENTTALMOD !FURUTA
      use tdepwgtsum_global !cABE 2016/11/24
      use partmod, only: itmxpt ! frtati 2021/10/05
*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param00.inc'
      include 'param.inc'

*-----------------------------------------------------------------------

      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
!$OMP THREADPRIVATE(/icomon/)
      common /wtsave/ oldwt
!$OMP THREADPRIVATE(/wtsave/)
      common /tlgeom/ iblz1,iblz2
!$OMP THREADPRIVATE(/tlgeom/)
      common /tlglat/ ilev1,ilev2,ilat1(5,10),ilat2(5,10)
!$OMP THREADPRIVATE(/tlglat/)

*-----------------------------------------------------------------------

      common /clustt/ nclsts, iclusts(nnn)
!$OMP THREADPRIVATE(/clustt/)
      common /clustw/ jclusts(0:8,nnn),  qclusts(0:12,nnn)
!$OMP THREADPRIVATE(/clustw/)
      common /cntcls/ jcount(3,nnn)
!$OMP THREADPRIVATE(/cntcls/)

*-----------------------------------------------------------------------

      common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                itout(itlmax), ittwo(itlmax)
      common /tall22/ itmcn(itlmax), itmtn(itlmax), itmtt(itlmax)
      common /tall36/ itdpo(itlmax)   ! S.Abe 2015/12/03
      common /tall37/ itcnt(9,itlmax)

      common /tall50/ itlmt(itlmax), ite2l(itlmax)   ! CCSE 2022.08.31
      common /tall51/ eletb(2000), nlete
      common /tall53/ itdfn(itlmax,2)

c T.Sato 2014/8/19 for detector resolution
      common /tall61/ rtdre(itlmax),rtdfa(itlmax)

      common /kmat1d/ idmn(0:kvlmax), idnm(kvmmax)

cABE 2016/03/04 for identification part=all
      common /tall65/ ipall(itlmax)

      common /tall70/ itrwgtsum(itlmax), itrgn1(itlmax),
     &                itrgm1(itlmax), itrncd(itlmax)

      common /tall00/ itnm, ital(itlmax), itals(itlmax), italm(itlmax)
      common /tall75/ itnlatcel(itlmax), itnlatmem(itlmax)

      common /tall82/ itcnth(9,itlmax)

      integer,save:: inilat
!$OMP THREADPRIVATE( inilat )
      data inilat /0/
      integer,allocatable,save:: nldtct(:)
!$OMP THREADPRIVATE( nldtct )
      integer,allocatable,save:: nlats(:,:)
!$OMP THREADPRIVATE( nlats )
      integer,allocatable,save:: nlatt(:,:)
!$OMP THREADPRIVATE( nlatt )
      integer,allocatable,save:: nlatu(:,:)
!$OMP THREADPRIVATE( nlatu )
      integer,save:: inilaterr
      data inilaterr /0/

*-----------------------------------------------------------------------

      dimension   lt(nl)
      dimension   kr(mr)
      dimension   eb(ne+1)
      dimension   tb(nt+1)
      dimension   tr(np,0:ne,nr,nt,2)
      dimension   trEVENT(np,nr,nt,itnlatmem(m))   ! S.Abe 2018/10/25, add itnlatmem(m)

      dimension ips(itmxpt) ! frtati 2021/10/05 6 -> itmxpt
      dimension facm(6)

*-----------------------------------------------------------------------

      double precision, save :: sweight(itlmax)
!$OMP THREADPRIVATE(sweight)

      common /jcomon/ nabov,nobch,nocas,nomax
!$OMP THREADPRIVATE(/jcomon/)

      common /regdc/  idrg(kvlmax), idgr(kvmmax)
      common /  tscmsg   / ktsc(kvlmax), mntsc, ntsc(kvlmax)
      common /  tscreg   / ntscell(kvlmax)
      common / etsminmax / etsmin, etsmax
      common / ptsminmax / ptsmin, ptsmax
      common / ctsminmax / ctsmin, ctsmax

      common /eparm/  esmax, esmin, emin(20)

      common / tsminmax  / tsmax
      common / etsexe  / dexc_ene
!$OMP THREADPRIVATE(/etsexe/)
      common /celdg/  rhog(kvlmax)

*-----------------------------------------------------------------------

      real*8 edep
      common/EPCONTedep/ edep
!$OMP THREADPRIVATE(/EPCONTedep/)

      integer iegsemi, iegsout
      common /egsemi/ iegsemi, iegsout

C for nonshared_tally option
      integer italsh
      common /talsh/ italsh

      common /dedxfac/ dedxfd
!$OMP THREADPRIVATE(/dedxfac/)

      common /cntmx/ ncntmx(3), nctmxsr(3)
!$OMP THREADPRIVATE(/cntmx/)

! sumover
      real(8) tr0(np,0:ne,nr,nt)

*-----------------------------------------------------------------------
*     allocate and initialze nlats,t,u
*-----------------------------------------------------------------------

! sumover
      tr0(:,:,:,:) = 0.0d0

      rng1=0.0 ! T.Sato 2024/04/02

      if( inilat .eq. 0 ) then

         inilat = 1

         latmem = 1
         do i = 1, itnm
            latmem = max(latmem, itnlatmem(m))
         enddo

         allocate( nldtct(itlmax) )
         nldtct(:)=0
         allocate( nlats(itlmax,latmem) )
         allocate( nlatt(itlmax,latmem) )
         allocate( nlatu(itlmax,latmem) )
         nlats(:,:)=0
         nlatt(:,:)=0
         nlatu(:,:)=0

      endif

*-----------------------------------------------------------------------
! T.Sato 2014/8/21, for considering detector resolution, ncol = 101
      if(ncol.eq.101) then
       if(rtdre(m).eq.0.0.and.rtdfa(m).eq.0.0) return ! detector resolution = 0
       if( .not. FIRSTsrc ) then
        do il = 1, nldtct(m)   ! S.Abe 2018/10/25
        do ir = 1, nr
         do it = 1, nt
          gaurntmp=gaurn(dummy)  ! T.Sato 2018/01/25
          do ip = 1, np
           if( trEVENT(ip,ir,it,il).gt.0.0d0 ) then

            heats = trEVENT(ip,ir,it,il)
            if(rtdre(m) .lt. 0.d0) then ! Ogawa 2018/11/05 if sigma<0, gauss + exp asymmetric resolution
             heats = usrdefres(heats)
            else
             heats = heats + gaurntmp*  ! T.Sato 2018/01/25
     &       sqrt(rtdre(m)**2+rtdfa(m)*heats)
            endif

            if(heats.lt.0.0) heats = 0.0
            trEVENT(ip,ir,it,il)=heats
           endif
          enddo
         enddo
        enddo
        enddo
       endif
       return
      endif
*----- End of revision on 2014/8/21 ------------------------------------

*-----------------------------------------------------------------------
*     ncol = 0, 4 and reg = weightsum
*-----------------------------------------------------------------------

      if( (ncol.eq.0 .or. ncol.eq.4) .and. itrwgtsum(m).eq.1 ) then

       if( .not. FIRSTsrc ) then

*-----------------------------------------------------------------------
*       check the event satisfy conditions or not
*-----------------------------------------------------------------------

        do il = 1, nldtct(m)   ! S.Abe 2018/10/25
        do it = 1, nt

         nreg = 0
         nlist = 0

         do j = 1, ncond_old

          if( itcond(m,1,j,1) .eq. 0 ) cycle

          nagree = 1

          do k = 1, nadd_old

           if( itcond(m,1,j,k) .eq. 0 ) exit

           nreg = nreg + 1
           edep = trEVENT(ipall(m),nreg,it,il)

           iope = itcond(m,2,j,k)
           eth = rteth(m,j,k)

           if( ( iope.eq.1 .and. edep.lt.eth ) .or.
     &         ( iope.eq.2 .and. edep.le.eth ) .or.
     &         ( iope.eq.3 .and. edep.eq.eth ) .or.
     &         ( iope.eq.4 .and. edep.ge.eth ) .or.
     &         ( iope.eq.5 .and. edep.gt.eth ) ) then
            nagree = nagree
           else
            nagree = 0
           endif

          enddo

          if( nagree .eq. 1 ) then
           nlist = itcond(m,3,j,1)
           exit
          endif

         enddo

*-----------------------------------------------------------------------
*        multiply efficiency
*-----------------------------------------------------------------------

         lchk = 0

         do ilst = 1, nlist_old
          if( itlist(m,ilst) .eq. nlist ) then
           lchk = ilst
           exit
          endif
         enddo

         do ip = 1, np
          heatwgtsum = 0.d0
          do ir = itrgn1(m)+1, nr
           trEVENT(ip,ir,it,il) = trEVENT(ip,ir,it,il)
     &                        * rteff(m,ir-itrgn1(m),lchk)
           heatwgtsum = heatwgtsum + trEVENT(ip,ir,it,il)
          enddo
          trEVENT(ip,1,it,il) = heatwgtsum
         enddo

*-----------------------------------------------------------------------

        enddo
        enddo

       endif

      endif

*-----------------------------------------------------------------------
* check of history counter
*-----------------------------------------------------------------------
         if ( ncol.eq.0 .or. ncol.eq.4 ) then
           do i = 1, 3
             if( itcnth(i,m) .eq. 1 ) then
               if( ncntmx(i) .lt. itcnth(i*2+2,m) .or.
     &             ncntmx(i) .gt. itcnth(i*2+3,m) ) then
                 trEVENT = 0.d0
               end if
             end if
           end do
         end if

*-----------------------------------------------------------------------
*     deposit : ncol = 0, 4 and ns = 4
*-----------------------------------------------------------------------

      if( ( ncol .eq. 0 .or. ncol .eq. 4 ) .and.
     &      itout(m) .ge. 2 ) then

         if( .not. FIRSTsrc ) then

            if( itrwgtsum(m) .eq. 1 ) then
             nr0 = 1
            else
             nr0 = nr
            endif


            do il = 1, nldtct(m)   ! S.Abe 2018/10/25
            do ip = 1, np
            do ir = 1, nr0
            do it = 1, nt

               if( itdpo(m) .eq. 0 ) then
                if( trEVENT(ipall(m),ir,it,il) .gt. 0.0 ) then
                  heats = trEVENT(ipall(m),ir,it,il)
                  ratio = trEVENT(ip,ir,it,il)
     &                 / trEVENT(ipall(m),ir,it,il)
                else
                  heats = trEVENT(ip,ir,it,il)
                  ratio = 1.0d0
                endif                  
               else
                  heats = trEVENT(ip,ir,it,il)
                  ratio = 1.0d0
               end if

               if( heats .gt. 0.0d0 .and.
     &             heats .ge. eb(1) .and. heats .lt. eb(ne+1) ) then

                  do ie = 1, ne

                     if( heats .ge. eb(ie) .and.

     &                   heats .lt. eb(ie+1) ) then
!$OMP CRITICAL (tr_tdepstregEVENT_1)
! sumover
                        tr0(ip,ie,ir,it) =
     &                         tr0(ip,ie,ir,it) + sweight(m)*ratio
                        tr(ip,ie,ir,it,1) = tr(ip,ie,ir,it,1)
     &                                    + sweight(m)*ratio   ! S.Abe 2015/12/03
!$OMP END CRITICAL (tr_tdepstregEVENT_1)
                        if(ip.eq.ipall(m).and.ratio.gt.0) then ! T.Sato 2022/09/09, part=all is always probability mode
!$OMP CRITICAL (tr_tdepstregEVENT_2)
                         tr(ip,ie,ir,it,2) = tr(ip,ie,ir,it,2)
     &                                    + 1.d0            ! Ogawa 2020/04/22  if (weight < 1.d0), tr(ip,ie,ir,it,2) should be count number
!$OMP END CRITICAL (tr_tdepstregEVENT_2)
                        else
!$OMP CRITICAL (tr_tdepstregEVENT_3)
                         tr(ip,ie,ir,it,2) = tr(ip,ie,ir,it,2)
     &                                    + (sweight(m)*ratio)**2   ! S.Abe 2015/12/03
!$OMP END CRITICAL (tr_tdepstregEVENT_3)
                        endif

                     end if

                  end do

               end if

            end do
            end do
            end do
            end do

! sumover
           call tdepstreg_sumover_ip(m,1, ipall(m), ratio,
     &                               np,  ne, nr,  nt, tr0)
!           endif
           tr0(:,:,:,:) = 0.0d0


            do il = 1, nldtct(m)   ! S.Abe 2018/10/25
            do ip = 1, np
            do ir = 1, nr
            do it = 1, nt
               trEVENT(ip,ir,it,il) = 0.0d0
            end do
            end do
            end do
            end do

         end if

*-----------------------------------------------------------------------
*        save initial weight and zero set
*-----------------------------------------------------------------------

         if( ncol .eq. 4 ) then

            sweight(m) = 1.d30 ! initialize

            nldtct(m) = 0   ! S.Abe 2018/10/25

         end if

      end if

*-----------------------------------------------------------------------
*        check of ncol
*-----------------------------------------------------------------------

         if( ncol .lt. 9 ) return

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
*        check region
*-----------------------------------------------------------------------

               jj = 0

           do ii = 1, nr

               call tregck(iblz1,ilev1,ilat1,mr,kr,jj,icc)

               if( icc .ne. 0 ) goto 501

           end do

           return

  501      continue

*-----------------------------------------------------------------------
*        check time
*-----------------------------------------------------------------------

               tin = abs(t(ibkt+no,ipomp+1))

               if( tin .lt. tb(1) ) return
               if( tin .ge. tb(nt+1) ) return

            do i = 1, nt

               if( tin .ge. tb(i) .and.
     &             tin .lt. tb(i+1) ) goto 50

            end do

   50          it = i

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
*        check lattice position
*-----------------------------------------------------------------------

            if( itnlatcel(m) .eq. 0 ) then

               ilatflag = -1
               nldtct(m) = 1
               il = nldtct(m)

            elseif( itnlatcel(m) .gt. 0 ) then

               if( ilev1 .le. 0 ) return

               do jl = 1, ilev1
                  if( ilat1(1,jl) .eq. itnlatcel(m) ) goto 1000
               enddo
               return

 1000          continue

               ilats = ilat1(3,jl)
               ilatt = ilat1(4,jl)
               ilatu = ilat1(5,jl)

               ilatflag = 0
               if( nldtct(m) .gt. 0 ) then
                  do il = 1, nldtct(m)
                     if( ilats .eq. nlats(m,il) .and.
     &                   ilatt .eq. nlatt(m,il) .and.
     &                   ilatu .eq. nlatu(m,il) ) then
                        ilatflag = -1
                        exit
                     endif
                  enddo
               endif

               if( ilatflag .eq. 0 ) then

                  il = nldtct(m) + 1
                  if( il .gt. itnlatmem(m) ) then
                     if( inilaterr .eq. 0 ) then
                        write(*,'('' **** Warning: nlatmem overflowed'',
     &                            '' in [t-deposit]'')')
                        write(*,'('' please increase nlatmem, or'',
     &                 '' some deposited events are ignored.'')')

                        inilaterr = -1
                     endif
                     return
                  endif

                  nlats(m,il) = ilats
                  nlatt(m,il) = ilatt
                  nlatu(m,il) = ilatu

               endif

            endif

*-----------------------------------------------------------------------
*           track length
*-----------------------------------------------------------------------

               tlngth = sqrt( ( xc(ibkxc+no,ipomp+1) -
     &                          x(ibkx+no,ipomp+1) )**2
     &                      + ( yc(ibkyc+no,ipomp+1) -
     &                          y(ibky+no,ipomp+1) )**2
     &                      + ( zc(ibkzc+no,ipomp+1) -
     &                          z(ibkz+no,ipomp+1) )**2 )

*-----------------------------------------------------------------------
*           initial and final energy
*-----------------------------------------------------------------------

                  ein = e(ibke+no,ipomp+1)
                  ecc = ec(ibkec+no,ipomp+1)

                  ecc = ein - dedxfd * ( ein - ecc )   ! S.Abe 2015/09/08

                  if( ncol .eq. 9 .or. ncol .eq. 11 ) ecc = 0.d0

*-----------------------------------------------------------------------
cKN Iwase 2014/08/22 for EGS

               if( ( ityp.eq.12 .or. ityp.eq.13 .or. ityp.eq.14 ) .and.
     &               iegsemi .ne. 0 ) then

                     call egs5edxde(ecc,tlngth,ein,mat,ityp)

                  if(ncol.eq.9 .or. ncol.eq.11) then
                     ecc = 0d0
                     if( ityp.eq.12 .or. ityp.eq.13 ) ein = ein + edep
                  endif

               endif

*-----------------------------------------------------------------------
*        dedxfnc is specified
*-----------------------------------------------------------------------

         if( itdfn(m,1) .ne. 0 ) then

               do i = 1, nlete

                  if( ein .ge. eletb(i) .and.
     &                ein .lt. eletb(i+1) ) goto 30

               end do

   30             ie1 = min( i, nlete )

               do i = nlete, 1, -1

                  if( ecc .ge. eletb(i) .and.
     &                ecc .lt. eletb(i+1) ) goto 40

               end do

   40             ie2 = max( i, 1 )

*-----------------------------------------------------------------------
*           rrl : track length and rrr : charge particles correction
*-----------------------------------------------------------------------

               rrl = tlngth
               rrr = tlngth
               rrt = 0.0d0

            if( mat .gt. 0 .and. jtyp .ne. 0 ) then

                  call rainge(ein,rng0,mat,ityp,ktyp,jtyp,rtyp)

                  rng1 = max( 0.0d0, rng0 - tlngth )

                  if( rng1 .gt. 1.d+30 ) return

               if( ein .gt. eletb(nlete+1) ) then

                  call rainge(eletb(nlete+1),rngm,mat,
     &                        ityp,ktyp,jtyp,rtyp)

                  rrr = max( 0.0d0, rrr - ( rng0 - rngm ) )
                  rrt = rrt + rng0 -rngm

               end if

            end if

               rrl0 = rrl
               rrr0 = rrr
               rrt0 = rrt

         end if

*-----------------------------------------------------------------------
*     check of region
*-----------------------------------------------------------------------

                  jj = 0
                  ein0 = ein
                  ecc0 = ecc

*-----------------------------------------------------------------------

      do 100 ii = 1, nr
               ein = ein0
               ecc = ecc0

*-----------------------------------------------------------------------

               call tregck(iblz1,ilev1,ilat1,mr,kr,jj,icc)

               if( icc .eq. 0 ) goto 100

               ir = ii

*-----------------------------------------------------------------------
*        check of particles, only charged particles
*-----------------------------------------------------------------------
            if( iegsemi .eq. 0 ) then

               if( jtyp .eq. 0 ) goto 200

            else

               if( jtyp .eq. 0 .and. ityp .ne. 14 ) goto 200

            end if

            call pcheck(m,np,ityp,ktyp,jtyp,ipn,ips)

               if( ipn .eq. 0 )  goto 200

*-----------------------------------------------------------------------
*        check of counter
*-----------------------------------------------------------------------

            do i = 1, 3

               if( itcnt(i,m) .ne. 0 ) then

                  if( ncnt(ibknct+i,no,ipomp+1) .lt. itcnt(i*2+2,m) .or.
     &                ncnt(ibknct+i,no,ipomp+1) .gt. itcnt(i*2+3,m) )
     &                 goto 200

               end if

            end do

*-----------------------------------------------------------------------
*        track structure mode
*-----------------------------------------------------------------------
                mat1 = 0
         if(  mntsc .gt. 0 .and.
     &           ntscell( idgr(iblz(ibkblz+no,ipomp+1)) ) .ne. 0 .and.

     &        (( ityp .eq. 12 .or. ityp .eq. 13 ) .and.
     &         ( e(ibke+no,ipomp+1) .le. etsmax .and.
     &           e(ibke+no,ipomp+1) .ge. 0.d0 )   .or.

     &         ( ityp .eq. 1 ) .and.
     &         ( e(ibke+no,ipomp+1) .le. ptsmax .and.
     &           e(ibke+no,ipomp+1) .ge. 0.d0 )    .or.

     &         ( ktyp .eq. 6000012 ) .and.
     &         ( e(ibke+no,ipomp+1) .le. ctsmax .and.
     &           e(ibke+no,ipomp+1) .ge. 0.d0 )    .or.

     &         ( ibryf(ityp,ktyp) .ge. 1 ) .and.
     &         ( e(ibke+no,ipomp+1) .ge. 0.d0 .and.
     &           e(ibke+no,ipomp+1) .le. tsmax*dble(ibryf(ityp,ktyp)) ))

     &     ) then

               mat1 = ntscell( idgr(iblz(ibkblz+no,ipomp+1)) )

               if(mat1 .ne.  0) then
                   if( ityp .eq. 12 .or. ityp .eq. 13 ) ttsmin = etsmin
                   if( ibryf(ityp,ktyp) .ge. 1)
     &                ttsmin = ibryf(ityp,ktyp) * emin(ityp)
               endif
               if( mat1 .eq.  1 ) then
                   if( ityp .eq.  1 )                   ttsmin = ptsmin
                   if( ktyp .eq. 6000012)               ttsmin = ctsmin
               endif

               if(ncol .eq. 10 .or. ncol .eq. 11)then
                  tlv = oldwt*(e(ibke+no,ipomp+1)-ec(ibkec+no,ipomp+1))
               elseif(ncol .eq. 13) then
                  tlv = oldwt * dexc_ene
               elseif(e(ibke+no,ipomp+1).le.ttsmin)then
                  tlv = oldwt*e(ibke+no,ipomp+1)
               elseif(ec(ibkec+no,ipomp+1).lt.ttsmin .and.
     &                     ibryf(ityp,ktyp) .ge. 1) then
                  tlv =      0.d0
               else
                  tlv = oldwt * dexc_ene
               endif

               if(itunt(m).eq.0.and.mat.gt.0)then
                  tlv = tlv / rhog(mat)
               endif

               do ip = 1, ipn

                  trEVENT(ips(ip),ir,it,il) =
     &            trEVENT(ips(ip),ir,it,il) + tlv

               end do

         if(sweight(m) .ne. 1.d30 .and. sweight(m) .ne. oldwt)
     &     write(*,9000)
         sweight(m) = oldwt ! T.Sato 2018/03/07

         if( ilatflag .eq. 0 .and. tlv .gt. 0.d0 ) then
            ilatflag = 1
            nldtct(m) = nldtct(m) + 1
         endif

         goto 200
         endif

*-----------------------------------------------------------------------
*        dedxfnc is specified
*-----------------------------------------------------------------------

         if( itdfn(m,1) .ne. 0 ) then

               if( ein .lt. eletb(1) ) goto 200
               if( ecc .ge. eletb(nlete+1) ) goto 200

*-----------------------------------------------------------------------

               rrl = rrl0
               rrr = rrr0
               rrt = rrt0

*-----------------------------------------------------------------------
*           tally
*-----------------------------------------------------------------------

            do 270 ie = ie1, ie2, -1

                     emie = max( ecc, eletb(ie) )
                     emiu = min( ein, eletb(ie+1) )

                  if( emiu * emie .gt. 0.d0 ) then
                     erg  = sqrt( emiu * emie )
                  else
                     erg  = emiu / 2.d0
                  end if

                  if( mat .gt. 0 ) then

                     call rainge(emie,rng,mat,ityp,ktyp,jtyp,rtyp)

                     rng = max( 0.0d0, rng - rng1 )
                     rrl = max( 0.0d0, rrr - rng )

                        rrt  = rrt + rrl

                     if( rrt .gt. tlngth ) then

                        rrl = rrl + ( tlngth - rrt )
                        rrt = tlngth

                     else if( ie .eq. ie2 .and.
     &                        ecc .gt. eletb(ie) ) then

                        rrl = rrl + ( tlngth - rrt )

                     end if

                  end if

*-----------------------------------------------------------------------
*              dedx value multiplyed user factor
*-----------------------------------------------------------------------

               call dedxas(erg,dedx,lmat,ityp,ktyp,jtyp,rtyp)

                  dedx = dedx / 10.d0
                  dtrk = rrl
                  eini = emiu
                  efin = emie

               if( itdfn(m,1) .eq. 1 ) then

                  call usrdfn1(ityp,ktyp,jtyp,rtyp,
     &                         dedx,dtrk,eini,efin,dhet)

               else if( itdfn(m,1) .eq. 2 ) then

                  call usrdfn2(ityp,ktyp,jtyp,rtyp,
     &                         dedx,dtrk,eini,efin,dhet)

               end if

*-----------------------------------------------------------------------

                     tlv = dhet  ! T.Sato 2018/01/20, delete oldwt

                  do ip = 1, ipn

                     trEVENT(ips(ip),ir,it,il) =
     &               trEVENT(ips(ip),ir,it,il) + tlv

                  end do

         if(sweight(m) .ne. 1.d30 .and. sweight(m) .ne. oldwt)
     &     write(*,9000)
         sweight(m) = oldwt ! always remember the last energy deposition event, T.Sato 2018/03/07

                  rrr = rrr - rrl

                  if( ilatflag .eq. 0 .and. tlv .gt. 0.d0 ) then
                     ilatflag = 1
                     nldtct(m) = nldtct(m) + 1
                  endif

  270       continue

*-----------------------------------------------------------------------
*        normal case
*-----------------------------------------------------------------------

         else

                     tlv = ( ein - ecc ) ! T.Sato 2018/01/20, delete oldwt

                  do ip = 1, ipn

                     trEVENT(ips(ip),ir,it,il) =
     &               trEVENT(ips(ip),ir,it,il) + tlv

                  end do

         if(sweight(m) .ne. 1.d30 .and. sweight(m) .ne. oldwt)
     &     write(*,9000)
         sweight(m) = oldwt ! always remember the last energy deposition event, T.Sato 2018/03/07

                  if( ilatflag .eq. 0 .and. tlv .gt. 0.d0 ) then
                     ilatflag = 1
                     nldtct(m) = nldtct(m) + 1
                  endif

         end if

*-----------------------------------------------------------------------

  200       continue

*-----------------------------------------------------------------------
*     dead particles of produced charged particle
*-----------------------------------------------------------------------

      if( ( ncol .eq. 13 .or. ncol .eq. 14 ) .and.
     &      mat .gt. 0 .and. nclsts .gt. 0 ) then

*-----------------------------------------------------------------------
*        check of counter
*-----------------------------------------------------------------------

            do i = 1, 3

               if( itcnt(i,m) .ne. 0 ) then

                  if( jcount(i,1) .lt. itcnt(i*2+2,m) .or.
     &                jcount(i,1) .gt. itcnt(i*2+3,m) ) goto 300

               end if

            end do

*-----------------------------------------------------------------------
*        dead particle and nucleus except for electron positron
*-----------------------------------------------------------------------

         do j = 1, nclsts

         if( jclusts(4,j) .lt. 0 .and. jclusts(5,j) .ne. 0 .and.
     &       qclusts(7,j) .gt. eletb(1) ) then


                  ipart = jclusts(3,j)
                  kpart = jclusts(7,j)
                  jpart = jclusts(5,j)
                  rpart = qclusts(5,j) * 1000.d0

*-----------------------------------------------------------------------

            call pcheck(m,np,ipart,kpart,jpart,ipn,ips)

               if( ipn .eq. 0 )  goto 285

*-----------------------------------------------------------------------

                  ein = qclusts(7,j)
                  ecc = 0.0d0

*-----------------------------------------------------------------------
*        dedxfnc is specified
*-----------------------------------------------------------------------

         if( itdfn(m,1) .ne. 0 ) then

               if( ein .lt. eletb(1) ) goto 285

               do i = 1, nlete

                  if( ein .ge. eletb(i) .and.
     &                ein .lt. eletb(i+1) ) goto 35

               end do

   35             ie1 = min( i, nlete )
                  ie2 = 1

*-----------------------------------------------------------------------
*           rrl : track length and rrr : charge particles correction
*-----------------------------------------------------------------------

                  call rainge(ein,rng0,mat,
     &                        ipart,kpart,jpart,rpart)

                  if( rng0 .gt. 1.d+30 ) goto 285

                  tlngth = rng0
                  rng1 = 0.0d0

                  rrl = tlngth
                  rrr = tlngth
                  rrt = 0.0d0

               if( ein .gt. eletb(nlete+1) ) then

                  call rainge(eletb(nlete+1),rngm,mat,
     &                        ipart,kpart,jpart,rpart)

                  rrr = max( 0.0d0, rrr - ( rng0 - rngm ) )
                  rrt = rrt + rng0 -rngm

               end if

*-----------------------------------------------------------------------
*           tally
*-----------------------------------------------------------------------

            do 275 ie = ie1, ie2, -1

                     emie = max( ecc, eletb(ie) )
                     emiu = min( ein, eletb(ie+1) )

                  if( emiu * emie .gt. 0.d0 ) then
                     erg  = sqrt( emiu * emie )
                  else
                     erg  = emiu / 2.d0
                  end if

                  if( mat .gt. 0 ) then

                     call rainge(emie,rng,mat,
     &                           ipart,kpart,jpart,rpart)

                     rng = max( 0.0d0, rng - rng1 )
                     rrl = max( 0.0d0, rrr - rng )

                        rrt  = rrt + rrl

                     if( rrt .gt. tlngth ) then

                        rrl = rrl + ( tlngth - rrt )
                        rrt = tlngth

                     else if( ie .eq. ie2 .and.
     &                        ecc .gt. eletb(ie) ) then

                        rrl = rrl + ( tlngth - rrt )

                     end if

                  end if

*-----------------------------------------------------------------------
*              dedx value multiplyed user factor
*-----------------------------------------------------------------------

               call dedxas(erg,dedx,lmat,ipart,kpart,jpart,rpart)

                  dedx = dedx / 10.d0
                  dtrk = rrl
                  eini = emiu
                  efin = emie

               if( itdfn(m,1) .eq. 1 ) then

                  call usrdfn1(ipart,kpart,jpart,rpart,
     &                         dedx,dtrk,eini,efin,dhet)

               else if( itdfn(m,1) .eq. 2 ) then

                  call usrdfn2(ipart,kpart,jpart,rpart,
     &                         dedx,dtrk,eini,efin,dhet)

               end if

*-----------------------------------------------------------------------
cKN 2019/05/27 ?????? should be one
                     tlv = dhet  ! T.Sato 2018/01/20, delete qclusts(8,j)

                  do ip = 1, ipn

                     trEVENT(ips(ip),ir,it,il) =
     &               trEVENT(ips(ip),ir,it,il) + tlv

                  end do

         if(sweight(m) .ne. 1.d30 .and. sweight(m) .ne. oldwt)
     &    write(*,9000)
         sweight(m) = oldwt ! always remember the last energy deposition event, T.Sato 2018/03/07

                  rrr = rrr - rrl

                  if( ilatflag .eq. 0 .and. tlv .gt. 0.d0 ) then
                     ilatflag = 1
                     nldtct(m) = nldtct(m) + 1
                  endif

  275       continue

*-----------------------------------------------------------------------
*        normal case
*-----------------------------------------------------------------------

         else

                     tlv = ( ein - ecc ) ! T.Sato 2018/01/20, delete oldwt

                  do ip = 1, ipn

                     trEVENT(ips(ip),ir,it,il) =
     &               trEVENT(ips(ip),ir,it,il) + tlv

                  end do

         if(sweight(m) .ne. 1.d30 .and. sweight(m) .ne. oldwt)
     &    write(*,9000)
         sweight(m) = oldwt ! always remember the last energy deposition event, T.Sato 2018/03/07

                  if( ilatflag .eq. 0 .and. tlv .gt. 0.d0 ) then
                     ilatflag = 1
                     nldtct(m) = nldtct(m) + 1
                  endif

         end if

*-----------------------------------------------------------------------

  285       continue

         end if
         end do

*-----------------------------------------------------------------------

  300    continue

      end if

*-----------------------------------------------------------------------

  100 continue

*-----------------------------------------------------------------------

      return

 9000 Format('[T-deposit] with output = deposit, particle weight is not
     &unique in one history. Result may be unreasonable.')

      end

************************************************************************
*                                                                      *
      subroutine tdepsttetEVENT(ncol,m,nl,lt,np,nr,mr,
     &     ne,nt,kr,eb,tb,tr,
     &     trEVENT,tr0)
*                                                                      *
*       Deposit tally in tetra mesh                                    *
*       output = deposit                                               *
*       Last Modified by T.Furuta on 2025/01/16                        *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*      ncol  ..... reaction type                                       *
*              0 : end of batch                                        *
*              4 : source particle                                     *
*                                                                      *
*              9 : termination by time cut-off                         *
*             10 : geometry boundary crossing                          *
*             11 : termination by energy cut-off                       *
*             12 : termination by escape or leakage                    *
*             13 : nuclear reaction (n,x)                              *
*             14 : nuclear reaction (n,n'x)                            *
*             15 : sequential transport only for tally                 *
*                                                                      *
************************************************************************
      use MMBANKMOD   !FURUTA
      use EVENTTALMOD !FURUTA
      use tdepwgtsum_global !cABE 2016/11/24
      use partmod, only: itmxpt ! frtati 2021/10/05
*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param00.inc'
      include 'param.inc'

*-----------------------------------------------------------------------

      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
!$OMP THREADPRIVATE(/icomon/)
      common /wtsave/ oldwt
!$OMP THREADPRIVATE(/wtsave/)
      common /tlgeom/ iblz1,iblz2
!$OMP THREADPRIVATE(/tlgeom/)
      common /tlglat/ ilev1,ilev2,ilat1(5,10),ilat2(5,10)
!$OMP THREADPRIVATE(/tlglat/)

*-----------------------------------------------------------------------

      common /clustt/ nclsts, iclusts(nnn)
!$OMP THREADPRIVATE(/clustt/)
      common /clustw/ jclusts(0:8,nnn),  qclusts(0:12,nnn)
!$OMP THREADPRIVATE(/clustw/)
      common /cntcls/ jcount(3,nnn)
!$OMP THREADPRIVATE(/cntcls/)

*-----------------------------------------------------------------------

      common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                itout(itlmax), ittwo(itlmax)
      common /tall22/ itmcn(itlmax), itmtn(itlmax), itmtt(itlmax)
      common /tall36/ itdpo(itlmax)   ! S.Abe 2015/12/03
      common /tall37/ itcnt(9,itlmax)

      common /tall50/ itlmt(itlmax), ite2l(itlmax)   ! CCSE 2022.08.31
      common /tall51/ eletb(2000), nlete
      common /tall53/ itdfn(itlmax,2)

c T.Sato 2014/8/19 for detector resolution
      common /tall61/ rtdre(itlmax),rtdfa(itlmax)

      common /kmat1d/ idmn(0:kvlmax), idnm(kvmmax)

cABE 2016/03/04 for identification part=all
      common /tall65/ ipall(itlmax)

      common /tall70/ itrwgtsum(itlmax), itrgn1(itlmax),
     &                itrgm1(itlmax), itrncd(itlmax)

      common /tall00/ itnm, ital(itlmax), itals(itlmax), italm(itlmax)
      common /tall75/ itnlatcel(itlmax), itnlatmem(itlmax)

      common /tall82/ itcnth(9,itlmax)

      integer,save:: inilat
!$OMP THREADPRIVATE( inilat )
      data inilat /0/
      integer,allocatable,save:: nldtct(:)
!$OMP THREADPRIVATE( nldtct )
      integer,allocatable,save:: nlats(:,:)
!$OMP THREADPRIVATE( nlats )
      integer,allocatable,save:: nlatt(:,:)
!$OMP THREADPRIVATE( nlatt )
      integer,allocatable,save:: nlatu(:,:)
!$OMP THREADPRIVATE( nlatu )
      integer,save:: inilaterr
      data inilaterr /0/

*-----------------------------------------------------------------------

      dimension   kr(mr)
      dimension   lt(nl)
      dimension   eb(ne+1)
      dimension   tb(nt+1)
      dimension   tr(np,0:ne,nr,nt,2)
      dimension   trEVENT(np,nr,nt)

      dimension ips(itmxpt) ! frtati 2021/10/05 6 -> itmxpt

*-----------------------------------------------------------------------

      double precision, save :: sweight(itlmax)
!$OMP THREADPRIVATE(sweight)

      common /jcomon/ nabov,nobch,nocas,nomax
!$OMP THREADPRIVATE(/jcomon/)

      common /regdc/  idrg(kvlmax), idgr(kvmmax)
      common /  tscmsg   / ktsc(kvlmax), mntsc, ntsc(kvlmax)
      common /  tscreg   / ntscell(kvlmax)
      common / etsminmax / etsmin, etsmax
      common / ptsminmax / ptsmin, ptsmax
      common / ctsminmax / ctsmin, ctsmax
      common / tsminmax  / tsmax
      common / etsexe  / dexc_ene
!$OMP THREADPRIVATE(/etsexe/)
      common /celdg/  rhog(kvlmax)

      common /eparm/  esmax, esmin, emin(20)

*-----------------------------------------------------------------------

      real*8 edep
      common/EPCONTedep/ edep
!$OMP THREADPRIVATE(/EPCONTedep/)

      integer iegsemi, iegsout
      common /egsemi/ iegsemi, iegsout

C for nonshared_tally option
      integer italsh
      common /talsh/ italsh

      common /dedxfac/ dedxfd
!$OMP THREADPRIVATE(/dedxfac/)

      integer iii0,kkk0
      common /itettal2/ iii0,kkk0
!$OMP THREADPRIVATE(/itettal2/)

      common /cntmx/ ncntmx(3), nctmxsr(3)
!$OMP THREADPRIVATE(/cntmx/)

! sumover
      real(8) tr0(np,0:ne,nr,nt)

! sumover
      tr0(:,:,:,:) = 0.0d0

*-----------------------------------------------------------------------

      rng1=0.0 ! T.Sato 2024/04/02

! T.Sato 2014/8/21, for considering detector resolution, ncol = 101
      if(ncol.eq.101) then
       if(rtdre(m).eq.0.0.and.rtdfa(m).eq.0.0) return ! detector resolution = 0
       if( .not. FIRSTsrc ) then
        do ir = 1, nr
         do it = 1, nt
          gaurntmp=gaurn(dummy)  ! T.Sato 2018/01/25
          do ip = 1, np
           if( trEVENT(ip,ir,it).gt.0.0d0 ) then

            heats = trEVENT(ip,ir,it)

            if(rtdre(m) .lt. 0.d0) then ! Ogawa 2018/11/05 if sigma<0, gauss + exp asymmetric resolution
             heats = usrdefres(heats)
            else
             heats = heats + gaurntmp*  ! T.Sato 2018/01/25
     &       sqrt(rtdre(m)**2+rtdfa(m)*heats)
            endif

            if(heats.lt.0.0) heats = 0.0
            trEVENT(ip,ir,it)=heats
           endif
          enddo
         enddo
        enddo
       endif
       return
      endif
*----- End of revision on 2014/8/21 ------------------------------------

*-----------------------------------------------------------------------
* check of history counter
*-----------------------------------------------------------------------
         if ( ncol.eq.0 .or. ncol.eq.4 ) then
           do i = 1, 3
             if( itcnth(i,m) .eq. 1 ) then
               if( ncntmx(i) .lt. itcnth(i*2+2,m) .or.
     &             ncntmx(i) .gt. itcnth(i*2+3,m) ) then
                 trEVENT = 0.d0
               end if
             end if
           end do
         end if

*-----------------------------------------------------------------------
*     deposit : ncol = 0, 4 and ns = 4
*-----------------------------------------------------------------------

      if( ( ncol .eq. 0 .or. ncol .eq. 4 ) .and.
     &      itout(m) .ge. 2 ) then

         if( .not. FIRSTsrc ) then

            if( itrwgtsum(m) .eq. 1 ) then
             nr0 = 1
            else
             nr0 = nr
            endif

            do ip = 1, np
            do ir = 1, nr0
            do it = 1, nt

               if( itdpo(m) .eq. 0 .and.
     &             trEVENT(ipall(m),ir,it) .gt. 0.0 ) then
                  heats = trEVENT(ipall(m),ir,it)
                  ratio = trEVENT(ip,ir,it)
     &                   / trEVENT(ipall(m),ir,it)
               else
                  heats = trEVENT(ip,ir,it)
                  ratio = 1.0d0
               end if

               if( heats .gt. 0.0d0 .and.
     &             heats .ge. eb(1) .and. heats .lt. eb(ne+1) ) then

                  do ie = 1, ne

                     if( heats .ge. eb(ie) .and.
     &                   heats .lt. eb(ie+1) ) then
!$OMP CRITICAL (tr_tdepsttetEVENT_1)
! sumover
                        tr0(ip,ie,ir,it) =
     &                           tr0(ip,ie,ir,it) + sweight(m)*ratio
                        tr(ip,ie,ir,it,1) = tr(ip,ie,ir,it,1)
     &                                    + sweight(m)*ratio   ! S.Abe 2015/12/03
!$OMP END CRITICAL (tr_tdepsttetEVENT_1)
                        if(ip.eq.ipall(m).and.ratio.gt.0) then ! T.Sato 2022/09/09, part=all is always probability mode
!$OMP CRITICAL (tr_tdepsttetEVENT_2)
                         tr(ip,ie,ir,it,2) = tr(ip,ie,ir,it,2)
     &                                    + 1.d0            ! Ogawa 2020/04/22  if (weight < 1.d0), tr(ip,ie,ir,it,2) should be count number
!$OMP END CRITICAL (tr_tdepsttetEVENT_2)
                        else
!$OMP CRITICAL (tr_tdepsttetEVENT_3)
                         tr(ip,ie,ir,it,2) = tr(ip,ie,ir,it,2)
     &                                    + (sweight(m)*ratio)**2   ! S.Abe 2015/12/03
!$OMP END CRITICAL (tr_tdepsttetEVENT_3)
                        endif

                     end if

                  end do

               end if

            end do
            end do
            end do

! sumover
            call tdepstreg_sumover_ip(m,1, ipall(m), ratio,
     &                                np,  ne, nr,  nt, tr0)
           tr0(:,:,:,:) = 0.0d0

            do ip = 1, np
            do ir = 1, nr
            do it = 1, nt
               trEVENT(ip,ir,it) = 0.0d0
            end do
            end do
            end do

         end if

*-----------------------------------------------------------------------
*        save initial weight and zero set
*-----------------------------------------------------------------------

         if( ncol .eq. 4 ) then

            sweight(m) = 1.d30 ! initialize

         end if

      end if

*-----------------------------------------------------------------------
*        check of ncol
*-----------------------------------------------------------------------

         if( ncol .lt. 9 ) return

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
*        check region
*-----------------------------------------------------------------------

               jj = 0

               ic=idgr(iblz1)
               itet=kkk0-10000
               ihelem=iii0

               call ttetck(ic,mr,kr,itet,ihelem,ir,icc)

               if( icc .ne. 1 )return

  501      continue

*-----------------------------------------------------------------------
*        check time
*-----------------------------------------------------------------------

               tin = abs(t(ibkt+no,ipomp+1))

               if( tin .lt. tb(1) ) return
               if( tin .ge. tb(nt+1) ) return

            do i = 1, nt

               if( tin .ge. tb(i) .and.
     &             tin .lt. tb(i+1) ) goto 50

            end do

   50          it = i

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
*           track length
*-----------------------------------------------------------------------

               tlngth = sqrt( ( xc(ibkxc+no,ipomp+1) -
     &                                        x(ibkx+no,ipomp+1) )**2
     &                      + ( yc(ibkyc+no,ipomp+1) -
     &                                        y(ibky+no,ipomp+1) )**2
     &                      + ( zc(ibkzc+no,ipomp+1) -
     &                                        z(ibkz+no,ipomp+1) )**2 )

*-----------------------------------------------------------------------
*           initial and final energy
*-----------------------------------------------------------------------

                  ein = e(ibke+no,ipomp+1)
                  ecc = ec(ibkec+no,ipomp+1)

                  ecc = ein - dedxfd * ( ein - ecc )   ! S.Abe 2015/09/08

                  if( ncol .eq. 9 .or. ncol .eq. 11 ) ecc = 0.d0

*-----------------------------------------------------------------------
cKN Iwase 2014/08/22 for EGS

               if( ( ityp.eq.12 .or. ityp.eq.13 .or. ityp.eq.14 ) .and.
     &               iegsemi .ne. 0 ) then

                     call egs5edxde(ecc,tlngth,ein,mat,ityp)

                  if(ncol.eq.9 .or. ncol.eq.11) then
                     ecc = 0d0
                     if( ityp.eq.12 .or. ityp.eq.13 ) ein = ein + edep
                  endif

               endif

*-----------------------------------------------------------------------
*        dedxfnc is specified
*-----------------------------------------------------------------------

         if( itdfn(m,1) .ne. 0 ) then

               do i = 1, nlete

                  if( ein .ge. eletb(i) .and.
     &                ein .lt. eletb(i+1) ) goto 30

               end do

   30             ie1 = min( i, nlete )

               do i = nlete, 1, -1

                  if( ecc .ge. eletb(i) .and.
     &                ecc .lt. eletb(i+1) ) goto 40

               end do

   40             ie2 = max( i, 1 )

*-----------------------------------------------------------------------
*           rrl : track length and rrr : charge particles correction
*-----------------------------------------------------------------------

               rrl = tlngth
               rrr = tlngth
               rrt = 0.0d0

            if( mat .gt. 0 .and. jtyp .ne. 0 ) then

                  call rainge(ein,rng0,mat,ityp,ktyp,jtyp,rtyp)

                  rng1 = max( 0.0d0, rng0 - tlngth )

                  if( rng1 .gt. 1.d+30 ) return

               if( ein .gt. eletb(nlete+1) ) then

                  call rainge(eletb(nlete+1),rngm,mat,
     &                        ityp,ktyp,jtyp,rtyp)

                  rrr = max( 0.0d0, rrr - ( rng0 - rngm ) )
                  rrt = rrt + rng0 -rngm

               end if

            end if

               rrl0 = rrl
               rrr0 = rrr
               rrt0 = rrt

         end if

*-----------------------------------------------------------------------
*        check of particles, only charged particles
*-----------------------------------------------------------------------
            if( iegsemi .eq. 0 ) then

               if( jtyp .eq. 0 ) goto 200

            else

               if( jtyp .eq. 0 .and. ityp .ne. 14 ) goto 200

            end if

            call pcheck(m,np,ityp,ktyp,jtyp,ipn,ips)

               if( ipn .eq. 0 )  goto 200

*-----------------------------------------------------------------------
*        check of counter
*-----------------------------------------------------------------------

            do i = 1, 3

               if( itcnt(i,m) .ne. 0 ) then

                  if( ncnt(ibknct+i,no,ipomp+1) .lt. itcnt(i*2+2,m) .or.
     &                ncnt(ibknct+i,no,ipomp+1) .gt. itcnt(i*2+3,m) )
     &                           goto 200

               end if

            end do

*-----------------------------------------------------------------------
*        track structure mode
*-----------------------------------------------------------------------
                mat1 = 0
            if(  mntsc .gt. 0 .and.
     &           ntscell( idgr(iblz(ibkblz+no,ipomp+1)) ) .ne. 0 .and.

     &        (( ityp .eq. 12 .or. ityp .eq. 13 ) .and.
     &         ( e(ibke+no,ipomp+1) .le. etsmax .and.
     &           e(ibke+no,ipomp+1) .ge. 0.d0 )   .or.

     &         ( ityp .eq. 1 ) .and.
     &         ( e(ibke+no,ipomp+1) .le. ptsmax .and.
     &           e(ibke+no,ipomp+1) .ge. 0.d0 )    .or.

     &         ( ktyp .eq. 6000012 ) .and.
     &         ( e(ibke+no,ipomp+1) .le. ctsmax .and.
     &           e(ibke+no,ipomp+1) .ge. 0.d0 )    .or.

     &         ( ibryf(ityp,ktyp) .ge. 1 ) .and.
     &         ( e(ibke+no,ipomp+1) .ge. 0.d0 .and.
     &           e(ibke+no,ipomp+1) .le. tsmax*dble(ibryf(ityp,ktyp)) ))

     &        ) mat1 = ntscell( idgr(iblz(ibkblz+no,ipomp+1)) )


         if(
     .         ( ( ktyp.eq.6000012            ) .and.
     .             mat1 .eq. 1                  .and.
     .             ( ec(ibkec+no,ipomp+1) .le. ctsmax )  .and.
     .             ( ec(ibkec+no,ipomp+1) .ge. 0.d0 ) )
     .   .or.
     .         ( ( ityp.eq.1                  ) .and.
     .             mat1 .eq. 1                  .and.
     .             ( ec(ibkec+no,ipomp+1) .le. ptsmax )  .and.
     .             ( ec(ibkec+no,ipomp+1) .ge. 0.d0 ) )
     .   .or.
     .         ( ( ityp.eq.12 .or. ityp.eq.13 ) .and.
     .             mat1 .ne. 0                  .and.
     .             ( ec(ibkec+no,ipomp+1) .le. etsmax )  .and.
     .             ( ec(ibkec+no,ipomp+1) .ge. 0.d0 ) )
     .   .or.
     .         ( (  ibryf(ityp,ktyp) .ge. 1   ) .and.
     .             mat1 .ne. 0                  .and.
     .             ( ec(ibkec+no,ipomp+1) .le. tsmax
     .                       *dble(ibryf(ityp,ktyp)) )   .and.
     .             ( ec(ibkec+no,ipomp+1) .ge. 0.d0 )  )
     .   ) then

               if(mat1 .ne.  0) then
                   if( ityp .eq. 12 .or. ityp .eq. 13 ) ttsmin = etsmin
                   if( ibryf(ityp,ktyp) .ge. 1)
     &                ttsmin = ibryf(ityp,ktyp) * emin(ityp)
               endif
               if( mat1 .eq.  1 ) then
                   if( ityp .eq.  1 )                   ttsmin = ptsmin
                   if( ktyp .eq. 6000012)               ttsmin = ctsmin
               endif

               if(ncol .eq. 10 .or. ncol .eq. 11)then
                  tlv=oldwt*(e(ibke+no,ipomp+1)-ec(ibkec+no,ipomp+1))
               elseif(ncol .eq. 13) then
                  tlv = oldwt * dexc_ene
               elseif(e(ibke+no,ipomp+1).le.ttsmin)then
                  tlv = oldwt * e(ibke+no,ipomp+1)
               elseif(ec(ibkec+no,ipomp+1).lt.ttsmin .and.
     &                     ibryf(ityp,ktyp) .ge. 1) then
                  tlv = 0.d0
               else
                  tlv = oldwt * dexc_ene
               endif

               if(itunt(m).eq.0.and.mat.gt.0)then
                  tlv = tlv / rhog(mat)
               endif

               do ip = 1, ipn

                  trEVENT(ips(ip),ir,it) =
     &            trEVENT(ips(ip),ir,it) + tlv

               end do

         if(sweight(m) .ne. 1.d30 .and. sweight(m) .ne. oldwt)
     &      write(*,9000)
         sweight(m) = oldwt ! T.Sato 2018/03/07  remember weight of this history
         goto 100
         endif

*-----------------------------------------------------------------------
*        dedxfnc is specified
*-----------------------------------------------------------------------

         if( itdfn(m,1) .ne. 0 ) then

               if( ein .lt. eletb(1) ) goto 200
               if( ecc .ge. eletb(nlete+1) ) goto 200

*-----------------------------------------------------------------------

               rrl = rrl0
               rrr = rrr0
               rrt = rrt0

*-----------------------------------------------------------------------
*           tally
*-----------------------------------------------------------------------

            do 270 ie = ie1, ie2, -1

                     emie = max( ecc, eletb(ie) )
                     emiu = min( ein, eletb(ie+1) )

                  if( emiu * emie .gt. 0.d0 ) then
                     erg  = sqrt( emiu * emie )
                  else
                     erg  = emiu / 2.d0
                  end if

                  if( mat .gt. 0 ) then

                     call rainge(emie,rng,mat,ityp,ktyp,jtyp,rtyp)

                     rng = max( 0.0d0, rng - rng1 )
                     rrl = max( 0.0d0, rrr - rng )

                        rrt  = rrt + rrl

                     if( rrt .gt. tlngth ) then

                        rrl = rrl + ( tlngth - rrt )
                        rrt = tlngth

                     else if( ie .eq. ie2 .and.
     &                        ecc .gt. eletb(ie) ) then

                        rrl = rrl + ( tlngth - rrt )

                     end if

                  end if

*-----------------------------------------------------------------------
*              dedx value multiplyed user factor
*-----------------------------------------------------------------------

               call dedxas(erg,dedx,lmat,ityp,ktyp,jtyp,rtyp)

                  dedx = dedx / 10.d0
                  dtrk = rrl
                  eini = emiu
                  efin = emie

               if( itdfn(m,1) .eq. 1 ) then

                  call usrdfn1(ityp,ktyp,jtyp,rtyp,
     &                         dedx,dtrk,eini,efin,dhet)

               else if( itdfn(m,1) .eq. 2 ) then

                  call usrdfn2(ityp,ktyp,jtyp,rtyp,
     &                         dedx,dtrk,eini,efin,dhet)

               end if

*-----------------------------------------------------------------------

                     tlv = dhet  ! T.Sato 2018/01/20, delete oldwt

                  do ip = 1, ipn

                     trEVENT(ips(ip),ir,it) =
     &               trEVENT(ips(ip),ir,it) + tlv

                  end do

                  if(sweight(m) .ne. 1.d30 .and. sweight(m) .ne. oldwt)
     &              write(*,9000)
                  sweight(m) = oldwt ! T.Sato 2018/03/07  remember weight of this history

                  rrr = rrr - rrl

  270       continue

*-----------------------------------------------------------------------
*        normal case
*-----------------------------------------------------------------------

         else

                     tlv = ( ein - ecc ) ! T.Sato 2018/01/20, delete oldwt

                  do ip = 1, ipn

                     trEVENT(ips(ip),ir,it) =
     &               trEVENT(ips(ip),ir,it) + tlv

                  end do

                  if(sweight(m) .ne. 1.d30 .and. sweight(m) .ne. oldwt)
     &              write(*,9000)
                  sweight(m) = oldwt ! T.Sato 2018/03/07  remember weight of this history

         end if

*-----------------------------------------------------------------------

  200       continue

*-----------------------------------------------------------------------
*     dead particles of produced charged particle
*-----------------------------------------------------------------------

      if( ( ncol .eq. 13 .or. ncol .eq. 14 ) .and.
     &      mat .gt. 0 .and. nclsts .gt. 0 ) then

*-----------------------------------------------------------------------
*        check of counter
*-----------------------------------------------------------------------

            do i = 1, 3

               if( itcnt(i,m) .ne. 0 ) then

                  if( jcount(i,1) .lt. itcnt(i*2+2,m) .or.
     &                jcount(i,1) .gt. itcnt(i*2+3,m) ) goto 300

               end if

            end do

*-----------------------------------------------------------------------
*        dead particle and nucleus except for electron positron
*-----------------------------------------------------------------------

         do j = 1, nclsts

         if( jclusts(4,j) .lt. 0 .and. jclusts(5,j) .ne. 0 .and.
     &       qclusts(7,j) .gt. eletb(1) ) then


                  ipart = jclusts(3,j)
                  kpart = jclusts(7,j)
                  jpart = jclusts(5,j)
                  rpart = qclusts(5,j) * 1000.d0

*-----------------------------------------------------------------------

            call pcheck(m,np,ipart,kpart,jpart,ipn,ips)

               if( ipn .eq. 0 )  goto 285

*-----------------------------------------------------------------------

                  ein = qclusts(7,j)
                  ecc = 0.0d0

*-----------------------------------------------------------------------
*        dedxfnc is specified
*-----------------------------------------------------------------------

         if( itdfn(m,1) .ne. 0 ) then

               if( ein .lt. eletb(1) ) goto 285

               do i = 1, nlete

                  if( ein .ge. eletb(i) .and.
     &                ein .lt. eletb(i+1) ) goto 35

               end do

   35             ie1 = min( i, nlete )
                  ie2 = 1

*-----------------------------------------------------------------------
*           rrl : track length and rrr : charge particles correction
*-----------------------------------------------------------------------

                  call rainge(ein,rng0,mat,
     &                        ipart,kpart,jpart,rpart)

                  if( rng0 .gt. 1.d+30 ) goto 285

                  tlngth = rng0
                  rng1 = 0.0d0

                  rrl = tlngth
                  rrr = tlngth
                  rrt = 0.0d0

               if( ein .gt. eletb(nlete+1) ) then

                  call rainge(eletb(nlete+1),rngm,mat,
     &                        ipart,kpart,jpart,rpart)

                  rrr = max( 0.0d0, rrr - ( rng0 - rngm ) )
                  rrt = rrt + rng0 -rngm

               end if

*-----------------------------------------------------------------------
*           tally
*-----------------------------------------------------------------------

            do 275 ie = ie1, ie2, -1

                     emie = max( ecc, eletb(ie) )
                     emiu = min( ein, eletb(ie+1) )

                  if( emiu * emie .gt. 0.d0 ) then
                     erg  = sqrt( emiu * emie )
                  else
                     erg  = emiu / 2.d0
                  end if

                  if( mat .gt. 0 ) then

                     call rainge(emie,rng,mat,
     &                           ipart,kpart,jpart,rpart)

                     rng = max( 0.0d0, rng - rng1 )
                     rrl = max( 0.0d0, rrr - rng )

                        rrt  = rrt + rrl

                     if( rrt .gt. tlngth ) then

                        rrl = rrl + ( tlngth - rrt )
                        rrt = tlngth

                     else if( ie .eq. ie2 .and.
     &                        ecc .gt. eletb(ie) ) then

                        rrl = rrl + ( tlngth - rrt )

                     end if

                  end if

*-----------------------------------------------------------------------
*              dedx value multiplyed user factor
*-----------------------------------------------------------------------

               call dedxas(erg,dedx,lmat,ipart,kpart,jpart,rpart)

                  dedx = dedx / 10.d0
                  dtrk = rrl
                  eini = emiu
                  efin = emie

               if( itdfn(m,1) .eq. 1 ) then

                  call usrdfn1(ipart,kpart,jpart,rpart,
     &                         dedx,dtrk,eini,efin,dhet)

               else if( itdfn(m,1) .eq. 2 ) then

                  call usrdfn2(ipart,kpart,jpart,rpart,
     &                         dedx,dtrk,eini,efin,dhet)

               end if

*-----------------------------------------------------------------------

                     tlv = dhet  ! T.Sato 2018/01/20, delete qclusts(8,j)

                  do ip = 1, ipn

                     trEVENT(ips(ip),ir,it) =
     &               trEVENT(ips(ip),ir,it) + tlv

                  end do

                  if(sweight(m) .ne. 1.d30 .and. sweight(m) .ne. oldwt)
     &              write(*,9000)
                  sweight(m) = oldwt ! T.Sato 2018/03/07  remember weight of this history

                  rrr = rrr - rrl

                  if( ilatflag .eq. 0 .and. tlv .gt. 0.d0 ) then
                     ilatflag = 1
                     nldtct(m) = nldtct(m) + 1
                  endif

  275       continue

*-----------------------------------------------------------------------
*        normal case
*-----------------------------------------------------------------------

         else

                     tlv = ( ein - ecc ) ! T.Sato 2018/01/20, delete oldwt

                  do ip = 1, ipn

                     trEVENT(ips(ip),ir,it) =
     &               trEVENT(ips(ip),ir,it) + tlv

                  end do

                  if(sweight(m) .ne. 1.d30 .and. sweight(m) .ne. oldwt)
     &              write(*,9000)
                  sweight(m) = oldwt ! T.Sato 2018/03/07  remember weight of this history

         end if

*-----------------------------------------------------------------------

  285       continue

         end if
         end do

*-----------------------------------------------------------------------

  300    continue

      end if

*-----------------------------------------------------------------------

  100 continue

*-----------------------------------------------------------------------

      return

 9000 Format('[T-deposit] with output = deposit, particle weight is not
     &unique in one history. Result may be unreasonable.')

      end

************************************************************************
*                                                                      *
      subroutine tdepstrzEVENT(ncol,m,nl,lt,np,nr,nz,ne,nt,rm,zm,eb,tb,
     &     tr,trEVENT,tr0)
*                                                                      *
*       Deposit tally in r-z scoring mesh                              *
*       output = deposit                                               *
*       last modified by T.Furuta on 2012/05/08                        *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*      ncol  ..... reaction type                                       *
*              0 : end of batch                                        *
*              4 : source particle                                     *
*                                                                      *
*              9 : termination by time cut-off                         *
*             10 : geometry boundary crossing                          *
*             11 : termination by energy cut-off                       *
*             12 : termination by escape or leakage                    *
*             13 : nuclear reaction (n,x)                              *
*             14 : nuclear reaction (n,n'x)                            *
*             15 : sequential transport only for tally                 *
*                                                                      *
************************************************************************
      use MMBANKMOD   !FURUTA
      use EVENTTALMOD !FURUTA
      use partmod, only: itmxpt ! frtati 2021/10/05
*-----------------------------------------------------------------------

      implicit double precision( a-h, o-z )

*-----------------------------------------------------------------------

      include 'param00.inc'
      include 'param.inc'

*-----------------------------------------------------------------------

      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
!$OMP THREADPRIVATE(/icomon/)
      common /wtsave/ oldwt
!$OMP THREADPRIVATE(/wtsave/)
      common /tlgeom/ iblz1,iblz2
!$OMP THREADPRIVATE(/tlgeom/)

*-----------------------------------------------------------------------

      common /clustt/ nclsts, iclusts(nnn)
!$OMP THREADPRIVATE(/clustt/)
      common /clustw/ jclusts(0:8,nnn),  qclusts(0:12,nnn)
!$OMP THREADPRIVATE(/clustw/)
      common /cntcls/ jcount(3,nnn)
!$OMP THREADPRIVATE(/cntcls/)

*-----------------------------------------------------------------------

      common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                itout(itlmax), ittwo(itlmax)
      common /tall08/ rtrx0(itlmax), rtry0(itlmax)
      common /tall22/ itmcn(itlmax), itmtn(itlmax), itmtt(itlmax)
      common /tall36/ itdpo(itlmax)   ! S.Abe 2015/09/10
      common /tall37/ itcnt(9,itlmax)
      common /tall46/ itmtr(itlmax,4), rtmtr(itlmax,13)

      common /tall50/ itlmt(itlmax), ite2l(itlmax)   ! CCSE 2022.08.31
      common /tall51/ eletb(2000), nlete
      common /tall53/ itdfn(itlmax,2)

c T.Sato 2014/8/19 for detector resolution
      common /tall61/ rtdre(itlmax),rtdfa(itlmax)

      common /kmat1d/ idmn(0:kvlmax), idnm(kvmmax)

cABE 2016/03/04 for identification part=all
      common /tall65/ ipall(itlmax)

      common /tall82/ itcnth(9,itlmax)

*-----------------------------------------------------------------------

      dimension   lt(nl)
      dimension rm(nr+1)
      dimension zm(nz+1)
      dimension eb(ne+1)
      dimension tb(nt+1)
      dimension tr(np,0:ne,nt,nr,nz,2)
      dimension trEVENT(np,nt,nr,nz)

*-----------------------------------------------------------------------

      dimension ud(3)

      data dmax  /1.0d+19/
      data dmax0 /1.0d+18/

      dimension ips(itmxpt) ! frtati 2021/10/05 6 -> itmxpt
      dimension facm(6)

*-----------------------------------------------------------------------

      common /paraj/ mstz(300), parz(300)

*-----------------------------------------------------------------------

      double precision, save :: sweight(itlmax)
!$OMP THREADPRIVATE(sweight)

      common /jcomon/ nabov,nobch,nocas,nomax
!$OMP THREADPRIVATE(/jcomon/)

*-----------------------------------------------------------------------

      real*8 edep
      common/EPCONTedep/ edep
!$OMP THREADPRIVATE(/EPCONTedep/)

      integer iegsemi, iegsout
      common /egsemi/ iegsemi, iegsout

C for nonshared_tally option
      integer italsh
      common /talsh/ italsh

*------------------------------------------------------------------------

      common /spred/ nspred, nwsprd, nedisp, itstep, ndedx

      common /regdc/  idrg(kvlmax), idgr(kvmmax)
      common /  tscmsg   / ktsc(kvlmax), mntsc, ntsc(kvlmax)
      common /  tscreg   / ntscell(kvlmax)
      common / etsminmax / etsmin, etsmax
      common / ptsminmax / ptsmin, ptsmax
      common / ctsminmax / ctsmin, ctsmax
      common / tsminmax  / tsmax ! T.Ogawa 2020/06/15
      common / etsexe  / dexc_ene
!$OMP THREADPRIVATE(/etsexe/)
      common /celdg/  rhog(kvlmax)

      common /eparm/  esmax, esmin, emin(20)

      common /cntmx/ ncntmx(3), nctmxsr(3)
!$OMP THREADPRIVATE(/cntmx/)

! sumover
      real(8) tr0(np,0:ne,nt,nr,nz)

      tr0(:,:,:,:,:) = 0.0d0

*-----------------------------------------------------------------------

      rng1=0.0 ! T.Sato 2024/04/02

*-----------------------------------------------------------------------
! T.Sato 2014/8/21, for considering detector resolution, ncol = 101
      if(ncol.eq.101) then
       if(rtdre(m).eq.0.0.and.rtdfa(m).eq.0.0) return ! detector resolution = 0
       if( .not. FIRSTsrc ) then
        do it = 1, nt
         do ir = 1, nr
          do iz = 1, nz
          gaurntmp=gaurn(dummy)  ! T.Sato 2018/01/25
           do ip = 1, np
            if( trEVENT(ip,it,ir,iz).gt.0.0d0 ) then
             heats = trEVENT(ip,it,ir,iz)

             if(rtdre(m) .lt. 0.d0) then ! Ogawa 2018/11/05 if sigma<0, gauss + exp asymmetric resolution
              heats = usrdefres(heat)
             else
              heats = heats + gaurntmp*  ! T.Sato 2018/01/25
     &        sqrt(rtdre(m)**2+rtdfa(m)*heats)
             endif

             if(heats.lt.0.0) heats = 0.0
             trEVENT(ip,it,ir,iz)=heats
            endif
           enddo
          enddo
         enddo
        enddo
       endif
       return
      endif
*----- End of revision on 2014/8/21 ------------------------------------

            small = parz(28) * 1.d-4 ! T.Sato 2023/11/01 avoid infinite loop

*-----------------------------------------------------------------------
* check of history counter
*-----------------------------------------------------------------------
         if ( ncol.eq.0 .or. ncol.eq.4 ) then
           do i = 1, 3
             if( itcnth(i,m) .eq. 1 ) then
               if( ncntmx(i) .lt. itcnth(i*2+2,m) .or.
     &             ncntmx(i) .gt. itcnth(i*2+3,m) ) then
                 trEVENT = 0.d0
               end if
             end if
           end do
         end if

*-----------------------------------------------------------------------
*     deposit : ncol = 0, 4 and ns = 4
*-----------------------------------------------------------------------

      if( ( ncol .eq. 0 .or. ncol .eq. 4 ) .and.
     &      itout(m) .ge. 2 ) then

         if( .not. FIRSTsrc ) then

            do ip = 1, np
            do it = 1, nt
            do ir = 1, nr
            do iz = 1, nz

               if( itdpo(m) .eq. 0 .and.
     &             trEVENT(ipall(m),it,ir,iz) .gt. 0.0 ) then
                  heats = trEVENT(ipall(m),it,ir,iz)
                  ratio = trEVENT(ip,it,ir,iz)
     &                   / trEVENT(ipall(m),it,ir,iz)
               else
                  heats = trEVENT(ip,it,ir,iz)
                  ratio = 1.0d0
               end if

               if( heats .gt. 0.0d0 .and.
     &             heats .ge. eb(1) .and. heats .lt. eb(ne+1) ) then

                  do ie = 1, ne

                     if( heats .ge. eb(ie) .and.
     &                   heats .lt. eb(ie+1) ) then
!$OMP CRITICAL (tr_tdepstrzEVENT_1)
! sumover
                        tr0(ip,ie,it,ir,iz) =
     &                          tr0(ip,ie,it,ir,iz) + sweight(m)*ratio
                        tr(ip,ie,it,ir,iz,1) = tr(ip,ie,it,ir,iz,1)
     &                                    + sweight(m)*ratio   ! S.Abe 2015/12/03
!$OMP END CRITICAL (tr_tdepstrzEVENT_1)
                        if(ip.eq.ipall(m).and.ratio.gt.0) then ! T.Sato 2022/09/09, part=all is always probability mode
!$OMP CRITICAL (tr_tdepstrzEVENT_2)
                         tr(ip,ie,it,ir,iz,2) = tr(ip,ie,it,ir,iz,2)
     &                                    + 1.d0            ! Ogawa 2020/04/22  if (weight < 1.d0), tr(ip,ie,ir,it,2) should be count number
!$OMP END CRITICAL (tr_tdepstrzEVENT_2)
                        else
!$OMP CRITICAL (tr_tdepstrzEVENT_3)
                         tr(ip,ie,it,ir,iz,2) = tr(ip,ie,it,ir,iz,2)
     &                                    + (sweight(m)*ratio)**2   ! S.Abe 2015/12/03
!$OMP END CRITICAL (tr_tdepstrzEVENT_3)

                        endif

                     end if

                  end do

               end if

            end do
            end do
            end do
            end do

! sumover
            call tdepstrz_sumover_ip(m,1, ipall(m), ratio,
     &                               np,  ne, nt, nr,  nz, tr0)
           tr0(:,:,:,:,:) = 0.0d0

            do ip = 1, np
            do it = 1, nt
            do ir = 1, nr
            do iz = 1, nz
               trEVENT(ip,it,ir,iz) = 0.0d0
            end do
            end do
            end do
            end do

         end if

*-----------------------------------------------------------------------
*        save initial weight and zero set
*-----------------------------------------------------------------------

         if( ncol .eq. 4 ) then

            sweight(m) = 1.d30 ! initialize

         end if

      end if

*-----------------------------------------------------------------------
*        check of ncol
*-----------------------------------------------------------------------

         if( ncol .lt. 9 ) return

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
*        check time
*-----------------------------------------------------------------------

               tin = abs(t(ibkt+no,ipomp+1))

               if( tin .lt. tb(1) ) return
               if( tin .ge. tb(nt+1) ) return

            do i = 1, nt

               if( tin .ge. tb(i) .and.
     &             tin .lt. tb(i+1) ) goto 501

            end do

  501          it = i

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
*        transform positions
*-----------------------------------------------------------------------

            call trnsxx(x(ibkx+no,ipomp+1),y(ibky+no,ipomp+1),
     &                  z(ibkz+no,ipomp+1),
     &                  xxa,yya,zza,itmtr(m,4))

            call trnsxx(xc(ibkxc+no,ipomp+1),yc(ibkyc+no,ipomp+1),
     &                  zc(ibkzc+no,ipomp+1),
     &                  xxc,yyc,zzc,itmtr(m,4))

*-----------------------------------------------------------------------
* add dummy step for EGS5 electrons at ncol=11
*     (stopped but still have dE)
cc H.Iwase 2014/8/21 add the change of H.Iwase 2014/1/9
*-----------------------------------------------------------------------

         if( iegsemi .ne. 0 .and.
     &     ( ityp .eq. 12 .or. ityp .eq. 13 ) .and. ncol .eq. 11 ) then

               dum = ( xxc - xxa )**2
     &             + ( yyc - yya )**2
     &             + ( zzc - zza )**2

            if( dum .le. small**2 ) then

               dumstep = small * 10.d0

               xdum = x(ibkx+no,ipomp+1) + u(ibku+no,ipomp+1)*dumstep
               ydum = y(ibky+no,ipomp+1) + v(ibkv+no,ipomp+1)*dumstep
               zdum = z(ibkz+no,ipomp+1) + w(ibkw+no,ipomp+1)*dumstep

               call trnsxx(xdum,ydum,zdum,xxc,yyc,zzc,itmtr(m,4))

            end if

         end if
*-----------------------------------------------------------------------
*        check z mesh and r mesh region
*-----------------------------------------------------------------------

               x0 = rtrx0(m)
               y0 = rtry0(m)

            if( zza-small .lt. zm(1) .and.
     &          zzc-small .lt. zm(1) ) return

            if( zza+small .ge. zm(nz+1) .and.
     &          zzc+small .ge. zm(nz+1) ) return

               dis0 = sqrt( ( xxa - x0 )**2
     &                    + ( yya - y0 )**2 )

               dis1 = sqrt( ( xxc - x0 )**2
     &                    + ( yyc - y0 )**2 )

            if( dis0 .lt. rm(1) .and.
     &          dis1 .lt. rm(1) ) return

            if( dis0 .ge. rm(nr+1) .and.
     &          dis1 .ge. rm(nr+1) ) then

               aa =   yyc - yya
               bb = - xxc + xxa
               cc = - aa * xxa - bb * yya

               if( aa**2 + bb**2 .ne. 0.0d0 ) then

                  dd = abs( aa * x0 + bb * y0 + cc )
     &               / sqrt( aa**2 + bb**2 )

                  if( dd .ge. rm(nr+1) ) return

               end if

            end if

*-----------------------------------------------------------------------
*        check of particles, only charged particles
*-----------------------------------------------------------------------

            if( iegsemi .eq. 0 ) then

               if( jtyp .eq. 0 ) goto 200
               if( abs( ec(ibkec+no,ipomp+1) - e(ibke+no,ipomp+1) )
     &             .eq. 0.0d0 )
     &             goto 200

             else

               if( jtyp .eq. 0 .and. ityp .ne. 14 ) goto 200
               if( ( abs( ec(ibkec+no,ipomp+1) - e(ibke+no,ipomp+1) )
     $                                        .eq. 0.0d0 )
     $              .and.
     $          ( ityp .ne. 12 .and. ityp .ne. 13 .and. ityp .ne. 14))
     $            goto 200

             end if

            call pcheck(m,np,ityp,ktyp,jtyp,ipn,ips)

               if( ipn .eq. 0 ) goto 200

*-----------------------------------------------------------------------
*        check of counter
*-----------------------------------------------------------------------

            do i = 1, 3

               if( itcnt(i,m) .ne. 0 ) then

                  if( ncnt(ibknct+i,no,ipomp+1) .lt. itcnt(i*2+2,m) .or.
     &                ncnt(ibknct+i,no,ipomp+1) .gt. itcnt(i*2+3,m) )
     &                     goto 200

               end if

            end do

*-----------------------------------------------------------------------
*        distance and unit vector ud(i)
*-----------------------------------------------------------------------

            dis = ( xxc - xxa )**2
     &          + ( yyc - yya )**2
     &          + ( zzc - zza )**2

                mat1 = 0
            if(  mntsc .gt. 0 .and.
     &           ntscell( idgr(iblz(ibkblz+no,ipomp+1)) ) .ne. 0 .and.

     &        (( ityp .eq. 12 .or. ityp .eq. 13 ) .and.
     &         ( e(ibke+no,ipomp+1) .le. etsmax .and.
     &           e(ibke+no,ipomp+1) .ge. 0.d0 )   .or.

     &         ( ityp .eq. 1 ) .and.
     &         ( e(ibke+no,ipomp+1) .le. ptsmax .and.
     &           e(ibke+no,ipomp+1) .ge. 0.d0 )    .or.

     &         ( ktyp .eq. 6000012 ) .and.
     &         ( e(ibke+no,ipomp+1) .le. ctsmax .and.
     &           e(ibke+no,ipomp+1) .ge. 0.d0 )    .or.

     &         ( ibryf(ityp,ktyp) .ge. 1 ) .and.
     &         ( e(ibke+no,ipomp+1) .ge. 0.d0 .and.
     &           e(ibke+no,ipomp+1) .le. tsmax*dble(ibryf(ityp,ktyp)) ))

     &        ) mat1 = ntscell( idgr(iblz(ibkblz+no,ipomp+1)) )


         if(
     .         ( ( ktyp.eq.6000012            ) .and.
     .             mat1 .eq. 1                  .and.
     .             ( ec(ibkec+no,ipomp+1) .le. ctsmax )  .and.
     .             ( ec(ibkec+no,ipomp+1) .ge. 0.d0 ) )
     .   .or.
     .         ( ( ityp.eq.1                  ) .and.
     .             mat1 .eq. 1                  .and.
     .             ( ec(ibkec+no,ipomp+1) .le. ptsmax )  .and.
     .             ( ec(ibkec+no,ipomp+1) .ge. 0.d0 ) )
     .   .or.
     .         ( ( ityp.eq.12 .or. ityp.eq.13 ) .and.
     .             mat1 .ne. 0                  .and.
     .             ( ec(ibkec+no,ipomp+1) .le. etsmax )  .and.
     .             ( ec(ibkec+no,ipomp+1) .ge. 0.d0 ) )
     .   .or.
     .         ( (  ibryf(ityp,ktyp) .ge. 1   ) .and.
     .             mat1 .ne. 0                  .and.
     .             ( ec(ibkec+no,ipomp+1) .le. tsmax
     .                       *dble(ibryf(ityp,ktyp)) )   .and.
     .             ( ec(ibkec+no,ipomp+1) .ge. 0.d0 )  )
     .   ) then
           else
            if( dis .le. small**2 ) return ! org
           endif

            dis = sqrt(dis)

            ud(1) = ( xxc - xxa ) / dis
            ud(2) = ( yyc - yya ) / dis
            ud(3) = ( zzc - zza ) / dis

*-----------------------------------------------------------------------
*     initial position, energy and initial range rng
*-----------------------------------------------------------------------

            tot = 0.0d0
            rng1 = 0.0d0 !FURUTA

            se  = e(ibke+no,ipomp+1)
            xpp = xxa
            ypp = yya
            zpp = zza

            if(  mat .gt. 0 )
     &          call rainge(se,rng,mat,ityp,ktyp,jtyp,rtyp)

*-----------------------------------------------------------------------
*     initial z-position
*-----------------------------------------------------------------------

               izm = 0
               izc = 0

            if( ud(3) .gt. 0.0d0 ) izk =  1
            if( ud(3) .eq. 0.0d0 ) izk =  0
            if( ud(3) .lt. 0.0d0 ) izk = -1

            if( izk .ge. 0 ) then

               do i = 1, nz + 1
                  if( zm(i) .gt. zpp + small ) then
                         izm = i - 1
                         izc = i - 1
                         goto 38
                  end if
               end do

            else

               izm = nz + 2
               izc = nz + 1

               do i = 1, nz + 1
                  if( zm(i) .ge. zpp - small ) then
                         izm = i
                         izc = i - 1
                         goto 38
                  end if
               end do

            end if

   38       continue

*-----------------------------------------------------------------------
*     initial r-position
*-----------------------------------------------------------------------

               uvec = ( xxa - x0 ) * ud(1) + ( yya - y0 ) * ud(2)

            if( uvec .ge. 0.0d0 ) then

               do i = 1, nr + 1
                  if( rm(i) .gt. dis0 + small ) then
                        irc = i - 1
                        goto 47
                  end if
               end do

            else

               do i = 1, nr + 1
                  if( rm(i) .ge. dis0 - small ) then
                        irc = i - 1
                        goto 47
                  end if
               end do

            end if

                  irc = nr + 1

   47          continue

*-----------------------------------------------------------------------
*     loop for finding mesh and booking upto total distance
*-----------------------------------------------------------------------

   50 continue

*-----------------------------------------------------------------------
*      calculation is finished
*-----------------------------------------------------------------------

         if( tot .ge. dis ) goto 200

*-----------------------------------------------------------------------
*        dr : distance to the nearest r mesh
*-----------------------------------------------------------------------

               dr = dmax

               a  = ud(1)**2 + ud(2)**2

         if( a .gt. 0.0d0 ) then

               b = 2.0d0 * ( ud(1) * ( xpp - x0 )
     &                     + ud(2) * ( ypp - y0 ) )

               ir1 = irc
               ir2 = irc + 1

              if( b .ge. 0.0 ) ir1 = irc + 1

            do 51 i = max( 1, ir1 ), min( nr + 1, ir2 )

               if( rm(i) .eq. 0.0d0 ) goto 51

               c = ( xpp - x0 )**2
     &           + ( ypp - y0 )**2 - rm(i)**2

               dchk = b**2 - 4.0d0 * a * c

               if( dchk .lt. 0.0d0 ) goto 51

               t1 = ( -b + sqrt(dchk) ) / ( a * 2.0 )
               t2 = ( -b - sqrt(dchk) ) / ( a * 2.0 )

               if( t1 .lt. small ) t1 = dmax
               if( t2 .lt. small ) t2 = dmax

               ttr = min(t1,t2)

               if( ttr .lt. dr ) then

                  dr  = ttr
                  idr = i

               end if

   51       continue

         end if

*-----------------------------------------------------------------------
*        dz : distance to the nearest z mesh
*-----------------------------------------------------------------------

               dz = dmax

         if( izk .ne. 0 ) then
            if( izm + izk .ge. 1 .and.
     &          izm + izk .le. nz + 1 ) then

               dz = abs( ( zm(izm+izk) - zpp ) / ud(3) )

            end if
         end if

*-----------------------------------------------------------------------
*        which boundary is the nearlist
*-----------------------------------------------------------------------

            if( dr .gt. dz ) then

               dd = dz
               jz = 1

            else

               dd = dr
               jz = 0

            end if

            if( dd .gt. dmax0 ) return

*-----------------------------------------------------------------------
*        recalculate rng for ncol = 11 and nedisp
*-----------------------------------------------------------------------

            if( mat .gt. 0 .and. jtyp .ne. 0 .and. nedisp .ne. 0 .and.
     &           (  ityp .ne. 12 .and. ityp .ne. 13 ) ) then


               call rainge(se,rng,mat,ityp,ktyp,jtyp,rtyp)

               if( ncol .eq. 11 )  dis = tot + rng

            end if

*-----------------------------------------------------------------------
*        propagate position upto the boundary or the final point
*-----------------------------------------------------------------------

            if( tot + dd .ge. dis - small ) then

               dd  = dis - tot
               tot = dis
               xpp = xxc
               ypp = yyc
               zpp = zzc

            else

               tot = tot + dd
               xpp = xpp + dd * ud(1)
               ypp = ypp + dd * ud(2)
               zpp = zpp + dd * ud(3)

            end if

*-----------------------------------------------------------------------
*        final energy ee and range rng1
*-----------------------------------------------------------------------

            if( mat .gt. 0 ) then

               delt = dd

               call ecol(ee,delt,se,rng,
     &                   mat,ityp,ktyp,jtyp,rtyp)

               rng1 = max( 0.0d0, rng - delt )

               if( rng1 .gt. 1.d+30 ) goto 200

               ee = max( 0.0d0, ee )
               ee = min( se, ee )

               if( iegsemi .ne. 0 .and. ityp .eq. 14 ) then
                  se = se + edep
               end if

cc H.Iwase 2014/8/21 change under the assumption itype=14 with ncol=11 does not come here
! for the egs5 electron last hinge ( e=ecut but still have edep (dE) )

               if( iegsemi .ne. 0 .and.
     &           ( ncol .eq. 9 .or. ncol .eq. 11 ) ) then


                  ee = 0d0
                  if( ityp.eq.12 .or. ityp.eq.13 ) se = se + edep

               end if
            else

               ee = se

            end if

*-----------------------------------------------------------------------
*        booking
*-----------------------------------------------------------------------

            iz = izc
            ir = irc

      if( ir .ge. 1 .and. ir .lt. nr + 1 .and.
     &    iz .ge. 1 .and. iz .lt. nz + 1 ) then


c Takeshi Kai (track structure mode)(2018/11/30)
         if(
     .         ( ( ktyp.eq.6000012            ) .and.
     .             mat1 .eq. 1                  .and.
     .             ( ec(ibkec+no,ipomp+1) .le. ctsmax )  .and.
     .             ( ec(ibkec+no,ipomp+1) .ge. 0.d0 ) )
     .   .or.
     .         ( ( ityp.eq.1                  ) .and.
     .             mat1 .eq. 1                  .and.
     .             ( ec(ibkec+no,ipomp+1) .le. ptsmax )  .and.
     .             ( ec(ibkec+no,ipomp+1) .ge. 0.d0 ) )
     .   .or.
     .         ( ( ityp.eq.12 .or. ityp.eq.13 ) .and.
     .             mat1 .ne. 0                  .and.
     .             ( ec(ibkec+no,ipomp+1) .le. etsmax )  .and.
     .             ( ec(ibkec+no,ipomp+1) .ge. 0.d0 ) )
     .   .or.
     .         ( (  ibryf(ityp,ktyp) .ge. 1   ) .and.
     .             mat1 .ne. 0                  .and.
     .             ( ec(ibkec+no,ipomp+1) .le. tsmax
     .                       *dble(ibryf(ityp,ktyp)) )   .and.
     .             ( ec(ibkec+no,ipomp+1) .ge. 0.d0 )  )
     .   ) then

                if(mat1 .ne.  0) then
                    if( ityp .eq. 12 .or. ityp .eq. 13 ) ttsmin = etsmin
                    if( ibryf(ityp,ktyp) .ge. 1)
     &                 ttsmin = ibryf(ityp,ktyp) * emin(ityp)
                endif
                if( mat1 .eq.  1 ) then
                    if( ityp .eq.  1 )                   ttsmin = ptsmin
                    if( ktyp .eq. 6000012)               ttsmin = ctsmin
                endif

                if(ncol .eq. 10 .or. ncol .eq. 11)then
                  tlv = oldwt*(e(ibke+no,ipomp+1)-ec(ibkec+no,ipomp+1))
                elseif(ncol .eq. 13) then
                  tlv = oldwt * dexc_ene
                elseif(e(ibke+no,ipomp+1).le.ttsmin)then
                  tlv = oldwt * e(ibke+no,ipomp+1)
                elseif(ec(ibkec+no,ipomp+1).lt.ttsmin .and.
     &                      ibryf(ityp,ktyp) .ge. 1) then
                  tlv = 0.d0
                else
                  tlv = oldwt * dexc_ene
                endif

                if(itunt(m).eq.0.and.mat.gt.0)then
                   tlv = tlv / rhog(mat)
                endif

                  do ip = 1, ipn

                     trEVENT(ips(ip),it,ir,iz) =
     &               trEVENT(ips(ip),it,ir,iz) + tlv

                  end do

         if(sweight(m) .ne. 1.d30 .and. sweight(m) .ne. oldwt)
     &      write(*,9000)
         sweight(m) = oldwt ! T.Sato 2018/03/07  remember weight of this history
         goto 200
         end if

*-----------------------------------------------------------------------
*        dedxfnc is specified
*-----------------------------------------------------------------------

         if( itdfn(m,1) .ne. 0 ) then

               if( se .lt. eletb(1) ) goto 200

               do i = 1, nlete

                  if( se .ge. eletb(i) .and.
     &                se .lt. eletb(i+1) ) goto 30

               end do

   30             ie1 = min( i, nlete )

               do i = nlete, 1, -1

                  if( ee .ge. eletb(i) .and.
     &                ee .lt. eletb(i+1) ) goto 40

               end do

   40             ie2 = max( i, 1 )

*-----------------------------------------------------------------------
*           rrl : track length and rrr : charge particles correction
*-----------------------------------------------------------------------

                  rrl = dd
                  rrr = dd
                  rrt = 0.0d0

            if( mat .gt. 0 ) then

               if( se .gt. eletb(nlete+1) ) then

                  call rainge(eletb(nlete+1),rngm,mat,
     &                        ityp,ktyp,jtyp,rtyp)

                  rrr = max( 0.0d0, rrr - ( rng - rngm ) )
                  rrt = rrt + rng -rngm

               end if

            end if

*-----------------------------------------------------------------------
*           tally
*-----------------------------------------------------------------------

            do 270 ie = ie1, ie2, -1

                     emie = max( ee, eletb(ie) )
                     emiu = min( se, eletb(ie+1) )

                  if( emiu * emie .gt. 0.d0 ) then
                     erg  = sqrt( emiu * emie )
                  else
                     erg  = emiu / 2.d0
                  end if

                  if( mat .gt. 0 ) then

                     call rainge(emie,rng2,mat,ityp,ktyp,jtyp,rtyp)

                     rng2 = max( 0.0d0, rng2 - rng1 )
                     rrl  = max( 0.0d0, rrr  - rng2 )

                        rrt  = rrt + rrl

                     if( rrt .gt. dd ) then

                        rrl = rrl + ( dd - rrt )
                        rrt = dd

                     else if( ie .eq. ie2 .and.
     &                        ee .gt. eletb(ie) ) then

                        rrl = rrl + ( dd - rrt )

                     end if

                  end if

*-----------------------------------------------------------------------
*              LET bin
*-----------------------------------------------------------------------

               call dedxas(erg,dedx,lmat,ityp,ktyp,jtyp,rtyp)

                  dedx = dedx / 10.d0
                  dtrk = rrl
                  eini = emiu
                  efin = emie

               if( itdfn(m,1) .eq. 1 ) then

                  call usrdfn1(ityp,ktyp,jtyp,rtyp,
     &                         dedx,dtrk,eini,efin,dhet)

               else if( itdfn(m,1) .eq. 2 ) then

                  call usrdfn2(ityp,ktyp,jtyp,rtyp,
     &                         dedx,dtrk,eini,efin,dhet)

               end if

*-----------------------------------------------------------------------

                     tlv = dhet ! T.Sato 2018/01/20, delete oldwt

                  do ip = 1, ipn

                     trEVENT(ips(ip),it,ir,iz) =
     &               trEVENT(ips(ip),it,ir,iz) + tlv

                  end do

                  if(sweight(m) .ne. 1.d30 .and. sweight(m) .ne. oldwt)
     &              write(*,9000)
                  sweight(m) = oldwt ! T.Sato 2018/03/07  remember weight of this history

                  rrr = rrr - rrl

  270       continue

*-----------------------------------------------------------------------
*        normal case
*-----------------------------------------------------------------------

         else

                     tlv = ( se - ee )  ! T.Sato 2018/01/20, delete oldwt

                  do ip = 1, ipn

                     trEVENT(ips(ip),it,ir,iz) =
     &               trEVENT(ips(ip),it,ir,iz) + tlv

                  end do

                  if(sweight(m) .ne. 1.d30 .and. sweight(m) .ne. oldwt)
     &              write(*,9000)
                  sweight(m) = oldwt ! T.Sato 2018/03/07  remember weight of this history

         end if

*-----------------------------------------------------------------------

      end if

*-----------------------------------------------------------------------
*     next position
*     for z-crossing jz = 1, r-crossing jz = 0
*-----------------------------------------------------------------------

            if( jz .eq. 1 ) then

                  izm = izm + izk
                  izc = izc + izk

            else if( jz .eq. 0 ) then

               if( idr .gt. irc ) then

                  irc = idr

               else

                  irc = idr - 1

               end if

            end if

            se  = ee
            rng = rng1

*-----------------------------------------------------------------------

      goto 50

*-----------------------------------------------------------------------

  200       continue

*-----------------------------------------------------------------------
*     dead particles of produced charged particle
*-----------------------------------------------------------------------

      if( ( ncol .eq. 13 .or. ncol .eq. 14 ) .and.
     &      mat .gt. 0 .and. nclsts .gt. 0 ) then

*-----------------------------------------------------------------------
*        check of counter
*-----------------------------------------------------------------------

            do i = 1, 3

               if( itcnt(i,m) .ne. 0 ) then

                  if( jcount(i,1) .lt. itcnt(i*2+2,m) .or.
     &                jcount(i,1) .gt. itcnt(i*2+3,m) ) return

               end if

            end do

*-----------------------------------------------------------------------
*        check z mesh and r mesh region
*-----------------------------------------------------------------------

            if(  zzc .lt. zm(1) ) return
            if(  zzc .ge. zm(nz+1) ) return
            if( dis1 .lt. rm(1) ) return
            if( dis1 .ge. rm(nr + 1) ) return

*-----------------------------------------------------------------------
*           z-position
*-----------------------------------------------------------------------

               do i = 1, nz

                  if( zzc .ge. zm(i) .and.
     &                zzc .lt. zm(i+1) ) goto 39

               end do

   39          iz = i

*-----------------------------------------------------------------------
*           r-position
*-----------------------------------------------------------------------

               do i = 1, nr

                  if( dis1 .ge. rm(i) .and.
     &                dis1 .lt. rm(i+1) ) goto 48

               end do

   48          ir = i

*-----------------------------------------------------------------------
*        dead particle and nucleus except for electron positron
*-----------------------------------------------------------------------

         do j = 1, nclsts

         if( jclusts(4,j) .lt. 0 .and. jclusts(5,j) .ne. 0 .and.
     &       qclusts(7,j) .gt. eletb(1) ) then


                  ipart = jclusts(3,j)
                  kpart = jclusts(7,j)
                  jpart = jclusts(5,j)
                  rpart = qclusts(5,j) * 1000.d0

*-----------------------------------------------------------------------

            call pcheck(m,np,ipart,kpart,jpart,ipn,ips)

               if( ipn .eq. 0 )  goto 285

*-----------------------------------------------------------------------

                  ein = qclusts(7,j)
                  ecc = 0.0d0

*-----------------------------------------------------------------------
*        dedxfnc is specified
*-----------------------------------------------------------------------

         if( itdfn(m,1) .ne. 0 ) then

               if( ein .lt. eletb(1) ) goto 285

               do i = 1, nlete

                  if( ein .ge. eletb(i) .and.
     &                ein .lt. eletb(i+1) ) goto 35

               end do

   35             ie1 = min( i, nlete )
                  ie2 = 1

*-----------------------------------------------------------------------
*           rrl : track length and rrr : charge particles correction
*-----------------------------------------------------------------------

                  call rainge(ein,rng0,mat,
     &                        ipart,kpart,jpart,rpart)

                  if( rng0 .gt. 1.d+30 ) goto 285

                  tlngth = rng0
                  rng1 = 0.0d0

                  rrl = tlngth
                  rrr = tlngth
                  rrt = 0.0d0

               if( ein .gt. eletb(nlete+1) ) then

                  call rainge(eletb(nlete+1),rngm,mat,
     &                        ipart,kpart,jpart,rpart)

                  rrr = max( 0.0d0, rrr - ( rng0 - rngm ) )
                  rrt = rrt + rng0 -rngm

               end if

*-----------------------------------------------------------------------
*           tally
*-----------------------------------------------------------------------

            do 275 ie = ie1, ie2, -1

                     emie = max( ecc, eletb(ie) )
                     emiu = min( ein, eletb(ie+1) )

                  if( emiu * emie .gt. 0.d0 ) then
                     erg  = sqrt( emiu * emie )
                  else
                     erg  = emiu / 2.d0
                  end if

                  if( mat .gt. 0 ) then

                     call rainge(emie,rng,mat,
     &                           ipart,kpart,jpart,rpart)

                     rng = max( 0.0d0, rng - rng1 )
                     rrl = max( 0.0d0, rrr - rng )

                        rrt  = rrt + rrl

                     if( rrt .gt. tlngth ) then

                        rrl = rrl + ( tlngth - rrt )
                        rrt = tlngth

                     else if( ie .eq. ie2 .and.
     &                        ecc .gt. eletb(ie) ) then

                        rrl = rrl + ( tlngth - rrt )

                     end if

                  end if

*-----------------------------------------------------------------------
*              LET bin
*-----------------------------------------------------------------------

               call dedxas(erg,dedx,lmat,ipart,kpart,jpart,rpart)

                  dedx = dedx / 10.d0
                  dtrk = rrl
                  eini = emiu
                  efin = emie

               if( itdfn(m,1) .eq. 1 ) then

                  call usrdfn1(ipart,kpart,jpart,rpart,
     &                         dedx,dtrk,eini,efin,dhet)

               else if( itdfn(m,1) .eq. 2 ) then

                  call usrdfn2(ipart,kpart,jpart,rpart,
     &                         dedx,dtrk,eini,efin,dhet)

               end if

*-----------------------------------------------------------------------

                     tlv = dhet  ! T.Sato 2018/01/20, delete qclusts(8,j)

                  do ip = 1, ipn

                     trEVENT(ips(ip),it,ir,iz) =
     &               trEVENT(ips(ip),it,ir,iz) + tlv

                  end do

                  if(sweight(m) .ne. 1.d30 .and. sweight(m) .ne. oldwt)
     &              write(*,9000)
                  sweight(m) = oldwt ! T.Sato 2018/03/07  remember weight of this history

                  rrr = rrr - rrl

  275       continue

*-----------------------------------------------------------------------
*        normal case
*-----------------------------------------------------------------------

         else

                     tlv = ( ein - ecc )  ! T.Sato 2018/01/20, delete qclusts(8,j)

                  do ip = 1, ipn

                     trEVENT(ips(ip),it,ir,iz) =
     &               trEVENT(ips(ip),it,ir,iz) + tlv

                  end do

                  if(sweight(m) .ne. 1.d30 .and. sweight(m) .ne. oldwt)
     &              write(*,9000)
                  sweight(m) = oldwt ! T.Sato 2018/03/07  remember weight of this history

                  endif

*-----------------------------------------------------------------------

  285       continue

         end if
         end do

      end if

*-----------------------------------------------------------------------

      return

 9000 Format('[T-deposit] with output = deposit, particle weight is not
     &unique in one history. Result may be unreasonable.')

      end

************************************************************************
*                                                                      *
      subroutine tdepstxyzEVENT(ncol,m,nl,lt,np,nx,ny,nz,ne,nt,
     &                     xm,ym,zm,eb,tb,tr,trEVENT,tr0)
*                                                                      *
*       Deposit tally in xyz scoring mesh                              *
*       output = deposit                                               *
*       last modified by T.Furuta on 2012/05/08                        *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*      ncol  ..... reaction type                                       *
*              0 : end of batch                                        *
*              4 : source particle                                     *
*                                                                      *
*              9 : termination by time cut-off                         *
*             10 : geometry boundary crossing                          *
*             11 : termination by energy cut-off                       *
*             12 : termination by escape or leakage                    *
*             13 : nuclear reaction (n,x)                              *
*             14 : nuclear reaction (n,n'x)                            *
*             15 : sequential transport only for tally                 *
*                                                                      *
************************************************************************
      use MMBANKMOD   !FURUTA
      use EVENTTALMOD !FURUTA
      use partmod, only: itmxpt ! frtati 2021/10/05
*-----------------------------------------------------------------------

      implicit double precision( a-h, o-z )

*-----------------------------------------------------------------------

      include 'param00.inc'
      include 'param.inc'

*-----------------------------------------------------------------------

      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
!$OMP THREADPRIVATE(/icomon/)
      common /wtsave/ oldwt
!$OMP THREADPRIVATE(/wtsave/)
      common /tlgeom/ iblz1,iblz2
!$OMP THREADPRIVATE(/tlgeom/)

*-----------------------------------------------------------------------

      common /clustt/ nclsts, iclusts(nnn)
!$OMP THREADPRIVATE(/clustt/)
      common /clustw/ jclusts(0:8,nnn),  qclusts(0:12,nnn)
!$OMP THREADPRIVATE(/clustw/)
      common /cntcls/ jcount(3,nnn)
!$OMP THREADPRIVATE(/cntcls/)

*-----------------------------------------------------------------------

      common /tall01/ itmsh(itlmax), itunt(itlmax), itspc(itlmax),
     &                itout(itlmax), ittwo(itlmax)
      common /tall22/ itmcn(itlmax), itmtn(itlmax), itmtt(itlmax)
      common /tall36/ itdpo(itlmax)   ! S.Abe 2015/12/03
      common /tall37/ itcnt(9,itlmax)
      common /tall46/ itmtr(itlmax,4), rtmtr(itlmax,13)

      common /tall50/ itlmt(itlmax), ite2l(itlmax)   ! CCSE 2022.08.31
      common /tall51/ eletb(2000), nlete
      common /tall53/ itdfn(itlmax,2)

c T.Sato 2014/8/19 for detector resolution
      common /tall61/ rtdre(itlmax),rtdfa(itlmax)

      common /kmat1d/ idmn(0:kvlmax), idnm(kvmmax)

cABE 2016/03/04 for identification part=all
      common /tall65/ ipall(itlmax)

      common /tall82/ itcnth(9,itlmax)

*-----------------------------------------------------------------------

      dimension   lt(nl)
      dimension xm(nx+1)
      dimension ym(ny+1)
      dimension zm(nz+1)
      dimension eb(ne+1)
      dimension tb(nt+1)
      dimension tr(np,0:ne,nt,nx*ny*nz,2)
      dimension trEVENT(np,nt,nx*ny*nz)

*-----------------------------------------------------------------------

      dimension ud(3)

      data dmax  /1.0d+19/
      data dmax0 /1.0d+18/

      dimension ips(itmxpt) ! frtati 2021/10/05 6 -> itmxpt
      dimension facm(6)

*-----------------------------------------------------------------------

      common /paraj/ mstz(300), parz(300)

*-----------------------------------------------------------------------

      double precision, save :: sweight(itlmax)
!$OMP THREADPRIVATE(sweight)

      common /jcomon/ nabov,nobch,nocas,nomax
!$OMP THREADPRIVATE(/jcomon/)

*-----------------------------------------------------------------------

      real*8 edep
      common/EPCONTedep/ edep
!$OMP THREADPRIVATE(/EPCONTedep/)

      integer iegsemi, iegsout
      common /egsemi/ iegsemi, iegsout

C for nonshared_tally option
      integer italsh
      common /talsh/ italsh

*------------------------------------------------------------------------

      common /spred/ nspred, nwsprd, nedisp, itstep, ndedx

      common /regdc/  idrg(kvlmax), idgr(kvmmax)
      common /  tscmsg   / ktsc(kvlmax), mntsc, ntsc(kvlmax)
      common /  tscreg   / ntscell(kvlmax)
      common / etsminmax / etsmin, etsmax
      common / ptsminmax / ptsmin, ptsmax
      common / ctsminmax / ctsmin, ctsmax
      common / tsminmax  / tsmax
      common / etsexe  / dexc_ene
!$OMP THREADPRIVATE(/etsexe/)
      common /celdg/  rhog(kvlmax)

      common /eparm/  esmax, esmin, emin(20)

      common /cntmx/ ncntmx(3), nctmxsr(3)
!$OMP THREADPRIVATE(/cntmx/)

! sumover
      real(8) tr0(np,0:ne,nt,nx*ny*nz)

*-----------------------------------------------------------------------

      icf(ix,iy,iz) = ix + ( iy - 1 ) * nx + ( iz - 1 ) * nx * ny

*-----------------------------------------------------------------------

! sumover
      tr0(:,:,:,:) = 0.0d0

*-----------------------------------------------------------------------

! T.Sato 2014/8/21, for considering detector resolution, ncol = 101
      if(ncol.eq.101) then
       if(rtdre(m).eq.0.0.and.rtdfa(m).eq.0.0) return ! detector resolution = 0
       if( .not. FIRSTsrc ) then
        do it = 1, nt
         do ix = 1, nx
          do iy = 1, ny
           do iz = 1, nz
            gaurntmp=gaurn(dummy)  ! T.Sato 2018/01/25
            do ip = 1, np
             if( trEVENT(ip,it,icf(ix,iy,iz)).gt.0.0d0 ) then
              heats = trEVENT(ip,it,icf(ix,iy,iz))
              if(rtdre(m) .lt. 0.d0) then ! Ogawa 2018/11/05 if sigma<0, userdefined resolution
               heats = usrdefres(heat)  ! default (example) is gauss + exp asymmetric
              else
               heats = heats + gaurntmp*  ! T.Sato 2018/01/25
     &         sqrt(rtdre(m)**2+rtdfa(m)*heats)
              endif
              if(heats.lt.0.0) heats = 0.0
              trEVENT(ip,it,icf(ix,iy,iz))=heats
             endif
            enddo
           enddo
          enddo
         enddo
        enddo
       endif
       return
      endif
*----- End of revision on 2014/8/21 ------------------------------------


         small = parz(28) * 1.d-4 ! T.Sato 2023/11/01 avoid infinite loop

*-----------------------------------------------------------------------
* check of history counter
*-----------------------------------------------------------------------
         if ( ncol.eq.0 .or. ncol.eq.4 ) then
           do i = 1, 3
             if( itcnth(i,m) .eq. 1 ) then
               if( ncntmx(i) .lt. itcnth(i*2+2,m) .or.
     &             ncntmx(i) .gt. itcnth(i*2+3,m) ) then
                 trEVENT = 0.d0
               end if
             end if
           end do
         end if

*-----------------------------------------------------------------------
*     deposit : ncol = 0, 4 and ns = 4
*-----------------------------------------------------------------------

      if( ( ncol .eq. 0 .or. ncol .eq. 4 ) .and.
     &      itout(m) .ge. 2 ) then

         if( .not. FIRSTsrc ) then

            do ip = 1, np
            do it = 1, nt
            do ix = 1, nx
            do iy = 1, ny
            do iz = 1, nz

               if( itdpo(m) .eq. 0 .and.
     &             trEVENT(ipall(m),it,icf(ix,iy,iz)) .gt. 0.0 ) then
                  heats = trEVENT(ipall(m),it,icf(ix,iy,iz))
                  ratio = trEVENT(ip,it,icf(ix,iy,iz))
     &                   / trEVENT(ipall(m),it,icf(ix,iy,iz))
               else
                  heats = trEVENT(ip,it,icf(ix,iy,iz))
                  ratio = 1.0d0
               end if

               if( heats .gt. 0.0d0 .and.
     &             heats .ge. eb(1) .and. heats .lt. eb(ne+1) ) then

                  do ie = 1, ne

                     if( heats .ge. eb(ie) .and.
     &                   heats .lt. eb(ie+1) ) then
!$OMP CRITICAL (tr_tdepstxyzEVENT_1)
! sumover
                        tr0(ip,ie,it,icf(ix,iy,iz)) = 
     &                  tr0(ip,ie,it,icf(ix,iy,iz)) + sweight(m)*ratio
                        tr(ip,ie,it,icf(ix,iy,iz),1) =
     &                  tr(ip,ie,it,icf(ix,iy,iz),1) + sweight(m)*ratio   ! S.Abe 2015/12/03
!$OMP END CRITICAL (tr_tdepstxyzEVENT_1)

                        if(ip.eq.ipall(m).and.ratio.gt.0) then ! T.Sato 2022/09/09, part=all is always probability mode
!$OMP CRITICAL (tr_tdepstxyzEVENT_2)
                         tr(ip,ie,it,icf(ix,iy,iz),2) =
     &                            tr(ip,ie,it,icf(ix,iy,iz),2) + 1.d0            ! Ogawa 2020/04/22  if (weight < 1.d0), tr(ip,ie,ir,it,2) should be count number
!$OMP END CRITICAL (tr_tdepstxyzEVENT_2)
                        else
!$OMP CRITICAL (tr_tdepstxyzEVENT_3)
                         tr(ip,ie,it,icf(ix,iy,iz),2) =
     &                   tr(ip,ie,it,icf(ix,iy,iz),2)+(sweight(m)*ratio)
     &                        **2   ! S.Abe 2015/12/03
!$OMP END CRITICAL (tr_tdepstxyzEVENT_3)
                        endif

                     end if

                  end do
               end if

            end do
            end do
            end do
            end do
            end do

! sumover
           call tdepstxyz_sumover_ip(m,1, ipall(m), ratio,
     &                               np,  ne, nt, nx, ny, nz, tr0)
           tr0(:,:,:,:) = 0.0d0

            do ip = 1, np
            do it = 1, nt
            do ix = 1, nx
            do iy = 1, ny
            do iz = 1, nz
               trEVENT(ip,it,icf(ix,iy,iz)) = 0.0d0
            end do
            end do
            end do
            end do
            end do

         end if

*-----------------------------------------------------------------------
*        save initial weight and zero set
*-----------------------------------------------------------------------

         if( ncol .eq. 4 ) then

            sweight(m) = 1.d30 ! initialize

         end if

      end if

*-----------------------------------------------------------------------
*        check of ncol
*-----------------------------------------------------------------------

         if( ncol .lt. 9 ) return

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
*        check time
*-----------------------------------------------------------------------

               tin = abs(t(ibkt+no,ipomp+1))

               if( tin .lt. tb(1) ) return
               if( tin .ge. tb(nt+1) ) return

            do i = 1, nt

               if( tin .ge. tb(i) .and.
     &             tin .lt. tb(i+1) ) goto 501

            end do

  501          it = i

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
*        transform positions
*-----------------------------------------------------------------------

            call trnsxx(x(ibkx+no,ipomp+1),y(ibky+no,ipomp+1),
     &                  z(ibkz+no,ipomp+1),
     &                  xxa,yya,zza,itmtr(m,4))

            call trnsxx(xc(ibkxc+no,ipomp+1),yc(ibkyc+no,ipomp+1),
     &                  zc(ibkzc+no,ipomp+1),
     &                  xxc,yyc,zzc,itmtr(m,4))

*-----------------------------------------------------------------------
* add dummy step for EGS5 electrons at ncol=11
*     (stopped but still have dE)
cc H.Iwase 2014/8/21 add the change of H.Iwase 2014/1/9
*-----------------------------------------------------------------------

         if( iegsemi .ne. 0 .and.
     &     ( ityp .eq. 12 .or. ityp .eq. 13 ) .and. ncol .eq. 11 ) then

               dum = ( xxc - xxa )**2
     &             + ( yyc - yya )**2
     &             + ( zzc - zza )**2

            if( dum .le. small**2 ) then

               dumstep = small * 10.d0

               xdum = x(ibkx+no,ipomp+1) + u(ibku+no,ipomp+1)*dumstep
               ydum = y(ibky+no,ipomp+1) + v(ibkv+no,ipomp+1)*dumstep
               zdum = z(ibkz+no,ipomp+1) + w(ibkw+no,ipomp+1)*dumstep

               call trnsxx(xdum,ydum,zdum,xxc,yyc,zzc,itmtr(m,4))

            end if

         end if
*-----------------------------------------------------------------------
*        check position : out of rainge
*-----------------------------------------------------------------------

            if( xxa-small .lt. xm(1) .and.
     &          xxc-small .lt. xm(1) ) return

            if( yya-small .lt. ym(1) .and.
     &          yyc-small .lt. ym(1) ) return

            if( zza-small .lt. zm(1) .and.
     &          zzc-small .lt. zm(1) ) return

            if( xxa+small .ge. xm(nx+1) .and.
     &          xxc+small .ge. xm(nx+1) ) return

            if( yya+small .ge. ym(ny+1) .and.
     &          yyc+small .ge. ym(ny+1) ) return

            if( zza+small .ge. zm(nz+1) .and.
     &          zzc+small .ge. zm(nz+1) ) return

*-----------------------------------------------------------------------
*        check of particles
*-----------------------------------------------------------------------

            if( iegsemi .eq. 0 ) then

               if( jtyp .eq. 0 ) goto 200
               if( abs( ec(ibkec+no,ipomp+1) - e(ibke+no,ipomp+1) )
     &                   .eq. 0.0d0 )
     &             goto 200

             else

               if( jtyp .eq. 0 .and. ityp .ne. 14 ) goto 200
               if( ( abs( ec(ibkec+no,ipomp+1) - e(ibke+no,ipomp+1) )
     $                      .eq. 0.0d0 )
     $              .and.
     $          ( ityp .ne. 12 .and. ityp .ne. 13 .and. ityp .ne. 14))
     $            goto 200

             end if

            call pcheck(m,np,ityp,ktyp,jtyp,ipn,ips)

               if( ipn .eq. 0 ) goto 200

*-----------------------------------------------------------------------
*        check of counter
*-----------------------------------------------------------------------

            do i = 1, 3

               if( itcnt(i,m) .ne. 0 ) then

                  if( ncnt(ibknct+i,no,ipomp+1) .lt. itcnt(i*2+2,m) .or.
     &                ncnt(ibknct+i,no,ipomp+1) .gt. itcnt(i*2+3,m) )
     &               goto 200

               end if

            end do

*-----------------------------------------------------------------------
*        distance and unit vector ud(i)
*-----------------------------------------------------------------------

            dis = ( xxc - xxa )**2
     &          + ( yyc - yya )**2
     &          + ( zzc - zza )**2

c Takeshi Kai (track structure mode)(2018/11/30)
                mat1 = 0
            if(  mntsc .gt. 0 .and.
     &           ntscell( idgr(iblz(ibkblz+no,ipomp+1)) ) .ne. 0 .and.

     &        (( ityp .eq. 12 .or. ityp .eq. 13 ) .and.
     &         ( e(ibke+no,ipomp+1) .le. etsmax .and.
     &           e(ibke+no,ipomp+1) .ge. 0.d0 )   .or.

     &         ( ityp .eq. 1 ) .and.
     &         ( e(ibke+no,ipomp+1) .le. ptsmax .and.
     &           e(ibke+no,ipomp+1) .ge. 0.d0 )    .or.

     &         ( ktyp .eq. 6000012 ) .and.
     &         ( e(ibke+no,ipomp+1) .le. ctsmax .and.
     &           e(ibke+no,ipomp+1) .ge. 0.d0 )    .or.

     &         ( ibryf(ityp,ktyp) .ge. 1 ) .and.
     &         ( e(ibke+no,ipomp+1) .ge. 0.d0 .and.
     &           e(ibke+no,ipomp+1) .le. tsmax*dble(ibryf(ityp,ktyp)) ))

     &        ) mat1 = ntscell( idgr(iblz(ibkblz+no,ipomp+1)) )


         if(
     .         ( ( ktyp.eq.6000012            ) .and.
     .             mat1 .eq. 1                  .and.
     .             ( ec(ibkec+no,ipomp+1) .le. ctsmax )  .and.
     .             ( ec(ibkec+no,ipomp+1) .ge. 0.d0 ) )
     .   .or.
     .         ( ( ityp.eq.1                  ) .and.
     .             mat1 .eq. 1                  .and.
     .             ( ec(ibkec+no,ipomp+1) .le. ptsmax )  .and.
     .             ( ec(ibkec+no,ipomp+1) .ge. 0.d0 ) )
     .   .or.
     .         ( ( ityp.eq.12 .or. ityp.eq.13 ) .and.
     .             mat1 .ne. 0                  .and.
     .             ( ec(ibkec+no,ipomp+1) .le. etsmax )  .and.
     .             ( ec(ibkec+no,ipomp+1) .ge. 0.d0 ) )
     .   .or.
     .         ( (  ibryf(ityp,ktyp) .ge. 1   ) .and.
     .             mat1 .ne. 0                  .and.
     .             ( ec(ibkec+no,ipomp+1) .le. tsmax
     .                       *dble(ibryf(ityp,ktyp)) )   .and.
     .             ( ec(ibkec+no,ipomp+1) .ge. 0.d0 )  )
     .   ) then
           else
            if( dis .le. small**2 ) return ! org
           endif

            dis = sqrt(dis)

            ud(1) = ( xxc - xxa ) / dis
            ud(2) = ( yyc - yya ) / dis
            ud(3) = ( zzc - zza ) / dis

*-----------------------------------------------------------------------
*     initial position, energy and initial range
*-----------------------------------------------------------------------

            tot = 0.0d0
            rng1 = 0.0d0 !FURUTA

            se  = e(ibke+no,ipomp+1)
            xpp = xxa
            ypp = yya
            zpp = zza

            if( mat .gt. 0 )
     &          call rainge(se,rng,mat,ityp,ktyp,jtyp,rtyp)

*-----------------------------------------------------------------------
*     initial xyz position
*-----------------------------------------------------------------------

               ixm = 0
               ixc = 0

            if( ud(1) .gt. 0.0d0 ) ixk =  1
            if( ud(1) .eq. 0.0d0 ) ixk =  0
            if( ud(1) .lt. 0.0d0 ) ixk = -1

            if( ixk .ge. 0 ) then

               do i = 1, nx + 1
                  if( xm(i) .gt. xpp + small ) then
                         ixm = i - 1
                         ixc = i - 1
                         goto 32
                  end if
               end do

            else

               ixm = nx + 2
               ixc = nx + 1

               do i = 1, nx + 1
                  if( xm(i) .ge. xpp - small ) then
                         ixm = i
                         ixc = i - 1
                         goto 32
                  end if
               end do

            end if

   32       continue

*-----------------------------------------------------------------------

               iym = 0
               iyc = 0

            if( ud(2) .gt. 0.0d0 ) iyk =  1
            if( ud(2) .eq. 0.0d0 ) iyk =  0
            if( ud(2) .lt. 0.0d0 ) iyk = -1

            if( iyk .ge. 0 ) then

               do i = 1, ny + 1
                  if( ym(i) .gt. ypp + small ) then
                         iym = i - 1
                         iyc = i - 1
                         goto 33
                  end if
               end do

            else

               iym = ny + 2
               iyc = ny + 1

               do i = 1, ny + 1
                  if( ym(i) .ge. ypp - small ) then
                         iym = i
                         iyc = i - 1
                         goto 33
                  end if
               end do

            end if

   33       continue

*-----------------------------------------------------------------------

               izm = 0
               izc = 0

            if( ud(3) .gt. 0.0d0 ) izk =  1
            if( ud(3) .eq. 0.0d0 ) izk =  0
            if( ud(3) .lt. 0.0d0 ) izk = -1

            if( izk .ge. 0 ) then

               do i = 1, nz + 1
                  if( zm(i) .gt. zpp + small ) then
                         izm = i - 1
                         izc = i - 1
                         goto 34
                  end if
               end do

            else

               izm = nz + 2
               izc = nz + 1

               do i = 1, nz + 1
                  if( zm(i) .ge. zpp - small ) then
                         izm = i
                         izc = i - 1
                         goto 34
                  end if
               end do

            end if

   34       continue

*-----------------------------------------------------------------------
*     loop for finding mesh and booking upto total distance
*-----------------------------------------------------------------------

   50 continue

*-----------------------------------------------------------------------
*     calculation is finished
*-----------------------------------------------------------------------

         if( tot .ge. dis ) goto 200

*-----------------------------------------------------------------------
*        dx,dy,dz : distance to the nearest xyz mesh
*-----------------------------------------------------------------------

            dx = dmax
            dy = dmax
            dz = dmax

         if( ixk .ne. 0 ) then
            if( ixm + ixk .ge. 1 .and.
     &         ixm + ixk .le. nx + 1 ) then

               dx = abs( ( xm(ixm+ixk) - xpp ) / ud(1) )

            end if
         end if

         if( iyk .ne. 0 ) then
            if( iym + iyk .ge. 1 .and.
     &          iym + iyk .le. ny + 1 ) then

               dy = abs( ( ym(iym+iyk) - ypp ) / ud(2) )

            end if
         end if

         if( izk .ne. 0 ) then
            if( izm + izk .ge. 1 .and.
     &          izm + izk .le. nz + 1 ) then

               dz = abs( ( zm(izm+izk) - zpp ) / ud(3) )

            end if
         end if

*-----------------------------------------------------------------------
*        which boundary is the nearlist
*-----------------------------------------------------------------------

               dd = dx
               jj = 1

            if( dy .lt. dd ) then

               dd = dy
               jj = 2

            end if

            if( dz .lt. dd ) then

               dd = dz
               jj = 3

            end if

               tt = dd

            if( dd .gt. dmax0 ) return

*-----------------------------------------------------------------------
*        recalculate rng for ncol = 11 and nedisp
*-----------------------------------------------------------------------

            if( mat .gt. 0 .and. jtyp .ne. 0 .and. nedisp .ne. 0 .and.
     &           (  ityp .ne. 12 .and. ityp .ne. 13 ) ) then


               call rainge(se,rng,mat,ityp,ktyp,jtyp,rtyp)

               if( ncol .eq. 11 )  dis = tot + rng

            end if

*-----------------------------------------------------------------------
*        propagate position upto the boundary or the final point
*-----------------------------------------------------------------------

            if( tot + dd .ge. dis - small ) then

               dd  = dis - tot
               tot = dis
               xpp = xxc
               ypp = yyc
               zpp = zzc

            else

               tot = tot + dd
               xpp = xpp + tt * ud(1)
               ypp = ypp + tt * ud(2)
               zpp = zpp + tt * ud(3)

            end if

*-----------------------------------------------------------------------
*        final energy ee and range rng1
*-----------------------------------------------------------------------

            if( mat .gt. 0 ) then

               delt = dd

               call ecol(ee,delt,se,rng,
     &                   mat,ityp,ktyp,jtyp,rtyp)

               rng1 = max( 0.0d0, rng - delt )

               if( rng1 .gt. 1.d+30 ) goto 200

               ee = max( 0.0d0, ee )
               ee = min( se, ee )

               if( iegsemi .ne. 0 .and. ityp .eq. 14 ) then
                  se = se + edep
               end if

cc H.Iwase 2014/8/21 change under the assumption itype=14 with ncol=11 does not come here
! for the egs5 electron last hinge ( e=ecut but still have edep (dE) )

               if( iegsemi .ne. 0 .and.
     &           ( ncol .eq. 9 .or. ncol .eq. 11 ) ) then

                  if( ityp .eq. 12 .or. ityp .eq. 13 ) then

                        de = se - ee
                        se = se + de
                        ee = 0.0

                  else if( ityp .eq. 14 ) then

                        ee = 0.d0

                  end if

               end if
            else

               ee = se

            end if

*-----------------------------------------------------------------------
*        booking
*-----------------------------------------------------------------------

      if( ixc .ge. 1 .and. ixc .lt. nx + 1 .and.
     &    iyc .ge. 1 .and. iyc .lt. ny + 1 .and.
     &    izc .ge. 1 .and. izc .lt. nz + 1 ) then


c Takeshi Kai (track structure mode)(2018/11/30)
         if(
     .         ( ( ktyp.eq.6000012            ) .and.
     .             mat1 .eq. 1                  .and.
     .             ( ec(ibkec+no,ipomp+1) .le. ctsmax )  .and.
     .             ( ec(ibkec+no,ipomp+1) .ge. 0.d0 ) )
     .   .or.
     .         ( ( ityp.eq.1                  ) .and.
     .             mat1 .eq. 1                  .and.
     .             ( ec(ibkec+no,ipomp+1) .le. ptsmax )  .and.
     .             ( ec(ibkec+no,ipomp+1) .ge. 0.d0 ) )
     .   .or.
     .         ( ( ityp.eq.12 .or. ityp.eq.13 ) .and.
     .             mat1 .ne. 0                  .and.
     .             ( ec(ibkec+no,ipomp+1) .le. etsmax )  .and.
     .             ( ec(ibkec+no,ipomp+1) .ge. 0.d0 ) )
     .   .or.
     .         ( (  ibryf(ityp,ktyp) .ge. 1   ) .and.
     .             mat1 .ne. 0                  .and.
     .             ( ec(ibkec+no,ipomp+1) .le. tsmax
     .                       *dble(ibryf(ityp,ktyp)) )   .and.
     .             ( ec(ibkec+no,ipomp+1) .ge. 0.d0 )  )
     .   ) then

               if(mat1 .ne.  0) then
                   if( ityp .eq. 12 .or. ityp .eq. 13 ) ttsmin = etsmin
                   if( ibryf(ityp,ktyp) .ge. 1)
     &                ttsmin = ibryf(ityp,ktyp) * emin(ityp)
               endif
               if( mat1 .eq.  1 ) then
                   if( ityp .eq.  1 )                   ttsmin = ptsmin
                   if( ktyp .eq. 6000012)               ttsmin = ctsmin
               endif

               if(ncol .eq. 10 .or. ncol .eq. 11)then
                  tlv = oldwt*(e(ibke+no,ipomp+1)-ec(ibkec+no,ipomp+1))
               elseif(ncol .eq. 13) then
                  tlv = oldwt * dexc_ene
               elseif(e(ibke+no,ipomp+1).le.ttsmin)then
                  tlv = oldwt * e(ibke+no,ipomp+1)
               elseif(ec(ibkec+no,ipomp+1).lt.ttsmin .and.
     &                     ibryf(ityp,ktyp) .ge. 1) then
                  tlv = 0.d0
               else
                  tlv = oldwt * dexc_ene
               endif

               if(itunt(m).eq.0.and.mat.gt.0)then
                  tlv = tlv / rhog(mat)
               endif

               do ip = 1, ipn

                  trEVENT(ips(ip),it,icf(ixc,iyc,izc)) =
     &            trEVENT(ips(ip),it,icf(ixc,iyc,izc)) + tlv

               end do

         if(sweight(m) .ne. 1.d30 .and. sweight(m) .ne. oldwt)
     &      write(*,9000)
         sweight(m) = oldwt ! T.Sato 2018/03/07  remember weight of this history
         goto 200
         end if


*-----------------------------------------------------------------------
*        dedxfnc is specified
*-----------------------------------------------------------------------

         if( itdfn(m,1) .ne. 0 ) then

            if( se .lt. eletb(1) ) goto 200

            do i = 1, nlete

               if( se .ge. eletb(i) .and.
     &             se .lt. eletb(i+1) ) goto 30

            end do

   30          ie1 = min( i, nlete )

            do i = nlete, 1, -1

               if( ee .ge. eletb(i) .and.
     &             ee .lt. eletb(i+1) ) goto 40

            end do

   40          ie2 = max( i, 1 )

*-----------------------------------------------------------------------
*           rrl : track length and rrr : charge particles correction
*-----------------------------------------------------------------------

                  rrl = dd
                  rrr = dd
                  rrt = 0.0d0

            if( mat .gt. 0 ) then

               if( se .gt. eletb(ne+1) ) then

                  call rainge(eletb(nlete+1),rngm,mat,
     &                        ityp,ktyp,jtyp,rtyp)

                  rrr = max( 0.0d0, rrr - ( rng - rngm ) )
                  rrt = rrt + rng -rngm

               end if

            end if

*-----------------------------------------------------------------------
*           tally
*-----------------------------------------------------------------------

            do 270 ie = ie1, ie2, -1

                     emie = max( ee, eletb(ie) )
                     emiu = min( se, eletb(ie+1) )

                  if( emiu * emie .gt. 0.d0 ) then
                     erg  = sqrt( emiu * emie )
                  else
                     erg  = emiu / 2.d0
                  end if

                  if( mat .gt. 0 ) then

                     call rainge(emie,rng2,mat,ityp,ktyp,jtyp,rtyp)

                     rng2 = max( 0.0d0, rng2 - rng1 )
                     rrl  = max( 0.0d0, rrr  - rng2 )

                        rrt  = rrt + rrl

                     if( rrt .gt. dd ) then

                        rrl = rrl + ( dd - rrt )
                        rrt = dd

                     else if( ie .eq. ie2 .and.
     &                        ee .gt. eletb(ie) ) then

                        rrl = rrl + ( dd - rrt )

                     end if

                  end if

*-----------------------------------------------------------------------
*              LET bin
*-----------------------------------------------------------------------

               call dedxas(erg,dedx,lmat,ityp,ktyp,jtyp,rtyp)

                  dedx = dedx / 10.d0
                  dtrk = rrl
                  eini = emiu
                  efin = emie

               if( itdfn(m,1) .eq. 1 ) then

                  call usrdfn1(ityp,ktyp,jtyp,rtyp,
     &                         dedx,dtrk,eini,efin,dhet)

               else if( itdfn(m,1) .eq. 2 ) then

                  call usrdfn2(ityp,ktyp,jtyp,rtyp,
     &                         dedx,dtrk,eini,efin,dhet)

               end if

*-----------------------------------------------------------------------

                  tlv = dhet  ! T.Sato 2018/01/20, delete oldwt

               do ip = 1, ipn

                  trEVENT(ips(ip),it,icf(ixc,iyc,izc)) =
     &            trEVENT(ips(ip),it,icf(ixc,iyc,izc)) + tlv

               end do

               if(sweight(m) .ne. 1.d30 .and. sweight(m) .ne. oldwt)
     &              write(*,9000)
               sweight(m) = oldwt ! T.Sato 2018/03/07  remember weight of this history

                  rrr = rrr - rrl

  270       continue

*-----------------------------------------------------------------------
*        normal case
*-----------------------------------------------------------------------

         else

                  tlv = ( se - ee ) ! T.Sato 2018/01/20, delete oldwt

               do ip = 1, ipn

                  trEVENT(ips(ip),it,icf(ixc,iyc,izc)) =
     &            trEVENT(ips(ip),it,icf(ixc,iyc,izc)) + tlv

               end do

               if(sweight(m) .ne. 1.d30 .and. sweight(m) .ne. oldwt)
     &              write(*,9000)
               sweight(m) = oldwt ! T.Sato 2018/03/07  remember weight of this history

         end if

*-----------------------------------------------------------------------

      end if

*-----------------------------------------------------------------------
*     next position
*-----------------------------------------------------------------------

            if( jj .eq. 1 ) then

               ixm = ixm + ixk
               ixc = ixc + ixk

            else if( jj .eq. 2 ) then

               iym = iym + iyk
               iyc = iyc + iyk

            else if( jj .eq. 3 ) then

               izm = izm + izk
               izc = izc + izk

            end if

               se  = ee
               rng = rng1

*-----------------------------------------------------------------------

      goto 50

*-----------------------------------------------------------------------

  200       continue

*-----------------------------------------------------------------------
*     dead particles of produced charged particle
*-----------------------------------------------------------------------

      if( ( ncol .eq. 13 .or. ncol .eq. 14 ) .and.
     &      mat .gt. 0 .and. nclsts .gt. 0 ) then

*-----------------------------------------------------------------------
*        check of counter
*-----------------------------------------------------------------------

            do i = 1, 3

               if( itcnt(i,m) .ne. 0 ) then

                  if( jcount(i,1) .lt. itcnt(i*2+2,m) .or.
     &                jcount(i,1) .gt. itcnt(i*2+3,m) ) return

               end if

            end do

*-----------------------------------------------------------------------
*        check position : out of rainge
*-----------------------------------------------------------------------

            if(  xxc .lt. xm(1) ) return
            if(  xxc .ge. xm(nx+1) ) return

            if(  yyc .lt. ym(1) ) return
            if(  yyc .ge. ym(ny+1) ) return

            if(  zzc .lt. zm(1) ) return
            if(  zzc .ge. zm(nz+1) ) return

*-----------------------------------------------------------------------
*           x-position
*-----------------------------------------------------------------------

               do i = 1, nx

                  if( xxc .ge. xm(i) .and.
     &                xxc .lt. xm(i+1) ) goto 36

               end do

   36          ix = i

*-----------------------------------------------------------------------
*           y-position
*-----------------------------------------------------------------------

               do i = 1, ny

                  if( yyc .ge. ym(i) .and.
     &                yyc .lt. ym(i+1) ) goto 37

               end do

   37          iy = i

*-----------------------------------------------------------------------
*           z-position
*-----------------------------------------------------------------------

               do i = 1, nz

                  if( zzc .ge. zm(i) .and.
     &                zzc .lt. zm(i+1) ) goto 38

               end do

   38          iz = i

*-----------------------------------------------------------------------
*        dead particle and nucleus except for electron positron
*-----------------------------------------------------------------------

         do j = 1, nclsts

         if( jclusts(4,j) .lt. 0 .and. jclusts(5,j) .ne. 0 .and.
     &       qclusts(7,j) .gt. eletb(1) ) then


                  ipart = jclusts(3,j)
                  kpart = jclusts(7,j)
                  jpart = jclusts(5,j)
                  rpart = qclusts(5,j) * 1000.d0

*-----------------------------------------------------------------------

            call pcheck(m,np,ipart,kpart,jpart,ipn,ips)

               if( ipn .eq. 0 )  goto 285

*-----------------------------------------------------------------------

                  ein = qclusts(7,j)
                  ecc = 0.0d0

*-----------------------------------------------------------------------
*        dedxfnc is specified
*-----------------------------------------------------------------------

         if( itdfn(m,1) .ne. 0 ) then

               if( ein .lt. eletb(1) ) goto 285

               do i = 1, nlete

                  if( ein .ge. eletb(i) .and.
     &                ein .lt. eletb(i+1) ) goto 35

               end do

   35             ie1 = min( i, nlete )
                  ie2 = 1

*-----------------------------------------------------------------------
*           rrl : track length and rrr : charge particles correction
*-----------------------------------------------------------------------

                  call rainge(ein,rng0,mat,
     &                        ipart,kpart,jpart,rpart)

                  if( rng0 .gt. 1.d+30 ) goto 285

                  tlngth = rng0
                  rng1 = 0.0d0

                  rrl = tlngth
                  rrr = tlngth
                  rrt = 0.0d0

               if( ein .gt. eletb(nlete+1) ) then

                  call rainge(eletb(nlete+1),rngm,mat,
     &                        ipart,kpart,jpart,rpart)

                  rrr = max( 0.0d0, rrr - ( rng0 - rngm ) )
                  rrt = rrt + rng0 -rngm

               end if

*-----------------------------------------------------------------------
*           tally
*-----------------------------------------------------------------------

            do 275 ie = ie1, ie2, -1

                     emie = max( ecc, eletb(ie) )
                     emiu = min( ein, eletb(ie+1) )

                  if( emiu * emie .gt. 0.d0 ) then
                     erg  = sqrt( emiu * emie )
                  else
                     erg  = emiu / 2.d0
                  end if

                  if( mat .gt. 0 ) then

                     call rainge(emie,rng,mat,
     &                           ipart,kpart,jpart,rpart)

                     rng = max( 0.0d0, rng - rng1 )
                     rrl = max( 0.0d0, rrr - rng )

                        rrt  = rrt + rrl

                     if( rrt .gt. tlngth ) then

                        rrl = rrl + ( tlngth - rrt )
                        rrt = tlngth

                     else if( ie .eq. ie2 .and.
     &                        ecc .gt. eletb(ie) ) then

                        rrl = rrl + ( tlngth - rrt )

                     end if

                  end if

*-----------------------------------------------------------------------
*              LET bin
*-----------------------------------------------------------------------

               call dedxas(erg,dedx,lmat,ipart,kpart,jpart,rpart)

                  dedx = dedx / 10.d0
                  dtrk = rrl
                  eini = emiu
                  efin = emie

               if( itdfn(m,1) .eq. 1 ) then

                  call usrdfn1(ipart,kpart,jpart,rpart,
     &                         dedx,dtrk,eini,efin,dhet)

               else if( itdfn(m,1) .eq. 2 ) then

                  call usrdfn2(ipart,kpart,jpart,rpart,
     &                         dedx,dtrk,eini,efin,dhet)

               end if

*-----------------------------------------------------------------------

                     tlv = dhet  ! T.Sato 2018/01/20, delete qclusts(8,j)

                  do ip = 1, ipn

                     trEVENT(ips(ip),it,icf(ix,iy,iz)) =
     &               trEVENT(ips(ip),it,icf(ix,iy,iz)) + tlv

                  end do

                  if(sweight(m) .ne. 1.d30 .and. sweight(m) .ne. oldwt)
     &              write(*,9000)
                  sweight(m) = oldwt ! T.Sato 2018/03/07  remember weight of this history

                  rrr = rrr - rrl

  275       continue

*-----------------------------------------------------------------------
*        normal case
*-----------------------------------------------------------------------

         else

                     tlv = ( ein - ecc )  ! T.Sato 2018/01/20, delete qclusts(8,j)

                  do ip = 1, ipn

                     trEVENT(ips(ip),it,icf(ix,iy,iz)) =
     &               trEVENT(ips(ip),it,icf(ix,iy,iz)) + tlv

                  end do

                  if(sweight(m) .ne. 1.d30 .and. sweight(m) .ne. oldwt)
     &              write(*,9000)
                  sweight(m) = oldwt ! T.Sato 2018/03/07  remember weight of this history

                  endif

*-----------------------------------------------------------------------

  285       continue

         end if
         end do

      end if

*-----------------------------------------------------------------------

      return

 9000 Format('[T-deposit] with output = deposit, particle weight is not
     &unique in one history. Result may be unreasonable.')
      end

************************************************************************
*                                                                      *
      subroutine tdpst2regEVENT(ncol,m,np,nr,mr,ne1,ne2,nt,kr,eb1,eb2,
     &     tb,tr,trEVENT,tr0)
*                                                                      *
*       Deposit2 tally in region mesh                                  *
*       last modified by T.Furuta on 2012/05/08                        *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*      ncol  ..... reaction type                                       *
*              0 : end of batch                                        *
*              4 : source particle                                     *
*                                                                      *
*              9 : termination by time cut-off                         *
*             10 : geometry boundary crossing                          *
*             11 : termination by energy cut-off                       *
*             12 : termination by escape or leakage                    *
*             13 : nuclear reaction (n,x)                              *
*             14 : nuclear reaction (n,n'x)                            *
*             15 : sequential transport only for tally                 *
*                                                                      *
************************************************************************
      use MMBANKMOD   !FURUTA
      use EVENTTALMOD !FURUTA
      use partmod, only: itmxpt ! frtati 2021/10/05
*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param00.inc'
      include 'param.inc'

*-----------------------------------------------------------------------

      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
!$OMP THREADPRIVATE(/icomon/)
      common /wtsave/ oldwt
!$OMP THREADPRIVATE(/wtsave/)
      common /tlgeom/ iblz1,iblz2
!$OMP THREADPRIVATE(/tlgeom/)
      common /tlglat/ ilev1,ilev2,ilat1(5,10),ilat2(5,10)
!$OMP THREADPRIVATE(/tlglat/)

*-----------------------------------------------------------------------

      common /clustt/ nclsts, iclusts(nnn)
!$OMP THREADPRIVATE(/clustt/)
      common /clustw/ jclusts(0:8,nnn),  qclusts(0:12,nnn)
!$OMP THREADPRIVATE(/clustw/)
      common /cntcls/ jcount(3,nnn)
!$OMP THREADPRIVATE(/cntcls/)

*-----------------------------------------------------------------------

      common /tall37/ itcnt(9,itlmax)

      common /tall50/ itlmt(itlmax), ite2l(itlmax)   ! CCSE 2022.08.31
      common /tall55/ itlmt2(itlmax)

      common /tall53/ itdfn(itlmax,2)
      common /tall54/ itdfn2(itlmax,2)

      common /tall51/ eletb(2000), nlete

c T.Sato 2023/08/19 for detector resolution
      common /tall611/rtdr2(itlmax,2),rtdf2(itlmax,2) ! T.Sato 2023/08/19 for [t-deposit2]

      common /tall82/ itcnth(9,itlmax)

      common /kmat1d/ idmn(0:kvlmax), idnm(kvmmax)

*-----------------------------------------------------------------------

      dimension   kr(mr)
      dimension   eb1(ne1+1)
      dimension   eb2(ne2+1)
      dimension   tb(nt+1)
      dimension   tr(np,0:ne1,0:ne2,nt,2)
      dimension   trEVENT(np,2,nt)

      dimension ips(itmxpt) ! frtati 2021/10/05 6 -> itmxpt

*-----------------------------------------------------------------------

      double precision, save :: sweight(itlmax)
!$OMP THREADPRIVATE(sweight)
      common /jcomon/ nabov,nobch,nocas,nomax
!$OMP THREADPRIVATE(/jcomon/)

*-----------------------------------------------------------------------
      common /dedxfac/ dedxfd
!$OMP THREADPRIVATE(/dedxfac/)

*-----------------------------------------------------------------------

      real*8 edep
      common/EPCONTedep/ edep
!$OMP THREADPRIVATE(/EPCONTedep/)

      integer iegsemi, iegsout
      common /egsemi/ iegsemi, iegsout

C for nonshared_tally option
      integer italsh
      common /talsh/ italsh

      common /cntmx/ ncntmx(3), nctmxsr(3)
!$OMP THREADPRIVATE(/cntmx/)

! sumover
      real(8) tr0(np,0:ne1,0:ne2,nt)

! sumover
      tr0(:,:,:,:) = 0.0d0


      rng1=0.0 ! T.Sato 2024/04/02

! T.Sato 2023/8/19, for considering detector resolution, ncol = 101
      if(ncol.eq.101) then
       if(rtdr2(m,1).eq.0.0.and.rtdf2(m,1).eq.0.0.and.
     &    rtdr2(m,2).eq.0.0.and.rtdf2(m,2).eq.0.0) return ! detector resolution = 0
       if( .not. FIRSTsrc ) then
        do it = 1, nt
         do ir = 1, 2 ! always 2
          gaurntmp=gaurn(dummy)
          do ip = 1, np
           if( trEVENT(ip,ir,it).gt.0.0d0 ) then
            heats = trEVENT(ip,ir,it)
            if(rtdr2(m,1) .lt. 0.d0) then ! Ogawa 2018/11/05 if sigma<0, gauss + exp asymmetric resolution
             heats = usrdefres(heat)
            else
             heats = heats + gaurntmp*  ! T.Sato 2018/01/25
     &       sqrt(rtdr2(m,ir)**2+rtdf2(m,ir)*heats)
            endif
            if(heats.lt.0.0) heats = 0.0
            trEVENT(ip,ir,it)=heats
           endif
          enddo
         enddo
        enddo
       endif
       return
      endif

*-----------------------------------------------------------------------
* check of history counter
*-----------------------------------------------------------------------
         if ( ncol.eq.0 .or. ncol.eq.4 ) then
           do i = 1, 3
             if( itcnth(i,m) .eq. 1 ) then
               if( ncntmx(i) .lt. itcnth(i*2+2,m) .or.
     &             ncntmx(i) .gt. itcnth(i*2+3,m) ) then
                 sweight(m) = 0.d0
               end if
             end if
           end do
         end if

*-----------------------------------------------------------------------
*     deposit : ncol = 0, 4 and ns = 4
*-----------------------------------------------------------------------

      if( ncol .eq. 0 .or. ncol .eq. 4 ) then

         if( .not. FIRSTsrc ) then

            do ip = 1, np
            do it = 1, nt

                  heat1 = trEVENT(ip,1,it)
                  heat2 = trEVENT(ip,2,it)
                  trEVENT(ip,1,it) = 0.0d0
                  trEVENT(ip,2,it) = 0.0d0

               if( ( heat1 .gt.  0.0d0 .or. heat2 .gt. 0.0d0 ) .and.
     &             heat1 .ge. eb1(1) .and. heat1 .lt. eb1(ne1+1) .and.
     &             heat2 .ge. eb2(1) .and. heat2 .lt. eb2(ne2+1) ) then

                  do ie1 = 1, ne1
                  do ie2 = 1, ne2

                     if( heat1 .ge. eb1(ie1) .and.
     &                   heat1 .lt. eb1(ie1+1) .and.
     &                   heat2 .ge. eb2(ie2) .and.
     &                   heat2 .lt. eb2(ie2+1) ) then
                      if (italsh .eq. 0 ) then
! sumover
                      tr0(ip,ie1,ie2,it) = tr0(ip,ie1,ie2,it)
     &                                    + sweight(m)
                        tr(ip,ie1,ie2,it,1) = tr(ip,ie1,ie2,it,1)
     &                                      + sweight(m)
                        tr(ip,ie1,ie2,it,2) = tr(ip,ie1,ie2,it,2)
     &                                      + 1
                      else
!$OMP CRITICAL (tr_tdpst2regEVENT)
! sumover
                      tr0(ip,ie1,ie2,it) = tr0(ip,ie1,ie2,it)
     &                                    + sweight(m)
                        tr(ip,ie1,ie2,it,1) = tr(ip,ie1,ie2,it,1)
     &                                      + sweight(m)
                        tr(ip,ie1,ie2,it,2) = tr(ip,ie1,ie2,it,2)
     &                                      + 1
!$OMP END CRITICAL (tr_tdpst2regEVENT)

                      end if

                     end if

                  end do
                  end do

               end if

            end do
            end do

! sumover
           call tdpst2reg_sumover(m, 0,
     &                            np,  ne1, ne2, nt, tr0)
           tr0(:,:,:,:) = 0.0d0


         end if

*-----------------------------------------------------------------------
*        save initial weight and zero set
*-----------------------------------------------------------------------

         if( ncol .eq. 4 ) then

            sweight(m) = 1.d30 ! initialize

         end if

      end if

*-----------------------------------------------------------------------
*        check of ncol
*-----------------------------------------------------------------------

         if( ncol .lt. 9 ) return

*-----------------------------------------------------------------------
*        check region
*-----------------------------------------------------------------------

               jj = 0

      do 1000 ii = 1, nr

               call tregck(iblz1,ilev1,ilat1,mr,kr,jj,icc)

               if( icc .eq. 0 ) goto 1000

               ir = ii

*-----------------------------------------------------------------------
*        ir = 1, 2 only
*-----------------------------------------------------------------------

               if( ir .eq. 1 ) then

                  itlmat = itlmt(m)
                  itdfun = itdfn(m,1)

                  i1 = 1

               else

                  itlmat = itlmt2(m)
                  itdfun = itdfn2(m,1)

                  i1 = 2

               end if

*-----------------------------------------------------------------------
*        check time
*-----------------------------------------------------------------------

               tin = abs(t(ibkt+no,ipomp+1))

               if( tin .lt. tb(1) ) return
               if( tin .ge. tb(nt+1) ) return

            do i = 1, nt

               if( tin .ge. tb(i) .and.
     &             tin .lt. tb(i+1) ) goto 50

            end do

   50          it = i

*-----------------------------------------------------------------------
*        material for LET
*-----------------------------------------------------------------------

               if( itlmat .gt. 0 ) then

                  lmat = idnm( itlmat )

               else if( itlmat .eq. 0 ) then

                  lmat = mat

               else

                  lmat = -idnm( -itlmt(m) )

               end if

*-----------------------------------------------------------------------
*           track length
*-----------------------------------------------------------------------

               tlngth = sqrt( ( xc(ibkxc+no,ipomp+1) -
     &                                    x(ibkx+no,ipomp+1) )**2
     &                      + ( yc(ibkyc+no,ipomp+1) -
     &                                    y(ibky+no,ipomp+1) )**2
     &                      + ( zc(ibkzc+no,ipomp+1) -
     &                                    z(ibkz+no,ipomp+1) )**2 )

*-----------------------------------------------------------------------
*           initial and final energy
*-----------------------------------------------------------------------

                  ein = e(ibke+no,ipomp+1)
                  ecc = ec(ibkec+no,ipomp+1)

                  ecc = ein - dedxfd * ( ein - ecc )

                  if( ncol .eq. 9 .or. ncol .eq. 11 ) ecc = 0.d0

*-----------------------------------------------------------------------
cKN Iwase 2014/08/22 for EGS

               if( ( ityp.eq.12 .or. ityp.eq.13 .or. ityp.eq.14 ) .and.
     &               iegsemi .ne. 0 ) then

                     call egs5edxde(ecc,tlngth,ein,mat,ityp)

                  if(ncol.eq.9 .or. ncol.eq.11) then
                     ecc = 0d0
                     if( ityp.eq.12 .or. ityp.eq.13 ) ein = ein + edep
                  endif

               endif

*-----------------------------------------------------------------------
*        dedxfnc is specified
*-----------------------------------------------------------------------

         if( itdfun .ne. 0 ) then

               do i = 1, nlete

                  if( ein .ge. eletb(i) .and.
     &                ein .lt. eletb(i+1) ) goto 30

               end do

   30             ie1 = min( i, nlete )

               do i = nlete, 1, -1

                  if( ecc .ge. eletb(i) .and.
     &                ecc .lt. eletb(i+1) ) goto 40

               end do

   40             ie2 = max( i, 1 )

*-----------------------------------------------------------------------
*           rrl : track length and rrr : charge particles correction
*-----------------------------------------------------------------------

               rrl = tlngth
               rrr = tlngth
               rrt = 0.0d0

            if( mat .gt. 0 .and. jtyp .ne. 0 ) then

                  call rainge(ein,rng0,mat,ityp,ktyp,jtyp,rtyp)

                  rng1 = max( 0.0d0, rng0 - tlngth )

                  if( rng1 .gt. 1.d+30 ) return

               if( ein .gt. eletb(nlete+1) ) then

                  call rainge(eletb(nlete+1),rngm,mat,
     &                        ityp,ktyp,jtyp,rtyp)

                  rrr = max( 0.0d0, rrr - ( rng0 - rngm ) )
                  rrt = rrt + rng0 -rngm

               end if

            end if

               rrl0 = rrl
               rrr0 = rrr
               rrt0 = rrt

         end if

*-----------------------------------------------------------------------
*        check of particles, only charged particles
*-----------------------------------------------------------------------
            if( iegsemi .eq. 0 ) then

               if( jtyp .eq. 0 ) goto 200

             else

               if( jtyp .eq. 0 .and. ityp .ne. 14 ) goto 200

             end if

            call pcheck(m,np,ityp,ktyp,jtyp,ipn,ips)

               if( ipn .eq. 0 )  goto 200

*-----------------------------------------------------------------------
*        check of counter
*-----------------------------------------------------------------------

            do i = 1, 3

               if( itcnt(i,m) .ne. 0 ) then

                  if( ncnt(ibknct+i,no,ipomp+1) .lt. itcnt(i*2+2,m) .or.
     &                ncnt(ibknct+i,no,ipomp+1) .gt. itcnt(i*2+3,m) )
     &                 goto 200

               end if

            end do

*-----------------------------------------------------------------------
*        dedxfnc is specified
*-----------------------------------------------------------------------

         if( itdfun .ne. 0 ) then

               if( ein .lt. eletb(1) ) goto 200
               if( ecc .ge. eletb(nlete+1) ) goto 200

*-----------------------------------------------------------------------

               rrl = rrl0
               rrr = rrr0
               rrt = rrt0

*-----------------------------------------------------------------------
*           tally
*-----------------------------------------------------------------------

            do 270 ie = ie1, ie2, -1

                     emie = max( ecc, eletb(ie) )
                     emiu = min( ein, eletb(ie+1) )

                  if( emiu * emie .gt. 0.d0 ) then
                     erg  = sqrt( emiu * emie )
                  else
                     erg  = emiu / 2.d0
                  end if

                  if( mat .gt. 0 ) then

                     call rainge(emie,rng,mat,ityp,ktyp,jtyp,rtyp)

                     rng = max( 0.0d0, rng - rng1 )
                     rrl = max( 0.0d0, rrr - rng )

                        rrt  = rrt + rrl

                     if( rrt .gt. tlngth ) then

                        rrl = rrl + ( tlngth - rrt )
                        rrt = tlngth

                     else if( ie .eq. ie2 .and.
     &                        ecc .gt. eletb(ie) ) then

                        rrl = rrl + ( tlngth - rrt )

                     end if

                  end if

*-----------------------------------------------------------------------
*              dedx value multiplyed user factor
*-----------------------------------------------------------------------

               call dedxas(erg,dedx,lmat,ityp,ktyp,jtyp,rtyp)

                  dedx = dedx / 10.d0
                  dtrk = rrl
                  eini = emiu
                  efin = emie

               if( itdfun .eq. 1 ) then

                  call usrdfn1(ityp,ktyp,jtyp,rtyp,
     &                         dedx,dtrk,eini,efin,dhet)

               else if( itdfun .eq. 2 ) then

                  call usrdfn2(ityp,ktyp,jtyp,rtyp,
     &                         dedx,dtrk,eini,efin,dhet)

               end if

*-----------------------------------------------------------------------

                     tlv = dhet  ! T.Sato 2018/01/20, delete oldwt

                  do ip = 1, ipn

                     trEVENT(ips(ip),i1,it) =
     &               trEVENT(ips(ip),i1,it) + tlv

                  end do

                  if(sweight(m) .ne. 1.d30 .and. sweight(m) .ne. oldwt)
     &              write(*,9000)
                  sweight(m) = oldwt ! T.Sato 2018/03/07  remember weight of this history

                  rrr = rrr - rrl

  270       continue

*-----------------------------------------------------------------------
*        normal case
*-----------------------------------------------------------------------

         else

                     tlv = ( ein - ecc ) ! T.Sato 2018/01/20, delete oldwt

                  do ip = 1, ipn

                     trEVENT(ips(ip),i1,it) =
     &               trEVENT(ips(ip),i1,it) + tlv

                  end do

                  if(sweight(m) .ne. 1.d30 .and. sweight(m) .ne. oldwt)
     &              write(*,9000)
                  sweight(m) = oldwt ! T.Sato 2018/03/07  remember weight of this history

         end if

*-----------------------------------------------------------------------

  200       continue

*-----------------------------------------------------------------------
*     dead particles of produced charged particle
*-----------------------------------------------------------------------

      if( ( ncol .eq. 13 .or. ncol .eq. 14 ) .and.
     &      mat .gt. 0 .and. nclsts .gt. 0 ) then

*-----------------------------------------------------------------------
*        check of counter
*-----------------------------------------------------------------------

            do i = 1, 3

               if( itcnt(i,m) .ne. 0 ) then

                  if( jcount(i,1) .lt. itcnt(i*2+2,m) .or.
     &                jcount(i,1) .gt. itcnt(i*2+3,m) ) goto 300

               end if

            end do

*-----------------------------------------------------------------------
*        dead particle and nucleus except for electron positron
*-----------------------------------------------------------------------

         do j = 1, nclsts

         if( jclusts(4,j) .lt. 0 .and. jclusts(5,j) .ne. 0 .and.
     &       qclusts(7,j) .gt. eletb(1) ) then


                  ipart = jclusts(3,j)
                  kpart = jclusts(7,j)
                  jpart = jclusts(5,j)
                  rpart = qclusts(5,j) * 1000.d0

*-----------------------------------------------------------------------

            call pcheck(m,np,ipart,kpart,jpart,ipn,ips)

               if( ipn .eq. 0 )  goto 285

*-----------------------------------------------------------------------

                  ein = qclusts(7,j)
                  ecc = 0.0d0

*-----------------------------------------------------------------------
*        dedxfnc is specified
*-----------------------------------------------------------------------

         if( itdfun .ne. 0 ) then

               if( ein .lt. eletb(1) ) goto 285

               do i = 1, nlete

                  if( ein .ge. eletb(i) .and.
     &                ein .lt. eletb(i+1) ) goto 35

               end do

   35             ie1 = min( i, nlete )
                  ie2 = 1

*-----------------------------------------------------------------------
*           rrl : track length and rrr : charge particles correction
*-----------------------------------------------------------------------

                  call rainge(ein,rng0,mat,
     &                        ipart,kpart,jpart,rpart)

                  if( rng0 .gt. 1.d+30 ) goto 285

                  tlngth = rng0
                  rng1 = 0.0d0

                  rrl = tlngth
                  rrr = tlngth
                  rrt = 0.0d0

               if( ein .gt. eletb(nlete+1) ) then

                  call rainge(eletb(nlete+1),rngm,mat,
     &                        ipart,kpart,jpart,rpart)

                  rrr = max( 0.0d0, rrr - ( rng0 - rngm ) )
                  rrt = rrt + rng0 -rngm

               end if

*-----------------------------------------------------------------------
*           tally
*-----------------------------------------------------------------------

            do 275 ie = ie1, ie2, -1

                     emie = max( ecc, eletb(ie) )
                     emiu = min( ein, eletb(ie+1) )

                  if( emiu * emie .gt. 0.d0 ) then
                     erg  = sqrt( emiu * emie )
                  else
                     erg  = emiu / 2.d0
                  end if

                  if( mat .gt. 0 ) then

                     call rainge(emie,rng,mat,
     &                           ipart,kpart,jpart,rpart)

                     rng = max( 0.0d0, rng - rng1 )
                     rrl = max( 0.0d0, rrr - rng )

                        rrt  = rrt + rrl

                     if( rrt .gt. tlngth ) then

                        rrl = rrl + ( tlngth - rrt )
                        rrt = tlngth

                     else if( ie .eq. ie2 .and.
     &                        ecc .gt. eletb(ie) ) then

                        rrl = rrl + ( tlngth - rrt )

                     end if

                  end if

*-----------------------------------------------------------------------
*              dedx value multiplyed user factor
*-----------------------------------------------------------------------

               call dedxas(erg,dedx,lmat,ipart,kpart,jpart,rpart)

                  dedx = dedx / 10.d0
                  dtrk = rrl
                  eini = emiu
                  efin = emie

               if( itdfun .eq. 1 ) then

                  call usrdfn1(ipart,kpart,jpart,rpart,
     &                         dedx,dtrk,eini,efin,dhet)

               else if( itdfun .eq. 2 ) then

                  call usrdfn2(ipart,kpart,jpart,rpart,
     &                         dedx,dtrk,eini,efin,dhet)

               end if

*-----------------------------------------------------------------------

                     tlv = dhet ! T.Sato 2018/01/20, delete qclusts(8,j)

                  do ip = 1, ipn

                     trEVENT(ips(ip),i1,it) =
     &               trEVENT(ips(ip),i1,it) + tlv

                  end do

                  if(sweight(m) .ne. 1.d30 .and. sweight(m) .ne. oldwt)
     &              write(*,9000)
                  sweight(m) = oldwt ! T.Sato 2018/03/07  remember weight of this history

                  rrr = rrr - rrl

  275       continue

*-----------------------------------------------------------------------
*        normal case
*-----------------------------------------------------------------------

         else

                     tlv = ( ein - ecc )  ! T.Sato 2018/01/20, delete qclusts(8,j)

                  do ip = 1, ipn

                     trEVENT(ips(ip),i1,it) =
     &               trEVENT(ips(ip),i1,it) + tlv

                  end do

                  if(sweight(m) .ne. 1.d30 .and. sweight(m) .ne. oldwt)
     &              write(*,9000)
                  sweight(m) = oldwt ! T.Sato 2018/03/07  remember weight of this history

         end if

*-----------------------------------------------------------------------

  285       continue

         end if
         end do

*-----------------------------------------------------------------------

  300    continue

      end if

*-----------------------------------------------------------------------

 1000 continue

*-----------------------------------------------------------------------

      return

 9000 Format('[T-deposit] with output = deposit, particle weight is not
     &unique in one history. Result may be unreasonable.')
      end


************************************************************************
*                                                                      *
      subroutine tdepstreg_sumover_ip(m,maxcas, ipall, ratio,
     &                   np,  ne, nr,  nt, tr0x)
*     tr_sum(np,0:ne,nr,nt,2)                                          *
*                                                                      *
************************************************************************

      use TALMOD
!$      use TALMOD0

      implicit double precision (a-h,o-z)
      dimension   tr0x(np,0:ne,nr,nt)

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

          call tdepstreg_sumover_sub_ip(maxcas, ipall, ratio,
     &      itaxs(m,iax),
     &      np,  ne, nr, nt, tr0x,
     &      itpan_sum(m,iax),itenm_sum(m,iax),
     &      itrgn_sum(m,iax),ittnm_sum(m,iax),
     &      tr_sum)

        endif

      enddo

      return
      end


************************************************************************
*                                                                      *
      subroutine tdepstreg_sumover_sub_ip(maxcas, ipall, ratio,
     &           itaxs_in,
     &           np,  ne, nr, nt, tr0,
     &           np_sum,  ne_sum, nr_sum, nt_sum,
     &           tr_sum)
*                                                                      *
************************************************************************

      implicit double precision (a-h,o-z)
      dimension   tr0(np,0:ne,nr,nt)
      dimension   tr_sum(np_sum,0:ne_sum,nr_sum,nt_sum,2)

      if(itaxs_in == 1) then  ! energ
        do it = 1,nt
        do ir = 1,nr
        do ip = 1,np

          tr0_sum = 0.0d0
          tr_sum(ip,0,ir,it,1) = tr_sum(ip,0,ir,it,1) + 
     &                         tr0(ip,0,ir,it) / maxcas
          if(ip.eq.ipall.and.ratio.gt.0) then
            tr_sum(ip,0,ir,it,2) = tr_sum(ip,0,ir,it,2) + 1.0d0
          else
            tr_sum(ip,0,ir,it,2) = tr_sum(ip,0,ir,it,2) +
     &           (tr0(ip,0,ir,it) / maxcas) ** 2
          endif

          tr0_sum = 0.0d0
          do ie = 1,ne
            tr0_sum = tr0_sum + tr0(ip,ie,ir,it) / maxcas
          end do
          tr_sum(ip,1,ir,it,1) =
     &           tr_sum(ip,1,ir,it,1) + tr0_sum
          if(ip.eq.ipall.and.ratio.gt.0) then
            tr_sum(ip,1,ir,it,2) = tr_sum(ip,1,ir,it,2) + 1.0d0
           else
            tr_sum(ip,1,ir,it,2) =
     &           tr_sum(ip,1,ir,it,2) + tr0_sum ** 2
          endif
        end do
        end do
        end do
!

      elseif(itaxs_in == 2) then    ! r
        do it = 1,nt
        do ie = 0,ne
        do ip = 1,np
          tr0_sum = 0.0d0
          do ir = 1,nr
            tr0_sum = tr0_sum + tr0(ip,ie,ir,it) / maxcas
          end do
          tr_sum(ip,ie,1,it,1) =
     &           tr_sum(ip,ie,1,it,1) + tr0_sum
          if(ip.eq.ipall.and.ratio.gt.0) then
            tr_sum(ip,ie,1,it,2) = tr_sum(ip,ie,1,it,2) + 1.0d0
          else
            tr_sum(ip,ie,1,it,2) =
     &           tr_sum(ip,ie,1,it,2) + tr0_sum ** 2
          endif
        end do
        end do
        end do

      elseif(itaxs_in == 11) then   ! time
        do ir = 1,nr
        do ie = 0,ne
        do ip = 1,np
          tr0_sum = 0.0d0
          do it = 1,nt
            tr0_sum = tr0_sum + tr0(ip,ie,ir,it) / maxcas
          end do
          tr_sum(ip,ie,ir,1,1) =
     &           tr_sum(ip,ie,ir,1,1) + tr0_sum
          if(ip.eq.ipall.and.ratio.gt.0) then
            tr_sum(ip,ie,ir,1,2) = tr_sum(ip,ie,ir,1,2) + 1.0d0
          else
            tr_sum(ip,ie,ir,1,2) =
     &           tr_sum(ip,ie,ir,1,2) + tr0_sum ** 2
          endif
        end do
        end do
        end do

       end if
       return
       end


************************************************************************
*                                                                      *
      subroutine tdepstrz_sumover_ip(m,maxcas, ipall, ratio,
     &                   np,  ne, nt, nr,  nz, tr0x)
*     tr_sum(np,0:ne,nt,nr,nz,2)                                          *
*                                                                      *
************************************************************************

      use TALMOD
!$      use TALMOD0

      implicit double precision (a-h,o-z)
      dimension   tr0x(np,0:ne,nt,nr,nz)

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

          call tdepstrz_sumover_sub_ip(maxcas, ipall, ratio,
     &      itaxs(m,iax),
     &      np,  ne, nt, nr, nz, tr0x,
     &      itpan_sum(m,iax),itenm_sum(m,iax),
     &      ittnm_sum(m,iax),itrnm_sum(m,iax),
     &      itznm_sum(m,iax),
     &      tr_sum)

        endif

      enddo

      return
      end


************************************************************************
*                                                                      *
      subroutine tdepstrz_sumover_sub_ip(maxcas, ipall, ratio,
     &           itaxs_in,
     &           np,  ne, nt, nr, nz, tr0,
     &           np_sum,  ne_sum, nt_sum, nr_sum, nz_sum, 
     &           tr_sum)
*                                                                      *
************************************************************************

      implicit double precision (a-h,o-z)
      dimension   tr0(np,0:ne,nt,nr,nz)
      dimension   tr_sum(np_sum,0:ne_sum,nt_sum,nr_sum,nz_sum,2)

      if(itaxs_in == 1) then  ! energ
        do iz = 1,nz
        do ir = 1,nr
        do it = 1,nt
        do ip = 1,np

          tr0_sum = 0.0d0
          tr_sum(ip,0,it,ir,iz,1) = tr_sum(ip,0,it,ir,iz,1) + 
     &                         tr0(ip,0,it,ir,iz) / maxcas
          if(ip.eq.ipall.and.ratio.gt.0) then
            tr_sum(ip,0,it,ir,iz,2) = tr_sum(ip,0,it,ir,iz,2) + 1.0d0
          else
            tr_sum(ip,0,it,ir,iz,2) = tr_sum(ip,0,it,ir,iz,2) +
     &           (tr0(ip,0,it,ir,iz) / maxcas) ** 2
          endif

          tr0_sum = 0.0d0
          do ie = 1,ne
            tr0_sum = tr0_sum + tr0(ip,ie,it,ir,iz) / maxcas
          end do
          tr_sum(ip,1,it,ir,iz,1) =
     &           tr_sum(ip,1,it,ir,iz,1) + tr0_sum
          if(ip.eq.ipall.and.ratio.gt.0) then
            tr_sum(ip,1,it,ir,iz,2) = tr_sum(ip,1,it,ir,iz,2) + 1.0d0
           else
            tr_sum(ip,1,it,ir,iz,2) =
     &           tr_sum(ip,1,it,ir,iz,2) + tr0_sum ** 2
          endif
        end do
        end do
        end do
        end do
!
      elseif(itaxs_in == 5) then    ! z
        do ir = 1,nr
        do it = 1,nt
        do ie = 0,ne
        do ip = 1,np
          tr0_sum = 0.0d0
          do iz = 1,nz
            tr0_sum = tr0_sum + tr0(ip,ie,it,ir,iz) / maxcas
          end do
          tr_sum(ip,ie,it,ir,1,1) =
     &           tr_sum(ip,ie,it,ir,1,1) + tr0_sum
          if(ip.eq.ipall.and.ratio.gt.0) then
            tr_sum(ip,ie,it,ir,1,2) = tr_sum(ip,ie,it,ir,1,2) + 1.0d0
          else
            tr_sum(ip,ie,it,ir,1,2) =
     &           tr_sum(ip,ie,it,ir,1,2) + tr0_sum ** 2
          endif
        end do
        end do
        end do
        end do

      elseif(itaxs_in == 6) then    ! r
        do iz = 1,nz
        do it = 1,nt
        do ie = 0,ne
        do ip = 1,np
          tr0_sum = 0.0d0
          do ir = 1,nr
            tr0_sum = tr0_sum + tr0(ip,ie,it,ir,iz) / maxcas
          end do
          tr_sum(ip,ie,it,1,iz,1) =
     &           tr_sum(ip,ie,it,1,iz,1) + tr0_sum
          if(ip.eq.ipall.and.ratio.gt.0) then
            tr_sum(ip,ie,it,1,iz,2) = tr_sum(ip,ie,it,1,iz,2) + 1.0d0
          else
            tr_sum(ip,ie,it,1,iz,2) =
     &           tr_sum(ip,ie,it,1,iz,2) + tr0_sum ** 2
          endif
        end do
        end do
        end do
        end do

      elseif(itaxs_in == 11) then   ! time
        do iz = 1,nz
        do ir = 1,nr
        do ie = 0,ne
        do ip = 1,np
          tr0_sum = 0.0d0
          do it = 1,nt
            tr0_sum = tr0_sum + tr0(ip,ie,it,ir,iz) / maxcas
          end do
          tr_sum(ip,ie,1,ir,iz,1) =
     &           tr_sum(ip,ie,1,ir,iz,1) + tr0_sum
          if(ip.eq.ipall.and.ratio.gt.0) then
            tr_sum(ip,ie,1,ir,iz,2) = tr_sum(ip,ie,1,ir,iz,2) + 1.0d0
          else
            tr_sum(ip,ie,1,ir,iz,2) =
     &           tr_sum(ip,ie,1,ir,iz,2) + tr0_sum ** 2
          endif
        end do
        end do
        end do
        end do

       end if
       return
       end

************************************************************************
*                                                                      *
      subroutine tdepstxyz_sumover_ip(m,maxcas, ipall, ratio,
     &                   np,  ne, nt, nx, ny, nz, tr0x)
*     tr_sum(np,0:ne,nt,nx,ny,nz,2)                                    *
*                                                                      *
************************************************************************

      use TALMOD
!$      use TALMOD0

      implicit double precision (a-h,o-z)
      dimension   tr0x(np,0:ne,nt,nx*ny*nz)

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

          call tdepstxyz_sumover_sub_ip(maxcas, ipall, ratio,
     &      itaxs(m,iax),
     &      np,  ne, nt, nx, ny, nz, tr0x,
     &      itpan_sum(m,iax),itenm_sum(m,iax),
     &      ittnm_sum(m,iax),itxnm_sum(m,iax),
     &      itynm_sum(m,iax),itznm_sum(m,iax),
     &      tr_sum)

        endif

      enddo

      return
      end

************************************************************************
*                                                                      *
      subroutine tdepstxyz_sumover_sub_ip(maxcas, ipall, ratio,
     &           itaxs_in,
     &           np,  ne, nt, nx, ny, nz, tr0,
     &           np_sum,  ne_sum, nt_sum, nx_sum, ny_sum, nz_sum, 
     &           tr_sum)
*                                                                      *
************************************************************************

      implicit double precision (a-h,o-z)
      dimension   tr0(np,0:ne,nt,nx*ny*nz)
      dimension   tr_sum(np_sum,0:ne_sum,nt_sum,nx_sum*ny_sum*nz_sum,2)

      icf(ix,iy,iz) = ix + ( iy - 1 ) * nx + ( iz - 1 ) * nx * ny
      icf_sum(ix,iy,iz) = ix + ( iy - 1 ) * nx_sum
     &                  + ( iz - 1 ) * nx_sum * ny_sum

      if(itaxs_in == 1) then  ! energ
        do iz = 1,nz
        do iy = 1,ny
        do ix = 1,nx
        do it = 1,nt
        do ip = 1,np

          tr0_sum = 0.0d0
          tr_sum(ip,0,it,icf_sum(ix,iy,iz),1) =
     &        tr_sum(ip,0,it,icf_sum(ix,iy,iz),1) + 
     &                         tr0(ip,0,it,icf(ix,iy,iz)) / maxcas
          if(ip.eq.ipall.and.ratio.gt.0) then
            tr_sum(ip,0,it,icf_sum(ix,iy,iz),2) =
     &                  tr_sum(ip,0,it,icf_sum(ix,iy,iz),2) + 1.0d0
          else
            tr_sum(ip,0,it,icf_sum(ix,iy,iz),2) =
     &            tr_sum(ip,0,it,icf_sum(ix,iy,iz),2) +
     &           (tr0(ip,0,it,icf(ix,iy,iz)) / maxcas) ** 2
          endif

          tr0_sum = 0.0d0
          do ie = 1,ne
            tr0_sum = tr0_sum + tr0(ip,ie,it,icf(ix,iy,iz)) / maxcas
          end do
          tr_sum(ip,1,it,icf_sum(ix,iy,iz),1) =
     &           tr_sum(ip,1,it,icf_sum(ix,iy,iz),1) + tr0_sum
          if(ip.eq.ipall.and.ratio.gt.0) then
            tr_sum(ip,1,it,icf_sum(ix,iy,iz),2) =
     &          tr_sum(ip,1,it,icf_sum(ix,iy,iz),2) + 1.0d0
           else
            tr_sum(ip,1,it,icf_sum(ix,iy,iz),2) =
     &           tr_sum(ip,1,it,icf_sum(ix,iy,iz),2) + tr0_sum ** 2
          endif
        end do
        end do
        end do
        end do
        end do
!
      elseif(itaxs_in == 3) then    ! x
        do iz = 1,nz
        do iy = 1,ny
        do it = 1,nt
        do ie = 0,ne
        do ip = 1,np
          tr0_sum = 0.0d0
          do ix = 1,nx
            tr0_sum = tr0_sum + tr0(ip,ie,it,icf(ix,iy,iz)) / maxcas
          end do
          tr_sum(ip,ie,it,icf_sum(1,iy,iz),1) =
     &           tr_sum(ip,ie,it,icf_sum(1,iy,iz),1) + tr0_sum
          if(ip.eq.ipall.and.ratio.gt.0) then
            tr_sum(ip,ie,it,icf_sum(1,iy,iz),2) =
     &          tr_sum(ip,ie,it,icf_sum(1,iy,iz),2) + 1.0d0
           else
            tr_sum(ip,ie,it,icf_sum(1,iy,iz),2) =
     &           tr_sum(ip,ie,it,icf_sum(1,iy,iz),2) + tr0_sum ** 2
          endif
        end do
        end do
        end do
        end do
        end do

      elseif(itaxs_in == 4) then    ! y
        do iz = 1,nz
        do ix = 1,nx
        do it = 1,nt
        do ie = 0,ne
        do ip = 1,np
          tr0_sum = 0.0d0
          do iy = 1,ny
            tr0_sum = tr0_sum + tr0(ip,ie,it,icf(ix,iy,iz)) / maxcas
          end do
          tr_sum(ip,ie,it,icf_sum(ix,1,iz),1) =
     &           tr_sum(ip,ie,it,icf_sum(ix,1,iz),1) + tr0_sum
          if(ip.eq.ipall.and.ratio.gt.0) then
            tr_sum(ip,ie,it,icf_sum(ix,1,iz),2) =
     &          tr_sum(ip,ie,it,icf_sum(ix,1,iz),2) + 1.0d0
           else
            tr_sum(ip,ie,it,icf_sum(ix,1,iz),2) =
     &           tr_sum(ip,ie,it,icf_sum(ix,1,iz),2) + tr0_sum ** 2
          endif
        end do
        end do
        end do
        end do
        end do

      elseif(itaxs_in == 5) then    ! z
        do iy = 1,ny
        do ix = 1,nx
        do it = 1,nt
        do ie = 0,ne
        do ip = 1,np
          tr0_sum = 0.0d0
          do iz = 1,nz
            tr0_sum = tr0_sum + tr0(ip,ie,it,icf(ix,iy,iz)) / maxcas
          end do
          tr_sum(ip,ie,it,icf_sum(ix,iy,1),1) =
     &           tr_sum(ip,ie,it,icf_sum(ix,iy,1),1) + tr0_sum
          if(ip.eq.ipall.and.ratio.gt.0) then
            tr_sum(ip,ie,it,icf_sum(ix,iy,1),2) =
     &          tr_sum(ip,ie,it,icf_sum(ix,iy,1),2) + 1.0d0
           else
            tr_sum(ip,ie,it,icf_sum(ix,iy,1),2) =
     &           tr_sum(ip,ie,it,icf_sum(ix,iy,1),2) + tr0_sum ** 2
          endif
        end do
        end do
        end do
        end do
        end do

      elseif(itaxs_in == 11) then    ! t
        do iz = 1,nz
        do iy = 1,ny
        do ix = 1,nx
        do ie = 0,ne
        do ip = 1,np
          tr0_sum = 0.0d0
          do it = 1,nt
            tr0_sum = tr0_sum + tr0(ip,ie,it,icf(ix,iy,iz)) / maxcas
          end do
          tr_sum(ip,ie,1,icf_sum(ix,iy,iz),1) =
     &           tr_sum(ip,ie,1,icf_sum(ix,iy,iz),1) + tr0_sum
          if(ip.eq.ipall.and.ratio.gt.0) then
            tr_sum(ip,ie,1,icf_sum(ix,iy,iz),2) =
     &          tr_sum(ip,ie,1,icf_sum(ix,iy,iz),2) + 1.0d0
           else
            tr_sum(ip,ie,1,icf_sum(ix,iy,iz),2) =
     &           tr_sum(ip,ie,1,icf_sum(ix,iy,iz),2) + tr0_sum ** 2
          endif
        end do
        end do
        end do
        end do
        end do

       end if
       return
       end
