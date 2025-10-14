************************************************************************
*                                                                      *
      subroutine update(ic)
*                                                                      *
*       update particle data                                           *
*       modified by K.Niita on 2002/02/06                              *
*                                                                      *
*                                                                      *
************************************************************************
      use MMBANKMOD !FURUTA
      implicit real*8 (a-h,o-z)

cKN 2015/03/19 nsos, nsosa, ibksos, iaksos
*-----------------------------------------------------------------------

      include 'param.inc'

*-----------------------------------------------------------------------

      common /bnkmem/ maxbnk, maxbn2, rtrckflp
      common /icomon/ no,mat,ityp,ktyp,jtyp,mtyp,rtyp,ctyp
!$OMP THREADPRIVATE(/icomon/)
      common /jcomon/ nabov,nobch,nocas,nomax
!$OMP THREADPRIVATE(/jcomon/)
      common /bnkmm/  mbmax, mbfin, mbtfin
      common /kcomon/ mbnk, ibnk, jbnk
!$OMP THREADPRIVATE(/kcomon/)
      common /kcomsi/ iibnk, jjbnk !FURUTA
!$OMP THREADPRIVATE(/kcomsi/)

*-----------------------------------------------------------------------

      dimension idas(1)
      equivalence ( das, idas )

      integer ifile                    !FURUTA
      integer,save:: iflsave(15)       !FURUTA
!$OMP THREADPRIVATE(iflsave)
      logical,save:: flock(15)=.false. !FURUTA

      integer,save:: iwarn=0 !FURUTA20151009

      integer,save::   name_leadsave , nty_leadsave  , nkf_leadsave  ,
     & iblz_leadsave , nmed_leadsave , nfcs_leadsave , nzst_leadsave ,
     & nsos_leadsave , ncnt1_leadsave, ncnt2_leadsave, ncnt3_leadsave,
     & itetpos_leadsave

      real*8, save::   e_leadsave    , t_leadsave    , x_leadsave    ,
     & y_leadsave    , z_leadsave    , u_leadsave    , v_leadsave    ,
     & w_leadsave    , wt_leadsave   , wtin_leadsave , wtnz_leadsave ,
     & xfcs_leadsave , spx_leadsave  , spy_leadsave  , spz_leadsave

!$OMP THREADPRIVATE(name_leadsave , nty_leadsave  , nkf_leadsave  )
!$OMP THREADPRIVATE(iblz_leadsave , nmed_leadsave , nfcs_leadsave )
!$OMP THREADPRIVATE(nzst_leadsave , nsos_leadsave , ncnt1_leadsave)
!$OMP THREADPRIVATE(ncnt2_leadsave, ncnt3_leadsave, itetpos_leadsave)
!$OMP THREADPRIVATE(e_leadsave    , t_leadsave    , x_leadsave    )
!$OMP THREADPRIVATE(y_leadsave    , z_leadsave    , u_leadsave    )
!$OMP THREADPRIVATE(v_leadsave    , w_leadsave    , wt_leadsave   )
!$OMP THREADPRIVATE(wtin_leadsave , wtnz_leadsave , xfcs_leadsave )
!$OMP THREADPRIVATE(spx_leadsave  , spy_leadsave  , spz_leadsave  )

      common /wwindp0/ wupn0, wsurvn0, mxspln0, mwhere0, mvoww0 ! T.Sato 2024/05/26
!$OMP THREADPRIVATE(/wwindp0/)

*-----------------------------------------------------------------------
*     ic = 1 : check remaining data
*     ic = 2 : normal update
*-----------------------------------------------------------------------

      Select Case(ic)
         Case(1)
*-----------------------------------------------------------------------
*           read from remaining memory
*-----------------------------------------------------------------------

            if( jbnk .gt. 0 ) then

              write(*,*)'ERROR in UPDATE' !FURUTA
              STOP                        !FURUTA

*-----------------------------------------------------------------------
*           read from temporary files
*-----------------------------------------------------------------------

            else

               if( ibnk .le. 0 ) then

                  no   = 0
                  mbnk = 0

                  write(6,'(/"**** Warning: ",
     &                       "something wrong in update ****")')

                  return

               end if

                  iop = iflsave(ibnk) + 39                              !FURUTA
                  if(iop.lt.40.or.iop.gt.54)then                        !FURUTA
                    no   = 0                                            !FURUTA
                    mbnk = 0                                            !FURUTA
                    write(6,'(/"**** Warning: ",
     &                   "something wrong in update ver.FURUTA****")')!FURUTA
                    return                                              !FURUTA
                  endif                                                 !FURUTA

                  rewind iop
                  read(iop) k

               do m = 1, k

                  read(iop)
     &            iblz(ibkblz+m,ipomp+1), name(ibknam+m,ipomp+1),
     &            nmed(ibknmd+m,ipomp+1),
     &            nty(ibknty+m,ipomp+1),  nkf(ibknkf+m,ipomp+1),
     &            nfcs(ibknfc+m,ipomp+1),
     &            ncnt(ibknct+1,m,ipomp+1),ncnt(ibknct+2,m,ipomp+1),
     &            ncnt(ibknct+3,m,ipomp+1),
     &            nzst(ibkzst+m,ipomp+1),nsos(ibksos+m,ipomp+1),
     &            itetpos(ibtetpos+m,ipomp+1) !FURUTA20170315

               end do

               do m = 1, k

                  read(iop)
     &            wt(ibkwt+m,ipomp+1), u(ibku+m,ipomp+1),
     &            v(ibkv+m,ipomp+1), w(ibkw+m,ipomp+1),
     &            e(ibke+m,ipomp+1),   x(ibkx+m,ipomp+1),
     &            y(ibky+m,ipomp+1), z(ibkz+m,ipomp+1),
     &            t(ibkt+m,ipomp+1),   wtin(ibkwin+m,ipomp+1),
     &            wtnz(ibkwnz+m,ipomp+1), xfcs(ibkxfc+m,ipomp+1),
     &            spx(ibkspx+m,ipomp+1),spy(ibkspy+m,ipomp+1),
     &            spz(ibkspz+m,ipomp+1)

               end do

                  close(iop)
!$OMP CRITICAL (flock_crit)
                  flock(iflsave(ibnk))=.false. !FURUTA
!$OMP END CRITICAL (flock_crit)
                  ibnk = ibnk - 1

!                  write(*,'("jbnk =",i3,", ibnk =",i3 )') jbnk, ibnk

                  if( ibnk .eq. 0 ) mbnk = 0

            end if

*-----------------------------------------------------------------------

               nomax = k
               no    = 1

       return

*-----------------------------------------------------------------------
*     ic = 2 : normal update
*-----------------------------------------------------------------------

         Case(2)

*-----------------------------------------------------------------------
*           nabov = 0
*-----------------------------------------------------------------------

            if( nabov .le. 0 ) return

*-----------------------------------------------------------------------
*     nabov > maxbnk - nomax
*-----------------------------------------------------------------------

      if( nomax + nabov .gt. maxbnk ) then

*-----------------------------------------------------------------------
*        delete information if number exceeds maxbnk
*-----------------------------------------------------------------------

         if( no .gt. maxbn2 ) then

                  k = nomax - no + 1

               do ll = 1, k

                  m = no + ll - 1

                  e(ibke+ll,ipomp+1)      = e(ibke+m,ipomp+1)
                  t(ibkt+ll,ipomp+1)      = t(ibkt+m,ipomp+1)
                  name(ibknam+ll,ipomp+1) = name(ibknam+m,ipomp+1)
                  nty(ibknty+ll,ipomp+1)  = nty(ibknty+m,ipomp+1)
                  nkf(ibknkf+ll,ipomp+1)  = nkf(ibknkf+m,ipomp+1)

                  x(ibkx+ll,ipomp+1) = x(ibkx+m,ipomp+1)
                  y(ibky+ll,ipomp+1) = y(ibky+m,ipomp+1)
                  z(ibkz+ll,ipomp+1) = z(ibkz+m,ipomp+1)

                  u(ibku+ll,ipomp+1) = u(ibku+m,ipomp+1)
                  v(ibkv+ll,ipomp+1) = v(ibkv+m,ipomp+1)
                  w(ibkw+ll,ipomp+1) = w(ibkw+m,ipomp+1)

                  wt(ibkwt+ll,ipomp+1)      = wt(ibkwt+m,ipomp+1)
                  iblz(ibkblz+ll,ipomp+1)   = iblz(ibkblz+m,ipomp+1)
                  nmed(ibknmd+ll,ipomp+1)   = nmed(ibknmd+m,ipomp+1)
                  wtin(ibkwin+ll,ipomp+1)   = wtin(ibkwin+m,ipomp+1)
                  wtnz(ibkwnz+ll,ipomp+1)   = wtnz(ibkwnz+m,ipomp+1)
                  nfcs(ibknfc+ll,ipomp+1)   = nfcs(ibknfc+m,ipomp+1)
                  xfcs(ibkxfc+ll,ipomp+1)   = xfcs(ibkxfc+m,ipomp+1)

                  spx(ibkspx+ll,ipomp+1)    = spx(ibkspx+m,ipomp+1)
                  spy(ibkspy+ll,ipomp+1)    = spy(ibkspy+m,ipomp+1)
                  spz(ibkspz+ll,ipomp+1)    = spz(ibkspz+m,ipomp+1)

                  nzst(ibkzst+ll,ipomp+1)   = nzst(ibkzst+m,ipomp+1)
                  nsos(ibksos+ll,ipomp+1)   = nsos(ibksos+m,ipomp+1)

                  ncnt(ibknct+1,ll,ipomp+1) = ncnt(ibknct+1,m,ipomp+1)
                  ncnt(ibknct+2,ll,ipomp+1) = ncnt(ibknct+2,m,ipomp+1)
                  ncnt(ibknct+3,ll,ipomp+1) = ncnt(ibknct+3,m,ipomp+1)
                itetpos(ibtetpos+ll,ipomp+1)=itetpos(ibtetpos+m,ipomp+1)

               end do

                  nomax = k
                  no    = 1

*-----------------------------------------------------------------------
*        write the information on the temporary bank or file
*-----------------------------------------------------------------------

         else

                  mbnk = 1

                  k = min( nomax - no, maxbn2 )

*-----------------------------------------------------------------------
*           write on remaining memory
*-----------------------------------------------------------------------
cKN:byFuruta,2011/02/21(K*20+1)
cFURUTA memory buffer move to MMBANKMOD

                  ibnk = ibnk + 1
                  iibnk = iibnk + 1
!$OMP CRITICAL (flock_crit)
                  do ifile=1,15                !FURUTA
                    if(.not.flock(ifile))exit  !FURUTA
                  enddo                        !FURUTA
                  flock(ifile)=.true.          !FURUTA
!$OMP END CRITICAL (flock_crit)

               if( ifile .gt. 15 ) then !FURUTA

                  write(6,'(/"**** Fatal ERROR: lack of bank",
     &                       " space ****"/
     &                       "     please increase maxbnk")')
                  call parastop( 837 )

               end if

                  iop = ifile + 39      !FURUTA
                  iflsave(ibnk) = ifile !FURUTA

                  open(iop,form='unformatted',status='scratch')

                  write(iop) k

               do ll = 1, k

                  m = no + ll

                  write(iop)
     &            iblz(ibkblz+m,ipomp+1), name(ibknam+m,ipomp+1),
     &            nmed(ibknmd+m,ipomp+1),
     &            nty(ibknty+m,ipomp+1),  nkf(ibknkf+m,ipomp+1),
     &            nfcs(ibknfc+m,ipomp+1),
     &            ncnt(ibknct+1,m,ipomp+1),ncnt(ibknct+2,m,ipomp+1),
     &            ncnt(ibknct+3,m,ipomp+1),
     &            nzst(ibkzst+m,ipomp+1),nsos(ibksos+m,ipomp+1),
     &            itetpos(ibtetpos+m,ipomp+1) !FURUTA20170315

               end do

               do ll = 1, k

                  m = no + ll

                  write(iop)
     &            wt(ibkwt+m,ipomp+1), u(ibku+m,ipomp+1),
     &            v(ibkv+m,ipomp+1), w(ibkw+m,ipomp+1),
     &            e(ibke+m,ipomp+1),   x(ibkx+m,ipomp+1),
     &            y(ibky+m,ipomp+1), z(ibkz+m,ipomp+1),
     &            t(ibkt+m,ipomp+1),   wtin(ibkwin+m,ipomp+1),
     &            wtnz(ibkwnz+m,ipomp+1), xfcs(ibkxfc+m,ipomp+1),
     &            spx(ibkspx+m,ipomp+1),spy(ibkspy+m,ipomp+1),
     &            spz(ibkspz+m,ipomp+1)

               end do

cFURUTA20151009-----
               iwarn=iwarn+1
               if(iwarn.eq.1.or.mod(iwarn,10).eq.0)then
                  write(*,'(/"**** Warning: Too many secondary",
     &                       " particles created ****")')
                  write(*,'("**** MAXBNK overflowed"
     &                 " thus HDD is used ",i5," times****")')iwarn
           write(*,'(" If you feel calculation is too slow,",
     &    " please increase MAXBNK or reconsider the setup"/)')
                 wupn0=1.0e10  ! no more particle split allowed after this warning, T.Sato 2024/5/26
               endif
c-------------------

!                  write(*,'("jbnk =",i3,", ibnk =",i3 )') jbnk, ibnk
*-----------------------------------------------------------------------

                  m  = no
                  ll = 1

                  e(ibke+ll,ipomp+1)      = e(ibke+m,ipomp+1)
                  t(ibkt+ll,ipomp+1)      = t(ibkt+m,ipomp+1)
                  name(ibknam+ll,ipomp+1) = name(ibknam+m,ipomp+1)
                  nty(ibknty+ll,ipomp+1)  = nty(ibknty+m,ipomp+1)
                  nkf(ibknkf+ll,ipomp+1)  = nkf(ibknkf+m,ipomp+1)

                  x(ibkx+ll,ipomp+1) = x(ibkx+m,ipomp+1)
                  y(ibky+ll,ipomp+1) = y(ibky+m,ipomp+1)
                  z(ibkz+ll,ipomp+1) = z(ibkz+m,ipomp+1)

                  u(ibku+ll,ipomp+1) = u(ibku+m,ipomp+1)
                  v(ibkv+ll,ipomp+1) = v(ibkv+m,ipomp+1)
                  w(ibkw+ll,ipomp+1) = w(ibkw+m,ipomp+1)

                  wt(ibkwt+ll,ipomp+1)      =   wt(ibkwt+m,ipomp+1)
                  iblz(ibkblz+ll,ipomp+1)   = iblz(ibkblz+m,ipomp+1)
                  nmed(ibknmd+ll,ipomp+1)   = nmed(ibknmd+m,ipomp+1)
                  wtin(ibkwin+ll,ipomp+1)   = wtin(ibkwin+m,ipomp+1)
                  wtnz(ibkwnz+ll,ipomp+1)   = wtnz(ibkwnz+m,ipomp+1)
                  nfcs(ibknfc+ll,ipomp+1)   = nfcs(ibknfc+m,ipomp+1)
                  xfcs(ibkxfc+ll,ipomp+1)   = xfcs(ibkxfc+m,ipomp+1)

                  spx(ibkspx+ll,ipomp+1)    = spx(ibkspx+m,ipomp+1)
                  spy(ibkspy+ll,ipomp+1)    = spy(ibkspy+m,ipomp+1)
                  spz(ibkspz+ll,ipomp+1)    = spz(ibkspz+m,ipomp+1)

                  nzst(ibkzst+ll,ipomp+1)   = nzst(ibkzst+m,ipomp+1)
                  nsos(ibksos+ll,ipomp+1)   = nsos(ibksos+m,ipomp+1)

                  ncnt(ibknct+1,ll,ipomp+1) = ncnt(ibknct+1,m,ipomp+1)
                  ncnt(ibknct+2,ll,ipomp+1) = ncnt(ibknct+2,m,ipomp+1)
                  ncnt(ibknct+3,ll,ipomp+1) = ncnt(ibknct+3,m,ipomp+1)
                itetpos(ibtetpos+ll,ipomp+1)=itetpos(ibtetpos+m,ipomp+1)

*-----------------------------------------------------------------------

                  ki = no + k
                  k  = nomax - ki

            if( k .gt. 0 ) then

               do kk = 1, k

                  m  = ki + kk
                  ll = kk + 1

                  e(ibke+ll,ipomp+1)      = e(ibke+m,ipomp+1)
                  t(ibkt+ll,ipomp+1)      = t(ibkt+m,ipomp+1)
                  name(ibknam+ll,ipomp+1) = name(ibknam+m,ipomp+1)
                  nty(ibknty+ll,ipomp+1)  = nty(ibknty+m,ipomp+1)
                  nkf(ibknkf+ll,ipomp+1)  = nkf(ibknkf+m,ipomp+1)

                  x(ibkx+ll,ipomp+1) = x(ibkx+m,ipomp+1)
                  y(ibky+ll,ipomp+1) = y(ibky+m,ipomp+1)
                  z(ibkz+ll,ipomp+1) = z(ibkz+m,ipomp+1)

                  u(ibku+ll,ipomp+1) = u(ibku+m,ipomp+1)
                  v(ibkv+ll,ipomp+1) = v(ibkv+m,ipomp+1)
                  w(ibkw+ll,ipomp+1) = w(ibkw+m,ipomp+1)

                  wt(ibkwt+ll,ipomp+1)      = wt(ibkwt+m,ipomp+1)
                  iblz(ibkblz+ll,ipomp+1)   = iblz(ibkblz+m,ipomp+1)
                  nmed(ibknmd+ll,ipomp+1)   = nmed(ibknmd+m,ipomp+1)
                  wtin(ibkwin+ll,ipomp+1)   = wtin(ibkwin+m,ipomp+1)
                  wtnz(ibkwnz+ll,ipomp+1)   = wtnz(ibkwnz+m,ipomp+1)
                  nfcs(ibknfc+ll,ipomp+1)   = nfcs(ibknfc+m,ipomp+1)
                  xfcs(ibkxfc+ll,ipomp+1)   = xfcs(ibkxfc+m,ipomp+1)

                  spx(ibkspx+ll,ipomp+1)    = spx(ibkspx+m,ipomp+1)
                  spy(ibkspy+ll,ipomp+1)    = spy(ibkspy+m,ipomp+1)
                  spz(ibkspz+ll,ipomp+1)    = spz(ibkspz+m,ipomp+1)

                  nzst(ibkzst+ll,ipomp+1)   = nzst(ibkzst+m,ipomp+1)
                  nsos(ibksos+ll,ipomp+1)   = nsos(ibksos+m,ipomp+1)

                  ncnt(ibknct+1,ll,ipomp+1) = ncnt(ibknct+1,m,ipomp+1)
                  ncnt(ibknct+2,ll,ipomp+1) = ncnt(ibknct+2,m,ipomp+1)
                  ncnt(ibknct+3,ll,ipomp+1) = ncnt(ibknct+3,m,ipomp+1)
               itetpos(ibtetpos+ll,ipomp+1)=itetpos(ibtetpos+m,ipomp+1)

               end do

                  nomax = k + 1
                  no    = 1

            else

                  nomax = 1
                  no    = 1

            end if

*-----------------------------------------------------------------------

         end if

      end if

*-----------------------------------------------------------------------
*     add the particles
*-----------------------------------------------------------------------

         do ll = 1, nabov

               m = nomax + ll

               e(ibke+m,ipomp+1)      = ea(iake+ll,ipomp+1)
               t(ibkt+m,ipomp+1)      = ta(iakt+ll,ipomp+1)
               name(ibknam+m,ipomp+1) = namea(iaknam+ll,ipomp+1)
               nty(ibknty+m,ipomp+1)  = ntya(iaknty+ll,ipomp+1)
               nkf(ibknkf+m,ipomp+1)  = nkfa(iaknkf+ll,ipomp+1)

               u(ibku+m,ipomp+1)    = ua(iaku+ll,ipomp+1)
               v(ibkv+m,ipomp+1)    = va(iakv+ll,ipomp+1)
               w(ibkw+m,ipomp+1)    = wa(iakw+ll,ipomp+1)
               wt(ibkwt+m,ipomp+1)  = wta(iakwt+ll,ipomp+1)

               x(ibkx+m,ipomp+1)    = xa(iakx+ll,ipomp+1)
               y(ibky+m,ipomp+1)    = ya(iaky+ll,ipomp+1)
               z(ibkz+m,ipomp+1)    = za(iakz+ll,ipomp+1)

               iblz(ibkblz+m,ipomp+1)   = iblza(iakblz+ll,ipomp+1)
               nmed(ibknmd+m,ipomp+1)   = nmeda(iaknmd+ll,ipomp+1)
               wtin(ibkwin+m,ipomp+1)   = wtina(iakwin+ll,ipomp+1)
               wtnz(ibkwnz+m,ipomp+1)   = wtnza(iakwnz+ll,ipomp+1)
               nfcs(ibknfc+m,ipomp+1)   = nfcsa(iaknfc+ll,ipomp+1)
               xfcs(ibkxfc+m,ipomp+1)   = xfcsa(iakxfc+ll,ipomp+1)

               spx(ibkspx+m,ipomp+1)    = spxa(iakspx+ll,ipomp+1)
               spy(ibkspy+m,ipomp+1)    = spya(iakspy+ll,ipomp+1)
               spz(ibkspz+m,ipomp+1)    = spza(iakspz+ll,ipomp+1)

               nzst(ibkzst+m,ipomp+1)   = nzsta(iakzst+ll,ipomp+1)
               nsos(ibksos+m,ipomp+1)   = nsosa(iaksos+ll,ipomp+1)

               ncnt(ibknct+1,m,ipomp+1) = ncnta(iaknct+1,ll,ipomp+1)
               ncnt(ibknct+2,m,ipomp+1) = ncnta(iaknct+2,ll,ipomp+1)
               ncnt(ibknct+3,m,ipomp+1) = ncnta(iaknct+3,ll,ipomp+1)
              itetpos(ibtetpos+m,ipomp+1)=itetposa(ibtetposa+ll,ipomp+1)

         end do

               nomax = m
               nabov = 0


*-----------------------------------------------------------------------

       return

*-----------------------------------------------------------------------
*     ic = 3 : recover leading particle
*-----------------------------------------------------------------------

         Case(3)

                no    = 1

                e(ibke+no,ipomp+1)           = e_leadsave
                t(ibkt+no,ipomp+1)           = t_leadsave
                name(ibknam+no,ipomp+1)      = name_leadsave
                nty(ibknty+no,ipomp+1)       = nty_leadsave
                nkf(ibknkf+no,ipomp+1)       = nkf_leadsave

                x(ibkx+no,ipomp+1)           = x_leadsave
                y(ibky+no,ipomp+1)           = y_leadsave
                z(ibkz+no,ipomp+1)           = z_leadsave

                u(ibku+no,ipomp+1)           = u_leadsave
                v(ibkv+no,ipomp+1)           = v_leadsave
                w(ibkw+no,ipomp+1)           = w_leadsave

                wt(ibkwt+no,ipomp+1)         = wt_leadsave
                iblz(ibkblz+no,ipomp+1)      = iblz_leadsave
                nmed(ibknmd+no,ipomp+1)      = nmed_leadsave
                wtin(ibkwin+no,ipomp+1)      = wtin_leadsave
                wtnz(ibkwnz+no,ipomp+1)      = wtnz_leadsave
                nfcs(ibknfc+no,ipomp+1)      = nfcs_leadsave
                xfcs(ibkxfc+no,ipomp+1)      = xfcs_leadsave

                spx(ibkspx+no,ipomp+1)       = spx_leadsave
                spy(ibkspy+no,ipomp+1)       = spy_leadsave
                spz(ibkspz+no,ipomp+1)       = spz_leadsave

                nzst(ibkzst+no,ipomp+1)      = nzst_leadsave
                nsos(ibksos+no,ipomp+1)      = nsos_leadsave

                ncnt(ibknct+1,no,ipomp+1)    = ncnt1_leadsave
                ncnt(ibknct+2,no,ipomp+1)    = ncnt2_leadsave
                ncnt(ibknct+3,no,ipomp+1)    = ncnt3_leadsave
                itetpos(ibtetpos+no,ipomp+1) = itetpos_leadsave

                no    = 0
                nomax = 1


       return

*-----------------------------------------------------------------------
*     ic = 4 : store leading particle
*-----------------------------------------------------------------------

         Case(4)

                  e_leadsave        = e(ibke+no,ipomp+1)
                  t_leadsave        = t(ibkt+no,ipomp+1)
                  name_leadsave     = name(ibknam+no,ipomp+1)
                  nty_leadsave      = nty(ibknty+no,ipomp+1)
                  nkf_leadsave      = nkf(ibknkf+no,ipomp+1)

                  x_leadsave        = x(ibkx+no,ipomp+1)
                  y_leadsave        = y(ibky+no,ipomp+1)
                  z_leadsave        = z(ibkz+no,ipomp+1)

                  u_leadsave        = u(ibku+no,ipomp+1)
                  v_leadsave        = v(ibkv+no,ipomp+1)
                  w_leadsave        = w(ibkw+no,ipomp+1)

                  wt_leadsave       = wt(ibkwt+no,ipomp+1)
                  iblz_leadsave     = iblz(ibkblz+no,ipomp+1)
                  nmed_leadsave     = nmed(ibknmd+no,ipomp+1)
                  wtin_leadsave     = wtin(ibkwin+no,ipomp+1)
                  wtnz_leadsave     = wtnz(ibkwnz+no,ipomp+1)
                  nfcs_leadsave     = nfcs(ibknfc+no,ipomp+1)
                  xfcs_leadsave     = xfcs(ibkxfc+no,ipomp+1)

                  spx_leadsave      = spx(ibkspx+no,ipomp+1)
                  spy_leadsave      = spy(ibkspy+no,ipomp+1)
                  spz_leadsave      = spz(ibkspz+no,ipomp+1)

                  nzst_leadsave     = nzst(ibkzst+no,ipomp+1)
                  nsos_leadsave     = nsos(ibksos+no,ipomp+1)

                  ncnt1_leadsave   = ncnt(ibknct+1,no,ipomp+1)
                  ncnt2_leadsave   = ncnt(ibknct+2,no,ipomp+1)
                  ncnt3_leadsave   = ncnt(ibknct+3,no,ipomp+1)
                  itetpos_leadsave = itetpos(ibtetpos+no,ipomp+1)

       return
*-----------------------------------------------------------------------

      End Select

      end

