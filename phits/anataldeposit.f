!***********************************************************************
!                                                                      *
      subroutine anatal_calc_anova_depstreg(m,ntf,
     &                             np,nei,ne,nt,nr,
     &                             sumfactor,weightRate,
     &                             trRES,tranatal,nfile,
     &                             resc2SUMTAL,resc3SUMTAL,
     &                             ierr)
!                                                                      *
!     setting of tranatal of t-deposit with reg mesh                   *
!     created by S.Hashimoto on 2015/1/20                              *
!                                                                      *
!***********************************************************************
      implicit none
!-----------------------------------------------------------------------
      include 'param.inc'
      include 'err.inc'
!-----------------------------------------------------------------------
      integer   istdevres
      integer   maxcasres
      real*8    rijklstres
      integer   irdrf
      common /res01/ istdevres,maxcasres,rijklstres,irdrf

      double precision resc2   ! W
      double precision resc3   ! N
      common /restart/ resc2(itlmax), resc3(itlmax)
!-----------------------------------------------------------------------
      character m_err*200
      integer   l_err, k_err
      common /error/ m_err, l_err, k_err
!-----------------------------------------------------------------------
      integer           m
      integer           ntf
      integer           np
      integer           ne
      integer           nt
      integer           nr
      integer           nei

      double precision  sumfactor     ! normalization factor
      double precision  weightRate(*) ! weighting rate
      double precision  trRES   (np,0:ne,nr,nt,2)
      integer           nfile
      double precision  tranatal(np,0:ne,nr,nt,2*nfile)
      double precision  resc2SUMTAL
      double precision  resc3SUMTAL

      double precision  cmin, cmax

      integer           ierr
!-----------------------------------------------------------------------
      integer           ir

      integer           ip
      integer           ie
      integer           it
      double precision  Xa, sigx

      integer   nr0
      integer   itrwgtsum
      integer   itrgn1
      integer   itrgm1
      integer   itrncd
      common /tall70/ itrwgtsum(itlmax), itrgn1(itlmax),
     &                itrgm1(itlmax), itrncd(itlmax)

!-----------------------------------------------------------------------
      integer           iat, iad, iaf

      iat(iad,iaf) = iad + (iaf-1) * 2
!-----------------------------------------------------------------------

      ierr = 0

      if( itrwgtsum(m) .eq. 1 ) then
       nr0 = 1
      else
       nr0 = nr
      endif

*-----------------------------------------------------------------------
      cmax = 0.0
      cmin = 1.e+33

      do ir = 1, nr0
       do it = 1, nt
        do ie = nei, ne
         do ip = 1, np

          if( trRES(ip,ie,ir,it,1) .gt. 0.d0 ) then

           call calc_deposit_stdev(m,Xa,sigx,
     &                             trRES(ip,ie,ir,it,1),
     &                             trRES(ip,ie,ir,it,2),
     &                             1.0d+0,ip)
           trRES(ip,ie,ir,it,1) = Xa
           trRES(ip,ie,ir,it,2) = sigx
          end if

          ! X_bar = Xj_bar
          tranatal(ip,ie,ir,it,iat(1,ntf)) =
     &       trRES(ip,ie,ir,it,1)

          ! (sig_xj)**2
          tranatal(ip,ie,ir,it,iat(2,ntf)) =
     &    (  trRES(ip,ie,ir,it,2)
     &     * trRES(ip,ie,ir,it,1) )**2

          ! sig_x
          tranatal(ip,ie,ir,it,iat(2,ntf)) = dsqrt(
     &    tranatal(ip,ie,ir,it,iat(2,ntf)))

          ! Sigma(xi wi)**2 = (sig_X**2 N(N-1)+N X_bar**2)(W/N)**2
          tranatal(ip,ie,ir,it,iat(2,ntf)) =
     &   (tranatal(ip,ie,ir,it,iat(2,ntf))**2 *
     &        resc3(m) * (resc3(m)-1.0d0) +
     &        resc3(m) *
     &    tranatal(ip,ie,ir,it,iat(1,ntf))**2) *
     &       (resc2(m)/resc3(m))**2

          ! Sigma xi wi = X_bar W
          tranatal(ip,ie,ir,it,iat(1,ntf)) =
     &    tranatal(ip,ie,ir,it,iat(1,ntf)) * resc2(m)

          if( tranatal(ip,ie,ir,it,iat(1,ntf)) > cmax )
     & cmax = tranatal(ip,ie,ir,it,iat(1,ntf))
          if( tranatal(ip,ie,ir,it,iat(1,ntf)) < cmin )
     & cmin = tranatal(ip,ie,ir,it,iat(1,ntf))

         end do   ! ip loop end
        end do    ! ie loop end
       end do     ! it loop end
      end do      ! ir loop end

! sumover
      call anatal_calc_deposit_anova_sub(m,ntf,np,
     &                               resc2(m),resc3(m) )

!-----------------------------------------------------------------------
      return

      end subroutine anatal_calc_anova_depstreg

!***********************************************************************
!                                                                      *
      subroutine anatal_calc_anova_depstrz(m,ntf,
     &                             np,nei,ne,nt,nr,nz,
     &                             sumfactor,weightRate,
     &                             trRES,tranatal,nfile,
     &                             resc2SUMTAL,resc3SUMTAL,
     &                             ierr)
!                                                                      *
!     setting of tranatal of t-deposit with r-z mesh                   *
!     created by S.Hashimoto on 2015/1/20                              *
!                                                                      *
!***********************************************************************
      implicit none
!-----------------------------------------------------------------------
      include 'param.inc'
      include 'err.inc'
!-----------------------------------------------------------------------
      integer   istdevres
      integer   maxcasres
      real*8    rijklstres
      integer   irdrf
      common /res01/ istdevres,maxcasres,rijklstres,irdrf

      double precision resc2   ! W
      double precision resc3   ! N
      common /restart/ resc2(itlmax), resc3(itlmax)
!-----------------------------------------------------------------------
      character m_err*200
      integer   l_err, k_err
      common /error/ m_err, l_err, k_err
!-----------------------------------------------------------------------
      integer           m
      integer           ntf
      integer           np
      integer           ne
      integer           nt
      integer           nr
      integer           nz
      integer           nei

      double precision  sumfactor     ! normalization factor
      double precision  weightRate(*) ! weighting rate
      double precision  trRES   (np,0:ne,nt,nr,nz,2)
      integer           nfile
      double precision  tranatal(np,0:ne,nt,nr,nz,2*nfile)

      double precision  resc2SUMTAL
      double precision  resc3SUMTAL

      double precision  cmin, cmax

      integer           ierr
!-----------------------------------------------------------------------
      integer           ir
      integer           iz

      integer           ip
      integer           ie
      integer           it
      double precision  Xa, sigx
*-----------------------------------------------------------------------
      integer           iat, iad, iaf

      iat(iad,iaf) = iad + (iaf-1) * 2
!-----------------------------------------------------------------------

      ierr = 0

*-----------------------------------------------------------------------
      cmax = 0.0
      cmin = 1.e+33

      do iz = 1, nz
       do ir = 1, nr
        do it = 1, nt
         do ie = nei, ne
          do ip = 1, np

           if( trRES(ip,ie,it,ir,iz,1) .gt. 0.d0 ) then
            call calc_deposit_stdev(m,Xa,sigx,
     &                              trRES(ip,ie,it,ir,iz,1),
     &                              trRES(ip,ie,it,ir,iz,2),
     &                              1.0d+0,ip)
            trRES(ip,ie,it,ir,iz,1) = Xa
            trRES(ip,ie,it,ir,iz,2) = sigx
           end if

           ! X_bar = Xj_bar
           tranatal(ip,ie,it,ir,iz,iat(1,ntf)) =
     &        trRES(ip,ie,it,ir,iz,1)

           ! (sig_xj)**2
           tranatal(ip,ie,it,ir,iz,iat(2,ntf)) =
     &     (  trRES(ip,ie,it,ir,iz,2)
     &      * trRES(ip,ie,it,ir,iz,1) )**2

           ! sig_x
           tranatal(ip,ie,it,ir,iz,iat(2,ntf)) = dsqrt(
     &     tranatal(ip,ie,it,ir,iz,iat(2,ntf)))

           ! Sigma(xi wi)**2 = (sig_X**2 N(N-1)+N X_bar**2)(W/N)**2
           tranatal(ip,ie,it,ir,iz,iat(2,ntf)) =
     &    (tranatal(ip,ie,it,ir,iz,iat(2,ntf))**2 *
     &      resc3(m) * (resc3(m)-1.0d0) +
     &      resc3(m) *
     &     tranatal(ip,ie,it,ir,iz,iat(1,ntf))**2) *
     &     (resc2(m)/resc3(m))**2

           ! Sigma xi wi = X_bar W
           tranatal(ip,ie,it,ir,iz,iat(1,ntf)) =
     &     tranatal(ip,ie,it,ir,iz,iat(1,ntf)) * resc2(m)

           if( tranatal(ip,ie,it,ir,iz,iat(1,ntf)) > cmax )
     &  cmax = tranatal(ip,ie,it,ir,iz,iat(1,ntf))
           if( tranatal(ip,ie,it,ir,iz,iat(1,ntf)) < cmin )
     &  cmin = tranatal(ip,ie,it,ir,iz,iat(1,ntf))

          end do   ! ip loop end
         end do    ! ie loop end
        end do     ! it loop end
       end do      ! ir loop end
      end do       ! iz loop end

! sumover
      call anatal_calc_deposit_anova_sub(m,ntf,np,
     &                               resc2(m),resc3(m) )

!-----------------------------------------------------------------------
      return

      end subroutine anatal_calc_anova_depstrz


!***********************************************************************
!                                                                      *
      subroutine anatal_calc_anova_depstxyz(m,ntf,
     &                             np,nei,ne,nt,nx,ny,nz,
     &                             sumfactor,weightRate,
     &                             trRES,tranatal,nfile,
     &                             resc2SUMTAL,resc3SUMTAL,
     &                             ierr)
!                                                                      *
!     setting of tranatal of t-deposit with xyz mesh                   *
!     create by S. Hashimoto on 2019.8.1                               *
!                                                                      *
!***********************************************************************
      implicit none
!-----------------------------------------------------------------------
      include 'param.inc'
      include 'err.inc'
!-----------------------------------------------------------------------
      integer   istdevres
      integer   maxcasres
      real*8    rijklstres
      integer   irdrf
      common /res01/ istdevres,maxcasres,rijklstres,irdrf

      double precision resc2   ! W
      double precision resc3   ! N
      common /restart/ resc2(itlmax), resc3(itlmax)
!-----------------------------------------------------------------------
      character m_err*200
      integer   l_err, k_err
      common /error/ m_err, l_err, k_err
!-----------------------------------------------------------------------
      integer           m
      integer           ntf
      integer           np
      integer           ne
      integer           nt
      integer           nx
      integer           ny
      integer           nz
      integer           nei

      double precision  sumfactor     ! normalization factor
      double precision  weightRate(*) ! weighting rate
      double precision  trRES   (np,0:ne,nt,nx*ny*nz,2)
      integer           nfile
      double precision  tranatal(np,0:ne,nt,nx*ny*nz,2*nfile)
      double precision  resc2SUMTAL
      double precision  resc3SUMTAL

      double precision  cmin, cmax

      integer           ierr
!-----------------------------------------------------------------------
      integer           ix
      integer           iy
      integer           iz

      integer           ip
      integer           ie
      integer           it
      double precision  Xa, sigx
!-----------------------------------------------------------------------
      integer           icf
      integer           iat, iad, iaf

      icf(ix,iy,iz) = ix + ( iy - 1 ) * nx + ( iz - 1 ) * nx * ny
      iat(iad,iaf) = iad + (iaf-1) * 2
!-----------------------------------------------------------------------

      ierr = 0

*-----------------------------------------------------------------------
      cmax = 0.0
      cmin = 1.e+33

      do iz = 1, nz
       do iy = 1, ny
        do ix = 1, nx
         do it = 1, nt
          do ie = nei, ne
           do ip = 1, np

            if( trRES(ip,ie,it,icf(ix,iy,iz),1) .gt. 0.d0 ) then
              call calc_deposit_stdev(m,Xa,sigx,
     &                                trRES(ip,ie,it,icf(ix,iy,iz),1),
     &                                trRES(ip,ie,it,icf(ix,iy,iz),2),
     &                                1.0d+0,ip)
              trRES(ip,ie,it,icf(ix,iy,iz),1) = Xa
              trRES(ip,ie,it,icf(ix,iy,iz),2) = sigx
            end if

            ! X_bar = Xj_bar
            tranatal(ip,ie,it,icf(ix,iy,iz),iat(1,ntf)) =
     &         trRES(ip,ie,it,icf(ix,iy,iz),1)

            ! (sig_xj)**2
            tranatal(ip,ie,it,icf(ix,iy,iz),iat(2,ntf)) =
     &      (  trRES(ip,ie,it,icf(ix,iy,iz),2)
     &       * trRES(ip,ie,it,icf(ix,iy,iz),1) )**2

            ! sig_x
            tranatal(ip,ie,it,icf(ix,iy,iz),iat(2,ntf)) = dsqrt(
     &      tranatal(ip,ie,it,icf(ix,iy,iz),iat(2,ntf)))

            ! Sigma(xi wi)**2 = (sig_X**2 N(N-1)+N X_bar**2)(W/N)**2
            tranatal(ip,ie,it,icf(ix,iy,iz),iat(2,ntf)) =
     &     (tranatal(ip,ie,it,icf(ix,iy,iz),iat(2,ntf))**2 *
     &      resc3(m) * (resc3(m)-1.0d0) +
     &      resc3(m) *
     &      tranatal(ip,ie,it,icf(ix,iy,iz),iat(1,ntf))**2) *
     &      (resc2(m)/resc3(m))**2

            ! Sigma xi wi = X_bar W
            tranatal(ip,ie,it,icf(ix,iy,iz),iat(1,ntf)) =
     &      tranatal(ip,ie,it,icf(ix,iy,iz),iat(1,ntf)) * resc2(m)

            if( tranatal(ip,ie,it,icf(ix,iy,iz),iat(1,ntf)) > cmax )
     &   cmax = tranatal(ip,ie,it,icf(ix,iy,iz),iat(1,ntf))
            if( tranatal(ip,ie,it,icf(ix,iy,iz),iat(1,ntf)) < cmin )
     &   cmin = tranatal(ip,ie,it,icf(ix,iy,iz),iat(1,ntf))

           end do   ! ip loop end
          end do    ! ie loop end
         end do     ! it loop end
        end do      ! ix loop end
       end do       ! iy loop end
      end do        ! iz loop end

! sumover
      call anatal_calc_deposit_anova_sub(m,ntf,np,
     &                               resc2(m),resc3(m) )

!-----------------------------------------------------------------------
      return

      end subroutine anatal_calc_anova_depstxyz
!***********************************************************************


************************************************************************
*                                                                      *
      subroutine anatal_pdepstreg(m,np,nr,mr,ne,nt,kr,eb,tb,tr,nfile,
     &                    weightRate,
     &                    nvl,ivl,rvl,
     &                    nx,ny,nz,xm,ym,zm,igsh,idasa,manatally)
*                                                                      *
*       anatally calculation and output of t-deposit with reg mesh     *
*                                                                      *
************************************************************************
      use sumtallymod , only : isumtally
      use sangelmod, only: itsans
      use anatallymod, only: anataldata,fg,fgaxs,fgaxs2,fgaxs3,delvol,
     &     cx_txt,lcx_txt,cy_txt,lcy_txt,cz_txt,lcz_txt
      use partmod ! frtati 2021/10/05

*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

      common /stat / istdev, irestart, ireschk

*-----------------------------------------------------------------------

      include 'param.inc'

*-----------------------------------------------------------------------

      double precision resc2   ! W
      double precision resc3   ! N
      common /restart/ resc2(itlmax), resc3(itlmax)

*-----------------------------------------------------------------------

      common /taliin/ rsouin, nzztin, nrgnin

*-----------------------------------------------------------------------

      common /verjam/ versn, lastr, iyeav, imonv, idayv

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

      common /tall20/ ittty(itlmax), ittnm(itlmax), ittrg(itlmax),
     &                rttmi(itlmax), rttma(itlmax), rttdl(itlmax)
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

      common /tall70/ itrwgtsum(itlmax), itrgn1(itlmax),
     &                itrgm1(itlmax), itrncd(itlmax)

      common /istcut/ ist_cut, ist_bat
*-----------------------------------------------------------------------

      common /volreg/ dvol(kvlmax)

*-----------------------------------------------------------------------

      common /tcntl/  icntl, inucr

*-----------------------------------------------------------------------

      common /cparm/  maxbch,maxcas
      common /jcomon/ nabov,nobch,nocas,nomax
!$OMP THREADPRIVATE(/jcomon/)
!$    integer ipomp,npomp
!$    integer OMP_GET_NUM_THREADS,OMP_GET_THREAD_NUM

      common /talout/ itall

      character fname*100, fnume*3

*-----------------------------------------------------------------------
      character gfnam*100
      common /tall60/ itger(itlmax)

      common /tall67/ itstd(itlmax), rtstd(itlmax)

*-----------------------------------------------------------------------

      dimension   kr(mr)
      dimension   eb(ne+1)
      dimension   tb(nt+1)
      dimension   tw(nt)
!      dimension   vl(nr,1,1)
      dimension   lr(nr)
      dimension   tr(np,0:ne,nr,nt,2*nfile)
!      dimension   weightRate(nfile), rdata(2,nfile)
      dimension   weightRate(nfile)
      dimension   ivl(nvl)
      dimension   rvl(nvl)

      dimension   xm(nx+1)
      dimension   ym(ny+1)
      dimension   zm(nz+1)
!      dimension   val(nr)

      integer,allocatable :: ixyz(:)

      integer irst,nrst, itmpdata
!      dimension   anatalrst(np,ne+1,1,nt+1,nr*1*1,1,2*nfile+3)
! sumover
!      dimension   anatalrst(np,ne+1,1,nt+1,(nr+1)*1*1,1,2*nfile+3)
      real(8),allocatable :: anatalrst(:,:,:,:,:,:,:)
      integer ianataldata,nanataldata
      integer, parameter :: ntaxis=9
      integer ij(ntaxis),nij(ntaxis),ibin(ntaxis)
      integer nijaxs,nijaxs2,nijaxs3,nfgmax
      character cij(ntaxis)*27,cijaxs*27,cijaxs2*27
      integer iDaxis
!      dimension   dlr(nr+1) !  integer lr(nr) ---> double dlr(nr)
      real(8),allocatable :: vl(:,:,:),rdata(:,:),val(:),dlr(:)

*-----------------------------------------------------------------------
      real(8) cfac
      real(8),parameter:: c2gy=1.602d-10,c2jperm3=1.602d-7

      dimension tott(np,2) ! frtati 2021/10/05 6 -> np

      character hsunit(0:5)*22
      data hsunit( 0) / 'Dose [Gy/source]      '/
      data hsunit( 1) / 'Dose [MeV/cm^3/source]'/
      data hsunit( 2) / 'Dose [MeV/source]     '/
      data hsunit( 3) / 'Number [1/source]     '/
      data hsunit( 4) / 'Number [1/nsec/source]'/
      data hsunit( 5) / 'Dose [J/m^3/source]   '/

      character cha*1
      data cha /"'"/

      character dc2*4
      character chp(itmxpt)*11  ! kitamura22/03/31
      character chq(itmxpt)*9   ! kitamura22/03/31

      character rpa*1
      data rpa /'}'/
      character yen*1

      dimension dt_one(1)
      data dt_one/1.0d0/
!-----------------------------------------------------------------------
      integer           iat, iad, iaf

      iat(iad,iaf) = iad + (iaf-1) * 2
!-----------------------------------------------------------------------
      yen  = char(92)
      allocate( ixyz(max((nx+1)*(ny+1),(ny+1)*(nz+1),(nx+1)*(nz+1))) )

      allocate (vl(nr,1,1),rdata(2,nfile),val(nr),dlr(nr+1))

      nrst = 2*nfile+3
      allocate (anatalrst(np,ne+1,1,nt+1,(nr+1)*1*1,1,nrst))

*-----------------------------------------------------------------------
      np_mxang = np
      if ( itmxang(m).ne.0 .and. iabs(itmxang(m)).lt.np ) then
        np_mxang = iabs(itmxang(m))
      end if

*-----------------------------------------------------------------------

         if( itrwgtsum(m) .eq. 1 ) then
            nr0 = 1
         else
            nr0 = nr
         endif

*-----------------------------------------------------------------------

         if( itout(m) .le. 1 ) then

            nei = 0

         else

            nei = 1

         end if

*-----------------------------------------------------------------------
*     for igshow
*-----------------------------------------------------------------------

         if( igsh .eq. 0 ) then

            npg = np
            neg = ne
            ntg = nt

         else

            npg = 1
            neg = nei
            ntg = 1

         end if

*-----------------------------------------------------------------------
*        particle name
*-----------------------------------------------------------------------

               call wtpname(m,np,chp,chq)

*-----------------------------------------------------------------------
*           itunt(m) = 0 : Gy/source
*                    = 1 : MeV/cm^3/source
*                    = 2 : MeV/source
*                    = 3 : 1/source
*                    = 4 : 1/nsec/source
*                    = 5 : J/m^3/source
*-----------------------------------------------------------------------
*        set mesh volume and name
*-----------------------------------------------------------------------

               call tregvl(nr,mr,kr,vl,lr,nvl,ivl,rvl)
               vl_sum = sum(vl(:,1,1))

         if( itrwgtsum(m) .eq. 1 ) then
            vlsum = 0.d0
            do ir = itrgn1(m)+1, nr
               vlsum = vlsum + vl(ir,1,1)
            enddo
            vl(1,1,1) = vlsum
         endif

ccse 2021.08 add (use anatal_rearrange sub.)
         do ir = 1, nr
            dlr(ir) = dble(lr(ir))
         end do
         dlr(nr+1) = dble(lr(nr))

*-----------------------------------------------------------------------
*        set volume ( unit = 2, 3, 4 ; vol = 1.0 )
*-----------------------------------------------------------------------

            if( itunt(m).ge.2 .and. itunt(m).le.4 )then

               do ir = 1, nr

                  vl(ir,1,1) = 1.0d0

               end do
               vl_sum = 1.0d0

            end if

*-----------------------------------------------------------------------
*           itunt(m) = 0 ; cfac = c2gy
*           itunt(m) = 5 ; cfac = c2jperm3
*-----------------------------------------------------------------------
            if( itunt(m) .eq. 0)then
              cfac = c2gy
            elseif( itunt(m) .eq. 5)then
              cfac = c2jperm3
            else
              cfac = 1.0d0
            endif

*-----------------------------------------------------------------------
*           itunt(m) = 4 : /nsec/source
*-----------------------------------------------------------------------

            if( itunt(m) .eq. 4 ) then

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
*        relative error and  unit conversion
*        c1 : normalization for source
*-----------------------------------------------------------------------

         if( igsh .eq. 0 ) then

                  cmax = 0.0
                  cmin = 1.e+33
                  dnon = 1.e-33
                  isdz = 0
                  stdm = 0.0
                  itstd(m) = 1

                if( isumtally(m)/=2 .or. itout(m)==1 ) then
                 if ( rsouin .gt. 0d0 ) then ! S.H. (2016.7.26)
                  c1 = 1.0d+0 / rsouin
                 else
                  c1 = 0d0
                 end if
                else          ! deposit, weighted average
                  c1 = 1.0d+0
                end if

      do ir=1,nr0
            do 100 it = 1, nt
            do 100 ie = nei, ne
            do 100 ip = 1, np
            do ntf = 1, nfile

             if( isumtally(m)/=2 .or. itout(m)==1 ) then

               if( tr(ip,ie,ir,it,iat(1,ntf)) .gt. 0.d0 ) then

                  call calc_deposit_stdev(m,Xa,sigx,
     &                            tr(ip,ie,ir,it,iat(1,ntf)),
     &                            tr(ip,ie,ir,it,iat(2,ntf)),
     &                            rtfac(m)/vl(ir,1,1)/tw(it),ip)
                  tr(ip,ie,ir,it,iat(1,ntf)) = cfac*Xa
                  tr(ip,ie,ir,it,iat(2,ntf)) = sigx
                  if( sigx .gt. stdm ) stdm = sigx

                  if( tr(ip,ie,ir,it,iat(1,ntf)) .gt. cmax )
     &         cmax = tr(ip,ie,ir,it,iat(1,ntf))

                  if( tr(ip,ie,ir,it,iat(1,ntf)) .lt. cmin )
     &         cmin = tr(ip,ie,ir,it,iat(1,ntf))

               else

                  isdz = 1
                  tr(ip,ie,ir,it,iat(2,ntf)) = 0.0

               end if
             else
               if( tr(ip,ie,ir,it,iat(1,ntf)) .gt. 0.d0 ) then
                 if( tr(ip,ie,ir,it,iat(1,ntf)) .gt. cmax
     &            .and. manatally .ne. 0 )  ! T.Sato 2021/01/30
     &        cmax = tr(ip,ie,ir,it,iat(1,ntf))
                 if( tr(ip,ie,ir,it,iat(1,ntf)) .lt. cmin
     &            .and. manatally .ne. 0 )  ! T.Sato 2021/01/30
     &        cmin = tr(ip,ie,ir,it,iat(1,ntf))
               else
                  isdz = 1
                  tr(ip,ie,ir,it,iat(2,ntf)) = 0.0
               end if
             end if

! sumover
              call pdepstreg_sumover_stdev_ntf(0,m,ntf,
     &             ip,ie,ir,it,
     &             nei,
     &             rtfac(m),cfac,tw(it),vl(ir,1,1),tw_sum,vl_sum)

            end do ! loop end ntf

  100       continue
      enddo  ! ir loop T.Sato 2021/04/15

!$OMP PARALLEL
!$OMP& private(ipomp,ir,it,ie,ip,ntf,rdata,answer,rerr)
!$OMP& private(iax,ir_a,it_a,ie_a,ido_ana,ioe)
!$      npomp=OMP_GET_NUM_THREADS()
!$      ipomp=OMP_GET_THREAD_NUM()
!$      write(*,'(''OpenMP PARALLEL PROCESS'',
!$   &     i4,''/'',i4)') ipOMP+1,npomp
!$OMP DO SCHEDULE(dynamic) reduction(max:cmax) reduction(min:cmin)
      do ir=1,nr0
            do 101 it = 1, nt
            do 101 ie = nei, ne
            do 101 ip = 1, np

            if( manatally .eq. 0 ) then ! T.Sato 2020/10/19

              do ntf = 1, nfile
                rdata(1,ntf)=tr(ip,ie,ir,it,iat(1,ntf))
                rdata(2,ntf)=tr(ip,ie,ir,it,iat(2,ntf))
              end do

              call usranatal(nfile,rdata,answer,rerr)

              tr(ip,ie,ir,it,1) = answer
              tr(ip,ie,ir,it,2) = rerr
             if ( itout(m) .le. 1 ) then ! ie=nei=ne=0
              anatalrst(ip, 1,1,it,ir,1,1) = answer
              anatalrst(ip, 1,1,it,ir,1,2) = rerr
             else
              anatalrst(ip,ie,1,it,ir,1,1) = answer
              anatalrst(ip,ie,1,it,ir,1,2) = rerr
             end if

              if(answer.gt.cmax) cmax = answer
              if(answer.gt.0d0 .and. answer.lt.cmin) cmin = answer

            else if( manatally .eq. 1 ) then

               do ntf = 1, nfile
                  rdata(1,ntf) = tr(ip,ie,ir,it,iat(1,ntf))
                  rdata(2,ntf) = tr(ip,ie,ir,it,iat(2,ntf))
               end do

               call anova(nfile,rdata,resc3(m)
     &              ,fmval,unca,uncerr,uratio)

              if ( itout(m) .le. 1 ) then ! ie=nei=ne=0
               anatalrst(ip,1,1,it,ir,1,1) ! mean
     &              = fmval
               if ( fmval .gt. 0d0 ) then
                  anatalrst(ip,1,1,it,ir,1,2) ! total uncertainty
     &                 = dsqrt(unca**2 + uncerr**2/resc3(m)) / fmval
                  anatalrst(ip,1,1,it,ir,1,3) ! syst. uncertainty
     &                 = unca / fmval
                  anatalrst(ip,1,1,it,ir,1,4) ! stat. uncertainty
     &                 = uncerr /dsqrt(resc3(m)) / fmval
                  anatalrst(ip,1,1,it,ir,1,5) ! squared unc. ratio
     &                 = uratio
               else
                  anatalrst(ip,1,1,it,ir,1,2) ! total uncertainty
     &                 = 0d0
                  anatalrst(ip,1,1,it,ir,1,3) ! syst. uncertainty
     &                 = 0d0
                  anatalrst(ip,1,1,it,ir,1,4) ! stat. uncertainty
     &                 = 0d0
                  anatalrst(ip,1,1,it,ir,1,5) ! squared unc. ratio
     &                 = 0d0
               end if
              else
               anatalrst(ip,ie,1,it,ir,1,1) ! mean
     &              = fmval
               if ( fmval .gt. 0d0 ) then
                  anatalrst(ip,ie,1,it,ir,1,2) ! total uncertainty
     &                 = dsqrt(unca**2 + uncerr**2/resc3(m)) / fmval
                  anatalrst(ip,ie,1,it,ir,1,3) ! syst. uncertainty
     &                 = unca / fmval
                  anatalrst(ip,ie,1,it,ir,1,4) ! stat. uncertainty
     &                 = uncerr /dsqrt(resc3(m)) / fmval
                  anatalrst(ip,ie,1,it,ir,1,5) ! squared unc. ratio
     &                 = uratio
               else
                  anatalrst(ip,ie,1,it,ir,1,2) ! total uncertainty
     &                 = 0d0
                  anatalrst(ip,ie,1,it,ir,1,3) ! syst. uncertainty
     &                 = 0d0
                  anatalrst(ip,ie,1,it,ir,1,4) ! stat. uncertainty
     &                 = 0d0
                  anatalrst(ip,ie,1,it,ir,1,5) ! squared unc. ratio
     &                 = 0d0
               end if
              end if

            else if( manatally .eq. 2 ) then

               do ntf = 1, nfile
                do ioe = 1, 2
                 if ( itout(m) .le. 1 ) then ! ie=nei=ne=0
                  anatalrst(ip,1,1,it,ir,1,iat(ioe,ntf))
     &              = tr(ip,ie,ir,it,iat(ioe,ntf))
                 else
                  anatalrst(ip,ie,1,it,ir,1,iat(ioe,ntf))
     &              = tr(ip,ie,ir,it,iat(ioe,ntf))
                 end if
                end do
               end do

            end if

!sumover
            do iax = 1,itaxn(m)

              ie_a = ie
              it_a = it
              ir_a = ir
              ido_ana = 0
              if(itaxs(m,iax) == 1 .and. ie == 1) then
                ie_a = ne + 1
                ido_ana = 1
              elseif(itaxs(m,iax) == 2 .and. ir == 1) then
                ir_a = nr + 1
                ido_ana = 1
              elseif(itaxs(m,iax) == 11 .and. it == 1) then
                it_a = nt + 1
                ido_ana = 1
              endif

              if(ido_ana == 1) then

                 call get_depstreg_tr_sum_data(m,iax,nfile,
     &                 ip,ie,ir,it,rdata)

                if( manatally .eq. 0 ) then ! T.Sato 2020/10/22

                   call usranatal(nfile,rdata,answer,rerr)

                   if ( itout(m) .le. 1 ) then ! ie=nei=ne=0
                     anatalrst(ip,1,1,it_a,ir_a,1,1) = answer
                     anatalrst(ip,1,1,it_a,ir_a,1,2) = rerr
                   else
                     anatalrst(ip,ie_a,1,it_a,ir_a,1,1) = answer
                     anatalrst(ip,ie_a,1,it_a,ir_a,1,2) = rerr
                   endif

                else if( manatally .eq. 1 ) then

                  call anova(nfile,rdata,resc3(m)
     &               ,fmval,unca,uncerr,uratio)

                  if ( itout(m) .le. 1 ) then ! ie=nei=ne=0

                    anatalrst(ip,1,1,it_a,ir_a,1,1) ! mean
     &                = fmval
                    if ( fmval .gt. 0.0d0 ) then
                      anatalrst(ip,1,1,it_a,ir_a,1,2) ! total uncertainty
     &                     = dsqrt(unca**2 + uncerr**2/resc3(m)) / fmval
                    else
                      anatalrst(ip,1,1,it_a,ir_a,1,2) ! total uncertainty
     &                   = 0.0d0
                    end if
                  else
                    anatalrst(ip,ie_a,1,it_a,ir_a,1,1) ! mean
     &                = fmval
                    if ( fmval .gt. 0.0d0 ) then
                      anatalrst(ip,ie_a,1,it_a,ir_a,1,2) ! total uncertainty
     &                     = dsqrt(unca**2 + uncerr**2/resc3(m)) / fmval
                    else
                      anatalrst(ip,ie_a,1,it_a,ir_a,1,2) ! total uncertainty
     &                   = 0.0d0
                    end if
                  endif

                else if( manatally .eq. 2 ) then
                  do ntf = 1, nfile
                    do ioe = 1, 2
                      if ( itout(m) .le. 1 ) then ! ie=nei=ne=0
                        anatalrst(ip,1,1,it_a,ir_a,1,iat(ioe,ntf))
     &                  = rdata(ioe,ntf)
                      else
                        anatalrst(ip,ie_a,1,it_a,ir_a,1,iat(ioe,ntf))
     &                  = rdata(ioe,ntf)
                      endif
                    end do
                  end do

                end if

              endif
            enddo

  101       continue
      enddo  ! ir loop T.Sato 2021/04/15
!$OMP END DO
!$OMP END PARALLEL

            if( istdev .eq. 2 .or. nobch .gt. 20 ) then
            if( rtstd(m) .gt. 0.d0 .and.
     &          stdm .gt. 0.d0 .and. isdz .eq. 0 .and.
     &          stdm .lt. rtstd(m) ) itstd(m) = 0
            end if

         end if   ! igsh

*-----------------------------------------------------------------------
*     out put unit = 31  : temporary number
*-----------------------------------------------------------------------

      do 900 iax = 1, itfln(m)

!OBINATA(2012.8.3): output *.err
        noe = 1
        if ( any( itaxs(m,iax) .eq. (/ 12, 13 /) )
     &      .and. ittwo(m) .ne. 4 .and. icntl .ne. 8 )
     &    noe = 2

        do 900 ioe = 1, noe

         if( igsh .ne. 0 .and.
     &     ( itaxs(m,iax) .lt. 7 .or. itrsh(m) .eq. 0 ) ) goto 900

         if( itall .eq. 2 .and. nobch .lt. maxbch .and.
     &       igsh .eq. 0 ) then

            write(fnume,'(i3.3)') nobch
!OBINATA(2012.8.3): output *.err
            if ( ioe .eq. 1 ) then
              fname = ctfln(m,iax)(1:itfll(m,iax))//'.'//fnume
            else
              call mk_2derrfn(ctfln(m,iax),fname,itfll(m,iax))
            end if

         else

!OBINATA(2012.8.3): output *.err
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

         call depstech(iot,m,iax,1)

*-----------------------------------------------------------------------
* manatally
*  0  :  user defined analysis
*  1  :  systematic uncertainty
*  2  :  c-value dependence
*-----------------------------------------------------------------------

*-----------------------------------------------------------------------
*        energy, reg, time axis
*-----------------------------------------------------------------------

         if( itaxs(m,iax) .eq. 1
     &        .or. itaxs(m,iax) .eq. 2
     &        .or. itaxs(m,iax) .eq. 11 ) then

          select case( itaxs(m,iax) )
           case ( 1 )
            iDaxis = 1   ! energy
           case ( 2 )
            iDaxis = 2   ! reg axis
           case ( 11 )
            iDaxis = 9   ! time
          end select
          if( manatally.eq.2 ) iDaxis = 51 ! S.H. 2021.8.15

C rearrange data from anatalrst to anataldata
          call anatal_rearrange_sum(0,0,1,
     &         iDaxis,ntaxis,anatalrst,
     &         np,     ne,       0,    nt,    nr0,  0,  0,   0, nrst,
     &         eb,ew,    1,1,  tb,tw, dlr,  dt_one,  dt_one,1,
     &         vl,
     &         nij,ibin,nijaxs,nijaxs2,nijaxs3,nfgmax,
     &         cij,cijaxs,cijaxs2,nanataldata,ierr)

          inum = 0

          do ianataldata = 1, nanataldata

           inum = inum + 1

C exchange ID number from ianataldata to ij2,...,ij8
           itmpdata = int(ianataldata-1)
           ij(8) = mod(itmpdata,nij(8))+1
           itmpdata = int((ianataldata-ij(8))/nij(8))
           ij(7) = mod(itmpdata,nij(7))+1
           itmpdata = int((itmpdata-ij(7)+1)/nij(7))
           ij(6) = mod(itmpdata,nij(6))+1
           itmpdata = int((itmpdata-ij(6)+1)/nij(6))
           ij(5) = mod(itmpdata,nij(5))+1
           itmpdata = int((itmpdata-ij(5)+1)/nij(5))
           ij(4) = mod(itmpdata,nij(4))+1
           itmpdata = int((itmpdata-ij(4)+1)/nij(4))
           ij(3) = mod(itmpdata,nij(3))+1
           itmpdata = int((itmpdata-ij(3)+1)/nij(3))
           ij(2) = mod(itmpdata,nij(2))+1

*-----------------------------------------------------------------------

           write(iot,'(/"#",78("-"))')

           if( inum .eq. 1 ) then
              write(iot,'( "#newpage:")')
           else
              write(iot,'( " newpage:")')
           end if

           if( lcx_txt(m) .gt. 0 ) then
              write(iot,'(/"x: ",200a1)')(cx_txt(m)(i:i),i=1,lcx_txt(m))
           else
            if( itaxl(m) .eq. 0 ) then
               write(iot,'(/"x: ",a27)') cijaxs
            else
               write(iot,'(/"x: ",200a1)') (itaxt(m)(i:i),i=1,itaxl(m))
            end if
           end if

           if( lcy_txt(m) .gt. 0 ) then
              write(iot,'( "y: ",200a1)')(cy_txt(m)(i:i),i=1,lcy_txt(m))
           else
            if( itayl(m) .eq. 0 ) then
               write(iot,'( "y: ",a22)') hsunit(itunt(m))
            else
               write(iot,'( "y: ",200a1)') (itayt(m)(i:i),i=1,itayl(m))
            end if
           end if

           if( (itaxs(m,iax) .eq. 1).and.
     &          (itety(m) .eq. 3 .or. itety(m) .eq. 5) ) then
              write(iot,'( "p: xlog ylog afac(0.8) form(0.9)")')

           else if( itaxs(m,iax) .eq. 2 ) then
              write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')

           else if( (itaxs(m,iax) .eq. 11) .and.
     &             (ittty(m) .eq. 3 .or. ittty(m) .eq. 5) ) then
              write(iot,'( "p: xlog ylog afac(0.8) form(0.9)")')

           else
              write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')

           end if

           if( itanl(m) .gt. 0 ) then
              write(iot,'( "p: ",200a1)')
     &             ( itang(m)(i:i),i = 1, itanl(m) )
           end if

           if( itsans(m) .gt. 0 ) then
              call write_sangel(iot,m,0)
           end if

           if ( manatally .eq. 0 ) then ! user defined analysis
            if ( iteps(m) .ne. 2 ) then
             if( itaxs(m,iax) .eq. 2 ) then ! reg axis
                write(iot,'( "h:   x      n",5x,"n",10x,
     &                        1000(a10,"),",a4," n   "))')
     &                         ( chp(i), chm(i) , i = 1, np_mxang ) ! frtati 2021/10/05 np -> np_mxang
             else
                write(iot,'( "h: n",12x,"x",12x,
     &                        1000(a10,"),hh0",a3," n "))')
     &                         ( chp(i), chl(i) , i = 1, np_mxang ) ! frtati 2021/10/05 np -> np_mxang
             end if

            else
             if( itaxs(m,iax) .eq. 2 ) then ! reg axis
                write(iot,'( "h:   x      n",5x,"n",10x,
     &             1000(a1,i1,a9,"),",a4
     &             ," ny",i0," dy",i0,"=[y",i0,"*y",i0,"] "))')
     &             ( chp(i)(1:1),i,chp(i)(2:10), chm(i)
     &             ,i+np_mxang,i,i,i+np_mxang, i = 1, np_mxang ) ! frtati 2021/10/05
             else
                write(iot,'( "h: n",12x,"x",12x,
     &             1000(a1,i1,a9,"),hh0",a3
     &             ," ny",i0," dy",i0,"=[y",i0,"*y",i0,"] "))')
     &             ( chp(i)(1:1),i,chp(i)(2:10), chl(i)
     &             ,i+np_mxang,i,i,i+np_mxang, i = 1, np_mxang ) ! frtati 2021/10/05
             end if
            end if

            if( itaxs(m,iax) .eq. 2 ) then ! reg axis
               write(iot,'( "#  num    reg     volume  ",
     &                      1000(a1,"  ",a8,"    r.err "))')
     &                       ( chb(i), chp(i)(3:10), i = 1, np ) ! frtati 2021/10/05
            else
               write(iot,'( "#  lower        upper  ",3x,
     &                      1000(a1,2x,a8,4x,"r.err "))')
     &                       ( chb(i), chp(i)(3:10), i = 1, np ) ! frtati 2021/10/05
            end if

           else if ( manatally .eq. 1 ) then ! systematic uncertainty
            if ( iteps(m) .ne. 2 ) then
             if( itaxs(m,iax) .eq. 2 ) then ! reg axis
                write(iot,'( "h:   x      n",5x,"n",10x,
     &                        1000(a10,"),",a4," n n n "))')
     &                         ( chp(i), chm(i) , i = 1, np_mxang ) ! frtati 2021/10/05 np -> np_mxang
             else
                write(iot,'( "h: n",12x,"x",12x,
     &                        1000(a10,"),hh0",a3," n n n "))')
     &                         ( chp(i), chl(i) , i = 1, np_mxang ) ! frtati 2021/10/05 np -> np_mxang
             end if

            else
             if( itaxs(m,iax) .eq. 2 ) then ! reg axis
                write(iot,'( "h:   x      n",5x,"n",10x,
     &             1000(a1,i1,a9,"),",a4
     &             ," ny",i0," dy",i0,"=[y",i0,"*y",i0,"] n n "))')
     &             ( chp(i)(1:1),i,chp(i)(2:10), chm(i)
     &             ,i+np_mxang,i,i,i+np_mxang, i = 1, np_mxang ) ! frtati 2021/10/05
             else
                write(iot,'( "h: n",12x,"x",12x,
     &             1000(a1,i1,a9,"),hh0",a3
     &             ," ny",i0," dy",i0,"=[y",i0,"*y",i0,"] n n "))')
     &             ( chp(i)(1:1),i,chp(i)(2:10), chl(i)
     &             ,i+np_mxang,i,i,i+np_mxang, i = 1, np_mxang ) ! frtati 2021/10/05
             end if
            end if

            if( itaxs(m,iax) .eq. 2 ) then ! reg axis
               write(iot,'( "#  num    reg     volume  ",
     &                    1000(a1,2x,a8,4x,"r.err(tot, syst, stat)"))')
     &                        ( chb(i), chp(i)(3:10), i = 1, np ) ! frtati 2021/10/05
            else
               write(iot,'( "#  lower        upper  ",3x,
     &                    1000(a1,2x,a8,4x,"r.err(tot, syst, stat)"))')
     &                        ( chb(i), chp(i)(3:10), i = 1, np ) ! frtati 2021/10/05
            end if

           else if ( manatally .eq. 2 ) then ! c-value dependence
            if ( iteps(m) .ne. 2 ) then
               write(iot,'( "h:   x",10x,
     &                    1000(a10,"),",a4," n   "))')
     &              ( chp(i), chm(i) , i = 1, np_mxang ) ! frtati 2021/10/05 np -> np_mxang

            else
               write(iot,'( "h:   x",10x,
     &              1000(a1,i1,a9,"),",a4
     &              ," ny",i0," dy",i0,"=[y",i0,"*y",i0,"] "))')
     &              ( chp(i)(1:1),i,chp(i)(2:10), chm(i)
     &              ,i+np_mxang,i,i,i+np_mxang, i = 1, np_mxang ) ! frtati 2021/10/05
            end if

            write(iot,'( "#  c-value   ",
     &                 1000(a1,"  ",a8,"    r.err "))')
     &                 ( chb(i), chp(i)(3:10), i = 1, np ) ! frtati 2021/10/05

           end if

           do ip = 1, np
              tott(ip,1) = 0.d+0
           end do

           voll = 0.0d0
           njaxs = nijaxs + 1

           do ijaxs = 1, nijaxs

            if( itrwgtsum(m) .eq. 1 ) then
             if( nr-itrgn1(m) .gt. 1 ) then
                ireg = 1000001
             else
                ireg = lr(nr)
             end if
            else
               ireg = lr(ijaxs)
            end if

            if ( manatally .eq. 0 ) then ! user defined analysis
             if( itaxs(m,iax) .eq. 2 ) then ! reg axis
                write(iot,'(i5,1x,i7,1pe13.4,1000(1pe13.4,0p3f8.4))')
     &               ijaxs, ireg, fgaxs(ijaxs),
     &               ((anataldata(ip,ianataldata,ijaxs,irst),irst=1,2)
     &               ,ip=1,np) ! frtati 2021/10/05
             else
                write(iot,'(1p2e13.4,1000(1pe13.4,0pf8.4))')
     &               fgaxs(ijaxs),fgaxs(ijaxs+1),
     &               ((anataldata(ip,ianataldata,ijaxs,irst),irst=1,2)
     &               ,ip=1,np) ! frtati 2021/10/05
             end if

            else if ( manatally .eq. 1 ) then ! systematic uncertainty
             if( itaxs(m,iax) .eq. 2 ) then ! reg axis
                write(iot,'(i5,1x,i7,1pe13.4,1000(1pe13.4,0p3f8.4))')
     &               ijaxs, ireg, fgaxs(ijaxs),
     &               ((anataldata(ip,ianataldata,ijaxs,irst),irst=1,4)
     &               ,ip=1,np) ! frtati 2021/10/05
             else
                write(iot,'(1p2e13.4,1000(1pe13.4,0p3f8.4))')
     &               fgaxs(ijaxs),fgaxs(ijaxs+1),
     &               ((anataldata(ip,ianataldata,ijaxs,irst),irst=1,4)
     &               ,ip=1,np) ! frtati 2021/10/05
             end if

            else if ( manatally .eq. 2 ) then ! c-value dependence
               write(iot,'(1pe13.4,1000(1pe13.4,0pf8.4))')
     &              fgaxs(ijaxs),
     &              ((anataldata(ip,ianataldata,ijaxs,irst),irst=1,2)
     &              ,ip=1,np) ! frtati 2021/10/05
            end if

            voll = voll + delvol(ianataldata,ijaxs)

            do ip = 1, np
               vn = anataldata(ip,ianataldata,ijaxs,1)
     &              * delvol(ianataldata,ijaxs)
               tott(ip,1) = tott(ip,1) + vn
            end do

           end do                 ! ijaxs = 1, nijaxs

           if( itaxs(m,iax) .eq. 1 ) then ! when axis = eng
            do ip = 1, np
               tott(ip,1) = tott(ip,1)
            end do

           else if( itaxs(m,iax) .eq. 2 ) then ! when axis = reg
            do ip = 1, np
             if( itunt(m) .eq. 1 .or. itunt(m) .eq. 5 ) then
                tott(ip,1) = tott(ip,1) / voll
             end if
            end do

           else if( itaxs(m,iax) .eq. 11 ) then ! when axis = time
            do ip = 1, np
             if( itunt(m) .eq. 4 ) then
                tott(ip,1) = tott(ip,1) / voll
             end if
            end do

           end if

! sumover
            do ip = 1, np
            do irst=1,2
              tott(ip,irst) = anataldata(ip,ianataldata,njaxs,irst)
            enddo
            enddo

           if( itaxs(m,iax).eq.1 .or. itaxs(m,iax).eq.11 ) then ! when axis = eng,t
            if( itunt(m) .ne. 0 ) then
             if ( manatally .eq. 0 ) then ! user defined analysis
                write(iot,'(/"#   sum over",14x,1000(1pe13.4,0pf8.4))')
     &               (tott(ip,1),tott(ip,2),ip=1,np) ! frtati 2021/10/05
             else if ( manatally .eq. 1 ) then ! systematic uncertainty
             write(iot,'(/"#   sum over",14x,1000(1pe13.4,0pf8.4,16x))')
     &               (tott(ip,1),tott(ip,2),ip=1,np) ! frtati 2021/10/05
             else if ( manatally .eq. 2 ) then ! c-value dependence
             end if
            end if

           else if( itaxs(m,iax).eq.2 ) then ! when axis = reg
            if( itunt(m) .eq. 1 .or. itunt(m) .eq. 5 ) then
             if ( manatally .eq. 0 ) then ! user defined analysis
                write(iot,'(/"#   sum over ",1pe13.4,
     &                     1000(1pe13.4,0pf8.4))')
     &                     voll, (tott(ip,1),tott(ip,2),ip=1,np) ! frtati 2021/10/05
             else if ( manatally .eq. 1 ) then ! systematic uncertainty
                write(iot,'(/"#   sum over ",1pe13.4,
     &                     1000(1pe13.4,0pf8.4,16x))')
     &                     voll, (tott(ip,1),tott(ip,2),ip=1,np) ! frtati 2021/10/05
             else if ( manatally .eq. 2 ) then ! c-value dependence
             end if

            else if( itunt(m) .ne. 0 ) then
             if ( manatally .eq. 0 ) then ! user defined analysis
                write(iot,'(/"#   sum over ",13x,
     &                     1000(1pe13.4,0pf8.4))')
     &                     (tott(ip,1),tott(ip,2),ip=1,np) ! frtati 2021/10/05
             else if ( manatally .eq. 1 ) then ! systematic uncertainty
                write(iot,'(/"#   sum over ",13x,
     &                     1000(1pe13.4,0pf8.4,16x))')
     &                     (tott(ip,1),tott(ip,2),ip=1,np) ! frtati 2021/10/05
             else if ( manatally .eq. 2 ) then ! c-value dependence
             end if
            end if
           end if

*-----------------------------------------------------------------------

           write(iot,'(/a1,"no. =",i3,a1)') cha, inum, cha

           write(iot,'("msuc: {",a1,"huge ",80a1)')
     &          yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

           write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
           write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

C write tally condition of each figure
           write(iot,'("wt: s(0.7)",/a1,"vspace{-3}")') yen
           itmp = 0
           do ijtmp = 2,8
            if (cij(ijtmp) .ne. 'F' ) then
             itmp = itmp + ij(ijtmp)
             if ( ibin(ijtmp) .eq. 0 ) then
                write(iot,'(a15,"&=&",i5)')
     &               cij(ijtmp), idnint(fg(itmp))
                itmp = itmp - ij(ijtmp) + nij(ijtmp)
             else
                write(iot,'(a15,2x,"(",i5,"-th bin)"
     &                             /1pe13.4,2x,"$--$",1pe13.4)')
     &               cij(ijtmp),ij(ijtmp),fg(itmp),fg(itmp+1)
                itmp = itmp - ij(ijtmp) + nij(ijtmp)+1
             end if
            end if
           end do

           write(iot,'("e:")')

          end do                  ! ianataldata = 1, nanataldata

          deallocate( fg )
          deallocate( fgaxs )
          deallocate( anataldata )
          deallocate( delvol )

*-----------------------------------------------------------------------
*        t-eng, eng-t axes (matrix)
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 12
     &       .or. itaxs(m,iax) .eq. 13 ) then

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

          select case( itaxs(m,iax) )
           case ( 12 )
            iDaxis = 35   ! t-eng
           case ( 13 )
            iDaxis = 36   ! eng-t
          end select

C rearrange data from anatalrst to anataldata
          call anatal_rearrange_sum(0,0,1,
     &         iDaxis,ntaxis,anatalrst,
     &         np,     ne,       0,    nt,    nr0,  0,  0,   0, nrst,
     &         eb,ew,    1,1,  tb,tw, dlr,  dt_one,  dt_one,1,
     &         vl,
     &         nij,ibin,nijaxs,nijaxs2,nijaxs3,nfgmax,
     &         cij,cijaxs,cijaxs2,nanataldata,ierr)

          inum = 0

          do ip = 1, np
          do ianataldata = 1, nanataldata

           inum = inum + 1

C exchange ID number from ianataldata to ij2,...,ij8
           itmpdata = int(ianataldata-1)
           ij(8) = mod(itmpdata,nij(8))+1
           itmpdata = int((ianataldata-ij(8))/nij(8))
           ij(7) = mod(itmpdata,nij(7))+1
           itmpdata = int((itmpdata-ij(7)+1)/nij(7))
           ij(6) = mod(itmpdata,nij(6))+1
           itmpdata = int((itmpdata-ij(6)+1)/nij(6))
           ij(5) = mod(itmpdata,nij(5))+1
           itmpdata = int((itmpdata-ij(5)+1)/nij(5))
           ij(4) = mod(itmpdata,nij(4))+1
           itmpdata = int((itmpdata-ij(4)+1)/nij(4))
           ij(3) = mod(itmpdata,nij(3))+1
           itmpdata = int((itmpdata-ij(3)+1)/nij(3))
           ij(2) = mod(itmpdata,nij(2))+1

*-----------------------------------------------------------------------

           write(iot,'(/"#",78("-"))')

           if( inum .eq. 1 .or. ip.gt.np_mxang ) then ! frtati 2021/10/05
              write(iot,'( "#newpage:")')
           else
              write(iot,'( " newpage:")')
           end if

           if ( ip.gt.np_mxang ) then
             write(iot,'( " SKIPPAGE:")')
             inum = inum - 1
           end if

*-----------------------------------------------------------------------

           if( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) then

              write(iot,'(/a1,"no. =",i3,a1)') cha, inum, cha

              write(iot,'("msuc: {",a1,"huge ",80a1)')
     &             yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

              write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
              write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

           if( lcx_txt(m) .gt. 0 ) then
              write(iot,'(/"x: ",200a1)')(cx_txt(m)(i:i),i=1,lcx_txt(m))
           else
              if( itaxl(m) .eq. 0 ) then
                 write(iot,'(/"x: ",a27)') cijaxs
              else
                 write(iot,'(/"x: ",200a1)')
     &                (itaxt(m)(i:i),i=1,itaxl(m))
              end if
           end if

           if( lcy_txt(m) .gt. 0 ) then
              write(iot,'( "y: ",200a1)')(cy_txt(m)(i:i),i=1,lcy_txt(m))
           else
              if( itayl(m) .eq. 0 ) then
                 write(iot,'("y: ",a27)') cijaxs2
              else
                 write(iot,'( "y: ",200a1)')
     &                (itayt(m)(i:i),i=1,itayl(m))
              end if
           end if

*-----------------------------------------------------------------------

              if( inum .eq. 1 ) then

                 form  = ( fgaxs2(nijaxs2+1) - fgaxs2(1) )
     &                / ( fgaxs(nijaxs+1) - fgaxs(1) )
                 xfac  = 0.9
                 afac  = 0.8
                 izlog = 1
                 inocm = 1
                 inolg = 1

                 if( itanl(m) .gt. 0 )
     &                call anset(itang(m),itanl(m),
     &                form,xfac,afac,cmin,cmax,izlog,inocm,inolg)

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
     &             ( ittwo(m) .ge. 6 .and. ittwo(m) .le. 7 ) ) .and.
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
     &              write(iot,'( "p: zlog")')

            end if

*-----------------------------------------------------------------------

            axs1min = fgaxs(1)
            axs1max = fgaxs(nijaxs+1)
            axs2min = fgaxs2(1)
            axs2max = fgaxs2(nijaxs2+1)
            axs1del = (fgaxs(nijaxs+1)-fgaxs(1))/dble(nijaxs)
            axs2del = (fgaxs2(nijaxs2+1)-fgaxs2(1))/dble(nijaxs2)

            write(iot,'( "p: xmin(",1p1g14.7,") xmax(",1p1g14.7,
     &                       ")")') axs1min, axs1max

            write(iot,'( "p: ymin(",1p1g14.7,") ymax(",1p1g14.7,
     &                       ")")') axs2min, axs2max

            if( itanl(m) .gt. 0 .and. ioe .eq. 1 ) then

               write(iot,'( "p: ",200a1)')
     &              ( itang(m)(i:i),i = 1, itanl(m) )

            else if( itanl(m) .gt. 0 .and. ioe .ne. 1 ) then

               call terrang(iot,m)

            end if

            if( itsans(m) .gt. 0 ) then
               call write_sangel(iot,m,0)
            end if

           end if

*-----------------------------------------------------------------------

           write(iot,'("#  ny = ",i3,"   nx = ",i3)')
     &          nijaxs2, nijaxs

           if( ittwo(m) .ne. 4 ) then
              write(iot,'( "# ( ( data(i1,i2), i1 = 1, n1 ),",
     &                         " i2 = n2, 1, -1 )")')
           end if

*-----------------------------------------------------------------------

           if( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) then

              write(iot,'(/a4," y = ",1p1g14.7," to ",
     &            1p1g14.7," by ",1p1g14.7," ; x = ",
     &            1p1g14.7," to ",1p1g14.7," by ",
     &            1p1g14.7," ;")') dc2,
     &            axs2max - axs2del/2d0, axs2min + axs2del/2d0, axs2del,
     &            axs1min + axs1del/2d0, axs1max - axs1del/2d0, axs1del

!OBINATA(2012.8.3): output *.err
              write(iot,'(1p10e11.3)')
     &       ( ( anataldata(ip,ianataldata,ijaxs+(ijaxs2-1)*nijaxs,ioe),
     &                ijaxs = 1, nijaxs ), ijaxs2 = nijaxs2, 1, -1 )

*-----------------------------------------------------------------------

           else if( ittwo(m) .eq. 4 ) then
            write(iot,'(/"# axis2      axis1    ",
     &                      "  number     r.err")')

            do ijaxs = 1, nijaxs
             do ijaxs2 = 1, nijaxs2

                write(iot,'(1p3e13.4,0pf8.4)')
     &               fgaxs2(ijaxs2)  + axs2del/2d0,
     &               fgaxs(ijaxs)  + axs1del/2d0,
     &             anataldata(ip,ianataldata,ijaxs+(ijaxs2-1)*nijaxs,1),
     &             anataldata(ip,ianataldata,ijaxs+(ijaxs2-1)*nijaxs,2)

             end do
            end do

*-----------------------------------------------------------------------

           else if( ittwo(m) .eq. 5 ) then

            write(iot,'("#   t = ",1p1g14.7," to ",
     &            1p1g14.7," by ",1p1g14.7/
     &                        "#   e = ",1p1g14.7," to ",
     &            1p1g14.7," by ",1p1g14.7)')
     &            axs2min + axs2del/2d0, axs2max - axs2del/2d0, axs2del,
     &            axs1min + axs1del/2d0, axs1max - axs1del/2d0, axs1del

            write(iot,'(/4x,a3,4x,1p1000e11.3)')
     &          'ax2/ax1',( fgaxs(ijaxs) + axs1del/2d0, ijaxs=1,nijaxs )

            do ijaxs2 = nijaxs2, 1, -1

!OBINATA(2012.8.3): output *.err
               write(iot,'(1p1000e11.3)')
     &              fgaxs2(ijaxs2)  + axs2del/2d0,
     &         ( anataldata(ip,ianataldata,ijaxs+(ijaxs2-1)*nijaxs,ioe),
     &              ijaxs = 1, nijaxs )

            end do

           end if

*-----------------------------------------------------------------------

      if( inolg .eq. 1 .and. ( ittwo(m) .eq. 3 .or. ittwo(m) .eq. 7 )
     &    .and. cmin .gt. 0.0 .and. cmax .gt. cmin ) then

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

      else

       if( lcz_txt(m) .gt. 0 ) then
          write(iot,'("y: ",200a1)')(cz_txt(m)(i:i),i=1,lcz_txt(m))
       else
        if( itazl(m) .eq. 0 ) then
            write(iot,'("y: ",a22)') hsunit(itunt(m))
        else
            write(iot,'("y: ",200a1)') (itazt(m)(i:i),i=1,itazl(m))
        end if
       end if

      end if

      if( izlog .eq. 0 ) write(iot,'( "p: ylin")')

      write(iot,'( "#",78("-"))')

      end if

*-----------------------------------------------------------------------

      if( inocm .eq. 1 ) then

C write tally condition of each figure
       write(iot,'("wt: s[c5]",/a1,"vspace{-3}")') yen
       itmp = 0
       do ijtmp = 2,8
        if (cij(ijtmp) .ne. 'F' ) then
         itmp = itmp + ij(ijtmp)
         if ( ibin(ijtmp) .eq. 0 ) then
            write(iot,'(a15,"&=&",i5)')
     &           cij(ijtmp), idnint(fg(itmp))
            itmp = itmp - ij(ijtmp) + nij(ijtmp)
         else
            write(iot,'(a15,2x,"(",i5,"-th bin)"
     &                             /1pe13.4,2x,"$--$",1pe13.4)')
     &           cij(ijtmp),ij(ijtmp),fg(itmp),fg(itmp+1)
            itmp = itmp - ij(ijtmp) + nij(ijtmp)+1
         end if
        end if
       end do
       write(iot,'("Particle = ",a8)') chq(ip)

       write(iot,'("e:")')

      end if

*-----------------------------------------------------------------------

      if( inolg .eq. 1 .and. ( ittwo(m) .eq. 3 .or. ittwo(m) .eq. 7 )
     &    .and. cmin .gt. 0.0 .and. cmax .gt. cmin ) then

         write(iot,'(/"#",78("-"))')
         write(iot,'("z: xorg[-1.03/0.05]"/
     &               "p: form[c1] xfac[c2] afac[c5] nosp notf")')

               write(iot,'( "p: xmin(",1p1g14.7,") xmax(",1p1g14.7,
     &                       ")")') axs1min, axs1max

               write(iot,'( "p: ymin(",1p1g14.7,") ymax(",1p1g14.7,
     &                       ")")') axs2min, axs2max

      end if

*-----------------------------------------------------------------------

          end do   ! ip = 1, np
          end do   ! ianataldata

          deallocate( fg )
          deallocate( fgaxs )
          deallocate( fgaxs2 )
          deallocate( fgaxs3 )
          deallocate( anataldata )
          deallocate( delvol )

*-----------------------------------------------------------------------

         end if          ! itaxis

*-----------------------------------------------------------------------

        close(iot)

        if( iteps(m) .ne. 0 ) then

          open(iot, file = fname, status = 'unknown' )
          call a_angel(idasa,fname)

        end if

*-----------------------------------------------------------------------

  900 continue     ! ioe, iax

*-----------------------------------------------------------------------
      deallocate( ixyz )

      deallocate (anatalrst)
      deallocate (vl,rdata,val,dlr)
      
      return
      end subroutine anatal_pdepstreg


************************************************************************
*                                                                      *
      subroutine anatal_pdepstrz(m,np,nr,nz,ne,nt,rm,zm,eb,tb,
     &                   tr,nfile,weightRate,idasa,manatally)
*                                                                      *
*       anatally calculation and output of t-deposit with r-z mesh     *
*                                                                      *
************************************************************************
      use sumtallymod , only : isumtally
      use sangelmod, only: itsans
      use anatallymod, only: anataldata,fg,fgaxs,fgaxs2,fgaxs3,delvol,
     &     cx_txt,lcx_txt,cy_txt,lcy_txt,cz_txt,lcz_txt
      use partmod ! frtati 2021/10/05

*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

      common /stat / istdev, irestart, ireschk

*-----------------------------------------------------------------------

      include 'param.inc'

*-----------------------------------------------------------------------

      double precision resc2   ! W
      double precision resc3   ! N
      common /restart/ resc2(itlmax), resc3(itlmax)

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
      common /tall10/ itaxn(itlmax), itaxs(itlmax,6)
      common /tall11/ itfln(itlmax), itfll(itlmax,6), ctfln(itlmax,6)
      character ctfln*100

      common /tall14/ ittll(itlmax), ittle(itlmax)
      character ittle*80
      common /tall15/ itanl(itlmax), itang(itlmax)
      character itang*200

      common /tall20/ ittty(itlmax), ittnm(itlmax), ittrg(itlmax),
     &                rttmi(itlmax), rttma(itlmax), rttdl(itlmax)
      common /tall27/ itaxl(itlmax), itaxt(itlmax),
     &                itayl(itlmax), itayt(itlmax),
     &                itazl(itlmax), itazt(itlmax)
      character itaxt*200, itayt*200, itazt*200

      common /tall21/ rtfac(itlmax)
      common /tall35/ iteps(itlmax)

      common /istcut/ ist_cut, ist_bat
*-----------------------------------------------------------------------

      common /tcntl/  icntl, inucr

*-----------------------------------------------------------------------

      common /cparm/  maxbch,maxcas
      common /jcomon/ nabov,nobch,nocas,nomax
!$OMP THREADPRIVATE(/jcomon/)
!$    integer ipomp,npomp
!$    integer OMP_GET_NUM_THREADS,OMP_GET_THREAD_NUM

      common /talout/ itall

      character fname*100, fnume*3

*-----------------------------------------------------------------------
      character gfnam*100
      common /tall60/ itger(itlmax)

      common /tall67/ itstd(itlmax), rtstd(itlmax)

*-----------------------------------------------------------------------

      dimension   rm(nr+1)
      dimension   zm(nz+1)
      dimension   eb(ne+1)
      dimension   tb(nt+1)
!      dimension   tw(nt)
      dimension   tr(np,0:ne,nt,nr,nz,2*nfile)
      dimension   trunc(np,0:ne,nt,nr,nz,5)
      integer irst,nrst, itmpdata
!      dimension   anatalrst(np,ne+1,1,nt+1,nr*nz,1,2*nfile+3)
! sumover
!      dimension   anatalrst(np,ne+1,1,nt+1,(nr+1)*(nz+1)*2,
!     &                      1,2*nfile+3)
      real(8),allocatable :: anatalrst(:,:,:,:,:,:,:)

      integer ianataldata,nanataldata
      integer ntaxis
      parameter (ntaxis=9) ! S.H. 2021.8.9
      integer ij(ntaxis),nij(ntaxis),ibin(ntaxis)
      integer nijaxs,nijaxs2,nijaxs3,nfgmax
      character cij(ntaxis)*15,cijaxs*15,cijaxs2*15
      integer iDaxis
!      dimension   weightRate(nfile), rdata(2,nfile)
      dimension   weightRate(nfile)

      real(8),allocatable :: tw(:),rdata(:,:)

*-----------------------------------------------------------------------
      real(8) cfac
      real(8),parameter:: c2gy=1.602d-10,c2jperm3=1.602d-7

      dimension tott(np,2) ! frtati 2021/10/05 6 -> np

      character hsunit(0:5)*22
      data hsunit( 0) / 'Dose [Gy/source]      '/
      data hsunit( 1) / 'Dose [MeV/cm^3/source]'/
      data hsunit( 2) / 'Dose [MeV/source]     '/
      data hsunit( 3) / 'Number [1/source]     '/
      data hsunit( 4) / 'Number [1/nsec/source]'/
      data hsunit( 5) / 'Dose [J/m^3/source]   '/

      character cha*1
      data cha /"'"/

      character dc2*4
      character chp(itmxpt)*11  ! kitamura22/03/31
      character chq(itmxpt)*9   ! kitamura22/03/31

*-----------------------------------------------------------------------

      data pi/3.14159265d+0/

      character rpa*1
      data rpa /'}'/
      character yen*1

      dimension dt_one(1)
      data dt_one/1.0d0/

*-----------------------------------------------------------------------
      real(8),allocatable :: vl(:,:,:) ! S.H. 2021.8.17

      real(8),allocatable :: vl_r(:,:),vl_z(:,:)

*-----------------------------------------------------------------------
      integer           iat, iad, iaf

      icf(ir,iz,ia) = ir + ( iz - 1 ) * (nr+1)
     &              + ( ia - 1 ) * (nr+1) * (nz+1)

      iat(iad,iaf) = iad + (iaf-1) * 2
*-----------------------------------------------------------------------

      yen  = char(92)
      igsh = 0

      allocate (tw(nt),rdata(2,nfile))

      nrst = 2*nfile+3
      allocate (anatalrst(np,ne+1,1,nt+1,(nr+1)*(nz+1)*1,1,nrst))


*-----------------------------------------------------------------------
      np_mxang = np
      if ( itmxang(m).ne.0 .and. iabs(itmxang(m)).lt.np ) then
        np_mxang = iabs(itmxang(m))
      end if

*-----------------------------------------------------------------------

         if( itout(m) .le. 1 ) then

            nei = 0

         else

            nei = 1

         end if

*-----------------------------------------------------------------------
*        particle name
*-----------------------------------------------------------------------

               call wtpname(m,np,chp,chq)

*-----------------------------------------------------------------------
*        set mesh volume ( unit = 2, 3, 4 ; vol = 1.0 )
*-----------------------------------------------------------------------
               allocate( vl_r(nz,1),vl_z(nr,1) )

               allocate( vl(nr,nz,1) )
               if( itunt(m).ge.2 .and. itunt(m).le.4 )then
                vl(1:nr,1:nz,1)=1.0d0
                vl_r(:,1)=1.0d0
                vl_z(:,1)=1.0d0
               else
                vl_r(:,1)=0.0d0
                vl_z(:,1)=0.0d0
                do iz=1,nz
                 do ir=1,nr
                  vl(ir,iz,1)=pi * ( rm(ir+1)**2 - rm(ir)**2 )
     &                 * ( zm(iz+1) - zm(iz) )
                  vl_r(iz,1) = vl_r(iz,1) + vl(ir,iz,1)
                  vl_z(ir,1) = vl_z(ir,1) + vl(ir,iz,1)
                 enddo
                enddo
               endif

*-----------------------------------------------------------------------
*           itunt(m) = 0 : cfac = c2gy
*           itunt(m) = 5 ; cfac = c2jperm3
*-----------------------------------------------------------------------
            if( itunt(m) .eq. 0 )then
              cfac=c2gy
            elseif( itunt(m) .eq. 5)then
              cfac = c2jperm3
            else
              cfac=1.0d0
            endif

*-----------------------------------------------------------------------
*           itunt(m) = 4 : /nsec/source
*-----------------------------------------------------------------------

            if( itunt(m) .eq. 4 ) then

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
*        c1 : normalization for source
*-----------------------------------------------------------------------

             if( isumtally(m)/=2 .or. itout(m)==1 ) then
              if ( rsouin .gt. 0d0 ) then ! S.H. (2016.7.26)
               c1 = 1.0d+0 / rsouin
              else
               c1 = 0d0
              end if
             else
               c1 = 1.0d+0
             end if

*-----------------------------------------------------------------------
*           itunt(m) = 0 : Gy/source
*                    = 1 : MeV/cm^3/source
*                    = 2 : MeV/source
*                    = 3 : 1/source
*                    = 4 : 1/nsec/source
*                    = 5 : J/m^3/source
*-----------------------------------------------------------------------
*        relative error and  unit conversion
*-----------------------------------------------------------------------

                  cmax = 0.0
                  cmin = 1.e+33
                  isdz = 0
                  stdm = 0.0
                  itstd(m) = 1

        do itmprz=1,nr*nz
            ir=(itmprz-1)/nz+1
            iz=itmprz-(ir-1)*nz
            do 100 it = 1, nt
            do 100 ie = nei, ne
            do 100 ip = 1, np
            do ntf = 1, nfile

             if( isumtally(m)/=2 .or. itout(m)==1 ) then   ! Weighted average

               if( tr(ip,ie,it,ir,iz,iat(1,ntf)) .gt. 0.d0 ) then

                  call calc_deposit_stdev(m,Xa,sigx,
     &                            tr(ip,ie,it,ir,iz,iat(1,ntf)),
     &                            tr(ip,ie,it,ir,iz,iat(2,ntf)),
     &                            rtfac(m)/vl(ir,iz,1)/tw(it),ip)   ! S.Abe 2015/12/03, add ip
                  tr(ip,ie,it,ir,iz,iat(1,ntf)) = cfac*Xa
                  tr(ip,ie,it,ir,iz,iat(2,ntf)) = sigx
                  if( sigx .gt. stdm ) stdm = sigx

                  if( tr(ip,ie,it,ir,iz,iat(1,ntf)) .gt. cmax )
     &         cmax = tr(ip,ie,it,ir,iz,iat(1,ntf))

                  if( tr(ip,ie,it,ir,iz,iat(1,ntf)) .lt. cmin )
     &         cmin = tr(ip,ie,it,ir,iz,iat(1,ntf))

               else

                  isdz = 1
                  tr(ip,ie,it,ir,iz,iat(2,ntf)) = 0.0

               end if

             else
               if( tr(ip,ie,it,ir,iz,iat(1,ntf)) .gt. 0.d0 ) then
                 if( tr(ip,ie,it,ir,iz,iat(1,ntf)) .gt. cmax
     &              .and. manatally .ne. 0 )  ! T.Sato 2021/01/30
     &        cmax = tr(ip,ie,it,ir,iz,iat(1,ntf))

                 if( tr(ip,ie,it,ir,iz,iat(1,ntf)) .lt. cmin
     &            .and. manatally .ne. 0 )  ! T.Sato 2021/01/30
     &        cmin = tr(ip,ie,it,ir,iz,iat(1,ntf))
               else
                 isdz = 1
                 tr(ip,ie,it,ir,iz,iat(2,ntf)) = 0.0
               end if
             end if

! sumover
              call pdepstrz_sumover_stdev_ntf(0,m,ntf,
     &             ip,ie,it,ir,iz,
     &             nei,
     &             rtfac(m),cfac,tw(it),vl(ir,iz,1),
     &             tw_sum,vl_r(iz,1),vl_z(ir,1))

            end do   ! loop end ntf
  100       continue
      enddo  ! rz loop T.Sato 2021/04/15

!$OMP PARALLEL
!$OMP& private(ipomp,ir,iz,it,ie,ip,ntf,itmprz,rdata,answer,rerr)
!$OMP& private(iax,ir_a,iz_a,it_a,ie_a,ido_ana,ioe)
!$      npomp=OMP_GET_NUM_THREADS()
!$      ipomp=OMP_GET_THREAD_NUM()
!$      write(*,'(''OpenMP PARALLEL PROCESS'',
!$   &     i4,''/'',i4)') ipOMP+1,npomp
!$OMP DO SCHEDULE(dynamic) reduction(max:cmax) reduction(min:cmin)
        do itmprz=1,nr*nz
            ir=(itmprz-1)/nz+1
            iz=itmprz-(ir-1)*nz
            do 101 it = 1, nt
            do 101 ie = nei, ne
            do 101 ip = 1, np

            if( manatally .eq. 0 ) then ! T.Sato 2020/10/19

              do ntf = 1, nfile
                rdata(1,ntf)=tr(ip,ie,it,ir,iz,iat(1,ntf))
                rdata(2,ntf)=tr(ip,ie,it,ir,iz,iat(2,ntf))
              end do

              call usranatal(nfile,rdata,answer,rerr)

              tr(ip,ie,it,ir,iz,1) = answer
              tr(ip,ie,it,ir,iz,2) = rerr
             if ( itout(m) .le. 1 ) then ! ie=nei=ne=0
              anatalrst(ip,1,1,it,icf(ir,iz,1),1,1) = answer
              anatalrst(ip,1,1,it,icf(ir,iz,1),1,2) = rerr
             else
              anatalrst(ip,ie,1,it,icf(ir,iz,1),1,1) = answer
              anatalrst(ip,ie,1,it,icf(ir,iz,1),1,2) = rerr
             end if

              if(answer.gt.cmax) cmax = answer
              if(answer.gt.0d0 .and. answer.lt.cmin) cmin = answer

            else if( manatally .eq. 1 ) then

               do ntf = 1, nfile
                  rdata(1,ntf)=tr(ip,ie,it,ir,iz,iat(1,ntf))
                  rdata(2,ntf)=tr(ip,ie,it,ir,iz,iat(2,ntf))
               end do

               call anova(nfile,rdata,resc3(m)
     &              ,fmval,unca,uncerr,uratio)

               trunc(ip,ie,it,ir,iz,1) ! mean
     &              = fmval
               if ( fmval .gt. 0d0 ) then
                  trunc(ip,ie,it,ir,iz,2) ! total uncertainty
     &                 = dsqrt(unca**2 + uncerr**2/resc3(m)) / fmval
                  trunc(ip,ie,it,ir,iz,3) ! syst. uncertainty
     &                 = unca / fmval
                  trunc(ip,ie,it,ir,iz,4) ! stat. uncertainty
     &                 = uncerr /dsqrt(resc3(m)) / fmval
                  trunc(ip,ie,it,ir,iz,5) ! squared unc. ratio
     &                 = uratio
               else
                  trunc(ip,ie,it,ir,iz,2) ! total uncertainty
     &                 = 0d0
                  trunc(ip,ie,it,ir,iz,3) ! syst. uncertainty
     &                 = 0d0
                  trunc(ip,ie,it,ir,iz,4) ! stat. uncertainty
     &                 = 0d0
                  trunc(ip,ie,it,ir,iz,5) ! squared unc. ratio
     &                 = 0d0
               end if

              if ( itout(m) .le. 1 ) then ! ie=nei=ne=0
               anatalrst(ip,1,1,it,icf(ir,iz,1),1,1) ! mean
     &              = fmval
               if ( fmval .gt. 0d0 ) then
                  anatalrst(ip,1,1,it,icf(ir,iz,1),1,2) ! total uncertainty
     &                 = dsqrt(unca**2 + uncerr**2/resc3(m)) / fmval
                  anatalrst(ip,1,1,it,icf(ir,iz,1),1,3) ! syst. uncertainty
     &                 = unca / fmval
                  anatalrst(ip,1,1,it,icf(ir,iz,1),1,4) ! stat. uncertainty
     &                 = uncerr /dsqrt(resc3(m)) / fmval
                  anatalrst(ip,1,1,it,icf(ir,iz,1),1,5) ! squared unc. ratio
     &                 = uratio
               else
                  anatalrst(ip,1,1,it,icf(ir,iz,1),1,2) ! total uncertainty
     &                 = 0d0
                  anatalrst(ip,1,1,it,icf(ir,iz,1),1,3) ! syst. uncertainty
     &                 = 0d0
                  anatalrst(ip,1,1,it,icf(ir,iz,1),1,4) ! stat. uncertainty
     &                 = 0d0
                  anatalrst(ip,1,1,it,icf(ir,iz,1),1,5) ! squared unc. ratio
     &                 = 0d0
               end if
              else
               anatalrst(ip,ie,1,it,icf(ir,iz,1),1,1) ! mean
     &              = fmval
               if ( fmval .gt. 0d0 ) then
                  anatalrst(ip,ie,1,it,icf(ir,iz,1),1,2) ! total uncertainty
     &                 = dsqrt(unca**2 + uncerr**2/resc3(m)) / fmval
                  anatalrst(ip,ie,1,it,icf(ir,iz,1),1,3) ! syst. uncertainty
     &                 = unca / fmval
                  anatalrst(ip,ie,1,it,icf(ir,iz,1),1,4) ! stat. uncertainty
     &                 = uncerr /dsqrt(resc3(m)) / fmval
                  anatalrst(ip,ie,1,it,icf(ir,iz,1),1,5) ! squared unc. ratio
     &                 = uratio
               else
                  anatalrst(ip,ie,1,it,icf(ir,iz,1),1,2) ! total uncertainty
     &                 = 0d0
                  anatalrst(ip,ie,1,it,icf(ir,iz,1),1,3) ! syst. uncertainty
     &                 = 0d0
                  anatalrst(ip,ie,1,it,icf(ir,iz,1),1,4) ! stat. uncertainty
     &                 = 0d0
                  anatalrst(ip,ie,1,it,icf(ir,iz,1),1,5) ! squared unc. ratio
     &                 = 0d0
               end if
              end if

            else if( manatally .eq. 2 ) then

               do ntf = 1, nfile
                do ioe = 1, 2
                 if ( itout(m) .le. 1 ) then ! ie=nei=ne=0
                  anatalrst(ip,1,1,it,icf(ir,iz,1),1,iat(ioe,ntf))
     &              = tr(ip,ie,it,ir,iz,iat(ioe,ntf))
                 else
                  anatalrst(ip,ie,1,it,icf(ir,iz,1),1,iat(ioe,ntf))
     &              = tr(ip,ie,it,ir,iz,iat(ioe,ntf))
                 end if
                end do
               end do

            end if

!sumover
            do iax = 1,itaxn(m)

              ie_a = ie
              it_a = it
              ir_a = ir
              iz_a = iz
              ido_ana = 0
              if(itaxs(m,iax) == 1 .and. ie == 1) then
                ie_a = ne + 1
                ido_ana = 1
              elseif(itaxs(m,iax) == 5 .and. iz == 1) then
                iz_a = nz + 1
                ido_ana = 1
              elseif(itaxs(m,iax) == 6 .and. ir == 1) then
                ir_a = nr + 1
                ido_ana = 1
              elseif(itaxs(m,iax) == 11 .and. it == 1) then
                it_a = nt + 1
                ido_ana = 1
              endif

              if(ido_ana == 1) then

                 call get_depstrz_tr_sum_data(m,iax,nfile,
     &                 ip,ie,it,ir,iz,rdata)

                if( manatally .eq. 0 ) then ! T.Sato 2020/10/22

                   call usranatal(nfile,rdata,answer,rerr)

                   if ( itout(m) .le. 1 ) then ! ie=nei=ne=0
                     anatalrst(ip,1,1,it_a,icf(ir_a,iz_a,1),1,1)
     &               = answer
                     anatalrst(ip,1,1,it_a,icf(ir_a,iz_a,1),1,2)
     &               = rerr
                   else
                     anatalrst(ip,ie_a,1,it_a,icf(ir_a,iz_a,1),1,1)
     &               = answer
                     anatalrst(ip,ie_a,1,it_a,icf(ir_a,iz_a,1),1,2)
     &               = rerr
                   endif

                else if( manatally .eq. 1 ) then

                  call anova(nfile,rdata,resc3(m)
     &               ,fmval,unca,uncerr,uratio)

                  if ( itout(m) .le. 1 ) then ! ie=nei=ne=0

                    anatalrst(ip,1,1,it_a,icf(ir_a,iz_a,1),1,1) ! mean
     &                = fmval
                    if ( fmval .gt. 0.0d0 ) then
                      anatalrst(ip,1,1,it_a,icf(ir_a,iz_a,1),1,2) ! total uncertainty
     &                     = dsqrt(unca**2 + uncerr**2/resc3(m)) / fmval
                    else
                      anatalrst(ip,1,1,it_a,icf(ir_a,iz_a,1),1,2) ! total uncertainty
     &                   = 0.0d0
                    end if
                  else
                    anatalrst(ip,ie_a,1,it_a,icf(ir_a,iz_a,1),1,1) ! mean
     &                = fmval
                    if ( fmval .gt. 0.0d0 ) then
                      anatalrst(ip,ie_a,1,it_a,icf(ir_a,iz_a,1),1,2) ! total uncertainty
     &                     = dsqrt(unca**2 + uncerr**2/resc3(m)) / fmval
                    else
                      anatalrst(ip,ie_a,1,it_a,icf(ir_a,iz_a,1),1,2) ! total uncertainty
     &                   = 0.0d0
                    end if
                  endif

                else if( manatally .eq. 2 ) then
                  do ntf = 1, nfile
                    do ioe = 1, 2
                      if ( itout(m) .le. 1 ) then ! ie=nei=ne=0
                        anatalrst(ip,1,1,it_a,icf(ir_a,iz_a,1),
     &                            1,iat(ioe,ntf))
     &                  = rdata(ioe,ntf)
                      else
                        anatalrst(ip,ie_a,1,it_a,icf(ir_a,iz_a,1),
     &                            1,iat(ioe,ntf))
     &                  = rdata(ioe,ntf)
                      endif
                    end do
                  end do

                end if

              endif
            enddo


  101       continue
      enddo  ! rz loop T.Sato 2021/04/15
!$OMP END DO
!$OMP END PARALLEL
            if( istdev .eq. 2 .or. nobch .gt. 20 ) then
            if( rtstd(m) .gt. 0.d0 .and.
     &          stdm .gt. 0.d0 .and. isdz .eq. 0 .and.
     &          stdm .lt. rtstd(m) ) itstd(m) = 0
            end if

*-----------------------------------------------------------------------
*     out put unit = 31  : temporary number
*-----------------------------------------------------------------------

      do iax = 1, itfln(m)

!OBINATA(2012.8.6): output *.err
        noe = 1
        if ( any( itaxs(m,iax) .eq. (/ 10, 12, 13 /) )
     &      .and. ittwo(m) .ne. 4 .and. icntl .ne. 8 )
     &    noe = 2

        do ioe = 1, noe

         if( itall .eq. 2 .and. nobch .lt. maxbch ) then

            write(fnume,'(i3.3)') nobch
!OBINATA(2012.8.6): output *.err
            if ( ioe .eq. 1 ) then
              fname = ctfln(m,iax)(1:itfll(m,iax))//'.'//fnume
            else
              call mk_2derrfn(ctfln(m,iax),fname,itfll(m,iax))
            end if

         else

!OBINATA(2012.8.6): output *.err
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

         call depstech(iot,m,iax,1)

*-----------------------------------------------------------------------
* manatally
*  0  :  user defined analysis
*  1  :  systematic uncertainty
*  2  :  c-value dependence
*-----------------------------------------------------------------------

*-----------------------------------------------------------------------
*        energy, r, z, time axes
*-----------------------------------------------------------------------

         if( itaxs(m,iax) .eq. 1
     &        .or. itaxs(m,iax) .eq. 6
     &        .or. itaxs(m,iax) .eq. 5
     &        .or. itaxs(m,iax) .eq. 11 ) then

          select case( itaxs(m,iax) )
           case ( 1 )
            iDaxis = 1   ! energy
           case ( 6 )
            iDaxis = 6   ! r axis
           case ( 5 )
            iDaxis = 7   ! z axis
           case ( 11 )
            iDaxis = 9   ! time
          end select
          if( manatally.eq.2 ) iDaxis = 51 ! S.H. 2021.8.19

C rearrange data from anatalrst to anataldata
          call anatal_rearrange_sum(0,0,2,
     &         iDaxis,ntaxis,anatalrst,
     &         np,     ne,      0,     nt,     nr, nz,  0,     0,nrst,
     &         eb,  1.0d0,  1,  1, tb, tw, rm, zm,  dt_one,  1,
     &         vl,
     &         nij,ibin,nijaxs,nijaxs2,nijaxs3,nfgmax,
     &         cij,cijaxs,cijaxs2,nanataldata,ierr)

          inum = 0

          do ianataldata = 1, nanataldata

           inum = inum + 1

C exchange ID number from ianataldata to ij2,...,ij8
           itmpdata = int(ianataldata-1)
           ij(8) = mod(itmpdata,nij(8))+1
           itmpdata = int((ianataldata-ij(8))/nij(8))
           ij(7) = mod(itmpdata,nij(7))+1
           itmpdata = int((itmpdata-ij(7)+1)/nij(7))
           ij(6) = mod(itmpdata,nij(6))+1
           itmpdata = int((itmpdata-ij(6)+1)/nij(6))
           ij(5) = mod(itmpdata,nij(5))+1
           itmpdata = int((itmpdata-ij(5)+1)/nij(5))
           ij(4) = mod(itmpdata,nij(4))+1
           itmpdata = int((itmpdata-ij(4)+1)/nij(4))
           ij(3) = mod(itmpdata,nij(3))+1
           itmpdata = int((itmpdata-ij(3)+1)/nij(3))
           ij(2) = mod(itmpdata,nij(2))+1

*-----------------------------------------------------------------------

           write(iot,'(/"#",78("-"))')

           if( inum .eq. 1 ) then
              write(iot,'( "#newpage:")')
           else
              write(iot,'( " newpage:")')
           end if

           if( lcx_txt(m) .gt. 0 ) then
              write(iot,'(/"x: ",200a1)')(cx_txt(m)(i:i),i=1,lcx_txt(m))
           else
            if( itaxl(m) .eq. 0 ) then
               write(iot,'(/"x: ",a15)') cijaxs
            else
               write(iot,'(/"x: ",200a1)') (itaxt(m)(i:i),i=1,itaxl(m))
            end if
           end if

           if( lcy_txt(m) .gt. 0 ) then
              write(iot,'( "y: ",200a1)')(cy_txt(m)(i:i),i=1,lcy_txt(m))
           else
            if( itayl(m) .eq. 0 ) then
               write(iot,'( "y: ",a22)') hsunit(itunt(m))
            else
               write(iot,'( "y: ",200a1)') (itayt(m)(i:i),i=1,itayl(m))
            end if
           end if

           if( itaxs(m,iax).eq.1 .and.
     &          ( itety(m) .eq. 3 .or. itety(m) .eq. 5 ) ) then
              write(iot,'( "p: xlog ylog afac(0.8) form(0.9)")')

           else if( itaxs(m,iax).eq.6 .and.
     &             ( itrty(m) .eq. 3 .or. itrty(m) .eq. 5 ) ) then
              write(iot,'( "p: xlog ylog afac(0.8) form(0.9)")')

           else if( itaxs(m,iax).eq.5 .and.
     &             ( itzty(m) .eq. 3 .or. itzty(m) .eq. 5 ) ) then
              write(iot,'( "p: xlog ylog afac(0.8) form(0.9)")')

           else if( itaxs(m,iax).eq.11 .and.
     &             ( ittty(m) .eq. 3 .or. ittty(m) .eq. 5 ) ) then
              write(iot,'( "p: xlog ylog afac(0.8) form(0.9)")')

           else
              write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')

           end if

           if( itanl(m) .gt. 0 ) then
              write(iot,'( "p: ",200a1)')
     &             ( itang(m)(i:i),i = 1, itanl(m) )
           end if

           if( itsans(m) .gt. 0 ) then
              call write_sangel(iot,m,0)
           end if

           if ( manatally .eq. 0 ) then ! user defined analysis
            if ( iteps(m) .ne. 2 ) then
               write(iot,'( "h: n",12x,"x",12x,
     &                    1000(a10,"),hh0",a3," n "))')
     &                     ( chp(i), chl(i) , i = 1, np_mxang ) ! frtati 2021/10/05 np -> np_mxang

            else
               write(iot,'( "h: n",12x,"x",12x,
     &         1000(a1,i1,a9,"),hh0",a3
     &         ," ny",i0," dy",i0,"=[y",i0,"*y",i0,"] "))')
     &         ( chp(i)(1:1),i,chp(i)(2:10), chl(i)
     &         ,i+np_mxang,i,i,i+np_mxang, i = 1, np_mxang ) ! frtati 2021/10/05

            end if

            write(iot,'( "#  lower        upper  ",3x,
     &                   1000(a1,"  ",a8,"    r.err "))')
     &                    ( chb(i), chp(i)(3:10), i = 1, np ) ! frtati 2021/10/05

           else if ( manatally .eq. 1 ) then ! systematic uncertainty
            if ( iteps(m) .ne. 2 ) then
               write(iot,'( "h: n",12x,"x",12x,
     &                    1000(a10,"),hh0",a3," n n n "))')
     &                     ( chp(i), chl(i) , i = 1, np_mxang ) ! frtati 2021/10/05 np -> np_mxang

            else
               write(iot,'( "h: n",12x,"x",12x,
     &         1000(a1,i1,a9,"),hh0",a3
     &         ," ny",i0," dy",i0,"=[y",i0,"*y",i0,"] n n "))')
     &         ( chp(i)(1:1),i,chp(i)(2:10), chl(i)
     &         ,i+np_mxang,i,i,i+np_mxang, i = 1, np_mxang ) ! frtati 2021/10/05

            end if

            write(iot,'( "#  lower        upper  ",3x,
     &                   1000(a1,2x,a8,4x,"r.err(tot, syst, stat)"))')
     &                    ( chb(i), chp(i)(3:10), i = 1, np ) ! frtati 2021/10/05

           else if ( manatally .eq. 2 ) then ! c-value dependence
            if ( iteps(m) .ne. 2 ) then
               write(iot,'( "h: n",12x,"x",12x,
     &                    1000(a10,"),hh0",a3," n "))')
     &                     ( chp(i), chl(i) , i = 1, np_mxang ) ! frtati 2021/10/05 np -> np_mxang

            else
               write(iot,'( "h: n",12x,"x",12x,
     &         1000(a1,i1,a9,"),hh0",a3
     &         ," ny",i0," dy",i0,"=[y",i0,"*y",i0,"] "))')
     &         ( chp(i)(1:1),i,chp(i)(2:10), chl(i)
     &         ,i+np_mxang,i,i,i+np_mxang, i = 1, np_mxang ) ! frtati 2021/10/05

            end if

            write(iot,'( "#  lower        upper  ",3x,
     &                   1000(a1,2x,a8,4x,"r.err "))')
     &                    ( chb(i), chp(i)(3:10), i = 1, np ) ! frtati 2021/10/05

           end if

           do ip = 1, np
              tott(ip,1) = 0.d+0
           end do

           voll = 0.0d0
           njaxs = nijaxs + 1

           do ijaxs = 1, nijaxs

            if ( manatally .eq. 0 ) then ! user defined analysis
               write(iot,'(1p2e13.4,1000(1pe13.4,0pf8.4))')
     &              fgaxs(ijaxs),fgaxs(ijaxs+1),
     &              ((anataldata(ip,ianataldata,ijaxs,irst),irst=1,2)
     &              ,ip=1,np) ! frtati 2021/10/05

            else if ( manatally .eq. 1 ) then ! systematic uncertainty
               write(iot,'(1p2e13.4,1000(1pe13.4,0p3f8.4))')
     &              fgaxs(ijaxs),fgaxs(ijaxs+1),
     &              ((anataldata(ip,ianataldata,ijaxs,irst),irst=1,4)
     &              ,ip=1,np) ! frtati 2021/10/05

            else if ( manatally .eq. 2 ) then ! c-value dependence
               write(iot,'(1p2e13.4,1000(1pe13.4,0pf8.4))')
     &              fgaxs(ijaxs),fgaxs(ijaxs+1),
     &              ((anataldata(ip,ianataldata,ijaxs,irst),irst=1,2)
     &              ,ip=1,np) ! frtati 2021/10/05

            end if

            voll = voll + delvol(ianataldata,ijaxs)

            do ip = 1, np
               vn = anataldata(ip,ianataldata,ijaxs,1)
     &              * delvol(ianataldata,ijaxs)
               tott(ip,1) = tott(ip,1) + vn
            end do

           end do

           if( itaxs(m,iax) .eq. 1 ) then ! when axis = eng
            do ip = 1, np
               tott(ip,1) = tott(ip,1)
            end do

           else if( itaxs(m,iax).eq.6 .or. itaxs(m,iax).eq.5 ) then ! when axis = r,z
            do ip = 1, np
             if( itunt(m) .eq. 1 .or. itunt(m) .eq. 5 ) then
                tott(ip,1) = tott(ip,1) / voll
             end if
            end do

           else if( itaxs(m,iax) .eq. 11 ) then ! when axis = t
            do ip = 1, np
             if( itunt(m) .eq. 4 ) then
                tott(ip,1) = tott(ip,1) / voll
             end if
            end do

           end if

! sumover
            do ip = 1, np
            do irst=1,2
              tott(ip,irst) = anataldata(ip,ianataldata,njaxs,irst)
            enddo
            enddo

           if( itaxs(m,iax).eq.1 .or. itaxs(m,iax).eq.11 ) then ! when axis = eng,t
            if( itunt(m) .ne. 0 ) then
             if ( manatally .eq. 0 ) then ! user defined analysis
                write(iot,'(/"#   sum over",14x,1000(1pe13.4,0pf8.4))')
     &               (tott(ip,1),tott(ip,2),ip=1,np) ! frtati 2021/10/05
             else if ( manatally .eq. 1 ) then ! systematic uncertainty
             write(iot,'(/"#   sum over",14x,1000(1pe13.4,0pf8.4,16x))')
     &                     (tott(ip,1),tott(ip,2),ip=1,np)
             else if ( manatally .eq. 2 ) then ! c-value dependence
                write(iot,'(/"#   sum over",14x,1000(1pe13.4,0pf8.4))')
     &               (tott(ip,1),tott(ip,2),ip=1,np) ! frtati 2021/10/05
             end if
            end if

           else if( itaxs(m,iax).eq.6 .or. itaxs(m,iax).eq.5 ) then ! when axis = r,z
            if( itunt(m) .eq. 1 .or. itunt(m) .eq. 5 ) then
             if ( manatally .eq. 0 ) then ! user defined analysis
                write(iot,'(/"#   sum over ",1pe13.4,
     &                     1000(1pe13.4,0pf8.4))')
     &                     voll, (tott(ip,1),tott(ip,2),ip=1,np) ! frtati 2021/10/05
             else if ( manatally .eq. 1 ) then ! systematic uncertainty
                write(iot,'(/"#   sum over ",1pe13.4,
     &                     1000(1pe13.4,0pf8.4,16x))')
     &                     voll, (tott(ip,1),tott(ip,2),ip=1,np) ! frtati 2021/10/05
             else if ( manatally .eq. 2 ) then ! c-value dependence
                write(iot,'(/"#   sum over ",1pe13.4,
     &                     1000(1pe13.4,0pf8.4))')
     &                     voll, (tott(ip,1),tott(ip,2),ip=1,np) ! frtati 2021/10/05
             end if

            else if( itunt(m) .ne. 0 ) then
             if ( manatally .eq. 0 ) then ! user defined analysis
                write(iot,'(/"#   sum over ",13x,
     &                     1000(1pe13.4,0pf8.4))')
     &                     (tott(ip,1),tott(ip,2),ip=1,np) ! frtati 2021/10/05
             else if ( manatally .eq. 1 ) then ! systematic uncertainty
                write(iot,'(/"#   sum over ",13x,
     &                     1000(1pe13.4,0pf8.4,16x))')
     &                     (tott(ip,1),tott(ip,2),ip=1,np) ! frtati 2021/10/05
             else if ( manatally .eq. 2 ) then ! c-value dependence
                write(iot,'(/"#   sum over ",13x,
     &                     1000(1pe13.4,0pf8.4))')
     &                     (tott(ip,1),tott(ip,2),ip=1,np) ! frtati 2021/10/05
             end if
            end if

           end if

*-----------------------------------------------------------------------

           write(iot,'(/a1,"no. =",i3,a1)') cha, inum, cha

           write(iot,'("msuc: {",a1,"huge ",80a1)')
     &          yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

           write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
           write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

C write tally condition of each figure
           write(iot,'("wt: s(0.7)",/a1,"vspace{-3}")') yen
           itmp = 0
           do ijtmp = 2,8
            if (cij(ijtmp) .ne. 'F' ) then
             itmp = itmp + ij(ijtmp)
             if ( ibin(ijtmp) .eq. 0 ) then
                write(iot,'(a15,"&=&",i5)')
     &               cij(ijtmp), idnint(fg(itmp))
                itmp = itmp - ij(ijtmp) + nij(ijtmp)
             else
                write(iot,'(a15,2x,"(",i5,"-th bin)"
     &                             /1pe13.4,2x,"$--$",1pe13.4)')
     &               cij(ijtmp),ij(ijtmp),fg(itmp),fg(itmp+1)
                itmp = itmp - ij(ijtmp) + nij(ijtmp)+1
             end if
            end if
           end do

           write(iot,'("e:")')

          end do                  ! ianataldata = 1, nanataldata

          deallocate( fg )
          deallocate( fgaxs )
          deallocate( anataldata )
          deallocate( delvol )

*-----------------------------------------------------------------------
*        rz , t-eng, eng-t axes (matrix)
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 10
     &           .or. itaxs(m,iax) .eq. 12
     &           .or. itaxs(m,iax) .eq. 13 ) then

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

          select case( itaxs(m,iax) )
           case ( 10 )
            iDaxis = 34   ! rz
           case ( 12 )
            iDaxis = 35   ! t-eng
           case ( 13 )
            iDaxis = 36   ! eng-t
          end select

C rearrange data from anatalrst to anataldata
          call anatal_rearrange_sum(0,0,2,
     &         iDaxis,ntaxis,anatalrst,
     &         np,     ne,      0,     nt,     nr, nz,  0,     0,nrst,
     &         eb,  1.0d0,  1,  1, tb, tw, rm, zm,  dt_one,  1,
     &         vl,
     &         nij,ibin,nijaxs,nijaxs2,nijaxs3,nfgmax,
     &         cij,cijaxs,cijaxs2,nanataldata,ierr)

          inum = 0

          do ip = 1, np
          do ianataldata = 1, nanataldata

           inum = inum + 1

C exchange ID number from ianataldata to ij2,...,ij8
           itmpdata = int(ianataldata-1)
           ij(8) = mod(itmpdata,nij(8))+1
           itmpdata = int((ianataldata-ij(8))/nij(8))
           ij(7) = mod(itmpdata,nij(7))+1
           itmpdata = int((itmpdata-ij(7)+1)/nij(7))
           ij(6) = mod(itmpdata,nij(6))+1
           itmpdata = int((itmpdata-ij(6)+1)/nij(6))
           ij(5) = mod(itmpdata,nij(5))+1
           itmpdata = int((itmpdata-ij(5)+1)/nij(5))
           ij(4) = mod(itmpdata,nij(4))+1
           itmpdata = int((itmpdata-ij(4)+1)/nij(4))
           ij(3) = mod(itmpdata,nij(3))+1
           itmpdata = int((itmpdata-ij(3)+1)/nij(3))
           ij(2) = mod(itmpdata,nij(2))+1

           write(iot,'(/"#",78("-"))')

           if( inum .eq. 1 .or. ip.gt.np_mxang ) then ! frtati 2021/10/05
              write(iot,'( "#newpage:")')
           else
              write(iot,'( " newpage:")')
           end if

           if ( ip.gt.np_mxang ) then
             write(iot,'( " SKIPPAGE:")')
             inum = inum - 1
           end if

*-----------------------------------------------------------------------

           if( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) then

              write(iot,'(/a1,"no. =",i3,a1)') cha, inum, cha

              write(iot,'("msuc: {",a1,"huge ",80a1)')
     &             yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

              write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
              write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

           if( lcx_txt(m) .gt. 0 ) then
              write(iot,'(/"x: ",200a1)')(cx_txt(m)(i:i),i=1,lcx_txt(m))
           else
              if( itaxl(m) .eq. 0 ) then
                 write(iot,'(/"x: ",a15)') cijaxs
              else
                 write(iot,'(/"x: ",200a1)')
     &                (itaxt(m)(i:i),i=1,itaxl(m))
              end if
           end if

           if( lcy_txt(m) .gt. 0 ) then
              write(iot,'( "y: ",200a1)')(cy_txt(m)(i:i),i=1,lcy_txt(m))
           else
              if( itayl(m) .eq. 0 ) then
                 write(iot,'("y: ",a15)') cijaxs2
              else
                 write(iot,'( "y: ",200a1)')
     &                (itayt(m)(i:i),i=1,itayl(m))
              end if
           end if

*-----------------------------------------------------------------------

              if( inum .eq. 1 ) then

                 form  = ( fgaxs2(nijaxs2+1) - fgaxs2(1) )
     &                / ( fgaxs(nijaxs+1) - fgaxs(1) )
                 xfac  = 0.9
                 afac  = 0.8
                 izlog = 1
                 inocm = 1
                 inolg = 1

                 if( itanl(m) .gt. 0 )
     &                call anset(itang(m),itanl(m),
     &                form,xfac,afac,cmin,cmax,izlog,inocm,inolg)

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
     &             ( ittwo(m) .ge. 6 .and. ittwo(m) .le. 7 ) ) .and.
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
     &                write(iot,'( "p: zlog")')

              end if

*-----------------------------------------------------------------------

              axs1min = fgaxs(1)
              axs1max = fgaxs(nijaxs+1)
              axs2min = fgaxs2(1)
              axs2max = fgaxs2(nijaxs2+1)
              axs1del = (fgaxs(nijaxs+1)-fgaxs(1))/dble(nijaxs)
              axs2del = (fgaxs2(nijaxs2+1)-fgaxs2(1))/dble(nijaxs2)

              write(iot,'( "p: xmin(",1p1g14.7,") xmax(",1p1g14.7,
     &                       ")")') axs1min, axs1max

              write(iot,'( "p: ymin(",1p1g14.7,") ymax(",1p1g14.7,
     &                       ")")') axs2min, axs2max

              if( itanl(m) .gt. 0 .and. ioe .eq. 1 ) then

                 write(iot,'( "p: ",200a1)')
     &                ( itang(m)(i:i),i = 1, itanl(m) )

              else if( itanl(m) .gt. 0 .and. ioe .ne. 1 ) then

                 call terrang(iot,m)

              end if

              if( itsans(m) .gt. 0 ) then
                 call write_sangel(iot,m,0)
              end if

           end if

*-----------------------------------------------------------------------

           write(iot,'("#  ny = ",i3,"   nx = ",i3)')
     &          nijaxs2, nijaxs

           if( ittwo(m) .ne. 4 ) then
              write(iot,'( "# ( ( data(i1,i2), i1 = 1, n1 ),",
     &                         " i2 = n2, 1, -1 )")')
           end if

*-----------------------------------------------------------------------

           if( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) then

              write(iot,'(/a4," y = ",1p1g14.7," to ",
     &            1p1g14.7," by ",1p1g14.7," ; x = ",
     &            1p1g14.7," to ",1p1g14.7," by ",
     &            1p1g14.7," ;")') dc2,
     &            axs2max - axs2del/2d0, axs2min + axs2del/2d0, axs2del,
     &            axs1min + axs1del/2d0, axs1max - axs1del/2d0, axs1del

              write(iot,'(1p10e11.3)')
     &       ( ( anataldata(ip,ianataldata,ijaxs+(ijaxs2-1)*nijaxs,ioe),
     &                ijaxs = 1, nijaxs ), ijaxs2 = nijaxs2, 1, -1 )

*-----------------------------------------------------------------------

           else if( ittwo(m) .eq. 4 ) then

              write(iot,'(/"# axis2      axis1    ",
     &                      "  number     r.err")')

              do ijaxs = 1, nijaxs
               do ijaxs2 = 1, nijaxs2

                  write(iot,'(1p3e11.3,0pf8.4)')
     &                 fgaxs2(ijaxs2)  + axs2del/2d0,
     &                 fgaxs(ijaxs)  + axs1del/2d0,
     &             anataldata(ip,ianataldata,ijaxs+(ijaxs2-1)*nijaxs,1),
     &             anataldata(ip,ianataldata,ijaxs+(ijaxs2-1)*nijaxs,2)

               end do
              end do

*-----------------------------------------------------------------------

           else if( ittwo(m) .eq. 5 ) then

              write(iot,'("#   t = ",1p1g14.7," to ",
     &            1p1g14.7," by ",1p1g14.7/
     &                        "#   e = ",1p1g14.7," to ",
     &            1p1g14.7," by ",1p1g14.7)')
     &            axs2min + axs2del/2d0, axs2max - axs2del/2d0, axs2del,
     &            axs1min + axs1del/2d0, axs1max - axs1del/2d0, axs1del

              write(iot,'(/4x,a3,4x,1p1000e11.3)')
     &          'ax2/ax1',( fgaxs(ijaxs) + axs1del/2d0, ijaxs=1,nijaxs )

              do ijaxs2 = nijaxs2, 1, -1

                 write(iot,'(1p1000e11.3)')
     &                fgaxs2(ijaxs2)  + axs2del/2d0,
     &         ( anataldata(ip,ianataldata,ijaxs+(ijaxs2-1)*nijaxs,ioe),
     &                ijaxs = 1, nijaxs )

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

      else

       if( lcz_txt(m) .gt. 0 ) then
          write(iot,'("y: ",200a1)')(cz_txt(m)(i:i),i=1,lcz_txt(m))
       else
        if( itazl(m) .eq. 0 ) then
            write(iot,'("y: ",a22)') hsunit(itunt(m))
        else
            write(iot,'("y: ",200a1)') (itazt(m)(i:i),i=1,itazl(m))
        end if
       end if

      end if

      if( izlog .eq. 0 ) write(iot,'( "p: ylin")')

      write(iot,'( "#",78("-"))')

      end if

*-----------------------------------------------------------------------

      if( inocm .eq. 1 ) then

C write tally condition of each figure
       write(iot,'("wt: s[c5]",/a1,"vspace{-3}")') yen
       itmp = 0
       do ijtmp = 2,8
        if (cij(ijtmp) .ne. 'F' ) then
         itmp = itmp + ij(ijtmp)
         if ( ibin(ijtmp) .eq. 0 ) then
            write(iot,'(a15,"&=&",i5)')
     &           cij(ijtmp), idnint(fg(itmp))
            itmp = itmp - ij(ijtmp) + nij(ijtmp)
         else
            write(iot,'(a15,2x,"(",i5,"-th bin)"
     &                             /1pe13.4,2x,"$--$",1pe13.4)')
     &           cij(ijtmp),ij(ijtmp),fg(itmp),fg(itmp+1)
            itmp = itmp - ij(ijtmp) + nij(ijtmp)+1
         end if
        end if
       end do
       write(iot,'("Particle = ",a8)') chq(ip)

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
     &                       ")")') axs1min, axs1max

               write(iot,'( "p: ymin(",1p1g14.7,") ymax(",1p1g14.7,
     &                       ")")') axs2min, axs2max

      end if

*-----------------------------------------------------------------------

          end do                  ! ianataldata = 1, nanataldata
          end do                  ! ip = 1, np

          deallocate( fg )
          deallocate( fgaxs )
          deallocate( fgaxs2 )
          deallocate( fgaxs3 )
          deallocate( anataldata )
          deallocate( delvol )

*-----------------------------------------------------------------------

         end if          ! itaxis

*-----------------------------------------------------------------------

         close(iot)

         if( iteps(m) .ne. 0 ) then

            open(iot, file = fname, status = 'unknown' )
            call a_angel(idasa,fname)

         end if

*-----------------------------------------------------------------------

      end do !ioe

      end do !iax

      deallocate( vl ) !FURUTA20190208
*-----------------------------------------------------------------------
      deallocate( vl_r,vl_z )
      deallocate (tw,rdata)

      return
      end subroutine anatal_pdepstrz


************************************************************************
*                                                                      *
      subroutine anatal_pdepstxyz(m,np,nl,lt,
     &                    nx,ny,nz,ne,nt,xm,ym,zm,eb,tb,tr,nfile,
     &                    weightRate,
     &                    igsh,idasa,manatally)
*                                                                      *
*       anatally calculation and output of t-deposit with xyz mesh     *
*                                                                      *
************************************************************************
      use sumtallymod , only : isumtally
      use sangelmod, only: itsans
      use anatallymod, only: anataldata,fg,fgaxs,fgaxs2,fgaxs3,delvol,
     &     cx_txt,lcx_txt,cy_txt,lcy_txt,cz_txt,lcz_txt
      use partmod ! frtati 2021/10/05

*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

      common /stat / istdev, irestart, ireschk

*-----------------------------------------------------------------------

      include 'param.inc'

*-----------------------------------------------------------------------

      double precision resc2   ! W
      double precision resc3   ! N
      common /restart/ resc2(itlmax), resc3(itlmax)

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

      common /tall20/ ittty(itlmax), ittnm(itlmax), ittrg(itlmax),
     &                rttmi(itlmax), rttma(itlmax), rttdl(itlmax)
      common /tall21/ rtfac(itlmax)
      common /tall28/ itgsh(itlmax)
      common /tall22/ itmcn(itlmax), itmtn(itlmax), itmtt(itlmax)
      common /tall34/ itvm(itlmax)
      common /tall35/ iteps(itlmax)
      common /tall38/ itres(itlmax)
      common /tall39/ rtwid(itlmax)

      common /tall46/ itmtr(itlmax,4), rtmtr(itlmax,13)
      common /tall49/ itglt(itlmax)

      common /istcut/ ist_cut, ist_bat
*-----------------------------------------------------------------------

      common /tcntl/  icntl, inucr

*-----------------------------------------------------------------------

      common /cparm/  maxbch,maxcas
      common /jcomon/ nabov,nobch,nocas,nomax
!$OMP THREADPRIVATE(/jcomon/)
!$    integer ipomp,npomp
!$    integer OMP_GET_NUM_THREADS,OMP_GET_THREAD_NUM

      common /talout/ itall

      character fname*100, fnume*3

*-----------------------------------------------------------------------
      character gfnam*100
      common /tall60/ itger(itlmax)

      common /tall67/ itstd(itlmax), rtstd(itlmax)

*-----------------------------------------------------------------------

      dimension   lt(nl)
      dimension   xm(nx+1)
      dimension   ym(ny+1)
      dimension   zm(nz+1)
      dimension   eb(ne+1)
      dimension   tb(nt+1)
!      dimension   tw(nt)
      dimension   tr(np,0:ne,nt,nx*ny*nz,2*nfile)
!      dimension   weightRate(nfile), rdata(2,nfile)
      dimension   weightRate(nfile)
      real(8),allocatable :: tw(:),rdata(:,:)

      integer irst,nrst, itmpdata
      real(8),allocatable:: anatalrst(:,:,:,:,:,:,:)
      integer ianataldata,nanataldata
      integer, parameter :: ntaxis=9  ! S.H. 2021.8.9
      integer ij(ntaxis),nij(ntaxis),ibin(ntaxis)
      integer nijaxs,nijaxs2,nijaxs3,nfgmax
      character cij(ntaxis)*32,cijaxs*32,cijaxs2*32
      integer iDaxis

      integer,allocatable :: ixyz(:)

*-----------------------------------------------------------------------
      real(8) cfac
      real(8),parameter:: c2gy=1.602d-10,c2jperm3=1.602d-7

      dimension tott(np,2) ! frtati 2021/10/05 6 -> np

      character hsunit(0:5)*22
      data hsunit( 0) / 'Dose [Gy/source]      '/
      data hsunit( 1) / 'Dose [MeV/cm^3/source]'/
      data hsunit( 2) / 'Dose [MeV/source]     '/
      data hsunit( 3) / 'Number [1/source]     '/
      data hsunit( 4) / 'Number [1/nsec/source]'/
      data hsunit( 5) / 'Dose [J/m^3/source]   '/

      character cha*1
      data cha /"'"/


      character dc2*4

      character chp(itmxpt)*11  ! kitamura22/03/31
      character chq(itmxpt)*9   ! kitamura22/03/31

*-----------------------------------------------------------------------

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

*-----------------------------------------------------------------------
      real(8),allocatable :: vl(:,:,:)

      real(8),allocatable :: vl_x(:,:),vl_y(:,:),vl_z(:,:)

!-----------------------------------------------------------------------
      integer           icf
      integer           iat, iad, iaf

      icf(ix,iy,iz) = ix + ( iy - 1 ) * nx 
     &              + ( iz - 1 ) * nx * ny
      icf2(ix,iy,iz) = ix + ( iy - 1 ) * (nx+1) 
     &              + ( iz - 1 ) * (nx+1) * (ny+1)
      iat(iad,iaf) = iad + (iaf-1) * 2
*-----------------------------------------------------------------------
      yen  = char(92)
      allocate( ixyz(max((nx+1)*(ny+1),(ny+1)*(nz+1),(nx+1)*(nz+1))) )

      allocate (tw(nt),rdata(2,nfile))
      nrst = 2*nfile+3
      allocate(anatalrst(np,ne+1,1,nt+1,(nx+1)*(ny+1)*(nz+1),
     &                   1,nrst)) ! S.H. 2022.5.20

*-----------------------------------------------------------------------
      np_mxang = np
      if ( itmxang(m).ne.0 .and. iabs(itmxang(m)).lt.np ) then
        np_mxang = iabs(itmxang(m))
      end if

*-----------------------------------------------------------------------

         if( itout(m) .le. 1 ) then

            nei = 0

         else

            nei = 1

         end if

*-----------------------------------------------------------------------
*     for igshow
*-----------------------------------------------------------------------

         if( igsh .eq. 0 ) then

            npg = np
            neg = ne
            ntg = nt

         else

            npg = 1
            neg = nei
            ntg = 1

         end if

*-----------------------------------------------------------------------
*        particle name
*-----------------------------------------------------------------------

               call wtpname(m,np,chp,chq)

*-----------------------------------------------------------------------
*        set mesh volume ( unit = 2, 3, 4 ; vol = 1.0 )
*-----------------------------------------------------------------------
               allocate( vl_x(ny,nz),vl_y(nx,nz),vl_z(nx,ny) )

               allocate( vl(nx,ny,nz) )
               if( itunt(m).ge.2 .and. itunt(m).le.4 )then
                vl(1:nx,1:ny,1:nz)=1.0d0
                vl_x(:,:)=1.0d0
                vl_y(:,:)=1.0d0
                vl_z(:,:)=1.0d0
               else
                vl_x(:,:)=0.0d0
                vl_y(:,:)=0.0d0
                vl_z(:,:)=0.0d0
                do iz=1,nz
                 do iy=1,ny
                  do ix=1,nx
                   vl(ix,iy,iz)=vls(nl,lt,itmcn(m),itvm(m),itmtr(m,4),
     &                  xm(ix),xm(ix+1),
     &                  ym(iy),ym(iy+1),
     &                  zm(iz),zm(iz+1))
                   vl_x(iy,iz) = vl_x(iy,iz) + vl(ix,iy,iz)
                   vl_y(ix,iz) = vl_y(ix,iz) + vl(ix,iy,iz)
                   vl_z(ix,iy) = vl_z(ix,iy) + vl(ix,iy,iz)
                  enddo
                 enddo
                enddo
               endif

*-----------------------------------------------------------------------
*           itunt(m) = 0 : cfac = c2gy
*-----------------------------------------------------------------------
            if( itunt(m) .eq. 0 )then
              cfac=c2gy
            elseif( itunt(m) .eq. 5)then
              cfac = c2jperm3
            else
              cfac=1.0d0
            endif

*-----------------------------------------------------------------------
*           itunt(m) = 4 : /nsec/source
*-----------------------------------------------------------------------

            if( itunt(m) .eq. 4 ) then

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
*           itunt(m) = 0 : Gy/source
*                    = 1 : MeV/cm^3/source
*                    = 2 : MeV/source
*                    = 3 : 1/source
*                    = 4 : 1/nsec/source
*-----------------------------------------------------------------------
*        relative error and  unit conversion
*        c1 : normalization for source
*-----------------------------------------------------------------------

                  cmax = 0.0
                  cmin = 1.e+33
                  isdz = 0
                  stdm = 0.0
                  itstd(m) = 1

         if( igsh .eq. 0 ) then

             if( isumtally(m)/=2 .or. itout(m)==1 ) then
              if ( rsouin .gt. 0d0 ) then ! S.H. (2016.7.26)
               c1 = 1.0d+0 / rsouin
              else
               c1 = 0d0
              end if
             else
               c1 = 1.0d+0
             end if

           do itmpxyz=1,nx*ny*nz
            iz=(itmpxyz-1)/(ny*nx)+1
            iy=(itmpxyz-1-(iz-1)*ny*nx)/nx+1
            ix=itmpxyz-(iy-1)*nx-(iz-1)*ny*nx
            do 100 it = 1, nt
            do 100 ie = nei, ne
            do 100 ip = 1, np
            do ntf = 1, nfile

             if( isumtally(m)/=2 .or. itout(m)==1 ) then   ! Weighted average

               if( tr(ip,ie,it,icf(ix,iy,iz),iat(1,ntf)).gt.0.d0 ) then

                  call calc_deposit_stdev(m,Xa,sigx,
     &                            tr(ip,ie,it,icf(ix,iy,iz),iat(1,ntf)),
     &                            tr(ip,ie,it,icf(ix,iy,iz),iat(2,ntf)),
     &                            rtfac(m)/vl(ix,iy,iz)/tw(it),ip)   ! S.Abe 2015/12/03, add ip
                  tr(ip,ie,it,icf(ix,iy,iz),iat(1,ntf)) = cfac*Xa
                  tr(ip,ie,it,icf(ix,iy,iz),iat(2,ntf)) = sigx
                  if( sigx .gt. stdm ) stdm = sigx

                  if( tr(ip,ie,it,icf(ix,iy,iz),iat(1,ntf)) .gt. cmax )
     &         cmax = tr(ip,ie,it,icf(ix,iy,iz),iat(1,ntf))

                  if( tr(ip,ie,it,icf(ix,iy,iz),iat(1,ntf)) .lt. cmin )
     &         cmin = tr(ip,ie,it,icf(ix,iy,iz),iat(1,ntf))

               else

                  isdz = 1
                  tr(ip,ie,it,icf(ix,iy,iz),iat(2,ntf)) = 0.d+0

               end if

             else
               if( tr(ip,ie,it,icf(ix,iy,iz),iat(1,ntf)).gt.0.d0 ) then
                 if( tr(ip,ie,it,icf(ix,iy,iz),iat(1,ntf)) .gt. cmax
     &            .and. manatally .ne. 0 )  ! T.Sato 2021/01/30
     &        cmax = tr(ip,ie,it,icf(ix,iy,iz),iat(1,ntf))
                 if( tr(ip,ie,it,icf(ix,iy,iz),iat(1,ntf)) .lt. cmin
     &            .and. manatally .ne. 0 )  ! T.Sato 2021/01/30
     &        cmin = tr(ip,ie,it,icf(ix,iy,iz),iat(1,ntf))
               else
                 isdz = 1
                 tr(ip,ie,it,icf(ix,iy,iz),iat(2,ntf)) = 0.d+0
               end if
             end if

! sumover
              call pdepstxyz_sumover_stdev_ntf(0,m,ntf,
     &             ip,ie,it,ix,iy,iz,
     &             nei,
     &             rtfac(m),cfac,tw(it),vl(ix,iy,iz),
     &             tw_sum,vl_x(iy,iz),vl_y(ix,iz),vl_z(ix,iy))

            end do     ! loop end ntf
  100       continue
      enddo  ! xyz loop T.Sato 2021/04/15

!$OMP PARALLEL
!$OMP& private(ipomp,ix,iy,iz,it,ie,ip,ntf,itmpxyz,rdata,answer,rerr)
!$OMP& private(iax,ix_a,iy_a,iz_a,it_a,ie_a,ido_ana,ioe)
!$      npomp=OMP_GET_NUM_THREADS()
!$      ipomp=OMP_GET_THREAD_NUM()
!$      write(*,'(''OpenMP PARALLEL PROCESS'',
!$   &     i4,''/'',i4)') ipOMP+1,npomp
!$OMP DO SCHEDULE(dynamic) reduction(max:cmax) reduction(min:cmin)
           do itmpxyz=1,nx*ny*nz
            iz=(itmpxyz-1)/(ny*nx)+1
            iy=(itmpxyz-1-(iz-1)*ny*nx)/nx+1
            ix=itmpxyz-(iy-1)*nx-(iz-1)*ny*nx
            do 101 it = 1, nt
            do 101 ie = nei, ne
            do 101 ip = 1, np

            if( manatally .eq. 0 ) then ! T.Sato 2020/10/19

               do ntf = 1, nfile
                  rdata(1,ntf)=tr(ip,ie,it,icf(ix,iy,iz),iat(1,ntf))
                  rdata(2,ntf)=tr(ip,ie,it,icf(ix,iy,iz),iat(2,ntf))
               end do

              call usranatal(nfile,rdata,answer,rerr)

              tr(ip,ie,it,icf(ix,iy,iz),1) = answer
              tr(ip,ie,it,icf(ix,iy,iz),2) = rerr
             if ( itout(m) .le. 1 ) then ! ie=nei=ne=0
              anatalrst(ip, 1,1,it,icf2(ix,iy,iz),1,1) = answer
              anatalrst(ip, 1,1,it,icf2(ix,iy,iz),1,2) = rerr
             else
              anatalrst(ip,ie,1,it,icf2(ix,iy,iz),1,1) = answer
              anatalrst(ip,ie,1,it,icf2(ix,iy,iz),1,2) = rerr
             end if

              if(answer.gt.cmax) cmax = answer
              if(answer.gt.0d0 .and. answer.lt.cmin) cmin = answer

            else if( manatally .eq. 1 ) then

               do ntf = 1, nfile
                  rdata(1,ntf)=tr(ip,ie,it,icf(ix,iy,iz),iat(1,ntf))
                  rdata(2,ntf)=tr(ip,ie,it,icf(ix,iy,iz),iat(2,ntf))
               end do

               call anova(nfile,rdata,resc3(m)
     &              ,fmval,unca,uncerr,uratio)

              if ( itout(m) .le. 1 ) then ! ie=nei=ne=0
               anatalrst(ip,1,1,it,icf2(ix,iy,iz),1,1) ! mean
     &              = fmval
               if ( fmval .gt. 0d0 ) then
                  anatalrst(ip,1,1,it,icf2(ix,iy,iz),1,2) ! total uncertainty
     &                 = dsqrt(unca**2 + uncerr**2/resc3(m)) / fmval
                  anatalrst(ip,1,1,it,icf2(ix,iy,iz),1,3) ! syst. uncertainty
     &                 = unca / fmval
                  anatalrst(ip,1,1,it,icf2(ix,iy,iz),1,4) ! stat. uncertainty
     &                 = uncerr /dsqrt(resc3(m)) / fmval
                  anatalrst(ip,1,1,it,icf2(ix,iy,iz),1,5) ! squared unc. ratio
     &                 = uratio
               else
                  anatalrst(ip,1,1,it,icf2(ix,iy,iz),1,2) ! total uncertainty
     &                 = 0d0
                  anatalrst(ip,1,1,it,icf2(ix,iy,iz),1,3) ! syst. uncertainty
     &                 = 0d0
                  anatalrst(ip,1,1,it,icf2(ix,iy,iz),1,4) ! stat. uncertainty
     &                 = 0d0
                  anatalrst(ip,1,1,it,icf2(ix,iy,iz),1,5) ! squared unc. ratio
     &                 = 0d0
               end if
              else
               anatalrst(ip,ie,1,it,ir,1,1) ! mean
     &              = fmval
               if ( fmval .gt. 0d0 ) then
                  anatalrst(ip,ie,1,it,icf2(ix,iy,iz),1,2) ! total uncertainty
     &                 = dsqrt(unca**2 + uncerr**2/resc3(m)) / fmval
                  anatalrst(ip,ie,1,it,icf2(ix,iy,iz),1,3) ! syst. uncertainty
     &                 = unca / fmval
                  anatalrst(ip,ie,1,it,icf2(ix,iy,iz),1,4) ! stat. uncertainty
     &                 = uncerr /dsqrt(resc3(m)) / fmval
                  anatalrst(ip,ie,1,it,icf2(ix,iy,iz),1,5) ! squared unc. ratio
     &                 = uratio
               else
                  anatalrst(ip,ie,1,it,icf2(ix,iy,iz),1,2) ! total uncertainty
     &                 = 0d0
                  anatalrst(ip,ie,1,it,icf2(ix,iy,iz),1,3) ! syst. uncertainty
     &                 = 0d0
                  anatalrst(ip,ie,1,it,icf2(ix,iy,iz),1,4) ! stat. uncertainty
     &                 = 0d0
                  anatalrst(ip,ie,1,it,icf2(ix,iy,iz),1,5) ! squared unc. ratio
     &                 = 0d0
               end if
              end if

            else if( manatally .eq. 2 ) then

               do ntf = 1, nfile
                do ioe = 1, 2
                 if ( itout(m) .le. 1 ) then ! ie=nei=ne=0
                  anatalrst(ip,1,1,it,icf2(ix,iy,iz),1,iat(ioe,ntf))
     &              = tr(ip,ie,it,icf(ix,iy,iz),iat(ioe,ntf))
                 else
                  anatalrst(ip,ie,1,it,icf2(ix,iy,iz),1,iat(ioe,ntf))
     &              = tr(ip,ie,it,icf(ix,iy,iz),iat(ioe,ntf))
                 end if
                end do
               end do

            end if

!sumover
            do iax = 1,itaxn(m)

              ie_a = ie
              it_a = it
              ix_a = ix
              iy_a = iy
              iz_a = iz
              ido_ana = 0
              if(itaxs(m,iax) == 1 .and. ie == 1) then
                ie_a = ne + 1
                ido_ana = 1
              elseif(itaxs(m,iax) == 3 .and. ix == 1) then
                ix_a = nx + 1
                ido_ana = 1
              elseif(itaxs(m,iax) == 4 .and. iy == 1) then
                iy_a = ny + 1
                ido_ana = 1
              elseif(itaxs(m,iax) == 5 .and. iz == 1) then
                iz_a = nz + 1
                ido_ana = 1
              elseif(itaxs(m,iax) == 11 .and. it == 1) then
                it_a = nt + 1
                ido_ana = 1
              endif

              if(ido_ana == 1) then

                 call get_depstxyz_tr_sum_data(m,iax,nfile,
     &                 ip,ie,it,ix,iy,iz,rdata)

                if( manatally .eq. 0 ) then ! T.Sato 2020/10/22

                   call usranatal(nfile,rdata,answer,rerr)

                   if ( itout(m) .le. 1 ) then ! ie=nei=ne=0
                     anatalrst(ip,1,1,it_a,icf2(ix_a,iy_a,iz_a),1,1)
     &               = answer
                     anatalrst(ip,1,1,it_a,icf2(ix_a,iy_a,iz_a),1,2)
     &               = rerr
                   else
                     anatalrst(ip,ie_a,1,it_a,icf2(ix_a,iy_a,iz_a),1,1)
     &               = answer
                     anatalrst(ip,ie_a,1,it_a,icf2(ix_a,iy_a,iz_a),1,2)
     &               = rerr
                   endif

                else if( manatally .eq. 1 ) then

                  call anova(nfile,rdata,resc3(m)
     &               ,fmval,unca,uncerr,uratio)

                  if ( itout(m) .le. 1 ) then ! ie=nei=ne=0

                    anatalrst(ip,1,1,it_a,icf2(ix_a,iy_a,iz_a),1,1) ! mean
     &                = fmval
                    if ( fmval .gt. 0.0d0 ) then
                      anatalrst(ip,1,1,it_a,icf2(ix_a,iy_a,iz_a),1,2) ! total uncertainty
     &                     = dsqrt(unca**2 + uncerr**2/resc3(m)) / fmval
                    else
                      anatalrst(ip,1,1,it_a,icf2(ix_a,iy_a,iz_a),1,2) ! total uncertainty
     &                   = 0.0d0
                    end if
                  else
                    anatalrst(ip,ie_a,1,it_a,icf2(ix_a,iy_a,iz_a),1,1) ! mean
     &                = fmval
                    if ( fmval .gt. 0.0d0 ) then
                      anatalrst(ip,ie_a,1,it_a,icf2(ix_a,iy_a,iz_a),1,2) ! total uncertainty
     &                     = dsqrt(unca**2 + uncerr**2/resc3(m)) / fmval
                    else
                      anatalrst(ip,ie_a,1,it_a,icf2(ix_a,iy_a,iz_a),1,2) ! total uncertainty
     &                   = 0.0d0
                    end if
                  endif

                else if( manatally .eq. 2 ) then
                  do ntf = 1, nfile
                    do ioe = 1, 2
                      if ( itout(m) .le. 1 ) then ! ie=nei=ne=0
                        anatalrst(ip,1,1,it_a,icf2(ix_a,iy_a,iz_a),
     &                            1,iat(ioe,ntf))
     &                  = rdata(ioe,ntf)
                      else
                        anatalrst(ip,ie_a,1,it_a,icf2(ix_a,iy_a,iz_a),
     &                            1,iat(ioe,ntf))
     &                  = rdata(ioe,ntf)
                      endif
                    end do
                  end do

                end if

              endif
            enddo

  101       continue
      enddo  ! xyz loop T.Sato 2021/04/15
!$OMP END DO
!$OMP END PARALLEL

            if( rtstd(m) .gt. 0.d0 .and.
     &          stdm .gt. 0.d0 .and. isdz .eq. 0 .and.
     &          stdm .lt. rtstd(m) ) itstd(m) = 0

         end if     ! igsh

*-----------------------------------------------------------------------
*     out put unit = 31  : temporary number
*-----------------------------------------------------------------------

      do 900 iax = 1, itfln(m)

!OBINATA(2012.8.6): output *.err
        noe = 1
        if ( any( itaxs(m,iax) .eq. (/ 7, 8, 9, 12, 13 /) )
     &      .and. ittwo(m) .ne. 4 .and. icntl .ne. 8 )
     &    noe = 2

        do 900 ioe = 1, noe

         if( igsh .ne. 0 .and.
     &     ( itaxs(m,iax) .lt. 7 .or.
     &       ittwo(m) .eq. 4 .or. ittwo(m) .eq. 5 ) ) goto 900

         if( itall .eq. 2 .and. nobch .lt. maxbch .and.
     &       igsh .eq. 0 ) then

            write(fnume,'(i3.3)') nobch
!OBINATA(2012.8.6): output *.err
            if ( ioe .eq. 1 ) then
              fname = ctfln(m,iax)(1:itfll(m,iax))//'.'//fnume
            else
              call mk_2derrfn(ctfln(m,iax),fname,itfll(m,iax))
            end if

         else

!OBINATA(2012.8.6): output *.err
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

         call depstech(iot,m,iax,1)

*-----------------------------------------------------------------------
* manatally
*  0  :  user defined analysis
*  1  :  systematic uncertainty
*  2  :  c-value dependence
*-----------------------------------------------------------------------

*-----------------------------------------------------------------------
*        energy, x, y, z, time axis
*-----------------------------------------------------------------------
         if( itaxs(m,iax) .eq. 1
     &        .or. itaxs(m,iax) .eq. 3
     &        .or. itaxs(m,iax) .eq. 4
     &        .or. itaxs(m,iax) .eq. 5
     &        .or. itaxs(m,iax) .eq. 11 ) then

          select case( itaxs(m,iax) )
           case ( 1 )
            iDaxis = 1   ! energy
           case ( 3 )
            iDaxis = 3   ! x axis
           case ( 4 )
            iDaxis = 4   ! y axis
           case ( 5 )
            iDaxis = 5   ! z axis
           case ( 11 )
            iDaxis = 9   ! time
          end select
          if( manatally.eq.2 ) iDaxis = 51 ! S.H. 2021.8.15

C rearrange data from anatalrst to anataldata
          call anatal_rearrange_sum(0,0,3,
     &         iDaxis,ntaxis,anatalrst,
     &         np,     ne,       0,    nt,    nx, ny, nz,   0, nrst,
     &         eb,ew,    1,1,  tb,tw, xm, ym, zm,1,
     &         vl,
     &         nij,ibin,nijaxs,nijaxs2,nijaxs3,nfgmax,
     &         cij,cijaxs,cijaxs2,nanataldata,ierr)

          inum = 0

          do ianataldata = 1, nanataldata

           inum = inum + 1

C exchange ID number from ianataldata to ij2,...,ij8
           itmpdata = int(ianataldata-1)
           ij(8) = mod(itmpdata,nij(8))+1
           itmpdata = int((ianataldata-ij(8))/nij(8))
           ij(7) = mod(itmpdata,nij(7))+1
           itmpdata = int((itmpdata-ij(7)+1)/nij(7))
           ij(6) = mod(itmpdata,nij(6))+1
           itmpdata = int((itmpdata-ij(6)+1)/nij(6))
           ij(5) = mod(itmpdata,nij(5))+1
           itmpdata = int((itmpdata-ij(5)+1)/nij(5))
           ij(4) = mod(itmpdata,nij(4))+1
           itmpdata = int((itmpdata-ij(4)+1)/nij(4))
           ij(3) = mod(itmpdata,nij(3))+1
           itmpdata = int((itmpdata-ij(3)+1)/nij(3))
           ij(2) = mod(itmpdata,nij(2))+1

*-----------------------------------------------------------------------

           write(iot,'(/"#",78("-"))')

           if( inum .eq. 1 ) then
              write(iot,'( "#newpage:")')
           else
              write(iot,'( " newpage:")')
           end if

           if( lcx_txt(m) .gt. 0 ) then
              write(iot,'(/"x: ",200a1)')(cx_txt(m)(i:i),i=1,lcx_txt(m))
           else
            if( itaxl(m) .eq. 0 ) then
               write(iot,'(/"x: ",a32)') cijaxs
            else
               write(iot,'(/"x: ",200a1)') (itaxt(m)(i:i),i=1,itaxl(m))
            end if
           end if

           if( lcy_txt(m) .gt. 0 ) then
              write(iot,'( "y: ",200a1)')(cy_txt(m)(i:i),i=1,lcy_txt(m))
           else
            if( itayl(m) .eq. 0 ) then
               write(iot,'( "y: ",a22)') hsunit(itunt(m))
            else
               write(iot,'( "y: ",200a1)') (itayt(m)(i:i),i=1,itayl(m))
            end if
           end if

           if(  itaxs(m,iax).eq.1 .and.
     &          (itety(m) .eq. 3 .or. itety(m) .eq. 5) ) then
              write(iot,'( "p: xlog ylog afac(0.8) form(0.9)")')

           else if(  itaxs(m,iax).eq.3 .and.
     &             (itxty(m) .eq. 3 .or. itxty(m) .eq. 5) ) then
              write(iot,'( "p: xlog ylog afac(0.8) form(0.9)")')

           else if(  itaxs(m,iax).eq.4 .and.
     &             (ityty(m) .eq. 3 .or. ityty(m) .eq. 5) ) then
              write(iot,'( "p: xlog ylog afac(0.8) form(0.9)")')


           else if( itaxs(m,iax).eq.5 .and.
     &             (itzty(m) .eq. 3 .or. itzty(m) .eq. 5) ) then
              write(iot,'( "p: xlog ylog afac(0.8) form(0.9)")')

           else if( itaxs(m,iax).eq.11 .and.
     &             (ittty(m) .eq. 3 .or. ittty(m) .eq. 5) ) then
              write(iot,'( "p: xlog ylog afac(0.8) form(0.9)")')

           else
              write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')

           end if

           if( itanl(m) .gt. 0 ) then
              write(iot,'( "p: ",200a1)')
     &             ( itang(m)(i:i),i = 1, itanl(m) )
           end if

           if( itsans(m) .gt. 0 ) then
              call write_sangel(iot,m,0)
           end if

           if ( manatally .eq. 0 ) then ! user defined analysis
            if ( iteps(m) .ne. 2 ) then
               write(iot,'( "h: n",12x,"x",12x,
     &                   1000(a10,"),hh0",a3," n "))')
     &                    ( chp(i), chl(i) , i = 1, np_mxang ) ! frtati 2021/10/05 np -> np_mxang
            else
               write(iot,'( "h: n",12x,"x",12x,
     &        1000(a1,i1,a9,"),hh0",a3
     &        ," ny",i0," dy",i0,"=[y",i0,"*y",i0,"] "))')
     &        ( chp(i)(1:1),i,chp(i)(2:10), chl(i)
     &        ,i+np_mxang,i,i,i+np_mxang, i = 1, np_mxang ) ! frtati 2021/10/05
            end if
            write(iot,'( "#  lower        upper  ",3x,
     &                  1000(a1,2x,a8,4x,"r.err "))')
     &                   ( chb(i), chp(i)(3:10), i = 1, np ) ! frtati 2021/10/05

           else if ( manatally .eq. 1 ) then ! systematic uncertainty
            if ( iteps(m) .ne. 2 ) then
               write(iot,'( "h: n",12x,"x",12x,
     &                   1000(a10,"),hh0",a3," n n n "))')
     &                    ( chp(i), chl(i) , i = 1, np_mxang ) ! frtati 2021/10/05 np -> np_mxang
            else
               write(iot,'( "h: n",12x,"x",12x,
     &        1000(a1,i1,a9,"),hh0",a3
     &        ," ny",i0," dy",i0,"=[y",i0,"*y",i0,"] n n "))')
     &        ( chp(i)(1:1),i,chp(i)(2:10), chl(i)
     &        ,i+np_mxang,i,i,i+np_mxang, i = 1, np_mxang ) ! frtati 2021/10/05
            end if
            write(iot,'( "#  lower        upper  ",3x,
     &                  1000(a1,2x,a8,4x,"r.err(tot, syst, stat)"))')
     &                   ( chb(i), chp(i)(3:10), i = 1, np ) ! frtati 2021/10/05

           else if ( manatally .eq. 2 ) then ! c-value dependence
            if ( iteps(m) .ne. 2 ) then
               write(iot,'( "h: n",12x,"x",12x,
     &                   1000(a10,"),hh0",a3," n "))')
     &                    ( chp(i), chl(i) , i = 1, np_mxang ) ! frtati 2021/10/05 np -> np_mxang
            else
               write(iot,'( "h: n",12x,"x",12x,
     &        1000(a1,i1,a9,"),hh0",a3
     &        ," ny",i0," dy",i0,"=[y",i0,"*y",i0,"] "))')
     &        ( chp(i)(1:1),i,chp(i)(2:10), chl(i)
     &        ,i+np_mxang,i,i,i+np_mxang, i = 1, np_mxang ) ! frtati 2021/10/05
            end if
            write(iot,'( "#  lower        upper  ",3x,
     &                  1000(a1,2x,a8,4x,"r.err "))')
     &                   ( chb(i), chp(i)(3:10), i = 1, np ) ! frtati 2021/10/05

           end if

           do ip = 1, np
              tott(ip,1) = 0.d+0
           end do

           voll = 0.0d0
           njaxs = nijaxs + 1

           do ijaxs = 1, nijaxs

            if ( manatally .eq. 0 ) then ! user defined analysis
               write(iot,'(1p2e13.4,1000(1pe13.4,0pf8.4))')
     &              fgaxs(ijaxs),fgaxs(ijaxs+1),
     &              ((anataldata(ip,ianataldata,ijaxs,irst),irst=1,2)
     &              ,ip=1,np) ! frtati 2021/10/05

            else if ( manatally .eq. 1 ) then ! systematic uncertainty
               write(iot,'(1p2e13.4,1000(1pe13.4,0p3f8.4))')
     &              fgaxs(ijaxs),fgaxs(ijaxs+1),
     &              ((anataldata(ip,ianataldata,ijaxs,irst),irst=1,4)
     &              ,ip=1,np) ! frtati 2021/10/05

            else if ( manatally .eq. 2 ) then ! c-value dependence
               write(iot,'(1p2e13.4,1000(1pe13.4,0pf8.4))')
     &              fgaxs(ijaxs),fgaxs(ijaxs),
     &              ((anataldata(ip,ianataldata,ijaxs,irst),irst=1,2)
     &              ,ip=1,np) ! frtati 2021/10/05

            end if

            voll = voll + delvol(ianataldata,ijaxs)

            do ip = 1, np
               vn = anataldata(ip,ianataldata,ijaxs,1)
     &              * delvol(ianataldata,ijaxs)
               tott(ip,1) = tott(ip,1) + vn
            end do

           end do                ! ijaxs = 1, nijaxs

           if( itaxs(m,iax) .eq. 1 ) then ! when axis = eng
            do ip = 1, np
               tott(ip,1) = tott(ip,1)
            end do

           else if( itaxs(m,iax) .eq. 3 .or. ! when axis = x, y, z
     &           itaxs(m,iax) .eq. 4 .or.
     &           itaxs(m,iax) .eq. 5 ) then
            do ip = 1, np
             if( itunt(m) .eq. 1 .or. itunt(m) .eq. 5 ) then
                tott(ip,1) = tott(ip,1) / voll
             end if
            end do

           else if( itaxs(m,iax) .eq. 11 ) then ! when axis = time
            do ip = 1, np
             if( itunt(m) .eq. 4 ) then
                tott(ip,1) = tott(ip,1) / voll
             end if
            end do

           end if

! sumover
            do ip = 1, np
            do irst=1,2
              tott(ip,irst) = anataldata(ip,ianataldata,njaxs,irst)
            enddo
            enddo

           if( itaxs(m,iax).eq.1 .or. itaxs(m,iax).eq.11 ) then ! when axis = eng,t
            if( itunt(m) .ne. 0 ) then
             if ( manatally .eq. 0 ) then ! user defined analysis
                write(iot,'(/"#   sum over",14x,1000(1pe13.4,0pf8.4))')
     &               (tott(ip,1),tott(ip,2),ip=1,np) ! frtati 2021/10/05
             else if ( manatally .eq. 1 ) then ! systematic uncertainty
             write(iot,'(/"#   sum over",14x,1000(1pe13.4,0pf8.4,16x))')
     &               (tott(ip,1),tott(ip,2),ip=1,np) ! frtati 2021/10/05
             else if ( manatally .eq. 2 ) then ! c-value dependence
                write(iot,'(/"#   sum over",14x,1000(1pe13.4,0pf8.4))')
     &               (tott(ip,1),tott(ip,2),ip=1,np) ! frtati 2021/10/05
             end if
            end if

           else if( itaxs(m,iax).eq.3
     &           .or. itaxs(m,iax).eq.4
     &           .or. itaxs(m,iax).eq.5 ) then ! when axis = x,y,z
            if( itunt(m) .eq. 1 .or. itunt(m) .eq. 5 ) then
             if ( manatally .eq. 0 ) then ! user defined analysis
                write(iot,'(/"#   sum over ",1pe13.4,
     &                    1000(1pe13.4,0pf8.4))')
     &                    voll, (tott(ip,1),tott(ip,2),ip=1,np) ! frtati 2021/10/05
             else if ( manatally .eq. 1 ) then ! systematic uncertainty
                write(iot,'(/"#   sum over ",1pe13.4,
     &                    1000(1pe13.4,0pf8.4,16x))')
     &                    voll, (tott(ip,1),tott(ip,2),ip=1,np) ! frtati 2021/10/05
             else if ( manatally .eq. 2 ) then ! c-value dependence
                write(iot,'(/"#   sum over ",1pe13.4,
     &                    1000(1pe13.4,0pf8.4))')
     &                    voll, (tott(ip,1),tott(ip,2),ip=1,np) ! frtati 2021/10/05
             end if

            else if( itunt(m) .ne. 0 ) then
             if ( manatally .eq. 0 ) then ! user defined analysis
                write(iot,'(/"#   sum over ",13x,1000(1pe13.4,0pf8.4))')
     &                          (tott(ip,1),tott(ip,2),ip=1,np) ! frtati 2021/10/05
             else if ( manatally .eq. 1 ) then ! systematic uncertainty
                write(iot,'(/"#   sum over ",13x,
     &                1000(1pe13.4,0pf8.4,16x))')
     &                          (tott(ip,1),tott(ip,2),ip=1,np) ! frtati 2021/10/05
             else if ( manatally .eq. 2 ) then ! c-value dependence
                write(iot,'(/"#   sum over ",13x,1000(1pe13.4,0pf8.4))')
     &                          (tott(ip,1),tott(ip,2),ip=1,np) ! frtati 2021/10/05
             end if
            end if

           end if

*-----------------------------------------------------------------------

           write(iot,'(/a1,"no. =",i3,a1)') cha, inum, cha

           write(iot,'("msuc: {",a1,"huge ",80a1)')
     &          yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

           write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
           write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

C write tally condition of each figure
           write(iot,'("wt: s(0.7)",/a1,"vspace{-3}")') yen
           itmp = 0
           do ijtmp = 2,8
            if (cij(ijtmp) .ne. 'F' ) then
             itmp = itmp + ij(ijtmp)
             if ( ibin(ijtmp) .eq. 0 ) then
                write(iot,'(a15,"&=&",i5)')
     &               cij(ijtmp), idnint(fg(itmp))
                itmp = itmp - ij(ijtmp) + nij(ijtmp)
             else
                write(iot,'(a15,2x,"(",i5,"-th bin)"
     &                            /1pe13.4,2x,"$--$",1pe13.4)')
     &               cij(ijtmp),ij(ijtmp),fg(itmp),fg(itmp+1)
                itmp = itmp - ij(ijtmp) + nij(ijtmp)+1
             end if
            end if
           end do

           write(iot,'("e:")')

          end do                 ! ianataldata = 1, nanataldata

          deallocate( fg )
          deallocate( fgaxs )
          deallocate( anataldata )
          deallocate( delvol )

*-----------------------------------------------------------------------
*        xy, yz, xz, t-eng, eng-t axes (matrix)
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 7
     &       .or. itaxs(m,iax) .eq. 8
     &       .or. itaxs(m,iax) .eq. 9
     &       .or. itaxs(m,iax) .eq. 12
     &       .or. itaxs(m,iax) .eq. 13 ) then

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

          select case( itaxs(m,iax) )
           case ( 7 )
            iDaxis = 31   ! xy
           case ( 8 )
            iDaxis = 32   ! yz
           case ( 9 )
            iDaxis = 33   ! xz
           case ( 12 )
            iDaxis = 35   ! t-eng
           case ( 13 )
            iDaxis = 36   ! eng-t
          end select

C rearrange data from anatalrst to anataldata
          call anatal_rearrange_sum(0,0,3,
     &         iDaxis,ntaxis,anatalrst,
     &         np,     ne,       0,    nt,    nx, ny, nz,   0, nrst,
     &         eb,ew,    1,1,  tb,tw, xm, ym, zm,1,
     &         vl,
     &         nij,ibin,nijaxs,nijaxs2,nijaxs3,nfgmax,
     &         cij,cijaxs,cijaxs2,nanataldata,ierr)

          inum = 0

          do ip = 1, npg
          do ianataldata = 1, nanataldata

           inum = inum + 1

C exchange ID number from ianataldata to ij1,...,ij6
           itmpdata = int(ianataldata-1)
           ij(8) = mod(itmpdata,nij(8))+1
           itmpdata = int((ianataldata-ij(8))/nij(8))
           ij(7) = mod(itmpdata,nij(7))+1
           itmpdata = int((itmpdata-ij(7)+1)/nij(7))
           ij(6) = mod(itmpdata,nij(6))+1
           itmpdata = int((itmpdata-ij(6)+1)/nij(6))
           ij(5) = mod(itmpdata,nij(5))+1
           itmpdata = int((itmpdata-ij(5)+1)/nij(5))
           ij(4) = mod(itmpdata,nij(4))+1
           itmpdata = int((itmpdata-ij(4)+1)/nij(4))
           ij(3) = mod(itmpdata,nij(3))+1
           itmpdata = int((itmpdata-ij(3)+1)/nij(3))
           ij(2) = mod(itmpdata,nij(2))+1

           write(iot,'(/"#",78("-"))')

           if( inum .eq. 1 .or. ip.gt.np_mxang ) then ! frtati 2021/10/05
              write(iot,'( "#newpage:")')
           else
              write(iot,'( " newpage:")')
           end if

           if ( ip.gt.np_mxang ) then
             write(iot,'( " SKIPPAGE:")')
             inum = inum - 1
           end if

*-----------------------------------------------------------------------

           if( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) then

              write(iot,'(/a1,"no. =",i3,a1)') cha, inum, cha

              write(iot,'("msuc: {",a1,"huge ",80a1)')
     &             yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

              write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
              write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

           if( lcx_txt(m) .gt. 0 ) then
              write(iot,'(/"x: ",200a1)')(cx_txt(m)(i:i),i=1,lcx_txt(m))
           else
              if( itaxl(m) .eq. 0 ) then
                 write(iot,'(/"x: ",a32)') cijaxs
              else
                 write(iot,'(/"x: ",200a1)')(itaxt(m)(i:i),i=1,itaxl(m))
              end if
           end if

           if( lcy_txt(m) .gt. 0 ) then
              write(iot,'( "y: ",200a1)')(cy_txt(m)(i:i),i=1,lcy_txt(m))
           else
              if( itayl(m) .eq. 0 ) then
                 write(iot,'("y: ",a32)') cijaxs2
              else
                 write(iot,'( "y: ",200a1)')
     &                (itayt(m)(i:i),i=1,itayl(m))
              end if
           end if

*-----------------------------------------------------------------------

              if( inum .eq. 1 ) then

                 form  = ( fgaxs2(nijaxs2+1) - fgaxs2(1) )
     &                / ( fgaxs(nijaxs+1) - fgaxs(1) )
                 xfac  = 0.9
                 afac  = 0.8
                 izlog = 1
                 inocm = 1
                 inolg = 1

                 if( itanl(m) .gt. 0 )
     &                call anset(itang(m),itanl(m),
     &                form,xfac,afac,cmin,cmax,izlog,inocm,inolg)

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
     &             ( ittwo(m) .ge. 6 .and. ittwo(m) .le. 7 ) ) .and.
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
     &                write(iot,'( "p: zlog")')

              end if

*-----------------------------------------------------------------------

              axs1min = fgaxs(1)
              axs1max = fgaxs(nijaxs+1)
              axs2min = fgaxs2(1)
              axs2max = fgaxs2(nijaxs2+1)
              axs1del = (fgaxs(nijaxs+1)-fgaxs(1))/dble(nijaxs)
              axs2del = (fgaxs2(nijaxs2+1)-fgaxs2(1))/dble(nijaxs2)

              write(iot,'( "p: xmin(",1p1g14.7,") xmax(",1p1g14.7,
     &                       ")")') axs1min, axs1max

              write(iot,'( "p: ymin(",1p1g14.7,") ymax(",1p1g14.7,
     &                       ")")') axs2min, axs2max

              if( itanl(m) .gt. 0 .and. ioe .eq. 1 ) then
                 write(iot,'( "p: ",200a1)')
     &                ( itang(m)(i:i),i = 1, itanl(m) )

              else if( itanl(m) .gt. 0 .and. ioe .ne. 1 ) then

                 call terrang(iot,m)

              end if

              if( itsans(m) .gt. 0 ) then
                 call write_sangel(iot,m,0)
              end if

           end if

*-----------------------------------------------------------------------

           write(iot,'("#  n2 = ",i3,"   n1 = ",i3)')
     &          nijaxs2, nijaxs

           if( ittwo(m) .ne. 4 ) then
              write(iot,'( "# ( ( data(i1,i2), i1 = 1, n1 ),",
     &                         " i2 = n2, 1, -1 )")')
           end if

*-----------------------------------------------------------------------

           if( ( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) .and.
     &          igsh .eq. 0 ) then

              write(iot,'(/a4," y = ",1p1g14.7," to ",
     &            1p1g14.7," by ",1p1g14.7," ; x = ",
     &            1p1g14.7," to ",1p1g14.7," by ",
     &            1p1g14.7," ;")') dc2,
     &            axs2max - axs2del/2d0, axs2min + axs2del/2d0, axs2del,
     &            axs1min + axs1del/2d0, axs1max - axs1del/2d0, axs1del

!OBINATA(2012.8.6): output *.err
              write(iot,'(1p10e11.3)')
     &       ( ( anataldata(ip,ianataldata,ijaxs+(ijaxs2-1)*nijaxs,ioe),
     &                ijaxs = 1, nijaxs ), ijaxs2 = nijaxs2, 1, -1 )

*-----------------------------------------------------------------------

           else if( ittwo(m) .eq. 4 ) then

              write(iot,'(/"# axis2      axis1    ",
     &                      "  Flux       r.err")')

              do ijaxs = 1, nijaxs
               do ijaxs2 = 1, nijaxs2

                  write(iot,'(1p3e13.4,0pf8.4)')
     &                 fgaxs2(ijaxs2)  + axs2del/2d0,
     &                 fgaxs(ijaxs)  + axs1del/2d0,
     &             anataldata(ip,ianataldata,ijaxs+(ijaxs2-1)*nijaxs,1),
     &             anataldata(ip,ianataldata,ijaxs+(ijaxs2-1)*nijaxs,2)

               end do
              end do

*-----------------------------------------------------------------------

           else if( ittwo(m) .eq. 5 ) then

              write(iot,'("#  ax2= ",1p1g14.7," to ",
     &            1p1g14.7," by ",1p1g14.7/
     &                        "#  ax1= ",1p1g14.7," to ",
     &            1p1g14.7," by ",1p1g14.7)')
     &            axs2min + axs2del/2d0, axs2max - axs2del/2d0, axs2del,
     &            axs1min + axs1del/2d0, axs1max - axs1del/2d0, axs1del

              write(iot,'(/4x,a3,4x,1p1000e11.3)')
     &          'ax2/ax1',( fgaxs(ijaxs) + axs1del/2d0, ijaxs=1,nijaxs )

              do ijaxs2 = nijaxs2, 1, -1

!OBINATA(2012.8.6): output *.err
                 write(iot,'(1p1000e11.3)')
     &                fgaxs2(ijaxs2)  + axs2del/2d0,
     &         ( anataldata(ip,ianataldata,ijaxs+(ijaxs2-1)*nijaxs,ioe),
     &                 ijaxs = 1, nijaxs )

              end do

           end if

*-----------------------------------------------------------------------
*        gshow
*-----------------------------------------------------------------------

         if ( any(itaxs(m,iax) .eq. (/ 7, 8, 9 /)) ) then ! xy, yz, xz

         if( itgsh(m) .ne. 0 .and.
     &     ( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) ) then

               write(iot,'(/"#",78("-"))')
               write(iot,'( "# gshow")')
               write(iot,'( "#",78("-"))')
               write(iot,'(/"p: legs[c5*0.875]")')

               if ( itaxs(m,iax) .eq. 7 ) then ! when axis = xy (x-y)
                  yval = ( fgaxs3(ij(5)) + fgaxs3(ij(5)+1) ) / 2.0d0
               else if ( itaxs(m,iax) .eq. 8 ) then ! when axis = yz (z-y)
                  yval = ( fgaxs3(ij(3)) + fgaxs3(ij(3)+1) ) / 2.0d0
               else if ( itaxs(m,iax) .eq. 9 ) then ! when axis = xz (z-x)
                  yval = ( fgaxs3(ij(4)) + fgaxs3(ij(4)+1) ) / 2.0d0
               end if
               none = 1
               iaxs = 3
               iuni = itgsh(m) * 2 - 1 + igsh
               ires = itres(m)
               widt = rtwid(m)

               call gshow(0,iot,iaxs,iuni,widt,ires,
     &                    inum,igsh,igser,gfnam,igfmn,
     &                nijaxs+1,nijaxs2+1,none,fgaxs,fgaxs2,yval,ixyz(1),
     &                    1,1,krr,vll,itmtr(m,4),itglt(m))

         end if

         end if   ! itaxs(m,iax) = 7(xy), 8(yz), 9(xz)

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

      else

       if( lcz_txt(m) .gt. 0 ) then
          write(iot,'("y: ",200a1)')(cz_txt(m)(i:i),i=1,lcz_txt(m))
       else
        if( itazl(m) .eq. 0 ) then
            write(iot,'("y: ",a22)') hsunit(itunt(m))
        else
            write(iot,'("y: ",200a1)') (itazt(m)(i:i),i=1,itazl(m))
        end if
       end if

      end if

      if( izlog .eq. 0 ) write(iot,'( "p: ylin")')

      write(iot,'( "#",78("-"))')

      end if

*-----------------------------------------------------------------------

      if( inocm .eq. 1 ) then

C write tally condition of each figure
       write(iot,'("wt: s[c5]",/a1,"vspace{-3}")') yen
       itmp = 0
       do ijtmp = 2,8
        if (cij(ijtmp) .ne. 'F' ) then
         itmp = itmp + ij(ijtmp)
         if ( ibin(ijtmp) .eq. 0 ) then
            write(iot,'(a15,"&=&",i5)')
     &           cij(ijtmp), idnint(fg(itmp))
            itmp = itmp - ij(ijtmp) + nij(ijtmp)
         else
            write(iot,'(a15,2x,"(",i5,"-th bin)"
     &                             /1pe13.4,2x,"$--$",1pe13.4)')
     &           cij(ijtmp),ij(ijtmp),fg(itmp),fg(itmp+1)
            itmp = itmp - ij(ijtmp) + nij(ijtmp)+1
         end if
        end if
       end do
       write(iot,'("Particle = ",a8)') chq(ip)

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
     &                       ")")') axs1min, axs1max

               write(iot,'( "p: ymin(",1p1g14.7,") ymax(",1p1g14.7,
     &                       ")")') axs2min, axs2max

      end if

*-----------------------------------------------------------------------

          end do   ! ianataldata = 1, nanataldata
          end do   ! ip = 1, npg

*-----------------------------------------------------------------------

         if ( any(itaxs(m,iax) .eq. (/ 7, 8, 9 /)) ) then ! xy, yz, xz

            if ( itbmp(m) .ne. 0 ) then

               numIndex = 4
               allocate( bmpfIType(1:numIndex) )
               allocate( bmpfIndex(1:numIndex) )
               if ( itaxs(m,iax) .eq. 7 ) then ! when axis = xy (x-y)
                  bmpfIType = (/ 't', 'e', 'p', 'z' /)
               else if ( itaxs(m,iax) .eq. 8 ) then ! when axis = yz (z-y)
                  bmpfIType = (/ 't', 'e', 'p', 'x' /)
               else if ( itaxs(m,iax) .eq. 9 ) then ! when axis = xz (z-x)
                  bmpfIType = (/ 't', 'e', 'p', 'y' /)
               end if
                  bmpfIndex = (/ ntg, (neg-nei+1), npg, nijaxs3 /)
                  bmpWidth  = nijaxs
                  bmpHeight = nijaxs2

            end if

            if ( itvtk(m) .ne. 0 ) then

               call open_file(isunit_vtk_meta_default,
     &                 "", isunit_vtk_meta, ios, .true.)
               call open_file(isunit_vtk_default,
     &                 "", isunit_vtk, ios, .true.)

               if ( itaxs(m,iax) .eq. 7 ) then ! when axis = xy (x-y)
                  iaxs = 1
               else if ( itaxs(m,iax) .eq. 8 ) then ! when axis = yz (z-y)
                  iaxs = 2
               else if ( itaxs(m,iax) .eq. 9 ) then ! when axis = xz (z-x)
                  iaxs = 3
               end if
               nparam = 2

               write(isunit_vtk_meta) ntg

               do it = 1, ntg

                  write(isunit_vtk_meta) iaxs
                  write(isunit_vtk_meta) nparam
                  write(isunit_vtk_meta) 'p', 'e'
                  write(isunit_vtk_meta) npg, (neg-nei+1)

                  if ( itaxs(m,iax) .eq. 7 ) then ! when axis = xy (x-y)
                     write(isunit_vtk_meta) nijaxs,nijaxs2,nijaxs3
                     write(isunit_vtk_meta)
     &                    ( fgaxs(ijaxs), ijaxs=1,nijaxs+1 )
                     write(isunit_vtk_meta)
     &                    ( fgaxs2(ijaxs2), ijaxs2=1,nijaxs2+1 )
                     write(isunit_vtk_meta)
     &                    ( fgaxs3(ijaxs3), ijaxs3=1,nijaxs3+1 )
                  else if ( itaxs(m,iax) .eq. 8 ) then ! when axis = yz (z-y)
                     write(isunit_vtk_meta) nijaxs3,nijaxs2,nijaxs
                     write(isunit_vtk_meta)
     &                    ( fgaxs3(ijaxs3), ijaxs3=1,nijaxs3+1 )
                     write(isunit_vtk_meta)
     &                    ( fgaxs2(ijaxs2), ijaxs2=1,nijaxs2+1 )
                     write(isunit_vtk_meta)
     &                    ( fgaxs(ijaxs), ijaxs=1,nijaxs+1 )
                  else if ( itaxs(m,iax) .eq. 9 ) then ! when axis = xz (z-x)
                     write(isunit_vtk_meta) nijaxs2,nijaxs3,nijaxs
                     write(isunit_vtk_meta)
     &                    ( fgaxs2(ijaxs2), ijaxs2=1,nijaxs2+1 )
                     write(isunit_vtk_meta)
     &                    ( fgaxs3(ijaxs3), ijaxs3=1,nijaxs3+1 )
                     write(isunit_vtk_meta)
     &                    ( fgaxs(ijaxs), ijaxs=1,nijaxs+1 )
                  end if

                  do ip = 1, npg
                  do ie = nei, neg

                    write(isunit_vtk_meta) ip, (ie-nei+1)

                   if ( itaxs(m,iax) .eq. 7 ) then ! when axis = xy (x-y)
                     write(isunit_vtk)
     &                    ( ( ( anataldata(ip,
     &                     1+((ijaxs3-1)+(it-1+(ie-1)*ntg)*nijaxs3),
     &                    ijaxs+(ijaxs2-1)*nijaxs,ioe),
     &                             ijaxs = 1, nijaxs ),
     &                             ijaxs2 = 1, nijaxs2 ),
     &                             ijaxs3 = 1, nijaxs3 )
                   else if ( itaxs(m,iax) .eq. 8 ) then ! when axis = yz (z-y)
                     write(isunit_vtk)
     &                    ( ( ( anataldata(ip,
     &                     1+((ijaxs3-1)+(it-1+(ie-1)*ntg)*nijaxs3),
     &                    ijaxs+(ijaxs2-1)*nijaxs,ioe),
     &                             ijaxs3 = 1, nijaxs3 ),
     &                             ijaxs2 = 1, nijaxs2 ),
     &                             ijaxs = 1, nijaxs )
                   else if ( itaxs(m,iax) .eq. 9 ) then ! when axis = xz (z-x)
                     write(isunit_vtk)
     &                    ( ( ( anataldata(ip,
     &                     1+((ijaxs3-1)+(it-1+(ie-1)*ntg)*nijaxs3),
     &                    ijaxs+(ijaxs2-1)*nijaxs,ioe),
     &                             ijaxs2 = 1, nijaxs2 ),
     &                             ijaxs3 = 1, nijaxs3 ),
     &                             ijaxs = 1, nijaxs )
                   end if

                  end do   ! ie = nei, neg
                  end do   ! ip = 1, npg

               end do      ! it = 1, ntg

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

             if ( itaxs(m,iax) .eq. 7 ) then ! when axis = xy (x-y)
               call vtk_set_cell_regmat(
     &                 isunit_vtk_rm,
     &                 isunit_vtk_geom, isunit_vtk_geom_meta,
     &                 nijaxs, nijaxs2, nijaxs3, fgaxs, fgaxs2, fgaxs3,
     &                 iaxs,iuni,ires, igser,
     &                 1,1,krr,vll,itmtr(m,4))
             else if ( itaxs(m,iax) .eq. 8 ) then ! when axis = yz (z-y)
               call vtk_set_cell_regmat(
     &                 isunit_vtk_rm,
     &                 isunit_vtk_geom, isunit_vtk_geom_meta,
     &                 nijaxs3, nijaxs2, nijaxs, fgaxs3, fgaxs2, fgaxs,
     &                 iaxs,iuni,ires, igser,
     &                 1,1,krr,vll,itmtr(m,4))
             else if ( itaxs(m,iax) .eq. 9 ) then ! when axis = xz (z-x)
               call vtk_set_cell_regmat(
     &                 isunit_vtk_rm,
     &                 isunit_vtk_geom, isunit_vtk_geom_meta,
     &                 nijaxs2, nijaxs3, nijaxs, fgaxs2, fgaxs3, fgaxs,
     &                 iaxs,iuni,ires, igser,
     &                 1,1,krr,vll,itmtr(m,4))
             end if   ! itaxs(m,iax)

            end if    ! itvtk(m)

         end if   ! itaxs(m,iax) = 7(xy), 8(yz), 9(xz)

*-----------------------------------------------------------------------

         deallocate( fg )
         deallocate( fgaxs )
         deallocate( fgaxs2 )
         deallocate( fgaxs3 )
         deallocate( anataldata )
         deallocate( delvol )

*-----------------------------------------------------------------------

         end if                    ! itaxis

*-----------------------------------------------------------------------

         close(iot)

         if( iteps(m) .ne. 0 ) then

            open(iot, file = fname, status = 'unknown' )
            call a_angel(idasa,fname)

         end if

*-----------------------------------------------------------------------

  900    continue     ! ioe, iax

         deallocate( ixyz )
         deallocate( vl ) !FURUTA20190208
*-----------------------------------------------------------------------
      deallocate(anatalrst)
      
      deallocate(vl_x,vl_y,vl_z)
      deallocate (tw,rdata)


      return
      end subroutine anatal_pdepstxyz


************************************************************************
*                                                                      *
      subroutine anatal_pdepsttet(m,np,nr,mr,ne,nt,eb,tb,tr,nfile,
     &                    weightRate,
     &                    nx,ny,nz,kr,xm,ym,zm,igsh,idasa,manatally)
*                                                                      *
*       anatally calculation and output of t-deposit with tet mesh     *
*                                                                      *
************************************************************************
      use sumtallymod , only : isumtally
      use sangelmod, only: itsans
      use anatallymod, only: anataldata,fg,fgaxs,fgaxs2,fgaxs3,delvol,
     &     cx_txt,lcx_txt,cy_txt,lcy_txt,cz_txt,lcz_txt
      use partmod ! frtati 2021/10/05

*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

*-----------------------------------------------------------------------

      include 'param.inc'

      common /stat / istdev, irestart, ireschk

*-----------------------------------------------------------------------

      double precision resc2   ! W
      double precision resc3   ! N
      common /restart/ resc2(itlmax), resc3(itlmax)

*-----------------------------------------------------------------------

      common /taliin/ rsouin, nzztin, nrgnin

*-----------------------------------------------------------------------

      common /verjam/ versn, lastr, iyeav, imonv, idayv

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

      common /tall20/ ittty(itlmax), ittnm(itlmax), ittrg(itlmax),
     &                rttmi(itlmax), rttma(itlmax), rttdl(itlmax)
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

      common /istcut/ ist_cut, ist_bat
*-----------------------------------------------------------------------

      common /volreg/ dvol(kvlmax)

*-----------------------------------------------------------------------

      common /tcntl/  icntl, inucr

*-----------------------------------------------------------------------

      common /cparm/  maxbch,maxcas
      common /jcomon/ nabov,nobch,nocas,nomax
!$OMP THREADPRIVATE(/jcomon/)
!$    integer ipomp,npomp
!$    integer OMP_GET_NUM_THREADS,OMP_GET_THREAD_NUM

      common /talout/ itall

      character fname*100, fnume*3

*-----------------------------------------------------------------------
      character gfnam*100
      common /tall60/ itger(itlmax)

      common /tall67/ itstd(itlmax), rtstd(itlmax)

*-----------------------------------------------------------------------

      dimension   kr(mr)
      dimension   eb(ne+1)
      dimension   tb(nt+1)
!      dimension   tw(nt)
      dimension   tr(np,0:ne,nr,nt,2*nfile)
!      dimension   weightRate(nfile) ,rdata(2,nfile)
      dimension   weightRate(nfile)

      dimension   xm(nx+1)
      dimension   ym(ny+1)
      dimension   zm(nz+1)

      integer,allocatable :: ixyz(:)

      integer irst,nrst, itmpdata
!      dimension   anatalrst(np,ne+1,1,nt+1,nr*1*1,1,2*nfile+3)
! sumover
!      dimension   anatalrst(np,ne+1,1,nt+1,(nr+1)*1*1,1,2*nfile+3)
      real(8),allocatable :: anatalrst(:,:,:,:,:,:,:)
      integer ianataldata,nanataldata
      integer, parameter :: ntaxis=9  ! S.H. 2021.8.9
      integer ij(ntaxis),nij(ntaxis),ibin(ntaxis)
      integer nijaxs,nijaxs2,nijaxs3,nfgmax
      character cij(ntaxis)*29,cijaxs*29,cijaxs2*29
      integer iDaxis
!      dimension   dlr(nr+1) !  integer lr(nr) ---> double dlr(nr)
      real(8),allocatable :: tw(:),rdata(:,:),dlr(:)
      integer,allocatable :: lr(:)        !FURUTA20190204
      real(8),allocatable :: vl(:,:,:),val(:) !FURUTA20190204
c Dont know why but necessary to avoid segmentation fault

*-----------------------------------------------------------------------
      real(8) cfac
      real(8),parameter:: c2gy=1.602d-10,c2jperm3=1.602d-7

      dimension tott(np,2) ! frtati 2021/10/05 6 -> np

      character hsunit(0:5)*22
      data hsunit( 0) / 'Dose [Gy/source]      '/
      data hsunit( 1) / 'Dose [MeV/cm^3/source]'/
      data hsunit( 2) / 'Dose [MeV/source]     '/
      data hsunit( 3) / 'Number [1/source]     '/
      data hsunit( 4) / 'Number [1/nsec/source]'/
      data hsunit( 5) / 'Dose [J/m^3/source]   '/

      character cha*1
      data cha /"'"/

      character dc2*4
      character chp(itmxpt)*11  ! kitamura22/03/31
      character chq(itmxpt)*9   ! kitamura22/03/31

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

      dimension dt_one(1)
      data dt_one/1.0d0/

      integer           iat, iad, iaf
      iat(iad,iaf) = iad + (iaf-1) * 2
*-----------------------------------------------------------------------
      yen  = char(92)
      allocate( ixyz(max((nx+1)*(ny+1),(ny+1)*(nz+1),(nx+1)*(nz+1))) )

      allocate (tw(nt),rdata(2,nfile),dlr(nr+1))

      nrst = 2*nfile+3
      allocate (anatalrst(np,ne+1,1,nt+1,(nr+1)*1*1,1,nrst))

*-----------------------------------------------------------------------
      np_mxang = np
      if ( itmxang(m).ne.0 .and. iabs(itmxang(m)).lt.np ) then
        np_mxang = iabs(itmxang(m))
      end if

*-----------------------------------------------------------------------

         if( itout(m) .le. 1 ) then

            nei = 0

         else

            nei = 1

         end if

*-----------------------------------------------------------------------
*     for igshow
*-----------------------------------------------------------------------

         if( igsh .eq. 0 ) then

            npg = np
            neg = ne
            ntg = nt

         else

            npg = 1
            neg = nei
            ntg = 1

         end if

*-----------------------------------------------------------------------
*        particle name
*-----------------------------------------------------------------------

               call wtpname(m,np,chp,chq)

*-----------------------------------------------------------------------
*           itunt(m) = 0 : Gy/source
*                    = 1 : MeV/cm^3/source
*                    = 2 : MeV/source
*                    = 3 : 1/source
*                    = 4 : 1/nsec/source
*                    = 5 : J/m^3/source
*-----------------------------------------------------------------------
*        set mesh volume and name
*-----------------------------------------------------------------------

               allocate ( lr(nr),vl(nr,1,1),val(nr) )
               call ttetvl(mr,kr,nr,vl,lr)

               vl_sum = sum(vl(:,1,1))

ccse 2021.08 add (use anatal_rearrange sub.)
         do ir = 1, nr
            dlr(ir) = dble(lr(ir))
         end do
         dlr(nr+1) = dble(lr(nr))

*-----------------------------------------------------------------------
*        set volume ( unit = 2, 3, 4 ; vol = 1.0 )
*-----------------------------------------------------------------------

            if( itunt(m).ge.2 .and. itunt(m).le.4 )then

               do ir = 1, nr

                  vl(ir,1,1) = 1.0d0

               end do
               vl_sum = 1.0d0

            end if

*-----------------------------------------------------------------------
*           itunt(m) = 0 ; cfac = c2gy
*           itunt(m) = 5 ; cfac = c2jperm3
*-----------------------------------------------------------------------
            if( itunt(m) .eq. 0)then
              cfac = c2gy
            elseif( itunt(m) .eq. 5)then
              cfac = c2jperm3
            else
              cfac = 1.0d0
            endif

*-----------------------------------------------------------------------
*           itunt(m) = 4 : /nsec/source
*-----------------------------------------------------------------------

            if( itunt(m) .eq. 4 ) then

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
*        relative error and  unit conversion
*        c1 : normalization for source
*-----------------------------------------------------------------------

         if( igsh .eq. 0 ) then

                  cmax = 0.0
                  cmin = 1.e+33
                  dnon = 1.e-33
                  isdz = 0
                  stdm = 0.0
                  itstd(m) = 1

                if( isumtally(m)/=2 .or. itout(m)==1 ) then
                 if ( rsouin .gt. 0d0 ) then ! S.H. (2016.7.26)
                  c1 = 1.0d+0 / rsouin
                 else
                  c1 = 0d0
                 end if
                else          ! deposit, weighted average
                  c1 = 1.0d+0
                end if

        do ir = 1, nr
            do 100 it = 1, nt
            do 100 ie = nei, ne
            do 100 ip = 1, np
            do ntf = 1, nfile

             if( isumtally(m)/=2 .or. itout(m)==1 ) then

               if( tr(ip,ie,ir,it,iat(1,ntf)) .gt. 0.d0 ) then

                  call calc_deposit_stdev(m,Xa,sigx,
     &                            tr(ip,ie,ir,it,iat(1,ntf)),
     &                            tr(ip,ie,ir,it,iat(2,ntf)),
     &                            rtfac(m)/vl(ir,1,1)/tw(it),ip)

                  tr(ip,ie,ir,it,iat(1,ntf)) = cfac*Xa
                  tr(ip,ie,ir,it,iat(2,ntf)) = sigx
                  if( sigx .gt. stdm ) stdm = sigx

                  if( tr(ip,ie,ir,it,iat(1,ntf)) .gt. cmax )
     &         cmax = tr(ip,ie,ir,it,iat(1,ntf))

                  if( tr(ip,ie,ir,it,iat(1,ntf)) .lt. cmin )
     &         cmin = tr(ip,ie,ir,it,iat(1,ntf))

               else

                  isdz = 1
                  tr(ip,ie,ir,it,iat(2,ntf)) = 0.0

               end if
             else
               if( tr(ip,ie,ir,it,iat(1,ntf)) .gt. 0.d0 ) then
                 if( tr(ip,ie,ir,it,iat(1,ntf)) .gt. cmax
     &            .and. manatally .ne. 0 )  ! T.Sato 2021/01/30
     &        cmax = tr(ip,ie,ir,it,iat(1,ntf))
                 if( tr(ip,ie,ir,it,iat(1,ntf)) .lt. cmin
     &            .and. manatally .ne. 0 )  ! T.Sato 2021/01/30
     &        cmin = tr(ip,ie,ir,it,iat(1,ntf))
               else
                  isdz = 1
                  tr(ip,ie,ir,it,iat(2,ntf)) = 0.0
               end if
             end if

! sumover
              call pdepstreg_sumover_stdev_ntf(0,m,ntf,
     &             ip,ie,ir,it,
     &             nei,
     &             rtfac(m),cfac,tw(it),vl(ir,1,1),tw_sum,vl_sum)

            end do   ! loop end ntf
  100       continue
      enddo  ! ir loop T.Sato 2021/04/15

!$OMP PARALLEL
!$OMP& private(ipomp,ir,it,ie,ip,ntf,rdata,answer,rerr)
!$OMP& private(iax,ir_a,it_a,ie_a,ido_ana,ioe)
!$      npomp=OMP_GET_NUM_THREADS()
!$      ipomp=OMP_GET_THREAD_NUM()
!$      write(*,'(''OpenMP PARALLEL PROCESS'',
!$   &     i4,''/'',i4)') ipOMP+1,npomp
!$OMP DO SCHEDULE(dynamic) reduction(max:cmax) reduction(min:cmin)
        do ir = 1, nr
            do 101 it = 1, nt
            do 101 ie = nei, ne
            do 101 ip = 1, np

            if( manatally .eq. 0 ) then ! T.Sato 2020/10/19

               do ntf = 1, nfile
                  rdata(1,ntf) = tr(ip,ie,ir,it,iat(1,ntf))
                  rdata(2,ntf) = tr(ip,ie,ir,it,iat(2,ntf))
               end do

              call usranatal(nfile,rdata,answer,rerr)

              tr(ip,ie,ir,it,1) = answer
              tr(ip,ie,ir,it,2) = rerr
             if ( itout(m) .le. 1 ) then ! ie=nei=ne=0
              anatalrst(ip, 1,1,it,ir,1,1) = answer
              anatalrst(ip, 1,1,it,ir,1,2) = rerr
             else
              anatalrst(ip,ie,1,it,ir,1,1) = answer
              anatalrst(ip,ie,1,it,ir,1,2) = rerr
             end if

              if(answer.gt.cmax) cmax = answer
              if(answer.gt.0d0 .and. answer.lt.cmin) cmin = answer

            else if( manatally .eq. 1 ) then

               do ntf = 1, nfile
                  rdata(1,ntf) = tr(ip,ie,ir,it,iat(1,ntf))
                  rdata(2,ntf) = tr(ip,ie,ir,it,iat(2,ntf))
               end do

               call anova(nfile,rdata,resc3(m)
     &              ,fmval,unca,uncerr,uratio)

              if ( itout(m) .le. 1 ) then ! ie=nei=ne=0
               anatalrst(ip,1,1,it,ir,1,1) ! mean
     &              = fmval
               if ( fmval .gt. 0d0 ) then
                  anatalrst(ip,1,1,it,ir,1,2) ! total uncertainty
     &                 = dsqrt(unca**2 + uncerr**2/resc3(m)) / fmval
                  anatalrst(ip,1,1,it,ir,1,3) ! syst. uncertainty
     &                 = unca / fmval
                  anatalrst(ip,1,1,it,ir,1,4) ! stat. uncertainty
     &                 = uncerr /dsqrt(resc3(m)) / fmval
                  anatalrst(ip,1,1,it,ir,1,5) ! squared unc. ratio
     &                 = uratio
               else
                  anatalrst(ip,1,1,it,ir,1,2) ! total uncertainty
     &                 = 0d0
                  anatalrst(ip,1,1,it,ir,1,3) ! syst. uncertainty
     &                 = 0d0
                  anatalrst(ip,1,1,it,ir,1,4) ! stat. uncertainty
     &                 = 0d0
                  anatalrst(ip,1,1,it,ir,1,5) ! squared unc. ratio
     &                 = 0d0
               end if
              else
               anatalrst(ip,ie,1,it,ir,1,1) ! mean
     &              = fmval
               if ( fmval .gt. 0d0 ) then
                  anatalrst(ip,ie,1,it,ir,1,2) ! total uncertainty
     &                 = dsqrt(unca**2 + uncerr**2/resc3(m)) / fmval
                  anatalrst(ip,ie,1,it,ir,1,3) ! syst. uncertainty
     &                 = unca / fmval
                  anatalrst(ip,ie,1,it,ir,1,4) ! stat. uncertainty
     &                 = uncerr /dsqrt(resc3(m)) / fmval
                  anatalrst(ip,ie,1,it,ir,1,5) ! squared unc. ratio
     &                 = uratio
               else
                  anatalrst(ip,ie,1,it,ir,1,2) ! total uncertainty
     &                 = 0d0
                  anatalrst(ip,ie,1,it,ir,1,3) ! syst. uncertainty
     &                 = 0d0
                  anatalrst(ip,ie,1,it,ir,1,4) ! stat. uncertainty
     &                 = 0d0
                  anatalrst(ip,ie,1,it,ir,1,5) ! squared unc. ratio
     &                 = 0d0
               end if
              end if

            else if( manatally .eq. 2 ) then

               do ntf = 1, nfile
                do ioe = 1, 2
                 if ( itout(m) .le. 1 ) then ! ie=nei=ne=0
                  anatalrst(ip,1,1,it,ir,1,iat(ioe,ntf))
     &              = tr(ip,ie,ir,it,iat(ioe,ntf))
                 else
                  anatalrst(ip,ie,1,it,ir,1,iat(ioe,ntf))
     &              = tr(ip,ie,ir,it,iat(ioe,ntf))
                 end if
                end do
               end do

            end if

!sumover
            do iax = 1,itaxn(m)

              ie_a = ie
              it_a = it
              ir_a = ir
              ido_ana = 0
              if(itaxs(m,iax) == 1 .and. ie == 1) then
                ie_a = ne + 1
                ido_ana = 1
              elseif(itaxs(m,iax) == 2 .and. ir == 1) then
                ir_a = nr + 1
                ido_ana = 1
              elseif(itaxs(m,iax) == 11 .and. it == 1) then
                it_a = nt + 1
                ido_ana = 1
              endif

              if(ido_ana == 1) then

                 call get_depstreg_tr_sum_data(m,iax,nfile,
     &                 ip,ie,ir,it,rdata)

                if( manatally .eq. 0 ) then ! T.Sato 2020/10/22

                   call usranatal(nfile,rdata,answer,rerr)

                   if ( itout(m) .le. 1 ) then ! ie=nei=ne=0
                     anatalrst(ip,1,1,it_a,ir_a,1,1) = answer
                     anatalrst(ip,1,1,it_a,ir_a,1,2) = rerr
                   else
                     anatalrst(ip,ie_a,1,it_a,ir_a,1,1) = answer
                     anatalrst(ip,ie_a,1,it_a,ir_a,1,2) = rerr
                   endif

                else if( manatally .eq. 1 ) then

                  call anova(nfile,rdata,resc3(m)
     &               ,fmval,unca,uncerr,uratio)

                  if ( itout(m) .le. 1 ) then ! ie=nei=ne=0

                    anatalrst(ip,1,1,it_a,ir_a,1,1) ! mean
     &                = fmval
                    if ( fmval .gt. 0.0d0 ) then
                      anatalrst(ip,1,1,it_a,ir_a,1,2) ! total uncertainty
     &                     = dsqrt(unca**2 + uncerr**2/resc3(m)) / fmval
                    else
                      anatalrst(ip,1,1,it_a,ir_a,1,2) ! total uncertainty
     &                   = 0.0d0
                    end if
                  else
                    anatalrst(ip,ie_a,1,it_a,ir_a,1,1) ! mean
     &                = fmval
                    if ( fmval .gt. 0.0d0 ) then
                      anatalrst(ip,ie_a,1,it_a,ir_a,1,2) ! total uncertainty
     &                     = dsqrt(unca**2 + uncerr**2/resc3(m)) / fmval
                    else
                      anatalrst(ip,ie_a,1,it_a,ir_a,1,2) ! total uncertainty
     &                   = 0.0d0
                    end if
                  endif

                else if( manatally .eq. 2 ) then
                  do ntf = 1, nfile
                    do ioe = 1, 2
                      if ( itout(m) .le. 1 ) then ! ie=nei=ne=0
                        anatalrst(ip,1,1,it_a,ir_a,1,iat(ioe,ntf))
     &                  = rdata(ioe,ntf)
                      else
                        anatalrst(ip,ie_a,1,it_a,ir_a,1,iat(ioe,ntf))
     &                  = rdata(ioe,ntf)
                      endif
                    end do
                  end do

                end if

              endif
            enddo


  101       continue
      enddo  ! ir loop T.Sato 2021/04/15
!$OMP END DO
!$OMP END PARALLEL

            if( rtstd(m) .gt. 0.d0 .and.
     &          stdm .gt. 0.d0 .and. isdz .eq. 0 .and.
     &          stdm .lt. rtstd(m) ) itstd(m) = 0

         end if     ! igsh

*-----------------------------------------------------------------------
*     out put unit = 31  : temporary number
*-----------------------------------------------------------------------

      do 900 iax = 1, itfln(m)

!OBINATA(2012.8.3): output *.err
        noe = 1
        if ( any( itaxs(m,iax) .eq. (/ 12, 13 /) )
     &      .and. ittwo(m) .ne. 4 .and. icntl .ne. 8 )
     &    noe = 2

        do 900 ioe = 1, noe

         if( igsh .ne. 0 .and.
     &     ( itaxs(m,iax) .lt. 7 .or. itrsh(m) .eq. 0 ) ) goto 900

         if( itall .eq. 2 .and. nobch .lt. maxbch .and.
     &       igsh .eq. 0 ) then

            write(fnume,'(i3.3)') nobch
!OBINATA(2012.8.3): output *.err
            if ( ioe .eq. 1 ) then
              fname = ctfln(m,iax)(1:itfll(m,iax))//'.'//fnume
            else
              call mk_2derrfn(ctfln(m,iax),fname,itfll(m,iax))
            end if

         else

!OBINATA(2012.8.3): output *.err
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

         call depstech(iot,m,iax,1)

*-----------------------------------------------------------------------
* manatally
*  0  :  user defined analysis
*  1  :  systematic uncertainty
*  2  :  c-value dependence
*-----------------------------------------------------------------------

*-----------------------------------------------------------------------
*        energy, tet, time axis
*-----------------------------------------------------------------------

         if( itaxs(m,iax) .eq. 1
     &        .or. itaxs(m,iax) .eq. 14
     &        .or. itaxs(m,iax) .eq. 11 ) then

          select case( itaxs(m,iax) )
           case ( 1 )
            iDaxis = 1   ! energy
           case ( 14 )
            iDaxis = 8   ! tet axis
           case ( 11 )
            iDaxis = 9   ! time
          end select
          if( manatally.eq.2 ) iDaxis = 51 ! S.H. 2021.8.15

C rearrange data from anatalrst to anataldata
          call anatal_rearrange_sum(0,0,4,
     &         iDaxis,ntaxis,anatalrst,
     &         np,     ne,       0,    nt,     nr,  0,  0,   0, nrst,
     &         eb,ew,    1,1,  tb,tw, dlr,  dt_one,  dt_one,1,
     &         vl,
     &         nij,ibin,nijaxs,nijaxs2,nijaxs3,nfgmax,
     &         cij,cijaxs,cijaxs2,nanataldata,ierr)

          inum = 0

          do ianataldata = 1, nanataldata

           inum = inum + 1

C exchange ID number from ianataldata to ij2,...,ij8
           itmpdata = int(ianataldata-1)
           ij(8) = mod(itmpdata,nij(8))+1
           itmpdata = int((ianataldata-ij(8))/nij(8))
           ij(7) = mod(itmpdata,nij(7))+1
           itmpdata = int((itmpdata-ij(7)+1)/nij(7))
           ij(6) = mod(itmpdata,nij(6))+1
           itmpdata = int((itmpdata-ij(6)+1)/nij(6))
           ij(5) = mod(itmpdata,nij(5))+1
           itmpdata = int((itmpdata-ij(5)+1)/nij(5))
           ij(4) = mod(itmpdata,nij(4))+1
           itmpdata = int((itmpdata-ij(4)+1)/nij(4))
           ij(3) = mod(itmpdata,nij(3))+1
           itmpdata = int((itmpdata-ij(3)+1)/nij(3))
           ij(2) = mod(itmpdata,nij(2))+1

*-----------------------------------------------------------------------

           write(iot,'(/"#",78("-"))')

           if( inum .eq. 1 ) then
              write(iot,'( "#newpage:")')
           else
              write(iot,'( " newpage:")')
           end if

           if( lcx_txt(m) .gt. 0 ) then
              write(iot,'(/"x: ",200a1)')(cx_txt(m)(i:i),i=1,lcx_txt(m))
           else
            if( itaxl(m) .eq. 0 ) then
               write(iot,'(/"x: ",a29)') cijaxs
            else
               write(iot,'(/"x: ",200a1)') (itaxt(m)(i:i),i=1,itaxl(m))
            end if
           end if

           if( lcy_txt(m) .gt. 0 ) then
              write(iot,'( "y: ",200a1)')(cy_txt(m)(i:i),i=1,lcy_txt(m))
           else
            if( itayl(m) .eq. 0 ) then
               write(iot,'( "y: ",a22)') hsunit(itunt(m))
            else
               write(iot,'( "y: ",200a1)') (itayt(m)(i:i),i=1,itayl(m))
            end if
           end if

           if( (itaxs(m,iax) .eq. 1) .and.
     &          (itety(m) .eq. 3 .or. itety(m) .eq. 5) ) then
              write(iot,'( "p: xlog ylog afac(0.8) form(0.9)")')

           else if( itaxs(m,iax) .eq. 14 ) then
              write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')

           else if( (itaxs(m,iax) .eq. 11) .and.
     &             (ittty(m) .eq. 3 .or. ittty(m) .eq. 5) ) then
              write(iot,'( "p: xlog ylog afac(0.8) form(0.9)")')

           else
              write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')

           end if

           if( itanl(m) .gt. 0 ) then
              write(iot,'( "p: ",200a1)')
     &             ( itang(m)(i:i),i = 1, itanl(m) )
           end if

           if( itsans(m) .gt. 0 ) then
              call write_sangel(iot,m,0)
           end if

           if ( manatally .eq. 0 ) then ! user defined analysis
            if ( iteps(m) .ne. 2 ) then
             if( itaxs(m,iax) .eq. 14 ) then ! tet axis
                write(iot,'( "h:   x      n",5x,"n",10x,
     &                          1000(a10,"),",a4," n   "))')
     &                           ( chp(i), chm(i) , i = 1, np_mxang ) ! frtati 2021/10/05 np -> np_mxang
             else
                write(iot,'( "h: n",12x,"x",12x,
     &                          1000(a10,"),hh0",a3," n "))')
     &                           ( chp(i), chl(i) , i = 1, np_mxang ) ! frtati 2021/10/05 np -> np_mxang
             end if

            else
             if( itaxs(m,iax) .eq. 14 ) then ! tet axis
                write(iot,'( "h:   x      n",5x,"n",10x,
     &               1000(a1,i1,a9,"),",a4
     &               ," ny",i0," dy",i0,"=[y",i0,"*y",i0,"] "))')
     &               ( chp(i)(1:1),i,chp(i)(2:10), chm(i)
     &               ,i+np_mxang,i,i,i+np_mxang, i = 1, np_mxang ) ! frtati 2021/10/05
             else
                write(iot,'( "h: n",12x,"x",12x,
     &               1000(a1,i1,a9,"),hh0",a3
     &               ," ny",i0," dy",i0,"=[y",i0,"*y",i0,"] "))')
     &               ( chp(i)(1:1),i,chp(i)(2:10), chl(i)
     &               ,i+np_mxang,i,i,i+np_mxang, i = 1, np_mxang ) ! frtati 2021/10/05
             end if
            end if

            if( itaxs(m,iax) .eq. 14 ) then ! tet axis
               write(iot,'( "#  num    tetra   volume  ",
     &                        1000(a1,2x,a8,4x,"r.err "))')
     &                         ( chb(i), chp(i)(3:10), i = 1, np ) ! frtati 2021/10/05
            else
               write(iot,'( "#  lower        upper  ",3x,
     &                       1000(a1,"  ",a8,"    r.err "))')
     &                        ( chb(i), chp(i)(3:10), i = 1, np ) ! frtati 2021/10/05
            end if

           else if ( manatally .eq. 1 ) then ! systematic uncertainty
            if ( iteps(m) .ne. 2 ) then
             if( itaxs(m,iax) .eq. 14 ) then ! tet axis
                write(iot,'( "h:   x      n",5x,"n",10x,
     &                          1000(a10,"),",a4," n n n "))')
     &                           ( chp(i), chm(i) , i = 1, np_mxang ) ! frtati 2021/10/05 np -> np_mxang
             else
                write(iot,'( "h: n",12x,"x",12x,
     &                          1000(a10,"),hh0",a3," n n n "))')
     &                           ( chp(i), chl(i) , i = 1, np_mxang ) ! frtati 2021/10/05 np -> np_mxang
             end if

            else
             if( itaxs(m,iax) .eq. 14 ) then ! tet axis
                write(iot,'( "h:   x      n",5x,"n",10x,
     &               1000(a1,i1,a9,"),",a4
     &               ," ny",i0," dy",i0,"=[y",i0,"*y",i0,"] n n "))')
     &               ( chp(i)(1:1),i,chp(i)(2:10), chm(i)
     &               ,i+np_mxang,i,i,i+np_mxang, i = 1, np_mxang ) ! frtati 2021/10/05
             else
                write(iot,'( "h: n",12x,"x",12x,
     &               1000(a1,i1,a9,"),hh0",a3
     &               ," ny",i0," dy",i0,"=[y",i0,"*y",i0,"] n n "))')
     &               ( chp(i)(1:1),i,chp(i)(2:10), chl(i)
     &               ,i+np_mxang,i,i,i+np_mxang, i = 1, np_mxang ) ! frtati 2021/10/05
             end if
            end if

            if( itaxs(m,iax) .eq. 14 ) then ! tet axis
               write(iot,'( "#  num    tetra   volume  ",
     &                    1000(a1,2x,a8,4x,"r.err(tot, syst, stat)"))')
     &                        ( chb(i), chp(i)(3:10), i = 1, np ) ! frtati 2021/10/05
            else
               write(iot,'( "#  lower        upper  ",3x,
     &                    1000(a1,2x,a8,4x,"r.err(tot, syst, stat)"))')
     &                        ( chb(i), chp(i)(3:10), i = 1, np ) ! frtati 2021/10/05
            end if

         else if ( manatally .eq. 2 ) then ! c-value dependence
          if ( iteps(m) .ne. 2 ) then
             write(iot,'( "h:   x",10x,
     &                  1000(a10,"),",a4," n   "))')
     &                  ( chp(i), chm(i) , i = 1, np_mxang ) ! frtati 2021/10/05 np -> np_mxang

          else
             write(iot,'( "h:   x",10x,
     &               1000(a1,i1,a9,"),",a4
     &               ," ny",i0," dy",i0,"=[y",i0,"*y",i0,"] "))')
     &               ( chp(i)(1:1),i,chp(i)(2:10), chm(i)
     &               ,i+np_mxang,i,i,i+np_mxang, i = 1, np_mxang ) ! frtati 2021/10/05
          end if

          write(iot,'( "#  c-value   ",
     &                        1000(a1,2x,a8,4x,"r.err "))')
     &                         ( chb(i), chp(i)(3:10), i = 1, np ) ! frtati 2021/10/05

         end if

         do ip = 1, np
            tott(ip,1) = 0.d+0
         end do

         voll = 0.0d0
         njaxs = nijaxs + 1

         do ijaxs = 1, nijaxs

          ireg = lr(ijaxs)

          if ( manatally .eq. 0 ) then ! user defined analysis
           if( itaxs(m,iax) .eq. 14 ) then
              write(iot,'(i5,1x,i7,1pe13.4,1000(1pe13.4,0p3f8.4))')
     &             ijaxs, ireg, fgaxs(ijaxs),
     &             ((anataldata(ip,ianataldata,ijaxs,irst),irst=1,2)
     &             ,ip=1,np) ! frtati 2021/10/05
           else
              write(iot,'(1p2e13.4,1000(1pe13.4,0pf8.4))')
     &             fgaxs(ijaxs),fgaxs(ijaxs+1),
     &             ((anataldata(ip,ianataldata,ijaxs,irst),irst=1,2)
     &             ,ip=1,np) ! frtati 2021/10/05
           end if

          else if ( manatally .eq. 1 ) then ! systematic uncertainty
           if( itaxs(m,iax) .eq. 14 ) then
              write(iot,'(i5,1x,i7,1pe13.4,1000(1pe13.4,0p3f8.4))')
     &             ijaxs, ireg, fgaxs(ijaxs),
     &             ((anataldata(ip,ianataldata,ijaxs,irst),irst=1,4)
     &             ,ip=1,np) ! frtati 2021/10/05
           else
              write(iot,'(1p2e13.4,1000(1pe13.4,0p3f8.4))')
     &             fgaxs(ijaxs),fgaxs(ijaxs+1),
     &             ((anataldata(ip,ianataldata,ijaxs,irst),irst=1,4)
     &             ,ip=1,np) ! frtati 2021/10/05
           end if

          else if ( manatally .eq. 2 ) then ! c-value dependence
             write(iot,'(1pe13.4,1000(1pe13.4,0pf8.4))')
     &            fgaxs(ijaxs),
     &            ((anataldata(ip,ianataldata,ijaxs,irst),irst=1,2)
     &            ,ip=1,np) ! frtati 2021/10/05

          end if

          voll = voll + delvol(ianataldata,ijaxs)

          do ip = 1, np
             vn = anataldata(ip,ianataldata,ijaxs,1)
     &            * delvol(ianataldata,ijaxs)
             tott(ip,1) = tott(ip,1) + vn
          end do

         end do                   ! ijaxs = 1, nijaxs

         if( itaxs(m,iax) .eq. 1 ) then ! when axis = eng
          do ip = 1, np
             tott(ip,1) = tott(ip,1)
          end do

         else if( itaxs(m,iax) .eq. 14 ) then ! when axis = tet
          do ip = 1, np
           if( itunt(m) .eq. 1 .or. itunt(m) .eq. 5 ) then
              tott(ip,1) = tott(ip,1) / voll
           end if
          end do

         else if( itaxs(m,iax) .eq. 11 ) then ! when axis = time
          do ip = 1, np
           if( itunt(m) .eq. 4 ) then
              tott(ip,1) = tott(ip,1) / voll
           end if
          end do

         end if

! sumover
            do ip = 1, np
            do irst=1,2
              tott(ip,irst) = anataldata(ip,ianataldata,njaxs,irst)
            enddo
            enddo

         if( itunt(m) .ne. 0 ) then
            write(iot,'(/"#   sum over",14x,1000(1pe13.4,0pf8.4))')
     &           (tott(ip,1),tott(ip,2),ip=1,np) ! frtati 2021/10/05
         end if

         if( itaxs(m,iax).eq.1 .or. itaxs(m,iax).eq.11 ) then ! when axis = eng,t
          if( itunt(m) .ne. 0 ) then
           if ( manatally .eq. 0 ) then ! user defined analysis
              write(iot,'(/"#   sum over",14x,1000(1pe13.4,0pf8.4))')
     &             (tott(ip,1),tott(ip,2),ip=1,np) ! frtati 2021/10/05
           else if ( manatally .eq. 1 ) then ! systematic uncertainty
           write(iot,'(/"#   sum over",14x,1000(1pe13.4,0pf8.4,16x))')
     &             (tott(ip,1),tott(ip,2),ip=1,np) ! frtati 2021/10/05
           else if ( manatally .eq. 2 ) then ! c-value dependence
           end if
          end if

         else if( itaxs(m,iax).eq.14 ) then ! when axis = tet
          if( itunt(m) .eq. 1 .or. itunt(m) .eq. 5 ) then
           if ( manatally .eq. 0 ) then ! user defined analysis
              write(iot,'(/"#   sum over ",1pe13.4,
     &                     1000(1pe13.4,0pf8.4))')
     &                     voll, (tott(ip,1),tott(ip,2),ip=1,np) ! frtati 2021/10/05
           else if ( manatally .eq. 1 ) then ! systematic uncertainty
              write(iot,'(/"#   sum over ",1pe13.4,
     &                     1000(1pe13.4,0pf8.4,16x))')
     &                     voll, (tott(ip,1),tott(ip,2),ip=1,np) ! frtati 2021/10/05
           else if ( manatally .eq. 2 ) then ! c-value dependence
           end if

          else if( itunt(m) .ne. 0 ) then
           if ( manatally .eq. 0 ) then ! user defined analysis
              write(iot,'(/"#   sum over ",13x,
     &                     1000(1pe13.4,0pf8.4))')
     &                     (tott(ip,1),tott(ip,2),ip=1,np) ! frtati 2021/10/05
           else if ( manatally .eq. 1 ) then ! systematic uncertainty
              write(iot,'(/"#   sum over ",13x,
     &                     1000(1pe13.4,0pf8.4,16x))')
     &                     (tott(ip,1),tott(ip,2),ip=1,np) ! frtati 2021/10/05
           else if ( manatally .eq. 2 ) then ! c-value dependence
           end if
          end if
         end if

*-----------------------------------------------------------------------

         write(iot,'(/a1,"no. =",i3,a1)') cha, inum, cha

         write(iot,'("msuc: {",a1,"huge ",80a1)')
     &        yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

         write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
         write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

C write tally condition of each figure
         write(iot,'("wt: s(0.7)",/a1,"vspace{-3}")') yen
         itmp = 0
         do ijtmp = 2,8
          if (cij(ijtmp) .ne. 'F' ) then
           itmp = itmp + ij(ijtmp)
           if ( ibin(ijtmp) .eq. 0 ) then
              write(iot,'(a15,"&=&",i5)')
     &             cij(ijtmp), idnint(fg(itmp))
              itmp = itmp - ij(ijtmp) + nij(ijtmp)
           else
              write(iot,'(a15,2x,"(",i5,"-th bin)"
     &                             /1pe13.4,2x,"$--$",1pe13.4)')
     &             cij(ijtmp),ij(ijtmp),fg(itmp),fg(itmp+1)
              itmp = itmp - ij(ijtmp) + nij(ijtmp)+1
           end if
          end if
         end do

         write(iot,'("e:")')

        end do                    ! ianataldata = 1, nanataldata

        deallocate( fg )
        deallocate( fgaxs )
        deallocate( anataldata )
        deallocate( delvol )

*-----------------------------------------------------------------------
*        t-eng, eng-t axis (matrix)
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 12
     &       .or. itaxs(m,iax) .eq. 13 ) then

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
          select case( itaxs(m,iax) )
           case ( 12 )
            iDaxis = 35   ! t-eng
           case ( 13 )
            iDaxis = 36   ! eng-t
          end select

C rearrange data from anatalrst to anataldata
          call anatal_rearrange_sum(0,0,4,
     &         iDaxis,ntaxis,anatalrst,
     &         np,     ne,       0,    nt,     nr,  0,  0,   0, nrst,
     &         eb,ew,    1,1,  tb,tw, dlr,  dt_one,  dt_one,1,
     &         vl,
     &         nij,ibin,nijaxs,nijaxs2,nijaxs3,nfgmax,
     &         cij,cijaxs,cijaxs2,nanataldata,ierr)

          inum = 0

          do ip = 1, np
          do ianataldata = 1, nanataldata

           inum = inum + 1

C exchange ID number from ianataldata to ij2,...,ij8
           itmpdata = int(ianataldata-1)
           ij(8) = mod(itmpdata,nij(8))+1
           itmpdata = int((ianataldata-ij(8))/nij(8))
           ij(7) = mod(itmpdata,nij(7))+1
           itmpdata = int((itmpdata-ij(7)+1)/nij(7))
           ij(6) = mod(itmpdata,nij(6))+1
           itmpdata = int((itmpdata-ij(6)+1)/nij(6))
           ij(5) = mod(itmpdata,nij(5))+1
           itmpdata = int((itmpdata-ij(5)+1)/nij(5))
           ij(4) = mod(itmpdata,nij(4))+1
           itmpdata = int((itmpdata-ij(4)+1)/nij(4))
           ij(3) = mod(itmpdata,nij(3))+1
           itmpdata = int((itmpdata-ij(3)+1)/nij(3))
           ij(2) = mod(itmpdata,nij(2))+1

*-----------------------------------------------------------------------

           write(iot,'(/"#",78("-"))')

           if( inum .eq. 1 .or. ip.gt.np_mxang ) then ! frtati 2021/10/05
              write(iot,'( "#newpage:")')
           else
              write(iot,'( " newpage:")')
           end if

           if ( ip.gt.np_mxang ) then
             write(iot,'( " SKIPPAGE:")')
             inum = inum - 1
           end if

*-----------------------------------------------------------------------

           if( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) then

              write(iot,'(/a1,"no. =",i3,a1)') cha, inum, cha

              write(iot,'("msuc: {",a1,"huge ",80a1)')
     &             yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

              write(iot,'("msdr: {",a1,"it plotted by ",
     &                  a1,"ANGEL ",a1,"version}")') yen, yen, yen
              write(iot,'("msdl: {",a1,"it calculated by ",
     &                  a1,"PHITS ",f5.2"}")') yen, yen, versn

           if( lcx_txt(m) .gt. 0 ) then
              write(iot,'(/"x: ",200a1)')(cx_txt(m)(i:i),i=1,lcx_txt(m))
           else
              if( itaxl(m) .eq. 0 ) then
                 write(iot,'(/"x: ",a29)') cijaxs
              else
                 write(iot,'(/"x: ",200a1)')
     &                (itaxt(m)(i:i),i=1,itaxl(m))
              end if
           end if

           if( lcy_txt(m) .gt. 0 ) then
              write(iot,'( "y: ",200a1)')(cy_txt(m)(i:i),i=1,lcy_txt(m))
           else
              if( itayl(m) .eq. 0 ) then
                 write(iot,'("y: ",a29)') cijaxs2
              else
                 write(iot,'( "y: ",200a1)')
     &                (itayt(m)(i:i),i=1,itayl(m))
              end if
           end if

*-----------------------------------------------------------------------

              if( inum .eq. 1 ) then

                 form  = ( fgaxs2(nijaxs2+1) - fgaxs2(1) )
     &                / ( fgaxs(nijaxs+1) - fgaxs(1) )
                 xfac  = 0.9
                 afac  = 0.8
                 izlog = 1
                 inocm = 1
                 inolg = 1

                 if( itanl(m) .gt. 0 )
     &                call anset(itang(m),itanl(m),
     &                form,xfac,afac,cmin,cmax,izlog,inocm,inolg)

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
     &             ( ittwo(m) .ge. 6 .and. ittwo(m) .le. 7 ) ) .and.
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

              axs1min = fgaxs(1)
              axs1max = fgaxs(nijaxs+1)
              axs2min = fgaxs2(1)
              axs2max = fgaxs2(nijaxs2+1)
              axs1del = (fgaxs(nijaxs+1)-fgaxs(1))/dble(nijaxs)
              axs2del = (fgaxs2(nijaxs2+1)-fgaxs2(1))/dble(nijaxs2)

              write(iot,'( "p: xmin(",1p1g14.7,") xmax(",1p1g14.7,
     &                       ")")') axs1min, axs1max

              write(iot,'( "p: ymin(",1p1g14.7,") ymax(",1p1g14.7,
     &                       ")")') axs2min, axs2max

              if( itanl(m) .gt. 0 .and. ioe .eq. 1 ) then
                 write(iot,'( "p: ",200a1)')
     &                ( itang(m)(i:i),i = 1, itanl(m) )

              else if( itanl(m) .gt. 0 .and. ioe .ne. 1 ) then

                 call terrang(iot,m)

              end if

              if( itsans(m) .gt. 0 ) then
                 call write_sangel(iot,m,0)
              end if

           end if

*-----------------------------------------------------------------------

           write(iot,'("#  ny = ",i3,"   nx = ",i3)')
     &          nijaxs2, nijaxs

           if( ittwo(m) .ne. 4 ) then
              write(iot,'( "# ( ( data(i1,i2), i1 = 1, n1 ),",
     &                         " i2 = n2, 1, -1 )")')
           end if

*-----------------------------------------------------------------------

           if( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) then

              write(iot,'(/a4," y = ",1p1g14.7," to ",
     &            1p1g14.7," by ",1p1g14.7," ; x = ",
     &            1p1g14.7," to ",1p1g14.7," by ",
     &            1p1g14.7," ;")') dc2,
     &            axs2max - axs2del/2d0, axs2min + axs2del/2d0, axs2del,
     &            axs1min + axs1del/2d0, axs1max - axs1del/2d0, axs1del

!OBINATA(2012.8.3): output *.err
              write(iot,'(1p10e11.3)')
     &       ( ( anataldata(ip,ianataldata,ijaxs+(ijaxs2-1)*nijaxs,ioe),
     &                ijaxs = 1, nijaxs ), ijaxs2 = nijaxs2, 1, -1 )

*-----------------------------------------------------------------------

           else if( ittwo(m) .eq. 4 ) then

               write(iot,'(/"# axis2      axis1    ",
     &                      "  number     r.err")')

               do ijaxs = 1, nijaxs
                do ijaxs2 = 1, nijaxs2

                   write(iot,'(1p3e13.4,0pf8.4)')
     &               fgaxs2(ijaxs2)  + axs2del/2d0,
     &               fgaxs(ijaxs)  + axs1del/2d0,
     &             anataldata(ip,ianataldata,ijaxs+(ijaxs2-1)*nijaxs,1),
     &             anataldata(ip,ianataldata,ijaxs+(ijaxs2-1)*nijaxs,2)

                end do
               end do

*-----------------------------------------------------------------------

            else if( ittwo(m) .eq. 5 ) then

               write(iot,'("#   t = ",1p1g14.7," to ",
     &            1p1g14.7," by ",1p1g14.7/
     &                        "#   e = ",1p1g14.7," to ",
     &            1p1g14.7," by ",1p1g14.7)')
     &            axs2min + axs2del/2d0, axs2max - axs2del/2d0, axs2del,
     &            axs1min + axs1del/2d0, axs1max - axs1del/2d0, axs1del

               write(iot,'(/4x,a3,4x,1p1000e11.3)')
     &          'ax2/ax1',( fgaxs(ijaxs) + axs1del/2d0, ijaxs=1,nijaxs )

               do ijaxs2 = nijaxs2, 1, -1

!OBINATA(2012.8.3): output *.err
                  write(iot,'(1p1000e11.3)')
     &                 fgaxs2(ijaxs2)  + axs2del/2d0,
     &         ( anataldata(ip,ianataldata,ijaxs+(ijaxs2-1)*nijaxs,ioe),
     &                 ijaxs = 1, nijaxs )

               end do

            end if

*-----------------------------------------------------------------------

      if( inolg .eq. 1 .and. ( ittwo(m) .eq. 3 .or. ittwo(m) .eq. 7 )
     &    .and. cmin .gt. 0.0 .and. cmax .gt. cmin ) then

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

      else

       if( lcz_txt(m) .gt. 0 ) then
          write(iot,'("y: ",200a1)')(cz_txt(m)(i:i),i=1,lcz_txt(m))
       else
        if( itazl(m) .eq. 0 ) then
            write(iot,'("y: ",a22)') hsunit(itunt(m))
        else
            write(iot,'("y: ",200a1)') (itazt(m)(i:i),i=1,itazl(m))
        end if
       end if

      end if

      if( izlog .eq. 0 ) write(iot,'( "p: ylin")')

      write(iot,'( "#",78("-"))')

      end if

*-----------------------------------------------------------------------

         if( inocm .eq. 1 ) then

          write(iot,'("wt: s[c5]",/a1,"vspace{-3}")') yen
          itmp = 0
          do ijtmp = 2,8
           if (cij(ijtmp) .ne. 'F' ) then
              itmp = itmp + ij(ijtmp)
              if ( ibin(ijtmp) .eq. 0 ) then
                 write(iot,'(a15,"&=&",i5)')
     &                cij(ijtmp), idnint(fg(itmp))
                 itmp = itmp - ij(ijtmp) + nij(ijtmp)
              else
                 write(iot,'(a15,2x,"(",i5,"-th bin)"
     &                             /1pe13.4,2x,"$--$",1pe13.4)')
     &                cij(ijtmp),ij(ijtmp),fg(itmp),fg(itmp+1)
                 itmp = itmp - ij(ijtmp) + nij(ijtmp)+1
              end if
             end if
            end do
            write(iot,'("Particle = ",a8)') chq(ip)

            write(iot,'("e:")')

           end if

*-----------------------------------------------------------------------

      if( inolg .eq. 1 .and. ( ittwo(m) .eq. 3 .or. ittwo(m) .eq. 7 )
     &    .and. cmin .gt. 0.0 .and. cmax .gt. cmin ) then

         write(iot,'(/"#",78("-"))')
         write(iot,'("z: xorg[-1.03/0.05]"/
     &               "p: form[c1] xfac[c2] afac[c5] nosp notf")')

               write(iot,'( "p: xmin(",1p1g14.7,") xmax(",1p1g14.7,
     &                       ")")') axs1min, axs1max

               write(iot,'( "p: ymin(",1p1g14.7,") ymax(",1p1g14.7,
     &                       ")")') axs2min, axs2max

      end if

*-----------------------------------------------------------------------

          end do   ! ip = 1, np
          end do   ! ianataldata

          deallocate( fg )
          deallocate( fgaxs )
          deallocate( fgaxs2 )
          deallocate( fgaxs3 )
          deallocate( anataldata )
          deallocate( delvol )

*-----------------------------------------------------------------------

         end if   ! itaxs(m,iax)

         close(iot)

         if( iteps(m) .ne. 0 ) then

            open(iot, file = fname, status = 'unknown' )
            call a_angel(idasa,fname)

         end if

*-----------------------------------------------------------------------

  900 continue     ! ioe, iax

*-----------------------------------------------------------------------
      deallocate( ixyz )
      deallocate( lr,vl,val ) !FURUTA20190204
      if(itfoam(m).eq.2) deallocate(xcm) !FURUTA20191028

      deallocate (anatalrst)
      deallocate (tw,rdata,dlr)

      return
      end subroutine anatal_pdepsttet



!***********************************************************************
!                                                                      *
      subroutine anatal_calc_deposit_anova_sub(m, ntf,np,
     &                                        resc2,resc3 )
!                                                                      *
!***********************************************************************

      use TALMOD
!$      use TALMOD0
      use anatallymod
! lanatalm_sum(m)     :  start address of one tally at nfile(m)=1
! ianatalm_sum(m,iax,nfile(m)) :  start address of one sumover
! manatalm_sum(m,iax,nfile(m)) :  one sumover length

      implicit none

C for nonshared_tally option
      integer italsh
      common /talsh/ italsh

      integer :: m, ntf, np

      double precision :: resc2,resc3

      real(8),pointer :: p_sum(:)
      real(8),pointer :: tran_sum(:)

      integer :: ln_sum,iax

      do iax =1,6

C for nonshared_tally option
!$       if(italsh .eq. 0) then
!$          call GET_TR_HEAD_POINTER0_SUM(p_sum,m,iax)
!$       else
            call GET_TR_HEAD_POINTER_SUM(p_sum,m,iax)
C for nonshared_tally option
!$       end if

        tran_sum => trANATAL_SUM(ianatalm_sum(m,iax,ntf):)
        ln_sum = manatalm_sum(m,iax,ntf)/2

        if(ln_sum > 0 ) then 
          call anatal_calc_deposit_anova_2trAN(m,np,p_sum,tran_sum,
     &       ln_sum, resc2,resc3)

        endif
      enddo

      return
      end

!***********************************************************************
!                                                                      *
! sumover subroutine                                                   *
!                                                                      *
!***********************************************************************


!***********************************************************************
!                                                                      *
      subroutine anatal_calc_deposit_anova_2trAN(m,np,
     &       tr_sum,trANATAL_SUM, ln_sum, resc2,resc3)
!                                                                      *
!***********************************************************************

      implicit none

      integer :: m,np,ln_sum
      real(8) :: tr_sum(ln_sum,2),trANATAL_SUM(ln_sum,2)
      real(8) :: resc2,resc3
      integer :: i,ip
      real(8) :: Xa,sigx

      ip = 0
      do i=1,ln_sum
             ip = ip + 1
             if(ip > np) then
               ip = 1
             endif
             if( tr_sum(i,1) .gt. 0.d0 ) then
              call calc_deposit_stdev(m,Xa,sigx,
     &                        tr_sum(i,1),
     &                        tr_sum(i,2),
     &                        1.0d+0,ip)
              tr_sum(i,1) = Xa
              tr_sum(i,2) = sigx
             end if
      enddo

! X_bar = Xj_bar
      trANATAL_SUM(1:ln_sum,1) = tr_sum(1:ln_sum,1)

! (sig_xj)**2
      trANATAL_SUM(1:ln_sum,2) = ( tr_sum(1:ln_sum,2)
     &                         *   tr_sum(1:ln_sum,1) ) **2

! sig_x
      trANATAL_SUM(1:ln_sum,2) =
     &                          dsqrt( trANATAL_SUM(1:ln_sum,2))

! Sigma(xi wi)**2 = (sig_X**2 N(N-1)+N X_bar**2)(W/N)**2
      trANATAL_SUM(1:ln_sum,2) = ( trANATAL_SUM(1:ln_sum,2)**2
     &                 * resc3 * (resc3-1.0d0)
     &                 + resc3 * trANATAL_SUM(1:ln_sum,1)**2 )
     &                 * (resc2/resc3)**2

! Sigma xi wi = X_bar W
      trANATAL_SUM(1:ln_sum,1) =
     &             trANATAL_SUM(1:ln_sum,1) * resc2

      return
      end

!***********************************************************************
!                                                                      *
! sumover subroutine                                                   *
!                                                                      *
!***********************************************************************

!***********************************************************************
!                                                                      *
      subroutine get_depstreg_tr_sum_data(m,iax,nfile,
     &           ip,ie,ir,it,rdata)
!                                                                      *
!***********************************************************************

      use TALMOD
!$      use TALMOD0
      use anatallymod, only:ianatalm_sum
! lanatalm_sum(m)     :  start address of one tally at nfile(m)=1
! ianatalm_sum(m,iax,nfile(m)) :  start address of one sumover
! manatalm_sum(m,iax,nfile(m)) :  one sumover length

      implicit double precision (a-h,o-z)


      integer :: m,iax,nfile,ip,ie,ir,it
      real(8) :: rdata(2,nfile)

      integer :: ntf

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

      real(8),pointer :: p_sum(:)

        do ntf=1,nfile

C for nonshared_tally option
!$       if(italsh .eq. 0) then
!$         p_sum => tr00_sum(ianatalm_sum(m,iax,ntf):)
!$       else
           p_sum => tr0_sum(ianatalm_sum(m,iax,ntf):)
C for nonshared_tally option
!$       end if

        call get_depstreg_tr_sum_data_sub(p_sum,
     &     itpan_sum(m,iax),itenm_sum(m,iax),
     &     itrgn_sum(m,iax),ittnm_sum(m,iax),
     &     ip,ie,ir,it,rdata(1,ntf))

      enddo

      return
      end

!***********************************************************************
!                                                                      *
      subroutine get_depstreg_tr_sum_data_sub(tr_sum,
     &          np,ne,nr,nt,ip,ie,ir,it,rdata)
!                                                                      *
!***********************************************************************

      implicit double precision (a-h,o-z)

      integer :: m,iax,np,ne,nr,nt,ip,ie,ir,it
      real(8) :: rdata(2)

      real(8) :: tr_sum(np,0:ne,nr,nt,2)
 
      do i=1,2
        rdata(i)= tr_sum(ip,ie,ir,it,i)
      enddo

      end

!***********************************************************************
!                                                                      *
      subroutine get_depstrz_tr_sum_data(m,iax,nfile,
     &           ip,ie,it,ir,iz,rdata)
!                                                                      *
!***********************************************************************

      use TALMOD
!$      use TALMOD0
      use anatallymod, only:ianatalm_sum
! lanatalm_sum(m)     :  start address of one tally at nfile(m)=1
! ianatalm_sum(m,iax,nfile(m)) :  start address of one sumover
! manatalm_sum(m,iax,nfile(m)) :  one sumover length

      implicit double precision (a-h,o-z)


      integer :: m,iax,nfile,ip,ie,it,ir,iz
      real(8) :: rdata(2,nfile)

      integer :: ntf

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

      real(8),pointer :: p_sum(:)

        do ntf=1,nfile

C for nonshared_tally option
!$       if(italsh .eq. 0) then
!$         p_sum => tr00_sum(ianatalm_sum(m,iax,ntf):)
!$       else
           p_sum => tr0_sum(ianatalm_sum(m,iax,ntf):)
C for nonshared_tally option
!$       end if

        call get_depstrz_tr_sum_data_sub(p_sum,
     &     itpan_sum(m,iax),itenm_sum(m,iax),ittnm_sum(m,iax),
     &     itrnm_sum(m,iax),itznm_sum(m,iax),
     &     ip,ie,it,ir,iz,rdata(1,ntf))

      enddo

      return
      end

!***********************************************************************
!                                                                      *
      subroutine get_depstrz_tr_sum_data_sub(tr_sum,
     &          np,ne,nt,nr,nz,ip,ie,it,ir,iz,rdata)
!                                                                      *
!***********************************************************************

      implicit double precision (a-h,o-z)

      integer :: m,iax,np,ne,nt,nr,nz,ip,ie,it,ir,iz
      real(8) :: rdata(2)

      real(8) :: tr_sum(np,0:ne,nt,nr,nz,2)

      do i=1,2
        rdata(i)= tr_sum(ip,ie,it,ir,iz,i)
      enddo

      end

!***********************************************************************
!                                                                      *
      subroutine get_depstxyz_tr_sum_data(m,iax,nfile,
     &           ip,ie,it,ix,iy,iz,rdata)
!                                                                      *
!***********************************************************************

      use TALMOD
!$      use TALMOD0
      use anatallymod, only:ianatalm_sum
! lanatalm_sum(m)     :  start address of one tally at nfile(m)=1
! ianatalm_sum(m,iax,nfile(m)) :  start address of one sumover
! manatalm_sum(m,iax,nfile(m)) :  one sumover length

      implicit double precision (a-h,o-z)


      integer :: m,iax,nfile,ip,ie,it,ix,iy,iz
      real(8) :: rdata(2,nfile)

      integer :: ntf

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

      real(8),pointer :: p_sum(:)

        do ntf=1,nfile

C for nonshared_tally option
!$       if(italsh .eq. 0) then
!$         p_sum => tr00_sum(ianatalm_sum(m,iax,ntf):)
!$       else
           p_sum => tr0_sum(ianatalm_sum(m,iax,ntf):)
C for nonshared_tally option
!$       end if

        call get_depstxyz_tr_sum_data_sub(p_sum,
     &     itpan_sum(m,iax),itenm_sum(m,iax),ittnm_sum(m,iax),
     &     itxnm_sum(m,iax),itynm_sum(m,iax),itznm_sum(m,iax),
     &     ip,ie,it,ix,iy,iz,rdata(1,ntf))

      enddo

      return
      end

!***********************************************************************
!                                                                      *
      subroutine get_depstxyz_tr_sum_data_sub(tr_sum,
     &          np,ne,nt,nx,ny,nz,ip,ie,it,ix,iy,iz,rdata)
!                                                                      *
!***********************************************************************

      implicit double precision (a-h,o-z)

      integer :: m,iax,np,ne,nt,nx,ny,nz,nm,ip,ie,it,ix,iy,iz
      real(8) :: rdata(2)

      real(8) :: tr_sum(np,0:ne,nt,nx*ny*nz,2)

      icf(ix,iy,iz) = ix + ( iy - 1 ) * nx + ( iz - 1 ) * nx * ny


      do i=1,2
        rdata(i)= tr_sum(ip,ie,it,icf(ix,iy,iz),i)
      enddo

      end
!***********************************************************************
