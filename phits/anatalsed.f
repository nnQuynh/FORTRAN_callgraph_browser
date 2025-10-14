!***********************************************************************
!                                                                      *
      subroutine anatal_calc_anova_sedreg(m,ntf,
     &                             np,ne,nr,
     &                             sumfactor,weightRate,
     &                             trRES,tranatal,nfile,
     &                             resc2SUMTAL,resc3SUMTAL,
     &                             ierr)
!                                                                      *
!     created by T.Miura on 2021/09/10                                 *
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
      integer           nr

      double precision  sumfactor     ! normalization factor
      double precision  weightRate(*) ! weighting rate
      double precision  trRES   (np,ne,nr,2)
      integer           nfile
      double precision  tranatal(np,ne,nr,2*nfile)
      double precision  resc2SUMTAL
      double precision  resc3SUMTAL

      double precision  cmin, cmax

      integer           ierr
!-----------------------------------------------------------------------
      integer           ir
      integer           ie
      integer           ip

      double precision  Xa, sigx
!-----------------------------------------------------------------------
      integer           iat, iad, iaf

      iat(iad,iaf) = iad + (iaf-1) * 2
!-----------------------------------------------------------------------

      ierr = 0

*-----------------------------------------------------------------------
      cmax = 0.0
      cmin = 1.e+33

      do ir = 1, nr
       do ie = 1, ne
        do ip = 1, np

             if( trRES(ip,ie,ir,1) .gt. 0.d0 ) then
              call calc_stdev(m,Xa,sigx,
     &                        trRES(ip,ie,ir,1),
     &                        trRES(ip,ie,ir,2),
     &                        1.0d+0)
              trRES(ip,ie,ir,1) = Xa
              trRES(ip,ie,ir,2) = sigx
             end if

              ! X_bar = Xj_bar
              tranatal(ip,ie,ir,iat(1,ntf)) =
     &           trRES(ip,ie,ir,1)

              ! (sig_xj)**2
              tranatal(ip,ie,ir,iat(2,ntf)) =
     &        (  trRES(ip,ie,ir,2)
     &         * trRES(ip,ie,ir,1) )**2

              ! sig_x
              tranatal(ip,ie,ir,iat(2,ntf)) = dsqrt(
     &        tranatal(ip,ie,ir,iat(2,ntf)))

              ! Sigma(xi wi)**2 = (sig_X**2 N(N-1)+N X_bar**2)(W/N)**2
              tranatal(ip,ie,ir,iat(2,ntf)) =
     &       (tranatal(ip,ie,ir,iat(2,ntf))**2 *
     &        resc3(m) * (resc3(m)-1.0d0) +
     &        resc3(m) *
     &        tranatal(ip,ie,ir,iat(1,ntf))**2) *
     &        (resc2(m)/resc3(m))**2

              ! Sigma xi wi = X_bar W
              tranatal(ip,ie,ir,iat(1,ntf)) =
     &        tranatal(ip,ie,ir,iat(1,ntf)) * resc2(m)

            if( tranatal(ip,ie,ir,iat(1,ntf)) > cmax )
     &            cmax = tranatal(ip,ie,ir,iat(1,ntf))
            if( tranatal(ip,ie,ir,iat(1,ntf)) < cmin )
     &            cmin = tranatal(ip,ie,ir,iat(1,ntf))

        end do     ! ip loop end
       end do      ! ie loop end
      end do       ! ir loop end

! sumover
      call anatal_calc_anova_sub(m,ntf,
     &                               resc2(m),resc3(m) )

!-----------------------------------------------------------------------
      return
      end subroutine anatal_calc_anova_sedreg



!***********************************************************************
!                                                                      *
      subroutine anatal_calc_anova_sedrz(m,ntf,
     &                             np,ne,nr,nz,
     &                             sumfactor,weightRate,
     &                             trRES,tranatal,nfile,
     &                             resc2SUMTAL,resc3SUMTAL,
     &                             ierr)
!                                                                      *
!     created by T.Miura on 2021/09/10                                 *
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
      integer           nr
      integer           nz

      double precision  sumfactor     ! normalization factor
      double precision  weightRate(*) ! weighting rate
      double precision  trRES   (np,ne,nr,nz,2)
      integer           nfile
      double precision  tranatal(np,ne,nr,nz,2*nfile)
      double precision  resc2SUMTAL
      double precision  resc3SUMTAL

      double precision  cmin, cmax

      integer           ierr
!-----------------------------------------------------------------------
      integer           iz
      integer           ir
      integer           ie
      integer           ip

      double precision  Xa, sigx
!-----------------------------------------------------------------------
      integer           iat, iad, iaf

      iat(iad,iaf) = iad + (iaf-1) * 2
!-----------------------------------------------------------------------

      ierr = 0

*-----------------------------------------------------------------------
      cmax = 0.0
      cmin = 1.e+33

      do iz = 1, nz
       do ir = 1, nr
        do ie = 1, ne
         do ip = 1, np

             if( trRES(ip,ie,ir,iz,1) .gt. 0.d0 ) then
              call calc_stdev(m,Xa,sigx,
     &                        trRES(ip,ie,ir,iz,1),
     &                        trRES(ip,ie,ir,iz,2),
     &                        1.0d+0)
              trRES(ip,ie,ir,iz,1) = Xa
              trRES(ip,ie,ir,iz,2) = sigx
             end if

              ! X_bar = Xj_bar
              tranatal(ip,ie,ir,iz,iat(1,ntf)) =
     &           trRES(ip,ie,ir,iz,1)

              ! (sig_xj)**2
              tranatal(ip,ie,ir,iz,iat(2,ntf)) =
     &        (  trRES(ip,ie,ir,iz,2)
     &         * trRES(ip,ie,ir,iz,1) )**2

              ! sig_x
              tranatal(ip,ie,ir,iz,iat(2,ntf)) = dsqrt(
     &        tranatal(ip,ie,ir,iz,iat(2,ntf)))

              ! Sigma(xi wi)**2 = (sig_X**2 N(N-1)+N X_bar**2)(W/N)**2
              tranatal(ip,ie,ir,iz,iat(2,ntf)) =
     &       (tranatal(ip,ie,ir,iz,iat(2,ntf))**2 *
     &        resc3(m) * (resc3(m)-1.0d0) +
     &        resc3(m) *
     &        tranatal(ip,ie,ir,iz,iat(1,ntf))**2) *
     &        (resc2(m)/resc3(m))**2

              ! Sigma xi wi = X_bar W
              tranatal(ip,ie,ir,iz,iat(1,ntf)) =
     &        tranatal(ip,ie,ir,iz,iat(1,ntf)) * resc2(m)

            if( tranatal(ip,ie,ir,iz,iat(1,ntf)) > cmax )
     &            cmax = tranatal(ip,ie,ir,iz,iat(1,ntf))
            if( tranatal(ip,ie,ir,iz,iat(1,ntf)) < cmin )
     &            cmin = tranatal(ip,ie,ir,iz,iat(1,ntf))

         end do     ! ip loop end
        end do      ! ie loop end
       end do       ! ir loop end
      end do        ! iz loop end

! sumover
      call anatal_calc_anova_sub(m,ntf,
     &                               resc2(m),resc3(m) )

!-----------------------------------------------------------------------
      return
      end subroutine anatal_calc_anova_sedrz



!***********************************************************************
!                                                                      *
      subroutine anatal_calc_anova_sedxyz(m,ntf,
     &                             np,ne,nx,ny,nz,
     &                             sumfactor,weightRate,
     &                             trRES,tranatal,nfile,
     &                             resc2SUMTAL,resc3SUMTAL,
     &                             ierr)
!                                                                      *
!     created by T.Miura on 2021/09/10                                 *
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
      integer           nx
      integer           ny
      integer           nz

      double precision  sumfactor     ! normalization factor
      double precision  weightRate(*) ! weighting rate
      double precision  trRES   (np,ne,nx*ny*nz,2)
      integer           nfile
      double precision  tranatal(np,ne,nx*ny*nz,2*nfile)
      double precision  resc2SUMTAL
      double precision  resc3SUMTAL

      double precision  cmin, cmax

      integer           ierr
!-----------------------------------------------------------------------
      integer           iz
      integer           iy
      integer           ix
      integer           ie
      integer           ip

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
         do ie = 1, ne
          do ip = 1, np

             if( trRES(ip,ie,icf(ix,iy,iz),1) .gt. 0.d0 ) then
              call calc_stdev(m,Xa,sigx,
     &                        trRES(ip,ie,icf(ix,iy,iz),1),
     &                        trRES(ip,ie,icf(ix,iy,iz),2),
     &                        1.0d+0)
              trRES(ip,ie,icf(ix,iy,iz),1) = Xa
              trRES(ip,ie,icf(ix,iy,iz),2) = sigx
             end if

              ! X_bar = Xj_bar
              tranatal(ip,ie,icf(ix,iy,iz),iat(1,ntf)) =
     &           trRES(ip,ie,icf(ix,iy,iz),1)

              ! (sig_xj)**2
              tranatal(ip,ie,icf(ix,iy,iz),iat(2,ntf)) =
     &        (  trRES(ip,ie,icf(ix,iy,iz),2)
     &         * trRES(ip,ie,icf(ix,iy,iz),1) )**2

              ! sig_x
              tranatal(ip,ie,icf(ix,iy,iz),iat(2,ntf)) = dsqrt(
     &        tranatal(ip,ie,icf(ix,iy,iz),iat(2,ntf)))

              ! Sigma(xi wi)**2 = (sig_X**2 N(N-1)+N X_bar**2)(W/N)**2
              tranatal(ip,ie,icf(ix,iy,iz),iat(2,ntf)) =
     &       (tranatal(ip,ie,icf(ix,iy,iz),iat(2,ntf))**2 *
     &        resc3(m) * (resc3(m)-1.0d0) +
     &        resc3(m) *
     &        tranatal(ip,ie,icf(ix,iy,iz),iat(1,ntf))**2) *
     &        (resc2(m)/resc3(m))**2

              ! Sigma xi wi = X_bar W
              tranatal(ip,ie,icf(ix,iy,iz),iat(1,ntf)) =
     &        tranatal(ip,ie,icf(ix,iy,iz),iat(1,ntf)) * resc2(m)

            if( tranatal(ip,ie,icf(ix,iy,iz),iat(1,ntf)) > cmax )
     &            cmax = tranatal(ip,ie,icf(ix,iy,iz),iat(1,ntf))
            if( tranatal(ip,ie,icf(ix,iy,iz),iat(1,ntf)) < cmin )
     &            cmin = tranatal(ip,ie,icf(ix,iy,iz),iat(1,ntf))

          end do    ! ip loop end
         end do     ! ie loop end
        end do      ! ix loop end
       end do       ! iy loop end
      end do        ! iz loop end

! sumover
      call anatal_calc_anova_sub(m,ntf,
     &                               resc2(m),resc3(m) )

!-----------------------------------------------------------------------
      return
      end subroutine anatal_calc_anova_sedxyz



************************************************************************
*                                                                      *
      subroutine anatal_psedreg(m,np,nr,mr,ne,kr,eb,tr,
     &                   nfile,weightRate,
     &                   nvl,ivl,rvl,
     &                   nx,ny,nz,xm,ym,zm,igsh,idasa,manatally)
*                                                                      *
*       output the SED tally in region mesh                            *
*       created by T.Miura on 2021/09/13                               *
*                                                                      *
************************************************************************
      use sangelmod, only: itsans
      use anatallymod, only: anataldata,fg,fgaxs,fgaxs2,fgaxs3,delvol,
     &     cx_txt,lcx_txt,cy_txt,lcy_txt,cz_txt,lcz_txt
      use partmod ! frtati 2021/10/05

*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

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

      common /tall57/ itsun(itlmax),rtdim(itlmax),rtucv(itlmax),
     &rtrho(itlmax),itmodel(itlmax) ! T.Sato 2022/08/14

      common /istcut/ ist_cut, ist_bat
*-----------------------------------------------------------------------

      common /volreg/ dvol(kvlmax)

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

      dimension   kr(mr)
      dimension   eb(ne+1)
!      dimension   ew(ne)
!      dimension   vl(nr)
!      dimension   lr(nr)
      dimension   tr(np,ne,nr,2*nfile)
      dimension   ivl(nvl)
      dimension   rvl(nvl)

      integer     nfile, manatally
!      dimension   weightRate(nfile), rdata(2,nfile)
      dimension   weightRate(nfile)

      dimension   xm(nx+1)
      dimension   ym(ny+1)
      dimension   zm(nz+1)
!      dimension   val(nr)

      integer,allocatable :: ixyz(:)

      integer irst,nrst, itmpdata
!      dimension   anatalrst(np,ne+1,1,1,nr*1*1,1,2*nfile+3)
! sumover
!      dimension   anatalrst(np,ne+1,1,1,(nr+1)*1*1,1,2*nfile+3)
      real(8),allocatable :: anatalrst(:,:,:,:,:,:,:)
      integer ianataldata,nanataldata
      integer, parameter :: ntaxis=9  ! S.H. 2021.8.9
      integer ij(ntaxis),nij(ntaxis),ibin(ntaxis)
      integer nijaxs,nijaxs2,nijaxs3,nfgmax
      character cij(ntaxis)*21,cijaxs*21,cijaxs2*21
      integer iDaxis
!      dimension   dlr(nr+1) !  integer lr(nr) ---> double dlr(nr)
      real(8),allocatable :: ew(:),vl(:),val(:),rdata(:,:), dlr(:)
      integer,allocatable :: lr(:)

*-----------------------------------------------------------------------

      dimension tott(6,2),totnorm(6) ! totnorm is added on 2017/11/21 for normalization

      character hsunit(24)*30

      data hsunit(1) / 'T.L. [cm/MeV/source]  '/
      data hsunit(2) / 'T.L. [cm/(keV/um)/source]'/
      data hsunit(3) / 'T.L. [cm/Gy/source]  '/
      data hsunit(4) / 'Dose [MeV/MeV/source]  '/
      data hsunit(5) / 'Dose [MeV/(keV/um)/source]'/
      data hsunit(6) / 'Dose [MeV/Gy/source]  '/
      data hsunit(7) / 'T.L. [cm/ln(MeV)/source]   '/
      data hsunit(8) / 'T.L. [cm/ln(keV/um)/source]'/
      data hsunit(9) / 'T.L. [cm/ln(Gy)/source]    '/
      data hsunit(10)/ 'Dose [MeV/ln(MeV)/source]   '/
      data hsunit(11)/ 'Dose [MeV/ln(keV/um)/source]'/
      data hsunit(12)/ 'Dose [MeV/ln(Gy)/source]    '/
      data hsunit(13)/ 'T.L. [cm/source]'/
      data hsunit(14)/ 'T.L. [cm/source]'/
      data hsunit(15)/ 'T.L. [cm/source]'/
      data hsunit(16)/ 'Dose [MeV/source]'/
      data hsunit(17)/ 'Dose [MeV/source]'/
      data hsunit(18)/ 'Dose [MeV/source]'/
      data hsunit(19)/ 'E*f(E) [dimensionless]'/
      data hsunit(20)/ 'y*f(y) [dimensionless]'/
      data hsunit(21)/ 'z*f(z) [dimensionless]'/
      data hsunit(22)/ 'E*d(E) [MeV]'/
      data hsunit(23)/ 'y*d(y) [keV/um]'/
      data hsunit(24)/ 'z*d(z) [Gy]'/

      character hxunit(3)*12
      data hxunit(1) / 'z [MeV]     '/
      data hxunit(2) / 'y [keV/um]'/
      data hxunit(3) / 'z [Gy]      '/

      character cha*1
      data cha /"'"/

      character chp(itmxpt)*11  ! kitamura22/03/31
      character chq(itmxpt)*9   ! kitamura22/03/31

      character rpa*1
      data rpa /'}'/
      character yen*1

      dimension dt_one(1)
      data dt_one /1.0d0/

      integer           iat, iad, iaf
      iat(iad,iaf) = iad + (iaf-1) * 2
!-----------------------------------------------------------------------
      yen  = char(92)
      allocate( ixyz(max((nx+1)*(ny+1),(ny+1)*(nz+1),(nx+1)*(nz+1))) )
      allocate (ew(ne),vl(nr),val(nr),rdata(2,nfile), dlr(nr+1))
      allocate (lr(nr))

      nrst = 2*nfile+3 ! S.H. 2021.8.15
      allocate (anatalrst(np,ne+1,1,1,(nr+1)*1*1,1,nrst))


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

         else

            npg = 1
            neg = 1

         end if

*-----------------------------------------------------------------------
*        particle name
*-----------------------------------------------------------------------

               call wtpname(m,np,chp,chq)

*-----------------------------------------------------------------------
*           itunt(m) = 1 : cm/z/source
*                    = 2 : MeV/z/source
*                    = 3 : cm/ln(z)/source
*                    = 4 : MeV/ln(z)/source
*                    = 5 : cm/source
*                    = 6 : MeV/source
*                    = 7 : y*f(y)
*                    = 8 : y*d(y)
*           itsun(m) = 1 : MeV
*                    = 2 : keV/um
*                    = 3 : Gy
*-----------------------------------------------------------------------
            if( itunt(m) .le. 2 .or. itunt(m) .eq. 7) then ! T.Sato 2017/11/21

               do i = 1, ne

                  ew(i) = ( eb(i+1) - eb(i) )

               end do
               ew_sum = eb(ne+1) - eb(1)


            else if( itunt(m) .le. 4 .or. itunt(m) .eq. 8) then ! T.Sato 2017/11/21

               do i = 1, ne

                  ew(i) = log( eb(i+1) / eb(i) )

               end do
               ew_sum = log( eb(ne+1) / eb(1) )


            else

               do i = 1, ne

                  ew(i) = 1.d0

               end do
               ew_sum = 1.0d0

            end if

*-----------------------------------------------------------------------
*        set mesh volume and name
*-----------------------------------------------------------------------

               call tregvl(nr,mr,kr,vl,lr,nvl,ivl,rvl)
               vl_sum = sum(vl(:))

*-----------------------------------------------------------------------
*        vol = 1.0
*-----------------------------------------------------------------------

               do ir = 1, nr

                  vl(ir) = 1.0d0

               end do
               vl_sum = 1.0d0

*-----------------------------------------------------------------------

         do ir = 1, nr
            dlr(ir) = dble(lr(ir))
         end do
         dlr(nr+1) = dble(lr(nr))

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

            do 100 ir = 1, nr
            do 100 ie = 1, ne
            do 100 ip = 1, np
            do ntf = 1, nfile

               if( tr(ip,ie,ir,iat(1,ntf)) .gt. 0.d0 ) then

                  call calc_stdev(m,Xa,sigx,
     &                            tr(ip,ie,ir,iat(1,ntf)),
     &                            tr(ip,ie,ir,iat(2,ntf)),
     &                            rtfac(m)/vl(ir)/ew(ie))

                  tr(ip,ie,ir,iat(1,ntf)) = Xa
                  tr(ip,ie,ir,iat(2,ntf)) = sigx
                  if( sigx .gt. stdm ) stdm = sigx

                  if( tr(ip,ie,ir,iat(1,ntf)) .gt. cmax )
     &                        cmax = tr(ip,ie,ir,iat(1,ntf))

                  if( tr(ip,ie,ir,iat(1,ntf)) .lt. cmin )
     &                        cmin = tr(ip,ie,ir,iat(1,ntf))

               else

                  isdz = 1
                  tr(ip,ie,ir,iat(2,ntf)) = 0.0

               end if

! sumover
               call psedreg_sumover_stdev_ntf(0,m,ntf,
     &             ip,ie,ir,
     &             rtfac(m),ew(ie),vl(ir),
     &             ew_sum,vl_sum)

            end do
  100       continue

!$OMP PARALLEL
!$OMP& private(ipomp,ir,ie,ip,ntf,isdz,rdata,answer,rerr)
!$OMP& private(iax,ir_a,ie_a,ido_ana,ioe)
!$      npomp=OMP_GET_NUM_THREADS()
!$      ipomp=OMP_GET_THREAD_NUM()
!$      write(*,'(''OpenMP PARALLEL PROCESS'',
!$   &     i4,''/'',i4)') ipOMP+1,npomp
!$OMP DO SCHEDULE(dynamic) reduction(max:cmax) reduction(min:cmin)
            do 101 ir = 1, nr
            do 101 ie = 1, ne
            do 101 ip = 1, np

            if( manatally .eq. 0 ) then ! T.Sato 2020/10/22

               do ntf = 1, nfile
                  rdata(1,ntf) = tr(ip,ie,ir,iat(1,ntf))
                  rdata(2,ntf) = tr(ip,ie,ir,iat(2,ntf))
               end do

              call usranatal(nfile,rdata,answer,rerr)

              anatalrst(ip,ie,1,1,ir,1,1) = answer
              anatalrst(ip,ie,1,1,ir,1,2) = rerr

              if(answer.gt.cmax) cmax = answer
              if(answer.gt.0d0 .and. answer.lt.cmin) cmin = answer

            else if( manatally .eq. 1 ) then

               do ntf = 1, nfile
                  rdata(1,ntf) = tr(ip,ie,ir,iat(1,ntf))
                  rdata(2,ntf) = tr(ip,ie,ir,iat(2,ntf))
               end do

               call anova(nfile,rdata,resc3(m)
     &              ,fmval,unca,uncerr,uratio)

               anatalrst(ip,ie,1,1,ir,1,1) ! mean
     &              = fmval
               if ( fmval .gt. 0d0 ) then
                  anatalrst(ip,ie,1,1,ir,1,2) ! total uncertainty
     &                 = dsqrt(unca**2 + uncerr**2/resc3(m)) / fmval
                  anatalrst(ip,ie,1,1,ir,1,3) ! syst. uncertainty
     &                 = unca / fmval
                  anatalrst(ip,ie,1,1,ir,1,4) ! stat. uncertainty
     &                 = uncerr /dsqrt(resc3(m)) / fmval
                  anatalrst(ip,ie,1,1,ir,1,5) ! squared unc. ratio
     &                 = uratio
               else
                  anatalrst(ip,ie,1,1,ir,1,2) ! total uncertainty
     &                 = 0d0
                  anatalrst(ip,ie,1,1,ir,1,3) ! syst. uncertainty
     &                 = 0d0
                  anatalrst(ip,ie,1,1,ir,1,4) ! stat. uncertainty
     &                 = 0d0
                  anatalrst(ip,ie,1,1,ir,1,5) ! squared unc. ratio
     &                 = 0d0
               end if

            else if( manatally .eq. 2 ) then

               do ntf = 1, nfile
                do ioe = 1, 2
                  anatalrst(ip,ie,1,1,ir,1,iat(ioe,ntf))
     &              = tr(ip,ie,ir,iat(ioe,ntf))
                end do
               end do

            end if

!sumover
            do iax = 1,itaxn(m)

              ie_a = ie
              ir_a = ir
              ido_ana = 0
              if(itaxs(m,iax) == 1 .and. ie == 1) then
                ie_a = ne + 1
                ido_ana = 1
              elseif(itaxs(m,iax) == 2 .and. ir == 1) then
                ir_a = nr + 1
                ido_ana = 1
              endif

              if(ido_ana == 1) then

                 call get_sedreg_tr_sum_data(m,iax,nfile,
     &                 ip,ie,ir,rdata)

                if( manatally .eq. 0 ) then ! T.Sato 2020/10/22

                   call usranatal(nfile,rdata,answer,rerr)

                   anatalrst(ip,ie_a,1,1,ir_a,1,1) = answer
                   anatalrst(ip,ie_a,1,1,ir_a,1,2) = rerr

                else if( manatally .eq. 1 ) then

                  call anova(nfile,rdata,resc3(m)
     &               ,fmval,unca,uncerr,uratio)

                  anatalrst(ip,ie_a,1,1,ir_a,1,1) ! mean
     &              = fmval
                 if ( fmval .gt. 0.0d0 ) then
                  anatalrst(ip,ie_a,1,1,ir_a,1,2) ! total uncertainty
     &                 = dsqrt(unca**2 + uncerr**2/resc3(m)) / fmval
                 else
                   anatalrst(ip,ie_a,1,1,ir_a,1,2) ! total uncertainty
     &                 = 0.0d0
                 end if

                else if( manatally .eq. 2 ) then
                  do ntf = 1, nfile
                    do ioe = 1, 2
                      anatalrst(ip,ie_a,1,1,ir_a,1,iat(ioe,ntf))
     &                = rdata(ioe,ntf)
                    end do
                  end do

                end if

              endif
            enddo

  101       continue
!$OMP END DO
!$OMP END PARALLEL


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

         if( ( itall .eq. 2 .and. nobch .lt. maxbch .and. igsh .eq. 0 )
     &        .or. ( itall .eq. 4 .and. igsh .eq. 0 ) ) then

            call mk_2dnumfn(ctfln(m,iax),fname,itfll(m,iax),
     &           nobch,maxbch,npe)

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

               call sedech(iot,m,iax,1)

*-----------------------------------------------------------------------
* manatally
*  0  :  user defined analysis
*  1  :  systematic uncertainty
*  2  :  c-value dependence
*-----------------------------------------------------------------------

*-----------------------------------------------------------------------
*        SED, reg axes
*-----------------------------------------------------------------------

         if( itaxs(m,iax) .eq. 1
     &    .or. itaxs(m,iax) .eq. 2 ) then

          select case( itaxs(m,iax) )
           case ( 1 )
            iDaxis = 1   ! sed axis
           case ( 2 )
            iDaxis = 2   ! reg axis
          end select
          if( manatally.eq.2 ) iDaxis = 51 ! S.H. 2021.8.15

C rearrange data from anatalrst to anataldata
          call anatal_rearrange_sum(0,0,1,
     &         iDaxis,ntaxis,anatalrst,
     &         np,     ne,      0,    0,     nr,  0,  0,   0, nrst,
     &         eb,ew,  1,1,  1,1, dlr,  dt_one,  dt_one,1,
     &         vl,
     &         nij,ibin,nijaxs,nijaxs2,nijaxs3,nfgmax,
     &         cij,cijaxs,cijaxs2,nanataldata,ierr)

          if ( itaxs(m,iax) .eq. 1 ) cijaxs = hxunit(itsun(m))

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
               write(iot,'(/"x: ",a21)') cijaxs
            else
               write(iot,'(/"x: ",200a1)') (itaxt(m)(i:i),i=1,itaxl(m))
            end if
           end if

           if( lcy_txt(m) .gt. 0 ) then
              write(iot,'( "y: ",200a1)')(cy_txt(m)(i:i),i=1,lcy_txt(m))
           else
            if( itayl(m) .eq. 0 ) then
               write(iot,'( "y: ",a30)') hsunit(3*(itunt(m)-1)+itsun(m))
            else
               write(iot,'( "y: ",200a1)') (itayt(m)(i:i),i=1,itayl(m))
            end if
           end if

         if( itaxs(m,iax) .eq. 1 ) then ! sed axis
               if( itety(m) .eq. 3 .or. itety(m) .eq. 5 ) then
                  write(iot,'( "p: xlog ylog afac(0.8) form(0.9)")')
               else
                  write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')
               end if

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

           if ( manatally .eq. 0 ) then ! user defined analysis
            if ( iteps(m) .ne. 2 ) then
             if( itaxs(m,iax) .eq. 1 ) then ! sed axis
               write(iot,'( "h: n",12x,"x",12x,
     &                    1000(a10,"),hh0",a3," n  "))')
     &                     ( chp(i), chl(i) , i = 1, np_mxang ) ! frtati 2021/10/05 np -> np_mxang

             else
                write(iot,'( "h:   x      n",5x,"n",10x,
     &                     1000(a10,"),",a4," n   "))')
     &                     ( chp(i), chm(i) , i = 1, np_mxang ) ! frtati 2021/10/05 np -> np_mxang
             end if

            else
             if( itaxs(m,iax) .eq. 1 ) then ! sed axis
               write(iot,'( "h: n",12x,"x",12x,
     &         1000(a1,i0,a9,"),hh0",a3
     &         ," ny",i0," dy",i0,"=[y",i0,"*y",i0,"] "))')
     &         ( chp(i)(1:1),i,chp(i)(2:10), chl(i)
     &         ,i+np_mxang,i,i,i+np_mxang, i = 1, np_mxang ) ! frtati 2021/10/05

             else
                write(iot,'( "h:   x      n",5x,"n",10x,
     &                     1000(a1,i0,a9,"),",a4
     &                    ," ny",i0," dy",i0,"=[y",i0,"*y",i0,"] "))')
     &                     ( chp(i)(1:1),i,chp(i)(2:10), chm(i)
     &                     ,i+np_mxang,i,i,i+np_mxang, i = 1, np_mxang ) ! frtati 2021/10/05
             end if
            end if

             if( itaxs(m,iax) .eq. 1 ) then ! sed axis
               write(iot,'( "#  lower        upper  ",3x,
     &                 1000(a1,2x,a8,4x,"r.err "))')
     &                 ( chb(i), chp(i)(3:10), i = 1, np ) ! frtati 2021/10/05

             else
               write(iot,'( "#  num    reg     volume  ",
     &                    1000(a1,2x,a8,4x,"r.err "))')
     &                    ( chb(i), chp(i)(3:10), i = 1, np ) ! frtati 2021/10/05
             end if

           else if ( manatally .eq. 1 ) then ! systematic uncertainty
            if ( iteps(m) .ne. 2 ) then
             if( itaxs(m,iax) .eq. 1 ) then ! sed axis
               write(iot,'( "h: n",12x,"x",12x,
     &                    1000(a10,"),hh0",a3," n n n "))')
     &                     ( chp(i), chl(i) , i = 1, np_mxang ) ! frtati 2021/10/05 np -> np_mxang

             else
                write(iot,'( "h:   x      n",5x,"n",10x,
     &                     1000(a10,"),",a4," n n n "))')
     &                     ( chp(i), chm(i) , i = 1, np_mxang ) ! frtati 2021/10/05 np -> np_mxang
             end if

            else
             if( itaxs(m,iax) .eq. 1 ) then ! sed axis
               write(iot,'( "h: n",12x,"x",12x,
     &         1000(a1,i0,a9,"),hh0",a3
     &         ," ny",i0," dy",i0,"=[y",i0,"*y",i0,"] n n "))')
     &         ( chp(i)(1:1),i,chp(i)(2:10), chl(i)
     &         ,i+np_mxang,i,i,i+np_mxang, i = 1, np_mxang ) ! frtati 2021/10/05

             else
                write(iot,'( "h:   x      n",5x,"n",10x,
     &                     1000(a1,i0,a9,"),",a4
     &                ," ny",i0," dy",i0,"=[y",i0,"*y",i0,"] n n "))')
     &                     ( chp(i)(1:1),i,chp(i)(2:10), chm(i)
     &                     ,i+np_mxang,i,i,i+np_mxang, i = 1, np_mxang ) ! frtati 2021/10/05
             end if
            end if

               write(iot,'( "#  num    reg     volume  ",
     &                    1000(a1,2x,a8,4x,"r.err(tot, syst, stat)"))')
     &                    ( chb(i), chp(i)(3:10), i = 1, np ) ! frtati 2021/10/05

           else if ( manatally .eq. 2 ) then ! c-value dependence
            if ( iteps(m) .ne. 2 ) then
                write(iot,'( "h:   x",10x,
     &                     1000(a10,"),",a4," n   "))')
     &                     ( chp(i), chm(i) , i = 1, np_mxang ) ! frtati 2021/10/05 np -> np_mxang

            else
                write(iot,'( "h:   x",10x,
     &                     1000(a1,i0,a9,"),",a4
     &                    ," ny",i0," dy",i0,"=[y",i0,"*y",i0,"] "))')
     &                     ( chp(i)(1:1),i,chp(i)(2:10), chm(i)
     &                     ,i+np_mxang,i,i,i+np_mxang, i = 1, np_mxang ) ! frtati 2021/10/05
            end if

               write(iot,'( "#  c-value   ",
     &                    1000(a1,2x,a8,4x,"r.err "))')
     &                    ( chb(i), chp(i)(3:10), i = 1, np ) ! frtati 2021/10/05
           end if

               do ip = 1, np
                     tott(ip,1) = 0.d+0
               end do

           voll = 0.0
           njaxs = nijaxs + 1

           do ijaxs = 1, nijaxs

      if( itaxs(m,iax) .eq. 1 ) then ! sed axis
! T.Sato 2017/11/21 calculate normalization factor
      if(itunt(m).ge.13) then ! yf(y) or yd(y) mode
       do ip=1,np
        totnorm(ip)=0.0
        do ie=1,ne
         if(eb(ie).eq.0.0) then ! lower boundary is 0
          totnorm(ip)=totnorm(ip)+anataldata(ip,ianataldata,ijaxs,1)
     &    *(fgaxs(ijaxs+1)-fgaxs(ijaxs))/(fgaxs(ijaxs+1)/2.0d0)
         else ! normal case
          totnorm(ip)=totnorm(ip)+anataldata(ip,ianataldata,ijaxs,1)
     &         *(fgaxs(ijaxs+1)-fgaxs(ijaxs))
     &     /sqrt(fgaxs(ijaxs+1)*fgaxs(ijaxs))
         endif
        enddo
        if(totnorm(ip).eq.0.0) totnorm(ip)=1.0d0
        anataldata(ip,ianataldata,ijaxs,1)
     &       = anataldata(ip,ianataldata,ijaxs,1) / totnorm(ip)
       enddo
      else ! conventional mode
       totnorm(:)=1.0d0
      endif
      end if

            if ( manatally .eq. 0 ) then ! user defined analysis
             if( itaxs(m,iax) .eq. 1 ) then ! sed axis
               write(iot,'(1p2e13.4,1000(1pe13.4,0pf8.4))')
     &              fgaxs(ijaxs),fgaxs(ijaxs+1),
     &              ((anataldata(ip,ianataldata,ijaxs,irst),irst=1,2)
     &              ,ip=1,np) ! frtati 2021/10/05

             else
                write(iot,'(i5,1x,i7,1pe13.4,1000(1pe13.4,0p3f8.4))')
     &               ijaxs, lr(ijaxs), fgaxs(ijaxs),
     &               ((anataldata(ip,ianataldata,ijaxs,irst),irst=1,2)
     &               ,ip=1,np) ! frtati 2021/10/05
             end if

            else if ( manatally .eq. 1 ) then ! systematic uncertainty
             if( itaxs(m,iax) .eq. 1 ) then ! sed axis
               write(iot,'(1p2e13.4,1000(1pe13.4,0p3f8.4))')
     &              fgaxs(ijaxs),fgaxs(ijaxs+1),
     &              ((anataldata(ip,ianataldata,ijaxs,irst),irst=1,4)
     &              ,ip=1,np) ! frtati 2021/10/05

             else
                write(iot,'(i5,1x,i7,1pe13.4,1000(1pe13.4,0p3f8.4))')
     &               ijaxs, lr(ijaxs), fgaxs(ijaxs),
     &               ((anataldata(ip,ianataldata,ijaxs,irst),irst=1,4)
     &               ,ip=1,np) ! frtati 2021/10/05
             end if

            else if ( manatally .eq. 2 ) then ! c-value dependence
                write(iot,'(1pe13.4,1000(1pe13.4,0pf8.4))')
     &               fgaxs(ijaxs),
     &               ((anataldata(ip,ianataldata,ijaxs,irst),irst=1,2)
     &               ,ip=1,np) ! frtati 2021/10/05

            end if

            voll = voll + delvol(ianataldata,ijaxs)

                  do ip = 1, np
               vn = anataldata(ip,ianataldata,ijaxs,1)
     &              * delvol(ianataldata,ijaxs)
                     tott(ip,1) = tott(ip,1) + vn
                  end do

               end do

               do ip = 1, np
             if( itaxs(m,iax) .eq. 1 ) then ! sed axis
                     if( itunt(m) .ne.  5 .and. itunt(m) .ne. 6 ) then
                        tott(ip,1) = tott(ip,1) / voll
                     end if

             else
                     tott(ip,1) = tott(ip,1)
             end if
               end do

! sumover
            do ip = 1, np
            do irst=1,2
              tott(ip,irst) = anataldata(ip,ianataldata,njaxs,irst)
            enddo
            enddo

             if( itaxs(m,iax) .eq. 1 ) then ! sed axis
               write(iot,'(/"#   sum over in cm or MeV ",
     &         1000(1pe13.4,0pf8.4))') (tott(ip,1),tott(ip,2),ip=1,np)   ! T.Sato 2017/11/21

             else
               write(iot,'(/"#   sum over ",   13x ,
     6               1000(1pe13.4,0pf8.4,16x))')
     &              (tott(ip,1),tott(ip,2),ip=1,np) ! frtati 2021/10/05
             end if

*-----------------------------------------------------------------------

           write(iot,'(/a1,"no. =",i3,a1)') cha, inum, cha

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

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
     &                     /1pe13.4,2x,"$--$",1pe13.4)')
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

         end if

            call prestart(m,iot) !OBINATA(2012.9.7)

            close(iot)

         if( iteps(m) .ne. 0 ) then

            open(iot, file = fname, status = 'unknown' )
            call a_angel(idasa,fname)

         end if

*-----------------------------------------------------------------------

  900 continue

*-----------------------------------------------------------------------
      deallocate( ixyz )
      deallocate (anatalrst)
      deallocate (ew,vl,val,rdata,dlr)
      deallocate (lr)

      return
      end subroutine anatal_psedreg



************************************************************************
*                                                                      *
      subroutine anatal_psedrz(m,np,nr,nz,ne,rm,zm,eb,
     &                  tr,nfile,weightRate,idasa,manatally)
*                                                                      *
*       output r-z scoring mesh SED tally                              *
*       created by T.Miura on 2021/09/13                               *
*                                                                      *
************************************************************************
      use sangelmod, only: itsans
      use anatallymod, only: anataldata,fg,fgaxs,fgaxs2,fgaxs3,delvol,
     &     cx_txt,lcx_txt,cy_txt,lcy_txt,cz_txt,lcz_txt
      use partmod ! frtati 2021/10/05

*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

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

      common /tall27/ itaxl(itlmax), itaxt(itlmax),
     &                itayl(itlmax), itayt(itlmax),
     &                itazl(itlmax), itazt(itlmax)
      character itaxt*200, itayt*200, itazt*200

      common /tall21/ rtfac(itlmax)
      common /tall35/ iteps(itlmax)

      common /tall57/ itsun(itlmax),rtdim(itlmax),rtucv(itlmax),
     &rtrho(itlmax),itmodel(itlmax) ! T.Sato 2022/08/14

      common /istcut/ ist_cut, ist_bat
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
      dimension   eb(ne+1)
!      dimension   ew(ne)
      dimension   tr(np,ne,nr,nz,2*nfile)
      integer irst,nrst, itmpdata
!      dimension   anatalrst(np,ne+1,1,1,nr*nz*1,1,2*nfile+3)
! sumover
!      dimension   anatalrst(np,ne+1,1,1,(nr+1)*(nz+1)*1,1,2*nfile+3)
      real(8),allocatable :: anatalrst(:,:,:,:,:,:,:)

      integer ianataldata,nanataldata
      integer, parameter :: ntaxis=9  ! S.H. 2021.8.9
      integer ij(ntaxis),nij(ntaxis),ibin(ntaxis)
      integer nijaxs,nijaxs2,nijaxs3,nfgmax
      character cij(ntaxis)*15,cijaxs*15,cijaxs2*15
      integer iDaxis

      integer     nfile, manatally
!      dimension   weightRate(nfile), rdata(2,nfile)
      dimension   weightRate(nfile)

*-----------------------------------------------------------------------

!      dimension vl(nr,nz,1) ! T.Sato 2017/11/21
      real(8),allocatable :: ew(:),rdata(:,:),vl(:,:,:)


      real(8),allocatable :: vl_r(:,:),vl_z(:,:)
      dimension dt_one(1)
      data dt_one /1.0d0/

      dimension tott(6,2),totnorm(6) ! totnorm is added on 2017/11/21 for normalization

      character hsunit(24)*30

      data hsunit(1) / 'T.L. [cm/MeV/source]  '/
      data hsunit(2) / 'T.L. [cm/(keV/um)/source]'/
      data hsunit(3) / 'T.L. [cm/Gy/source]  '/
      data hsunit(4) / 'Dose [MeV/MeV/source]  '/
      data hsunit(5) / 'Dose [MeV/(keV/um)/source]'/
      data hsunit(6) / 'Dose [MeV/Gy/source]  '/
      data hsunit(7) / 'T.L. [cm/ln(MeV)/source]   '/
      data hsunit(8) / 'T.L. [cm/ln(keV/um)/source]'/
      data hsunit(9) / 'T.L. [cm/ln(Gy)/source]    '/
      data hsunit(10)/ 'Dose [MeV/ln(MeV)/source]   '/
      data hsunit(11)/ 'Dose [MeV/ln(keV/um)/source]'/
      data hsunit(12)/ 'Dose [MeV/ln(Gy)/source]    '/
      data hsunit(13)/ 'T.L. [cm/source]'/
      data hsunit(14)/ 'T.L. [cm/source]'/
      data hsunit(15)/ 'T.L. [cm/source]'/
      data hsunit(16)/ 'Dose [MeV/source]'/
      data hsunit(17)/ 'Dose [MeV/source]'/
      data hsunit(18)/ 'Dose [MeV/source]'/
      data hsunit(19)/ 'E*f(E) [dimensionless]'/
      data hsunit(20)/ 'y*f(y) [dimensionless]'/
      data hsunit(21)/ 'z*f(z) [dimensionless]'/
      data hsunit(22)/ 'E*d(E) [MeV]'/
      data hsunit(23)/ 'y*d(y) [keV/um]'/
      data hsunit(24)/ 'z*d(z) [Gy]'/

      character hxunit(3)*12
      data hxunit(1) / 'z [MeV]     '/
      data hxunit(2) / 'y [keV/um]'/
      data hxunit(3) / 'z [Gy]      '/

      character dc2*4
      character chp(itmxpt)*11  ! kitamura22/03/31
      character chq(itmxpt)*9   ! kitamura22/03/31


*-----------------------------------------------------------------------

      data pi/3.14159265d+0/

      character rpa*1
      data rpa /'}'/
      character yen*1

*-----------------------------------------------------------------------
      integer           iat, iad, iaf
      iat(iad,iaf) = iad + (iaf-1) * 2

*-----------------------------------------------------------------------
      integer           icf,ir,iz
      icf(ir,iz) = ir + ( iz - 1 ) * (nr+1)

*-----------------------------------------------------------------------
*        set mesh volume ( unit = 1-3 ; vol = 1.0 )
*-----------------------------------------------------------------------
      allocate (ew(ne),rdata(2,nfile),vl(nr,nz,1))
! sumover
      allocate (vl_r(nz,1),vl_z(nr,1))
      vl_r(:,:) = 1.0d0
      vl_z(:,:) = 1.0d0

      do ir=1,nr  ! T.Sato 2017/11/21
       do iz=1,nz
                 vl(ir,iz,1) = 1.0d0
       enddo
      enddo

*-----------------------------------------------------------------------

      yen  = char(92)
      igsh = 0

      nrst = 2*nfile+3 ! S.H. 2021.8.15
      allocate (anatalrst(np,ne+1,1,1,(nr+1)*(nz+1)*1,1,nrst))


*-----------------------------------------------------------------------
      np_mxang = np
      if ( itmxang(m).ne.0 .and. iabs(itmxang(m)).lt.np ) then
        np_mxang = iabs(itmxang(m))
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
*           itunt(m) = 1 : cm/z/source
*                    = 2 : MeV/z/source
*                    = 3 : cm/ln(z)/source
*                    = 4 : MeV/ln(z)/source
*                    = 5 : cm/source
*                    = 6 : MeV/source
*                    = 7 : y*f(y)
*                    = 8 : y*d(y)
*           itsun(m) = 1 : MeV
*                    = 2 : keV/um
*                    = 3 : Gy
*-----------------------------------------------------------------------
            if( itunt(m) .le. 2 .or. itunt(m) .eq. 7) then ! T.Sato 2017/11/21

               do i = 1, ne

                  ew(i) = ( eb(i+1) - eb(i) )

               end do
               ew_sum = eb(ne+1) - eb(1)

            else if( itunt(m) .le. 4 .or. itunt(m) .eq. 8) then ! T.Sato 2017/11/21

               do i = 1, ne

                  ew(i) = log( eb(i+1) / eb(i) )

               end do
               ew_sum = log( eb(ne+1) / eb(1) )

            else

               do i = 1, ne

                  ew(i) = 1.d0

               end do
               ew_sum = 1.0d0

            end if

*-----------------------------------------------------------------------
*        relative error and  unit conversion
*-----------------------------------------------------------------------

                  cmax = 0.0
                  cmin = 1.e+33
                  isdz = 0
                  stdm = 0.0
                  itstd(m) = 1

            do 100 iz = 1, nz
            do 100 ir = 1, nr
            do 100 ie = 1, ne
            do 100 ip = 1, np
            do ntf = 1, nfile

               if( tr(ip,ie,ir,iz,iat(1,ntf)) .gt. 0.d0 ) then

                  call calc_stdev(m,Xa,sigx,
     &                            tr(ip,ie,ir,iz,iat(1,ntf)),
     &                            tr(ip,ie,ir,iz,iat(2,ntf)),
     &                            rtfac(m)/vl(ir,iz,1)/ew(ie))

                  tr(ip,ie,ir,iz,iat(1,ntf)) = Xa
                  tr(ip,ie,ir,iz,iat(2,ntf)) = sigx
                  if( sigx .gt. stdm ) stdm = sigx

                  if( tr(ip,ie,ir,iz,iat(1,ntf)) .gt. cmax )
     &                 cmax = tr(ip,ie,ir,iz,iat(1,ntf))

                  if( tr(ip,ie,ir,iz,iat(1,ntf)) .lt. cmin )
     &                 cmin = tr(ip,ie,ir,iz,iat(1,ntf))

               else

                  isdz = 1
                  tr(ip,ie,ir,iz,iat(2,ntf)) = 0.0

               end if

! sumover
               call psedrz_sumover_stdev_ntf(0,m,ntf,
     &               ip,ie,ir,iz,
     &               rtfac(m),ew(ie),vl(ir,iz,1),
     &               ew_sum,vl_r(iz,1),vl_z(ir,1))

            end do
  100       continue

!$OMP PARALLEL
!$OMP& private(ipomp,ir,iz,ie,ip,ntf,itmprz,rdata,answer,rerr)
!$OMP& private(iax,ir_a,iz_a,ie_a,ido_ana,ioe)
!$      npomp=OMP_GET_NUM_THREADS()
!$      ipomp=OMP_GET_THREAD_NUM()
!$      write(*,'(''OpenMP PARALLEL PROCESS'',
!$   &     i4,''/'',i4)') ipOMP+1,npomp
!$OMP DO SCHEDULE(dynamic) reduction(max:cmax) reduction(min:cmin)
        do itmprz=1,nr*nz
            ir=(itmprz-1)/nz+1
            iz=itmprz-(ir-1)*nz
            do 101 ie = 1, ne
            do 101 ip = 1, np

            if( manatally .eq. 0 ) then ! T.Sato 2020/10/22

               do ntf = 1, nfile
                  rdata(1,ntf)=tr(ip,ie,ir,iz,iat(1,ntf))
                  rdata(2,ntf)=tr(ip,ie,ir,iz,iat(2,ntf))
               end do

              call usranatal(nfile,rdata,answer,rerr)

              anatalrst(ip,ie,1,1,icf(ir,iz),1,1) = answer
              anatalrst(ip,ie,1,1,icf(ir,iz),1,2) = rerr

              if(answer.gt.cmax) cmax = answer
              if(answer.gt.0d0 .and. answer.lt.cmin) cmin = answer

            else if( manatally .eq. 1 ) then

               do ntf = 1, nfile
                  rdata(1,ntf)=tr(ip,ie,ir,iz,iat(1,ntf))
                  rdata(2,ntf)=tr(ip,ie,ir,iz,iat(2,ntf))
               end do

               call anova(nfile,rdata,resc3(m)
     &              ,fmval,unca,uncerr,uratio)

               anatalrst(ip,ie,1,1,icf(ir,iz),1,1) ! mean
     &              = fmval
               if ( fmval .gt. 0d0 ) then
                  anatalrst(ip,ie,1,1,icf(ir,iz),1,2) ! total uncertainty
     &                 = dsqrt(unca**2 + uncerr**2/resc3(m)) / fmval
                  anatalrst(ip,ie,1,1,icf(ir,iz),1,3) ! syst. uncertainty
     &                 = unca / fmval
                  anatalrst(ip,ie,1,1,icf(ir,iz),1,4) ! stat. uncertainty
     &                 = uncerr /dsqrt(resc3(m)) / fmval
                  anatalrst(ip,ie,1,1,icf(ir,iz),1,5) ! squared unc. ratio
     &                 = uratio
               else
                  anatalrst(ip,ie,1,1,icf(ir,iz),1,2) ! total uncertainty
     &                 = 0d0
                  anatalrst(ip,ie,1,1,icf(ir,iz),1,3) ! syst. uncertainty
     &                 = 0d0
                  anatalrst(ip,ie,1,1,icf(ir,iz),1,4) ! stat. uncertainty
     &                 = 0d0
                  anatalrst(ip,ie,1,1,icf(ir,iz),1,5) ! squared unc. ratio
     &                 = 0d0
               end if

            else if( manatally .eq. 2 ) then

               do ntf = 1, nfile
                do ioe = 1, 2
                  anatalrst(ip,ie,1,1,icf(ir,iz),1,iat(ioe,ntf))
     &              = tr(ip,ie,ir,iz,iat(ioe,ntf))
                end do
               end do

            end if

! sumover
            do iax = 1,itaxn(m)

              ie_a = ie
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
              endif

              if(ido_ana == 1) then

                 call get_sedrz_tr_sum_data(m,iax,nfile,
     &                 ip,ie,ir,iz,rdata)

                if( manatally .eq. 0 ) then ! T.Sato 2020/10/22

                   call usranatal(nfile,rdata,answer,rerr)

                   anatalrst(ip,ie_a,1,1,icf(ir_a,iz_a),1,1)
     &             = answer
                   anatalrst(ip,ie_a,1,1,icf(ir_a,iz_a),1,2)
     &             = rerr

                else if( manatally .eq. 1 ) then

                  call anova(nfile,rdata,resc3(m)
     &               ,fmval,unca,uncerr,uratio)

                  anatalrst(ip,ie_a,1,1,icf(ir_a,iz_a),1,1) ! mean
     &              = fmval
                 if ( fmval .gt. 0.0d0 ) then
                   anatalrst(ip,ie_a,1,1,icf(ir_a,iz_a),1,2) ! total uncertainty
     &                 = dsqrt(unca**2 + uncerr**2/resc3(m)) / fmval
                 else
                   anatalrst(ip,ie_a,1,1,icf(ir_a,iz_a),1,2) ! total uncertainty
     &                 = 0.0d0
                 end if

                else if( manatally .eq. 2 ) then
                  do ntf = 1, nfile
                    do ioe = 1, 2
                      anatalrst(ip,ie_a,1,1,1,icf(ir_a,iz_a),
     &                          iat(ioe,ntf)) = rdata(ioe,ntf)
                    end do
                  end do

                end if

              endif
            enddo

  101       continue
      enddo  ! rz loop T.Sato 2021/04/15
!$OMP END DO
!$OMP END PARALLEL


            if( nobch .gt. ist_bat ) then
            if( rtstd(m) .gt. 0.d0 .and.
     &          stdm .gt. 0.d0 .and. isdz .eq. 0 .and.
     &          stdm .lt. rtstd(m) ) itstd(m) = 0
            end if

*-----------------------------------------------------------------------
*     out put unit = 31  : temporary number
*-----------------------------------------------------------------------

      do iax = 1, itfln(m)

!OBINATA(2012.9.7): output *.err
        noe = 1
        if ( any( itaxs(m,iax) .eq. (/ 7, 8, 9, 10 /) )
     &      .and. ittwo(m) .ne. 4 .and. icntl .ne. 8 )
     &    noe = 2

        do ioe = 1, noe

         if( itall .eq. 2 .and. nobch .lt. maxbch ) then

            write(fnume,'(i3.3)') nobch
!OBINATA(2012.9.7): output *.err
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

!OBINATA(2012.9.7): output *.err
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

               call sedech(iot,m,iax,1)

*-----------------------------------------------------------------------
* manatally
*  0  :  user defined analysis
*  1  :  systematic uncertainty
*  2  :  c-value dependence
*-----------------------------------------------------------------------

*-----------------------------------------------------------------------
*        SED, r, z axes
*-----------------------------------------------------------------------

         if( itaxs(m,iax) .eq. 1 .or.
     &       itaxs(m,iax) .eq.  6 .or.
     &       itaxs(m,iax) .eq.  5 ) then

          select case( itaxs(m,iax) )
           case ( 1 )
            iDaxis = 1   ! sed axis
           case ( 6 )
            iDaxis = 6   ! r(-z) axis
           case ( 5 )
            iDaxis = 7   ! (r-)z axis
          end select
          if( manatally.eq.2 ) iDaxis = 51 ! S.H. 2021.8.15

C rearrange data from anatalrst to anataldata
          call anatal_rearrange_sum(0,0,2,
     &         iDaxis,ntaxis,anatalrst,
     &         np,     ne,      0,    0,     nr, nz,  0,    0, nrst,
     &         eb,ew,    1,1,  1,1,  rm, zm, dt_one,1,
     &         vl,
     &         nij,ibin,nijaxs,nijaxs2,nijaxs3,nfgmax,
     &         cij,cijaxs,cijaxs2,nanataldata,ierr)

          if ( itaxs(m,iax) .eq. 1 ) cijaxs = hxunit(itsun(m))

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
               write(iot,'( "y: ",a30)') hsunit(3*(itunt(m)-1)+itsun(m))
            else
               write(iot,'( "y: ",200a1)') (itayt(m)(i:i),i=1,itayl(m))
            end if
           end if

           if( itaxs(m,iax) .eq. 1 ) then ! sed axis
               if( itety(m) .eq. 3 .or. itety(m) .eq. 5 ) then
                  write(iot,'( "p: xlog ylog afac(0.8) form(0.9)")')
               else
                  write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')
               end if

           else if( itaxs(m,iax) .eq. 6 ) then ! r axis
               if( itrty(m) .eq. 3 .or. itrty(m) .eq. 5 ) then
                  write(iot,'( "p: xlog ylog afac(0.8) form(0.9)")')
               else
                  write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')
               end if

           else if( itaxs(m,iax) .eq. 5 ) then ! z axis
               if( itzty(m) .eq. 3 .or. itzty(m) .eq. 5 ) then
                  write(iot,'( "p: xlog ylog afac(0.8) form(0.9)")')
               else
                  write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')
               end if

           end if

               if( itanl(m) .gt. 0 ) then

                  write(iot,'( "p: ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               end if

               if( itsans(m) .gt. 0 ) then

                  call write_sangel(iot,m,0)

               end if

           if ( manatally .eq. 0 ) then ! user defined analysis
            if ( iteps(m) .ne. 2 ) then
               write(iot,'( "h: n",12x,"x",12x,
     &                    1000(a10,"),hh0",a3," n "))')
     &                    ( chp(i), chl(i) , i = 1, np_mxang ) ! frtati 2021/10/05 np -> np_mxang

            else
               write(iot,'( "h: n",12x,"x",12x,
     &                    1000(a1,i0,a9,"),hh0",a3
     &                    ," ny",i0," dy",i0,"=[y",i0,"*y",i0,"] "))')
     &                    ( chp(i)(1:1),i,chp(i)(2:10), chl(i)
     &                    ,i+np_mxang,i,i,i+np_mxang, i = 1, np_mxang ) ! frtati 2021/10/05
            end if

            write(iot,'( "#  lower        upper  ",3x,
     &                 1000(a1,2x,a8,4x,"r.err "))')
     &                 ( chb(i), chp(i)(3:10), i = 1, np ) ! frtati 2021/10/05

           else if ( manatally .eq. 1 ) then ! systematic uncertainty
            if ( iteps(m) .ne. 2 ) then
               write(iot,'( "h: n",12x,"x",12x,
     &                    1000(a10,"),hh0",a3," n n n "))')
     &                    ( chp(i), chl(i) , i = 1, np_mxang ) ! frtati 2021/10/05 np -> np_mxang

            else
               write(iot,'( "h: n",12x,"x",12x,
     &                    1000(a1,i0,a9,"),hh0",a3
     &                ," ny",i0," dy",i0,"=[y",i0,"*y",i0,"] n n "))')
     &                    ( chp(i)(1:1),i,chp(i)(2:10), chl(i)
     &                    ,i+np_mxang,i,i,i+np_mxang, i = 1, np_mxang ) ! frtati 2021/10/05
            end if

            write(iot,'( "#  lower        upper  ",3x,
     &                 1000(a1,2x,a8,4x,"r.err(tot, syst, stat)"))')
     &                 ( chb(i), chp(i)(3:10), i = 1, np ) ! frtati 2021/10/05

           else if ( manatally .eq. 2 ) then ! c-value dependence
            if ( iteps(m) .ne. 2 ) then
               write(iot,'( "h: n",12x,"x",12x,
     &                    1000(a10,"),hh0",a3," n "))')
     &                    ( chp(i), chl(i) , i = 1, np_mxang ) ! frtati 2021/10/05 np -> np_mxang

            else
               write(iot,'( "h: n",12x,"x",12x,
     &                    1000(a1,i0,a9,"),hh0",a3
     &                    ," ny",i0," dy",i0,"=[y",i0,"*y",i0,"] "))')
     &                    ( chp(i)(1:1),i,chp(i)(2:10), chl(i)
     &                    ,i+np_mxang,i,i,i+np_mxang, i = 1, np_mxang ) ! frtati 2021/10/05 np -> np_mxang

            end if

            write(iot,'( "#  lower        upper  ",3x,
     &                 1000(a1,2x,a8,4x,"r.err "))')
     &                 ( chb(i), chp(i)(3:10), i = 1, np ) ! frtati 2021/10/05

         end if

               do ip = 1, np
                     tott(ip,1) = 0.d+0
               end do

           voll = 0.0d0
           njaxs = nijaxs + 1

           do ijaxs = 1, nijaxs

      if( itaxs(m,iax) .eq. 1 ) then ! sed axis
! T.Sato 2017/11/21 calculate normalization factor
      if(itunt(m).ge.13) then ! yf(y) or yd(y) mode
       do ip=1,np
        totnorm(ip)=0.0
        do ie=1,ne
         if(eb(ie).eq.0.0) then ! lower boundary is 0
          totnorm(ip)=totnorm(ip)+anataldata(ip,ianataldata,ijaxs,1)
     &    *(fgaxs(ijaxs+1)-fgaxs(ijaxs))/(fgaxs(ijaxs+1)/2.0d0)
         else ! normal case
          totnorm(ip)=totnorm(ip)+anataldata(ip,ianataldata,ijaxs,1)
     &         *(fgaxs(ijaxs+1)-fgaxs(ijaxs))
     &     /sqrt(fgaxs(ijaxs+1)*fgaxs(ijaxs))
         endif
        enddo
        if(totnorm(ip).eq.0.0) totnorm(ip)=1.0d0
        anataldata(ip,ianataldata,ijaxs,1)
     &       = anataldata(ip,ianataldata,ijaxs,1) / totnorm(ip)
       enddo
      else ! conventional mode
       totnorm(:)=1.0d0
      endif
      endif

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

               end do

               do ip = 1, np
             if( itaxs(m,iax) .eq. 1 ) then ! sed axis
                     if( itunt(m) .ne.  5 .and. itunt(m) .ne. 6 ) then
                     tott(ip,1) = tott(ip,1) / voll
                     end if
             else
                     tott(ip,1) = tott(ip,1)
             end if
               end do

! sumover
            do ip = 1, np
            do irst=1,2
              tott(ip,irst) = anataldata(ip,ianataldata,njaxs,irst)
            enddo
            enddo

           if( itaxs(m,iax) .eq. 1 ) then ! sed axis
               write(iot,'(/"#   sum over in cm or MeV ",
     &         1000(1pe13.4,0pf8.4))') (tott(ip,1),tott(ip,2),ip=1,np)   ! T.Sato 2017/11/21

           else
               write(iot,'(/"#   sum over ",   13x ,
     6               1000(1pe13.4,0pf8.4,16x))')
     &              (tott(ip,1),tott(ip,2),ip=1,np) ! frtati 2021/10/05
           end if

*-----------------------------------------------------------------------

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

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
     &                     /1pe13.4,2x,"$--$",1pe13.4)')
     &                     cij(ijtmp),ij(ijtmp),fg(itmp),fg(itmp+1)
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

          select case( itaxs(m,iax) )
           case ( 10 )
            iDaxis = 34   ! rz
          end select

C rearrange data from anatalrst to anataldata
          call anatal_rearrange_sum(0,0,2,
     &         iDaxis,ntaxis,anatalrst,
     &         np,     ne,      0,    0,     nr, nz,  0,    0, nrst,
     &         eb,ew,    1,1,  1,1,  rm, zm, dt_one,1,
     &         vl,
     &         nij,ibin,nijaxs,nijaxs2,nijaxs3,nfgmax,
     &         cij,cijaxs,cijaxs2,nanataldata,ierr)

               inum = 0

            do ip = 1, np
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

*-----------------------------------------------------------------------

                  write(iot,'(/"#",78("-"))')

               if( inum .eq. 1 .or. ip.gt.np_mxang ) then ! frtati 2021/10/05
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

           if ( np.gt.np_mxang ) then
             write(iot,'( " SKIPPAGE:")')
             inum = inum - 1
           end if

*-----------------------------------------------------------------------

            if( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) then

              write(iot,'(/a1,"no. =",i3,a1)') cha, inum, cha

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

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
               write(iot,'(/"x: ",200a1)') (itaxt(m)(i:i),i=1,itaxl(m))
            end if
           end if

           if( lcy_txt(m) .gt. 0 ) then
              write(iot,'( "y: ",200a1)')(cy_txt(m)(i:i),i=1,lcy_txt(m))
           else
            if( itayl(m) .eq. 0 ) then
               write(iot,'("y: ",a15)') cijaxs2
            else
               write(iot,'( "y: ",200a1)') (itayt(m)(i:i),i=1,itayl(m))
            end if
           end if

*-----------------------------------------------------------------------

               if( inum .eq. 1 ) then

             form  = ( fgaxs2(nijaxs2+1) - fgaxs2(1) )
     &              / ( fgaxs(nijaxs+1) - fgaxs(1) )
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
     &            ( itang(m)(i:i),i = 1, itanl(m) )

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
     &                   " i2 = n2, 1, -1 )")')
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
     &                      "  SED        r.err")')

            do ijaxs = 1, nijaxs
             do ijaxs2 = 1, nijaxs2

                write(iot,'(1p3e11.3,0pf8.4)')
     &               fgaxs2(ijaxs2)  + axs2del/2d0,
     &               fgaxs(ijaxs)  + axs1del/2d0,
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

!OBINATA(2012.7.2): output *.err
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
     &"p: xfac[c2*0.05] form[c1/0.05] notc noxt noyt"/
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
           write(iot,'("y: ",a30)') hsunit(3*(itunt(m)-1)+itsun(m))
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
     &                 /1pe13.4,2x,"$--$",1pe13.4)')
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
          end do   ! do ip = 1, np

*-----------------------------------------------------------------------

          deallocate( fg )
          deallocate( fgaxs )
          deallocate( fgaxs2 )
          deallocate( fgaxs3 )
          deallocate( anataldata )
          deallocate( delvol )

*-----------------------------------------------------------------------

         end if

*-----------------------------------------------------------------------

            call prestart(m,iot) !OBINATA(2012.9.7)

            close(iot)

         if( iteps(m) .ne. 0 ) then

            open(iot, file = fname, status = 'unknown' )
            call a_angel(idasa,fname)

         end if

        end do

      end do

*-----------------------------------------------------------------------
      deallocate (vl_r,vl_z)
      deallocate (anatalrst)
      deallocate (ew,rdata,vl)

      return
      end subroutine anatal_psedrz



************************************************************************
*                                                                      *
      subroutine anatal_psedxyz(m,np,
     &                   nx,ny,nz,ne,xm,ym,zm,eb,tr,
     &                   nfile,weightRate,
     &                   igsh,idasa,manatally)
*                                                                      *
*       output xyz scoring mesh SED tally                              *
*       created by T.Miura on 2021/09/13                               *
*                                                                      *
************************************************************************
      use sangelmod, only: itsans
      use anatallymod, only: anataldata,fg,fgaxs,fgaxs2,fgaxs3,delvol,
     &     cx_txt,lcx_txt,cy_txt,lcy_txt,cz_txt,lcz_txt
      use partmod ! frtati 2021/10/05

*-----------------------------------------------------------------------

      implicit double precision (a-h,o-z)

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

      common /tall21/ rtfac(itlmax)
      common /tall28/ itgsh(itlmax)
      common /tall22/ itmcn(itlmax), itmtn(itlmax), itmtt(itlmax)
      common /tall34/ itvm(itlmax)
      common /tall35/ iteps(itlmax)
      common /tall38/ itres(itlmax)
      common /tall39/ rtwid(itlmax)

      common /tall46/ itmtr(itlmax,4), rtmtr(itlmax,13)
      common /tall49/ itglt(itlmax)

      common /tall57/ itsun(itlmax),rtdim(itlmax),rtucv(itlmax),
     &rtrho(itlmax),itmodel(itlmax) ! T.Sato 2022/08/14

      common /istcut/ ist_cut, ist_bat
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

      dimension   xm(nx+1)
      dimension   ym(ny+1)
      dimension   zm(nz+1)
      dimension   eb(ne+1)
!      dimension   ew(ne)
      dimension   tr(np,ne,nx*ny*nz,2*nfile)
      integer irst,nrst, itmpdata
!      dimension   anatalrst(np,ne+1,1,1,nx*ny*nz,1,2*nfile+3)
! sumover
!      dimension   anatalrst(np,ne+1,1,1,(nx+1)*(ny+1)*(nz+1),
!     &                      1,2*nfile+3)
      real(8),allocatable :: anatalrst(:,:,:,:,:,:,:)
      integer ianataldata,nanataldata
      integer, parameter :: ntaxis=9  ! S.H. 2021.8.9
      integer ij(ntaxis),nij(ntaxis),ibin(ntaxis)
      integer nijaxs,nijaxs2,nijaxs3,nfgmax
      character cij(ntaxis)*15,cijaxs*15,cijaxs2*15
      integer iDaxis

      integer     nfile, manatally
!      dimension   weightRate(nfile), rdata(2,nfile)
      dimension   weightRate(nfile)
      real(8),allocatable :: ew(:),rdata(:,:)

      integer,allocatable :: ixyz(:)
      real(8),allocatable :: vl(:,:,:)
      real(8),allocatable :: vl_x(:,:),vl_y(:,:),vl_z(:,:)

*-----------------------------------------------------------------------

      dimension tott(6,2),totnorm(6) ! totnorm is added on 2017/11/21 for normalization

      character hsunit(24)*30

      data hsunit(1) / 'T.L. [cm/MeV/source]  '/
      data hsunit(2) / 'T.L. [cm/(keV/um)/source]'/
      data hsunit(3) / 'T.L. [cm/Gy/source]  '/
      data hsunit(4) / 'Dose [MeV/MeV/source]  '/
      data hsunit(5) / 'Dose [MeV/(keV/um)/source]'/
      data hsunit(6) / 'Dose [MeV/Gy/source]  '/
      data hsunit(7) / 'T.L. [cm/ln(MeV)/source]   '/
      data hsunit(8) / 'T.L. [cm/ln(keV/um)/source]'/
      data hsunit(9) / 'T.L. [cm/ln(Gy)/source]    '/
      data hsunit(10)/ 'Dose [MeV/ln(MeV)/source]   '/
      data hsunit(11)/ 'Dose [MeV/ln(keV/um)/source]'/
      data hsunit(12)/ 'Dose [MeV/ln(Gy)/source]    '/
      data hsunit(13)/ 'T.L. [cm/source]'/
      data hsunit(14)/ 'T.L. [cm/source]'/
      data hsunit(15)/ 'T.L. [cm/source]'/
      data hsunit(16)/ 'Dose [MeV/source]'/
      data hsunit(17)/ 'Dose [MeV/source]'/
      data hsunit(18)/ 'Dose [MeV/source]'/
      data hsunit(19)/ 'E*f(E) [dimensionless]'/
      data hsunit(20)/ 'y*f(y) [dimensionless]'/
      data hsunit(21)/ 'z*f(z) [dimensionless]'/
      data hsunit(22)/ 'E*d(E) [MeV]'/
      data hsunit(23)/ 'y*d(y) [keV/um]'/
      data hsunit(24)/ 'z*d(z) [Gy]'/

      character hxunit(3)*12
      data hxunit(1) / 'z [MeV]     '/
      data hxunit(2) / 'y [keV/um]'/
      data hxunit(3) / 'z [Gy]      '/

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
      integer           iat, iad, iaf
      iat(iad,iaf) = iad + (iaf-1) * 2

*-----------------------------------------------------------------------

      icf(ix,iy,iz) = ix + ( iy - 1 ) * nx
     &              + ( iz - 1 ) * nx * ny
      icf2(ix,iy,iz) = ix + ( iy - 1 ) * (nx+1)
     &              + ( iz - 1 ) * (nx+1) * (ny+1)

*-----------------------------------------------------------------------
*        set mesh volume ( unit = 1-3 ; vol = 1.0 )
*-----------------------------------------------------------------------
      allocate (ew(ne),rdata(2,nfile))
! sumover
      allocate( vl_x(ny,nz),vl_y(nx,nz),vl_z(nx,ny) )
      vl_x(:,:) = 1.0d0
      vl_y(:,:) = 1.0d0
      vl_z(:,:) = 1.0d0

      allocate( vl(nx,ny,nz) )
      do ix=1,nx  ! T.Sato 2017/11/21
       do iy=1,ny
        do iz=1,nz
         vl(ix,iy,iz) = 1.0d0
        enddo
       enddo
      enddo

*-----------------------------------------------------------------------

      yen  = char(92)
      allocate( ixyz(max((nx+1)*(ny+1),(ny+1)*(nz+1),(nx+1)*(nz+1))) )

      nrst = 2*nfile+3 ! S.H. 2021.8.15
      allocate (anatalrst(np,ne+1,1,1,(nx+1)*(ny+1)*(nz+1),1,nrst))


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

         else

            npg = 1
            neg = 1

         end if

*-----------------------------------------------------------------------
*        particle name
*-----------------------------------------------------------------------

               call wtpname(m,np,chp,chq)

*-----------------------------------------------------------------------
*           itunt(m) = 1 : cm/z/source
*                    = 2 : MeV/z/source
*                    = 3 : cm/ln(z)/source
*                    = 4 : MeV/ln(z)/source
*                    = 5 : cm/source
*                    = 6 : MeV/source
*                    = 7 : y*f(y)
*                    = 8 : y*d(y)
*           itsun(m) = 1 : MeV
*                    = 2 : keV/um
*                    = 3 : Gy
*-----------------------------------------------------------------------
            if( itunt(m) .le. 2 .or. itunt(m) .eq. 7) then ! T.Sato 2017/11/21

               do i = 1, ne

                  ew(i) = ( eb(i+1) - eb(i) )

               end do
               ew_sum = eb(ne+1) - eb(1)

            else if( itunt(m) .le. 4 .or. itunt(m) .eq. 8) then ! T.Sato 2017/11/21

               do i = 1, ne

                  ew(i) = log( eb(i+1) / eb(i) )

               end do
               ew_sum = log( eb(ne+1) / eb(1) )

            else

               do i = 1, ne

                  ew(i) = 1.d0

               end do
              ew_sum = 1.0d0

            end if

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

           do itmpxyz=1,nx*ny*nz
            iz=(itmpxyz-1)/(ny*nx)+1
            iy=(itmpxyz-1-(iz-1)*ny*nx)/nx+1
            ix=itmpxyz-(iy-1)*nx-(iz-1)*ny*nx
            do 100 ie = 1, ne
            do 100 ip = 1, np
            do ntf = 1, nfile

               if( tr(ip,ie,icf(ix,iy,iz),iat(1,ntf)) .gt. 0.d0 ) then

                  call calc_stdev(m,Xa,sigx,
     &                            tr(ip,ie,icf(ix,iy,iz),iat(1,ntf)),
     &                            tr(ip,ie,icf(ix,iy,iz),iat(2,ntf)),
     &                            rtfac(m)/vl(ix,iy,iz)/ew(ie))

                  tr(ip,ie,icf(ix,iy,iz),iat(1,ntf)) = Xa
                  tr(ip,ie,icf(ix,iy,iz),iat(2,ntf)) = sigx
                  if( sigx .gt. stdm ) stdm = sigx

                  if( tr(ip,ie,icf(ix,iy,iz),iat(1,ntf)) .gt. cmax )
     &                     cmax = tr(ip,ie,icf(ix,iy,iz),iat(1,ntf))

                  if( tr(ip,ie,icf(ix,iy,iz),iat(1,ntf)) .lt. cmin )
     &                     cmin = tr(ip,ie,icf(ix,iy,iz),iat(1,ntf))

               else

                  isdz = 1
                  tr(ip,ie,icf(ix,iy,iz),iat(2,ntf)) = 0.d+0

               end if

! sumover
               call psedxyz_sumover_stdev_ntf(0,m,ntf,
     %               ip,ie,ix,iy,iz,
     &               rtfac(m),ew(ie),vl(ix,iy,iz),
     &               ew_sum,vl_x(iy,iz),vl_y(ix,iz),vl_z(ix,iy))

            end do
  100       continue
      enddo  ! xyz loop T.Sato 2021/04/15

!$OMP PARALLEL
!$OMP& private(ipomp,ix,iy,iz,ie,ip,ntf,itmpxyz,rdata,answer,rerr)
!$OMP& private(iax,ix_a,iy_a,iz_a,ie_a,ido_ana,ioe)
!$      npomp=OMP_GET_NUM_THREADS()
!$      ipomp=OMP_GET_THREAD_NUM()
!$      write(*,'(''OpenMP PARALLEL PROCESS'',
!$   &     i4,''/'',i4)') ipOMP+1,npomp
!$OMP DO SCHEDULE(dynamic) reduction(max:cmax) reduction(min:cmin)
           do itmpxyz=1,nx*ny*nz
            iz=(itmpxyz-1)/(ny*nx)+1
            iy=(itmpxyz-1-(iz-1)*ny*nx)/nx+1
            ix=itmpxyz-(iy-1)*nx-(iz-1)*ny*nx
            do 101 ie = 1, ne
            do 101 ip = 1, np

            if( manatally .eq. 0 ) then ! T.Sato 2020/10/22

               do ntf = 1, nfile
                  rdata(1,ntf)=tr(ip,ie,icf(ix,iy,iz),iat(1,ntf))
                  rdata(2,ntf)=tr(ip,ie,icf(ix,iy,iz),iat(2,ntf))
               end do

              call usranatal(nfile,rdata,answer,rerr)

              anatalrst(ip,ie,1,1,icf2(ix,iy,iz),1,1) = answer
              anatalrst(ip,ie,1,1,icf2(ix,iy,iz),1,2) = rerr

              if(answer.gt.cmax) cmax = answer
              if(answer.gt.0d0 .and. answer.lt.cmin) cmin = answer

            else if( manatally .eq. 1 ) then

               do ntf = 1, nfile
                  rdata(1,ntf)=tr(ip,ie,icf(ix,iy,iz),iat(1,ntf))
                  rdata(2,ntf)=tr(ip,ie,icf(ix,iy,iz),iat(2,ntf))
               end do

               call anova(nfile,rdata,resc3(m)
     &              ,fmval,unca,uncerr,uratio)

               anatalrst(ip,ie,1,1,icf2(ix,iy,iz),1,1) ! mean
     &              = fmval
               if ( fmval .gt. 0d0 ) then
                  anatalrst(ip,ie,1,1,icf2(ix,iy,iz),1,2) ! total uncertainty
     &                 = dsqrt(unca**2 + uncerr**2/resc3(m)) / fmval
                  anatalrst(ip,ie,1,1,icf2(ix,iy,iz),1,3) ! syst. uncertainty
     &                 = unca / fmval
                  anatalrst(ip,ie,1,1,icf2(ix,iy,iz),1,4) ! stat. uncertainty
     &                 = uncerr /dsqrt(resc3(m)) / fmval
                  anatalrst(ip,ie,1,1,icf2(ix,iy,iz),1,5) ! squared unc. ratio
     &                 = uratio
               else
                  anatalrst(ip,ie,1,1,icf2(ix,iy,iz),1,2) ! total uncertainty
     &                 = 0d0
                  anatalrst(ip,ie,1,1,icf2(ix,iy,iz),1,3) ! syst. uncertainty
     &                 = 0d0
                  anatalrst(ip,ie,1,1,icf2(ix,iy,iz),1,4) ! stat. uncertainty
     &                 = 0d0
                  anatalrst(ip,ie,1,1,icf2(ix,iy,iz),1,5) ! squared unc. ratio
     &                 = 0d0
               end if

            else if( manatally .eq. 2 ) then

               do ntf = 1, nfile
                do ioe = 1, 2
                  anatalrst(ip,ie,1,1,icf2(ix,iy,iz),1,iat(ioe,ntf))
     &              = tr(ip,ie,icf(ix,iy,iz),iat(ioe,ntf))
                end do
               end do

            end if

!sumover
            do iax = 1,itaxn(m)

              ie_a = ie
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
              endif

              if(ido_ana == 1) then

                 call get_sedxyz_tr_sum_data(m,iax,nfile,
     &                 ip,ie,ix,iy,iz,rdata)

                if( manatally .eq. 0 ) then ! T.Sato 2020/10/22

                   call usranatal(nfile,rdata,answer,rerr)

                   anatalrst(ip,ie_a,1,1,icf2(ix_a,iy_a,iz_a),1,1)
     &             = answer
                   anatalrst(ip,ie_a,1,1,icf2(ix_a,iy_a,iz_a),1,2)
     &             = rerr

                else if( manatally .eq. 1 ) then

                  call anova(nfile,rdata,resc3(m)
     &               ,fmval,unca,uncerr,uratio)

                  anatalrst(ip,ie_a,1,1,icf2(ix_a,iy_a,iz_a),1,1) ! mean
     &              = fmval
                 if ( fmval .gt. 0.0d0 ) then
                   anatalrst(ip,ie_a,1,1,icf2(ix_a,iy_a,iz_a),1,2) ! total uncertainty
     &                 = dsqrt(unca**2 + uncerr**2/resc3(m)) / fmval
                 else
                   anatalrst(ip,ie_a,1,1,icf2(iz_a,iy_a,iz_a),1,2) ! total uncertainty
     &                 = 0.0d0
                 end if

                else if( manatally .eq. 2 ) then
                  do ntf = 1, nfile
                    do ioe = 1, 2
                      anatalrst(ip,ie_a,1,1,icf2(ix_a,iy_a,iz_a),
     &                          1,iat(ioe,ntf)) = rdata(ioe,ntf)
                    end do
                  end do

                end if

              endif
            enddo
  101       continue
      enddo  ! xyz loop T.Sato 2021/04/15
!$OMP END DO
!$OMP END PARALLEL


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

!OBINATA(2012.9.7): output *.err
        noe = 1
        if ( any( itaxs(m,iax) .eq. (/ 7, 8, 9, 10 /) )
     &      .and. ittwo(m) .ne. 4 .and. icntl .ne. 8 )
     &    noe = 2

        do 900 ioe = 1, noe

         if( igsh .ne. 0 .and.
     &     ( itaxs(m,iax) .lt. 7 .or.
     &       ittwo(m) .eq. 4 .or. ittwo(m) .eq. 5 ) ) goto 900

         if( itall .eq. 2 .and. nobch .lt. maxbch .and.
     &       igsh .eq. 0 ) then

            write(fnume,'(i3.3)') nobch
!OBINATA(2012.9.7): output *.err
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

!OBINATA(2012.9.7): output *.err
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

               call sedech(iot,m,iax,1)

*-----------------------------------------------------------------------
* manatally
*  0  :  user defined analysis
*  1  :  systematic uncertainty
*  2  :  c-value dependence
*-----------------------------------------------------------------------

*-----------------------------------------------------------------------
*        SED, x, y, z axes
*-----------------------------------------------------------------------

         if( itaxs(m,iax) .eq.  1 .or.
     &       itaxs(m,iax) .eq.  3 .or.
     &       itaxs(m,iax) .eq.  4 .or.
     &       itaxs(m,iax) .eq.  5 ) then

          select case( itaxs(m,iax) )
           case ( 1 )
            iDaxis = 1   ! sed axis
           case ( 3 )
            iDaxis = 3   ! x axis
           case ( 4 )
            iDaxis = 4   ! y axis
           case ( 5 )
            iDaxis = 5   ! z axis
          end select
          if( manatally.eq.2 ) iDaxis = 51 ! S.H. 2021.8.15

C rearrange data from anatalrst to anataldata
          call anatal_rearrange_sum(0,0,3,
     &         iDaxis,ntaxis,anatalrst,
     &         np,     ne,      0,    0,     nx, ny, nz,    0, nrst,
     &         eb,ew,    1,1,  1,1,  xm, ym, zm,1,
     &         vl,
     &         nij,ibin,nijaxs,nijaxs2,nijaxs3,nfgmax,
     &         cij,cijaxs,cijaxs2,nanataldata,ierr)

          if ( itaxs(m,iax) .eq. 1 ) cijaxs = hxunit(itsun(m))

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
               write(iot,'( "y: ",a30)') hsunit(3*(itunt(m)-1)+itsun(m))
            else
               write(iot,'( "y: ",200a1)') (itayt(m)(i:i),i=1,itayl(m))
            end if
           end if

           if( itaxs(m,iax) .eq. 1 ) then ! sed axis
               if( itety(m) .eq. 3 .or. itety(m) .eq. 5 ) then
                  write(iot,'( "p: xlog ylog afac(0.8) form(0.9)")')
               else
                  write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')
               end if

           else if( itaxs(m,iax) .eq. 3 ) then ! x axis
               if( itxty(m) .eq. 3 .or. itxty(m) .eq. 5 ) then
                  write(iot,'( "p: xlog ylog afac(0.8) form(0.9)")')
               else
                  write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')
               end if

           else if( itaxs(m,iax) .eq. 4 ) then ! y axis
               if( ityty(m) .eq. 3 .or. ityty(m) .eq. 5 ) then
                  write(iot,'( "p: xlog ylog afac(0.8) form(0.9)")')
               else
                  write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')
               end if

           else if( itaxs(m,iax) .eq. 5 ) then ! z axis
               if( itzty(m) .eq. 3 .or. itzty(m) .eq. 5 ) then
                  write(iot,'( "p: xlog ylog afac(0.8) form(0.9)")')
               else
                  write(iot,'( "p: xlin ylog afac(0.8) form(0.9)")')
               end if

           end if

               if( itanl(m) .gt. 0 ) then

                  write(iot,'( "p: ",200a1)')
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               end if

               if( itsans(m) .gt. 0 ) then

                  call write_sangel(iot,m,0)

               end if

           if ( manatally .eq. 0 ) then ! user defined analysis
            if ( iteps(m) .ne. 2 ) then
               write(iot,'( "h: n",12x,"x",12x,
     &                    1000(a10,"),hh0",a3," n "))')
     &                    ( chp(i), chl(i) , i = 1, np_mxang ) ! frtati 2021/10/05 np -> np_mxang

            else
               write(iot,'( "h: n",12x,"x",12x,
     &                    1000(a1,i0,a9,"),hh0",a3
     &                    ," ny",i0," dy",i0,"=[y",i0,"*y",i0,"] "))')
     &                    ( chp(i)(1:1),i,chp(i)(2:10), chl(i)
     &                    ,i+np_mxang,i,i,i+np_mxang, i = 1, np_mxang ) ! frtati 2021/10/05
            end if

            write(iot,'( "#  lower        upper  ",3x,
     &                 1000(a1,2x,a8,4x,"r.err "))')
     &                 ( chb(i), chp(i)(3:10), i = 1, np ) ! frtati 2021/10/05

           else if ( manatally .eq. 1 ) then ! systematic uncertainty
            if ( iteps(m) .ne. 2 ) then
               write(iot,'( "h: n",12x,"x",12x,
     &                    1000(a10,"),hh0",a3," n n n "))')
     &                    ( chp(i), chl(i) , i = 1, np_mxang ) ! frtati 2021/10/05 np -> np_mxang

            else
               write(iot,'( "h: n",12x,"x",12x,
     &                    1000(a1,i0,a9,"),hh0",a3
     &                ," ny",i0," dy",i0,"=[y",i0,"*y",i0,"] n n "))')
     &                    ( chp(i)(1:1),i,chp(i)(2:10), chl(i)
     &                    ,i+np_mxang,i,i,i+np_mxang, i = 1, np_mxang ) ! frtati 2021/10/05
            end if

            write(iot,'( "#  lower        upper  ",3x,
     &                 1000(a1,2x,a8,4x,"r.err(tot, syst, stat)"))')
     &                 ( chb(i), chp(i)(3:10), i = 1, np ) ! frtati 2021/10/05

           else if ( manatally .eq. 2 ) then ! c-value dependence
            if ( iteps(m) .ne. 2 ) then
               write(iot,'( "h: n",12x,"x",12x,
     &                    1000(a10,"),hh0",a3," n "))')
     &                    ( chp(i), chl(i) , i = 1, np_mxang ) ! frtati 2021/10/05 np -> np_mxang

            else
               write(iot,'( "h: n",12x,"x",12x,
     &                    1000(a1,i0,a9,"),hh0",a3
     &                    ," ny",i0," dy",i0,"=[y",i0,"*y",i0,"] "))')
     &                    ( chp(i)(1:1),i,chp(i)(2:10), chl(i)
     &                    ,i+np_mxang,i,i,i+np_mxang, i = 1, np_mxang ) ! frtati 2021/10/05 np -> np_mxang

            end if

            write(iot,'( "#  lower        upper  ",3x,
     &                 1000(a1,2x,a8,4x,"r.err "))')
     &                 ( chb(i), chp(i)(3:10), i = 1, np ) ! frtati 2021/10/05

         end if

               do ip = 1, np
                     tott(ip,1) = 0.d+0
               end do

           voll = 0.0d0
           njaxs = nijaxs + 1


           do ijaxs = 1, nijaxs

      if( itaxs(m,iax) .eq. 1 ) then ! sed axis
! T.Sato 2017/11/21 calculate normalization factor
      if(itunt(m).ge.13) then ! yf(y) or yd(y) mode
       do ip=1,np
        totnorm(ip)=0.0
        do ie=1,ne
         if(eb(ie).eq.0.0) then ! lower boundary is 0
          totnorm(ip)=totnorm(ip)+anataldata(ip,ianataldata,ijaxs,1)
     &    *(fgaxs(ijaxs+1)-fgaxs(ijaxs))/(fgaxs(ijaxs+1)/2.0d0)
         else ! normal case
          totnorm(ip)=totnorm(ip)+anataldata(ip,ianataldata,ijaxs,1)
     &         *(fgaxs(ijaxs+1)-fgaxs(ijaxs))
     &     /sqrt(fgaxs(ijaxs+1)*fgaxs(ijaxs))
         endif
        enddo
        if(totnorm(ip).eq.0.0) totnorm(ip)=1.0d0
        anataldata(ip,ianataldata,ijaxs,1)
     &       = anataldata(ip,ianataldata,ijaxs,1) / totnorm(ip)
       enddo
      else ! conventional mode
       totnorm(:)=1.0d0
      endif
      end if

            if ( manatally .eq. 0 ) then ! user defined analysis
               write(iot,'(1p2e13.4,1000(1pe13.4,0p3f8.4))')
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

               end do

               do ip = 1, np
             if( itaxs(m,iax) .eq. 1 ) then ! let axis
                     if( itunt(m) .ne.  5 .and. itunt(m) .ne. 6 ) then
                     tott(ip,1) = tott(ip,1) / voll
                     end if

             else
                        tott(ip,1) = tott(ip,1)
             end if
               end do

! sumover
            do ip = 1, np
            do irst=1,2
              tott(ip,irst) = anataldata(ip,ianataldata,njaxs,irst)
            enddo
            enddo

             if( itaxs(m,iax) .eq. 1 ) then ! sed axis
               write(iot,'(/"#   sum over in cm or MeV ",
     &         1000(1pe13.4,0pf8.4))') (tott(ip,1),tott(ip,2),ip=1,np)   ! T.Sato 2017/11/21

             else
               write(iot,'(/"#   sum over ",   13x ,
     &               1000(1pe13.4,0pf8.4,16x))')
     &              (tott(ip,1),tott(ip,2),ip=1,np) ! frtati 2021/10/05
             end if

*-----------------------------------------------------------------------

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

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
     &                     /1pe13.4,2x,"$--$",1pe13.4)')
     &                     cij(ijtmp),ij(ijtmp),fg(itmp),fg(itmp+1)
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
*        xy, yz, xz axes (matrix)
*-----------------------------------------------------------------------

         else if( itaxs(m,iax) .eq. 7
     &         .or. itaxs(m,iax) .eq. 8
     &         .or. itaxs(m,iax) .eq. 9 ) then

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
            iDaxis = 33   ! xz axes
          end select

C rearrange data from anatalrst to anataldata
          call anatal_rearrange_sum(0,0,3,
     &         iDaxis,ntaxis,anatalrst,
     &         np,     ne,      0,    0,     nx, ny, nz,    0, nrst,
     &         eb,ew,    1,1,  1,1,  xm, ym, zm,1,
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

*-----------------------------------------------------------------------

                  write(iot,'(/"#",78("-"))')

           if( inum .eq. 1 .or. ip.gt.np_mxang ) then ! frtati 2021/10/05
                  write(iot,'( "#newpage:")')
               else
                  write(iot,'( " newpage:")')
               end if

           if ( np.gt.np_mxang ) then
             write(iot,'( " SKIPPAGE:")')
             inum = inum - 1
           end if

*-----------------------------------------------------------------------

            if( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) then

              write(iot,'(/a1,"no. =",i3,a1)') cha, inum, cha

               write(iot,'("msuc: {",a1,"huge ",80a1)')
     &                  yen, ( ittle(m)(i:i), i = 1, ittll(m) ), rpa

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
     &            ( itang(m)(i:i),i = 1, itanl(m) )

               else if( itanl(m) .gt. 0 .and. ioe .ne. 1 ) then

                   call terrang(iot,m)

               end if

               if( itsans(m) .gt. 0 ) then

                  call write_sangel(iot,m,0)

               end if

            end if

*-----------------------------------------------------------------------

             write(iot,'("#  n2 = ",i3,"   n1 = ",i3)')
     &            nijaxs2, nijaxs

            if( ittwo(m) .ne. 4 ) then

                write(iot,'( "# ( ( data(i1,i2), i1 = 1, n1 ),",
     &                         " i2 = n2, 1, -1 )")')

            end if

*-----------------------------------------------------------------------

            if( ( ittwo(m) .le. 3 .or. ittwo(m) .ge. 6 ) .and.
     &            igsh .eq. 0 ) then

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
     &                    "  SED          r.err")')

              do ijaxs = 1, nijaxs
               do ijaxs2 = 1, nijaxs2

                  write(iot,'(1p3e11.3,0pf8.4)')
     &                 fgaxs2(ijaxs2)  + axs2del/2d0,
     &                 fgaxs(ijaxs)  + axs1del/2d0,
     &           anataldata(ip,ianataldata,ijaxs+(ijaxs2-1)*nijaxs,1),
     &           anataldata(ip,ianataldata,ijaxs+(ijaxs2-1)*nijaxs,2)

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

                  write(iot,'(1p1000e11.3)')
     &                fgaxs2(ijaxs2)  + axs2del/2d0,
     &         ( anataldata(ip,ianataldata,ijaxs+(ijaxs2-1)*nijaxs,ioe),
     &                ijaxs = 1, nijaxs )

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

               if ( itaxs(m,iax) .eq. 7 ) then ! when axis = xy (x-y)
                  yval = ( fgaxs3(ij(5)) + fgaxs3(ij(5)+1) ) / 2.0d0
               else if ( itaxs(m,iax) .eq. 8 ) then ! when axis = yz (z-y)
                  yval = ( fgaxs3(ij(3)) + fgaxs3(ij(3)+1) ) / 2.0d0
               else if ( itaxs(m,iax) .eq. 9 ) then ! when axis = xz (z-x)
                  yval = ( fgaxs3(ij(4)) + fgaxs3(ij(4)+1) ) / 2.0d0
               end if
               none = 1
               iaxs = 1
               iuni = itgsh(m) * 2 - 1 + igsh
               ires = itres(m)
               widt = rtwid(m)

               call gshow(0,iot,iaxs,iuni,widt,ires,
     &                    inum,igsh,igser,gfnam,igfmn,
     &                nijaxs+1,nijaxs2+1,none,fgaxs,fgaxs2,yval,ixyz(1),
     &                    1,1,krr,vll,itmtr(m,4),itglt(m))

         end if

*-----------------------------------------------------------------------

      if( inolg .eq. 1 .and. ( ittwo(m) .eq. 3 .or. ittwo(m) .eq. 7 )
     &    .and. cmin .gt. 0.0 .and. cmax .gt. cmin .and.
     &    igsh .eq. 0 ) then

      write(iot,'(/"#",78("-"))')
      write(iot,'(
     &"z: xorg(1.03)"/
     &"p: xfac[c2*0.05] form[c1/0.05] notc noxt noyt"/
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
           write(iot,'("y: ",a30)') hsunit(3*(itunt(m)-1)+itsun(m))
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
     &                 /1pe13.4,2x,"$--$",1pe13.4)')
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
          end do   ! do ip = 1, np

*-----------------------------------------------------------------------

          deallocate( fg )
          deallocate( fgaxs )
          deallocate( fgaxs2 )
          deallocate( fgaxs3 )
          deallocate( anataldata )
          deallocate( delvol )

*-----------------------------------------------------------------------

         end if

            call prestart(m,iot) !OBINATA(2012.9.7)

            close(iot)

         if( iteps(m) .ne. 0 ) then

            open(iot, file = fname, status = 'unknown' )
            call a_angel(idasa,fname)

         end if

  900    continue

*-----------------------------------------------------------------------
         deallocate( ixyz )
         deallocate( vl )

         deallocate( vl_x,vl_y,vl_z )
         deallocate (anatalrst)
         deallocate (ew,rdata)

         
      return
      end subroutine anatal_psedxyz

!***********************************************************************
!                                                                      *
! sumover subroutine                                                   *
!                                                                      *
!***********************************************************************

!***********************************************************************
!                                                                      *
      subroutine get_sedreg_tr_sum_data(m,iax,nfile,
     &           ip,ie,ir,rdata)
!                                                                      *
!***********************************************************************

      use TALMOD
!$      use TALMOD0
      use anatallymod, only:ianatalm_sum
! lanatalm_sum(m)     :  start address of one tally at nfile(m)=1
! ianatalm_sum(m,iax,nfile(m)) :  start address of one sumover
! manatalm_sum(m,iax,nfile(m)) :  one sumover length

      implicit double precision (a-h,o-z)


      integer :: m,iax,nfile,ip,ie,ir
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

        call get_sedreg_tr_sum_data_sub(p_sum,
     &     itpan_sum(m,iax),itenm_sum(m,iax),
     &     itrgn_sum(m,iax),
     &     ip,ie,ir,rdata(1,ntf))

      enddo

      return
      end

!***********************************************************************
!                                                                      *
      subroutine get_sedreg_tr_sum_data_sub(tr_sum,
     &          np,ne,nr,ip,ie,ir,rdata)
!                                                                      *
!***********************************************************************

      implicit double precision (a-h,o-z)

      integer :: m,iax,np,ne,nr,ip,ie,ir
      real(8) :: rdata(2)

      real(8) :: tr_sum(np,ne,nr,2)
 
      do i=1,2
        rdata(i)= tr_sum(ip,ie,ir,i)
      enddo

      end

!***********************************************************************
!                                                                      *
      subroutine get_sedrz_tr_sum_data(m,iax,nfile,
     &           ip,ie,ir,iz,rdata)
!                                                                      *
!***********************************************************************

      use TALMOD
!$      use TALMOD0
      use anatallymod, only:ianatalm_sum
! lanatalm_sum(m)     :  start address of one tally at nfile(m)=1
! ianatalm_sum(m,iax,nfile(m)) :  start address of one sumover
! manatalm_sum(m,iax,nfile(m)) :  one sumover length

      implicit double precision (a-h,o-z)


      integer :: m,iax,nfile,ip,ie,ir,iz
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

        call get_sedrz_tr_sum_data_sub(p_sum,
     &     itpan_sum(m,iax),itenm_sum(m,iax),
     &     itrnm_sum(m,iax),itznm_sum(m,iax),
     &     ip,ie,ir,iz,rdata(1,ntf))

      enddo

      return
      end

!***********************************************************************
!                                                                      *
      subroutine get_sedrz_tr_sum_data_sub(tr_sum,
     &          np,ne,nr,nz,ip,ie,ir,iz,rdata)
!                                                                      *
!***********************************************************************

      implicit double precision (a-h,o-z)

      integer :: m,iax,np,ne,nr,nz,ip,ie,ir,iz
      real(8) :: rdata(2)

      real(8) :: tr_sum(np,ne,nr,nz,2)

      do i=1,2
        rdata(i)= tr_sum(ip,ie,ir,iz,i)
      enddo

      end

!***********************************************************************
!                                                                      *
      subroutine get_sedxyz_tr_sum_data(m,iax,nfile,
     &           ip,ie,ix,iy,iz,rdata)
!                                                                      *
!***********************************************************************

      use TALMOD
!$      use TALMOD0
      use anatallymod, only:ianatalm_sum
! lanatalm_sum(m)     :  start address of one tally at nfile(m)=1
! ianatalm_sum(m,iax,nfile(m)) :  start address of one sumover
! manatalm_sum(m,iax,nfile(m)) :  one sumover length

      implicit double precision (a-h,o-z)


      integer :: m,iax,nfile,ip,ie,ix,iy,iz
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

        call get_sedxyz_tr_sum_data_sub(p_sum,
     &     itpan_sum(m,iax),itenm_sum(m,iax),
     &     itxnm_sum(m,iax),itynm_sum(m,iax),itznm_sum(m,iax),
     &     ip,ie,ix,iy,iz,rdata(1,ntf))

      enddo

      return
      end

!***********************************************************************
!                                                                      *
      subroutine get_sedxyz_tr_sum_data_sub(tr_sum,
     &          np,ne,nx,ny,nz,ip,ie,ix,iy,iz,rdata)
!                                                                      *
!***********************************************************************

      implicit double precision (a-h,o-z)

      integer :: m,iax,np,ne,nx,ny,nz,nm,ip,ie,ix,iy,iz
      real(8) :: rdata(2)

      real(8) :: tr_sum(np,ne,nx*ny*nz,2)

      icf(ix,iy,iz) = ix + ( iy - 1 ) * nx + ( iz - 1 ) * nx * ny


      do i=1,2
        rdata(i)= tr_sum(ip,ie,icf(ix,iy,iz),i)
      enddo

      end


